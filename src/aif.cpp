// aif: FracMinHash + Bloom-filter tier sketches for ANI / alignment-fraction
// estimation of (short) contigs against a sketched reference.
//
// Design (see fmh_bloom_filter_theory.pdf and aif_benchmark_report.pdf):
//   * reference: canonical k-mers are FracMinHash-selected into disjoint hash-range
//     tiers, each a cache-blocked Bloom filter with its own FPR (hedged tiers);
//     exact FPRs from fill ratios; spurious-match rate calibrated empirically;
//     k chosen per reference (adaptive k).
//   * query (raw contigs): every selected k-mer is probed; a multi-level HMM
//     (unaligned + aligned identity levels) over the probe hits gives the aligned
//     fraction and the FP-corrected aligned containment theta; ANI = theta^(1/k)
//     with bias correction, or the posterior-mean ANI for heterogeneous contigs.
//
// v1.1 speed work: rolling canonical k-mers into flat buffers, vectorised
// selection hashing, branchless compaction, prefetched filter tests, binned HMM
// (B probes per step, emissions from power tables), buffered FASTA reading,
// O(1) spurious-rate calibration during the build pass, optional threads.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread -o aif aif.cpp
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif
#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#define AIF_AVX512 1
#endif

namespace {

const uint64_t SEED_SEL = 0x9E3779B97F4A7C15ULL;
const uint64_t SEED_BF = 0xC2B2AE3D27D4EB4FULL;
const uint64_t SEED_BIT = 0x165667B19E3779F9ULL;
const uint32_t MAGIC = 0x33464941;  // "AIF3" (packed low-first k-mer encoding)
const int MAXT = 4;                 // maximum number of tiers

inline uint64_t fmix64(uint64_t k) {
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53ULL;
  k ^= k >> 33;
  return k;
}
inline uint64_t mulhi(uint64_t a, uint64_t b) { return (uint64_t)(((unsigned __int128)a * b) >> 64); }

int8_t CODE[256];
inline double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
void init_code() {
  for (int i = 0; i < 256; ++i) CODE[i] = -1;
  CODE['A'] = CODE['a'] = 0;
  CODE['C'] = CODE['c'] = 1;
  CODE['G'] = CODE['g'] = 2;
  CODE['T'] = CODE['t'] = 3;
}

// ---------------------------------------------------------------- FASTA (buffered)
class Fasta {
  // memory-mapped input (regular files); falls back to buffered stdio otherwise
  const char* map = nullptr;
  size_t map_len = 0, mpos = 0, released = 0;
  int fd = -1;
  FILE* f = nullptr;
  std::vector<char> buf;
  size_t beg = 0, end = 0;
  bool eof = false;
  bool have_header = false;
  std::string header;

  bool refill() {
    if (eof) return false;
    if (beg > 0) {
      std::memmove(buf.data(), buf.data() + beg, end - beg);
      end -= beg;
      beg = 0;
    }
    if (end == buf.size()) buf.resize(buf.size() * 2);
    size_t got = fread(buf.data() + end, 1, buf.size() - end, f);
    if (got == 0) eof = true;
    end += got;
    return got > 0;
  }
  // reads one line (without '\n'); returns false at EOF with nothing read
  bool line(const char*& p, size_t& n) {
    for (;;) {
      const char* s = buf.data() + beg;
      const char* nl = (const char*)memchr(s, '\n', end - beg);
      if (nl) {
        p = s;
        n = (size_t)(nl - s);
        beg += n + 1;
        if (n && p[n - 1] == '\r') --n;
        return true;
      }
      if (!refill()) {
        if (end > beg) {  // last line without newline
          p = buf.data() + beg;
          n = end - beg;
          beg = end;
          if (n && p[n - 1] == '\r') --n;
          return true;
        }
        return false;
      }
    }
  }

 public:
  explicit Fasta(const char* path) {
    fd = open(path, O_RDONLY);
    struct stat st;
    if (fd >= 0 && fstat(fd, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) {
      void* p = mmap(nullptr, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
      if (p != MAP_FAILED) {
        map = (const char*)p;
        map_len = (size_t)st.st_size;
        madvise(p, map_len, MADV_SEQUENTIAL);
        return;
      }
    }
    if (fd >= 0) close(fd);
    fd = -1;
    buf.resize(1 << 22);
    f = fopen(path, "rb");
    if (!f) throw std::runtime_error(std::string("cannot open ") + path);
  }
  ~Fasta() {
    if (map) munmap((void*)map, map_len);
    if (fd >= 0) close(fd);
    if (f) fclose(f);
  }
  bool next_mapped(std::string& name, std::string& seq);
  bool next(std::string& name, std::string& seq) {
    if (map) return next_mapped(name, seq);
    seq.clear();
    const char* p;
    size_t n;
    if (!have_header) {
      for (;;) {
        if (!line(p, n)) return false;
        if (n && p[0] == '>') {
          header.assign(p + 1, n - 1);
          break;
        }
      }
    }
    have_header = false;
    size_t e = header.find_first_of(" \t");
    name.assign(header, 0, e);
    while (line(p, n)) {
      if (n && p[0] == '>') {
        header.assign(p + 1, n - 1);
        have_header = true;
        break;
      }
      seq.append(p, n);
    }
    return true;
  }
};

bool Fasta::next_mapped(std::string& name, std::string& seq) {
  seq.clear();
  const char* p = map + mpos;
  const char* end = map + map_len;
  while (p < end && *p != '>') {  // skip to the next header
    const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
    p = nl ? nl + 1 : end;
  }
  if (p >= end) {
    mpos = map_len;
    return false;
  }
  const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
  const char* he = nl ? nl : end;
  const char* hs = p + 1;
  const char* ws = hs;
  while (ws < he && *ws != ' ' && *ws != '\t' && *ws != '\r') ++ws;
  name.assign(hs, (size_t)(ws - hs));
  p = nl ? nl + 1 : end;
  // sequence lines until the next header
  const char* nxt = p;
  while (nxt < end) {
    const char* q = (const char*)memchr(nxt, '>', (size_t)(end - nxt));
    if (!q) {
      nxt = end;
      break;
    }
    if (q == map || q[-1] == '\n') {
      nxt = q;
      break;
    }
    nxt = q + 1;
  }
  seq.reserve((size_t)(nxt - p));
  while (p < nxt) {
    const char* e2 = (const char*)memchr(p, '\n', (size_t)(nxt - p));
    const char* le = e2 ? e2 : nxt;
    size_t len = (size_t)(le - p);
    if (len && p[len - 1] == '\r') --len;
    seq.append(p, len);
    p = e2 ? e2 + 1 : nxt;
  }
  mpos = (size_t)(nxt - map);
  // release consumed pages so the resident set stays small on large inputs
  const size_t pg = 1 << 22;
  if (mpos - released >= pg) {
    size_t upto = (mpos / 4096) * 4096;
    madvise((void*)(map + released), upto - released, MADV_DONTNEED);
    released = upto;
  }
  return true;
}

[[maybe_unused]] std::vector<std::string> load_fasta(const char* fasta) {
  std::vector<std::string> seqs;
  std::string name, seq;
  Fasta fa(fasta);
  while (fa.next(name, seq)) seqs.push_back(seq);
  return seqs;
}

// ---------------------------------------------------------------- packed 2-bit sequences
// Encoding (AIF3): code = (ascii >> 1) & 3, i.e. A=0 C=1 T=2 G=3; complement = code ^ 2.
// Base i is stored at bits 2*(i%32) of word i/32 ("low-first"); the forward k-mer at
// start i has base i in its lowest two bits.
struct Packed {
  std::vector<uint64_t> w;     // words (+2 zero words of slack)
  std::vector<uint8_t> badk;   // per start: 1 if the k-mer overlaps a non-ACGT base (empty if none)
  std::vector<uint32_t> bad;   // positions of non-ACGT bases
  size_t L = 0;
  void set_k(int k) {          // rebuild badk for a given k
    badk.clear();
    if (bad.empty() || L < (size_t)k) return;
    badk.assign(L - k + 1 + 64, 0);
    for (uint32_t b : bad) {
      size_t lo = b + 1 >= (uint32_t)k ? b + 1 - k : 0, hi = std::min<size_t>(b, L - k);
      for (size_t i = lo; i <= hi; ++i) badk[i] = 1;
    }
  }
};

inline uint64_t zero_bytes(uint64_t y) {  // 0x80 in every byte of y that is zero (exact)
  const uint64_t m = 0x7F7F7F7F7F7F7F7FULL;
  return ~(((y & m) + m) | y | m);
}

void pack_seq(const uint8_t* s, size_t L, int k, Packed& P) {
  P.L = L;
  const size_t nw = (L + 31) / 32;
  P.w.assign(nw + 2, 0);
  P.badk.clear();
  std::vector<uint32_t>& bad = P.bad;
  bad.clear();
  for (size_t w = 0; w < nw; ++w) {
    const size_t i0 = w * 32;
    uint64_t acc = 0;
    if (i0 + 32 <= L) {
      for (int q = 0; q < 4; ++q) {
        uint64_t v;
        std::memcpy(&v, s + i0 + 8 * q, 8);
        const uint64_t lc = v | 0x2020202020202020ULL;
        const uint64_t ok = zero_bytes(lc ^ 0x6161616161616161ULL) | zero_bytes(lc ^ 0x6363636363636363ULL) |
                            zero_bytes(lc ^ 0x6767676767676767ULL) | zero_bytes(lc ^ 0x7474747474747474ULL);
        if (__builtin_expect(ok != 0x8080808080808080ULL, 0))
          for (int b = 0; b < 8; ++b)
            if (CODE[s[i0 + 8 * q + b]] < 0) bad.push_back((uint32_t)(i0 + 8 * q + b));
        const uint64_t c = (v >> 1) & 0x0303030303030303ULL;
#if defined(__BMI2__)
        acc |= _pext_u64(c, 0x0303030303030303ULL) << (16 * q);
#else
        uint64_t pk = 0;
        for (int b = 0; b < 8; ++b) pk |= ((c >> (8 * b)) & 3) << (2 * b);
        acc |= pk << (16 * q);
#endif
      }
    } else {
      for (size_t i = i0; i < L; ++i) {
        if (CODE[s[i]] < 0) bad.push_back((uint32_t)i);
        acc |= (uint64_t)((s[i] >> 1) & 3) << (2 * (i - i0));
      }
    }
    P.w[w] = acc;
  }
  P.set_k(k);
}

inline uint64_t kmer_at(const uint64_t* w, size_t i, uint64_t mask) {
  const size_t q = i >> 5, r = 2 * (i & 31);
  return ((w[q] >> r) | ((w[q + 1] << 1) << (63 - r))) & mask;
}
inline uint64_t rc_kmer(uint64_t f, int k) {
  uint64_t x = f ^ 0xAAAAAAAAAAAAAAAAULL;
  x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
  x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
  x = __builtin_bswap64(x);
  return x >> (64 - 2 * k);
}

#ifdef AIF_AVX512
inline __m512i fmix64_v(__m512i k) {
  k = _mm512_xor_si512(k, _mm512_srli_epi64(k, 33));
  k = _mm512_mullo_epi64(k, _mm512_set1_epi64((long long)0xff51afd7ed558ccdULL));
  k = _mm512_xor_si512(k, _mm512_srli_epi64(k, 33));
  k = _mm512_mullo_epi64(k, _mm512_set1_epi64((long long)0xc4ceb9fe1a85ec53ULL));
  k = _mm512_xor_si512(k, _mm512_srli_epi64(k, 33));
  return k;
}
#endif

// selected (h(x) < hmax) canonical k-mers of a packed sequence, in position order;
// out arrays need L - k + 1 + 8 entries
size_t select_kmers(const Packed& P, int k, uint64_t hmax, size_t start, size_t end, uint64_t* out_x, uint64_t* out_h,
                    uint32_t* out_pos) {
  if (P.L < (size_t)k) return 0;
  const size_t nst = std::min(P.L - k + 1, end);  // start must be a multiple of 32
  const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const uint64_t* w = P.w.data();
  const uint8_t* badk = P.badk.empty() ? nullptr : P.badk.data();
  size_t m = 0;
#ifdef AIF_AVX512
  alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                              23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                              39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                              55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
  const __m512i bsw = _mm512_load_si512((const void*)BSW);
  const __m512i vmask = _mm512_set1_epi64((long long)mask);
  const __m512i vcomp = _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL);
  const __m512i m33 = _mm512_set1_epi64((long long)0x3333333333333333ULL);
  const __m512i m0f = _mm512_set1_epi64((long long)0x0F0F0F0F0F0F0F0FULL);
  const __m512i vseed = _mm512_set1_epi64((long long)SEED_SEL);
  const __m512i vhmax = _mm512_set1_epi64((long long)hmax);
  const __m512i lane = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  const __m512i v63 = _mm512_set1_epi64(63);
  const __m128i rcs = _mm_cvtsi32_si128(64 - 2 * k);
  alignas(64) uint64_t tc[32], th[32];
  for (size_t base = start; base < nst; base += 32) {
    const __m512i W0 = _mm512_set1_epi64((long long)w[base >> 5]);
    const __m512i W1 = _mm512_set1_epi64((long long)(w[(base >> 5) + 1] << 1));
    uint32_t selm = 0;
    for (int r0 = 0; r0 < 32; r0 += 8) {
      const __m512i rv = _mm512_add_epi64(lane, _mm512_set1_epi64(r0));
      const __m512i sh = _mm512_add_epi64(rv, rv);
      __m512i f = _mm512_or_si512(_mm512_srlv_epi64(W0, sh), _mm512_sllv_epi64(W1, _mm512_sub_epi64(v63, sh)));
      f = _mm512_and_si512(f, vmask);
      __m512i x = _mm512_xor_si512(f, vcomp);
      x = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 2), m33), _mm512_slli_epi64(_mm512_and_si512(x, m33), 2));
      x = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(x, m0f), 4));
      x = _mm512_shuffle_epi8(x, bsw);
      x = _mm512_srl_epi64(x, rcs);
      const __m512i c = _mm512_min_epu64(f, x);
      const __m512i h = fmix64_v(_mm512_xor_si512(c, vseed));
      _mm512_store_si512((void*)(tc + r0), c);
      _mm512_store_si512((void*)(th + r0), h);
      selm |= (uint32_t)_mm512_cmplt_epu64_mask(h, vhmax) << r0;
    }
    if (base + 32 > nst) selm &= (nst - base >= 32) ? 0xFFFFFFFFu : ((1u << (nst - base)) - 1);
    if (badk && selm)
      for (uint32_t q = selm; q; q &= q - 1) {
        int j = __builtin_ctz(q);
        if (badk[base + j]) selm &= ~(1u << j);
      }
    while (selm) {
      int j = __builtin_ctz(selm);
      selm &= selm - 1;
      out_x[m] = tc[j];
      out_h[m] = th[j];
      out_pos[m] = (uint32_t)(base + j);
      ++m;
    }
  }
#else
  for (size_t i = start; i < nst; ++i) {
    uint64_t f = kmer_at(w, i, mask), r = rc_kmer(f, k);
    uint64_t c = f < r ? f : r, h = fmix64(c ^ SEED_SEL);
    out_x[m] = c;
    out_h[m] = h;
    out_pos[m] = (uint32_t)i;
    m += (h < hmax) & (badk == nullptr || !badk[i]);
  }
#endif
  return m;
}

#ifdef AIF_AVX512
// positions of selected k-mers only (compressed 16 at a time); out needs L - k + 1 + 32 entries
size_t select_positions(const Packed& P, int k, uint64_t hmax, uint32_t* out_pos, size_t start = 0,
                        size_t end = ~(size_t)0) {
  if (P.L < (size_t)k) return 0;
  const size_t nst = std::min(P.L - k + 1, end);  // start must be a multiple of 32
  const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const uint64_t* w = P.w.data();
  const uint8_t* badk = P.badk.empty() ? nullptr : P.badk.data();
  alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                              23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                              39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                              55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
  const __m512i bsw = _mm512_load_si512((const void*)BSW);
  const __m512i vmask = _mm512_set1_epi64((long long)mask);
  const __m512i vcomp = _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL);
  const __m512i m33 = _mm512_set1_epi64((long long)0x3333333333333333ULL);
  const __m512i m0f = _mm512_set1_epi64((long long)0x0F0F0F0F0F0F0F0FULL);
  const __m512i vseed = _mm512_set1_epi64((long long)SEED_SEL);
  const __m512i vhmax = _mm512_set1_epi64((long long)hmax);
  const __m512i lane = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  const __m512i v63 = _mm512_set1_epi64(63);
  const __m128i rcs = _mm_cvtsi32_si128(64 - 2 * k);
  const __m512i iota16 = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
  size_t m = 0;
  for (size_t base = start; base < nst; base += 32) {
    const __m512i W0 = _mm512_set1_epi64((long long)w[base >> 5]);
    const __m512i W1 = _mm512_set1_epi64((long long)(w[(base >> 5) + 1] << 1));
    uint32_t selm = 0;
    for (int r0 = 0; r0 < 32; r0 += 8) {
      const __m512i rv = _mm512_add_epi64(lane, _mm512_set1_epi64(r0));
      const __m512i sh = _mm512_add_epi64(rv, rv);
      __m512i f = _mm512_or_si512(_mm512_srlv_epi64(W0, sh), _mm512_sllv_epi64(W1, _mm512_sub_epi64(v63, sh)));
      f = _mm512_and_si512(f, vmask);
      __m512i x = _mm512_xor_si512(f, vcomp);
      x = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 2), m33), _mm512_slli_epi64(_mm512_and_si512(x, m33), 2));
      x = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(x, m0f), 4));
      x = _mm512_srl_epi64(_mm512_shuffle_epi8(x, bsw), rcs);
      const __m512i h = fmix64_v(_mm512_xor_si512(_mm512_min_epu64(f, x), vseed));
      selm |= (uint32_t)_mm512_cmplt_epu64_mask(h, vhmax) << r0;
    }
    if (base + 32 > nst) selm &= (1u << (nst - base)) - 1;
    if (badk && selm)
      for (uint32_t q = selm; q; q &= q - 1) {
        int j = __builtin_ctz(q);
        if (badk[base + j]) selm &= ~(1u << j);
      }
    if (!selm) continue;
    const __m512i pb = _mm512_add_epi32(_mm512_set1_epi32((int)base), iota16);
    const __m512i pb16 = _mm512_add_epi32(pb, _mm512_set1_epi32(16));
    const __mmask16 lo = (__mmask16)(selm & 0xFFFF), hi = (__mmask16)(selm >> 16);
    _mm512_storeu_si512((void*)(out_pos + m), _mm512_maskz_compress_epi32(lo, pb));
    m += (size_t)__builtin_popcount(lo);
    _mm512_storeu_si512((void*)(out_pos + m), _mm512_maskz_compress_epi32(hi, pb16));
    m += (size_t)__builtin_popcount(hi);
  }
  return m;
}
#endif

// ---------------------------------------------------------------- Bloom tier
struct Tier {
  uint64_t lo = 0, hi = 0;  // selection-hash range [lo, hi)
  uint32_t nh = 1;          // hash functions
  uint64_t nblocks = 1;     // 512-bit blocks
  uint64_t n_ins = 0;       // insertions (with multiplicity)
  double eps = 0;           // exact FPR of a random non-member
  std::vector<uint64_t> bits;

  double s() const { return (double)(hi - lo) / 18446744073709551616.0; }
  inline uint64_t block_index(uint64_t h) const { return ((h >> 32) * nblocks) >> 32; }
  inline uint64_t* block_of(uint64_t x, uint64_t& h) {
    h = fmix64(x ^ SEED_BF);
    return &bits[block_index(h) * 8];
  }
  inline const uint64_t* block_of(uint64_t x, uint64_t& h) const {
    h = fmix64(x ^ SEED_BF);
    return &bits[block_index(h) * 8];
  }
  inline void set_at(uint64_t* blk, uint64_t h) const {
    uint64_t h2 = fmix64(h ^ SEED_BIT);
    for (uint32_t j = 0; j < nh; ++j) {
      if (j && j % 7 == 0) h2 = fmix64(h2 + j);
      uint32_t bp = (uint32_t)((h2 >> (9 * (j % 7))) & 511);
      blk[bp >> 6] |= 1ULL << (bp & 63);
    }
  }
  // branch-free membership test (same bit positions as set_at)
  inline bool contains_at(const uint64_t* blk, uint64_t h) const {
    uint64_t h2 = fmix64(h ^ SEED_BIT);
    uint64_t acc = 1;
    const uint32_t n1 = nh < 7 ? nh : 7;
    for (uint32_t j = 0; j < n1; ++j) {
      uint32_t bp = (uint32_t)((h2 >> (9 * j)) & 511);
      acc &= blk[bp >> 6] >> (bp & 63);
    }
    for (uint32_t j = 7; j < nh; ++j) {
      if (j % 7 == 0) h2 = fmix64(h2 + j);
      uint32_t bp = (uint32_t)((h2 >> (9 * (j % 7))) & 511);
      acc &= blk[bp >> 6] >> (bp & 63);
    }
    return acc & 1;
  }
  bool contains(uint64_t x) const {
    uint64_t h;
    const uint64_t* blk = block_of(x, h);
    return contains_at(blk, h);
  }
  void compute_eps() {
    double acc = 0;
    for (uint64_t b = 0; b < nblocks; ++b) {
      int ones = 0;
      for (int w = 0; w < 8; ++w) ones += __builtin_popcountll(bits[b * 8 + w]);
      acc += std::pow(ones / 512.0, (double)nh);
    }
    eps = acc / (double)nblocks;
  }
};

// exact Bloom cost (bits/element) minimised over integer h
double c_bf(double eps) {
  double best = 1e300;
  for (int h = 1; h <= 30; ++h) best = std::min(best, -h / std::log(1 - std::pow(eps, 1.0 / h)));
  return best;
}

struct Sketch {
  int k = 21;
  uint64_t total_bases = 0;
  uint64_t n_kmers = 0;
  double budget = 0;
  double sigma_emp = -1;  // empirically calibrated spurious-match rate (<0: unknown)
  std::vector<Tier> tiers;

  double sigma() const { return -std::expm1(-2.0 * (double)n_kmers / std::pow(4.0, k)); }
  double sigma_used() const { return sigma_emp >= 0 ? sigma_emp : sigma(); }
  size_t filter_bytes() const {
    size_t b = 0;
    for (auto& t : tiers) b += t.bits.size() * 8;
    return b;
  }
  inline int tier_of(uint64_t hs) const {
    int t = 0;
    while (t + 1 < (int)tiers.size() && hs >= tiers[t].hi) ++t;
    return t;
  }
  void save(const char* path) const {
    FILE* f = fopen(path, "wb");
    if (!f) throw std::runtime_error("cannot write sketch");
    uint32_t nt = tiers.size();
    fwrite(&MAGIC, 4, 1, f);
    fwrite(&k, 4, 1, f);
    fwrite(&total_bases, 8, 1, f);
    fwrite(&n_kmers, 8, 1, f);
    fwrite(&budget, 8, 1, f);
    fwrite(&sigma_emp, 8, 1, f);
    fwrite(&nt, 4, 1, f);
    for (auto& t : tiers) {
      fwrite(&t.lo, 8, 1, f);
      fwrite(&t.hi, 8, 1, f);
      fwrite(&t.nh, 4, 1, f);
      fwrite(&t.nblocks, 8, 1, f);
      fwrite(&t.n_ins, 8, 1, f);
      fwrite(&t.eps, 8, 1, f);
      fwrite(t.bits.data(), 8, t.bits.size(), f);
    }
    fclose(f);
  }
  void load(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) throw std::runtime_error("cannot read sketch");
    uint32_t magic, nt;
    auto rd = [&](void* p, size_t sz) {
      if (fread(p, sz, 1, f) != 1) throw std::runtime_error("truncated sketch");
    };
    rd(&magic, 4);
    if (magic != MAGIC) throw std::runtime_error("not an aif sketch");
    rd(&k, 4);
    rd(&total_bases, 8);
    rd(&n_kmers, 8);
    rd(&budget, 8);
    rd(&sigma_emp, 8);
    rd(&nt, 4);
    if (nt == 0 || nt > (uint32_t)MAXT) throw std::runtime_error("bad tier count");
    tiers.resize(nt);
    for (auto& t : tiers) {
      rd(&t.lo, 8);
      rd(&t.hi, 8);
      rd(&t.nh, 4);
      rd(&t.nblocks, 8);
      rd(&t.n_ins, 8);
      rd(&t.eps, 8);
      t.bits.resize(t.nblocks * 8);
      rd(t.bits.data(), 8 * t.bits.size());
    }
    fclose(f);
  }
};

// ---------------------------------------------------------------- sketching
struct TierSpec {
  double eps, frac;
};

std::vector<TierSpec> parse_tiers(const std::string& s) {
  std::vector<TierSpec> v;
  size_t i = 0;
  while (i < s.size()) {
    size_t j = s.find(',', i);
    if (j == std::string::npos) j = s.size();
    std::string tok = s.substr(i, j - i);
    size_t c = tok.find(':');
    if (c == std::string::npos) throw std::runtime_error("tier spec must be eps:frac[,eps:frac]");
    v.push_back({atof(tok.substr(0, c).c_str()), atof(tok.substr(c + 1).c_str())});
    i = j + 1;
  }
  if (v.empty() || v.size() > (size_t)MAXT) throw std::runtime_error("1 to 4 tiers supported");
  return v;
}

#ifdef AIF_AVX512
// insert the k-mers at the given (selected) positions: 8 per step, gather-OR-scatter per hash
// function; lanes whose target words collide fall back to scalar updates for that round
void insert_positions_v(Sketch& sk, const Packed& pk, int k, const uint32_t* pos32, size_t m) {
  const int T = (int)sk.tiers.size();
  Tier& t0 = sk.tiers[0];
  Tier& t1 = sk.tiers[T - 1];
  long long* base0 = (long long*)t0.bits.data();
  const long long delta1 = (long long)((long long*)t1.bits.data() - base0);
  const long long* wbase64 = (const long long*)pk.w.data();
  const uint64_t kmask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
  alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                              23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                              39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                              55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
  const __m512i bsw = _mm512_load_si512((const void*)BSW);
  const __m512i vmask = _mm512_set1_epi64((long long)kmask), vcomp = _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL);
  const __m512i m33 = _mm512_set1_epi64((long long)0x3333333333333333ULL), m0f = _mm512_set1_epi64((long long)0x0F0F0F0F0F0F0F0FULL);
  const __m512i vsel = _mm512_set1_epi64((long long)SEED_SEL), v63 = _mm512_set1_epi64(63), v31 = _mm512_set1_epi64(31);
  const __m128i rcs = _mm_cvtsi32_si128(64 - 2 * k);
  const __m512i vhi0 = _mm512_set1_epi64((long long)t0.hi);
  const __m512i vnb0 = _mm512_set1_epi64((long long)t0.nblocks), vnb1 = _mm512_set1_epi64((long long)t1.nblocks);
  const __m512i voff1 = _mm512_set1_epi64(delta1), zero = _mm512_setzero_si512();
  const __m512i sbf = _mm512_set1_epi64((long long)SEED_BF), sbit = _mm512_set1_epi64((long long)SEED_BIT);
  const __m512i m511 = _mm512_set1_epi64(511), m63 = _mm512_set1_epi64(63), one = _mm512_set1_epi64(1);
  const uint32_t nh0 = t0.nh, nh1 = t1.nh, nhmax = std::max(nh0, nh1);
  alignas(64) long long ix[8], bv[8];
  for (size_t j = 0; j < m; j += 8) {
    const __mmask8 tail = (m - j >= 8) ? (__mmask8)0xFF : (__mmask8)((1u << (m - j)) - 1);
    const __m512i pos = _mm512_cvtepu32_epi64(_mm256_maskz_loadu_epi32(tail, pos32 + j));
    const __m512i wi = _mm512_srli_epi64(pos, 5);
    const __m512i A = _mm512_mask_i64gather_epi64(zero, tail, wi, wbase64, 8);
    const __m512i Bw = _mm512_slli_epi64(_mm512_mask_i64gather_epi64(zero, tail, _mm512_add_epi64(wi, one), wbase64, 8), 1);
    const __m512i sh = _mm512_slli_epi64(_mm512_and_si512(pos, v31), 1);
    __m512i f = _mm512_and_si512(_mm512_or_si512(_mm512_srlv_epi64(A, sh), _mm512_sllv_epi64(Bw, _mm512_sub_epi64(v63, sh))), vmask);
    __m512i r = _mm512_xor_si512(f, vcomp);
    r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 2), m33), _mm512_slli_epi64(_mm512_and_si512(r, m33), 2));
    r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(r, m0f), 4));
    r = _mm512_srl_epi64(_mm512_shuffle_epi8(r, bsw), rcs);
    const __m512i x = _mm512_min_epu64(f, r);
    const __m512i h = fmix64_v(_mm512_xor_si512(x, vsel));
    const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_cmpge_epu64_mask(h, vhi0) & tail) : (__mmask8)0;
    t1.n_ins += (uint64_t)__builtin_popcount((unsigned)m1);
    t0.n_ins += (uint64_t)__builtin_popcount((unsigned)(tail & ~m1));
    if (T == 1) t0.n_ins -= (uint64_t)__builtin_popcount((unsigned)m1);
    const __m512i bh = fmix64_v(_mm512_xor_si512(x, sbf));
    const __m512i nb = _mm512_mask_blend_epi64(m1, vnb0, vnb1);
    const __m512i blk = _mm512_srli_epi64(_mm512_mul_epu32(_mm512_srli_epi64(bh, 32), nb), 32);
    const __m512i wb = _mm512_add_epi64(_mm512_slli_epi64(blk, 3), _mm512_mask_blend_epi64(m1, zero, voff1));
    __m512i h2 = fmix64_v(_mm512_xor_si512(bh, sbit));
    for (uint32_t q = 0; q < nhmax; ++q) {
      if (q && q % 7 == 0) h2 = fmix64_v(_mm512_add_epi64(h2, _mm512_set1_epi64(q)));
      const __m512i bp = _mm512_and_si512(_mm512_srl_epi64(h2, _mm_cvtsi32_si128(9 * (q % 7))), m511);
      const __m512i idx = _mm512_add_epi64(wb, _mm512_srli_epi64(bp, 6));
      const __m512i bitv = _mm512_sllv_epi64(one, _mm512_and_si512(bp, m63));
      const __mmask8 am = (__mmask8)(((q < nh0) ? (__mmask8)(~m1) : 0) | ((q < nh1) ? m1 : 0)) & tail;
      const __m512i cf = _mm512_maskz_conflict_epi64(am, idx);
      if (_mm512_mask_test_epi64_mask(am, cf, cf)) {  // colliding lanes: scalar updates
        _mm512_store_si512((void*)ix, idx);
        _mm512_store_si512((void*)bv, bitv);
        for (int l = 0; l < 8; ++l)
          if ((am >> l) & 1) base0[ix[l]] |= bv[l];
      } else {
        const __m512i g = _mm512_mask_i64gather_epi64(zero, am, idx, base0, 8);
        _mm512_mask_i64scatter_epi64(base0, am, idx, _mm512_or_si512(g, bitv), 8);
      }
    }
  }
}
#endif

Sketch build_sketch(std::vector<Packed>& packs, int k, double budget, const std::vector<TierSpec>& spec,
                    size_t calib_target = 1000000, uint32_t max_hashes = 6) {
  Sketch sk;
  sk.k = k;
  sk.budget = budget;
  for (auto& pk : packs) {
    sk.total_bases += pk.L;
    if (pk.L >= (size_t)k) sk.n_kmers += pk.L - k + 1;
  }
  double fsum = 0;
  for (auto& t : spec) fsum += t.frac;
  std::vector<double> s(spec.size());
  double stot = 0;
  for (size_t i = 0; i < spec.size(); ++i) {
    s[i] = budget * (spec[i].frac / fsum) / c_bf(spec[i].eps);
    stot += s[i];
  }
  if (stot > 1.0)
    for (auto& x : s) x /= stot;
  sk.tiers.resize(spec.size());
  long double lo = 0;
  const long double TWO64 = 18446744073709551616.0L;
  for (size_t i = 0; i < spec.size(); ++i) {
    Tier& t = sk.tiers[i];
    t.lo = (uint64_t)std::min(lo, TWO64 - 1);
    lo += (long double)s[i] * TWO64;
    t.hi = (uint64_t)std::min(lo, TWO64 - 1);
    double bits = budget * (spec[i].frac / fsum) * (double)sk.total_bases;
    t.nblocks = std::max<uint64_t>(1, (uint64_t)(bits / 512.0));
    t.bits.assign(t.nblocks * 8, 0);
    double m = (double)t.nblocks * 512.0, n = std::max<double>(1.0, t.s() * (double)sk.n_kmers);
    double best = 2;
    for (uint32_t h = 1; h <= std::max<uint32_t>(1, max_hashes); ++h) {
      double hh = (double)h;
      double fpr = std::pow(1 - std::exp(-hh * n / m), hh);
      if (fpr < best) {
        best = fpr;
        t.nh = h;
      }
    }
  }
  const uint64_t hmax = sk.tiers.back().hi;
  const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
  // calibration: single-substitution neighbours of random reference k-mers
  std::vector<uint64_t> cal_x;
  std::vector<uint8_t> cal_t;
  uint64_t rng = 0x2545F4914F6CDD1DULL;
  auto next = [&]() {
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return rng;
  };
  const size_t CH = 1 << 16;
  std::vector<uint64_t> sx(CH + 8), shv(CH + 8), bh(CH + 8);
  std::vector<uint32_t> spos(CH + 40);
  cal_x.reserve(calib_target / 4 + 64);
  cal_t.reserve(calib_target / 4 + 64);
  std::vector<uint64_t*> blk(CH + 8);
  std::vector<uint8_t> bt(CH + 8);
  for (auto& pk : packs) {
    if (pk.L < (size_t)k) continue;
    pk.set_k(k);
    const size_t nst = pk.L - k + 1;
    for (size_t c0 = 0; c0 < nst; c0 += CH) {
#ifdef AIF_AVX512
      if (sk.tiers.size() <= 2) {
        const size_t m = select_positions(pk, k, hmax, spos.data(), c0, c0 + CH);
        insert_positions_v(sk, pk, k, spos.data(), m);
        continue;
      }
#endif
      size_t m = select_kmers(pk, k, hmax, c0, c0 + CH, sx.data(), shv.data(), spos.data());
      for (size_t j = 0; j < m; ++j) {
        int t = sk.tier_of(shv[j]);
        bt[j] = (uint8_t)t;
        blk[j] = sk.tiers[t].block_of(sx[j], bh[j]);
        __builtin_prefetch(blk[j], 1);
      }
      for (size_t j = 0; j < m; ++j) {
        sk.tiers[bt[j]].set_at(blk[j], bh[j]);
        ++sk.tiers[bt[j]].n_ins;
      }
    }
    // stratified calibration samples, 8 independent xorshift streams (identical in the
    // AVX-512 and scalar paths, so sketches do not depend on the instruction set)
    if (nst < 64) continue;
    const size_t G = std::min<size_t>(nst / 8, std::max<size_t>(1, (size_t)((double)calib_target * (double)nst /
                                                   std::max<double>(1.0, (double)sk.n_kmers) / 8.0)));
    const uint64_t stride = std::max<uint64_t>(1, nst / (8 * G));
    const uint64_t* w = pk.w.data();
    alignas(64) uint64_t st[8];
    for (int l = 0; l < 8; ++l) st[l] = fmix64(next() + (uint64_t)l) | 1;
    auto emit = [&](uint64_t pos, uint64_t x, uint64_t h) {
      if (h >= hmax || pos >= nst) return;
      if (!pk.badk.empty() && pk.badk[pos]) return;
      cal_x.push_back(x);
      cal_t.push_back((uint8_t)sk.tier_of(h));
    };
#ifdef AIF_AVX512
    {
      alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                                  23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                                  39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                                  55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
      const __m512i bsw = _mm512_load_si512((const void*)BSW);
      const __m512i vmask = _mm512_set1_epi64((long long)mask), vcomp = _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL);
      const __m512i m33 = _mm512_set1_epi64((long long)0x3333333333333333ULL), m0f = _mm512_set1_epi64((long long)0x0F0F0F0F0F0F0F0FULL);
      const __m512i vsel = _mm512_set1_epi64((long long)SEED_SEL), vhmax = _mm512_set1_epi64((long long)hmax);
      const __m512i v63 = _mm512_set1_epi64(63), v31 = _mm512_set1_epi64(31), one = _mm512_set1_epi64(1);
      const __m512i lo32 = _mm512_set1_epi64(0xFFFFFFFFLL), vk = _mm512_set1_epi64(k), v3 = _mm512_set1_epi64(3);
      const __m512i vstride = _mm512_set1_epi64((long long)stride), vnst1 = _mm512_set1_epi64((long long)(nst - 1));
      const __m512i lane = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
      const __m128i rcs = _mm_cvtsi32_si128(64 - 2 * k);
      __m512i S = _mm512_load_si512((const void*)st);
      auto xs = [](__m512i v) {
        v = _mm512_xor_si512(v, _mm512_slli_epi64(v, 13));
        v = _mm512_xor_si512(v, _mm512_srli_epi64(v, 7));
        return _mm512_xor_si512(v, _mm512_slli_epi64(v, 17));
      };
      alignas(64) uint64_t ox[8], oh[8], op[8];
      for (size_t g = 0; g < G; ++g) {
        S = xs(S);
        const __m512i rr = S;
        S = xs(S);
        const __m512i r2 = S;
        __m512i pos = _mm512_add_epi64(_mm512_mullo_epi64(_mm512_add_epi64(_mm512_set1_epi64((long long)(8 * g)), lane), vstride),
                                       _mm512_srli_epi64(_mm512_mul_epu32(_mm512_srli_epi64(rr, 32), vstride), 32));
        pos = _mm512_min_epu64(pos, vnst1);
        const __m512i off2 = _mm512_slli_epi64(_mm512_srli_epi64(_mm512_mul_epu32(_mm512_and_si512(r2, lo32), vk), 32), 1);
        const __m512i sub = _mm512_add_epi64(one, _mm512_srli_epi64(_mm512_mul_epu32(_mm512_srli_epi64(r2, 32), v3), 32));
        const __m512i wi = _mm512_srli_epi64(pos, 5);
        const __m512i A = _mm512_i64gather_epi64(wi, (const long long*)w, 8);
        const __m512i Bw = _mm512_slli_epi64(_mm512_i64gather_epi64(_mm512_add_epi64(wi, one), (const long long*)w, 8), 1);
        const __m512i sh = _mm512_slli_epi64(_mm512_and_si512(pos, v31), 1);
        __m512i f = _mm512_and_si512(_mm512_or_si512(_mm512_srlv_epi64(A, sh), _mm512_sllv_epi64(Bw, _mm512_sub_epi64(v63, sh))), vmask);
        f = _mm512_xor_si512(f, _mm512_sllv_epi64(sub, off2));
        __m512i r = _mm512_xor_si512(f, vcomp);
        r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 2), m33), _mm512_slli_epi64(_mm512_and_si512(r, m33), 2));
        r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(r, m0f), 4));
        r = _mm512_srl_epi64(_mm512_shuffle_epi8(r, bsw), rcs);
        const __m512i x = _mm512_min_epu64(f, r);
        const __m512i h = fmix64_v(_mm512_xor_si512(x, vsel));
        __mmask8 sel = _mm512_cmplt_epu64_mask(h, vhmax);
        if (!sel) continue;
        _mm512_store_si512((void*)ox, x);
        _mm512_store_si512((void*)oh, h);
        _mm512_store_si512((void*)op, pos);
        for (; sel; sel &= (__mmask8)(sel - 1)) {
          int l = __builtin_ctz((unsigned)sel);
          emit(op[l], ox[l], oh[l]);
        }
      }
    }
#else
    for (size_t g = 0; g < G; ++g)
      for (int l = 0; l < 8; ++l) {
        auto xs1 = [](uint64_t v) {
          v ^= v << 13;
          v ^= v >> 7;
          return v ^ (v << 17);
        };
        st[l] = xs1(st[l]);
        const uint64_t rr = st[l];
        st[l] = xs1(st[l]);
        const uint64_t r2 = st[l];
        uint64_t pos = (8 * g + l) * stride + (((rr >> 32) * (stride & 0xFFFFFFFFULL)) >> 32);
        pos = std::min<uint64_t>(pos, nst - 1);
        const int off = (int)(((r2 & 0xFFFFFFFFULL) * (uint64_t)k) >> 32);
        const uint64_t sub = 1 + (((r2 >> 32) * 3) >> 32);
        const uint64_t f = kmer_at(w, pos, mask) ^ (sub << (2 * off));
        const uint64_t rc = rc_kmer(f, k);
        const uint64_t x = f < rc ? f : rc;
        emit(pos, x, fmix64(x ^ SEED_SEL));
      }
#endif
  }
  for (auto& t : sk.tiers) t.compute_eps();
  // spurious rate: inverse-variance combination of per-tier estimates
  double N[MAXT] = {0}, Y[MAXT] = {0};
  for (size_t i = 0; i < cal_x.size(); ++i) {
    N[cal_t[i]] += 1;
    Y[cal_t[i]] += sk.tiers[cal_t[i]].contains(cal_x[i]);
  }
  double num = 0, den = 0;
  for (size_t t = 0; t < sk.tiers.size(); ++t) {
    if (N[t] < 1) continue;
    double e = sk.tiers[t].eps;
    double est = (Y[t] / N[t] - e) / (1 - e);
    double w = N[t] * (1 - e) * (1 - e) / std::max(e * (1 - e), 1e-9);
    num += w * est;
    den += w;
  }
  sk.sigma_emp = den > 0 ? std::max(0.0, num / den) : -1;
  return sk;
}

// ---------------------------------------------------------------- estimation
struct Probe {
  uint32_t pos;
  uint8_t tier;
  uint8_t hit;
};

struct Params {
  double rho_u = 3e-6;       // per-base switch rate aligned <-> unaligned
  double rho_h = 1e-4;       // per-base switch rate between aligned identity levels
  int window = 1000;         // window (bases) of the heterogeneity test
  int het_mode = 1;          // 0 = pooled, 1 = auto (test -> posterior mean), 2 = always posterior mean
  double het_z = 3.0;        // z threshold of the heterogeneity test
  double lr_min = 10.0;      // log-likelihood ratio (vs. all-unaligned) needed to report
  double ident_c = 8.0;      // keep aligned levels with theta >= ident_c * min effective FPR
  double grid_step = 0.005;  // ANI spacing of aligned levels
  double floor_frac = 2.0;   // flag below_floor if this much aligned mass is on the 2 lowest levels (>1 = off)
  int bin = 32;              // probes per HMM step (main pass)
  int bin_het = 8;           // probes per HMM step when re-run for the posterior-mean ANI
  bool refine = true;        // exact per-probe posteriors around alignment boundaries
  double refine_lo = 0.05;   // refine bins with posterior in (lo, 1-lo) and 0.5-crossings (lo >= 0.5: crossings only)
  int fast_win = 32;         // probes per window of the all-aligned fast-path check (0 = off)
  double fast_llr = 10.0;    // every window must favour an aligned level over U by this (nats)
  bool bias_correct = true;
};

struct Result {
  double ani = NAN, af = 0, theta = 0, lr = 0;
  size_t probes = 0, hits = 0;
  double aligned_probes = 0;
  int windows = 0;
  bool het = false;
  double low_frac = 0;
  const char* status = "ok";
};

double Bcoef(double p, int k) {
  if (p <= 1e-9) return 0;
  double q = 1 - std::pow(1 - p, k);
  return std::max(0.0, q * (2.0 * k - 1 + 2.0 / p) - 2.0 * k);
}

struct alignas(64) Work {  // per-thread scratch buffers, reused across contigs
  Packed pk;
  std::vector<float> Ef, Ff, bf, wf, Eu, Fu, Es;
  std::vector<uint64_t> canon, hs, bh;
  std::vector<uint32_t> kpos, idx;
  std::vector<const uint64_t*> blk;
  std::vector<Probe> pr;
  std::vector<uint32_t> bstart, bfirst, blast;
  std::vector<double> ab;  // per-bin aligned posterior
  std::vector<float> Bb, Bbu;           // stored backward messages per bin
  std::vector<uint8_t> refined;         // per bin: per-probe posteriors available in ap
  std::vector<double> ap;               // per-probe aligned posterior (refined bins only)
  std::vector<float> rfa, rfu, rb, rw;  // refinement scratch
  std::vector<uint8_t> bcnt;  // per bin: hits and misses per tier
  std::vector<double> E, F, a, b, w, conf;
};

struct Estimator {
  const Sketch& sk;
  Params P;
  int T;
  std::vector<double> epsp, stier, gth, gani;
  size_t SPA = 16;                 // aligned states padded to a multiple of 16
  std::vector<float> PHa, PMa;     // [t][c][aligned state] = (hit prob)^c, (miss prob)^c
  std::vector<double> PHu, PMu;    // [t][c] for the unaligned state
  std::vector<float> gaf;          // padded aligned-level ANIs
  std::vector<float> LH, LM;       // [t][aligned state]: log((hit|A)/(hit|U)), log((miss|A)/(miss|U))
  std::vector<float> PE;           // per-probe emissions [t*2+hit][aligned state]
  float PEu[2 * MAXT];             // per-probe emission of U
  float PEs[2 * MAXT];             // row sums of PE
  double logU[2 * MAXT];
  static const size_t QMAX = 1 << 14;
  std::vector<double> tab_su, tab_h;

  Estimator(const Sketch& s, const Params& p) : sk(s), P(p), T((int)s.tiers.size()) {
    double sg = sk.sigma_used();
    for (auto& t : sk.tiers) {
      epsp.push_back(t.eps + sg - t.eps * sg);
      stier.push_back(t.s());
    }
    double emin = 1;
    for (double e : epsp) emin = std::min(emin, e);
    double amin = std::max(std::pow(std::min(0.9, P.ident_c * emin), 1.0 / sk.k), 0.70);
    for (double a = 0.9995; a >= amin - 1e-12; a -= P.grid_step) {
      gani.push_back(a);
      gth.push_back(std::min(std::max(std::pow(a, sk.k), 1e-7), 0.9995));
    }
    const size_t m = gth.size(), B = (size_t)std::max(1, std::max(P.bin, P.bin_het));
    SPA = (m + 15) / 16 * 16;
    PHa.assign((size_t)T * (B + 1) * SPA, 0.0f);
    PMa.assign((size_t)T * (B + 1) * SPA, 0.0f);
    PHu.assign((size_t)T * (B + 1), 1.0);
    PMu.assign((size_t)T * (B + 1), 1.0);
    gaf.assign(SPA, 0.0f);
    for (size_t j = 0; j < m; ++j) gaf[j] = (float)gani[j];
    PE.assign((size_t)T * 2 * SPA, 0.0f);
    for (int t = 0; t < T; ++t) {
      double e = epsp[t];
      PEu[2 * t] = (float)(1 - e);
      PEu[2 * t + 1] = (float)e;
      PEs[2 * t] = PEs[2 * t + 1] = 0.0f;
      for (size_t j = 0; j < m; ++j) {
        double ph = gth[j] + (1 - gth[j]) * e;
        PE[((size_t)t * 2 + 0) * SPA + j] = (float)(1 - ph);
        PE[((size_t)t * 2 + 1) * SPA + j] = (float)ph;
        PEs[2 * t] += (float)(1 - ph);
        PEs[2 * t + 1] += (float)ph;
      }
    }
    LH.assign((size_t)T * SPA, -1e30f);
    LM.assign((size_t)T * SPA, -1e30f);
    for (int t = 0; t < T; ++t) {
      double e = epsp[t];
      for (size_t j = 0; j < m; ++j) {
        double ph = gth[j] + (1 - gth[j]) * e, hp = 1, mp = 1;
        LH[(size_t)t * SPA + j] = (float)std::log(ph / e);
        LM[(size_t)t * SPA + j] = (float)std::log((1 - ph) / (1 - e));
        for (size_t c = 0; c <= B; ++c) {
          PHa[((size_t)t * (B + 1) + c) * SPA + j] = (float)hp;
          PMa[((size_t)t * (B + 1) + c) * SPA + j] = (float)mp;
          hp *= ph;
          mp *= 1 - ph;
        }
      }
      double hp = 1, mp = 1;
      for (size_t c = 0; c <= B; ++c) {
        PHu[(size_t)t * (B + 1) + c] = hp;
        PMu[(size_t)t * (B + 1) + c] = mp;
        hp *= e;
        mp *= 1 - e;
      }
      logU[2 * t] = std::log(1 - e);
      logU[2 * t + 1] = std::log(e);
    }
    tab_su.resize(QMAX);
    tab_h.resize(QMAX);
    for (size_t d = 0; d < QMAX; ++d) {
      tab_su[d] = std::exp(-P.rho_u * (double)d);
      tab_h[d] = -std::expm1(-P.rho_h * (double)d);
    }
  }
  inline double su(uint32_t d) const { return d < QMAX ? tab_su[d] : std::exp(-P.rho_u * d); }
  inline double hh(uint32_t d) const { return d < QMAX ? tab_h[d] : -std::expm1(-P.rho_h * d); }

  double combine(const double* N, const double* Y, double theta_ref, double* var_out = nullptr) const {
    double wsum = 0, acc = 0;
    double th = std::min(std::max(theta_ref, 1e-4), 1.0);
    for (int t = 0; t < T; ++t) {
      if (N[t] < 1e-9) continue;
      double e = epsp[t];
      double c = (Y[t] / N[t] - e) / (1 - e);
      double v = (th * (1 - th) * std::max(1 - stier[t], 0.0) + (1 - th) * e / (1 - e)) / N[t];
      v = std::max(v, 1e-12);
      acc += c / v;
      wsum += 1 / v;
    }
    if (wsum == 0) return 0;
    if (var_out) *var_out = 1 / wsum;
    return acc / wsum;
  }

  double theta_range(const std::vector<Probe>& pr, const double* w, double theta_ref, size_t lo, size_t hi,
                     double* var = nullptr, double* neff = nullptr) const {
    double N[MAXT] = {0}, Y[MAXT] = {0};
    double ne = 0;
    for (size_t i = lo; i < hi; ++i) {
      N[pr[i].tier] += w[i];
      Y[pr[i].tier] += w[i] * pr[i].hit;
      ne += w[i];
    }
    if (neff) *neff = ne;
    return combine(N, Y, theta_ref, var);
  }

  double ani_from_theta(double th, double var_theta, double nkmers) const {
    int k = sk.k;
    th = std::min(std::max(th, 1e-12), 1.0);
    double ani = std::pow(th, 1.0 / k);
    if (P.bias_correct) {
      double p = 1 - ani;
      double vproc = nkmers > 0 ? th * Bcoef(p, k) / nkmers : 0;
      double v = var_theta + vproc;
      ani *= 1 + (k - 1.0) / (2.0 * k * k) * std::min(v / (th * th), 1.0);
    }
    return std::min(ani, 1.0);
  }

  // probe all selected k-mers of seq into W.pr
  void probe(const std::string& seq, Work& W) const {
    const size_t L = seq.size();
    const int k = sk.k;
    W.pr.clear();
    if (L < (size_t)k) return;
    const size_t nst = L - k + 1;
    if (W.canon.size() < nst + 8) {
      size_t cap = nst + 8;
      W.canon.resize(cap);
      W.kpos.resize(cap);
      W.hs.resize(cap);
      W.blk.resize(cap);
      W.bh.resize(cap);
    }
    pack_seq((const uint8_t*)seq.data(), L, k, W.pk);
#ifdef AIF_AVX512
    if (T <= 2) {
      if (W.kpos.size() < nst + 32) W.kpos.resize(nst + 32);
      const size_t m = select_positions(W.pk, k, sk.tiers.back().hi, W.kpos.data());
      W.pr.resize(m);
      const long long* wbase64 = (const long long*)W.pk.w.data();
      const uint64_t kmask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
      alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                                  23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                                  39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                                  55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
      const __m512i bsw = _mm512_load_si512((const void*)BSW);
      const __m512i vmask = _mm512_set1_epi64((long long)kmask), vcomp = _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL);
      const __m512i m33 = _mm512_set1_epi64((long long)0x3333333333333333ULL), m0f = _mm512_set1_epi64((long long)0x0F0F0F0F0F0F0F0FULL);
      const __m512i vsel = _mm512_set1_epi64((long long)SEED_SEL), v63 = _mm512_set1_epi64(63), v31 = _mm512_set1_epi64(31);
      const __m128i rcs = _mm_cvtsi32_si128(64 - 2 * k);
      const Tier& t0 = sk.tiers[0];
      const Tier& t1 = sk.tiers[T - 1];
      const long long* base0 = (const long long*)t0.bits.data();
      const long long delta1 = (long long)((const long long*)t1.bits.data() - base0);
      const __m512i vhi0 = _mm512_set1_epi64((long long)t0.hi);
      const __m512i vnb0 = _mm512_set1_epi64((long long)t0.nblocks), vnb1 = _mm512_set1_epi64((long long)t1.nblocks);
      const __m512i voff1 = _mm512_set1_epi64(delta1), zero = _mm512_setzero_si512();
      const __m512i sbf = _mm512_set1_epi64((long long)SEED_BF), sbit = _mm512_set1_epi64((long long)SEED_BIT);
      const __m512i m511 = _mm512_set1_epi64(511), m63 = _mm512_set1_epi64(63), one = _mm512_set1_epi64(1);
      const uint32_t nh0 = t0.nh, nh1 = t1.nh, nhmax = std::max(nh0, nh1);
      for (size_t j = 0; j < m; j += 8) {
        const __mmask8 tail = (m - j >= 8) ? (__mmask8)0xFF : (__mmask8)((1u << (m - j)) - 1);
        const __m512i pos = _mm512_cvtepu32_epi64(_mm256_maskz_loadu_epi32(tail, W.kpos.data() + j));
        const __m512i wi = _mm512_srli_epi64(pos, 5);
        const __m512i A = _mm512_mask_i64gather_epi64(zero, tail, wi, wbase64, 8);
        const __m512i Bw = _mm512_slli_epi64(_mm512_mask_i64gather_epi64(zero, tail, _mm512_add_epi64(wi, one), wbase64, 8), 1);
        const __m512i sh = _mm512_slli_epi64(_mm512_and_si512(pos, v31), 1);
        __m512i f = _mm512_or_si512(_mm512_srlv_epi64(A, sh), _mm512_sllv_epi64(Bw, _mm512_sub_epi64(v63, sh)));
        f = _mm512_and_si512(f, vmask);
        __m512i r = _mm512_xor_si512(f, vcomp);
        r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 2), m33), _mm512_slli_epi64(_mm512_and_si512(r, m33), 2));
        r = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(r, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(r, m0f), 4));
        r = _mm512_srl_epi64(_mm512_shuffle_epi8(r, bsw), rcs);
        const __m512i x = _mm512_min_epu64(f, r);
        const __m512i h = fmix64_v(_mm512_xor_si512(x, vsel));
        const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_cmpge_epu64_mask(h, vhi0) & tail) : (__mmask8)0;
        const __m512i bh = fmix64_v(_mm512_xor_si512(x, sbf));
        const __m512i nb = _mm512_mask_blend_epi64(m1, vnb0, vnb1);
        const __m512i blk = _mm512_srli_epi64(_mm512_mul_epu32(_mm512_srli_epi64(bh, 32), nb), 32);
        const __m512i wb = _mm512_add_epi64(_mm512_slli_epi64(blk, 3), _mm512_mask_blend_epi64(m1, zero, voff1));
        __m512i h2 = fmix64_v(_mm512_xor_si512(bh, sbit));
        __m512i acc = one;
        for (uint32_t q = 0; q < nhmax; ++q) {  // independent gathers (no early exit: it would serialise them)
          if (q && q % 7 == 0) h2 = fmix64_v(_mm512_add_epi64(h2, _mm512_set1_epi64(q)));
          const __m512i bp = _mm512_and_si512(_mm512_srl_epi64(h2, _mm_cvtsi32_si128(9 * (q % 7))), m511);
          const __m512i idx = _mm512_add_epi64(wb, _mm512_srli_epi64(bp, 6));
          const __mmask8 am = (__mmask8)(((q < nh0) ? (__mmask8)(~m1) : 0) | ((q < nh1) ? m1 : 0)) & tail;
          const __m512i g = _mm512_mask_i64gather_epi64(zero, am, idx, base0, 8);
          acc = _mm512_mask_and_epi64(acc, am, acc, _mm512_srlv_epi64(g, _mm512_and_si512(bp, m63)));
        }
        const __mmask8 hit = _mm512_test_epi64_mask(acc, one) & tail;
        const size_t lim = std::min<size_t>(8, m - j);
        for (size_t l = 0; l < lim; ++l)
          W.pr[j + l] = {W.kpos[j + l], (uint8_t)((m1 >> l) & 1), (uint8_t)((hit >> l) & 1)};
      }
      return;
    }
#endif
    const size_t m =
        select_kmers(W.pk, k, sk.tiers.back().hi, 0, nst, W.canon.data(), W.hs.data(), W.kpos.data());
    W.pr.resize(m);
    const uint64_t* hs = W.hs.data();
    for (size_t j = 0; j < m; ++j) {
      int t = sk.tier_of(hs[j]);
      W.blk[j] = sk.tiers[t].block_of(W.canon[j], W.bh[j]);
      __builtin_prefetch(W.blk[j]);
      W.pr[j] = {W.kpos[j], (uint8_t)t, 0};
    }
    for (size_t j = 0; j < m; ++j) W.pr[j].hit = sk.tiers[W.pr[j].tier].contains_at(W.blk[j], W.bh[j]);
  }

  // binned forward-backward over states {U, A_1..A_m} in single precision; aligned states
  // are padded to a multiple of 16 (zero emission). Fills W.a (per-probe aligned posterior).
  // exact per-probe forward-backward over bins with uncertain posterior (alignment boundaries),
  // conditioned on the binned forward/backward messages at the region edges
  void refine(Work& W) const {
    const size_t nb = W.bstart.size(), SP = SPA, m = gth.size();
    const float pi = 1.0f / (float)m;
    const std::vector<Probe>& pr = W.pr;
    W.refined.assign(nb, 0);
    for (size_t bi = 0; bi < nb; ++bi) {
      const double a = W.ab[bi];
      bool unc = a > P.refine_lo && a < 1 - P.refine_lo;
      if (bi + 1 < nb && (a - 0.5) * (W.ab[bi + 1] - 0.5) < 0) unc = true;
      if (unc)
        for (size_t q = (bi ? bi - 1 : 0); q <= std::min(nb - 1, bi + 1); ++q) W.refined[q] = 1;
    }
    W.ap.resize(pr.size());
    size_t bi = 0;
    while (bi < nb) {
      if (!W.refined[bi]) {
        ++bi;
        continue;
      }
      size_t b0 = bi;
      while (bi < nb && W.refined[bi]) ++bi;
      const size_t b1 = bi - 1;
      const size_t i0 = W.bfirst[b0], i1 = W.bfirst[b1 + 1], nr = i1 - i0;
      W.rfa.resize(nr * SP);
      W.rfu.resize(nr);
      W.rb.resize(SP);
      W.rw.resize(SP);
      // forward over the region's probes
      float pu, pa_sum;
      const float* pa = nullptr;
      uint32_t ppos = 0;
      for (size_t r = 0; r < nr; ++r) {
        const Probe& p = pr[i0 + r];
        const float* e = &PE[((size_t)p.tier * 2 + p.hit) * SP];
        const float eu = PEu[2 * p.tier + p.hit];
        float* fc = &W.rfa[r * SP];
        float fu, iz;
        if (r == 0 && b0 == 0) {  // contig start: initial distribution
          fu = 0.5f * eu;
          float z = fu + 0.5f * pi * PEs[2 * p.tier + p.hit];
          iz = 1.0f / z;
          const float c0 = 0.5f * pi * iz;
          for (size_t j = 0; j < SP; ++j) fc[j] = c0 * e[j];
        } else {
          if (r == 0) {
            pa = &W.Ff[(b0 - 1) * SP];
            pu = W.Fu[b0 - 1];
            ppos = W.bstart[b0 - 1];
          }
          const uint32_t d = p.pos - ppos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
          pa_sum = std::max(0.0f, 1.0f - pu);
          const float base = pu * u * pi + stay * h * pi * pa_sum, c1 = stay * (1.0f - h);
          fu = (pu * stay + pa_sum * u) * eu;
          float dot = 0;
          for (size_t j = 0; j < SP; ++j) dot += pa[j] * e[j];
          const float z = fu + base * PEs[2 * p.tier + p.hit] + c1 * dot;
          iz = 1.0f / z;
          const float bz = base * iz, cz = c1 * iz;
          for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * pa[j]) * e[j];
        }
        W.rfu[r] = fu * iz;
        pa = fc;
        pu = W.rfu[r];
        ppos = p.pos;
      }
      // backward message entering the region's last probe
      float* bb = W.rb.data();
      float* ww = W.rw.data();
      float bu = 1.0f;
      if (b1 + 1 < nb) {
        const float* e = &W.Ef[(b1 + 1) * SP];
        const float* bn = &W.Bb[(b1 + 1) * SP];
        const float wu = W.Eu[b1 + 1] * W.Bbu[b1 + 1];
        float sw = 0;
        for (size_t j = 0; j < SP; ++j) {
          ww[j] = e[j] * bn[j];
          sw += ww[j];
        }
        const uint32_t d = W.bstart[b1 + 1] - pr[i1 - 1].pos;
        const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * pi;
        const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
        const float zz = b0v + (float)m * base + c1 * sw, iz = 1.0f / zz;
        for (size_t j = 0; j < SP; ++j) bb[j] = (base + c1 * ww[j]) * iz;
        bu = b0v * iz;
      } else {
        for (size_t j = 0; j < SP; ++j) bb[j] = 1.0f;
      }
      // backward over the region with per-probe posteriors
      for (size_t r = nr; r-- > 0;) {
        const float* fa = &W.rfa[r * SP];
        float tot = W.rfu[r] * bu, gu = tot;
        for (size_t j = 0; j < SP; ++j) tot += fa[j] * bb[j];
        W.ap[i0 + r] = 1.0 - (double)gu / (double)tot;
        if (r == 0) break;
        const Probe& p = pr[i0 + r];
        const float* e = &PE[((size_t)p.tier * 2 + p.hit) * SP];
        const float wu = PEu[2 * p.tier + p.hit] * bu;
        float sw = 0;
        for (size_t j = 0; j < SP; ++j) {
          ww[j] = e[j] * bb[j];
          sw += ww[j];
        }
        const uint32_t d = p.pos - pr[i0 + r - 1].pos;
        const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * pi;
        const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
        const float zz = b0v + (float)m * base + c1 * sw, iz = 1.0f / zz;
        for (size_t j = 0; j < SP; ++j) bb[j] = (base + c1 * ww[j]) * iz;
        bu = b0v * iz;
      }
      // per-bin mean posterior for bin-level uses (heterogeneity windows)
      for (size_t q = b0; q <= b1; ++q) {
        double acc = 0;
        for (uint32_t i = W.bfirst[q]; i < W.bfirst[q + 1]; ++i) acc += W.ap[i];
        W.ab[q] = acc / (double)(W.bfirst[q + 1] - W.bfirst[q]);
      }
    }
  }

  // bins of B consecutive probes: start/last position, per-tier hit/miss counts; returns llU
  double make_bins(Work& W, size_t& hits, int binsz) const {
    const std::vector<Probe>& pr = W.pr;
    const size_t n = pr.size(), B = (size_t)std::max(1, binsz), nb = (n + B - 1) / B;
    W.bstart.resize(nb);
    W.blast.resize(nb);
    W.bfirst.resize(nb + 1);
    W.bcnt.assign(nb * 2 * T, 0);
    double llU = 0;
    hits = 0;
    for (size_t bi = 0; bi < nb; ++bi) {
      size_t f = bi * B, l = std::min(n, f + B);
      W.bstart[bi] = pr[f].pos;
      W.blast[bi] = pr[l - 1].pos;
      W.bfirst[bi] = (uint32_t)f;
      uint8_t* c = &W.bcnt[bi * 2 * T];
      for (size_t i = f; i < l; ++i) ++c[2 * pr[i].tier + pr[i].hit];
      for (int t = 0; t < T; ++t) {
        llU += c[2 * t] * logU[2 * t] + c[2 * t + 1] * logU[2 * t + 1];
        hits += c[2 * t + 1];
      }
    }
    W.bfirst[nb] = (uint32_t)n;
    return llU;
  }

  void hmm(Work& W, double& ll, double& pm_ani, double& low_frac) const {
    const size_t m = gth.size(), SP = SPA;
    const size_t B = (size_t)std::max(1, std::max(P.bin, P.bin_het));  // power-table stride
    const size_t nb = W.bstart.size();
    const float pi = 1.0f / (float)m;
    // emissions
    W.Ef.resize(nb * SP);
    W.Eu.resize(nb);
    W.Es.resize(nb);
    for (size_t bi = 0; bi < nb; ++bi) {
      float* __restrict e = &W.Ef[bi * SP];
      const uint8_t* c = &W.bcnt[bi * 2 * T];
      const float* h0 = &PHa[((size_t)0 * (B + 1) + c[1]) * SP];
      const float* m0 = &PMa[((size_t)0 * (B + 1) + c[0]) * SP];
      for (size_t j = 0; j < SP; ++j) e[j] = h0[j] * m0[j];
      double eu = PHu[c[1]] * PMu[c[0]];
      for (int t = 1; t < T; ++t) {
        const float* ht = &PHa[((size_t)t * (B + 1) + c[2 * t + 1]) * SP];
        const float* mt = &PMa[((size_t)t * (B + 1) + c[2 * t]) * SP];
        for (size_t j = 0; j < SP; ++j) e[j] *= ht[j] * mt[j];
        eu *= PHu[t * (B + 1) + c[2 * t + 1]] * PMu[t * (B + 1) + c[2 * t]];
      }
      W.Eu[bi] = (float)eu;
      float se = 0;
      for (size_t j = 0; j < SP; ++j) se += e[j];
      W.Es[bi] = se;
    }
    // forward
    W.Ff.resize(nb * SP);
    W.Fu.resize(nb);
    ll = 0;
    double zprod = 1.0;
    int zcnt = 0;
    {
      const float* e = &W.Ef[0];
      float* f = &W.Ff[0];
      float fu = 0.5f * W.Eu[0], z = fu;
      for (size_t j = 0; j < SP; ++j) {
        f[j] = 0.5f * pi * e[j];
        z += f[j];
      }
      const float iz = 1.0f / z;
      for (size_t j = 0; j < SP; ++j) f[j] *= iz;
      W.Fu[0] = fu * iz;
      ll += std::log((double)z);
    }
    for (size_t bi = 1; bi < nb; ++bi) {
      uint32_t d = W.bstart[bi] - W.bstart[bi - 1];
      const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
      const float* __restrict fp = &W.Ff[(bi - 1) * SP];
      float* __restrict fc = &W.Ff[bi * SP];
      const float* __restrict e = &W.Ef[bi * SP];
      const float fpu = W.Fu[bi - 1];
      const float SA = std::max(0.0f, 1.0f - fpu);  // previous step is normalised
      const float base = fpu * u * pi + stay * h * pi * SA;
      const float c1 = stay * (1.0f - h);
      const float fu = (fpu * stay + SA * u) * W.Eu[bi];
      float dot = 0;
      for (size_t j = 0; j < SP; ++j) dot += fp[j] * e[j];
      const float z = fu + base * W.Es[bi] + c1 * dot;
      const float iz = 1.0f / z, bz = base * iz, cz = c1 * iz;
      for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * fp[j]) * e[j];
      W.Fu[bi] = fu * iz;
      zprod *= (double)z;
      if (++zcnt == 16) {
        ll += std::log(zprod);
        zprod = 1.0;
        zcnt = 0;
      }
    }
    ll += std::log(zprod);
    // backward with posteriors
    W.ab.resize(nb);
    W.bf.assign(SP, 1.0f);
    W.wf.resize(SP);
    float bu = 1.0f;
    double pm_num = 0, pm_den = 0, low = 0;
    const size_t nlow = std::min<size_t>(2, m);
    auto post = [&](size_t bi, const float* fa, float fu, const float* ba, float bub) {
      float pu = fu * bub, tot = pu, num = 0, lw = 0;
      for (size_t j = 0; j < SP; ++j) {
        float g = fa[j] * ba[j];
        tot += g;
        num += g * gaf[j];
      }
      for (size_t j = m - nlow; j < m; ++j) lw += fa[j] * ba[j];
      double it = 1.0 / (double)tot, al = 1.0 - (double)pu * it;
      double cnt = (double)(W.bfirst[bi + 1] - W.bfirst[bi]);
      W.ab[bi] = al;
      pm_num += cnt * (double)num * it;
      pm_den += cnt * al;
      low += cnt * (double)lw * it;
    };
    W.Bb.resize(nb * SP);
    W.Bbu.resize(nb);
    std::fill(W.Bb.begin() + (nb - 1) * SP, W.Bb.begin() + nb * SP, 1.0f);
    W.Bbu[nb - 1] = 1.0f;
    post(nb - 1, &W.Ff[(nb - 1) * SP], W.Fu[nb - 1], W.bf.data(), bu);
    for (size_t bi = nb - 1; bi > 0; --bi) {
      uint32_t d = W.bstart[bi] - W.bstart[bi - 1];
      const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
      const float* __restrict e = &W.Ef[bi * SP];
      float* __restrict bb = W.bf.data();
      float* __restrict ww = W.wf.data();
      const float wu = W.Eu[bi] * bu;
      float sw = 0;
      for (size_t j = 0; j < SP; ++j) {
        ww[j] = e[j] * bb[j];
        sw += ww[j];
      }
      const float SW = sw * pi;
      const float base = u * wu + stay * h * SW;
      const float c1 = stay * (1.0f - h);
      const float b0 = stay * wu + u * SW;
      const float zz = b0 + (float)m * base + c1 * sw;  // any positive scale works
      const float iz = 1.0f / zz, bz = base * iz, cz = c1 * iz;
      for (size_t j = 0; j < SP; ++j) bb[j] = bz + cz * ww[j];
      bu = b0 * iz;
      std::memcpy(&W.Bb[(bi - 1) * SP], bb, SP * sizeof(float));
      W.Bbu[bi - 1] = bu;
      post(bi - 1, &W.Ff[(bi - 1) * SP], W.Fu[bi - 1], bb, bu);
    }
    pm_ani = pm_den > 0 ? pm_num / pm_den : NAN;
    low_frac = pm_den > 0 ? low / pm_den : 0;
  }

  // true if every window of ~fast_win probes (whole bins) favours some aligned level over U
  // by fast_llr nats; returns the summed log-likelihood ratio in llr_out
  bool all_aligned(const Work& W, double& llr_out) const {
    const size_t nb = W.bstart.size(), SP = SPA;
    const size_t wb = std::max<size_t>(1, (size_t)P.fast_win / (size_t)std::max(1, P.bin));
    if (P.fast_win == 0 || nb < wb || SP > 256) return false;
    llr_out = 0;
    alignas(64) float acc[256];
    for (size_t lo = 0; lo < nb; lo += wb) {
      size_t hi = std::min(nb, lo + wb);
      if (hi - lo < wb && lo > 0) lo = nb - wb, hi = nb;  // short tail: re-check the last full window
      uint32_t cnt[2 * MAXT] = {0};
      for (size_t bi = lo; bi < hi; ++bi)
        for (int q = 0; q < 2 * T; ++q) cnt[q] += W.bcnt[bi * 2 * T + q];
      for (size_t j = 0; j < SP; ++j) acc[j] = 0;
      for (int t = 0; t < T; ++t) {
        const float ch = (float)cnt[2 * t + 1], cm = (float)cnt[2 * t];
        const float* lh = &LH[(size_t)t * SP];
        const float* lm = &LM[(size_t)t * SP];
        for (size_t j = 0; j < SP; ++j) acc[j] += ch * lh[j] + cm * lm[j];
      }
      float best = -1e30f;
      for (size_t j = 0; j < SP; ++j) best = std::max(best, acc[j]);
      if (best < (float)P.fast_llr) return false;
      llr_out += best;
      if (hi == nb) break;
    }
    return true;
  }

  Result estimate(const std::string& seq, Work& W) const {
    Result R;
    const int k = sk.k;
    probe(seq, W);
    const std::vector<Probe>& pr = W.pr;
    R.probes = pr.size();
    if (pr.size() < 3) {
      for (auto& p : pr) R.hits += p.hit;
      R.status = "too_few_probes";
      return R;
    }
    const double llU = make_bins(W, R.hits, P.bin);
    const size_t nb = W.bstart.size();
    double ll = 0, pm_ani = NAN, low_frac = 0, fast_lr = 0;
    const bool fast = all_aligned(W, fast_lr);
    if (fast) {  // confidently aligned everywhere: the posterior is ~1 in every bin
      W.ab.assign(nb, 1.0);
      R.lr = fast_lr;
    } else {
      hmm(W, ll, pm_ani, low_frac);
      R.lr = ll - llU;
      if (P.refine && P.bin > 1) refine(W);
      else W.refined.assign(nb, 0);
    }
    if (fast) W.refined.assign(nb, 0);
    if (R.lr < P.lr_min) {
      R.status = "no_homology";
      return R;
    }
    R.low_frac = low_frac;
    const double* ab = W.ab.data();
    // pass 1 over bins: posterior-weighted tier counts, aligned-start integral, segments
    double N[MAXT] = {0}, Y[MAXT] = {0}, neff = 0;
    const size_t L = seq.size();
    const size_t nstarts = L - k + 1;
    // per-bin posterior at its first/last probe and integral over its span (exact per-probe
    // trapezoids in refined bins, constant otherwise)
    auto afirst = [&](size_t b) { return W.refined[b] ? W.ap[W.bfirst[b]] : ab[b]; };
    auto alast = [&](size_t b) { return W.refined[b] ? W.ap[W.bfirst[b + 1] - 1] : ab[b]; };
    double Sx = afirst(0) * (double)(W.bstart[0] + 1);
    int segs = 0;
    double prev = 0;
    for (size_t bi = 0; bi < nb; ++bi) {
      const uint8_t* c = &W.bcnt[bi * 2 * T];
      if (W.refined[bi]) {
        for (uint32_t i = W.bfirst[bi]; i < W.bfirst[bi + 1]; ++i) {
          const double a = W.ap[i];
          N[pr[i].tier] += a;
          Y[pr[i].tier] += a * pr[i].hit;
          neff += a;
          if (i > W.bfirst[bi]) Sx += 0.5 * (W.ap[i - 1] + a) * (double)(pr[i].pos - pr[i - 1].pos);
          segs += (a > 0.5) & (i == 0 || prev <= 0.5);
          prev = a;
        }
      } else {
        const double a = ab[bi];
        for (int t = 0; t < T; ++t) {
          N[t] += a * (c[2 * t] + c[2 * t + 1]);
          Y[t] += a * c[2 * t + 1];
        }
        neff += a * (double)(W.bfirst[bi + 1] - W.bfirst[bi]);
        Sx += a * (double)(W.blast[bi] - W.bstart[bi]);
        segs += (a > 0.5) & (bi == 0 || prev <= 0.5);
        prev = a;
      }
      if (bi + 1 < nb) Sx += 0.5 * (alast(bi) + afirst(bi + 1)) * (double)(W.bstart[bi + 1] - W.blast[bi]);
    }
    Sx += alast(nb - 1) * (double)(nstarts - 1 - W.blast[nb - 1]);
    double var = 0;
    const double th0 = combine(N, Y, 0.3);
    const double th = combine(N, Y, th0, &var);
    R.theta = th;
    R.aligned_probes = neff;
    R.af = std::min(1.0, std::max(0.0, (Sx + (double)(k - 1) * segs) / (double)L));
    const double nk_aligned = std::max(1.0, R.af * L - (k - 1));
    if (th <= 0) {
      R.status = "no_homology";
      return R;
    }
    const double ani_pooled = ani_from_theta(th, var, nk_aligned);
    R.ani = ani_pooled;
    if (P.het_mode == 0 || P.window <= 0) return R;
    // pass 2: heterogeneity test over windows (whole bins) of confidently aligned bins
    const double bproc = std::max(th, 1e-6) * Bcoef(1 - ani_pooled, k);
    double Q = 0;
    int nw = 0;
    size_t i = 0;
    while (i < nb) {
      size_t j = i;
      const uint32_t start = W.bstart[i];
      while (j < nb && W.bstart[j] < start + (uint32_t)P.window) ++j;
      if (j < nb && W.blast[nb - 1] < W.bstart[j] + (uint32_t)P.window / 2) j = nb;
      double Nw[MAXT] = {0}, Yw[MAXT] = {0}, ne = 0, nall = 0;
      for (size_t q = i; q < j; ++q) {
        const uint8_t* c = &W.bcnt[q * 2 * T];
        const double cnt = (double)(W.bfirst[q + 1] - W.bfirst[q]);
        nall += cnt;
        if (ab[q] <= 0.9) continue;
        for (int t = 0; t < T; ++t) {
          Nw[t] += c[2 * t] + c[2 * t + 1];
          Yw[t] += c[2 * t + 1];
        }
        ne += cnt;
      }
      double v = 0;
      const double tw = combine(Nw, Yw, th, &v);
      const double len = ne / std::max(1.0, nall) * (double)(W.blast[j - 1] - start + 1);
      if (len > 100 && ne >= 10) {
        Q += (tw - th) * (tw - th) / std::max(v + bproc / len, 1e-12);
        ++nw;
      }
      i = j;
    }
    R.windows = nw;
    bool het = P.het_mode == 2;
    if (P.het_mode == 1 && nw >= 2) {
      const double df = (double)nw - 1;
      het = Q > df + P.het_z * std::sqrt(2 * df);
    }
    if (het && (fast || P.bin_het != P.bin)) {  // posterior-mean ANI from a finer-binned HMM
      size_t hh = 0;
      make_bins(W, hh, P.bin_het);
      hmm(W, ll, pm_ani, low_frac);
      R.low_frac = low_frac;
    }
    if (het && pm_ani == pm_ani) {
      R.het = true;
      R.ani = pm_ani;
    }
    if (R.low_frac > P.floor_frac) R.status = "below_floor";
    return R;
  }
};

double composition_sf2(const std::vector<Packed>& packs) {
  double cnt[4] = {0, 0, 0, 0};
  const uint64_t lo = 0x5555555555555555ULL;
  for (auto& pk : packs) {
    const size_t nw = (pk.L + 31) / 32;
    for (size_t i = 0; i < nw; ++i) {
      const uint64_t x = pk.w[i], hb = (x >> 1) & lo, lb = x & lo;
      cnt[1] += __builtin_popcountll(lb & ~hb);
      cnt[2] += __builtin_popcountll(hb & ~lb);
      cnt[3] += __builtin_popcountll(hb & lb);
    }
    cnt[0] += (double)pk.L - pk.bad.size();
  }
  cnt[0] -= cnt[1] + cnt[2] + cnt[3];
  double tot = cnt[0] + cnt[1] + cnt[2] + cnt[3], sf2 = 0;
  for (double c : cnt) sf2 += (c / tot) * (c / tot);
  return std::max(sf2, 0.25);
}

void usage() {
  fprintf(stderr,
          "aif sketch <ref.fa> -o <out.aif> [-k auto|N] [--kmin 17] [--kmax 21] [--sigma-target 0.003]\n"
          "           [--budget 0.88] [--tiers 0.005:0.6,0.1:0.4] [--max-hashes 6]\n"
          "aif query <ref.aif> <contigs.fa> [-t threads] [--bin 32] [--bin-het 8] [--fast-win 32] [--fast-llr 10] [--rho-u 3e-6] [--rho-h 1e-4] [--ident-c 8]\n"
          "          [--lr-min 10] [--het 0|1|2 (1)] [--window 1000] [--het-z 3] [--grid-step 0.005]\n"
          "          [--floor-frac 2] [--no-bias-correct]\n"
          "output: contig ani af theta probes hits aligned_probes status lr het low_frac\n"
          "aif info <ref.aif>\n");
  exit(1);
}

void format_result(std::string& out, const std::string& name, const Result& r) {
  char buf[512];
  int n = snprintf(buf, sizeof buf, "\t%.6f\t%.6f\t%.6f\t%zu\t%zu\t%.1f\t%s\t%.1f\t%d\t%.3f\n", r.ani, r.af, r.theta,
                   r.probes, r.hits, r.aligned_probes, r.status, r.lr, (int)r.het, r.low_frac);
  out.append(name);
  out.append(buf, (size_t)n);
}

int cmd_query(int argc, char** argv) {
  if (argc < 4) usage();
  Sketch sk;
  sk.load(argv[2]);
  Params P;
  int threads = 1;
  for (int i = 4; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-t") threads = std::max(1, atoi(argv[++i]));
    else if (a == "--bin") P.bin = std::max(1, std::min(64, atoi(argv[++i])));
    else if (a == "--bin-het") P.bin_het = std::max(1, std::min(64, atoi(argv[++i])));
    else if (a == "--fast-win") P.fast_win = std::max(0, atoi(argv[++i]));
    else if (a == "--no-refine") P.refine = false;
    else if (a == "--refine-lo") P.refine_lo = atof(argv[++i]);
    else if (a == "--fast-llr") P.fast_llr = atof(argv[++i]);
    else if (a == "--rho-u") P.rho_u = atof(argv[++i]);
    else if (a == "--rho-h") P.rho_h = atof(argv[++i]);
    else if (a == "--window") P.window = atoi(argv[++i]);
    else if (a == "--het") P.het_mode = atoi(argv[++i]);
    else if (a == "--het-z") P.het_z = atof(argv[++i]);
    else if (a == "--lr-min") P.lr_min = atof(argv[++i]);
    else if (a == "--ident-c") P.ident_c = atof(argv[++i]);
    else if (a == "--grid-step") P.grid_step = atof(argv[++i]);
    else if (a == "--floor-frac") P.floor_frac = atof(argv[++i]);
    else if (a == "--no-bias-correct") P.bias_correct = false;
    else usage();
  }
  Estimator est(sk, P);
#if defined(__x86_64__)
  _mm_setcsr(_mm_getcsr() | 0x8040);  // flush-to-zero / denormals-are-zero
#endif
  Fasta fa(argv[3]);
  fputs("contig\tani\taf\ttheta\tprobes\thits\taligned_probes\tstatus\tlr\thet\tlow_frac\n", stdout);
  std::string out;
  if (threads == 1) {
    Work W;
    std::string name, seq;
    while (fa.next(name, seq)) {
      format_result(out, name, est.estimate(seq, W));
      if (out.size() > (1 << 20)) {
        fwrite(out.data(), 1, out.size(), stdout);
        out.clear();
      }
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
  }
  // multi-threaded pipeline: the main thread reads batch i+1 while persistent workers
  // process batch i (double buffering); results are printed in input order
  struct Batch {
    std::vector<std::string> names, seqs;
    std::vector<Result> res;
    size_t n = 0;
  };
  Batch bufs[2];
  const size_t BATCH_BASES = (size_t)4 << 20;
  auto fill = [&](Batch& bt) {
    bt.n = 0;
    size_t bases = 0;
    while (bases < BATCH_BASES) {
      if (bt.n == bt.names.size()) {
        bt.names.emplace_back();
        bt.seqs.emplace_back();
      }
      if (!fa.next(bt.names[bt.n], bt.seqs[bt.n])) break;
      bases += bt.seqs[bt.n].size();
      ++bt.n;
    }
    if (bt.res.size() < bt.n) bt.res.resize(bt.n);
    return bt.n > 0;
  };
  std::mutex mu;
  std::condition_variable cv_work, cv_done;
  Batch* cur = nullptr;
  std::atomic<size_t> next_idx(0);
  int active = 0;
  uint64_t gen = 0;
  bool quit = false;
  std::vector<Work> works(threads);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&, t]() {
#if defined(__x86_64__)
      _mm_setcsr(_mm_getcsr() | 0x8040);
#endif
      uint64_t seen = 0;
      for (;;) {
        Batch* bt;
        {
          std::unique_lock<std::mutex> lk(mu);
          cv_work.wait(lk, [&] { return quit || gen != seen; });
          if (quit) return;
          seen = gen;
          bt = cur;
        }
        const size_t CHUNK = 16;  // contigs per grab (less contention / false sharing)
        for (size_t i0; (i0 = next_idx.fetch_add(CHUNK)) < bt->n;)
          for (size_t i = i0; i < std::min(bt->n, i0 + CHUNK); ++i) bt->res[i] = est.estimate(bt->seqs[i], works[t]);
        {
          std::lock_guard<std::mutex> lk(mu);
          if (--active == 0) cv_done.notify_one();
        }
      }
    });
  int which = 0;
  bool have = fill(bufs[which]);
  while (have) {
    {
      std::lock_guard<std::mutex> lk(mu);
      cur = &bufs[which];
      next_idx = 0;
      active = threads;
      ++gen;
    }
    cv_work.notify_all();
    const bool have_next = fill(bufs[1 - which]);  // read ahead while workers run
    {
      std::unique_lock<std::mutex> lk(mu);
      cv_done.wait(lk, [&] { return active == 0; });
    }
    out.clear();
    Batch& bt = bufs[which];
    for (size_t i = 0; i < bt.n; ++i) format_result(out, bt.names[i], bt.res[i]);
    fwrite(out.data(), 1, out.size(), stdout);
    which = 1 - which;
    have = have_next;
  }
  {
    std::lock_guard<std::mutex> lk(mu);
    quit = true;
  }
  cv_work.notify_all();
  for (auto& th : pool) th.join();
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  init_code();
  if (argc < 3) usage();
  std::string cmd = argv[1];
  try {
    if (cmd == "sketch") {
      const char* in = argv[2];
      std::string out;
      int k = 0, kmin = 17, kmax = 21;
      double sigma_target = 0.003, budget = 0.88;
      uint32_t max_hashes = 6;
      std::string tiers = "0.005:0.6,0.1:0.4";
      for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-o") out = argv[++i];
        else if (a == "-k") k = (std::string(argv[i + 1]) == "auto") ? (++i, 0) : atoi(argv[++i]);
        else if (a == "--kmin") kmin = atoi(argv[++i]);
        else if (a == "--kmax") kmax = atoi(argv[++i]);
        else if (a == "--sigma-target") sigma_target = atof(argv[++i]);
        else if (a == "--budget") budget = atof(argv[++i]);
        else if (a == "--tiers") tiers = argv[++i];
        else if (a == "--max-hashes") max_hashes = (uint32_t)std::max(1, atoi(argv[++i]));
        else usage();
      }
      if (out.empty() || (k != 0 && (k < 5 || k > 31))) usage();
      auto spec = parse_tiers(tiers);
      std::vector<Packed> packs;
      {
        Fasta fa(in);
        std::string name, seq;
        while (fa.next(name, seq)) {
          packs.emplace_back();
          pack_seq((const uint8_t*)seq.data(), seq.size(), k ? k : kmin, packs.back());
        }
      }
      Sketch sk;
      if (k != 0) {
        sk = build_sketch(packs, k, budget, spec, 1000000, max_hashes);
      } else {
        // adaptive k: smallest k >= kmin whose calibrated spurious rate is <= sigma_target;
        // start at the k predicted from composition, then verify
        double sf2 = composition_sf2(packs), npos = 0;
        for (auto& pk : packs) npos += (double)pk.L;
        int kk = kmin;
        while (kk < kmax && 4.0 * 2.0 * npos * std::pow(sf2, kk) > sigma_target) ++kk;
        sk = build_sketch(packs, kk, budget, spec, 1000000, max_hashes);
        while (sk.sigma_emp > sigma_target && kk < kmax) {
          int next = kk + 1;
          while (next < kmax && sk.sigma_emp * std::pow(sf2, next - kk) > sigma_target) ++next;
          kk = next;
          sk = Sketch();
          sk = build_sketch(packs, kk, budget, spec, 1000000, max_hashes);
        }
      }
      sk.save(out.c_str());
      fprintf(stderr, "k=%d bases=%lu bytes=%zu bits/base=%.4f sigma_analytic=%.3g sigma_emp=%.3g\n", sk.k,
              (unsigned long)sk.total_bases, sk.filter_bytes(), 8.0 * sk.filter_bytes() / sk.total_bases,
              sk.sigma(), sk.sigma_emp);
      for (auto& t : sk.tiers)
        fprintf(stderr, "  tier s=%.5f nh=%u n=%lu blocks=%lu eps=%.5f\n", t.s(), t.nh, (unsigned long)t.n_ins,
                (unsigned long)t.nblocks, t.eps);
    } else if (cmd == "query") {
      return cmd_query(argc, argv);
    } else if (cmd == "info") {
      Sketch sk;
      sk.load(argv[2]);
      printf("k=%d bases=%lu kmers=%lu budget=%.4f filter_bytes=%zu bits/base=%.4f sigma_analytic=%.3g sigma_emp=%.3g\n",
             sk.k, (unsigned long)sk.total_bases, (unsigned long)sk.n_kmers, sk.budget, sk.filter_bytes(),
             8.0 * sk.filter_bytes() / sk.total_bases, sk.sigma(), sk.sigma_emp);
      for (auto& t : sk.tiers)
        printf("tier s=%.5f nh=%u n=%lu blocks=%lu eps=%.5f\n", t.s(), t.nh, (unsigned long)t.n_ins,
               (unsigned long)t.nblocks, t.eps);
    } else {
      usage();
    }
  } catch (std::exception& e) {
    fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
