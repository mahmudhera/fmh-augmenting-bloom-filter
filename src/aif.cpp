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
// v1.2 (real genomes): long contigs are probed in 64 kb windows (cache-resident buffers,
// per-thread memory independent of contig length); threads use a streaming pipeline
// (chunks of contigs, ordered writer, bounded outstanding input) instead of batches with
// barriers; --het-fine-max N limits the fine-binned heterogeneity re-run to contigs with
// <= N probes (faster on Mbp-scale contigs, slightly less accurate genome-level ANI).
//
// v1.3 (single-thread speed): forward/backward recursions interleaved, with look-ahead
// normalisers that keep horizontal reductions off the dependency chain (HMM and
// refinement; rounding-level output changes, <= 1e-6); cached emissions for fine bins;
// GFNI reverse complements; AVX-512 2-bit packing and mask-popcount bin counts; batched,
// branch-free calibration. Sketches are byte-identical to v1.1/v1.2.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread -o aif aif.cpp
#include <algorithm>
#include <cstddef>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
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
  bool next_mapped(std::string& name, std::string& seq, bool append);
  // memory-mapped input only: the next record as a name and a raw span [b, e) of the map
  // (sequence lines with their line breaks); release_upto drops mapped pages before q
  bool mapped() const { return map != nullptr; }
  bool next_span(std::string& name, const char*& b, const char*& e, bool unused = false);
  void release_upto(const char* q);
  // next record; append = true appends the sequence to seq instead of replacing it
  bool next(std::string& name, std::string& seq, bool append = false) {
    if (map) return next_mapped(name, seq, append);
    if (!append) seq.clear();
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

// append [p, e) to seq without line breaks ('\n', '\r'); AVX-512 VBMI2 compresses 64
// bytes per step, otherwise one memchr per line
void append_stripped(std::string& seq, const char* p, const char* e) {
  const size_t n = (size_t)(e - p), at = seq.size();
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  seq.resize(at + n + 64);  // slack for full-width stores
  char* d = &seq[at];
  char* const d0 = d;
  const __m512i nl = _mm512_set1_epi8('\n'), cr = _mm512_set1_epi8('\r');
  size_t i = 0;
  for (; i + 64 <= n; i += 64) {
    const __m512i v = _mm512_loadu_si512((const void*)(p + i));
    const __mmask64 keep = _mm512_cmpneq_epi8_mask(v, nl) & _mm512_cmpneq_epi8_mask(v, cr);
    _mm512_storeu_si512((void*)d, _mm512_maskz_compress_epi8(keep, v));
    d += __builtin_popcountll(keep);
  }
  if (i < n) {
    const __mmask64 tail = (n - i == 64) ? ~0ULL : ((1ULL << (n - i)) - 1);
    const __m512i v = _mm512_maskz_loadu_epi8(tail, (const void*)(p + i));
    const __mmask64 keep = _mm512_cmpneq_epi8_mask(v, nl) & _mm512_cmpneq_epi8_mask(v, cr) & tail;
    _mm512_storeu_si512((void*)d, _mm512_maskz_compress_epi8(keep, v));
    d += __builtin_popcountll(keep);
  }
  seq.resize(at + (size_t)(d - d0));
#else
  seq.reserve(at + n);
  while (p < e) {
    const char* e2 = (const char*)memchr(p, '\n', (size_t)(e - p));
    const char* le = e2 ? e2 : e;
    size_t len = (size_t)(le - p);
    if (len && p[len - 1] == '\r') --len;
    seq.append(p, len);
    p = e2 ? e2 + 1 : e;
  }
#endif
}

bool Fasta::next_mapped(std::string& name, std::string& seq, bool append) {
  if (!append) seq.clear();
  const char *b, *e;
  if (!next_span(name, b, e, false)) return false;
  append_stripped(seq, b, e);
  // release consumed pages so the resident set stays small on large inputs
  const size_t pg = 1 << 22;
  if (mpos - released >= pg) release_upto(map + mpos);
  return true;
}

void Fasta::release_upto(const char* q) {
  const size_t upto = ((size_t)(q - map) / 4096) * 4096;
  if (upto > released) {
    madvise((void*)(map + released), upto - released, MADV_DONTNEED);
    released = upto;
  }
}

bool Fasta::next_span(std::string& name, const char*& b, const char*& e, bool) {
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
  b = p;
  e = nxt;
  mpos = (size_t)(nxt - map);
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
  size_t w0 = 0;
#if defined(__AVX512BW__)
  // 64 bases per step: code = (byte >> 1) & 3; pairs, quads and bytes are merged with
  // multiply-adds, and vpmovdb keeps one byte (4 codes) per dword -> two packed words
  {
    const __m512i lcbit = _mm512_set1_epi8(0x20), m3 = _mm512_set1_epi8(3);
    const __m512i ca = _mm512_set1_epi8('a'), cc = _mm512_set1_epi8('c'), cg = _mm512_set1_epi8('g'),
                  ct = _mm512_set1_epi8('t');
    const __m512i w14 = _mm512_set1_epi16(0x0401), w116 = _mm512_set1_epi32(0x00100001);
    for (; (w0 + 2) * 32 <= L; w0 += 2) {
      const size_t i0 = w0 * 32;
      const __m512i v = _mm512_loadu_si512((const void*)(s + i0));
      const __m512i lc = _mm512_or_si512(v, lcbit);
      const __mmask64 ok = _mm512_cmpeq_epi8_mask(lc, ca) | _mm512_cmpeq_epi8_mask(lc, cc) |
                           _mm512_cmpeq_epi8_mask(lc, cg) | _mm512_cmpeq_epi8_mask(lc, ct);
      if (__builtin_expect(ok != ~0ULL, 0))
        for (int b = 0; b < 64; ++b)
          if (CODE[s[i0 + b]] < 0) bad.push_back((uint32_t)(i0 + b));
      const __m512i c = _mm512_and_si512(_mm512_srli_epi16(v, 1), m3);
      const __m512i p2 = _mm512_maddubs_epi16(c, w14);   // c0 + 4 c1 per 16 bits
      const __m512i p4 = _mm512_madd_epi16(p2, w116);    // q0 + 16 q1 per 32 bits
      _mm_storeu_si128((__m128i*)&P.w[w0], _mm512_cvtepi32_epi8(p4));
    }
  }
#endif
  for (size_t w = w0; w < nw; ++w) {
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
// Complement every 2-bit base and reverse the 4 bases within each byte (the per-byte part
// of a reverse complement; a byte shuffle and a shift finish it). With GFNI this is one
// affine bit-matrix op: output bit pairs (7:6, 5:4, 3:2, 1:0) take input pairs (1:0, 3:2,
// 5:4, 7:6), and 0xAA flips each base's high bit (A<->T, C<->G).
inline __m512i comp_rev_bytes(__m512i x) {
#if defined(__GFNI__)
  return _mm512_gf2p8affine_epi64_epi8(x, _mm512_set1_epi64((long long)0x4080102004080102ULL), 0xAA);
#else
  x = _mm512_xor_si512(x, _mm512_set1_epi64((long long)0xAAAAAAAAAAAAAAAAULL));
  x = _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 2), m33), _mm512_slli_epi64(_mm512_and_si512(x, m33), 2));
  return _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi64(x, 4), m0f), _mm512_slli_epi64(_mm512_and_si512(x, m0f), 4));
#endif
}

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
      __m512i x = comp_rev_bytes(f);
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
      __m512i x = comp_rev_bytes(f);
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
  const __m512i vmask = _mm512_set1_epi64((long long)kmask);
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
    __m512i r = comp_rev_bytes(f);
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
      const __m512i vmask = _mm512_set1_epi64((long long)mask);
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
      // branch-free emission (no non-ACGT k-mers, <= 2 tiers): selected lanes are compressed
      // straight into cal_x/cal_t in lane order, exactly as the scalar emit() loop would
      const bool fast_emit = pk.badk.empty() && sk.tiers.size() <= 2;
      const __m512i vhi0 = _mm512_set1_epi64((long long)sk.tiers[0].hi);
      const size_t base_n = cal_x.size();
      size_t ne = 0;
      if (fast_emit) {
        cal_x.resize(base_n + 8 * G + 8);
        cal_t.resize(base_n + 8 * G + 8);
      }
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
        __m512i r = comp_rev_bytes(f);
        r = _mm512_srl_epi64(_mm512_shuffle_epi8(r, bsw), rcs);
        const __m512i x = _mm512_min_epu64(f, r);
        const __m512i h = fmix64_v(_mm512_xor_si512(x, vsel));
        __mmask8 sel = _mm512_cmplt_epu64_mask(h, vhmax);
        if (fast_emit) {
          const unsigned tm = sk.tiers.size() == 2 ? (unsigned)_mm512_mask_cmpge_epu64_mask(sel, h, vhi0) : 0u;
          _mm512_storeu_si512((void*)(cal_x.data() + base_n + ne), _mm512_maskz_compress_epi64(sel, x));
          const uint64_t tb = _pdep_u64(_pext_u32(tm, (unsigned)sel), 0x0101010101010101ULL);
          std::memcpy(cal_t.data() + base_n + ne, &tb, 8);
          ne += (size_t)__builtin_popcount((unsigned)sel);
          continue;
        }
        if (!sel) continue;
        _mm512_store_si512((void*)ox, x);
        _mm512_store_si512((void*)oh, h);
        _mm512_store_si512((void*)op, pos);
        for (; sel; sel &= (__mmask8)(sel - 1)) {
          int l = __builtin_ctz((unsigned)sel);
          emit(op[l], ox[l], oh[l]);
        }
      }
      if (fast_emit) {
        cal_x.resize(base_n + ne);
        cal_t.resize(base_n + ne);
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
  {  // batched membership tests: block addresses and prefetches first, then the bit tests
    const size_t BT = 64;
    const uint64_t* cb[BT];
    uint64_t ch[BT];
    for (size_t i0 = 0; i0 < cal_x.size(); i0 += BT) {
      const size_t n = std::min(BT, cal_x.size() - i0);
      for (size_t j = 0; j < n; ++j) {
        cb[j] = sk.tiers[cal_t[i0 + j]].block_of(cal_x[i0 + j], ch[j]);
        __builtin_prefetch(cb[j]);
      }
      for (size_t j = 0; j < n; ++j) {
        const int t = cal_t[i0 + j];
        N[t] += 1;
        Y[t] += sk.tiers[t].contains_at(cb[j], ch[j]);
      }
    }
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
  size_t het_fine_max = ~(size_t)0;  // re-run at bin_het only for contigs with <= this many probes (speed option)
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
  std::vector<const float*> ep;         // per-bin emission vectors (into Ef or the emission cache)
  std::vector<uint8_t> refined;         // per bin: per-probe posteriors available in ap
  std::vector<double> ap;               // per-probe aligned posterior (refined bins only)
  std::vector<float> rfa, rfu, rb, rw;  // refinement scratch
  std::vector<float> rbb, rbu;          // refinement: per-probe backward messages
  std::vector<uint8_t> bcnt;  // per bin: hits and misses per tier
  std::vector<double> E, F, a, b, w, conf;
};

#ifdef AIF_PROF  // per-stage cycle counters (make prof): probe, bins, hmm, refine, integrate, het-test, het-rerun
#include <x86intrin.h>
std::atomic<uint64_t> g_prof[10];
#define PSTART uint64_t _pt = __rdtsc()
#define PSTAMP(i)                       \
  do {                                  \
    const uint64_t _pn = __rdtsc();     \
    g_prof[i] += _pn - _pt;             \
    _pt = _pn;                          \
  } while (0)
#else
#define PSTART
#define PSTAMP(i)
#endif

void prof_report() {
#ifdef AIF_PROF
  static const char* nm[8] = {"probe", "bins", "hmm", "refine", "integrate", "het-test", "het-rerun", "-"};
  uint64_t tot = 0;
  for (int i = 0; i < 7; ++i) tot += g_prof[i];
  for (int i = 0; i < 7; ++i)
    fprintf(stderr, "prof %-10s %8.3f Gcyc %5.1f%%\n", nm[i], g_prof[i] / 1e9, 100.0 * g_prof[i] / std::max<uint64_t>(tot, 1));
  fprintf(stderr, "prof   of which pack %.3f Gcyc, select %.3f Gcyc\n", g_prof[7] / 1e9, g_prof[8] / 1e9);
  fprintf(stderr, "prof reader (FASTA parsing, threaded mode) %.3f Gcyc\n", g_prof[9] / 1e9);
#endif
}

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
  float PED[4 * MAXT * MAXT];      // <PE[a], PE[b]> for probe outcome pairs (refine look-ahead)
  double logU[2 * MAXT];
  static const size_t QMAX = 1 << 14;
  std::vector<double> tab_su, tab_h;
  // emission cache for fine (bin_het) bins: one vector per hit/miss count combination
  // (computed exactly as in hmm()); ecix maps the mixed-radix counts to a slot
  size_t EB = 0;                   // counts per tier and outcome are < EB + 1
  std::vector<int32_t> ecix;
  std::vector<float> ecE, ecU, ecS;
  inline size_t ec_key(const uint8_t* c) const {
    size_t key = 0;
    for (int q = 0; q < 2 * T; ++q) key = key * (EB + 1) + c[q];
    return key;
  }

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
    for (int a = 0; a < 2 * T; ++a)
      for (int b = 0; b < 2 * T; ++b) {
        float d = 0;
        for (size_t j = 0; j < SPA; ++j) d += PE[(size_t)a * SPA + j] * PE[(size_t)b * SPA + j];
        PED[a * 2 * T + b] = d;
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
    if (P.bin_het > 0 && P.bin_het < P.bin && T <= 2) build_emission_cache((size_t)P.bin_het, B);
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

  // probe all selected k-mers of seq into W.pr. Long contigs are probed in windows of
  // PROBE_WIN k-mer starts so the packed sequence and position buffers stay in L1/L2 and
  // the per-thread scratch does not grow with contig length.
  static constexpr size_t PROBE_WIN = (size_t)1 << 16;
  void probe(std::string_view seq, Work& W) const {
    const size_t L = seq.size();
    const int k = sk.k;
    W.pr.clear();
    if (L < (size_t)k) return;
    const size_t nst = L - k + 1;
    for (size_t off = 0; off < nst; off += PROBE_WIN) {
      const size_t n = std::min(PROBE_WIN, nst - off);
      probe_window((const uint8_t*)seq.data() + off, n + k - 1, (uint32_t)off, W);
    }
  }

  // probe the k-mers of s[0, L) and append them to W.pr with positions shifted by off
  void probe_window(const uint8_t* s, size_t L, uint32_t off, Work& W) const {
    const int k = sk.k;
    const size_t nst = L - k + 1;
    const size_t base = W.pr.size();
#ifdef AIF_PROF
    const uint64_t _p0 = __rdtsc();
#endif
    pack_seq(s, L, k, W.pk);
#ifdef AIF_PROF
    const uint64_t _p1 = __rdtsc();
    g_prof[7] += _p1 - _p0;  // pack (sub-stage of probe)
#endif
#ifdef AIF_AVX512
    if (T <= 2) {
      if (W.kpos.size() < nst + 32) W.kpos.resize(nst + 32);
      const size_t m = select_positions(W.pk, k, sk.tiers.back().hi, W.kpos.data());
#ifdef AIF_PROF
      g_prof[6 + 2] += __rdtsc() - _p1;  // select (sub-stage of probe)
#endif
      W.pr.resize(base + m);
      const long long* wbase64 = (const long long*)W.pk.w.data();
      const uint64_t kmask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
      alignas(64) static const uint8_t BSW[64] = {7,  6,  5,  4,  3,  2,  1,  0,  15, 14, 13, 12, 11, 10, 9,  8,
                                                  23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24,
                                                  39, 38, 37, 36, 35, 34, 33, 32, 47, 46, 45, 44, 43, 42, 41, 40,
                                                  55, 54, 53, 52, 51, 50, 49, 48, 63, 62, 61, 60, 59, 58, 57, 56};
      const __m512i bsw = _mm512_load_si512((const void*)BSW);
      const __m512i vmask = _mm512_set1_epi64((long long)kmask);
      const __m512i v63 = _mm512_set1_epi64(63), v31 = _mm512_set1_epi64(31), vsel = _mm512_set1_epi64((long long)SEED_SEL);
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
        __m512i r = comp_rev_bytes(f);
        r = _mm512_srl_epi64(_mm512_shuffle_epi8(r, bsw), rcs);
        const __m512i x = _mm512_min_epu64(f, r);
        const __m512i hs = fmix64_v(_mm512_xor_si512(x, vsel));  // selection hash: tier
        const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_cmpge_epu64_mask(hs, vhi0) & tail) : (__mmask8)0;
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
        Probe* out = W.pr.data() + base + j;
        for (size_t l = 0; l < lim; ++l)
          out[l] = {W.kpos[j + l] + off, (uint8_t)((m1 >> l) & 1), (uint8_t)((hit >> l) & 1)};
      }
      return;
    }
#endif
    if (W.canon.size() < nst + 8) {
      const size_t cap = nst + 8;
      W.canon.resize(cap);
      W.kpos.resize(cap);
      W.hs.resize(cap);
      W.blk.resize(cap);
      W.bh.resize(cap);
    }
    const size_t m =
        select_kmers(W.pk, k, sk.tiers.back().hi, 0, nst, W.canon.data(), W.hs.data(), W.kpos.data());
    W.pr.resize(base + m);
    Probe* out = W.pr.data() + base;
    const uint64_t* hs = W.hs.data();
    for (size_t j = 0; j < m; ++j) {
      int t = sk.tier_of(hs[j]);
      W.blk[j] = sk.tiers[t].block_of(W.canon[j], W.bh[j]);
      __builtin_prefetch(W.blk[j]);
      out[j] = {W.kpos[j] + off, (uint8_t)t, 0};
    }
    for (size_t j = 0; j < m; ++j) out[j].hit = sk.tiers[out[j].tier].contains_at(W.blk[j], W.bh[j]);
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
      W.rbb.resize(nr * SP);
      W.rbu.resize(nr);
      // Forward and backward recursions over the region's probes are independent: run them
      // interleaved, with look-ahead dots (as in hmm(); <e[r-1], e[r]> is a table lookup since
      // a probe's emission is one of 2T vectors), then the per-probe posteriors.
      const int K = 2 * T;
      auto kind = [&](size_t r) { return 2 * (int)pr[i0 + r].tier + (int)pr[i0 + r].hit; };
      float fbz, fcz;            // F[r] = (fbz + fcz F[r-1]) e[r]
      const float* fm1 = nullptr;  // forward message before the region (F[-1])
      float pu;
      uint32_t ppos;
      {  // forward, probe 0
        const int k0 = kind(0);
        const float* e = &PE[(size_t)k0 * SP];
        const float eu = PEu[k0];
        float* fc = &W.rfa[0];
        float fu, iz;
        if (b0 == 0) {  // contig start: initial distribution
          fu = 0.5f * eu;
          const float z = fu + 0.5f * pi * PEs[k0];
          iz = 1.0f / z;
          fbz = 0.5f * pi * iz;
          fcz = 0.0f;
          for (size_t j = 0; j < SP; ++j) fc[j] = fbz * e[j];
        } else {
          fm1 = &W.Ff[(b0 - 1) * SP];
          const float pu0 = W.Fu[b0 - 1];
          const uint32_t d = pr[i0].pos - W.bstart[b0 - 1];
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
          const float pa_sum = std::max(0.0f, 1.0f - pu0);
          const float base = pu0 * u * pi + stay * h * pi * pa_sum, c1 = stay * (1.0f - h);
          fu = (pu0 * stay + pa_sum * u) * eu;
          float dot = 0;
          for (size_t j = 0; j < SP; ++j) dot += fm1[j] * e[j];
          const float z = fu + base * PEs[k0] + c1 * dot;
          iz = 1.0f / z;
          fbz = base * iz;
          fcz = c1 * iz;
          for (size_t j = 0; j < SP; ++j) fc[j] = (fbz + fcz * fm1[j]) * e[j];
        }
        W.rfu[0] = fu * iz;
        pu = W.rfu[0];
        ppos = pr[i0].pos;
      }
      float bbz, bcz, bu;          // B[r] = bbz + bcz (e[r+1] B[r+1])
      const float* bnext_e = nullptr;  // e and B of the step after the region's last probe
      const float* bnext_b = nullptr;
      {  // backward message entering the region's last probe
        float* bo = &W.rbb[(nr - 1) * SP];
        if (b1 + 1 < nb) {
          bnext_e = W.ep[b1 + 1];
          bnext_b = &W.Bb[(b1 + 1) * SP];
          const float wu = W.Eu[b1 + 1] * W.Bbu[b1 + 1];
          float sw = 0;
          for (size_t j = 0; j < SP; ++j) sw += bnext_e[j] * bnext_b[j];
          const uint32_t d = W.bstart[b1 + 1] - pr[i1 - 1].pos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * pi;
          const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
          const float zz = b0v + (float)m * base + c1 * sw, iz = 1.0f / zz;
          bbz = base * iz;
          bcz = c1 * iz;
          for (size_t j = 0; j < SP; ++j) bo[j] = bbz + bcz * (bnext_e[j] * bnext_b[j]);
          bu = b0v * iz;
        } else {
          for (size_t j = 0; j < SP; ++j) bo[j] = 1.0f;
          bbz = 1.0f;
          bcz = 0.0f;
          bu = 1.0f;
        }
        W.rbu[nr - 1] = bu;
      }
      for (size_t i = 1; i < nr; ++i) {
        {  // forward step into probe i
          const size_t r = i;
          const int kr = kind(r), kp = kind(r - 1);
          const float* __restrict e = &PE[(size_t)kr * SP];
          const float* __restrict ep1 = &PE[(size_t)kp * SP];
          const float* __restrict fprev = &W.rfa[(r - 1) * SP];
          const float* __restrict fpp = r >= 2 ? &W.rfa[(r - 2) * SP] : fm1;
          float s2 = 0;
          if (fpp)
            for (size_t j = 0; j < SP; ++j) s2 += fpp[j] * (ep1[j] * e[j]);
          const float dot = fbz * PED[kp * K + kr] + fcz * s2;
          const Probe& p = pr[i0 + r];
          const uint32_t d = p.pos - ppos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
          const float pa_sum = std::max(0.0f, 1.0f - pu);
          const float base = pu * u * pi + stay * h * pi * pa_sum, c1 = stay * (1.0f - h);
          const float fu = (pu * stay + pa_sum * u) * PEu[kr];
          const float z = fu + base * PEs[kr] + c1 * dot;
          const float iz = 1.0f / z, bz = base * iz, cz = c1 * iz;
          float* __restrict fc = &W.rfa[r * SP];
          for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * fprev[j]) * e[j];
          fbz = bz;
          fcz = cz;
          W.rfu[r] = fu * iz;
          pu = W.rfu[r];
          ppos = p.pos;
        }
        {  // backward step out of probe r = nr - i (gives B[r-1])
          const size_t r = nr - i;
          const int kr = kind(r);
          const float* __restrict e = &PE[(size_t)kr * SP];
          const float* __restrict bb = &W.rbb[r * SP];
          const float* __restrict en = r + 1 < nr ? &PE[(size_t)kind(r + 1) * SP] : bnext_e;
          const float* __restrict bn = r + 1 < nr ? &W.rbb[(r + 1) * SP] : bnext_b;
          float t2 = 0;
          if (en)
            for (size_t j = 0; j < SP; ++j) t2 += e[j] * (en[j] * bn[j]);
          const float sw = bbz * PEs[kr] + bcz * t2;
          const float wu = PEu[kr] * bu;
          const uint32_t d = pr[i0 + r].pos - pr[i0 + r - 1].pos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * pi;
          const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
          const float zz = b0v + (float)m * base + c1 * sw, iz = 1.0f / zz;
          const float bz = base * iz, cz = c1 * iz;
          float* __restrict bo = &W.rbb[(r - 1) * SP];
          for (size_t j = 0; j < SP; ++j) bo[j] = bz + cz * (e[j] * bb[j]);
          bbz = bz;
          bcz = cz;
          bu = b0v * iz;
          W.rbu[r - 1] = bu;
        }
      }
      for (size_t r = 0; r < nr; ++r) {  // per-probe aligned posterior
        const float* __restrict fa = &W.rfa[r * SP];
        const float* __restrict bb = &W.rbb[r * SP];
        float gu = W.rfu[r] * W.rbu[r], tot = gu;
        for (size_t j = 0; j < SP; ++j) tot += fa[j] * bb[j];
        W.ap[i0 + r] = 1.0 - (double)gu / (double)tot;
      }
      // per-bin mean posterior for bin-level uses (heterogeneity windows)
      for (size_t q = b0; q <= b1; ++q) {
        double acc = 0;
        for (uint32_t i = W.bfirst[q]; i < W.bfirst[q + 1]; ++i) acc += W.ap[i];
        W.ab[q] = acc / (double)(W.bfirst[q + 1] - W.bfirst[q]);
      }
    }
  }

  // emission vector of a bin with per-tier miss/hit counts c (B = power-table stride)
  void emission(const uint8_t* c, size_t B, float* __restrict e, float& eu_out, float& es_out) const {
    const size_t SP = SPA;
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
    eu_out = (float)eu;
    float se = 0;
    for (size_t j = 0; j < SP; ++j) se += e[j];
    es_out = se;
  }

  void build_emission_cache(size_t Bh, size_t B) {
    EB = Bh;
    size_t nkey = 1;
    for (int q = 0; q < 2 * T; ++q) nkey *= EB + 1;
    ecix.assign(nkey, -1);
    uint8_t c[2 * MAXT] = {0};
    for (size_t key = 0; key < nkey; ++key) {
      size_t r = key, sum = 0;
      for (int q = 2 * T - 1; q >= 0; --q) {
        c[q] = (uint8_t)(r % (EB + 1));
        r /= EB + 1;
        sum += c[q];
      }
      if (sum == 0 || sum > EB) continue;
      const size_t slot = ecU.size();
      ecix[key] = (int32_t)slot;
      ecE.resize((slot + 1) * SPA);
      ecU.push_back(0);
      ecS.push_back(0);
      emission(c, B, &ecE[slot * SPA], ecU[slot], ecS[slot]);
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
#ifdef AIF_AVX512
    // counts of full bins from tier/hit bit masks: 8 probes (8-byte structs) per load
    const bool vec = T <= 2 && sizeof(Probe) == 8 && B % 8 == 0 && B <= 64;
    static_assert(offsetof(Probe, tier) == 4 && offsetof(Probe, hit) == 5, "Probe layout");
    const __m512i tbit = _mm512_set1_epi64(1LL << 32), hbit = _mm512_set1_epi64(1LL << 40);
#endif
    for (size_t bi = 0; bi < nb; ++bi) {
      size_t f = bi * B, l = std::min(n, f + B);
      W.bstart[bi] = pr[f].pos;
      W.blast[bi] = pr[l - 1].pos;
      W.bfirst[bi] = (uint32_t)f;
      uint8_t* c = &W.bcnt[bi * 2 * T];
#ifdef AIF_AVX512
      if (vec && l - f == B) {
        uint64_t tm = 0, hm = 0;
        for (size_t g = 0; g < B; g += 8) {
          const __m512i v = _mm512_loadu_si512((const void*)&pr[f + g]);
          tm |= (uint64_t)_mm512_test_epi64_mask(v, tbit) << g;
          hm |= (uint64_t)_mm512_test_epi64_mask(v, hbit) << g;
        }
        const int h1 = __builtin_popcountll(hm & tm), h0 = __builtin_popcountll(hm & ~tm);
        const int n1 = __builtin_popcountll(tm);
        c[1] = (uint8_t)h0;
        c[0] = (uint8_t)(B - n1 - h0);
        if (T == 2) {
          c[3] = (uint8_t)h1;
          c[2] = (uint8_t)(n1 - h1);
        }
      } else
#endif
        for (size_t i = f; i < l; ++i) ++c[2 * pr[i].tier + pr[i].hit];
      for (int t = 0; t < T; ++t) {
        llU += c[2 * t] * logU[2 * t] + c[2 * t + 1] * logU[2 * t + 1];
        hits += c[2 * t + 1];
      }
    }
    W.bfirst[nb] = (uint32_t)n;
    return llU;
  }

  void hmm(Work& W, double& ll, double& pm_ani, double& low_frac, bool fine = false) const {
    const size_t m = gth.size(), SP = SPA;
    const size_t B = (size_t)std::max(1, std::max(P.bin, P.bin_het));  // power-table stride
    const size_t nb = W.bstart.size();
    const float pi = 1.0f / (float)m;
    // emissions: fine bins (no refinement follows) take vectors from the cache, others are
    // computed into Ef (which refine() reads)
    W.ep.resize(nb);
    W.Eu.resize(nb);
    W.Es.resize(nb);
    const bool cached = fine && EB > 0;
    if (!cached) W.Ef.resize(nb * SP);
    for (size_t bi = 0; bi < nb; ++bi) {
      const uint8_t* c = &W.bcnt[bi * 2 * T];
      if (cached) {
        const int32_t slot = ecix[ec_key(c)];
        W.ep[bi] = &ecE[(size_t)slot * SP];
        W.Eu[bi] = ecU[slot];
        W.Es[bi] = ecS[slot];
      } else {
        emission(c, B, &W.Ef[bi * SP], W.Eu[bi], W.Es[bi]);
        W.ep[bi] = &W.Ef[bi * SP];
      }
    }
    const float* const* EP = W.ep.data();
    // forward / backward. Both recursions are latency-bound, so (1) they run interleaved in
    // one loop and (2) each step's normalising dot product is obtained from scalars of the
    // previous step and sums over older messages (look-ahead), keeping the horizontal
    // reductions off the dependency chain:
    //   F[b-1] = (bz + cz F[b-2]) e[b-1]  =>  <F[b-1], e[b]> = bz <e[b-1],e[b]> + cz <F[b-2] e[b-1], e[b]>
    //   B[b]   = bz' + cz' e[b+1] B[b+1]  =>  <e[b], B[b]>   = bz' Es[b] + cz' <e[b] e[b+1], B[b+1]>
    W.Ff.resize(nb * SP);
    W.Fu.resize(nb);
    W.ab.resize(nb);
    W.Bb.resize(nb * SP);
    W.Bbu.resize(nb);
    std::fill(W.Bb.begin() + (nb - 1) * SP, W.Bb.begin() + nb * SP, 1.0f);
    W.Bbu[nb - 1] = 1.0f;
    ll = 0;
    double zprod = 1.0;
    int zcnt = 0;
    float fbz, fcz;  // forward: F[b-1] = (fbz + fcz F[b-2]) e[b-1]
    {
      const float* e = EP[0];
      float* f = &W.Ff[0];
      const float fu = 0.5f * W.Eu[0], z = fu + 0.5f * pi * W.Es[0];
      const float iz = 1.0f / z;
      fbz = 0.5f * pi * iz;
      fcz = 0.0f;
      for (size_t j = 0; j < SP; ++j) f[j] = fbz * e[j];
      W.Fu[0] = fu * iz;
      ll += std::log((double)z);
    }
    float bu = 1.0f, bbz = 1.0f, bcz = 0.0f;  // backward: B[b] = bbz + bcz e[b+1] B[b+1]
    for (size_t i = 1; i < nb; ++i) {
      {  // forward step into bin i
        const size_t bi = i;
        uint32_t d = W.bstart[bi] - W.bstart[bi - 1];
        const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
        const float* __restrict fp = &W.Ff[(bi - 1) * SP];
        float* __restrict fc = &W.Ff[bi * SP];
        const float* __restrict e = EP[bi];
        const float* __restrict e1 = EP[bi - 1];
        float s1 = 0, s2 = 0;
        if (bi >= 2) {
          const float* __restrict fpp = &W.Ff[(bi - 2) * SP];
          for (size_t j = 0; j < SP; ++j) {
            const float ee = e1[j] * e[j];
            s1 += ee;
            s2 += fpp[j] * ee;
          }
        } else {
          for (size_t j = 0; j < SP; ++j) s1 += e1[j] * e[j];
        }
        const float dot = fbz * s1 + fcz * s2;
        const float fpu = W.Fu[bi - 1];
        const float SA = std::max(0.0f, 1.0f - fpu);  // previous step is normalised
        const float base = fpu * u * pi + stay * h * pi * SA;
        const float c1 = stay * (1.0f - h);
        const float fu = (fpu * stay + SA * u) * W.Eu[bi];
        const float z = fu + base * W.Es[bi] + c1 * dot;
        const float iz = 1.0f / z, bz = base * iz, cz = c1 * iz;
        for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * fp[j]) * e[j];
        fbz = bz;
        fcz = cz;
        W.Fu[bi] = fu * iz;
        zprod *= (double)z;
        if (++zcnt == 16) {
          ll += std::log(zprod);
          zprod = 1.0;
          zcnt = 0;
        }
      }
      {  // backward step out of bin nb - i
        const size_t bi = nb - i;
        uint32_t d = W.bstart[bi] - W.bstart[bi - 1];
        const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
        const float* __restrict e = EP[bi];
        const float* __restrict bb = &W.Bb[bi * SP];
        float* __restrict bo = &W.Bb[(bi - 1) * SP];
        float t2 = 0;
        if (bi + 1 < nb) {
          const float* __restrict e1 = EP[bi + 1];
          const float* __restrict bn = &W.Bb[(bi + 1) * SP];
          for (size_t j = 0; j < SP; ++j) t2 += e[j] * e1[j] * bn[j];
        }
        const float sw = bbz * W.Es[bi] + bcz * t2;
        const float wu = W.Eu[bi] * bu;
        const float SW = sw * pi;
        const float base = u * wu + stay * h * SW;
        const float c1 = stay * (1.0f - h);
        const float b0 = stay * wu + u * SW;
        const float zz = b0 + (float)m * base + c1 * sw;  // any positive scale works
        const float iz = 1.0f / zz, bz = base * iz, cz = c1 * iz;
        for (size_t j = 0; j < SP; ++j) bo[j] = bz + cz * (e[j] * bb[j]);
        bbz = bz;
        bcz = cz;
        bu = b0 * iz;
        W.Bbu[bi - 1] = bu;
      }
    }
    ll += std::log(zprod);
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
    for (size_t bi = nb; bi-- > 0;) post(bi, &W.Ff[bi * SP], W.Fu[bi], &W.Bb[bi * SP], W.Bbu[bi]);
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

  Result estimate(std::string_view seq, Work& W) const {
    Result R;
    const int k = sk.k;
    PSTART;
    probe(seq, W);
    PSTAMP(0);
    const std::vector<Probe>& pr = W.pr;
    R.probes = pr.size();
    if (pr.size() < 3) {
      for (auto& p : pr) R.hits += p.hit;
      R.status = "too_few_probes";
      return R;
    }
    const double llU = make_bins(W, R.hits, P.bin);
    PSTAMP(1);
    const size_t nb = W.bstart.size();
    double ll = 0, pm_ani = NAN, low_frac = 0, fast_lr = 0;
    const bool fast = all_aligned(W, fast_lr);
    if (fast) {  // confidently aligned everywhere: the posterior is ~1 in every bin
      W.ab.assign(nb, 1.0);
      R.lr = fast_lr;
    } else {
      hmm(W, ll, pm_ani, low_frac);
      R.lr = ll - llU;
      PSTAMP(2);
      if (P.refine && P.bin > 1) refine(W);
      else W.refined.assign(nb, 0);
      PSTAMP(3);
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
    PSTAMP(4);
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
    PSTAMP(5);
    if (het) {  // posterior-mean ANI from a finer-binned HMM (long contigs: the main bins suffice)
      if (P.bin_het != P.bin && pr.size() <= P.het_fine_max) {
        size_t hh = 0;
        make_bins(W, hh, P.bin_het);
        hmm(W, ll, pm_ani, low_frac, true);
        R.low_frac = low_frac;
      } else if (fast) {  // the fast path skipped the HMM; the main bins are still in W
        hmm(W, ll, pm_ani, low_frac);
        R.low_frac = low_frac;
      }
    }
    if (het && pm_ani == pm_ani) {
      R.het = true;
      R.ani = pm_ani;
    }
    if (R.low_frac > P.floor_frac) R.status = "below_floor";
    PSTAMP(6);
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
          "aif query <ref.aif> <contigs.fa> [-t threads] [--bin 32] [--bin-het 8] [--het-fine-max N (off)] [--fast-win 32] [--fast-llr 10] [--rho-u 3e-6] [--rho-h 1e-4] [--ident-c 8]\n"
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
    else if (a == "--het-fine-max") P.het_fine_max = (size_t)std::max(0.0, atof(argv[++i]));
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
    prof_report();
    return 0;
  }
  // multi-threaded streaming pipeline. The main thread reads contigs into chunks (one long
  // contig, or a group of short ones) and queues them; workers estimate whole chunks and
  // format their output; a writer thread prints chunks in input order. Outstanding input is
  // capped (MAX_OUT bases), so memory stays bounded and no thread waits on a batch barrier.
  // Memory-mapped input: chunks carry raw spans of the map and workers strip line breaks
  // into their own reused buffers, so the reader only locates records. Other input (pipes):
  // the reader copies sequences into the chunk.
  struct Chunk {
    size_t id = 0, bytes = 0;
    std::vector<std::string> names;
    std::string data;
    std::vector<size_t> off{0};
    std::vector<std::pair<const char*, const char*>> spans;
    std::string out;
  };
  const bool mapped = fa.mapped();
  std::vector<std::string> sbuf(threads);
  const size_t CHUNK_BASES = (size_t)1 << 18, CHUNK_CONTIGS = 64, MAX_OUT = (size_t)128 << 20;
  std::mutex mu;
  std::condition_variable cv_work, cv_done, cv_space;
  std::deque<std::unique_ptr<Chunk>> todo;
  std::map<size_t, std::unique_ptr<Chunk>> done;
  std::vector<std::unique_ptr<Chunk>> spare;  // written chunks, recycled (no fresh-page faults)
  size_t outstanding = 0, produced = 0;
  bool reading = true;
  std::vector<Work> works(threads);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&, t]() {
#if defined(__x86_64__)
      _mm_setcsr(_mm_getcsr() | 0x8040);
#endif
      for (;;) {
        std::unique_ptr<Chunk> c;
        {
          std::unique_lock<std::mutex> lk(mu);
          cv_work.wait(lk, [&] { return !todo.empty() || !reading; });
          if (todo.empty()) return;
          c = std::move(todo.front());
          todo.pop_front();
        }
        if (mapped) {
          for (size_t i = 0; i < c->spans.size(); ++i) {
            sbuf[t].clear();
            append_stripped(sbuf[t], c->spans[i].first, c->spans[i].second);
            format_result(c->out, c->names[i], est.estimate(sbuf[t], works[t]));
          }
        } else {
          for (size_t i = 0; i + 1 < c->off.size(); ++i)
            format_result(c->out, c->names[i],
                          est.estimate(std::string_view(c->data.data() + c->off[i], c->off[i + 1] - c->off[i]), works[t]));
        }
        {
          std::lock_guard<std::mutex> lk(mu);
          const size_t id = c->id;
          done.emplace(id, std::move(c));
        }
        cv_done.notify_one();
      }
    });
  std::thread writer([&]() {
    size_t next = 0;
    for (;;) {
      std::unique_ptr<Chunk> c;
      {
        std::unique_lock<std::mutex> lk(mu);
        cv_done.wait(lk, [&] { return done.count(next) || (!reading && next == produced); });
        if (!done.count(next)) return;
        c = std::move(done[next]);
        done.erase(next);
      }
      fwrite(c->out.data(), 1, c->out.size(), stdout);
      // chunks are written in input order, so every mapped page before this chunk's end is done
      if (mapped && !c->spans.empty()) fa.release_upto(c->spans.back().second);
      {
        std::lock_guard<std::mutex> lk(mu);
        outstanding -= c->bytes;
        spare.push_back(std::move(c));
      }
      cv_space.notify_one();
      ++next;
    }
  });
  std::string name;
  for (;;) {
    std::unique_ptr<Chunk> c;
    {
      std::lock_guard<std::mutex> lk(mu);
      if (!spare.empty()) {
        c = std::move(spare.back());
        spare.pop_back();
      }
    }
    if (c) {
      c->names.clear();
      c->data.clear();
      c->off.assign(1, 0);
      c->spans.clear();
      c->out.clear();
      if (c->data.capacity() > 2 * CHUNK_BASES) std::string().swap(c->data);  // drop a long contig's buffer
    } else {
      c = std::make_unique<Chunk>();
    }
    c->id = produced;
    c->bytes = 0;
    while (mapped && c->bytes < CHUNK_BASES && c->names.size() < CHUNK_CONTIGS) {
      const char *b, *e;
      if (!fa.next_span(name, b, e)) break;
      c->names.push_back(name);
      c->spans.emplace_back(b, e);
      c->bytes += (size_t)(e - b);
    }
    while (!mapped && c->data.size() < CHUNK_BASES && c->names.size() < CHUNK_CONTIGS) {
      const size_t before = c->data.size();
#ifdef AIF_PROF
      const uint64_t _r0 = __rdtsc();
      const bool got = fa.next(name, c->data, true);
      g_prof[9] += __rdtsc() - _r0;
      if (!got) break;
#else
      if (!fa.next(name, c->data, true)) break;
#endif
      c->names.push_back(name);
      c->off.push_back(c->data.size());
      if (before == 0 && c->data.size() >= CHUNK_BASES) break;  // a long contig travels alone
    }
    if (c->names.empty()) break;
    if (!mapped) c->bytes = c->data.size();
    const size_t nb = c->bytes;
    {
      std::unique_lock<std::mutex> lk(mu);
      cv_space.wait(lk, [&] { return outstanding == 0 || outstanding + nb <= MAX_OUT; });
      outstanding += nb;
      todo.push_back(std::move(c));
      ++produced;
    }
    cv_work.notify_one();
  }
  {
    std::lock_guard<std::mutex> lk(mu);
    reading = false;
  }
  cv_work.notify_all();
  cv_done.notify_all();
  for (auto& th : pool) th.join();
  cv_done.notify_all();
  writer.join();
  prof_report();
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
