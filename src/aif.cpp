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
// v1.4 (uncertainty): per-contig confidence intervals for ANI and AF. Statistical SEs come
// from Theorem 1 (FP-corrected containment, conditioned on the realised substitutions) and from
// the posterior uncertainty of alignment boundaries; calibrated intervals (default) use
// stratified split-conformal multipliers fitted on the real tuning set (src/ci_table.inc,
// bench_real/calibration.py). Optional AF-model variants: a local identity prior
// (--ident-window) and overdispersed bin emissions (--overdisp); both off by default.
//
// v1.3 (single-thread speed): forward/backward recursions interleaved, with look-ahead
// normalisers that keep horizontal reductions off the dependency chain (HMM and
// refinement; rounding-level output changes, <= 1e-6); cached emissions for fine bins;
// GFNI reverse complements; AVX-512 2-bit packing and mask-popcount bin counts; batched,
// branch-free calibration. Sketches are byte-identical to v1.1/v1.2.
//
// v1.5 (database mode): "aif db build/search/info". References share one k and tier layout,
// so a contig is prepared once (selection, filter hashes) and each candidate reference costs
// only its filter tests and the HMM; candidates come from a bit-sliced screening index.
// Exhaustive search is byte-identical to pairwise queries against "-k K" sketches.
//
// v1.6 (GTDB scale, identity maps): gzip input; --segments writes aligned/unaligned segments and
// windowed identity with calibrated intervals; db search adds a pre-HMM local-score filter, caps
// screening keys per contig, and orders hits by the model's likelihood ratio.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread -o aif aif.cpp -lz
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
#include <unordered_map>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
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
  gzFile gz = nullptr;  // gzip-compressed input
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
    size_t got;
    if (gz) {
      const int r = gzread(gz, buf.data() + end, (unsigned)std::min<size_t>(buf.size() - end, 1u << 30));
      if (r < 0) throw std::runtime_error("gzip read error");
      got = (size_t)r;
    } else {
      got = fread(buf.data() + end, 1, buf.size() - end, f);
    }
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
    unsigned char magic[2] = {0, 0};
    if (fd >= 0 && fstat(fd, &st) == 0 && S_ISREG(st.st_mode) && st.st_size >= 2 && pread(fd, magic, 2, 0) == 2 &&
        magic[0] == 0x1f && magic[1] == 0x8b) {  // gzip: buffered reading through zlib
      close(fd);
      fd = -1;
      buf.resize(1 << 22);
      gz = gzopen(path, "rb");
      if (!gz) throw std::runtime_error(std::string("cannot open ") + path);
      gzbuffer(gz, 1 << 20);
      return;
    }
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
    if (gz) gzclose(gz);
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
  const uint64_t* bp = nullptr;  // filter words in a memory-mapped database (else bits)
  inline const uint64_t* data() const { return bp ? bp : bits.data(); }

  double s() const { return (double)(hi - lo) / 18446744073709551616.0; }
  inline uint64_t block_index(uint64_t h) const { return ((h >> 32) * nblocks) >> 32; }
  inline uint64_t* block_of(uint64_t x, uint64_t& h) {
    h = fmix64(x ^ SEED_BF);
    return &bits[block_index(h) * 8];
  }
  inline const uint64_t* block_of(uint64_t x, uint64_t& h) const {
    h = fmix64(x ^ SEED_BF);
    return data() + block_index(h) * 8;
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
  inline bool contains_at(const uint64_t* blk, uint64_t h) const { return contains_h2(blk, fmix64(h ^ SEED_BIT)); }
  inline bool contains_h2(const uint64_t* blk, uint64_t h2) const {
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
    for (auto& t : tiers) b += t.nblocks * 64;
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
  double ident_window = 0.0;   // local identity prior: aligned levels within this ANI of the local identity (0 = off)
  int ident_win_probes = 256;  // window (probes) of the local identity profile
  int ident_radius = 3;        // local identity = median over windows within +-radius (>= 10000: whole contig)
  double overdisp = 0.0;       // scale of the k-mer-overlap overdispersion of aligned bin counts (0 = binomial)
  double ci_level = 0.95;      // two-sided level of the reported intervals (0.90, 0.95 or 0.99 calibrated)
  bool ci_stat = false;        // intervals from the statistical SE alone (normal quantile), no calibration
  std::string ci_table;        // file with a calibration table (default: built in, from the real tuning set)
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
  double ani_se = NAN, af_se = NAN;        // standard errors (statistical + model discrepancy)
  double ani_lo = NAN, ani_hi = NAN, af_lo = NAN, af_hi = NAN;  // confidence intervals (P.ci_level)
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
  std::vector<uint16_t> wj;             // local identity prior: allowed levels 0..wj-1 per probe window
  size_t wprobes = 0;                   // probes per window (0: prior off)
  std::vector<std::pair<double, double>> wtheta;  // scratch: (window containment, probes)
  std::vector<double> wth, wn;          // scratch: per-window containment and probe counts
  std::vector<uint16_t> bjc;            // allowed level count per bin of the current binning
  std::vector<std::pair<double, double>> hw;  // heterogeneity windows: (containment, aligned length)
  double se_pooled = 0;                 // statistical SE of the pooled ANI (for the heterogeneous CI)
  std::vector<double> cxs, cas;         // scratch for the AF interval
  double afvar = 0;                     // AF boundary variance (bp^2), accumulated by refine()
  std::vector<uint8_t> refined;         // per bin: per-probe posteriors available in ap
  std::vector<double> ap;               // per-probe aligned posterior (refined bins only)
  std::vector<float> rfa, rfu, rb, rw;  // refinement scratch
  std::vector<float> rbb, rbu;          // refinement: per-probe backward messages
  std::vector<uint8_t> bcnt;  // per bin: hits and misses per tier
  std::vector<int32_t> bslot; // per bin: emission-cache slot (fine bins only)
  std::vector<double> E, F, a, b, w, conf;
  // database search (db_prepare / db_hits): reference-independent per-probe hashes
  std::vector<uint32_t> pbh;        // Bloom block hash (high 32 bits) per probe
  std::vector<uint64_t> ph2;        // bit-position hash per probe
  std::vector<uint64_t> skey;       // screening keys of the contig (block hashes)
  std::vector<uint64_t> skh;        // their selection hashes (for capping the number of keys)
  std::vector<float> psw;           // per-probe switching cost (local-score filter)
  bool want_seg = false;            // keep per-probe aligned posteriors for segments()
  std::vector<float> pa;            // per-probe aligned posterior (main pass)
  std::vector<uint16_t> scnt16;     // screening counts (per group, short-term)
  std::vector<uint32_t> scnt;       // screening counts per reference
  std::vector<std::pair<float, uint32_t>> cand;  // (screening containment, reference)
  std::vector<std::pair<float, uint32_t>> cand2; // (probe containment, reference)
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

// level mask value inside templates on M (local identity prior): constant 1 when the prior is off
#define MKV(p, j) (M ? (p)[j] : 1.0f)

// Built-in calibration table for the confidence intervals: conformal multipliers [level 90/95/99%]
// [stratum][ANI, AF], fitted on the real tuning set by bench_real/calibration.py (see there).
static const double CI_DEFAULT[3][24][2] = {
#include "ci_table.inc"
};

struct Estimator {
  const Sketch& sk;
  Params P;
  int T;
  std::vector<double> epsp, stier, gth, gani;
  size_t SPA = 16;                 // aligned states padded to a multiple of 16
  std::vector<float> PHa, PMa;     // [t][c][aligned state] = (hit prob)^c, (miss prob)^c
  // Overdispersed aligned emissions: overlapping k-mers share bases, so probe survivals within a
  // bin are positively correlated under the simple mutation model. A bin's tier-t hit count is
  // modelled as beta-binomial with the intra-bin correlation implied by that overlap structure
  // (density s_t, k, identity), scaled by P.overdisp. BBa[t][hits][misses][state] holds the
  // probability of one ordered outcome (the binomial coefficient cancels against U's).
  std::vector<float> BBa;
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
  const double* tsu = nullptr;     // transition tables (own, or shared with a prototype Estimator)
  const double* thh = nullptr;
  // emission cache for fine (bin_het) bins: one vector per hit/miss count combination
  // (computed exactly as in hmm()); ecix maps the mixed-radix counts to a slot
  size_t EB = 0;                   // counts per tier and outcome are < EB + 1
  std::vector<int32_t> ecix;
  std::vector<float> ecE, ecU, ecS, ecC;  // ecC: prefix sums of ecE over levels (masked row sums)
  std::vector<float> MT;                  // MT[c*SP + j] = (j < c): prefix masks, c = 0..m
  std::vector<float> PIc;                 // 1 / c
  std::vector<float> PEcum, PEDcum;       // prefix sums over levels of PE rows and of PE row products
  inline size_t ec_key(const uint8_t* c) const {
    size_t key = 0;
    for (int q = 0; q < 2 * T; ++q) key = key * (EB + 1) + c[q];
    return key;
  }

  Estimator(const Estimator&) = delete;
  Estimator& operator=(const Estimator&) = delete;
  // proto: an Estimator with the same Params whose distance tables are shared (database mode)
  Estimator(const Sketch& s, const Params& p, const Estimator* proto = nullptr) : sk(s), P(p), T((int)s.tiers.size()) {
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
    if (P.overdisp > 0) {
      BBa.assign((size_t)T * (B + 1) * (B + 1) * SPA, 0.0f);
      for (int t = 0; t < T; ++t) {
        const double e = epsp[t], st = stier[t];
        for (size_t j = 0; j < m; ++j) {
          const double th = gth[j], q = gani[j], pi = th + (1 - th) * e;
          double sum = 0;
          if (th < 1)
            for (int d = 1; d < sk.k; ++d) sum += (std::pow(q, d) - th) / (1 - th);
          const double Dm1 = std::max(0.0, 2.0 * st * sum * P.overdisp);  // design effect - 1 (survivals)
          for (size_t c = 0; c <= B; ++c)
            for (size_t mm = 0; c + mm <= B; ++mm) {
              const size_t n = c + mm;
              double lp;
              const double rho = n >= 2 ? (1 - e) * (1 - e) * th * (1 - th) * Dm1 / ((double)(n - 1) * pi * (1 - pi)) : 0;
              if (rho < 1e-7) {
                lp = (double)c * std::log(pi) + (double)mm * std::log(1 - pi);
              } else {
                const double A = 1.0 / std::min(rho, 0.5) - 1.0, a = pi * A, b = (1 - pi) * A;
                lp = std::lgamma(a + c) - std::lgamma(a) + std::lgamma(b + mm) - std::lgamma(b) - std::lgamma(A + n) + std::lgamma(A);
              }
              BBa[(((size_t)t * (B + 1) + c) * (B + 1) + mm) * SPA + j] = (float)std::exp(lp);
            }
        }
      }
    }
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
    PIc.assign(m + 1, 0.0f);
    for (size_t c = 1; c <= m; ++c) PIc[c] = 1.0f / (float)c;
    MT.assign((m + 1) * SPA, 0.0f);
    for (size_t c = 0; c <= m; ++c)
      for (size_t j = 0; j < c; ++j) MT[c * SPA + j] = 1.0f;
    PEcum.assign((size_t)2 * T * SPA, 0.0f);
    PEDcum.assign((size_t)4 * T * T * SPA, 0.0f);
    for (int a = 0; a < 2 * T; ++a) {
      float acc = 0;
      for (size_t j = 0; j < SPA; ++j) PEcum[(size_t)a * SPA + j] = (acc += PE[(size_t)a * SPA + j]);
      for (int b = 0; b < 2 * T; ++b) {
        float d = 0;
        for (size_t j = 0; j < SPA; ++j)
          PEDcum[((size_t)a * 2 * T + b) * SPA + j] = (d += PE[(size_t)a * SPA + j] * PE[(size_t)b * SPA + j]);
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
    if (P.bin_het > 0 && P.bin_het < P.bin && T <= 2) build_emission_cache((size_t)P.bin_het, B);
    load_ci_table();
    if (proto) {
      tsu = proto->tsu;
      thh = proto->thh;
    } else {
      tab_su.resize(QMAX);
      tab_h.resize(QMAX);
      for (size_t d = 0; d < QMAX; ++d) {
        tab_su[d] = std::exp(-P.rho_u * (double)d);
        tab_h[d] = -std::expm1(-P.rho_h * (double)d);
      }
      tsu = tab_su.data();
      thh = tab_h.data();
    }
  }
  inline double su(uint32_t d) const { return d < QMAX ? tsu[d] : std::exp(-P.rho_u * d); }
  inline double hh(uint32_t d) const { return d < QMAX ? thh[d] : -std::expm1(-P.rho_h * d); }

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
      const long long* base0 = (const long long*)t0.data();
      const long long delta1 = (long long)((const long long*)t1.data() - base0);
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
    if (W.wprobes)
      refine_t<true>(W);
    else
      refine_t<false>(W);
  }
  template <bool M>
  void refine_t(Work& W) const {
    W.afvar = 0;
    const size_t nb = W.bstart.size(), SP = SPA, m = gth.size();
    const float* MTp = MT.data();
    // masked per-probe emission sums / dots for a probe kind and an allowed level count
    auto pes = [&](int kk, size_t c) { return (!M || c >= m) ? PEs[kk] : PEcum[(size_t)kk * SP + c - 1]; };
    auto ped = [&](int ka, int kb, size_t c) { return PEDcum[((size_t)ka * 2 * T + kb) * SP + (M ? c : m) - 1]; };
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
      auto kind = [&](size_t r) { return 2 * (int)pr[i0 + r].tier + (int)pr[i0 + r].hit; };
      float fbz, fcz;            // F[r] = (fbz + fcz F[r-1]) e[r]
      const float* fm1 = nullptr;  // forward message before the region (F[-1])
      float pu;
      uint32_t ppos;
      {  // forward, probe 0
        const int k0 = kind(0);
        const float* e = &PE[(size_t)k0 * SP];
        const float eu = PEu[k0];
        const size_t jc0 = M ? jc_at(W, i0) : m;
        const float* mk0 = MTp + jc0 * SP;
        const float pi = PIc[jc0];
        float* fc = &W.rfa[0];
        float fu, iz;
        if (b0 == 0) {  // contig start: initial distribution
          fu = 0.5f * eu;
          const float z = fu + 0.5f * pi * pes(k0, jc0);
          iz = 1.0f / z;
          fbz = 0.5f * pi * iz;
          fcz = 0.0f;
          for (size_t j = 0; j < SP; ++j) fc[j] = fbz * e[j] * MKV(mk0, j);
        } else {
          fm1 = &W.Ff[(b0 - 1) * SP];
          const float pu0 = W.Fu[b0 - 1];
          const uint32_t d = pr[i0].pos - W.bstart[b0 - 1];
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
          const float pa_sum = std::max(0.0f, 1.0f - pu0);
          const float base = pu0 * u * pi + stay * h * pi * pa_sum, c1 = stay * (1.0f - h);
          fu = (pu0 * stay + pa_sum * u) * eu;
          float dot = 0;
          for (size_t j = 0; j < SP; ++j) dot += fm1[j] * e[j] * MKV(mk0, j);
          const float z = fu + base * pes(k0, jc0) + c1 * dot;
          iz = 1.0f / z;
          fbz = base * iz;
          fcz = c1 * iz;
          for (size_t j = 0; j < SP; ++j) fc[j] = (fbz + fcz * fm1[j]) * e[j] * MKV(mk0, j);
        }
        W.rfu[0] = fu * iz;
        pu = W.rfu[0];
        ppos = pr[i0].pos;
      }
      float bbz, bcz, bu;          // B[r] = bbz + bcz (e[r+1] B[r+1])
      const float* bnext_e = nullptr;  // e and B of the step after the region's last probe
      const float* bnext_b = nullptr;
      size_t jcn = m;                  // allowed levels of that step
      {  // backward message entering the region's last probe
        float* bo = &W.rbb[(nr - 1) * SP];
        if (b1 + 1 < nb) {
          bnext_e = W.ep[b1 + 1];
          bnext_b = &W.Bb[(b1 + 1) * SP];
          jcn = M ? W.bjc[b1 + 1] : m;
          const float* mkn = MTp + jcn * SP;
          const float wu = W.Eu[b1 + 1] * W.Bbu[b1 + 1];
          float sw = 0;
          for (size_t j = 0; j < SP; ++j) sw += bnext_e[j] * MKV(mkn, j) * bnext_b[j];
          const uint32_t d = W.bstart[b1 + 1] - pr[i1 - 1].pos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * PIc[jcn];
          const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
          const float zz = b0v + (float)jcn * base + c1 * sw, iz = 1.0f / zz;
          bbz = base * iz;
          bcz = c1 * iz;
          for (size_t j = 0; j < SP; ++j) bo[j] = bbz + bcz * (bnext_e[j] * MKV(mkn, j) * bnext_b[j]);
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
          const size_t jcr = M ? jc_at(W, i0 + r) : m, jcm = M ? std::min(jcr, jc_at(W, i0 + r - 1)) : m;
          const float* __restrict mmin = MTp + jcm * SP;
          const float* __restrict mkr = MTp + jcr * SP;
          const float pi = PIc[jcr];
          float s2 = 0;
          if (fpp)
            for (size_t j = 0; j < SP; ++j) s2 += fpp[j] * (ep1[j] * e[j]) * MKV(mmin, j);
          const float dot = fbz * ped(kp, kr, jcm) + fcz * s2;
          const Probe& p = pr[i0 + r];
          const uint32_t d = p.pos - ppos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d);
          const float pa_sum = std::max(0.0f, 1.0f - pu);
          const float base = pu * u * pi + stay * h * pi * pa_sum, c1 = stay * (1.0f - h);
          const float fu = (pu * stay + pa_sum * u) * PEu[kr];
          const float z = fu + base * pes(kr, jcr) + c1 * dot;
          const float iz = 1.0f / z, bz = base * iz, cz = c1 * iz;
          float* __restrict fc = &W.rfa[r * SP];
          for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * fprev[j]) * e[j] * MKV(mkr, j);
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
          const size_t jcr = M ? jc_at(W, i0 + r) : m, jcx = M ? (r + 1 < nr ? jc_at(W, i0 + r + 1) : jcn) : m;
          const float* __restrict mmin = MTp + std::min(jcr, jcx) * SP;
          const float* __restrict mkr = MTp + jcr * SP;
          const float pi = PIc[jcr];
          float t2 = 0;
          if (en)
            for (size_t j = 0; j < SP; ++j) t2 += e[j] * (en[j] * bn[j]) * MKV(mmin, j);
          const float sw = bbz * pes(kr, jcr) + bcz * t2;
          const float wu = PEu[kr] * bu;
          const uint32_t d = pr[i0 + r].pos - pr[i0 + r - 1].pos;
          const float stay = (float)su(d), u = 1.0f - stay, h = (float)hh(d), SW = sw * pi;
          const float base = u * wu + stay * h * SW, c1 = stay * (1.0f - h), b0v = stay * wu + u * SW;
          const float zz = b0v + (float)jcr * base + c1 * sw, iz = 1.0f / zz;
          const float bz = base * iz, cz = c1 * iz;
          float* __restrict bo = &W.rbb[(r - 1) * SP];
          for (size_t j = 0; j < SP; ++j) bo[j] = bz + cz * (e[j] * MKV(mkr, j) * bb[j]);
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
        W.ap[i0 + r] = tot > 0 ? 1.0 - (double)gu / (double)tot : 0.0;
      }
      // AF boundary variance of this region (confidence interval), from the per-probe posterior
      W.afvar += boundary_var(W, i0, i1, [&](size_t i) { return W.ap[i]; }, b0 ? W.ab[b0 - 1] : 0.0,
                              b1 + 1 < nb ? W.ab[b1 + 1] : 0.0);
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
    if (!BBa.empty()) {
      const float* r0 = &BBa[(((size_t)0 * (B + 1) + c[1]) * (B + 1) + c[0]) * SP];
      for (size_t j = 0; j < SP; ++j) e[j] = r0[j];
      double eu = PHu[c[1]] * PMu[c[0]];
      for (int t = 1; t < T; ++t) {
        const float* rt = &BBa[(((size_t)t * (B + 1) + c[2 * t + 1]) * (B + 1) + c[2 * t]) * SP];
        for (size_t j = 0; j < SP; ++j) e[j] *= rt[j];
        eu *= PHu[t * (B + 1) + c[2 * t + 1]] * PMu[t * (B + 1) + c[2 * t]];
      }
      eu_out = (float)eu;
      float se = 0;
      for (size_t j = 0; j < SP; ++j) se += e[j];
      es_out = se;
      return;
    }
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
    ecC.resize(ecE.size());
    for (size_t sl = 0; sl < ecU.size(); ++sl) {
      float acc = 0;
      for (size_t j = 0; j < SPA; ++j) ecC[sl * SPA + j] = (acc += ecE[sl * SPA + j]);
    }
  }

  // Local identity prior. Uniform jumps between identity levels let an unaligned island be
  // absorbed as "aligned at ~85%" inside a high-identity region. Identity varies smoothly
  // along a genome except at real breakpoints, so the aligned levels allowed at a position are
  // restricted to within ident_window of the LOCAL identity: the probe-weighted median
  // FP-corrected containment of the aligned-looking windows (containment above the lowest
  // level's) within +-ident_radius windows of ident_win_probes probes. Islands (containment
  // ~0) and isolated conserved genes do not move the median; a long divergent tract sets its
  // own local identity and stays aligned. Fills W.wj (allowed level count per window).
  void set_local_prior(Work& W) const {
    const size_t m = gth.size(), nb = W.bstart.size();
    W.wprobes = 0;
    W.wj.clear();
    if (P.ident_window <= 0 || nb == 0) return;
    const size_t B = (size_t)std::max(1, P.bin);
    const size_t wbins = std::max<size_t>(1, ((size_t)P.ident_win_probes + B - 1) / B);
    const size_t nw = (nb + wbins - 1) / wbins;
    W.wprobes = wbins * B;
    W.wth.assign(nw, -1.0);
    W.wn.assign(nw, 0.0);
    for (size_t w = 0; w < nw; ++w) {
      double N[MAXT] = {0}, Y[MAXT] = {0}, tot = 0;
      for (size_t bi = w * wbins; bi < std::min(nb, (w + 1) * wbins); ++bi)
        for (int t = 0; t < T; ++t) {
          N[t] += W.bcnt[bi * 2 * T + 2 * t] + W.bcnt[bi * 2 * T + 2 * t + 1];
          Y[t] += W.bcnt[bi * 2 * T + 2 * t + 1];
        }
      for (int t = 0; t < T; ++t) tot += N[t];
      const double th = combine(N, Y, 0.5);
      if (th >= gth.back() && th == th) {
        W.wth[w] = std::min(th, 1.0);
        W.wn[w] = tot;
      }
    }
    W.wj.assign(nw, (uint16_t)m);
    const size_t R = (size_t)P.ident_radius;
    for (size_t w = 0; w < nw; ++w) {
      const size_t lo = (R >= nw || w < R) ? 0 : w - R, hi = (R >= nw) ? nw : std::min(nw, w + R + 1);
      auto& v = W.wtheta;
      v.clear();
      double wsum = 0;
      for (size_t q = lo; q < hi; ++q)
        if (W.wth[q] >= 0) {
          v.emplace_back(W.wth[q], W.wn[q]);
          wsum += W.wn[q];
        }
      if (v.empty()) continue;
      std::sort(v.begin(), v.end());
      double acc = 0, med = v.back().first;
      for (auto& x : v)
        if ((acc += x.second) >= 0.5 * wsum) {
          med = x.first;
          break;
        }
      const double floor = std::pow(med, 1.0 / sk.k) - P.ident_window;
      size_t jc = 0;
      while (jc < m && gani[jc] >= floor) ++jc;
      W.wj[w] = (uint16_t)std::max<size_t>(jc, std::min<size_t>(m, 2));
    }
  }
  // allowed level count at probe index i (m when the prior is off)
  inline size_t jc_at(const Work& W, size_t i) const {
    return W.wprobes ? W.wj[std::min(W.wj.size() - 1, i / W.wprobes)] : gth.size();
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
    if (EB > 0 && (size_t)binsz == EB) {  // fine bins: emission-cache slots (hmm reads them)
      W.bslot.resize(nb);
      for (size_t bi = 0; bi < nb; ++bi) W.bslot[bi] = ecix[ec_key(&W.bcnt[bi * 2 * T])];
    }
    return llU;
  }

  void hmm(Work& W, double& ll, double& pm_ani, double& low_frac, bool fine = false) const {
    if (W.wprobes)
      hmm_t<true>(W, ll, pm_ani, low_frac, fine);
    else
      hmm_t<false>(W, ll, pm_ani, low_frac, fine);
  }
  // M: local identity prior active (per-bin level masks); M = false compiles the masks away
  template <bool M>
  void hmm_t(Work& W, double& ll, double& pm_ani, double& low_frac, bool fine) const {
    const size_t m = gth.size(), SP = SPA;
    const size_t B = (size_t)std::max(1, std::max(P.bin, P.bin_het));  // power-table stride
    const size_t nb = W.bstart.size();
    // local identity prior: allowed aligned levels per bin (all m when off); jumps are uniform
    // over the levels allowed at the destination bin
    W.bjc.resize(nb);
    if (M)
      for (size_t bi = 0; bi < nb; ++bi) W.bjc[bi] = (uint16_t)jc_at(W, W.bfirst[bi]);
    else
      std::fill(W.bjc.begin(), W.bjc.end(), (uint16_t)m);
    const uint16_t* JC = W.bjc.data();
    const float* MTp = MT.data();
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
        const int32_t slot = W.bslot[bi];
        W.ep[bi] = &ecE[(size_t)slot * SP];
        W.Eu[bi] = ecU[slot];
        W.Es[bi] = (M && JC[bi] < m) ? ecC[(size_t)slot * SP + JC[bi] - 1] : ecS[slot];
      } else {
        emission(c, B, &W.Ef[bi * SP], W.Eu[bi], W.Es[bi]);
        W.ep[bi] = &W.Ef[bi * SP];
        if (M && JC[bi] < m) {
          float sm = 0;
          for (size_t j = 0; j < JC[bi]; ++j) sm += W.Ef[bi * SP + j];
          W.Es[bi] = sm;
        }
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
      const float* mk0 = MTp + (size_t)JC[0] * SP;
      const float pi = M ? PIc[JC[0]] : PIc[m];
      float* f = &W.Ff[0];
      const float fu = 0.5f * W.Eu[0], z = fu + 0.5f * pi * W.Es[0];
      const float iz = 1.0f / z;
      fbz = 0.5f * pi * iz;
      fcz = 0.0f;
      for (size_t j = 0; j < SP; ++j) f[j] = fbz * e[j] * MKV(mk0, j);
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
        const float* __restrict mkb = MTp + (size_t)JC[bi] * SP;
        const float* __restrict mmin = MTp + (size_t)std::min(JC[bi], JC[bi - 1]) * SP;
        const float pi = M ? PIc[JC[bi]] : PIc[m];
        float s1 = 0, s2 = 0;
        if (bi >= 2) {
          const float* __restrict fpp = &W.Ff[(bi - 2) * SP];
          for (size_t j = 0; j < SP; ++j) {
            const float ee = e1[j] * e[j] * MKV(mmin, j);
            s1 += ee;
            s2 += fpp[j] * ee;
          }
        } else {
          for (size_t j = 0; j < SP; ++j) s1 += e1[j] * e[j] * MKV(mmin, j);
        }
        const float dot = fbz * s1 + fcz * s2;
        const float fpu = W.Fu[bi - 1];
        const float SA = std::max(0.0f, 1.0f - fpu);  // previous step is normalised
        const float base = fpu * u * pi + stay * h * pi * SA;
        const float c1 = stay * (1.0f - h);
        const float fu = (fpu * stay + SA * u) * W.Eu[bi];
        const float z = fu + base * W.Es[bi] + c1 * dot;
        const float iz = 1.0f / z, bz = base * iz, cz = c1 * iz;
        for (size_t j = 0; j < SP; ++j) fc[j] = (bz + cz * fp[j]) * e[j] * MKV(mkb, j);
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
        const float* __restrict mkb = MTp + (size_t)JC[bi] * SP;
        const float pi = M ? PIc[JC[bi]] : PIc[m];
        const float mcb = M ? (float)JC[bi] : (float)m;
        float t2 = 0;
        if (bi + 1 < nb) {
          const float* __restrict e1 = EP[bi + 1];
          const float* __restrict bn = &W.Bb[(bi + 1) * SP];
          const float* __restrict mmin = MTp + (size_t)std::min(JC[bi], JC[bi + 1]) * SP;
          for (size_t j = 0; j < SP; ++j) t2 += e[j] * e1[j] * bn[j] * MKV(mmin, j);
        }
        const float sw = bbz * W.Es[bi] + bcz * t2;
        const float wu = W.Eu[bi] * bu;
        const float SW = sw * pi;
        const float base = u * wu + stay * h * SW;
        const float c1 = stay * (1.0f - h);
        const float b0 = stay * wu + u * SW;
        const float zz = b0 + mcb * base + c1 * sw;  // any positive scale works
        const float iz = 1.0f / zz, bz = base * iz, cz = c1 * iz;
        for (size_t j = 0; j < SP; ++j) bo[j] = bz + cz * (e[j] * MKV(mkb, j) * bb[j]);
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
      if (!(tot > 0)) {  // degenerate (underflow): treat the bin as unaligned
        W.ab[bi] = 0;
        return;
      }
      double it = 1.0 / (double)tot, al = 1.0 - (double)pu * it;
      double cnt = (double)(W.bfirst[bi + 1] - W.bfirst[bi]);
      if (!fine) W.ab[bi] = al;  // per-bin posteriors are not used after the fine re-run
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
      size_t jc = 0;
      for (size_t bi = lo; bi < hi; ++bi) jc = std::max(jc, jc_at(W, W.bfirst[bi]));
      for (size_t j = 0; j < std::min(jc, SP); ++j) best = std::max(best, acc[j]);
      if (best < (float)P.fast_llr) return false;
      llr_out += best;
      if (hi == nb) break;
    }
    return true;
  }

  // ---- confidence intervals
  static constexpr double CI_ANI_FLOOR = 5e-4, CI_AF_FLOOR = 1e-2;  // SE floors: 0.05 pp ANI, 1 pp (+1/L) AF
  static constexpr int CI_NSTRATA = 24, CI_NLEVELS = 3;  // 4 identity x 3 length x 2 AF classes; 90/95/99%
  std::vector<double> ci_tab;  // [level][stratum][ani, af] conformal multipliers
  static int ci_stratum(double ani, double af, size_t L) {
    const int ia = ani >= 0.99 ? 0 : ani >= 0.97 ? 1 : ani >= 0.93 ? 2 : 3;
    const int il = L < 1000 ? 0 : L < 3000 ? 1 : 2;
    const int ic = af >= 0.999 ? 0 : 1;
    return (ia * 3 + il) * 2 + ic;
  }
  static int ci_level_index(double lv) { return lv < 0.925 ? 0 : lv < 0.97 ? 1 : 2; }
  void load_ci_table() {
    ci_tab.assign((size_t)CI_NLEVELS * CI_NSTRATA * 2, 0.0);
    if (!P.ci_table.empty()) {
      FILE* f = fopen(P.ci_table.c_str(), "r");
      if (!f) throw std::runtime_error("cannot open " + P.ci_table);
      int lv, st;
      double qa, qf;
      char line[256];
      while (fgets(line, sizeof line, f))
        if (sscanf(line, "%d %d %lf %lf", &lv, &st, &qa, &qf) == 4 && lv >= 0 && lv < CI_NLEVELS && st >= 0 && st < CI_NSTRATA) {
          ci_tab[((size_t)lv * CI_NSTRATA + st) * 2] = qa;
          ci_tab[((size_t)lv * CI_NSTRATA + st) * 2 + 1] = qf;
        }
      fclose(f);
      return;
    }
    for (int lv = 0; lv < CI_NLEVELS; ++lv)
      for (int st = 0; st < CI_NSTRATA; ++st) {
        ci_tab[((size_t)lv * CI_NSTRATA + st) * 2] = CI_DEFAULT[lv][st][0];
        ci_tab[((size_t)lv * CI_NSTRATA + st) * 2 + 1] = CI_DEFAULT[lv][st][1];
      }
  }
  // z quantile of the two-sided level (Acklam's rational approximation of the normal quantile)
  static double zq(double level) {
    const double p = 0.5 + 0.5 * level, q = std::sqrt(-2 * std::log(1 - p));
    return q - (2.515517 + 0.802853 * q + 0.010328 * q * q) / (1 + 1.432788 * q + 0.189269 * q * q + 0.001308 * q * q * q);
  }
  // Calibrated intervals (stratified split-conformal): half-width = q * sigma, where sigma is
  // the statistical SE with a small floor and q is the conformal quantile of |error| / sigma
  // in the contig's stratum (estimated identity x length x fully/partially aligned), fitted on
  // the real tuning set against MUMmer truth (bench_real/calibration.py). The ANI statistical
  // SE is set aside in R.ani_se until both estimates are final (finish_ci).
  void set_ani_ci(Result& R, double se_stat) const { R.ani_se = se_stat; }
  void finish_ci(Result& R, size_t L) const {
    const double sa = std::sqrt(R.ani_se * R.ani_se + CI_ANI_FLOOR * CI_ANI_FLOOR);
    const double sf = std::sqrt(R.af_se * R.af_se + std::pow(1.0 / (double)std::max<size_t>(L, 1) + CI_AF_FLOOR, 2));
    double qa, qf;
    if (P.ci_stat || ci_tab.empty()) {
      qa = qf = zq(P.ci_level);
    } else {
      const int st = ci_stratum(R.ani, R.af, L), lv = ci_level_index(P.ci_level);
      qa = ci_tab[((size_t)lv * CI_NSTRATA + st) * 2 + 0];
      qf = ci_tab[((size_t)lv * CI_NSTRATA + st) * 2 + 1];
    }
    const double z = zq(0.95);
    R.ani_lo = std::max(0.0, R.ani - qa * sa);
    R.ani_hi = std::min(1.0, R.ani + qa * sa);
    R.af_lo = std::max(0.0, R.af - qf * sf);
    R.af_hi = std::min(1.0, R.af + qf * sf);
    R.ani_se = std::min(0.5, qa * sa / z);  // reported as the equivalent normal SE of the interval
    R.af_se = std::min(0.5, qf * sf / z);
  }
  // Theorem 1, conditioned on the contig's realised substitutions: probe-level variance of the
  // FP-corrected containment (FMH sampling + false positives) plus the mutation-process term
  // minus the part explained by the realised substitution count, k^2 theta^2 p / ((1-p) n).
  void ani_ci_pooled(Result& R, Work& W, double th, double var, double nk) const {
    const int k = sk.k;
    th = std::min(std::max(th, 1e-9), 1.0);
    const double p = std::max(1e-9, 1 - std::pow(th, 1.0 / k));
    const double cond = std::max(0.0, Bcoef(p, k) - (double)k * k * th * p / (1 - p));
    const double vth = var + (nk > 0 ? th * cond / nk : 0);
    W.se_pooled = R.ani * std::sqrt(vth) / (th * k);
    set_ani_ci(R, W.se_pooled);
  }
  // heterogeneous contigs: the posterior-mean ANI's error is dominated by the spread of local
  // identity, estimated from the heterogeneity-test windows (cluster-robust)
  void ani_ci_het(Result& R, const Work& W) const {
    double sw = 0, s1 = 0, s2 = 0, sw2 = 0;
    for (auto& x : W.hw) {
      const double a = std::pow(std::min(std::max(x.first, 1e-9), 1.0), 1.0 / sk.k);
      sw += x.second;
      s1 += x.second * a;
      s2 += x.second * a * a;
      sw2 += x.second * x.second;
    }
    double se = W.se_pooled;
    if (W.hw.size() >= 2 && sw > 0) {
      const double mu = s1 / sw, var = std::max(0.0, s2 / sw - mu * mu), neff = sw * sw / sw2;
      se = std::max(se, std::sqrt(var / std::max(1.0, neff - 1)));
    }
    set_ani_ci(R, se);
  }
  // AF: the aligned length is uncertain only where the per-probe aligned posterior is
  // intermediate. Each such region is one alignment boundary whose position has survival
  // function a(x) (aligned on the left) or 1 - a(x); a region with the same state on both
  // sides is split at its extremum into two boundaries. Var(AF) = sum of boundary variances / L^2.
  // Variance (bp^2) of the aligned length contributed by uncertain probes in [p0, p1): each run
  // of probes with intermediate aligned posterior is one alignment boundary whose position has
  // survival a(x) (aligned on the left) or 1 - a(x); a run with the same state on both sides is
  // split at its extremum into two boundaries. aleft/aright: posterior just outside the range.
  template <class Post>
  double boundary_var(const Work& W, size_t p0, size_t p1, Post post, double aleft, double aright) const {
    const std::vector<Probe>& pr = W.pr;
    const size_t n = pr.size();
    const double lo = 0.02, hi = 0.98;
    std::vector<double>& xs = const_cast<Work&>(W).cxs;
    std::vector<double>& as = const_cast<Work&>(W).cas;
    auto bvar = [&](size_t s0, size_t e, bool down) {
      double m1 = 0, m2 = 0;
      for (size_t q = s0 + 1; q <= e; ++q) {
        const double dx = xs[q] - xs[q - 1];
        const double S0 = down ? as[q - 1] : 1 - as[q - 1], S1 = down ? as[q] : 1 - as[q];
        const double xm = 0.5 * (xs[q] + xs[q - 1]) - xs[s0];
        m1 += 0.5 * (S0 + S1) * dx;
        m2 += (S0 + S1) * xm * dx;
      }
      return std::max(0.0, m2 - m1 * m1);
    };
    double var = 0, prev = aleft;
    for (size_t i = p0; i < p1;) {
      const double a = post(i);
      if (a <= lo || a >= hi) {
        prev = a;
        ++i;
        continue;
      }
      xs.clear();
      as.clear();
      xs.push_back(i ? (double)pr[i - 1].pos : (double)pr[i].pos);
      as.push_back(prev >= 0.5 ? 1.0 : 0.0);
      while (i < p1) {
        const double ai = post(i);
        if (ai <= lo || ai >= hi) break;
        xs.push_back((double)pr[i].pos);
        as.push_back(ai);
        ++i;
      }
      const double ar = i < p1 ? post(i) : aright;
      xs.push_back(i < n ? (double)pr[std::min(i, n - 1)].pos : (double)pr[n - 1].pos);
      as.push_back(ar >= 0.5 ? 1.0 : 0.0);
      const size_t e = xs.size() - 1;
      const bool lA = as.front() >= 0.5, rA = as.back() >= 0.5;
      if (lA != rA) {
        var += bvar(0, e, lA);
      } else {
        size_t ext = 1;
        for (size_t q = 1; q < e; ++q)
          if (lA ? as[q] < as[ext] : as[q] > as[ext]) ext = q;
        var += bvar(0, ext, lA) + bvar(ext, e, !lA);
      }
      prev = ar;
    }
    return var;
  }

  void af_ci(Result& R, const Work& W, size_t L) const {
    double var = W.afvar;  // refined regions (refine() accumulates it); fast path: 0
    if (!P.refine && !W.refined.empty() && !W.pr.empty()) {  // no refinement: bin-level posteriors
      const size_t nb = W.bstart.size();
      size_t b = 0;
      var = boundary_var(W, 0, W.pr.size(), [&](size_t i) {
        while (b + 1 < nb && W.bfirst[b + 1] <= i) ++b;
        return W.ab[b];
      }, 0.0, 0.0);
    }
    R.af_se = L > 0 ? std::sqrt(var) / (double)L : 0;  // statistical SE; finish_ci calibrates
  }

  Result estimate(std::string_view seq, Work& W) const {
#ifdef AIF_PROF
    const uint64_t _p0 = __rdtsc();
#endif
    probe(seq, W);
#ifdef AIF_PROF
    g_prof[0] += __rdtsc() - _p0;
#endif
    return estimate_probed(seq.size(), W);
  }

  // everything after probing: W.pr holds the contig's probes with this reference's hits
  Result estimate_probed(size_t L, Work& W) const {
    Result R;
    const int k = sk.k;
    W.afvar = 0;
    PSTART;
    const std::vector<Probe>& pr = W.pr;
    R.probes = pr.size();
    if (pr.size() < 3) {
      for (auto& p : pr) R.hits += p.hit;
      R.status = "too_few_probes";
      return R;
    }
    const double llU = make_bins(W, R.hits, P.bin);
    set_local_prior(W);
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
    if (fast) {
      W.refined.assign(nb, 0);
      W.afvar = 0;
    }
    if (R.lr < P.lr_min) {
      R.status = "no_homology";
      return R;
    }
    R.low_frac = low_frac;
    const double* ab = W.ab.data();
    // pass 1 over bins: posterior-weighted tier counts, aligned-start integral, segments
    double N[MAXT] = {0}, Y[MAXT] = {0}, neff = 0;
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
    if (W.want_seg) {
      W.pa.resize(pr.size());
      for (size_t bi = 0; bi < nb; ++bi)
        for (uint32_t i = W.bfirst[bi]; i < W.bfirst[bi + 1]; ++i) W.pa[i] = (float)(W.refined[bi] ? W.ap[i] : ab[bi]);
    }
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
    if (P.het_mode == 0 || P.window <= 0) {
      finish_ci(R, L);
      return R;
    }
    ani_ci_pooled(R, W, th, var, nk_aligned);
    af_ci(R, W, L);
    W.hw.clear();
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
        W.hw.emplace_back(tw, len);
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
      ani_ci_het(R, W);
    }
    if (R.low_frac > P.floor_frac) R.status = "below_floor";
    finish_ci(R, L);
    PSTAMP(6);
    return R;
  }

  // Segments and identity map of a contig with a reported estimate (requires W.want_seg during
  // estimate_probed): maximal runs of probes with aligned posterior > 0.5 (aligned) or <= 0.5
  // (unaligned) tile the contig, boundaries halfway between the flanking probes. Aligned runs are
  // split into windows of about win bases, each with its FP-corrected containment, ANI and a
  // statistical 95% interval (delta method on Theorem 1's variance). Rows: prefix, start, end
  // (0-based, half-open), aligned, ani, ani_lo, ani_hi, probes, hits.
  // split-conformal 95 percent multipliers of the window intervals, fitted on the real tuning set
  // (bench_real/idmap.py calibrate): [window 1 kb, 3 kb][estimated ANI >= 0.99, 0.97-0.99, 0.93-0.97, < 0.93]
  static constexpr double SEG_Q[2][4] = {{2.6187, 3.1523, 4.4475, 3.3788}, {3.3400, 5.0861, 5.8977, 4.5768}};
  void segments(const Work& W, size_t L, int win, const std::string& prefix, std::string& out, bool ci_stat = false) const {
    const std::vector<Probe>& pr = W.pr;
    const size_t m = pr.size();
    if (m == 0 || W.pa.size() != m) return;
    const int k = sk.k;
    char buf[256];
    auto emit = [&](size_t b, size_t e, int al, double ani, double lo, double hi, size_t np, size_t nh) {
      int n = snprintf(buf, sizeof buf, "\t%zu\t%zu\t%d\t%.6f\t%.6f\t%.6f\t%zu\t%zu\n", b, e, al, ani, lo, hi, np, nh);
      out.append(prefix);
      out.append(buf, (size_t)n);
    };
    auto cut = [&](size_t i) -> size_t {  // boundary between probe i-1 and probe i
      if (i == 0) return 0;
      if (i >= m) return L;
      return ((size_t)pr[i - 1].pos + k + (size_t)pr[i].pos) / 2;
    };
    size_t i = 0;
    while (i < m) {
      const bool al = W.pa[i] > 0.5f;
      size_t j = i;
      while (j < m && (W.pa[j] > 0.5f) == al) ++j;
      const size_t b = cut(i), e = cut(j);
      if (!al) {
        size_t nh = 0;
        for (size_t q = i; q < j; ++q) nh += pr[q].hit;
        emit(b, e, 0, NAN, NAN, NAN, j - i, nh);
      } else {  // windows of ~win bases inside the aligned run
        const size_t span = e - b, nwin = std::max<size_t>(1, (span + (size_t)win / 2) / (size_t)std::max(1, win));
        size_t q = i;
        for (size_t w = 0; w < nwin; ++w) {
          const size_t wb = b + span * w / nwin, we = b + span * (w + 1) / nwin;
          double N[MAXT] = {0}, Y[MAXT] = {0};
          size_t np = 0, nh = 0;
          while (q < j && (size_t)pr[q].pos + k <= we + k / 2) {
            N[pr[q].tier] += 1;
            Y[pr[q].tier] += pr[q].hit;
            ++np;
            nh += pr[q].hit;
            ++q;
          }
          if (w + 1 == nwin)
            while (q < j) {
              N[pr[q].tier] += 1;
              Y[pr[q].tier] += pr[q].hit;
              ++np;
              nh += pr[q].hit;
              ++q;
            }
          double ani = NAN, lo = NAN, hi = NAN;
          if (np >= 3) {
            double v = 0;
            const double th0 = combine(N, Y, 0.3);
            const double th = combine(N, Y, std::max(th0, 1e-3), &v);
            if (th > 0) {
              const double nk = std::max(1.0, (double)(we - wb) - (k - 1));
              ani = ani_from_theta(th, v, nk);
              const double se = std::sqrt(std::max(v + th * Bcoef(1 - ani, k) / nk, 0.0)) * ani / (k * th);
              lo = std::max(0.0, ani - 1.959964 * se);
              hi = std::min(1.0, ani + 1.959964 * se);
              if (!ci_stat) {  // calibrated: conformal multiplier on the (clipped) statistical SE
                const int wi = win <= 1732 ? 0 : 1, ci = ani >= 0.99 ? 0 : ani >= 0.97 ? 1 : ani >= 0.93 ? 2 : 3;
                const double se2 = std::max((hi - lo) / (2 * 1.959964), 1e-4);
                lo = std::max(0.0, ani - SEG_Q[wi][ci] * se2);
                hi = std::min(1.0, ani + SEG_Q[wi][ci] * se2);
              }
            } else {
              ani = lo = hi = 0;
            }
          }
          emit(wb, we, 1, ani, lo, hi, np, nh);
        }
      }
      i = j;
    }
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
          "          intervals: [--ci-level 0.95 (0.90/0.95/0.99)] [--ci-stat (statistical SE only)] [--ci-table FILE]\n"
          "          AF model options: [--ident-window D (0: off) --ident-radius 3] (local identity prior)\n"
          "                            [--overdisp L (0: off)] (overdispersed bin emissions, experimental)\n"
          "          [--segments FILE --seg-window 1000 --seg-ci-stat] (aligned/unaligned segments and windowed identity map\n"
          "           per contig, with calibrated 95%% window intervals; --seg-ci-stat: statistical intervals)\n"
          "output: contig ani af theta probes hits aligned_probes status lr het low_frac ani_lo ani_hi af_lo af_hi ani_se af_se\n"
          "aif info <ref.aif>\n"
          "database mode (one k and tier layout for all references; results equal pairwise query with -k K sketches):\n"
          "aif db build -o <db.aifd> [refs.fa ...] [--list refs.txt (path [name] per line)] [-k 20] [-t threads]\n"
          "             [--budget 0.88] [--tiers 0.005:0.6,0.1:0.4] [--max-hashes 6]\n"
          "             [--screen-frac 0.015625] [--screen-bits 8] [--screen-hashes 3] [--screen-group-ratio 1.5] (screening index)\n"
          "aif db search <db.aifd> <contigs.fa> [-t threads] [-n 1 (hits per contig; 0: every candidate, any status)]\n"
          "             [--max-cand 16 (evaluated by the model)] [--prescreen K (0: off; rank K screened references by probe containment)]\n"
          "             [--screen-z 2] [--screen-min 2] [--screen-max-keys 128] [--no-screen (every reference)] [--sort lr|aniaf|ani|af]\n"
          "             [--vmin 5 (skip the model for candidates whose local-homology score is below this; -1: off)]\n"
          "             [query options]\n"
          "output: contig reference rank <query columns> screen_count screen_keys screen_containment candidates\n"
          "aif db info <db.aifd> [--refs]\n");
  exit(1);
}

void format_result(std::string& out, const std::string& name, const Result& r, const char* extra = nullptr) {
  char buf[512];
  int n = snprintf(buf, sizeof buf, "\t%.6f\t%.6f\t%.6f\t%zu\t%zu\t%.1f\t%s\t%.1f\t%d\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                   r.ani, r.af, r.theta, r.probes, r.hits, r.aligned_probes, r.status, r.lr, (int)r.het, r.low_frac,
                   r.ani_lo, r.ani_hi, r.af_lo, r.af_hi, r.ani_se, r.af_se);
  out.append(name);
  if (extra) {
    out.append(buf, (size_t)n - 1);
    out.append(extra);
    out.push_back('\n');
  } else {
    out.append(buf, (size_t)n);
  }
}

// query options shared by "query" and "db search"; returns false if argv[i] is not one
std::string g_seg_path;  // --segments: per-contig segments / identity-map windows (query and db search)
int g_seg_win = 1000;
bool g_seg_ci_stat = false;  // statistical window intervals only (no calibration)
bool parse_query_opt(int& i, int argc, char** argv, Params& P, int& threads) {
  (void)argc;
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
  else if (a == "--ident-window") P.ident_window = atof(argv[++i]);
  else if (a == "--ident-win-probes") P.ident_win_probes = std::max(8, atoi(argv[++i]));
  else if (a == "--ident-radius") P.ident_radius = std::max(0, atoi(argv[++i]));
  else if (a == "--overdisp") P.overdisp = std::max(0.0, atof(argv[++i]));
  else if (a == "--ci-level") P.ci_level = std::min(0.9999, std::max(0.5, atof(argv[++i])));
  else if (a == "--ci-stat") P.ci_stat = true;
  else if (a == "--ci-table") P.ci_table = argv[++i];
  else if (a == "--grid-step") P.grid_step = atof(argv[++i]);
  else if (a == "--floor-frac") P.floor_frac = atof(argv[++i]);
  else if (a == "--no-bias-correct") P.bias_correct = false;
  else if (a == "--segments") g_seg_path = argv[++i];
  else if (a == "--seg-window") g_seg_win = std::max(100, atoi(argv[++i]));
  else if (a == "--seg-ci-stat") g_seg_ci_stat = true;
  else return false;
  return true;
}

// Runs proc(name, seq, Work&, out, out2) over every record of fa and prints the outputs in input
// order (out to stdout, out2 to f2 when given): single-threaded, or the multi-threaded streaming
// pipeline described below.
template <class Proc>
void run_pipeline(Fasta& fa, int threads, FILE* f2, Proc proc) {
  std::string out, out2;
  if (threads == 1) {
    Work W;
    std::string name, seq;
    while (fa.next(name, seq)) {
      proc(name, std::string_view(seq), W, out, out2);
      if (out.size() > (1 << 20)) {
        fwrite(out.data(), 1, out.size(), stdout);
        out.clear();
      }
      if (out2.size() > (1 << 20)) {
        if (f2) fwrite(out2.data(), 1, out2.size(), f2);
        out2.clear();
      }
    }
    fwrite(out.data(), 1, out.size(), stdout);
    if (f2) fwrite(out2.data(), 1, out2.size(), f2);
    prof_report();
    return;
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
    std::string out, out2;
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
            proc(c->names[i], std::string_view(sbuf[t]), works[t], c->out, c->out2);
          }
        } else {
          for (size_t i = 0; i + 1 < c->off.size(); ++i)
            proc(c->names[i], std::string_view(c->data.data() + c->off[i], c->off[i + 1] - c->off[i]), works[t], c->out, c->out2);
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
      if (f2) fwrite(c->out2.data(), 1, c->out2.size(), f2);
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
      c->out2.clear();
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
}

// ---------------------------------------------------------------- database mode
// A database holds the sketches of many references built with ONE k and ONE tier layout
// (budget and tier spec). The tier hash ranges depend only on the layout, so the selected
// k-mers of a contig, their tiers and their filter hashes are the same for every reference:
// a contig is prepared once (packing, selection, block and bit-position hashes; db_prepare),
// and each candidate reference then costs only its filter tests (db_hits) and the HMM.
// Results equal pairwise "aif query" against "aif sketch -k K" sketches of the references.
//
// Candidates come from a screening index, a bit-sliced Bloom filter (as in COBS): k-mers with
// selection hash below h_scr (a fixed fraction of all k-mers, a subset of the probes) are
// inserted with scr_h hash functions into per-reference Bloom columns of a common length.
// References are interleaved in groups of 512 (sorted by size, so a group's columns fill
// alike); one 64-byte row answers a hash function for 512 references. A contig's screening
// count per reference is compared with its null expectation (from the column's fill ratio)
// and the references with the largest z are evaluated with the full model.
//
// File ("AIFD"): 64-byte header, 64-byte aligned blobs (filters, screening rows), metadata,
// trailer (metadata offset and length, magic). The file is memory-mapped; a reference's
// filters are paged in on first use.
const uint32_t DB_MAGIC = 0x44464941;  // "AIFD"
const uint32_t DB_VERSION = 1;
const uint32_t DB_GROUP = 512;
const uint64_t SEED_SCR = 0x94D049BB133111EBULL;
static_assert(sizeof(Probe) == 8 && offsetof(Probe, tier) == 4 && offsetof(Probe, hit) == 5, "Probe layout");

struct DbGroup {
  uint32_t d0 = 0, nd = 0;  // references [d0, d0 + nd)
  uint64_t rows = 0, off = 0;
  size_t words() const { return (nd + 63) / 64; }  // 64-bit words per row
};

// row of screening hash function i for a key (Kirsch-Mitzenmacher double hashing)
inline uint64_t scr_row(uint64_t g, uint32_t i, uint64_t rows) {
  const uint32_t a = (uint32_t)g, b = (uint32_t)(g >> 32) | 1u;
  return ((uint64_t)(uint32_t)(a + i * b) * rows) >> 32;
}

struct Db {
  int k = 20;
  double budget = 0;
  std::vector<std::pair<uint64_t, uint64_t>> tier_rng;
  double scr_frac = 0, scr_bits = 0;
  uint64_t h_scr = 0;
  uint32_t scr_h = 1;
  std::vector<std::string> names;
  std::vector<Sketch> refs;
  std::vector<double> scr_f;     // per reference: screening false-positive rate (fill^scr_h)
  std::vector<uint64_t> scr_n;   // per reference: screening k-mers inserted
  std::vector<DbGroup> groups;
  const uint8_t* map = nullptr;
  size_t map_len = 0;

  Db() = default;
  Db(const Db&) = delete;
  Db& operator=(const Db&) = delete;
  ~Db() {
    if (map) munmap((void*)map, map_len);
  }
  size_t filter_bytes() const {
    size_t b = 0;
    for (auto& r : refs) b += r.filter_bytes();
    return b;
  }
  size_t screen_bytes() const {
    size_t b = 0;
    for (auto& g : groups) b += g.rows * g.words() * 8;
    return b;
  }
  uint64_t total_bases() const {
    uint64_t b = 0;
    for (auto& r : refs) b += r.total_bases;
    return b;
  }
  const uint64_t* rows_of(const DbGroup& g) const { return (const uint64_t*)(map + g.off); }

  void load(const char* path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) throw std::runtime_error(std::string("cannot open database ") + path);
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size < 64 + 24) {
      close(fd);
      throw std::runtime_error("not an aif database (too short)");
    }
    map_len = (size_t)st.st_size;
    void* m = mmap(nullptr, map_len, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) throw std::runtime_error("cannot map database");
    map = (const uint8_t*)m;
    uint32_t magic, ver;
    std::memcpy(&magic, map, 4);
    std::memcpy(&ver, map + 4, 4);
    if (magic != DB_MAGIC) throw std::runtime_error("not an aif database");
    if (ver != DB_VERSION) throw std::runtime_error("unsupported database version");
    uint64_t moff, mlen;
    uint32_t tmagic;
    std::memcpy(&moff, map + map_len - 24, 8);
    std::memcpy(&mlen, map + map_len - 16, 8);
    std::memcpy(&tmagic, map + map_len - 8, 4);
    if (tmagic != DB_MAGIC || moff < 64 || mlen > map_len || moff > map_len - 24 - mlen)
      throw std::runtime_error("truncated or corrupt database (trailer)");
    const uint8_t* p = map + moff;
    const uint8_t* const e = p + mlen;
    auto rd = [&](void* dst, size_t n) {
      if ((size_t)(e - p) < n) throw std::runtime_error("corrupt database (metadata)");
      std::memcpy(dst, p, n);
      p += n;
    };
    auto blob_ok = [&](uint64_t off, uint64_t len) { return off >= 64 && off % 64 == 0 && off <= moff && len <= moff - off; };
    uint32_t nt, nr, ng;
    rd(&k, 4);
    rd(&budget, 8);
    rd(&nt, 4);
    if (k < 5 || k > 31 || nt == 0 || nt > (uint32_t)MAXT) throw std::runtime_error("corrupt database (k or tiers)");
    tier_rng.resize(nt);
    for (auto& r : tier_rng) {
      rd(&r.first, 8);
      rd(&r.second, 8);
    }
    rd(&scr_frac, 8);
    rd(&scr_bits, 8);
    rd(&h_scr, 8);
    rd(&scr_h, 4);
    if (scr_h == 0 || scr_h > 8) throw std::runtime_error("corrupt database (screening)");
    rd(&nr, 4);
    names.resize(nr);
    refs.resize(nr);
    scr_f.resize(nr);
    scr_n.resize(nr);
    for (uint32_t d = 0; d < nr; ++d) {
      uint32_t len;
      rd(&len, 4);
      if (len > (1u << 20)) throw std::runtime_error("corrupt database (name)");
      names[d].resize(len);
      rd(&names[d][0], len);
      Sketch& sk = refs[d];
      sk.k = k;
      sk.budget = budget;
      rd(&sk.total_bases, 8);
      rd(&sk.n_kmers, 8);
      rd(&sk.sigma_emp, 8);
      rd(&scr_f[d], 8);
      rd(&scr_n[d], 8);
      sk.tiers.resize(nt);
      for (uint32_t t = 0; t < nt; ++t) {
        Tier& T = sk.tiers[t];
        T.lo = tier_rng[t].first;
        T.hi = tier_rng[t].second;
        uint64_t off;
        rd(&T.nh, 4);
        rd(&T.nblocks, 8);
        rd(&T.n_ins, 8);
        rd(&T.eps, 8);
        rd(&off, 8);
        if (T.nh == 0 || T.nh > 30 || T.nblocks == 0 || T.nblocks >= (1ULL << 32) || !blob_ok(off, T.nblocks * 64))
          throw std::runtime_error("corrupt database (filter of " + names[d] + ")");
        T.bp = (const uint64_t*)(map + off);
      }
    }
    rd(&ng, 4);
    groups.resize(ng);
    uint32_t next = 0;
    for (auto& g : groups) {
      rd(&g.d0, 4);
      rd(&g.nd, 4);
      rd(&g.rows, 8);
      rd(&g.off, 8);
      if (g.d0 != next || g.nd == 0 || g.nd > DB_GROUP || g.rows == 0 || g.rows >= (1ULL << 32) || !blob_ok(g.off, g.rows * g.words() * 8))
        throw std::runtime_error("corrupt database (screening group)");
      next += g.nd;
    }
    if (next != nr) throw std::runtime_error("corrupt database (groups do not cover the references)");
    // the screening rows are read for every contig: ask for them up front
    for (auto& g : groups) madvise((void*)(map + g.off), g.rows * g.words() * 8, MADV_WILLNEED);
  }
};

int cmd_db_build(int argc, char** argv) {
  std::string out, list;
  int k = 20, threads = 1;
  double budget = 0.88, scr_frac = 1.0 / 64, scr_bits = 8, grp_ratio = 1.5;
  uint32_t max_hashes = 6, scr_h = 3;
  std::string tiers = "0.005:0.6,0.1:0.4";
  struct RefIn {
    std::string path, name;
    uint64_t fsize = 0;
  };
  std::vector<RefIn> in;
  for (int i = 3; i < argc; ++i) {
    std::string a = argv[i];
    auto val = [&]() -> const char* {
      if (i + 1 >= argc) usage();
      return argv[++i];
    };
    if (a == "-o") out = val();
    else if (a == "-k") k = atoi(val());
    else if (a == "-t") threads = std::max(1, atoi(val()));
    else if (a == "--list") list = val();
    else if (a == "--budget") budget = atof(val());
    else if (a == "--tiers") tiers = val();
    else if (a == "--max-hashes") max_hashes = (uint32_t)std::max(1, atoi(val()));
    else if (a == "--screen-frac") scr_frac = atof(val());
    else if (a == "--screen-bits") scr_bits = atof(val());
    else if (a == "--screen-hashes") scr_h = (uint32_t)std::max(1, std::min(8, atoi(val())));
    else if (a == "--screen-group-ratio") grp_ratio = std::max(1.0, atof(val()));
    else if (!a.empty() && a[0] == '-') usage();
    else in.push_back({a, "", 0});
  }
  if (!list.empty()) {
    FILE* f = fopen(list.c_str(), "r");
    if (!f) throw std::runtime_error("cannot read " + list);
    char line[8192];
    while (fgets(line, sizeof line, f)) {
      std::string l(line);
      while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
      if (l.empty() || l[0] == '#') continue;
      size_t tab = l.find_first_of("\t ");
      if (tab == std::string::npos) in.push_back({l, "", 0});
      else in.push_back({l.substr(0, tab), l.substr(l.find_first_not_of("\t ", tab)), 0});
    }
    fclose(f);
  }
  if (out.empty() || in.empty() || k < 5 || k > 31) usage();
  auto spec = parse_tiers(tiers);
  {  // sizes (in parallel: on a cold disk the per-file reads dominate) and default names
    std::atomic<size_t> nx{0};
    std::vector<std::string> err(in.size());
    auto sizer = [&]() {
      for (size_t i; (i = nx.fetch_add(1)) < in.size();) {
        RefIn& r = in[i];
        struct stat st;
        if (stat(r.path.c_str(), &st) != 0) {
          err[i] = "cannot stat " + r.path;
          continue;
        }
        r.fsize = (uint64_t)st.st_size;
        // gzip: the uncompressed size from the trailer (ISIZE, exact below 4 GiB), else 4x the file
        int fd = open(r.path.c_str(), O_RDONLY);
        unsigned char m[2] = {0, 0}, z[4];
        if (fd >= 0 && st.st_size >= 18 && pread(fd, m, 2, 0) == 2 && m[0] == 0x1f && m[1] == 0x8b) {
          uint64_t isz = 0;
          if (pread(fd, z, 4, st.st_size - 4) == 4) isz = (uint64_t)z[0] | ((uint64_t)z[1] << 8) | ((uint64_t)z[2] << 16) | ((uint64_t)z[3] << 24);
          r.fsize = isz > r.fsize ? isz : 4 * r.fsize;
        }
        if (fd >= 0) close(fd);
        if (r.name.empty()) {
          std::string bn = r.path.substr(r.path.find_last_of('/') + 1);
          for (const char* ext : {".gz", ".fasta", ".fna", ".fa", ".fas", ".ffn"})
            if (bn.size() > strlen(ext) && bn.compare(bn.size() - strlen(ext), strlen(ext), ext) == 0) bn.resize(bn.size() - strlen(ext));
          r.name = bn;
        }
      }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < std::max(threads, 16); ++t) pool.emplace_back(sizer);
    for (auto& th : pool) th.join();
    for (auto& e : err)
      if (!e.empty()) throw std::runtime_error(e);
  }
  // references sorted by size: a screening group's columns then have similar fill
  std::stable_sort(in.begin(), in.end(), [](const RefIn& a, const RefIn& b) { return a.fsize < b.fsize; });
  const double t0 = now();
  // the tier layout (hash ranges) is the same for every reference: take it from a tiny sketch
  std::vector<std::pair<uint64_t, uint64_t>> rng;
  {
    std::vector<Packed> none;
    Sketch probe = build_sketch(none, k, budget, spec, 0, max_hashes);
    for (auto& t : probe.tiers) rng.emplace_back(t.lo, t.hi);
  }
  const long double TWO64 = 18446744073709551616.0L;
  const uint64_t hmax = rng.back().second;
  uint64_t h_scr = (uint64_t)std::min<long double>((long double)scr_frac * TWO64, (long double)hmax);
  if (scr_frac <= 0 || scr_bits <= 0) throw std::runtime_error("--screen-frac and --screen-bits must be positive");
  FILE* f = fopen(out.c_str(), "wb");
  if (!f) throw std::runtime_error("cannot write " + out);
  uint64_t pos = 0;
  auto put = [&](const void* p, size_t n) {
    if (n && fwrite(p, 1, n, f) != n) throw std::runtime_error("write failed (disk full?)");
    pos += n;
  };
  auto align64 = [&]() {
    static const char z[64] = {0};
    if (pos % 64) put(z, 64 - pos % 64);
  };
  {
    char hdr[64] = {0};
    std::memcpy(hdr, &DB_MAGIC, 4);
    std::memcpy(hdr + 4, &DB_VERSION, 4);
    put(hdr, 64);
  }
  struct RefMeta {
    uint64_t total_bases, n_kmers;
    double sigma_emp, scr_f;
    uint64_t scr_n;
    std::vector<uint32_t> nh;
    std::vector<uint64_t> nblocks, n_ins, off;
    std::vector<double> eps;
  };
  std::vector<RefMeta> meta(in.size());
  std::vector<DbGroup> groups;
  uint64_t fbytes = 0, sbytes = 0;
  // screening groups: consecutive (size-sorted) references, at most DB_GROUP per group and at most
  // grp_ratio between the largest and smallest file, so the common column length fits them all
  for (size_t g0 = 0, nd = 0; g0 < in.size(); g0 += nd) {
    nd = 1;
    while (g0 + nd < in.size() && nd < DB_GROUP && (double)in[g0 + nd].fsize <= grp_ratio * (double)std::max<uint64_t>(1, in[g0].fsize)) ++nd;
    double meansize = 0;
    for (size_t d = 0; d < nd; ++d) meansize += (double)in[g0 + d].fsize / (double)nd;
    // rows from the group's mean file size (file bytes ~ bases); larger members get a somewhat
    // higher screening false-positive rate, which is measured per reference (column fill)
    const uint64_t rows = std::max<uint64_t>(64, (uint64_t)std::ceil(scr_bits * scr_frac * meansize));
    if (rows >= (1ULL << 32)) throw std::runtime_error("reference too large for the screening index");
    const size_t words = (nd + 63) / 64;
    std::vector<uint64_t> mat(rows * words, 0);
    std::vector<Sketch> gs(nd);
    std::vector<uint64_t> gn(nd, 0);
    std::atomic<size_t> next{0};
    std::vector<std::string> errs(nd);
    auto work = [&]() {
      const size_t CH = 1 << 16;
      std::vector<uint64_t> sx(CH + 8), shv(CH + 8);
      std::vector<uint32_t> spos(CH + 40);
      for (;;) {
        const size_t d = next.fetch_add(1);
        if (d >= nd) return;
        try {
          std::vector<Packed> packs;
          {
            Fasta fa(in[g0 + d].path.c_str());
            std::string name, seq;
            while (fa.next(name, seq)) {
              packs.emplace_back();
              pack_seq((const uint8_t*)seq.data(), seq.size(), k, packs.back());
            }
          }
          gs[d] = build_sketch(packs, k, budget, spec, 1000000, max_hashes);  // also sets badk for k
          uint64_t n = 0;
          for (auto& pk : packs) {
            if (pk.L < (size_t)k) continue;
            const size_t nst = pk.L - k + 1;
            for (size_t c0 = 0; c0 < nst; c0 += CH) {
              const size_t m = select_kmers(pk, k, h_scr, c0, c0 + CH, sx.data(), shv.data(), spos.data());
              for (size_t j = 0; j < m; ++j) {
                const uint64_t g = fmix64(fmix64(sx[j] ^ SEED_BF) ^ SEED_SCR);
                for (uint32_t i = 0; i < scr_h; ++i)
                  __atomic_fetch_or(&mat[scr_row(g, i, rows) * words + (d >> 6)], 1ULL << (d & 63), __ATOMIC_RELAXED);
              }
              n += m;
            }
          }
          gn[d] = n;
        } catch (std::exception& e) {
          errs[d] = e.what();
        }
      }
    };
    std::vector<std::thread> pool;
    for (int t = 1; t < threads; ++t) pool.emplace_back(work);
    work();
    for (auto& th : pool) th.join();
    for (size_t d = 0; d < nd; ++d)
      if (!errs[d].empty()) throw std::runtime_error(in[g0 + d].path + ": " + errs[d]);
    // column fill -> screening false-positive rate of each reference
    std::vector<uint64_t> ones(DB_GROUP, 0);
    for (uint64_t r = 0; r < rows; ++r)
      for (size_t w = 0; w < words; ++w)
        for (uint64_t b = mat[r * words + w]; b; b &= b - 1) ++ones[w * 64 + __builtin_ctzll(b)];
    for (size_t d = 0; d < nd; ++d) {
      const Sketch& sk = gs[d];
      RefMeta& M = meta[g0 + d];
      M.total_bases = sk.total_bases;
      M.n_kmers = sk.n_kmers;
      M.sigma_emp = sk.sigma_emp;
      M.scr_n = gn[d];
      M.scr_f = std::pow((double)ones[d] / (double)rows, (double)scr_h);
      for (auto& t : sk.tiers) {
        if (t.lo != rng[&t - &sk.tiers[0]].first || t.hi != rng[&t - &sk.tiers[0]].second)
          throw std::runtime_error("internal: tier layout differs between references");
        align64();
        M.nh.push_back(t.nh);
        M.nblocks.push_back(t.nblocks);
        M.n_ins.push_back(t.n_ins);
        M.eps.push_back(t.eps);
        M.off.push_back(pos);
        put(t.bits.data(), t.bits.size() * 8);
        fbytes += t.bits.size() * 8;
      }
      gs[d] = Sketch();
    }
    align64();
    groups.push_back({(uint32_t)g0, (uint32_t)nd, rows, pos});
    put(mat.data(), mat.size() * 8);
    sbytes += mat.size() * 8;
    fprintf(stderr, "\r%zu / %zu references", g0 + nd, in.size());
  }
  fprintf(stderr, "\n");
  std::string md;
  auto w = [&](const void* p, size_t n) { md.append((const char*)p, n); };
  const uint32_t nt = (uint32_t)rng.size(), nr = (uint32_t)in.size(), ng = (uint32_t)groups.size();
  w(&k, 4);
  w(&budget, 8);
  w(&nt, 4);
  for (auto& r : rng) {
    w(&r.first, 8);
    w(&r.second, 8);
  }
  w(&scr_frac, 8);
  w(&scr_bits, 8);
  w(&h_scr, 8);
  w(&scr_h, 4);
  w(&nr, 4);
  uint64_t tb = 0;
  for (size_t d = 0; d < in.size(); ++d) {
    const RefMeta& M = meta[d];
    const uint32_t len = (uint32_t)in[d].name.size();
    w(&len, 4);
    w(in[d].name.data(), len);
    w(&M.total_bases, 8);
    w(&M.n_kmers, 8);
    w(&M.sigma_emp, 8);
    w(&M.scr_f, 8);
    w(&M.scr_n, 8);
    for (uint32_t t = 0; t < nt; ++t) {
      w(&M.nh[t], 4);
      w(&M.nblocks[t], 8);
      w(&M.n_ins[t], 8);
      w(&M.eps[t], 8);
      w(&M.off[t], 8);
    }
    tb += M.total_bases;
  }
  w(&ng, 4);
  for (auto& g : groups) {
    w(&g.d0, 4);
    w(&g.nd, 4);
    w(&g.rows, 8);
    w(&g.off, 8);
  }
  const uint64_t moff = pos, mlen = md.size();
  put(md.data(), md.size());
  char tr[24] = {0};
  std::memcpy(tr, &moff, 8);
  std::memcpy(tr + 8, &mlen, 8);
  std::memcpy(tr + 16, &DB_MAGIC, 4);
  put(tr, 24);
  if (fclose(f) != 0) throw std::runtime_error("write failed");
  fprintf(stderr, "references=%u bases=%lu k=%d filters=%.1f MB (%.4f bits/base) screening=%.1f MB (%.4f bits/base) time=%.1f s\n",
          nr, (unsigned long)tb, k, fbytes / 1e6, 8.0 * fbytes / std::max<double>(1, tb), sbytes / 1e6,
          8.0 * sbytes / std::max<double>(1, tb), now() - t0);
  return 0;
}

// Reference-independent probe preparation for a database search: positions and tiers into
// W.pr (hits cleared), the block hash (high 32 bits) and bit-position hash of every probe into
// W.pbh / W.ph2, and the screening keys (probes with selection hash < h_scr) into W.skey.
void db_prepare_window(const Db& D, const uint8_t* s, size_t L, uint32_t off, Work& W) {
  const int k = D.k;
  const size_t nst = L - k + 1, base = W.pr.size();
  const int T = (int)D.tier_rng.size();
  const uint64_t hmax = D.tier_rng.back().second;
  pack_seq(s, L, k, W.pk);
#ifdef AIF_AVX512
  if (T <= 2) {
    if (W.kpos.size() < nst + 32) W.kpos.resize(nst + 32);
    const size_t m = select_positions(W.pk, k, hmax, W.kpos.data());
    W.pr.resize(base + m);
    W.pbh.resize(base + m);
    W.ph2.resize(base + m);
    size_t ns = W.skey.size();
    W.skey.resize(ns + m);
    W.skh.resize(ns + m);
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
    const __m512i vhi0 = _mm512_set1_epi64((long long)D.tier_rng[0].second), vscr = _mm512_set1_epi64((long long)D.h_scr);
    const __m512i sbf = _mm512_set1_epi64((long long)SEED_BF), sbit = _mm512_set1_epi64((long long)SEED_BIT);
    const __m512i zero = _mm512_setzero_si512(), one = _mm512_set1_epi64(1), lane = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
    const __m512i tierb = _mm512_set1_epi64(1LL << 32), voff = _mm512_set1_epi64((long long)off);
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
      const __m512i hs = fmix64_v(_mm512_xor_si512(x, vsel));
      const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_cmpge_epu64_mask(hs, vhi0) & tail) : (__mmask8)0;
      const __m512i bh = fmix64_v(_mm512_xor_si512(x, sbf));
      const __m512i h2 = fmix64_v(_mm512_xor_si512(bh, sbit));
      _mm256_mask_storeu_epi32(W.pbh.data() + base + j, tail, _mm512_cvtepi64_epi32(_mm512_srli_epi64(bh, 32)));
      _mm512_mask_storeu_epi64(W.ph2.data() + base + j, tail, h2);
      const __mmask8 ms = _mm512_cmplt_epu64_mask(hs, vscr) & tail;
      _mm512_mask_compressstoreu_epi64(W.skey.data() + ns, ms, bh);
      _mm512_mask_compressstoreu_epi64(W.skh.data() + ns, ms, hs);
      ns += (size_t)__builtin_popcount((unsigned)ms);
      // Probe {pos + off, tier, hit = 0} as one 64-bit lane
      const __m512i pv = _mm512_mask_or_epi64(_mm512_add_epi64(pos, voff), m1, _mm512_add_epi64(pos, voff), tierb);
      _mm512_mask_storeu_epi64((void*)(W.pr.data() + base + j), tail, pv);
      (void)lane;
    }
    W.skey.resize(ns);
    W.skh.resize(ns);
    return;
  }
#endif
  if (W.canon.size() < nst + 8) {
    W.canon.resize(nst + 8);
    W.kpos.resize(nst + 8);
    W.hs.resize(nst + 8);
  }
  const size_t m = select_kmers(W.pk, k, hmax, 0, nst, W.canon.data(), W.hs.data(), W.kpos.data());
  W.pr.resize(base + m);
  W.pbh.resize(base + m);
  W.ph2.resize(base + m);
  for (size_t j = 0; j < m; ++j) {
    const uint64_t hs = W.hs[j];
    int t = 0;
    while (t + 1 < T && hs >= D.tier_rng[t].second) ++t;
    const uint64_t bh = fmix64(W.canon[j] ^ SEED_BF);
    W.pbh[base + j] = (uint32_t)(bh >> 32);
    W.ph2[base + j] = fmix64(bh ^ SEED_BIT);
    W.pr[base + j] = {W.kpos[j] + off, (uint8_t)t, 0};
    if (hs < D.h_scr) {
      W.skey.push_back(bh);
      W.skh.push_back(hs);
    }
  }
}

void db_prepare(const Db& D, std::string_view seq, Work& W) {
  W.pr.clear();
  W.pbh.clear();
  W.ph2.clear();
  W.skey.clear();
  W.skh.clear();
  const size_t L = seq.size();
  if (L < (size_t)D.k) return;
  const size_t nst = L - D.k + 1;
  for (size_t off = 0; off < nst; off += Estimator::PROBE_WIN) {
    const size_t n = std::min(Estimator::PROBE_WIN, nst - off);
    db_prepare_window(D, (const uint8_t*)seq.data() + off, n + D.k - 1, (uint32_t)off, W);
  }
}

// filter tests of the prepared probes against one reference: sets W.pr[i].hit; adds the hits of
// each tier to hits[] when given
void db_hits(const Sketch& R, Work& W, uint64_t* hits = nullptr) {
  const size_t m = W.pr.size();
  const int T = (int)R.tiers.size();
  uint64_t hc[MAXT] = {0};
#ifdef AIF_AVX512
  if (T <= 2) {
    const Tier& t0 = R.tiers[0];
    const Tier& t1 = R.tiers[T - 1];
    const long long* base0 = (const long long*)t0.data();
    const long long delta1 = (long long)((const long long*)t1.data() - base0);
    const __m512i vnb0 = _mm512_set1_epi64((long long)t0.nblocks), vnb1 = _mm512_set1_epi64((long long)t1.nblocks);
    const __m512i voff1 = _mm512_set1_epi64(delta1), zero = _mm512_setzero_si512();
    const __m512i m511 = _mm512_set1_epi64(511), m63 = _mm512_set1_epi64(63), one = _mm512_set1_epi64(1);
    const __m512i tierm = _mm512_set1_epi64(0xFFLL << 32), hitclr = _mm512_set1_epi64(~(0xFFLL << 40)), hitb = _mm512_set1_epi64(1LL << 40);
    const uint32_t nh0 = t0.nh, nh1 = t1.nh, nhmax = std::max(nh0, nh1);
    long long* prw = (long long*)W.pr.data();
    for (size_t j = 0; j < m; j += 8) {
      const __mmask8 tail = (m - j >= 8) ? (__mmask8)0xFF : (__mmask8)((1u << (m - j)) - 1);
      __m512i pv = _mm512_maskz_loadu_epi64(tail, prw + j);
      const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_test_epi64_mask(pv, tierm) & tail) : (__mmask8)0;
      const __m512i bh32 = _mm512_cvtepu32_epi64(_mm256_maskz_loadu_epi32(tail, W.pbh.data() + j));
      const __m512i nb = _mm512_mask_blend_epi64(m1, vnb0, vnb1);
      const __m512i blk = _mm512_srli_epi64(_mm512_mul_epu32(bh32, nb), 32);
      const __m512i wb = _mm512_add_epi64(_mm512_slli_epi64(blk, 3), _mm512_mask_blend_epi64(m1, zero, voff1));
      __m512i h2 = _mm512_maskz_loadu_epi64(tail, W.ph2.data() + j);
      __m512i acc = one;
      for (uint32_t q = 0; q < nhmax; ++q) {
        if (q && q % 7 == 0) h2 = fmix64_v(_mm512_add_epi64(h2, _mm512_set1_epi64(q)));
        const __m512i bp = _mm512_and_si512(_mm512_srl_epi64(h2, _mm_cvtsi32_si128(9 * (q % 7))), m511);
        const __m512i idx = _mm512_add_epi64(wb, _mm512_srli_epi64(bp, 6));
        const __mmask8 am = (__mmask8)(((q < nh0) ? (__mmask8)(~m1) : 0) | ((q < nh1) ? m1 : 0)) & tail;
        const __m512i g = _mm512_mask_i64gather_epi64(zero, am, idx, base0, 8);
        acc = _mm512_mask_and_epi64(acc, am, acc, _mm512_srlv_epi64(g, _mm512_and_si512(bp, m63)));
      }
      const __mmask8 hit = _mm512_test_epi64_mask(acc, one) & tail;
      pv = _mm512_or_si512(_mm512_and_si512(pv, hitclr), _mm512_maskz_mov_epi64(hit, hitb));
      _mm512_mask_storeu_epi64(prw + j, tail, pv);
      hc[T - 1] += (uint64_t)__builtin_popcount((unsigned)(hit & m1));
      hc[0] += (uint64_t)__builtin_popcount((unsigned)(hit & ~m1));
    }
    if (hits)
      for (int t = 0; t < T; ++t) hits[t] = hc[t];
    return;
  }
#endif
  for (size_t j = 0; j < m; ++j) {
    const Tier& t = R.tiers[W.pr[j].tier];
    const uint64_t* blk = t.data() + (((uint64_t)W.pbh[j] * t.nblocks) >> 32) * 8;
    W.pr[j].hit = t.contains_h2(blk, W.ph2[j]);
    hc[W.pr[j].tier] += W.pr[j].hit;
  }
  if (hits)
    for (int t = 0; t < T; ++t) hits[t] = hc[t];
}

// hits per tier of the prepared probes against one reference, without storing them (candidate
// ranking); returns the FP-corrected containment of the probes
double db_containment(const Sketch& R, const Work& W, const size_t* ntier) {
  const size_t m = W.pr.size();
  const int T = (int)R.tiers.size();
  uint64_t h[MAXT] = {0};
#ifdef AIF_AVX512
  if (T <= 2) {
    const Tier& t0 = R.tiers[0];
    const Tier& t1 = R.tiers[T - 1];
    const long long* base0 = (const long long*)t0.data();
    const long long delta1 = (long long)((const long long*)t1.data() - base0);
    const __m512i vnb0 = _mm512_set1_epi64((long long)t0.nblocks), vnb1 = _mm512_set1_epi64((long long)t1.nblocks);
    const __m512i voff1 = _mm512_set1_epi64(delta1), zero = _mm512_setzero_si512();
    const __m512i m511 = _mm512_set1_epi64(511), m63 = _mm512_set1_epi64(63), one = _mm512_set1_epi64(1);
    const __m512i tierm = _mm512_set1_epi64(0xFFLL << 32);
    const uint32_t nh0 = t0.nh, nh1 = t1.nh, nhmax = std::max(nh0, nh1);
    const long long* prw = (const long long*)W.pr.data();
    for (size_t j = 0; j < m; j += 8) {
      const __mmask8 tail = (m - j >= 8) ? (__mmask8)0xFF : (__mmask8)((1u << (m - j)) - 1);
      const __m512i pv = _mm512_maskz_loadu_epi64(tail, prw + j);
      const __mmask8 m1 = (T == 2) ? (__mmask8)(_mm512_test_epi64_mask(pv, tierm) & tail) : (__mmask8)0;
      const __m512i bh32 = _mm512_cvtepu32_epi64(_mm256_maskz_loadu_epi32(tail, W.pbh.data() + j));
      const __m512i nb = _mm512_mask_blend_epi64(m1, vnb0, vnb1);
      const __m512i blk = _mm512_srli_epi64(_mm512_mul_epu32(bh32, nb), 32);
      const __m512i wb = _mm512_add_epi64(_mm512_slli_epi64(blk, 3), _mm512_mask_blend_epi64(m1, zero, voff1));
      __m512i h2 = _mm512_maskz_loadu_epi64(tail, W.ph2.data() + j);
      __m512i acc = one;
      for (uint32_t q = 0; q < nhmax; ++q) {
        if (q && q % 7 == 0) h2 = fmix64_v(_mm512_add_epi64(h2, _mm512_set1_epi64(q)));
        const __m512i bp = _mm512_and_si512(_mm512_srl_epi64(h2, _mm_cvtsi32_si128(9 * (q % 7))), m511);
        const __m512i idx = _mm512_add_epi64(wb, _mm512_srli_epi64(bp, 6));
        const __mmask8 am = (__mmask8)(((q < nh0) ? (__mmask8)(~m1) : 0) | ((q < nh1) ? m1 : 0)) & tail;
        const __m512i g = _mm512_mask_i64gather_epi64(zero, am, idx, base0, 8);
        acc = _mm512_mask_and_epi64(acc, am, acc, _mm512_srlv_epi64(g, _mm512_and_si512(bp, m63)));
      }
      const __mmask8 hit = _mm512_test_epi64_mask(acc, one) & tail;
      h[T - 1] += (uint64_t)__builtin_popcount((unsigned)(hit & m1));
      h[0] += (uint64_t)__builtin_popcount((unsigned)(hit & ~m1));
    }
  } else
#endif
  {
    for (size_t j = 0; j < m; ++j) {
      const Tier& t = R.tiers[W.pr[j].tier];
      const uint64_t* blk = t.data() + (((uint64_t)W.pbh[j] * t.nblocks) >> 32) * 8;
      h[W.pr[j].tier] += t.contains_h2(blk, W.ph2[j]);
    }
  }
  const double sg = R.sigma_used();
  double num = 0, den = 0;
  for (int t = 0; t < T; ++t) {
    const double e = R.tiers[t].eps + sg - R.tiers[t].eps * sg;
    num += (double)h[t] - (double)ntier[t] * e;
    den += (double)ntier[t] * (1 - e);
  }
  return den > 0 ? num / den : 0.0;
}

// screening counts of the contig's keys for every reference into W.scnt16 (one count per
// reference). Contigs with more than max_keys keys use the max_keys keys of smallest selection
// hash (a consistent subsample: containment estimates stay unbiased; long contigs have signal to
// spare). Returns the number of keys used.
size_t db_screen(const Db& D, Work& W, size_t max_keys) {
  const size_t nr = D.refs.size();
  size_t ns = W.skey.size();
  if (ns > max_keys) {
    std::vector<std::pair<uint64_t, uint64_t>> kh(ns);
    for (size_t i = 0; i < ns; ++i) kh[i] = {W.skh[i], W.skey[i]};
    std::nth_element(kh.begin(), kh.begin() + (max_keys - 1), kh.end());
    for (size_t i = 0; i < max_keys; ++i) W.skey[i] = kh[i].second;
    ns = max_keys;
  }
  W.scnt16.assign(nr + 64, 0);
  if (!ns) return 0;
  const uint32_t H = D.scr_h;
  for (const DbGroup& g : D.groups) {
    const uint64_t* rows = D.rows_of(g);
    const size_t words = g.words();
    uint16_t* c16 = W.scnt16.data() + g.d0;
    const size_t PF = 8;
    for (size_t i = 0; i < std::min(ns, PF); ++i) {
      const uint64_t gk = fmix64(W.skey[i] ^ SEED_SCR);
      for (uint32_t h = 0; h < H; ++h) __builtin_prefetch(rows + scr_row(gk, h, g.rows) * words);
    }
    for (size_t i = 0; i < ns; ++i) {
      if (i + PF < ns) {
        const uint64_t gk = fmix64(W.skey[i + PF] ^ SEED_SCR);
        for (uint32_t h = 0; h < H; ++h) __builtin_prefetch(rows + scr_row(gk, h, g.rows) * words);
      }
      const uint64_t gk = fmix64(W.skey[i] ^ SEED_SCR);
      uint64_t b[8];
      const uint64_t* r0 = rows + scr_row(gk, 0, g.rows) * words;
      for (size_t q = 0; q < words; ++q) b[q] = r0[q];
      for (uint32_t h = 1; h < H; ++h) {
        const uint64_t* rh = rows + scr_row(gk, h, g.rows) * words;
        for (size_t q = 0; q < words; ++q) b[q] &= rh[q];
      }
#ifdef AIF_AVX512
      const __m512i one16 = _mm512_set1_epi16(1);
      for (size_t q = 0; q < words; ++q) {
        if (!b[q]) continue;
        __m512i lo = _mm512_loadu_si512((const void*)(c16 + q * 64));
        __m512i hi = _mm512_loadu_si512((const void*)(c16 + q * 64 + 32));
        lo = _mm512_mask_add_epi16(lo, (__mmask32)(uint32_t)b[q], lo, one16);
        hi = _mm512_mask_add_epi16(hi, (__mmask32)(uint32_t)(b[q] >> 32), hi, one16);
        _mm512_storeu_si512((void*)(c16 + q * 64), lo);
        _mm512_storeu_si512((void*)(c16 + q * 64 + 32), hi);
      }
#else
      for (size_t q = 0; q < words; ++q)
        for (uint64_t x = b[q]; x; x &= x - 1) ++c16[q * 64 + __builtin_ctzll(x)];
#endif
    }
  }
  return ns;
}

// Local-homology score of the prepared probes against one reference (hits set by db_hits): the
// best path of a two-state (unaligned / aligned at identity a) model, maximised over a few
// identity levels, as a log-likelihood ratio against "all unaligned". Probe scores are the HMM's
// emission log-ratios; switches cost the HMM's switching probability over the probe gap; starting
// and ending aligned is free. A cheap stand-in for the HMM's likelihood ratio, used to skip the
// model for candidates that cannot reach it (calibrated empirically, see db_bench.py).
// per-probe switching cost (log switch probability over the gap to the previous probe, relative
// to staying); reference-independent, computed once per contig
void db_switch_costs(Work& W, double rho_u) {
  W.psw.resize(W.pr.size());
  uint32_t prev = W.pr.empty() ? 0 : W.pr[0].pos;
  for (size_t i = 0; i < W.pr.size(); ++i) {
    const double d = (double)(W.pr[i].pos - prev) + 1.0;
    prev = W.pr[i].pos;
    W.psw[i] = (float)(std::log(-std::expm1(-rho_u * d)) + rho_u * d);
  }
}

double db_vscore(const Sketch& R, const Work& W) {
  static const double LEV[8] = {0.995, 0.98, 0.96, 0.93, 0.90, 0.87, 0.84, 0.81};
  const int T = (int)R.tiers.size(), k = R.k;
  const double sg = R.sigma_used();
  alignas(32) float S[2 * MAXT][8];  // score per probe kind (2 * tier + hit) and level
  for (int j = 0; j < 8; ++j) {
    const double th = std::pow(LEV[j], k);
    for (int t = 0; t < T; ++t) {
      const double e = R.tiers[t].eps + sg - R.tiers[t].eps * sg;
      S[2 * t][j] = (float)std::log(1 - th);
      S[2 * t + 1][j] = (float)std::log((th + (1 - th) * e) / e);
    }
  }
  alignas(32) float vu[8] = {0}, va[8] = {0};
  const Probe* pr = W.pr.data();
  const float* sw = W.psw.data();
  for (size_t i = 0, m = W.pr.size(); i < m; ++i) {
    const float* sc = S[2 * pr[i].tier + pr[i].hit];
    const float c = sw[i];
    for (int j = 0; j < 8; ++j) {
      const float nu = std::max(vu[j], va[j] + c);
      const float na = std::max(va[j], vu[j] + c) + sc[j];
      vu[j] = nu;
      va[j] = na;
    }
  }
  float best = 0;
  for (int j = 0; j < 8; ++j) best = std::max(best, std::max(vu[j], va[j]));
  return best;
}

// Lazily constructed estimators, one per reference (shared by all threads; construction is
// race-free: a losing thread discards its copy).
struct DbEstimators {
  const Db& D;
  Params P;
  std::unique_ptr<std::atomic<Estimator*>[]> slot;
  std::atomic<size_t> built{0};
  DbEstimators(const Db& d, const Params& p) : D(d), P(p), slot(new std::atomic<Estimator*>[d.refs.size()]) {
    for (size_t i = 0; i < D.refs.size(); ++i) slot[i].store(nullptr);
    if (!D.refs.empty()) get(0);  // prototype for the shared distance tables
  }
  ~DbEstimators() {
    for (size_t i = 0; i < D.refs.size(); ++i) delete slot[i].load();
  }
  const Estimator& get(uint32_t d) {
    Estimator* e = slot[d].load(std::memory_order_acquire);
    if (e) return *e;
    Estimator* fresh = new Estimator(D.refs[d], P, d ? slot[0].load(std::memory_order_acquire) : nullptr);
    Estimator* expected = nullptr;
    if (slot[d].compare_exchange_strong(expected, fresh, std::memory_order_acq_rel)) {
      built.fetch_add(1);
      return *fresh;
    }
    delete fresh;
    return *expected;
  }
};

int cmd_db_search(int argc, char** argv) {
  if (argc < 5) usage();
  Db D;
  D.load(argv[3]);
  Params P;
  // hit order: 3 = the model's log-likelihood ratio (default: how well a reference explains the
  // contig's probe hits; robust to AF over-estimates on divergent hits), 0 = ani*af, 1 = ani, 2 = af
  int threads = 1, nout = 1, sort_by = 3;
  size_t max_cand = 16, prescreen = 0, max_keys = 128;  // screening keys per contig (smallest selection hashes)
  // pre-HMM local-score filter: on the 2,005-genome benchmark the lowest score of any reported
  // hit was 9.1 (370,757 hits), and a cutoff of 5 skipped 99% of the candidates that report nothing
  double vmin = 5.0;
  bool vreport = false;
  double zmin = 2.0;
  uint32_t cmin = 2;
  bool screen = true;
  for (int i = 5; i < argc; ++i) {
    std::string a = argv[i];
    auto val = [&]() -> const char* {
      if (i + 1 >= argc) usage();
      return argv[++i];
    };
    if (a == "-n") nout = std::max(0, atoi(val()));
    else if (a == "--max-cand") max_cand = (size_t)std::max(1, atoi(val()));
    else if (a == "--prescreen") prescreen = (size_t)std::max(0, atoi(val()));
    else if (a == "--screen-max-keys") max_keys = (size_t)std::min(65535, std::max(16, atoi(val())));
    else if (a == "--vmin") vmin = atof(val());
    else if (a == "--vreport") vreport = true;
    else if (a == "--screen-z") zmin = atof(val());
    else if (a == "--screen-min") cmin = (uint32_t)std::max(0, atoi(val()));
    else if (a == "--no-screen") screen = false;
    else if (a == "--sort") {
      std::string v = val();
      sort_by = v == "ani" ? 1 : v == "af" ? 2 : v == "aniaf" ? 0 : v == "lr" ? 3 : -1;
      if (sort_by < 0) usage();
    } else if (!parse_query_opt(i, argc, argv, P, threads)) usage();
  }
#if defined(__x86_64__)
  _mm_setcsr(_mm_getcsr() | 0x8040);
#endif
  DbEstimators ests(D, P);
  Fasta fa(argv[4]);
  fputs("contig\treference\trank\tani\taf\ttheta\tprobes\thits\taligned_probes\tstatus\tlr\thet\tlow_frac\tani_lo\tani_hi"
        "\taf_lo\taf_hi\tani_se\taf_se\tscreen_count\tscreen_keys\tscreen_containment\tcandidates\n",
        stdout);
  const uint32_t nr = (uint32_t)D.refs.size();
  struct Hit {
    Result R;
    uint32_t d;
    float z;
    double score;
  };
  FILE* fseg = nullptr;
  if (!g_seg_path.empty()) {
    fseg = fopen(g_seg_path.c_str(), "w");
    if (!fseg) throw std::runtime_error("cannot write " + g_seg_path);
    fputs("contig\treference\tstart\tend\taligned\tani\tani_lo\tani_hi\tprobes\thits\n", fseg);
  }
  run_pipeline(fa, threads, fseg, [&](const std::string& name, std::string_view seq, Work& W, std::string& out, std::string& out2) {
    thread_local std::vector<Hit> hits;
    W.want_seg = false;
    hits.clear();
    db_prepare(D, seq, W);
    size_t ns = W.skey.size();
    W.cand.clear();
    if (screen) {
      ns = db_screen(D, W, max_keys);
      const uint16_t* cnt = W.scnt16.data();
      const uint16_t c0 = (uint16_t)std::max<uint32_t>(1, cmin);
      for (uint32_t d0 = 0; d0 < nr; d0 += 32) {
#ifdef AIF_AVX512
        uint32_t msk = (uint32_t)_mm512_cmpge_epu16_mask(_mm512_loadu_si512((const void*)(cnt + d0)), _mm512_set1_epi16((short)c0));
        if (nr - d0 < 32) msk &= (1u << (nr - d0)) - 1;
#else
        uint32_t msk = 0;
        for (uint32_t q = 0; q < 32 && d0 + q < nr; ++q) msk |= (uint32_t)(cnt[d0 + q] >= c0) << q;
#endif
        for (; msk; msk &= msk - 1) {
        const uint32_t d = d0 + (uint32_t)__builtin_ctz(msk);
        const uint32_t c = cnt[d];
        if (c < cmin || c == 0) continue;
        const double f = D.scr_f[d], mu = (double)ns * f, sd = std::sqrt(std::max((double)ns * f * (1 - f), 1e-12));
        const double z = ((double)c - mu) / sd;
        // significance by z; ranking by the FP-corrected screening containment (z would favour
        // references with a low screening FPR on long contigs)
        if (z >= zmin) W.cand.emplace_back((float)(((double)c / (double)ns - f) / (1 - f)), d);
        }
      }
      auto better = [](const std::pair<float, uint32_t>& a, const std::pair<float, uint32_t>& b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;
      };
      const size_t K1 = std::max(max_cand, prescreen);
      if (W.cand.size() > K1) {
        std::partial_sort(W.cand.begin(), W.cand.begin() + K1, W.cand.end(), better);
        W.cand.resize(K1);
      } else {
        std::sort(W.cand.begin(), W.cand.end(), better);
      }
      // second stage: rank the screened references by the FP-corrected containment of ALL probes
      // (filter tests only), and keep the best max_cand for the full model
      if (W.cand.size() > max_cand) {
        size_t nt[MAXT] = {0};
        for (const Probe& p : W.pr) ++nt[p.tier];
        W.cand2.clear();
        for (auto& cd : W.cand) W.cand2.emplace_back((float)db_containment(D.refs[cd.second], W, nt), cd.second);
        std::partial_sort(W.cand2.begin(), W.cand2.begin() + max_cand, W.cand2.end(), better);
        std::unordered_map<uint32_t, float> zof;
        W.cand2.resize(max_cand);
        for (auto& cd : W.cand) zof[cd.second] = cd.first;
        W.cand.clear();
        for (auto& cd : W.cand2) W.cand.emplace_back(zof[cd.second], cd.second);
      }
    } else {
      W.scnt16.assign(nr + 64, 0);
      for (uint32_t d = 0; d < nr; ++d) W.cand.emplace_back(0.0f, d);
    }
    size_t nskip = 0;
    bool have_sw = false;
    size_t ntier[MAXT] = {0};
    for (const Probe& p : W.pr) ++ntier[p.tier];
    for (auto& cd : W.cand) {
      const uint32_t d = cd.second;
      uint64_t hc[MAXT] = {0};
      db_hits(D.refs[d], W, hc);
      double vs = 0;
      // the local score is computed only for candidates without a clear global hit excess (z < 4)
      bool weak = vreport;
      if (vmin > 0 && !weak) {
        const Sketch& R = D.refs[d];
        const double sg = R.sigma_used();
        double ex = 0, var = 0;
        for (size_t t = 0; t < R.tiers.size(); ++t) {
          const double e = R.tiers[t].eps + sg - R.tiers[t].eps * sg;
          ex += (double)hc[t] - (double)ntier[t] * e;
          var += (double)ntier[t] * e * (1 - e);
        }
        weak = ex < 4.0 * std::sqrt(std::max(var, 1e-9));
      }
      if (weak) {
        if (!have_sw) {
          db_switch_costs(W, P.rho_u);
          have_sw = true;
        }
        vs = db_vscore(D.refs[d], W);
        if (vs < vmin) {  // cannot reach the model's reporting threshold: skip the HMM
          ++nskip;
          if (vreport) {
            Result R;
            R.probes = W.pr.size();
            R.status = "skipped";
            char ex[160];
            snprintf(ex, sizeof ex, "\t%u\t%zu\t%.4f\t%zu\t%.2f", W.scnt16[d], ns, cd.first, W.cand.size(), vs);
            format_result(out, name + "\t" + D.names[d] + "\t0", R, ex);
          }
          continue;
        }
      }
      const Estimator& E = ests.get(d);
      Result R = E.estimate_probed(seq.size(), W);
      if (vreport) {
        char ex[160];
        snprintf(ex, sizeof ex, "\t%u\t%zu\t%.4f\t%zu\t%.2f", W.scnt16[d], ns, cd.first, W.cand.size(), vs);
        format_result(out, name + "\t" + D.names[d] + "\t0", R, ex);
        continue;
      }
      const bool ok = std::strcmp(R.status, "ok") == 0;
      if (!ok && nout > 0) continue;
      const double sc = !ok ? -1e300 : sort_by == 1 ? R.ani : sort_by == 2 ? R.af : sort_by == 3 ? R.lr : R.ani * R.af;
      hits.push_back({R, d, cd.first, sc});
    }
    (void)nskip;
    if (vreport) return;
    if (nout > 0) {
      std::stable_sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        return a.score != b.score ? a.score > b.score : a.d < b.d;
      });
      if (hits.size() > (size_t)nout) hits.resize(nout);
    }
    char extra[128];
    if (hits.empty()) {
      Result R;
      R.probes = W.pr.size();
      R.status = W.cand.empty() ? "no_candidate" : "no_homology";
      snprintf(extra, sizeof extra, "\t0\t%zu\tnan\t%zu", ns, W.cand.size());
      format_result(out, name + "\t*\t0", R, extra);
      return;
    }
    if (fseg && !hits.empty() && std::strcmp(hits[0].R.status, "ok") == 0) {  // segments of the top hit
      const uint32_t d = hits[0].d;
      db_hits(D.refs[d], W);
      W.want_seg = true;
      const Estimator& E = ests.get(d);
      E.estimate_probed(seq.size(), W);
      E.segments(W, seq.size(), g_seg_win, name + "\t" + D.names[d], out2, g_seg_ci_stat);
      W.want_seg = false;
    }
    for (size_t r = 0; r < hits.size(); ++r) {
      const Hit& h = hits[r];
      snprintf(extra, sizeof extra, "\t%u\t%zu\t%.4f\t%zu", (unsigned)W.scnt16[h.d], ns, screen ? h.z : NAN, W.cand.size());
      format_result(out, name + "\t" + D.names[h.d] + "\t" + std::to_string(r + 1), h.R, extra);
    }
  });
  if (fseg) fclose(fseg);
  fprintf(stderr, "estimators built: %zu of %u references\n", ests.built.load(), nr);
  return 0;
}

int cmd_db_info(int argc, char** argv) {
  if (argc < 4) usage();
  Db D;
  D.load(argv[3]);
  const bool refs = argc > 4 && std::string(argv[4]) == "--refs";
  const uint64_t tb = D.total_bases();
  printf("references=%zu bases=%lu k=%d budget=%.4f tiers=%zu filter_bytes=%zu (%.4f bits/base) screen_bytes=%zu (%.4f bits/base)"
         " screen_frac=%.5f screen_bits=%.2f screen_hashes=%u groups=%zu file_bytes=%zu\n",
         D.refs.size(), (unsigned long)tb, D.k, D.budget, D.tier_rng.size(), D.filter_bytes(),
         8.0 * D.filter_bytes() / std::max<double>(1, tb), D.screen_bytes(), 8.0 * D.screen_bytes() / std::max<double>(1, tb),
         D.scr_frac, D.scr_bits, D.scr_h, D.groups.size(), D.map_len);
  if (refs) {
    printf("reference\tbases\tsigma_emp\tscreen_kmers\tscreen_fpr");
    for (size_t t = 0; t < D.tier_rng.size(); ++t) printf("\tnh%zu\teps%zu", t, t);
    printf("\n");
    for (size_t d = 0; d < D.refs.size(); ++d) {
      const Sketch& sk = D.refs[d];
      printf("%s\t%lu\t%.3g\t%lu\t%.4f", D.names[d].c_str(), (unsigned long)sk.total_bases, sk.sigma_emp,
             (unsigned long)D.scr_n[d], D.scr_f[d]);
      for (auto& t : sk.tiers) printf("\t%u\t%.5f", t.nh, t.eps);
      printf("\n");
    }
  }
  return 0;
}

int cmd_query(int argc, char** argv) {
  if (argc < 4) usage();
  Sketch sk;
  sk.load(argv[2]);
  Params P;
  int threads = 1;
  for (int i = 4; i < argc; ++i)
    if (!parse_query_opt(i, argc, argv, P, threads)) usage();
  Estimator est(sk, P);
#if defined(__x86_64__)
  _mm_setcsr(_mm_getcsr() | 0x8040);  // flush-to-zero / denormals-are-zero
#endif
  Fasta fa(argv[3]);
  fputs("contig\tani\taf\ttheta\tprobes\thits\taligned_probes\tstatus\tlr\thet\tlow_frac\tani_lo\tani_hi\taf_lo\taf_hi\tani_se\taf_se\n",
        stdout);
  FILE* fseg = nullptr;
  if (!g_seg_path.empty()) {
    fseg = fopen(g_seg_path.c_str(), "w");
    if (!fseg) throw std::runtime_error("cannot write " + g_seg_path);
    fputs("contig\tstart\tend\taligned\tani\tani_lo\tani_hi\tprobes\thits\n", fseg);
  }
  run_pipeline(fa, threads, fseg, [&](const std::string& name, std::string_view seq, Work& W, std::string& out, std::string& out2) {
    W.want_seg = fseg != nullptr;
    const Result R = est.estimate(seq, W);
    format_result(out, name, R);
    if (fseg && std::strcmp(R.status, "ok") == 0) est.segments(W, seq.size(), g_seg_win, name, out2, g_seg_ci_stat);
  });
  if (fseg) fclose(fseg);
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
    } else if (cmd == "db") {
      std::string sub = argv[2];
      if (sub == "build") return cmd_db_build(argc, argv);
      if (sub == "search") return cmd_db_search(argc, argv);
      if (sub == "info") return cmd_db_info(argc, argv);
      usage();
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
