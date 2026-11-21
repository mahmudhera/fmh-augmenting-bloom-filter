// aif: FracMinHash + Bloom-filter tier sketches for ANI / alignment-fraction
// estimation of (short) contigs against a sketched reference.
//
// Design (see fmh_bloom_filter_theory.pdf, Sections 3-7):
//   * reference side: canonical k-mers are FracMinHash-selected into disjoint
//     hash-range tiers; each tier is a blocked Bloom filter with its own FPR
//     (hedged tiers, D8). Total size is set by a bits-per-base budget.
//   * query side (raw contig): every selected k-mer is probed (D2). A two-state
//     HMM (aligned / unaligned) over the probe hits, fitted by EM, yields the
//     aligned containment theta (FP- and spurious-match corrected, Theorem 1)
//     and the posterior aligned fraction (D12, Proposition 5).
//   * ANI = theta^(1/k) with the second-order bias correction (Corollary 1);
//     optional window-wise estimation reduces heterogeneity bias (Prop. 9).
//
// Build: g++ -O3 -std=c++17 -o aif aif.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const uint64_t SEED_SEL = 0x9E3779B97F4A7C15ULL;
const uint64_t SEED_BF = 0xC2B2AE3D27D4EB4FULL;
const uint64_t SEED_BIT = 0x165667B19E3779F9ULL;
const uint32_t MAGIC = 0x32464941;  // "AIF2"

inline uint64_t fmix64(uint64_t k) {
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53ULL;
  k ^= k >> 33;
  return k;
}
inline uint64_t mulhi(uint64_t a, uint64_t b) {
  return (uint64_t)(((unsigned __int128)a * b) >> 64);
}

int8_t CODE[256];
double T_PROBE = 0, T_HMM = 0;
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

// ---------------------------------------------------------------- FASTA
struct Fasta {
  FILE* f = nullptr;
  char* line = nullptr;
  size_t cap = 0;
  std::string pending;
  bool have_pending = false;
  explicit Fasta(const char* path) {
    f = fopen(path, "r");
    if (!f) throw std::runtime_error(std::string("cannot open ") + path);
  }
  ~Fasta() {
    if (f) fclose(f);
    free(line);
  }
  bool next(std::string& name, std::string& seq) {
    seq.clear();
    if (!have_pending) {
      ssize_t n;
      while ((n = getline(&line, &cap, f)) > 0) {
        if (line[0] == '>') {
          pending.assign(line + 1, n - 1);
          have_pending = true;
          break;
        }
      }
      if (!have_pending) return false;
    }
    size_t e = pending.find_first_of(" \t\r\n");
    name = pending.substr(0, e);
    have_pending = false;
    ssize_t n;
    while ((n = getline(&line, &cap, f)) > 0) {
      if (line[0] == '>') {
        pending.assign(line + 1, n - 1);
        have_pending = true;
        break;
      }
      while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r')) --n;
      seq.append(line, n);
    }
    return true;
  }
};

inline uint64_t revcomp_kmer(uint64_t x, int k) {
  x = ~x;  // complement (A<->T, C<->G under 0..3 coding)
  x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
  x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
  x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
  x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
  x = (x >> 32) | (x << 32);
  return x >> (64 - 2 * k);
}

template <class F>
void for_each_kmer(const std::string& s, int k, F&& f) {
  const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const int shift = 2 * (k - 1);
  uint64_t fw = 0, rv = 0;
  int len = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    int c = CODE[(uint8_t)s[i]];
    if (c < 0) {
      len = 0;
      fw = rv = 0;
      continue;
    }
    fw = ((fw << 2) | (uint64_t)c) & mask;
    rv = (rv >> 2) | ((uint64_t)(3 - c) << shift);
    if (++len >= k) f((uint32_t)(i + 1 - k), fw < rv ? fw : rv);
  }
}

// ---------------------------------------------------------------- Bloom tier
struct Tier {
  uint64_t lo = 0, hi = 0;  // selection-hash range [lo, hi)
  uint32_t nh = 1;          // hash functions
  uint64_t nblocks = 1;     // 512-bit blocks
  uint64_t n_ins = 0;       // insertions (with multiplicity)
  double eps = 0;           // exact FPR of a random non-member
  std::vector<uint64_t> bits;

  double s() const { return (double)(hi - lo) / 18446744073709551616.0; }
  inline bool selects(uint64_t hs) const { return hs >= lo && hs < hi; }
  inline void positions(uint64_t x, uint64_t& block, uint32_t* bitpos) const {
    uint64_t h = fmix64(x ^ SEED_BF);
    block = mulhi(h, nblocks);
    uint64_t h2 = fmix64(h ^ SEED_BIT);
    for (uint32_t j = 0; j < nh; ++j) {
      if (j && j % 7 == 0) h2 = fmix64(h2 + j);
      bitpos[j] = (uint32_t)((h2 >> (9 * (j % 7))) & 511);
    }
  }
  void insert(uint64_t x) {
    uint64_t b;
    uint32_t p[64];
    positions(x, b, p);
    uint64_t* blk = &bits[b * 8];
    for (uint32_t j = 0; j < nh; ++j) blk[p[j] >> 6] |= 1ULL << (p[j] & 63);
  }
  inline const uint64_t* block_of(uint64_t x, uint64_t& h) const {
    h = fmix64(x ^ SEED_BF);
    return &bits[mulhi(h, nblocks) * 8];
  }
  // early-exit membership test given the block pointer and its hash
  inline bool contains_at(const uint64_t* blk, uint64_t h) const {
    uint64_t h2 = fmix64(h ^ SEED_BIT);
    for (uint32_t j = 0; j < nh; ++j) {
      if (j && j % 7 == 0) h2 = fmix64(h2 + j);
      uint32_t bp = (uint32_t)((h2 >> (9 * (j % 7))) & 511);
      if (!(blk[bp >> 6] >> (bp & 63) & 1)) return false;
    }
    return true;
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
  uint64_t n_kmers = 0;   // k-mer positions in the reference (~ distinct k-mers)
  double budget = 0;      // bits per base requested
  double sigma_emp = -1;  // empirically calibrated spurious-match rate (<0: unknown)
  std::vector<Tier> tiers;

  double sigma() const {  // analytic spurious rate of a uniform random canonical k-mer
    return -std::expm1(-2.0 * (double)n_kmers / std::pow(4.0, k));
  }
  double sigma_used() const { return sigma_emp >= 0 ? sigma_emp : sigma(); }
  size_t filter_bytes() const {
    size_t b = 0;
    for (auto& t : tiers) b += t.bits.size() * 8;
    return b;
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
  return v;
}

std::vector<std::string> load_fasta(const char* fasta) {
  std::vector<std::string> seqs;
  std::string name, seq;
  Fasta fa(fasta);
  while (fa.next(name, seq)) seqs.push_back(seq);
  return seqs;
}

Sketch build_sketch(const std::vector<std::string>& seqs, int k, double budget, const std::vector<TierSpec>& spec) {
  Sketch sk;
  sk.k = k;
  sk.budget = budget;
  for (auto& sq : seqs) {
    sk.total_bases += sq.size();
    if (sq.size() >= (size_t)k) sk.n_kmers += sq.size() - k + 1;
  }
  double fsum = 0;
  for (auto& t : spec) fsum += t.frac;
  // sampling rates from the integer-h Bloom cost at each target FPR
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
  }
  // size each tier from its expected number of selections (s * k-mer positions)
  for (auto& t : sk.tiers) {
    double m = (double)t.nblocks * 512.0, n = std::max<double>(1.0, t.s() * (double)sk.n_kmers);
    double best = 2;
    for (uint32_t h = 1; h <= 40; ++h) {
      double hh = (double)h;
      double fpr = std::pow(1 - std::exp(-hh * n / m), hh);
      if (fpr < best) {
        best = fpr;
        t.nh = h;
      }
    }
  }
  // single pass: insert selected k-mers
  for (auto& sq : seqs)
    for_each_kmer(sq, k, [&](uint32_t, uint64_t x) {
      uint64_t hs = fmix64(x ^ SEED_SEL);
      for (auto& t : sk.tiers)
        if (t.selects(hs)) {
          t.insert(x);
          ++t.n_ins;
          break;
        }
    });
  for (auto& t : sk.tiers) t.compute_eps();
  // empirical spurious-match calibration: probe single-substitution neighbours of
  // reference k-mers (captures composition bias and near-repeats)
  {
    uint64_t rng = 0x2545F4914F6CDD1DULL;
    auto next = [&]() {
      rng ^= rng << 13;
      rng ^= rng >> 7;
      rng ^= rng << 17;
      return rng;
    };
    const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    std::vector<double> N(sk.tiers.size(), 0), Y(sk.tiers.size(), 0);
    size_t target = 400000, tries = 0;
    std::vector<size_t> lens;
    size_t totlen = 0;
    for (auto& sq : seqs) {
      lens.push_back(sq.size());
      totlen += sq.size();
    }
    for (size_t it = 0; it < target && tries < 4 * target; ++tries) {
      size_t r = next() % totlen, si = 0;
      while (r >= lens[si]) r -= lens[si++];
      const std::string& sq = seqs[si];
      if (r + k > sq.size()) continue;
      uint64_t fw = 0;
      bool ok = true;
      for (int j = 0; j < k; ++j) {
        int c = CODE[(uint8_t)sq[r + j]];
        if (c < 0) { ok = false; break; }
        fw = (fw << 2) | (uint64_t)c;
      }
      if (!ok) continue;
      ++it;
      int off = (int)(next() % k);
      int sh = 2 * (k - 1 - off);
      uint64_t old = (fw >> sh) & 3, nw = (old + 1 + next() % 3) & 3;
      fw = (fw & ~(3ULL << sh)) | (nw << sh);
      fw &= mask;
      uint64_t rv = revcomp_kmer(fw, k);
      uint64_t x = fw < rv ? fw : rv;
      uint64_t hs = fmix64(x ^ SEED_SEL);
      for (size_t t = 0; t < sk.tiers.size(); ++t)
        if (sk.tiers[t].selects(hs)) {
          N[t] += 1;
          Y[t] += sk.tiers[t].contains(x);
          break;
        }
    }
    double num = 0, den = 0;
    for (size_t t = 0; t < sk.tiers.size(); ++t) {
      num += Y[t] - N[t] * sk.tiers[t].eps;
      den += N[t] * (1 - sk.tiers[t].eps);
    }
    sk.sigma_emp = den > 0 ? std::max(0.0, num / den) : -1;
  }
  return sk;
}

// ---------------------------------------------------------------- estimation
struct Probe {
  uint32_t pos;
  uint8_t tier;
  uint8_t hit;
};

struct Params {
  double rho_u = 3e-6;     // per-base switch rate aligned <-> unaligned
  double rho_h = 1e-4;     // per-base switch rate between aligned identity levels
  int window = 1000;       // window (bases) of the heterogeneity test
  int het_mode = 1;        // 0 = pooled, 1 = auto (test -> posterior mean), 2 = always posterior mean
  double het_z = 3.0;      // z threshold of the heterogeneity test
  double lr_min = 10.0;    // log-likelihood ratio (vs. all-unaligned) needed to report
  double ident_c = 8.0;    // keep aligned levels with theta >= ident_c * min effective FPR
  double grid_step = 0.005;  // ANI spacing of aligned levels
  double floor_frac = 2.0;   // flag below_floor if this much aligned mass sits on the 2 lowest levels (>1 = off)
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

// process-variance coefficient B(p,k) of Blanca et al. (leading term)
double Bcoef(double p, int k) {
  if (p <= 1e-9) return 0;
  double q = 1 - std::pow(1 - p, k);
  return std::max(0.0, q * (2.0 * k - 1 + 2.0 / p) - 2.0 * k);
}

struct Estimator {
  const Sketch& sk;
  Params P;
  std::vector<double> epsp;   // effective FP rate per tier (FP + spurious)
  std::vector<double> stier;  // sampling rate per tier
  std::vector<double> gth;    // aligned-state containment levels
  std::vector<double> gani;   // their ANI values
  std::vector<double> etab;   // emission table
  std::vector<double> logU;   // log U-emission per (tier, hit)
  static const size_t QMAX = 1 << 16;
  std::vector<double> tab_su, tab_h;  // exp(-rho_u d), 1-exp(-rho_h d)

  Estimator(const Sketch& s, const Params& p) : sk(s), P(p) {
    double sg = sk.sigma_used();
    for (auto& t : sk.tiers) {
      epsp.push_back(t.eps + sg - t.eps * sg);
      stier.push_back(t.s());
    }
    // aligned identity levels: identifiable ones only (theta well above the noise floor)
    double emin = 1;
    for (double e : epsp) emin = std::min(emin, e);
    double amin = std::pow(std::min(0.9, P.ident_c * emin), 1.0 / sk.k);
    amin = std::max(amin, 0.70);
    for (double a = 0.9995; a >= amin - 1e-12; a -= P.grid_step) {
      gani.push_back(a);
      gth.push_back(std::min(std::max(std::pow(a, sk.k), 1e-7), 0.9995));
    }
    // emission table per (tier, hit): [U, A_1..A_m]
    {
      size_t S = gth.size() + 1;
      etab.assign(sk.tiers.size() * 2 * S, 0);
      for (size_t t = 0; t < sk.tiers.size(); ++t)
        for (int hit = 0; hit < 2; ++hit) {
          double* dst = &etab[(t * 2 + hit) * S];
          double e = epsp[t];
          dst[0] = hit ? e : 1 - e;
          for (size_t j = 0; j < gth.size(); ++j) {
            double t1 = gth[j] + (1 - gth[j]) * e;
            dst[1 + j] = hit ? t1 : 1 - t1;
          }
        }
    }
    for (size_t t = 0; t < sk.tiers.size(); ++t)
      for (int hit = 0; hit < 2; ++hit) logU.push_back(std::log(hit ? epsp[t] : 1 - epsp[t]));
    tab_su.resize(QMAX);
    tab_h.resize(QMAX);
    for (size_t d = 0; d < QMAX; ++d) {
      tab_su[d] = std::exp(-P.rho_u * (double)d);
      tab_h[d] = -std::expm1(-P.rho_h * (double)d);
    }
  }
  inline double su(uint32_t d) const { return d < QMAX ? tab_su[d] : std::exp(-P.rho_u * d); }
  inline double hh(uint32_t d) const { return d < QMAX ? tab_h[d] : -std::expm1(-P.rho_h * d); }

  // tier-combined, inverse-variance weighted FP-corrected containment
  double combine(const double* N, const double* Y, double theta_ref, double* var_out = nullptr) const {
    double wsum = 0, acc = 0;
    double th = std::min(std::max(theta_ref, 1e-4), 1.0);
    for (size_t t = 0; t < sk.tiers.size(); ++t) {
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

  double theta_range(const std::vector<Probe>& pr, const std::vector<double>& w, double theta_ref, size_t lo,
                     size_t hi, double* var = nullptr, double* neff = nullptr) const {
    double N[8] = {0}, Y[8] = {0};
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

  // forward-backward over states {U, A_1..A_m}; returns aligned posterior and log-likelihoods
  void hmm(const std::vector<Probe>& pr, std::vector<double>& aligned, double& ll, double& llU,
           double* pm_ani = nullptr, double* low_frac = nullptr) const {
    const size_t n = pr.size(), m = gth.size(), S = m + 1;
    std::vector<double> F(n * S);
    const double pi = 1.0 / m;
    ll = 0;
    llU = 0;
    double zprod = 1.0;
    int zcnt = 0;
    auto row = [&](const Probe& p) { return &etab[(p.tier * 2 + p.hit) * S]; };
    {
      const double* e = row(pr[0]);
      double* f = &F[0];
      f[0] = 0.5 * e[0];
      double z = f[0];
      for (size_t j = 1; j < S; ++j) {
        f[j] = 0.5 * pi * e[j];
        z += f[j];
      }
      double iz = 1 / z;
      for (size_t j = 0; j < S; ++j) f[j] *= iz;
      ll += std::log(z);
      llU += logU[pr[0].tier * 2 + pr[0].hit];
    }
    for (size_t i = 1; i < n; ++i) {
      uint32_t d = pr[i].pos - pr[i - 1].pos;
      const double stay = su(d), u = 1 - stay, h = hh(d);
      const double* fp = &F[(i - 1) * S];
      double* __restrict fc = &F[i * S];
      const double* __restrict e = row(pr[i]);
      double SA = 0;
      for (size_t j = 1; j < S; ++j) SA += fp[j];
      const double base = fp[0] * u * pi + stay * h * pi * SA;
      const double c1 = stay * (1 - h);
      double z0 = (fp[0] * stay + SA * u) * e[0];
      fc[0] = z0;
      double z = 0;
      for (size_t j = 1; j < S; ++j) {
        double v = (base + c1 * fp[j]) * e[j];
        fc[j] = v;
        z += v;
      }
      z += z0;
      const double iz = 1 / z;
      for (size_t j = 0; j < S; ++j) fc[j] *= iz;
      zprod *= z;
      if (++zcnt == 32) {
        ll += std::log(zprod);
        zprod = 1.0;
        zcnt = 0;
      }
      llU += logU[pr[i].tier * 2 + pr[i].hit];
    }
    ll += std::log(zprod);
    aligned.assign(n, 0);
    std::vector<double> b(S, 1.0), w(S);
    double pm_num = 0, pm_den = 0, low = 0;
    const size_t nlow = std::min<size_t>(2, m);
    {
      const double* fl = &F[(n - 1) * S];
      aligned[n - 1] = 1 - fl[0];
      for (size_t j = 0; j < m; ++j) {
        pm_num += fl[1 + j] * gani[j];
        pm_den += fl[1 + j];
      }
      for (size_t j = m - nlow; j < m; ++j) low += fl[1 + j];
    }
    for (size_t i = n - 1; i > 0; --i) {
      uint32_t d = pr[i].pos - pr[i - 1].pos;
      const double stay = su(d), u = 1 - stay, h = hh(d);
      const double* __restrict e = row(pr[i]);
      double* __restrict bb = b.data();
      double* __restrict ww = w.data();
      ww[0] = e[0] * bb[0];
      double SW = 0;
      for (size_t j = 1; j < S; ++j) {
        ww[j] = e[j] * bb[j];
        SW += ww[j];
      }
      SW *= pi;
      const double base = u * ww[0] + stay * h * SW;
      const double c1 = stay * (1 - h);
      double b0 = stay * ww[0] + u * SW;
      double zz = b0;
      for (size_t j = 1; j < S; ++j) {
        double v = base + c1 * ww[j];
        bb[j] = v;
        zz += v;
      }
      bb[0] = b0;
      const double iz = 1 / zz;
      const double* fp = &F[(i - 1) * S];
      double pu = fp[0] * bb[0] * iz, tot = pu, num = 0;
      bb[0] *= iz;
      for (size_t j = 1; j < S; ++j) {
        bb[j] *= iz;
        double g = fp[j] * bb[j];
        tot += g;
        num += g * gani[j - 1];
      }
      double itot = 1 / tot;
      aligned[i - 1] = 1 - pu * itot;
      pm_num += num * itot;
      pm_den += 1 - pu * itot;
      double lw = 0;
      for (size_t j = S - nlow; j < S; ++j) lw += fp[j] * bb[j];
      low += lw * itot;
    }
    if (pm_ani) *pm_ani = pm_den > 0 ? pm_num / pm_den : NAN;
    if (low_frac) *low_frac = pm_den > 0 ? low / pm_den : 0;
  }

  Result estimate(const std::string& seq) const {
    Result R;
    const int k = sk.k;
    double t0 = now();
    std::vector<Probe> pr;
    pr.reserve(seq.size() / 2 + 8);
    // phase 1: select k-mers and prefetch their filter blocks; phase 2: test membership
    const uint64_t hmax = sk.tiers.back().hi;
    struct Sel {
      const uint64_t* blk;
      uint64_t h;
    };
    std::vector<Sel> sel;
    sel.reserve(seq.size() / 4 + 8);
    for_each_kmer(seq, k, [&](uint32_t pos, uint64_t x) {
      uint64_t hs = fmix64(x ^ SEED_SEL);
      if (hs >= hmax) return;
      for (size_t t = 0; t < sk.tiers.size(); ++t)
        if (sk.tiers[t].selects(hs)) {
          uint64_t h;
          const uint64_t* blk = sk.tiers[t].block_of(x, h);
          __builtin_prefetch(blk);
          sel.push_back({blk, h});
          pr.push_back({pos, (uint8_t)t, 0});
          break;
        }
    });
    for (size_t i = 0; i < pr.size(); ++i) pr[i].hit = sk.tiers[pr[i].tier].contains_at(sel[i].blk, sel[i].h);
    R.probes = pr.size();
    for (auto& p : pr) R.hits += p.hit;
    if (pr.size() < 3) {
      R.status = "too_few_probes";
      return R;
    }
    const size_t n = pr.size();
    std::vector<double> a;
    double ll, llU, pm_ani = NAN, low_frac = 0;
    double t1 = now();
    T_PROBE += t1 - t0;
    hmm(pr, a, ll, llU, &pm_ani, &low_frac);
    T_HMM += now() - t1;
    R.lr = ll - llU;
    if (R.lr < P.lr_min) {
      R.status = "no_homology";
      return R;
    }
    R.low_frac = low_frac;
    // pooled aligned containment
    double var = 0, neff = 0;
    std::vector<double> ones(n, 1.0);
    double th0 = theta_range(pr, a, 0.3, 0, n);
    double th = theta_range(pr, a, th0, 0, n, &var, &neff);
    R.theta = th;
    R.aligned_probes = neff;
    // aligned fraction: posterior-expected aligned k-mer starts + (k-1) per aligned segment
    const size_t L = seq.size();
    const size_t nstarts = L - k + 1;
    double Sx = a[0] * (pr[0].pos + 1);
    for (size_t i = 1; i < n; ++i) Sx += 0.5 * (a[i - 1] + a[i]) * (double)(pr[i].pos - pr[i - 1].pos);
    Sx += a[n - 1] * (double)(nstarts - 1 - pr[n - 1].pos);
    int segs = 0;
    for (size_t i = 0; i < n; ++i)
      if (a[i] > 0.5 && (i == 0 || a[i - 1] <= 0.5)) ++segs;
    R.af = std::min(1.0, std::max(0.0, (Sx + (double)(k - 1) * segs) / (double)L));
    const double nk_aligned = std::max(1.0, R.af * L - (k - 1));
    if (th <= 0) {
      R.status = "no_homology";
      return R;
    }
    const double ani_pooled = ani_from_theta(th, var, nk_aligned);
    R.ani = ani_pooled;
    if (P.het_mode == 0 || P.window <= 0) return R;
    // heterogeneity test on windows of confidently aligned probes
    std::vector<double> conf(n);
    for (size_t q = 0; q < n; ++q) conf[q] = a[q] > 0.9 ? 1.0 : 0.0;
    std::vector<double> wc, wv;
    size_t i = 0;
    while (i < n) {
      size_t j = i;
      uint32_t start = pr[i].pos;
      while (j < n && pr[j].pos < start + (uint32_t)P.window) ++j;
      if (j < n && pr[n - 1].pos < pr[j].pos + (uint32_t)P.window / 2) j = n;
      double v = 0, ne = 0;
      double tw = theta_range(pr, conf, th, i, j, &v, &ne);
      double frac = ne / std::max<double>(1, (double)(j - i));
      double len = frac * (double)(pr[j - 1].pos - start + 1);
      if (len > 100 && ne >= 10) {
        double vproc = std::max(th, 1e-6) * Bcoef(1 - ani_pooled, k) / len;
        wc.push_back(tw);
        wv.push_back(v + vproc);
      }
      i = j;
    }
    R.windows = (int)wc.size();
    bool het = P.het_mode == 2;
    if (P.het_mode == 1 && wc.size() >= 2) {
      double Q = 0;
      for (size_t w = 0; w < wc.size(); ++w) Q += (wc[w] - th) * (wc[w] - th) / std::max(wv[w], 1e-12);
      double df = (double)wc.size() - 1;
      het = Q > df + P.het_z * std::sqrt(2 * df);
    }
    if (het && pm_ani == pm_ani) {
      R.het = true;
      R.ani = pm_ani;
    }
    if (R.low_frac > P.floor_frac) R.status = "below_floor";
    return R;
  }
};

// sum of squared base frequencies (composition factor of spurious matches per extra base)
double composition_sf2(const std::vector<std::string>& seqs) {
  double cnt[4] = {0, 0, 0, 0};
  for (auto& seq : seqs)
    for (char c : seq) {
      int x = CODE[(uint8_t)c];
      if (x >= 0) cnt[x] += 1;
    }
  double tot = cnt[0] + cnt[1] + cnt[2] + cnt[3], sf2 = 0;
  for (double c : cnt) sf2 += (c / tot) * (c / tot);
  return std::max(sf2, 0.25);
}

void usage() {
  fprintf(stderr,
          "aif sketch <ref.fa> -o <out.aif> [-k auto|N] [--kmin 17] [--kmax 21] [--sigma-target 0.003]\n"
          "           [--budget 0.88] [--tiers 0.005:0.6,0.1:0.4]\n"
          "aif query <ref.aif> <contigs.fa> [--rho-u 3e-6] [--rho-h 1e-4] [--ident-c 8] [--lr-min 10]\n"
          "          [--het 0|1|2 (1)] [--window 1000] [--het-z 3] [--grid-step 0.005] [--floor-frac 2]\n"
          "          [--no-bias-correct]\n"
          "output: contig ani af theta probes hits aligned_probes status lr het low_frac\n"
          "aif info <ref.aif>\n");
  exit(1);
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
      double sigma_target = 0.003;
      double budget = 0.88;
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
        else usage();
      }
      if (out.empty() || (k != 0 && (k < 5 || k > 31))) usage();
      auto spec = parse_tiers(tiers);
      Sketch sk;
      const std::vector<std::string> seqs = load_fasta(in);
      if (k != 0) {
        sk = build_sketch(seqs, k, budget, spec);
      } else {
        // adaptive k: smallest k >= kmin whose calibrated spurious rate is <= sigma_target
        double sf2 = composition_sf2(seqs);
        // start at the k predicted from composition (empirical sigma is ~2-5x the
        // composition-aware estimate 2 n (sum f^2)^k; factor 4 used), then verify
        double npos = 0;
        for (auto& sq : seqs) npos += (double)sq.size();
        int kk = kmin;
        while (kk < kmax && 4.0 * 2.0 * npos * std::pow(sf2, kk) > sigma_target) ++kk;
        sk = build_sketch(seqs, kk, budget, spec);
        while (sk.sigma_emp > sigma_target && kk < kmax) {
          int next = kk + 1;
          while (next < kmax && sk.sigma_emp * std::pow(sf2, next - kk) > sigma_target) ++next;
          kk = next;
          sk = Sketch();
          sk = build_sketch(seqs, kk, budget, spec);
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
      if (argc < 4) usage();
      Sketch sk;
      sk.load(argv[2]);
      Params P;
      for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--rho-u") P.rho_u = atof(argv[++i]);
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
      Fasta fa(argv[3]);
      std::string name, seq;
      printf("contig\tani\taf\ttheta\tprobes\thits\taligned_probes\tstatus\tlr\thet\tlow_frac\n");
      double tq0 = now();
      while (fa.next(name, seq)) {
        Result r = est.estimate(seq);
        printf("%s\t%.6f\t%.6f\t%.6f\t%zu\t%zu\t%.1f\t%s\t%.1f\t%d\t%.3f\n", name.c_str(), r.ani, r.af, r.theta,
               r.probes, r.hits, r.aligned_probes, r.status, r.lr, (int)r.het, r.low_frac);
      }
      if (getenv("AIF_PROFILE"))
        fprintf(stderr, "total %.3f probe %.3f hmm %.3f states %zu\n", now() - tq0, T_PROBE, T_HMM, est.gth.size() + 1);
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
