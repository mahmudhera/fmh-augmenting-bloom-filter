#!/usr/bin/env python3
"""Closed forms and Monte Carlo checks of the idealised model for the theory write-up (aif_theory.pdf).

Writes results_theory/numbers.json and report_theory/gen/*.tbl. Everything here evaluates formulas of the
theory or simulates the idealised model (independent selection and false positives, simple substitution model);
checks against the aif binary itself are in scripts/theory_sims.py.

Usage: /usr/bin/python3 -I scripts/theory_numbers.py
"""
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results_theory")
GEN = os.path.join(ROOT, "report_theory", "gen")
os.makedirs(RES, exist_ok=True)
os.makedirs(GEN, exist_ok=True)

# ---------------------------------------------------------------- defaults of aif v1.7
K = 20                      # database k (pairwise: adaptive 17-21)
BUDGET = 0.88               # bits per reference base
SPEC = [(0.005, 0.6), (0.10, 0.4)]   # (target FPR, budget share) per tier
MAX_H = 6
BLOCK = 512
RHO_U, RHO_H = 3e-6, 1e-4
GRID = 0.005
IDENT_C = 8.0
LR_MIN = 10.0
SIGMA = 2.3e-4              # calibrated spurious rate of the random 5 Mbp test genome at k = 20 (aif info)
EPS_MEAS = (0.00635, 0.10156)  # measured tier FPRs of that sketch (aif info)


def c_bf(eps, hmax=30):
    """exact Bloom cost (bits/element) minimised over integer h (as aif's c_bf)"""
    return min(-h / math.log(1 - eps ** (1.0 / h)) for h in range(1, hmax + 1))


def c_bf_real(eps):
    return math.log2(1 / eps) / math.log(2)


def layout(budget=BUDGET, spec=SPEC, max_h=MAX_H):
    """per tier: selection rate s, bits per element c, hash count h, elements per block, standard and blocked FPR"""
    fs = sum(f for _, f in spec)
    out = []
    for e, f in spec:
        c = c_bf(e)
        s = budget * (f / fs) / c
        best, hb = 2.0, 1
        for h in range(1, max_h + 1):
            fpr = (1 - math.exp(-h / c)) ** h
            if fpr < best:
                best, hb = fpr, h
        lam = BLOCK / c
        out.append({"target": e, "share": f / fs, "c": c, "s": s, "h": hb, "per_block": lam,
                    "eps_std": best, "eps_blocked": blocked_fpr(lam, hb)})
    return out


def blocked_fpr(lam, h, B=BLOCK):
    """expected FPR of a blocked Bloom filter: Poisson(lam) elements per block, each setting h of B bits uniformly
    (with replacement); a non-member tests h uniform positions of a uniform block"""
    # distribution of the number r of distinct positions among h uniform draws from B
    pr = np.zeros(h + 1)
    pr[1] = 1.0
    for _ in range(h - 1):
        nx = np.zeros(h + 1)
        for r in range(1, h + 1):
            if pr[r]:
                nx[r] += pr[r] * r / B
                if r + 1 <= h:
                    nx[r + 1] += pr[r] * (B - r) / B
        pr = nx
    jmax = int(lam + 12 * math.sqrt(lam) + 20)
    tot = 0.0
    for j in range(jmax + 1):
        pj = math.exp(-lam + j * math.log(lam) - math.lgamma(j + 1))
        cond = 0.0
        for r in range(1, h + 1):
            if pr[r]:
                cond += pr[r] * sum((-1) ** i * math.comb(r, i) * (1 - i / B) ** (j * h) for i in range(r + 1))
        tot += pj * cond
    return tot


def kl(a, b):
    """KL(Bern(a) || Bern(b))"""
    return a * math.log(a / b) + (1 - a) * math.log((1 - a) / (1 - b))


def levels(k, eps_min):
    amin = max(min(0.9, IDENT_C * eps_min) ** (1.0 / k), 0.70)
    out, a = [], 0.9995
    while a >= amin - 1e-12:
        out.append(a)
        a -= GRID
    return out


def Bcoef(p, k):
    if p <= 1e-12:
        return 0.0
    q = 1 - (1 - p) ** k
    return max(0.0, q * (2 * k - 1 + 2 / p) - 2 * k)


def var_mut_exact(n, p, k):
    """Var(N_mut) of the simple mutation model (Blanca et al.), n k-mers"""
    th = (1 - p) ** k
    q = 1 - th
    if p <= 0:
        return 0.0
    return n * th * Bcoef(p, k) + k * (k - 1) * th * th + 2 * th / (p * p) * ((1 + (k - 1) * th) * p - q)


def v_tier(th, s, e):
    """sampling + false-positive variance of one tier's FP-corrected containment, per expected probe count (x n s)"""
    return th * (1 - th) * (1 - s) + (1 - th) * e / (1 - e)


def tier_eps(lay, sigma=SIGMA, measured=True):
    eps = EPS_MEAS if measured else tuple(t["eps_blocked"] for t in lay)
    return [e + sigma - e * sigma for e in eps]


def main():
    R = {}
    lay = layout()
    epsp = tier_eps(lay)
    s = [t["s"] for t in lay]
    stot = sum(s)
    R["layout"] = lay
    R["eps_prime"] = epsp
    R["s_total"] = stot
    lev = levels(K, min(epsp))
    R["levels"] = {"m": len(lev), "a_min": lev[-1], "a_max": lev[0]}

    # ---- Table: tier layout
    rows = []
    for i, t in enumerate(lay):
        rows.append(f"{i + 1}@{t['target']:g}@{100 * t['share']:.0f}%@{t['c']:.2f}@{t['s']:.4f}@{t['h']}@{t['per_block']:.1f}"
                    f"@{t['eps_std']:.4f}@{t['eps_blocked']:.4f}@{EPS_MEAS[i]:.4f}")
    write_tbl("tab_layout.tbl", "c c c c c c c c c c", "c c c c c c c c c c",
              ["tier@target@budget@bits per@selection@hash@elements@FPR@FPR@FPR",
               "@FPR@share@element $c$@rate $s sub t$@functions@per block@standard@blocked@measured"], rows)

    # ---- spurious matches of a random genome: analytic vs calibrated (aif info of data/theory/ref_k*.aif)
    sig = {}
    for k in range(17, 22):
        p = os.path.join(ROOT, "data", "theory", f"ref_k{k}.aif")
        sig[k] = {"analytic": -math.expm1(-2 * 4999981 / 4.0 ** k)}
    R["sigma_random"] = sig
    # calibration SE: neighbours sampled = selected share of ~1e6 draws; tier-wise binomial SE, inverse-variance combined
    ncal = 1e6
    w = 0.0
    for t, e in zip(lay, EPS_MEAS):
        nt = ncal * t["s"]
        v = e * (1 - e) / nt / (1 - e) ** 2
        w += 1 / v
    R["sigma_cal_se"] = math.sqrt(1 / w)

    # ---- hedge efficiency: J of the two-tier default vs the best single tier and exact 64-bit FMH
    hedge = []
    for a in [0.80, 0.83, 0.85, 0.88, 0.90, 0.93, 0.95, 0.97, 0.99, 0.995, 0.999]:
        th = a ** K
        # combined non-process variance per unit n: 1 / sum_t (s_t / v_t)
        vc = 1.0 / sum(st / v_tier(th, st, e) for st, e in zip(s, epsp))
        J_h = vc * BUDGET / (th * (1 - th))
        # best single tier: minimise J(eps) = c(eps) (1 - s + eps/((1-eps) th)) / (1 - s) ... (s small: 1 + phi)
        best = (1e300, None)
        dbest = (0.0, None)
        for le in np.linspace(math.log(1e-5), math.log(0.9), 900):
            e = math.exp(le)
            c = c_bf(e)
            st = BUDGET / c
            if st >= 1:
                continue
            J = c * v_tier(th, st, e) / (th * (1 - th))
            if J < best[0]:
                best = (J, e)
            ep = e + SIGMA - e * SIGMA
            D = st * kl(th + (1 - th) * ep, ep)  # expected log-likelihood ratio per base (detection, boundaries)
            if D > dbest[0]:
                dbest = (D, e)
        D_h = sum(st * kl(th + (1 - th) * e, e) for st, e in zip(s, epsp))
        st0 = BUDGET / 64
        J_fmh = 64 * v_tier(th, st0, 0.0) / (th * (1 - th))
        hedge.append({"ani": a, "theta": th, "J_hedge": J_h, "J_best": best[0], "eps_best": best[1], "J_fmh": J_fmh,
                      "eff_vs_best": best[0] / J_h, "gain_vs_fmh": J_fmh / J_h,
                      "D_hedge": D_h, "D_best": dbest[0], "eps_dbest": dbest[1], "D_eff": D_h / dbest[0]})
    R["hedge"] = hedge
    rows = [f"{100 * h['ani']:.1f}@{h['theta']:.3g}@{h['J_hedge']:.1f}@{h['gain_vs_fmh']:.1f}@{h['eps_best']:.3f}@{h['eff_vs_best']:.2f}"
            f"@{h['eps_dbest']:.4f}@{h['D_eff']:.2f}" for h in hedge]
    write_tbl("tab_hedge.tbl", "c c c c c c c c", "n n n n n n n n",
              ["ANI@$theta = a sup k$@$J$, aif@gain over@containment@aif share@detection@aif share",
               "(%)@@(two tiers)@exact FMH@best FPR@of best@best FPR@of best"], rows)

    # ---- ANI standard deviation of a fully aligned contig: aif tiers vs exact FMH at equal bits vs alignment
    sdt = []
    for L in [1000, 2000, 5000, 20000]:
        for a in [0.85, 0.90, 0.95, 0.99]:
            p = 1 - a
            th = a ** K
            n = L - K + 1
            vproc = var_mut_exact(n, p, K) / n ** 2
            vc = 1.0 / sum(n * st / v_tier(th, st, e) for st, e in zip(s, epsp))
            st0 = BUDGET / 64
            vf = v_tier(th, st0, 0.0) / (n * st0)
            fac = (a / (K * th)) ** 2
            sd_aif = math.sqrt(fac * (vproc + vc))
            sd_fmh = math.sqrt(fac * (vproc + vf))
            sd_floor = math.sqrt(fac * vproc)
            sd_aln = math.sqrt(p * (1 - p) / L)
            sdt.append({"L": L, "ani": a, "aif": sd_aif, "fmh_equal_bits": sd_fmh, "containment_floor": sd_floor,
                        "alignment": sd_aln})
    R["sd_table"] = sdt
    rows = [f"{d['L'] // 1000}@{100 * d['ani']:.0f}@{100 * d['aif']:.2f}@{100 * d['fmh_equal_bits']:.2f}@"
            f"{100 * d['containment_floor']:.2f}@{100 * d['alignment']:.2f}" for d in sdt]
    write_tbl("tab_sd.tbl", "c c c c c c", "n n n n n n",
              ["length@ANI@aif@exact FMH,@all $k$-mers@alignment",
               "(kb)@(%)@(two tiers)@equal bits@(floor)@(oracle)"], rows)

    # ---- conditional (realised-identity) variance: share of process variance explained by the substitution count
    cond = []
    for a in [0.80, 0.85, 0.90, 0.95, 0.97, 0.99, 0.999]:
        p = 1 - a
        th = a ** K
        B = Bcoef(p, K)
        expl = K * K * th * p / (1 - p)
        cond.append({"ani": a, "B": B, "explained": expl, "share": min(1.0, expl / B) if B > 0 else 1.0})
    R["conditional"] = cond

    # ---- detection: expected LLR per base and minimum length
    m = len(lev)
    det = []
    for a in [0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.95, 0.97, 0.99]:
        th = a ** K
        best = None
        for aj in lev:  # the best-fitting level of the grid (the true identity may lie below it)
            tj = aj ** K
            r_, v_ = 0.0, 0.0
            for st, e in zip(s, epsp):
                pi = th + (1 - th) * e
                pj = tj + (1 - tj) * e
                lh, lm = math.log(pj / e), math.log((1 - pj) / (1 - e))
                r_ += st * (pi * lh + (1 - pi) * lm)
                v_ += st * pi * (1 - pi) * (lh - lm) ** 2
            if best is None or r_ > best[0]:
                best = (r_, v_)
        rate, var = best
        if rate <= 0:
            det.append({"ani": a, "theta": th, "rate": rate, "var": var, "L50": float("inf"), "L95": float("inf"),
                        "in_grid": a >= lev[-1] - 1e-9})
            continue
        off = math.log(2 * m)
        L50 = (LR_MIN + off) / rate
        # 95% detection: L r - off - 1.645 sqrt(L v) = LR_MIN
        A, Bq, C = rate, -1.645 * math.sqrt(var), -(LR_MIN + off)
        x = (-Bq + math.sqrt(Bq * Bq - 4 * A * C)) / (2 * A)
        det.append({"ani": a, "theta": th, "rate": rate, "var": var, "L50": L50, "L95": x * x, "in_grid": a >= lev[-1] - 1e-9})
    R["detection"] = det
    R["detection_offset"] = math.log(2 * m)
    fmt = lambda x: "\\(if" if not math.isfinite(x) else f"{x:,.0f}"
    mi = lambda x: f"{1000 * x:.3g}".replace("-", "\\(mi")
    rows = [f"{100 * d['ani']:.0f}@{d['theta']:.3g}@{mi(d['rate'])}@{fmt(d['L50'])}@{fmt(d['L95'])}" for d in det]
    write_tbl("tab_detect.tbl", "c c c c c", "n n n n n",
              ["ANI (%)@$theta$@nats per kb@50% detected@95% detected", "@@@from (bp)@from (bp)"], rows)

    # ---- boundary resolution: approximation and exact single-change-point posterior (MC)
    bnd = []
    rng = np.random.default_rng(7)
    wt = [st / stot for st in s]
    ebar = sum(w * e for w, e in zip(wt, epsp))
    for a in [0.85, 0.90, 0.95, 0.99]:
        th = a ** K
        DA = sum(w * kl(th + (1 - th) * e, e) for w, e in zip(wt, epsp))   # per probe, aligned data
        DU = sum(w * kl(e, th + (1 - th) * e) for w, e in zip(wt, epsp))   # per probe, unaligned data
        mc = boundary_mc(rng, th, epsp, wt, reps=4000)
        bnd.append({"ani": a, "theta": th, "DA": DA, "DU": DU, "scale_bp": (1 / DA + 1 / DU) / stot,
                    "mc_bias_bp": mc[0] / stot, "mc_sd_bp": mc[1] / stot, "mc_post_sd_bp": mc[2] / stot})
    R["boundary"] = bnd
    rows = [f"{100 * b['ani']:.0f}@{b['DA']:.3f}@{b['DU']:.3f}@{b['scale_bp']:.0f}@{b['mc_bias_bp']:+.1f}@{b['mc_sd_bp']:.0f}@{b['mc_post_sd_bp']:.0f}"
            for b in bnd]
    write_tbl("tab_boundary_ideal.tbl", "c c c c c c c", "n n n n n n n",
              ["ANI@$D sub A$@$D sub U$@scale (16)@mean error@SD of error@mean posterior",
               "(%)@(nats/probe)@(nats/probe)@(bp)@(bp)@(bp)@SD (bp)"], rows)

    # ---- absorption of unaligned stretches by the lowest identity levels
    ab = absorption(lev, s, epsp)
    R["absorption"] = ab
    rows = [f"{x['a']:.3f}@{x['theta']:.4f}@{1000 * x['kl_rate']:.3f}@{x['G_int']:,.0f}@{x['G_edge']:,.0f}" for x in ab["per_level"][-6:]]
    write_tbl("tab_absorb.tbl", "c c c c c", "n n n n n",
              ["level $a sub j$@$theta sub j$@nats per kb@internal gap@flank",
               "@@(U data)@break-even (bp)@break-even (bp)"], rows)

    # ---- power-mean bias of the pooled ANI on mosaic contigs
    pm = []
    for (a1, a2) in [(0.99, 0.95), (0.99, 0.90), (0.97, 0.93), (0.995, 0.85)]:
        for w1 in (0.5, 0.8):
            mean = w1 * a1 + (1 - w1) * a2
            pooled = (w1 * a1 ** K + (1 - w1) * a2 ** K) ** (1 / K)
            var = w1 * (a1 - mean) ** 2 + (1 - w1) * (a2 - mean) ** 2
            approx = (K - 1) / 2 * var / mean
            pm.append({"a1": a1, "a2": a2, "w1": w1, "mean": mean, "pooled": pooled, "bias": pooled - mean, "approx": approx})
    R["power_mean"] = pm
    rows = [f"{100 * x['a1']:.1f} / {100 * x['a2']:.1f}@{100 * x['w1']:.0f}%@{100 * x['mean']:.2f}@{100 * x['pooled']:.2f}"
            f"@{100 * x['bias']:+.2f}@{100 * x['approx']:+.2f}" for x in pm]
    write_tbl("tab_powermean.tbl", "c c c c c c", "c n n n n n",
              ["identities@share of@true ANI@pooled@bias@approx.",
               "(%)@first part@(%)@$M sub k$ (%)@(pp)@(20) (pp)"], rows)

    # ---- screening recall at database scale
    R["screen"] = screening()

    # ---- indels: k-mer survival and the two alignment conventions
    ind = []
    for p, r in [(0.05, 0.0), (0.05, 0.005), (0.05, 0.01), (0.02, 0.002), (0.0, 0.01)]:
        th = (1 - p) ** K * (1 - r) ** K
        ind.append({"p": p, "r": r, "aif": th ** (1 / K), "excl": 1 - p, "incl": 1 - (p + r) / (1 + r / 2)})
    R["indels"] = ind

    # ---- Monte Carlo check of Theorem 1 with two tiers and the inverse-variance combination (idealised model)
    R["mc_thm1"] = mc_thm1(np.random.default_rng(11), s, epsp)
    rows = [f"{x['L']:,}@{100 * x['ani']:.0f}@{x['mean']:.4f}@{x['theta']:.4f}@{x['sd_mc']:.4f}@{x['sd_thy']:.4f}"
            f"@{x['ani_sd_mc'] * 100:.2f}@{x['ani_sd_thy'] * 100:.2f}@{100 * x['ani_bias_mc']:+.2f}" for x in R["mc_thm1"]]
    write_tbl("tab_mc.tbl", "c c c c c c c c c", "n n n n n n n n n",
              ["length@ANI@mean of@$theta$@SD of $Chat$@SD of $Chat$@SD of ANI@SD of ANI@ANI bias",
               "(bp)@(%)@$Chat$ (MC)@@(MC)@(theory)@(MC, pp)@(theory, pp)@(MC, pp)"], rows)

    json.dump(R, open(os.path.join(RES, "numbers.json"), "w"), indent=1, default=float)
    print("wrote", os.path.join(RES, "numbers.json"))


def boundary_mc(rng, th, epsp, wt, reps=4000, nA=400, nU=400):
    """exact posterior of one change point (aligned -> unaligned) under a uniform prior, probe units.
    returns (mean error of the posterior mean, SD of that error, mean posterior SD)"""
    T = len(epsp)
    errs, psd = [], []
    piA = [th + (1 - th) * e for e in epsp]
    for _ in range(reps):
        tier = rng.choice(T, size=nA + nU, p=wt)
        pa = np.array(piA)[tier]
        pu = np.array(epsp)[tier]
        x = np.empty(nA + nU, dtype=bool)
        x[:nA] = rng.random(nA) < pa[:nA]
        x[nA:] = rng.random(nU) < pu[nA:]
        la = np.where(x, np.log(pa), np.log(1 - pa))
        lu = np.where(x, np.log(pu), np.log(1 - pu))
        # change point b: probes [0,b) aligned, [b, n) unaligned
        ca = np.concatenate([[0.0], np.cumsum(la)])
        cu = np.concatenate([[0.0], np.cumsum(lu)])
        ll = ca + (cu[-1] - cu)
        ll -= ll.max()
        w = np.exp(ll)
        w /= w.sum()
        b = np.arange(nA + nU + 1)
        mb = (w * b).sum()
        errs.append(mb - nA)
        psd.append(math.sqrt(max(0.0, (w * b * b).sum() - mb * mb)))
    e = np.array(errs)
    return float(e.mean()), float(e.std()), float(np.mean(psd))


def absorption(lev, s, epsp):
    """posterior odds 'aligned at level j' : 'unaligned' for a stretch of U data inside an alignment (internal)
    or at a contig end (flank), with typical (expected) log-likelihood ratios"""
    m = len(lev)
    k_int = RHO_H ** 2 / (m * RHO_U ** 2)
    k_edge = RHO_H / (m * RHO_U)
    per = []
    for a in lev:
        th = a ** K
        rate = sum(st * kl(e, th + (1 - th) * e) for st, e in zip(s, epsp))  # E_U log(PU/Pj) per base
        per.append({"a": a, "theta": th, "kl_rate": rate,
                    "G_int": (math.log(k_int)) / (rate + RHO_H),
                    "G_edge": math.log(k_edge) / (rate + RHO_H) if k_edge > 1 else 0.0})

    def odds(G, kap):
        return kap * sum(math.exp(-G * (x["kl_rate"] + RHO_H)) for x in per)

    curve = []
    for G in [100, 200, 500, 1000, 1500, 2000, 3000, 4000, 6000, 8000]:
        oi, oe = odds(G, k_int), odds(G, k_edge)
        curve.append({"G": G, "odds_int": oi, "p_int": oi / (1 + oi), "odds_edge": oe, "p_edge": oe / (1 + oe)})
    # break-even lengths of the summed odds
    def solve(kap):
        lo, hi = 0.0, 1e6
        if odds(0, kap) <= 1:
            return 0.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if odds(mid, kap) > 1:
                lo = mid
            else:
                hi = mid
        return lo
    return {"m": m, "kappa_int": k_int, "kappa_edge": k_edge, "per_level": per, "curve": curve,
            "G_star_int": solve(k_int), "G_star_edge": solve(k_edge)}


def screening(R_refs=199923, max_keys=128, scr_frac=1 / 64, zmin=2.0, cmin=2, top=16):
    """recall of the screening stage: the true reference must pass c >= cmin and z >= zmin, and rank in the top
    `top` by FP-corrected screening containment among R_refs unrelated references (whose FPRs follow the GTDB
    database, data/theory/gtdb_dbinfo.tsv). Ties are broken at random. Idealised: keys independent."""
    from scipy import stats
    p = os.path.join(ROOT, "data", "theory", "gtdb_dbinfo.tsv")
    f = np.array([float(l.split("\t")[4]) for l in open(p).readlines()[2:]])
    qs = np.quantile(f, (np.arange(40) + 0.5) / 40)
    nper = len(f) / 40 * R_refs / len(f)
    fmed = float(np.median(f))
    out = []
    for L in [500, 1000, 2000, 3000, 5000, 10000]:
        ns = min(max_keys, int(round((L - K + 1) * scr_frac)))
        for a in [0.80, 0.85, 0.90, 0.95, 0.99]:
            th = a ** K
            ptrue = th + (1 - th) * fmed
            pmf = stats.binom.pmf(np.arange(ns + 1), ns, ptrue)
            thr = max(cmin, math.ceil(ns * fmed + zmin * math.sqrt(ns * fmed * (1 - fmed)) - 1e-9))
            pass_p = pmf[thr:].sum()
            rec = 0.0
            recK = 0.0
            nfalse = 0.0
            for c in range(thr, ns + 1):
                x = (c / ns - fmed) / (1 - fmed)
                # nulls ranking above: containment > x (strictly) and passing the tests; ties: half
                lam_gt, lam_eq = 0.0, 0.0
                for fq in qs:
                    tq = max(cmin, math.ceil(ns * fq + zmin * math.sqrt(ns * fq * (1 - fq)) - 1e-9))
                    cx = ns * (x * (1 - fq) + fq)  # count with the same containment
                    gt = max(tq, math.floor(cx + 1e-9) + 1)
                    lam_gt += nper * stats.binom.sf(gt - 1, ns, fq)
                    if abs(cx - round(cx)) < 1e-9 and round(cx) >= tq:
                        lam_eq += nper * stats.binom.pmf(round(cx), ns, fq)
                # rank = Poisson(lam_gt) + Binomial-ish share of ties ~ Poisson(lam_eq / 2)
                prob_top = stats.poisson.cdf(top - 1, lam_gt + 0.5 * lam_eq)
                rec += pmf[c] * prob_top
                recK += pmf[c] * stats.poisson.cdf(1024 - 1, lam_gt + 0.5 * lam_eq)
            for fq in qs:
                tq = max(cmin, math.ceil(ns * fq + zmin * math.sqrt(ns * fq * (1 - fq)) - 1e-9))
                nfalse += nper * stats.binom.sf(tq - 1, ns, fq)
            out.append({"L": L, "keys": ns, "ani": a, "pass": float(pass_p), "recall_top16": float(rec),
                        "recall_top1024": float(recK), "null_candidates": float(nfalse)})
    rows = []
    for L in sorted({x["L"] for x in out}):
        r = [x for x in out if x["L"] == L and x["ani"] >= 0.85]
        rows.append(f"{L:,}@{r[0]['keys']}@{r[0]['null_candidates']:,.0f}@" +
                    "@".join(f"{100 * x['pass']:.0f} / {100 * x['recall_top16']:.0f} / {100 * x['recall_top1024']:.0f}" for x in r))
    write_tbl("tab_screen.tbl", "c c c c c c c", "n n n c c c c",
              ["contig@screening@unrelated refs.@true reference: passes / top 16 / top 1,024 (%)",
               "(bp)@keys@passing@85% ANI@90% ANI@95% ANI@99% ANI"], rows, span_first=(3, 6))
    return {"fpr_quantiles": [float(x) for x in np.quantile(f, [0.01, 0.1, 0.5, 0.9, 0.99])], "rows": out}


def mc_thm1(rng, s, epsp, reps=20000):
    """two tiers, simple mutation model, independent selection; aif's two-step inverse-variance combination"""
    out = []
    T = len(s)
    for L in [2000, 5000]:
        n = L - K + 1
        for a in [0.85, 0.90, 0.95, 0.99]:
            p = 1 - a
            th = a ** K
            est, anis = [], []
            for r0 in range(0, reps, 500):
                rr = min(500, reps - r0)
                mut = rng.random((rr, L)) < p
                cs = np.concatenate([np.zeros((rr, 1), dtype=np.int32), np.cumsum(mut, axis=1, dtype=np.int32)], axis=1)
                surv = (cs[:, K:] - cs[:, :n]) == 0  # k-mer i survives iff no substitution in [i, i+K)
                u = rng.random((rr, n))
                tier = np.full((rr, n), -1, dtype=np.int8)
                lo = 0.0
                for t in range(T):
                    tier[(u >= lo) & (u < lo + s[t])] = t
                    lo += s[t]
                fp = rng.random((rr, n))
                for i in range(rr):
                    N = np.zeros(T)
                    Y = np.zeros(T)
                    for t in range(T):
                        sel = tier[i] == t
                        N[t] = sel.sum()
                        hit = surv[i, sel] | (fp[i, sel] < epsp[t])
                        Y[t] = hit.sum()
                    c0, _ = combine(N, Y, 0.3, s, epsp)
                    c, v = combine(N, Y, c0, s, epsp)
                    est.append(c)
                    anis.append(ani_bc(c, v, n, p_hat=None))
            e = np.array(est)
            an = np.array(anis)
            vproc = var_mut_exact(n, p, K) / n ** 2
            vc = 1.0 / sum(n * st / v_tier(th, st, ep) for st, ep in zip(s, epsp))
            vt = vproc + vc
            out.append({"L": L, "ani": a, "theta": th, "mean": float(e.mean()), "sd_mc": float(e.std()),
                        "sd_thy": math.sqrt(vt), "ani_sd_mc": float(np.nanstd(an)),
                        "ani_sd_thy": math.sqrt(vt) * a / (K * th), "ani_bias_mc": float(np.nanmean(an) - a)})
            print("mc", L, a, out[-1])
    return out


def combine(N, Y, theta_ref, s, epsp):
    th = min(max(theta_ref, 1e-4), 1.0)
    acc = wsum = 0.0
    for t in range(len(N)):
        if N[t] < 1e-9:
            continue
        e = epsp[t]
        c = (Y[t] / N[t] - e) / (1 - e)
        v = max((th * (1 - th) * max(1 - s[t], 0) + (1 - th) * e / (1 - e)) / N[t], 1e-12)
        acc += c / v
        wsum += 1 / v
    if wsum == 0:
        return 0.0, 0.0
    return acc / wsum, 1 / wsum


def ani_bc(th, var_theta, nk, p_hat=None):
    """aif's ANI with the delta-method bias correction"""
    if th <= 0:
        return float("nan")
    th = min(th, 1.0)
    ani = th ** (1 / K)
    p = 1 - ani
    v = var_theta + (th * Bcoef(p, K) / nk if nk > 0 else 0)
    ani *= 1 + (K - 1) / (2 * K * K) * min(v / (th * th), 1.0)
    return min(ani, 1.0)


def write_tbl(name, head_fmt, body_fmt, head_rows, rows, span_first=None):
    """groff tbl: header rows centred, body numeric; '@' separates columns"""
    ncol = len(head_fmt.split())
    lines = [".TS", "center box tab(@);"]
    if span_first:
        a, b = span_first
        lines.append(" ".join("c" if i <= a else ("s" if i <= b else "c") for i in range(ncol)))
        lines += [head_fmt] * (len(head_rows) - 1)
    else:
        lines += [head_fmt] * len(head_rows)
    lines += [body_fmt + "."]
    out = lines + head_rows + ["_"] + rows + [".TE"]
    open(os.path.join(GEN, name), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
