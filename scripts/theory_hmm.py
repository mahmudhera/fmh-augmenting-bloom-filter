#!/usr/bin/env python3
"""Exact per-probe forward-backward of aif's segment HMM (unbinned, numpy) on simulated probe sequences, for the
boundary-resolution analysis of aif_theory.pdf: aligned stretch at identity a, then unrelated sequence.

Compares, on the same simulated probes:
  two-state  : unaligned + the single true level (the change-point model of Section 5.4)
  all levels : the full model (unaligned + the 28 levels of the default grid at k = 20)
and AF read as all aligned mass (v1.6) or as mass within 0.05 of the true identity (v1.7's rule without windows).

Usage: /usr/bin/python3 -I scripts/theory_hmm.py   -> results_theory/hmm_boundary.json
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_numbers as tn  # noqa: E402

ROOT = tn.ROOT
K = tn.K


def simulate(rng, a, LA=5000, LU=5000, s=None):
    """probe positions, tiers and hits of a contig: LA bases aligned at identity a (simple substitution model,
    overlapping k-mers share substitutions), then LU unrelated bases"""
    L = LA + LU
    n = L - K + 1
    mut = rng.random(LA) < 1 - a
    cs = np.concatenate([[0], np.cumsum(mut)])
    u = rng.random(n)
    tier = np.full(n, -1)
    lo = 0.0
    for t, st in enumerate(s):
        tier[(u >= lo) & (u < lo + st)] = t
        lo += st
    pos = np.flatnonzero(tier >= 0)
    tr = tier[pos]
    eps = np.array(tn.EPS_MEAS) + tn.SIGMA
    hit = rng.random(len(pos)) < eps[tr]
    al = pos + K <= LA  # k-mer entirely inside the aligned stretch
    surv = np.zeros(len(pos), dtype=bool)
    pa = pos[al]
    surv[al] = (cs[pa + K] - cs[pa]) == 0
    hit |= surv
    return pos, tr, hit


def fb(pos, tr, hit, levels, epsp, rho_u=tn.RHO_U, rho_h=tn.RHO_H):
    """scaled forward-backward; returns P(aligned) per probe and P(level j) per probe"""
    m = len(levels)
    th = np.array([a ** K for a in levels])
    epsp = np.asarray(epsp)
    n = len(pos)
    e_u = np.where(hit, epsp[tr], 1 - epsp[tr])
    ph = th[None, :] + (1 - th[None, :]) * epsp[tr][:, None]
    e_a = np.where(hit[:, None], ph, 1 - ph)
    d = np.diff(pos).astype(float)
    stay = np.exp(-rho_u * d)
    u = 1 - stay
    h = -np.expm1(-rho_h * d)
    F = np.empty((n, m + 1))
    f = np.concatenate([[0.5 * e_u[0]], 0.5 / m * e_a[0]])
    f /= f.sum()
    F[0] = f
    for i in range(1, n):
        fu, fa = f[0], f[1:]
        FA = fa.sum()
        nu = (fu * stay[i - 1] + FA * u[i - 1]) * e_u[i]
        na = (fu * u[i - 1] / m + FA * stay[i - 1] * h[i - 1] / m + fa * stay[i - 1] * (1 - h[i - 1])) * e_a[i]
        f = np.concatenate([[nu], na])
        f /= f.sum()
        F[i] = f
    B = np.empty((n, m + 1))
    b = np.ones(m + 1)
    B[-1] = b
    for i in range(n - 1, 0, -1):
        bu, ba = b[0], b[1:]
        w = e_a[i] * ba
        W = w.sum()
        nbu = stay[i - 1] * e_u[i] * bu + u[i - 1] / m * W
        nba = u[i - 1] * e_u[i] * bu + stay[i - 1] * h[i - 1] / m * W + stay[i - 1] * (1 - h[i - 1]) * w
        b = np.concatenate([[nbu], nba])
        b /= b.sum()
        B[i - 1] = b
    G = F * B
    G /= G.sum(axis=1, keepdims=True)
    return 1 - G[:, 0], G[:, 1:]


def af_bases(pos, pa, L):
    """aligned bases from per-probe posteriors: trapezoids between probes, constant to the ends, + (k-1) per run"""
    x = pos.astype(float)
    S = pa[0] * (x[0] + 1) + np.sum(0.5 * (pa[1:] + pa[:-1]) * np.diff(x)) + pa[-1] * (L - K - x[-1])
    runs = int(np.sum((pa > 0.5) & np.concatenate([[True], pa[:-1] <= 0.5])))
    return S + (K - 1) * runs


def main():
    lay = tn.layout()
    s = [t["s"] for t in lay]
    epsp = tn.tier_eps(lay)
    lev = tn.levels(K, min(epsp))
    rng = np.random.default_rng(31)
    out = []
    for a in [0.90, 0.95, 0.99]:
        j0 = int(np.argmin(np.abs(np.array(lev) - a)))
        keep = np.array(lev) >= lev[j0] - 0.05 - 1e-9
        e2, e_all, e_loc = [], [], []
        for r in range(300):
            pos, tr, hit = simulate(rng, a, s=s)
            p2, _ = fb(pos, tr, hit, [lev[j0]], epsp)
            pall, G = fb(pos, tr, hit, lev, epsp)
            ploc = G[:, keep].sum(axis=1)
            for arr, lst in ((p2, e2), (pall, e_all), (ploc, e_loc)):
                lst.append(af_bases(pos, arr, 10000) - 5000)
        res = {"ani": a}
        for nm, lst in (("two_state", e2), ("all_levels", e_all), ("local", e_loc)):
            x = np.array(lst)
            res[nm] = {"mean": float(x.mean()), "sd": float(x.std()), "q": [float(v) for v in np.percentile(x, [1, 10, 50, 90, 99])]}
        out.append(res)
        print(json.dumps(res))
    json.dump(out, open(os.path.join(ROOT, "results_theory", "hmm_boundary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
