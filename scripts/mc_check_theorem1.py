#!/usr/bin/env python3
"""Monte Carlo sanity check of the closed-form variance in Theorem 1 of the report.

This is NOT an implementation of any sketching method. It simulates the
idealized model used in the theory (simple substitution model, independent
hash-based selection with probability s, independent Bloom-filter false
positives with probability eps) and compares the empirical mean/variance of the
FP-corrected containment estimator and of the ANI estimator with the closed
forms derived in the report.

Usage: python3 mc_check_theorem1.py  (needs numpy)
"""
import math
import numpy as np


def B_term(p, k):
    q = 1 - (1 - p) ** k
    return q * (2 * k - 1 + 2 / p) - 2 * k


def var_nmut_exact(n, p, k):
    """Exact Var(N_mut) for n k-mers (Blanca, Harris, Koslicki, Medvedev 2022)."""
    q = 1 - (1 - p) ** k
    return (n * (1 - q) * (q * (2 * k + 2 / p - 1) - 2 * k)
            + k * (k - 1) * (1 - q) ** 2
            + 2 * (1 - q) / p ** 2 * ((1 + (k - 1) * (1 - q)) * p - q))


def theory(n, p, k, s, eps):
    th = (1 - p) ** k
    v_proc = var_nmut_exact(n, p, k) / n ** 2
    v_samp = th * (1 - th) * (1 - s) / (n * s)
    v_fp = (1 - th) * eps / ((1 - eps) * n * s)
    vC = v_proc + v_samp + v_fp
    vA = vC * (1 - p) ** 2 / (k ** 2 * th ** 2)
    biasA = -(k - 1) / (2 * k ** 2) * vC / th ** 2 * (1 - p)
    return vC, vA, biasA


def simulate(n, p, k, s, eps, reps, rng, batch=500):
    L = n + k - 1
    C_hat, empty = [], 0
    done = 0
    while done < reps:
        b = min(batch, reps - done)
        mut = rng.random((b, L)) < p
        cs = np.concatenate([np.zeros((b, 1), dtype=np.int32), np.cumsum(mut, axis=1, dtype=np.int32)], axis=1)
        U = (cs[:, k:] - cs[:, :-k]) == 0          # k-mer i unmutated, shape (b, n)
        S = rng.random((b, n)) < s                  # FMH / tier selection
        F = rng.random((b, n)) < eps                # AMQ false positive on probe
        Y = (S & (U | F)).sum(axis=1)
        N = S.sum(axis=1)
        ok = N > 0
        empty += int((~ok).sum())
        C_hat.append(((Y[ok] / N[ok]) - eps) / (1 - eps))
        done += b
    return np.concatenate(C_hat), empty


def main():
    rng = np.random.default_rng(20261121)
    k = 21
    cases = [  # (n, p, s, eps)
        (5000, 0.05, 1.0, 0.0),
        (5000, 0.05, 0.02, 0.0),
        (5000, 0.05, 0.02, 0.25),
        (5000, 0.05, 0.005, 0.01),
        (20000, 0.01, 0.01, 0.10),
        (10000, 0.10, 0.05, 0.001),
        (5000, 0.02, 0.10, 0.30),
    ]
    print("n,p,s,eps,mean_C_sim,theta,var_C_sim,var_C_thy,ratio_C,var_ANI_sim,var_ANI_thy,ratio_ANI,bias_ANI_sim,bias_ANI_thy,frac_empty")
    for (n, p, s, eps) in cases:
        reps = 6000
        C, empty = simulate(n, p, k, s, eps, reps, rng)
        th = (1 - p) ** k
        A = np.clip(C, 1e-12, None) ** (1 / k)
        vC, vA, bA = theory(n, p, k, s, eps)
        print(f"{n},{p},{s},{eps},{C.mean():.4f},{th:.4f},{C.var():.3e},{vC:.3e},{C.var()/vC:.3f},"
              f"{A.var():.3e},{vA:.3e},{A.var()/vA:.3f},{A.mean()-(1-p):.2e},{bA:.2e},{empty/reps:.4f}")


if __name__ == "__main__":
    main()
