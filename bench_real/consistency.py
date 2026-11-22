#!/usr/bin/env python3
"""Is aif's error lower than skani's consistently, not just on average?

Compares an aif configuration with skani (default) on the real sets in every cell of:
set, set x contig length, set x true ANI, set x true AF, and per genome.
  ANI: paired MAE on homologous contigs (true AF >= 0.2) that BOTH tools report
       (coverage differs, so unpaired MAEs compare different contigs); a 95% bootstrap
       interval of the paired difference is given per cell.
  AF:  MAE over all contigs of the cell (unreported = AF 0), with a paired bootstrap interval.
Cells with fewer than MIN contigs for a metric are reported but not counted.

Usage: python3 bench_real/consistency.py [aif_config] [--json out.json]
"""
import csv
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaluate_real as ev  # noqa: E402

SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur", "negative"]
MIN = 20
LEN = [(500, 1000, "0.5-1k"), (1000, 2000, "1-2k"), (2000, 5000, "2-5k"), (5000, 20000, "5-20k"), (20000, 10**12, ">=20k")]
ANI = [(0, 0.90, "<90"), (0.90, 0.95, "90-95"), (0.95, 0.98, "95-98"), (0.98, 0.995, "98-99.5"), (0.995, 1.01, ">=99.5")]
AFB = [(0, 0.05, "<0.05"), (0.05, 0.5, "0.05-0.5"), (0.5, 0.9, "0.5-0.9"), (0.9, 1.01, ">=0.9")]


def boot(diffs, rng, n=1000):
    if len(diffs) < 2:
        return (float("nan"), float("nan"))
    m = len(diffs)
    means = []
    for _ in range(n):
        s = 0.0
        for _ in range(m):
            s += diffs[rng.randrange(m)]
        means.append(s / m)
    means.sort()
    return (means[int(0.025 * n)], means[int(0.975 * n) - 1])


def cell(names, T, A, S, rng):
    """Paired ANI and AF comparison on a list of contigs; differences are aif - skani (pp)."""
    hom = [c for c in names if T[c][1] >= 0.2 and A[c][2] and S[c][2]]
    da = [(abs(A[c][0] - T[c][2]) - abs(S[c][0] - T[c][2])) * 100 for c in hom]
    df = [(abs(A[c][1] - T[c][1]) - abs(S[c][1] - T[c][1])) * 100 for c in names]
    r = {"n": len(names), "n_paired": len(hom)}
    if hom:
        r["ani_aif"] = sum(abs(A[c][0] - T[c][2]) for c in hom) / len(hom) * 100
        r["ani_skani"] = sum(abs(S[c][0] - T[c][2]) for c in hom) / len(hom) * 100
        r["ani_diff"] = sum(da) / len(da)
        r["ani_ci"] = boot(da, rng) if len(hom) >= MIN else (float("nan"), float("nan"))
    if names:
        r["af_aif"] = sum(abs(A[c][1] - T[c][1]) for c in names) / len(names) * 100
        r["af_skani"] = sum(abs(S[c][1] - T[c][1]) for c in names) / len(names) * 100
        r["af_diff"] = sum(df) / len(df)
        r["af_ci"] = boot(df, rng) if len(names) >= MIN else (float("nan"), float("nan"))
    hom_all = [c for c in names if T[c][1] >= 0.2]
    if hom_all:
        r["cov_aif"] = sum(1 for c in hom_all if A[c][2]) / len(hom_all)
        r["cov_skani"] = sum(1 for c in hom_all if S[c][2]) / len(hom_all)
    return r


def main():
    argv = sys.argv[1:]
    cfg = argv[0] if argv and not argv[0].startswith("--") else "aif"
    jout = argv[argv.index("--json") + 1] if "--json" in argv else None
    rng = random.Random(1)
    out = {"config": cfg, "cells": []}
    for s in SETS:
        T = ev.load_truth(s)
        A = ev.load_pred(cfg, s)
        S = ev.load_pred("skani", s)
        names = sorted(T)
        groups = [("set", "all", names)]
        groups += [("length", lab, [c for c in names if lo <= T[c][0] < hi]) for lo, hi, lab in LEN]
        groups += [("true_ani", lab, [c for c in names if T[c][1] >= 0.2 and lo <= T[c][2] < hi]) for lo, hi, lab in ANI]
        groups += [("true_af", lab, [c for c in names if lo <= T[c][1] < hi]) for lo, hi, lab in AFB]
        gids = sorted(set(c.split("__")[0] for c in names))
        groups += [("genome", g, [c for c in names if c.startswith(g + "__")]) for g in gids]
        for kind, lab, cs in groups:
            if not cs:
                continue
            r = cell(cs, T, A, S, rng)
            r.update(set=s, kind=kind, label=lab)
            out["cells"].append(r)
    # summary per kind
    summ = {}
    for kind in ("set", "length", "true_ani", "true_af", "genome"):
        cs = [c for c in out["cells"] if c["kind"] == kind]
        a = [c for c in cs if c["n_paired"] >= MIN and "ani_diff" in c]
        f = [c for c in cs if c["n"] >= MIN and "af_diff" in c]
        summ[kind] = {
            "ani_cells": len(a), "ani_aif_better": sum(1 for c in a if c["ani_diff"] < 0),
            "ani_aif_sig_worse": sum(1 for c in a if c["ani_ci"][0] > 0),
            "af_cells": len(f), "af_aif_better": sum(1 for c in f if c["af_diff"] < 0),
            "af_aif_sig_worse": sum(1 for c in f if c["af_ci"][0] > 0),
        }
    out["summary"] = summ
    print(f"aif config: {cfg} vs skani (default). Cells with >= {MIN} contigs. 'sig. worse' = 95% CI of aif-skani above 0.")
    print(f"{'stratification':12s} {'ANI cells':>9s} {'aif better':>10s} {'sig. worse':>10s} | {'AF cells':>8s} {'aif better':>10s} {'sig. worse':>10s}")
    for kind, m in summ.items():
        print(f"{kind:12s} {m['ani_cells']:9d} {m['ani_aif_better']:10d} {m['ani_aif_sig_worse']:10d} | {m['af_cells']:8d} "
              f"{m['af_aif_better']:10d} {m['af_aif_sig_worse']:10d}")
    print("\ncells where aif's error is higher (ANI paired / AF), with 95% CI of the difference (pp):")
    for c in out["cells"]:
        bad = []
        if c["n_paired"] >= MIN and c.get("ani_diff", -1) > 0:
            bad.append(f"ANI {c['ani_aif']:.3f} vs {c['ani_skani']:.3f} (n={c['n_paired']}, diff {c['ani_diff']:+.3f} "
                       f"[{c['ani_ci'][0]:+.3f},{c['ani_ci'][1]:+.3f}])")
        if c["n"] >= MIN and c.get("af_diff", -1) > 0:
            bad.append(f"AF {c['af_aif']:.2f} vs {c['af_skani']:.2f} (n={c['n']}, diff {c['af_diff']:+.2f} "
                       f"[{c['af_ci'][0]:+.2f},{c['af_ci'][1]:+.2f}])")
        if bad:
            print(f"  {c['set']:9s} {c['kind']:8s} {c['label']:10s} " + "; ".join(bad))
    if jout:
        with open(jout, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
