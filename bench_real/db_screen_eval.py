#!/usr/bin/env python3
"""Screening recall of aif db search variants against the exhaustive search (--no-screen).

For every contig with an exhaustive hit (status ok): found = the screened search's top hit is the
exhaustive best reference, or scores (ani*af) within 0.001 of it; reported = any screened hit.
Also: same-score recall split by contig length and by the exhaustive best hit's AF.
Usage: python3 bench_real/db_screen_eval.py var_32_8_3 var_64_8_3 ...  (dirs under data/db/results)
"""
import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "data", "db", "results")
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]


def top(path):
    H = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = r["contig"].split()[0]
        if r["status"] == "ok" and r["rank"] == "1":
            H[c] = (r["reference"], float(r["ani"]) * float(r["af"]), float(r["af"]), int(r["candidates"]))
        H.setdefault(c, None)
    return H


def main():
    out = {}
    times = {}
    tp = os.path.join(RES, "var_times.txt")
    if os.path.exists(tp):
        for l in open(tp):
            p = l.split()
            if len(p) >= 4 and p[3] == "s":
                times.setdefault(p[0], {})[p[1]] = (float(p[2]), int(p[4]))
    X = {s: top(os.path.join(RES, "aif_exh", s + ".tsv")) for s in SETS}
    for v in sys.argv[1:]:
        n = found = rep = 0
        ncand = []
        bylen = {}
        byaf = {}
        for s in SETS:
            V = top(os.path.join(RES, v, s + ".tsv"))
            L = {r["contig"]: int(r["length"]) for r in csv.DictReader(open(os.path.join(ROOT, "data", "real", "sets", s, "truth.tsv")), delimiter="\t")}
            for c, x in X[s].items():
                if x is None:
                    continue
                y = V.get(c)
                n += 1
                ok = y is not None and (y[0] == x[0] or y[1] >= x[1] - 1e-3)
                found += ok
                rep += y is not None
                if y is not None:
                    ncand.append(y[3])
                lb = "<1k" if L[c] < 1000 else "1-3k" if L[c] < 3000 else "3-10k" if L[c] < 10000 else ">=10k"
                ab = "af<0.2" if x[2] < 0.2 else "af<0.9" if x[2] < 0.9 else "af>=0.9"
                for d, k in ((bylen, lb), (byaf, ab)):
                    e = d.setdefault(k, [0, 0])
                    e[0] += ok
                    e[1] += 1
        t = times.get(v.replace("var_", ""), {})
        out[v] = {"n": n, "recall": found / n, "reported": rep / n, "mean_candidates": sum(ncand) / max(1, len(ncand)),
                  "by_length": {k: a / b for k, (a, b) in bylen.items()}, "by_af": {k: a / b for k, (a, b) in byaf.items()},
                  "time_s_t8": sum(x[0] for x in t.values()), "max_kb": max([x[1] for x in t.values()] or [0])}
        o = out[v]
        print(f"{v:12s} n={n} recall={o['recall']:.4f} reported={o['reported']:.4f} cand={o['mean_candidates']:.1f} "
              f"time(t8, 5 sets)={o['time_s_t8']:.1f}s  len:{ {k: round(x, 4) for k, x in o['by_length'].items()} } "
              f"af:{ {k: round(x, 4) for k, x in o['by_af'].items()} }")
    json.dump(out, open(os.path.join(RES, "screen_variants.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
