#!/usr/bin/env python3
"""Score real-data predictions against MUMmer truth.

Usage: python3 bench_real/evaluate_real.py [--configs a,b,...] [--json out.json] [--quiet]

Per set, and pooled over all sets:
  ANI  on homologous contigs (true AF >= 0.2): coverage (fraction reported), MAE and bias
       on reported contigs, and paired MAE on contigs every compared config reports.
  AF   on all contigs, counting an unreported contig as AF 0 (a tool that reports nothing
       claims no homology); also AF MAE on reported contigs only.
  FP   false homology calls: fraction of non-homologous contigs (true AF < 0.05) that a tool
       reports with an aligned fraction >= 0.2 (reports of a small aligned fraction on such
       contigs are scored by the AF error instead).
Strata: contig length, true ANI and true AF bins (pooled over the five within-species sets).
"""
import csv
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
SETS = os.environ.get("REAL_SETS", "bcereus,ecoli,mag2328,paer,saur,negative").split(",")
HOM = 0.2
LEN_BINS = [(500, 1000), (1000, 2000), (2000, 5000), (5000, 20000), (20000, 10**12)]
ANI_BINS = [(0, 0.90), (0.90, 0.95), (0.95, 0.98), (0.98, 0.995), (0.995, 1.01)]
AF_BINS = [(0, 0.05), (0.05, 0.5), (0.5, 0.9), (0.9, 1.01)]


def mean(x):
    return sum(x) / len(x) if x else float("nan")


def median(x):
    s = sorted(x)
    n = len(s)
    return float("nan") if n == 0 else (s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]))


def load_truth(s):
    return {r["contig"]: (int(r["length"]), float(r["true_af"]), float(r["true_ani"]))
            for r in csv.DictReader(open(os.path.join(R, "sets", s, "truth.tsv")), delimiter="\t")}


def load_pred(cfg, s):
    return {r["contig"]: (float(r["ani"]), float(r["af"]), int(r["reported"]))
            for r in csv.DictReader(open(os.path.join(R, "results", cfg, s + ".pred.tsv")), delimiter="\t")}


def score(contigs, truth, preds, cfgs):
    """contigs: list of names. Returns {cfg: metrics}."""
    hom = [c for c in contigs if truth[c][1] >= HOM]
    non = [c for c in contigs if truth[c][1] < 0.05]
    both = [c for c in hom if all(preds[k][c][2] for k in cfgs)]
    out = {}
    for k in cfgs:
        P = preds[k]
        rep = [c for c in hom if P[c][2]]
        ea = [abs(P[c][0] - truth[c][2]) * 100 for c in rep]
        ba = [(P[c][0] - truth[c][2]) * 100 for c in rep]
        ef = [abs(P[c][1] - truth[c][1]) * 100 for c in contigs]
        bf = [(P[c][1] - truth[c][1]) * 100 for c in contigs]
        efr = [abs(P[c][1] - truth[c][1]) * 100 for c in contigs if P[c][2]]
        out[k] = {
            "n": len(contigs), "n_hom": len(hom), "coverage": len(rep) / len(hom) if hom else float("nan"),
            "ani_mae": mean(ea), "ani_medae": median(ea), "ani_bias": mean(ba),
            "ani_p90": sorted(ea)[int(0.9 * (len(ea) - 1))] if ea else float("nan"),
            "paired_n": len(both), "paired_ani_mae": mean([abs(P[c][0] - truth[c][2]) * 100 for c in both]),
            "af_mae": mean(ef), "af_bias": mean(bf), "af_mae_reported": mean(efr),
            "paired_af_mae": mean([abs(P[c][1] - truth[c][1]) * 100 for c in both]),
            "fp_rate": (sum(1 for c in non if P[c][2] and P[c][1] >= HOM) / len(non)) if non else float("nan"),
            "any_report_rate": (sum(1 for c in non if P[c][2]) / len(non)) if non else float("nan"),
            "n_non": len(non),
        }
    return out


def main():
    argv = sys.argv[1:]
    cfgs = ["skani", "skani_relaxed", "skani_slow", "aif", "aif_0.44", "aif_rt", "aif_0.44_rt", "aif_fast", "aif_lp"]
    if "--configs" in argv:
        cfgs = argv[argv.index("--configs") + 1].split(",")
    jout = argv[argv.index("--json") + 1] if "--json" in argv else None
    quiet = "--quiet" in argv
    truth, preds = {}, {k: {} for k in cfgs}
    members = {}
    for s in SETS:
        t = load_truth(s)
        truth.update(t)
        members[s] = sorted(t)
        for k in cfgs:
            preds[k].update(load_pred(k, s))
    within = [c for s in SETS if s != "negative" for c in members[s]]
    res = {"configs": cfgs, "sets": {}, "strata": {}}
    for s in SETS:
        res["sets"][s] = score(members[s], truth, preds, cfgs)
    res["sets"]["all_within"] = score(within, truth, preds, cfgs)
    res["sets"]["all"] = score([c for s in SETS for c in members[s]], truth, preds, cfgs)
    # macro over the five within-species sets
    res["macro"] = {k: {m: mean([res["sets"][s][k][m] for s in SETS if s != "negative"])
                        for m in ("ani_mae", "paired_ani_mae", "af_mae", "coverage")} for k in cfgs}
    for name, bins, key in (("length", LEN_BINS, lambda c: truth[c][0]),
                            ("true_ani", ANI_BINS, lambda c: truth[c][2] if truth[c][1] >= HOM else -1),
                            ("true_af", AF_BINS, lambda c: truth[c][1])):
        res["strata"][name] = {}
        for lo, hi in bins:
            cs = [c for c in within if lo <= key(c) < hi]
            if cs:
                res["strata"][name][f"{lo}-{hi}"] = score(cs, truth, preds, cfgs)
    if not quiet:
        def table(title, d):
            print(f"\n== {title}")
            print(f"{'':16s}" + "".join(f"{k:>34s}" for k in cfgs))
            print(f"{'':16s}" + "".join(f"{'cov  aniMAE paired  afMAE   FP':>34s}" for _ in cfgs))
            for row, m in d.items():
                line = f"{row[:16]:16s}"
                for k in cfgs:
                    x = m[k]
                    line += f"   {x['coverage']:5.3f} {x['ani_mae']:6.3f} {x['paired_ani_mae']:6.3f} {x['af_mae']:6.2f} {x['fp_rate']:5.3f}"
                print(line)
        table("per set (ANI on true AF>=0.2; AF on all contigs, unreported = 0)", res["sets"])
        for name in res["strata"]:
            table(f"within-species contigs by {name}", res["strata"][name])
        print("\nmacro over 5 within-species sets:")
        for k in cfgs:
            m = res["macro"][k]
            print(f"  {k:14s} cov={m['coverage']:.3f} ANI MAE={m['ani_mae']:.3f} paired={m['paired_ani_mae']:.3f} AF MAE={m['af_mae']:.2f}")
    if jout:
        with open(jout, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
