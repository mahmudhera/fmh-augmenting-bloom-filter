#!/usr/bin/env python3
"""Score tool predictions against benchmark truth.

Usage: python3 evaluate.py <truth.tsv> <label>=<predictions.tsv> [<label>=<predictions.tsv> ...]
         [--json out.json] [--quiet]

predictions.tsv columns: contig, ani (fraction), af (fraction of query aligned).
For every condition and tool it reports coverage (fraction of contigs with an
output), mean absolute error (MAE) of ANI and AF in percentage points, and mean
signed error (bias). "paired" metrics use only contigs reported by ALL tools.
"""
import csv
import json
import sys
from collections import OrderedDict, defaultdict


def load_pred(path):
    pred = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            name = row["contig"].split()[0]
            pred[name] = (float(row["ani"]), float(row["af"]))
    return pred


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    jout = None
    if "--json" in sys.argv:
        jout = sys.argv[sys.argv.index("--json") + 1]
        args = [a for a in args if a != jout]
    quiet = "--quiet" in sys.argv
    truth_path = args[0]
    tools = OrderedDict()
    for spec in args[1:]:
        label, path = spec.split("=", 1)
        tools[label] = load_pred(path)
    truth = list(csv.DictReader(open(truth_path), delimiter="\t"))
    conds = OrderedDict()
    for t in truth:
        conds.setdefault(t["condition"], []).append(t)
    res = {"conditions": OrderedDict(), "overall": {}}
    for cname, rows in conds.items():
        entry = {"scenario": rows[0]["scenario"], "n": len(rows), "tools": {}}
        both = [t for t in rows if all(t["contig"] in p for p in tools.values())]
        for label, pred in tools.items():
            ea, ef, ba, bf = [], [], [], []
            for t in rows:
                if t["contig"] in pred:
                    a, f = pred[t["contig"]]
                    ea.append(abs(a - float(t["true_ani"])) * 100)
                    ef.append(abs(f - float(t["true_af"])) * 100)
                    ba.append((a - float(t["true_ani"])) * 100)
                    bf.append((f - float(t["true_af"])) * 100)
            pa = [abs(pred[t["contig"]][0] - float(t["true_ani"])) * 100 for t in both]
            pf = [abs(pred[t["contig"]][1] - float(t["true_af"])) * 100 for t in both]
            entry["tools"][label] = {
                "coverage": len(ea) / len(rows), "ani_mae": mean(ea), "af_mae": mean(ef),
                "ani_bias": mean(ba), "af_bias": mean(bf),
                "paired_n": len(both), "paired_ani_mae": mean(pa), "paired_af_mae": mean(pf)}
        res["conditions"][cname] = entry
    for label in tools:
        cs = [c["tools"][label] for c in res["conditions"].values()]
        def macro(key):
            v = [c[key] for c in cs if c[key] == c[key]]
            return mean(v)
        allrows = [t for rows in conds.values() for t in rows]
        pred = tools[label]
        rep = [t for t in allrows if t["contig"] in pred]
        res["overall"][label] = {
            "coverage": len(rep) / len(allrows),
            "ani_mae_micro": mean([abs(pred[t["contig"]][0] - float(t["true_ani"])) * 100 for t in rep]),
            "af_mae_micro": mean([abs(pred[t["contig"]][1] - float(t["true_af"])) * 100 for t in rep]),
            "ani_mae_macro": macro("ani_mae"), "af_mae_macro": macro("af_mae"),
            "paired_ani_mae_macro": macro("paired_ani_mae"), "paired_af_mae_macro": macro("paired_af_mae"),
        }
    if not quiet:
        labels = list(tools)
        hdr = f"{'condition':22s}" + "".join(f" | {l[:10]:>10s} cov  aniMAE  afMAE" for l in labels)
        print(hdr)
        print("-" * len(hdr))
        for cname, e in res["conditions"].items():
            line = f"{cname:22s}"
            for l in labels:
                m = e["tools"][l]
                line += f" | {'':>10s} {m['coverage']:4.2f} {m['ani_mae']:7.3f} {m['af_mae']:6.2f}"
            print(line)
        print("-" * len(hdr))
        for l in labels:
            o = res["overall"][l]
            print(f"{l:12s} coverage={o['coverage']:.3f}  ANI MAE micro={o['ani_mae_micro']:.3f} macro={o['ani_mae_macro']:.3f} "
                  f"paired-macro={o['paired_ani_mae_macro']:.3f} | AF MAE micro={o['af_mae_micro']:.2f} "
                  f"macro={o['af_mae_macro']:.2f} paired-macro={o['paired_af_mae_macro']:.2f}")
    if jout:
        with open(jout, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
