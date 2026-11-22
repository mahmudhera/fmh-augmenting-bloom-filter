#!/usr/bin/env python3
"""Held-out simulated accuracy of the AF-model options (default, aif-RT, local identity prior)
on the simulated main and fresh benchmarks. Writes data/real/results/afopt_sim.json.

Usage: python3 bench_real/afopt_heldout.py
"""
import csv
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
AIF = os.path.join(ROOT, "bin", "aif")
G = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]
OPTS = {"default": [], "aif_rt": ["--rho-u", "2e-5"], "aif_lp": ["--ident-window", "0.08", "--ident-radius", "3"]}


def main():
    res = {}
    for bench in ("bench", "bench2"):
        truth = list(csv.DictReader(open(os.path.join(ROOT, "data", bench, "truth.tsv")), delimiter="\t"))
        res[bench] = {}
        for name, args in OPTS.items():
            pred = {}
            for g in G:
                out = subprocess.run([AIF, "query", os.path.join(ROOT, "data/results/main/aif", g + ".aif"),
                                      os.path.join(ROOT, "data", bench, g + ".contigs.fa"), "-t", "8"] + args,
                                     capture_output=True, text=True, check=True).stdout
                for line in out.splitlines()[1:]:
                    p = line.split("\t")
                    if p[7] == "ok":
                        pred[p[0].split()[0]] = (float(p[1]), float(p[2]))
            conds = {}
            for t in truth:
                if t["contig"] in pred:
                    a, f = pred[t["contig"]]
                    conds.setdefault(t["condition"], []).append((abs(a - float(t["true_ani"])) * 100, abs(f - float(t["true_af"])) * 100))
            allr = [x for v in conds.values() for x in v]
            res[bench][name] = {
                "ani_macro": sum(sum(x[0] for x in v) / len(v) for v in conds.values()) / len(conds),
                "af_macro": sum(sum(x[1] for x in v) / len(v) for v in conds.values()) / len(conds),
                "ani_micro": sum(x[0] for x in allr) / len(allr), "af_micro": sum(x[1] for x in allr) / len(allr)}
            print(bench, name, {k: round(v, 3) for k, v in res[bench][name].items()})
    json.dump(res, open(os.path.join(ROOT, "data", "real", "results", "afopt_sim.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
