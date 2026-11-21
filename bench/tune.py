#!/usr/bin/env python3
"""Parallel parameter search for aif on the TUNING set (never the main benchmark).

Usage: python3 tune.py <tune_dir> <skani_tune_predictions.tsv> <out.tsv> [--quick]
"""
import itertools
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
AIF = os.path.join(ROOT, "bin", "aif")
GENOMES = sorted(f[:-3] for f in os.listdir(os.path.join(ROOT, "data", "genomes")) if f.endswith(".fa"))


def sketch(k, tiers, budget, cache):
    d = os.path.join(cache, f"k{k}_{tiers.replace(':', '-').replace(',', '_')}_b{budget}")
    os.makedirs(d, exist_ok=True)
    for g in GENOMES:
        out = os.path.join(d, g + ".aif")
        if not os.path.exists(out):
            subprocess.run([AIF, "sketch", os.path.join(ROOT, "data", "genomes", g + ".fa"), "-o", out,
                            "-k", str(k), "--tiers", tiers, "--budget", str(budget)],
                           check=True, stderr=subprocess.DEVNULL)
    return d


def run(cfg, tune_dir, skani_pred, cache):
    k, tiers, budget, qargs = cfg
    d = sketch(k, tiers, budget, cache)
    with tempfile.TemporaryDirectory() as tmp:
        pred = os.path.join(tmp, "pred.tsv")
        with open(pred, "w") as out:
            out.write("contig\tani\taf\n")
            for g in GENOMES:
                r = subprocess.run([AIF, "query", os.path.join(d, g + ".aif"), os.path.join(tune_dir, g + ".contigs.fa")]
                                   + qargs, check=True, capture_output=True, text=True)
                for line in r.stdout.splitlines()[1:]:
                    f = line.split("\t")
                    if f[7] == "ok":
                        out.write(f"{f[0]}\t{f[1]}\t{f[2]}\n")
        js = os.path.join(tmp, "e.json")
        subprocess.run([sys.executable, os.path.join(HERE, "evaluate.py"), os.path.join(tune_dir, "truth.tsv"),
                        f"skani={skani_pred}", f"aif={pred}", "--json", js, "--quiet"], check=True)
        res = json.load(open(js))
    o = res["overall"]["aif"]
    s = res["overall"]["skani"]
    worst = max((c["tools"]["aif"]["ani_mae"] - c["tools"]["skani"]["ani_mae"]) for c in res["conditions"].values()
                if c["tools"]["aif"]["ani_mae"] == c["tools"]["aif"]["ani_mae"] and c["tools"]["skani"]["ani_mae"] == c["tools"]["skani"]["ani_mae"])
    score = o["ani_mae_macro"] / s["ani_mae_macro"] + o["af_mae_macro"] / s["af_mae_macro"] + 2 * (1 - o["coverage"])
    return cfg, o, score, worst


def main():
    tune_dir, skani_pred, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    cache = os.path.join(ROOT, "data", "tune_sketch")
    ks = [17, 18, 19]
    tiers = ["0.005:0.5,0.1:0.5", "0.005:0.6,0.1:0.4", "0.004:0.5,0.08:0.5", "0.007:0.5,0.15:0.5",
             "0.01:1", "0.007:1", "0.01:0.6,0.2:0.4"]
    qgrid = []
    for ru, rh, ic, lr in itertools.product([3e-6, 1e-5, 3e-5], [3e-5, 1e-4, 3e-4], [4, 6, 8], [6, 10]):
        qgrid.append(["--rho-u", str(ru), "--rho-h", str(rh), "--ident-c", str(ic), "--lr-min", str(lr)])
    if os.environ.get("ROUND") == "1":
        ks = [17, 18, 19]
        tiers = ["0.003:0.5,0.05:0.5", "0.002:0.4,0.03:0.6", "0.005:0.5,0.1:0.5", "0.003:0.3,0.05:0.7",
                 "0.001:0.5,0.03:0.5", "0.01:1", "0.02:1"]
        qgrid = []
        for ru, rh, ic in itertools.product([1e-5, 1e-4, 1e-3], [1e-4, 1e-3], [2, 4, 8]):
            qgrid.append(["--rho-u", str(ru), "--rho-h", str(rh), "--ident-c", str(ic)])
    if "--quick" in sys.argv:
        ks, tiers = [17], tiers[:2]
    cfgs = [(k, t, 0.85, q) for k in ks for t in tiers for q in qgrid]
    # build sketches first (serially per config dir, parallel across dirs)
    with ThreadPoolExecutor(32) as ex:
        list(ex.map(lambda kt: sketch(kt[0], kt[1], 0.85, cache), set((c[0], c[1]) for c in cfgs)))
    with ThreadPoolExecutor(48) as ex:
        results = list(ex.map(lambda c: run(c, tune_dir, skani_pred, cache), cfgs))
    results.sort(key=lambda r: r[2])
    with open(outp, "w") as fh:
        fh.write("score\tk\ttiers\tqargs\tcoverage\tani_micro\tani_macro\taf_micro\taf_macro\tworst_cond_vs_skani\n")
        for cfg, o, score, worst in results:
            fh.write(f"{score:.4f}\t{cfg[0]}\t{cfg[1]}\t{' '.join(cfg[3])}\t{o['coverage']:.3f}\t{o['ani_mae_micro']:.3f}\t"
                     f"{o['ani_mae_macro']:.3f}\t{o['af_mae_micro']:.3f}\t{o['af_mae_macro']:.3f}\t{worst:.3f}\n")
    for cfg, o, score, worst in results[:15]:
        print(f"{score:.3f} k={cfg[0]} tiers={cfg[1]} {' '.join(cfg[3])} cov={o['coverage']:.3f} "
              f"ani={o['ani_mae_micro']:.3f}/{o['ani_mae_macro']:.3f} af={o['af_mae_micro']:.2f}/{o['af_mae_macro']:.2f} worst={worst:+.2f}")


if __name__ == "__main__":
    main()
