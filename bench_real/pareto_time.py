#!/usr/bin/env python3
"""Pareto benchmark timing: sketch + query wall time (1 thread, core-pinned) of every aif budget and
skani compression factor on the five within-species real sets; 3 interleaved repetitions, median.
Writes data/real/results/pareto_time.json."""
import json
import os
import statistics
import subprocess
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "real")
AIF = os.path.join(ROOT, "bin", "aif")
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
SETS = [("bcereus", "BA52G"), ("ecoli", "MG1655"), ("mag2328", "REFMAG"), ("paer", "PAO1"), ("saur", "NCTC8325")]
CFG = [("aif", b) for b in ("0.11", "0.22", "0.44", "0.88", "1.76")] + [("skani", c) for c in ("500", "250", "125", "60", "30")]


def run(tool, par, s, ref, td):
    q = os.path.join(R, "sets", s, f"{ref}.contigs.fa")
    rfa = os.path.join(R, "refs", ref + ".fa")
    pin = ["taskset", "-c", "3"]
    t0 = time.time()
    if tool == "aif":
        sk = os.path.join(td, "r.aif")
        subprocess.run(pin + [AIF, "sketch", rfa, "-o", sk, "--budget", par], stderr=subprocess.DEVNULL, check=True)
        subprocess.run(pin + [AIF, "query", sk, q, "-t", "1"], stdout=subprocess.DEVNULL, check=True)
    else:
        subprocess.run(pin + [SKANI, "dist", "-t", "1", "--qi", "-q", q, "-r", rfa, "-c", par, "-o", os.path.join(td, "o.tsv")],
                       stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=True)
    return time.time() - t0


def main():
    T = {f"{t}_{p}": [] for t, p in CFG}
    with tempfile.TemporaryDirectory() as td:
        for rep in range(3):
            for t, p in CFG:
                T[f"{t}_{p}"].append(sum(run(t, p, s, ref, td) for s, ref in SETS))
    out = {k: {"median_s": statistics.median(v), "runs": v} for k, v in T.items()}
    json.dump(out, open(os.path.join(R, "results", "pareto_time.json"), "w"), indent=1)
    for k, v in out.items():
        print(k, round(v["median_s"], 2))


if __name__ == "__main__":
    main()
