#!/usr/bin/env python3
"""Grid search of aif query parameters on the REAL tuning set (data/real_tune).

Uses the sketches written by run_tools.py (REAL_DIR=data/real_tune ... aif). Each grid
point runs aif query on every tuning set and is scored with evaluate_real.score.
Objective: macro ANI MAE / ANI MAE(default) + macro AF MAE / AF MAE(default).

Usage: python3 bench_real/tune_real.py "<grid spec>" [out.tsv]
  grid spec: "--rho-u=3e-6,1e-5 --ident-c=8,16"   (Cartesian product)
"""
import itertools
import math
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("REAL_DIR", os.path.join(ROOT, "data", "real_tune"))
os.environ.setdefault("REAL_SETS", "bcereus,ecoli,mag2328,paer,saur")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaluate_real as ev  # noqa: E402

R = os.environ["REAL_DIR"]
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SETS = os.environ["REAL_SETS"].split(",")


def jobs_for(args):
    out = []
    for s in SETS:
        refs = sorted(set(l.split("\t")[4] for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]))
        out += [(s, r, args) for r in refs]
    return out


def run(job):
    s, ref, args = job
    sk = os.path.join(R, "results", "aif", f"{s}.{ref}.aif")
    q = os.path.join(R, "sets", s, f"{ref}.contigs.fa")
    res = subprocess.run([AIF, "query", sk, q, "-t", "1"] + list(args), capture_output=True, text=True, check=True).stdout
    pred = {}
    for line in res.splitlines()[1:]:
        p = line.split("\t")
        ok = p[7] == "ok"
        pred[p[0]] = (float(p[1]) if ok else math.nan, float(p[2]) if ok else 0.0, int(ok))
    return s, pred


def score(args, pool, truth, members):
    preds = {}
    for s, p in pool.map(run, jobs_for(args)):
        preds.update(p)
    m = {s: ev.score(members[s], truth, {"x": preds}, ["x"])["x"] for s in SETS}
    ani = sum(m[s]["ani_mae"] for s in SETS) / len(SETS)
    af = sum(m[s]["af_mae"] for s in SETS) / len(SETS)
    cov = sum(m[s]["coverage"] for s in SETS) / len(SETS)
    paired = sum(m[s]["ani_medae"] for s in SETS) / len(SETS)
    return ani, af, cov, paired


def main():
    spec = sys.argv[1] if len(sys.argv) > 1 else ""
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(R, "results", "tune_real.tsv")
    axes = []
    for tok in spec.split():
        k, v = tok.split("=", 1)
        axes.append([(k, x) for x in v.split(",")])
    truth, members = {}, {}
    for s in SETS:
        t = ev.load_truth(s)
        truth.update(t)
        members[s] = sorted(t)
    with Pool(32) as pool:
        base = score([], pool, truth, members)
        rows = []
        for combo in itertools.product(*axes) if axes else [()]:
            args = [x for kv in combo for x in kv]
            ani, af, cov, med = score(args, pool, truth, members)
            obj = ani / base[0] + af / base[1]
            rows.append((obj, " ".join(args), ani, af, cov, med))
            print(f"{obj:.4f}  ani={ani:.4f} af={af:.3f} cov={cov:.4f} medAE={med:.4f}  {' '.join(args) or '(default)'}", flush=True)
    rows.sort()
    with open(out, "a") as fh:
        for r in rows:
            fh.write(f"{r[0]:.4f}\t{r[1] or 'default'}\t{r[2]:.4f}\t{r[3]:.3f}\t{r[4]:.4f}\t{r[5]:.4f}\n")
    print("best:", rows[0])


if __name__ == "__main__":
    main()
