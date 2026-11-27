#!/usr/bin/env python3
"""Held-out check of --af-local (AF over aligned mass within D of the local identity; chosen on the UHGG tuning set
T, bench_real/uhgg_af.py) on aif's earlier benchmarks, all with their own MUMmer truth:
  simulated   data/bench, data/bench2 (4 genomes x conditions; macro over conditions)
  real        data/real (the 5 within-species sets of the real-genomes report) and data/real_fresh (fresh drafts)
Contig ANI MAE (truth AF >= 0.2, reported), AF MAE and bias (every contig; unreported = AF 0), macro over sets, and
AF bias by truth-ANI class. Writes results_real/af_local_check.json.

  python3 -I bench_real/af_local_check.py [--tune] [binary] [name=opt,opt,...]   (default: bin/aif_dev, new=--af-local,0.05)
  --tune: the simulated (data/tune) and real (data/real_tune) TUNING sets instead of the held-out ones
"""
import csv
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]
BINS = [(0, 0.85, "<0.85"), (0.85, 0.90, "0.85-0.90"), (0.90, 0.95, "0.90-0.95"), (0.95, 0.99, "0.95-0.99"),
        (0.99, 1.01, ">=0.99")]


def query(aif, sketch, fa, args):
    out = subprocess.run([aif, "query", sketch, fa, "-t", "4"] + args, capture_output=True, text=True, check=True).stdout
    pred = {}
    for line in out.splitlines()[1:]:
        p = line.split("\t")
        ok = p[7] == "ok"
        pred[p[0].split()[0]] = (float(p[1]) if ok else math.nan, float(p[2]) if ok else 0.0, ok)
    return pred


def summarise(rows):
    """rows: (group, true_ani, true_af, ani, af, ok)"""
    out = {}
    groups = sorted({r[0] for r in rows})
    for g in groups:
        rs = [r for r in rows if r[0] == g]
        a = [abs(r[3] - r[1]) * 100 for r in rs if r[2] >= 0.2 and r[5]]
        f = [(r[4] - r[2]) * 100 for r in rs]
        out[g] = {"n": len(rs), "ani_mae": sum(a) / max(1, len(a)), "af_mae": sum(map(abs, f)) / len(f),
                  "af_bias": sum(f) / len(f)}
    macro = {k: sum(v[k] for v in out.values()) / len(out) for k in ("ani_mae", "af_mae", "af_bias")}
    byani = {}
    for lo, hi, lab in BINS:
        f = [(r[4] - r[2]) * 100 for r in rows if lo <= r[1] < hi and r[2] > 0]
        if f:
            byani[lab] = (len(f), round(sum(f) / len(f), 2), round(sum(map(abs, f)) / len(f), 2))
    return {"groups": out, "macro": macro, "af_by_true_ani": byani}


def sim(aif, args, bench, ex):
    truth = list(csv.DictReader(open(os.path.join(ROOT, "data", bench, "truth.tsv")), delimiter="\t"))
    pred = {}
    for p in ex.map(lambda g: query(aif, os.path.join(ROOT, "data/results/main/aif", g + ".aif"),
                                    os.path.join(ROOT, "data", bench, g + ".contigs.fa"), args), G):
        pred.update(p)
    rows = []
    for t in truth:
        a, f, ok = pred.get(t["contig"], (math.nan, 0.0, False))
        rows.append((t["condition"], float(t["true_ani"]), float(t["true_af"]), a, f, ok))
    return summarise(rows)


def real(aif, args, root, ex):
    R = os.path.join(ROOT, "data", root)
    rows = []
    for s in SETS:
        d = os.path.join(R, "sets", s)
        if not os.path.isdir(d):
            continue
        refs = sorted({r["reference"] for r in csv.DictReader(open(os.path.join(d, "genomes.tsv")), delimiter="\t")})
        jobs = []
        for ref in refs:
            sk = os.path.join(R, "results", "aif", f"{s}.{ref}.aif")
            if not os.path.exists(sk):
                sk = os.path.join(ROOT, "data", "real", "results", "aif", f"{s}.{ref}.aif")
            jobs.append((sk, os.path.join(d, f"{ref}.contigs.fa")))
        pred = {}
        for p in ex.map(lambda j: query(aif, j[0], j[1], args), jobs):
            pred.update(p)
        for t in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t"):
            a, f, ok = pred.get(t["contig"], (math.nan, 0.0, False))
            rows.append((s, float(t["true_ani"]), float(t["true_af"]), a, f, ok))
    return summarise(rows)


def main():
    """argv: [--tune] [binary] variant specs name=opt,opt,... (default: new=--af-local,0.05)"""
    a = sys.argv[1:]
    tune = bool(a) and a[0] == "--tune"
    a = a[1:] if tune else a
    aif = a[0] if a and "=" not in a[0] else os.path.join(ROOT, "bin", "aif_dev")
    specs = [x for x in a if "=" in x] or ["new=--af-local,0.05"]
    V = [("default", [])] + [(x.split("=", 1)[0], [y for y in x.split("=", 1)[1].split(",") if y]) for x in specs]
    sims, reals = (("tune",), ("real_tune",)) if tune else (("bench", "bench2"), ("real", "real_fresh"))
    res = {"binary": aif, "variants": dict(V)}
    with ThreadPoolExecutor(12) as ex:
        for name, args in V:
            res[name] = {b: sim(aif, args, b, ex) for b in sims}
            res[name].update({r: real(aif, args, r, ex) for r in reals})
    for b in sims + reals:
        for name, _ in V:
            m = res[name][b]["macro"]
            g = res[name][b]["groups"]
            het = " ".join(f"{c}:{g[c]['af_bias']:+.1f}" for c in ("hetblock_10k_ani95", "hetgamma_10k_ani95") if c in g)
            print(f"{b:10s} {name:10s} ANI MAE {m['ani_mae']:.3f}  AF MAE {m['af_mae']:.3f}  AF bias {m['af_bias']:+.3f}  {het}  "
                  f"by true ANI {res[name][b]['af_by_true_ani']}")
    out = os.path.join(ROOT, "results_real", "af_local_check" + ("_tune" if tune else "") + ".json")
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
