#!/usr/bin/env python3
"""Genome-level ANI and AF of aif v1.6 (--af-local 0) and v1.7 on the held-out real sets: per-contig output aggregated as
in bench_real/genome_level.py (AF = sum(af L) / sum(L), ANI = aligned-length-weighted mean of contig ANI), against
MUMmer ANIm (1-to-1) and the many-to-one AF of the per-contig truth. Writes results_real/af_local_genome.json."""
import csv
import json
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import af_local_check as A  # noqa: E402

ROOT = A.ROOT
AIF = os.path.join(ROOT, "bin", "aif")


def run(root, args, ex):
    R = os.path.join(ROOT, "data", root)
    out = {}
    for s in A.SETS:
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
        for p in ex.map(lambda j: A.query(AIF, j[0], j[1], args), jobs):
            pred.update(p)
        acc, qaf = {}, {}
        for t in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t"):
            g, ln = t["contig"].split("__")[0], int(t["length"])
            a = acc.setdefault(g, [0.0, 0.0, 0.0])
            q = qaf.setdefault(g, [0.0, 0.0])
            q[0] += ln
            q[1] += ln * float(t["true_af"])
            a[0] += ln
            ani, af, ok = pred.get(t["contig"], (math.nan, 0.0, False))
            if ok:
                a[1] += af * ln
                a[2] += af * ln * ani
        gt = {r["gid"]: (float(r["anim_ani"]), float(r["anim_af"])) for r in
              csv.DictReader(open(os.path.join(d, "genome_truth.tsv")), delimiter="\t")}
        ea, ef, eq = [], [], []
        for g, a in acc.items():
            if g not in gt or a[1] == 0:
                continue
            ani, af = a[2] / a[1], a[1] / a[0]
            if gt[g][0] == gt[g][0]:
                ea.append(abs(ani - gt[g][0]) * 100)
            ef.append(abs(af - gt[g][1]) * 100)
            eq.append(abs(af - qaf[g][1] / qaf[g][0]) * 100)
        out[s] = {"n": len(ef), "ani_mae": sum(ea) / len(ea), "af_mae_anim": sum(ef) / len(ef), "af_mae_q": sum(eq) / len(eq)}
    out["macro"] = {k: sum(v[k] for s, v in out.items()) / len(out) for k in ("ani_mae", "af_mae_anim", "af_mae_q")}
    return out


def main():
    res = {}
    with ThreadPoolExecutor(12) as ex:
        for name, args in (("v16", ["--af-local", "0"]), ("v17", [])):
            res[name] = {r: run(r, args, ex) for r in ("real", "real_fresh")}
    for r in ("real", "real_fresh"):
        for name in ("v16", "v17"):
            m = res[name][r]["macro"]
            print(f"{r:10s} {name} genome ANI MAE {m['ani_mae']:.3f} pp  AF MAE vs ANIm {m['af_mae_anim']:.2f}  vs many-to-one {m['af_mae_q']:.2f}")
    json.dump(res, open(os.path.join(ROOT, "results_real", "af_local_genome.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
