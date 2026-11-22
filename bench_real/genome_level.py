#!/usr/bin/env python3
"""Genome-level ANI / AF on the real sets, against MUMmer ANIm (dnadiff-style 1-to-1).

Compared:
  skani, skani_slow       skani dist in genome mode (query genome vs reference)
  fastANI                 fastANI (AF = mapped / total 3 kb fragments)
  aif, aif_rt (contig agg.) aif per-contig output aggregated: AF = sum(af_i L_i) / sum(L_i),
                          ANI = aligned-length-weighted mean of contig ANI
  skani --qi (contig agg.) the same aggregation of skani's per-contig output

AF is scored against both the 1-to-1 ANIm AF (af_mae) and the many-to-one AF of the
per-contig truth (af_mae_q; the quantity aif and skani --qi estimate per contig).

Usage: python3 bench_real/genome_level.py   -> data/real/results/genome_level.{tsv,json}
"""
import csv
import json
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
FASTANI = os.environ.get("FASTANI", "/usr/bin/fastANI")
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]


def genomes(s):
    return [l.rstrip("\n").split("\t") for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]]


def run_skani(s, ref, gids, extra, tmp):
    lst = os.path.join(tmp, "l.txt")
    with open(lst, "w") as fh:
        for g in gids:
            fh.write(os.path.join(R, "sets", s, "genomes", g + ".fa") + "\n")
    out = os.path.join(tmp, "o.tsv")
    subprocess.run([SKANI, "dist", "-t", "8", "--ql", lst, "-r", os.path.join(R, "refs", ref + ".fa"), "-o", out] + extra,
                   check=True, capture_output=True)
    res = {}
    for r in list(open(out))[1:]:
        p = r.rstrip("\n").split("\t")
        res[os.path.basename(p[1])[:-3]] = (float(p[2]) / 100, float(p[4]) / 100)
    return res


def run_fastani(s, ref, gids, tmp):
    lst = os.path.join(tmp, "l.txt")
    with open(lst, "w") as fh:
        for g in gids:
            fh.write(os.path.join(R, "sets", s, "genomes", g + ".fa") + "\n")
    out = os.path.join(tmp, "f.tsv")
    subprocess.run([FASTANI, "--ql", lst, "-r", os.path.join(R, "refs", ref + ".fa"), "-o", out, "-t", "8"],
                   check=True, capture_output=True)
    res = {}
    for r in open(out):
        p = r.split()
        res[os.path.basename(p[0])[:-3]] = (float(p[2]) / 100, int(p[3]) / int(p[4]))
    return res


def aggregate(cfg, s):
    L = {r["contig"]: int(r["length"]) for r in csv.DictReader(open(os.path.join(R, "sets", s, "truth.tsv")), delimiter="\t")}
    acc = {}
    for r in csv.DictReader(open(os.path.join(R, "results", cfg, s + ".pred.tsv")), delimiter="\t"):
        g = r["contig"].split("__")[0]
        a = acc.setdefault(g, [0.0, 0.0, 0.0])
        ln = L[r["contig"]]
        a[0] += ln
        if r["reported"] == "1":
            al = float(r["af"]) * ln
            a[1] += al
            a[2] += al * float(r["ani"])
    return {g: (a[2] / a[1] if a[1] else float("nan"), a[1] / a[0]) for g, a in acc.items()}


def main():
    rows = []
    for s in SETS:
        gs = genomes(s)
        truth = {r["gid"]: (float(r["anim_ani"]), float(r["anim_af"]))
                 for r in csv.DictReader(open(os.path.join(R, "sets", s, "genome_truth.tsv")), delimiter="\t")}
        qaf = {}  # many-to-one AF (nucmer --maxmatch, delta-filter -q), from the per-contig truth
        for r in csv.DictReader(open(os.path.join(R, "sets", s, "truth.tsv")), delimiter="\t"):
            a = qaf.setdefault(r["contig"].split("__")[0], [0.0, 0.0])
            a[0] += int(r["length"])
            a[1] += int(r["length"]) * float(r["true_af"])
        by_ref = {}
        for g in gs:
            by_ref.setdefault(g[4], []).append(g[0])
        tools = {"skani": {}, "skani_slow": {}, "fastANI": {}}
        with tempfile.TemporaryDirectory() as tmp:
            for ref, gids in by_ref.items():
                tools["skani"].update(run_skani(s, ref, gids, [], tmp))
                tools["skani_slow"].update(run_skani(s, ref, gids, ["--slow"], tmp))
                tools["fastANI"].update(run_fastani(s, ref, gids, tmp))
        tools["aif_agg"] = aggregate("aif", s)
        tools["aif_rt_agg"] = aggregate("aif_rt", s)
        tools["aif_fast_agg"] = aggregate("aif_fast", s)
        tools["skani_qi_agg"] = aggregate("skani", s)
        tools["skani_slow_qi_agg"] = aggregate("skani_slow", s)
        for g in gs:
            gid = g[0]
            row = {"set": s, "gid": gid, "stratum": g[3], "n50": int(g[7]), "contigs": int(g[5]),
                   "anim_ani": truth[gid][0], "anim_af": truth[gid][1], "q_af": qaf[gid][1] / qaf[gid][0]}
            for t, d in tools.items():
                a, f = d.get(gid, (float("nan"), float("nan")))
                row[t + "_ani"], row[t + "_af"] = a, f
            rows.append(row)
        print(s, "done", file=sys.stderr)
    keys = list(rows[0].keys())
    with open(os.path.join(R, "results", "genome_level.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.5f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    tl = ["skani", "skani_slow", "fastANI", "skani_qi_agg", "skani_slow_qi_agg", "aif_agg", "aif_rt_agg", "aif_fast_agg"]
    summ = {}
    for s in SETS + ["all"]:
        rr = [r for r in rows if s == "all" or r["set"] == s]
        summ[s] = {}
        for t in tl:
            ea = [abs(r[t + "_ani"] - r["anim_ani"]) * 100 for r in rr if r[t + "_ani"] == r[t + "_ani"]]
            ef = [abs(r[t + "_af"] - r["anim_af"]) * 100 for r in rr if r[t + "_af"] == r[t + "_af"]]
            eq = [abs(r[t + "_af"] - r["q_af"]) * 100 for r in rr if r[t + "_af"] == r[t + "_af"]]
            summ[s][t] = {"n": len(rr), "reported": len(ea), "ani_mae": sum(ea) / len(ea) if ea else float("nan"),
                          "af_mae": sum(ef) / len(ef) if ef else float("nan"),
                          "af_mae_q": sum(eq) / len(eq) if eq else float("nan"),
                          "ani_max": max(ea) if ea else float("nan")}
    json.dump(summ, open(os.path.join(R, "results", "genome_level.json"), "w"), indent=1)
    print(f"{'set':10s}" + "".join(f"{t:>22s}" for t in tl))
    for s, d in summ.items():
        print(f"{s:10s}" + "".join(f"  {d[t]['ani_mae']:6.3f} {d[t]['af_mae']:5.2f} {d[t]['af_mae_q']:5.2f}" for t in tl))


if __name__ == "__main__":
    main()
