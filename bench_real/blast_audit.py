#!/usr/bin/env python3
"""blastn audit of the MUMmer truth on contigs where aif and the truth disagree.

Samples 150 'disputed' contigs (aif AF exceeds the MUMmer AF by > 0.2) and 150 controls
from the five within-species sets (seed 5), aligns them to their reference with
blastn -task blastn (word size 11, e-value 1e-5), and compares the union of blastn HSPs
on the contig (a generous upper bound on alignable sequence) with nucmer and aif.

Usage: python3 bench_real/blast_audit.py [aif_config]  -> data/real/results/blast_audit.{tsv,json}
"""
import csv
import json
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "real")
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]


def union(iv):
    tot, cs, ce = 0, None, None
    for a, b in sorted(iv):
        if cs is None or a > ce + 1:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = a, b
        else:
            ce = max(ce, b)
    return tot + (ce - cs + 1 if cs is not None else 0)


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "aif"
    rng = random.Random(5)
    T, P, ref_of, set_of = {}, {}, {}, {}
    for s in SETS:
        g = {l.split("\t")[0]: l.split("\t")[4] for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]}
        for r in csv.DictReader(open(os.path.join(R, "sets", s, "truth.tsv")), delimiter="\t"):
            T[r["contig"]] = r
            ref_of[r["contig"]] = g[r["contig"].split("__")[0]]
            set_of[r["contig"]] = s
        for r in csv.DictReader(open(os.path.join(R, "results", cfg, s + ".pred.tsv")), delimiter="\t"):
            P[r["contig"]] = r
    names = sorted(T)
    dis = [c for c in names if float(P[c]["af"]) - float(T[c]["true_af"]) > 0.2]
    ctl = [c for c in names if float(P[c]["af"]) - float(T[c]["true_af"]) <= 0.2]
    sel = [(c, "disputed") for c in rng.sample(dis, 150)] + [(c, "control") for c in rng.sample(ctl, 150)]
    want = {c for c, _ in sel}
    seqs = {}
    for s in SETS:
        for gid in sorted({c.split("__")[0] for c in want if set_of[c] == s}):
            name = None
            for line in open(os.path.join(R, "sets", s, "genomes", gid + ".fa")):
                if line.startswith(">"):
                    name = line[1:].strip()
                elif name in want:
                    seqs[name] = seqs.get(name, "") + line.strip()
    hsps = {}
    with tempfile.TemporaryDirectory() as td:
        for ref in sorted({ref_of[c] for c in want}):
            db = os.path.join(td, ref)
            subprocess.run(["makeblastdb", "-in", os.path.join(R, "refs", ref + ".fa"), "-dbtype", "nucl", "-out", db],
                           check=True, capture_output=True)
            q = os.path.join(td, ref + ".q.fa")
            with open(q, "w") as fh:
                for c in sorted(want):
                    if ref_of[c] == ref:
                        fh.write(f">{c}\n{seqs[c]}\n")
            out = subprocess.run(["blastn", "-task", "blastn", "-query", q, "-db", db, "-evalue", "1e-5",
                                  "-outfmt", "6 qseqid pident length qstart qend", "-num_threads", "16",
                                  "-max_target_seqs", "50", "-max_hsps", "100"],
                                 check=True, capture_output=True, text=True).stdout
            for line in out.splitlines():
                p = line.split("\t")
                a, b = int(p[3]), int(p[4])
                hsps.setdefault(p[0], []).append((min(a, b), max(a, b)))
    rows = []
    for c, grp in sel:
        L = int(T[c]["length"])
        rows.append({"contig": c, "group": grp, "set": set_of[c], "length": L, "nucmer_af": float(T[c]["true_af"]),
                     "blastn_af": union(hsps.get(c, [])) / L, "aif_af": float(P[c]["af"]),
                     "aif_ani": float(P[c]["ani"])})
    with open(os.path.join(R, "results", "blast_audit.tsv"), "w") as fh:
        keys = list(rows[0])
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.5f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    summ = {"n_disputed_total": len(dis), "n_contigs": len(names)}
    for grp in ("disputed", "control"):
        rr = [r for r in rows if r["group"] == grp]
        m = lambda k: sum(r[k] for r in rr) / len(rr)  # noqa: E731
        summ[grp] = {"n": len(rr), "nucmer_af": m("nucmer_af"), "blastn_af": m("blastn_af"), "aif_af": m("aif_af"),
                     "aif_vs_nucmer_mae": sum(abs(r["aif_af"] - r["nucmer_af"]) for r in rr) / len(rr) * 100,
                     "aif_vs_blastn_mae": sum(abs(r["aif_af"] - r["blastn_af"]) for r in rr) / len(rr) * 100,
                     "blastn_minus_nucmer": (m("blastn_af") - m("nucmer_af")) * 100,
                     "within_10pp_of_blastn": sum(1 for r in rr if abs(r["aif_af"] - r["blastn_af"]) < 0.1)}
    json.dump(summ, open(os.path.join(R, "results", "blast_audit.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
