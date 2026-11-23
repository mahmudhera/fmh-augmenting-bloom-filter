#!/usr/bin/env python3
"""Fresh real-data VALIDATION set for the frozen v1.6 defaults (seed 4242).

Disjoint from the evaluation set (data/real), the tuning set (data/real_tune) and every genome of the
2,005-genome database (data/db/meta/db_accessions.tsv; matched on the accession number):
  bcereus  24 D2 drafts (2 at 85-90% ANIm, all that remain unused; 8 at 90-95; 7 each at 95-98 and
           >= 99.5%)                                                       vs B. anthracis 52-G
  ecoli    24 D3 drafts                                                    vs E. coli K-12
  paer     16 GenBank contig-level P. aeruginosa drafts                    vs PAO1
  saur     16 GenBank contig-level S. aureus drafts                        vs NCTC 8325
(The bin-2328 MAGs are all used by the evaluation and tuning sets.)
Writes data/real_fresh/lists/<set>.tsv; raw/, meta/ and the references link to data/real.
Then: REAL_DIR=data/real_fresh REAL_SETS=bcereus,ecoli,paer,saur bench_real/download.sh, prepare.py, truth.py.
"""
import json
import os
import random
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(ROOT, "data", "real")
TUNE = os.path.join(ROOT, "data", "real_tune")
FRESH = os.path.join(ROOT, "data", "real_fresh")
DATASETS = os.environ.get("DATASETS", "/home/mrahmanh/miniconda3/envs/meta-sim/bin/datasets")
SEED = 4242

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_datasets import acc_of, bin_of  # noqa: E402


def num(a):
    return a.split("_")[1].split(".")[0] if a and "_" in a else a


def used_accs():
    used = set()
    for base in (MAIN, TUNE):
        for s in ["bcereus", "ecoli", "paer", "saur", "negative"]:
            f = os.path.join(base, "lists", s + ".tsv")
            if os.path.exists(f):
                for l in list(open(f))[1:]:
                    used.add(num(l.split("\t")[3].split(".")[0]))
    for l in list(open(os.path.join(ROOT, "data", "db", "meta", "db_accessions.tsv")))[1:]:
        used.add(num(l.split("\t")[0]))
    return used


def ncbi_drafts(taxon, n, rng, used):
    cmd = [DATASETS, "summary", "genome", "taxon", taxon, "--assembly-level", "contig",
           "--assembly-source", "GenBank", "--limit", "600", "--as-json-lines"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    recs = []
    for line in out.splitlines():
        if line.startswith("{"):
            d = json.loads(line)
            nc = int(d.get("assembly_stats", {}).get("number_of_contigs", 0) or 0)
            if (d.get("organism", {}).get("organism_name", "").startswith(taxon) and 20 <= nc <= 2000
                    and num(d["accession"]) not in used):
                recs.append((d["accession"], nc))
    return rng.sample(sorted(recs), n)


def main():
    rng = random.Random(SEED)
    os.makedirs(os.path.join(FRESH, "lists"), exist_ok=True)
    for d in ("raw", "meta"):
        if not os.path.exists(os.path.join(FRESH, d)):
            os.symlink(os.path.join(MAIN, d), os.path.join(FRESH, d))
    os.makedirs(os.path.join(FRESH, "refs"), exist_ok=True)
    for r in ("BA52G", "MG1655", "PAO1", "NCTC8325"):
        if not os.path.exists(os.path.join(FRESH, "refs", r + ".fa")):
            os.symlink(os.path.join(MAIN, "refs", r + ".fa"), os.path.join(FRESH, "refs", r + ".fa"))
    used = used_accs()
    rows = {}
    by = {}
    for line in open(os.path.join(MAIN, "meta", "D2_anim.txt")):
        p = line.rstrip("\n").split("\t")
        a = acc_of(p[1]) if len(p) > 2 else None
        if a and num(a) not in used:
            by.setdefault(bin_of(float(p[2])), []).append((p[1], a, float(p[2])))
    bc = []
    # only 2 unused drafts remain at 85-90%: take both; 8 at 90-95%, 7 each at 95-98% and >= 99.5%
    for b, n in (("85-90", 2), ("90-95", 8), ("95-98", 7), (">=99.5", 7)):
        for name, a, ani in sorted(rng.sample(sorted(by[b]), min(n, len(by[b])))):
            bc.append(("bcereus", name.replace(".LargeContigs.fna", ""), "ncbi", a, f"{b}|{ani:.2f}", "BA52G"))
    rows["bcereus"] = bc
    d3 = sorted(set(l.strip().split("/")[-1] for l in open(os.path.join(MAIN, "meta", "d3_l.txt"))
                    if acc_of(l) and num(acc_of(l)) not in used))
    rows["ecoli"] = [("ecoli", n.replace(".LargeContigs.fna", ""), "ncbi", acc_of(n), "D3", "MG1655")
                     for n in sorted(rng.sample(d3, 24))]
    rows["paer"] = [("paer", a, "ncbi", a, f"contigs={nc}", "PAO1") for a, nc in ncbi_drafts("Pseudomonas aeruginosa", 16, rng, used)]
    rows["saur"] = [("saur", a, "ncbi", a, f"contigs={nc}", "NCTC8325") for a, nc in ncbi_drafts("Staphylococcus aureus", 16, rng, used)]
    for s, rr in rows.items():
        with open(os.path.join(FRESH, "lists", s + ".tsv"), "w") as fh:
            fh.write("set\tgenome\tsource\taccession\tstratum\treference\n")
            for r in rr:
                fh.write("\t".join(r) + "\n")
        print(s, len(rr), file=sys.stderr)


if __name__ == "__main__":
    main()
