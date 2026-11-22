#!/usr/bin/env python3
"""Real-data TUNING set, disjoint from the evaluation genomes (seed 777).

Same sources and references as select_datasets.py, none of the evaluation genomes:
  bcereus  24 D2 drafts (6 each in 85-90, 90-95, 95-98, >=99.5; all 98-99.5 drafts are used
           for evaluation)            vs B. anthracis 52-G
  ecoli    20 D3 drafts               vs E. coli K-12
  mag2328  20 MAGs of bin 2328        vs the most contiguous MAG among the remaining ones
  paer     10 P. aeruginosa drafts    vs PAO1
  saur     10 S. aureus drafts        vs NCTC 8325
Writes data/real_tune/lists/<set>.tsv; data/real_tune/{raw,meta} and the four fixed
references link to data/real.
"""
import json
import os
import random
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(ROOT, "data", "real")
TUNE = os.path.join(ROOT, "data", "real_tune")
DATASETS = os.environ.get("DATASETS", "/home/mrahmanh/miniconda3/envs/meta-sim/bin/datasets")
SEED = 777

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_datasets import acc_of, bin_of  # noqa: E402


def main_accs():
    used = set()
    for s in ["bcereus", "ecoli", "paer", "saur", "negative"]:
        for l in list(open(os.path.join(MAIN, "lists", s + ".tsv")))[1:]:
            used.add(l.split("\t")[3].split(".")[0])
    return used


def ncbi_drafts(taxon, n, rng, used):
    cmd = [DATASETS, "summary", "genome", "taxon", taxon, "--assembly-level", "contig",
           "--assembly-source", "GenBank", "--limit", "400", "--as-json-lines"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    recs = []
    for line in out.splitlines():
        if line.startswith("{"):
            d = json.loads(line)
            nc = int(d.get("assembly_stats", {}).get("number_of_contigs", 0) or 0)
            if (d.get("organism", {}).get("organism_name", "").startswith(taxon) and 20 <= nc <= 2000
                    and d["accession"].split(".")[0] not in used):
                recs.append((d["accession"], nc))
    return rng.sample(sorted(recs), n)


def main():
    rng = random.Random(SEED)
    os.makedirs(os.path.join(TUNE, "lists"), exist_ok=True)
    for d in ("raw", "meta"):
        if not os.path.exists(os.path.join(TUNE, d)):
            os.symlink(os.path.join(MAIN, d), os.path.join(TUNE, d))
    os.makedirs(os.path.join(TUNE, "refs"), exist_ok=True)  # REFMAG.fa here is the tuning set's own
    for r in ("BA52G", "MG1655", "PAO1", "NCTC8325"):
        if not os.path.exists(os.path.join(TUNE, "refs", r + ".fa")):
            os.symlink(os.path.join(MAIN, "refs", r + ".fa"), os.path.join(TUNE, "refs", r + ".fa"))
    used = main_accs()
    rows = {}
    by = {}
    for line in open(os.path.join(MAIN, "meta", "D2_anim.txt")):
        p = line.rstrip("\n").split("\t")
        a = acc_of(p[1]) if len(p) > 2 else None
        if a and a not in used:
            by.setdefault(bin_of(float(p[2])), []).append((p[1], a, float(p[2])))
    bc = []
    for b in ["85-90", "90-95", "95-98", ">=99.5"]:
        for name, a, ani in sorted(rng.sample(sorted(by[b]), 6)):
            bc.append(("bcereus", name.replace(".LargeContigs.fna", ""), "ncbi", a, f"{b}|{ani:.2f}", "BA52G"))
    rows["bcereus"] = bc
    d3 = sorted(set(l.strip().split("/")[-1] for l in open(os.path.join(MAIN, "meta", "d3_l.txt"))
                    if acc_of(l) and acc_of(l) not in used))
    rows["ecoli"] = [("ecoli", n.replace(".LargeContigs.fna", ""), "ncbi", acc_of(n), "D3", "MG1655")
                     for n in sorted(rng.sample(d3, 20))]
    main_mags = set(l.split("\t")[1] for l in list(open(os.path.join(MAIN, "sets", "mag2328", "genomes.tsv")))[1:])
    ref_mag = open(os.path.join(MAIN, "refs", "REFMAG.fa")).readline()  # header of the evaluation reference
    mags = []
    for l in open(os.path.join(MAIN, "meta", "2328_l.txt")):
        m = l.strip()
        nm = os.path.basename(m).replace(".fa", "")
        if m and nm not in main_mags:
            mags.append(("mag2328", nm, "skani-test", m, "SGB2328", "REFMAG"))
    # drop the evaluation reference MAG (identified by its first contig name)
    first = {}
    for r in mags:
        with open(os.path.join(MAIN, "raw", "mag2328", os.path.basename(r[3]))) as fh:
            first[r[1]] = fh.readline()
    rows["mag2328"] = [r for r in mags if first[r[1]].split()[0] != ref_mag.split()[0]]
    rows["paer"] = [("paer", a, "ncbi", a, f"contigs={nc}", "PAO1") for a, nc in ncbi_drafts("Pseudomonas aeruginosa", 10, rng, used)]
    rows["saur"] = [("saur", a, "ncbi", a, f"contigs={nc}", "NCTC8325") for a, nc in ncbi_drafts("Staphylococcus aureus", 10, rng, used)]
    for s, rr in rows.items():
        with open(os.path.join(TUNE, "lists", s + ".tsv"), "w") as fh:
            fh.write("set\tgenome\tsource\taccession\tstratum\treference\n")
            for r in rr:
                fh.write("\t".join(r) + "\n")
        print(s, len(rr), file=sys.stderr)


if __name__ == "__main__":
    main()
