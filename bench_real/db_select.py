#!/usr/bin/env python3
"""Choose the genomes of the database-mode benchmark (data/db/meta/db_accessions.tsv).

  near neighbours: 100 random RefSeq complete genomes of each query species (E. coli, S. aureus,
                   P. aeruginosa, B. cereus group), excluding every assembly used as a query
                   (matched on the accession number, GCA_/GCF_ alike)
  background:      1600 random RefSeq reference bacterial genomes (complete), same exclusion,
                   in a random order; the first 0/400/800/1600 form the nested scaling databases
  references:      the five pairwise references (BA52G, MG1655, PAO1, NCTC8325, REFMAG) are added
                   by db_download.sh from data/real/refs
Inputs: NCBI datasets "summary genome taxon ... --report ids_only" listings in data/db/meta.
"""
import json
import os
import random

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = os.path.join(ROOT, "data", "db", "meta")
SPECIES = [("ecoli", "Escherichia_coli"), ("saur", "Staphylococcus_aureus"), ("paer", "Pseudomonas_aeruginosa"),
           ("bcereus", "Bacillus_cereus_group")]
N_SP, N_BG = 100, 1600


def ids(path):
    return [json.loads(l)["accession"] for l in open(path) if l.strip()]


def num(acc):
    return acc.split("_")[1].split(".")[0]


def main():
    excl = set()
    lists = os.path.join(ROOT, "data", "real", "lists")
    for f in os.listdir(lists):
        for l in list(open(os.path.join(lists, f)))[1:]:
            p = l.rstrip("\n").split("\t")
            if p[2] == "ncbi":
                excl.add(num(p[3]))
    # the real tuning set's queries too (data/real_tune)
    tl = os.path.join(ROOT, "data", "real_tune", "lists")
    if os.path.isdir(tl):
        for f in os.listdir(tl):
            for l in list(open(os.path.join(tl, f)))[1:]:
                p = l.rstrip("\n").split("\t")
                if p[2] == "ncbi":
                    excl.add(num(p[3]))
    rng = random.Random(2026)
    chosen, seen = [], set()
    for grp, fn in SPECIES:
        pool = sorted(a for a in set(ids(os.path.join(M, fn + ".ids.jsonl"))) if num(a) not in excl)
        pick = rng.sample(pool, N_SP)
        for a in pick:
            chosen.append((a, grp, -1))
            seen.add(num(a))
    pool = sorted(a for a in set(ids(os.path.join(M, "bg_ids.jsonl"))) if num(a) not in excl and num(a) not in seen)
    for i, a in enumerate(rng.sample(pool, N_BG)):
        chosen.append((a, "background", i))
    with open(os.path.join(M, "db_accessions.tsv"), "w") as fh:
        fh.write("accession\tgroup\tbg_order\n")
        for a, g, i in chosen:
            fh.write(f"{a}\t{g}\t{i}\n")
    print(f"excluded query accessions: {len(excl)}; chosen: {len(chosen)}")


if __name__ == "__main__":
    main()
