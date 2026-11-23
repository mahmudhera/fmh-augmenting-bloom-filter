#!/usr/bin/env python3
"""GTDB contig-classification benchmark: choose the query genomes.

300 GTDB non-representative genomes from 300 distinct species (bacteria and archaea in proportion
to their species counts, at least 15 archaeal), each a draft NCBI assembly (contig or scaffold
level) with CheckM2 completeness >= 95% and contamination <= 2%, from species whose representative
is in the database. Their contigs >= 500 bp are the queries; the species representative (in the
database) is the expected top hit. Seed 2611. GTDB_SET=2 selects a fresh validation set (seed 4611)
sharing no species with the first.
Writes data/gtdb/meta/queries.tsv (accession, ncbi accession, species, genus, representative).
"""
import csv
import gzip
import os
import random

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
N, N_AR_MIN = 300, 15
SET = os.environ.get("GTDB_SET", "")  # "" = first query set; e.g. "2" = fresh validation set
SEED = 2611 if not SET else 2611 + 1000 * int(SET)


def rows(fn):
    with gzip.open(os.path.join(G, fn), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r


def main():
    rng = random.Random(SEED)
    # species of earlier query sets are excluded (a fresh set shares no species with them)
    used_sp = set()
    for prev in [""] + [str(i) for i in range(2, int(SET or 1))]:
        f = os.path.join(G, "meta", f"queries{prev}.tsv")
        if prev != SET and os.path.exists(f):
            used_sp |= {r["species"] for r in csv.DictReader(open(f), delimiter="\t")}
    pools = {}
    for dom, fn in (("bac", "bac120_metadata.tsv.gz"), ("ar", "ar53_metadata.tsv.gz")):
        by_sp = {}
        for r in rows(fn):
            if r["gtdb_representative"] != "f" or r["ncbi_assembly_level"] not in ("Contig", "Scaffold"):
                continue
            try:
                if float(r["checkm2_completeness"]) < 95 or float(r["checkm2_contamination"]) > 2:
                    continue
            except ValueError:
                continue
            if r["ncbi_genome_category"] not in ("none", "", "na") and "derived from metagenome" in r["ncbi_genome_category"]:
                continue
            sp = r["gtdb_taxonomy"].split(";")[-1]
            if sp[3:] in used_sp:
                continue
            by_sp.setdefault(sp, []).append(r)
        pools[dom] = by_sp
    nsp = {d: len(v) for d, v in pools.items()}
    n_ar = max(N_AR_MIN, round(N * nsp["ar"] / (nsp["ar"] + nsp["bac"])))
    out = []
    for dom, n in (("ar", n_ar), ("bac", N - n_ar)):
        species = sorted(pools[dom])
        for sp in rng.sample(species, n):
            r = rng.choice(sorted(pools[dom][sp], key=lambda x: x["accession"]))
            tax = r["gtdb_taxonomy"].split(";")
            out.append((r["accession"], r["accession"][3:], sp[3:], tax[-2][3:], r["gtdb_genome_representative"], dom,
                        r["contig_count"], r["genome_size"]))
    with open(os.path.join(G, "meta", f"queries{SET}.tsv"), "w") as fh:
        fh.write("accession\tncbi_accession\tspecies\tgenus\trepresentative\tdomain\tcontigs\tgenome_size\n")
        for o in out:
            fh.write("\t".join(o) + "\n")
    print(f"set '{SET}' (seed {SEED}, {len(used_sp)} species excluded): species with eligible non-representatives: "
          f"bacteria {nsp['bac']}, archaea {nsp['ar']}; chose {len(out)} ({n_ar} archaeal)")


if __name__ == "__main__":
    main()
