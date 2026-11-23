#!/usr/bin/env python3
"""GTDB classification benchmark: query contigs (>= 500 bp) of the downloaded query genomes.
Writes data/gtdb/queries/contigs.fa (contig names <gid>__<n>) and data/gtdb/meta/query_contigs.tsv
(contig, gid, accession, species, genus, representative, domain, length). Run with python3 -I.
"""
import csv
import glob
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
SET = os.environ.get("GTDB_SET", "")


def main():
    Q = list(csv.DictReader(open(os.path.join(G, "meta", f"queries{SET}.tsv")), delimiter="\t"))
    qd = os.path.join(G, f"queries{SET}")
    os.makedirs(qd, exist_ok=True)
    fo = open(os.path.join(qd, "contigs.fa"), "w")
    mo = open(os.path.join(G, "meta", f"query_contigs{SET}.tsv"), "w")
    mo.write("contig\tgid\taccession\tspecies\tgenus\trepresentative\tdomain\tlength\n")
    ng = nc = nb = 0
    for i, q in enumerate(Q):
        fs = glob.glob(os.path.join(G, f"qraw{SET}", "ncbi_dataset", "data", q["ncbi_accession"], "*.fna"))
        if not fs:
            continue
        gid = f"gq{i:03d}"
        os.makedirs(os.path.join(qd, "genomes"), exist_ok=True)
        go = open(os.path.join(qd, "genomes", gid + ".fa"), "w")
        seqs, name, buf = [], None, []
        for line in open(fs[0]):
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(buf))
                name, buf = line, []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs.append("".join(buf))
        k = 0
        for s in seqs:
            if len(s) < 500:
                continue
            c = f"{gid}__{k}"
            k += 1
            fo.write(f">{c}\n{s}\n")
            go.write(f">{c}\n{s}\n")
            mo.write(f"{c}\t{gid}\t{q['accession']}\t{q['species']}\t{q['genus']}\t{q['representative']}\t{q['domain']}\t{len(s)}\n")
            nc += 1
            nb += len(s)
        go.close()
        ng += 1
    print(f"{ng} genomes, {nc} contigs >= 500 bp, {nb / 1e6:.1f} Mbp")


if __name__ == "__main__":
    main()
