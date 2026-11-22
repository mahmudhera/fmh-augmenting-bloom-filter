#!/usr/bin/env python3
"""Build the query contig files for the real-data benchmark.

For every genome in data/real/lists/<set>.tsv: keep contigs >= MIN_LEN bp, rename them
<gid>__<n> (gid = short genome id), and write
  data/real/sets/<set>/genomes/<gid>.fa      one file per query genome (alignment truth)
  data/real/sets/<set>/<ref>.contigs.fa      all query genomes of the set against <ref>
  data/real/sets/<set>/genomes.tsv           gid, genome, accession, stratum, reference, contigs, bases, n50
For mag2328 the reference is the MAG with the highest contig N50 (written to
data/real/refs/REFMAG.fa) and 48 other MAGs are drawn at random (seed 20261121).
When NCBI returned several versions of an accession, the newest is used.
"""
import glob
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
MIN_LEN = int(os.environ.get("MIN_LEN", "500"))
N_MAG = int(os.environ.get("N_MAG", "48"))
SEED = 20261121
PREFIX = {"bcereus": "bc", "ecoli": "ec", "mag2328": "mg", "paer": "pa", "saur": "sa", "negative": "ng"}


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        yield name, "".join(seq)


def n50(lengths):
    s, half = 0, sum(lengths) / 2
    for x in sorted(lengths, reverse=True):
        s += x
        if s >= half:
            return x
    return 0


def ncbi_file(acc):
    if "." in acc:
        return os.path.join(R, "raw", "ncbi", acc + ".fna")
    fs = glob.glob(os.path.join(R, "raw", "ncbi", acc + ".*.fna"))
    return max(fs, key=lambda f: int(f.rsplit(".", 2)[-2]))


def write_fa(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def main():
    rng = random.Random(SEED)
    for set_ in ["bcereus", "ecoli", "mag2328", "paer", "saur", "negative"]:
        if not os.path.exists(os.path.join(R, "lists", set_ + ".tsv")):
            continue
        rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(R, "lists", set_ + ".tsv"))][1:]
        if set_ == "mag2328":
            stats = []
            for r in rows:
                L = [len(s) for _, s in read_fasta(os.path.join(R, "raw", "mag2328", os.path.basename(r[3])))]
                stats.append((n50(L), sum(L), r))
            stats.sort(key=lambda x: (-x[0], -x[1], x[2][1]))
            ref = stats[0][2]
            with open(os.path.join(R, "refs", "REFMAG.fa"), "w") as fh:
                for name, s in read_fasta(os.path.join(R, "raw", "mag2328", os.path.basename(ref[3]))):
                    write_fa(fh, name, s)
            print(f"mag2328 reference: {ref[1]} (N50 {stats[0][0]}, {stats[0][1]} bp)", file=sys.stderr)
            rows = sorted(rng.sample([x[2] for x in stats[1:]], N_MAG), key=lambda r: r[1])
        out = os.path.join(R, "sets", set_)
        os.makedirs(os.path.join(out, "genomes"), exist_ok=True)
        combined = {}
        meta = []
        for i, r in enumerate(rows):
            gid = f"{PREFIX[set_]}{i:03d}"
            src = (os.path.join(R, "raw", "mag2328", os.path.basename(r[3])) if r[2] == "skani-test"
                   else ncbi_file(r[3]))
            ref = r[5]
            if ref not in combined:
                combined[ref] = open(os.path.join(out, f"{ref}.contigs.fa"), "w")
            L = []
            with open(os.path.join(out, "genomes", gid + ".fa"), "w") as gf:
                for j, (name, s) in enumerate(read_fasta(src)):
                    if len(s) < MIN_LEN:
                        continue
                    cname = f"{gid}__{j}"
                    write_fa(gf, cname, s)
                    write_fa(combined[ref], cname, s)
                    L.append(len(s))
            meta.append((gid, r[1], os.path.basename(src), r[4], ref, len(L), sum(L), n50(L)))
        for f in combined.values():
            f.close()
        with open(os.path.join(out, "genomes.tsv"), "w") as fh:
            fh.write("gid\tgenome\tfile\tstratum\treference\tcontigs\tbases\tn50\n")
            for m in meta:
                fh.write("\t".join(map(str, m)) + "\n")
        nc = sum(m[5] for m in meta)
        nb = sum(m[6] for m in meta)
        print(f"{set_}: {len(meta)} genomes, {nc} contigs >= {MIN_LEN} bp, {nb/1e6:.1f} Mbp", file=sys.stderr)


if __name__ == "__main__":
    main()
