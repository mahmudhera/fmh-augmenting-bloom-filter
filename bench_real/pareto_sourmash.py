#!/usr/bin/env python3
"""Pareto baseline: exact FracMinHash containment ANI (sourmash) per contig on the real sets.

For each set and scaled value S: the reference is sketched at k = 21, scaled S; every contig is
sketched alone; ANI = C^(1/k) from the containment C of the contig's sketch in the reference's (sourmash's point estimate) (reported when
they share >= 1 hash). sourmash estimates no AF; prediction files carry AF = 1 when reported.
Sketch size: 64 bits per stored hash (sourmash), i.e. 64 / S bits per reference base.
Writes data/real/results/sm_S<S>/<set>.pred.tsv and timing to data/real/results/pareto_sourmash.json.
Run with a Python that has sourmash (e.g. the binchicken env): python -I bench_real/pareto_sourmash.py
"""
import json
import os
import sys
import time

from sourmash import MinHash

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "real")
SETS = [("bcereus", "BA52G"), ("ecoli", "MG1655"), ("mag2328", "REFMAG"), ("paer", "PAO1"), ("saur", "NCTC8325"),
        ("negative", None)]
SCALED = [int(x) for x in (sys.argv[1:] or ["18", "36", "73", "145", "290", "580"])]
K = 21


def fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    out = {}
    tp = os.path.join(R, "results", "pareto_sourmash.json")
    if os.path.exists(tp):
        out = json.load(open(tp))
    for S in SCALED:
        d = os.path.join(R, "results", f"sm_S{S}")
        os.makedirs(d, exist_ok=True)
        tot = 0.0
        for s, _ in SETS:
            refs = sorted(set(l.split("\t")[4] for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]))
            rows = []
            for ref in refs:
                t0 = time.time()
                mr = MinHash(n=0, ksize=K, scaled=S)
                for _, seq in fasta(os.path.join(R, "refs", ref + ".fa")):
                    mr.add_sequence(seq, force=True)
                for c, seq in fasta(os.path.join(R, "sets", s, f"{ref}.contigs.fa")):
                    mc = MinHash(n=0, ksize=K, scaled=S)
                    mc.add_sequence(seq, force=True)
                    if len(mc) == 0 or mc.count_common(mr) == 0:
                        rows.append((c, float("nan"), 0.0, 0))
                        continue
                    a = mc.contained_by(mr) ** (1.0 / K)  # containment ANI point estimate
                    rows.append((c, a, 1.0, 1))
                tot += time.time() - t0
            with open(os.path.join(d, f"{s}.pred.tsv"), "w") as fh:
                fh.write("contig\tani\taf\treported\n")
                for c, a, f, r in rows:
                    fh.write(f"{c}\t{a:.6f}\t{f:.6f}\t{r}\n")
        out[str(S)] = {"time_s": tot, "bits_per_base": 64.0 / S}
        print(f"S={S}: {tot:.1f} s", flush=True)
        json.dump(out, open(tp, "w"), indent=1)


if __name__ == "__main__":
    main()
