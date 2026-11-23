#!/usr/bin/env python3
"""Pareto baseline: alignment with minimap2 (asm20 preset) per contig on the real sets.

AF = union of aligned contig intervals over primary alignments / contig length; ANI = 1 - gap-
compressed divergence (de tag), weighted by aligned contig length over primary alignments
(matches / block length would count every base of a long indel inside one alignment as a
difference, where nucmer splits the alignment; secondary alignments to paralogs bias it low). Reported when a
contig has any alignment. Index: the reference (built on the fly; size reported from a saved
.mmi). Writes data/real/results/minimap2/<set>.pred.tsv and data/real/results/pareto_minimap2.json.
"""
import json
import os
import subprocess
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "real")
MM2 = os.environ.get("MINIMAP2", "minimap2")
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur", "negative"]


def names(fa):
    return [l[1:].split()[0] for l in open(fa) if l.startswith(">")]


def union(iv):
    iv.sort()
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)


def main():
    d = os.path.join(R, "results", "minimap2")
    os.makedirs(d, exist_ok=True)
    tot_t, idx_bytes, bases = 0.0, 0, 0
    for s in SETS:
        refs = sorted(set(l.split("\t")[4] for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]))
        rows = []
        for ref in refs:
            rfa = os.path.join(R, "refs", ref + ".fa")
            q = os.path.join(R, "sets", s, f"{ref}.contigs.fa")
            paf = os.path.join(d, f"{s}.{ref}.paf")
            if not os.path.exists(paf) or "--rerun" in os.sys.argv:
                tf = tempfile.NamedTemporaryFile(delete=False).name
                with open(paf, "w") as fh:
                    subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf, "taskset", "-c", "2", MM2, "-c", "-x", "asm20",
                                    "-t", "1", rfa, q], stdout=fh, stderr=subprocess.DEVNULL, check=True)
                tot_t += float(open(tf).read().split()[-2])
                os.unlink(tf)
            mmi = os.path.join(d, ref + ".mmi")
            if not os.path.exists(mmi):
                subprocess.run([MM2, "-x", "asm20", "-d", mmi, rfa], stderr=subprocess.DEVNULL, check=True)
            if s != "negative":
                idx_bytes += os.path.getsize(mmi)
                bases += sum(len(l.strip()) for l in open(rfa) if not l.startswith(">"))
            iv, mt, bl, L = {}, {}, {}, {}
            for line in open(paf):
                p = line.rstrip("\n").split("\t")
                if "tp:A:P" not in p[12:]:
                    continue
                c = p[0]
                L[c] = int(p[1])
                iv.setdefault(c, []).append((int(p[2]), int(p[3])))
                de = next((float(x[5:]) for x in p[12:] if x.startswith("de:f:")), None)
                if de is None:
                    continue
                ql = int(p[3]) - int(p[2])
                mt[c] = mt.get(c, 0) + (1 - de) * ql
                bl[c] = bl.get(c, 0) + ql
            for c in names(q):
                if c in iv:
                    rows.append((c, mt[c] / bl[c], min(1.0, union(iv[c]) / L[c]), 1))
                else:
                    rows.append((c, float("nan"), 0.0, 0))
        with open(os.path.join(d, f"{s}.pred.tsv"), "w") as fh:
            fh.write("contig\tani\taf\treported\n")
            for c, a, f, r in rows:
                fh.write(f"{c}\t{a:.6f}\t{f:.6f}\t{r}\n")
    jp = os.path.join(R, "results", "pareto_minimap2.json")
    if tot_t == 0 and os.path.exists(jp):  # alignments reused: keep the measured time
        tot_t = json.load(open(jp))["time_s"]
    json.dump({"time_s": tot_t, "index_bits_per_base": 8.0 * idx_bytes / max(1, bases)}, open(jp, "w"), indent=1)
    print(f"minimap2: {tot_t:.1f} s, index {8.0 * idx_bytes / max(1, bases):.1f} bits/base")


if __name__ == "__main__":
    main()
