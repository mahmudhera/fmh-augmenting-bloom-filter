#!/usr/bin/env python3
"""Generate benchmark contigs with exactly known ANI and AF from real genomes.

Each contig is cut from a reference genome, optionally reverse-complemented,
mutated (substitutions, optional indels, optional rate heterogeneity), and
optionally padded with foreign sequence from a different genome (AF < 1) or
assembled from two distant loci (chimera).

Ground truth (per contig):
  true_ani = 1 - substitutions / aligned_columns   (indel columns excluded,
             the usual alignment-based ANI convention)
  true_af  = homologous contig bases / contig length (inserted bases inside a
             homologous segment count as homologous, as in alignment spans)

Usage: python3 make_scenarios.py <genome_dir> <out_dir> [contigs_per_genome] [seed]
Writes <out_dir>/<genome>.contigs.fa, <out_dir>/truth.tsv, <out_dir>/conditions.tsv
"""
import os
import sys

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")

# (name, scenario tag, length, target ANI, AF, indel events / substitution rate,
#  heterogeneity profile, chimera)
CONDITIONS = [
    # S1/S5: length sweep at 95% ANI
    ("len01k_ani95", "S1/S5 length", 1000, 0.95, 1.0, 0.0, None, False),
    ("len02k_ani95", "S1/S5 length", 2000, 0.95, 1.0, 0.0, None, False),
    ("len03k_ani95", "S1/S5 length", 3000, 0.95, 1.0, 0.0, None, False),
    ("len05k_ani95", "S1 length", 5000, 0.95, 1.0, 0.0, None, False),
    ("len10k_ani95", "S1 length", 10000, 0.95, 1.0, 0.0, None, False),
    ("len20k_ani95", "S1 length", 20000, 0.95, 1.0, 0.0, None, False),
    # S2/S3: ANI sweep at 5 kb
    ("ani800_5k", "S3 distant", 5000, 0.80, 1.0, 0.0, None, False),
    ("ani850_5k", "S3 distant", 5000, 0.85, 1.0, 0.0, None, False),
    ("ani900_5k", "S3 distant", 5000, 0.90, 1.0, 0.0, None, False),
    ("ani930_5k", "ANI sweep", 5000, 0.93, 1.0, 0.0, None, False),
    ("ani970_5k", "ANI sweep", 5000, 0.97, 1.0, 0.0, None, False),
    ("ani990_5k", "S2 strain", 5000, 0.99, 1.0, 0.0, None, False),
    ("ani995_5k", "S2 strain", 5000, 0.995, 1.0, 0.0, None, False),
    ("ani999_5k", "S2 strain", 5000, 0.999, 1.0, 0.0, None, False),
    # S4: partial overlap
    ("af20_10k_ani95", "S4 partial AF", 10000, 0.95, 0.2, 0.0, None, False),
    ("af40_10k_ani95", "S4 partial AF", 10000, 0.95, 0.4, 0.0, None, False),
    ("af60_10k_ani95", "S4 partial AF", 10000, 0.95, 0.6, 0.0, None, False),
    ("af80_10k_ani95", "S4 partial AF", 10000, 0.95, 0.8, 0.0, None, False),
    ("af30_10k_ani99", "S4 partial AF", 10000, 0.99, 0.3, 0.0, None, False),
    ("af60_10k_ani99", "S4 partial AF", 10000, 0.99, 0.6, 0.0, None, False),
    # S8: indels
    ("indel10_5k_ani95", "S8 indels", 5000, 0.95, 1.0, 0.10, None, False),
    ("indel20_5k_ani95", "S8 indels", 5000, 0.95, 1.0, 0.20, None, False),
    ("indel10_5k_ani99", "S8 indels", 5000, 0.99, 1.0, 0.10, None, False),
    # S9: heterogeneous rates (mean ANI 95%)
    ("hetblock_10k_ani95", "S9 heterogeneity", 10000, 0.95, 1.0, 0.0, "blocks", False),
    ("hetgamma_10k_ani95", "S9 heterogeneity", 10000, 0.95, 1.0, 0.0, "gamma", False),
    # S11: chimera (two distant loci, both homologous)
    ("chimera_5k_ani95", "S11 chimera", 5000, 0.95, 1.0, 0.0, None, True),
    # combined short/distant and long/distant
    ("len02k_ani90", "S3/S5 short+distant", 2000, 0.90, 1.0, 0.0, None, False),
    ("len20k_ani85", "S3 long+distant", 20000, 0.85, 1.0, 0.0, None, False),
]


def read_fasta(path):
    seq = []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                continue
            seq.append(line.strip().upper())
    return b"".join(seq)


def revcomp(s):
    return s.translate(COMP)[::-1]


def rate_profile(n, p, hetero, rng):
    if hetero is None:
        return np.full(n, p)
    if hetero == "blocks":           # 1 kb blocks, alternately identical and at 2p
        block = 1000
        prof = np.zeros(n)
        first = rng.integers(2)
        for b in range(0, n, block):
            if ((b // block) + first) % 2 == 1:
                prof[b:b + block] = 2 * p
        return prof
    if hetero == "gamma":            # 500 bp blocks, rate ~ Gamma(mean p, CV 1)
        block = 500
        prof = np.empty(n)
        for b in range(0, n, block):
            prof[b:b + block] = min(rng.gamma(1.0, p), 0.5)
        return prof
    raise ValueError(hetero)


def mutate(seg, prof, indel_rate, rng):
    """Return (mutated bytes, substitutions, aligned columns)."""
    arr = np.frombuffer(seg, dtype=np.uint8).copy()
    n = len(arr)
    sub = rng.random(n) < prof
    if sub.any():
        idx = np.nonzero(sub)[0]
        cur = arr[idx]
        # choose a different base uniformly
        code = np.searchsorted(BASES, cur)
        code = np.where(code > 3, 0, code)
        new = (code + rng.integers(1, 4, size=len(idx))) % 4
        arr[idx] = BASES[new]
    n_sub_total = int(sub.sum())
    if indel_rate <= 0:
        return arr.tobytes(), n_sub_total, n
    # indel events (positions sorted); geometric lengths, mean 2, max 10
    ev = np.nonzero(rng.random(n) < indel_rate)[0]
    out = []
    aligned = 0
    nsub = 0
    pos = 0
    for e in ev:
        if e < pos:
            continue
        out.append(arr[pos:e])
        aligned += e - pos
        nsub += int(sub[pos:e].sum())
        d = int(min(rng.geometric(0.5), 10))
        if rng.random() < 0.5:        # deletion of d reference bases
            pos = min(e + d, n)
        else:                         # insertion of d random bases before e
            out.append(BASES[rng.integers(0, 4, size=d)])
            pos = e
    out.append(arr[pos:])
    aligned += n - pos
    nsub += int(sub[pos:].sum())
    return np.concatenate(out).tobytes(), nsub, aligned


def sample_segment(genome, length, rng, avoid=None, min_dist=1_000_000):
    G = len(genome)
    for _ in range(1000):
        st = int(rng.integers(0, G - length))
        if avoid is not None:
            d = abs(st - avoid)
            d = min(d, G - d)
            if d < min_dist:
                continue
        seg = genome[st:st + length]
        if b"N" not in seg:
            return st, seg
    raise RuntimeError("could not sample segment")


def main():
    gdir, out, = sys.argv[1], sys.argv[2]
    per = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 20261121
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(seed)
    names = sorted(f[:-3] for f in os.listdir(gdir) if f.endswith(".fa"))
    genomes = {g: read_fasta(os.path.join(gdir, g + ".fa")) for g in names}
    truth = open(os.path.join(out, "truth.tsv"), "w")
    truth.write("contig\tcondition\tscenario\tgenome\tlength\ttrue_ani\ttrue_af\tn_sub\taligned_cols\thomolog_len\n")
    with open(os.path.join(out, "conditions.tsv"), "w") as fh:
        fh.write("condition\tscenario\tlength\ttarget_ani\taf\tindel_frac\thetero\tchimera\n")
        for c in CONDITIONS:
            fh.write("\t".join(map(str, c)) + "\n")
    for gi, g in enumerate(names):
        genome = genomes[g]
        donor = genomes[names[(gi + 1) % len(names)]]
        with open(os.path.join(out, f"{g}.contigs.fa"), "w") as fa:
            for (cname, scen, L, ani, af, indel_frac, hetero, chim) in CONDITIONS:
                p = 1 - ani
                for i in range(per):
                    name = f"{cname}__{g}__{i:03d}"
                    if chim:
                        half = L // 2
                        st1, s1 = sample_segment(genome, half, rng)
                        _, s2 = sample_segment(genome, L - half, rng, avoid=st1)
                        parts = []
                        nsub = ncol = 0
                        for s in (s1, s2):
                            if rng.random() < 0.5:
                                s = revcomp(s)
                            m, ns, nc = mutate(s, rate_profile(len(s), p, hetero, rng), indel_frac * p, rng)
                            parts.append(m)
                            nsub += ns
                            ncol += nc
                        seq = b"".join(parts)
                        hom = len(seq)
                    else:
                        Lh = int(round(af * L))
                        _, seg = sample_segment(genome, Lh, rng)
                        if rng.random() < 0.5:
                            seg = revcomp(seg)
                        m, nsub, ncol = mutate(seg, rate_profile(Lh, p, hetero, rng), indel_frac * p, rng)
                        hom = len(m)
                        Lf = L - Lh
                        if Lf > 0:
                            _, foreign = sample_segment(donor, Lf, rng)
                            mode = rng.integers(3)
                            if mode == 0:
                                seq = foreign + m
                            elif mode == 1:
                                seq = m + foreign
                            else:
                                a = Lf // 2
                                seq = foreign[:a] + m + foreign[a:]
                        else:
                            seq = m
                    true_ani = 1 - nsub / ncol
                    true_af = hom / len(seq)
                    fa.write(f">{name}\n")
                    for j in range(0, len(seq), 80):
                        fa.write(seq[j:j + 80].decode() + "\n")
                    truth.write(f"{name}\t{cname}\t{scen}\t{g}\t{len(seq)}\t{true_ani:.6f}\t{true_af:.6f}\t{nsub}\t{ncol}\t{hom}\n")
    truth.close()


if __name__ == "__main__":
    main()
