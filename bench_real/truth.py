#!/usr/bin/env python3
"""Alignment-based (MUMmer4) truth for the real-data benchmark.

For every query genome, nucmer against the set's reference, then
  per contig   nucmer --maxmatch (anchors may be repeated in the reference, so contigs that are
               multi-copy elements such as IS sequences still align), delta-filter -q (best
               reference placement for every query base):
               true_af  = union of aligned query intervals / contig length
               true_ani = aligned-length-weighted %IDY of those alignments (indels count
                          as differences, as in ANIm / dnadiff)
  per genome   nucmer with default seeds and delta-filter -1 (1-to-1), as dnadiff:
               length-weighted identity and the fraction of query bases aligned
               (genome-level ANIm and AF)
Writes data/real/sets/<set>/truth.tsv and genome_truth.tsv.

Usage: python3 bench_real/truth.py [sets...]   (MUMMER=<bin dir>, JOBS=<n>)
"""
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
MUMMER = os.environ.get("MUMMER", os.path.join(ROOT, "data", "envs", "mummer", "bin"))
JOBS = int(os.environ.get("JOBS", "48"))
SETS = os.environ.get("REAL_SETS", "bcereus,ecoli,mag2328,paer,saur,negative").split(",")


def contig_lengths(path):
    L, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            L[name] = 0
        else:
            L[name] += len(line.strip())
    return L


def coords(delta, flt):
    """[(qname, qstart, qend, len_q, idy)] from delta-filter <flt> | show-coords."""
    f = subprocess.run([os.path.join(MUMMER, "delta-filter"), flt, delta], capture_output=True, check=True, text=True).stdout
    with tempfile.NamedTemporaryFile("w", suffix=".delta", delete=False) as fh:
        fh.write(f)
        tmp = fh.name
    out = subprocess.run([os.path.join(MUMMER, "show-coords"), "-T", "-H", "-l", tmp], capture_output=True, check=True, text=True).stdout
    os.unlink(tmp)
    rows = []
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 11:
            continue
        s2, e2 = int(p[2]), int(p[3])
        rows.append((p[10], min(s2, e2), max(s2, e2), int(p[5]), float(p[6]), int(p[4])))
    return rows


def union_len(iv):
    iv.sort()
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce + 1:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        tot += ce - cs + 1
    return tot


def one(job):
    set_, gid, ref = job
    q = os.path.join(R, "sets", set_, "genomes", gid + ".fa")
    with tempfile.TemporaryDirectory() as td:
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre + "m", ref, q],
                       check=True, capture_output=True)
        subprocess.run([os.path.join(MUMMER, "nucmer"), "-t", "1", "--prefix", pre, ref, q],
                       check=True, capture_output=True)
        cq = coords(pre + "m.delta", "-q")
        c1 = coords(pre + ".delta", "-1")
    L = contig_lengths(q)
    iv, wid, wl, na = {}, {}, {}, {}
    for name, s, e, lq, idy, _ in cq:
        iv.setdefault(name, []).append((s, e))
        wid[name] = wid.get(name, 0.0) + idy * lq
        wl[name] = wl.get(name, 0) + lq
        na[name] = na.get(name, 0) + 1
    rows = []
    for name, ln in L.items():
        af = min(1.0, union_len(iv.get(name, [])) / ln)
        ani = wid[name] / wl[name] / 100 if wl.get(name) else float("nan")
        rows.append((name, ln, af, ani, wl.get(name, 0), na.get(name, 0)))
    gl = sum(L.values())
    g_iv = {}
    gw = gwl = 0.0
    for name, s, e, lq, idy, _ in c1:
        g_iv.setdefault(name, []).append((s, e))
        gw += idy * lq
        gwl += lq
    g_al = sum(union_len(v) for v in g_iv.values())
    grow = (gid, gl, g_al / gl, gw / gwl / 100 if gwl else float("nan"))
    return set_, rows, grow


def main():
    sets = sys.argv[1:] or SETS
    jobs = []
    for s in sets:
        for line in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]:
            p = line.rstrip("\n").split("\t")
            jobs.append((s, p[0], os.path.join(R, "refs", p[4] + ".fa")))
    res = {s: ([], []) for s in sets}
    with Pool(JOBS) as pool:
        for k, (s, rows, grow) in enumerate(pool.imap_unordered(one, jobs, chunksize=1)):
            res[s][0].extend(rows)
            res[s][1].append(grow)
            if (k + 1) % 20 == 0:
                print(f"{k + 1}/{len(jobs)}", file=sys.stderr, flush=True)
    for s, (rows, grows) in res.items():
        rows.sort()
        grows.sort()
        with open(os.path.join(R, "sets", s, "truth.tsv"), "w") as fh:
            fh.write("contig\tlength\ttrue_af\ttrue_ani\taligned_bp\tn_aln\n")
            for r in rows:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.5f}\t{r[3]:.5f}\t{r[4]}\t{r[5]}\n")
        with open(os.path.join(R, "sets", s, "genome_truth.tsv"), "w") as fh:
            fh.write("gid\tbases\tanim_af\tanim_ani\n")
            for g in grows:
                fh.write(f"{g[0]}\t{g[1]}\t{g[2]:.5f}\t{g[3]:.5f}\n")
        print(s, len(rows), "contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
