#!/usr/bin/env python3
"""Run aif and skani configurations on every real set (accuracy runs; timing is in time_real.py).

For each set and reference: aif sketch + aif query (1 thread), skani dist --qi (1 thread).
Raw outputs go to data/real/results/<config>/<set>.<ref>.tsv; the common prediction file
data/real/results/<config>/<set>.pred.tsv has columns contig, ani (fraction, nan when the
tool reports nothing), af (fraction of the contig aligned; 0 when nothing is reported),
reported (0/1). Every query contig appears in it.

Usage: python3 bench_real/run_tools.py [config ...]
"""
import math
import os
import subprocess
import sys
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
SETS = os.environ.get("REAL_SETS", "bcereus,ecoli,mag2328,paer,saur,negative").split(",")
CONFIGS = {
    "aif": ("aif", [], []),
    "aif_0.44": ("aif", ["--budget", "0.44"], []),
    # real-data setting chosen on data/real_tune (bench_real/tune_real.py): lowest AF error
    # among settings whose ANI error is no worse than the default's
    "aif_rt": ("aif", [], ["--rho-u", "2e-5"]),
    "aif_0.44_rt": ("aif", ["--budget", "0.44"], ["--rho-u", "2e-5"]),
    # speed option: no fine-binned heterogeneity re-run for contigs with > 16384 probes
    "aif_fast": ("aif", [], ["--het-fine-max", "16384"]),
    # AF-oriented local identity prior (v1.4 option; tuned on data/real_tune)
    "aif_lp": ("aif", [], ["--ident-window", "0.08", "--ident-radius", "3"]),
    "skani": ("skani", [], None),
    "skani_relaxed": ("skani", ["-s", "70", "--min-af", "5"], None),
    "skani_slow": ("skani", ["--slow"], None),
}


def names(fa):
    return [l[1:].split()[0] for l in open(fa) if l.startswith(">")]


def timed(cmd, stdout=None):
    tf = stdout + ".time" if stdout else None
    with open(stdout, "w") if stdout else open(os.devnull, "w") as out:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
    t, m = open(tf).read().split()[-2:]
    return float(t), int(m)


def run(job):
    cfg, set_, ref = job
    tool, sargs, qargs = CONFIGS[cfg]
    out = os.path.join(R, "results", cfg)
    q = os.path.join(R, "sets", set_, f"{ref}.contigs.fa")
    rfa = os.path.join(R, "refs", ref + ".fa")
    raw = os.path.join(out, f"{set_}.{ref}.tsv")
    pred = {}
    if tool == "aif":
        sk = os.path.join(out, f"{set_}.{ref}.aif")
        st, sm = timed([AIF, "sketch", rfa, "-o", sk] + sargs, os.path.join(out, f"{set_}.{ref}.sketch.log"))
        qt, qm = timed([AIF, "query", sk, q, "-t", "1"] + qargs, raw)
        for line in list(open(raw))[1:]:
            p = line.rstrip("\n").split("\t")
            ok = p[7] == "ok"
            pred[p[0]] = (float(p[1]) if ok else math.nan, float(p[2]) if ok else 0.0, int(ok))
        res = (st, sm, qt, qm, os.path.getsize(sk))
    else:
        qt, qm = timed([SKANI, "dist", "-t", "1", "--qi", "-q", q, "-r", rfa, "-o", raw] + sargs, raw + ".stdout")
        for line in list(open(raw))[1:]:
            p = line.rstrip("\n").split("\t")
            pred[p[6].split()[0]] = (float(p[2]) / 100, float(p[4]) / 100, 1)
        res = (0.0, 0, qt, qm, 0)
    rows = [(c,) + pred.get(c, (math.nan, 0.0, 0)) for c in names(q)]
    return cfg, set_, ref, rows, res


def main():
    cfgs = sys.argv[1:] or list(CONFIGS)
    jobs = []
    for cfg in cfgs:
        os.makedirs(os.path.join(R, "results", cfg), exist_ok=True)
        for s in SETS:
            refs = sorted(set(l.split("\t")[4] for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]))
            jobs += [(cfg, s, r) for r in refs]
    out = {}
    with Pool(min(16, len(jobs))) as pool:
        for cfg, s, ref, rows, res in pool.imap_unordered(run, jobs):
            out.setdefault((cfg, s), []).append((ref, rows, res))
    for (cfg, s), parts in sorted(out.items()):
        d = os.path.join(R, "results", cfg)
        with open(os.path.join(d, f"{s}.pred.tsv"), "w") as fh:
            fh.write("contig\tani\taf\treported\n")
            for _, rows, _ in sorted(parts):
                for c, a, f, rep in rows:
                    fh.write(f"{c}\t{a:.6f}\t{f:.6f}\t{rep}\n")
        with open(os.path.join(d, f"{s}.resources.tsv"), "w") as fh:
            fh.write("reference\tsketch_s\tsketch_kb\tquery_s\tquery_kb\tsketch_bytes\n")
            for ref, _, r in sorted(parts):
                fh.write(f"{ref}\t" + "\t".join(map(str, r)) + "\n")
        n = sum(len(p[1]) for p in parts)
        rep = sum(x[3] for p in parts for x in p[1])
        print(f"{cfg:14s} {s:9s} contigs={n} reported={rep}", file=sys.stderr)


if __name__ == "__main__":
    main()
