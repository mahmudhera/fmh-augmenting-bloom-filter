#!/usr/bin/env python3
"""Interleaved, core-pinned timing of skani and aif on the real sets.

Per set (one reference each): skani dist --qi (sketch + query, 1 thread), aif sketch and
aif query (1 thread; also with the --het-fine-max 16384 speed option); then 8-thread
runs. Median wall time over rounds after a warm-up round; peak RSS from one GNU time
run. Also records sketch sizes (bits per reference base).
Writes data/real/results/timing_real.tsv and sketch_sizes.tsv.

Usage: python3 bench_real/time_real.py [rounds]
"""
import os
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
PAIRS = [("bcereus", "BA52G"), ("ecoli", "MG1655"), ("mag2328", "REFMAG"), ("paer", "PAO1"), ("saur", "NCTC8325")]
C1 = ["taskset", "-c", "7"]
C8 = ["taskset", "-c", "8-15"]


def wall(cmd):
    t = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return time.perf_counter() - t


def rss_kb(cmd):
    with tempfile.NamedTemporaryFile("r") as tf:
        subprocess.run(["/usr/bin/time", "-f", "%M", "-o", tf.name] + cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        return int(open(tf.name).read().split()[-1])


def fasta_stats(path):
    n = b = 0
    for line in open(path):
        if line.startswith(">"):
            n += 1
        else:
            b += len(line) - 1
    return n, b


def main():
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    tmp = tempfile.mkdtemp()
    rows, sizes = [], []
    for s, ref in PAIRS:
        rfa = os.path.join(R, "refs", ref + ".fa")
        q = os.path.join(R, "sets", s, f"{ref}.contigs.fa")
        sk = os.path.join(tmp, ref + ".aif")
        subprocess.run([AIF, "sketch", rfa, "-o", sk], check=True, stderr=subprocess.DEVNULL)
        skd = os.path.join(tmp, ref + "_skani")
        subprocess.run([SKANI, "sketch", "-t", "1", rfa, "-o", skd], check=True, capture_output=True)
        rbases = fasta_stats(rfa)[1]
        skb = sum(os.path.getsize(os.path.join(skd, x)) for x in os.listdir(skd))
        sizes.append((ref, rbases, os.path.getsize(sk), skb))
        nq, bq = fasta_stats(q)
        cmds = {
            ("skani", "sketch+query", 1): C1 + [SKANI, "dist", "-t", "1", "--qi", "-q", q, "-r", rfa, "-o", os.devnull],
            ("aif", "sketch", 1): C1 + [AIF, "sketch", rfa, "-o", os.path.join(tmp, "x.aif")],
            ("aif", "query", 1): C1 + [AIF, "query", sk, q, "-t", "1"],
            ("aif", "query-fast", 1): C1 + [AIF, "query", sk, q, "-t", "1", "--het-fine-max", "16384"],
            ("skani", "sketch+query", 8): C8 + [SKANI, "dist", "-t", "8", "--qi", "-q", q, "-r", rfa, "-o", os.devnull],
            ("aif", "query", 8): C8 + [AIF, "query", sk, q, "-t", "8"],
            ("aif", "query-fast", 8): C8 + [AIF, "query", sk, q, "-t", "8", "--het-fine-max", "16384"],
        }
        for c in cmds.values():
            wall(c)
        times = {k: [] for k in cmds}
        for _ in range(rounds):
            for k, c in cmds.items():
                times[k].append(wall(c))
        for k, v in times.items():
            rows.append((s, ref, nq, bq) + k + (statistics.median(v), min(v), rss_kb(cmds[k])))
        print(s, {f"{k[0]}:{k[1]}:{k[2]}": f"{statistics.median(v):.2f}s" for k, v in times.items()}, flush=True)
    out = os.path.join(R, "results")
    with open(os.path.join(out, "timing_real.tsv"), "w") as fh:
        fh.write("set\treference\tcontigs\tbases\ttool\tphase\tthreads\tmedian_s\tmin_s\tmax_rss_kb\n")
        for r in rows:
            fh.write("\t".join(map(str, r[:7])) + f"\t{r[7]:.4f}\t{r[8]:.4f}\t{r[9]}\n")
    with open(os.path.join(out, "sketch_sizes.tsv"), "w") as fh:
        fh.write("reference\tbases\taif_bytes\tskani_bytes\taif_bits_per_base\tskani_bits_per_base\n")
        for ref, b, a, k in sizes:
            fh.write(f"{ref}\t{b}\t{a}\t{k}\t{a * 8 / b:.4f}\t{k * 8 / b:.4f}\n")


if __name__ == "__main__":
    main()
