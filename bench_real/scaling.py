#!/usr/bin/env python3
"""Thread scaling of aif query on the largest real set (B. cereus group, 449 Mbp of contigs).

Compares two aif binaries (e.g. v1.1 built from git and the current one) at 1-16 threads,
interleaved, median of rounds; peak RSS from one GNU time run.

Usage: python3 bench_real/scaling.py <aif_old> <aif_new> [rounds]  -> data/real/results/scaling.tsv
"""
import os
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "real")
CORES = {1: "8", 2: "8-9", 4: "8-11", 8: "8-15", 16: "8-23"}


def wall(cmd):
    t = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return time.perf_counter() - t


def rss_kb(cmd):
    with tempfile.NamedTemporaryFile("r") as tf:
        subprocess.run(["/usr/bin/time", "-f", "%M", "-o", tf.name] + cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        return int(open(tf.name).read().split()[-1])


def main():
    old, new = sys.argv[1], sys.argv[2]
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    tmp = tempfile.mkdtemp()
    sk = os.path.join(tmp, "BA52G.aif")
    subprocess.run([new, "sketch", os.path.join(R, "refs", "BA52G.fa"), "-o", sk], check=True, stderr=subprocess.DEVNULL)
    q = os.path.join(R, "sets", "bcereus", "BA52G.contigs.fa")
    cmds = {(tag, t): ["taskset", "-c", c, b, "query", sk, q, "-t", str(t)]
            for t, c in CORES.items() for tag, b in (("v1.1", old), ("v1.2", new))}
    for c in cmds.values():
        wall(c)
    times = {k: [] for k in cmds}
    for _ in range(rounds):
        for k, c in cmds.items():
            times[k].append(wall(c))
    with open(os.path.join(R, "results", "scaling.tsv"), "w") as fh:
        fh.write("version\tthreads\tmedian_s\tmin_s\tmax_rss_kb\n")
        for (tag, t), v in sorted(times.items(), key=lambda x: (x[0][0], x[0][1])):
            m = rss_kb(cmds[(tag, t)])
            fh.write(f"{tag}\t{t}\t{statistics.median(v):.4f}\t{min(v):.4f}\t{m}\n")
            print(f"{tag} t={t:2d} median {statistics.median(v):.3f} s  rss {m/1024:.0f} MB", flush=True)


if __name__ == "__main__":
    main()
