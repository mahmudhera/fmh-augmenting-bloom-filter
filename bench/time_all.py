#!/usr/bin/env python3
"""Interleaved, pinned timing of skani and aif (median wall time, peak RSS).

Every round runs each command once (alternating tools, same core), after warm-up
rounds, so both tools see the same CPU frequency state. Wall time is measured with
perf_counter (ms resolution); peak RSS comes from one extra GNU time run per command.

Usage: python3 bench/time_all.py [rounds] [aif_binary] [out_dir]
Writes <out_dir>/timing.tsv and <out_dir>/throughput.tsv.
"""
import os
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
GENOMES = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]


def wall(cmd):
    t = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return time.perf_counter() - t


def rss_kb(cmd):
    with tempfile.NamedTemporaryFile("r") as tf:
        subprocess.run(["/usr/bin/time", "-f", "%M", "-o", tf.name] + cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        return int(open(tf.name).read().split()[-1])


def bench(cmds, rounds):
    for _ in range(3):
        for c in cmds.values():
            wall(c)
    times = {k: [] for k in cmds}
    for _ in range(rounds):
        for k, c in cmds.items():
            times[k].append(wall(c))
    return {k: (statistics.median(v), min(v), rss_kb(cmds[k])) for k, v in times.items()}


def main():
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    aif = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "bin", "aif")
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "data", "results", "main")
    os.makedirs(out, exist_ok=True)
    tmp = tempfile.mkdtemp()
    core = ["taskset", "-c", "7"]
    rows = []
    for g in GENOMES:
        ref = os.path.join(ROOT, "data", "genomes", g + ".fa")
        q = os.path.join(ROOT, "data", "bench", g + ".contigs.fa")
        sk = os.path.join(tmp, g + ".aif")
        subprocess.run([aif, "sketch", ref, "-o", sk], check=True, stderr=subprocess.DEVNULL)
        res = bench({
            "skani_total": core + [SKANI, "dist", "-t", "1", "--qi", "-q", q, "-r", ref, "-o", os.devnull],
            "aif_sketch": core + [aif, "sketch", ref, "-o", os.path.join(tmp, "x.aif")],
            "aif_query": core + [aif, "query", sk, q],
        }, rounds)
        for k, (med, mn, mem) in res.items():
            tool, phase = k.split("_", 1)
            rows.append((tool, g, phase, med, mn, mem))
        print(g, {k: f"{v[0]*1000:.0f} ms" for k, v in res.items()}, flush=True)
    with open(os.path.join(out, "timing.tsv"), "w") as fh:
        fh.write("tool\tgenome\tphase\tmedian_s\tmin_s\tmax_rss_kb\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]:.4f}\t{r[5]}\n")
    # throughput: 2,800 contigs (20.3 Mbp) against E. coli, 1 and 8 threads
    q = os.path.join(ROOT, "data", "throughput", "NC_000913.3.contigs.fa")
    ref = os.path.join(ROOT, "data", "genomes", "NC_000913.3.fa")
    sk = os.path.join(tmp, "thr.aif")
    subprocess.run([aif, "sketch", ref, "-o", sk], check=True, stderr=subprocess.DEVNULL)
    nbases = sum(len(l.strip()) for l in open(q) if not l.startswith(">"))
    ncontigs = sum(1 for l in open(q) if l.startswith(">"))
    cores8 = ["taskset", "-c", "8-15"]
    cmds = {
        "skani:sketch+query:1": core + [SKANI, "dist", "-t", "1", "--qi", "-q", q, "-r", ref, "-o", os.devnull],
        "aif:sketch:1": core + [aif, "sketch", ref, "-o", os.path.join(tmp, "y.aif")],
        "aif:query:1": core + [aif, "query", sk, q],
        "skani:sketch+query:8": cores8 + [SKANI, "dist", "-t", "8", "--qi", "-q", q, "-r", ref, "-o", os.devnull],
    }
    mt = subprocess.run([aif, "query", sk, q, "-t", "2"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
    if mt:
        cmds["aif:query:8"] = cores8 + [aif, "query", sk, q, "-t", "8"]
    res = bench(cmds, rounds)
    # 4x larger query set (81 Mbp): amortises skani's reference sketching
    q4 = os.path.join(ROOT, "data", "throughput", "NC_000913.3.x4.contigs.fa")
    cmds4 = {
        "skani:sketch+query:1": core + [SKANI, "dist", "-t", "1", "--qi", "-q", q4, "-r", ref, "-o", os.devnull],
        "aif:query:1": core + [aif, "query", sk, q4],
        "skani:sketch+query:8": cores8 + [SKANI, "dist", "-t", "8", "--qi", "-q", q4, "-r", ref, "-o", os.devnull],
    }
    if mt:
        cmds4["aif:query:8"] = cores8 + [aif, "query", sk, q4, "-t", "8"]
    res4 = bench(cmds4, max(3, rounds // 2)) if os.path.exists(q4) else {}
    with open(os.path.join(out, "throughput.tsv"), "w") as fh:
        fh.write("set\ttool\tphase\tthreads\tcontigs\tbases\tmedian_s\tmin_s\tmax_rss_kb\n")
        for label, rr, mult in (("x1", res, 1), ("x4", res4, 4)):
            for k, (med, mn, mem) in rr.items():
                tool, phase, thr = k.split(":")
                fh.write(f"{label}\t{tool}\t{phase}\t{thr}\t{ncontigs*mult}\t{nbases*mult}\t{med:.4f}\t{mn:.4f}\t{mem}\n")
                print(f"{label} {k:24s} median {med*1000:7.1f} ms  min {mn*1000:7.1f} ms  rss {mem} KB")


if __name__ == "__main__":
    main()
