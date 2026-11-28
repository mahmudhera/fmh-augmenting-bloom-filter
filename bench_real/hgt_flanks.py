#!/usr/bin/env python3
"""Flank placement for the HGT atlas: is a shared block embedded in its host genome, or a foreign contig?

Contamination and mis-binning put whole contigs of one organism into the assembly of another. A block shared with
partner P that sits on such a contig is not transfer. For every block of >= 1 kb aligned bases (results_hgt/TAG/
blocks.tsv.gz), the windows immediately left and right of it on its contig (the flanks; by construction not
shared with P at >= 99%) are searched against all representatives without the 99% filter, and the closest relative of
each flank is kept. hgt_analyze.py (HGT_FLANKS) then rejects a block in a draft or metagenome-assembled genome Q when
it has no flank (the block is the whole contig) or when a flank's closest relative is nearer to P's lineage than to
Q's (the contig belongs to P's side of the tree).

Stages (data in data/hgt/SET/flanks/):
  extract SET TAG   flank windows of TAG's blocks into FASTA parts (gzip level 1)
  search SET        aif db search of the parts -> fhits/partNNN.tsv (filter stage on the fly)
  filter            stdin -> stdout: per window, the 10 best non-self hits (ANI >= 0.80, >= 500 aligned bp) by AF x ANI
  calibrate SET TAG flank evidence of TAG's cross-genus blocks by genome category -> results_hgt/TAG/flank_calibration.json
                    (complete genomes cannot carry foreign contigs, so their evidence is the background of real transfer)

Usage: /usr/bin/python3 -I bench_real/hgt_flanks.py extract all all; ... search all
"""
import collections
import csv
import gzip
import json
import multiprocessing
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
DB = os.path.join(ROOT, "data", "gtdb", "aif", "gtdb_0.88.aifd")
AIF = os.path.join(ROOT, "bin", "aif")
W, WMIN = 5000, 1000
PAIR_BP = 1000
THREADS = os.environ.get("THREADS", "48")
NCAND = os.environ.get("HGT_NCAND", "128")
ANI_KEEP, BP_KEEP, TOP = 0.80, 500, 10
PART_BP = 4e9


def read_genome(path):
    seqs, buf = [], []
    with gzip.open(path, "rt") as fi:
        for line in fi:
            if line[0] == ">":
                if buf:
                    seqs.append("".join(buf))
                    buf = []
                seqs.append(None)
            else:
                buf.append(line.strip())
    if buf:
        seqs.append("".join(buf))
    return [s for s in seqs if s is not None]


def write_part(job):
    p, items, outd = job
    path = os.path.join(outd, f"part{p:03d}.fa.gz")
    n = 0
    with gzip.open(path + ".tmp", "wt", compresslevel=1) as fo:
        for acc, gpath, fl in items:
            seqs = read_genome(gpath)
            for ci, st, ln in sorted(fl):
                t = seqs[ci][st:st + ln]
                if len(t) == ln:
                    fo.write(f">{acc}|{ci}|{st}|{ln}\n{t}\n")
                    n += 1
    os.replace(path + ".tmp", path)
    return p, n


def extract(setname, tag):
    base = os.path.join(HD, setname)
    outd = os.path.join(base, "flanks", "win")
    os.makedirs(outd, exist_ok=True)
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    blocks = collections.defaultdict(list)
    for b in csv.DictReader(gzip.open(os.path.join(ROOT, "results_hgt", tag, "blocks.tsv.gz"), "rt"), delimiter="\t"):
        if float(b["aligned_bp"]) >= PAIR_BP:
            blocks[b["query"]].append((int(b["contig"]), int(b["start"]), int(b["end"])))
    clen = {}
    for l in open(os.path.join(base, "contigs.tsv")):
        acc, ci, nm, ln, nw = l.split("\t", 4)
        if acc in blocks:
            clen[(acc, int(ci))] = int(ln)
    flanks = {}
    nblk = noflank = 0
    for acc, bl in blocks.items():
        fl = set()
        for ci, s0, e in bl:
            nblk += 1
            ln = clen[(acc, ci)]
            k = 0
            if s0 >= W:
                fl.add((ci, s0 - W, W))
                k += 1
            if ln - e >= WMIN:
                fl.add((ci, e, min(W, ln - e)))
                k += 1
            noflank += k == 0
        if fl:
            flanks[acc] = fl
    tot = sum(x[2] for v in flanks.values() for x in v)
    nparts = max(1, round(tot / PART_BP))
    parts = [[] for _ in range(nparts)]
    load = [0] * nparts
    for acc in sorted(flanks, key=lambda a: -sum(x[2] for x in flanks[a])):
        i = load.index(min(load))
        parts[i].append((acc, G[acc]["path"], flanks[acc]))
        load[i] += sum(x[2] for x in flanks[acc])
    info = {"tag": tag, "blocks": nblk, "blocks_without_flank": noflank, "genomes": len(flanks),
            "flank_windows": sum(len(v) for v in flanks.values()), "flank_bp": tot, "parts": nparts}
    print(json.dumps(info), flush=True)
    jobs = [(p, parts[p], outd) for p in range(nparts) if not os.path.exists(os.path.join(outd, f"part{p:03d}.fa.gz"))]
    with multiprocessing.get_context("fork").Pool(min(32, max(1, len(jobs)))) as pool:
        for p, n in pool.imap_unordered(write_part, jobs):
            print(f"part {p} written ({n} windows)", flush=True)
    json.dump(info, open(os.path.join(base, "flanks", "extract.json"), "w"), indent=1)


def filt():
    """stdin: aif db search output; stdout: window, ref, ani, af for the TOP best hits of each window by AF x ANI"""
    out = sys.stdout
    cur, buf = None, []

    def flush():
        buf.sort(key=lambda x: -x[0])
        for sc, ref, ani, af in buf[:TOP]:
            out.write(f"{cur}\t{ref}\t{ani}\t{af}\n")
    for line in sys.stdin:
        p = line.split("\t", 12)
        if p[0] == "contig" or len(p) < 12 or p[9] != "ok":
            continue
        q = p[0].split("|")
        if q[0] == p[1]:
            continue
        ani, af = float(p[3]), float(p[4])
        if ani < ANI_KEEP or af * float(q[3]) < BP_KEEP:
            continue
        if p[0] != cur:
            if cur is not None:
                flush()
            cur, buf = p[0], []
        buf.append((af * ani, p[1], p[3], p[4]))
    if cur is not None:
        flush()


def search(setname):
    base = os.path.join(HD, setname, "flanks")
    wd, hd = os.path.join(base, "win"), os.path.join(base, "fhits")
    os.makedirs(hd, exist_ok=True)
    tlog = os.path.join(base, "search_time.json")
    times = json.load(open(tlog)) if os.path.exists(tlog) else {}
    for f in sorted(os.listdir(wd)):
        if not f.endswith(".fa.gz"):
            continue
        out = os.path.join(hd, f.split(".")[0] + ".tsv")
        if os.path.exists(out):
            continue
        t0 = time.time()
        tf = out + ".time"
        cmd = (f"/usr/bin/time -f '%e %M' -o {tf} {AIF} db search {DB} {os.path.join(wd, f)} -t {THREADS} "
               f"-n {NCAND} --max-cand {NCAND} | /usr/bin/python3 -I {os.path.abspath(__file__)} filter > {out}.tmp")
        subprocess.run(cmd, shell=True, check=True)
        os.replace(out + ".tmp", out)
        e, m = open(tf).read().split()
        times[f] = {"wall": time.time() - t0, "aif_s": float(e), "kb": int(m),
                    "bytes": os.path.getsize(os.path.join(wd, f))}
        json.dump(times, open(tlog, "w"), indent=1)
        print(f, times[f], flush=True)


def calibrate(setname, tag):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hgt_analyze import lca_class
    base = os.path.join(HD, setname)
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    for g in G.values():
        g["tax"] = [x[3:] for x in g["taxonomy"].split(";")]
    comp = {a for a, g in G.items() if g["category"] == "isolate" and g["assembly_level"] in ("Complete Genome", "Chromosome")}
    FL = {}
    hd = os.path.join(base, "flanks", "fhits")
    for f in sorted(os.listdir(hd)):
        if f.endswith(".tsv"):
            for l in open(os.path.join(hd, f)):
                w, ref, ani, af = l.rstrip("\n").split("\t")
                a, ci, st, ln = w.split("|")
                FL.setdefault((a, int(ci), int(st)), (ref, float(ani), float(af)))  # best hit first
    B = []
    for b in csv.DictReader(gzip.open(os.path.join(ROOT, "results_hgt", tag, "blocks.tsv.gz"), "rt"), delimiter="\t"):
        if b["distance"] != "same genus" and float(b["aligned_bp"]) >= PAIR_BP:
            B.append((b["query"], b["partner"], int(b["contig"]), int(b["start"]), int(b["end"])))
    need = {b[0] for b in B}
    clen = {}
    for l in open(os.path.join(base, "contigs.tsv")):
        acc, ci, nm, ln, nw = l.split("\t", 4)
        if acc in need:
            clen[(acc, int(ci))] = int(ln)

    def cat(a):
        return "complete" if a in comp else "draft isolate" if G[a]["category"] == "isolate" else "MAG or SAG"
    order = ["partner itself", "partner lineage", "neither", "host lineage", "unplaced"]
    res = collections.defaultdict(collections.Counter)
    pb = collections.defaultdict(collections.Counter)
    for q, p, ci, s0, e in B:
        c = cat(q)
        keys = ([(q, ci, s0 - W)] if s0 >= W else []) + ([(q, ci, e)] if clen[(q, ci)] - e >= WMIN else [])
        if not keys:
            res[c]["no flank"] += 1
            continue
        worst = "unplaced"
        for k in keys:
            h = FL.get(k)
            if not h:
                continue
            ref, ani, af = h
            if ref == p:
                s_ = "partner itself"
                pb[c][("ani>=0.97" if ani >= 0.97 else "ani<0.97") + ("|af>=0.8" if af >= 0.8 else "|af<0.8")] += 1
            else:
                dq, dp = lca_class(G[q]["tax"], G[ref]["tax"]), lca_class(G[p]["tax"], G[ref]["tax"])
                s_ = "partner lineage" if dp < dq else "host lineage" if dq < dp else "neither"
            if order.index(s_) < order.index(worst):
                worst = s_
        res[c][worst] += 1
    out = {c: {"blocks": sum(v.values()), **dict(v), "partner_itself_flanks": dict(pb[c])} for c, v in res.items()}
    json.dump(out, open(os.path.join(ROOT, "results_hgt", tag, "flank_calibration.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


def main():
    st = sys.argv[1] if len(sys.argv) > 1 else ""
    if st == "extract":
        extract(sys.argv[2], sys.argv[3])
    elif st == "search":
        search(sys.argv[2])
    elif st == "filter":
        filt()
    elif st == "calibrate":
        calibrate(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
