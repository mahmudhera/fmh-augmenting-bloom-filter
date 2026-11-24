#!/usr/bin/env python3
"""Kraken2 baseline on the same GTDB R232 species representatives as aif's database.

  python3 bench_real/kraken2_gtdb.py prepare   # GTDB taxonomy -> NCBI-style taxdump; library .fna files and
                                               # seqid2taxid.map from data/gtdb/meta/reps.list (no masking)
  python3 bench_real/kraken2_gtdb.py build     # kraken2-build --build --fast-build --block-size 16 MB, capped at
                                               # aif's database size (--max-db-size = size of gtdb_0.88.aifd), timed
  python3 bench_real/kraken2_gtdb.py build_step3  # resume: step 3 only, with --fast-build (-F), reusing the
                                                 # capacity estimate of a previous build's step 2 (build.log)
The reported database came from a deterministic build stopped in step 3, then build_step3 (see build.json).
  python3 bench_real/kraken2_gtdb.py classify  # both GTDB query sets (GTDB_SET "" and 2), timed
  python3 bench_real/kraken2_gtdb.py evaluate  # species / genus accuracy as in gtdb_bench.py

Kraken2 assigns each contig to a taxon (LCA of its k-mer hits); a contig counts as correct at species
(genus) level when the assigned taxon is that species (or a species of that genus / the genus itself).
"""
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
K = os.path.join(G, "kraken2_db")
THREADS = os.environ.get("THREADS", "48")
RANKS = [("d", "superkingdom"), ("p", "phylum"), ("c", "class"), ("o", "order"), ("f", "family"), ("g", "genus"), ("s", "species")]
NPART = 64


def timed(cmd, out=None, err=None):
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=fo, stderr=fe, check=True)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def taxonomy():
    """{ncbi accession: [d__X, p__Y, ..., s__Z]} for every GTDB genome"""
    T = {}
    for fn in ("bac120_taxonomy.tsv.gz", "ar53_taxonomy.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            for line in fh:
                acc, tax = line.rstrip("\n").split("\t")
                T[acc[3:]] = tax.split(";")
    return T


def write_part(job):
    i, items = job
    out = os.path.join(K, "library", "gtdb", f"part{i:02d}.fna")
    mp = []
    with open(out, "w") as fo:
        for path, taxid in items:
            with gzip.open(path, "rt") as fi:
                for line in fi:
                    if line.startswith(">"):
                        sid = line[1:].split()[0]
                        mp.append(f"{sid}\t{taxid}\n")
                        fo.write(f">{sid}\n")
                    else:
                        fo.write(line)
    return mp


def prepare():
    T = taxonomy()
    os.makedirs(os.path.join(K, "taxonomy"), exist_ok=True)
    os.makedirs(os.path.join(K, "library", "gtdb"), exist_ok=True)
    ids = {"root": 1}
    nodes = [(1, 1, "no rank")]
    names = [(1, "root")]

    def node(path):  # path: tuple of rank strings from the top
        if path in ids:
            return ids[path]
        parent = node(path[:-1]) if len(path) > 1 else 1
        tid = len(ids) + 1
        ids[path] = tid
        nodes.append((tid, parent, dict(RANKS)[path[-1][0]]))
        names.append((tid, path[-1]))
        return tid
    reps = [l.rstrip("\n").split("\t") for l in open(os.path.join(G, "meta", "reps.list"))]
    items = []
    for path, acc in reps:
        items.append((path, node(tuple(T[acc]))))
    for p in list(ids):  # every GTDB taxon of the query genomes too (taxonomy lookups at evaluation)
        pass
    with open(os.path.join(K, "taxonomy", "nodes.dmp"), "w") as fh:
        for tid, par, rank in nodes:
            fh.write(f"{tid}\t|\t{par}\t|\t{rank}\t|\t\t|\n")
    with open(os.path.join(K, "taxonomy", "names.dmp"), "w") as fh:
        for tid, name in names:
            fh.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    jobs = [(i, items[i::NPART]) for i in range(NPART)]
    with Pool(int(THREADS)) as pool, open(os.path.join(K, "seqid2taxid.map"), "w") as fm:
        for k, mp in enumerate(pool.imap_unordered(write_part, jobs)):
            fm.writelines(mp)
            print(f"part {k + 1}/{NPART}", flush=True)
    json.dump({"taxa": len(ids), "genomes": len(items)}, open(os.path.join(K, "prepare.json"), "w"))
    print(f"{len(items)} genomes, {len(ids)} taxa")


def build():
    cap = os.path.getsize(os.path.join(G, "aif", "gtdb_0.88.aifd"))
    # --fast-build and 16 MB read blocks, as in build_step3 (the default deterministic step 3 is ~50x slower here)
    block = int(os.environ.get("K2_BLOCK", 16 * 1024 * 1024))
    t, m = timed(["kraken2-build", "--build", "--db", K, "--threads", THREADS, "--max-db-size", str(cap), "--fast-build",
                  "--block-size", str(block)], out=os.path.join(K, "build.log"), err=os.path.join(K, "build.err"))
    size = sum(os.path.getsize(os.path.join(K, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    r = {"build_s": t, "build_kb": m, "bytes": size, "max_db_size": cap, "threads": int(THREADS), "fast_build": True, "block_size": block}
    json.dump(r, open(os.path.join(K, "build.json"), "w"), indent=1)
    print(r)


def build_step3():
    """Kraken2's default (deterministic) step 3 parallelises within one sequence at a time and ran at ~2 MB/s on
    GTDB's many short contigs (29.7 GB of 618 GB in 4.2 h), so it is rerun with Kraken2's --fast-build mode
    (build_db -F: parallel over sequences; minimizer/LCA pairs may vary between runs). The arguments are those
    that build_kraken2_db.sh passed (capacity from its step 2, -M from --max-db-size)."""
    import re
    log = open(os.path.join(K, "build.log")).read()
    req = int(re.search(r"Estimated hash table requirement: (\d+) bytes", log).group(1))
    h, mi, sec = re.search(r"Capacity estimation complete. \[(?:(\d+)h)?(?:(\d+)m)?([\d.]+)s\]", log).groups()
    step2 = int(h or 0) * 3600 + int(mi or 0) * 60 + float(sec)
    cap = os.path.getsize(os.path.join(G, "aif", "gtdb_0.88.aifd"))
    # kraken2-build --block-size (default 16 KB): with 16 KB blocks the fast build's threads queue on the input lock
    # (~26 MB/s); 16 MB blocks remove that bottleneck. Sub-block size as kraken2-build sets it: block / threads.
    BLOCK = int(os.environ.get("K2_BLOCK", 16 * 1024 * 1024))
    files = sorted(os.path.join("library", "gtdb", f) for f in os.listdir(os.path.join(K, "library", "gtdb")) if f.endswith(".fna"))
    for f in ("hash.k2d.tmp", "taxo.k2d.tmp", "opts.k2d.tmp"):
        if os.path.exists(os.path.join(K, f)):
            os.unlink(os.path.join(K, f))
    cmd = (f"cat {' '.join(files)} | /usr/lib/kraken2/build_db -k 35 -l 31 -S 1111111111111111101010101010101 -H hash.k2d.tmp -t taxo.k2d.tmp "
           f"-o opts.k2d.tmp -n taxonomy/ -m seqid2taxid.map -c {req // 4} -p {THREADS} -M {cap // 4} -B {BLOCK} -b {-(-BLOCK // int(THREADS))} -r 0 -F")
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(os.path.join(K, "build3.err"), "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf, "bash", "-c", cmd], cwd=K, stderr=fe, stdout=fe, check=True)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    for f in ("taxo.k2d", "opts.k2d", "hash.k2d"):
        os.rename(os.path.join(K, f + ".tmp"), os.path.join(K, f))
    size = sum(os.path.getsize(os.path.join(K, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    r = {"build_s": step2 + float(t), "step2_estimate_s": step2, "step3_s": float(t), "build_kb": int(m), "bytes": size,
         "max_db_size": cap, "required_bytes_uncapped": req, "threads": int(THREADS), "fast_build": True, "block_size": BLOCK,
         "note": "step 2 (capacity estimate) from the first build; step 3 rerun with --fast-build after the default "
                 "deterministic step 3 processed 29.7 GB of 618 GB in 4.2 h (and --fast-build with the default 16 KB blocks ran "
                 "at 26 MB/s); peak memory is that of step 3"}
    json.dump(r, open(os.path.join(K, "build.json"), "w"), indent=1)
    print(r)


def classify():
    out = {}
    tp = os.path.join(K, "classify_time.json")
    if os.path.exists(tp):
        out = json.load(open(tp))
    for s in ("", "2"):
        q = os.path.join(G, f"queries{s}", "contigs.fa")
        key = f"set{s or '1'}"
        if key in out:
            continue
        subprocess.run(["cat", os.path.join(K, "hash.k2d")], stdout=subprocess.DEVNULL, check=True)  # warm cache
        t, m = timed(["kraken2", "--db", K, "--threads", os.environ.get("SEARCH_THREADS", "16"), "--output",
                      os.path.join(K, f"out_{key}.tsv"), "--report", os.path.join(K, f"report_{key}.txt"), q])
        out[key] = {"s": t, "kb": m}
        print(key, out[key], flush=True)
        json.dump(out, open(tp, "w"), indent=1)


def evaluate():
    nodes, names = {}, {}
    for line in open(os.path.join(K, "taxonomy", "nodes.dmp")):
        p = [x.strip() for x in line.split("|")]
        nodes[int(p[0])] = (int(p[1]), p[2])
    for line in open(os.path.join(K, "taxonomy", "names.dmp")):
        p = [x.strip() for x in line.split("|")]
        names[int(p[0])] = p[1]

    def lineage(t):
        L = {}
        while t != 1 and t in nodes:
            L[nodes[t][1]] = names[t]
            t = nodes[t][0]
        return L
    res = {}
    LB = [(500, 1000, "0.5-1k"), (1000, 2000, "1-2k"), (2000, 5000, "2-5k"), (5000, 20000, "5-20k"), (20000, 10 ** 12, ">=20k")]
    for s in ("", "2"):
        key = f"set{s or '1'}"
        meta = {r["contig"]: r for r in csv.DictReader(open(os.path.join(G, "meta", f"query_contigs{s}.tsv")), delimiter="\t")}
        A = {}
        for line in open(os.path.join(K, f"out_{key}.tsv")):
            p = line.rstrip("\n").split("\t")
            A[p[1]] = int(p[2]) if p[0] == "C" else 0
        recs = []
        for c, m in meta.items():
            t = A.get(c, 0)
            L = lineage(t) if t else {}
            recs.append({"L": int(m["length"]), "dom": m["domain"], "rep": int(t != 0), "at_species": int("species" in L),
                         "sp": int(L.get("species", "")[3:] == m["species"]), "gen": int(L.get("genus", "")[3:] == m["genus"])})

        def summ(rs):
            n = max(1, len(rs))
            rep = [r for r in rs if r["rep"]]
            sp_calls = [r for r in rs if r["at_species"]]
            return {"n": len(rs), "classified": len(rep) / n, "species_acc": sum(r["sp"] for r in rs) / n, "genus_acc": sum(r["gen"] for r in rs) / n,
                    "species_calls": len(sp_calls) / n, "species_prec": sum(r["sp"] for r in sp_calls) / max(1, len(sp_calls))}
        d = {"all": summ(recs), "bac": summ([r for r in recs if r["dom"] == "bac"]), "ar": summ([r for r in recs if r["dom"] == "ar"])}
        for lo, hi, lab in LB:
            d[lab] = summ([r for r in recs if lo <= r["L"] < hi])
        res[key] = d
        print(key, {k: round(v, 4) for k, v in d["all"].items()}, {lab: round(d[lab]["species_acc"], 3) for _, _, lab in LB})
    for f in ("build.json", "classify_time.json", "prepare.json"):
        p = os.path.join(K, f)
        if os.path.exists(p):
            res[f[:-5]] = json.load(open(p))
    json.dump(res, open(os.path.join(G, "kraken2_eval.json"), "w"), indent=1)


if __name__ == "__main__":
    {"prepare": prepare, "build": build, "build_step3": build_step3, "classify": classify, "evaluate": evaluate}[sys.argv[1]]()
