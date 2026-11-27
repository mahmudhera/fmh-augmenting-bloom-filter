#!/usr/bin/env python3
"""UHGG contig benchmark: aif against fast and alignment-based tools, with MUMmer as the reference standard.

Question: on real UHGG v3.0 contigs, sampled where the screen's decisions are made, how close are the fast
tools to alignment, and what would each cost on the whole catalogue (55.7 M contigs, 707 GB)?

Design (fixed before any tool was run on the sample):
  population  held-out UHGG genomes (not the 3,000-genome pilot), contigs >= 500 bp
  stage 1     8 of the 128 screen parts at random (a hash partition of the genomes), then 1,000 MAGs and 300
              isolates at random among their held-out genomes (isolates over-sampled); every contig of a drawn
              genome forms the pool
  pool calls  aif, skani and Kraken2 are run on the whole pool (sourmash too, once its database is built)
  stage 2     contigs drawn within strata: call (aif flag class; label-conflict reference; flagged by skani or
              Kraken2 but not aif; aif near miss; ambiguous best hit; host only; no species-level hit) x length
              class (< 2 kb, 2-10 kb, >= 10 kb); nested subsets S1 < S2 < S3 take the first n of a seeded random
              order in every cell, so all three share the design
  context     host genus, self reference and host-equivalent references of each genome come from the screen and
              are the same for every tool and for the truth: only contig-level hits differ
  weights     inverse inclusion probabilities of both stages, for catalogue-level estimates
  tools       aif (v1.7 defaults: v1.6 plus the AF fix of bench_real/uhgg_af.py, same hits and ANI; v1.6 results in
              results_uhgg/compare/v16; GTDB R232 0.88 database, top 5), skani 0.4.0 (search --qi), sourmash
              (branchwater manysearch, k = 21, scaled 64 = 1 bit per reference base as aif, containment >= 0.05),
              minimap2 -x asm20 (whole GTDB as a split index), LexicMap (-n 100; S1 only, as its S1 run took 6.3 h,
              scored on S1 with S1's own weights), Kraken2 (database capped at aif's size),
              CAT (DIAMOND against GTDB R232 proteins); all at 48 threads with their database in the page cache
              where it fits in memory
  truth       candidates: every hit of every tool plus minimap2's exhaustive search of all of GTDB R232; each
              (contig, reference) pair aligned with nucmer --maxmatch, delta-filter -q; AF = union of aligned
              contig bases / contig length, ANI = length-weighted identity (the convention of every earlier aif
              benchmark, bench_real/truth.py). Alternatives reported: delta-filter -1 (1-to-1, as dnadiff) and AF
              counted only over alignments of >= 90% identity
  cost        time and peak memory per subset size; build cost per tool; projection to the catalogue by a fit
              of time against query bases, checked by projecting aif and comparing with its measured full run

  python3 -I bench_real/uhgg_compare.py build <mm2|lexicmap|sourmash|cat>
  python3 -I bench_real/uhgg_compare.py pool                       # stage 1 genomes and the contig pool
  python3 -I bench_real/uhgg_compare.py strata                     # stage 2 strata and nested subsets S1 < S2 < S3
  python3 -I bench_real/uhgg_compare.py run <pool|S1|S2|S3> <tool,...>   # timed runs (CPUS=0-47 when the node is quiet)
  python3 -I bench_real/uhgg_compare.py seqmap                     # GTDB sequence -> genome (for minimap2's PAF)
  python3 -I bench_real/uhgg_compare.py truth [S3]                 # MUMmer for every candidate pair
  python3 -I bench_real/uhgg_compare.py evaluate [S3]              # -> results_uhgg/compare/eval_S3.json
  python3 -I bench_real/uhgg_compare.py project                    # -> results_uhgg/compare/cost.json
"""
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
import zlib
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
C = os.path.join(ROOT, "data", "compare")
ENV = os.path.join(ROOT, "data", "envs", "compare", "bin")
REPS = os.path.join(G, "meta", "reps.list")
U = os.path.join(ROOT, "data", "uhgg3")
RESU = os.path.join(ROOT, "results_uhgg")
OUT = os.path.join(RESU, "compare")
AIF = os.path.join(ROOT, "bin", "aif")
AIFDB = os.path.join(G, "aif", "gtdb_0.88.aifd")
SKANI = "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani"  # 0.4.0, as the GTDB sketch database
SKDB = os.path.join(G, "skani_db")
K2 = os.path.join(G, "kraken2_db")
K2BIN = "/usr/bin/kraken2"  # 2.1.3, as the database
MM2I = os.path.join(C, "mm2", "gtdb_r232.asm20.I25G.mmi")
LMI = os.path.join(C, "lexicmap", "gtdb_r232.lmi")
SMZ = os.path.join(C, "sourmash", "gtdb_r232.k21.s64.sig.zip")
SMR = os.path.join(C, "sourmash", "gtdb_r232.k21.s64.rocksdb")
CATDB = os.path.join(C, "cat_gtdb", "db")
NPART, SEED, POOL_PARTS, N_MAG, N_ISO = 128, 20261125, 8, 1000, 300
THREADS = os.environ.get("THREADS", "48")
CPUS = os.environ.get("CPUS", "")  # e.g. 0-47 for the timed runs; empty while builds share the node


def timed(cmd, out=None, err=None, cpus=None, cwd=None, env=None):
    """run cmd under /usr/bin/time; returns (seconds, peak KiB)"""
    tf = tempfile.NamedTemporaryFile(delete=False, dir=C).name
    pre = ["taskset", "-c", cpus] if cpus else []
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + pre + cmd, stdout=fo, stderr=fe, check=True,
                       cwd=cwd, env=env)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def du(path):
    if os.path.isfile(path):
        return os.path.getsize(path)
    return sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(path) for f in fs)


def save_build(tool, rec):
    p = os.path.join(C, f"build_{tool}.json")
    json.dump(rec, open(p, "w"), indent=1)
    print(tool, rec, flush=True)


def ref_paths():
    return [l.split("\t")[0] for l in open(REPS)]


def build_mm2():
    """all GTDB R232 representatives as one gzip stream; minimap2 -x asm20 split index (-I 25G: the build peaks at
    about 8.5 bytes per base, so 25 Gbp parts fit beside other jobs on a 503 GB node)"""
    d = os.path.join(C, "mm2")
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, "gtdb_r232.fa.gz")
    if not os.path.exists(fa):
        with open(fa + ".tmp", "wb") as fo:
            for p in ref_paths():
                with open(p, "rb") as fi:
                    while True:
                        b = fi.read(1 << 24)
                        if not b:
                            break
                        fo.write(b)
        os.replace(fa + ".tmp", fa)
    mmi = os.path.join(d, "gtdb_r232.asm20.I25G.mmi")
    t, m = timed([os.path.join(ENV, "minimap2"), "-x", "asm20", "-I", "25G", "-t", "16", "-d", mmi + ".tmp", fa],
                 err=os.path.join(d, "index.log"), cpus="0-15")
    os.replace(mmi + ".tmp", mmi)
    save_build("minimap2", {"build_s": t, "build_kb": m, "bytes": du(mmi), "threads": 16, "split": "25G"})


def build_lexicmap():
    d = os.path.join(C, "lexicmap")
    os.makedirs(d, exist_ok=True)
    lst = os.path.join(d, "files.txt")
    with open(lst, "w") as fh:
        fh.write("\n".join(ref_paths()) + "\n")
    out = os.path.join(d, "gtdb_r232.lmi")
    t, m = timed([os.path.join(ENV, "lexicmap"), "index", "-S", "-X", lst, "-O", out, "-j", "32", "--force",
                  "-N", r"^(GC[AF]_\d+\.\d+)", "-G", os.path.join(d, "skipped.tsv")],
                 err=os.path.join(d, "index.log"), cpus="16-47")
    save_build("lexicmap", {"build_s": t, "build_kb": m, "bytes": du(out), "threads": 32})


def build_sourmash():
    """k = 21 (as the earlier sourmash baseline), scaled 64: 64-bit hashes at 1 per 64 bases = 1 bit per reference
    base, the budget of aif's 0.88 database (82.5 GB over 654.7 Gbp)"""
    d = os.path.join(C, "sourmash")
    os.makedirs(d, exist_ok=True)
    lst = os.path.join(d, "refs.csv")
    with open(lst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "genome_filename", "protein_filename"])
        for l in open(REPS):
            p, acc = l.rstrip("\n").split("\t")
            w.writerow([acc, p, ""])
    z = os.path.join(d, "gtdb_r232.k21.s64.sig.zip")
    t, m = timed([os.path.join(ENV, "sourmash"), "scripts", "manysketch", lst, "-p", "k=21,scaled=64,dna", "-c", "16",
                  "-o", z], err=os.path.join(d, "sketch.log"), cpus="48-63")
    rec = {"build_s": t, "build_kb": m, "bytes": du(z), "threads": 16}
    rdb = os.path.join(d, "gtdb_r232.k21.s64.rocksdb")
    t2, m2 = timed([os.path.join(ENV, "sourmash"), "scripts", "index", z, "-k", "21", "-s", "64", "-o", rdb, "-c", "16",
                    "--internal-storage"], err=os.path.join(d, "index.log"), cpus="48-63")
    rec.update({"index_s": t2, "index_kb": m2, "index_bytes": du(rdb)})
    save_build("sourmash", rec)


def build_sourmash_index():
    """branchwater 0.9.15 'index' does not read the 84 GB (zip64) sketch zip: it falls back to reading it as a
    path list and fails after 2 h 15 min. The 199,923 sketches are unpacked unchanged (stored, not recompressed)
    and indexed from a path list instead. The manysketch run's own timer was lost with that failure; its wall time
    is taken from file times (refs.csv written just before it started; the zip at its end)."""
    import zipfile
    d = os.path.join(C, "sourmash")
    z = os.path.join(d, "gtdb_r232.k21.s64.sig.zip")
    sd = os.path.join(d, "sigs")
    lst = os.path.join(d, "sigs.txt")
    t0 = time.time()
    if not os.path.exists(lst):
        os.makedirs(sd, exist_ok=True)
        with zipfile.ZipFile(z) as zf, open(lst + ".tmp", "w") as fh:
            for n in zf.namelist():
                if n.endswith(".sig.gz"):
                    out = os.path.join(sd, os.path.basename(n))
                    if os.path.exists(out) and os.path.getsize(out) == zf.getinfo(n).file_size:  # restart
                        fh.write(out + "\n")
                        continue
                    with zf.open(n) as fi, open(out, "wb") as fo:
                        shutil.copyfileobj(fi, fo, 1 << 20)
                    fh.write(out + "\n")
        os.replace(lst + ".tmp", lst)
    ex = time.time() - t0
    rdb = SMR
    t2, m2 = timed([os.path.join(ENV, "sourmash"), "scripts", "index", lst, "-k", "21", "-s", "64", "-o", rdb, "-c", "16",
                    "--internal-storage"], err=os.path.join(d, "index.log"), cpus="48-63")
    sk = os.path.getmtime(z) - os.path.getmtime(os.path.join(d, "refs.csv"))
    save_build("sourmash", {"build_s": round(sk, 2), "build_kb": None, "build_from_mtimes": True, "bytes": du(z),
                            "threads": 16, "unpack_s": round(ex, 2), "index_s": t2, "index_kb": m2,
                            "index_bytes": du(rdb), "failed_zip_index_s": 8086})


CAT_LOCAL = """
import os, sys
import CAT_pack.io.download as D
local = sys.argv.pop(1)
orig = D.download_singleton
def ds(url, path, report=lambda *_: None):
    if path.name.endswith("gtdb_proteins_aa_reps.tar.gz"):
        if not path.exists():
            os.symlink(local, path)
        n = path.stat().st_size
        report(path, n, n)
        return
    orig(url, path, report)
D.download_singleton = ds
from CAT_pack.cli import app
sys.argv[0] = "CAT_pack"
sys.exit(app())
"""


def cat_download():
    """'CAT_pack download --db gtdb' with the 132 GB protein archive taken from a resumed local download (CAT's own
    download has no resume); CAT's checksum test and processing run unchanged"""
    d = os.path.join(C, "cat_gtdb")
    local = os.path.join(d, "proteins_dl", "gtdb_proteins_aa_reps_r232.tar.gz")
    env = dict(os.environ, PATH=ENV + ":" + os.environ["PATH"])
    t, m = timed([os.path.join(ENV, "python3"), "-I", "-c", CAT_LOCAL, local, "download", "--db", "gtdb", "--out_dir",
                  os.path.join(d, "dl"), "-n", "32"], err=os.path.join(d, "download_local.log"), cpus="16-47", env=env)
    json.dump({"process_s": t, "process_kb": m}, open(os.path.join(d, "download_local.json"), "w"))


def build_cat():
    """CAT_pack prepare on the GTDB R232 protein representatives fetched by 'CAT_pack download --db gtdb'"""
    d = os.path.join(C, "cat_gtdb")
    dl = os.path.join(d, "dl")
    if not os.path.isdir(dl) or not any(f.endswith("accession2taxid.txt.gz") for f in os.listdir(dl)):
        cat_download()
    pick = lambda suf: sorted(os.path.join(dl, f) for f in os.listdir(dl) if f.endswith(suf))[-1]
    db = os.path.join(d, "db")
    env = dict(os.environ, PATH=ENV + ":" + os.environ["PATH"])
    t, m = timed([os.path.join(ENV, "CAT_pack"), "prepare", "--db_fasta", pick(".fa.gz"), "--names", pick("names.dmp"),
                  "--nodes", pick("nodes.dmp"), "--acc2tax", pick("accession2taxid.txt.gz"), "-d", db, "-n", "32",
                  "--common_prefix", "gtdb_r232"], err=os.path.join(d, "prepare.log"), cpus="16-47", env=env)
    save_build("cat", {"build_s": t, "build_kb": m, "bytes": du(db), "threads": 32})


def part_of(g):
    return zlib.crc32(g.encode()) % NPART


def extract(job):
    """contigs of the chosen genomes from one screen part -> (fasta chunk path, contig rows)"""
    p, gs, d = job
    rows = []
    for line in open(os.path.join(U, "fa", f"part{p:03d}.contigs.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[1] in gs:
            rows.append(f[:4])
    keep = {r[0] for r in rows}
    o = os.path.join(d, f"part{p:03d}.fa")
    with open(os.path.join(U, "fa", f"part{p:03d}.fa")) as fi, open(o, "w") as fo:
        w = False
        for line in fi:
            if line[0] == ">":
                w = line[1:].split()[0] in keep
            if w:
                fo.write(line)
    return o, rows


def sample_pool():
    """stage 1: POOL_PARTS screen parts at random (parts are a hash partition of the genomes), then N_MAG MAGs and
    N_ISO isolates at random among their held-out genomes; every contig of a drawn genome joins the pool"""
    d = os.path.join(C, "pool")
    os.makedirs(d, exist_ok=True)
    pilot = set(open(os.path.join(U, "pilot", "genomes.txt")).read().split())
    meta = {}
    with open(os.path.join(U, "meta", "genomes-all_metadata.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["Genome"]] = r
    scr = {}
    with gzip.open(os.path.join(RESU, "genomes.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            scr[r["genome"]] = r
    rng = random.Random(SEED)
    parts = sorted(rng.sample(range(NPART), POOL_PARTS))
    elig = sorted(g for g in scr if g not in pilot and part_of(g) in parts)
    pick, pi = [], {}
    for t, n in (("MAG", N_MAG), ("Isolate", N_ISO)):
        gs = [g for g in elig if meta[g]["Genome_type"] == t]
        ch = sorted(rng.sample(gs, min(n, len(gs))))
        for g in ch:
            pi[g] = POOL_PARTS / NPART * len(ch) / len(gs)
        pick += ch
        print(f"stage 1: {len(ch)} of {len(gs)} held-out {t}s in parts {parts}", flush=True)
    with open(os.path.join(d, "genomes.tsv"), "w") as fh:
        fh.write("genome\ttype\tpi1\tcompleteness\tcontamination\thost_genus\thost_src\tself_ref\tequiv_xgenus_refs\n")
        for g in pick:
            r, m = scr[g], meta[g]
            fh.write(f"{g}\t{m['Genome_type']}\t{pi[g]:.6g}\t{m['Completeness']}\t{m['Contamination']}\t{r['host_genus']}\t"
                     f"{r['host_src']}\t{r['self_ref']}\t{r['equiv_xgenus_refs']}\n")
    gs = set(pick)
    with Pool(len(parts)) as pool:
        outs = pool.map(extract, [(p, {g for g in gs if part_of(g) == p}, d) for p in parts])
    with open(os.path.join(d, "pool.fa"), "w") as fo, open(os.path.join(d, "contigs.tsv"), "w") as ft:
        ft.write("contig\tgenome\tlength\tgc\n")
        for o, rows in outs:
            for r in rows:
                ft.write("\t".join(r) + "\n")
            with open(o) as fi:
                for line in fi:
                    fo.write(line)
            os.remove(o)
    n = sum(len(r) for _, r in outs)
    bp = sum(int(x[2]) for _, r in outs for x in r)
    print(f"pool: {len(pick)} genomes, {n} contigs, {bp / 1e9:.2f} Gbp", flush=True)


def warm(paths):
    for p in paths:
        fs = [os.path.join(dd, f) for dd, _, ff in os.walk(p) for f in ff] if os.path.isdir(p) else [p]
        for f in fs:
            subprocess.run(["cat", f], stdout=subprocess.DEVNULL, check=True)


def run_tool(tool, fa, d, warm_db=True):
    """one timed run of a tool on a FASTA of contigs; outputs in d/<tool>.*; returns {s, kb, ...}"""
    os.makedirs(d, exist_ok=True)
    o = os.path.join(d, tool)
    cp = CPUS or None
    env = dict(os.environ, PATH=ENV + ":" + os.environ["PATH"])
    if tool == "aif":
        warm_db and warm([AIFDB])
        t, m = timed([AIF, "db", "search", AIFDB, fa, "-n", "5", "-t", THREADS], out=o + ".tsv", err=o + ".log", cpus=cp)
    elif tool == "skani":
        warm_db and warm([SKDB])
        t, m = timed([SKANI, "search", "-d", SKDB, "--qi", fa, "-t", THREADS, "-o", o + ".tsv"], err=o + ".log", cpus=cp)
    elif tool == "kraken2":
        warm_db and warm([os.path.join(K2, f) for f in ("hash.k2d", "opts.k2d", "taxo.k2d")])
        t, m = timed([K2BIN, "--db", K2, "--threads", THREADS, "--output", o + ".tsv", "--report", o + ".report", fa],
                     err=o + ".log", cpus=cp)
    elif tool == "sourmash":
        warm_db and warm([SMR])
        qz = o + ".query.sig.zip"
        lst = o + ".query.csv"
        with open(lst, "w") as fh:
            fh.write(f"name,genome_filename,protein_filename\nq,{fa},\n")
        t1, m1 = timed([os.path.join(ENV, "sourmash"), "scripts", "manysketch", lst, "-s", "-p", "k=21,scaled=64,dna",
                        "-c", THREADS, "-o", qz], err=o + ".sketch.log", cpus=cp, env=env)
        t2, m2 = timed([os.path.join(ENV, "sourmash"), "scripts", "manysearch", qz, SMR, "-k", "21", "-s", "64", "-t", "0.05",
                        "-c", THREADS, "-o", o + ".csv", "--ignore-abundance", "-N"], err=o + ".log", cpus=cp, env=env)
        return {"s": t1 + t2, "kb": max(m1, m2), "sketch_s": t1, "search_s": t2}
    elif tool == "minimap2":
        sp = os.path.join(d, "mm2split")
        t, m = timed([os.path.join(ENV, "minimap2"), "-c", "-x", "asm20", "-N", "50", "-t", THREADS, "--split-prefix", sp,
                      MM2I, fa], out=o + ".paf", err=o + ".log", cpus=cp)
    elif tool == "lexicmap":
        # seeds (344 GB) and metadata in the page cache, as aif's database is; the 161 GB of genome sequences
        # (read only to align candidates) do not fit beside them in 503 GB of RAM and stay on disk
        # -n 100: keep the top 100 genomes per contig in chaining (LexicMap's advice: >= 100); only the top 20 are
        # scored, and the default (all genomes) aligns every one of thousands of seed-matched genomes per contig
        warm_db and warm([os.path.join(LMI, f) for f in os.listdir(LMI) if f != "genomes"])
        t, m = timed([os.path.join(ENV, "lexicmap"), "search", "-d", LMI, fa, "-o", o + ".tsv", "-j", THREADS, "-n", "100"],
                     err=o + ".log", cpus=cp, env=env)
    elif tool == "cat":
        # CAT_pack 7 keeps the DIAMOND database and the taxonomy in one directory
        warm_db and warm([CATDB])
        t, m = timed([os.path.join(ENV, "CAT_pack"), "contigs", "-c", fa, "-d", CATDB, "-n", THREADS, "-o", o],
                     err=o + ".log", cpus=cp, env=env)   # CAT_pack 7 contigs has no --force (add_names does)
        subprocess.run([os.path.join(ENV, "CAT_pack"), "add_names", "-i", o + ".contig2classification.txt", "-o",
                        o + ".names.txt", "-d", CATDB, "--only_official", "--force"], env=env, check=True,
                       capture_output=True)
    else:
        raise ValueError(tool)
    return {"s": t, "kb": m}


NOTAX = ("?", "?", "?", "?")
ANI_SP, AF_MIN = 0.95, 0.5
CATS = ["anchored", "accessory", "mobile", "noncoding", "label", "skani_only", "kraken2_only", "near_miss", "ambiguous",
        "host", "weak"]
LENS = [(500, 2000, "<2kb"), (2000, 10000, "2-10kb"), (10000, 10 ** 12, ">=10kb")]
SIZES = {"S1": 1000, "S2": 3000, "S3": 10000}


def gtdb_tax():
    """{accession: (phylum, family, genus, species)} for GTDB R232 representatives"""
    T = {}
    for fn in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            rd = csv.reader(fh, delimiter="\t")
            h = next(rd)
            ia, it, ir = h.index("accession"), h.index("gtdb_taxonomy"), h.index("gtdb_representative")
            for r in rd:
                if r[ir] == "t":
                    x = [y[3:] for y in r[it].split(";")]
                    T[r[ia][3:]] = (x[1], x[4], x[5], x[6])
    return T


def acc_of(path):
    b = os.path.basename(path)
    return b[: b.index("_genomic")] if "_genomic" in b else b.split(".fna")[0]


def load_aif(path):
    H = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            H.setdefault(f[0], [])
            if f[1] != "*" and f[9] == "ok":
                H[f[0]].append((f[1], float(f[3]), float(f[4]), float(f[13]), float(f[14]), float(f[15]), float(f[16])))
    return H


def load_skani(path):
    H = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            H.setdefault(p[6].split()[0], []).append((acc_of(p[0]), float(p[2]) / 100, float(p[4]) / 100))
    for c in H:  # skani's ranking (ANI); ties by AF, then name
        H[c].sort(key=lambda x: (-x[1], -x[2], x[0]))
    return H


def k2_genus(path):
    """{contig: Kraken2 genus or None} (GTDB taxonomy of the Kraken2 database)"""
    nodes, names = {}, {}
    for line in open(os.path.join(K2, "taxonomy", "nodes.dmp")):
        f = [x.strip() for x in line.split("|")]
        nodes[int(f[0])] = (int(f[1]), f[2])
    for line in open(os.path.join(K2, "taxonomy", "names.dmp")):
        f = [x.strip() for x in line.split("|")]
        names[int(f[0])] = f[1]
    memo = {}

    def genus(t):
        if t not in memo:
            g, x = None, t
            while x != 1 and x in nodes:
                if nodes[x][1] == "genus":
                    g = names[x][3:]
                    break
                x = nodes[x][0]
            memo[t] = g
        return memo[t]
    out = {}
    for line in open(path):
        f = line.split("\t", 3)
        out[f[1]] = genus(int(f[2])) if f[0] == "C" else None
    return out


def contexts():
    """{genome: (host genus, self reference, cross-genus host-equivalent references, type, pi1)} of the pool"""
    X = {}
    for r in csv.DictReader(open(os.path.join(C, "pool", "genomes.tsv")), delimiter="\t"):
        X[r["genome"]] = (r["host_genus"] or None, r["self_ref"] or None,
                          set(x for x in r["equiv_xgenus_refs"].split(";") if x), r["type"], float(r["pi1"]))
    return X


def foreign(hits, ctx, TX):
    """the screen's foreign-contig rule on a tool's ranked hits [(ref, ani, af or None, ...)], top 5"""
    host, sref, E = ctx[:3]
    hs = [h for h in hits if h[0] != sref][:5]
    if not hs or not host:
        return False
    t = hs[0]
    if not (t[1] >= ANI_SP and (t[2] is None or t[2] >= AF_MIN) and TX.get(t[0], NOTAX)[2] != host):
        return False
    return not any(h[1] >= ANI_SP and (h[0] in E or TX.get(h[0], NOTAX)[2] == host) for h in hs)


def category(c, flag, aif_h, sk_h, k2g, ctx, TX, label_refs):
    if flag:
        return "anchored" if flag["anchored"] == "1" else ("mobile" if flag["cls"] == "mobile" else
                                                           "noncoding" if flag["cls"] == "noncoding" else "accessory")
    host, sref, E = ctx[:3]
    hs = [h for h in aif_h if h[0] != sref][:5]
    if hs and hs[0][0] in label_refs:
        return "label"
    if foreign(sk_h, ctx, TX):
        return "skani_only"
    if host and k2g and k2g != host and k2g not in {TX.get(r, NOTAX)[2] for r in E | ({sref} if sref else set())}:
        return "kraken2_only"
    if not hs:
        return "host" if aif_h else "weak"
    t = hs[0]
    tg = TX.get(t[0], NOTAX)[2]
    hosty = lambda h: h[0] in E or TX.get(h[0], NOTAX)[2] == host
    if host and tg != host and t[0] not in E and not any(h[1] >= ANI_SP and hosty(h) for h in hs):
        if (0.93 <= t[1] < ANI_SP and t[2] >= 0.3) or (t[1] >= ANI_SP and 0.3 <= t[2] < AF_MIN):
            return "near_miss"
    if len(hs) > 1 and hs[0][1] >= 0.90 and TX.get(hs[1][0], NOTAX)[2] != tg and hs[0][1] - hs[1][1] <= 0.005:
        return "ambiguous"
    if t[1] >= ANI_SP and hosty(t):
        return "host"
    return "weak"


def waterfill(avail, target):
    """per-cell sizes: min(avail, L), with L the smallest level whose total reaches the target"""
    lo, hi = 0, max(avail.values())
    while lo < hi:
        mid = (lo + hi) // 2
        if sum(min(a, mid) for a in avail.values()) >= target:
            hi = mid
        else:
            lo = mid + 1
    return {k: min(a, lo) for k, a in avail.items()}


def stage2_cells():
    """stage-2 cells {(category, length class): pool contigs in one seeded random order}, rows {contig: (genome,
    length, category, length class)} and the genome contexts; subsets take prefixes of every cell's order"""
    TX = gtdb_tax()
    X = contexts()
    label_refs = {r["reference"] for r in csv.DictReader(open(os.path.join(RESU, "label_inconsistent_refs.tsv")),
                                                          delimiter="\t")}
    pd_ = os.path.join(C, "runs", "pool")
    A = load_aif(os.path.join(pd_, "aif.tsv"))
    S = load_skani(os.path.join(pd_, "skani.tsv"))
    K = k2_genus(os.path.join(pd_, "kraken2.tsv"))
    F = {}
    for p in sorted({part_of(g) for g in X}):
        for r in csv.DictReader(open(os.path.join(U, "screen", f"part{p:03d}.flagged.tsv")), delimiter="\t"):
            if r["genome"] in X:
                F[r["contig"]] = r
    cells = {}
    rows = {}
    for r in csv.DictReader(open(os.path.join(C, "pool", "contigs.tsv")), delimiter="\t"):
        c, g, ln = r["contig"], r["genome"], int(r["length"])
        if ln < 500:  # listed in the screen's contig table, never searched
            continue
        cat = category(c, F.get(c), A.get(c, []), S.get(c, []), K.get(c), X[g], TX, label_refs)
        lc = next(n for a, b, n in LENS if a <= ln < b)
        cells.setdefault((cat, lc), []).append(c)
        rows[c] = (g, ln, cat, lc)
    rng = random.Random(SEED + 1)
    for k in sorted(cells):
        cells[k].sort()
        rng.shuffle(cells[k])
    return cells, rows, X


def strata():
    """stage 2: categories x length classes; nested subsets S1 < S2 < S3 from one seeded order per cell"""
    cells, rows, X = stage2_cells()
    avail = {k: len(v) for k, v in cells.items()}
    alloc = {s: waterfill(avail, n) for s, n in SIZES.items()}
    d = os.path.join(C, "sample")
    os.makedirs(d, exist_ok=True)
    member = {}
    for s in SIZES:
        for k, v in cells.items():
            for c in v[: alloc[s][k]]:
                member.setdefault(c, set()).add(s)
    with open(os.path.join(d, "sample.tsv"), "w") as fh:
        fh.write("contig\tgenome\ttype\tlength\tcategory\tlength_class\tN_cell\tn_cell_S3\tpi1\tpi2_S3\tS1\tS2\tS3\n")
        for c in sorted(member):
            g, ln, cat, lc = rows[c]
            k = (cat, lc)
            fh.write(f"{c}\t{g}\t{X[g][3]}\t{ln}\t{cat}\t{lc}\t{avail[k]}\t{alloc['S3'][k]}\t{X[g][4]:.6g}\t"
                     f"{alloc['S3'][k] / avail[k]:.6g}\t" + "\t".join(str(int(s in member[c])) for s in SIZES) + "\n")
    outs = {s: open(os.path.join(d, f"{s}.fa"), "w") for s in SIZES}
    with open(os.path.join(C, "pool", "pool.fa")) as fi:
        w = None
        for line in fi:
            if line[0] == ">":
                w = member.get(line[1:].split()[0])
            if w:
                for s in w:
                    outs[s].write(line)
    for f in outs.values():
        f.close()
    summ = {"cells": {f"{a}|{b}": {"N": avail[(a, b)], **{s: alloc[s][(a, b)] for s in SIZES}} for a, b in sorted(cells)}}
    for s in SIZES:
        summ[s] = {"contigs": sum(alloc[s].values()),
                   "bp": sum(rows[c][1] for c in member if s in member[c])}
    json.dump(summ, open(os.path.join(d, "design.json"), "w"), indent=1)
    for k in sorted(cells):
        print(f"{k[0]:13s} {k[1]:7s} N={avail[k]:7d} " + " ".join(f"{s}={alloc[s][k]}" for s in SIZES))
    print({s: summ[s] for s in SIZES}, flush=True)


def _headers(path):
    acc = acc_of(path)
    out = subprocess.run(f"zcat {path} | grep '>'", shell=True, capture_output=True, text=True).stdout
    return [(l[1:].split()[0], acc) for l in out.splitlines()]


def seqmap():
    """GTDB sequence ID -> genome accession, for minimap2's PAF. From the Kraken2 build's seqid2taxid.map (same
    first header tokens, same reps.list): its taxids are GTDB species nodes, one per species representative."""
    p = os.path.join(C, "mm2", "seq2genome.tsv")
    if os.path.exists(p):
        return {a: b for a, b in (l.rstrip("\n").split("\t") for l in open(p))}
    K = os.path.join(G, "kraken2_db")
    name = {}
    for line in open(os.path.join(K, "taxonomy", "names.dmp")):
        f = line.split("\t|\t")
        if f[1].startswith("s__"):
            name[f[0]] = f[1]
    sp2acc = {}
    TX = {}
    for fn in ("bac120_taxonomy.tsv.gz", "ar53_taxonomy.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            for line in fh:
                acc, tax = line.rstrip("\n").split("\t")
                TX[acc[3:]] = tax.split(";")[-1]
    for line in open(REPS):
        acc = line.rstrip("\n").split("\t")[1]
        assert TX[acc] not in sp2acc
        sp2acc[TX[acc]] = acc
    tid2acc = {t: sp2acc[n] for t, n in name.items()}
    assert len(tid2acc) == len(sp2acc), (len(tid2acc), len(sp2acc))
    M = {}
    with open(os.path.join(K, "seqid2taxid.map")) as fi, open(p + ".tmp", "w") as fo:
        for line in fi:
            sid, tid = line.rstrip("\n").split("\t")
            M[sid] = tid2acc[tid]
            fo.write(f"{sid}\t{M[sid]}\n")
    os.replace(p + ".tmp", p)
    return M


def union(iv):
    tot, ce = 0, -1
    for a, b in sorted(iv):
        if b <= ce:
            continue
        tot += b - max(a, ce)
        ce = b
    return tot


def load_minimap2(path, S2G):
    """{contig: [(genome, ani, af)]}: per genome, alignments kept greedily by matches unless half of their contig span
    is already covered (as delta-filter -q); AF = covered / length, ANI = 1 - de weighted by contig span; ranked by
    matching bases"""
    A = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            if de is None:
                continue
            g = S2G.get(f[5], f[5])
            A.setdefault(f[0], {}).setdefault(g, []).append((int(f[9]), int(f[2]), int(f[3]), de, int(f[1])))
    H = {}
    for c, byg in A.items():
        out = []
        for g, al in byg.items():
            al.sort(reverse=True)
            kept, iv = [], []
            for m, a, b, de, ql in al:
                if union(iv + [(a, b)]) - union(iv) >= 0.5 * (b - a):
                    kept.append((a, b, de, m))
                    iv.append((a, b))
            span = sum(b - a for a, b, _, _ in kept)
            ani = sum((1 - de) * (b - a) for a, b, de, _ in kept) / span
            out.append((g, ani, min(1.0, union([(a, b) for a, b, _, _ in kept]) / al[0][4]), sum(m for *_, m in kept)))
        out.sort(key=lambda x: (-x[3], x[0]))
        H[c] = [x[:3] for x in out]
    return H


def load_lexicmap(path):
    """{contig: [(genome, ani, af)]}: AF = qcovGnm, ANI = pident weighted by aligned length; LexicMap's genome order"""
    H, W = {}, {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c, g = f[0], f[3]
            if c not in H:
                H[c] = []
            if (c, g) not in W:
                W[(c, g)] = [0.0, 0, float(f[5]) / 100]
                H[c].append(g)
            W[(c, g)][0] += float(f[10]) / 100 * int(f[9])
            W[(c, g)][1] += int(f[9])
    return {c: [(g, W[(c, g)][0] / W[(c, g)][1], min(1.0, W[(c, g)][2])) for g in gs] for c, gs in H.items()}


def load_sourmash(path):
    """{contig: [(genome, containment ANI, None)]} ranked by containment (sourmash estimates no AF)"""
    H = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ani = r.get("query_containment_ani") or ""
            H.setdefault(r["query_name"].split()[0], []).append(
                (r["match_name"].split()[0], float(ani) if ani else 0.0, None, float(r["containment"])))
    for c in H:
        H[c].sort(key=lambda x: (-x[3], x[0]))
        H[c] = [x[:3] for x in H[c]]
    return H


def load_cat(path):
    """{contig: CAT genus or None} from CAT_pack 7 add_names --only_official: the genus column holds "Name: score",
    "Name*: score" (suggestive, kept) or "no support"; GTDB names without the g__ prefix, as gtdb_tax()"""
    out, ig = {}, None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                ig = f.index("genus")
                continue
            x = f[ig] if len(f) > ig else ""
            out[f[0]] = None if not x or x == "no support" else x.rsplit(":", 1)[0].rstrip("*").strip()
    return out


TOPK = {"aif": 5, "skani": 10, "sourmash": 10, "minimap2": 20, "lexicmap": 20}


HITF = {"aif": "aif.tsv", "skani": "skani.tsv", "sourmash": "sourmash.csv", "minimap2": "minimap2.paf",
        "lexicmap": "lexicmap.tsv"}


def hit_tools(stage):
    """the TOPK tools that have hits for a stage (LexicMap ran on S1 only)"""
    return [t for t in TOPK if os.path.exists(os.path.join(C, "runs", stage, HITF[t]))]


def load_hits(stage, tool, S2G=None):
    f = os.path.join(C, "runs", stage, HITF[tool])
    if tool == "aif":
        return load_aif(f)
    if tool == "skani":
        return load_skani(f)
    if tool == "sourmash":
        return load_sourmash(f)
    if tool == "minimap2":
        return load_minimap2(f, S2G if S2G is not None else seqmap())
    if tool == "lexicmap":
        return load_lexicmap(f)
    raise ValueError(tool)


def truth_job(job):
    """nucmer --maxmatch of a set of contigs against one reference; three conventions"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import truth as T
    T.MUMMER = os.path.join(ROOT, "data", "envs", "mummer", "bin")
    ref, rpath, qfa, L = job
    rows = []
    with tempfile.TemporaryDirectory(dir=os.path.join(C, "tmp")) as td:
        rfa = os.path.join(td, "ref.fa")
        with gzip.open(rpath, "rb") as fi, open(rfa, "wb") as fo:
            fo.write(fi.read())
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(T.MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre, rfa, qfa], check=True,
                       capture_output=True)
        conv = {"q": T.coords(pre + ".delta", "-q"), "one": T.coords(pre + ".delta", "-1")}
    for name, ln in L.items():
        v = []
        for key, rows_ in (("q", conv["q"]), ("one", conv["one"]), ("q90", [r for r in conv["q"] if r[4] >= 90])):
            iv = [(s - 1, e) for n, s, e, lq, idy, _ in rows_ if n == name]
            w = [(lq, idy) for n, s, e, lq, idy, _ in rows_ if n == name]
            af = min(1.0, union(iv) / ln)
            ani = sum(a * b for a, b in w) / sum(a for a, _ in w) / 100 if w else float("nan")
            v += [af, ani]
        rows.append((name, ref, ln, *v))
    os.unlink(qfa)
    return rows


def truth(stage="S3"):
    """candidates: top hits of every tool (TOPK) on the stage; MUMmer for each (contig, reference) pair"""
    src = os.environ.get("HITS_FROM", stage)  # e.g. pool: the same per-contig hits, restricted to the stage below
    tools = os.environ.get("TOOLS", ",".join(hit_tools(src))).split(",")
    S2G = seqmap() if "minimap2" in tools else None
    cand = {}
    for tool in tools:
        k = TOPK[tool]
        H = load_hits(src, tool, S2G)
        for c, hs in H.items():
            for h in hs[:k]:
                cand.setdefault(h[0], set()).add(c)
    paths = {acc: p for p, acc in (l.rstrip("\n").split("\t") for l in open(REPS))}
    seqs, name = {}, None
    for line in open(os.path.join(C, "sample", f"{stage}.fa")):
        if line[0] == ">":
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    d = os.path.join(C, "truth")
    os.makedirs(os.path.join(C, "tmp"), exist_ok=True)
    os.makedirs(d, exist_ok=True)
    tp = os.path.join(d, f"{stage}.pairs.tsv")
    done = set()
    if os.path.exists(tp):
        with open(tp) as fh:
            next(fh)
            done = {tuple(l.split("\t", 2)[:2]) for l in fh}
    jobs = []
    for ref, cs in sorted(cand.items(), key=lambda x: -len(x[1])):
        todo = sorted(c for c in cs if (c, ref) not in done and c in seqs)
        if not todo or ref not in paths:
            continue
        qfa = os.path.join(C, "tmp", f"tq_{ref}.fa")
        with open(qfa, "w") as fh:
            for c in todo:
                fh.write(f">{c}\n{seqs[c]}\n")
        jobs.append((ref, paths[ref], qfa, {c: len(seqs[c]) for c in todo}))
    print(f"truth {stage}: {sum(len(v) for v in cand.values())} pairs, {len(cand)} references, {len(jobs)} jobs", flush=True)
    new = not os.path.exists(tp)
    with open(tp, "a") as fh, Pool(int(os.environ.get("JOBS", "56"))) as pool:
        if new:
            fh.write("contig\treference\tlength\taf_q\tani_q\taf_one\tani_one\taf_q90\tani_q90\n")
        for i, rows in enumerate(pool.imap_unordered(truth_job, jobs, chunksize=1)):
            for r in rows:
                fh.write("\t".join([r[0], r[1], str(r[2])] + [f"{x:.5f}" for x in r[3:]]) + "\n")
            fh.flush()
            if (i + 1) % 500 == 0:
                print(f"  {i + 1}/{len(jobs)}", flush=True)


def load_truth(stage, conv="q"):
    """{contig: [(ref, ani, af)]} ranked by identical aligned bases (ANI x AF), aligned pairs only"""
    T = {}
    ia, ib = {"q": (3, 4), "one": (5, 6), "q90": (7, 8)}[conv]
    p = os.path.join(C, "truth", f"{stage}.pairs.tsv")
    if not os.path.exists(p):  # S1 < S2 < S3: a subset's pairs are S3's, restricted by the caller to its contigs
        p = os.path.join(C, "truth", "S3.pairs.tsv")
    with open(p) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            af, ani = float(f[ia]), float(f[ib])
            if af > 0 and ani == ani:
                T.setdefault(f[0], []).append((f[1], ani, af))
    for c in T:
        T[c].sort(key=lambda x: (-x[1] * x[2], x[0]))
    return T


def wmean(xs, ws):
    sw = sum(ws)
    return sum(x * w for x, w in zip(xs, ws)) / sw if sw else float("nan")


def evaluate(stage="S3"):
    """estimation, search and decision accuracy against MUMmer; unweighted (sample) and weighted (catalogue)"""
    TX = gtdb_tax()
    X = contexts()
    S = {r["contig"]: r for r in csv.DictReader(open(os.path.join(C, "sample", "sample.tsv")), delimiter="\t")
         if r[stage] == "1"}
    # each subset is a simple random sample within its cell (a prefix of one seeded order): pi2 = n_cell / N_cell
    ncell = {}
    for r in S.values():
        k = (r["category"], r["length_class"])
        ncell[k] = ncell.get(k, 0) + 1
    wt = {c: 1 / (float(r["pi1"]) * ncell[(r["category"], r["length_class"])] / int(r["N_cell"])) for c, r in S.items()}
    src = os.environ.get("HITS_FROM", stage)  # pool holds the same per-contig hits for every S3 contig
    avail = hit_tools(src)
    tools = [t for t in os.environ.get("TOOLS", ",".join(avail)).split(",") if t in TOPK]
    S2G = seqmap() if "minimap2" in tools else None
    HT = {t: load_hits(src, t, S2G) for t in tools}
    d = os.path.join(C, "runs", src)
    GEN = {"kraken2": k2_genus(os.path.join(d, "kraken2.tsv"))}
    if os.path.exists(os.path.join(d, "cat.names.txt")):
        GEN["cat"] = load_cat(os.path.join(d, "cat.names.txt"))
    res = {"stage": stage, "contigs": len(S), "hits_from": src, "tools": tools}
    for conv in ("q", "one", "q90"):
        TR = load_truth(stage, conv)
        R = res[conv] = {}
        tcall = {c: foreign(TR.get(c, []), X[S[c]["genome"]], TX) for c in S}
        sp = {c for c in S if TR.get(c) and TR[c][0][1] >= ANI_SP and TR[c][0][2] >= AF_MIN}
        R["truth"] = {"foreign": sum(tcall.values()), "species_hit": len(sp),
                      "foreign_w": wmean([float(tcall[c]) for c in S], [wt[c] for c in S])}
        for tool in tools + list(GEN):
            out = {}
            if tool in HT:
                H = HT[tool]
                # (1) estimation on the tool's top hit (and on all its top-k hits), truth ANI >= 0.90
                for lab, k in (("top1", 1), ("topk", TOPK[tool])):
                    ea, ef, ws, lens, cov_a, cov_f = [], [], [], [], [], []
                    for c in S:
                        tr = {r: (a, f) for r, a, f in TR.get(c, [])}
                        sref = X[S[c]["genome"]][1]
                        for h in [h for h in H.get(c, []) if h[0] != sref][:k]:
                            if h[0] in tr and tr[h[0]][0] >= 0.90:
                                ea.append(h[1] - tr[h[0]][0])
                                ef.append(None if h[2] is None else h[2] - tr[h[0]][1])
                                ws.append(wt[c])
                                lens.append(S[c]["length_class"])
                                if tool == "aif":
                                    cov_a.append(h[3] <= tr[h[0]][0] <= h[4])
                                    cov_f.append(h[5] <= tr[h[0]][1] <= h[6])
                    e = {"pairs": len(ea), "ani_mae_pp": 100 * sum(map(abs, ea)) / max(1, len(ea)),
                         "ani_bias_pp": 100 * sum(ea) / max(1, len(ea)),
                         "ani_mae_pp_w": 100 * wmean([abs(x) for x in ea], ws),
                         "ani_within1pp": sum(abs(x) <= 0.01 for x in ea) / max(1, len(ea))}
                    fa = [x for x in ef if x is not None]
                    if fa:
                        e.update({"af_mae_pp": 100 * sum(map(abs, fa)) / len(fa), "af_bias_pp": 100 * sum(fa) / len(fa),
                                  "af_mae_pp_w": 100 * wmean([abs(x) for x in ef if x is not None],
                                                             [w for x, w in zip(ef, ws) if x is not None])})
                    by = {}
                    for _, _, ln in LENS:
                        sel = [i for i, l in enumerate(lens) if l == ln]
                        if sel:
                            by[ln] = {"pairs": len(sel), "ani_mae_pp": 100 * sum(abs(ea[i]) for i in sel) / len(sel)}
                            f2 = [ef[i] for i in sel if ef[i] is not None]
                            if f2:
                                by[ln]["af_mae_pp"] = 100 * sum(map(abs, f2)) / len(f2)
                    e["by_length"] = by
                    if cov_a:
                        e["ci_cover_ani"], e["ci_cover_af"] = sum(cov_a) / len(cov_a), sum(cov_f) / len(cov_f)
                    out[lab] = e
                # (2) search: contigs with a species-level hit (truth): reported, same top hit, same genus
                rep = top = gen = 0
                for c in sp:
                    sref = X[S[c]["genome"]][1]
                    hs = [h for h in H.get(c, []) if h[0] != sref]
                    tb = [r for r in TR[c] if r[0] != sref]
                    if not tb:
                        continue
                    if hs:
                        rep += 1
                        best = tb[0][1] * tb[0][2]
                        tv = {r: a * f for r, a, f in tb}
                        top += tv.get(hs[0][0], 0) >= 0.98 * best
                        gen += TX.get(hs[0][0], NOTAX)[2] == TX.get(tb[0][0], NOTAX)[2]
                out["search"] = {"n": len(sp), "reported": rep / max(1, len(sp)), "top_hit": top / max(1, len(sp)),
                                 "top_genus": gen / max(1, len(sp))}
                call = {c: foreign(H.get(c, []), X[S[c]["genome"]], TX) for c in S}
            else:
                G2 = GEN[tool]

                def gcall(c):
                    host, sref, E = X[S[c]["genome"]][:3]
                    g = G2.get(c)
                    return bool(host and g and g != host and g not in {TX.get(r, NOTAX)[2] for r in E | ({sref} if sref else set())})
                call = {c: gcall(c) for c in S}
            # (3) decision: the screen's foreign call against the truth's
            tp = sum(call[c] and tcall[c] for c in S)
            fp = sum(call[c] and not tcall[c] for c in S)
            fn = sum(not call[c] and tcall[c] for c in S)
            tn = len(S) - tp - fp - fn
            W = lambda f: sum(wt[c] for c in S if f(c))
            wtp, wfp, wfn = W(lambda c: call[c] and tcall[c]), W(lambda c: call[c] and not tcall[c]), W(lambda c: not call[c] and tcall[c])
            out["decision"] = {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
                               "precision": tp / max(1, tp + fp), "recall": tp / max(1, tp + fn),
                               "f1": 2 * tp / max(1, 2 * tp + fp + fn),
                               "precision_w": wtp / max(1e-12, wtp + wfp), "recall_w": wtp / max(1e-12, wtp + wfn),
                               "flag_rate_w": wmean([float(call[c]) for c in S], [wt[c] for c in S])}
            byc = {}
            for cat in CATS:
                cs = [c for c in S if S[c]["category"] == cat]
                if cs:
                    byc[cat] = {"n": len(cs), "agree": sum(call[c] == tcall[c] for c in cs) / len(cs),
                                "flag": sum(call[c] for c in cs) / len(cs), "truth_flag": sum(tcall[c] for c in cs) / len(cs)}
            out["decision"]["by_category"] = byc
            R[tool] = out
        # (4) paired: aif against each tool on the (contig, reference) pairs both report, truth ANI >= 0.90, so the
        # two error rates come from the same pairs (the per-tool rows above use each tool's own reported pairs)
        if "aif" in HT:
            P = R["paired"] = {}
            for tool in [t for t in tools if t != "aif"]:
                rows = []
                for c in S:
                    sref = X[S[c]["genome"]][1]
                    tr = {r: (a, f) for r, a, f in TR.get(c, [])}
                    ka = {h[0]: h for h in HT["aif"].get(c, []) if h[0] != sref}
                    kt = {h[0]: h for h in HT[tool].get(c, []) if h[0] != sref}
                    for ref in set(ka) & set(kt):
                        if ref in tr and tr[ref][0] >= 0.90:
                            fa, ft = ka[ref][2] - tr[ref][1], (None if kt[ref][2] is None else kt[ref][2] - tr[ref][1])
                            rows.append((S[c]["length_class"], ka[ref][1] - tr[ref][0], fa, kt[ref][1] - tr[ref][0], ft))

                def summ(sel):
                    n = len(sel)
                    if not n:
                        return {"pairs": 0}
                    e = {"pairs": n, "ani_mae_pp_aif": 100 * sum(abs(r[1]) for r in sel) / n,
                         "ani_mae_pp_tool": 100 * sum(abs(r[3]) for r in sel) / n,
                         "ani_bias_pp_aif": 100 * sum(r[1] for r in sel) / n, "ani_bias_pp_tool": 100 * sum(r[3] for r in sel) / n,
                         "ani_aif_closer": sum(abs(r[1]) < abs(r[3]) for r in sel) / n,
                         "ani_tool_closer": sum(abs(r[3]) < abs(r[1]) for r in sel) / n}
                    if all(r[4] is not None for r in sel):
                        e.update({"af_mae_pp_aif": 100 * sum(abs(r[2]) for r in sel) / n,
                                  "af_mae_pp_tool": 100 * sum(abs(r[4]) for r in sel) / n,
                                  "af_bias_pp_aif": 100 * sum(r[2] for r in sel) / n, "af_bias_pp_tool": 100 * sum(r[4] for r in sel) / n,
                                  "af_aif_closer": sum(abs(r[2]) < abs(r[4]) for r in sel) / n,
                                  "af_tool_closer": sum(abs(r[4]) < abs(r[2]) for r in sel) / n})
                    return e
                P[tool] = summ(rows)
                P[tool]["by_length"] = {ln: summ([r for r in rows if r[0] == ln]) for _, _, ln in LENS}
    os.makedirs(OUT, exist_ok=True)
    part = "" if src == stage and tools == avail else ".partial"
    json.dump(res, open(os.path.join(OUT, f"eval_{stage}{part}.json"), "w"), indent=1)
    for tool in tools + list(GEN):
        o = res["q"][tool]
        t1 = o.get("top1", {})
        print(f"{tool:9s} ANI MAE {t1.get('ani_mae_pp', float('nan')):.3f} AF MAE {t1.get('af_mae_pp', float('nan')):.2f} "
              f"search {o.get('search', {})} decision P {o['decision']['precision']:.3f} R {o['decision']['recall']:.3f}",
              flush=True)


def project():
    """time against query bases per tool; projection to the catalogue (707 GB) from the pool run (a random sample
    of genomes), with the fixed cost from the S1 < S2 < S3 fit; aif's projection is checked against its full run"""
    import numpy as np
    fp = os.path.join(C, "full_bp.json")
    if not os.path.exists(fp):
        tot = n = 0
        for p in range(NPART):
            for line in open(os.path.join(U, "fa", f"part{p:03d}.contigs.tsv")):
                ln = int(line.split("\t", 3)[2])
                if ln >= 500:
                    tot += ln
                    n += 1
        json.dump({"bp": tot, "contigs": n}, open(fp, "w"))
    full_bp = json.load(open(fp))["bp"]
    full_s_aif = sum(v["s"] for v in json.load(open(os.path.join(U, "search_time.json"))).values() if isinstance(v, dict) and "s" in v)
    bp = {}
    for st in ("S1", "S2", "S3"):
        bp[st] = sum(len(l) - 1 for l in open(os.path.join(C, "sample", f"{st}.fa")) if l[0] != ">")
    bp["pool"] = sum(int(r["length"]) for r in csv.DictReader(open(os.path.join(C, "pool", "contigs.tsv")), delimiter="\t")
                     if int(r["length"]) >= 500)
    T = {st: json.load(open(os.path.join(C, "runs", st, "time.json"))) for st in bp
         if os.path.exists(os.path.join(C, "runs", st, "time.json"))}
    builds = {}
    for f in os.listdir(C):
        if f.startswith("build_") and f.endswith(".json"):
            builds[f[6:-5]] = json.load(open(os.path.join(C, f)))
    gb = json.load(open(os.path.join(G, "results", "build.json")))
    builds.setdefault("aif", {"build_s": gb["aif_0.88"]["build_s"], "bytes": gb["aif_0.88"]["bytes"], "threads": 32})
    builds.setdefault("skani", {"build_s": gb["skani"]["build_s"], "bytes": gb["skani"]["bytes"], "threads": 32})
    out = {"bp": bp, "full_bp": full_bp, "aif_full_measured_s": full_s_aif, "builds": builds, "tools": {}}
    tools = sorted({t for st in T for t in T[st]})
    for t in tools:
        pts = [(bp[st], T[st][t]["s"], T[st][t]["kb"]) for st in ("S1", "S2", "S3") if st in T and t in T[st]]
        rec = {"points": {st: T[st][t] for st in T if t in T[st]}}
        if len(pts) >= 2:
            x, y = np.array([p[0] for p in pts], float), np.array([p[1] for p in pts], float)
            b, a = np.polyfit(x, y, 1)
            rec["fit_fixed_s"], rec["fit_s_per_gbp"] = max(0.0, float(a)), float(b) * 1e9
        if "pool" in T and t in T["pool"]:
            fixed = rec.get("fit_fixed_s", 0.0)
            rate = max(1e-12, (T["pool"][t]["s"] - fixed) / bp["pool"])
            rec["pool_s_per_gbp"] = rate * 1e9
            rec["projected_full_s"] = fixed + rate * full_bp
            rec["basis"] = "pool rate, S1-S3 fixed cost"
        elif len(pts) >= 2:  # no pool run: the S1-S3 slope (S3 over-represents flagged categories)
            rec["projected_full_s"] = rec["fit_fixed_s"] + max(0.0, rec["fit_s_per_gbp"]) * full_bp / 1e9
            rec["basis"] = "S1-S3 fit"
        elif pts:  # one subset: all of its time taken as per-base cost (no fixed cost)
            rec["s_per_gbp"] = pts[0][1] / pts[0][0] * 1e9
            rec["projected_full_s"] = pts[0][1] / pts[0][0] * full_bp
            rec["basis"] = f"{next(st for st in ('S1', 'S2', 'S3') if st in T and t in T[st])} rate, no fixed cost"
        if rec["points"]:
            rec["peak_gb"] = max(v["kb"] for v in rec["points"].values()) / 2 ** 20
        out["tools"][t] = rec
    json.dump(out, open(os.path.join(OUT, "cost.json"), "w"), indent=1)
    for t, r in out["tools"].items():
        print(f"{t:9s} fixed {r.get('fit_fixed_s', float('nan')):8.1f} s  pool {r.get('pool_s_per_gbp', float('nan')):9.1f} s/Gbp  "
              f"projected {r.get('projected_full_s', float('nan')) / 3600:8.2f} h  peak {r.get('peak_gb', float('nan')):.1f} GB")
    print(f"aif measured on the full catalogue: {full_s_aif / 3600:.2f} h", flush=True)


def run(stage, tools):
    """stage: pool | S1 | S2 | S3; appends {tool: {s, kb}} to data/compare/runs/<stage>/time.json"""
    fa = os.path.join(C, "pool", "pool.fa") if stage == "pool" else os.path.join(C, "sample", f"{stage}.fa")
    d = os.path.join(C, "runs", stage)
    os.makedirs(d, exist_ok=True)
    tp = os.path.join(d, "time.json")
    T = json.load(open(tp)) if os.path.exists(tp) else {}
    for tool in tools:
        if tool in T and not os.environ.get("RERUN"):
            continue
        r = run_tool(tool, fa, d)
        r.update({"threads": int(THREADS), "cpus": CPUS or "unpinned"})
        T[tool] = r
        json.dump(T, open(tp, "w"), indent=1)
        print(stage, tool, r, flush=True)


if __name__ == "__main__":
    os.makedirs(C, exist_ok=True)
    cmd = sys.argv[1]
    if cmd == "build":
        {"mm2": build_mm2, "lexicmap": build_lexicmap, "sourmash": build_sourmash, "sourmash_index": build_sourmash_index,
         "cat": build_cat}[sys.argv[2]]()
    elif cmd == "pool":
        sample_pool()
    elif cmd == "strata":
        strata()
    elif cmd == "truth":
        truth(sys.argv[2] if len(sys.argv) > 2 else "S3")
    elif cmd == "evaluate":
        evaluate(sys.argv[2] if len(sys.argv) > 2 else "S3")
    elif cmd == "project":
        project()
    elif cmd == "seqmap":
        print(len(seqmap()), "sequences")
    elif cmd == "run":
        run(sys.argv[2], sys.argv[3].split(","))
    else:
        raise SystemExit(__doc__)
