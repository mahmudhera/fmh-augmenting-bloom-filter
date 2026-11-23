#!/usr/bin/env python3
"""GTDB-scale database and contig-classification benchmark (aif db vs skani search).

Database: all GTDB species representatives of the latest release (bench_real/gtdb_download.sh).
Queries: contigs >= 500 bp of 294 non-representative draft genomes from distinct species
(bench_real/gtdb_select.py, gtdb_prepare.py); each query's species representative is in the
database. Truth: GTDB taxonomy (species/genus of the top hit) and MUMmer for every
(contig, reference) pair among each tool's top 3 hits.

  python3 bench_real/gtdb_bench.py list       # reference list from the unpacked tarball
  python3 bench_real/gtdb_bench.py build      # aif db build (budgets 0.88, 0.22) and skani sketch, timed
  python3 bench_real/gtdb_bench.py search     # aif (both databases) and skani search, timed
  python3 bench_real/gtdb_bench.py truth      # MUMmer truth for top-hit pairs
  python3 bench_real/gtdb_bench.py evaluate   # -> data/gtdb/results/gtdb_eval.json
"""
import csv
import glob
import gzip
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
RES = os.path.join(G, "results")
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
MUMMER = os.environ.get("MUMMER", os.path.join(ROOT, "data", "envs", "mummer", "bin"))
THREADS = int(os.environ.get("THREADS", "48"))
BUDGETS = ["0.88", "0.22"]
TOP = 3
Q = os.path.join(G, "queries", "contigs.fa")


def timed(cmd, out=None, err=None):
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=fo, stderr=fe, check=True)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def taxonomy():
    """{ncbi accession: (species, genus, domain)} for every GTDB genome."""
    T = {}
    for fn, dom in (("bac120_metadata.tsv.gz", "bac"), ("ar53_metadata.tsv.gz", "ar")):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tax = r["gtdb_taxonomy"].split(";")
                T[r["accession"][3:]] = (tax[-1][3:], tax[-2][3:], dom)
    return T


def acc_of(path):
    b = os.path.basename(path)
    return b[: b.index("_genomic")] if "_genomic" in b else b.split(".fna")[0]


def list_refs():
    files = sorted(glob.glob(os.path.join(G, "reps", "*", "database", "GC*", "*", "*", "*", "*_genomic.fna.gz")))
    with open(os.path.join(G, "meta", "reps.list"), "w") as fh, open(os.path.join(G, "meta", "reps.skani.list"), "w") as fs:
        for p in files:
            fh.write(f"{p}\t{acc_of(p)}\n")
            fs.write(p + "\n")
    print(f"{len(files)} representative genomes")


def build():
    os.makedirs(os.path.join(G, "aif"), exist_ok=True)
    out = {}
    bp = os.path.join(RES, "build.json")
    os.makedirs(RES, exist_ok=True)
    if os.path.exists(bp):
        out = json.load(open(bp))
    lst = os.path.join(G, "meta", "reps.list")
    for b in BUDGETS:
        key = f"aif_{b}"
        if key in out:
            continue
        db = os.path.join(G, "aif", f"gtdb_{b}.aifd")
        t, m = timed([AIF, "db", "build", "--list", lst, "-o", db, "-k", "20", "-t", str(THREADS), "--budget", b], err=db + ".log")
        info = subprocess.run([AIF, "db", "info", db], capture_output=True, text=True, check=True).stdout.split()
        kv = dict(x.split("=", 1) for x in info if "=" in x)
        out[key] = {"build_s": t, "build_kb": m, "bytes": os.path.getsize(db), "bases": int(kv["bases"]), "refs": int(kv["references"]),
                    "filter_bytes": int(kv["filter_bytes"]), "screen_bytes": int(kv["screen_bytes"])}
        print(key, out[key], flush=True)
        json.dump(out, open(bp, "w"), indent=1)
    if "skani" not in out:
        sdb = os.path.join(G, "skani_db")
        subprocess.run(["rm", "-rf", sdb], check=True)
        t, m = timed([SKANI, "sketch", "-l", os.path.join(G, "meta", "reps.skani.list"), "-o", sdb, "-t", str(THREADS)],
                     err=sdb + ".log")
        sz = sum(os.path.getsize(os.path.join(sdb, f)) for f in os.listdir(sdb))
        out["skani"] = {"build_s": t, "build_kb": m, "bytes": sz}
        print("skani", out["skani"], flush=True)
        json.dump(out, open(bp, "w"), indent=1)


def search():
    os.makedirs(RES, exist_ok=True)
    tp = os.path.join(RES, "search_time.json")
    out = json.load(open(tp)) if os.path.exists(tp) else {}
    th = os.environ.get("SEARCH_THREADS", "16")
    runs = [(f"aif_{b}", [AIF, "db", "search", os.path.join(G, "aif", f"gtdb_{b}.aifd"), Q, "-n", "5", "-t", th]) for b in BUDGETS]
    runs.append(("skani", [SKANI, "search", "-d", os.path.join(G, "skani_db"), "--qi", Q, "-t", th]))
    for key, cmd in runs:
        if key in out:
            continue
        o = os.path.join(RES, key + ".tsv")
        if key == "skani":
            cmd = cmd + ["-o", o]
            t, m = timed(cmd, err=o + ".log")
        else:
            t, m = timed(cmd, out=o, err=o + ".log")
        out[key] = {"s": t, "kb": m, "threads": int(th)}
        print(key, out[key], flush=True)
        json.dump(out, open(tp, "w"), indent=1)


def load_aif(path):
    H = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = r["contig"].split()[0]
        H.setdefault(c, [])
        if r["status"] == "ok" and r["reference"] != "*":
            H[c].append((r["reference"], float(r["ani"]), float(r["af"]), float(r["ani_lo"]), float(r["ani_hi"]),
                         float(r["af_lo"]), float(r["af_hi"])))
    return H


def load_skani(path):
    H = {}
    if not os.path.exists(path):
        return H
    for line in list(open(path))[1:]:
        p = line.rstrip("\n").split("\t")
        H.setdefault(p[6].split()[0], []).append((acc_of(p[0]), float(p[2]) / 100, float(p[4]) / 100))
    for c in H:
        H[c].sort(key=lambda x: -x[1])
    return H


def truth_one(job):
    ref, rpath, qfa = job
    import truth as T
    T.MUMMER = MUMMER
    with tempfile.TemporaryDirectory() as td:
        rfa = os.path.join(td, "ref.fa")
        with gzip.open(rpath, "rb") as fi, open(rfa, "wb") as fo:
            fo.write(fi.read())
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre, rfa, qfa], check=True, capture_output=True)
        cq = T.coords(pre + ".delta", "-q")
    L = T.contig_lengths(qfa)
    iv, wid, wl = {}, {}, {}
    for name, s, e, lq, idy, _ in cq:
        iv.setdefault(name, []).append((s, e))
        wid[name] = wid.get(name, 0.0) + idy * lq
        wl[name] = wl.get(name, 0) + lq
    rows = []
    for name, ln in L.items():
        af = min(1.0, T.union_len(iv.get(name, [])) / ln)
        ani = wid[name] / wl[name] / 100 if wl.get(name) else float("nan")
        rows.append((name, ref, ln, af, ani))
    os.unlink(qfa)
    return rows


def truth():
    pairs = {}
    for f in [f"aif_{b}" for b in BUDGETS] + ["skani"]:
        H = load_skani(os.path.join(RES, f + ".tsv")) if f == "skani" else load_aif(os.path.join(RES, f + ".tsv"))
        for c, hs in H.items():
            for h in hs[:TOP]:
                pairs.setdefault(h[0], set()).add(c)
    meta = {r["contig"]: r for r in csv.DictReader(open(os.path.join(G, "meta", "query_contigs.tsv")), delimiter="\t")}
    for c, m in meta.items():  # the query's own species representative, always
        pairs.setdefault(m["representative"][3:], set()).add(c)
    tp = os.path.join(RES, "truth_pairs.tsv")
    done = set()
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp), delimiter="\t"):
            done.add((r["contig"], r["reference"]))
    paths = {acc: p for p, acc in (l.rstrip("\n").split("\t") for l in open(os.path.join(G, "meta", "reps.list")))}
    seqs, name, buf = {}, None, []
    for line in open(Q):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    tmpd = os.path.join(G, "tmp", "truthq")
    os.makedirs(tmpd, exist_ok=True)
    jobs = []
    for ref, cs in sorted(pairs.items(), key=lambda x: -len(x[1])):
        todo = sorted(c for c in cs if (c, ref) not in done)
        if not todo or ref not in paths:
            continue
        qfa = os.path.join(tmpd, ref + ".fa")
        with open(qfa, "w") as fh:
            for c in todo:
                fh.write(f">{c}\n{seqs[c]}\n")
        jobs.append((ref, paths[ref], qfa))
    print(f"truth: {sum(len(v) for v in pairs.values())} pairs over {len(pairs)} references; {len(jobs)} nucmer jobs", flush=True)
    new = not os.path.exists(tp)
    with open(tp, "a") as fh, Pool(int(os.environ.get("JOBS", "56"))) as pool:
        if new:
            fh.write("contig\treference\tlength\ttrue_af\ttrue_ani\n")
        for k, rows in enumerate(pool.imap_unordered(truth_one, jobs, chunksize=1)):
            for c, ref, ln, af, ani in rows:
                fh.write(f"{c}\t{ref}\t{ln}\t{af:.5f}\t{ani:.5f}\n")
            fh.flush()
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(jobs)}", flush=True)


def evaluate():
    TX = taxonomy()
    meta = {r["contig"]: r for r in csv.DictReader(open(os.path.join(G, "meta", "query_contigs.tsv")), delimiter="\t")}
    TR = {}
    tp = os.path.join(RES, "truth_pairs.tsv")
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp), delimiter="\t"):
            TR[(r["contig"], r["reference"])] = (float(r["true_ani"]), float(r["true_af"]))
    LB = [(500, 1000, "0.5-1k"), (1000, 2000, "1-2k"), (2000, 5000, "2-5k"), (5000, 20000, "5-20k"), (20000, 10 ** 12, ">=20k")]
    res = {"tools": {}}
    tools = [(f"aif_{b}", load_aif(os.path.join(RES, f"aif_{b}.tsv"))) for b in BUDGETS]
    tools.append(("skani", load_skani(os.path.join(RES, "skani.tsv"))))
    mean = lambda v: sum(v) / len(v) if v else float("nan")  # noqa: E731
    for name, H in tools:
        recs = []
        for c, m in meta.items():
            L = int(m["length"])
            own = m["representative"][3:]
            tself = TR.get((c, own))
            homol = tself is not None and tself[1] >= 0.2  # aligns to its own species representative
            hs = H.get(c, [])
            if not hs:
                recs.append({"L": L, "dom": m["domain"], "rep": 0, "sp": 0, "gen": 0, "homol": homol})
                continue
            h = hs[0]
            sp, gen, _ = TX.get(h[0], ("?", "?", "?"))
            r = {"L": L, "dom": m["domain"], "rep": 1, "sp": int(sp == m["species"]), "gen": int(gen == m["genus"]), "homol": homol,
                 "ani": h[1], "af": h[2]}
            t = TR.get((c, h[0]))
            if t is not None:
                r["t_ani"], r["t_af"] = t
                if len(h) > 3:
                    r["ci_a"] = h[3] - 1e-9 <= t[0] <= h[4] + 1e-9 if t[1] >= 0.2 and t[0] == t[0] else None
                    r["ci_f"] = h[5] - 1e-9 <= t[1] <= h[6] + 1e-9
            recs.append(r)

        def summ(rs):
            n = len(rs)
            rep = [r for r in rs if r["rep"]]
            hom = [r for r in rs if r["homol"]]
            ea = [abs(r["ani"] - r["t_ani"]) for r in rep if "t_ani" in r and r["t_af"] >= 0.2 and r["t_ani"] == r["t_ani"]]
            ef = [abs(r["af"] - r["t_af"]) for r in rep if "t_af" in r]
            ca = [r["ci_a"] for r in rep if r.get("ci_a") is not None]
            cf = [r["ci_f"] for r in rep if r.get("ci_f") is not None]
            return {"n": n, "reported": len(rep) / max(1, n), "species_acc": sum(r["sp"] for r in rs) / max(1, n),
                    "genus_acc": sum(r["gen"] for r in rs) / max(1, n),
                    "species_prec": sum(r["sp"] for r in rep) / max(1, len(rep)),
                    "n_homol": len(hom), "species_acc_homol": sum(r["sp"] for r in hom) / max(1, len(hom)),
                    "ani_mae_pp": 100 * mean(ea), "n_ani": len(ea), "af_mae_pp": 100 * mean(ef), "n_af": len(ef),
                    "ci_cov_ani": mean(ca) if ca else None, "ci_cov_af": mean(cf) if cf else None}

        d = {"all": summ(recs), "bac": summ([r for r in recs if r["dom"] == "bac"]), "ar": summ([r for r in recs if r["dom"] == "ar"])}
        for lo, hi, lab in LB:
            d[lab] = summ([r for r in recs if lo <= r["L"] < hi])
        res["tools"][name] = d
    for f in ("build.json", "search_time.json"):
        p = os.path.join(RES, f)
        if os.path.exists(p):
            res[f[:-5]] = json.load(open(p))
    json.dump(res, open(os.path.join(RES, "gtdb_eval.json"), "w"), indent=1)
    keys = ["all", "bac", "ar"] + [x[2] for x in LB]
    print(f"{'tool':10s} {'subset':8s} {'n':>6s} {'rep':>6s} {'sp acc':>7s} {'gen acc':>7s} {'sp prec':>7s} {'sp|hom':>7s} {'ANI MAE':>8s} {'AF MAE':>7s} {'CI A/F':>11s}")
    for name, d in res["tools"].items():
        for k in keys:
            s = d[k]
            ci = f"{s['ci_cov_ani']:.3f}/{s['ci_cov_af']:.3f}" if s["ci_cov_ani"] is not None else "-"
            print(f"{name:10s} {k:8s} {s['n']:6d} {s['reported']:6.3f} {s['species_acc']:7.3f} {s['genus_acc']:7.3f} {s['species_prec']:7.3f} "
                  f"{s['species_acc_homol']:7.3f} {s['ani_mae_pp']:8.3f} {s['af_mae_pp']:7.2f} {ci:>11s}")


if __name__ == "__main__":
    {"list": list_refs, "build": build, "search": search, "truth": truth, "evaluate": evaluate}[sys.argv[1]]()
