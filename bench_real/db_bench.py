#!/usr/bin/env python3
"""Database-mode benchmark: aif db vs skani search on the real contig sets.

Database: 2,005 genomes (bench_real/db_select.py, db_download.sh) = the five pairwise references,
100 near neighbours of each query species and 1,600 background reference genomes. Nested
scaling databases keep the 405 references + near neighbours and add 0/400/800/1600 background
genomes (db405, db805, db1205, db2005). Queries: the contigs of the five within-species real sets
(the query assemblies are excluded from the database).

  python3 bench_real/db_bench.py build      # aif db build / skani sketch for every size (8 threads), timed
  python3 bench_real/db_bench.py search     # accuracy runs on db2005: aif -n 5 (default and size-matched),
                                            # aif --no-screen -n 1, skani search, skani dist --rl (all 2,005)
  python3 bench_real/db_bench.py truth      # MUMmer truth for every (contig, reference) among both tools'
                                            # top-3 hits and the exhaustive aif best hit
  python3 bench_real/db_bench.py time       # search time / memory, every size, 1 and 8 threads (warm cache)
  python3 bench_real/db_bench.py evaluate   # -> data/db/results/db_eval.json
"""
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
B = os.path.join(ROOT, "data", "db")
RES = os.path.join(B, "results")
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
MUMMER = os.environ.get("MUMMER", os.path.join(ROOT, "data", "envs", "mummer", "bin"))
SETS = [("bcereus", "BA52G"), ("ecoli", "MG1655"), ("mag2328", "REFMAG"), ("paer", "PAO1"), ("saur", "NCTC8325")]
SIZES = [0, 400, 800, 1600]  # background genomes added to the 405 references + near neighbours
TOP = 3


def dbname(nbg):
    return f"db{405 + nbg}"


def ref_rows():
    """[(path, name, group, bg_order)] in database order: references, near neighbours, background."""
    meta = {r["accession"]: (r["group"], int(r["bg_order"])) for r in csv.DictReader(open(os.path.join(B, "meta", "db_accessions.tsv")), delimiter="\t")}
    rows = []
    for line in open(os.path.join(B, "meta", "refs.tsv")):
        path, name = line.rstrip("\n").split("\t")
        g, o = meta.get(name, ("reference", -1))
        rows.append((path, name, g, o))
    return rows


def subset(nbg):
    return [r for r in ref_rows() if r[2] != "background" or r[3] < nbg]


def timed(cmd, out=None, err=None):
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=fo, stderr=fe, check=True)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def dir_bytes(d):
    return sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))


def queries():
    return [(s, os.path.join(ROOT, "data", "real", "sets", s, f"{r}.contigs.fa")) for s, r in SETS]


def all_queries():
    p = os.path.join(B, "queries", "all5.fa")
    if not os.path.exists(p):
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            for _, q in queries():
                fh.write(open(q).read())
    return p


# ------------------------------------------------------------------ build
def build():
    os.makedirs(os.path.join(B, "aif"), exist_ok=True)
    os.makedirs(os.path.join(B, "skani"), exist_ok=True)
    out = {}
    for nbg in SIZES:
        name = dbname(nbg)
        rows = subset(nbg)
        lst = os.path.join(B, "meta", f"{name}.list")
        with open(lst, "w") as fh:
            for p, n, _, _ in rows:
                fh.write(f"{p}\t{n}\n")
        slst = os.path.join(B, "meta", f"{name}.skani.list")
        with open(slst, "w") as fh:
            for p, _, _, _ in rows:
                fh.write(p + "\n")
        adb = os.path.join(B, "aif", name + ".aifd")
        ta, ma = timed([AIF, "db", "build", "--list", lst, "-o", adb, "-k", "20", "-t", "8"], err=adb + ".log")
        sdb = os.path.join(B, "skani", name)
        subprocess.run(["rm", "-rf", sdb], check=True)
        ts, ms = timed([SKANI, "sketch", "-l", slst, "-o", sdb, "-t", "8"], err=sdb + ".log")
        info = subprocess.run([AIF, "db", "info", adb], capture_output=True, text=True, check=True).stdout.split()
        kv = dict(x.split("=", 1) for x in info if "=" in x)
        bases = int(kv["bases"])
        out[name] = {"refs": len(rows), "bases": bases,
                     "aif": {"build_s": ta, "build_kb": ma, "bytes": os.path.getsize(adb),
                             "filter_bytes": int(kv["filter_bytes"]), "screen_bytes": int(kv["screen_bytes"])},
                     "skani": {"build_s": ts, "build_kb": ms, "bytes": dir_bytes(sdb)}}
        o = out[name]
        print(f"{name}: {len(rows)} refs, {bases / 1e9:.2f} Gbp | aif {ta:.1f} s {ma / 1e6:.2f} GB, {o['aif']['bytes'] / 1e6:.0f} MB "
              f"({8 * o['aif']['bytes'] / bases:.3f} bits/base) | skani {ts:.1f} s {ms / 1e6:.2f} GB, {o['skani']['bytes'] / 1e6:.0f} MB "
              f"({8 * o['skani']['bytes'] / bases:.3f} bits/base)", flush=True)
    # size-matched variant of the full database: filter budget 0.74 + screening ~ skani's bits/base
    name = dbname(SIZES[-1])
    adb = os.path.join(B, "aif", name + "_b074.aifd")
    ta, ma = timed([AIF, "db", "build", "--list", os.path.join(B, "meta", f"{name}.list"), "-o", adb, "-k", "20", "-t", "8",
                    "--budget", "0.74"], err=adb + ".log")
    out[name]["aif_b074"] = {"build_s": ta, "build_kb": ma, "bytes": os.path.getsize(adb)}
    print(f"{name} budget 0.74: {ta:.1f} s, {os.path.getsize(adb) / 1e6:.0f} MB ({8 * os.path.getsize(adb) / out[name]['bases']:.3f} bits/base)")
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(os.path.join(RES, "build.json"), "w"), indent=1)


# ------------------------------------------------------------------ search (accuracy runs)
def search():
    name = dbname(SIZES[-1])
    adb = os.path.join(B, "aif", name + ".aifd")
    sdb = os.path.join(B, "skani", name)
    for d in ("aif", "aif_b074", "aif_exh", "skani", "skani_dist"):
        os.makedirs(os.path.join(RES, d), exist_ok=True)
    slst = os.path.join(B, "meta", f"{name}.skani.list")
    for s, q in queries():
        if "--dist-only" in sys.argv:
            out = os.path.join(RES, "skani_dist", s + ".tsv")
            if not os.path.exists(out):
                subprocess.run([SKANI, "dist", "--qi", "-q", q, "--rl", slst, "-t", "16", "-o", out], stderr=subprocess.DEVNULL, check=True)
            print("searched (skani dist)", s, flush=True)
            continue
        subprocess.run([AIF, "db", "search", adb, q, "-n", "5", "-t", "16"], stdout=open(os.path.join(RES, "aif", s + ".tsv"), "w"),
                       stderr=subprocess.DEVNULL, check=True)
        subprocess.run([AIF, "db", "search", adb.replace(".aifd", "_b074.aifd"), q, "-n", "5", "-t", "16"],
                       stdout=open(os.path.join(RES, "aif_b074", s + ".tsv"), "w"), stderr=subprocess.DEVNULL, check=True)
        subprocess.run([AIF, "db", "search", adb, q, "--no-screen", "-n", "1", "-t", "48"],
                       stdout=open(os.path.join(RES, "aif_exh", s + ".tsv"), "w"), stderr=subprocess.DEVNULL, check=True)
        subprocess.run([SKANI, "search", "-d", sdb, "--qi", q, "-t", "16", "-o", os.path.join(RES, "skani", s + ".tsv")],
                       stderr=subprocess.DEVNULL, check=True)
        subprocess.run([SKANI, "dist", "--qi", "-q", q, "--rl", slst, "-t", "16", "-o", os.path.join(RES, "skani_dist", s + ".tsv")],
                       stderr=subprocess.DEVNULL, check=True)
        print("searched", s, flush=True)


def ref_name_of_path():
    return {p: n for p, n, _, _ in ref_rows()}


def load_aif(path):
    """{contig: [(ref, ani, af, rank_score, ci...)] in rank order} (status ok rows only)."""
    H = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = r["contig"].split()[0]
        H.setdefault(c, [])
        if r["status"] == "ok" and r["reference"] != "*":
            H[c].append((r["reference"], float(r["ani"]), float(r["af"]), float(r["ani_lo"]), float(r["ani_hi"]),
                         float(r["af_lo"]), float(r["af_hi"])))
    return H


def load_skani(path):
    """{contig: [(ref, ani, af)] sorted by ANI (skani's ranking)}."""
    pn = ref_name_of_path()
    H = {}
    if not os.path.exists(path):
        return H
    for line in list(open(path))[1:]:
        p = line.rstrip("\n").split("\t")
        c = p[6].split()[0]
        H.setdefault(c, []).append((pn.get(p[0], os.path.basename(p[0])), float(p[2]) / 100, float(p[4]) / 100))
    for c in H:
        H[c].sort(key=lambda x: -x[1])
    return H


# ------------------------------------------------------------------ truth
def contig_seqs(path):
    S, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                S[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        S[name] = "".join(buf)
    return S


def truth_one(job):
    ref, rpath, qfa = job
    import truth as T  # bench_real/truth.py helpers (coords, union_len)
    T.MUMMER = MUMMER
    with tempfile.TemporaryDirectory() as td:
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre, rpath, qfa],
                       check=True, capture_output=True)
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
    for s, q in queries():
        A = load_aif(os.path.join(RES, "aif", s + ".tsv"))
        A2 = load_aif(os.path.join(RES, "aif_b074", s + ".tsv"))
        X = load_aif(os.path.join(RES, "aif_exh", s + ".tsv"))
        S = load_skani(os.path.join(RES, "skani", s + ".tsv"))
        S2 = load_skani(os.path.join(RES, "skani_dist", s + ".tsv"))
        for H in (A, A2, S, S2):
            for c, hs in H.items():
                for h in hs[:TOP]:
                    pairs.setdefault(h[0], set()).add((s, c))
        for c, hs in X.items():
            for h in hs[:1]:
                pairs.setdefault(h[0], set()).add((s, c))
    done = set()
    tp = os.path.join(RES, "truth_pairs.tsv")
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp), delimiter="\t"):
            done.add((r["contig"], r["reference"]))
    seqs = {}
    for s, q in queries():
        seqs.update(contig_seqs(q))
    paths = {n: p for p, n, _, _ in ref_rows()}
    jobs = []
    tmpd = os.path.join(B, "tmp", "truthq")
    os.makedirs(tmpd, exist_ok=True)
    for ref, cs in sorted(pairs.items(), key=lambda x: -len(x[1])):
        todo = sorted(c for _, c in cs if (c, ref) not in done)
        if not todo:
            continue
        qfa = os.path.join(tmpd, ref + ".fa")
        with open(qfa, "w") as fh:
            for c in todo:
                fh.write(f">{c}\n{seqs[c]}\n")
        jobs.append((ref, paths[ref], qfa))
    print(f"truth: {sum(len(v) for v in pairs.values())} pairs over {len(pairs)} references; {len(jobs)} nucmer jobs", flush=True)
    new = not os.path.exists(tp)
    with open(tp, "a") as fh, Pool(int(os.environ.get("JOBS", "48"))) as pool:
        if new:
            fh.write("contig\treference\tlength\ttrue_af\ttrue_ani\n")
        for k, rows in enumerate(pool.imap_unordered(truth_one, jobs, chunksize=1)):
            for c, ref, ln, af, ani in rows:
                fh.write(f"{c}\t{ref}\t{ln}\t{af:.5f}\t{ani:.5f}\n")
            fh.flush()
            if (k + 1) % 50 == 0:
                print(f"  {k + 1}/{len(jobs)}", flush=True)


# ------------------------------------------------------------------ timing
def time_runs():
    q = all_queries()
    out = {}
    tp = os.path.join(RES, "time.json")
    if os.path.exists(tp):
        out = json.load(open(tp))
    for nbg in SIZES:
        name = dbname(nbg)
        adb = os.path.join(B, "aif", name + ".aifd")
        sdb = os.path.join(B, "skani", name)
        for th in (8, 1):
            for tool, cmd in (("aif", [AIF, "db", "search", adb, q, "-t", str(th)]),
                              ("skani", [SKANI, "search", "-d", sdb, "--qi", q, "-t", str(th)]),
                              ("skani_keep", [SKANI, "search", "-d", sdb, "--qi", q, "-t", str(th), "--keep-refs"]),
                              ("skani_dist", [SKANI, "dist", "--qi", "-q", q, "--rl", os.path.join(B, "meta", f"{name}.skani.list"),
                                              "-t", str(th)])):
                key = f"{name}.{tool}.t{th}"
                if key in out or (tool == "skani_dist" and th == 1):  # skani dist re-sketches every reference per run
                    continue
                if th == 1 and nbg not in (SIZES[0], SIZES[-1]):  # single-thread runs: smallest and largest database
                    continue
                cpu = f"0-{th - 1}" if th > 1 else "2"
                if th > 1:
                    timed(["taskset", "-c", cpu] + cmd)  # warm the page cache (8-thread runs come first)
                t, m = timed(["taskset", "-c", cpu] + cmd)
                out[key] = {"s": t, "kb": m}
                print(f"{key}: {t:.1f} s, {m / 1e6:.2f} GB", flush=True)
                json.dump(out, open(tp, "w"), indent=1)


# ------------------------------------------------------------------ evaluation
def evaluate():
    TR = {}
    for r in csv.DictReader(open(os.path.join(RES, "truth_pairs.tsv")), delimiter="\t"):
        TR[(r["contig"], r["reference"])] = (float(r["true_ani"]), float(r["true_af"]))
    res = {"sets": {}}
    pooled = {}
    for s, q in queries():
        own = dict(SETS)[s]
        T0 = {r["contig"]: (float(r["true_af"]), float(r["true_ani"]), int(r["length"]))
              for r in csv.DictReader(open(os.path.join(ROOT, "data", "real", "sets", s, "truth.tsv")), delimiter="\t")}
        A = load_aif(os.path.join(RES, "aif", s + ".tsv"))
        A2 = load_aif(os.path.join(RES, "aif_b074", s + ".tsv"))
        X = load_aif(os.path.join(RES, "aif_exh", s + ".tsv"))
        S = load_skani(os.path.join(RES, "skani", s + ".tsv"))
        S2 = load_skani(os.path.join(RES, "skani_dist", s + ".tsv"))
        names = sorted(T0)
        hom = [c for c in names if T0[c][0] >= 0.2]  # homologous to the pairwise reference (in the database)

        def score(c, ref):
            t = TR.get((c, ref))
            if t is None or t[1] <= 0 or t[0] != t[0]:
                return 0.0
            return t[0] * t[1]

        def cands_of(c):
            return set(h[0] for H in (A, A2, S, S2) for h in H.get(c, [])[:TOP]) | set(h[0] for h in X.get(c, [])[:1]) | {own}

        def best(c):
            return max(score(c, r) for r in cands_of(c))

        def best_ani(c):  # highest true ANI among candidate references aligned over >= half the contig
            v = [TR[(c, r)][0] for r in cands_of(c) if (c, r) in TR and TR[(c, r)][1] >= 0.5 and TR[(c, r)][0] == TR[(c, r)][0]]
            return max(v) if v else None

        m = {}
        for tool, H in (("aif", A), ("aif_b074", A2), ("skani", S), ("skani_dist", S2), ("aif_exhaustive", X)):
            rep = [c for c in hom if H.get(c)]
            e_ani, e_af, e_af_all, reg, isbest, cov_ani, cov_af, reg_ani = [], [], [], [], [], [], [], []
            for c in hom:
                hs = H.get(c, [])
                if not hs:
                    e_af_all.append(T0[c][0])  # nothing reported: AF 0 vs the best available
                    continue
                h = hs[0]
                t = TR.get((c, h[0]))
                if t is None:
                    continue
                ta, tf = t
                e_af.append(abs(h[2] - tf))
                e_af_all.append(abs(h[2] - tf))
                if tf >= 0.2 and ta == ta:
                    e_ani.append(abs(h[1] - ta))
                    if len(h) > 3:
                        cov_ani.append(h[3] - 1e-9 <= ta <= h[4] + 1e-9)
                        cov_af.append(h[5] - 1e-9 <= tf <= h[6] + 1e-9)
                b = best(c)
                sc = score(c, h[0])
                reg.append(b - sc)
                isbest.append(sc >= b - 1e-3)
                ba = best_ani(c)
                if ba is not None and tf >= 0.2 and ta == ta:
                    reg_ani.append(max(0.0, ba - ta))
            mean = lambda v: sum(v) / len(v) if v else float("nan")  # noqa: E731
            m[tool] = {"n_hom": len(hom), "coverage": len(rep) / max(1, len(hom)),
                       "ani_mae_pp": 100 * mean(e_ani), "n_ani": len(e_ani), "af_mae_pp": 100 * mean(e_af),
                       "af_mae_all_pp": 100 * mean(e_af_all),
                       "regret_pp": 100 * mean(reg), "top1_best": mean(isbest), "regret_ani_pp": 100 * mean(reg_ani),
                       "ci_cov_ani": mean(cov_ani) if cov_ani else None, "ci_cov_af": mean(cov_af) if cov_af else None}
            for k, v in (("e_ani", e_ani), ("e_af", e_af), ("e_af_all", e_af_all), ("reg", reg), ("isbest", isbest),
                         ("cov_ani", cov_ani), ("cov_af", cov_af), ("reg_ani", reg_ani)):
                pooled.setdefault(tool, {}).setdefault(k, []).extend(v)
            pooled[tool].setdefault("n_hom", 0)
            pooled[tool]["n_hom"] += len(hom)
            pooled[tool].setdefault("n_rep", 0)
            pooled[tool]["n_rep"] += len(rep)
        # paired ANI error on contigs both tools report with a homologous top-1
        for key, SS in (("paired", S), ("paired_dist", S2)):
            both = []
            for c in hom:
                if A.get(c) and SS.get(c):
                    ta_ = TR.get((c, A[c][0][0]))
                    ts_ = TR.get((c, SS[c][0][0]))
                    if ta_ and ts_ and ta_[1] >= 0.2 and ts_[1] >= 0.2:
                        both.append((abs(A[c][0][1] - ta_[0]), abs(SS[c][0][1] - ts_[0])))
            m[key] = {"n": len(both), "aif_ani_mae_pp": 100 * sum(x for x, _ in both) / max(1, len(both)),
                      "skani_ani_mae_pp": 100 * sum(y for _, y in both) / max(1, len(both))}
            pooled.setdefault(key, []).extend(both)
        # screening recall: exhaustive best hit found by the screened search (same reference or same score)
        rec = [c for c in names if X.get(c)]
        found = sum(1 for c in rec if A.get(c) and (A[c][0][0] == X[c][0][0] or
                                                    A[c][0][1] * A[c][0][2] >= X[c][0][1] * X[c][0][2] - 1e-3))
        m["screen_recall"] = found / max(1, len(rec))
        m["screen_recall_n"] = len(rec)
        res["sets"][s] = m
    P = {}
    mean = lambda v: sum(v) / len(v) if v else float("nan")  # noqa: E731
    for tool in ("aif", "aif_b074", "skani", "skani_dist", "aif_exhaustive"):
        d = pooled[tool]
        P[tool] = {"n_hom": d["n_hom"], "coverage": d["n_rep"] / d["n_hom"], "ani_mae_pp": 100 * mean(d["e_ani"]),
                   "n_ani": len(d["e_ani"]), "af_mae_pp": 100 * mean(d["e_af"]), "af_mae_all_pp": 100 * mean(d["e_af_all"]),
                   "regret_pp": 100 * mean(d["reg"]), "top1_best": mean(d["isbest"]), "regret_ani_pp": 100 * mean(d["reg_ani"]),
                   "ci_cov_ani": mean(d["cov_ani"]) if d["cov_ani"] else None,
                   "ci_cov_af": mean(d["cov_af"]) if d["cov_af"] else None}
    for key in ("paired", "paired_dist"):
        both = pooled[key]
        P[key] = {"n": len(both), "aif_ani_mae_pp": 100 * sum(x for x, _ in both) / max(1, len(both)),
                  "skani_ani_mae_pp": 100 * sum(y for _, y in both) / max(1, len(both))}
    tot = sum(res["sets"][s]["screen_recall_n"] for s, _ in SETS)
    P["screen_recall"] = sum(res["sets"][s]["screen_recall"] * res["sets"][s]["screen_recall_n"] for s, _ in SETS) / max(1, tot)
    res["pooled"] = P
    for f in ("build.json", "time.json"):
        p = os.path.join(RES, f)
        if os.path.exists(p):
            res[f[:-5]] = json.load(open(p))
    json.dump(res, open(os.path.join(RES, "db_eval.json"), "w"), indent=1)
    hdr = f"{'set':9s} {'tool':15s} {'cov':>6s} {'ANI MAE':>8s} {'AF MAE':>7s} {'AF all':>7s} {'regret':>7s} {'top1best':>8s} {'CI ANI':>6s} {'CI AF':>6s} {'ANIreg':>6s}"
    print(hdr)
    for s in [x for x, _ in SETS] + ["pooled"]:
        m = res["sets"][s] if s != "pooled" else P
        for tool in ("aif", "aif_b074", "aif_exhaustive", "skani", "skani_dist"):
            r = m[tool]
            ca = f"{r['ci_cov_ani']:.3f}" if r.get("ci_cov_ani") is not None else "-"
            cf = f"{r['ci_cov_af']:.3f}" if r.get("ci_cov_af") is not None else "-"
            print(f"{s:9s} {tool:15s} {r['coverage']:6.3f} {r['ani_mae_pp']:8.3f} {r['af_mae_pp']:7.2f} {r['af_mae_all_pp']:7.2f} "
                  f"{r['regret_pp']:7.3f} {r['top1_best']:8.3f} {ca:>6s} {cf:>6s} {r['regret_ani_pp']:6.3f}")
        pr, pd = m["paired"], m["paired_dist"]
        print(f"{s:9s} paired ANI MAE vs search (n={pr['n']}): aif {pr['aif_ani_mae_pp']:.3f} skani {pr['skani_ani_mae_pp']:.3f}; "
              f"vs dist (n={pd['n']}): aif {pd['aif_ani_mae_pp']:.3f} skani {pd['skani_ani_mae_pp']:.3f}; screen recall {m['screen_recall']:.4f}")


if __name__ == "__main__":
    {"build": build, "search": search, "truth": truth, "time": time_runs, "evaluate": evaluate}[sys.argv[1]]()
