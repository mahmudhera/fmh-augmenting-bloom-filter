#!/usr/bin/env python3
"""Validation of the HGT atlas calls (bench_real/hgt_atlas.py) against alignment.

  precision SET [N]   N sampled cross-genus calls (>= 99%, >= 500 bp), each window aligned to its partner genome with
                      MUMmer (nucmer, delta-filter -q, show-coords): confirmed when >= 500 query bases align at >= 99%
  recall SET          a universe U of 2,000 genomes of the set and 300 query genomes R in U: minimap2 (asm5, base-level,
                      all secondary alignments kept) maps R's windows to U; (window, partner) pairs with >= 500 bases at
                      >= 99% identity are compared with aif's calls restricted to partners in U

Usage: /usr/bin/python3 -I bench_real/hgt_validate.py precision complete 400
       /usr/bin/python3 -I bench_real/hgt_validate.py recall complete
"""
import collections
import csv
import gzip
import json
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
MUM = os.path.join(ROOT, "data", "envs", "mummer", "bin")
MM2 = "/usr/bin/minimap2"
W, WMIN = 5000, 1000
ANI_MIN, BP_MIN = 0.99, 500
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def genomes():
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    for g in G.values():
        g["tax"] = [x[3:] for x in g["taxonomy"].split(";")]
    return G


def read_fasta_gz(path):
    seqs, buf = [], []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l[0] == ">":
                if buf:
                    seqs.append("".join(buf))
                    buf = []
                seqs.append(None)
            else:
                buf.append(l.strip())
    if buf:
        seqs.append("".join(buf))
    return [s for s in seqs if s is not None]


def calls(setname):
    """aif calls: (window, partner) -> (ani, aligned bp)"""
    out = {}
    hd = os.path.join(HD, setname, "hits")
    for f in sorted(os.listdir(hd)):
        if f.endswith(".tsv"):
            for l in open(os.path.join(hd, f)):
                w, ref, rank, ani, af, theta, lr = l.rstrip("\n").split("\t")
                al = float(af) * int(w.split("|")[3])
                if float(ani) >= ANI_MIN and al >= BP_MIN:
                    out[(w, ref)] = (float(ani), al)
    return out


def nucmer_job(job):
    """align one window to one partner genome; aligned query bases at >= 99% and their identity"""
    w, ref, qpath, rpath = job
    acc, ci, st, ln = w.split("|")
    seq = read_fasta_gz(qpath)[int(ci)][int(st):int(st) + int(ln)]
    td = tempfile.mkdtemp(dir=os.path.join(HD, "tmp"))
    try:
        q = os.path.join(td, "q.fa")
        open(q, "w").write(f">w\n{seq}\n")
        r = os.path.join(td, "r.fa")
        with gzip.open(rpath, "rb") as fi, open(r, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        subprocess.run([os.path.join(MUM, "nucmer"), "--maxmatch", "-p", os.path.join(td, "o"), r, q],
                       capture_output=True, check=True)
        flt = subprocess.run([os.path.join(MUM, "delta-filter"), "-q", os.path.join(td, "o.delta")],
                             capture_output=True, check=True, text=True).stdout
        open(os.path.join(td, "f.delta"), "w").write(flt)
        co = subprocess.run([os.path.join(MUM, "show-coords"), "-rclTH", os.path.join(td, "f.delta")],
                            capture_output=True, check=True, text=True).stdout
        cov = [0] * len(seq)
        al = idw = 0.0
        for l in co.splitlines():
            p = l.split("\t")
            s2, e2, l2, idy = int(p[2]), int(p[3]), int(p[5]), float(p[6])
            lo, hi = min(s2, e2), max(s2, e2)
            al += l2
            idw += l2 * idy
            if idy >= 99.0:
                for i in range(lo - 1, hi):
                    cov[i] = 1
        return w, ref, sum(cov), (idw / al / 100 if al else None), al
    finally:
        shutil.rmtree(td, ignore_errors=True)


def precision(setname, n):
    from hgt_analyze import lca_class
    G = genomes()
    C = calls(setname)
    xs = [(w, r) for (w, r) in C if lca_class(G[w.split("|")[0]]["tax"], G[r]["tax"]) >= 1]
    random.seed(20261127)
    samp = random.sample(xs, min(n, len(xs)))
    os.makedirs(os.path.join(HD, "tmp"), exist_ok=True)
    jobs = [(w, r, G[w.split("|")[0]]["path"], G[r]["path"]) for w, r in samp]
    res = []
    with multiprocessing.get_context("fork").Pool(32) as pool:
        for x in pool.imap_unordered(nucmer_job, jobs):
            res.append(x)
    rows = []
    for w, r, cov99, idy, al in res:
        ani, alb = C[(w, r)]
        rows.append({"window": w, "partner": r, "aif_ani": ani, "aif_bp": alb, "mum_bp99": cov99, "mum_identity": idy,
                     "mum_aligned": al, "distance": lca_class(G[w.split("|")[0]]["tax"], G[r]["tax"])})
    conf = [x for x in rows if x["mum_bp99"] >= BP_MIN]
    near = [x for x in rows if x["mum_bp99"] < BP_MIN and x["mum_identity"] and x["mum_identity"] >= 0.98]
    S = {"n": len(rows), "confirmed": len(conf), "precision": len(conf) / len(rows),
         "near_miss_98": len(near), "unaligned": sum(1 for x in rows if not x["mum_aligned"]),
         "ani_mae_pp": 100 * sum(abs(x["aif_ani"] - x["mum_identity"]) for x in conf if x["mum_identity"]) / max(1, len(conf)),
         "bp_median_ratio": sorted(x["aif_bp"] / max(1, x["mum_bp99"]) for x in conf)[len(conf) // 2] if conf else None}
    out = os.path.join(ROOT, "results_hgt", setname)
    os.makedirs(out, exist_ok=True)
    json.dump({"summary": S, "rows": rows}, open(os.path.join(out, "validate_precision.json"), "w"), indent=1)
    print(json.dumps(S, indent=1))


def recall(setname, nu=2000, nr=300):
    from hgt_analyze import lca_class
    G = genomes()
    base = os.path.join(HD, setname)
    members = sorted({l.split("\t")[0] for l in open(os.path.join(base, "contigs.tsv"))})
    random.seed(7)
    U = random.sample(members, min(nu, len(members)))
    R = set(random.sample(U, min(nr, len(U))))
    vd = os.path.join(HD, "validate_" + setname)
    os.makedirs(vd, exist_ok=True)
    uf, rf = os.path.join(vd, "U.fa"), os.path.join(vd, "R_windows.fa")
    if not os.path.exists(uf):
        with open(uf + ".tmp", "w") as fo:
            for a in U:
                for ci, s in enumerate(read_fasta_gz(G[a]["path"])):
                    fo.write(f">{a}|{ci}\n{s}\n")
        os.replace(uf + ".tmp", uf)
    if not os.path.exists(rf):
        with open(rf + ".tmp", "w") as fo:
            for a in sorted(R):
                for ci, s in enumerate(read_fasta_gz(G[a]["path"])):
                    for st in range(0, len(s), W):
                        t = s[st:st + W]
                        if len(t) >= WMIN:
                            fo.write(f">{a}|{ci}|{st}|{len(t)}\n{t}\n")
        os.replace(rf + ".tmp", rf)
    paf = os.path.join(vd, "mm2.paf")
    tm = {}
    if not os.path.exists(paf):
        t0 = time.time()
        subprocess.run(f"/usr/bin/time -f '%e %M' -o {vd}/mm2.time {MM2} -x asm5 -c -N 5000 -p 0.0 --secondary=yes "
                       f"-t 48 -I 100G {uf} {rf} > {paf}.tmp 2> {vd}/mm2.log", shell=True, check=True)
        os.replace(paf + ".tmp", paf)
        tm["mm2_wall"] = time.time() - t0
    e, m = open(os.path.join(vd, "mm2.time")).read().split()
    tm["mm2_s"], tm["mm2_kb"] = float(e), int(m)
    # minimap2 calls: per (window, partner genome) the query bases covered by alignments at >= 99% identity
    cov = collections.defaultdict(list)
    for l in open(paf):
        p = l.split("\t")
        q, t = p[0], p[5].split("|")[0]
        if t == q.split("|")[0]:
            continue
        matches, blen = int(p[9]), int(p[10])
        if blen and matches / blen >= ANI_MIN:
            cov[(q, t)].append((int(p[2]), int(p[3])))
    mm = set()
    for k, iv in cov.items():
        iv.sort()
        tot, cs, ce = 0, -1, -1
        for s, e in iv:
            if s > ce:
                tot += max(0, ce - cs)
                cs, ce = s, e
            else:
                ce = max(ce, e)
        tot += max(0, ce - cs)
        if tot >= BP_MIN:
            mm.add(k)
    C = calls(setname)
    Uset = set(U)
    aif = {(w, r) for (w, r) in C if w.split("|")[0] in R and r in Uset}
    both = aif & mm
    def by_dist(S_):
        c = collections.Counter(lca_class(G[w.split("|")[0]]["tax"], G[r]["tax"]) >= 1 for w, r in S_)
        return {"cross_genus": c[True], "same_genus": c[False]}
    S = {"U": len(U), "R": len(R), "minimap2_calls": len(mm), "aif_calls": len(aif), "both": len(both),
         "recall_vs_minimap2": len(both) / max(1, len(mm)), "aif_also_in_minimap2": len(both) / max(1, len(aif)),
         "mm2_by_distance": by_dist(mm), "both_by_distance": by_dist(both), "aif_by_distance": by_dist(aif), **tm}
    xg_mm = {k for k in mm if lca_class(G[k[0].split("|")[0]]["tax"], G[k[1]]["tax"]) >= 1}
    S["recall_cross_genus"] = len(xg_mm & aif) / max(1, len(xg_mm))
    random.seed(11)
    S["sample_mm2_only"] = random.sample(sorted(mm - aif), min(60, len(mm - aif)))
    S["sample_aif_only"] = random.sample(sorted(aif - mm), min(60, len(aif - mm)))
    out = os.path.join(ROOT, "results_hgt", setname)
    os.makedirs(out, exist_ok=True)
    json.dump(S, open(os.path.join(out, "validate_recall.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in S.items() if not k.startswith("sample")}, indent=1))


def adjudicate(setname):
    """MUMmer on the sampled disagreements of the recall comparison"""
    G = genomes()
    out = os.path.join(ROOT, "results_hgt", setname)
    S = json.load(open(os.path.join(out, "validate_recall.json")))
    os.makedirs(os.path.join(HD, "tmp"), exist_ok=True)
    res = {}
    for key in ("sample_mm2_only", "sample_aif_only"):
        jobs = [(w, r, G[w.split("|")[0]]["path"], G[r]["path"]) for w, r in S[key]]
        with multiprocessing.get_context("fork").Pool(32) as pool:
            rr = list(pool.imap_unordered(nucmer_job, jobs))
        res[key] = {"n": len(rr), "mummer_confirms": sum(1 for x in rr if x[2] >= BP_MIN),
                    "rows": [{"window": x[0], "partner": x[1], "mum_bp99": x[2], "mum_identity": x[3]} for x in rr]}
        print(key, res[key]["n"], res[key]["mummer_confirms"])
    json.dump(res, open(os.path.join(out, "validate_adjudicate.json"), "w"), indent=1)


def mm2_calls(vd):
    """minimap2 (window, partner) pairs with >= BP_MIN query bases at >= 99% identity (5 kb windows)"""
    cov = collections.defaultdict(list)
    for l in open(os.path.join(vd, "mm2.paf")):
        p = l.split("\t")
        q, t = p[0], p[5].split("|")[0]
        if t == q.split("|")[0]:
            continue
        matches, blen = int(p[9]), int(p[10])
        if blen and matches / blen >= ANI_MIN:
            cov[(q, t)].append((int(p[2]), int(p[3])))
    out = {}
    for k, iv in cov.items():
        iv.sort()
        tot, cs, ce = 0, -1, -1
        for s_, e in iv:
            if s_ > ce:
                tot += max(0, ce - cs)
                cs, ce = s_, e
            else:
                ce = max(ce, e)
        tot += max(0, ce - cs)
        if tot >= BP_MIN:
            out[k] = tot
    return out


def variants(setname, specs=None):
    """aif with smaller windows and larger candidate caps on the recall genomes, scored against minimap2 per 5 kb tile"""
    from hgt_analyze import lca_class
    G = genomes()
    vd = os.path.join(HD, "validate_" + setname)
    R = sorted({l[1:].split("|")[0] for l in open(os.path.join(vd, "R_windows.fa")) if l[0] == ">"})
    U = {l[1:].split("|")[0] for l in open(os.path.join(vd, "U.fa")) if l[0] == ">"}
    mm = mm2_calls(vd)
    db = os.path.join(ROOT, "data", "gtdb", "aif", "gtdb_0.88.aifd")
    aif = os.path.join(ROOT, "bin", "aif")
    res = json.load(open(os.path.join(ROOT, "results_hgt", setname, "validate_variants.json"))) \
        if os.path.exists(os.path.join(ROOT, "results_hgt", setname, "validate_variants.json")) else {}
    todo = [(5000, 32, 0), (5000, 256, 0), (2000, 32, 0), (1000, 32, 0), (1000, 64, 0), (2000, 64, 0)]
    if specs:  # e.g. w5000_c64_m99: window 5000, 64 candidates, --af-min-ident 0.99
        todo = []
        for sp in specs:
            d = {x[0]: int(x[1:]) for x in sp.split("_")}
            todo.append((d["w"], d["c"], d.get("m", 0)))
    for wsz, cap, mi in todo:
        key = f"w{wsz}_c{cap}" + (f"_m{mi}" if mi else "")
        if key in res:
            continue
        fa = os.path.join(vd, f"R_w{wsz}.fa")
        if not os.path.exists(fa):
            with open(fa + ".tmp", "w") as fo:
                for a in R:
                    for ci, s_ in enumerate(read_fasta_gz(G[a]["path"])):
                        for st in range(0, len(s_), wsz):
                            t = s_[st:st + wsz]
                            if len(t) >= min(WMIN, wsz // 2):
                                fo.write(f">{a}|{ci}|{st}|{len(t)}\n{t}\n")
            os.replace(fa + ".tmp", fa)
        out = os.path.join(vd, f"aif_{key}.tsv")
        t0 = time.time()
        extra = f"--af-min-ident {mi / 100:.2f} " if mi else ""
        subprocess.run(f"/usr/bin/time -f '%e %M' -o {out}.time {aif} db search {db} {fa} -t 48 -n {cap} --max-cand {cap} {extra}"
                       f"| /usr/bin/python3 -I {os.path.join(ROOT, 'bench_real', 'hgt_atlas.py')} filter > {out}",
                       shell=True, check=True)
        e, m = open(out + ".time").read().split()
        # aggregate sub-window calls to 5 kb tiles: aligned bases at >= 99% per (tile, partner)
        agg = collections.defaultdict(float)
        for l in open(out):
            w, ref, rank, ani, af, theta, lr = l.rstrip("\n").split("\t")
            a, ci, st, ln = w.split("|")
            if (mi or float(ani) >= ANI_MIN) and ref in U:
                tile = f"{a}|{ci}|{int(st) // W * W}"
                agg[(tile, ref)] += float(af) * int(ln)
        calls_ = {k for k, v in agg.items() if v >= BP_MIN}
        mmk = {(w.rsplit("|", 1)[0], r) for w, r in mm}
        both = calls_ & mmk
        xg = lambda S_: {k for k in S_ if lca_class(G[k[0].split("|")[0]]["tax"], G[k[1]]["tax"]) >= 1}
        res[key] = {"window": wsz, "cap": cap, "aif_s": float(e), "kb": int(m), "wall": time.time() - t0,
                    "aif_calls": len(calls_), "mm2_calls": len(mmk), "both": len(both),
                    "recall": len(both) / len(mmk), "agree": len(both) / max(1, len(calls_)),
                    "recall_cross_genus": len(xg(both)) / max(1, len(xg(mmk))), "aif_cross_genus": len(xg(calls_)),
                    "bases": os.path.getsize(fa)}
        print(key, res[key], flush=True)
        os.makedirs(os.path.join(ROOT, "results_hgt", setname), exist_ok=True)
        json.dump(res, open(os.path.join(ROOT, "results_hgt", setname, "validate_variants.json"), "w"), indent=1)


def blast_job(job):
    """BLAST (megablast) of one window against one partner genome: longest HSP at >= 99% identity, and query bases
    covered by HSPs at >= 99% identity (Smillie et al.: >= 99% over >= 500 bp)"""
    w, ref, qpath, rpath = job
    acc, ci, st, ln = w.split("|")
    seq = read_fasta_gz(qpath)[int(ci)][int(st):int(st) + int(ln)]
    td = tempfile.mkdtemp(dir=os.path.join(HD, "tmp"))
    try:
        q = os.path.join(td, "q.fa")
        open(q, "w").write(f">w\n{seq}\n")
        r = os.path.join(td, "r.fa")
        with gzip.open(rpath, "rb") as fi, open(r, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        o = subprocess.run(["/usr/bin/blastn", "-task", "megablast", "-query", q, "-subject", r, "-evalue", "1e-20",
                            "-outfmt", "6 pident length qstart qend", "-max_hsps", "1000"],
                           capture_output=True, check=True, text=True).stdout
        best, cov = 0, [0] * len(seq)
        for l in o.splitlines():
            pi, le, qs, qe = l.split("\t")
            if float(pi) >= 99.0:
                best = max(best, int(le))
                for i in range(min(int(qs), int(qe)) - 1, max(int(qs), int(qe))):
                    cov[i] = 1
        return w, ref, best, sum(cov)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def blast_check(setname):
    """BLAST adjudication of the precision sample and of the recall disagreements"""
    G = genomes()
    out = os.path.join(ROOT, "results_hgt", setname)
    os.makedirs(os.path.join(HD, "tmp"), exist_ok=True)
    sets = {"precision": [(r["window"], r["partner"]) for r in json.load(open(os.path.join(out, "validate_precision.json")))["rows"]]}
    S = json.load(open(os.path.join(out, "validate_recall.json")))
    sets["mm2_only"] = [tuple(x) for x in S["sample_mm2_only"]]
    sets["aif_only"] = [tuple(x) for x in S["sample_aif_only"]]
    res = {}
    for k, lst in sets.items():
        jobs = [(w, r, G[w.split("|")[0]]["path"], G[r]["path"]) for w, r in lst]
        with multiprocessing.get_context("fork").Pool(32) as pool:
            rr = list(pool.imap_unordered(blast_job, jobs))
        res[k] = {"n": len(rr), "hsp500": sum(1 for x in rr if x[2] >= BP_MIN), "cov500": sum(1 for x in rr if x[3] >= BP_MIN),
                  "rows": [{"window": x[0], "partner": x[1], "best_hsp99": x[2], "cov99": x[3]} for x in rr]}
        print(k, res[k]["n"], "single HSP >= 500 bp at >= 99%:", res[k]["hsp500"], " covered >= 500 bp:", res[k]["cov500"], flush=True)
    json.dump(res, open(os.path.join(out, "validate_blast.json"), "w"), indent=1)


def blast_chunk(job):
    q, db, out = job
    if not os.path.exists(out):
        subprocess.run(["/usr/bin/blastn", "-task", "megablast", "-query", q, "-db", db, "-evalue", "1e-20",
                        "-max_target_seqs", "5000", "-outfmt", "6 qseqid sseqid pident length qstart qend",
                        "-perc_identity", "95", "-num_threads", "1", "-out", out + ".tmp"], check=True)
        os.replace(out + ".tmp", out)
    return out


def blast_truth(setname, ngen=60):
    """BLAST (megablast) truth for the windows of ngen recall genomes against the universe U: (window, partner) with
    an HSP at >= 99% identity over >= 500 bp (Smillie et al.); then minimap2 and every aif variant scored against it"""
    from hgt_analyze import lca_class
    G = genomes()
    vd = os.path.join(HD, "validate_" + setname)
    bd = os.path.join(vd, "blast")
    os.makedirs(bd, exist_ok=True)
    R = sorted({l[1:].split("|")[0] for l in open(os.path.join(vd, "R_windows.fa")) if l[0] == ">"})
    random.seed(5)
    Rb = set(random.sample(R, ngen))
    db = os.path.join(bd, "U")
    if not os.path.exists(db + ".ndb") and not os.path.exists(db + ".nal"):
        with open(os.path.join(bd, "U.fa"), "w") as fo:
            for l in open(os.path.join(vd, "U.fa")):
                fo.write(l.replace("|", "__") if l[0] == ">" else l)
        subprocess.run(["/usr/bin/makeblastdb", "-in", os.path.join(bd, "U.fa"), "-dbtype", "nucl", "-out", db,
                        "-max_file_sz", "4GB"], check=True, capture_output=True)
    # query chunks: the 5 kb windows of the chosen genomes
    recs = []
    keep = False
    for l in open(os.path.join(vd, "R_windows.fa")):
        if l[0] == ">":
            keep = l[1:].split("|")[0] in Rb
            if keep:
                recs.append([l, ""])
        elif keep:
            recs[-1][1] += l
    nch = 48
    jobs = []
    for c in range(nch):
        qf = os.path.join(bd, f"q{c:02d}.fa")
        if not os.path.exists(qf):
            with open(qf, "w") as fo:
                for i in range(c, len(recs), nch):
                    fo.write(recs[i][0].replace("|", "__") + recs[i][1])
        jobs.append((qf, db, os.path.join(bd, f"q{c:02d}.tsv")))
    t0 = time.time()
    with multiprocessing.get_context("fork").Pool(nch) as pool:
        outs = pool.map(blast_chunk, jobs)
    tb = time.time() - t0
    truth = set()
    for o in outs:
        for l in open(o):
            q, sid, pi, le, qs, qe = l.rstrip("\n").split("\t")
            q = q.replace("__", "|")
            t = sid.split("__")[0]
            if t != q.split("|")[0] and float(pi) >= 99.0 and int(le) >= BP_MIN:
                truth.add((q.rsplit("|", 1)[0], t))
    U = {l[1:].split("|")[0] for l in open(os.path.join(vd, "U.fa")) if l[0] == ">"}
    xg = lambda S_: {k for k in S_ if lca_class(G[k[0].split("|")[0]]["tax"], G[k[1]]["tax"]) >= 1}
    def score(calls_):
        calls_ = {k for k in calls_ if k[0].split("|")[0] in Rb}
        both = calls_ & truth
        return {"calls": len(calls_), "both": len(both), "recall": len(both) / max(1, len(truth)),
                "precision": len(both) / max(1, len(calls_)), "recall_cross_genus": len(xg(both)) / max(1, len(xg(truth))),
                "precision_cross_genus": len(xg(both)) / max(1, len(xg(calls_)))}
    res = {"genomes": len(Rb), "windows": len(recs), "truth": len(truth), "truth_cross_genus": len(xg(truth)),
           "blast_wall_s": tb}
    mm = mm2_calls(vd)
    res["minimap2_summed"] = score({(w.rsplit("|", 1)[0], r) for w, r in mm})
    single = set()
    for l in open(os.path.join(vd, "mm2.paf")):
        p_ = l.split("\t")
        q, t = p_[0], p_[5].split("|")[0]
        if t != q.split("|")[0] and int(p_[10]) >= BP_MIN and int(p_[9]) / int(p_[10]) >= ANI_MIN:
            single.add((q.rsplit("|", 1)[0], t))
    res["minimap2_single"] = score(single)
    for f in sorted(os.listdir(vd)):
        if f.startswith("aif_w") and f.endswith(".tsv"):
            agg = collections.defaultdict(float)
            for l in open(os.path.join(vd, f)):
                w, ref, rank, ani, af, theta, lr = l.rstrip("\n").split("\t")
                a, ci, st, ln = w.split("|")
                if ("_m" in f or float(ani) >= ANI_MIN) and ref in U:
                    agg[(f"{a}|{ci}|{int(st) // W * W}", ref)] += float(af) * int(ln)
            res[f[4:-4]] = score({k for k, v in agg.items() if v >= BP_MIN})
    out = os.path.join(ROOT, "results_hgt", setname)
    json.dump(res, open(os.path.join(out, "validate_blast_truth.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "precision":
        precision(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 400)
    elif sys.argv[1] == "recall":
        recall(sys.argv[2])
    elif sys.argv[1] == "blast_truth":
        blast_truth(sys.argv[2])
    elif sys.argv[1] == "blast":
        blast_check(sys.argv[2])
    elif sys.argv[1] == "variants":
        variants(sys.argv[2], sys.argv[3:] or None)
    elif sys.argv[1] == "adjudicate":
        adjudicate(sys.argv[2])
    else:
        raise SystemExit(__doc__)
