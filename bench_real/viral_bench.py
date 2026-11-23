#!/usr/bin/env python3
"""Viral genome clustering (vOTUs) benchmark: aif vs skani vs vclust, against BLAST + CheckV-style ANI/AF.

Data sets (bench_real/viral_download.sh), 10,000 sequences each (seed 31):
  inphared  random complete phage genomes from INPHARED (GenBank, not dereplicated: pairs at every identity)
  gpd       random metagenomic viral contigs >= 10 kb of the Gut Phage Database (dereplicated at 95% ANI:
            mostly near-threshold rejections and partial overlaps)
Truth: all-vs-all blastn megablast (evalue 1e-5), then ANI/AF per pair as CheckV's anicalc.py: alignments
  sorted by bitscore, those overlapping already-kept ones on the query by > 10% dropped; ANI = alignment-length
  weighted identity; qcov / tcov = merged aligned fraction of query / target. Pairs with both directions keep
  the direction whose query is the shorter sequence, so AF = aligned fraction of the shorter sequence.
vOTUs: greedy centroid clustering as CheckV's aniclust.py (longest first; a sequence joins the first
  centroid with ANI >= 95% and AF(shorter) >= 85%), applied to every tool's pair table.
Tools (16 threads): aif db build -i (every sequence a reference) + db search all-vs-all (-n 100 --max-cand 100); skani triangle -i --sparse
  (default, and --small-genomes); vclust prefilter (--min-ident 0.8) + align.

  python3 bench_real/viral_bench.py prepare | blast | tools | evaluate [dataset ...]
  python3 bench_real/viral_bench.py prepare_tune; blast tune; tune   # choose --af-min-ident on a disjoint
                                         # INPHARED tuning set (3,000 genomes, seed 32)
"""
import csv
import gzip
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V = os.path.join(ROOT, "data", "viral")
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
SKANI = os.environ.get("SKANI", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani")
VCLUST = os.environ.get("VCLUST", os.path.join(ROOT, "data", "envs", "vclust", "bin", "vclust"))
TH = os.environ.get("THREADS", "16")
N = 10000
DATASETS = ["inphared", "gpd"]
TAUS = ["0.86", "0.88", "0.90", "0.92"]


def timed(cmd, out=None, err=None):
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=fo, stderr=fe, check=True)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def prepare():
    rng = random.Random(31)
    src = {"inphared": os.path.join(V, "inphared", "14Apr2025_genomes_excluding_refseq.fa.gz"),
           "gpd": os.path.join(V, "raw", "GPD_sequences.fa.gz")}
    for ds, path in src.items():
        names = [n for n, _ in fasta(path)]
        keep = set(rng.sample(names, N))
        d = os.path.join(V, ds)
        os.makedirs(d, exist_ok=True)
        L = {}
        with open(os.path.join(d, "seqs.fa"), "w") as fh:
            for n, s in fasta(path):
                if n in keep and n not in L:
                    s = s.upper()
                    fh.write(f">{n}\n{s}\n")
                    L[n] = len(s)
        with open(os.path.join(d, "lengths.tsv"), "w") as fh:
            for n, l in L.items():
                fh.write(f"{n}\t{l}\n")
        print(f"{ds}: {len(L)} sequences, {sum(L.values()) / 1e6:.1f} Mbp")


def prepare_tune():
    """3,000 INPHARED genomes disjoint from the test sample (seed 32): tuning set for --af-min-ident."""
    test = set(lengths("inphared"))
    path = os.path.join(V, "inphared", "14Apr2025_genomes_excluding_refseq.fa.gz")
    names = [n for n, _ in fasta(path) if n not in test]
    keep = set(random.Random(32).sample(sorted(set(names)), 3000))
    d = os.path.join(V, "tune")
    os.makedirs(d, exist_ok=True)
    L = {}
    with open(os.path.join(d, "seqs.fa"), "w") as fh:
        for n, s in fasta(path):
            if n in keep and n not in L:
                fh.write(f">{n}\n{s.upper()}\n")
                L[n] = len(s)
    with open(os.path.join(d, "lengths.tsv"), "w") as fh:
        for n, l in L.items():
            fh.write(f"{n}\t{l}\n")
    print(f"tune: {len(L)} sequences, {sum(L.values()) / 1e6:.1f} Mbp")


def tune():
    """aif with --af-min-ident in TAUS on the tuning set; choose the value with the best vOTU-edge F1."""
    d = os.path.join(V, "tune")
    fa = os.path.join(d, "seqs.fa")
    db = os.path.join(d, "aif.aifd")
    if not os.path.exists(db):
        subprocess.run([AIF, "db", "build", "-i", "-o", db, "-k", "20", "-t", TH, fa], stderr=subprocess.DEVNULL, check=True)
    L = lengths("tune")
    T = {(r["short"], r["long"]): (float(r["ani"]), float(r["af_short"])) for r in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t")}
    te = {k for k, (a, f) in T.items() if a >= 0.95 and f >= 0.85}
    res = {}
    for tau in ["0"] + TAUS:
        o = os.path.join(d, f"aif_t{tau}.tsv")
        subprocess.run([AIF, "db", "search", db, fa, "-t", TH, "-n", "100", "--max-cand", "100", "--af-min-ident", tau],
                       stdout=open(o, "w"), stderr=subprocess.DEVNULL, check=True)
        P = load_aif_file(o, L)
        near = [k for k, (a, f) in T.items() if (0.93 <= a < 0.97 and f >= 0.5) or (a >= 0.93 and 0.75 <= f < 0.95)]
        rep = [k for k in near if k in P]
        allp = [k for k, (a, f) in T.items() if a >= 0.8 and f >= 0.1 and k in P]
        edges = {k for k, (a, f) in P.items() if a >= 0.95 and f >= 0.85}
        tp = len(edges & te)
        res[tau] = {"af_mae_near_pp": 100 * sum(abs(P[k][1] - T[k][1]) for k in rep) / max(1, len(rep)),
                    "af_mae_all_pp": 100 * sum(abs(P[k][1] - T[k][1]) for k in allp) / max(1, len(allp)),
                    "edge_precision": tp / max(1, len(edges)), "edge_recall": tp / max(1, len(te)),
                    "edge_f1": 2 * tp / max(1, len(edges) + len(te))}
        print(tau, {k: round(v, 4) for k, v in res[tau].items()}, flush=True)
    # rule (fixed on the tuning set): highest edge F1; candidates within 0.005 of it are tied and the
    # tie is broken by the near-threshold AF error
    f1max = max(res[t]["edge_f1"] for t in TAUS)
    best = min([t for t in TAUS if res[t]["edge_f1"] >= f1max - 0.005], key=lambda t: res[t]["af_mae_near_pp"])
    res["chosen"] = best
    res["rule"] = "highest edge F1; ties within 0.005 broken by near-threshold AF MAE"
    json.dump(res, open(os.path.join(V, "tune_af_min_ident.json"), "w"), indent=1)
    print("chosen --af-min-ident", best)


def load_aif_file(path, L):
    P = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "ok" or r["reference"] == "*":
            continue
        q, t = r["contig"].split()[0], r["reference"]
        if q == t or (L[q], q) > (L[t], t):
            continue
        P[(q, t)] = (float(r["ani"]), float(r["af"]))
    return P


def lengths(ds):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(os.path.join(V, ds, "lengths.tsv"))}


def blast(datasets):
    out = {}
    for ds in datasets:
        d = os.path.join(V, ds)
        fa = os.path.join(d, "seqs.fa")
        if os.path.exists(os.path.join(d, "blast.done")):
            anicalc(ds)
            continue
        t1, m1 = timed(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", os.path.join(d, "blastdb")])
        t2, m2 = timed(["blastn", "-task", "megablast", "-query", fa, "-db", os.path.join(d, "blastdb"), "-evalue", "1e-5",
                        "-max_target_seqs", "25000", "-num_threads", os.environ.get("BLAST_THREADS", "56"), "-mt_mode", "1",
                        "-outfmt", "6 std qlen slen", "-out", os.path.join(d, "blast.tsv")])
        out[ds] = {"makeblastdb_s": t1, "blastn_s": t2, "blastn_kb": m2, "threads": int(os.environ.get("BLAST_THREADS", "56"))}
        print(ds, out[ds], flush=True)
        p = os.path.join(V, "results_blast_time.json")
        old = json.load(open(p)) if os.path.exists(p) else {}
        old.update({ds: out[ds]})
        json.dump(old, open(p, "w"), indent=1)
        open(os.path.join(d, "blast.done"), "w").write("ok\n")
        anicalc(ds)


def merged_len(iv):
    iv.sort()
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce + 1:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs + 1 if cs is not None else 0)


def anicalc(ds):
    """CheckV anicalc.py logic on blast.tsv -> truth.tsv (qname tname ani qcov tcov), shorter-as-query direction."""
    d = os.path.join(V, ds)
    L = lengths(ds)
    pairs = defaultdict(list)
    for line in open(os.path.join(d, "blast.tsv")):
        p = line.rstrip("\n").split("\t")
        q, t = p[0], p[1]
        if q == t:
            continue
        pairs[(q, t)].append((float(p[11]), int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[3]), float(p[2])))
    best = {}
    for (q, t), alns in pairs.items():
        alns.sort(key=lambda a: -a[0])
        cov = bytearray(L[q] + 2)
        keep = []
        for bs, qs, qe, ts, te, ln, pid in alns:
            a, b = min(qs, qe), max(qs, qe)
            covered = cov[a:b + 1].count(1)
            if covered / (b - a + 1) <= 0.10:
                cov[a:b + 1] = b"\x01" * (b - a + 1)
                keep.append((qs, qe, ts, te, ln, pid))
        ani = sum(k[4] * k[5] for k in keep) / sum(k[4] for k in keep)
        qcov = merged_len([(min(k[0], k[1]), max(k[0], k[1])) for k in keep]) / L[q]
        tcov = merged_len([(min(k[2], k[3]), max(k[2], k[3])) for k in keep]) / L[t]
        key = tuple(sorted((q, t)))
        short = q if (L[q], q) < (L[t], t) else t
        rec = (q, t, ani / 100, min(1.0, qcov), min(1.0, tcov))
        # keep the direction whose query is the shorter sequence (AF = its aligned fraction)
        if key not in best or (q == short and best[key][0] != short):
            best[key] = rec
    with open(os.path.join(d, "truth.tsv"), "w") as fh:
        fh.write("short\tlong\tani\taf_short\taf_long\n")
        for key, (q, t, ani, qc, tc) in best.items():
            s, l, afs, afl = (q, t, qc, tc) if (L[q], q) <= (L[t], t) else (t, q, tc, qc)
            fh.write(f"{s}\t{l}\t{ani:.5f}\t{afs:.5f}\t{afl:.5f}\n")
    print(ds, "truth pairs:", len(best), flush=True)


def tools(datasets):
    tp = os.path.join(V, "results_tool_time.json")
    out = json.load(open(tp)) if os.path.exists(tp) else {}
    for ds in datasets:
        d = os.path.join(V, ds)
        fa = os.path.join(d, "seqs.fa")
        r = {}
        db = os.path.join(d, "aif.aifd")
        r["aif_build"] = timed([AIF, "db", "build", "-i", "-o", db, "-k", "20", "-t", TH, fa], err=db + ".log")
        r["aif_search"] = timed([AIF, "db", "search", db, fa, "-t", TH, "-n", "100", "--max-cand", "100"],
                                out=os.path.join(d, "aif.tsv"), err=os.path.join(d, "aif.err"))
        r["aif_db_bytes"] = os.path.getsize(db)
        tau = json.load(open(os.path.join(V, "tune_af_min_ident.json")))["chosen"]
        r["aif_t_search"] = timed([AIF, "db", "search", db, fa, "-t", TH, "-n", "100", "--max-cand", "100", "--af-min-ident", tau],
                                  out=os.path.join(d, "aif_t.tsv"), err=os.path.join(d, "aif_t.err"))
        r["af_min_ident"] = tau
        for name, extra in (("skani", []), ("skani_small", ["--small-genomes"])):
            r[name] = timed([SKANI, "triangle", "-i", fa, "--sparse", "-t", TH, "-o", os.path.join(d, name + ".tsv")] + extra,
                            err=os.path.join(d, name + ".log"))
        fl = os.path.join(d, "vclust_fltr.txt")
        r["vclust_prefilter"] = timed([VCLUST, "prefilter", "-i", fa, "-o", fl, "--min-ident", "0.8", "-t", TH],
                                      err=os.path.join(d, "vclust.log"))
        r["vclust_align"] = timed([VCLUST, "align", "-i", fa, "-o", os.path.join(d, "vclust_ani.tsv"), "--filter", fl, "-t", TH],
                                  err=os.path.join(d, "vclust_align.log"))
        out[ds] = r
        print(ds, r, flush=True)
        json.dump(out, open(tp, "w"), indent=1)


def load_tool(ds, tool, L):
    """{(short, long): (ani, af_short)}"""
    d = os.path.join(V, ds)
    P = {}

    def put(q, t, ani, afq, aft):
        if q == t:
            return
        s, l, afs = (q, t, afq) if (L[q], q) <= (L[t], t) else (t, q, aft)
        key = (s, l)
        if afs is None:
            return
        if key not in P or P[key][1] < afs:
            P[key] = (ani, afs)
    if tool in ("aif", "aif_t"):
        # aif's AF is the query's aligned fraction: use rows whose query is the shorter sequence
        return load_aif_file(os.path.join(d, tool + ".tsv"), L)
    elif tool.startswith("skani"):
        for line in list(open(os.path.join(d, tool + ".tsv")))[1:]:
            p = line.rstrip("\n").split("\t")
            # Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
            rn, qn = p[5].split()[0], p[6].split()[0]
            put(qn, rn, float(p[2]) / 100, float(p[4]) / 100, float(p[3]) / 100)
    elif tool == "vclust":
        ids = {}
        for line in list(open(os.path.join(d, "vclust_ani.ids.tsv")))[1:]:
            p = line.rstrip("\n").split("\t")
            ids[p[0]] = p[1]
        for r in csv.DictReader(open(os.path.join(d, "vclust_ani.tsv")), delimiter="\t"):
            q, t = r["query"], r["reference"]
            ani = float(r["ani"])
            put(q, t, ani, float(r["qcov"]), float(r["rcov"]))
    return P


def greedy_clusters(names, L, edges):
    """CheckV aniclust: longest first; each sequence joins the first (longest) centroid it has an edge to."""
    adj = defaultdict(set)
    for s, l in edges:
        adj[s].add(l)
        adj[l].add(s)
    order = sorted(names, key=lambda n: (-L[n], n))
    cl = {}
    for n in order:
        if n in cl:
            continue
        cl[n] = n
        for m in adj[n]:
            if m not in cl:
                cl[m] = n
    return cl


def ari(a, b):
    from math import comb
    names = list(a)
    cont = defaultdict(int)
    for n in names:
        cont[(a[n], b[n])] += 1
    sa, sb = defaultdict(int), defaultdict(int)
    for (x, y), c in cont.items():
        sa[x] += c
        sb[y] += c
    idx = sum(comb(c, 2) for c in cont.values())
    ea = sum(comb(c, 2) for c in sa.values())
    eb = sum(comb(c, 2) for c in sb.values())
    tot = comb(len(names), 2)
    exp = ea * eb / tot
    mx = (ea + eb) / 2
    return (idx - exp) / (mx - exp) if mx != exp else 1.0


def evaluate(datasets):
    res = {}
    TT = json.load(open(os.path.join(V, "results_tool_time.json"))) if os.path.exists(os.path.join(V, "results_tool_time.json")) else {}
    BT = json.load(open(os.path.join(V, "results_blast_time.json"))) if os.path.exists(os.path.join(V, "results_blast_time.json")) else {}
    ANIB = [(0.8, 0.9, "80-90"), (0.9, 0.95, "90-95"), (0.95, 0.98, "95-98"), (0.98, 1.01, ">=98")]
    for ds in datasets:
        L = lengths(ds)
        T = {}
        for r in csv.DictReader(open(os.path.join(V, ds, "truth.tsv")), delimiter="\t"):
            T[(r["short"], r["long"])] = (float(r["ani"]), float(r["af_short"]))
        true_edges = {k for k, (a, f) in T.items() if a >= 0.95 and f >= 0.85}
        truth_cl = greedy_clusters(list(L), L, true_edges)
        res[ds] = {"n_seq": len(L), "mbp": sum(L.values()) / 1e6, "truth_pairs": len(T), "truth_edges": len(true_edges),
                   "truth_votus": len(set(truth_cl.values())), "tools": {}}
        for tool in ("aif", "aif_t", "skani", "skani_small", "vclust"):
            try:
                P = load_tool(ds, tool, L)
            except FileNotFoundError:
                continue
            ev = {}
            # pair accuracy on truth pairs with ANI >= 80% and AF(shorter) >= 0.1
            for lo, hi, lab in ANIB + [(0.8, 1.01, "all")]:
                ks = [k for k, (a, f) in T.items() if lo <= a < hi and f >= 0.1]
                rep = [k for k in ks if k in P]
                ea = [abs(P[k][0] - T[k][0]) for k in rep]
                ef = [abs(P[k][1] - T[k][1]) for k in rep]
                ev[lab] = {"n": len(ks), "reported": len(rep) / max(1, len(ks)),
                           "ani_mae_pp": 100 * sum(ea) / max(1, len(ea)), "af_mae_pp": 100 * sum(ef) / max(1, len(ef))}
            # near the vOTU thresholds: truth ANI 93-97% or AF 75-95%
            near = [k for k, (a, f) in T.items() if (0.93 <= a < 0.97 and f >= 0.5) or (a >= 0.93 and 0.75 <= f < 0.95)]
            rep = [k for k in near if k in P]
            ev["near_threshold"] = {"n": len(near), "reported": len(rep) / max(1, len(near)),
                                    "ani_mae_pp": 100 * sum(abs(P[k][0] - T[k][0]) for k in rep) / max(1, len(rep)),
                                    "af_mae_pp": 100 * sum(abs(P[k][1] - T[k][1]) for k in rep) / max(1, len(rep))}
            edges = {k for k, (a, f) in P.items() if a >= 0.95 and f >= 0.85}
            tp = len(edges & true_edges)
            ev["edges"] = {"predicted": len(edges), "tp": tp, "precision": tp / max(1, len(edges)), "recall": tp / max(1, len(true_edges)),
                           "f1": 2 * tp / max(1, len(edges) + len(true_edges))}
            cl = greedy_clusters(list(L), L, edges)
            ev["clusters"] = {"n_votus": len(set(cl.values())), "ari": ari(truth_cl, cl),
                              "identical_cluster_frac": None}
            # fraction of truth vOTUs (size >= 2) reproduced exactly
            tg = defaultdict(set)
            pg = defaultdict(set)
            for n in L:
                tg[truth_cl[n]].add(n)
                pg[cl[n]].add(n)
            pset = {frozenset(v) for v in pg.values()}
            big = [frozenset(v) for v in tg.values() if len(v) >= 2]
            ev["clusters"]["identical_cluster_frac"] = sum(1 for v in big if v in pset) / max(1, len(big))
            ev["clusters"]["n_truth_multi"] = len(big)
            res[ds]["tools"][tool] = ev
        res[ds]["time"] = {"tools": TT.get(ds), "blast": BT.get(ds)}
    p = os.path.join(V, "viral_eval.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    old.update(res)
    json.dump(old, open(p, "w"), indent=1)
    for ds, r in res.items():
        print(f"== {ds}: {r['n_seq']} seqs, {r['mbp']:.0f} Mbp, truth pairs {r['truth_pairs']}, edges {r['truth_edges']}, vOTUs {r['truth_votus']}")
        for tool, ev in r["tools"].items():
            a, e, c, nt = ev["all"], ev["edges"], ev["clusters"], ev["near_threshold"]
            print(f"  {tool:12s} pairs rep {a['reported']:.3f} ANI MAE {a['ani_mae_pp']:.3f} AF MAE {a['af_mae_pp']:.2f} | near-thr rep {nt['reported']:.3f} "
                  f"ANI {nt['ani_mae_pp']:.3f} AF {nt['af_mae_pp']:.2f} | edges P {e['precision']:.3f} R {e['recall']:.3f} F1 {e['f1']:.3f} | "
                  f"vOTUs {c['n_votus']} ARI {c['ari']:.4f} identical {c['identical_cluster_frac']:.3f}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    ds = sys.argv[2:] or DATASETS
    {"prepare": lambda: prepare(), "prepare_tune": lambda: prepare_tune(), "tune": lambda: tune(), "blast": lambda: blast(ds),
     "tools": lambda: tools(ds), "evaluate": lambda: evaluate(ds)}[cmd]()
