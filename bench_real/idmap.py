#!/usr/bin/env python3
"""Identity maps along contigs: aif --segments vs fastANI --visualize, against MUMmer alignments.

Real sets (within-species): every query genome against its reference.
  truth    nucmer --maxmatch + delta-filter -q (as bench_real/truth.py): aligned blocks per contig;
           per base: aligned or not; identity of a window = 1 - (substitutions + indel bases in the
           window, show-snps) / aligned bases in the window (block-level %IDY would average identity
           over whole alignment blocks, often tens of kb).
  aif      aif query --segments (1 kb windows, and 3 kb windows for a like-for-like comparison):
           aligned/unaligned segments tiling each contig; per aligned window: ANI and 95% interval.
  fastANI  --visualize: 3 kb fragments with identity (mapped fragments only; global coordinates are
           mapped back to contigs by cumulative offsets).
Metrics
  identity windows with truth AF >= 0.9: MAE vs truth identity, Pearson r, CI coverage (aif),
           by true identity class;
  islands  per base, unaligned (truth: not covered by any block) vs aif unaligned segments and
           fastANI's unmapped 3 kb fragments: precision / recall of unaligned bases, on contigs the
           tool reports;
  breakpoints true aligned/unaligned boundaries (unaligned stretches >= 300 bp, >= 300 bp from contig
           ends): distance to the nearest predicted boundary.
  python3 bench_real/idmap.py run        # truth coordinates, aif segments, fastANI
  REAL_DIR=data/real_tune python3 bench_real/idmap.py calibrate   # window-interval multipliers from the
                                         # tuning set (split conformal, 95%) -> data/real/results/idmap_ci.json
  python3 bench_real/idmap.py evaluate   # -> data/real/results/idmap.json (statistical and calibrated coverage)
"""
import bisect
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("REAL_DIR", os.path.join(ROOT, "data", "real"))
OUT = os.path.join(R, "idmap")
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
FASTANI = os.environ.get("FASTANI", "fastANI")
MUMMER = os.environ.get("MUMMER", os.path.join(ROOT, "data", "envs", "mummer", "bin"))
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]


def genomes(s):
    rows = [l.rstrip("\n").split("\t") for l in list(open(os.path.join(R, "sets", s, "genomes.tsv")))[1:]]
    return [(r[0], r[4]) for r in rows]


def contig_lengths(fa):
    L, order, name = {}, [], None
    for line in open(fa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            L[name] = 0
            order.append(name)
        else:
            L[name] += len(line.strip())
    return L, order


def one_genome(job):
    s, gid, ref = job
    import truth as T
    T.MUMMER = MUMMER
    q = os.path.join(R, "sets", s, "genomes", gid + ".fa")
    rfa = os.path.join(R, "refs", ref + ".fa")
    with tempfile.TemporaryDirectory() as td:
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre, rfa, q], check=True, capture_output=True)
        cq = T.coords(pre + ".delta", "-q")
        fq = subprocess.run([os.path.join(MUMMER, "delta-filter"), "-q", pre + ".delta"], capture_output=True, check=True, text=True).stdout
        open(pre + ".q.delta", "w").write(fq)
        snps = subprocess.run([os.path.join(MUMMER, "show-snps"), "-T", "-H", "-q", pre + ".q.delta"], capture_output=True, check=True,
                              text=True).stdout
        fo = os.path.join(td, "fa")
        subprocess.run([FASTANI, "-q", q, "-r", rfa, "--visualize", "-o", fo], check=True, capture_output=True)
        vis = open(fo + ".visual").read() if os.path.exists(fo + ".visual") else ""
    with open(os.path.join(OUT, s, gid + ".coords.tsv"), "w") as fh:
        for name, a, b, lq, idy, _ in cq:
            fh.write(f"{name}\t{a}\t{b}\t{idy}\n")
    with open(os.path.join(OUT, s, gid + ".snps.tsv"), "w") as fh:  # contig, 0-based position of a difference
        for line in snps.splitlines():
            p = line.split("\t")
            if len(p) >= 12:
                fh.write(f"{p[-1]}\t{int(p[3]) - 1}\n")
    L, order = contig_lengths(q)
    off, acc = [], 0
    for c in order:
        off.append(acc)
        acc += L[c]
    with open(os.path.join(OUT, s, gid + ".fastani.tsv"), "w") as fh:
        for line in vis.splitlines():
            p = line.split("\t")
            if len(p) < 9:
                continue
            g0, g1, idy = int(p[6]), int(p[7]), float(p[2])
            i = bisect.bisect_right(off, g0) - 1
            fh.write(f"{order[i]}\t{g0 - off[i]}\t{g1 - off[i] + 1}\t{idy}\n")
    return s, gid


def run():
    jobs = []
    for s in SETS:
        os.makedirs(os.path.join(OUT, s), exist_ok=True)
        jobs += [(s, g, ref) for g, ref in genomes(s) if not os.path.exists(os.path.join(OUT, s, g + ".snps.tsv"))]
    with Pool(int(os.environ.get("JOBS", "32"))) as pool:
        for k, _ in enumerate(pool.imap_unordered(one_genome, jobs, chunksize=1)):
            if (k + 1) % 50 == 0:
                print(f"{k + 1}/{len(jobs)}", flush=True)
    for s in SETS:
        ref = genomes(s)[0][1]
        sk = os.path.join(R, "results", "aif", f"{s}.{ref}.aif")
        q = os.path.join(R, "sets", s, f"{ref}.contigs.fa")
        for w in (1000, 3000):
            subprocess.run([AIF, "query", sk, q, "-t", "16", "--segments", os.path.join(OUT, s, f"aif_seg{w}.tsv"), "--seg-window", str(w), "--seg-ci-stat"],
                           stdout=subprocess.DEVNULL, check=True)
        print("aif segments", s, flush=True)


def load_diffs(s):
    """{contig: sorted 0-based positions of alignment differences (substitutions, indel bases)}"""
    d = {}
    for g, _ in genomes(s):
        for line in open(os.path.join(OUT, s, g + ".snps.tsv")):
            c, pos = line.split("\t")
            d.setdefault(c, []).append(int(pos))
    for v in d.values():
        v.sort()
    return d


def est_class(a):
    return "ge99" if a >= 0.99 else "97-99" if a >= 0.97 else "93-97" if a >= 0.93 else "lt93"


def window_rows(w):
    """(set, contig, start, end, ani, lo, hi) of aif's aligned windows of size w, with truth (AF, identity)"""
    out = []
    for s in SETS:
        blocks, cov = {}, {}
        for g, _ in genomes(s):
            for line in open(os.path.join(OUT, s, g + ".coords.tsv")):
                c, a, b, idy = line.split("\t")
                blocks.setdefault(c, []).append((int(a) - 1, int(b), float(idy) / 100))
        cov = {c: merge([(a, b) for a, b, _ in v]) for c, v in blocks.items()}
        diffs = load_diffs(s)
        for line in list(open(os.path.join(OUT, s, f"aif_seg{w}.tsv")))[1:]:
            p = line.rstrip("\n").split("\t")
            if p[3] != "1" or p[4] == "nan":
                continue
            c, a, b = p[0], int(p[1]), int(p[2])
            alb = overlap(cov.get(c, []), a, b)
            al = alb / max(1, b - a)
            d = diffs.get(c, [])
            nd = bisect.bisect_left(d, b) - bisect.bisect_left(d, a)
            tid = 1 - nd / alb if alb > 0 else float("nan")
            out.append((s, c, a, b, float(p[4]), float(p[5]), float(p[6]), al, tid))
    return out


def calibrate():
    import math
    res = {}
    for w in (1000, 3000):
        sc = {}
        for s, c, a, b, ani, lo, hi, al, tid in window_rows(w):
            if al < 0.9 or tid != tid:
                continue
            se = max((hi - lo) / (2 * 1.959964), 1e-4)
            sc.setdefault(est_class(ani), []).append(abs(ani - tid) / se)
        res[str(w)] = {}
        for k, v in sc.items():
            v.sort()
            n = len(v)
            res[str(w)][k] = {"n": n, "q95": v[min(n, math.ceil((n + 1) * 0.95)) - 1]}
    json.dump(res, open(os.path.join(ROOT, "data", "real", "results", "idmap_ci.json"), "w"), indent=1)
    print(res)


def merge(iv):
    iv = sorted(iv)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def overlap(iv, a, b):
    """total overlap of merged intervals iv with [a, b)"""
    t = 0
    for x, y in iv:
        if y <= a:
            continue
        if x >= b:
            break
        t += min(b, y) - max(a, x)
    return t


def evaluate():
    res = {}
    IDB = [(0, 0.9, "<90"), (0.9, 0.95, "90-95"), (0.95, 0.98, "95-98"), (0.98, 0.995, "98-99.5"), (0.995, 1.01, ">=99.5")]
    agg = {}
    for s in SETS:
        L = {}
        blocks = {}  # contig -> [(a, b, idy)] (0-based half-open)
        fa = {}      # contig -> [(a, b, idy)]
        for g, _ in genomes(s):
            Lg, _ = contig_lengths(os.path.join(R, "sets", s, "genomes", g + ".fa"))
            L.update(Lg)
            for line in open(os.path.join(OUT, s, g + ".coords.tsv")):
                c, a, b, idy = line.split("\t")
                blocks.setdefault(c, []).append((int(a) - 1, int(b), float(idy) / 100))
            for line in open(os.path.join(OUT, s, g + ".fastani.tsv")):
                c, a, b, idy = line.split("\t")
                fa.setdefault(c, []).append((int(a), int(b), float(idy) / 100))
        cov = {c: merge([(a, b) for a, b, _ in v]) for c, v in blocks.items()}
        diffs = load_diffs(s)

        def true_window(c, a, b):
            alb = overlap(cov.get(c, []), a, b)
            al = alb / max(1, b - a)
            d = diffs.get(c, [])
            nd = bisect.bisect_left(d, b) - bisect.bisect_left(d, a)
            return al, (1 - nd / alb if alb > 0 else float("nan"))

        def add(tool, key, est, tru, lo=None, hi=None):
            d = agg.setdefault((tool, key), {"e": [], "x": [], "y": [], "cov": []})
            d["e"].append(abs(est - tru))
            d["x"].append(est)
            d["y"].append(tru)
            if lo is not None:
                d["cov"].append(lo - 1e-9 <= tru <= hi + 1e-9)

        def idclass(v):
            for lo_, hi_, lab in IDB:
                if lo_ <= v < hi_:
                    return lab
            return IDB[-1][2]

        # aif windows
        cal = {}
        cp = os.path.join(ROOT, "data", "real", "results", "idmap_ci.json")
        if os.path.exists(cp):
            cal = json.load(open(cp))
        seg = {}
        for w in (1000, 3000):
            for line in list(open(os.path.join(OUT, s, f"aif_seg{w}.tsv")))[1:]:
                p = line.rstrip("\n").split("\t")
                c, a, b, al = p[0], int(p[1]), int(p[2]), int(p[3])
                if w == 1000:
                    seg.setdefault(c, []).append((a, b, al))
                if not al or p[4] == "nan":
                    continue
                tal, tid = true_window(c, a, b)
                if tal >= 0.9 and tid == tid:
                    ani, lo, hi = float(p[4]), float(p[5]), float(p[6])
                    for key in (s, "all", "all|" + idclass(tid)):
                        add(f"aif_{w}", key, ani, tid, lo, hi)
                    q = cal.get(str(w), {}).get(est_class(ani), {}).get("q95")
                    if q is not None:
                        se = max((hi - lo) / (2 * 1.959964), 1e-4)
                        for key in (s, "all", "all|" + idclass(tid)):
                            add(f"aif_{w}_cal", key, ani, tid, ani - q * se, ani + q * se)
        # fastANI fragments
        for c, v in fa.items():
            for a, b, idy in v:
                tal, tid = true_window(c, a, b)
                if tal >= 0.9 and tid == tid:
                    for key in (s, "all", "all|" + idclass(tid)):
                        add("fastani_3000", key, idy, tid)
        # islands: per-base unaligned calls on contigs each tool reports
        for tool in ("aif", "fastani"):
            tp = fp = fn = 0
            for c, Lc in L.items():
                if tool == "aif":
                    if c not in seg:
                        continue
                    pred_un = merge([(a, b) for a, b, al in seg[c] if not al])
                else:
                    if c not in fa:
                        continue
                    # fastANI: 3 kb fragments of the contig; unmapped fragments are unaligned
                    nfr = Lc // 3000
                    mapped = set(a // 3000 for a, _, _ in fa[c])
                    pred_un = merge([(j * 3000, (j + 1) * 3000) for j in range(nfr) if j not in mapped])
                    Lc = nfr * 3000  # bases beyond the last full fragment are not assessed
                true_un = []
                last = 0
                for a, b in cov.get(c, []):
                    if a > last:
                        true_un.append((last, min(a, Lc)))
                    last = max(last, b)
                if last < Lc:
                    true_un.append((last, Lc))
                true_un = [(a, b) for a, b in true_un if b > a]
                inter = sum(overlap(pred_un, a, b) for a, b in true_un)
                tp += inter
                fp += sum(b - a for a, b in pred_un if b <= Lc) - inter
                fn += sum(b - a for a, b in true_un) - inter
            d = agg.setdefault((tool, "islands|" + s), {})
            d.update(tp=tp, fp=fp, fn=fn)
            g = agg.setdefault((tool, "islands|all"), {"tp": 0, "fp": 0, "fn": 0})
            g["tp"] += tp
            g["fp"] += fp
            g["fn"] += fn
        # breakpoints (aif only; fastANI resolution is 3 kb fragments)
        dist = agg.setdefault(("aif", "bp|all"), {"d": []})
        dists = agg.setdefault(("aif", "bp|" + s), {"d": []})
        for c, Lc in L.items():
            if c not in seg:
                continue
            cv = cov.get(c, [])
            tb = []
            for (a0, b0), (a1, b1) in zip(cv, cv[1:]):
                if a1 - b0 >= 300:
                    tb += [b0, a1]
            if cv and cv[0][0] >= 300:
                tb.append(cv[0][0])
            if cv and Lc - cv[-1][1] >= 300:
                tb.append(cv[-1][1])
            tb = [x for x in tb if 300 <= x <= Lc - 300]
            pb = sorted(set(b for a, b, al in seg[c][:-1]))
            for x in tb:
                if not pb:
                    dd = float("inf")
                else:
                    i = bisect.bisect_left(pb, x)
                    dd = min(abs(pb[j] - x) for j in (i - 1, i) if 0 <= j < len(pb))
                dist["d"].append(dd)
                dists["d"].append(dd)
    for (tool, key), d in agg.items():
        if "e" in d:
            n = len(d["e"])
            mx, my = sum(d["x"]) / n, sum(d["y"]) / n
            sxy = sum((x - mx) * (y - my) for x, y in zip(d["x"], d["y"]))
            sxx = sum((x - mx) ** 2 for x in d["x"])
            syy = sum((y - my) ** 2 for y in d["y"])
            res.setdefault(tool, {})[key] = {"n": n, "mae_pp": 100 * sum(d["e"]) / n,
                                             "r": sxy / (sxx * syy) ** 0.5 if sxx > 0 and syy > 0 else float("nan"),
                                             "ci_cov": sum(d["cov"]) / len(d["cov"]) if d["cov"] else None}
        elif "tp" in d:
            tp, fp, fn = d["tp"], d["fp"], d["fn"]
            res.setdefault(tool, {})[key] = {"precision": tp / max(1, tp + fp), "recall": tp / max(1, tp + fn), "true_bp": tp + fn}
        else:
            v = sorted(d["d"])
            n = len(v)
            res.setdefault(tool, {})[key] = {"n": n, "median_bp": v[n // 2] if n else None,
                                             "within_100": sum(x <= 100 for x in v) / max(1, n), "within_250": sum(x <= 250 for x in v) / max(1, n),
                                             "within_500": sum(x <= 500 for x in v) / max(1, n)}
    json.dump(res, open(os.path.join(R, "results", "idmap.json"), "w"), indent=1)
    for tool in sorted(res):
        for key in sorted(res[tool]):
            print(tool, key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[tool][key].items()})


if __name__ == "__main__":
    {"run": run, "evaluate": evaluate, "calibrate": calibrate}[sys.argv[1]]()
