#!/usr/bin/env python3
"""AF bias of aif on UHGG contigs: diagnosis and tuning on a set disjoint from the benchmark sample.

The benchmark (uhgg_compare.py) found aif's top-hit AF biased high (+5 pp against MUMmer delta-filter -q on S3).
Every choice here is made on the tuning set T; S1-S3 are scored once, after.

  T         the next contigs of every stage-2 cell's seeded order after S3's (same design, disjoint from S1-S3);
            cells that S3 exhausts are absent
  variants  aif db search with extra options, on T (runs/T_<name>/aif.tsv; the top hit's segments in aif.seg.tsv)
  truth     MUMmer for every variant's top-5 hits, appended to truth/T.pairs.tsv (uhgg_compare.truth)
  diag      base-level comparison of aif's aligned segments (top hit) with MUMmer's -q alignments: aif-only bases
            split into edge extensions, bridged gaps and islands; MUMmer-only bases; aif's window identity there

  python3 -I bench_real/uhgg_af.py tune [3000]
  python3 -I bench_real/uhgg_af.py run <name> [aif db search options...]
  python3 -I bench_real/uhgg_af.py truth <name>
  python3 -I bench_real/uhgg_af.py score <name> [<name> ...]
  python3 -I bench_real/uhgg_af.py diag <name>
"""
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uhgg_compare as U  # noqa: E402

C = U.C
SD = os.path.join(C, "sample")
NT = 3000


def tune(n=NT):
    cells, rows, X = U.stage2_cells()
    avail = {k: len(v) for k, v in cells.items()}
    a3 = U.waterfill(avail, U.SIZES["S3"])
    # the refactored cells must reproduce the benchmark's S3 exactly
    s3 = {r["contig"] for r in csv.DictReader(open(os.path.join(SD, "sample.tsv")), delimiter="\t") if r["S3"] == "1"}
    mine = {c for k, v in cells.items() for c in v[: a3[k]]}
    assert mine == s3, (len(mine), len(s3), len(mine ^ s3))
    rem = {k: avail[k] - a3[k] for k in avail if avail[k] > a3[k]}
    at = U.waterfill(rem, n)
    T = {c: k for k in rem for c in cells[k][a3[k]: a3[k] + at[k]]}
    assert not (set(T) & s3)
    with open(os.path.join(SD, "T.tsv"), "w") as fh:
        fh.write("contig\tgenome\tlength\tcategory\tlength_class\tN_cell\tn_cell_T\n")
        for c in sorted(T):
            g, ln, cat, lc = rows[c]
            fh.write(f"{c}\t{g}\t{ln}\t{cat}\t{lc}\t{avail[T[c]]}\t{at[T[c]]}\n")
    w = None
    with open(os.path.join(C, "pool", "pool.fa")) as fi, open(os.path.join(SD, "T.fa"), "w") as fo:
        for line in fi:
            if line[0] == ">":
                w = line[1:].split()[0] in T
            if w:
                fo.write(line)
    print(f"T: {len(T)} contigs in {len(rem)} cells (S3 reproduced: {len(s3)})")
    for k in sorted(rem):
        print(f"  {k[0]:13s} {k[1]:7s} avail-after-S3 {rem[k]:7d} T {at[k]}")


def run(name, extra, stage="T"):
    d = os.path.join(C, "runs", f"{stage}_{name}")
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(SD, f"{stage}.fa")
    U.warm([U.AIFDB])
    aif = os.environ.get("AIF_BIN", U.AIF)
    cmd = [aif, "db", "search", U.AIFDB, fa, "-n", "5", "-t", U.THREADS, "--segments", os.path.join(d, "aif.seg.tsv")]
    t, m = U.timed(cmd + extra, out=os.path.join(d, "aif.tsv"), err=os.path.join(d, "aif.log"), cpus=U.CPUS or None)
    json.dump({"binary": aif, "args": extra, "s": t, "kb": m}, open(os.path.join(d, "run.json"), "w"))
    print(f"{stage}_{name}: {t:.1f} s, {m / 2**20:.1f} GB, args {extra}")


def truth(name, stage="T"):
    os.environ["HITS_FROM"] = f"{stage}_{name}"
    os.environ["TOOLS"] = "aif"
    U.truth(stage)


def members(stage):
    if stage == "T":
        return {r["contig"]: r for r in csv.DictReader(open(os.path.join(SD, "T.tsv")), delimiter="\t")}
    return {r["contig"]: r for r in csv.DictReader(open(os.path.join(SD, "sample.tsv")), delimiter="\t") if r[stage] == "1"}


def score(names, stage="T", quiet=False):
    """top-1 and top-5 AF / ANI error of each variant against MUMmer (truth ANI >= 0.90, self reference excluded)"""
    S = members(stage)
    X = U.contexts()
    out = {}
    for conv in ("q", "one", "q90"):
        TR = U.load_truth(stage, conv)
        for name in names:
            H = U.load_aif(os.path.join(C, "runs", f"{stage}_{name}", "aif.tsv"))
            for lab, k in (("top1", 1), ("top5", 5)):
                rows = []
                for c in S:
                    tr = {r: (a, f) for r, a, f in TR.get(c, [])}
                    sref = X[S[c]["genome"]][1]
                    for h in [h for h in H.get(c, []) if h[0] != sref][:k]:
                        if h[0] in tr and tr[h[0]][0] >= 0.90:
                            rows.append((S[c]["length_class"], h[1] - tr[h[0]][0], h[2] - tr[h[0]][1], tr[h[0]][1],
                                         h[5] <= tr[h[0]][1] <= h[6], h[3] <= tr[h[0]][0] <= h[4]))
                n = max(1, len(rows))
                e = {"pairs": len(rows), "af_mae": 100 * sum(abs(r[2]) for r in rows) / n,
                     "af_bias": 100 * sum(r[2] for r in rows) / n, "ani_mae": 100 * sum(abs(r[1]) for r in rows) / n,
                     "ani_bias": 100 * sum(r[1] for r in rows) / n, "af_cover": sum(r[4] for r in rows) / n,
                     "ani_cover": sum(r[5] for r in rows) / n}
                for _, _, ln in U.LENS:
                    sel = [r for r in rows if r[0] == ln]
                    if sel:
                        e[ln] = (len(sel), round(100 * sum(abs(r[2]) for r in sel) / len(sel), 2),
                                 round(100 * sum(r[2] for r in sel) / len(sel), 2))
                # bias by truth AF
                for lo, hi in ((0, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 1.01)):
                    sel = [r for r in rows if lo <= r[3] < hi]
                    if sel:
                        e[f"tAF{lo}-{min(hi, 1)}"] = (len(sel), round(100 * sum(r[2] for r in sel) / len(sel), 2))
                out[(conv, name, lab)] = e
    if not quiet:
        for (conv, name, lab), e in out.items():
            if lab != "top1":
                continue
            print(f"{conv:4s} {name:14s} n {e['pairs']:5d} AF mae {e['af_mae']:5.2f} bias {e['af_bias']:+5.2f} cover "
                  f"{e['af_cover']:.3f} | ANI mae {e['ani_mae']:.3f} bias {e['ani_bias']:+.3f} cover {e['ani_cover']:.3f} | "
                  + " ".join(f"{k}:{v}" for k, v in e.items() if isinstance(v, tuple)))
    return out


def diag_job(job):
    ref, rpath, qfa = job
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import truth as T
    T.MUMMER = os.path.join(U.ROOT, "data", "envs", "mummer", "bin")
    with tempfile.TemporaryDirectory(dir=os.path.join(C, "tmp")) as td:
        rfa = os.path.join(td, "ref.fa")
        with gzip.open(rpath, "rb") as fi, open(rfa, "wb") as fo:
            fo.write(fi.read())
        pre = os.path.join(td, "x")
        subprocess.run([os.path.join(T.MUMMER, "nucmer"), "-t", "1", "--maxmatch", "--prefix", pre, rfa, qfa], check=True,
                       capture_output=True)
        rows = T.coords(pre + ".delta", "-q")
    os.unlink(qfa)
    return ref, [(n, s - 1, e, idy) for n, s, e, lq, idy, _ in rows]


def ivunion(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def ivlen(iv):
    return sum(b - a for a, b in iv)


def ivinter(A, B):
    i = j = 0
    out = []
    while i < len(A) and j < len(B):
        a, b = max(A[i][0], B[j][0]), min(A[i][1], B[j][1])
        if a < b:
            out.append([a, b])
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def ivminus(A, B):
    out = []
    for a, b in A:
        cur = a
        for c, d in B:
            if d <= cur or c >= b:
                continue
            if c > cur:
                out.append([cur, c])
            cur = max(cur, d)
        if cur < b:
            out.append([cur, b])
    return out


def diag(name, stage="T"):
    d = os.path.join(C, "runs", f"{stage}_{name}")
    S = members(stage)
    X = U.contexts()
    H = U.load_aif(os.path.join(d, "aif.tsv"))
    seg = {}
    for line in open(os.path.join(d, "aif.seg.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] == "contig":
            continue
        seg.setdefault((f[0], f[1]), []).append((int(f[2]), int(f[3]), int(f[4]), float(f[5]), int(f[8]), int(f[9])))
    # the segments file holds the top hit over all references (self included); keep pairs whose reference is not
    # the genome's self reference
    pairs = {(c, r) for (c, r) in seg if c in S and r != X[S[c]["genome"]][1]}
    cache = os.path.join(d, "diag_coords.json")
    if os.path.exists(cache):
        CO = {tuple(k.split("\t")): v for k, v in json.load(open(cache)).items()}
    else:
        seqs, nm = {}, None
        for line in open(os.path.join(SD, f"{stage}.fa")):
            if line[0] == ">":
                nm = line[1:].split()[0]
                seqs[nm] = []
            else:
                seqs[nm].append(line.strip())
        byref = {}
        for c, r in pairs:
            byref.setdefault(r, []).append(c)
        paths = {acc: p for p, acc in (l.rstrip("\n").split("\t") for l in open(U.REPS))}
        os.makedirs(os.path.join(C, "tmp"), exist_ok=True)
        jobs = []
        for r, cs in byref.items():
            q = os.path.join(C, "tmp", f"dq_{r}.fa")
            with open(q, "w") as fh:
                for c in cs:
                    fh.write(f">{c}\n{''.join(seqs[c])}\n")
            jobs.append((r, paths[r], q))
        CO = {p: [] for p in pairs}
        with Pool(int(os.environ.get("JOBS", "56"))) as pool:
            for r, rows in pool.imap_unordered(diag_job, jobs, chunksize=1):
                for n, a, b, idy in rows:
                    if (n, r) in CO:
                        CO[(n, r)].append((a, b, idy))
        json.dump({"\t".join(k): v for k, v in CO.items()}, open(cache, "w"))
    tot = {"pairs": 0, "bases": 0, "aif": 0, "mum": 0, "both": 0, "edge": 0, "gap": 0, "island": 0, "mum_only": 0}
    idbin = {}  # aif window identity of aif-only vs shared aligned bases
    perlen = {}
    for (c, r), co in CO.items():
        if not co:
            continue  # MUMmer found nothing: the pair is outside the estimation set (truth ANI unknown)
        L = int(S[c]["length"])
        M = ivunion([(a, b) for a, b, _ in co])
        ani_t = sum((b - a) * i for a, b, i in co) / sum(b - a for a, b, _ in co) / 100
        if ani_t < 0.90:
            continue
        A = ivunion([(s, e) for s, e, al, *_ in seg[(c, r)] if al])
        both = ivinter(A, M)
        aonly = ivminus(A, M)
        tot["pairs"] += 1
        tot["bases"] += L
        tot["aif"] += ivlen(A)
        tot["mum"] += ivlen(M)
        tot["both"] += ivlen(both)
        tot["mum_only"] += ivlen(ivminus(M, A))
        lc = S[c]["length_class"]
        pl = perlen.setdefault(lc, {"pairs": 0, "bases": 0, "aonly": 0, "monly": 0})
        pl["pairs"] += 1
        pl["bases"] += L
        pl["aonly"] += ivlen(aonly)
        pl["monly"] += ivlen(ivminus(M, A))
        for a, b in aonly:
            # island: the aif segment containing it has no MUMmer overlap; gap: MUMmer on both sides within the
            # same aif segment; edge: otherwise (extension past one MUMmer end, or at a contig end)
            sg = next(x for x in A if x[0] <= a and b <= x[1])
            inside = ivinter([sg], M)
            if not inside:
                tot["island"] += b - a
            elif any(x[1] <= a for x in inside) and any(x[0] >= b for x in inside):
                tot["gap"] += b - a
            else:
                tot["edge"] += b - a
        for s, e, al, ani, np_, nh in seg[(c, r)]:
            if not al or ani != ani:
                continue
            ov = ivlen(ivinter([[s, e]], M))
            key = "<0.86" if ani < 0.86 else "0.86-0.90" if ani < 0.90 else "0.90-0.95" if ani < 0.95 else ">=0.95"
            v = idbin.setdefault(key, [0, 0])
            v[0] += ov
            v[1] += (e - s) - ov
    B = max(1, tot["bases"])
    print(f"{stage}_{name}: {tot['pairs']} top-hit pairs (truth ANI >= 0.90), {tot['bases']} contig bases")
    print(f"  aif aligned {tot['aif'] / B:.4f}  MUMmer aligned {tot['mum'] / B:.4f}  both {tot['both'] / B:.4f}")
    print(f"  aif-only {(tot['aif'] - tot['both']) / B:.4f} = edge {tot['edge'] / B:.4f} + bridged gap {tot['gap'] / B:.4f}"
          f" + island {tot['island'] / B:.4f};  MUMmer-only {tot['mum_only'] / B:.4f}")
    for lc, v in sorted(perlen.items()):
        print(f"  {lc:7s} pairs {v['pairs']:5d} aif-only {v['aonly'] / v['bases']:.4f} MUMmer-only {v['monly'] / v['bases']:.4f}")
    print("  aif window identity: bases shared with MUMmer / aif-only")
    for key in ("<0.86", "0.86-0.90", "0.90-0.95", ">=0.95"):
        if key in idbin:
            sh, ao = idbin[key]
            print(f"    {key:9s} shared {sh / B:.4f} aif-only {ao / B:.4f} (aif-only share {ao / max(1, sh + ao):.2f})")
    return tot


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "tune":
        tune(int(sys.argv[2]) if len(sys.argv) > 2 else NT)
    elif cmd == "run":
        run(sys.argv[2], sys.argv[3:], os.environ.get("STAGE", "T"))
    elif cmd == "truth":
        truth(sys.argv[2], os.environ.get("STAGE", "T"))
    elif cmd == "score":
        score(sys.argv[2:], os.environ.get("STAGE", "T"))
    elif cmd == "diag":
        diag(sys.argv[2], os.environ.get("STAGE", "T"))
    else:
        raise SystemExit(__doc__)
