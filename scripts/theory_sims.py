#!/usr/bin/env python3
"""Checks of the theory against the aif binary (aif_theory.pdf).

Contigs are simulated from a random 5 Mbp reference (data/theory/ref.fa, iid uniform bases, sketched with -k 20) so
that the idealised model holds except for what the experiment varies. Each experiment writes contigs and truth to
data/theory/<exp>/, runs bin/aif query, and stores a summary in results_theory/sims.json.

Usage: /usr/bin/python3 -I scripts/theory_sims.py [exp ...]   (exp: null detect boundary absorb het indel sigma)
"""
import json
import math
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "theory")
RES = os.path.join(ROOT, "results_theory", "sims.json")
AIF = os.path.join(ROOT, "bin", "aif")
SK = os.path.join(D, "ref_k20.aif")
THREADS = "16"
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
COLS = "contig ani af theta probes hits aligned_probes status lr het low_frac ani_lo ani_hi af_lo af_hi ani_se af_se".split()


def load_ref():
    seq = "".join(l.strip() for l in open(os.path.join(D, "ref.fa")) if not l.startswith(">"))
    return np.frombuffer(seq.encode(), dtype=np.uint8).copy()


def code(a):
    return np.searchsorted(ACGT, a)


def mutate(rng, seg, p):
    """substitutions at rate p (each to one of the other three bases); returns (seq, substitutions)"""
    seg = seg.copy()
    m = rng.random(len(seg)) < p
    c = code(seg[m])
    seg[m] = ACGT[(c + rng.integers(1, 4, m.sum())) % 4]
    return seg, int(m.sum())


def randseq(rng, n):
    return ACGT[rng.integers(0, 4, n)]


def write_fa(path, recs):
    with open(path, "w") as f:
        for name, s in recs:
            f.write(f">{name}\n{s.tobytes().decode()}\n")


def run(fa, *opts):
    out = subprocess.run([AIF, "query", SK, fa, "-t", THREADS, *opts], capture_output=True, text=True, check=True).stdout
    rows = {}
    for line in out.splitlines():
        p = line.split("\t")
        if p[0] == "contig" or len(p) < len(COLS):
            continue
        r = dict(zip(COLS, p))
        for c in COLS[1:]:
            if c != "status":
                try:
                    r[c] = float(r[c])
                except ValueError:
                    r[c] = float("nan")
        rows[r["contig"]] = r
    return rows


def segment(rng, ref, L):
    s = int(rng.integers(0, len(ref) - L))
    return ref[s:s + L]


def exp_null(rng, ref):
    """unrelated contigs: how often does the likelihood ratio reach lambda? (bound: exp(-lambda))"""
    out = {}
    d = os.path.join(D, "null")
    os.makedirs(d, exist_ok=True)
    for L, n in [(1000, 40000), (5000, 20000), (20000, 5000)]:
        fa = os.path.join(d, f"null_{L}.fa")
        write_fa(fa, [(f"n{i}", randseq(rng, L)) for i in range(n)])
        rows = run(fa)
        lr = np.array([r["lr"] for r in rows.values()])
        rep = sum(r["status"] == "ok" for r in rows.values())
        out[L] = {"n": n, "reported": int(rep), "frac_ge": {str(l): float((lr >= l).mean()) for l in range(0, 11)},
                  "max_lr": float(np.nanmax(lr)), "status": {s: int(sum(r["status"] == s for r in rows.values()))
                                                            for s in {r["status"] for r in rows.values()}}}
        print("null", L, out[L])
    return out


def exp_detect(rng, ref):
    out = []
    d = os.path.join(D, "detect")
    os.makedirs(d, exist_ok=True)
    recs, truth = [], {}
    for L in [500, 1000, 2000, 5000, 10000, 20000]:
        for a in [0.84, 0.86, 0.88, 0.90, 0.92, 0.95]:
            for i in range(200):
                s, _ = mutate(rng, segment(rng, ref, L), 1 - a)
                nm = f"L{L}_a{a}_{i}"
                recs.append((nm, s))
                truth[nm] = (L, a)
    fa = os.path.join(d, "detect.fa")
    write_fa(fa, recs)
    rows = run(fa)
    for L in [500, 1000, 2000, 5000, 10000, 20000]:
        for a in [0.84, 0.86, 0.88, 0.90, 0.92, 0.95]:
            rr = [rows[n] for n, t in truth.items() if t == (L, a)]
            lr = np.array([r["lr"] for r in rr])
            out.append({"L": L, "ani": a, "reported": float(np.mean([r["status"] in ("ok", "below_floor") for r in rr])),
                        "lr_median": float(np.median(lr))})
            print("detect", out[-1])
    return out


def exp_boundary(rng, ref):
    """5 kb aligned at identity a, then 5 kb unrelated: AF error in bases"""
    out = []
    d = os.path.join(D, "boundary")
    os.makedirs(d, exist_ok=True)
    recs, truth = [], {}
    for a in [0.85, 0.90, 0.95, 0.99]:
        for i in range(400):
            s, _ = mutate(rng, segment(rng, ref, 5000), 1 - a)
            nm = f"a{a}_{i}"
            recs.append((nm, np.concatenate([s, randseq(rng, 5000)])))
            truth[nm] = a
    fa = os.path.join(d, "boundary.fa")
    write_fa(fa, recs)
    for tag, opts in [("v17", ()), ("v16", ("--af-local", "0"))]:
        rows = run(fa, *opts)
        for a in [0.85, 0.90, 0.95, 0.99]:
            err = np.array([(rows[n]["af"] - 0.5) * 10000 for n, t in truth.items() if t == a and rows[n]["status"] != "no_homology"])
            out.append({"mode": tag, "ani": a, "n": int(len(err)), "mean_bp": float(err.mean()), "sd_bp": float(err.std()),
                        "median_bp": float(np.median(err))})
            print("boundary", out[-1])
    return out


def exp_absorb(rng, ref, a=0.97, A=3000):
    """aligned A + unrelated gap G + aligned A (internal), and aligned A + unrelated flank G (edge)"""
    out = []
    d = os.path.join(D, "absorb")
    os.makedirs(d, exist_ok=True)
    Gs = [100, 200, 500, 1000, 1500, 2000, 3000, 5000]
    recs, truth = [], {}
    for G in Gs:
        for i in range(200):
            seg = segment(rng, ref, 2 * A)
            s1, _ = mutate(rng, seg[:A], 1 - a)
            s2, _ = mutate(rng, seg[A:], 1 - a)
            nm = f"int_G{G}_{i}"
            recs.append((nm, np.concatenate([s1, randseq(rng, G), s2])))
            truth[nm] = ("int", G, 2 * A)
            s3, _ = mutate(rng, segment(rng, ref, A), 1 - a)
            nm = f"edge_G{G}_{i}"
            recs.append((nm, np.concatenate([s3, randseq(rng, G)])))
            truth[nm] = ("edge", G, A)
    fa = os.path.join(d, "absorb.fa")
    write_fa(fa, recs)
    for tag, opts in [("v17", ()), ("v16", ("--af-local", "0"))]:
        rows = run(fa, *opts)
        for kind in ("int", "edge"):
            for G in Gs:
                v = [(rows[n]["af"] * (t[2] + G) - t[2]) / G for n, t in truth.items() if t[0] == kind and t[1] == G]
                out.append({"mode": tag, "kind": kind, "G": G, "absorbed_share": float(np.mean(v)),
                            "median": float(np.median(v))})
                print("absorb", out[-1])
    return out


def exp_het(rng, ref, L=10000):
    out = []
    d = os.path.join(D, "het")
    os.makedirs(d, exist_ok=True)
    cases = [(0.99, 0.95), (0.99, 0.90), (0.97, 0.93), (0.995, 0.85), (0.95, 0.95)]
    recs, truth = [], {}
    for (a1, a2) in cases:
        for i in range(200):
            seg = segment(rng, ref, L)
            s1, m1 = mutate(rng, seg[:L // 2], 1 - a1)
            s2, m2 = mutate(rng, seg[L // 2:], 1 - a2)
            nm = f"h{a1}_{a2}_{i}"
            recs.append((nm, np.concatenate([s1, s2])))
            truth[nm] = ((a1, a2), 1 - (m1 + m2) / L)
    fa = os.path.join(d, "het.fa")
    write_fa(fa, recs)
    r_auto, r_pool = run(fa), run(fa, "--het", "0")
    K = 20
    for c in cases:
        nms = [n for n, t in truth.items() if t[0] == c]
        tru = np.array([truth[n][1] for n in nms])
        auto = np.array([r_auto[n]["ani"] for n in nms])
        pool = np.array([r_pool[n]["ani"] for n in nms])
        het = np.mean([r_auto[n]["het"] for n in nms])
        pm = (0.5 * c[0] ** K + 0.5 * c[1] ** K) ** (1 / K)
        out.append({"a1": c[0], "a2": c[1], "true": float(tru.mean()), "pooled": float(np.nanmean(pool)),
                    "auto": float(np.nanmean(auto)), "het_called": float(het), "power_mean": pm,
                    "pooled_err_pp": float(100 * np.nanmean(pool - tru)), "auto_err_pp": float(100 * np.nanmean(auto - tru)),
                    "auto_mae_pp": float(100 * np.nanmean(np.abs(auto - tru)))})
        print("het", out[-1])
    return out


def exp_indel(rng, ref, L=10000, p=0.03):
    out = []
    d = os.path.join(D, "indel")
    os.makedirs(d, exist_ok=True)
    recs, truth = [], {}
    for r in [0.0, 0.002, 0.005, 0.01]:
        for i in range(200):
            seg = segment(rng, ref, L + 200)
            s, m = mutate(rng, seg, p)
            # single-base indels at rate r per base (half insertions, half deletions)
            ev = rng.random(len(s)) < r
            ins = ev & (rng.random(len(s)) < 0.5)
            dele = ev & ~ins
            parts, j = [], 0
            idx = np.flatnonzero(ev)
            for q in idx:
                parts.append(s[j:q])
                if ins[q]:
                    parts.append(s[q:q + 1])
                    parts.append(randseq(rng, 1))
                j = q + 1
                if ins[q]:
                    continue
            parts.append(s[j:])
            out_s = np.concatenate(parts)[:L]
            frac = len(out_s) / len(s)
            n_ins, n_del = int(ins.sum() * frac), int(dele.sum() * frac)
            subs = m * frac
            cols = L + n_del  # alignment columns: query bases plus deleted reference bases
            ex = 1 - subs / (cols - n_ins - n_del)
            inc = 1 - (subs + n_ins + n_del) / cols
            nm = f"r{r}_{i}"
            recs.append((nm, out_s))
            truth[nm] = (r, ex, inc)
    fa = os.path.join(d, "indel.fa")
    write_fa(fa, recs)
    rows = run(fa)
    for r in [0.0, 0.002, 0.005, 0.01]:
        nms = [n for n, t in truth.items() if t[0] == r]
        an = np.array([rows[n]["ani"] for n in nms])
        ex = np.array([truth[n][1] for n in nms])
        inc = np.array([truth[n][2] for n in nms])
        out.append({"r": r, "aif": float(np.nanmean(an)), "excl": float(ex.mean()), "incl": float(inc.mean()),
                    "pred": float((1 - p) * (1 - r)), "af": float(np.nanmean([rows[n]["af"] for n in nms]))})
        print("indel", out[-1])
    return out


def exp_sigma():
    """calibrated spurious-match rate vs the random-sequence value 1 - exp(-2N/4^k)"""
    out = []
    gen = [("random 5 Mbp", os.path.join(D, "ref.fa"))] + [
        (g.replace(".fa", ""), os.path.join(ROOT, "data", "genomes", g)) for g in sorted(os.listdir(os.path.join(ROOT, "data", "genomes")))]
    for name, fa in gen:
        for k in range(17, 22):
            sk = os.path.join(D, "sig", f"{os.path.basename(fa)}.k{k}.aif")
            os.makedirs(os.path.dirname(sk), exist_ok=True)
            if not os.path.exists(sk):
                subprocess.run([AIF, "sketch", fa, "-o", sk, "-k", str(k)], capture_output=True, check=True)
            info = subprocess.run([AIF, "info", sk], capture_output=True, text=True, check=True).stdout
            kv = dict(x.split("=", 1) for x in info.split() if "=" in x)
            out.append({"genome": name, "k": k, "bases": int(kv["bases"]), "analytic": float(kv["sigma_analytic"]),
                        "emp": float(kv["sigma_emp"])})
            print("sigma", out[-1])
    return out


def prepare():
    """the random 5 Mbp reference (seed 20261127) and its sketches at k = 17..21"""
    os.makedirs(D, exist_ok=True)
    fa = os.path.join(D, "ref.fa")
    if not os.path.exists(fa):
        r = np.random.default_rng(20261127)
        s = ACGT[r.integers(0, 4, 5_000_000)].tobytes().decode()
        with open(fa, "w") as f:
            f.write(">ref\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + "\n")
    for k in range(17, 22):
        sk = os.path.join(D, f"ref_k{k}.aif")
        if not os.path.exists(sk):
            subprocess.run([AIF, "sketch", fa, "-o", sk, "-k", str(k)], capture_output=True, check=True)


def main():
    if sys.argv[1:] == ["prepare"]:
        prepare()
        return
    exps = sys.argv[1:] or ["null", "detect", "boundary", "absorb", "het", "indel", "sigma"]
    res = json.load(open(RES)) if os.path.exists(RES) else {}
    ref = load_ref()
    for e in exps:
        rng = np.random.default_rng({"null": 1, "detect": 2, "boundary": 3, "absorb": 4, "het": 5, "indel": 6, "sigma": 7}[e])
        res[e] = exp_sigma() if e == "sigma" else globals()[f"exp_{e}"](rng, ref)
        json.dump(res, open(RES, "w"), indent=1)
    print("wrote", RES)


if __name__ == "__main__":
    main()
