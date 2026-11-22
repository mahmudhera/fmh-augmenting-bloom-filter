#!/usr/bin/env python3
"""Calibration of aif's per-contig confidence intervals (ANI and AF).

Statistical intervals (aif --ci-stat): half-width = z * sigma, where sigma is the statistical SE
(Theorem 1 conditioned on the realised substitutions for ANI; posterior boundary variance for
AF) with a small floor. Calibrated intervals (aif default): stratified split-conformal. The
nonconformity score is |estimate - truth| / sigma; for each stratum (estimated identity
>= 0.99 / 0.97 / 0.93 / lower x length < 1 / 3 kb / longer x AF estimate >= 0.999 or not) and
level (90, 95, 99%), the multiplier q is the conformal quantile of the scores on the REAL TUNING
SET (MUMmer truth), ceil((n + 1) level)-th smallest. Strata with fewer than MIN_N contigs fall
back to the pooled (identity x AF class), then (identity), then all contigs.

  python3 bench_real/calibration.py fit        # fit on data/real_tune, write src/ci_table.inc
  python3 bench_real/calibration.py evaluate   # coverage of statistical and calibrated intervals
                                               # on all sets -> data/real/results/calibration.json
"""
import csv
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
REAL_SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]
SIM_G = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]
LEVELS = [0.90, 0.95, 0.99]
MIN_N = 50
OUT = os.path.join(ROOT, "data", "real", "results", "calibration")


def stratum(ani, af, L):
    ia = 0 if ani >= 0.99 else 1 if ani >= 0.97 else 2 if ani >= 0.93 else 3
    il = 0 if L < 1000 else 1 if L < 3000 else 2
    ic = 0 if af >= 0.999 else 1
    return (ia * 3 + il) * 2 + ic


def datasets():
    out = []
    for tag, rdir in (("real_tune", "data/real_tune"), ("real_eval", "data/real")):
        for s in REAL_SETS:
            refs = sorted(set(l.split("\t")[4] for l in list(open(f"{rdir}/sets/{s}/genomes.tsv"))[1:]))
            for ref in refs:
                out.append((tag, s, f"{rdir}/results/aif/{s}.{ref}.aif", f"{rdir}/sets/{s}/{ref}.contigs.fa"))
    for tag, d in (("sim_tune", "data/tune"), ("sim_main", "data/bench"), ("sim_fresh", "data/bench2")):
        for g in SIM_G:
            out.append((tag, g, f"data/results/main/aif/{g}.aif", f"{d}/{g}.contigs.fa"))
    return out


def truth_all():
    T = {}
    for tag, rdir in (("real_tune", "data/real_tune"), ("real_eval", "data/real")):
        for s in REAL_SETS:
            for r in csv.DictReader(open(f"{rdir}/sets/{s}/truth.tsv"), delimiter="\t"):
                T[(tag, r["contig"])] = (float(r["true_ani"]), float(r["true_af"]), int(r["length"]), s)
    for tag, d in (("sim_tune", "data/tune"), ("sim_main", "data/bench"), ("sim_fresh", "data/bench2")):
        for r in csv.DictReader(open(f"{d}/truth.tsv"), delimiter="\t"):
            T[(tag, r["contig"])] = (float(r["true_ani"]), float(r["true_af"]), int(r["length"]), r["condition"])
    return T


def run_all(mode, extra):
    os.makedirs(OUT, exist_ok=True)

    def one(job):
        tag, s, sk, q = job
        out = f"{OUT}/{mode}.{tag}.{s}.{os.path.basename(sk)}.tsv"
        with open(out, "w") as fh:
            subprocess.run([AIF, "query", sk, q, "-t", "4"] + extra, stdout=fh, check=True)
        rows = []
        for r in csv.DictReader(open(out), delimiter="\t"):
            if r["status"] == "ok":
                rows.append((tag, r["contig"].split()[0], float(r["ani"]), float(r["af"]), float(r["ani_lo"]),
                             float(r["ani_hi"]), float(r["af_lo"]), float(r["af_hi"]), float(r["ani_se"]), float(r["af_se"])))
        return rows

    with ThreadPoolExecutor(8) as ex:
        return [r for rows in ex.map(one, datasets()) for r in rows]


def conformal_q(scores, level):
    s = sorted(scores)
    n = len(s)
    k = min(n, int(math.ceil((n + 1) * level)))
    return s[k - 1]


def fit():
    T = truth_all()
    rows = [r for r in run_all("stat", ["--ci-stat", "--ci-level", "0.95"]) if r[0] == "real_tune"]
    z95 = 1.959964
    recs = []
    for tag, c, ani, af, alo, ahi, flo, fhi, ase, fse in rows:
        ta, tf, L, _ = T[(tag, c)]
        st = stratum(ani, af, L)
        sa = ase  # --ci-stat at 95%: reported SE equals sigma (statistical SE with floor)
        sf = fse
        s_ani = abs(ani - ta) / sa if (tf >= 0.2 and ta == ta) else None
        s_af = abs(af - tf) / sf
        recs.append((st, s_ani, s_af))
    table = {}
    for li, lv in enumerate(LEVELS):
        for st in range(24):
            ia, rest = divmod(st, 6)
            il, ic = divmod(rest, 2)
            pools = [lambda r: r[0] == st,
                     lambda r: r[0] // 6 == ia and r[0] % 2 == ic,
                     lambda r: r[0] // 6 == ia,
                     lambda r: True]
            q = []
            for j in (1, 2):
                for pool in pools:
                    sc = [r[j] for r in recs if pool(r) and r[j] is not None]
                    if len(sc) >= MIN_N:
                        q.append(conformal_q(sc, lv))
                        break
            table[(li, st)] = q
    with open(os.path.join(ROOT, "src", "ci_table.inc"), "w") as fh:
        fh.write("// generated by bench_real/calibration.py fit: conformal multipliers fitted on the real tuning set\n")
        for li in range(len(LEVELS)):
            fh.write("  {" + ", ".join("{%.4f, %.4f}" % tuple(table[(li, st)]) for st in range(24)) + "},\n")
    with open(os.path.join(ROOT, "data", "real", "results", "ci_table.tsv"), "w") as fh:
        fh.write("level\tstratum\tq_ani\tq_af\n")
        for (li, st), q in sorted(table.items()):
            fh.write(f"{li}\t{st}\t{q[0]:.4f}\t{q[1]:.4f}\n")
    print("wrote src/ci_table.inc; rebuild aif (make aif)")


def evaluate():
    T = truth_all()
    res = {}
    for mode, extra in (("stat", ["--ci-stat"]), ("calibrated", [])):
        res[mode] = {}
        for lv in LEVELS:
            rows = run_all(f"{mode}{lv}", extra + ["--ci-level", str(lv)])
            per = {}
            for tag, c, ani, af, alo, ahi, flo, fhi, ase, fse in rows:
                ta, tf, L, grp = T[(tag, c)]
                d = per.setdefault(tag, {"ani": [0, 0], "af": [0, 0], "grp": {}, "strata": {}, "wa": [], "wf": []})
                g = d["grp"].setdefault(grp, {"ani": [0, 0], "af": [0, 0]})
                st = d["strata"].setdefault(stratum(ani, af, L), {"ani": [0, 0], "af": [0, 0], "wa": [], "wf": []})
                if tf >= 0.2 and ta == ta:
                    hit = alo - 1e-9 <= ta <= ahi + 1e-9
                    for x in (d, g, st):
                        x["ani"][0] += hit
                        x["ani"][1] += 1
                    d["wa"].append((ahi - alo) / 2 * 100)
                    st["wa"].append((ahi - alo) / 2 * 100)
                hit = flo - 1e-9 <= tf <= fhi + 1e-9
                for x in (d, g, st):
                    x["af"][0] += hit
                    x["af"][1] += 1
                d["wf"].append((fhi - flo) / 2 * 100)
                st["wf"].append((fhi - flo) / 2 * 100)
            out = {}
            for tag, d in per.items():
                med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")  # noqa: E731
                out[tag] = {
                    "ani": d["ani"][0] / d["ani"][1], "af": d["af"][0] / d["af"][1], "n_ani": d["ani"][1], "n_af": d["af"][1],
                    "halfwidth_ani_med": med(d["wa"]), "halfwidth_af_med": med(d["wf"]),
                    "groups": {k: {"ani": v["ani"][0] / max(1, v["ani"][1]), "af": v["af"][0] / max(1, v["af"][1]),
                                   "n_ani": v["ani"][1], "n_af": v["af"][1]} for k, v in d["grp"].items()},
                    "strata": {str(k): {"ani": v["ani"][0] / max(1, v["ani"][1]), "af": v["af"][0] / max(1, v["af"][1]),
                                        "n_ani": v["ani"][1], "n_af": v["af"][1], "hw_ani": med(v["wa"]), "hw_af": med(v["wf"])}
                               for k, v in d["strata"].items()}}
            res[mode][str(lv)] = out
            for tag in ["real_tune", "real_eval", "sim_tune", "sim_main", "sim_fresh"]:
                o = out[tag]
                print(f"{mode:10s} {lv:.2f} {tag:9s} ANI cov {o['ani']:.3f} (n={o['n_ani']}, median half-width {o['halfwidth_ani_med']:.3f} pp) "
                      f"| AF cov {o['af']:.3f} (n={o['n_af']}, median half-width {o['halfwidth_af_med']:.2f} pp)", flush=True)
    json.dump(res, open(os.path.join(ROOT, "data", "real", "results", "calibration.json"), "w"), indent=1)


if __name__ == "__main__":
    {"fit": fit, "evaluate": evaluate}[sys.argv[1] if len(sys.argv) > 1 else "evaluate"]()
