#!/usr/bin/env python3
"""Build tables (tbl) and figures (EPS) for the aif-vs-skani benchmark report.

Reads data/results/{main,fresh}/..., data/bench*/truth.tsv and writes
report_bench/gen/*.tbl, report_bench/gen/numbers.txt, report_bench/figs/*.eps.
Usage: python3 scripts/benchmark_report.py   (run from the repository root)
"""
import csv
import json
import math
import os
import subprocess
import tempfile
from collections import OrderedDict, defaultdict

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "data", "results", "main")
GEN = os.path.join(ROOT, "report_bench", "gen")
FIG = os.path.join(ROOT, "report_bench", "figs")
SKANI = "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani"
GENOMES = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]
GNAME = {"NC_000913.3": "E. coli K-12", "NC_000964.3": "B. subtilis 168",
         "NC_002516.2": "P. aeruginosa PAO1", "NC_007795.1": "S. aureus NCTC 8325"}
C_AIF, C_SK, INK, INK2, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d9d8d4"


def tbl(path, fmt_header, fmt_body, header_rows, rows, opts="center box tab(@)"):
    lines = [".TS", opts + ";"] + fmt_header + [fmt_body + "."]
    for i, hr in enumerate(header_rows):
        fmt = fmt_header[i].split() if i < len(fmt_header) else []
        if fmt and len(fmt) == len(hr):
            hr = [h for h, f in zip(hr, fmt) if f != "s"]
        lines.append("@".join(hr))
    lines.append("_")
    for r in rows:
        lines.append("_" if r == "_" else "@".join(r))
    lines.append(".TE")
    open(path, "w").write("\n".join(lines) + "\n")


def f(x, d=2):
    return "\\(em" if x != x else f"{x:.{d}f}"


def load_eval(path):
    return json.load(open(path))


def load_pred(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["contig"].split()[0]] = (float(r["ani"]), float(r["af"]))
    return out


def theta(p, k):
    return (1 - p) ** k


def varn(n, p, k):
    q = 1 - (1 - p) ** k
    return (n * (1 - q) * (q * (2 * k + 2 / p - 1) - 2 * k) + k * (k - 1) * (1 - q) ** 2
            + 2 * (1 - q) / p ** 2 * ((1 + (k - 1) * (1 - q)) * p - q))


def sketch_info(path):
    out = subprocess.run([os.path.join(ROOT, "bin", "aif"), "info", path], capture_output=True, text=True).stdout
    info = {"tiers": []}
    for line in out.splitlines():
        kv = dict(t.split("=", 1) for t in line.split() if "=" in t)
        if line.startswith("k="):
            info.update(k=int(kv["k"]), bases=int(kv["bases"]), bytes=int(kv["filter_bytes"]),
                        sigma=float(kv["sigma_emp"]), bpb=float(kv["bits/base"]))
        elif line.startswith("tier"):
            info["tiers"].append((float(kv["s"]), float(kv["eps"]), int(kv["nh"])))
    return info


def pred_sd(L, ani, info):
    p = 1 - ani
    k = info["k"]
    n = L - k + 1
    th = theta(p, k)
    inv = 0
    for s, e, _ in info["tiers"]:
        ep = e + info["sigma"] - e * info["sigma"]
        v = (th * (1 - th) * (1 - s) + (1 - th) * ep / (1 - ep)) / (n * s)
        inv += 1 / v
    vc = varn(n, p, k) / n ** 2 + 1 / inv
    return vc * (1 - p) ** 2 / (k ** 2 * th ** 2)


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    nums = OrderedDict()
    truth = list(csv.DictReader(open(os.path.join(ROOT, "data", "bench", "truth.tsv")), delimiter="\t"))
    T = {t["contig"]: t for t in truth}
    conds = OrderedDict()
    for t in truth:
        conds.setdefault(t["condition"], []).append(t)
    cinfo = {r["condition"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "bench", "conditions.tsv")), delimiter="\t")}

    # ---------- conditions table
    rows = []
    for c, r in cinfo.items():
        extra = []
        if float(r["indel_frac"]) > 0:
            extra.append(f"indels {float(r['indel_frac'])*100:.0f}% of subs")
        if r["hetero"] != "None":
            extra.append({"blocks": "1 kb blocks at p=0 / 2p", "gamma": "500 bp blocks, rate ~ Gamma(CV 1)"}[r["hetero"]])
        if r["chimera"] == "True":
            extra.append("two loci $>=$ 1 Mbp apart")
        rows.append([c.replace("_", "\\(ul"), r["scenario"], f"{int(r['length'])/1000:g}", f"{100*float(r['target_ani']):g}",
                     f"{float(r['af']):g}", "; ".join(extra) or "\\(em"])
    tbl(f"{GEN}/tab_conditions.tbl", ["l l c c c l"], "l l n n n l",
        [["condition", "scenario (theory report)", "kb", "ANI %", "AF", "other"]], rows)

    # ---------- evaluations
    E = {}
    for d in os.listdir(R):
        p = os.path.join(R, d, "eval.json")
        if os.path.exists(p):
            E[d] = load_eval(p)
    cmp3 = load_eval(os.path.join(R, "compare3.json"))
    fresh = load_eval(os.path.join(ROOT, "data", "results", "fresh", "eval.json"))

    def ov(name):
        e = E[name]
        return e["overall"][name] if name in e["overall"] else e["overall"]["skani"]

    sk, ai = cmp3["overall"]["skani"], cmp3["overall"]["aif"]
    skr = cmp3["overall"]["skani_relaxed"]
    fsk, fai = fresh["overall"]["skani"], fresh["overall"]["aif"]

    # ---------- per-genome sketch info + skani sketch sizes
    infos = {g: sketch_info(os.path.join(R, "aif", g + ".aif")) for g in GENOMES}
    skbytes = {}
    with tempfile.TemporaryDirectory() as tmp:
        for g in GENOMES:
            d = os.path.join(tmp, g)
            subprocess.run([SKANI, "sketch", "-t", "1", os.path.join(ROOT, "data", "genomes", g + ".fa"), "-o", d],
                           capture_output=True, check=True)
            skbytes[g] = sum(os.path.getsize(os.path.join(d, x)) for x in ("sketches.db2", "markers.bin2", "markers.bin", "index.db2"))
    aif_bytes = {g: os.path.getsize(os.path.join(R, "aif", g + ".aif")) for g in GENOMES}
    tot_bases = sum(infos[g]["bases"] for g in GENOMES)
    nums["skani_bpb"] = sum(skbytes.values()) * 8 / tot_bases
    nums["aif_bpb"] = sum(aif_bytes.values()) * 8 / tot_bases

    # ---------- timing
    tim = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(R, "timing.tsv")), delimiter="\t"):
        tim[r["genome"]][(r["tool"], r["phase"])] = (float(r["median_s"]), int(r["max_rss_kb"]))
    thr = list(csv.DictReader(open(os.path.join(R, "throughput.tsv")), delimiter="\t"))

    # ---------- headline table
    sk_time = sum(tim[g][("skani", "total")][0] for g in GENOMES)
    ai_time = sum(tim[g][("aif", "sketch")][0] + tim[g][("aif", "query")][0] for g in GENOMES)
    sk_mem = max(tim[g][("skani", "total")][1] for g in GENOMES) / 1024
    ai_mem = max(max(tim[g][("aif", "sketch")][1], tim[g][("aif", "query")][1]) for g in GENOMES) / 1024
    nums.update(sk_time=sk_time, ai_time=ai_time, sk_mem=sk_mem, ai_mem=ai_mem)
    rows = [
        ["sketch size (bits / reference base)", f(nums["skani_bpb"], 3), f(nums["aif_bpb"], 3), "\\(em"],
        ["coverage (contigs reported)", f(sk["coverage"], 3), f(ai["coverage"], 3), f(fsk["coverage"], 3) + " / " + f(fai["coverage"], 3)],
        ["ANI MAE, all contigs (pp)", f(sk["ani_mae_micro"], 3), f(ai["ani_mae_micro"], 3), f(fsk["ani_mae_micro"], 3) + " / " + f(fai["ani_mae_micro"], 3)],
        ["ANI MAE, mean over conditions (pp)", f(sk["ani_mae_macro"], 3), f(ai["ani_mae_macro"], 3), f(fsk["ani_mae_macro"], 3) + " / " + f(fai["ani_mae_macro"], 3)],
        ["ANI MAE, paired contigs (pp)", f(sk["paired_ani_mae_macro"], 3), f(ai["paired_ani_mae_macro"], 3), f(fsk["paired_ani_mae_macro"], 3) + " / " + f(fai["paired_ani_mae_macro"], 3)],
        ["AF MAE, all contigs (pp)", f(sk["af_mae_micro"]), f(ai["af_mae_micro"]), f(fsk["af_mae_micro"]) + " / " + f(fai["af_mae_micro"])],
        ["AF MAE, mean over conditions (pp)", f(sk["af_mae_macro"]), f(ai["af_mae_macro"]), f(fsk["af_mae_macro"]) + " / " + f(fai["af_mae_macro"])],
        ["AF MAE, paired contigs (pp)", f(sk["paired_af_mae_macro"]), f(ai["paired_af_mae_macro"]), f(fsk["paired_af_mae_macro"]) + " / " + f(fai["paired_af_mae_macro"])],
        ["wall time, 4 references x 700 contigs (s)", f(sk_time), f(ai_time), "\\(em"],
        ["peak memory (MB)", f(sk_mem, 1), f(ai_mem, 1), "\\(em"],
    ]
    tbl(f"{GEN}/tab_headline.tbl", ["l c c c"], "l n n c",
        [["metric", "skani 0.4.0 (default)", "aif (this work)", "fresh set: skani / aif"]], rows)

    # ---------- per-condition table + win counts
    rows = []
    win_ani = win_af = 0
    nboth = 0
    for c, e in cmp3["conditions"].items():
        a, s = e["tools"]["aif"], e["tools"]["skani"]
        ok = a["ani_mae"] == a["ani_mae"] and s["ani_mae"] == s["ani_mae"]
        if ok:
            nboth += 1
            win_ani += a["ani_mae"] < s["ani_mae"]
            win_af += a["af_mae"] < s["af_mae"]
        rows.append([c.replace("_", "\\(ul"), f(s["coverage"]), f(a["coverage"]), f(s["ani_mae"], 3), "\\f[B]" + f(a["ani_mae"], 3) + "\\f[]",
                     f(s["af_mae"]), "\\f[B]" + f(a["af_mae"]) + "\\f[]", f(s["ani_bias"], 2), f(a["ani_bias"], 2)])
    nums.update(win_ani=win_ani, win_af=win_af, nboth=nboth)
    fw_ani = sum(1 for e in fresh["conditions"].values() if e["tools"]["aif"]["ani_mae"] < e["tools"]["skani"]["ani_mae"])
    fw_af = sum(1 for e in fresh["conditions"].values() if e["tools"]["aif"]["af_mae"] < e["tools"]["skani"]["af_mae"])
    nums.update(fresh_win_ani=fw_ani, fresh_win_af=fw_af, fresh_n=len(fresh["conditions"]))
    tbl(f"{GEN}/tab_percond.tbl", ["l c s c s c s c s", "l c c c c c c c c"], "l n n n n n n n n",
        [["condition", "coverage", "", "ANI MAE (pp)", "", "AF MAE (pp)", "", "ANI bias (pp)", ""],
         ["", "skani", "aif", "skani", "aif", "skani", "aif", "skani", "aif"]], rows)

    # ---------- ablations
    abl = [("aif", "final: adaptive k, tiers 0.005:0.6 + 0.1:0.4, auto heterogeneity"),
           ("abl_k17", "fixed k = 17"), ("abl_k21", "fixed k = 21"),
           ("abl_pooled", "pooled ANI only (no heterogeneity mode)"),
           ("abl_single_tier", "single tier, FPR 0.01 (no hedging)"),
           ("abl_v1_k21_hedge", "k = 21, tiers 0.003:0.55 + 0.58:0.45 (first design)"),
           ("abl_nobias", "no delta-method bias correction")]
    rows = []
    for name, desc in abl:
        o = ov(name)
        rows.append([desc, f(o["coverage"], 3), f(o["ani_mae_micro"], 3), f(o["ani_mae_macro"], 3), f(o["af_mae_micro"]), f(o["af_mae_macro"])])
    tbl(f"{GEN}/tab_ablation.tbl", ["l c c s c s", "l c c c c c"], "l n n n n n",
        [["variant (all at 0.88 bits/base)", "coverage", "ANI MAE (pp)", "", "AF MAE (pp)", ""],
         ["", "", "all", "cond. mean", "all", "cond. mean"]], rows)

    # ---------- budget sweep
    skb = {}
    for line in open(os.path.join(R, "skani_budgets.txt")):
        name, v = line.split()
        skb[name] = float(v.split("=")[1])
    sweep = []
    for name, b in [("budget_0.064", 0.064), ("budget_0.22", 0.22), ("budget_0.44", 0.44), ("aif", 0.88), ("budget_1.76", 1.76)]:
        o = ov(name)
        sweep.append(("aif", b, o))
    for name, key in [("skani", "skani_default"), ("skani_medium", "skani--medium"), ("skani_slow", "skani--slow")]:
        sweep.append((name, skb[key], ov(name)))
    rows = []
    for name, b, o in sweep:
        label = {"aif": "aif", "skani": "skani default (c=125)", "skani_medium": "skani --medium (c=70)",
                 "skani_slow": "skani --slow (c=30)"}[name]
        rows.append([label, f(b, 3), f(o["coverage"], 3), f(o["ani_mae_micro"], 3), f(o["ani_mae_macro"], 3), f(o["af_mae_micro"]), f(o["af_mae_macro"])])
    tbl(f"{GEN}/tab_budget.tbl", ["l c c c s c s", "l c c c c c c"], "l n n n n n n",
        [["tool", "bits/base", "coverage", "ANI MAE (pp)", "", "AF MAE (pp)", ""],
         ["", "", "", "all", "cond. mean", "all", "cond. mean"]], rows)
    nums["skani_slow_bpb"] = skb["skani--slow"]

    # ---------- resources table
    rows = []
    for g in GENOMES:
        i = infos[g]
        st, sm = tim[g][("skani", "total")]
        at, am = tim[g][("aif", "sketch")]
        qt, qm = tim[g][("aif", "query")]
        rows.append([GNAME[g], f"{i['bases']/1e6:.2f}", str(i["k"]), f"{i['sigma']:.4f}", f"{skbytes[g]/1000:.0f}", f"{aif_bytes[g]/1000:.0f}",
                     f(st), f"{sm/1024:.1f}", f(at), f"{am/1024:.1f}", f(qt), f"{qm/1024:.1f}"])
    tbl(f"{GEN}/tab_resources.tbl", ["l c c c c s c s c s c s", "l c c c c c c c c c c c"], "l n n n n n n n n n n n",
        [["reference", "Mbp", "aif", "aif", "sketch (KB)", "", "skani total", "", "aif sketch", "", "aif query", ""],
         ["", "", "$k$", "$sigma$", "skani", "aif", "s", "MB", "s", "MB", "s", "MB"]], rows)
    rows = []
    for r in thr:
        mbps = int(r["bases"]) / 1e6 / float(r["median_s"])
        rows.append([r["tool"], r["phase"], r["threads"], r["contigs"], f"{int(r['bases'])/1e6:.1f}", f(float(r["median_s"]), 3),
                     f"{int(r['max_rss_kb'])/1024:.1f}", f"{mbps:.0f}" if r["phase"] != "sketch" else "\\(em"])
    tbl(f"{GEN}/tab_throughput.tbl", ["l l c c c c c c"], "l l n n n n n n",
        [["tool", "phase", "threads", "contigs", "Mbp", "median s", "peak MB", "Mbp/s"]], rows)
    # speed: v1.0 vs v1.1 vs skani (single thread)
    t10 = list(csv.DictReader(open(os.path.join(ROOT, "data", "results", "timing_v10", "throughput.tsv")), delimiter="\t"))
    tim10 = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(ROOT, "data", "results", "timing_v10", "timing.tsv")), delimiter="\t"):
        tim10[r["genome"]][(r["tool"], r["phase"])] = float(r["median_s"])
    def getthr(rows_, st, tool, phase, threads="1"):
        for r in rows_:
            if r.get("set", "x1") == st and r["tool"] == tool and r["phase"] == phase and r["threads"] == threads:
                return float(r["median_s"])
        return float("nan")
    srows = []
    sk_ref = sum(tim[g][("skani", "total")][0] for g in GENOMES)
    a10_s = sum(tim10[g][("aif", "sketch")] for g in GENOMES); a10_q = sum(tim10[g][("aif", "query")] for g in GENOMES)
    a11_s = sum(tim[g][("aif", "sketch")][0] for g in GENOMES); a11_q = sum(tim[g][("aif", "query")][0] for g in GENOMES)
    srows.append(["4 references x 700 contigs: sketch", "\\(em", f(a10_s, 3), f(a11_s, 3), f(a10_s / a11_s, 1) + "\\(mu"])
    srows.append(["4 references x 700 contigs: query", "\\(em", f(a10_q, 3), f(a11_q, 3), f(a10_q / a11_q, 1) + "\\(mu"])
    srows.append(["4 references x 700 contigs: total", f(sk_ref, 3), f(a10_s + a10_q, 3), f(a11_s + a11_q, 3), f((a10_s + a10_q) / (a11_s + a11_q), 1) + "\\(mu"])
    for st, lab in (("x1", "2,800 contigs (20 Mbp)"), ("x4", "11,200 contigs (81 Mbp)")):
        sk1 = getthr(thr, st, "skani", "sketch+query")
        q10 = getthr(t10, st, "aif", "query"); q11 = getthr(thr, st, "aif", "query")
        s10 = getthr(t10, "x1", "aif", "sketch"); s11 = getthr(thr, "x1", "aif", "sketch")
        srows.append([lab + ": query", "\\(em", f(q10, 3), f(q11, 3), f(q10 / q11, 1) + "\\(mu"])
        srows.append([lab + ": sketch + query", f(sk1, 3), f(s10 + q10, 3), f(s11 + q11, 3), f((s10 + q10) / (s11 + q11), 1) + "\\(mu"])
        sk8 = getthr(thr, st, "skani", "sketch+query", "8"); q8 = getthr(thr, st, "aif", "query", "8")
        srows.append([lab + ": 8 threads (aif query only)", f(sk8, 3), "\\(em", f(q8, 3), "\\(em"])
    tbl(f"{GEN}/tab_speed.tbl", ["l c c c c"], "l n n n n",
        [["workload (seconds; single thread unless noted)", "skani", "aif v1.0", "aif v1.1", "v1.1 speed-up"]], srows)
    nums.update(a10_q=a10_q, a11_q=a11_q, a10_s=a10_s, a11_s=a11_s, sk_ref=sk_ref)

    # ---------- tier parameters
    rows = []
    for g in GENOMES:
        i = infos[g]
        for j, (s, e, nh) in enumerate(i["tiers"]):
            rows.append([GNAME[g] if j == 0 else "", str(i["k"]) if j == 0 else "", f"tier {j+1}", f"{s:.4f}", f"{e:.4f}", str(nh)])
    tbl(f"{GEN}/tab_tiers.tbl", ["l c c c c c"], "l n l n n n",
        [["reference", "$k$", "tier", "sampling $s$", "exact FPR $epsilon$", "hashes $h$"]], rows)

    # ---------- theory vs observed (model-relative RMSE, AF=1 substitution-only conditions)
    paif = load_pred(os.path.join(R, "aif", "predictions.tsv"))
    psk = load_pred(os.path.join(R, "skani", "predictions.tsv"))
    rows = []
    for c in ["len01k_ani95", "len02k_ani95", "len05k_ani95", "len10k_ani95", "len20k_ani95",
              "ani900_5k", "ani930_5k", "ani970_5k", "ani990_5k", "ani995_5k", "ani999_5k"]:
        r = cinfo[c]
        L, ani = int(r["length"]), float(r["target_ani"])
        pv = sum(pred_sd(L, ani, infos[g]) for g in GENOMES) / len(GENOMES)
        ea = [(paif[t["contig"]][0] - ani) ** 2 for t in conds[c] if t["contig"] in paif]
        es = [(psk[t["contig"]][0] - ani) ** 2 for t in conds[c] if t["contig"] in psk]
        rows.append([c.replace("_", "\\(ul"), f(100 * math.sqrt(pv), 3), f(100 * math.sqrt(sum(ea) / len(ea)), 3) if ea else "\\(em",
                     f(100 * math.sqrt(sum(es) / len(es)), 3) if es else "\\(em", f(100 * math.sqrt(ani * (1 - ani) / L), 3)])
    tbl(f"{GEN}/tab_theory.tbl", ["l c c c c"], "l n n n n",
        [["condition", "Theorem 1 sd", "aif RMSE", "skani RMSE", "alignment oracle sd"]], rows)

    # ---------- figures: per-condition dot plots
    names = list(cmp3["conditions"].keys())
    for metric, fname, xlabel in [("ani_mae", "fig_percond_ani.eps", "ANI mean absolute error (percentage points, log scale)"),
                                  ("af_mae", "fig_percond_af.eps", "AF mean absolute error (percentage points, log scale)")]:
        fig, ax = plt.subplots(figsize=(6.4, 6.6))
        ys = list(range(len(names)))[::-1]
        for y, c in zip(ys, names):
            a = cmp3["conditions"][c]["tools"]["aif"][metric]
            s = cmp3["conditions"][c]["tools"]["skani"][metric]
            ca = cmp3["conditions"][c]["tools"]["aif"]["coverage"]
            cs = cmp3["conditions"][c]["tools"]["skani"]["coverage"]
            vals = [v for v in (a, s) if v == v]
            if len(vals) == 2:
                ax.plot([max(a, 0.005), max(s, 0.005)], [y, y], color=GRID, lw=2, zorder=1)
            if s == s:
                ax.plot([max(s, 0.005)], [y], "o", ms=6.5, color=C_SK, mfc=C_SK if cs >= 0.5 else "white", mew=1.4, zorder=2)
            if a == a:
                ax.plot([max(a, 0.005)], [y], "s", ms=6, color=C_AIF, mfc=C_AIF if ca >= 0.5 else "white", mew=1.4, zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels(names, fontsize=6.8, color=INK)
        ax.set_xscale("log")
        ax.set_xlabel(xlabel, fontsize=8, color=INK2)
        ax.grid(True, axis="x", which="major", color=GRID, lw=0.5)
        ax.tick_params(axis="x", labelsize=7, colors=INK2)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.plot([], [], "o", color=C_SK, label="skani (default)")
        ax.plot([], [], "s", color=C_AIF, label="aif")
        ax.plot([], [], "o", color=INK2, mfc="white", label="hollow: coverage < 50%")
        fig.legend(fontsize=7, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.55, 0.0))
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        fig.savefig(os.path.join(FIG, fname))
        plt.close(fig)

    # ---------- figure: budget sweep
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, key, lab in [(axes[0], "ani_mae_micro", "ANI MAE, all contigs (pp)"), (axes[1], "af_mae_micro", "AF MAE, all contigs (pp)")]:
        xs = [b for n, b, o in sweep if n == "aif"]
        ys = [o[key] for n, b, o in sweep if n == "aif"]
        ax.plot(xs, ys, "-s", color=C_AIF, lw=1.8, ms=5, label="aif")
        sx = [b for n, b, o in sweep if n != "aif"]
        sy = [o[key] for n, b, o in sweep if n != "aif"]
        ax.plot(sx, sy, "--o", color=C_SK, lw=1.8, ms=5.5, label="skani presets")
        for n, b, o in sweep:
            if n != "aif":
                ax.annotate({"skani": "default", "skani_medium": "medium", "skani_slow": "slow"}[n], (b, o[key]),
                            xytext=(4, 3), textcoords="offset points", fontsize=6.5, color=INK)
        ax.axvline(nums["skani_bpb"], color=GRID, lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("sketch size (bits per reference base)", fontsize=8, color=INK2)
        ax.set_title(lab, fontsize=8.5, color=INK)
        ax.grid(True, which="major", color=GRID, lw=0.5)
        ax.tick_params(labelsize=7, colors=INK2)
        for sp in ax.spines.values():
            sp.set_color(GRID)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_budget.eps"))
    plt.close(fig)

    # ---------- figure: error scatter vs true ANI
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, (lab, pred, col) in zip(axes, [("skani (default)", psk, C_SK), ("aif", paif, C_AIF)]):
        xs, ys = [], []
        for t in truth:
            if t["contig"] in pred:
                xs.append(100 * float(t["true_ani"]))
                ys.append(100 * (pred[t["contig"]][0] - float(t["true_ani"])))
        ax.scatter(xs, ys, s=3, color=col, alpha=1.0, linewidths=0)
        ax.axhline(0, color=INK2, lw=0.7)
        ax.set_ylim(-6, 6)
        ax.set_title(f"{lab}: {len(xs)} contigs reported", fontsize=8.5, color=INK)
        ax.set_xlabel("true ANI (%)", fontsize=8, color=INK2)
        ax.grid(True, color=GRID, lw=0.5)
        ax.tick_params(labelsize=7, colors=INK2)
        for sp in ax.spines.values():
            sp.set_color(GRID)
    axes[0].set_ylabel("estimate − truth (pp)", fontsize=8, color=INK2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_scatter.eps"))
    plt.close(fig)

    # ---------- development history (evaluated on the main benchmark during development)
    hist = load_eval(os.path.join(R, "history.json"))
    desc = OrderedDict([
        ("v1", "k=21, tiers 0.003:0.55 + 0.58:0.45, 2-state HMM fitted by EM"),
        ("v2", "k=17, one tier (0.02), 15-level HMM, calibrated sigma"),
        ("v3", "two tiers, identifiable fine grid, posterior-mean ANI"),
        ("tuned_k17", "parameters from the tuning set, fixed k=17"),
        ("final", "+ adaptive k per reference (v1.0)"),
        ("v11", "v1.1: speed rewrite (binned HMM + refinement)"),
    ])
    rows = []
    for key, d in desc.items():
        o = hist["overall"][key]
        rows.append([key.replace("_", "\\(ul"), d, f(o["coverage"], 3), f(o["ani_mae_micro"], 3), f(o["ani_mae_macro"], 3),
                     f(o["af_mae_micro"]), f(o["af_mae_macro"])])
    tbl(f"{GEN}/tab_history.tbl", ["l l c c s c s", "l l c c c c c"], "l lw(1.45i) n n n n n",
        [["version", "change", "coverage", "ANI MAE (pp)", "", "AF MAE (pp)", ""],
         ["", "", "", "all", "cond. mean", "all", "cond. mean"]], rows)

    # ---------- numbers
    for key, val in [("sk", sk), ("ai", ai), ("skr", skr), ("fsk", fsk), ("fai", fai)]:
        for m in ("coverage", "ani_mae_micro", "ani_mae_macro", "af_mae_micro", "af_mae_macro", "paired_ani_mae_macro", "paired_af_mae_macro"):
            nums[f"{key}_{m}"] = val[m]
    for name in ["skani_medium", "skani_slow", "budget_0.44", "budget_0.22", "budget_0.064"]:
        o = ov(name)
        for m in ("coverage", "ani_mae_micro", "af_mae_micro"):
            nums[f"{name}_{m}"] = o[m]
    with open(os.path.join(GEN, "numbers.txt"), "w") as fh:
        for k_, v in nums.items():
            fh.write(f"{k_} = {v}\n")
    for k_, v in nums.items():
        print(f"{k_} = {v}")


if __name__ == "__main__":
    main()
