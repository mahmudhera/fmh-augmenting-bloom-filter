#!/usr/bin/env python3
"""Tables, figures and numbers for the viral clustering report (report_viral/).

Reads data/viral (bench_real/viral_bench.py outputs) and writes report_viral/gen/*.tbl,
report_viral/gen/numbers.txt and report_viral/figs/*.eps; copies summaries to results_viral/.
Usage: /usr/bin/python3 scripts/viral_report.py
"""
import csv
import json
import os
import sys
from collections import OrderedDict

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench_real"))
import viral_bench as vb  # noqa: E402

V = os.path.join(ROOT, "data", "viral")
GEN = os.path.join(ROOT, "report_viral", "gen")
FIG = os.path.join(ROOT, "report_viral", "figs")
DS = [("inphared", "INPHARED"), ("gpd", "GPD")]
TOOLS = [("aif", "aif (default AF)"), ("aif_t", "aif (AF at identity \\(>= 88%)"), ("skani", "skani"),
         ("skani_small", "skani --small-genomes"), ("vclust", "vclust (alignment)")]
PLAIN = {"aif": "aif (default AF)", "aif_t": "aif (AF at identity \u2265 88%)", "skani": "skani", "skani_small": "skani --small-genomes",
         "vclust": "vclust"}
COL = {"aif": "#8fb8ea", "aif_t": "#2a78d6", "skani": "#eb6834", "skani_small": "#eda100", "vclust": "#1baf7a"}
MK = {"aif": "o", "aif_t": "o", "skani": "s", "skani_small": "s", "vclust": "D"}
INK2, GRID = "#52514e", "#d9d8d4"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2})


def tbl(path, fmt_header, fmt_body, header_rows, rows, opts="center box tab(@)"):
    lines = [".TS", opts + ";"] + fmt_header + [fmt_body + "."]
    for i, hr in enumerate(header_rows):
        fmt = fmt_header[i].split() if i < len(fmt_header) else []
        if fmt and len(fmt) == len(hr):
            hr = [h for h, f_ in zip(hr, fmt) if f_ != "s"]
        lines.append("@".join(hr))
    lines.append("_")
    for r in rows:
        lines.append("_" if r == "_" else "@".join(r))
    lines.append(".TE")
    open(os.path.join(GEN, path), "w").write("\n".join(lines) + "\n")


def f(x, d=2):
    return "\\(em" if x is None or x != x else f"{x:.{d}f}"


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    E = json.load(open(os.path.join(V, "viral_eval.json")))
    TT = json.load(open(os.path.join(V, "results_tool_time.json")))
    BT = json.load(open(os.path.join(V, "results_blast_time.json")))
    TU = json.load(open(os.path.join(V, "tune_af_min_ident.json")))
    nums = OrderedDict()
    # (1) data
    rows = []
    for ds, lab in DS:
        r = E[ds]
        L = vb.lengths(ds)
        med = sorted(L.values())[len(L) // 2]
        T = list(csv.DictReader(open(os.path.join(V, ds, "truth.tsv")), delimiter="\t"))
        n80 = sum(1 for x in T if float(x["ani"]) >= 0.8 and float(x["af_short"]) >= 0.1)
        multi = r["tools"]["aif"]["clusters"]["n_truth_multi"]
        rows.append([lab, f"{r['n_seq']:,}", f(r["mbp"], 0), f(med / 1e3, 1), f"{r['truth_pairs']:,}", f"{n80:,}", f"{r['truth_edges']:,}",
                     f"{r['truth_votus']:,}", f"{multi:,}", f(BT[ds]["blastn_s"], 0)])
        nums[f"data_{ds}"] = (r["n_seq"], r["mbp"], med, r["truth_pairs"], n80, r["truth_edges"], r["truth_votus"], multi, BT[ds]["blastn_s"])
    tbl("tab_data.tbl", ["l c c c c c c c c c", "^ ^ ^ ^ ^ c c c c c"], "l n n n n n n n n n",
        [["data set", "sequences", "Mbp", "median", "BLAST", "pairs", "vOTU", "vOTUs", "multi-member", "BLAST"],
         ["", "", "", "(kb)", "pairs", "\\(>= 80% ANI", "edges", "", "vOTUs", "(s, 56 thr)"]], rows)
    # (2) pair accuracy
    rows = []
    for ds, lab in DS:
        for i, (t, tl) in enumerate(TOOLS):
            ev = E[ds]["tools"][t]
            a, nt = ev["all"], ev["near_threshold"]
            rows.append([lab if i == 0 else "", tl, f(a["reported"], 3), f(a["ani_mae_pp"], 2), f(a["af_mae_pp"], 2), f(nt["reported"], 3),
                         f(nt["ani_mae_pp"], 3), f(nt["af_mae_pp"], 2)])
            nums[f"pairs_{ds}_{t}"] = (a["reported"], a["ani_mae_pp"], a["af_mae_pp"], nt["reported"], nt["ani_mae_pp"], nt["af_mae_pp"], nt["n"], a["n"])
        rows.append("_")
    tbl("tab_pairs.tbl", ["l l c s s c s s", "^ ^ c c c c c c"], "l l n n n n n n",
        [["data set", "tool", "pairs \\(>= 80% ANI, AF \\(>= 0.1", "", "", "near the thresholds", "", ""],
         ["", "", "reported", "ANI MAE", "AF MAE", "reported", "ANI MAE", "AF MAE"]], rows[:-1])
    # (3) ANI error by identity stratum
    rows = []
    for ds, lab in DS:
        for i, (t, tl) in enumerate(TOOLS):
            ev = E[ds]["tools"][t]
            r = [lab if i == 0 else "", tl]
            for st in ("80-90", "90-95", "95-98", ">=98"):
                r += [f(ev[st]["reported"], 3), f(ev[st]["ani_mae_pp"], 2)]
            rows.append(r)
        rows.append("_")
    tbl("tab_strata.tbl", ["l l c s c s c s c s", "^ ^ c c c c c c c c"], "l l n n n n n n n n",
        [["data set", "tool", "ANI 80\\(en90%", "", "90\\(en95%", "", "95\\(en98%", "", "\\(>= 98%", ""],
         ["", "", "rep.", "MAE", "rep.", "MAE", "rep.", "MAE", "rep.", "MAE"]], rows[:-1])
    # (4) vOTU edges and clusters
    rows = []
    for ds, lab in DS:
        for i, (t, tl) in enumerate(TOOLS):
            ev = E[ds]["tools"][t]
            e, c = ev["edges"], ev["clusters"]
            rows.append([lab if i == 0 else "", tl, f"{e['predicted']:,}", f(e["precision"], 3), f(e["recall"], 3), f(e["f1"], 3), f"{c['n_votus']:,}",
                         f(c["ari"], 3), f(c["identical_cluster_frac"], 3)])
            nums[f"votu_{ds}_{t}"] = (e["precision"], e["recall"], e["f1"], c["n_votus"], c["ari"], c["identical_cluster_frac"])
        rows.append("_")
    tbl("tab_votu.tbl", ["l l c c c c c c c", "^ ^ c c c c c c c"], "l l n n n n n n n",
        [["data set", "tool", "vOTU edges", "precision", "recall", "F1", "vOTUs", "ARI", "identical"],
         ["", "", "predicted", "", "", "", "", "", "clusters"]], rows[:-1])
    # (5) resources
    rows = []
    for ds, lab in DS:
        t = TT[ds]
        recs = [("aif (default AF)", t["aif_build"], t["aif_search"], t["aif_db_bytes"]),
                ("aif (AF at identity \\(>= 88%)", t["aif_build"], t["aif_t_search"], t["aif_db_bytes"]),
                ("skani", None, t["skani"], None), ("skani --small-genomes", None, t["skani_small"], None),
                ("vclust (alignment)", t["vclust_prefilter"], t["vclust_align"], None)]
        for i, (tl, b, s_, db) in enumerate(recs):
            tot = (b[0] if b else 0) + s_[0]
            mem = max(b[1] if b else 0, s_[1]) / 1e6
            rows.append([lab if i == 0 else "", tl, f(b[0], 1) if b else "\\(em", f(s_[0], 1), f(tot, 1), f(mem, 2), f(db / 1e6, 0) if db else "\\(em"])
            nums[f"time_{ds}_{tl.split()[0]}{'_t' if 'identity' in tl else ''}{'_small' if 'small' in tl else ''}"] = (tot, mem)
        rows.append("_")
    tbl("tab_time.tbl", ["l l c c c c c", "^ ^ c c c c c"], "l l n n n n n",
        [["data set", "tool", "index / prefilter", "search / align", "total", "peak memory", "database"],
         ["", "", "(s)", "(s)", "(s)", "(GB)", "(MB)"]], rows[:-1])
    # (6) tuning grid
    rows = []
    for tau in ["0", "0.86", "0.88", "0.90", "0.92"]:
        r = TU[tau]
        rows.append([("off" if tau == "0" else tau) + (" (chosen)" if tau == TU["chosen"] else ""), f(r["af_mae_near_pp"], 2), f(r["af_mae_all_pp"], 2),
                     f(r["edge_precision"], 3), f(r["edge_recall"], 3), f(r["edge_f1"], 4)])
    tbl("tab_tune.tbl", ["l c c c c c", "^ c c c c c"], "l n n n n n",
        [["--af-min-ident", "AF MAE", "AF MAE", "edge", "edge", "edge"], ["", "near thresholds", "all pairs", "precision", "recall", "F1"]], rows)
    nums["tune"] = TU
    # figure 1: AF bias against truth AF
    bins = [(0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 0.95), (0.95, 1.01)]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5), sharey=True)
    hs = []
    for ax, (ds, lab) in zip(axes, DS):
        L = vb.lengths(ds)
        T = {(r["short"], r["long"]): (float(r["ani"]), float(r["af_short"])) for r in csv.DictReader(open(os.path.join(V, ds, "truth.tsv")), delimiter="\t")}
        for t, tl in TOOLS:
            P = vb.load_tool(ds, t, L)
            xs, ys = [], []
            for lo, hi in bins:
                ks = [k for k, (a, f_) in T.items() if a >= 0.8 and lo <= f_ < hi and k in P]
                if len(ks) < 30:
                    continue
                xs.append((lo + min(hi, 1.0)) / 2)
                ys.append(100 * sum(P[k][1] - T[k][1] for k in ks) / len(ks))
            h, = ax.plot(xs, ys, "-" + MK[t], color=COL[t], lw=1.4, ms=3.5, mfc=COL[t] if t != "aif" else "white")
            if ds == "inphared":
                hs.append((h, PLAIN[t]))
        ax.axhline(0, color=INK2, lw=0.6)
        ax.axvline(0.85, color=INK2, lw=0.5, ls=":")
        ax.set_title(lab, fontsize=8)
        ax.set_xlabel("BLAST AF of the shorter sequence", fontsize=7.5)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("AF bias, tool \u2212 BLAST (pp)", fontsize=7.5)
    fig.legend([h for h, _ in hs], [l_ for _, l_ in hs], loc="upper center", ncol=5, frameon=False, fontsize=6.3)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(FIG, "fig_afbias.eps"))
    plt.close(fig)
    # figure 2: clustering agreement against run time
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4))
    for ax, key, ylab in ((axes[0], "f1", "vOTU-edge F1"), (axes[1], "ari", "ARI vs BLAST vOTUs")):
        for ds, lab in DS:
            t = TT[ds]
            tm = {"aif": t["aif_build"][0] + t["aif_search"][0], "aif_t": t["aif_build"][0] + t["aif_t_search"][0], "skani": t["skani"][0],
                  "skani_small": t["skani_small"][0], "vclust": t["vclust_prefilter"][0] + t["vclust_align"][0]}
            for tl, _ in TOOLS:
                ev = E[ds]["tools"][tl]
                y = ev["edges"]["f1"] if key == "f1" else ev["clusters"]["ari"]
                ax.plot([tm[tl]], [y], MK[tl], color=COL[tl], ms=6 if ds == "inphared" else 5, mfc=COL[tl] if ds == "inphared" else "white", mew=1.3)
        ax.set_xscale("log")
        ax.set_xlabel("all-vs-all wall time, 16 threads (s)", fontsize=7.5)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.grid(True, which="both", color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    hh = [plt.Line2D([], [], color=COL[t], marker=MK[t], ls="none", ms=5) for t, _ in TOOLS]
    hh += [plt.Line2D([], [], color=INK2, marker="o", ls="none", ms=5), plt.Line2D([], [], color=INK2, marker="o", mfc="white", ls="none", ms=5)]
    fig.legend(hh, [PLAIN[t] for t, _ in TOOLS] + ["INPHARED (filled)", "GPD (open)"], loc="upper center", ncol=4, frameon=False, fontsize=6.3)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(os.path.join(FIG, "fig_time.eps"))
    plt.close(fig)
    with open(os.path.join(GEN, "numbers.txt"), "w") as fh:
        for k, v in nums.items():
            fh.write(f"{k}\t{v}\n")
    out = os.path.join(ROOT, "results_viral")
    os.makedirs(out, exist_ok=True)
    for fn in ("viral_eval.json", "results_tool_time.json", "results_blast_time.json", "tune_af_min_ident.json"):
        open(os.path.join(out, fn), "w").write(open(os.path.join(V, fn)).read())
    for ds in ("inphared", "gpd", "tune"):
        open(os.path.join(out, f"{ds}_ids.txt"), "w").write("".join(n + "\n" for n in vb.lengths(ds)))
    print("\n".join(f"{k} = {v}" for k, v in nums.items() if k != "tune"))


if __name__ == "__main__":
    main()
