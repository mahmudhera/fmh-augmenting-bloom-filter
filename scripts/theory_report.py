#!/usr/bin/env python3
"""Figures and simulation tables of aif_theory.pdf (report_theory/), from results_theory/{numbers,sims,hmm_boundary}.json.

Usage: /usr/bin/python3 -I scripts/theory_report.py   (after theory_numbers.py, theory_sims.py, theory_hmm.py)
"""
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_numbers as tn  # noqa: E402

ROOT = tn.ROOT
RES = os.path.join(ROOT, "results_theory")
GEN = os.path.join(ROOT, "report_theory", "gen")
FIG = os.path.join(ROOT, "report_theory", "figs")
os.makedirs(FIG, exist_ok=True)
C1, C2, C_TH = "#2a78d6", "#eb6834", "#52514e"
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
INK2, GRID = "#52514e", "#d9d8d4"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2, "legend.frameon": False})


def axes_style(ax):
    ax.grid(True, color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def sgn(x, d=0):
    v = f"{abs(x):.{d}f}"
    return ("\\(mi" + v if x < 0 and float(v) != 0 else "+" + v)


def sci(x):
    if x >= 0.01:
        return f"{x:.2g}"
    e = math.floor(math.log10(x))
    return f"{x / 10 ** e:.1f}\\(mu10\\u\\(mi{-e}\\d"


def tbl(name, fmt_header, fmt_body, header_rows, rows):
    lines = [".TS", "center box tab(@);"] + fmt_header + [fmt_body + "."] + header_rows + ["_"] + rows + [".TE"]
    open(os.path.join(GEN, name), "w").write("\n".join(lines) + "\n")


def main():
    N = json.load(open(os.path.join(RES, "numbers.json")))
    S = json.load(open(os.path.join(RES, "sims.json")))
    H = json.load(open(os.path.join(RES, "hmm_boundary.json")))
    nums = {}

    # ---- Figure: detection, theory vs aif
    det = {round(d["ani"], 2): d for d in N["detection"]}
    off = N["detection_offset"]
    fig, ax = plt.subplots(figsize=(6.2, 2.6))
    anis = [0.84, 0.86, 0.88, 0.90, 0.92, 0.95]
    for c, a in zip(SEQ, anis):
        d = det[a]
        xs = [10 ** (2.3 + i * 0.01) for i in range(230)]
        ys = [stats.norm.cdf((x * d["rate"] - off - tn.LR_MIN) / math.sqrt(x * d["var"])) for x in xs]
        ax.plot(xs, [100 * y for y in ys], color=c, lw=1.4)
        pts = [x for x in S["detect"] if abs(x["ani"] - a) < 1e-9]
        ax.plot([p["L"] for p in pts], [100 * p["reported"] for p in pts], "o", color=c, ms=4, mec="white", mew=0.5,
                label=f"{100 * a:.0f}%")
    ax.set_xscale("log")
    ax.set_xlim(200, 25000)
    ax.set_ylim(-3, 103)
    ax.set_xlabel("contig length (bp)")
    ax.set_ylabel("contigs reported (%)")
    ax.legend(title="ANI", loc="lower right", fontsize=7, title_fontsize=7, ncol=2)
    axes_style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_detect.eps"))
    plt.close(fig)

    # ---- Figure: absorption of unaligned stretches
    A = N["absorption"]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.4), sharey=True)
    for ax, kind, key, title in ((axes[0], "int", "p_int", "inside an alignment (3 kb + gap + 3 kb)"),
                                 (axes[1], "edge", "p_edge", "at a contig end (3 kb + flank)")):
        kap = A["kappa_int"] if kind == "int" else A["kappa_edge"]
        xs = [10 ** (1.9 + i * 0.01) for i in range(185)]
        ys = []
        for G in xs:
            o = kap * sum(math.exp(-G * (x["kl_rate"] + tn.RHO_H)) for x in A["per_level"])
            ys.append(100 * o / (1 + o))
        ax.plot(xs, ys, color=C_TH, lw=1.2, ls="--", label="theory (17)")
        for mode, col, lab in (("v16", C1, "aif, all aligned mass (v1.6)"), ("v17", C2, "aif, local-identity AF (v1.7)")):
            pts = [x for x in S["absorb"] if x["mode"] == mode and x["kind"] == kind]
            ax.plot([p["G"] for p in pts], [100 * p["absorbed_share"] for p in pts], "o-", color=col, ms=4, lw=1.0,
                    mec="white", mew=0.5, label=lab)
        ax.axhline(0, color=INK2, lw=0.5)
        ax.set_xscale("log")
        ax.set_xlabel("length of the unaligned stretch (bp)")
        ax.set_title(title, fontsize=8, color=INK2)
        axes_style(ax)
    axes[0].set_ylabel("counted as aligned (%)")
    axes[1].legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_absorb.eps"))
    plt.close(fig)

    # ---- Table: null likelihood ratios
    nl = S["null"]
    rows = []
    for lam in ["1", "2", "3", "5", "7", "10"]:
        rows.append(f"{lam}@{sci(math.exp(-int(lam)))}@" + "@".join(
            (sci(nl[L]['frac_ge'][lam]) if nl[L]['frac_ge'][lam] > 0 else f"0 (< {sci(1 / nl[L]['n'])})") for L in ("1000", "5000", "20000")))
    tbl("tab_null.tbl", ["c c c s s", "c c c c c"], "n n n n n",
        ["$lambda$@bound@observed share of unrelated contigs with $LR >= e sup lambda$",
         "@$e sup {- lambda}$@1 kb (40,000)@5 kb (20,000)@20 kb (5,000)"], rows)
    nums["null_max_lr"] = max(nl[L]["max_lr"] for L in nl)
    nums["null_n"] = sum(nl[L]["n"] for L in nl)

    # ---- Table: boundary resolution (idealised, exact HMM, aif)
    bnd = {round(b["ani"], 2): b for b in N["boundary"]}
    hb = {round(h["ani"], 2): h for h in H}
    sb = {(b["mode"], round(b["ani"], 2)): b for b in S["boundary"]}
    rows = []
    for a in (0.90, 0.95, 0.99):
        b, h = bnd[a], hb[a]
        fm = lambda m, s: f"{sgn(m)} / {s:.0f}"
        rows.append(f"{100 * a:.0f}@{b['scale_bp']:.0f}@{b['mc_sd_bp']:.0f}@{fm(h['two_state']['mean'], h['two_state']['sd'])}"
                    f"@{fm(h['all_levels']['mean'], h['all_levels']['sd'])}@{fm(h['local']['mean'], h['local']['sd'])}"
                    f"@{fm(sb[('v16', a)]['mean_bp'], sb[('v16', a)]['sd_bp'])}@{fm(sb[('v17', a)]['mean_bp'], sb[('v17', a)]['sd_bp'])}")
    tbl("tab_boundary.tbl", ["c c c c s s c s", "c c c c c c c c", "c c c c c c c c"], "n n n c c c c c",
        ["ANI@scale@independent@substitution model, exact HMM@aif",
         "(%)@(16)@probes, SD@two states@all levels,@floor rule,@all mass@v1.7",
         "@(bp)@(bp)@(true level)@all mass@no windows@(v1.6)@"], rows)

    # ---- Table: spurious-match rates
    sg = S["sigma"]
    names = {"random 5 Mbp": "random sequence", "NC_000913.3": "\\f[I]E. coli\\f[R] K-12", "NC_000964.3": "\\f[I]B. subtilis\\f[R] 168",
             "NC_002516.2": "\\f[I]P. aeruginosa\\f[R] PAO1", "NC_007795.1": "\\f[I]S. aureus\\f[R] NCTC 8325"}
    rows = []
    for g in names:
        r = [x for x in sg if x["genome"] == g]
        rows.append(f"{names[g]}@{r[0]['bases'] / 1e6:.2f}@" + "@".join(
            f"{1e4 * x['analytic']:.2f} / {1e4 * x['emp']:.1f}" for x in sorted(r, key=lambda x: x["k"])))
    tbl("tab_sigma.tbl", ["l c c s s s s", "l c c c c c c"], "l n c c c c c",
        ["genome@Mbp@$sigma$ (\\(mu10\\u\\(mi4\\d): random-sequence value / calibrated",
         "@@$k=17$@18@19@20@21"], rows)

    # ---- Table: heterogeneous contigs
    rows = []
    for x in S["het"]:
        rows.append(f"{100 * x['a1']:.1f} / {100 * x['a2']:.1f}@{100 * x['true']:.2f}@{100 * x['power_mean']:.2f}@{100 * x['pooled']:.2f}"
                    f"@{100 * x['auto']:.2f}@{100 * x['het_called']:.1f}")
    tbl("tab_het.tbl", ["c c c c c c", "c c c c c c"], "c n n n n n",
        ["halves (ANI %)@true ANI@$M sub k$ (20)@pooled, aif@reported, aif@heterogeneous",
         "@(%)@(%)@\\f[CR]--het 0\\f[R] (%)@(default) (%)@called (%)"], rows)

    # ---- Table: indels
    rows = []
    for x in S["indel"]:
        rows.append(f"{100 * x['r']:.1f}@{100 * x['excl']:.2f}@{100 * x['incl']:.2f}@{100 * x['pred']:.2f}@{100 * x['aif']:.2f}@{100 * x['af']:.1f}")
    tbl("tab_indel.tbl", ["c c c c c c", "c c c c c c"], "n n n n n n",
        ["indel rate@ANI, indel columns@ANI, indels as@$(1-p)(1-r)$@aif ANI@aif AF",
         "$r$ (%)@excluded (%)@mismatches (%)@(%)@(%)@(%)"], rows)

    json.dump(nums, open(os.path.join(RES, "report_numbers.json"), "w"), indent=1)
    print("figures and tables written")


if __name__ == "__main__":
    main()
