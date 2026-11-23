#!/usr/bin/env python3
"""Generated items for the method-and-results overview (report_overview/): a method schematic and a
results-at-a-glance table, every number read from the committed result files of the other reports.
Usage: /usr/bin/python3 scripts/overview_report.py
"""
import json
import os

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "report_overview", "gen")
FIG = os.path.join(ROOT, "report_overview", "figs")
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "serif", "font.size": 7.5})


def J(*p):
    return json.load(open(os.path.join(ROOT, *p)))


def numbers(path):
    d = {}
    for line in open(os.path.join(ROOT, path)):
        if "=" in line and "\t" not in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
        elif "\t" in line:
            k, v = line.rstrip("\n").split("\t", 1)
            d[k] = v
    return d


def f(x, d=2):
    return "\\(em" if x is None else f"{float(x):.{d}f}"


def schematic():
    fig, ax = plt.subplots(figsize=(6.8, 3.3))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    def box(x, y, w, h, txt, fc="#eef4fc", ec="#2a78d6"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=0.9))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=6.6, color=INK, linespacing=1.25)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=INK2))
    # reference side
    ax.text(2, 47.5, "Reference (sketch once)", fontsize=7.5, weight="bold", color=INK)
    box(2, 37, 18, 7, "reference genome(s)\n(FASTA / gzip)", fc="#f3f2ef", ec=INK2)
    box(25, 37, 20, 7, "FracMinHash selection\nh(x) < s  (s about 0.12)")
    box(50, 36.5, 23, 8, "two Bloom tiers (hash ranges)\nFPR 0.005 and 0.1\n0.88 bits per base\nexact FPR from fill ratio")
    box(78, 36.5, 20, 8, "database (v1.5-1.6)\none k for all references\nbit-sliced screening\nmemory-mapped", fc="#eaf7f1", ec="#1baf7a")
    arrow(20.5, 40.5, 24.5, 40.5)
    arrow(45.5, 40.5, 49.5, 40.5)
    arrow(73.5, 40.5, 77.5, 40.5)
    # query side
    ax.text(2, 31.5, "Query contig (raw sequence)", fontsize=7.5, weight="bold", color=INK)
    box(2, 21, 18, 7, "contig\n(no sketch needed)", fc="#f3f2ef", ec=INK2)
    box(25, 21, 20, 7, "selected k-mers in order\n(every k-mer with h < s)")
    box(50, 21, 23, 7, "probe the tiers: hit / miss\nper position along the contig")
    box(78, 21, 20, 7, "screening -> candidates\n-> local-score filter", fc="#eaf7f1", ec="#1baf7a")
    arrow(20.5, 24.5, 24.5, 24.5)
    arrow(45.5, 24.5, 49.5, 24.5)
    arrow(61.5, 36.5, 61.5, 28.5)
    arrow(88, 36.5, 88, 28.5)
    # model and outputs
    box(14, 2, 31, 13, "multi-level HMM over binned hits\nstates: unaligned + aligned at\nidentity levels $a$ = 0.81 ... 1.0\nhit probability $a^k + (1 - a^k)\\,\\varepsilon'$")
    box(51, 2, 47, 13, "outputs per contig (and per reference hit)\n"
                          "AF: aligned posterior (optionally at identity $\\geq T$)\n"
                          "ANI $= \\hat{C}^{1/k}$,  $\\hat{C} = (Y/N - \\varepsilon')/(1 - \\varepsilon')$  (Theorem 1)\n"
                          "calibrated 95% intervals for ANI and AF\n"
                          "identity map: segments + windowed identity", fc="#fdf3ea", ec="#eb6834")
    arrow(61.5, 20.5, 38, 15.5)
    arrow(88, 20.5, 80, 15.5)
    arrow(45.5, 8.5, 50.5, 8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_method.eps"))
    plt.close(fig)


def glance():
    nb = numbers("report_bench/gen/numbers.txt")
    real = J("results_real", "db", "eval_fixed_k.json")["sets"]["all_within"]
    pareto = J("results_real", "pareto", "eval_pareto.json")["sets"]["all_within"]
    db = J("results_real", "db", "db_eval.json")["pooled"]
    dbf = J("results_real", "fresh", "db_eval_fresh.json")["pooled"]
    g1 = J("results_real", "gtdb", "gtdb_eval.json")
    g2 = J("results_real", "fresh", "gtdb_eval_set2.json")
    con = J("results_real", "gtdb", "contam.json")
    idm = J("results_real", "idmap", "idmap.json")
    vir = J("results_viral", "viral_eval.json")
    st = J("results_real", "gtdb", "search_time.json")
    rows = [
        ["simulated contigs (2,800)", "ANI / AF MAE (pp)", f"{f(nb['ai_ani_mae_micro'], 2)} / {f(nb['ai_af_mae_micro'], 2)}",
         "skani", f"{f(nb['sk_ani_mae_micro'], 2)} / {f(nb['sk_af_mae_micro'], 2)}"],
        ["real contigs, pairwise", "ANI / AF MAE (pp)", f"{f(real['aif']['ani_mae'], 2)} / {f(real['aif']['af_mae'], 2)}",
         "skani", f"{f(real['skani']['ani_mae'], 2)} / {f(real['skani']['af_mae'], 2)}"],
        ["", "coverage (homologous)", f(real["aif"]["coverage"], 3), "skani", f(real["skani"]["coverage"], 3)],
        ["exact FracMinHash, 0.88 bits/base", "ANI MAE (pp)", f(pareto["aif"]["ani_mae"], 2), "sourmash", f(pareto["sm_S73"]["ani_mae"], 2)],
        ["alignment (27 bits/base)", "ANI / AF MAE (pp)", f"{f(pareto['aif']['ani_mae'], 2)} / {f(pareto['aif']['af_mae'], 2)}",
         "minimap2", f"{f(pareto['minimap2']['ani_mae'], 2)} / {f(pareto['minimap2']['af_mae'], 2)}"],
        ["2,005-genome DB", "top-hit ANI / AF MAE", f"{f(db['aif']['ani_mae_pp'], 2)} / {f(db['aif']['af_mae_pp'], 2)}",
         "skani search", f"{f(db['skani']['ani_mae_pp'], 2)} / {f(db['skani']['af_mae_pp'], 2)}"],
        ["2,005-genome DB, fresh", "top-hit ANI / AF MAE", f"{f(dbf['aif']['ani_mae_pp'], 2)} / {f(dbf['aif']['af_mae_pp'], 2)}",
         "skani search", f"{f(dbf['skani']['ani_mae_pp'], 2)} / {f(dbf['skani']['af_mae_pp'], 2)}"],
        ["GTDB, query set 1", "species / genus", f"{f(g1['tools']['aif_0.88']['all']['species_acc'], 3)} / "
         f"{f(g1['tools']['aif_0.88']['all']['genus_acc'], 3)}", "skani", f"{f(g1['tools']['skani']['all']['species_acc'], 3)} / {f(g1['tools']['skani']['all']['genus_acc'], 3)}"],
        ["GTDB, fresh set 2", "species / genus", f"{f(g2['aif_0.88']['all']['species_acc'] if 'aif_0.88' in g2 else g2['tools']['aif_0.88']['all']['species_acc'], 3)} / "
         f"{f(g2['tools']['aif_0.88']['all']['genus_acc'], 3)}", "skani", f"{f(g2['tools']['skani']['all']['species_acc'], 3)} / {f(g2['tools']['skani']['all']['genus_acc'], 3)}"],
        ["GTDB search, 16 thr", "wall time (s)", f(st["aif_0.88"]["s"], 1), "skani", f(st["skani"]["s"], 0)],
        ["contamination (genus rule)", "recall / precision", f"{f(con['aif|genus|default|all']['recall_contigs'], 3)} / "
         f"{f(con['aif|genus|default|all']['precision_contigs'], 3)}", "skani", f"{f(con['skani|genus|default|all']['recall_contigs'], 3)} / "
         f"{f(con['skani|genus|default|all']['precision_contigs'], 3)}"],
        ["identity maps, 3 kb", "identity MAE (pp)", f(idm["aif_3000"]["all"]["mae_pp"], 2), "fastANI", f(idm["fastani_3000"]["all"]["mae_pp"], 2)],
        ["unaligned islands", "precision / recall", f"{f(idm['aif']['islands|all']['precision'], 3)} / {f(idm['aif']['islands|all']['recall'], 3)}",
         "fastANI", f"{f(idm['fastani']['islands|all']['precision'], 3)} / {f(idm['fastani']['islands|all']['recall'], 3)}"],
        ["vOTUs, INPHARED / GPD", "edge F1 (aif: AF \\(>= 88%)", f"{f(vir['inphared']['tools']['aif_t']['edges']['f1'], 3)} / {f(vir['gpd']['tools']['aif_t']['edges']['f1'], 3)}",
         "skani small", f"{f(vir['inphared']['tools']['skani_small']['edges']['f1'], 3)} / {f(vir['gpd']['tools']['skani_small']['edges']['f1'], 3)}"],
        ["", "", "", "vclust", f"{f(vir['inphared']['tools']['vclust']['edges']['f1'], 3)} / {f(vir['gpd']['tools']['vclust']['edges']['f1'], 3)}"],
    ]
    lines = [".TS", "center box tab(@);", "l l c l c", "l l n l n.", "experiment@metric@aif@comparator@value", "_"]
    lines += ["@".join(r) for r in rows]
    lines.append(".TE")
    open(os.path.join(GEN, "tab_glance.tbl"), "w").write("\n".join(lines) + "\n")


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    schematic()
    glance()
    print("wrote", os.listdir(GEN), os.listdir(FIG))


if __name__ == "__main__":
    main()
