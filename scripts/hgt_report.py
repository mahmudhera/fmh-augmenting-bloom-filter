#!/usr/bin/env python3
"""Tables and figures of the HGT atlas report (report_hgt/ -> aif_hgt_atlas_report.pdf).

Reads results_hgt/complete_complete/{summary,ecology,annotation}.json, results_hgt/complete/validate_*.json,
data/hgt/complete/search_time*.json and the validation files in data/hgt/validate_complete/.
Usage: /usr/bin/python3 -I scripts/hgt_report.py
"""
import collections
import csv
import gzip
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench_real"))
from hgt_analyze import lca_class  # noqa: E402

RES = os.path.join(ROOT, "results_hgt")
HD = os.path.join(ROOT, "data", "hgt")
GEN = os.path.join(ROOT, "report_hgt", "gen")
FIG = os.path.join(ROOT, "report_hgt", "figs")
os.makedirs(GEN, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
C1, C2 = "#2a78d6", "#eb6834"
INK2, GRID = "#52514e", "#d9d8d4"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2, "legend.frameon": False})
DIST = ["same genus", "same family", "same order", "same class", "same phylum", "same domain", "different domain"]
DLAB = ["genus", "family", "order", "class", "phylum", "domain", "none"]
HABS = ["gut", "host_other", "food", "engineered", "marine", "sediment", "freshwater", "subsurface", "extreme", "soil",
        "plant"]
HLAB = {"gut": "gut", "host_other": "host, other sites", "food": "food", "engineered": "engineered", "marine": "marine",
        "sediment": "sediment", "freshwater": "freshwater", "subsurface": "subsurface", "extreme": "extreme",
        "soil": "soil", "plant": "plant"}


def axes_style(ax):
    ax.grid(True, color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def tbl(name, fmt_header, fmt_body, header_rows, rows):
    lines = [".TS", "center box tab(@);"] + fmt_header + [fmt_body + "."] + header_rows + ["_"] + rows + [".TE"]
    open(os.path.join(GEN, name), "w").write("\n".join(lines) + "\n")


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def genomes():
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    for g in G.values():
        g["tax"] = [x[3:] for x in g["taxonomy"].split(";")]
    return G


def pair_level(G):
    """genome-pair recall/precision against BLAST (pairs with an HSP >= 99% over >= 1 kb), cross-genus"""
    vd = os.path.join(HD, "validate_complete")
    U = {l[1:].split("|")[0] for l in open(os.path.join(vd, "U.fa")) if l[0] == ">"}
    Rb = {l[1:].split("|")[0] for l in open(os.path.join(vd, "blast_sub_windows.fa")) if l[0] == ">"}
    best = collections.defaultdict(int)
    for f in os.listdir(os.path.join(vd, "blast")):
        if f.endswith(".tsv"):
            for l in open(os.path.join(vd, "blast", f)):
                q, sid, pi, le, qs, qe = l.rstrip("\n").split("\t")
                a, t = q.split("__")[0], sid.split("__")[0]
                if t != a and float(pi) >= 99:
                    best[(a, t)] = max(best[(a, t)], int(le))
    T = {k for k, v in best.items() if v >= 1000 and lca_class(G[k[0]]["tax"], G[k[1]]["tax"]) >= 1}
    out = {}
    for f in sorted(os.listdir(vd)):
        if f.startswith("aif_w5000") and f.endswith(".tsv") and "_m" not in f:
            agg = collections.defaultdict(float)
            for l in open(os.path.join(vd, f)):
                w, ref, rank, ani, af, th, lr = l.rstrip("\n").split("\t")
                a, ci, st, ln = w.split("|")
                if a in Rb and ref in U and float(ani) >= 0.99 and float(af) * int(ln) >= 500:
                    agg[(a, ref)] = max(agg[(a, ref)], float(af) * int(ln))
            C = {k for k, v in agg.items() if v >= 1000 and lca_class(G[k[0]]["tax"], G[k[1]]["tax"]) >= 1}
            out[f[4:-4]] = {"truth": len(T), "calls": len(C), "both": len(C & T), "recall": len(C & T) / max(1, len(T)),
                            "precision": len(C & T) / max(1, len(C))}
    return out


FULL = "all_flankcontig_reciprocal"     # full atlas: flank on the contig + found in both directions
FULLU = "all_flankuniform_reciprocal"  # sensitivity: and no flank placed on the partner's side, in any genome
SUBLAB = {"isolate": "isolates", "isolate_xstudy": "isolates, other projects", "mag": "MAGs",
          "mag_xstudy": "MAGs, other projects"}
SUBHEAD = {"isolate": ("isolates", ""), "isolate_xstudy": ("isolates,", "other projects"), "mag": ("MAGs", ""),
           "mag_xstudy": ("MAGs,", "other projects")}


def full_atlas(G, N, A0):
    """tables and figure of Section 7 (all 199,923 representatives)"""
    F = {}
    SF = json.load(open(os.path.join(RES, FULL, "summary.json")))
    SU = json.load(open(os.path.join(RES, FULLU, "summary.json")))
    S1 = json.load(open(os.path.join(RES, "all", "summary.json")))
    CAL = json.load(open(os.path.join(RES, "all", "flank_calibration.json")))
    AF = json.load(open(os.path.join(RES, FULL, "annotation.json")))
    ST = json.load(open(os.path.join(HD, "all", "search_time.json")))
    FST = json.load(open(os.path.join(HD, "all", "flanks", "search_time.json")))
    FEX = json.load(open(os.path.join(HD, "all", "flanks", "extract.json")))
    F["search"] = {"genomes": SF["genomes"], "windows": SF["windows"], "parts": len(ST),
                   "aif_s": sum(v["aif_s"] for v in ST.values()), "peak_gb": max(v["kb"] for v in ST.values()) / 1e6,
                   "flank_s": sum(v["aif_s"] for v in FST.values()), **{"flank_" + k: v for k, v in FEX.items()}}
    # flank evidence by genome category
    cats = [("complete", "complete isolates"), ("draft isolate", "draft isolates"), ("MAG or SAG", "MAGs and SAGs")]
    ev = [("no flank", "none: the block is a whole contig"), ("host lineage", "nearer the host's lineage"),
          ("neither", "equally distant"), ("unplaced", "no relative at $>=$ 80%"),
          ("partner lineage", "nearer the partner's lineage"), ("partner itself", "the partner genome itself")]
    rows = [f"{lab}@" + "@".join(pct(CAL[c].get(k, 0) / CAL[c]["blocks"]) for c, _ in cats) for k, lab in ev]
    rows.append("blocks@" + "@".join(f"{CAL[c]['blocks']:,}" for c, _ in cats))
    tbl("tab_full_flank.tbl", ["l c c c"], "l n n n", ["closest relative of the flanks@" + "@".join(l for _, l in cats)],
        rows)
    F["calibration"] = CAL
    # extent by category, main rule and sensitivity
    order = ["complete isolate", "draft isolate", "MAG", "SAG"]
    rows = []
    for c in order:
        a, u = SF["by_category"][c], SU["by_category"][c]
        rows.append(f"{c}s@{a['genomes']:,}@{pct(a['with_cross_genus'] / a['genomes'])}@"
                    f"{pct(a['with_cross_phylum'] / a['genomes'])}@{pct(u['with_cross_genus'] / u['genomes'])}@"
                    f"{pct(u['with_cross_phylum'] / u['genomes'])}")
    rows.append(f"all@{SF['genomes']:,}@{pct(SF['genomes_with_cross_genus'] / SF['genomes'])}@"
                f"{pct(SF['genomes_with_cross_phylum'] / SF['genomes'])}@{pct(SU['genomes_with_cross_genus'] / SU['genomes'])}@"
                f"{pct(SU['genomes_with_cross_phylum'] / SU['genomes'])}")
    tbl("tab_full_cat.tbl", ["c c c s c s", "c c c c c c"], "l n n n n n",
        ["genomes@number@main rule (%)@with partner test (%)", "@@other genus@other phylum@other genus@other phylum"], rows)
    # pairs by distance
    rows = [f"{d}@{S1['pairs_by_distance'][d]:,}@{SF['pairs_by_distance'][d]:,}@{SU['pairs_by_distance'][d]:,}" for d in DIST]
    tbl("tab_full_dist.tbl", ["c c c c", "c c c c"], "l n n n",
        ["lowest common taxon@one direction@both directions@both directions,", "@@(main rule)@partner test"], rows)
    F["summary"] = {k: SF[k] for k in SF if k not in ("ecology",)}
    F["summary_uniform"] = {k: SU[k] for k in SU if k not in ("ecology",)}
    F["summary_one_direction"] = {k: S1[k] for k in ("pairs_by_distance", "cross_genus_pairs", "genomes_with_cross_genus",
                                                     "genomes_with_cross_phylum", "anchored_pairs", "windows_with_hits",
                                                     "windows_saturated")}
    # phyla
    xg, xp = set(), set()
    for r in csv.DictReader(gzip.open(os.path.join(RES, FULL, "pairs.tsv.gz"), "rt"), delimiter="\t"):
        if r["distance"] != "same genus":
            xg.add(r["query"])
        if r["distance"] in ("same domain", "different domain"):
            xp.add(r["query"])
    ph = collections.defaultdict(collections.Counter)
    for a, g in G.items():
        c = ph[g["tax"][1]]
        c["n"] += 1
        c["iso"] += g["category"] == "isolate"
        c["xg"] += a in xg
        c["xp"] += a in xp
    top = sorted(ph.items(), key=lambda x: -x[1]["n"])
    sel = [x for x in top if x[1]["n"] >= 1000][:16]
    sel += [x for x in top if x[1]["n"] >= 300 and x[1]["iso"] / x[1]["n"] < 0.01 and x not in sel][:4]
    rows = [f"{p}@{c['n']:,}@{pct(c['iso'] / c['n'])}@{pct(c['xg'] / c['n'])}@"
            f"{pct(c['xp'] / c['n'])}" for p, c in sel]
    tbl("tab_full_phyla.tbl", ["c c c c c"], "l n n n n",
        ["phylum@genomes@isolates (%)@other genus (%)@other phylum (%)"], rows)
    F["phyla"] = {p: dict(c) for p, c in top}
    # ecology subsets
    E = {}
    for k in ("isolate", "isolate_xstudy", "mag", "mag_xstudy"):
        f = os.path.join(RES, FULL, f"ecology_{k}.json")
        if os.path.exists(f):
            E[k] = json.load(open(f))
    Ecc = json.load(open(os.path.join(RES, "complete_complete", "ecology.json")))
    rows = []
    def eco_row(lab, e):
        r = [lab, f"{e['genomes']:,}", f"{e['exchanging_pairs']:,}"]
        for d in DIST[:-1]:
            v = e["strata"][d]
            r.append(f"{v['ratio']:.1f}" if v.get("ratio") else "\\(em")
        ci = e["mh_ci"]["cross_genus"]
        r.append(f"{e['mh_odds_cross_genus']:.1f} ({ci[0]:.1f}\\(en{ci[1]:.1f})")
        return "@".join(r)
    rows.append(eco_row("complete genomes", Ecc))
    for k, e in E.items():
        rows.append(eco_row(SUBLAB[k], e))
    tbl("tab_full_eco.tbl", ["c c c c s s s s s c", "c c c c c c c c c c"], "l n n n n n n n n c",
        ["subset@genomes@exchanging@same / different habitat, by lowest common taxon@MH odds ratio",
         "@@pairs@genus@family@order@class@phylum@domain@between genera"], rows)
    F["ecology"] = {k: {kk: e[kk] for kk in ("genomes", "habitats", "exchanging_pairs", "mh_odds_cross_genus",
                                             "mh_odds_cross_family", "mh_odds_same_genus", "mh_ci")}
                    | {"ratios": {d: e["strata"][d].get("ratio") for d in DIST}} for k, e in E.items()}
    hp = ["gut|gut", "host_other|host_other", "engineered|engineered", "gut|host_other", "food|engineered",
          "gut|food", "gut|engineered", "marine|marine", "freshwater|freshwater", "soil|soil", "plant|plant",
          "sediment|sediment", "extreme|extreme"]
    cols = [("complete", Ecc)] + [(k, e) for k, e in E.items()]
    rows = []
    for h in hp:
        a, b = h.split("|")
        lab = HLAB[a] if a == b else f"{HLAB[a]} \\(em {HLAB[b]}"
        cells = []
        for _, e in cols:
            v = e["habitat_pairs_cross_order"].get(h) or e["habitat_pairs_cross_order"].get(f"{b}|{a}")
            cells.append(f"{v['ratio']:.2f}" if v and v["expected"] >= 2 else "\\(em")
        rows.append(lab + "@" + "@".join(cells))
    tbl("tab_full_hab.tbl", ["c c " + " ".join("c" for _ in cols[1:])] * 2, "l n " + " ".join("n" for _ in cols[1:]),
        ["habitats@complete@" + "@".join(SUBHEAD[k][0] for k, _ in cols[1:]),
         "@genomes@" + "@".join(SUBHEAD[k][1] for k, _ in cols[1:])], rows)
    F["habitats"] = {k: {h: e["habitat_pairs_cross_order"].get(h) for h in hp} for k, e in cols}
    if "isolate" in E:
        heatmap(E["isolate"]["habitat_pairs_cross_order"], "fig_full_habitats.eps")
    if os.path.exists(os.path.join(RES, FULL, "magcheck.json")):
        provenance(F, E, eco_row)
    # gene classes, full atlas against the complete-genome atlas
    c0 = {t["class"]: t for t in A0["classes"]}
    rows = []
    for t in AF["classes"]:
        z = c0[t["class"]]
        rows.append(f"{t['class']}@{pct(t['shared_share'])}@{pct(t['background_share'])}@{t['odds_ratio']:.2f}@"
                    f"{z['odds_ratio']:.2f}")
    tbl("tab_full_func.tbl", ["c c s c c", "c c c c c"], "l n n n n",
        ["gene class@full atlas, genes (%)@odds ratio@odds ratio", "@shared@background@full atlas@complete genomes"], rows)
    F["annotation"] = {k: AF[k] for k in ("regions", "region_bp", "genomes", "genes_shared", "genes_background",
                                          "by_distance")}
    F["annotation_classes"] = AF["classes"]
    F["annotation_top_kos"] = AF["top_kos"][:25]
    N["full"] = F


def provenance(F, E, eco_row):
    """Section 7.5: MAG source samples, co-sample exchange, controls, and the gene content of shared DNA"""
    PR = json.load(open(os.path.join(RES, "provenance.json")))
    MC = json.load(open(os.path.join(RES, FULL, "magcheck.json")))

    def load(tag, k):
        return json.load(open(os.path.join(RES, tag, f"ecology_{k}.json")))
    ks = load(FULL, "mag_known_xsample")
    # co-sample pairs against all other pairs (known-source MAGs), strata above the genus
    def sig(x):
        return f"{x:,.0f}" if x >= 100 else f"{x:.1f}" if x >= 1 else f"{x:.2f}"
    cs = ks["cosample"]
    rows = ["pairs sharing a source sample@" + "@".join(f"{cs[d]['pairs']:,}" for d in DIST[1:]),
            "\\h'0.8m'exchanging@" + "@".join(f"{cs[d]['exchanging']:,}" for d in DIST[1:]),
            "exchanging per million pairs:@" + "@" * 5,
            "\\h'0.8m'pairs sharing a sample@" + "@".join(sig(1e6 * cs[d]["rate"]) for d in DIST[1:]),
            "\\h'0.8m'all other pairs@" + "@".join(sig(1e6 * cs[d]["rate_other_pairs"]) for d in DIST[1:]),
            "\\h'0.8m'ratio@" + "@".join(f"{cs[d]['rate'] / cs[d]['rate_other_pairs']:.0f}" for d in DIST[1:])]
    tbl("tab_full_cosample.tbl", ["c c s s s s s", "c c c c c c c"], "l n n n n n n",
        ["MAGs with a known source@lowest common taxon", "@family@order@class@phylum@domain@none"], rows)
    # ecology under provenance and content controls
    M = FULL + "_mobile"
    sets = [("MAGs", E["mag"]), ("\\h'0.8m'known source", load(FULL, "mag_known")),
            ("\\h'1.6m'other samples", ks), ("\\h'1.6m'other projects", load(FULL, "mag_known_xstudy")),
            ("\\h'1.6m'other projects, samples", load(FULL, "mag_known_xstudy_xsample")),
            ("mobile gene: isolates", load(M, "isolate")),
            ("\\h'1.6m'other projects", load(M, "isolate_xstudy")),
            ("mobile gene: MAGs", load(M, "mag")),
            ("\\h'1.6m'all controls", load(M, "mag_known_xstudy_xsample"))]
    tbl("tab_full_prov.tbl", ["c c c c s s s s s c", "c c c c c c c c c c"], "l n n n n n n n n c",
        ["subset@genomes@exchanging@same / different habitat, by lowest common taxon@MH odds ratio",
         "@@pairs@genus@family@order@class@phylum@domain@between genera"], [eco_row(l, e) for l, e in sets])
    # gene content of the shared regions by kind of pair (between genera)
    sh, bg, ps = MC["shared"]["between genera"], MC["background"], MC["pairs"]
    kinds = [("background, isolates", bg["isolate"], None), ("background, MAGs", bg["MAG"], None),
             ("isolates", sh["isolates"], ps["isolates"]),
             ("MAGs, other study, other habitat", sh["MAG, other study, different habitat"],
              ps["MAG, other study, different habitat"]),
             ("MAGs, other study, habitat unknown", sh["MAG, other study, habitat unknown"],
              ps["MAG, other study, habitat unknown"]),
             ("MAGs, other study, same habitat", sh["MAG, other study, same habitat"],
              ps["MAG, other study, same habitat"]),
             ("MAGs, same study", sh["MAG, same study"], ps["MAG, same study"])]
    rows = []
    for lab, v, p in kinds:
        r = [lab, f"{v['genes']:,}", pct(v["mobile"]), pct(v["ribosome"], 2), pct(v["resistance"])]
        r += ([f"{p['pairs']:,}", pct(p["with_mobile"] / p["pairs"]), pct(p["with_ribosome"] / p["pairs"])] if p
              else ["\\(em"] * 3)
        rows.append("@".join(r))
    tbl("tab_full_content.tbl", ["c c c s s c c s", "c c c c c c c c"], "l n n n n n n n",
        ["genome pair@genes@share of genes (%)@pairs@share of pairs (%)",
         "@@mobile@ribosomal@resistance@@mobile gene@ribosomal gene"], rows)
    F["provenance"] = {"sources": PR, "cosample": cs, "mag_pairs_provenance": MC["mag_pairs_provenance"],
                       "content": {"background": bg, "between_genera": sh,
                                   "same_habitat_between_families": MC["same_habitat_between_families"],
                                   "pairs": ps},
                       "ecology": {lab.replace("\\h'0.8m'", "").replace("\\h'1.6m'", "  "):
                                   {"genomes": e["genomes"], "exchanging_pairs": e["exchanging_pairs"],
                                    "mh": e["mh_odds_cross_genus"], "mh_ci": e["mh_ci"]["cross_genus"],
                                    "ratios": {d: e["strata"][d].get("ratio") for d in DIST}}
                                   for lab, e in sets}}


def heatmap(M, name):
    """observed / expected exchanging pairs for every pair of habitats (log colour), observed counts printed"""
    n = len(HABS)
    Z = np.full((n, n), np.nan)
    O = np.zeros((n, n))
    for i, a in enumerate(HABS):
        for j, b in enumerate(HABS):
            if j > i:
                continue
            v = M.get(f"{a}|{b}") or M.get(f"{b}|{a}")
            if v and v["expected"] >= 2:
                Z[i, j] = math.log2(max(v["observed"], 0.5) / v["expected"])
                O[i, j] = v["observed"]
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#e6e5e1", "#f4a582", "#c0451a"])
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(4.9, 4.0))
    im = ax.imshow(Z, cmap=cmap, vmin=-4.5, vmax=4.5)
    for i in range(n):
        for j in range(n):
            if j <= i and not np.isnan(Z[i, j]):
                ax.text(j, i, f"{int(O[i, j])}", ha="center", va="center", fontsize=5.5,
                        color="white" if abs(Z[i, j]) > 3 else "#1a1a19")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([HLAB[h] for h in HABS], rotation=55, ha="right", fontsize=7)
    ax.set_yticklabels([HLAB[h] for h in HABS], fontsize=7)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.75, ticks=[-4, -2, 0, 2, 4])
    cb.ax.set_yticklabels(["1/16", "1/4", "1", "4", "16"], fontsize=7)
    cb.set_label("observed / expected", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name))
    plt.close(fig)


def main():
    G = genomes()
    S = json.load(open(os.path.join(RES, "complete_complete", "summary.json")))
    S32 = json.load(open(os.path.join(RES, "complete_complete_c32", "summary.json")))
    SA = json.load(open(os.path.join(RES, "complete", "summary.json")))
    E = json.load(open(os.path.join(RES, "complete_complete", "ecology.json")))
    A = json.load(open(os.path.join(RES, "complete_complete", "annotation.json")))
    VP = json.load(open(os.path.join(RES, "complete", "validate_precision.json")))
    VB = json.load(open(os.path.join(RES, "complete", "validate_blast_truth.json")))
    VV = json.load(open(os.path.join(RES, "complete", "validate_variants.json")))
    ST = json.load(open(os.path.join(HD, "complete", "search_time.json")))
    ST32 = json.load(open(os.path.join(HD, "complete", "search_time_c32.json")))
    PL = pair_level(G)
    N = {}

    # ---- validation table: window level (BLAST truth, 60 genomes) and genome-pair level, with time for 300 genomes
    rows = []
    order = [("minimap2_single", "minimap2 (asm5, every secondary)", None),
             ("w5000_c32", "aif, 5 kb windows, 32 candidates", "w5000_c32"),
             ("w5000_c128", "aif, 5 kb, 128 candidates (atlas)", "w5000_c128"),
             ("w5000_c256", "aif, 5 kb, 256 candidates", "w5000_c256"),
             ("w2000_c64", "aif, 2 kb, 64 candidates", "w2000_c64"),
             ("w1000_c64", "aif, 1 kb, 64 candidates", "w1000_c64"),
             ("w5000_c128_m99", "aif, 5 kb, 128, AF at $>=$ 99%", "w5000_c128_m99")]
    for key, lab, vk in order:
        v = VB[key]
        t = f"{VV[vk]['aif_s']:.0f}" if vk else "390"
        pl = PL.get(key)
        plt_ = f"{pct(pl['recall'], 0)} / {pct(pl['precision'], 0)}" if pl else "\\(em"
        rows.append(f"{lab}@{t}@{pct(v['recall_cross_genus'], 0)} / {pct(v['precision_cross_genus'], 0)}@"
                    f"{pct(v['recall'], 0)} / {pct(v['precision'], 0)}@{plt_}")
    tbl("tab_valid.tbl", ["l c c s c", "l c c c c"], "l n c c c",
        ["method@time (s)@windows: recall / precision (%)@genome pairs (%)",
         "@300 genomes@cross-genus@all@cross-genus"], rows)
    N["blast_truth"] = {k: VB[k] for k in ("genomes", "windows", "truth", "truth_cross_genus")}
    N["pair_level"] = PL
    N["precision"] = VP["summary"]
    N["precision_blast"] = {k: VP["blast"][k] for k in ("hsp500", "cov500", "n")}

    # ---- extent table by distance (complete-complete, cap 128 vs cap 32)
    rows = []
    for d in DIST:
        rows.append(f"{d}@{S['pairs_by_distance'][d]:,}@{S32['pairs_by_distance'][d]:,}@{SA['pairs_by_distance'][d]:,}")
    tbl("tab_dist.tbl", ["c c c c", "c c c c"], "l n n n",
        ["lowest common taxon@complete partners@complete partners@any partner",
         "@(128 candidates)@(32 candidates)@(128 candidates)"], rows)

    # ---- ecology strata table
    rows = []
    for d in DIST[:-1]:
        v = E["strata"][d]
        ci = v.get("ratio_ci") or [float("nan")] * 2
        rows.append(f"{d}@{v['same_exchanging']:,.0f} / {v['same_pairs']:,.0f}@{100 * v['rate_same']:.3f}@"
                    f"{v['diff_exchanging']:,.0f} / {v['diff_pairs']:,.0f}@{100 * v['rate_diff']:.3f}@"
                    f"{v['ratio']:.2f} ({ci[0]:.2f}\\(en{ci[1]:.2f})")
    tbl("tab_strata.tbl", ["c c s c s c", "c c c c c c"], "l n n n n n",
        ["lowest common@same habitat@different habitat@ratio",
         "taxon@exchanging / pairs@rate (%)@exchanging / pairs@rate (%)@(95% CI)"], rows)
    N["mh"] = {"cross_genus": [E["mh_odds_cross_genus"], *E["mh_ci"]["cross_genus"]],
               "cross_family": [E["mh_odds_cross_family"], *E["mh_ci"]["cross_family"]],
               "same_genus": [E["mh_odds_same_genus"], *E["mh_ci"]["same_genus"]]}

    # ---- figure: exchange rate by distance and habitat sameness
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    xs = list(range(6))
    for key, col, lab in (("rate_same", C1, "same habitat"), ("rate_diff", C2, "different habitat")):
        ys = [E["strata"][d][key] for d in DIST[:-1]]
        ax.plot(xs, ys, "o-", color=col, ms=4, lw=1.2, mec="white", mew=0.5, label=lab)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(DLAB[:-1])
    ax.set_xlabel("lowest common taxon of the pair")
    ax.set_ylabel("pairs sharing $\\geq$1 kb at $\\geq$99%")
    ax.legend(loc="upper right", fontsize=7)
    axes_style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_strata.eps"))
    plt.close(fig)

    # ---- figure: habitat-pair enrichment heatmap (cross-order)
    M = E["habitat_pairs_cross_order"]
    heatmap(M, "fig_habitats.eps")
    N["habitat_order"] = {k: v for k, v in M.items()}

    # ---- function table and figure
    rows = []
    for t in A["classes"]:
        rows.append(f"{t['class']}@{t['shared']:,}@{pct(t['shared_share'])}@{t['background']:,}@{pct(t['background_share'])}"
                    f"@{t['odds_ratio']:.2f}")
    tbl("tab_func.tbl", ["c c s c s c", "c c c c c c"], "l n n n n n",
        ["gene class@shared blocks@background windows@odds", "@genes@%@genes@%@ratio"], rows)
    rows = []
    for t in A["top_kos"][:18]:
        nm = t["name"].split("; ", 1)[-1] if "; " in t["name"] else t["name"]
        nm = nm.replace("[EC:", "(EC ").replace("]", ")")[:72]
        rows.append(f"{t['ko']}@{nm}@{t['shared']}@{t['background']}")
    tbl("tab_kos.tbl", ["l l c c"], "l l n n", ["KEGG@function@shared@background"], rows)
    byd = A["by_distance"]
    N["func_by_distance"] = {d: {k: v for k, v in byd[d].items()} for d in byd}
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    cl = [t for t in A["classes"] if t["class"] not in ("no KEGG ortholog",)]
    cl.sort(key=lambda t: t["odds_ratio"])
    ax.barh(range(len(cl)), [math.log2(t["odds_ratio"]) for t in cl], color=[C2 if t["odds_ratio"] > 1 else C1 for t in cl],
            height=0.62)
    ax.set_yticks(range(len(cl)))
    ax.set_yticklabels([t["class"] for t in cl], fontsize=7)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xticklabels(["1/8", "1/4", "1/2", "1", "2", "4", "8"])
    ax.set_xlabel("odds ratio, shared blocks against background windows")
    axes_style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_func.eps"))
    plt.close(fig)

    # ---- AMR / metal share by habitat of the query genome (from the annotation files)
    import hgt_annotate as ha
    names, brite = ha.kegg()
    ko = {}
    for l in open(os.path.join(HD, "annot_complete_complete", "ko.tsv")):
        q, s_ = l.split("\t")[:2]
        if q not in ko and "~" in s_:
            ko[q] = s_.split("~")[1]
    by = collections.defaultdict(collections.Counter)
    for l in open(os.path.join(HD, "annot_complete_complete", "proteins.faa")):
        if l[0] != ">":
            continue
        gid = l[1:].split()[0]
        p = gid.rsplit("_", 1)[0].split("|")
        by[(p[0], G[p[1]]["habitat"])][ha.classify(ko.get(gid), names, brite)] += 1
    rows = []
    hab_func = {}
    for h in HABS:
        R_, B_ = by[("R", h)], by[("B", h)]
        nr, nb = sum(R_.values()), sum(B_.values())
        if nr < 300:
            continue
        hab_func[h] = {"genes": nr, "amr": R_["antimicrobial resistance"] / nr, "amr_bg": B_["antimicrobial resistance"] / nb,
                       "metal": R_["metal resistance"] / nr, "transposase": R_["transposase"] / nr,
                       "conj": R_["conjugation / T4SS"] / nr}
        x = hab_func[h]
        rows.append(f"{HLAB[h]}@{nr:,}@{pct(x['amr'])}@{pct(x['amr_bg'])}@{pct(x['metal'])}@{pct(x['transposase'])}"
                    f"@{pct(x['conj'])}")
    tbl("tab_habfunc.tbl", ["c c c c c c c", "c c c c c c c"], "l n n n n n n",
        ["habitat of the@genes on@resistance@resistance,@metal@transposase@conjugation",
         "query genome@shared blocks@(%)@background (%)@resistance (%)@(%)@(%)"], rows)
    N["hab_func"] = hab_func

    # ---- hubs and cross-phylum pairs
    rows = [f"{g}@{n}" for g, n in S["top_genera_by_partner_genera"][:12]]
    tbl("tab_hubs.tbl", ["l c"], "l n", ["genus@partner genera"], [r.replace("_", "\\_") if False else r for r in rows])
    rows = [f"{a}@{n:,}" for a, n in S["cross_phylum_pairs_top"][:10]]
    tbl("tab_phyla.tbl", ["l c"], "l n", ["phyla@genome pairs"], rows)

    # ---- timing
    N["time"] = {"c128_s": sum(v["aif_s"] for v in ST.values()), "c32_s": sum(v["aif_s"] for v in ST32.values()),
                 "bytes": sum(v["bytes"] for v in ST.values()), "peak_gb": max(v["kb"] for v in ST.values()) / 1e6}
    N["summary"] = {k: S[k] for k in S if k not in ("ecology",)}
    N["summary_c32"] = {k: S32[k] for k in ("cross_genus_pairs", "genomes_with_cross_genus", "genomes_with_cross_family",
                                            "genomes_with_cross_phylum", "genus_pairs", "genera_involved")}
    N["summary_any"] = {k: SA[k] for k in ("cross_genus_pairs", "genomes_with_cross_genus", "genomes_with_cross_family",
                                           "genomes_with_cross_phylum", "genus_pairs", "genera_involved",
                                           "cross_genus_partner_category")}
    N["annotation"] = {k: A[k] for k in ("regions", "region_bp", "genomes", "genes_shared", "genes_background")}
    N["ecology"] = {k: E[k] for k in ("genomes", "habitats", "exchanging_pairs")}
    if os.path.exists(os.path.join(RES, FULL, "annotation.json")):
        full_atlas(G, N, A)
    json.dump(N, open(os.path.join(GEN, "numbers.json"), "w"), indent=1, default=float)
    print("tables and figures written")


if __name__ == "__main__":
    main()
