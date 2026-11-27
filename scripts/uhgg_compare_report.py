#!/usr/bin/env python3
"""Tables, figure and numbers for the tool comparison and the v1.7 AF and interval changes in the UHGG report.

Reads results_uhgg/compare (eval_S3.json, eval_S1.json, cost.json, v16/, af_tuning_T.json, af_diag_T_base.json),
results_real (af_local_check.json, calibration.json, calibration_v16.json) and the Kraken2 build record; writes
report_uhgg/gen/tab_cmp_*.tbl, gen/numbers_cmp.txt and figs/fig_cost.eps.
Usage: /usr/bin/python3 scripts/uhgg_compare_report.py
"""
import json
import os

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CMP = os.path.join(ROOT, "results_uhgg", "compare")
RR = os.path.join(ROOT, "results_real")
GEN = os.path.join(ROOT, "report_uhgg", "gen")
FIG = os.path.join(ROOT, "report_uhgg", "figs")
INK, INK2, GRID, HI, MID = "#1f1e1c", "#52514e", "#d9d8d4", "#2a78d6", "#8a8984"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2})
NAME = {"aif": "aif v1.7", "minimap2": "minimap2", "sourmash": "sourmash", "skani": "skani", "kraken2": "Kraken2",
        "cat": "CAT", "lexicmap": "LexicMap"}
DASH = "\\(em"


def J(*p):
    return json.load(open(os.path.join(*p)))


def tbl(path, fmt, header_rows, rows, opts="center box tab(@)"):
    lines = [".TS", opts + ";"] + fmt[:-1] + [fmt[-1] + "."]
    for i, hr in enumerate(header_rows):
        cols = fmt[i].split() if i < len(fmt) - 1 else []
        if len(cols) == len(hr):  # tbl wants no entry for a spanned ("s") column
            hr = [h for h, c in zip(hr, cols) if c != "s"]
        lines.append("@".join(hr))
    lines.append("_")
    for r in rows:
        lines.append("_" if r == "_" else "@".join(r))
    lines.append(".TE")
    open(os.path.join(GEN, path), "w").write("\n".join(lines) + "\n")


def f(x, d=2):
    return DASH if x is None or x != x else f"{x:.{d}f}"


def pct(x, d=1):
    return DASH if x is None or x != x else f"{100 * x:.{d}f}"


def dur(s):
    if s is None:
        return DASH
    if s < 600:
        return f"{s:.1f} s" if s < 10 else f"{s:.0f} s"
    if s < 2 * 3600:
        return f"{s / 60:.0f} min"
    if s < 3 * 86400:
        return f"{s / 3600:.1f} h"
    if s < 2 * 365 * 86400:
        return f"{s / 86400:.0f} d"
    return f"{s / (365 * 86400):.0f} y"


def accuracy_rows(E, tools, extra_v16=None):
    rows = []
    for t in tools:
        o = E["q"][t]
        t1, s, d = o.get("top1", {}), o.get("search", {}), o["decision"]
        rows.append([NAME[t], f(t1.get("ani_mae_pp")), f(t1.get("af_mae_pp")), pct(t1.get("ani_within1pp"), 0),
                     pct(s.get("reported")), pct(s.get("top_hit")), f(d["precision"], 3), f(d["recall"], 3),
                     f(d["f1"], 3), f(d["precision_w"], 3), f(d["recall_w"], 3)])
        if t == "aif" and extra_v16 is not None:
            o = extra_v16["q"]["aif"]
            t1, s, d = o["top1"], o["search"], o["decision"]
            rows.append(["aif v1.6", f(t1.get("ani_mae_pp")), f(t1.get("af_mae_pp")), pct(t1.get("ani_within1pp"), 0),
                         pct(s.get("reported")), pct(s.get("top_hit")), f(d["precision"], 3), f(d["recall"], 3),
                         f(d["f1"], 3), f(d["precision_w"], 3), f(d["recall_w"], 3)])
    return rows


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    E3, E1 = J(CMP, "eval_S3.json"), J(CMP, "eval_S1.json")
    V3, V1 = J(CMP, "v16", "eval_S3.json"), J(CMP, "v16", "eval_S1.json")
    cost = J(CMP, "cost.json")
    N = {}
    hdr = [["Tool", "ANI", "AF", "ANI", "has a", "top hit", "foreign call", "", "", "catalogue-weighted", ""],
           ["", "MAE (pp)", "MAE (pp)", "\\(<= 1 pp (%)", "hit (%)", "best (%)", "P", "R", "F1", "P", "R"]]
    fmt = ["l c c c c c c s s c s", "l c c c c c c c c c c", "l n n n n n n n n n n"]
    tools3 = ["aif", "minimap2", "sourmash", "skani", "kraken2", "cat"]
    tbl("tab_cmp_acc.tbl", fmt, hdr, accuracy_rows(E3, tools3, V3))
    tools1 = ["aif", "lexicmap", "minimap2", "sourmash", "skani", "kraken2", "cat"]
    tbl("tab_cmp_s1.tbl", fmt, hdr, accuracy_rows(E1, tools1, V1))

    # by convention: aif v1.6 / v1.7 and the AF tools, top hit
    rows = []
    for conv, lab in (("q", "delta-filter \\-q"), ("one", "delta-filter \\-1"), ("q90", "\\(>= 90% identity")):
        r = [lab]
        for E, t in ((V3, "aif"), (E3, "aif"), (E3, "skani"), (E3, "minimap2")):
            r.append(f(E[conv][t]["top1"].get("af_mae_pp")))
        p = E3[conv]["paired"]["skani"]
        r += [f(p["af_mae_pp_aif"]), f(p["af_mae_pp_tool"])]
        r += [f(E3[conv]["aif"]["top1"]["ani_mae_pp"], 3), f(E3[conv]["skani"]["top1"]["ani_mae_pp"], 3),
              f(E3[conv]["minimap2"]["top1"]["ani_mae_pp"], 3)]
        rows.append(r)
    tbl("tab_cmp_conv.tbl", ["l c s s s c s c s s", "l c c c c c c c c c", "l n n n n n n n n n"],
        [["MUMmer convention", "AF MAE, top hit (pp)", "", "", "", "AF MAE, shared pairs", "", "ANI MAE, top hit (pp)", "", ""],
         ["", "aif v1.6", "aif v1.7", "skani", "minimap2", "aif v1.7", "skani", "aif", "skani", "minimap2"]], rows)

    # cost
    k2b = J(ROOT, "data", "gtdb", "kraken2_db", "build.json")
    B = dict(cost["builds"])
    B["kraken2"] = {"build_s": k2b["build_s"], "bytes": k2b["bytes"], "threads": k2b["threads"]}
    T = cost["tools"]
    rows = []
    for t in ["aif", "kraken2", "skani", "sourmash", "minimap2", "cat", "lexicmap"]:
        b, r = B[t], T[t]
        bs = b.get("build_s_active") or b["build_s"]
        star = "\\(dg" if b.get("build_s_active") else ""
        pts = r["points"]
        s3 = pts.get("S3", {}).get("s")
        pool = pts.get("pool", {}).get("s")
        proj = r["projected_full_s"]
        basis = {"pool rate, S1-S3 fixed cost": "pool", "S1-S3 fit": "S1\\(enS3", "S1 rate, no fixed cost": "S1"}[r["basis"]]
        rows.append([NAME[t].replace(" v1.7", ""), f"{b['bytes'] / 1e9:.1f}" if b["bytes"] < 1e12 else f"{b['bytes'] / 1e12:.2f} TB",
                     dur(bs) + star, str(b.get("threads", "")), dur(s3) if s3 else f"{dur(pts['S1']['s'])} (S1)",
                     dur(pool), dur(proj), basis, f"{r['peak_gb']:.0f}"])
        N[f"proj_h_{t}"] = proj / 3600
        N[f"peak_gb_{t}"] = r["peak_gb"]
    tbl("tab_cmp_cost.tbl", ["l c c c c c c c c", "l c c c c c c c c", "l n n n n n n c n"],
        [["Tool", "index (GB)", "build", "build", "S3 search", "pool search", "full catalogue", "basis", "search peak"],
         ["", "", "", "threads", "(114 Mbp)", "(3.64 Gbp)", "(706 Gbp)", "", "memory (GB)"]], rows)
    N["aif_full_measured_h"] = cost["aif_full_measured_s"] / 3600
    N["full_gbp"] = cost["full_bp"] / 1e9

    # AF fix: held-out and tuning
    ch = J(RR, "af_local_check.json")
    tune = J(CMP, "af_tuning_T.json")["scores"]
    rows = [["UHGG tuning set T (top hit, \\-q)", "tuning", f(tune["q|base|top1"]["af_mae"]), f(tune["q|v17|top1"]["af_mae"]),
             f(tune["q|base|top1"]["af_bias"]), f(tune["q|v17|top1"]["af_bias"])]]
    for E0, E, lab in ((V3, E3, "UHGG S3 (top hit, \\-q)"), (V1, E1, "UHGG S1 (top hit, \\-q)")):
        a, b = E0["q"]["aif"]["top1"], E["q"]["aif"]["top1"]
        rows.append([lab, "held out", f(a["af_mae_pp"]), f(b["af_mae_pp"]), f(a["af_bias_pp"]), f(b["af_bias_pp"])])
    for key, lab in (("real", "real genomes, 5 sets"), ("real_fresh", "real genomes, fresh drafts"),
                     ("bench", "simulated, main"), ("bench2", "simulated, fresh")):
        a, b = ch["v16"][key]["macro"], ch["default"][key]["macro"]
        rows.append([lab, "held out", f(a["af_mae"]), f(b["af_mae"]), f(a["af_bias"]), f(b["af_bias"])])
        N[f"af_{key}_v16"], N[f"af_{key}_v17"] = a["af_mae"], b["af_mae"]
    hb = [(ch[v]["bench"]["groups"]["hetblock_10k_ani95"]["af_bias"], ch[v]["bench"]["groups"]["hetgamma_10k_ani95"]["af_bias"])
          for v in ("v16", "default")]
    N["hetblock_v16"], N["hetgamma_v16"], N["hetblock_v17"], N["hetgamma_v17"] = hb[0][0], hb[0][1], hb[1][0], hb[1][1]
    tbl("tab_cmp_affix.tbl", ["l c c s c s", "l c c c c c", "l c n n n n"],
        [["Data set", "", "AF MAE (pp)", "", "AF bias (pp)", ""], ["", "role", "v1.6", "v1.7", "v1.6", "v1.7"]], rows)

    # the AF models that framed the choice, on the three tuning sets
    low = J(CMP, "af_lowani_T.json")
    tch = J(RR, "af_local_check_tune.json")
    order = [("base", "v16", "v1.6"), ("ami0.89", "flat", "flat: identity \\(>= 0.89"),
             ("al0.05", "local", "local threshold only"), ("v17", "default", "v1.7 (local + evidence)")]
    rows = []
    for name, tk, lab in order:
        e = tune[f"q|{name}|top1"]
        lo = [low[name][b] for b in ("0.80-0.85", "0.85-0.88", "0.88-0.90")]
        lo_mae = sum(x["af_mae"] * x["n"] for x in lo) / sum(x["n"] for x in lo)
        st, rt = tch[tk]["tune"], tch[tk]["real_tune"]
        rows.append([lab, f(e["af_mae"]), f(e["af_bias"]), f(lo_mae), f(rt["macro"]["af_mae"]), f(st["macro"]["af_mae"]),
                     f(st["groups"]["hetblock_10k_ani95"]["af_bias"], 1)])
    tbl("tab_cmp_afvar.tbl", ["l c s s c c s", "l c c c c c c", "l n n n n n n"],
        [["AF model", "UHGG tuning set T", "", "", "real tuning", "simulated tuning", ""],
         ["", "AF MAE", "AF bias", "AF MAE at ANI 0.80\\(en0.90", "AF MAE", "AF MAE", "hetblock bias"]], rows)

    # intervals: coverage of nominal 95%
    c16, c17 = J(RR, "calibration_v16.json")["calibrated"]["0.95"], J(RR, "calibration.json")["calibrated"]["0.95"]
    rows = []
    for E0, E, lab in ((V3, E3, "UHGG S3 top hits (db search)"), (V1, E1, "UHGG S1 top hits (db search)")):
        a, b = E0["q"]["aif"]["top1"], E["q"]["aif"]["top1"]
        rows.append([lab, pct(a["ci_cover_ani"]), pct(b["ci_cover_ani"]), pct(a["ci_cover_af"]), pct(b["ci_cover_af"]), DASH, DASH])
    for key, lab in (("real_eval", "real genomes, 5 sets (query)"), ("real_fresh", "real genomes, fresh drafts (query)"),
                     ("sim_main", "simulated, main (query)")):
        a = c16.get(key)
        b = c17[key]
        rows.append([lab, pct(a["ani"]) if a else DASH, pct(b["ani"]), pct(a["af"]) if a else DASH, pct(b["af"]),
                     f(a["halfwidth_af_med"]) if a else DASH, f(b["halfwidth_af_med"])])
    tbl("tab_cmp_ci.tbl", ["l c s c s c s", "l c c c c c c", "l n n n n n n"],
        [["Held-out set", "ANI cover (%)", "", "AF cover (%)", "", "AF half-width (pp)", ""],
         ["", "v1.6", "v1.7", "v1.6", "v1.7", "v1.6", "v1.7"]], rows)

    # real-genome report (Section 17): per-set AF of v1.6 and v1.7, and pairwise interval coverage
    RGEN = os.path.join(ROOT, "report_real", "gen")
    SETN = {"bcereus": "B. cereus", "ecoli": "E. coli", "mag2328": "MAGs (bin 2328)", "paer": "P. aeruginosa",
            "saur": "S. aureus"}
    rows = []
    for key, lab in (("real", "main"), ("real_fresh", "fresh")):
        for g in sorted(ch["default"][key]["groups"]):
            a, b = ch["v16"][key]["groups"][g], ch["default"][key]["groups"][g]
            rows.append([f"{lab}: \\f[I]{SETN[g]}\\f[R]" if g != "mag2328" else f"{lab}: {SETN[g]}", str(b["n"]),
                         f(a["af_mae"]), f(b["af_mae"]), f(a["af_bias"]), f(b["af_bias"]), f(b["ani_mae"], 3)])
        a, b = ch["v16"][key]["macro"], ch["default"][key]["macro"]
        rows.append([f"{lab}: macro", "", f(a["af_mae"]), f(b["af_mae"]), f(a["af_bias"]), f(b["af_bias"]), f(b["ani_mae"], 3)])
        rows.append("_")
    for key, lab in (("bench", "simulated main"), ("bench2", "simulated fresh")):
        a, b = ch["v16"][key], ch["default"][key]
        rows.append([f"{lab}: macro", "", f(a["macro"]["af_mae"]), f(b["macro"]["af_mae"]), f(a["macro"]["af_bias"]),
                     f(b["macro"]["af_bias"]), f(b["macro"]["ani_mae"], 3)])
        for c in ("hetblock_10k_ani95", "hetgamma_10k_ani95"):
            rows.append([f"\\0\\0{c.split('_')[0]}", str(b["groups"][c]["n"]), f(a["groups"][c]["af_mae"]),
                         f(b["groups"][c]["af_mae"]), f(a["groups"][c]["af_bias"]), f(b["groups"][c]["af_bias"]),
                         f(b["groups"][c]["ani_mae"], 3)])
    gen0 = GEN
    globals()["GEN"] = RGEN
    tbl("tab_v17_af.tbl", ["l c c s c s c", "l c c c c c c", "l n n n n n n"],
        [["Set", "contigs", "AF MAE (pp)", "", "AF bias (pp)", "", "ANI MAE (pp)"],
         ["", "", "v1.6", "v1.7", "v1.6", "v1.7", "v1.6 = v1.7"]], rows)
    rows = []
    c16a, c17a = J(RR, "calibration_v16.json")["calibrated"], J(RR, "calibration.json")["calibrated"]
    for key, lab in (("real_eval", "real, main"), ("real_fresh", "real, fresh"), ("sim_main", "simulated, main"),
                     ("sim_fresh", "simulated, fresh")):
        r = [lab]
        for lv in ("0.9", "0.95", "0.99"):
            a, b = c16a[lv].get(key), c17a[lv][key]
            r += [f"{pct(a['ani'])}/{pct(a['af'])}" if a else DASH, f"{pct(b['ani'])}/{pct(b['af'])}"]
        rows.append(r)
    tbl("tab_v17_ci.tbl", ["l c s c s c s", "l c c c c c c", "l c c c c c c"],
        [["Held-out set", "90%", "", "95%", "", "99%", ""], ["", "v1.6", "v1.7", "v1.6", "v1.7", "v1.6", "v1.7"]], rows)
    globals()["GEN"] = gen0

    # base-level diagnosis
    dg = J(CMP, "af_diag_T_base.json")
    for k, v in dg["frac"].items():
        N[f"diag_{k}"] = v
    N["diag_pairs"] = dg["pairs"]
    for k, v in dg["by_aif_window_identity"].items():
        N[f"diag_id_{k}_aifonly_share"] = v["aif_only"] / max(1e-12, v["aif_only"] + v["shared"])

    # figure: projected full-catalogue time against foreign-call F1 (S3) and top-hit AF MAE
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.5))
    pts = []
    for t in ["aif", "kraken2", "skani", "sourmash", "minimap2", "cat", "lexicmap"]:
        E = E1 if t == "lexicmap" else E3
        pts.append((t, T[t]["projected_full_s"] / 3600, E["q"][t]["decision"]["f1"], E["q"][t].get("top1", {}).get("af_mae_pp")))
    # label offsets (points) and alignment per panel, set by eye against collisions
    OFF = [{"aif": (0, 7, "center"), "kraken2": (0, 7, "center"), "skani": (7, -3, "left"), "sourmash": (6, -10, "left"),
            "minimap2": (6, -10, "left"), "cat": (6, 2, "left"), "lexicmap": (-6, -10, "right")},
           {"aif": (6, -10, "left"), "skani": (6, 2, "left"), "minimap2": (6, 2, "left"), "lexicmap": (-8, 6, "right")}]
    for i, (ax, j, ylab) in enumerate(((axs[0], 2, "foreign-contig F1 (S3)"), (axs[1], 3, "top-hit AF MAE (pp, S3)"))):
        for p in pts:
            if p[j] is None:
                continue
            c = HI if p[0] == "aif" else MID
            ax.scatter([p[1]], [p[j]], s=26, color=c, edgecolor="white", linewidth=0.8, zorder=3)
            lab = NAME[p[0]] + (" (S1)" if p[0] == "lexicmap" else "")
            dx, dy, ha = OFF[i][p[0]]
            ax.annotate(lab, (p[1], p[j]), textcoords="offset points", xytext=(dx, dy), ha=ha, fontsize=7, color=INK)
        ax.set_xscale("log")
        ax.set_xlabel("projected time, full catalogue (h, 48 threads)")
        ax.set_ylabel(ylab)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xlim(0.5, 3e6)
    axs[0].set_ylim(0.35, 0.95)
    axs[1].set_ylim(0, 6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_cost.eps"))
    with open(os.path.join(GEN, "numbers_cmp.txt"), "w") as fh:
        for k, v in sorted(N.items()):
            fh.write(f"{k}\t{v}\n")
    print("wrote", ", ".join(sorted(x for x in os.listdir(GEN) if x.startswith("tab_cmp"))), "fig_cost.eps, numbers_cmp.txt")


if __name__ == "__main__":
    main()
