#!/usr/bin/env python3
"""Tables, figures and numbers for the UHGG contamination report (report_uhgg/).

Reads results_uhgg/summary.json and kraken2_check.json (bench_real/uhgg_summary.py) and writes report_uhgg/gen/*.tbl,
report_uhgg/gen/numbers.txt and report_uhgg/figs/*.eps.
Usage: /usr/bin/python3 scripts/uhgg_report.py            (UHGG_PILOT=1: the pilot's results_uhgg/pilot)
"""
import json
import os
from collections import OrderedDict

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results_uhgg", "pilot") if os.environ.get("UHGG_PILOT") else os.path.join(ROOT, "results_uhgg")
GEN = os.path.join(ROOT, "report_uhgg", "gen")
FIG = os.path.join(ROOT, "report_uhgg", "figs")
LB = ["0.5\\(en2", "2\\(en5", "5\\(en20", "\\(>= 20"]
LBP = ["0.5–2", "2–5", "5–20", "≥20"]
C_MAG, C_ISO, C_BG = "#2a78d6", "#eb6834", "#9a9893"
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


def pct(x, d=1):
    return "\\(em" if x is None else f"{100 * x:.{d}f}"


def ci(w, d=1):
    return f"{100 * w[0]:.{d}f} ({100 * w[1]:.{d}f}\\(en{100 * w[2]:.{d}f})"


def axes_style(ax):
    ax.grid(True, color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def extra_numbers(label_refs):
    """numbers quoted in the text: MAGs >= 5% contaminated against CheckM2, representatives >= 1%, and the flagged
    contigs whose best hit is a label-inconsistent reference (flagged_contigs.tsv.gz and the metadata, if present)"""
    import csv
    import gzip
    import statistics
    out = {}
    meta = os.path.join(ROOT, "data", "uhgg3", "meta", "genomes-all_metadata.tsv")
    if os.path.exists(meta):
        cm = {r["Genome"]: (r["Genome_type"], float(r["Contamination"])) for r in csv.DictReader(open(meta), delimiter="\t")}
        g5 = [(float(r["contam"]), cm[r["genome"]][1]) for r in csv.DictReader(gzip.open(os.path.join(RES, "genomes.tsv.gz"), "rt"),
                                                                                 delimiter="\t")
              if cm[r["genome"]][0] == "MAG" and float(r["contam"]) >= 0.05]
        out["mags_ge5"] = {"n": len(g5), "checkm_median": statistics.median(c for _, c in g5),
                           "aif_median": statistics.median(a for a, _ in g5), "checkm_lt2": sum(c < 2 for _, c in g5)}
    reps = list(csv.DictReader(open(os.path.join(RES, "species_reps_contaminated.tsv")), delimiter="\t"))
    out["reps_ge1"] = sum(float(r["aif_contam"]) >= 0.01 for r in reps)
    fc = os.path.join(RES, "flagged_contigs.tsv.gz")
    if os.path.exists(fc):
        bad = {r[0] for r in label_refs}
        n = nb = bp = bpb = abp = abpb = 0
        for r in csv.DictReader(gzip.open(fc, "rt"), delimiter="\t"):
            L, hit = int(r["length"]), r["ref"] in bad
            n += 1
            nb += hit
            bp += L
            bpb += L * hit
            if r["anchored"] == "1":
                abp += L
                abpb += L * hit
        out["label_resid"] = {"contigs": n, "share": nb / n, "gbp": bpb / 1e9, "anchored_share": abpb / max(1, abp)}
    return out


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    S = json.load(open(os.path.join(RES, "summary.json")))
    kp = os.path.join(RES, "kraken2_check.json")
    KC = json.load(open(kp)) if os.path.exists(kp) else {}
    nums = OrderedDict()
    g = S["groups"]
    # (1) data and run
    sr = S.get("search", {})
    pr = S.get("prepare", {})
    rows = []
    for k, lab in (("isolates", "isolates"), ("mags", "MAGs"), ("hq_mags", "\\0\\0high-quality MAGs"), ("reps", "species representatives")):
        d = g[k]
        hs = d["host_src"]
        rows.append([lab, f"{d['n']:,}", f(d["bp"] / 1e9, 1), pct(d["hit_share"]), pct(d["species_hit_share"]),
                     pct(hs.get("self", 0) / d["n"]), pct(hs.get("species", 0) / d["n"]), pct((hs.get("weak", 0) + hs.get("none", 0)) / d["n"])])
    tbl("tab_data.tbl", ["l c c c s c s s", "^ ^ ^ c c c c c"], "l n n n n n n n",
        [["genomes", "number", "Gbp", "bases with a GTDB hit (%)", "", "host genus from (% of genomes)", "", ""],
         ["", "", "", "any", "ANI \\(>= 95%", "self", "species", "weak"]], rows)
    nums["data"] = {k: (g[k]["n"], g[k]["bp"]) for k in ("all", "isolates", "mags", "hq_mags", "reps")}
    nums["search"] = sr
    nums["prepare"] = pr
    # (2) flagged contigs by class
    fc = S["flagged_contigs"]
    nums["flagged"] = (fc, S["flagged_anchored"], S["shared"], S["contamination_sources"])
    # (3) contig profiles
    P = S["contig_profiles"]
    rows = []
    for t, tl in (("MAG", "MAGs"), ("Isolate", "isolates")):
        for b in range(4):
            fl, mo, bg = P[f"{t}_flag_chrom_{b}"], P[f"{t}_flag_mobile_{b}"], P[f"{t}_bg_{b}"]
            rows.append([tl if b == 0 else "", LB[b], f"{fl['n']:,}", f(fl["core_per_kb"], 4), f(bg["core_per_kb"], 4),
                         f"{fl['core_per_kb'] / bg['core_per_kb']:.2f}" if bg["core_per_kb"] else "\\(em", pct(fl["hypo_frac"], 0),
                         pct(bg["hypo_frac"], 0), f"{mo['n']:,}"])
        rows.append("_")
    tbl("tab_profile.tbl", ["l l c c s s c s c", "^ ^ c c c c c c c"], "l l n n n n n n n",
        [["genomes", "contig (kb)", "flagged", "core genes per kb", "", "", "hypothetical CDS (%)", "", "flagged"],
         ["", "", "non-mobile", "flagged", "unflagged", "ratio", "flagged", "unflagged", "mobile"]], rows[:-1])
    # (4) core mixture
    cm = S["core_mixture"]
    rows = []
    for b, d in enumerate(cm["bins"]):
        rows.append([LB[b], f(d["mag_flag_bp"] / 1e6, 1), f(d["core_per_kb_iso_flag"], 4), f(d["core_per_kb_mag_flag"], 4),
                     f(d["core_per_kb_mag_bg"], 4), pct(d["share"], 0) + f" ({pct(d['ci'][0], 0)}\\(en{pct(d['ci'][1], 0)})"])
    rows.append(["all", f(sum(d["mag_flag_bp"] for d in cm["bins"]) / 1e6, 1), "", "", "",
                 pct(cm["overall"], 0) + f" ({pct(cm['overall_ci'][0], 0)}\\(en{pct(cm['overall_ci'][1], 0)})"])
    tbl("tab_mixture.tbl", ["l c c s s c", "^ c c c c c"], "l n n n n n",
        [["contig (kb)", "flagged MAG", "core genes per kb", "", "", "contamination-like"],
         ["", "Mbp", "isolate flagged", "MAG flagged", "MAG unflagged", "share (%), 95% CI"]], rows)
    nums["mixture"] = cm
    # (5) per-genome rates
    rows = []
    for k, lab in (("isolates", "isolates"), ("mags", "MAGs"), ("hq_mags", "\\0\\0high-quality"), ("mq_mags", "\\0\\0medium-quality"),
                   ("reps_iso", "reps, isolate"), ("reps_mag", "reps, MAG"), ("heldout_isolates", "held-out isolates"),
                   ("heldout_mags", "held-out MAGs"), ("pilot_isolates", "pilot isolates"), ("pilot_mags", "pilot MAGs")):
        d = g[k]
        if not d["n"]:
            continue
        rows.append([lab, f"{d['n']:,}", pct(d["shared_any"][0]), pct(d["foreign_ge1"][0]), ci(d["contam_any"]),
                     ci(d["contam_ge1"], 2), pct(d["contam_ge5"][0], 2), pct(d["strict_any"][0])])
        nums[f"rates_{k}"] = {x: d[x] for x in ("n", "shared_any", "mobile_any", "foreign_ge1", "foreign_ge5", "contam_any", "contam_ge1",
                                                "contam_ge5", "strict_any", "contam_mean", "anchored_share", "foreign_share", "mobile_share")}
    tbl("tab_rates.tbl", ["l c c s c s s c", "^ ^ c c c c c c"], "l n n n n n n n",
        [["genomes", "n", "shared (%)", "", "core-anchored contamination (%, 95% CI)", "", "", "\\(>= 2 core"],
         ["", "", "any", "non-mob. \\(>= 1%", "any", "\\(>= 1%", "\\(>= 5%", "genes (%)"]], rows)
    # (6) species-matched
    sm = S["species_matched"]
    rows = []
    for k, lab in (("shared_any", "any shared sequence"), ("foreign_ge1", "non-mobile shared \\(>= 1%"), ("contam_any", "contamination, any"),
                   ("contam_ge1", "contamination \\(>= 1%"), ("contam_ge5", "contamination \\(>= 5%")):
        d = sm[k]
        rows.append([lab, pct(d["iso_rate"]), pct(d["mag_rate"]), f"{100 * d['mh_rd']:+.1f}", f"{d['mag_higher']} / {d['iso_higher']} / {d['tie']}",
                     f"{d['sign_p']:.1g}"])
    tbl("tab_matched.tbl", ["l c c c c c", "^ c c c c c"], "l n n n n n",
        [["per-genome call", "isolates", "MAGs", "difference", "species: MAG higher /", "sign test"],
         ["", "(%)", "(%)", "(pp)", "isolate higher / tie", "p"]], rows)
    nums["matched"] = sm
    # (7) against CheckM2 contamination and completeness
    rows = []
    for d in S["by_checkm_contam"]:
        rows.append([f"{d['lo']:g}\\(en{d['hi']:g}", f"{d['n']:,}", ci(d["contam_any"]), ci(d["contam_ge1"], 2), pct(d["foreign_ge5"][0])])
    tbl("tab_checkm.tbl", ["l c c s c", "l ^ c c c"], "l n n n n",
        [["CheckM2 contamination", "MAGs", "core-anchored contamination (%)", "", "non-mobile shared"],
         ["(MGnify, %)", "", "any", "\\(>= 1%", "\\(>= 5% (%)"]], rows)
    nums["checkm"] = (S["by_checkm_contam"], S["spearman_checkm"], S["by_completeness"], S["by_contigs"], S["hq_contaminated"])
    # (8) validation
    ss, gc, rc = S["same_sample"], S["gc"], S["recurrence_pilot"]
    rows = []
    for k, lab in (("anchored", "core-anchored"), ("accessory", "accessory"), ("mobile", "mobile"), ("noncoding", "non-coding")):
        for t, tl in (("MAG", "MAGs"), ("Isolate", "isolates")):
            d = rc.get(f"{k}_{t}")
            if d:
                rows.append([lab, tl, f"{d['n']:,}", pct(d["recurrent"]), pct(d["same_study"]), pct(d["sporadic"])])
    tbl("tab_recur.tbl", ["l l c c s c", "^ ^ ^ c c c"], "l l n n n n",
        [["flagged contigs", "genomes", "tested", "in a conspecific genome (%)", "", "absent"],
         ["", "", "", "other study", "same study only", "(%)"]], rows)
    nums["recurrence"] = rc
    nums["same_sample"] = ss
    nums["gc"] = gc
    rows = []
    for k, lab in (("contaminated", "core-anchored contamination"), ("shared_only", "shared sequence only")):
        d = ss[k]
        rows.append([lab, f"{d['mags']:,}", pct(d["observed"]), pct(d["expected"]), f(d["ratio"], 2), f(d["z"], 1)])
    tbl("tab_sample.tbl", ["l c c s c c", "^ ^ c c ^ ^"], "l n n n n n",
        [["MAGs with", "n", "source genus also binned from the sample (%)", "", "ratio", "z"], ["", "", "observed", "expected", "", ""]], rows)
    if KC:
        rows = []
        for k, lab in (("flag_anchored_MAG", "core-anchored, MAGs"), ("flag_accessory_MAG", "accessory, MAGs"), ("flag_mobile_MAG", "mobile, MAGs"),
                       ("flag_noncoding_MAG", "non-coding, MAGs"), ("flag_anchored_Isolate", "core-anchored, isolates"),
                       ("flag_accessory_Isolate", "accessory, isolates"), ("flag_mobile_Isolate", "mobile, isolates"),
                       ("flag_noncoding_Isolate", "non-coding, isolates"), ("_", ""), ("control_MAG", "control, MAGs"),
                       ("control_Isolate", "control, isolates"), ("control_Isolate_noEnt", "\\0\\0host not \\f[I]Enterobacteriaceae\\f[R]")):
            if k == "_":
                rows.append("_")
                continue
            d = KC.get(k)
            if not d:
                continue
            sf = "\\(em" if k.startswith("control") else pct(d.get("source_family", 0))
            rows.append([lab, f"{d['n']:,}", pct(d.get("source", 0)), sf, pct(d.get("host", 0)),
                         pct(d.get("host_family", 0)), pct(d.get("other", 0)), pct(d.get("above", 0) + d.get("unclassified", 0))])
        tbl("tab_kraken.tbl", ["l c c s s s s s", "^ ^ c s c s c c", "^ ^ c c c c c c"], "l n n n n n n n",
            [["contigs", "n", "Kraken2 assignment (%)", "", "", "", "", ""],
             ["", "", "source", "", "host", "", "other", "higher"],
             ["", "", "genus", "family", "genus", "family", "genus", "or none"]], rows)
        nums["kraken"] = KC
    # (9) anchoring sensitivity
    rows = [[pct(d["level"], 1), pct(d["p_core1"], 0), pct(d["p_core2"], 0)] for d in S["anchor_sensitivity"]]
    tbl("tab_sens.tbl", ["c c c"], "n n n", [["contamination (%)", "anchored, \\(>= 1 core gene (%)", "strict, \\(>= 2 (%)"]], rows)
    nums["sensitivity"] = S["anchor_sensitivity"]
    # (10) sources and pairs
    def read_tsv(fn):
        L = [l.rstrip("\n").split("\t") for l in open(os.path.join(RES, fn))]
        return L[0], L[1:]
    _, src = read_tsv("contamination_sources.tsv")
    rows = [[r[0].replace("_", "\\_"), r[1].replace("_", "\\_"), r[2].replace("_", "\\_"), f(float(r[3]) / 1e6, 2), pct(float(r[4])), r[5]] for r in src[:12]]
    tbl("tab_sources.tbl", ["l l l c c c"], "l l l n n n", [["source genus", "family", "phylum", "Mbp", "share (%)", "MAGs"]], rows)
    _, prs = read_tsv("shared_pairs.tsv")
    rows = [[r[0].replace("_", "\\_"), r[1].replace("_", "\\_"), r[2]] for r in prs[:15]]
    tbl("tab_pairs.tbl", ["l l c"], "l l n", [["host genus", "genus sharing the sequence", "genomes"]], rows)
    # (11a) the pre-specified rule against the final rule (host-equivalent references)
    p0 = os.path.join(RES, "prespecified", "summary.json")
    if os.path.exists(p0):
        S0 = json.load(open(p0))
        rows = []
        for k, lab in (("isolates", "isolates"), ("mags", "MAGs"), ("hq_mags", "\\0\\0high-quality MAGs")):
            a, b = S0["groups"][k], g[k]
            rows.append([lab, f"{b['n']:,}", pct(a["contam_any"][0]), pct(b["contam_any"][0]), pct(a["contam_ge1"][0], 2),
                         pct(b["contam_ge1"][0], 2), pct(a["contam_ge5"][0], 2), pct(b["contam_ge5"][0], 2)])
        a, b = S0["species_matched"]["contam_any"], S["species_matched"]["contam_any"]
        rows.append(["MAG \\(mi isolate, same species (pp)", "", f"{100 * a['mh_rd']:+.1f}", f"{100 * b['mh_rd']:+.1f}", "", "", "", ""])
        a, b = S0["core_mixture"], S["core_mixture"]
        rows.append(["core-gene mixture, all lengths (%)", "", pct(a["overall"], 0), pct(b["overall"], 0), "", "", "", ""])
        tbl("tab_rule.tbl", ["l c c s c s c s", "^ ^ c c c c c c"], "l n n n n n n n",
            [["", "n", "contamination, any (%)", "", "\\(>= 1% (%)", "", "\\(>= 5% (%)", ""],
             ["", "", "pre-specified", "final", "pre-specified", "final", "pre-specified", "final"]], rows)
        nums["prespecified"] = {"groups": {k: {x: S0["groups"][k][x] for x in ("n", "contam_any", "contam_ge1", "contam_ge5")}
                                           for k in ("isolates", "mags", "hq_mags")},
                                "mixture": S0["core_mixture"]["overall"], "matched": S0["species_matched"]["contam_any"],
                                "flagged": S0["flagged_contigs"]}
    # (11b) GTDB references conspecific with UHGG genomes of another genus label
    _, lr = read_tsv("label_inconsistent_refs.tsv")
    rows = [[r[0], r[1].replace("_", "\\_"), r[3].replace("_", "\\_"), r[4], f"{int(r[6]):,}",
             r[7].split(";")[0].split(":")[0].replace("_", "\\_")] for r in lr]
    tbl("tab_labels.tbl", ["l l l c c c", "l ^ ^ c c c"], "l l l n n l",
        [["GTDB R232", "GTDB species", "NCBI species", "CheckM2", "UHGG", "most frequent"],
         ["representative", "", "", "compl. (%)", "genomes", "host genus"]], rows)
    nums["labels"] = S["label_inconsistent"]
    # (11c) contamination against MAG fragmentation and completeness
    rows = []
    for d in S["by_contigs"]:
        lab = f"{d['lo']}\\(en{d['hi'] - 1}" if d["hi"] < 10 ** 8 else f"\\(>= {d['lo']}"
        rows.append(["contigs" if d is S["by_contigs"][0] else "", lab, f"{d['n']:,}", ci(d["contam_ge1"], 2), pct(d["foreign_ge5"][0])])
    rows.append("_")
    for d in S["by_completeness"]:
        rows.append(["completeness (%)" if d is S["by_completeness"][0] else "", f"{d['lo']:g}\\(en{d['hi']:g}", f"{d['n']:,}",
                     ci(d["contam_ge1"], 2), pct(d["foreign_ge5"][0])])
    tbl("tab_frag.tbl", ["l l c c c", "^ ^ ^ c c"], "l l n n n",
        [["MAGs by", "", "n", "contamination \\(>= 1%", "non-mobile shared"], ["", "", "", "(%, 95% CI)", "\\(>= 5% (%)"]], rows)
    # (11) AMR
    nums["amr"] = S["amr_reps"]
    nums["studies"] = S.get("studies")
    nums["euk"] = S["euk_known_issue"]
    nums["low_hit"] = S["low_hit"]
    nums["host_vs_mgnify"] = S["host_vs_mgnify"]
    nums["extra"] = extra_numbers(lr)
    # figure 1: distribution of the per-genome fractions
    E = S["hist"]["edges"]
    labs = ["0", "<0.01", "0.01–0.1", "0.1–0.3", "0.3–1", "1–2", "2–5", "5–10", "10–20", "≥20"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4), sharey=True)
    for ax, key, title in ((axes[0], "foreign", "non-mobile shared sequence"), (axes[1], "contam", "core-anchored contamination")):
        for j, (t, col, lab) in enumerate((("isolates", C_ISO, "isolates"), ("mags", C_MAG, "MAGs"))):
            h = S["hist"][f"{key}_{t}"]
            n = sum(h)
            ax.bar([i + (j - 0.5) * 0.38 for i in range(len(h))], [100 * x / n for x in h], width=0.36, color=col, label=lab)
        ax.set_yscale("log")
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=45, fontsize=6.3, ha="right")
        ax.set_xlabel("per-genome fraction (%)", fontsize=7.5)
        ax.set_title(title, fontsize=8)
        axes_style(ax)
    axes[0].set_ylabel("genomes (%)", fontsize=7.5)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_hist.eps"))
    plt.close(fig)
    # figure 2: core-gene deficit and the CheckM2 gradient
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4))
    ax = axes[0]
    xs = range(4)
    for j, (key, col, lab) in enumerate((("MAG_bg", C_BG, "MAG, unflagged"), ("MAG_flag_chrom", C_MAG, "MAG, flagged non-mobile"),
                                          ("Isolate_flag_chrom", C_ISO, "isolate, flagged non-mobile"))):
        ax.bar([x + (j - 1) * 0.27 for x in xs], [P[f"{key}_{b}"]["core_per_kb"] for b in xs], width=0.25, color=col, label=lab)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(LBP)
    ax.set_xlabel("contig length (kb)", fontsize=7.5)
    ax.set_ylabel("core genes per kb", fontsize=7.5)
    ax.legend(frameon=False, fontsize=6.3)
    axes_style(ax)
    ax = axes[1]
    cb = S["by_checkm_contam"]
    x = [(d["lo"] + d["hi"]) / 2 for d in cb]
    y = [100 * d["contam_any"][0] for d in cb]
    lo = [100 * (d["contam_any"][0] - d["contam_any"][1]) for d in cb]
    hi = [100 * (d["contam_any"][2] - d["contam_any"][0]) for d in cb]
    ax.errorbar(x, y, yerr=[lo, hi], fmt="o-", color=C_MAG, ms=4, lw=1.2, capsize=2)
    ax.set_xlabel("CheckM2 contamination reported by MGnify (%)", fontsize=7.5)
    ax.set_ylabel("MAGs with core-anchored\ncontamination (%)", fontsize=7.5)
    axes_style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_core.eps"))
    plt.close(fig)
    with open(os.path.join(GEN, "numbers.txt"), "w") as fh:
        for k, v in nums.items():
            fh.write(f"{k}\t{json.dumps(v)}\n")
    print(open(os.path.join(GEN, "numbers.txt")).read()[:6000])


if __name__ == "__main__":
    main()
