#!/usr/bin/env python3
"""Tables (tbl), figures (EPS) and key numbers for the real-genome report.

Reads data/real/... (bench_real/ outputs) and writes report_real/gen/*.tbl,
report_real/gen/numbers.txt and report_real/figs/*.eps. Also re-runs a few tuning points
(final binary) on the real and simulated tuning sets for the tuning table.
Usage: python3 scripts/real_report.py   (from the repository root)
"""
import csv
import json
import math
import os
import subprocess
import sys
from collections import OrderedDict

import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench_real"))
import evaluate_real as ev  # noqa: E402

R = os.path.join(ROOT, "data", "real")
RT = os.path.join(ROOT, "data", "real_tune")
GEN = os.path.join(ROOT, "report_real", "gen")
FIG = os.path.join(ROOT, "report_real", "figs")
AIF = os.path.join(ROOT, "bin", "aif")
SETS = ["bcereus", "ecoli", "mag2328", "paer", "saur"]
SNAME = {"bcereus": "B. cereus group (D2)", "ecoli": "E. coli (D3)", "mag2328": "MAG bin 2328",
         "paer": "P. aeruginosa", "saur": "S. aureus", "negative": "negative controls"}
SSHORT = {"bcereus": "B. cereus", "ecoli": "E. coli", "mag2328": "MAGs 2328", "paer": "P. aerug.", "saur": "S. aureus"}
REFNAME = {"BA52G": "B. anthracis 52-G chr.", "MG1655": "E. coli K-12 MG1655", "REFMAG": "best MAG of bin 2328",
           "PAO1": "P. aeruginosa PAO1", "NCTC8325": "S. aureus NCTC 8325"}
# entity colours (dataviz reference palette, light mode): fixed per tool across all figures
COL = {"aif": "#2a78d6", "skani": "#eb6834", "aif_rt": "#1baf7a", "skani_slow": "#eda100"}
LAB = {"aif": "aif", "skani": "skani", "aif_rt": "aif-RT", "skani_slow": "skani slow"}  # "slow" = skani --slow
CMAP = {"aif": "Blues", "skani": "Oranges", "aif_rt": "Greens", "skani_slow": "YlOrBr"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2})


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
    open(os.path.join(GEN, path), "w").write("\n".join(lines) + "\n")


def f(x, d=2):
    return "\\(em" if x is None or x != x else f"{x:.{d}f}"


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    nums = OrderedDict()
    E = json.load(open(os.path.join(R, "results", "eval_main.json")))
    S = E["sets"]

    # ---------------- data table
    rows = []
    tot_g = tot_c = tot_b = 0
    for s in SETS + ["negative"]:
        g = tsv(os.path.join(R, "sets", s, "genomes.tsv"))
        t = tsv(os.path.join(R, "sets", s, "truth.tsv"))
        gt = tsv(os.path.join(R, "sets", s, "genome_truth.tsv"))
        L = sorted(int(r["length"]) for r in t)
        hom = sum(1 for r in t if float(r["true_af"]) >= 0.2) / len(t)
        anis = [float(r["anim_ani"]) * 100 for r in gt]
        refs = sorted(set(r["reference"] for r in g))
        nc = sum(int(r["contigs"]) for r in g)
        nb = sum(int(r["bases"]) for r in g)
        tot_g += len(g)
        tot_c += nc
        tot_b += nb
        rows.append([SNAME[s], ", ".join(REFNAME.get(x, x) for x in refs) if s != "negative" else "BA52G / MG1655 (crossed)",
                     str(len(g)), f"{nc:,}", f"{nb/1e6:.0f}", f"{L[len(L)//2]/1000:.1f}",
                     f"{min(anis):.1f}\\(en{max(anis):.1f}", f"{hom*100:.0f}"])
    rows.append("_")
    rows.append(["total", "", str(tot_g), f"{tot_c:,}", f"{tot_b/1e6:.0f}", "", "", ""])
    nums.update(n_genomes=tot_g, n_contigs=tot_c, n_mbp=round(tot_b / 1e6))
    tbl("tab_data.tbl", ["l l c c c c c c", "l l c c c c c c"], "l l n n n n c n",
        [["set", "reference", "genomes", "contigs", "Mbp", "median", "genome ANIm", "homol."],
         ["", "", "", "\\(>= 500 bp", "", "len. (kb)", "range (%)", "(%)"]], rows)

    # ---------------- truth validation: our ANIm vs published D2 ANIm (skani-test D2_anim.txt)
    pub = {}
    for line in open(os.path.join(R, "meta", "D2_anim.txt")):
        p = line.rstrip("\n").split("\t")
        if len(p) > 2:
            pub[p[1].replace(".LargeContigs.fna", "")] = float(p[2])
    gmap = {r["gid"]: r["genome"] for r in tsv(os.path.join(R, "sets", "bcereus", "genomes.tsv"))}
    ours = {r["gid"]: float(r["anim_ani"]) * 100 for r in tsv(os.path.join(R, "sets", "bcereus", "genome_truth.tsv"))}
    d = [ours[g] - pub[gmap[g]] for g in ours if gmap[g] in pub]
    nums.update(anim_vs_pub_n=len(d), anim_vs_pub_mae=sum(abs(x) for x in d) / len(d),
                anim_vs_pub_max=max(abs(x) for x in d), anim_vs_pub_bias=sum(d) / len(d))

    # ---------------- sketch sizes
    sz = tsv(os.path.join(R, "results", "sketch_sizes.tsv"))
    tb = sum(int(r["bases"]) for r in sz)
    nums["aif_bpb"] = sum(int(r["aif_bytes"]) for r in sz) * 8 / tb
    nums["skani_bpb"] = sum(int(r["skani_bytes"]) for r in sz) * 8 / tb
    rows = [[REFNAME[r["reference"]], f"{int(r['bases'])/1e6:.2f}", f"{int(r['skani_bytes'])/1024:.0f}",
             f"{int(r['aif_bytes'])/1024:.0f}", f(float(r["skani_bits_per_base"]), 3), f(float(r["aif_bits_per_base"]), 3)]
            for r in sz]
    tbl("tab_sizes.tbl", ["l c c c c c"], "l n n n n n",
        [["reference", "Mbp", "skani KiB", "aif KiB", "skani bits/base", "aif bits/base"]], rows)

    # ---------------- timing
    tim = {}
    for r in tsv(os.path.join(R, "results", "timing_real.tsv")):
        tim[(r["set"], r["tool"], r["phase"], int(r["threads"]))] = (float(r["median_s"]), int(r["max_rss_kb"]) / 1024)
    t_sk1 = sum(tim[(s, "skani", "sketch+query", 1)][0] for s in SETS)
    t_ai1 = sum(tim[(s, "aif", "sketch", 1)][0] + tim[(s, "aif", "query", 1)][0] for s in SETS)
    t_af1 = sum(tim[(s, "aif", "sketch", 1)][0] + tim[(s, "aif", "query-fast", 1)][0] for s in SETS)
    t_sk8 = sum(tim[(s, "skani", "sketch+query", 8)][0] for s in SETS)
    t_ai8 = sum(tim[(s, "aif", "sketch", 1)][0] + tim[(s, "aif", "query", 8)][0] for s in SETS)
    t_af8 = sum(tim[(s, "aif", "sketch", 1)][0] + tim[(s, "aif", "query-fast", 8)][0] for s in SETS)
    nums.update(t_sk1=t_sk1, t_ai1=t_ai1, t_af1=t_af1, t_sk8=t_sk8, t_ai8=t_ai8, t_af8=t_af8)
    rows = []
    for s in SETS:
        nb = sum(int(r["bases"]) for r in tsv(os.path.join(R, "sets", s, "genomes.tsv")))
        sk1, m_sk1 = tim[(s, "skani", "sketch+query", 1)]
        a_s = tim[(s, "aif", "sketch", 1)][0]
        aq1, m_aq1 = tim[(s, "aif", "query", 1)]
        af1 = tim[(s, "aif", "query-fast", 1)][0]
        sk8, m_sk8 = tim[(s, "skani", "sketch+query", 8)]
        aq8, m_aq8 = tim[(s, "aif", "query", 8)]
        af8 = tim[(s, "aif", "query-fast", 8)][0]
        rows.append([SSHORT[s], f"{nb/1e6:.0f}", f(sk1), f(a_s + aq1), f(a_s + af1), f(sk8), f(a_s + aq8), f(a_s + af8),
                     f"{m_sk1:.0f} / {m_sk8:.0f}", f"{m_aq1:.0f} / {m_aq8:.0f}"])
    rows.append("_")
    rows.append(["total", "", f(t_sk1), f(t_ai1), f(t_af1), f(t_sk8), f(t_ai8), f(t_af8), "", ""])
    tbl("tab_speed.tbl", ["l c c s s c s s c c", "^ ^ c c c c c c c c"], "l n n n n n n n c c",
        [["set", "Mbp", "1 thread (s)", "", "", "8 threads (s)", "", "", "peak RSS (MB)", ""],
         ["", "", "skani", "aif", "aif fast", "skani", "aif", "aif fast", "skani 1t / 8t", "aif 1t / 8t"]], rows)
    sc = tsv(os.path.join(R, "results", "scaling.tsv"))
    scd = {(r["version"], int(r["threads"])): (float(r["median_s"]), int(r["max_rss_kb"]) / 1024) for r in sc}
    thr = sorted(set(int(r["threads"]) for r in sc))
    vers = sorted(set(r["version"] for r in sc))
    rows = []
    for v in vers:
        rows.append([v] + [f"{scd[(v, t)][0]:.2f} ({scd[(v, t)][1]:.0f})" for t in thr])
    tbl("tab_scaling.tbl", ["l " + " ".join("c" for _ in thr)], "l " + " ".join("c" for _ in thr),
        [["aif query"] + [f"{t} thread{'s' if t > 1 else ''}" for t in thr]], rows)
    for v in vers:
        for t in thr:
            nums[f"sc_{v}_{t}"] = scd[(v, t)][0]
            nums[f"sc_{v}_{t}_rss"] = scd[(v, t)][1]

    # ---------------- genome level
    G = json.load(open(os.path.join(R, "results", "genome_level.json")))
    gtools = [("skani", "skani (genome)"), ("skani_slow", "skani --slow (genome)"), ("fastANI", "fastANI"),
              ("skani_qi_agg", "skani --qi (agg.)"), ("aif_agg", "aif (agg.)"), ("aif_rt_agg", "aif-RT (agg.)"),
              ("aif_fast_agg", "aif fast (agg.)")]
    rows = []
    for t, lab in gtools:
        rows.append([lab] + [f(G[s][t]["ani_mae"], 3) for s in SETS] + [f(G["all"][t]["ani_mae"], 3),
                                                                        f(G["all"][t]["af_mae"], 2), f(G["all"][t]["af_mae_q"], 2)])
    tbl("tab_genome.tbl", ["l c s s s s s c c", "^ c c c c c c c c"], "l n n n n n n n n",
        [["method", "genome ANI MAE vs ANIm (pp)", "", "", "", "", "", "AF MAE (pp)", ""],
         ["", "B. cer.", "E. coli", "MAGs", "P. aer.", "S. aur.", "all 224", "1-to-1", "many-to-1"]], rows)
    for t, _ in gtools:
        nums[f"g_{t}_ani"] = G["all"][t]["ani_mae"]
        nums[f"g_{t}_af"] = G["all"][t]["af_mae"]
        nums[f"g_{t}_afq"] = G["all"][t]["af_mae_q"]

    # ---------------- headline and per-set tables (contig level)
    C4 = ["skani", "skani_slow", "aif", "aif_rt"]
    W = S["all_within"]
    M = E["macro"]
    for k in C4 + ["aif_0.44", "aif_fast", "skani_relaxed"]:
        nums[f"{k}_cov"] = W[k]["coverage"]
        nums[f"{k}_ani"] = W[k]["ani_mae"]
        nums[f"{k}_pani"] = W[k]["paired_ani_mae"]
        nums[f"{k}_af"] = W[k]["af_mae"]
        nums[f"{k}_fp"] = W[k]["fp_rate"]
        nums[f"{k}_mani"] = M[k]["ani_mae"]
        nums[f"{k}_maf"] = M[k]["af_mae"]
    nums["paired_n"] = W["aif"]["paired_n"]
    nums["n_hom"] = W["aif"]["n_hom"]
    nums["n_within"] = W["aif"]["n"]
    nums["n_non"] = W["aif"]["n_non"]
    wins = {k: sum(1 for s in SETS if S[s]["aif"]["ani_mae"] < S[s][k]["ani_mae"] and S[s]["aif"]["af_mae"] < S[s][k]["af_mae"])
            for k in ("skani", "skani_slow")}
    nums.update(wins_skani=wins["skani"], wins_slow=wins["skani_slow"])
    hrows = [
        ["sketch size (bits / reference base)", f(nums["skani_bpb"], 3), "3.33", f(nums["aif_bpb"], 3), f(nums["aif_bpb"], 3)],
        ["homologous contigs reported (coverage)"] + [f(W[k]["coverage"], 3) for k in C4],
        ["ANI MAE, reported homologous contigs (pp)"] + [f(W[k]["ani_mae"], 3) for k in C4],
        [f"ANI MAE, contigs reported by all ({W['aif']['paired_n']:,}) (pp)"] + [f(W[k]["paired_ani_mae"], 3) for k in C4],
        ["ANI MAE, mean over the 5 sets (pp)"] + [f(M[k]["ani_mae"], 3) for k in C4],
        ["AF MAE, all contigs (pp)"] + [f(W[k]["af_mae"], 2) for k in C4],
        ["AF MAE, mean over the 5 sets (pp)"] + [f(M[k]["af_mae"], 2) for k in C4],
        ["false homology calls (% of non-homologous)"] + [f(W[k]["fp_rate"] * 100, 2) for k in C4],
        ["genome ANI MAE vs ANIm, 224 genomes (pp)", f(G["all"]["skani"]["ani_mae"], 3), f(G["all"]["skani_slow"]["ani_mae"], 3),
         f(G["all"]["aif_agg"]["ani_mae"], 3), f(G["all"]["aif_rt_agg"]["ani_mae"], 3)],
        ["wall time, 5 sets, 1 / 8 threads (s)", f"{t_sk1:.2f} / {t_sk8:.2f}", "\\(em", f"{t_ai1:.2f} / {t_ai8:.2f}", "\\(em"],
    ]
    tbl("tab_headline.tbl", ["l c c c c"], "l n n n n",
        [["contig level (33,464 contigs, 5 within-species sets)", "skani", "skani --slow", "aif", "aif-RT"]], hrows)
    rows = []
    for s in SETS + ["negative"]:
        r = [SNAME[s], f"{S[s]['aif']['n']:,}"]
        for k in C4:
            x = S[s][k]
            if s == "negative":
                r += ["\\(em", f(x["af_mae"], 2)]
            else:
                r += [f(x["ani_mae"], 3), f(x["af_mae"], 2)]
        rows.append(r)
    tbl("tab_sets.tbl", ["l c c s c s c s c s", "^ ^ c c c c c c c c"], "l n n n n n n n n n",
        [["set", "contigs", "skani", "", "skani --slow", "", "aif", "", "aif-RT", ""],
         ["", "", "ANI", "AF", "ANI", "AF", "ANI", "AF", "ANI", "AF"]], rows)
    rows = []
    for s in SETS:
        rows.append([SNAME[s]] + [f(S[s][k]["coverage"], 3) for k in C4] + [f(S[s][k]["paired_ani_mae"], 3) for k in C4])
    tbl("tab_cov.tbl", ["l c s s s c s s s", "^ c c c c c c c c"], "l n n n n n n n n",
        [["set", "coverage of homologous contigs", "", "", "", "paired ANI MAE (pp)", "", "", ""],
         ["", "skani", "--slow", "aif", "aif-RT", "skani", "--slow", "aif", "aif-RT"]], rows)

    # ---------------- strata tables
    def strata_table(name, labels, fname, ani_only=False):
        rows = []
        for (b, m), lab in zip(E["strata"][name].items(), labels):
            r = [lab, f"{m['aif']['n']:,}"]
            for k in C4:
                r += [f(m[k]["coverage"], 2) if not ani_only else f(m[k]["coverage"], 2), f(m[k]["ani_mae"], 2), f(m[k]["af_mae"], 1)]
            rows.append(r)
        tbl(fname, ["l c c s s c s s c s s c s s", "^ ^ c c c c c c c c c c c c"], "l n n n n n n n n n n n n n",
            [["stratum", "contigs", "skani", "", "", "skani --slow", "", "", "aif", "", "", "aif-RT", "", ""],
             ["", "", "cov", "ANI", "AF", "cov", "ANI", "AF", "cov", "ANI", "AF", "cov", "ANI", "AF"]], rows)
    strata_table("length", ["0.5\\(en1 kb", "1\\(en2 kb", "2\\(en5 kb", "5\\(en20 kb", "\\(>= 20 kb"], "tab_len.tbl")
    strata_table("true_ani", ["< 90%", "90\\(en95%", "95\\(en98%", "98\\(en99.5%", "\\(>= 99.5%"], "tab_ani.tbl")
    strata_table("true_af", ["< 0.05", "0.05\\(en0.5", "0.5\\(en0.9", "\\(>= 0.9"], "tab_af.tbl")
    bias_rows = []
    for (b, m), lab in zip(E["strata"]["true_ani"].items(), ["< 90%", "90\\(en95%", "95\\(en98%", "98\\(en99.5%", "\\(>= 99.5%"]):
        bias_rows.append([lab] + [f"{m[k]['ani_bias']:+.2f}" for k in C4] + [f"{m[k]['af_bias']:+.1f}" for k in C4])
    tbl("tab_bias.tbl", ["l c s s s c s s s", "^ c c c c c c c c"], "l n n n n n n n n",
        [["true ANI", "ANI bias (pp)", "", "", "", "AF bias (pp)", "", "", ""],
         ["", "skani", "--slow", "aif", "aif-RT", "skani", "--slow", "aif", "aif-RT"]], bias_rows)
    st = E["strata"]["length"]
    first = list(st.values())[0]
    nums.update(short_cov_aif=first["aif"]["coverage"], short_cov_sk=first["skani"]["coverage"],
                short_cov_slow=first["skani_slow"]["coverage"], short_af_aif=first["aif"]["af_mae"],
                short_af_sk=first["skani"]["af_mae"], short_pani_aif=first["aif"]["paired_ani_mae"],
                short_pani_sk=first["skani"]["paired_ani_mae"])
    sa = E["strata"]["true_af"]
    part = list(sa.values())[1]
    nums.update(part_af_aif=part["aif"]["af_mae"], part_af_sk=part["skani"]["af_mae"], part_af_slow=part["skani_slow"]["af_mae"],
                part_af_rt=part["aif_rt"]["af_mae"], part_bias_aif=part["aif"]["af_bias"])
    low = list(E["strata"]["true_ani"].values())[0]
    nums.update(low_ani_aif=low["aif"]["ani_mae"], low_ani_sk=low["skani"]["ani_mae"], low_cov_aif=low["aif"]["coverage"],
                low_cov_sk=low["skani"]["coverage"], low_bias_aif=low["aif"]["ani_bias"], low_bias_sk=low["skani"]["ani_bias"])

    # ---------------- truth sensitivity: --mumreference truth vs --maxmatch truth
    truth_mm, truth_mr, members = {}, {}, {}
    preds = {k: {} for k in ("skani", "aif")}
    for s in SETS:
        t = ev.load_truth(s)
        truth_mm.update(t)
        members[s] = sorted(t)
        for r in tsv(os.path.join(R, "truth_mumref", s + ".truth.tsv")):
            truth_mr[r["contig"]] = (int(r["length"]), float(r["true_af"]), float(r["true_ani"]))
        for k in preds:
            preds[k].update(ev.load_pred(k, s))
    allc = [c for s in SETS for c in members[s]]
    a = ev.score(allc, truth_mr, preds, ["skani", "aif"])
    b = ev.score(allc, truth_mm, preds, ["skani", "aif"])
    changed = sum(1 for c in allc if truth_mm[c][1] - truth_mr[c][1] > 0.05)
    nums["truth_changed"] = changed
    rows = []
    for lab, sc_ in (("nucmer default seeds (unique anchors)", a), ("nucmer --maxmatch (used)", b)):
        rows.append([lab] + [f(sc_[k]["ani_mae"], 3) for k in ("skani", "aif")] + [f(sc_[k]["af_mae"], 2) for k in ("skani", "aif")]
                    + [f(sc_[k]["any_report_rate"] * 100, 1) for k in ("skani", "aif")])
    tbl("tab_truth.tbl", ["l c s c s c s", "^ c c c c c c"], "l n n n n n n",
        [["per-contig truth", "ANI MAE (pp)", "", "AF MAE (pp)", "", "reports on AF<0.05 (%)", ""],
         ["", "skani", "aif", "skani", "aif", "skani", "aif"]], rows)
    nums.update(mr_aif_af=a["aif"]["af_mae"], mm_aif_af=b["aif"]["af_mae"], mr_aif_any=a["aif"]["any_report_rate"],
                mm_aif_any=b["aif"]["any_report_rate"])
    BA = json.load(open(os.path.join(R, "results", "blast_audit.json")))
    rows = []
    for g in ("disputed", "control"):
        x = BA[g]
        rows.append([g + f" ({x['n']})", f(x["nucmer_af"], 3), f(x["blastn_af"], 3), f(x["aif_af"], 3),
                     f(x["aif_vs_nucmer_mae"], 1), f(x["aif_vs_blastn_mae"], 1)])
    tbl("tab_audit.tbl", ["l c c c c c"], "l n n n n n",
        [["sample", "nucmer AF", "blastn AF", "aif AF", "|aif\\(minnucmer| (pp)", "|aif\\(minblastn| (pp)"]], rows)
    nums.update(audit_dis_total=BA["n_disputed_total"], audit_dis_nuc=BA["disputed"]["nucmer_af"],
                audit_dis_bl=BA["disputed"]["blastn_af"], audit_dis_aif=BA["disputed"]["aif_af"])

    # ---------------- tuning table (re-run with the final binary)
    def real_tune(args):
        os.environ["REAL_DIR"] = RT
        truth, mem = {}, {}
        pr = {}
        for s in SETS:
            tt = {r["contig"]: (int(r["length"]), float(r["true_af"]), float(r["true_ani"]))
                  for r in tsv(os.path.join(RT, "sets", s, "truth.tsv"))}
            truth.update(tt)
            mem[s] = sorted(tt)
            refs = sorted(set(r["reference"] for r in tsv(os.path.join(RT, "sets", s, "genomes.tsv"))))
            for ref in refs:
                out = subprocess.run([AIF, "query", os.path.join(RT, "results", "aif", f"{s}.{ref}.aif"),
                                      os.path.join(RT, "sets", s, f"{ref}.contigs.fa"), "-t", "8"] + args,
                                     capture_output=True, text=True, check=True).stdout
                for line in out.splitlines()[1:]:
                    p = line.split("\t")
                    ok = p[7] == "ok"
                    pr[p[0]] = (float(p[1]) if ok else math.nan, float(p[2]) if ok else 0.0, int(ok))
        m = {s: ev.score(mem[s], truth, {"x": pr}, ["x"])["x"] for s in SETS}
        return (sum(m[s]["ani_mae"] for s in SETS) / 5, sum(m[s]["af_mae"] for s in SETS) / 5)

    def sim_tune(args):
        sys.path.insert(0, os.path.join(ROOT, "bench"))
        truth = tsv(os.path.join(ROOT, "data", "tune", "truth.tsv"))
        pr = {}
        for g in ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NC_007795.1"]:
            sk = os.path.join(ROOT, "data", "results", "main", "aif", g + ".aif")
            out = subprocess.run([AIF, "query", sk, os.path.join(ROOT, "data", "tune", g + ".contigs.fa"), "-t", "8"] + args,
                                 capture_output=True, text=True, check=True).stdout
            for line in out.splitlines()[1:]:
                p = line.split("\t")
                if p[7] == "ok":
                    pr[p[0].split()[0]] = (float(p[1]), float(p[2]))
        rep = [t for t in truth if t["contig"] in pr]
        ani = sum(abs(pr[t["contig"]][0] - float(t["true_ani"])) for t in rep) / len(rep) * 100
        af = sum(abs(pr[t["contig"]][1] - float(t["true_af"])) for t in rep) / len(rep) * 100
        return ani, af

    rows = []
    for ru in ["1e-6", "3e-6", "1e-5", "2e-5", "3e-5", "1e-4"]:
        args = ["--rho-u", ru]
        ra, rf = real_tune(args)
        sa_, sf = sim_tune(args)
        lab = ru + (" (default)" if ru == "3e-6" else " (aif-RT)" if ru == "2e-5" else "")
        rows.append([lab, f(ra, 3), f(rf, 2), f(sa_, 3), f(sf, 2)])
        nums[f"tune_{ru}"] = (ra, rf, sa_, sf)
    tbl("tab_tune.tbl", ["l c s c s", "^ c c c c"], "l n n n n",
        [["rho_u", "real tuning set (5 sets)", "", "simulated tuning set", ""],
         ["", "ANI MAE", "AF MAE", "ANI MAE", "AF MAE"]], rows)

    # ---------------- figures
    # (1) per-set grouped bars: ANI MAE and AF MAE
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))
    for ax, key, lab in ((axes[0], "ani_mae", "ANI MAE, homologous contigs (pp)"), (axes[1], "af_mae", "AF MAE, all contigs (pp)")):
        w = 0.19
        for i, k in enumerate(C4):
            xs = [j + (i - 1.5) * (w + 0.02) for j in range(len(SETS))]
            ax.bar(xs, [S[s][k][key] for s in SETS], width=w, color=COL[k], label=LAB[k], edgecolor="white", linewidth=0.6)
        ax.set_xticks(range(len(SETS)))
        ax.set_xticklabels([SSHORT[s] for s in SETS], fontsize=7, color=INK)
        ax.set_title(lab, fontsize=8, color=INK)
        ax.grid(True, axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].legend(frameon=False, fontsize=6.8, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_sets.eps"))
    plt.close(fig)

    # (2) by contig length: coverage, ANI MAE, AF MAE
    labs = ["0.5–1k", "1–2k", "2–5k", "5–20k", "≥20k"]
    fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.2))
    for ax, key, lab in ((axes[0], "coverage", "coverage (homologous)"), (axes[1], "paired_ani_mae", "paired ANI MAE (pp)"),
                         (axes[2], "af_mae", "AF MAE (pp)")):
        for k, mk in zip(C4, ["o", "D", "s", "^"]):
            ys = [m[k][key] for m in st.values()]
            ax.plot(range(5), ys, "-" + mk, color=COL[k], lw=1.6, ms=4, label=LAB[k])
        ax.set_xticks(range(5))
        ax.set_xticklabels(labs, fontsize=6.5)
        ax.set_title(lab, fontsize=8, color=INK)
        ax.grid(True, color=GRID, lw=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[2].legend(frameon=False, fontsize=6.5)
    axes[1].set_xlabel("contig length (bp)", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_length.eps"))
    plt.close(fig)

    # (3) AF: estimate vs truth, small multiples (one tool per panel)
    truth_all = {}
    for s in SETS:
        truth_all.update(ev.load_truth(s))
    P = {k: {} for k in C4}
    for k in C4:
        for s in SETS:
            P[k].update(ev.load_pred(k, s))
    names = sorted(truth_all)
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 1.95), sharey=True)
    for ax, k in zip(axes, C4):
        xs = [truth_all[c][1] for c in names]
        ys = [P[k][c][1] for c in names]
        ax.hexbin(xs, ys, gridsize=40, extent=(0, 1, 0, 1), bins="log", cmap=CMAP[k], mincnt=1, linewidths=0.1)
        ax.plot([0, 1], [0, 1], color=INK2, lw=0.6)
        ax.set_title(f"{LAB[k]} (MAE {W[k]['af_mae']:.2f} pp)", fontsize=7.5, color=INK)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("true AF (nucmer)", fontsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("estimated AF", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_af_scatter.eps"))
    plt.close(fig)

    # (4) ANI signed error vs true ANI (homologous contigs), skani vs aif
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.3), sharey=True)
    hom = [c for c in names if truth_all[c][1] >= 0.2]
    for ax, k in zip(axes, ["skani", "aif"]):
        xs = [truth_all[c][2] * 100 for c in hom if P[k][c][2]]
        ys = [(P[k][c][0] - truth_all[c][2]) * 100 for c in hom if P[k][c][2]]
        ax.hexbin(xs, ys, gridsize=(50, 30), extent=(80, 100, -6, 6), bins="log", cmap=CMAP[k], mincnt=1, linewidths=0.1)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_xlim(80, 100.2)
        ax.set_ylim(-6, 6)
        ax.set_title(f"{LAB[k]}: {len(xs):,} of {len(hom):,} homologous contigs", fontsize=7.5, color=INK)
        ax.set_xlabel("true contig ANI (%)", fontsize=7)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("estimate − truth (pp)", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_ani_err.eps"))
    plt.close(fig)

    # (5) thread scaling
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    vstyle = {"v1.1": (INK2, "o", ":"), "v1.2": ("#8fb8ea", "^", "-"), "v1.3": (COL["aif"], "s", "-")}
    for v in vers:
        col, mk, ls = vstyle.get(v, (INK2, "o", "-"))
        ax.plot(thr, [scd[(v, t)][0] for t in thr], ls + mk, color=col, lw=1.6, ms=4, label=f"aif {v}")
    sk1 = tim[("bcereus", "skani", "sketch+query", 1)][0]
    sk8 = tim[("bcereus", "skani", "sketch+query", 8)][0]
    ax.plot([1, 8], [sk1, sk8], "--o", color=COL["skani"], lw=1.4, ms=4, label="skani (sketch+query)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(thr)
    ax.set_xticklabels([str(t) for t in thr])
    yt = [0.3, 0.5, 1, 2, 4]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{y:g}" for y in yt])
    ax.minorticks_off()
    ax.set_xlabel("threads", fontsize=7.5)
    ax.set_ylabel("wall time (s), 449 Mbp", fontsize=7.5)
    ax.grid(True, which="both", color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_scaling.eps"))
    plt.close(fig)

    # ================= v1.4: confidence intervals =================
    CAL = json.load(open(os.path.join(R, "results", "calibration.json")))
    lvls = ["0.9", "0.95", "0.99"]
    DSN = {"real_tune": "real, tuning (fit)", "real_eval": "real, held out", "sim_main": "simulated main", "sim_fresh": "simulated fresh"}
    rows = []
    for tag in ["real_tune", "real_eval", "sim_main", "sim_fresh"]:
        for mode, lab in (("stat", "statistical"), ("calibrated", "calibrated")):
            r = [DSN[tag], lab]
            r += [f"{CAL[mode][l][tag]['ani'] * 100:.1f}" for l in lvls]
            r += [f"{CAL[mode][l][tag]['af'] * 100:.1f}" for l in lvls]
            r += [f"{CAL[mode]['0.95'][tag]['halfwidth_ani_med']:.2f}", f"{CAL[mode]['0.95'][tag]['halfwidth_af_med']:.2f}"]
            rows.append(r)
        if tag != "sim_fresh":
            rows.append("_")
    tbl("tab_ci.tbl", ["l l c s s c s s c s", "^ ^ c c c c c c c c"], "l l n n n n n n n n",
        [["data", "intervals", "ANI coverage (%)", "", "", "AF coverage (%)", "", "", "median half-width", ""],
         ["", "", "90%", "95%", "99%", "90%", "95%", "99%", "ANI (pp)", "AF (pp)"]], rows)
    for tag in ["real_eval", "sim_main"]:
        for mode in ("stat", "calibrated"):
            for l in lvls:
                nums[f"ci_{mode}_{tag}_{l}_ani"] = CAL[mode][l][tag]["ani"]
                nums[f"ci_{mode}_{tag}_{l}_af"] = CAL[mode][l][tag]["af"]
            nums[f"ci_{mode}_{tag}_hw_ani"] = CAL[mode]["0.95"][tag]["halfwidth_ani_med"]
            nums[f"ci_{mode}_{tag}_hw_af"] = CAL[mode]["0.95"][tag]["halfwidth_af_med"]
    # statistical coverage by simulated condition family
    grp = CAL["stat"]["0.95"]["sim_main"]["groups"]
    fam = OrderedDict([("substitutions only, AF = 1", lambda c: c.startswith(("len", "ani"))),
                       ("partial AF", lambda c: c.startswith("af")),
                       ("indels", lambda c: c.startswith("indel")),
                       ("rate heterogeneity", lambda c: c.startswith("het")),
                       ("chimeras", lambda c: c.startswith("chimera"))])
    rows = []
    for name, sel in fam.items():
        cs = [g for c, g in grp.items() if sel(c)]
        na = sum(g["n_ani"] for g in cs)
        nf = sum(g["n_af"] for g in cs)
        ca = sum(g["ani"] * g["n_ani"] for g in cs) / max(1, na)
        cf = sum(g["af"] * g["n_af"] for g in cs) / max(1, nf)
        rows.append([name, str(len(cs)), f"{na:,}", f"{ca * 100:.1f}", f"{cf * 100:.1f}"])
        nums[f"ci_simfam_{name.split()[0]}_ani"] = ca
        nums[f"ci_simfam_{name.split()[0]}_af"] = cf
    tbl("tab_ci_sim.tbl", ["l c c c c"], "l n n n n",
        [["simulated conditions (main)", "conditions", "contigs", "ANI coverage (%)", "AF coverage (%)"]], rows)
    # calibrated, per stratum (held-out real)
    st = CAL["calibrated"]["0.95"]["real_eval"]["strata"]
    IDL = ["\\(>= 99%", "97\\(en99%", "93\\(en97%", "< 93%"]
    LNL = ["< 1 kb", "1\\(en3 kb", "\\(>= 3 kb"]
    rows = []
    for k in sorted(st, key=lambda x: int(x)):
        v = st[k]
        if v["n_af"] < 100:
            continue
        ia, rest = divmod(int(k), 6)
        il, ic = divmod(rest, 2)
        rows.append([IDL[ia], LNL[il], "full" if ic == 0 else "partial", f"{v['n_af']:,}", f"{v['ani'] * 100:.1f}",
                     f"{v['af'] * 100:.1f}", f"{v['hw_ani']:.2f}", f"{v['hw_af']:.2f}"])
    tbl("tab_ci_strata.tbl", ["l l l c c c c c"], "l l l n n n n n",
        [["identity", "length", "AF", "contigs", "ANI cov. (%)", "AF cov. (%)", "ANI half-w. (pp)", "AF half-w. (pp)"]], rows)
    # calibration figure
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    xs = [90, 95, 99]
    series = [("calibrated, real held out", "calibrated", "real_eval", COL["aif"], "s", "-"),
              ("statistical, real held out", "stat", "real_eval", COL["skani"], "o", "--"),
              ("statistical, simulated main", "stat", "sim_main", COL["aif_rt"], "^", "-")]
    for ax, key, title in ((axes[0], "ani", "ANI intervals"), (axes[1], "af", "AF intervals")):
        ax.plot([88, 100], [88, 100], color=GRID, lw=1)
        ax.set_xlim(88, 100)
        for lab, mode, tag, col, mk, ls in series:
            ax.plot(xs, [CAL[mode][l][tag][key] * 100 for l in lvls], ls + mk, color=col, lw=1.6, ms=5, label=lab)
        ax.set_xticks(xs)
        ax.set_xlabel("nominal coverage (%)", fontsize=7.5)
        ax.set_title(title, fontsize=8, color=INK)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("empirical coverage (%)", fontsize=7.5)
    axes[0].set_ylim(60, 101)
    axes[1].set_ylim(60, 101)
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_calibration.eps"))
    plt.close(fig)

    # ================= v1.4: AF model grid and options =================
    grid = tsv(os.path.join(R, "results", "af_model_grid.tsv"))
    rows = []
    for g in grid:
        rows.append([g["setting"], f(float(g["real_ani_mae"]), 3), f(float(g["real_af_mae"]), 2), f(float(g["genome_ani_mae"]), 3),
                     f(float(g["sim_ani_mae"]), 3), f(float(g["sim_af_mae"]), 2)])
    tbl("tab_afgrid.tbl", ["l c s c c s", "^ c c c c c"], "l n n n n n",
        [["setting (tuning sets)", "real contigs", "", "real genomes", "simulated", ""],
         ["", "ANI MAE", "AF MAE", "ANI MAE", "ANI MAE", "AF MAE"]], rows)
    SIMO = json.load(open(os.path.join(R, "results", "afopt_sim.json")))
    rows = []
    for k, lab in (("aif", "default"), ("aif_rt", "aif-RT (--rho-u 2e-5)"), ("aif_lp", "local prior (0.08, R = 3)")):
        gk = {"aif": "aif_agg", "aif_rt": "aif_rt_agg", "aif_lp": "aif_lp_agg"}[k]
        sk_ = {"aif": "default", "aif_rt": "aif_rt", "aif_lp": "aif_lp"}[k]
        rows.append([lab, f(M[k]["ani_mae"], 3), f(M[k]["af_mae"], 2), f(G["all"][gk]["ani_mae"], 3), f(G["all"][gk]["af_mae"], 2),
                     f(SIMO["bench"][sk_]["ani_macro"], 3), f(SIMO["bench"][sk_]["af_macro"], 2)])
        nums[f"afopt_{k}"] = (M[k]["ani_mae"], M[k]["af_mae"], G["all"][gk]["ani_mae"], G["all"][gk]["af_mae"],
                              SIMO["bench"][sk_]["ani_macro"], SIMO["bench"][sk_]["af_macro"])
    tbl("tab_afopt.tbl", ["l c s c s c s", "^ c c c c c c"], "l n n n n n n",
        [["held-out evaluation", "real contigs (macro)", "", "real genomes (224)", "", "simulated main", ""],
         ["", "ANI MAE", "AF MAE", "ANI MAE", "AF MAE", "ANI MAE", "AF MAE"]], rows)
    # trade-off figure
    fam_style = [("default", INK, "*", True), ("switching rate", COL["skani"], "o", True), ("contig window", COL["aif"], "s", False),
                 ("local prior", COL["aif"], "s", True), ("overdispersion", COL["aif_rt"], "^", True)]
    short = {"default": "default", "switching rate 1e-5": "1e-5", "switching rate 2e-5 (aif-RT)": "2e-5 (RT)",
             "switching rate 3e-5": "3e-5", "contig window 0.06": "W6", "contig window 0.08": "W8", "contig window 0.10": "W10",
             "local prior 0.08, R=1": "L8/1", "local prior 0.08, R=3": "L8/3", "local prior 0.10, R=1": "L10/1",
             "local prior 0.10, R=3": "L10/3", "overdispersion 1": "od", "overdispersion 1, rate 1e-5": "od+1e-5"}
    right_labels = {"default", "switching rate 1e-5", "switching rate 2e-5 (aif-RT)", "switching rate 3e-5", "overdispersion 1",
                    "overdispersion 1, rate 1e-5", "contig window 0.06"}
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, ykey, ylab in ((axes[0], "genome_ani_mae", "genome-level ANI MAE (pp)"), (axes[1], "sim_af_mae", "simulated AF MAE (pp)")):
        for famname, col, mk, filled in fam_style:
            pts = [g for g in grid if g["setting"].startswith(famname)]
            if not pts:
                continue
            ax.plot([float(g["real_af_mae"]) for g in pts], [float(g[ykey]) for g in pts], mk, color=col,
                    mfc=col if filled else "white", mew=1.3, ms=7 if mk == "*" else 5.5, label=famname, linestyle="none")
            for g in pts:
                if ykey == "sim_af_mae" and g["setting"] not in right_labels:
                    continue
                ax.annotate(short.get(g["setting"], ""), (float(g["real_af_mae"]), float(g[ykey])), xytext=(3, 3),
                            textcoords="offset points", fontsize=5.8, color=INK2)
        ax.set_xlabel("real contig AF MAE (pp)", fontsize=7.5)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].legend(frameon=False, fontsize=6.3, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_af_tradeoff.eps"))
    plt.close(fig)

    db_section(nums)

    with open(os.path.join(GEN, "numbers.txt"), "w") as fh:
        for k, v in nums.items():
            fh.write(f"{k}\t{v}\n")
    print("\n".join(f"{k} = {v}" for k, v in nums.items()))


def db_section(nums):
    """Database mode (v1.5): tables 21-25 and figure 8, from data/db/results (bench_real/db_bench.py)."""
    DB = os.path.join(ROOT, "data", "db", "results")
    if not os.path.exists(os.path.join(DB, "db_eval.json")):
        print("database results missing; skipping the database section")
        return
    E = json.load(open(os.path.join(DB, "db_eval.json")))
    K = json.load(open(os.path.join(R, "results", "eval_k.json")))
    SG = json.load(open(os.path.join(DB, "screen_grid.json")))
    BLD, TIM = E.get("build", {}), E.get("time", {})
    # (21) cost of one shared k (pairwise, real sets, pooled within-species)
    rows = []
    for k, lab in (("aif", "adaptive k (17\\(en20), default"), ("aif_k19", "k = 19 for every reference"),
                   ("aif_k20", "k = 20 for every reference (database default)"), ("skani", "skani (default)")):
        w = K["sets"]["all_within"][k]
        rows.append([lab, f(w["coverage"], 3), f(w["ani_mae"], 3), f(w["paired_ani_mae"], 3), f(w["af_mae"], 2)])
        nums[f"k_{k}"] = (w["coverage"], w["ani_mae"], w["paired_ani_mae"], w["af_mae"])
    tbl("tab_db_k.tbl", ["l c c c c"], "l n n n n",
        [["pairwise, 33,464 within-species contigs", "coverage", "ANI MAE", "paired ANI", "AF MAE"]], rows)
    # (22) database build and size
    rows = []
    for name in ("db405", "db805", "db1205", "db2005"):
        b = BLD[name]
        a, s = b["aif"], b["skani"]
        rows.append([str(b["refs"]), f(b["bases"] / 1e9, 2), f(a["build_s"], 1), f(a["build_kb"] / 1e6, 2), f(a["bytes"] / 1e6, 0),
                     f(8 * a["bytes"] / b["bases"], 3), f(s["build_s"], 1), f(s["build_kb"] / 1e6, 2), f(s["bytes"] / 1e6, 0),
                     f(8 * s["bytes"] / b["bases"], 3)])
    b = BLD["db2005"]
    a = b["aif_b074"]
    rows.append(["2005*", f(b["bases"] / 1e9, 2), f(a["build_s"], 1), f(a["build_kb"] / 1e6, 2), f(a["bytes"] / 1e6, 0),
                 f(8 * a["bytes"] / b["bases"], 3), "", "", "", ""])
    nums["db2005_aif_bpb"] = 8 * b["aif"]["bytes"] / b["bases"]
    nums["db2005_skani_bpb"] = 8 * b["skani"]["bytes"] / b["bases"]
    nums["db2005_aif_screen_bpb"] = 8 * b["aif"]["screen_bytes"] / b["bases"]
    nums["db2005_build"] = (b["aif"]["build_s"], b["skani"]["build_s"])
    tbl("tab_db_build.tbl", ["c c c s s s c s s s", "^ ^ c c c c c c c c"], "n n n n n n n n n n",
        [["references", "Gbp", "aif db build", "", "", "", "skani sketch", "", "", ""],
         ["", "", "time (s)", "peak (GB)", "MB", "bits/base", "time (s)", "peak (GB)", "MB", "bits/base"]], rows)
    # (23) accuracy on db2005
    P = E["pooled"]
    labs = (("aif", "aif db search"), ("aif_b074", "aif, size-matched database"), ("aif_exhaustive", "aif, exhaustive (no screening)"),
            ("skani", "skani search"), ("skani_dist", "skani dist (all references)"))
    rows = []
    for k, lab in labs:
        r = P[k]
        rows.append([lab, f(r["coverage"], 3), f(r["ani_mae_pp"], 3), f(r["af_mae_pp"], 2), f(r["af_mae_all_pp"], 2),
                     f(r["regret_pp"], 2), f(r["top1_best"], 3), f(r["regret_ani_pp"], 3),
                     f(r["ci_cov_ani"], 3) if r.get("ci_cov_ani") is not None else "\\(em",
                     f(r["ci_cov_af"], 3) if r.get("ci_cov_af") is not None else "\\(em"])
        nums[f"db_{k}"] = (r["coverage"], r["ani_mae_pp"], r["af_mae_pp"], r["af_mae_all_pp"], r["regret_pp"], r["top1_best"],
                           r["regret_ani_pp"], r.get("ci_cov_ani"), r.get("ci_cov_af"), r["n_ani"])
    nums["db_paired"] = (P["paired"]["n"], P["paired"]["aif_ani_mae_pp"], P["paired"]["skani_ani_mae_pp"])
    nums["db_paired_dist"] = (P["paired_dist"]["n"], P["paired_dist"]["aif_ani_mae_pp"], P["paired_dist"]["skani_ani_mae_pp"])
    nums["db_screen_recall"] = P["screen_recall"]
    nums["db_n_hom"] = P["aif"]["n_hom"]
    tbl("tab_db_acc.tbl", ["l c c c c c s c c s", "^ ^ ^ ^ ^ c c c c c"], "l n n n n n n n n n",
        [["top hit (2,005 genomes)", "coverage", "ANI MAE", "AF MAE", "AF MAE*", "regret (ANI\\(muAF)", "", "ANI", "CI coverage", ""],
         ["", "", "", "", "", "pp", "best", "regret", "ANI", "AF"]], rows)
    rows = []
    for s in SETS:
        m = E["sets"][s]
        rows.append([SSHORT[s], str(m["aif"]["n_hom"]), f(m["aif"]["coverage"], 3), f(m["skani"]["coverage"], 3),
                     f(m["skani_dist"]["coverage"], 3), f(m["paired_dist"]["aif_ani_mae_pp"], 3), f(m["paired_dist"]["skani_ani_mae_pp"], 3),
                     f(m["aif"]["af_mae_all_pp"], 2), f(m["skani"]["af_mae_all_pp"], 2), f(m["skani_dist"]["af_mae_all_pp"], 2),
                     f(m["aif"]["regret_pp"], 2), f(m["skani_dist"]["regret_pp"], 2), f(m["screen_recall"], 3)])
    tbl("tab_db_sets.tbl", ["l c c s s c s c s s c s c", "^ ^ c c c c c c c c c c ^"], "l n n n n n n n n n n n n",
        [["set", "homol.", "coverage", "", "", "paired ANI MAE", "", "AF MAE*", "", "", "regret", "", "screen"],
         ["", "contigs", "aif", "search", "dist", "aif", "dist", "aif", "search", "dist", "aif", "dist", "recall"]], rows,
        opts="center box tab(@) expand")
    # (24) screening configurations
    rows = []
    for c in SG["configs"]:
        rows.append([c["frac"], str(c["bits"]), str(c["hashes"]), f(c["screen_bits_per_base"], 3), c["ranking"],
                     str(c["max_cand"]) if not c["prescreen"] else f"{c['prescreen']}\\(->{c['max_cand']}", f(c["recall"], 4),
                     f(c["reported"], 4), f(c.get("time_s_t8"), 1)])
    tbl("tab_db_screen.tbl", ["c c c c l c c c c"], "n n n n l n n n n",
        [["fraction", "bits", "hashes", "bits/base", "ranking", "evaluated", "recall", "any hit", "time (s)"]], rows)
    # (25) search time and memory
    rows = []
    for name in ("db405", "db805", "db1205", "db2005"):
        r = [str(BLD[name]["refs"])]
        for tool in ("aif", "skani", "skani_keep"):
            for th in (1, 8):
                t = TIM.get(f"{name}.{tool}.t{th}")
                r.append(f(t["s"], 1) if t else "\\(em")
            t = TIM.get(f"{name}.{tool}.t8")
            r.append(f(t["kb"] / 1e6, 2) if t else "\\(em")
        t = TIM.get(f"{name}.skani_dist.t8")
        r += [f(t["s"], 1) if t else "\\(em", f(t["kb"] / 1e6, 2) if t else "\\(em"]
        rows.append(r)
        for tool in ("aif", "skani", "skani_keep", "skani_dist"):
            for th in (1, 8):
                t = TIM.get(f"{name}.{tool}.t{th}")
                if t:
                    nums[f"time_{name}_{tool}_t{th}"] = (t["s"], t["kb"] / 1e6)
    tbl("tab_db_time.tbl", ["c c s s c s s c s s c s", "^ c c c c c c c c c c c"], "n n n n n n n n n n n n",
        [["refs", "aif db search", "", "", "skani search", "", "", "search --keep-refs", "", "", "skani dist", ""],
         ["", "1 thr", "8 thr", "GB", "1 thr", "8 thr", "GB", "1 thr", "8 thr", "GB", "8 thr", "GB"]], rows)
    # figure 8: search time and memory against database size
    names_ = ("db405", "db805", "db1205", "db2005")
    sizes = [BLD[n]["refs"] for n in names_]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7))
    styles = (("aif", "aif db search", COL["aif"], "o"), ("skani", "skani search", COL["skani"], "s"),
              ("skani_keep", "skani search --keep-refs", "#eda100", "^"), ("skani_dist", "skani dist", INK2, "D"))
    handles = []
    for tool, lab, col, mk in styles:
        for th, ls in ((8, "-"), (1, "--")):
            pts = [(BLD[n]["refs"], TIM[f"{n}.{tool}.t{th}"]["s"]) for n in names_ if f"{n}.{tool}.t{th}" in TIM]
            if not pts:
                continue
            h, = axes[0].plot([x for x, _ in pts], [y for _, y in pts], ls + mk, color=col, lw=1.5, ms=4,
                              mfc=col if th == 8 else "white")
            if th == 8:
                handles.append((h, lab))
        pts = [(BLD[n]["refs"], TIM[f"{n}.{tool}.t8"]["kb"] / 1e6) for n in names_ if f"{n}.{tool}.t8" in TIM]
        axes[1].plot([x for x, _ in pts], [y for _, y in pts], "-" + mk, color=col, lw=1.5, ms=4)
    axes[0].set_yscale("log")
    yt = [3, 5, 10, 20, 50, 100, 200]
    axes[0].set_yticks(yt)
    axes[0].set_yticklabels([str(y) for y in yt])
    axes[0].minorticks_off()
    axes[0].set_ylabel("wall time (s), 33,464 contigs", fontsize=7.5)
    axes[1].set_ylabel("peak memory (GB), 8 threads", fontsize=7.5)
    for ax in axes:
        ax.set_xlabel("references in the database", fontsize=7.5)
        ax.set_xticks(sizes)
        ax.grid(True, which="both", color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.legend([h for h, _ in handles], [l for _, l in handles], loc="upper center", ncol=4, frameon=False, fontsize=6.8)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(FIG, "fig_db_scaling.eps"))
    plt.close(fig)
    # committed summaries
    out = os.path.join(ROOT, "results_real", "db")
    os.makedirs(out, exist_ok=True)
    for src, dst in ((os.path.join(DB, "db_eval.json"), "db_eval.json"), (os.path.join(DB, "screen_grid.json"), "screen_grid.json"),
                     (os.path.join(R, "results", "eval_k.json"), "eval_fixed_k.json")):
        open(os.path.join(out, dst), "w").write(open(src).read())


if __name__ == "__main__":
    main()
