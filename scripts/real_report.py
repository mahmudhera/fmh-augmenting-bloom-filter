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
    pareto_section(nums)
    idmap_section(nums)
    gtdb_section(nums)
    fresh_section(nums)
    tools_section(nums)

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


def pareto_section(nums):
    """Pareto benchmark: error against sketch size and time (table 27, figure 9)."""
    E = json.load(open(os.path.join(R, "results", "eval_pareto.json")))
    SZ = json.load(open(os.path.join(R, "results", "pareto_sizes.json")))["bits_per_base"]
    PT = json.load(open(os.path.join(R, "results", "pareto_time.json")))
    SMT = json.load(open(os.path.join(R, "results", "pareto_sourmash.json")))
    MMT = json.load(open(os.path.join(R, "results", "pareto_minimap2.json")))
    W = E["sets"]["all_within"]
    tkey = {"aif_b011": "aif_0.11", "aif_b022": "aif_0.22", "aif_0.44": "aif_0.44", "aif": "aif_0.88", "aif_b176": "aif_1.76",
            "skani_c500": "skani_500", "skani_c250": "skani_250", "skani": "skani_125", "skani_c60": "skani_60", "skani_slow": "skani_30"}

    def tsec(k):
        if k in tkey:
            return PT[tkey[k]]["median_s"]
        if k.startswith("sm_S"):
            return SMT[k[4:]]["time_s"]
        if k == "minimap2":
            return MMT["time_s"]
        return None
    fam = [("aif", ["aif_b011", "aif_b022", "aif_0.44", "aif", "aif_b176"], lambda k: "budget " + {"aif_b011": "0.11", "aif_b022": "0.22",
                                                                                                    "aif_0.44": "0.44", "aif": "0.88", "aif_b176": "1.76"}[k]),
           ("skani", ["skani_c500", "skani_c250", "skani", "skani_c60", "skani_slow"], lambda k: "c = " + {"skani_c500": "500", "skani_c250": "250",
                                                                                                          "skani": "125", "skani_c60": "60", "skani_slow": "30"}[k]),
           ("FracMinHash (sourmash)", ["sm_S580", "sm_S290", "sm_S145", "sm_S73", "sm_S36", "sm_S18"], lambda k: "scaled " + k[4:]),
           ("minimap2 (alignment)", ["minimap2"], lambda k: "asm20")]
    rows = []
    for fname, ks, lab in fam:
        for k in ks:
            w = W[k]
            t = tsec(k)
            rows.append([fname if k == ks[0] else "", lab(k), f(SZ[k], 3 if SZ[k] < 10 else 1), f(w["coverage"], 3), f(w["ani_mae"], 3),
                         f(w["paired_ani_mae"], 3), "\\(em" if k.startswith("sm_") else f(w["af_mae"], 2), f(t, 1) if t is not None else "\\(em"])
            nums[f"pareto_{k}"] = (SZ[k], w["coverage"], w["ani_mae"], w["paired_ani_mae"], w["af_mae"], t)
        rows.append("_")
    rows = rows[:-1]
    tbl("tab_pareto.tbl", ["l l c c c c c c"], "l l n n n n n n",
        [["method", "setting", "bits/base", "coverage", "ANI MAE", "paired ANI", "AF MAE", "time (s)"]], rows)
    # figure 9
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.5))
    sty = [("aif", ["aif_b011", "aif_b022", "aif_0.44", "aif", "aif_b176"], COL["aif"], "o"),
           ("skani", ["skani_c500", "skani_c250", "skani", "skani_c60", "skani_slow"], COL["skani"], "s"),
           ("FracMinHash (sourmash)", ["sm_S580", "sm_S290", "sm_S145", "sm_S73", "sm_S36", "sm_S18"], INK2, "^"),
           ("minimap2", ["minimap2"], "#1baf7a", "D")]
    hs = []
    for lab, ks, col, mk in sty:
        x = [SZ[k] for k in ks]
        h, = axes[0].plot(x, [W[k]["ani_mae"] for k in ks], "-" + mk, color=col, lw=1.4, ms=4)
        hs.append((h, lab))
        if not ks[0].startswith("sm_"):
            axes[1].plot(x, [W[k]["af_mae"] for k in ks], "-" + mk, color=col, lw=1.4, ms=4)
        tt = [tsec(k) for k in ks]
        axes[2].plot(tt, [W[k]["ani_mae"] for k in ks], "-" + mk, color=col, lw=1.4, ms=4)
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_xlabel("sketch / index size (bits per reference base)", fontsize=7)
    axes[1].set_xlabel("sketch / index size (bits per reference base)", fontsize=7)
    axes[2].set_xlabel("sketch + query time, 1 thread (s)", fontsize=7)
    axes[0].set_ylabel("contig ANI MAE (pp)", fontsize=7.5)
    axes[1].set_ylabel("contig AF MAE (pp)", fontsize=7.5)
    axes[2].set_ylabel("contig ANI MAE (pp)", fontsize=7.5)
    axes[1].set_yscale("log")
    axes[1].set_yticks([1, 2, 3, 5, 10, 20, 30])
    axes[1].set_yticklabels(["1", "2", "3", "5", "10", "20", "30"])
    axes[1].minorticks_off()
    fig.legend([h for h, _ in hs], [l for _, l in hs], loc="upper center", ncol=4, frameon=False, fontsize=6.8)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(FIG, "fig_pareto.eps"))
    plt.close(fig)


def idmap_section(nums):
    """Identity maps along contigs (tables 28-29, figure 10)."""
    D = json.load(open(os.path.join(R, "results", "idmap.json")))
    classes = ["<90", "90-95", "95-98", "98-99.5", ">=99.5"]
    rows = []
    for key, lab in (("aif_1000", "aif, 1 kb windows"), ("aif_3000", "aif, 3 kb windows"), ("fastani_3000", "fastANI, 3 kb fragments")):
        d = D[key]
        cal = D.get(key + "_cal", {}).get("all", {})
        rows.append([lab, str(d["all"]["n"]), f(d["all"]["mae_pp"], 3), f(d["all"]["r"], 3)] +
                    [f(d["all|" + c]["mae_pp"], 2) for c in classes] +
                    [f(d["all"]["ci_cov"], 3) if d["all"]["ci_cov"] is not None else "\\(em", f(cal.get("ci_cov"), 3) if cal else "\\(em"])
        nums[f"idmap_{key}"] = (d["all"]["n"], d["all"]["mae_pp"], d["all"]["r"], d["all"]["ci_cov"], cal.get("ci_cov"))
    tbl("tab_idmap.tbl", ["l c c c c s s s s c s", "^ ^ ^ ^ c c c c c c c"], "l n n n n n n n n n n",
        [["windows with true AF \\(>= 0.9", "n", "MAE", "r", "MAE by true identity (%)", "", "", "", "", "95% CI coverage", ""],
         ["", "", "", "", "<90", "90\\(en95", "95\\(en98", "98\\(en99.5", "\\(>=99.5", "stat.", "calib."]], rows)
    rows = []
    for st in SETS + ["all"]:
        a_, f_ = D["aif"]["islands|" + st], D["fastani"]["islands|" + st]
        b_ = D["aif"]["bp|" + st]
        rows.append([SSHORT.get(st, "all sets"), f(a_["true_bp"] / 1e6, 1), f(a_["precision"], 3), f(a_["recall"], 3), f(f_["precision"], 3),
                     f(f_["recall"], 3), str(b_["n"]), str(b_["median_bp"]), f(b_["within_250"], 3), f(b_["within_500"], 3)])
    nums["idmap_islands"] = (D["aif"]["islands|all"], D["fastani"]["islands|all"], D["aif"]["bp|all"])
    tbl("tab_islands.tbl", ["l c c s c s c c c c", "^ ^ c c c c c c c c"], "l n n n n n n n n n",
        [["set", "unaligned", "aif", "", "fastANI", "", "break-", "median", "within", "within"],
         ["", "Mbp", "precision", "recall", "precision", "recall", "points", "dist. (bp)", "250 bp", "500 bp"]], rows)
    # figure 10: identity map of one mostly aligned contig with identity variation and an island
    Q = json.load(open(os.path.join(R, "results", "idmap_ci.json")))["1000"]

    def load_set(st):
        seg = {}
        for line in list(open(os.path.join(R, "idmap", st, "aif_seg1000.tsv")))[1:]:
            p = line.rstrip("\n").split("\t")
            seg.setdefault(p[0], []).append((int(p[1]), int(p[2]), int(p[3]), float(p[4]), float(p[5]), float(p[6])))
        return seg
    best, bs, bset = None, -1, None
    for st in ("ecoli", "paer"):
        seg = load_set(st)
        nfa = {}
        for fn in os.listdir(os.path.join(R, "idmap", st)):
            if fn.endswith(".fastani.tsv"):
                for l in open(os.path.join(R, "idmap", st, fn)):
                    nfa[l.split("\t")[0]] = nfa.get(l.split("\t")[0], 0) + 1
        for c, v in seg.items():
            L = v[-1][1]
            if L < 60000 or L > 250000 or nfa.get(c, 0) < 12:
                continue
            al = sum(x[1] - x[0] for x in v if x[2] == 1) / L
            isl = sum(1 for x in v if x[2] == 0 and x[1] - x[0] >= 2000)
            ids = sorted(x[3] for x in v if x[2] == 1 and x[3] == x[3])
            if al < 0.75 or isl < 1 or len(ids) < 20 or ids[len(ids) // 2] < 0.95:
                continue
            spread = ids[int(0.9 * len(ids))] - ids[int(0.1 * len(ids))]
            sc = 100 * spread + min(isl, 3)
            if sc > bs:
                best, bs, bset = c, sc, st
    c, st = best, bset
    seg = load_set(st)[c]
    gid = c.split("__")[0]
    blocks = [(int(x[1]) - 1, int(x[2]), float(x[3]) / 100) for x in (l.split("\t") for l in open(os.path.join(R, "idmap", st, gid + ".coords.tsv")))
              if x[0] == c]
    fa = [(int(x[1]), int(x[2]), float(x[3]) / 100) for x in (l.split("\t") for l in open(os.path.join(R, "idmap", st, gid + ".fastani.tsv"))) if x[0] == c]
    fig, ax = plt.subplots(figsize=(6.8, 2.5))
    for a0, b0, al, ani, lo, hi in seg:
        if not al:
            ax.axvspan(a0 / 1e3, b0 / 1e3, color="#e8e7e3", lw=0)
    xs, ys, lo_, hi_, tx, ty = [], [], [], [], [], []
    for a0, b0, al, ani, lo, hi in seg:
        w = sum(max(0, min(b0, y) - max(a0, x)) for x, y, _ in blocks)
        tid = sum(max(0, min(b0, y) - max(a0, x)) * i for x, y, i in blocks) / w if w else float("nan")
        if not al:  # break the lines across unaligned segments
            for v in (xs, ys, lo_, hi_, tx, ty):
                v.append(float("nan"))
        if al and ani == ani:
            k = "ge99" if ani >= 0.99 else "97-99" if ani >= 0.97 else "93-97" if ani >= 0.93 else "lt93"
            se = max((hi - lo) / (2 * 1.959964), 1e-4)
            xs.append((a0 + b0) / 2e3)
            ys.append(100 * ani)
            lo_.append(100 * max(0, ani - Q[k]["q95"] * se))
            hi_.append(100 * min(1, ani + Q[k]["q95"] * se))
            tx.append((a0 + b0) / 2e3)
            ty.append(100 * tid)
    ax.fill_between(xs, lo_, hi_, color="#c9def5", lw=0, label="aif calibrated 95% interval")
    ax.plot(tx, ty, "-", color=INK, lw=1.0, label="MUMmer identity (same windows)")
    ax.plot(xs, ys, "-", color=COL["aif"], lw=1.2, label="aif, 1 kb windows")
    for a0, b0, idy in fa:
        ax.plot([a0 / 1e3, b0 / 1e3], [100 * idy] * 2, color=COL["skani"], lw=2.0, solid_capstyle="butt")
    ax.plot([], [], color=COL["skani"], lw=2.0, label="fastANI, 3 kb fragments")
    ax.fill_between([], [], color="#e8e7e3", label="aif: unaligned")
    ax.set_xlabel(f"position on contig (kb); {SSHORT[st]} draft {gid} vs its reference", fontsize=7.5)
    ax.set_ylabel("identity (%)", fontsize=7.5)
    allv = [v for v in ty + ys + [100 * x[2] for x in fa] if v == v]
    ax.set_ylim(max(80, min(allv) - 1), 100.3)
    ax.set_xlim(0, seg[-1][1] / 1e3)
    ax.grid(True, color=GRID, lw=0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=6.3, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_idmap.eps"))
    plt.close(fig)
    nums["idmap_example"] = c


def gtdb_section(nums):
    """GTDB-scale database and contig classification (tables 27-28, 32; figure 9)."""
    GR = os.path.join(ROOT, "data", "gtdb", "results")
    if not os.path.exists(os.path.join(GR, "gtdb_eval.json")):
        print("GTDB results missing; skipping")
        return
    E = json.load(open(os.path.join(GR, "gtdb_eval.json")))
    B, S = E.get("build", {}), E.get("search_time", {})
    LAB = {"aif_0.88": "aif (budget 0.88)", "aif_0.22": "aif (budget 0.22)", "skani": "skani 0.4.0"}
    rows = []
    for k in ("aif_0.88", "aif_0.22", "skani"):
        b = B.get(k)
        st = S.get(k)
        if not b:
            continue
        bases = B["aif_0.88"]["bases"]
        rows.append([LAB[k], f(b["build_s"] / 60, 1), f(b["build_kb"] / 1e6, 1), f(b["bytes"] / 1e9, 1), f(8 * b["bytes"] / bases, 3),
                     f(st["s"] / 60, 1) if st else "\\(em", f(st["kb"] / 1e6, 1) if st else "\\(em"])
        nums[f"gtdb_build_{k}"] = (b["build_s"], b["build_kb"], b["bytes"], 8 * b["bytes"] / bases, st["s"] if st else None, st["kb"] if st else None)
    nums["gtdb_refs"] = (B["aif_0.88"]["refs"], B["aif_0.88"]["bases"])
    tbl("tab_gtdb_build.tbl", ["l c c c c c c", "^ c c c c c c"], "l n n n n n n",
        [["GTDB R232 representatives", "build", "build peak", "database", "bits per", "search", "search peak"],
         ["", "(min)", "(GB)", "(GB)", "base", "(min)", "(GB)"]], rows)
    rows = []
    for k in ("aif_0.88", "aif_0.22", "skani"):
        if k not in E["tools"]:
            continue
        a = E["tools"][k]["all"]
        rows.append([LAB[k], f(a["reported"], 3), f(a["species_acc"], 3), f(a["genus_acc"], 3), f(a["species_prec"], 3),
                     f(a["species_acc_homol"], 3), f(a["call95_sens"], 3), f(a["call95_prec"], 3), f(a["ani_mae_pp"], 3), f(a["af_mae_pp"], 2)])
        nums[f"gtdb_{k}"] = a
        nums[f"gtdb_{k}_ar"] = E["tools"][k]["ar"]
    tbl("tab_gtdb_acc.tbl", ["l c c s c c c s c c", "^ ^ c c c c c c c c"], "l n n n n n n n n n",
        [["tool", "any hit", "top hit correct", "", "species", "species,", "ANI \\(>= 95% call", "", "ANI MAE", "AF MAE"],
         ["", "", "species", "genus", "precision", "homol.", "sens.", "prec.", "(pp)", "(pp)"]], rows)
    # figure 9: species-level accuracy and ANI error by contig length
    bins = ["0.5-1k", "1-2k", "2-5k", "5-20k", ">=20k"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4))
    sty = {"aif_0.88": (COL["aif"], "o", "-"), "aif_0.22": (COL["aif"], "o", "--"), "skani": (COL["skani"], "s", "-")}
    for k, (col, mk, ls) in sty.items():
        if k not in E["tools"]:
            continue
        d = E["tools"][k]
        axes[0].plot(range(len(bins)), [d[b]["species_acc"] for b in bins], ls + mk, color=col, lw=1.5, ms=4, mfc=col if ls == "-" else "white",
                     label=LAB[k])
        axes[1].plot(range(len(bins)), [d[b]["ani_mae_pp"] for b in bins], ls + mk, color=col, lw=1.5, ms=4, mfc=col if ls == "-" else "white")
    for ax in axes:
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels(["0.5\u20131", "1\u20132", "2\u20135", "5\u201320", "\u226520"])
        ax.set_xlabel("contig length (kb)", fontsize=7.5)
        ax.grid(True, color=GRID, lw=0.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("top hit = own species", fontsize=7.5)
    axes[1].set_ylabel("top-hit ANI MAE (pp)", fontsize=7.5)
    axes[0].legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_gtdb.eps"))
    plt.close(fig)
    # contamination (table 32)
    cp = os.path.join(GR, "contam.json")
    if os.path.exists(cp):
        C = json.load(open(cp))
        rows = []
        for tool in ("aif", "skani"):
            for rank in ("species", "genus"):
                for rule in ("default", "strict"):
                    r = [("aif" if tool == "aif" else "skani") if rank == "species" and rule == "default" else "", rank, rule]
                    for lv in ("0.05", "0.1", "0.2", "all"):
                        v = C[f"{tool}|{rank}|{rule}|{lv}"]
                        r += [f(v["precision_bases"], 3), f(v["recall_bases"], 3)]
                    rows.append(r)
                    nums[f"contam_{tool}_{rank}_{rule}"] = C[f"{tool}|{rank}|{rule}|all"]
            rows.append("_")
        rows = rows[:-1]
        tbl("tab_contam.tbl", ["l l l c s c s c s c s", "^ ^ ^ c c c c c c c c"], "l l l n n n n n n n n",
            [["tool", "taxon", "rule", "5% contamination", "", "10%", "", "20%", "", "all bins", ""],
             ["", "", "", "prec.", "recall", "prec.", "recall", "prec.", "recall", "prec.", "recall"]], rows)
    out = os.path.join(ROOT, "results_real", "gtdb")
    os.makedirs(out, exist_ok=True)
    for fn in ("gtdb_eval.json", "contam.json", "build.json", "search_time.json"):
        if os.path.exists(os.path.join(GR, fn)):
            open(os.path.join(out, fn), "w").write(open(os.path.join(GR, fn)).read())


def fresh_section(nums):
    """Validation of the frozen v1.6 defaults on fresh query sets (tables 33-34)."""
    DBF = os.path.join(ROOT, "data", "db", "results_fresh", "db_eval.json")
    G2 = os.path.join(ROOT, "data", "gtdb", "results2", "gtdb_eval.json")
    if not (os.path.exists(DBF) and os.path.exists(G2)):
        print("fresh validation results missing; skipping")
        return
    E0 = json.load(open(os.path.join(ROOT, "data", "db", "results", "db_eval.json")))["pooled"]
    E1 = json.load(open(DBF))["pooled"]
    rows = []
    for k, lab in (("aif", "aif db search"), ("aif_b074", "aif, size-matched"), ("aif_exhaustive", "aif, exhaustive"),
                   ("skani", "skani search"), ("skani_dist", "skani dist")):
        r0, r1 = E0[k], E1[k]
        rows.append([lab, f(r0["coverage"], 3), f(r1["coverage"], 3), f(r0["ani_mae_pp"], 3), f(r1["ani_mae_pp"], 3), f(r0["af_mae_pp"], 2),
                     f(r1["af_mae_pp"], 2), f(r0["regret_ani_pp"], 3), f(r1["regret_ani_pp"], 3), f(r0["top1_best"], 3), f(r1["top1_best"], 3)])
        nums[f"fresh_db_{k}"] = (r1["coverage"], r1["ani_mae_pp"], r1["af_mae_pp"], r1["regret_ani_pp"], r1["top1_best"], r1["regret_pp"],
                                 r1.get("ci_cov_ani"), r1.get("ci_cov_af"), r1["n_hom"])
    nums["fresh_db_paired"] = (E1["paired"], E1["paired_dist"], E1["screen_recall"])
    tbl("tab_fresh_db.tbl", ["l c s c s c s c s c s", "^ c c c c c c c c c c"], "l n n n n n n n n n n",
        [["2,005-genome database", "coverage", "", "ANI MAE", "", "AF MAE", "", "ANI regret", "", "top hit best", ""],
         ["", "orig.", "fresh", "orig.", "fresh", "orig.", "fresh", "orig.", "fresh", "orig.", "fresh"]], rows)
    G0 = json.load(open(os.path.join(ROOT, "data", "gtdb", "results", "gtdb_eval.json")))["tools"]
    G1 = json.load(open(G2))["tools"]
    rows = []
    LAB = {"aif_0.88": "aif (budget 0.88)", "aif_0.22": "aif (budget 0.22)", "skani": "skani"}
    for k in ("aif_0.88", "aif_0.22", "skani"):
        a0, a1 = G0[k]["all"], G1[k]["all"]
        pair = lambda key, d=3: f"{f(a0[key], d)} / {f(a1[key], d)}"  # noqa: E731
        rows.append([LAB[k], pair("species_acc"), pair("genus_acc"), pair("call95_sens"), pair("call95_prec"), pair("ani_mae_pp"),
                     pair("af_mae_pp", 2)])
        nums[f"fresh_gtdb_{k}"] = a1
        nums[f"fresh_gtdb_{k}_ar"] = G1[k]["ar"]
    tbl("tab_fresh_gtdb.tbl", ["l c c c c c c", "^ c c c c c c"], "l c c c c c c",
        [["GTDB R232", "species", "genus", "\\(>= 95% call", "\\(>= 95% call", "ANI MAE", "AF MAE"],
         ["(set 1 / set 2)", "correct", "correct", "sensitivity", "precision", "(pp)", "(pp)"]], rows)
    for k in ("aif_0.88", "skani"):
        bins = ["0.5-1k", "1-2k", "2-5k", "5-20k", ">=20k"]
        nums[f"fresh_gtdb_bylen_{k}"] = [round(G1[k][b]["species_acc"], 3) for b in bins]
    cp = os.path.join(ROOT, "data", "gtdb", "results2", "contam.json")
    if os.path.exists(cp):
        C = json.load(open(cp))
        nums["fresh_contam"] = {t: C[f"{t}|genus|default|all"] for t in ("aif", "skani")}
    S2 = os.path.join(ROOT, "data", "gtdb", "results2", "search_time.json")
    if os.path.exists(S2):
        nums["fresh_gtdb_time"] = json.load(open(S2))
    out = os.path.join(ROOT, "results_real", "fresh")
    os.makedirs(out, exist_ok=True)
    for src, dst in ((DBF, "db_eval_fresh.json"), (G2, "gtdb_eval_set2.json"), (cp, "contam_set2.json"), (S2, "gtdb_search_time_set2.json")):
        if os.path.exists(src):
            open(os.path.join(out, dst), "w").write(open(src).read())
    for fn in ("queries2.tsv", "query_contigs2.tsv"):
        open(os.path.join(out, "gtdb_" + fn), "w").write(open(os.path.join(ROOT, "data", "gtdb", "meta", fn)).read())
    for s_ in ("bcereus", "ecoli", "paer", "saur"):
        open(os.path.join(out, f"real_fresh_{s_}.tsv"), "w").write(open(os.path.join(ROOT, "data", "real_fresh", "lists", s_ + ".tsv")).read())


def tools_section(nums):
    """Kraken2 on GTDB R232 (tables 35-36, figure 12); CheckM2 and GUNC on simulated bins (table 37, figure 13)."""
    GD = os.path.join(ROOT, "data", "gtdb")
    out = os.path.join(ROOT, "results_real", "tools")
    os.makedirs(out, exist_ok=True)
    COLX = {"kraken2": "#008300", "checkm2": "#e87ba4", "gunc": "#6250d6"}  # palette slots 6, 5, 7
    KP = os.path.join(GD, "kraken2_eval.json")
    G = [json.load(open(os.path.join(GD, d, "gtdb_eval.json"))) for d in ("results", "results2")]
    if os.path.exists(KP):
        K = json.load(open(KP))
        open(os.path.join(out, "kraken2_eval.json"), "w").write(open(KP).read())
        B, S1, S2 = G[0]["build"], G[0]["search_time"], json.load(open(os.path.join(GD, "results2", "search_time.json")))
        kb, kt = K["build"], K.get("classify_time", {})
        rows = []
        for k, lab in (("aif_0.88", "aif (budget 0.88)"), ("skani", "skani 0.4.0")):
            b = B[k]
            rows.append([lab, f(b["bytes"] / 1e9, 1), f(b["build_s"] / 60, 1), f(b["build_kb"] / 1e6, 1), "32", f(S1[k]["s"], 1),
                         f(S2[k]["s"], 1), f(max(S1[k]["kb"], S2[k]["kb"]) / 1e6, 1)])
        rows.append(["Kraken2 2.1.3", f(kb["bytes"] / 1e9, 1), f(kb["build_s"] / 60, 1), f(kb["build_kb"] / 1e6, 1), str(kb["threads"]),
                     f(kt["set1"]["s"], 1) if "set1" in kt else "\\(em", f(kt["set2"]["s"], 1) if "set2" in kt else "\\(em",
                     f(max(v["kb"] for v in kt.values()) / 1e6, 1) if kt else "\\(em"])
        nums["k2_build"] = kb
        nums["k2_time"] = kt
        tbl("tab_k2_build.tbl", ["l c c c c c c c", "^ c c c c c c c"], "l n n n n n n n",
            [["GTDB R232 (199,923)", "database", "build", "build peak", "build", "search (s)", "search (s)", "search peak"],
             ["", "(GB)", "(min)", "(GB)", "threads", "set 1", "set 2", "(GB)"]], rows)
        rows = []
        pr = lambda a, b, d=3: f"{f(a, d)} / {f(b, d)}"  # noqa: E731
        for k, lab in (("aif_0.88", "aif"), ("skani", "skani")):
            a, b = G[0]["tools"][k], G[1]["tools"][k]
            rows.append([lab, pr(a["all"]["reported"], b["all"]["reported"]), pr(a["all"]["species_acc"], b["all"]["species_acc"]),
                         pr(a["all"]["genus_acc"], b["all"]["genus_acc"]), pr(a["all"]["call95_sens"], b["all"]["call95_sens"]),
                         pr(a["all"]["call95_prec"], b["all"]["call95_prec"]), pr(a["0.5-1k"]["species_acc"], b["0.5-1k"]["species_acc"]),
                         pr(a[">=20k"]["species_acc"], b[">=20k"]["species_acc"])])
        a, b = K["set1"], K["set2"]
        rows.append(["Kraken2", pr(a["all"]["classified"], b["all"]["classified"]), pr(a["all"]["species_acc"], b["all"]["species_acc"]),
                     pr(a["all"]["genus_acc"], b["all"]["genus_acc"]), pr(a["all"]["species_acc"], b["all"]["species_acc"]),
                     pr(a["all"]["species_prec"], b["all"]["species_prec"]), pr(a["0.5-1k"]["species_acc"], b["0.5-1k"]["species_acc"]),
                     pr(a[">=20k"]["species_acc"], b[">=20k"]["species_acc"])])
        nums["k2_set1"], nums["k2_set2"] = K["set1"], K["set2"]
        tbl("tab_k2_acc.tbl", ["l c c c c c c c", "^ c c c c c c c"], "l c c c c c c c",
            [["set 1 / set 2", "assigned", "species", "genus", "species call", "species call", "species,", "species,"],
             ["", "", "correct", "correct", "sensitivity", "precision", "0.5\\(en1 kb", "\\(>= 20 kb"]], rows)
        # figure 12: species accuracy by contig length, aif / skani / Kraken2, both query sets
        bins = ["0.5-1k", "1-2k", "2-5k", "5-20k", ">=20k"]
        fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4), sharey=True)
        for ax, gi, key, title in ((axes[0], 0, "set1", "query set 1 (294 genomes)"), (axes[1], 1, "set2", "query set 2 (298 genomes)")):
            x = list(range(len(bins)))
            for k, col, mk, lab in (("aif_0.88", COL["aif"], "o", "aif (top hit)"), ("skani", COL["skani"], "s", "skani (top hit)")):
                ax.plot(x, [100 * G[gi]["tools"][k][b]["species_acc"] for b in bins], color=col, marker=mk, ms=4, lw=1.4, label=lab)
            ax.plot(x, [100 * K[key][b]["species_acc"] for b in bins], color=COLX["kraken2"], marker="^", ms=4, lw=1.4, label="Kraken2")
            ax.set_xticks(x)
            ax.set_xticklabels(["0.5\u20131", "1\u20132", "2\u20135", "5\u201320", "\u226520"])
            ax.set_xlabel("contig length (kb)")
            ax.set_title(title, fontsize=8, color=INK)
            ax.grid(axis="y", color=GRID, lw=0.5)
            ax.set_ylim(0, 100)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        axes[0].set_ylabel("species correct (%)")
        axes[0].legend(frameon=False, fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig_k2.eps"))
        plt.close(fig)
    else:
        print("Kraken2 results missing; skipping tables 35-36")
    CP = os.path.join(GD, "results", "contam_tools.json")
    BP = os.path.join(GD, "results", "contam_tools_bins.tsv")
    if not os.path.exists(CP):
        print("bin-tool results missing; skipping table 37")
        return
    C = json.load(open(CP))
    for src, dst in ((CP, "contam_tools.json"), (BP, "contam_tools_bins.tsv")):
        if os.path.exists(src):
            open(os.path.join(out, dst), "w").write(open(src).read())
    T = C.get("seq_time", {})
    TIM = {"aif": T.get("aif"), "skani": T.get("skani"), "checkm2": C.get("checkm2_time"), "gunc": C.get("gunc_time")}
    rows = []
    LABT = {"aif": "aif", "skani": "skani", "checkm2": "CheckM2 1.1.0", "gunc": "GUNC 1.1.1"}
    for k in ("aif", "skani", "checkm2", "gunc"):
        if k not in C:
            continue
        d = C[k]
        tm = TIM.get(k)
        rows.append([LABT[k], f(d["0.0"]["mae_pp"], 2), f(d["0.0"]["flag_5pct"], 2)] + [f(d[lv]["mae_pp"], 2) for lv in ("0.05", "0.1", "0.2")] +
                    [f(d[lv]["flag_5pct"], 2) for lv in ("0.05", "0.1", "0.2")] + [f(tm["s"], 0) if tm else "\\(em", f(tm["kb"] / 1e6, 1) if tm else "\\(em"])
        nums[f"bintool_{k}"] = d
        if k == "gunc":
            rows.append(["GUNC, own call", "\\(em", f(d["0.0"]["gunc_fail"], 2), "\\(em", "\\(em", "\\(em"] +
                        [f(d[lv]["gunc_fail"], 2) for lv in ("0.05", "0.1", "0.2")] + ["", ""])
    if os.path.exists(BP):  # detection at a lower threshold, near GUNC's own false-call rate on controls
        rb = tsv(BP)
        for k in ("aif", "skani"):
            fl = lambda lv: sum(float(r[k]) >= 0.03 for r in rb if r["level"] == lv) / max(1, sum(r["level"] == lv for r in rb))  # noqa: E731
            d3 = {lv: fl(lv) for lv in ("0.0", "0.05", "0.1", "0.2")}
            nums[f"bintool_{k}_3pct"] = d3
            rows.insert([r[0] for r in rows].index(LABT["skani"]) + 1 if k == "skani" else 1,
                        [f"{LABT[k]}, flag at 3%", "\\(em", f(d3["0.0"], 2), "\\(em", "\\(em", "\\(em"] + [f(d3[lv], 2) for lv in ("0.05", "0.1", "0.2")] + ["", ""])
    nums["bintool_time"] = TIM
    tbl("tab_bintools.tbl", ["l c s c s s c s s c c", "^ c c c c c c c c c c"], "l n n n n n n n n n n",
        [["tool", "controls", "", "MAE (pp) at", "", "", "flagged at", "", "", "time", "peak"],
         ["", "MAE (pp)", "flagged", "5%", "10%", "20%", "5%", "10%", "20%", "(s)", "(GB)"]], rows)
    if not os.path.exists(BP):
        return
    rs = tsv(BP)
    tools = [t for t in ("aif", "skani", "checkm2", "gunc") if t in rs[0]]
    fig, axes = plt.subplots(1, len(tools), figsize=(6.6, 1.9), sharey=True)
    colt = {"aif": COL["aif"], "skani": COL["skani"], "checkm2": COLX["checkm2"], "gunc": COLX["gunc"]}
    for ax, t in zip(axes, tools):
        xs = [100 * float(r["true"]) for r in rs if r[t] != "NA"]
        ys = [100 * float(r[t]) for r in rs if r[t] != "NA"]
        ax.plot([0, 30], [0, 30], color=INK2, lw=0.6, ls="--")
        ax.axhline(5, color=GRID, lw=0.6)
        ax.scatter(xs, ys, s=4, color=colt[t], lw=0)
        ax.set_title(LABT[t], fontsize=8, color=INK)
        ax.set_xlim(-1, 30)
        ax.set_ylim(-1, 40)
        ax.set_xlabel("true (%)")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("estimated contamination (%)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_bintools.eps"))
    plt.close(fig)


if __name__ == "__main__":
    main()
