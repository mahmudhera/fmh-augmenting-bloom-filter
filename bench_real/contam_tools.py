#!/usr/bin/env python3
"""Bin-level contamination estimates: aif and skani (contig classification against GTDB R232) vs CheckM2 and GUNC.

Bins are those of gtdb_contam.py (host genome + 5/10/20% of another query genome's contigs, seed 7). The
gene-based tools are expensive, so a deterministic subset is used: the first 100 bins of each level (hosts in
sorted order) and 100 uncontaminated controls (the first 100 hosts alone).
  true contamination = contaminant bases / host bases (0 for controls)
  aif / skani        = flagged bases / unflagged bases, a contig being flagged when its top hit's genus differs
                       from the bin's base-weighted majority genus (contigs without a hit are not flagged)
  CheckM2            = 'Contamination' / 100 (redundant marker genes)
  GUNC               = 'contamination_portion' at the level of maximal clade separation; GUNC's own call:
                       pass.GUNC == False (DIAMOND 2.2.8 with GUNC_SKIP_DIAMOND_VERSION_CHECK; GUNC pins 2.1.24)
Metrics: MAE of the estimate (pp), detection at a 5% contamination threshold (the MIMAG high-quality limit)
for every level, and false calls on controls.
Run time: aif and skani search all 400 bins' contigs as one file (seqtime; database in the page cache) at the
thread count of CheckM2 and GUNC; the gene-based tools are timed on the bin directory.
  python3 bench_real/contam_tools.py bins | checkm2 | gunc | seqtime | evaluate   (GTDB_SET selects the query set)
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gtdb_bench as gb  # noqa: E402
import gtdb_contam as gc  # noqa: E402

G = gb.G
SET = gb.SET
D = os.path.join(G, f"contam_bins{SET}")
CHECKM2 = os.environ.get("CHECKM2", "/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/checkm2")
CHECKM2_DB = os.environ.get("CHECKM2_DB", "/usr1/mrahmanh/meta-co-assembler/data/checkm2/CheckM2_database/uniref100.KO.1.dmnd")
GUNC = os.environ.get("GUNC", os.path.join(ROOT, "data", "envs", "gunc", "bin", "gunc"))
GUNC_DB = os.environ.get("GUNC_DB", os.path.join(ROOT, "data", "envs", "gunc_db", "gunc_db_gtdb214.dmnd"))
TOOLBIN = os.environ.get("TOOLBIN", "/home/mrahmanh/miniconda3/envs/binchicken/bin")  # diamond, prodigal
TH = os.environ.get("THREADS", "32")
NPER = 100


def timed(cmd, out=None, err=None, env=None):
    tf = tempfile.NamedTemporaryFile(delete=False).name
    with open(out or os.devnull, "w") as fo, open(err or os.devnull, "w") as fe:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf] + cmd, stdout=fo, stderr=fe, check=True, env=env)
    t, m = open(tf).read().split()[-2:]
    os.unlink(tf)
    return float(t), int(m)


def subset():
    by_g = gc.genomes()
    bins = gc.make_bins(by_g)
    out = []
    for lvl in gc.LEVELS:
        lb = [b for b in bins if b[0] == lvl][:NPER]
        for lv, h, c, host, add in lb:
            out.append((f"b{int(lv * 100):02d}_{h}_{c}", lv, h, c, host, add))
    for h in sorted(by_g)[:NPER]:
        out.append((f"ctrl_{h}", 0.0, h, "", by_g[h], []))
    return out


def bins():
    os.makedirs(os.path.join(D, "fa"), exist_ok=True)
    seqs = {}
    name = None
    for line in open(gb.Q):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    with open(os.path.join(D, "bins.tsv"), "w") as fh:
        fh.write("bin\tlevel\thost\tcontaminant\thost_bp\tcont_bp\tcont_contigs\n")
        for b, lv, h, c, host, add in subset():
            with open(os.path.join(D, "fa", b + ".fa"), "w") as fo:
                for ctg, _ in host + add:
                    fo.write(f">{ctg}\n{''.join(seqs[ctg])}\n")
            fh.write(f"{b}\t{lv}\t{h}\t{c}\t{sum(L for _, L in host)}\t{sum(L for _, L in add)}\t{','.join(x for x, _ in add)}\n")
    print(len(subset()), "bins written to", D)


def env():
    e = dict(os.environ)
    e["PATH"] = TOOLBIN + ":" + e["PATH"]
    return e


def checkm2():
    out = os.path.join(D, "checkm2")
    e = env()
    e["PATH"] = os.path.dirname(CHECKM2) + ":" + e["PATH"]  # CheckM2's own environment (prodigal, diamond)
    t, m = timed([CHECKM2, "predict", "--input", os.path.join(D, "fa"), "-x", "fa", "--output-directory", out, "--threads", TH,
                  "--database_path", CHECKM2_DB, "--force"], err=os.path.join(D, "checkm2.log"), env=e)
    json.dump({"s": t, "kb": m, "threads": int(TH)}, open(os.path.join(D, "checkm2_time.json"), "w"))
    print("checkm2", t, m)


def gunc():
    out = os.path.join(D, "gunc")
    tmp = os.path.join(D, "gunc_tmp")
    os.makedirs(out, exist_ok=True)
    os.makedirs(tmp, exist_ok=True)
    e = env()
    # GUNC 1.1.1 pins DIAMOND 2.1.24; the installed DIAMOND is 2.2.8, so GUNC's documented override is used
    e["GUNC_SKIP_DIAMOND_VERSION_CHECK"] = "1"
    t, m = timed([GUNC, "run", "--input_dir", os.path.join(D, "fa"), "--file_suffix", ".fa", "--out_dir", out, "--threads", TH,
                  "--db_file", GUNC_DB, "--temp_dir", tmp], err=os.path.join(D, "gunc.log"), env=e)
    json.dump({"s": t, "kb": m, "threads": int(TH)}, open(os.path.join(D, "gunc_time.json"), "w"))
    print("gunc", t, m)


def seqtime():
    """aif and skani on every bin's contigs (one FASTA, contig names prefixed by bin), timed, warm cache"""
    fa = os.path.join(D, "bins_all.fa")
    with open(fa, "w") as fo:
        for fn in sorted(os.listdir(os.path.join(D, "fa"))):
            b = fn[:-3]
            for line in open(os.path.join(D, "fa", fn)):
                fo.write(f">{b}__{line[1:]}" if line.startswith(">") else line)
    tp = os.path.join(D, "seq_time.json")
    out = json.load(open(tp)) if os.path.exists(tp) else {}
    runs = [("aif", [gb.AIF, "db", "search", os.path.join(G, "aif", "gtdb_0.88.aifd"), fa, "-n", "5", "-t", TH],
             [os.path.join(G, "aif", "gtdb_0.88.aifd")]),
            ("skani", [gb.SKANI, "search", "-d", os.path.join(G, "skani_db"), "--qi", fa, "-t", TH, "-o", os.path.join(D, "skani_bins.tsv")],
             [os.path.join(G, "skani_db", x) for x in os.listdir(os.path.join(G, "skani_db"))])]
    for key, cmd, warm in runs:
        if key in out:
            continue
        for w in warm:
            subprocess.run(["cat", w], stdout=subprocess.DEVNULL, check=True)
        o = os.path.join(D, f"{key}_bins.tsv")
        t, m = timed(cmd, out=None if key == "skani" else o, err=o + ".log")
        out[key] = {"s": t, "kb": m, "threads": int(TH)}
        print(key, out[key], flush=True)
        json.dump(out, open(tp, "w"), indent=1)


def evaluate():
    TX = gb.taxonomy()
    rows = list(csv.DictReader(open(os.path.join(D, "bins.tsv")), delimiter="\t"))
    contigs = {r["bin"]: [l.split()[0][1:] for l in open(os.path.join(D, "fa", r["bin"] + ".fa")) if l.startswith(">")] for r in rows}
    L = {r["contig"]: int(r["length"]) for r in csv.DictReader(open(os.path.join(G, "meta", f"query_contigs{SET}.tsv")), delimiter="\t")}
    est = {}
    for tool, H in (("aif", gb.load_aif(os.path.join(gb.RES, "aif_0.88.tsv"))), ("skani", gb.load_skani(os.path.join(gb.RES, "skani.tsv")))):
        top = {c: TX.get(hs[0][0], ("?", "?", "?"))[1] for c, hs in H.items() if hs}
        for r in rows:
            cs = contigs[r["bin"]]
            votes = {}
            for c in cs:
                if c in top:
                    votes[top[c]] = votes.get(top[c], 0) + L[c]
            maj = max(votes, key=votes.get) if votes else None
            fl = sum(L[c] for c in cs if c in top and top[c] != maj)
            tot = sum(L[c] for c in cs)
            est.setdefault(tool, {})[r["bin"]] = fl / max(1, tot - fl)
    cm = os.path.join(D, "checkm2", "quality_report.tsv")
    if os.path.exists(cm):
        for r in csv.DictReader(open(cm), delimiter="\t"):
            est.setdefault("checkm2", {})[r["Name"]] = float(r["Contamination"]) / 100
    gcall = {}
    gd = os.path.join(D, "gunc")
    gf = [f for f in os.listdir(gd) if f.endswith("maxCSS_level.tsv")] if os.path.isdir(gd) else []
    if gf:
        for r in csv.DictReader(open(os.path.join(gd, gf[0])), delimiter="\t"):
            est.setdefault("gunc", {})[r["genome"]] = float(r["contamination_portion"])
            gcall[r["genome"]] = r["pass.GUNC"] == "False"
    with open(os.path.join(gb.RES, "contam_tools_bins.tsv"), "w") as fh:
        tools = sorted(est)
        fh.write("bin\tlevel\ttrue\t" + "\t".join(tools) + "\tgunc_fail\n")
        for r in rows:
            b = r["bin"]
            fh.write(f"{b}\t{r['level']}\t{float(r['cont_bp']) / float(r['host_bp']):.5f}\t" +
                     "\t".join(f"{est[t][b]:.5f}" if b in est[t] else "NA" for t in tools) + f"\t{int(gcall.get(b, False))}\n")
    res = {}
    for tool, E in est.items():
        d = {}
        for lv in ("0.0", "0.05", "0.1", "0.2", "contaminated"):
            rs = [r for r in rows if (r["level"] == lv if lv != "contaminated" else float(r["level"]) > 0) and r["bin"] in E]
            if not rs:
                continue
            tru = {r["bin"]: float(r["cont_bp"]) / float(r["host_bp"]) for r in rs}
            mae = sum(abs(E[b] - t) for b, t in tru.items()) / len(rs)
            det = sum(1 for b in tru if E[b] >= 0.05) / len(rs)
            d[lv] = {"n": len(rs), "mae_pp": 100 * mae, "flag_5pct": det, "mean_est_pct": 100 * sum(E[b] for b in tru) / len(rs),
                     "mean_true_pct": 100 * sum(tru.values()) / len(rs)}
            if tool == "gunc":
                d[lv]["gunc_fail"] = sum(1 for b in tru if gcall.get(b)) / len(rs)
        res[tool] = d
    for f in ("checkm2_time.json", "gunc_time.json", "seq_time.json"):
        p = os.path.join(D, f)
        if os.path.exists(p):
            res[f[:-5]] = json.load(open(p))
    json.dump(res, open(os.path.join(gb.RES, "contam_tools.json"), "w"), indent=1)
    for tool, d in res.items():
        if "time" in tool:
            print(tool, d)
            continue
        print(tool, {lv: {k: round(v, 3) for k, v in x.items()} for lv, x in d.items()})


if __name__ == "__main__":
    {"bins": bins, "checkm2": checkm2, "gunc": gunc, "seqtime": seqtime, "evaluate": evaluate}[sys.argv[1]]()
