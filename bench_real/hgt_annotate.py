#!/usr/bin/env python3
"""What moves between genera: genes on shared blocks against genes on background windows of the same genomes.

Regions: the union, per query contig, of its cross-genus blocks (results_hgt/TAG/blocks.tsv.gz), each tagged with
the most distant partner. Background: random 5 kb windows of the same genomes outside any region, matched in
total length per genome. Genes: prodigal -p meta. Function: DIAMOND against CheckM2's UniRef100 database with KEGG
orthologs (best hit, e-value <= 1e-10, >= 30% identity over >= 50% of the query), classes from the KEGG hierarchy.

Usage: /usr/bin/python3 -I bench_real/hgt_annotate.py TAG     (TAG e.g. complete_complete)
"""
import collections
import csv
import gzip
import json
import math
import multiprocessing
import os
import random
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
PRODIGAL = os.path.join(ROOT, "data", "envs", "compare", "bin", "prodigal")
DIAMOND = os.path.join(ROOT, "data", "envs", "compare", "bin", "diamond")
KODB = os.environ.get("CHECKM2_DB", "/usr1/mrahmanh/meta-co-assembler/data/checkm2/CheckM2_database/uniref100.KO.1.dmnd")
W = 5000
DORDER = ["same genus", "same family", "same order", "same class", "same phylum", "same domain", "different domain"]

CLASS_RULES = [  # (class, BRITE ids, name regex); first match wins
    ("transposase", (), r"transposase|insertion element|\bIS\d"),
    ("integrase/recombinase", (), r"integrase|recombinase|resolvase|invertase|excisionase"),
    ("antimicrobial resistance", ("01504",), r"beta-lactamase|tetracycline resistance|aminoglycoside|chloramphenicol|"
                                            r"vancomycin|macrolide|streptogramin|lincosamide|sulfonamide|trimethoprim|"
                                            r"fosfomycin|multidrug|bleomycin|quinolone"),
    ("conjugation / T4SS", (), r"conjugal|conjugative|type IV secretion|relaxase|mobilization|\btra[A-Z]\b|\btrb[A-Z]\b|"
                              r"\bvir[BD]\d"),
    ("phage", ("03200",), r"phage|terminase|portal protein|capsid|tail fiber|tail protein|holin|baseplate"),
    ("defense (RM, CRISPR, TA)", ("02048",), r"restriction|CRISPR|\bcas\d|toxin-antitoxin|antitoxin|abortive"),
    ("metal resistance", (), r"mercur|arsen|cadmium|copper resistance|silver|cobalt-zinc-cadmium|czc|chromate|"
                             r"tellurite|heavy metal"),
    ("toxins", ("02042",), r"toxin"),
    ("transporters", ("02000",), r"transporter|permease|ABC|efflux"),
    ("ribosome / translation", ("03011", "03012", "03016"), r"ribosomal protein|translation"),
    ("DNA replication / repair", ("03032", "03400"), r"DNA polymerase|helicase|repair"),
]


def kegg():
    names = {}
    for l in open(os.path.join(HD, "kegg", "ko_list.tsv")):
        k, n = l.rstrip("\n").split("\t", 1)
        names[k] = n
    brite = collections.defaultdict(set)
    d = json.load(open(os.path.join(HD, "kegg", "ko00001.json")))

    def walk(node, ctx):
        nm = node["name"]
        m = re.search(r"\[BR:ko(\d+)\]", nm)
        c2 = ctx + ([m.group(1)] if m else [])
        if nm.startswith("K") and re.match(r"K\d{5}", nm):
            for b in c2:
                brite[nm[:6]].add(b)
        for ch in node.get("children", []):
            walk(ch, c2)
    walk(d, [])
    return names, brite


def classify(ko, names, brite):
    if ko is None:
        return "no KEGG ortholog"
    n = names.get(ko, "")
    for cls, br, rx in CLASS_RULES:
        if (set(br) & brite.get(ko, set())) or re.search(rx, n, re.I):
            return cls
    return "other KEGG ortholog"


def read_fasta_gz(path):
    seqs, buf = [], []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l[0] == ">":
                if buf:
                    seqs.append("".join(buf))
                    buf = []
                seqs.append(None)
            else:
                buf.append(l.strip())
    if buf:
        seqs.append("".join(buf))
    return [s for s in seqs if s is not None]


def extract(job):
    acc, path, regions, nbg, seed = job
    seqs = read_fasta_gz(path)
    out, bg = [], []
    occupied = collections.defaultdict(list)
    for ci, s0, s1, dmax, ani in regions:
        out.append(f">R|{acc}|{ci}|{s0}|{s1}|{dmax}|{ani:.4f}\n{seqs[ci][s0:s1]}\n")
        occupied[ci].append((s0, s1))
    rnd = random.Random(seed)
    tries = 0
    while len(bg) < nbg and tries < 200 * nbg:
        tries += 1
        ci = rnd.randrange(len(seqs))
        if len(seqs[ci]) < W:
            continue
        st = rnd.randrange(0, len(seqs[ci]) - W + 1)
        if any(st < e and st + W > s for s, e in occupied[ci]):
            continue
        occupied[ci].append((st, st + W))
        bg.append(f">B|{acc}|{ci}|{st}|{st + W}\n{seqs[ci][st:st + W]}\n")
    return "".join(out), "".join(bg)


def prodigal(fa):
    faa = fa.replace(".fa", ".faa")
    if not os.path.exists(faa):
        subprocess.run([PRODIGAL, "-p", "meta", "-q", "-i", fa, "-a", faa + ".tmp", "-o", "/dev/null"], check=True)
        os.replace(faa + ".tmp", faa)
    return faa


def main():
    tag = sys.argv[1]
    res = os.path.join(ROOT, "results_hgt", tag)
    wd = os.path.join(HD, "annot_" + tag)
    os.makedirs(wd, exist_ok=True)
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    # union of cross-genus blocks per query contig
    iv = collections.defaultdict(list)
    for b in csv.DictReader(gzip.open(os.path.join(res, "blocks.tsv.gz"), "rt"), delimiter="\t"):
        if b["distance"] == "same genus":
            continue
        iv[(b["query"], int(b["contig"]))].append((int(b["start"]), int(b["end"]), DORDER.index(b["distance"]),
                                                    float(b["max_ani"])))
    regions = collections.defaultdict(list)
    for (acc, ci), lst in iv.items():
        lst.sort()
        cur = None
        for s, e, d, a in lst:
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
                cur[2] = max(cur[2], d)
                cur[3] = max(cur[3], a)
            else:
                if cur:
                    regions[acc].append((ci, *cur))
                cur = [s, e, d, a]
        regions[acc].append((ci, *cur))
    jobs = []
    for i, (acc, rg) in enumerate(sorted(regions.items())):
        bp = sum(r[2] - r[1] for r in rg)
        jobs.append((acc, G[acc]["path"], rg, max(1, round(bp / W)), i))
    nch = 48
    chunks = [[] for _ in range(nch)]
    for i, j in enumerate(jobs):
        chunks[i % nch].append(j)
    fas = [os.path.join(wd, f"c{c:02d}.fa") for c in range(nch)]
    if not all(os.path.exists(f) for f in fas):
        with multiprocessing.get_context("fork").Pool(32) as pool:
            outs = [open(f, "w") for f in fas]
            for i, (r, b) in enumerate(pool.imap(extract, jobs, chunksize=4)):
                outs[i % nch].write(r)
                outs[i % nch].write(b)
            for o in outs:
                o.close()
    with multiprocessing.get_context("fork").Pool(48) as pool:
        faas = pool.map(prodigal, fas)
    allp = os.path.join(wd, "proteins.faa")
    if not os.path.exists(allp):
        with open(allp + ".tmp", "w") as fo:
            for f in faas:
                fo.write(open(f).read())
        os.replace(allp + ".tmp", allp)
    hits = os.path.join(wd, "ko.tsv")
    if not os.path.exists(hits):
        subprocess.run([DIAMOND, "blastp", "--db", KODB, "-q", allp, "-o", hits + ".tmp", "--max-target-seqs", "1",
                        "--evalue", "1e-10", "--id", "30", "--query-cover", "50", "-p", "48", "--outfmt", "6", "qseqid",
                        "sseqid", "pident", "length", "qlen", "evalue"], check=True, capture_output=True)
        os.replace(hits + ".tmp", hits)
    ko = {}
    for l in open(hits):
        q, s = l.split("\t")[:2]
        if q not in ko and "~" in s:
            ko[q] = s.split("~")[1]
    names, brite = kegg()
    genes = []  # (kind R/B, region id, dmax, ko, class)
    for l in open(allp):
        if l[0] != ">":
            continue
        gid = l[1:].split()[0]
        reg = gid.rsplit("_", 1)[0]
        p = reg.split("|")
        k = ko.get(gid)
        genes.append((p[0], reg, int(p[5]) if p[0] == "R" else -1, k, classify(k, names, brite)))
    R = [g for g in genes if g[0] == "R"]
    B = [g for g in genes if g[0] == "B"]
    cls = [c for c, _, _ in CLASS_RULES] + ["other KEGG ortholog", "no KEGG ortholog"]
    cr, cb = collections.Counter(g[4] for g in R), collections.Counter(g[4] for g in B)
    table = []
    for c in cls:
        a, b = cr[c], cb[c]
        orr = ((a + 0.5) / (len(R) - a + 0.5)) / ((b + 0.5) / (len(B) - b + 0.5))
        table.append({"class": c, "shared": a, "shared_share": a / len(R), "background": b, "background_share": b / len(B),
                      "odds_ratio": orr})
    # by distance of the farthest partner
    byd = {}
    for d in range(1, 7):
        gd = [g for g in R if g[2] == d]
        if gd:
            c = collections.Counter(g[4] for g in gd)
            byd[DORDER[d]] = {"genes": len(gd), **{k: c[k] / len(gd) for k in cls}}
    topko = collections.Counter(g[3] for g in R if g[3])
    bko = collections.Counter(g[3] for g in B if g[3])
    tops = [{"ko": k, "name": names.get(k, ""), "shared": n, "background": bko[k]} for k, n in topko.most_common(40)]
    nreg = sum(len(v) for v in regions.values())
    S = {"regions": nreg, "region_bp": sum(r[2] - r[1] for v in regions.values() for r in v), "genomes": len(regions),
         "genes_shared": len(R), "genes_background": len(B), "classes": table, "by_distance": byd, "top_kos": tops}
    json.dump(S, open(os.path.join(res, "annotation.json"), "w"), indent=1)
    for t in table:
        print(f"{t['class']:28s} {t['shared']:7d} {100 * t['shared_share']:5.1f}%  {t['background']:7d} "
              f"{100 * t['background_share']:5.1f}%  OR {t['odds_ratio']:.2f}")
    for t in tops[:25]:
        print(t)


if __name__ == "__main__":
    main()
