#!/usr/bin/env python3
"""Analysis of the HGT atlas (bench_real/hgt_atlas.py): shared windows -> blocks, genome pairs, network,
recency, replicons, and the ecology-versus-phylogeny test of Smillie et al. (2011).

Definitions
  call      a window of genome Q aligned to another species' representative P at ANI >= ANI_MIN over >= BP_MIN bases
  block     consecutive called windows of one contig of Q with the same partner P
  pair      (Q, P) with >= 1 block of >= PAIR_BP aligned bases; pairs whose calls cover >= SAME_FRAC of Q's windows
            are near-identical genomes or label conflicts and are set aside
  distance  rank of the lowest common GTDB taxon of Q and P

Usage: /usr/bin/python3 -I bench_real/hgt_analyze.py SET   -> results_hgt/SET/{summary.json, pairs.tsv.gz, blocks.tsv.gz}
"""
import collections
import csv
import gzip
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
ANI_MIN = float(os.environ.get("HGT_ANI", "0.99"))
BP_MIN = float(os.environ.get("HGT_BP", "500"))
PAIR_BP = float(os.environ.get("HGT_PAIR_BP", "1000"))
SAME_FRAC = 0.30
NCAND = int(os.environ.get("HGT_NCAND", "128"))  # search cap; a window with NCAND - 1 kept hits may have lost partners
PARTNERS = os.environ.get("HGT_PARTNERS", "all")
# reciprocal anchoring (full atlas): keep a pair only if both genomes carry a shared block of >= PAIR_BP on a contig with
# native sequence (a window not shared with that partner) or are complete genomes - contamination is a whole foreign contig
RECIPROCAL = os.environ.get("HGT_RECIPROCAL", "0") == "1"
# flank placement (bench_real/hgt_flanks.py) replaces the geometric anchoring; only blocks of >= PAIR_BP are tested and kept.
#   contig   a block in a draft or metagenome-assembled genome needs a flank window on its contig (it is not a whole contig)
#   uniform  contig, and for blocks between genera in ANY genome (complete ones too): no flank whose closest relative is
#            nearer the partner's lineage than the host's
#   lenient  contig, and the partner test for draft and metagenome-assembled genomes only
#   strict   lenient, and >= 1 flank placed on the host's side
FLANKS = os.environ.get("HGT_FLANKS", "0") == "1"
FLANK_MODE = os.environ.get("HGT_FLANK_MODE", "contig")
WMIN = 1000
W = 5000
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
DIST = ["same genus", "same family", "same order", "same class", "same phylum", "same domain", "different domain"]
HABS = ["gut", "host_other", "food", "engineered", "marine", "plant", "extreme", "subsurface", "freshwater", "soil",
        "sediment"]


def tax(t):
    return [x[3:] for x in t.split(";")]


def lca_class(a, b):
    """index into DIST of the lowest shared rank (species never shared between representatives)"""
    for i, r in enumerate(range(5, -1, -1)):  # genus, family, ..., domain
        if a[r] == b[r]:
            return i
    return 6


def load(setname):
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    for g in G.values():
        g["tax"] = tax(g["taxonomy"])
    base = os.path.join(HD, setname)
    contigs = collections.defaultdict(dict)  # acc -> ci -> (len, nwin)
    for l in open(os.path.join(base, "contigs.tsv")):
        acc, ci, nm, ln, nw = l.rstrip("\n").split("\t")
        contigs[acc][int(ci)] = (int(ln), int(nw))
    return G, contigs, base


def read_hits(base):
    hd = os.path.join(base, "hits")
    for f in sorted(os.listdir(hd)):
        if not f.endswith(".tsv"):
            continue
        for l in open(os.path.join(hd, f)):
            w, ref, rank, ani, af, theta, lr = l.rstrip("\n").split("\t")
            yield w, ref, int(rank), float(ani), float(af), float(lr)


def load_flanks(base):
    """flank window (acc, ci, start) -> refs of its best hits, best first (hgt_flanks.py filter output)"""
    fl = collections.defaultdict(list)
    hd = os.path.join(base, "flanks", "fhits")
    for f in sorted(os.listdir(hd)):
        if f.endswith(".tsv"):
            for l in open(os.path.join(hd, f)):
                w, ref, ani, af = l.rstrip("\n").split("\t")
                a, ci, st, ln = w.split("|")
                fl[(a, int(ci), int(st))].append(ref)
    return fl


def partner_ok(g):
    if PARTNERS == "all":
        return True
    if PARTNERS == "isolate":
        return g["category"] == "isolate"
    if PARTNERS == "complete":
        return g["category"] == "isolate" and g["assembly_level"] in ("Complete Genome", "Chromosome")
    raise SystemExit(f"HGT_PARTNERS={PARTNERS}?")


def main():
    setname = sys.argv[1]
    G, contigs, base = load(setname)
    out = os.path.join(ROOT, "results_hgt", setname + ("" if PARTNERS == "all" else "_" + PARTNERS)
                       + (("_flank" + ("" if FLANK_MODE == "lenient" else FLANK_MODE)) if FLANKS else "")
                       + ("_reciprocal" if RECIPROCAL else ""))
    os.makedirs(out, exist_ok=True)
    nwin = {a: sum(v[1] for v in c.values()) for a, c in contigs.items()}
    # chromosome = the largest replicon and any >= 1 Mbp; the rest are extrachromosomal (complete genomes only)
    chrom = {}
    for a, c in contigs.items():
        big = max(v[0] for v in c.values())
        for ci, (ln, _) in c.items():
            chrom[(a, ci)] = ln == big or ln >= 1_000_000
    # ---- calls
    calls = collections.defaultdict(list)  # (Q, P) -> [(ci, start, aligned, ani)]
    nhits = collections.Counter()
    ani_hist = collections.Counter()
    for w, ref, rank, ani, af, lr in read_hits(base):
        q, ci, st, ln = w.split("|")
        aligned = af * int(ln)
        nhits[w] += 1
        if ani >= ANI_MIN and aligned >= BP_MIN and ref in G and partner_ok(G[ref]):
            calls[(q, ref)].append((int(ci), int(st), aligned, ani))
    sat = sum(1 for v in nhits.values() if v >= NCAND - 1)
    # ---- pairs and blocks
    pairs, blocks, same_entity = [], [], []
    complete = {a for a, g in G.items() if g["assembly_level"] in ("Complete Genome", "Chromosome")
                and g["category"] == "isolate"}
    anchored = set()
    FL = load_flanks(base) if FLANKS else None
    fstat = collections.Counter()
    kept_blocks = set()  # (q, p, ci, s0) of blocks that pass the flank test (FLANKS mode, draft/MAG genomes)

    def flank_status(q, p, key, skip_partner=False):
        refs = [r for r in FL.get(key, ()) if not (skip_partner and r == p)]
        if not refs:
            return "unplaced"
        tr = G[refs[0]]["tax"]
        dq = lca_class(G[q]["tax"], tr)
        dp = -1 if refs[0] == p else lca_class(G[p]["tax"], tr)
        return "partner" if dp < dq else "host" if dq < dp else "neutral"

    def block_test(q, p, b):
        ci, s0, s1 = b[0], b[1], b[2]
        keys = []
        if s0 >= W:
            keys.append((q, ci, s0 - W))
        if contigs[q][ci][0] - (s1 + W) >= WMIN:
            keys.append((q, ci, s1 + W))
        if not keys:
            return "no flank", "no flank"
        st = [flank_status(q, p, k) for k in keys]
        st2 = [flank_status(q, p, k, True) for k in keys]
        v = "partner" if "partner" in st else "host" if "host" in st else "unplaced"
        v2 = "partner" if "partner" in st2 else "host" if "host" in st2 else "unplaced"
        return v, v2
    for (q, p), cs in calls.items():
        frac = len(cs) / max(1, nwin.get(q, 1))
        if frac >= SAME_FRAC:
            same_entity.append((q, p, frac))
            continue
        cs.sort()
        cur = None
        bl = []
        for ci, st, al, ani in cs:
            if cur and cur[0] == ci and st == cur[2] + W:
                cur[2] = st
                cur[3] += al
                cur[4] = max(cur[4], ani)
                cur[5].append(ani)
            else:
                if cur:
                    bl.append(cur)
                cur = [ci, st, st, al, ani, [ani]]
        bl.append(cur)
        tot = sum(b[3] for b in bl)
        big = max(b[3] for b in bl)
        if big < PAIR_BP:
            continue
        if FLANKS:
            cat = "complete" if q in complete else "draft"
            dd = lca_class(G[q]["tax"], G[p]["tax"])
            ok = False
            for b in bl:
                if b[3] < PAIR_BP:
                    continue
                v, v2 = block_test(q, p, b)
                if dd >= 1:
                    fstat[(cat, v)] += 1
                    fstat[(cat + "_ignoring_partner", v2)] += 1
                if q in complete:
                    passed = not (FLANK_MODE == "uniform" and dd >= 1 and v == "partner")
                elif v == "no flank":
                    passed = False
                elif dd >= 1 and FLANK_MODE in ("uniform", "lenient", "strict") and v == "partner":
                    passed = False
                elif dd >= 1 and FLANK_MODE == "strict" and v != "host":
                    passed = False
                else:
                    passed = True
                if passed:
                    kept_blocks.add((q, p, b[0], b[1]))
                    ok = True
            if ok:
                anchored.add((q, p))
        else:
            percontig = collections.Counter(c[0] for c in cs)
            if q in complete or any(b[3] >= PAIR_BP and contigs[q][b[0]][1] > percontig[b[0]] for b in bl):
                anchored.add((q, p))
        d = lca_class(G[q]["tax"], G[p]["tax"])
        pairs.append((q, p, d, len(cs), tot, max(c[3] for c in cs), len(bl)))
        for ci, s0, s1, al, mx, anis in bl:
            ln = contigs[q][ci][0]
            embedded = s0 > 0 and s1 + W < ln  # native sequence on both sides within the contig
            blocks.append((q, p, d, ci, s0, min(s1 + W, ln), al, sum(anis) / len(anis), mx,
                           int(chrom.get((q, ci), True)), int(embedded)))
    if RECIPROCAL:
        keep = {(q, p) for (q, p) in anchored if (p, q) in anchored}
        S0 = {"pairs_before_reciprocal": len(pairs), "anchored_pairs": len(anchored), "reciprocal_pairs": len(keep)}
        pairs = [r for r in pairs if (r[0], r[1]) in keep]
        blocks = [b for b in blocks if (b[0], b[1]) in keep]
    elif FLANKS:
        keep = anchored
        S0 = {"pairs_before_flanks": len(pairs), "anchored_pairs": len(anchored)}
        pairs = [r for r in pairs if (r[0], r[1]) in keep]
        blocks = [b for b in blocks if (b[0], b[1]) in keep]
    else:
        keep = None
        S0 = {"anchored_pairs": len(anchored)}
    if FLANKS:
        blocks = [b for b in blocks if (b[0], b[1], b[3], b[4]) in kept_blocks]
        S0["flank_mode"] = FLANK_MODE
        S0["flank_block_status"] = {f"{a}|{b}": n for (a, b), n in sorted(fstat.items())}
    with gzip.open(os.path.join(out, "pairs.tsv.gz"), "wt") as f:
        f.write("query\tpartner\tdistance\twindows\taligned_bp\tmax_ani\tblocks\n")
        for r in pairs:
            f.write(f"{r[0]}\t{r[1]}\t{DIST[r[2]]}\t{r[3]}\t{r[4]:.0f}\t{r[5]:.4f}\t{r[6]}\n")
    with gzip.open(os.path.join(out, "blocks.tsv.gz"), "wt") as f:
        f.write("query\tpartner\tdistance\tcontig\tstart\tend\taligned_bp\tmean_ani\tmax_ani\tchromosome\tembedded\n")
        for r in blocks:
            f.write(f"{r[0]}\t{r[1]}\t{DIST[r[2]]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.0f}\t{r[7]:.4f}\t{r[8]:.4f}\t{r[9]}\t{r[10]}\n")
    # ---- summary
    S = {"set": setname, "partners": PARTNERS, "ani_min": ANI_MIN, "bp_min": BP_MIN, "pair_bp": PAIR_BP, "genomes": len(contigs),
         "windows": sum(nwin.values()), "windows_with_hits": len(nhits), "windows_saturated": sat,
         "same_entity_pairs": len(same_entity), "reciprocal": RECIPROCAL, **S0}
    Q = set(contigs)
    by_d = collections.Counter(r[2] for r in pairs)
    S["pairs_by_distance"] = {DIST[i]: by_d[i] for i in range(7)}
    xg = [r for r in pairs if r[2] >= 1]
    S["cross_genus_pairs"] = len(xg)
    S["genomes_with_cross_genus"] = len({r[0] for r in xg})
    S["genomes_with_cross_family"] = len({r[0] for r in pairs if r[2] >= 2})
    S["genomes_with_cross_phylum"] = len({r[0] for r in pairs if r[2] >= 5})
    gp = {(G[r[0]]["tax"][5], G[r[1]]["tax"][5]) for r in xg}
    S["genus_pairs"] = len({tuple(sorted(x)) for x in gp})
    S["genera_involved"] = len({x for p in gp for x in p})
    xb = [b for b in blocks if b[2] >= 1]
    S["cross_genus_blocks"] = len(xb)
    S["cross_genus_block_bp"] = sum(b[6] for b in xb)
    # share of each query genome's windows that carry cross-genus shared sequence (union over partners)
    shared_w = collections.defaultdict(set)
    se = {(a, b) for a, b, _ in same_entity}
    for (q, p), cs in calls.items():
        if (RECIPROCAL or FLANKS) and (q, p) not in keep:
            continue
        if lca_class(G[q]["tax"], G[p]["tax"]) >= 1 and (q, p) not in se:
            for ci, st, al, ani in cs:
                shared_w[q].add((ci, st))
    fr = sorted(len(shared_w[q]) / nwin[q] for q in Q)
    S["shared_window_frac_quantiles"] = {str(x): fr[min(len(fr) - 1, int(x * len(fr)))] for x in (0.5, 0.9, 0.99)}
    S["shared_window_frac_mean"] = sum(fr) / len(fr)
    S["blocks_extrachromosomal_share"] = sum(1 for b in xb if not b[9]) / max(1, len(xb))
    xbc = [b for b in xb if b[0] in complete]  # replicons are known only in complete genomes
    S["blocks_extrachromosomal_share_complete"] = sum(1 for b in xbc if not b[9]) / max(1, len(xbc))
    S["cross_genus_blocks_complete"] = len(xbc)
    # genomes with cross-genus / cross-phylum sharing by category of the query genome
    def gcat(a):
        g = G[a]
        return "complete isolate" if a in complete else "draft isolate" if g["category"] == "isolate" else g["category"]
    den_c = collections.Counter(gcat(a) for a in contigs)
    xgq = {r[0] for r in xg}
    xpq = {r[0] for r in pairs if r[2] >= 5}
    S["by_category"] = {c: {"genomes": n, "with_cross_genus": sum(1 for a in xgq if gcat(a) == c),
                            "with_cross_phylum": sum(1 for a in xpq if gcat(a) == c)} for c, n in den_c.items()}
    S["blocks_embedded_share"] = sum(1 for b in xb if b[10]) / max(1, len(xb))
    lens = sorted(b[5] - b[4] for b in xb)
    S["block_span_quantiles"] = {str(x): lens[min(len(lens) - 1, int(x * len(lens)))] for x in (0.1, 0.5, 0.9, 0.99)} if lens else {}
    # recency: identity of cross-genus blocks
    S["block_ani_bins"] = dict(collections.Counter(
        ("100" if b[7] >= 0.9995 else "99.9-99.95" if b[7] >= 0.999 else "99.5-99.9" if b[7] >= 0.995 else "99-99.5")
        for b in xb))
    # partner categories
    S["cross_genus_partner_category"] = dict(collections.Counter(G[r[1]]["category"] for r in xg))
    # hubs
    deg = collections.Counter()
    for a, b in {tuple(sorted(x)) for x in gp}:
        deg[a] += 1
        deg[b] += 1
    S["top_genera_by_partner_genera"] = sorted(deg.items(), key=lambda x: (-x[1], x[0]))[:25]
    phy = collections.Counter(tuple(sorted((G[r[0]]["tax"][1], G[r[1]]["tax"][1]))) for r in pairs if r[2] >= 5)
    S["cross_phylum_pairs_top"] = [[f"{a} / {b}", n] for (a, b), n in sorted(phy.items(), key=lambda x: (-x[1], x[0]))[:20]]
    S["ecology"] = ecology(G, Q, pairs)
    json.dump(S, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in S.items() if k not in ("ecology", "top_genera_by_partner_genera", "cross_phylum_pairs_top")}, indent=1))
    print(json.dumps(S["ecology"]["by_distance"], indent=1))


def ecology(G, Q, pairs):
    """HGT rate of directed pairs (Q in the query set, P any representative), both with a habitat label, by the
    rank of their lowest common taxon and by habitat sameness. Denominators count all such pairs exactly."""
    lab = {a: g["habitat"] for a, g in G.items() if g["habitat"] in HABS and partner_ok(g)}
    Ql = [a for a in Q if a in lab]
    Al = list(lab)
    # pairs sharing the rank-r taxon (r = genus .. domain), split by same habitat; self pairs removed
    def share_counts(r):
        nq, na = collections.Counter(), collections.Counter()
        nqh, nah = collections.Counter(), collections.Counter()
        for a in Ql:
            t = tuple(G[a]["tax"][:r + 1])
            nq[t] += 1
            nqh[(t, lab[a])] += 1
        for a in Al:
            t = tuple(G[a]["tax"][:r + 1])
            na[t] += 1
            nah[(t, lab[a])] += 1
        allp = sum(n * na[t] for t, n in nq.items()) - len(Ql)
        samep = sum(n * nah[k] for k, n in nqh.items()) - len(Ql)
        return allp, samep
    tot_all = len(Ql) * len(Al) - len(Ql)
    hq = collections.Counter(lab[a] for a in Ql)
    ha = collections.Counter(lab[a] for a in Al)
    tot_same = sum(hq[h] * ha[h] for h in hq) - len(Ql)
    share = {r: share_counts(r) for r in range(0, 6)}  # r = 0 domain ... 5 genus
    den = {}
    prev_all, prev_same = tot_all, tot_same
    # DIST index i: 0 same genus (share r=5), 1 same family not genus (share 4 - share 5), ..., 5 same domain, 6 diff domain
    for i in range(7):
        r = 5 - i
        cur_all, cur_same = share[r] if r >= 0 else (0, 0)
        nxt_all, nxt_same = share[r + 1] if 0 <= r + 1 <= 5 else (0, 0)
        if i == 0:
            a_, s_ = cur_all, cur_same
        elif i == 6:
            a_, s_ = tot_all - share[0][0], tot_same - share[0][1]
        else:
            a_, s_ = cur_all - nxt_all, cur_same - nxt_same
        den[i] = (a_, s_)
    num = collections.Counter()
    for q, p, d, *_ in pairs:
        if q in lab and p in lab:
            num[(d, lab[q] == lab[p])] += 1
    res = {}
    for i in range(7):
        a_, s_ = den[i]
        d_ = a_ - s_
        ns, nd = num[(i, True)], num[(i, False)]
        rs = ns / s_ if s_ else float("nan")
        rd = nd / d_ if d_ else float("nan")
        res[DIST[i]] = {"same_habitat_pairs": s_, "same_habitat_hgt": ns, "rate_same": rs,
                        "diff_habitat_pairs": d_, "diff_habitat_hgt": nd, "rate_diff": rd,
                        "ratio": rs / rd if rd and rd == rd and rd > 0 else None}
    # per habitat (cross-genus): rate within habitat, and habitat-pair matrix of observed / expected
    mat = collections.Counter()
    for q, p, d, *_ in pairs:
        if d >= 1 and q in lab and p in lab:
            mat[(lab[q], lab[p])] += 1
    within = {}
    for h in HABS:
        within[h] = {"query_genomes": hq[h], "partners": ha[h], "cross_genus_pairs": mat[(h, h)]}
    return {"labelled_queries": len(Ql), "labelled_all": len(Al), "by_distance": res, "within_habitat": within,
            "habitat_matrix": {f"{a}|{b}": n for (a, b), n in mat.items()}, "hq": dict(hq), "ha": dict(ha)}


if __name__ == "__main__":
    main()
