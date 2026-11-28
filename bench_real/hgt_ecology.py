#!/usr/bin/env python3
"""Ecology versus phylogeny in the HGT atlas (Smillie et al. 2011, Nature 480:241), on unordered genome pairs.

Universe: genomes of the query set that pass the partner filter (e.g. complete isolates) and have a habitat label.
A pair exchanges if either genome has a block of >= 1 kb at >= 99% with the other (results_hgt/TAG/pairs.tsv.gz).
Pairs are stratified by the rank of their lowest common GTDB taxon; within each stratum the exchange rate of
same-habitat pairs is compared with that of different-habitat pairs. The Mantel-Haenszel odds ratio summarises the
strata above the genus; a bootstrap over genomes (pairs share genomes) gives its interval. Habitat-pair enrichment:
observed exchanging pairs over the number expected from the stratum rates.

HGT_ECO_XSTUDY=1 leaves out pairs of genomes from the same NCBI BioProject (numerators and denominators), so that
genomes assembled from the same samples cannot share contigs by cross-binning or contamination. RefSeq assemblies all
carry the umbrella project PRJNA224116 and count as their own project (data/hgt/bioproject.tsv, from GTDB metadata).

HGT_ECO_XSAMPLE=1 also leaves out pairs of MAGs whose source samples overlap (data/hgt/provenance/sources.tsv, from
hgt_provenance.py): different projects reassemble the same public metagenomes, which BioProjects do not reveal.
HGT_ECO_KNOWN=1 restricts the universe to MAGs with a known source sample. With XSAMPLE, results also report the
exchange rate of co-sample pairs against all other pairs in each stratum ("cosample").

Usage: /usr/bin/python3 -I bench_real/hgt_ecology.py TAG SET [PARTNERS]   e.g. complete_complete complete complete
       PARTNERS: complete | all | isolate | mag
"""
import collections
import csv
import gzip
import json
import os
import random
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
DIST = ["same genus", "same family", "same order", "same class", "same phylum", "same domain", "different domain"]
HABS = ["gut", "host_other", "food", "engineered", "marine", "plant", "extreme", "subsurface", "freshwater", "soil",
        "sediment"]


def bioprojects():
    """representative accession -> NCBI BioProject (data/hgt/bioproject.tsv, built from the GTDB metadata once)"""
    path = os.path.join(HD, "bioproject.tsv")
    if not os.path.exists(path):
        with open(path + ".tmp", "w") as fo:
            for f in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
                with gzip.open(os.path.join(ROOT, "data", "gtdb", f), "rt") as fi:
                    h = fi.readline().rstrip("\n").split("\t")
                    ia, ib, ir = h.index("accession"), h.index("ncbi_bioproject"), h.index("gtdb_representative")
                    for l in fi:
                        x = l.rstrip("\n").split("\t")
                        if x[ir] in ("t", "True", "TRUE"):
                            fo.write(f"{x[ia][3:]}\t{x[ib]}\n")
        os.replace(path + ".tmp", path)
    return dict(l.rstrip("\n").split("\t") for l in open(path))


def main():
    tag, setname = sys.argv[1], sys.argv[2]
    partners = sys.argv[3] if len(sys.argv) > 3 else "complete"
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    members = {l.split("\t")[0] for l in open(os.path.join(HD, setname, "contigs.tsv"))}
    def ok(g):
        if partners == "complete":
            return g["category"] == "isolate" and g["assembly_level"] in ("Complete Genome", "Chromosome")
        if partners == "isolate":
            return g["category"] == "isolate"
        if partners == "mag":
            return g["category"] != "isolate"
        return True
    xsample = os.environ.get("HGT_ECO_XSAMPLE") == "1"
    known = os.environ.get("HGT_ECO_KNOWN") == "1"
    SRC = {}
    if xsample or known:
        for l in open(os.path.join(HD, "provenance", "sources.tsv")):
            a, s = l.rstrip("\n").split("\t")
            SRC[a] = s.split(",")
    U = sorted(a for a in members if ok(G[a]) and G[a]["habitat"] in HABS and (not known or a in SRC))
    n = len(U)
    idx = {a: i for i, a in enumerate(U)}
    tax = [[x[3:] for x in G[a]["taxonomy"].split(";")] for a in U]
    hab = [G[a]["habitat"] for a in U]
    xstudy = os.environ.get("HGT_ECO_XSTUDY") == "1"
    if xstudy:
        bpm = bioprojects()
        grp = [bpm[a] if bpm.get(a) not in (None, "", "none", "PRJNA224116") else "_" + a for a in U]
    else:
        grp = ["_" + a for a in U]
    # integer codes: taxon prefix per rank (0 domain .. 5 genus), habitat, group
    T = np.zeros((6, n), dtype=np.int64)
    for r in range(6):
        ids = {}
        T[r] = [ids.setdefault(tuple(t[:r + 1]), len(ids)) for t in tax]
    H = np.array([HABS.index(h) for h in hab], dtype=np.int64)
    gid = {}
    GID = np.array([gid.setdefault(g, len(gid)) for g in grp], dtype=np.int64)

    def lca_v(I, J):
        d = np.full(len(I), 6, dtype=np.int64)
        for r in range(6):
            d[T[r][I] == T[r][J]] = 5 - r
        return d
    # co-sample pairs (unordered, packed i * n + j with i < j)
    cos = np.zeros(0, dtype=np.int64)
    if xsample:
        by = collections.defaultdict(set)
        for a in U:
            for s in SRC.get(a, ()):
                by[s].add(idx[a])
        chunks = []
        for L in by.values():
            if len(L) > 1:
                L = np.array(sorted(L), dtype=np.int64)
                ii, jj = np.triu_indices(len(L), 1)
                chunks.append(L[ii] * n + L[jj])
        if chunks:
            cos = np.unique(np.concatenate(chunks))
    CI_, CJ_ = cos // max(n, 1), cos % max(n, 1)
    cos_d, cos_same = lca_v(CI_, CJ_), H[CI_] == H[CJ_]
    # extra exclusions: co-sample pairs not already within one group
    keep = GID[CI_] != GID[CJ_]
    XI, XJ, Xd, Xs = CI_[keep], CJ_[keep], cos_d[keep], cos_same[keep]
    xcode = Xd * 2 + Xs
    xkey = np.unique(cos[keep])
    ex_all = set()
    for r in csv.DictReader(gzip.open(os.path.join(ROOT, "results_hgt", tag, "pairs.tsv.gz"), "rt"), delimiter="\t"):
        if r["query"] in idx and r["partner"] in idx:
            i, j = idx[r["query"]], idx[r["partner"]]
            ex_all.add((min(i, j), max(i, j)))
    ex_arr = np.array(sorted(ex_all), dtype=np.int64).reshape(-1, 2)
    ex_key = ex_arr[:, 0] * n + ex_arr[:, 1]
    ex_ok = GID[ex_arr[:, 0]] != GID[ex_arr[:, 1]]
    ex_ok &= ~np.isin(ex_key, xkey)
    ex = {(int(i), int(j)) for i, j in ex_arr[ex_ok]}
    cosample = None
    if xsample:
        cex = np.isin(cos, ex_key)
        cosample = {}
        for d in range(7):
            m = cos_d == d
            cosample[DIST[d]] = {"pairs": int(m.sum()), "exchanging": int((cex & m).sum()),
                                 "pairs_same_habitat": int((m & cos_same).sum()),
                                 "exchanging_same_habitat": int((cex & m & cos_same).sum())}

    def lca(i, j):
        for d, r in enumerate(range(5, -1, -1)):
            if tax[i][r] == tax[j][r]:
                return d
        return 6
    exl = [(i, j, lca(i, j), hab[i] == hab[j], hab[i], hab[j]) for i, j in ex]

    def tables(w):
        """weighted pair counts per (stratum, same habitat) and exchanging pairs; w = genome weights (bootstrap)"""
        den = collections.Counter()
        # pairs sharing the rank-r taxon (r = 5 genus .. 0 domain), all and same-habitat, unordered, no self pairs
        share = {}
        def cpairs(keyf):
            """weighted unordered pairs of genomes with equal key, minus those within one group (BioProject)"""
            n, n2, ng, ng2 = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
            for i, x in enumerate(w):
                if not x:
                    continue
                k = keyf(i)
                n[k] += x
                n2[k] += x * x
                if xstudy:
                    ng[(k, grp[i])] += x
                    ng2[(k, grp[i])] += x * x
            return sum((n[k] ** 2 - n2[k]) / 2 for k in n) - sum((ng[k] ** 2 - ng2[k]) / 2 for k in ng)
        tot2 = cpairs(lambda i: 0)
        same_tot = cpairs(lambda i: hab[i])
        for r in range(6):
            share[r] = (cpairs(lambda i: tuple(tax[i][:r + 1])), cpairs(lambda i: (tuple(tax[i][:r + 1]), hab[i])))
        for d in range(7):
            r = 5 - d
            if d == 0:
                a_, s_ = share[5]
            elif d == 6:
                a_, s_ = tot2 - share[0][0], same_tot - share[0][1]
            else:
                a_, s_ = share[r][0] - share[r + 1][0], share[r][1] - share[r + 1][1]
            den[(d, True)] = s_
            den[(d, False)] = a_ - s_
        if len(XI):
            wa = np.asarray(w, dtype=float)
            sub = np.bincount(xcode, weights=wa[XI] * wa[XJ], minlength=14)
            for d in range(7):
                den[(d, True)] -= sub[2 * d + 1]
                den[(d, False)] -= sub[2 * d]
        num = collections.Counter()
        for i, j, d, same, _, _ in exl:
            num[(d, same)] += w[i] * w[j]
        return den, num

    def mh(den, num, strata):
        nume = deno = 0.0
        for d in strata:
            a, b = num[(d, True)], den[(d, True)] - num[(d, True)]
            c, e = num[(d, False)], den[(d, False)] - num[(d, False)]
            n = a + b + c + e
            if n <= 0:
                continue
            nume += a * e / n
            deno += b * c / n
        return nume / deno if deno else float("nan")

    w1 = [1] * len(U)
    den, num = tables(w1)
    res = {"genomes": len(U), "habitats": dict(collections.Counter(hab)), "exchanging_pairs": len(ex), "strata": {}}
    if xsample:
        res["cosample_pairs"] = int(len(cos))
        res["cosample_pairs_excluded_beyond_project"] = int(len(XI))
        res["genomes_with_source"] = sum(1 for a in U if a in SRC)
        # exchange rate of co-sample pairs against the pairs analysed (not co-sample, and not same-project if XSTUDY)
        for d in range(7):
            c = cosample[DIST[d]]
            allp = den[(d, True)] + den[(d, False)]
            alle = num[(d, True)] + num[(d, False)]
            c["rate"] = c["exchanging"] / c["pairs"] if c["pairs"] else None
            c["rate_same_habitat"] = (c["exchanging_same_habitat"] / c["pairs_same_habitat"]
                                      if c["pairs_same_habitat"] else None)
            c["rate_other_pairs"] = alle / allp if allp else None
            c["rate_other_same_habitat"] = num[(d, True)] / den[(d, True)] if den[(d, True)] else None
        res["cosample"] = cosample
    for d in range(7):
        s_, df_ = den[(d, True)], den[(d, False)]
        res["strata"][DIST[d]] = {"same_pairs": s_, "same_exchanging": num[(d, True)], "diff_pairs": df_,
                                  "diff_exchanging": num[(d, False)],
                                  "rate_same": num[(d, True)] / s_ if s_ else None,
                                  "rate_diff": num[(d, False)] / df_ if df_ else None}
    cross = [1, 2, 3, 4, 5, 6]
    res["mh_odds_cross_genus"] = mh(den, num, cross)
    res["mh_odds_cross_family"] = mh(den, num, [2, 3, 4, 5, 6])
    res["mh_odds_same_genus"] = mh(den, num, [0])
    rng = random.Random(20261127)
    boots = collections.defaultdict(list)
    for b in range(200):
        w = [0] * len(U)
        for _ in range(len(U)):
            w[rng.randrange(len(U))] += 1
        dB, nB = tables(w)
        boots["cross_genus"].append(mh(dB, nB, cross))
        boots["cross_family"].append(mh(dB, nB, [2, 3, 4, 5, 6]))
        boots["same_genus"].append(mh(dB, nB, [0]))
        for d in range(7):
            rs = nB[(d, True)] / dB[(d, True)] if dB[(d, True)] else float("nan")
            rd = nB[(d, False)] / dB[(d, False)] if dB[(d, False)] else float("nan")
            boots[f"ratio_{d}"].append(rs / rd if rd and rd == rd and rd > 0 else float("nan"))
    def ci(v):
        v = sorted(x for x in v if x == x)
        return [v[int(0.025 * len(v))], v[int(0.975 * len(v)) - 1]] if v else None
    res["mh_ci"] = {k: ci(boots[k]) for k in ("cross_genus", "cross_family", "same_genus")}
    for d in range(7):
        st = res["strata"][DIST[d]]
        st["ratio"] = st["rate_same"] / st["rate_diff"] if st["rate_same"] is not None and st["rate_diff"] else None
        st["ratio_ci"] = ci(boots[f"ratio_{d}"])
    cnt, cntg = {}, {}
    for r in range(-1, 6):  # r = -1: no shared taxon required
        c, cg = collections.Counter(), collections.Counter()
        for i in range(len(U)):
            t = tuple(tax[i][:r + 1])
            c[(t, hab[i])] += 1
            if xstudy:
                cg[(t, hab[i], grp[i])] += 1
        cnt[r], cntg[r] = c, cg

    def pairs_between(h1, h2, r):
        """unordered pairs (h1, h2) sharing the rank-r taxon, minus pairs within one group"""
        tot = 0
        for (t, h), n in cnt[r].items():
            if h == h1:
                m = cnt[r].get((t, h2), 0)
                tot += n * (n - 1) / 2 if h1 == h2 else n * m
        for (t, h, g), n in cntg[r].items():
            if h == h1:
                m = cntg[r].get((t, h2, g), 0)
                tot -= n * (n - 1) / 2 if h1 == h2 else n * m
        return tot

    xh = collections.Counter()
    for a_, b_, d_ in zip(H[XI].tolist(), H[XJ].tolist(), Xd.tolist()):
        xh[(min(a_, b_), max(a_, b_), d_)] += 1

    def habitat_matrix(dmin):
        """observed exchanging pairs of each habitat pair over those expected from the stratum rates (strata >= dmin)"""
        strata = list(range(dmin, 7))
        rate = {d: (num[(d, True)] + num[(d, False)]) / (den[(d, True)] + den[(d, False)]) for d in strata}
        obs = collections.Counter()
        for i, j, d, same, ha, hb in exl:
            if d >= dmin:
                obs[tuple(sorted((ha, hb)))] += 1
        mat = {}
        for a_i, h1 in enumerate(HABS):
            for h2 in HABS[a_i:]:
                allp = pairs_between(h1, h2, -1)
                sh = {r: pairs_between(h1, h2, r) for r in range(6)}
                e = 0.0
                k1, k2 = sorted((HABS.index(h1), HABS.index(h2)))
                for d in strata:
                    r = 5 - d
                    n_ = allp - sh[0] if d == 6 else sh[r] - sh[r + 1]
                    n_ -= xh.get((k1, k2, d), 0)
                    e += n_ * rate[d]
                o = obs[tuple(sorted((h1, h2)))]
                mat[f"{h1}|{h2}"] = {"observed": o, "expected": e, "ratio": o / e if e else None}
        return mat
    res["habitat_pairs_cross_genus"] = habitat_matrix(1)
    res["habitat_pairs_cross_family"] = habitat_matrix(2)
    res["habitat_pairs_cross_order"] = habitat_matrix(3)
    name = ("ecology" + ("" if partners in ("complete", "all") else "_" + partners) + ("_known" if known else "")
            + ("_xstudy" if xstudy else "") + ("_xsample" if xsample else ""))
    res["xstudy"], res["xsample"], res["known_source_only"] = xstudy, xsample, known
    out = os.path.join(ROOT, "results_hgt", tag, name + ".json")
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if not k.startswith("habitat_pairs")}, indent=1))


if __name__ == "__main__":
    main()
