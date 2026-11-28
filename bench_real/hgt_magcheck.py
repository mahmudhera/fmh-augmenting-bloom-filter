#!/usr/bin/env python3
"""Is DNA shared by MAGs transfer or assembly artifact? Gene content of shared blocks by the kind of genome pair.

Transferred DNA is rich in mobile-element genes (transposases, integrases, conjugation, phage) and poor in
housekeeping genes; a fragment of another genome put into a MAG by chimeric assembly or mis-binning carries that
genome's genes at random, so its content approaches the background. Genes of the shared regions
(hgt_annotate.py, data/hgt/annot_TAG/) are assigned to the cross-genus blocks they lie in, and each block to the kind
of its genome pair:
  isolates         both genomes isolates
  MAG, same study  at least one MAG, and the two share a BioProject or a source sample (hgt_provenance.py)
  MAG, other study, same / different habitat
Background genes come from random windows of the same genomes, split by the carrier's category (isolate / MAG).

Pair-level content filter: results_hgt/TAG_mobile/pairs.tsv.gz keeps the cross-genus pairs whose shared regions carry
at least one mobile-element gene (either direction), for hgt_ecology.py; the same filter applies to isolates and MAGs.

Also: the share of exchanging MAG pairs (both with known sources) that have a common source sample or BioProject, by
distance ("mag_pairs_provenance").

Usage: /usr/bin/python3 -I bench_real/hgt_magcheck.py TAG    -> results_hgt/TAG/magcheck.json, results_hgt/TAG_mobile/
"""
import collections
import csv
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hgt_annotate import CLASS_RULES, DORDER, HD, ROOT, classify, kegg  # noqa: E402

MOBILE = {"transposase", "integrase/recombinase", "conjugation / T4SS", "phage"}
HABS = {"gut", "host_other", "food", "engineered", "marine", "plant", "extreme", "subsurface", "freshwater", "soil",
        "sediment"}


def main():
    tag = sys.argv[1]
    res = os.path.join(ROOT, "results_hgt", tag)
    wd = os.path.join(HD, "annot_" + tag)
    G = {r["accession"]: r for r in csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t")}
    bp = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(HD, "bioproject.tsv")))
    SRC = {}
    for l in open(os.path.join(HD, "provenance", "sources.tsv")):
        a, s = l.rstrip("\n").split("\t")
        SRC[a] = set(s.split(","))

    def kind(q, p):
        iq, ip = G[q]["category"] == "isolate", G[p]["category"] == "isolate"
        if iq and ip:
            return "isolates"
        same_bp = bp.get(q) == bp.get(p) and bp.get(q) not in (None, "", "none", "PRJNA224116")
        if same_bp or (SRC.get(q, set()) & SRC.get(p, set())):
            return "MAG, same study"
        hq, hp = G[q]["habitat"], G[p]["habitat"]
        if hq in HABS and hp in HABS:
            return "MAG, other study, " + ("same habitat" if hq == hp else "different habitat")
        return "MAG, other study, habitat unknown"
    blocks = collections.defaultdict(list)
    for b in csv.DictReader(gzip.open(os.path.join(res, "blocks.tsv.gz"), "rt"), delimiter="\t"):
        if b["distance"] != "same genus":
            blocks[(b["query"], int(b["contig"]))].append((int(b["start"]), int(b["end"]), b["partner"],
                                                           DORDER.index(b["distance"])))
    ko = {}
    for l in open(os.path.join(wd, "ko.tsv")):
        q, s = l.split("\t")[:2]
        if q not in ko and "~" in s:
            ko[q] = s.split("~")[1]
    names, brite = kegg()
    # counts[(dmin, kind)][class]; background by carrier category
    cnt = collections.defaultdict(collections.Counter)
    bg = collections.defaultdict(collections.Counter)
    pair_content = collections.defaultdict(collections.Counter)
    cnt_h = collections.defaultdict(collections.Counter)
    for l in open(os.path.join(wd, "proteins.faa")):
        if l[0] != ">":
            continue
        f = l[1:].split(" # ")
        gid = f[0]
        p = gid.rsplit("_", 1)[0].split("|")
        c = classify(ko.get(gid), names, brite)
        acc = p[1]
        if p[0] == "B":
            bg["isolate" if G[acc]["category"] == "isolate" else "MAG"][c] += 1
            continue
        ci, s0 = int(p[2]), int(p[3])
        mid = s0 + (int(f[1]) + int(f[2])) // 2
        seen = set()
        for st, en, partner, d in blocks.get((acc, ci), ()):
            if st <= mid < en:
                pk = (min(acc, partner), max(acc, partner))
                pc = pair_content[pk]
                pc["genes"] += 1
                pc["mobile"] += c in MOBILE
                pc["ribosome"] += c == "ribosome / translation"
                k = kind(acc, partner)
                for dmin in (1, 2, 5):
                    if d >= dmin and (dmin, k) not in seen:
                        seen.add((dmin, k))
                        cnt[(dmin, k)][c] += 1
                # same-habitat pairs between families, by habitat: isolates against MAGs of other studies
                h = G[acc]["habitat"]
                if d >= 2 and h == G[partner]["habitat"] and h in HABS and k in ("isolates",
                                                                                  "MAG, other study, same habitat"):
                    if (h, k) not in seen:
                        seen.add((h, k))
                        cnt_h[(h, k)][c] += 1
    cls = [c for c, _, _ in CLASS_RULES] + ["other KEGG ortholog", "no KEGG ortholog"]

    def summ(C):
        n = sum(C.values())
        mob = sum(C[c] for c in MOBILE)
        return {"genes": n, "mobile": mob / n if n else None,
                "ribosome": C["ribosome / translation"] / n if n else None,
                "resistance": (C["antimicrobial resistance"] + C["metal resistance"]) / n if n else None,
                "no_ko": C["no KEGG ortholog"] / n if n else None, "classes": {c: C[c] for c in cls}}
    out = {"background": {k: summ(v) for k, v in bg.items()}, "shared": {}}
    lab = {1: "between genera", 2: "between families", 5: "between phyla"}
    for (dmin, k), C in sorted(cnt.items()):
        out["shared"].setdefault(lab[dmin], {})[k] = summ(C)
    out["same_habitat_between_families"] = {}
    for (h, k), C in sorted(cnt_h.items()):
        out["same_habitat_between_families"].setdefault(h, {})["isolates" if k == "isolates" else "MAGs"] = summ(C)
    # pair-level support: share of cross-genus pairs (by kind) with a mobile gene, with a ribosomal gene
    ps = collections.defaultdict(collections.Counter)
    for (a, b), pc in pair_content.items():
        k = kind(a, b)
        ps[k]["pairs"] += 1
        ps[k]["with_mobile"] += pc["mobile"] > 0
        ps[k]["with_ribosome"] += pc["ribosome"] > 0
    out["pairs"] = {k: dict(v) for k, v in ps.items()}
    od = os.path.join(ROOT, "results_hgt", tag + "_mobile")
    os.makedirs(od, exist_ok=True)
    nk = 0
    opath = os.path.join(od, "pairs.tsv.gz")
    with gzip.open(opath + ".tmp", "wt") as fo, gzip.open(os.path.join(res, "pairs.tsv.gz"), "rt") as fi:
        fo.write(fi.readline())
        for l in fi:
            x = l.split("\t", 3)
            if x[2] == "same genus":
                continue
            pc = pair_content.get((min(x[0], x[1]), max(x[0], x[1])))
            if pc and pc["mobile"] > 0:
                fo.write(l)
                nk += 1
    os.replace(opath + ".tmp", opath)
    out["mobile_filter_rows"] = nk
    # exchanging pairs of MAGs with known sources, by distance: share with a common source sample / BioProject
    cs = collections.defaultdict(collections.Counter)
    seen = set()
    for r in csv.DictReader(gzip.open(os.path.join(res, "pairs.tsv.gz"), "rt"), delimiter="\t"):
        q, p = r["query"], r["partner"]
        k = (min(q, p), max(q, p))
        if k in seen or q not in SRC or p not in SRC:
            continue
        seen.add(k)
        c = cs[r["distance"]]
        a = bool(SRC[q] & SRC[p])
        b = bp.get(q) == bp.get(p) and bp.get(q) not in (None, "", "none", "PRJNA224116")
        c["pairs"] += 1
        c["cosample"] += a
        c["same_project"] += b
        c["cosample_other_project"] += a and not b
    out["mag_pairs_provenance"] = {d: dict(cs[d]) for d in DORDER}
    json.dump(out, open(os.path.join(res, "magcheck.json"), "w"), indent=1)
    for h, D in out["same_habitat_between_families"].items():
        print(f"same habitat {h:11s} " + "  ".join(
            f"{k}: genes {v['genes']:6d} mobile {100 * v['mobile']:4.1f}% ribosome {100 * v['ribosome']:4.2f}%"
            for k, v in sorted(D.items())))
    for k, v in sorted(ps.items()):
        print(f"pairs {k:36s} {v['pairs']:7d}  with mobile gene {100 * v['with_mobile'] / v['pairs']:5.1f}%  "
              f"with ribosomal gene {100 * v['with_ribosome'] / v['pairs']:5.1f}%")
    for k, v in out["background"].items():
        print(f"background {k:8s} genes {v['genes']:8d}  mobile {100 * v['mobile']:5.1f}%  ribosome "
              f"{100 * v['ribosome']:5.2f}%  resistance {100 * v['resistance']:5.2f}%  no KO {100 * v['no_ko']:5.1f}%")
    for L, D in out["shared"].items():
        print(L)
        for k, v in D.items():
            print(f"  {k:36s} genes {v['genes']:8d}  mobile {100 * v['mobile']:5.1f}%  ribosome "
                  f"{100 * v['ribosome']:5.2f}%  resistance {100 * v['resistance']:5.2f}%  no KO {100 * v['no_ko']:5.1f}%")


if __name__ == "__main__":
    main()
