#!/usr/bin/env python3
"""Contig-level contamination screening of simulated bins, from the GTDB classification runs.

Every contig was classified on its own (gtdb_bench.py search), so bins can be assembled after the
fact. For each of the 294 query genomes (host) and each contamination level (5, 10, 20% of host
bases), a bin is the host's contigs plus contigs of one other random query genome (contaminant),
added in random order until the level is reached, skipping contigs that would overshoot it by more than
25% (seed 7); hosts whose contaminant has no small-enough contig are skipped. The bin's taxon is the base-weighted
majority of its contigs' top-hit taxa; a contig is flagged when its top hit has another taxon
(species or genus rule). Contigs without a hit are not flagged ("unknown"); the "strict" rule also
flags them. Precision and recall of contaminant bases and contigs, overall and by contig length.
Writes data/gtdb/results/contam.json.
"""
import csv
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gtdb_bench as gb  # noqa: E402

G = gb.G
LEVELS = [0.05, 0.10, 0.20]
LB = [(500, 2000, "0.5-2k"), (2000, 10000, "2-10k"), (10000, 10 ** 12, ">=10k")]


def main():
    TX = gb.taxonomy()
    meta = list(csv.DictReader(open(os.path.join(G, "meta", f"query_contigs{gb.SET}.tsv")), delimiter="\t"))
    by_g = {}
    for m in meta:
        by_g.setdefault(m["gid"], []).append((m["contig"], int(m["length"])))
    gids = sorted(by_g)
    tools = {"aif": gb.load_aif(os.path.join(gb.RES, "aif_0.88.tsv")), "skani": gb.load_skani(os.path.join(gb.RES, "skani.tsv"))}
    rng = random.Random(7)
    bins = []
    for lvl in LEVELS:
        for h in gids:
            c = rng.choice([g for g in gids if g != h])
            hb = sum(L for _, L in by_g[h])
            cont = by_g[c][:]
            rng.shuffle(cont)
            add, tot = [], 0
            for x in cont:  # contigs that keep the bin within 1.25x the target level (no overshoot by one huge contig)
                if tot >= lvl * hb:
                    break
                if tot + x[1] > 1.25 * lvl * hb:
                    continue
                add.append(x)
                tot += x[1]
            if not add:
                continue
            bins.append((lvl, h, c, by_g[h], add))
    res = {}
    for tname, H in tools.items():
        top = {}
        for c, hs in H.items():
            if hs:
                sp, gen, _ = TX.get(hs[0][0], ("?", "?", "?"))
                top[c] = (sp, gen)
        for rank in ("species", "genus"):
            ri = 0 if rank == "species" else 1
            for rule in ("default", "strict"):
                for lvl in LEVELS + ["all"]:
                    st = {"tp_b": 0, "fp_b": 0, "fn_b": 0, "tp_c": 0, "fp_c": 0, "fn_c": 0}
                    bylen = {lab: {"tp": 0, "fn": 0, "fp": 0} for _, _, lab in LB}
                    for lv, h, cg, host, add in bins:
                        if lvl != "all" and lv != lvl:
                            continue
                        votes = {}
                        for c, L in host + add:
                            if c in top:
                                votes[top[c][ri]] = votes.get(top[c][ri], 0) + L
                        maj = max(votes, key=votes.get) if votes else None
                        for c, L, is_cont in [(c, L, 0) for c, L in host] + [(c, L, 1) for c, L in add]:
                            flag = (c in top and top[c][ri] != maj) or (rule == "strict" and c not in top)
                            lab = next(lb for lo, hi, lb in LB if lo <= L < hi)
                            if flag and is_cont:
                                st["tp_b"] += L
                                st["tp_c"] += 1
                                bylen[lab]["tp"] += 1
                            elif flag:
                                st["fp_b"] += L
                                st["fp_c"] += 1
                                bylen[lab]["fp"] += 1
                            elif is_cont:
                                st["fn_b"] += L
                                st["fn_c"] += 1
                                bylen[lab]["fn"] += 1
                    pr = lambda tp, fp: tp / max(1, tp + fp)  # noqa: E731
                    out = {"precision_bases": pr(st["tp_b"], st["fp_b"]), "recall_bases": st["tp_b"] / max(1, st["tp_b"] + st["fn_b"]),
                           "precision_contigs": pr(st["tp_c"], st["fp_c"]), "recall_contigs": st["tp_c"] / max(1, st["tp_c"] + st["fn_c"]),
                           "recall_by_length": {k: v["tp"] / max(1, v["tp"] + v["fn"]) for k, v in bylen.items()},
                           "precision_by_length": {k: v["tp"] / max(1, v["tp"] + v["fp"]) for k, v in bylen.items()}}
                    res[f"{tname}|{rank}|{rule}|{lvl}"] = out
    res["n_bins"] = {str(l): sum(1 for b in bins if b[0] == l) for l in LEVELS}
    json.dump(res, open(os.path.join(gb.RES, "contam.json"), "w"), indent=1)
    print("bins per level:", res["n_bins"])
    for k, v in res.items():
        if k.endswith("|all"):
            print(f"{k:32s} bases P={v['precision_bases']:.3f} R={v['recall_bases']:.3f} | contigs P={v['precision_contigs']:.3f} "
                  f"R={v['recall_contigs']:.3f} | R by length {({a: round(b, 3) for a, b in v['recall_by_length'].items()})}")


if __name__ == "__main__":
    main()
