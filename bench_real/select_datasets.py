#!/usr/bin/env python3
"""Choose the real genome sets for the real-data benchmark (deterministic, seed 20261121).

Inputs are the skani-test metadata files (github.com/bluenote-1577/skani-test) copied
to data/real/meta: D2_anim.txt (fastANI/skani D2: B. cereus group vs B. anthracis 52-G),
d3_l.txt (fastANI/skani D3: E. coli assemblies), 2328_l.txt (Pasolli et al. species bin
2328, 195 MAGs).

Writes data/real/lists/<set>.tsv with columns: set, genome, source, accession_or_path,
stratum, reference.

Sets
  bcereus   80 draft assemblies from D2, stratified by published ANIm: 18 each in
            85-90, 90-95, 95-98, >=99.5 and all 8 in 98-99.5; vs B. anthracis 52-G chromosome
  ecoli     48 draft assemblies from D3 vs E. coli K-12 MG1655
  mag2328   48 MAGs of bin 2328 vs the most contiguous MAG of the bin (chosen later)
  paer      24 P. aeruginosa contig-level assemblies (66% GC) vs PAO1
  saur      24 S. aureus contig-level assemblies (33% GC) vs NCTC 8325
  negative  8 E. coli drafts vs B. anthracis and 8 B. cereus-group drafts vs E. coli
"""
import json
import os
import random
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "data", "real", "meta")
OUT = os.path.join(ROOT, "data", "real", "lists")
DATASETS = os.environ.get("DATASETS", "/home/mrahmanh/miniconda3/envs/meta-sim/bin/datasets")
SEED = 20261121


def acc_of(name):
    m = re.search(r"(GC[AF]_\d+)", name)
    return m.group(1) if m else None


def bin_of(a):
    if a < 90:
        return "85-90"
    if a < 95:
        return "90-95"
    if a < 98:
        return "95-98"
    if a < 99.5:
        return "98-99.5"
    return ">=99.5"


def ncbi_drafts(taxon, n, rng):
    """n random GenBank contig-level assemblies of a species (from the first 400 listed)."""
    cmd = [DATASETS, "summary", "genome", "taxon", taxon, "--assembly-level", "contig",
           "--assembly-source", "GenBank", "--limit", "400", "--as-json-lines"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    recs = []
    for line in out.splitlines():
        if not line.startswith("{"):
            continue
        d = json.loads(line)
        st = d.get("assembly_stats", {})
        nc = int(st.get("number_of_contigs", 0) or 0)
        if d.get("organism", {}).get("organism_name", "").startswith(taxon) and 20 <= nc <= 2000:
            recs.append((d["accession"], nc))
    recs.sort()
    return rng.sample(recs, n)


def main():
    rng = random.Random(SEED)
    os.makedirs(OUT, exist_ok=True)
    rows = {}
    # --- bcereus: stratified by the published ANIm value (D2_anim.txt column 3)
    by = {}
    for line in open(os.path.join(META, "D2_anim.txt")):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        a = acc_of(p[1])
        if a:
            by.setdefault(bin_of(float(p[2])), []).append((p[1], a, float(p[2])))
    bc = []
    quota = {"85-90": 18, "90-95": 18, "95-98": 18, "98-99.5": 8, ">=99.5": 18}  # 98-99.5 has only 8 drafts
    for b, q in quota.items():
        for name, a, ani in sorted(rng.sample(sorted(by[b]), q)):
            bc.append(("bcereus", name.replace(".LargeContigs.fna", ""), "ncbi", a, f"{b}|{ani:.2f}", "BA52G"))
    rows["bcereus"] = bc
    # --- ecoli: random D3 drafts
    d3 = sorted(set(l.strip().split("/")[-1] for l in open(os.path.join(META, "d3_l.txt")) if acc_of(l)))
    ec_pick = rng.sample(d3, 56)
    rows["ecoli"] = [("ecoli", n.replace(".LargeContigs.fna", ""), "ncbi", acc_of(n), "D3", "MG1655")
                     for n in sorted(ec_pick[:48])]
    # --- mag2328: random MAGs (the reference MAG is chosen after download, by N50)
    mags = sorted(l.strip() for l in open(os.path.join(META, "2328_l.txt")) if l.strip())
    rows["mag2328"] = [("mag2328", os.path.basename(m).replace(".fa", ""), "skani-test", m, "SGB2328", "REFMAG")
                       for m in mags]  # all 195; prepare.py picks the reference and 48 queries
    # --- composition extremes
    rows["paer"] = [("paer", a, "ncbi", a, f"contigs={nc}", "PAO1") for a, nc in ncbi_drafts("Pseudomonas aeruginosa", 24, rng)]
    rows["saur"] = [("saur", a, "ncbi", a, f"contigs={nc}", "NCTC8325") for a, nc in ncbi_drafts("Staphylococcus aureus", 24, rng)]
    # --- negative controls: other-genus drafts (E. coli not in the ecoli set; B. cereus group not in bcereus)
    used = set(r[3] for r in bc)
    rest_bc = sorted(set((name, a) for lst in by.values() for (name, a, _) in lst if a not in used))
    neg = [("negative", n.replace(".LargeContigs.fna", ""), "ncbi", acc_of(n), "Ecoli_vs_BA52G", "BA52G")
           for n in sorted(ec_pick[48:56])]
    neg += [("negative", n.replace(".LargeContigs.fna", ""), "ncbi", a, "Bcereus_vs_MG1655", "MG1655")
            for n, a in sorted(rng.sample(rest_bc, 8))]
    rows["negative"] = neg
    for s, rr in rows.items():
        with open(os.path.join(OUT, s + ".tsv"), "w") as fh:
            fh.write("set\tgenome\tsource\taccession\tstratum\treference\n")
            for r in rr:
                fh.write("\t".join(r) + "\n")
        print(s, len(rr), file=sys.stderr)


if __name__ == "__main__":
    main()
