#!/usr/bin/env python3
"""Contig-level screen of the UHGG v3.0 human-gut catalogue (286,621 genomes) against GTDB R232 with aif.

  bench_real/uhgg_download.sh                  # genomes (Prokka GFF3 with ##FASTA) -> data/uhgg3/gff
  python3 -I bench_real/uhgg_screen.py prepare # assemblies -> NPART FASTA parts (genomes kept whole), contig table
  python3 -I bench_real/uhgg_screen.py search  # aif db search (GTDB R232, budget 0.88, top 5 hits) per part, timed
  python3 -I bench_real/uhgg_screen.py genomes # per-genome foreign-sequence calls (rule below) and flagged contigs

Rule (fixed on the 3,000-genome pilot and the simulated GTDB bins before the full screen):
  self reference  a GTDB representative that is the top hit at ANI >= 0.999 for >= 50% of a genome's bases (the
                  genome, or a near-clone, is itself in GTDB R232); its hits are set aside
  host genus      the self reference's genus; otherwise the base-weighted majority genus of top hits at ANI >= 0.95
                  (weak: of all top hits, when species-level hits cover < 20% of the genome)
  host-equivalent a reference whose hits at ANI >= 0.95 cover >= 50% of the genome's bases (GTDB's own species
  references      criterion, ANI >= 95% and AF >= 0.5): conspecific with the genome whatever its genus label. Added
                  after the held-out screen: GTDB R232 files some E. coli assemblies under placeholder genera
                  (ECMA0423, G047199095, G047199115) and some chimeric Enterobacter/Escherichia assemblies under the
                  wrong genus, which made isolates of these genera look contaminated (results_uhgg/prespecified/)
  foreign contig  best remaining hit at ANI >= 0.95 over AF >= 0.5 to another genus, and no hit at ANI >= 0.95 to the
                  host genus or to a host-equivalent reference among the top 5 (strict variant: AF >= 0.8): sequence
                  shared at species-level identity with another genus, by contamination or by recent horizontal transfer
  class           mobile: a mobile-hallmark CDS (Prokka product: transposase, integrase, recombinase, conjugation,
                  relaxase, phage, plasmid ...; MOBKW), or MGnify's mobilome (>= 50% of CDS tagged mobile_element_type,
                  or >= 50% of length in plasmid/prophage regions; species representatives only);
                  core: otherwise, carries a core-gene CDS (ribosomal proteins, RNA polymerase, EF-Tu/Ts/G, IF-1/2/3,
                  gyrase, RecA, Ffh, SecY, aminoacyl-tRNA ligases; CORE); accessory: other coding contigs; non-coding
  chromosomal     core + accessory (the first pilot's 'chromosomal'); foreign = chromosomal bases / remaining bases
  contamination   core-anchored: the non-mobile foreign contigs of a source species (GTDB species of the best hit) in a
                  genome, when they carry >= 1 core-gene CDS (strict: >= 2); contam = anchored bases / remaining bases.
                  On the pilot, foreign contigs in MAGs carried 8x fewer core genes per kb than unflagged contigs of the
                  same length and 70% hypothetical proteins: most cross-genus sequence is accessory (mobile) DNA, not
                  misbinned genome, so per-genome contamination is anchored on core genes

UHGG_PILOT=N restricts prepare/search to the first N downloaded genomes (metadata order) in data/uhgg3/pilot.
"""
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import zlib
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
U = os.environ.get("UHGG_DIR", os.path.join(ROOT, "data", "uhgg3"))
PILOT = int(os.environ.get("UHGG_PILOT", "0"))
W = os.path.join(U, "pilot") if PILOT else U
G = os.path.join(ROOT, "data", "gtdb")
AIF = os.environ.get("AIF", os.path.join(ROOT, "bin", "aif"))
DB = os.path.join(G, "aif", "gtdb_0.88.aifd")
THREADS = os.environ.get("THREADS", "48")
NPART = 8 if PILOT else 128
MINLEN = 500
# ribosomal-protein CDS by Prokka product name (excludes modifying enzymes such as "ribosomal protein L11 methyltransferase")
RIBO = re.compile(r"product=(?:30S|50S) ribosomal protein [LS][0-9A-Za-z]*(?:;|$)")
# MGnify's mobilome (mobile_element_type, plasmid / prophage regions) and AMR annotations exist on species representatives
# only; every genome has Prokka products, so mobile hallmark genes and core genes are also counted from product names
MOBKW = re.compile(r"product=[^;]*(?:transpos|integrase|recombinase|resolvase|excisionase|conjugal|conjugative|relaxase|"
                   r"mobilization|(?<![a-z])phage(?! shock)|terminase|capsid|portal protein|\bholin\b|baseplate|tail fiber|"
                   r"tape measure|plasmid|antirestriction|type iv secretion)", re.I)
CORE = re.compile(r"product=(?:(?:30S|50S) ribosomal protein [LS][0-9A-Za-z]*|DNA-directed RNA polymerase subunit (?:alpha|beta)"
                  r"[^;]*|Elongation factor (?:Tu|Ts|G)|Translation initiation factor IF-[123]|DNA gyrase subunit [AB]|"
                  r"Protein RecA|Signal recognition particle protein|Protein translocase subunit SecY|[A-Za-z()/-]+--tRNA ligase"
                  r"[^;]*)(?:;|$)")
CONTIG_COLS = ["contig", "genome", "length", "gc", "screened", "cds", "mge_cds", "ribo_cds", "plasmid_bp", "prophage_bp",
               "amr_cds", "mobkw_cds", "hypo_cds", "core_cds"]
ANI_SP, AF_MIN, AF_STRICT, SELF_ANI, SELF_COV, SP_COV, EQ_COV = 0.95, 0.5, 0.8, 0.999, 0.5, 0.2, 0.5


def meta():
    return list(csv.DictReader(open(os.path.join(U, "meta", "genomes-all_metadata.tsv")), delimiter="\t"))


def part_of(genome):
    return zlib.crc32(genome.encode()) % NPART


def extract(job):
    """one part: FASTA of its genomes' contigs (>= MINLEN) and a contig table (CONTIG_COLS: GC, MGnify mobilome
    annotations: CDS tagged mobile_element_type, plasmid and prophage regions; ribosomal-protein CDS; AMR CDS; CDS whose
    Prokka product is a mobile hallmark (MOBKW), hypothetical, or a core gene (CORE)); truncated gzip files are reported"""
    p, genomes = job
    fa = os.path.join(W, "fa", f"part{p:03d}.fa")
    bad, n = [], 0
    with open(fa + ".tmp", "w") as fo, open(os.path.join(W, "fa", f"part{p:03d}.contigs.tsv"), "w") as ft:
        for g in genomes:
            ann = {}  # contig -> [CDS, mobile_element_type CDS, ribosomal-protein CDS, plasmid bp, prophage bp, AMR CDS,
            #                     mobile-hallmark CDS, hypothetical CDS, core CDS]
            try:
                with gzip.open(os.path.join(U, "gff", g + ".gff.gz"), "rt") as fh:
                    for line in fh:
                        if line.startswith("##FASTA"):
                            break
                        if line.startswith("#"):
                            continue
                        f = line.rstrip("\n").split("\t")
                        if len(f) < 9:
                            continue
                        a = ann.setdefault(f[0], [0] * 9)
                        if f[2] == "CDS":
                            a[0] += 1
                            a[1] += "mobile_element_type=" in f[8]
                            a[2] += RIBO.search(f[8]) is not None
                            a[5] += "element_type=AMR" in f[8]
                            a[6] += MOBKW.search(f[8]) is not None
                            a[7] += "product=hypothetical protein" in f[8]
                            a[8] += CORE.search(f[8]) is not None
                        elif f[2] in ("plasmid", "prophage"):
                            a[3 if f[2] == "plasmid" else 4] += int(f[4]) - int(f[3]) + 1
                    seqs, name = {}, None
                    for line in fh:
                        if line.startswith(">"):
                            name = line[1:].split()[0]
                            seqs[name] = []
                        else:
                            seqs[name].append(line.strip())
            except (EOFError, OSError, zlib.error) as e:
                bad.append((g, str(e)))
                continue
            for c, parts in seqs.items():
                s = "".join(parts).upper()
                gc = (s.count("G") + s.count("C")) / max(1, len(s) - s.count("N"))
                a = ann.get(c, [0] * 9)
                ft.write(f"{c}\t{g}\t{len(s)}\t{gc:.4f}\t{int(len(s) >= MINLEN)}\t" + "\t".join(map(str, a)) + "\n")
                if len(s) >= MINLEN:
                    fo.write(f">{c}\n{s}\n")
                    n += 1
    os.replace(fa + ".tmp", fa)
    return p, n, bad


def prepare():
    M = meta()
    have = set(f[:-7] for f in os.listdir(os.path.join(U, "gff")) if f.endswith(".gff.gz"))
    gs = [r["Genome"] for r in M if r["Genome"] in have]
    if PILOT:  # the pilot's genome list is fixed by its first prepare (downloads were still running)
        pl = os.path.join(W, "genomes.txt")
        if os.path.exists(pl):
            gs = [x.strip() for x in open(pl)]
        else:
            gs = gs[:PILOT]
            os.makedirs(W, exist_ok=True)
            open(pl, "w").write("\n".join(gs) + "\n")
    missing = len(M) - len(have & {r["Genome"] for r in M})
    os.makedirs(os.path.join(W, "fa"), exist_ok=True)
    jobs = [(p, sorted(g for g in gs if part_of(g) == p)) for p in range(NPART)]
    bad, n = [], 0
    with Pool(int(THREADS)) as pool:
        for p, k, b in pool.imap_unordered(extract, jobs):
            n += k
            bad += b
    json.dump({"genomes": len(gs), "missing_downloads": missing, "contigs_screened": n, "bad_gzip": bad},
              open(os.path.join(W, "prepare.json"), "w"), indent=1)
    print(len(gs), "genomes,", n, "contigs >=", MINLEN, "bp;", len(bad), "truncated files;", missing, "not downloaded")


def search():
    """aif db search per part (database read into the page cache first); wall time and peak memory per part"""
    out = os.path.join(W, "aif")
    os.makedirs(out, exist_ok=True)
    for _ in range(2):  # a second read moves the pages to the active list, out of reach of streaming writes
        subprocess.run(["cat", DB], stdout=subprocess.DEVNULL, check=True)
    tp = os.path.join(W, "search_time.json")
    T = json.load(open(tp)) if os.path.exists(tp) else {}
    for p in range(NPART):
        key = f"part{p:03d}"
        if key in T:
            continue
        fa = os.path.join(W, "fa", key + ".fa")
        tf = tempfile.NamedTemporaryFile(delete=False).name
        with open(os.path.join(out, key + ".tsv.tmp"), "w") as fo, open(os.path.join(out, key + ".log"), "w") as fe:
            subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf, AIF, "db", "search", DB, fa, "-n", "5", "-t", THREADS],
                           stdout=fo, stderr=fe, check=True)
        os.replace(os.path.join(out, key + ".tsv.tmp"), os.path.join(out, key + ".tsv"))
        t, m = open(tf).read().split()[-2:]
        os.unlink(tf)
        T[key] = {"s": float(t), "kb": int(m), "bytes": os.path.getsize(fa)}
        json.dump(T, open(tp, "w"), indent=1)
        print(key, T[key], flush=True)
    tot = sum(v["s"] for v in T.values())
    print(f"total {tot:.0f} s for {sum(v['bytes'] for v in T.values()) / 1e9:.1f} GB of FASTA at {THREADS} threads")


def taxonomy():
    """{GTDB accession: (domain, phylum, class, order, family, genus, species)}"""
    T = {}
    for fn in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                T[r["accession"][3:]] = tuple(x[3:] for x in r["gtdb_taxonomy"].split(";"))
    return T


TX = {}
NOTAX = ("?",) * 7
GENOME_COLS = ["genome", "contigs", "total_bp", "screened_bp", "hit_bp", "species_hit_bp", "self_ref", "host_genus",
               "host_family", "host_src", "foreign_n", "foreign_bp", "mobile_n", "mobile_bp", "core_n", "core_bp",
               "accessory_n", "accessory_bp", "noncoding_bp", "chrom_bp", "chrom_strict_bp", "chrom_xfam_bp", "anchored_bp",
               "anchored_strict_bp", "anchored_species", "anchored_core", "amr_on_foreign", "top_source", "top_anchor",
               "foreign", "contam", "equiv_n", "equiv_xgenus", "equiv_xgenus_refs"]
FLAG_COLS = ["contig", "genome", "length", "gc", "cls", "strict", "ref", "ani", "af", "source_species", "source_genus",
             "source_family", "source_phylum", "host_genus", "host_family", "cds", "mge_cds", "mobkw_cds", "hypo_cds",
             "core_cds", "ribo_cds", "amr_cds", "anchored"]


def screen_part(p):
    key = f"part{p:03d}"
    C = {}
    for line in open(os.path.join(W, "fa", key + ".contigs.tsv")):
        f = line.rstrip("\n").split("\t")
        C[f[0]] = (f[1], int(f[2]), float(f[3]), int(f[4])) + tuple(int(x) for x in f[5:14])
    H = {}
    with open(os.path.join(W, "aif", key + ".tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t", 10)
            if f[1] != "*" and f[9] == "ok":
                H.setdefault(f[0], []).append((f[1], float(f[3]), float(f[4])))
    byg = {}
    for c, v in C.items():
        byg.setdefault(v[0], []).append(c)
    rows, flags = [], []
    for g in sorted(byg):
        cs = byg[g]
        L = {c: C[c][1] for c in cs}
        tot = sum(L.values())
        cov = {}
        for c in cs:
            if c in H and H[c][0][1] >= SELF_ANI:
                cov[H[c][0][0]] = cov.get(H[c][0][0], 0) + L[c]
        sref = max(cov, key=cov.get) if cov else None
        if sref and cov[sref] < SELF_COV * tot:
            sref = None
        if sref:
            host, src = TX.get(sref, NOTAX), "self"
        else:
            v, w = {}, {}
            for c in cs:
                if c in H:
                    t = TX.get(H[c][0][0], NOTAX)
                    w[t[:6]] = w.get(t[:6], 0) + L[c]
                    if H[c][0][1] >= ANI_SP:
                        v[t[:6]] = v.get(t[:6], 0) + L[c]
            if v and sum(v.values()) >= SP_COV * tot:
                host, src = max(v, key=v.get), "species"
            elif w:
                host, src = max(w, key=w.get), "weak"
            else:
                host, src = None, "none"
        hg = host[5] if host else None
        eq = {}  # host-equivalent references: hits at ANI >= 0.95 cover >= EQ_COV of the genome (any genus label)
        for c in cs:
            for h in H.get(c, []):
                if h[1] >= ANI_SP and h[0] != sref:
                    eq[h[0]] = eq.get(h[0], 0) + L[c] * h[2]
        E = {r for r, v in eq.items() if v >= EQ_COV * tot}
        exg = sorted(r for r in E if TX.get(r, NOTAX)[5] != hg)
        acc = dict(foreign_n=0, foreign_bp=0, mobile_n=0, mobile_bp=0, core_n=0, core_bp=0, accessory_n=0, accessory_bp=0,
                   noncoding_bp=0, chrom_bp=0, chrom_strict_bp=0, chrom_xfam_bp=0, amr=0)
        srcbp, grp, fl = {}, {}, []
        hit_bp = sum(L[c] for c in cs if c in H)
        sp_bp = sum(L[c] for c in cs if c in H and H[c][0][1] >= ANI_SP)
        for c in cs:
            hs = [h for h in H.get(c, []) if h[0] != sref]
            if not hs or hg is None:
                continue
            t = hs[0]
            tt = TX.get(t[0], NOTAX)
            if not (t[1] >= ANI_SP and t[2] >= AF_MIN and tt[5] != hg):
                continue
            if any(h[1] >= ANI_SP and (h[0] in E or TX.get(h[0], NOTAX)[5] == hg) for h in hs):
                continue
            _, ln, gc, _, cds, mge, ribo, pl, pp, amr, mobkw, hypo, core = C[c]
            if (cds and mge / cds >= 0.5) or pl + pp >= 0.5 * ln or mobkw:
                cls = "mobile"
            elif cds == 0:
                cls = "noncoding"
            else:
                cls = "core" if core else "accessory"
            strict = int(t[2] >= AF_STRICT)
            acc["foreign_n"] += 1
            acc["foreign_bp"] += ln
            if cls == "noncoding":
                acc["noncoding_bp"] += ln
            else:
                acc[cls + "_n"] += 1
                acc[cls + "_bp"] += ln
            if cls in ("core", "accessory"):
                acc["chrom_bp"] += ln
                acc["chrom_strict_bp"] += ln * strict
                acc["chrom_xfam_bp"] += ln * (tt[4] != host[4])
                srcbp[tt[5]] = srcbp.get(tt[5], 0) + ln
                gs_ = grp.setdefault(tt[6], [0, 0])
                gs_[0] += ln
                gs_[1] += core
            acc["amr"] += amr
            fl.append([c, g, ln, f"{gc:.4f}", cls, strict, t[0], f"{t[1]:.4f}", f"{t[2]:.3f}", tt[6], tt[5], tt[4], tt[1], hg,
                       host[4], cds, mge, mobkw, hypo, core, ribo, amr, tt[6]])
        anch = {k: v for k, v in grp.items() if v[1] >= 1}
        for r in fl:  # anchored: a core or accessory contig of a source species whose contigs carry core genes
            r[-1] = int(r[4] in ("core", "accessory") and r[-1] in anch)
        flags += fl
        abp = sum(v[0] for v in anch.values())
        sbp = sum(v[0] for v in anch.values() if v[1] >= 2)
        top = max(srcbp, key=srcbp.get) if srcbp else ""
        topa = max(anch, key=lambda k: anch[k][0]) if anch else ""
        rows.append([g, len(cs), tot, sum(L[c] for c in cs if C[c][3]), hit_bp, sp_bp, sref or "", hg or "", host[4] if host else "",
                     src, acc["foreign_n"], acc["foreign_bp"], acc["mobile_n"], acc["mobile_bp"], acc["core_n"], acc["core_bp"],
                     acc["accessory_n"], acc["accessory_bp"], acc["noncoding_bp"], acc["chrom_bp"], acc["chrom_strict_bp"],
                     acc["chrom_xfam_bp"], abp, sbp, len(anch), sum(v[1] for v in anch.values()), acc["amr"], top, topa,
                     f"{acc['chrom_bp'] / max(1, tot - acc['chrom_bp']):.5f}", f"{abp / max(1, tot - abp):.5f}", len(E), len(exg),
                     ";".join(exg)])
    out = os.path.join(W, "screen")
    with open(os.path.join(out, key + ".genomes.tsv"), "w") as fh:
        fh.write("\t".join(GENOME_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out, key + ".flagged.tsv"), "w") as fh:
        fh.write("\t".join(FLAG_COLS) + "\n")
        for r in flags:
            fh.write("\t".join(map(str, r)) + "\n")
    return p, len(rows), len(flags)


def genomes():
    TX.update(taxonomy())  # inherited by the forked workers
    os.makedirs(os.path.join(W, "screen"), exist_ok=True)
    parts = [p for p in range(NPART) if os.path.exists(os.path.join(W, "aif", f"part{p:03d}.tsv"))]
    import multiprocessing  # fork, so the workers share TX (Python >= 3.14 defaults to forkserver)
    with multiprocessing.get_context("fork").Pool(min(len(parts), int(THREADS))) as pool:
        n = [r for r in pool.imap_unordered(screen_part, parts)]
    print(sum(x[1] for x in n), "genomes,", sum(x[2] for x in n), "flagged contigs in", len(parts), "parts")


if __name__ == "__main__":
    {"prepare": prepare, "search": search, "genomes": genomes}[sys.argv[1]]()
