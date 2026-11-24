#!/usr/bin/env python3
"""Summaries and validation of the UHGG v3.0 contig screen (bench_real/uhgg_screen.py); system python (numpy, scipy).

  /usr/bin/python3 -I bench_real/uhgg_summary.py stats   # catalogue-wide results -> results_uhgg/summary.json + tables
  /usr/bin/python3 -I bench_real/uhgg_summary.py kraken  # Kraken2 (GTDB R232 reps, the build of kraken2_gtdb.py) on
                                                         # flagged contigs and length-matched controls from the same genomes
UHGG_PILOT=1 reads the pilot (data/uhgg3/pilot) and writes results_uhgg/pilot.

What the screen measures: contigs that share species-level identity (ANI >= 95%, AF >= 50%) with another genus
('foreign'). The analyses separate the two sources of such sequence:
  core-gene deficit  a misbinned fragment of another genome carries core genes at the genome-wide rate; mobile and
                     accessory DNA carries almost none. Comparing core genes per kb of flagged non-mobile contigs with
                     unflagged contigs of the same length class (MAGs), against flagged contigs of isolates (no binning)
                     as the transferred-DNA reference, estimates the contamination share of flagged bases (bootstrap
                     over the 128 parts)
  core anchoring     per genome, contamination is called only for source species whose contigs carry a core gene;
                     its sensitivity at 0.5-10% contamination is simulated from unflagged contigs of high-quality MAGs
  isolates           pure cultures as the control, overall and within species (>= 3 isolates and >= 3 MAGs)
  same sample        is the source genus also binned from the same sample (MGnify's GTDB labels)? Expected rate from
                     the samples of the same study with the same number of MAGs that hold a MAG of the host genus
  GC                 flagged contigs vs length-matched unflagged contigs of the same genomes
  recurrence         pilot only: is a flagged contig present in other genomes of the same species (other studies)?
  Kraken2            genus of each flagged contig and its control (minimizer LCA; another algorithm, same references)
"""
import collections
import csv
import glob
import gzip
import json
import math
import multiprocessing
import os
import random
import subprocess
import sys
import zlib

import numpy as np
from scipy import stats as sst

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
U = os.path.join(ROOT, "data", "uhgg3")
PILOT = bool(os.environ.get("UHGG_PILOT"))
W = os.path.join(U, "pilot") if PILOT else U
OUT = os.path.join(ROOT, "results_uhgg", "pilot") if PILOT else os.path.join(ROOT, "results_uhgg")
G = os.path.join(ROOT, "data", "gtdb")
K2 = os.path.join(G, "kraken2_db")
THREADS = os.environ.get("THREADS", "48")
NPART = 8 if PILOT else 128
GC_MIN_LEN = 5000
LBINS = [(500, 2000), (2000, 5000), (5000, 20000), (20000, 10 ** 9)]
EUK = ["MGYG000906601", "MGYG000840688", "MGYG000826516", "MGYG000825234", "MGYG000903611"]  # known_issues.txt
GINT = ["contigs", "total_bp", "screened_bp", "hit_bp", "species_hit_bp", "foreign_n", "foreign_bp", "mobile_n", "mobile_bp",
        "core_n", "core_bp", "accessory_n", "accessory_bp", "noncoding_bp", "chrom_bp", "chrom_strict_bp", "chrom_xfam_bp",
        "anchored_bp", "anchored_strict_bp", "anchored_species", "anchored_core", "amr_on_foreign", "equiv_n", "equiv_xgenus"]
FINT = ["length", "strict", "cds", "mge_cds", "mobkw_cds", "hypo_cds", "core_cds", "ribo_cds", "amr_cds", "anchored"]


def _np(o):
    """JSON fallback for numpy scalars and arrays."""
    if hasattr(o, "tolist"):
        return o.tolist()
    raise TypeError(type(o).__name__)


def wilson(k, n, z=1.96):
    if n == 0:
        return [0.0, 0.0, 0.0]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [p, max(0.0, c - h), min(1.0, c + h)]


def genus_of(lineage):
    x = lineage.split(";")
    return x[5][3:] if len(x) > 5 else ""


def species_of(lineage):
    x = lineage.split(";")
    return x[6][3:] if len(x) > 6 else ""


def lbin(L):
    return next(i for i, (lo, hi) in enumerate(LBINS) if lo <= L < hi)


def gtdb_tax():
    """{accession: (phylum, family, genus, species)}, {genus: (phylum, family)}, {accession: GC fraction}"""
    T, GC = {}, {}
    for fn in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            rd = csv.reader(fh, delimiter="\t")
            h = next(rd)
            ia, it, ig, ir = h.index("accession"), h.index("gtdb_taxonomy"), h.index("gc_percentage"), h.index("gtdb_representative")
            for r in rd:
                if r[ir] != "t":
                    continue
                x = [y[3:] for y in r[it].split(";")]
                T[r[ia][3:]] = (x[1], x[4], x[5], x[6])
                GC[r[ia][3:]] = float(r[ig]) / 100
    GEN = {t[2]: (t[0], t[1]) for t in T.values()}
    return T, GEN, GC


def ncbi_meta(accs):
    """{accession: (NCBI species, CheckM2 completeness, CheckM2 contamination)} for the given GTDB accessions"""
    out = {}
    for fn in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(os.path.join(G, fn), "rt") as fh:
            rd = csv.reader(fh, delimiter="\t")
            h = next(rd)
            ia, inc, icc, ict = h.index("accession"), h.index("ncbi_taxonomy"), h.index("checkm2_completeness"), h.index("checkm2_contamination")
            for r in rd:
                a = r[ia][3:]
                if a in accs:
                    out[a] = (r[inc].split(";")[-1][3:], r[icc], r[ict])
    return out


def load():
    M = {r["Genome"]: r for r in csv.DictReader(open(os.path.join(U, "meta", "genomes-all_metadata.tsv")), delimiter="\t")}
    R = {}
    for f in sorted(glob.glob(os.path.join(W, "screen", "part*.genomes.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            m = M[r["genome"]]
            for k in GINT:
                r[k] = int(r[k])
            t = r["total_bp"]
            r["foreign"], r["contam"] = float(r["foreign"]), float(r["contam"])
            r["strict"] = r["chrom_strict_bp"] / max(1, t - r["chrom_strict_bp"])
            r["contam_strict"] = r["anchored_strict_bp"] / max(1, t - r["anchored_strict_bp"])
            r["type"] = m["Genome_type"]
            r["rep"] = m["Genome"] == m["Species_rep"]
            r["compl"], r["contam_checkm"] = float(m["Completeness"]), float(m["Contamination"])
            r["hq"] = r["type"] == "MAG" and r["compl"] >= 90 and r["contam_checkm"] <= 5
            R[r["genome"]] = r
    F = []
    for f in sorted(glob.glob(os.path.join(W, "screen", "part*.flagged.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            for k in FINT:
                r[k] = int(r[k])
            r["gc"], r["ani"], r["af"] = float(r["gc"]), float(r["ani"]), float(r["af"])
            F.append(r)
    pilot = set(open(os.path.join(U, "pilot", "genomes.txt")).read().split())
    return M, R, F, pilot


def gstats(rows):
    n = len(rows)
    if n == 0:
        return {"n": 0}
    tot = sum(r["total_bp"] for r in rows)

    def frac(f):
        return wilson(sum(1 for r in rows if f(r)), n)
    f = np.array([r["foreign"] for r in rows])
    c = np.array([r["contam"] for r in rows])
    return {"n": n, "bp": tot,
            "shared_any": frac(lambda r: r["foreign_n"] > 0), "mobile_any": frac(lambda r: r["mobile_n"] > 0),
            "core_contig_any": frac(lambda r: r["core_n"] > 0),
            "foreign_ge1": frac(lambda r: r["foreign"] >= 0.01), "foreign_ge5": frac(lambda r: r["foreign"] >= 0.05),
            "contam_any": frac(lambda r: r["anchored_bp"] > 0), "contam_ge1": frac(lambda r: r["contam"] >= 0.01),
            "contam_ge5": frac(lambda r: r["contam"] >= 0.05), "strict_any": frac(lambda r: r["anchored_strict_bp"] > 0),
            "strict_ge1": frac(lambda r: r["contam_strict"] >= 0.01),
            "foreign_mean": float(f.mean()), "contam_mean": float(c.mean()), "contam_p99": float(np.quantile(c, 0.99)),
            "foreign_share": sum(r["foreign_bp"] for r in rows) / max(1, tot),
            "mobile_share": sum(r["mobile_bp"] for r in rows) / max(1, tot),
            "chrom_share": sum(r["chrom_bp"] for r in rows) / max(1, tot),
            "anchored_share": sum(r["anchored_bp"] for r in rows) / max(1, tot),
            "hit_share": sum(r["hit_bp"] for r in rows) / max(1, tot),
            "species_hit_share": sum(r["species_hit_bp"] for r in rows) / max(1, tot),
            "host_src": dict(collections.Counter(r["host_src"] for r in rows))}


def species_matched(R, M):
    """isolates vs MAGs of the same species (host called from a self or species-level reference)"""
    sp = collections.defaultdict(lambda: {"Isolate": [], "MAG": []})
    for g, r in R.items():
        if r["host_src"] in ("self", "species"):
            sp[M[g]["Species_rep"]][r["type"]].append(r)
    el = {k: v for k, v in sp.items() if len(v["Isolate"]) >= 3 and len(v["MAG"]) >= 3}
    out = {"species": len(el), "isolates": sum(len(v["Isolate"]) for v in el.values()),
           "mags": sum(len(v["MAG"]) for v in el.values())}
    tests = {"shared_any": lambda r: r["foreign_n"] > 0, "foreign_ge1": lambda r: r["foreign"] >= 0.01,
             "foreign_ge5": lambda r: r["foreign"] >= 0.05, "contam_any": lambda r: r["anchored_bp"] > 0,
             "contam_ge1": lambda r: r["contam"] >= 0.01, "contam_ge5": lambda r: r["contam"] >= 0.05}
    for name, f in tests.items():
        num = den = pi = pm = 0.0
        up = down = tie = 0
        for v in el.values():
            a = np.array([f(r) for r in v["Isolate"]], float)
            b = np.array([f(r) for r in v["MAG"]], float)
            w = len(a) * len(b) / (len(a) + len(b))
            num += w * (b.mean() - a.mean())
            den += w
            pi += w * a.mean()
            pm += w * b.mean()
            up += b.mean() > a.mean()
            down += b.mean() < a.mean()
            tie += b.mean() == a.mean()
        p = sst.binomtest(up, up + down).pvalue if up + down else 1.0
        den = den or 1.0
        out[name] = {"mh_rd": num / den, "iso_rate": pi / den, "mag_rate": pm / den, "mag_higher": up, "iso_higher": down,
                     "tie": tie, "sign_p": p}
    return out


def same_sample(R, M):
    """Is the source genus of a MAG's foreign sequence also binned (another MAG) from the same sample?"""
    bins = [(2, 3), (4, 6), (7, 10), (11, 20), (21, 10 ** 9)]

    def cbin(k):
        return next(i for i, (lo, hi) in enumerate(bins) if lo <= k <= hi)
    samp = collections.defaultdict(list)
    for g, m in M.items():
        if m["Genome_type"] == "MAG":
            samp[m["Sample_accession"]].append(g)
    S, pair, single = {}, collections.Counter(), collections.Counter()
    for s, gs in samp.items():
        if len(gs) < 2:
            continue
        st_ = (M[gs[0]]["Study_accession"], cbin(len(gs)))
        gen = {genus_of(M[g]["Lineage"]) for g in gs} - {""}
        S[s] = st_
        for a in gen:
            single[(st_, a)] += 1
            for b in gen:
                if a != b:
                    pair[(st_, a, b)] += 1
    res = {"samples_with_2plus": len(S)}
    units = {"contaminated": lambda r: r["anchored_bp"] > 0 and r["top_anchor"].split(" ")[0],
             "shared_only": lambda r: r["anchored_bp"] == 0 and r["chrom_bp"] > 0 and r["top_source"]}
    for lab, fx in units.items():
        O = E = V = 0.0
        n = 0
        for g, r in R.items():
            if r["type"] != "MAG" or r["host_src"] not in ("self", "species"):
                continue
            x = fx(r)
            s = M[g]["Sample_accession"]
            if not x or s not in S:
                continue
            h = genus_of(M[g]["Lineage"])
            if not h or x == h:
                continue
            st_ = S[s]
            nh = single[(st_, h)] - 1
            if nh < 1:
                continue
            o = int(x in {genus_of(M[y]["Lineage"]) for y in samp[s] if y != g})
            e = (pair[(st_, h, x)] - o) / nh
            n += 1
            O += o
            E += e
            V += e * (1 - e)
        res[lab] = {"mags": n, "observed": O / max(1, n), "expected": E / max(1, n), "ratio": O / max(E, 1e-9),
                    "z": (O - E) / math.sqrt(max(V, 1e-9))}
    return res


FL = {}  # contig -> flagged class (set before forking the scan workers)
TYPE = {}  # genome -> Genome_type
SAMPLE_HQ = set()


def scan_part(p):
    """aggregates per (type, length bin, group) of one part: contigs, bp, CDS, core, hypothetical, mobile-hallmark CDS;
    and the (length, core) lists of the sampled high-quality MAGs (anchoring sensitivity)"""
    agg = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    donors = collections.defaultdict(list)
    for line in open(os.path.join(W, "fa", f"part{p:03d}.contigs.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[4] != "1":
            continue
        L, cds = int(f[2]), int(f[5])
        cls = FL.get(f[0])
        grp = "flag_chrom" if cls in ("core", "accessory") else ("flag_mobile" if cls == "mobile" else
                                                                ("bg" if cls is None and cds else None))
        if grp:
            a = agg[(TYPE[f[1]], lbin(L), grp)]
            a[0] += 1
            a[1] += L
            a[2] += cds
            a[3] += int(f[13])
            a[4] += int(f[12])
            a[5] += int(f[11])
        if f[1] in SAMPLE_HQ and cls is None and cds:
            donors[f[1]].append((L, int(f[13])))
    return p, {k: v for k, v in agg.items()}, dict(donors)


def core_mixture(parts, B=300, seed=3):
    """share of flagged non-mobile MAG bases that is contamination-like, from core genes per kb (bootstrap over parts)"""
    def est(sel):
        A = collections.defaultdict(lambda: np.zeros(6))
        for p in sel:
            for k, v in parts[p].items():
                A[k] += v
        out, num, den = {}, 0.0, 0.0
        for b in range(len(LBINS)):
            def d(t, g):
                v = A[(t, b, g)]
                return v[3] / max(1, v[1]) * 1000, v[1]
            (di, _), (dm, Lm), (db, _) = d("Isolate", "flag_chrom"), d("MAG", "flag_chrom"), d("MAG", "bg")
            sh = min(1.0, max(0.0, (dm - di) / max(1e-12, db - di)))
            out[b] = (sh, Lm, di, dm, db)
            num += sh * Lm
            den += Lm
        return out, num / max(1, den)
    keys = sorted(parts)
    point, overall = est(keys)
    rng = random.Random(seed)
    boots = [est([rng.choice(keys) for _ in keys]) for _ in range(B)]
    res = {"overall": overall, "overall_ci": [float(np.quantile([x[1] for x in boots], q)) for q in (0.025, 0.975)], "bins": []}
    for b, (lo, hi) in enumerate(LBINS):
        sh, Lm, di, dm, db = point[b]
        bs = [x[0][b][0] for x in boots]
        res["bins"].append({"lo": lo, "hi": hi, "share": sh, "ci": [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))],
                            "mag_flag_bp": Lm, "core_per_kb_iso_flag": di, "core_per_kb_mag_flag": dm, "core_per_kb_mag_bg": db})
    return res


def anchor_sensitivity(donors, sizes, levels=(0.005, 0.01, 0.02, 0.05, 0.1), n=4000, seed=9):
    """P(>= 1 / >= 2 core genes) among random unflagged contigs of a high-quality MAG adding up to a contamination level"""
    rng = random.Random(seed)
    ds = [v for v in donors.values() if len(v) >= 10]
    out = []
    for lv in levels:
        k1 = k2 = 0
        for _ in range(n):
            d = rng.choice(ds)[:]
            rng.shuffle(d)
            target = lv * rng.choice(sizes)
            tot = core = 0
            for L, c in d:
                if tot >= target:
                    break
                if tot + L > 1.25 * target and tot > 0:
                    continue
                tot += L
                core += c
            k1 += core >= 1
            k2 += core >= 2
        out.append({"level": lv, "p_core1": k1 / n, "p_core2": k2 / n})
    return out


def contig_table(genomes):
    """every contig of the given genomes from the prepare tables: {genome: [(contig, length, gc, screened, cds, core)]}"""
    out = collections.defaultdict(list)
    parts = sorted({zlib.crc32(g.encode()) % NPART for g in genomes})
    for p in parts:
        for line in open(os.path.join(W, "fa", f"part{p:03d}.contigs.tsv")):
            f = line.rstrip("\n").split("\t")
            if f[1] in genomes:
                out[f[1]].append((f[0], int(f[2]), float(f[3]), int(f[4]), int(f[5]), int(f[13]), int(f[10])))
    return out


def controls(F, CT=None, seed=11):
    """one unflagged contig per flagged coding contig: same genome, nearest length (screened)"""
    flagged = {r["contig"] for r in F}
    want = [r for r in F if r["cls"] != "noncoding"]
    CT = CT or contig_table({r["genome"] for r in want})
    rng = random.Random(seed)
    used = set()
    pairs = []
    for r in want:
        cand = [c for c in CT[r["genome"]] if c[3] and c[0] not in flagged and c[0] not in used]
        if not cand:
            continue
        best = min(abs(math.log(c[1] / r["length"])) for c in cand)
        near = [c for c in cand if abs(math.log(c[1] / r["length"])) <= best + 1e-9]
        c = rng.choice(near)
        used.add(c[0])
        pairs.append((r, c))
    return pairs


def gc_test(pairs, M, GCR):
    out = {}
    for lab, sel in (("anchored", lambda r: r["anchored"]), ("shared_chrom", lambda r: r["cls"] in ("core", "accessory") and not r["anchored"]),
                     ("mobile", lambda r: r["cls"] == "mobile")):
        for t in ("Isolate", "MAG"):
            rows = [(r, c) for r, c in pairs if sel(r) and M[r["genome"]]["Genome_type"] == t and r["length"] >= GC_MIN_LEN
                    and c[1] >= GC_MIN_LEN / 2]
            if not rows:
                continue
            dh_f = np.array([abs(r["gc"] - float(M[r["genome"]]["GC_content"]) / 100) for r, _ in rows])
            dh_c = np.array([abs(c[2] - float(M[r["genome"]]["GC_content"]) / 100) for r, c in rows])
            inf = []
            for r, _ in rows:
                hg = float(M[r["genome"]]["GC_content"]) / 100
                sg = GCR.get(r["ref"])
                if sg is not None and abs(sg - hg) >= 0.05:
                    inf.append(abs(r["gc"] - sg) < abs(r["gc"] - hg))
            out[f"{lab}_{t}"] = {"pairs": len(rows), "flag_dev_median": float(np.median(dh_f)), "ctrl_dev_median": float(np.median(dh_c)),
                                 "flag_dev_ge5": float((dh_f >= 0.05).mean()), "ctrl_dev_ge5": float((dh_c >= 0.05).mean()),
                                 "informative": len(inf), "closer_to_source": wilson(sum(inf), len(inf))}
    return out


def recurrence(F, M):
    """pilot test (data/uhgg3/consp_pilot): flagged contigs of pilot genomes searched against <= 300 genomes of their species"""
    cls = {r["contig"]: r for r in F}
    res = collections.defaultdict(collections.Counter)
    for d in glob.glob(os.path.join(U, "consp_pilot", "MGYG*", "hits.tsv")):
        car, seen = collections.defaultdict(set), set()
        for line in open(d):
            f = line.rstrip("\n").split("\t")
            if f[0] == "contig":
                continue
            seen.add(f[0])
            if f[1] != "*" and f[9] == "ok" and float(f[3]) >= 0.95 and float(f[4]) >= 0.5:
                car[f[0]].add(f[1])
        for c in seen:
            if c not in cls:
                continue
            r = cls[c]
            g, st = r["genome"], M[r["genome"]]["Study_accession"]
            ot = [x for x in car[c] if x != g]
            k = "recurrent" if any(M[x]["Study_accession"] != st for x in ot) else ("same_study" if ot else "sporadic")
            lab = "anchored" if r["anchored"] else r["cls"]
            res[f"{lab}_{r['type'] if 'type' in r else M[g]['Genome_type']}"][k] += 1
    return {k: {"n": sum(v.values()), **{x: v[x] / sum(v.values()) for x in ("recurrent", "same_study", "sporadic")}} for k, v in res.items()}


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(f"{x:.5g}" if isinstance(x, float) else str(x) for x in r) + "\n")


def stats():
    os.makedirs(OUT, exist_ok=True)
    M, R, F, pilot = load()
    T, GEN, GCR = gtdb_tax()
    S = {}
    rows = list(R.values())
    mags = [r for r in rows if r["type"] == "MAG"]
    iso = [r for r in rows if r["type"] == "Isolate"]
    known = ("self", "species")
    grp = {"all": rows, "isolates": iso, "mags": mags, "hq_mags": [r for r in mags if r["hq"]],
           "mq_mags": [r for r in mags if not r["hq"]],
           "reps": [r for r in rows if r["rep"]], "reps_mag": [r for r in mags if r["rep"]], "reps_iso": [r for r in iso if r["rep"]],
           "isolates_known": [r for r in iso if r["host_src"] in known], "mags_known": [r for r in mags if r["host_src"] in known],
           "mags_weak": [r for r in mags if r["host_src"] not in known],
           "pilot_isolates": [r for r in iso if r["genome"] in pilot], "pilot_mags": [r for r in mags if r["genome"] in pilot],
           "heldout_isolates": [r for r in iso if r["genome"] not in pilot], "heldout_mags": [r for r in mags if r["genome"] not in pilot]}
    S["groups"] = {k: gstats(v) for k, v in grp.items()}
    S["species_matched"] = species_matched(R, M)
    # whole-catalogue contig scan: core-gene mixture and anchoring sensitivity
    FL.update({r["contig"]: r["cls"] for r in F})
    TYPE.update({g: m["Genome_type"] for g, m in M.items()})
    hq = sorted(r["genome"] for r in mags if r["hq"] and r["foreign_n"] == 0)
    SAMPLE_HQ.update(random.Random(4).sample(hq, min(3000, len(hq))))
    parts, donors = {}, {}
    with multiprocessing.get_context("fork").Pool(min(NPART, int(THREADS))) as pool:
        for p, agg, dn in pool.imap_unordered(scan_part, range(NPART)):
            parts[p] = {k: np.array(v, float) for k, v in agg.items()}
            donors.update(dn)
    S["core_mixture"] = core_mixture(parts)
    S["anchor_sensitivity"] = anchor_sensitivity(donors, [r["total_bp"] for r in mags])
    prof = {}
    for t in ("Isolate", "MAG"):
        for g_ in ("flag_chrom", "flag_mobile", "bg"):
            for b in range(len(LBINS)):
                v = sum((parts[p].get((t, b, g_), np.zeros(6)) for p in parts), np.zeros(6))
                prof[f"{t}_{g_}_{b}"] = {"n": int(v[0]), "bp": int(v[1]), "core_per_kb": v[3] / max(1, v[1]) * 1000,
                                         "hypo_frac": v[4] / max(1, v[2]), "mobkw_per_kb": v[5] / max(1, v[1]) * 1000}
    S["contig_profiles"] = prof
    # MGnify quality columns (CheckM2 contamination, completeness) vs aif
    cb = []
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5.01)):
        sel = [r for r in mags if lo <= r["contam_checkm"] < hi]
        cb.append({"lo": lo, "hi": min(hi, 5), "n": len(sel), "contam_any": wilson(sum(r["anchored_bp"] > 0 for r in sel), len(sel)),
                   "contam_ge1": wilson(sum(r["contam"] >= 0.01 for r in sel), len(sel)),
                   "foreign_ge5": wilson(sum(r["foreign"] >= 0.05 for r in sel), len(sel))})
    S["by_checkm_contam"] = cb
    kb = []
    for lo, hi in ((50, 70), (70, 90), (90, 100.01)):
        sel = [r for r in mags if lo <= r["compl"] < hi]
        kb.append({"lo": lo, "hi": min(hi, 100), "n": len(sel), "contam_ge1": wilson(sum(r["contam"] >= 0.01 for r in sel), len(sel)),
                   "foreign_ge5": wilson(sum(r["foreign"] >= 0.05 for r in sel), len(sel))})
    S["by_completeness"] = kb
    nb = []
    for lo, hi in ((0, 50), (50, 150), (150, 300), (300, 10 ** 9)):
        sel = [r for r in mags if lo <= int(M[r["genome"]]["N_contigs"]) < hi]
        nb.append({"lo": lo, "hi": hi, "n": len(sel), "contam_ge1": wilson(sum(r["contam"] >= 0.01 for r in sel), len(sel)),
                   "foreign_ge5": wilson(sum(r["foreign"] >= 0.05 for r in sel), len(sel))})
    S["by_contigs"] = nb
    mk = [r for r in mags if r["host_src"] in known]
    S["spearman_checkm"] = {"contam": float(sst.spearmanr([r["contam_checkm"] for r in mk], [r["contam"] for r in mk])[0]),
                            "foreign": float(sst.spearmanr([r["contam_checkm"] for r in mk], [r["foreign"] for r in mk])[0])}
    hqc = [r for r in mags if r["hq"] and r["contam"] >= 0.01]
    S["hq_contaminated"] = {"ge1": len(hqc), "ge5": sum(r["contam"] >= 0.05 for r in hqc), "of_hq": sum(r["hq"] for r in mags),
                            "checkm_median": float(np.median([r["contam_checkm"] for r in hqc] or [0])),
                            "aif_median": float(np.median([r["contam"] for r in hqc] or [0]))}
    ag = [genus_of(M[r["genome"]]["Lineage"]) == r["host_genus"] for r in rows if r["host_src"] in known]
    agw = [genus_of(M[r["genome"]]["Lineage"]) == r["host_genus"] for r in rows if r["host_src"] == "weak"]
    S["host_vs_mgnify"] = {"known": wilson(sum(ag), len(ag)), "weak": wilson(sum(agw), len(agw))}
    # GTDB label inconsistencies: references conspecific with a genome (ANI >= 95% over >= 50% of it) under another genus
    xr = collections.defaultdict(lambda: [0, collections.Counter()])
    for r in rows:
        if r["equiv_xgenus"]:
            for a in r["equiv_xgenus_refs"].split(";"):
                xr[a][0] += 1
                xr[a][1][r["host_genus"]] += 1
    NB = ncbi_meta(set(xr))
    xl = sorted(xr.items(), key=lambda kv: -kv[1][0])
    write_tsv(os.path.join(OUT, "label_inconsistent_refs.tsv"),
              ["reference", "gtdb_species", "gtdb_family", "ncbi_species", "checkm2_completeness", "checkm2_contamination", "genomes",
               "host_genera"],
              [[a, T.get(a, ("?",) * 4)[3], T.get(a, ("?",) * 4)[1], *NB.get(a, ("?", "?", "?")), v[0],
                ";".join(f"{g}:{n}" for g, n in v[1].most_common(5))] for a, v in xl])
    fam = collections.Counter()
    for a, v in xr.items():
        fam[T.get(a, ("?",) * 4)[1]] += v[0]
    S["label_inconsistent"] = {"genomes": {t: sum(1 for r in rows if r["type"] == t and r["equiv_xgenus"]) for t in ("Isolate", "MAG")},
                               "genomes_with_equiv": {t: sum(1 for r in rows if r["type"] == t and r["equiv_n"]) for t in ("Isolate", "MAG")},
                               "refs": len(xr), "refs_ge10": sum(1 for v in xr.values() if v[0] >= 10),
                               "families": fam.most_common(10)}
    # contamination sources (anchored) and shared-sequence genus pairs (everything else)
    anc = [r for r in F if r["anchored"] and R[r["genome"]]["type"] == "MAG"]
    sbp, sg = collections.Counter(), collections.defaultdict(set)
    for r in anc:
        sbp[r["source_genus"]] += r["length"]
        sg[r["source_genus"]].add(r["genome"])
    tb = sum(sbp.values())
    write_tsv(os.path.join(OUT, "contamination_sources.tsv"), ["source_genus", "family", "phylum", "bp", "share", "mags"],
              [(k, GEN.get(k, ("?", "?"))[1], GEN.get(k, ("?", "?"))[0], v, v / max(1, tb), len(sg[k])) for k, v in sbp.most_common(40)])
    sh = [r for r in F if not r["anchored"]]
    pairs = collections.defaultdict(set)
    xp = xf = tbp = 0
    for r in sh:
        pairs[(r["host_genus"], r["source_genus"])].add(r["genome"])
        tbp += r["length"]
        xp += r["length"] * (GEN.get(r["host_genus"], ("?",))[0] != r["source_phylum"])
        xf += r["length"] * (r["host_family"] != r["source_family"])
    S["shared"] = {"contigs": len(sh), "bp": tbp, "genomes": len({r["genome"] for r in sh}), "genus_pairs": len(pairs),
                   "cross_family_share": xf / max(1, tbp), "cross_phylum_share": xp / max(1, tbp),
                   "by_class": dict(collections.Counter(r["cls"] for r in sh)),
                   "genera_connected": len({a for a, _ in pairs} | {b for _, b in pairs})}
    S["contamination_sources"] = {"contigs": len(anc), "bp": tb, "genera": len(sbp),
                                  "cross_family_share": sum(r["length"] for r in anc if r["host_family"] != r["source_family"]) / max(1, tb)}
    write_tsv(os.path.join(OUT, "shared_pairs.tsv"), ["host_genus", "source_genus", "genomes"],
              [(a, b, len(v)) for (a, b), v in sorted(pairs.items(), key=lambda x: -len(x[1]))[:60]])
    # species representatives
    reps = sorted((r for r in rows if r["rep"] and r["contam"] >= 0.01), key=lambda r: -r["contam"])
    write_tsv(os.path.join(OUT, "species_reps_contaminated.tsv"),
              ["genome", "type", "lineage", "completeness", "checkm2_contamination", "aif_contam", "aif_contam_strict", "anchored_species",
               "anchored_core", "top_anchor", "host_src"],
              [(r["genome"], r["type"], M[r["genome"]]["Lineage"], r["compl"], r["contam_checkm"], r["contam"], r["contam_strict"],
                r["anchored_species"], r["anchored_core"], r["top_anchor"], r["host_src"]) for r in reps])
    # AMR (AMRFinderPlus, species representatives only)
    repg = {r["genome"] for r in rows if r["rep"]}
    CT = contig_table(repg)
    fr = [r for r in F if r["genome"] in repg]
    S["amr_reps"] = {"reps": len(repg), "amr_genes": sum(c[6] for g in repg for c in CT[g]),
                     **{f"on_{k}": sum(r["amr_cds"] for r in fr if (r["cls"] == k and not r["anchored"])) for k in ("mobile", "core", "accessory", "noncoding")},
                     "on_anchored": sum(r["amr_cds"] for r in fr if r["anchored"]),
                     "reps_with_shared_amr": len({r["genome"] for r in fr if r["amr_cds"] and not r["anchored"]}),
                     "reps_with_anchored_amr": len({r["genome"] for r in fr if r["amr_cds"] and r["anchored"]})}
    # studies (>= 1000 MAGs)
    bys = collections.defaultdict(list)
    for r in mags:
        bys[M[r["genome"]]["Study_accession"]].append(r["contam"] >= 0.01)
    big = {k: v for k, v in bys.items() if len(v) >= 1000}
    if len(big) > 1:
        chi = sst.chi2_contingency([[sum(v), len(v) - sum(v)] for v in big.values()])
        S["studies"] = {"n": len(big), "rates": sorted([[k, len(v), sum(v) / len(v)] for k, v in big.items()], key=lambda x: -x[2]),
                        "chi2_p": float(chi[1])}
    S["euk_known_issue"] = [[g, R[g]["host_src"], R[g]["hit_bp"] / max(1, R[g]["total_bp"])] for g in EUK if g in R]
    hf = np.array([r["hit_bp"] / max(1, r["total_bp"]) for r in rows])
    S["low_hit"] = {"lt10": int((hf < 0.10).sum()), "lt50": int((hf < 0.5).sum())}
    S["same_sample"] = same_sample(R, M)
    S["gc"] = gc_test(controls(F), M, GCR)
    S["recurrence_pilot"] = recurrence(F, M)
    edges = [0, 1e-4, 0.001, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 10]
    S["hist"] = {"edges": edges, **{f"{k}_{lab}": np.histogram([r[k] for r in sel], edges)[0].tolist()
                                    for k in ("foreign", "contam") for lab, sel in (("isolates", iso), ("mags", mags))}}
    tp = os.path.join(W, "search_time.json")
    if os.path.exists(tp):
        t = json.load(open(tp))
        S["search"] = {"parts": len(t), "s": sum(v["s"] for v in t.values()), "max_kb": max(v["kb"] for v in t.values()),
                       "bytes": sum(v["bytes"] for v in t.values()), "threads": int(THREADS)}
    pp = os.path.join(W, "prepare.json")
    if os.path.exists(pp):
        S["prepare"] = json.load(open(pp))
    S["flagged_contigs"] = dict(collections.Counter(r["cls"] for r in F))
    S["flagged_anchored"] = sum(r["anchored"] for r in F)
    json.dump(S, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=_np)
    cols = ["genome", "contigs", "total_bp", "self_ref", "host_genus", "host_src", "foreign_bp", "mobile_bp", "core_bp", "accessory_bp",
            "chrom_bp", "anchored_bp", "anchored_strict_bp", "anchored_species", "top_source", "top_anchor", "foreign", "contam",
            "equiv_n", "equiv_xgenus_refs"]
    with gzip.open(os.path.join(OUT, "genomes.tsv.gz"), "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in sorted(R):
            r = R[g]
            fh.write("\t".join(f"{r[c]:.5f}" if c in ("foreign", "contam") else str(r[c]) for c in cols) + "\n")
    with gzip.open(os.path.join(OUT, "flagged_contigs.tsv.gz"), "wt") as fh:
        k = ["contig", "genome", "length", "gc", "cls", "anchored", "strict", "ref", "ani", "af", "source_species", "source_genus",
             "host_genus", "cds", "mobkw_cds", "hypo_cds", "core_cds", "amr_cds"]
        fh.write("\t".join(k) + "\n")
        for r in F:
            fh.write("\t".join(str(r[c]) for c in k) + "\n")
    g = S["groups"]
    for k in ("isolates", "mags", "hq_mags", "reps", "heldout_isolates", "heldout_mags"):
        d = g[k]
        if d["n"]:
            print(f"{k:17s} n={d['n']:7d} shared={d['shared_any'][0]:.4f} foreign>=1%={d['foreign_ge1'][0]:.4f} "
                  f"contam any={d['contam_any'][0]:.4f} >=1%={d['contam_ge1'][0]:.4f} >=5%={d['contam_ge5'][0]:.4f} "
                  f"mean={d['contam_mean']:.5f}")
    print("core mixture", json.dumps(S["core_mixture"], default=_np))
    print("anchor sensitivity", S["anchor_sensitivity"])
    print("species matched", {k: v for k, v in S["species_matched"].items()})
    print("same sample", S["same_sample"])
    print("gc", S["gc"])
    print("recurrence", S["recurrence_pilot"])
    print("amr", S["amr_reps"])


WANT = {}


def extract_part(job):
    """the wanted contigs of one FASTA part (WANT is inherited from the parent by fork)"""
    p, kd = job
    out = os.path.join(kd, f"part{p:03d}.fa")
    keep = False
    with open(os.path.join(W, "fa", f"part{p:03d}.fa")) as fi, open(out, "w") as fo:
        for line in fi:
            if line[0] == ">":
                keep = line[1:].split()[0] in WANT
            if keep:
                fo.write(line)
    return out


def kraken():
    """Kraken2 genus of flagged contigs and of their length-matched controls"""
    os.makedirs(OUT, exist_ok=True)
    M, R, F, _ = load()
    pairs = controls(F)
    kd = os.path.join(W, "kraken2")
    os.makedirs(kd, exist_ok=True)
    want = {r["contig"]: ("flag", r) for r in F}
    for r, c in pairs:
        want[c[0]] = ("ctrl", r)
    fa = os.path.join(kd, "contigs.fa")
    if not os.path.exists(fa):
        WANT.update(want)
        parts = sorted({zlib.crc32(r["genome"].encode()) % NPART for r in F})
        with multiprocessing.get_context("fork").Pool(min(len(parts), int(THREADS))) as pool:
            outs = list(pool.imap_unordered(extract_part, [(p, kd) for p in parts]))
        with open(fa + ".tmp", "w") as fo:
            for o in sorted(outs):
                with open(o) as fi:
                    for line in fi:
                        fo.write(line)
                os.remove(o)
        os.replace(fa + ".tmp", fa)
    out = os.path.join(kd, "out.tsv")
    if not os.path.exists(out):
        subprocess.run(["cat", os.path.join(K2, "hash.k2d")], stdout=subprocess.DEVNULL, check=True)
        tf = os.path.join(kd, "time.txt")
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", tf, "kraken2", "--db", K2, "--threads", THREADS,
                        "--output", out + ".tmp", "--report", os.path.join(kd, "report.txt"), fa], check=True,
                       stderr=open(os.path.join(kd, "kraken2.log"), "w"))
        os.replace(out + ".tmp", out)
    nodes, names = {}, {}
    for line in open(os.path.join(K2, "taxonomy", "nodes.dmp")):
        p = [x.strip() for x in line.split("|")]
        nodes[int(p[0])] = (int(p[1]), p[2])
    for line in open(os.path.join(K2, "taxonomy", "names.dmp")):
        p = [x.strip() for x in line.split("|")]
        names[int(p[0])] = p[1]

    def lineage(t):
        L = {}
        while t != 1 and t in nodes:
            L[nodes[t][1]] = names[t][3:]
            t = nodes[t][0]
        return L
    A = {}
    for line in open(out):
        p = line.split("\t", 3)
        A[p[1]] = lineage(int(p[2])) if p[0] == "C" else {}
    groups = collections.defaultdict(list)
    for cid, (kind, r) in want.items():
        L = A.get(cid, {})
        gk, fk = L.get("genus"), L.get("family")
        if kind == "flag":
            lab = "anchored" if r["anchored"] else r["cls"]
            v = ("source" if gk == r["source_genus"] else "host" if gk == r["host_genus"] else "other" if gk
                 else "source_family" if fk and fk == r["source_family"] and fk != r["host_family"]
                 else "host_family" if fk and fk == r["host_family"] else "above" if L else "unclassified")
            key = f"flag_{lab}_{M[r['genome']]['Genome_type']}"
        else:
            v = ("host" if gk == r["host_genus"] else "source" if gk == r["source_genus"] else "other" if gk
                 else "host_family" if fk and fk == r["host_family"] else "above" if L else "unclassified")
            key = f"control_{M[r['genome']]['Genome_type']}"
        groups[key].append(v)
        if r["host_family"] != "Enterobacteriaceae":  # GTDB R232 genus labels conflict there (methods)
            groups[key + "_noEnt"].append(v)
    res = {}
    for k, v in groups.items():
        c = collections.Counter(v)
        res[k] = {"n": len(v), **{x: c[x] / len(v) for x in c}}
    tf = os.path.join(kd, "time.txt")
    if os.path.exists(tf):
        s, kb = open(tf).read().split()[-2:]
        res["time"] = {"s": float(s), "kb": int(kb), "contigs": len(want)}
    json.dump(res, open(os.path.join(OUT, "kraken2_check.json"), "w"), indent=1, default=_np)
    for k in sorted(res):
        print(k, {x: (round(y, 4) if isinstance(y, float) else y) for x, y in res[k].items()})


if __name__ == "__main__":
    {"stats": stats, "kraken": kraken}[sys.argv[1]]()
