#!/usr/bin/env python3
"""Global atlas of recent horizontal gene transfer across GTDB R232 species representatives with aif.

Every representative genome is cut into 5 kb windows, and each window is searched against the GTDB database
(aif db search, all 199,923 representatives). A window that aligns to ANOTHER species' representative at >= 99%
identity over >= 500 bp carries recently shared sequence (Smillie et al. 2011 used >= 99% over >= 500 bp with BLAST
on 2,235 genomes). Windows keep their contig coordinates, so shared blocks can be told apart from whole foreign
contigs (contamination) by their native flanks.

Stages (data in data/hgt/, results in results_hgt/):
  select              genome table: taxonomy, category (isolate / MAG / SAG), assembly level, habitat class
  windows SET         5 kb windows of the set's genomes into FASTA parts (SET: complete | isolates | all | a list)
  search SET          aif db search of every part, filtered on the fly (filter stage) -> hits/partNNN.tsv
  filter              stdin -> stdout: non-self hits with status ok, ANI >= 0.97 and >= 500 aligned bp

Usage: /usr/bin/python3 -I bench_real/hgt_atlas.py <stage> [SET]
"""
import csv
import gzip
import json
import multiprocessing
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "data", "gtdb")
REPS = os.path.join(G, "reps", "gtdb_genomes_reps_r232", "database")
DB = os.path.join(G, "aif", "gtdb_0.88.aifd")
HD = os.path.join(ROOT, "data", "hgt")
RES = os.path.join(ROOT, "results_hgt")
AIF = os.path.join(ROOT, "bin", "aif")
W = 5000          # window length
WMIN = 1000       # shortest window kept (contig ends)
THREADS = os.environ.get("THREADS", "48")
NCAND = os.environ.get("HGT_NCAND", "32")   # references evaluated by the model per window (and hits reported)
ANI_KEEP, BP_KEEP = 0.97, 500
GZ = os.environ.get("HGT_GZ", "0") == "1"  # write window parts gzip-compressed (level 1); aif reads them directly

# habitat classes from free-text isolation sources (first matching rule wins)
HABITAT_RULES = [
    ("gut", r"rumen|ruminant|caec|cecum|\bgut\b|\bfec|\bfaec|stool|intestin|\bdung|manure|gastrointestinal|digestive|"
            r"cloaca|dropping|\bcolon\b|meconium|rect(al|um)\b|methane producing sheep"),
    ("host_other", r"\b(oral|saliva|dental|plaque|tongue|tooth|teeth|gingiv|subgingiv|skin|sputum|nasal|nasopharyn|throat|"
                   r"lung|bronch|blood|urin|vagin|cervi|wound|abscess|clinical|patient|hospital|cerebrospinal|biopsy|"
                   r"pus\b|ear\b|eye|bile|peritone|catheter|human|homo sapiens|tissue|lesion|infection)"),
    ("food", r"food|milk|cheese|ferment|kimchi|meat|dairy|yogurt|yoghurt|wine|beer|sausage|sourdough|vinegar|pickle"),
    ("engineered", r"wastewater|waste water|sludge|bioreactor|digest|compost|biogas|sewage|anammox|reactor|landfill|"
                   r"treatment plant|biofilter"),
    ("marine", r"marine|\bsea\b|ocean|seawater|sea water|coral|hydrothermal|black smoker|\bvent\b|cold seep|saline water|"
               r"sponge|estuar|salt marsh|mangrove|seagrass|tidal|beach|macroalg|seaweed|pyropia|kelp|\bbay\b"),
    ("plant", r"rhizosph|\broot|leaf|leaves|phyllosph|plant|nodule|seed|\bstem|flower|fruit|endophyt|rhizoplane|tuber|\bwood"),
    ("extreme", r"hot spring|thermal|hypersaline|acid mine|soda lake|salt|brine|geyser|volcan|glacier|\bice\b|permafrost|"
                r"desert|stromatolite|microbial mat"),
    ("subsurface", r"groundwater|aquifer|subsurface|fracture fluid|borehole|tunnel|aspo|olkiluoto|\bmine\b|\brock|"
                   r"deep biosphere|cave|\bwell\b|oil reservoir|oil field"),
    ("freshwater", r"freshwater|fresh water|lake|river|pond|spring|stream|reservoir|drinking water|creek|wetland"),
    ("soil", r"soil|paddy|agricultur|grassland|peat|forest"),
    ("sediment", r"sediment|\bmud\b"),
]
HAB_RE = [(h, re.compile(p)) for h, p in HABITAT_RULES]
NONE = {"", "none", "missing", "not collected", "not applicable", "unknown", "n/a", "na", "metagenome", "not provided",
        "not available", "other"}


def habitat(src):
    s = src.strip().lower()
    if s in NONE or s.startswith("not provided") or s.startswith("missing"):
        return "unknown"
    for h, rx in HAB_RE:
        if rx.search(s):
            return h
    return "other"


def select():
    os.makedirs(HD, exist_ok=True)
    paths = {}
    for d, _, fs in os.walk(REPS):
        for f in fs:
            if f.endswith("_genomic.fna.gz"):
                paths[f[:-len("_genomic.fna.gz")]] = os.path.join(d, f)
    out = []
    for f in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        for r in csv.DictReader(gzip.open(os.path.join(G, f), "rt"), delimiter="\t"):
            if r["gtdb_representative"] != "t":
                continue
            acc = r["accession"][3:]
            cat = {"none": "isolate", "derived from single cell": "SAG"}.get(r["ncbi_genome_category"], "MAG")
            out.append([acc, r["gtdb_taxonomy"], cat, r["ncbi_assembly_level"], r["genome_size"], r["contig_count"],
                        r["checkm2_completeness"], r["checkm2_contamination"], habitat(r["ncbi_isolation_source"]),
                        r["ncbi_isolation_source"].replace("\t", " ").strip(), paths.get(acc, "")])
    with open(os.path.join(HD, "genomes.tsv"), "w") as f:
        f.write("accession\ttaxonomy\tcategory\tassembly_level\tgenome_size\tcontigs\tcheckm2_completeness\t"
                "checkm2_contamination\thabitat\tisolation_source\tpath\n")
        for r in out:
            f.write("\t".join(r) + "\n")
    missing = sum(1 for r in out if not r[-1])
    print(f"{len(out)} representatives, {missing} without a genome file")


def load_genomes():
    return list(csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t"))


def set_members(name):
    G_ = load_genomes()
    if name == "complete":
        return [g for g in G_ if g["category"] == "isolate" and g["assembly_level"] in ("Complete Genome", "Chromosome")]
    if name == "isolates":
        return [g for g in G_ if g["category"] == "isolate"]
    if name == "all":
        return G_
    if os.path.exists(name):
        keep = {l.strip() for l in open(name) if l.strip()}
        return [g for g in G_ if g["accession"] in keep]
    raise SystemExit(f"unknown set {name}")


def window_genome(g):
    """windows of one genome: list of (name, seq); contig table rows"""
    acc = g["accession"]
    recs, contigs = [], []
    seqs, buf = [], []
    with gzip.open(g["path"], "rt") as fi:
        for line in fi:
            if line[0] == ">":
                if buf:
                    seqs.append("".join(buf))
                    buf = []
                contigs.append(line[1:].split()[0])
            else:
                buf.append(line.strip())
        if buf:
            seqs.append("".join(buf))
    rows = []
    for ci, s in enumerate(seqs):
        nwin = 0
        for st in range(0, len(s), W):
            t = s[st:st + W]
            if len(t) < WMIN:
                continue
            recs.append(f">{acc}|{ci}|{st}|{len(t)}\n{t}\n")
            nwin += 1
        rows.append(f"{acc}\t{ci}\t{contigs[ci] if ci < len(contigs) else ci}\t{len(s)}\t{nwin}\n")
    return "".join(recs), rows


def write_part(job):
    p, gs, outd = job
    path = os.path.join(outd, f"part{p:03d}.fa" + (".gz" if GZ else ""))
    rows = []
    with (gzip.open(path + ".tmp", "wt", compresslevel=1) if GZ else open(path + ".tmp", "w")) as fo:
        for g in gs:
            txt, r = window_genome(g)
            fo.write(txt)
            rows += r
    os.replace(path + ".tmp", path)
    return p, rows


def windows(name, nparts):
    gs = [g for g in set_members(name) if g["path"]]
    d = os.path.join(HD, os.path.basename(name), "win")
    os.makedirs(d, exist_ok=True)
    gs.sort(key=lambda g: g["accession"])
    # balance parts by genome size
    parts = [[] for _ in range(nparts)]
    load = [0] * nparts
    for g in sorted(gs, key=lambda g: -int(g["genome_size"] or 0)):
        i = load.index(min(load))
        parts[i].append(g)
        load[i] += int(g["genome_size"] or 0)
    ext = ".fa.gz" if GZ else ".fa"
    jobs = [(p, parts[p], d) for p in range(nparts) if not os.path.exists(os.path.join(d, f"part{p:03d}{ext}"))]
    rows = []
    with multiprocessing.get_context("fork").Pool(min(32, max(1, len(jobs)))) as pool:
        for p, r in pool.imap_unordered(write_part, jobs):
            rows += r
            print(f"part {p} written", flush=True)
    ct = os.path.join(HD, os.path.basename(name), "contigs.tsv")
    with open(ct, "a") as f:
        for r in rows:
            f.write(r)
    print(f"{len(gs)} genomes, {sum(load) / 1e9:.2f} Gbp in {nparts} parts")


def filt():
    """stdin: aif db search output; stdout: non-self hits worth keeping"""
    out = sys.stdout
    for line in sys.stdin:
        p = line.split("\t", 12)
        if p[0] == "contig" or len(p) < 12 or p[9] != "ok":
            continue
        q = p[0].split("|")
        if q[0] == p[1]:
            continue
        ani, af = float(p[3]), float(p[4])
        if ani < ANI_KEEP or af * float(q[3]) < BP_KEEP:
            continue
        out.write("\t".join(p[:6]) + "\t" + p[10] + "\n")


def search(name):
    base = os.path.join(HD, os.path.basename(name))
    wd, hd = os.path.join(base, "win"), os.path.join(base, "hits")
    os.makedirs(hd, exist_ok=True)
    tlog = os.path.join(base, "search_time.json")
    times = json.load(open(tlog)) if os.path.exists(tlog) else {}
    for f in sorted(os.listdir(wd)):
        if not (f.endswith(".fa") or f.endswith(".fa.gz")):
            continue
        out = os.path.join(hd, f.split(".")[0] + ".tsv")
        if os.path.exists(out):
            continue
        t0 = time.time()
        nbytes = os.path.getsize(os.path.join(wd, f))
        tf = out + ".time"
        cmd = (f"/usr/bin/time -f '%e %M' -o {tf} {AIF} db search {DB} {os.path.join(wd, f)} -t {THREADS} "
               f"-n {NCAND} --max-cand {NCAND} | /usr/bin/python3 -I {os.path.abspath(__file__)} filter > {out}.tmp")
        subprocess.run(cmd, shell=True, check=True)
        os.replace(out + ".tmp", out)
        if os.environ.get("HGT_DELETE") == "1":  # parts are regenerable; free the space once searched
            os.remove(os.path.join(wd, f))
        e, m = open(tf).read().split()
        times[f] = {"wall": time.time() - t0, "aif_s": float(e), "kb": int(m), "bytes": nbytes}
        json.dump(times, open(tlog, "w"), indent=1)
        print(f, times[f], flush=True)


def main():
    st = sys.argv[1]
    if st == "select":
        select()
    elif st == "windows":
        windows(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 64)
    elif st == "search":
        search(sys.argv[2])
    elif st == "filter":
        filt()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
