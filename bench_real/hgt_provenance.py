#!/usr/bin/env python3
"""Sample provenance of the metagenome-assembled and single-cell representatives in the HGT atlas.

Two MAGs assembled from the same metagenome can share a contig by cross-binning or chimeric assembly, and different
projects often reassemble the same public samples. A pair of MAGs is "co-sample" when their source samples overlap.
Every MAG has its own NCBI BioSample; its source is named in the record's derived-from attributes, either as
BioSamples of the metagenome or as SRA runs/experiments/samples. This script fetches those records, resolves SRA
accessions to BioSamples, and writes each MAG's set of source BioSamples (hgt_ecology.py, HGT_ECO_XSAMPLE).

Stages (data in data/hgt/provenance/; raw NCBI records are kept gzipped in raw/ and parsed as data only):
  fetch     BioSample records of all non-isolate representatives (E-utilities, <= 3 requests/s, no API key)
  parse     derived-from accessions per MAG -> derived.tsv
  resolve   SRA run/experiment/sample accessions -> BioSample (SRA runinfo) -> sra_map.tsv
  resolve_samples  SRA sample accessions without runs at NCBI (mostly ENA ERS) -> BioSample -> sample_map.tsv
  sources   per MAG: source BioSamples -> sources.tsv (an unresolved SRA accession is kept as its own key);
            coverage statistics -> results_hgt/provenance.json

Usage: /usr/bin/python3 -I bench_real/hgt_provenance.py <stage>
"""
import collections
import csv
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HD = os.path.join(ROOT, "data", "hgt")
PD = os.path.join(HD, "provenance")
RAW = os.path.join(PD, "raw")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "hgt_provenance"
BATCH = 200
HABS = {"gut", "host_other", "food", "engineered", "marine", "plant", "extreme", "subsurface", "freshwater", "soil",
        "sediment"}
BS_RE = re.compile(r"\bSAM(?:N|EA|D)\d{6,}\b")
SRA_RE = re.compile(r"\b[SED]R[RXS]\d{6,}\b")
_last = [0.0]


def eutil(name, data, tries=6):
    """POST to an E-utility, at most ~3 requests per second, retrying on network and server errors"""
    data = dict(data, tool=TOOL)
    body = urllib.parse.urlencode(data).encode()
    for k in range(tries):
        wait = 0.36 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(urllib.request.Request(EU + name, data=body), timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            code = getattr(e, "code", None)
            if code is not None and code < 500 and code != 429:
                raise
            time.sleep(2 * (k + 1))
    raise RuntimeError(f"{name} failed after {tries} tries")


def search_ids(db, accs):
    term = " OR ".join(f"{a}[accn]" for a in accs)
    x = eutil("esearch.fcgi", {"db": db, "term": term, "retmax": str(20 * len(accs))})
    return re.findall(rb"<Id>(\d+)</Id>", x)


def genomes():
    G = list(csv.DictReader(open(os.path.join(HD, "genomes.tsv")), delimiter="\t"))
    return [g for g in G if g["category"] != "isolate"]


def biosample_of():
    want = {g["accession"] for g in genomes()}
    out = {}
    for f in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"):
        with gzip.open(os.path.join(ROOT, "data", "gtdb", f), "rt") as fi:
            h = fi.readline().rstrip("\n").split("\t")
            ia, ib = h.index("accession"), h.index("ncbi_biosample")
            for l in fi:
                x = l.rstrip("\n").split("\t")
                a = x[ia][3:]
                if a in want:
                    out[a] = x[ib]
    return out


def fetch():
    os.makedirs(RAW, exist_ok=True)
    G = genomes()
    bs = biosample_of()
    # habitat-labelled genomes first: the ecology analysis needs those
    order = sorted(G, key=lambda g: (g["habitat"] not in HABS, g["accession"]))
    accs = []
    seen = set()
    for g in order:
        b = bs.get(g["accession"], "")
        if BS_RE.fullmatch(b) and b not in seen:
            seen.add(b)
            accs.append(b)
    json.dump({"genomes": len(G), "with_biosample": sum(1 for g in G if BS_RE.fullmatch(bs.get(g["accession"], ""))),
               "biosamples": len(accs)}, open(os.path.join(PD, "fetch.json"), "w"), indent=1)
    nb = (len(accs) + BATCH - 1) // BATCH
    t0 = time.time()
    for k in range(nb):
        path = os.path.join(RAW, f"bs{k:04d}.xml.gz")
        if os.path.exists(path):
            continue
        chunk = accs[k * BATCH:(k + 1) * BATCH]
        ids = search_ids("biosample", chunk)
        x = eutil("efetch.fcgi", {"db": "biosample", "id": b",".join(ids).decode(), "retmode": "xml"}) if ids else b""
        with gzip.open(path + ".tmp", "wb") as fo:
            fo.write(x)
        os.replace(path + ".tmp", path)
        if k % 20 == 0:
            print(f"batch {k + 1}/{nb}, {len(ids)} records, {time.time() - t0:.0f} s", flush=True)
    print("fetch done", flush=True)


def parse():
    G = genomes()
    bs = biosample_of()
    by_bs = collections.defaultdict(list)
    for g in G:
        by_bs[bs.get(g["accession"], "")].append(g["accession"])
    rec = {}
    names = collections.Counter()
    for f in sorted(os.listdir(RAW)):
        if not f.endswith(".xml.gz"):
            continue
        x = gzip.open(os.path.join(RAW, f)).read()
        if not x.strip():
            continue
        for el in ET.fromstring(x).iter("BioSample"):
            acc = el.get("accession")
            src, how = set(), set()
            for at in el.iter("Attribute"):
                nm = (at.get("attribute_name") or "") + "|" + (at.get("harmonized_name") or "")
                if "deriv" in nm.lower() or "source material" in nm.lower():
                    t = at.text or ""
                    found = set(BS_RE.findall(t)) | set(SRA_RE.findall(t))
                    if found:
                        names[at.get("attribute_name")] += 1
                        how.add(at.get("attribute_name"))
                    src |= found
            src.discard(acc)
            rec[acc] = (sorted(src), sorted(how))
    with open(os.path.join(PD, "derived.tsv"), "w") as fo:
        n = nsrc = 0
        for b, gl in by_bs.items():
            if b not in rec:
                continue
            s, how = rec[b]
            for a in gl:
                fo.write(f"{a}\t{b}\t{','.join(s)}\t{';'.join(how)}\n")
                n += 1
                nsrc += bool(s)
    print(json.dumps({"records": len(rec), "genomes_with_record": n, "with_derived_from": nsrc,
                      "attribute_names": names.most_common(10)}, indent=1))


def resolve():
    sra = set()
    for l in open(os.path.join(PD, "derived.tsv")):
        a, b, s, how = l.rstrip("\n").split("\t")
        sra |= {x for x in s.split(",") if SRA_RE.fullmatch(x)}
    path = os.path.join(PD, "sra_map.tsv")
    done = set()
    if os.path.exists(path):
        done = {l.split("\t")[0] for l in open(path)}
    todo = sorted(sra - done)
    print(f"{len(sra)} SRA accessions, {len(todo)} to resolve", flush=True)
    t0 = time.time()
    with open(path, "a") as fo:
        for k in range(0, len(todo), BATCH):
            chunk = todo[k:k + BATCH]
            ids = search_ids("sra", chunk)
            got = collections.defaultdict(set)
            for j in range(0, len(ids), 500):
                x = eutil("efetch.fcgi", {"db": "sra", "id": b",".join(ids[j:j + 500]).decode(),
                                          "rettype": "runinfo", "retmode": "text"})
                for r in csv.DictReader(io.StringIO(x.decode("utf-8", "replace"))):
                    b = r.get("BioSample", "")
                    if not BS_RE.fullmatch(b or ""):
                        continue
                    for col in ("Run", "Experiment", "Sample"):
                        if r.get(col):
                            got[r[col]].add(b)
            for a in chunk:
                fo.write(f"{a}\t{','.join(sorted(got.get(a, ())))}\n")
            fo.flush()
            if (k // BATCH) % 20 == 0:
                print(f"{k + len(chunk)}/{len(todo)}, {time.time() - t0:.0f} s", flush=True)
    print("resolve done", flush=True)


def resolve_samples():
    """SRA sample accessions without runs at NCBI (mostly ENA ERS) -> BioSample, by the record's own SRA id"""
    todo = []
    for l in open(os.path.join(PD, "sra_map.tsv")):
        a, b = l.rstrip("\n").split("\t")
        if not b and re.fullmatch(r"[SED]RS\d+", a):
            todo.append(a)
    path = os.path.join(PD, "sample_map.tsv")
    done = {l.split("\t")[0] for l in open(path)} if os.path.exists(path) else set()
    todo = sorted(set(todo) - done)
    print(f"{len(todo)} sample accessions to resolve", flush=True)
    with open(path, "a") as fo:
        for k in range(0, len(todo), 100):
            chunk = todo[k:k + 100]
            ids = search_ids("biosample", chunk)
            got = collections.defaultdict(set)
            if ids:
                x = eutil("efetch.fcgi", {"db": "biosample", "id": b",".join(ids).decode(), "retmode": "xml"})
                for el in ET.fromstring(x).iter("BioSample"):
                    for i in el.iter("Id"):
                        if i.get("db") == "SRA" and i.text:
                            got[i.text.strip()].add(el.get("accession"))
            for a in chunk:
                fo.write(f"{a}\t{','.join(sorted(got.get(a, ())))}\n")
            fo.flush()
    print("resolve_samples done", flush=True)


def sources():
    smap = {}
    for f in ("sra_map.tsv", "sample_map.tsv"):
        if os.path.exists(os.path.join(PD, f)):
            for l in open(os.path.join(PD, f)):
                a, b = l.rstrip("\n").split("\t")
                b = [x for x in b.split(",") if x]
                if b or a not in smap:
                    smap[a] = b
    G = {g["accession"]: g for g in genomes()}
    st = collections.Counter()
    per_src = collections.Counter()
    nsrc = collections.Counter()
    with open(os.path.join(PD, "sources.tsv"), "w") as fo:
        for l in open(os.path.join(PD, "derived.tsv")):
            a, b, s, how = l.rstrip("\n").split("\t")
            out, unres = set(), 0
            for x in (s.split(",") if s else []):
                if BS_RE.fullmatch(x):
                    out.add(x)
                elif x in smap and smap[x]:
                    out.update(smap[x])
                else:
                    out.add(x)  # unresolved: the SRA accession itself is still a sample key shared between MAGs
                    unres += 1
            out.discard(b)
            lab = "hab" if G[a]["habitat"] in HABS else "nohab"
            st[(lab, "record")] += 1
            st[(lab, "with_source")] += bool(out)
            st[(lab, "with_unresolved_key")] += unres > 0
            if out:
                nsrc[min(len(out), 10)] += 1
                for x in out:
                    per_src[x] += 1
                fo.write(f"{a}\t{','.join(sorted(out))}\n")
    G_hab = sum(1 for g in G.values() if g["habitat"] in HABS)
    res = {"genomes": len(G), "habitat_labelled": G_hab,
           "stats": {f"{k[0]}|{k[1]}": v for k, v in sorted(st.items())},
           "sources_per_genome(10=10+)": dict(sorted(nsrc.items())),
           "source_samples": len(per_src),
           "genomes_per_source_top": per_src.most_common(10),
           "sources_with_2plus_genomes": sum(1 for v in per_src.values() if v > 1)}
    json.dump(res, open(os.path.join(ROOT, "results_hgt", "provenance.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


def main():
    st = sys.argv[1] if len(sys.argv) > 1 else ""
    f = {"fetch": fetch, "parse": parse, "resolve": resolve, "resolve_samples": resolve_samples,
         "sources": sources}.get(st)
    if not f:
        raise SystemExit(__doc__)
    os.makedirs(PD, exist_ok=True)
    f()


if __name__ == "__main__":
    main()
