#!/usr/bin/env bash
# Download the four reference genomes used by the benchmark from NCBI (nuccore efetch).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/genomes
for acc in NC_000913.3 NC_000964.3 NC_002516.2 NC_007795.1; do
  [ -s data/genomes/$acc.fa ] && continue
  curl -sS -m 300 -o data/genomes/$acc.fa \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text"
  sleep 1
done
