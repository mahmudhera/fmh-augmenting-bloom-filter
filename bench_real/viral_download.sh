#!/usr/bin/env bash
# Gut Phage Database (Camarillo-Guerrero et al. 2021): sequences and metadata from EBI, in parallel
# byte ranges (single streams are slow). Output: data/viral/raw/
set -uo pipefail
cd "$(dirname "$0")/.."
B=https://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database
D=data/viral/raw
mkdir -p $D/parts
pget() {  # url out
  local url=$1 out=$2
  local size=$(curl -sSIL "$url" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1)
  local ch=$(( (size + 31) / 32 ))
  rm -f $D/parts/*
  seq 0 31 | xargs -P 16 -I{} bash -c 'a=$(( {} * '$ch' )); b=$(( a + '$ch' - 1 )); [ $b -ge '$size' ] && b='$((size-1))'; [ $a -lt '$size' ] && curl -sS -L --retry 8 -r $a-$b -o '$D'/parts/$(printf %02d {}) "'$url'"'
  cat $D/parts/* > $out
  [ $(stat -c %s $out) -eq $size ] && echo "$out ok ($size bytes)" || echo "$out SIZE MISMATCH"
  rm -f $D/parts/*
}
pget $B/GPD_metadata.tsv $D/GPD_metadata.tsv
pget $B/README.txt $D/README.txt
pget $B/GPD_sequences.fa.gz $D/GPD_sequences.fa.gz
gzip -t $D/GPD_sequences.fa.gz && echo "gzip integrity ok"
# INPHARED (Cook et al. 2021): complete phage genomes from GenBank (excluding RefSeq duplicates), with metadata
I=https://millardlab-inphared.s3.climb.ac.uk
mkdir -p data/viral/inphared
for f in 14Apr2025_genomes_excluding_refseq.fa.gz 14Apr2025_data_excluding_refseq.tsv.gz; do
  [ -s data/viral/inphared/$f ] || curl -sS -L --retry 5 -o data/viral/inphared/$f $I/$f
  gzip -t data/viral/inphared/$f && echo "$f ok"
done
