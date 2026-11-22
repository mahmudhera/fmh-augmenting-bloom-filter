#!/usr/bin/env bash
# Download the real-data genomes chosen by select_datasets.py.
#   NCBI assemblies  -> data/real/raw/ncbi/<accession>.fna   (NCBI datasets CLI)
#   bin-2328 MAGs    -> data/real/raw/mag2328/<name>.fa       (github.com/bluenote-1577/skani-test)
#   references       -> data/real/refs/{BA52G,MG1655,PAO1,NCTC8325}.fa
set -euo pipefail
cd "$(dirname "$0")/.."
DATASETS=${DATASETS:-/home/mrahmanh/miniconda3/envs/meta-sim/bin/datasets}
R=${REAL_DIR:-data/real}
mkdir -p $R/raw/ncbi $R/raw/mag2328 $R/refs $R/tmp

# references: B. anthracis 52-G chromosome (the D2 reference), and the three simulation references
[ -s $R/refs/BA52G.fa ] || curl -sS -m 300 \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CM002395.1&rettype=fasta&retmode=text" \
  -o $R/refs/BA52G.fa
cp --update=none data/genomes/NC_000913.3.fa $R/refs/MG1655.fa
cp --update=none data/genomes/NC_002516.2.fa $R/refs/PAO1.fa
cp --update=none data/genomes/NC_007795.1.fa $R/refs/NCTC8325.fa

# NCBI assemblies, one batch download
for s in bcereus ecoli paer saur negative; do if [ -f $R/lists/$s.tsv ]; then tail -n +2 $R/lists/$s.tsv; fi; done | awk -F'\t' '$3=="ncbi"{print $4}' | sort -u > $R/tmp/accs.txt
todo=$R/tmp/todo.txt; : > $todo
while read -r a; do ls $R/raw/ncbi/"$a"*.fna >/dev/null 2>&1 || echo "$a" >> $todo; done < $R/tmp/accs.txt
if [ -s $todo ]; then
  "$DATASETS" download genome accession --inputfile $todo --include genome --filename $R/tmp/ncbi.zip --no-progressbar
  rm -rf $R/tmp/ncbi && mkdir -p $R/tmp/ncbi && (cd $R/tmp/ncbi && unzip -q ../ncbi.zip)
  for d in $R/tmp/ncbi/ncbi_dataset/data/GC*; do
    acc=$(basename "$d"); f=$(ls "$d"/*.fna | head -1)
    cp "$f" $R/raw/ncbi/"$acc".fna
  done
fi
echo "NCBI assemblies: $(ls $R/raw/ncbi | wc -l) of $(wc -l < $R/tmp/accs.txt) requested"

# bin-2328 MAGs from the skani-test repository
tail -n +2 $R/lists/mag2328.tsv | cut -f4 | while read -r p; do
  o=$R/raw/mag2328/$(basename "$p")
  [ -s "$o" ] || curl -sS -m 120 -o "$o" "https://raw.githubusercontent.com/bluenote-1577/skani-test/main/$p"
done
echo "MAGs: $(ls $R/raw/mag2328 | wc -l)"
