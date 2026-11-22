#!/usr/bin/env bash
# Download the database-mode benchmark genomes (bench_real/db_select.py) with the NCBI datasets
# CLI (dehydrated package + rehydrate) into data/db/genomes/<accession>.fna, add the five
# pairwise references, and write data/db/meta/refs.tsv (path, name) for aif db build / skani.
set -euo pipefail
cd "$(dirname "$0")/.."
DATASETS=${DATASETS:-/home/mrahmanh/miniconda3/envs/meta-sim/bin/datasets}
B=data/db
mkdir -p $B/genomes $B/tmp
tail -n +2 $B/meta/db_accessions.tsv | cut -f1 > $B/tmp/all.txt
: > $B/tmp/todo.txt
while read -r a; do [ -s $B/genomes/$a.fna ] || echo "$a" >> $B/tmp/todo.txt; done < $B/tmp/all.txt
if [ -s $B/tmp/todo.txt ]; then
  rm -rf $B/tmp/pkg && mkdir -p $B/tmp/pkg
  "$DATASETS" download genome accession --inputfile $B/tmp/todo.txt --include genome --dehydrated \
    --filename $B/tmp/pkg/pkg.zip --no-progressbar
  (cd $B/tmp/pkg && unzip -q -o pkg.zip)
  "$DATASETS" rehydrate --directory $B/tmp/pkg --no-progressbar --max-workers 16
  for d in $B/tmp/pkg/ncbi_dataset/data/GC*; do
    acc=$(basename "$d"); f=$(ls "$d"/*.fna 2>/dev/null | head -1 || true)
    [ -n "$f" ] && mv "$f" $B/genomes/"$acc".fna
  done
  rm -rf $B/tmp/pkg
fi
for r in BA52G MG1655 PAO1 NCTC8325 REFMAG; do cp --update=none data/real/refs/$r.fa $B/genomes/$r.fna; done
{ for r in BA52G MG1655 PAO1 NCTC8325 REFMAG; do echo -e "$PWD/$B/genomes/$r.fna\t$r"; done
  while read -r a; do [ -s $B/genomes/$a.fna ] && echo -e "$PWD/$B/genomes/$a.fna\t$a"; done < $B/tmp/all.txt; } > $B/meta/refs.tsv
echo "genomes: $(wc -l < $B/meta/refs.tsv) of $(( $(wc -l < $B/tmp/all.txt) + 5 ))"
