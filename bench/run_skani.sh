#!/usr/bin/env bash
# Run skani (default parameters unless extra args given) on every genome's contig set.
# Usage: bench/run_skani.sh <outdir> [extra skani dist args...]
set -euo pipefail
SKANI=${SKANI:-/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani}
out=$1; shift
mkdir -p "$out"
: > "$out/resources.tsv"
bench=${BENCH:-data/bench}
for ref in data/genomes/*.fa; do
  g=$(basename "$ref" .fa)
  /usr/bin/time -f "%e\t%M" -o "$out/$g.time" \
    "$SKANI" dist -t 1 --qi -q "$bench/$g.contigs.fa" -r "$ref" -o "$out/$g.tsv" "$@" 2> "$out/$g.log"
  printf "%s\t%s\n" "$g" "$(cat "$out/$g.time")" >> "$out/resources.tsv"
done
# common format: contig, ani (fraction), af_query (fraction)
for f in "$out"/NC_*.tsv; do tail -n +2 "$f"; done | awk -F'\t' 'BEGIN{OFS="\t"; print "contig","ani","af"} {split($7,a," "); print a[1], $3/100, $5/100}' > "$out/predictions.tsv"
cat "$out/resources.tsv"
