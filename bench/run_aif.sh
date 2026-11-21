#!/usr/bin/env bash
# Sketch each reference with aif and query its contig set; record time/memory/sketch size.
# Usage: bench/run_aif.sh <outdir> "<sketch args>" "<query args>" [bench_dir]
set -euo pipefail
AIF=${AIF:-bin/aif}
out=$1; sargs=${2:-}; qargs=${3:-}; bench=${4:-data/bench}
mkdir -p "$out"
printf "genome\tsketch_s\tsketch_kb\tquery_s\tquery_kb\tsketch_bytes\tbases\n" > "$out/resources.tsv"
for ref in data/genomes/*.fa; do
  g=$(basename "$ref" .fa)
  /usr/bin/time -f "%e\t%M" -o "$out/$g.stime" $AIF sketch "$ref" -o "$out/$g.aif" $sargs 2> "$out/$g.slog"
  /usr/bin/time -f "%e\t%M" -o "$out/$g.qtime" $AIF query "$out/$g.aif" "$bench/$g.contigs.fa" $qargs > "$out/$g.tsv"
  bases=$(grep -v '>' "$ref" | tr -d '\n' | wc -c)
  printf "%s\t%s\t%s\t%s\t%s\n" "$g" "$(cat "$out/$g.stime")" "$(cat "$out/$g.qtime")" "$(stat -c %s "$out/$g.aif")" "$bases" >> "$out/resources.tsv"
done
{ printf "contig\tani\taf\n"; for f in "$out"/NC_*.tsv; do awk -F'\t' 'NR>1 && $8=="ok"{print $1"\t"$2"\t"$3}' "$f"; done; } > "$out/predictions.tsv"
