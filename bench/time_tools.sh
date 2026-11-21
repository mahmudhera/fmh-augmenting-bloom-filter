#!/usr/bin/env bash
# Median-of-N wall time and peak RSS, single thread, sequential runs.
# Usage: bench/time_tools.sh [N]   -> data/results/main/timing.tsv
set -euo pipefail
cd "$(dirname "$0")/.."
N=${1:-5}
SKANI=${SKANI:-/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani}
T=$(mktemp -d)
out=data/results/main/timing.tsv
printf "tool\tgenome\tphase\tmedian_s\tmax_rss_kb\n" > $out
med() { sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}'; }
for ref in data/genomes/*.fa; do
  g=$(basename $ref .fa); q=data/bench/$g.contigs.fa
  for phase in skani_total aif_sketch aif_query; do
    : > $T/t; : > $T/m
    for i in $(seq $N); do
      case $phase in
        skani_total) /usr/bin/time -f "%e %M" -o $T/r $SKANI dist -t 1 --qi -q $q -r $ref -o $T/o >/dev/null 2>&1 ;;
        aif_sketch)  /usr/bin/time -f "%e %M" -o $T/r bin/aif sketch $ref -o $T/s.aif >/dev/null 2>&1 ;;
        aif_query)   /usr/bin/time -f "%e %M" -o $T/r bin/aif query $T/s.aif $q > $T/o 2>/dev/null ;;
      esac
      awk '{print $1}' $T/r >> $T/t; awk '{print $2}' $T/r >> $T/m
    done
    printf "%s\t%s\t%s\t%s\t%s\n" "${phase%%_*}" "$g" "${phase#*_}" "$(med < $T/t)" "$(sort -n $T/m | tail -1)" >> $out
  done
done
rm -rf $T
column -t $out
