#!/usr/bin/env bash
# Query throughput on a 2,800-contig set against E. coli (median of 3; single thread).
set -euo pipefail
cd "$(dirname "$0")/.."
SKANI=${SKANI:-/home/mrahmanh/miniconda3/envs/meta-co-assembler-env/bin/skani}
[ -s data/throughput/NC_000913.3.contigs.fa ] || /usr/bin/python3 bench/make_scenarios.py data/genomes data/throughput 100 4242
q=data/throughput/NC_000913.3.contigs.fa; ref=data/genomes/NC_000913.3.fa
T=$(mktemp -d); out=data/results/main/throughput.tsv
n=$(grep -c '>' $q); bp=$(grep -v '>' $q | tr -d '\n' | wc -c)
med() { sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}'; }
run3() { : > $T/t; : > $T/m; for i in 1 2 3; do /usr/bin/time -f "%e %M" -o $T/r "$@" > /dev/null 2>&1; awk '{print $1}' $T/r >> $T/t; awk '{print $2}' $T/r >> $T/m; done; echo "$(med < $T/t) $(sort -n $T/m | tail -1)"; }
bin/aif sketch $ref -o $T/r.aif 2>/dev/null
read sk_t sk_m < <(run3 $SKANI dist -t 1 --qi -q $q -r $ref -o $T/o)
read as_t as_m < <(run3 bin/aif sketch $ref -o $T/r2.aif)
read aq_t aq_m < <(run3 bin/aif query $T/r.aif $q)
printf "tool\tphase\tcontigs\tbases\tmedian_s\tmax_rss_kb\n" > $out
printf "skani\tsketch+query\t%s\t%s\t%s\t%s\n" $n $bp $sk_t $sk_m >> $out
printf "aif\tsketch\t%s\t%s\t%s\t%s\n" $n $bp $as_t $as_m >> $out
printf "aif\tquery\t%s\t%s\t%s\t%s\n" $n $bp $aq_t $aq_m >> $out
rm -rf $T; column -t $out
