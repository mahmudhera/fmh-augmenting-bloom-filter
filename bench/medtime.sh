#!/usr/bin/env bash
# medtime.sh N cmd... : median wall seconds and max RSS (KB) over N runs after 3 warm-up runs
# (stdout of cmd discarded). For tool-vs-tool comparisons prefer bench/time_all.py (interleaved).
N=$1; shift
T=$(mktemp); M=$(mktemp); X=$(mktemp)
for i in 1 2 3; do "$@" > /dev/null 2>/dev/null; done  # warm-up runs (CPU frequency ramp)
for i in $(seq $N); do
  /usr/bin/time -f "%e %M" -o $X "$@" > /dev/null 2>/dev/null
  awk '{print $1}' $X >> $T; awk '{print $2}' $X >> $M
done
echo "$(sort -n $T | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}') s  $(sort -n $M | tail -1) KB"
rm -f $T $M $X
