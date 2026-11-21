#!/usr/bin/env bash
# Interleaved timing (same core, alternating tools, warm-up first) so both tools see the
# same CPU frequency state. Usage: bench/time_compare.sh ROUNDS label1 "cmd1" label2 "cmd2" ...
R=$1; shift
labels=(); cmds=()
while [ $# -gt 0 ]; do labels+=("$1"); cmds+=("$2"); shift 2; done
T=$(mktemp -d)
for c in "${cmds[@]}"; do for i in 1 2 3; do bash -c "$c" > /dev/null 2>&1; done; done
for r in $(seq $R); do
  for i in "${!cmds[@]}"; do
    /usr/bin/time -f "%e %M" -o $T/x bash -c "${cmds[$i]}" > /dev/null 2>&1
    cat $T/x >> $T/t$i
  done
done
for i in "${!cmds[@]}"; do
  med=$(awk '{print $1}' $T/t$i | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  mn=$(awk '{print $1}' $T/t$i | sort -n | head -1)
  mem=$(awk '{print $2}' $T/t$i | sort -n | tail -1)
  printf "%-28s median %6.3f s  min %6.3f s  peak %6d KB\n" "${labels[$i]}" "$med" "$mn" "$mem"
done
rm -rf $T
