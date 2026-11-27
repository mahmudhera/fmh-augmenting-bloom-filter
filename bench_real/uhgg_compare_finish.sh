#!/bin/bash
# What ran after uhgg_timed_runs.sh (2026-11-26/27): (1) CAT S1-S3 again (the first attempt passed CAT_pack 7 a
# contigs --force it no longer has); (2) MUMmer truth for every tool's new candidates (untimed, all cores), LexicMap's
# from S1; (3) evaluate S3 (all tools but LexicMap) and S1 (all tools); (4) a timed minimap2 pool run, because on
# S1-S3 its mapping time is per-index-part overhead, flat in query size, sampling its --split-prefix temp size;
# (5) project.
cd "$(dirname "$0")/.." || exit 1
L=data/compare/runs/timed.log
PY="/usr/bin/python3 -I bench_real/uhgg_compare.py"
step(){ echo "$(date -u +%FT%TZ) start $*" >> $L; "$@" >> $L 2>&1; local rc=$?; echo "$(date -u +%FT%TZ) end rc=$rc $*" >> $L; return $rc; }
export CPUS=0-47 THREADS=48
for s in S1 S2 S3; do step $PY run $s cat; done
unset CPUS THREADS
step env JOBS=60 $PY truth S3
step env JOBS=60 HITS_FROM=S1 TOOLS=lexicmap $PY truth S3
step $PY evaluate S3
step $PY evaluate S1
step $PY project
export CPUS=0-47 THREADS=48
D=data/compare/runs/pool; T=$D/mm2split_du.tsv; : > $T
( while :; do b=$(find $D -maxdepth 1 -name 'mm2split*' -printf '%s\n' 2>/dev/null | awk '{s+=$1} END {print s+0}'); echo "$(date +%s)	$b" >> $T; sleep 60; done ) &
SAMP=$!
step $PY run pool minimap2
kill $SAMP
echo "$(date -u +%FT%TZ) mm2split peak bytes $(sort -k2,2n $T | tail -1 | cut -f2)" >> $L
step $PY project
