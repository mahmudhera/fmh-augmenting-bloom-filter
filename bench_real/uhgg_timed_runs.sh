#!/bin/bash
# Timed runs of the UHGG contig comparison on a quiet node (48 threads pinned to cores 0-47), tool-major so each
# database is warmed once. Waits for the sourmash index build. Every step is logged; a failed step does not stop
# the rest. LexicMap goes on to S2 only if S1 took <= 2 h, and to S3 only if S2 took <= 3 h (S3 has ~3.5x S2's bases);
# otherwise its larger subsets are left to a decision on the measured times.
cd "$(dirname "$0")/.." || exit 1
L=data/compare/runs/timed.log
mkdir -p data/compare/runs
export CPUS=0-47 THREADS=48
PY="/usr/bin/python3 -I bench_real/uhgg_compare.py"
building(){ for p in $(pgrep -u "$USER" -x python3); do tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE "uhgg_compare.py build " && return 0; done; return 1; }
step(){ echo "$(date -u +%FT%TZ) start $*" >> $L; "$@" >> $L 2>&1; local rc=$?; echo "$(date -u +%FT%TZ) end rc=$rc $*" >> $L; return $rc; }
secs(){ /usr/bin/python3 -I -c "import json,sys; print(int(json.load(open(sys.argv[1])).get(sys.argv[2],{}).get('s',10**9)))" \
        "data/compare/runs/$1/time.json" "$2" 2>/dev/null || echo 1000000000; }

while building; do sleep 60; done
echo "$(date -u +%FT%TZ) builds done; timed runs begin (CPUS=$CPUS THREADS=$THREADS)" >> $L

for t in aif skani kraken2; do
  for s in S1 S2 S3; do step $PY run $s $t; done
  step env RERUN=1 $PY run pool $t
done
if [ -f data/compare/build_sourmash.json ]; then
  for s in S1 S2 S3 pool; do step $PY run $s sourmash; done
else
  echo "$(date -u +%FT%TZ) skip sourmash: no build_sourmash.json" >> $L
fi
for s in S1 S2 S3; do step $PY run $s cat; done
step $PY run S1 lexicmap
if [ "$(secs S1 lexicmap)" -le 7200 ]; then
  step $PY run S2 lexicmap
  if [ "$(secs S2 lexicmap)" -le 10800 ]; then step $PY run S3 lexicmap; else echo "$(date -u +%FT%TZ) defer lexicmap S3" >> $L; fi
else
  echo "$(date -u +%FT%TZ) defer lexicmap S2/S3 (S1 took $(secs S1 lexicmap) s)" >> $L
fi
for s in S1 S2 S3; do step $PY run $s minimap2; done
echo "$(date -u +%FT%TZ) timed runs done" >> $L
