#!/usr/bin/env bash
# Full UHGG v3.0 screen after the download: fetch anything missing, extract (re-fetching truncated files), search,
# per-genome calls. Logs to data/uhgg3/run.log.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data/uhgg3
export THREADS=${THREADS:-48}
bench_real/uhgg_download.sh                      # second pass: only genomes still missing
python3 -I bench_real/uhgg_screen.py prepare
bad=$(python3 -I -c "import json; print(' '.join(g for g, _ in json.load(open('$D/prepare.json'))['bad_gzip']))")
if [ -n "$bad" ]; then
  echo "re-fetching truncated files: $bad"
  for g in $bad; do rm -f "$D/gff/$g.gff.gz"; done
  bench_real/uhgg_download.sh
  python3 -I bench_real/uhgg_screen.py prepare
fi
python3 -I bench_real/uhgg_screen.py search
python3 -I bench_real/uhgg_screen.py genomes
echo UHGG_DONE
