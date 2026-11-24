#!/usr/bin/env bash
# Download every genome of the UHGG v3.0 human-gut catalogue (MGnify, 286,621 genomes; Prokka GFF3 with the
# assembly in its ##FASTA section, about 240 GB of gff.gz) for the contig-level contamination screen.
#   bench_real/uhgg_download.sh            (UHGG_DIR, default data/uhgg3; UHGG_JOBS parallel connections, default 16)
# Each connection fetches 50 genomes (connection reuse; one TLS handshake per file trips the server's rate limit
# at 32 connections). An interrupted transfer can leave a truncated file: uhgg_prepare.py checks every gzip.
# Metadata (genomes-all_metadata.tsv, added_genomes.tsv, removed_genomes.tsv, README) go to $UHGG_DIR/meta.
# Re-running skips genomes already downloaded.
set -euo pipefail
cd "$(dirname "$0")/.."
D=${UHGG_DIR:-data/uhgg3}
J=${UHGG_JOBS:-16}
B=https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v3.0
mkdir -p "$D/meta" "$D/gff"
for f in README_v3.0.txt known_issues.txt study_list.txt added_genomes.tsv removed_genomes.tsv genomes-all_metadata.tsv; do
  [ -s "$D/meta/$f" ] || curl -s -f --retry 5 -o "$D/meta/$f" "$B/$f"
done
python3 -I - "$D" <<'PY' > "$D/urls.txt"
import csv, os, sys
d = sys.argv[1]
for r in csv.DictReader(open(os.path.join(d, "meta", "genomes-all_metadata.tsv")), delimiter="\t"):
    if not os.path.exists(os.path.join(d, "gff", r["Genome"] + ".gff.gz")):
        print(r["FTP_download"].replace("ftp://", "https://"))
PY
echo "$(wc -l < "$D/urls.txt") genomes to fetch"
[ -s "$D/urls.txt" ] || exit 0
cd "$D/gff"
split -l 50 ../urls.txt ../batch_
fetch() {
  args=()
  while read -r u; do args+=(-o "$(basename "$u")" "$u"); done < "$1"
  curl -s -f --retry 8 --retry-delay 10 --retry-all-errors --remove-on-error "${args[@]}" || echo "batch $1: curl exit $?"
  rm -f "$1"
}
export -f fetch
ls ../batch_* | xargs -P "$J" -I{} bash -c 'fetch "$@"' _ {}
echo "done: $(ls | grep -c '\.gff\.gz$') files"
