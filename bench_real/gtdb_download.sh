#!/usr/bin/env bash
# Download the GTDB species-representative genomes (latest release) in parallel 1 GB byte ranges
# from the Aalborg mirror, verify the MD5, and unpack to data/gtdb/genomes (one .fna.gz per genome).
set -uo pipefail
cd "$(dirname "$0")/.."
M=${GTDB_MIRROR:-https://data.gtdb.aau.ecogenomic.org/releases/latest}
U=$M/genomic_files_reps/gtdb_genomes_reps.tar.gz
D=data/gtdb
mkdir -p $D/parts
for f in VERSION.txt MD5SUM.txt bac120_metadata.tsv.gz ar53_metadata.tsv.gz bac120_taxonomy.tsv.gz ar53_taxonomy.tsv.gz; do
  [ -s $D/$f ] || curl -sS -L --retry 5 -o $D/$f "$M/$f"
done
SIZE=$(curl -sSIL "$U" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1)
CH=$((1<<30)); N=$(( (SIZE + CH - 1) / CH ))
echo "size $SIZE bytes, $N parts"
get() {  # part index
  local i=$1 a=$(( $1 * CH )) b=$(( ($1 + 1) * CH - 1 )); [ $b -ge $SIZE ] && b=$((SIZE - 1))
  local want=$((b - a + 1)) p=$D/parts/$(printf %05d $i)
  for t in 1 2 3 4 5 6; do
    [ -f $p ] && [ $(stat -c %s $p) -eq $want ] && return 0
    curl -sS -L --retry 5 -m 1800 -r $a-$b -o $p.tmp "$U" && mv $p.tmp $p
  done
  echo "part $i failed" >&2; return 1
}
export -f get; export D U CH SIZE
seq 0 $((N - 1)) | xargs -P 16 -I{} bash -c 'get {}'
cat $D/parts/* > $D/gtdb_genomes_reps.tar.gz
md5=$(md5sum $D/gtdb_genomes_reps.tar.gz | cut -d' ' -f1)
grep -q "$md5" $D/MD5SUM.txt && echo "MD5 ok $md5" || echo "MD5 MISMATCH $md5"
