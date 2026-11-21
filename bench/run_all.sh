#!/usr/bin/env bash
# Full comparison on the main benchmark: skani presets, aif final, ablations, budget sweep.
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/results/main
mkdir -p $R
PY=${PY:-/usr/bin/python3}
ev() { $PY bench/evaluate.py data/bench/truth.tsv skani=data/results/main/skani/predictions.tsv "$1=$R/$1/predictions.tsv" --json $R/$1/eval.json --quiet; }
# skani presets (sketch size measured separately)
bench/run_skani.sh $R/skani >/dev/null
bench/run_skani.sh $R/skani_medium --medium >/dev/null
bench/run_skani.sh $R/skani_slow --slow >/dev/null
for t in skani skani_medium skani_slow; do ev $t; done
# aif final + ablations
declare -A S Q
S[aif]="";                                   Q[aif]=""
S[abl_k17]="-k 17";                          Q[abl_k17]=""
S[abl_k21]="-k 21";                          Q[abl_k21]=""
S[abl_pooled]="";                            Q[abl_pooled]="--het 0"
S[abl_single_tier]="--tiers 0.01:1";         Q[abl_single_tier]=""
S[abl_nobias]="";                            Q[abl_nobias]="--no-bias-correct"
S[abl_v1_k21_hedge]="-k 21 --tiers 0.003:0.55,0.58:0.45"; Q[abl_v1_k21_hedge]=""
for b in 0.064 0.22 0.44 1.76; do S[budget_$b]="--budget $b"; Q[budget_$b]=""; done
for name in "${!S[@]}"; do
  bench/run_aif.sh $R/$name "${S[$name]}" "${Q[$name]}" >/dev/null 2>&1 &
done
wait
for name in "${!S[@]}"; do ev $name; done
echo done
