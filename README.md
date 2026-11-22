# FracMinHash + Bloom filters for ANI / AF on short contigs

Three reports:

1. **Theory**: [`fmh_bloom_filter_theory.pdf`](fmh_bloom_filter_theory.pdf) (27 pages). Variance framework, 14 designs, 16 scenarios.
2. **Implementation + simulated benchmark**: [`aif_benchmark_report.pdf`](aif_benchmark_report.pdf) (15 pages). The `aif` tool (v1.1) compared with skani 0.4.0 at equal sketch budget, including the v1.1 speed rewrite.
3. **Real genomes**: [`aif_real_genomes_report.pdf`](aif_real_genomes_report.pdf) (18 pages). 240 real assemblies and MAGs from the skani paper's data sets, with MUMmer ground truth; the v1.2/v1.3 speed work; v1.4 calibrated confidence intervals; and an analysis of the AF model.

## Real genomes (skani paper data; 35,362 contigs >= 500 bp; MUMmer truth; no aif tuning on these data)

Sets: B. cereus group drafts (D2) vs B. anthracis 52-G; E. coli drafts (D3) vs K-12; 48 MAGs of Pasolli bin 2328 vs
the bin's best MAG; P. aeruginosa and S. aureus drafts; cross-genus negative controls. Pooled over the 33,464
within-species contigs:

| | skani 0.4.0 | skani `--slow` | **aif** | aif-RT (`--rho-u 2e-5`) |
|---|---|---|---|---|
| sketch size (bits / reference base) | 0.903 | 3.33 | 0.880 | 0.880 |
| homologous contigs reported | 90.1% | 95.6% | **99.2%** | 99.3% |
| ANI MAE, reported homologous contigs (pp) | 0.677 | 0.611 | **0.502** | 0.496 |
| ANI MAE, contigs reported by all (pp) | 0.675 | 0.532 | **0.392** | 0.387 |
| AF MAE, all contigs (pp) | 9.35 | 4.29 | **2.63** | 1.85 |
| AF MAE on 0.5-1 kb contigs (pp) | 39.2 | 21.3 | **4.3** | 3.9 |
| genome ANI MAE vs ANIm, 224 genomes (pp) | 0.237 | 0.259 | **0.190** | 0.263 |
| wall time, 5 sets, 1 / 8 threads (s), aif v1.4 | 3.27 / 2.72 | | **2.92 / 0.86** | |

aif beats both skani presets on ANI and AF on every set. aif-RT is one setting chosen on a disjoint real tuning set
(it helps contig AF, costs genome-level ANI on divergent genomes and simulated AF; the default is unchanged).
fastANI's genome ANI MAE on the same genomes is 0.247 pp. Broken down by set, length, divergence, aligned
fraction and genome (`bench_real/consistency.py`), aif's paired ANI error is lower than skani's in every cell.

**Confidence intervals (v1.4).** Every estimate has an ANI and an AF interval (`ani_lo ani_hi af_lo af_hi`).
- **Statistical intervals** (`--ci-stat`) come from the estimator's variance: Theorem 1, conditioned on the
  realised substitutions, for ANI; posterior boundary uncertainty for AF. They are calibrated when truth
  follows the model: 96.5% coverage at 95% on simulated substitution-only contigs.
- **Default intervals** are calibrated against MUMmer truth by stratified split-conformal prediction, fitted on
  the real tuning set only. On held-out real contigs they cover 94.4% (ANI) and 94.2% (AF) at 95%, and
  98.6% / 99.0% at 99%. Median half-widths are 0.87 pp (ANI) and 2.6 pp (AF); at >= 99% identity on contigs
  >= 3 kb they are 0.44 / 0.54 pp. skani reports no uncertainty.

**AF model (v1.4 analysis).** aif over-estimates AF on partially aligned, divergent contigs because unaligned
islands get absorbed as low-identity alignment.
- **Local identity prior** (`--ident-window 0.08 --ident-radius 3`, an option): lowers held-out contig AF error
  by 13% (2.45 -> 2.14 pp) but raises genome-level ANI error (0.190 -> 0.251 pp).
- **Rejected:** model averaging over switching rates, and overdispersed emissions.
- **Conclusion:** no tested correction improves both contig AF and genome ANI, so the defaults are unchanged.

## Headline benchmark result (held-out set, 2,800 contigs, 28 conditions, 4 bacterial references)

| | skani 0.4.0 (default) | **aif v1.1** |
|---|---|---|
| sketch size (bits / reference base) | 0.900 | 0.880 |
| contigs reported | 86.8% | **95.1%** |
| ANI MAE, all contigs / mean over conditions (pp) | 1.117 / 1.448 | **0.357 / 0.502** |
| AF MAE, all contigs / mean over conditions (pp) | 7.43 / 9.98 | **0.39 / 0.51** |
| conditions where aif has lower ANI and AF error | | **28 / 28** |
| peak memory | 19.0 MB | **15.5 MB** |

A freshly generated confirmation benchmark reproduces this: ANI 0.366 vs 1.206 pp, AF 0.55 vs 7.46 pp, 28/28 conditions.
aif at half the budget (0.44 bits/base) still beats skani's densest preset (`--slow`, 3.33 bits/base).

## Speed (single thread unless noted; medians of interleaved, core-pinned runs)

Real genomes, aif v1.4 (with confidence intervals) sketch + query vs skani `dist --qi` (sketch + query):

| set (contigs vs reference) | Mbp | skani 1 thread | **aif 1 thread** | skani 8 threads | **aif 8 threads** |
|---|---|---|---|---|---|
| B. cereus group vs B. anthracis 52-G | 449 | 1.34 s | **1.33 s** | 1.16 s | **0.37 s** |
| E. coli vs K-12 | 243 | 0.75 s | **0.66 s** | 0.62 s | **0.17 s** |
| MAGs of bin 2328 vs best MAG | 116 | 0.39 s | **0.29 s** | 0.31 s | **0.09 s** |
| P. aeruginosa vs PAO1 | 165 | 0.56 s | **0.45 s** | 0.44 s | **0.14 s** |
| S. aureus vs NCTC 8325 | 68 | 0.23 s | **0.18 s** | 0.19 s | **0.08 s** |

v1.3 (no intervals) was 4% faster (B. cereus 1.28 s).

Simulated contigs:

| workload | skani (sketch + query) | aif v1.0 | aif v1.1 | **aif v1.3** |
|---|---|---|---|---|
| 20 Mbp of contigs vs E. coli, sketch + query | 0.095 s | 0.373 s | 0.084 s | **0.067 s** |
| 81 Mbp of contigs vs E. coli, sketch + query | 0.300 s | 1.319 s | 0.259 s | **0.194 s** |
| 81 Mbp, 8 threads (aif: query only) | 0.452 s | n/a | 0.122 s | **0.062 s** |
| 4 references x 700 contigs, total | **0.138 s** | 0.545 s | 0.178 s | 0.155 s |

v1.1 queries are 5.2x faster than v1.0 and sketching is 1.9x faster, at equal accuracy. skani is still faster when a
few hundred contigs meet a newly sketched reference (aif's sketch build costs more); aif is faster on larger
query sets and scales with threads. The fast paths need AVX-512; a scalar fallback gives byte-identical sketches.

**v1.2 (real assemblies).** Draft assemblies have Mbp-scale contigs, on which v1.1's batch-and-barrier thread
pipeline did not scale (449 Mbp B. cereus set: 1.9 s at 1 thread, 3.5 s at 8 threads, 1.6 GB). v1.2 probes long
contigs in cache-sized windows and uses a streaming pipeline (chunked work queue, ordered writer, workers strip FASTA
line breaks from the memory map), with output byte-identical to v1.1. `--het-fine-max 16384` is an optional switch,
about 20% faster on Mbp-scale contigs, with slightly less accurate genome-level ANI.

**v1.3 (single thread).** After v1.2, aif was still up to 31% slower than skani at 1 thread on real sets. v1.3 is
29% faster than v1.2 there (B. cereus set: 1.19 s at 1 thread, 0.34 s at 8, 0.29 s at 16). Changes:
- HMM forward/backward recursions are interleaved, and their per-step normalisers come from look-ahead scalars, so
  horizontal reductions leave the dependency chain (same for per-probe refinement).
- Fine bins use cached emission vectors.
- Reverse complements use GFNI (`vgf2p8affineqb`).
- AVX-512 2-bit packing; mask-popcount bin counts; batched, branch-free calibration.

Sketches are byte-identical. Estimates differ from v1.2 only by floating-point rounding (at most 1e-6, in 0.3% of rows),
and every accuracy metric is unchanged. When a few hundred contigs meet a newly sketched reference, the sketch build
(18-40 ms, built twice for GC-rich references by the adaptive-k search) still makes skani slightly faster end to end.
Simulated-set timings are in `results/main/*_v12.tsv` and `*_v13.tsv`.

## The `aif` tool

`src/aif.cpp` is a single dependency-free C++17 file.

* **Sketch.** Canonical k-mers are FracMinHash-selected into two hash-range tiers, each stored as a blocked Bloom filter (FPR 0.005 and 0.1; 0.88 bits/base by default). The exact FPR comes from the fill ratio. The spurious-match rate is calibrated empirically, and k is chosen per reference (17-21) to keep it at or below 0.003.
* **Query** (raw contigs). Every selected k-mer is probed. A multi-level HMM (one unaligned state plus aligned identity levels) gives the aligned fraction and the FP-corrected aligned containment, so ANI = C^(1/k) with bias correction. Heterogeneous contigs use posterior-mean ANI instead.
* **v1.1 speed**: 2-bit packing and AVX-512 k-mer selection, gather-based Bloom tests, a binned single-precision HMM with exact per-probe refinement around alignment boundaries, mmap input, vectorised sketch insertion, and `-t N` threads with ordered output.
* **v1.2/v1.3**: windowed probing of long contigs, a streaming multithreaded pipeline, latency-optimised HMM recursions, GFNI/AVX-512 packing and reverse complements (see above); `make prof` builds a per-stage cycle-counter binary.
* **v1.4**: per-contig confidence intervals (statistical, or conformal-calibrated from the table in `src/ci_table.inc`, refit with `bench_real/calibration.py fit`); optional local identity prior (`--ident-window`) and overdispersed emissions (`--overdisp`).

```
make aif
bin/aif sketch ref.fa -o ref.aif          # -k auto --budget 0.88 --tiers 0.005:0.6,0.1:0.4
bin/aif query ref.aif contigs.fa -t 8 > out.tsv  # contig ani af theta probes hits aligned_probes status lr het low_frac
                                                #   ani_lo ani_hi af_lo af_hi ani_se af_se   (95% intervals; --ci-level, --ci-stat)
```

## Layout

```
fmh_bloom_filter_theory.pdf, aif_benchmark_report.pdf, aif_real_genomes_report.pdf
src/aif.cpp, src/ci_table.inc  the tool (+ built-in confidence-interval calibration table)
bench/                         simulated benchmark: download_genomes.sh, make_scenarios.py, run_skani.sh, run_aif.sh,
                               run_all.sh, evaluate.py, tune.py, time_all.py (interleaved timing), time_compare.sh, medtime.sh
bench_real/                    real genomes: select_datasets.py, select_tune.py, download.sh, prepare.py, truth.py (MUMmer),
                               run_tools.py, evaluate_real.py, genome_level.py, blast_audit.py, tune_real.py,
                               time_real.py, scaling.py, consistency.py, calibration.py, af_model_grid.py, afopt_heldout.py
results/, results_real/        committed summaries: evaluations (JSON), predictions, truth, timings, tuning tables
scripts/                       theory tables/figures, Monte Carlo check, report tables/figures (benchmark_report.py, real_report.py)
report/, report_bench/,        groff sources of the three PDFs (custom macros; built with groff + ghostscript)
report_real/
data/  (git-ignored)           genomes, assemblies, generated contigs, raw run outputs; regenerate with the scripts
```

## Reproduce

```
make aif && bench/download_genomes.sh
python3 bench/make_scenarios.py data/genomes data/bench 25            # main benchmark
python3 bench/make_scenarios.py data/genomes data/tune 10 777         # tuning set
python3 bench/make_scenarios.py data/genomes data/bench2 25 99991     # fresh confirmation set
bench/run_all.sh && python3 bench/time_all.py 11
python3 scripts/benchmark_report.py && report_bench/build.sh
```

Real genomes (needs the NCBI `datasets` CLI, MUMmer 4, blastn, fastANI and skani; paths are overridable with
`DATASETS=`, `MUMMER=`, `SKANI=`, `FASTANI=`):

```
python3 -I bench_real/select_datasets.py && bench_real/download.sh && python3 -I bench_real/prepare.py
python3 -I bench_real/truth.py && python3 -I bench_real/run_tools.py
python3 -I bench_real/evaluate_real.py --json data/real/results/eval_main.json
python3 -I bench_real/genome_level.py && python3 -I bench_real/blast_audit.py && python3 -I bench_real/time_real.py 7
python3 -I bench_real/consistency.py aif && python3 -I bench_real/calibration.py evaluate
python3 -I bench_real/af_model_grid.py && python3 -I bench_real/afopt_heldout.py
python3 scripts/real_report.py && report_real/build.sh
```

The scripts expect skani at the path in `bench/run_skani.sh` (override with `SKANI=...`), plus numpy and matplotlib for the Python parts.
