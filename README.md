# FracMinHash + Bloom filters for ANI / AF on short contigs

Two deliverables:

1. **Theory**: [`fmh_bloom_filter_theory.pdf`](fmh_bloom_filter_theory.pdf) (27 pages). Variance framework, 14 designs, 16 scenarios.
2. **Implementation + benchmark**: [`aif_benchmark_report.pdf`](aif_benchmark_report.pdf) (15 pages). The `aif` tool (v1.1) compared with skani 0.4.0 at equal sketch budget, including the v1.1 speed rewrite.

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

| workload | skani (sketch + query) | aif v1.0 | **aif v1.1** |
|---|---|---|---|
| 20 Mbp of contigs vs E. coli, sketch + query | 0.095 s | 0.373 s | **0.084 s** |
| 81 Mbp of contigs vs E. coli, sketch + query | 0.300 s | 1.319 s | **0.259 s** |
| 81 Mbp, 8 threads (aif: query only) | 0.452 s | n/a | **0.122 s** |
| 4 references x 700 contigs, total | **0.136 s** | 0.545 s | 0.178 s |

v1.1 queries are 5.2x faster than v1.0 and sketching is 1.9x faster, at equal accuracy. skani is still faster when a
few hundred contigs meet a newly sketched reference (aif's sketch build costs more); aif is faster on larger
query sets and scales with threads. The fast paths need AVX-512; a scalar fallback gives byte-identical sketches.

## The `aif` tool

`src/aif.cpp` is a single dependency-free C++17 file.

* **Sketch.** Canonical k-mers are FracMinHash-selected into two hash-range tiers, each stored as a blocked Bloom filter (FPR 0.005 and 0.1; 0.88 bits/base by default). The exact FPR comes from the fill ratio. The spurious-match rate is calibrated empirically, and k is chosen per reference (17-21) to keep it at or below 0.003.
* **Query** (raw contigs). Every selected k-mer is probed. A multi-level HMM (one unaligned state plus aligned identity levels) gives the aligned fraction and the FP-corrected aligned containment, so ANI = C^(1/k) with bias correction. Heterogeneous contigs use posterior-mean ANI instead.
* **v1.1 speed**: 2-bit packing and AVX-512 k-mer selection, gather-based Bloom tests, a binned single-precision HMM with exact per-probe refinement around alignment boundaries, mmap input, vectorised sketch insertion, and `-t N` threads with ordered output.

```
make aif
bin/aif sketch ref.fa -o ref.aif          # -k auto --budget 0.88 --tiers 0.005:0.6,0.1:0.4
bin/aif query ref.aif contigs.fa -t 8 > out.tsv  # contig ani af theta probes hits aligned_probes status lr het low_frac
```

## Layout

```
fmh_bloom_filter_theory.pdf, aif_benchmark_report.pdf
src/aif.cpp                    the tool
bench/                         download_genomes.sh, make_scenarios.py, run_skani.sh, run_aif.sh, run_all.sh,
                               evaluate.py, tune.py, time_all.py (interleaved timing), time_compare.sh, medtime.sh
results/                       committed summaries: evaluations (JSON), predictions, timings, tuning tables
scripts/                       theory tables/figures, Monte Carlo check, benchmark-report tables/figures
report/, report_bench/         groff sources of the two PDFs (custom macros; built with groff + ghostscript)
data/  (git-ignored)           genomes, generated contigs, raw run outputs; regenerate with bench/ scripts
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

The scripts expect skani at the path in `bench/run_skani.sh` (override with `SKANI=...`), plus numpy and matplotlib for the Python parts.
