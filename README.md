# FracMinHash + Bloom filters for ANI / AF on short contigs

**Start here:** [`aif_method_and_results.pdf`](aif_method_and_results.pdf) (15 pages): the method and every result so far,
with pointers to the detailed reports below.

Five detailed reports:

1. **Theory**: [`fmh_bloom_filter_theory.pdf`](fmh_bloom_filter_theory.pdf) (27 pages). Variance framework, 14 designs, 16 scenarios.
2. **Implementation + simulated benchmark**: [`aif_benchmark_report.pdf`](aif_benchmark_report.pdf) (15 pages). The `aif` tool (v1.1) compared with skani 0.4.0 at equal sketch budget, including the v1.1 speed rewrite.
3. **Real genomes**: [`aif_real_genomes_report.pdf`](aif_real_genomes_report.pdf). 240 real assemblies and MAGs from the skani paper's data sets, with MUMmer ground truth; the v1.2/v1.3 speed work; v1.4 calibrated confidence intervals; an analysis of the AF model; the database mode (v1.5) against 2,005 genomes and (v1.6) all 199,923 GTDB species representatives; a Pareto sweep of sketch sizes; identity maps; contamination screening; comparisons with Kraken2, CheckM2 and GUNC; a v1.7 addendum (AF near the local identity, refitted intervals).
4. **Viral clustering**: [`aif_viral_clustering_report.pdf`](aif_viral_clustering_report.pdf) (8 pages). vOTU clustering (95% ANI over 85% AF) of 2 x 10,000 phage genomes/contigs (INPHARED, Gut Phage Database) against BLAST/CheckV truth, vs skani and vclust.
5. **Human gut catalogue**: [`aif_uhgg_contamination_report.pdf`](aif_uhgg_contamination_report.pdf). All 286,621 UHGG v3.0 genomes (276,463 MAGs, 10,158 isolates) screened contig by contig against GTDB R232: cross-genus shared sequence (mostly transferred DNA) and core-anchored contamination, checked against isolates, CheckM2, Kraken2 and recurrence in other strains. Also: a seven-tool benchmark on the screen's contigs against MUMmer (minimap2, LexicMap, sourmash, skani, Kraken2, CAT), with each tool's cost projected to the whole catalogue, and the v1.7 AF fix and interval tables that it led to.

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

## GTDB scale (v1.6): contig classification against all 199,923 GTDB R232 species representatives

Database: every GTDB R232 species representative (189,801 bacterial + 10,122 archaeal genomes, 655 Gbp; read
from the release's gzip files). Queries: 23,561 contigs >= 500 bp of 294 held-out GTDB non-representative draft
genomes from 294 species (15 archaeal); each query's own species representative is in the database. Truth:
GTDB taxonomy, and MUMmer for top-hit ANI/AF (72,135 contig-reference pairs).

| | skani 0.4.0 | **aif** (budget 0.88) | aif (budget 0.22) |
|---|---|---|---|
| database size (GB; bits / reference base) | 82.7 (1.011) | 82.5 (1.008) | 28.5 (0.348) |
| build, 32 threads, genomes in page cache (min; peak GB) | 10.8; 15.6 | **6.7; 1.9** | 14.7; 1.7 |
| search, 23,561 contigs, 16 threads, warm (s) | 252 | **11.5** | 7.4 |
| contigs with a hit | 71.6% | **93.5%** | 93.1% |
| top hit = query's species / genus | 61.8% / 69.0% | **78.3% / 88.5%** | 77.8% / 88.1% |
| species, contigs homologous to their own representative | 72.3% | **91.4%** | 90.9% |
| species calls at top-hit ANI >= 95%: sensitivity / precision | 60.9% / 87.5% | **76.6% / 86.4%** | 75.8% / 86.2% |
| species accuracy, 0.5-1 / 1-2 / 5-20 / >=20 kb contigs | 2.7 / 18.5 / 76.1 / 90.5% | **45.7 / 62.6 / 83.1 / 92.5%** | 44.9 / 61.8 / 82.9 / 92.3% |
| top-hit ANI MAE / AF MAE (pp) | 0.60 / 3.34 | **0.49 / 2.58** | 0.58 / 3.84 |

Contamination screening (simulated bins: each query genome plus 5/10/20% of another's contigs, 861 bins; contig
flagged when its top-hit genus differs from the bin's majority): aif flags 90.8% of contaminant contigs (75.5% of
those of 0.5-2 kb) at 75% precision, and 98.8% of contaminant bases at 95% precision; skani flags 67.8% of
contigs (13.3% of 0.5-2 kb) at 82% precision, and 97.6% of bases at 93% precision. Calibrated intervals cover 91.5% (ANI) / 91.9% (AF)
on GTDB's diverse taxa (95% nominal; calibrated on five bacterial species).

**Held-out confirmation (frozen defaults, git tag `v1.6-freeze`).** Database-mode defaults chosen while looking at
the benchmarks (screening setup, pre-HMM filter, key cap, likelihood-ratio hit order) were frozen and re-tested
on fresh queries: GTDB query set 2 (298 genomes from 300 unseen species, 19,555 contigs): aif species / genus
accuracy 78.6% / 90.1% (skani 66.1% / 73.3%), ANI>=95% calls 76.6% at 85.8% precision (skani 65.0% at 88.8%),
search 12.8 s vs 212 s; 80 fresh drafts against the 2,005-genome database: top-hit ANI / AF MAE 0.23 / 0.77 pp
(skani search 0.41 / 2.26, dist 0.52 / 2.93), coverage 99.9% (80.2%, 92.4%). Report Section 15.

**Kraken2 on the same references (report Section 16).** Kraken2 2.1.3 built from the same 199,923 representatives
(GTDB taxonomy as an NCBI-style taxdump), capped at aif's database size (`--max-db-size` 82.5 GB, which keeps
about 8% of its minimizers; uncapped it estimated 990 GB), `--fast-build`. Query set 1 / set 2:

| | skani 0.4.0 | Kraken2 2.1.3 | **aif** (budget 0.88) |
|---|---|---|---|
| build (min; threads; peak GB) | 10.8; 32; 15.6 | 233; 48; 88.7 | **6.7; 32; 1.9** |
| search, 16 threads, warm (s; peak GB) | 252 / 212; 78 | 73.7 / 82.1; 83 | **11.5 / 12.8**; 65 |
| contigs assigned | 71.6% / 75.5% | **97.5% / 98.7%** | 93.5% / 94.6% |
| species correct | 61.8% / 66.1% | 76.8% / 78.1% | **78.3% / 78.6%** |
| genus correct | 69.0% / 73.3% | **90.1% / 92.1%** | 88.5% / 90.1% |
| species calls: sensitivity / precision (aif, skani: ANI >= 95%; Kraken2: species rank) | 60.9% / 87.5% (set 1) | 76.8% / 84.0% (set 1) | 76.6% / 86.4% (set 1) |
| species correct, 0.5-1 kb contigs | 2.7% / 3.0% | 44.2% / 39.9% | **45.7% / 40.7%** |
| ANI and AF per hit | yes | no | yes, with intervals |

**Bin quality: CheckM2 and GUNC (report Section 16).** 400 simulated bins (100 each at 5/10/20% contamination
and 100 uncontaminated controls), 24 threads; contamination = contaminant bases / host bases; aif and skani
estimate it as flagged / unflagged bases (genus rule). GUNC 1.1.1 ran with DIAMOND 2.2.8 via
`GUNC_SKIP_DIAMOND_VERSION_CHECK` (it pins 2.1.24).

| | controls flagged | MAE at 5 / 10 / 20% (pp) | flagged at 5 / 10 / 20% | time, 400 bins |
|---|---|---|---|---|
| **aif** (estimate >= 5%) | 2% | **0.55 / 0.63 / 0.84** | 86 / 98 / 99% | **11 s** |
| aif (estimate >= 3%) | 6% | | 95 / 100 / 99% | |
| skani (>= 5%) | 4% | 0.89 / 1.03 / 1.35 | 77 / 96 / 99% | 206 s |
| CheckM2 1.1.0 (>= 5%) | 0% | 2.87 / 5.23 / 10.70 | 23 / 61 / 70% | 23.5 min |
| GUNC 1.1.1, own call | 5% | 0.69 / 1.33 / 3.86 | 94 / 100 / 99% | 6.6 h (upper bound: shared host) |

Every host species has a representative in the database, which favours aif and skani; the gene-based tools need
no close reference (and CheckM2 also estimates completeness). MMseqs2 taxonomy and CAT (about 700 million GTDB
proteins) were not run.

## Database mode (v1.5-1.6): contigs against 2,005 genomes

`aif db build` puts many references in one memory-mapped file with one k (20) and one tier layout. A contig is then
hashed once, and each candidate reference costs only its filter tests and the HMM. Candidates come from a
bit-sliced (COBS-style) screening index at 0.155 bits/base; v1.6 adds a pre-HMM local-score filter (output
unchanged), caps screening keys per contig at 128, and orders hits by the model's likelihood ratio.
Exhaustive search (`--no-screen`) is byte-identical to pairwise `aif query` against `-k 20` sketches.

Benchmark database: the 5 pairwise references, 100 RefSeq complete genomes each of E. coli, S. aureus,
P. aeruginosa and the B. cereus group, and 1,600 other RefSeq reference genomes (8.62 Gbp). All query
assemblies are excluded. Queries are the 33,464 within-species contigs; truth is MUMmer for every
(contig, reference) pair among the tools' top 3 hits. Pooled over the 26,084 homologous contigs (aif v1.6):

| | skani search | skani dist (2,005 FASTA) | **aif db search** | aif, size-matched DB |
|---|---|---|---|---|
| database (bits / reference base) | 0.896 | (FASTA) | 1.035 | 0.895 |
| homologous contigs with a hit | 77.0% | 94.5% | **99.9%** | 99.9% |
| top-hit ANI MAE (pp) | 0.513 | 0.690 | **0.310** | 0.320 |
| paired top-hit ANI MAE vs aif (pp) | 0.513 vs 0.246 | 0.690 vs 0.286 | | |
| top-hit AF MAE (pp) | 3.03 | 4.15 | **1.12** | 1.22 |
| top hit is the best reference | 78.3% | 74.9% | **88.9%** | 88.8% |
| regret, true ANI x AF / true ANI (pp) | 1.81 / 0.15 | 2.24 / 0.25 | **0.28 / 0.06** | 0.27 / 0.06 |
| search time, 33,464 contigs, 1 / 8 threads (s) | 264.5 / 37.9 (`--keep-refs` 83.8 / 15.9) | - / 118.5 | **32.9 / 4.4** | |
| peak memory, 8 threads (GB) | 1.32 (1.87) | 1.62 | 1.77 | |

Calibrated intervals keep 94.0% (ANI) / 94.5% (AF) coverage at 95% in database mode; screening recall against
exhaustive search is 97.5%.

## Pareto sweep (v1.6): sketch size vs accuracy vs time, real sets, pooled within-species

| method | bits / base | ANI MAE (pp) | AF MAE (pp) | sketch + query, 1 thread (s) |
|---|---|---|---|---|
| **aif**, budget 0.11 / 0.22 / 0.44 / 0.88 / 1.76 | 0.11-1.76 | 0.67 / 0.62 / 0.56 / 0.50 / 0.46 | 6.3 / 4.9 / 3.6 / 2.6 / 1.7 | 1.4 / 1.6 / 2.0 / 2.9 / 5.0 |
| skani, c = 500 / 250 / 125 / 60 / 30 | 0.26-3.27 | 0.84 / 0.74 / 0.68 / 0.63 / 0.61 | 26.9 / 15.4 / 9.4 / 6.0 / 4.3 | 2.4 / 2.7 / 3.2 / 4.4 / 7.1 |
| exact FracMinHash (sourmash), scaled 580 ... 18 | 0.11-3.56 | 0.97 ... 0.71 | (no AF) | 21-184 |
| minimap2 asm20 (alignment) | 27.4 (index) | 0.23 | 1.78 | 622 |

Exact FracMinHash at 32x aif's smallest sketch is still less accurate than it (the bits-per-sample theorem).

## Identity maps (v1.6): `--segments`

`aif query ... --segments FILE` (and `aif db search`) writes aligned/unaligned segments tiling each contig and
windowed identity (`--seg-window 1000`) with calibrated 95% intervals. Against MUMmer (224 genomes; window truth
from alignment differences): 3 kb window identity MAE 0.55 pp (fastANI fragments 0.78 pp; r = 0.97 vs 0.86;
0.035 vs 0.235 pp at >= 99.5% identity, fastANI better below 90%); 1 kb windows 0.63 pp; calibrated window
intervals cover 94.7% (1 kb) / 94.4% (3 kb) of held-out windows; unaligned islands found at 97% precision, 81%
recall (fastANI 56% / 51%); breakpoints within a median 171 bp (92% within 500 bp).

## Viral clustering (vOTUs; 10,000 sequences per data set; truth: BLAST megablast + CheckV anicalc/aniclust)

| INPHARED / GPD | vOTU-edge F1 | ARI vs BLAST vOTUs | near-threshold ANI / AF MAE (pp) | all-vs-all time, 16 thr (s) |
|---|---|---|---|---|
| aif (default AF) | 0.926 / 0.861 | 0.891 / 0.775 | 0.57 / 6.3 ; 0.63 / 7.2 | 8.0 / 4.9 |
| **aif `--af-min-ident 0.88`** | **0.938 / 0.902** | **0.896 / 0.807** | 0.57 / 2.5 ; 0.63 / 3.2 | 7.7 / 4.8 |
| skani | 0.915 / 0.833 | 0.888 / 0.723 | 0.99 / 6.9 ; 1.10 / 6.5 | 3.1 / 2.4 |
| skani `--small-genomes` | 0.917 / 0.857 | 0.890 / 0.750 | 0.89 / 2.4 ; 0.94 / 3.0 | 4.0 / 1.7 |
| vclust (alignment) | 0.970 / 0.948 | 0.937 / 0.893 | 0.34 / 1.0 ; 0.55 / 1.5 | 74 / 35 |

aif's default AF counts divergent homologous stretches that megablast does not align (the vOTU convention);
`--af-min-ident T` reports the AF at identity >= T (T = 0.88 chosen on a disjoint tuning set). Best among the
sketch tools; vclust's alignment is closest to BLAST at ~10x aif's time. `aif db build -i` makes every FASTA record a reference.

## Human gut catalogue (UHGG v3.0, 286,621 genomes, every contig >= 500 bp against GTDB R232)

Search (aif v1.7): 55.7 M contigs (707 GB) in 2.2 h at 48 threads, 89 GB peak memory. A contig is *foreign* when its best hit
is >= 95% ANI over AF >= 0.5 to another genus and nothing of the host (genus, or a reference conspecific with the
genome) is among its top 5 hits.

| | isolates (10,158) | MAGs (276,463) | high-quality MAGs (125,972) |
|---|---|---|---|
| any foreign (shared) contig | 74.4% | 63.9% | 60.3% |
| core-anchored contamination, any | 1.5% | 9.6% | 6.7% |
| contamination >= 1% | 0.47% | 2.68% | 1.73% |
| contamination >= 5% | 0.02% | 0.20% | 0.08% |

- **Shared sequence is mostly transferred DNA**: non-mobile foreign contigs of MAGs carry 0.29-0.70x the core-gene
  density of unflagged contigs of the same length and 71-74% hypothetical proteins; isolates (never binned) carry
  shared sequence as often as MAGs. 1.73 M contigs link 3,117 genera (46% of pairs across families).
- **Contamination is therefore core-anchored**: called only when a source species' foreign contigs carry a core gene
  (ribosomal proteins, RNA polymerase, EF-Tu/G, gyrase, RecA, tRNA ligases, ...). Within 256 species, MAGs are more
  contaminated than isolates in 207 and less in 21 (p = 1e-39); the rate rises from 3.9% to 32.6% across MGnify's
  CheckM2 contamination bins (0-1% ... 4-5%). A core-gene mixture puts 45% (44-46%) of non-mobile shared MAG
  sequence as misbinned fragments, so per-genome calls are conservative.
- **Kraken2 agrees on the genera**: on the same GTDB R232 references, it puts 95.0% of anchored MAG contigs in the
  source genus aif found (2.0% in the host genus) and 83.6% of length-matched unflagged control contigs in the
  host genus. This checks the search, not the contamination/transfer split, which rests on core genes.
- **GTDB R232 label problem**: five *E. coli*/*Shigella* assemblies are filed under other genera (ECMA0423,
  G047199095, G047199115, "Enterobacter flexneri"); with the pre-specified rule, isolates looked more contaminated
  than MAGs (20.7%). References conspecific with the genome by GTDB's criterion now count as the host (the one
  change after the held-out screen; `results_uhgg/prespecified/` keeps the first version).

### Seven tools on the screen's contigs (aif v1.7; report Section 7)

Stratified sample of the screen's contigs (S3: 10,005 contigs, 114 Mbp), scored against MUMmer alignments of every tool's top hits (`delta-filter -q`). Each tool was timed at 48 pinned threads, and the time projected to the 706 Gbp catalogue from a pool of 3.64 Gbp (`bench_real/uhgg_compare.py`, `results_uhgg/compare/`).

| | top-hit ANI MAE (pp) | AF MAE (pp) | top hit = best ref. | foreign-contig F1 | full catalogue | search memory |
|---|---|---|---|---|---|---|
| **aif v1.7** | 0.75 | 3.58 | **89.5%** | 0.871 | **1.8 h** (2.2 h measured) | 82 GB |
| minimap2 2.31 (asm20) | **0.32** | **2.73** | 87.9% | **0.881** | 10 days | 177 GB |
| sourmash (manysearch) | 1.23 | - | 89.0% | 0.844 | 19 h | 5 GB |
| skani 0.4.0 | 1.08 | 5.40 | 46.3% | 0.578 | 2.5 h | 65 GB |
| Kraken2 2.1.3 | - | - | - | 0.609 | 1.4 h | 81 GB |
| CAT_pack 7 (DIAMOND) | - | - | - | 0.411 | ~100 days | 53 GB |
| LexicMap 0.10.1 (S1 only) | 1.22 | 1.41 | 86.1% | 0.848 | decades (disk-bound) | 3 GB |

minimap2 is the most accurate tool. aif makes the screen's decisions nearly as well, at the cost of a sketch tool. The benchmark exposed aif v1.6's AF over-call (5.13 pp, bias +5.0), which v1.7 fixes (see "The `aif` tool" below).

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
* **v1.5**: database mode (`aif db build/search/info`): shared-k two-stage probing, bit-sliced screening index, memory-mapped database, lazily built per-reference estimators, top-n hits per contig.
* **v1.6**: `db build -i` (each FASTA record a reference); `--af-min-ident` (identity-thresholded AF, e.g. for BLAST-convention vOTU AF); gzip input; `--segments` identity maps (aligned/unaligned segments, windowed identity with calibrated intervals); database search scaled to GTDB (pre-HMM local-score filter, capped screening keys, vectorised candidate scan, parallel build sizing), hits ordered by likelihood ratio.
* **v1.7**: AF counts the aligned posterior mass within `--af-local` (0.05) of the local identity. Unaligned stretches inside an alignment were absorbed by the lowest identity levels, whose expected hit rates are close to the filters' false-positive rates, so AF read about 5 pp high on real contigs. Windows with their own evidence of alignment keep divergent blocks. Hits and ANI are unchanged, and `--af-local 0` reproduces v1.6 exactly. Numbers in the reports above are v1.6, except the gut-catalogue report (screen rerun with v1.7; contamination rates moved by at most 0.23 pp) and the real-genome report's v1.7 addendum. The setting was chosen on tuning sets (`bench_real/uhgg_af.py`, `bench_real/af_local_check.py --tune`). Held-out AF MAE: UHGG contigs against GTDB R232 5.13 → 3.58 pp; real-genome sets 2.45 → 1.97 and fresh drafts 2.53 → 1.97. Simulated benchmarks get worse, 5.30 → 5.83 and 4.95 → 5.51, mostly on heterogeneous-identity contigs (`hetblock` −5.8 pp). Intervals were refitted for the new AF. `db search` has its own conformal table (`src/ci_table_db.inc`, fitted on the UHGG tuning set); `query` keeps the real-tuning-set table (`src/ci_table.inc`). Held-out 95% coverage, ANI / AF: UHGG S3 top hits 0.917 / 0.868 → 0.951 / 0.942; real genomes 0.946 / 0.943, fresh drafts 0.948 / 0.957. The median AF half-width on real genomes fell from 2.60 to 1.98 pp. Genome-level estimates aggregated from contigs: ANI error is unchanged (0.140 → 0.141 pp), and AF error against ANIm halves (1.87 → 0.97 pp; `bench_real/af_local_genome.py`).
* **v1.4**: per-contig confidence intervals (statistical, or conformal-calibrated from the table in `src/ci_table.inc`, refit with `bench_real/calibration.py fit`); optional local identity prior (`--ident-window`) and overdispersed emissions (`--overdisp`).

```
make aif
bin/aif sketch ref.fa -o ref.aif          # -k auto --budget 0.88 --tiers 0.005:0.6,0.1:0.4
bin/aif query ref.aif contigs.fa -t 8 > out.tsv  # contig ani af theta probes hits aligned_probes status lr het low_frac
                                                #   ani_lo ani_hi af_lo af_hi ani_se af_se   (95% intervals; --ci-level, --ci-stat)
bin/aif db build --list refs.txt -o refs.aifd -t 8   # refs.txt: "path [name]" per line; -k 20 --budget 0.88
bin/aif db search refs.aifd contigs.fa -t 8 -n 5 > hits.tsv  # contig reference rank <query columns> screen_count ...
bin/aif query ref.aif contigs.fa --segments seg.tsv   # + segments / windowed identity map (contig start end aligned ani ani_lo ani_hi probes hits)
bin/aif db info refs.aifd --refs
```

## Layout

```
fmh_bloom_filter_theory.pdf, aif_benchmark_report.pdf, aif_real_genomes_report.pdf
src/aif.cpp, src/ci_table.inc  the tool (+ built-in confidence-interval calibration table)
bench/                         simulated benchmark: download_genomes.sh, make_scenarios.py, run_skani.sh, run_aif.sh,
                               run_all.sh, evaluate.py, tune.py, time_all.py (interleaved timing), time_compare.sh, medtime.sh
bench_real/                    real genomes: select_datasets.py, select_tune.py, download.sh, prepare.py, truth.py (MUMmer),
                               run_tools.py, evaluate_real.py, genome_level.py, blast_audit.py, tune_real.py,
                               time_real.py, scaling.py, consistency.py, calibration.py, af_model_grid.py, afopt_heldout.py,
                               database mode: db_select.py, db_download.sh, db_bench.py, db_screen_eval.py;
                               GTDB: gtdb_download.sh, gtdb_select.py, gtdb_prepare.py, gtdb_bench.py, gtdb_contam.py;
                               other tools: kraken2_gtdb.py (Kraken2), contam_tools.py (CheckM2, GUNC);
                               pareto_sourmash.py, pareto_minimap2.py, pareto_time.py; idmap.py;
                               viral clustering: viral_download.sh, viral_bench.py;
                               gut catalogue: uhgg_download.sh, uhgg_screen.py, uhgg_run.sh, uhgg_summary.py
results/, results_real/,       committed summaries: evaluations (JSON), predictions, truth, timings, tuning tables
results_viral/, results_uhgg/  (the full UHGG flagged-contig lists, 54-57 MB, are not committed; `uhgg_summary.py stats` rewrites them)
scripts/                       theory tables/figures, Monte Carlo check, report tables/figures (benchmark_report.py, real_report.py)
report/, report_bench/,        groff sources of the three PDFs (custom macros; built with groff + ghostscript)
report_real/, report_viral/,
report_uhgg/, report_overview/
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
python3 -I bench_real/db_select.py && bench_real/db_download.sh      # database mode
for s in build search truth time evaluate; do python3 -I bench_real/db_bench.py $s; done
bench_real/gtdb_download.sh && python3 -I bench_real/gtdb_select.py   # then NCBI datasets download of the queries
python3 -I bench_real/gtdb_prepare.py && for s in list build search truth evaluate; do python3 -I bench_real/gtdb_bench.py $s; done
python3 -I bench_real/gtdb_contam.py
python3 -I bench_real/select_fresh.py; REAL_DIR=data/real_fresh bench_real/download.sh   # fresh validation queries
REAL_DIR=data/real_fresh python3 -I bench_real/prepare.py; REAL_DIR=data/real_fresh REAL_SETS=bcereus,ecoli,paer,saur python3 -I bench_real/truth.py
DB_QUERY_DIR=data/real_fresh DB_SETS=bcereus,ecoli,paer,saur DB_RES=data/db/results_fresh python3 -I bench_real/db_bench.py search  # + truth, evaluate
GTDB_SET=2 python3 -I bench_real/gtdb_select.py   # + download, gtdb_prepare.py, gtdb_bench.py search/truth/evaluate, gtdb_contam.py
for s in prepare build classify evaluate; do python3 -I bench_real/kraken2_gtdb.py $s; done   # Kraken2 2.1.3
for s in bins checkm2 gunc seqtime evaluate; do THREADS=24 python3 -I bench_real/contam_tools.py $s; done  # CheckM2, GUNC
python3 -I bench_real/run_tools.py aif_b011 aif_b022 aif_b176 skani_c60 skani_c250 skani_c500
<python with sourmash> -I bench_real/pareto_sourmash.py; python3 -I bench_real/pareto_minimap2.py; python3 -I bench_real/pareto_time.py
python3 -I bench_real/idmap.py run; REAL_DIR=data/real_tune python3 -I bench_real/idmap.py run
REAL_DIR=data/real_tune python3 -I bench_real/idmap.py calibrate; python3 -I bench_real/idmap.py evaluate
python3 scripts/real_report.py && report_real/build.sh
bench_real/viral_download.sh && for s in prepare blast prepare_tune; do python3 -I bench_real/viral_bench.py $s; done
python3 -I bench_real/viral_bench.py blast tune && for s in tune tools evaluate; do python3 -I bench_real/viral_bench.py $s; done
/usr/bin/python3 scripts/viral_report.py && report_viral/build.sh
bench_real/uhgg_run.sh        # UHGG v3.0: download, prepare, search, per-genome calls (UHGG_PILOT=3000 for the pilot)
/usr/bin/python3 -I bench_real/uhgg_summary.py stats && /usr/bin/python3 -I bench_real/uhgg_summary.py kraken
/usr/bin/python3 scripts/uhgg_report.py && report_uhgg/build.sh
/usr/bin/python3 scripts/overview_report.py && report_overview/build.sh   # method-and-results overview
```

The scripts expect skani at the path in `bench/run_skani.sh` (override with `SKANI=...`), plus numpy and matplotlib for the Python parts.
