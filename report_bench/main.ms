.\" Build: report_bench/build.sh
.EQ
delim $$
gsize 11
define eps 'epsilon'
define Chat '{ C hat }'
.EN
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
aif: FracMinHash + Bloom-Filter Sketches
for Contig ANI and Alignment Fraction
.ft R
.ps 14
.vs 18p
Implementation and Benchmark against skani
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Implementation report v1.1 \(em November 21, 2026
Implementation, experiments and typesetting by Claude (Anthropic) for the project owner
.ps \n(PS
.vs \n(VSp
.sp 1.2v
.in 0.5i
.ll 6.0i
.ce 1
\f[B]Abstract\f[R]
.sp 0.3v
.fi
.ad b
.ps 10
.vs 12.4p
We implemented the strategy that the theory report identified as best for short contigs queried
against a sketched reference. It is a hedged two-tier FracMinHash/Bloom-filter sketch probed densely by the
raw contig, and it estimates ANI and AF jointly with a multi-level hidden Markov model over the probe hits.
We compared it with skani 0.4.0 at the same sketch budget: 0.88 against 0.90 bits per reference base.
The benchmark has 2,800 contigs over 28 conditions covering the report's scenarios: length 1\(en20 kb,
ANI 80\(en99.9%, partial overlap, indels, rate heterogeneity and chimeras. Contigs were cut from four real
bacterial genomes with exactly known ANI and AF.
On the held-out benchmark, aif (version 1.1) reports more contigs (95.1% vs 86.8%) and reduces the ANI mean
absolute error from 1.12 to 0.36 percentage points (pp) over all contigs and from 1.45 to 0.50 pp averaged over
conditions. It reduces the AF error from 7.43 to 0.39 pp (all contigs) and from 9.98 to 0.51 pp (condition average).
aif's ANI and AF errors are lower than skani's in all 28 conditions, and a freshly generated
confirmation benchmark reproduces this (28/28; ANI 0.37 vs 1.21 pp, AF 0.55 vs 7.46 pp).
aif also beats skani's densest preset (\f[CR]--slow\f[R], 3.8\(mu the space) on both errors, even at half its
own budget.
A speed rewrite (v1.1: AVX-512 k-mer extraction and Bloom tests, a binned single-precision HMM with exact
boundary refinement) made queries 5.2\(mu and sketching 1.9\(mu faster than v1.0 at equal accuracy.
aif now uses less memory than skani (at most 15.5 vs 19.0 MB) and is faster on larger query sets: 0.26 vs
0.30 s for 81 Mbp of contigs single-threaded, and 0.12 vs 0.45 s with 8 threads. For a single reference
with only 700 contigs skani remains faster end to end (0.14 vs 0.18 s over the four references) because
aif's sketch construction costs more.
The observed errors match Theorem 1 of the theory report to within 0\(en17% (median 4%).
.ps \n(PS
.vs \n(VSp
.in 0
.ll 6.5i
.bp
.ce 1
\f[B]\s(13Contents\s0\f[R]
.sp 0.8v
.so gen/toc.ms
.nr BODY 1
.pn 1
.bp
.\" =====================================================================
.H1 "Summary of results"
.LP
The task was to beat skani on contig ANI and AF error at the same sketch budget, over scenarios that
cover the theory report, and to record time and memory.
Table 1 gives the headline. Every number comes from the held-out main benchmark, except the last column,
which is a second, freshly generated benchmark used only once, for confirmation.
.TC 1
Headline comparison at equal budget. MAE: mean absolute error in percentage points (pp) of ANI or AF.
\(lqAll contigs\(rq pools every reported contig; \(lqmean over conditions\(rq averages the 28 per-condition MAEs;
\(lqpaired\(rq restricts each condition to contigs that both tools report. Times are single-threaded medians.
.TCE
.so gen/tab_headline.tbl
.IP \(bu 0.2i
\f[B]Targets met.\f[R] skani's ANI MAE (1.117 pp all contigs; 1.448 pp condition mean) and AF MAE (7.43; 9.98 pp)
were the targets. aif reaches 0.357 / 0.502 pp for ANI (3.1\(mu / 2.9\(mu lower) and 0.39 / 0.51 pp for AF
(19\(mu / 20\(mu lower), while reporting 8.3 percentage points more contigs.
.IP \(bu 0.2i
\f[B]Uniform improvement.\f[R] aif's errors are lower for both ANI and AF in 28 of 28 conditions, on both
the main and the fresh benchmark (Section 5.1, Figures 1\(en2).
.IP \(bu 0.2i
\f[B]Budget efficiency.\f[R] aif at 0.44 bits per base (half of skani's budget) still beats every skani preset,
including \f[CR]--slow\f[R] at 3.33 bits per base, on both errors (Section 5.5, Figure 4).
.IP \(bu 0.2i
\f[B]Resources.\f[R] aif's peak memory is lower (at most 15.5 vs 19.0 MB). It is faster than skani on query sets of
tens of Mbp (single-threaded and multithreaded) and slower only when a few hundred contigs meet a newly
sketched reference (Sections 5.6 and 6).
.IP \(bu 0.2i
\f[B]Theory confirmed.\f[R] For substitution-only contigs, Theorem 1 of the theory report predicts aif's
RMSE to within 0\(en17%, median 4% (Table 6).
.\" =====================================================================
.H1 "Benchmark design"
.H2 "Genomes and ground truth"
.LP
Four complete bacterial genomes spanning GC content serve as references and as templates for contigs:
\f[I]E. coli\f[R] K-12 MG1655 (NC_000913.3, 4.64 Mbp, GC 51%), \f[I]B. subtilis\f[R] 168 (NC_000964.3, 4.22 Mbp, 44%),
\f[I]P. aeruginosa\f[R] PAO1 (NC_002516.2, 6.26 Mbp, 67%) and \f[I]S. aureus\f[R] NCTC 8325 (NC_007795.1, 2.82 Mbp, 33%).
For each contig, a segment of the reference is sampled and reverse-complemented with probability 1/2.
It is then mutated with substitutions, plus indels or heterogeneous rates in some conditions, and optionally
padded with foreign sequence from another of the four genomes or built from two distant loci.
Truth is exact:
$roman ANI = 1 - $ substitutions / aligned columns, with indel columns excluded as in alignment-based ANI; and
$roman AF =$ homologous contig bases / contig length, with inserted bases inside homologous segments counted as homologous.
.H2 "Conditions"
.LP
Table 2 lists the 28 conditions. Each has 25 contigs per genome, 100 in total; each genome's 700 contigs are
queried against that genome. The scenario column maps conditions to the scenarios of the theory report.
.TC 2
Benchmark conditions (2,800 contigs). Lengths are before indels; the AF column is the homologous fraction.
.TCE
.so gen/tab_conditions.tbl
.H2 "Tools, metrics and protocol"
.LP
\f[B]skani\f[R] 0.4.0 runs with default parameters (c = 125, m = 1000, screen 80%, min-AF 15%),
\f[CR]-t 1 --qi\f[R], raw reference FASTA. On individual contigs skani does not apply its learned ANI
regression, so it reports its raw estimate.
\f[B]aif\f[R] sketches the reference (\f[CR]aif sketch\f[R]) and queries the raw contigs (\f[CR]aif query\f[R]),
single-threaded. Both tools receive identical inputs.
skani's sketch for these genomes occupies 0.8995 bits per reference base (seeds plus markers);
aif's sketch files occupy 0.8801 bits per base, header included.
.LP
Metrics are coverage (the fraction of contigs with an output), MAE of ANI and AF in pp, and mean signed error (bias).
Contigs a tool does not report are excluded from its MAE, so coverage must always be read alongside MAE.
\(lqPaired\(rq metrics remove that confound by using only contigs reported by both tools.
.LP
\f[B]Protocol.\f[R] Numeric parameters were tuned only on a separate \f[I]tuning set\f[R] (same 28 conditions,
10 contigs per genome, different random seed). The main benchmark was used during development to evaluate
successive design versions (Section 4), which is a mild form of test-set exposure. To guard against it, a
third \f[I]fresh benchmark\f[R] (25 contigs per genome, new seed) was generated after all decisions were frozen
and evaluated once (Table 1, last column).
Time and peak resident memory come from GNU time, as medians of 5 sequential runs on a 64-core Linux host
(503 GB RAM). Single-run timings on this host vary by up to 2\(mu.
.\" =====================================================================
.H1 "The implemented strategy"
.H2 "What was implemented, and why"
.LP
For a raw contig queried against a sketched reference, the theory report predicted the following gains:
(i) many cheap, noisy probes instead of a few exact ones (Theorem 3), (ii) hedged FPR tiers when the ANI is
unknown (D8), and (iii) AF from the query's own coordinates (D12, Proposition 5). It also predicted that,
for short contigs, reference-side anchors and island keys (D3/D5) cost more than the filter, and that their
benefits (colinearity, chimeras, multi-reference localisation) are secondary in this setting (scenario S1).
aif therefore implements D2 + D8 + D12, with the FP-corrected estimator of Theorem 1, the bias correction of
Corollary 1, the heterogeneity remedy of Proposition 9, and a per-reference choice of $k$ that controls
spurious matches (Proposition 6).
Anchors, island keys and island fingerprints (D6) were left for later (Section 7).
.H2 "Reference sketch"
.IP \(bu 0.2i
\f[B]Selection.\f[R] Every canonical $k$-mer is hashed with a 64-bit mixer. Tiers occupy disjoint ranges of
the hash, so tier $t$ holds $roman FMH sub {s sub t}$-style consistent samples, and no $k$-mer is stored twice.
.IP \(bu 0.2i
\f[B]Filters.\f[R] Each tier is a cache-blocked Bloom filter (512-bit blocks) with an independent hash family.
Its share of the budget fixes its size, and its sampling rate is $s sub t = beta w sub t / c( eps sub t )$, using
the exact integer-$h$ Bloom cost. Final tiers: $eps$ = 0.005 with 60% of the bits, and $eps$ = 0.1 with 40%.
The exact FPR of the built filter, the mean over blocks of $( roman fill ) sup h$, is stored in the sketch (Table 3).
.IP \(bu 0.2i
\f[B]Spurious-match calibration.\f[R] About 1,000,000 reference $k$-mers at stratified random positions receive one random
substitution and are probed. The excess positive rate over the known FPR estimates $sigma$, the probability that a non-shared query
$k$-mer occurs elsewhere in the reference. The calibrated $sigma$ was 4\(en28\(mu the uniform-sequence estimate
(for example 0.020 vs 0.0007 for \f[I]P. aeruginosa\f[R] at $k = 17$) because of composition bias and near-repeats.
.IP \(bu 0.2i
\f[B]Adaptive $k$.\f[R] Choose the smallest $k >= 17$ with calibrated $sigma <= 0.003$. The starting $k$ is
predicted from base composition ($sigma approx 4 cdot 2n ( sum f sub b sup 2 ) sup k$) and then verified.
This gave $k = 17$ for \f[I]E. coli\f[R] and \f[I]B. subtilis\f[R] and 18 for \f[I]S. aureus\f[R], each with a single build,
and 20 for \f[I]P. aeruginosa\f[R], whose calibrated $sigma$ at the predicted $k = 19$ (0.0032) sat just above the target.
.TC 3
Built sketches: chosen $k$, sampling rate and exact FPR per tier.
.TCE
.so gen/tab_tiers.tbl
.H2 "Query"
.IP \(bu 0.2i
\f[B]Probing.\f[R] Every selected $k$-mer of the contig is probed in its tier's filter, and its filter block is
prefetched during hashing. This yields a sequence of (position, tier, hit) probes, about 12% of all $k$-mer positions.
.IP \(bu 0.2i
\f[B]HMM.\f[R] One \(lqunaligned\(rq state $U$ plus aligned states $A sub j$ on an ANI grid with 0.005 spacing,
from 99.95% down to an \f[I]identifiability floor\f[R]: only levels with $theta sub j = roman ANI sub j sup k >= 8 cdot min sub t eps sub t prime$
are kept, where $eps prime = eps + sigma - eps sigma$.
Emissions are Bernoulli: $Pr ( roman hit | U ) = eps sub t prime$ and $Pr ( roman hit | A sub j ) = theta sub j + (1 - theta sub j ) eps sub t prime$.
Switching rates are per base: $rho sub u$ = 3\(mu10\u\s-2\-6\s+2\d between aligned and unaligned, and $rho sub h$ = 10\u\s-2\-4\s+2\d between levels.
A scaled forward\(enbackward pass gives the aligned posterior of every probe and the posterior mean ANI.
In v1.1 the pass runs over bins of 32 consecutive probes (emissions from precomputed power tables), and an
exact per-probe pass is then run only around alignment boundaries, i.e. bins with posterior between 0.05 and
0.95 or at a 0.5 crossing, conditioned on the binned messages at the region edges (Section 6).
.IP \(bu 0.2i
\f[B]Reporting.\f[R] A contig is reported if the HMM log-likelihood exceeds that of the all-$U$ model by at least 10 nats.
.IP \(bu 0.2i
\f[B]AF\f[R] = (posterior-expected number of aligned $k$-mer starts + $(k - 1)$ per aligned segment) / contig length.
.IP \(bu 0.2i
\f[B]ANI.\f[R] The pooled estimate is $Chat sup {1/k}$, where $Chat$ is the posterior-weighted, FP- and $sigma$-corrected
containment (5) of the theory report, combined over tiers by inverse variance and bias-corrected with (8).
If a $chi sup 2$-type test over 1 kb windows of confidently aligned probes rejects homogeneity ($z > 3$),
the posterior-mean ANI is reported instead, which removes the power-mean bias of Proposition 9; in v1.1 it
comes from a finer HMM pass (bins of 8 probes) run only for such contigs.
.LP
All operations are linear in sequence length. The query costs $O(n sub probe cdot m)$ for $m approx$ 20\(en40
HMM levels, and memory is the sketch plus one contig.
The implementation is a single dependency-free C++17 file (\f[CR]src/aif.cpp\f[R]) with AVX-512 fast paths and a
portable scalar fallback that produces byte-identical sketches.
.\" =====================================================================
.H1 "Development and parameter tuning"
.LP
Table 4 traces the design through six versions, each scored on the main benchmark when it was built.
.TC 4
Design versions, scored on the main benchmark (2,800 contigs) at the time they were built.
.TCE
.so gen/tab_history.tbl
.IP v1 0.5i
Two tiers at $k = 21$ with a 2-state HMM (one aligned containment fitted by EM) already halved skani's ANI error
for contigs at 90% ANI and above. It failed in three ways. Under rate heterogeneity, one aligned level called
the divergent blocks \(lqunaligned\(rq (AF error 48 pp). At 80\(en85% ANI almost nothing was reported. And EM
made queries slow.
.IP v2 0.5i
Fifteen fixed aligned levels and a single pass fixed heterogeneity and speed. Calibrating $sigma$ exposed the
composition effect. With one tier at FPR 0.02, foreign sequence in partial-overlap contigs was absorbed into
the lowest aligned levels (AF +10\(en15 pp), because those levels are statistically almost indistinguishable from $U$.
.IP v3 0.5i
Restoring two tiers, so that low-FPR hits carry strong evidence, keeping only identifiable levels and adding
posterior-mean ANI under heterogeneity cut the AF error to 1.9 pp.
.IP tuned 0.5i
Two grid searches on the tuning set (378 + 1,134 configurations over $k$, tier FPRs and split, $rho sub u$,
$rho sub h$, the identifiability factor and the reporting threshold) found that $rho sub h$ must be well below
the earlier value and that the identifiability factor should be about 8. This step brought the largest gain.
.IP final 0.5i
Adaptive $k$ raised coverage from 94.2% to 95.5% (v1.0), mostly 85% ANI contigs on the high- and low-GC genomes,
and also lowered both errors. Flagging estimates that sit on the identifiability floor was tested and rejected,
because it cost 6 points of coverage for a negligible gain in ANI error.
.LP
Replacing EM with a single forward\(enbackward pass, and then low-level optimisation (a single build pass,
precomputed emission and transition tables, vectorisable loops, prefetching), reduced query time about 5\(mu.
.IP v1.1 0.5i
The speed rewrite of Section 6. Its accuracy-relevant decisions (bin size 32 with a bin-8 re-run for
heterogeneous contigs, refinement threshold 0.05, grid step 0.005, at most 6 Bloom hashes) were taken on the
tuning set; on the main benchmark it is as accurate as v1.0 or slightly better (Table 4).
.\" =====================================================================
.H1 "Results"
.H2 "Per-condition accuracy"
.LP
Table 5 and Figures 1\(en2 give every condition. aif's ANI and AF errors are lower than skani's in all 28.
The largest relative gains are on short contigs, where skani reports 32% of 1 kb contigs at 2.0 pp error
and aif reports all of them at 0.62 pp. They are also large for partial overlap, where skani's ANI error is 0.3\(en3.3 pp (mostly variance, bias up to
$-$0.85 pp) against 0.1\(en0.5 pp for aif, and for chimeras, where skani's AF is 13 pp off.
Under rate heterogeneity skani is biased high (+1.3 to +1.9 pp); aif's heterogeneity mode keeps the bias within $+-$0.3 pp.
At 80% ANI neither tool works: skani reports 7% of contigs with a 7.1 pp error, and aif reports 2% with a 4.3 pp error.
.TC 5
Per-condition results on the main benchmark (aif in bold).
.TCE
.so gen/tab_percond.tbl
.ne 7.2i
.PSPIC -C figs/fig_percond_ani.eps 5.6i
.FC 1
ANI mean absolute error per condition (log scale). Hollow markers: the tool reports fewer than half of
the condition's contigs.
.FE
.ne 7.2i
.PSPIC -C figs/fig_percond_af.eps 5.6i
.FC 2
AF mean absolute error per condition (log scale). For contigs fully contained in the reference,
aif's AF error is typically below 0.1 pp.
.FE
.H2 "Error distribution"
.LP
Figure 3 shows the signed error of every reported contig against its true ANI.
skani's errors spread by $+-$2\(en5 pp at 85\(en95% ANI. aif's are concentrated within about $+-$1 pp, with a
few outliers at 85% ANI and in partial-overlap contigs.
.ne 3.3i
.PSPIC -C figs/fig_scatter.eps 6.4i
.FC 3
Estimate minus truth (pp) for every reported contig. aif reports 2,662 contigs and skani 2,431.
.FE
.H2 "Agreement with the theory"
.LP
For substitution-only contigs at AF = 1, Table 6 compares aif's RMSE about the model ANI with the standard
deviation predicted by Theorem 1. The prediction uses each sketch's actual $k$, tier rates, exact FPRs and
calibrated $sigma$, averaged over the four references. Agreement is within 0\(en17% (median 4%) across lengths of
1\(en20 kb and ANIs of 90\(en99.9%. aif's error is 1.1\(en1.8\(mu that of the alignment oracle, against 2.8\(en4.2\(mu for skani.
.TC 6
Predicted (Theorem 1) vs observed error, in pp of ANI, measured about the condition's model ANI.
.TCE
.so gen/tab_theory.tbl
.H2 "Ablations"
.TC 7
Ablations on the main benchmark, all at 0.88 bits per base.
.TCE
.so gen/tab_ablation.tbl
.LP
Each component matters (Table 7). Hedged tiers versus a single tier cut the AF error from 1.13 to 0.39 pp
(all contigs) and from 1.31 to 0.51 pp (condition mean). The heterogeneity mode versus pooled-only ANI cuts the
all-contig ANI error from 0.45 to 0.36 pp.
Adaptive $k$ versus fixed $k = 17$ raises coverage from 94.1% to 95.1% and lowers the all-contig AF error from
0.50 to 0.39 pp; fixed $k = 17$ has a lower condition-mean ANI error only because it reports fewer hard
contigs. Fixed $k = 21$ is worse on every metric.
The first design (v1 tiers, under the final HMM) remains clearly worse.
With about 100 or more probes per contig, the delta-method bias correction is negligible.
.H2 "Accuracy versus budget"
.LP
Table 8 and Figure 4 sweep aif's budget and compare skani's presets, which densify its seeds.
aif at 0.44 bits per base, half of skani's default, beats skani \f[CR]--slow\f[R] (3.33 bits per base) on both errors.
At 0.22 bits per base it still beats default skani.
At 0.064 bits per base, the budget of plain FracMinHash at $s$ = 1/1000, aif's AF error approaches skani's and
its coverage drops to 74%, as the theory report predicts for very sparse sketches.
At 1.76 bits per base, ANI keeps improving (paired ANI MAE 0.30 vs 0.36 pp), but AF does not (paired 0.47 vs 0.39 pp).
The HMM treats overlapping $k$-mer probes as independent, so at double the density it overcounts evidence and
overreacts to mutation clusters. The parameters were tuned at 0.88, and tempering the emissions by probe density
should restore monotonicity (Section 7).
.TC 8
Accuracy versus sketch size.
.TCE
.so gen/tab_budget.tbl
.ne 3.2i
.PSPIC -C figs/fig_budget.eps 6.4i
.FC 4
ANI and AF MAE versus sketch size. The vertical guide marks skani's default budget.
.FE
.H2 "Time and memory"
.TC 9
Per-reference resources (v1.1; medians of 11 interleaved runs; 700 contigs, about 3.6 Mbp, per reference).
.TCE
.so gen/tab_resources.tbl
.TC 10
Throughput against \f[I]E. coli\f[R] on 2,800 contigs (20.3 Mbp) and 11,200 contigs (81 Mbp); skani times include
its sketching of the reference, aif's sketch (one per reference) is listed separately.
.TCE
.so gen/tab_throughput.tbl
.LP
aif's sketch files are 2% smaller than skani's, and its peak memory is lower in every single-threaded run
(Table 9). For one reference and 700 short contigs skani is still faster end to end (0.027\(en0.040 s against
0.037\(en0.062 s), because aif's sketch construction (0.020\(en0.044 s) costs more than skani's; aif's queries
alone take 0.017\(en0.019 s. On larger query sets aif is faster (Table 10): 0.084 vs 0.095 s for 20 Mbp and
0.26 vs 0.30 s for 81 Mbp single-threaded, and 0.12 vs 0.45 s with 8 threads, where skani gets slower and aif
scales. aif probes about 12% of all contig $k$-mers, roughly 15\(mu more than skani's 1/125 seeds; Section 6
explains how that work was made cheap. Timings are medians of interleaved, core-pinned runs after warm-up
(\f[CR]bench/time_all.py\f[R]); single runs on this host vary by up to 2\(mu with CPU frequency state.
.\" =====================================================================
.H1 "Speed optimisation (v1.1)"
.LP
v1.0 spent its query time in two places: scalar rolling-hash probing (about 40%) and a 30-state HMM run over
every probe (about 50%). v1.1 rewrites both, without changing the statistical model:
.IP (1) 0.3i
\f[B]Packed sequences and vectorised selection.\f[R] Contigs are packed to 2 bits per base with SWAR and BMI2
(8 bases per step; non-ACGT bases flagged). Canonical $k$-mers are extracted directly from the packed words,
8 positions at a time with AVX-512: variable shifts, a bit-reversal reverse complement, a vectorised fmix64
selection hash and mask-compression of the selected positions. There is no per-base dependency chain.
.IP (2) 0.3i
\f[B]Vectorised Bloom tests.\f[R] For 8 selected positions at a time the $k$-mer and both hashes are recomputed from
the packed words, and the filter bits are fetched with gathers across both tiers at once. An early exit on rejected
lanes was tried and removed: it serialises the gathers and was slower.
.IP (3) 0.3i
\f[B]Binned HMM with exact boundary refinement.\f[R] The forward\(enbackward pass runs in single precision over bins of
32 probes, with emissions from power tables, aligned states padded to a multiple of 16, and a normaliser computed
from one dot product (the aligned mass of a normalised step is $1 - f sub U$), so each step carries one reduction.
An exact per-probe pass then refines only the uncertain bins around alignment boundaries, using the stored binned
messages as boundary conditions. This recovers AF precision: the refined AF error is lower than either the coarse
binned pass or v1.0's per-probe HMM. Heterogeneous contigs get a bin-8 re-run for their posterior-mean ANI, and
contigs that every 32-probe window confidently places as aligned skip the HMM altogether.
.IP (4) 0.3i
\f[B]Input, sketching and threads.\f[R] FASTA is read through \f[CR]mmap\f[R], releasing consumed pages. The sketch build
packs each reference once, inserts with vectorised gather\(enOR\(enscatter (lanes that collide on the same word fall back to
scalar updates), and draws its calibration samples from 8 vectorised random streams on stratified positions; the
scalar fallback produces byte-identical sketches. The low-FPR tier is capped at 6 hash functions (its exact FPR
rises from about 0.0053 to 0.0058). Queries can run multithreaded with persistent workers and double-buffered batches,
and the output is identical for any thread count.
.LP
Table 11 compares v1.0, v1.1 and skani on the same workloads. Table 12 shows where v1.1 spends its time.
.TC 11
Speed of aif v1.0 and v1.1 against skani (medians of interleaved, core-pinned runs; skani includes sketching its
reference).
.TCE
.so gen/tab_speed.tbl
.TC 12
Where v1.1 spends its query time (81 Mbp of contigs against \f[I]E. coli\f[R], single thread; ms, from instrumented
runs, so the stages sum only approximately to the end-to-end time).
.TCE
.TS
center box tab(@);
l n l.
\f[B]stage\f[R]@\f[B]ms\f[R]@\f[B]v1.0 counterpart\f[R]
_
read FASTA (mmap)@26@buffered stdio, about 80 ms
pack to 2-bit words@23@\(em
select $k$-mers (AVX-512)@53@scalar rolling hash, about 510 ms in total for probing
Bloom tests (gathers)@59@
bins + fast-path check@10@\(em
binned HMM (B = 32)@18@per-probe HMM, about 650 ms
boundary refinement@23@\(em
heterogeneity re-run, post-processing@about 25@about 25 ms
.TE
.\" =====================================================================
.H1 "Discussion"
.H2 "Why aif wins"
.IP (1) 0.3i
\f[B]Information per stored bit.\f[R] At the same size, a Bloom tier at FPR 0.005\(en0.1 lets the raw query probe
about 15\(mu more $k$-mers than skani's seeds. Theorem 1 says that, below the knee $s approx 1/k$, error is governed by
the number of probes, and Table 6 confirms it quantitatively.
.IP (2) 0.3i
\f[B]AF from dense query-side coordinates.\f[R] The HMM locates alignment boundaries to within tens of bases,
because hits arrive every few bases. Chain ends from sparse seeds are hundreds of bases off. Separating AF from
ANI also removes the partial-overlap bias in ANI (Proposition 3).
.IP (3) 0.3i
\f[B]Calibrated noise.\f[R] Exact filter FPRs and the empirically calibrated $sigma$ make the containment estimator
unbiased, and adaptive $k$ keeps $sigma$ small on GC-skewed references.
.IP (4) 0.3i
\f[B]Heterogeneity.\f[R] Multi-level states and posterior-mean ANI remove most of the power-mean bias that every
containment-based estimator, skani's included, suffers.
.H2 "Limitations and threats to validity"
.IP \(bu 0.2i
\f[B]Asymmetric.\f[R] aif needs the query's raw sequence. skani also compares sketch against sketch.
For all-versus-all use, the theory report recommends $b$-bit FracMinHash (D9), which is not implemented here.
.IP \(bu 0.2i
\f[B]Speed.\f[R] For a few hundred contigs against a newly sketched reference, skani is still about 1.3\(mu faster
end to end because aif's sketch construction is slower; sketches are built once per reference in database use.
On larger query sets aif is faster, and its multithreaded query scales where skani's does not.
The fast paths need AVX-512; the scalar fallback is correct but several times slower.
.IP \(bu 0.2i
\f[B]Low ANI.\f[R] Below about 85% ANI, the identifiability floor (Section 3.3) limits detection on 5 kb contigs,
and both tools fail at 80%.
.IP \(bu 0.2i
\f[B]Indels.\f[R] aif counts an indel as about one substitution, so its ANI is biased low by about the indel rate:
$-$0.5 pp when indel events are 10% of substitutions and $-$1.0 pp at 20%. skani is less affected ($-$0.2 and $-$0.6 pp),
and aif's MAE is only slightly lower there (0.98 vs 1.09 pp at 20%). Reference anchors would allow the
span-discrepancy correction of the theory report (Section 5.6).
.IP \(bu 0.2i
\f[B]Colinearity.\f[R] Without reference positions, aif cannot flag chimeras or rearrangements; it reports the
union AF, which is the truth definition used here.
.IP \(bu 0.2i
\f[B]Simulated data.\f[R] Contigs are real genome segments with simulated evolution and a single reference per
query. Real genome pairs, with alignment-based ground truth, and multi-reference databases remain to be tested.
skani was designed for genome-scale comparisons, and its contig mode uses neither its learned ANI model nor
its database screening advantages.
.IP \(bu 0.2i
\f[B]Tuning exposure.\f[R] Design choices were partly guided by main-benchmark failures. The fresh benchmark
(Table 1) shows no sign of overfitting.
.\" =====================================================================
.H1 "Next steps"
.IP 1. 0.3i
Temper the HMM emissions by probe density, so that parameters stay optimal at any budget and accuracy improves
monotonically with space.
.IP 2. 0.3i
Add sparse reference anchors and island keys (D3/D5) for chimera and rearrangement detection, indel correction
and multi-reference localisation; then evaluate scenarios S6, S10\(enS13 of the theory report.
.IP 3. 0.3i
Implement island fingerprints (D6) for strain-level comparisons (ANI $>=$ 99.9%).
.IP 4. 0.3i
Test on real data: genome pairs across species and genera with alignment-based ANI/AF, and MAG contigs (experiments E2\(enE3).
.IP 5. 0.3i
Add a database mode with a screening tier, a faster (vectorised) sketch build, an AVX2 fast path, and a $b$-bit
FracMinHash variant for symmetric use.
.\" =====================================================================
.H1 "Reproducibility"
.CS
make aif                               # build bin/aif (C++17, no dependencies)
bench/download_genomes.sh              # 4 reference genomes from NCBI
python3 bench/make_scenarios.py data/genomes data/bench 25            # main benchmark
python3 bench/make_scenarios.py data/genomes data/tune 10 777         # tuning set
python3 bench/make_scenarios.py data/genomes data/bench2 25 99991     # fresh benchmark
bench/run_all.sh                       # skani presets, aif, ablations, budget sweep
python3 bench/time_all.py 11                                           # timing (interleaved, pinned)
python3 bench/tune.py data/tune <skani_tune_predictions.tsv> out.tsv  # parameter search
python3 scripts/benchmark_report.py; report_bench/build.sh             # this report

aif sketch ref.fa -o ref.aif           # defaults: -k auto, --budget 0.88, --tiers 0.005:0.6,0.1:0.4
aif query ref.aif contigs.fa [-t N] > out.tsv  # columns: contig ani af theta probes hits ... status
.CE
.LP
Benchmark data are regenerated from the scripts; summaries of all runs (evaluations, predictions,
timings) are committed under \f[CR]results/\f[R].
