.\" Build: report_overview/build.sh   (generated items: /usr/bin/python3 scripts/overview_report.py;
.\" tables and figures of the other reports are included from ../report_*/)
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
aif: Method and Results
.ft R
.ps 14
.vs 18p
Contig ANI and Alignment Fraction from FracMinHash-Sampled
Bloom-Filter Sketches, from Pairwise Queries to GTDB Scale
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Overview of method and results, aif v1.6 \(em November 23, 2026
Experiments and typesetting by Claude (Anthropic) for the project owner
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
aif estimates the average nucleotide identity (ANI) and the alignment fraction (AF) of short sequences
(contigs, MAG fragments, viral genomes) against reference genomes from compact sketches. The reference is
stored as FracMinHash-selected k-mers in two Bloom-filter tiers with deliberately high false-positive rates
(0.5% and 10%, 0.88 bits per reference base). The query is used raw: every selected k-mer is probed in
order, and a multi-level hidden Markov model over the hits gives the aligned fraction, a false-positive-corrected
containment and hence ANI, calibrated confidence intervals, and a segmentation of the contig into aligned
and unaligned regions with local identity. The theory behind the design shows that at a fixed budget many
noisy samples beat fewer exact ones (the bits-per-sample theorem), and gives the estimator's exact variance.
.sp 0.4v
This document collects the method and every result so far. At equal sketch size aif is more accurate than
skani 0.4.0 on simulated contigs (ANI 0.36 vs 1.12 pp, AF 0.39 vs 7.43 pp) and on 33,464 real contigs
(0.50 vs 0.68 pp, 2.6 vs 9.4 pp), and at 0.11\(en1.76 bits per base it is on the size\(enaccuracy front, ahead
of exact FracMinHash at up to 32\(mu the size. In database mode it searches all 199,923 GTDB R232 species
representatives (82.5 GB, the size of skani's database) 22\(mu faster than skani, and assigns contigs of held-out
genomes to the right species 78\(en79% of the time against 62\(en66%, confirmed on a fresh set of unseen species with
frozen defaults. It also gives the best vOTU clustering among sketch tools and identity maps more accurate than
fastANI's. Alignment (minimap2, vclust) remains more accurate at a much larger cost; the method's weak points
are contigs below 90% identity and the AF definition at that boundary.
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
.H1 "Introduction"
.LP
Comparing a genome fragment with a reference is a basic operation in microbial genomics: placing contigs and
MAGs, screening bins for contamination, dereplicating genomes, clustering viruses into vOTUs. Two numbers
matter, the identity of the aligned part (ANI) and the fraction of the fragment that aligns (AF). Alignment
gives both but needs the reference sequence. Sketch-based tools such as skani, fastANI and Mash compress the
reference, but their accuracy degrades on short and partially aligned fragments, where AF in particular is
hard to estimate.
.LP
aif ("ANI and AF") is built on three observations from the accompanying theory report
(\f[CR]fmh_bloom_filter_theory.pdf\f[R]):
.IP (1) 0.3i
\f[B]Bits per sample.\f[R] For a fixed number of bits per reference base, storing many FracMinHash samples in a
filter with a deliberately high false-positive rate (FPR) gives a lower-variance containment estimate than
storing fewer samples exactly. The optimal FPR is around 0.5\(en10%, far above usual practice.
.IP (2) 0.3i
\f[B]Exact variance.\f[R] The false-positive-corrected containment estimator is unbiased, and its variance
decomposes into mutation, sampling and false-positive terms.
.IP (3) 0.3i
\f[B]AF needs query positions only.\f[R] Unordered containment cannot separate AF from ANI, but the positions of
hits along the query make AF identifiable; positions on the reference are not needed.
.LP
The project has four detailed reports: the theory, a simulated benchmark of aif v1.1
(\f[CR]aif_benchmark_report.pdf\f[R]), real genomes, database mode and GTDB scale
(\f[CR]aif_real_genomes_report.pdf\f[R]), and viral clustering (\f[CR]aif_viral_clustering_report.pdf\f[R]).
This overview summarises the method (Section 2) and the results (Section 3); Table 1 lists the headline
numbers.
.TC 1
Results at a glance (aif v1.6 unless noted; MAE in percentage points (pp)). Details in Section 3 and the
detailed reports.
.TCE
.ps 9
.vs 11p
.so gen/tab_glance.tbl
.ps \n(PS
.vs \n(VSp
.\" =====================================================================
.H1 "Method"
.ne 3.3i
.PSPIC -C figs/fig_method.eps 6.5i
.FC 1
Overview of aif. The reference is sketched once; contigs are used raw. In database mode a contig is prepared
once and screened against all references, and only plausible candidates reach the model.
.FE
.H2 "Sketch"
.LP
Canonical k-mers of the reference are selected by FracMinHash: a k-mer x is kept when its 64-bit hash h(x)
falls below a threshold s\(md2\u64\d, about 12% of k-mers at the default budget. The selected hash range is
split into two tiers, each a cache-blocked Bloom filter (512-bit blocks) with its own target FPR: 0.5% for the
lower range, 10% for the upper. The tier boundaries follow from the budget (0.88 bits per reference base by
default) and the tier specification; they do not depend on the genome. Each tier's exact FPR \(*e is computed
from its fill ratio. The rate \(*s of spurious matches (non-homologous k-mers that occur by chance in a
genome-sized reference) is calibrated empirically from single-substitution neighbours of reference k-mers. For
pairwise use k is chosen per reference (17\(en21) so that \(*s \(<= 0.003; databases use one k (20).
The hedged design keeps a low-FPR tier for divergent contigs, where few k-mers survive, and a high-FPR tier
that maximises samples per bit.
.H2 "Containment and ANI"
.LP
Every query k-mer with h(x) < s is probed in its tier. With N probes and Y hits, the false-positive-corrected
containment is
.DS
	\f[I]C\f[R]  =  (Y/N \(mi \(*e\(fm) / (1 \(mi \(*e\(fm),        \(*e\(fm = \(*e + \(*s \(mi \(*e\(*s ,
.DE
computed per tier and combined with inverse-variance weights. Under the simple mutation model with ideal
hashing, Theorem 1 of the theory report gives E[\f[I]C\f[R]] = \(*h (the true containment of the aligned
part) exactly, and
.DS
	Var(\f[I]C\f[R])  =  Var(N\dmut\u)/n\u2\d  +  \(*h(1\(mi\(*h)(1\(mis)/(ns)  +  (1\(mi\(*h)\(*e\(fm/((1\(mi\(*e\(fm)ns) ,
.DE
the mutation-process, sampling and false-positive terms for n aligned k-mer positions. ANI is
\f[I]C\f[R]\u1/k\d with a delta-method bias correction. Theorem 3 compares an exact 64-bit FracMinHash tier
(sampling rate \(*b/64 at \(*b bits per base) with a filter tier at FPR \(*e (rate \(*b/c(\(*e)): the filter
divides the non-process variance by 64/J(\(*e), J(\(*e) = c(\(*e)[1 + \(*e/((1\(mi\(*e)\(*h)]. The gain is
4\(en23\(mu depending on identity, and the optimal FPR is near \(*h/ln(1/\(*h), which is the reason for the two tiers.
.H2 "Segmentation, AF and heterogeneous identity"
.LP
Containment alone cannot separate "half the contig aligned at high identity" from "all of it at lower identity"
(Proposition 3). aif therefore reads the probe outcomes in contig order with a hidden Markov model:
.IP \(bu 0.2i
\f[B]States.\f[R] One unaligned state (hits only by false positives, probability \(*e\(fm) and aligned
states at identity levels a from about 0.81 to 0.9995 in steps of 0.005 (hit probability
a\uk\d + (1\(mia\uk\d)\(*e\(fm).
.IP \(bu 0.2i
\f[B]Transitions.\f[R] Switching between aligned and unaligned costs a per-base rate (3\(mu10\u\(mi6\d), and
switching between identity levels a smaller one (10\u\(mi4\d).
.IP \(bu 0.2i
\f[B]Inference.\f[R] Forward\(enbackward over bins of 32 probes, with emissions from precomputed power tables, and
exact per-probe refinement around aligned/unaligned boundaries. AF is the posterior-weighted aligned length.
.IP \(bu 0.2i
\f[B]Heterogeneous contigs.\f[R] A test over 1 kb windows detects contigs whose identity varies (recombination,
mosaics). For them ANI is the posterior mean over identity levels, from a finer re-run (bins of 8 probes).
.LP
A contig is reported when the log-likelihood ratio against "all unaligned" is at least 10.
.H2 "Confidence intervals"
.LP
Statistical standard errors come from Theorem 1, conditioned on the realised substitutions (for ANI), and from
the posterior uncertainty of alignment boundaries (for AF). Real contigs violate the model (indels, repeats,
mosaics), so the default intervals are calibrated: stratified split-conformal multipliers (identity
\(mu length \(mu AF class, 24 strata) fitted on a disjoint real tuning set against MUMmer truth.
.H2 "Database mode"
.LP
\f[CR]aif db build\f[R] stores many references in one memory-mapped file with one k and one tier layout. Since
tier boundaries do not depend on the genome, a contig's selected k-mers and their filter hashes are the same
for every reference. A contig is therefore prepared once, and each candidate reference costs only its filter
tests and the HMM. An exhaustive search is byte-identical to pairwise queries.
.IP \(bu 0.2i
\f[B]Screening.\f[R] Candidates come from a bit-sliced Bloom filter, as in COBS. k-mers in the lowest 1/64 of
the hash range are inserted with 3 hash functions into one column per reference (8 bits per k-mer; 0.13\(en0.16
bits per base). References are sorted by size into groups of up to 512, so one 64-byte row answers a hash
function for 512 references.
.IP \(bu 0.2i
\f[B]Candidates.\f[R] A reference is a candidate when its shared-key count is significant (z \(>= 2), and
candidates are ranked by FP-corrected screening containment. Long contigs use their 128 keys of smallest hash.
.IP \(bu 0.2i
\f[B]Pre-HMM filter.\f[R] Up to 16 top candidates are evaluated. A two-state Viterbi score skips the model for
those that cannot reach the reporting threshold, without changing the output.
.IP \(bu 0.2i
\f[B]Output.\f[R] Hits are ordered by the model's likelihood ratio, and estimators are built lazily per
reference.
.H2 "Identity maps and identity-thresholded AF"
.LP
With \f[CR]--segments\f[R], aif writes the aligned and unaligned segments that tile each contig, and windowed
identity (about 1 kb) with calibrated 95% intervals. With \f[CR]--af-min-ident T\f[R], AF counts only regions
aligned at identity \(>= T, using the HMM's posterior over identity levels. This matches AF definitions based
on seed-and-extend alignment, such as megablast's for vOTUs.
.H2 "Implementation"
.LP
aif is one dependency-light C++17 file (zlib for gzip input), vectorised with AVX-512: k-mer packing and
selection, GFNI reverse complements, gather-based filter tests, and HMM recursions arranged to keep reductions off
the dependency chain. Threads use a streaming pipeline with ordered output.
.\" =====================================================================
.H1 "Results"
.H2 "Simulated contigs"
.LP
On 2,800 contigs simulated from four bacterial references under 28 conditions (identity, length, AF,
indels), aif v1.1 reports 95.1% of contigs (skani 86.8%) with 0.36 pp ANI and 0.39 pp AF error (skani 1.12 and
7.43 pp), and has lower ANI and AF error in all 28 conditions. A freshly generated set reproduces this
(Table 2).
.TC 2
Simulated benchmark (held-out set and a fresh confirmation set; aif v1.1, equal sketch size).
.TCE
.so ../report_bench/gen/tab_headline.tbl
.H2 "Real genomes, pairwise"
.LP
The real benchmark uses 240 assemblies (35,362 contigs \(>= 500 bp) from the skani paper's data:
B. cereus-group drafts against B. anthracis, E. coli drafts against K-12, MAGs of one species bin, P. aeruginosa
and S. aureus drafts, and negative controls. Truth is MUMmer alignment of every contig. No aif parameter was tuned
on these data. aif is more accurate than skani on every set for both ANI and AF (Table 3), and also than
skani's 3.8\(mu larger \f[CR]--slow\f[R] preset. On 0.5\(en1 kb contigs it reports 95% of homologous contigs
against 23%. Aggregated per genome, its contig estimates give the most accurate genome-level ANI (0.19 pp
against ANIm; skani 0.24, fastANI 0.25). Single-threaded it matches skani's speed on every set; at 8 threads it
is 2.4\(en3.6\(mu faster.
.TC 3
Real contigs, pooled over the five within-species sets (MUMmer truth; aif-RT: one setting chosen on a disjoint
real tuning set).
.TCE
.so ../report_real/gen/tab_headline.tbl
.H2 "Calibrated intervals"
.LP
Statistical intervals cover the truth at nominal rates on simulated substitution-only contigs, but are too
narrow against alignment truth on real contigs. Conformally calibrated intervals cover 94.4% (ANI) and 94.2%
(AF) of held-out real contigs at 95% (Table 4, Figure 2). In database mode they keep 94\(en96% coverage; on
GTDB's diverse taxa coverage drops to 92\(en93%, which measures the shift away from the five calibration species.
.TC 4
Coverage of statistical and calibrated intervals (median half-widths in pp).
.TCE
.so ../report_real/gen/tab_ci.tbl
.ne 3.0i
.PSPIC -C ../report_real/figs/fig_calibration.eps 6.4i
.FC 2
Empirical against nominal coverage of ANI and AF intervals.
.FE
.H2 "Size, accuracy and time"
.LP
Each method was swept across sketch sizes on the real sets (Table 5, Figure 3). aif is more accurate than
skani and than exact FracMinHash (sourmash) at every size, and faster up to 0.88 bits per base. Exact
FracMinHash at 3.56 bits per base (0.71 pp ANI error) is less accurate than aif at 0.11 bits per base
(0.67 pp): the bits-per-sample theorem on real data. Alignment (minimap2) reaches 0.23 pp ANI and 1.78 pp AF
error, with an index of 27 bits per base and 210\(mu aif's run time. aif's AF error at 1.76 bits per base
(1.73 pp) matches it.
.ne 4.6i
.TC 5
Pairwise accuracy against sketch (or index) size and time (pooled within-species contigs).
.TCE
.so ../report_real/gen/tab_pareto.tbl
.ne 3.0i
.PSPIC -C ../report_real/figs/fig_pareto.eps 6.6i
.FC 3
Error against sketch size and time.
.FE
.H2 "Database search: 2,005 genomes"
.LP
Against a database of 2,005 complete genomes (8.6 Gbp: the pairwise references, 100 genomes of each query
species and 1,600 other RefSeq reference genomes), aif's top hit has 0.31 pp ANI and 1.12 pp AF error over
26,084 homologous contigs. skani search has 0.51 / 3.03 pp, and skani dist 0.69 / 4.15 pp (Table 6). aif
reports 99.9% of homologous contigs (77.0% and 94.5%), its top hit is the best evaluated reference for 88.9%
(78.3%, 74.9%), and it searches 33,464 contigs in 4.4 s at 8 threads (skani search 37.9 s). A size-matched aif
database (0.895 bits per base, the size of skani's) is almost as accurate.
.TC 6
Top-hit accuracy against the 2,005-genome database (MAE and regret in pp).
.TCE
.so ../report_real/gen/tab_db_acc.tbl
.H2 "GTDB scale"
.LP
The GTDB database holds all 199,923 species representatives of release R232 (655 Gbp). With the genomes in
the page cache aif builds it in 6.7 min with 1.9 GB of memory (skani: 10.8 min, 15.6 GB), and the databases are
the same size (82.5 vs 82.7 GB). On 23,561 contigs of 294 held-out non-representative genomes, aif's search
takes 11.5 s at 16 threads against skani's 252 s (Table 7).
.LP
aif's top hit has the query's species for 78.3% of contigs and its genus for 88.5%; skani's for 61.6% and 68.9%
(Table 8, Figure 4). aif is more accurate at every contig length, and under an equal calling rule (ANI \(>= 95%)
it calls 76.6% of contigs correctly at 86.4% precision (skani 60.7% at 87.2%). A 28.5 GB aif database
(budget 0.22) is nearly as accurate.
.TC 7
Building and searching the GTDB R232 database.
.TCE
.so ../report_real/gen/tab_gtdb_build.tbl
.TC 8
Contig classification against GTDB R232 (query set 1).
.TCE
.so ../report_real/gen/tab_gtdb_acc.tbl
.ne 2.8i
.PSPIC -C ../report_real/figs/fig_gtdb.eps 6.4i
.FC 4
Species-level accuracy and top-hit ANI error against contig length (GTDB R232, query set 1).
.FE
.H2 "Held-out confirmation with frozen defaults"
.LP
Several database defaults were chosen while looking at the benchmark results: the screening setup, the
pre-HMM filter cutoff, the key cap and the hit order. They were frozen (git tag \f[CR]v1.6-freeze\f[R]), and
two fresh query sets were drawn and scored once.
.IP \(bu 0.2i
\f[B]2,005-genome database, 80 new drafts.\f[R] aif's top-hit errors are 0.23 pp ANI and 0.77 pp AF, against
0.41 / 2.26 pp for skani search and 0.52 / 2.93 pp for skani dist. Coverage is 99.9%, against 80.2% and 92.4%.
.IP \(bu 0.2i
\f[B]GTDB query set 2.\f[R] 298 genomes from 300 species not used before (Table 9). aif's species accuracy is
78.6% (skani 65.8%) and its genus accuracy 90.1% (73.2%).
.LP
Every conclusion holds on unseen genomes and species.
.TC 9
GTDB R232: query set 1 and the fresh set 2 (disjoint species), frozen defaults.
.TCE
.so ../report_real/gen/tab_fresh_gtdb.tbl
.H2 "Contamination screening"
.LP
Bins were simulated from the GTDB queries: each genome plus 5\(en20% of another genome's contigs. A contig is
flagged when its top hit's genus differs from the bin's majority. aif flags 91% (set 1) and 93% (set 2) of
contaminant contigs at 75\(en77% precision; skani flags 68% and 70% at 81\(en84%. On contigs of 0.5\(en2 kb
aif's recall is 76\(en84% and skani's 13\(en15%, because skani reports nothing for most short contigs.
.H2 "Identity maps"
.LP
Against MUMmer alignments of 224 query genomes, aif's window identity at 3 kb has 0.55 pp error and
correlation 0.97 with the truth. fastANI's 3 kb fragments have 0.78 pp and 0.86 (Table 10, Figure 5). fastANI
is better below 90% identity. Calibrated window intervals cover 94.4\(en94.7% of held-out windows. aif's
unaligned segments find unaligned islands at 97% precision and 81% recall (fastANI 56% / 51%), and
aligned/unaligned breakpoints within a median 171 bp.
.TC 10
Window identity against MUMmer (windows with true AF \(>= 0.9).
.TCE
.so ../report_real/gen/tab_idmap.tbl
.ne 3.0i
.PSPIC -C ../report_real/figs/fig_idmap.eps 6.6i
.FC 5
Identity map of a 162 kb P. aeruginosa draft contig against PAO1.
.FE
.H2 "Viral clustering"
.LP
vOTU clustering uses \(>= 95% ANI over \(>= 85% AF, with BLAST and CheckV's ANI/AF as the standard. On
10,000 INPHARED phage genomes and 10,000 Gut Phage Database contigs (Table 11), aif's ANI is the most accurate
of the sketch tools near the thresholds (0.57\(en0.63 pp; skani 0.89\(en1.10 pp). Its default AF over-estimates
BLAST's on partial overlaps, because megablast does not align the divergent stretches the HMM counts. With AF
counted at identity \(>= 0.88 (chosen on a disjoint tuning set), aif has the best vOTU agreement among sketch
tools. The alignment-based vclust is best overall, at 7\(en10\(mu aif's run time.
.TC 11
vOTU edges and clusters against BLAST/CheckV (aif with \f[CR]--af-min-ident 0.88\f[R] in the second aif row).
.TCE
.so ../report_viral/gen/tab_votu.tbl
.\" =====================================================================
.H1 "Limitations and open items"
.IP \(bu 0.2i
\f[B]Divergent contigs.\f[R] Below about 90% identity aif misses pairs, over-estimates ANI by about 2.6 pp
against MUMmer, and sometimes over-calls AF. Positional anchors on the reference would allow colinearity checks.
.IP \(bu 0.2i
\f[B]AF definition.\f[R] aif's AF measures homology down to its \(ap84% detection floor. Matching an alignment
convention (MUMmer, megablast) is a choice, made explicit by \f[CR]--af-min-ident\f[R].
.IP \(bu 0.2i
\f[B]Alignment is more accurate.\f[R] minimap2 on bacterial contigs and vclust on viral genomes are more
accurate, at a much larger index or run time. aif's niche is accuracy per stored bit, uncertainty and
segmentation at database scale.
.IP \(bu 0.2i
\f[B]Baselines.\f[R] Contig classification is compared with skani only. Kraken2, MMseqs2 taxonomy and CAT, and
GUNC or CheckM2 for contamination, are the next comparisons (in progress).
.IP \(bu 0.2i
\f[B]Biology.\f[R] All applications so far are benchmarks; no biological finding yet.
.IP \(bu 0.2i
\f[B]Engineering.\f[R] The fast paths need AVX-512, and timings come from one host. There is no packaging yet.
Search memory is roughly the database size (memory-mapped).
.IP \(bu 0.2i
\f[B]Calibration scope.\f[R] Intervals were calibrated on five bacterial species against MUMmer; coverage is
92\(en93% on GTDB's diverse taxa.
.\" =====================================================================
.H1 "Reports, code and reproducibility"
.LP
Everything is in the repository \f[CR]fmh-augmenting-bloom-filter\f[R]: the tool (\f[CR]src/aif.cpp\f[R]),
benchmark scripts (\f[CR]bench/\f[R], \f[CR]bench_real/\f[R]), report generators (\f[CR]scripts/\f[R]) and
committed result summaries (\f[CR]results/\f[R], \f[CR]results_real/\f[R], \f[CR]results_viral/\f[R]). The
detailed reports give full methods, data selection, all tables, and the commands to reproduce them:
.IP \(bu 0.2i
\f[CR]fmh_bloom_filter_theory.pdf\f[R]: variance framework, designs, theorems and proofs;
.IP \(bu 0.2i
\f[CR]aif_benchmark_report.pdf\f[R]: aif v1.1 on simulated contigs;
.IP \(bu 0.2i
\f[CR]aif_real_genomes_report.pdf\f[R]: real genomes, speed, intervals, database mode, GTDB, contamination,
size\(enaccuracy front, identity maps, held-out confirmation;
.IP \(bu 0.2i
\f[CR]aif_viral_clustering_report.pdf\f[R]: vOTU clustering.
.LP
Software versions: aif v1.6 (git tag \f[CR]v1.6-freeze\f[R] for the database defaults), skani 0.4.0, sourmash
4.8.14, minimap2 2.26, fastANI, vclust 1.3.1, MUMmer 4, BLAST+ 2.12. Timings were measured on one Xeon Gold
5416S host.
