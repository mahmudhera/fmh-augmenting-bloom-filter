.\" Build: report_real/build.sh   (tables/figures: python3 scripts/real_report.py)
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
aif on Real Genomes
.ft R
.ps 14
.vs 18p
Contig ANI and Alignment Fraction on Draft Assemblies and MAGs
from the skani Benchmark Data
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Real-genome report, aif v1.5 \(em November 22, 2026
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
The earlier benchmark used contigs simulated from four reference genomes. Here aif is tested on real genome pairs.
The data are 240 real assemblies with 35,362 contigs of at least 500 bp (1.13 Gbp), drawn from the data sets of the
skani paper. They are the B. cereus group (D2) and E. coli (D3) draft assemblies, metagenome-assembled genomes of
species bin 2328, GC-extreme P. aeruginosa and S. aureus drafts, and cross-genus negative controls.
The ground truth is MUMmer alignment of every contig to its reference.
No aif parameter was tuned on these data. At equal sketch size (0.88 vs 0.90 bits per reference base) aif is better
than skani 0.4.0 on every set, for both ANI and AF.
Pooled over the 33,464 within-species contigs, aif reports 99.2% of homologous contigs (skani 90.1%).
Its ANI mean absolute error is 0.50 pp (skani 0.68 pp), or 0.39 vs 0.68 pp on the 23,488 contigs that every
configuration reports. Its AF error is 2.63 pp (skani 9.35 pp).
aif also beats skani's 3.8\(mu larger \f[CR]--slow\f[R] preset on every set. At half its budget, aif still beats
skani's default in pooled ANI (0.56 vs 0.68 pp) and in AF on every set.
On 500\(en1000 bp contigs aif reports 95% of homologous contigs, against 23% for skani, with an AF error of 4.3 vs 39 pp.
Aggregating aif's contig estimates also gives the most accurate genome-level ANI against ANIm (0.19 pp vs 0.24 pp for
skani and 0.25 pp for fastANI).
One setting change, chosen on a disjoint real tuning set (aif-RT), lowers the held-out AF error to 1.85 pp.
Real assemblies also exposed speed problems. Version 1.2 made the thread pipeline scale on Mbp-scale contigs.
Version 1.3 then closed the single-thread gap with:
GFNI reverse complements; vectorised 2-bit packing and bin counting; and HMM recursions restructured so that their
reductions leave the dependency chain.
Single-threaded, v1.3 processed the five sets in 2.81 s against skani's 3.28 s, faster on every set; at 8 threads,
0.89 s against 2.66 s. Estimates moved by at most 10\u\s-2\(mi6\s0\d.
Version 1.4 adds per-contig confidence intervals for ANI and AF. Statistical intervals from Theorem 1 cover the
truth at nominal rates when the model holds (96.5% at 95% on simulated substitution-only contigs).
Intervals calibrated by stratified split-conformal prediction on the real tuning set cover 94.4% (ANI) and 94.2% (AF)
of held-out real contigs at 95%, with median half-widths of 0.87 and 2.6 pp. They cost about 4% in run time, and aif
still matches or beats skani single-threaded on every set (2.92 vs 3.27 s).
We also traced aif's AF over-estimate on partially aligned contigs to two mechanisms, and tested principled
corrections. One of them, a local identity prior, is offered as an option. Each correction trades contig AF
against genome-level ANI, so the defaults are unchanged.
Version 1.5 adds a database mode. All references share one k, so a contig is hashed once and each candidate
reference costs only its filter tests and the HMM. A bit-sliced screening index chooses the candidates, and an
exhaustive search is byte-identical to pairwise queries. Against 2,005 complete genomes (8.6 Gbp), aif reports
99.9% of homologous contigs, against 77% for skani search and 95% for skani dist. Its top-hit ANI error is
0.32 pp (skani 0.51 and 0.69 pp), its AF error 1.3 pp (3.0 and 4.2 pp), and its top hit is the best reference for
90% of contigs (77% and 74%). At 8 threads it searches 33,464 contigs in 5.2 s, against 37.9 s for skani search.
Its database is 1.04 bits per base against skani's 0.90, and a size-matched aif database is as accurate.
We also report where aif errs and how the truth was audited.
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
The earlier benchmark (\f[CR]aif_benchmark_report.pdf\f[R]) showed that aif beats skani at equal sketch budget on
contigs simulated from four reference genomes. This report asks whether that holds on real genome pairs.
Real assemblies differ from simulations in several ways: their divergence varies along the genome; they contain
mosaic and accessory segments, repeats and mobile elements; contig lengths are very uneven; and MAGs add
fragmentation and contamination.
Table 1 gives the headline. aif's numbers come from parameters tuned on the simulated tuning set only. aif-RT
changes one setting, chosen on a disjoint real tuning set (Section 6).
.TC 1
Headline, pooled over the five within-species sets. ANI is scored on homologous contigs (true AF \(>= 0.2)
that the tool reports. AF is scored on all contigs, with an unreported contig counted as AF 0. The false-homology
column counts contigs with true AF < 0.05 that are reported with AF \(>= 0.2.
.TCE
.so gen/tab_headline.tbl
.LP
The main findings:
.IP (1) 0.3i
\f[B]Contig ANI.\f[R] aif has the lowest ANI error on every set (Table 4). Pooled, it is 0.50 pp against 0.68
(skani) and 0.61 (skani \f[CR]--slow\f[R]), while reporting more contigs. On contigs that every configuration
reports, the error is 0.39 vs 0.68 pp.
.IP (2) 0.3i
\f[B]Contig AF.\f[R] aif's AF error is 2.6 pp against 9.4 (skani) and 4.3 (skani \f[CR]--slow\f[R]).
skani under-estimates AF, worst on short contigs. aif slightly over-estimates AF on partially aligned, divergent
contigs (Section 5).
.IP (3) 0.3i
\f[B]Short contigs.\f[R] On 0.5\(en1 kb contigs, skani reports 23% of homologous contigs and aif 95%. AF error
is 39 vs 4.3 pp, and paired ANI error 0.92 vs 0.49 pp.
.IP (4) 0.3i
\f[B]Genome level.\f[R] Aggregated per-contig aif output matches ANIm best (0.190 pp MAE over 224 genomes).
skani genome mode gives 0.237 pp and fastANI 0.247 pp. skani \f[CR]--slow\f[R] has the best genome AF (Section 4).
.IP (5) 0.3i
\f[B]Speed.\f[R] Single-threaded, aif v1.4 (sketch + query, with intervals) processes the five sets in 2.92 s
against skani's 3.27 s, matching or beating it on every set; at 8 threads, 0.86 s against 2.72 s. Version 1.1 was
20\(en55% slower than skani single-threaded and did not scale with threads on these inputs (Section 7).
.IP (6) 0.3i
\f[B]Uncertainty.\f[R] Every estimate now has a confidence interval. Calibrated intervals cover held-out real
contigs at the nominal rate (94.4% ANI, 94.2% AF at 95%); skani reports no uncertainty (Section 8).
.IP (7) 0.3i
\f[B]AF model.\f[R] Uniform identity jumps let unaligned islands be absorbed as low-identity alignment. A local
identity prior fixes this and lowers held-out contig AF error by 13%, but raises genome-level ANI error. No
correction we tested improves both, so the defaults are unchanged (Section 9).
.IP (8) 0.3i
\f[B]Database mode.\f[R] Against 2,005 genomes, aif's top hit has 0.32 pp ANI and 1.3 pp AF error, against
0.51 / 3.0 pp for skani search and 0.69 / 4.2 pp for skani dist. aif covers 99.9% of homologous contigs
(77% and 95%), and is 7\(mu faster than skani search, 3\(mu faster than its \f[CR]--keep-refs\f[R] mode. It uses
more memory (2.0 vs 1.3 GB) and a 16% larger database; a size-matched database is as accurate (Section 10).
.\" =====================================================================
.H1 "Data and ground truth"
.H2 "Genome sets"
.LP
All sets come from, or follow, the skani paper's benchmark data (Shaw and Yu, Nature Methods 2023; data and scripts in
github.com/bluenote-1577/skani-test).
.IP \(bu 0.2i
\f[B]B. cereus group (D2).\f[R] This is the fastANI/skani D2 set: draft assemblies queried against B. anthracis 52-G.
We drew 80 draft assemblies, stratified by the published ANIm value. There are 18 each in 85\(en90, 90\(en95 and
95\(en98%, all 8 in 98\(en99.5%, and 18 at \(>= 99.5%.
The reference is the 52-G chromosome (NZ_CM002395.1).
.IP \(bu 0.2i
\f[B]E. coli (D3).\f[R] 48 random draft assemblies from the fastANI/skani D3 list, against E. coli K-12 MG1655.
.IP \(bu 0.2i
\f[B]MAG bin 2328.\f[R] These are 48 random metagenome-assembled genomes of species-level bin 2328 from Pasolli et al.
(the 195-MAG set of skani's clustering figure). The reference is the most contiguous MAG of the bin (N50 597 kb),
excluded from the queries.
.IP \(bu 0.2i
\f[B]GC extremes.\f[R] 24 GenBank contig-level P. aeruginosa assemblies (66% GC) against PAO1, and 24 S. aureus
assemblies (33% GC) against NCTC 8325. These references are also used in the simulation benchmark.
.IP \(bu 0.2i
\f[B]Negative controls.\f[R] 8 E. coli drafts against B. anthracis, and 8 B. cereus-group drafts against E. coli.
.LP
Contigs shorter than 500 bp are removed, and every remaining contig is a query. Assemblies come from NCBI.
When an accession has several versions, the newest is used. Selection is deterministic (seed 20261121;
\f[CR]bench_real/select_datasets.py\f[R]). Table 2 summarises the sets.
Divergence is real: B. cereus-group genomes span 86\(en100% ANIm, and one P. aeruginosa assembly is at 93.6%
(a P. paraeruginosa-like strain).
.TC 2
Real genome sets. "Homol." is the fraction of contigs with true AF \(>= 0.2. "Genome ANIm" is MUMmer 1-to-1 identity.
.TCE
.so gen/tab_data.tbl
.H2 "Alignment-based truth"
.LP
Every query genome is aligned to its reference with MUMmer 4 (\f[CR]bench_real/truth.py\f[R]).
.IP \(bu 0.2i
\f[B]Per contig:\f[R] \f[CR]nucmer --maxmatch\f[R], then \f[CR]delta-filter -q\f[R], which keeps the best
reference placement of every query base.
True AF is the union of aligned query intervals divided by contig length.
True ANI is the length-weighted identity of those alignments; indels count as differences, as in ANIm.
.IP \(bu 0.2i
\f[B]Per genome:\f[R] \f[CR]nucmer\f[R] with default seeds and \f[CR]delta-filter -1\f[R], as in
\f[CR]dnadiff\f[R], giving genome-level ANIm and the 1-to-1 aligned fraction.
.LP
\f[CR]--maxmatch\f[R] matters for per-contig truth. With default (unique) anchors, contigs that are entirely a
multi-copy element have no unique seed and stay unaligned, even though their sequence is in the reference.
IS-element contigs of 1.3\(en1.4 kb in S. aureus and P. aeruginosa are examples (Section 5.1).
As a check, our genome-level ANIm reproduces the published D2 ANIm values (skani-test \f[CR]D2_anim.txt\f[R]).
The median absolute difference is 0.05 pp and the mean 0.21 pp over the 80 genomes.
The three largest differences (2\(en4.7 pp) are assemblies whose GenBank record has since been replaced by a
new version.
.H2 "Tools, budgets and metrics"
.LP
skani 0.4.0 runs as \f[CR]skani dist --qi\f[R] (per-contig mode) with default parameters (c = 125), with
\f[CR]--slow\f[R] (c = 30), and with relaxed reporting thresholds (\f[CR]-s 70 --min-af 5\f[R]; nearly identical to
default, so omitted from most tables).
aif runs with its v1.1 defaults. All accuracy numbers are v1.4 runs. Default estimates differ from v1.1\(en1.3 by
at most 10\u\s-2\(mi6\s0\d (Section 7.3).
Sketch sizes are those of the earlier report (Table 3). aif uses 0.880 bits per reference base and skani 0.903.
.TC 3
Sketch size per reference (skani: all files of \f[CR]skani sketch\f[R]; aif: the \f[CR].aif\f[R] file).
.TCE
.so gen/tab_sizes.tbl
.LP
The metrics are:
.IP \(bu 0.2i
\f[B]ANI\f[R] is scored on homologous contigs (true AF \(>= 0.2). Coverage is the fraction of them that a tool
reports. The mean absolute error (MAE, in percentage points) is over reported contigs. The paired MAE uses only
contigs that every compared configuration reports (23,488 contigs).
.IP \(bu 0.2i
\f[B]AF\f[R] is scored on all contigs. A tool that reports nothing for a contig is taken to claim no
homology (AF = 0). Without this rule, a tool could lower its AF error by staying silent.
.IP \(bu 0.2i
\f[B]False homology calls\f[R] are non-homologous contigs (true AF < 0.05) reported with AF \(>= 0.2.
.LP
Pooled numbers cover the five within-species sets (33,464 contigs); macro numbers average over the five sets.
.\" =====================================================================
.H1 "Contig-level accuracy"
.H2 "Per set"
.LP
aif has lower ANI and AF error than both skani presets on all five sets (Table 4, Figure 1). At half its budget
(0.44 bits/base) aif still has the lower AF error on every set and the lower ANI error on three. It is worse on
E. coli (0.72 vs 0.62 pp) and level on P. aeruginosa. The AF gap is largest
on the B. cereus group and S. aureus: skani's errors there are 14.0 and 14.7 pp, aif's 3.8 and 1.8 pp.
aif also covers more homologous contigs everywhere, and its paired ANI error is lower everywhere (Table 5).
.TC 4
ANI MAE (pp, homologous contigs) and AF MAE (pp, all contigs) per set. The negative controls have no homologous
contigs, so only AF is shown.
.TCE
.so gen/tab_sets.tbl
.TC 5
Coverage of homologous contigs and paired ANI MAE (contigs reported by all configurations).
.TCE
.so gen/tab_cov.tbl
.ne 3.0i
.PSPIC -C figs/fig_sets.eps 6.4i
.FC 1
ANI and AF mean absolute error per set ("skani slow" is skani \f[CR]--slow\f[R], 3.8\(mu aif's sketch size).
.FE
.H2 "By contig length"
.LP
skani's contig mode needs enough seeds. On 0.5\(en1 kb contigs it reports 23% of homologous contigs (58% with
\f[CR]--slow\f[R]), and its AF error is 39 pp (Table 6, Figure 2).
aif reports 95% of them with a 4.3 pp AF error, and its paired ANI error there (0.49 pp) is about half of skani's.
At 20 kb and longer the tools converge in coverage. aif's ANI error is still about half of skani's
(0.24 vs 0.43 pp), and its AF error is lower than skani's but slightly above skani \f[CR]--slow\f[R]'s.
.TC 6
Accuracy by contig length (cov = coverage of homologous contigs; ANI and AF MAE in pp).
.TCE
.so gen/tab_len.tbl
.ne 2.7i
.PSPIC -C figs/fig_length.eps 6.4i
.FC 2
Coverage, paired ANI MAE and AF MAE by contig length.
.FE
.H2 "By divergence and by aligned fraction"
.LP
Table 7 stratifies homologous contigs by true ANI. Table 8 stratifies all contigs by true AF, and Table 9 gives the
signed errors.
aif's ANI error is lowest in every stratum from 90% upward (for example 1.24 vs 1.91 pp at 90\(en95%).
Below 90% all tools are inaccurate (MAE about 2.6\(en2.7 pp). Here aif reports 81% of contigs and skani 36%, and
aif over-estimates ANI by 2.6 pp on average.
The likely reason is a difference in what is called aligned. At this divergence nucmer aligns through
highly divergent stretches and counts them in the identity. aif's HMM classifies the most divergent stretches as
unaligned, so its identity describes the better-conserved part.
.LP
For AF, skani under-estimates at every level (\(mi5 to \(mi47 pp). It reports AF only where its sparse chains
reach, and its 15% reporting threshold zeroes partial contigs.
aif is nearly unbiased on fully aligned contigs (+0.6 pp at AF \(>= 0.9). It over-estimates partially aligned
contigs by about 7 pp: an unaligned island between aligned stretches is absorbed into a low-identity aligned
state when evidence is weak. aif-RT cuts this bias to +2.8\(en3.6 pp (Section 6).
.TC 7
Accuracy by true contig ANI (homologous contigs).
.TCE
.so gen/tab_ani.tbl
.TC 8
Accuracy by true AF. Contigs with AF < 0.05 contribute only AF error.
.TCE
.so gen/tab_af.tbl
.TC 9
Mean signed error (estimate minus truth) by true ANI.
.TCE
.so gen/tab_bias.tbl
.ne 2.4i
.PSPIC -C figs/fig_af_scatter.eps 6.4i
.FC 3
Estimated versus true AF for all 33,464 within-species contigs (hexagonal bins, log density). Points on the
bottom edge are contigs that the tool did not report (AF 0).
.FE
.ne 2.7i
.PSPIC -C figs/fig_ani_err.eps 6.4i
.FC 4
ANI error versus true ANI for reported homologous contigs (log density). skani's upper edge is its
100% cap. Both tools' errors grow with divergence; aif's are concentrated within about \(+-1 pp down to 95%.
.FE
.H2 "Non-homologous contigs and negative controls"
.LP
Of 6,404 within-species contigs with true AF < 0.05, aif reports 14 with AF \(>= 0.2 (0.22%) and skani 2 (0.03%).
All 14 of aif's calls are weak: likelihood ratios 10\(en40, with ANI 87\(en92%. In two of them skani \f[CR]--slow\f[R]
also finds 28% AF, so some may be divergent homology that nucmer misses.
aif reports a small AF (< 0.05) on a further 4% of these contigs. These are real short homologies, such as IS
fragments or rRNA, which the AF error already scores.
On the 1,898 cross-genus negative controls, aif reports 4 contigs (AF MAE 0.41 pp) and skani 13 (0.18 pp).
Nearly all of the negative-control truth is AF 0, and both tools are essentially silent.
.\" =====================================================================
.H1 "Genome-level ANI and AF"
.LP
skani and fastANI are designed for genome-versus-genome ANI. aif estimates per contig, but a genome-level value
follows by aggregation:
.IP \(bu 0.2i
genome AF = \(*S af\d\s-2i\s0\u L\d\s-2i\s0\u / \(*S L\d\s-2i\s0\u;
.IP \(bu 0.2i
genome ANI = the aligned-length-weighted mean of contig ANI.
.LP
Table 10 compares this aggregate (and the same aggregation of skani \f[CR]--qi\f[R]) with skani genome mode and
fastANI. The truth is the dnadiff-style ANIm of Section 2.2.
aif's aggregate has the lowest genome ANI error overall: 0.190 pp, against 0.237 for skani, 0.259 for
skani \f[CR]--slow\f[R] and 0.247 for fastANI. On the divergent B. cereus set all methods are within
0.34\(en0.38 pp, except aif-RT (0.54) and aif fast (0.42).
skani \f[CR]--slow\f[R] gives the best genome AF, at 3.8\(mu the space. aif's aggregated AF (2.0\(en2.4 pp)
is better than skani default's (3.6\(en4.0 pp) and fastANI's.
aif-RT has a better genome AF (1.0 pp many-to-one) but a worse genome ANI on the divergent B. cereus set.
Its higher switching rate drops the most divergent stretches of long contigs from the aligned set.
For genome-level ANI the default is therefore the better setting.
.ne 2.3i
.TC 10
Genome-level accuracy against MUMmer ANIm. AF is scored against the 1-to-1 aligned fraction (dnadiff) and the
many-to-one fraction (the per-contig truth aggregated). "agg." means aggregated per-contig output.
.TCE
.so gen/tab_genome.tbl
.\" =====================================================================
.H1 "How reliable is the truth?"
.H2 "Repeats and seeding"
.LP
Our first truth used nucmer's default unique-anchor seeding. With it, aif appeared to report 6.4% of
non-homologous contigs. Inspection showed that many were 1.3\(en1.4 kb contigs with 100% of their k-mers present
in the reference: whole IS elements, multi-copy in the reference, which nucmer could not seed.
Re-aligning with \f[CR]--maxmatch\f[R] raised the true AF by more than 5 pp for 302 contigs. It lowered aif's AF
error from 2.93 to 2.63 pp and left skani's essentially unchanged (Table 11).
Under either truth, aif is ahead by a wide margin.
.TC 11
Effect of the per-contig truth's seeding on the scores (pooled within-species contigs).
.TCE
.so gen/tab_truth.tbl
.H2 "A blastn audit"
.LP
In 887 contigs (2.7%), aif's AF exceeds the MUMmer AF by more than 0.2.
We re-aligned 150 of them, and 150 controls, with the more sensitive blastn (\f[CR]-task blastn\f[R], word
size 11). We took the union of all HSPs on each contig, a generous upper bound on alignable sequence
(\f[CR]bench_real/blast_audit.py\f[R]; Table 12).
On controls, nucmer, blastn and aif agree to within a few points.
On disputed contigs, blastn finds 13 pp more homology than nucmer (0.61 vs 0.48), but aif still reports much more
(0.87). Both explanations hold, then. nucmer under-aligns divergent sequence, which accounts for about a third of
the disagreement. aif really does over-call AF on a few percent of contigs, mostly those it places at 85\(en92% identity.
Dispersed repeat k-mers can mimic a low-identity alignment, and aif's sketch has no positions to check
colinearity.
.TC 12
blastn audit of the MUMmer truth (mean AF; mean absolute difference in pp).
.TCE
.so gen/tab_audit.tbl
.H2 "Other caveats"
.LP
The truth's identity counts indels as differences, as ANIm does, and skani's ANI model was trained to reproduce
MUMmer identity. If anything, MUMmer truth favours skani. Each query is compared with one reference; database
search is not tested.
.\" =====================================================================
.H1 "Real-data tuning (aif-RT)"
.LP
To ask whether real data call for different settings without tuning on the evaluation genomes, we built a disjoint
real tuning set (seed 777; \f[CR]bench_real/select_tune.py\f[R]). It draws from the same sources and references,
with none of the evaluation genomes: 24 B. cereus-group, 20 E. coli, 20 MAG, 10 P. aeruginosa and 10 S. aureus
assemblies, 12,416 contigs in all. The MAG set has its own reference MAG.
On it, aif v1.1 already scores like the evaluation set (ANI 0.486, AF 2.62 pp; skani 0.622 and 9.72 pp).
.LP
A grid over query-time parameters (\f[CR]bench_real/tune_real.py\f[R]) showed one strong effect: the
aligned/unaligned switching prior \f[CR]--rho-u\f[R] (Table 13).
Real contigs are more mosaic than the simulation, and a higher switching rate lets the HMM close unaligned
islands. The other parameters (\f[CR]--rho-h\f[R], \f[CR]--ident-c\f[R], \f[CR]--lr-min\f[R], \f[CR]--het-z\f[R])
gave no improvement without costing coverage or ANI.
The rule was fixed before looking at held-out data: the lowest AF error among settings whose ANI error is no worse
than the default's. It selects \f[CR]--rho-u 2e-5\f[R] (aif-RT).
On the held-out evaluation set, aif-RT lowers the AF error from 2.63 to 1.85 pp pooled (2.45 to 1.71 macro), and on
partially aligned contigs from 14.4 to 9.5 pp. ANI changes from 0.502 to 0.496 pp.
.LP
The cost shows on the simulated tuning set (Table 13, right). There AF error rises from 0.45 to 0.72 pp, mostly
on uniformly divergent (80\(en90%) contigs, whose construction truth is AF = 1 even where no aligner would align.
The setting also worsens genome-level ANI on divergent genomes (Section 4).
We therefore keep the simulation-tuned default. We recommend \f[CR]--rho-u 2e-5\f[R] when per-contig AF on
real assemblies is the main output.
.TC 13
Effect of the switching prior rho_u on the two tuning sets (MAE in pp; real set: macro over 5 sets; simulated set:
all reported contigs).
.TCE
.so gen/tab_tune.tbl
.\" =====================================================================
.H1 "Speed and memory on real data (v1.2\(en1.4)"
.H2 "What real assemblies exposed"
.LP
The simulated contigs were at most 20 kb. Draft assemblies contain contigs up to 5 Mbp, and that broke two
assumptions of v1.1.
.IP (1) 0.3i
\f[B]Thread pipeline.\f[R] v1.1 read 4 Mbp batches, handed out 16 contigs per grab, and waited at a barrier after
each batch. A batch with a few Mbp-scale contigs kept one thread busy while the others idled. Every thread also
grew its scratch buffers to the largest contig. On the 449 Mbp B. cereus set, v1.1 took 1.9 s at 1 thread and
3.5 s at 8 threads, using 1.6 GB at 8 threads (Table 15).
.IP (2) 0.3i
\f[B]Cache behaviour.\f[R] Probing a 5 Mbp contig in one pass streams megabytes of packed sequence and positions
through L2. Per-stage cycle counters (\f[CR]make prof\f[R]) showed probing at 6.9 cycles per base on \(>= 500 kb
contigs, against 4.2 on 20\(en500 kb contigs.
.H2 "Changes in v1.2"
.IP \(bu 0.2i
\f[B]Windowed probing.\f[R] Contigs are probed in windows of 65,536 k-mer starts. Output is unchanged, and the
per-thread scratch no longer depends on contig length.
On \(>= 500 kb contigs the probing cost fell from 1.27 to 0.69 Gcycles.
.IP \(bu 0.2i
\f[B]Streaming pipeline.\f[R] The reader queues chunks: one long contig, or up to 64 short ones totalling 256 kb.
Workers take chunks as they free up, and a writer thread prints them in input order. Outstanding input is capped
at 128 MB.
For memory-mapped input the reader only locates records. Workers strip line breaks (AVX-512 VBMI2
\f[CR]vpcompressb\f[R]) into reused buffers. Written chunks are recycled, and mapped pages are released behind
the writer.
.IP \(bu 0.2i
\f[B]Optional speed switch.\f[R] \f[CR]--het-fine-max N\f[R] skips the fine-binned (8-probe) heterogeneity
re-run for contigs with more than N probes. The default is off.
With N = 16384, the B. cereus set runs 22% faster. Contig-level ANI is essentially unchanged (Table 10 row
"aif fast"), but genome-level ANI on divergent genomes is worse (0.415 vs 0.373 pp). Accuracy comes first, so the
switch is not the default.
.LP
v1.2 output is byte-identical to v1.1 at any thread count. This was checked on all real sets, the simulated
benchmark, and piped input. On the simulated throughput sets v1.2 is as fast as or faster than v1.1 (20 Mbp at
8 threads: 33 vs 49 ms).
.H2 "Matching skani single-threaded (v1.3)"
.LP
After v1.2, aif was still 31% slower than skani at 1 thread on the B. cereus set and 20% slower on E. coli.
Per-stage counters on B. cereus showed probing at 50% of the time, the fine-binned heterogeneity re-run at 29%,
and refinement at 10%. v1.3 changes:
.IP \(bu 0.2i
\f[B]HMM recursions.\f[R] Each forward/backward step was latency-bound by about 70 cycles: a horizontal dot
product over the identity states, then a division, both feeding the next step. Two changes break this chain.
The forward and backward recursions are now interleaved in one loop. Each step's normaliser is also computed
ahead of time from the previous step's two scalars and sums over older messages, for example
<F[b\(mi1], e[b]> = bz <e[b\(mi1], e[b]> + cz <F[b\(mi2] e[b\(mi1], e[b]>. This removes the reductions from
the chain. The same applies to per-probe refinement, where <e[r\(mi1], e[r]> is a table lookup.
Only floating-point rounding changes. 0.3% of output rows differ, all in the sixth decimal (10\u\s-2\(mi6\s0\d).
.IP \(bu 0.2i
\f[B]Emission cache.\f[R] An 8-probe bin's emission vector depends only on its four hit/miss counts (495
combinations), so these vectors are precomputed per sketch (exact).
.IP \(bu 0.2i
\f[B]GFNI reverse complement.\f[R] Complementing the bases and reversing them within each byte is one
\f[CR]vgf2p8affineqb\f[R] with a constant bit matrix, replacing 13 vector operations, in selection, probing and
sketching (bit-identical).
.IP \(bu 0.2i
\f[B]Vectorised packing and bin counting.\f[R] ASCII-to-2-bit packing does 64 bases per step with two
multiply-adds and \f[CR]vpmovdb\f[R]. Bin counts come from tier/hit bit masks and popcounts (both exact).
Calibration membership tests are batched and prefetched, and calibration samples are written branch-free.
Sketches are byte-identical.
.IP \(bu 0.2i
\f[B]Tried and dropped.\f[R] We tried passing canonical k-mers or tier bits from selection to probing; the
compressed 64-bit stores cost more than the recomputation they saved (net +3%). A 0.01 ANI grid (16 states) was 9%
faster but raised simulated AF error from 0.79 to 0.90 pp. Refining only at posterior crossings raised AF error.
.LP
Accuracy is unchanged. Every per-set metric of Tables 4\(en9 moved by less than 3\(mu10\u\s-2\(mi7\s0\d, and the
simulated main and fresh benchmarks reproduce v1.1 to the printed digit. Step by step, aif's single-thread time on
the B. cereus set, relative to skani, went from 1.55\(mu (v1.1) through 1.16\(mu (HMM changes), 1.09\(mu
(refinement, GFNI) and 1.01\(mu (packing) to 0.95\(mu (bin counting).
.H2 "Timing"
.LP
Table 14 gives medians of 7 interleaved, core-pinned rounds per command, for v1.4 (with confidence intervals).
Single-threaded, aif (sketch + query) takes 0.74\(en0.99\(mu skani's time. The closest set is B. cereus, with its
Mbp-scale contigs (1.33 vs 1.34 s). v1.3, which computes no intervals, was at 0.73\(en0.95\(mu; the intervals cost
about 4%.
With \f[CR]--het-fine-max 16384\f[R] the ratio is 0.67\(en0.79\(mu.
At 8 threads aif is 2.4\(en3.6\(mu faster than skani on every set.
aif's sketch takes 14\(en41 ms per reference.
.TC 14
Wall time (s) per set: skani \f[CR]dist --qi\f[R] (sketch + query) against aif v1.4 sketch + query, at 1 and 8
threads. "aif fast" adds \f[CR]--het-fine-max 16384\f[R]. Peak RSS is for the query step.
.TCE
.so gen/tab_speed.tbl
.TC 15
Thread scaling of aif query on the B. cereus set (449 Mbp, 10,355 contigs): median seconds (peak RSS, MB), all
three versions interleaved in one run.
.TCE
.so gen/tab_scaling.tbl
.ne 2.6i
.PSPIC -C figs/fig_scaling.eps 3.6i
.FC 5
Thread scaling on the B. cereus set: aif v1.1, v1.2 and v1.3, and skani at 1 and 8 threads.
.FE
.LP
\f[B]Small single-reference workloads.\f[R] On the simulated sets of the earlier report (a few hundred contigs
per newly sketched reference), v1.3 is faster than skani on 20 Mbp of contigs (0.067 vs 0.094 s) and 81 Mbp
(0.19 vs 0.30 s).
For 700 contigs per reference, skani is still slightly faster over four references (0.138 vs 0.155 s), because
building aif's sketch (18\(en40 ms) dominates. aif ties on E. coli and B. subtilis and loses on S. aureus and
P. aeruginosa. For high-GC P. aeruginosa, the adaptive-k search builds the sketch twice (k = 19 just misses the
spurious-rate target). Sketches are built once per reference when they are reused.
.LP
\f[B]Memory.\f[R] Single-threaded, aif uses about as much memory as skani or less (19\(en64 MB against 17\(en74 MB).
At 8 threads its peak RSS is higher (106\(en473 MB against 21\(en89 MB). That figure includes up to 128 MB of
memory-mapped input pages, which are shared page cache. The rest is per-thread HMM state for Mbp-scale contigs,
which holds forward, backward and emission values per 8-probe bin. Recomputing emissions in the backward pass,
or checkpointing, would reduce it.
.\" =====================================================================
.H1 "Calibrated confidence intervals (v1.4)"
.LP
skani, fastANI and alignment-based tools report point estimates only. aif's estimators have known sampling
distributions, so every estimate can carry an interval. Output now has six more columns:
\f[CR]ani_lo ani_hi af_lo af_hi ani_se af_se\f[R].
.H2 "Statistical intervals"
.IP \(bu 0.2i
\f[B]ANI.\f[R] The variance of the FP-corrected containment is the probe-level term (FracMinHash sampling and
false positives) plus Theorem 1's mutation-process term. The truth users care about is this contig's realized
identity, not the mutation rate that generated it, so the process term is conditioned on the realized number of
substitutions. This removes the part explained by the substitution count, k\u\s-22\s0\d\(*h\u\s-22\s0\dp/((1\(mip)n), which is
almost all of it at high identity, where k-mer survival is nearly determined by the count. The variance is
carried to ANI by the delta method. For heterogeneous contigs, the spread of the heterogeneity-test windows
gives a cluster-robust standard error, used when larger.
.IP \(bu 0.2i
\f[B]AF.\f[R] The aligned length is uncertain only where the per-probe aligned posterior is intermediate. Each
such run is one alignment boundary whose position has the posterior as its survival function; its variance comes
from the first two moments. The variances are summed over boundaries during refinement, at negligible cost.
.LP
Where the model holds, these intervals are calibrated (Table 17). On simulated substitution-only contigs (16
conditions) they cover 96.5% (ANI) and 99.5% (AF) at nominal 95%. They fail where the model does: indels are not
in it, and ANI coverage there falls to 62% because aif's ANI is biased low by about the indel rate.
On real data with MUMmer truth they cover only 75% (ANI) and 82% (AF) at 95%. Real truth differs from the model
through indels, alignment-definition differences and mosaic structure, not through sampling error.
.TC 16
Coverage of 90/95/99% intervals and median 95% half-widths. Statistical: \f[CR]--ci-stat\f[R]. Calibrated:
default, fitted on the real tuning set only. Real truth is MUMmer; simulated truth is the construction values.
.TCE
.so gen/tab_ci.tbl
.TC 17
Statistical-interval coverage at nominal 95% on the simulated main benchmark, by condition family.
.TCE
.so gen/tab_ci_sim.tbl
.H2 "Calibration against alignment truth"
.LP
To give intervals that hold against alignment truth, aif uses stratified (Mondrian) split-conformal prediction:
.IP \(bu 0.2i
The score is |estimate \(mi truth| / \(*s, where \(*s is the statistical SE with a small floor (0.05 pp for ANI;
1 pp + 1/L for AF).
.IP \(bu 0.2i
Contigs are stratified by estimated identity (\(>= 99, 97, 93% or lower), length (< 1 kb, 1\(en3 kb, longer) and
AF class (fully or partially aligned), giving 24 strata.
.IP \(bu 0.2i
On the real tuning set only, each stratum's conformal quantile at 90, 95 and 99% is the multiplier on \(*s.
Strata with fewer than 50 contigs fall back to coarser pools.
.LP
The table is built into aif (\f[CR]src/ci_table.inc\f[R], from \f[CR]bench_real/calibration.py\f[R]); other
truth definitions can supply their own table (\f[CR]--ci-table\f[R]).
On held-out real contigs the calibrated intervals cover 89.7/94.4/98.6% (ANI) and 89.0/94.2/99.0% (AF) at
90/95/99% (Table 16, Figure 6). Coverage holds within strata too: 88\(en98% at 95% (Table 18). The intervals are
informative where the estimates are good. At \(>= 99% identity on contigs \(>= 3 kb, the median half-widths are
0.44 pp for ANI and 0.54 pp for AF. Where estimates are poor, the intervals say so: below 93% identity they widen
to several pp of ANI, and short, divergent, partially aligned contigs get AF intervals spanning [0, 1].
On simulated data the calibrated intervals over-cover (> 99%), because alignment-definition discrepancies are
absent there; \f[CR]--ci-stat\f[R] is the right choice when truth follows the model.
The intervals cost about 4% in run time (Section 7.4).
.TC 18
Calibrated 95% intervals on held-out real contigs, by stratum (strata with \(>= 100 contigs).
.TCE
.so gen/tab_ci_strata.tbl
.ne 3.0i
.PSPIC -C figs/fig_calibration.eps 6.4i
.FC 6
Empirical against nominal coverage of ANI and AF intervals. Calibrated intervals on held-out real contigs lie on
the diagonal; statistical intervals are calibrated on simulated contigs (ANI) but too narrow against alignment truth.
.FE
.\" =====================================================================
.H1 "Toward a principled AF model"
.LP
aif's main residual weakness on real data is AF over-estimation on partially aligned, divergent contigs: +7 pp
bias, and in the consistency analysis three set \(mu aligned-fraction cells where skani is closer.
aif-RT (Section 6) narrows it by raising the aligned/unaligned switching rate, but costs simulated AF and
genome-level ANI. We looked for a principled correction, chose all settings on the tuning sets only
(Table 19, Figure 7), and evaluated once on held-out data (Table 20).
.IP (1) 0.3i
\f[B]Model averaging over switching rates\f[R], weighting a "stable" and a "mosaic" HMM per contig by
marginal likelihood, with the mixing weight estimated by empirical Bayes. The fitted weight goes to 1 even on
simulated, uniformly divergent contigs. The marginal likelihood prefers the mosaic model where no mosaic exists,
so the emission model is misspecified and averaging cannot fix it.
.IP (2) 0.3i
\f[B]Identity-coupled prior.\f[R] Jumps between identity levels are uniform and about 33\(mu cheaper than
aligned/unaligned switches. An unaligned island inside a 98%-identity region is therefore cheaply explained as
"aligned at 85%". Restricting the aligned levels to within \(*D of the contig's dominant identity (a
probe-weighted median over aligned-looking windows, robust to islands and conserved genes) fixes this. With
\(*D = 0.08 it lowers real AF error from 2.62 to 2.25 pp and also improves simulated AF and ANI.
However, it raises genome-level ANI error from 0.164 to 0.227 pp. Long divergent tracts, which alignment keeps,
are now forced out of the aligned set, and the remaining identity is over-estimated. A local version, which takes
the identity from \(+-R neighbouring windows, keeps long tracts and halves that cost (\(*D = 0.10, R = 1:
0.174 pp), but its AF gain shrinks accordingly.
.IP (3) 0.3i
\f[B]Overdispersed emissions.\f[R] Overlapping k-mers share bases, so probe hits cluster. A beta-binomial bin
emission with the overdispersion implied by the k-mer overlap structure should tolerate the resulting gaps. It
improves simulated AF at the default switching rate, but worsens real ANI, AF and genome-level ANI, and does not
let higher switching rates leave simulated contigs intact. Hit clustering in this form is not the missing piece.
.LP
On these data contig AF and genome-level ANI pull against each other. The regions that inflate aif's AF
(absorbed islands) and the regions alignment counts at 85\(en90% identity (long divergent tracts) produce
similar sparse-hit patterns near the identifiability floor. We therefore keep the defaults, which give the best
genome-level ANI, and offer the local prior (\f[CR]--ident-window 0.08 --ident-radius 3\f[R]) for per-contig
AF. Held out, it lowers contig AF error by 13% (2.45 to 2.14 pp macro) and genome AF error from 2.39 to 1.87 pp,
leaves simulated accuracy unchanged, and raises genome-level ANI error to 0.251 pp, close to fastANI's 0.247 and
above skani genome mode's 0.237. A model with positional anchors on the reference, which could check colinearity,
is the clean way out.
.TC 19
AF-model settings on the tuning sets (MAE in pp; real: macro over 5 sets; genome level: aggregated aif output vs
ANIm; simulated: macro over conditions).
.TCE
.so gen/tab_afgrid.tbl
.ne 3.0i
.PSPIC -C figs/fig_af_tradeoff.eps 6.4i
.FC 7
Tuning-set trade-offs. Left: real contig AF error against genome-level ANI error; no setting improves both over the
default. Right: real contig AF error against simulated AF error. W: per-contig identity window (\(*D = 0.06\(en0.10);
L\(*D/R: local identity prior; od: overdispersed emissions.
.FE
.TC 20
Held-out evaluation of the default, aif-RT and the local identity prior.
.TCE
.so gen/tab_afopt.tbl
.\" =====================================================================
.H1 "Database mode (v1.5)"
.LP
Earlier sections give every query a single reference. In practice contigs are often searched against thousands of
references, as with skani's \f[CR]search\f[R] command. Version 1.5 adds \f[CR]aif db build\f[R],
\f[CR]aif db search\f[R] and \f[CR]aif db info\f[R]. This section describes the design, checks that it is exact,
and compares it with skani on a 2,005-genome database.
.H2 "Design"
.IP \(bu 0.2i
\f[B]One k, one tier layout.\f[R] Tier hash ranges depend only on the budget and the tier specification, not on
the genome. With one k for all references, a contig's selected k-mers, their tiers and their Bloom hashes are
therefore the same for every reference. A contig is prepared once: packing, selection, tier, block hash and
bit-position hash. Each candidate reference then costs only its filter tests and the HMM.
The per-reference filters are exactly those of \f[CR]aif sketch -k K\f[R]. An exhaustive search
(\f[CR]--no-screen -n 0\f[R]) is byte-identical to pairwise \f[CR]aif query\f[R] runs against \f[CR]-k 20\f[R]
sketches, for all 33,464 real contigs against their own and the cross-species references, at 1 and 8 threads.
.IP \(bu 0.2i
\f[B]Screening index.\f[R] Candidates come from a bit-sliced Bloom filter, as in COBS (Bingmann et al., 2019).
k-mers whose selection hash lies in the lowest 1/64 of the hash range (a subset of the probes) are inserted with
3 hash functions into one Bloom column per reference, at 8 bits per screening k-mer. References are sorted by size
into groups of at most 512, with at most a 1.5\(mu size range, and the group's columns share one length. One
64-byte row then answers a hash function for 512 references. A reference's column fill gives its exact screening
false-positive rate f.
A reference is a candidate when its shared-key count c over the contig's n keys has z \(>= 2 against
Binomial(n, f), with c \(>= 2.
Candidates are ranked by FP-corrected screening containment, (c/n \(mi f)/(1 \(mi f), and the best 16 are run
through the full model.
Ranking by z instead favours references with a low f on long contigs and loses 4\(en5 points of recall (Table 25).
The output is the top n hits per contig, ordered by ANI\(muAF (or by ANI or AF).
.IP \(bu 0.2i
\f[B]Lazy, memory-mapped.\f[R] The database file is memory-mapped, so a reference's filters are paged in on first
use. Per-reference estimators are built on first use and share their distance tables.
.H2 "Cost of one shared k"
.LP
Pairwise aif chooses k per reference (17\(en20 on these references) to cap the spurious-match rate.
A database needs one k; 20 suits GC-rich genomes up to about 10 Mbp.
On the pairwise benchmark, k = 20 raises the pooled ANI error from 0.50 to 0.52 pp and the AF error from 2.63 to
2.80 pp (Table 21). Both remain far below skani's.
.TC 21
Pairwise accuracy with adaptive and fixed k (real sets, pooled within-species; MAE in pp).
.TCE
.so gen/tab_db_k.tbl
.H2 "Benchmark database"
.LP
The database holds 2,005 complete genomes (8.62 Gbp; \f[CR]bench_real/db_select.py\f[R], seed 2026):
.IP \(bu 0.2i
the five pairwise references;
.IP \(bu 0.2i
100 random RefSeq complete genomes each of E. coli, S. aureus, P. aeruginosa and the B. cereus group;
.IP \(bu 0.2i
1,600 random RefSeq reference genomes of other bacteria.
.LP
Every assembly used as a query, in the evaluation or tuning sets, is excluded. Nested databases of 405, 805, 1,205
and 2,005 genomes keep the references and near neighbours and add background genomes, for the scaling runs.
The queries are the contigs of the five within-species sets; 26,084 of them are homologous to their pairwise
reference (true AF \(>= 0.2).
.LP
skani 0.4.0 runs in two database modes with default parameters: \f[CR]skani search --qi\f[R] against a
\f[CR]skani sketch\f[R] database, and \f[CR]skani dist --qi\f[R] against the list of all 2,005 FASTA files.
\f[CR]skani search\f[R] screens with an inverted marker index (c = 1000) and drops many short contigs that
\f[CR]skani dist\f[R] reports. \f[CR]skani dist\f[R], however, re-sketches every reference on each run.
.LP
\f[B]Build and size (Table 22).\f[R] The aif database is 1.035 bits per reference base: 0.880 for the filters and
0.155 for the screening index. skani's is 0.896. skani sketches 3\(en4\(mu faster: 3.0 s against 10.9 s for 2,005
genomes at 8 threads, both from the page cache. For an equal-size comparison we also built a size-matched aif
database, with filter budget 0.74, giving 0.895 bits per base (965 MB, the same as skani's).
.TC 22
Database construction at 8 threads. Size is the database file(s); bits/base is per reference base.
* Size-matched aif database (filter budget 0.74).
.TCE
.so gen/tab_db_build.tbl
.H2 "Accuracy of the top hit"
.LP
For every contig we take the top hit of each tool. The aif top hit is ordered by ANI\(muAF; skani's is its
highest-ANI hit. We also take the top three hits of each tool and the exhaustive aif best hit.
Every such (contig, reference) pair is aligned with MUMmer, exactly as in Section 2.2. That gives 145,591 pairs
over 580 references.
Each tool's top-hit ANI and AF are scored against the truth of that pair. Assignment quality is the regret:
the best true ANI\(muAF among all evaluated references minus that of the tool's top hit. We also report the
fraction of contigs whose top hit is the best (within 0.1 pp), and an ANI-only regret: the highest true ANI among
references aligned over at least half the contig, minus the true ANI of the top hit.
.LP
Results are in Tables 23 and 24. aif's database search reports 99.9% of homologous contigs, against 77.0%
for \f[CR]skani search\f[R] and 94.5% for \f[CR]skani dist\f[R].
Its top-hit ANI error is 0.32 pp, against 0.51 and 0.69 pp. On the contigs both tools report, it is 0.26 vs
0.51 pp (search) and 0.30 vs 0.69 pp (dist).
AF error is 1.30 pp against 3.03 and 4.15 pp; counting unreported contigs as AF 0, it is 1.36 against 23.6 and
8.7 pp.
aif's top hit is the best reference for 90.0% of contigs (skani 77.2% and 73.9%). Its regret is 0.20 pp
of ANI\(muAF, against 1.94 and 2.35 pp. On the ANI-only criterion the regret is 0.10 pp, against 0.15 and 0.24 pp.
The size-matched database is almost as accurate: 0.33 pp ANI and 1.38 pp AF.
The calibrated intervals keep their coverage in database mode, 93.8% (ANI) and 94.1% (AF) at 95%, although the
references differ from the pairwise ones they were calibrated on.
.TC 23
Top-hit accuracy against the 2,005-genome database, pooled over the 26,084 homologous contigs of the five sets
(MAE and regret in pp; AF MAE*: unreported contigs count as AF 0; best: fraction whose top hit is the best
evaluated reference; CI: 95% calibrated intervals).
.TCE
.so gen/tab_db_acc.tbl
.TC 24
Per-set comparison with skani search and skani dist. Paired ANI MAE: contigs that aif and skani dist both
report; regret in pp of ANI\(muAF.
.TCE
.so gen/tab_db_sets.tbl
.H2 "Screening"
.LP
Screening recall is the fraction of contigs whose exhaustive best hit (by estimated ANI\(muAF) is also the
screened top hit, or scores within 0.001 of it. It is 96.2% at the default setting (Table 25).
The misses cost little: exhaustive search lowers the pooled regret only from 0.20 to 0.18 pp, and leaves ANI and
AF error unchanged (Table 23). Most misses are in the B. cereus group, where many references are similarly
close mosaics. There, one reference has a higher identity over part of the contig and another a lower identity
over all of it. Containment and ANI\(muAF rank these two differently.
Table 25 compares screening densities. Doubling the screening fraction (0.31 bits/base) raises recall by
1.4 points. A second stage that ranks 32 screened references by the containment of all probes before the HMM
adds little.
.TC 25
Screening configurations on the 2,005-genome database (5 sets, 8 threads; recall against the exhaustive search;
"any hit": contigs with an exhaustive hit that get any hit; time: search time of the five sets). The default is
1/64, 8 bits, 3 hash functions, containment ranking, 16 evaluated.
.TCE
.so gen/tab_db_screen.tbl
.H2 "Speed and scaling"
.LP
aif searches the 33,464 contigs against the 2,005 genomes in 5.2 s at 8 threads and in 34.9 s on one thread
(Table 26, Figure 8). \f[CR]skani search\f[R] takes 37.9 and 264.5 s (7.3\(en7.6\(mu slower).
With \f[CR]--keep-refs\f[R], which keeps passing reference sketches in memory, it takes 15.9 and 83.8 s
(2.4\(en3.1\(mu slower).
\f[CR]skani dist\f[R], which also sketches the references on every run, takes 118.5 s at 8 threads.
From 405 to 2,005 genomes (5\(mu the references), aif's time grows 1.3\(mu. Screening is linear in the
number of reference groups, and noise candidates are drawn from a larger pool.
skani search's time grows 1.1\(mu over the same range.
.LP
Memory is aif's weak point: 2.0 GB at 8 threads against 1.3 GB (skani search) and 1.9 GB
(\f[CR]--keep-refs\f[R]). The peak counts the touched pages of the memory-mapped database (1.1 GB) and the
per-reference estimators. These take about 0.35 MB each and are built for nearly every reference, because
weak noise candidates eventually reach most of them. A cheaper rejection of noise candidates before the HMM,
and estimators shared between references with equal filter error rates, would lower both time and memory.
.TC 26
Search time (s) and peak memory (GB, 8 threads) for the 33,464 contigs of the five sets against the nested
databases (warm page cache; skani dist sketches the references in every run).
.TCE
.so gen/tab_db_time.tbl
.ne 3.0i
.PSPIC -C figs/fig_db_scaling.eps 6.4i
.FC 8
Database search time (left; solid, filled: 8 threads; dashed, open: 1 thread, measured for the smallest and
largest database) and peak memory (right, 8 threads) against database size.
.FE
.\" =====================================================================
.H1 "Discussion"
.H2 "Why the simulation results carry over"
.LP
The advantages measured on simulations come from three properties of aif's design, and none depends on
simulated data:
.IP \(bu 0.2i
dense probing per stored bit;
.IP \(bu 0.2i
query-side segmentation for AF;
.IP \(bu 0.2i
calibrated false-positive noise.
.LP
Real data add mosaic structure, repeats and uneven contig lengths. aif's AF errors grow on these contigs, but far
less than skani's, which are dominated by short and partial contigs. The largest real-data gains are on short
contigs. These are exactly the regime the project targeted: MAG and draft-assembly fragments of 0.5\(en2 kb, where
skani's seed chains are too sparse to report anything.
.H2 "Limitations"
.IP \(bu 0.2i
\f[B]Divergent contigs (< 90% ANI).\f[R] aif over-estimates ANI by about 2.6 pp against MUMmer and sometimes
over-calls AF. Positional anchors (theory report, designs D3/D5) would allow colinearity checks and catch
dispersed-repeat matches. The corrections tested in Section 9 trade contig AF against genome-level ANI.
.IP \(bu 0.2i
\f[B]Calibration scope.\f[R] The calibrated intervals are calibrated to MUMmer truth on these genome sets.
Exchangeability with new data (other clades, assemblers, or truth definitions) is an assumption; the stratified
table can be refitted (\f[CR]--ci-table\f[R]). Statistical intervals do not account for indels.
.IP \(bu 0.2i
\f[B]Truth.\f[R] The truth is MUMmer (blastn-audited). Alignment-based AF is itself uncertain at 85\(en90%
identity. At 85\(en90% the tools and aligners disagree by several points.
.IP \(bu 0.2i
\f[B]Database scale.\f[R] Database mode is tested on 2,005 genomes, not GTDB scale (tens of thousands). The
screening index is linear in the number of references per screening key. An inverted index or a reference
clustering step would be needed beyond about 10\u5\d references. The database uses one k (20), which costs a little
accuracy against adaptive k. Genome-to-genome screening, skani's other use case, is not tested.
.IP \(bu 0.2i
\f[B]Genome level.\f[R] skani \f[CR]--slow\f[R] gives better genome-level AF; aif's per-contig AF over-estimate
carries into the aggregate.
.IP \(bu 0.2i
\f[B]Hardware.\f[R] aif's speed relies on AVX-512 (with VBMI2 and GFNI used where present); the scalar build is
correct but slower. Timings come from one Xeon Gold 5416S host, interleaved to control frequency drift.
.IP \(bu 0.2i
\f[B]Fresh sketches for few contigs.\f[R] When a few hundred contigs meet a newly sketched reference, building
the sketch dominates, and skani remains about 10% faster end to end.
.H2 "Next steps"
.IP 1. 0.3i
Positional anchors for colinearity, then a dispersion test on low-identity calls to cut AF over-calls.
.IP 2. 0.3i
An identity-aware switching prior, so that aif-RT's AF gain comes without its cost on divergent genomes.
.IP 3. 0.3i
Lower per-thread HMM memory, a cheaper adaptive-k search (one build per reference), and pruning of identity
levels for the heterogeneity re-run.
.IP 4. 0.3i
Indel-aware statistical intervals, and conformal tables for other truth definitions (BLAST-based AF, as used for
viral clustering).
.IP 5. 0.3i
Database mode at GTDB scale: an inverted screening index or reference clusters, and per-k preparation so that
references keep their adaptive k.
.\" =====================================================================
.H1 "Reproducibility"
.CS
make aif                                     # aif v1.5 (bin/aif); make prof: stage counters
python3 -I bench_real/select_datasets.py     # genome lists (skani-test metadata)
bench_real/download.sh                       # NCBI datasets + skani-test MAGs
python3 -I bench_real/prepare.py             # contigs >= 500 bp, per-set query files
python3 -I bench_real/truth.py               # MUMmer truth (MUMMER=<bin dir>)
python3 -I bench_real/run_tools.py           # skani presets and aif configurations
python3 -I bench_real/evaluate_real.py --json data/real/results/eval_main.json
python3 -I bench_real/genome_level.py        # skani / fastANI genome mode vs ANIm
python3 -I bench_real/blast_audit.py         # blastn audit
python3 -I bench_real/time_real.py 7         # timing; scaling.py v1.1=<bin> v1.3=bin/aif
python3 -I bench_real/consistency.py aif     # per-stratum aif vs skani errors
python3 -I bench_real/calibration.py fit     # CI table from the tuning set
python3 -I bench_real/calibration.py evaluate  # CI coverage
python3 -I bench_real/af_model_grid.py       # AF-model settings (tuning sets)
python3 -I bench_real/afopt_heldout.py       # AF options, simulated held out
python3 -I bench_real/db_select.py; bench_real/db_download.sh   # 2,005-genome database (NCBI datasets)
python3 -I bench_real/db_bench.py build      # then: search, truth, time, evaluate (database mode)
python3 -I bench_real/db_screen_eval.py <dirs>   # screening recall of search variants vs exhaustive
python3 -I bench_real/select_tune.py         # disjoint real tuning set, then with
  REAL_DIR=data/real_tune REAL_SETS=bcereus,ecoli,mag2328,paer,saur: download, prepare,
  truth, run_tools.py aif, tune_real.py "--rho-u=1e-5,2e-5,3e-5 ..."
/usr/bin/python3 scripts/real_report.py; report_real/build.sh   # this report
.CE
.LP
Genome lists, per-contig truth, predictions, evaluations and timings are committed under
\f[CR]\%results_real/\f[R]; downloaded assemblies and alignments are regenerated by the scripts.
