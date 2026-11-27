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
Real-genome report, aif v1.6, with a v1.7 addendum (Section 17) \(em November 27, 2026
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
We test aif on real genome pairs: 240 assemblies with 35,362 contigs of at least 500 bp (1.13 Gbp) from the
skani paper's data sets, namely B. cereus group and E. coli drafts, MAGs of species bin 2328, P. aeruginosa and
S. aureus drafts, and negative controls. The truth is MUMmer alignment, and no aif parameter was tuned on these
data. At equal sketch size (0.88 vs 0.90 bits per reference base) aif beats skani 0.4.0 on every set. Pooled
over 33,464 within-species contigs, it reports 99.2% of homologous contigs (skani 90.1%), with 0.50 vs 0.68 pp
ANI error and 2.6 vs 9.4 pp AF error. On 0.5\(en1 kb contigs it reports 95% against 23%. Aggregated, its contig
estimates give the most accurate genome-level ANI (0.19 pp; skani 0.24, fastANI 0.25).
Speed work (v1.2\(en1.3) makes aif faster than skani single-threaded. Calibrated per-contig intervals (v1.4)
cover 94% of held-out contigs at 95%.
A database mode (v1.5\(en1.6) shares one k across references, so each contig is hashed once and each candidate
costs only filter tests and the HMM. Against 2,005 genomes, aif's top-hit errors are 0.31 pp ANI and 1.1 pp AF,
against 0.51/3.0 pp for skani search, and it is 8\(mu faster.
Against all 199,923 GTDB R232 species representatives (655 Gbp), aif builds an 82.5 GB database in 6.7 min
(skani: 82.7 GB, 10.8 min) and classifies the contigs of 294 held-out genomes in 11.5 s (skani: 252 s). Its top
hit has the right species for 78% of contigs and the right genus for 89% (skani 62% and 69%). It is more
accurate at every contig length, and flags 91% of contaminant contigs in simulated bins (skani 68%).
Across sketch sizes from 0.11 to 3.6 bits per base aif is on the Pareto front. Exact FracMinHash containment at
32\(mu the size is less accurate than aif's smallest sketch, as the bits-per-sample theorem predicts.
aif also reports identity maps along contigs. These are more accurate than fastANI's fragment identities
(0.55 vs 0.78 pp at 3 kb), with calibrated window intervals, and they find unaligned islands at 97% precision and
81% recall. Database defaults were frozen and re-validated on fresh queries and unseen species, with
unchanged conclusions. Kraken2, given the same references and database size, names the species about as often
as aif (77\(en78% vs 78\(en79%) and the genus more often (90\(en92% vs 89\(en90%), but it reports no ANI or AF,
searches 6\(mu slower and builds 35\(mu slower. On 400 simulated bins aif estimates contamination more
accurately than GUNC and CheckM2 (0.67 vs 1.96 and 6.27 pp error). It detects contamination as well as GUNC's
own call, and runs 2,100\(mu faster.
We also report where aif errs, notably AF over-estimates on divergent partial contigs.
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
\f[B]Database mode.\f[R] Against 2,005 genomes, aif's top hit has 0.31 pp ANI and 1.1 pp AF error, against
0.51 / 3.0 pp for skani search and 0.69 / 4.2 pp for skani dist. aif covers 99.9% of homologous contigs
(77% and 95%), and is 8\(mu faster than skani search, 3.6\(mu faster than its \f[CR]--keep-refs\f[R] mode. It uses
somewhat more memory (1.8 vs 1.3 GB) and a 16% larger database; a size-matched database is as accurate (Section 10).
.IP (9) 0.3i
\f[B]GTDB scale.\f[R] Against all 199,923 GTDB R232 representatives, aif classifies contigs of held-out genomes
to the right species 78% of the time (skani 62%), genus 89% (69%), and 22\(mu faster, with an equal-size
database. With an equal species-calling rule (ANI \(>= 95%) it calls 77% of contigs correctly at 86% precision,
against 61% at 88% for skani (Section 11). Contaminant contigs in simulated bins are flagged at genus level
for 91% (skani 68%, Section 12).
.IP (10) 0.3i
\f[B]Size\(enaccuracy front.\f[R] Across sketch sizes from 0.11 to 3.6 bits per base, aif is more accurate than
skani and exact FracMinHash at every size, and faster up to 0.88 bits per base. Alignment (minimap2) is more
accurate still, with a 31\(mu larger index and 210\(mu the run time (Section 13).
.IP (11) 0.3i
\f[B]Identity maps.\f[R] aif's segmentation gives windowed identity along contigs (0.55 pp error at 3 kb; fastANI
0.78 pp), calibrated window intervals (94.4\(en94.7% coverage), unaligned islands (97% precision, 81% recall) and
breakpoints (median 171 bp) (Section 14).
.IP (12) 0.3i
\f[B]Held-out confirmation.\f[R] With the database defaults frozen, fresh query sets reproduce every conclusion: on
80 new drafts against 2,005 genomes (ANI 0.23 vs 0.41 pp for skani search), and on 298 genomes of 300 unseen GTDB
species (species accuracy 78.6% vs 66.1%) (Section 15).
.IP (13) 0.3i
\f[B]Other tools.\f[R] With the same 199,923 references in an equal-size database, Kraken2 gets the species right
for 76.8% and 78.1% of contigs (aif 78.3 and 78.6%) and the genus for 90.1% and 92.1% (aif 88.5 and 90.1%). At
equal sensitivity aif's species calls are 2.4\(en2.7 points more precise. Kraken2 searches 6.4\(mu slower, takes
35\(mu longer to build, and reports no ANI or AF. On 400 simulated bins aif's contamination error is
0.55\(en0.84 pp; GUNC's is 0.69\(en3.86 pp and CheckM2's 2.9\(en10.7 pp. With a 3% limit aif detects
contamination as often as GUNC's own call (95\(en100% of bins, 6% of controls), and searches all 400 bins in
11 s, against 23.5 min for CheckM2 and 6.6 h for GUNC (Section 16).
.IP (14) 0.3i
\f[B]v1.7.\f[R] A benchmark on gut contigs (\f[CR]aif_uhgg_contamination_report.pdf\f[R], Section 7) showed that
aif over-called AF on real contigs: unaligned stretches inside an alignment were absorbed by the lowest identity
levels. v1.7 counts AF only near the local identity, which lowers the AF error on every real set here (macro
2.45 to 1.97 pp; fresh drafts 2.53 to 1.97) and raises it on simulated contigs of heterogeneous identity. ANI and
the hits are unchanged, and the intervals were refitted (Section 17). The rest of this report is v1.6.
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
.H1 "Database mode (v1.5\(en1.6)"
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
Long contigs use their 128 keys of smallest selection hash. This is a consistent subsample, so containment stays
unbiased, and it changed the top hit of 0.9% of contigs.
.IP \(bu 0.2i
\f[B]Pre-HMM filter (v1.6).\f[R] A two-state Viterbi score is computed for candidates without a clear global
hit excess (z < 4). It uses the HMM's emission log-ratios at eight identity levels and its switching
penalties, and costs a fraction of the HMM. A candidate scoring below 5 nats skips the HMM.
On the 2,005-genome benchmark no reported hit among 370,757 scored below 9.1, and 99% of the candidates that
report nothing were skipped. Output is identical; the main saving is estimator memory at GTDB scale.
.IP \(bu 0.2i
\f[B]Hit order (v1.6).\f[R] Hits are ordered by the model's log-likelihood ratio, i.e. how well each
reference explains the contig's probe hits (alternatives: ANI\(muAF, ANI, AF). Ordering by ANI\(muAF
lets aif's AF over-estimate on divergent hits outrank a closer reference. On GTDB, a 111 kb contig at
99.9% ANI and 0.90 AF to its own species lost to another species at 93.1% ANI and 0.99 AF. The likelihood
ratio lowers the top-hit ANI and AF errors and the ANI-only regret below (v1.5, ordered by ANI\(muAF:
0.32 pp, 1.30 pp, 0.10 pp), at a small cost in the ANI\(muAF regret (0.20 \(-> 0.28 pp).
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
For every contig we take the top hit of each tool. The aif top hit is the first by likelihood ratio; skani's is its
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
Its top-hit ANI error is 0.31 pp, against 0.51 and 0.69 pp. On the contigs both tools report, it is 0.25 vs
0.51 pp (search) and 0.29 vs 0.69 pp (dist).
AF error is 1.12 pp against 3.03 and 4.15 pp; counting unreported contigs as AF 0, it is 1.19 against 23.6 and
8.7 pp.
aif's top hit is the best reference for 88.9% of contigs (skani 78.3% and 74.9%). Its regret is 0.28 pp
of ANI\(muAF, against 1.81 and 2.24 pp. On the ANI-only criterion the regret is 0.06 pp, against 0.15 and 0.25 pp.
The size-matched database is almost as accurate: 0.32 pp ANI and 1.22 pp AF.
The calibrated intervals keep their coverage in database mode, 94.0% (ANI) and 94.5% (AF) at 95%, although the
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
Screening recall is the fraction of contigs whose exhaustive best hit is also the screened top hit, or scores
within 0.001 of it. It is 97.5% with the v1.6 defaults. The grid in Table 25 was measured in v1.5, with
ANI\(muAF ordering and no key cap. The misses cost little: exhaustive search lowers the pooled regret only from
0.28 to 0.27 pp, and leaves ANI and AF error unchanged (Table 23). Most misses are in the B. cereus group, where many references are similarly
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
With the v1.6 defaults, aif searches the 33,464 contigs against the 2,005 genomes in 4.4 s at 8 threads and in
32.9 s on one thread (Table 26, Figure 8). \f[CR]skani search\f[R] takes 37.9 and 264.5 s (8.0\(en8.6\(mu
slower). With \f[CR]--keep-refs\f[R], which keeps passing reference sketches in memory, it takes 15.9 and
83.8 s (2.5\(en3.6\(mu slower).
\f[CR]skani dist\f[R], which also sketches the references on every run, takes 118.5 s at 8 threads.
From 405 to 2,005 genomes (5\(mu the references), aif's time grows 1.1\(mu, as does skani search's.
In v1.5 aif's time grew 1.3\(mu (5.2 s at 2,005 genomes), because noise candidates went through the HMM.
.LP
aif uses more memory: 1.8 GB at 8 threads against 1.3 GB (skani search) and 1.9 GB
(\f[CR]--keep-refs\f[R]). The peak counts the touched pages of the memory-mapped database (1.1 GB) and the
per-reference estimators, about 0.35 MB each. The pre-HMM filter builds estimators only for references that
pass it: 748 rather than 2,002 of the 2,005 (v1.5: 2.0 GB peak).
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
.H1 "GTDB-scale database: contig classification"
.LP
The database is all 199,923 species representatives of GTDB release R232 (April 2026): 189,801
bacterial and 10,122 archaeal genomes, 655 Gbp. They come from the release tarball (192 GB, MD5-verified),
and aif reads the gzip files directly (v1.6).
.LP
\f[B]Queries.\f[R] The queries are 294 GTDB non-representative genomes from 294 distinct species, 15 of them
archaeal (\f[CR]bench_real/gtdb_select.py\f[R], seed 2611). Each is an NCBI draft assembly (contig or
scaffold level, not metagenome-derived) with CheckM2 completeness \(>= 95% and contamination \(<= 2%. Six
further genomes, ENA submissions, are not served by NCBI datasets and were dropped. The 23,561 contigs of at
least 500 bp (1.32 Gbp) are classified one by one.
Each query's own species representative is in the database, and the query genome itself is not.
.LP
\f[B]Truth.\f[R] A contig's top hit is correct at species (genus) level when the hit's GTDB species (genus)
is the query genome's. Top-hit ANI and AF are scored against MUMmer, as in Section 10, for every
(contig, reference) pair among each tool's top three hits and the contig's own species representative.
Contigs that align to their own species representative over at least 20% of their length form the
"homologous" subset. The others are plasmids, mobile elements and accessory regions, whose best match may
legitimately lie in another species.
.LP
\f[B]Tools.\f[R]
.IP \(bu 0.2i
aif at budget 0.88: default screening, hits ordered by the model's likelihood ratio;
.IP \(bu 0.2i
a smaller aif database, at filter budget 0.22;
.IP \(bu 0.2i
skani 0.4.0: \f[CR]skani sketch\f[R], then \f[CR]skani search --qi\f[R], with hits ordered by ANI.
.LP
Builds use 32 threads; searches use 16 threads, with the database in the page cache.
.LP
\f[B]Build and size (Table 27).\f[R] The aif database is 82.5 GB, 1.008 bits per base: 72.0 GB of filters
and 10.5 GB of screening index. skani's is 82.7 GB, 1.011 bits per base. With the genomes in the page cache,
aif builds it in 6.7 min with 1.9 GB of memory, against 10.8 min and 15.6 GB for skani. From a cold disk
aif's first build took 30 min, mostly reading 181 GB of gzip files. The 0.22-budget aif database is 28.5 GB
(0.348 bits per base), and its screening index is 37% of it.
.LP
\f[B]Search.\f[R] With the databases in the page cache, aif searches the 23,561 contigs (1.32 Gbp) against the 199,923
genomes in 11.5 s at 16 threads, and in 7.4 s with the 0.22 database. skani search takes 252 s, 22\(en34\(mu
longer. The first runs, with parts of the databases on disk, took 13\(en77 s for aif and 852 s for skani.
Peak resident memory is 65 GB for aif and 78 GB for skani. For aif it is almost all touched pages of the
memory-mapped database file, which the page cache shares. The 19,726 estimators built (one per reference that
passed the pre-HMM filter) add about 7 GB.
.LP
\f[B]Classification (Table 28, Figure 9).\f[R]
.IP \(bu 0.2i
\f[B]Coverage and accuracy.\f[R] aif finds a hit for 93.5% of contigs, skani for 71.6%. aif's top hit has the
query's species for 78.3% of all contigs and its genus for 88.5%; skani's for 61.8% and 69.0%.
.IP \(bu 0.2i
\f[B]Homologous contigs.\f[R] On contigs that align to their own species representative, aif's top hit is the
right species for 91.4% and skani's for 72.3%.
.IP \(bu 0.2i
\f[B]By length.\f[R] aif is more accurate at every contig length: 45.7 vs 2.7% at 0.5\(en1 kb,
62.6 vs 18.5% at 1\(en2 kb, 83.1 vs 76.1% at 5\(en20 kb, and 92.5 vs 90.5% at \(>= 20 kb.
.IP \(bu 0.2i
\f[B]Equal calling rule.\f[R] Species calls are made only when the top-hit ANI is \(>= 95%, GTDB's species
boundary. Under this rule aif calls 76.6% of contigs correctly at 86.4% precision, and skani 60.9% at 87.5%.
.IP \(bu 0.2i
\f[B]Archaea.\f[R] On the 1,015 archaeal contigs, aif's species accuracy is 87.4% against 69.1%.
.IP \(bu 0.2i
\f[B]ANI and AF.\f[R] Against MUMmer, aif's top-hit ANI error is 0.49 pp and skani's 0.60 pp, over a larger set
of contigs for aif. AF error is 2.58 vs 3.34 pp.
.IP \(bu 0.2i
\f[B]Smaller database.\f[R] The 0.22-budget database, a third of skani's size, still classifies 77.8% of
contigs to species (88.1% to genus), with 0.58 pp ANI error.
.LP
\f[B]Interval coverage under shift.\f[R] The calibrated intervals cover 91.5% (ANI) and 91.9% (AF) of top-hit
pairs here, against 95% nominal. They were calibrated on five bacterial species, and these queries span the
whole tree, so the gap measures distribution shift. Refitting on a taxonomically broad tuning set
(\f[CR]--ci-table\f[R]) is the remedy.
.LP
\f[B]Where species assignment fails.\f[R] For 300 contigs whose own representative was not among the reported
hits, aif's candidate list was examined. In 245 the representative was not a candidate at all. These were
mostly short contigs whose screening containment is tied (0.87\(en1.0) across many genomes, such as conserved
and mobile elements, or contigs about equally distant from many references. Neither a larger candidate list
nor taxonomy alone can resolve them, and they account for most of the gap between all contigs (78.3%) and
homologous ones (91.4%).
.TC 27
Building and searching the GTDB R232 database (199,923 genomes, 655 Gbp; build: 32 threads; search: the
23,561 query contigs, 16 threads). Peak memory of the search includes the memory-mapped database pages.
.TCE
.so gen/tab_gtdb_build.tbl
.TC 28
Contig classification against GTDB R232 (correct: the top hit's GTDB species or genus is the query's;
precision: among contigs with a hit; homol.: contigs aligning to their own species representative; ANI \(>= 95% call:
species assigned only when the top-hit ANI is \(>= 95% (sensitivity over all contigs, precision over calls);
ANI/AF: top hit against MUMmer).
.TCE
.so gen/tab_gtdb_acc.tbl
.ne 2.8i
.PSPIC -C figs/fig_gtdb.eps 6.4i
.FC 9
Species-level classification accuracy (left) and top-hit ANI error (right) against contig length,
GTDB R232.
.FE
.\" =====================================================================
.H1 "Contamination screening of bins"
.LP
Every contig was classified on its own, so contaminated bins can be assembled after the fact
(\f[CR]bench_real/gtdb_contam.py\f[R]). For each of the 294 query genomes (the host) and each contamination
level (5, 10 and 20% of host bases), contigs of one other random query genome are added until the level is
reached, skipping contigs that would overshoot it by more than 25%: 861 bins (hosts whose contaminant has no
small-enough contig are skipped). The bin's taxon is the base-weighted majority of its contigs' top-hit taxa. A contig is
flagged when its top hit has another species, or with the genus rule another genus. Contigs without a hit
are left unflagged, or flagged under the strict rule.
.LP
The genus rule gives the best precision for both tools (Table 29). Pooled over the three levels, aif flags
98.8% of contaminant bases at 95.1% precision, and skani 97.6% at 93.1%. By contig count aif flags 90.8% of
contaminant contigs at 75.3% precision, and skani 67.8% at 81.8%. On 0.5\(en2 kb contaminant contigs the recall is 75.5% against 13.3%:
skani reports nothing for most short contigs, so they cannot be flagged. Their bases are few, which is why the
base-level recalls are close. Flagging every contig without a hit (strict rule) recovers almost all contaminant
contigs for both tools, at a large cost in precision. The species rule flags host contigs whose best match lies
in another species (conserved or mobile elements), and its precision is lower for both tools (aif 76%,
skani 71% by bases).
.TC 29
Contaminant detection in 861 simulated bins: precision and recall of contaminant bases, by contamination
level (taxon: species or genus rule; strict: contigs without a hit are flagged too).
.TCE
.so gen/tab_contam.tbl
.\" =====================================================================
.H1 "Size, accuracy and time: a Pareto comparison"
.LP
The earlier comparisons use one sketch size per tool. Here each method is swept across sizes on the five
real sets (pairwise; pooled within-species contigs; time is sketch plus query on one thread):
.IP \(bu 0.2i
\f[B]aif\f[R] at budgets 0.11\(en1.76 bits per reference base;
.IP \(bu 0.2i
\f[B]skani\f[R] at compression factors c = 500 to 30, measured at 0.26\(en3.27 bits per base;
.IP \(bu 0.2i
\f[B]exact FracMinHash\f[R]: containment ANI with sourmash 4.8 (k = 21), at scaled values 580 to 18.
Hashes are stored at 64 bits each, so this spans 0.11\(en3.56 bits per base; sourmash estimates no AF;
.IP \(bu 0.2i
\f[B]alignment\f[R] with minimap2 (asm20 preset, identity from the gap-compressed divergence) as the
accuracy ceiling.
.LP
Results are in Table 30 and Figure 10.
.IP (1) 0.3i
\f[B]aif is on the front at every size.\f[R] At 0.22 bits/base its ANI error (0.62 pp) is below
skani's default at 0.84 bits/base (0.68 pp). At 0.44 bits/base (0.56 pp) it beats skani
\f[CR]--slow\f[R] at 3.27 bits/base (0.61 pp), with 7.4\(mu less space. At 0.44 bits/base its AF error
(3.6 pp) is below skani's at every size tested (4.3 pp at best).
.IP (2) 0.3i
\f[B]Filtered samples beat exact samples.\f[R] Exact FracMinHash at 0.88 bits/base has an ANI error of 0.77 pp,
worse than aif at 0.11 bits/base (0.67 pp). Even at 3.56 bits/base, 32\(mu aif's smallest sketch, it reaches
only 0.71 pp. This is the bits-per-sample principle (Theorem 3) on real data: at equal space, many
samples stored with 0.5\(en10% false positives carry more information than fewer samples stored exactly.
Containment alone cannot separate AF from ANI (Proposition 3), and sourmash reports no AF.
.IP (3) 0.3i
\f[B]Speed.\f[R] Up to 0.88 bits/base aif is also the fastest: 1.35\(en2.91 s for the five sets, against
2.42\(en3.22 s for skani at comparable sizes. sourmash (Python driving Rust) takes 21\(en184 s.
.IP (4) 0.3i
\f[B]Alignment.\f[R] minimap2 reaches 0.23 pp ANI and 1.78 pp AF error. Its index is 27 bits/base (31\(mu aif's
default sketch) and its run time 622 s (210\(mu aif's). At 1.76 bits/base aif's AF error (1.73 pp) matches
minimap2's, and its ANI error is twice as high.
.ne 4.0i
.TC 30
Pairwise accuracy against sketch (or index) size and time on the real sets (pooled within-species
contigs; MAE in pp; paired: contigs that every configuration reports; time: sketch + query of the five
sets, 1 thread, median of 3).
.TCE
.so gen/tab_pareto.tbl
.ne 3.0i
.PSPIC -C figs/fig_pareto.eps 6.6i
.FC 10
Error against sketch size (left, middle) and against time (right). Sketch size is bits per reference base:
aif budget, skani sketch files, 64 bits per sourmash hash, minimap2 index. sourmash has no AF estimate.
.FE
.\" =====================================================================
.H1 "Identity maps along contigs"
.LP
aif's HMM segments each contig, and v1.6 exposes the segmentation. With \f[CR]--segments FILE\f[R],
\f[CR]aif query\f[R] and \f[CR]aif db search\f[R] write aligned and unaligned segments that tile each contig.
Each aligned run is split into windows of about 1 kb (\f[CR]--seg-window\f[R]), and each window gets an
FP-corrected identity with a 95% interval. No other sketch tool reports positions. The closest comparator
is fastANI's \f[CR]--visualize\f[R], which lists the identity of each mapped 3 kb fragment.
.LP
\f[B]Truth.\f[R] Each query genome is aligned to its reference with MUMmer, as in Section 2.2. A window's true
identity is 1 \(mi (substitutions + indel bases in the window) / (aligned bases in the window), from
\f[CR]show-snps\f[R]. Alignment-block identity would average over whole blocks, often tens of kb long.
A base is unaligned when no alignment covers it. Windows are scored when at least 90% of their bases align.
.LP
\f[B]Intervals.\f[R] Statistical intervals (Theorem 1 with the delta method) cover 89% of held-out 1 kb
windows and 78% of 3 kb windows. Split-conformal multipliers, fitted on the real tuning set by window size and
estimated identity, give 94.7% and 94.4% at 95%. They are built into aif; \f[CR]--seg-ci-stat\f[R] reports
the statistical intervals.
.LP
\f[B]Results\f[R] (Tables 31\(en32, Figure 11; 224 genomes, 833,183 windows of 1 kb):
.IP \(bu 0.2i
\f[B]Window identity.\f[R] At fastANI's 3 kb resolution, aif's error is 0.55 pp against fastANI's 0.78 pp, and its
correlation with the truth is 0.97 against 0.86. aif is better above 90% identity, by up to 7\(mu at
\(>= 99.5% (0.035 vs 0.235 pp). fastANI is better below 90% (2.2 vs 2.7 pp), where aif's detection floor
biases identity upward. At 1 kb aif's error is 0.63 pp.
.IP \(bu 0.2i
\f[B]Islands.\f[R] Of 139 Mbp of unaligned sequence on reported contigs, aif's unaligned segments recover 81% at
97% precision. fastANI's unmapped fragments recover 51% at 56% precision.
.IP \(bu 0.2i
\f[B]Breakpoints.\f[R] Across 67,853 aligned/unaligned boundaries, the median distance to aif's nearest
segment boundary is 171 bp, and 92% are within 500 bp. The median is 59\(en70 bp on E. coli and
P. aeruginosa, against fastANI's 3 kb fragments.
.LP
Recall is lowest on the B. cereus group (73%). There, 85\(en90% regions sit between "aligned at low identity"
and "unaligned", the same ambiguity as in the AF model (Section 9).
.TC 31
Window identity against MUMmer (windows with true AF \(>= 0.9; MAE in pp, r: Pearson correlation; 95% CI
coverage of aif's statistical and calibrated window intervals on held-out data).
.TCE
.so gen/tab_idmap.tbl
.TC 32
Unaligned islands (per base, on contigs each tool reports) and aligned/unaligned breakpoints (aif).
.TCE
.so gen/tab_islands.tbl
.ne 3.0i
.PSPIC -C figs/fig_idmap.eps 6.6i
.FC 11
Identity map of a 162 kb P. aeruginosa draft contig (the P. paraeruginosa-like strain, 93.6% genome ANIm)
against PAO1: aif's 1 kb windows with calibrated 95% intervals, MUMmer identity of the same windows,
fastANI's 3 kb fragments, and aif's unaligned segments (shaded).
.FE
.\" =====================================================================
.H1 "Validation of the frozen defaults on fresh queries"
.LP
Several database-mode defaults were chosen while looking at the results reported above:
.IP \(bu 0.2i
the screening configuration (1/64 \(mu 8 bits \(mu 3 hashes, containment ranking, 16 candidates), on the
2,005-genome benchmark (Table 25);
.IP \(bu 0.2i
the pre-HMM filter cutoff (\f[CR]--vmin 5\f[R]) and the 128-key cap, on the same benchmark;
.IP \(bu 0.2i
hit ordering by likelihood ratio, after inspecting the GTDB results (Section 10.1).
.LP
Other choices were made properly on separate data: the shared k (planned before the database benchmark), the
interval calibration (real tuning set), the identity-map multipliers (tuning set) and \f[CR]--af-min-ident\f[R]
(viral report, tuning set). The defaults that were not were therefore frozen
(git tag \f[CR]v1.6-freeze\f[R], binary MD5 \f[CR]be055d32\f[R]). Two fresh query sets were then drawn and
scored once:
.IP \(bu 0.2i
\f[B]Database benchmark, fresh queries\f[R] (\f[CR]bench_real/select_fresh.py\f[R], seed 4242). There are 80 new
draft assemblies: 24 B. cereus group (D2, including the only 2 unused drafts at 85\(en90% ANIm), 24 E. coli (D3),
and 16 each of P. aeruginosa and S. aureus. All are disjoint from the evaluation, tuning and database genomes.
They give 9,617 contigs, 7,088 of them homologous to their pairwise reference, searched against the same
2,005-genome database.
.IP \(bu 0.2i
\f[B]GTDB, query set 2\f[R] (\f[CR]GTDB_SET=2\f[R], seed 4611). There are 298 non-representative genomes (15 archaeal)
from 300 species not used in set 1, giving 19,555 contigs (1.33 Gbp), searched against the same GTDB R232
database.
.LP
\f[B]Results\f[R] (Tables 33\(en34). Every conclusion holds.
.IP \(bu 0.2i
\f[B]Database benchmark.\f[R] On the fresh queries aif reports 99.9% of homologous contigs (skani search 80.2%,
dist 92.4%). Its top-hit errors are 0.23 pp ANI and 0.77 pp AF, against 0.41 / 2.26 pp (search) and
0.52 / 2.93 pp (dist). On paired contigs the ANI error is 0.17 vs 0.41 pp (search) and 0.19 vs 0.52 pp (dist).
Its ANI regret is 0.08 pp (skani 0.20 and 0.27), and its top hit is the best reference for 84.9% of contigs
(66.3% and 62.4%). Screening recall is 98.0%, and interval coverage 95.0% (ANI) and 95.8% (AF).
.IP \(bu 0.2i
\f[B]GTDB.\f[R] On set 2 aif's top hit has the right species for 78.6% of contigs (set 1: 78.3%) and the right
genus for 90.1% (88.5%); skani's for 66.1% and 73.3%. aif is again more accurate at every contig length:
40.7 vs 3.0% at 0.5\(en1 kb, and 93.6 vs 91.9% at \(>= 20 kb. Under the ANI \(>= 95% calling rule it
calls 76.6% correctly at 85.8% precision (skani 65.0% at 88.8%). Top-hit errors are 0.49 pp ANI and 2.32 pp AF
(skani 0.60, 3.20). Its margin over skani is somewhat smaller than on set 1 (12.5 vs 16.5 points of species
accuracy). Searches took 12.8 s (aif) and 212 s (skani).
.IP \(bu 0.2i
\f[B]Contamination.\f[R] In the bins simulated from set 2 (genus rule), aif flags 93.0% of contaminant contigs
(skani 70.3%), and 98.2% of contaminant bases at 97.0% precision (skani 97.4% at 95.2%).
.LP
The tuned defaults were therefore not fitted to the benchmark queries: on unseen genomes and unseen species
they give the same results, or slightly better.
.TC 33
2,005-genome database: original queries (Section 10, 26,084 homologous contigs) and fresh queries (7,088),
frozen defaults (MAE and regret in pp; best: fraction whose top hit is the best evaluated reference).
.TCE
.so gen/tab_fresh_db.tbl
.TC 34
GTDB R232: query set 1 (Section 11) and the fresh set 2 (disjoint species), frozen defaults.
.TCE
.so gen/tab_fresh_gtdb.tbl
.\" =====================================================================
.H1 "Comparison with a taxonomic classifier and bin-quality tools"
.LP
Sections 11 and 12 compare aif only with skani. Two other kinds of tool answer related questions. k-mer
taxonomic classifiers (Kraken2) assign contigs to taxa. Gene-based bin-quality tools (CheckM2, GUNC) estimate
how contaminated a bin is. Neither kind reports ANI or AF, so they are compared on classification and
contamination only.
.H2 "Kraken2 on GTDB R232"
.LP
Kraken2 2.1.3 was given exactly aif's references: the 199,923 GTDB R232 representatives, with GTDB's
taxonomy written as an NCBI-style taxonomy of 247,369 taxa (\f[CR]bench_real/kraken2_gtdb.py\f[R]). Its
defaults were kept (k = 35, minimizer length 31, no masking). The database was capped at aif's size
(\f[CR]--max-db-size\f[R] 82.5 GB); Kraken2 then keeps a hash-chosen subset of minimizers.
Uncapped, Kraken2 estimated its hash table at 990 GB, so the capped database keeps about 8% of the
minimizers (those with the largest hash values).
Kraken2's default, deterministic build parallelises within one sequence at a time. On GTDB's many short contigs
it read 2 MB/s (30 GB of the 618 GB library in 4.2 h, at 37 cores), so the build was rerun with
\f[CR]--fast-build\f[R] and 16 MB read blocks (\f[CR]--block-size\f[R]; the 16 KB default left the threads
waiting on the input lock). The fast build is not bit-for-bit reproducible between runs.
Kraken2 assigns each contig to the lowest common ancestor of its k-mer hits. A contig counts as correct at
species (genus) level when the assigned taxon is the query's species (or lies within the query's genus).
Assignments above species rank are Kraken2's way of abstaining, so its species-rank assignments are set
against aif's ANI \(>= 95% calls (Table 36). Searches used 16 threads with the hash table in the page cache, as
for aif and skani.
.LP
\f[B]Results (Tables 35\(en36, Figure 12).\f[R]
.IP \(bu 0.2i
\f[B]Species and genus.\f[R] At the same database size, aif and Kraken2 name the species about equally often.
aif's top hit has the right species for 78.3% of contigs on set 1 and 78.6% on set 2; Kraken2's assignment for
76.8% and 78.1%. Kraken2 assigns more contigs (97.5 and 98.7%, against aif's 93.5 and 94.6%) and has the right
genus more often (90.1 and 92.1%, against 88.5 and 90.1%). When a contig's k-mers are shared by several species of
one genus, its lowest-common-ancestor rule falls back to the genus; aif's top hit then names one of the species.
skani trails both (61.8 and 66.1% species).
.IP \(bu 0.2i
\f[B]Species calls.\f[R] Kraken2's species-rank assignments recover the right species for 76.8 and 78.1% of
contigs at 84.0 and 83.1% precision. aif's ANI \(>= 95% calls recover 76.6% at 86.4 and 85.8% precision. At
about the same sensitivity, aif's calls are 2.4 and 2.7 points more precise.
.IP \(bu 0.2i
\f[B]Contig length.\f[R] aif's species accuracy is level with or above Kraken2's at every length: 45.7 vs 44.2%
(set 1) and 40.7 vs 39.9% (set 2) at 0.5\(en1 kb, and 92.5 vs 92.0% and 93.6 vs 93.3% at \(>= 20 kb. Both keep
most of their accuracy on short contigs, where skani reports almost nothing.
.IP \(bu 0.2i
\f[B]Cost.\f[R] Kraken2's build took 3.9 h at 48 threads: 1.9 h to estimate the hash table size and 2.0 h to
fill it, with 88.7 GB peak memory. aif's took 6.7 min at 32 threads (30 min from a cold disk) and 1.9 GB.
Kraken2 loads the whole hash table into memory, so each search peaks at 83 GB. It took 73.7 and 82.1 s, 6.4
times aif's 11.5 and 12.8 s.
.LP
Kraken2 reports a taxon and no measure of similarity. aif reports the same labels about as accurately,
together with ANI and AF and their intervals for every hit, from a database it builds 35 times faster. Kraken2
kept 8% of its minimizers to fit aif's size; an uncapped Kraken2 database would be about 12 times larger and
probably more accurate.
.TC 35
Kraken2, aif and skani on GTDB R232 (199,923 representatives): database size, build and search cost. Searches
at 16 threads with the database in the page cache; peak memory includes mapped database pages.
.TCE
.so gen/tab_k2_build.tbl
.TC 36
Classification on GTDB R232, query set 1 / set 2. Species call: aif's ANI \(>= 95% calls; Kraken2's
assignments at species rank (sensitivity: correct calls over all contigs). Last two columns: species correct
among contigs of 0.5\(en1 kb and of \(>= 20 kb.
.TCE
.so gen/tab_k2_acc.tbl
.ne 3.0i
.PSPIC -C figs/fig_k2.eps 6.4i
.FC 12
Species accuracy by contig length on GTDB R232: aif's and skani's top hits and Kraken2's assignments, for
query set 1 and the fresh set 2.
.FE
.H2 "CheckM2 and GUNC on simulated bins"
.LP
The gene-based tools take minutes per bin, so they were run on a fixed subset of the bins of Section 12:
the first 100 bins of each contamination level, and 100 uncontaminated controls (the first 100 host genomes
alone), 400 bins in all (\f[CR]bench_real/contam_tools.py\f[R]). True contamination is contaminant bases
divided by host bases.
.IP \(bu 0.2i
CheckM2 1.1.0 (default models, its DIAMOND database) estimates contamination from redundant marker genes.
.IP \(bu 0.2i
GUNC 1.1.1 (GTDB r214 database) reports the fraction of genes from non-dominant clades at the taxonomic level
of maximal clade separation, and its own pass or fail call. GUNC pins DIAMOND 2.1.24; the installed DIAMOND
2.2.8 was used through GUNC's documented override (\f[CR]GUNC_SKIP_DIAMOND_VERSION_CHECK\f[R]).
.IP \(bu 0.2i
aif and skani estimate contamination as flagged over unflagged bases under the genus rule of Section 12.
.LP
All four tools ran at 24 threads. A bin counts as flagged when its estimate is \(>= 5%, the MIMAG
high-quality limit.
.LP
\f[B]Results (Table 37, Figure 13).\f[R]
.IP \(bu 0.2i
\f[B]Size of the contamination.\f[R] aif's estimates track the truth: the error is 0.55, 0.63 and 0.84 pp at 5,
10 and 20% contamination, and 0.51 pp on controls. skani's errors are 0.90\(en1.35 pp. GUNC's contamination
portion is close at 5% (0.69 pp) but falls short as contamination grows (1.33 and 3.86 pp; a mean of 17.0%
at a true 20.7%). CheckM2 under-estimates at every level, with errors of 2.9, 5.2 and 10.7 pp and a mean of
7.7% against a true 12.1%.
.IP \(bu 0.2i
\f[B]Detection.\f[R] At the 5% limit aif flags 86, 98 and 99% of the bins at the three levels, and 2% of the
controls. CheckM2 flags 23, 61 and 70%, and no controls. GUNC's own call rests on clade separation rather than
on the size of its estimate. It fails 94, 100 and 99% of the contaminated bins and 5% of the controls. A fixed
5% limit is a poor detector for bins whose true contamination is about 5%. With a 3% limit aif flags 95,
100 and 99% of the bins and 6% of the controls, level with GUNC's call.
.IP \(bu 0.2i
\f[B]Cost.\f[R] For the 400 bins at 24 threads, CheckM2 took 23.5 min and GUNC 6.6 h, each with 16.6 GB of
peak memory. GUNC's DIAMOND search shared the machine with the Kraken2 build for most of its run, so its time
is an upper bound.
For comparison, aif searched the contigs of all 400 bins (36,614 contigs, 1.83 Gbp, written as one file) in
11.0 s at 24 threads, with the database in the page cache and 46.6 GB peak memory (mapped database pages).
That is 128 times faster than CheckM2 and about 2,100 times faster than GUNC. skani took 206 s. The top hits
of these timed runs are identical, contig for contig, to those of the full-set searches used for the estimates.
.TC 37
Contamination estimates on 400 simulated bins (100 per level of 5, 10 and 20% contamination, and 100
uncontaminated controls). MAE against the true contamination; flagged: estimate \(>= 5% (rows \(lqflag at
3%\(rq: \(>= 3%); GUNC, own call: pass.GUNC false. Time and peak memory for all 400 bins at 24 threads.
.TCE
.so gen/tab_bintools.tbl
.ne 2.6i
.PSPIC -C figs/fig_bintools.eps 6.4i
.FC 13
Estimated against true contamination for the 400 bins (dashed: equality; grey line: the 5% limit).
.FE
.LP
The gene-based tools do things not tested here. They need no reference genome close to the bin, and CheckM2
also estimates completeness. aif and skani judge contamination by the taxa of the contigs' top hits, so they
depend on the reference database. Here every host species has a representative in GTDB R232; a host without a
relative of its own genus in the database would be split among other genera and its contigs flagged. GUNC's
database (GTDB r214) is older than the R232 database used by aif and skani.
.H2 "Tools not run"
.LP
MMseqs2 taxonomy and CAT/BAT classify contigs by protein search against a reference protein database. For
GTDB R232 that means predicting, storing and indexing about 700 million proteins from the 199,923
representatives, and a translated search for every contig. This was beyond the scope of this comparison.
These tools target a different question: placing contigs of novel lineages at higher ranks, where nucleotide
identity carries little signal. Their comparison with aif on such contigs is left open. CAT was later run on gut
contigs in the UHGG report (Section 7): with the GTDB R232 proteins (a 227 GB DIAMOND database), its genus calls
gave a foreign-contig F1 of 0.41, and the catalogue would take about 100 days.
.\" =====================================================================
.H1 "Addendum: aif v1.7, AF near the local identity"
.LP
The contig benchmark on the UHGG gut catalogue (\f[CR]aif_uhgg_contamination_report.pdf\f[R], Section 7) found aif
v1.6's top-hit AF 5.1 pp above MUMmer's. Base by base, three quarters of the excess were gaps bridged between two
MUMmer alignments, at aif window identities near its lowest levels. Those levels (about 85%) expect a hit rate close
to the filters' false-positive rates, so an unaligned stretch inside an alignment costs almost nothing as
"aligned", while two switches cost about 14 nats. The same mechanism was behind the AF over-calls on divergent
contigs analysed in Section 9.
.LP
v1.7 counts AF over the aligned posterior mass within 0.05 of the local identity: the median FP-corrected
containment of aligned-looking windows of 64 probes, three windows either side. A window that favours its best
aligned level over unaligned by 3 nats, or a neighbour of such a window, sets its floor at its own identity, so
divergent blocks keep their AF. The HMM, the hits and ANI are unchanged; \f[CR]\-\-af-local 0\f[R] reproduces
v1.6. The settings were chosen on tuning sets only: a gut tuning set disjoint from the UHGG benchmark, the real
tuning set of Section 6 and the simulated tuning set.
.LP
Table 38 gives the held-out effect on the data of this report. The AF error falls on all nine real sets, by a
fifth on average (main 2.45 to 1.97 pp, fresh 2.53 to 1.97), and the bias drops from +1.5\(en1.7 to +0.3\(en0.7 pp.
On simulated contigs it rises from 5.30 to 5.83 pp. Uniformly mutated contigs at 90\(en99% identity lose 0.5\(en1 pp,
and contigs of heterogeneous identity lose more: 5.8 pp for alternating 100%/90% blocks (\f[CR]hetblock\f[R]) and
4.4 pp for gamma-distributed rates (\f[CR]hetgamma\f[R]), where a short divergent block next to conserved ones
does not always show enough evidence of its own. The \f[CR]hetblock\f[R] failure was first seen on the held-out
simulated benchmark, and it shaped the evidence rule.
.LP
\f[B]Genome level.\f[R] Aggregated from contigs as in Section 4, genome-level ANI is unchanged (error 0.140 to
0.141 pp, macro over the five main sets; fresh drafts 0.145 both), and the genome-level AF error against ANIm
halves (1.87 to 0.97 pp; fresh 2.10 to 1.17; \f[CR]bench_real/af_local_genome.py\f[R]). The identity prior of
Section 9 changed which regions the HMM aligns and so paid in genome-level ANI; v1.7 leaves the HMM and contig ANI
alone and moves only the AF weights, which the AF gain makes more accurate.
.LP
\f[B]Intervals.\f[R] The conformal tables were refitted for the new AF (Table 39). Pairwise queries keep a table
fitted on the real tuning set, and \f[CR]aif db search\f[R] now has its own, fitted on gut contigs against GTDB R232,
where the pairwise table under-covered (92% of ANI and 87% of AF at a nominal 95%; 95.1% and 94.2% with the new
table). On the real sets the coverage stays at its nominal level and the median AF half-width falls from 2.60 to
1.98 pp, because AF is more accurate. Simulated contigs remain over-covered, as before.
.TC 38
AF error of aif v1.6 and v1.7 on the held-out real and simulated data of this report (every contig; macro over sets or
conditions; ANI is identical in both versions).
.TCE
.so gen/tab_v17_af.tbl
.TC 39
Coverage of the calibrated intervals (ANI/AF, %) at three levels, v1.6 and v1.7 (pairwise query).
.TCE
.so gen/tab_v17_ci.tbl
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
\f[B]Database scale.\f[R] At 200,000 references each screening key reads about 1,200 random cache lines.
This is affordable with at most 128 keys per contig, but the cost grows linearly with the database; an inverted
index or reference clusters would be needed well beyond GTDB. Searches assume the database fits in the page
cache (82.5 GB for GTDB at the default budget). The database uses one k (20), which costs a little accuracy
against adaptive k. Interval coverage drops to 92% on GTDB's diverse taxa (calibration covers five bacterial
species). Genome-to-genome screening, skani's other use case, is not tested.
.IP \(bu 0.2i
\f[B]Genome level.\f[R] skani \f[CR]--slow\f[R] gives better genome-level AF; aif's per-contig AF over-estimate
carries into the aggregate.
.IP \(bu 0.2i
\f[B]Hardware.\f[R] aif's speed relies on AVX-512 (with VBMI2 and GFNI used where present); the scalar build is
correct but slower. Timings come from one Xeon Gold 5416S host, interleaved to control frequency drift.
.IP \(bu 0.2i
\f[B]Fresh sketches for few contigs.\f[R] When a few hundred contigs meet a newly sketched reference, building
the sketch dominates, and skani remains about 10% faster end to end.
.IP \(bu 0.2i
\f[B]Other tools.\f[R] Kraken2 was compared at one size, aif's, where it keeps 8% of its minimizers; protein-based
classifiers (MMseqs2, CAT) were not run. The bin comparison uses simulated bins whose host species are all in the
database, which favours the reference-based estimates of aif and skani over the gene-based tools.
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
Beyond GTDB scale: an inverted screening index or reference clusters, and per-k preparation so that references
keep their adaptive k. A taxonomically broad calibration set for the intervals.
.IP 6. 0.3i
Application studies on the new outputs: identity maps across a pangenome (recombination tracts and islands),
and contamination screening of real MAG catalogues against GTDB.
.IP 7. 0.3i
Novel lineages: contigs and bins from genera absent from the database, compared with protein-based classifiers
and with GUNC and CheckM2, which do not need close references.
.\" =====================================================================
.H1 "Reproducibility"
.CS
make aif                                     # aif v1.6 (bin/aif; zlib); make prof: stage counters
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
bench_real/gtdb_download.sh                  # GTDB R232 representatives (parallel ranges, MD5)
python3 -I bench_real/gtdb_select.py; (NCBI datasets download); python3 -I bench_real/gtdb_prepare.py
python3 -I bench_real/gtdb_bench.py list     # then: build, search, truth, evaluate
python3 -I bench_real/gtdb_contam.py         # contamination screening of simulated bins
python3 -I bench_real/kraken2_gtdb.py prepare   # then: build, classify, evaluate (Kraken2 on GTDB R232)
python3 -I bench_real/contam_tools.py bins     # then: checkm2, gunc, seqtime, evaluate (bin-quality tools)
python3 -I bench_real/run_tools.py aif_b011 aif_b022 aif_b176 skani_c60 skani_c250 skani_c500
<sourmash python> -I bench_real/pareto_sourmash.py; python3 -I bench_real/pareto_minimap2.py
python3 -I bench_real/pareto_time.py         # Pareto timing
python3 -I bench_real/idmap.py run; REAL_DIR=data/real_tune ... run, calibrate; idmap.py evaluate
python3 -I bench_real/select_tune.py         # disjoint real tuning set, then with
  REAL_DIR=data/real_tune REAL_SETS=bcereus,ecoli,mag2328,paer,saur: download, prepare,
  truth, run_tools.py aif, tune_real.py "--rho-u=1e-5,2e-5,3e-5 ..."
/usr/bin/python3 scripts/real_report.py; report_real/build.sh   # this report
.CE
.LP
Genome lists, per-contig truth, predictions, evaluations and timings are committed under
\f[CR]\%results_real/\f[R]; downloaded assemblies and alignments are regenerated by the scripts.
