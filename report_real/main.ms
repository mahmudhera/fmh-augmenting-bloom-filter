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
Real-genome report, aif v1.2 \(em November 22, 2026
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
No aif parameter was tuned on these data. At equal sketch size (0.88 vs 0.90 bits per reference base) aif v1.1/v1.2
is better than skani 0.4.0 on every set, for both ANI and AF.
Pooled over the 33,464 within-species contigs, aif reports 99.2% of homologous contigs (skani 90.1%).
Its ANI mean absolute error is 0.50 pp (skani 0.68 pp), or 0.39 vs 0.68 pp on the 23,488 contigs that every
configuration reports. Its AF error is 2.63 pp (skani 9.35 pp).
aif also beats skani's 3.8\(mu larger \f[CR]--slow\f[R] preset on every set. At half its budget, aif still beats
skani's default in pooled ANI (0.56 vs 0.68 pp) and in AF on every set.
On 500\(en1000 bp contigs aif reports 95% of homologous contigs, against 23% for skani, with an AF error of 4.3 vs 39 pp.
Aggregating aif's contig estimates also gives the most accurate genome-level ANI against ANIm (0.19 pp vs 0.24 pp for
skani and 0.25 pp for fastANI).
One setting change, chosen on a disjoint real tuning set (aif-RT), lowers the held-out AF error to 1.85 pp.
Real assemblies exposed a speed problem: v1.1's multithreaded pipeline did not scale on Mbp-scale contigs.
Version 1.2 fixes this with cache-sized probe windows and a streaming thread pipeline. It is 7\(mu faster than v1.1
at 8 threads, gives byte-identical output, and is 2.6\(mu faster than skani at 8 threads over the five sets.
Single-threaded, it is within 20% of skani over the five sets, and an optional switch closes that gap.
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
\f[B]Speed.\f[R] At 8 threads, aif v1.2 processes the five sets in 1.05 s against skani's 2.68 s.
Single-threaded it takes 3.88 vs 3.27 s, or 3.20 s with \f[CR]--het-fine-max 16384\f[R].
Version 1.1 did not scale with threads on these inputs (Section 7).
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
aif runs with its v1.1 defaults. The same aif output is produced by v1.2 (Section 7); all accuracy numbers are
v1.2 runs.
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
.H1 "Speed and memory on real data (v1.2)"
.H2 "What real assemblies exposed"
.LP
The simulated contigs were at most 20 kb. Draft assemblies contain contigs up to 5 Mbp, and that broke two
assumptions of v1.1.
.IP (1) 0.3i
\f[B]Thread pipeline.\f[R] v1.1 read 4 Mbp batches, handed out 16 contigs per grab, and waited at a barrier after
each batch. A batch with a few Mbp-scale contigs kept one thread busy while the others idled. Every thread also
grew its scratch buffers to the largest contig. On the 449 Mbp B. cereus set, v1.1 took 3.6 s at 1 thread and
3.6 s at 8 threads, using 1.95 GB at 8 threads (Table 15).
.IP (2) 0.3i
\f[B]Cache behaviour.\f[R] Probing a 5 Mbp contig in one pass streams megabytes of packed sequence and positions
through L2. Per-stage cycle counters (\f[CR]make prof\f[R]) showed probing at 6.9 cycles per base on \(>= 500 kb
contigs, against 4.2 on 20\(en500 kb contigs.
.H2 "Changes"
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
.H2 "Timing"
.LP
Table 14 gives medians of 7 interleaved, core-pinned rounds per command.
At 8 threads aif is 1.9\(en3.1\(mu faster than skani on every set; over the five sets it takes 1.05 s against 2.68 s.
Single-threaded, aif is within 7% of skani on the three sets with shorter contigs. It is 20\(en31% slower on the
B. cereus and E. coli sets, where Mbp-scale contigs trigger the costly fine-binned heterogeneity pass.
With \f[CR]--het-fine-max 16384\f[R], aif is faster than skani on three sets and within 3% on the other two.
aif's sketch takes 17\(en45 ms per reference.
.TC 14
Wall time (s) per set: skani \f[CR]dist --qi\f[R] (sketch + query) against aif sketch + query, at 1 and 8 threads.
"aif fast" adds \f[CR]--het-fine-max 16384\f[R]. Peak RSS is for the query step.
.TCE
.so gen/tab_speed.tbl
.TC 15
Thread scaling of aif query on the B. cereus set (449 Mbp, 10,355 contigs): median seconds (peak RSS, MB).
.TCE
.so gen/tab_scaling.tbl
.ne 2.6i
.PSPIC -C figs/fig_scaling.eps 3.6i
.FC 5
Thread scaling on the B. cereus set: aif v1.1, v1.2, and skani at 1 and 8 threads.
.FE
.LP
\f[B]Memory.\f[R] Single-threaded, aif uses as much memory as skani or less (19\(en69 MB against 17\(en74 MB).
At 8 threads its peak RSS is higher (115\(en534 MB against 22\(en88 MB). That figure includes up to 128 MB of
memory-mapped input pages, which are shared page cache. The rest is per-thread HMM state for Mbp-scale contigs,
which holds forward, backward and emission values per 8-probe bin. Recomputing emissions in the backward pass,
or checkpointing, would reduce it.
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
dispersed-repeat matches.
.IP \(bu 0.2i
\f[B]Truth.\f[R] The truth is MUMmer (blastn-audited). Alignment-based AF is itself uncertain at 85\(en90%
identity. At 85\(en90% the tools and aligners disagree by several points.
.IP \(bu 0.2i
\f[B]Pairwise only.\f[R] Every query has a single reference. Multi-reference database search, and skani's
genome-to-genome screening use cases, are not tested.
.IP \(bu 0.2i
\f[B]Genome level.\f[R] skani \f[CR]--slow\f[R] gives better genome-level AF; aif's per-contig AF over-estimate
carries into the aggregate.
.IP \(bu 0.2i
\f[B]Hardware.\f[R] aif's speed relies on AVX-512. Timings come from one Xeon Gold 5416S host, interleaved to
control frequency drift.
.H2 "Next steps"
.IP 1. 0.3i
Positional anchors for colinearity, then a dispersion test on low-identity calls to cut AF over-calls.
.IP 2. 0.3i
An identity-aware switching prior, so that aif-RT's AF gain comes without its cost on divergent genomes.
.IP 3. 0.3i
Lower per-thread HMM memory, and pruning of identity levels for the heterogeneity re-run.
.IP 4. 0.3i
A database mode (many references per query) evaluated on GTDB-scale data, as in the skani paper.
.\" =====================================================================
.H1 "Reproducibility"
.CS
make aif                                     # aif v1.2 (bin/aif); make prof: stage counters
python3 -I bench_real/select_datasets.py     # genome lists (skani-test metadata)
bench_real/download.sh                       # NCBI datasets + skani-test MAGs
python3 -I bench_real/prepare.py             # contigs >= 500 bp, per-set query files
python3 -I bench_real/truth.py               # MUMmer truth (MUMMER=<bin dir>)
python3 -I bench_real/run_tools.py           # skani presets and aif configurations
python3 -I bench_real/evaluate_real.py --json data/real/results/eval_main.json
python3 -I bench_real/genome_level.py        # skani / fastANI genome mode vs ANIm
python3 -I bench_real/blast_audit.py         # blastn audit
python3 -I bench_real/time_real.py 7         # timing; scaling.py <aif_v1.1> bin/aif
python3 -I bench_real/select_tune.py         # disjoint real tuning set, then with
  REAL_DIR=data/real_tune REAL_SETS=bcereus,ecoli,mag2328,paer,saur: download, prepare,
  truth, run_tools.py aif, tune_real.py "--rho-u=1e-5,2e-5,3e-5 ..."
/usr/bin/python3 scripts/real_report.py; report_real/build.sh   # this report
.CE
.LP
Genome lists, per-contig truth, predictions, evaluations and timings are committed under
\f[CR]\%results_real/\f[R]; downloaded assemblies and alignments are regenerated by the scripts.
