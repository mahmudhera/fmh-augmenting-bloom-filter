.\" Build: report_viral/build.sh   (tables/figures: /usr/bin/python3 scripts/viral_report.py)
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
aif for Viral Genome Clustering
.ft R
.ps 14
.vs 18p
vOTU-level ANI and Alignment Fraction of Phage Genomes
against BLAST, skani and vclust
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Viral clustering report, aif v1.6 \(em November 23, 2026
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
Viral genomes and contigs are grouped into viral operational taxonomic units (vOTUs) when they share at least
95% ANI over at least 85% of the shorter sequence. The standard truth is all-vs-all BLAST followed by CheckV's
\f[CR]anicalc\f[R]/\f[CR]aniclust\f[R], and the AF threshold decides many pairs. We compared aif's database mode
with skani (default and \f[CR]--small-genomes\f[R]) and the alignment-based vclust. The data are two samples
of 10,000 sequences: complete phage genomes from INPHARED (600 Mbp, 55,751 same-vOTU pairs) and gut
phage contigs from the Gut Phage Database (374 Mbp, 8,311 same-vOTU pairs).
aif's ANI is the most accurate of the sketch-based tools near the vOTU thresholds (0.57\(en0.63 pp against
0.89\(en1.10 pp for skani), and at 95\(en98% identity it beats vclust as well. With its default AF, however, it
over-estimates BLAST's AF on partially aligned pairs by 6\(en12 pp. Megablast does not align the divergent
stretches that aif's HMM counts as homologous.
We therefore added an identity-thresholded AF, the fraction aligned at identity \(>= T, and chose T = 0.88 on a
disjoint tuning set. On the held-out samples it lowers the near-threshold AF error from 6.3\(en7.2 to 2.5\(en3.2 pp.
It then gives the best vOTU agreement among the sketch-based tools: edge F1 0.938 / 0.902 and adjusted Rand index
(ARI) 0.896 / 0.807 (INPHARED / GPD), against 0.917 / 0.857 and 0.890 / 0.750 for skani's small-genome mode.
vclust's alignments remain closest to BLAST (F1 0.970 / 0.948, ARI 0.937 / 0.893), at 7\(en10\(mu aif's run time
and twice its memory. skani is the fastest, at 2\(en4 s per 10,000 sequences against 5\(en8 s for aif.
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
This report tests aif on the clustering of viral genomes into vOTUs, the last application experiment of the
project plan. The real-genome report (\f[CR]aif_real_genomes_report.pdf\f[R]) covers bacterial contigs,
database mode and GTDB scale.
.IP (1) 0.3i
\f[B]ANI near the thresholds.\f[R] On pairs near the vOTU boundaries, aif's ANI error is 0.57 pp (INPHARED) and
0.63 pp (GPD). skani's is 0.99 and 1.10 pp (\f[CR]--small-genomes\f[R]: 0.89, 0.94), and vclust's 0.34 and
0.55 pp (Table 2). At 95\(en98% true identity aif has the lowest error of all four tools on INPHARED
(0.48 pp, vclust 0.80, skani 0.88\(en0.89). At 80\(en90% it has the highest, near its detection floor
(Table 3).
.IP (2) 0.3i
\f[B]AF definition.\f[R] With its default AF, aif over-estimates BLAST's AF by 5\(en12 pp on partially aligned
pairs (Figure 1). The model counts divergent homologous stretches that megablast, with its 28-mer seeds, does not
align. The option \f[CR]--af-min-ident T\f[R] reports the fraction aligned at identity \(>= T. T = 0.88 was
chosen on a disjoint tuning set of 3,000 INPHARED genomes (Table 6). On the held-out data it cuts the
near-threshold AF error from 6.3 to 2.5 pp (INPHARED) and from 7.2 to 3.2 pp (GPD). Mid-range AFs are then
under-estimated, but near the 85% threshold the bias is small.
.IP (3) 0.3i
\f[B]vOTUs.\f[R] With the thresholded AF, aif's vOTU-edge F1 is 0.938 (INPHARED) and 0.902 (GPD), and its
clusters agree with BLAST's at ARI 0.896 and 0.807. That beats skani (F1 0.915 / 0.833, ARI 0.888 / 0.723)
and skani \f[CR]--small-genomes\f[R] (0.917 / 0.857; 0.890 / 0.750). vclust is best (F1 0.970 / 0.948, ARI 0.937 /
0.893) (Table 4).
.IP (4) 0.3i
\f[B]Cost.\f[R] At 16 threads, all-vs-all comparison of 10,000 genomes takes aif 5\(en8 s (database build
included), skani 2\(en4 s and vclust 35\(en74 s. Peak memory is about 2 GB for aif, 0.1\(en0.4 GB for skani and
4\(en4.6 GB for vclust (Table 5, Figure 2). BLAST, the truth, takes 39\(en139 s at 56 threads.
.\" =====================================================================
.H1 "Background"
.LP
The MIUViG standard defines species-level vOTUs by \(>= 95% ANI over \(>= 85% of the shorter sequence (AF).
In practice the truth comes from all-vs-all \f[CR]blastn\f[R] (megablast) and CheckV's \f[CR]anicalc.py\f[R]
and \f[CR]aniclust.py\f[R] scripts. ANI is the alignment-length-weighted identity of the non-overlapping
alignments, AF the merged aligned fraction, and the clustering a greedy centroid one. Because viral contigs are
short and fragmented, the AF threshold decides many pairs. AF is also where aif's design differs most from
skani's (\f[CR]aif_real_genomes_report.pdf\f[R], Sections 3 and 9).
.LP
The tools compared:
.IP \(bu 0.2i
\f[B]aif v1.6\f[R], database mode. \f[CR]aif db build -i\f[R], new in this report, makes every FASTA record a
reference. \f[CR]aif db search\f[R] then runs all sequences against the database, with up to 100 hits per
sequence and 100 candidates. The AF of a pair is taken from the search whose query is the shorter sequence.
.IP \(bu 0.2i
\f[B]skani 0.4.0\f[R] \f[CR]triangle -i --sparse\f[R], with defaults (c = 125), and with
\f[CR]--small-genomes\f[R] (c = 30, 200 markers), which skani recommends for viruses.
.IP \(bu 0.2i
\f[B]vclust 1.3.1\f[R], an alignment-based tool (Kmer-db prefilter, then LZ-ANI alignment): \f[CR]prefilter
--min-ident 0.8\f[R] and \f[CR]align\f[R]. It reports ANI and query/reference coverage for every prefiltered pair.
.LP
Every tool's pair table goes through the same greedy centroid clustering (ANI \(>= 95%, AF \(>= 85%), as in
\f[CR]aniclust.py\f[R]. Differences in vOTUs therefore come only from the pair estimates.
.\" =====================================================================
.H1 "Data and truth"
.IP \(bu 0.2i
\f[B]INPHARED\f[R] (Cook et al., 2021; release of 14 April 2025): all complete phage genomes from GenBank,
excluding RefSeq duplicates. There are 28,665 genomes; we drew 10,000 at random (seed 31). They are not
dereplicated, so the sample has pairs at every identity.
.IP \(bu 0.2i
\f[B]Gut Phage Database\f[R] (GPD; Camarillo-Guerrero et al., 2021): 142,809 viral contigs of at least 10 kb
from gut metagenomes, 71% of them genome fragments by CheckV. We drew 10,000 at random. GPD was dereplicated
at 95% identity, so the sample is rich in near-threshold rejections and partial overlaps. Under the
CheckV-style definition it still contains 8,311 same-vOTU pairs.
.IP \(bu 0.2i
\f[B]Tuning set\f[R]: 3,000 further INPHARED genomes, disjoint from the test sample (seed 32), used only to
choose \f[CR]--af-min-ident\f[R].
.LP
\f[B]Truth.\f[R] Each sample is searched all-vs-all with \f[CR]blastn -task megablast -evalue 1e-5\f[R]
(BLAST+ 2.12, 56 threads, \f[CR]-mt_mode 1\f[R]). ANI and AF follow CheckV's \f[CR]anicalc.py\f[R]:
.IP \(bu 0.2i
alignments are sorted by bitscore, and those overlapping kept ones on the query by more than 10% are dropped;
.IP \(bu 0.2i
ANI is the alignment-length-weighted identity of the kept alignments;
.IP \(bu 0.2i
AF is the merged aligned fraction, of the shorter sequence.
.LP
A pair is a vOTU edge when ANI \(>= 95% and AF \(>= 85%. "Near the thresholds" means true ANI 93\(en97%
(with AF \(>= 0.5), or ANI \(>= 93% with AF 75\(en95%. All metrics compare each tool's estimate for a pair
with BLAST's.
.TC 1
Data sets (BLAST pairs: pairs with any megablast alignment; pairs \(>= 80% ANI: with AF \(>= 0.1; vOTU edges:
ANI \(>= 95% and AF \(>= 85%).
.TCE
.so gen/tab_data.tbl
.\" =====================================================================
.H1 "Results"
.H2 "Pair estimates"
.LP
Tables 2 and 3 compare pair estimates with BLAST. Over all pairs with \(>= 80% true ANI:
.IP \(bu 0.2i
\f[B]Coverage.\f[R] vclust, which aligns every prefiltered pair, reports almost all (99\(en100%). aif reports
90% (INPHARED) and 77% (GPD), skani \f[CR]--small-genomes\f[R] 95% and 79%, and skani's default 88% and 62%.
.IP \(bu 0.2i
\f[B]80\(en90% identity.\f[R] Most pairs any tool misses are here. aif's coverage falls to 50\(en69% and its
error rises to 3.4\(en4.3 pp, its weakest regime, as on bacterial contigs.
.IP \(bu 0.2i
\f[B]Above 90% identity.\f[R] aif's ANI error is lower than both skani modes in every stratum. It is also lower
than vclust at 95\(en98% on INPHARED (0.48 vs 0.80 pp), and on GPD above 95% (0.47\(en0.70 vs 1.6\(en2.0 pp).
vclust is best at \(>= 98% on INPHARED (0.06 vs 0.09 pp) and at 80\(en90%.
.LP
Near the vOTU thresholds, aif's ANI error is 0.57 pp (INPHARED) and 0.63 pp (GPD), against 0.89\(en1.10 pp
for skani and 0.34\(en0.55 pp for vclust.
.TC 2
Pair estimates against BLAST/anicalc (MAE in pp; reported: fraction of the truth pairs the tool reports).
.TCE
.so gen/tab_pairs.tbl
.TC 3
ANI error by true identity (pairs with AF \(>= 0.1; rep.: fraction reported; MAE in pp).
.TCE
.so gen/tab_strata.tbl
.H2 "The AF definition"
.LP
AF is where the tools differ most (Figure 1). Relative to BLAST's megablast-based AF:
.IP \(bu 0.2i
\f[B]aif's default AF\f[R] is too high by 5\(en12 pp for pairs that align over 30\(en85% of the shorter
sequence, and by 2\(en6 pp near the 85% threshold;
.IP \(bu 0.2i
\f[B]vclust\f[R] also over-estimates partial overlaps (+8 pp at AF 0.1\(en0.3), but is unbiased near the threshold;
.IP \(bu 0.2i
\f[B]skani\f[R] under-estimates, strongly so in its default mode (\(mi10 to \(mi23 pp at AF 0.5\(en0.85).
.LP
The cause is definitional. aif's HMM counts a stretch as aligned when its probes look homologous, down to its
\(ap84% detection floor. Megablast needs 28 bp exact seeds and leaves 85\(en90% stretches unaligned, so they do
not count toward BLAST's AF. On bacterial contigs, MUMmer's sensitive alignment is the truth and this choice
matters less; for vOTUs the convention is BLAST.
.LP
\f[B]Identity-thresholded AF.\f[R] aif v1.6 adds \f[CR]--af-min-ident T\f[R]. The HMM already holds each bin's
posterior over identity levels, so the AF is scaled by the share of the aligned posterior mass at levels
\(>= T. This is the fraction of the sequence aligned at identity \(>= T, and the AF interval is scaled alike.
On the tuning set (Table 6), T = 0.86 changes little, because aif's lowest level is already \(ap0.84. T = 0.88
halves the near-threshold AF error (5.5 \(-> 2.2 pp) and raises edge F1 from 0.933 to 0.942. T = 0.90 has an
F1 only 0.002 higher but twice the near-threshold error, and T = 0.92 loses recall. The rule, fixed before the
test sets were scored, was highest F1, with ties within 0.005 broken by the near-threshold AF error; it chose
T = 0.88. On the held-out samples the near-threshold AF error drops from 6.3 to 2.5 pp (INPHARED) and from 7.2
to 3.2 pp (GPD), close to skani \f[CR]--small-genomes\f[R] (2.4, 3.0 pp). Mid-range AFs are now
under-estimated (Figure 1), which matters little for the 85% threshold.
.ne 3.0i
.PSPIC -C figs/fig_afbias.eps 6.4i
.FC 1
Mean AF error against BLAST's AF of the shorter sequence, for pairs with \(>= 80% ANI (vertical line: the 85% vOTU
threshold).
.FE
.H2 "vOTU edges and clusters"
.LP
Table 4 scores the vOTU decisions.
.IP \(bu 0.2i
\f[B]aif, default AF.\f[R] It finds 97\(en99% of true edges, but its AF excess adds false ones (precision
0.88 and 0.76).
.IP \(bu 0.2i
\f[B]aif, thresholded AF.\f[R] Precision rises to 0.91 and 0.86 at a small cost in recall. F1 rises to
0.938 (INPHARED) and 0.902 (GPD), and the cluster agreement to ARI 0.896 and 0.807. That beats both skani modes on
both data sets.
.IP \(bu 0.2i
\f[B]skani.\f[R] The default misses many edges (recall 0.90 and 0.81), because its sketches are too sparse for
short phages. \f[CR]--small-genomes\f[R] recovers recall but loses precision.
.IP \(bu 0.2i
\f[B]vclust.\f[R] Its alignments are closest to BLAST's (F1 0.970 / 0.948, ARI 0.937 / 0.893), and it reproduces
the most BLAST vOTUs exactly.
.LP
The dereplicated GPD sample is the harder one for every tool: its same-vOTU pairs are mostly partial overlaps
near the AF threshold.
.TC 4
vOTU edges (ANI \(>= 95%, AF \(>= 85%) and greedy-centroid vOTUs against BLAST/anicalc (identical clusters: the
fraction of multi-member BLAST vOTUs reproduced exactly).
.TCE
.so gen/tab_votu.tbl
.H2 "Speed and memory"
.LP
Each tool ran all-vs-all on 10,000 sequences at 16 threads (Table 5, Figure 2).
.IP \(bu 0.2i
\f[B]skani\f[R] is the fastest: 1.7\(en4.0 s with 0.1\(en0.4 GB.
.IP \(bu 0.2i
\f[B]aif\f[R] takes 4.8\(en8.0 s, of which about 2 s builds the database (48 and 77 MB, 1.03 bits per base),
with about 2 GB of memory. Its search builds one estimator per reference, here all 10,000.
.IP \(bu 0.2i
\f[B]vclust\f[R] takes 35\(en74 s and 4.0\(en4.6 GB, mostly in its prefilter.
.IP \(bu 0.2i
\f[B]BLAST\f[R], the truth, takes 39\(en139 s at 56 threads.
.LP
aif's run time is dominated by the HMM on homologous pairs. Phages of the same genus give many candidates per
query, and up to 100 are evaluated.
.TC 5
Resources for all-vs-all comparison of 10,000 sequences (16 threads; aif: database build, then search; skani:
\f[CR]triangle\f[R]; vclust: prefilter, then align).
.TCE
.so gen/tab_time.tbl
.ne 3.0i
.PSPIC -C figs/fig_time.eps 6.4i
.FC 2
vOTU-edge F1 (left) and cluster agreement (right) against all-vs-all wall time.
.FE
.TC 6
Choosing \f[CR]--af-min-ident\f[R] on the disjoint INPHARED tuning set (3,000 genomes; AF MAE in pp).
.TCE
.so gen/tab_tune.tbl
.\" =====================================================================
.H1 "Discussion"
.H2 "What the experiment shows"
.LP
Alignment remains the reference for vOTU clustering, and vclust shows it can be cheap: about 1 minute per
10,000 genomes. Among the sketch-based tools, aif with the thresholded AF gives the best vOTU agreement on both
data sets, and the most accurate ANI near the 95% threshold. It costs about twice skani's run time and a tenth
to a seventh of vclust's. aif's niches here are compact databases, at 1.03 bits per base, and calibrated intervals for ANI and
AF. With those, borderline vOTU calls can be flagged rather than decided silently. A natural extension is to
report, for each edge, the probability that ANI \(>= 95% and AF \(>= 85%.
.H2 "Limitations"
.IP \(bu 0.2i
\f[B]AF convention.\f[R] The vOTU standard uses BLAST's AF. aif's default AF measures homology at lower
identity and needs the identity threshold to match it. T = 0.88 was chosen on INPHARED only. GPD, a different
source, confirms the gain, but other data may prefer a slightly different T.
.IP \(bu 0.2i
\f[B]Scale.\f[R] The samples hold 10,000 sequences each. IMG/VR-scale clustering (millions of sequences) was not
attempted, and IMG/VR itself needs a JGI login. At that scale aif's per-reference estimators and screening
would dominate memory, and vclust's sparse design is the benchmark to beat.
.IP \(bu 0.2i
\f[B]Detection floor.\f[R] Below about 90% identity aif reports fewer pairs and over-estimates ANI. This is
irrelevant to the 95% vOTU threshold, but matters for genus-level clustering.
.IP \(bu 0.2i
\f[B]Clustering.\f[R] All tools were scored with the same greedy centroid clustering. vclust also offers other
algorithms (e.g. Leiden), which were not evaluated.
.IP \(bu 0.2i
\f[B]Interval calibration.\f[R] aif's intervals were calibrated on bacterial contigs against MUMmer. They were
not recalibrated against BLAST on viral pairs.
.\" =====================================================================
.H1 "Reproducibility"
.CS
make aif                                          # aif v1.6 (adds db build -i and --af-min-ident)
bench_real/viral_download.sh                      # GPD (EBI) and INPHARED (14 Apr 2025)
python3 -I bench_real/viral_bench.py prepare      # 10,000-sequence samples (seed 31)
python3 -I bench_real/viral_bench.py blast        # BLAST all-vs-all + CheckV-style anicalc
python3 -I bench_real/viral_bench.py prepare_tune # disjoint tuning sample (seed 32); then:
python3 -I bench_real/viral_bench.py blast tune; python3 -I bench_real/viral_bench.py tune
python3 -I bench_real/viral_bench.py tools        # aif, skani (x2), vclust, timed (THREADS=16)
python3 -I bench_real/viral_bench.py evaluate     # -> data/viral/viral_eval.json
/usr/bin/python3 scripts/viral_report.py; report_viral/build.sh   # this report
.CE
.LP
vclust 1.3.1 was installed with pip into \f[CR]data/envs/vclust\f[R]. Committed summaries are in
\f[CR]results_viral/\f[R] (evaluations, timings, tuning grid and the identifiers of the sampled sequences).
