.\" Build: report_uhgg/build.sh   (tables/figures: /usr/bin/python3 scripts/uhgg_report.py)
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
aif on the Human Gut Genome Catalogue
.ft R
.ps 14
.vs 18p
Cross-Genus Shared Sequence and Core-Anchored Contamination
in 286,621 UHGG v3.0 Genomes against GTDB R232
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Biological application report, aif v1.7 \(em November 27, 2026
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
We searched every contig of at least 500 bp in the 286,621 genomes of the Unified Human Gastrointestinal Genome
catalogue v3.0 (276,463 MAGs and 10,158 isolates) against the 199,923 species representatives of GTDB R232
with aif. The search took 2.2 hours on 48 threads. 64% of MAGs and 74% of isolates carry at least one contig
that shares species-level identity (ANI \(>= 95%, AF \(>= 0.5) with another genus. In all, 1.73 million contigs
(7.8 Gbp) link 3,117 genera in 96,558 genus pairs, 46% of them across families. Most of this is not
contamination. These contigs carry 0.3\(en0.7\(mu the core-gene density of other contigs of the same length, and
71\(en74% of their genes are hypothetical. Isolates, which are never binned, carry them as often as MAGs. We
therefore call contamination only when a source species' foreign contigs carry a core gene. With this rule,
9.6% of MAGs have some contamination, 2.7% at least 1% and 0.20% at least 5%. For isolates the figures are
1.5%, 0.47% and 0.02%. Within the 256 species that have both, MAGs are the more contaminated in 207 species and
isolates in 21 (sign test p = 1\(mu10\(mi39). The rate rises with MGnify's CheckM2 contamination (from 3.9% to
32.6%) and with MAG fragmentation. Kraken2, run on the same references, puts 95% of the anchored MAG contigs in the source genus and 84% of
length-matched control contigs in the host genus. Among high-quality MAGs, 2,185 (1.7%) carry at least 1%
contamination from a source with a GTDB relative. The rule was fixed on a 3,000-genome pilot, with one later
change: GTDB R232 files five \f[I]Escherichia coli\f[R] or \f[I]Shigella\f[R] assemblies under other genera,
which made Enterobacteriaceae isolates look contaminated. References that are conspecific with a genome by GTDB's
own criterion now count as the host. A core-gene mixture model attributes 45% of the non-mobile shared MAG
sequence to misbinned genome fragments, so the per-genome calls are conservative.
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
This report applies aif to a real biological question: how much of the human-gut genome catalogue is shared
between genera, and how much of that is contamination. The method and benchmarks are in the real-genome report
(\f[CR]aif_real_genomes_report.pdf\f[R]); this screen uses its GTDB R232 database with aif v1.7, the frozen
v1.6 defaults plus the AF fix of Section 7.5.
.IP (1) 0.3i
\f[B]Shared sequence is common and mostly not contamination.\f[R] 64% of MAGs and 74% of isolates have at least
one contig at species-level identity (ANI \(>= 95%, AF \(>= 0.5) to another genus. The non-mobile ones in
MAGs carry 0.29\(en0.70\(mu the core-gene density of unflagged contigs of the same length, and 71\(en74% of their genes
are hypothetical (Table 4). Mobile-element hallmarks mark another 212,503 contigs.
.IP (2) 0.3i
\f[B]A mixture of two kinds.\f[R] Isolates show what transferred DNA looks like: their flagged contigs carry
almost no core genes. Measured against isolates and against unflagged MAG contigs, 45% (44\(en46%) of the
non-mobile flagged MAG bases look like misbinned genome fragments. That share is 68% for contigs under 2 kb and
25% for those of 20 kb or more (Table 5).
.IP (3) 0.3i
\f[B]Core-anchored contamination.\f[R] Contamination is called only when the foreign contigs of a source species
carry a core gene. 9.6% of MAGs then have some, 2.7% at least 1% and 0.20% at least 5%. For isolates the figures
are 1.5%, 0.47% and 0.02% (Table 6). Within species, MAGs exceed isolates by 7.8 pp, and in 207 of the 228
species that differ (Table 7).
.IP (4) 0.3i
\f[B]Agreement with independent evidence.\f[R] Among MAGs, the rate rises from 3.9% to 32.6% across MGnify's
CheckM2 contamination bins, and from 0.55% to 4.9% (at \(>= 1%) as MAGs become more fragmented (Tables 8 and 9).
Kraken2 on the same references puts 95.0% of anchored MAG contigs in the source genus (2.0% in the host genus)
and 83.6% of length-matched control contigs in the host genus (Table 11).
In the pilot, 94% of anchored MAG contigs were found in no other genome of the species; 76% of
accessory and about 50% of isolate contigs were not (Table 12).
.IP (5) 0.3i
\f[B]Sources.\f[R] The contaminating sequence comes mostly from abundant gut genera: \f[I]Bacteroides\f[R]
(17%), \f[I]Phocaeicola\f[R] (9%) and \f[I]Roseburia\f[R] (5%), with 889 source genera in all (Table 10). 549 MAGs
are at least 5% contaminated by aif's count, at a median CheckM2 contamination of 3.2%. Among them, 143 have
CheckM2 below 2%.
.IP (6) 0.3i
\f[B]A reference-label problem.\f[R] The pre-specified rule judged "another genus" by GTDB's labels. GTDB R232
files five \f[I]E. coli\f[R] or \f[I]Shigella\f[R] assemblies under other genera, so with that rule the isolates
looked more contaminated than the MAGs (20.7% vs 9.8%). Treating references that are conspecific with the genome
(ANI \(>= 95% over \(>= 50% of it) as the host fixes this. The change touched only Enterobacteriaceae and nothing in
the pilot (Section 2.5).
.IP (7) 0.3i
\f[B]Cost.\f[R] The search took 2.2 h on 48 threads for 707 GB of contigs (35 genomes per second), with a peak
of 89 GB of memory, mostly the memory-mapped database. The per-genome calls then took minutes.
.IP (8) 0.3i
\f[B]Against other tools.\f[R] On 10,005 sampled contigs scored against MUMmer, aif's foreign-contig calls (F1
0.871) come close to minimap2's (0.881), which would need about 10 days for the catalogue against aif's 2 h.
sourmash (0.844), Kraken2 (0.609), skani (0.578) and CAT (0.411) trail. minimap2 estimates ANI and AF best (0.32
and 2.73 pp of error, against aif's 0.75 and 3.58). The benchmark exposed an AF over-call in aif v1.6 (+5 pp),
fixed in v1.7 with recalibrated intervals (Section 7).
.so methods.ms
.\" =====================================================================
.H1 "What the shared sequence is"
.H2 "Contig profiles"
.LP
Table 4 compares the flagged contigs (foreign, not mobile) with unflagged contigs of the same length and type of
genome.
.IP \(bu 0.2i
\f[B]Core genes.\f[R] Flagged MAG contigs carry 0.0081\(en0.0166 core genes per kb, against 0.024\(en0.027 for
unflagged contigs. The deficit grows with length, from 0.70\(mu the unflagged rate under 2 kb to 0.29\(mu at 20 kb or
more. Flagged isolate contigs carry 0.06\(en0.24\(mu the rate of their unflagged contigs.
.IP \(bu 0.2i
\f[B]Hypothetical proteins.\f[R] 71\(en74% of the genes on flagged MAG contigs are "hypothetical protein", against
48\(en55% on unflagged contigs; isolates show the same contrast (71\(en78% vs 33\(en60%).
.IP \(bu 0.2i
\f[B]Mobile elements.\f[R] Another 186,774 flagged MAG contigs and 25,729 isolate contigs carry a mobile-element
hallmark (transposase, integrase, conjugation, phage or plasmid genes).
.LP
Misbinned genome fragments should look like unflagged contigs: same core-gene density, same share of named
genes. The flagged contigs look instead like the accessory genome: islands and integrative elements, whose
cargo is mostly hypothetical genes. Long flagged contigs are the most depleted, as expected for large
integrative elements.
.TC 4
Flagged (foreign, non-mobile) and unflagged contigs by length: core genes per kb and their ratio, and the share
of hypothetical CDS. The last column counts the flagged contigs classed as mobile.
.TCE
.so gen/tab_profile.tbl
.ne 3.0i
.PSPIC -C figs/fig_core.eps 6.4i
.FC 1
Left: core genes per kb of unflagged MAG contigs and of flagged non-mobile contigs of MAGs and isolates. Right:
MAGs with core-anchored contamination by MGnify's CheckM2 contamination (95% Wilson intervals).
.FE
.H2 "A core-gene mixture"
.LP
If the flagged non-mobile sequence of MAGs is a mixture of misbinned fragments and transferred DNA, its core-gene
density lies between that of genome fragments (unflagged MAG contigs) and that of transferred DNA (the flagged
contigs of isolates, which are not binned). The share of fragments is therefore
.CS
share = (MAG flagged \(mi isolate flagged) / (MAG unflagged \(mi isolate flagged)        (core genes per kb)
.CE
.LP
computed per length class and weighted by flagged bases. The 95% intervals come from 300 bootstrap resamples of
the 128 parts (Table 5). The share is 45% (44\(en46%) overall. It is 68% for contigs under 2 kb and 56% at
2\(en5 kb, and falls to 35% and 25% for longer contigs. Short flagged contigs are the ones a binner places by
coverage and composition alone, and mostly look like misbinned fragments. Long ones are mostly transferred
elements. This estimate is aggregate, not per genome. It assumes that isolates' transferred DNA is
representative of MAGs', and that binned fragments carry core genes at the unflagged rate.
.TC 5
Core-gene mixture of flagged non-mobile MAG sequence (bootstrap over the 128 parts).
.TCE
.so gen/tab_mixture.tbl
.\" =====================================================================
.H1 "Per-genome contamination"
.H2 "Rates"
.LP
Table 6 and Figure 2 give the per-genome calls.
.IP \(bu 0.2i
\f[B]MAGs.\f[R] 9.6% have some core-anchored contamination, 2.68% at least 1% and 0.20% at least 5%. For
high-quality MAGs (completeness \(>= 90%) the figures are 6.7%, 1.73% and 0.08%; for medium-quality ones 12.1%, 3.47%
and 0.30%.
.IP \(bu 0.2i
\f[B]Isolates.\f[R] 1.5%, 0.47% and 0.02%.
.IP \(bu 0.2i
\f[B]Pilot and held-out.\f[R] The held-out MAGs (9.7%, 2.69%, 0.20%) come out somewhat above the pilot MAGs
(6.7%, 1.58%, 0.09%), and so do the held-out isolates (1.6% vs 0.6%). The pilot took the first 3,000 genomes of
the metadata table, 2,294 of them species representatives. It was a tuning set, not a random sample.
.IP \(bu 0.2i
\f[B]Shared sequence\f[R], by contrast, is as common in isolates as in MAGs (74% vs 64% with any; 20% vs 26% with
at least 1% non-mobile).
.LP
In all, 120,242 contigs (610 Mbp) are anchored contamination: 0.091% of MAG bases, against 1.2% of bases implied
by CheckM2's mean MAG contamination. aif sees only contaminants whose species, or a species-level relative, is
in GTDB, and anchoring misses about half of low-level contamination (Section 5.4). The mixture estimate of
Section 3.2 (about 2.8 Gbp, 0.42% of MAG bases) is the better aggregate figure.
.TC 6
Per-genome calls (shared: any foreign contig, and non-mobile foreign bases \(>= 1%; contamination: any anchored
contig, and anchored bases over the remaining bases \(>= 1% and \(>= 5%; last column: anchoring with at least 2
core genes per source species; reps: species representatives).
.TCE
.so gen/tab_rates.tbl
.ne 3.0i
.PSPIC -C figs/fig_hist.eps 6.4i
.FC 2
Distribution of the per-genome fractions of non-mobile shared sequence (left) and core-anchored contamination
(right), isolates and MAGs (log scale).
.FE
.H2 "Isolates and MAGs of the same species"
.LP
Overall rates mix species. Table 7 compares isolates and MAGs within the 256 species that have at least 3 of
each (8,896 isolates, 124,780 MAGs).
.IP \(bu 0.2i
\f[B]Contamination.\f[R] MAGs have more in 207 species and isolates in 21, with 28 ties. The weighted difference is
+7.8 pp (any), +2.2 pp (\(>= 1%) and +0.1 pp (\(>= 5%), with sign-test p-values of 10\(mi39 to 10\(mi19.
.IP \(bu 0.2i
\f[B]Shared sequence.\f[R] Any shared sequence is as common in isolates as in MAGs of the same species (124 vs 125
species, p = 1.0). Isolates even lead in the weighted difference (\(mi10.1 pp), perhaps because their more
complete assemblies keep more plasmids and islands. Non-mobile shared sequence above 1% is only modestly higher in MAGs (+2.5 pp).
.LP
Binning should add contamination but not transferred DNA, and that is the pattern.
.TC 7
Isolates and MAGs of the same species (256 species with \(>= 3 of each): pooled rates, Mantel\(enHaenszel
weighted difference, and the count of species where each is higher.
.TCE
.so gen/tab_matched.tbl
.H2 "CheckM2, fragmentation and completeness"
.LP
MGnify's CheckM2 contamination is a model-based estimate from marker genes and owes nothing to identity with
other genomes. Across its five bins, 0\(en1% to 4\(en5%, the share of MAGs with core-anchored contamination rises
from 3.9% to 32.6%, and at \(>= 1% from 0.71% to 14.5% (Table 8, Figure 1). The Spearman correlation of the two
estimates over all MAGs with a species-level host is 0.25. They measure different things: CheckM2 sees any
foreign marker genes, aif only sequence with a GTDB relative at \(>= 95% ANI. Still, the gradient is steep and
monotone.
.LP
Contamination also rises with fragmentation (Table 9): from 0.55% at \(>= 1% for MAGs in under 50 contigs to 4.9%
for those in 300 or more. It is lower in high-quality MAGs (1.73%) than in MAGs of 50\(en90% completeness
(3.4\(en3.5%). Short contigs carry less coverage and composition signal, and are the ones binners misplace.
.TC 8
Core-anchored contamination of MAGs by MGnify's CheckM2 contamination (95% Wilson intervals).
.TCE
.so gen/tab_checkm.tbl
.TC 9
Contamination of MAGs by number of contigs and by CheckM2 completeness.
.TCE
.so gen/tab_frag.tbl
.H2 "Sources"
.LP
The anchored contamination of MAGs (119,105 contigs, 604 Mbp) comes from 889 source genera, 24% of the bases
from another family (Table 10). The leading sources are the most abundant genera of the gut:
\f[I]Bacteroides\f[R] (17% of bases, in 2,232 MAGs), \f[I]Phocaeicola\f[R] (9%), \f[I]Roseburia\f[R] (5%) and
\f[I]Blautia_A\f[R] (4%). A misbinned contig is most likely to come from an abundant co-assembled organism.
.LP
The most contaminated MAGs are clear cases. A 17.4% contamination of the \f[I]Sarcina perfringens_B\f[R]
species representative MGYG000833322 (CheckM2: 1.2%) comes from \f[I]Clostridium butyricum\f[R] and carries 52
core genes. Overall, 75 of the 6,461 species representatives (1.2%) have at least 1% contamination
(\f[CR]results_uhgg/species_reps_contaminated.tsv\f[R]). 549 MAGs are at least 5% contaminated by aif's count,
at a median CheckM2 contamination of 3.2% (aif's median: 6.4%), and 143 of them have CheckM2 below 2%.
.TC 10
Sources of core-anchored contamination in MAGs, by genus (share of anchored bases; MAGs: with any contamination
from the genus).
.TCE
.so gen/tab_sources.tbl
.\" =====================================================================
.H1 "Independent checks"
.H2 "Kraken2"
.LP
Kraken2 (minimizer LCA) is another algorithm on the same GTDB R232 references (Section 16 of the real-genome
report). It assigned every flagged contig and, for each coding flagged contig, a length-matched unflagged contig
of the same genome as a control (Table 11).
.IP \(bu 0.2i
\f[B]Anchored contigs.\f[R] Kraken2 puts 95.0% of the anchored MAG contigs in the source genus that aif found,
and 2.0% in the host genus. For isolates it is 82.9% and 7.1% (1,137 contigs).
.IP \(bu 0.2i
\f[B]Other flagged contigs.\f[R] Accessory, mobile and non-coding MAG contigs go to the source genus 67\(en74% of
the time and to a third genus 11\(en17% of the time. Sequence present in several genera has no single right genus,
and LCA moves it up the tree or to another genus.
.IP \(bu 0.2i
\f[B]Controls.\f[R] 83.6% of MAG control contigs go to the host genus and 1.7% to the source genus. Isolate
controls agree with the host genus only 45.2% of the time, because 76% of them come from
\f[I]Enterobacteriaceae\f[R] hosts, where Kraken2 meets the same conflicting labels (Section 2.5). For example,
it places \f[I]E. coli\f[R] contigs in \f[I]Escherichia\f[R] when the genome's host vote went to ECMA0423. With
\f[I]Enterobacteriaceae\f[R] hosts excluded, agreement is 72.8%.
.LP
The check shows that aif's genus assignments are not an artifact of sketching, Bloom filters or the HMM. An
exact-k-mer classifier puts the same contigs in the same genus. It does not separate contamination from
transfer, since both put source-genus sequence on the contig. That separation rests on the core-gene evidence.
.ne 3.2i
.TC 11
Kraken2 assignment of flagged contigs (by class) and of length-matched unflagged control contigs from the same
genomes. Source and host: the genera found by aif. Source family: the source's family, when it differs from
the host's.
.TCE
.so gen/tab_kraken.tbl
.H2 "Recurrence in other strains"
.LP
In the pilot, each flagged contig was searched against up to 300 other genomes of its species (Table 12; the
search used aif v1.6, and the counts are of the contigs that v1.7 still flags).
Transferred DNA in a species' gene pool should recur in other strains, while contamination is specific to one
assembly.
.IP \(bu 0.2i
94% of anchored MAG contigs occur in no other genome of the species (88 contigs).
.IP \(bu 0.2i
76% of accessory MAG contigs and 76% of mobile MAG contigs occur in no other genome; 22% recur in another study.
.IP \(bu 0.2i
About half of the accessory, mobile and non-coding isolate contigs recur (44\(en45%).
.LP
Anchored contigs thus behave most like one-off contamination. Accessory and mobile contigs are partly shared gene
pool. The pilot's anchored set is small, so this check is indicative.
.TC 12
Recurrence of flagged pilot contigs in up to 300 other genomes of the species (\(>= 95% ANI, AF \(>= 0.5): in a genome
of another study, of the same study only, or absent.
.TCE
.so gen/tab_recur.tbl
.H2 "Same sample and GC: two checks that do not discriminate"
.LP
\f[B]Same sample.\f[R] A contaminant's source should live in the same sample, and is often binned there too. For
MAGs with anchored contamination, the source genus was also binned from the same sample in 43.8% of cases,
against 32.2% expected from other samples of the same study (1.36\(mu, z = 44). MAGs with shared sequence only show
almost the same enrichment (1.28\(mu, z = 76): horizontal transfer is also local, between organisms that live
together (Table 13). The test confirms co-occurrence, but cannot separate the two explanations.
.LP
\f[B]GC.\f[R] Where host and source differ by at least 5 pp of GC, only 23% of anchored MAG contigs are closer
to the source's GC (4,447 contigs), against 98% of the contaminant contigs of the simulated bins. Binners
cluster on composition, so the contaminants that survive binning are those with host-like GC. Isolates, which
are not binned, show 64% (45 contigs). GC therefore cannot validate contamination in MAGs.
.TC 13
Is the source genus also binned as another MAG of the same sample? (Expected: other samples of the same study
with the same number of MAGs and a MAG of the host genus.)
.TCE
.so gen/tab_sample.tbl
.H2 "Sensitivity of anchoring"
.LP
How often does real contamination carry a core gene? We drew random contigs from 3,000 high-quality MAGs
without flagged contigs until they added up to a fixed share of a random MAG's size, 4,000 draws per level
(Table 14). At least one core gene is
present at 44% of 0.5% contaminations, at 49% of 1% and 63% of 2%, and at 84% and 96% at 5% and 10%. Anchoring
therefore misses about half of the contaminations near 1%, but few above 5%. Requiring two core genes (strict)
costs a further 17\(en22 pp at levels up to 2%.
.TC 14
Probability that a simulated contamination carries at least one (two) core genes.
.TCE
.so gen/tab_sens.tbl
.\" =====================================================================
.H1 "Shared sequence and resistance genes"
.LP
Apart from contamination, the screen maps recent cross-genus exchange in the gut. 1.73 million contigs (7.8 Gbp)
in 182,554 genomes share \(>= 95% identity with another genus. They link 3,117 genera in 96,558 directed genus pairs,
46% of them across families and 21% across phyla (Table 15). The densest exchange is within the
\f[I]Bacteroidales\f[R] (\f[I]Bacteroides\f[R], \f[I]Phocaeicola\f[R], \f[I]Parabacteroides\f[R],
\f[I]Prevotella\f[R]) and within the \f[I]Lachnospiraceae\f[R]. Across phyla, \f[I]Bacteroidota\f[R] and
\f[I]Bacillota\f[R] exchange most: 3,530 \f[I]Phocaeicola\f[R] genomes share sequence with \f[I]Roseburia\f[R],
3,209 with \f[I]Blautia_A\f[R] and 2,448 with \f[I]Faecalibacterium\f[R].
.LP
\f[B]Resistance genes.\f[R] Of the 5,314 AMRFinderPlus genes on species representatives, 370 (7.0%) lie on
shared contigs: 222 on mobile ones and 148 on accessory ones. Only 2 lie on anchored contamination. 219
representatives carry a resistance gene on sequence shared with another genus.
.LP
\f[B]Residual label effects.\f[R] Some \f[I]Escherichia\f[R] pairs in Table 15 (with G047199115 and
"\f[I]Enterobacter\f[R]") still come from the mislabelled references of Table 3. They passed the 50% rule
because each covers less than half of a given genome. 1.7% of the flagged contigs (0.13 Gbp) and 0.6% of the
anchored bases have one of those seven references as their best hit.
.TC 15
Most frequent genus pairs of shared sequence (number of host-genus genomes with shared sequence from the genus).
.TCE
.so gen/tab_pairs.tbl
.\" =====================================================================
.H1 "Accuracy and cost against other tools"
.H2 "Cost of the screen"
.LP
The search ran the 128 parts one at a time on 48 threads, with the 82.5 GB GTDB database in the page cache. Each
part took 56\(en76 s, 8,082 s in all for 55.7 million contigs and 707 GB of sequence: 6,900 contigs or
35 genomes per second (v1.6: 7,566 s). Peak memory was 89 GB, mostly the memory-mapped database. The per-genome calls take a
few minutes over all parts. Kraken2 assigned the 3.58 million flagged and control contigs (18.6 Gbp) in 224 s
on 48 threads with its database in the page cache, 108 s of it classifying (peak 87 GB, mostly its 77 GB hash
table; 137 s and 85 s in the v1.6 run, on this shared host): 83 Mbp/s overall and 173 Mbp/s while classifying,
against aif's 88 Mbp/s. Kraken2 returns a taxon
for each contig, while aif also returns ANI and AF for each of the top hits.
.H2 "A benchmark on the screen's own contigs"
.LP
How close are aif's calls to alignment, and what would the alternatives cost on the whole catalogue? We sampled
contigs from the screen and scored seven tools against MUMmer (\f[CR]bench_real/uhgg_compare.py\f[R]).
.IP \(bu 0.2i
\f[B]Sample.\f[R] 8 of the 128 parts at random, then 1,000 MAGs and 300 isolates at random among their held-out
genomes; every contig of at least 500 bp of a drawn genome forms the pool (237,338 contigs, 3.64 Gbp). The pool
is stratified by the screen's call (anchored, accessory, mobile, non-coding, label conflict, flagged by skani or
Kraken2 only, near miss, ambiguous best hit, host only, no species-level hit) and by length (< 2 kb, 2\(en10 kb,
\(>= 10 kb). Nested subsets S1 \(sb S2 \(sb S3 of 1,018, 3,013 and 10,005 contigs (11, 33 and 114 Mbp) take
the first contigs of one seeded order in every cell. Inverse inclusion probabilities give catalogue-level
("weighted") estimates.
.IP \(bu 0.2i
\f[B]Truth.\f[R] Every top hit of every tool (aif 5, skani and sourmash 10, minimap2 and LexicMap 20 per
contig) was aligned with \f[CR]nucmer \-\-maxmatch\f[R]: 105,950 contig\(enreference pairs. AF is the union of
aligned contig bases over the contig length, ANI the length-weighted identity. The primary convention keeps
\f[CR]delta-filter \-q\f[R] alignments, as every earlier aif benchmark; the one-to-one filter
(\f[CR]\-1\f[R], as dnadiff) and alignments of \(>= 90% identity are reported too.
.IP \(bu 0.2i
\f[B]Measures.\f[R] Estimation: each tool's top hit other than the genome's self reference, on pairs with true
ANI \(>= 0.90. Search: for the 7,268 contigs whose best reference is a species-level match, does the tool report a
hit, and is its top hit that reference? Decision: the screen's foreign-contig rule applied to each tool's hits
(Kraken2 and CAT: their genus against the host genus), scored against the same rule applied to MUMmer.
.IP \(bu 0.2i
\f[B]Tools.\f[R] aif v1.7 (top 5; Section 7.5); skani 0.4.0 (\f[CR]search \-\-qi\f[R]); sourmash branchwater
\f[CR]manysearch\f[R] (k = 21, scaled 64, containment \(>= 0.05); minimap2 2.31 (\f[CR]\-x asm20 \-c \-N 50\f[R],
all of GTDB as a 25 Gbp split index); LexicMap 0.10.1 (\f[CR]\-n 100\f[R]); Kraken2 2.1.3 (database capped at
aif's size); CAT_pack 7 (DIAMOND against the GTDB R232 proteins). All ran at 48 threads pinned to the same
cores, with their databases in the page cache where they fit. LexicMap ran on S1 only, which took 6.3 h.
.IP \(bu 0.2i
\f[B]Cost projection.\f[R] Time is taken as a fixed cost plus a rate per base. The fixed cost comes from the
S1\(enS3 fit and the rate from the pool, a random sample of genomes; CAT and LexicMap, without a pool run, are
projected from the subsets. Projected this way, aif's full screen comes to 1.75 h, against 2.24 h measured.
.H2 "Accuracy"
.LP
Table 16 gives the S3 results and Table 17 those of S1, the only subset with LexicMap.
.IP \(bu 0.2i
\f[B]minimap2 is the most accurate.\f[R] Its top-hit ANI error is 0.32 pp against aif's 0.75, its AF error 2.73
pp against 3.58, and its foreign-call F1 0.881 against 0.871. LexicMap has the lowest AF error (1.4 pp on S1)
and the highest weighted precision (0.95), at a lower recall (0.81; weighted 0.85).
.IP \(bu 0.2i
\f[B]aif is close behind on the decisions\f[R], and finds the best reference for 89.5% of contigs, more than any
other tool. sourmash comes next (F1 0.844, ANI error 1.23 pp, no AF). skani reports a hit for only 66.5% of the
contigs (F1 0.578); Kraken2 (0.609) flags 5.4% of the catalogue against a true 3.1%, and CAT (0.411) assigns a
genus to only 37% of contigs.
.IP \(bu 0.2i
\f[B]Conventions.\f[R] The ranking holds under all three (Table 18). On the pairs that both report, aif v1.7's
AF error is below skani's under each convention; v1.6's was above it under the one-to-one and 90%-identity
conventions.
.LP
The strata come from aif's own screen calls and the foreign rule was fixed on aif's pilot. The weights correct
the sampling, but the design centres on the decisions aif makes. The truth pairs are limited to the tools' top
hits, which include minimap2's search of all of GTDB.
.TC 16
Accuracy on S3 (10,005 contigs; MUMmer, \f[CR]delta-filter \-q\f[R]). MAE: mean absolute error of the top hit
(true ANI \(>= 0.90). Has a hit, top hit best: of the 7,268 contigs with a species-level best reference. Foreign
call: precision, recall and F1 of the screen's rule, unweighted and weighted to the catalogue (true foreign
rate 3.1%).
.TCE
.so gen/tab_cmp_acc.tbl
.TC 17
Accuracy on S1 (1,018 contigs), the subset that LexicMap was run on; columns as in Table 16.
.TCE
.so gen/tab_cmp_s1.tbl
.TC 18
AF and ANI error of the top hit by MUMmer convention (S3), and AF error on the (contig, reference) pairs that aif
and skani both report.
.TCE
.so gen/tab_cmp_conv.tbl
.H2 "Cost"
.LP
Table 19 and Figure 3 give the cost of each tool, measured on S1\(enS3 and the pool and projected to the 706 Gbp
of the catalogue.
.IP \(bu 0.2i
\f[B]aif, Kraken2 and skani\f[R] would each take 1.4\(en2.5 h on 48 threads, with 65\(en82 GB of memory.
sourmash would take 19 h, in 5 GB.
.IP \(bu 0.2i
\f[B]minimap2\f[R] would take 236 h (10 days) at a 177 GB peak, after building a 2.08 TB index in 12.9 h (347 GB
of memory). That is feasible on one node, at about 105 times aif's measured 2.2 h. On S1\(enS3 its time is
almost all index loading (about 2.6 h per pass over the 27 index parts) and hardly changes with query size, so the
per-base rate came from the pool run. Its temporary split files would need about 0.5 TB.
.IP \(bu 0.2i
\f[B]CAT\f[R] would take about 100 days (175 from the S2\(enS3 slope alone), most of it in DIAMOND.
.IP \(bu 0.2i
\f[B]LexicMap\f[R] took 6.3 h for the 11 Mbp of S1, almost all of it waiting on random reads of its 161 GB genome
store, which does not fit in memory beside its 344 GB of seeds; the node's disks are hard drives. Even 100
times faster on solid-state storage, the catalogue would take months.
.LP
aif thus makes the screen's decisions nearly as well as an aligner, at the cost of a sketch tool, and reports
ANI, AF and their intervals for every hit.
.TC 19
Cost: index size and build time (\(dg: wall time minus pauses, on a shared disk, so an upper bound), search time
at 48 threads on S3 and the pool, the projection to the full catalogue and the peak memory of the search.
.TCE
.so gen/tab_cmp_cost.tbl
.ne 3.0i
.PSPIC -C figs/fig_cost.eps 6.4i
.FC 3
Projected time for the full catalogue against the foreign-contig F1 (left) and the top-hit AF error (right) on
S3; LexicMap on S1. sourmash, Kraken2 and CAT report no AF.
.FE
.H2 "AF over-call and the v1.7 fix"
.LP
The benchmark first ran aif v1.6, whose top-hit AF read 5.1 pp high on S3 (Table 21). We diagnosed and fixed it
on a separate tuning set, T: 3,002 contigs that come next in the same strata after S3's and share no contig with
S1\(enS3 (\f[CR]bench_real/uhgg_af.py\f[R]).
.IP \(bu 0.2i
\f[B]Where the excess lies.\f[R] Base by base over 2,550 top-hit pairs, aif called 84.3% of the contig bases
aligned and MUMmer 81.1%. Three quarters of aif's extra 3.6% are gaps bridged between two MUMmer alignments (2.7%)
and a quarter extensions past an alignment's end (0.9%); there are no isolated islands. They sit where aif's own
identity estimate is low: in windows aif puts below 86% identity, 80% of the aligned bases have no MUMmer
alignment, against 1% at 95% or more.
.IP \(bu 0.2i
\f[B]Mechanism.\f[R] The lowest identity levels of aif's HMM (about 85%) expect a hit rate close to the Bloom
filters' false-positive rates. An unaligned stretch inside an alignment then costs almost nothing as "aligned at
85%", while leaving and re-entering the aligned state costs about 14 nats (two switches between bins of 32 probes,
about 300 bp apart). The posterior there follows the
switching prior, not the probes.
.IP \(bu 0.2i
\f[B]Fix.\f[R] AF now counts the aligned posterior mass at levels within 0.05 of the local identity: the median
FP-corrected containment of the aligned-looking windows of 64 probes within three windows either side. A window
that favours its best aligned level over unaligned by 3 nats, or has such a neighbour, lowers the floor to its
own identity, so that divergent blocks keep their AF. The HMM is unchanged, and so are the hits and ANI;
\f[CR]\-\-af-local 0\f[R] reproduces v1.6.
.LP
Table 20 shows why the local form was chosen. A flat threshold (AF over identity \(>= 0.89) removes the bias on
high-identity pairs, but cuts AF by 33 pp at true ANI 0.80\(en0.88, where the whole contig sits at the low
levels. The local threshold alone keeps those, but cuts real divergent blocks between conserved ones (\(mi34 pp on
the simulated tuning set's contigs of alternating 100% and 90% identity), which the evidence rule restores. v1.7 trades a little of the
local threshold's gain on T for the lowest error on the real tuning set and a far smaller loss on heterogeneous
contigs. All settings were chosen on these tuning sets.
.LP
On held-out data (Table 21) the AF error falls by 30% on S3 (5.13 to 3.58 pp) and by 20% on the real-genome
sets of the earlier report (2.45 to 1.97 pp; fresh drafts 2.53 to 1.97). There, genome-level ANI aggregated from
contigs is unchanged and genome-level AF error halves (Section 17 of that report). It rises on the simulated benchmarks
(5.30 to 5.83 pp): on contigs of heterogeneous identity, \f[CR]hetblock\f[R] loses 5.8 pp and \f[CR]hetgamma\f[R]
4.4 pp, and uniformly mutated contigs at 90\(en99% identity 0.5\(en1 pp. That failure was first seen on the
held-out simulated set, and it shaped the evidence rule. The identity maps of \f[CR]\-\-segments\f[R] still
come from the HMM's aligned segments and are unchanged.
.TC 20
AF models on the three tuning sets (pp): UHGG tuning set T (top hit, \f[CR]delta-filter \-q\f[R]; true ANI \(>= 0.90,
and the top-5 pairs at true ANI 0.80\(en0.90), the real and simulated tuning sets (every contig, macro over sets or
conditions), and the AF bias on simulated contigs of alternating 100% and 90% identity (\f[CR]hetblock\f[R]).
.TCE
.so gen/tab_cmp_afvar.tbl
.TC 21
AF error of aif v1.6 and v1.7 (UHGG: top hit, true ANI \(>= 0.90; real and simulated: every contig, macro over sets
or conditions as in the real-genome report).
.TCE
.so gen/tab_cmp_affix.tbl
.H2 "Intervals"
.LP
aif's intervals are split-conformal: the statistical SE times a multiplier fitted per stratum on tuning data.
Fitted on pairwise queries of the real tuning set, they covered 92% of the true ANI and 87% of the true AF on S3
at a nominal 95%. Database search compares a contig with references at all distances and its errors differ, and
one pooled table over-covered pairwise queries (96\(en97%) while still under-covering S3 (94% and 92%). v1.7
therefore has two built-in tables: pairwise queries use one refitted on the real tuning set, and database search
one fitted on T (9,880 hits). Table 22 gives the held-out coverage: 95.1% (ANI) and 94.2% (AF) for S3's top hits,
and 94.6% and 94.3% on the real genomes, whose AF intervals became narrower (median half-width 2.60 to 1.98 pp)
because AF is more accurate.
.TC 22
Coverage of the 95% intervals on held-out data (UHGG: top hit, \f[CR]delta-filter \-q\f[R], true ANI \(>= 0.90;
real and simulated: as in the real-genome report) and the median AF half-width.
.TCE
.so gen/tab_cmp_ci.tbl
.\" =====================================================================
.H1 "Discussion"
.H2 "What the screen shows"
.LP
Species-level identity with another genus is common in gut genomes. Taken at face value it would put
contamination in two-thirds of the catalogue, isolates included. Core genes separate the two meanings. Misbinned
fragments carry core genes at the genome-wide rate. Transferred DNA carries almost none. Anchoring on them gives
contamination rates that are low for isolates and higher for MAGs, within species as well as overall. These
rates rise steeply with CheckM2's own estimate and with fragmentation. Kraken2 puts the anchored contigs in the
same source genus, and in the pilot they are specific to one assembly. Beyond a CheckM2-style estimate, aif names the contaminating contigs and their source species.
.LP
The screen also shows a hazard of identity-based methods: reference labels. Five mislabelled
\f[I]E. coli\f[R]/\f[I]Shigella\f[R] assemblies in GTDB R232 were enough to make 20% of all isolates look
contaminated. The fix, judging sameness by sequence (ANI \(>= 95% over \(>= 50%), uses GTDB's own species criterion
and needed no tuning. It is now part of the rule.
.H2 "Limitations"
.IP \(bu 0.2i
\f[B]Detection scope.\f[R] Only contaminants with a species-level GTDB relative are seen. Contamination by
novel species, or by the host's own genus, is invisible to this screen.
.IP \(bu 0.2i
\f[B]Anchoring is conservative.\f[R] About half of 1% contaminations carry no core gene (Section 5.4), and the
mixture estimate puts about 4.7\(mu more contamination-like sequence in MAGs than the per-genome calls.
.IP \(bu 0.2i
\f[B]Gene names.\f[R] Core and mobile genes come from Prokka product names. Mobile elements without a named
hallmark are counted as accessory, and divergent core genes named "hypothetical protein" are missed.
.IP \(bu 0.2i
\f[B]One rule change after the held-out screen.\f[R] The host-equivalent references (Section 2.5) were added
after seeing the held-out isolates. Both versions are reported, and the change is confined to Enterobacteriaceae.
.IP \(bu 0.2i
\f[B]The pilot\f[R] was the first 3,000 genomes in metadata order, not a random sample, and the recurrence check
used only the pilot.
.IP \(bu 0.2i
\f[B]Benchmark design.\f[R] The tool comparison samples contigs by aif's own screen calls, and the foreign rule
was fixed on aif's pilot (with v1.6). The weights correct the sampling but not this focus. The AF fix was found
with this benchmark; it was tuned on a separate set, but its failure on heterogeneous simulated contigs was first
seen on held-out data.
.IP \(bu 0.2i
\f[B]No truth set.\f[R] There is no gold standard for contamination in real MAGs. The checks are independent
but indirect. The simulated GTDB bins of the real-genome report remain the only labelled test.
.\" =====================================================================
.H1 "Reproducibility"
.CS
bench_real/uhgg_download.sh                       # 286,621 Prokka GFF3 (with assemblies) from MGnify
python3 -I bench_real/uhgg_screen.py prepare      # 128 FASTA parts, contig table (CDS classes)
python3 -I bench_real/uhgg_screen.py search       # aif db search vs GTDB R232 (48 threads, top 5)
python3 -I bench_real/uhgg_screen.py genomes      # per-genome calls  (bench_real/uhgg_run.sh: all four)
/usr/bin/python3 -I bench_real/uhgg_summary.py stats    # -> results_uhgg/summary.json and tables
/usr/bin/python3 -I bench_real/uhgg_summary.py kraken   # Kraken2 check -> results_uhgg/kraken2_check.json
python3 -I bench_real/uhgg_compare.py pool|strata|run|truth|evaluate|project   # tool benchmark (Section 7)
python3 -I bench_real/uhgg_af.py tune|run|truth|score|diag      # AF diagnosis and tuning set T
python3 -I bench_real/af_local_check.py [\-\-tune]             # v1.7 AF on the earlier benchmarks
python3 bench_real/calibration.py fit|evaluate                 # interval tables (query, db search)
/usr/bin/python3 scripts/uhgg_report.py; /usr/bin/python3 scripts/uhgg_compare_report.py
report_uhgg/build.sh                                            # this report
.CE
.LP
The pilot is \f[CR]UHGG_PILOT=3000\f[R] for the screen and \f[CR]UHGG_PILOT=1\f[R] for the summary. Committed
results (\f[CR]results_uhgg/\f[R]): the summary, per-genome calls (\f[CR]genomes.tsv.gz\f[R]), contamination
sources, shared genus pairs, contaminated species representatives, the label-inconsistent references, and the
pre-specified rule's summaries (\f[CR]prespecified/\f[R]).
