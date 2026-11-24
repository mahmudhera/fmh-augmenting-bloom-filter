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
Biological application report, aif v1.6 \(em November 24, 2026
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
with aif. The search took 2.1 hours on 48 threads. 64% of MAGs and 75% of isolates carry at least one contig
that shares species-level identity (ANI \(>= 95%, AF \(>= 0.5) with another genus. In all, 1.75 million contigs
(7.9 Gbp) link 3,156 genera in 98,362 genus pairs, 46% of them across families. Most of this is not
contamination. These contigs carry 0.3\(en0.7\(mu the core-gene density of other contigs of the same length, and
71\(en74% of their genes are hypothetical. Isolates, which are never binned, carry them as often as MAGs. We
therefore call contamination only when a source species' foreign contigs carry a core gene. With this rule,
9.7% of MAGs have some contamination, 2.7% at least 1% and 0.20% at least 5%. For isolates the figures are
1.5%, 0.48% and 0.02%. Within the 256 species that have both, MAGs are the more contaminated in 205 species and
isolates in 22 (sign test p = 2\(mu10\(mi38). The rate rises with MGnify's CheckM2 contamination (from 4.0% to
32.7%) and with MAG fragmentation. Kraken2, run on the same references, puts 95% of the anchored MAG contigs in the source genus and 84% of
length-matched control contigs in the host genus. Among high-quality MAGs, 2,212 (1.8%) carry at least 1%
contamination from a source with a GTDB relative. The rule was fixed on a 3,000-genome pilot, with one later
change: GTDB R232 files five \f[I]Escherichia coli\f[R] or \f[I]Shigella\f[R] assemblies under other genera,
which made Enterobacteriaceae isolates look contaminated. References that are conspecific with a genome by GTDB's
own criterion now count as the host. A core-gene mixture model attributes 44% of the non-mobile shared MAG
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
(\f[CR]aif_real_genomes_report.pdf\f[R]); this screen uses its GTDB R232 database and frozen defaults.
.IP (1) 0.3i
\f[B]Shared sequence is common and mostly not contamination.\f[R] 64% of MAGs and 75% of isolates have at least
one contig at species-level identity (ANI \(>= 95%, AF \(>= 0.5) to another genus. The non-mobile ones in
MAGs carry 0.29\(en0.69\(mu the core-gene density of unflagged contigs of the same length, and 71\(en74% of their genes
are hypothetical (Table 4). Mobile-element hallmarks mark another 218,464 contigs.
.IP (2) 0.3i
\f[B]A mixture of two kinds.\f[R] Isolates show what transferred DNA looks like: their flagged contigs carry
almost no core genes. Measured against isolates and against unflagged MAG contigs, 44% (43\(en45%) of the
non-mobile flagged MAG bases look like misbinned genome fragments. That share is 68% for contigs under 2 kb and
25% for those of 20 kb or more (Table 5).
.IP (3) 0.3i
\f[B]Core-anchored contamination.\f[R] Contamination is called only when the foreign contigs of a source species
carry a core gene. 9.7% of MAGs then have some, 2.7% at least 1% and 0.20% at least 5%. For isolates the figures
are 1.5%, 0.48% and 0.02% (Table 6). Within species, MAGs exceed isolates by 7.8 pp, and in 205 of the 227
species that differ (Table 7).
.IP (4) 0.3i
\f[B]Agreement with independent evidence.\f[R] Among MAGs, the rate rises from 4.0% to 32.7% across MGnify's
CheckM2 contamination bins, and from 0.55% to 4.9% (at \(>= 1%) as MAGs become more fragmented (Tables 8 and 9).
Kraken2 on the same references puts 94.8% of anchored MAG contigs in the source genus (2.1% in the host genus)
and 84.0% of length-matched control contigs in the host genus (Table 11).
In the pilot, 93% of anchored MAG contigs were found in no other genome of the species; 76% of
accessory and about 50% of isolate contigs were not (Table 12).
.IP (5) 0.3i
\f[B]Sources.\f[R] The contaminating sequence comes mostly from abundant gut genera: \f[I]Bacteroides\f[R]
(17%), \f[I]Phocaeicola\f[R] (9%) and \f[I]Roseburia\f[R] (5%), with 901 source genera in all (Table 10). 556 MAGs
are at least 5% contaminated by aif's count, at a median CheckM2 contamination of 3.2%. Among them, 144 have
CheckM2 below 2%.
.IP (6) 0.3i
\f[B]A reference-label problem.\f[R] The pre-specified rule judged "another genus" by GTDB's labels. GTDB R232
files five \f[I]E. coli\f[R] or \f[I]Shigella\f[R] assemblies under other genera, so with that rule the isolates
looked more contaminated than the MAGs (20.7% vs 9.9%). Treating references that are conspecific with the genome
(ANI \(>= 95% over \(>= 50% of it) as the host fixes this. The change touched only Enterobacteriaceae and nothing in
the pilot (Section 2.5).
.IP (7) 0.3i
\f[B]Cost.\f[R] The search took 2.1 h on 48 threads for 707 GB of contigs (38 genomes per second), with a peak
of 89 GB of memory, mostly the memory-mapped database. The per-genome calls then took minutes.
.so methods.ms
.\" =====================================================================
.H1 "What the shared sequence is"
.H2 "Contig profiles"
.LP
Table 4 compares the flagged contigs (foreign, not mobile) with unflagged contigs of the same length and type of
genome.
.IP \(bu 0.2i
\f[B]Core genes.\f[R] Flagged MAG contigs carry 0.0080\(en0.0165 core genes per kb, against 0.024\(en0.027 for
unflagged contigs. The deficit grows with length, from 0.69\(mu the unflagged rate under 2 kb to 0.29\(mu at 20 kb or
more. Flagged isolate contigs carry 0.06\(en0.26\(mu the rate of their unflagged contigs.
.IP \(bu 0.2i
\f[B]Hypothetical proteins.\f[R] 71\(en74% of the genes on flagged MAG contigs are "hypothetical protein", against
48\(en55% on unflagged contigs; isolates show the same contrast (72\(en79% vs 33\(en60%).
.IP \(bu 0.2i
\f[B]Mobile elements.\f[R] Another 192,188 flagged MAG contigs and 26,276 isolate contigs carry a mobile-element
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
the 128 parts (Table 5). The share is 44% (43\(en45%) overall. It is 68% for contigs under 2 kb and 55% at
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
\f[B]MAGs.\f[R] 9.7% have some core-anchored contamination, 2.71% at least 1% and 0.20% at least 5%. For
high-quality MAGs (completeness \(>= 90%) the figures are 6.8%, 1.76% and 0.08%; for medium-quality ones 12.2%, 3.51%
and 0.30%.
.IP \(bu 0.2i
\f[B]Isolates.\f[R] 1.5%, 0.48% and 0.02%.
.IP \(bu 0.2i
\f[B]Pilot and held-out.\f[R] The held-out MAGs (9.8%, 2.72%, 0.20%) come out somewhat above the pilot MAGs
(6.9%, 1.58%, 0.09%), and so do the held-out isolates (1.6% vs 0.6%). The pilot took the first 3,000 genomes of
the metadata table, 2,294 of them species representatives. It was a tuning set, not a random sample.
.IP \(bu 0.2i
\f[B]Shared sequence\f[R], by contrast, is as common in isolates as in MAGs (75% vs 64% with any; 19% vs 26% with
at least 1% non-mobile).
.LP
In all, 121,372 contigs (615 Mbp) are anchored contamination: 0.092% of MAG bases, against 1.2% of bases implied
by CheckM2's mean MAG contamination. aif sees only contaminants whose species, or a species-level relative, is
in GTDB, and anchoring misses about half of low-level contamination (Section 5.4). The mixture estimate of
Section 3.2 (about 2.8 Gbp, 0.43% of MAG bases) is the better aggregate figure.
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
\f[B]Contamination.\f[R] MAGs have more in 205 species and isolates in 22, with 29 ties. The weighted difference is
+7.8 pp (any), +2.2 pp (\(>= 1%) and +0.1 pp (\(>= 5%), with sign-test p-values of 10\(mi38 to 10\(mi19.
.IP \(bu 0.2i
\f[B]Shared sequence.\f[R] Any shared sequence is as common in isolates as in MAGs of the same species (123 vs 126
species, p = 0.9). Isolates even lead in the weighted difference (\(mi10.1 pp), perhaps because their more
complete assemblies keep more plasmids and islands. Non-mobile shared sequence above 1% is only modestly higher in MAGs (+3.5 pp).
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
from 4.0% to 32.7%, and at \(>= 1% from 0.72% to 14.6% (Table 8, Figure 1). The Spearman correlation of the two
estimates over all MAGs with a species-level host is 0.25. They measure different things: CheckM2 sees any
foreign marker genes, aif only sequence with a GTDB relative at \(>= 95% ANI. Still, the gradient is steep and
monotone.
.LP
Contamination also rises with fragmentation (Table 9): from 0.55% at \(>= 1% for MAGs in under 50 contigs to 4.9%
for those in 300 or more. It is lower in high-quality MAGs (1.76%) than in MAGs of 50\(en90% completeness
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
The anchored contamination of MAGs (120,242 contigs, 610 Mbp) comes from 901 source genera, 24% of the bases
from another family (Table 10). The leading sources are the most abundant genera of the gut:
\f[I]Bacteroides\f[R] (17% of bases, in 2,236 MAGs), \f[I]Phocaeicola\f[R] (9%), \f[I]Roseburia\f[R] (5%) and
\f[I]Blautia_A\f[R] (4%). A misbinned contig is most likely to come from an abundant co-assembled organism.
.LP
The most contaminated MAGs are clear cases. A 17.4% contamination of the \f[I]Sarcina perfringens_B\f[R]
species representative MGYG000833322 (CheckM2: 1.2%) comes from \f[I]Clostridium butyricum\f[R] and carries 52
core genes. Overall, 76 of the 6,461 species representatives (1.2%) have at least 1% contamination
(\f[CR]results_uhgg/species_reps_contaminated.tsv\f[R]). 556 MAGs are at least 5% contaminated by aif's count,
at a median CheckM2 contamination of 3.2% (aif's median: 6.4%), and 144 of them have CheckM2 below 2%.
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
\f[B]Anchored contigs.\f[R] Kraken2 puts 94.8% of the anchored MAG contigs in the source genus that aif found,
and 2.1% in the host genus. For isolates it is 83.2% and 6.8% (1,130 contigs).
.IP \(bu 0.2i
\f[B]Other flagged contigs.\f[R] Accessory, mobile and non-coding MAG contigs go to the source genus 67\(en74% of
the time and to a third genus 11\(en18% of the time. Sequence present in several genera has no single right genus,
and LCA moves it up the tree or to another genus.
.IP \(bu 0.2i
\f[B]Controls.\f[R] 84.0% of MAG control contigs go to the host genus and 1.6% to the source genus. Isolate
controls agree with the host genus only 45.3% of the time, because 73% of them come from
\f[I]Enterobacteriaceae\f[R] hosts, where Kraken2 meets the same conflicting labels (Section 2.5). For example,
it places \f[I]E. coli\f[R] contigs in \f[I]Escherichia\f[R] when the genome's host vote went to ECMA0423. With
\f[I]Enterobacteriaceae\f[R] hosts excluded, agreement is 73.1%.
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
In the pilot, each flagged contig was searched against up to 300 other genomes of its species (Table 12).
Transferred DNA in a species' gene pool should recur in other strains, while contamination is specific to one
assembly.
.IP \(bu 0.2i
93% of anchored MAG contigs occur in no other genome of the species (90 contigs).
.IP \(bu 0.2i
76% of accessory MAG contigs and 76% of mobile MAG contigs occur in no other genome; 22% recur in another study.
.IP \(bu 0.2i
About half of the accessory, mobile and non-coding isolate contigs recur (45\(en46%).
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
MAGs with anchored contamination, the source genus was also binned from the same sample in 43.7% of cases,
against 32.2% expected from other samples of the same study (1.36\(mu, z = 44). MAGs with shared sequence only show
almost the same enrichment (1.28\(mu, z = 76): horizontal transfer is also local, between organisms that live
together (Table 13). The test confirms co-occurrence, but cannot separate the two explanations.
.LP
\f[B]GC.\f[R] Where host and source differ by at least 5 pp of GC, only 22% of anchored MAG contigs are closer
to the source's GC (4,509 contigs), against 98% of the contaminant contigs of the simulated bins. Binners
cluster on composition, so the contaminants that survive binning are those with host-like GC. Isolates, which
are not binned, show 64% (47 contigs). GC therefore cannot validate contamination in MAGs.
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
present at 44% of 0.5% contaminations, at 51% of 1% and 62% of 2%, and at 84% and 96% at 5% and 10%. Anchoring
therefore misses about half of the contaminations near 1%, but few above 5%. Requiring two core genes (strict)
costs a further 10\(en20 pp at low levels.
.TC 14
Probability that a simulated contamination carries at least one (two) core genes.
.TCE
.so gen/tab_sens.tbl
.\" =====================================================================
.H1 "Shared sequence and resistance genes"
.LP
Apart from contamination, the screen maps recent cross-genus exchange in the gut. 1.75 million contigs (7.9 Gbp)
in 183,878 genomes share \(>= 95% identity with another genus. They link 3,156 genera in 98,362 directed genus pairs,
46% of them across families and 21% across phyla (Table 15). The densest exchange is within the
\f[I]Bacteroidales\f[R] (\f[I]Bacteroides\f[R], \f[I]Phocaeicola\f[R], \f[I]Parabacteroides\f[R],
\f[I]Prevotella\f[R]) and within the \f[I]Lachnospiraceae\f[R]. Across phyla, \f[I]Bacteroidota\f[R] and
\f[I]Bacillota\f[R] exchange most: 3,563 \f[I]Phocaeicola\f[R] genomes share sequence with \f[I]Roseburia\f[R],
3,237 with \f[I]Blautia_A\f[R] and 2,461 with \f[I]Faecalibacterium\f[R].
.LP
\f[B]Resistance genes.\f[R] Of the 5,314 AMRFinderPlus genes on species representatives, 382 (7.2%) lie on
shared contigs: 231 on mobile ones and 151 on accessory ones. Only 2 lie on anchored contamination. 222
representatives carry a resistance gene on sequence shared with another genus.
.LP
\f[B]Residual label effects.\f[R] Some \f[I]Escherichia\f[R] pairs in Table 15 (with G047199115 and
"\f[I]Enterobacter\f[R]") still come from the mislabelled references of Table 3. They passed the 50% rule
because each covers less than half of a given genome. 1.6% of the flagged contigs (0.13 Gbp) and 0.6% of the
anchored bases have one of those nine references as their best hit.
.TC 15
Most frequent genus pairs of shared sequence (number of host-genus genomes with shared sequence from the genus).
.TCE
.so gen/tab_pairs.tbl
.\" =====================================================================
.H1 "Cost"
.LP
The search ran the 128 parts one at a time on 48 threads, with the 82.5 GB GTDB database in the page cache. Each
part took 57\(en62 s, 7,566 s in all for 55.7 million contigs and 707 GB of sequence: 7,400 contigs or
38 genomes per second. Peak memory was 89 GB, mostly the memory-mapped database. The per-genome calls take a
few minutes over all parts. Kraken2 assigned the 3.64 million flagged and control contigs (19.0 Gbp) in 137 s
on 48 threads (peak 87 GB, mostly its 77 GB hash table): 139 Mbp/s against aif's 93 Mbp/s. Kraken2 returns a taxon
for each contig, while aif also returns ANI and AF for each of the top hits.
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
mixture estimate puts about 4.6\(mu more contamination-like sequence in MAGs than the per-genome calls.
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
/usr/bin/python3 scripts/uhgg_report.py; report_uhgg/build.sh   # this report
.CE
.LP
The pilot is \f[CR]UHGG_PILOT=3000\f[R] for the screen and \f[CR]UHGG_PILOT=1\f[R] for the summary. Committed
results (\f[CR]results_uhgg/\f[R]): the summary, per-genome calls (\f[CR]genomes.tsv.gz\f[R]), contamination
sources, shared genus pairs, contaminated species representatives, the label-inconsistent references, and the
pre-specified rule's summaries (\f[CR]prespecified/\f[R]).
