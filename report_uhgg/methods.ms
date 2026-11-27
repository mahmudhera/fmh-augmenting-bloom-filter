.\" Methods sections of the UHGG report (included by main.ms)
.H1 "Data and screen"
.H2 "The catalogue"
.LP
The Unified Human Gastrointestinal Genome catalogue v3.0 (MGnify, released September 2026) holds 286,621
genomes: 276,463 metagenome-assembled genomes (MAGs) and 10,158 isolate genomes, in 6,461 species. Every genome
passed MGnify's quality control: CheckM2 completeness \(>= 50%, contamination \(<= 5%, and GUNC (genomes that
GUNC flagged as chimeric were removed). Each genome comes as a Prokka GFF3 file with the assembly appended.
Two properties of the annotation matter here:
.IP \(bu 0.2i
every genome has Prokka gene calls with product names;
.IP \(bu 0.2i
MGnify's mobilome annotation (\f[CR]mobile_element_type\f[R] on CDS, plasmid and prophage regions) and its
AMRFinderPlus genes exist only for the 6,461 species representatives.
.LP
Mobile and core genes were therefore identified from Prokka product names, which every genome has:
.IP \(bu 0.2i
\f[B]mobile hallmarks\f[R]: transposases, integrases, recombinases and resolvases, excisionases, conjugation and
type IV secretion proteins, relaxases and mobilization proteins, phage structural proteins (terminase, capsid,
portal, holin, baseplate, tail fibre, tape measure), plasmid proteins and antirestriction proteins;
.IP \(bu 0.2i
\f[B]core genes\f[R]: ribosomal proteins, RNA polymerase \(*a/\(*b/\(*b\(fm, EF-Tu/Ts/G, IF-1/2/3, gyrase A/B,
RecA, Ffh, SecY and the aminoacyl-tRNA ligases. These are single-copy, vertically inherited genes that every
bacterium carries and that rarely move between genera.
.H2 "The search"
.LP
Each assembly's contigs of at least 500 bp (55.7 million contigs, 707 GB) went into one of 128 FASTA parts, with
genomes kept whole. Each part was searched with \f[CR]aif db search\f[R] against the GTDB R232 database of 199,923 species representatives
(82.5 GB; Section 11 of the real-genome report), with the top 5 hits per contig and aif v1.7: the frozen v1.6
defaults plus the AF fix of Section 7.5. The screen was first run with v1.6 and rerun with v1.7; the rule was
developed with v1.6 on the pilot. The rerun moved the contamination rates of Table 6 by at most 0.23 pp and the
shared-sequence rates by at most 1.2 pp (isolates with \(>= 1% non-mobile shared sequence, 18.7% \(-> 19.8%).
Searches ran on 48 threads, with the database in the page cache. Table 1 summarizes the catalogue and the search.
.TC 1
The catalogue and its hits to GTDB R232 (host genus from: self, the genome's own GTDB copy; species, the majority
genus of species-level hits; weak, of all hits).
.TCE
.so gen/tab_data.tbl
.H2 "Host genus and foreign contigs"
.LP
The rule was fixed on a 3,000-genome pilot (the first 3,000 genomes of the metadata table) and on the simulated
GTDB bins of the real-genome report, before the rest of the catalogue was scored:
.IP \(bu 0.2i
\f[B]Self reference.\f[R] Many UHGG genomes, or near-clones of them, are GTDB representatives. A reference that is the
top hit at ANI \(>= 99.9% over \(>= 50% of a genome's bases is the genome itself; it defines the host
genus, and its hits are set aside (whatever contamination the genome has, its copy in GTDB has too).
.IP \(bu 0.2i
\f[B]Host genus\f[R] otherwise: the base-weighted majority genus of the top hits at ANI \(>= 95%, when those cover
\(>= 20% of the genome ("species" host); else the majority of all top hits ("weak").
.IP \(bu 0.2i
\f[B]Host-equivalent references\f[R]: references whose hits at ANI \(>= 95% cover \(>= 50% of the genome's bases.
By GTDB's own species criterion (ANI \(>= 95%, AF \(>= 0.5) they are conspecific with the genome, whatever genus
label they carry. This rule was added after the held-out screen (Section 2.5).
.IP \(bu 0.2i
\f[B]Foreign contig\f[R]: the best remaining hit has ANI \(>= 95% and AF \(>= 0.5 to another genus, and no hit at
ANI \(>= 95% to the host genus or to a host-equivalent reference is among the top 5. On the simulated bins this
rule recalls 95% of contaminant bases with 0.2% of host bases flagged.
.IP \(bu 0.2i
\f[B]Class\f[R]: \f[I]mobile\f[R] if the contig carries a mobile-hallmark CDS (or, for species representatives,
MGnify's mobilome covers half of it); \f[I]non-coding\f[R] if it has no CDS; otherwise \f[I]core\f[R] if it carries
a core gene and \f[I]accessory\f[R] if not. Core plus accessory is called \f[I]non-mobile\f[R].
.LP
A foreign contig shares species-level identity with another genus. That can mean two things: the contig was
misbinned from another organism (contamination), or the sequence itself moved between genera recently enough to
stay \(>= 95% identical (horizontal transfer). The analyses below separate the two.
.H2 "From the first rule to core anchoring"
.LP
The pilot was first scored with a simpler rule, which counted all non-mobile foreign bases as contamination
(mobile meaning MGnify's mobilome only). It put 3.3% of MAGs and 0.8% of isolates above 5% contamination.
Before the full screen, three checks on the pilot showed that this over-counts:
.IP \(bu 0.2i
the flagged contigs of MAGs were 69% hypothetical proteins, and their named genes were integrases (XerC/XerD-type),
Tn916 transposases, primases, topoisomerase III, LexA and efflux pumps: the cargo of integrative and conjugative
elements;
.IP \(bu 0.2i
where host and source genomes differ by \(>= 5 pp GC, only 19% of the flagged contigs (\(>= 5 kb) of MAGs were
closer to the source's GC, against 98% of the contaminant contigs of the simulated bins;
.IP \(bu 0.2i
single MAGs carried flagged contigs from six or more genera, each one with the host's GC.
.LP
Misbinned genome fragments should carry core genes at the genome-wide rate; transferred DNA carries almost none.
The revised rule therefore calls contamination only when it is \f[B]core-anchored\f[R]. Per genome, the
non-mobile foreign contigs are grouped by the GTDB species of their best hit. A group counts as contamination
when its contigs carry at least one core gene (strict: two), and contamination is the anchored bases over the
remaining bases. Everything else is reported as shared sequence. The rule, the hallmark and core lists, and every
threshold were fixed on the pilot before the held-out 283,621 genomes were scored, with the one exception
described next.
.H2 "Reference labels: the one change after the held-out screen"
.LP
With the pre-specified rule, the held-out isolates came out \f[I]more\f[R] contaminated than the MAGs: 20.7% had
core-anchored contamination and 4.0% had at least 5%, against 0.6% and 0% of the pilot isolates (Table 2).
Pure cultures are the screen's negative control, so this had to be explained before anything else was read. The
excess came almost entirely from isolates of \f[I]Escherichia\f[R] and \f[I]Enterobacter\f[R] (the pilot,
taken in metadata order, held 13 Enterobacter and no Escherichia isolates), and their top hits showed why:
.IP \(bu 0.2i
GTDB R232 files five \f[I]E. coli\f[R] or \f[I]Shigella flexneri\f[R] assemblies (NCBI taxonomy) as
representatives of placeholder genera (ECMA0423, G047199095, G047199115) or of \f[I]Enterobacter\f[R]. All five
are single-contig "complete" genomes of normal size (4.4\(en4.7 Mbp) with CheckM2 completeness of only
41\(en64%. That pattern fits uncorrected long-read assemblies, whose frameshifted genes would also mislead
gene-based placement; we did not verify this;
.IP \(bu 0.2i
two \f[I]Enterobacter hormaechei\f[R] assemblies are filed under \f[I]Escherichia\f[R]
("\f[I]Escherichia boydii\f[R]", "\f[I]Escherichia hormaechei_A\f[R]"), and \f[I]E. coli\f[R] contigs match
them at 97\(en98% ANI: they look chimeric.
.LP
An \f[I]E. coli\f[R] genome thus has species-level hits under four genus labels. Its host genus is whichever
label wins the vote, and its contigs that best match another label, often its core-gene-rich conserved
regions, became "foreign". The fix judges sameness by sequence rather than label: a reference that is
conspecific with the genome by GTDB's criterion is treated as the host (host-equivalent references, above).
The threshold is GTDB's AF of 0.5, not a tuned value. The fix changed nothing on the pilot (11,979 flagged contigs
before and after; 12,200 with v1.6) and removed 113,809 of 1,960,522 flagged contigs in the full screen (v1.6:
118,517 of 1,993,817). Only seven references were ever host-equivalent under another genus label, all in
\f[I]Enterobacteriaceae\f[R]: the five above and the two correctly labelled \f[I]E. coli\f[R] representatives,
for genomes whose host vote went to a placeholder genus (Table 3). With v1.6, whose AF read high, the two
\f[I]E. hormaechei\f[R] assemblies also covered half of seven genomes. Both versions
of every summary are kept (\f[CR]results_uhgg/prespecified/\f[R]), and Table 2 compares them. The final rule
brings the isolates back to the pilot's level (1.5% any, 0.02% \(>= 5%) and barely moves the MAGs (9.8% \(-> 9.6%).
.TC 2
Pre-specified and final rule (final: with host-equivalent references), all genomes. Same species: the
Mantel\(enHaenszel weighted difference in the share with any contamination (Section 4.2); core-gene mixture: Section 3.2.
.TCE
.so gen/tab_rule.tbl
.TC 3
GTDB R232 representatives that were host-equivalent (ANI \(>= 95% over \(>= 50%) to UHGG genomes of another host
genus, with the number of such genomes and their most frequent host genus.
.TCE
.so gen/tab_labels.tbl
.H2 "Validation"
.LP
None of the following uses aif's own estimates of contamination:
.IP \(bu 0.2i
\f[B]Core-gene mixture.\f[R] For non-mobile flagged contigs of MAGs, core genes per kb are compared with
unflagged contigs of the same length class (the rate of genuine genome fragments) and with flagged contigs of
isolates (transferred DNA, since isolates are not binned). The contamination-like share of flagged bases is
(MAG flagged \(mi isolate flagged) / (MAG unflagged \(mi isolate flagged), per length class, with 95% intervals
from a bootstrap over the 128 parts.
.IP \(bu 0.2i
\f[B]Isolates.\f[R] Pure cultures are not binned and should show little contamination; they are compared with
MAGs overall and within species (species with at least 3 isolates and 3 MAGs; Mantel\(enHaenszel weighted
difference and a sign test over species).
.IP \(bu 0.2i
\f[B]CheckM2.\f[R] MGnify's CheckM2 contamination is a marker- and model-based estimate, independent of sequence
identity to other genomes.
.IP \(bu 0.2i
\f[B]Kraken2\f[R] (GTDB R232 representatives, built as in Section 16 of the real-genome report) assigns every
flagged contig and a length-matched unflagged contig of the same genome (minimizer LCA: another algorithm on the
same references).
.IP \(bu 0.2i
\f[B]Recurrence\f[R] (pilot genomes): each flagged contig was searched against up to 300 other genomes of its
species (578 species with \(>= 10). Transferred DNA in a species' gene pool recurs in other strains; contamination
is specific to one assembly.
.IP \(bu 0.2i
\f[B]Same sample\f[R]: is the source genus also binned as another MAG of the same sample (by MGnify's own
labels)? The expectation comes from the other samples of the same study with the same number of MAGs and a MAG
of the host genus.
.IP \(bu 0.2i
\f[B]GC\f[R]: deviation of flagged contigs from their genome's GC, against length-matched controls.
