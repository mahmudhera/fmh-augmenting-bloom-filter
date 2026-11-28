.\" Build: report_hgt/build.sh   (tables and figures: /usr/bin/python3 -I scripts/hgt_report.py)
.EQ
delim $$
gsize 11
.EN
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
Recent Gene Exchange Across the Bacterial Tree
.ft R
.ps 14
.vs 18p
An aif Atlas of DNA Shared at \(>= 99% Identity
Among the 199,923 Species Representatives of GTDB R232
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Biological application report, aif v1.7 \(em November 28, 2026
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
Horizontal gene transfer is recorded in genomes as stretches of nearly identical DNA shared by distant lineages.
Finding them needs local alignment of every genome against every other, which limited the landmark study of
Smillie et al. (2011) to 2,235 genomes. aif makes the search a database query: every genome is cut into 5 kb
windows, and each window is searched against all 199,923 GTDB R232 species representatives. For the 10,296
complete isolate genomes (44.5 Gbp, 8.9 million windows) this took 39 minutes on 48 threads, against an
estimated 12 days for minimap2. Against BLAST, the standard of the earlier study (\(>= 99% identity over
\(>= 500 bp), aif's calls of shared DNA between genera are 98% precise at the level of genome pairs and recover 76%
of BLAST's pairs; MUMmer confirms 86% of sampled windows outright and places most of the rest at 98\(en99%.
.sp 0.3v
Among complete genomes, with both sides complete so that contamination cannot mimic transfer, 22% share
\(>= 1 kb at \(>= 99% identity with another genus, 14% with another family and 5.6% with another phylum, linking
1,116 genera in 6,528 genus pairs. The shared blocks carry metal resistance (odds ratio 6.4 against background
windows of the same genomes), transposases (4.6), integrases (4.1), antimicrobial resistance (2.6) and conjugation
genes (1.9), and few housekeeping genes (ribosome 0.21). Resistance genes make up 8.2% of genes on blocks shared
across phyla. At a fixed taxonomic distance, genomes from the same habitat are twice as likely to share DNA
(Mantel\(enHaenszel odds ratio 2.0, 95% CI 1.5\(en2.5), and the effect is concentrated in one connected group of
habitats: gut, other host sites, food and engineered systems such as wastewater. Between orders, gut genomes share
DNA 14 times more often than expected from taxonomic distance, engineered-system genomes 20 times, and food with
engineered systems 20 times, while soil, plant and freshwater genomes exchange 5 to 17 times less than expected.
Recent exchange between distant lineages is, in the sequenced record, mainly a feature of human-associated and
human-made habitats.
.sp 0.3v
The full atlas covers all 199,923 representatives (655 Gbp, 141 million windows) in 8 hours. 80% of them are
metagenome-assembled, so a shared block counts only if its contig extends beyond it and both genomes carry it. Then
23% of all representatives share \(>= 1 kb at \(>= 99% with another genus and 7.5% with another phylum; reduced-genome
symbionts (Patescibacteriota, Nanobdellota) share least. Among 14,881 habitat-labelled isolates the habitat effect is
stronger than among complete genomes (odds ratio 4.8, 95% CI 4.1\(en5.6), with the same human-associated cluster (gut
genomes 28 times the expected exchange between orders), and it survives a filter on the content of the shared DNA.
Among metagenome-assembled genomes the same-habitat excess grows with taxonomic distance. Tracing 59,700 MAGs to
their source samples shows why: MAGs from one sample share DNA at rates that barely fall with distance, and MAGs from
other studies of the same habitat share DNA with three times the housekeeping genes of isolate transfers, as
fragments of co-occurring organisms would. MAG pairs from different habitats look like transfer. Transfer ecology
from MAGs needs a check of what is shared, not only of where the genomes came from.
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
.H1 "Introduction"
.LP
Bacteria and archaea acquire genes from distant relatives: resistance to antibiotics and metals, new metabolism,
mobile elements that spread themselves. Transfer leaves a trace that does not depend on any gene model: a stretch of
DNA that two genomes share at far higher identity than their divergence allows. Smillie et al. used exactly this
trace, BLAST hits of \(>= 99% identity over \(>= 500 bp between genomes with less than 97% 16S identity, and found
among 2,235 genomes that shared ecology predicts recent exchange better than geography or phylogeny, with
human-associated bacteria exchanging 25 times more than others [1].
.LP
That study is 15 years old and its scale has not been repeated, because the search is all against all. Every
genome must be aligned locally against every other, and the cost grows with the square of the collection. GTDB R232
has 199,923 species representatives (655 Gbp). Whole-genome sketches (Mash, skani) give one identity per genome
pair and miss a 5 kb island in a 5 Mbp genome. Aligners find islands but cannot hold or scan the collection: our
earlier benchmark projected 236 h for minimap2 to map a 706 Gbp catalogue against GTDB.
.LP
aif answers the local question directly. A 5 kb window of a genome is a short contig, and aif's database search
returns its identity (ANI) and aligned fraction (AF) against every reference that shares sequence with it, in
milliseconds and from a 1 bit per base sketch. This report uses that to build a recent-transfer atlas:
.IP (1) 0.3i
\f[B]Method and validation\f[R] (Sections 2\(en3): windows, calls, blocks and genome pairs, checked against BLAST,
MUMmer and minimap2.
.IP (2) 0.3i
\f[B]A complete-genome atlas\f[R] (Sections 4\(en6): all 10,296 complete isolate genomes among the
representatives, searched against all 199,923 representatives and analysed among themselves, where contamination
cannot pose as transfer. Extent, what moves, and ecology against phylogeny.
.IP (3) 0.3i
\f[B]The full atlas\f[R] (Section 7): all 199,923 representatives, 80% of them metagenome-assembled, with a test
against contamination calibrated on the complete genomes, and the ecology retested on 14,881 isolates.
.H1 "Data and method"
.H2 "Genomes and habitats"
.LP
GTDB R232 lists 199,923 species representatives: 38,398 isolates (10,296 of them complete genomes or
chromosome-level assemblies), 160,761 MAGs and 764 single-cell genomes. Taxonomy is GTDB's. Habitats come from
NCBI's free-text isolation source, mapped by keyword rules (\f[CR]bench_real/hgt_atlas.py\f[R]) to 11 classes:
gut (any host), other host sites (oral, skin, respiratory, urogenital, blood, clinical), food, engineered systems
(wastewater, sludge, bioreactors, digesters, compost), marine, sediment, freshwater, subsurface, extreme (hot
springs, hypersaline, acid mine drainage, permafrost), soil and plant. 4,014 of the complete genomes receive a
habitat; the others have no usable isolation source (5,067) or one outside these classes (1,215).
.H2 "Windows, search and calls"
.LP
Every replicon is cut into consecutive 5 kb windows (a final piece of 1\(en5 kb is kept). Each window is searched
with \f[CR]aif db search\f[R] against the GTDB R232 database of all 199,923 representatives (0.88 bits per base,
$k = 20$), evaluating up to 128 candidate references per window. A \f[I]call\f[R] is a window aligned to the
representative of another species at ANI \(>= 99% over \(>= 500 aligned bases, the criterion of Smillie et al.
.IP \(bu 0.2i
\f[B]Blocks.\f[R] Consecutive called windows of one replicon with the same partner form a block.
.IP \(bu 0.2i
\f[B]Pairs.\f[R] Two genomes exchange when one of them has a block of \(>= 1 kb aligned to the other. Pairs whose
calls cover \(>= 30% of a genome are near-identical genomes or label conflicts, not transfer, and are set aside
(4 cases).
.IP \(bu 0.2i
\f[B]Distance.\f[R] The rank of the pair's lowest common GTDB taxon: same genus, same family (different genera),
and so on up to different domains.
.IP \(bu 0.2i
\f[B]Contamination.\f[R] A contaminated assembly carries a whole foreign contig, which would look like transfer
from the contaminant's lineage. Restricting both genomes to complete assemblies removes this: in a complete genome,
shared DNA is part of a chromosome or plasmid. Draft and metagenome-assembled partners are analysed only in the full
atlas, with a flank test and reciprocal anchoring (Section 7.1). Unrestricted partners show why: of 223 bacteria\(enarchaea "pairs"
with any partner, 162 involved a MAG, and the examples we examined were contaminant contigs (an archaeal
\f[I]Methanothrix\f[R] block of 23 kb at 99.99% inside a bacterial MAG; a \f[I]Cupriavidus\f[R] window, a common
reagent contaminant, in a draft archaeal assembly).
.H2 "Ecology against phylogeny"
.LP
Among the 4,014 habitat-labelled complete genomes (8.05 million unordered pairs), pairs are stratified by the rank
of their lowest common taxon, and in each stratum the exchange rate of same-habitat pairs is compared with that of
different-habitat pairs. The Mantel\(enHaenszel odds ratio summarises the strata. Pairs share genomes, so intervals
come from a bootstrap over genomes (200 replicates). For each pair of habitats, the number of exchanging pairs
is compared with the number expected from the stratum rates, which removes the effect of the habitats' different
taxonomic composition (\f[CR]bench_real/hgt_ecology.py\f[R]).
.H2 "Genes on shared blocks"
.LP
Blocks shared with another genus are merged per replicon into 9,492 regions (101 Mbp). Genes are called with
prodigal and assigned KEGG orthologs by DIAMOND against CheckM2's UniRef100 database (best hit, e-value
\(<= 10\u\(mi10\d, \(>= 30% identity over \(>= 50% of the protein). As background, random 5 kb windows of the same
genomes outside the regions are matched in total length per genome (111,632 genes on regions, 114,610 in background).
Genes are grouped into classes by the KEGG hierarchy and names (\f[CR]bench_real/hgt_annotate.py\f[R]).
.H1 "Validation and cost"
.LP
Three questions: are the calls real, how many are missed, and what would the alternatives cost.
.IP \(bu 0.2i
\f[B]Precision.\f[R] 400 calls between genera were drawn at random and each window aligned to its partner genome
with MUMmer and with BLAST. MUMmer confirms 342 (85.5%) at the strict criterion (\(>= 500 bases at \(>= 99%), BLAST
345 (86.3%). Of the rest, 39 align at 98\(en99%, below the threshold but clearly shared, and one window does not
align. aif's identity differs from MUMmer's by 0.23 pp on average.
.IP \(bu 0.2i
\f[B]Recall.\f[R] BLAST (megablast) was run for every window of 60 random complete genomes against a universe of
2,000 complete genomes, giving 629 (window, partner) pairs with an HSP of \(>= 99% over \(>= 500 bp, 227 of them
between genera. Table 1 scores aif and minimap2 against it. At the level that the analyses use, genome pairs
sharing \(>= 1 kb between genera, aif with 128 candidates per window finds 76% of BLAST's pairs at 98% precision.
At the window level it finds 52% at 94% precision.
.IP \(bu 0.2i
\f[B]What is missed.\f[R] Of BLAST's cross-genus windows, aif with 256 candidates evaluated 91%; most misses are windows
where a short stretch at \(>= 99% lies within a longer alignment at lower identity, for which aif reports the
identity of the whole aligned part (often 96\(en98%). These mosaic windows are a definitional difference rather
than missed homology: BLAST itself merges such stretches into a single lower-identity HSP when they are adjacent.
.IP \(bu 0.2i
\f[B]minimap2\f[R] recovers 87% of BLAST's windows but its extra calls are mostly not confirmed: BLAST confirms 16
of 60 sampled minimap2-only calls (and 14 of 60 aif-only calls). Its assembly preset cuts alignments at divergent
stretches and reports the identical pieces, which inflates "\(>= 99%" calls. Its index of the 2,000 genomes took
51 GB, so GTDB would have to be searched in about 80 parts.
.IP \(bu 0.2i
\f[B]Design.\f[R] Smaller windows raise recall but lose precision (1 kb: 53% and 47% between genera), because a 1 kb
window's identity is less certain. More candidates raise recall at a cost in time; 128 is the compromise used here.
Counting only aligned mass at levels \(>= 99% (\f[CR]--af-min-ident 0.99\f[R]) adds many unconfirmed calls.
.ne 2.4i
.TC 1
Validation against BLAST (HSP \(>= 99% over \(>= 500 bp). Windows: 52,055 windows of 60 complete genomes against
2,000 complete genomes (227 cross-genus pairs). Genome pairs: pairs sharing \(>= 1 kb between genera (55). Time:
the 300-genome validation set (1.3 Gbp) against all 199,923 representatives for aif, against the 2,000 genomes for
minimap2, both on 48 threads.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_valid.tbl
.ps \n(PS
.vs \n(VSp
.LP
\f[B]Cost.\f[R] The complete-genome search (44.5 Gbp, 8.9 million windows, against all 199,923 representatives)
took 2,346 s on 48 threads, with a peak of 94.5 GB, almost all of it the memory-mapped database. minimap2 took
390 s to map 1.3 Gbp of windows against 8.6 Gbp; by linear scaling, the same search against all representatives
would take about 12 days, in line with the 236 h projected for minimap2 at GTDB scale in the gut-catalogue report.
BLAST took 316 s for 0.26 Gbp against 8.6 Gbp on 48 processes, which scales to about 48 days. aif's search costs
the same per window whatever the size of the collection behind the screen.
.H1 "How much DNA moves between genera"
.LP
Table 2 counts exchanging genome pairs by their lowest common taxon. Most sharing is between species of one genus
(36,264 directed pairs), but distant sharing is common: 7,726 pairs between genera of one family, 7,633 between
orders of one class, and 3,332 between phyla.
.IP \(bu 0.2i
\f[B]Genomes.\f[R] 2,281 of the 10,296 complete genomes (22.2%) share \(>= 1 kb at \(>= 99% with another genus,
1,454 (14.1%) with another family and 575 (5.6%) with another phylum.
.IP \(bu 0.2i
\f[B]Network.\f[R] These pairs link 1,116 genera in 6,528 genus pairs. The best-connected genera are
\f[I]Enterobacter\f[R] (155 partner genera), \f[I]Aquipseudomonas\f[R] (137), \f[I]Providencia\f[R] (132),
\f[I]Proteus\f[R] (115), \f[I]Paracoccus\f[R] (114), \f[I]Streptococcus\f[R] (111), \f[I]Acinetobacter\f[R]
(110), \f[I]Vibrio\f[R] (107), \f[I]Corynebacterium\f[R] and \f[I]Klebsiella\f[R] (105 each).
.IP \(bu 0.2i
\f[B]Blocks.\f[R] 45,884 blocks between genera, 190 Mbp in all. 32% lie on plasmids or other secondary replicons,
the rest on chromosomes. 15% are identical (\(>= 99.95%), so transfer is ongoing. Half are a single window
(\(<= 5 kb); 10% span \(>= 10 kb and 1% \(>= 35 kb, the size of integrative and conjugative elements.
.IP \(bu 0.2i
\f[B]Across phyla.\f[R] The most connected phylum pairs are Actinomycetota\(enPseudomonadota (1,500 genome pairs),
Bacillota\(enBacillota_I (387), Actinomycetota\(enBacillota (312) and Bacteroidota\(enPseudomonadota (250).
.IP \(bu 0.2i
\f[B]Across domains.\f[R] 19 pairs of complete genomes share DNA between bacteria and archaea, and 18 of them
involve one \f[I]Methanobacterium\f[R] genome (GCA_000499765.1). Its shared block, 1.0\(en2.5 kb at 99.1\(en99.9%,
carries a transposase (KEGG K07484) next to a site-specific recombinase and is found in ten genera of thermophilic
or anaerobic bacteria (among them \f[I]Fervidobacterium\f[R], \f[I]Petrotoga\f[R], \f[I]Defluviitoga\f[R],
\f[I]Thermoanaerobacterium\f[R], \f[I]Mahella\f[R] and \f[I]Dehalobacter\f[R]) from hot springs, oil fields and
digesters: an insertion sequence moving among organisms that share those habitats, across the domain boundary.
.LP
The candidate cap matters for distant partners: with 32 candidates per window instead of 128, the same search
finds 2,146 genomes with cross-genus sharing (\(mi6%) but only 1,257 cross-phylum pairs (\(mi62%), because close
relatives of the query fill the slots. The counts here are therefore lower bounds, most of all for the most
distant and most promiscuous elements.
.ne 2.1i
.TC 2
Exchanging genome pairs (directed: a block of \(>= 1 kb at \(>= 99% in the first genome) by lowest common taxon.
\(lqAny partner\(rq includes draft and metagenome-assembled partners, which contamination can mimic (Section 2.2).
.TCE
.ps 9
.vs 11p
.so gen/tab_dist.tbl
.ps \n(PS
.vs \n(VSp
.H1 "What moves"
.LP
Genes on shared blocks differ sharply from genes on background windows of the same genomes (Table 3, Figure 1).
.IP \(bu 0.2i
\f[B]Enriched:\f[R] metal resistance (3.2% of genes against 0.5%; odds ratio 6.4: mercury \f[I]mer\f[R] operons,
\f[I]czc\f[R] cobalt\(enzinc\(encadmium efflux, copper and zinc ATPases), transposases (4.6), integrases and
recombinases (4.1), antimicrobial resistance (2.6) and conjugation machinery (1.9).
.IP \(bu 0.2i
\f[B]Depleted:\f[R] ribosomal and translation genes (0.21), transporters (0.59), and KEGG-annotated genes in
general (0.43); genes without any KEGG ortholog are enriched (1.6), as expected of mobile elements.
.LP
Individual orthologs make the point (Table 4): the IS26 transposase occurs 586 times on shared blocks and 21 times
in background, mercuric reductase \f[I]merA\f[R] 285 against 10, and the tetracycline resistance gene
\f[I]tetM\f[R] 225 against 4, the signature of Tn916-type conjugative transposons that cross phyla. Resistance
genes become more prominent with distance: they are 8.2% of genes on blocks shared between phyla, against 2.8%
overall. By habitat (Table 5), resistance genes make up 4.4\(en4.7% of the shared genes of gut, other host and food
genomes, against 1.4\(en1.8% for freshwater, soil and plant genomes and 0.2% for extreme habitats, while metal
resistance is highest in marine and subsurface genomes.
.ne 3.0i
.PSPIC -C figs/fig_func.eps 4.6i
.FC 1
Gene classes on blocks shared with another genus against background windows of the same genomes (odds ratio,
log scale). Complete genomes, 111,632 and 114,610 genes.
.FE
.ne 2.9i
.TC 3
Gene classes on shared blocks and on background windows (genes with a KEGG ortholog assigned to the first matching
class; transposases first).
.TCE
.ps 9
.vs 11p
.so gen/tab_func.tbl
.ps \n(PS
.vs \n(VSp
.ne 3.2i
.TC 4
The most frequent KEGG orthologs on shared blocks, with their count in the background windows.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_kos.tbl
.ps \n(PS
.vs \n(VSp
.ne 2.4i
.TC 5
Genes on shared blocks by the habitat of the genome that carries them: resistance (with its share in that
habitat's background windows), metal resistance, transposases and conjugation genes.
.TCE
.ps 9
.vs 11p
.so gen/tab_habfunc.tbl
.ps \n(PS
.vs \n(VSp
.H1 "Ecology against phylogeny"
.LP
Exchange falls steeply with taxonomic distance, 5 to 10 times per rank from genus to order (Figure 2). At every
distance, pairs from the same habitat exchange more (Table 6): the ratio of rates is 1.8 for pairs within a genus or
a family, and 2.2\(en2.3 for pairs that share only an order, a class or the domain, with bootstrap intervals that
exclude 1. The exception is pairs that share only a phylum, 1.25 (0.48\(en2.72). Summarised over the strata above the genus,
the Mantel\(enHaenszel odds ratio is 2.0 (1.5\(en2.5); between families, 2.1 (1.5\(en2.8). Shared habitat
predicts recent exchange independently of relatedness, as Smillie et al. found on 2,235 genomes.
.ne 3.0i
.PSPIC -C figs/fig_strata.eps 4.6i
.FC 2
Share of genome pairs that exchange (\(>= 1 kb at \(>= 99%), by the rank of their lowest common taxon, for pairs
from the same habitat and from different habitats (4,014 habitat-labelled complete genomes).
.FE
.ne 2.2i
.TC 6
Exchange by taxonomic distance and habitat (unordered pairs of habitat-labelled complete genomes). Intervals:
bootstrap over genomes.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_strata.tbl
.ps \n(PS
.vs \n(VSp
.LP
The average hides where the effect lies. Figure 3 compares, for every pair of habitats, the number of
exchanging pairs between different orders with the number expected from their taxonomic distances:
.IP \(bu 0.2i
\f[B]A human-associated and engineered cluster.\f[R] Within habitats, gut genomes exchange 14 times more than
expected, engineered-system genomes 20 times and other host sites 5.5 times. Between habitats, food with engineered
systems 20 times, gut with food 8.7, host sites with engineered systems 7.9, gut with engineered systems 7.7 and gut
with other host sites 4.4. Wastewater and other engineered systems connect the human body, food and the environment:
they are enriched with every member of the cluster, and with sediment (4.5 times between genera).
.IP \(bu 0.2i
\f[B]Natural environments.\f[R] Soil, plant and freshwater genomes exchange 5 to 17 times less than expected
among themselves (0.16, 0.20, 0.06), and marine genomes about as expected (1.1). Pairs across natural habitats are
also depleted (marine with soil 0.17).
.LP
The pattern does not come from shared plasmids within one family: it is stronger between orders than between genera
(gut 14 against 8.9 times; Table 7). It resembles Smillie et al.'s finding of
25-fold more exchange among human-associated bacteria, and adds the engineered systems that connect them.
.ne 4.4i
.PSPIC -C figs/fig_habitats.eps 4.9i
.FC 3
Exchange between different orders for each pair of habitats: observed exchanging pairs (numbers) over the number
expected from the pairs' taxonomic distances (colour, log scale). Blank: fewer than 2 expected.
.FE
.ne 1.6i
.TC 7
Exchange within the human-associated and engineered habitats, observed over expected, as more closely related pairs
are excluded.
.TCE
.ps 9
.vs 11p
.TS
center box tab(@);
l c c c.
habitats@between genera@between families@between orders
_
gut@8.9@11.8@14.2
engineered@20.5@22.8@20.4
other host sites@3.5@5.0@5.5
food \(em engineered@19.4@23.1@20.4
gut \(em food@5.6@8.6@8.7
gut \(em engineered@5.6@6.9@7.7
soil@0.50@0.18@0.16
freshwater@0.29@0.04@0.06
.TE
.ps \n(PS
.vs \n(VSp
.H2 "Caveats"
.IP \(bu 0.2i
\f[B]What is sequenced.\f[R] Complete genomes over-represent pathogens and model organisms. The 1,016
human-associated and engineered complete genomes are taxonomically varied (largest families Mycobacteriaceae 87,
Enterobacteriaceae 41, Campylobacteraceae 38, Streptococcaceae 32, Lactobacillaceae 31), and the comparison within
taxonomic strata removes the effect of composition on the expected counts. It does not remove the possibility that
sequenced human-associated isolates were selected for mobile elements (for example, resistant clinical isolates). Natural environments are also more diverse, so a transfer partner is less likely to have been sequenced.
.IP \(bu 0.2i
\f[B]Habitats.\f[R] Labels come from free text. Half the complete genomes have none, and the classes are broad.
.IP \(bu 0.2i
\f[B]Sensitivity.\f[R] About a quarter of BLAST's exchanging pairs are missed (Section 3), mostly short or mosaic
shared stretches, and widely shared elements can exceed 128 candidates per window. Counts are lower bounds; the
comparisons between habitats assume that misses do not depend on habitat.
.IP \(bu 0.2i
\f[B]Scope.\f[R] One genome per species, so exchange within species is not measured; \(>= 99% identity means recent
but is not dated; the direction of transfer is not inferred.
.H1 "The full atlas"
.LP
The same search was run for every representative: 199,923 genomes, 655 Gbp, 141.3 million windows. It took 28,764 s
(8.0 h) on 48 threads in 128 parts, with a peak of 112 GB, again almost all of it the memory-mapped database. 7.6
million windows (5.4%) align to another species at \(>= 97%. 80% of the representatives are metagenome-assembled
(MAGs), which adds uncultured lineages and 71,605 habitat-labelled genomes, but also contamination.
.H2 "Contamination and mis-binning"
.LP
A contig placed in the wrong bin carries its organism's DNA into another genome, where it looks like transfer. To see
how often, the 5 kb windows on either side of every cross-genus block (its flanks) were searched against all
representatives at \(>= 80% identity (1.38 million windows, 6.7 Gbp, 5 minutes of search;
\f[CR]bench_real/hgt_flanks.py\f[R]). The closest relative of a flank shows whose genome the contig belongs to.
Table 8 compares the genome categories:
.IP \(bu 0.2i
\f[B]Complete genomes\f[R] cannot carry a foreign contig, so they show what real transfer looks like. Almost every
block (99.3%) has a flank, and the closest relative of the flanks is on the host's side of the tree for 52% of
blocks. For 16% it is on the partner's side: the moved element continues beyond its \(>= 99% stretch at lower
identity.
.IP \(bu 0.2i
\f[B]Draft isolates and MAGs.\f[R] 30% and 46% of blocks are entire contigs. Assemblies break at repeats, so
plasmids, phages and insertion sequences often form contigs of their own, and so do contaminant contigs. The two
cannot be told apart, and such blocks are dropped. Among blocks with flanks, the partner's side is closest for 19%
of draft isolates' blocks and 24% of MAGs', against 16% in complete genomes: an excess of 3 to 8 points that may be
contamination.
.LP
The full atlas therefore counts a block only if its contig extends beyond it by at least one window, and a pair of
genomes only if \f[I]both\f[R] carry such a block with the other, so that contamination of either genome alone
cannot create a pair. The second condition costs recall: 68% of the cross-genus pairs between complete genomes
are found from both sides. As a sensitivity test, blocks with a flank on the partner's side were also rejected in
every genome, complete ones included. This removes 37% of the cross-genus pairs; the number of genomes with
cross-genus sharing falls by 25% among complete genomes and by 38% among MAGs.
.ne 1.9i
.TC 8
Closest relative of the flanks of cross-genus blocks (\(>= 1 kb), by the category of the genome that carries the
block (% of blocks; a block with flanks on both sides takes the flank nearer the partner).
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_flank.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Extent"
.LP
46,467 representatives (23.2%) share \(>= 1 kb at \(>= 99% identity with another genus, 33,135 (16.6%) with another
family and 14,909 (7.5%) with another phylum. They form 417,376 directed cross-genus pairs that link 12,915 genera in
77,228 genus pairs, through 612,754 blocks (3.55 Gbp counted once per partner). 17% of the blocks are identical
(\(>= 99.95%), and in complete genomes 27% lie on secondary replicons.
.IP \(bu 0.2i
\f[B]By genome category\f[R] (Table 9): complete isolates 30.6%, draft isolates 32.0%, MAGs 21.3%, single-cell
genomes 9.9%. Binning loses plasmids and other mobile DNA, whose coverage and composition differ from the
chromosome, and single-cell genomes are incomplete. The MAG figures are therefore lower bounds, but they still add
34,250 genomes with recent cross-genus sharing.
.IP \(bu 0.2i
\f[B]By distance\f[R] (Table 10): as in the complete genomes, sharing falls with distance, with a rise at the class
level (101,076 pairs). 802 directed pairs cross the domain boundary, 668 of them between MAGs; the partner test
leaves 262.
.IP \(bu 0.2i
\f[B]By phylum\f[R] (Table 11): Bacillota share most (37% of genomes with another genus), then Bacillota_I,
Acidobacteriota and Pseudomonadota (27\(en28%). Lineages with small genomes, many of them symbionts, share least:
Patescibacteriota (the
candidate phyla radiation) 6.6%, Nanobdellota (DPANN archaea) 4.9%, Micrarchaeota 6.1%, Chlamydiota 6.4% and
Omnitrophota 8.0%. Archaea share less than bacteria overall (10.4% against 23.9%).
.IP \(bu 0.2i
\f[B]Hubs.\f[R] The best-connected genera now include many gut Bacillota:
\f[I]Blautia_A\f[R] (442 partner genera), \f[I]Oscillibacter\f[R] (396), \f[I]Lawsonibacter\f[R] (389),
\f[I]Roseburia\f[R] (357), \f[I]Enterocloster\f[R] (349) and \f[I]Faecalibacterium\f[R] (333), together with
\f[I]Streptococcus\f[R] (485), \f[I]Aquipseudomonas\f[R] (447), \f[I]Collinsella\f[R] (408) and
\f[I]Paracoccus\f[R] (400). The most connected phyla are Bacillota and Bacillota_I (10,122 genome pairs),
Actinomycetota and Bacillota (6,494), and Actinomycetota and Pseudomonadota (5,948).
.ne 1.7i
.TC 9
Genomes sharing \(>= 1 kb at \(>= 99% with another genus or phylum, by category. Main rule: block on a contig that
extends beyond it, found from both genomes. Partner test: also no flank nearer the partner's lineage, in any genome.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_cat.tbl
.ps \n(PS
.vs \n(VSp
.ne 1.9i
.TC 10
Exchanging genome pairs (directed) by lowest common taxon in the full atlas. One direction: block of \(>= 1 kb in
the first genome, without contamination control; both directions: the main rule.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_dist.tbl
.ps \n(PS
.vs \n(VSp
.ne 3.0i
.TC 11
Genomes sharing \(>= 1 kb at \(>= 99% with another genus or phylum, by phylum (the 16 largest phyla and the largest
phyla with less than 1% isolates).
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_phyla.tbl
.ps \n(PS
.vs \n(VSp
.H2 "What moves"
.LP
Blocks shared with another genus were merged per contig into 151,196 regions (1.41 Gbp in 46,467 genomes), with
1.50 million genes, against 1.54 million genes in background windows of the same genomes. The enrichments of the
complete-genome atlas recur (Table 12): metal resistance (odds ratio 5.4), conjugation machinery (3.8), integrases
(3.6), transposases (3.6) and resistance genes (1.9), while ribosomal and translation genes are depleted (0.36). Most
odds ratios are somewhat smaller than among complete genomes; conjugation is the exception (3.8 against 1.9), and the
type IV coupling protein VirD4 is among the most frequent genes on shared blocks (4,461 against 901 in background).
Other frequent genes are transposases, the mercury-resistance regulator \f[I]merR\f[R] (3,146 against 224), zinc,
cadmium and copper ATPases (\f[I]zntA\f[R], \f[I]copA\f[R]) and \f[I]czc\f[R] efflux. Resistance genes again
become more common with distance, from 1.0% of genes on blocks shared within a family to 2.6% between phyla.
.ne 2.9i
.TC 12
Gene classes on blocks shared with another genus against background windows of the same genomes, full atlas, with
the odds ratios of the complete-genome atlas (Table 3).
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_func.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Ecology: isolates and MAGs"
.LP
The full atlas has 14,881 habitat-labelled isolates (complete and draft), 3.7 times the complete genomes of
Section 6 and 6.7 times the genomes of Smillie et al., and 71,605 habitat-labelled MAGs and single-cell genomes. The
two give different answers (Table 13).
.IP \(bu 0.2i
\f[B]Isolates\f[R] confirm Section 6 more strongly. Same-habitat pairs exchange 2.6 to 9.6 times more than
different-habitat pairs at every distance, and the Mantel\(enHaenszel odds ratio between genera is 4.8 (4.1\(en5.6).
The habitat matrix (Figure 4, Table 14) shows the same cluster: between orders, gut genomes exchange 28 times more
than expected, other host sites and engineered systems 13 times each, gut with other host sites 7.2 times and food
with engineered systems 6.4 times, while soil (0.16), freshwater (0.11) and plant genomes (0.05) exchange far less
than expected. Leaving out pairs from the same sequencing project changes little (odds ratio 4.6, 3.9\(en5.3).
.IP \(bu 0.2i
\f[B]MAGs\f[R] show a same-habitat excess that grows with distance beyond the family: 2.7 times within families,
13 times within phyla, 43 times between phyla and 80 times between domains. Natural habitats now look enriched as
well (between orders, soil 13 times, sediment 19 and extreme habitats 12). Leaving out pairs from the same project
removes 23% of the exchanging pairs and lowers these 3- to 4-fold (soil 3.5, sediment 4.9, extreme 3.2), and the
excess between phyla to 31 times and between domains to 46, but the rise with distance remains.
.LP
.ne 1.9i
.TC 13
Exchange by habitat in each subset (MAGs include single-cell genomes): ratio of the exchange rates of same- and different-habitat pairs at each
lowest common taxon, and the Mantel\(enHaenszel odds ratio over the strata above the genus (95% bootstrap interval).
\(lqOther projects\(rq leaves out pairs from the same NCBI BioProject.
.TCE
.ps 7.5
.vs 9.5p
.so gen/tab_full_eco.tbl
.ps \n(PS
.vs \n(VSp
.LP
A same-habitat excess that is largest where transfer is rarest is not what transfer produces; Section 7.5 traces it
to its sources. Section 6 rests on isolates, and the larger isolate set confirms it. Isolates and MAGs are not
pooled: gut genomes are 93% MAGs, whose lower sharing rates (Table 9) would pass for a habitat effect.
.br
.ne 2.6i
.TC 14
Exchange between different orders within and between habitats, observed over expected, in each subset.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_hab.tbl
.ps \n(PS
.vs \n(VSp
.br
.ne 4.4i
.PSPIC -C figs/fig_full_habitats.eps 4.9i
.FC 4
Exchange between different orders for each pair of habitats among the 14,881 habitat-labelled isolates of the full
atlas: observed exchanging pairs (numbers) over the number expected from the pairs' taxonomic distances (colour,
log scale). Blank: fewer than 2 expected.
.FE
.H2 "What MAGs share"
.LP
Two artifacts could make MAGs of one habitat share DNA. MAGs binned from the same sample can share a contig
(cross-binning, contamination), and different projects reassemble the same public metagenomes, which BioProjects do
not reveal. And a MAG contig can carry a fragment of another organism of its community, joined to it by chimeric
assembly or binned with it, while that organism's representative comes from another study of the same habitat. The
first needs a shared sample; the second only a shared habitat.
.LP
\f[B]Sample provenance.\f[R] Every MAG has its own BioSample, whose record names the source as the metagenome's
BioSamples or as SRA runs (\f[CR]hgt_provenance.py\f[R]; SRA accessions are resolved through SRA run tables). Of the
161,525 MAG and single-cell representatives, 147,569 records name a source, and 59,700 of the 71,605 habitat-labelled
MAGs (83%) resolve to one or more of 27,503 source samples (9,173 MAGs to ten or more: co-assemblies). The 2.2
million pairs of these MAGs with a common source sample exchange far more than other pairs, and their rate barely
falls with distance: 43-fold from the family to between domains, against 2,500-fold for other pairs, so their excess
grows from 7 to 382 times (Table 15). Almost all of them are also from one project, however. Across the atlas, 0.1 to
1.2% of exchanging MAG pairs share a sample but not a project, against 19 to 57% that share a project. Leaving out
co-sample pairs on top of same-project pairs changes nothing (odds ratio 4.2 either way, Table 16), and the rise with
distance remains (8 times within phyla, 30 times between phyla).
.br
.ne 1.9i
.TC 15
Exchange of MAG pairs with a common source sample, against all other pairs of the same lowest common taxon, among
the 59,700 habitat-labelled MAGs with a known source.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_cosample.tbl
.ps \n(PS
.vs \n(VSp
.LP
\f[B]Gene content.\f[R] Transferred DNA is rich in mobile-element genes and poor in housekeeping genes; a fragment
of another genome carries that genome's genes at random. Table 17 sorts the genes of the shared regions (Section
7.3) by the kind of pair that shares them (\f[CR]hgt_magcheck.py\f[R]). Isolate pairs, and MAG pairs from other
studies and different habitats, carry the same signature of transfer: 7% mobile-element genes and 0.6\(en0.7%
ribosomal or translation genes, against 3.5\(en4.0% in background windows. MAG pairs from one study, and MAG pairs
from other studies of the same habitat, carry half the mobile-element genes and about three times the ribosomal
genes. With isolate transfers and background windows as the two references, the ribosomal genes say that 36\(en45%
of their shared genes come from random fragments of genomes. The excess holds in each habitat: between families, gut MAG pairs
carry 1.75% ribosomal genes against 0.61% for gut isolate pairs, marine MAG pairs 2.6% against 0.4%. The second
artifact fits: organisms of different habitats cannot meet in one assembly, and their shared DNA looks like transfer.
Genuine exchange within a community also needs co-occurrence, but it would carry mobile-element genes.
.br
.ne 2.4i
.TC 16
Exchange by habitat among MAGs under provenance and content controls, as in Table 13. \(lqKnown source\(rq:
MAGs whose source samples are known; \(lqother samples\(rq and \(lqother projects\(rq leave out pairs with a common
source sample or BioProject. Mobile gene: only exchanging pairs whose shared regions carry a mobile-element
gene (no same-genus pairs: their regions were not annotated); \(lqall controls\(rq: known source, other projects and
samples.
.TCE
.ps 7
.vs 9p
.so gen/tab_full_prov.tbl
.ps \n(PS
.vs \n(VSp
.LP
\f[B]A content filter.\f[R] Keeping only the pairs whose shared DNA carries a mobile-element gene keeps 59% of
isolate pairs and 22\(en28% of the MAG pairs rich in fragments (Table 17). Among isolates the habitat effect
survives (odds ratio 3.5, 3.0\(en4.3; between orders gut 21 times the expected exchange, soil 0.18). Among MAGs with
all controls it is 2.7 (2.1\(en3.6) on 1,918 pairs, but it still grows with distance (22 times between phyla), because
the filter enriches transfer without removing fragments. MAG-based measures of transfer ecology are therefore upper
bounds until the junctions of shared blocks are checked in the source reads. Observed-over-expected ratios of MAGs
and isolates (Table 14) also cannot be compared directly: gut genomes are 35% of the MAGs and 11% of the isolates, so
gut pairs set much of the MAGs' expected rates and their ratio has little room to rise.
.br
.ne 2.2i
.TC 17
Gene content of the shared regions between genera by kind of genome pair, against background windows of the same
genomes. \(lqMAGs\(rq: pairs with at least one MAG; study: BioProject or source sample. Mobile: transposases,
integrases and recombinases, conjugation, phage. Resistance: antimicrobial and metal. Pairs: unordered pairs with
annotated shared genes.
.TCE
.ps 8.5
.vs 10.5p
.so gen/tab_full_content.tbl
.ps \n(PS
.vs \n(VSp
.H1 "Reproduction"
.LP
.EQ
delim off
.EN
.nf
.ft CR
.ps 8.5
.vs 10.5p
/usr/bin/python3 -I bench_real/hgt_atlas.py select                 # genome table, habitats
/usr/bin/python3 -I bench_real/hgt_atlas.py windows complete 64    # 5 kb windows
HGT_NCAND=128 /usr/bin/python3 -I bench_real/hgt_atlas.py search complete
HGT_PARTNERS=complete /usr/bin/python3 -I bench_real/hgt_analyze.py complete
/usr/bin/python3 -I bench_real/hgt_ecology.py complete_complete complete complete
/usr/bin/python3 -I bench_real/hgt_annotate.py complete_complete
/usr/bin/python3 -I bench_real/hgt_validate.py precision complete 400   # also: recall, variants,
/usr/bin/python3 -I bench_real/hgt_validate.py blast_truth complete     #   blast, adjudicate
# full atlas (Section 7); P="/usr/bin/python3 -I", T=all_flankcontig_reciprocal
HGT_GZ=1 $P bench_real/hgt_atlas.py windows all 128
HGT_NCAND=128 HGT_DELETE=1 $P bench_real/hgt_atlas.py search all
$P bench_real/hgt_analyze.py all                     # one direction, no control
$P bench_real/hgt_flanks.py extract all all          # then: search all; calibrate all all
HGT_FLANKS=1 HGT_FLANK_MODE=contig HGT_RECIPROCAL=1 $P bench_real/hgt_analyze.py all
HGT_FLANKS=1 HGT_FLANK_MODE=uniform HGT_RECIPROCAL=1 $P bench_real/hgt_analyze.py all
$P bench_real/hgt_ecology.py $T all isolate          # and: mag
HGT_ECO_XSTUDY=1 $P bench_real/hgt_ecology.py $T all isolate   # and: mag
$P bench_real/hgt_annotate.py $T
$P bench_real/hgt_provenance.py fetch    # then: parse, resolve, resolve_samples, sources
HGT_ECO_KNOWN=1 HGT_ECO_XSTUDY=1 HGT_ECO_XSAMPLE=1 $P bench_real/hgt_ecology.py $T all mag
$P bench_real/hgt_magcheck.py $T
$P bench_real/hgt_ecology.py ${T}_mobile all isolate   # and: mag
/usr/bin/python3 -I scripts/hgt_report.py; report_hgt/build.sh
.ps \n(PS
.vs \n(VSp
.ft R
.fi
.EQ
delim $$
.EN
.LP
Results are in \f[CR]results_hgt/\f[R]; windows, hits and validation files in \f[CR]data/hgt/\f[R]
(git-ignored).
.HU "References"
.RF 1
C. S. Smillie, M. B. Smith, J. Friedman, O. X. Cordero, L. A. David, E. J. Alm. Ecology drives a global network of
gene exchange connecting the human microbiome. \f[I]Nature\f[R] 480, 241\(en244, 2011.
.RE
