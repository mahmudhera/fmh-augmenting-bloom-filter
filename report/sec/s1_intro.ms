.H1 "Introduction and problem statement"
.LP
Average nucleotide identity (ANI) between two genomes is the mean identity of their
alignable regions. The usual 95% boundary between prokaryotic species is defined in
these terms [37, 38, 25].
Sketching methods replace alignment with comparisons between small samples of $k$-mers.
With FracMinHash (FMH) [4, 5, 6] a sequence is reduced to the $k$-mers whose hash
falls below a fraction $s$ of the hash range.
The containment of a query $A$ in a reference $B$, $C = | K sub A \(ca K sub B | / | K sub A |$,
is estimated from the two sketches.
Under a simple substitution model $Ex [C] = (1-p) sup k$, so ANI is estimated as
$ANIh = Chat sup {1/k}$ [6, 7].
This estimator works well for genome-sized inputs. Two things go wrong on short
contigs (a few kb), which are typical of metagenome-assembled genomes (MAGs), plasmids
and mobile elements:
.IP (i)
\f[B]Variance.\f[R] A 5 kb contig at $s = 1/1000$ contributes about five sampled $k$-mers.
Its ANI standard deviation is about 2.8 percentage points (pp) at 95% ANI
(Section 3), so the estimate cannot even separate species.
.IP (ii)
\f[B]Confounding with alignment fraction (AF).\f[R] If only a fraction $alpha$ of the contig is
homologous to $B$, containment estimates $alpha (1-p) sup k$ and the ANI estimate is biased by the
factor $alpha sup {1/k}$. For example, $alpha = 0.5$ turns 95% into 91.9%.
An unordered set comparison cannot undo this (Proposition 3).
.LP
\f[B]The proposal under study.\f[R]
Keep the FMH-sampled $k$-mers of the reference together with their \f[I]positions\f[R]
(we call them \f[I]anchors\f[R]). The stretches between consecutive anchors are
\f[I]islands\f[R]. Record the island $k$-mers in a Bloom filter (BF) [1, 13],
either all of them or a second FMH subsample with a larger threshold $s sub 2 > s$.
The questions are:
.IP Q1. 0.4i
Is this space-efficient?
.IP Q2. 0.4i
Is it linear (in time and in space)?
.IP Q3. 0.4i
What ANI estimator does it support, and by how much does it improve on plain FMH,
especially for short contigs?
.IP Q4. 0.4i
Can it estimate AF?
.IP Q5. 0.4i
Which nearby ideas are worth considering, and how do they compare?
.LP
\f[B]Summary of answers.\f[R]
.IP Q1: 0.4i
Not if the filter holds all island $k$-mers. Any approximate membership query (AMQ)
structure with useful false-positive rate (FPR) then costs more than the 2-bit-packed
genome (Theorem 2). With a second FMH threshold the cost becomes
$s sub a (b sub a + c sub pos ) + s sub 2 c( eps )$ bits per base,
and it can be set equal to the plain FMH budget (Table 6).
.IP Q2: 0.4i
Yes. Construction and querying are single streaming passes, linear in sequence length,
and space is linear in reference length with an explicit, tunable constant (Section 4.3).
.IP Q3: 0.4i
Probing every tier-2 query $k$-mer against the reference filter gives an unbiased
FP-corrected containment estimator whose variance is given exactly by Theorem 1.
At the same number of bits it reduces sampling variance by $R sup * = 64 / J sup *$,
between 4\(mu and 23\(mu depending on ANI (Theorem 3, Table 8).
The gain comes from spending about 1\(en10 bits per sampled $k$-mer instead of 64,
not from the Bloom filter as such.
It requires the query to enumerate its $k$-mers, which is the normal situation for
a contig queried against a sketched reference database.
.IP Q4: 0.4i
Yes. Positions make AF identifiable, and tier-2 hits localise alignment boundaries
with error of about $1/(s sub 2 theta )$ bases. Estimating AF of the query needs positions
only on the query side (Proposition 5); positions on the reference side add colinearity,
localisation and chimera detection.
.IP Q5: 0.4i
Fourteen designs are compared (Table 3). The most promising relatives are
island-keyed localisation, which enables small $k$ at low ANI;
island fingerprints, which win at strain level;
hedged multi-FPR tiers for unknown divergence;
and $b$-bit FMH for symmetric (sketch-vs-sketch) use.
.LP
\f[B]Scope.\f[R] This document is purely analytical; no sketching software has been implemented.
The only computation is (a) evaluation of the closed forms that produce the tables and figures,
and (b) one small Monte Carlo simulation of the idealised model,
used to check the central variance formula (Table 5).
All numerical examples use $k=21$, a 5 Mbp reference genome and, as the baseline,
FMH at $s = 1/1000$ with 64-bit hashes ($beta sub 0$ = 0.064 bits per reference base) unless stated otherwise.
.ne 3.2i
.sp 0.4v
.PSPIC -C figs/fig_schematic.eps 6.4i
.FC 1
The anchored island filter (AIF), shown in the setting of Section 5.
Reference $B$ stores FMH anchors $a sub j$ (blue) with positions; the stretches between
consecutive anchors are islands $I sub j$. Tier-2 $k$-mers ($h(x) < s sub 2$, orange) are
inserted into one Bloom filter under keys $(x, j)$, where $j$ is the island index.
The query contig $A$ is raw. Its anchors are matched and chained to $B$'s anchors, which gives
the island that each query position maps to. Its tier-2 $k$-mers are then probed with key
(y, island). Hits (filled) are dense in the aligned span and are false positives only
outside it (hatched), which yields AF; their density inside the span yields containment and ANI.
An anchor destroyed by a mutation (\(mu) merges two islands for the query.
.FE
