.H1 "Departures from the model, and cost"
.H2 "Indels"
.LP
The substitution model has no indels. An insertion of $l$ bases into the contig removes $k + l - 1$ of its
$k$-mers from the reference's set, and a deletion removes the $k - 1$ that span the junction.
.TH "Proposition 15" "indels act as mismatches"
With substitutions at rate $p$ and single-base indels at rate $r$ per base, $theta approx (1-p) sup k (1-r) sup k$,
so aif's ANI estimates $(1 - p)(1 - r) approx 1 - p - r$: the identity that counts each indel column as a
mismatch (as MUMmer's percent identity does [16]), not the one that excludes indel columns. AF is unaffected by
short indels. A long insertion is an unaligned stretch (Section 5.4); a long deletion is invisible to AF.
.TQ
.LP
Table 15 confirms this on 10 kb contigs at 97% substitution identity: with 1% indels aif reads 96.12%, against
96.04% for the mismatch convention and 96.99% for the indel-free one, and AF stays at 99.7\(en99.8%. The
convention matters when comparing tools: against a truth that excludes indel columns, aif's ANI is low by about the
indel rate.
.ne 1.6i
.TC 15
Contigs with substitutions ($p$ = 3%) and single-base indels at rate $r$ (half insertions), 200 per row: the two
alignment conventions, Proposition 15's prediction, and aif.
.TCE
.ps 9
.vs 11p
.so gen/tab_indel.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Other departures"
.IP \(bu 0.2i
\f[B]Repeats and low-complexity sequence\f[R] raise the spurious-match rate (Table 4), which aif calibrates and
counters with a larger $k$. A reference whose null hit rate exceeds the calibrated one breaks Theorem 4's
condition.
.IP \(bu 0.2i
\f[B]Clustered substitutions\f[R] (mutation hot spots, rate variation along the genome) kill fewer $k$-mers per
substitution than independent ones, so $theta > (1-p) sup k$ and ANI reads high; at the extreme of substitutions
in pairs closer than $k$, the deficit $1 - theta$ is nearly halved. They also inflate the variance beyond (9). Both
are absorbed, on average, by the calibrated intervals, but not removed from the point estimate.
.IP \(bu 0.2i
\f[B]Rearrangements and chimeras\f[R] are invisible: aif uses no reference positions, so a contig assembled from two
distant loci of the same genome is fully aligned. This was the price of dropping the anchors of the first design.
.IP \(bu 0.2i
\f[B]Composition\f[R] does not bias selection (hashing is uniform), but GC-skewed genomes have more spurious
matches at a given $k$ (\f[I]P. aeruginosa\f[R] in Table 4).
.H2 "Cost"
.LP
All operations are linear and streaming.
.IP \(bu 0.2i
\f[B]Sketch:\f[R] one pass over the reference, $s L sub B$ insertions of one cache line each, $beta L sub B$ bits.
.IP \(bu 0.2i
\f[B]Query:\f[R] one pass over the contig to hash and select its $k$-mers, independent of the number of
references. Per reference: $s L$ filter tests of one cache line each, and the HMM in $O(( s L / 32 ) m)$ for the
binned pass, plus $O(m)$ per probe in the refined bins around boundaries and, for heterogeneous contigs, a pass with
bins of 8 probes.
.IP \(bu 0.2i
\f[B]Database:\f[R] the screen reads $3 n sub s$ rows of 64 bytes per group of 512 references, about 17,700 rows for
a 1 kb contig against GTDB's 394 groups, then evaluates at most 16 references. The screening index adds 0.13 bits
per base (8 bits per key at 1/64 of $k$-mers), so a database costs 1.0 bits per reference base in all.
