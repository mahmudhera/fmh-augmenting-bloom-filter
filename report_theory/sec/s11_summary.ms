.H1 "What the theory says, and what it suggests"
.H2 "Summary"
.IP \(bu 0.2i
\f[B]Containment and ANI\f[R] are unbiased up to $O(1/( n s theta ))$, with a closed-form variance that simulation
confirms to within 3%. At the same bits, aif's standard deviation is about half that of exact FracMinHash and
1.1\(en1.6 times the floor of comparing every $k$-mer. The two-tier hedge keeps detection working down to about
85% ANI at the price of 26\(en30% of the attainable containment efficiency at 99% ANI and above.
.IP \(bu 0.2i
\f[B]AF\f[R] is identifiable from query positions alone, but its resolution is limited by the model's lowest
identity levels, which absorb unaligned stretches of up to about 1.5 kb inside alignments and 0.5 kb at contig
ends. v1.7's local-identity AF removes long absorbed stretches and leaves short ones.
.IP \(bu 0.2i
\f[B]Detection\f[R] has a parameter-free false-report bound per reference, $e sup { - lambda }$, which the data
beat by two orders of magnitude, and a predictable sensitivity: half of contigs are detected from 172 bp at 95%
ANI, 854 bp at 90%, 3.5 kb at 86%, and none below about 82% at $k = 20$.
.IP \(bu 0.2i
\f[B]Mosaic contigs\f[R] would read high by a power-mean bias; the heterogeneity test and the posterior-mean ANI
remove it unless part of the contig lies below the level grid.
.IP \(bu 0.2i
\f[B]Intervals\f[R] are calibrated per stratum with a finite-sample guarantee that holds only for contigs resembling
the calibration set.
.IP \(bu 0.2i
\f[B]Database search\f[R] is exactly pairwise search on the references it evaluates; at GTDB scale, the screen
decides which those are, and it misses short divergent contigs.
.H2 "Improvements the theory points to"
.IP (1) 0.3i
\f[B]A prior on identity for level switches.\f[R] Absorption is driven by the uniform redraw over 28 levels, which
gives the near-unaligned levels the same prior as the level of the alignment. By (17), dividing the redraw prior of
the lowest eight levels by 10 would shorten the break-even length from 1.5 kb to 0.8 kb inside alignments and
from 0.5 kb to 0.2 kb at contig ends, in the model itself rather than through window heuristics. It should also sharpen boundaries, since the same low
levels blur them (Table 8, column \(lqall levels\(rq).
.IP (2) 0.3i
\f[B]Detection floor versus AF bias.\f[R] The grid floor $theta sub m = 8 eps prime sub min$ sets both the detection
floor (82% at $k = 20$) and the number of near-unaligned levels. Moving it is a trade, which (17) and (19) quantify.
.IP (3) 0.3i
\f[B]Screening for short contigs.\f[R] Re-ranking the top 1,024 screened references by all-probe containment raises
recall at 1 kb and 95% ANI from 72% to 96% (Table 14). Denser screening keys would do the same at a cost in bits.
.IP (4) 0.3i
\f[B]A database-scale threshold.\f[R] A family-wise guarantee over $R$ references needs $lambda >= ln R + ln (1/ alpha )$;
the empirical null tail suggests a smaller increase suffices, which a database null experiment could settle.
.IP (5) 0.3i
\f[B]Pre-model filter.\f[R] Replacing the heuristic with an upper bound on $ln LR$ would make it provably
output-preserving.
.IP (6) 0.3i
\f[B]Reference positions.\f[R] Rearrangements and chimeras need positions on the reference, at the cost analysed
in the first theory report.
