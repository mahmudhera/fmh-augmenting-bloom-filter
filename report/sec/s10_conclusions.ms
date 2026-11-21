.H1 "Limitations, open questions and experimental plan"
.H2 "Limitations of this analysis"
.IP \(bu 0.2i
Idealised hashing and the SMM. Real genomes have repeated and low-complexity $k$-mers, indels and
non-uniform mutation. Only indels (Section 5.6) and heterogeneity (Section 6.4) are treated, and only to first order.
.IP \(bu 0.2i
The ANI variances use the delta method, which understates spread below about 30\(en100 probes (Table 5).
.IP \(bu 0.2i
Chaining is modelled only through matched-anchor counts (Proposition 4). Chaining errors, rearrangements and
orientation are not modelled.
.IP \(bu 0.2i
Spurious matches use a uniform random-sequence model (Table 15), which is optimistic for real genomes.
.IP \(bu 0.2i
Components are optimised separately. The joint allocation of a budget between anchors (localisation), the filter
tier (estimation) and island fingerprints is not solved.
.H2 "Open questions"
.IP 1. 0.3i
Joint budget allocation $( s sub a , s sub 2 , eps , f)$ for a given contig-length and ANI distribution.
.IP 2. 0.3i
Exact, non-asymptotic confidence intervals for the FP-corrected estimator (5), extending [6].
.IP 3. 0.3i
An error theory for HMM segmentation of AF, and the best tier-2 rate when AF rather than ANI is the target.
.IP 4. 0.3i
How much of the repeat-induced bias on real assemblies does island keying remove?
.IP 5. 0.3i
Multi-$k$ tiers versus a single $k$ with localisation, at low ANI.
.IP 6. 0.3i
Island fingerprints under sequencing errors, and the multi-resolution (Merkle) variant.
.H2 "Experimental plan for the implementation phase"
.IP E1. 0.4i
\f[B]Controlled simulation on real genomes.\f[R] Mutate real bacterial genomes (substitutions, plus indels in a second
arm). Cut contigs of 1\(en50 kb, including partial overlaps and chimeras. Compare D0, D2, D5, D6 and D8 at equal
budgets on bias, variance and CI coverage, against Theorem 1 and Tables 10 and 21\(en22.
.IP E2. 0.4i
\f[B]Real genome pairs.\f[R] Species- and genus-level pairs with alignment-based ANI and AF as ground truth;
baselines sourmash (FMH), skani and FastANI.
.IP E3. 0.4i
\f[B]MAG contigs.\f[R] Contig-level assignment and ANI on real metagenome assemblies.
.IP E4. 0.4i
\f[B]Ablations.\f[R] Island versus global keys; FPR choices and hedging; $k sub e$ (21 versus 10\(en12); FMH versus syncmer anchors.
.IP E5. 0.4i
\f[B]Resources.\f[R] Sketch size, build and query time, and peak memory; check the linearity claims of Table 7.
.H1 "Conclusions"
.LP
\f[B]Q1 (space).\f[R] A Bloom filter of all island $k$-mers is never space-efficient: it costs more than the 2-bit
genome for any useful FPR (Theorem 2). With a second FMH threshold the island tier is sketch-sized and can match
the plain FMH budget exactly.
For short contigs the anchors, which must be dense enough to land in the contig, are the larger cost.
.LP
\f[B]Q2 (linearity).\f[R] Yes: single-pass, linear-time construction and querying, and linear space with a tunable constant.
.LP
\f[B]Q3 (ANI).\f[R] Augmenting FMH with a filter does improve ANI estimates for short contigs, through one mechanism:
many cheap, noisy probes instead of a few expensive, exact ones. The exact variance is given by Theorem 1.
At equal space the gain is $R sup * = 64 / J sup *$, from 4\(mu to 23\(mu in variance, with an unexpectedly high
optimal FPR at species and strain level, about 3\(mu in standard deviation for a 5 kb contig at 95% ANI.
Reaching the containment floor is not sketch-sized (0.3\(en1.2 bits per base).
The gain requires the query to enumerate its $k$-mers; BF-versus-BF comparison fails for short contigs.
.LP
\f[B]Q4 (AF).\f[R] Containment cannot separate AF from ANI (Proposition 3), and ignoring AF costs several pp.
Positions make AF identifiable, query-side positions suffice for AF of the query (Proposition 5), and tier-2 hits
locate boundaries to about $1/(s sub 2 theta )$ bases.
Reference-side anchors add colinearity (chimeras), localisation and indel detection.
.LP
\f[B]Q5 (close ideas).\f[R] The positional part of the proposal pays off most through localisation: island keys
suppress spurious and paralogous matches, allow small $k$ at low ANI and fix the reversed asymmetry.
Its most powerful relative is the island fingerprint, which tests each island as a single long unit and is the best
design per bit at ANI $>=$ 99.9%.
The recommended system combines these as AIF-ANI (Section 7) and is ready for the experiments above.
