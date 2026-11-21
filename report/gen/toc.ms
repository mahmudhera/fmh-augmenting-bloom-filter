.nf
.ta 0.45i 6.5iR
.lc .
.ps 10
.vs 11.9p
.sp 0.25v
\f[B]1.	Introduction and problem statement\f[]1
.sp 0.25v
\f[B]2.	Model, notation and building blocks\f[]2
\h'0.45i'2.1\h'|1.05i'Sequences, k-mers and the simple mutation model2
\h'0.45i'2.2\h'|1.05i'FracMinHash and the consistency property3
\h'0.45i'2.3\h'|1.05i'Approximate membership query (AMQ) structures3
\h'0.45i'2.4\h'|1.05i'Designs studied4
.sp 0.25v
\f[B]3.	A unified variance framework\f[]4
\h'0.45i'3.1\h'|1.05i'Probing a reference AMQ with the query's sampled k-mers4
\h'0.45i'3.2\h'|1.05i'Main result4
\h'0.45i'3.3\h'|1.05i'Efficiency relative to alignment5
\h'0.45i'3.4\h'|1.05i'Monte Carlo check of Theorem 16
\h'0.45i'3.5\h'|1.05i'Why short contigs fail: a diagnosis6
.sp 0.25v
\f[B]4.	Space, linearity and the bits-per-sample principle\f[]7
\h'0.45i'4.1\h'|1.05i'Cost accounting7
\h'0.45i'4.2\h'|1.05i'Dense filters are never space-efficient7
\h'0.45i'4.3\h'|1.05i'Linearity7
\h'0.45i'4.4\h'|1.05i'The bits-per-sample principle8
\h'0.45i'4.5\h'|1.05i'Unknown divergence: hedged tiers9
\h'0.45i'4.6\h'|1.05i'Short contigs at an equal budget9
\h'0.45i'4.7\h'|1.05i'The asymmetry principle10
.sp 0.25v
\f[B]5.	Positions, islands and alignment fraction\f[]11
\h'0.45i'5.1\h'|1.05i'Anchors and islands11
\h'0.45i'5.2\h'|1.05i'Why containment alone cannot give alignment fraction12
\h'0.45i'5.3\h'|1.05i'Alignment fraction from positions12
\h'0.45i'5.4\h'|1.05i'ANI over the aligned span13
\h'0.45i'5.5\h'|1.05i'Island-keyed filters and localisation13
\h'0.45i'5.6\h'|1.05i'Indels14
\h'0.45i'5.7\h'|1.05i'Choice of anchors14
.sp 0.25v
\f[B]6.	Information per stored bit: what to test\f[]14
\h'0.45i'6.1\h'|1.05i'Choosing k14
\h'0.45i'6.2\h'|1.05i'Island fingerprints: testing the island itself15
\h'0.45i'6.3\h'|1.05i'Dense positional match vectors16
\h'0.45i'6.4\h'|1.05i'Rate heterogeneity16
.sp 0.25v
\f[B]7.	A concrete proposal: the AIF-ANI estimator\f[]16
\h'0.45i'7.1\h'|1.05i'Reference sketch16
\h'0.45i'7.2\h'|1.05i'Query algorithm (raw contig)17
\h'0.45i'7.3\h'|1.05i'Starting parameters17
\h'0.45i'7.4\h'|1.05i'Predictions to test17
.sp 0.25v
\f[B]8.	Scenario analysis\f[]18
.sp 0.25v
\f[B]9.	Related work and close ideas\f[]20
.sp 0.25v
\f[B]10.	Limitations, open questions and experimental plan\f[]21
\h'0.45i'10.1\h'|1.05i'Limitations of this analysis21
\h'0.45i'10.2\h'|1.05i'Open questions21
\h'0.45i'10.3\h'|1.05i'Experimental plan for the implementation phase21
.sp 0.25v
\f[B]11.	Conclusions\f[]22
.sp 0.25v
\f[B]	References\f[]22
.sp 0.25v
\f[B]A.	Proof of Theorem 1\f[]24
.sp 0.25v
\f[B]B.	Bloom-filter cardinality variance (Proposition 2c)\f[]24
.sp 0.25v
\f[B]C.	Additional tables\f[]25
.lc
.ps \n(PS
.vs \n(VSp
.fi
