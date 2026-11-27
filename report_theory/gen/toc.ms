.nf
.ta 0.45i 6.5iR
.lc .
.ps 10
.vs 11.9p
.sp 0.25v
\f[B]1.	Introduction\f[]1
.sp 0.25v
\f[B]2.	Model and notation\f[]2
\h'0.45i'2.1\h'|1.05i'Contigs, segments and the estimands2
\h'0.45i'2.2\h'|1.05i'Hashing and consistency2
.sp 0.25v
\f[B]3.	The sketch\f[]3
\h'0.45i'3.1\h'|1.05i'Tier layout3
\h'0.45i'3.2\h'|1.05i'The false-positive rate of a blocked filter is known exactly3
\h'0.45i'3.3\h'|1.05i'Spurious matches4
.sp 0.25v
\f[B]4.	Containment and ANI\f[]5
\h'0.45i'4.1\h'|1.05i'One tier5
\h'0.45i'4.2\h'|1.05i'Combining tiers5
\h'0.45i'4.3\h'|1.05i'Bits per sample, and what the hedge costs6
\h'0.45i'4.4\h'|1.05i'ANI and its bias correction7
\h'0.45i'4.5\h'|1.05i'Variance given the realised identity8
.sp 0.25v
\f[B]5.	Alignment fraction\f[]8
\h'0.45i'5.1\h'|1.05i'Counts cannot separate AF from ANI8
\h'0.45i'5.2\h'|1.05i'The segment HMM8
\h'0.45i'5.3\h'|1.05i'How precisely a boundary is placed9
\h'0.45i'5.4\h'|1.05i'Low levels absorb unaligned stretches10
\h'0.45i'5.5\h'|1.05i'What the local-identity AF of v1.7 can remove11
\h'0.45i'5.6\h'|1.05i'Identity-thresholded AF11
.sp 0.25v
\f[B]6.	Detection and false reports\f[]12
\h'0.45i'6.1\h'|1.05i'A bound on false reports that needs no tuning12
\h'0.45i'6.2\h'|1.05i'The shortest contig that can be detected12
.sp 0.25v
\f[B]7.	Contigs with heterogeneous identity\f[]14
\h'0.45i'7.1\h'|1.05i'Pooled containment over-states mosaic identity14
\h'0.45i'7.2\h'|1.05i'The heterogeneity test and the posterior-mean ANI15
.sp 0.25v
\f[B]8.	Uncertainty\f[]15
\h'0.45i'8.1\h'|1.05i'Statistical intervals15
\h'0.45i'8.2\h'|1.05i'Calibrated intervals16
.sp 0.25v
\f[B]9.	Database search\f[]16
\h'0.45i'9.1\h'|1.05i'Exactness16
\h'0.45i'9.2\h'|1.05i'The screening index17
\h'0.45i'9.3\h'|1.05i'Recall at database scale17
\h'0.45i'9.4\h'|1.05i'The pre-model filter18
.sp 0.25v
\f[B]10.	Departures from the model, and cost\f[]18
\h'0.45i'10.1\h'|1.05i'Indels18
\h'0.45i'10.2\h'|1.05i'Other departures18
\h'0.45i'10.3\h'|1.05i'Cost19
.sp 0.25v
\f[B]11.	What the theory says, and what it suggests\f[]19
\h'0.45i'11.1\h'|1.05i'Summary19
\h'0.45i'11.2\h'|1.05i'Improvements the theory points to19
.sp 0.25v
\f[B]	References\f[]20
.sp 0.25v
\f[B]A.	Simulation details\f[]20
.sp 0.25v
\f[B]B.	Reproduction\f[]21
.lc
.ps \n(PS
.vs \n(VSp
.fi
