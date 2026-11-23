.nf
.ta 0.45i 6.5iR
.lc .
.ps 10
.vs 11.9p
.sp 0.25v
\f[B]1.	Summary of results\f[]1
.sp 0.25v
\f[B]2.	Data and ground truth\f[]2
\h'0.45i'2.1\h'|1.05i'Genome sets2
\h'0.45i'2.2\h'|1.05i'Alignment-based truth3
\h'0.45i'2.3\h'|1.05i'Tools, budgets and metrics3
.sp 0.25v
\f[B]3.	Contig-level accuracy\f[]3
\h'0.45i'3.1\h'|1.05i'Per set3
\h'0.45i'3.2\h'|1.05i'By contig length4
\h'0.45i'3.3\h'|1.05i'By divergence and by aligned fraction5
\h'0.45i'3.4\h'|1.05i'Non-homologous contigs and negative controls6
.sp 0.25v
\f[B]4.	Genome-level ANI and AF\f[]7
.sp 0.25v
\f[B]5.	How reliable is the truth?\f[]7
\h'0.45i'5.1\h'|1.05i'Repeats and seeding7
\h'0.45i'5.2\h'|1.05i'A blastn audit7
\h'0.45i'5.3\h'|1.05i'Other caveats8
.sp 0.25v
\f[B]6.	Real-data tuning (aif-RT)\f[]8
.sp 0.25v
\f[B]7.	Speed and memory on real data (v1.2\(en1.4)\f[]8
\h'0.45i'7.1\h'|1.05i'What real assemblies exposed9
\h'0.45i'7.2\h'|1.05i'Changes in v1.29
\h'0.45i'7.3\h'|1.05i'Matching skani single-threaded (v1.3)9
\h'0.45i'7.4\h'|1.05i'Timing10
.sp 0.25v
\f[B]8.	Calibrated confidence intervals (v1.4)\f[]11
\h'0.45i'8.1\h'|1.05i'Statistical intervals11
\h'0.45i'8.2\h'|1.05i'Calibration against alignment truth12
.sp 0.25v
\f[B]9.	Toward a principled AF model\f[]13
.sp 0.25v
\f[B]10.	Database mode (v1.5\(en1.6)\f[]14
\h'0.45i'10.1\h'|1.05i'Design15
\h'0.45i'10.2\h'|1.05i'Cost of one shared k15
\h'0.45i'10.3\h'|1.05i'Benchmark database16
\h'0.45i'10.4\h'|1.05i'Accuracy of the top hit16
\h'0.45i'10.5\h'|1.05i'Screening17
\h'0.45i'10.6\h'|1.05i'Speed and scaling17
.sp 0.25v
\f[B]11.	GTDB-scale database: contig classification\f[]18
.sp 0.25v
\f[B]12.	Contamination screening of bins\f[]20
.sp 0.25v
\f[B]13.	Size, accuracy and time: a Pareto comparison\f[]21
.sp 0.25v
\f[B]14.	Identity maps along contigs\f[]22
.sp 0.25v
\f[B]15.	Discussion\f[]24
\h'0.45i'15.1\h'|1.05i'Why the simulation results carry over24
\h'0.45i'15.2\h'|1.05i'Limitations24
\h'0.45i'15.3\h'|1.05i'Next steps24
.sp 0.25v
\f[B]16.	Reproducibility\f[]25
.lc
.ps \n(PS
.vs \n(VSp
.fi
