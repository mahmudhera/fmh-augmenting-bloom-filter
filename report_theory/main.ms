.\" Build: report_theory/build.sh   (numbers: scripts/theory_numbers.py, theory_sims.py, theory_hmm.py, theory_report.py)
.EQ
delim $$
gsize 11
define eps 'epsilon'
define Chat '{ C hat }'
define ANIh '{ roman "ANI" } hat'
define Var 'roman "Var"'
define Ex 'roman "E"'
define Pr 'roman "Pr"'
define KL 'roman "KL"'
define LR 'roman "LR"'
.EN
.nr BODY 0
.sp 1.0i
.ce 3
.ps 18
.vs 23p
.ft B
The Theory of aif
.ft R
.ps 14
.vs 18p
Tiered Bloom-Filter FracMinHash Sketches for
Contig ANI and Alignment Fraction
.ps 13
.vs 17p
.sp 0.8v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Theory report, aif v1.7 \(em November 27, 2026
Analysis, simulations and typesetting by Claude (Anthropic) for the project owner
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
aif estimates the identity (ANI) and the aligned fraction (AF) of a contig against a reference from a sketch
of the reference alone: FracMinHash-selected $k$-mers stored in two Bloom-filter tiers with false-positive
rates (FPRs) of 0.6% and 10%, read along the contig by a hidden Markov model (HMM). This report gives the theory
of the method as implemented in v1.7 and checks every result twice, against Monte Carlo simulation of the
idealised model and against the aif binary on contigs simulated from a random 5 Mbp reference.
(1) The tier layout does not depend on the genome, which makes database search exactly equal to pairwise
queries. The FPR of a blocked filter is known exactly from its fill (0.0064 predicted, 0.0063 measured).
(2) The FP-corrected containment is unbiased with a closed-form variance (Theorem 1); the variance agrees with
simulation to within 3%. At 0.88 bits per base the two tiers give 2.7\(en8.6\(mu less sampling variance than exact
FracMinHash, but only 26\(en81% of what the best single tier would give, because the FPR that is best for
containment is the worst for detection (Theorem 2).
(3) Containment counts alone cannot separate AF from ANI. In the HMM, AF is limited less by the change-point
information of the two states than by the lowest identity levels, which look almost unaligned: an unaligned
stretch inside an alignment is counted as aligned until it is about 1.5 kb long, and one at a contig end until
about 0.5 kb (Theorem 3; aif: 86% of a 1 kb gap absorbed, against 85% predicted). This is the mechanism behind
the +5 pp AF bias that v1.7 corrected, and it bounds what that correction can do.
(4) Whatever the HMM's parameters, an unrelated contig reaches a likelihood ratio (LR) of $e sup lambda$ with
probability at most $e sup { - lambda }$ (Theorem 4); none of 65,000 unrelated contigs reached the reporting
threshold. The expected LR predicts the shortest detectable contig at each identity: 854 bp at 90% ANI,
3.5 kb at 86%, none below about 82% at $k = 20$, matching aif.
(5) Pooled containment over-states the ANI of mosaic contigs by a power-mean bias, up to 1.8 pp for halves at
99% and 90%, which aif's heterogeneity test removes.
(6) At GTDB scale the database screen, not the model, limits short contigs: with 15 screening keys a 1 kb contig
at 95% ANI must outrank about 13,000 unrelated references that pass by chance, and does so 72% of the time.
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
.so sec/s1_intro.ms
.so sec/s2_model.ms
.so sec/s3_sketch.ms
.so sec/s4_containment.ms
.so sec/s5_af.ms
.so sec/s6_detect.ms
.so sec/s7_het.ms
.so sec/s8_ci.ms
.so sec/s9_db.ms
.so sec/s10_departures.ms
.so sec/s11_summary.ms
.so sec/refs.ms
.so sec/appendix.ms
