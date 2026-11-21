.\" Build: see ../Makefile  (groff -s -t -e -k -Tps -M tmac -m fmhdoc main.ms)
.EQ
delim $$
gsize 11
define eps 'epsilon'
define Chat '{ C hat }'
define ANIh '{ roman "ANI" } hat'
define Var 'roman "Var"'
define Ex 'roman "E"'
define Pr 'roman "Pr"'
.EN
.\" ---------------------------------------------------------------- title page
.nr BODY 0
.sp 1.1i
.ce 3
.ps 18
.vs 23p
.ft B
Augmenting FracMinHash with Bloom Filters
and Positional Anchors for ANI and
Alignment-Fraction Estimation on Short Contigs
.ft R
.ps 13
.vs 17p
.sp 0.6v
.ce 1
A theoretical framework (pre-implementation analysis)
.sp 1.0v
.ps 10.5
.vs 13p
.ce 3
Project: \f[CR]fmh-augmenting-bloom-filter\f[R]
Working draft v0.1 \(em November 21, 2026
Analysis and typesetting drafted with Claude (Anthropic) for the project owner
.ps \n(PS
.vs \n(VSp
.sp 1.4v
.in 0.5i
.ll 6.0i
.ce 1
\f[B]Abstract\f[R]
.sp 0.3v
.fi
.ad b
.ps 10
.vs 12.4p
Containment-based ANI estimation with FracMinHash (FMH),
$ANIh = Chat sup {1/k}$, degrades sharply on short contigs, and it cannot
separate ANI from the alignment fraction (AF).
We build a single variance framework that splits the error of any
FMH-style estimator into three parts: mutation-process variance,
query-side sampling variance, and false-positive noise from approximate
membership query (AMQ) structures such as Bloom filters (BFs).
The resulting closed form (Theorem 1, checked by Monte Carlo) shows that at the usual
scale $s = 1/1000$ with $k=21$ the sampling term is about $1/(ks) approx 48$ times the
irreducible term. That is the root cause of the short-contig problem.
Inside this framework we analyse the proposal to store FMH anchor
positions plus a Bloom filter over the \(lqislands\(rq between anchors,
together with 13 related designs, across 16 scenarios.
(1) A filter that holds \f[I]all\f[R] island $k$-mers is never space-efficient:
any AMQ with false-positive rate (FPR) $eps < 1/4$ costs more than the 2-bit sequence itself (Theorem 2).
(2) A subsampled island tier (a second FMH threshold $s sub 2 > s$) is linear in
time and space. At an equal bit budget it reduces sampling variance by a factor
$R sup * = 64 / J sup *$, which ranges from 4.1 (ANI 80%) through 11.3 (95%) to 23 (99.9%).
The optimal FPR is surprisingly high ($eps sup * approx 0.5$ for ANI $>=$ 95%; Theorem 3).
For a 5 kb contig at 95% ANI the standard deviation falls from 2.84 to 0.91 percentage points.
(3) BF augmentation is intrinsically asymmetric: the query must enumerate its $k$-mers.
BF-versus-BF comparison is unusable for short contigs.
(4) Positions are what make AF identifiable. Containment alone gives
$ANIh -> alpha sup {1/k} (1-p)$, and tier-2 hits sharpen AF boundaries from about $1/(s theta )$
to about $1/(s sub 2 theta )$ bases.
(5) Keying the filter by island identity localises probes. Spurious matches fall by a factor
equal to the number of islands, which admits small $k$ at low ANI (up to a 2.9\(mu information gain
at 80%) and rescues the reversed asymmetry (contigs in the database).
(6) For ANI $>=$ 99.9%, a related idea, storing a short fingerprint per island, beats every
$k$-mer design per stored bit.
We give a concrete estimator (AIF-ANI), a regime map and an experimental plan for the implementation phase.
.ps \n(PS
.vs \n(VSp
.in 0
.ll 6.5i
.\" ---------------------------------------------------------------- contents
.bp
.ce 1
\f[B]\s(13Contents\s0\f[R]
.sp 0.8v
.so gen/toc.ms
.\" ---------------------------------------------------------------- body
.nr BODY 1
.pn 1
.bp
.so sec/s1_intro.ms
.so sec/s2_model.ms
.so sec/s3_variance.ms
.so sec/s4_space.ms
.so sec/s5_positions.ms
.so sec/s6_information.ms
.so sec/s7_estimator.ms
.so sec/s8_scenarios.ms
.so sec/s9_related.ms
.so sec/s10_conclusions.ms
.so sec/refs.ms
.so sec/appendix.ms
