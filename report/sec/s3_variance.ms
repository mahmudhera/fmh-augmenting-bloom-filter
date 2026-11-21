.H1 "A unified variance framework"
.H2 "Probing a reference AMQ with the query's sampled k-mers"
.LP
Suppose the reference side keeps an AMQ (for example a BF) of $roman FMH sub s (K sub B )$ with FPR $eps$
and no false negatives. Suppose the query can enumerate $roman FMH sub s (K sub A )$, which is true whenever
the query is a raw sequence or carries its sampled hashes explicitly.
Every query $k$-mer with $h(x) < s$ is probed. Let $N$ be the number of probes and $Y$ the number of positive answers.
A shared $k$-mer always answers positive. By (4), a non-shared one answers positive only through a false positive,
with probability $eps$. The natural estimator is
.DS
	$Chat ~=~ { Y/N - eps } over { 1 - eps } .$	(5)
.DE
With $eps = 0$ this is the FMH containment estimator. With $s = 1$ it is the \(lqcontainment MinHash\(rq
idea of probing a BF of the larger set [8], applied to every $k$-mer.
The design space spanned by $(s, eps )$ therefore contains plain FMH (D0), the dense BF (D1) and the
two-tier filter (D2) as special cases.
.H2 "Main result"
.TH "Theorem 1" "variance decomposition"
Assume the SMM, ideal hashing (independent selection with probability $s$, independent false positives
with probability $eps$), and condition on $N >= 1$. Then
.DS
	$Ex [ Chat | M, N ] = C sub M ,~~~~~~ "so " Ex [ Chat ] = theta " exactly,"$	(6)
.DE
.DS
	$Var ( Chat ) ~=~ { Var ( N sub mut ) } over { n sup 2 } ~+~ { theta (1- theta )(1-s) } over {ns} ~+~ { (1- theta ) eps } over { (1- eps ) ns } ~+~ O left ( 1 over { (ns) sup 2 } right ) ,$	(7)
.DE
where $M$ denotes the mutation process. The three terms are the \f[I]process\f[R], \f[I]sampling\f[R] and
\f[I]false-positive\f[R] variance, and $Var ( N sub mut )/ n sup 2 approx theta B(p,k) / n$ by (2)\(en(3).
.TQ
.LP
The proof (Appendix A) applies the law of total variance over the mutation process, the hash selection
and the false-positive draws. Given $M$ and $N$, the number of shared probes is hypergeometric and the
false positives are binomial, and the FP correction in (5) makes the conditional mean exactly $C sub M$.
Three remarks follow.
(i) The finite-population factor $(1-s)$ makes the sampling term vanish at $s=1$; the false-positive
term does not.
(ii) Spurious (non-homologous) matches enter exactly like false positives, with $eps$ replaced by
$eps prime = eps + sigma - eps sigma$. The difference is that $eps$ is known from the filter's fill
ratio while $sigma$ is not, so $sigma$ must be made negligible (Section 5.5).
(iii) The structure of (7) holds for any tier: if the query probes $k$-mers sampled at rate $s sub 2$
against an AMQ of $roman FMH sub {s sub 2} (K sub B )$, replace $s$ by $s sub 2$.
.TH "Corollary 1" "ANI"
By the delta method, with $ANIh = Chat sup {1/k}$,
.DS
	$Var ( ANIh ) ~ approx ~ { (1-p) sup 2 } over { k sup 2 theta sup 2 } ~ Var ( Chat ) ,~~~~~~ Ex [ ANIh ] - (1-p) ~ approx ~ - ~ { k - 1 } over { 2 k sup 2 } ~ { Var ( Chat ) } over { theta sup 2 } ~ (1-p) .$	(8)
.DE
In the regime $kp << 1$, using $B approx k sup 2 p$ and $q approx kp$,
.DS
	$Var ( ANIh ) ~ approx ~ p over n left [ 1 ~+~ { 1 - s + phi } over { k s } right ] ,~~~~~~ phi = { eps } over { (1 - eps ) theta } .$	(9)
.DE
.TQ
.H2 "Efficiency relative to alignment"
.LP
An alignment of the homologous region observes all $L sub A approx n$ positions. Under the SMM the
maximum-likelihood estimate $M/L$ has variance $p(1-p)/L$, which is the Cram\['e]r\(enRao bound.
We define the efficiency $eta = [p(1-p)/n] / Var ( ANIh )$.
By (9), when $kp << 1$, $eta approx [1 + (1-s+ phi )/(ks)] sup -1$, so the
\f[I]variance-inflation factor\f[R] of exact FMH is
.DS
	$roman VIF (s) ~=~ 1 ~+~ { 1-s } over { ks } ~~~ approx ~~~ 1 over { ks } ~~~ (s << 1/k).$	(10)
.DE
Table 4 and Figure 2 evaluate the exact expressions (7)\(en(8).
.TC 4
Efficiency $eta$ of exact ($eps = 0$) FMH containment-based ANI as a function of the query-side
sampling rate $s$ ($k=21$, large $n$). The VIF column is the small-$p$ factor (10).
The row $s=1$ is the \f[I]containment floor\f[R]: the best that any containment-based estimator can achieve.
.TCE
.so gen/tab_eff.tbl
.LP
Three conclusions follow.
.IP (a)
\f[B]Sampling dominates at practical scales.\f[R]
At $s = 1/1000$ the VIF is 48.6 and $eta$ is between 0.012 and 0.020 across 95\(en99.9% ANI.
Exact FMH at the default scale extracts about 2% of the information an alignment of the same
sequence would give.
.IP (b)
\f[B]The knee is at $s approx 1/k$.\f[R]
Below $1/k$ the selected $k$-mers rarely overlap, each one is a nearly independent test, and variance
scales as $1/(ns)$. Above $1/k$ neighbouring tests overlap and share mutations, and returns diminish.
.IP (c)
\f[B]Even the floor is imperfect at moderate divergence.\f[R]
With every $k$-mer compared, $eta$ is 0.99 at 99.9% ANI, but only 0.69 at 95%, 0.44 at 90% and 0.12 at 80%.
A $k$-mer test reports only \(lqat least one substitution in this window\(rq, and that loses information
once windows typically contain substitutions ($kp$ of order 1). Recovering it needs smaller $k$ or positional
information (Section 6).
.LP
.ne 3.1i
.PSPIC -C figs/fig_eff.eps 5.0i
.FC 2
Efficiency of exact FMH containment-based ANI versus sampling rate $s$ ($k=21$), from (7)\(en(8).
The dotted curve is the small-$p$ approximation $1/(1+(1-s)/(ks))$.
Vertical guides mark the default $s=1/1000$ and the knee $s = 1/k$.
.FE
.H2 "Monte Carlo check of Theorem 1"
.LP
Table 5 compares (7)\(en(8) with a simulation of the idealised model (6,000 replicates per row; script
\f[CR]\%scripts/mc_check_theorem1.py\f[R]). The simulation draws substitutions, marks surviving $k$-mers,
selects each with probability $s$, and turns each non-shared probe positive with probability $eps$.
The variance of $Chat$ agrees to within 2.5%, which is within Monte Carlo error (about 1.8% relative
standard error).
The ANI variance agrees to within 15% whenever the expected number of probes $ns$ is at least about 100.
With about 25 probes the delta method understates it by 1.7\(mu, because $Chat sup {1/k}$ is strongly
nonlinear near 0. In that regime confidence intervals should come from inverting the exact
distribution of $Chat$, as in [6], rather than from (8).
.TC 5
Monte Carlo check of Theorem 1 ($k=21$). \(lqthy\(rq is (7) with the exact (2); ANI variance from (8).
.TCE
.so gen/tab_mc.tbl
.H2 "Why short contigs fail: a diagnosis"
.LP
For a contig of length $L sub A$ the information in FMH comparisons is governed by the expected number of
probes $n s$, not by $L sub A$ alone.
Take a 5 kb contig at 95% ANI with $s = 1/1000$.
It yields $ns approx 5$ probes and an empty sketch with probability $e sup -5 approx 0.7%$.
The ANI standard deviation is 2.84 pp (Table 10), and the downward delta-method bias from (8) is about $-0.9$ pp.
The intrinsic limits at this length are 0.37 pp for the containment floor and 0.31 pp for an alignment oracle.
No sketch can beat these, but D0 is about 8\(mu above them.
Theorem 1 points to three remedies, which organise the rest of this document:
.IP R1. 0.4i
\f[B]More probes per stored bit\f[R]: lower the cost per sampled $k$-mer and accept FP noise (Section 4).
.IP R2. 0.4i
\f[B]Probe only where the contig aligns\f[R]: remove AF confounding and spurious hits using positions (Section 5).
.IP R3. 0.4i
\f[B]More information per probe\f[R]: choose $k$, or test longer units, to match the divergence (Section 6).
