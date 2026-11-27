.H1 "Containment and ANI"
.H2 "One tier"
.LP
Restrict attention for now to a contig that is aligned throughout at identity $a$. In tier $t$ it has $N sub t$
probes, of which $Y sub t$ hit. A shared $k$-mer always hits; a non-shared one hits with probability
$eps prime sub t$. The FP-corrected containment of the tier is
.DS
	$Chat sub t ~=~ { Y sub t / N sub t - eps prime sub t } over { 1 - eps prime sub t } .$	(8)
.DE
Let $C sub M$ be the realised containment, the share of the contig's $k$-mers left unmutated by the mutation
process $M$.
.TH "Theorem 1" "unbiasedness and variance [first theory report]"
Under the simple substitution model and ideal hashing, conditioning on $N sub t >= 1$,
$Ex [ Chat sub t | M, N sub t ] = C sub M$, so $Ex [ Chat sub t ] = theta$ exactly, and
.DS
	$Var ( Chat sub t ) ~=~ { Var ( N sub mut ) } over { n sup 2 } ~+~ { theta (1 - theta )(1 - s sub t ) + (1 - theta ) eps prime sub t / (1 - eps prime sub t ) } over { n s sub t } ~+~ O left ( 1 over { (n s sub t ) sup 2 } right ) .$	(9)
.DE
The first term is the \f[I]process\f[R] variance of the mutation process ($Var ( N sub mut )$ is given exactly
by Blanca et al. [7], about $n theta k sup 2 p$ for $kp << 1$). The second is the \f[I]sampling\f[R] and
\f[I]false-positive\f[R] variance, written $v sub t / (n s sub t )$ below.
.TQ
.LP
The proof is that of the first theory report (Appendix A there), with $eps$ replaced by $eps prime$ and spurious
matches treated like false positives: given $M$ and $N sub t$, the shared probes are hypergeometric and the false
positives binomial, and the affine correction (8) makes the conditional mean $C sub M$.
.H2 "Combining tiers"
.LP
Given $M$, the tiers are independent (disjoint hash ranges, independent false positives), and each $Chat sub t$
has conditional mean $C sub M$. aif combines them with inverse-variance weights:
.DS
	$Chat ~=~ { sum from t w sub t Chat sub t } over { sum from t w sub t } ,~~~~~~ w sub t = { N sub t } over { v sub t ( theta hat ) } ,~~~~~~ v sub t ( theta ) = theta (1 - theta )(1 - s sub t ) + (1 - theta ) { eps prime sub t } over { 1 - eps prime sub t } .$	(10)
.DE
.TH "Proposition 4" "combined estimator"
With weights computed at a fixed $theta$, $Chat$ is conditionally unbiased for $C sub M$, and
$Var ( Chat ) = Var ( C sub M ) + Ex [ ( sum w sub t ) sup -1 ]$, which is the smallest variance of any
fixed-weight average of the tiers. aif evaluates the weights at a plug-in $theta hat$ (two steps, starting from
0.3); this makes $Chat$ biased by $O(1 / ( n s theta ))$.
.TQ
.LP
\f[I]Proof.\f[R] Conditional on $M$ and the $N sub t$, the $Chat sub t$ are independent with means $C sub M$ and
variances $v sub t / N sub t$ to first order; the weighted mean is unbiased, its conditional variance
$( sum w sub t ) sup -1$ is minimal by the Gauss\(enMarkov argument, and the law of total variance adds
$Var ( C sub M )$. Plug-in weights correlate with the estimates through $theta hat$. \(sq
.LP
Table 5 checks Theorem 1 and Proposition 4 by Monte Carlo (20,000 contigs per row, two tiers, aif's
estimator). The standard deviation of $Chat$ agrees with (9)\(en(10) to within 3%, and the ANI standard deviation
to within 4%. The plug-in bias is visible only where probes are few and $theta$ small: +3% relative at 2 kb and
85% ANI, +1.4% at 5 kb, and at most 0.4% from 90% ANI up.
.ne 1.9i
.TC 5
Monte Carlo check of Theorem 1 with two tiers and aif's combination (10), and of the ANI estimator (13).
Simple substitution model, $k = 20$, default tiers, 20,000 contigs per row.
.TCE
.ps 9
.vs 11p
.so gen/tab_mc.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Bits per sample, and what the hedge costs"
.LP
Below the knee $s << 1/k$, accuracy is bought with probes, and a probe costs $c( eps )$ bits of reference. The
first theory report compared a single filter tier with exact 64-bit FracMinHash at the same budget:
.TH "Theorem 2" "bits per sample [first theory report]"
At a budget of $beta$ bits per base, a single tier at FPR $eps$ has non-process variance
$theta (1 - theta ) J( eps ) / ( n beta )$ with
.DS
	$J( eps ) ~=~ c( eps ) left [ 1 - s + { eps } over { (1 - eps ) theta } right ] ,$	(11)
.DE
against $J = 64$ for exact FracMinHash. The minimising FPR grows with $theta$, from about 0.002 at 80% ANI to
0.5 and above from 95% ANI up.
.TQ
.LP
An FPR that is best for containment is not best for everything else. Detection and segmentation (Sections 5 and 6)
depend on how much each probe tells aligned from unaligned, the expected log-likelihood ratio per base,
.DS
	$D( eps ) ~=~ s( eps ) ~ KL left ( pi ( a ) ~||~ eps prime right ) ,~~~~~~ s( eps ) = beta / c( eps ) ,$	(12)
.DE
which is largest at FPRs of 0.0001\(en0.008: a hit must be much more likely in an alignment than outside it.
The default design hedges with two tiers. Table 6 measures what that costs against the best single tier for each
criterion at each identity:
.IP \(bu 0.2i
\f[B]Against exact FracMinHash\f[R], the two tiers have 2.7\(en8.6 times less sampling variance at the same
bits (2.7\(mu at 80% ANI, 7.8\(mu at 95%, 8.6\(mu at 99.9%).
.IP \(bu 0.2i
\f[B]Against the best single tier\f[R], they keep 56\(en81% of the attainable containment efficiency between 85% and
95% ANI but only 26\(en30% at 99% and above, where a single tier at an FPR near 0.6 would be best. For detection they keep 40% at 85%
ANI, 57% at 90% and over 90% at 99%.
.LP
No single FPR serves both criteria across identities: containment wants a high FPR at high identity, detection a
low FPR at low identity. The two-tier default gives up some variance at high identity, where it matters least in
absolute terms (Table 7), to keep detection and AF working down to about 85% ANI.
.ne 2.6i
.TC 6
Efficiency of the default two tiers at $beta = 0.88$ against exact FracMinHash and against the best single tier,
for containment variance ($J$, lower is better) and for detection information $D$ (12). \(lqaif share of best\(rq
is $J sub best / J sub aif$ and $D sub aif / D sub best$.
.TCE
.ps 9
.vs 11p
.so gen/tab_hedge.tbl
.ps \n(PS
.vs \n(VSp
.H2 "ANI and its bias correction"
.LP
ANI is $Chat sup {1/k}$, corrected for the curvature of the $k$-th root by the second-order delta method:
.DS
	$ANIh ~=~ Chat sup {1/k} left [ 1 + { k - 1 } over { 2 k sup 2 } ~ min left ( V over { Chat sup 2 } , 1 right ) right ] ,~~~~~~ Var ( ANIh ) ~ approx ~ { a sup 2 } over { k sup 2 theta sup 2 } ~ V ,$	(13)
.DE
where $V$ is (9)\(en(10) evaluated at the estimate. The relative variance $V / theta sup 2$ grows as $theta$
falls, which is why ANI accuracy degrades below 90%. With very few probes the delta method under-states the variance,
and intervals from the exact distribution of the containment are preferable [5]. In Table 5 the corrected estimator's bias is at most 0.04 pp
from 90% ANI up and +0.28 pp at 2 kb and 85%, where the plug-in bias of Proposition 4 dominates. Table 7 puts the
standard deviation in context: at the same bits, aif's two tiers roughly halve the error of exact FracMinHash, and
they sit 1.1\(en1.6 times above the floor set by the mutation process itself (every $k$-mer compared).
.ne 2.9i
.TC 7
Standard deviation of ANI (pp) for a fully aligned contig, from (9)\(en(13): aif's two tiers, exact 64-bit
FracMinHash at the same 0.88 bits per base ($s = 0.0138$), every $k$-mer compared (the containment floor), and
an alignment of the same contig ($sqrt { p(1-p)/L }$).
.TCE
.ps 9
.vs 11p
.so gen/tab_sd.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Variance given the realised identity"
.LP
Alignment-based truth, such as MUMmer's, is the contig's \f[I]realised\f[R] identity $1 - M/L$, where $M$ is its
number of substitutions, not the rate $p$ that generated them. Intervals for that estimand should not include
the variance of $M$ itself.
.TH "Proposition 5" "conditional variance"
To first order in $1/L$,
.DS
	$Var ( Chat | M ) ~ approx ~ Var ( Chat ) ~-~ { k sup 2 theta sup 2 p } over { (1 - p) n } .$	(14)
.DE
.TQ
.LP
\f[I]Proof.\f[R] $Ex [ C sub M | M ] approx (1 - M/L) sup k$, and $M/L$ has variance $p(1-p)/L$; the delta method
gives $Var ( Ex [ C sub M | M ]) approx k sup 2 theta sup 2 p / ((1-p) L)$, and the law of total variance removes it
from the process term. \(sq
.LP
What remains of the process variance is the arrangement of the substitutions (how many $k$-mers each one kills,
which depends on their spacing). The realised count explains most of the process variance at high identity and
little at low: 94% of it at 99% ANI, 71% at 95%, 46% at 90% and 27% at 85%. aif's statistical ANI intervals use
(14) with $n$ replaced by the aligned $k$-mer count.
