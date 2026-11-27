.H1 "Contigs with heterogeneous identity"
.H2 "Pooled containment over-states mosaic identity"
.LP
Recombinant and mosaic contigs contain aligned segments at different identities. ANI (1) is their
length-weighted mean, but containment pools $k$-mer survival, which is a convex function of identity.
.TH "Proposition 10" "power-mean bias"
Let the aligned part consist of segments with identities $a sub i$ and length shares $w sub i$. The pooled
containment converges to $sum w sub i a sub i sup k$, so the pooled ANI converges to the power mean
$M sub k = ( sum w sub i a sub i sup k ) sup {1/k}$, which exceeds the true ANI $M sub 1 = sum w sub i a sub i$
unless all $a sub i$ are equal. For a small spread,
.DS
	$M sub k - M sub 1 ~ approx ~ { (k - 1) ~ Var sub w (a) } over { 2 M sub 1 } .$	(20)
.DE
.TQ
.LP
\f[I]Proof.\f[R] The power-mean inequality gives $M sub k >= M sub 1$ for $k >= 1$; (20) is the second-order
expansion of $x sup k$ around $M sub 1$. \(sq
.LP
The bias is large for strongly mixed contigs (Table 12): halves at 99% and 90% read 96.3% instead of 94.5%, +1.8
pp, and halves at 99.5% and 85% read 96.3% instead of 92.3%. Expansion (20) is accurate to within 0.05 pp up to
a 4 pp spread and over-states the bias beyond. aif run without its heterogeneity handling (\f[CR]--het 0\f[R])
returns the power mean to within 0.03 pp on simulated mosaics (Table 13), except where part of the contig lies
near the bottom of the grid.
.ne 2.3i
.TC 12
Power-mean bias (Proposition 10) for two-part contigs at $k = 20$.
.TCE
.ps 9
.vs 11p
.so gen/tab_powermean.tbl
.ps \n(PS
.vs \n(VSp
.H2 "The heterogeneity test and the posterior-mean ANI"
.LP
aif tests for heterogeneity over windows of about 1 kb made of confidently aligned bins (posterior above 0.9).
With $W$ windows, window containments $theta hat sub w$ and the pooled $theta hat$, the statistic is
.DS
	$Q ~=~ sum from w { ( theta hat sub w - theta hat ) sup 2 } over { V sub w + theta hat B(p,k) / l sub w } ,$	(21)
.DE
where $V sub w$ is the window's sampling and false-positive variance (10), $theta hat B(p,k) / l sub w$ its process
variance (Theorem 1, $B(p,k) approx k sup 2 p$ [7]) and $l sub w$ its aligned length. Under homogeneity $Q$ is
approximately $chi sup 2$ with $W - 1$ degrees of freedom, and the contig is called heterogeneous when
$Q > W - 1 + 3 sqrt { 2(W-1) }$, the normal approximation at $z = 3$. Because $chi sup 2$ is right-skewed, the
actual size is larger than the nominal 0.13%: 1.0% with 10 windows, 0.69% with 20, 0.34% with 100. On
homogeneous 10 kb contigs at 95% aif called 1.0%, exactly the $chi sup 2$ value, so the variance in (21) is right.
For mosaics with halves 4 pp or more apart the test fired on 94\(en100% of contigs (Table 13).
.LP
When the test fires, aif reports the posterior mean of the level identity over aligned probes, from a finer pass
(bins of 8 probes). This estimates $M sub 1$ directly, since the level posterior tracks each segment's identity.
On the simulated mosaics its error is \(mi0.14 to \(mi0.21 pp, against +0.37 to +1.80 pp for the pooled estimate.
It fails where a segment's identity lies below the grid (halves at 99.5% and 85%: +1.7 pp), because the posterior
cannot put mass below the lowest level.
.ne 2.0i
.TC 13
Mosaic contigs (10 kb, two halves, 200 per row) against the random reference: true ANI, power mean (20), aif's
pooled estimate, aif's reported ANI (heterogeneity test and posterior mean), and the share called heterogeneous.
.TCE
.ps 9
.vs 11p
.so gen/tab_het.tbl
.ps \n(PS
.vs \n(VSp
