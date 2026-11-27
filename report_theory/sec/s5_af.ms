.H1 "Alignment fraction"
.H2 "Counts cannot separate AF from ANI"
.TH "Proposition 6" "non-identifiability [first theory report]"
Let a fraction $alpha$ of the contig be aligned at containment $theta$ and the rest unaligned. In the sparse limit
(probes more than $k$ apart) each probe hits independently with probability
$eps prime + (1 - eps prime ) alpha theta$, so the hit counts depend on $( alpha , theta )$ only through
$alpha theta$. The containment estimator converges to $alpha theta$, and its ANI to $alpha sup {1/k} a$.
.TQ
.LP
At $alpha = 0.5$ and $k = 20$ a true 95% is read as 91.8%. AF is identifiable only from the \f[I]positions\f[R]
of the hits along the contig: aligned stretches have hits at rate about $s theta$ per base, unaligned ones at
$s eps prime$. Positions on the reference are not needed.
.H2 "The segment HMM"
.LP
aif reads the probes in contig order with a hidden Markov model [11]. The hidden state at a probe is $U$ (unaligned) or $A sub j$ (aligned at
identity $a sub j$). Emissions are independent given the state: a tier-$t$ probe hits with probability
$eps prime sub t$ in $U$ and $pi sub {tj} = theta sub j + (1 - theta sub j ) eps prime sub t$ in $A sub j$, with
$theta sub j = a sub j sup k$. Between probes $d$ bases apart the state changes as follows:
.IP \(bu 0.2i
from $U$, it moves to aligned with probability $1 - e sup { - rho sub u d }$, landing on a level drawn uniformly
from the $m$ levels;
.IP \(bu 0.2i
from $A sub i$, it moves to $U$ with probability $1 - e sup { - rho sub u d }$; otherwise it redraws its level
uniformly with probability $1 - e sup { - rho sub h d }$ and keeps it otherwise.
.LP
The initial state is $U$ or aligned with probability 1/2 each. The level grid runs from 0.9995 down in steps of
0.005 to the identity at which $theta sub j$ falls to $8 min sub t eps prime sub t$; with the default tiers at
$k = 20$ this gives $m = 28$ levels down to 0.8645. Two implementation details do not change the model: the
forward\(enbackward pass runs over bins of 32 probes (the state is constant within a bin), and bins whose
aligned posterior is between 0.05 and 0.95 are recomputed probe by probe.
.LP
From the posterior probability $gamma (x)$ that $k$-mer start $x$ is aligned, AF is
.DS
	$roman AF hat ~=~ 1 over L left [ int gamma (x) ~ dx ~+~ (k - 1) R right ] ,$	(15)
.DE
where the integral interpolates $gamma$ between probes and $R$ is the number of aligned runs (each run of
$k$-mer starts covers $k - 1$ more bases). Containment is (10) with every probe weighted by its aligned
posterior. The contig is reported when $LR$, the ratio of the HMM's marginal likelihood to that of "all
unaligned", is at least $e sup 10$.
.H2 "How precisely a boundary is placed"
.LP
Consider one boundary, from aligned at containment $theta$ to unaligned, with many probes on each side. Per probe,
the aligned side is distinguished from the unaligned side by the divergences
$D sub A = sum sub t w sub t KL ( pi sub t || eps prime sub t )$ and $D sub U = sum sub t w sub t KL ( eps prime sub t || pi sub t )$,
where $w sub t = s sub t / s$ is the tier's share of probes.
.TH "Proposition 7" "change-point resolution"
With independent probes and a known level, the maximum-likelihood boundary is the argmax of a two-sided random walk
whose steps have drift $- D sub A$ on the aligned side and $- D sub U$ on the unaligned side (Hinkley [12]).
Its error, and the spread of the posterior, are of order
.DS
	$Delta ~ approx ~ 1 over s left ( 1 over { D sub A } + 1 over { D sub U } right ) ~~ roman bases .$	(16)
.DE
.TQ
.LP
At 95% ANI, $D sub A = 0.67$ and $D sub U = 0.32$ nats per probe, so $Delta approx 38$ bp; at 90%, 181 bp; at 99%,
9 bp. The exact posterior of one change point (Monte Carlo, 4,000 boundaries per identity) has error standard
deviations of 39, 196 and 8 bp, and mean posterior standard deviations about 10% smaller (Table 8): in the
idealised model the boundary posterior is both sharp and calibrated.
.LP
aif does not reach this resolution. Table 8 decomposes the gap with an exact per-probe forward\(enbackward of the
full model (no bins) on probes simulated under the substitution model, 300 contigs of 5 kb aligned + 5 kb
unrelated sequence per identity:
.IP (i) 0.3i
\f[B]Runs of misses.\f[R] One substitution kills $k$ consecutive $k$-mers, about 2.4 probes, so misses come in
runs and the boundary is less certain: with the true level known, the error SD is 86 bp at 95% ANI and 36 bp at
99% (against 39 and 8 bp for independent probes).
.IP (ii) 0.3i
\f[B]Low levels near the boundary.\f[R] With all 28 levels, the stretch just past the boundary can be explained as
aligned at a low level instead of unaligned, at almost no cost (Section 5.4). This adds a positive tail: SD 147 bp
at 95% and 113 bp at 99% when AF counts all aligned mass (v1.6).
.IP (iii) 0.3i
\f[B]The v1.7 floor\f[R] (count only levels within 0.05 of the identity) removes that tail but also removes the
end of the true alignment whenever a chance run of misses makes it look like a low level: mean \(mi113 bp at 95%
without aif's per-window evidence rule. aif's rule keeps most of it (mean \(mi10 bp).
.LP
aif's own errors (Table 8, last two columns) match the exact full model, and they do not change when the bins are
reduced to one probe (SD 171 bp at 95% with single-probe bins, 163 bp with 32). \f[B]AF resolution is limited by the
model, not by its binned implementation\f[R], and mostly by the low identity levels, which the next section
analyses.
.ne 2.2i
.TC 8
AF error at one boundary (5 kb aligned + 5 kb unrelated), in bp: mean / SD. \(lqindependent probes\(rq: exact
posterior of one change point, known level. \(lqsubstitution model\(rq: exact per-probe forward\(enbackward on
probes simulated with overlapping $k$-mers (two states with the true level; all 28 levels with AF from all
aligned mass; all levels with mass counted only within 0.05 of the true identity). aif: 400 contigs per row,
against the random reference. At 90% the floor lies below the grid, so the floor rule changes nothing.
.TCE
.ps 9
.vs 11p
.so gen/tab_boundary.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Low levels absorb unaligned stretches"
.LP
The lowest levels expect hits barely more often than the FPRs: at $k = 20$ the lowest level, 86.45%, has
$theta = 0.054$, so a tier-2 probe hits with probability 0.150 there against 0.102 when unaligned. Unaligned
sequence explained as aligned at such a level costs little in emissions. What decides between the two
explanations is mostly the switching prior.
.TH "Theorem 3" "absorption"
Let a stretch of $G$ bases with no homology lie either between two aligned stretches at a common level (an
\f[I]internal gap\f[R]) or at a contig end next to one (a \f[I]flank\f[R]). In the limit of small switching
probabilities per probe gap, the posterior odds that the stretch is aligned at some level $j$ throughout rather
than unaligned throughout are
.DS
	$O(G) ~=~ kappa ~ e sup { - rho sub h G } sum from j e sup { - LAMBDA sub j } ,~~~~~~ kappa sub { roman "int" } = { rho sub h sup 2 } over { m rho sub u sup 2 } ,~~~~ kappa sub { roman "edge" } = { rho sub h } over { m rho sub u } ,$	(17)
.DE
where $LAMBDA sub j$ is the log-likelihood ratio of $U$ against $A sub j$ on the stretch's probes. Under unaligned
data $LAMBDA sub j / G -> r sub j = sum sub t s sub t KL ( eps prime sub t || pi sub {tj} )$, while
$Ex [ e sup { - LAMBDA sub j } ] = 1$.
.TQ
.LP
\f[I]Proof.\f[R] Compare the two path families, which agree outside the stretch. Leaving level $j sub 0$ for $U$
costs $1 - e sup { - rho sub u d } approx rho sub u d$, and returning costs $rho sub u d / m$ because the level is
redrawn. Moving to level $j$ costs $rho sub h d / m$, and returning to $j sub 0$ the same. Staying costs
$e sup { - rho sub u d }$ per step in $U$ and $e sup { - ( rho sub u + rho sub h ) d }$ in $A sub j$. The ratio
of prior path probabilities is therefore $kappa sub { roman "int" } e sup { - rho sub h G }$ for an internal gap. At a flank only
one switch is needed (or, at the contig start, the initial state replaces one), giving $kappa sub { roman "edge" }$. The
emission ratio is $e sup { - LAMBDA sub j }$; its expectation under $U$ is 1 because it is a likelihood ratio. \(sq
.LP
With the defaults, $kappa sub { roman "int" } = 40$ and $kappa sub { roman "edge" } = 1.2$: \f[B]the prior favours "aligned at a low
level" over "unaligned" for an internal gap\f[R], because one $U$ excursion costs two switches at rate $rho sub u$
plus a level redraw, while a level excursion costs two redraws at rate $rho sub h = 33 rho sub u$ spread over 28
levels. Only the emissions push back, at $r sub j$ = 2.7 nats per kb for the lowest level and 8.2 for the eighth
lowest (Table 9). The typical odds fall to 1 at $G sup * = 1.5$ kb for internal gaps and 0.5 kb at contig ends.
(The mean odds do not fall, since $Ex [ e sup { - LAMBDA } ] = 1$: some long gaps are absorbed by chance.)
.LP
Figure 1 compares (17) with aif. With AF counting all aligned mass (v1.6), aif counts 86% of a 1 kb internal gap
as aligned against 85% predicted, 63% of a 1.5 kb gap against 51%, and 14% of a 3 kb gap against 1%; at contig ends,
57% of a 500 bp flank against 51%. The theory reproduces the break-even lengths and slightly under-predicts the
tail, as expected from ignoring partial-absorption paths. This mechanism is the +5 pp AF bias found on real
contigs (gut contigs against GTDB, aif v1.6): there, three quarters of the excess came from bridged internal gaps
and one quarter from contig ends, in line with internal gaps being absorbed up to three times longer.
.ne 1.6i
.TC 9
Emission cost of explaining unaligned sequence by the six lowest levels ($k = 20$, default tiers), and the gap
length at which each level alone stops being favoured over $U$, $ln kappa / ( r sub j + rho sub h )$.
.TCE
.ps 9
.vs 11p
.so gen/tab_absorb.tbl
.ps \n(PS
.vs \n(VSp
.ne 2.9i
.PSPIC -C figs/fig_absorb.eps 6.4i
.FC 1
Share of an unaligned stretch that is counted as aligned: Theorem 3 (dashed; typical odds with $LAMBDA sub j$ at
its mean) and aif on 200 contigs per point, with AF from all aligned mass (v1.6) and with v1.7's local-identity AF.
Aligned parts are at 97% ANI. Negative values: true aligned bases removed.
.FE
.H2 "What the local-identity AF of v1.7 can remove"
.LP
v1.7 keeps the HMM and changes what AF counts: aligned posterior mass at levels within $delta = 0.05$ of the
\f[I]local identity\f[R], the probe-weighted median identity of the aligned-looking windows of 64 probes within
three windows on either side. A window whose best level beats $U$ by at least 3 nats on its own lowers its floor
to its own identity minus $delta$, as do its two neighbours, so that genuinely divergent blocks keep their AF.
.TH "Proposition 8" "scope of the local-identity AF"
(a) For a homogeneous alignment at identity $a >= a sub m + delta$, the level posterior concentrates at $a$ as
the aligned length grows, so the mass below $a - delta$ vanishes and AF is unchanged in the limit.
(b) An absorbed stretch carries its mass on the levels with the smallest $r sub j$, the lowest ones. They lie
below the floor whenever the local identity exceeds $a sub m + delta$ (0.915 at $k = 20$), so absorbed mass is
removed unless the stretch's window has its own 3 nats of evidence.
(c) A stretch that shares a window, or the neighbour of one, with its aligned flanks inherits their evidence and
stays counted. Conversely, a divergent block more than $delta$ below the local identity loses its AF if it is
too short to reach 3 nats on its own.
.TQ
.LP
Figure 1 shows the reach of (c): aif v1.7 still counts 82% of a 100 bp internal gap, 55% of 500 bp and 31% of
1 kb (v1.6: 100%, 99% and 86%), removes nearly all of gaps from 1.5 kb, and slightly over-corrects at 2 kb
(\(mi7% of the gap, 145 bp of true alignment). At contig ends the residual is 68% of 100 bp and 6% of 500 bp. The
correction is therefore complete for long insertions and partial for short ones, whose bases are few. For a
divergent block, 3 nats need about $3 / r(a)$ bases of own evidence: 360 bp at 88% ANI and 740 bp at 86% (Table
11), so blocks near the bottom of the grid inside a high-identity contig can lose AF.
.H2 "Identity-thresholded AF"
.LP
With \f[CR]--af-min-ident\f[R] $T$, AF is multiplied by the share of aligned posterior mass at levels
$>= T$. For homogeneous segments the level posterior concentrates at each segment's identity, so the estimate
converges to the length share of segments with $a sub i >= T$, the AF that a seed-and-extend aligner with an
identity floor reports. Being a global share, it assumes each segment's level posterior is concentrated; a contig
whose segments straddle $T$ within the posterior's spread is split between them.
