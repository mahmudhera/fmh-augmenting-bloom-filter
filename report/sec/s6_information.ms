.H1 "Information per stored bit: what to test"
.LP
Section 4 maximised the \f[I]number\f[R] of probes per bit. This section asks how much each probe tells us,
and whether something other than a $k$-mer should be tested.
.H2 "Choosing k"
.TH "Proposition 7" "optimal test length"
In the sparse regime a $k$-mer test is Bernoulli($e sup { - lambda k }$). Its Fisher information about $lambda$ is
.DS
	$I(k) ~=~ { k sup 2 e sup { - lambda k } } over { 1 - e sup { - lambda k } } ~=~ { g( lambda k ) } over { lambda sup 2 } ,~~~~~~ g(u) = { u sup 2 } over { e sup u - 1 } .$	(14)
.DE
The cost per test does not depend on $k$, so information per stored bit is maximised at $k sup * = u sup * / lambda$.
Here $u sup * approx 1.594$ is the positive root of $2 (1 - e sup -u ) = u$, and $g( u sup * ) approx 0.648$.
.TQ
.TC 16
Optimal $k$ and the information retained by $k=21$, together with the best admissible $k$ when spurious matches must
stay below $10 sup -3$ (global probes on 5 Mbp: $k >= 17$; island-local probes with $l = 1000$: $k >= 10$).
\(lqrel. info\(rq is $g( lambda k ) / g( u sup * )$.
.TCE
.so gen/tab_optk.tbl
.LP
$k = 21$ is within 8% of optimal at 90\(en95% ANI, so the standard choice is right for species-level work.
Elsewhere it is not.
At 99% it keeps only 29% of the attainable information per probe, and at 99.9% only 3%: almost every 21-mer
survives and the test says little. Section 6.2 addresses this.
At 80% it keeps 32%: almost every 21-mer is destroyed. The fix, $k approx 7$\(en10, is ruled out for global
probes by spurious matches ($k >= 17$), but island-local probes allow $k >= 10$.
\f[B]Localisation (D7) is worth 2.9\(mu more information per probe at 80% ANI and 5.7\(mu at 75%.\f[R]
A smaller $k$ also raises $theta$, which raises the admissible FPR in Theorem 3, a second gain.
The dense floor improves as well (Table 17): at 80% ANI its efficiency rises from 0.12 at $k=21$ to 0.46 at $k=10$.
Tiers at two values of $k$ (for example 21 and 12) could be stored and chosen at query time, in the spirit of
multi-resolution containment [36]. Small $k$ requires masking low-complexity sequence, where the
no-repeat assumption fails.
.TC 17
Efficiency of the containment floor ($s = 1$, $eps = 0$) as a function of $k$.
.TCE
.so gen/tab_floork.tbl
.H2 "Island fingerprints: testing the island itself"
.LP
Proposition 7 says that at high ANI the ideal test is \f[I]long\f[R], hundreds to thousands of bases.
Islands already provide long, content-defined, position-anchored units.
Design D6 stores, for each island, an $f$-bit hash $phi sub j$ of its interior: the $l sub j = pi sub {j+1} - pi sub j - k sub a$
bases between the two anchor $k$-mers.
A query with consecutive matched anchors ($j prime = j + 1$) and $DELTA = 0$ hashes its own interior and compares.
The two are identical with probability $e sup { - lambda l sub j }$; a fingerprint collision (probability $2 sup -f$)
acts as a false positive.
.TH "Proposition 8" "island fingerprints"
An island test carries Fisher information $l sup 2 e sup { - lambda l } / (1 - e sup { - lambda l } )$ about $lambda$ (for $f -> inf$).
Per stored bit (anchor plus fingerprint) this is maximised at mean island length $l sup * = u sup * / lambda$.
With anchors at rate $s sub a$ and both anchors surviving with probability $theta sub a sup 2$,
.DS
	$eta sub D6 ~=~ p over {1-p} ~ s sub a left [ I( k sub a ) ~+~ theta sub a sup 2 ~ Ex sub l left ( { l sup 2 e sup { - lambda l } } over { 1 - e sup { - lambda l } } right ) right ] ,$	(15)
.DE
where the expectation is over geometric island lengths. The first term is the anchors' own survival tests.
.TQ
.TC 18
Island fingerprints (D6, $k sub a = 21$, $f = 8$, anchors as 32-bit fingerprints plus positions) versus D0 and D2.
Columns 2\(en4 are at the D0 budget; the last three columns tune the island length for information per bit.
\(lqper-bit gain\(rq compares $eta / beta$ with D2's.
.TCE
.so gen/tab_islandfp.tbl
.LP
At 99.9% ANI and the D0 budget, D6 reaches $eta = 0.46$, against 0.33 for the tuned Bloom tier (D2) and
0.02 for plain FMH (D0). Tuned for information per bit, D6 at 99.95% gives $eta = 0.25$ for only 0.017 bits per
base, 2.9\(mu D2's efficiency per bit.
At 99.5% and below, islands are rarely identical and D6 loses to D2.
\f[B]D6 is the natural generalisation of the island idea: instead of testing the island's $k$-mers, test the
island, i.e. a position-anchored \(lq$k$-mer\(rq whose length adapts to the sequence.\f[R]
It is the design of choice at strain level (ANI $>=$ 99.9%), where 21-mer tests are nearly always positive.
A multi-resolution (Merkle [35]) variant stores fingerprints of each island and its halves, quarters, and so on
down to length $l sub 0$. It costs about $2 f / l sub 0$ bits per base, keeps near-optimal units for every ANI
above about $1 - 1.6 / l sub 0$, and with polynomial rolling hashes the query computes every level in one pass.
Sequencing or assembly errors make islands differ exactly as mutations do, so D6 measures the combined
per-base difference.
.H2 "Dense positional match vectors"
.LP
If the query probes \f[I]every\f[R] $k$-mer (D11), the hit vector $Z sub 1 ... Z sub n$ has run structure.
An isolated substitution produces a run of exactly $k$ misses. Substitutions closer than $k+1$ merge, and a run
of length $r$ contains at least $\(lc r/k \(rc$ substitutions. When $kp << 1$, the number of runs equals the number
of substitutions with high probability, so the realised mutation count is recovered almost exactly.
False positives are repairable: an interior run shorter than $k$ is impossible without one, so an isolated hit
that splits a run can be detected and removed.
D11 is not a sketch, however. By Theorem 2 it needs at least 2 bits per base, so one might as well store the
sequence. And the containment floor is already nearly efficient when $kp << 1$ (Tables 4 and 17), which caps
the gain at $1/ eta sub floor$: 1.4\(mu at 95% and 2.3\(mu at 90% for $k=21$.
D11 belongs in a second stage, applied to retrieved candidate regions, not in the sketch.
.H2 "Rate heterogeneity"
.TH "Proposition 9" "power-mean bias"
Let aligned regions have rates $p sub i$, weighted by length. Containment converges to $Ex [ (1 - p sub i ) sup k ]$, and
.DS
	$Chat sup {1/k} ~ -> ~ left ( Ex [ (1 - p sub i ) sup k ] right ) sup {1/k} ~ >= ~ 1 - Ex [ p sub i ] ,~~~~~~ roman bias ~ approx ~ { (k-1) Var ( p sub i ) } over { 2 (1 - p bar ) } .$	(16)
.DE
.TQ
.LP
The inequality is the power-mean inequality and the bias follows from a second-order expansion.
In Table 19 the bias at a true mean ANI of 95% ranges from +0.6 to +2.2 pp. A genome with half its aligned
length recently transferred (identical) and half at 90% identity is reported as 97.2%. Larger $k$ is worse.
Estimating ANI per window along the query, transforming each window to ANI and averaging by length removes
this bias when windows are homogeneous. It needs query-side positions only.
The price is the per-window small-sample bias of (8), which grows like $1/(W s sub 2 )$; windows with at least
30\(en50 probes ($W >= 40 / s sub 2$, about 2 kb at $s sub 2 = 1/50$) balance the two.
Within a short contig heterogeneity is limited, so this matters mainly for genome-scale ANI.
FastANI's fragment-based identity [25] is the alignment analogue.
.TC 19
ANI overestimation by containment under heterogeneous rates (true mean ANI 95%).
.TCE
.so gen/tab_jensen.tbl
