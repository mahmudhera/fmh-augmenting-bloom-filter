.H1 "Positions, islands and alignment fraction"
.H2 "Anchors and islands"
.LP
Let the reference store anchors $"{" ( f(x), pi sub x ) : x \(mo K sub B sup {( k sub a )} , h sub a (x) < s sub a "}"$.
Here $f$ is a hash or fingerprint and $pi sub x$ a position, sorted by position.
Consecutive anchors $pi sub j < pi sub {j+1}$ delimit island $I sub j$: the $k$-mers starting strictly between them.
Under ideal hashing, selection events of distinct $k$-mers are independent, so anchor starts form a
Bernoulli($s sub a$) process. Gaps are geometric with mean $1/s sub a$ and coefficient of variation close to 1,
so long islands are common: $Pr ( roman gap > 3/ s sub a ) approx e sup -3 approx 5%$.
Islands are \f[I]content-defined\f[R]: their boundaries depend only on local sequence, so an indel shifts
coordinates but not island identities. This is the principle behind content-defined chunking in data
deduplication [34].
.LP
A reference anchor survives in $A$ if its $k sub a$-mer is unmutated, with probability $theta sub a = (1-p) sup {k sub a}$.
.TH "Proposition 4" "anchoring feasibility"
Matched anchors in the aligned part of $A$ form, approximately, a Bernoulli process of rate $s sub a theta sub a$
per base, and their number is approximately Poisson($alpha L sub A s sub a theta sub a$).
Having at least $m sub 0$ matched anchors with probability at least $1 - delta$ requires
$alpha L sub A >= lambda ( m sub 0 , delta ) / ( s sub a theta sub a )$, with $lambda (3, 0.05) = 6.30$.
.TQ
.TC 12
Minimum contig length (kb) for at least 3 matched anchors with 95% probability ($alpha = 1$), anchor storage
cost, and the global probability that a random anchor $k sub a$-mer occurs somewhere in a 5 Mbp genome.
.TCE
.so gen/tab_anchor.tbl
.LP
Default FMH anchors ($s sub a = 1/1000$, $k sub a = 21$) need 18.5 kb at 95% ANI and 57 kb at 90%, so they cannot
anchor short contigs. Denser, shorter anchors ($s sub a = 1/200$, $k sub a = 15$) bring this down to 2.7 kb at
95% and 6.1 kb at 90%, for 0.21 bits per base. At 80% even these need 36 kb.
Pushing to $k sub a = 12$ makes 30% of random anchors match somewhere, so chaining must reject many stray
anchors. \f[B]Anchoring short contigs needs denser and shorter anchors than FMH's defaults; short contigs at
low ANI are fundamentally hard for any anchor-based method.\f[R]
.H2 "Why containment alone cannot give alignment fraction"
.TH "Proposition 3" "non-identifiability"
Let a fraction $alpha$ of $A$ be homologous to $B$ (SMM, rate $p$), with the rest absent from $B$ and no
spurious matches. In the sparse limit (selected $k$-mers more than $k$ apart),
$Y | N ~ \(ap ~ roman Bin (N, (1- eps ) alpha theta (p) + eps )$.
The likelihood of any unordered sketch comparison therefore depends on $( alpha , p)$ only through
$alpha theta (p)$. ANI and AF are not separately identifiable, and
.DS
	$Chat sup {1/k} ~ -> ~ alpha sup {1/k} (1-p) .$	(13)
.DE
.TQ
.LP
\f[I]Proof.\f[R] Each selected $k$-mer independently falls in the homologous part with probability $alpha$ and
is then shared with probability $theta$. \(sq
.LP
The bias is large (Table 13). At $alpha = 0.5$ a true 95% is reported as 91.9%, which crosses the species
boundary. Larger $k$ reduces the bias but does not remove it.
This is the formal content of the observation that FMH without positions cannot estimate AF.
.TC 13
ANI reported by the naive containment estimator when only a fraction $alpha$ of the contig aligns
(bias in pp in parentheses).
.TCE
.so gen/tab_afbias.tbl
.H2 "Alignment fraction from positions"
.LP
There are two routes.
\f[B]Reference-side anchors (D3).\f[R] Chain matched anchors colinearly (in both orientations) and take each
chain's span from its first to its last matched anchor. Each boundary is then misplaced by about $1/(s sub a theta sub a )$ bases.
\f[B]Query-side coordinates (D12).\f[R] AF is a property of the query's coordinates, which suggests the following.
.TH "Proposition 5" "AF needs query positions only"
(a) The AF of $A$ in $B$ depends only on which positions of $A$ have homologs in $B$. It can therefore be
estimated from the query's own tier-2 positions and their hit indicators, with no reference positions.
(b) In an aligned stretch, hits form approximately a Bernoulli process of rate $s sub 2 theta sub 1$ per base,
with $theta sub 1 = theta + (1- theta ) eps$. Outside, the rate is $s sub 2 theta sub 0$ with $theta sub 0 = eps + sigma$.
If $theta sub 0 << theta sub 1$, the distance from a true boundary to the nearest hit on its aligned side is
geometric with mean about $1/(s sub 2 theta sub 1 )$, and a two-state segmentation (HMM or change point)
locates boundaries to within this order.
.TQ
.TC 14
Mean boundary localisation error, $1/(s theta )$ bases ($k=21$). Rows with $s = s sub a$ apply to chain ends,
and rows with $s = s sub 2$ to tier-2 refinement.
.TCE
.so gen/tab_afres.tbl
.LP
For a 5 kb contig at 95% ANI, anchors or $k$-mers sampled at 1/1000 place a boundary to within about 2.9 kb,
which is useless. Tier-2 sampling at 1/50 gives 147 bp (about 3% of AF per boundary), and 1/20 gives 59 bp.
The two routes combine naturally: chain with anchors, then extend each chain end by probing tier-2
$k$-mers against the flanking islands. Query-side segmentation alone has three limitations.
It cannot check colinearity, so chimeric or rearranged contigs look fully aligned (Scenario S11).
It cannot give the AF of $B$ with respect to $A$.
And it needs $theta sub 1 >> theta sub 0$: at 80% ANI, $theta approx 0.009$, so $eps$ and $sigma$ must be
below about $10 sup -3$, and boundaries blur to about 5.5 kb at $s sub 2 = 1/50$.
.H2 "ANI over the aligned span"
.LP
Restricting (5) to probes inside the estimated span gives $Chat sub aln$, which removes the $alpha sup {1/k}$
factor of (13). Theorem 1 applies with $n$ replaced by $alpha n$.
Boundary handling adds a relative bias of order $1/( N sub aln theta )$ in $C$, because the stretch between the
last hit and the true boundary consists of misses that are dropped. For the 5 kb example this is about 0.1 pp
of ANI, and an HMM posterior or an explicit correction removes it.
This bias is negligible next to the 3 pp from ignoring AF at $alpha = 0.5$.
.H2 "Island-keyed filters and localisation"
.LP
The proposal inserts each island $k$-mer $x$ under the key $(x, j)$, where $j$ is its island index (D4/D5).
To use the key, the query needs $j$.
Suppose consecutive matched anchors pair $a sub t$ with $b sub j$ and $a sub {t+1}$ with $b sub {j prime}$.
Every query $k$-mer strictly between $a sub t$ and $a sub {t+1}$ maps between $b sub j$ and $b sub {j prime}$,
that is, into islands $j, ..., j prime - 1$.
If $j prime = j+1$ (no anchor lost to mutation), \f[I]the island is determined exactly, whatever indels it contains\f[R].
If $j prime > j + 1$, the query interpolates coordinates and probes the one or two islands its position maps to;
probing all $j prime - j$ islands would multiply the FPR by that number.
Between consecutive matched anchors, an average of $(1 - theta sub a )/ theta sub a$ reference anchors are lost:
1.9 at 95% ANI with $k sub a = 21$ and 1.2 with $k sub a = 15$.
.TH "Proposition 6" "localisation"
Model a non-homologous probe as a uniformly random $k$-mer. Its spurious-match probability is
$sigma sub g = 1 - (1 - 4 sup -k ) sup { L sub B } approx L sub B / 4 sup k$ under a global key, and
$sigma sub l approx l sub j / 4 sup k$ under island key $j$ (island length $l sub j$).
Localisation divides $sigma$ by about $L sub B / l approx s sub a L sub B$, the number of islands.
The per-probe FPR $eps$ is unchanged.
.TQ
.TC 15
Spurious-match probability per probe: island-local versus global keys (uniform random sequence model;
real genomes are worse because of low-complexity sequence and repeats).
.TCE
.so gen/tab_spur.tbl
.LP
Localisation has four consequences.
(i) \f[B]Small $k$ becomes admissible.\f[R] $sigma <= 10 sup -3$ needs $k >= 17$ with global keys for a 5 Mbp genome
($k >= 21$ for 3 Gbp), but only $k >= 10$ with island keys ($l = 1000$). Section 6.1 shows this matters at low ANI (D7).
(ii) \f[B]Repeats and paralogs are suppressed.\f[R] Global probes hit any copy; island-keyed probes hit only the
chained locus. This matters because assemblies break at repeats, so short contigs are enriched near them.
(iii) \f[B]The reversed asymmetry is rescued\f[R] (Table 11, last column).
(iv) \f[B]No extra space.\f[R] Keys change but the element count does not. The cost is the anchors, plus
the need to match them first.
The one new failure mode is a false negative from mis-localisation. It can occur only when an anchor is
lost and the interpolated position falls on the wrong side of an island boundary. Probing both islands
within the local indel uncertainty of a boundary avoids it, at twice the FPR for those few probes.
.H2 "Indels"
.LP
For consecutive matched anchors, the span discrepancy
$DELTA = ( pi sub {t+1} sup A - pi sub t sup A ) - ( pi sub {j prime} sup B - pi sub j sup B )$
is the net indel length between them, which is information an unordered sketch cannot see.
A deletion from $A$ relative to $B$ removes $k - 1$ shared $k$-mers, and an insertion of $d$ bases removes $k+d-1$.
Containment therefore counts each indel event as about $1 + (d-1)/k$ substitutions. Because ANI is
conventionally computed over aligned columns, $k$-mer ANI is biased low by about the indel event rate $rho$.
A single event always gives $DELTA != 0$; cancellation needs two or more events of opposite net length.
So the fraction of chained spans with $DELTA != 0$ estimates $1 - e sup { - rho l bar }$, which yields $rho hat$.
The tier-2 probes inside those spans can then be down-weighted, or the expected number of indel-killed
$k$-mers removed from the denominator of $Chat$.
Island keys are unaffected by indels within an island. Fixed-offset keys such as $(x, \(lf pi / W \(rf )$ would
break, which is why content-defined islands are the right granularity.
.H2 "Choice of anchors"
.LP
FMH anchors have geometric gaps, and a single long island can hold most of a short contig.
Anchors only \f[I]localise\f[R]; they need not give unbiased set estimates.
They can therefore come from a window-guaranteed scheme (D13), such as minimisers [28, 29],
syncmers [30] (better conserved under mutation) or strobemers [31] (robust to indels), which bound gap lengths.
The estimation tier keeps FMH's consistency property (4), which is what makes $Chat$ unbiased.
Decoupling \(lqwhere\(rq from \(lqhow much\(rq this way removes the bias that minimiser-based Jaccard estimators
suffer [33], while keeping their spacing guarantees.
