.H1 "Model, notation and building blocks"
.H2 "Sequences, k-mers and the simple mutation model"
.TC 1
Notation.
.TCE
.TS
center box tab(@);
l l.
\f[B]symbol\f[R]@\f[B]meaning\f[R]
_
$A, B, L sub A , L sub B , n$@query (contig), reference, lengths; $n = L sub A - k + 1$
$k, k sub a , k sub e$@$k$-mer length; anchor $k$; estimation $k$ (default $k sub a = k sub e = k = 21$)
$p, lambda , theta , q$@substitution rate; $lambda = - ln (1-p)$; $theta = (1-p) sup k$; $q = 1- theta$
$C, C sub M , Chat$@containment of $A$ in $B$; its realised value; estimator
$s, s sub a , s sub 2$@FMH scale (baseline); anchor rate; tier-2 (island) sampling rate
$eps , c( eps ), h$@AMQ false-positive rate; bits per element; number of hash functions
$beta$@sketch bits per reference base ($beta sub 0 = 64 times 1/1000 = 0.064$)
$alpha$@alignment fraction (AF) of $A$ in $B$
$sigma$@probability of a spurious (non-homologous) match per probe
$eta$@efficiency, i.e. alignment-oracle variance divided by estimator variance
$B(p,k)$, $phi$@process-variance coefficient (3); $phi = eps / ((1- eps ) theta )$
$J( eps )$@bits per unit of information, $c( eps )(1+ phi )$ (Theorem 3)
.TE
.LP
Let $A$ (the query, typically a contig) and $B$ (the reference) have lengths $L sub A$ and $L sub B$.
Let $K sub A$ and $K sub B$ be their $k$-mer sets, and $n = L sub A - k + 1$.
Following [7, 6] we assume no repeated $k$-mers and, unless stated otherwise, no
spurious matches; Section 5.5 removes the second assumption.
In the \f[I]simple mutation model\f[R] (SMM), the homologous part of $B$ arises from $A$ by
independent substitutions at rate $p$ per base. We write
.DS
	$theta = (1-p) sup k ,~~~~ q = 1 - theta ,~~~~ lambda = - ln (1-p),~~~~ roman ANI = 1-p .$	(1)
.DE
A $k$-mer of $A$ is shared with $B$ (\(lqunmutated\(rq) with probability $theta$.
The realised containment $C sub M = 1 - N sub mut / n$, where $N sub mut$ is the number of mutated $k$-mers,
has $Ex [ C sub M ] = theta$, and its variance is known exactly [7]:
.DS
	$Var ( N sub mut ) = n theta B(p,k) + k(k-1) theta sup 2 + { 2 theta } over { p sup 2 } left [ (1 + (k-1) theta ) p - q right ] ,$	(2)
.DE
.DS
	$B(p,k) ~=~ q (2k - 1 + 2/p ) - 2k ~~ approx ~~ k sup 2 p ~~~ ( kp << 1 ).$	(3)
.DE
We call $theta B(p,k)/n$ the \f[I]process variance\f[R] of $C sub M$. It is a property of the mutation
process and is present even if every $k$-mer is compared.
.H2 "FracMinHash and the consistency property"
.LP
Let $h$ map $k$-mers to $[0,1)$, idealised as uniform and independent.
$roman FMH sub s (K) = "{" x \(mo K : h(x) < s "}"$.
The FMH containment estimator is $Chat sub s = | roman FMH sub s (K sub A ) \(ca roman FMH sub s (K sub B )| / | roman FMH sub s (K sub A ) |$.
The property everything below relies on is \f[I]consistency\f[R]: a shared $k$-mer is selected in
both sequences or in neither. Consequently, for a query $k$-mer $x$ with $h(x) < s$,
.DS
	$x \(mo K sub B ~~ <=> ~~ x \(mo roman FMH sub s (K sub B ) ,$	(4)
.DE
so testing membership in the \f[I]sampled\f[R] reference set is as good as testing membership in all of $K sub B$.
Equation (4) is what makes it legitimate to store only a subsample of $B$ in a filter and still
obtain an unbiased containment estimate. Minimiser/winnowing schemes [28, 29] lack it,
which is why their Jaccard estimators are biased [33].
.H2 "Approximate membership query (AMQ) structures"
.LP
A Bloom filter with $m$ bits and $h$ hash functions storing $N$ items has FPR
$eps approx (1 - e sup { - hN/m } ) sup h$ and no false negatives [1, 13].
Its cost is
$c sub BF ( eps ) = m/N = - h / ln (1 - eps sup {1/h} )$ bits per element; minimised over real $h$
this is $log sub 2 e cdot log sub 2 (1/ eps ) approx 1.44 log sub 2 (1/ eps )$.
Two facts matter for estimation.
First, the FPR of a \f[I]given\f[R] filter is known exactly: if a fraction $rho$ of its bits are set, a
non-member is reported present with probability $rho sup h$, independently across non-members
(for ideal hashes).
Second, any AMQ needs at least $N log sub 2 (1/ eps )$ bits [12].
Other AMQs come closer to this bound: Golomb-coded sets [15], xor and binary fuse filters [16, 17],
ribbon filters [40], and quotient or cuckoo filters [41, 18].
Note also that a FracMinHash sketch whose hashes are truncated to $b$ bits is itself an AMQ:
a non-member collides with one of the $N$ stored fingerprints with probability $eps approx N 2 sup -b$,
and the sorted fingerprints can be Elias\(enFano coded in about $log sub 2 (1/ eps ) + 2$ bits each.
Table 2 lists the costs used throughout.
.TC 2
Bits per stored element $c( eps )$ for AMQ structures at FPR $eps$.
\(lqBF (best integer $h$)\(rq is the exact Bloom cost minimised over integer $h >= 1$;
binary fuse filters [17] require $eps = 2 sup -f$.
The last column is the factor by which a BF stores more elements than a 64-bit FMH in the same space.
.TCE
.so gen/tab_cost.tbl
.H2 "Designs studied"
.LP
Table 3 fixes names for the 14 designs analysed. D4 and D5 are the proposal in its dense
and subsampled forms; the others are the close relatives that Q5 asks about.
.TC 3
Designs. \(lqQuery side\(rq states what the query must supply: a raw sequence can produce any of these.
.TCE
.TS
center box tab(@);
l lw(3.4i) l.
\f[B]id\f[R]@\f[B]reference stores\f[R]@\f[B]query side\f[R]
_
D0@FMH $s$ with 64-bit hashes (baseline; sourmash-like)@sketch
D1@D0 + dense global BF over all $k$-mers of $B$@all $k$-mers
D2@two-tier FMH: global AMQ (BF) over $roman FMH sub {s sub 2} (K sub B )$, $s sub 2 > s$@tier-2 $k$-mers
D3@positional FMH: anchors ($roman FMH sub {s sub a}$) with positions; chaining@anchors + positions
D4@D3 + BF over \f[I]all\f[R] island $k$-mers, keyed $(x, j)$ by island $j$ (dense proposal)@raw
D5@D3 + BF over island $k$-mers with $h(x) < s sub 2$, keyed $(x,j)$ (\f[B]AIF\f[R], subsampled proposal)@raw
D6@D3 + one $f$-bit fingerprint per island (content-defined chunk hash)@raw
D7@D5 with small estimation $k sub e$ (10\(en15), made possible by localisation@raw
D8@several tiers with different FPRs (hedging against unknown ANI)@tier $k$-mers
D9@$b$-bit (truncated, Elias\(enFano coded) FMH: a symmetric AMQ@sketch
D10@BF only, compared BF-vs-BF (symmetric)@BF
D11@dense positional match vector, run-length estimator@all $k$-mers
D12@D2 + query-side segmentation of hits (AF without reference positions)@raw
D13@window-guaranteed anchors (minimisers/syncmers) + FMH estimation tier@raw
.TE
