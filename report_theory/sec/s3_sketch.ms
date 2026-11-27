.H1 "The sketch"
.H2 "Tier layout"
.LP
A sketch has a budget of $beta$ bits per reference base and a list of tiers, each with a target FPR $eps sub t sup *$
and a share $f sub t$ of the budget (default: $eps sup * = 0.005$ with 60%, $eps sup * = 0.1$ with 40%). Tier
$t$ is a Bloom filter of $beta f sub t L sub B$ bits holding the $k$-mers whose selection hash lies in
$[ S sub {t-1} , S sub t )$, where
.DS
	$s sub t ~=~ S sub t - S sub {t-1} ~=~ { beta f sub t } over { c( eps sub t sup * ) } ,~~~~~~ c( eps ) ~=~ min from {h >= 1} ~ { - h } over { ln (1 - eps sup {1/h} ) } ,$	(3)
.DE
so that the filter holds $s sub t L sub B$ elements at $c( eps sub t sup * )$ bits each. The number of hash
functions $h sub t$ minimises the standard Bloom FPR $(1 - e sup { - h/c } ) sup h$, capped at 6. Table 3 lists
the default layout: 12.1% of $k$-mers are selected, 4.8% into the first tier and 7.3% into the second.
.TH "Proposition 1" "the layout does not depend on the reference"
The tier boundaries $S sub t$ depend only on $( beta , f sub t , eps sub t sup * )$, and the selection hash only on
the canonical $k$-mer. Hence the set of probes of a contig, their positions and their tiers are the same for every
reference sketched with the same $k$, budget and tier list. Only the hit indicators depend on the reference.
.TQ
.LP
The proof is immediate from (3). Two consequences are used later: a contig is prepared once for a whole database
(Section 9), and the sketch size is exactly linear in genome size with no genome-dependent tuning.
.H2 "The false-positive rate of a blocked filter is known exactly"
.LP
Each tier is a \f[I]blocked\f[R] Bloom filter [1, 9]: an element selects one 512-bit block (one cache line) by one
hash and sets $h$ bits inside it by others. A query costs one cache miss instead of $h$.
.TH "Proposition 2" "exact FPR"
Let block $b$ of a built filter have a fraction $rho sub b$ of its bits set, and let there be $N sub b$ blocks.
Under ideal hashing a $k$-mer that is not in the filter tests positive with probability
.DS
	$eps ~=~ 1 over { N sub b } sum from b ~ rho sub b sup h ,$	(4)
.DE
independently for distinct non-members. Over random filters with $c$ bits per element, the block loads are
Poisson with mean $lambda = 512/c$, and
.DS
	$Ex [ eps ] ~=~ sum from {j >= 0} e sup { - lambda } { lambda sup j } over { j! } ~ q sub h (j) ,~~~~~~ q sub h (j) ~ approx ~ left ( 1 - e sup { - jh/512 } right ) sup h ,$	(5)
.DE
where $q sub h (j)$ is the probability that $h$ uniform positions are all set after $jh$ uniform insertions into
512 bits (computed exactly from the occupancy distribution).
.TQ
.LP
\f[I]Proof.\f[R] A non-member picks a uniform block and $h$ uniform positions within it, independent of the
filter's contents; given the filter the probability is (4). Distinct non-members use independent hashes. The
number of elements per block is binomial, which is Poisson in the limit of many blocks. \(sq
.LP
aif computes (4) from the filter after building it, so the FPR used by every estimator is exact for that filter,
not an approximation. Blocking costs accuracy. The block loads vary, so the FPR is that of a filter with Poisson
loads, which (5) puts at 0.0064 for the first tier against 0.0054 for an unblocked filter with the same bits (an
18% penalty), and 0.1012 against 0.1000 for the second. The measured values on the random reference are 0.0063 and
0.1016 (Table 3). The cap of six hash functions costs the first tier a little more: its target 0.005 would need
$h approx 8$.
.ne 1.6i
.TC 3
Default tier layout ($beta = 0.88$ bits per base). FPRs: \(lqstandard\(rq is an unblocked Bloom filter with the
same bits, \(lqblocked\(rq is (5), \(lqmeasured\(rq is (4) on the random 5 Mbp reference at $k = 20$.
.TCE
.ps 9
.vs 11p
.so gen/tab_layout.tbl
.ps \n(PS
.vs \n(VSp
.H2 "Spurious matches"
.LP
A query $k$-mer that is not shared with the homologous region can still occur elsewhere in $B$: a repeat, a
paralog, or chance. Let $sigma$ be the probability of this for a mutated $k$-mer. A non-shared probe then hits
with the \f[I]effective FPR\f[R]
.DS
	$eps prime sub t ~=~ eps sub t + sigma - eps sub t sigma .$	(6)
.DE
.TH "Proposition 3" "spurious-match rate and its calibration"
(a) If the $N$ $k$-mers of $B$ are independent and uniform, $sigma = 1 - exp ( - 2 N / 4 sup k )$ (both
strands). (b) Let a calibration sample consist of reference $k$-mers with one base substituted at a uniform offset,
kept if selected. If $N sub t$ of them fall in tier $t$ and $Y sub t$ hit, then
$( Y sub t / N sub t - eps sub t )/(1 - eps sub t )$ is unbiased for $sigma sub 1$, the chance-presence rate of
single-substitution neighbours, and the inverse-variance combination over tiers has standard error
.DS
	$roman SE ( sigma hat ) ~ approx ~ left [ sum from t { N sub t (1 - eps sub t ) } over { eps sub t } right ] sup { -1/2 } ,$	(7)
.DE
about $3.5 times 10 sup -4$ for a 5 Mbp genome (10\u6\d draws, 12% selected).
.TQ
.LP
\f[I]Proof.\f[R] (a) A given $k$-mer equals one of $N$ independent uniform $k$-mers or their reverse complements
with probability $1 - (1 - 4 sup -k ) sup {2N}$. (b) A neighbour hits if it is present ($sigma sub 1$) or, if not,
by a false positive ($eps sub t$), so $Ex [ Y sub t / N sub t ] = eps sub t + sigma sub 1 (1 - eps sub t )$;
the binomial variance of $Y sub t$ gives (7) for small $sigma$. \(sq
.LP
The neighbour rate $sigma sub 1$ is the relevant one at low divergence, where a mutated $k$-mer usually differs
from its homolog by one substitution. It is conservative at higher divergence: neighbours of repeat-family
$k$-mers are themselves enriched in the family, while $k$-mers with several substitutions are not. aif clips the
estimate at zero.
.LP
Real genomes are far from random (Table 4). At $k = 17$ their calibrated rates are 5\(en30 times the random
value: 0.0027 for \f[I]E. coli\f[R] and 0.021 for \f[I]P. aeruginosa\f[R], whose 66% GC content concentrates
its $k$-mer spectrum. Pairwise sketches therefore choose $k$ per reference: the smallest $k$ between 17 and 21 with
$sigma hat <= 0.003$ (17 for \f[I]E. coli\f[R], 20 for \f[I]P. aeruginosa\f[R]). The GTDB database uses $k = 20$
throughout; there the median calibrated rate is $4.7 times 10 sup -4$, the 90th percentile $1.6 times 10 sup -3$,
and 2.3% of the 199,923 references exceed 0.003.
.LP
$sigma$ matters in two places. Left out of (6), it biases containment up by about $(1 - theta ) sigma$, which is
$a (1 - theta ) sigma / (k theta )$ in ANI: 0.10 pp at 90% and 0.32 pp at 85% for $sigma = 0.003$. And it raises
the lowest identity the HMM can represent (Section 5.2), because the level grid stops where $theta$ falls to
eight times the smallest effective FPR.
.ne 1.7i
.TC 4
Spurious-match rate $sigma$: the random-sequence value of Proposition 3(a) and aif's calibrated estimate.
Values below the calibration's standard error (about $3.5 times 10 sup -4$) are noise; 0.0 is an estimate clipped
at zero.
.TCE
.ps 9
.vs 11p
.so gen/tab_sigma.tbl
.ps \n(PS
.vs \n(VSp
