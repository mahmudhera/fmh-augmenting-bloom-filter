.H1 "Space, linearity and the bits-per-sample principle"
.H2 "Cost accounting"
.LP
All designs are compared by their size in bits per reference base, $beta$.
A tier that samples at rate $s sub t$ and spends $c sub t$ bits per element costs $s sub t c sub t$.
Anchors with positions cost $s sub a ( b sub a + c sub pos )$. Here $b sub a$ is the anchor hash or
fingerprint width (64, or about 32 when chaining filters stray fingerprint matches), and
$c sub pos approx 2 + log sub 2 (1/ s sub a )$ is the Elias\(enFano cost of a position delta.
Table 6 evaluates the designs for a 5 Mbp genome.
.TC 6
Space per reference base and for a 5 Mbp genome. D0 is the baseline (0.064 bits/base).
\(lqD5-lite\(rq is the proposal tuned to the baseline budget.
.TCE
.so gen/tab_space.tbl
.H2 "Dense filters are never space-efficient"
.TH "Theorem 2" "dense AMQs versus the sequence"
Let $B$ have $N = L sub B - k + 1$ distinct $k$-mers. Any AMQ representing $K sub B$ with FPR $eps$
uses at least $N log sub 2 (1/ eps )$ bits [12]. For large $L sub B$ this exceeds the $2 L sub B$ bits of the
2-bit-packed sequence whenever $eps < 1/4$.
A Bloom filter, with any number of hash functions, exceeds 2 bits per element whenever
$eps < 1 - e sup {-1/2} approx 0.39$.
The same holds for an island filter that holds all island $k$-mers (D4), which stores $(1 - s sub a ) N$ elements.
.TQ
.LP
\f[I]Proof.\f[R] The bound is that of Carter et al. [12] for a universe ($4 sup k$) much larger than $N$;
$N log sub 2 (1/ eps ) > 2 L sub B$ reduces to $log sub 2 (1/ eps ) > 2$ as $N / L sub B -> 1$.
For the BF, $c sub BF ( eps ) = - h / ln (1 - eps sup {1/h} )$ equals 2 at $eps = 0.393$ for $h=1$
and at larger $eps$ for no $h >= 2$. \(sq
.LP
So the answer to Q1 for the dense proposal is \f[B]no\f[R].
At 1% FPR, D4 costs 9.6 bits per base (6 MB for a 5 Mbp genome).
That is 4.8\(mu the genome itself and 150\(mu the baseline sketch, and the genome would also permit
exact alignment, i.e. exact ANI and AF.
A dense filter can still be useful as a query-time accelerator, but not as a sketch.
The island tier must be subsampled. With a second FMH threshold, D5 costs 0.10 bits per base
($s sub a = 1/1000$, $s sub 2 = 1/50$, $eps = 0.25$), and \(lqD5-lite\(rq matches the baseline exactly
(Table 6).
.H2 "Linearity"
.LP
Every operation is a single streaming pass (Table 7).
A rolling hash gives each $k$-mer's hash in $O(1)$.
The island index of a $k$-mer is just the number of anchors seen so far, so island keys cost one counter.
A tier-2 insertion or probe touches $h$ filter bits; with blocked Bloom filters [15] that is one cache line,
and at the high-FPR optimum of Theorem 3 (ANI $>=$ 95%) $h = 1$.
Space is linear in $L sub B$ with the explicit constant $beta$. This is the same scaling as FMH and
unlike fixed-size MinHash [2, 3].
The filter must be sized before insertion ($m = c( eps ) s sub 2 L sub B$), which is easy because
$L sub B$ is known; partitioned filters handle unknown sizes.
So the answer to Q2 is \f[B]yes\f[R]: linear time and linear space.
.TC 7
Time and space complexity ($h$ = filter hash functions, $m sub a$ = matched anchors).
.TCE
.TS
center box tab(@);
l l l.
\f[B]operation\f[R]@\f[B]D0 (FMH)\f[R]@\f[B]added by D2 / D5\f[R]
_
build reference@$THETA ( L sub B )$ hashing, $THETA (s L sub B )$ output@$THETA (s sub 2 h L sub B )$ filter writes; island counter
space (bits)@$64 s L sub B$@$s sub 2 c( eps ) L sub B$ + anchors $s sub a (b sub a + c sub pos ) L sub B$
query a contig@$THETA ( L sub A )$ hashing, $THETA (s L sub A )$ lookups@$THETA (s sub 2 h L sub A )$ probes; chaining $O( m sub a log m sub a )$
merge two references@hash-set union@bitwise OR (same $m$, $h$, hashes)
compare two sketches@$THETA (s(L sub A + L sub B ))$ merge@not supported by a BF (Section 4.7); use D9
.TE
.H2 "The bits-per-sample principle"
.LP
Theorem 1 says that, below the knee, accuracy is bought with \f[I]probes\f[R]. A probe costs
$c( eps )$ bits of reference storage, and probes with smaller $c$ are noisier. The next result fixes the
best trade-off.
.TH "Theorem 3" "optimal false-positive rate at fixed budget"
Fix $beta$ bits per reference base for the comparison tier, and assume the sampling-dominated
regime ($s << 1/k$, so $1 - s approx 1$).
An exact 64-bit FMH tier samples at $s = beta / 64$; an AMQ tier at FPR $eps$ samples at
$s sub 2 = beta / c( eps )$. The sampling-plus-false-positive variance in (7) equals
$theta (1- theta ) J / (n beta )$ with
.DS
	$J( eps ) ~=~ c( eps ) left ( 1 ~+~ { eps } over { (1 - eps ) theta } right )$	(11)
.DE
($J = 64$ for the exact tier). The AMQ tier therefore reduces the non-process variance by $R( eps ) = 64 / J( eps )$.
For a Bloom filter with real-valued $h$ ($c = log sub 2 (1/ eps ) / ln 2$) the minimiser satisfies
.DS
	$eps ln (1/ eps ) ~=~ (1 - eps ) [ (1 - eps ) theta + eps ] ,~~~~~~ "and " eps sup * ~ \(ap ~ theta / ln (1/ theta ) " as " theta -> 0 .$	(12)
.DE
.TQ
.LP
\f[I]Proof.\f[R] Substitute $s sub 2 = beta / c$ into the last two terms of (7):
$(1- theta )[ theta + eps /(1- eps )] c / (n beta ) = theta (1 - theta ) J / (n beta )$.
For (12), set $d/d eps$ of $ln (1/ eps ) ( 1 + eps / ((1- eps ) theta ))$ to zero.
For small $theta$ the right side of (12) tends to $theta$, and $eps ln (1/ eps ) = theta$ gives
$eps sup * approx theta / ln (1 / theta )$. \(sq
.TC 8
Optimal FPR and variance-reduction factor $R = 64/J$ at equal bits ($k=21$), computed with the exact
integer-$h$ Bloom cost. $R$ is the factor by which the sampling-plus-FP variance falls relative to 64-bit FMH
in the same space; columns 7\(en9 use other AMQs (Table 2).
.TCE
.so gen/tab_epsopt.tbl
.LP
Table 8 and Figure 3 give the main quantitative answer to Q3.
.IP \(bu 0.2i
At ANI $>=$ 95% the optimal FPR is \f[I]high\f[R], $eps sup * approx$ 0.5\(en0.6.
The filter then costs only 1\(en1.5 bits per element ($h = 1$), and the variance falls 11\(en23\(mu.
Even a conventional $eps = 0.01$ gives 6.5\(mu.
.IP \(bu 0.2i
At 90% ANI, $eps sup * approx 0.05$ and $R sup * = 6.9$. At 80%, $eps sup * approx 0.0017$ and
$R sup * = 4.1$, because a false positive now competes directly with the small signal $theta approx 0.009$.
.IP \(bu 0.2i
Many noisy probes beat a few exact ones. This is the FracMinHash counterpart of $b$-bit minwise
hashing [9], where keeping 1\(en2 bits of each MinHash value is optimal at high similarity.
At high ANI non-shared $k$-mers are rare, so false positives rarely matter and the number of probes dominates.
.IP \(bu 0.2i
The Bloom filter itself is not essential. Binary fuse filters do better ($R sup *$ up to 28), and the
information-theoretic bound caps the gain at $R sup * <= 37$ (99.9%), 16.4 (95%) and 5.9 (80%).
$b$-bit FMH with Elias\(enFano coding is competitive only where small $eps$ is needed, but unlike a BF it is
symmetric (Section 4.7).
.LP
Caveats. (i) The query must enumerate its tier-2 $k$-mers.
(ii) The filter's hash functions must be independent of the selection hash; because $h(x) < s sub 2$ is
constrained, its bits must not be reused for filter positions.
(iii) The optimum is broad (Figure 3), so $eps$ need not be tuned finely.
(iv) The principle applies only while sampling dominates ($s sub 2 <= 1/k$); beyond that the process
variance is the floor.
.ne 3.2i
.sp 0.3v
.PSPIC -C figs/fig_eps.eps 5.0i
.FC 3
Sampling-plus-FP variance of a Bloom-filter tier relative to 64-bit FMH at the same bit budget, i.e.
$J( eps )/64 = 1/R( eps )$, against the FPR. Dots mark the optimum $eps sup *$ for each ANI ($k=21$).
Curves below the dotted line are improvements.
.FE
.H2 "Unknown divergence: hedged tiers"
.LP
$eps sup *$ depends on $theta$, which is unknown until the comparison is made, and a database serves
comparisons at every ANI.
Consider tiers $t = 1, ..., T$ on disjoint hash ranges $[ s sub {t-1} , s sub t )$, so no $k$-mer is stored twice.
Tier $t$ receives a fraction $w sub t$ of the budget at FPR $eps sub t$. Information adds across tiers, so
$1 / J sub design = sum from t w sub t / J( eps sub t )$, and the relative efficiency at $theta$ is
$J sup * ( theta ) / J sub design ( theta )$. Per-tier estimates are combined with plug-in inverse-variance weights.
Table 9 gives minimax designs over 80\(en99.9% ANI (D8).
.TC 9
Robust (minimax) FPR choices: efficiency relative to the $theta$-specific optimum of Table 8, at the same budget.
.TCE
.so gen/tab_hedge.tbl
.LP
A single tier at $eps = 0.044$ guarantees 40% of the attainable gain everywhere.
A two-tier split, with 55% of the bits at $eps = 0.0033$ for distant matches and 45% at $eps = 0.58$
for close ones, guarantees 57%.
Both are better choices than either \(lqtextbook\(rq FPR.
.H2 "Short contigs at an equal budget"
.LP
Table 10 and Figure 4 turn Theorem 1 into ANI standard deviations for contigs of 1\(en100 kb at 95% ANI.
Appendix C (Tables 21\(en22) repeats this at 99% and 90%.
All D2 variants use exactly the D0 budget of 0.064 bits per base.
.TC 10
Standard deviation of the ANI estimate (pp) at 95% ANI ($k=21$), from (7)\(en(8) with the exact process variance.
\(dg: fewer than 30 expected probes, where (8) understates the true spread (Table 5).
.TCE
.so gen/tab_std_p05.tbl
.LP
At 5 kb the standard deviation falls from 2.84 pp (D0) to 1.17 pp ($eps = 0.01$), 1.02 pp (robust
$eps = 0.044$) and 0.91 pp ($eps sup *$). At 1 kb it falls from 6.4 pp to 2.1 pp.
The gain is about 3\(mu in standard deviation, slightly below $sqrt { R sup * }$ because the process variance
is not reduced.
Getting within 1.5\(mu of the floor's variance needs much more:
$beta approx$ 0.28 bits per base at 99% ANI, 0.59 at 95% and 1.17 at 90%.
That is 4\(en18\(mu the D0 budget and a sizeable fraction of the 2-bit genome.
\f[B]In the sketch regime ($beta <= 0.1$) the realistic gain is about 3\(mu in standard deviation at the same
space, or the same precision with about $R sup *$ times fewer bits.\f[R]
.ne 3.2i
.sp 0.3v
.PSPIC -C figs/fig_std.eps 6.4i
.FC 4
ANI standard deviation versus contig length at 95% (left) and 99% (right) ANI, all sketches at the
D0 budget (0.064 bits per base). The containment floor and alignment oracle are intrinsic limits.
.FE
.H2 "The asymmetry principle"
.LP
The gain above depends on the query enumerating its $k$-mers. When both sides are sketches, or when
the roles are swapped, the picture changes.
.TH "Proposition 2" "orientation of the filter"
Let the contig contribute $N = s sub 2 n$ tier-2 $k$-mers and the genome $n sub B = s sub 2 L sub B$, with $L sub B >> L sub A$.
On the count scale, the false-positive variance is:
(a) $N (1- theta ) eps / (1- eps )$ when contig $k$-mers probe a genome filter;
(b) about $s sub 2 L sub B eps (1 - eps )$ when genome $k$-mers probe a contig filter, which scales with the
genome rather than the contig;
(c) about $0.443 ( n sub B + N ) / log sub 2 (1/ eps )$ for BF\(caBF comparison, which estimates $| A \(ca B |$
by inclusion\(enexclusion from fill ratios [44, 19, 20].
.TQ
.LP
Part (c) follows from the linear-counting variance $(m/h sup 2 )(e sup a - 1 - a)$ at $a = hn/m = ln 2$
(Appendix B). Table 11 compares signal-to-noise ratios.
.TC 11
Signal-to-noise ratio for the shared-$k$-mer count at 95% ANI ($s sub 2 = 1/50$, $eps = 0.01$, 5 Mbp genome).
\(lqLocalized\(rq probes only the genome region matched by anchors (about $2 L sub A$).
.TCE
.so gen/tab_sym.tbl
.LP
For a 5 kb contig the correct orientation has SNR 7.2. Symmetric BF\(caBF has 0.42, which is useless;
it needs contigs of roughly 50\(en100 kb to become informative (SNR 4 at 50 kb, 8 at 100 kb). The reversed orientation has 1.07.
Three rules follow.
(1) Put the AMQ on the large side and enumerate the small side.
(2) If the database holds contigs and the query is a genome, \f[I]localise first\f[R]: use anchors to find
the matching genome region, then probe only there (SNR 7.0; Section 5.5).
(3) For symmetric all-versus-all use, store $b$-bit FMH fingerprints (D9). Both sides stay enumerable at
$log sub 2 (1/ eps ) + 2$ bits per element, and Theorem 1 applies unchanged.
