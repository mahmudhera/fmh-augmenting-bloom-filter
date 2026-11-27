.H1 "Database search"
.H2 "Exactness"
.TH "Proposition 12" "database search equals pairwise query"
For a database built with $k$, budget and tier list, the estimates of a contig against any evaluated reference
equal those of a pairwise query against that reference's sketch built with the same $k$, budget and tier list.
.TQ
.LP
\f[I]Proof.\f[R] By Proposition 1 the probes are the same for every reference. Each reference's filters are built
by the same procedure as a pairwise sketch (same seeds; block count from its own size), its FPRs are computed by
(4) and its $sigma$ by the same calibration. The estimators depend only on the probes, their hits and these
rates. \(sq
.LP
An exhaustive search is byte-identical to pairwise queries (real-genome report, Section 10). The database
therefore changes only \f[I]which\f[R] references are evaluated, which is decided by a screen.
.H2 "The screening index"
.LP
The screen uses $k$-mers whose selection hash lies below $s sub scr = 1/64$, a subset of the first tier's probes.
Each reference is a column of a bit-sliced Bloom filter (as in COBS [10]) with about 8 bits per key and 3 hash
functions; references of similar size share rows, so one 64-byte row answers one hash function for 512
references. A reference's screening FPR $f sub d$ is its column's fill cubed (Proposition 2). Over GTDB R232 the
median is 0.028, and 98% of references lie between 0.024 and 0.030.
.LP
A contig contributes $n sub s = min (128, ~ approx L/64)$ keys, those of smallest hash. Reference $d$ shares
$c sub d$ of them and passes when $c sub d >= 2$ and
.DS
	$z sub d ~=~ { c sub d - n sub s f sub d } over { sqrt { n sub s f sub d (1 - f sub d ) } } ~>=~ 2 .$	(23)
.DE
Passing references are ranked by the FP-corrected screening containment $(c sub d / n sub s - f sub d )/(1 - f sub d )$,
and the top 16 go to the model.
.TH "Proposition 13" "ranking by z favours low-FPR references"
For a reference with containment $theta$, $Ex [ z sub d ] approx sqrt { n sub s } ~ theta sqrt { (1 - f sub d ) / f sub d }$.
Ranking by $z$ therefore prefers, at equal containment, the reference with the lower screening FPR, and the
preference grows with $sqrt { n sub s }$; halving $f$ is worth a factor 1.4 in containment. The FP-corrected
containment is unbiased for $theta$ whatever $f sub d$.
.TQ
.LP
\f[I]Proof.\f[R] $c sub d ~ \(ap ~ roman Bin ( n sub s , theta + (1 - theta ) f sub d )$, so
$Ex [ c sub d ] - n sub s f sub d = n sub s theta (1 - f sub d )$. \(sq
.LP
This is why aif uses $z$ only as a significance test and ranks by containment (an earlier version ranked by $z$,
which promoted references with low screening FPRs on long contigs).
.H2 "Recall at database scale"
.LP
The screen sees only $n sub s approx L/64$ keys, so short contigs have few: 15 at 1 kb. With 199,923 references,
many unrelated ones reach two or three shared keys by chance.
.TH "Proposition 14" "screening recall"
Let $R$ unrelated references have screening FPRs $f sub d$. The expected number that pass the screen is
$sum sub d Pr ( roman Bin ( n sub s , f sub d ) >= tau sub d )$, with $tau sub d$ the smallest count meeting (23)
and $c >= 2$. The true reference, with containment $theta$, is among the top $K$ by screening containment with
probability
.DS
	$sum from c Pr ( c sub true = c ) ~ Pr left ( roman Poisson left ( LAMBDA sub > (c) + half LAMBDA sub = (c) right ) <= K - 1 right ) ,$	(24)
.DE
where $LAMBDA sub > (c)$ and $LAMBDA sub = (c)$ are the expected numbers of passing unrelated references whose
containment exceeds or ties that of $c$ shared keys (ties broken at random).
.TQ
.LP
\f[I]Proof.\f[R] Unrelated references' counts are independent binomials; the number exceeding a threshold is a sum
of independent rare Bernoulli variables, hence Poisson in the limit. \(sq
.LP
Table 14 evaluates (24) with GTDB's screening FPRs. A 1 kb contig has 13,154 unrelated references passing the
screen. At 95% ANI the true reference expects 5.6 of its 15 keys and reaches the top 16 in 72% of cases; at 90%
ANI it expects 2.2 keys and makes it in 6%. Recall improves quickly with length (2 kb: 98% at 95%; 5 kb: 72% at
90%; 10 kb: 95% at 90%) and is negligible at 85% for any length shown. \f[B]For short divergent contigs the screen,
not the model, is the limit of database search.\f[R] Re-ranking a larger screened set by the FP-corrected
containment of all probes (\f[CR]--prescreen\f[R] $K$; about eight times as many probes as keys) recovers most of
it: with $K = 1024$, 96% at 1 kb and 95% ANI, and 38% at 90%. Two caveats apply. Related references (other strains
of the species, or the genus) also compete, but they are useful candidates rather than noise. And keys are not
independent of the probes the model uses, which is why Theorem 4 needs its union bound in database mode.
.ne 2.2i
.TC 14
Screening recall at GTDB scale (199,923 references, screening FPRs as measured; Proposition 14): keys per contig,
expected unrelated references passing, and for the true reference the probability of passing the screen, of
being among the 16 evaluated, and of being among the top 1,024 (as with \f[CR]--prescreen 1024\f[R]).
.TCE
.ps 9
.vs 11p
.so gen/tab_screen.tbl
.ps \n(PS
.vs \n(VSp
.H2 "The pre-model filter"
.LP
Candidates without a clear excess of hits are first scored by a two-state Viterbi path at eight fixed identities,
with the HMM's emission ratios and switching costs, and skipped when the score is below 5. This is a heuristic,
not a bound: the Viterbi score maximises over paths at levels that are not the HMM's grid and carries no level
prior, while $ln LR$ sums over paths. Neither bounds the other. Its threshold was set where the database benchmark
showed no change in output; a provable filter would need an upper bound on $ln LR$, for example the Viterbi score
over the HMM's own grid plus the log number of paths.
