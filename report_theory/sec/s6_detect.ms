.H1 "Detection and false reports"
.H2 "A bound on false reports that needs no tuning"
.LP
A contig is reported against a reference when $LR >= e sup lambda$ with $lambda = 10$. The threshold can be
justified without simulation.
.TH "Theorem 4" "false reports"
Suppose a contig has no homology with the reference, and its probes hit independently with probabilities
$q sub t <= eps prime sub t$. Then, for any switching rates, level grid, initial distribution and binning of the
HMM,
.DS
	$Pr ( LR >= e sup lambda ) ~<=~ e sup { - lambda } .$	(18)
.DE
.TQ
.LP
\f[I]Proof.\f[R] $LR = P sub H (x) / P sub U (x)$, where $P sub H$ is the HMM's marginal likelihood, a mixture over
state paths of products of per-probe emissions, and $P sub U$ is that of "all unaligned". For one path and one
probe in state $A sub j$,
$Ex sub q [ P sub j (x) / P sub U (x) ] = 1 + (q - eps prime )( pi sub j - eps prime ) / ( eps prime (1 - eps prime ))$,
which is at most 1 because $pi sub j >= eps prime >= q$; probes in $U$ contribute exactly 1. Independence across
probes and linearity over paths give $Ex sub q [ LR ] <= 1$, and Markov's inequality gives (18). \(sq
.LP
Three remarks. (i) The bound is the reason $sigma$ is calibrated conservatively (Section 3.3): it holds whenever
the true null hit rate does not exceed the assumed one. A reference whose repeats raise the null rate above
$eps prime$ voids it. (ii) aif skips the HMM when every 32-probe window favours some level by 10 nats (a fast path
for confidently aligned contigs). Its score is a maximum over levels rather than a likelihood ratio; by a union
bound over levels, an unrelated contig takes it with probability at most $m e sup -10 = 1.3 times 10 sup -3$ if
it has one window, and with vanishing probability if it has several. (iii) The bound is per reference. When
candidates are chosen from $R$ references using the same probes (database search, Section 9), a guarantee for the
whole search needs a union bound, $R e sup { - lambda }$, which is vacuous at $lambda = 10$ for GTDB's 199,923
references: a family-wise error of $10 sup -4$ would need $lambda >= ln R + 9.2 = 21.4$. (iv) The likelihood of the first
$i$ probes is the marginal of the same model, so $LR$ over prefixes is a non-negative supermartingale under the null,
and Ville's inequality [13] extends (18) to the maximum over prefixes.
.LP
In practice the bound is far from tight (Table 10). Of 65,000 unrelated contigs of 1\(en20 kb, none was reported,
and the largest log-likelihood ratio was 7.7. Among 1 kb contigs, $LR >= e sup 3$ occurred at a rate of
$4.7 times 10 sup -4$ against the bound's 0.050: the empirical tail runs about 100 times below
$e sup { - lambda }$, because most unrelated contigs have $LR << 1$. At that rate the 16 candidates of a database search
would produce a false report for about one contig in $10 sup 5$, but selection of the candidates by their chance
excess of hits makes this an estimate, not a guarantee.
.ne 1.8i
.TC 10
Unrelated contigs (independent uniform sequence) against the random reference: share with $LR >= e sup lambda$,
and the bound (18). \(lq0 (< x)\(rq: none observed among $1/x$ contigs.
.TCE
.ps 9
.vs 11p
.so gen/tab_null.tbl
.ps \n(PS
.vs \n(VSp
.H2 "The shortest contig that can be detected"
.LP
For a contig aligned throughout at identity $a$, the log-likelihood ratio grows linearly with length. Let
.DS
	$r(a) ~=~ max from j ~ sum from t s sub t ~ Ex sub { pi sub t (a) } left [ log { P sub j (x) } over { P sub U (x) } right ] ~~~ roman "nats per base,"$	(19)
.DE
the expected evidence per base for the best-fitting level, and $v(a)$ the corresponding variance.
.TH "Proposition 9" "detection"
To leading order, $ln LR approx L r(a) - ln (2m)$ with standard deviation $sqrt { L v(a) }$, so a contig of length
$L$ is reported with probability about $PHI left ( ( L r(a) - ln (2m) - lambda ) / sqrt { L v(a) } right )$.
If the true identity lies below the grid, $r(a)$ can be negative, and the contig is not detected at any length.
.TQ
.LP
\f[I]Proof sketch.\f[R] The path "aligned at level $j$ throughout" has prior probability $1/(2m)$ times stay
probabilities that are $1 - O( rho L )$, and its log-likelihood ratio is a sum of independent per-probe terms with
mean $L r(a)$ when $j$ is the best level. The marginal likelihood adds the neighbouring levels and paths, which
changes $ln LR$ by $O(1)$. \(sq
.LP
Table 11 and Figure 2 give the prediction and aif's outcome on 200 contigs per cell. Half of the contigs are
detected from 172 bp at 95% ANI, 854 bp at 90%, 1.7 kb at 88%, 3.5 kb at 86% and 9.4 kb at 84%; aif's detection
rates follow the predicted curves at every identity tested, and its median $ln LR$ matches $L r(a) - ln (2m)$
closely (78.6 observed against 78.1 predicted at 5 kb and 90%; 36.4 against 36.3 at 10 kb and 86%).
.LP
Below about 82% ANI, $r(a)$ turns negative at $k = 20$: the lowest level (86.45%) expects more hits than such a
contig produces, and every level fits worse than $U$. This is a hard floor of the default model, set by the grid
rule $theta sub m = 8 eps prime sub min$. Extending the grid lowers the floor but adds levels that look even more
like $U$, which raises $kappa$'s sum in (17), so detection at low identity and AF bias trade against each other
through the same parameter. The floor also moves with $k$: at $k = 17$ the lowest level is 84.45%.
.ne 2.0i
.TC 11
Expected evidence per base (19) and the predicted length at which 50% and 95% of contigs aligned at the given
identity are reported ($LR >= e sup 10$; $k = 20$, default tiers). \(if: never.
.TCE
.ps 9
.vs 11p
.so gen/tab_detect.tbl
.ps \n(PS
.vs \n(VSp
.ne 2.8i
.PSPIC -C figs/fig_detect.eps 6.2i
.FC 2
Share of contigs reported against contig length: Proposition 9 (lines) and aif (points, 200 contigs each),
contigs aligned throughout at the given identity.
.FE
