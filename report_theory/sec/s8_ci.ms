.H1 "Uncertainty"
.H2 "Statistical intervals"
.LP
\f[B]ANI.\f[R] The statistical standard error is (13) with the conditional variance (14), evaluated at the
estimate, over the aligned $k$-mers. For heterogeneous contigs it is replaced by the spread of the window
identities of the test (21), divided by the effective number of windows, when that is larger.
.LP
\f[B]AF.\f[R] The aligned length is uncertain only where the per-probe aligned posterior $gamma$ is between 0 and 1,
that is, around boundaries. aif treats each run of intermediate posteriors, bracketed by confident states, as one
boundary whose position $X$ has survival function $Pr (X > x) = gamma (x)$ (aligned on the left) or
$1 - gamma (x)$; a run with the same state on both sides is split at its extremum into two boundaries. Then
.DS
	$Var ( X ) ~=~ 2 int (x - x sub 0 ) ~ Pr (X > x) ~ dx ~-~ left ( int Pr (X > x) ~ dx right ) sup 2 ,~~~~~~ Var ( roman AF hat ) ~=~ 1 over { L sup 2 } sum from { roman boundaries } Var ( X ) ,$	(22)
.DE
the usual moments of a non-negative variable from its survival function, summed over boundaries treated as
independent. In the idealised model the posterior is calibrated: its standard deviation is within about 10% of the
actual error (Table 8, independent probes). On real contigs both the boundary model and the substitution model are
wrong in ways that make posteriors too narrow: runs of misses, the low levels of Section 5.4, indels and
rearrangements. The statistical intervals are therefore only the input to calibration.
.H2 "Calibrated intervals"
.LP
The default intervals are $theta hat \(+- q sub g sigma hat$, where $sigma hat$ is the statistical standard
error with a small floor (0.05 pp for ANI, $1/L$ + 1 pp for AF) and $q sub g$ is a multiplier for the contig's
stratum $g$: four estimated-identity classes \(mu three length classes \(mu fully or partly aligned, 24 strata.
The multipliers are split-conformal quantiles fitted on contigs with alignment truth.
.TH "Proposition 11" "coverage per stratum"
Let the calibration set contain $n sub g$ contigs in stratum $g$ with scores
$R sub i = | roman error sub i | / sigma hat sub i$, and let $q sub g$ be the
$left ceiling (n sub g + 1)(1 - alpha ) right ceiling$-th smallest. If a new contig in stratum $g$ is exchangeable
with them, its interval covers the truth with probability at least $1 - alpha$, and at most
$1 - alpha + 1/(n sub g + 1)$ when scores are almost surely distinct.
.TQ
.LP
\f[I]Proof.\f[R] Under exchangeability the rank of the new score among the $n sub g + 1$ is uniform, and coverage
is the event that it is at most $left ceiling (n sub g + 1)(1 - alpha ) right ceiling$ (split-conformal
prediction [14, 15]; stratifying on features of the contig itself keeps exchangeability within strata). \(sq
.LP
The guarantee holds as far as exchangeability does, and it is the assumption that fails in practice. Intervals
fitted on pairwise queries of real genomes (held-out coverage 94.4% for ANI and 94.2% for AF at 95%) covered only
92% and 87% (aif v1.6) when used for database search of gut contigs against GTDB, a different population of
errors. v1.7
therefore fits a separate table for database search, on gut contigs disjoint from the evaluation set, which covers
95.1% and 94.2% there. The theory says what to watch: coverage is guaranteed only for contigs resembling the
calibration set, stratum by stratum.
