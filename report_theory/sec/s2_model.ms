.H1 "Model and notation"
.H2 "Contigs, segments and the estimands"
.LP
A contig $A$ of length $L$ has $n = L - k + 1$ $k$-mer starts. It is a concatenation of segments. An
\f[I]aligned\f[R] segment $i$ of length $l sub i$ is homologous to part of the reference $B$ and differs from
it by independent substitutions at rate $p sub i = 1 - a sub i$ (the simple substitution model [7]). An
\f[I]unaligned\f[R] segment has no homolog in $B$. The two estimands are
.DS
	$roman AF = { sum sub { i ~ roman aligned } l sub i } over L ,~~~~~~ roman ANI = { sum sub { i ~ roman aligned } l sub i a sub i } over { sum sub { i ~ roman aligned } l sub i } ,$	(1)
.DE
which are the quantities an alignment reports (AF as the aligned span, ANI as matches over aligned columns;
Section 10 treats indels). For one aligned segment at identity $a$ a $k$-mer survives unmutated with probability
.DS
	$theta = a sup k ,$	(2)
.DE
the \f[I]containment\f[R] of the aligned part. At $k = 20$, $theta$ is 0.82 at 99% ANI, 0.36 at 95%, 0.12 at
90% and 0.039 at 85%: the signal fades geometrically with divergence, which sets most of the limits below.
.H2 "Hashing and consistency"
.LP
Results are stated under \f[I]ideal hashing\f[R]: the selection hash of a canonical $k$-mer is uniform on
$[0,1)$ and independent across distinct $k$-mers, and the filter hashes are independent of it and of each
other. aif uses 64-bit mixing functions with separate seeds for selection, block and bit positions, so the
assumption is the usual one for sketching. Selection is by FracMinHash: $x$ is selected when $h(x) < S$, and
its tier is the hash sub-range that contains $h(x)$. The property everything rests on is
\f[I]consistency\f[R]: a $k$-mer shared by $A$ and $B$ is selected, and put in the same tier, in both or in
neither. Testing a selected query $k$-mer against the reference's filter of selected $k$-mers is therefore as good
as testing it against all of $B$, up to false positives.
.ne 3.2i
.TC 2
Notation (default values at $k = 20$).
.TCE
.TS
center box tab(@);
l lw(4.6i).
\f[B]symbol\f[R]@\f[B]meaning\f[R]
_
$k$, $L$, $n$@$k$-mer length (20); contig length; $k$-mer starts, $n = L - k + 1$
$a$, $p$, $theta$@identity, substitution rate $1 - a$, containment $a sup k$
$beta$, $f sub t$@sketch budget in bits per reference base (0.88); budget share of tier $t$ (0.6, 0.4)
$s sub t$, $s$@selection rate of tier $t$ (0.048, 0.073); total $s = sum s sub t$ (0.121 probes per base)
$eps sub t$, $sigma$, $eps prime sub t$@filter FPR; spurious-match rate; effective FPR $eps + sigma - eps sigma$ (0.0066, 0.102)
$c( eps )$@bits per element of a Bloom filter at FPR $eps$
$N sub t$, $Y sub t$@probes and hits in tier $t$
$pi sub t (a)$@hit probability in tier $t$ when aligned at identity $a$: $theta + (1 - theta ) eps prime sub t$
$a sub j$, $m$@aligned levels of the HMM, 0.9995 down to 0.8645 in steps of 0.005; $m = 28$
$rho sub u$, $rho sub h$@per-base switching rates aligned\(<>unaligned ($3 times 10 sup -6$) and between levels ($10 sup -4$)
$LR$@likelihood ratio of the HMM against "all unaligned"; reporting threshold $e sup 10$
$D sub A$, $D sub U$@per-probe KL divergences $KL ( P sub A || P sub U )$ and $KL ( P sub U || P sub A )$
.TE
