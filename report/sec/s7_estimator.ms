.H1 "A concrete proposal: the AIF-ANI estimator"
.LP
This section assembles the pieces that survive the analysis into one estimator built on the
\f[I]anchored island filter\f[R] (AIF, design D5 with optional D6/D7/D8 components).
It is stated precisely so that the implementation phase can test it.
.H2 "Reference sketch"
.IP (a)
\f[B]Anchors.\f[R] Canonical $k sub a$-mers with $h sub a (x) < s sub a$ (or syncmers; Section 5.7), stored as $b sub a$-bit
fingerprints with Elias\(enFano-coded positions, in a hash map from fingerprint to (position, island index).
.IP (b)
\f[B]Island filter.\f[R] One AMQ (Bloom or binary fuse) over keys $(x, j)$ for canonical $k sub e$-mers with
$h sub e (x) < s sub 2$. $h sub e$ is independent of $h sub a$ and of the filter's own hashes, and $j$ is the island index.
Optionally there are two hash-range tiers at different FPRs (D8, Table 9).
.IP (c)
\f[B]Island fingerprints\f[R] (optional, D6). An $f$-bit hash of each island interior.
.IP (d)
\f[B]Metadata.\f[R] $L sub B$, parameters, and each tier's fill ratio, from which its exact FPR follows.
.LP
Cost: $beta = s sub a ( b sub a + c sub pos + f ) + sum from t s sub {2,t} c( eps sub t )$ bits per base,
for example 0.064 (D5-lite) or 0.10 (D5) (Table 6).
.H2 "Query algorithm (raw contig)"
.CS
1. ANCHOR  compute A's anchors; look up fingerprints in the anchor map;
           chain colinear matches per orientation (O(m log m), gap-tolerant LIS);
           keep chains with >= m0 anchors.  None? -> fall back to D12 (global tier)
           or report "insufficient anchors" with the Prop. 4 feasibility bound.
2. PROBE   for each tier-2 query k-mer y at position t inside a chain span:
             islands J(t) <- exact if the flanking matched anchors are consecutive
                              in B, else interpolated (1-2 islands);
             Z_t <- [ (y, j) in filter for some j in J(t) ];
             eps_t <- 1 - (1 - eps)^|J(t)|
3. EXTEND  probe tier-2 k-mers beyond each chain end against the flanking island;
           segment hits with a 2-state HMM -> aligned span(s) S of A
4. AF      AF_hat = |S| / L_A
5. C       C_aln = ( sum_{t in S} Z_t - sum eps_t ) / sum_{t in S} (1 - eps_t)
6. ANI     ANI_hat = C_aln^(1/k_e), bias-corrected with (8);
           CI by inverting Theorem 1 (exact inversion if N_aln*theta < 30)
7. REFINE  (optional) indel correction from span discrepancies (Sec. 5.6);
           window averaging for heterogeneity (Sec. 6.4); island-fingerprint
           estimate (Sec. 6.2) combined by inverse variance at high ANI
.CE
.LP
Step 5 generalises (5) to probe-specific FPRs. Since $Ex [ Z sub t ] = U sub t (1 - eps sub t ) + eps sub t$, the estimator
is a $(1 - eps sub t )$-weighted mean of the shared indicators $U sub t$, and so is unbiased for the containment of the span.
Every step is linear in $L sub A$ apart from chaining, which is $O( m log m )$ in the number of matched anchors.
.H2 "Starting parameters"
.TS
center box tab(@);
l l l.
\f[B]parameter\f[R]@\f[B]suggested start\f[R]@\f[B]reason\f[R]
_
$k sub a$, $s sub a$@15, 1/200 to 1/500 (contigs 2\(en10 kb)@Prop. 4, Table 12
$k sub e$@21 (species level); 10\(en12 with island keys below 90%@Prop. 7, Tables 16\(en17
$s sub 2$@1/20 to 1/70, depending on budget@Theorem 3, Table 10
$eps$@two tiers (0.003, 0.58) or one tier at 0.044@Table 9
$f$ (island fingerprint)@8 bits if strain level matters@Prop. 8, Table 18
budget@0.064 (= D0) to 0.10 bits per base@Table 6
.TE
.H2 "Predictions to test"
.IP P1. 0.4i
At equal space, the ANI standard deviation on 5 kb contigs at 95% ANI falls about 3\(mu (2.84 to about 1.0 pp) relative to D0.
.IP P2. 0.4i
AF boundaries are located to within about 150 bp at $s sub 2 = 1/50$ (95% ANI), and the $alpha sup {1/k}$ bias disappears.
.IP P3. 0.4i
At 80% ANI, island-local $k sub e = 10$\(en12 gives about 2.9\(mu more information per probe than $k = 21$, provided anchors exist.
.IP P4. 0.4i
At ANI $>=$ 99.9%, island fingerprints reach $eta approx 0.46$ at the D0 budget (D0 has 0.02).
.IP P5. 0.4i
BF-versus-BF comparison fails below about 50 kb, whereas orientation-correct probing does not.
.IP P6. 0.4i
Monte Carlo variances match Theorem 1 once real genomes (repeats, low complexity, indels) are used, up to the
spurious-match term $sigma$, which island keys should keep near zero.
