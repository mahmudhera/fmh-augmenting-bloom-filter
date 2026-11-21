.H1 "Scenario analysis"
.LP
Sixteen scenarios cover the situations in which the proposal would be used.
Each states the setting, applies the results above and ends with a verdict.
Table 20 at the end summarises them.
.H3 "S1. Species-level short contig (95% ANI, 5 kb; raw contig queried against sketched references)."
Plain FMH gives about 5 probes, a 0.7% chance of an empty sketch, a standard deviation of 2.84 pp and a bias of
about $-0.9$ pp. At the same 0.064 bits per base, a Bloom tier gives 0.91\(en1.17 pp (Table 10).
For ANI alone no reference positions are needed: D2 with the robust FPR of 0.044 (1.02 pp) is enough.
Island keys (D5) need anchors dense enough to land at least 3 matches in the contig. Proposition 4 requires
$s sub a approx 1/200$ with $k sub a = 15$, which costs 0.13\(en0.21 bits per base (16- to 32-bit fingerprints).
\f[B]For short contigs the anchors, not the filter, dominate the cost of the proposal.\f[R]
These anchors are not wasted: they are themselves an FMH tier at rate $s sub a$, and they are 5\(mu denser than D0.
\f[I]Verdict:\f[R] D2 (+ D12 for AF) at the D0 budget gives about 3\(mu lower spread. D5 is worth its extra
anchor cost when localisation is needed (S3, S6, S10, S11).
.H3 "S2. Strain-level short contig (99.9% ANI, 5 kb)."
The oracle standard deviation is 0.045 pp; the contig holds only about 5 substitutions, so the relative error is
at least 45% for any method. Using the asymptotic efficiencies of Table 18, D0 ($eta = 0.020$) gives about
0.32 pp, D2 with $eps sup * = 0.63$ ($eta = 0.33$) about 0.078 pp, and island fingerprints at the D0 budget
($eta = 0.46$) about 0.066 pp. At 99.9%, 21-mer tests keep only 3% of the attainable information per probe
(Table 16). Anchors survive almost always ($theta sub a approx 0.98$), and 800 bp islands are identical with
probability 0.45, close to the most informative regime.
\f[I]Verdict:\f[R] D6, or D2 at high FPR. This is the regime where the island idea shines, in its fingerprint form.
.H3 "S3. Distant homology (80\(en85% ANI)."
At 80%, $theta approx 0.009$. The optimal FPR is small (0.0017) and the Bloom tier gains only 4.1\(mu.
The containment floor itself is poor ($eta = 0.12$ at $k=21$), and 21-mers keep 32% of the information per probe.
Island-local $k sub e = 10$ recovers 2.9\(mu per probe and raises the floor to 0.46, but anchoring needs about
36 kb at 80% even with $k sub a = 15$, $s sub a = 1/200$ (Table 12).
\f[I]Verdict:\f[R] for contigs of at least 20\(en40 kb, D5 + D7 (island-local small $k$) is the most informative
sketch. For short contigs at 80% every nucleotide sketch is weak, and protein-level (AAI) or alignment methods are
appropriate.
.H3 "S4. Partial overlap (AF < 1)."
Containment reports $alpha sup {1/k} (1-p)$: $-3.1$ pp at $alpha = 0.5$ and $-9.9$ pp at $alpha = 0.1$ (Table 13).
This is often larger than the variance problem. Query-side segmentation (D12) or chaining plus tier-2 extension
(D5) estimates AF with about 150 bp boundary error (95% ANI, $s sub 2 = 1/50$), and restricting $Chat$ to the
aligned span removes the bias, leaving about 0.1 pp.
\f[I]Verdict:\f[R] the most important practical fix for contigs. It needs positions on the query side only,
plus reference anchors if colinearity matters.
.H3 "S5. Very short contigs (1\(en2 kb)."
At 1 kb and 95% ANI the oracle standard deviation is already 0.69 pp. D0 gives 6.4 pp (and fewer than one
expected probe), and the best Bloom tier 2.1 pp. Anchoring a 1 kb contig needs $s sub a approx 1/74$, about
0.33 bits per base even with 16-bit fingerprints.
\f[I]Verdict:\f[R] sketches give about \(+-2 pp at best. For sub-pp precision, retrieve candidate reference
regions (with D0/D2 as a filter) and align, or run D11 on the retrieved sequence.
.H3 "S6. Reversed asymmetry (contigs in the database, genome as the query)."
Probing genome $k$-mers against a contig's filter generates about $s sub 2 L sub B eps$ false positives:
1,000 at $eps = 0.01$, against a signal of 34 for a 5 kb contig. The SNR is 1.07 (Table 11).
Matching the genome's anchors to the contig's anchors first, and probing only the matched region, restores SNR 7.0.
\f[I]Verdict:\f[R] island keys and anchors are \f[I]necessary\f[R] here. This is the strongest argument for the
positional part of the proposal.
.H3 "S7. Symmetric all-versus-all comparison of sketches."
Neither side can enumerate $k$-mers. BF\(caBF inclusion\(enexclusion has SNR 0.42 at 5 kb and needs roughly 50\(en100 kb
to become usable (Table 11).
\f[I]Verdict:\f[R] use $b$-bit FMH (D9): enumerable fingerprints at $log sub 2 (1/ eps ) + 2$ bits each, analysed by
Theorem 1, with $R$ from the \(lq$b$-bit FMH\(rq column of Table 8 (4.8\(en11.9).
.H3 "S8. Indels."
Each indel event destroys $k - 1$ or $k + d - 1$ query $k$-mers and is counted as roughly one substitution, so
$k$-mer ANI runs low by about the indel rate. Span discrepancies $DELTA$ between matched anchors expose net indel
lengths (Section 5.6). Island keys are unaffected by indels inside islands; fixed-offset keys are not.
\f[I]Verdict:\f[R] positions allow an indel-corrected ANI that unordered sketches cannot provide.
.H3 "S9. Heterogeneous divergence (recombination, hypervariable regions)."
Containment overestimates ANI by the power-mean bias: +0.6 to +2.2 pp at 95% (Table 19), more for larger $k$.
Window-wise estimation along the query removes it, with windows of at least 30\(en50 probes.
\f[I]Verdict:\f[R] relevant mainly to genome-scale ANI; it needs query-side positions only.
.H3 "S10. Repeats and paralogs."
Global probes match any copy of a repeated $k$-mer, inflating containment and AF in repeat-adjacent contigs,
which assemblies produce in excess. Island-keyed probes match only the chained locus (Proposition 6).
\f[I]Verdict:\f[R] D5 > D2. The size of the effect on real assemblies is an open empirical question (E4).
.H3 "S11. Chimeric or misassembled contigs."
A contig whose halves come from distant loci, or from different genomes, looks fully aligned to set-based methods
and to query-side segmentation (D12). Chaining yields two discordant chains.
\f[I]Verdict:\f[R] only reference-side positions (D3/D5) detect this. Report per-chain AF and ANI.
.H3 "S12. Large eukaryotic references (Gbp)."
Global spurious matches at $k=21$ reach $6.8 times 10 sup -4$ per probe for 3 Gbp under the random model,
and repeats make it much worse. That biases $Chat$ upward unless subtracted, and it rules out small $k$.
The filter scales linearly: 37.5 MB per 3 Gbp at 0.10 bits per base.
\f[I]Verdict:\f[R] island keys are essential, and anchors should be processed per chromosome or contig.
.H3 "S13. Database scale (about 10\u5\d references)."
At the D0 budget the space is the same as D0 (about 4 GB for $10 sup 5$ genomes of 5 Mbp). Probing every
reference's filter costs $O( roman { #refs } times s sub 2 L sub A h )$, so screen first with the anchor or D0 tier
through an inverted index and refine the top candidates with tier-2 probes.
Bit-sliced multi-reference filters (BIGSI [22], COBS [42], IBF/HIBF [23, 43]) answer one probe for every
reference at once, but their keys must be reference-independent, which rules out island keys.
\f[I]Verdict:\f[R] a two-stage design (screen, then refine) keeps queries sublinear in practice.
.H3 "S14. Unknown divergence."
The optimal FPR ranges over three orders of magnitude across 80\(en99.9% ANI.
\f[I]Verdict:\f[R] hedged tiers (D8). Two tiers at FPR 0.0033 and 0.58 guarantee 57% of the attainable gain
everywhere, and one tier at 0.044 guarantees 40% (Table 9).
.H3 "S15. Sequencing or assembly errors."
An error rate $e$ adds to the divergence ($p sub eff approx p + e$). Errors concentrate at contig ends, where
query-side hit density can flag and trim them. D6 is the most sensitive design, because one error spoils an
island. At 99.9% ANI, $e = 10 sup -4$ inflates $p$ by 10%.
\f[I]Verdict:\f[R] report error-inclusive ANI, or trim low-density ends, before strain-level calls.
.H3 "S16. Streaming construction, multi-contig references and merging."
Everything is single-pass and linear (Table 7). Filters with identical parameters merge by bitwise OR. For
multi-contig references (draft genomes) island indices must be globally unique, e.g. (contig id, $j$).
Positions are per sequence, so merging positioned references needs a concatenated coordinate space.
\f[I]Verdict:\f[R] no obstacle.
.ne 3.4i
.TC 20
Scenario summary. \(lqGain\(rq is relative to D0 at equal space unless stated otherwise.
.TCE
.TS
center box tab(@);
l lw(1.25i) lw(1.55i) lw(2.05i).
\f[B]#\f[R]@\f[B]situation\f[R]@\f[B]recommended\f[R]@\f[B]expected effect / binding constraint\f[R]
_
S1@95%, 5 kb@D2 (+D12); D5 if localising@SD 2.84 \(-> 0.9\(en1.2 pp; anchors cost 0.13\(en0.21 b/base
S2@99.9%, 5 kb@D6 or D2 ($eps approx 0.6$)@$eta$: 0.02 \(-> 0.33\(en0.46; few substitutions per contig
S3@80\(en85%@D2 (small $eps$) + D7@4\(mu (D2) + 2.9\(mu (local $k$); anchors need 15\(en40 kb
S4@AF < 1@D12 / D5 + span-restricted $Chat$@removes $-3$ pp ($alpha = 0.5$) bias; 150 bp boundaries
S5@1\(en2 kb@D2, then retrieve and align@SD about 2 pp floor for sketches
S6@contigs in database@D5 with localisation@SNR 1.07 \(-> 7.0
S7@sketch vs sketch@D9 ($b$-bit FMH)@BF\(caBF unusable below about 50 kb
S8@indels@D3/D5 span discrepancies@indel-corrected ANI
S9@rate heterogeneity@windowed estimation@removes +0.6 to +2.2 pp bias
S10@repeats@D5 island keys@suppresses paralog hits
S11@chimeras@D3/D5 chaining@detects discordant chains
S12@Gbp references@D5 (keys essential)@$sigma$: $7 times 10 sup -4$ \(-> about 0
S13@$10 sup 5$ references@screen (D0/D3) + refine (D2/D5)@same space as D0; sublinear queries
S14@unknown ANI@D8 hedged tiers@$>=$ 57% of optimum everywhere
S15@errors@trim ends; error-aware D6@$e$ inflates $p$ at strain level
S16@streaming / merging@any@linear; OR-merge
.TE
