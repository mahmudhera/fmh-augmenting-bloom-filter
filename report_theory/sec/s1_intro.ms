.H1 "Introduction"
.LP
Two numbers describe how a contig relates to a reference genome: the identity of the part that aligns (ANI) and
the fraction of the contig that aligns (AF). Alignment gives both but needs the reference sequence. Sketching tools
(Mash [2], sourmash and FracMinHash [3, 4], FastANI [8], skani [6]) compare compressed representations instead,
but lose accuracy on short contigs, and AF in particular. aif gives both from a sketch of the reference that costs 0.88 bits per base, about 1/18 of the 2-bit sequence:
.IP \(bu 0.2i
\f[B]Sketch.\f[R] The canonical $k$-mers of the reference whose selection hash falls below a threshold (about
12% of them) are inserted into two cache-blocked Bloom filters, one per hash range, with FPRs of about 0.6% and
10%.
.IP \(bu 0.2i
\f[B]Query.\f[R] Every selected $k$-mer of the raw contig is tested in its tier. The hits and misses, in contig
order, are read by an HMM with one unaligned state and 28 aligned states at identities from 86.5% to 99.95%
($k = 20$).
.IP \(bu 0.2i
\f[B]Estimates.\f[R] AF is the posterior aligned length. Containment is the FP-corrected hit rate over the
aligned part, and ANI is its $k$-th root with a bias correction (or, for mosaic contigs, the posterior mean
identity). A contig is reported when the HMM's likelihood ratio against "all unaligned" is at least $e sup 10$.
.IP \(bu 0.2i
\f[B]Database search.\f[R] Candidates come from a bit-sliced screening index. Up to 16 references per contig
reach the model.
.LP
The project's first theory document (\f[CR]fmh_bloom_filter_theory.pdf\f[R], November 21) came before any code.
It analysed 14 candidate designs, including reference-side anchors and island filters, and found two principles
that the implementation kept: storing many samples at a high FPR beats storing fewer exactly (bits per sample),
and AF needs positions on the query only. aif uses neither anchors nor islands. It is the two-tier global filter
of that analysis, hedged across FPRs, with query-side segmentation by an HMM. Its defaults were then tuned on
data (the benchmark and real-genome reports).
.LP
This report is the theory of the method as built. It states the estimators, what they are unbiased for, how
much they vary, where AF and detection are limited, what the database screen loses, and how the model departs
from real sequence. Each result is checked in two ways (Table 1):
.IP (a) 0.3i
\f[B]Idealised model.\f[R] Monte Carlo of the simple substitution model with independent selection and false
positives (\f[CR]scripts/theory_numbers.py\f[R], \f[CR]scripts/theory_hmm.py\f[R]).
.IP (b) 0.3i
\f[B]The aif binary.\f[R] Contigs cut from a random 5 Mbp reference (independent uniform bases), mutated as
stated, and queried against its default sketch at $k = 20$ (\f[CR]scripts/theory_sims.py\f[R]). The random
reference removes repeats, so differences from theory come from the method, not the genome.
.LP
All numbers use aif v1.7 defaults at $k = 20$, the database $k$. Pairwise queries choose $k$ between 17 and
21 per reference (Section 3.3); the formulas hold for any $k$.
.ne 3i
.TC 1
Results of this report and how each was checked. MC: Monte Carlo of the idealised model.
.TCE
.ps 9
.vs 11p
.TS
center box tab(@);
l lw(2.75i) lw(1.95i).
\f[B]result\f[R]@\f[B]statement\f[R]@\f[B]check\f[R]
_
Proposition 1@tier layout and probes independent of the reference@code; database = pairwise, byte-identical
Proposition 2@exact FPR of a blocked Bloom filter from its fill@predicted 0.0064 / 0.1012, measured 0.0063 / 0.1016
Proposition 3@spurious-match rate and its calibration@random and four real genomes, $k$ = 17\(en21 (Table 4)
Theorem 1@FP-corrected containment: unbiased, closed-form variance@MC: SD within 3% (Table 5)
Proposition 4@inverse-variance combination of tiers@MC: bias $<=$ 3% relative at 2 kb, 85% ANI
Theorem 2@bits per sample; regret of the two-tier hedge@closed form (Table 6)
Proposition 5@variance given the realised identity@used by the statistical intervals
Proposition 6@counts cannot separate AF from ANI@closed form
Proposition 7@change-point resolution of a boundary@exact posterior; aif (Table 8)
Theorem 3@absorption of unaligned stretches by low levels@aif v1.6 and v1.7 (Figure 1)
Proposition 8@what the local-identity AF of v1.7 can remove@aif (Figure 1)
Theorem 4@false reports: $Pr (LR >= e sup lambda ) <= e sup { - lambda }$@aif: 0 of 65,000 unrelated contigs (Table 10)
Proposition 9@expected LR and the shortest detectable contig@aif: 36 length\(muANI cells (Figure 2)
Proposition 10@power-mean bias of pooled ANI on mosaics@aif \f[CR]--het 0\f[R] (Table 13)
Proposition 11@split-conformal coverage per stratum@measured coverage (earlier reports)
Proposition 12@database search equals pairwise query@byte-identical outputs (real-genome report)
Proposition 13@z-ranking favours low-FPR references@closed form
Proposition 14@screening recall at database scale@closed form with GTDB's screening FPRs (Table 14)
Proposition 15@indels act as mismatches@aif (Table 15)
.TE
.ps \n(PS
.vs \n(VSp
