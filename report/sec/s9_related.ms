.H1 "Related work and close ideas"
.LP
\f[B]Sketch-based similarity.\f[R]
MinHash [2] and Mash [3] estimate Jaccard similarity from fixed-size sketches. FracMinHash [4, 5] samples
consistently at a fixed fraction, so containment is estimable, and Hera et al. [6] derived its mutation-rate
confidence intervals using the $k$-mer statistics of Blanca et al. [7]. Theorem 1 adds the false-positive
term to that variance and makes the dependence on the query-side sampling rate explicit.
Sylph [27] estimates containment ANI for low-coverage metagenomes with a statistical model; CMash [36]
estimates containment at several $k$ from one structure.
.LP
\f[B]Bloom filters in containment estimation.\f[R]
The nearest precedent for the Bloom tier is the containment MinHash of Koslicki and Zabeti [8]: put the larger
set in a Bloom filter, sketch the smaller one, and count sketch elements found in the filter.
D2 differs in three ways. It subsamples the filtered set consistently (FMH), which keeps the filter
sketch-sized (equation (4)). It corrects for the known FPR exactly. And it chooses the FPR by the
bits-per-sample criterion of Theorem 3, which favours far higher FPRs than are usual.
Set-cardinality and intersection estimation from filter fill ratios goes back to linear counting [44] and
appears in chemoinformatics [19] and distributed systems [20]; Proposition 2(c) explains why it fails for small
query sets.
.LP
\f[B]Bits per sample.\f[R]
$b$-bit minwise hashing [9] showed that storing a few bits of each MinHash value beats storing many full values
in the same space. BinDash [10] and Dashing / Dashing 2 [11, 45] apply related ideas to genomes.
Theorem 3 is the FracMinHash analogue, for asymmetric membership probing.
.LP
\f[B]AMQ structures.\f[R]
Bloom filters [1, 13], compressed Bloom filters [14], blocked filters and Golomb-coded sets [15], cuckoo [18],
quotient [41], xor [16], binary fuse [17] and ribbon [40] filters trade space, speed and FPR differently;
the bound of Carter et al. [12] underlies Theorem 2.
Sequence-search indexes built from Bloom filters (SBT [21], BIGSI [22], COBS [42], Raptor/IBF [23], HIBF [43])
answer containment queries over many datasets. They target presence/absence thresholds rather than calibrated
ANI, and their keys are global $k$-mers, not island-keyed.
.LP
\f[B]Positional sketches and ANI/AF.\f[R]
MashMap [24] maps with minimiser sketches and estimates per-segment identity. FastANI [25] maps 3 kb query
fragments and reports ANI as mean fragment identity and AF as the fraction mapped. skani [26] chains sparse
sampled seeds, computes ANI on chained regions and reports AF, and was designed for fragmented MAGs.
The AIF shares skani's \(lqchain first, estimate on aligned regions\(rq structure. What it adds is a calibrated
probabilistic layer, an FP-corrected filter tier with exact variance and bits-per-sample optimisation, plus
localisation-enabled small $k$, island fingerprints for strain level, and hedged FPR tiers.
.LP
\f[B]Anchors and chunking.\f[R]
Minimisers [28], winnowing [29], syncmers [30], strobemers [31] and minimiser-space de Bruijn graphs [32]
provide anchors with spacing or robustness guarantees; minimiser Jaccard estimators are biased [33], which
motivates D13's split between localising anchors and an FMH estimation tier.
Island fingerprints are content-defined chunking [34] combined with Merkle-style hashing [35].
Spaced seeds [39] are a further way to make tests robust to dense substitutions at low identity.
