# FracMinHash + Bloom filters for ANI / AF on short contigs

**Status:** theoretical analysis only. No sketching method has been implemented yet.

**Deliverable:** [`fmh_bloom_filter_theory.pdf`](fmh_bloom_filter_theory.pdf) (27 pages)

## Question

Does augmenting a FracMinHash (FMH) sketch with Bloom filters improve containment-based ANI estimation
(ANI = C^(1/k)) on short contigs, and can it estimate alignment fraction (AF)?
The specific proposal: keep FMH k-mers with their positions ("anchors") and store the k-mers of the
"islands" between anchors in a Bloom filter. The island k-mers can be subsampled with a second FMH threshold.

## Headline results

| Question | Answer |
|---|---|
| Space-efficient? | Not with **all** island k-mers. Any AMQ with FPR < 1/4 costs more than the 2-bit genome (Theorem 2). With a second FMH threshold it can match the plain-FMH budget. For short contigs the dense anchors dominate the cost. |
| Linear? | Yes. Single streaming pass; linear time; linear space with a tunable constant. |
| Better ANI? | Yes, in the asymmetric setting (the query enumerates its k-mers). At equal bits, sampling variance falls by 4x (80% ANI) to 23x (99.9%), with an unexpectedly high optimal FPR (Theorem 3). For a 5 kb contig at 95% ANI the s.d. goes from 2.84 to 0.91 pp. |
| AF? | Containment cannot separate AF from ANI (Prop. 3). Positions fix this, and query-side positions suffice (Prop. 5). Tier-2 hits locate boundaries to about 1/(s2*theta) bp. |
| Close ideas | Island-keyed localisation (allows small k at low ANI; fixes contigs-in-database), island fingerprints (best per bit at >= 99.9% ANI), hedged multi-FPR tiers, b-bit FMH for symmetric use. |

The report gives a unified variance framework (Theorem 1, checked by Monte Carlo), 14 designs, 16 scenarios,
a concrete estimator (AIF-ANI) and an experimental plan.

## Layout

```
fmh_bloom_filter_theory.pdf   the report
report/                       groff sources (main.ms, sec/*.ms, tmac/fmhdoc.tmac), generated tables (gen/), figures (figs/)
scripts/compute_tables.py     evaluates the closed-form results -> tables + EPS figures
scripts/mc_check_theorem1.py  small Monte Carlo of the idealised model, used to check Theorem 1 (not an implementation)
Makefile                      `make pdf` (groff + ghostscript only) or `make` (also regenerates tables; needs numpy + matplotlib)
```
