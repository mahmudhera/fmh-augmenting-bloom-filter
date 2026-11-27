.HA A "Simulation details"
.LP
\f[B]Reference.\f[R] 5,000,000 independent uniform bases (seed 20261127), sketched with \f[CR]aif sketch -k 20\f[R]
and default settings ($beta = 0.88$, tiers \f[CR]0.005:0.6,0.1:0.4\f[R], at most 6 hash functions). Its calibrated
$sigma$ is $2.3 times 10 sup -4$ and its tier FPRs are 0.00635 and 0.10156; the formulas use these values.
.LP
\f[B]Contigs.\f[R] Substrings at uniform positions, each base substituted independently with probability $p$ by
one of the other three bases; unrelated sequence is independent uniform bases from another stream. Queries ran with
aif v1.7 defaults, or with \f[CR]--af-local 0\f[R] (v1.6 AF), \f[CR]--het 0\f[R] (pooled ANI) or
\f[CR]--bin 1\f[R] (no binning), as stated.
.IP \(bu 0.2i
Null: 40,000 contigs of 1 kb, 20,000 of 5 kb and 5,000 of 20 kb.
.IP \(bu 0.2i
Detection: lengths 0.5\(en20 kb \(mu identities 84\(en95%, 200 contigs per cell.
.IP \(bu 0.2i
Boundaries: 5 kb aligned + 5 kb unrelated, 400 contigs per identity.
.IP \(bu 0.2i
Absorption: 3 kb + gap + 3 kb and 3 kb + flank, aligned parts at 97%, 8 lengths, 200 contigs each.
.IP \(bu 0.2i
Mosaics: 10 kb, two halves, 200 contigs per case. Indels: 10 kb, $p = 3%$, single-base indels, 200 per rate.
.LP
\f[B]Idealised Monte Carlo.\f[R] Theorem 1: 20,000 contigs per row, substitutions on the sequence, $k$-mer survival,
independent selection and false positives, aif's two-step combination. One change point: 4,000 boundaries per
identity, 400 probes on each side, exact posterior under a uniform prior. Exact HMM: 300 contigs per identity
(5 kb + 5 kb), probes simulated with overlapping $k$-mers, per-probe scaled forward\(enbackward with aif's states,
emissions and transitions (no bins). Screening: (24) with GTDB R232's per-reference screening FPRs
(\f[CR]aif db info --refs\f[R]) in 40 quantile classes.
.HA B "Reproduction"
.LP
.nf
.ft CR
.ps 8.5
.vs 10.5p
/usr/bin/python3 -I scripts/theory_sims.py prepare   # random reference and sketches
/usr/bin/python3 -I scripts/theory_numbers.py        # closed forms and Monte Carlo
/usr/bin/python3 -I scripts/theory_sims.py           # aif on simulated contigs
/usr/bin/python3 -I scripts/theory_hmm.py            # exact HMM at boundaries
/usr/bin/python3 -I scripts/theory_report.py         # figures and tables
report_theory/build.sh                               # -> aif_theory.pdf
.ps \n(PS
.vs \n(VSp
.ft R
.fi
.LP
Results go to \f[CR]results_theory/\f[R] (\f[CR]numbers.json\f[R], \f[CR]sims.json\f[R],
\f[CR]hmm_boundary.json\f[R]). \f[CR]theory_numbers.py\f[R] reads the GTDB screening FPRs from \f[CR]data/theory/gtdb_dbinfo.tsv\f[R]
(\f[CR]aif db info data/gtdb/aif/gtdb_0.88.aifd --refs\f[R]); the other steps need only the aif binary.
