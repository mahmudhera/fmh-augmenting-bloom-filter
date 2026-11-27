#!/usr/bin/env python3
"""AF bias of aif variants by truth-ANI class on T (top-5 hits, self reference excluded), incl. pairs below 0.90"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uhgg_af as A  # noqa: E402
import uhgg_compare as U  # noqa: E402

stage = os.environ.get("STAGE", "T")
conv = os.environ.get("CONV", "q")
S = A.members(stage)
X = U.contexts()
TR = U.load_truth(stage, conv)
BINS = [(0.0, 0.80, "<0.80"), (0.80, 0.85, "0.80-0.85"), (0.85, 0.88, "0.85-0.88"), (0.88, 0.90, "0.88-0.90"),
        (0.90, 0.95, "0.90-0.95"), (0.95, 1.01, ">=0.95")]
for name in sys.argv[1:]:
    H = U.load_aif(os.path.join(U.C, "runs", f"{stage}_{name}", "aif.tsv"))
    out = {}
    for c in S:
        tr = {r: (a, f) for r, a, f in TR.get(c, [])}
        sref = X[S[c]["genome"]][1]
        for h in [h for h in H.get(c, []) if h[0] != sref][:5]:
            if h[0] in tr:
                b = next(lab for lo, hi, lab in BINS if lo <= tr[h[0]][0] < hi)
                out.setdefault(b, []).append(h[2] - tr[h[0]][1])
    print(f"{name:12s} " + "  ".join(f"{lab}: n {len(out[lab])} bias {100 * sum(out[lab]) / len(out[lab]):+.1f} "
                                     f"mae {100 * sum(map(abs, out[lab])) / len(out[lab]):.1f}"
                                     for _, _, lab in BINS if lab in out))
