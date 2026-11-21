#!/usr/bin/env python3
"""Evaluate the closed-form results of the report and emit tbl tables + EPS figures.

Nothing here implements a sketch. Every number is a direct evaluation of a
formula stated in the report (Theorems/Propositions referenced in comments).

Usage: python3 compute_tables.py <out_gen_dir> <out_fig_dir> [mc_csv]
"""
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("ps")
import matplotlib.pyplot as plt  # noqa: E402

K = 21                      # estimation k
LB = 5_000_000              # reference genome length (bacterial)
S0 = 1 / 1000               # baseline FracMinHash scale
HBITS = 64                  # bits per stored hash in baseline
BETA0 = HBITS * S0          # baseline bits per reference base (0.064)
LN2 = math.log(2)

# palette (validated, light mode): blue, orange, aqua, yellow, magenta
C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"


# ----------------------------------------------------------------------------
# Core formulas (Section 3)
# ----------------------------------------------------------------------------
def theta(p, k=K):
    return (1 - p) ** k


def Bterm(p, k=K):
    q = 1 - theta(p, k)
    return q * (2 * k - 1 + 2 / p) - 2 * k


def var_nmut_exact(n, p, k=K):
    q = 1 - theta(p, k)
    return (n * (1 - q) * (q * (2 * k + 2 / p - 1) - 2 * k)
            + k * (k - 1) * (1 - q) ** 2
            + 2 * (1 - q) / p ** 2 * ((1 + (k - 1) * (1 - q)) * p - q))


def var_C(n, p, s, eps, k=K):
    """Theorem 1: process + sampling + false-positive variance of C-hat."""
    th = theta(p, k)
    return (var_nmut_exact(n, p, k) / n ** 2
            + th * (1 - th) * (1 - s) / (n * s)
            + (1 - th) * eps / ((1 - eps) * n * s))


def var_ani(n, p, s, eps, k=K):
    th = theta(p, k)
    return var_C(n, p, s, eps, k) * (1 - p) ** 2 / (k ** 2 * th ** 2)


def eta(p, s, eps, k=K, n=10 ** 7):
    """Asymptotic efficiency vs. the alignment oracle p(1-p)/n."""
    return p * (1 - p) / n / var_ani(n, p, s, eps, k)


# ----------------------------------------------------------------------------
# AMQ cost models (Section 4)
# ----------------------------------------------------------------------------
def c_bf(e):          # Bloom filter, best integer number of hash functions h >= 1
    """bits/element m/n = -h / ln(1 - e^(1/h)), minimized over integer h."""
    return min(-h / math.log(1 - e ** (1 / h)) for h in range(1, 40))


def h_bf(e):
    return min(range(1, 40), key=lambda h: -h / math.log(1 - e ** (1 / h)))


def c_bf1(e):         # Bloom filter with one hash function (= bitmap)
    return 1 / (-math.log(1 - e))


def c_ef(e):          # b-bit FracMinHash fingerprints, Elias-Fano coded
    return math.log2(1 / e) + 2


def c_lb(e):          # information-theoretic AMQ lower bound
    return math.log2(1 / e)


def c_fuse(f):        # binary fuse filter (3-wise), f-bit fingerprints, eps = 2^-f
    return 1.125 * f


EPS_GRID = np.logspace(-6, math.log10(0.7), 6000)


def J(e, th, cfun=c_bf):
    return cfun(e) * (1 + e / ((1 - e) * th))


def J_vec(th, cfun=c_bf):
    return np.array([J(e, th, cfun) for e in EPS_GRID])


def eps_opt(th, cfun=c_bf):
    jj = J_vec(th, cfun)
    i = int(np.argmin(jj))
    return EPS_GRID[i], jj[i]


def fuse_opt(th):
    best = None
    for f in range(1, 25):
        e = 2.0 ** -f
        v = c_fuse(f) * (1 + e / ((1 - e) * th))
        if best is None or v < best[1]:
            best = (f, v)
    return best


# ----------------------------------------------------------------------------
# helpers for tbl output
# ----------------------------------------------------------------------------
def tbl(path, fmt_header, fmt_body, header_rows, rows, opts="center box tab(@)"):
    lines = [".TS", opts + ";"]
    lines += fmt_header
    lines.append(fmt_body + ".")
    for i, hr in enumerate(header_rows):
        # tbl: a spanned ('s') column takes no data entry -> drop its placeholder
        fmt = fmt_header[i].split() if i < len(fmt_header) else []
        if fmt and len(fmt) == len(hr):
            hr = [h for h, f in zip(hr, fmt) if f != "s"]
        lines.append("@".join(hr))
    lines.append("_")
    for r in rows:
        if r == "_":
            lines.append("_")
        else:
            lines.append("@".join(r))
    lines.append(".TE")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def f2(x, d=2):
    return f"{x:.{d}f}"


def sci(x, d=1):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, d) >= 10:
        m /= 10
        e += 1
    return f"{m:.{d}f}\\(mu10\\u\\s-2{e}\\s+2\\d"


def kb(L):
    return f"{L/1000:g}"


# ----------------------------------------------------------------------------
def main():
    gen, figs = sys.argv[1], sys.argv[2]
    mc_csv = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(gen, exist_ok=True)
    os.makedirs(figs, exist_ok=True)
    nums = {}   # scalar values quoted in the text

    # ---------------- Table: Monte Carlo check of Theorem 1 ----------------
    if mc_csv and os.path.exists(mc_csv):
        rows = []
        with open(mc_csv) as fh:
            hdr = fh.readline().strip().split(",")
            for line in fh:
                d = dict(zip(hdr, line.strip().split(",")))
                rows.append([d["n"], d["p"], d["s"], d["eps"],
                             f2(float(d["mean_C_sim"]), 4) + " / " + f2(float(d["theta"]), 4),
                             sci(float(d["var_C_sim"]), 2), sci(float(d["var_C_thy"]), 2), f2(float(d["ratio_C"]), 3),
                             sci(float(d["var_ANI_sim"]), 2), sci(float(d["var_ANI_thy"]), 2), f2(float(d["ratio_ANI"]), 3)])
        tbl(f"{gen}/tab_mc.tbl",
            ["c c c c c c s s c s s", "c c c c c c c c c c c"],
            "n n n n c n n n n n n",
            [["", "", "", "", "mean $C hat$", "Var $C hat$", "", "", "Var ANI", "", ""],
             ["$n$", "$p$", "$s$", "$epsilon$", "sim / $theta$", "sim", "theory", "ratio", "sim", "theory", "ratio"]],
            rows, opts="center box tab(@)")

    # ---------------- Table: efficiency vs s (Theorem 1, Cor. 1) ----------------
    ps = [0.001, 0.005, 0.01, 0.05, 0.10, 0.20]
    ss = [1, 1 / 2, 1 / 10, 1 / 20, 1 / 50, 1 / 200, 1 / 1000, 1 / 5000]
    rows = []
    for s in ss:
        r = [("1" if s == 1 else f"1/{round(1/s)}"), f2(1 + (1 - s) / (K * s), 1)]
        for p in ps:
            r.append(f2(eta(p, s, 0.0), 3))
        rows.append(r)
    tbl(f"{gen}/tab_eff.tbl",
        ["c c c s s s s s", "c c c c c c c c"],
        "c n n n n n n n",
        [["", "small-$p$", "efficiency $eta$ at ANI =", "", "", "", "", ""],
         ["$s$", "VIF", *[f"{100*(1-p):g}%" for p in ps]]],
        rows)
    nums["eta_floor_p05"] = eta(0.05, 1, 0)
    nums["eta_floor_p10"] = eta(0.10, 1, 0)
    nums["eta_s1000_p05"] = eta(0.05, S0, 0)

    # ---------------- Table: AMQ cost models ----------------
    eps_list = [0.5, 0.3, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001]
    rows = []
    for e in eps_list:
        f = -math.log2(e)
        fuse = f2(c_fuse(round(f)), 2) if abs(f - round(f)) < 1e-9 else "\\(em"
        rows.append([f"{e:g}", f2(c_lb(e)), f2(c_bf(e)), f2(c_bf1(e)), fuse, f2(c_ef(e)),
                     f2(HBITS / c_bf(e), 1)])
    tbl(f"{gen}/tab_cost.tbl",
        ["c c c c c c c"],
        "n n n n n n n",
        [["FPR $epsilon$", "lower bound", "BF (best integer $h$)", "BF ($h=1$)", "binary fuse", "$b$-bit FMH (EF)", "64 / $c sub BF$"]],
        rows)

    # ---------------- Table: optimal eps per ANI (Theorem 3) ----------------
    anis = [0.999, 0.995, 0.99, 0.98, 0.97, 0.95, 0.90, 0.85, 0.80]
    rows = []
    for a in anis:
        p = 1 - a
        th = theta(p)
        e_bf, j_bf = eps_opt(th, c_bf)
        e_lb, j_lb = eps_opt(th, c_lb)
        e_ef, j_ef = eps_opt(th, c_ef)
        f_fu, j_fu = fuse_opt(th)
        j01 = J(0.01, th, c_bf)
        rows.append([f"{100*a:g}", f2(th, 4), f"{e_bf:.3g}", f2(c_bf(e_bf), 2), f2(HBITS / j_bf, 1),
                     f2(HBITS / j01, 1), f2(HBITS / j_fu, 1) + f" ($f$={f_fu})", f2(HBITS / j_ef, 1), f2(HBITS / j_lb, 1)])
        if abs(a - 0.95) < 1e-9:
            nums["eps_opt_95"] = e_bf
            nums["R_opt_95"] = HBITS / j_bf
            nums["R_01_95"] = HBITS / j01
        if abs(a - 0.80) < 1e-9:
            nums["eps_opt_80"] = e_bf
            nums["R_opt_80"] = HBITS / j_bf
    tbl(f"{gen}/tab_epsopt.tbl",
        ["c c c c c c c c c", "c c c c c c c c c"],
        "n n n n n n l n n",
        [["ANI", "", "BF", "BF", "BF", "BF", "binary fuse", "$b$-bit FMH", "lower"],
         ["(%)", "$theta$", "$epsilon sup *$", "bits/elt", "$R sup *$", "$R$ at $epsilon$=0.01", "$R sup *$", "$R sup *$", "bound $R sup *$"]],
        rows)

    # ---------------- robust (minimax) single-tier and two-tier designs ----------------
    ani_grid = np.linspace(0.80, 0.999, 60)
    th_grid = theta(1 - ani_grid)
    Jstar = np.array([eps_opt(t, c_bf)[1] for t in th_grid])
    e_cand = np.logspace(-4.5, math.log10(0.7), 160)
    Jmat = np.array([[J(e, t) for t in th_grid] for e in e_cand])     # (E, T)
    eff1 = Jstar[None, :] / Jmat
    worst1 = eff1.min(axis=1)
    i_rob = int(np.argmax(worst1))
    eps_rob = e_cand[i_rob]
    nums["eps_rob"] = eps_rob
    nums["worst_rob"] = worst1[i_rob]
    best2 = (0, None)
    ws = np.linspace(0.0, 1.0, 41)
    for i in range(0, len(e_cand), 2):
        for j in range(i, len(e_cand), 2):
            info = ws[:, None] / Jmat[i][None, :] + (1 - ws[:, None]) / Jmat[j][None, :]
            eff = Jstar[None, :] * info
            wmin = eff.min(axis=1)
            kidx = int(np.argmax(wmin))
            if wmin[kidx] > best2[0]:
                best2 = (wmin[kidx], (e_cand[i], e_cand[j], ws[kidx]))
    worst2, (e1, e2, w2) = best2
    nums.update(dict(hedge_e1=e1, hedge_e2=e2, hedge_w=w2, hedge_worst=worst2))

    def eff_at(design, a):
        th = theta(1 - a)
        js = eps_opt(th, c_bf)[1]
        if design[0] == "single":
            return js / J(design[1], th)
        _, ea, eb, w = design
        return js * (w / J(ea, th) + (1 - w) / J(eb, th))

    designs = [("single fixed $epsilon$ = 0.01", ("single", 0.01)),
               ("single fixed $epsilon$ = 0.25", ("single", 0.25)),
               (f"single minimax $epsilon$ = {eps_rob:.3g}", ("single", eps_rob)),
               (f"two-tier: {w2:.2f} at $epsilon$={e1:.2g} + {1-w2:.2f} at $epsilon$={e2:.2g}", ("two", e1, e2, w2))]
    rows = []
    for name, d in designs:
        effs = [eff_at(d, a) for a in ani_grid]
        rows.append([name, f2(min(effs), 2)] + [f2(eff_at(d, a), 2) for a in [0.80, 0.85, 0.90, 0.95, 0.99, 0.999]])
    tbl(f"{gen}/tab_hedge.tbl",
        ["l c c s s s s s", "l c c c c c c c"],
        "l n n n n n n n",
        [["design (all at the same bit budget)", "worst case", "relative efficiency at ANI =", "", "", "", "", ""],
         ["", "80\\(en99.9%", "80%", "85%", "90%", "95%", "99%", "99.9%"]],
        rows)

    # ---------------- Table + Figure: std(ANI) vs contig length ----------------
    Ls = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
    blocks = []
    figdata = {}
    for p in [0.01, 0.05, 0.10]:
        th = theta(p)
        e_star, _ = eps_opt(th, c_bf)
        s_star = BETA0 / c_bf(e_star)
        s_01 = BETA0 / c_bf(0.01)
        s_rob = BETA0 / c_bf(eps_rob)
        designs = [
            ("oracle", None, None),
            ("floor", 1.0, 0.0),
            ("D0", S0, 0.0),
            ("D2-0.01", s_01, 0.01),
            ("D2-rob", s_rob, eps_rob),
            ("D2-opt", s_star, e_star),
        ]
        rows = []
        curves = {d[0]: [] for d in designs}
        for L in Ls:
            n = L - K + 1
            r = [kb(L)]
            for name, s, e in designs:
                if name == "oracle":
                    sd = math.sqrt(p * (1 - p) / L)
                    mark = ""
                else:
                    sd = math.sqrt(var_ani(n, p, s, e))
                    mark = "\\(dg" if n * s < 30 else ""
                curves[name].append(100 * sd)
                r.append(f2(100 * sd, 2) + mark)
            rows.append(r)
        blocks.append((p, th, e_star, s_star, s_01, s_rob, rows))
        figdata[p] = curves
    for (p, th, e_star, s_star, s_01, s_rob, rows) in blocks:
        tbl(f"{gen}/tab_std_p{int(round(p*100)):02d}.tbl",
            ["c c c c c c c", "c c c c c c c", "c c c c c c c"],
            "n n n n n n n",
            [["contig", "alignment", "floor", "D0: FMH", "D2: BF tier", "D2: BF tier", "D2: BF tier"],
             ["(kb)", "oracle", "($s$=1)", "$s$=1/1000", "$epsilon$=0.01", f"$epsilon$={eps_rob:.2g} (robust)", f"$epsilon sup *$={e_star:.2g}"],
             ["", "", "", "", f"$s sub 2$=1/{round(1/s_01)}", f"$s sub 2$=1/{round(1/s_rob)}", f"$s sub 2$=1/{round(1/s_star)}"]],
            rows)
        if abs(p - 0.05) < 1e-9:
            nums["s_star_95"] = s_star
            nums["s_01"] = s_01
    # example numbers for text (5 kb, 95%)
    n5 = 5000 - K + 1
    nums["sd5k_D0"] = 100 * math.sqrt(var_ani(n5, 0.05, S0, 0))
    nums["sd5k_D2"] = 100 * math.sqrt(var_ani(n5, 0.05, nums["s_star_95"], nums["eps_opt_95"]))
    nums["sd5k_or"] = 100 * math.sqrt(0.05 * 0.95 / 5000)
    nums["sd5k_floor"] = 100 * math.sqrt(var_ani(n5, 0.05, 1, 0))

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharex=True)
    styles = [("oracle", "alignment oracle", INK2, ":", None),
              ("floor", "containment floor (s=1)", INK, "-", None),
              ("D0", "D0: FMH s=1/1000", C1, "-", "o"),
              ("D2-0.01", "D2: BF tier, eps=0.01", C2, "--", "s"),
              ("D2-opt", "D2: BF tier, eps*", C3, "-.", "^")]
    for ax, p in zip(axes, [0.05, 0.01]):
        cur = figdata[p]
        x = np.array(Ls) / 1000
        for key, lab, col, ls, mk in styles:
            ax.plot(x, cur[key], color=col, ls=ls, lw=1.6, marker=mk, ms=4, label=lab)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"ANI = {100*(1-p):g}%  (k = 21, equal bit budget)", fontsize=8.5, color=INK)
        ax.set_xlabel("contig length (kb)", fontsize=8, color=INK2)
        ax.grid(True, which="major", color=GRID, lw=0.5)
        ax.tick_params(labelsize=7, colors=INK2)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.set_xticks([1, 2, 5, 10, 20, 50, 100])
        ax.set_xticklabels(["1", "2", "5", "10", "20", "50", "100"])
    axes[0].set_ylabel("std. dev. of ANI estimate (pp)", fontsize=8, color=INK2)
    # direct labels on left panel at x=100
    for key, lab, col, ls, mk in styles:
        y = figdata[0.05][key][-1]
        axes[0].annotate({"oracle": "oracle", "floor": "floor", "D0": "D0", "D2-0.01": "D2 (0.01)",
                          "D2-opt": "D2 (eps*)"}[key], (100, y), xytext=(3, 0), textcoords="offset points",
                         fontsize=6.5, color=INK, va="center")
    axes[0].set_xlim(0.8, 260)
    axes[1].set_xlim(0.8, 260)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(f"{figs}/fig_std.eps")
    plt.close(fig)

    # ---------------- Figure: J(eps)/64 vs eps ----------------
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    sel = [(0.995, C1, "-"), (0.99, C2, "--"), (0.95, C3, "-."), (0.90, C4, ":"), (0.80, C5, "-")]
    xs = np.logspace(-5, math.log10(0.7), 400)
    end_labels = []
    for a, col, ls in sel:
        th = theta(1 - a)
        ys = [J(e, th) / HBITS for e in xs]
        ax.plot(xs, ys, color=col, ls=ls, lw=1.8, label=f"ANI {100*a:g}%")
        e_b, j_b = eps_opt(th)
        ax.plot([e_b], [j_b / HBITS], marker="o", ms=5, color=col, mec="#fcfcfb", mew=1.2)
        end_labels.append([ys[-1], f"{100*a:g}%"])
    # spread right-end labels so they do not collide (log-space minimum gap)
    end_labels.sort()
    for i in range(1, len(end_labels)):
        if end_labels[i][0] < end_labels[i - 1][0] * 1.22:
            end_labels[i][0] = end_labels[i - 1][0] * 1.22
    for y, lab in end_labels:
        ax.text(0.78, y, lab, fontsize=7, color=INK, va="center", ha="left")
    ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
    ax.text(1e-5, 1.06, "64-bit FMH (no FP)", fontsize=7, color=INK2, ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(5e-6, 2.2)
    ax.set_xlabel("Bloom-filter false-positive rate  eps", fontsize=8, color=INK2)
    ax.set_ylabel("sampling variance relative to\n64-bit FMH at equal bits", fontsize=8, color=INK2)
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.tick_params(labelsize=7, colors=INK2)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=3)
    fig.tight_layout()
    fig.savefig(f"{figs}/fig_eps.eps")
    plt.close(fig)

    # ---------------- Figure: efficiency vs s ----------------
    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    sgrid = np.logspace(-4, 0, 200)
    for p, col, ls in [(0.001, C1, "-"), (0.01, C2, "--"), (0.05, C3, "-."), (0.10, C4, ":")]:
        ys = [eta(p, s, 0.0) for s in sgrid]
        ax.plot(sgrid, ys, color=col, ls=ls, lw=1.8, label=f"ANI {100*(1-p):g}%")
        ax.annotate(f"{100*(1-p):g}%", (1, ys[-1]), xytext=(3, 0), textcoords="offset points",
                    fontsize=7, color=INK, va="center")
    ax.plot(sgrid, [1 / (1 + (1 - s) / (K * s)) for s in sgrid], color=INK2, lw=0.9, ls=(0, (1, 2)),
            label="small-p approx. 1/(1+(1-s)/(ks))")
    ax.axvline(1 / 1000, color=GRID, lw=1)
    ax.text(1 / 1000, 1.02, " s = 1/1000", fontsize=7, color=INK2, va="bottom")
    ax.axvline(1 / K, color=GRID, lw=1)
    ax.text(1 / K, 1.02, " s = 1/k", fontsize=7, color=INK2, va="bottom")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.12)
    ax.set_xlim(1e-4, 2.5)
    ax.set_xlabel("query-side sampling rate s", fontsize=8, color=INK2)
    ax.set_ylabel("efficiency vs. alignment", fontsize=8, color=INK2)
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.tick_params(labelsize=7, colors=INK2)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.legend(fontsize=6.5, frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.62))
    fig.tight_layout()
    fig.savefig(f"{figs}/fig_eff.eps")
    plt.close(fig)

    # ---------------- Table: space accounting ----------------
    def cpos(s):
        return 2 + math.log2(1 / s)
    sa_lite = 1 / 2000
    anchor_bits_lite = sa_lite * (32 + cpos(sa_lite))
    s2_lite = (BETA0 - anchor_bits_lite) / c_bf(0.25)
    nums["s2_lite"] = s2_lite
    space = [
        ("2-bit packed genome (exact; alignment possible)", 2.0),
        ("D0: FMH, $s$=1/1000, 64-bit hashes", BETA0),
        ("D3: D0 + Elias\\(enFano positions", S0 * (64 + cpos(S0))),
        ("D1: global BF of all $k$-mers, $epsilon$=0.01", c_bf(0.01)),
        ("D1: global BF of all $k$-mers, $epsilon$=0.25", c_bf(0.25)),
        ("D4: anchors 1/1000 + island BF of ALL island $k$-mers, $epsilon$=0.01",
         (1 - S0) * c_bf(0.01) + S0 * (32 + cpos(S0))),
        ("D2: BF tier $s sub 2$=1/50, $epsilon$=0.25 (no anchors)", c_bf(0.25) / 50),
        ("D5: anchors 1/1000 (32-bit fp + pos.) + island BF $s sub 2$=1/50, $epsilon$=0.25",
         S0 * (32 + cpos(S0)) + c_bf(0.25) / 50),
        (f"D5-lite: anchors 1/2000 + island BF $s sub 2$=1/{round(1/s2_lite)}, $epsilon$=0.25 (= D0 budget)",
         anchor_bits_lite + s2_lite * c_bf(0.25)),
        ("D6: island fingerprints, anchors 1/800 + 8-bit island hash", (1 / 800) * (32 + cpos(1 / 800) + 8)),
    ]
    rows = []
    for name, b in space:
        mb = b * LB / 8 / 1e6
        rows.append([name, f2(b, 3), (f"{mb:.2f} MB" if mb >= 1 else f"{mb*1000:.0f} KB"), f2(b / BETA0, 1)])
    tbl(f"{gen}/tab_space.tbl",
        ["l c c c"],
        "l n n n",
        [["design", "bits / base", "5 Mbp genome", "\\(mu D0"]],
        rows)

    # ---------------- Table: naive ANI bias under AF < 1 (Prop. 3) ----------------
    rows = []
    for af in [1.0, 0.9, 0.75, 0.5, 0.25, 0.1]:
        r = [f"{af:g}"]
        for k in [21, 31]:
            for a in [0.99, 0.95]:
                est = af ** (1 / k) * a
                r.append(f"{100*est:.2f} ({100*(est-a):+.2f})")
        rows.append(r)
    tbl(f"{gen}/tab_afbias.tbl",
        ["c c s c s", "c c c c c"],
        "n n n n n",
        [["", "$k$ = 21", "", "$k$ = 31", ""],
         ["AF $alpha$", "true 99%", "true 95%", "true 99%", "true 95%"]],
        rows)

    # ---------------- Table: anchor feasibility (Prop. 4) ----------------
    from math import exp
    lo, hi = 0.0, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - exp(-mid) * (1 + mid + mid ** 2 / 2) < 0.95:
            lo = mid
        else:
            hi = mid
    lam3 = hi
    nums["lam3"] = lam3
    rows = []
    for (sa, ka) in [(1 / 1000, 21), (1 / 200, 21), (1 / 200, 15), (1 / 50, 15), (1 / 50, 12)]:
        r = [f"1/{round(1/sa)}", str(ka), f2(sa * (32 + cpos(sa)), 3)]
        for a in [0.99, 0.95, 0.90, 0.85, 0.80]:
            th_a = (a) ** ka
            Lmin = lam3 / (sa * th_a)
            r.append(f"{Lmin/1000:.1f}" if Lmin < 1e6 else ">1000")
        r.append(sci(LB / 4 ** ka, 1))
        rows.append(r)
    tbl(f"{gen}/tab_anchor.tbl",
        ["c c c c s s s s c", "c c c c c c c c c"],
        "n n n n n n n n n",
        [["anchor", "", "anchor", "min. contig length (kb) at ANI =", "", "", "", "", "global spurious"],
         ["rate $s sub a$", "$k sub a$", "bits/base", "99%", "95%", "90%", "85%", "80%", "anchor rate"]],
        rows)

    # ---------------- Table: optimal k (Prop. 6) ----------------
    def g(u):
        return u * u / math.expm1(u)
    u_star = 1.5936242600400401
    gmax = g(u_star)
    nums["u_star"] = u_star

    def kmin(N, sig=1e-3):
        return math.ceil(math.log(N / sig, 4))
    kg, kl = kmin(LB), kmin(1000)
    nums["kmin_global"], nums["kmin_local"] = kg, kl
    rows = []
    for p in [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25]:
        lam = -math.log(1 - p)
        kstar = u_star / lam
        r21 = g(21 * lam) / gmax
        keg = max(round(kstar), kg)
        kel = max(round(kstar), kl)
        rows.append([f"{100*(1-p):g}", f"{kstar:.0f}", f2(r21, 3), str(keg), f2(g(keg * lam) / gmax, 3),
                     str(kel), f2(g(kel * lam) / gmax, 3), f2(g(kel * lam) / g(21 * lam), 2)])
    tbl(f"{gen}/tab_optk.tbl",
        ["c c c c s c s c", "c c c c c c c c"],
        "n n n n n n n n",
        [["ANI", "", "rel. info", f"global probe ($k >=$ {kg})", "", f"island-local probe ($k >=$ {kl})", "", "local $k$ vs."],
         ["(%)", "$k sup *$", "at $k$=21", "$k$", "rel. info", "$k$", "rel. info", "$k$=21"]],
        rows)

    # ---------------- Table: dense-floor efficiency vs k ----------------
    rows = []
    for k in [10, 12, 15, 21, 31]:
        r = [str(k)]
        for a in [0.999, 0.99, 0.95, 0.90, 0.85, 0.80]:
            r.append(f2(eta(1 - a, 1.0, 0.0, k=k), 3))
        rows.append(r)
    tbl(f"{gen}/tab_floork.tbl",
        ["c c s s s s s", "c c c c c c c"],
        "n n n n n n n",
        [["", "efficiency of the containment floor ($s$ = 1) at ANI =", "", "", "", "", ""],
         ["$k$", "99.9%", "99%", "95%", "90%", "85%", "80%"]],
        rows)

    # ---------------- Table: spurious matches ----------------
    rows = []
    for k in [9, 11, 13, 15, 17, 19, 21, 25, 31]:
        r = [str(k)]
        for N in [200, 1000, LB, 3e9]:
            r.append(sci(-math.expm1(-N / 4 ** k), 1))
        rows.append(r)
    tbl(f"{gen}/tab_spur.tbl",
        ["c c c c c", "c c c c c"],
        "n n n n n",
        [["", "island, $l$ = 200", "island, $l$ = 1000", "genome 5 Mbp", "genome 3 Gbp"],
         ["$k$", "(local probe)", "(local probe)", "(global probe)", "(global probe)"]],
        rows)

    # ---------------- Table: AF boundary resolution (Prop. 5) ----------------
    rows = []
    for s in [1 / 1000, 1 / 200, 1 / 50, 1 / 20, 1 / 5]:
        r = [f"1/{round(1/s)}"]
        for a in [0.99, 0.95, 0.90, 0.85]:
            th = theta(1 - a)
            r.append(f"{1/(s*th):,.0f}")
        rows.append(r)
    tbl(f"{gen}/tab_afres.tbl",
        ["c c s s s", "c c c c c"],
        "n n n n n",
        [["sampling rate", "mean boundary error (bp) at ANI =", "", "", ""],
         ["", "99%", "95%", "90%", "85%"]],
        rows)

    # ---------------- Table: Jensen / heterogeneity bias (Prop. 7) ----------------
    rng = np.random.default_rng(7)
    het = []

    def mix(fr, pr):
        fr = np.array(fr)
        pr = np.array(pr)
        return 1 - (fr * pr).sum(), lambda k: (fr * (1 - pr) ** k).sum() ** (1 / k)
    het.append(("two blocks: 50% at $p$=0, 50% at $p$=0.10 (recent HGT)", *mix([0.5, 0.5], [0.0, 0.10])))
    het.append(("two blocks: 50% at $p$=0.02, 50% at $p$=0.08", *mix([0.5, 0.5], [0.02, 0.08])))
    het.append(("90% at $p$=0.03, 10% hypervariable at $p$=0.23", *mix([0.9, 0.1], [0.03, 0.23])))
    for cv in [0.5, 1.0]:
        shape = 1 / cv ** 2
        scale = 0.05 / shape
        draws = np.clip(rng.gamma(shape, scale, 2_000_000), 0, 0.75)
        het.append((f"region rates ~ Gamma(mean 0.05, CV {cv:g})", 1 - draws.mean(),
                    (lambda d: (lambda k: np.mean((1 - d) ** k) ** (1 / k)))(draws)))
    rows = []
    for name, true_ani, fn in het:
        rows.append([name, f2(100 * true_ani, 2), f2(100 * fn(21), 2), f"{100*(fn(21)-true_ani):+.2f}",
                     f2(100 * fn(31), 2), f"{100*(fn(31)-true_ani):+.2f}"])
    tbl(f"{gen}/tab_jensen.tbl",
        ["l c c s c s", "l c c c c c"],
        "l n n n n n",
        [["rate profile (aligned regions)", "mean ANI", "$k$ = 21", "", "$k$ = 31", ""],
         ["", "(%)", "est. (%)", "bias (pp)", "est. (%)", "bias (pp)"]],
        rows)

    # ---------------- Table: island fingerprints (Prop. 8) ----------------
    def eta_if(p, sa, ka=21, f=8):
        lam = -math.log(1 - p)
        tha = (1 - p) ** ka
        anchor_info = ka ** 2 * math.exp(-lam * ka) / -math.expm1(-lam * ka)
        gmaxlen = int(60 / sa)
        gs = np.arange(1, gmaxlen + 1)
        pg = sa * (1 - sa) ** (gs - 1)
        ell = np.maximum(gs - ka, 0).astype(float)
        e2 = 2.0 ** -f
        t = np.exp(-lam * ell)
        tobs = t + (1 - t) * e2
        d = (1 - e2) * ell * t
        with np.errstate(divide="ignore", invalid="ignore"):
            info = np.where(ell > 0, d ** 2 / (tobs * (1 - tobs)), 0.0)
        I_lam = sa * (anchor_info + tha ** 2 * float((pg * info).sum()))
        return I_lam * p / (1 - p)

    def bits_if(sa, f=8):
        return sa * (32 + cpos(sa) + f)
    # equal-budget s_a
    lo, hi = 1e-5, 1e-1
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if bits_if(mid) > BETA0:
            hi = mid
        else:
            lo = mid
    sa_eq = lo
    nums["sa_eq"] = sa_eq
    rows = []
    for p in [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]:
        th = theta(p)
        e_star, _ = eps_opt(th, c_bf)
        s2 = BETA0 / c_bf(e_star)
        lam = -math.log(1 - p)
        # island length tuned to the information optimum u*/lambda, kept in the sparse regime
        sa_t = min(1 / 200, 1 / (u_star / lam + 21))
        rows.append([f"{100*(1-p):g}", f2(eta(p, S0, 0.0), 3), f2(eta(p, s2, e_star), 3),
                     f2(eta_if(p, sa_eq), 3), f"1/{round(1/sa_t)}", f2(eta_if(p, sa_t), 3), f2(bits_if(sa_t), 3),
                     f2((eta_if(p, sa_t) / bits_if(sa_t)) / (eta(p, s2, e_star) / BETA0), 2)])
    tbl(f"{gen}/tab_islandfp.tbl",
        ["c c c c c s s s", "c c c c c c c c"],
        "n n n n n n n n",
        [["ANI", "D0", "D2 (BF, $epsilon sup *$)", "D6 at D0 budget", "D6, mean island length $u sup * / lambda$ (max. information per bit)", "", "", ""],
         ["(%)", "$eta$", "$eta$", f"$eta$ ($s sub a$=1/{round(1/sa_eq)})", "$s sub a$", "$eta$", "bits/base", "per-bit gain vs D2"]],
        rows)

    # ---------------- Table: symmetric vs asymmetric SNR (Prop. 2 / Sec. 4.5) ----------------
    s2 = 1 / 50
    e = 0.01
    th = theta(0.05)
    nB = s2 * LB
    rows = []
    for LA in [2000, 5000, 20000, 100000]:
        N = s2 * LA
        sig = N * th
        v_samp = N * th * (1 - th) * (1 - s2)
        sd_asym = math.sqrt(v_samp + N * (1 - th) * e / (1 - e))
        sd_sym = math.sqrt(v_samp + 0.4427 * (nB + N) / math.log2(1 / e))
        sd_rev = math.sqrt(v_samp + s2 * LB * e * (1 - e))
        sd_rev_loc = math.sqrt(v_samp + s2 * 2 * LA * e * (1 - e))
        rows.append([kb(LA), f"{sig:.0f}", f2(sig / sd_asym, 1), f2(sig / sd_sym, 2), f2(sig / sd_rev, 2), f2(sig / sd_rev_loc, 1)])
    tbl(f"{gen}/tab_sym.tbl",
        ["c c c c c c", "c c c c c c"],
        "n n n n n n",
        [["contig", "signal", "SNR: contig probes", "SNR: BF\\(ca BF", "SNR: genome probes", "SNR: genome probes"],
         ["(kb)", "$s sub 2 n theta$", "genome BF", "(symmetric)", "contig BF (global)", "contig BF (localized)"]],
        rows)

    # ---------------- Figure 1: schematic ----------------
    draw_schematic(f"{figs}/fig_schematic.eps")

    # ---------------- bits/base needed to get within 1.5x of the s=1 floor ----------------
    for p in [0.01, 0.05, 0.10]:
        th = theta(p)
        q = 1 - th
        B = Bterm(p)
        best = None
        for e in EPS_GRID:
            phi = e / ((1 - e) * th)
            s_need = q * (1 + phi) / (0.5 * B + q)
            if s_need > 1:
                continue
            beta = c_bf(e) * s_need
            if best is None or beta < best[0]:
                best = (beta, e, s_need)
        nums[f"floor15_beta_p{int(round(p*100)):02d}"] = best[0]
        nums[f"floor15_eps_p{int(round(p*100)):02d}"] = best[1]
        nums[f"floor15_s_p{int(round(p*100)):02d}"] = best[2]
    nums["spur_3G_k21"] = -math.expm1(-3e9 / 4 ** 21)
    nums["pempty_5k"] = (1 - S0) ** (5000 - K + 1)
    nums["sd1k_D2opt_95"] = 100 * math.sqrt(var_ani(1000 - K + 1, 0.05, BETA0 / c_bf(nums["eps_opt_95"]), nums["eps_opt_95"]))

    # ---------------- numbers for text ----------------
    with open(f"{gen}/numbers.txt", "w") as fh:
        for kname, v in sorted(nums.items()):
            fh.write(f"{kname} = {v}\n")
    for kname, v in sorted(nums.items()):
        print(f"{kname} = {v}")


def draw_schematic(path):
    import matplotlib.patches as mp
    fig, ax = plt.subplots(figsize=(6.6, 2.75))
    ax.set_xlim(-2, 102)
    ax.set_ylim(-1.2, 7.4)
    ax.axis("off")
    yB, yA = 5.0, 1.6
    # reference B
    ax.add_patch(mp.Rectangle((0, yB - 0.18), 100, 0.36, fc="#efeeea", ec=INK2, lw=0.6))
    ax.text(-1.5, yB, "B", ha="right", va="center", fontsize=10, color=INK, style="italic")
    anchorsB = [6, 19, 27, 44, 58, 71, 86, 95]
    for j, x in enumerate(anchorsB):
        ax.add_patch(mp.Rectangle((x - 0.7, yB - 0.42), 1.4, 0.84, fc=C1, ec=INK, lw=0.5))
        ax.text(x, yB + 0.62, f"$a_{{{j+1}}}$", ha="center", va="bottom", fontsize=7, color=INK)
    for j in range(len(anchorsB) - 1):
        xm = (anchorsB[j] + anchorsB[j + 1]) / 2
        ax.annotate("", xy=(anchorsB[j] + 0.9, yB + 1.45), xytext=(anchorsB[j + 1] - 0.9, yB + 1.45),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.6))
        ax.text(xm, yB + 1.55, f"$I_{{{j+1}}}$", ha="center", va="bottom", fontsize=7, color=INK2)
    rng = np.random.default_rng(3)
    tierB = np.sort(rng.uniform(0, 100, 34))
    tierB = [t for t in tierB if min(abs(t - a) for a in anchorsB) > 1.2]
    ax.plot(tierB, [yB - 0.75] * len(tierB), ls="none", marker="o", ms=4.2, color=C2)
    ax.text(101, yB - 0.75, "tier-2 $k$-mers\n($h<s_2$) $\\rightarrow$\nBF keys $(x, j)$", fontsize=6.3,
            color=INK, va="center", ha="left")
    # query A: covers B[30, 90] aligned + 12 units unaligned prefix
    offset = -14.0   # A coordinate = B coordinate + offset in aligned region
    a_start, a_end = 16, 76
    ax.add_patch(mp.Rectangle((a_start, yA - 0.18), a_end - a_start, 0.36, fc="#efeeea", ec=INK2, lw=0.6))
    ax.add_patch(mp.Rectangle((a_start, yA - 0.18), 14, 0.36, fc="#ffffff", ec=INK2, lw=0.6, hatch="////"))
    ax.text(a_start - 1.5, yA, "A", ha="right", va="center", fontsize=10, color=INK, style="italic")
    ax.text(a_start + 7, yA - 0.95, "unaligned", ha="center", fontsize=6.5, color=INK2)
    # aligned region maps B[44-? ...]: A x in [30,76] <-> B x - offset in [44, 90]
    mutated = {71}
    for x in anchorsB:
        xa = x + offset
        if 30 <= xa <= a_end:
            if x in mutated:
                ax.plot([xa], [yA], marker="x", ms=7, color=C5, mew=1.6)
                ax.text(xa, yA - 0.95, "anchor\nmutated", ha="center", fontsize=6, color=INK2)
            else:
                ax.add_patch(mp.Rectangle((xa - 0.7, yA - 0.42), 1.4, 0.84, fc=C1, ec=INK, lw=0.5))
                ax.plot([xa, x], [yA + 0.45, yB - 0.45], color=C1, lw=0.8, alpha=0.8)
    tierA = np.sort(rng.uniform(a_start, a_end, 30))
    for t in tierA:
        hit = (t >= 30) and (rng.random() < 0.6)
        ax.plot([t], [yA + 0.75], marker="o", ms=4.2, color=C2, mfc=(C2 if hit else "white"), mew=1.0)
    ax.text(a_end + 1.0, yA + 0.75, "query tier-2 $k$-mers: filled = BF hit,\nhollow = miss "
            "(probe key = ($y$, island of $B$))", fontsize=6.3, color=INK, va="center")
    ax.annotate("", xy=(30, yA - 1.45), xytext=(a_end, yA - 1.45), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.6))
    ax.text((30 + a_end) / 2, yA - 1.4, "aligned span (AF = 46/60)", ha="center", va="bottom", fontsize=6.5, color=INK2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
