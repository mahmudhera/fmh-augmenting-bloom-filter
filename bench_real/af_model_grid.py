#!/usr/bin/env python3
"""AF-model settings on the tuning sets (real data/real_tune and simulated data/tune).

For each setting (aif query options) reports, on the TUNING sets only:
  real contig ANI MAE and AF MAE (macro over the 5 sets), aggregated genome-level ANI MAE vs ANIm,
  and simulated ANI / AF MAE (macro over conditions). Writes data/real/results/af_model_grid.tsv.

Usage: python3 bench_real/af_model_grid.py   (settings below; simulated sketches from data/results/main/aif)
"""
import csv, math, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
R = 'data/real_tune'; SETS = ['bcereus', 'ecoli', 'mag2328', 'paer', 'saur']
SIMD = 'data/results/main/aif'
jobs = []
for s in SETS:
    refs = sorted(set(l.split('\t')[4] for l in list(open(f'{R}/sets/{s}/genomes.tsv'))[1:]))
    jobs += [(s, r) for r in refs]
L, T, GT = {}, {}, {}
for s in SETS:
    for r in csv.DictReader(open(f'{R}/sets/{s}/truth.tsv'), delimiter='\t'):
        L[r['contig']] = int(r['length']); T[r['contig']] = (float(r['true_af']), float(r['true_ani']), s)
    for r in csv.DictReader(open(f'{R}/sets/{s}/genome_truth.tsv'), delimiter='\t'):
        GT[(s, r['gid'])] = float(r['anim_ani'])
simtruth = {r['contig']: (float(r['true_af']), float(r['true_ani']), r['condition']) for r in csv.DictReader(open('data/tune/truth.tsv'), delimiter='\t')}
def run_real(a):
    (s, ref), args = a
    return s, subprocess.run(['bin/aif', 'query', f'{R}/results/aif/{s}.{ref}.aif', f'{R}/sets/{s}/{ref}.contigs.fa', '-t', '2'] + args, capture_output=True, text=True, check=True).stdout
def run_sim(a):
    g, args = a
    return subprocess.run(['bin/aif', 'query', f'{SIMD}/{g}.aif', f'data/tune/{g}.contigs.fa', '-t', '2'] + args, capture_output=True, text=True, check=True).stdout
def evaluate(args, ex):
    pred = {}
    for s, out in ex.map(run_real, [(j, args) for j in jobs]):
        for line in out.splitlines()[1:]:
            p = line.split('\t'); ok = p[7] == 'ok'
            pred[p[0]] = (float(p[1]) if ok else math.nan, float(p[2]) if ok else 0.0, ok)
    ani, af = [], []
    for s in SETS:
        cs = [c for c in T if T[c][2] == s]
        a = [abs(pred[c][0] - T[c][1]) * 100 for c in cs if T[c][0] >= 0.2 and pred[c][2]]
        f = [abs(pred[c][1] - T[c][0]) * 100 for c in cs]
        ani.append(sum(a) / len(a)); af.append(sum(f) / len(f))
    acc = {}
    for c, (a, f, ok) in pred.items():
        if not ok: continue
        g = (T[c][2], c.split('__')[0]); x = acc.setdefault(g, [0.0, 0.0]); al = f * L[c]; x[0] += al; x[1] += al * a
    gerr = [abs(x[1] / x[0] - GT[g]) * 100 for g, x in acc.items() if x[0] > 0]
    sp = {}
    for out in ex.map(run_sim, [(g, args) for g in ['NC_000913.3', 'NC_000964.3', 'NC_002516.2', 'NC_007795.1']]):
        for line in out.splitlines()[1:]:
            p = line.split('\t')
            if p[7] == 'ok': sp[p[0].split()[0]] = (float(p[1]), float(p[2]))
    conds = {}
    for c, (taf, tani, cond) in simtruth.items():
        if c in sp: conds.setdefault(cond, []).append((abs(sp[c][0] - tani) * 100, abs(sp[c][1] - taf) * 100))
    sa = sum(sum(x[0] for x in v) / len(v) for v in conds.values()) / len(conds)
    sf = sum(sum(x[1] for x in v) / len(v) for v in conds.values()) / len(conds)
    return sum(ani) / 5, sum(af) / 5, sum(gerr) / len(gerr), sa, sf

SETTINGS = [
    ("default", ""),
    ("switching rate 1e-5", "--rho-u 1e-5"), ("switching rate 2e-5 (aif-RT)", "--rho-u 2e-5"), ("switching rate 3e-5", "--rho-u 3e-5"),
    ("contig window 0.06", "--ident-window 0.06 --ident-radius 10000"), ("contig window 0.08", "--ident-window 0.08 --ident-radius 10000"),
    ("contig window 0.10", "--ident-window 0.10 --ident-radius 10000"),
    ("local prior 0.08, R=1", "--ident-window 0.08 --ident-radius 1"), ("local prior 0.08, R=3", "--ident-window 0.08 --ident-radius 3"),
    ("local prior 0.10, R=1", "--ident-window 0.10 --ident-radius 1"), ("local prior 0.10, R=3", "--ident-window 0.10 --ident-radius 3"),
    ("overdispersion 1", "--overdisp 1"), ("overdispersion 1, rate 1e-5", "--overdisp 1 --rho-u 1e-5"),
]
grid = [a.split() for _, a in SETTINGS]
rows = []
with ThreadPoolExecutor(16) as ex:
    for args in grid:
        r = evaluate(args, ex)
        rows.append(r)
        print(f"{' '.join(args):40s} real ANI {r[0]:.4f} AF {r[1]:.3f} genomeANI {r[2]:.4f} | sim ANI {r[3]:.4f} AF {r[4]:.3f}", flush=True)
with open('data/real/results/af_model_grid.tsv', 'w') as fh:
    fh.write('setting\toptions\treal_ani_mae\treal_af_mae\tgenome_ani_mae\tsim_ani_mae\tsim_af_mae\n')
    for (name, a), r in zip(SETTINGS, rows):
        fh.write(f"{name}\t{a}\t" + "\t".join(f"{x:.4f}" for x in r) + "\n")
