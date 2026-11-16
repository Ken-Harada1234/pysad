# -*- coding: utf-8 -*-
"""gs の出力 (.gs.txt) の数値から PNG を作る。py -3.12 tools/make_figs.py"""
import os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EX = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'


def lines(rel, kind='gs'):
    return open(os.path.join(EX, rel + '.%s.txt' % kind), encoding='utf-8', errors='replace').read().split('\n')


def save(fig, rel):
    fig.tight_layout()
    fig.savefig(os.path.join(EX, rel), dpi=80)
    plt.close(fig)
    print(rel, os.path.getsize(os.path.join(EX, rel)))


# ---- T1: optics (OpticsPlot data)
L = lines('T1_intro/FFS_draw_OpticsPlot_PF')
i0 = [k for k, l in enumerate(L) if l.startswith('PLOTDATA')][0]
rows = []
for l in L[i0 + 2:]:
    m = re.match(r'\s*(%s)\s+(%s)\s+(%s)\s+(%s)\s+(\S+)\s*$' % (num, num, num, num), l)
    if m:
        rows.append([float(m.group(k)) for k in range(1, 5)])
a = np.array(rows)
fig, ax = plt.subplots(figsize=(8, 3.2))
ax.plot(a[:, 0], a[:, 1], label='BX'); ax.plot(a[:, 0], a[:, 2], label='BY')
ax.plot(a[:, 0], a[:, 3] * 10, label='EX x10')
ax.set_xlabel('s [m]'); ax.set_ylabel('BX, BY [m],  EX x10 [m]'); ax.set_title('PF 2.5 GeV (v9_16_asym): optics from the gs output'); ax.legend(ncol=3, loc='upper right')
ax.set_xlim(0, a[:, 0].max())
save(fig, 'T1_intro/FFS_draw_OpticsPlot_PF_optics.png')

# ---- T1: tune vs dp/p
L = lines('T1_intro/fn_PolynomialFit_PF')
d = []
for l in L:
    m = re.match(r'\s*(%s)\s+(%s)\s+(%s)\s+(%s)\s*$' % (num, num, num, num), l)
    if m and abs(float(m.group(1)) - round(float(m.group(1)))) < 1e-9 and 1 <= float(m.group(1)) <= 41:
        d.append([float(m.group(k)) for k in range(2, 5)])
d = np.array(d[:41])
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.plot(d[:, 0] * 100, d[:, 1], 'o-', ms=3, label='nu_x')
ax.plot(d[:, 0] * 100, d[:, 2] + 4.3, 'o-', ms=3, label='nu_y + 4.3')
ax.set_xlabel('dp/p [%]'); ax.set_ylabel('tune'); ax.set_title('PF tune vs dp/p (design sextupoles)'); ax.legend()
save(fig, 'T1_intro/fn_PolynomialFit_PF_tune_vs_dp.png')

# ---- T2: dispersion at PM*
L = lines('T2_samples/calc_dispersion_PF')
rows = []
for l in L:
    m = re.match(r'\s*(\d+)\.0+\s+(%s)\s+(%s)\s*$' % (num, num), l)
    if m:
        rows.append([float(m.group(2)), float(m.group(3))])
a = np.array(rows)
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(a[:, 0], a[:, 1], 'o-', ms=3)
ax.set_xlabel('s [m]'); ax.set_ylabel('EX at PM* [m]'); ax.set_title('PF dispersion at the 65 monitors (gs)')
save(fig, 'T2_samples/calc_dispersion_PF_ex.png')

# ---- T2: mom_dep
L = lines('T2_samples/FFS_mom_dep_ts_PF')
blocks = {0: [], 1: []}; b = -1
for l in L:
    if l.startswith('--- with the design'): b = 0
    elif l.startswith('--- without'): b = 1
    elif b >= 0:
        m = re.match(r'\s*(%s)\s+(%s)\s+(%s)\s*$' % (num, num, num), l)
        if m: blocks[b].append([float(m.group(k)) for k in (1, 2, 3)])
fig, axs = plt.subplots(1, 2, figsize=(8, 3.4))
for k, ttl in ((0, 'design sextupoles'), (1, 'S* 0 (natural)')):
    a = np.array(blocks[k])
    axs[k].plot(a[:, 0] * 100, a[:, 1], 'o-', ms=3, label='nu_x'); axs[k].plot(a[:, 0] * 100, a[:, 2], 's-', ms=3, label='nu_y')
    axs[k].set_title(ttl); axs[k].set_xlabel('dp/p [%]'); axs[k].legend()
axs[0].set_ylabel('tune')
save(fig, 'T2_samples/FFS_mom_dep_ts_PF_tune_vs_dp.png')

# ---- T2: amp dep
L = lines('T2_samples/FFS_amp_dep_ts_PF')
X, Y = [], []; mode = None
for l in L:
    if l.startswith('--- horizontal'): mode = X
    elif l.startswith('--- vertical'): mode = Y
    elif mode is not None:
        m = re.match(r'\s*(%s)\s+(%s)\s+(%s)\s*$' % (num, num, num), l)
        if m: mode.append([float(m.group(k)) for k in (1, 2, 3)])
nx0, ny0 = 9.5975304909, 5.2789602116
X = np.array(X); Y = np.array(Y)
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.plot(X[:, 0] ** 2, X[:, 1] - nx0, 'o-', label='dnu_x (x amplitude)'); ax.plot(X[:, 0] ** 2, X[:, 2] - ny0, 's-', label='dnu_y (x amplitude)')
ax.plot(Y[:, 0] ** 2, Y[:, 1] - nx0, 'o--', label='dnu_x (y amplitude)'); ax.plot(Y[:, 0] ** 2, Y[:, 2] - ny0, 's--', label='dnu_y (y amplitude)')
ax.set_xlabel('amplitude^2 [sigma^2]'); ax.set_ylabel('tune shift'); ax.legend(fontsize=7); ax.set_title('PF amplitude-dependent tune shift (gs)')
save(fig, 'T2_samples/FFS_amp_dep_ts_PF_dnu_vs_A2.png')

# ---- T2: phase space
L = lines('T2_samples/phase_PF')
rec = {}
for l in L:
    m = re.match(r'amp\s+(%s) sigma\s+turn\s+(%s)\s+x \[mm\]\s+(%s)\s+px \[mrad\]\s+(%s)' % (num, num, num, num), l)
    if m:
        rec.setdefault(int(float(m.group(1))), []).append([float(m.group(3)), float(m.group(4))])
fig, ax = plt.subplots(figsize=(5, 4))
for k, v in sorted(rec.items()):
    v = np.array(v); ax.plot(v[:, 0], v[:, 1], '.', ms=5, label='%d sigma' % k)
ax.set_xlabel('x [mm]'); ax.set_ylabel("x' [mrad]"); ax.legend(); ax.set_title('PF one-turn x-x\' at the ring start (60 turns, gs)')
save(fig, 'T2_samples/phase_PF_xpx.png')

# ---- T2: da_xy
L = lines('T2_samples/fn_DynamicApertureSurvey_da_xy_PF')
blk = {}; cur = None
for l in L:
    m = re.match(r'dp/p =\s+(%s)' % num, l)
    if m: cur = round(float(m.group(1)), 1); blk[cur] = []; continue
    m = re.match(r'\s*(%s)\s+(%s)\s*$' % (num, num), l)
    if m and cur is not None: blk[cur].append([float(m.group(1)), float(m.group(2))])
fig, ax = plt.subplots(figsize=(5, 4.6))
for k, v in sorted(blk.items()):
    v = np.array(v); v = np.vstack([v, v[:1]]); ax.plot(v[:, 0], v[:, 1], 'o-', ms=3, label='dp = %+.1f sigma_E' % k)
ax.set_xlabel('x [sigma_x]'); ax.set_ylabel('y [sigma_y]'); ax.set_aspect('equal'); ax.legend(fontsize=7); ax.set_title('PF DA in x-y, 1000 turns (gs)')
save(fig, 'T2_samples/fn_DynamicApertureSurvey_da_xy_PF_da.png')

# ---- T2: topdrawer file -> PNG, and the .td text
L = lines('T2_samples/FFS_draw_topdrawer_PF')
i0 = [k for k, l in enumerate(L) if l.strip() == 'TDBEGIN'][0]; i1 = [k for k, l in enumerate(L) if l.strip() == 'TDEND'][0]
body = [l.strip() for l in L[i0 + 1:i1] if 'Print[' not in l and 'Do[' not in l and l.strip()]
open(os.path.join(EX, 'T2_samples/FFS_draw_topdrawer_PF_fort12.td'), 'w', encoding='utf-8', newline='\n').write('\n'.join(body) + '\n')
curves = [[]]; names = []
for l in body:
    m = re.match(r'(%s)\s+(%s)$' % (num, num), l)
    if m: curves[-1].append([float(m.group(1)), float(m.group(2))]); continue
    m = re.match(r"JOIN 1 TEXT '(.*)'", l)
    if m: names.append(m.group(1)); curves.append([])
fig, ax = plt.subplots(figsize=(8, 3.2))
for c, n in zip(curves, names):
    c = np.array(c); ax.plot(c[:, 0], c[:, 1] * (10 if n == 'EX' else 1), label=n + (' x10' if n == 'EX' else ''))
ax.set_xlabel('s [m]'); ax.legend(ncol=3); ax.set_title('the hand-written TopDrawer data, redrawn with matplotlib')
save(fig, 'T2_samples/FFS_draw_topdrawer_PF_plot.png')
