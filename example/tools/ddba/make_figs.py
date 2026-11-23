# -*- coding: utf-8 -*-
"""DDBA_NN の εx・チューン・周長・自然クロマティシティの図 (example/ddba/DDBA_overview.png)。
結果は example/ddba/results_k64.json (pysad k64) と chrom_k64.json。段階ごとに色を変える。"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'tools', 'figures'))
import figsave  # 折れ線系は SVG、そうでなければ PNG
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE)); D = os.path.join(EX, 'ddba')
sys.path.insert(0, HERE); import ddba_info as I
plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'DejaVu Sans']
R = json.load(open(os.path.join(D, 'results_k64.json'))); C = json.load(open(os.path.join(D, 'chrom_k64.json')))
cols = {1: '#4c78a8', 2: '#f58518', 3: '#54a24b', 4: '#e45756', 5: '#72b7b2', 6: '#b279a2', 7: '#9d755d', 8: '#bab0ac', 9: '#222222'}
ns = [n for n in range(1, 29) if R[str(n)].get('C')]
col = [cols[I.I[n]['STAGE']] for n in ns]
fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
a = ax[0, 0]; a.scatter(ns, [R[str(n)]['ex'] * 1e9 * (2.5 / (R[str(n)]['P'] / 1e9)) ** 2 for n in ns], c=col, s=36)
a.set_ylabel(r'$\varepsilon_x$ (scaled to 2.5 GeV) [nm rad]'); a.set_yscale('log'); a.set_title('エミッタンス (pysad k64。エネルギーの違いは E² で 2.5 GeV に換算)')
a = ax[0, 1]; a.scatter(ns, [R[str(n)]['C'] for n in ns], c=col, s=36); a.set_ylabel('周長 [m]'); a.set_title('周長')
a = ax[1, 0]
a.scatter(ns, [R[str(n)]['nx'] for n in ns], c=col, s=36, marker='o', label=r'$\nu_x$'); a.scatter(ns, [R[str(n)]['ny'] for n in ns], c=col, s=36, marker='s', label=r'$\nu_y$')
a.set_ylabel('チューン'); a.legend(); a.set_title('チューン (丸 νx、四角 νy)')
a = ax[1, 1]; nc = [n for n in ns if C.get(str(n))]
a.scatter(nc, [C[str(n)]['nx'] for n in nc], c=[cols[I.I[n]['STAGE']] for n in nc], s=36, marker='o'); a.scatter(nc, [C[str(n)]['ny'] for n in nc], c=[cols[I.I[n]['STAGE']] for n in nc], s=36, marker='s')
a.set_ylabel('自然クロマティシティ (6 極 K2 = 0)'); a.set_title('自然クロマティシティ (丸 x、四角 y)')
for x in ax.flat: x.set_xlabel('DDBA_NN'); x.grid(alpha=.3)
h = [plt.Line2D([], [], marker='o', ls='', color=cols[k], label='%d. %s' % (k, v.split(' (')[0])) for k, v in I.STAGES.items()]
fig.legend(handles=h, loc='lower center', ncol=3, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.1, 1, 1)); figsave.save_figure(fig, os.path.join(D, 'DDBA_overview.png'), dpi=130, bbox_inches=None)
