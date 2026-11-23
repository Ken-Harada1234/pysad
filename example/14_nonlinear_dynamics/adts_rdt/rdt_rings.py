# -*- coding: utf-8 -*-
"""共鳴駆動項 (6 極の 1 次)・ADTS (6 極の 2 次 + 8 極の 1 次)・クロマティシティを 3 つのリングで計算する Python の例。
リング: (1) 小さな FODO (example/00_lattice/FODO_ring.sad に SF・SD の K2 = 4, −6 を入れたもの、3 GeV、105.6 m)
        (2) ASTRID2 (example/astrid/lattice/ASTRID2_v01_01.sad、0.58 GeV、45.7 m)
        (3) 大きなリング A (pysad/tests/kek-ls_v03_12.lat、3 GeV、603 m、3563 要素、誤差なし)
条件: δ = 0 の設計光学 (CALC)、結合・垂直分散なし、6 極は厚いレンズの解析式 (method='analytic')、数値積分 (nslice=16) と照合。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_rings.py            出力: 標準出力 (rdt_rings.pysad.txt に保存)、図 rdt_rings_elements.png"""
import os, sys, math
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, ROOT)
from pysad import resdrive as R, compat
from pysad.parser import load_sad_file, load_sad_string
compat.set_mode('v2')
HERE = os.path.dirname(os.path.abspath(__file__))
fodo_src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
RINGS = [('小さな FODO (3 GeV、105.6 m)', lambda: load_sad_string(fodo_src).build_lattice('RING')),
         ('ASTRID2 (0.58 GeV、45.7 m)', lambda: load_sad_file(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')).build_lattice('RING')),
         ('大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)', lambda: load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat')).build_lattice('RING'))]
fig_data = {}
for name, mk in RINGS:
    lat = mk()
    rows = lat.twiss_table()
    print('=' * 100)
    print(name, '  素子数', len(lat.elements), '  周長 %.3f m' % rows[-1].s)
    wp = R.working_point_check(rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi)
    print('作動点の確認: νx = %.6f νy = %.6f; 整数共鳴まで x %.4f y %.4f; 半整数まで x %.4f y %.4f; 3 次までで最も近い共鳴線 %s' % (
        rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi, wp['dist_x_int'], wp['dist_y_int'], wp['dist_x_half'], wp['dist_y_half'],
        ['%d·νx%+d·νy: %.4f' % t for t in wp['nearest_lines'][:3]]))
    h = R.resonance_terms(lat, rows=rows, bend_cubic=False)      # この図は 6 極だけ (既定は偏向磁石の項つき)
    print(h)
    h2 = R.resonance_terms(lat, method='numeric', nslice=16, rows=rows, bend_cubic=False)
    print('  数値積分 (nslice=16) との差の最大 |Δh| / |h| = %.2e' % max(abs(h.total[k] - h2.total[k]) / abs(h.total[k]) for k in h.names))
    for nm in ('h30000', 'h10110'):
        print('  要素ごとの寄与の大きい 5 つ (%s): ' % nm + ', '.join('%s(s=%.1f) %.3f' % (n, s, abs(v)) for n, s, v in h.element_table(nm, 5)))
    a = R.adts(lat, rows=rows)
    print(a)
    an = R.adts(lat, method='numeric', nslice=16, rows=rows)
    print('  数値積分 (nslice=16) の 6 極: %s' % ', '.join('%.6e' % an.sext[k] for k in ('xx', 'xy', 'yy')))
    print('  6 極の調和ごとの寄与 (k=(kx,ky) → (xx, xy, yy)):')
    for kv, v in sorted(a.by_harmonic.items()):
        if any(abs(t) > 0 for t in v) and (kv[0] > 0 or (kv[0] == 0 and kv[1] > 0)):
            print('     k = %-8s  νx 側の共鳴 %d·νx%+d·νy:  %s' % (kv, kv[0], kv[1], ', '.join('%.4e' % (2 * t if False else t) for t in v)))
    # ∂ν/∂δ・∂²ν/∂δ²・∂³ν/∂δ³ (CALC の δ 依存の数値微分) と 6 極の 1 次の寄与
    c = R.chromaticity_numeric(lat, order=3, dmax=2e-3)
    sx, sy = R.chromaticity_sextupole(lat, rows=rows)
    print('クロマティシティ (CALC の δ = ±2e-3 の範囲の多項式フィット): ∂ν/∂δ = (%.4f, %.4f)  ∂²ν/∂δ² = (%.2f, %.2f)  ∂³ν/∂δ³ = (%.0f, %.0f)' % (
        c['d1'] + c['d2'] + c['d3']))
    print('  6 極の 1 次の寄与 (h11001/π, h00111/π) = (%.4f, %.4f)  (1 次の分散だけ。残りは自然クロマティシティ・高次)' % (sx, sy))
    # δ ≠ 0 の ADTS (δ のオフモメンタムの光学で): 一次の δ 依存 ∂(∂ν/∂J)/∂δ の中心差分
    dd = 1e-3
    ap, am = R.adts(lat, delta=dd), R.adts(lat, delta=-dd)
    print('  ADTS の δ 依存 (6 極 2 次、δ = ±1e-3 の光学): ∂(dνx/dJx)/∂δ = %.4e  ∂(dνx/dJy)/∂δ = %.4e  ∂(dνy/dJy)/∂δ = %.4e' % tuple(
        (ap.sext[k] - am.sext[k]) / (2 * dd) for k in ('xx', 'xy', 'yy')))
    fig_data[name] = (h, a, rows)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(RINGS), 2, figsize=(11, 3.2 * len(RINGS)))
    for i, (name, (h, a, rows)) in enumerate(fig_data.items()):
        k = h.names.index('h30000')
        ax[i, 0].bar(h.elem_s, np.abs(h.per_element[:, k]), width=max(1.0, rows[-1].s / 400))
        ax[i, 0].set_title('%s: |h30000| per sextupole' % ('FODO', 'ASTRID2', 'ring A (kek-ls_v03_12)')[i], fontsize=9)
        ax[i, 0].set_xlabel('s [m]')
        for j, lab in enumerate(('dnux/dJx', 'dnux/dJy', 'dnuy/dJy')):
            ax[i, 1].plot(a.elem_s, np.cumsum(a.per_element[j]), label=lab)
        ax[i, 1].set_title('cumulative sextupole 2nd-order ADTS contribution', fontsize=9)
        ax[i, 1].set_xlabel('s [m]'); ax[i, 1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'rdt_rings_elements.png'), dpi=110)
    print('図: rdt_rings_elements.png')
except Exception as e:
    print('図は作らなかった:', e)
