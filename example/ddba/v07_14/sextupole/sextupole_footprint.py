# -*- coding: utf-8 -*-
"""チューンの振幅依存 (追跡 + NAFF) と、解析式の ADTS・高次の ν(J) の比較、DA の端のチューンずれ (R7)。
使い方: py -3.12 sextupole_footprint.py [NAME ...] [--plot-only]    出力: results/footprint.json、results/footprint.out.txt、fig_footprint.png
条件: RF・放射なし、リング先頭から px = py = 0 で出発、1024 周 (前半・後半 512 周のチューン = FrequencyMap の nux1・nuy1)。
  x の枝: x = 0.2〜6 mm (16 点)、y = 0.2 mm。y の枝: y = 0.2〜5 mm (16 点)、x = 0.2 mm。生き残った点だけ。
小振幅の ν(x) = ν0 + a x² + b x⁴ をフィットして (x < 2 mm の点)、a を解析式の ∂νx/∂Jx / (2βx) と比べる (J = x²/(2β))。b ≠ 0 は ADTS の線形の近似を超える (J² の項: 5 次の写像の項)。
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
C.lower_priority()
from pysad import sxtopt as S
from pysad.freqmap import frequency_map

PLOT_ONLY = '--plot-only' in sys.argv
names = [a for a in sys.argv[1:] if not a.startswith('--')] or ['final', 'rev9', 'retune1', 'opt_full_1', 'opt_full_k125', 'opt_adts', 'v07_12']
CAND = json.load(open(os.path.join(C.OUT, 'candidates.json'), encoding='utf-8'))
src0 = C.load_source()
RDT = json.load(open(os.path.join(C.OUT, 'rdt_table.json'), encoding='utf-8'))
out = {}
LOG = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


for n in ([] if PLOT_ONLY else names):
    k = np.array(CAND['candidates'][n]['k'])
    src = S.replace_k2_in_source(src0, CAND['family_names'], k)
    st = C.make_track_setup(src, rad=False, rfsw=False, alost=0.05)
    bx, by = st.rows[0][1], st.rows[0][3]
    xs = np.linspace(0.2e-3, 6e-3, 16)
    ys = np.linspace(0.2e-3, 5e-3, 16)
    rx = frequency_map(st, 1, x=xs, y=[0.2e-3], nturn=1024, workers=1, warn=False)
    ry = frequency_map(st, 1, x=[0.2e-3], y=ys, nturn=1024, workers=1, warn=False)
    ok_x = (rx.lost_turn == 0) & np.isfinite(rx.nux1)
    ok_y = (ry.lost_turn == 0) & np.isfinite(ry.nuy1)
    nu0 = st.tunes_full()
    d = dict(x=xs.tolist(), nux_x=rx.nux1.tolist(), nuy_x=rx.nuy1.tolist(), ok_x=ok_x.tolist(), y=ys.tolist(), nux_y=ry.nux1.tolist(), nuy_y=ry.nuy1.tolist(), ok_y=ok_y.tolist(),
             nu0=[float(v % 1) for v in nu0], beta=[bx, by])
    sel = ok_x & (xs < 2.5e-3)
    if sel.sum() >= 4:
        A = np.vstack([xs[sel] ** 2, xs[sel] ** 4]).T
        a, b = np.linalg.lstsq(A, rx.nux1[sel] - rx.nux1[sel][0] + 0 * xs[sel], rcond=None)[0]
        # nux1[sel][0] は x = 0.2 mm の値 ≈ ν0 + a x0²。ν0 を未知数として入れ直す
        A2 = np.vstack([np.ones(sel.sum()), xs[sel] ** 2, xs[sel] ** 4]).T
        c0, a, b = np.linalg.lstsq(A2, rx.nux1[sel], rcond=None)[0]
        an = RDT[n]['adts_nform'][0] if n in RDT and 'adts_nform' in RDT[n] else None
        d['fit_x'] = dict(nu0=float(c0), a=float(a), b=float(b), a_adts=(an / (2 * bx) if an is not None else None))
    out[n] = d
    last_x = xs[ok_x][-1] if ok_x.any() else 0.0
    P('%-14s 生き残った x の最大 %.1f mm (x の枝)、y の最大 %.1f mm (y の枝)、νx の x² 係数 a = %s (解析の ADTS から %s)、x⁴ 係数 b = %s' % (
        n, last_x * 1e3, (ys[ok_y][-1] if ok_y.any() else 0) * 1e3,
        '%.1f' % d['fit_x']['a'] if 'fit_x' in d else '-', '%.1f' % d['fit_x']['a_adts'] if 'fit_x' in d and d['fit_x']['a_adts'] is not None else '-',
        '%.3g' % d['fit_x']['b'] if 'fit_x' in d else '-'))
if not PLOT_ONLY:
    json.dump(out, open(os.path.join(C.OUT, 'footprint.json'), 'w', encoding='utf-8'), ensure_ascii=False)
    open(os.path.join(C.OUT, 'footprint.out.txt'), 'w', encoding='utf-8').write('\n'.join(LOG) + '\n')
def make_plot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    show = [n for n in ('v07_12', 'retune1', 'final', 'rev9', 'opt_full_1', 'opt_full_k125', 'opt_adts') if n in out]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for n in show:
        d = out[n]
        ok = np.array(d['ok_x'])
        x = np.array(d['x'])[ok] * 1e3
        nux = np.array(d['nux_x'])[ok]
        nuy = np.array(d['nuy_x'])[ok]
        dnx = ((nux - nux[0] + 0.5) % 1) - 0.5
        dny = ((nuy - nuy[0] + 0.5) % 1) - 0.5
        axs[0].plot(x, dnx, 'o-', ms=3, lw=1, label=n)
        okm = np.abs(dny) < 0.05
        axs[1].plot(x[okm], dny[okm], 'o-', ms=3, lw=1, label=n)
        ok = np.array(d['ok_y'])
        y = np.array(d['y'])[ok] * 1e3
        nuy = np.array(d['nuy_y'])[ok]
        dyy = ((nuy - nuy[0] + 0.5) % 1) - 0.5
        axs[2].plot(y, dyy, 'o-', ms=3, lw=1, label=n)
    for ax, t, xl in zip(axs, ('Δνx (x の枝、y = 0.2 mm)', 'Δνy (x の枝、y = 0.2 mm、|Δ| < 0.05 の点)', 'Δνy (y の枝、x = 0.2 mm)'), ('x [mm]', 'x [mm]', 'y [mm]')):
        ax.set_title(t, fontsize=9)
        ax.set_xlabel(xl)
        ax.set_ylabel('チューンのずれ (最小振幅の点から)')
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.suptitle('チューンの振幅依存 (追跡 + NAFF、1024 周、RF・放射なし、リング先頭から)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_footprint.png'), dpi=120)


if PLOT_ONLY:
    make_plot(json.load(open(os.path.join(C.OUT, 'footprint.json'), encoding='utf-8')))
else:
    make_plot(out)
