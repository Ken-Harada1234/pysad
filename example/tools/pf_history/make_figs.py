# -*- coding: utf-8 -*-
"""PF の履歴ラティス 4 本の図を作る (pysad k64 の Twiss から)。
  example/pf/lattice/PF_history_optics.png : リング 1 周の βx βy ηx と、偏向磁石・四極の並び
  example/pf/lattice/PF_history_cell.png   : 標準セル (偏向磁石 2 個分、8.59 m) の拡大
  example/pf/lattice/PF_history_layout.png : リングの配置図 (長直線部 0.55 pt、その他の基準軌道 1.8 pt、偏向磁石 3 pt。作図ルールによる)
使い方: py -3.12 example/tools/pf_history/make_figs.py"""
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pf_lib
OUT = os.path.join(pf_lib.EX, 'pf', 'lattice')
CB, CY, CE = '#0072B2', '#D55E00', '#009E73'
T = lambda e: type(e).__name__
TITLE = {'pf_1981_460nmrad': '1981  pf_1981_460nmrad', 'pf_1986_128nmrad': '1986  pf_1986_128nmrad',
         'pf_1997_036nmrad': '1997  pf_1997_036nmrad', 'pf_2005_036nmrad': '2005  pf_2005_036nmrad'}
for _f in ('Meiryo', 'Yu Gothic', 'Noto Sans CJK JP', 'IPAGothic', 'MS Gothic'):
    if any(_f == f.name for f in matplotlib.font_manager.fontManager.ttflist):
        matplotlib.rcParams['font.family'] = _f
        break
matplotlib.rcParams['axes.unicode_minus'] = False


def is_mag(e):
    t = T(e)
    return (t == 'Bend' and getattr(e, 'ANGLE', 0) != 0) or (t == 'Quadrupole' and e.L > 0 and e.K1 != 0) or (t == 'Sextupole' and e.L > 0 and e.K2 != 0)


def strip(ax, els, rows, y0, h, smin=None, smax=None):
    """磁石の並び: 偏向磁石は黒、収束四極 (K1>0) は上、発散四極は下、6 極は緑の細い箱、空洞は赤い三角"""
    for i, e in enumerate(els):
        s0 = rows[i].s; s1 = s0 + e.L
        if smin is not None and (s1 < smin or s0 > smax):
            continue
        t = T(e)
        if t == 'Bend' and getattr(e, 'ANGLE', 0) != 0:
            ax.add_patch(Rectangle((s0, y0 - h / 2), e.L, h, fc='#333333', ec='none'))
        elif t == 'Quadrupole' and e.L > 0 and e.K1 != 0:
            ax.add_patch(Rectangle((s0, y0 if e.K1 > 0 else y0 - h * 0.9), e.L, h * 0.9, fc='#4C78A8' if e.K1 > 0 else '#E45756', ec='none'))
        elif t == 'Sextupole' and e.L > 0 and e.K2 != 0:
            ax.add_patch(Rectangle((s0, y0 - h * 0.25), e.L, h * 0.5, fc=CE, ec='none'))
        elif t == 'Cavity':
            ax.plot([s0], [y0], marker='v', color='#B22222', ms=5)


def twiss_curve(ffs, i0, i1, nsub=6):
    """要素 i0..i1 (0 始まり) の中を細かく刻んだ s, βx, βy, ηx"""
    S, BX, BY, EX = [], [], [], []
    for i in range(i0, i1 + 1):
        e = ffs.els[i]; r = ffs.rows[i]
        n = nsub if (e.L > 0.05) else 1
        for k in range(n):
            f = k / n
            rr = r if f == 0 else ffs.row_frac(i + 1, f)
            S.append(r.s + f * e.L); BX.append(rr.betx); BY.append(rr.bety); EX.append(rr.etx)
    r = ffs.rows[i1 + 1]; S.append(r.s); BX.append(r.betx); BY.append(r.bety); EX.append(r.etx)
    return np.array(S), np.array(BX), np.array(BY), np.array(EX)


def main():
    data = {}
    for lat, yr, nm in pf_lib.LATS:
        data[lat] = pf_lib.run(lat, 'k64')[0]
    # ---- リング 1 周
    fig, axs = plt.subplots(4, 1, figsize=(10, 11.5), sharex=True)
    for ax, (lat, yr, nm) in zip(axs, pf_lib.LATS):
        ffs = data[lat]; n = len(ffs.els)
        S, BX, BY, EX = twiss_curve(ffs, 0, n - 1, nsub=4)
        ax.plot(S, BX, color=CB, lw=1.1, label='βx'); ax.plot(S, BY, color=CY, lw=1.1, label='βy'); ax.plot(S, EX * 10, color=CE, lw=1.1, label='ηx × 10')
        strip(ax, ffs.els, ffs.rows, -3.0, 2.2)
        ax.set_ylim(-5.5, 32); ax.set_xlim(0, ffs.rows[-1].s)
        ax.set_ylabel('βx, βy [m],  ηx×10 [m]')
        ax.set_title('%s   (C = %.3f m,  νx = %.4f,  νy = %.4f)' % (TITLE[lat], ffs.rows[-1].s, ffs.rows[-1].mux / 2 / math.pi, ffs.rows[-1].muy / 2 / math.pi), fontsize=10, loc='left')
        ax.grid(alpha=0.25)
        if lat == pf_lib.LATS[0][0]:
            ax.legend(ncol=3, loc='upper right', fontsize=9)
    axs[-1].set_xlabel('s [m]')
    fig.suptitle('PF 2.5 GeV の履歴ラティス: リング 1 周の β 関数と分散 (pysad、設計値・誤差なし)\n下の帯: 黒 = 偏向磁石、青 = 収束四極、赤 = 発散四極、緑 = 6 極', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(OUT, 'PF_history_optics.png'), dpi=100); plt.close(fig)
    # ---- 標準セル
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.5))
    for ax, (lat, yr, nm) in zip(axs.flat, pf_lib.LATS):
        ffs = data[lat]
        bi = [i for i, e in enumerate(ffs.els) if T(e) == 'Bend' and getattr(e, 'ANGLE', 0) != 0]
        i0 = bi[5]; i1 = bi[7]                     # 6 番目の偏向磁石の入口から、2 個分 (8.59 m)
        S, BX, BY, EX = twiss_curve(ffs, i0, i1 - 1, nsub=8)
        ax.plot(S, BX, color=CB, lw=1.4, label='βx'); ax.plot(S, BY, color=CY, lw=1.4, label='βy'); ax.plot(S, EX * 10, color=CE, lw=1.4, label='ηx × 10')
        strip(ax, ffs.els, ffs.rows, -1.0, 0.8, ffs.rows[i0].s, ffs.rows[i1].s)
        ax.set_xlim(ffs.rows[i0].s, ffs.rows[i1].s); ax.set_ylim(-2.2, max(BY.max(), BX.max()) * 1.12)
        mx = (ffs.rows[i1].mux - ffs.rows[i0].mux) * 180 / math.pi; my = (ffs.rows[i1].muy - ffs.rows[i0].muy) * 180 / math.pi
        ax.set_title('%s: 偏向磁石 2 個分 (Δs = %.2f m)  Δμx = %.1f°  Δμy = %.1f°' % (TITLE[lat], ffs.rows[i1].s - ffs.rows[i0].s, mx, my), fontsize=9, loc='left')
        ax.set_xlabel('s [m]'); ax.grid(alpha=0.25)
        if lat == pf_lib.LATS[0][0]:
            ax.legend(ncol=3, loc='upper right', fontsize=8)
    fig.suptitle('標準セルの拡大 (偏向磁石 2 個分)。帯の色はリング 1 周の図と同じ', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(OUT, 'PF_history_cell.png'), dpi=100); plt.close(fig)
    # ---- 配置図 (基準軌道を、偏向磁石の ANGLE で積分する。長直線部は細線)
    fig, axs = plt.subplots(2, 2, figsize=(10, 9.4))
    for ax, (lat, yr, nm) in zip(axs.flat, pf_lib.LATS):
        ffs = data[lat]; els = ffs.els; rows = ffs.rows
        x = y = th = 0.0
        mags = [i for i, e in enumerate(els) if is_mag(e)]
        C = rows[-1].s
        longs = []                                 # 長直線部: 磁石のない区間が 4.9 m 以上のもの
        for k, i in enumerate(mags):
            j = mags[(k + 1) % len(mags)]
            s0 = rows[i + 1].s; s1 = rows[j].s + (C if j <= i else 0)
            if s1 - s0 >= 4.9 and abs((s0 + s1) / 2 % C - (0.0 if (s0 + s1) / 2 % C < C / 4 or (s0 + s1) / 2 % C > 3 * C / 4 else C / 2)) < 5:
                longs.append((s0, s1))

        def is_long(s):
            return any((a <= s <= b) or (a <= s + C <= b) for a, b in longs)
        segs = []
        for i, e in enumerate(els):
            L = e.L
            if L <= 0:
                continue
            ang = getattr(e, 'ANGLE', 0) if T(e) == 'Bend' else 0
            n = 8 if ang != 0 else 1
            kind = 'b' if ang != 0 else ('l' if is_long(rows[i].s + L / 2) else 'o')
            for k in range(n):
                dl = L / n; dth = ang / n
                th1 = th + dth * 0.5
                x1 = x + dl * math.cos(th1); y1 = y + dl * math.sin(th1)
                segs.append((x, y, x1, y1, kind)); x, y = x1, y1; th += dth
        for (a, b, c, d, kd) in segs:
            lw, col = (3.0, '#333333') if kd == 'b' else ((0.55, '#555555') if kd == 'l' else (1.8, '#777777'))
            ax.plot([a, c], [b, d], lw=lw, color=col, solid_capstyle='butt')
        ax.set_aspect('equal'); ax.set_title(TITLE[lat] + '   (周回の閉じ方: %.3f m)' % math.hypot(x, y), fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle('リングの配置図 (偏向磁石 3 pt、長直線部 0.55 pt、その他の基準軌道 1.8 pt)', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(OUT, 'PF_history_layout.png'), dpi=100); plt.close(fig)
    for f in ('optics', 'cell', 'layout'):
        print(f, os.path.getsize(os.path.join(OUT, 'PF_history_%s.png' % f)))


if __name__ == '__main__':
    main()
