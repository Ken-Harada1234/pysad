# -*- coding: utf-8 -*-
"""example/kek_ls の図 (PNG) を、gs の出力 (NAME.gs.txt) の数値から作る。使い方 (リポジトリ直下から): py -3.12 example/tools/kek_ls_figs.py
作図のルール (D:\\work\\CLAUDE.md): ラティス配置図は、長直線部を細線 (0.55 pt)、その他の基準軌道 1.8 pt、偏向磁石 3 pt。"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = ['Yu Gothic', 'Meiryo', 'MS Gothic', 'sans-serif']
matplotlib.rcParams['axes.unicode_minus'] = False
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'kek_ls')
CB, CY, CE = '#1f5fa8', '#c0392b', '#2a8a3e'   # βx、βy、ηx


def lines(rel, kind='gs'):
    return open(os.path.join(EX, rel + '.%s.txt' % kind), encoding='utf-8', errors='replace').read().split('\n')


def save(fig, rel):
    fig.tight_layout()
    fig.savefig(os.path.join(EX, rel), dpi=90)
    plt.close(fig)
    print(rel, os.path.getsize(os.path.join(EX, rel)))


def rows(rel, tag, kind='gs'):
    out = []
    for l in lines(rel, kind):
        if l.startswith(tag + ' '):
            t = l.split()
            out.append(t[1:])
    return out


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float('nan')


# ---- 1. 標準セルの配置とオプティクス
def cell_figure():
    R = [r for r in rows('lattice/KEKLS_cell_optics', 'E') if len(r) >= 10]
    name = [r[0] for r in R]; typ = [r[1] for r in R]
    S = np.array([fnum(r[2]) for r in R]); L = np.array([fnum(r[3]) for r in R]); ANG = np.array([fnum(r[4]) for r in R])
    K1 = np.array([fnum(r[5]) for r in R]); K2 = np.array([fnum(r[6]) for r in R])
    BX = np.array([fnum(r[7]) for r in R]); BY = np.array([fnum(r[8]) for r in R]); EXX = np.array([fnum(r[9]) for r in R])
    C = S[-1] + L[-1]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.2), gridspec_kw={'height_ratios': [1, 1.5]}, sharex=True)
    # 配置: 基準軌道を s 方向に引く。長直線部 (DRID1) は 0.55 pt、その他のドリフト・要素 1.8 pt、偏向磁石 3 pt
    for i in range(len(R)):
        s0, s1 = S[i], S[i] + L[i]
        if typ[i] == 'DRIFT':
            lw = 0.55 if name[i].startswith('DRID1') else 1.8
            a1.plot([s0, s1], [0, 0], color='k', lw=lw, solid_capstyle='butt')
        elif typ[i] == 'BEND':
            if abs(ANG[i]) > 0:
                a1.plot([s0, s1], [0, 0], color='#7f3fbf', lw=3, solid_capstyle='butt')
                if K1[i] != 0:
                    a1.add_patch(plt.Rectangle((s0, -0.05), L[i], 0.1, fc='none', ec='#7f3fbf', lw=0.6))
        else:
            a1.plot([s0, s1], [0, 0], color='k', lw=1.8, solid_capstyle='butt')
        if typ[i] == 'QUAD':
            sgn = 1 if K1[i] > 0 else -1
            a1.add_patch(plt.Rectangle((s0, 0 if sgn > 0 else -0.45), max(L[i], 0.05), 0.45, fc='#e67e22' if sgn > 0 else '#2c7fb8', ec='none'))
        elif typ[i] == 'SEXT':
            a1.add_patch(plt.Rectangle((s0, 0.0), max(L[i], 0.05), 0.7, fc='#2a8a3e', ec='none'))
        elif typ[i] == 'MONI':
            a1.plot([s0], [-0.65], marker='v', ms=3, color='gray')
    a1.set_ylim(-0.9, 1.0); a1.set_yticks([])
    a1.set_title('KEK-LS v03_68 標準セル (ARCA, %.3f m): 配置 (橙: 収束四極、青: 発散四極、緑: 六極、紫太線: 偏向磁石、細線: 長直線部)' % C, fontsize=9)
    a2.plot(S, BX, color=CB, lw=1.3, label='βx'); a2.plot(S, BY, color=CY, lw=1.3, label='βy'); a2.plot(S, EXX * 20, color=CE, lw=1.3, label='ηx × 20')
    a2.set_xlabel('s [m]'); a2.set_ylabel('βx, βy [m],  ηx×20 [m]'); a2.legend(ncol=3, loc='upper center'); a2.set_xlim(0, C)
    a2.set_title('標準セルのオプティクス (gs の出力)', fontsize=9)
    save(fig, 'lattice/KEKLS_cell_layout_optics.png')


# ---- 2. リングの配置 (x-z) とオプティクス
def ring_figure():
    R = rows('display/KEKLS_optics_table', 'E')
    name = [r[0] for r in R]
    S = np.array([fnum(r[1]) for r in R]); BX = np.array([fnum(r[2]) for r in R]); BY = np.array([fnum(r[3]) for r in R]); EXX = np.array([fnum(r[4]) for r in R])
    NX = np.array([fnum(r[5]) for r in R]); NY = np.array([fnum(r[6]) for r in R]); GX = np.array([fnum(r[7]) for r in R]); GZ = np.array([fnum(r[8]) for r in R])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={'width_ratios': [1, 1.6]})
    for i in range(len(R) - 1):
        lw = 0.55 if name[i].startswith('DRID1') else 1.8
        a1.plot([GX[i], GX[i + 1]], [GZ[i], GZ[i + 1]], color='k', lw=lw, solid_capstyle='butt')
    a1.set_aspect('equal'); a1.set_xlabel('X [m]'); a1.set_ylabel('Y [m]')
    a1.set_title('RING の配置 (20 セル。長直線部は細線。入射セルは S = 28.5 m から)', fontsize=9)
    a2.plot(S, BX, color=CB, lw=0.9, label='βx'); a2.plot(S, BY, color=CY, lw=0.9, label='βy'); a2.plot(S, EXX * 50, color=CE, lw=0.9, label='ηx × 50')
    a2.set_xlabel('s [m]'); a2.set_ylabel('[m]'); a2.legend(ncol=3, loc='upper right'); a2.set_xlim(0, S.max())
    a2.set_title('RING のオプティクス (gs の出力)', fontsize=9)
    save(fig, 'display/KEKLS_ring_layout_optics.png')


# ---- 3. 運動量依存のチューン
def chrom_figure():
    r1 = np.array([[fnum(x) for x in r] for r in rows('display/KEKLS_chromaticity', 'T1')])
    r0 = np.array([[fnum(x) for x in r] for r in rows('display/KEKLS_chromaticity', 'T0')])
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, a, t in ((axs[0], r1, '設計の六極あり'), (axs[1], r0, '六極なし (K2 = 0)')):
        ax.plot(a[:, 0] * 100, a[:, 1] - a[0, 1] * 0 - 48.58, 'o-', color=CB, label='νx − 48.58')
        ax.plot(a[:, 0] * 100, a[:, 2] - 17.62, 's-', color=CY, label='νy − 17.62')
        ax.set_xlabel('δ [%]'); ax.set_ylabel('Δν'); ax.set_title('RINGC のチューン対運動量: ' + t, fontsize=9); ax.grid(alpha=.3); ax.legend()
    save(fig, 'display/KEKLS_chromaticity.png')


# ---- 4. 振幅依存チューンシフト
def adts_figure():
    ax_ = np.array([[fnum(x) for x in r[:3]] for r in rows('display/KEKLS_amp_dep_ts', 'AX')])
    ay_ = np.array([[fnum(x) for x in r[:3]] for r in rows('display/KEKLS_amp_dep_ts', 'AY')])
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, a, t in ((axs[0], ax_, '水平振幅 (n σx)'), (axs[1], ay_, '垂直振幅 (n σy)')):
        ok = np.abs(a[:, 1]) < 1e3
        ax.plot(a[ok, 0], a[ok, 1] - 48.58, 'o-', color=CB, label='Δνx'); ax.plot(a[ok, 0], a[ok, 2] - 17.62, 's-', color=CY, label='Δνy')
        ax.set_xlabel('n [σ]'); ax.set_ylabel('Δν'); ax.set_title('振幅依存チューンシフト: ' + t, fontsize=9); ax.grid(alpha=.3); ax.legend()
    save(fig, 'display/KEKLS_amp_dep_ts.png')


# ---- 5. 六極のスキャン (dtux・dtuy の地図)
def sext_figure():
    a = np.array([[fnum(x) for x in r] for r in rows('display/KEKLS_sextupole_tune_scan', 'S')])
    ds = np.unique(a[:, 0]); n = len(ds)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for k, (ax, t) in enumerate(zip(axs, ('Σ|Δνx|', 'Σ|Δνy|'))):
        z = a[:, 2 + k].reshape(n, n)  # 行: dsd、列: dsf
        im = ax.imshow(z.T, origin='lower', extent=[ds[0] - 0.75, ds[-1] + 0.75, ds[0] - 0.75, ds[-1] + 0.75], aspect='auto', cmap='viridis')
        for i in range(n):
            for j in range(n):
                ax.text(ds[i], ds[j], '%.2f' % z[i, j], ha='center', va='center', fontsize=6, color='w')
        ax.set_xlabel('ΔK2(SDH)'); ax.set_ylabel('ΔK2(SFH)'); ax.set_title(t + ' (δ = ±3% の和)', fontsize=9)
        fig.colorbar(im, ax=ax)
    save(fig, 'display/KEKLS_sextupole_tune_scan.png')


# ---- 6. セルのチューンサーベイ
def survey_figure(kind='gs'):
    P = np.array([[fnum(x) for x in r] for r in rows('design/KEKLS_tune_survey_cell', 'POINT', kind)])
    tx = np.unique(P[:, 0]); ty = np.unique(P[:, 1])
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
    for ax, col, t, cm in ((axs[0], 6, 'εx [pm rad]', 'viridis'), (axs[1], 4, 'ξx', 'coolwarm'), (axs[2], 5, 'ξy', 'coolwarm')):
        z = P[:, col].reshape(len(tx), len(ty))
        if col == 6:
            z = np.clip(z, 0, 450)
        else:
            z = np.clip(z, -1.5, 1.5)
        im = ax.pcolormesh(np.append(tx - 0.075, tx[-1] + 0.075), np.append(ty - 0.075, ty[-1] + 0.075), z.T, cmap=cm, shading='flat')
        fig.colorbar(im, ax=ax); ax.set_xlabel('νx (リング)'); ax.set_ylabel('νy (リング)'); ax.set_title(t + ' (%s)' % kind, fontsize=9)
        for xv in (48.5, 49.0): ax.axvline(xv, color='k', lw=0.4, ls=':')
        for yv in (17.5, 18.0): ax.axhline(yv, color='k', lw=0.4, ls=':')
    save(fig, 'design/KEKLS_tune_survey_cell_%s.png' % kind)


if __name__ == '__main__':
    todo = sys.argv[1:] or ['cell', 'ring', 'chrom', 'adts', 'sext', 'survey']
    for t in todo:
        try:
            {'cell': cell_figure, 'ring': ring_figure, 'chrom': chrom_figure, 'adts': adts_figure, 'sext': sext_figure,
             'survey': lambda: (survey_figure('gs'), survey_figure('pysad'))}[t]()
        except Exception as e:
            print(t, 'failed:', repr(e))
