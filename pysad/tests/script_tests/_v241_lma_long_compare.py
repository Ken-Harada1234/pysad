# -*- coding: utf-8 -*-
"""v2.4.1: 長時間追跡 (scripts/v241_lma_long_*) の gs・pysad (k64・v2) の比較表。有る組だけ出す。
使い方: py -3.12 pysad/tests/script_tests/_v241_lma_long_compare.py > pysad/tests/script_tests/scripts/v241_lma_long_summary.out.txt"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _v241_lma_long_analyze as A
SP = os.path.join(HERE, 'scripts')
RINGS = {'kekls': 'KEK-LS (3 GeV、603.12 m、3563 要素)', 'astrid2': 'ASTRID2 (0.58 GeV、45.704 m、199 要素)'}


def J(d, tw, cod):
    X = d['X'] - cod[None, None, :]
    x, px, y, py, z, dp = [X[:, :, j] for j in range(6)]
    xb, pb = x - tw['etx'] * dp, px - tw['etpx'] * dp
    Jx = (1 + tw['alfx'] ** 2) / tw['betx'] * xb ** 2 + 2 * tw['alfx'] * xb * pb + tw['betx'] * pb ** 2
    Jy = (1 + tw['alfy'] ** 2) / tw['bety'] * y ** 2 + 2 * tw['alfy'] * y * py + tw['bety'] * py ** 2
    return Jx, Jy


def f(ring, mode, n, t=10000):
    return os.path.join(SP, 'v241_lma_long_%s_%s_n%d_t%d' % (ring, mode, n, t))


def main():
    for ring in RINGS:
        tw, cod = A.twiss(ring)
        print('== %s' % RINGS[ring])
        # NOFLUC
        b = f(ring, 'nofluc', 10)
        if os.path.exists(b + '.gs.txt'):
            g = A.load(b + '.gs.txt')
            em = g['emit']
            print('NOFLUC 10 粒子 10000 周: Emittance[] の減衰率 α (振幅/周) x,y,z = %s' % ['%.4e' % v for v in em['alpha']])
            for nm, fn in (('pysad k64', '.pysadk64.txt'), ('pysad v2', '.pysad.txt')):
                if os.path.exists(b + fn):
                    o = A.load(b + fn)
                    d = np.abs(g['X'] - o['X'])
                    print('  gs − %s: 座標の最大の差 n=250 %.1e、n=10000 %.1e' % (nm, d[0].max(), d[-1].max()))
            for nm, fn in (('gs', '.gs.txt'), ('pysad k64', '.pysadk64.txt'), ('pysad v2', '.pysad.txt')):
                if os.path.exists(b + fn):
                    d = A.load(b + fn)
                    Jx, Jy = J(d, tw, cod)
                    n = d['n']
                    r = []
                    for Jq in (Jx, Jy):
                        s = np.array([-np.polyfit(n, np.log(Jq[:, i]), 1)[0] / 2 for i in range(Jq.shape[1])])
                        r.append((s.mean(), s.std(ddof=1) / np.sqrt(len(s))))
                    print('  %-10s 追跡の α_x = %.3e ± %.1e (Emittance[] の %.3f 倍)、α_y = %.3e ± %.1e (%.3f 倍)' % (
                        nm, r[0][0], r[0][1], r[0][0] / em['alpha'][0], r[1][0], r[1][1], r[1][0] / em['alpha'][1]))
        # FLUC
        for n in (40, 120):
            b = f(ring, 'fluc', n)
            if not (os.path.exists(b + '.pysad.txt') or os.path.exists(b + '.pysadk64.txt')):
                continue
            ref = None
            for fn in ('.gs.txt', '.pysadk64.txt', '.pysad.txt'):
                if os.path.exists(b + fn):
                    ref = A.load(b + fn)
                    break
            em = ref['emit']
            ax, ay, az = em['alpha']
            print('FLUC %d 粒子 10000 周: Emittance[] 平衡 εx %.3e σδ %.3e σz %.3e、α_x %.3e α_z %.3e' % (n, em['eps'][0], em['sigd'], em['sigz'], ax, az))
            g = A.load(b + '.gs.txt') if os.path.exists(b + '.gs.txt') else None
            for nm, fn in (('pysad k64', '.pysadk64.txt'), ('pysad v2', '.pysad.txt')):
                if g is not None and os.path.exists(b + fn):
                    o = A.load(b + fn)
                    d = np.abs(g['X'] - o['X'])
                    print('  gs − %s: 座標の最大の差 n=250 %.1e、n=10000 %.1e (x, px, y, py, z, δ の最大)' % (nm, d[0].max(), d[-1].max()))
            for nm, fn in (('gs', '.gs.txt'), ('pysad k64', '.pysadk64.txt'), ('pysad v2', '.pysad.txt')):
                if os.path.exists(b + fn):
                    d = A.load(b + fn)
                    s = A.stats(d, tw, cod)
                    t = s['n']
                    fa = lambda a: 1 - np.exp(-2 * a * t)
                    m = t >= 1500
                    r_sd = s['sd'] ** 2 / (em['sigd'] ** 2 * fa(az))
                    r_ex = s['ex'] / (em['eps'][0] * fa(ax))
                    # 標準誤差: 時間相関 (同期振動の周期 ~250 周) を無視しない見積もり = 独立な点を 1500 周おきの 6 点とみなす
                    ind = np.array([np.where(t >= v)[0][0] for v in (1500, 3000, 4500, 6000, 7500, 9000)])
                    print('  %-10s 生き残り %d、t≥1500 の平均 σδ²/予想 = %.3f (6 点の標準誤差 %.3f)、εx/予想 = %.3f (%.3f)、σz²/予想 %.3f、εy 最大 %.1e' % (
                        nm, s['nalive'][-1], r_sd[m].mean(), r_sd[ind].std(ddof=1) / np.sqrt(6), r_ex[m].mean(), r_ex[ind].std(ddof=1) / np.sqrt(6),
                        (s['sz'] ** 2 / (em['sigz'] ** 2 * fa(az)))[m].mean(), s['ey'].max()))


if __name__ == '__main__':
    main()
