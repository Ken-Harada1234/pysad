# -*- coding: utf-8 -*-
"""v2.3.0 Q1: 周波数マップの例 (pysad/freqmap.py)。1 つの位置 (長直線部) から (x, y) の格子を追跡し、NAFF で前半・後半のチューンを出して
D = log10 sqrt(Δνx² + Δνy²) を図にする (チューン平面と初期条件平面)。拡散の大きい領域が、チューン図の共鳴線 (次数 4 以下) の上に載ることを確かめる。
使い方 (リポジトリ直下から): py -3.12 example/v2_3/Q1_fmap_ring.py astrid2   または   kekls   [workers]
条件 (どちらも v2 モード、RFSW OFF、RAD OFF、δ = 0、px = py = 0・z = 0 から出発): 追跡 NTURN = 1024 周 (前半・後半 512 周ずつ)、NAFF の窓 = Hann (次数 1)、
アパーチャ LOSSAMPL = 0.1 m (実質なし。失われる点は力学的な不安定)、格子 41 × 31 (x が速く変わる)。
出力: Q1_fmap_<名前>.table.txt (表)、.png (図)、標準出力 (.out.txt)。
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), loc=2, x=(0.0005, 0.018), y=(0.0005, 0.014),
                    label='ASTRID2 (0.58 GeV、周長 45.7 m、6 セル)、位置 2 (長直線部の入口、βx=3.9 m、βy=1.6 m)'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), loc=3, x=(0.0002, 0.006), y=(0.0002, 0.008),
                  label='kek-ls_v03_12.lat (3 GeV、周長 603 m、20 セル)、位置 3 (長直線部の入口、βx=4.6 m、βy=6.2 m)'),
}
NX, NY, NTURN = 41, 31, 1024


def main(name, workers=1):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.freqmap import frequency_map, resonance_check, line_distance, resonance_lines
    compat.set_mode('v2')
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)
    print(cs['label'])
    nux, nuy = st.tunes_full()
    print('作動点 (CALC): nu_x = %.6f, nu_y = %.6f' % (nux, nuy))
    print('作動点に近い共鳴線 (次数 4 以下。距離はチューン平面の距離。整数・半整数は必ず入れた):')
    for (m, n, p, d, o, nm) in resonance_check(nux, nuy, 4):
        print('   %2d nu_x %+2d nu_y = %4d   距離 %.4f   %s' % (m, n, p, d, nm))
    t0 = time.time()
    r = frequency_map(st, cs['loc'], x=np.linspace(*cs['x'], NX), y=np.linspace(*cs['y'], NY), nturn=NTURN, workers=workers, progress=1)
    print('追跡と NAFF: %.0f s' % (time.time() - t0))
    print('条件:', r.cond)
    r.save(os.path.join(HERE, 'Q1_fmap_%s.table.txt' % name))
    r.plot(os.path.join(HERE, 'Q1_fmap_%s.png' % name), order=4)
    ok = np.isfinite(r.D)
    print('点 %d、失われた点 %d、D の範囲 %.2f 〜 %.2f' % (len(r), (r.lost_turn > 0).sum(), np.nanmin(r.D), np.nanmax(r.D)))
    # 共鳴線との一致: 点 (nu_x1, nu_y1) から、次数 order 以下の共鳴線までの最短距離 (チューン平面の距離) で、近い (< EPS) 点と遠い点を比べる
    EPS = 0.003
    for order in (3, 4, 5):
        dist = line_distance(r.nux1[ok], r.nuy1[ok], order)
        D = r.D[ok]
        near = dist < EPS
        print('次数 %d 以下の共鳴線から %.3f 以内 (近い) %d 点、それより遠い %d 点: D の中央値 %.2f / %.2f、D > -5 の割合 %.3f / %.3f、D > -4 の点の %.2f が近い' % (
            order, EPS, near.sum(), (~near).sum(), np.median(D[near]) if near.any() else np.nan, np.median(D[~near]) if (~near).any() else np.nan,
            np.mean(D[near] > -5) if near.any() else np.nan, np.mean(D[~near] > -5) if (~near).any() else np.nan,
            np.mean(near[D > -4]) if (D > -4).any() else np.nan))
    # D が大きい点の最も近い共鳴線
    from collections import Counter
    cnt = Counter()
    hiidx = np.nonzero(ok & (r.D > -4.0))[0]
    for i in hiidx:
        best = None
        for (m, n, p) in resonance_lines(5, ((r.nux1[i] - 0.05, r.nux1[i] + 0.05), (r.nuy1[i] - 0.05, r.nuy1[i] + 0.05))):
            d = abs(m * r.nux1[i] + n * r.nuy1[i] - p) / np.hypot(m, n)
            if best is None or d < best[0]:
                best = (d, m, n, p)
        if best:
            cnt[(best[1], best[2], best[3])] += 1
    print('D > -4 の点 (%d 点) の最も近い共鳴線 (次数 5 以下) の上位:' % len(hiidx))
    for k, v in cnt.most_common(8):
        print('   %d nu_x %+d nu_y = %d (整数部を除いた端数の平面で)  %d 点' % (k[0], k[1], k[2], v))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 1)
