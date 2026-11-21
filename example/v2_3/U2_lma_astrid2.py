# -*- coding: utf-8 -*-
"""v2.3.0 U2: ASTRID2 の局所運動量アクセプタンス (1 セル、pysad の LMA) を K1 の定義 'log' (K64) と 'lin' (B' 一様) で比べる。
使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/v2_3/U2_lma_astrid2.py [workers]   (標準出力を U2_lma_astrid2.out.txt に保存した)
条件は Q1_lma_ring.py と同じ (RF・放射つき、アパーチャ ±10 mm、check_every = 1、既定の追跡周回数)。"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)


def main(workers=4):
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity, select_locations
    from pysad.lma import momentum_acceptance
    compat.set_mode('v2')
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    res = {}
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.010)
        st = st.copy(check_every=1)
        per = find_periodicity(st.els, st.rows)
        locs = select_locations(st, per['starts'][0], per['starts'][1] - 1, 1, 0.0)
        t0 = time.time()
        r = momentum_acceptance(st, locs, workers=workers)
        res[d] = r
        print('%s: %d 位置、%.0f s' % (d, len(locs), time.time() - t0))
    a, b = res['log'], res['lin']
    dp = np.abs(np.array(b.dp_plus) / np.array(a.dp_plus) - 1)
    dm = np.abs(np.array(b.dp_minus) / np.array(a.dp_minus) - 1)
    print('δ+ : log 最小 %.5f 最大 %.5f / lin 最小 %.5f 最大 %.5f、相対差の最大 %.2e' % (min(a.dp_plus), max(a.dp_plus), min(b.dp_plus), max(b.dp_plus), dp.max()))
    print('δ- : log 最小 %.5f 最大 %.5f / lin 最小 %.5f 最大 %.5f、相対差の最大 %.2e' % (min(a.dp_minus), max(a.dp_minus), min(b.dp_minus), max(b.dp_minus), dm.max()))
    print('LMA の刻み (dpmax/…) 程度の差か: 差の最大 δ+ %.2e δ- %.2e (絶対)' % (np.max(np.abs(np.array(b.dp_plus) - np.array(a.dp_plus))), np.max(np.abs(np.array(b.dp_minus) - np.array(a.dp_minus)))))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
