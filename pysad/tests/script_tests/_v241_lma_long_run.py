# -*- coding: utf-8 -*-
"""v2.4.1: 10^4 周級の長時間追跡 (放射・量子ゆらぎ FLUC あり/なし) の gs・pysad 比較。
閉軌道 (RADCOD) から出発する N 粒子 (fluc: 振幅 0 から FLUC で広がる。nofluc: 初期振幅を与えて放射減衰だけを見る) を、位置 1 から interval 周ずつ TrackParticles で追跡して
座標を印刷する。台本は gs・pysad の台本評価器に同じものを渡す。解析は _v241_lma_long_analyze.py (座標 → σδ・εx・εy・減衰、Emittance[] の値との比較)。
使い方: py -3.12 pysad/tests/script_tests/_v241_lma_long_run.py RING MODE N NTURN INTERVAL ID [--nogs] [--nopy]
  RING = astrid2 | kekls (kek-ls_v03_12.lat、3 GeV、周長 603.12 m、全周)   MODE = fluc | nofluc
出力: pysad/tests/script_tests/scripts/v241_lma_long_<ring>_<mode>_n<N>_t<NTURN>.{sad,pysad.txt,gs.txt}
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _v241_lma_tools as T

PATH = {'astrid2': T.ASTRID2, 'kekls': T.KEKLS12}
# nofluc の初期振幅 (x, y [m], δ): 小振幅 (線形に近い)
A0 = {'astrid2': (1e-3, 0.5e-3, 3e-3), 'kekls': (0.5e-3, 0.3e-3, 3e-3)}


def make_script(ring, mode, n, nturn, interval, seed=1):
    src = T.read(PATH[ring])
    rng = np.random.default_rng(seed)
    X = np.zeros((6, n))
    if mode == 'nofluc':
        ax, ay, ad = A0[ring]
        # 位相をばらけさせた (x, y, δ) の初期の偏り: 振幅 ax, ay, ad の 0.2〜1 倍
        f = rng.uniform(0.2, 1.0, (3, n))
        X[0], X[2], X[5] = ax * f[0] * rng.choice([-1, 1], n), ay * f[1] * rng.choice([-1, 1], n), ad * f[2] * rng.choice([-1, 1], n)
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    st, _ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.1)
    cod = st.closed_orbit([1])[0]
    cols = ['{' + ','.join('%.17g' % v for v in X[k]) + '}' for k in range(6)]
    cols.append('{' + ','.join('1' for _ in range(n)) + '}')
    fl = 'FLUC' if mode == 'fluc' else 'NOFLUC'
    nc = nturn // interval
    L = ['FFS;', 'USE RING; CELL; CALC; RFSW; RAD; RADCOD; %s;' % fl, 'emit = Emittance[];', 'Print["EMIT ", emit];',
         'LOSSAMPL = 0.1; $FORM = "20.14"; NPARA = 1;',
         'COD0 = {%s};' % ','.join('%.17g' % v for v in cod),
         'beam = {%s};' % ','.join(cols),
         'b = {beam[[1]] + COD0[[1]], beam[[2]] + COD0[[2]], beam[[3]] + COD0[[3]], beam[[4]] + COD0[[4]], beam[[5]] + COD0[[5]], beam[[6]] + COD0[[6]], beam[[7]]};',
         'Do[r = TrackParticles[{1, b}, 1, 1, %d]; b = r[[2]];' % interval,
         '  Do[Print["S ", k*%d, " ", i, " ", b[[1,i]], " ", b[[2,i]], " ", b[[3,i]], " ", b[[4,i]], " ", b[[5,i]], " ", b[[6,i]], " ", b[[7,i]]], {i, 1, %d}], {k, 1, %d}];' % (interval, n, nc),
         'Exit;']
    return src + '\n' + '\n'.join(L) + '\n'


def main():
    a = sys.argv[1:]
    ring, mode, n, nturn, interval, uid = a[0], a[1], int(a[2]), int(a[3]), int(a[4]), a[5]
    tag = 'v241_lma_long_%s_%s_n%d_t%d' % (ring, mode, n, nturn)
    script = make_script(ring, mode, n, nturn, interval)
    sad = os.path.join(T.SP, tag + '.sad')
    open(sad, 'w', encoding='utf-8').write(script)
    if '--nopy' not in a:
        _t, tpy = T.run_pysad_text(script, os.path.join(T.SP, tag + '.pysad.txt'))
        print('pysad: %.0f s' % tpy, flush=True)
    if '--nogs' not in a:
        gsname = os.path.join(T.SP, '_gs_lma_%s_%s.sad' % (uid, tag))
        open(gsname, 'w', encoding='utf-8').write(script)
        _g, tgs, left = T.run_gs(gsname, os.path.join(T.SP, tag + '.gs.txt'))
        if os.path.exists(gsname):
            os.remove(gsname)
        print('gs: %.0f s、自分の台本の sad1.exe の残り: %s' % (tgs, left or 'なし'), flush=True)


if __name__ == '__main__':
    main()
