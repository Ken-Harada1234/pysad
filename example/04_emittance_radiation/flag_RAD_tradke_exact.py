# -*- coding: utf-8 -*-
"""放射反作用の角度の項 ddpx の符号 (B111)。K64 の tradke (Emittance[]・閉軌道) は ddpx = +0.5 dpx dg、tradkn (追跡) は -0.5 dpx dg (dg < 0)。
独立な基準: 放射反作用つきのローレンツ力の数値積分 (pysad/tests/script_tests/_p2_indep.py の bend_with_radiation: 一様磁場の扇形の偏向磁石。連続な放射)。
扇形の偏向磁石 1 台 (3 GeV、L = 1.2 m、ANGLE = 0.2618 rad、入口の軌道 x = 2 mm・px = 1 mm・δ = -1e-4) の、放射 ON と OFF の出口の差 Δx・Δpx。
放射の強さ k は、pysad の v2 の Δδ に合わせる。結果 (出力は flag_RAD_tradke_exact.out.txt): v2 の tradke は厳密な積分と 5 桁一致 (Δpx で 6e-6)、k64 の tradke は Δpx が 2.2e-4 ずれる
(符号が逆の ddpx の分。偏向角が小さい 0.05 rad の磁石では 1.6e-3)。追跡 (TrackParticles、tradkn) は k64・v2 とも厳密な積分と 5 桁一致。"""
import io, math, os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'pysad', 'tests', 'script_tests'))
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
from pysad import ffscmd, compat, emit6d
from pysad.elements import Bend
import _p2_indep as X

TH, L, X0, PX0, G0 = 0.2618, 1.2, 2e-3, 1e-3, -1e-4
DIV = int(sys.argv[1]) if len(sys.argv) > 1 else 20       # EMITDIV (放射の刻みの倍率)。既定の 1 では K64 の tradke の誤差が 9%


def emit_end(mode, radcod):
    compat.set_mode(mode)
    el = Bend(name='B', L=L, ANGLE=TH, E1=0.0, E2=0.0)
    emit6d.EMITDIV[:] = [DIV, DIV, DIV]          # _Phys が作られるとき EMITDIV を読むので、先に設定する
    ph = emit6d._Phys(3e9, 100.0, amass=emit6d.ELMASS, charge=1.0)
    st = emit6d.tturne([el], np.array([X0, PX0, 0, 0, 0, G0], float), ph, trf0=0.0, rfsw=False, radcod=radcod, want_db=False, want_jac=False)
    return st.cod


def track_end(mode, rad):
    compat.set_mode(mode)
    src = f'''MOMENTUM = 3GEV;
BEND B = (L={L} ANGLE={TH});
LINE RING = (B);
FFS;
USE RING;
{'RAD; NOFLUC;' if rad else 'NORAD;'}
r = TrackParticles[{{1, {{{{{X0}}},{{{PX0}}},{{0}},{{0}},{{0}},{{{G0}}},{{1}}}}}}, 1, 1, 1];
Print["R ", Flatten[Take[r[[2]],6]]];
Exit;
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    l = [x for x in buf.getvalue().splitlines() if x.startswith('R ')][0]
    return np.array([float(v) for v in l[2:].replace('{', ' ').replace('}', ' ').replace(',', ' ').split()])


off = emit_end('v2', False)
dg = emit_end('v2', True)[5] - G0
k = X.radiation_reaction_k_for_dg(X0, PX0, G0, TH, L, dg)
xc0, pxc0, _, _ = X.bend_with_radiation(X0, PX0, G0, TH, L, 0.0)
xe, pxe, _, _ = X.bend_with_radiation(X0, PX0, G0, TH, L, k)
print('EMITDIV = %d。放射による Δδ = %.4e (放射の強さ k を合わせた)。放射なしの出口: pysad x %.6e px %.6e、厳密 x %.6e px %.6e' % (DIV, dg, off[0], off[1], xc0, pxc0))
print('厳密な積分 (連続な放射の反作用):      Δx %.5e  Δpx %.5e' % (xe - xc0, pxe - pxc0))
for mode in ('k64', 'v2'):
    on = emit_end(mode, True)
    print('Emittance[] の tradke (%s):   Δx %.5e  Δpx %.5e   (厳密との比 Δpx: %.6f)' % (mode, on[0] - off[0], on[1] - off[1], (on[1] - off[1]) / (pxe - pxc0)))
t0, t1 = track_end('v2', False), track_end('v2', True)
print('追跡の tradkn (v2):  Δx %.5e  Δpx %.5e   (厳密との比 Δpx: %.6f)' % (t1[0] - t0[0], t1[1] - t0[1], (t1[1] - t0[1]) / (pxe - pxc0)))
t0, t1 = track_end('k64', False), track_end('k64', True)
print('追跡の tradkn (k64): Δx %.5e  Δpx %.5e   (厳密との比 Δpx: %.6f)' % (t1[0] - t0[0], t1[1] - t0[1], (t1[1] - t0[1]) / (pxe - pxc0)))
