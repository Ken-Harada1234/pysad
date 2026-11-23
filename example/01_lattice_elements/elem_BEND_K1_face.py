# -*- coding: utf-8 -*-
"""K1≠0 の BEND の端の模型。K64 の tbende・tbendi は、双極磁場を端面から、四極 (K1) の磁場を半径方向の入口・出口の面から始める
('radial': 既定、k64 は常にこれ)。v2 のオプション 'face' は四極の磁場も端面から始める (硬い端の磁石)。
端の角 e ≠ 0 のとき、2 つの模型は半径方向の面と端面の間の楔形の領域 (厚さ x tan e) の四極の蹴りが違う: 半径方向の面で Δpx = +K1 tan(e) x²、Δpy = −K1 tan(e) x y。
独立な基準: ローレンツ力の数値積分 (DOP853、rtol 1e-13) による硬い端の厳密な写像 (pysad/tests/_p2_indep.py の exact_map)。
BEND 1 台 (3 GeV、L = 1.2 m、ANGLE = 0.2618 rad、E1 = E2 = 0.5·ANGLE (矩形))、軌道 (x, px, δ) の (x, px) の出口。
切り替え: 環境変数 PYSAD_BEND_K1_EDGE=face、または pysad.compat.set_bend_k1_edge('face') (v2 だけ)。出力は elem_BEND_K1_face.out.txt。"""
import io, os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'pysad', 'tests'))
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
from pysad import ffscmd, compat, lattice, track6n
import _p2_indep as X

L, TH, E, K1 = 1.2, 0.2617993877991494, 0.5, -0.05
src = f'''MOMENTUM = 3GEV;
BEND BE = (L={L} ANGLE={TH} E1={E} E2={E} K1={K1});
LINE X = (BE);
FFS;
USE X; INS; NORAD; NOFLUC; CALC;
Exit;
'''
compat.set_mode('v2')
ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
el = [x for x in ffs.lat.elements if x.name == 'BE'][0]
print('K1 = %g、端の角 e = %g·ANGLE。出口の (x, px) の厳密な写像との差 (最大)' % (K1, E))
print('%-26s %-10s %-12s %-12s' % ('軌道 (x0, px0, δ)', '模型', 'CALC', '追跡'))
for c in [(1e-3, 0, 0), (3e-3, 0, 0), (0, 5e-4, 0), (2e-3, 3e-4, 0), (1e-3, 0, 5e-3)]:
    for edge in ('radial', 'face'):
        compat.set_bend_k1_edge(edge)
        cc, _m = lattice._bend_element_cod_matrix(el, np.array([c[0], c[1], 0, 0, 0, c[2]], float))
        t = np.asarray(track6n.track_particles_fast([el], np.array([[c[0], c[1], 0, 0, 0, c[2]]], float), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * L)[0])
        cl = [c[0], c[1], 0, 0, 0, c[2]]
        ex = X.exact_map(cl, el.K1, E * TH, E * TH, quad=edge)
        exo = X.exact_map(cl, el.K1, E * TH, E * TH, quad='face' if edge == 'radial' else 'radial')
        print('%-26s %-10s %-12.1e %-12.1e (別の模型の厳密な写像との差 %.1e)' % (str(c), edge, np.abs(cc[:2] - ex).max(), np.abs(t[:2] - ex).max(), np.abs(ex - exo).max()))
compat.set_bend_k1_edge(None)
