# -*- coding: utf-8 -*-
"""v2.4.0 F1: ROTATE した偏向磁石を RDT・ADTS (resdrive) と ν(J, δ) の級数 (ampnf) に入れた検証。
使い方: py -3.12 example/v2_4/F1_rotated_bend.py fodo|astrid  (出力は標準出力。保存は F1_*.out.txt)
条件: k64 (模型は K1 = 0 の BEND か ASTRID2 の複合機能 BEND。pysad の追跡は k64 = gs)、RF・放射なし、結合した光学の固有モード (4×4)、1024 周の NAFF。
fodo: 3 GeV の 24 セル FODO、6 極 K2 = 4・−6。全 24 台を θ = 0.05・0.1・0.2 傾ける。
astrid: ASTRID2 (0.58 GeV、45.7 m、複合機能 BEND)。全 BEND の ROTATE を θ = 1・5・20 mrad にする (据え付けのロールの典型 〜 大きめ)。"""
import math
import os
import sys
import time
import dataclasses
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import ampnf as A, compat, resdrive as R
from pysad.lattice import Lattice
from pysad.parser import load_sad_file, load_sad_string
compat.set_mode('k64')
P = lambda *a: print(*a, flush=True)


def roll_all_bends(lat, th):
    out = []
    for el in lat.elements:
        if type(el).__name__ == 'Bend' and float(el.ANGLE) != 0.0 and float(el.L) > 0:
            el = dataclasses.replace(el, ROTATE=th)
        out.append(el)
    return Lattice(out, name=lat.name)


def fodo():
    src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
        'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
    return load_sad_string(src).build_lattice('RING'), 3e9, 3e-7, (0.2, 0.1, 0.05)


def astrid():
    lat = load_sad_file(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')).build_lattice('RING')
    return lat, 0.58e9, 3e-9, (20e-3, 5e-3, 1e-3)


def main(name):
    lat0, p0, j, thetas = {'fodo': fodo, 'astrid': astrid}[name]()
    opt = dict(D=5, nslice=16, bend_cubic=True, kin=True, fringe=True, nslice_bend=4)
    P('# %s: 全 BEND の ROTATE = θ。adts().total = (ロールなしの解析式) + .rotate (結合した光学の模型の J¹ の差)。' % name)
    P('# ampnf = 模型 D=5 の J¹ (全部入り)、追跡 = pysad の追跡 (全素子) + NAFF のフィット (正規形の J・φ から出発、J = %g の格子 8 点)。単位 [1/(m rad)]' % j)
    P('θ [rad], 結合, 項, xx, xy, yy')
    base = R.adts(lat0)
    r0 = A.amplitude_tune(lat0, **opt)
    tr0 = A.tracked_adts(lat0, r0, j=j, momentum_eV=p0)
    P('0 | 0 | 解析 (ロールなし) | %.4f %.4f %.4f' % tuple(base.total[k] for k in ('xx', 'xy', 'yy')))
    P('0 | 0 | ampnf | %.4f %.4f %.4f' % tuple(r0.adts[k] for k in ('xx', 'xy', 'yy')))
    P('0 | 0 | 追跡 | %.4f %.4f %.4f' % tuple(tr0[k] for k in ('xx', 'xy', 'yy')))
    for th in thetas:
        t0 = time.time()
        lat = roll_all_bends(lat0, th)
        a = R.adts(lat)
        r = A.amplitude_tune(lat, **opt)
        tr = A.tracked_adts(lat, r, j=j, momentum_eV=p0)
        P('%g | %.3g | 解析 adts().total | %.4f %.4f %.4f' % (th, r.res['coupling'], *[a.total[k] for k in ('xx', 'xy', 'yy')]))
        P('%g | | .rotate (傾きの分) | %.4f %.4f %.4f   (θ²: %.4f %.4f %.4f)' % (th, *[a.rotate[k] for k in ('xx', 'xy', 'yy')], *[a.rotate[k] / th ** 2 for k in ('xx', 'xy', 'yy')]))
        P('%g | | ampnf D=5 | %.4f %.4f %.4f' % (th, *[r.adts[k] for k in ('xx', 'xy', 'yy')]))
        P('%g | | 追跡 + NAFF | %.4f %.4f %.4f' % (th, *[tr[k] for k in ('xx', 'xy', 'yy')]))
        P('%g | | 差 (解析 − 追跡) / 追跡 | %.2e %.2e %.2e   (傾きの分を入れない: %.2e %.2e %.2e)  [%.0f s]' % (
            th, *[(a.total[k] - tr[k]) / tr[k] for k in ('xx', 'xy', 'yy')], *[(a.total[k] - a.rotate[k] - tr[k]) / tr[k] for k in ('xx', 'xy', 'yy')], time.time() - t0))
        res = A.one_turn_series(lat, D=3, nslice=16, bend_cubic=True, nslice_bend=8)
        h = A.rdt_series(res)
        sk = ('h10101', 'h10011', 'h20010', 'h02100', 'h11010')
        P('%g | | スキューの h (|h|): ' % th + '  '.join('%s %.3e' % (k, abs(h[k])) for k in sk) + '  (h30000 %.3e)' % abs(h['h30000']))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fodo')
