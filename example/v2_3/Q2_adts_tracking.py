# -*- coding: utf-8 -*-
"""v2.3.0 Q2: pysad の追跡 (track6n) と NAFF で ADTS を求め、解析値 (pysad/resdrive.py) と比べる。gs (SAD) での同じ比較は Q2_adts_gs.py。
条件: 互換モード k64 (gs と同じ追跡)、RF・放射なし、1024 周、リング先頭から。6 極の強さを λ 倍 (λ = +1, −1, 0)、8 極を μ 倍 (μ = ±1) の追跡から
  6 極の 2 次の項 C = (f(1,1) + f(−1,1) − 2 f(0,1))/2 、8 極の 1 次の項 B = (f(0,1) − f(0,−1))/2
を取り出す。(λ の奇数次 (6 極 1 次 × 偏向磁石などの非線形) と 6 極以外の非線形を消すため。f は NAFF のチューンの中心差分から、振幅 b と 2b の外挿)
使い方 (リポジトリ直下から): py -3.12 example/v2_3/Q2_adts_tracking.py fodo|astrid|kekls"""
import os, sys, time
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, compat
from pysad.parser import load_sad_string
compat.set_mode('k64')
sys.path.insert(0, os.path.join(ROOT, 'example', 'v2_3'))
import Q2_adts_gs as G
key = sys.argv[1] if len(sys.argv) > 1 else 'fodo'
nturn = 1024
src, line, p, jb = G.read_ring(key)
lat = load_sad_string(src).build_lattice(line)
rows = lat.twiss_table()
r0 = rows[0]
tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
sx = [e for e in lat.elements if type(e).__name__ == 'Sextupole']
oc = [e for e in lat.elements if type(e).__name__ == 'Octupole']
k2, k3 = [e.K2 for e in sx], [e.K3 for e in oc]
A = R.adts(lat, rows=rows)        # 解析値は強さを変える前に求める


def f(lam, mu):
    for e, k in zip(sx, k2):
        e.K2 = k * lam
    for e, k in zip(oc, k3):
        e.K3 = k * mu
    t = R.adts_tracking(lat.elements, p, tw, jref=(jb, jb), nturn=nturn)
    return np.array([t['xx'], t['xy'], t['yy']])


t0 = time.time()
res = {c: f(*c) for c in [(1, 1), (-1, 1), (0, 1), (0, -1)] if oc or c[1] == 1}
C = (res[(1, 1)] + res[(-1, 1)] - 2 * res[(0, 1)]) / 2
S = np.array([A.sext['xx'], A.sext['xy'], A.sext['yy']])
print('リング %s: エネルギー %.3g GeV、周長 %.3f m、νx=%.6f νy=%.6f、J の基準 %.1e m rad、%d 周、k64 (経過 %.0f s)' % (
    key, p / 1e9, rows[-1].s, rows[-1].mux / 6.283185307179586, rows[-1].muy / 6.283185307179586, jb, nturn, time.time() - t0))
print('6 極の 2 次  追跡 C   = %s' % np.array2string(C, precision=5))
print('            解析     = %s    差/|解析| = %s' % (np.array2string(S, precision=5), np.array2string((C - S) / np.abs(S), formatter={'float_kind': lambda v: '%.1e' % v})))
if oc:
    B = (res[(0, 1)] - res[(0, -1)]) / 2
    O = np.array([A.octu['xx'], A.octu['xy'], A.octu['yy']])
    print('8 極の 1 次  追跡 B   = %s' % np.array2string(B, precision=5))
    print('            解析     = %s    差/|解析| = %s' % (np.array2string(O, precision=5), np.array2string((B - O) / np.abs(O), formatter={'float_kind': lambda v: '%.1e' % v})))
print('全部入り (λ=1, μ=1) の追跡 = %s ; 解析の和 (6 極 2 次 + 8 極 1 次) = %s ; 差は 6 極×8 極・偏向磁石などの交差項' % (
    np.array2string(res[(1, 1)], precision=5), np.array2string(S + (O if oc else 0), precision=5)))
print('(0, 1) の追跡 f(0,1) = %s  (6 極なし。8 極 1 次 + 偏向磁石・端部の非線形)' % np.array2string(res[(0, 1)], precision=5))
