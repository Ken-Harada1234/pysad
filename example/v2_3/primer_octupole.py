# -*- coding: utf-8 -*-
"""入門 §5: 8 極の 1 次の ADTS の式 (薄い 8 極で ∂νx/∂Jx = K3L βx²/(16π)、∂νx/∂Jy = −K3L βx βy/(8π)、∂νy/∂Jy = K3L βy²/(16π)) を、
FODO リング (3 GeV、105.6 m、6 極なし) のリング先頭に薄い 8 極 (MULT、L = 0、K3L = 1e4 [1/m³]) を 1 個置いて確かめる。
  (a) 式 (リング先頭の Twiss の βx、βy) (b) pysad の解析 adts().octu (c) 追跡 + NAFF + フィット (8 極あり − なし、Jmax = 1e-7、1024 周、order = 2)
使い方 (リポジトリ直下から): py -3.12 example/v2_3/primer_octupole.py   (所要 約 20 s)
"""
import math
import numpy as np
import primer_common as C
from pysad import adtsdirect as D, resdrive as R
from pysad.parser import load_sad_string

K3L = 1e4
src0 = open(C.os.path.join(C.ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
src = src0.replace('LINE RING = (ORIGIN RF 24*CELL);', 'MULT OC = (L=0 K3=%g);\nLINE RING = (ORIGIN OC RF 24*CELL);' % K3L)
assert 'OC' in src
res = {}
for tag, s in (('なし', src0), ('あり', src)):
    lat = load_sad_string(s).build_lattice('RING')
    rows = lat.twiss_table()
    r0 = rows[0]
    bx, by = float(r0.betx), float(r0.bety)
    a = R.adts(lat, rows=rows)
    t = D.direct_adts(lat, C.P_FODO, 1e-7, nturn=1024, order=2, rows=rows)
    res[tag] = dict(bx=bx, by=by, octu=np.array([a.octu['xx'], a.octu['xy'], a.octu['yy']]), tr=np.array([t['xx'], t['xy_mean'], t['yy']]))
bx, by = res['あり']['bx'], res['あり']['by']
form = np.array([K3L * bx ** 2 / (16 * math.pi), -K3L * bx * by / (8 * math.pi), K3L * by ** 2 / (16 * math.pi)])
tr = res['あり']['tr'] - res['なし']['tr']
print('FODO リング (3 GeV、105.6 m、6 極なし)、リング先頭に薄い 8 極 K3L = %g [1/m³]。βx = %.4f m、βy = %.4f m (単位 1/(m rad)、(xx, xy, yy))' % (K3L, bx, by))
print('(a) 式                     %12.3f %12.3f %12.3f' % tuple(form))
print('(b) pysad 解析 adts().octu %12.3f %12.3f %12.3f   差/|式| = %.1e %.1e %.1e' % (*res['あり']['octu'], *((res['あり']['octu'] - form) / np.abs(form))))
print('(c) 追跡 (あり − なし)     %12.3f %12.3f %12.3f   差/|式| = %.1e %.1e %.1e' % (*tr, *((tr - form) / np.abs(form))))
