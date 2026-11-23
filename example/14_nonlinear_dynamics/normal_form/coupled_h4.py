# -*- coding: utf-8 -*-
"""結合した光学 (ROTATE した磁石) の共鳴駆動項の拡張 (4 次の h・δ≠0・ref・DX/DY/DROTATE。`pysad/ampnf_coupled.py`) の検証。
使い方: PYTHONUTF8=1 py -3.12 example/14_nonlinear_dynamics/normal_form/coupled_h4.py fodo|astrid   (出力は標準出力。保存は coupled_h4_*.out.txt)
条件: v2 モード (DX・DROTATE のある BEND の tchge は v2 の符号。k64 では BEND の DX が追跡と 4e-8 違う)、RF・放射なし、追跡の写像は ring_map3 (1 周の写像の 3 次までのテイラー展開、
結合したモードの正規化、閉軌道のまわり)。表の値は最大の項に対する相対差 (最大の項の 1e-3 倍以上の項)。ADTS は h3・h4 の正規形 (nform.adts_from_map) と
ampnf の正規形 (D=5) と追跡 + NAFF のフィット (ampnf.tracked_adts) [1/(m rad)]。
fodo: 3 GeV の 24 セル FODO (K2 = 4・−6)、先頭の 1 台を 50 mrad 傾ける。
astrid: ASTRID2 (0.58 GeV、45.7 m、複合機能 BEND)、全 BEND の ROTATE = 20 mrad。"""
import copy
import dataclasses
import math
import os
import sys
import time

import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from pysad import ampnf as A, ampnf_coupled as C, compat, nform, rdttrack as T
from pysad.lattice import Lattice
from pysad.parser import load_sad_file, load_sad_string
compat.set_mode('v2')
P = lambda *a: print(*a, flush=True)
OPT = dict(nslice=16, bend_cubic=True, kin=True, fringe=True)


def fodo():
    src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
        'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
    ANG, EDG = 0.2617993877991494, 0.1308996938995747
    src = src.replace('BEND  B  = (', 'BEND BR = (ROTATE=0.05 L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g)\nBEND  B  = (' % (ANG, EDG, EDG))
    src = src.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE CELLR = (SF DS QF D2 BR D2 QD DS SD D1);\nLINE RING = (ORIGIN RF 1*CELLR 23*CELL);')
    return load_sad_string(src).build_lattice('RING'), 3e9, 4, 3e-7


def astrid():
    lat = load_sad_file(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')).build_lattice('RING')
    out = []
    for el in lat.elements:
        if type(el).__name__ == 'Bend' and float(el.ANGLE) != 0.0 and float(el.L) > 0:
            el = dataclasses.replace(el, ROTATE=20e-3)
        out.append(el)
    return Lattice(out, name=lat.name), 0.58e9, 4, 3e-9


def with_attrs(lat, spec):
    els = [copy.copy(e) for e in lat.elements]
    for tn, lst in spec.items():
        idx = [i for i, e in enumerate(els) if ((type(e).__name__ == 'Bend' and float(e.ANGLE) != 0.0 and float(e.L) > 0) if tn == 'BendLive' else type(e).__name__ == tn)]
        for k, attrs in lst:
            for a, v in attrs.items():
                setattr(els[idx[k % len(idx)]], a, v)
    return Lattice(els, name=lat.name)


def rel(h, ht, frac=1e-3):
    m = max(abs(v) for v in ht.values())
    return max(abs(h[k] - ht[k]) for k in ht if abs(ht[k]) > frac * m) / m, m


def case(label, lat, p0, D, delta=0.0, start=None, adts=False, jn=3e-7, amp3=2e-4):
    t0 = time.time()
    rm = C.taylor_map(lat, delta=delta, start=start, D=D, **OPT, nslice_bend=4)
    t1 = time.time() - t0
    els = list(lat.elements)
    if start is not None:
        els = els[start:] + els[:start]
    rmt = T.ring_map3(els, p0, rm.A, delta=delta, amp=1e-5, amp3=amp3)
    h3, (h4, info) = C.rdt3(rm), C.rdt4(rm)
    t3 = T.rdt_from_map(rmt)
    t4, _ = T.rdt4_from_map(rmt)
    e3, m3 = rel(h3, t3)
    e4, m4 = rel(h4, t4)
    P('%-34s D=%d  閉軌道 %.2e m (追跡との差 %.1e)  結合 %.3f | h3 max %.3g 差 %.1e | h4 max %.4g 差 %.1e | 模型 %.0f s' % (
        label, rm.D, np.abs(rm.orbit).max(), np.abs(rm.orbit - rmt.orbit).max(), rm.coupling, m3, e3, m4, e4, t1))
    if adts:
        a_map, _ = nform.adts_from_map(rm)
        r = A.amplitude_tune(lat, D=5, nslice=16, bend_cubic=True, kin=True, fringe=True, nslice_bend=4)
        a_tr = A.tracked_adts(lat, r, j=jn, momentum_eV=p0)
        a_t, _ = nform.adts_from_map(rmt)
        P('   ADTS xx xy yy: h3,h4 の正規形 %s | ampnf D=5 %s | 追跡の写像 %s | 追跡 + NAFF %s' % tuple(
            ' '.join('%.3f' % d[k] for k in ('xx', 'xy', 'yy')) for d in (a_map, r.adts, a_t, a_tr)))
        P('   相対差: 正規形 − ampnf %s | 正規形 − 追跡 + NAFF %s' % tuple(' '.join('%.1e' % (a_map[k] / d[k] - 1) for k in ('xx', 'xy', 'yy')) for d in (r.adts, a_tr)))
    return rm


def main(name):
    lat, p0, D, jn = {'fodo': fodo, 'astrid': astrid}[name]()
    P('# %s (v2 モード)。3 次・4 次の h を、模型の 1 周の写像 (ampnf_coupled) と追跡の写像のテイラー展開 (rdttrack.ring_map3) から出して比べる' % name)
    case('傾き (ROTATE) だけ', lat, p0, D, adts=True, jn=jn, amp3=2e-4 if name == 'fodo' else 1e-4)
    case('δ = +1e-3', lat, p0, 5, delta=1e-3, amp3=2e-4 if name == 'fodo' else 1e-4)
    case('δ = −1e-3', lat, p0, 5, delta=-1e-3, amp3=2e-4 if name == 'fodo' else 1e-4)
    case('始点 = 素子 30 (1 周の写像の始点)', lat, p0, D, start=30, amp3=2e-4 if name == 'fodo' else 1e-4)
    spec = dict(Quadrupole=[(2, dict(DX=1.5e-4, DY=-1e-4)), (7, dict(DX=-2e-4, DY=1.2e-4)), (12, dict(DX=1e-4, DY=2e-4))],
                Sextupole=[(3, dict(DX=2e-4, DY=-1e-4)), (10, dict(DX=-1e-4, DY=3e-4))], BendLive=[(5, dict(DX=2e-4, DY=1e-4))])
    case('DX・DY (四極 3・6 極 2・BEND 1 台)', with_attrs(lat, spec), p0, D, amp3=2e-4 if name == 'fodo' else 1e-4)
    case('BEND の DROTATE (3 台)', with_attrs(lat, dict(BendLive=[(1, dict(DROTATE=3e-3)), (7, dict(DROTATE=-2e-3)), (13, dict(DROTATE=4e-3))])), p0, D, amp3=2e-4 if name == 'fodo' else 1e-4)
    spec2 = dict(spec, BendLive=[(5, dict(DX=2e-4, DY=1e-4, DROTATE=3e-3))])
    case('DX・DY + DROTATE + δ = 5e-4', with_attrs(lat, spec2), p0, 5, delta=5e-4, amp3=2e-4 if name == 'fodo' else 1e-4)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fodo')
