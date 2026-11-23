# -*- coding: utf-8 -*-
"""ずれて傾いたソレノイド (BZ = 2 T、長さ 0.8 m の DRIFT 2 本) のあるリングの Emittance[] の POL の量 (SpinTune・n0・srot) を、追跡の 1 周のスピン行列と比べる。
ソレノイドの領域をリングの先頭に足す (GEO=1 の SOL に DX DY DPX DPY: 設計軌道の位置・向きをソレノイドの系で書く)。Emittance[] は NORFSW NORAD で POL ON、
追跡は track6 (POL + RAD、放射の損失は 1e-9 倍にしてスピンの歳差だけ見る。RFSW OFF)。
「新」= v2 の Emittance[] (運動量の向きの磁場の積分 B∥ を向きの変化から作る B⊥ と分けて足す)、「K64 式」= spin.SOL_BPAR = False (K64 の bsi を z 軸の向きに足す式)。
リング: 小さな FODO (0.5 GeV、周長 45.6 m)、ASTRID2 (0.58 GeV、45.7 m)、ASTRID2 + MPW、PF リング (2.5 GeV、187 m)、大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)。
使い方: py -3.12 example/16_spin_polarization/spin_sol_tilt.py (約 1 分)。結果は spin_sol_tilt.out.txt。"""
import io
import math
import os
import sys

import numpy as np

import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
from pysad import compat, emit6d, ffscmd, spin, track6

compat.set_mode('v2')
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
FODO = '''MOMENTUM = 0.5 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELLQ = (QF D2 B D2 QD D1);
MARK MKS = ();
LINE RING = (MKS 12*CELLQ);
'''
RINGS = [('小さな FODO (0.5 GeV、45.6 m)', None, 'RING'),
         ('ASTRID2 (0.58 GeV、45.7 m)', 'example/astrid/lattice/ASTRID2_v01_01.sad', 'RING'),
         ('ASTRID2 + MPW', 'example/astrid/lattice/ASTRID2_v01_02_id.sad', 'RINGMPW'),
         ('PF リング (2.5 GeV、187 m)', 'example/pf/lattice/v9_16_asym.sad', 'RING'),
         ('大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)', 'pysad/tests/kek-ls_v03_12.lat', 'RING')]
GEOS = [('軸上', ''), ('傾き 3 mrad・ずれ 2 mm', 'DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002'), ('傾き 10 mrad・ずれ 4 mm', 'DX=0.004 DY=-0.002 DPX=0.01 DPY=0.006')]


def run(base, line, geo, flag):
    spin.SOL_BPAR = flag
    extra = f'''
SOL SA = (L=0 BZ=2 BOUND=1 GEO=1 {geo}) SB = (L=0 BZ=0 BOUND=1);
DRIFT DSS = (L=0.4);
LINE SOLR = (SA DSS DSS SB);
LINE RINGS = (SOLR {line});
FFS;
USE RINGS; NORFSW; NORAD; NOFLUC; POL;
Exit;
'''
    ffs, _ = ffscmd.run_sad_string(base + extra, out=io.StringIO())
    ffs.run_emit(quiet=True)
    sp = ffs.emit['spin']
    els = ffs.lat.elements
    emit6d.CAVE *= 1e-9
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _o, s1 = track6.track_particles(els, np.zeros((1, 6)), 1, ffs.momentum(), 0.0, False, True, False, spins=[np.array(s0, float)], radpol=False,
                                        circ=emit6d._emit_circ(els))
        M[:, j] = s1[0]
    emit6d.CAVE /= 1e-9
    nu_tr = math.acos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))) / (2 * math.pi)
    w, v = np.linalg.eig(M)
    k = int(np.argmin(abs(w - 1)))
    n0 = np.real(v[:, k])
    n0 /= np.linalg.norm(n0)
    if n0 @ sp['sps'][:, 0] < 0:
        n0 = -n0
    return abs(sp['nup']), nu_tr, math.acos(min(1.0, n0 @ sp['sps'][:, 0])), np.abs(sp['srot1'] - M).max()


print('リング / ソレノイドの置き方 / 式 : SpinTune (Emittance[])  SpinTune (追跡)  差    n0 の向きの差 [rad]  |srot1 − 追跡の行列|')
for name, path, line in RINGS:
    base = FODO if path is None else open(os.path.join(ROOT, path), encoding='utf-8', errors='replace').read()
    for gname, geo in GEOS:
        for flag in (True, False):
            nu_em, nu_tr, ang, dm = run(base if path else FODO.replace('MKS', 'MKS'), line, geo, flag)
            print('%-44s %-22s %-6s %.9f  %.9f  %+.2e  %.2e  %.2e' % (name, gname, '新' if flag else 'K64式', nu_em, nu_tr, nu_em - nu_tr, ang, dm), flush=True)
