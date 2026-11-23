# -*- coding: utf-8 -*-
"""K64 の tcod (6D 閉軌道) が収束しないときの Emittance[] (B163)。

条件: DDBA_15 のラティス (5 GeV、RF 7 MV、6 極 600 台) のノーマルセルだけのリング RINGM (1 周 656.4 m)、ON RAD FLUC COD RFSW RADCOD。
RFSW ON の最初の Newton 法で、運動量の誤差が大きい (dpthre 超え) ために横方向の修正が 0 にされ、その方向が以後使われ続けて r が下がらず、
tcod は fact < fmin で終わる。K64 の temit は "???-Emittance[]-closed orbit not found." を出し、収束していない cod から続行する (Stable = 0)。
- k64: 同じに続行する。gs の Emittance[] (2026-11-22 の gs、同じ台本) と全桁一致 (定数は gs の出力)。
- v2: 放射込みの一周写像の 6D 固定点を Newton 法で解く。固定点の残差が 1e-12 以下で、Stable = 1。
6 極を 0 にした同じリングは両モードで収束し gs と一致する。"""
import io
import os
import re

import numpy as np
import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
LAT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), 'example', 'ddba', 'lattice', 'DDBA_15_v05_21.sad')
NOSX = 'Do[If[LINE["TYPENAME",i]=="SEXT",LINE["K2",i]=0],{i,1,LINE["LENGTH"]}];\n'
# gs (K64) の Emittance[] (RINGM、6 極つき): EnergyLossU0, EquilibriumPosition, MomentumCompaction, MomentumSpread, Emittances[1], Stable
GS_SEXT = dict(U0=2412450.636949, dz=0.0301668200, alpha=1.340803290e-4, sigE=0.0010139795, ex=6.27329883e-10, stable=0.0)


def _run(extra=''):
    if not os.path.exists(LAT):
        pytest.skip('DDBA_15 のラティスがない')
    src = open(LAT, encoding='utf-8').read()
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\nON LOG;\nON RAD FLUC COD RFSW RADCOD;\nFFS;\n$FORM = "15.10"; PageWidth = 1999;\nUSE RINGM;\n' + extra +
                          'cell calc emit;\nem=Emittance[];\nPrint["EMRES ",em];\nExit;\n', out=buf)
    t = buf.getvalue()
    line = ''.join(t[t.index('EMRES'):].split(chr(10))[:12]).replace(' ', '')
    def g(k):
        return float(re.search(k + r'->([-\d.E+]+)', line).group(1))
    ex = float(re.search(r'Emittances->\{([-\d.E+]+)', line).group(1))
    return dict(U0=g('EnergyLossU0'), dz=g('EquilibriumPosition'), alpha=g('MomentumCompaction'),
                sigE=g('MomentumSpread'), ex=ex, stable=g('Stable')), t


def test_k64_continues_from_unconverged_orbit_like_gs(k64_mode):
    r, t = _run()
    assert 'closed orbit not found' in t
    assert r['stable'] == 0.0
    for k in ('U0', 'dz', 'alpha', 'sigE', 'ex'):
        assert r[k] == pytest.approx(GS_SEXT[k], rel=2e-9), k


def test_v2_solves_fixed_point(v2_mode):
    r, t = _run()
    assert 'closed orbit not found' not in t
    assert r['stable'] == 1.0
    # k64 (収束していない軌道 + K64 の他の癖) と εx が 1〜3e-3 違う。固定点の側が正しい (docs/30 B163)
    assert 1e-3 < abs(r['ex'] / GS_SEXT['ex'] - 1.0) < 3e-3


def test_no_sextupole_both_modes_converge(k64_mode):
    r, t = _run(NOSX)
    assert 'closed orbit not found' not in t
    assert r['stable'] == 1.0
    assert r['ex'] == pytest.approx(7.74473888e-10, rel=1e-9)     # gs (6 極なし)
    compat.set_mode('v2')
    try:
        r2, _ = _run(NOSX)
    finally:
        compat.set_mode('k64')
    assert r2['ex'] == pytest.approx(r['ex'], rel=2e-3)
