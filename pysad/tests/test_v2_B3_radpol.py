# -*- coding: utf-8 -*-
"""担当 B3 (2): K64 の RADPOL 追跡 (photons.f の sprot の放射偏極、dst = (PST − st)·cfpd·(p p0/h²)²·ucde/bt) の緩和が Sokolov-Ternov の値と合うか。

A2 の記録: ASTRID2 (0.58 GeV) で P=1 の偏極の緩和が (1 − PST)/τ の 1.58 倍速い (−4.17e-6 /s 対 −2.65e-6 /s)。
結論 (docs/history/v2next_B3.md §3): **K64 の式の誤りではない**。
  (a) 式: tradk の 1 刻みの ucde = uc·de = ħ r_e γ^5 ds/ρ³ で、cfpd·ucde = (5√3/8) r_e ƛ_C γ^5 ds/ρ³ (ST の率そのもの。de = (2/3) r_e γ³ ds/ρ²、uc = (3/2) ƛ_C γ²/ρ)。
      dst = (PST − st)·率 を n の向きに足すので、P=1 からの緩和は (1 − PST)×率、P=0 (n に垂直) からの立ち上がりは PST×率。
  (b) 3 GeV の FODO (τ = 143 s、1 刻みの増分 ~1e-12): 追跡の P=1 からの緩和と P=0 からの立ち上がりが、式の −(1−PST)·sdamp、+PST·sdamp と 1e-4 以内で一致する (K64 の式どおり)。
  (c) ASTRID2 (τ = 28815 s): 1 刻みの増分が (1−PST)·cfpd·ucde = 6e-17 で、倍精度の |s|≈1 の刻み (ulp = 1.1e-16、下側 5.5e-17) の半分より大きい。
      s = 1.0 に −6.4e-17 を足すと 1 − 1.1e-16 に丸まる (1.7 倍)。ASTRID2 の 1 周 6336 刻みで、この丸めが 1.58 倍の緩和を作る (追跡の値は 1 − |s| = 1 − P の増分だけ)。
      P=0.5 から始めると (増分 3.5e-16 = 3 ulp) 丸めが目立たず、式と合う。
v2・k64 とも追跡は変えない。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, track6

PST = 8.0 / (5.0 * math.sqrt(3.0))

_FODO = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
USE RINGQ; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''


@pytest.fixture(scope='module')
def ring():
    compat.set_mode('v2')
    try:
        ffs, _ = ffscmd.run_sad_string(_FODO, out=io.StringIO())
        p = ffs.momentum()
        els = ffs.lat.elements
        e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
    finally:
        compat.set_mode('k64')
    return dict(p=p, els=els, e=e, sp=e['spin'], circ=emit6d._emit_circ(els), cod=np.array(e['orbit_entrance']), trf0=e.get('trf0', 0.0))


def _track(r, s0, nt, radpol=True):
    _o, s_ = track6.track_particles(r['els'], np.array([r['cod']]), nt, r['p'], r['trf0'], True, True, False, spins=[s0], radpol=radpol, circ=r['circ'])
    return np.array(s_[0])


@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_radpol_relaxation_equals_sokolov_ternov(ring, mode):
    """P=1 (n の向き) からの緩和 −(1−PST)·sdamp、P=0 (n に垂直) からの立ち上がり +PST·sdamp /周 を、追跡が 1e-4 で再現する (3 GeV FODO、200 周)"""
    compat.set_mode(mode)
    try:
        sd = ring['sp']['sdamp']
        n0 = np.array([0.0, 1.0, 0.0])
        s = _track(ring, n0, 200)
        assert abs((s @ n0 - 1.0) / 200 / (-(1 - PST) * sd) - 1) < 1e-4
        s = _track(ring, np.array([1.0, 0.0, 0.0]), 200)
        assert abs(s @ n0 / 200 / (PST * sd) - 1) < 1e-4
        s = _track(ring, n0, 200, radpol=False)
        assert abs(1 - s @ n0) < 1e-14                    # RADPOL なしでは偏極は変わらない
    finally:
        compat.set_mode('k64')


def test_radpol_increment_is_below_double_precision_for_low_energy_rings():
    """ASTRID2 の 1 刻みの緩和 (1−PST)·cfpd·ucde = 6.4e-17 を 1.0 に足すと 1 − 1.1e-16 に丸まる (倍精度の下側の ulp = 1.1e-16、半分 5.5e-17 より大きい)。
    増えた量の比 (1.73) が、ASTRID2 の 6336 刻み/周のうち効く刻みの混ざり方で 1.58 になる"""
    inc = 6.4e-17
    s = np.float64(1.0)
    n = 100000
    for _ in range(n):
        s = s - np.float64(inc)
    assert abs((1.0 - s) / (n * inc) - 1.73) < 0.02
    s = np.float64(0.5)                        # |s|<1/2 で ulp は半分、増分が 3 ulp 以上なら相対 <20% の誤差
    big = 3.5e-16
    for _ in range(n):
        s = s - np.float64(big)
    assert abs((0.5 - s) / (n * big) - 1.0) < 0.2
