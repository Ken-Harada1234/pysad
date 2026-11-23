# -*- coding: utf-8 -*-
"""担当 A2 (4): ずれた BEND (DX) の CALC (tbende の tchge) と追跡 (tbend の tbshift) の差。

K64 の tchge.f (tbrote) の出口は、中点の座標の (dx, dz) を出口の座標へ回すとき s0 (= sin(phig)) の項の符号が入口と逆で、
出口のドリフト tdrife1(-dz1) が逆向きになる (dz1 = dz c0 - dx s0。入口は dz c0 + dx s0)。DX だけのとき、出口のドリフトの長さが
2 dx sin(θ/2) だけ違い (DX=2 mm、15 度の BEND で 5.2e-4 m)、1 要素の x-px・y-py の行列と z が ドリフト 5e-4 m の分ずれる。
v2 は出口も入口と同じ回し方 (R(phig)(dx, dz)) にした。

根拠 (追跡も SAD も使わない): 硬い端の厳密な幾何 (磁石を、弦に垂直な中点の x 方向へ dx だけ平行移動した帯の中の円弧。
名目の出口の平面で座標を測る)。矩形磁石 (E1=E2=0.5) は面に沿ったずれなので、光学・軌道が何も変わらない (厳密)。
扇形 (E=0) では x-px 行列が DX=2 mm で 8.6e-6 変わり、出口の軌道は 6.76e-5, 1.12e-4。
条件: 3 GeV、L=1.2 m、15 度。SAD の ANGLE>0 は x の負の側へ曲がる (厳密な幾何は鏡像にして比べる)。
"""
import io
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from pysad import compat, ffscmd, lattice, track6n

TH = 0.2617993877991494
L = 1.2


def _exact(xin, ain, dX, e):
    """x-z 平面の厳密な幾何 (曲率 1/rho で +x 側に曲がる向き、入口 (0,0) で +z)。dX: 弦に垂直な方向への磁石の平行移動。
    e: 端の角の比 (E1=E2=e*ANGLE)。戻り: 名目の出口の平面上の (x, 傾き)"""
    rho = L / TH
    e1 = e * TH
    dv = dX * np.array([math.cos(TH / 2), -math.sin(TH / 2)])
    n1 = np.array([math.sin(e1), math.cos(e1)])
    P1 = dv.copy()
    Pe = np.array([rho - rho * math.cos(TH), rho * math.sin(TH)])
    n2 = np.array([math.sin(TH - e1), math.cos(TH - e1)])
    P2 = Pe + dv
    p0 = np.array([xin, 0.0])
    d0 = np.array([math.sin(ain), math.cos(ain)])
    q1 = p0 + np.dot(P1 - p0, n1) / np.dot(d0, n1) * d0
    phi1 = math.atan2(d0[0], d0[1])

    def pos(phi):
        return q1 + rho * np.array([math.cos(phi1) - math.cos(phi), math.sin(phi) - math.sin(phi1)])
    phi2 = brentq(lambda phi: np.dot(pos(phi) - P2, n2), phi1, phi1 + 0.6, xtol=1e-15)
    q2 = pos(phi2)
    d2 = np.array([math.sin(phi2), math.cos(phi2)])
    te = np.array([math.sin(TH), math.cos(TH)])
    q3 = q2 + np.dot(Pe - q2, te) / np.dot(d2, te) * d2
    return float(np.dot(q3 - Pe, [math.cos(TH), -math.sin(TH)])), phi2 - TH


def _exact_sad(dx, e, h=1e-7):
    """SAD の座標 (ANGLE>0 は -x 側へ曲がる) に直した (出口の軌道 (x, px)、x-px の 2x2 行列)。鏡像: x -> -x、dX -> -dx"""
    x0, a0 = _exact(0.0, 0.0, -dx, e)
    M = np.zeros((2, 2))
    for j in range(2):
        d = np.zeros(2)
        d[j] = h
        p = _exact(d[0], d[1], -dx, e)
        m = _exact(-d[0], -d[1], -dx, e)
        M[:, j] = [(p[0] - m[0]) / (2 * h), (p[1] - m[1]) / (2 * h)]
    return np.array([-x0, -a0]), M


def _bend(dx, e):
    src = f'''MOMENTUM = 3GEV;
BEND BE = (L={L} ANGLE={TH} E1={e} E2={e} DX={dx});
LINE X = (BE);
FFS;
USE X;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return [x for x in ffs.lat.elements if x.name == 'BE'][0], ffs


def _calc(el):
    return lattice._bend_element_cod_matrix(el, np.zeros(6))


@pytest.mark.parametrize('e', [0.5, 0.0, 0.25])
def test_calc_shifted_bend_matches_exact_geometry_v2(v2_mode, e):
    """ずれた BEND (DX=2 mm) の CALC の 1 要素の出口の軌道と行列 (ずれなしとの差) が、厳密な幾何と 2e-8 で一致する"""
    be, _ = _bend(2e-3, e)
    b0, _ = _bend(0.0, e)
    co, Mc = _calc(be)
    _co0, Mc0 = _calc(b0)
    orb, Mx = _exact_sad(2e-3, e)
    _o0, Mx0 = _exact_sad(0.0, e)
    assert abs(co[0] - orb[0]) < 2e-8 and abs(co[1] - orb[1]) < 2e-8, (co[:2], orb)
    assert np.abs((Mc - Mc0)[:2, :2] - (Mx - Mx0)).max() < 2e-8
    assert abs(co[4]) < 2e-5        # 行路: 厳密な幾何で矩形は 0、扇形は -8.9e-6 (K64 は +5.2e-4)


@pytest.mark.parametrize('e', [0.5, 0.0])
def test_calc_shifted_bend_k64_has_the_exit_drift_error(k64_mode, e):
    """k64: 出口のドリフトが逆向きで、x-px の行列が 2 dx sin(θ/2) = 5.2e-4 だけずれる (記録)"""
    be, _ = _bend(2e-3, e)
    b0, _ = _bend(0.0, e)
    _co, Mc = _calc(be)
    _co0, Mc0 = _calc(b0)
    _orb, Mx = _exact_sad(2e-3, e)
    _o0, Mx0 = _exact_sad(0.0, e)
    d = ((Mc - Mc0)[:2, :2] - (Mx - Mx0))[0, 1]
    assert abs(abs(d) - 2 * 2e-3 * math.sin(TH / 2)) < 4e-5, d


def test_calc_shifted_bend_equals_tracking_element_v2(v2_mode):
    """CALC の 1 要素 (cod と x-px、y-py の行列) が、追跡 (tbshift + tbendcore) の写像とその数値微分と 2e-8 で一致する (矩形・扇形)"""
    for e in (0.5, 0.0):
        be, ffs = _bend(2e-3, e)
        p = ffs.momentum()

        def trk(c):
            return track6n.track_particles_fast([be], np.array([c]), 1, p, 0.0, False, False, False, circ=2 * L)[0]
        co, Mc = _calc(be)
        c0 = np.zeros(6)
        M = np.zeros((6, 6))
        for j in range(6):
            d = np.zeros(6)
            d[j] = 1e-6
            M[:, j] = (trk(c0 + d) - trk(c0 - d)) / 2e-6
        assert np.abs(co[:4] - trk(c0)[:4]).max() < 5e-9
        assert np.abs(Mc[:4, :4] - M[:4, :4]).max() < 2e-8, (e, np.abs(Mc[:4, :4] - M[:4, :4]).max())


_RING = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=%(e)s E2=%(e)s);
BEND  BE = (L=1.2 ANGLE=0.2617993877991494 E1=%(e)s E2=%(e)s DX=%(dx)s);
MARK  ORIGIN = ();
LINE CELLQ = (QF D2 B D2 QD DS D1);
LINE CELLE = (QF D2 BE D2 QD DS D1);
LINE RING = (ORIGIN CELLE 23*CELLQ);
FFS;
USE RING; CELL; CALC;
$FORM = "16.10";
Print["NU ", Twiss["NX","$$$"]/(2*Pi), " ", Twiss["NY","$$$"]/(2*Pi)];
Exit;
'''


def _tracked_tunes(els, p):
    circ = sum(x.L for x in els)

    def f(c):
        return np.asarray(track6n.track_particles_fast(els, np.array([list(c) + [0.0, 0.0]]), 1, p, 0.0, False, False, False, circ=circ))[0][:4]

    def jac(c, h=1e-6):
        J = np.zeros((4, 4))
        for i in range(4):
            d = np.zeros(4)
            d[i] = h
            J[:, i] = (f(c + d) - f(c - d)) / (2 * h)
        return J
    v = np.zeros(4)
    for _ in range(8):
        v = v - np.linalg.solve(jac(v) - np.eye(4), f(v) - v)
    ev = np.linalg.eigvals(jac(v))
    return sorted(set(round(abs(math.atan2(z.imag, z.real)) / 2 / math.pi, 10) for z in ev))


def _tune_diff(dx, e):
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(_RING % dict(dx=dx, e=e), out=buf)
    t = [ln for ln in buf.getvalue().splitlines() if ln.startswith('NU ')][0].split()
    calc = sorted(min(float(x) % 1.0, 1 - float(x) % 1.0) for x in (t[1], t[2]))
    trk = _tracked_tunes(ffs.lat.elements, ffs.momentum())
    return max(abs(a - b) for a, b in zip(calc, trk))


def test_ring_tune_with_one_shifted_rectangular_bend_v2(v2_mode):
    """矩形磁石 (E=0.5) 1 台を DX=2 mm ずらしたリング: CALC のチューン = 追跡の一周写像の固有値 (1e-9)。厳密な幾何でも光学は変わらない"""
    assert _tune_diff(2e-3, 0.5) < 1e-9
    assert _tune_diff(-3e-3, 0.5) < 1e-9


def test_ring_tune_with_one_shifted_rectangular_bend_k64(k64_mode):
    """k64: CALC のチューンが 2.7e-5 (2 mm) 動く (追跡は動かない)"""
    d = _tune_diff(2e-3, 0.5)
    assert 2e-5 < d < 4e-5, d


def test_ring_tune_with_shifted_sector_bend_v2(v2_mode):
    """扇形 (E=0) では厳密な幾何でも光学が動く (8.6e-6 の行列の差) が、CALC と追跡が 1e-8 で一致する"""
    assert _tune_diff(2e-3, 0.0) < 1e-8


# --- 大きなリング A (3 GeV、周長 603 m、kek-ls_v03_12.lat): K1=0 の偏向磁石 30 台に DX ---

def _ring_a_diff(mode):
    import os
    lat = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'kek-ls_v03_12.lat')
    if not os.path.exists(lat):
        pytest.skip('kek-ls_v03_12.lat がない')
    old = compat._mode
    compat.set_mode(mode)
    try:
        with open(lat, encoding='utf-8') as fh:
            src = fh.read()
        ffs, _ = ffscmd.run_sad_string(src + '\nFFS;\nUSE RING;\nCELL; CALC;\nExit;\n', out=io.StringIO())
        els = ffs.lat.elements
        rng = np.random.default_rng(1)
        bends = [e for e in els if type(e).__name__ == 'Bend' and e.ANGLE != 0 and e.K1 == 0]
        for i in rng.choice(len(bends), 30, replace=False):
            bends[i].DX = float(rng.uniform(-5e-4, 5e-4))
        ffs.calc(quiet=True)
        n = len(els) + 1
        calc = sorted(min((ffs.twiss_value(k, n) / (2 * math.pi)) % 1.0, 1 - (ffs.twiss_value(k, n) / (2 * math.pi)) % 1.0) for k in ('NX', 'NY'))
        trk = _tracked_tunes(els, ffs.momentum())
        return max(abs(a - b) for a, b in zip(calc, trk))
    finally:
        compat.set_mode(old)


def test_large_ring_a_shifted_bends_calc_tune_matches_tracking_v2(v2_mode):
    """大きなリング A に DX (±0.5 mm、30 台) を入れた CALC のチューンが、追跡の一周写像と 1e-7 以内。K64 は 2 倍以上ずれる"""
    d2 = _ring_a_diff('v2')
    dk = _ring_a_diff('k64')
    assert d2 < 1e-7, d2
    assert dk > 2 * d2, (dk, d2)
