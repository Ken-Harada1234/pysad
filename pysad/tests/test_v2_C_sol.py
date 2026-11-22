# -*- coding: utf-8 -*-
"""担当 C: ソレノイド領域 (BOUND≠0 の SOL の間) の中の BEND・MULT・CAVI・APERT (docs/21 B69、docs/46)、MULT のスキュー成分。

K64 の tsol.f は領域内の BEND を tdrift(al, bzs, phiy, phix)、MULT を tmulti1、CAVI を tmultiacc (どれも bz つき) で進めるが、
gs (この K64 のビルド) は領域内の偏向をほぼ無視する (B69)。pysad の v1 は未対応 (NotImplementedError) だった。v2 は実装する。

根拠 (gs は根拠にならないので、独立な数値計算と極限):
  - 一様磁場 B = (Bx, By, bz) と多極磁場 By + i Bx = Σ (K_n + i SK_n)/L (x+iy)^n/n! の中のローレンツ力 dp/ds = (p × B)/pz を
    scipy の solve_ivp で積分した結果と、追跡のカーネルが一致 (力学的運動量で比較)。
  - 放射損失 (一様な双極 + ソレノイド): Δδ = -(2/3) r_e γ³ (K0/L)² L。
  - bz → 0 の領域内の要素は、領域の外の普通の要素と同じ光学・エミッタンス (リングの CALC・Emittance[])。
  - CALC (INS) の軌道が追跡と一致 (領域内の BEND の蹴りが光学にも入る)。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, emit6d as _e, track6n as tn, track6n_ext as x

P = 1.0e9
RE = 2.8179403e-15


def _elem(defs):
    import pysad.ffscmd as F
    src = 'MOMENTUM = 1GEV;\n' + defs + '\nMARK ST=();\nLINE BL=(ST E1);\nFFS;\nUSE BL;\nExit;\n'
    holder = {}
    orig = F.FFS.use

    def use(self, *a, **k):
        r = orig(self, *a, **k)
        holder['f'] = self
        return r
    F.FFS.use = use
    try:
        F.run_sad_string(src, out=io.StringIO())
    finally:
        F.FFS.use = orig
    return [e for e in holder['f'].els if e.name == 'E1'][0]


class _Col:
    def __init__(self):
        self.items = []

    def _add(self, k, p, name, cnt, sm):
        self.items.append((k, p))


def _trk(el, c, bz, rad=False, rfsw=False):
    ph = _e._Phys(P, 100.0)
    col = _Col()
    x._compile_sol_v2(col, el, type(el).__name__, bz, ph, 2)
    k, p = col.items[0]
    S = tn.make_state(np.array([c], float), ph)
    x.step_kind(k, p, S, 1, ph, 0.0, rfsw, rad, False, tn.make_rs(1))
    return tn.state_out(S, 1)[0], ph


def _kin(v, bz):
    return np.array([v[0], v[1] + bz * v[2] / 2, v[2], v[3] - bz * v[0] / 2])


def _ode(c, L, bz, K, SK=None):
    x0, px, y0, py = c[:4]
    pn = 1.0 + c[5]

    def f(s, u):
        X, PX, Y, PY = u
        pz = math.sqrt(pn * pn - PX * PX - PY * PY)
        w = complex(X, Y)
        b = 0j
        t = 1 + 0j
        fct = 1.0
        for n in range(len(K)):
            if n > 0:
                fct *= n
            b += complex(K[n], SK[n] if SK else 0.0) / L * t / fct
            t *= w
        By, Bx = b.real, b.imag
        return [PX / pz, (PY * bz - pz * By) / pz, PY / pz, (pz * Bx - PX * bz) / pz]
    return solve_ivp(f, [0, L], [x0, px, y0, py], rtol=1e-13, atol=1e-16).y[:, -1]


C0 = np.array([1e-3, 2e-3, -1.5e-3, 1e-3, 0.0, 0.01])
CASES = [
    ('BEND E1=(L=0.5 K0=0.02 ROTATE=0.4);', [0.02 * math.cos(0.4)], [0.02 * math.sin(0.4)], 1e-14),
    ('BEND E1=(L=0.5 ANGLE=0.02);', [0.02], None, 1e-14),
    ('MULT E1=(L=0.5 K0=0.02 DISFRIN=1);', [0.02], None, 1e-14),
    ('MULT E1=(L=0.5 K1=0.3 DISFRIN=1);', [0, 0.3], None, 1e-10),
    ('MULT E1=(L=0.5 K1=-0.3 DISFRIN=1);', [0, -0.3], None, 1e-10),
    ('MULT E1=(L=0.5 SK1=0.2 DISFRIN=1);', [0, 0], [0, 0.2], 1e-10),
    ('MULT E1=(L=0.5 K0=0.01 K1=0.3 SK0=0.005 DISFRIN=1);', [0.01, 0.3], [0.005, 0], 5e-8),
    ('MULT E1=(L=0.5 SK2=30 K1=0.1 SK1=0.05 DISFRIN=1 EPS=0.01);', [0, 0.1, 0], [0, 0.05, 30], 1e-8),
    ('MULT E1=(L=0.5 K0=0.01 K1=0.3 K2=30 K3=2000 DISFRIN=1 EPS=0.01);', [0.01, 0.3, 30, 2000], None, 1e-8),
]


@pytest.mark.parametrize('bz', [0.3, 3.0])
@pytest.mark.parametrize('case', range(len(CASES)))
def test_sol_region_element_matches_lorentz_force(v2_mode, bz, case):
    defs, K, SK, tol = CASES[case]
    o, _ = _trk(_elem(defs), C0, bz)
    r = _ode(list(_kin(C0, bz)) + [0.0, C0[5]], 0.5, bz, K, SK)
    k1 = _kin(o, bz)
    assert max(abs(r[0] - k1[0]), abs(r[1] - k1[1]), abs(r[2] - k1[2]), abs(r[3] - k1[3])) < tol


@pytest.mark.parametrize('bz', [0.0, 0.4, 2.0])
def test_sol_region_bend_radiation_loss(v2_mode, bz):
    for defs in ('BEND E1=(L=0.5 K0=0.02);', 'MULT E1=(L=0.5 K0=0.02 DISFRIN=1);'):
        if bz == 0.0 and defs.startswith('MULT'):
            continue
        o, ph = _trk(_elem(defs), np.zeros(6), bz, rad=True)
        exp = -(2.0 / 3.0) * RE * ph.h0**3 * (0.02 / 0.5)**2 * 0.5
        assert abs(o[5] / exp - 1.0) < 2e-4


def test_sol_region_cavity_and_rf_limits(v2_mode):
    c = np.array([1e-3, 2e-3, -1.5e-3, 1e-3, 3e-3, 0.01])
    for defs in ('CAVI E1=(L=0.5 VOLT=1e6 FREQ=5e8 PHI=0.3);', 'CAVI E1=(L=0 VOLT=1e6 FREQ=5e8 PHI=0.3);'):
        el = _elem(defs)
        a, _ = _trk(el, c, 0.0, rfsw=True)             # bz=0: 普通の空洞 (tcav)
        b, _ = _trk(el, c, 1e-9, rfsw=True)            # bz→0 の領域内の式 (tmultiacc の骨組み)
        assert np.abs(a - b).max() < 1e-10
        assert abs(a[5] - c[5]) > 1e-4                 # 加速されている
    # 電圧 0 の空洞は領域内のドリフト (ソレノイドの厳密解) と同じ
    el = _elem('CAVI E1=(L=0.5 VOLT=0 FREQ=5e8);')
    d = _elem('DRIFT E1=(L=0.5);')
    a, ph = _trk(el, c, 0.7, rfsw=True)
    S = tn.make_state(np.array([c], float), ph)
    x.k_sol_drift(S, 1, 0.5, 0.7, False, False, tn.make_rs(1), ph.p0, ph.anrad, ph.cuc)
    assert np.abs(a[:4] - tn.state_out(S, 1)[0][:4]).max() < 1e-14


RING = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) DS=(L=0.3);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
BEND BK = (L=0.2 K0=1e-4 DISFRIN=1) BKM=(L=0.2 K0=-1e-4 DISFRIN=1);
MULT MS = (L=0.2 K1=0.05 K2=1 DISFRIN=1);
CAVI RFC = (VOLT=1.5e6 HARM=64);
CAVI RF2 = (L=0.3 VOLT=2e5 HARM=64 DISFRIN=1);
SOL SB = (BZ=%s BOUND=1 GEO=1) SE = (BOUND=1 GEO=1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE INS = (%s DS %s DS %s);
LINE RING = (ST 8*CELL1 INS RFC);
FFS;
USE RING;
CELL; CALC;
$FORM="20.13";
PageWidth=1999;
Print["NU ",Twiss["NX","$$$"]/2/Pi," ",Twiss["NY","$$$"]/2/Pi," ",Twiss["DX","$$$"]," ",Twiss["DY","$$$"]];
e=Emittance[];
em=Emittances/.e; dr=DampingRate/.e; tu=Tunes/.e;
Print["EM ",em[[1]]," ",em[[3]]," ",EnergyLossU0/.e," ",dr[[1]]," ",dr[[2]]," ",dr[[3]]," ",tu[[1]]," ",tu[[2]]," ",tu[[3]]];
Exit;
'''


def _ring(bz, ins, sol):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(RING % (bz, 'SB' if sol else '', ins, 'SE' if sol else ''), out=buf)
    return {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
            if l.split() and l.split()[0] in ('NU', 'EM')}


@pytest.mark.parametrize('ins', ['BK DS BKM', 'MS', 'RF2'])
def test_sol_region_with_zero_field_equals_ordinary_elements_in_ring(v2_mode, ins):
    a = _ring('1e-9', ins, False)
    b = _ring('1e-9', ins, True)
    for k in ('NU', 'EM'):
        for u, v in zip(a[k], b[k]):
            assert abs(u - v) <= 2e-7 * abs(u) + 1e-12


@pytest.mark.parametrize('ins', ['BK DS BKM', 'MS', 'RF2'])
def test_sol_region_with_field_ring_is_consistent(v2_mode, ins):
    """bz = 3 T: 完走し、Emittance[] の横のチューンが CALC と一致、減衰率の和が 2 U0/E (Robinson)"""
    o = _ring('3', ins, True)
    assert all(math.isfinite(v) for v in o['NU'] + o['EM'])
    assert abs(-sum(o['EM'][3:6]) / (2 * o['EM'][2] / 3.0e9) - 1.0) < 2e-3
    fx = o['NU'][0] - round(o['NU'][0])
    fy = o['NU'][1] - round(o['NU'][1])
    # 結合があるので CALC のモードのチューンと 6D のチューンは同じ組 (小さな差は RF・放射)
    assert min(abs(abs(o['EM'][6]) - abs(fx)), abs(abs(o['EM'][6]) - abs(fy))) < 5e-3
    assert min(abs(abs(o['EM'][7]) - abs(fx)), abs(abs(o['EM'][7]) - abs(fy))) < 5e-3


LINE = r'''MOMENTUM = 3GEV;
DRIFT DS=(L=0.3);
BEND BK = (L=0.2 K0=1e-3) BKM=(L=0.2 K0=-1e-3);
MULT MK = (L=0.2 K0=5e-4 K1=0.2 K2=5 SK1=0.05);
APERT AP = (DX1=-0.01 DX2=0.01 DY1=-0.01 DY2=0.01);
SOL SB = (BZ=2 BOUND=1 GEO=1) SE = (BOUND=1 GEO=1);
MARK ST = ();
LINE BL = (ST SB DS BK DS MK AP DS BKM DS SE DS);
FFS;
USE BL;
INS; CALC;
$FORM="20.14";
PageWidth=1999;
Print["TW ",Twiss["DX","$$$"]," ",Twiss["DPX","$$$"]," ",Twiss["DY","$$$"]," ",Twiss["DPY","$$$"]];
NORAD; LOSSMAP;
r=TrackParticles[{1,{{0,0.012},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1}}},"$$$",1,1];
Print["TR ",r[[2,1,1]]," ",r[[2,2,1]]," ",r[[2,3,1]]," ",r[[2,4,1]]];
Print["LS ",r[[2,7,1]]," ",r[[2,7,2]]," ",r[[2,8,2]]," ",LINE["POSITION","AP"]];
Exit;
'''


def _line():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(LINE, out=buf)
    return {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
            if l.split() and l.split()[0] in ('TW', 'TR', 'LS')}


def test_sol_region_calc_orbit_equals_tracking_and_apert_drops(v2_mode):
    o = _line()
    assert abs(o['TW'][0]) > 1e-5                           # 領域内の蹴りで軌道が動く
    for a, b in zip(o['TW'], o['TR']):
        assert abs(a - b) < 1e-12
    # APERT (領域内): 12 mm の粒子は ±10 mm の絞りで落ち、落ちた要素は APERT
    assert o['LS'][0] == 1.0 and o['LS'][1] == 0.0 and o['LS'][2] == o['LS'][3]


def test_sol_region_bend_is_not_tracked_in_k64_mode(k64_mode):
    with pytest.raises(NotImplementedError):
        _line()


# ---------------------------------------------------------------- MULT のスキュー成分 (pysad の不具合の修正。両モード、参照は gs)
SK = r'''MOMENTUM = 1GEV;
DRIFT D1 = (L=0.5);
MULT M1 = (L=0.5 K1=0.1 SK1=0.2 SK2=30 K2=10 SK0=0.001 DISFRIN=1);
MULT M2 = (L=0 SK1=0.05 SK3=500);
MARK ST = ();
LINE BL = (ST D1 M1 D1 M2 D1);
FFS;
USE BL;
INS; CALC;
$FORM="18.14";
PageWidth=1999;
NORAD; NORFSW;
r=TrackParticles[{1,{{1e-3},{2e-3},{-1.5e-3},{1e-3},{0},{0.01},{1}}},6,1,1];
Print["TR ",r[[2,1,1]]," ",r[[2,2,1]]," ",r[[2,3,1]]," ",r[[2,4,1]]];
Print["TW ",Twiss["BX",6]," ",Twiss["BY",6]," ",Twiss["R1",6]," ",Twiss["R2",6]," ",Twiss["DX",6]," ",Twiss["DY",6]];
Exit;
'''
SK_TR = [.00363483140233, .00161126848062, .00107981471487, .00267838027652]
SK_TW = [2.98217220173114, 3.59817004483586, -.15058637850085, .11736876681750, 1.045384629890E-5, 7.551615120748E-4]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_mult_skew_components_match_gs(mode):
    """MULT の SK0 SK1 SK2 SK3 (ANGLE=0): 以前は黙って無視していた (lattice._mult_params)。参照は gs (2026-11-17)。
    k64 は gs と 2e-13 で一致。v2 の本体の積分は 4 次の合成 (B226) で、K64 の Strang 分割 (2 次) の誤りの分だけ gs と違う:
    この MULT (L = 0.5、K1 = 0.1、軌道の傾き 1〜3 mrad) では TR の差が 1.1e-11 以下 (EPS を 1 → 0.1 → 0.01 と小さくすると差が
    1.1e-11 → 1.5e-12 → 1.6e-13 と刻みの 2 乗に比例して減る = K64 側の積分誤り)。v2 で quadbody.set_scheme(0) にすると gs と
    2e-13 で一致する (差の原因が本体の積分の方式だけであることの確認)。TW (光学) も同じで、k64・v2 (K64 の方式) は 1e-9、v2 の既定は
    BX・BY で 2.5e-9・1.2e-8 (相対 3.3e-9。EPS に比例して減る) の差"""
    from pysad.ffscmd import run_sad_string
    from pysad import quadbody
    old = compat._mode
    old_sch = quadbody.set_scheme(0 if mode == 'v2' else None)       # v2 で本体を K64 の方式に固定した場合
    results = {}
    try:
        for tag, sch in (('fixed', 0 if mode == 'v2' else None), ('default', None)):
            quadbody.set_scheme(sch)
            compat.set_mode(mode)
            buf = io.StringIO()
            run_sad_string(SK, out=buf)
            results[tag] = {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
                            if l.split() and l.split()[0] in ('TW', 'TR')}
    finally:
        compat.set_mode(old)
        quadbody.set_scheme(old_sch)
    for tag, o in results.items():
        tol = 2e-13 if (mode == 'k64' or tag == 'fixed') else 2e-11      # v2 の既定: K64 の Strang の誤り (実測 1.1e-11)
        for a, b in zip(o['TR'], SK_TR):
            assert abs(a - b) < tol
        tolw = 1e-9 if tol == 2e-13 else 1e-8                            # v2 の既定: CALC の本体 (B226) の差 (BY で相対 3.3e-9)
        for a, b in zip(o['TW'], SK_TW):
            assert abs(a - b) < tolw * max(abs(b), 1.0)
    if mode == 'v2':                                                  # 4 次の合成が実際に効いている (K64 の方式と違う)
        assert max(abs(a - b) for a, b in zip(results['default']['TR'], SK_TR)) > 1e-12
