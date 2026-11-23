# -*- coding: utf-8 -*-
"""v2.4.0 D2: 硬い端の極 (ANGLE=0 の BEND、K0・AE) の端の幾何 (B220・B221) と、CALC の分散 (B222: 誤りではないと確認)。

独立な基準は `_d2_indep.py` (硬い端の面を実験室の直線座標で厳密に積分する写像)。
- B220: 放射の「場の中の行程」を、面までの直進の行程 ℓ (px の 1 次・p/pz の因子を含む厳密な式) で引く (v2)。K64 は tan(e)·x だけ。
- B221: 追跡 tsteer の端 (v2) を CALC (tsteee) と同じ厳密な硬い端の鎖にする。K64 の追跡は x tan(e)/ρ の薄い蹴りだけ。
- B222: K64 の CALC の分散 (EX・EPX) と行列の δ 列は、追跡の閉軌道の δ 微分・追跡のヤコビアンと一致する (誤りではない。η' の発生 K0 は正準運動量では
  後続のドリフトの T[x, δ] = −L px の中にある)。
"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat, lattice as lat
from pysad.perfkern import edge_free, edge_len
from . import _d2_indep as I
from . import _p3_indep as X

pytestmark = pytest.mark.filterwarnings('ignore')


def _mode(m):
    old = compat._mode
    compat.set_mode(m)
    return old


@pytest.mark.parametrize('sgn', [-1.0, 1.0])
def test_edge_free_gradient_and_linear_limit(sgn):
    args = (0.3, 0.1, 0.004, -0.002, 0.03, 0.01, 0.002)
    ell, g = edge_free(*args, sgn)
    names = [0, 1, 2, 3, 4, 5]
    idx_map = {0: 2, 1: 4, 2: 3, 3: 5, 5: 6}      # 座標 (x, px, y, py, δ) -> 引数の位置
    for k, pos in idx_map.items():
        a = list(args)
        a[pos] += 1e-7
        b = list(args)
        b[pos] -= 1e-7
        num = (edge_len(*a, sgn) - edge_len(*b, sgn)) / 2e-7
        assert abs(num - g[k]) < 1e-7 * max(1.0, abs(num)), (k, num, g[k])
    assert g[4] == 0.0
    # tan e -> 0: 線形の式 tan e · x
    l0, _ = edge_free(1e-9, 0.0, 0.004, 0.0, 0.03, 0.0, 0.0, sgn)
    assert abs(l0 / (1e-9 * 0.004) - 1.0) < 1e-3


@pytest.mark.parametrize('K0,ae1,ae2,c', [(0.06, 0.3, 0.3, [0.01, 0.03, 0, 0, 0, 0.004]), (-0.12, 0.3, 0.2, [-0.02, -0.06, 0, 0, 0, 0]),
                                          (0.006, 0.3, 0.3, [0.001, 0.003, 0, 0, 0, 0.01])])
@pytest.mark.parametrize('fast', [False, True])
def test_steer_tracking_edge_is_exact_hard_edge_v2_and_not_k64(K0, ae1, ae2, c, fast):
    """v2: 追跡 tsteer の 1 要素の写像が硬い端の厳密な写像と 1e-13 以内 (python 版・numba 版とも)。k64 は K64 のまま (薄い蹴り: K0=0.06 で 1e-5)"""
    el = 'BEND SA = (L=0.2 ANGLE=0 K0=%g AE1=%g AE2=%g FRINGE=0);' % (K0, ae1, ae2)
    r = I.pole_exact(c[0], c[1], c[5], 0.2, K0, ae1, ae2)
    old = _mode('v2')
    try:
        o = X.one_element_track(el, c, False, fast)
        assert abs(o[0] - r[0]) < 1e-13 and abs(o[1] - r[1]) < 1e-13, (o[:2], r[:2])
        compat.set_mode('k64')
        o = X.one_element_track(el, c, False, fast)
        if K0 == 0.06:
            assert abs(o[1] - r[1]) > 2e-6
    finally:
        compat.set_mode(old)


def test_calc_steer_map_is_exact_hard_edge_both_modes():
    """CALC (tsteee の端の鎖) の写像は k64・v2 とも硬い端の厳密な写像と一致 (x・px・z で 1e-12)"""
    from types import SimpleNamespace as NS
    for mode in ('k64', 'v2'):
        old = _mode(mode)
        try:
            for K0, ae, c in [(0.06, 0.3, [0.02, 0.06, 0, 0, 0, 0.01]), (-0.12, 0.3, [0.01, -0.03, 0, 0, 0, 0])]:
                el = NS(L=0.2, K0=K0, E1=ae, E2=ae, DX=0.0, DY=0.0)
                out, M = lat._tsteee_cod_matrix(el, np.array(c, float))
                r = I.pole_exact(c[0], c[1], c[5], 0.2, K0, ae, ae)
                assert abs(out[0] - r[0]) < 1e-12 and abs(out[1] - r[1]) < 1e-12 and abs(out[4] + r[4]) < 1e-12
        finally:
            compat.set_mode(old)


@pytest.mark.parametrize('ae', [0.1, 0.3])
def test_emit_radiation_jacobian_matches_exact_pole_v2(ae):
    """v2: 1 つの極の放射つき写像 (x・px・δ・z) とそのヤコビアン (T + D、δ の行を含む) が、硬い端の厳密な放射つき写像 (RK4) の数値微分と一致する。
    δ の行 (dδ/dx・dδ/dpx) の誤差は線形の行程の式だと 1% (8e-9)、厳密な行程で 3e-11 以下 (基準の RK4・数値微分の雑音)"""
    from types import SimpleNamespace as NS
    from pysad import emit6d
    from pysad.emit6d import _State, _Phys, _tsteee
    old = _mode('v2')
    old_div = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [20.0] * 3
    try:
        E = 1e9
        L, K0 = 0.2, 0.06
        c = [-0.01, 0.03, 0, 0, 0, 0.0]

        def calc(cc):
            st = _State(np.array([cc[0], cc[1], 0, 0, cc[4], cc[5]], float), _Phys(E, 100.0))
            st.radcod = True
            _tsteee(st, NS(L=L, K0=K0, E1=ae, E2=ae, DX=0.0, DY=0.0, F1=0, FB1=0, FB2=0, FRINGE=0, DISFRIN=0), True)
            return st

        def exact(cc):
            y = I.pole_radiative([cc[0], cc[1], cc[5], -cc[4]], L, K0, ae, ae, E)
            return np.array([y[0], y[1], y[2], -y[3]])
        st = calc(c)
        ex = exact(c)
        assert abs(st.cod[0] - ex[0]) < 1e-11 and abs(st.cod[1] - ex[1]) < 1e-11
        assert abs(st.cod[5] - ex[2]) < 5e-12                  # 放射による δ (2.5e-7)
        idx = [0, 1, 5, 4]
        Jc = (st.T + st.D)[np.ix_(idx, idx)]
        Je = I.numeric_jacobian(lambda v: exact(v), c, [0, 1, 5, 4])
        # exact の戻りの並びは (x, px, δ, z)、calc の idx は (x, px, δ, z) に合わせた
        assert abs(Jc[2, 0] - Je[2, 0]) < 3e-12 and abs(Jc[2, 1] - Je[2, 1]) < 3e-11, (Jc[2], Je[2])
        assert np.abs(Jc[:2, :3] - Je[:2, :3]).max() < 1e-8
    finally:
        emit6d.EMITDIV[:] = old_div
        compat.set_mode(old)


_FODO = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) DS = (L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND SA = (L=0.2 ANGLE=0 K0=%g AE1=%g AE2=%g);
BEND SB = (L=0.2 ANGLE=0 K0=%g AE1=%g AE2=%g);
MARK ST = ();
LINE CELLA = (QF D2 B D2 QD DS SA D1);
LINE CELLB = (QF D2 B D2 QD DS SB D1);
LINE CC = (CELLA CELLB);
LINE RING = (ST 6*CC);
FFS;
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
Exit;
'''


@pytest.mark.parametrize('K0', [0.0015, 0.003, 0.006])
def test_steer_tracking_tune_equals_calc_v2_ae_part_exact(K0):
    """FODO + ステアリング 2 台 (AE1 = AE2 = 0.3): 追跡の一周写像の固有値と CALC のチューンの差のうち AE に依る分が v2 で 2e-7 以下 (以前 AE = 0.3 で
    6e-6・2.5e-5・1e-4 (K0 = 0.0015・0.003・0.006)。追跡の薄い蹴りの K0² tan²e のヤコビアンの誤り)。AE に依らない K0² の差 (K0 = 0.006 で 2.3e-6、
    0.012 で 9e-6。原因は AE と無関係で、この作業では未調査) は残る。k64 は K64 の差のまま (>1e-5)"""
    old = _mode('v2')
    try:
        d = {}
        for ae in (0.0, 0.3):
            nu, calc = X.tracking_and_calc_tunes(_FODO % (K0, ae, ae, -K0, ae, ae))
            d[ae] = np.array(nu) - np.array(calc)
        assert np.abs(d[0.3] - d[0.0]).max() < 2e-7, d
        assert np.abs(d[0.3]).max() < 5e-6
        compat.set_mode('k64')
        nu, calc = X.tracking_and_calc_tunes(_FODO % (K0, 0.3, 0.3, -K0, 0.3, 0.3))
        assert max(abs(np.array(nu) - np.array(calc))) > 1e-5
    finally:
        compat.set_mode(old)


_WIG = '''MOMENTUM = 1GEV;
DRIFT D1=(L=0.6) DW=(L=0.3);
QUAD QF=(L=0.3 K1=0.7 %(fr)s) QD=(L=0.3 K1=-0.7 %(fr)s);
BEND B=(L=1.0 ANGLE=0.39269908169872414 E1=0 E2=0);
BEND P1=(L=0.2 ANGLE=0 K0=0.06 AE1=%(ae)g AE2=%(ae)g) P2=(L=0.4 ANGLE=0 K0=-0.12 AE1=%(ae)g AE2=%(ae)g);
LINE CELL=(QF D1 B D1 QD D1);
LINE WIG=(DW P1 DW P2 DW P1 DW);
CAVI RF=(L=0 VOLT=2E6 HARM=40);
LINE RING=(RF WIG 16*CELL);
FFS;
USE RING;
CELL; CALC; %(sw)s
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi)];
Exit;
'''


@pytest.mark.parametrize('ae,tol_damp', [(0.1, 3e-5), (0.3, 3e-5)])
def test_wiggler_ring_damping_and_tune_match_exact_v2(ae, tol_damp):
    """16 セル (1 GeV、周長 32 m 余り) + 3 極 (K0 = 0.06・−0.12・0.06 rad) の挿入光源、AE1 = AE2 = AE: v2 の Emittance[] の x の減衰 (ξ = 2 p0 dr/(2 E))
    が、硬い端の厳密な 4 次元写像 (一周写像の固有値の絶対値) と相対 3e-5 以内、チューンは 1e-5 以内 (AE = 0 でも 4e-6 残る: 基準の刻みの雑音) (四極のフリンジ tqfrie を切った比較。K64・pysad の
    四極にはハードエッジの非線形フリンジが既定で入り、AE があると軌道が四極の中心から外れるので、独立な基準にはそれが無い。
    以前の「AE² に比例する +1.0%・2.0% (B143 の残り)」は、(1) 四極のフリンジの違い (チューン 8e-5・9e-4)、(2) 基準の U0 が RF なしの軌道上の値で
    閉軌道のエネルギーの沈み込み (δ = −8.8e-4) の分 0.2% 違っていた、の 2 つで、端の模型の誤りではなかった)"""
    import os
    from pysad import ffscmd, emit6d
    from .rad_indep_d2 import Ring4X
    nb = 16
    th = 2 * math.pi / nb
    old = _mode('v2')
    old_div = list(emit6d.EMITDIV)
    try:
        src = _WIG % dict(fr='DISFRIN=1', ae=ae, sw='NORFSW; NORAD; NOFLUC;')
        buf = io.StringIO()
        ffs, _ = ffscmd.run_sad_string(src, out=buf)
        calc_nu = float(re.search(r'CALC tunes\s+([-\d.eE+]+)', buf.getvalue()).group(1)) % 1.0
        calc_nu = min(calc_nu, 1 - calc_nu)
        ffs, _ = ffscmd.run_sad_string(src.replace('NORFSW; NORAD; NOFLUC;', 'RFSW; RAD; NOFLUC;'), out=io.StringIO())
        emit6d.EMITDIV[:] = [20.0] * 3
        em = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
        dr = np.abs(em['damping_rate'])
        E = 1e9
        ring = Ring4X(E, ae, th, nb)
        ev = ring.solve()
        ph = np.angle(ev) / 2 / math.pi
        xm = [z for z, f in zip(ev, ph) if f > 0.15][0]
        dec_exact = 1 - abs(xm)
        dec_v2 = ffs.momentum() * dr[0] / E            # 一周あたりの x の減衰 (2 p dr/2/E)
        assert abs(dec_v2 / dec_exact - 1.0) < tol_damp, (dec_v2, dec_exact)
        nu_ex = [f for f in ph if f > 0.15][0]
        assert abs(em['tunes'][0] - nu_ex) < 1e-5, (em["tunes"][0], nu_ex)
    finally:
        emit6d.EMITDIV[:] = old_div
        compat.set_mode(old)


def test_calc_dispersion_equals_closed_orbit_derivative_and_delta_column_equals_tracking():
    """B222: K64 の CALC の分散 (EX・EPX) は、極のある (K0 = 0.06・−0.12・0.06 rad、AE = 0.2) リングの全要素で閉軌道の δ 微分と一致する (1e-8)。
    K0 の極の η' の発生 (K0 δ) は、正準運動量では極の後の px_cod = ∓K0 に比例するドリフトの T[x, δ] = −L px の中に入っている。
    行列の δ 列も追跡のヤコビアンと一致する (ステアリング・L=0 の薄いキック・MULT の K0・ANGLE≠0 + K0 の BEND)"""
    from pysad import ffscmd
    for mode in ('k64', 'v2'):
        old = _mode(mode)
        try:
            src = _WIG % dict(fr='', ae=0.2, sw='NORFSW; NORAD; NOFLUC;')
            src = src.replace('CAVI RF=(L=0 VOLT=2E6 HARM=40);', '').replace('LINE RING=(RF WIG', 'LINE RING=(WIG')
            ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
            els = ffs.lat.elements
            rows = ffs.rows if hasattr(ffs, 'rows') else ffs._last_rows
            h = 1e-6
            _, _, ep = lat._qcod_periodic(els, h)
            _, _, em = lat._qcod_periodic(els, -h)
            for k in range(0, 9):
                dcod = (ep[k][1] - em[k][1]) / (2 * h)
                assert abs(dcod[0] - rows[k].etx) < 1e-8 and abs(dcod[1] - rows[k].etpx) < 1e-8, (k, dcod[:2], rows[k].etx, rows[k].etpx)
            # δ 列
            for eld, c0 in [('BEND SA = (L=0.2 ANGLE=0 K0=0.06);', [0.001, 0.02, 0, 0, 0, 0]),
                            ('MULT SA = (L=0.2 K0=0.06);', [0.001, 0.02, 0, 0, 0, 0]),
                            ('BEND SA = (L=0.2 ANGLE=0.1 K0=0.06);', [0.001, 0.02, 0, 0, 0, 0]),
                            ('BEND SA = (L=0.0 ANGLE=0 K0=0.002);', [0.001, 0.02, 0, 0, 0, 0])]:
                src = 'MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n' % eld
                f1, _ = X.run(src)
                el = [e for e in f1.lat.elements if e.name == 'SA'][0]
                c0 = np.array(c0, float)
                J = np.zeros(4)
                for i, j in enumerate((0, 1, 2, 3)):
                    cp = c0.copy()
                    cm = c0.copy()
                    cp[5] += h
                    cm[5] -= h
                    J[i] = (X.one_element_track(eld, cp, False, False)[j] - X.one_element_track(eld, cm, False, False)[j]) / (2 * h)
                Mc, _ = lat._el_cod_step(el, c0, 0.0)
                assert np.abs(J - Mc[:4, 5]).max() < 1e-9, (mode, eld, J, Mc[:4, 5])
        finally:
            compat.set_mode(old)


# ---------------------------------------------------------------- B221 の続き: F1・FB1・FB2 つきの tsteer の端
@pytest.mark.parametrize('extra', [' F1=0.04 FRINGE=3', ' FB1=0.03 FB2=0.02 F1=0.04 FRINGE=3'])
@pytest.mark.parametrize('fast', [False, True])
def test_steer_tracking_with_f1_edge_equals_calc_v2(extra, fast):
    """v2: F1・FB1・FB2 (tblfre の断片) のある tsteer も端を CALC (tsteee) と同じ順序 (tblfre → 厳密な端の鎖 → 断片) の厳密な鎖にした。
    追跡と CALC のチューンの差のうち AE に依る分が 2e-7 以下 (K0 = 0.003、AE1 = AE2 = 0.3)。FRINGE = −1・−2 (片側だけ) は K64 の CALC (tsteee) が
    反対側の端を飛ばすのに追跡 (tsteer) は入れるので、AE = 0 でも y のチューンが 2e-5〜3e-4 違う (K64 の不整合、未対応)"""
    old = _mode('v2')
    try:
        d = {}
        for ae in (0.0, 0.3):
            src = _FODO.replace('AE2=%g)', 'AE2=%g' + extra + ')') % (0.003, ae, ae, -0.003, ae, ae)
            nu, calc = X.tracking_and_calc_tunes(src, fast=fast)
            d[ae] = np.array(nu) - np.array(calc)
        assert np.abs(d[0.3] - d[0.0]).max() < 2e-7, d
    finally:
        compat.set_mode(old)


# ---------------------------------------------------------------- B223: MULT の端の面
@pytest.mark.parametrize('L,K0,K1,a1,a2,x,px', [(0.2, 0.06, 0.0, 0.3, 0.3, 0.01, 0.03), (0.2, 0.06, 0.0, 0.3, 0.0, 0.01, 0.03),
                                               (0.2, 0.06, 0.0, 0.0, 0.3, 0.01, 0.03), (0.2, 0.06, 0.05, 0.3, 0.3, 0.01, 0.03),
                                               (0.2, 0.06, 0.3, 0.3, 0.3, 0.01, 0.03)])
def test_mult_edge_is_exact_hard_edge_v2_calc_and_tracking(L, K0, K1, a1, a2, x, px):
    """v2: MULT (K0 のみ・勾配 K1 つき、DISFRIN = 1) の CALC の写像が硬い端の厳密な写像 (K1 = 0 は解析解、K1 ≠ 0 は `_p3_indep.pole_k1` の ODE) と
    K1 = 0 で 1e-15、K1 ≠ 0 で 3e-7 以内 (K1 ≠ 0 の差は本体の積分 (tsolque・tqente) の差で AE = 0 でも同じ)。追跡 (python の tmulti) と CALC は 1e-14 で一致。
    k64 は K64 の薄い蹴り (K0 = 0.06・AE = 0.3 で 1e-5)"""
    from types import SimpleNamespace as NS
    kn = {0: K0}
    if K1:
        kn[1] = K1
    e = NS(L=L, KnL=kn, ROTATION=0.0, DROTATE=0.0, EPS=0.0, DX=0.0, DY=0.0, DISFRIN=1.0, FRINGE=0.0, K0FR=0.0, FB1=0.0, FB2=0.0, AE1=a1, AE2=a2)
    c = [x, px, 0, 0, 0, 0.0]
    xe, pe, _W = X.pole_k1(x, math.asin(px), 0.0, L, K0, K1, a1, a2)
    old = _mode('v2')
    try:
        m = lat._mult_cod_map(e)(np.array(c, float))
        tol = 1e-14 if K1 == 0.0 else 3e-7
        assert abs(m[0] - xe) < tol and abs(m[1] - math.sin(pe)) < tol, (m[:2], xe, math.sin(pe))
        el = 'MULT SA = (L=%g K0=%g%s AE1=%g AE2=%g DISFRIN=1);' % (L, K0, (' K1=%g' % K1) if K1 else '', a1, a2)
        o = X.one_element_track(el, c, False, False)
        assert abs(o[0] - m[0]) < 1e-14 and abs(o[1] - m[1]) < 1e-14, (o[:2], m[:2])
        compat.set_mode('k64')
        m64 = lat._mult_cod_map(e)(np.array(c, float))
        if K1 == 0.0 and a1 == a2 == 0.3:
            assert abs(m64[1] - math.sin(pe)) > 3e-6
    finally:
        compat.set_mode(old)


@pytest.mark.parametrize('ae,tol_damp,tol_nu', [(0.1, 3e-5, 1e-5), (0.3, 2e-4, 1e-5)])
def test_wiggler_ring_with_mult_poles_matches_exact_v2(ae, tol_damp, tol_nu):
    """MULT で作った 3 極 (K0 = 0.06・−0.12・0.06、AE1 = AE2 = AE、四極のフリンジなし) のリング: v2 の減衰・チューンが硬い端の厳密解と一致
    (以前 (K64 の薄い蹴り) は AE = 0.1・0.3 でチューンが 2.9e-4・2.5e-3、減衰が 2.6e-4・3.2e-3 違った)"""
    from pysad import ffscmd, emit6d
    from .rad_indep_d2 import Ring4X
    th = 2 * math.pi / 16
    old = _mode('v2')
    old_div = list(emit6d.EMITDIV)
    try:
        src = _WIG % dict(fr='DISFRIN=1', ae=ae, sw='RFSW; RAD; NOFLUC;')
        src = src.replace('BEND P1=(L=0.2 ANGLE=0 K0', 'MULT P1=(L=0.2 K0').replace('P2=(L=0.4 ANGLE=0 K0', 'P2=(L=0.4 K0')
        ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
        emit6d.EMITDIV[:] = [20.0] * 3
        em = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
        ring = Ring4X(1e9, ae, th, 16)
        ev = ring.solve()
        ph = np.angle(ev) / 2 / math.pi
        xm = [z for z, f in zip(ev, ph) if f > 0.15][0]
        dec_v2 = ffs.momentum() * abs(em['damping_rate'][0]) / 1e9
        assert abs(dec_v2 / (1 - abs(xm)) - 1.0) < tol_damp
        assert abs(em['tunes'][0] - [f for f in ph if f > 0.15][0]) < tol_nu
    finally:
        emit6d.EMITDIV[:] = old_div
        compat.set_mode(old)


# ---------------------------------------------------------------- 追跡と CALC の AE に依らない K0² の差の原因
def test_tracking_calc_quad_jacobian_differs_at_orbit_offset_and_steer_does_not():
    """追跡と CALC のチューンの AE に依らない K0² の差 (K0 = 0.006 で 2.3e-6、0.012 で 9e-6) の原因: 四極の本体の積分器。軌道が x = 3 mm・px = 2 mrad ずれた
    四極 1 台 (DISFRIN = 1、K1 = −0.6) のヤコビアンが、追跡と CALC (どちらも K64 の tsolqu・tsolque) で 5e-7 違う (厳密な写像との差: CALC 1.8e-7、
    追跡 7e-7。v2 は 4 次の合成で CALC 5e-10・追跡 7.9e-9 に減った。B226)。ステアリング (AE = 0・0.3)・扇形の偏向磁石・ドリフトは 1e-9 以下で一致する"""
    old = _mode('v2')
    try:
        def diff(eld):
            f1, _ = X.run('MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n' % eld)
            el = [e for e in f1.lat.elements if e.name == 'SA'][0]
            c0 = np.array([0.003, 0.002, 0, 0, 0, 0.0])
            h = 1e-6
            J = np.zeros((4, 4))
            for j in range(4):
                cp, cm = c0.copy(), c0.copy()
                cp[j] += h
                cm[j] -= h
                J[:, j] = (X.one_element_track(eld, cp, False, False)[:4] - X.one_element_track(eld, cm, False, False)[:4]) / (2 * h)
            M, _ = lat._el_cod_step(el, c0, 0.0)
            return np.abs(J - M[:4, :4]).max()
        compat.set_mode('k64')
        assert diff('QUAD SA=(L=0.3 K1=-0.6 DISFRIN=1);') > 4e-7      # K64 の四極の積分器の誤り (B226) は k64 で再現のまま
        compat.set_mode('v2')
        assert diff('QUAD SA=(L=0.3 K1=-0.6 DISFRIN=1);') < 3e-8      # v2 は 4 次の合成 (追跡 7.1e-7 → 7.9e-9、CALC 1.8e-7 → 5e-10。test_v241_quad_integrator.py)
        for eld in ('BEND SA=(L=0.2 ANGLE=0 K0=0.06 AE1=0.3 AE2=0.3);', 'BEND SA=(L=0.2 ANGLE=0 K0=0.06);',
                    'BEND SA=(L=1.2 ANGLE=0.26179938779914913 E1=0.1308996938995747 E2=0.1308996938995747);', 'DRIFT SA=(L=1);'):
            assert diff(eld) < 1e-9, eld
    finally:
        compat.set_mode(old)


# ---------------------------------------------------------------- B224: MULT の端の縦収束と tradke のヤコビアン
def test_mult_edge_vertical_focusing_matches_bend_chain_v2():
    """v2 (B224): MULT (K0 のみ) の端の写像が縦方向の収束 (Maxwell) も持つ: y・py が ANGLE=0 の BEND (厳密な端の鎖) と 2e-5 以内
    (縦収束を入れる前は 8e-4。K64 の薄い蹴りの縦成分 Δpy = fe·Im(cx) を足している。残りは薄い蹴りと厳密な Maxwell の違い)"""
    from types import SimpleNamespace as NS
    old = _mode('v2')
    try:
        L, K0, ae = 0.2, 0.06, 0.3
        c = np.array([0.01, 0.03, 0.004, 0.002, 0, 0.0])
        cb, _ = lat._tsteee_cod_matrix(NS(L=L, K0=K0, E1=ae, E2=ae, DX=0.0, DY=0.0), c)
        e = NS(L=L, KnL={0: K0}, ROTATION=0.0, DROTATE=0.0, EPS=0.0, DX=0.0, DY=0.0, DISFRIN=1.0, FRINGE=0.0, K0FR=0.0, FB1=0.0, FB2=0.0, AE1=ae, AE2=ae)
        cm = lat._mult_cod_map(e)(c)
        assert np.abs(cm[:4] - cb[:4]).max() < 2e-5, (cm[:4], cb[:4])
    finally:
        compat.set_mode(old)


@pytest.mark.parametrize('K1', [-0.01, 0.01])
def test_mult_edge_radiation_jacobian_equals_numeric_derivative_v2(K1):
    """v2 (B224): K1 ≠ 0 の MULT (K1 < 0 は面の傾きが実験室系の y) の AE のある端 (厳密な硬い端の写像はシンプレクティックでない: K1 = −0.01・K0 = 0.005・AE = 0.3 で
    3e-2) を含む 1 要素の放射つき写像のヤコビアンの δ の行が数値微分と一致する (以前はシンプレクティックの逆行列の式 (tradke_db の tinv・tr2i) のため
    y の微分が 2〜3 倍違い、K1 < 0 の Emittance[] の ΔJx が −0.011 だった。v2 は厳密な逆行列)"""
    from types import SimpleNamespace as NS
    from pysad import emit6d
    from pysad.emit6d import _State, _Phys, _tmulte
    old = _mode('v2')
    old_div = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [20.0] * 3
    try:
        c = [0.002, 0.001, 0.001, 0.0005, 0, 0.0]

        def run(cc):
            st = _State(np.array(cc, float), _Phys(1e9, 100.0))
            st.radcod = True
            el = NS(L=0.2, KnL={0: 0.005, 1: K1}, ROTATION=0.0, DROTATE=0.0, EPS=0.0, DX=0.0, DY=0.0, DISFRIN=1.0, FRINGE=0.0, K0FR=0.0, FB1=0.0, FB2=0.0,
                    AE1=0.3, AE2=0.3, F1=0, F2=0, F1K1F=0, F1K1B=0, F2K1F=0, F2K1B=0, ANGLE=0.0, VOLT=0.0)
            _tmulte(st, el, True)
            return st
        st = run(c)
        idx = [0, 1, 2, 3, 5]
        for j, i in enumerate(idx):
            h = 1e-6
            cp, cm = list(c), list(c)
            cp[i] += h
            cm[i] -= h
            num = (run(cp).cod[5] - run(cm).cod[5]) / (2 * h)
            ana = (st.T + st.D)[5, i]
            assert abs(ana - num) < 2e-11 + 1e-6 * abs(num) if i != 5 else abs(ana - num) < 1e-9, (i, ana, num)
    finally:
        emit6d.EMITDIV[:] = old_div
        compat.set_mode(old)


# ---------------------------------------------------------------- B225: 片側のフリンジ (FRINGE = −1・−2) のステアリング
@pytest.mark.parametrize('extra', [' FRINGE=-1 F1=0.05', ' FRINGE=-2 F1=0.05'])
@pytest.mark.parametrize('fast', [False, True])
def test_steer_one_sided_fringe_tracking_equals_calc_v2(extra, fast):
    """v2 (B225): FRINGE = −1 (入口だけ)・−2 (出口だけ) のステアリング (AE1 = AE2 = 0.3、K0 = 0.003) は、追跡も CALC (tsteee) と同じく反対側の端 (面・縦収束) を飛ばす。
    追跡と CALC のチューンの差が AE の有無によらず四極の積分器の差 (6e-7) だけになる (K64 の追跡は両方の端を入れるので AE = 0 でも y のチューンが
    2e-5・3e-4 違い、AE = 0.3 では 1e-3 違った)"""
    old = _mode('v2')
    try:
        for ae in (0.0, 0.3):
            src = _FODO.replace('AE2=%g)', 'AE2=%g' + extra + ')') % (0.003, ae, ae, -0.003, ae, ae)
            nu, calc = X.tracking_and_calc_tunes(src, fast=fast)
            assert np.abs(np.array(nu) - np.array(calc)).max() < 1e-6, (ae, nu, calc)
    finally:
        compat.set_mode(old)
