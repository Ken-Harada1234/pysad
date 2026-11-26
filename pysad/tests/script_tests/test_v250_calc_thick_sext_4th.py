# -*- coding: utf-8 -*-
"""v2.5.0: CALC・Emittance[] の厚い多極 (SEXT・OCT・DECA・DODECA、L ≠ 0) の積分を追跡 (tthin) と同じ 4 次の力勾配型の分割にする (docs/30 B326)。

K64 の CALC の厚い多極 (tthine) は L/6・キック/2・2L/3・キック/2・L/6 (2 次、補正キックなし)。追跡の tthin は L(1/2−α)・キック/2・αL・補正キック・αL・キック/2・L(1/2−α)
(4 次、α = 1/√12)。このため動いた軌道のまわりの線形化 (行列) が追跡のヤコビアンとずれ、3 GeV の 8 セルリング (K2 = 2、L = 0.2 m の SF) の dν/dδ が水平で 4.4e-4 (相対)
(L² に比例) ずれた。v2 は CALC の COD の写像・解析ヤコビアン (連鎖律の積) を追跡と同じ 4 次の分割にする (lattice._tthine_cod の v2 の枝、emit6d._jac_thin4、
elkern._tthine_map・_jac_thin4)。k64 は従来のまま (ビット単位で同じ)、L = 0 は変えない。

基準 (pysad のコードを使わない): 運動学的な運動方程式 x' = px/pz、y' = py/pz、z' = 1 − pr/pz、px' = −(K/L) Re (x − iy)^n/n!、py' = −(K/L) Im (x − iy)^n/n! とその変分方程式
(6×6 のヤコビアン) を solve_ivp (DOP853、rtol 1e-13) で積分する。δ は一定、運動量は正準 (1+δ で規格化しない)、dv 項は含めない。
限界: 試験の条件は 1 つの軌道 (x = 3 mm、y = −2 mm、px = 1 mrad、py = 0.5 mrad、δ = 2e-3) と L = 0.1〜0.4 m、強さ K2 = 2〜50 など。
ずれ・ロール (tchge) と非線形フリンジ (ttfrin) は 4 次の本体の前後の別の段なので、ここでは対称性と numba 版・python 版の一致だけを見る。
MULT (tmulte の tmulke の分割) と BEND の 6 極成分は対象外 (従来のまま)。放射・スピンの経路 (tradke) は本体の写像と別 (本体だけが変わる)。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, elements as E, emit6d as e, lattice as lat, track6
from pysad.ffscmd import run_sad_string
from . import test_v250_chromaticity_bend as CB

pytestmark = pytest.mark.filterwarnings('ignore')

C0 = np.array([3e-3, 1e-3, -2e-3, 5e-4, 0.0, 2e-3])
S6 = np.kron(np.eye(3), np.array([[0.0, 1.0], [-1.0, 0.0]]))
# (kord, 強さ K [/m^kord、ak = K·L]) の組。L = 0.2 で誤差が見える程度
CASES = [(2, 2.0), (2, 50.0), (3, 500.0), (4, 1e4), (5, 1e5)]
CLASSES = {2: (E.Sextupole, 'K2'), 3: (E.Octupole, 'K3'), 4: (E.Decapole, 'K4'), 5: (E.Dodecapole, 'K5')}


def _mode(m):
    old = compat._mode
    compat.set_mode(m)
    return old


def ref(c, L, ak, n):
    """独立な基準: 厚い多極の厳密な運動方程式の積分。戻り (COD、6×6 ヤコビアン)。dv 項なし"""
    aki = ak / math.factorial(n)

    def rhs(s, v):
        x, px, y, py, z, d = v[:6]
        M = v[6:].reshape(6, 6)
        pr = 1 + d
        pz = math.sqrt(pr * pr - px * px - py * py)
        a = px * px + py * py
        w = complex(x, -y)
        cx = w ** n
        g = n * w ** (n - 1)
        f = np.array([px / pz, -aki / L * cx.real, py / pz, -aki / L * cx.imag, 1 - pr / pz, 0.0])
        A = np.zeros((6, 6))
        A[0, 1] = 1 / pz + px * px / pz ** 3
        A[0, 3] = px * py / pz ** 3
        A[0, 5] = -px * pr / pz ** 3
        A[2, 1] = A[0, 3]
        A[2, 3] = 1 / pz + py * py / pz ** 3
        A[2, 5] = -py * pr / pz ** 3
        A[4, 1] = -pr * px / pz ** 3
        A[4, 3] = -pr * py / pz ** 3
        A[4, 5] = a / pz ** 3
        A[1, 0] = -aki / L * g.real
        A[1, 2] = -aki / L * g.imag
        A[3, 0] = -aki / L * g.imag
        A[3, 2] = aki / L * g.real
        return np.concatenate([f, (A @ M).ravel()])
    s = solve_ivp(rhs, (0, L), np.concatenate([c, np.eye(6).ravel()]), method='DOP853', rtol=1e-13, atol=1e-16)
    return s.y[:6, -1], s.y[6:, -1].reshape(6, 6)


def _err(m, c, L, ak, n):
    old = _mode(m)
    try:
        xe, Je = ref(c, L, ak, n)
        xc = lat._tthine_cod(c, L, ak, n)
        Jc = e._jac_thin(c, L, ak, n, 0.0)
    finally:
        compat.set_mode(old)
    return np.abs(xc - xe).max(), np.abs(Jc - Je).max()


# (a) 厳密な積分との一致と、L の冪 (4 次: 1 刻みの誤差 L^5、2 次: L^3)
@pytest.mark.parametrize('n,K', CASES)
def test_v2_matches_exact_integration_and_error_scales_L5(n, K):
    errs = {m: {L: _err(m, C0, L, K * L, n) for L in (0.4, 0.2, 0.1)} for m in ('v2', 'k64')}
    for L in (0.4, 0.2, 0.1):
        assert errs['v2'][L][1] < errs['k64'][L][1] / 20.0, (n, K, L, errs)       # 観測: 260〜17000 倍小さい
    # 誤差の L の冪: v2 は L を半分にすると 1/32 (16 以上 64 以下)、k64 は 1/8 (6〜10)
    r2 = errs['v2'][0.2][1] / errs['v2'][0.1][1]
    r64 = errs['k64'][0.2][1] / errs['k64'][0.1][1]
    assert 16.0 < r2 < 64.0, (n, K, r2)
    assert 6.0 < r64 < 10.0, (n, K, r64)
    assert 16.0 < errs['v2'][0.4][1] / errs['v2'][0.2][1] < 64.0


def test_v2_cod_error_scales_L5_sextupole():
    e1 = {L: _err('v2', C0, L, 50.0 * L, 2)[0] for L in (0.4, 0.2, 0.1)}
    assert 16.0 < e1[0.4] / e1[0.2] < 64.0 and 16.0 < e1[0.2] / e1[0.1] < 64.0, e1
    k1 = {L: _err('k64', C0, L, 50.0 * L, 2)[0] for L in (0.4, 0.2, 0.1)}
    assert e1[0.2] < k1[0.2] / 100.0


# 追跡 (tthin) の 1 刻みの写像との一致 (刻みの数: tthin も 1 刻み。CALC の積分の構造が同じ)
@pytest.mark.parametrize('n', [2, 3, 4, 5])
@pytest.mark.parametrize('L', [0.1, 0.4])
def test_v2_cod_equals_tracking_single_step(n, L):
    cls, kn = CLASSES[n]
    K = {2: 40.0, 3: 900.0, 4: 2e4, 5: 4e5}[n]
    el = cls(L=L, **{kn: K})
    el.DISFRIN = 1.0                 # ttfrin (フリンジ) は別の段
    ph = e._Phys(3e9, 100.0)
    c = np.array([1.2e-3, 2e-4, -8e-4, 1e-4, 0.0, 2.5e-3])
    old = _mode('v2')
    try:
        p = track6._Particle(c, ph)
        track6._tthin_sext(p, el, False)
        out_t = p.out()
        out_c = lat._tthine_cod(c, L, K * L, n)
        out_c[4] -= ph.dv(c[5]) * L
        assert np.abs(out_t - out_c).max() < 1e-14
        # numba 版 (elkern.k_thin) も同じ (dv 項は k_thin の use_dv=False では入らない)
        r = lat._thin_step_fast(L, K * L, n, 0.0, 0.0, 0.0, c, False)
        assert np.abs(r[0] - (out_c + np.array([0, 0, 0, 0, ph.dv(c[5]) * L, 0]))).max() < 1e-14
    finally:
        compat.set_mode(old)
    old = _mode('k64')
    try:
        o2 = lat._tthine_cod(c, L, K * L, n)
        o2[4] -= ph.dv(c[5]) * L
    finally:
        compat.set_mode(old)
    assert np.abs(out_t - o2).max() > 1e-16


# (b) シンプレクティック性。numba 版・python 版の一致
@pytest.mark.parametrize('n,K', CASES)
def test_v2_jacobian_symplectic_and_numba_equals_python(n, K):
    rng = np.random.default_rng(n)
    old = _mode('v2')
    try:
        for _ in range(3):
            c = C0 * rng.uniform(0.3, 1.7, 6)
            c[4] = 0.0
            for L in (0.15, 0.4):
                ak = K * L
                J = e._jac_thin(c, L, ak, n, 0.7)          # dv 項 (J[4,5] += 0.7 L) も入れる
                assert np.abs(J.T @ S6 @ J - S6).max() < 1e-12
                # 数値微分 (COD の写像) との一致
                f = lambda cc: lat._tthine_cod(cc, L, ak, n)
                Jn = np.zeros((6, 6))
                h = 1e-7
                for j in range(6):
                    d = np.zeros(6)
                    d[j] = h
                    Jn[:, j] = (f(c + d) - f(c - d)) / (2 * h)
                Jn[4, 5] += 0.7 * L
                assert np.abs(J - Jn).max() < 1e-8
                # numba (フリンジ付き・なし)
                for fr in (False, True):
                    r = lat._thin_step_fast(L, ak, n, 0.0, 0.0, 0.0, c, fr)
                    assert r is not None
                    assert np.abs(r[1].T @ S6 @ r[1] - S6).max() < 1e-12, (n, L, fr)
                    # python 版 (_thin_cod_matrix) と一致
                    cod_p, M_p = lat._thin_cod_matrix(L, ak, n, 0.0, 0.0, 0.0, c, fr)
                    assert np.abs(r[1] - M_p).max() < 1e-12, (n, L, fr)
                    assert np.abs(r[0] - cod_p).max() < 1e-10, (n, L, fr)       # dv 項の有無 (z だけ) の差を除く
    finally:
        compat.set_mode(old)


def test_v2_rotated_thick_sextupole_matrix_symplectic():
    """ロール (tchge) つきの厚い 6 極の行列 (python 版・numba 版) もシンプレクティック"""
    old = _mode('v2')
    try:
        c = C0.copy()
        L, ak = 0.3, 6.0
        cod_p, M_p = lat._thin_cod_matrix(L, ak, 2, 0.0, 0.0, 0.3, c, True)
        r = lat._thin_step_fast(L, ak, 2, 0.0, 0.0, 0.3, c, True)
        assert np.abs(M_p.T @ S6 @ M_p - S6).max() < 1e-12
        assert np.abs(r[1] - M_p).max() < 1e-12
    finally:
        compat.set_mode(old)


# (c) k64 のビット同一: 変更前 (ブランチの親) で計算した値 (float.hex)
K64_COD2 = ['0x1.4a0fbb7e1bb03p-10', '0x1.982a7de82c177p-13', '-0x1.94a5d33f35976p-11', '0x1.72d11c9a5af2bp-14', '-0x1.eb05e433ab334p-28', '0x1.47ae147ae147bp-9']
K64_COD3 = ['0x1.4a4358d34ed2fp-10', '0x1.a3781f369177fp-13', '-0x1.93c06b7e162a8p-11', '0x1.a2d9747b1cebep-14', '-0x1.005dfdd4bb453p-27', '0x1.47ae147ae147bp-9']
K64_CHRO_SEXT = ['0x1.a92b48e15a60fp+1', '-0x1.3ff6b5fe44b7cp+4']                  # _sext_ring(2.0, 0.2) の chromaticity()
K64_CHRO_CASE = ['-0x1.5289fffe3ff45p+4', '-0x1.72b5c5bba0c1bp+0']                  # CASES['bend_k1_edge_sext']
V2_CHRO_THIN = ['0x1.b5477183d0999p+1', '-0x1.45d392782034dp+4']                    # v2、_sext_ring(2.0, 0.0) (薄い 6 極。変更前と同じ)
C_BIT = np.array([1.2e-3, 2e-4, -8e-4, 1e-4, 0.0, 2.5e-3])


def test_k64_bit_identical_to_before():
    old = _mode('k64')
    try:
        assert [float(v).hex() for v in lat._tthine_cod(C_BIT, 0.3, 40.0 * 0.3, 2)] == K64_COD2
        assert [float(v).hex() for v in lat._tthine_cod(C_BIT, 0.3, 900.0 * 0.3, 3)] == K64_COD3
        # python のヤコビアンと numba 版は、k64 では従来の 2 次の積 (4 次の枝に入らない)
        for n, K in ((2, 40.0), (3, 900.0), (4, 2e4), (5, 4e5)):
            Jp = e._jac_thin(C_BIT, 0.3, K * 0.3, n, 0.3)
            r = lat._thin_step_fast(0.3, K * 0.3, n, 0.0, 0.0, 0.0, C_BIT, False)
            assert np.abs(r[1] - e._jac_thin(C_BIT, 0.3, K * 0.3, n, 0.0)).max() < 1e-13
            assert Jp.shape == (6, 6)
        assert [float(v).hex() for v in CB._sext_ring(2.0, 0.2).lat.chromaticity()] == K64_CHRO_SEXT
        assert [float(v).hex() for v in CB._ring('bend_k1_edge_sext').lat.chromaticity()] == K64_CHRO_CASE
    finally:
        compat.set_mode(old)


def _old_2nd_order(c, L, ak, n):
    """変更前の 2 次の分割 (独立な書き直し): ドリフト L/6、キック aki/2、ドリフト 2L/3、キック aki/2、ドリフト L/6"""
    aki = ak / math.factorial(n)

    def kick(c):
        c = c.copy()
        cx = complex(c[0], -c[2]) ** n
        c[1] -= aki * 0.5 * cx.real
        c[3] -= aki * 0.5 * cx.imag
        return c
    c = lat._exact_drift_cod(c, L / 6.0)
    c = kick(c)
    c = lat._exact_drift_cod(c, L / 1.5)
    c = kick(c)
    return lat._exact_drift_cod(c, L / 6.0)


@pytest.mark.parametrize('n', [2, 3, 4, 5])
def test_k64_is_second_order_split(n):
    old = _mode('k64')
    try:
        K = {2: 40.0, 3: 900.0, 4: 2e4, 5: 4e5}[n]
        assert np.array_equal(lat._tthine_cod(C_BIT, 0.3, K * 0.3, n), _old_2nd_order(C_BIT, 0.3, K * 0.3, n))
    finally:
        compat.set_mode(old)


# (d) L = 0 (薄い多極) は v2 でも変えない
@pytest.mark.parametrize('n,K', [(2, 40.0), (3, 900.0), (4, 2e4), (5, 4e5)])
def test_thin_L0_unchanged_in_v2(n, K):
    res = {}
    for m in ('v2', 'k64'):
        old = _mode(m)
        try:
            res[m] = (lat._tthine_cod(C_BIT, 0.0, K, n), e._jac_thin(C_BIT, 0.0, K, n, 0.0))
        finally:
            compat.set_mode(old)
    assert np.array_equal(res['v2'][0], res['k64'][0])
    assert np.array_equal(res['v2'][1], res['k64'][1])
    assert not lat._thin4_active(0.0, K, n)


def test_thin_sextupole_ring_chromaticity_unchanged_in_v2():
    old = _mode('v2')
    try:
        assert [float(v).hex() for v in CB._sext_ring(2.0, 0.0).lat.chromaticity()] == V2_CHRO_THIN
    finally:
        compat.set_mode(old)


# 係数が追跡と同じ (独立に作らない)
def test_coefficients_equal_tracking_and_derivation():
    assert (lat._TT_ALPHA, lat._TT_BETA, lat._TT_GAMMA, lat._TT_ALPHA1) == (track6._TT_ALPHA, track6._TT_BETA, track6._TT_GAMMA, track6._TT_ALPHA1)
    from pysad import elkern as ek
    assert (ek._TT_ALPHA, ek._TT_BETA, ek._TT_GAMMA) == (track6._TT_ALPHA, track6._TT_BETA, track6._TT_GAMMA)
    assert abs(lat._TT_ALPHA - 1.0 / math.sqrt(12.0)) < 1e-16
    assert abs(lat._TT_BETA - (1.0 / 6.0 - lat._TT_ALPHA / 2.0)) < 1e-16
    assert abs(lat._TT_GAMMA - (1.0 / 40.0 - 1.0 / 24.0 / math.sqrt(3.0))) < 1e-16


# (a') リング: dν/dδ が追跡と、薄い 6 極に細かく分けたリングの追跡 (独立な厳密解) と 1e-5 (相対) 以内。k64 は離れる
@pytest.mark.parametrize('K2,L', [(2.0, 0.2), (4.0, 0.4), (4.0, 0.2)])
def test_ring_dnu_ddelta_v2_equals_tracking_and_thin_slice_limit(K2, L):
    res = {}
    for m in ('v2', 'k64'):
        old = _mode(m)
        try:
            res[m] = np.array(CB._sext_ring(K2, L).lat.chromaticity())
        finally:
            compat.set_mode(old)
    old = _mode('v2')
    try:
        trk = CB._trk(CB._sext_ring(K2, L))
        trk_slices = CB._trk(CB._sext_ring(K2, L, nslice=40))     # 追跡の Jacobian (CALC の薄い 6 極の差は含まない)
    finally:
        compat.set_mode(old)
    assert np.allclose(res['v2'], trk, rtol=1e-5, atol=1e-5)
    assert np.allclose(res['v2'], trk_slices, rtol=1e-5, atol=2e-5)
    assert abs(res['k64'][0] / trk[0] - 1) > 1e-4            # k64: 水平で 4.4e-4 (K2 = 2、L = 0.2)〜1.6e-3
    assert np.abs(res['v2'] - trk).max() < np.abs(res['k64'] - trk).max() / 20.0


def test_ring_thick_sextupole_v2_differs_from_k64_by_second_order_error():
    """CALC (FFS の run_sad_string) で厚い 6 極 (L = 0.4、K2 = 4) の dν/dδ が v2 と k64 で 1e-3 以上違う (k64 の 2 次の分割の誤差。水平 0.0115)"""
    src = ('MOMENTUM = 3GEV;\nDRIFT D1=(L=0.8) D2=(L=0.3);\nBEND B=(L=1.2 ANGLE=%.15g);\n'
           'QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);\nSEXT SF=(L=0.4 K2=4);\n'
           'MARK ST=();\nLINE CELL=(QF D1 B D2 SF QD D1 B D2);\nLINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);\n'
           'FFS;\nUSE RING;\nCELL;\nCALC;\n' % CB.ANG)
    vals = {}
    for m in ('v2', 'k64'):
        old = _mode(m)
        try:
            buf = io.StringIO()
            ffs, _ = run_sad_string(src, out=buf)
            vals[m] = np.array(ffs.lat.chromaticity())
            assert np.all(np.isfinite(vals[m]))
        finally:
            compat.set_mode(old)
    assert np.abs(vals['v2'] - vals['k64']).max() > 1e-3        # 4 次への変更で dν/dδ が動く (L = 0.4、K2 = 4: 水平 0.0115)
