# -*- coding: utf-8 -*-
"""v2.4.0 担当 D: B112 の β の因子 (docs/30 B167) と BEND の T[z, δ] (B164 = B139)。

B167 (放射の β の因子): K64 の放射は β = 1 の式。v2 は厳密な放射のスペクトル (Schott) から決めた係数 (pysad/radbeta.py) を掛ける。
  独立な基準: (1) 円軌道の遅延放射場の時間領域の数値積分から求めた第 n 高調波のパワー (`_rad_beta_indep.schott_pn_direct`) が Schott の閉じた式と 1e-10 で一致、
  (2) 全放射パワーが Larmor の (2/3) β³γ⁴ (Σ n S_n)、(3) 古典式 U0 = (2/3) r mc² β³γ⁴ Σθ²/L と Emittance[]・RADINT の U0。
  εx の β の因子 (r2/β⁴) は Emittance[] と RADINT が互いに一致する (どちらも同じ理論の実装。独立な写像との比較は未完: docs/31)。
B164 (T[z, δ]): `tbendebody0` の dω/dp の厳密な式を、写像の数値微分 (4 次の中心差分) と比べる。"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, lattice as L, radbeta

from . import _rad_beta_indep as R

MP = 0.93827208816e9
RP = 1.5346982e-18      # 陽子の古典半径 [m]


# ======================================================================= B164 (B139)
def _jac_fd(c, phi0, rhob, drhob=0.0, h=2e-4):
    J = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = h
        f = lambda d: L._tbendebody0_cod(c + d, phi0, rhob, drhob)
        J[:, j] = (-f(2 * e) + 8 * f(e) - 8 * f(-e) + f(-2 * e)) / (12 * h)
    return J


@pytest.mark.parametrize('c', [np.array([1e-3, 1e-2, 2e-3, 1e-3, 0, 3e-3]), np.array([0, 3e-2, 1e-2, 5e-3, 0, 0.0])])
def test_bendebody0_matrix_exact_v2(c):
    """v2: tbendebody0 の解析行列が写像の数値微分と 1e-10 で一致 (K64 は T[z, δ] ∼ L (px² + py²)、T[y, δ] ∼ L py (px² + py²) の誤差が 1e-4〜1e-3)"""
    phi0 = 0.05
    rhob = 1.2 / phi0
    J = _jac_fd(c, phi0, rhob)
    assert np.abs(L._tbendebody0_mat(c, phi0, rhob, 0.0, True) - J).max() < 1e-10
    err_k64 = np.abs(L._tbendebody0_mat(c, phi0, rhob, 0.0, False) - J).max()
    assert err_k64 > 5e-5
    compat.set_mode('k64')           # 既定の引数は compat に従う: k64 は K64 の式のまま
    try:
        assert np.abs(L._tbendebody0_mat(c, phi0, rhob) - J).max() == pytest.approx(err_k64, rel=1e-12)
    finally:
        compat.set_mode(None)


def test_bendebody0_numba_kernel_matches_python():
    from pysad import elkern
    phi0 = 0.05
    rhob = 1.2 / phi0
    c = np.array([1e-3, 1e-2, 2e-3, 1e-3, 0, 3e-3])
    for flag in (False, True):
        cc = c.copy()
        t = elkern._tbendebody0(cc, phi0, rhob, 0.0, 0.0, flag)
        assert np.abs(t - L._tbendebody0_mat(c, phi0, rhob, 0.0, flag)).max() < 1e-12


# ======================================================================= B167: 厳密なスペクトル
@pytest.mark.parametrize('beta', [0.6, 0.8])
def test_schott_equals_time_domain_radiation(beta):
    """第 n 高調波のパワー: 遅延放射場の時間領域の数値積分 (リエナール・ウィーヘルト、観測方向 ψ で積分) = Schott の閉じた式 (β n S_n)。1e-9 以内"""
    for n in (1, 2, 4):
        d = R.schott_pn_direct(n, beta, nt=1500, npsi=150)
        s = beta * n * radbeta.schott_s(n, beta)
        assert d == pytest.approx(s, rel=1e-9)


def test_schott_total_power_is_larmor():
    """Σ n S_n = (2/3) β³γ⁴ (古典の全放射パワー。β = 0.8、γ = 5/3)"""
    beta = 0.8
    g = 1 / math.sqrt(1 - beta**2)
    tot = sum(n * radbeta.schott_s(n, beta) for n in range(1, 120))
    assert tot == pytest.approx(2.0 / 3.0 * beta**3 * g**4, rel=1e-9)


def test_r2_table_matches_exact_spectrum():
    """r2(γ) の表が Schott の式の直接計算と一致 (γ = 1.5: 離散和、γ = 4: 連続の n の積分。線形補間の誤差 2e-4)"""
    for g in (1.5, 4.0):
        assert radbeta.r2_of_gamma(g) == pytest.approx(radbeta.r2_exact(g), rel=2e-4)
    assert radbeta.r2_of_gamma(1.5) == pytest.approx(0.23774, abs=2e-4)       # 記録: γ = 1.5・2.5・4 で 0.238・0.646・0.850
    assert radbeta.r2_of_gamma(100.0) == pytest.approx(1 - 2.5291 / 1e4, abs=1e-5)


def test_ab_limits_and_derivatives():
    a, b, da, db = radbeta.ab(1e5)
    assert abs(a - 1) < 1e-9 and abs(b - 1) < 1e-9 and abs(da) < 1e-9 and abs(db) < 1e-9
    for g in (1.3, 1.77, 5.0, 50.0):
        h = 1e-7
        a, b, da, db = radbeta.ab(g)
        assert (math.log(radbeta.ab(g + h)[0]) - math.log(radbeta.ab(g - h)[0])) / (2 * h) == pytest.approx(da, rel=1e-5, abs=1e-8)
        assert (math.log(radbeta.ab(g + h)[1]) - math.log(radbeta.ab(g - h)[1])) / (2 * h) == pytest.approx(db, rel=1e-5, abs=1e-8)
        beta2 = 1 - 1 / g**2
        assert a * b == pytest.approx(beta2, rel=1e-12)                              # 平均の損失 β² (δ の損失)
        assert a * b * b == pytest.approx(radbeta.r2_of_gamma(g) / beta2, rel=1e-9)  # N<u²> の係数 r2/β²


# ======================================================================= B167: 陽子のリング
_PC = 3.0e9
_TH = 2 * math.pi / 24
_SRC = '''MOMENTUM = 3.0 GEV; MASS = 0.93827208816 GEV;
DRIFT D1=(L=1.0) D2=(L=0.5);
QUAD QF=(L=0.3 K1=0.7) QD=(L=0.3 K1=-0.6);
BEND B=(L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g);
LINE CELL=(QF D2 B D2 QD D1);
LINE RING=(24*CELL);
FFS;
USE RING; CELL; CALC; RAD; NOFLUC;
Exit;
''' % (_TH, _TH / 2, _TH / 2)


def _ring():
    ffs, _ = ffscmd.run_sad_string(_SRC, out=io.StringIO())
    return ffs


def _emit(ffs, z=1.0):
    return emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False, amass=MP, charge=z)


def _classical_u0():
    g = math.sqrt(1 + (_PC / MP) ** 2)
    b = _PC / MP / g
    rho = 1.2 / _TH
    return (2.0 / 3.0) * RP * MP * b**3 * g**4 * (24 * _TH / rho), b, g


def test_proton_u0_matches_classical_v2():
    """3 GeV/c の陽子 (γ = 3.35、β = 0.954) の 1 周のエネルギー損失 = 古典式 (2/3) r mc² β³γ⁴ Σθ²/L。v2 の Emittance[] は 1e-6 以内、
    β の因子を切ると 1/β³ = 1.15 倍 (超相対論の式)"""
    u0, b, _ = _classical_u0()
    compat.set_mode('v2')
    try:
        ffs = _ring()
        compat.set_radiation_beta(True)
        e1 = _emit(ffs)
        compat.set_radiation_beta(False)
        e0 = _emit(ffs)
    finally:
        compat.set_radiation_beta(None)
        compat.set_mode(None)
    assert e1['energy_loss_eV'] == pytest.approx(u0, rel=2e-6)
    assert e0['energy_loss_eV'] == pytest.approx(u0 / b**3, rel=2e-6)


def test_proton_damping_and_emittance_factors_v2():
    """電荷 z = 1000 (放射を 1e6 倍に強めて数値的に条件のよい状態にする。古典半径 ∝ z²): δ の損失 (減衰率) の β の因子 = β²、εx の因子 = r2/β⁴"""
    u0, b, g = _classical_u0()
    compat.set_mode('v2')
    try:
        ffs = _ring()
        compat.set_radiation_beta(True)
        e1 = _emit(ffs, 1000.0)
        compat.set_radiation_beta(False)
        e0 = _emit(ffs, 1000.0)
    finally:
        compat.set_radiation_beta(None)
        compat.set_mode(None)
    assert e1['damping_rate'][0] / e0['damping_rate'][0] == pytest.approx(b**2, rel=2e-3)
    assert e1['emittance_x'] / e0['emittance_x'] == pytest.approx(radbeta.r2_of_gamma(g) / b**4, rel=3e-3)


def test_proton_radint_factors_v2():
    """RADINT (v2): 陽子の COEFF・CO4 に粒子の古典半径・質量を使い (K64 は電子の値)、β の因子 (U0: β³、εx: r2/β⁴) を掛ける。U0 は古典式と 2e-4 (出力 4 桁)"""
    u0, b, g = _classical_u0()

    def run(flag):
        compat.set_mode('v2')
        compat.set_radiation_beta(flag)
        try:
            buf = io.StringIO()
            ffscmd.run_sad_string(_SRC.replace('Exit;', 'RADINT;\nExit;'), out=buf)
        finally:
            compat.set_radiation_beta(None)
            compat.set_mode(None)
        t = buf.getvalue()
        loss = float(re.search(r'ENERGY LOSS PER TURN \(kV\)\s*=\s*([-0-9.Ee+]+)', t).group(1)) * 1e3
        em = float(re.search(r'EMITTANCE U\s*=\s*([-0-9.Ee+]+)', t).group(1))
        return loss, em
    l1, e1 = run(True)
    l0, e0 = run(False)
    assert l1 == pytest.approx(u0, rel=2e-4)
    assert l0 == pytest.approx(u0 / b**3, rel=2e-4)
    assert e1 / e0 == pytest.approx(radbeta.r2_of_gamma(g) / b**4, rel=1e-3)


def test_k64_radiation_unchanged():
    """k64: β の因子は入らない (anrad は正)"""
    compat.set_mode('k64')
    try:
        ph = emit6d._Phys(_PC, 100.0, amass=MP)
        assert ph.anrad > 0
        e = _emit(_ring())
    finally:
        compat.set_mode(None)
    assert e['energy_loss_eV'] > 0
