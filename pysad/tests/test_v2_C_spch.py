# -*- coding: utf-8 -*-
"""v2 (担当 C): 空間電荷の追跡 — SPAC フラグ + DRIFT の RADIUS、PSPAC フラグ + SPCH 要素 (pysad/spch.py)。

K64 の動作と v2 の根拠は docs/history/v2_changes_C_spch.md。要点:
- K64 の SPAC (spkick.f) は静止系の場に γ が 1 つ足りず、力が 1/γ 倍に弱い (gs で確かめた: 5 MeV、1e11、半径 1 mm の一様な
  ビームの 1 m のドリフトで rms x は 0.62 mm (y は 0.73 mm)。解析式は 1.64 mm、力を 1/γ にした解析式は 0.65 mm)。
- v2 は、静止系の 3 次元のポアソン方程式 (円筒の壁) を解き、正しい係数のキックを入れる。

v2 の値の根拠は解析式と、独立な数値計算 (包絡方程式の数値積分、自前の行列での線形レンズのリング):
  一様な円柱 (ペルビアンス K = 2 r_e λ/(β²γ³))、Bassetti-Erskine、球対称のガウス分布、円筒の中の点電荷 (ベッセル級数)、
  ずれた線電荷の像の力 2λ x0/(R²−x0²)、縦の場の g 因子 1+2 ln(R/a)、KV ビームの包絡方程式、Brillouin 流、Laslett のチューンシフト。
条件: 電子 5 MeV (γ = 9.84)、1e9〜1e11 個、ビーム半径 1 mm 程度、真空ダクトの半径 3 mm〜3 cm。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import erf, j0, j1, jn_zeros, wofz
from scipy.stats import norm, qmc

from pysad import compat, spch
from pysad.ffscmd import run_sad_string

RE = 2.8179403205e-15
ME = 0.51099895e6
P5 = 5e6


def bg(pev):
    p = pev / ME
    g = math.sqrt(1.0 + p * p)
    return p / g, g


def sobol(n, d, seed):
    return qmc.Sobol(d, seed=seed).random(n)


def gauss(u):
    return norm.ppf(u * (1 - 2e-9) + 1e-9)


def setup(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out


def track(I, beam, cmd="TrackParticles[{1, beam0}]"):
    """beam (n,6) を SAD の TrackParticles で追う。戻り: (座標 (n,6), alive (n,))"""
    n = beam.shape[0]
    I.globals['beam0'] = [list(map(float, beam[:, j])) for j in range(6)] + [[1.0] * n]
    I.run("r0 = " + cmd + "[[2]];", 'ffs')
    r = I.globals['r0']
    return np.array(r[:6]).T, np.array(r[6])


DRIFT_LINE = """MOMENTUM = 5MEV;
DRIFT D1 = (L=%s RADIUS=%s);
MARK ORIGIN = ();
LINE TL = (ORIGIN D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = %s;
"""


def uniform_beam(n, a, lb, seed=3, eps=0.0):
    """半径 a の一様な円 (KV 分布: 4 次元の球面の一様分布)、z は長さ lb の一様分布"""
    u = sobol(n, 5, seed)
    v = gauss(u[:, :4])
    v /= np.sqrt((v * v).sum(axis=1))[:, None]
    b = np.zeros((n, 6))
    b[:, 0] = a * v[:, 0]
    b[:, 2] = a * v[:, 2]
    if eps:
        b[:, 1] = eps / a * v[:, 1]
        b[:, 3] = eps / a * v[:, 3]
    b[:, 4] = lb * (u[:, 4] - 0.5)
    return b


# ============================================================================ ソルバー (解析式との比較)
def test_uniform_cylinder_perveance_gradient():
    """長い一様な円柱 (半径 a、線密度 λ): E'_x = 2λ x/a² (→ x'' = K x/a²、K = 2 r_c λ/(β²γ³))、軸上の電位 λ(1+2 ln(R/a))。
    壁の半径 R を変えても横の場は変わらない (ガウスの法則)"""
    n, a, lb = 32768, 1e-3, 1.0
    u = sobol(n, 3, 1)
    r = a * np.sqrt(u[:, 0])
    x, y, z = r * np.cos(2 * np.pi * u[:, 1]), r * np.sin(2 * np.pi * u[:, 1]), lb * (u[:, 2] - 0.5)
    lam = n / lb
    core = np.abs(z) < 0.25 * lb
    slopes = []
    for R in (0.02, 0.003):
        f = spch.Solver3D().solve(x, y, z, R)
        _, gx, gy, _ = f.grad(x, y, z, selfsub=True)
        sl = -(gx[core] * x[core] + gy[core] * y[core]).sum() / (x[core] ** 2 + y[core] ** 2).sum()
        slopes.append(sl)
        assert sl == pytest.approx(2 * lam / a ** 2, rel=5e-3)
        u0 = f.grad(np.zeros(1), np.zeros(1), np.zeros(1))[0][0]
        assert u0 == pytest.approx(lam * (1 + 2 * math.log(R / a)), rel=2e-3)
    assert slopes[0] == pytest.approx(slopes[1], rel=1e-4)


def bassetti_erskine(x, y, sx, sy):
    """単位線密度の 2 次元ガウス分布の場 (Ex, Ey) (ガウス単位: 遠方で 2/r)"""
    s = math.sqrt(2 * (sx * sx - sy * sy))
    xa, ya = np.abs(x), np.abs(y)
    w1 = wofz((xa + 1j * ya) / s)
    w2 = wofz((xa * sy / sx + 1j * ya * sx / sy) / s)
    f = 2 * math.sqrt(math.pi) / s * (w1 - np.exp(-xa * xa / 2 / sx / sx - ya * ya / 2 / sy / sy) * w2)
    return np.sign(x) * f.imag, np.sign(y) * f.real


def test_gaussian_matches_bassetti_erskine():
    """長い楕円のガウスビーム (σx:σy = 10:3): 横の場 = λ(z) × Bassetti-Erskine。中心の勾配は 2λ/(σx(σx+σy))"""
    n = 65536
    g = gauss(sobol(n, 3, 2))
    sx, sy, sz, R = 1e-3, 0.3e-3, 2.0, 0.05
    x, y, z = sx * g[:, 0], sy * g[:, 1], sz * g[:, 2]
    f = spch.Solver3D().solve(x, y, z, R)
    _, gx, gy, _ = f.grad(x, y, z, selfsub=True)
    lam = n * np.exp(-z * z / 2 / sz / sz) / math.sqrt(2 * math.pi) / sz
    ex, ey = bassetti_erskine(x, y, sx, sy)
    m = (np.abs(x) < 3 * sx) & (np.abs(y) < 3 * sy) & (np.abs(z) < 2 * sz)
    for gg, ee in ((gx, ex), (gy, ey)):
        ref = lam[m] * ee[m]
        assert (-gg[m] * ref).sum() / (ref ** 2).sum() == pytest.approx(1.0, abs=5e-3)
        assert np.sqrt(((-gg[m] - ref) ** 2).mean() / (ref ** 2).mean()) < 0.015
    d = 1e-6
    _, cx, cy, _ = f.grad(np.array([d, -d, 0, 0]), np.array([0, 0, d, -d]), np.zeros(4))
    lam0 = n / math.sqrt(2 * math.pi) / sz
    assert -(cx[0] - cx[1]) / (2 * d) == pytest.approx(2 * lam0 / (sx * (sx + sy)), rel=0.02)
    assert -(cy[2] - cy[3]) / (2 * d) == pytest.approx(2 * lam0 / (sy * (sx + sy)), rel=0.02)


def test_short_bunch_sphere_3d():
    """静止系で球対称のガウス分布 (壁は遠い): E_r = N [erf(r/√2σ) − √(2/π) (r/σ) exp(−r²/2σ²)] / r² (3 次元の解)"""
    n, sg = 65536, 1e-3
    g = gauss(sobol(n, 3, 2))
    x, y, z = sg * g[:, 0], sg * g[:, 1], sg * g[:, 2]
    f = spch.Solver3D().solve(x, y, z, 0.1)
    _, gx, gy, gz = f.grad(x, y, z, selfsub=True)
    r = np.sqrt(x * x + y * y + z * z)
    er = n * (erf(r / math.sqrt(2) / sg) - math.sqrt(2 / math.pi) * r / sg * np.exp(-r * r / 2 / sg / sg)) / r ** 2
    en = -(gx * x + gy * y + gz * z) / r
    m = r < 3 * sg
    assert (en[m] * er[m]).sum() / (er[m] ** 2).sum() == pytest.approx(1.0, abs=5e-3)
    assert np.sqrt(((en - er)[m] ** 2).mean() / (er[m] ** 2).mean()) < 0.015


def test_point_charge_in_pipe_potential():
    """短いバンチ (σ = R/40) の円筒の壁による遮蔽: 軸上の電位 = 円筒の中の点電荷の解
    φ(z) = (2/R) Σ_j exp(−j_0j |z|/R) / (j_0j J1(j_0j)²) (周期の像を足したもの)。z = 3.5σ で自由空間の 1/z より 4.5% 小さい"""
    n, sg, R = 32768, 1e-4, 0.004
    g = gauss(sobol(n, 3, 1))
    x, y, z = sg * g[:, 0], sg * g[:, 1], sg * g[:, 2]
    f = spch.Solver3D(nz=32).solve(x, y, z, R)
    lz = f.phi.shape[2] * f.hz
    jz = jn_zeros(0, 3000)

    def pipe(zz):
        return 2.0 / R * np.sum(np.exp(-jz * abs(zz) / R) / (jz * j1(jz) ** 2))
    zt = 3.5 * sg
    u = f.grad(np.zeros(1), np.zeros(1), np.array([zt]))[0][0] / n
    ref = sum(pipe(zt - m * lz) for m in range(-3, 4)) * erf(3.5 / math.sqrt(2))
    assert u == pytest.approx(ref, rel=3e-3)
    assert ref < 0.97 / zt                      # 壁の遮蔽が効いている条件であること


def test_pipe_image_force_offset_beam():
    """軸から (x0, y0) ずれた細い長いビームが壁の像から受ける力: E = 2λ (x0, y0)/(R² − x0² − y0²)"""
    n, a, lb, R = 16384, 1e-4, 1.0, 0.02
    u = sobol(n, 3, 1)
    r = a * np.sqrt(u[:, 0])
    x, y, z = r * np.cos(2 * np.pi * u[:, 1]), r * np.sin(2 * np.pi * u[:, 1]), lb * (u[:, 2] - 0.5)
    zt = np.linspace(-0.2, 0.2, 5)
    for x0, y0 in ((0.008, 0.0), (0.006, 0.006), (0.0, 0.015)):
        f = spch.Solver3D().solve(x + x0, y + y0, z, R)
        f0 = spch.Solver3D().solve(x + x0, y + y0, z, 1e3)          # 壁なし (像が無視できる)
        _, gx, gy, _ = f.grad(0 * zt + x0, 0 * zt + y0, zt)
        _, hx, hy, _ = f0.grad(0 * zt + x0, 0 * zt + y0, zt)
        c = 2 * n / lb / (R * R - x0 * x0 - y0 * y0)
        assert -(gx - hx).mean() == pytest.approx(c * x0, rel=1e-4, abs=1e-5 * c * R)
        assert -(gy - hy).mean() == pytest.approx(c * y0, rel=1e-4, abs=1e-5 * c * R)


def test_longitudinal_field_g_factor():
    """長いバンチ (一様な円柱 a、ガウスの λ(z)、壁 R) の縦の場: 軸上の電位 λ(z)(1+2 ln(R/a)) → E_z = −λ'(z) (1+2 ln(R/a))、
    r = a では −λ'(z)·2 ln(R/a) (g 因子)"""
    n, a, R, sz = 65536, 1e-3, 0.02, 1.0
    u = sobol(n, 3, 2)
    r = a * np.sqrt(u[:, 0])
    x, y, z = r * np.cos(2 * np.pi * u[:, 1]), r * np.sin(2 * np.pi * u[:, 1]), sz * gauss(u[:, 2])
    f = spch.Solver3D().solve(x, y, z, R)
    zt = np.array([-2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0]) * sz
    dl = -zt / sz ** 2 * n * np.exp(-zt * zt / 2 / sz / sz) / math.sqrt(2 * math.pi) / sz
    gz = f.grad(0 * zt, 0 * zt, zt)[3]
    assert np.allclose(gz / (dl * (1 + 2 * math.log(R / a))), 1.0, atol=0.01)
    gz = f.grad(0 * zt + 0.999 * a, 0 * zt, zt)[3]
    assert np.allclose(gz / (dl * 2 * math.log(R / a)), 1.0, atol=0.01)


# ============================================================================ キックの性質
def _canon_kick(v, field, mode, rcw, al=1.0):
    p0 = P5 / ME
    h0 = math.sqrt(1 + p0 * p0)
    b0 = p0 / h0
    S = np.zeros((1, 11))
    pr = 1 + v[5]
    S[0, :6] = [v[0], v[1] / pr, v[2], v[3] / pr, v[4], v[5]]
    p = p0 * pr
    zeta = np.array([h0 * b0 * v[4] / (p / math.sqrt(1 + p * p))])
    if mode == 3:
        spch.kick3d(S, 1, field, al, rcw, p0, h0, b0, h0, zeta, selfsub=False)
    else:
        spch.kick2d(S, 1, field, al, rcw, p0, h0, b0, h0, selfsub=False)
    pr = 1 + S[0, 5]
    return np.array([S[0, 0], S[0, 1] * pr, S[0, 2], S[0, 3] * pr, S[0, 4], S[0, 5]])


@pytest.mark.parametrize('mode', [3, 2])
def test_frozen_field_kick_is_symplectic(mode):
    """凍結した電位のキックは (x, px, y, py, z, δ) でシンプレクティック (ヤコビ行列 J で Jᵀ S J = S)。
    キックの大きさは J−I の最大の要素が 1 以上 (強いキック) で確かめる"""
    n = 8192
    g = gauss(sobol(n, 5, 5))
    sx, sy, sz = 1e-3, 0.6e-3, 0.02
    x, y, z = sx * g[:, 0], sy * g[:, 2], sz * g[:, 4]
    _, g0 = bg(P5)
    if mode == 3:
        field = spch.Solver3D().solve(x, y, g0 * z, 0.02)
    else:
        field = spch.Solver2D(64, 64, 1, 1.5e-4, 1e-4, 1e-3).solve(x, y, z)
    rcw = RE * 3e10 / n
    sm = np.zeros((6, 6))
    for i in range(3):
        sm[2 * i, 2 * i + 1] = 1.0
        sm[2 * i + 1, 2 * i] = -1.0
    hs = np.array([1e-8, 1e-7, 1e-8, 1e-7, 1e-7, 1e-7])
    for v0 in ([0.3e-3, 1e-4, -0.2e-3, 2e-4, 0.007, 1e-3], [1.1e-3, -3e-4, 0.5e-3, 1e-4, -0.02, -5e-3]):
        v0 = np.array(v0)
        J = np.zeros((6, 6))
        for j in range(6):
            e = np.zeros(6)
            e[j] = hs[j]
            J[:, j] = (_canon_kick(v0 + e, field, mode, rcw) - _canon_kick(v0 - e, field, mode, rcw)) / (2 * hs[j])
        assert np.abs(J - np.eye(6)).max() > 0.5
        assert np.abs(J.T @ sm @ J - sm).max() < 1e-7
        assert abs(np.linalg.det(J) - 1.0) < 1e-7


def test_total_momentum_nearly_conserved_in_free_space():
    """壁が遠いとき、バンチ全体の運動量の和は変わらない (メッシュの誤差の範囲: 力の絶対値の和の 1e-3 以下)"""
    n = 8192
    g = gauss(sobol(n, 5, 5))
    p0 = P5 / ME
    h0 = math.sqrt(1 + p0 * p0)
    S = np.zeros((n, 11))
    S[:, 0], S[:, 2], S[:, 4] = 1e-3 * g[:, 0], 0.6e-3 * g[:, 2], 0.02 * g[:, 4]
    b0, g0, zeta = spch.bunch_frame(S, n, p0)
    f = spch.Solver3D().solve(S[:, 0], S[:, 2], zeta, 10.0)
    spch.kick3d(S, n, f, 1.0, RE * 3e10 / n, p0, h0, b0, g0, zeta)
    for col in (1, 3, 5):
        assert abs(S[:, col].sum()) < 1e-3 * np.abs(S[:, col]).sum()
        assert np.abs(S[:, col]).sum() > 0.0


# ============================================================================ SAD スクリプト (TrackParticles)
def test_spac_kv_envelope_in_drift(v2_mode):
    """SPAC: KV ビーム (半径 1 mm、エミッタンス 2e-6 m、1e11 個/0.1 m、5 MeV) の 1 m のドリフト。包絡方程式
    a'' = K/a + ε²/a³ (K = 2 r_e λ/(β²γ³) = 5.98e-6) の数値積分 a(1 m) = 3.996 mm と、追跡の 2·rms x が 0.5% で合う。
    NOSPAC は a'' = ε²/a³ (2.236 mm)"""
    n, a, eps, lb, npb = 4096, 1e-3, 2e-6, 0.1, 1e11
    b, g = bg(P5)
    K = 2 * RE * (npb / lb) / (b * b * g ** 3)
    beam = uniform_beam(n, a, lb, eps=eps)
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '0.02', '1e11'))
    core = np.abs(beam[:, 4]) < 0.3 * lb
    I.run("SPAC;", 'ffs')
    c, alive = track(I, beam)
    assert alive.sum() == n
    ode = solve_ivp(lambda t, y: [y[1], K / y[0] + eps ** 2 / y[0] ** 3], [0, 1.0], [a, 0], rtol=1e-10, atol=1e-14)
    assert ode.y[0, -1] == pytest.approx(3.9958e-3, rel=1e-3)
    assert 2 * c[core, 0].std() == pytest.approx(ode.y[0, -1], rel=5e-3)
    assert 2 * c[core, 2].std() == pytest.approx(ode.y[0, -1], rel=5e-3)
    assert 2 * c[core, 1].std() == pytest.approx(ode.y[1, -1], rel=0.02)          # 包絡の傾き a'
    I.run("NOSPAC;", 'ffs')
    c0, _ = track(I, beam)
    ode0 = solve_ivp(lambda t, y: [y[1], eps ** 2 / y[0] ** 3], [0, 1.0], [a, 0], rtol=1e-10, atol=1e-14)
    assert 2 * c0[:, 0].std() == pytest.approx(ode0.y[0, -1], rel=1e-3)


def test_spac_zero_charge_and_gamma_scaling(v2_mode):
    """極限: PBUNCH = 0 ではドリフトと同じ。短いドリフト (1 キック) の傾き px/x は K L/a² で、K ∝ 1/(β²γ³)
    (5 MeV と 50 MeV で 1e3 倍違う: γ⁻² の打ち消し (電場と磁場) と γ⁻¹ (慣性))"""
    n, a, lb = 4096, 1e-3, 0.1
    beam = uniform_beam(n, a, lb, eps=1e-6)
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '0.02', '0'))
    I.run("NOSPAC;", 'ffs')
    c0, _ = track(I, beam)
    I.run("SPAC;", 'ffs')
    c1, _ = track(I, beam)
    assert np.abs(c1 - c0).max() < 1e-15
    beam = uniform_beam(n, a, lb)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    sl = {}
    for mom, pev in (('5MEV', 5e6), ('50MEV', 50e6)):
        ffs, I, out = setup((DRIFT_LINE % ('0.05', '0.02', '1e11')).replace('5MEV', mom))
        I.run("SPAC;", 'ffs')
        c, _ = track(I, beam)
        b, g = bg(pev)
        K = 2 * RE * (1e11 / lb) / (b * b * g ** 3)
        sl[mom] = (c[core, 1] * beam[core, 0]).sum() / (beam[core, 0] ** 2).sum()
        assert sl[mom] == pytest.approx(K * 0.05 / a ** 2, rel=0.01)
    assert sl['5MEV'] / sl['50MEV'] == pytest.approx(1.0e3, rel=0.03)


def test_spac_k64_mode_unchanged(k64_mode):
    """k64 モード: SPAC は今までどおり追跡に入らない (NOSPAC とビット単位で同じ)"""
    beam = uniform_beam(512, 1e-3, 0.1, eps=1e-6)
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '0.02', '1e11'))
    I.run("NOSPAC;", 'ffs')
    c0, a0 = track(I, beam)
    I.run("SPAC; PSPAC;", 'ffs')
    c1, a1 = track(I, beam)
    assert np.array_equal(c0, c1) and np.array_equal(a0, a1)


def test_spac_negative_radius_is_aperture_only(v2_mode):
    """RADIUS<0: 空間電荷なし、半径 |RADIUS| の円の絞りだけ。RADIUS>0 は空間電荷と絞り。NOSPAC では絞りも入らない"""
    n = 2048
    beam = uniform_beam(n, 1e-3, 0.1)
    inside = np.hypot(beam[:, 0], beam[:, 2]) <= 0.6e-3
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '-0.6e-3', '1e11'))
    I.run("SPAC;", 'ffs')
    c, alive = track(I, beam)
    assert np.array_equal(alive > 0, inside)
    assert np.allclose(c[inside, 0], beam[inside, 0], atol=1e-15) and np.abs(c[inside, 1]).max() == 0.0
    I.run("NOSPAC;", 'ffs')
    c, alive = track(I, beam)
    assert alive.sum() == n
    # FFS で RADIUS を変えられる (正にすると空間電荷が入る)
    I.run('SPAC; Element["RADIUS", "D1"] = 0.02;', 'ffs')
    c, alive = track(I, beam)
    assert alive.sum() == n and c[:, 0].std() > 2.5 * beam[:, 0].std()


def test_spac_zero_bunch_length_error(v2_mode):
    """すべての粒子の z が同じ (線密度が決まらない) とき、NORFSW のときは、理由を書いたエラーにする"""
    beam = uniform_beam(256, 1e-3, 0.1)
    beam[:, 4] = 0.0
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '0.02', '1e11'))
    I.run("SPAC;", 'ffs')
    I.globals['beam0'] = [list(map(float, beam[:, j])) for j in range(6)] + [[1.0] * 256]
    I.run("r1 = TrackParticles[{1, beam0}];", 'ffs')
    assert 'バンチ長が 0' in out.getvalue()
    # NORFSW (区間ごとに z を 0 にする) では、空間電荷を入れる前にエラーにする。RADIUS<0 (絞りだけ) は NORFSW でも走る
    beam = uniform_beam(256, 1e-3, 0.1)
    I.run("NORFSW;", 'ffs')
    I.globals['beam0'] = [list(map(float, beam[:, j])) for j in range(6)] + [[1.0] * 256]
    I.run("r2 = TrackParticles[{1, beam0}];", 'ffs')
    assert 'RFSW を ON' in out.getvalue() and 'r2' not in I.globals
    I.run('Element["RADIUS", "D1"] = -0.01; r3 = TrackParticles[{1, beam0}];', 'ffs')
    assert sum(I.globals['r3'][1][6]) == 256


def test_spac_solenoid_brillouin_flow(v2_mode):
    """SOL 領域の中のドリフト (spdrift_solenoid): 場の外で作ったビーム (正準角運動量 0) の包絡は
    a'' = −k_L² a + K/a (k_L = bz/2)。k_L² = K/a² (Brillouin 流、BZ = 0.0816 T) なら半径は変わらない。
    0.7 倍の磁場では 1.81 mm に広がる (空間電荷が無ければ 0.14 mm に絞られる)"""
    n, a, lb, npb = 4096, 1e-3, 0.1, 1e11
    b, g = bg(P5)
    K = 2 * RE * (npb / lb) / (b * b * g ** 3)
    brho = P5 / 299792458.0
    kl = math.sqrt(K) / a
    beam = uniform_beam(n, a, lb)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    for fac in (1.0, 0.7):
        bz = 2 * kl * fac * brho
        src = """MOMENTUM = 5MEV;
DRIFT D1 = (L=1.0 RADIUS=0.02);
SOL SA = (L=0 BZ=%.12g BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
MARK ORIGIN = ();
LINE TL = (ORIGIN SA D1 SB);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
""" % bz
        ffs, I, out = setup(src)
        I.run("SPAC;", 'ffs')
        c, alive = track(I, beam)
        assert alive.sum() == n
        ode = solve_ivp(lambda t, y: [y[1], -(kl * fac) ** 2 * y[0] + K / y[0]], [0, 1.0], [a, 0], rtol=1e-10, atol=1e-14)
        rr = math.sqrt(2.0 * (c[core, 0] ** 2 + c[core, 2] ** 2).mean())
        assert rr == pytest.approx(ode.y[0, -1], rel=0.01)
        if fac == 1.0:
            assert rr == pytest.approx(a, rel=0.01)
        else:
            assert ode.y[0, -1] == pytest.approx(1.8087e-3, rel=1e-3)
        I.run("NOSPAC;", 'ffs')
        c0, _ = track(I, beam)
        ode0 = solve_ivp(lambda t, y: [y[1], -(kl * fac) ** 2 * y[0]], [0, 1.0], [a, 0], rtol=1e-10, atol=1e-14)
        assert math.sqrt(2.0 * (c0[:, 0] ** 2 + c0[:, 2] ** 2).mean()) == pytest.approx(abs(ode0.y[0, -1]), rel=2e-3)
        assert abs(ode0.y[0, -1]) < (0.8 if fac == 1.0 else 0.1) * ode.y[0, -1]      # 空間電荷の有無で大きく違う条件であること


# ---- リングのチューンシフト
LQ, LD, KQ, NC = 0.2, 0.8, 1.2, 4
RING = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.8 RADIUS=0.03);
QUAD QF = (L=0.2 K1=1.2) QD = (L=0.2 K1=-1.2);
SPCH SP = ();
MARK ORIGIN = ();
LINE CELL = (QF D1 SP QD D1 SP);
LINE RING = (ORIGIN 4*CELL);
FFS;
RFSW; NORAD; NOFLUC;
USE RING; CELL; CALC;
SPACSTEP = 0.1;
PSPACNX = 64; PSPACNY = 64; PSPACDX = 4e-4; PSPACDY = 4e-4;
"""


def _mq(k, l):
    if k > 0:
        w = math.sqrt(k)
        return np.array([[math.cos(w * l), math.sin(w * l) / w], [-w * math.sin(w * l), math.cos(w * l)]])
    w = math.sqrt(-k)
    return np.array([[math.cosh(w * l), math.sinh(w * l) / w], [w * math.sinh(w * l), math.cosh(w * l)]])


def _md(l):
    return np.array([[1.0, l], [0.0, 1.0]])


def _twiss(m):
    mu = math.acos(0.5 * (m[0, 0] + m[1, 1]))
    s = math.sin(mu)
    if m[0, 1] < 0:
        mu, s = 2 * math.pi - mu, -s
    return mu, m[0, 1] / s, (m[0, 0] - m[1, 1]) / (2 * s)


def _ring_model(kp, eps, tw, kicks):
    """自前の 2x2 行列のリング。空間電荷は、ガウスビームの中心の勾配 kp/(σx(σx+σy)) の薄い線形レンズ (σ は裸の光学から)。
    kicks='drift': ドリフトを 0.1 m に割った間 (SPAC)、'spch': 各ドリフトの後に、前のレンズからの距離の分 (PSPAC)。
    戻り: (νx, νy) の小数部 (acos)"""
    (bx, ax), (by, ay) = tw
    mx, my = np.eye(2), np.eye(2)
    sgx = eps * np.array([[bx, -ax], [-ax, (1 + ax * ax) / bx]])
    sgy = eps * np.array([[by, -ay], [-ay, (1 + ay * ay) / by]])

    def app(a, b):
        nonlocal mx, my, sgx, sgy
        mx, my, sgx, sgy = a @ mx, b @ my, a @ sgx @ a.T, b @ sgy @ b.T

    def lens(al):
        nonlocal mx, my
        sx, sy = math.sqrt(sgx[0, 0]), math.sqrt(sgy[0, 0])
        mx = np.array([[1, 0], [kp * al / (sx * (sx + sy)), 1]]) @ mx
        my = np.array([[1, 0], [kp * al / (sy * (sx + sy)), 1]]) @ my
    for _ in range(NC):
        for sgn in (1, -1):
            k = KQ / LQ * sgn
            app(_mq(k, LQ), _mq(-k, LQ))
            if kicks == 'drift':
                nd = round(LD / 0.1)
                app(_md(0.05), _md(0.05))
                for i in range(nd):
                    if i > 0:
                        app(_md(0.1), _md(0.1))
                    lens(0.1)
                app(_md(0.05), _md(0.05))
            else:
                app(_md(LD), _md(LD))
                if kicks == 'spch':
                    lens(LQ + LD)
    return [math.acos(0.5 * np.trace(m)) / 2 / math.pi for m in (mx, my)]


@pytest.mark.parametrize('flag', ['SPAC', 'PSPAC'])
def test_ring_tune_shift_laslett(v2_mode, flag):
    """4 セルの FODO リング (周長 8 m、5 MeV、ν = 0.754)、整合したガウスビーム (ε = 1e-6 m、σz = 5 cm、2.5e8 個)。
    中心の小振幅の粒子の 1 周の行列からチューンを求め、空間電荷のチューンシフト (Laslett: Δν = −(1/4π)∮β K_sc ds、
    K_sc = 2 r_e λ/(β²γ³)/(σx(σx+σy))) を、自前の行列に線形レンズを入れたリングと比べる (Δν = −0.003、5% で一致)。
    SPAC はドリフトの中だけ (6.4 m)、PSPAC は SPCH 要素に全周 (8 m) の分をまとめて入れる"""
    b, g = bg(P5)
    eps, sz, npb, n = 1e-6, 0.05, 2.5e8, 4096
    k = KQ / LQ
    cx = _md(LD) @ _mq(-k, LQ) @ _md(LD) @ _mq(k, LQ)
    cy = _md(LD) @ _mq(k, LQ) @ _md(LD) @ _mq(-k, LQ)
    (mux, bx, ax), (muy, by, ay) = _twiss(cx), _twiss(cy)
    v = gauss(sobol(n, 5, 5))
    beam = np.zeros((n, 6))
    beam[:, 0] = math.sqrt(eps * bx) * v[:, 0]
    beam[:, 1] = math.sqrt(eps / bx) * (v[:, 1] - ax * v[:, 0])
    beam[:, 2] = math.sqrt(eps * by) * v[:, 2]
    beam[:, 3] = math.sqrt(eps / by) * (v[:, 3] - ay * v[:, 2])
    beam[:, 4] = sz * v[:, 4]
    ffs, I, out = setup(RING)
    I.globals['PBUNCH'] = npb * (n + 1) / n
    d = 1e-6

    def tunes():
        p = np.zeros((8, 6))
        for i in range(8):                       # 粒子同士の力を避けるため、試験粒子は 1 個ずつ
            pr = np.zeros((1, 6))
            pr[0, i // 2] = d if i % 2 == 0 else -d
            c, alive = track(I, np.vstack([beam, pr]), "TrackParticles[{1, beam0}, 1, 1, 1]")
            assert alive.sum() == n + 1
            p[i] = c[n]
        res = []
        for o in (0, 2):
            r = 2 * o                             # 試験粒子の行: x± (0,1)、px± (2,3)、y± (4,5)、py± (6,7)
            m = np.array([[p[r, o] - p[r + 1, o], p[r + 2, o] - p[r + 3, o]],
                          [p[r, o + 1] - p[r + 1, o + 1], p[r + 2, o + 1] - p[r + 3, o + 1]]]) / (2 * d)
            assert np.linalg.det(m) == pytest.approx(1.0, abs=2e-3)
            res.append(math.acos(0.5 * np.trace(m)) / 2 / math.pi)
        return res
    I.run("NOSPAC; NOPSPAC;", 'ffs')
    nu0 = tunes()
    bare = _ring_model(0.0, eps, ((bx, ax), (by, ay)), 'none')
    assert nu0 == pytest.approx(bare, abs=1e-8)          # 裸のチューン (自前の行列と追跡)
    assert bare[0] == pytest.approx(1 - 0.754366, abs=1e-6)
    I.run(flag + ";", 'ffs')
    nu1 = tunes()
    kp = 2 * RE * npb / math.sqrt(2 * math.pi) / sz / (b * b * g ** 3)
    mod = _ring_model(kp, eps, ((bx, ax), (by, ay)), 'drift' if flag == 'SPAC' else 'spch')
    for i in range(2):
        dnu_model = mod[i] - bare[i]                     # acos の小数部は 1−ν なので、正がチューンの低下
        assert dnu_model == pytest.approx(0.0030 if flag == 'SPAC' else 0.0038, abs=3e-4)
        assert nu1[i] - nu0[i] == pytest.approx(dnu_model, rel=0.05)


def test_pspac_spch_kick_uniform_beam(v2_mode):
    """PSPAC + SPCH: 一様なビーム (半径 1 mm) の 1 キックの傾き px/x = K L/a²。L は前の SPCH (または区間の始め) からの距離で、
    区間の終わりでは最後の SPCH からの残りを入れる。PSPACNZ を与えたとき (K64 と同じ固定メッシュ) と自動のメッシュが合う"""
    n, a, lb, npb = 8192, 1e-3, 0.1, 1e11
    b, g = bg(P5)
    K = 2 * RE * (npb / lb) / (b * b * g ** 3)
    beam = uniform_beam(n, a, lb)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    src = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.3) D2 = (L=0.2);
SPCH SP1 = ();
MARK ORIGIN = ();
LINE TL = (ORIGIN SP1 D1 D2);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
PSPACNX = 128; PSPACNY = 128; PSPACDX = 1e-4; PSPACDY = 1e-4;
"""
    ffs, I, out = setup(src)
    I.run("PSPAC;", 'ffs')
    res = []
    for extra in ("", "PSPACNZ = 40; PSPACDZ = 0.004;"):
        if extra:
            I.run(extra, 'ffs')
        c, alive = track(I, beam)
        assert alive.sum() == n
        assert np.abs(c[:, 0] - beam[:, 0]).max() == 0.0             # キックは区間の終わり (0.5 m の分)
        for o in (0, 2):
            sl = (c[core, o + 1] * beam[core, o]).sum() / (beam[core, o] ** 2).sum()
            assert sl == pytest.approx(K * 0.5 / a ** 2, rel=0.01)
            res.append(sl)
        assert np.abs(c[core, 5]).max() < 2e-5                       # 中央部は λ' = 0: エネルギーの蹴りはほぼ 0
    assert res[0] == pytest.approx(res[2], rel=0.01)
    # SPCH が区間の途中にあるとき: SPCH の位置までの分と、残りの分に分かれる (和は同じ長さ)
    src2 = src.replace("LINE TL = (ORIGIN SP1 D1 D2)", "LINE TL = (ORIGIN D1 SP1 D2)")
    ffs, I, out = setup(src2)
    I.run("PSPAC;", 'ffs')
    c, _ = track(I, beam)
    sl1 = K * 0.3 / a ** 2
    x1 = 1.0 + sl1 * 0.2                                            # 0.3 m の分のキックの後、0.2 m のドリフトで x が伸びる
    assert (c[core, 0] * beam[core, 0]).sum() / (beam[core, 0] ** 2).sum() == pytest.approx(x1, rel=3e-3)
    assert (c[core, 1] * beam[core, 0]).sum() / (beam[core, 0] ** 2).sum() == pytest.approx(sl1 + K * 0.2 / a ** 2 / x1, rel=0.02)


def test_pspac_longitudinal_kick(v2_mode):
    """PSPAC の縦のキック: Δγ = −r_e N L/γ² · (dλ/dz) Φ(x,y)/N。ガウスの λ(z)、一様な円柱 a、矩形の境界 (半幅 b):
    軸上の Φ/N = 1 + 2 ln(1.08 b/a) (正方形の導体の箱は半径 1.08 b の円筒に当たる)"""
    n, a, sz, npb = 16384, 1e-3, 0.02, 1e11
    b, g = bg(P5)
    u = sobol(n, 3, 4)
    r = a * np.sqrt(u[:, 0])
    beam = np.zeros((n, 6))
    beam[:, 0], beam[:, 2], beam[:, 4] = r * np.cos(2 * np.pi * u[:, 1]), r * np.sin(2 * np.pi * u[:, 1]), sz * gauss(u[:, 2])
    src = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.5);
SPCH SP1 = ();
MARK ORIGIN = ();
LINE TL = (ORIGIN SP1 D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
PSPACNX = 128; PSPACNY = 128; PSPACDX = 1e-4; PSPACDY = 1e-4;
PSPAC;
"""
    ffs, I, out = setup(src)
    c, _ = track(I, beam)
    hb = 64 * 1e-4
    gfac = 1 + 2 * math.log(1.0787 * hb / a)
    near = (np.hypot(beam[:, 0], beam[:, 2]) < 0.15 * a) & (np.abs(beam[:, 4]) < 2 * sz)
    z = beam[near, 4]
    dlam = -z / sz ** 2 * npb * np.exp(-z * z / 2 / sz / sz) / math.sqrt(2 * math.pi) / sz
    dgam = -RE * 0.5 / g ** 2 * dlam * gfac
    ddel = dgam / (b * b * g)                       # Δδ = Δγ/(β²γ)
    assert (c[near, 5] * ddel).sum() / (ddel ** 2).sum() == pytest.approx(1.0, abs=0.02)
    assert np.abs(ddel).max() > 1e-5


# ============================================================================ 使い方の例 (docs/history/v2_changes_C_spch.md と同じスクリプト)
EXAMPLE_LINE = r"""MOMENTUM = 5MEV;
! 輸送路: 真空ダクト半径 2 cm のドリフト 1 m (RADIUS>0 のドリフトに空間電荷が入る)
DRIFT D1 = (L=1.0 RADIUS=0.02);
MARK  ORIGIN = ();
LINE TL = (ORIGIN D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
$FORM = "12.6";
np = 2000;
PBUNCH = 1e11;                   ! 1 バンチの実粒子数 (マクロ粒子 1 個は PBUNCH/np 個)
fr[x_] := x - Floor[x];
a = 1e-3;                        ! 半径 1 mm の一様な円、長さ 0.1 m の一様なバンチ
rr = Table[a*Sqrt[fr[i*0.6180339887+0.5]], {i, np}];
th = Table[2*Pi*fr[i*0.7548776662], {i, np}];
x0 = rr*Cos[th]; y0 = rr*Sin[th];
z0 = Table[0.1*(fr[i*0.5698402910]-0.5), {i, np}];
zero = Table[0, {i, np}];
beam = {x0, zero, y0, zero, z0, zero, Table[1, {i, np}]};
rms[v_] := Sqrt[Plus@@(v^2)/Length[v]];
NOSPAC;
r = TrackParticles[{1, beam}][[2]];
Print["RES NOSPAC rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
SPAC;
r = TrackParticles[{1, beam}][[2]];
Print["RES SPAC   rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
"""

EXAMPLE_RING = r"""MOMENTUM = 5MEV;
! 小さなリング: 4 セルの FODO (周長 8 m)。SPAC はドリフト (RADIUS>0) の中、PSPAC は SPCH 要素の位置にまとめて入れる
DRIFT D1 = (L=0.8 RADIUS=0.03);
QUAD  QF = (L=0.2 K1=1.2) QD = (L=0.2 K1=-1.2);
SPCH  SP = ();
MARK  ORIGIN = ();
LINE CELL = (QF D1 SP QD D1 SP);
LINE RING = (ORIGIN 4*CELL);
FFS;
RFSW; NORAD; NOFLUC;
USE RING; CELL; CALC;
$FORM = "12.6";
np = 1000;
PBUNCH = 2.5e9;
SPACSTEP = 0.1;                  ! SPAC のキックの間隔 [m] (既定 0.05)
PSPACNZ = 64; PSPACDZ = 0.008; PSPACNX = 64; PSPACNY = 64; PSPACDX = 4e-4; PSPACDY = 4e-4;      ! PSPAC の横のメッシュ (±12.8 mm)
fr[x_] := x - Floor[x];
gs[k_, c_] := Table[Sqrt[-2*Log[fr[i*k+0.5]+1e-12]]*Cos[2*Pi*fr[i*c]], {i, np}];   ! ガウス分布 (乱数を使わない)
eps = 1e-6; bx = Twiss["BX", 1]; ax = Twiss["AX", 1]; by = Twiss["BY", 1]; ay = Twiss["AY", 1];
u1 = gs[0.6180339887, 0.7548776662]; u2 = gs[0.5698402910, 0.8191725134];
u3 = gs[0.6710436067, 0.4502835858]; u4 = gs[0.8566748839, 0.7338918566];
x0 = Sqrt[eps*bx]*u1; px0 = Sqrt[eps/bx]*(u2 - ax*u1);
y0 = Sqrt[eps*by]*u3; py0 = Sqrt[eps/by]*(u4 - ay*u3);
z0 = 0.05*gs[0.7071067812, 0.3247179572];
zero = Table[0, {i, np}];
beam = {x0, px0, y0, py0, z0, zero, Table[1, {i, np}]};
rms[v_] := Sqrt[Plus@@(v^2)/Length[v]];
Print["RES initial rms x=", rms[x0], " y=", rms[y0]];
NOSPAC; NOPSPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES NOSPAC 3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
SPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES SPAC   3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
NOSPAC; PSPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES PSPAC  3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
PSPACNTURNCALC = 1;              ! 1 周目の電位を凍結して使い回す
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES PSPAC frozen rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
"""


def _res(out):
    """Print した 'RES <名前> ... rms x= .. y= .. [alive= ..]' の行を {名前: (x, y, alive)} にする"""
    d = {}
    for l in out.getvalue().splitlines():
        if l.startswith('RES '):
            name = l[4:l.index('rms x=')].strip()
            v = l[l.index('rms x=') + 6:].replace('y=', ' ').replace('alive=', ' ').split()
            d[name] = tuple(float(t) for t in v)
    return d


def test_example_transport_line_script(v2_mode):
    """輸送路の例 (SAD スクリプトだけで走る): 5 MeV、1e11 個、半径 1 mm の一様なビーム、1 m のドリフト (RADIUS = 2 cm)。
    解析式 (包絡方程式 a'' = K/a) は rms = a/2 = 1.639 mm。gs (K64) は x 0.623 mm、y 0.725 mm (力が 1/γ 倍に弱い。
    力を K/γ にした解析式は 0.645 mm)"""
    ffs, I, out = setup(EXAMPLE_LINE)
    r = _res(out)
    assert r['NOSPAC'] == pytest.approx((5.002783e-4, 4.997272e-4, 2000.0), rel=1e-6)
    assert r['SPAC'][0] == pytest.approx(1.6387e-3, rel=0.02)
    assert r['SPAC'][1] == pytest.approx(1.6387e-3, rel=0.02)
    assert r['SPAC'][2] == 2000.0


def test_example_ring_script(v2_mode):
    """小さなリングの例 (4 セルの FODO、5 MeV、2.5e9 個、3 周): SPAC (ドリフトの中の 3 次元の解) と PSPAC (SPCH の位置の
    2 次元の解) が合い、gs (K64) の PSPAC の値 (rms x 1.776 mm、y 1.042 mm。凍結は 1.779、1.041) とも合う。
    gs の SPAC は x 1.750 mm、y 1.034 mm で、空間電荷が弱い (NOSPAC は 1.750、1.024)"""
    ffs, I, out = setup(EXAMPLE_RING)
    r = _res(out)
    assert r['NOSPAC 3 turns'] == pytest.approx((1.750e-3, 1.024e-3, 1000.0), rel=5e-4)
    for name in ('SPAC   3 turns', 'PSPAC  3 turns', 'PSPAC frozen'):
        assert r[name][0] == pytest.approx(1.776e-3, rel=2e-3)
        assert r[name][1] == pytest.approx(1.042e-3, rel=2e-3)
        assert r[name][2] == 1000.0


def test_element_radius_readable(v2_mode):
    """DRIFT の RADIUS は Element["RADIUS", 名前] で読める (gs と同じ)"""
    ffs, I, out = setup(DRIFT_LINE % ('1.0', '0.02', '1e11') + 'Print["RES r rms x=", Element["RADIUS", "D1"], " y=", 0];' + chr(10))
    assert _res(out)['r'][0] == pytest.approx(0.02)
