# -*- coding: utf-8 -*-
"""v2.4.0 担当 D3: B167 の残り (RF の電圧の β の因子・PhotonList の β の因子・独立な (x, px, δ, z) 写像との εx・σδ)、
B146 の ψ (SOL 領域内の四極の本体の B∥)、SOL 領域内の TRPT の加速 (POL)、`Bend.transfer_matrix` の K0≠0、B145 の確認。

独立な基準 (pysad の光学・放射・追跡のコードを使わない):
  - 同期位相: 1 周でエネルギーを失う粒子の δ が閉じる (空洞の仕事 = 損失) ことを追跡 (RF の 1 周、RAD) で確かめる。
  - 光子のエネルギーの総和 = 粒子が失ったエネルギー (エネルギー保存)。
  - `_rad_beta_indep.Ring4B`: (x, px, δ, z) の放射つき 1 周写像 (RK4、空洞はエネルギー kick) + 量子励起の Lyapunov 方程式 (厳密なスペクトルの r2)。
  - `_spin_maxwell_indep2.make_field` の Maxwell 場 (四極のソフトなフリンジ) + 一様なソレノイド場の Thomas-BMT の積分。
  - 一様な E_z を持つ空洞 + ソレノイドの運動方程式・Thomas-BMT (E×β の項を含む) の積分。
陽子の放射は 1 周 1e-17 と弱く 1 周写像の微分が数値的に条件が悪いので、電荷 z = 1000 (古典半径 ∝ z²、空洞の仕事も z 倍) で放射を強め、独立な写像は
放射をさらに z_I = 1.7e5 に強めて解く (平衡の εx・σδ は z によらない)。RF は同期位相の安定な側 (PHI = π) に置く。"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, radbeta, track6

from . import _rad_beta_indep as RB
from . import _spin_maxwell_indep2 as M

MP = 0.93827208816e9
RP = 1.5346982e-18
LAM = 1.054571817e-34 * 299792458 / 1.602176634e-19 / MP        # ħ/(mc) [m]
TH = 2 * math.pi / 24


def _ring(pc_gev, volt=None, phi=0.0, e1e2=True):
    cav = ('CAVI RF=(L=0 VOLT=%.8E HARM=24 PHI=%.15g);\nMARK ST=();\nLINE RING=(ST RF 24*CELL);' % (volt, phi)) if volt is not None else 'LINE RING=(24*CELL);'
    e = TH / 2 if e1e2 else 0.0
    src = f'''MOMENTUM = {pc_gev} GEV;
DRIFT D1=(L=1.0) D2=(L=0.5);
QUAD QF=(L=0.3 K1=0.7) QD=(L=0.3 K1=-0.6);
BEND B=(L=1.2 ANGLE={TH:.16g} E1={e:.16g} E2={e:.16g});
LINE CELL=(QF D2 B D2 QD D1);
{cav}
FFS;
USE RING;
Exit;
'''
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return ffs


@pytest.fixture
def v2():
    compat.set_mode('v2')
    yield
    compat.set_mode(None)


# ======================================================================= B232 (RF の電圧の β の因子・電荷)
def test_rf_synchronous_phase_closes_delta_for_proton(v2):
    """陽子 (3 GeV/c、β = 0.954、z = 1000、V = 20 V) の同期位相 (Emittance[]、RADCOD OFF の trf0) で、RF + RAD の追跡 1 周の δ が閉じる
    (空洞の仕事 z V sinφ = 1 周のエネルギー損失)。K64 は損失を Δp c (運動量) で電圧と比べるので同期位相が 1/β = 5% ずれ、
    電荷も掛けない (z = 1000 で 1000 倍ずれる)。v2: 閉じ方 1e-6 以下 (損失 U0/pc に対する比)。β だけ入れない場合は 0.05"""
    ffs = _ring(3.0, 20.0)
    pc = 3e9

    def closure():
        e = emit6d.emittance6d(ffs.lat, pc, charge=1000.0, amass=MP, rfsw=True, radcod=False)
        c0 = np.array(e['orbit_entrance'], dtype=float)
        c0[4] = 0.0
        o = track6.track_particles(ffs.lat.elements, c0[None, :], 1, pc, e['trf0'], True, True, False, charge=1000.0, amass=MP, circ=91.2)
        return (o[0, 5] - c0[5]) / (e['energy_loss_eV'] / pc)
    assert abs(closure()) < 1e-6
    orig = emit6d._u0_energy
    try:
        emit6d._u0_energy = lambda st, ph: st.u0 * ph.pgev / abs(ph.charge)          # β の因子なし (電荷だけ)
        assert closure() == pytest.approx(0.05, rel=0.05)
        emit6d._u0_energy = lambda st, ph: st.u0 * ph.pgev                           # K64 の式
        assert abs(closure()) > 100.0
    finally:
        emit6d._u0_energy = orig


# ======================================================================= B233 (PhotonList の β の因子)
@pytest.mark.parametrize('pc,z', [(6.0, 100.0), (3.0, 20.0)])
def test_photon_energies_sum_to_energy_loss_for_proton(v2, pc, z):
    """陽子の追跡 (FLUC + PHOTONS、24 セル 1 周) の光子のエネルギーの総和 = 粒子が失ったエネルギー (エネルギー保存、dE = v dp)。
    以前は光子のエネルギー = dp m γ で 1/β² 倍 (γ = 6.5 で 2.5%、3.35 で 10%) 違った"""
    ffs = _ring(3.0, e1e2=True)
    n = 30
    sink = []
    o = track6.track_particles(ffs.lat.elements, np.zeros((n, 6)), 1, pc * 1e9, 0.0, False, True, True, seed=3, charge=z, amass=MP,
                               circ=91.2, photons=sink, photon_geo=ffs.ensure_geo())
    p0 = pc * 1e9 / MP
    E = lambda d: MP * np.sqrt(1 + (p0 * (1 + d)) ** 2)
    loss = E(0.0) - E(o[:, 5])
    pl = np.array(sink).reshape(-1, 12)
    assert len(pl) > 500
    assert pl[:, 0].sum() / loss.sum() == pytest.approx(1.0, abs=6e-3)


# ======================================================================= B167: 独立な (x, px, δ, z) 写像との εx・σδ
def _independent(pc_gev, zi, gain_eV, ds):
    elems = []
    for _ in range(24):
        elems += [dict(type='quad', L=0.3, k1=0.7 / 0.3), dict(type='drift', L=0.5), dict(type='bend', L=1.2, h=TH / 1.2),
                  dict(type='drift', L=0.5), dict(type='quad', L=0.3, k1=-0.6 / 0.3), dict(type='drift', L=1.0)]
    r4 = RB.Ring4B(pc_gev * 1e9 / MP, zi * zi * RP, LAM, elems, radbeta.r2_of_gamma, gain_eV / MP, 24, ds=ds)
    r4.solve()
    S = r4.equilibrium()
    return r4.mode_emittance(0.3, 0.5), math.sqrt(S[2, 2]), r4


@pytest.mark.parametrize('pc,tol_e,tol_s', [(3.0, 4e-3, 1e-3), (1.5, 1e-2, 3e-3)])
def test_proton_emittance_and_sigma_delta_match_independent_map(v2, pc, tol_e, tol_s):
    """3 GeV/c の陽子 (γ = 3.35、β = 0.954)、24 セル FODO (扇形)、RF あり (h = 24): Emittance[] (z = 1000、EMITDIV = 4) の εx・σδ が、独立な
    (x, px, δ, z) 写像 (z_I = 1.7e5、刻み 0.02 m) と εx 0.4%・σδ 0.1% で一致 (刻み 0.0125 m では 4e-5・8e-6。1.5 GeV/c (β = 0.848) も 6e-5・4e-4。ここでは 1.5 GeV/c を εx 1%・σδ 0.3% で確認)。
    εx の因子 r2/β⁴ (B167) の独立な確認 (従来は Emittance[] と RADINT が互いに一致するだけだった)"""
    z = 1000.0
    e0 = emit6d.emittance6d(_ring(pc, e1e2=False).lat, pc * 1e9, charge=z, amass=MP, rfsw=False, radcod=False)
    u0 = e0['energy_loss_eV']
    gain = 1e5 * u0                                                   # z V = 1e5 U0 (同期位相ほぼ 0)
    ffs = _ring(pc, gain / z, phi=math.pi, e1e2=False)
    e = emit6d.emittance6d(ffs.lat, pc * 1e9, charge=z, amass=MP, rfsw=True, radcod=False)
    zi = 1.7e5
    exm, sd, r4 = _independent(pc, zi, gain, 0.02)
    assert r4.phis == pytest.approx(math.asin(u0 * (zi / z) ** 2 / gain), abs=0.2)       # 安定な側の同期位相 (0 < φs < π/2)
    assert e['emittance_x'] == pytest.approx(exm, rel=tol_e)
    assert e['momentum_spread'] == pytest.approx(sd, rel=tol_s)


# ======================================================================= B234 (SOL 領域内の四極の本体の B∥)
def _sol_quad_setup(disfrin, p=0.3, bzT=3.0, el='QUAD E1 = (L=0.4 K1=0.3 DISFRIN=%s);'):
    from .test_v2_B1_sol import _run
    src = f'''MOMENTUM = {p} GEV;
{el % disfrin}
DRIFT D = (L=0.05);
SOL SA = (L=0 BZ={bzT} BOUND=1 GEO=1 DISFRIN=1) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA D E1 D SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _ = _run(src)
    return ffs


def _sol_refs(p, bzT, c, k=0.75):
    """Maxwell (ソフトなフリンジ、w = 4e-3・8e-3 m の w² 外挿) と硬い端 (B_s なし) の Thomas-BMT のスピン回転行列"""
    from .test_v2_B1_sol import CL, ME, _tbmt
    bz = bzT * CL / (p * 1e9)
    sharp = lambda x, y, s: (np.array([k * y, k * x, bz]) if 0.05 < s < 0.45 else np.array([0, 0, bz]))

    def fld(w):
        B = M.make_field(0.4, {1: k}, w, 2)
        return lambda x, y, s: B(x, y, s - 0.05) + np.array([0, 0, bz])
    ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
    R = lambda f: np.array([_tbmt(ck, s0, 0.5, f, p * 1e9 / ME)[6:] for s0 in np.eye(3)]).T
    return R(sharp), (4 * R(fld(4e-3)) - R(fld(8e-3))) / 3


@pytest.mark.parametrize('disfrin', ['0', '1'])
def test_sol_quad_body_bpar_tracking_and_emittance(v2, disfrin):
    """領域内の四極 (L = 0.4 m、K1 = 0.3、0.3 GeV、BZ = 3 T) のスピン (座標 3 倍の _P[1]): DISFRIN=0 (端のフリンジ + ψ あり) は Maxwell 場の積分と、
    DISFRIN=1 (端なし) は硬い端 (B_s なし) の積分と、追跡 1e-6・Emittance[] 5e-7 以内で一致する。以前 (本体の ∫B·dr = Δψ が無い) は
    DISFRIN=0 が Maxwell と 5.5e-6 (座標の 2 乗に比例)、DISFRIN=1 が硬い端と 5.5e-6 ずれていた (端の ψ だけで本体が無い = 和が −Δψ)"""
    from .test_v2_B1_sol import CL, ME, _P, _track
    ffs = _sol_quad_setup(disfrin)
    p, bzT = 0.3, 3.0
    c = _P[1].copy()
    c[:4] *= 3
    Rs, Rm = _sol_refs(p, bzT, c)
    ref = Rm if disfrin == '0' else Rs
    assert np.abs(Rs - Rm).max() > 5e-6
    cols = []
    for s0 in np.eye(3):
        _b, s = _track(ffs, c[None, :], 1, len(ffs.els) + 1, [s0])
        cols.append(s[0])
    assert np.abs(np.array(cols).T - ref).max() < 1e-6
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne(ffs.els, np.array(c, float), emit6d._Phys(p * 1e9, 100.0), rfsw=False, radcod=False)
        Re = st.srot[:, :3].copy()
        compat.set_spin_bpar(False)
        st0 = emit6d.tturne(ffs.els, np.array(c, float), emit6d._Phys(p * 1e9, 100.0), rfsw=False, radcod=False)
        Re0 = st0.srot[:, :3].copy()
    finally:
        compat.set_spin_bpar(None)
        emit6d.POL_MODE = False
    assert np.abs(Re - ref).max() < 5e-7
    assert np.abs(Re0 - ref).max() > 3e-6                # B165 を切ると (本体の Δψ なし = K64 の動作) 5.5e-6 ずれる


# ======================================================================= B238 (SOL 領域内の TRPT の加速、POL)
def test_sol_cavity_trpt_acceleration_spin_matches_integration(v2):
    """SOL 領域内の CAVI (L = 0.4 m、2 MV、0.02 GeV の電子で γ が 39.2 → 45.0、BZ = 0.1 T) の TRPT の加速 + POL (スピン追跡) を、一様な E_z (ソフトな端
    w = 2e-3 m の Maxwell 場: E_r = −r E'/2) と一様な B_z の中の運動方程式・Thomas-BMT (E×β の項も) の積分と比べる。以前は SADError (POL + TRPT + RFSW の
    SOL 領域内の加速は未対応)。差 8e-6 (γ の 15% の加速)。V = 0 の差は 4e-8 (積分の床)。加速に比例して増える残りの原因は未確認 (K64 の薄い空洞のモデルと
    厳密な E の場の違い。端のフリンジ tcavfrin の模型を含む)"""
    from scipy.integrate import solve_ivp
    from .test_v2_B1_sol import CL, ME, _run, _track, GS
    p, bzT, V, Lc = 0.02, 0.1, 3e6, 0.4
    src = f'''MOMENTUM = {p} GEV;
CAVI E1 = (L={Lc} VOLT={V} FREQ=1E6 PHI=-1.5707963267948966);
DRIFT D = (L=0.05);
SOL SA = (L=0 BZ={bzT} BOUND=1 GEO=1 DISFRIN=1) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA D E1 D SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _ = _run(src)
    ffs.flags['TRPT'] = True
    ffs.flags['RFSW'] = True
    bz = bzT * CL / (p * 1e9)
    u0 = p * 1e9 / ME
    eg = V / ME / Lc
    w = 2e-3
    Gk = lambda k, z: M._Gk(k, z - 0.05, Lc, w)

    def rhs(t, y):
        r, u, S = y[:3], y[3:6], y[6:].reshape(3, 3)
        g = math.sqrt(1 + u @ u)
        beta = u / g
        rr2 = r[0] ** 2 + r[1] ** 2
        E = np.array([-eg * Gk(1, r[2]) / 2 * r[0], -eg * Gk(1, r[2]) / 2 * r[1], eg * (Gk(0, r[2]) - rr2 * Gk(2, r[2]) / 4)])
        B = np.array([0, 0, bz * u0])
        Bpar = (B @ beta) / (beta @ beta) * beta
        Om = ((1 + GS * g) * (B - Bpar) + (1 + GS) * Bpar) / g - (GS + 1 / (g + 1)) * np.cross(beta, E)
        return np.concatenate([beta, E + np.cross(beta, B), np.cross(S, Om[None, :]).ravel()])
    worst = 0.0
    for c in (np.array([1e-3, 2e-4, -1e-3, 1e-4, 0, 1e-3]), np.array([-2e-3, -3e-4, 5e-4, 4e-4, 0, -5e-3])):
        ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0]]
        pt = 1 + c[5]
        uv = np.array([ck[1], ck[3], math.sqrt(pt * pt - ck[1] ** 2 - ck[3] ** 2)]) * u0
        ev = lambda t, y: y[2] - (0.1 + Lc)
        ev.terminal = True
        sol = solve_ivp(rhs, (0, 3), np.concatenate([[ck[0], ck[2], 0.0], uv, np.eye(3).ravel()]), events=ev, rtol=1e-11, atol=1e-14,
                        method='DOP853', max_step=w / 3)
        Sref = sol.y_events[0][0][6:].reshape(3, 3)
        for j, s0 in enumerate(np.eye(3)):
            _b, s = _track(ffs, c[None, :], 1, len(ffs.els) + 1, [s0])
            worst = max(worst, float(np.abs(s[0] - Sref[j]).max()))
    assert worst < 2e-5


# ======================================================================= B237 (Bend.transfer_matrix の K0 ≠ 0)
def _track_jac(els, p):
    h = 1e-6
    J = np.zeros((6, 6))
    f = lambda c: track6.track_particles(els, np.array([c], float), 1, p, 0.0, False, False, False, circ=50.0)[0]
    for j in range(6):
        d = np.zeros(6)
        d[j] = h
        J[:, j] = (f(d) - f(-d)) / (2 * h)
    return J


@pytest.mark.parametrize('decl,tol', [('BEND M1=(L=1 ANGLE=0.1 K0=0.004 E1=0.02 E2=0.03)', 2e-9),
                                      ('BEND M1=(L=1 ANGLE=-0.1 K0=0.01 E1=0.0 E2=0.0)', 2e-9),
                                      ('BEND M1=(L=1 ANGLE=0.1 K0=-0.01 K1=0.3 E1=0.05 E2=0.05)', 3e-5)])
def test_bend_transfer_matrix_with_k0_matches_tracking_jacobian(v2, decl, tol):
    """ANGLE ≠ 0 かつ K0 ≠ 0 の BEND の transfer_matrix (4×4 と δ の列) が追跡の写像の COD = 0 での微分と一致する。K0 なしの式 (h = ANGLE/L) は
    場の角 ANGLE + K0 を使わないので分散 (x の δ 列) が K0/ANGLE (K0 = 0.004・ANGLE = 0.1 で 4%) ずれた。K1 ≠ 0 は追跡の tbendi (数値積分の刻み) との差 8e-6"""
    src = f'MOMENTUM = 1 GEV;\n{decl};\nMARK ST=();\nLINE RING=(ST M1);\nFFS;\nUSE RING;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    els = [e for e in ffs.lat.elements if e.name == 'M1']
    J = _track_jac(els, ffs.momentum())
    Mt = np.array(els[0].transfer_matrix(0.0))
    assert np.abs(J[:4, :4] - Mt[:4, :4]).max() < tol
    assert np.abs(J[:4, 5] - Mt[:4, 5]).max() < tol
    compat.set_mode('k64')
    try:
        Mk = np.array(els[0].transfer_matrix(0.0))
    finally:
        compat.set_mode('v2')
    assert np.abs(J[:4, 5] - Mk[:4, 5]).max() > 1e-3          # k64 は K64 の (K0 なしの) 式のまま


# ======================================================================= B145 (F1 の fb1 ≠ fb2): 基準軌道の取り方
def test_f1_unequal_fringe_is_design_orbit_convention():
    """ANGLE = 0 の BEND (K0 = 0.003、L = 0.2 m、3 GeV) の F1 (磁場が直線で 0 → B0 に立ち上がる長さ fb1 = F1 + FB1、fb2 = F1 + FB2。立ち上がりの中心が要素の端)。
    x・px を厳密な Maxwell 場の運動 (`_p3_indep.maxwell_steer`) と比べる: px は 1e-10 で一致し、x の差は δ に依らない定数 (K0/al)(fb1² − fb2²)/24
    (fb1 = fb2 では 1e-11)。厳密な入口の変位は −dxfr1/(1+δ)、K64 は dxfr1 δ/(1+δ) (δ = 0 の粒子は動かさない) なので、差は設計軌道 (δ = 0) の変位 dxfr1 そのもの =
    『フリンジが設計軌道を動かさない』という基準軌道の取り方 (リングの幾何を決める硬い端の模型と一致させる)。誤りではないので直さない"""
    from . import _p3_indep as X
    al, K0 = 0.2, 0.003
    B0 = K0 / al
    diffs = []
    for fb1, fb2, d in ((0.04, 0.04, 0.0), (0.04, 0.04, 0.01), (0.1, 0.06, 0.02), (0.1, 0.06, -0.02), (0.3, 0.0001, 0.01)):
        x, px, _y, _py = X.maxwell_steer(fb1, fb2, al, K0, 0.0, d)
        xnom = x - (fb2 / 2) * px / math.sqrt((1 + d) ** 2 - px * px)
        c = np.zeros(6)
        c[5] = d
        k = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r F1=%r FB1=0 FB2=%r FRINGE=3);' % (al, K0, fb1, fb2 - fb1), c)
        assert k[1] == pytest.approx(px, abs=1e-10)
        diffs.append((k[0] - xnom) / (B0 * (fb1 ** 2 - fb2 ** 2) / 24) if fb1 != fb2 else k[0] - xnom)
    assert abs(diffs[0]) < 1e-10 and abs(diffs[1]) < 1e-10
    for r in diffs[2:]:
        assert r == pytest.approx(1.0, abs=2e-5)
