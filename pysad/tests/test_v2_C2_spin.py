# -*- coding: utf-8 -*-
"""担当 C2: POL + RAD のスピン追跡 (track6 の経路) での SOL 領域 (端のずれ・回転・F1、領域内の DRIFT・QUAD・BEND・MULT・CAVI)、
UND、PHSROT。v1 と k64 モードは、領域内は DRIFT・QUAD だけで、UND・PHSROT は NotImplementedError だった。

根拠 (gs を使わない):
  - Thomas-BMT の式 dS/dl = S × [(1 + aγ) B⊥ + (1 + a) B∥]/|p| と運動方程式 dp/dl = p × B/|p| を、同じ一様磁場・四極磁場の中で
    数値積分 (solve_ivp) したスピンと一致。
  - 偏向磁石: スピンと運動量のなす角 = aγ × 偏向角 (ν_s = aγ)。アンジュレータの 1 極でも同じ。
  - 平面リングにソレノイド 1 台: 1 周のスピンの回転角 μ は cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1 + a) ∫Bz dl/(Bρ)。
  - 軌道 (座標) は高速追跡 (track6n) と同じ (FLUC OFF)。

K64 の式との違い (v2 で直した。記録は docs/history/v2_changes_C2.md):
  - sprot の縦磁場の項 bsi = bz·al は p0 で規格化したまま z 軸の向きに足される: δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならず、
    軌道が傾いていると運動量に垂直な成分が二重に入る。双極 + ソレノイド (tdrift) の bsi は符号・大きさとも合わない (3.6e-3 rad の差)。
  - tsolqur (領域内の四極、RAD): 縦磁場の積分が 1 スライス分足りない (1 スライスなら半分)。放射の基準の向きが正準運動量のままで、
    軸を外れた粒子の損失が大きく出る (試験の粒子で 2.2 倍)。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat

GS = 0.001159652181280002            # a = (g − 2)/2
ME = 0.51099895e6
CL = 299792458.0


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, I = run_sad_string(src, out=buf)
    return ffs


def _track(ffs, c, k0, k1, spin=None, rad=True):
    """spin が None なら高速追跡 (POL OFF)、あれば POL ON (track6 の経路)。戻り値 (座標 (n, 6), スピン (n, 3) か None)"""
    c = np.asarray(c, float)
    n = len(c)
    cols = [list(c[:, j]) for j in range(6)]
    ffs.flags['RAD'] = rad
    ffs.flags['FLUC'] = False
    if spin is None:
        ffs.flags['POL'] = False
        out = ffs.fn_track([[float(k0), cols + [[1.0] * n]], float(k1)])
        return np.array(out[1][:6]).T, None
    ffs.flags['POL'] = True
    try:
        sp = np.asarray(spin, float)
        out = ffs.fn_track([[float(k0), cols + [[math.atan2(s[2], s[0]) for s in sp], [s[1] for s in sp], [1.0] * n]], float(k1)])
    finally:
        ffs.flags['POL'] = False
    r = out[1]
    so = [[math.sqrt(max(0.0, 1 - y * y)) * math.cos(a), y, math.sqrt(max(0.0, 1 - y * y)) * math.sin(a)] for a, y in zip(r[6], r[7])]
    return np.array(r[:6]).T, np.array(so)


def _tbmt(c, s, L, field, p0mc):
    """運動方程式と Thomas-BMT の数値積分。c = (x, 力学的 px, y, 力学的 py, z, δ)、field(x, y, s) → B (1/m、p0 で規格化)"""
    pt = 1 + c[5]
    pv = np.array([c[1], c[3], math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)])
    gam = math.sqrt(1 + (p0mc * pt) ** 2)

    def rhs(l, u):
        r, p, S = u[:3], u[3:6], u[6:]
        B = field(r[0], r[1], r[2])
        ph = p / pt
        Bl = np.dot(B, ph) * ph
        Om = ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt
        return np.concatenate([ph, np.cross(p, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: u[2] - L
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2 * L + 1), np.concatenate([[c[0], c[2], 0.0], pv, s]), events=ev, rtol=1e-12, atol=1e-15, method='DOP853')
    return sol.y_events[0][0]


_BASE = '''MOMENTUM = 1 GEV;
%s
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1 DISFRIN=1) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA E1 SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
_P = np.array([[0, 0, 0, 0, 0, 0], [1e-3, 2e-4, -1e-3, 1e-4, 0, 1e-3], [-2e-3, -3e-4, 5e-4, 4e-4, 0, -5e-3]], float)
_S0 = ([0, 0, 1], [1, 0, 0], [0, 1, 0])
_CASES = {
    'drift': ('DRIFT E1 = (L=0.4);', 1.0, lambda bz: (lambda x, y, s: np.array([0.0, 0.0, bz])), 5e-7),
    'bend': ('BEND E1 = (L=0.4 K0=0.02 ROTATE=0.3);', 0.0,
             lambda bz: (lambda x, y, s: np.array([0.05 * math.sin(0.3), 0.05 * math.cos(0.3), bz])), 3e-6),
    'bend_sol': ('BEND E1 = (L=0.4 K0=0.02 ROTATE=0.3);', 1.0,
                 lambda bz: (lambda x, y, s: np.array([0.05 * math.sin(0.3), 0.05 * math.cos(0.3), bz])), 1e-5),
    'quad_sol': ('QUAD E1 = (L=0.4 K1=0.3 DISFRIN=1);', 1.0, lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.75 * x, bz])), 3e-6),
    'mult_sol': ('MULT E1 = (L=0.4 K0=0.02 K1=0.3 DISFRIN=1);', 1.0,
                 lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.05 + 0.75 * x, bz])), 5e-5),
}


@pytest.mark.parametrize('case', list(_CASES))
def test_spin_in_sol_region_matches_thomas_bmt_integration(v2_mode, case):
    """領域内の DRIFT・BEND (bz = 0 と 1 T)・QUAD・MULT (K0 + K1) の中のスピン。1 GeV (aγ = 2.27)、0.4 m、軸上と軸を外れた粒子
    (δ = 1e-3、−5e-3)、3 つの向きのスピン。残る差は放射の刻みの粗さ (歳差を刻みごとの 1 回の回転にまとめる) と放射による δ の変化"""
    el, bzT, fld, tol = _CASES[case]
    ffs = _run(_BASE % (el, bzT))
    bz = bzT * CL / 1e9
    f = fld(bz)
    worst = 0.0
    for s0 in _S0:
        b, s = _track(ffs, _P, 1, len(ffs.els) + 1, spin=[s0] * len(_P))
        for c, sv in zip(_P, s):
            ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
            u = _tbmt(ck, s0, 0.4, f, 1e9 / ME)
            worst = max(worst, float(np.abs(sv - u[6:]).max()))
            assert abs(np.linalg.norm(sv) - 1.0) < 1e-9
    assert worst < tol, (case, worst)


def test_spin_precession_in_dipole_is_a_gamma_times_bend_angle(v2_mode):
    """領域内の BEND (bz = 0): スピンと運動量のなす角 = aγ × 偏向角 (ν_s = aγ)。偏向角 0.02 と 0.3 rad"""
    for el in ('BEND E1 = (L=0.4 K0=0.02);', 'BEND E1 = (L=0.4 K0=-0.3);'):
        ffs = _run(_BASE % (el, 0.0))
        b, s = _track(ffs, np.zeros((1, 6)), 1, len(ffs.els) + 1, spin=[[0, 0, 1]])
        th = -b[0, 1]                                        # 偏向角 (px の変化)
        assert abs(th) > 5e-3
        gam = math.sqrt(1 + (1e9 / ME) ** 2)
        ang_s = math.atan2(s[0, 0], s[0, 2])                 # スピンの向き (z から x へ)
        ang_p = math.atan2(b[0, 1], math.sqrt(1 - b[0, 1] ** 2))
        assert abs((ang_s - ang_p) / (GS * gam * (ang_p)) - 1) < 2e-4, (el, ang_s, ang_p)


@pytest.mark.parametrize('pole,y0', [(0.25, 0.0), (0.25, 2e-3), (5, 2e-3), (5.25, 3e-3)])
def test_spin_in_undulator_matches_thomas_bmt_integration(v2_mode, pole, y0):
    """UND (FBX = 1 T、1 GeV、0.5 m): 磁場 B = rot A (By = F0 cosh(kz y) cos(kz s)、Bs = −F0 sinh(kz y) sin(kz s)) の中の
    Thomas-BMT の積分と比べる。周期の途中で終わる UND (POLE = 0.25、5.25) の出口はポテンシャルの飛びが薄い蹴りになるので、
    スピンもその分 (1 + aγ) × (向きの変化) だけ回す。中点面を外れた粒子 (y0≠0) は縦磁場 Bs の歳差 (1 + a) ∫Bs ds も入る"""
    L = 0.5
    ffs = _run('MOMENTUM = 1 GEV;\nUND E1 = (L=%s FBX=1.0 POLE=%s SLICE=2000);\nMARK ST = (BX=5 BY=5);\nLINE BL = (ST E1);\nFFS;\nUSE BL;\n'
               'INS; NORAD; NOFLUC; CALC;\nExit;\n' % (L, pole))
    F0 = 1.0 * CL / 1e9
    kz = 2 * math.pi * pole / L
    fld = lambda x, y, s: np.array([0.0, F0 * math.cosh(kz * y) * math.cos(kz * s), -F0 * math.sinh(kz * y) * math.sin(kz * s)])
    c = np.array([[1e-3, 1e-4, y0, 2e-4, 0, 1e-3]])
    pt = 1 + c[0, 5]
    gam = math.sqrt(1 + (1e9 / ME * pt) ** 2)
    for s0 in _S0:
        b, s = _track(ffs, c, 1, 3, spin=[s0])
        u = _tbmt(c[0], s0, L, fld, 1e9 / ME)
        ax = F0 / kz * math.cosh(kz * u[1]) * math.sin(kz * L)          # 出口のポテンシャル
        p0v = u[3:6] / pt
        p1 = np.array([p0v[0] + ax / pt, p0v[1], 0.0])
        p1[2] = math.sqrt(1 - p1[0] ** 2 - p1[1] ** 2)
        axis = np.cross(p0v, p1)
        sn = np.linalg.norm(axis)
        S = u[6:]
        if sn > 1e-14:
            axis = axis / sn
            a = (1 + GS * gam) * math.asin(sn)
            S = S * math.cos(a) + np.cross(axis, S) * math.sin(a) + axis * np.dot(axis, S) * (1 - math.cos(a))
        assert np.abs(s[0] - S).max() < 2e-6, (pole, y0, s0)


_RING = '''MOMENTUM = 0.5 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA DS DS SB);
LINE RINGQ = (ORIGIN 12*CELLQ CELLS 11*CELLQ);
FFS;
USE RINGQ;
CELL; CALC;
Exit;
'''


@pytest.mark.parametrize('bzT', [0.0, 1.0, 2.5])
def test_ring_spin_tune_with_solenoid(v2_mode, bzT):
    """平面リング (0.5 GeV、aγ = 1.135) の 1 周のスピンの回転: BZ = 0 なら回転角 2π aγ (ν_s = aγ)、長さ 1 m のソレノイドを入れると
    cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1 + a) bz L。閉軌道 (設計軌道) の上の粒子を 3 つの向きのスピンで 1 周追う"""
    ffs = _run(_RING % bzT)
    n = len(ffs.els) + 1
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        b, s = _track(ffs, np.zeros((1, 6)), 1, n, spin=[s0])
        M[:, j] = s[0]
    assert np.max(np.abs(M @ M.T - np.eye(3))) < 1e-8
    cmu = (np.trace(M) - 1.0) / 2.0
    gam = math.sqrt(1 + (0.5e9 / ME) ** 2)
    phi = (1 + GS) * (bzT * CL / 0.5e9) * 1.0
    half = math.cos(math.pi * GS * gam) * math.cos(phi / 2)
    assert abs(cmu - (2 * half * half - 1)) < 2e-5, (cmu, 2 * half * half - 1)


_ALL = '''MOMENTUM = 3 GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
BEND BK = (L=0.4 K0=0.002 ROTATE=0.3);
MULT MK = (L=0.3 K0=-0.001 K1=0.2 K2=3 ROTATE=0.2);
QUAD QF = (L=0.3 K1=0.3);
CAVI CA = (L=0.3 VOLT=1E6 FREQ=5E8 PHI=0.3);
UND U1 = (L=0.5 FBX=1.0 POLE=5 SLICE=200);
PHSROT PR = (AX=0.5 BX=4 PX=0.31 AY=-0.3 BY=2 PY=0.22 AZ=0 BZ=10 PZ=0.01 D11=1e-4 D22=1e-4 D33=2e-4 D44=2e-4 D55=1e-4 D66=1e-4);
SOL SA = (L=0 BZ=3 BOUND=1 GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002 F1=0.1) SM=(L=0 BZ=1.5 F1=0.05) SB = (L=0 BZ=0 BOUND=1 F1=0.1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST DR SA DR2 BK DR2 SM MK DR2 QF CA SB DR U1 DR PR);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
_P4 = np.array([[0, 0, 0, 0, 0, 0], [1e-3, 2e-4, -1e-3, 1e-4, 1e-3, 1e-3], [-2e-3, -3e-4, 5e-4, 4e-4, -2e-3, -5e-3],
                [3e-3, 1e-3, 2e-3, -1e-3, 0, 2e-2]], float)


@pytest.mark.parametrize('rfsw', [False, True])
def test_spin_path_follows_the_same_orbit_as_fast_tracking(v2_mode, rfsw):
    """track6 の経路 (POL + RAD) の座標は高速追跡 (track6n) と同じ: ずれて傾いた SOL (F1 つき)、領域内の BEND・MULT・QUAD・CAVI・SOL、
    UND (磁場つき)、PHSROT (減衰つき)。スピンの大きさは 1 のまま"""
    ffs = _run(_ALL)
    ffs.flags['RFSW'] = rfsw
    n = len(ffs.els) + 1
    a, _ = _track(ffs, _P4, 1, n)
    b, s = _track(ffs, _P4, 1, n, spin=[[0, 0, 1]] * len(_P4))
    assert np.max(np.abs(a - b)) < 1e-10
    assert np.max(np.abs(np.linalg.norm(s, axis=1) - 1.0)) < 1e-9
    assert abs(s[0, 2] - 1.0) > 1e-4                         # スピンは動いている


def test_spin_frame_rotation_at_tilted_solenoid_is_the_geometric_rotation(v2_mode):
    """SOL の端のスピンの回転 (trots の R) は、幾何の軸の向きの変化そのもの: R = G1ᵀ G0 (G は x・y・s の軸を列に持つ行列)。
    磁場のない傾いた領域を通ると、スピンは元に戻る (追跡)"""
    from pysad import geo as G, sgeo, emit6d as _e
    src = _ALL.replace('BZ=3 ', 'BZ=%s ').replace('BZ=1.5 ', 'BZ=0 ').replace(
        'LINE BL = (ST DR SA DR2 BK DR2 SM MK DR2 QF CA SB DR U1 DR PR);', 'LINE BL = (ST DR SA DR2 DR2 SB DR);')
    ffs = _run(src % '3')
    i = [k for k, e in enumerate(ffs.els) if e.name == 'SA'][0]
    ke = [k for k, e in enumerate(ffs.els) if e.name == 'SB'][0]
    ent, ext = sgeo.region_frames(ffs.els, i, _e._Phys(ffs.momentum(), 10.0))
    res = G.compute_geo(ffs.els, warn=lambda s: None, momentum_eV=ffs.momentum())
    g0, g1, g2, g3 = res.geo[i][:, :3], res.geo[i + 1][:, :3], res.geo[ke][:, :3], res.geo[ke + 1][:, :3]
    assert np.max(np.abs(ent[0] - g1.T @ g0)) < 1e-12
    assert np.max(np.abs(ext[0] - g3.T @ g2)) < 1e-12
    ffs = _run(src % '0')
    for s0 in ([0, 0, 1], [1, 0, 0], [0.6, 0.8, 0]):
        b, s = _track(ffs, _P4[:2], 1, len(ffs.els) + 1, spin=[s0] * 2)
        assert np.max(np.abs(s - np.array(s0))) < 1e-9


def test_k64_mode_spin_path_is_unchanged(k64_mode):
    """k64 モード: 領域内の BEND・UND・PHSROT は従来どおり未対応。領域内の DRIFT・QUAD のスピンは K64 の式のまま"""
    ffs = _run(_BASE % ('DRIFT E1 = (L=0.4);', 1.0))
    b, s = _track(ffs, _P, 1, len(ffs.els) + 1, spin=[[1, 0, 0]] * len(_P))
    bz = CL / 1e9
    u = _tbmt([_P[2, 0], _P[2, 1] + 0.5 * bz * _P[2, 2], _P[2, 2], _P[2, 3] - 0.5 * bz * _P[2, 0], 0, _P[2, 5]], [1, 0, 0], 0.4,
              lambda x, y, z: np.array([0.0, 0.0, bz]), 1e9 / ME)
    assert 3e-4 < np.abs(s[2] - u[6:]).max() < 1e-3          # K64 の式: δ = −5e-3 の粒子で 6e-4 の差 (v2 は 1e-7)
    for el in ('BEND E1 = (L=0.4 K0=0.02);', ):
        ffs = _run(_BASE % (el, 1.0)) if False else None
    with pytest.raises(Exception):
        ffs = _run('MOMENTUM = 1 GEV;\nUND E1 = (L=0.25 POLE=1 SLICE=40);\nMARK ST = (BX=5 BY=5);\nLINE BL = (ST E1);\nFFS;\nUSE BL;\n'
                   'INS; NORAD; NOFLUC; CALC;\nExit;\n')
        _track(ffs, np.zeros((1, 6)), 1, 3, spin=[[0, 0, 1]])
