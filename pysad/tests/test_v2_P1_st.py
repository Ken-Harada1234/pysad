# -*- coding: utf-8 -*-
"""v2next 担当 P1 項目 2 (B116 の続き): n が傾く・非平面のリングの Polarization の緩和率。

v2.1.0 (A2) の緩和率 Σ r (b·n)² は、n と磁場の向き b が平行 (平面リング) で Sokolov-Ternov の厳密な式と同じだが、n が b から傾くと合わない。
厳密な式 (Baier-Katkov の反転の率): 1/τ_ST = (5√3/8) r_e ƛ_C c γ⁵ ⟨|ρ|⁻³ (1 − 2/9 (n·v)²)⟩ (v = 速度の向き)、
P_eq = −(8/5√3) ⟨ρ⁻³ b·(n − ∂n/∂δ)⟩ / ⟨|ρ|⁻³ (1 − 2/9 (n·v)² + 11/18 |∂n/∂δ|²)⟩ (Derbenev-Kondratenko)。反転は磁場が n に垂直 (垂直の偏向、水平磁場の挿入光源) でも起きる。
v2 (P1): 緩和率を Σ r (1 − 2/9 (n·v)²) (spin.ST_EXACT)、分子に運動論的偏極 b·∂n/∂δ (spin.KINETIC)、脱偏極率を放射の刻みごとのスピン軸の場で積分 (spin.STEP_RATE)。
独立な計算: pysad/tests/_p1_dk.py (閉軌道に沿った DK の式の数値積分。要素を刻み、軌道とスピンは追跡 (track6)、磁場は要素の定義から。放射積分・tradke・gintd・SLIM は使わない)。
追跡 (v2) は ROTATE のある BEND のスピンの座標系を回す (K64 は回さない。docs/history/v2next_P1.md §2.3)。
"""
import io
import math
import os

import sys

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, emit6d, ffscmd, spin, track6

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _p1_dk as DK

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
GS = 0.001159652181280002
ME = 0.51099895e6

_HEAD = '''MOMENTUM = %s GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DV=(L=0.3);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=QDK);
QUAD  QFE = (L=0.3 K1=0.7 DY=0.003) QDE = (L=0.3 K1=QDK DY=-0.003);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND  VB1 = (L=0.4 ANGLE=%g ROTATE=90 DEG) VB2=(L=0.4 ANGLE=%g ROTATE=90 DEG);
SOL SA = (L=0 BZ=%g BOUND=1) SB = (L=0 BZ=0 BOUND=1);
CAVI RFC = (VOLT=0.2e6 HARM=64);
MARK ST = ();
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLE1 = (QFE D2 B D2 QD D1);
LINE CELLE2 = (QF D2 B D2 QDE D1);
LINE VDOG = (VB1 DV VB2);
LINE SOLR = (SA DV DV SB);
'''


def _fodo(ring_line, E=0.5, vb=0.0, bz=0.0, qd=-0.6):
    src = _HEAD.replace('QDK', str(qd)) % (E, vb, -vb, bz) + ring_line + '\nFFS;\nUSE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return ffs


def _both(ffs):
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    sp = e['spin']
    dk = DK.DK(list(ffs.lat.elements), ffs.momentum())
    try:
        r = dk.run()
    finally:
        DK.restore_cave()
    return sp, r


def test_planar_ring_is_unchanged_by_the_exact_relaxation_rate_v2(v2_mode):
    """平面リング: Σ r (1 − 2/9 (n·v)²) = Σ r (b·n)² = Σ r (A2 と同じ)、P = 8/(5√3)、独立な DK の積分と τ が 1e-5 で合う"""
    ffs = _fodo('LINE RING = (ST 12*CELLQ RFC);')
    sp, r = _both(ffs)
    assert abs(sp['sdamp_rel'] / sp['sdamp_rel_b2'] - 1) < 1e-9
    assert abs(sp['taup'] / r['tau_st'] - 1) < 1e-5
    assert abs(sp['equpol'][0] - 8 / (5 * math.sqrt(3))) < 1e-7
    assert abs(r['P_dk'] - 8 / (5 * math.sqrt(3))) < 1e-9


def test_vertical_bends_flip_spins_although_b_is_perpendicular_to_n_v2(v2_mode):
    """垂直の偏向 (ROTATE=90° の BEND の組、軌道が鉛直方向に ±0.1 rad 曲がる、0.5 GeV の FODO): n は鉛直のまま、垂直の偏向の磁場は n に垂直 (b·n = 0)。
    Σ r (b·n)² (A2) は垂直の偏向の放射を数えず τ が 8% 長い。厳密な式は反転の率に数えて、DK の積分と τ が 1e-5、P が 2e-3 で合う"""
    ffs = _fodo('LINE RING = (ST 6*CELLQ VDOG 6*CELLQ RFC);', vb=0.1, qd=-0.55)
    sp, r = _both(ffs)
    assert abs(sp['taup'] / r['tau_st'] - 1) < 1e-4, (sp['taup'], r['tau_st'])
    assert abs(sp['equpol'][0] / r['P_dk'] - 1) < 2e-3
    assert sp['sdamp_rel'] / sp['sdamp_rel_b2'] > 1.05              # 垂直の偏向の分
    assert r['tau_b2'] / r['tau_st'] > 1.05                         # A2 の式なら τ が 8% 長い


def test_tilted_n_with_solenoid_relaxation_and_kinetic_polarization_v2(v2_mode):
    """ソレノイド (BZ = 1 T、0.8 m) のある 0.5 GeV の 12 セル FODO: n が傾く (n·v ≠ 0)。τ が DK の積分と 1e-5 (A2 の式は 3% 長い)、
    P が 3% 以内 (運動論的偏極 b·∂n/∂δ が分子の 1.7%。残る差は脱偏極率: Emittance[] の D と追跡の D の違い、記録 §5)"""
    ffs = _fodo('LINE RING = (ST SOLR 12*CELLQ RFC);', bz=1.0, qd=-0.55)
    sp, r = _both(ffs)
    assert abs(sp['taup'] / r['tau_st'] - 1) < 1e-4
    assert abs(sp['equpol'][0] / r['P_dk'] - 1) < 4e-2, (sp['equpol'][0], r['P_dk'])
    assert abs(r['kin_over_num']) > 0.01 and abs(sp['kinetic']) > 0.005 * abs(sp['sdamp'])
    assert r['tau_b2'] / r['tau_st'] > 1.02
    old_k, old_s = spin.KINETIC, spin.ST_EXACT
    spin.KINETIC = False
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
        assert abs(e['spin']['equpol'][0] / sp['equpol'][0] - 1) > 0.005
    finally:
        spin.KINETIC = old_k
    spin.ST_EXACT = False
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
        assert abs(e['spin']['taup'] / r['tau_st'] - 1) > 0.02          # A2 の緩和率 (b·n)² なし
    finally:
        spin.ST_EXACT = old_s


def test_vertical_dispersion_ring_relaxation_time_v2(v2_mode):
    """四極の縦ずれ (±3 mm の 2 台) で垂直の分散・閉軌道 (y 6 mm) のある 0.5 GeV の FODO: τ が DK の積分と 2e-3 (P は 3% 以内: 脱偏極率は
    Emittance[] の D と追跡の D が違うので (記録 §5) 独立な値が 3% の幅を持つ)"""
    ffs = _fodo('LINE RING = (ST 2*CELLQ CELLE1 4*CELLQ CELLE2 4*CELLQ RFC);')
    sp, r = _both(ffs)
    assert abs(sp['taup'] / r['tau_st'] - 1) < 2e-3
    assert abs(sp['equpol'][0] / r['P_dk'] - 1) < 0.03


@pytest.mark.parametrize('chi,ptol', [(0.3, 8e-3), (math.pi / 2, 4e-3)])
def test_rotated_insertion_device_polarization_v2(v2_mode, chi, ptol):
    """ASTRID2 + MPW (0.58 GeV、11 周期の K0 の極) の MPW を ROTATE = 0.3 rad / 90° (水平磁場の挿入光源 = 垂直の波打ち) にしたリング:
    τ が DK の積分と 2e-3、P が 0.8% / 0.4% で合う。90° では磁場が n に垂直で、A2 の Σ r (b·n)² は MPW の放射を数えず τ が 32% 短い (反転は起きる)"""
    path = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
    if not os.path.exists(path):
        pytest.skip('lattice がない')
    src = open(path, encoding='utf-8').read() + '\nFFS;\nUSE RINGMPW;\nCELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    for e in ffs.lat.elements:
        if e.name.startswith('BMPW'):
            e.ROTATE = chi
    ffs.calc(quiet=True)
    sp, r = _both(ffs)
    assert abs(sp['taup'] / r['tau_st'] - 1) < 2e-3, (sp['taup'], r['tau_st'])
    assert abs(sp['equpol'][0] / r['P_dk'] - 1) < ptol, (sp['equpol'][0], r['P_dk'])
    if chi > 1.0:
        assert r['tau_b2'] / r['tau_st'] > 1.2


def test_k64_mode_keeps_the_signed_relaxation_rate(k64_mode):
    """k64: 緩和率は K64 のまま (符号つきの sdamp。sdamp_rel_b2 や運動論的偏極は出ない)。垂直の偏向のリングで τ は垂直の偏向の放射を数えない"""
    ffs = _fodo('LINE RING = (ST 6*CELLQ VDOG 6*CELLQ RFC);', vb=0.1, qd=-0.55)
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    sp = e['spin']
    assert 'sdamp_rel_b2' not in sp and 'kinetic' not in sp
    dk = DK.DK(list(ffs.lat.elements), ffs.momentum())
    try:
        r = dk.run()
    finally:
        DK.restore_cave()
    assert r['tau_st'] / sp['taup'] < 0.95         # k64 の τ (符号つきの sdamp から) は垂直の偏向の分だけ長い


# ---------------------------------------------------------------- 追跡の ROTATE (v2: スピンの座標系を回す)
def _tbmt(c, s, L, field, p0mc):
    pt = 1 + c[5]
    pv = np.array([c[1], c[3], math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)])
    gam = math.sqrt(1 + (p0mc * pt) ** 2)

    def rhs(l, u):
        r, pp, S = u[:3], u[3:6], u[6:]
        B = field(r[0], r[1], r[2])
        ph = pp / pt
        Bl = np.dot(B, ph) * ph
        Om = ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt
        return np.concatenate([ph, np.cross(pp, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: u[2] - L
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2 * L + 1), np.concatenate([[c[0], c[2], 0.0], pv, s]), events=ev, rtol=1e-12, atol=1e-15, method='DOP853')
    return sol.y_events[0][0]


def _rotated_k0_bend_spin(chi, c0):
    k0, L, E = 0.0640435, 0.0875, 0.58
    src = f'MOMENTUM = {E} GEV;\nDRIFT D1=(L=1);\nBEND E1 = (L={L} K0={k0} ROTATE={chi});\nMARK ST = ();\nLINE RING = (ST E1 D1);\nFFS;\nUSE RING;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    els = [e for e in ffs.lat.elements if e.name == 'E1']
    p = ffs.momentum()
    field = lambda x, y, s: np.array([k0 / L * math.sin(chi), k0 / L * math.cos(chi), 0.0])
    Rt = np.zeros((3, 3))
    Rx = np.zeros((3, 3))
    for j in range(3):
        s0 = np.zeros(3)
        s0[j] = 1.0
        _o, sp = track6.track_particles(els, np.array([c0], float), 1, p, 0.0, False, True, False, spins=[s0], radpol=False, circ=100.0)
        Rt[:, j] = sp[0]
        Rx[:, j] = _tbmt(c0, s0, L, field, p / ME)[6:]
    # 厳密な積分は実験室系 (回さない座標) なので、軌道が軸上 (c0 = 0) のときだけ座標系の違いが出ない
    return Rt, Rx


def test_tracking_spin_of_rotated_bend_precesses_about_the_rotated_axis_v2(v2_mode):
    """ROTATE = 90° の K0 の BEND (磁場は水平): 追跡のスピンは水平軸 (x) のまわりに歳差する (厳密な Thomas-BMT の積分と 1e-9。軸上)"""
    Rt, Rx = _rotated_k0_bend_spin(math.pi / 2, np.zeros(6))
    assert np.abs(Rt - Rx).max() < 1e-8, np.abs(Rt - Rx).max()
    assert abs(Rt[1, 2] - 0.1479) < 1e-3                # x のまわりの回転 (角 (1+aγ) K0 − K0 = aγ K0)


def test_tracking_spin_of_rotated_bend_k64_precesses_about_the_unrotated_axis(k64_mode):
    """k64 (K64 の tbshift・tent は座標だけ回してスピンを回さない): 同じ BEND でスピンが鉛直軸のまわりに歳差する (厳密な値と 0.1 以上違う)"""
    Rt, Rx = _rotated_k0_bend_spin(math.pi / 2, np.zeros(6))
    assert np.abs(Rt - Rx).max() > 0.1
    assert abs(Rt[0, 2]) > 0.14                          # y のまわりの回転
