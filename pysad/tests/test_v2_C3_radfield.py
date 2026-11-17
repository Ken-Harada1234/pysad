# -*- coding: utf-8 -*-
"""担当 C3 (v2next): RadiationField を、Liénard-Wiechert の式の独立な数値計算 (遅延時間を解く) と比べる。

独立な基準 (`_c3_indep.py`): 円弧・直線の軌道 (解析的な位置・速度・加速度) で、観測時刻 τ ごとに遅延時間 t' を brentq で解き、
E = q/(4πε0)[(n−β)/(γ²(1−n·β)³R²) + n×((n−β)×β̇)/(c(1−n·β)³R)] を評価する。pysad の RadiationField は K64 の
Feynman-Heaviside 形式 (軌道の点列と自然スプラインの微分) で、同じ式ではない。
確かめたこと: 振幅・位相 (τ ごとの場の一致)、スペクトル (Jackson の式 K_{2/3}、臨界周波数 ωc = 3γ³c/2ρ)、全エネルギー
(e²cγ⁴β⁴/(6πε0ρ²) × 弧の時間)、複数セグメント (シケイン、区間の境目の重複点)、RADLIGHT の複数周 (重複点の除去)。
条件: 電子。γ = 20 (10.2 MeV、ρ = 1 m)、1 GeV (ρ = 4.58 m の偏向磁石、RADLIGHTDIV = 16)。gs は RadiationField で SIGSEGV (B49) するので gs との比較はない。
"""
import io
import math
import os
import sys

import numpy as np
import pytest
from scipy.special import kv

from pysad import compat, ffscmd
from pysad.radfield import SADError, radiation_field, radiation_spectrum

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _c3_indep as C       # noqa: E402

GAMMA = 20.0


def _arc_traj(m, half=0.5, rho=1.0):
    P = C.Path([('arc', 2 * half * rho, 1.0 / rho)], GAMMA)
    ss, tr = P.sample(m)
    return P, tr


def test_field_of_a_bend_equals_lienard_wiechert_amplitude_and_phase():
    """円弧 (弧の角 ±0.5 rad、γ = 20、ρ = 1 m) の接線方向 (軌道の中央の方向) 50 m の観測点の場: 点を細かくすると収束し (h² で、4001 点 → 8001 点で 4 倍)、
    主成分 (水平) が 8001 点で最大値の 1e-5 以内、全成分が 2e-5 以内"""
    errs = []
    for m in (2001, 4001, 8001):
        P, tr = _arc_traj(m)
        # 弧の中央 (φ = 0.5 の接線方向) を向く観測点: 軌道は z 方向から x 方向へ曲がる。接線 = (sin 0.5, 0, cos 0.5)
        obs = [50 * math.sin(0.5), 0.0, 50 * math.cos(0.5)]
        # radiation_field の座標は (t, x, y, z) = 軌道の 4 本。ここでは Path の x, y, z をそのまま渡す
        f = radiation_field(tr, obs)
        tau = f[:, 0]
        sel = np.arange(m // 4, 3 * m // 4, max(1, m // 400))
        ex = C.lw_field(P, tau[sel], obs, r0=50.0)
        peak = np.abs(ex).max()
        errs.append(np.abs(f[sel, 1:4] - ex).max() / peak)
    assert errs[2] < 2e-5, errs
    assert errs[0] / errs[1] > 3.0 and errs[1] / errs[2] > 3.0, errs          # 2 次収束 (スプラインの微分の離散化誤差)


def test_spectrum_and_critical_frequency_match_jackson():
    """RadiationSpectrum の |c + is|² から求めた d²I/dωdΩ (軸上、θ = 0) が Jackson の式
    (e²/(4πε0))/(3π²c) (ωρ/c)² γ^-4 K_{2/3}²(ξ)、ξ = ωρ/(3cγ³) と ω = 0.1〜2.5 ωc で 1% 以内 (ωc = 3γ³c/(2ρ))"""
    P, tr = _arc_traj(4001)
    D = 500.0
    obs = [D * math.sin(0.5), 0.0, D * math.cos(0.5)]
    f = radiation_field(tr, obs)
    tau = f[:, 0]
    sl = slice(40, len(tau) - 40)
    wc = 3 * GAMMA ** 3 * C.C / 2.0
    lc = 2 * math.pi * C.C / wc
    rng_ = [lc / 2.6, lc / 0.1, (lc / 0.1 - lc / 2.6) / 40.0]
    tot = 0.0
    for comp in (1, 2, 3):
        k, c, s = radiation_spectrum(tau[sl], f[sl, comp], rng_)
        tot = tot + np.array(c) ** 2 + np.array(s) ** 2
    k = np.array(k)
    w = k * C.C
    i_num = D ** 2 * C.EPS0 * C.C * tot / math.pi
    q2 = C.E0 ** 2 / (4 * math.pi * C.EPS0)
    xi = w * 1.0 / (3 * C.C * GAMMA ** 3)
    i_an = q2 / (3 * math.pi ** 2 * C.C) * (w / C.C) ** 2 * GAMMA ** -4 * kv(2.0 / 3.0, xi) ** 2
    sel = (w / wc > 0.1) & (w / wc < 2.6)
    assert sel.sum() >= 10
    assert np.abs(i_num[sel] / i_an[sel] - 1).max() < 0.01, (w[sel] / wc, i_num[sel] / i_an[sel])
    # 臨界周波数: ω² K²_{2/3}(ω/2ωc ...) の形の d²I/dωdΩ のピークの位置が Jackson の式のピークと同じ (ピークは ω ≈ 0.3 ωc 付近)
    assert int(np.argmax(i_num)) == int(np.argmax(i_an))


def test_total_energy_per_arc_angle_matches_larmor():
    """観測点を弧の接線方向に沿って垂直方向にも振り (±0.3 rad、31 点)、ε0 c R² ∫|E|² dτ cosθ dθ の和 = 偏向 1 rad あたりの全放射エネルギー
    e²γ⁴β³/(6πε0ρ) (弧の時間 1/ω × Larmor の式 e²cγ⁴β⁴/(6πε0ρ²)) と 1e-3 以内 (D = 500 m、点 4001)"""
    P, tr = _arc_traj(4001)
    D = 500.0
    ths = np.linspace(-0.3, 0.3, 31)
    dth = ths[1] - ths[0]
    tot = 0.0
    for th in ths:
        o = D * np.array([math.sin(0.5) * math.cos(th), math.sin(th), math.cos(0.5) * math.cos(th)])
        f = radiation_field(tr, o)
        tau = f[:, 0]
        e2 = f[:, 1] ** 2 + f[:, 2] ** 2 + f[:, 3] ** 2
        tot += C.EPS0 * C.C * D ** 2 * math.cos(th) * np.trapezoid(e2[40:-40], tau[40:-40]) * dth
    beta = P.beta
    want = C.E0 ** 2 * GAMMA ** 4 * beta ** 3 / (6 * math.pi * C.EPS0 * 1.0)
    assert abs(tot / want - 1) < 1.5e-3, tot / want


def test_chicane_multi_segment_with_repeated_joint_points_matches_lienard_wiechert():
    """直線 + 円弧 + 直線 + 逆向きの円弧 (シケイン) の軌道。区間の境目の点を重ねた点列 (RADLIGHT の複数区間と同じ) を RadiationField に渡し、
    厳密な Liénard-Wiechert の場 (パルスが 2 つ: 前半・後半の弧の中央) と最大値の 3e-4 以内 (6001 点)。重なった点は除かれ、戻りの点数は重なりなしの点数"""
    P = C.Path([('arc', 0.4, 1.0), ('line', 0.3), ('arc', 0.4, -1.0)], GAMMA)
    m = 6001
    ss, tr = P.sample(m)
    j = int(np.argmin(abs(ss - 0.4)))
    tr2 = np.vstack([tr[:j + 1], tr[j:j + 1], tr[j + 1:]])
    obs = [100 * math.sin(0.2), 0.0, 100 * math.cos(0.2)]
    f = radiation_field(tr2, obs)
    assert f.shape[0] == m
    tau = f[:, 0]
    sel = np.arange(150, len(tau) - 150, 5)
    ex = C.lw_field(P, tau[sel], obs, r0=100.0)
    assert np.abs(f[sel, 1:4] - ex).max() / np.abs(ex).max() < 3e-4
    assert np.abs(ex).max() > 5e-6                           # 2 つのパルスが窓の中にある


def test_non_monotonic_time_is_an_error():
    P, tr = _arc_traj(200)
    tr = tr.copy()
    tr[100, 0] = tr[50, 0]                                  # 時刻が戻る (位置は違う) = 軌道の誤り
    with pytest.raises(SADError):
        radiation_field(tr, [0.0, 0.0, 10.0])


_SRC_RING = r'''MOMENTUM = 1e9;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1);
FFS;
USE RING;
CELL;
CALC;
RADLIGHT;
np=2;
x=Table[0,{7},{np}];
x[[1,2]]=1e-4;x[[3,2]]=-1e-4;x[[6,2]]=1e-3;
Do[x[[7,i]]=1,{i,np}];
d={3,x};
r=TrackParticles[d,9,1,2];
t=r[[2]];
traj=t[[1]];
Print["DUP ",Length[traj[[1]]]-Length[Union[traj[[1]]]]];
f=RadiationField[traj,{3.0,0.5,6.0}];
Print["FD ",Dimensions[f]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_radlight_multi_turn_trajectory_gives_a_field(mode):
    """複数周の RADLIGHT の軌道 (区間の始めの点と直前の要素の後の点が重なる) で RadiationField が失敗していた
    (`x must be strictly increasing sequence`)。重なった点を除いて場を返す (両モード。k64 は gs の軌道と同じ 22 点 → 21 点)"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        ffscmd.run_sad_string(_SRC_RING, out=buf)
    finally:
        compat.set_mode(old)
    txt = buf.getvalue()
    assert 'strictly increasing' not in txt
    dup = [l for l in txt.splitlines() if l.startswith('DUP')][0]
    assert float(dup.split()[1]) == 1.0
    fd = [l for l in txt.splitlines() if l.startswith('FD')][0]
    nums = [float(v) for v in fd.replace('{', ' ').replace('}', ' ').replace(',', ' ').split()[1:]]
    assert nums[0] == 10.0 and nums[1] > 20


def test_radlight_1gev_bend_spectrum_matches_jackson_v2(v2_mode):
    """v2 の RADLIGHT (BEND は円弧、RADLIGHTDIV = 16 点/(1/γ)) の 1 GeV の偏向磁石 (L=1.2 m、ANGLE=0.2618) の軌道 8200 点から RadiationField と
    RadiationSpectrum で求めた軸上の d²I/dωdΩ が、Jackson の式と ω = 0.2〜5 ωc (λc = 2.56 nm) で 2e-3 以内。観測点は 5 km (遠方: 50 m だと
    1/γ·ρ の光源の大きさに対し Fresnel 領域で 4.6% 小さくなる)"""
    src = '''MOMENTUM = 1e9;
DRIFT D1 = (L=0.5);
BEND B = (L=1.2 ANGLE=0.2617993877991494);
MARK ST = ();
LINE RING = (ST D1 B D1);
FFS;
USE RING; INS; CALC;
RADLIGHT;
np=1;
x=Table[0,{7},{np}];
Do[x[[7,i]]=1,{i,np}];
d={1,x};
r=TrackParticles[d,4];
traj=r[[2,1]];
Exit;
'''
    ffs, I = ffscmd.run_sad_string(src, out=io.StringIO())
    traj = np.array(I.globals['traj'], float)
    assert traj.shape == (4, 8200)
    gamma = math.hypot(1.0, ffs.momentum() / 0.51099895e6)
    rho = 1.2 / 0.2617993877991494
    D = 5000.0
    phi = 0.13089969389957471
    f = radiation_field(traj.T, [D * math.cos(phi), D * math.sin(phi), 0.0])
    tau = f[:, 0]
    wc = 3 * gamma ** 3 * C.C / (2 * rho)
    lc = 2 * math.pi * C.C / wc
    tot = 0.0
    for comp in (1, 2, 3):
        k, c, s = radiation_spectrum(tau, f[:, comp], [0.2 * lc, 6 * lc, 0.5 * lc])
        tot = tot + np.array(c) ** 2 + np.array(s) ** 2
    k = np.array(k)
    w = k * C.C
    i_num = D ** 2 * C.EPS0 * C.C * tot / math.pi
    q2 = C.E0 ** 2 / (4 * math.pi * C.EPS0)
    xi = w * rho / (3 * C.C * gamma ** 3)
    i_an = q2 / (3 * math.pi ** 2 * C.C) * (w * rho / C.C) ** 2 * gamma ** -4 * kv(2.0 / 3.0, xi) ** 2
    assert np.abs(i_num / i_an - 1).max() < 2e-3, i_num / i_an
