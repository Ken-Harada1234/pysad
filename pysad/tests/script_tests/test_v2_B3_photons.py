# -*- coding: utf-8 -*-
"""担当 B3 (3): PHOTONS (TrackParticles の光子の記録 PhotonList)。BEAMBEAM 以外の要素 (BEND・QUAD) の光子 (photons.f の tphrec・tphotonlist)。

根拠:
  (a) gs (K64) との照合 (k64 モード): 乱数の引き方 (tdusrnpl・radangle1 の tgauss) と tphrec の位置・方向・偏光の式まで直訳した。BEND (K1=0・K1≠0・ROTATE・DX・DY・DROTATE) と QUAD の光子
      358 個 (example/07_spin_solenoid_special/flag_PHOTONS_bend.sad、4 粒子 × 2 セル) が、gs と光子数・位置 (gx gy gz)・方向・偏光・粒子番号・要素番号まで一致する。
      エネルギー列 (E) だけ 3e-8 (相対) 違う (原因は gs の gammab と pysad の h0 の定義の違いと推定。確かめていない)。
  (b) 解析式 (gs を使わない、v2): 曲率半径 ρ の偏向磁石で、1 粒子 1 回通過あたりの光子数 = (5/(2√3)) α γ θ、臨界エネルギー u_c = (3/2) ħc γ³/ρ、
      平均エネルギー <u> = (8/(15√3)) u_c、<u²> = (11/27) u_c²、スペクトル (3/(5π)) ∫_y^∞ K_{5/3}(z) dz (y = u/u_c)。
  (c) 位置と方向: 光子の位置は偏向磁石の円弧上、方向は円弧の接線から ~1/γ の角度の中。
v2 の修正: K64 は偏光の回転で c1 = (cos·ppx1 + sin·ppy1)/pp と、すでに規格化した ppx1・ppy1 を曲げ角 pp (≈ 数 mrad) でもう一度割る。xi1 が ±1 を超え 1e3〜1e5 になる
    (gs の PhotonList の 8 列目が −3735 など)。Stokes パラメータの大きさは (xi1² + xi2² + xi3²) ≤ 1 のはずで、v2 は割らない。
条件: 3 GeV、BEND L=1.2 m・ANGLE=15 度 (ρ=4.58 m)、24 セルの FODO。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
HBARC = 1.973269804e-7          # eV m
ALPHA = 1.0 / 137.035999

_RING = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST=();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST %d*CELL1);
FFS;
USE RING; CELL; CALC;
RAD; FLUC; PHOTONS;
SeedRandom[%d];
n=%d;
x0=Table[1e-5*Sin[0.37*i],{i,n}]; px0=Table[0,{i,n}]; y0=Table[0,{i,n}];
py0=Table[0,{i,n}]; z0=Table[0.0,{i,n}]; dp0=Table[0.0,{i,n}]; fl0=Table[1,{i,n}];
r=TrackParticles[{1,{x0,px0,y0,py0,z0,dp0,fl0}},%d,1];
Exit;
'''


def _run(src):
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    return ffs, np.array(I.globals['PhotonList'], dtype=float).reshape(-1, 12)


def _ring(ncell, seed, npart):
    return _run(_RING % (ncell, seed, npart, 6 * ncell + 1))


def _parse_gs(path):
    rows_a, rows_b = [], []
    for ln in open(path, encoding='utf-8', errors='ignore'):
        if ln.startswith('PA'):
            rows_a.append([float(v) for v in ln.split()[1:]])
        elif ln.startswith('PB'):
            rows_b.append([float(v) for v in ln.split()[1:]])
    return np.hstack([np.array(rows_a), np.array(rows_b)])


def test_photon_list_matches_gs_for_bend_and_quad(k64_mode):
    """K1≠0・ROTATE・DX・DY・DROTATE つきの BEND と QUAD の光子 358 個が gs と一致する (E 列だけ 3e-8)"""
    ref = _parse_gs(os.path.join(ROOT, 'example', '07_spin_solenoid_special', 'flag_PHOTONS_bend.gs.txt'))
    src = open(os.path.join(ROOT, 'example', '07_spin_solenoid_special', 'flag_PHOTONS_bend.sad'), encoding='utf-8').read()
    _ffs, ph = _run(src)
    assert ph.shape == ref.shape == (358, 12)
    assert np.abs(ph[:, 0] / ref[:, 0] - 1).max() < 5e-8
    assert np.abs(ph[:, 1:4] - ref[:, 1:4]).max() < 1e-9          # 位置 (GX GY GZ [m])
    assert np.abs(ph[:, 4:7] - ref[:, 4:7]).max() < 2e-10         # 方向余弦 (gs の表示は 10 桁)
    assert (np.abs(ph[:, 7:10] - ref[:, 7:10]) / (np.abs(ref[:, 7:10]) + 1e-3)).max() < 5e-8        # 偏光 (K64 の xi1 は 1e6 のオーダーの値を含む)
    assert np.array_equal(ph[:, 10:], ref[:, 10:])                # 粒子番号・要素番号
    assert set(ph[:, 11].astype(int)) == {4, 10, 16, 20, 22, 28, 34}       # 20 は QUAD の光子


def test_stokes_xi1_is_unbounded_in_k64(k64_mode):
    _ffs, ph = _ring(3, 7, 4)
    assert np.abs(ph[:, 7]).max() > 1e3            # K64 の xi1 は ±1 を超える (pp で二重に割る)


def test_stokes_parameters_are_bounded_v2(v2_mode):
    _ffs, ph = _ring(3, 7, 4)
    assert np.abs(ph[:, 7:10]).max() <= 1.0 + 1e-12
    s2 = ph[:, 7] ** 2 + ph[:, 8] ** 2 + ph[:, 9] ** 2
    assert np.all(s2 <= 1.0 + 1e-9) and np.all(s2 > 0.9)         # ほぼ完全偏光 (radangle1 の xi2² + xi3² = 1)


@pytest.fixture(scope='module')
def _stat_run():
    compat.set_mode('v2')
    try:
        ffs, ph = _ring(24, 11, 100)
    finally:
        compat.set_mode('k64')
    return ffs, ph


def test_photon_number_energy_and_spectrum_v2(_stat_run):
    """光子数 N = (5/(2√3)) α γ θ × 磁石数 × 粒子数 (Poisson の ±3σ)、<u>/u_c = 8/(15√3)、<u²>/u_c² = 11/27、スペクトルの KS 検定"""
    _ffs, ph = _stat_run
    gam = math.sqrt(1 + (3e9 / 0.51099895e6) ** 2)
    theta = 0.2617993877991494
    rho = 1.2 / theta
    uc = 1.5 * HBARC * gam ** 3 / rho
    nexp = 5 / (2 * math.sqrt(3)) * ALPHA * gam * theta * 24 * 100
    assert abs(len(ph) / nexp - 1) < 3.0 / math.sqrt(nexp)
    u = ph[:, 0] / uc
    # 平均エネルギー (放射で粒子のエネルギーが下がる効果は <1e-4)。標準誤差 = std/√N
    assert abs(u.mean() - 8 / (15 * math.sqrt(3))) < 3.5 * u.std() / math.sqrt(len(u))
    assert abs((u ** 2).mean() / (11 / 27) - 1) < 0.05
    # KS: pdf(y) = (3/(5π)) ∫_y^∞ K_{5/3}(z) dz。対数の格子 (s = ln y) で積分する (y^{-2/3} の特異性を避ける)
    from scipy.integrate import cumulative_trapezoid
    from scipy.special import kv
    s_ = np.linspace(-30.0, math.log(60.0), 20001)
    y_ = np.exp(s_)
    g_ = kv(5.0 / 3.0, y_) * y_                                    # K_{5/3}(z) dz = K y ds
    f5 = cumulative_trapezoid(g_[::-1], -s_[::-1], initial=0.0)[::-1]       # ∫_y^{60} K_{5/3}
    f5 = f5 + (1.5 * math.gamma(5.0 / 3.0) * 2 ** (2.0 / 3.0)) * 0.0       # (60 より先は < 1e-25)
    pdf = 3.0 / (5.0 * math.pi) * f5
    cdf = cumulative_trapezoid(pdf * y_, s_, initial=0.0)
    cdf = cdf + 3.0 / (5.0 * math.pi) * 1.5 * math.gamma(2.0 / 3.0) * 2 ** (-1.0 / 3.0) * 0.0     # (y < e^-30 の寄与は ~1e-9)
    assert abs(cdf[-1] - 1.0) < 2e-3                              # 規格化 (<y> = 8/(15√3) は上の検算)
    ys = np.sort(u)
    emp = np.arange(1, len(ys) + 1) / len(ys)
    d = np.abs(np.interp(np.log(ys), s_, cdf / cdf[-1]) - emp).max()
    assert d < 1.63 / math.sqrt(len(ys))                           # KS の 1% 点 (1.63/√N)


def test_photon_position_and_direction_v2(_stat_run):
    """光子の位置は円弧 (幾何の部分要素) の上 (粒子の軌道のずれ 2 mm 以内)、方向は円弧の接線から ~1/γ (u > 0.5 u_c の光子で rms 0.3〜3 /γ)"""
    ffs, ph = _stat_run
    g = ffs.ensure_geo()
    gam = math.sqrt(1 + (3e9 / 0.51099895e6) ** 2)
    from scipy.optimize import minimize_scalar
    uc = 1.5 * HBARC * gam ** 3 / (1.2 / 0.2617993877991494)
    sel = ph[ph[:, 0] > 0.5 * uc][::11]               # 高いエネルギーの光子 (角度 ~1/γ。低いエネルギーは (u_c/u)^{1/3}/γ に広がる)
    frs = np.linspace(0.0, 1.0, 41)
    dmax = 0.0
    ang = []
    for r in sel:
        l = int(r[11])
        dist = lambda f: float(np.linalg.norm(g.geo_frac(l, f)[:, 3] - r[1:4]))
        f0 = frs[int(np.argmin([dist(f) for f in frs]))]
        res = minimize_scalar(dist, bounds=(max(0.0, f0 - 0.025), min(1.0, f0 + 0.025)), method='bounded', options={'xatol': 1e-7})
        dmax = max(dmax, res.fun)
        tang = g.geo_frac(l, res.x)[:, 2]
        ang.append(math.acos(min(1.0, float(r[4:7] @ tang))))
    assert dmax < 2e-3          # 粒子の軌道のずれ (分散 × 放射で減った δ ≈ 1e-3、RF なし) の範囲
    rms = math.sqrt(np.mean(np.array(ang) ** 2))
    assert 0.3 / gam < rms < 3.0 / gam


def test_photons_with_beambeam_and_bend_are_recorded(v2_mode):
    """D1: BEAMBEAM と BEND を同じラティスに置いた PHOTONS (以前は SADError で未対応)。BEND と BEAMBEAM の両方の光子が PhotonList に入る
    (gs との一致は test_v2_D1_emit.py の k64 の試験)"""
    src = _RING % (2, 1, 2, 13)
    src = src.replace('MARK ST=();', 'MARK ST=();\nBEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11 BSTRL=1);').replace('LINE CELL1 = (QF D2 B D2 QD D1);', 'LINE CELL1 = (QF D2 B BB D2 QD D1);')
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    assert 'BEAMBEAM と BEND' not in buf.getvalue()
    pl = I.globals['PhotonList']
    assert len(pl) > 0 and len({int(r[11]) for r in pl}) >= 2
