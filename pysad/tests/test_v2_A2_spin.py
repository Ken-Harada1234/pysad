# -*- coding: utf-8 -*-
"""担当 A2 (1): example/07 の flag_POL_RADPOL で、500 周後のスピン角が gs・k64 (0.8794 rad) と v2 (1.1516 rad) で違う原因の切り分け。

結論 (docs/history/v2next_A2.md §1): 原因は A の B56・B58・B64・K1 つき BEND の修正ではなく、C2 の N8 (sprot の刻みの角)。
K64 の sprot は歳差の角が (1+aγ)·2 sin(θ/2) になり、放射の分割 (EPS) の 1 刻みの偏向角 θ について θ²/24 だけ小さい。
v2 は厳密な θ。根拠 (gs を使わない): 要素ごとに「スピンの回転 = a γ × 運動量の向きの変化 + (向きの変化 − 基準軌道の回転)」
(Thomas-BMT: スピンと運動量のなす角は a γ × 偏向角、ν_s = aγ) の和が、追跡の 1 周のスピン角と合う。
条件: 3 GeV、24 セルの平面 FODO リング (周長 105.6 m、偏向磁石 L=1.2 m 15 度)、RAD NOFLUC POL、RF 4 MV、閉軌道の近く (z = zeq)。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat

GS = 0.001159652181280002
ME = 0.51099895e6
CL = 299792458.0


def _ring(e_gev=3.0, k1=0.0, eps=None, ncell=24, sol=None):
    ep = '' if eps is None else ' EPS=%g' % eps
    cells = 'LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);\nLINE RING = (ORIGIN RF %d*CELL);\n' % ncell
    return f'''MOMENTUM = {e_gev}GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE={2 * math.pi / ncell} E1={math.pi / ncell / 2} E2={math.pi / ncell / 2} K1={k1}{ep});
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
{cells}FFS;
USE RING;
CELL; CALC;
Exit;
'''


def _make(src):
    from pysad.ffscmd import run_sad_string
    ffs, I = run_sad_string(src, out=io.StringIO())
    ffs.flags['RAD'] = True
    ffs.flags['FLUC'] = False
    ffs.flags['POL'] = True
    ffs.flags['RFSW'] = True
    return ffs


def _per_element(ffs, e_gev, z0):
    """1 周を要素ごとに追い、(1 周のスピン角、要素ごとの残差の和)。残差 = 追跡の Δ(スピン角) − [aγ·θ_lab + (θ_lab − θ0)]、
    θ_lab = 基準軌道の回転 θ0 − Δx' (SAD の ANGLE>0 は x' を減らす向き)、γ は要素の前後の δ の平均のエネルギー"""
    gam0 = math.sqrt(1 + (e_gev * 1e9 / ME) ** 2)
    n = len(ffs.els) + 1
    prev = [[0.0], [0.0], [0.0], [0.0], [z0], [0.0], [0.0], [0.0], [1.0]]
    full = ffs.fn_track([[1.0, [list(c) for c in prev]], 1.0])[1]
    res = 0.0
    for i in range(1, n):
        out = ffs.fn_track([[float(i), [list(c) for c in prev]], float(i + 1)])[1]
        cur = [[out[k][0]] for k in range(9)]
        d = (cur[6][0] - prev[6][0] + math.pi) % (2 * math.pi) - math.pi
        th0 = getattr(ffs.els[i - 1], 'ANGLE', 0.0) or 0.0
        dp0, dp1 = prev[5][0], cur[5][0]
        dxp = math.asin(cur[1][0] / (1 + dp1)) - math.asin(prev[1][0] / (1 + dp0))
        thl = th0 - dxp
        g = gam0 * (1 + 0.5 * (dp0 + dp1))
        res += d - (GS * g * thl + (thl - th0))
        prev = cur
    return full[6][0], res


ZEQ = 0.038315657


def test_k64_vs_v2_turn_angle_is_slice_angle_error(v2_mode):
    """example/07 のリング: 1 周のスピン角の v2 − k64 は 5.6e-4 rad (500 周で 0.28 rad)。放射の分割 EPS を細かくすると両者が 1e-6 で一致する
    (分割 θ²/24 の誤差)。原因は C2 の N8 で、A の B56・B58・B64 ではない (y=0 の閉軌道の近くでは端の縦磁場の項は 0)"""
    def turn(mode, eps):
        compat.set_mode(mode)
        ffs = _make(_ring(eps=eps))
        return ffs.fn_track([[1.0, [[0.0], [0.0], [0.0], [0.0], [ZEQ], [0.0], [0.0], [0.0], [1.0]]], 1.0])[1][6][0]
    d0 = turn('v2', None) - turn('k64', None)
    assert 4e-4 < d0 < 7e-4, d0
    d1 = turn('v2', 1e-5) - turn('k64', 1e-5)
    assert abs(d1) < 2e-6, d1
    assert abs(d1) < abs(d0) / 100


def test_per_element_spin_is_a_gamma_times_deflection_v2(v2_mode):
    """v2: 1 周の要素ごとの Thomas-BMT の和 (aγ × 運動量の向きの変化) が、追跡のスピン角と 2e-5 rad 以内"""
    full, res = _per_element(_make(_ring()), 3.0, ZEQ)
    assert abs(res) < 2e-5, res
    # 閉軌道の近く: 1 周の角は 2π aγ (mod 2π) にエネルギーの損失・加速 (δ の平均) の分を足したもの
    gam = math.sqrt(1 + (3e9 / ME) ** 2)
    base = (2 * math.pi * GS * gam + math.pi) % (2 * math.pi) - math.pi
    assert 0.005 < full - base < 0.02


def test_per_element_spin_k64_has_the_slice_error(k64_mode):
    """k64: 同じ和が −5.5e-4 rad ずれる (N8)。k64 の動作の記録"""
    full, res = _per_element(_make(_ring()), 3.0, ZEQ)
    assert 3e-4 < abs(res) < 8e-4, res


@pytest.mark.parametrize('k1', [0.05, -0.05])
def test_per_element_spin_with_k1_bend_v2(v2_mode, k1):
    """K1 つき BEND (tbendi の経路): 同じ和が 5e-5 rad 以内 (端の縦磁場の項は y=0 では 0)"""
    full, res = _per_element(_make(_ring(k1=k1)), 3.0, ZEQ)
    assert abs(res) < 5e-5, res


_RING_SK = '''MOMENTUM = 0.5 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 K1=%s E1=0.1308996938995747 E2=0.1308996938995747);
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


@pytest.mark.parametrize('k1,bzT', [(0.0, 1.0), (0.05, 0.0), (0.05, 1.0), (-0.05, 2.5)])
def test_ring_spin_tune_with_k1_bends_and_solenoid_v2(v2_mode, k1, bzT):
    """K1 つき BEND (tbendi の経路) と 1 m のソレノイド (BZ=bzT) のある平面リング (0.5 GeV、aγ=1.135) の 1 周のスピン回転:
    cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1+a) bz L (閉軌道の上のスピンの 3x3 行列の跡から)。解析式と 5e-5 以内
    (K1 の大きさは勾配 0.05/m²。軌道が閉軌道から外れる分は無い: 設計軌道)"""
    from .test_v2_C2_spin import _run as _runc2, _track as _trackc2
    ffs = _runc2(_RING_SK % (k1, bzT))
    n = len(ffs.els) + 1
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _b, s = _trackc2(ffs, np.zeros((1, 6)), 1, n, spin=[s0])
        M[:, j] = s[0]
    assert np.max(np.abs(M @ M.T - np.eye(3))) < 1e-8
    cmu = (np.trace(M) - 1.0) / 2.0
    gam = math.sqrt(1 + (0.5e9 / ME) ** 2)
    phi = (1 + GS) * (bzT * CL / 0.5e9) * 1.0
    half = math.cos(math.pi * GS * gam) * math.cos(phi / 2)
    assert abs(cmu - (2 * half * half - 1)) < 5e-5, (cmu, 2 * half * half - 1)
