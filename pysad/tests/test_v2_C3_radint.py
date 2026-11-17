# -*- coding: utf-8 -*-
"""担当 C3 (v2next): A が「要判断」にした K64 の疑い B106・B107・B109 の独立な確認。

- B107 (RADINT の I4 の磁場勾配ベクトル): 傾いた・Y 方向の分散のある場所の K1≠0 の BEND (CKL1Y)、軌道が中心からずれた四極 (QUADRU)
  独立な基準: (a) 結合リングの物理的な I4 (A2 の 5x5 行列、`_a2_indep.coupled_i4`)、(b) 放射の 1 次の摂動論 (6x6 の固有ベクトル、`coupled_partition`)、
  (c) 平面リングの 3x3 行列の放射積分に四極の中の軌道のずれを足したもの (`_c3_indep.planar_quad_offset_jx`)、
  (d) RAD つき 6 次元追跡の 1 周写像の固有値の減衰 (RADINT・Emittance[] のどちらも使わない)。
- B109 (tbgrote の trans(5,6)): 写像の微分は、写像そのものの有限差分と一致しなければならない。
- B106 (tradke のビームの広がりの項の b): 力学的運動量は px + bzh y。`test_v2_A_spin.py` に、b = bzh の式を入れた。
条件: 3 GeV、24 セルの FODO (周長 105.6 m、RF 1.8 MV)。skew 四極 SQ (L=0.2、K1 = ks、ROTATE=45 度) 4 台。
"""
import io
import math
import os
import sys

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, radint, track6n
from pysad import lattice as _lat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _a2_indep as A       # noqa: E402
import _c3_indep as C       # noqa: E402


def _ring(ks=0.0, kb=0.0, dxf=0.0, dxd=0.0, dyf=0.0, dyd=0.0, coupled=True):
    line = '4*(CELLS 5*CELLQ)' if coupled else '24*CELLQ'
    return f'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7 DX={dxf} DY={dyf}) QD = (L=0.3 K1=-0.6 DX={dxd} DY={dyd});
QUAD  SQ = (L=0.2 K1={ks} ROTATE=0.7853981633974483);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 K1={kb} E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.8E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE CELLS = (QF D2 B D2 SQ QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF {line});
FFS;
USE RINGQ;
CELL; CALC; RFSW; RAD; NOFLUC;
Exit;
'''


def _radint_j(ffs):
    res, _m = radint.compute(ffs, ffs.ensure_rows())
    ri = [float(res['RI'][i]) for i in range(1, 12)]
    du, dv = 1 - ri[3] / ri[1], 1 - ri[6] / ri[1]
    return ri, (du, dv, 4 - du - dv)


# ---------------------------------------------------------------------------------------------
# B107: K1≠0 の偏向磁石 + Y 方向の分散 (傾いた四極の結合)
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize('kb', [-0.1, 0.1])
def test_b107_bend_k1_with_vertical_dispersion_v2(v2_mode, kb):
    """K1≠0 の偏向磁石で、結合で Y 方向の分散がある: I4U + I4V = 独立な物理的 I4 (1e-6)。DU・DV・Je が 1 次の摂動論と 3e-3 以内。
    (K64 は y 成分の勾配 CKL1Y = −cos(θ)·DKL1 が θ = 0 で 0 にならず、2 k η_y の項が余計に入る: 1〜2% ずれる)"""
    ffs, _ = ffscmd.run_sad_string(_ring(0.05, kb), out=io.StringIO())
    ri, (du, dv, de) = _radint_j(ffs)
    els = ffs.lat.elements
    i4, _i2 = A.coupled_i4(els)
    assert abs((ri[3] + ri[6]) / i4 - 1) < 1e-6, (ri[3], ri[6], i4)
    modes, _u0 = A.coupled_partition(els, ffs.momentum())
    ms = sorted(modes, key=lambda t: -t[0])                  # ν の大きい順: 横 2 つ (x 様と y 様の順は K1 で入れ替わる)、z 様 (νs)
    tr_own = sorted([ms[0][1], ms[1][1]])
    tr_got = sorted([du, dv])
    assert abs(tr_got[0] - tr_own[0]) < 3e-3 and abs(tr_got[1] - tr_own[1]) < 3e-3 and abs(de - ms[2][1]) < 3e-3, ((du, dv, de), ms)


def test_b107_bend_k1_k64_keeps_the_old_i4v(k64_mode):
    ffs, _ = ffscmd.run_sad_string(_ring(0.05, 0.1), out=io.StringIO())
    ri, _ = _radint_j(ffs)
    i4, _ = A.coupled_i4(ffs.lat.elements)
    assert abs((ri[3] + ri[6]) / i4 - 1) > 5e-3                # K64 は I4U + I4V が物理的な I4 から 1.6% ずれる (記録: v2next_C3.md)


# ---------------------------------------------------------------------------------------------
# B107b: 四極の中の軌道のずれ (QUADRU)。要素の DX・DY を引く + K1<0 の勾配の符号
# ---------------------------------------------------------------------------------------------

R_ZERO_DC = 1.0474767979479964 / 0.45593683933259485        # 経路長が変わらない (δ のずれが出ない) QF : QD のずれの比 (planar の計算)


@pytest.mark.parametrize('dxf,dxd', [(1e-4, 0.0), (0.0, 1e-4), (-1e-4, 0.0), (1e-4, 1e-4 * R_ZERO_DC), (-1e-4, -1e-4 * R_ZERO_DC)])
def test_b107b_quad_offset_jx_matches_planar_integrals_v2(v2_mode, dxf, dxd):
    """四極の据え付けのずれで軌道が四極の中心からずれる: RADINT の Jx が、平面リングの 3x3 の放射積分 (h = k x_o、I4 の 2 k² x_o η) と 3e-4 以内。
    (K64 は閉軌道だけを使い四極自身のずれを無視し、K1<0 の I4 の項の符号も逆: QD のずれで 2e-3、QF + QD で 5e-3 ずれる)"""
    ffs, _ = ffscmd.run_sad_string(_ring(0.0, 0.0, dxf, dxd, coupled=False), out=io.StringIO())
    _ri, (du, _dv, _de) = _radint_j(ffs)
    want = C.planar_quad_offset_jx(dxf, dxd)
    assert abs(du - want) < 3e-4, (du, want)


def test_b107b_quad_offset_zero_dc_matches_tracking_and_emittance_v2(v2_mode):
    """経路長が変わらない組 (δ のずれが出ない) の ΔJx が、RAD つき 6 次元追跡の 1 周写像の固有値の減衰と 3% 以内で一致する (+ΔJx = 7.8e-3)。
    Emittance[] も 1.5% 以内"""
    j0 = None
    out = {}
    for dx in (0.0, 1e-4):
        ffs, _ = ffscmd.run_sad_string(_ring(0.0, 0.0, dx, dx * R_ZERO_DC, coupled=False), out=io.StringIO())
        _ri, (du, _dv, _de) = _radint_j(ffs)
        jt = _tracking_j(ffs)
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
        je = [-2 * d * ffs.momentum() / e['energy_loss_eV'] for d in e['damping_rate'][:3]]
        out[dx] = (du, jt[0], je[0])
    d_radint = out[1e-4][0] - out[0.0][0]
    d_track = out[1e-4][1] - out[0.0][1]
    d_emit = out[1e-4][2] - out[0.0][2]
    assert abs(d_track / 7.8e-3 - 1) < 0.05
    assert abs(d_radint / d_track - 1) < 0.03, (d_radint, d_track)
    assert abs(d_emit / d_track - 1) < 0.03, (d_emit, d_track)


def _tracking_j(ffs):
    """RAD つき 6 次元追跡 (track6n、RAD の追跡 tradk を使う) の 1 周写像の固有値の減衰から分配数 (Jx, Jy, Jz)"""
    els = ffs.lat.elements
    p = ffs.momentum()
    circ = sum(x.L for x in els)
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True)
    trf0 = e['trf0']

    def f(c):
        return np.asarray(track6n.track_particles_fast(els, np.array([c]), 1, p, trf0, True, True, False, circ=circ))[0]

    def jac(c, h=1e-7):
        J = np.zeros((6, 6))
        for i in range(6):
            d = np.zeros(6)
            d[i] = h
            J[:, i] = (f(c + d) - f(c - d)) / (2 * h)
        return J
    c0 = np.array(e['orbit_entrance'], float)
    for _ in range(4):
        c0 = c0 - np.linalg.solve(jac(c0) - np.eye(6), f(c0) - c0)
    ev = np.linalg.eigvals(jac(c0))
    fold = lambda x: min(x % 1.0, 1 - x % 1.0)
    mags = []
    for t in e['tunes']:
        z = min((z for z in ev if z.imag >= 0), key=lambda z: abs(fold(math.atan2(z.imag, z.real) / 2 / math.pi) - fold(t)))
        mags.append(abs(z))
    return [2 * p * (1 - m) / e['energy_loss_eV'] for m in mags]


@pytest.mark.parametrize('dyf,dyd,want', [(1e-4, 0.0, 0.00042), (0.0, 1e-4, -0.00042), (1e-4, -2e-4, 0.00121)])
def test_b107b_vertical_offset_changes_jy_like_tracking_v2(v2_mode, dyf, dyd, want):
    """Y 方向のずれの四極 + skew 四極 (結合、η_y ≠ 0): RADINT の ΔJy (y のずれ 0 との差) が、6 次元追跡の ΔJy (+0.00042、−0.00042、+0.00121) と 1e-4 以内。
    (K64 の y 成分は符号が逆で、QF のずれで ΔJy の符号が逆になる)"""
    ffs0, _ = ffscmd.run_sad_string(_ring(0.05), out=io.StringIO())
    ffs1, _ = ffscmd.run_sad_string(_ring(0.05, 0.0, 0.0, 0.0, dyf, dyd), out=io.StringIO())
    d_radint = _radint_j(ffs1)[1][1] - _radint_j(ffs0)[1][1]
    d_track = _tracking_j(ffs1)[1] - _tracking_j(ffs0)[1]
    assert abs(d_track - want) < 1.5e-4, (d_track, want)
    assert abs(d_radint - d_track) < 1.2e-4, (d_radint, d_track)


def test_b107b_k64_ignores_quad_displacement(k64_mode):
    """K64: 四極自身のずれ DX は RADINT に入らない (閉軌道だけ)。ずれのある四極だけのリングの J が、ずれのないリングと軌道の分だけ違う"""
    ffs, _ = ffscmd.run_sad_string(_ring(0.0, 0.0, 1e-4, 1e-4 * R_ZERO_DC, coupled=False), out=io.StringIO())
    _ri, (du, _dv, _de) = _radint_j(ffs)
    assert abs(du - C.planar_quad_offset_jx(1e-4, 1e-4 * R_ZERO_DC)) > 3e-3


# ---------------------------------------------------------------------------------------------
# B109: tbgrote の trans(5,6)
# ---------------------------------------------------------------------------------------------

def _fd_jac(f, c0, h=1e-7):
    J = np.zeros((6, 6))
    for i in range(6):
        d = np.zeros(6)
        d[i] = h
        J[:, i] = (f(c0 + d) - f(c0 - d)) / (2 * h)
    return J


_TCHGE_CASES = [(1e-3, 0.0, 2e-3, 0.0, 1e-3, 0.0, 0.0),        # dx0, dy, dz, theta, dtheta, alg, phig
                (1e-3, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1),
                (1e-3, 2e-3, 3e-3, 0.2, 2e-3, 0.3, 0.15)]


@pytest.mark.parametrize('case', _TCHGE_CASES)
@pytest.mark.parametrize('ent', [True, False])
def test_b109_tchge_jacobian_equals_finite_difference_v2(v2_mode, case, ent):
    """tchge (据え付けの誤差の座標変換) の解析的な 6x6 行列 (δ 列の z 行 (5,6) を含む) = 写像の有限差分。k64 は (5,6) が px²+py² のオーダーで違う"""
    dx0, dy, dz, th, dth, alg, phig = case
    c0 = np.array([1e-3, 3e-3, -2e-3, 2e-3, 0.0, 1e-3])
    f = lambda c: np.array(_lat._tchge_cod(c.copy(), dx0, dy, dz, th, dth, alg, phig, ent))
    J = _lat._tchge_jac(c0.copy(), dx0, dy, dz, th, dth, alg, phig, ent)
    assert np.abs(J - _fd_jac(f, c0)).max() < 2e-9


@pytest.mark.parametrize('case', _TCHGE_CASES[:2])
def test_b109_k64_56_element_keeps_the_k64_value(k64_mode, case):
    dx0, dy, dz, th, dth, alg, phig = case
    c0 = np.array([1e-3, 3e-3, -2e-3, 2e-3, 0.0, 1e-3])
    f = lambda c: np.array(_lat._tchge_cod(c.copy(), dx0, dy, dz, th, dth, alg, phig, False))
    d56 = (_lat._tchge_jac(c0.copy(), dx0, dy, dz, th, dth, alg, phig, False) - _fd_jac(f, c0))[4, 5]
    assert abs(d56 + (3e-3 ** 2 + 2e-3 ** 2)) < 2e-6           # K64: (5,6) が −(px² + py²) だけ違う (dpzi (pr+1)/pzi)
