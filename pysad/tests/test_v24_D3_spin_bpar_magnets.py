# -*- coding: utf-8 -*-
"""v2.4.0 担当 D3: 磁石の本体の運動量の向きの磁場 (1+a)B∥ ds の範囲拡大 (docs/21 B230・B231。B165 の続き)。

独立な基準 `_spin_maxwell_indep2.maxwell_R`: Maxwell 方程式を満たす直線の磁石 (QUAD・MULT・SEXT・ステアリング) の中の Thomas-BMT の数値積分
(フリンジ幅 w = 8e-3 と 4e-3 m の結果を w² で外挿)。入口の面の座標 (x, px, y, py) = (2e-3, 3e-3, 2e-3, 3e-3)、0.5 GeV、L = 0.4 m。
z 軸のまわりの回転 R01・R10 が運動量の向きの磁場の積分 (1+a)∫B∥ds に対応する。
  - QUAD・MULT: K64 (と B165 の前の v2) は端の ψ だけで本体 Δψ が無く、R01 が 6.5e-6〜1.0e-5 ずれる → v2 で 1e-8 以下。
  - SEXT・ステアリング (K64 は端の ψ も本体の B∥ も入れない = 3 つの和が 0): 元から 1e-7 以内で一致。
  - K64 の追跡の QUAD は刻みが 1 つ (ndiv = 1) のとき入口の +ψ を出口の −ψ が上書きして落とす (v2 は足す)。
  - 追跡の BEND の刻み (1 刻みの偏向角) を 5e-3 rad 以下にして dR/dpy の R12・R21 の誤差 (刻みに比例する 1 次の誤差) を 1e-2 → 6e-4。"""
import io

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, track6

from . import _spin_maxwell_indep2 as M

E, L = 0.5, 0.4
C = (2e-3, 3e-3, 2e-3, 3e-3)


def _els(decl):
    src = f'MOMENTUM = {E} GEV;\n{decl};\nMARK ST = ();\nLINE RING = (ST M1);\nFFS;\nUSE RING;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return [e for e in ffs.lat.elements if e.name == 'M1'], ffs.momentum()


def _ref(mult, **kw):
    r8, _ = M.maxwell_R(E * 1e9, L, mult, C, w=8e-3, **kw)
    r4, _ = M.maxwell_R(E * 1e9, L, mult, C, w=4e-3, **kw)
    return (4 * r4 - r8) / 3


def _track(els, p):
    cols = []
    for j in range(3):
        _o, s = track6.track_particles(els, np.array([[C[0], C[1], C[2], C[3], 0, 0.0]]), 1, p, 0.0, False, True, False,
                                       spins=[np.eye(3)[j]], radpol=False, circ=100.0)
        cols.append(s[0])
    return np.array(cols).T


def _emit(els, p):
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne(els, np.array([C[0], C[1], C[2], C[3], 0, 0.0]), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
        return st.srot[:, :3].copy()
    finally:
        emit6d.POL_MODE = False


CASES = {
    'quad1': ('QUAD M1=(L=0.4 K1=0.4)', {1: 1.0}),
    'quad8': ('QUAD M1=(L=0.4 K1=3.2)', {1: 8.0}),
    'mult': ('MULT M1=(L=0.4 K1=0.4 K2=16)', {1: 1.0, 2: 40.0}),
    'mult_k0': ('MULT M1=(L=0.4 K0=0.002 K1=0.2 K2=12)', {0: 0.005, 1: 0.5, 2: 30.0}),
    # v2.4.0 D3 (続き): ずれ・ROTATE・K1 < 0 (−90° 回転)・F1/FRINGE つき (F1 の断片の放射は ψ を足さない: 本体の Δψ だけで和が 0)
    'quad_dxdy': ('QUAD M1=(L=0.4 K1=0.4 DX=0.5E-3 DY=-0.3E-3)', {1: 1.0}, dict(dx=0.5e-3, dy=-0.3e-3)),
    'quad_rot': ('QUAD M1=(L=0.4 K1=0.4 ROTATE=0.3)', {1: 1.0}, dict(rot=0.3)),
    'quad_neg': ('QUAD M1=(L=0.4 K1=-0.4)', {1: -1.0}),
    'quad_f1': ('QUAD M1=(L=0.4 K1=0.4 F1=0.03 FRINGE=3)', {1: 1.0}),
    'mult_rot': ('MULT M1=(L=0.4 K1=0.4 K2=16 DX=0.5E-3 ROTATE=0.2)', {1: 1.0, 2: 40.0}, dict(dx=0.5e-3, rot=0.2)),
    'mult_f1': ('MULT M1=(L=0.4 K1=0.4 K2=16 F1=0.03 FRINGE=3)', {1: 1.0, 2: 40.0}),
}


@pytest.fixture(scope='module')
def refs():
    return {k: _ref(v[1], **(v[2] if len(v) > 2 else {})) for k, v in CASES.items()}


@pytest.mark.parametrize('name', list(CASES))
def test_body_bpar_tracking_and_emittance(name, refs):
    els, p = _els(CASES[name][0])
    R0 = refs[name]
    try:
        compat.set_mode('v2')
        compat.set_spin_bpar(False)
        d_old = _track(els, p) - R0
        compat.set_spin_bpar(True)
        d_trk = _track(els, p) - R0
        d_em = _emit(els, p) - R0
        compat.set_mode('k64')
        d_k64 = _track(els, p) - R0
    finally:
        compat.set_spin_bpar(None)
        compat.set_mode(None)
    z_old = max(abs(d_old[0, 1]), abs(d_old[1, 0]))
    assert z_old > 4e-6                                    # 本体の B∥ が無いと z 軸まわりが 1e-5 ずれる
    assert abs(d_k64[0, 1]) > 4e-6                         # k64 は変えない
    tol = 3e-8 if name in ('quad1', 'mult') else 1.5e-7
    big = 1e-5 if name != 'mult_rot' else 1e-6
    assert abs(d_trk[0, 1]) < tol and abs(d_trk[1, 0]) < tol
    assert abs(d_em[0, 1]) < tol and abs(d_em[1, 0]) < tol
    # 残る誤差 (R02・R12 など、B⊥ に対応) は軌道の違い (K64 の端の模型と Maxwell のフリンジ) 由来で、B∥ とは別
    assert np.abs(d_trk).max() < 1e-5
    if name == 'mult_rot':
        assert np.abs(d_trk).max() < 1e-6          # MULT の ROTATE: 追跡もスピンを回す (B241。以前は全体が 5.7e-4 ずれた)


@pytest.mark.parametrize('kind', ['sext', 'steer'])
def test_sext_steering_already_consistent(kind):
    """SEXT (L ≠ 0) とステアリングは K64 が端の ψ も本体の B∥ も入れない (3 つの和 0)。追跡は元から基準と合う"""
    if kind == 'sext':
        els, p = _els('SEXT M1=(L=0.4 K2=24)')
        mult = {2: 60.0}
    else:
        els, p = _els('BEND M1=(L=0.4 ANGLE=0 K0=0.003)')
        mult = {0: 0.003 / L}
    R0 = _ref(mult)
    try:
        compat.set_mode('v2')
        d = _track(els, p) - R0
    finally:
        compat.set_mode(None)
    assert np.abs(d).max() < 1e-7


def test_quad_single_slice_entrance_psi_not_overwritten():
    """ndiv = 1 (弱い四極) の追跡: K64 は入口の +ψ を出口の −ψ で上書きする。v2 は足す (B165 を切っても入口の ψ は残る)"""
    els, p = _els('QUAD M1=(L=0.4 K1=0.4)')
    R0 = _ref({1: 1.0})
    try:
        compat.set_mode('v2')
        compat.set_spin_bpar(False)
        d_v2 = _track(els, p) - R0
        compat.set_mode('k64')
        d_k = _track(els, p) - R0
    finally:
        compat.set_spin_bpar(None)
        compat.set_mode(None)
    # 入口 ψ 4e-6 の分だけ v2 (本体なし) の誤差が k64 より小さい (k64 −1.0e-5、v2 本体なし −6.2e-6 = −Δψ)
    assert abs(d_k[0, 1]) - abs(d_v2[0, 1]) > 3e-6


def test_bend_default_slicing_converges():
    """追跡の BEND (15°、L = 1.2 m、0.5 GeV) の dR/dpy の R12・R21: 刻みの偏向角を 5e-3 rad 以下にして 1e-2 → 1e-3 以下"""
    from . import test_v24_D_spin_bpar as T
    els, p = T._setup()
    try:
        compat.set_mode('v2')
        compat.set_spin_bpar(True)
        d = T._lin(lambda py: T._R_track(els, p, py)) - T.REF
    finally:
        compat.set_spin_bpar(None)
        compat.set_mode(None)
    assert abs(d[1, 2]) < 1.5e-3 and abs(d[2, 1]) < 1.5e-3 and abs(d[1, 0]) < 1e-3 and abs(d[0, 1]) < 2e-3


def test_sext_emittance_spin_and_radiation_match_tracking():
    """B240: Emittance[] の 6 極 (L ≠ 0) は K64 の tthine が放射しないので spin (srot) が入らなかった (追跡は放射・歳差する)。v2 は追跡と同じ 1 刻みの tradke を入れる。
    スピン回転行列は Maxwell 場の積分と 5e-6 以内 (残りは Emittance[] が 6 極のフリンジ ttfrin を入れないための軌道の差)、k64 は 3.5e-4 ずれたまま。
    座標 (3, 2) mm を通る 1 刻みの δ の損失 U0 (st.u0) は追跡の損失と 1e-5 以内"""
    els, p = _els('SEXT M1=(L=0.4 K2=24)')
    R0 = _ref({2: 60.0})
    try:
        compat.set_mode('v2')
        d2 = _emit(els, p) - R0
        compat.set_mode('k64')
        d1 = _emit(els, p) - R0
        compat.set_mode('v2')
        els2, p2 = _els('SEXT M1=(L=0.4 K2=9.6)')
        c = np.array([3e-3, 0, 2e-3, 0, 0, 0.0])
        o = track6.track_particles(els2, c[None, :], 1, p2, 0.0, False, True, False, circ=100.0)
        st = emit6d.tturne(els2, c, emit6d._Phys(p2, 100.0), rfsw=False, radcod=False)
    finally:
        compat.set_mode(None)
    assert np.abs(d2).max() < 5e-6
    assert np.abs(d1).max() > 3e-4
    assert st.u0 == pytest.approx(-o[0, 5], rel=1e-5)


@pytest.mark.parametrize('dp', [0.02, -0.03])
def test_body_bpar_with_momentum_offset(dp):
    """B242: 端の縦磁場の積分 ψ (p0 で規格化) を粒子の運動量 (1+δ) で規格化し直す。K64 は割らないので、本体の Δψ/(1+δ) との和が δΔψ だけ 0 にならない
    (δ = 0.02 の QUAD で R01 に 1.2e-7 = δ × 6e-6)。v2: 追跡・Emittance[] とも Maxwell 場の積分 (運動量 (1+δ) の粒子) と 3e-8 以内。k=1 の QUAD、座標 (2, 3, 2, 3) mm"""
    els, p = _els('QUAD M1=(L=0.4 K1=0.4)')
    r8, _ = M.maxwell_R(E * 1e9, L, {1: 1.0}, C, dp=dp, w=8e-3)
    r4, _ = M.maxwell_R(E * 1e9, L, {1: 1.0}, C, dp=dp, w=4e-3)
    R0 = (4 * r4 - r8) / 3

    def trk():
        cols = []
        for j in range(3):
            _o, s_ = track6.track_particles(els, np.array([[C[0], C[1], C[2], C[3], 0, dp]]), 1, p, 0.0, False, True, False,
                                            spins=[np.eye(3)[j]], radpol=False, circ=100.0)
            cols.append(s_[0])
        return np.array(cols).T

    def emi():
        emit6d.POL_MODE = True
        try:
            st = emit6d.tturne(els, np.array([C[0], C[1], C[2], C[3], 0, dp]), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
            return st.srot[:, :3].copy()
        finally:
            emit6d.POL_MODE = False
    try:
        compat.set_mode('v2')
        dt, de = trk() - R0, emi() - R0
    finally:
        compat.set_mode(None)
    assert abs(dt[0, 1]) < 3e-8 and abs(dt[1, 0]) < 3e-8
    assert abs(de[0, 1]) < 3e-8 and abs(de[1, 0]) < 3e-8
