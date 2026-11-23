# -*- coding: utf-8 -*-
"""v2.4.0 (放射・エミッタンス) B115 の残り: ANGLE=0・K0≠0 の極 (挿入光源) の RADINT の εx と、Emittance[] との差の原因。

独立な計算は `_rad_indep.py` (直線の基準軌道の中の磁場として I1..I5 を積分。4 次元の放射つき写像 + 空洞を RK4 で作り固有値から Jx)。
結論: (1) RADINT は CALC の分散を seed にするが、K64 の CALC は極の分散の発生 (η' += K0 δ) を入れない。v2 は極の分散を 1 周の周期解で補正し、
極の中は直線の基準軌道の場として積分する (扇形近似の h² の収束・I4 の h³ の項を入れない)。独立な積分と I1..I5 が 5e-5 以内。
(2) Emittance[] は正しい: RF なしでは RADINT と一致し (εx 1e-4)、RF ありの Jx の +0.7% (εx −1.4%) は、独立な 4 次元の写像 (空洞つき) が再現する。
条件: ASTRID2 + MPW (0.58 GeV、周長 45.7 m) と、16 セルの小リング (1 GeV、極の K0 = 0.06・−0.12・0.06 rad)。"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, radint

from . import _rad_indep as R

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
MPW = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
CQ = 55 / (32 * math.sqrt(3)) * 3.8615926796e-13


def _ffs(src, ring):
    ffs, _ = ffscmd.run_sad_string(src + f'\nFFS;\nUSE {ring};\nCELL; CALC; RFSW; RAD; NOFLUC;\nExit;\n', out=io.StringIO())
    return ffs


def _radint_I(ffs):
    res, _m = radint.compute(ffs, ffs.ensure_rows())
    return np.array([float(res['RI'][i]) for i in range(1, 6)])


def _emit(ffs, rf, rc, div=20.0):
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [div] * 3
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=rf, radcod=rc)
    finally:
        emit6d.EMITDIV[:] = old
    dr = np.abs(e['damping_rate'])
    return e, 2 * ffs.momentum() * dr[0] / e['energy_loss_eV']


def _ex_from_I(I, p):
    g = p / 0.51099895e6
    jx = 1 - I[3] / I[1]
    return jx, CQ * g * g * I[4] / (jx * I[1])


def _mpw_ffs():
    if not os.path.exists(MPW):
        pytest.skip('ASTRID2 + MPW のラティスがない')
    with open(MPW, encoding='utf-8') as fh:
        return _ffs(fh.read(), 'RINGMPW')


def test_radint_poles_equal_independent_integrals_v2(v2_mode):
    """RADINT の I1..I5 (極の分散を補正) が、独立な積分 (直線の基準軌道の極) と 5e-5 以内。補正なしは I1 が 4e-3、I5 が 2.1e-2 ずれた"""
    ffs = _mpw_ffs()
    I, _C = R.integrals_pole_rich(ffs.lat.elements)
    ri = _radint_I(ffs)
    for a, b in zip(I, ri):
        assert abs(b / a - 1) < 5e-5, (list(I), list(ri))


def test_radint_ex_equals_emittance_without_rf_v2(v2_mode):
    """RADINT の Jx・εx が Emittance[] (RFSW OFF、EMITDIV=20) と Jx 2e-4・εx 1e-4 以内。独立な積分とも一致"""
    ffs = _mpw_ffs()
    p = ffs.momentum()
    jx_r, ex_r = _ex_from_I(_radint_I(ffs), p)
    e, jx_e = _emit(ffs, False, False)
    assert abs(jx_r - jx_e) < 2e-4
    assert abs(ex_r / e['emittance_x'] - 1) < 1e-4
    Ii, _ = R.integrals_pole_rich(ffs.lat.elements)
    jx_i, ex_i = _ex_from_I(Ii, p)
    assert abs(ex_r / ex_i - 1) < 5e-5 and abs(jx_r - jx_i) < 1e-4


def test_radint_k64_ignores_poles(k64_mode):
    """k64 は極を積分に入れない (K64 のまま)。MPW なしの ASTRID2 と I2 が同じ"""
    ffs = _mpw_ffs()
    assert abs(_radint_I(ffs)[1] - 3.870436) < 1e-4


def test_rf_on_jx_shift_is_physical_4d_map(v2_mode):
    """16 セルの小リング + 3 極の挿入光源: Emittance[] の Jx は RF なし 0.6890 → RF あり 0.7040 (RF の位置の分散)。
    独立な 4 次元 (x, px, δ, z) の放射つき写像 + 空洞 (RK4) の固有値から Jx = 2 E (1−|λx|)/U0 = 0.7038 (差 3e-4)。チューン (x 0.38954、z 0.03404) も 1e-4 で一致"""
    src, elems = R.toy_ring(16, 0.06)
    ffs = _ffs(src, 'RING')
    _e0, jx_off = _emit(ffs, False, False)
    e1, jx_on = _emit(ffs, True, True)
    r4 = R.Ring4(1e9, elems, V=2e6, harm=40)
    ev = r4.solve()
    ph = np.angle(ev) / 2 / math.pi
    lam = [e for e, f in zip(ev, ph) if f > 0.2][0]
    jx_map = 2 * (1 - abs(lam)) / r4.u0
    assert abs(jx_on - jx_map) < 1e-3, (jx_on, jx_map)
    assert abs(jx_on - jx_off) > 0.01
    tx = [abs(f) for f in ph if abs(f) > 0.2][0]
    tz = [abs(f) for f in ph if 0 < abs(f) < 0.2][0]
    assert abs(tx - abs(e1['tunes'][0])) < 1e-4
    assert abs(tz - abs(e1['tunes'][2])) < 1e-4
