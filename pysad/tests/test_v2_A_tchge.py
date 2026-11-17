# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): tchge (据え付け誤差の座標変換) — B09 (入口の回転が一周行列に入らない) と B56 (tbrote の srot の逐次更新)。

根拠:
- B09: 転送行列は軌道の写像の微分でなければならない。v2 の `_tchge_jac` は `_tchge_cod` の数値微分と一致する。
  さらに、ロールとずれのある BEND・六極を含むリングで、CALC のチューンが、別のコードである追跡 (TrackParticles 相当) の
  一周写像の固有値から求めたチューンと一致する (K64 は 0.02〜0.18 ずれる)。
- B56: tbrote の rr は回転行列から単位行列を引いたもの。スピンの 3x3 行列への作用は行列積 (I+rr)·srot で、直交行列は直交行列に移る。
  K64 は行を逐次更新するので、直交性が崩れる。

条件: 3 GeV、24 セルの FODO リング (周長 105.6 m)、矩形の偏向磁石。
"""
import io
import math

import numpy as np

from pysad import ffscmd, track6n, lattice as L, spin


def _num_jac(f, c, h=1e-7):
    J = np.zeros((6, 6))
    for i in range(6):
        d = np.zeros(6)
        d[i] = h
        J[:, i] = (f(c + d) - f(c - d)) / (2 * h)
    return J


ARGS = [  # dx, dy, dz, theta, dtheta, alg, phig
    (2e-3, -1e-3, 0.0, 0.3, 0.0, 0.0, 0.0),        # 薄い多極・空洞の形 (ロールとずれ)
    (2e-3, 1e-3, 0.0, 0.3, 0.0, 0.6, 0.13),        # BEND (ロール ROTATE とずれ、中点まわり)
    (1e-3, 2e-3, 5e-4, -0.2, 1e-3, 0.6, 0.13),     # BEND (DZ・DROTATE も)
]
COD = np.array([1.2e-3, 3e-4, -8e-4, 2e-4, 1e-3, 2e-3])


def test_tchge_entrance_matrix_is_jacobian_v2(v2_mode):
    for a in ARGS:
        J = L._tchge_jac(COD, *a, True)
        Jn = _num_jac(lambda c: L._tchge_cod(c, *a, True), COD)
        assert np.abs(J[:4, :4] - Jn[:4, :4]).max() < 1e-8, a
        assert np.abs(J[:4, 5] - Jn[:4, 5]).max() < 1e-8


def test_tchge_entrance_matrix_k64_has_no_rotation(k64_mode):
    a = ARGS[0]
    J = L._tchge_jac(COD, *a, True)
    Jn = _num_jac(lambda c: L._tchge_cod(c, *a, True), COD)
    assert np.abs(J[:4, :4] - np.eye(4)).max() < 1e-12          # K64: 回転が落ちて単位行列
    assert np.abs(J[:4, :4] - Jn[:4, :4]).max() > 0.2


def test_tchge_exit_matrix_is_jacobian_both_modes():
    """出口は K64 でも回転が入る (v2 で変えていない)"""
    from pysad import compat
    for m in ('k64', 'v2'):
        old = compat._mode
        compat.set_mode(m)
        try:
            for a in ARGS:
                J = L._tchge_jac(COD, *a, False)
                Jn = _num_jac(lambda c: L._tchge_cod(c, *a, False), COD)
                assert np.abs(J[:4, :4] - Jn[:4, :4]).max() < 1e-8
        finally:
            compat.set_mode(old)


def _ring(kind):
    extra = {
        'bend': 'BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5) '
                'BR = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5 ROTATE=0.3 DY=1E-3);\nSEXT SX=(L=0 K2=0);',
        'sext': 'BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5) '
                'BR = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);\nSEXT SX=(L=0 K2=8 ROTATE=0.4 DX=3E-3 DY=-2E-3);'}[kind]
    return '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
%s
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4E6 FREQ=5.0E8 PHI=0);
LINE CELL = (QF D2 B D2 QD DS D1 DS);
LINE CELLR = (QF D2 BR D2 QD DS SX D1 DS);
LINE RING = (ORIGIN RF 23*CELL CELLR);
FFS;
USE RING; CELL; CALC;
Print[{Twiss["NX","$$$"],Twiss["NY","$$$"]}/2/Pi];
''' % extra


def _tunes(kind):
    buf = io.StringIO()
    ffs, _I = ffscmd.run_sad_string(_ring(kind), out=buf)
    o = [l for l in buf.getvalue().splitlines() if l.startswith('{')]
    nu = [float(x) for x in o[0].strip('{}').split(',')]
    calc = sorted(min(x % 1, 1 - x % 1) for x in nu)
    els = ffs.lat.elements
    p = ffs.momentum()

    def f(v):
        return np.asarray(track6n.track_particles_fast(els, np.array([list(v) + [0.0, 0.0]]), 1, p, 0.0,
                                                       False, False, False))[0][:4]

    def jac(v, h=1e-6):
        J = np.zeros((4, 4))
        for i in range(4):
            d = np.zeros(4)
            d[i] = h
            J[:, i] = (f(v + d) - f(v - d)) / (2 * h)
        return J
    v = np.zeros(4)
    for _ in range(6):
        v = v - np.linalg.solve(jac(v) - np.eye(4), f(v) - v)
    ev = np.linalg.eigvals(jac(v))
    trk = sorted(set(round(abs(math.atan2(e.imag, e.real)) / 2 / math.pi, 9) for e in ev))
    return calc, trk


def test_calc_tunes_match_tracking_v2(v2_mode):
    """ロール + ずれのある六極・BEND を含むリング: CALC のチューン = 追跡の一周写像のチューン"""
    # BEND は 2e-5: CALC (tbende) と追跡 (tbend) は、ずれた BEND の 2 次の幾何の項が違う (ロールなしの DX=2 mm でも 2.9e-5 の差が K64 にある)
    for kind, tol in (('sext', 1e-7), ('bend', 2e-5)):
        calc, trk = _tunes(kind)
        assert len(trk) == 2
        assert max(abs(a - b) for a, b in zip(calc, trk)) < tol, (kind, calc, trk)


def test_calc_tunes_k64_differ_from_tracking(k64_mode):
    for kind in ('sext', 'bend'):
        calc, trk = _tunes(kind)
        assert max(abs(a - b) for a, b in zip(calc, trk)) > 5e-3, (kind, calc, trk)


# ------------------------------------------------------------------ B56: tbrote の srot
class _St:
    def __init__(self):
        self.srot = spin.srotinit()
        rng = np.random.default_rng(3)
        # 直交行列 (適当な回転) と微分の列
        q, _r = np.linalg.qr(rng.normal(size=(3, 3)))
        self.srot[:, :3] = q
        self.srot[:, 3:] = rng.normal(size=(3, 6))


def test_tbrote_rr_is_rotation_minus_identity():
    for ent in (True, False):
        rr = spin.tbrote_rr(0.6, 0.13, 0.02, 0.01, ent)
        R = np.eye(3) + rr
        assert np.abs(R @ R.T - np.eye(3)).max() < 1e-14
        assert abs(np.linalg.det(R) - 1.0) < 1e-14


def test_tbrote_spin_is_matrix_product_v2(v2_mode):
    st = _St()
    s0 = st.srot.copy()
    spin.tbrote_spin(st, 0.6, 0.13, 0.02, 0.01, True)
    R = np.eye(3) + spin.tbrote_rr(0.6, 0.13, 0.02, 0.01, True)
    assert np.abs(st.srot - R @ s0).max() < 1e-15
    q = st.srot[:, :3]
    assert np.abs(q @ q.T - np.eye(3)).max() < 1e-14            # 直交性が保たれる
    # 入口と出口の回転は互いに逆
    spin.tbrote_spin(st, 0.6, 0.13, 0.02, 0.01, False)
    assert np.abs(st.srot - s0).max() < 1e-14


def test_tbrote_spin_k64_breaks_orthogonality(k64_mode):
    st = _St()
    spin.tbrote_spin(st, 0.6, 0.13, 0.02, 0.01, True)
    q = st.srot[:, :3]
    assert np.abs(q @ q.T - np.eye(3)).max() > 1e-5
