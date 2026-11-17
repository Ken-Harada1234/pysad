# -*- coding: utf-8 -*-
"""v2 の修正 B34: BEAMBEAM の傾き回転 (K64 は cost=sin(傾き角)、sint=傾き角)。

根拠: 強ビームの包絡を角 t だけ回したもの (Σ' = R Σ R^T) の蹴りは、粒子を -t 回し、回していない包絡で蹴り、+t 回したものに等しい
(回転に対する共変性)。K64 の式は cost²+sint² ≠ 1 で回転にならず、この等式を満たさない。
"""
import math

import numpy as np

from pysad import beambeam, compat


def _params(t):
    sx, sy, sp, sq = 2.0e-4, 2.0e-5, 1.0e-4, 5.0e-5
    env = np.diag([sx**2, sp**2, sy**2, sq**2, 6.0e-3**2, 1.0e-3**2])
    c, s = math.cos(t), math.sin(t)
    R = np.eye(6)
    R[0, 0] = R[1, 1] = R[2, 2] = R[3, 3] = c
    R[0, 2] = R[1, 3] = -s
    R[2, 0] = R[3, 1] = s
    env = R @ env @ R.T
    p = {'SLICE': 3, 'NP': 1.0e10, 'SIGZ': 6.0e-3}
    for i in range(6):
        for j in range(i, 6):
            p['R%d%d' % (i + 1, j + 1)] = env[i, j]
    return p, R


def _kick(p, x):
    bb = beambeam.setup(p, 4.0e9)
    S = np.zeros((1, 9))
    S[0, :6] = x
    beambeam.apply(S, 1, bb)
    return S[0, :6].copy()


X0 = np.array([1.5e-4, 2.0e-5, 3.0e-5, -1.0e-5, 1.0e-3, 0.0])


def _expected(t):
    p0, _ = _params(0.0)
    _, R = _params(t)
    return R @ _kick(p0, R.T @ X0)


def test_b34_tilt_rotation_is_covariant(v2_mode):
    t = 0.3
    p, _ = _params(t)
    got = _kick(p, X0)
    exp = _expected(t)
    d0 = np.abs(exp - X0)[[1, 3]]
    assert d0.max() > 1e-8                      # 蹴りがある
    assert np.allclose(got[:4], exp[:4], rtol=0, atol=1e-9 * d0.max() + 1e-18)
    assert abs(got[5] - exp[5]) < 1e-12


def test_b34_k64_mode_is_not_a_rotation(k64_mode):
    t = 0.3
    p, _ = _params(t)
    got = _kick(p, X0)
    exp = _expected(t)
    assert abs(got[0] - exp[0]) > 1e-6         # K64: 座標がつぶれる (cost=sin t, sint=t)
