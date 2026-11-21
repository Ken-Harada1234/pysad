# -*- coding: utf-8 -*-
"""v2.4.0: BEAMBEAM の傾き (B34) の独立な確認。強ビームの x-y 結合のある (傾いた) 2 次元ガウス分布の場を、極座標の数値積分
(試験粒子のまわりで 1/|u| の特異性が消える: F = −∫dθ û ∫0^∞ ρ(r+uû) du) で作り、`beambeam.apply` の蹴りと比べる。
条件: SLICE=1 (z=0 の 1 枚、試験粒子も z=0)、強ビームは 4 GeV の電子 1e10 個、σx 2e-4・σy 2e-5・σz 6 mm の包絡を x-y 面内で角 t 回したもの (px・py も同じ回転)。
v2 は 1e-4 以内で一致 (傾き 0・0.3・0.7・−0.5 rad)。k64 は傾き 0 では一致し、傾きがあると一致しない (B34)。
独立な数値積分は SAD・pysad のコードを使わない。
"""
import math

import numpy as np
import pytest

from pysad import beambeam

RE = 2.8179403262e-15
ME = 0.51099895e6
SX, SY = 2.0e-4, 2.0e-5


def _env(t):
    e = np.diag([SX**2, 1.0e-4**2, SY**2, 5.0e-5**2, 6.0e-3**2, 1.0e-3**2])
    c, s = math.cos(t), math.sin(t)
    R = np.eye(6)
    R[0, 0] = R[1, 1] = R[2, 2] = R[3, 3] = c
    R[0, 2] = R[1, 3] = -s
    R[2, 0] = R[3, 1] = s
    return R @ e @ R.T


def _kick_pysad(t, x, y):
    env = _env(t)
    p = {'SLICE': 1, 'NP': 1.0e10, 'SIGZ': 6.0e-3}
    for i in range(6):
        for j in range(i, 6):
            p['R%d%d' % (i + 1, j + 1)] = env[i, j]
    bb = beambeam.setup(p, 4.0e9)
    S = np.zeros((1, 9))
    S[0, :6] = [x, 0.0, y, 0.0, 0.0, 0.0]
    beambeam.apply(S, 1, bb)
    return S[0, 1], S[0, 3]


def _field_numeric(t, x, y):
    """単位電荷の傾いた 2 次元ガウス分布 (共分散 = _env(t) の x-y 部分) の 2 次元クーロン場 ∫ρ(r')(r−r')/|r−r'|² d²r'"""
    cov = _env(t)[np.ix_([0, 2], [0, 2])]
    ci = np.linalg.inv(cov)
    norm = 1.0 / (2 * math.pi * math.sqrt(np.linalg.det(cov)))
    r = np.array([x, y])
    xg, wg = np.polynomial.legendre.leggauss(400)
    th = (xg + 1) * math.pi
    wth = wg * math.pi
    umax = 12.0 * SX
    ug, wu = np.polynomial.legendre.leggauss(600)
    u = (ug + 1) * 0.5 * umax
    wu = wu * 0.5 * umax
    F = np.zeros(2)
    for a, wa in zip(th, wth):
        uh = np.array([math.cos(a), math.sin(a)])
        pts = r[None, :] + u[:, None] * uh[None, :]
        q = np.einsum('ni,ij,nj->n', pts, ci, pts)
        rho = norm * np.exp(-0.5 * q)
        F += -uh * float(np.sum(rho * wu)) * wa
    return F


@pytest.mark.parametrize('t', [0.0, 0.3, 0.7, -0.5])
def test_tilted_beambeam_kick_equals_numeric_field_v2(v2_mode, t):
    x, y = 3.0e-4, 4.0e-5
    kx, ky = _kick_pysad(t, x, y)
    F = _field_numeric(t, x, y)
    ce = 1.0e10 * RE / (4.0e9 / ME)
    exp = -2.0 * ce * F
    assert abs(kx / exp[0] - 1) < 2e-4 or abs(kx / exp[0] + 1) < 2e-4
    assert abs(ky / exp[1] - 1) < 2e-4 or abs(ky / exp[1] + 1) < 2e-4
    s = np.sign(kx / exp[0])
    assert abs(kx - s * exp[0]) < 2e-4 * abs(exp[0]) + 1e-15 and abs(ky - s * exp[1]) < 2e-4 * abs(exp[1]) + 1e-15


def test_untilted_matches_in_k64_but_tilted_does_not(k64_mode):
    x, y = 3.0e-4, 4.0e-5
    ce = 1.0e10 * RE / (4.0e9 / ME)
    kx0, ky0 = _kick_pysad(0.0, x, y)
    e0 = -2.0 * ce * _field_numeric(0.0, x, y)
    assert abs(abs(kx0 / e0[0]) - 1) < 2e-4
    kx, ky = _kick_pysad(0.3, x, y)
    e1 = -2.0 * ce * _field_numeric(0.3, x, y)
    assert abs(abs(kx / e1[0]) - 1) > 1e-2 or abs(abs(ky / e1[1]) - 1) > 1e-2
