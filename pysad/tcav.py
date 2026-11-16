"""
TCAVI 要素 (横偏向空洞・クラブ空洞)。K64 の ttcav.f (追跡)、ttcave.f (6 次元の写像と行列)。
RFSW OFF では長さ L のドリフト (tturn、ttcave とも)。追跡は RANKICK / RANPHASE (tgauss) 未対応。

ttcav.f の TENT / TEXIT (ずれ DX DY と回転 ROTATE) を入口と出口に置く。TENT の回転は
  x=c·x−s·y; y=(s·x'+y)/c (|c|≥|s|) の順 (更新した x を使う) で、そのまま写す。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import track6n as _n
from .emit6d import CVELOC

NJ = dict(cache=True, error_model='numpy')


def params(el):
    p = el.params
    g = lambda k: float(p.get(k, 0.0) or 0.0)
    return dict(al=float(el.L), ak=g('K0'), harm=g('HARM'), phi=g('PHI'), freq=g('FREQ'), dx=g('DX'), dy=g('DY'),
                theta=g('ROTATE'), rank=g('RANKICK'), ranp=g('RANPHASE'), rad=g('RAD') == 0.0)


@njit(**NJ)
def _tent(S, n, dx, dy, theta):
    if theta != 0.0:
        cost = math.cos(theta)
        sint = math.sin(theta)
    else:
        cost = 1.0
        sint = 0.0
    for i in range(n):
        if dx != 0.0:
            S[i, 0] -= dx
        if dy != 0.0:
            S[i, 2] -= dy
        if theta != 0.0:
            if abs(cost) >= abs(sint):
                S[i, 0] = cost * S[i, 0] - sint * S[i, 2]
                S[i, 2] = (sint * S[i, 0] + S[i, 2]) / cost
                S[i, 1] = cost * S[i, 1] - sint * S[i, 3]
                S[i, 3] = (sint * S[i, 1] + S[i, 3]) / cost
            else:
                xi = S[i, 0]
                S[i, 0] = xi * cost - S[i, 2] * sint
                S[i, 2] = xi * sint + S[i, 2] * cost
                pxi = S[i, 1]
                S[i, 1] = pxi * cost - S[i, 3] * sint
                S[i, 3] = pxi * sint + S[i, 3] * cost


@njit(**NJ)
def _texit(S, n, dx, dy, theta):
    if theta != 0.0:
        cost = math.cos(theta)
        sint = math.sin(theta)
    else:
        cost = 1.0
        sint = 0.0
    for i in range(n):
        if theta != 0.0:
            if abs(cost) >= abs(sint):
                S[i, 0] = cost * S[i, 0] + sint * S[i, 2]
                S[i, 2] = (-sint * S[i, 0] + S[i, 2]) / cost
                S[i, 1] = cost * S[i, 1] + sint * S[i, 3]
                S[i, 3] = (-sint * S[i, 1] + S[i, 3]) / cost
            else:
                xi = S[i, 0]
                S[i, 0] = xi * cost + S[i, 2] * sint
                S[i, 2] = -xi * sint + S[i, 2] * cost
                pxi = S[i, 1]
                S[i, 1] = pxi * cost + S[i, 3] * sint
                S[i, 3] = -pxi * sint + S[i, 3] * cost
        if dy != 0.0:
            S[i, 2] += dy
        if dx != 0.0:
            S[i, 0] += dx


@njit(**NJ)
def k_tcav(S, n, al, ak, w, phic, p0, h0, krad, fluc, rs, anrad, cuc):
    """ttcav.f 本体 (TENT/TEXIT は呼び出し側)。w: 波数、phic = (PHI + 周回の位相)·sign(charge)"""
    if al != 0.0:
        _n.k_drift(S, n, al * 0.5)
        if krad:
            for i in range(n):
                _n.setr0_1(i, S)
    v = ak * p0 * w
    for i in range(n):
        dp1r = S[i, 5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        if S[i, 6] > 0.1:
            h1 = _n.p2h(p1)
        else:
            h1 = p1r * h0 / (1.0 - S[i, 6] + S[i, 10])
        t = min(1.0e4, max(-1.0e4, -S[i, 4] * h1 / p1))
        ph = w * t + phic
        dh = v * S[i, 0] * math.cos(ph)
        S[i, 1] = S[i, 1] - ak * math.sin(ph) / p1r
        h2 = h1 + dh
        a = max(-1.0, dh * (h1 + h2) / p1**2)
        dpr = a * (0.5 - a * (0.125 - a * 0.0625))
        dpr = (dpr**2 + a) / (2.0 + 2.0 * dpr)
        dpr = (dpr**2 + a) / (2.0 + 2.0 * dpr)
        dp2r = dp1r + p1r * dpr
        p2r = 1.0 + dp2r
        S[i, 5] = dp2r
        S[i, 6] = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + S[i, 10]
        S[i, 1] = S[i, 1] * p1r / p2r
        S[i, 3] = S[i, 3] * p1r / p2r
        S[i, 4] = -t * p2r * p0 / h2
    if al != 0.0:
        _n.k_drift(S, n, al * 0.5)
        if krad:
            _n.tradk_all(S, n, al, 0.0, fluc, rs, p0, anrad, cuc)


def compile_tcav(el):
    pr = params(el)
    if pr['rank'] != 0.0 or pr['ranp'] != 0.0:
        raise NotImplementedError('track6n: TCAVI の RANKICK / RANPHASE (tgauss) は未対応')
    return pr


def apply_tcav(S, n, pr, ph, rad, fluc, rs, kturn):
    """tturn.f icTCAV (RFSW ON のとき)。ph: _Phys"""
    harm = pr['harm']
    w = 2.0 * math.pi * pr['freq'] / CVELOC if harm == 0.0 else ph.omega0 * harm / CVELOC
    harmf = harm - int(harm)
    phase = pr['phi'] + harmf * (kturn - 1) * 2.0 * math.pi
    phic = phase * (1.0 if ph.charge >= 0 else -1.0)
    _tent(S, n, pr['dx'], pr['dy'], pr['theta'])
    k_tcav(S, n, pr['al'], pr['ak'], w, phic, ph.p0, ph.h0, bool(rad and pr['rad']), fluc, rs, ph.anrad, ph.cuc)
    _texit(S, n, pr['dx'], pr['dy'], pr['theta'])


def tcav_kick(pr, ph):
    """ttcave.f のキック: (COD 写像 f(c), 解析行列 J(c))。RFSW ON 用 (t=-z/v1)"""
    harm = pr['harm']
    w = 2.0 * math.pi * pr['freq'] / CVELOC if harm == 0.0 else ph.omega0 * harm / CVELOC
    ak = pr['ak']
    p0, h0 = ph.p0, ph.h0
    v = ak * p0 * w
    phic = pr['phi'] * (1.0 if ph.charge >= 0 else -1.0)

    def _core(c):
        p1 = p0 * (1.0 + c[5])
        h1 = p1 * math.sqrt(1.0 + 1.0 / p1**2)
        h1 = p1 + 1.0 / (h1 + p1)
        dh1 = p0 * (p0 + p1) / (h1 + h0) * c[5]
        v1 = p1 / h1
        t = -c[4] / v1
        phii = w * t + phic
        cosp = math.cos(phii)
        sinp = math.sin(phii)
        dh = v * c[0] * cosp
        dh2 = dh1 + dh
        h2 = h0 + dh2
        p2 = h2 * math.sqrt(1.0 - 1.0 / h2**2)
        pf = (h2 + h0) / (p0 + p2) * dh2 / p0
        p2 = p0 * (1.0 + pf)
        v2 = p2 / h2
        return p1, h1, v1, t, cosp, sinp, dh, h2, p2, pf, v2

    def f(c):
        p1, h1, v1, t, cosp, sinp, dh, h2, p2, pf, v2 = _core(c)
        out = np.array(c, dtype=float)
        out[5] = pf
        out[4] = -t * v2
        out[1] = c[1] - ak * sinp
        return out

    def jac(c):
        p1, h1, v1, t, cosp, sinp, dh, h2, p2, pf, v2 = _core(c)
        a = -v * w * c[0] * sinp
        J = np.eye(6)
        J[1, 4] = ak * w / v1 * cosp
        J[1, 5] = ak * w * cosp * t * p0 / p1 / h1**2
        J[4, 0] = -t * p0 / p2 / h2**2 * ak * w * cosp
        J[4, 4] = (p2 + a * t / p2 / h2) / h2 / v1
        J[4, 5] = t * p0 / h1**2 / h2**2 / p1 / p2 * (dh * (h2 * (h2 + h1) + p1**2) + a * t)
        J[5, 0] = ak * w * cosp / v2
        J[5, 4] = -a / v1 / v2 / p0
        J[5, 5] = (p1 - a * t / p1 / h1) / h1 / v2
        return J
    return f, jac
