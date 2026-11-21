# -*- coding: utf-8 -*-
"""v2.3.0 U2 の試験・例の補助: ドリフト・四極 (厚い)・扇形偏向磁石 (複合機能、K1 の定義 'log'|'lin') のリングの独立な数値積分。
pysad・SAD のコードを使わない。_u1_indep.py (実験室系の円筒座標のローレンツ力) の sector_map を偏向磁石に使い、ドリフトと四極も厳密な運動方程式で積分する。
Bρ0 = 1 の単位、運動量 p = 1 + δ、ux・uy は方向余弦 (SAD の px/p・py/p)。
 edge=True なら扇形磁石の端の縦方向の集束のインパルス (_vedge) を入れる (K64・elegant の硬い端の標準の模型に相当)。
 ドリフト: dx/dz = ux/uz、dy/dz = uy/uz、四極 (B_y = K1 x、B_x = K1 y): dux/dz = -K1 x/p、duy/dz = +K1 y/p。端は硬い端・フリンジなし。
リングの要素の並び ring = [('D', L), ('Q', L, K1), ('B', model, rho0, theta0, K1), ...]。閉軌道 (δ ごと) の周りの線形写像から (νx, νy)。
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

import _u1_indep as U1


def _drift(s, L):
    x, ux, y, uy = s
    uz = math.sqrt(1.0 - ux * ux - uy * uy)
    return np.array([x + L * ux / uz, ux, y + L * uy / uz, uy])


def _quad(s, L, k1, p, rtol=1e-13, atol=1e-15):
    def rhs(z, v):
        x, ux, y, uy = v
        uz = math.sqrt(1.0 - ux * ux - uy * uy)
        return [ux / uz, -k1 * x / p, uy / uz, k1 * y / p]
    sol = solve_ivp(rhs, (0.0, L), list(s), method='DOP853', rtol=rtol, atol=atol)
    return sol.y[:, -1]


def _vedge(s, f, rho0, p, sgn):
    """扇形磁石の端 (硬い端の標準の縦方向の集束): 端面 (磁石の中心軌道に垂直) の縁の磁場 B_z = y ∂B_y/∂z (∇×B = 0) を軌道の傾き ux/uz で横切るインパルス。
    入口 Δuy = −(ux/uz) g y / p、出口 +。g は端の位置の中央面の磁場 B_y (Bρ0 = 1)。y の 1 次まで (y³ は落とす)。"""
    x, ux, y, uy = s
    uz = math.sqrt(1.0 - ux * ux - uy * uy)
    g = f(rho0 + x, 0.0)[0]
    return np.array([x, ux, y, uy - sgn * (ux / uz) * g * y / p])


def one_turn(ring, s, p, fields, edge=False):
    for e in ring:
        if e[0] == 'D':
            s = _drift(s, e[1])
        elif e[0] == 'Q':
            s = _quad(s, e[1], e[2], p)
        else:
            _, model, rho0, th, k1 = e
            f = fields.get((model, rho0, k1))
            if f is None:
                f = fields[(model, rho0, k1)] = U1.make_field(model, rho0, k1)
            if edge:
                s = _vedge(s, f, rho0, p, +1.0)
            s = U1.sector_map(f, rho0, th, p, s)
            if edge:
                s = _vedge(s, f, rho0, p, -1.0)
    return s


def closed_orbit(ring, delta, fields, edge=False):
    p = 1.0 + delta
    # 水平の閉軌道 (y = 0)。初期値は分散 ≈ 数 m × δ 程度から始め、fsolve で求める
    def g(v):
        s = one_turn(ring, np.array([v[0], v[1], 0.0, 0.0]), p, fields, edge)
        return [s[0] - v[0], s[1] - v[1]]
    v = fsolve(g, [0.0, 0.0], xtol=1e-14)
    return np.array([v[0], v[1], 0.0, 0.0])


def tunes(ring, delta, eps=1e-7, fields=None, edge=False):
    fields = {} if fields is None else fields
    p = 1.0 + delta
    s0 = closed_orbit(ring, delta, fields, edge)
    m = np.zeros((4, 4))
    for j in range(4):
        d = np.zeros(4)
        d[j] = eps
        m[:, j] = (one_turn(ring, s0 + d, p, fields, edge) - one_turn(ring, s0 - d, p, fields, edge)) / (2 * eps)
    out = []
    for a, b in ((0, 1), (2, 3)):
        tr = 0.5 * (m[a, a] + m[b, b])
        out.append(math.acos(tr) / (2 * math.pi))      # 端数 (整数部分は呼び出し側で決める)
    return out[0], out[1]
