# -*- coding: utf-8 -*-
"""B165 の独立な基準 (その 2、pysad のコードを使わない): 直線の磁石 (QUAD・SEXT・MULT・ステアリング) の Maxwell 場の中の Thomas-BMT の数値積分。

場は スカラーポテンシャル Φ = Σ_n Φ_n、Φ_n = Σ_j c_j G^(2j)(z) r^(2j) P_n(x, y) (c_0 = 1、c_(j+1) = −c_j/(4 (j+1)(j+m+1))、m = n+1)、
P_n = k_n Im(z^m)/m! (z = x + i y。By + i Bx = Σ k_n z^n/n!)、G(z) = f(z) f(L − z)、f(d) = (1 + tanh(d/w))/2 (フリンジ幅 w)。
j = 0, 1, 2 まで入れる (∇²Φ = 0 は j の和で z についての微分の項が打ち消す)。B = ∇Φ (curl B = 0・div B = 0)。
運動方程式は dr/dl = p/pt、dp/dl = p × B/pt (B は Bρ0 で規格化、p は p0 で規格化)、スピンは dS/dl = S × Ω、Ω = [(1 + aγ) B⊥ + (1 + a) B∥]/pt
(符号は _spin_maxwell_indep.py と同じ)。入口の面 z = 0 での座標 (x0, y0, px0, py0) から出口の面 z = L の外の直線まで積分し、スピン回転行列 R (3x3) と出口の面の座標を返す。"""
import math

import numpy as np
from scipy.integrate import solve_ivp

GS = 0.001159652181280002
ME = 0.51099895e6


def _tanh_polys(kmax):
    """tanh の k 階微分を t = tanh の多項式として (numpy.poly1d の係数)"""
    polys = [np.poly1d([1.0, 0.0])]
    for _ in range(kmax):
        p = polys[-1]
        polys.append(p.deriv() * np.poly1d([-1.0, 0.0, 1.0]))
    return polys


_TP = _tanh_polys(6)


def _fk(k, d, w):
    t = math.tanh(d / w)
    if k == 0:
        return 0.5 * (1.0 + t)
    return 0.5 * _TP[k](t) / w ** k


def _Gk(k, z, L, w):
    """G(z) = f(z) f(L − z) の k 階微分"""
    s = 0.0
    for i in range(k + 1):
        s += math.comb(k, i) * _fk(i, z, w) * (-1) ** (k - i) * _fk(k - i, L - z, w)
    return s


def make_field(L, mult, w, jmax=2):
    """mult: {n: k_n} (n = 0: 二極 By = k0、1: 四極、2: 六極 …)。B(x, y, z) を返す"""
    terms = []
    for n, kn in mult.items():
        m = n + 1
        a = kn / math.factorial(m)
        c = [1.0]
        for j in range(jmax):
            c.append(-c[-1] / (4.0 * (j + 1) * (j + m + 1)))
        terms.append((m, a, c))

    def B(x, y, z):
        Gs = [_Gk(k, z, L, w) for k in range(2 * jmax + 2)]
        bx = by = bz = 0.0
        zc = complex(x, y)
        r2 = x * x + y * y
        for (m, a, c) in terms:
            zm = zc ** m
            P = a * zm.imag
            Px = a * m * (zc ** (m - 1)).imag
            Py = a * m * (zc ** (m - 1)).real
            for j in range(jmax + 1):
                g, g1 = Gs[2 * j], Gs[2 * j + 1]
                rj = r2 ** j
                if j > 0:
                    rjm = r2 ** (j - 1)
                    bx += c[j] * g * (2 * j * x * rjm * P + rj * Px)
                    by += c[j] * g * (2 * j * y * rjm * P + rj * Py)
                else:
                    bx += c[j] * g * Px
                    by += c[j] * g * Py
                bz += c[j] * g1 * rj * P
        return np.array([bx, by, bz])
    return B


def maxwell_R(Ep, L, mult, coord, dp=0.0, w=4e-3, zs=None, jmax=2, a_eff=None, dx=0.0, dy=0.0, rot=0.0, z0=0.0):
    """Ep: 運動量 eV/c。coord = (x0, px0, y0, py0) (px, py は p0 で規格化。入口の面 z = 0)。R (3x3) と出口の面 z = L での (x, px, y, py)"""
    x0, px0, y0, py0 = coord
    if zs is None:
        zs = 8.0 * w
    pt = 1.0 + dp
    gam = math.sqrt(1.0 + (Ep / ME * pt) ** 2)
    B0 = make_field(L, mult, w, jmax)
    if dx != 0.0 or dy != 0.0 or rot != 0.0:
        # 磁石の座標 (ずらして rot だけ回した: x' = cos(rot)(x − dx) − sin(rot)(y − dy)、y' = sin(rot)(x − dx) + cos(rot)(y − dy)) の場を実験室系に戻す
        cr, sr = math.cos(rot), math.sin(rot)

        def B(x, y, z):
            xm, ym = x - dx, y - dy
            bx, by, bz = B0(cr * xm - sr * ym, sr * xm + cr * ym, z - z0)
            return np.array([cr * bx + sr * by, -sr * bx + cr * by, bz])
    else:
        B = B0
    pz0 = math.sqrt(pt * pt - px0 * px0 - py0 * py0)

    def rhs(l, u):
        r, pp, S = u[:3], u[3:6], u[6:].reshape(3, 3)
        Bv = B(r[0], r[1], r[2])
        ph = pp / pt
        Bl = (Bv @ ph) * ph
        Om = ((1 + GS * gam) * (Bv - Bl) + (1 + GS) * Bl) / pt
        dS = np.cross(S, Om[None, :])
        return np.concatenate([ph, np.cross(pp, Bv) / pt, dS.ravel()])
    ev = lambda l, u: u[2] - (L + zs)
    ev.terminal = True
    r0 = np.array([x0 - px0 / pz0 * zs, y0 - py0 / pz0 * zs, -zs])
    u0 = np.concatenate([r0, [px0, py0, pz0], np.eye(3).ravel()])
    sol = solve_ivp(rhs, (0, 4 * L + 8 * zs + 1), u0, events=ev, rtol=1e-12, atol=1e-15, method='DOP853', max_step=w / 3)
    uf = sol.y_events[0][0]
    R = uf[6:].reshape(3, 3).T
    # 出口の面 z = L の外の直線で (x, y) を z = L まで戻す
    pp = uf[3:6]
    dz = (L + zs) - L
    xf = uf[0] - pp[0] / pp[2] * dz
    yf = uf[1] - pp[1] / pp[2] * dz
    return R, np.array([xf, pp[0], yf, pp[1]])
