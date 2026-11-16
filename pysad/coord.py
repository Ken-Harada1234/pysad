"""
COORD 要素 (座標変換: ずれ DX DY DZ と回転 CHI1 CHI2 CHI3)。K64 の tcoord.f (追跡)、tcoorde.f / qcoord.f (6 次元の写像と行列)。
幾何 (LINE["GX"] など) は geo.py `_coord`。

DIR=0 (既定): 入口で DX DY DZ だけずらしてから回転。DIR≠0: 逆変換 (回転してから出口で戻す)。
追跡は px,py が斜率 (px/(1+δ))、6 次元 (emit/Twiss) は正準運動量で、z を持つ y 軸まわりの平面 z=0 へ戻す。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

NJ = dict(cache=True, error_model='numpy')
_KEYS = ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')


def params(el):
    """要素から (dx, dy, dz, chi1, chi2, chi3, dir0) (dir0: DIR==0 → 順変換)"""
    p = el.params
    g = lambda k: float(p.get(k, 0.0) or 0.0)
    return g('DX'), g('DY'), g('DZ'), g('CHI1'), g('CHI2'), g('CHI3'), g('DIR') == 0.0


@njit(**NJ)
def _rot(dx, dy, dz, chi1, chi2, chi3, direct):
    """(r11..r33, dx0, dy0, dz0, dx1, dy1, dz1)"""
    c1 = math.cos(chi1); s1 = math.sin(chi1)
    c2 = math.cos(chi2); s2 = math.sin(chi2)
    c3 = math.cos(chi3); s3 = math.sin(chi3)
    if direct:
        dx0 = dx; dy0 = dy; dz0 = dz
        dx1 = 0.0; dy1 = 0.0; dz1 = 0.0
        r11 = c1 * c3 + s1 * s2 * s3
        r12 = -c2 * s3
        r13 = s1 * c3 - c1 * s2 * s3
        r21 = -s1 * s2 * c3 + c1 * s3
        r22 = c2 * c3
        r23 = c1 * s2 * c3 + s1 * s3
        r31 = -s1 * c2
        r32 = -s2
        r33 = c1 * c2
    else:
        dx0 = 0.0; dy0 = 0.0; dz0 = 0.0
        dx1 = dx; dy1 = -dy; dz1 = -dz
        r11 = c1 * c3 + s1 * s2 * s3
        r21 = c2 * s3
        r31 = -s1 * c3 + c1 * s2 * s3
        r12 = s1 * s2 * c3 - c1 * s3
        r22 = c2 * c3
        r32 = c1 * s2 * c3 + s1 * s3
        r13 = s1 * c2
        r23 = -s2
        r33 = c1 * c2
    return r11, r12, r13, r21, r22, r23, r31, r32, r33, dx0, dy0, dz0, dx1, dy1, dz1


@njit(**NJ)
def _sqrt1m(a):
    """mathfun sqrt1(-a) 相当: sqrt(1-a)-1 (a が小さいときは級数)"""
    if a < 8e-6:
        return -a * (0.5 + a * (0.125 + a * 0.0625))
    return -a / (1.0 + math.sqrt(max(1e-200, 1.0 - a)))


@njit(**NJ)
def k_coord(S, n, dx, dy, dz, chi1, chi2, chi3, direct):
    """tcoord.f (追跡: S の x,px,y,py,z。px,py は斜率)"""
    (r11, r12, r13, r21, r22, r23, r31, r32, r33, dx0, dy0, dz0, dx1, dy1, dz1) = _rot(dx, dy, dz, chi1, chi2, chi3, direct)
    for i in range(n):
        pxi = S[i, 1]
        pyi = S[i, 3]
        a = pxi * pxi + pyi * pyi
        dpz = _sqrt1m(a)
        pzi = 1.0 + dpz
        zi = pzi * S[i, 4] - dz0
        xi = S[i, 0] + pxi * S[i, 4] - dx0
        yi = S[i, 2] + pyi * S[i, 4] - dy0
        xf = r11 * xi + r12 * yi + r13 * zi
        yf = r21 * xi + r22 * yi + r23 * zi
        zf = r31 * xi + r32 * yi + r33 * zi + dz1
        px = r11 * pxi + r12 * pyi + r13 * pzi
        py = r21 * pxi + r22 * pyi + r23 * pzi
        pzf = r31 * pxi + r32 * pyi + r33 * pzi
        z = zf / pzf
        S[i, 1] = px
        S[i, 3] = py
        S[i, 4] = z
        S[i, 0] = xf - px * z + dx1
        S[i, 2] = yf - py * z + dy1


def coord_map(c, prm):
    """tcoorde の COD 写像と 6×6 の行列 (T, c_out)。c = (x, px, y, py, z, δ) (px,py は正準運動量)"""
    dx, dy, dz, chi1, chi2, chi3, direct = prm
    (r11, r12, r13, r21, r22, r23, r31, r32, r33, dx0, dy0, dz0, dx1, dy1, dz1) = _rot(dx, dy, dz, chi1, chi2, chi3, direct)
    p = 1.0 + c[5]
    z = c[4]
    pxi = c[1]
    pyi = c[3]
    a = pxi**2 + pyi**2
    pzi = p * math.sqrt(max(1.0 - a / p**2, 1e-200))
    zi = pzi / p * z - dz0
    xi = c[0] - dx0 + pxi / p * z
    yi = c[2] - dy0 + pyi / p * z
    T1 = np.eye(6)
    T1[0, 1] = z / p
    T1[0, 4] = pxi / p
    T1[0, 5] = -pxi * z / p**2
    T1[2, 3] = z / p
    T1[2, 4] = pyi / p
    T1[2, 5] = -pyi * z / p**2
    T1[4, 1] = -pxi / p / pzi * z
    T1[4, 3] = -pyi / p / pzi * z
    T1[4, 4] = pzi / p
    T1[4, 5] = (pxi**2 + pyi**2) / p**2 / pzi * z
    T1[5, 1] = -pxi / pzi
    T1[5, 3] = -pyi / pzi
    T1[5, 5] = p / pzi
    xf = r11 * xi + r12 * yi + r13 * zi
    yf = r21 * xi + r22 * yi + r23 * zi
    zf = r31 * xi + r32 * yi + r33 * zi + dz1
    pxf = r11 * pxi + r12 * pyi + r13 * pzi
    pyf = r21 * pxi + r22 * pyi + r23 * pzi
    pzf = r31 * pxi + r32 * pyi + r33 * pzi
    R = np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])
    A = T1.copy()
    for rows in ((0, 2, 4), (1, 3, 5)):            # (x, y, z) と (px, py, pz) の行を R で混ぜる
        A[list(rows), :] = R @ T1[list(rows), :]
    T2 = np.eye(6)
    T2[0, 1] = -zf / pzf
    T2[0, 4] = -pxf / pzf
    T2[0, 5] = pxf / pzf**2 * zf
    T2[2, 3] = -zf / pzf
    T2[2, 4] = -pyf / pzf
    T2[2, 5] = pyf / pzf**2 * zf
    T2[4, 1] = pxf / p / pzf * zf
    T2[4, 3] = pyf / p / pzf * zf
    T2[4, 4] = p / pzf
    T2[4, 5] = -(pxf**2 + pyf**2) / p / pzf**2 * zf
    T2[5, 1] = pxf / p
    T2[5, 3] = pyf / p
    T2[5, 5] = pzf / p
    T = T2 @ A
    out = np.array([xf - pxf / pzf * zf + dx1, pxf, yf - pyf / pzf * zf + dy1, pyf, p / pzf * zf, c[5]])
    return T, out
