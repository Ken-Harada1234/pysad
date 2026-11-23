# -*- coding: utf-8 -*-
"""領域と座標の切替 (09_solenoid.md §8、11_misalignment_and_coord.md §7) の式の確認。

K64 の tsol.f の `trots` (SOL 領域の端)、`tsolrot` (要素の DX・DY・DZ・CHI1・CHI2・CHI3) の枝を numpy に写し、
(a) trots: 直線を回して新しい面まで延ばす (直線と平面の交点を 3 次元で独立に解く) と一致するか、
(b) tsolrot (bz = 0): 要素の中央を支点とした剛体回転 + 面への戻し (3 次元) と一致するか (入口・出口)、
(c) tsolrot (bz != 0、傾いたソレノイド): 一様磁場の中のらせん軌道の ODE と一致するか。
(c) では、Fortran の式をそのまま写すと、位置・運動量の更新の dcosphi の符号が面の条件 (s の式) と食い違う結果になる
(xsincos の dc = cos - 1 のとき)。移し間違いか K64 の不整合かは K64 を実行していないので未確認。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
rng = np.random.default_rng(7)


def xsincos(x):
    s = math.sin(x)
    return s, x - s, math.cos(x), math.cos(x) - 1.0     # tfloor.f: s, xs = x - sin x, c, dc = cos x - 1


def pz_of(px, py):
    return math.sqrt(1 - px * px - py * py)


# ------------------------------------------------------------------ (a) trots
def trots(x, px, y, py, z, R, dx, dy, dz, ent, dv=0.0):
    pzi = pz_of(px, py)
    if ent:
        xi, yi = x, y
    else:
        xi, yi = x - dx, y - dy
    xf = R[0, 0] * xi + R[0, 1] * yi
    yf = R[1, 0] * xi + R[1, 1] * yi
    zf = R[2, 0] * xi + R[2, 1] * yi
    pxf = R[0, 0] * px + R[0, 1] * py + R[0, 2] * pzi
    pyf = R[1, 0] * px + R[1, 1] * py + R[1, 2] * pzi
    pzf = R[2, 0] * px + R[2, 1] * py + R[2, 2] * pzi
    if ent:
        return xf - pxf / pzf * zf + dx, pxf, yf - pyf / pzf * zf + dy, pyf, z + zf / pzf + dz
    return xf - pxf / pzf * zf, pxf, yf - pyf / pzf * zf, pyf, z + zf / pzf + dz - dv * dz


def rot(chi1, chi2, chi3):
    c1, s1, c2, s2, c3, s3 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2), math.cos(chi3), math.sin(chi3)
    R1 = np.array([[c1, 0, s1], [0, 1, 0], [-s1, 0, c1]])
    R2 = np.array([[1, 0, 0], [0, c2, s2], [0, -s2, c2]])
    R3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
    return R3 @ R2 @ R1


def line_plane(r, u, plane_z):
    """点 r を通り向き u の直線と、平面 z = plane_z の交点と、r からの媒介変数 λ (r + λ u) を返す"""
    lam = (plane_z - r[2]) / u[2]
    return r + lam * u, lam


print('== (a) trots (SOL 領域の端): 直線を回して新しい面 z = 0 まで延ばす ==')
mx = 0.0
for k in range(5):
    x, y = rng.normal(0, 2e-3, 2)
    px, py = rng.normal(0, 3e-2, 2)
    R = rot(*rng.normal(0, 0.05, 3))
    dx, dy, dz = rng.normal(0, 1e-3, 3)
    # 独立: r = (x, y, 0)、方向 u = (px, py, pz)。回して直線と面 z = 0 の交点を求める。z は面までの距離 -λ (面が前方にあれば z は増える)
    r = R @ np.array([x, y, 0.0])
    u = R @ np.array([px, py, pz_of(px, py)])
    pt, lam = line_plane(r, u, 0.0)
    ref = np.array([pt[0] + dx, u[0], pt[1] + dy, u[1], 0.0 - lam + dz])
    got = np.array(trots(x, px, y, py, 0.0, R, dx, dy, dz, True))
    mx = max(mx, np.abs(ref - got).max())
print('  入口: Fortran の式と 3 次元の直線・平面の交点の差の最大 = %.1e (z の増分は交点までの媒介変数 -λ + dz)' % mx)


# ------------------------------------------------------------------ (b) tsolrot bz = 0
def tsolrot0(x, px, y, py, z, alg, dx, dy, dz, chi1, chi2, ent):
    schi1, _, cchi1, dcchi1 = xsincos(chi1)
    schi2, _, cchi2, dcchi2 = xsincos(chi2)
    if ent:
        s0s = -alg
        x -= dx
        y -= dy
        pz0 = pz_of(px, py)
        pz1 = schi1 * px + cchi1 * pz0
        px = cchi1 * px - schi1 * pz0
        py = cchi2 * py - schi2 * pz1
        ds1 = schi1 * x + dcchi1 * s0s - cchi1 * dz
        x = cchi1 * x - schi1 * (s0s - dz)
        ds2 = schi2 * y + cchi2 * ds1 + dcchi2 * s0s
        y = cchi2 * y - schi2 * (ds1 + s0s)
        pz2 = pz_of(px, py)
        a = ds2 / pz2
        return x - a * px, px, y - a * py, py, z + a
    s0s = alg
    pz0 = pz_of(px, py)
    pz1 = -schi2 * py + cchi2 * pz0
    pyi = cchi2 * py + schi2 * pz0
    pz2 = -schi1 * px + cchi1 * pz1
    pxi = cchi1 * px + schi1 * pz1
    y1 = cchi2 * y - schi2 * s0s
    ds1 = -schi2 * y - dcchi2 * s0s
    ds2 = -schi1 * x + cchi1 * ds1 - dcchi1 * s0s + dz
    x1 = cchi1 * x + schi1 * (ds1 - s0s)
    a = ds2 / pz2
    return x1 - a * pxi + dx, pxi, y1 - a * pyi + dy, pyi, z + a


def _rmat(chi1, chi2):
    c1, s1, c2, s2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    R1p = np.array([[c1, 0, -s1], [0, 1, 0], [s1, 0, c1]])        # tsolrot の第 1 の回転 (x' = c1 x - s1 z)
    R2p = np.array([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])        # 第 2 (y' = c2 y - s2 z)
    return R2p @ R1p


def pivot_ent(x, px, y, py, z, alg, dx, dy, dz, chi1, chi2):
    """要素の中央 (入口の面から前方 alg) を支点として、要素を (dx, dy, dz) ずらし (chi1, chi2) 回す。
    入口の面は支点から見て z = -alg。粒子の位置 (x - dx, y - dy, -alg - dz) を回し、方向も回して、面 z = -alg まで延ばす"""
    R = _rmat(chi1, chi2)
    v = np.array([x - dx, y - dy, -alg - dz])
    u = np.array([px, py, pz_of(px, py)])
    r2 = R @ v
    u2 = R @ u
    pt, lam = line_plane(r2, u2, -alg)
    return np.array([pt[0], u2[0], pt[1], u2[1], z - lam])


def pivot_ext(x, px, y, py, z, alg, dx, dy, dz, chi1, chi2):
    """出口: 要素の座標で (x, y, +alg) (出口の面は支点から前方 alg) の粒子を、逆回転 R^T して (dx, dy, dz) だけずらし、面 z = +alg まで延ばす"""
    RT = _rmat(chi1, chi2).T
    v = RT @ np.array([x, y, alg]) + np.array([dx, dy, dz])
    u = RT @ np.array([px, py, pz_of(px, py)])
    pt, lam = line_plane(v, u, alg)
    return np.array([pt[0], u[0], pt[1], u[1], z - lam])


print('\n== (b) tsolrot (bz = 0): 要素の中央を支点とする剛体回転 + 面への戻し ==')
for ent, fn in ((True, pivot_ent), (False, pivot_ext)):
    mx = 0.0
    for k in range(5):
        x, y = rng.normal(0, 2e-3, 2)
        px, py = rng.normal(0, 3e-2, 2)
        dx, dy, dz = rng.normal(0, 1e-3, 3)
        chi1, chi2 = rng.normal(0, 0.03, 2)
        got = np.array(tsolrot0(x, px, y, py, 0.0, 0.4 if ent else -0.4, dx, dy, dz, chi1, chi2, ent))   # 出口は tmulti が alg - al (= -alg) を渡す
        ref = fn(x, px, y, py, 0.0, 0.4, dx, dy, dz, chi1, chi2)
        mx = max(mx, np.abs(got - ref).max())
    print('  %s: Fortran の枝と 3 次元の構成の差の最大 = %.1e' % ('入口' if ent else '出口', mx))
print('  (z は、面まで前に進んだ分だけ増える。ent では z + ds2/pz2、(a) の z + zf/pzf と同じ形)')


# ------------------------------------------------------------------ (c) tsolrot bz != 0 (ent)
def tsolrot_bz_ent(x, px, y, py, g, bz, alg, dx, dy, dz, chi1, chi2, dc_sign=1):
    s0s = -alg
    schi1, _, cchi1, dcchi1 = xsincos(chi1)
    schi2, _, cchi2, dcchi2 = xsincos(chi2)
    b = abs(bz)
    sg = math.copysign(1, bz)
    phix, phiy, phiz = -schi1 * sg, cchi1 * (-schi2) * sg, cchi1 * cchi2 * sg
    dphizsq = phix**2 + phiy**2
    pr = 1 + g
    px = px + bz * y / pr * .5
    py = py - bz * x / pr * .5
    x -= dx
    y -= dy
    ds1 = schi1 * x + dcchi1 * s0s - cchi1 * dz
    x = cchi1 * x - schi1 * (s0s - dz)
    ds2 = -(schi2 * y + cchi2 * ds1 + dcchi2 * s0s)
    y = cchi2 * y - schi2 * (ds1 + s0s)
    pz0 = pz_of(px, py)
    pz1 = schi1 * px + cchi1 * pz0
    px = cchi1 * px - schi1 * pz0
    py = cchi2 * py - schi2 * pz1
    start = (x, y, s0s - ds2, px, py)
    alb = pr / b
    dpz0 = pz_of(px, py) - 1
    pz0 = 1 + dpz0
    dpl = px * phix + py * phiy + dpz0 * phiz
    pl = phiz + dpl
    plx, ply, plz = pl * phix, pl * phiy, pl * phiz
    ptx, pty, ptz = px - plx, py - ply, dpz0 - dpl * phiz + dphizsq
    pbx, pby, pbz = pty * phiz - ptz * phiy, ptz * phix - ptx * phiz, ptx * phiy - pty * phix
    ps2 = ds2 / alb
    phi = math.asin(ps2 / pz0)
    for _ in range(10):
        sn, xsn, cs, dc = xsincos(phi)
        dcs = dc * dc_sign
        s = plz * xsn + pz0 * sn - pbz * dcs
        dphi = (ps2 - s) / (pz0 + ptz * dcs + pbz * sn)
        phi0 = phi
        phi += dphi
        if phi0 == phi or abs(dphi / phi) < 3e-16:
            break
    sn, xsn, cs, dc = xsincos(phi)
    dcs = dc * dc_sign
    out = dict(x=x + (plx * phi + ptx * sn + pbx * dcs) * alb, y=y + (ply * phi + pty * sn + pby * dcs) * alb,
               px=px - ptx * dcs + pbx * sn, py=py - pty * dcs + pby * sn, l=phi * alb)
    return out, dict(start=start, phi=(phix, phiy, phiz), bp=b / pr, s0s=s0s, ds2=ds2)


def helix_ode(info):
    x, y, z, px, py = info['start']
    phi = np.array(info['phi'])
    bp = info['bp']

    def rhs(l, v):
        return list(v[3:]) + list(bp * np.cross(v[3:], phi))       # dr/dl = u、du/dl = b_p u x phi
    ev = lambda l, v: v[2] - info['s0s']
    ev.terminal = True
    sgn = 1 if info['ds2'] > 0 else -1
    sol = solve_ivp(rhs, (0, sgn * 10), [x, y, z, px, py, pz_of(px, py)], events=ev, rtol=1e-13, atol=1e-16, method='DOP853')
    v = sol.y_events[0][0]
    return dict(x=v[0], y=v[1], px=v[3], py=v[4], l=sol.t_events[0][0])


print('\n== (c) tsolrot (bz != 0、傾いたソレノイド): らせん軌道の ODE との比較 ==')
case = dict(x=0.002, px=0.01, y=-0.001, py=0.004, g=0.01, bz=1.3, alg=0.5, dx=0.0005, dy=-0.0003, dz=0.01, chi1=0.03, chi2=-0.02)
for name, dcs in (('dcosphi = cos - 1 (xsincos の戻り値のまま)', 1), ('dcosphi = 1 - cos (符号を反転)', -1)):
    out, info = tsolrot_bz_ent(dc_sign=dcs, **case)
    ref = helix_ode(info)
    print('  %s:' % name)
    for k in ('x', 'y', 'px', 'py', 'l'):
        print('     %-2s: Fortran の式 = %.12e、ODE = %.12e、差 = %+.2e' % (k, out[k], ref[k], out[k] - ref[k]))
print('  読み: 面の条件の式 (s = plz*xsinphi + pz0*sinphi - pbz*dcosphi) は dcosphi = cos - 1 で ODE と一致する (l の差 1e-17)。')
print('        位置・運動量の更新は dcosphi = 1 - cos のとき ODE に近い。同じ dcosphi で両方は成り立たない。')
