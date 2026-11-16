"""
GaussianCoulomb[x, y, sigx, sigy] / GaussianCoulombU[x, y, sigx, sigy]。K64 の tspac.f (bbkick1、cwerf)、tintrb.f (twspu)、bint.f (rombint、polint)、gamma.f (gamma0log)。

GaussianCoulomb: ガウス分布のバンチが (x, y) に作る力 f_x + i f_y (bbkick1。Yokoya の cwerf)。
GaussianCoulombU: そのポテンシャル (twspu。ロンバーグ積分、大きい引数は漸近展開)。
GaussianCoulombFitted: twspfu (21×61×61 の表を spline1 で作り 3 次のスプラインで補間。最初の呼び出しで表を作る (数秒))。gs と相対 1e-9 以内で一致 (K64 は spline1 の境界値に ktaloc の記憶を読むが 0 とみなして一致)。
"""
from __future__ import annotations

import cmath
import math
import numpy as np

_RTPI = math.sqrt(math.pi)
_EULER = 0.57721566490153286061


def cwerf(z: complex) -> complex:
    """tspac.f CWERF (Gautschi の方法。w(z) = exp(-z²) erfc(-iz))"""
    const = 1.1283791670955
    xx, yy = z.real, z.imag
    x, y = abs(xx), abs(yy)
    h2 = 0.0
    if y < 4.29 and x < 5.33:
        s = (1.0 - y / 4.29) * math.sqrt(1.0 - x * x / 28.41)
        h = 1.6 * s
        h2 = 2.0 * h
        nc = 6 + int(23.0 * s)
        nu = 9 + int(21.0 * s)
        alamda = h2 ** nc
        b = alamda == 0.0
    else:
        h = 0.0
        nc = 0
        nu = 8
        alamda = 0.0
        b = True
    r1 = r2 = s1 = s2 = 0.0
    n = nu + 1
    while True:
        n -= 1
        fn = n + 1
        t1 = y + h + fn * r1
        t2 = x - fn * r2
        c = 0.5 / (t1 * t1 + t2 * t2)
        r1 = c * t1
        r2 = c * t2
        if h > 0.0 and n <= nc:
            t1 = alamda + s1
            s1 = r1 * t1 - r2 * s2
            s2 = r2 * t1 + r1 * s2
            alamda = alamda / h2
        if n <= 0:
            break
    if b:
        rs1, rs2 = r1, r2
    else:
        rs1, rs2 = s1, s2
    rs1 = const * rs1
    if y == 0.0:
        rs1 = math.exp(-x * x)
    w = complex(rs1, const * rs2)
    if yy >= 0.0:
        return w.conjugate() if xx < 0.0 else w
    w = 2.0 * cmath.exp(-complex(x, y) ** 2) - w
    return w.conjugate() if xx > 0.0 else w


def bbkick1(dx: float, dy: float, sx: float, sy: float) -> complex:
    """tspac.f bbkick1: ガウスバンチの力 (px + i py)。sx < sy は x と y を入れ替える"""
    if sx < sy:
        return 1j * bbkick1(dy, dx, sy, sx).conjugate()
    ex = 0.5 * ((dx / sx) ** 2 + (dy / sy) ** 2)
    if ex <= 1e-12:
        t21 = -2.0 / (sx * (sx + sy))
        t43 = -2.0 / (sy * (sx + sy))
        return complex(t21 * dx, t43 * dy)
    g = math.exp(-ex) if ex <= 100.0 else 0.0
    if abs(sx - sy) / sx <= 1e-5:
        c1 = 2.0 * (1.0 - g)
        ddf = 1.0 / complex(dx, -dy)
        return complex(-c1 * ddf.real, -c1 * ddf.imag)
    r1 = sy / sx
    r2 = sx / sy
    c0 = math.sqrt(2.0 * (sx * sx - sy * sy))
    c1 = 2.0 * _RTPI / c0
    z1 = complex(abs(dx), abs(dy)) / c0
    z2 = complex(r1 * abs(dx), r2 * abs(dy)) / c0
    ddf = cwerf(z1) - g * cwerf(z2)
    re, im = ddf.real, ddf.imag
    if dx == 0.0:
        im = 0.0
    if dy == 0.0:
        re = 0.0
    px = -c1 * im
    py = -c1 * re
    if dx < 0.0:
        px = -px
    if dy < 0.0:
        py = -py
    return complex(px, py)


def _gamma0ser(x):
    if x <= 0.0:
        return float('nan') if x == 0.0 else -float('inf')
    ap = 1.0
    s = x
    d = x
    for _ in range(300):
        d = -d * x * ap
        ap += 1.0
        d = d / ap ** 2
        s += d
        if abs(d) < abs(s) * 1e-13:
            break
    return s


def _gamma0cf(x):
    fpmin = 1e-30
    a1 = -1.0
    b = x - a1
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, 301):
        an = -float(i * i)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        dl = d * c
        h *= dl
        if abs(dl - 1.0) < 1e-13:
            break
    return math.exp(-x) * h


def gamma0log(x):
    """gamma.f gamma0log (オイラー定数込み)"""
    if x < 0.0:
        return float('-inf')
    if x < 1.0:
        return _gamma0ser(x)
    return _gamma0cf(x) + math.log(x) + _EULER


def _polint(xa, ya, x):
    """bint.f polint (Neville)。xa, ya は 0 始まり、長さ n"""
    n = len(xa)
    ns = 1
    dif = abs(x - xa[0])
    for i in range(1, n + 1):
        dift = abs(x - xa[i - 1])
        if dift < dif:
            ns = i
            dif = dift
    c = list(ya)
    d = list(ya)
    y = ya[ns - 1]
    ns -= 1
    dy = 0.0
    for m in range(1, n):
        for i in range(1, n - m + 1):
            ho = xa[i - 1] - x
            hp = xa[i + m - 1] - x
            w = c[i] - d[i - 1]
            den = w / (ho - hp)
            d[i - 1] = hp * den
            c[i - 1] = ho * den
        if 2 * ns < n - m:
            dy = c[ns]
        else:
            dy = d[ns - 1]
            ns -= 1
        y += dy
    return y, dy


def rombint(f, a, b, epslon, epsabs):
    """bint.f rombint (k=5 の多項式外挿)"""
    k = 5
    itmax = 30
    ba = b - a
    s0 = ba * (f(a) + f(b)) * 0.5
    s20 = s0
    h = [0.0] * (itmax + 2)
    y = [0.0] * (itmax + 2)
    h[1] = 1.0
    y[1] = s0
    n = 1
    it = 2
    while True:
        xstep = ba / n
        x = a - xstep * 0.5
        s = 0.0
        for _ in range(n):
            x += xstep
            s += f(x)
        s *= xstep
        s1 = (s + s0) * 0.5
        h[it] = h[it - 1] * 0.25
        y[it] = s1
        if it < k:
            s2 = (2.0 * s + s0) / 3.0
            if abs(s2 - s20) < max(epslon * abs(s2), epsabs):
                return s2
            s20 = s2
        else:
            s2, ds2 = _polint(h[it - k + 1:it + 1], y[it - k + 1:it + 1], 0.0)
            if abs(ds2) < max(epslon * abs(s2), epsabs):
                return s2
            if it >= itmax:
                return s2
        s0 = s1
        n *= 2
        it += 1


def twspu(x, y, sigx, sigy, epslon, epsabs):
    """tintrb.f twspu: ガウスバンチのポテンシャル U"""
    if sigx < sigy:
        return twspu(y, x, sigy, sigx, epslon, epsabs)
    sig2 = sigx * (sigx + sigy)
    ax = x / sig2
    ay = y / sig2
    a = 0.5 * x * ax
    b = 0.5 * y * ay
    r = sigy / sigx
    if r == 1.0:
        u = 2.0 * (a + b)
        return 0.0 if u == 0.0 else gamma0log(u)
    u0 = 2.0 * (a + b / r)
    if u0 <= min(epslon, epsabs * 10.0):
        return 2.0 * (a + b / r) + math.log(0.5 + 0.5 * r)

    def f1(v):
        t = (1.0 - r) * (1.0 - v) + r
        tsq = t * t
        f = a + b / tsq
        w = v * (1.0 + t)
        xx = w * f
        if xx < 0.001:
            return f * (1.0 - xx / 2.0 * (1.0 - xx / 3.0 * (1.0 - xx / 4.0 * (1.0 - xx / 5.0 * (1.0 - xx / 6.0)))))
        return (1.0 - math.exp(-xx)) / w
    if a + b < 30.0:
        return 2.0 * rombint(f1, 0.0, 1.0, epslon, epsabs) + math.log(0.5 + 0.5 * r)
    r1 = 1.0 - r
    u = a + b
    c = r1 / u
    d = (a - b) / u
    d2 = d * d
    f0 = 1.0 - d2
    f1_ = f0 - d2
    f2 = 3.0 - 4.0 * d2
    f3 = 4.0 * f0 * f1_ - f2
    return (gamma0log(2.0 * u)
            - c * (64.0 * d + c * (-48.0 * f1_ + c * (-80.0 * d * f2 + c * (210.0 * f3
            + c * (756.0 * d * (f2 + 2.0 * f3) + c * (3465.0 * (2.0 * (d * f2) ** 2 - 1.0)
            + c * (19305.0 * d * (1.0 - 4.0 * f1_ * (f3 + f1_))))))))) / 128.0)


# ---------------------------------------------------------------------------
# GaussianCoulombFitted: tintrb.f twspfu・twspfu0・twspfuinit と spline.f spline1 (mode1, mode2 付き)
# ---------------------------------------------------------------------------

def spline1m(y, ddy, mode1, mode2):
    """spline.f spline1(np,y,ddy,work,mode1,mode2): dx=1 の自然スプライン。ddy は入出力 (mode1=1・2 は ddy[0]、mode2=1・2 は ddy[np-1] を境界の値として読む)。
    ddy は y''/6。配列は 0 始まり"""
    n = len(y)
    work = np.zeros(n)
    if mode1 == 0:
        ddy[0] = 0.0
        if n <= 2:
            ddy[1] = 0.0
            return ddy
        ddy[1] = y[2] - 2.0 * y[1] + y[0]
        work[1] = 4.0
    elif mode1 == 1:
        work[0] = 0.5
        work[1] = 3.5
        ddy[0] = 0.5 * (y[1] - y[0] - ddy[0])
        ddy[1] = y[2] - 2.0 * y[1] + y[0] - ddy[0]
    elif mode1 == 2:
        ddy[1] = y[2] - 2.0 * y[1] + y[0] - ddy[0]
        work[1] = 4.0
    for i in range(1, n - 2):                      # Fortran i=2..np-2
        work[i + 1] = 4.0 - 1.0 / work[i]
        ddy[i + 1] = y[i + 2] - 2.0 * y[i + 1] + y[i] - ddy[i] / work[i]
    if mode2 == 0:
        ddy[n - 1] = 0.0
        ddy[n - 2] = ddy[n - 2] / work[n - 2]
    elif mode2 == 1:
        ddy[n - 1] = ddy[n - 1] - (y[n - 1] - y[n - 2])
        ddy[n - 2] = (ddy[n - 2] - ddy[n - 1] * 0.5) / (work[n - 2] - 0.5)
        ddy[n - 1] = (ddy[n - 1] - ddy[n - 2]) * 0.5
    elif mode2 == 2:
        ddy[n - 2] = (ddy[n - 2] - ddy[n - 1]) / work[n - 2]
    for i in range(n - 3, 0, -1):                  # Fortran i=np-2..2
        ddy[i] = (ddy[i] - ddy[i + 1]) / work[i]
    if mode1 == 1:
        ddy[0] = ddy[0] - ddy[1] * 0.5
    return ddy


_NR, _NX, _NY = 20, 60, 60
_XM, _YM, _RM = 15.0, 30.0, 5.0
_XSTEP = _XM / _NX
_RSTEP = _RM / _NR ** 2
_YSTEP = _YM / _NY
_TABLES = {}


def twspfuinit():
    """tintrb.f twspfuinit: (nx+1, ny+1, nr+1) の表 u、uxx、uyy、uxxyy、urr、uxxrr、uyyrr、uxxyyrr。K64 は ktaloc の記憶 (0 とみなす) を spline1 の境界値に読む"""
    if _TABLES:
        return _TABLES
    nr, nx, ny = _NR, _NX, _NY
    shape = (nx + 1, ny + 1, nr + 1)
    u = np.zeros(shape)
    uxx = np.zeros(shape)
    uyy = np.zeros(shape)
    uxxyy = np.zeros(shape)
    urr = np.zeros(shape)
    uxxrr = np.zeros(shape)
    uyyrr = np.zeros(shape)
    uxxyyrr = np.zeros(shape)
    for n in range(nr + 1):
        r = 10.0 ** (-n ** 2 * _RSTEP)
        ystepr = r * _YSTEP
        for i in range(nx + 1):
            x = _XSTEP * i
            sv = np.array([twspu(x, j * ystepr, 1.0, r, 1.0e-8, 1.0e-11) for j in range(ny + 1)])
            u[i, :, n] = sv
            dds = np.zeros(ny + 1)
            dds[ny] = -bbkick1(x, ny * ystepr, 1.0, r).imag * ystepr
            spline1m(sv, dds, 1, 1)
            uyy[i, :, n] = dds
        sj = np.zeros(ny + 1)
        for j in range(ny + 1):
            ddx = uxx[:, j, n].copy()
            ddx[0] = 0.0
            ddx[nx] = 0.0                 # K64 は dds に境界値を入れるが spline1 には uxx(0,j,n) を渡す (ktaloc の記憶 = 0)
            spline1m(u[:, j, n].copy(), ddx, 1, 1)
            uxx[:, j, n] = ddx
            sj[j] = ddx[0]
        dds = np.zeros(ny + 1)
        spline1m(sj, dds, 1, 0)
        uxxyy[0, :, n] = dds
        sj2 = uxx[nx, :, n].copy()
        dds2 = np.zeros(ny + 1)
        spline1m(sj2, dds2, 0, 0)
        for j in range(ny + 1):
            uxxyy[nx, j, n] = dds2[j]
            ddv = uxxyy[:, j, n].copy()
            spline1m(uyy[:, j, n].copy(), ddv, 2, 2)
            uxxyy[:, j, n] = ddv
    for i in range(nx + 1):
        for j in range(ny + 1):
            for src, dst in ((u, urr), (uxx, uxxrr), (uyy, uyyrr), (uxxyy, uxxyyrr)):
                d = np.zeros(nr + 1)
                spline1m(src[i, j, :].copy(), d, 1, 0)
                dst[i, j, :] = d
    _TABLES.update(u=u, uxx=uxx, uyy=uyy, uxxyy=uxxyy, urr=urr, uxxrr=uxxrr, uyyrr=uyyrr, uxxyyrr=uxxyyrr)
    return _TABLES


def twspfu0(x, y, sigx, sigy):
    """tintrb.f twspfu0: (fx, fy, fu, fxx, fyy, fxy)。表の外は bbkick1 と twspu (fxx = fyy = fxy = 0)"""
    nr, nx, ny = _NR, _NX, _NY
    r = sigy / sigx
    rl = math.sqrt(-math.log10(r) / _RSTEP)
    n = int(rl)
    if n >= nr:
        return _twspfu_out(x, y, sigx, sigy)
    dxs = _XSTEP * sigx
    dys = _YSTEP * sigy
    ax = min(abs(x) / dxs, _XM * 100.0)
    i = int(ax)
    if i >= nx:
        return _twspfu_out(x, y, sigx, sigy)
    ay = min(abs(y) / dys, _YM * 100.0)
    j = int(ay)
    if j >= ny:
        return _twspfu_out(x, y, sigx, sigy)
    T = twspfuinit()
    u, uxx, uyy, uxxyy = T['u'], T['uxx'], T['uyy'], T['uxxyy']
    urr, uxxrr, uyyrr, uxxyyrr = T['urr'], T['uxxrr'], T['uyyrr'], T['uxxyyrr']
    i1 = i + 1
    j1 = j + 1
    n1 = n + 1
    bax = ax - i
    bay = ay - j
    br = rl - n
    aax = 1.0 - bax
    aay = 1.0 - bay
    ar = 1.0 - br
    ar2 = -br * (ar + 1.0)
    br2 = -ar * (br + 1.0)
    aax2 = -bax * (aax + 1.0)
    bax2 = -aax * (bax + 1.0)
    aay2 = -bay * (aay + 1.0)
    bay2 = -aay * (bay + 1.0)

    def cmb(a, arr):
        return ar * (a[arr[0], arr[1], n] + ar2 * arr_rr[arr[0], arr[1], n]) + br * (a[arr[0], arr[1], n1] + br2 * arr_rr[arr[0], arr[1], n1])

    def blend(a, rr, ii, jj):
        return ar * (a[ii, jj, n] + ar2 * rr[ii, jj, n]) + br * (a[ii, jj, n1] + br2 * rr[ii, jj, n1])

    u0 = blend(u, urr, i, j)
    u1 = blend(u, urr, i1, j)
    u2 = blend(u, urr, i, j1)
    u3 = blend(u, urr, i1, j1)
    uxx0 = blend(uxx, uxxrr, i, j)
    uxx1 = blend(uxx, uxxrr, i1, j)
    uxx2 = blend(uxx, uxxrr, i, j1)
    uxx3 = blend(uxx, uxxrr, i1, j1)
    uyy0 = blend(uyy, uyyrr, i, j)
    uyy1 = blend(uyy, uyyrr, i1, j)
    uyy2 = blend(uyy, uyyrr, i, j1)
    uyy3 = blend(uyy, uyyrr, i1, j1)
    uxxyy0 = blend(uxxyy, uxxyyrr, i, j)
    uxxyy1 = blend(uxxyy, uxxyyrr, i1, j)
    uxxyy2 = blend(uxxyy, uxxyyrr, i, j1)
    uxxyy3 = blend(uxxyy, uxxyyrr, i1, j1)
    up = aax * (u0 + aax2 * uxx0) + bax * (u1 + bax2 * uxx1)
    uq = aax * (u2 + aax2 * uxx2) + bax * (u3 + bax2 * uxx3)
    uyyp = aax * (uyy0 + aax2 * uxxyy0) + bax * (uyy1 + bax2 * uxxyy1)
    uyyq = aax * (uyy2 + aax2 * uxxyy2) + bax * (uyy3 + bax2 * uxxyy3)
    fu = aay * (up + aay2 * uyyp) + bay * (uq + bay2 * uyyq)
    aax23 = 3.0 * aax2 + 2.0
    bax23 = 3.0 * bax2 + 2.0
    aay23 = 3.0 * aay2 + 2.0
    bay23 = 3.0 * bay2 + 2.0
    fy = ((uq - up) - aay23 * uyyp + bay23 * uyyq) / dys
    fyy = -6.0 * (aay * uyyp + bay * uyyq) / dys ** 2
    ur = aay * (u0 + aay2 * uyy0) + bay * (u2 + bay2 * uyy2)
    us = aay * (u1 + aay2 * uyy1) + bay * (u3 + bay2 * uyy3)
    uxxr = aay * (uxx0 + aay2 * uxxyy0) + bay * (uxx2 + bay2 * uxxyy2)
    uxxs = aay * (uxx1 + aay2 * uxxyy1) + bay * (uxx3 + bay2 * uxxyy3)
    fx = ((us - ur) - aax23 * uxxr + bax23 * uxxs) / dxs
    fxx = -6.0 * (aax * uxxr + bax * uxxs) / dxs ** 2
    fxy = (aax23 * (uxx2 - uxx0 + bay23 * uxxyy2 - aay23 * uxxyy0)
           - bax23 * (uxx3 - uxx1 + bay23 * uxxyy3 - aay23 * uxxyy1)
           + ((u2 - u0 + bay23 * uyy2 - aay23 * uyy0)
              - (u3 - u1 + bay23 * uyy3 - aay23 * uyy1))) / dxs / dys
    if x > 0.0:
        fx = -fx
        fxy = -fxy
    if y > 0.0:
        fy = -fy
        fxy = -fxy
    return fx, fy, fu, fxx, fyy, fxy


def _twspfu_out(x, y, sigx, sigy):
    z = bbkick1(x, y, sigx, sigy)
    return z.real, z.imag, twspu(x, y, sigx, sigy, 1.0e-8, 1.0e-11), 0.0, 0.0, 0.0


def twspfu(x, y, sigx, sigy):
    """tintrb.f twspfu: (fx, fy, fu, fxx, fyy, fxy)。sigx < sigy は x と y を入れ替えて twspfu0 を呼ぶ (K64 は fxx・fyy の位置を入れ替えない)"""
    if sigx < sigy:
        a, b, c, d, e, f = twspfu0(y, x, sigy, sigx)
        return b, a, c, d, e, f
    return twspfu0(x, y, sigx, sigy)
