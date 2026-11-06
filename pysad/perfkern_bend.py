"""
perfkern_bend.py — emit6d の性能用 numba カーネル (11/05)

- apply_tdb: _State.apply (tmultr5 + tmulbs) の in-place 版 (T ← t·T, D ← t·D, B ← t·B·t^T)
- bendecorr / bendebody: K1≠0 BEND の tbendecorr / tbendebody (tbende.f 直訳。emit6d.py の旧 Python 版と同じ演算順)
  戻りは (t, 出口 cod)。呼び出し側 (emit6d._tbendecorr/_tbendebody) が apply_tdb する。
  参照値: pysad/tests/data_bend_kernels.npz (numba 化前の Python 版の出力。test_perf_kernels.py)
"""
import math

import numpy as np
from numba import njit

NJ = dict(cache=True, error_model='numpy', nogil=True)


@njit(**NJ)
def apply_tdb(t, T, D, B, want_db):
    tmp = np.empty((6, 6))
    for i in range(6):
        for j in range(6):
            s = 0.0
            for k in range(6):
                s += t[i, k] * T[k, j]
            tmp[i, j] = s
    for i in range(6):
        for j in range(6):
            T[i, j] = tmp[i, j]
    if want_db:
        for i in range(6):
            for j in range(6):
                s = 0.0
                for k in range(6):
                    s += t[i, k] * D[k, j]
                tmp[i, j] = s
        for i in range(6):
            for j in range(6):
                D[i, j] = tmp[i, j]
        for i in range(6):
            for j in range(6):
                s = 0.0
                for k in range(6):
                    s += t[i, k] * B[k, j]
                tmp[i, j] = s
        for i in range(6):
            for j in range(6):
                s = 0.0
                for k in range(6):
                    s += tmp[i, k] * t[j, k]
                B[i, j] = s


@njit(**NJ)
def _sqrt1(x):
    """mathfun sqrt1: sqrt(1+x)-1 (emit6d._sqrt1 と同じ)"""
    if abs(x) < 8e-6:
        return x * (0.5 + x * (-0.125 + x * 0.0625))
    if abs(x) < 3.4e-3:
        return x * (0.5 + x * (-0.125 + x * (0.0625 + x * (
            -.0390625 + x * (.02734375 - .0205078125 * x)))))
    return x / (1.0 + math.sqrt(max(1e-200, 1.0 + x)))


@njit(**NJ)
def bendecorr(cod, ak1, al, rho0, psqmax):
    """tbende.f tbendecorr: (t, 出口 cod)"""
    phi0 = al / rho0
    p = 1.0 + cod[5]
    xi = cod[0]
    pxi = min(p, max(-p, cod[1])) / p
    yi = cod[2]
    pyi = min(p, max(-p, cod[3])) / p
    s = min(psqmax, pxi**2 + pyi**2)
    dpzi = _sqrt1(-s)
    pzi = 1.0 + dpzi
    a = 1.0 / (pzi - phi0 * pxi)
    xf = (-al * dpzi * pxi + pzi * xi) * a
    pxf = pxi + phi0 * dpzi
    yf = yi + (-al * dpzi + xf * phi0) * pyi / pzi
    zf = cod[4] + (0.5 * s + (1.0 + xi / rho0) * dpzi / pzi) * al
    xr = xf / rho0
    pxf = pxf - ak1 * xf * xr * (0.5 - xr * (2.0 - xr) / 12.0) / p
    dk = -ak1 * xr * (1.0 - xr * (1.0 - 2.0 * xr))
    t = np.eye(6)
    t[0, 0] = pzi * a
    t[0, 1] = ((al + phi0 * xi) * (pzi**2 + pxi**2)
               - al * (phi0 * pxi**3 + pzi**3)) * a**2 / pzi / p
    t[0, 3] = (al * (1.0 - phi0 * pxi) + phi0 * xi) * pxi * pyi * a**2 / pzi / p
    t[0, 5] = pxi / pzi * a * (al * (-a * (1.0 - phi0 * pxi) + pzi**2)
                               - a * phi0 * xi)
    t[1, 0] = dk * t[0, 0]
    t[1, 1] = 1.0 - phi0 * pxi / pzi + dk * t[0, 1]
    t[1, 3] = -phi0 * pyi / pzi + dk * t[0, 3]
    t[1, 5] = -phi0 / pzi * dpzi + dk * t[0, 5]
    t[2, 0] = phi0 * pyi * a
    t[2, 1] = (-al * phi0 * (pzi**2 + pxi**2)
               + (al + phi0 * xi) * (pxi + phi0 * pzi)) * a**2 * pyi / pzi / p
    t[2, 3] = (-al * (pzi**2 * (pzi - phi0 * pxi) + phi0 * pxi * pyi**2)
               + (pyi**2 + pzi * (pzi - phi0 * pxi)) * (al + phi0 * xi)
               ) * a**2 / p / pzi
    t[2, 5] = -(al * (dpzi * (1.0 + pzi) * phi0 * pxi - pzi**3)
                + (al + phi0 * xi)) * a**2 * pyi / pzi / p
    t[4, 0] = -t[0, 0] * t[1, 5]
    t[4, 1] = -t[0, 1] * t[1, 5] + t[1, 1] * t[0, 5]
    t[4, 3] = t[2, 5]
    t[4, 5] = (-s / p + (1.0 + xi / rho0) * s / pzi**3) * al
    c = cod.copy()
    c[0] = xf
    c[1] = min(p, max(-p, pxf * p))
    c[2] = yf
    c[4] = zf
    return t, c


@njit(**NJ)
def bendebody(cod, al, phin, ak1, rhob, rho0, h0, h1emit, dvemit, psqmax):
    """tbende.f tbendebody (tbinit=.true. で毎回 tbendiinit/tbendeinit): (t, 出口 cod)"""
    dp = cod[5]
    p = 1.0 + dp
    xi = cod[0]
    pxi = min(p, max(-p, cod[1])) / p
    yi = cod[2]
    pyi = min(p, max(-p, cod[3])) / p
    zi = cod[4]
    # tbendiinit
    rhoe = rhob * p
    rhosq = rho0 * rhoe
    drhop = (rhoe - rho0) / rhosq
    akk = ak1 / al
    akxsq = -1.0 / rhosq - akk / p
    akysq = akk / p
    if akxsq > 0.0:
        akx = math.sqrt(akxsq)
        phix = akx * al
        dcx = 2.0 * math.sinh(0.5 * phix)**2
        spx = math.sinh(phix)
        xspx = phix - spx
        aksx = akx * spx
        dcxkx = dcx / akxsq
        sxkx = spx / akx
        xsxkx = xspx / akx / akxsq
    elif akxsq < 0.0:
        akx = math.sqrt(-akxsq)
        phix = akx * al
        dcx = -2.0 * math.sin(0.5 * phix)**2
        spx = math.sin(phix)
        xspx = phix - spx
        aksx = -akx * spx
        dcxkx = dcx / akxsq
        sxkx = spx / akx
        xsxkx = -xspx / akx / akxsq
    else:
        akx = 0.0
        phix = 0.0
        dcx = 0.0
        spx = 0.0
        aksx = 0.0
        dcxkx = 0.5 * al**2
        sxkx = al
        xsxkx = -1.0 / 6.0 * al**3
    if akysq > 0.0:
        aky = math.sqrt(akysq)
        phiy = aky * al
        dcy = 2.0 * math.sinh(0.5 * phiy)**2
        spy = math.sinh(phiy)
        aksy = aky * spy
        syky = spy / aky
    elif akysq < 0.0:
        aky = math.sqrt(-akysq)
        phiy = aky * al
        dcy = -2.0 * math.sin(0.5 * phiy)**2
        spy = math.sin(phiy)
        aksy = -aky * spy
        syky = spy / aky
    else:
        aky = 0.0
        phiy = 0.0
        dcy = 0.0
        spy = 0.0
        aksy = 0.0
        syky = al
    # tbendeinit
    cx = 1.0 + dcx
    cy = 1.0 + dcy
    drhopp = 1.0 / rhoe / p
    akxsqp = -akxsq / p
    phixp = -0.5 * phix / p
    if akxsq > 0.0:
        dcxp = spx * phixp
        dcxkxp = dcxp / akxsq + dcxkx / p
        sxkxp = 0.5 * sxkx / p + cx * phixp / akx
        xsxkxp = (-dcx * phixp / akxsq - 1.5 * xsxkx / p) / akx
        aksxp = -0.5 * aksx / p + akx * cx * phixp
    elif akxsq < 0.0:
        dcxp = -spx * phixp
        dcxkxp = dcxp / akxsq + dcxkx / p
        sxkxp = 0.5 * sxkx / p + cx * phixp / akx
        xsxkxp = (-dcx * phixp / akxsq - 1.5 * xsxkx / p) / akx
        aksxp = -0.5 * aksx / p - akx * cx * phixp
    else:
        dcxp = 0.0
        dcxkxp = 0.0
        sxkxp = 0.0
        xsxkxp = 0.0
        aksxp = 0.0
    akysqp = -akysq / p
    phiyp = -0.5 * phiy / p
    if akysq > 0.0:
        dcyp = phiyp * spy
        sykyp = 0.5 * syky / p + cy * phiyp / aky
        aksyp = -0.5 * aksy / p + aky * cy * phiyp
    elif akysq < 0.0:
        dcyp = -phiyp * spy
        sykyp = 0.5 * syky / p + cy * phiyp / aky
        aksyp = -0.5 * aksy / p - aky * cy * phiyp
    else:
        dcyp = 0.0
        sykyp = 0.0
        aksyp = 0.0
    # tbendebody
    xiksq = drhop + akxsq * xi
    dxf = dcxkx * xiksq + sxkx * pxi
    dpxf = sxkx * xiksq + dcx * pxi
    dyf = dcy * yi + syky * pyi
    dpyf = aksy * yi + dcy * pyi
    hi = (0.5 * (pxi**2 + pyi**2 - akxsq * xi**2 - akysq * yi**2)
          - drhop * xi)
    hip = (0.5 * (-2.0 * (pxi**2 + pyi**2) / p - akxsqp * xi**2
                  - akysqp * yi**2) - drhopp * xi)
    dxpx = dxf * pxi + xi * dpxf + dxf * dpxf + dyf * pyi + yi * dpyf + dyf * dpyf
    zf = (zi - 0.25 * (dxpx - (5.0 / rho0 - 1.0 / rhoe)
                       * (drhop * xsxkx - sxkx * xi - dcxkx * pxi))
          - hi * al * 0.5)
    xf = xi + dxf
    pxf = pxi + dpxf
    yf = yi + dyf
    pyf = pyi + dpyf
    t = np.eye(6)
    t[0, 0] = cx
    t[0, 1] = sxkx / p
    t[0, 5] = (drhopp * dcxkx + drhop * dcxkxp
               + xi * dcxp + (sxkxp - sxkx / p) * pxi)
    t[1, 0] = aksx * p
    t[1, 1] = cx
    t[1, 5] = (p * (drhopp * sxkx + drhop * sxkxp) + drhop * sxkx
               + (p * aksxp + aksx) * xi + dcxp * pxi)
    t[2, 2] = cy
    t[2, 3] = syky / p
    t[2, 5] = yi * dcyp + (sykyp - syky / p) * pyi
    t[3, 2] = aksy * p
    t[3, 3] = cy
    t[3, 5] = (aksyp * p + aksy) * yi + pyi * dcyp
    t[4, 0] = -t[0, 0] * t[1, 5] + t[1, 0] * t[0, 5]
    t[4, 1] = -t[0, 1] * t[1, 5] + t[1, 1] * t[0, 5]
    t[4, 2] = -t[2, 2] * t[3, 5] + t[3, 2] * t[2, 5]
    t[4, 3] = -t[2, 3] * t[3, 5] + t[3, 3] * t[2, 5]
    t[4, 5] = (-0.25 * (-dxpx / p
                        + t[0, 5] * pxf + t[1, 5] / p * xf
                        + t[2, 5] * pyf + t[3, 5] / p * yf
                        - 1.0 / rhoe / p * (drhop * xsxkx - sxkx * xi
                                            - dcxkx * pxi)
                        - (5.0 / rho0 - 1.0 / rhoe)
                        * (drhopp * xsxkx + drhop * xsxkxp
                           - sxkxp * xi - dcxkxp * pxi
                           + dcxkx * pxi / p))
               - hip * al * 0.5 + h0 / h1emit**3 * al)
    c = cod.copy()
    c[0] = xf
    c[1] = min(p, max(-p, pxf * p))
    c[2] = yf
    c[3] = min(p, max(-p, pyf * p))
    c[4] = zf - dvemit * al
    return t, c
