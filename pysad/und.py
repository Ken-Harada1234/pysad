"""
UND 要素 (アンジュレータ、K64 の undulator.f: K.Ohmi と E.Forest、KEK Report 92-14)。追跡だけ。

K64 では UND に磁場のキーワードが無い: initb1.f は `kytbl(kwFBX,icUND)=ky_FBX_UND` と書くが kwFBX は未定義の変数 (=0) で、
`FBX` `FBY` `BX` `BY` はどれも "Undefined keyword" になる。そのため F0 = G0 = 0 で、undulator は磁場なしの (正準運動量の
近似の) ドリフトを刻みごとに行う。Twiss・Emittance[] の写像は無い (恒等、長さは幾何にだけ入る)。
直訳は F0, G0 を保ったままにしてあるが、利用者は設定できない。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

NJ = dict(cache=True, error_model='numpy')
_BAD = ('FBX', 'FBY', 'BX', 'BY')


def compile_und(el, pgev):
    """undinit。戻り値は k_und の引数 (nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds)"""
    p = el.params
    g = lambda k: float(p.get(k, 0.0) or 0.0)
    for k in _BAD:
        if k in p:
            raise NotImplementedError(f'UND: キーワード {k} は K64 に無い (Undefined keyword)')
    length = float(el.L)
    nsli = g('SLICE')
    npole = g('POLE')
    f0 = 0.0 * pgev            # p_in(2)/pgev*c: FBX が設定できないので 0
    g0 = 0.0
    dph = g('PHI') * math.pi / 180.0
    with np.errstate(all='ignore'):
        kz = 2.0 * math.pi * npole / length if length != 0.0 else float('nan')
        kx = g('KX')
        ky = math.sqrt(kz * kz + kx * kx)
        qy = g('QY')
        qx = math.sqrt(kz * kz + qy * qy)
        ds = length / nsli if nsli != 0.0 else float('inf')
    return (int(nsli + 0.01), f0, g0, dph, kz, kx, ky, qx, qy, ds)


@njit(**NJ)
def k_und(S, n, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, p0, h0):
    for i in range(n):
        pn = 1.0 + S[i, 5]
        S[i, 1] = S[i, 1] * pn
        S[i, 3] = S[i, 3] * pn
    for i in range(n):
        x = S[i, 0]
        px = S[i, 1]
        y = S[i, 2]
        py = S[i, 3]
        z = S[i, 4]
        g = S[i, 5]
        for j in range(1, nslice + 1):
            kx1 = Kx * x
            kx2 = Ky * y
            kx3 = Qx * x
            kx4 = Qy * y
            sinkz = math.sin(Kz * j * ds)
            sinkzp = math.sin(Kz * j * ds + dph)
            ff = (F0 / Kz * math.cos(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Qx / Kz * math.sinh(kx3) * math.sin(kx4) * sinkzp)
            gg = (-G0 / Kz * math.cos(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Ky / Kz * math.sinh(kx2) * math.sin(kx1) * sinkz)
            fx = (-Kx * F0 / Kz * math.sin(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Kz * math.cosh(kx3) * math.sin(kx4) * sinkzp)
            fy = (Ky * F0 / Kz * math.cos(kx1) * math.sinh(kx2) * sinkz
                  - Qy * Qy * G0 / Qx / Kz * math.sinh(kx3) * math.cos(kx4) * sinkzp)
            gx = (-G0 * Qx / Kz * math.cos(kx4) * math.sinh(kx3) * sinkzp
                  + Kx * Kx * F0 / Ky / Kz * math.sinh(kx2) * math.cos(kx1) * sinkz)
            gy = (G0 * Qy / Kz * math.sin(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Kz * math.cosh(kx2) * math.sin(kx1) * sinkz)
            dpds = ds / (1.0 + g)
            fx = dpds * fx
            gx = dpds * gx
            fy = dpds * fy
            gy = dpds * gy
            det1a = 1.0 - fx - gy + fx * gy - fy * gx
            fgx = dpds * (ff * fx + gg * gx)
            fgy = dpds * (ff * fy + gg * gy)
            tmp = fx
            fx = (1.0 - gy) / det1a
            gx = gx / det1a
            fy = fy / det1a
            gy = (1.0 - tmp) / det1a
            tmp = px - fgx
            px = fx * tmp + gx * (py - fgy)
            py = fy * tmp + gy * (py - fgy)
            fgx = px - ff
            fgy = py - gg
            x = x + dpds * fgx
            y = y + dpds * fgy
            z = z - 0.5 * (fgx * fgx + fgy * fgy) * (dpds * dpds) / ds
        S[i, 0] = x
        S[i, 1] = px
        S[i, 2] = y
        S[i, 3] = py
        S[i, 4] = z
    for i in range(n):
        pn = 1.0 + S[i, 5]
        S[i, 1] = S[i, 1] / pn
        S[i, 3] = S[i, 3] / pn
        h1 = math.sqrt(1.0 + (pn * p0)**2)
        S[i, 6] = -S[i, 5] * (1.0 + pn) / h1 / (h1 + pn * h0) + S[i, 10]
