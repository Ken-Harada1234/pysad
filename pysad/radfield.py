"""
RadiationField / RadiationSpectrum。K64 の tlinit.f (tlfield、tlcalc、tlspect) の直訳。
RadiationField[traj, obs]: traj = {{t1..tm},{x1..xm},{y1..ym},{z1..zm}} (TrackParticles の RADLIGHT の軌道。GEO 座標)、obs = 観測点 {x,y,z}。
  Feynman-Heaviside の式 (E = (mu0 e CHARGE/4π)(c²n/R² + R/c d(c²n/R²)/dt + d²n/dt²)) の時間領域の場。微分は自然スプライン (K64 の spline: y''/6)。
  戻り: {{τ..},{Ex..},{Ey..},{Ez..},{Hx..},{Hy..},{Hz..},{Sx..},{Sy..},{Sz..}}。
RadiationSpectrum[{τ, f}, {λ1, λ2, dλ}]: 場 f(τ) の cos・sin 積分の和 (区分線形近似の解析解)。戻り: {{k..},{c..},{s..}}。
RADLIGHT での軌道の記録は ffscmd の `_RadLight` (TrackParticles が {結果, 軌道} を返す。gs と一致: test_radlight.py)。未対応: RadiationSpectrum の NPARA による fork 並列 (逐次で計算)。
"""
from __future__ import annotations

import math

import numpy as np

from .sadscript import SADError, is_num

_C = 299792458.0
_ELEMCH = 1.602176634e-19


def _natural_ddy(x, y):
    """spline.f mode 0: 自然スプライン、戻りは y''(i)/6"""
    n = len(x)
    if n <= 2:
        return np.zeros(n)
    from scipy.interpolate import CubicSpline
    return CubicSpline(x, y, bc_type='natural')(x, 2) / 6.0


def radiation_field(traj, obs, charge=-1.0, e=_ELEMCH):
    """tlcalc。traj: (m, 4) = [t, x, y, z]。戻り: (m, 10)"""
    xt = np.array(traj, dtype=float)
    m = xt.shape[0]
    tx, ty, tz = [float(v) for v in obs]
    c = _C
    z0 = 4e-7 * 2.0 * math.asin(1.0) * c
    r0 = math.sqrt(tx * tx + ty * ty + tz * tz)
    f = np.zeros((m, 10))
    for i in range(m):
        dxi = tx - xt[i, 1]
        dyi = ty - xt[i, 2]
        dzi = tz - xt[i, 3]
        ri = math.sqrt(dxi * dxi + dyi * dyi + dzi * dzi)
        dri = -(xt[i, 1] * (tx + dxi) + xt[i, 2] * (ty + dyi) + xt[i, 3] * (tz + dzi)) / (r0 + ri)
        xt[i, 0] = xt[i, 0] + dri / c
        xt[i, 1] = dxi / ri
        xt[i, 2] = dyi / ri
        xt[i, 3] = dzi / ri
        f[i, 7] = ri / c
    t = xt[:, 0]
    for i in (1, 2, 3):
        f[:, i] = _natural_ddy(t, xt[:, i])
    for i in (4, 5, 6):
        f[:, 9] = xt[:, i - 3] / f[:, 7] ** 2
        f[:, i] = _natural_ddy(t, f[:, 9])
        f1 = f[m - 2, i]
        for j in range(m - 1):
            dt = t[j + 1] - t[j]
            f[j, i] = f[j, 9] + f[j, 7] * ((f[j + 1, 9] - f[j, 9]) / dt - dt * (2.0 * f[j, i] + f[j + 1, i]))
        dt = t[m - 1] - t[m - 2]
        f[m - 1, i] = f[m - 1, 9] + f[m - 1, 7] * ((f[m - 1, 9] - f[m - 2, 9]) / dt + dt * (f1 + 2.0 * f[m - 1, i]))
    coeff = 1e-7 * e * charge
    out = np.zeros((m, 10))
    for i in range(m):
        out[i, 0] = xt[i, 0]
        out[i, 1] = coeff * (6.0 * f[i, 1] + f[i, 4])
        out[i, 2] = coeff * (6.0 * f[i, 2] + f[i, 5])
        out[i, 3] = coeff * (6.0 * f[i, 3] + f[i, 6])
        out[i, 4] = (xt[i, 2] * out[i, 3] - xt[i, 3] * out[i, 2]) / z0
        out[i, 5] = (xt[i, 3] * out[i, 1] - xt[i, 1] * out[i, 3]) / z0
        out[i, 6] = (xt[i, 1] * out[i, 2] - xt[i, 2] * out[i, 1]) / z0
        out[i, 7] = out[i, 2] * out[i, 6] - out[i, 3] * out[i, 5]
        out[i, 8] = out[i, 3] * out[i, 4] - out[i, 1] * out[i, 6]
        out[i, 9] = out[i, 1] * out[i, 5] - out[i, 2] * out[i, 4]
    return out


def radiation_spectrum(tau, f, rng):
    """tlspect (逐次)。rng = {λ1, λ2, dλ}。戻り: (k, c, s)"""
    t = np.asarray(tau, dtype=float)
    y = np.asarray(f, dtype=float)
    m = len(t)
    al0, al1, dal = [float(v) for v in rng]
    na = int((al1 - al0) / dal + 2.5)
    ak0 = 2.0 * math.pi / al1
    ak1 = 2.0 * math.pi / al0
    dak = (ak1 - ak0) / (na - 1)
    slope = (y[1:] - y[:-1]) / (t[1:] - t[:-1])
    da = slope[:-1] - slope[1:]           # a0 - a1 (i = 2..m-1)
    ks, cs, ss = [], [], []
    for k in range(na):
        ak = ak0 + dak * k
        w = ak * _C
        x = w * t[1:m - 1]
        c = float(np.sum(np.cos(x) * da))
        s = float(np.sum(np.sin(x) * da))
        a0 = slope[m - 2]
        xm = w * t[m - 1]
        c = c / w + math.cos(xm) * a0 / w + y[m - 1] * math.sin(xm)
        s = s / w + math.sin(xm) * a0 / w - y[m - 1] * math.cos(xm)
        a1 = slope[0]
        x1 = w * t[0]
        c = (c - math.cos(x1) * a1 / w - y[0] * math.sin(x1)) / w
        s = (s - math.sin(x1) * a1 / w + y[0] * math.cos(x1)) / w
        ks.append(ak)
        cs.append(c)
        ss.append(s)
    return ks, cs, ss


def install(I):
    B = I.builtins

    def _rf(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in RadiationField')
        tr, ob = a
        if not isinstance(tr, list) or len(tr) != 4 or not all(isinstance(r, list) for r in tr):
            raise SADError('General::wrongleng: "#1","4"')
        if not isinstance(ob, list) or len(ob) != 3:
            raise SADError('General::wrongleng: "#2","3"')
        ch = I.globals.get('CHARGE', -1.0)
        ch = float(ch) if is_num(ch) else -1.0
        out = radiation_field(np.array(tr, dtype=float).T, ob, charge=ch)
        return [[float(v) for v in out[:, j]] for j in range(10)]
    B['RadiationField'] = _rf

    def _rs(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in RadiationSpectrum')
        p, r = a
        if not isinstance(p, list) or len(p) != 2:
            raise SADError('General::wrongleng: "#1","2"')
        if len(p[0]) != len(p[1]):
            raise SADError('General::equalleng: "#1[[1]] and #[[2]]"')
        if not isinstance(r, list) or len(r) != 3:
            raise SADError('General::wrongleng: "#","3"')
        return [list(x) for x in radiation_spectrum(p[0], p[1], r)]
    B['RadiationSpectrum'] = _rs
