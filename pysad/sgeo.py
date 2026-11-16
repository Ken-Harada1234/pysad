"""
SOL 領域 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで) の幾何。K64 の tsgeo.f (tsgeo、tsgeo1)。

GEO≠0 の SOL (領域の始端か終端のどちらか一方) を基準に、SOL のずれ (DX DY) と向き (DPX DPY = CHI1 CHI2) で与えた軌道を
領域の中で追い (ドリフト・四極・領域内の SOL の写像)、その軌道に沿って幾何 geo(k+1..ke+1) と pos を作る。
GEO が始端なら前向き (dir=+1)、終端なら後ろ向きに追って全体を始端の幾何につなぎ直す。

対応: DRIFT、QUAD (K1≠0、ずれ・F1・FRINGE なし、tracking と同じ制限)、L=0 の要素、領域内の SOL (BOUND=0)。
未対応 (NotImplementedError): 領域内の BEND・MULT・CAVI。SOL に DX DY DPX DPY がない (軌道が軸上) ときは直進と同じ結果になるので、その場合は呼び出し側 (geo.py) が従来の直進で済ませる。
"""
from __future__ import annotations

import math

import numpy as np

from . import track6n as _n
from . import track6n_ext as _x
from . import geo as _g

_RS = np.array([1], dtype=np.int64)


def _val(el, key, default=0.0):
    v = getattr(el, key, None)
    if v is None:
        return default
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


class SolRegion:
    """tsgeo の結果: ke (領域の終点、1 始まり)、SOL ごとの CHI1..3・DX DY DZ (K64 が要素に書く値。追跡の tsconv が使う)"""
    def __init__(self, ke):
        self.ke = ke
        self.sol_values = {}


def _tsfrin(x, px, y, py, z, dbz, ph):
    """tsfrin (g=0): (x, px, y, py, z) → 同じ。z は ds の累積用に呼び出し側が 0 から始める"""
    S = _n.make_state(np.array([[x, px, y, py, z, 0.0]]), ph)
    _x.tsfrin1(0, S, dbz)
    return float(S[0, 0]), float(S[0, 1]), float(S[0, 2]), float(S[0, 3]), float(S[0, 4])


def _quad_step(el, bz, xi, pxi, yi, pyi, dir_, f, ph):
    """tsgeo1 の QUAD: tquad (bz*dir) を 1 粒子に使う"""
    k, pr = _x.compile_sol_quad(el, bz)
    L, ak, ct, st, fringe, _bz, radok = pr
    S = _n.make_state(np.array([[xi, (pxi - f * yi) * dir_, yi, (pyi + f * xi) * dir_, 0.0, 0.0]]), ph)
    _x.k_sol_quad(S, 1, L, ak, ct, st, fringe, bz * dir_, False, False, _RS, ph.p0, ph.h0, ph.anrad, ph.cuc)
    xf, pxf, yf, pyf, zf = (float(S[0, 0]), float(S[0, 1]), float(S[0, 2]), float(S[0, 3]), float(S[0, 4]))
    return xf, pxf * dir_ + f * yf, yf, pyf * dir_ - f * xf, -zf


def _bz_after(elements, k, ke, ph):
    """領域内の各要素の直後の bz (K64 tfbzs(l) = l 番目の要素を通ったあとの領域の bz。0 始まりの添字 l: k..ke)"""
    out = {}
    b = 0.0
    for l in range(k, ke + 1):
        el = elements[l]
        if type(el).__name__ == 'Solenoid':
            b = _x._t6._tfbzs(el, ph) if l < ke else 0.0
        out[l] = b
    return out


def tsgeo(elements, k1, geo, pos, ph, warn):
    """k1: SOL の位置 (1 始まり)。戻り: (SolRegion, ke1) か None (K64 の sol=.false.: 領域として扱わない)。geo (nlat,3,4)、pos (nlat) を直接更新する。
    geo[k-1] = K64 の geo(:,:,k)"""
    nlat = len(elements) + 1
    k = k1 - 1                                           # 0 始まりの要素の添字 (K64 の k に対応するのは k+1)
    cmp1 = elements[k]
    kg = 0
    ke = None
    if _val(cmp1, 'BOUND') != 0.0:
        if _val(cmp1, 'GEO') != 0.0:
            kg = k1
        for i in range(k1 + 1, nlat):
            e = elements[i - 1]
            if type(e).__name__ == 'Solenoid' and _val(e, 'BOUND') != 0.0:
                if _val(e, 'GEO') != 0.0:
                    if kg != 0:
                        warn(' Duplicated GEO of Solenoid ' + cmp1.name)
                        return None
                    kg = i
                ke = i
                break
    if ke is None:
        warn(' Missing BOUND of Solenoid ' + cmp1.name)
        return None
    if kg == 0:
        warn(' Missing GEO of Solenoid ' + cmp1.name)
        return None
    cmp2 = elements[ke - 1]
    ke1 = ke + 1
    if kg == k1:
        kk1, kk2, cmpp, cmpe, dir_, idir = k1, ke, cmp1, cmp2, 1.0, 1
    else:
        kk1, kk2, cmpp, cmpe, dir_, idir = ke, k1, cmp2, cmp1, -1.0, -1
    chi1 = _val(cmpp, 'DPX')
    chi2 = _val(cmpp, 'DPY')
    xi = _val(cmpp, 'DX')
    yi = _val(cmpp, 'DY')
    pxi = -math.sin(chi1) * math.cos(chi2)
    pyi = -math.sin(chi2)
    bza = _bz_after(elements, k1 - 1, ke - 1, ph)         # 0 始まり
    reg = SolRegion(ke)
    gk = geo[k1 - 1]
    if idir > 0:
        pos0 = 0.0
        bzs = bza[kk1 - 1]
        chi3 = _g.tfchi(gk, 3)
        geo[k1][:, 0:3] = _g.tfderotgeo(gk, (chi1, chi2, chi3))
        geo[k1][:, 3] = gk[:, 3] + geo[k1][:, 0:3] @ np.array([-xi, -yi, 0.0])
        reg.sol_values[k1] = tuple(_g.tgrot(gk, geo[k1]))
        pos[k1] = pos[k1 - 1]
    else:
        geos = gk.copy()
        pos0 = pos[k1 - 1]
        bzs = bza[kk1 - 2]
        geo[ke - 1] = _g.GEOINI
        geo[ke1 - 1][:, 0:3] = _g.tfchitogeo((chi1, chi2, 0.0))
        reg.sol_values[ke] = tuple(_g.tgrot(geo[ke - 1][:, 0:3], geo[ke1 - 1]))
        pos[ke - 1] = 0.0
        pos[ke1 - 1] = 0.0
        geo[ke1 - 1][:, 3] = 0.0
        geo[ke - 1][:, 3] = geo[ke - 1][:, 0:3] @ np.array([-xi, -yi, 0.0])
    ds = 0.0
    if _val(cmpp, 'DISFRIN') == 0.0:
        xi, pxi, yi, pyi, ds = _tsfrin(xi, pxi, yi, pyi, ds, bzs, ph)
        ds = -ds * dir_
    pxi = pxi + yi * bzs * 0.5
    pyi = pyi - xi * bzs * 0.5
    state = [xi, pxi, yi, pyi, ds, bzs]
    for i in range(kk1 + idir, kk2 + idir, idir):
        _tsgeo1(elements, i, state, geo, pos, idir, dir_, bza, ph)
    xi, pxi, yi, pyi, ds, bzs = state
    reg.sol_values[('end', ke if idir > 0 else k1)] = (xi, yi, ds)
    if idir > 0:
        chi3 = _g.tfchi(geo[ke1 - 1][:, 0:3], 3)
        geo[ke1 - 1][:, 0:3] = _g.tfderotgeo(geo[ke1 - 1][:, 0:3], (0.0, 0.0, chi3))
        reg.sol_values[ke] = tuple(_g.tgrot(geo[ke - 1][:, 0:3], geo[ke1 - 1]))
    else:
        a = _g.tgrot(geo[k1 - 1], geo[ke1 - 1])
        geot = _g.tfderotgeo(geos, (0.0, 0.0, a[2])) @ geo[k1 - 1][:, 0:3].T
        dpos = pos0 - pos[k1 - 1]
        pos[k1 - 1:ke1] += dpos
        dg4 = geos[:, 3] - geot @ geo[k1 - 1][:, 3]
        for i in range(k1 + 1, ke1 + 1):
            geo[i - 1] = geot @ geo[i - 1]
            geo[i - 1][:, 3] += dg4
        geo[k1 - 1] = geos
        reg.sol_values[k1] = tuple(_g.tgrot(geos, geo[k1]))
    return reg, ke1


def _tsgeo1(elements, i, state, geo, pos, idir, dir_, bza, ph):
    """tsgeo1 (i: 1 始まり)"""
    xi, pxi, yi, pyi, ds, bzs = state
    i0 = i + (1 - idir) // 2
    i1 = 2 * i + 1 - i0
    el = elements[i - 1]
    tn = type(el).__name__
    al = float(el.L)
    dl = dx = dy = 0.0
    if tn == 'Drift' and getattr(el, '_sad_type', 'DRIFT') == 'DRIFT':
        pzi = 1.0 - (pxi ** 2 + pyi ** 2) / (1.0 + math.sqrt((1.0 - pyi) * (1.0 + pyi) - pxi ** 2))
        phi = al * bzs / pzi * dir_
        ak = bzs
        sinphi = math.sin(phi)
        a22 = math.cos(phi)
        if ak == 0.0:
            a12 = al / pzi * dir_
            a14 = 0.0
        else:
            a12 = sinphi / ak
            a14 = a12 * sinphi / (1.0 + a22) if a22 >= 0.0 else (1.0 - a22) / ak
        a24 = sinphi
        dxx = a12 * pxi + a14 * pyi
        pxf = a22 * pxi + a24 * pyi
        dyy = -a14 * pxi + a12 * pyi
        pyf = -a24 * pxi + a22 * pyi
        dl = (pxi ** 2 + pyi ** 2) / pzi / (1.0 + pzi) * al
        xi = xi + dxx
        yi = yi + dyy
    elif tn == 'Quadrupole':
        f = 0.5 * bzs
        if al != 0.0 and float(getattr(el, 'K1', 0.0) or 0.0) != 0.0:
            xf, pxf, yf, pyf, dl = _quad_step(el, bzs, xi, pxi, yi, pyi, dir_, f, ph)
            dx, dy = xf - xi, yf - yi
            xi, yi = xf, yf
        else:
            raise NotImplementedError('tsgeo: 領域内の L=0 / K1=0 の四極')
    elif tn == 'Solenoid':
        if _val(el, 'BOUND') == 0.0:
            bzs0 = bzs
            bzs = bza[(i - 2) if idir < 0 else (i - 1)]
            db = bzs - bzs0
            if _val(el, 'DISFRIN') == 0.0:
                pxf0 = pxi - yi * bzs0 * 0.5
                pyf0 = pyi + xi * bzs0 * 0.5
                xf, pxf0, yf, pyf0, dl = _tsfrin(xi, pxf0, yi, pyf0, 0.0, db, ph)
                xi, yi = xf, yf
                dl = -dl
                pxf = pxf0 + yi * bzs * 0.5
                pyf = pyf0 - xi * bzs * 0.5
            else:
                pxf = pxi + yi * db * 0.5
                pyf = pyi - xi * db * 0.5
        else:
            pxf = pxi - yi * bzs * 0.5
            pyf = pyi + xi * bzs * 0.5
            if _val(el, 'DISFRIN') == 0.0:
                xf, pxf, yf, pyf, dl = _tsfrin(xi, pxf, yi, pyf, 0.0, -bzs, ph)
                xi, yi = xf, yf
                dl = -dl
            geo[i1 - 1] = _g.tforbitgeo(geo[i0 - 1], (xi, pxf, yi, pyf))
            pxi, pyi = pxf, pyf
            pos[i1 - 1] = pos[i0 - 1] + dl * dir_
            ds = ds + dl
            state[:] = [xi, pxi, yi, pyi, ds, bzs]
            return
        al = 0.0
    elif tn in ('Bend', 'Multipole', 'Cavity') or (tn == 'Drift' and getattr(el, '_sad_type', 'DRIFT') != 'DRIFT' and al != 0.0):
        raise NotImplementedError(f'tsgeo: 領域内の {tn} は未対応 (SOL にずれ・傾きがあるときの幾何)')
    else:
        al = 0.0
        pxf, pyf = pxi, pyi
    geo[i1 - 1][:, 3] = geo[i0 - 1][:, 3] + geo[i0 - 1][:, 2] * al * dir_
    geo[i1 - 1][:, 0:3] = geo[i0 - 1][:, 0:3]
    pxi, pyi = pxf, pyf
    pos[i1 - 1] = pos[i0 - 1] + (al + dl) * dir_
    ds = ds + dl
    state[:] = [xi, pxi, yi, pyi, ds, bzs]
