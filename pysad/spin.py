"""
スピン (POL フラグ) の直訳: photons.f の tradke のスピン部、tspin の spnorm / srequpol / spdepol、tchge.f の srot 回転。

状態 (temw / temit の srot(3,9), gintd(3)):
    srot[:, 0:3]  スピン歳差の 1 周行列 (srotinit で単位行列)
    srot[:, 3:9]  (x,px,y,py,z,pz) の初期座標に対する 1 周のスピン微分 (脱偏極ベクトル)
    gintd[3]      sflc·anp·(bt·h1·p/al1)² のスピン軸成分の積算 (Sokolov-Ternov の傾き)。1/sdamp が偏極時間
K64 の癖 (docs/21 B55〜B57) と v2 (compat.v2()) の扱い (docs/history/v2_changes_A.md):
  - tadddsg の呼び出しは dgz に dgx を渡す (dsg は spequ でしか使わず、spequ は呼ばれないので無影響。ここでは dsg を持たない)
  - tbrote の srot 更新は列ごとに行を逐次更新する (rr·srot の行列積ではない: 2 行目が更新後の 1 行目を使う)。v2 は行列積 (B56)
  - serot の c2 = (−s·c1' + c2)/c と、tradke の末尾の srot(3,:) = (−sp·srot(1,:)' + srot(3,:))/cp は、更新後の値 (') を使うが、
    展開すると −s·c1 + c·c2 で、正しい回転をその場で計算する書き方 (誤りではない。v2 でも同じ)
  - 平衡偏極度: K64 の spdepol は NaN (cerfc の桁あふれ)。v2 は桁あふれを直したうえで、線形のスピン軌道結合の脱偏極率
    (slim_depolarization) を使う (B55)
"""
from __future__ import annotations

import math

import numpy as np

GSPIN = 0.001159652181280002
CL = 1.0 + GSPIN
PST = 8.0 * math.sqrt(3.0) / 15.0
CFPD = 5.0 * math.sqrt(3.0) / 8.0
CEZLIM = 18.0
M_2PI = 2.0 * math.pi
M_SQRTPI = math.sqrt(math.pi)


def sflc():
    from .emit6d import ELRADI, FINEST
    return 0.75 * (ELRADI / FINEST) ** 2


def srotinit():
    s = np.zeros((3, 9))
    s[0, 0] = s[1, 1] = s[2, 2] = 1.0
    return s


def xsincos(x):
    """tfloor.f xsincos: sin, cos と 1−cos の精度の良い形 (s, c, dc)"""
    s = math.sin(x)
    c = math.cos(x)
    dc = -2.0 * math.sin(0.5 * x) ** 2
    return s, c, dc


def _pxy2dpz(px, py):
    a = px * px + py * py
    return -a / (1.0 + math.sqrt(max(0.0, 1.0 - a)))


# ---------------------------------------------------------------------------
# tradke のスピン部 (photons.f 1521-1616)
# ---------------------------------------------------------------------------

def tradke_spin(st, c):
    """c: _State.tradke のローカル量の dict。st.srot / st.gintd を更新する。transr は更新前の値 (st.transr)"""
    pr = c['pr']
    px, py, pz = c['px'], c['py'], c['pz']
    pxi, pyi, pz00, pz0, pxr0 = c['pxi'], c['pyi'], c['pz00'], c['pz0'], c['pxr0']
    cp, sp = c['cp'], c['sp']
    xpx, xpy, xpz, xpa = c['xpx'], c['xpy'], c['xpz'], c['xpa']
    bzh, bzhr0, al = c['bzh'], c['bzhr0'], c['al']
    h1, p, p0, al1, anp = c['h1'], c['p'], c['p0'], c['al1'], c['anp']
    from .emit6d import _tinv6
    transr = st.transr
    tr2 = st.T @ _tinv6(transr)
    if bzh != 0.0:
        tr2[1, :] = tr2[1, :] + bzh * tr2[2, :]
        tr2[3, :] = tr2[3, :] - bzh * tr2[0, :]
    ddpz = (tr2[5, :] * pr - tr2[1, :] * px - tr2[3, :] * py) / pz
    dpxi = np.array([0.0, 1.0, bzhr0, 0.0, 0.0, 0.0])
    dpyi = np.array([-bzhr0, 0.0, 0.0, 1.0, 0.0, 0.0])
    dpz00 = (-pxi * dpxi - pyi * dpyi) / pz00
    dpz00[5] += pr / pz00
    dpxr0 = cp * dpxi + sp * dpz00
    dpz0 = (-pxr0 * dpxr0 - pyi * dpyi) / pz0
    dpz0[5] += pr / pz0
    dxpx = tr2[3, :] * pz0 + py * dpz0 - ddpz * pyi - pz * dpyi
    dxpy = ddpz * pxr0 + pz * dpxr0 - tr2[1, :] * pz0 - px * dpz0
    dxpz = tr2[1, :] * pyi + px * dpyi - tr2[3, :] * pxr0 - py * dpxr0
    dh1r = p * p0 / h1 ** 2

    bsir0 = st.bsir0
    xpzb = xpz + (bsir0 + bzh * 2.0 * al) * pr ** 2
    dxpzb = dxpz.copy()
    dxpzb[5] += 2.0 * (bsir0 + bzh * 2.0 * al) * pr
    bd0 = getattr(st, 'bsi_d0', None)
    if bd0 is not None and (bd0.any() or st.bsi_d1.any()):
        # v2 (B58): bsir0 (端の縦磁場の積分 ψ(x, y)) の軌道についての微分。K64 は bsir0 を定数として扱い、
        # スピン軌道結合 (脱偏極ベクトル) に端の縦磁場の y・x 依存が入らない (追跡のスピンの軌道依存と BEND で 6 倍違う)
        dxpzb = dxpzb + (bd0 + st.bsi_d1 @ tr2) * pr ** 2
    pxh = (pxr0 + px) / pr * 0.5
    pyh = (pyi + py) / pr * 0.5
    dpxh = (tr2[1, :] + dpxr0) / pr * 0.5
    dpxh[5] -= pxh / pr
    dpyh = (tr2[3, :] + dpyi) / pr * 0.5
    dpyh[5] -= pyh / pr
    pzh = 1.0 + _pxy2dpz(pxh, pyh)
    dpzh = -(pxh * dpxh + pyh * dpyh) / pzh
    bp = (xpx * pxh + xpy * pyh + xpzb * pzh) / pr ** 2
    dbp = (dxpx * pxh + xpx * dpxh + dxpy * pyh + xpy * dpyh + dxpzb * pzh + xpzb * dpzh) / pr ** 2
    dbp[5] -= 2.0 * bp / pr
    blx, bly, blz = bp * pxh, bp * pyh, bp * pzh
    btx = xpx / pr ** 2 - blx
    bty = xpy / pr ** 2 - bly
    btz = xpzb / pr ** 2 - blz
    dblx = dbp * pxh + bp * dpxh
    dbly = dbp * pyh + bp * dpyh
    dblz = dbp * pzh + bp * dpzh
    dbtx = dxpx / pr ** 2 - dblx
    dbtx[5] -= 2.0 * xpx / pr ** 3
    dbty = dxpy / pr ** 2 - dbly
    dbty[5] -= 2.0 * xpy / pr ** 3
    dbtz = dxpzb / pr - dblz                       # K64 のまま (pr**2 ではなく pr)
    dbtz[5] -= 2.0 * xpzb / pr ** 3
    ct = h1 * GSPIN
    dct = ct * dh1r
    ct += 1.0
    gx = ct * btx + CL * blx
    gy = ct * bty + CL * bly
    gz = ct * btz + CL * blz
    dgx = ct * dbtx + CL * dblx
    dgx[5] += dct * btx
    dgy = ct * dbty + CL * dbly
    dgy[5] += dct * bty
    dgz = ct * dbtz + CL * dblz
    dgz[5] += dct * btz
    srot = st.srot
    srot[0, 3:9] += dgx @ transr
    srot[1, 3:9] += dgy @ transr
    srot[2, 3:9] += dgz @ transr
    g = math.sqrt(gx * gx + gy * gy + gz * gz)
    if g != 0.0:
        bt = math.sqrt(btx * btx + bty * bty + btz * btz)
        sinu, cosu, dcosu = xsincos(g)
        sinu = sinu / g
        sx, sy, sz = srot[0].copy(), srot[1].copy(), srot[2].copy()
        sw = (sx * gx + sy * gy + sz * gz) / g
        st.gintd = st.gintd + sw[0:3] * sflc() * anp * (bt * h1 * p / al1) ** 2
        sw = -sw * dcosu / g
        sux = sy * gz - sz * gy
        suy = sz * gx - sx * gz
        suz = sx * gy - sy * gx
        srot[0] = cosu * sx + sinu * sux + sw * gx
        srot[1] = cosu * sy + sinu * suy + sw * gy
        srot[2] = cosu * sz + sinu * suz + sw * gz
    srot[0] = cp * srot[0] + sp * srot[2]
    srot[2] = (-sp * srot[0] + srot[2]) / cp       # 更新後の 1 行目を使う (K64 のまま)


# ---------------------------------------------------------------------------
# tchge.f: 座標回転と tbrote の srot 更新
# ---------------------------------------------------------------------------

def tchge_rotate_spin(st, th, ent):
    """tchge: ent → srot(1,:)= c·sx − s·srot(2,:); srot(2,:)= s·sx + c·srot(2,:)。出口は逆回転"""
    if th == 0.0:
        return
    ct, sn = math.cos(th), math.sin(th)
    srot = st.srot
    sx = srot[0].copy()
    if ent:
        srot[0] = ct * sx - sn * srot[1]
        srot[1] = sn * sx + ct * srot[1]
    else:
        srot[0] = ct * sx + sn * srot[1]
        srot[1] = -sn * sx + ct * srot[1]


def tbrote_rr(alg, phig, dtheta, dchi2, ent):
    """tbrote の rr (単位行列からの差。対角は cos−1 の形)"""
    from .lattice import _xsincos
    if phig == 0.0:
        c0, s0 = 1.0, 0.0
    else:
        c0, s0 = math.cos(phig), math.sin(phig)
    s2, _, c2, dc2 = _xsincos(dchi2)
    st, _, ct, dct = _xsincos(dtheta)
    rr = np.zeros((3, 3))
    rr[0, 0] = dct * c0 ** 2 + dc2 * s0 ** 2 + c0 * s0 * s2 * st
    rr[1, 1] = c2 * dct + dc2
    rr[2, 2] = dc2 * c0 ** 2 + dct * s0 ** 2 - c0 * s0 * s2 * st
    if ent:
        rr[0, 1] = s0 * s2 - c0 * c2 * st
        rr[0, 2] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[1, 0] = -ct * s0 * s2 + c0 * st
        rr[1, 2] = c0 * ct * s2 + s0 * st
        rr[2, 0] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[2, 1] = -c0 * s2 - c2 * s0 * st
    else:
        rr[1, 0] = s0 * s2 - c0 * c2 * st
        rr[2, 0] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[0, 1] = -ct * s0 * s2 + c0 * st
        rr[2, 1] = c0 * ct * s2 + s0 * st
        rr[0, 2] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[1, 2] = -c0 * s2 - c2 * s0 * st
    return rr


def tbrote_spin(st, alg, phig, dtheta, dchi2, ent):
    rr = tbrote_rr(alg, phig, dtheta, dchi2, ent)
    srot = st.srot
    from . import compat
    if compat.v2():
        # v2 (B56): 行列積 (I + rr)·srot。回転 (直交行列) なので srot(:,1:3) の直交性が保たれる
        st.srot[:, :] = srot + rr @ srot
        return
    for i in range(9):                            # K64: 行を逐次更新 (2 行目以降は更新後の値を使う)
        srot[0, i] = rr[0, :] @ srot[:, i] + srot[0, i]
        srot[1, i] = rr[1, :] @ srot[:, i] + srot[1, i]
        srot[2, i] = rr[2, :] @ srot[:, i] + srot[2, i]


def tchge_spin(st, dx0, dy, dz, theta, dtheta, alg, phig, ent):
    """tchge の calpol 部 (cod の更新は lattice._tchge_cod が別にやる)"""
    th = theta if phig != 0.0 else theta + dtheta
    if ent:
        tchge_rotate_spin(st, th, True)
        if th != 0.0:
            dx = math.cos(th) * dx0 - math.sin(th) * dy
        else:
            dx = dx0
        if dx != 0.0 or dz != 0.0 or dtheta != 0.0:
            tbrote_spin(st, alg, phig, dtheta, 0.0, True)
    else:
        if th != 0.0:
            dx = math.cos(th) * dx0 - math.sin(th) * dy
        else:
            dx = dx0
        if dx != 0.0 or dz != 0.0 or dtheta != 0.0:
            tbrote_spin(st, alg, phig, dtheta, 0.0, False)
        tchge_rotate_spin(st, th, False)


# ---------------------------------------------------------------------------
# spnorm
# ---------------------------------------------------------------------------

def _abs3(a, b, c):
    return math.hypot(a, math.hypot(b, c))


def spnorm(srot, gintd):
    """spnorm: srot(:,1:3) を直交化し、スピン歳差の不変軸 sps(:,1) と歳差角 smu、減衰率 sdamp = n·gintd を返す。srot は更新される"""
    from .emit6d import _tsvd_tsolva
    srot[:, 1] = srot[:, 1] / np.linalg.norm(srot[:, 1])
    s = srot[0, 0] * srot[0, 1] + srot[1, 0] * srot[1, 1] + srot[2, 0] * srot[2, 1]
    srot[:, 0] = srot[:, 0] - s * srot[:, 1]
    s = _abs3(srot[0, 0], srot[1, 0], srot[2, 0])
    srot[:, 0] = srot[:, 0] / s
    srot[:, 2] = np.cross(srot[:, 0], srot[:, 1])
    sps = np.zeros((3, 3))
    sps[0, 0] = srot[1, 2] - srot[2, 1]
    sps[1, 0] = srot[2, 0] - srot[0, 2]
    sps[2, 0] = srot[0, 1] - srot[1, 0]
    s = _abs3(sps[0, 0], sps[1, 0], sps[2, 0])
    smin = 1e-4
    if s < smin:
        # teigen: 実固有値に対応する固有ベクトル
        w, v = np.linalg.eig(srot[:, 0:3])
        for i in range(3):
            if abs(w[i].imag) < 1e-12:
                vv = np.real(v[:, i])
                sps[:, 0] = vv
                s = _abs3(vv[0], vv[1], vv[2])
                break
    sps[:, 0] = sps[:, 0] / s
    dr = sps[:, 0] - srot[:, 0] * sps[0, 0] - srot[:, 1] * sps[1, 0] - srot[:, 2] * sps[2, 0]
    a = srot[:, 0:3].copy()
    a[0, 0] -= 1.0
    a[1, 1] -= 1.0
    a[2, 2] -= 1.0
    dsps = np.asarray(_tsvd_tsolva(a, dr, 1e-8))
    sps[:, 0] = sps[:, 0] + dsps
    sps[:, 0] = sps[:, 0] / np.linalg.norm(sps[:, 0])
    if abs(min(sps[:, 0])) > abs(max(sps[:, 0])):
        sps[:, 0] = -sps[:, 0]
    sps[:, 1] = 0.0
    if abs(sps[0, 0]) > abs(sps[1, 0]):
        sps[1, 1] = 1.0
        s = sps[1, 0]
    else:
        sps[0, 1] = 1.0
        s = sps[0, 0]
    sps[:, 1] = sps[:, 1] - s * sps[:, 0]
    sps[:, 1] = sps[:, 1] / np.linalg.norm(sps[:, 1])
    sps[:, 2] = np.cross(sps[:, 0], sps[:, 1])
    spsa1 = srot[:, 0] * sps[0, 1] + srot[:, 1] * sps[1, 1] + srot[:, 2] * sps[2, 1]
    cm = float(spsa1 @ sps[:, 1])
    sm = float(spsa1 @ sps[:, 2])
    smu = math.atan2(-sm, cm)
    sdamp = float(sps[:, 0] @ gintd)
    return sps, smu, sdamp


# ---------------------------------------------------------------------------
# srequpol / spdepol
# ---------------------------------------------------------------------------

def _serot(xx, xy, yy, c1, c2, d1, d2, e1, e2):
    tx = 0.5 * math.atan2(2.0 * xy, xx - yy)
    c = math.cos(tx)
    s = math.sin(tx)
    em1 = c ** 2 * xx + s ** 2 * yy + 2.0 * c * s * xy
    em2 = c ** 2 * yy + s ** 2 * xx - 2.0 * c * s * xy
    c1 = c * c1 + s * c2
    c2 = (-s * c1 + c2) / c                       # 更新後の c1 (K64 のまま)
    d1 = c * d1 + s * d2
    d2 = (-s * d1 + d2) / c
    e1 = c * e1 + s * e2
    e2 = (-s * e1 + e2) / c
    return c1, c2, d1, d2, e1, e2, em1, em2


def _cexp1(z):
    z = complex(z)
    if abs(z) < 2e-4:
        return z * (1.0 + 0.5 * z * (1.0 + z / 3.0 * (1.0 + z / 4.0 * (1.0 + z / 5.0))))
    return complex(np.exp(z) - 1.0)


def _cerfc(z):
    """gamma.f cerfc: |z|>18 は漸近展開、そうでなければ erfc"""
    z = complex(z)
    if abs(z) > CEZLIM:
        x = 0.5 / z ** 2
        base = 1.0 if z.real == 0.0 else 1.0 - math.copysign(1.0, z.real)
        with np.errstate(all='ignore'):
            return base + (1.0 - x * (1.0 - 3.0 * x * (1.0 - 5.0 * x * (1.0 - 7.0 * x * (1.0 - 9.0 * x)))))\
                * np.exp(-0.5 / x) / z / M_SQRTPI
    from scipy.special import erfc
    return complex(erfc(z))


def _mat(rows):
    return np.array(rows, dtype=complex)


def _indn(i, i1, i2, idx, imx, is_):
    ind = [0] * 13
    ind[2 * i - 2] = i1
    ind[2 * i - 1] = i2
    ind[i + 5] = idx
    ind[i + 8] = imx
    ind[12] = is_
    return tuple(ind)


def _erfcx(z):
    """exp(z²)·erfc(z) (複素数、桁あふれなし) = w(i z) (Faddeeva 関数)"""
    from scipy.special import wofz
    return complex(wofz(1j * complex(z)))


def _spintrm(rm, dx, amx, ams, ssprd):
    from . import compat
    _v2 = compat.v2()
    for ind in list(rm.keys()):
        mat = rm[ind]
        cphi = complex(-(ind[6] * dx[0] + ind[7] * dx[1] + ind[8] * dx[2]),
                       ind[9] * amx[0] + ind[10] * amx[1] + ind[11] * amx[2] + ind[12] * ams)
        with np.errstate(all='ignore'):
            c1 = _cexp1(cphi)
            cm0 = -mat / c1 if c1 != 0 else -mat * complex(float('inf'), 0.0)
            a = math.sqrt(2.0) * abs(ind[12]) * ssprd
            if a == 0.0:
                cm = cm0
            else:
                cb = -_cexp1(-cphi) / a
                if _v2:
                    # v2 (B55): exp(cb²)·erfc(cb) を桁あふれしない形 (Faddeeva 関数 w(i·cb)) で計算する。
                    # K64 は |cb|>18 で exp(cb²-cφ) → 0 と erfc の漸近式の exp(-cb²) → ∞ を別々に計算して NaN になる
                    cm = -mat * M_2PI * M_SQRTPI * np.exp(-cphi) * _erfcx(cb) / a
                else:
                    cm = -mat * M_2PI * M_SQRTPI * np.exp(cb ** 2 - cphi) * _cerfc(cb) / a
        rm[ind] = cm0
        k0 = ind[:6] + (0,) * 7
        if k0 in rm:
            rm[k0] = rm[k0] + cm
        else:
            rm[k0] = cm.copy() if hasattr(cm, 'copy') else cm


def _spdotrm(rma, rmb, rmc):
    for ia, ma in rma.items():
        for ib, mb in rmb.items():
            k = tuple(x + y for x, y in zip(ia, ib))
            t = ma @ mb
            if k in rmc:
                rmc[k] = rmc[k] + t
            else:
                rmc[k] = t


def spdepol(gxr, gxi, gyr, gyi, gzr, gzi, e1, e2, dx, ssprd, amx, ams):
    """spdepol: 2 次までの摂動 (K64 は spcalcres を k=2 だけ呼ぶ)。rmd(3,3,1) = rm1 を返す"""
    rm1 = {}
    cI = 1j
    for i in (1, 2, 3):
        gx = complex(gxr[i - 1], gxi[i - 1])
        gxc = gx.conjugate()
        gy = complex(gyr[i - 1] + gzi[i - 1], gzr[i - 1] - gyi[i - 1])
        gyc = gy.conjugate()
        gz = complex(gyr[i - 1] - gzi[i - 1], gzr[i - 1] + gyi[i - 1])
        gzc = gz.conjugate()
        m1 = 0.25 * gzc * np.array([0, -cI, 1, cI, 0, 0, -1, 0, 0], dtype=complex).reshape((3, 3), order='F')
        m2 = 0.5 * gxc * np.array([0, 0, 0, 0, 0, -1, 0, 1, 0], dtype=complex).reshape((3, 3), order='F')
        m3 = 0.25 * gy * np.array([0, cI, 1, -cI, 0, 0, -1, 0, 0], dtype=complex).reshape((3, 3), order='F')
        rm1[_indn(i, 0, 1, 1, -1, -1)] = m1
        rm1[_indn(i, 0, 1, 1, -1, 0)] = m2
        rm1[_indn(i, 0, 1, 1, -1, 1)] = m3
        rm1[_indn(i, 1, 0, 1, 1, -1)] = m3.conj()
        rm1[_indn(i, 1, 0, 1, 1, 0)] = m2.conj()
        rm1[_indn(i, 1, 0, 1, 1, 1)] = m1.conj()
    rmi1 = {k: v.copy() for k, v in rm1.items()}
    _spintrm(rmi1, dx, amx, ams, ssprd)
    # spcalcres(m=2)
    rm2 = {}
    _spdotrm(rm1, rm1, rm2)
    for k in rm2:
        rm2[k] = rm2[k] / 2.0
    rmi2 = {k: v.copy() for k, v in rm2.items()}
    _spdotrm(rm1, rmi1, rmi2)
    _spintrm(rmi2, dx, amx, ams, ssprd)
    rmd = np.zeros((3, 3))
    z = np.zeros((3, 3), dtype=complex)
    for i in (1, 2, 3):
        a20 = rmi2.get(_indn(i, 2, 0, 0, 0, 0), z)
        a21 = rmi2.get(_indn(i, 1, 1, 0, 0, 0), z)
        se12 = (e1[i - 1] + e2[i - 1]) * 0.5
        de12 = (e1[i - 1] - e2[i - 1]) * 0.5
        rmd = rmd + se12 * np.real(a21) + 2.0 * de12 * np.real(a20)
    return rmd


def srequpol(srot, sps, params, demit, sdamp, r):
    """srequpol: 平衡偏極 equpol(3) と脱偏極行列 rm1(3,3)。
    params: dict(emit=(ex,ey,ez), damp=(3), tunes=(3), smu)。demit は規格化モード座標のビーム行列 (21 個の packed)"""
    smu = params['smu']
    drot = sps.T @ srot[:, 3:9] @ r
    c1, c2, c3, c4, c5, c6 = drot[0]
    d1, d2, d3, d4, d5, d6 = drot[1]
    e1, e2, e3, e4, e5, e6 = drot[2]
    c1, c2, d1, d2, e1, e2, dex1, dex2 = _serot(demit[0], demit[1], demit[2], c1, c2, d1, d2, e1, e2)
    c3, c4, d3, d4, e3, e4, dey1, dey2 = _serot(demit[5], demit[8], demit[9], c3, c4, d3, d4, e3, e4)
    c5, c6, d5, d6, e5, e6, dez1, dez2 = _serot(demit[14], demit[19], demit[20], c5, c6, d5, d6, e5, e6)
    emit1 = params['emit']
    damp = np.abs(np.asarray(params['damp'], dtype=float))
    ssprd = math.sqrt(float(np.dot(drot[0] ** 2, [emit1[0], emit1[0], emit1[1], emit1[1], emit1[2], emit1[2]])))
    amu = np.asarray(params['tunes'], dtype=float) * M_2PI
    rm1 = spdepol([c1, c3, c5], [c2, c4, c6], [d1, d3, d5], [d2, d4, d6], [e1, e3, e5], [e2, e4, e6],
                  [dex1, dey1, dez1], [dex2, dey2, dez2], damp, ssprd, amu, smu)
    with np.errstate(all='ignore'):
        equpol = np.array([PST * sdamp / (sdamp - rm1[0, 0]), 0.0, 0.0])
    return equpol, rm1


def _cross_mat(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def slim_depolarization(els, ph, cod, trf0=0.0, rfsw=False, radcod=False):
    """線形のスピン軌道結合 (SLIM / Derbenev-Kondratenko の脱偏極項) による 1 周あたりの脱偏極率 (v2 の Polarization が使う)。

    要素 i ごとに emit6d.tturne (POL) で、軌道の行列 T_i (減衰込み)、その要素の放射の拡散行列 B_i (要素の出口での共分散)、
    スピンの回転 R_i、スピン軌道結合 D_i (軌道のずれ u による余分の回転ベクトル D_i u。K64 の回転の向きは ds = s × g なので
    s' = (1 − [D_i u]×) R_i s。追跡 (track6 のスピン) の 1 周の応答と 2〜6% で一致することを確かめた) を取る。
    スピンの周期解 n(u) = n0 + N u (N は 3x6) は  N_{i+1} = (R_i N_i + [n0_{i+1}]× D_i) T_i^{-1}、N_N = N_0 を満たす。
    放射で軌道が ξ (共分散 B_i) だけ跳ぶと、スピンは変わらず n が N ξ だけ動くので、偏極は 1 周あたり
        rate = (1/2) Σ_i tr(N⊥ B_i N⊥ᵀ)   (N⊥ は n0 に垂直な成分)
    だけ減る。これは DK の式の (11/18)|γ ∂n/∂γ|² の項と同じ (B_i の δδ 成分 = (55/24√3) r_e ƛ γ^5/|ρ|³ ds)。
    戻り dict(rate, n0, N0 (入口の ∂n/∂u)、rev (平衡分布での n の広がり (1/2) tr(N⊥ Σ N⊥ᵀ) を呼び出し側が作るための N0))。
    限界: 軌道の振幅について 1 次 (シンクロトロン側帯などの高次の共鳴は入らない)、運動論的偏極 (b·γ∂n/∂γ) の項は入らない"""
    from . import emit6d as _e6
    old = _e6.POL_MODE
    _e6.POL_MODE = True
    secs = []
    c = np.array(cod, dtype=float)
    try:
        for el in els:
            st = _e6.tturne([el], c, ph, trf0=trf0, rfsw=rfsw, radcod=radcod)
            secs.append((st.T + st.D, st.B.copy(), st.srot[:, 0:3].copy(), st.srot[:, 3:9].copy()))
            c = st.cod.copy()
    finally:
        _e6.POL_MODE = old
    rtot = np.eye(3)
    for _T, _B, R, _D in secs:
        rtot = R @ rtot
    w, v = np.linalg.eig(rtot)
    k = int(np.argmin(np.abs(w - 1.0)))
    n0 = np.real(v[:, k])
    n0 = n0 / np.linalg.norm(n0)
    if abs(min(n0)) > abs(max(n0)):
        n0 = -n0
    tinv = [np.linalg.inv(T) for T, _B, _R, _D in secs]

    def prop(N):
        n = n0
        out = []
        for (T, _B, R, D), ti in zip(secs, tinv):
            n = R @ n
            N = (R @ N + _cross_mat(n) @ D) @ ti
            out.append((N, n))
        return N, out
    cvec, _ = prop(np.zeros((3, 6)))
    a = np.zeros((18, 18))
    for j in range(18):
        e = np.zeros(18)
        e[j] = 1.0
        a[:, j] = (prop(e.reshape(3, 6))[0] - cvec).ravel()
    try:
        nn0 = np.linalg.solve(np.eye(18) - a, cvec.ravel()).reshape(3, 6)
    except np.linalg.LinAlgError:
        nn0 = np.linalg.lstsq(np.eye(18) - a, cvec.ravel(), rcond=None)[0].reshape(3, 6)
    _, out = prop(nn0)
    rate = 0.0
    for (_T, B, _R, _D), (N, n) in zip(secs, out):
        npp = (np.eye(3) - np.outer(n, n)) @ N
        rate += 0.5 * float(np.trace(npp @ B @ npp.T))
    return dict(rate=rate, n0=n0, N0=nn0)


def spin_results(st, em, ph, els=None, cod=None, trf0=0.0, rfsw=False, radcod=False) -> dict:
    """temit.f tecalc の calpol 部: srot1 = srot(:,1:3) → spnorm → srequpol。偏極時間 = 1/sdamp/revolution frequency [s]。
    v2 (B55) で els・cod が渡されたとき: 平衡偏極度は線形のスピン軌道結合 (slim_depolarization) の脱偏極率から
    P = (8/5√3)·sdamp/(sdamp + rate) とする (K64 の spdepol は NaN、桁あふれを直しても追跡と 3 桁以上合わない)"""
    from . import compat
    srot = st.srot
    srot1 = srot[:, 0:3].copy()
    sps, smu, sdamp = spnorm(srot, st.gintd)
    params = dict(emit=(em['emittance_x'], em['emittance_y'], em['emittance_z']),
                  damp=em['damping_rate'], tunes=em['tunes'], smu=smu)
    equpol, spm = srequpol(srot, sps, params, em['beamn'], sdamp, em['normal_r'])
    frev = ph.omega0 / M_2PI
    with np.errstate(all='ignore'):
        taup = 1.0 / sdamp / frev if sdamp != 0.0 else float('inf')
    res = dict(srot=srot.copy(), srot1=srot1, sps=sps, smu=smu, sdamp=sdamp, nup=smu / M_2PI,
               equpol=equpol, spm=spm, taup=taup, gintd=st.gintd.copy())
    if compat.v2() and els is not None and cod is not None:
        res['equpol_k64'] = equpol                      # K64 の spdepol の式 (桁あふれだけ直したもの)。参考値
        try:
            sl = slim_depolarization(els, ph, cod, trf0=trf0, rfsw=rfsw, radcod=radcod)
            rate = sl['rate']
            with np.errstate(all='ignore'):
                pol = PST * sdamp / (sdamp + rate) if (sdamp + rate) != 0.0 else float('nan')
            npp = (np.eye(3) - np.outer(sl['n0'], sl['n0'])) @ sl['N0']
            res.update(equpol=np.array([pol, 0.0, 0.0]), depol_rate=rate, spin_orbit=sl['N0'],
                       spin_spread=0.5 * float(np.trace(npp @ em['beam_matrix'] @ npp.T)),
                       taudep=(1.0 / rate / frev if rate > 0.0 else float('inf')))
        except Exception as exc:            # 求まらないときは NaN (K64 と同じ表示) にして理由を残す
            res.update(equpol=np.array([float('nan'), 0.0, 0.0]), depol_error=str(exc))
    return res
