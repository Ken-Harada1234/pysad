"""
BEAMBEAM 要素 (強弱ビームのビーム・ビーム相互作用) の追跡。K64 の bbstrhl.f (beambeam, bbinit)、bbtool.f (cerrf, gauinv)、
matrix.f (包絡行列) の直訳。ビームストラールング (BSTRL>0) は乱数 (tran) の引き方まで直訳。光子の出力 (photons) は未対応。ルミノシティは LuminosityMonitor に渡す。光学 (Twiss) には効かない。

K64 の癖 (直訳で残す):
  - 傾きの回転 (tilt) は sint=xangle(5) (= 傾き角そのもの)、cost=xangle(6) (= sin(傾き角)) を使う (bbinit が書くのは
    xangle(6)=sin、xangle(7)=cos)。tilt≠0 (強ビームに x-y の結合がある) のときだけ効く。
  - cerrf の 3.2 は単精度の定数 (H=1/(3.2*Q))。
  - 丸いビーム (|σx²-σy²|/σx² ≤ 1e-5) で x=y=0 の粒子は 0/0 が NaN になり、limitnan で px=+1e4 になる。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import compat
from .emit6d import ELMASS, ELRADI, FINEST, CVELOC
from .track6n import tran, lcg, pxy2dpz

NJ = dict(cache=True, error_model='numpy')
LUM_CB = None            # TrackParticles が設定する関数 f(lum) (LuminosityMonitor)
NSLIMAX = 500
SQRPI = math.sqrt(math.pi)
_H32 = float(np.float32(3.2))

# 要素キーワード → p_in の位置 (1 始まり、initb1.f ky_*_BEAM)
_KEYS = ['AX', 'BX', 'AY', 'BY', 'R1', 'R2', 'R3', 'R4', 'EX', 'EPX', 'EY', 'EPY', 'ZX', 'ZPX', 'ZY', 'ZPY',
         'DX', 'DPX', 'DY', 'DPY', 'XANGLE', 'EMIX', 'EMIY', 'DP', 'AZ', 'DZ', 'SIGZ', 'SLICE', 'NP', 'STURN',
         'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'R22', 'R23', 'R24', 'R25', 'R26', 'R33', 'R34', 'R35', 'R36',
         'R44', 'R45', 'R46', 'R55', 'R56', 'R66', 'COUPLE', 'BSTRL']


_HBAR = float(np.float32(6.58211928e-16))      # besseltab の hbar [eV s] は単精度の定数のまま
LBTAB = 1000
_TABLE = None
_TABLE_V2 = None


def _chebev(a, b, c, m, x):
    d = 0.0
    dd = 0.0
    y = (2.0 * x - a - b) / (b - a)
    y2 = 2.0 * y
    for j in range(m, 1, -1):
        sv = d
        d = y2 * d - dd + c[j - 1]
        dd = sv
    return y * d - dd + 0.5 * c[0]


_C1 = [-1.142022680371168, 6.5165112670737e-3, 3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13]
_C2 = [1.843740587300905, -7.68528408447867e-2, 1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10,
       -1.702e-13, -1.49e-15]


def _beschb(x):
    xx = 8.0 * x * x - 1.0
    gam1 = _chebev(-1.0, 1.0, _C1, 7, xx)
    gam2 = _chebev(-1.0, 1.0, _C2, 8, xx)
    return gam1, gam2, gam2 - x * gam1, gam2 + x * gam1


def _bessk(x, xnu):
    """bbtool.f bessik の rk (修正ベッセル関数 K_xnu(x))。Numerical Recipes 版の直訳 (EPS=1e-10 は単精度の定数)"""
    EPS = float(np.float32(1.0e-10))
    FPMIN = float(np.float32(1.0e-30))
    MAXIT = 10000
    XMIN = 2.0
    nl = int(xnu + 0.5)
    xmu = xnu - nl
    xmu2 = xmu * xmu
    xi = 1.0 / x
    xi2 = 2.0 * xi
    h = xnu * xi
    if h < FPMIN:
        h = FPMIN
    b = xi2 * xnu
    d = 0.0
    c = h
    for _ in range(MAXIT):
        b = b + xi2
        d = 1.0 / (b + d)
        c = b + 1.0 / c
        de = c * d
        h = de * h
        if abs(de - 1.0) < EPS:
            break
    if x < XMIN:
        x2 = 0.5 * x
        pimu = math.pi * xmu
        fact = 1.0 if abs(pimu) < EPS else pimu / math.sin(pimu)
        d = -math.log(x2)
        e = xmu * d
        fact2 = 1.0 if abs(e) < EPS else math.sinh(e) / e
        gam1, gam2, gampl, gammi = _beschb(xmu)
        ff = fact * (gam1 * math.cosh(e) + gam2 * fact2 * d)
        sm = ff
        e = math.exp(e)
        p = 0.5 * e / gampl
        q = 0.5 / (e * gammi)
        c = 1.0
        d = x2 * x2
        sum1 = p
        for i in range(1, MAXIT + 1):
            ff = (i * ff + p + q) / (i * i - xmu2)
            c = c * d / i
            p = p / (i - xmu)
            q = q / (i + xmu)
            de = c * ff
            sm = sm + de
            sum1 = sum1 + c * (p - i * ff)
            if abs(de) < abs(sm) * EPS:
                break
        rkmu = sm
        rk1 = sum1 * xi2
    else:
        b = 2.0 * (1.0 + x)
        d = 1.0 / b
        delh = d
        h = delh
        q1 = 0.0
        q2 = 1.0
        a1 = 0.25 - xmu2
        c = a1
        q = c
        a = -a1
        sm = 1.0 + q * delh
        for i in range(2, MAXIT + 1):
            a = a - 2 * (i - 1)
            c = -a * c / i
            qnew = (q1 - b * q2) / a
            q1 = q2
            q2 = qnew
            q = q + c * qnew
            b = b + 2.0
            d = 1.0 / (b + a * d)
            delh = (b * d - 1.0) * delh
            h = h + delh
            dels = q * delh
            sm = sm + dels
            if abs(dels / sm) < EPS:
                break
        h = a1 * h
        rkmu = math.sqrt(math.pi / (2.0 * x)) * math.exp(-x) / sm
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi
    for i in range(1, nl + 1):
        rktemp = (xmu + i) * xi2 * rk1 + rkmu
        rkmu = rk1
        rk1 = rktemp
    return rkmu


def _powi32(x, n):
    """単精度の x**n (gfortran の整数べき。2 進展開)"""
    x = np.float32(x)
    r = np.float32(1.0)
    while n:
        if n & 1:
            r = np.float32(r * x)
        x = np.float32(x * x)
        n >>= 1
    return float(r)


def besseltab():
    """bbstrhl.f mkbesseltab: (xi, SwI) の 1 始まりの表 (長さ 1002)。光子のエネルギー分布 S(x) の積分の逆引き。
    1.1**40、1.02、5./3. は単精度の定数のまま"""
    global _TABLE, _TABLE_V2
    if compat.v2():
        if _TABLE_V2 is None:
            _TABLE_V2 = _besseltab_v2()
        return _TABLE_V2
    if _TABLE is not None:
        return _TABLE
    xnu = float(np.float32(5.0) / np.float32(3.0))
    x = _powi32(1.1, 40)
    c102 = float(np.float32(1.02))
    xi = np.zeros(LBTAB + 2)
    sw = np.zeros(LBTAB + 2)
    swi = np.zeros(LBTAB + 2)
    rk0 = 0.0
    for i in range(1, LBTAB + 1):
        x0 = x
        x = x / c102
        x2 = (x + x0) * 0.5
        dx6 = (x0 - x) / 6.0
        xi[i] = math.log(x)
        rk2 = _bessk(x2, xnu)
        rk = _bessk(x, xnu)
        if i > 1:
            sw[i] = sw[i - 1] + (rk0 + 4 * rk2 + rk) * dx6
            swi[i] = swi[i - 1] + sw[i] * (x0 - x)
        rk0 = rk
    top = swi[LBTAB]
    for i in range(1, LBTAB + 1):
        swi[i] = swi[i] / top
    _TABLE = (xi, swi)
    return _TABLE


def _besseltab_v2():
    """v2 (B335): 光子エネルギーの逆引き表。K64 (`mkbesseltab`) は累積 SwI を Sw(x_i)·(x_{i-1} − x_i) の矩形則 (区間の下端の Sw) で作り、
    x_min = 1.1^40/1.02^1000 = 1.1e-7 より下の光子 (全体の 0.6%) を捨てて SwI(1000) = 1 に規格化する。
    v2 は累積を台形則 (区間の両端の平均) で作り、x_min より下の裾は小 x の漸近形 Sw ∝ x^(-2/3) から 3·x_min·Sw(x_min) と見積って全体に入れる
    (SwI(1000) < 1)。u が SwI(1000) より大きい光子は、表の外側 (ξ → 0、1001 番目の点を ln ξ が 40 小さい点に置く) で、エネルギーはほぼ 0。
    格子 (x = 1.1^40/1.02^i) と Sw の Simpson 積分は K64 と同じ。理論値 <ξ> = 0.30792、<ξ²> = 0.40741 に対し、表から (線形補間込みで) 求めたモーメントは K64 の表が +1.60%・+2.62%、この表が +0.016%・+0.046%"""
    xnu = float(np.float32(5.0) / np.float32(3.0))
    x = _powi32(1.1, 40)
    c102 = float(np.float32(1.02))
    xi = np.zeros(LBTAB + 2)
    sw = np.zeros(LBTAB + 2)
    swi = np.zeros(LBTAB + 2)
    xs = np.zeros(LBTAB + 2)
    rk0 = 0.0
    for i in range(1, LBTAB + 1):
        x0 = x
        x = x / c102
        x2 = (x + x0) * 0.5
        dx6 = (x0 - x) / 6.0
        xs[i] = x
        xi[i] = math.log(x)
        rk2 = _bessk(x2, xnu)
        rk = _bessk(x, xnu)
        if i > 1:
            sw[i] = sw[i - 1] + (rk0 + 4 * rk2 + rk) * dx6
            swi[i] = swi[i - 1] + 0.5 * (sw[i] + sw[i - 1]) * (x0 - x)
        rk0 = rk
    total = swi[LBTAB] + 3.0 * xs[LBTAB] * sw[LBTAB]
    for i in range(1, LBTAB + 1):
        swi[i] = swi[i] / total
    xi[LBTAB + 1] = xi[LBTAB] - 40.0
    swi[LBTAB + 1] = 1.0
    return (xi, swi)


def gauinv(p0):
    """bbtool.f GAUINV (最大相対誤差 6.4e-4 の近似。そのまま)"""
    PP1, QQ2, QQ3, QQ4, QQ5 = 0.334624883253, 0.090230446775, 0.049905685242, 0.027852994157, 0.015645650215
    A3, A2, A1, A0 = 4.5585614e+01, 2.1635544e+00, 2.7724523e+00, 2.5050240e+00
    B4, B3, B2, B1, B0 = 4.0314354e+02, -2.7713713e+02, 7.9731883e+01, -1.4946512e+01, 2.2157257e+00
    C4, C3, C2, C1, C0 = 4.1394487e+03, -1.5585873e+03, 2.4648581e+02, -2.4719139e+01, 2.4335936e+00
    D4, D3, D2, D1, D0 = 4.0895693e+04, -8.5400893e+03, 7.4942805e+02, -4.1028898e+01, 2.6346872e+00
    E4, E3, E2, E1, E0 = 3.9399134e+05, -4.6004775e+04, 2.2566998e+03, -6.8317697e+01, 2.8224654e+00
    F0, F1, F2 = -8.1807613e-02, -2.8358733e+00, 1.4902469e+00
    p = p0 - 0.5
    p1 = abs(p)
    if p1 < PP1:
        p2 = p * p
        return (((A3 * p2 + A2) * p2 + A1) * p2 + A0) * p
    q = 0.5 - p1
    if q > QQ2:
        g = (((B4 * q + B3) * q + B2) * q + B1) * q + B0
    elif q > QQ3:
        g = (((C4 * q + C3) * q + C2) * q + C1) * q + C0
    elif q > QQ4:
        g = (((D4 * q + D3) * q + D2) * q + D1) * q + D0
    elif q > QQ5:
        g = (((E4 * q + E3) * q + E2) * q + E1) * q + E0
    elif q > 0.0:
        t = math.sqrt(-2.0 * math.log(q))
        g = t + F0 + F1 / (F2 + t)
    else:
        return 0.0
    return -g if p < 0.0 else g


def _norm_axis(env, vec):
    """matrix.f benv_norm_axis: x-y の結合を消す回転 (5×5 包絡とその重心を回す)。傾き角を返す"""
    t = math.atan2(2.0 * env[0, 2], env[0, 0] - env[2, 2]) * 0.5
    c, s = math.cos(t), math.sin(t)
    u = np.zeros((5, 5))
    u[0, 0] = u[1, 1] = u[2, 2] = u[3, 3] = c
    u[4, 4] = 1.0
    u[0, 2] = u[1, 3] = s
    u[2, 0] = u[3, 1] = -s
    return t, u @ env @ u.T, u @ vec


def m_B(ax, bx, ay, by, az, bz):
    """matrix.f m_B_set_6d"""
    B = np.zeros((6, 6))
    sb = math.sqrt(bx)
    B[0, 0], B[1, 0], B[1, 1] = sb, -ax / sb, 1.0 / sb
    sb = math.sqrt(by)
    B[2, 2], B[3, 2], B[3, 3] = sb, -ay / sb, 1.0 / sb
    sb = math.sqrt(bz)
    B[4, 4], B[5, 4], B[5, 5] = sb, -az / sb, 1.0 / sb
    return B


def m_R(r1, r2, r3, r4):
    """matrix.f m_R_set_6d"""
    rmu = math.sqrt(1.0 - r1 * r4 + r2 * r3)
    R = np.zeros((6, 6))
    R[0, 0] = R[1, 1] = R[2, 2] = R[3, 3] = rmu
    R[4, 4] = R[5, 5] = 1.0
    R[2, 0], R[2, 1], R[3, 0], R[3, 1] = -r1, -r2, -r3, -r4
    R[0, 2], R[0, 3], R[1, 2], R[1, 3] = r4, -r2, -r3, r1
    return R


def m_H(ex, epx, ey, epy, zx, zpx, zy, zpy):
    """matrix.f m_H_set_6d"""
    detx = zx * epx - zpx * ex
    dety = zy * epy - zpy * ey
    a2 = 1.0 - detx - dety
    if not a2 > 0.0:
        raise ValueError(f'det Hx + det Hy > 1 ({detx}, {dety})')
    a = math.sqrt(a2)
    ai = 1.0 / (1.0 + a)
    H = np.zeros((6, 6))
    H[0, 0] = H[1, 1] = 1.0 - detx * ai
    H[2, 2] = H[3, 3] = 1.0 - dety * ai
    H[4, 4] = H[5, 5] = a
    H[0, 4], H[0, 5], H[1, 4], H[1, 5] = zx, ex, zpx, epx
    H[2, 4], H[2, 5], H[3, 4], H[3, 5] = zy, ey, zpy, epy
    H[4, 0], H[4, 1], H[5, 0], H[5, 1] = -epx, ex, zpx, -zx
    H[4, 2], H[4, 3], H[5, 2], H[5, 3] = -epy, ey, zpy, -zy
    H[0, 2] = (ex * zpy - epy * zx) * ai
    H[0, 3] = (ey * zx - ex * zy) * ai
    H[1, 2] = (epx * zpy - epy * zpx) * ai
    H[1, 3] = (ey * zpx - epx * zy) * ai
    H[2, 0] = (ey * zpx - epx * zy) * ai
    H[2, 1] = (-ey * zx + ex * zy) * ai
    H[3, 0] = (epy * zpx - epx * zpy) * ai
    H[3, 1] = (ex * zpy - epy * zx) * ai
    return H


def setup(params: dict, momentum_eV: float, photon=None):
    """bbinit。戻り値は追跡カーネルに渡す辞書 (nslice=0 なら相互作用なし)"""
    p = {k: float(params.get(k, 0.0) or 0.0) for k in _KEYS}
    nsli = int(p['SLICE'])
    if nsli > NSLIMAX:
        raise ValueError(f'BEAMBEAM: SLICE={nsli} は最大 {NSLIMAX} を超える')
    out = dict(nslice=max(nsli, 0))
    if nsli <= 0:
        return out
    if all(p[k] > 0.0 for k in ('R11', 'R22', 'R33', 'R44', 'R55', 'R66')):
        env = np.zeros((6, 6))
        k = 0
        for i in range(6):
            for j in range(i, 6):
                env[i, j] = env[j, i] = p[('R%d%d' % (i + 1, j + 1))]
                k += 1
    else:
        emiz2 = p['DP']**2 * p['SIGZ']**2 - p['AZ']**2
        if emiz2 <= 0.0:
            raise ValueError('BEAMBEAM: 縦のエミッタンスが 0 以下 (DP^2 SIGZ^2 - AZ^2)')
        emiz = math.sqrt(emiz2)
        bz = p['SIGZ']**2 / emiz
        az = -p['AZ'] / emiz
        B = m_B(p['AX'], p['BX'], p['AY'], p['BY'], az, bz)
        R = m_R(p['R1'], p['R2'], p['R3'], p['R4'])
        H = m_H(p['EX'], p['EPX'], p['EY'], p['EPY'], p['ZX'], p['ZPX'], p['ZY'], p['ZPY'])
        E = np.diag([p['EMIX'], p['EMIX'], p['EMIY'], p['EMIY'], p['DP'] * p['SIGZ'], p['DP'] * p['SIGZ']])
        HRB = H @ R @ B
        env = HRB @ E @ HRB.T
    xang = p['XANGLE']
    if xang != 0.0:
        cosxi = 1.0 / math.cos(xang)
        tanx = math.tan(xang)
        C = np.zeros((6, 6))
        C[0, 0] = C[2, 2] = C[5, 5] = 1.0
        C[1, 1] = C[3, 3] = C[4, 4] = cosxi
        C[0, 4] = tanx
        C[5, 1] = -tanx
        env = C @ env @ C.T
    gambet = momentum_eV / ELMASS
    gamma = math.sqrt(gambet**2 + 1.0)
    sinx, cosx, tanx = math.sin(xang), math.cos(xang), math.tan(xang)
    ce = p['NP'] * ELRADI / gamma
    cod = np.array([p['DX'], p['DPX'] / cosx, p['DY'], p['DPY'] / cosx, p['DZ'] / cosx])
    # z スライス: x,px,y,py,pz の 5 次元 (z を固定した条件付きの包絡)
    inv = np.linalg.inv(env)
    sub = np.zeros((5, 5))
    sub[:4, :4] = inv[:4, :4]
    vsub = np.zeros(5)
    for i in range(4):
        vsub[i] = inv[4, i]
        sub[4, i] = inv[5, i]
        sub[i, 4] = inv[i, 5]
    vsub[4] = inv[4, 5]
    sub[4, 4] = inv[5, 5]
    if np.min(np.linalg.eigvalsh(0.5 * (sub + sub.T))) < 0.0:
        raise ValueError('BEAMBEAM: 5 次元の包絡に 0 か負の固有値がある')
    benv5 = np.linalg.inv(sub)
    vcen = -(benv5 @ vsub)
    tilt, benv5, vcen = _norm_axis(benv5, vcen)
    xangle = np.array([xang, sinx, cosx, tanx, tilt, math.sin(tilt), math.cos(tilt), 0.0])
    # スライスの中心 (σz 単位、前から後ろへ符号を反転済み)
    zs = np.zeros(NSLIMAX * 2 + 2)           # 1 始まりの添字をそのまま使う
    for i in range(1, nsli):
        zs[i + NSLIMAX] = gauinv(1.0 / nsli * i)
        if i == 1:
            zs[i] = -math.exp(-zs[i + NSLIMAX]**2 / 2.0) / math.sqrt(2.0 * math.pi) * nsli
        else:
            zs[i] = (math.exp(-zs[i + NSLIMAX - 1]**2 / 2.0) - math.exp(-zs[i + NSLIMAX]**2 / 2.0)) \
                / math.sqrt(2.0 * math.pi) * nsli
    if nsli == 1:
        zs[1] = 0.0
    else:
        zs[nsli] = math.exp(-zs[nsli + NSLIMAX - 1]**2 / 2.0) / math.sqrt(2.0 * math.pi) * nsli
    for i in range(1, nsli + 1):
        zs[i] = -zs[i]
    bstrl = int(p['BSTRL']) > 0
    cnbs = 5.0 * math.sqrt(3.0) * FINEST * gamma / 6.0
    cubs = 1.5 * _HBAR * CVELOC * gamma**2 / ELMASS
    zb = zs[NSLIMAX:NSLIMAX + nsli + 1].copy()           # スライスの境界 (σz 単位)。zs[NSLIMAX] は 100 (符号を反転しない)
    zb[1:nsli] = -zb[1:nsli]                          # mkbesseltab 直前の符号反転 (zslice(i+nslimax) も反転される)
    zb[0] = 100.0
    zb[nsli] = -100.0
    dsz = np.zeros(nsli)
    sg = p['SIGZ']
    for k in range(1, nsli + 1):
        dsz[k - 1] = (zb[k - 1] - zb[k]) * sg * 0.5
        if k == 1:
            dsz[k - 1] = (5.0 - zb[k]) * sg * 0.5
        if k == nsli:
            dsz[k - 1] = (zb[k - 1] + 5.0) * sg * 0.5
    out.update(bstrl=bstrl, cnbs=cnbs, cubs=cubs, dsz=dsz, photon=photon if bstrl else None, momentum=momentum_eV)
    out.update(gamma=gamma, gambet=gambet, ce=ce, xangle=xangle, cod=cod, benv5=benv5, vcen=vcen,
               zc=zs[1:nsli + 1].copy(), sigz=p['SIGZ'])
    return out


@njit(**NJ)
def cerrf(xx, yy):
    """bbtool.f CERRF (複素誤差関数 w(z) の Koelbig 版)。(wx, wy)"""
    CC = 1.1283791670955
    XLIM = 5.33
    YLIM = 4.29
    if np.isnan(xx) or np.isnan(yy):
        return 0.0, 0.0
    x = abs(xx)
    y = abs(yy)
    rx = np.zeros(34)
    ry = np.zeros(34)
    if y < YLIM and x < XLIM:
        q = (1.0 - y / YLIM) * math.sqrt(1.0 - (x / XLIM)**2)
        h = 1.0 / (_H32 * q)
        nc = 7 + int(23.0 * q)
        xl = h**(1 - nc)
        xh = y + 0.5 / h
        yh = x
        nu = 10 + int(21.0 * q)
        rx[nu + 1] = 0.0
        ry[nu + 1] = 0.0
        for n in range(nu, 0, -1):
            tx = xh + n * rx[n + 1]
            ty = yh - n * ry[n + 1]
            tn = tx * tx + ty * ty
            rx[n] = 0.5 * tx / tn
            ry[n] = 0.5 * ty / tn
        sx = 0.0
        sy = 0.0
        for n in range(nc, 0, -1):
            saux = sx + xl
            sx = rx[n] * saux - ry[n] * sy
            sy = rx[n] * sy + ry[n] * saux
            xl = h * xl
        wx = CC * sx
        wy = CC * sy
    else:
        xh = y
        yh = x
        r1 = 0.0
        r2 = 0.0
        for n in range(9, 0, -1):
            tx = xh + n * r1
            ty = yh - n * r2
            tn = tx * tx + ty * ty
            r1 = 0.5 * tx / tn
            r2 = 0.5 * ty / tn
        wx = CC * r1
        wy = CC * r2
    if y == 0.0:
        wx = math.exp(-x**2)
    if yy < 0.0:
        wx = 2.0 * math.exp(y * y - x * x) * math.cos(2.0 * x * y) - wx
        wy = -2.0 * math.exp(y * y - x * x) * math.sin(2.0 * x * y) - wy
        if xx > 0.0:
            wy = -wy
    else:
        if xx < 0.0:
            wy = -wy
    return wx, wy


@njit(**NJ)
def _limitnan(a, lim):
    if np.isnan(a):
        return lim
    if abs(a) <= lim:
        return a
    return lim if a > 0.0 else -lim


@njit(**NJ)
def bb_kick(S, n, nsli, gamma, gambet, ce, xangle, cod, benv5, vcen, zc, sigz, tilt, icross, deform, codf,
            bstrl, cnbs, cubs, dszs, xitab, switab, rs, info, photon, prec):
    """beambeam の本体。S[:, (x, px, y, py, z, g, dv, ..., dvfs)] (px, py は px/(1+δ))"""
    sinx, cosx, tanx = xangle[1], xangle[2], xangle[3]
    sint, cost = xangle[4], xangle[5]            # K64 のまま (sint = 傾き角、cost = sin(傾き角))
    for i in range(n):
        x = S[i, 0]
        y = S[i, 2]
        z = S[i, 4]
        g = S[i, 5]
        pn = 1.0 + g
        px = S[i, 1] * pn
        py = S[i, 3] * pn
        if icross:
            pxy2 = px**2 + py**2
            pn = 1.0 + g
            pz = math.sqrt(pn**2 - pxy2)
            h1 = pxy2 / (pn + pz)
            g = g - tanx * px + tanx * tanx * h1
            px = (px - tanx * h1) / cosx
            py = py / cosx
            pxy2 = px**2 + py**2
            pz = 1.0 / math.sqrt((1.0 + g)**2 - pxy2)
            xx = x
            x = tanx * z + (1.0 + sinx * px * pz) * xx
            y = y + sinx * py * xx * pz
            z = z / cosx - tanx / cosx * h1 * xx * pz
        if codf:
            x = x - cod[0]
            px = px + cod[1]
            y = y - cod[2]
            py = py + cod[3]
            z = z + cod[4]
        if tilt:
            xt = cost * x + sint * y
            y = cost * y - sint * x
            x = xt
            xt = cost * px + sint * py
            py = cost * py - sint * px
            px = xt
        for s in range(nsli):
            zcs = zc[s] * sigz
            pxy2 = px**2 + py**2
            sz = (z - zcs) / 2.0
            x = x + px * sz
            y = y + py * sz
            g = g - pxy2 * 0.25
            pns = 1.0 + g
            if deform:
                x = x - zcs * vcen[0]
                px = px - zcs * vcen[1]
                y = y - zcs * vcen[2]
                py = py - zcs * vcen[3]
            sig11x = benv5[0, 0] + 2.0 * benv5[1, 0] * sz + benv5[1, 1] * sz**2
            sig12x = benv5[1, 0] + benv5[1, 1] * sz
            sig11y = benv5[2, 2] + 2.0 * benv5[3, 2] * sz + benv5[3, 3] * sz**2
            sig12y = benv5[3, 2] + benv5[3, 3] * sz
            asp = math.sqrt(sig11y / sig11x)
            xx2 = x * x / sig11x
            yy2 = y * y / sig11y
            if xx2 < 36.0 and yy2 < 36.0:
                info[2] += math.exp(-xx2 * 0.5) * math.exp(-yy2 * 0.5) / (2.0 * math.pi) / math.sqrt(sig11x * sig11y)
            cn = ce / nsli
            hi = -x**2 / 2.0 / sig11x - y**2 / 2.0 / sig11y
            if abs((sig11x - sig11y) / sig11x) <= 1.0e-5:
                c1 = cn
                a2 = 2.0 * (1.0 - math.exp(hi)) / (x**2 + y**2)
                fx0 = c1 * a2 * x
                fy0 = c1 * a2 * y
                gxy = c1 * sig12x / sig11x * math.exp(hi)
            elif sig11x > sig11y:
                sq = math.sqrt(2.0 * (sig11x - sig11y))
                c1 = ce * 2.0 * SQRPI / sq / nsli
                d1x = abs(x) / sq
                d1y = abs(y) / sq
                w1y, w1x = cerrf(d1x, d1y)
                w2y, w2x = cerrf(d1x * asp, d1y / asp)
                a2 = math.exp(hi)
                fx0 = c1 * (w1x - a2 * w2x)
                fy0 = c1 * (w1y - a2 * w2y)
                if x < 0.0:
                    fx0 = -fx0
                if y < 0.0:
                    fy0 = -fy0
                cdu = -1.0 / (sig11x - sig11y) / 2.0
                fxy = x * fx0 + y * fy0
                dux = cdu * (fxy + 2.0 * cn * (a2 * asp - 1.0))
                duy = -cdu * (fxy + 2.0 * cn * (a2 / asp - 1.0))
                gxy = sig12x * dux + sig12y * duy
            else:
                sq = math.sqrt(2.0 * (sig11y - sig11x))
                c1 = cn * 2.0 * SQRPI / sq
                d1x = abs(x) / sq
                d1y = abs(y) / sq
                w1x, w1y = cerrf(d1y, d1x)
                w2x, w2y = cerrf(d1y / asp, d1x * asp)
                a2 = math.exp(hi)
                fx0 = c1 * (w1x - a2 * w2x)
                fy0 = c1 * (w1y - a2 * w2y)
                if x < 0.0:
                    fx0 = -fx0
                if y < 0.0:
                    fy0 = -fy0
                cdu = -1.0 / (sig11x - sig11y) / 2.0
                fxy = x * fx0 + y * fy0
                dux = cdu * (fxy + 2.0 * cn * (a2 * asp - 1.0))
                duy = -cdu * (fxy + 2.0 * cn * (a2 / asp - 1.0))
                gxy = sig12x * dux + sig12y * duy
            px = px - fx0
            py = py - fy0
            g = g - gxy
            if bstrl:
                dsz = dszs[s]
                rhoi2 = (fx0 * fx0 + fy0 * fy0) / dsz / dsz
                rhoi = math.sqrt(rhoi2)
                rnbs = cnbs * rhoi * dsz
                if rnbs > 0.5:                       # K64: 'Nbs is too large' で stop
                    info[0] = rnbs
                    return
                xr = tran(rs)
                if xr < rnbs:
                    xr = xr / rnbs
                    jl = 0
                    ju = 1001
                    while ju - jl > 1:
                        jm = (ju + jl) // 2
                        if (switab[1000] > switab[1]) == (xr > switab[jm]):
                            jl = jm
                        else:
                            ju = jm
                    j = jl
                    xij = xitab[j] + (xitab[j + 1] - xitab[j]) * (xr - switab[j]) / (switab[j + 1] - switab[j])
                    delg = cubs * rhoi * math.exp(xij)
                    if photon:
                        kk = int(info[3])
                        prec[kk, 0] = i + 1
                        prec[kk, 1] = x
                        prec[kk, 2] = px / pns
                        prec[kk, 3] = y
                        prec[kk, 4] = py / pns
                        prec[kk, 5] = delg
                        prec[kk, 6] = sz
                        info[3] += 1.0
                        rs[0] = lcg(lcg(rs[0]))      # radangle1 の tgauss 2 回 (Box-Muller の 1 組)
                    g = g - delg
                    info[1] += 1.0
            if deform:
                x = x + zcs * vcen[0]
                px = px + zcs * vcen[1]
                y = y + zcs * vcen[2]
                py = py + zcs * vcen[3]
            pxy2 = px**2 + py**2
            x = x - px * sz
            y = y - py * sz
            g = g + pxy2 * 0.25
        if tilt:
            xt = cost * x - sint * y
            y = cost * y + sint * x
            x = xt
            xt = cost * px - sint * py
            py = cost * py + sint * px
            px = xt
        if codf:
            x = x + cod[0]
            px = px - cod[1]
            y = y + cod[2]
            py = py - cod[3]
            z = z - cod[4]
        if icross:
            pxy2 = px**2 + py**2
            pn = 1.0 + g
            pz = math.sqrt(pn**2 - pxy2)
            h1 = pxy2 / (pn + pz)
            pz = 1.0 / pz
            x = (x - sinx * z) / (1.0 + (sinx * (px + sinx * h1)) * pz)
            y = y - sinx * py * pz * x
            z = z * cosx + sinx * cosx * h1 * pz * x
            h1 = h1 * cosx**2
            xx = px
            px = xx * cosx + tanx * h1
            py = py * cosx
            g = g + sinx * xx
        x = _limitnan(x, 1.0e4)
        px = _limitnan(px, 1.0e4)
        y = _limitnan(y, 1.0e4)
        py = _limitnan(py, 1.0e4)
        z = _limitnan(z, 1.0e10)
        g = _limitnan(g, 1.0e4)
        pn = 1.0 + g
        S[i, 0] = x
        S[i, 1] = px / pn
        S[i, 2] = y
        S[i, 3] = py / pn
        S[i, 4] = z
        S[i, 5] = g
        h1 = math.sqrt(1.0 + (pn * gambet)**2)
        S[i, 6] = -g * (1.0 + pn) / h1 / (h1 + pn * gamma) + S[i, 10]


PHOTON_SINK = None        # TrackParticles (PHOTONS) が設定するリスト。PhotonList の行を足す
PHOTON_OUT = None         # K64 の write(*,*) x,pxi,y,pyi,delg に当たる出力関数 f(str)


def _fort(v):
    """gfortran の list-directed real*8 (1PE24.16E3 に近い)"""
    t = '%.16E' % v
    m, e = t.split('E')
    return m + 'E' + e[0] + '%03d' % int(e[1:])


def _photon_rows(rec, bb):
    """tphotonconv (rho=0: 直線) と tphotonlist の 1 光子ぶん。rec: (k, x, pxi, y, pyi, dp, ds)"""
    l, gv = bb['photon']
    k, xi, pxi, yi, pyi, dp, ds = rec
    pzi = 1.0 + pxy2dpz(pxi, pyi)
    x1, x2, x3 = gv[0, 0], gv[1, 0], gv[2, 0]
    y1, y2, y3 = gv[0, 1], gv[1, 1], gv[2, 1]
    z1, z2, z3 = gv[0, 2], gv[1, 2], gv[2, 2]
    gx = gv[0, 3] + xi * x1 + yi * y1 + ds * z1
    gy = gv[1, 3] + xi * x2 + yi * y2 + ds * z2
    gz = gv[2, 3] + xi * x3 + yi * y3 + ds * z3
    dpx = dp * (pzi * z1 + pxi * x1 + pyi * y1)
    dpy = dp * (pzi * z2 + pxi * x2 + pyi * y2)
    dpz = dp * (pzi * z3 + pxi * x3 + pyi * y3)
    chi = 0.0 if x3 == 0.0 else 2.0 * math.atan2(x3, -y3)       # tsetphotongeo
    xi3, xi1, xi2 = 1.0, 0.0, 0.0
    xi3a = math.cos(chi) * xi3 + math.sin(chi) * xi1
    xi1 = -math.sin(chi) * xi3 + math.cos(chi) * xi1
    xi3 = xi3a
    dn = math.sqrt(dpx**2 + dpy**2 + dpz**2)
    return [dn * bb['momentum'], float(gx), float(gy), float(gz), dpx / dn, dpy / dn, dpz / dn, xi1, xi2, xi3,
            float(k), float(l)]


def apply(S, n, bb, rs=None):
    """bb = setup() の戻り値。rs: 乱数の状態 (int64 の 1 要素の配列、BSTRL>0 のとき)。
    戻り値は (ルミノシティ Σ dLum / スライス数 / 粒子数、放射した光子の数)"""
    ns = bb['nslice']
    if ns <= 0 or n <= 0:
        return 0.0, 0.0
    xa = bb['xangle']
    tilt_on = xa[4] != 0.0
    if compat.v2() and tilt_on:
        # B34 の修正: 傾きの回転は (cos t, sin t) を使う (bbinit は xangle(6)=sin、(7)=cos を書くが、K64 の beambeam は
        # sint=xangle(5)=t、cost=xangle(6)=sin t を読む)。bb_kick は xangle[4] を sint、xangle[5] を cost として読む。
        xa = xa.copy()
        xa[4] = bb['xangle'][5]
        xa[5] = bb['xangle'][6]
    cod = bb['cod']
    info = np.zeros(4)
    photon = bool(bb['bstrl']) and bb.get('photon') is not None and PHOTON_SINK is not None
    prec = np.zeros((n * ns if photon else 1, 7))
    if bb['bstrl']:
        xit, swit = besseltab()
        if rs is None:
            rs = np.array([1], dtype=np.int64)
        dsz = bb['dsz']
    else:
        xit = swit = np.zeros(1)
        dsz = np.zeros(1)
        rs = np.zeros(1, dtype=np.int64)
    bb_kick(S, n, ns, bb['gamma'], bb['gambet'], bb['ce'], xa, cod, bb['benv5'], bb['vcen'], bb['zc'], bb['sigz'],
            tilt_on, xa[0] != 0.0, bool(np.any(bb['vcen'][:4] != 0.0)), bool(np.any(cod[:4] != 0.0)),
            bool(bb['bstrl']), bb.get('cnbs', 0.0), bb.get('cubs', 0.0), dsz, xit, swit, rs, info, photon, prec)
    if info[0] > 0.0:
        raise ValueError(f'BEAMBEAM: Nbs is too large ({info[0]}) (K64 は stop する)')
    if photon:
        for r in prec[:int(info[3])]:
            if PHOTON_OUT is not None:
                PHOTON_OUT(''.join(' ' + _fort(v).rjust(25) for v in (r[1], r[2], r[3], r[4], r[5])) + chr(10))
            PHOTON_SINK.append(_photon_rows(r, bb))
    with np.errstate(all='ignore'):
        lum = float(info[2] / ns / n)
    if not math.isfinite(lum):
        lum = 0.0
    if LUM_CB is not None:
        LUM_CB(lum)                     # K64: symbol LuminosityMonitor に値を入れ、LuminosityMonitor[lum] を呼ぶ
    return lum, info[1]
