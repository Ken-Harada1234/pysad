"""
track6n.py — track6.py (SAD tturn の直訳) の numba 版: 粒子配列に対するカーネル

track6.py と同じ演算を同じ順序で行う (結果は同一になるはず)。要素は compile_lattice() で
(種類コード, パラメータ) の表にし、1 周は turn() が要素ごとに njit カーネルを呼ぶ。
対応: DRIFT, QUAD (tquad: 回転・ttfrin4・tqlfr・tsolqu/tsolqur), SEXT/OCT/DECA/DODECA (tthin),
BEND (tbendcore/tbendr, K1≠0 は tbendi), CAVI (tcav), ANGLE=0 BEND (tsteer: L=0 薄いキック / ドリフト),
MARK/MONI。放射 (tradk, tradkn) と FLUC (tradkfn: tdusrn/drndsr の乱数は SAD の LCG を int64 で)。
MULT とソレノイドは未対応 (track_particles_fast は track6 にフォールバックする)。
粒子の損失判定 (tapert) と配列の詰め替え (SAD と同じ: 落ちた粒子を末尾の生きている粒子と入れ替える) も kernel。
"""
from __future__ import annotations

import math
from typing import List

import numpy as np
from numba import njit

from . import emit6d as _e
from . import lattice as _lat
from . import track6 as _t6
from .elements import Drift, Bend, Quadrupole

CVELOC = 299792458.0
GMIN = -0.999999
CAVE = 8.0 / 15.0 / math.sqrt(3.0)
TT_ALPHA = 2.88675134594812882e-1
TT_BETA = 2.23290993692602255e-2
TT_GAMMA = 9.43738783765593145e-4
TT_ALPHA1 = 0.5 - TT_ALPHA
RBH = 0.5 + 1.0 / math.sqrt(12.0)
RBL = 1.0 / 6.0 / RBH
FACT = np.array([1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0])

NJ = dict(cache=True, error_model='numpy', nogil=True)

# 要素の種類
K_DRIFT, K_QUAD, K_THIN, K_BEND, K_CAV, K_STEER0, K_NONE = 0, 1, 2, 3, 4, 5, 6
K_STEER = 14           # ANGLE=0 の BEND (L>0, K0≠0): tsteer (K_MULT.. = 10..13 は track6n_ext)


# gfortran (-O1, x86-64 の maxsd/minsd) の MAX/MIN: a > b ? a : b。NaN が混じると後ろの引数を返す
# (Python の max/min は前の引数を返すので、発散した粒子の挙動が変わる)。NAN_MODE で切り替え
NAN_MODE = 1


@njit(**NJ)
def fmax(a, b):
    if NAN_MODE == 1:
        return a if a > b else b
    return b if b > a else a


@njit(**NJ)
def fmin(a, b):
    if NAN_MODE == 1:
        return a if a < b else b
    return b if b < a else a


# ============================================================================ mathfun
@njit(**NJ)
def gfloor(v):
    """max(gmin, g+dg) (tradkfn/tradkn)。SAD (resync 検証、scratch_rng2) では g+dg が NaN のとき gmin が返る
    (maxsd のオペランド順が fmax(GMIN, v) = 'b if not a > b' と逆)。NaN の粒子が g=NaN のままだと次の要素の tdusrn が
    anp=NaN で 10000 回ループして乱数を 20001 個消費してしまう"""
    return v if v > GMIN else GMIN


@njit(**NJ)
def sqrt1(x):
    """mathfun sqrt1: sqrt(1+x)-1 (SAD と同じ級数の分岐と max(xmin, 1+x) のガード。
    1+x<0 でも NaN にならず有限値を返す = 落ちる粒子の挙動が SAD と同じになる)"""
    ax = abs(x)
    if ax < 8e-6:
        return x * (0.5 + x * (-0.125 + x * 0.0625))
    if ax < 3.4e-3:
        return x * (0.5 + x * (-0.125 + x * (0.0625 + x * (-.0390625 + x * (.02734375 - .0205078125 * x)))))
    return x / (1.0 + math.sqrt(fmax(1e-200, 1.0 + x)))


@njit(**NJ)
def pxy2dpz(px, py):
    """mathfun pxy2dpz: sqrt(1-px^2-py^2)-1 (同上)"""
    x = px * px + py * py
    if x < 8e-6:
        return -x * (0.5 + x * (0.125 + x * 0.0625))
    return -x / (1.0 + math.sqrt(fmax(1e-200, 1.0 - x)))


@njit(**NJ)
def p2h(p):
    if p > 1e3:
        p2 = 1.0 / (p * p)
        return p * (1.0 + p2 * (0.5 - p2 * 0.125))
    return math.hypot(1.0, p)


@njit(**NJ)
def asinz(x):
    """mathfun asinz: |x|<0.0032 は級数、それ以外は asin(min(1, max(-1, x))) (tradkfn, tbendcore が使う)"""
    if abs(x) < 0.0032:
        x2 = x * x
        return x + x * x2 * (1.0 / 6.0 + x2 * 3.0 / 40.0)
    if x != x:
        # SAD (gfortran -Ofast) では max(-1, NaN) が -1 になり asin(-1) = -π/2 が返る (scratch_loss/ftest で確認)。
        # 発散した粒子の tradkfn で theta<0 → anp<0 → 光子なし (乱数 1 個) となり、落ちるまでの乱数消費が SAD と揃う
        return -0.5 * math.pi
    return math.asin(fmin(1.0, fmax(-1.0, x)))


asinz_e = asinz
asinz_l = asinz


@njit(**NJ)
def norm2_3(a, b, c):
    """gfortran の norm2([a,b,c]) (インライン展開: scale/ssq 法)。tradkfn の ppa に使う。
    成分に ±inf が 2 つ以上あると (inf/inf)^2 で NaN になる (発散した粒子で SAD と同じ結果にするため必要:
    NaN → asinz = -π/2 → anp<0 → 光子なし。math.sqrt(a²+b²+c²) だと inf → +π/2 → 光子が出て乱数消費がずれる)"""
    scale = 0.0
    ssq = 1.0
    for x in (a, b, c):
        if x != 0.0:
            ax = abs(x)
            if scale < ax:
                v = scale / ax
                ssq = 1.0 + ssq * v * v
                scale = ax
            else:
                v = ax / scale
                ssq += v * v
    return scale * math.sqrt(ssq)


@njit(**NJ)
def dvfun(p0, h0, dp):
    """_Phys.dv (dvfs=0。dvfs は呼び出し側で足す)"""
    p1 = p0 * (1.0 + dp)
    h1 = p2h(p1)
    return -(p1 + p0) / h1 / (p0 * h1 + p1 * h0) * dp


# ============================================================================ 乱数 (sadrandom の int64 版)
@njit(**NJ)
def lcg(seed):
    s = (48828125 * seed) % 4294967296
    if s >= 2147483648:
        s -= 4294967296
    return s


@njit(**NJ)
def tran(rs):
    rs[0] = lcg(rs[0])
    return (rs[0] + 2147483648.0) * (1.0 / 4294967296.0)


@njit(**NJ)
def drndsr(rs):
    p = tran(rs)
    if p < 0.658:
        p2 = p * p
        return (((0.41841680 * p2 + 0.14175015) * p2 + 0.30528148) * p2 + 0.53520052) * p2 * p
    q = 1.0 - p
    if q > 0.0297:
        return (((-0.32809490 * q + 0.20653610) * q + 0.011920005)
                / (((q + 0.88765651) * q + 0.19270658) * q + 0.0033141692))
    t = -math.log(q)
    return (674.99117 * t + 148.32871) / ((t - 225.45670) * t - 692.23986) + t


@njit(**NJ)
def tdusrn(rs, anp):
    p0 = math.exp(-anp)
    pn = p0
    p = tran(rs)
    r1 = 0.0
    r2 = 0.0
    an = 0.0
    dph = 0.0
    x = 1.0
    for _ in range(10000):
        if p <= pn:
            break
        an += 1.0
        x = x * anp / an
        pn = pn + x * p0
        dph1 = drndsr(rs)
        dph += dph1
        r = tran(rs)
        r1 += r * dph1
        r2 += r * r * dph1
        if pn > 0.9999:
            break
    if r1 != 0.0:
        r2 = r2 / r1 * 0.5
        r1 = r1 / dph
    return dph, r1, r2, an


# ============================================================================ 放射 (tradk)
@njit(**NJ)
def tradk_core1(i, S, al, phi0, fluc, rs, p0, anrad, cuc):
    """S: (n, 11) [x, px, y, py, z, g, dv, pxr0, pyr0, zr0, dvfs]。基準点は更新しない"""
    if al == 0.0:
        return
    x = S[i, 0]; px = S[i, 1]; y = S[i, 2]; py = S[i, 3]; z = S[i, 4]; g = S[i, 5]
    pxr0 = S[i, 7]; pyr0 = S[i, 8]; zr0 = S[i, 9]
    cphi0 = math.cos(phi0)
    sphi0 = math.sin(phi0)
    dpz0 = pxy2dpz(pxr0, pyr0)
    px0 = cphi0 * pxr0 + sphi0 * (1.0 + dpz0)
    dpz0 = cphi0 * dpz0 - sphi0 * pxr0
    dpx = px - px0
    dpy = py - pyr0
    dpz = pxy2dpz(px, py)
    if fluc:
        dpz0 = pxy2dpz(px0, pyr0)
    ppx = py * dpz0 - dpz * pyr0 + dpy
    ppy = dpz * px0 - px * dpz0 - dpx
    ppz = px * pyr0 - py * px0
    ppa = norm2_3(ppx, ppy, ppz)
    theta = asinz(ppa)
    pr = 1.0 + g
    pp = p0 * pr
    h1 = p2h(pp)
    al1 = al - z + zr0
    anp = anrad * h1 * theta
    if fluc:
        dph, r1, r2, an = tdusrn(rs, anp)
        if an == 0.0:
            return
        rho = al1 / theta
        uc = cuc * h1 * h1 * h1 / p0 / rho
        de = dph * uc
        dg = -de / (1.0 + 2.0 * de)
        g = gfloor(g + dg)
        ddpx = -r1 * dpx * dg
        ddpy = -r1 * dpy * dg
        x += r2 * ddpx * al1
        y += r2 * ddpy * al1
    else:
        uc = cuc * h1 * h1 * h1 / p0 * theta / al1
        de = CAVE * anp * uc
        dg = -de / (1.0 + 2.0 * de)
        g = gfloor(g + dg)
        ddpx = -0.5 * dpx * dg
        ddpy = -0.5 * dpy * dg
        x += ddpx * al1 / 3.0
        y += ddpy * al1 / 3.0
    px += ddpx
    py += ddpy
    pr = 1.0 + g
    p2 = p0 * pr
    h2 = p2h(p2)
    S[i, 6] = -g * (1.0 + pr) / h2 / (h2 + p2) + S[i, 10]
    z = z * p2 / h2 * h1 / pp
    S[i, 0] = x; S[i, 1] = px; S[i, 2] = y; S[i, 3] = py; S[i, 4] = z; S[i, 5] = g


@njit(**NJ)
def setr0_1(i, S):
    S[i, 7] = S[i, 1]
    S[i, 8] = S[i, 3]
    S[i, 9] = S[i, 4]


@njit(**NJ)
def tradk_all(S, n, al, phi0, fluc, rs, p0, anrad, cuc):
    """tradk: 全粒子 (SAD と同じ順で乱数を引く) + setr0"""
    for i in range(n):
        tradk_core1(i, S, al, phi0, fluc, rs, p0, anrad, cuc)
        setr0_1(i, S)


# ============================================================================ ドリフト
@njit(**NJ)
def drift1(i, S, al):
    px = S[i, 1]; py = S[i, 3]
    dpz = pxy2dpz(px, py)
    al1 = al / (1.0 + dpz)
    S[i, 0] += px * al1
    S[i, 2] += py * al1
    S[i, 4] += dpz * al1 - S[i, 6] * al


@njit(**NJ)
def k_drift(S, n, al):
    for i in range(n):
        drift1(i, S, al)


# ============================================================================ 四極
@njit(**NJ)
def ttfrin4_1(i, S, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    aki = ak / al / 4.0 / (1.0 + S[i, 5])
    xi = S[i, 0]; yi = S[i, 2]
    a = aki * xi * xi
    b = aki * yi * yi
    ab = a - b
    t = ab * ab / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    S[i, 0] = xi + dx1
    S[i, 2] = yi + dy1
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h * h
    px = S[i, 1]; py = S[i, 3]
    px1 = (px * yy + py * h) / f
    py1 = (py * xx - px * h) / f
    S[i, 4] = S[i, 4] - px1 * (dx1 + t * xi) - py1 * (dy1 + t * yi)
    S[i, 1] = px1
    S[i, 3] = py1


@njit(**NJ)
def tqlfr1(i, S, a1, b1):
    pr = 1.0 + S[i, 5]
    a = a1 / pr
    ea = math.exp(a)
    b = b1 / pr
    x = S[i, 0]; px = S[i, 1]; y = S[i, 2]; py = S[i, 3]
    pxf = px / ea
    pyf = py * ea
    S[i, 4] = (S[i, 4] - (a * x + b * (1.0 + 0.5 * a) * pxf) * px
               + (a * y + b * (1.0 - 0.5 * a) * pyf) * py)
    S[i, 0] = ea * x + b * px
    S[i, 2] = y / ea - b * py
    S[i, 1] = pxf
    S[i, 3] = pyf


@njit(**NJ)
def rot1(i, S, ct, st):
    # tsolrot は chi3 /= 0 のときだけ回す。無条件に掛けると x=inf, y=inf の粒子で 0*inf=NaN になり SAD (x=inf のまま) と変わる
    if st == 0.0 and ct == 1.0:
        return
    x = S[i, 0]; px = S[i, 1]; y = S[i, 2]; py = S[i, 3]
    S[i, 0] = ct * x - st * y
    S[i, 2] = st * x + ct * y
    S[i, 1] = ct * px - st * py
    S[i, 3] = st * px + ct * py


# mathfun (tfloor.f): 級数に切り替えるしきい値
_EPSIEEE = 5e-17
XTHS1 = (_EPSIEEE * 362880.0 / 6.0) ** (1.0 / 6.0)
XTHS2 = (_EPSIEEE * 6227020800.0 / 6e9) ** (1.0 / 10.0)
XTHC1 = (_EPSIEEE * 40320.0 / 2.0) ** (1.0 / 6.0)
XTHC2 = (_EPSIEEE * 479001600.0 / 2.0) ** (1.0 / 10.0)
XTHCS1 = min(XTHS1, XTHC1)
XTHCS2 = min(XTHS2, XTHC2)
# gfortran -Ofast は定数での割り算を逆数の掛け算にする (-freciprocal-math): x/5040 は x*(1/5040) と丸めが違う
R5040 = 1.0 / 5040.0
R720 = 1.0 / 720.0
R39916800 = 1.0 / 39916800.0
R3628800 = 1.0 / 3628800.0


@njit(**NJ)
def xsincos(x):
    """mathfun xsincos の直訳: 小さい |x| は級数 (s=x-xs, dc, c=1+dc)、それ以外は sin/cos。
    戻り (s, xs, c, dc)。ライブラリ関数だけで書くと SAD と下位ビットが変わり、カオス的な粒子で増幅される"""
    if abs(x) < XTHCS1:
        x2 = x * x
        xs = x * x2 * (1.0 / 6.0 - x2 * (1.0 / 120.0 - x2 * R5040))
        dc = -x2 * (0.5 - x2 * (1.0 / 24.0 - x2 * R720))
        return x - xs, xs, 1.0 + dc, dc
    if abs(x) < XTHCS2:
        x2 = x * x
        xs = x * x2 * (1.0 / 6.0 - x2 * (1.0 / 120.0 - x2 * (1.0 / 5040.0 - x2 * (1.0 / 362880.0 - x2 * R39916800))))
        dc = -x2 * (0.5 - x2 * (1.0 / 24.0 - x2 * (1.0 / 720.0 - x2 * (1.0 / 40320.0 - x2 * R3628800))))
        return x - xs, xs, 1.0 + dc, dc
    s = math.sin(x)
    c = math.cos(x)
    if c > 0.0:
        dc = -s * s / (1.0 + c)
    else:
        dc = c - 1.0
    return s, x - s, c, dc


@njit(**NJ)
def xsincosh(x):
    """mathfun xsincosh の直訳。戻り (sh, xsh, ch, dch)"""
    if abs(x) < XTHCS1:
        x2 = x * x
        xsh = -x * x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 * R5040))
        dch = x2 * (0.5 + x2 * (1.0 / 24.0 + x2 * R720))
        return x - xsh, xsh, 1.0 + dch, dch
    if abs(x) < XTHCS2:
        x2 = x * x
        xsh = -x * x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 * (1.0 / 5040.0 + x2 * (1.0 / 362880.0 + x2 * R39916800))))
        dch = x2 * (0.5 + x2 * (1.0 / 24.0 + x2 * (1.0 / 720.0 + x2 * (1.0 / 40320.0 + x2 * R3628800))))
        return x - xsh, xsh, 1.0 + dch, dch
    sh = math.sinh(x)
    ch = math.cosh(x)
    return sh, x - sh, ch, sh * sh / (1.0 + ch)


@njit(**NJ)
def tzsetparam0(dp, aln0, akk):
    pr = 1.0 + dp
    akkp = akk / pr
    w1 = math.sqrt(akkp)
    phi1 = fmin(100.0, fmax(-100.0, aln0 * w1))
    s1, xs1, c1, dc1 = xsincos(phi1)
    sh2, xsh2, ch2, dch2 = xsincosh(phi1)
    return w1, s1, xs1, dc1, sh2, xsh2, dch2


@njit(**NJ)
def tsolqu1(i, S, al, ak, eps0, qum):
    """tsolqu / tsolqum (bz=0, ak0=0, ak>0)"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    ndiv = 1 + int(abs(al * ak / eps))
    aln0 = al / ndiv
    akk = ak / al
    w1, s1, xs1, dc1, sh2, xsh2, dch2 = tzsetparam0(S[i, 5], aln0, akk)
    x = S[i, 0]; px = S[i, 1]; y = S[i, 2]; py = S[i, 3]; z = S[i, 4]
    dv = S[i, 6]
    for n in range(ndiv):
        pxi = px; pyi = py
        ap = pxi * pxi + pyi * pyi
        dpz = sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * aln0
        if n == 0:
            r *= 0.5
        x += pxi * r
        y += pyi * r
        xi = x; yi = y
        u1 = xi * dc1 + pxi * s1 / w1
        u2 = -xi * w1 * s1 + pxi * dc1
        v1 = yi * dch2 + pyi * sh2 / w1
        v2 = yi * w1 * sh2 + pyi * dch2
        x += u1
        px = pxi + u2
        y += v1
        py = pyi + v2
        if qum:
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi * xi * xs1 - yi * yi * xsh2) + (ap - dv) * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2))
        else:
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi * xi * xs1 - yi * yi * xsh2) + ap * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2)
                 - dv * aln0)
    ap = px * px + py * py
    dpz = sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0 * 0.5
    x += px * r
    y += py * r
    z -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
    S[i, 0] = x; S[i, 1] = px; S[i, 2] = y; S[i, 3] = py; S[i, 4] = z


@njit(**NJ)
def ndivrad(phi, ak, bz, eps, h0):
    eps1 = 1.0 if eps <= 0.0 else eps
    aka = (abs(phi) + math.hypot(ak, bz) * 0.02) / eps1
    return int(1.0 + h0 / 100.0 * aka)


@njit(**NJ)
def tsolqur_all(S, n, al, ak, eps0, fluc, rs, p0, h0, anrad, cuc):
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    ndiv = 1 + int(abs(al) * ak / eps)
    ndiv = fmin(1000, fmax(ndiv, ndivrad(0.0, ak, 0.0, eps0, h0)))
    aln = al / ndiv
    akk = ak / al
    for m in range(ndiv):
        alr = aln * 0.5 if m == 0 else aln
        for i in range(n):
            w1, s1, xs1, dc1, sh2, xsh2, dch2 = tzsetparam0(S[i, 5], aln, akk)
            px = S[i, 1]; py = S[i, 3]
            ap = px * px + py * py
            dpz = sqrt1(-ap)
            r = -dpz / (1.0 + dpz) * alr
            S[i, 0] += px * r
            S[i, 2] += py * r
            S[i, 4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
            xi = S[i, 0]; yi = S[i, 2]
            b = px; d = py
            u1 = xi * dc1 + b * s1 / w1
            u2 = -xi * w1 * s1 + b * dc1
            v1 = yi * dch2 + d * sh2 / w1
            v2 = yi * w1 * sh2 + d * dch2
            S[i, 0] += u1
            S[i, 1] += u2
            S[i, 2] += v1
            S[i, 3] += v2
            S[i, 4] = (S[i, 4] - 0.25 * (w1 * (xi * xi * xs1 - yi * yi * xsh2) + (b * b + d * d) * aln
                                         + u1 * (u2 + b) + xi * b * dc1 + v1 * (v2 + d) + yi * d * dch2)
                       - S[i, 6] * aln)
        tradk_all(S, n, aln, 0.0, fluc, rs, p0, anrad, cuc)
    for i in range(n):
        px = S[i, 1]; py = S[i, 3]
        ap = px * px + py * py
        dpz = sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * aln * 0.5
        S[i, 0] += px * r
        S[i, 2] += py * r
        S[i, 4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r


@njit(**NJ)
def k_quad(S, n, L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0, krad,
           fluc, rs, p0, h0, anrad, cuc):
    if K1L == 0.0:
        for i in range(n):
            drift1(i, S, L)
        return
    ak = abs(K1L)
    for i in range(n):
        S[i, 0] -= dx
        S[i, 2] -= dy
        rot1(i, S, ct, st)
        if krad:
            setr0_1(i, S)
        if fringe and mfring > -4 and mfring != 2:
            ttfrin4_1(i, S, ak, L)
        if mfring == 1 or mfring == 3:
            tqlfr1(i, S, f1in, f2in)
    if krad:
        if f1in != 0.0:
            tradk_all(S, n, f1in, 0.0, fluc, rs, p0, anrad, cuc)
        for i in range(n):
            setr0_1(i, S)
        tsolqur_all(S, n, L, ak, eps0, fluc, rs, p0, h0, anrad, cuc)
    else:
        for i in range(n):
            tsolqu1(i, S, L, ak, eps0, False)
    for i in range(n):
        if mfring == 2 or mfring == 3:
            tqlfr1(i, S, -f1out, f2out)
        if fringe and mfring > -4 and mfring != 1:
            ttfrin4_1(i, S, -ak, L)
    if krad and f1out != 0.0:
        tradk_all(S, n, f1out, 0.0, fluc, rs, p0, anrad, cuc)
    for i in range(n):
        rot1(i, S, ct, -st)
        S[i, 0] += dx
        S[i, 2] += dy


# ============================================================================ 薄い多極 (tthin)
@njit(**NJ)
def tent1(i, S, dx, dy, theta, cost, sint):
    if dx != 0.0:
        S[i, 0] -= dx
    if dy != 0.0:
        S[i, 2] -= dy
    if theta != 0.0:
        if abs(cost) >= abs(sint):
            S[i, 0] = cost * S[i, 0] - sint * S[i, 2]
            S[i, 2] = (sint * S[i, 0] + S[i, 2]) / cost
            S[i, 1] = cost * S[i, 1] - sint * S[i, 3]
            S[i, 3] = (sint * S[i, 1] + S[i, 3]) / cost
        else:
            xi = S[i, 0]
            S[i, 0] = xi * cost - S[i, 2] * sint
            S[i, 2] = xi * sint + S[i, 2] * cost
            pxi = S[i, 1]
            S[i, 1] = pxi * cost - S[i, 3] * sint
            S[i, 3] = pxi * sint + S[i, 3] * cost


@njit(**NJ)
def texit1(i, S, dx, dy, theta, cost, sint):
    if theta != 0.0:
        if abs(cost) >= abs(sint):
            S[i, 0] = cost * S[i, 0] + sint * S[i, 2]
            S[i, 2] = (-sint * S[i, 0] + S[i, 2]) / cost
            S[i, 1] = cost * S[i, 1] + sint * S[i, 3]
            S[i, 3] = (-sint * S[i, 1] + S[i, 3]) / cost
        else:
            xi = S[i, 0]
            S[i, 0] = xi * cost + S[i, 2] * sint
            S[i, 2] = -xi * sint + S[i, 2] * cost
            pxi = S[i, 1]
            S[i, 1] = pxi * cost + S[i, 3] * sint
            S[i, 3] = -pxi * sint + S[i, 3] * cost
    if dy != 0.0:
        S[i, 2] += dy
    if dx != 0.0:
        S[i, 0] += dx


@njit(**NJ)
def ttfrin6_1(i, S, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    aki = ak / al / 24.0 / (1.0 + S[i, 5])
    cx = complex(S[i, 0], S[i, 2])
    cp = complex(S[i, 1], -S[i, 3])
    cz1 = cx * cx
    cz = cz1 * cx
    a = aki * cz.imag * 2.0
    ca = -aki * cx * (cz / 4.0 - cz.conjugate())
    cx1 = cx + ca
    d = 1.0 + a * a - (aki * 3.0)**2 * (cz.real * cz.real + cz.imag * cz.imag)
    cp1 = (complex(1.0, a) * cp - (aki * 3.0) * cz1 * (cx * cp).conjugate()) / d
    S[i, 4] = S[i, 4] - (ca * cp1).real
    S[i, 0] = cx1.real
    S[i, 2] = cx1.imag
    S[i, 1] = cp1.real
    S[i, 3] = -cp1.imag


@njit(**NJ)
def ttfrin_gen_1(i, S, nord, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    kord = nord // 2
    akk = ak / al / FACT[kord] / 4.0
    aki = akk / (1.0 + S[i, 5])
    cx = complex(S[i, 0], S[i, 2])
    cp = complex(S[i, 1], -S[i, 3])
    cz1 = cx**(kord - 1)
    cz = cz1 * cx
    a = aki * cz.imag * 2.0
    ca = -aki * cx * (cz / (kord + 1) - cz.conjugate())
    cx1 = cx + ca
    d = 1.0 + a * a - (aki * kord)**2 * (cz.real * cz.real + cz.imag * cz.imag)
    cp1 = (complex(1.0, a) * cp - (aki * kord) * cz1 * (cx * cp).conjugate()) / d
    S[i, 4] = S[i, 4] - (ca * cp1).real
    S[i, 0] = cx1.real
    S[i, 2] = cx1.imag
    S[i, 1] = cp1.real
    S[i, 3] = -cp1.imag


@njit(**NJ)
def ttfrin_n_1(i, S, nord, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    if nord == 4:
        ttfrin4_1(i, S, ak, al)
    elif nord == 6:
        ttfrin6_1(i, S, ak, al)
    elif nord == 2:
        aki = ak / al / (1.0 + S[i, 5])
        S[i, 0] += 0.5 * aki * S[i, 2]**2
        S[i, 3] -= aki * S[i, 1] * S[i, 2]
        S[i, 4] -= 0.5 * aki * S[i, 1] * S[i, 2]**2
    else:
        ttfrin_gen_1(i, S, nord, ak, al)


@njit(**NJ)
def k_thin(S, n, L, ak, kord, dx, dy, theta, cost, sint, fringe, krad, fluc, rs, p0, anrad, cuc):
    """tthin (SEXT kord=2 と一般 kord)。ak = 積分強度"""
    if ak == 0.0:
        for i in range(n):
            drift1(i, S, L)
        return
    ala = L * TT_ALPHA1
    alb = L * TT_ALPHA
    for i in range(n):
        tent1(i, S, dx, dy, theta, cost, sint)
        if krad:
            setr0_1(i, S)
        if fringe:
            ttfrin_n_1(i, S, 2 * kord + 2, ak, L)
        akf = ak / FACT[kord]
        if L != 0.0:
            drift1(i, S, ala)
            akf *= 0.5
        # kick 1
        aki = akf / (1.0 + S[i, 5])
        if kord == 2:
            x = S[i, 0]; y = S[i, 2]
            S[i, 1] -= aki * (x - y) * (x + y)
            S[i, 3] += 2.0 * aki * x * y
        else:
            cx = complex(S[i, 0], -S[i, 2])**kord
            S[i, 1] -= aki * cx.real
            S[i, 3] -= aki * cx.imag
        if L != 0.0:
            if kord == 1:
                f1 = TT_BETA * L
                f3 = TT_GAMMA * 2.0 * L * L
                f4 = 0.5 * f1
                f5 = f3
                f2 = 0.0
            else:
                f1 = TT_BETA * kord * L
                f2 = 2.0 * TT_GAMMA * kord * (kord - 1) * L * L
                f3 = TT_GAMMA * kord * (kord + 1) * L * L
                f4 = 0.5 * TT_BETA * L
                f5 = 2.0 * TT_GAMMA * kord * L * L
            px = S[i, 1]; py = S[i, 3]
            dpz = pxy2dpz(px, py)
            al1 = alb / (1.0 + dpz)
            xi = S[i, 0] + px * al1
            yi = S[i, 2] + py * al1
            zi = S[i, 4] + dpz * al1 - S[i, 6] * alb
            aki = akf / (1.0 + S[i, 5]) * 2.0
            r = xi * xi + yi * yi
            if kord == 1:
                rcx1 = (xi - yi) * (xi + yi)
                S[i, 1] += aki * aki * (f1 - aki * f3) * xi
                S[i, 3] += aki * aki * (f1 + aki * f3) * yi
                zi += aki * aki * (f4 * r - f5 * aki * rcx1)
            elif kord == 2:
                cx = complex(xi, -yi)**2
                rcx1 = xi * cx.real + yi * cx.imag
                S[i, 1] += aki * aki * (f1 * r * xi - aki * (f2 * xi * rcx1 + f3 * r * cx.real))
                S[i, 3] += aki * aki * (f1 * r * yi - aki * (f2 * yi * rcx1 + f3 * r * cx.imag))
                zi += aki * aki * r * (f4 * r - f5 * aki * rcx1)
            else:
                rk1 = r**(kord - 2)
                rk = r * rk1
                cx = complex(xi, -yi)**kord
                rcx1 = xi * cx.real + yi * cx.imag
                S[i, 1] += aki * aki * (f1 * rk * xi - aki * (f2 * rk1 * xi * rcx1 + f3 * rk * cx.real))
                S[i, 3] += aki * aki * (f1 * rk * yi - aki * (f2 * rk1 * yi * rcx1 + f3 * rk * cx.imag))
                zi += aki * aki * rk * (f4 * r - f5 * aki * rcx1)
            px = S[i, 1]; py = S[i, 3]
            dpz = pxy2dpz(px, py)
            al1 = alb / (1.0 + dpz)
            S[i, 0] = xi + px * al1
            S[i, 2] = yi + py * al1
            S[i, 4] = zi + dpz * al1 - S[i, 6] * alb
            # kick 2
            aki = akf / (1.0 + S[i, 5])
            if kord == 2:
                x = S[i, 0]; y = S[i, 2]
                S[i, 1] -= aki * (x - y) * (x + y)
                S[i, 3] += 2.0 * aki * x * y
            else:
                cx = complex(S[i, 0], -S[i, 2])**kord
                S[i, 1] -= aki * cx.real
                S[i, 3] -= aki * cx.imag
            drift1(i, S, ala)
        if fringe:
            ttfrin_n_1(i, S, 2 * kord + 2, -ak, L)
        if krad:
            tradk_core1(i, S, L, 0.0, fluc, rs, p0, anrad, cuc)
            setr0_1(i, S)
        texit1(i, S, dx, dy, theta, cost, sint)


# ============================================================================ BEND
@njit(**NJ)
def tbendcore1(i, S, al, phi0, cosp1, sinp1, cosp2, sinp2, mfring, fringe,
               cosw, sinw, sqwh, sinwp1, rhob, rho0, drhob, fb1, fb2):
    if (mfring > 0 or mfring == -1) and fb1 != 0.0:
        dxfr1 = fb1 * fb1 / rhob / 24.0
        dyfr1 = fb1 / (rhob * rhob) / 6.0
        dzfr1 = dxfr1 * sinp1
        dyfra1 = 4.0 * dyfr1 / (fb1 * fb1) if fringe else 0.0
    else:
        dxfr1 = 0.0; dyfr1 = 0.0; dyfra1 = 0.0; dzfr1 = 0.0
    if (mfring > 0 or mfring == -2) and fb2 != 0.0:
        dxfr2 = fb2 * fb2 / rhob / 24.0
        dyfr2 = fb2 / (rhob * rhob) / 6.0
        dzfr2 = dxfr2 * sinp2
        dyfra2 = 4.0 * dyfr2 / (fb2 * fb2) if fringe else 0.0
    else:
        dxfr2 = 0.0; dyfr2 = 0.0; dyfra2 = 0.0; dzfr2 = 0.0
    if cosp1 * cosp2 > 0.0:
        dcosp = (sinp2 - sinp1) * (sinp2 + sinp1) / (cosp1 + cosp2)
    else:
        dcosp = cosp1 - cosp2
    sinp2wp1 = sinp2 + sinwp1
    sinp2p1 = sinp2 + sinp1
    cosp1p2 = cosp1 * cosp2
    smin = 1e-4
    dp = S[i, 5]
    pp = 1.0 + dp
    pinv = 1.0 / pp
    rhoe = rhob * pp
    pxi = S[i, 1]
    pyi = S[i, 3]
    dpzi = pxy2dpz(pxi, pyi)
    pzi = 1.0 + dpzi
    sp1 = sinp1 / pzi
    x1 = S[i, 0] / (cosp1 - pxi * sp1)
    dz1 = x1 * sp1
    y1 = S[i, 2] + pyi * dz1
    px1 = pxi * cosp1 + pzi * sinp1
    x1 = x1 + dxfr1 * dp * pinv
    py1 = pyi + (dyfr1 - dyfra1 * y1 * y1) * y1 * pinv * pinv
    z1 = S[i, 4] - dz1 + (dxfr1 * px1 + (0.5 * dyfr1 - 0.25 * dyfra1 * y1 * y1) * y1 * y1 * pinv) * pinv - dzfr1
    pv1sqi = 1.0 / fmax(smin, 1.0 - px1 * px1)
    fa = y1 / rhoe * math.sqrt(pv1sqi)
    f = (1.0 - (y1 / rhob)**2 / 6.0) * fa
    ff = 0.25 * (f + fa) * y1 * pv1sqi
    x2 = x1 + ff
    py2 = py1 - px1 * f
    z2 = z1 - px1 * ff
    dph2 = sqrt1(-py2 * py2)
    ph2 = 1.0 + dph2
    dpx2 = pxi * cosp1 + (dpzi - dph2) * sinp1
    pz2 = 1.0 + pxy2dpz(px1, py2)
    drho1 = drhob + rhob * dp
    drho = drho1 + rhoe * dph2
    t2 = (px1 + ph2 * sinp1) / (pz2 + ph2 * cosp1)
    dpx3a = (x2 * sinw - drho * sinp2wp1) / rhoe
    dpx3 = dpx3a - dpx2 * (cosw - sinw * t2)
    px3 = ph2 * sinp2 + dpx3
    sqwht2 = sqwh + sinw * t2
    px1px3 = ph2 * sinp2p1 + dpx3a + dpx2 * sqwht2
    dpz3 = pxy2dpz(px3, py2)
    pz3 = 1.0 + dpz3
    dpz32 = px1px3 * (px1 - px3) / (pz2 + pz3)
    pz3ph2 = pz3 + ph2 * cosp2
    pz2px1 = pz2 * cosp1 + px1 * sinp1
    t3 = (px3 + ph2 * sinp2) / pz3ph2
    t2t3 = ((ph2 * sinp2 + px1px3) / pz3ph2
            + (ph2 * sinp1 + px1 * (dpz32 - ph2 * dcosp) / pz3ph2) / (pz2 + ph2 * cosp1))
    t4 = (cosp2 + t3 * sinp2) * pz2px1
    x3 = (x2 * (cosw - rho0 / rhoe * t3 * sinw / ph2)
          + (rho0 * (cosw * t2t3 + sinw * (1.0 - t2 * t3)) * dpx2
             - drho * (-sinp2wp1 * rho0 / rhoe * t3 - dpz32 - sqwh * pz2 - sinw * px1)) / ph2)
    da = asinz_e(
        (dpx2 * (sinp1 * (t2 * (pz3 * cosp2 + px3 * sinp2) - px1 * cosp2)
                 - sinp2 * (t3 * pz2px1 - px3 * cosp1)
                 + cosp1p2 * dpz32
                 + t4 * sqwht2)
         + dpx3a * t4) / (ph2 * ph2))
    phi0a = phi0 + da
    y3 = y1 + py2 * rhoe * phi0a
    z3 = z2 - phi0 * drho1 - da * rhoe - S[i, 6] * al
    pv2sqi = 1.0 / fmax(smin, 1.0 - px3 * px3)
    fa = y3 / rhoe * math.sqrt(pv2sqi)
    f = (1.0 - (y3 / rhob)**2 / 6.0) * fa
    ff = 0.25 * (f + fa) * y3 * pv2sqi
    py4 = py2 - px3 * f
    z4 = z3 - px3 * ff
    x4 = x3 - ff - dxfr2 * dp * pinv
    py4 = py4 + (dyfr2 - dyfra2 * y3 * y3) * y3 * pinv * pinv
    z4 = z4 + (dxfr2 * px3 + (0.5 * dyfr2 - 0.25 * dyfra2 * y3 * y3) * y3 * y3 * pinv) * pinv - dzfr2
    dpz4 = pxy2dpz(px3, py4)
    pxo = -cosp2 * dpx3 + sinp2 * (dpz4 - dpz3 - dpx3 * t3)
    dz4 = x4 * sinp2 / (cosp2 * (1.0 + dpz4) + sinp2 * px3)
    S[i, 0] = x4 * cosp2 + pxo * dz4
    S[i, 1] = pxo
    S[i, 2] = y3 + py4 * dz4
    S[i, 3] = py4
    S[i, 4] = z4 - dz4


@njit(**NJ)
def tbshift1(i, S, dx, dy, phi0, cost, sint, ent):
    if dx != 0.0 or dy != 0.0:
        phih = phi0 * 0.5
        sph = math.sin(phih)
        dcph = -2.0 * math.sin(phih * 0.5)**2
        st1 = sint if ent else -sint
        dxa = dx * cost - dy * st1
        dx1 = dx + dcph * cost * dxa
        dy1 = dy - dcph * st1 * dxa
        ds = dxa * sph
        if ent:
            al = ds / (1.0 + pxy2dpz(S[i, 1], S[i, 3]))
            x1 = S[i, 0] + S[i, 1] * al - dx1
            S[i, 2] = S[i, 2] + S[i, 3] * al - dy1
            S[i, 4] = S[i, 4] - al
            S[i, 0] = x1 * cost - S[i, 2] * sint
            S[i, 2] = x1 * sint + S[i, 2] * cost
            px1 = S[i, 1]
            S[i, 1] = px1 * cost - S[i, 3] * sint
            S[i, 3] = px1 * sint + S[i, 3] * cost
        else:
            x1 = S[i, 0] * cost - S[i, 2] * sint
            S[i, 2] = S[i, 0] * sint + S[i, 2] * cost
            px1 = S[i, 1]
            S[i, 1] = px1 * cost - S[i, 3] * sint
            S[i, 3] = px1 * sint + S[i, 3] * cost
            al = ds / (1.0 + pxy2dpz(S[i, 1], S[i, 3]))
            S[i, 0] = x1 + S[i, 1] * al - dx1
            S[i, 2] = S[i, 2] + S[i, 3] * al - dy1
            S[i, 4] = S[i, 4] - al
    elif sint != 0.0 or cost != 1.0:
        y1 = S[i, 2]
        S[i, 2] = S[i, 0] * sint + y1 * cost
        S[i, 0] = S[i, 0] * cost - y1 * sint
        px1 = S[i, 1]
        S[i, 1] = px1 * cost - S[i, 3] * sint
        S[i, 3] = px1 * sint + S[i, 3] * cost


@njit(**NJ)
def tbrot1(i, S, alg, phig, dtheta, dchi2, ent):
    c0 = math.cos(phig)
    s0 = math.sin(phig)
    s2, xs2, c2, dc2 = xsincos(dchi2)
    st, xst, ct, dct = xsincos(dtheta)
    dr11 = dct * c0 * c0 + dc2 * s0 * s0 + c0 * s0 * s2 * st
    dr22 = c2 * dct + dc2
    dr33 = dc2 * c0 * c0 + dct * s0 * s0 - c0 * s0 * s2 * st
    r11 = dr11 + 1.0; r22 = dr22 + 1.0; r33 = dr33 + 1.0
    if ent:
        r12 = s0 * s2 - c0 * c2 * st
        r13 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r21 = -ct * s0 * s2 + c0 * st
        r23 = c0 * ct * s2 + s0 * st
        r31 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r32 = -c0 * s2 - c2 * s0 * st
    else:
        r21 = s0 * s2 - c0 * c2 * st
        r31 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r12 = -ct * s0 * s2 + c0 * st
        r32 = c0 * ct * s2 + s0 * st
        r13 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r23 = -c0 * s2 - c2 * s0 * st
    dl = alg if phig == 0.0 else alg * s0 / phig
    xl0 = dl * s0
    zi = -dl * c0
    dxz = r13 * zi
    dyz = r23 * zi
    dzz = dr33 * zi
    xi = S[i, 0] + xl0
    yi = S[i, 2]
    pxi = S[i, 1]
    pyi = S[i, 3]
    pzi = 1.0 + pxy2dpz(pxi, pyi)
    xf = dr11 * xi + r12 * yi + dxz + S[i, 0]
    yf = r21 * xi + dr22 * yi + dyz + yi
    zf = r31 * xi + r32 * yi + dzz
    pxf = r11 * pxi + r12 * pyi + r13 * pzi
    pyf = r21 * pxi + r22 * pyi + r23 * pzi
    pzf = r31 * pxi + r32 * pyi + r33 * pzi
    S[i, 1] = pxf
    S[i, 3] = pyf
    S[i, 0] = xf - pxf / pzf * zf
    S[i, 2] = yf - pyf / pzf * zf
    S[i, 4] = S[i, 4] + zf / pzf


@njit(**NJ)
def tbendal(n, ndiv, f1r, f2r, aln):
    if 0 < n <= ndiv:
        return aln, aln
    if n == -1:
        return RBL * f1r, f1r
    if n == 0:
        return RBH * f1r, f1r
    if n == ndiv + 1:
        return RBH * f2r, f2r
    return RBL * f2r, f2r


@njit(**NJ)
def k_bend0(S, n, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2,
            cosw, sinw, sqwh, sinwp1, fb1, fb2, mfring, fringe, krad, first, last, eps0,
            dx, dy, cost, sint, dtheta, dchi2, alg, phig, fluc, rs, p0, h0, anrad, cuc):
    """tbend (K1=0): tbshift/tbrot → tbendcore (放射ありは tbendr のスライス) → 戻し"""
    for i in range(n):
        tbshift1(i, S, dx, dy, phi0, cost, sint, True)
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg, phig, dtheta, dchi2, True)
    rhob = al / phib
    rho0 = al / phi0
    drhob = -rho0 * ak0 / phib
    n1 = 1; n2 = 0; ndiv = 1
    f1r = 0.0; f2r = 0.0
    if krad:
        for i in range(n):
            setr0_1(i, S)
        if first and fb1 != 0.0:
            n1 = -1
            f1r = fb1 * 0.5
        if last and fb2 != 0.0:
            n2 = 2
            f2r = fb2 * 0.5
        phir = phib * (al - f1r - f2r) / al
        ndiv = fmin(1024, ndivrad(phir, 0.0, 0.0, eps0, h0))
    n2 = ndiv + n2
    if n1 == n2:
        for i in range(n):
            tbendcore1(i, S, al, phi0, cosp1, sinp1, cosp2, sinp2, mfring, fringe,
                       cosw, sinw, sqwh, sinwp1, rhob, rho0, drhob, fb1, fb2)
        if krad:
            tradk_all(S, n, al, phi0, fluc, rs, p0, anrad, cuc)
    else:
        aln = (al - f1r - f2r) / ndiv
        coswn = 0.0; sinwn = 0.0; sqwhn = 0.0; sinwpn = 0.0
        for m in range(n1, n2 + 1):
            mfr1 = 0
            psi1n = 0.0; psi2n = 0.0
            cosp1n = 1.0; sinp1n = 0.0; cosp2n = 1.0; sinp2n = 0.0
            alx, alr = tbendal(m, ndiv, f1r, f2r, aln)
            if m == n1:
                if mfring > 0 or mfring == -1:
                    mfr1 = -1
                psi1n = psi1; cosp1n = cosp1; sinp1n = sinp1
            elif m == n2:
                if mfring > 0 or mfring == -2:
                    mfr1 = -2
                psi2n = psi2; cosp2n = cosp2; sinp2n = sinp2
            phi0n = alx / al * phi0
            if m <= 2 or m >= ndiv:
                wn = phi0n - psi1n - psi2n
                sinwn, xsinwn, coswn, sqwhn = xsincos(wn)
                sinwpn = math.sin(phi0n - psi2n)
            for i in range(n):
                tbendcore1(i, S, alx, phi0n, cosp1n, sinp1n, cosp2n, sinp2n, mfr1, fringe,
                           coswn, sinwn, -sqwhn, sinwpn, rhob, rho0, drhob, fb1, fb2)
            tradk_all(S, n, alr, phi0n, fluc, rs, p0, anrad, cuc)
    for i in range(n):
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg - al, phig - phi0, dtheta, dchi2, False)
        tbshift1(i, S, -dx, -dy, -phi0, cost, -sint, False)


@njit(**NJ)
def k_bendi(S, n, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak, dx, dy, theta, dtheta, cost, sint,
            dchi2, alg, phig, fb1, fb2, mfring, enarad, fringe, eps0, first, last,
            fluc, rs, p0, h0, anrad, cuc):
    """tbendi (K1≠0 の BEND)。粒子外側・スライス内側、最後の tradk は全粒子分まとめて"""
    tanp1 = sinp1 / cosp1
    tanp2 = sinp2 / cosp2
    rhob = al / phib
    rho0 = al / phi0
    eps = 1e-3 if eps0 == 0.0 else 1e-3 * eps0
    ndiv0 = 1 + int(abs(phi0 / eps))
    krad = enarad and al != 0.0
    f1r = 0.0; f2r = 0.0
    n1 = 1; n2 = 0
    alc = al
    ndiv = ndiv0
    if krad:
        if first and fb1 != 0.0:
            f1r = 0.5 * fb1
            n1 = -1
        if last and fb2 != 0.0:
            f2r = 0.5 * fb2
            n2 = 2
        alc = al - f1r - f2r
        ndiv = fmax(ndiv0, ndivrad(phib * alc / al, ak * alc / al, 0.0, eps0, h0))
    ndiv = fmin(1000, ndiv)
    n2 = ndiv + n2
    af = 1.0 if fringe else 0.0
    aln = alc / ndiv
    alr_e = 0.0; phi_e = 0.0
    for i in range(n):
        tbshift1(i, S, dx, dy, phi0, cost, sint, True)
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg, phig, dtheta, dchi2, True)
        if krad:
            setr0_1(i, S)
        if fringe and mfring > -4 and mfring != 2:
            ttfrin4_1(i, S, ak, al)
        if fb1 != 0.0:
            dxfr1 = fb1 * fb1 / rhob / 24.0
            dyfr1 = fb1 / (rhob * rhob) / 6.0
            dyfra1 = 4.0 * dyfr1 / (fb1 * fb1)
            dp = S[i, 5]
            pr = 1.0 + dp
            y = S[i, 2]
            S[i, 0] = S[i, 0] + dxfr1 * dp / pr
            S[i, 3] = S[i, 3] + (dyfr1 - dyfra1 * y * y) * y / (pr * pr)
            S[i, 4] = S[i, 4] + (dxfr1 * S[i, 1] + (0.5 * dyfr1 - 0.25 * dyfra1 * y * y) * y * y / pr) / pr
        m_dp = S[i, 5]
        m_p = 1.0 + m_dp
        rhoe = rhob * m_p            # bendib の rhoe: 入口で 1 回だけ決める (放射で dp が変わっても更新しない)
        yi = S[i, 2]
        f = yi / rhoe
        fpx = af * S[i, 1]
        m_pyi = S[i, 3] - (tanp1 + fpx) * f
        ff = af * yi * f * 0.5
        m_xi = S[i, 0] + ff
        m_zi = S[i, 4] - ff * fpx
        m_pxi = S[i, 1] + tanp1 * m_xi / rhoe
        m_yi = yi
        akx0 = 0.0; alx0 = 0.0; phix0 = 0.0
        # tbendiinit の値
        rhosq = 0.0; drhop = 0.0; akk = 0.0; akxsq = 0.0; akysq = 0.0
        dcx = 0.0; dcxkx = 0.0; sxkx = 0.0; xsxkx = 0.0; dcy = 0.0; aksy = 0.0; syky = 0.0
        for m in range(n1, n2 + 1):
            alx, alr = tbendal(m, ndiv, f1r, f2r, aln)
            akxn = ak * alx / al
            phixn = phi0 * alx / al
            # --- tbendiinit(m, akxn, alx)  (m.p は更新、rhoe はそのまま)
            m_p = 1.0 + m_dp
            rhosq = rho0 * rhoe
            drhop = (rhoe - rho0) / rhosq
            akk = akxn / alx
            akxsq = -1.0 / rhosq - akk / m_p
            akysq = akk / m_p
            if akxsq > 0.0:
                akx = math.sqrt(akxsq)
                phix = akx * alx
                dcx = 2.0 * math.sinh(0.5 * phix)**2
                spx = math.sinh(phix)
                xspx = phix - spx
                dcxkx = dcx / akxsq
                sxkx = spx / akx
                xsxkx = xspx / akx / akxsq
            elif akxsq < 0.0:
                akx = math.sqrt(-akxsq)
                phix = akx * alx
                dcx = -2.0 * math.sin(0.5 * phix)**2
                spx = math.sin(phix)
                xspx = phix - spx
                dcxkx = dcx / akxsq
                sxkx = spx / akx
                xsxkx = -xspx / akx / akxsq
            else:
                dcx = 0.0
                dcxkx = 0.5 * alx * alx
                sxkx = alx
                xsxkx = -1.0 / 6.0 * alx**3
            if akysq > 0.0:
                aky = math.sqrt(akysq)
                phiy = aky * alx
                dcy = 2.0 * math.sinh(0.5 * phiy)**2
                spy = math.sinh(phiy)
                aksy = aky * spy
                syky = spy / aky
            elif akysq < 0.0:
                aky = math.sqrt(-akysq)
                phiy = aky * alx
                dcy = -2.0 * math.sin(0.5 * phiy)**2
                spy = math.sin(phiy)
                aksy = -aky * spy
                syky = spy / aky
            else:
                dcy = 0.0
                aksy = 0.0
                syky = alx
            # --- tbendicorr(m, (akxn+akx0)/2, (alx+alx0)/2, (phixn+phix0)/2)
            ak1c = (akxn + akx0) * 0.5
            alc1 = (alx + alx0) * 0.5
            phic = (phixn + phix0) * 0.5
            s = fmin(0.9999, m_pxi * m_pxi + m_pyi * m_pyi)
            dpzi = sqrt1(-s)
            pzi = 1.0 + dpzi
            m_zi = m_zi + (0.5 * s + (1.0 + m_xi / rho0) * dpzi / pzi) * alc1
            m_xi = (-alc1 * dpzi * m_pxi + pzi * m_xi) / (pzi - phic * m_pxi)
            m_pxi = m_pxi + phic * dpzi
            m_yi = m_yi + (-alc1 * dpzi + m_xi * phic) * m_pyi / pzi
            xr = m_xi / rho0
            m_pxi = m_pxi - ak1c * m_xi * xr * (0.5 - xr * (2.0 - xr) / 12.0) / m_p
            akx0 = akxn; alx0 = alx; phix0 = phixn
            # --- tbendibody(m, alx0)
            xiksq = drhop + akxsq * m_xi
            dxf = dcxkx * xiksq + sxkx * m_pxi
            dpxf = sxkx * xiksq + dcx * m_pxi
            dyf = dcy * m_yi + syky * m_pyi
            dpyf = aksy * m_yi + dcy * m_pyi
            hi = (0.5 * (m_pxi * m_pxi + m_pyi * m_pyi - akxsq * m_xi * m_xi - akysq * m_yi * m_yi)
                  - drhop * m_xi)
            m_zi = (m_zi - 0.25 * (dxf * m_pxi + m_xi * dpxf + dxf * dpxf
                                   + dyf * m_pyi + m_yi * dpyf + dyf * dpyf
                                   - (5.0 / rho0 - 1.0 / rhoe)
                                   * (drhop * xsxkx - sxkx * m_xi - dcxkx * m_pxi))
                    - hi * alx0 * 0.5)
            m_xi += dxf
            m_pxi += dpxf
            m_yi += dyf
            m_pyi += dpyf
            if krad and m != n2:
                S[i, 0] = m_xi; S[i, 1] = m_pxi; S[i, 2] = m_yi; S[i, 3] = m_pyi; S[i, 4] = m_zi; S[i, 5] = m_dp
                tradk_core1(i, S, alr, phixn, fluc, rs, p0, anrad, cuc)
                m_xi = S[i, 0]; m_pxi = S[i, 1]; m_yi = S[i, 2]; m_pyi = S[i, 3]; m_zi = S[i, 4]; m_dp = S[i, 5]
                setr0_1(i, S)
            alr_e = alr
            phi_e = phixn
        # tbendicorr(m, akx0/2, alx0/2, phix0/2)
        ak1c = akx0 * 0.5
        alc1 = alx0 * 0.5
        phic = phix0 * 0.5
        s = fmin(0.9999, m_pxi * m_pxi + m_pyi * m_pyi)
        dpzi = sqrt1(-s)
        pzi = 1.0 + dpzi
        m_zi = m_zi + (0.5 * s + (1.0 + m_xi / rho0) * dpzi / pzi) * alc1
        m_xi = (-alc1 * dpzi * m_pxi + pzi * m_xi) / (pzi - phic * m_pxi)
        m_pxi = m_pxi + phic * dpzi
        m_yi = m_yi + (-alc1 * dpzi + m_xi * phic) * m_pyi / pzi
        xr = m_xi / rho0
        m_pxi = m_pxi - ak1c * m_xi * xr * (0.5 - xr * (2.0 - xr) / 12.0) / m_p
        m_zi = m_zi - S[i, 6] * al
        px = m_pxi + tanp2 * m_xi / rhoe
        f = m_yi / rhoe
        fpx = af * px
        S[i, 1] = px
        S[i, 2] = m_yi
        S[i, 3] = m_pyi - (tanp2 - fpx) * f
        ff = af * m_yi * f * 0.5
        S[i, 0] = m_xi - ff
        S[i, 4] = m_zi + ff * fpx
        S[i, 5] = m_dp
        if fb2 != 0.0:
            dxfr2 = fb2 * fb2 / rhob / 24.0
            dyfr2 = fb2 / (rhob * rhob) / 6.0
            dyfra2 = 4.0 * dyfr2 / (fb2 * fb2)
            dp = S[i, 5]
            pr = 1.0 + dp
            y = S[i, 2]
            S[i, 0] = S[i, 0] - dxfr2 * dp / pr
            S[i, 3] = S[i, 3] + (dyfr2 - dyfra2 * y * y) * y / (pr * pr)
            S[i, 4] = S[i, 4] - (dxfr2 * S[i, 1] - (0.5 * dyfr2 - 0.25 * dyfra2 * y * y) * y * y / pr) / pr
        if fringe and mfring > -4 and mfring != 1:
            ttfrin4_1(i, S, -ak, al)
    if krad:
        tradk_all(S, n, alr_e, phi_e, fluc, rs, p0, anrad, cuc)
    for i in range(n):
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg - al, phig - phi0, dtheta, dchi2, False)
        tbshift1(i, S, -dx, -dy, -phi0, cost, -sint, False)


# ============================================================================ 空洞 (tcav)
@njit(**NJ)
def cav_drift1(i, S, alx):
    px = S[i, 1]; py = S[i, 3]
    a = px * px + py * py
    dpz = sqrt1(-a)
    al1 = alx / (1.0 + dpz)
    S[i, 0] += px * al1
    S[i, 2] += py * al1
    S[i, 4] += dpz * al1 - S[i, 6] * alx


@njit(**NJ)
def tcavfrin1(i, S, al, v, w, dphis, offset, p0, h0):
    vf = v / al / p0 * 0.5
    s0 = math.sin(dphis) - offset
    dp1r = S[i, 5]
    p1r = 1.0 + dp1r
    p1 = p0 * p1r
    if S[i, 6] > 0.1:
        h1 = p2h(p1)
    else:
        h1 = p1r * h0 / (1.0 - S[i, 6] + S[i, 10])
    t = -S[i, 4] * h1 / p1
    phw = w * t
    dpt = vf * (math.sin(phw - dphis) + s0) / p1r
    S[i, 1] += S[i, 0] * dpt
    S[i, 3] += S[i, 2] * dpt
    dh = -w * p0 * vf * (S[i, 0]**2 + S[i, 2]**2) * 0.5 * math.cos(phw - dphis)
    h2 = h1 + dh
    a = fmax(dh * (h1 + h2) / (p1 * p1), -1.0)
    dpr = sqrt1(a)
    dp2r = fmax(dp1r + p1r * dpr, 1e-10 - 1.0)
    p2r = 1.0 + dp2r
    S[i, 5] = dp2r
    S[i, 4] = S[i, 4] * p2r / p1r * h1 / h2
    S[i, 6] = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + S[i, 10]


@njit(**NJ)
def k_cav(S, n, al, v, w, wi, dphis, offset, offset1, tlim, fringe, mfring, p0, h0):
    """tcav (vnominal=0, ndiv=1)。v = vc/amass·|charge|、w = 2π f/c、dphis = phis − phic、offset は rad で切替済み"""
    wsn = 1.0
    vcorr = v * (w * (0.5 / p0 + 0.5 / p0))**2 / 4.0
    v20a = vcorr
    v02a = vcorr
    for i in range(n):
        if al != 0.0 and fringe and mfring >= 0 and mfring != 2:
            tcavfrin1(i, S, al, v, w, dphis, offset, p0, h0)
        if al != 0.0:
            cav_drift1(i, S, al * wsn * 0.5)
        dp1r = S[i, 5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        if S[i, 6] > 0.1:
            h1 = p2h(p1)
        else:
            h1 = p1r * h0 / (1.0 - S[i, 6] + S[i, 10])
        t = fmin(tlim, fmax(-tlim, -S[i, 4] * h1 / p1))
        phh = 0.5 * w * t
        dh = fmax(1.0 + 3.83e-12 - h1,
                 (v + v20a * S[i, 0]**2 + v02a * S[i, 2]**2)
                 * (-2.0 * math.sin(phh) * math.cos(phh - dphis) + offset) * wsn)
        h2 = h1 + dh
        a = fmax(dh * (h1 + h2) / (p1 * p1), -1.0)
        dpr = sqrt1(a)
        dp2r = dp1r + p1r * dpr
        p2r = 1.0 + dp2r
        S[i, 5] = dp2r
        S[i, 6] = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + S[i, 10]
        av = -(math.cos(2.0 * phh - dphis) * wi - offset1 * t) * wsn / p0
        dpx = 2.0 * v20a * S[i, 0] * av
        dpy = 2.0 * v02a * S[i, 2] * av
        S[i, 1] = (S[i, 1] * p1r + dpx) / p2r
        S[i, 3] = (S[i, 3] * p1r + dpy) / p2r
        S[i, 4] = -t * p2r * p0 / h2
        if al != 0.0:
            cav_drift1(i, S, al * wsn * 0.5)
            if fringe and mfring >= 0 and mfring != 1:
                tcavfrin1(i, S, al, -v, w, dphis, offset, p0, h0)


@njit(**NJ)
def k_steer0(S, n, phib, dx, dy, theta, cost, sint):
    """tthin(nord=2, ak=K0): px −= K0/(1+δ) (K0>0 は px を減らす)。ずれ・ロールは TENT/TEXIT (tthin の tent/texit)"""
    for i in range(n):
        tent1(i, S, dx, dy, theta, cost, sint)
        S[i, 1] -= phib / (1.0 + S[i, 5])
        texit1(i, S, dx, dy, theta, cost, sint)


XTHASX = (2.0**-52 * 13312.0 / 231.0 / 6.0)**0.1       # tfloor.f xthasx


@njit(**NJ)
def asinx(x):
    """mathfun asinx = asin(x) − x (小さい x は級数)"""
    if abs(x) < XTHASX:
        x2 = x * x
        return x * x2 * (1.0 / 6.0 + x2 * (3.0 / 40.0 + x2 * (5.0 / 112.0 + x2 * (35.0 / 1152.0 + x2 * 63.0 / 2816.0))))
    return math.asin(fmin(1.0, fmax(-1.0, x))) - x


@njit(**NJ)
def k_steer(S, n, al, phib, dx, dy, theta, cost, sint, tanp1, krad, eps0, fluc, rs, p0, h0, anrad, cuc):
    """tsteer.f (track6._tsteer_ens と同じ。F1/FB1/FB2 なし)。スライス外側・粒子内側、各スライスの後に全粒子 tradk。
    入口エッジは tanp1 (最初のスライスだけ)。出口の tanp2 は fb2≠0 のときだけなので 0 のまま"""
    rhob = al / phib
    for i in range(n):
        tent1(i, S, dx, dy, theta, cost, sint)
    ndiv = 1
    if krad:
        for i in range(n):
            setr0_1(i, S)
        ndiv = min(1000, max(1, ndivrad(phib, 0.0, 0.0, eps0, h0)))
    aln = al / ndiv
    for m in range(1, ndiv + 1):
        for i in range(n):
            dp = S[i, 5]
            pr = 1.0 + dp
            rhoe = rhob * pr
            pxi = S[i, 1] - tanp1 * S[i, 0] / rhoe
            S[i, 3] += (tanp1 / rhoe) * S[i, 2]
            dpv1 = sqrt1(-pxi**2)
            pv1 = 1.0 + dpv1
            fa = S[i, 2] / rhoe / pv1
            f = (1.0 - (S[i, 2] / rhob)**2 / 6.0) * fa
            ff = (0.5 - (S[i, 2] / rhob)**2 / 24.0) * S[i, 2] * fa / pv1**2
            pyi = S[i, 3] + pxi * f
            xi = S[i, 0] - ff
            S[i, 4] += ff * pxi
            dpz1 = pxy2dpz(pxi, pyi)
            pz1 = 1.0 + dpz1
            dpx = aln / rhoe
            pxf = pxi + dpx
            dpz2 = pxy2dpz(pxf, pyi)
            pz2 = 1.0 + dpz2
            d = pxf * pz1 + pxi * pz2
            if d == 0.0 or pxi * pxf < 0.0:
                phsq = pxi**2 + pz1**2
                u = (dpx + pxf * dpz1 - pxi * dpz2) / phsq
                w = (pxf * dpz1 - pxi * dpz2 + dpx * pyi**2) / phsq
            else:
                u = dpx * (pxf + pxi) / d
                w = -dpx * (pxf * dpz1 + pxi * dpz2) / d
            dl = rhoe * (asinx(u) + w)
            S[i, 0] = xi + aln * (pxf + pxi) / (pz2 + pz1)
            S[i, 2] += pyi * (aln + dl)
            dpv2 = sqrt1(-pxf**2)
            pv2 = 1.0 + dpv2
            fa = S[i, 2] / rhoe / pv2
            f = (1.0 - (S[i, 2] / rhob)**2 / 6.0) * fa
            ff = (0.5 - (S[i, 2] / rhob)**2 / 24.0) * S[i, 2] * fa / pv2**2
            S[i, 3] = pyi - pxf * f
            S[i, 0] += ff
            S[i, 1] = pxf
            S[i, 4] -= (dl + S[i, 6] * aln + ff * pxf)
        tanp1 = 0.0
        if krad:
            tradk_all(S, n, aln, 0.0, fluc, rs, p0, anrad, cuc)
    for i in range(n):
        texit1(i, S, dx, dy, theta, cost, sint)


# ============================================================================ 損失判定 (tapert) と詰め替え
@njit(**NJ)
def k_apert(S, n, kid, lost_turn, kturn, alost, zlost, log, nlog, lcur):
    """limittrack + tapert (ax=ay=0)。落ちた粒子は末尾の生きている粒子と入れ替え (SAD と同じ)。
    kid[i] は粒子の ID (スロット→ID)、lost_turn[ID] に kturn を記録。戻り: 新しい n
    log (m,6) / nlog[0]: MAIN の TRACK 用の損失報告 (tapert_report_dropped: スロット順に ID, l_track=lcur, x, y, z, kturn)。
    log の行数が 0 なら記録しない (DA / TrackParticles は報告を出さない)"""
    xlimit = alost * 3.0
    drop = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        for k in (0, 2):
            v = S[i, k]
            if v != v:
                S[i, k] = xlimit
            elif abs(v) > xlimit:
                S[i, k] = math.copysign(xlimit, v)
        for k in (1, 3):
            v = S[i, k]
            if v != v:
                S[i, k] = 0.9999
            elif abs(v) > 0.9999:
                S[i, k] = math.copysign(0.9999, v)
        v = S[i, 4]
        if v != v:
            S[i, 4] = 1e10
        elif abs(v) > 1e10:
            S[i, 4] = math.copysign(1e10, v)
        x = S[i, 0]; y = S[i, 2]; z = S[i, 4]
        ok = (-alost < x) and (x < alost) and (-alost < y) and (y < alost) and (abs(z) <= zlost)
        if not ok:
            drop[i] = True
            lost_turn[kid[i]] = kturn
            if nlog[0] < log.shape[0]:
                r = nlog[0]
                log[r, 0] = kid[i]
                log[r, 1] = lcur
                log[r, 2] = x
                log[r, 3] = y
                log[r, 4] = z
                log[r, 5] = kturn
                nlog[0] = r + 1
    i = 0
    np_ = n
    while i < np_:
        if drop[i]:
            while i < np_ - 1 and drop[np_ - 1]:
                np_ -= 1
            if not drop[np_ - 1] and np_ - 1 != i:
                for k in range(S.shape[1]):
                    tmp = S[i, k]; S[i, k] = S[np_ - 1, k]; S[np_ - 1, k] = tmp
                tk = kid[i]; kid[i] = kid[np_ - 1]; kid[np_ - 1] = tk
                drop[i] = False
                drop[np_ - 1] = True
            np_ -= 1
        i += 1
    return np_


# ============================================================================ ラティスの表
# ソレノイド領域の la (tapert) の扱い: 0 = 通常、1 = 領域の始点 SOL、2 = 領域内、3 = 終点 SOL (詳細は track6n_ext.py)
SM_NONE, SM_START, SM_IN, SM_END = 0, 1, 2, 3


def _ext():
    """MULT / ソレノイド領域のカーネル (track6n_ext.py)。必要になったときに import する
    (track6n_ext が track6n のヘルパーを import するので、循環 import を避けるため遅延させる)"""
    from . import track6n_ext
    return track6n_ext


class Prog:
    """要素ごとの (種類, パラメータ) の列。MULT (K_MULT) とソレノイド領域 (K_SOL_*) は track6n_ext で compile する"""
    def __init__(self, elements, ph):
        self.kinds = []
        self.params = []
        self.names = []
        # tturn の la (tapert までの要素数) を減らす要素か: SAD では DRIFT (go to 1011) と MARK/MONI (exit) は
        # la を減らさない。parser が L=0 の BEND (VD/HD 等) を Drift にしたものは _sad_type='BEND' で見分ける
        self.la_count = []
        self.sol_mode = []              # ソレノイド領域内の要素の la の扱い (SM_*)
        i = 0
        nel = len(elements)
        while i < nel:
            el = elements[i]
            tn = type(el).__name__
            if tn == 'Solenoid':
                if _t6._g(el, 'BOUND') == 0.0:
                    raise NotImplementedError('track6n: BOUND=0 の SOL (領域外)')
                i = _ext().compile_sol_region(self, elements, i, ph) + 1
                continue
            if tn == 'Multipole':
                k, p = _ext().compile_mult(el)
                cnt = True
            else:
                k, p = self._compile(el, ph)
                if tn in ('Marker', 'Monitor'):
                    cnt = False
                elif isinstance(el, Drift):
                    cnt = getattr(el, '_sad_type', 'DRIFT') != 'DRIFT'
                else:
                    cnt = True
            self._add(k, p, el.name, cnt, SM_NONE)
            i += 1

    def _add(self, k, p, name, cnt, sm):
        self.kinds.append(k)
        self.params.append(p)
        self.names.append(name)
        self.la_count.append(cnt)
        self.sol_mode.append(sm)

    @staticmethod
    def _compile(el, ph):
        g = _t6._g
        tn = type(el).__name__
        if isinstance(el, Drift):
            return K_DRIFT, (el.L,)
        if isinstance(el, Quadrupole):
            L = el.L
            K1L = g(el, 'K1') * L
            if L == 0.0:
                # tquad.f: al=0 → tthin(nord=4, ak=K1 (積分値), dx, dy, theta=ROT, krad=.false., fringe=.false.)
                ak = g(el, 'K1')
                if ak == 0.0:
                    return K_DRIFT, (0.0,)
                theta = g(el, 'ROTATION')
                ct, st = (math.cos(theta), math.sin(theta)) if theta != 0.0 else (1.0, 0.0)
                return K_THIN, (0.0, ak, 1, g(el, 'DX'), g(el, 'DY'), theta, ct, st, False, False)
            th2 = g(el, 'ROTATION') + (0.5 * math.pi if K1L < 0 else 0.0)
            f1in, f2in, f1out, f2out = _t6._tsetfringep(el, K1L / L) if L else (0.0, 0.0, 0.0, 0.0)
            return K_QUAD, (L, K1L, g(el, 'DX'), g(el, 'DY'), math.cos(th2), math.sin(th2),
                            g(el, 'DISFRIN') == 0.0, int(round(g(el, 'FRINGE'))), f1in, f2in, f1out, f2out,
                            g(el, 'EPS'), g(el, 'RAD') == 0.0 and L != 0.0)
        if tn in _t6._THIN_KEY:
            kn, kord = _t6._THIN_KEY[tn]
            L = el.L
            K = g(el, kn)
            ak = K * L if L != 0.0 else K
            theta = g(el, 'ROTATION')
            ct, st = (math.cos(theta), math.sin(theta)) if theta != 0.0 else (1.0, 0.0)
            return K_THIN, (L, ak, kord, g(el, 'DX'), g(el, 'DY'), theta, ct, st, g(el, 'DISFRIN') == 0.0,
                            g(el, 'RAD') == 0.0 and L != 0.0)
        if isinstance(el, Bend):
            al0 = el.L
            phi0 = el.ANGLE
            ak0 = g(el, 'K0')
            if phi0 == 0.0:
                if g(el, 'K1') != 0.0:
                    # tbend.f: K0=0 なら tquad (theta+dtheta、フリンジ係数 0) → K_QUAD、K0≠0 なら tmulti → track6 へ
                    kind, proxy = _t6._bend_k1_proxy(el)
                    if kind != 'quad':
                        raise NotImplementedError('track6n: ANGLE=0 で K1≠0 かつ K0≠0 の BEND (tmulti)')
                    L = proxy.L
                    K1L = proxy.K1 * L
                    if L == 0.0 or K1L == 0.0:
                        raise NotImplementedError('track6n: ANGLE=0 K1≠0 の BEND の L=0')
                    th2 = proxy.ROTATION + (0.5 * math.pi if K1L < 0 else 0.0)
                    return K_QUAD, (L, K1L, proxy.DX, proxy.DY, math.cos(th2), math.sin(th2),
                                    proxy.DISFRIN == 0.0, 0, 0.0, 0.0, 0.0, 0.0, proxy.EPS, proxy.RAD == 0.0 and L != 0.0)
                if ak0 == 0.0:
                    return K_DRIFT, (al0,)
                # tbend.f: tsteer / tthin には theta+dtheta (ROT + DROT) を渡す
                dx, dy = g(el, 'DX'), g(el, 'DY')
                theta = _lat._bend_theta(el) + g(el, 'DROTATE')
                cost, sint = (math.cos(theta), math.sin(theta)) if theta != 0.0 else (1.0, 0.0)
                if al0 == 0.0:
                    return K_STEER0, (ak0, dx, dy, theta, cost, sint)
                if g(el, 'F1') != 0.0 or g(el, 'FB1') != 0.0 or g(el, 'FB2') != 0.0:
                    raise NotImplementedError('track6n: ANGLE=0 の BEND (tsteer) の F1/FB1/FB2 フリンジは未対応')
                psi1 = g(el, 'E1')
                return K_STEER, (al0, -ak0, dx, dy, theta, cost, sint, math.sin(psi1) / math.cos(psi1),
                                 g(el, 'EPS'), g(el, 'RAD') == 0.0)
            if al0 == 0.0:
                raise NotImplementedError('track6n: L=0 で ANGLE≠0 の BEND')
            psi1 = g(el, 'E1')
            psi2 = g(el, 'E2')
            fb1, fb2, mfring, fringe = _lat._bend_fringe_params(el)
            if mfring <= 0:
                if mfring != -1:
                    fb1 = 0.0
                if mfring != -2:
                    fb2 = 0.0
            w = phi0 - psi1 - psi2
            al = al0
            if fb1 != 0.0 or fb2 != 0.0:
                al = al - ((phi0 * fb1)**2 + (phi0 * fb2)**2) / al / 48.0 * math.sin(0.5 * w) / math.sin(0.5 * phi0)
            cosp1, sinp1 = math.cos(psi1), math.sin(psi1)
            cosp2, sinp2 = math.cos(psi2), math.sin(psi2)
            cosw, sinw = math.cos(w), math.sin(w)
            sqwh = sinw**2 / (1.0 + cosw) if cosw >= 0.0 else 1.0 - cosw
            sinwp1 = math.sin(phi0 - psi2)
            dx, dy = g(el, 'DX'), g(el, 'DY')
            theta = _lat._bend_theta(el)        # ROT = 設計ロール ROTATE + ROTATION
            cost, sint = math.cos(theta), math.sin(theta)
            dtheta = g(el, 'DROTATE')
            dchi2 = g(el, 'CHI2')
            alg, phig = 0.5 * al0, 0.5 * phi0
            phib = phi0 + ak0
            ak1 = g(el, 'K1') * al0
            common = (al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
                      fb1, fb2, mfring, fringe, g(el, 'EPS'), dx, dy, cost, sint, dtheta, dchi2, alg, phig,
                      ak1, theta, g(el, 'RAD') == 0.0 and al0 != 0.0)
            return K_BEND, common
        if tn == 'Cavity':
            al = el.L
            vc = g(el, 'VOLT') + g(el, 'DVOLT')
            harm = g(el, 'HARM') or g(el, '_harm')
            freq = g(el, 'FREQ')
            w = 2.0 * math.pi * freq / CVELOC if harm == 0.0 else ph.omega0 * harm / CVELOC
            wi = 1.0 / w if w != 0.0 else 0.0
            phic = g(el, 'PHI') * ph.charge
            v = vc / ph.amass * abs(ph.charge)
            return K_CAV, (al, v, w, wi, phic, g(el, 'DISFRIN') == 0.0, int(round(g(el, 'FRINGE'))))
        if tn in ('Multipole', 'Solenoid', 'TCavity'):
            # Prog.__init__ が MULT/SOL を先に処理するので、ここに来るのは TCavity (と直接 _compile を呼んだ場合) だけ。
            # L=0 の MULT (薄い多極) や SOL を黙って無視しないよう、L==0 の判定より前に置く
            raise NotImplementedError(f'track6n: 要素 {tn} は未対応')
        if tn in ('Marker', 'Monitor'):
            return K_NONE, ()
        if el.L == 0.0 and tn in ('Drift', 'Sextupole', 'Octupole', 'Decapole', 'Dodecapole', 'Cavity'):
            return K_NONE, ()
        raise NotImplementedError(f'track6n: 要素 {tn} は未対応')


TRACE = None


_NOLOG = np.zeros((0, 6))                      # 損失報告を記録しないときの log (行数 0)
_NONLOG = np.zeros(1, dtype=np.int64)


def turn(prog: Prog, S, n, ph, trf0, rfsw, rad, fluc, rs, apert=None):
    """1 周。apert=(kid, lost_turn, kturn, alost, zlost[, log, nlog]) を渡すと tapert を最初の要素・15 要素ごと・周末に行い、
    新しい n を返す。ソレノイド領域の la の扱いは SAD のまま:
      tturn1 の主ループは始点 SOL k で tapert 判定 → tsol → (1020) la-1、その後 l=k+1..ke も主ループを通り
      (`sol=l<ke; go to 1020`) 各 l で「la<=0 なら tapert, la=15」→ la-1 (座標は変わらないので判定は領域の出口で)。
      tsol 内のループ (l=k+1..ke) は各要素の前に「la<=0 なら tapert, la=15、そうでなければ la-1」(tsol.f:604-615)"""
    p0, h0, anrad, cuc = ph.p0, ph.h0, ph.anrad, ph.cuc
    la = 0
    if apert is not None:
        kid, lost_turn, kturn, alost, zlost = apert[:5]
        log, nlog = apert[5:7] if len(apert) >= 7 else (_NOLOG, _NONLOG)
        nlat_ = len(prog.kinds) + 1      # SAD の nlat (要素数 + 1)。l_track の周回またぎの値 (前の周の最後 = nlat) に使う
    nreg = 0                             # 領域内の要素数 (k+1..ke) の数え上げ (主ループの後処理用)
    for l_el, (kind, pr, cnt, sm) in enumerate(zip(prog.kinds, prog.params, prog.la_count, prog.sol_mode)):
        if apert is not None:
            if sm == SM_NONE or sm == SM_START:
                if la <= 0:
                    # tapert の l_track は直前に処理した要素の番号 (tturn.f:243 の l_track=l は tapert の後)。周の最初は前の周の最後 (nlat)
                    n = k_apert(S, n, kid, lost_turn, kturn, alost, zlost, log, nlog, l_el if l_el >= 1 else nlat_)
                    if n <= 0:
                        return 0
                    la = 15
                if cnt and sm == SM_NONE:
                    la -= 1
                if sm == SM_START:
                    nreg = 0
            else:                        # tsol 内 (tsol.f:604-615)
                nreg += 1
                if la <= 0:
                    n = k_apert(S, n, kid, lost_turn, kturn, alost, zlost, log, nlog, l_el)
                    if n <= 0:
                        return 0
                    la = 15
                else:
                    la -= 1
        if kind == K_DRIFT:
            k_drift(S, n, pr[0])
        elif kind == K_QUAD:
            (L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0, radok) = pr
            k_quad(S, n, L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0,
                   rad and radok, fluc, rs, p0, h0, anrad, cuc)
        elif kind == K_THIN:
            (L, ak, kord, dx, dy, theta, ct, st, fringe, radok) = pr
            k_thin(S, n, L, ak, kord, dx, dy, theta, ct, st, fringe, rad and radok, fluc, rs, p0, anrad, cuc)
        elif kind == K_BEND:
            (al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
             fb1, fb2, mfring, fringe, eps0, dx, dy, cost, sint, dtheta, dchi2, alg, phig, ak1, theta, radok) = pr
            krad = rad and radok
            if ak1 != 0.0:
                k_bendi(S, n, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak1, dx, dy, theta, dtheta, cost, sint,
                        dchi2, alg, phig, fb1, fb2, mfring, krad, fringe, eps0, True, True,
                        fluc, rs, p0, h0, anrad, cuc)
            else:
                k_bend0(S, n, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
                        fb1, fb2, mfring, fringe, krad, True, True, eps0, dx, dy, cost, sint, dtheta, dchi2, alg, phig,
                        fluc, rs, p0, h0, anrad, cuc)
        elif kind == K_CAV:
            (al, v, w, wi, phic, fringe, mfring) = pr
            phis = w * trf0
            dphis = phis - phic
            if rad:
                offset = math.sin(dphis)
                offset1 = 0.0
            else:
                offset = math.sin(dphis) - math.sin(phis)
                offset1 = math.sin(phis)
            tlim = 1e4 if rfsw else 0.0
            k_cav(S, n, al, v, w, wi, dphis, offset, offset1, tlim, fringe, mfring, p0, h0)
        elif kind == K_STEER0:
            k_steer0(S, n, pr[0], pr[1], pr[2], pr[3], pr[4], pr[5])
        elif kind == K_STEER:
            (al, phib, dx, dy, theta, cost, sint, tanp1, eps0, radok) = pr
            k_steer(S, n, al, phib, dx, dy, theta, cost, sint, tanp1, rad and radok, eps0, fluc, rs, p0, h0, anrad, cuc)
        elif kind > K_NONE:
            # MULT / ソレノイド領域 (K_MULT, K_SOL_FRIN, K_SOL_DRIFT, K_SOL_QUAD)
            _ext().step_kind(kind, pr, S, n, ph, trf0, rfsw, rad, fluc, rs)
        if TRACE is not None:            # デバッグ用: 要素ごとのフック (切り分けスクリプトが設定する)
            TRACE(l_el, prog, S, n, la, rs)
        if apert is not None and sm == SM_END:
            la -= 1                      # 主ループ: 始点 SOL の 1020
            for _ in range(nreg):        # 主ループ: l=k+1..ke (sol=l<ke; go to 1020)
                if la <= 0:
                    n = k_apert(S, n, kid, lost_turn, kturn, alost, zlost, log, nlog, l_el + 1)
                    if n <= 0:
                        return 0
                    la = 15
                la -= 1
    if apert is not None:
        n = k_apert(S, n, kid, lost_turn, kturn, alost, zlost, log, nlog, nlat_)
    return n


def make_state(coords, ph):
    """物理座標 (np, 6) → S (np, 11)。最後の列は dvfs (ph.dvfs: tsetdvfs = RFSW ON のとき FSHIFT。tconvm/tradk/tcav の dv に足される)"""
    c = np.array(coords, dtype=float)
    npart = c.shape[0]
    S = np.zeros((npart, 11))
    S[:, 10] = ph.dvfs
    S[:, 0] = c[:, 0]
    S[:, 2] = c[:, 2]
    S[:, 4] = c[:, 4]
    S[:, 5] = c[:, 5]
    pr = 1.0 + c[:, 5]
    S[:, 1] = c[:, 1] / pr
    S[:, 3] = c[:, 3] / pr
    for i in range(npart):
        S[i, 6] = ph.dv(S[i, 5])
    S[:, 7] = S[:, 1]
    S[:, 8] = S[:, 3]
    S[:, 9] = S[:, 4]
    return S


def state_out(S, n):
    c = S[:n, :6].copy()
    pr = 1.0 + c[:, 5]
    c[:, 1] *= pr
    c[:, 3] *= pr
    return c


MAXTURN_LOST = 2**62            # lost_turn の「落ちていない」印


class TrackRun:
    """TrackParticles (tftrack.f) 用: tturn を区間 (ls..ld) ごとに呼び、tapert (LOSSAMPL/LOSSDZ、最初の要素・15 要素ごと・周末) で
    落ちた粒子を SAD と同じく末尾へ寄せて追跡から外す。初期フラグ 0 の粒子は追跡しない (座標はそのまま)。
    結果は粒子の番号順に戻し、alive は落ちていなければ 1、落ちた粒子の座標は落ちた時点の値 (limittrack で |x|,|y|≤3·LOSSAMPL、
    |px|,|py|≤0.9999、|z|≤1e10 に抑えたもの)。Prog が未対応の要素を含むときは NotImplementedError (呼び出し側が track6 に落とす)"""

    def __init__(self, elements, coords, flags, momentum_eV, alost=1.0, zlost=100.0, kturn0=1, trf0=0.0, rfsw=True,
                 rad=False, fluc=False, fshift=0.0, circ=None, charge=1.0, amass=None):
        self.elements = list(elements)
        if circ is None:
            circ = sum(el.L for el in self.elements)
        kw = {} if amass is None else {'amass': amass}
        ph = _e._Phys(momentum_eV, circ, charge=charge, **kw)
        ph.dvfs = fshift if rfsw else 0.0
        self.ph = ph
        self.trf0, self.rfsw, self.rad = trf0, rfsw, rad
        self.fl = bool(rad and fluc)
        self.alost, self.zlost = alost, zlost
        self.kturn = kturn0
        self.coords0 = np.array(coords, dtype=float).reshape(-1, 6)
        npart = self.coords0.shape[0]
        self.alive_in = [float(f) != 0.0 for f in (flags if flags is not None else [1.0] * npart)]
        idx = [i for i in range(npart) if self.alive_in[i]]
        self.kid = np.array(idx, dtype=np.int64)
        self.lost_turn = np.full(max(npart, 1), MAXTURN_LOST, dtype=np.int64)
        self.n = len(idx)
        self.S = make_state(self.coords0[idx], ph) if idx else None
        self.rs = np.array([_t6.RNG.seed], dtype=np.int64)
        Prog(self.elements, ph)                  # 未対応の要素なら NotImplementedError (全体で 1 回確かめる)
        self._progs = {}

    def segment(self, a, b):
        """tturn(ls=a, ld=b): 要素 a..b−1 を 1 回通す (tapert 付き)"""
        if self.n <= 0 or b <= a:
            return
        key = (a, b)
        prog = self._progs.get(key)
        if prog is None:
            prog = self._progs[key] = Prog(self.elements[a - 1:b - 1], self.ph)
        with np.errstate(all='ignore'):
            self.n = turn(prog, self.S, self.n, self.ph, self.trf0, self.rfsw, self.rad, self.fl, self.rs,
                          apert=(self.kid, self.lost_turn, self.kturn, self.alost, self.zlost))
        if self.n > 0 and not self.rfsw:
            self.S[:self.n, 4] = 0.0             # tturn.f:142-145: RFSW OFF は区間ごとに z=0

    def next_turn(self):
        self.kturn += 1

    def result(self):
        """(座標 (npart,6) 粒子の番号順, alive のリスト)"""
        _t6.RNG.seed = int(self.rs[0])
        out = self.coords0.copy()
        alive = [1.0 if a else 0.0 for a in self.alive_in]
        if self.S is not None and len(self.kid):
            c = state_out(self.S, self.S.shape[0])
            for s in range(len(self.kid)):
                i = int(self.kid[s])
                out[i] = c[s]
                alive[i] = 1.0 if self.lost_turn[i] == MAXTURN_LOST else 0.0
        return out, alive


def track_particles_fast(elements, coords, nturn, momentum_eV, trf0=0.0, rfsw=True, rad=False, fluc=False,
                         seed=None, charge=1.0, amass=None, progress=0, fshift=0.0, circ=None):
    """track6.track_particles と同じ仕様 (MULT/SOL があれば track6 にフォールバック)。
    circ: 全周長 (omega0 = 2πc/circ、空洞の周波数)。elements が区間のときに渡す。None なら elements の長さの和"""
    import time
    if circ is None:
        circ = sum(el.L for el in elements)
    kw = {} if amass is None else {'amass': amass}
    ph = _e._Phys(momentum_eV, circ, charge=charge, **kw)
    ph.dvfs = fshift if rfsw else 0.0
    try:
        prog = Prog(elements, ph)
    except NotImplementedError:
        return _t6.track_particles(elements, coords, nturn, momentum_eV, trf0, rfsw, rad, fluc, seed, charge, amass, progress, fshift, circ)
    if seed is not None:
        _t6.RNG.restore(seed)
    rs = np.array([_t6.RNG.seed], dtype=np.int64)
    fl = bool(rad and fluc)
    S = make_state(coords, ph)
    n = S.shape[0]
    t0 = time.time()
    for it in range(nturn):
        turn(prog, S, n, ph, trf0, rfsw, rad, fl, rs)
        if not np.all(np.isfinite(S[:n, :6])):
            raise FloatingPointError(f'track6n: {it + 1} 周目で発散')
        if not rfsw:
            S[:n, 4] = 0.0
        if progress and ((it + 1) % progress == 0 or it + 1 == nturn):
            print(f'    [track6n] turn {it + 1}/{nturn} ({time.time() - t0:.1f} s)', flush=True)
    _t6.RNG.seed = int(rs[0])
    return state_out(S, n)
