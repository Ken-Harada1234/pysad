"""
elkern.py — 要素ごとの 6D/4D 追跡カーネル (numba) (2026-11-05, docs/history/work_perf2.md)

emit6d.tturne (6D: 放射・dv 項・D/B 行列) と lattice._el_cod_step (4D の CALC: 放射なし・dv なし・z 行の組み直しなし) の
要素ごとの Python (piece の連続) を 1 回の numba 呼び出しにまとめたもの。式はすべて既存の Python 版
(lattice.py の _tbedge*/_tbendebody0*/_tsteee_body/_tqente_cod/_tsolque_body_cod/_tqfrie_cod/_tqlfre_cod/_tthine_cod/
 _exact_drift_cod/_tchge_cod/_tbrote_jac、emit6d.py の _jac_*/tradke/_tquade/_tbende_core/_tsteee/_tthine_gen) の直訳。

状態: cod(6), T(6,6), D(6,6), B(6,6), codr0(6), transr(6,6), S(3) = [u0, h1emit, dvemit]。
PH(11) = [p0, h0, anrad, cuc, crad, brhoz, emidiv, emidib, emidiq, dvfs, pgev] (emit6d._Phys)。
フラグ: want_jac (T を作る), want_db (D・B を作る), symp_fix (piece ごとに z 行をシンプレクティック条件で組み直す: 6D のみ),
       use_dv (z −= dv(δ)·L の項: 6D のみ), enarad (放射: tradke), radcod (放射で cod を変える)。
各カーネルは状態の複製で計算し、成功 (戻り値 0) のときだけ書き戻す。0 以外 (解析ヤコビアンが作れない等) は呼び出し側が
従来の Python 経路で計算し直す。
"""
import math

import numpy as np
from numba import njit

from .perfkern import tradke_db as _tradke_db, edge_free as _edge_free
from .radbeta import ab as _rb_ab
from . import quadbody as _qb
from .perfkern_bend import apply_tdb as _apply_tdb, bendecorr as _bendecorr, bendebody as _bendebody

NJ = dict(cache=True, error_model='numpy', nogil=True)

# emit6d の定数
CAVE = 8.0 / 15.0 / math.sqrt(3.0)
CUU = 11.0 / 27.0
GMIN = -0.999999
EPSRAD = 1e-6
EPSBEND = 1e-3
PSQMAX = 0.9999
XTHASX = (2.0**-52 * 13312.0 / 231.0 / 6.0)**0.1     # tfloor.f xthasx

# PH の並び
P0, H0, ANRAD, CUC, CRAD, BRHOZ, EMIDIV, EMIDIB, EMIDIQ, DVFS, PGEV = range(11)
RSGN = 11            # tradke の ddpx の符号 (K64 = +1、v2 = −1。v2next A1 項目 7: K64 の追跡 tradkn と逆)
EXACTD = 12          # tradke の D・B の更新に終点の座標に関する微分を使う (v2 = 1。v2next P2 項目 4)
QBSCH = 13           # 四極の本体の積分の方式 (quadbody.scheme(): 0 = K64 の Strang 分割、1 = 4 次の合成。B226)
# S の並び
U0, H1EMIT, DVEMIT = range(3)


def phys_array(ph):
    """emit6d._Phys → PH 配列"""
    from . import compat as _cmp
    return np.array([ph.p0, ph.h0, ph.anrad, ph.cuc, ph.crad, ph.brhoz, ph.emidiv, ph.emidib, ph.emidiq,
                     ph.dvfs, ph.pgev, -1.0 if _cmp.v2() else 1.0, 1.0 if _cmp.v2() else 0.0, float(_qb.scheme())], dtype=float)


# ---------------------------------------------------------------------------
# 数学関数 (mathfun)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _p2h(p):
    if p > 1e3:
        p2 = 1.0 / p**2
        return p * (1.0 + p2 * (0.5 - p2 * 0.125))
    return math.hypot(1.0, p)


@njit(**NJ)
def _pxy2dpz(px, py):
    """emit6d._pxy2dpz: sqrt(1-px^2-py^2)-1"""
    x = px**2 + py**2
    if x < 8e-6:
        return -x * (0.5 + x * (0.125 + x * 0.0625))
    return -x / (1.0 + math.sqrt(max(1e-200, 1.0 - x)))


@njit(**NJ)
def _pxy2dpz3(pxi, pyi, p):
    """lattice._pxy2dpz: sqrt(p^2 - pxi^2 - pyi^2) - p"""
    arg = p**2 - pxi**2 - pyi**2
    return math.sqrt(max(0.0, arg)) - p


@njit(**NJ)
def _asinz(x):
    return math.asin(max(-1.0, min(1.0, x)))


@njit(**NJ)
def _asinx(x):
    """mathfun asinx = asin(x) − x"""
    if abs(x) < XTHASX:
        x2 = x * x
        return x * x2 * (1.0 / 6.0 + x2 * (3.0 / 40.0 + x2 * (5.0 / 112.0 + x2 * (35.0 / 1152.0 + x2 * 63.0 / 2816.0))))
    return math.asin(min(1.0, max(-1.0, x))) - x


@njit(**NJ)
def _xsincos(x):
    """mathfun xsincos: (sin, x−sin, cos, cos−1)"""
    sn = math.sin(x)
    return sn, x - sn, math.cos(x), -2.0 * math.sin(0.5 * x)**2


@njit(**NJ)
def _dv(PH, dp):
    """_Phys.dv: dvemit = -(p1+p0)/h1/(p0 h1+p1 h0) dp + dvfs"""
    p0 = PH[P0]
    p1 = p0 * (1.0 + dp)
    h1 = _p2h(p1)
    return -(p1 + p0) / h1 / (p0 * h1 + p1 * PH[H0]) * dp + PH[DVFS]


@njit(**NJ)
def _dvz(PH, dp):
    h1 = _p2h(PH[P0] * (1.0 + dp))
    return PH[H0] / (h1 * h1 * h1)


@njit(**NJ)
def _mm(a, b):
    c = np.empty((6, 6))
    for i in range(6):
        for j in range(6):
            s = 0.0
            for k in range(6):
                s += a[i, k] * b[k, j]
            c[i, j] = s
    return c


@njit(**NJ)
def _piece(J, T, D, B, want_jac, want_db, symp_fix):
    """_State.piece のヤコビアン部分: z 行の組み直し (symp_fix) → apply"""
    if not want_jac:
        return
    if symp_fix:
        for j in range(4):
            J[4, j] = (-J[0, j] * J[1, 5] + J[1, j] * J[0, 5]
                       - J[2, j] * J[3, 5] + J[3, j] * J[2, 5])
    _apply_tdb(J, T, D, B, want_db)


@njit(**NJ)
def _setr0(cod, T, codr0, transr):
    for i in range(6):
        codr0[i] = cod[i]
        for j in range(6):
            transr[i, j] = T[i, j]


# ---------------------------------------------------------------------------
# tradke (photons.f, bzh=0, calpol=F)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _tradke(cod, T, D, B, codr0, transr, S, PH, al, phir0, want_jac, want_db, radcod, tin=0.0, tout=0.0):
    if not want_jac and not radcod:
        return
    gi = codr0[5]
    pr = 1.0 + gi
    th = math.tan(0.5 * phir0)
    sp = 2.0 * th / (1.0 + th**2)
    cp = 1.0 - th * sp
    pxi = codr0[1]
    pyi = codr0[3]
    pz00 = pr * (1.0 + _pxy2dpz(pxi / pr, pyi / pr))
    pxr0 = cp * pxi + sp * pz00
    pz0 = -sp * pxi + cp * pz00
    px = cod[1]
    py = cod[3]
    pz = pr * (1.0 + _pxy2dpz(px / pr, py / pr))
    dpx = px - pxr0
    dpy = py - pyi
    xpx = (py * pz0 - pz * pyi)
    xpy = (pz * pxr0 - px * pz0)
    xpz = (px * pyi - py * pxr0)
    xpa = math.sqrt(xpx**2 + xpy**2 + xpz**2) / pr**2
    theta = _asinz(xpa)
    p = PH[P0] * pr
    h1 = _p2h(p)
    al1 = al - cod[4] + codr0[4]
    gin = np.zeros(6)
    gout = np.zeros(6)
    if tin != 0.0 or tout != 0.0:
        # v2 (B220): 端の面までの直進の行程を厳密に (perfkern.edge_free。tan e → 0 で線形の式 tin·x に一致)
        l1, gin = _edge_free(tin, 0.0, codr0[0], codr0[2], codr0[1], codr0[3], codr0[5], -1.0)
        l2, gout = _edge_free(tout, 0.0, cod[0], cod[2], cod[1], cod[3], cod[5], 1.0)
        al1 = al1 - l1 - l2
    anrad = abs(PH[ANRAD])
    anp = anrad * h1 * theta
    uct = PH[CUC] * h1**3 / PH[P0] / al1
    dla = 0.0
    dlb = 0.0
    if PH[ANRAD] < 0.0:
        # v2 (B167): 厳密な放射のスペクトル (β < 1)。光子数に a、光子エネルギーに b (radbeta.ab)
        fa, fb, dla, dlb = _rb_ab(h1)
        anp *= fa
        uct *= fb
        anrad *= fa
    uc = uct * theta
    dg = -CAVE * anp * uc
    S[U0] -= dg
    ddpx = PH[RSGN] * 0.5 * dpx * dg      # v2: 符号を反転 (K64 の追跡 tradkn は ddpx = −0.5 dpx dg)
    ddpy = PH[RSGN] * 0.5 * dpy * dg
    c1 = al1 / pr / 3.0
    if radcod:
        g = max(GMIN, gi + dg)
        pr1 = 1.0 + g
        cod[0] += ddpx * c1
        cod[2] += ddpy * c1
        cod[1] = px * pr1 / pr + ddpx
        cod[3] = py * pr1 / pr + ddpy
        cod[5] = g
        p2 = PH[P0] * pr1
        h2 = _p2h(p2)
        cod[4] = cod[4] * p2 / h2 * h1 / p
        S[H1EMIT] = _p2h(p2)
        S[DVEMIT] = _dv(PH, g)
    if not want_db:
        _setr0(cod, T, codr0, transr)
        return
    ok, Dn, Bn = _tradke_db(T, transr, D, B, pr, px, py, pxi, pyi, pz, pz00, pz0, pxr0, cp, sp,
                            xpx, xpy, xpz, xpa, anp, uc, uct, al1, dg, ddpx, ddpy, c1, h1, p, PH[P0], anrad, CAVE, CUU,
                            0.0, 0.0, tin, tout, False, PH[RSGN], PH[EXACTD] != 0.0, 0.0, 0.0, dla, dlb, gin, gout)
    if ok:
        for i in range(6):
            for j in range(6):
                D[i, j] = Dn[i, j]
                B[i, j] = Bn[i, j]
    _setr0(cod, T, codr0, transr)


# ---------------------------------------------------------------------------
# ドリフト・薄い多極 (tdrife0, tthine)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _jac_drift(c, L, dvzL):
    """emit6d._jac_drift。戻り (ok, J)"""
    pr = 1.0 + c[5]
    px = c[1]
    py = c[3]
    pz2 = pr * pr - px * px - py * py
    J = np.eye(6)
    if pz2 <= 1e-12:
        return False, J
    k = L / (pz2 * math.sqrt(pz2))
    J[0, 1] = k * (pz2 + px * px)
    J[0, 3] = k * px * py
    J[0, 5] = -k * px * pr
    J[2, 1] = J[0, 3]
    J[2, 3] = k * (pz2 + py * py)
    J[2, 5] = -k * py * pr
    J[4, 1] = -k * pr * px
    J[4, 3] = -k * pr * py
    J[4, 5] = k * (px * px + py * py) + dvzL
    return True, J


@njit(**NJ)
def _exact_drift(c, L):
    """lattice._exact_drift_cod (in-place)"""
    pr = 1.0 + c[5]
    a = c[1]**2 + c[3]**2
    pzi = pr * math.sqrt(max(0.0, 1.0 - a / pr**2))
    c[0] += c[1] / pzi * L
    c[2] += c[3] / pzi * L
    c[4] -= a / (pr + pzi) / pzi * L


@njit(**NJ)
def _drift_map(c, L, PH, use_dv):
    """emit6d._drift_map (in-place): 厳密ドリフト + dv 項"""
    pr = 1.0 + c[5]
    a = c[1]**2 + c[3]**2
    dpz = -a / pr / (1.0 + math.sqrt(max(0.0, 1.0 - a / pr**2)))
    pzi = pr + dpz
    ale = L / pzi
    c[0] += c[1] * ale
    c[2] += c[3] * ale
    c[4] += dpz * ale
    if use_dv:
        c[4] -= _dv(PH, c[5]) * L


@njit(**NJ)
def _cpowi(a, b, n):
    """(a + i b)^n (n ≥ 0 の整数)。CPython の complex ** int (c_powu: 2 進べき乗) と同じ掛け算の順。
    numba の complex ** int は exp/log 経由で b=0 でも 1e-38 程度の虚部の雑音を作るので使わない"""
    rr = 1.0
    ri = 0.0
    pr = a
    pi = b
    mask = 1
    while mask > 0 and n >= mask:
        if n & mask:
            rr, ri = rr * pr - ri * pi, rr * pi + ri * pr
        mask <<= 1
        pr, pi = pr * pr - pi * pi, pr * pi + pi * pr
    return rr, ri


@njit(**NJ)
def _jac_kick(c, aki, kord):
    gr, gi = _cpowi(c[0], -c[2], kord - 1)
    gr *= kord
    gi *= kord
    J = np.eye(6)
    J[1, 0] = -aki * gr
    J[1, 2] = -aki * gi
    J[3, 0] = -aki * gi
    J[3, 2] = aki * gr
    return J


@njit(**NJ)
def _kick(c, a, kord):
    cr, ci = _cpowi(c[0], -c[2], kord)
    c[1] -= a * cr
    c[3] -= a * ci


@njit(**NJ)
def _factorial(k):
    f = 1.0
    for i in range(2, k + 1):
        f *= i
    return f


_TT_ALPHA = 2.88675134594812882e-1       # 4 次の力勾配型の分割の係数 (lattice._TT_*、K64 tquad.f tthin)
_TT_BETA = 2.23290993692602255e-2
_TT_GAMMA = 9.43738783765593145e-4
_TT_ALPHA1 = 0.5 - _TT_ALPHA


@njit(**NJ)
def _thin4_corr(x, y, pr, aki, kord, L, want_jac):
    """lattice._thin4_corr (v2 の B326: 4 次の分割の中央の補正キック)。戻り (dp (3,) = (Δpx, Δpy, Δz)、Jd (3, 3) = (x, y, δ) での微分)"""
    k = kord
    f1 = _TT_BETA * k * L
    f2 = 2.0 * _TT_GAMMA * k * (k - 1) * L**2
    f3 = _TT_GAMMA * k * (k + 1) * L**2
    f4 = 0.5 * _TT_BETA * L
    f5 = 2.0 * _TT_GAMMA * k * L**2
    r = x * x + y * y
    rk1 = 1.0
    for _ in range(k - 2):
        rk1 *= r
    rk = r * rk1
    cr, ci = _cpowi(x, -y, k)
    R = x * cr + y * ci
    B1 = aki * aki / pr
    B2 = aki**3 / pr**2
    Sx = f2 * rk1 * x * R + f3 * rk * cr
    Sy = f2 * rk1 * y * R + f3 * rk * ci
    d3 = aki**2 / pr**2
    d4 = aki**3 / pr**3
    dp = np.zeros(3)
    Jd = np.zeros((3, 3))
    dp[0] = B1 * f1 * rk * x - B2 * Sx
    dp[1] = B1 * f1 * rk * y - B2 * Sy
    dp[2] = d3 * f4 * rk * r - d4 * f5 * rk * R
    if not want_jac:
        return dp, Jd
    gr, gi = _cpowi(x, -y, k - 1)
    gr *= k
    gi *= k
    drk_x = 2.0 * (k - 1) * rk1 * x
    drk_y = 2.0 * (k - 1) * rk1 * y
    r3 = 0.0
    if k >= 3:
        r3 = 1.0
        for _ in range(k - 3):
            r3 *= r
    drk1_x = 2.0 * (k - 2) * r3 * x
    drk1_y = 2.0 * (k - 2) * r3 * y
    R_x = cr + x * gr + y * gi
    R_y = ci + x * gi - y * gr
    Sx_x = f2 * (drk1_x * x * R + rk1 * R + rk1 * x * R_x) + f3 * (drk_x * cr + rk * gr)
    Sx_y = f2 * (drk1_y * x * R + rk1 * x * R_y) + f3 * (drk_y * cr + rk * gi)
    Sy_x = f2 * (drk1_x * y * R + rk1 * y * R_x) + f3 * (drk_x * ci + rk * gi)
    Sy_y = f2 * (drk1_y * y * R + rk1 * R + rk1 * y * R_y) + f3 * (drk_y * ci - rk * gr)
    A = B1 * f1
    Jd[0, 0] = A * (drk_x * x + rk) - B2 * Sx_x
    Jd[0, 1] = A * drk_y * x - B2 * Sx_y
    Jd[0, 2] = -A / pr * rk * x + 2.0 * B2 / pr * Sx
    Jd[1, 0] = A * drk_x * y - B2 * Sy_x
    Jd[1, 1] = A * (drk_y * y + rk) - B2 * Sy_y
    Jd[1, 2] = -A / pr * rk * y + 2.0 * B2 / pr * Sy
    Jd[2, 0] = d3 * f4 * (drk_x * r + rk * 2.0 * x) - d4 * f5 * (drk_x * R + rk * R_x)
    Jd[2, 1] = d3 * f4 * (drk_y * r + rk * 2.0 * y) - d4 * f5 * (drk_y * R + rk * R_y)
    Jd[2, 2] = -2.0 * d3 / pr * f4 * rk * r + 3.0 * d4 / pr * f5 * rk * R
    return dp, Jd


@njit(**NJ)
def _jac_thin4(c, L, aki, kord, dvzL):
    """emit6d._jac_thin4 (v2 の B326: 4 次の分割の COD まわりの微分)。戻り (ok, J)"""
    ala = L * _TT_ALPHA1
    alb = L * _TT_ALPHA
    ok, J = _jac_drift(c, ala, 0.0)
    if not ok:
        return False, J
    c1 = c.copy()
    _exact_drift(c1, ala)
    J = _mm(_jac_kick(c1, aki * 0.5, kord), J)
    _kick(c1, aki * 0.5, kord)
    ok, J2 = _jac_drift(c1, alb, 0.0)
    if not ok:
        return False, J
    J = _mm(J2, J)
    _exact_drift(c1, alb)
    dp, Jd = _thin4_corr(c1[0], c1[2], 1.0 + c1[5], aki, kord, L, True)
    Jc = np.eye(6)
    for i in range(3):
        row = 1 if i == 0 else (3 if i == 1 else 4)
        Jc[row, 0] += Jd[i, 0]
        Jc[row, 2] += Jd[i, 1]
        Jc[row, 5] += Jd[i, 2]
    J = _mm(Jc, J)
    c1[1] += dp[0]
    c1[3] += dp[1]
    c1[4] += dp[2]
    ok, J3 = _jac_drift(c1, alb, 0.0)
    if not ok:
        return False, J
    J = _mm(J3, J)
    _exact_drift(c1, alb)
    J = _mm(_jac_kick(c1, aki * 0.5, kord), J)
    _kick(c1, aki * 0.5, kord)
    ok, J4 = _jac_drift(c1, ala, 0.0)
    if not ok:
        return False, J
    J = _mm(J4, J)
    J[4, 5] += dvzL
    return True, J


@njit(**NJ)
def _jac_thin(c, L, ak, kord, dvzL, v2=False):
    """emit6d._jac_thin。戻り (ok, J)。v2=True: L≠0・kord ≥ 2 は 4 次の分割 (B326)"""
    if ak == 0.0:
        return _jac_drift(c, L, dvzL)
    aki = ak / _factorial(kord)
    if L == 0.0:
        return True, _jac_kick(c, aki, kord)
    if v2 and kord >= 2:
        return _jac_thin4(c, L, aki, kord, dvzL)
    ok, J = _jac_drift(c, L / 6.0, 0.0)
    if not ok:
        return False, J
    c1 = c.copy()
    _exact_drift(c1, L / 6.0)
    J = _mm(_jac_kick(c1, aki * 0.5, kord), J)
    _kick(c1, aki * 0.5, kord)
    ok, J2 = _jac_drift(c1, L / 1.5, 0.0)
    if not ok:
        return False, J
    J = _mm(J2, J)
    _exact_drift(c1, L / 1.5)
    J = _mm(_jac_kick(c1, aki * 0.5, kord), J)
    _kick(c1, aki * 0.5, kord)
    ok, J3 = _jac_drift(c1, L / 6.0, 0.0)
    if not ok:
        return False, J
    J = _mm(J3, J)
    J[4, 5] += dvzL
    return True, J


@njit(**NJ)
def _tthine_map(c, L, ak, kord, v2=False):
    """lattice._tthine_cod (in-place)。v2=True: L≠0・kord ≥ 2 は 4 次の分割 (B326)"""
    if ak == 0.0:
        _exact_drift(c, L)
        return
    aki = ak / _factorial(kord)
    if L == 0.0:
        _kick(c, aki, kord)
        return
    if v2 and kord >= 2:
        _exact_drift(c, L * _TT_ALPHA1)
        _kick(c, aki * 0.5, kord)
        _exact_drift(c, L * _TT_ALPHA)
        dp, Jd = _thin4_corr(c[0], c[2], 1.0 + c[5], aki, kord, L, False)
        c[1] += dp[0]
        c[3] += dp[1]
        c[4] += dp[2]
        _exact_drift(c, L * _TT_ALPHA)
        _kick(c, aki * 0.5, kord)
        _exact_drift(c, L * _TT_ALPHA1)
        return
    _exact_drift(c, L / 6.0)
    _kick(c, aki * 0.5, kord)
    _exact_drift(c, L / 1.5)
    _kick(c, aki * 0.5, kord)
    _exact_drift(c, L / 6.0)


@njit(**NJ)
def _rot4(c, ct, st):
    """(x,px,y,py) を θ で回す: x' = c x − s y, y' = s x + c y (px, py も)"""
    x = c[0]
    px = c[1]
    c[0] = ct * x - st * c[2]
    c[2] = st * x + ct * c[2]
    c[1] = ct * px - st * c[3]
    c[3] = st * px + ct * c[3]


@njit(**NJ)
def _rotmat(ct, st, ent):
    """_tchge_jac の回転行列 (ent: x' = c x − s y、出口: 逆)"""
    J = np.eye(6)
    if ent:
        J[0, 0] = ct; J[0, 2] = -st; J[1, 1] = ct; J[1, 3] = -st
        J[2, 0] = st; J[2, 2] = ct; J[3, 1] = st; J[3, 3] = ct
    else:
        J[0, 0] = ct; J[0, 2] = st; J[1, 1] = ct; J[1, 3] = st
        J[2, 0] = -st; J[2, 2] = ct; J[3, 1] = -st; J[3, 3] = ct
    return J


@njit(**NJ)
def _quirk56(c, v2=False):
    """tbgrote (dr=0) の trans: 単位行列で (5,6) だけ dpzi·(pr+1)/pzi (lattice._tchge_jac の simple 分岐)"""
    J = np.eye(6)
    if v2:
        return J            # v2 (B109): dr=0 の写像は恒等なので (5,6) = 0
    pr = 1.0 + c[5]
    a = (c[1] / pr)**2 + (c[3] / pr)**2
    dpzi = (math.sqrt(1.0 - a) - 1.0) * pr
    J[4, 5] = dpzi * (pr + 1.0) / (pr + dpzi)
    return J


@njit(**NJ)
def _tchge_simple(c, T, D, B, dx0, dy, th, ent, want_jac, want_db, symp_fix, v2=False):
    """tchge (alg=phig=dz=dtheta=0): 薄い多極・四極のずれ dx0,dy とロール th。lattice._tchge_cod/_tchge_jac と同じ
    (入口: 回す → ずらす。th≠0 かつ dxr≠0 のとき行列には回転が入らない SAD の癖。出口: 戻す → 回転)"""
    if th != 0.0:
        ct = math.cos(th)
        st = math.sin(th)
        dxr = ct * dx0 - st * dy
    else:
        ct = 1.0
        st = 0.0
        dxr = dx0
    if ent:
        if th != 0.0:
            _rot4(c, ct, st)
            c[2] -= st * dx0 + ct * dy
            if dxr != 0.0:
                J = _quirk56(c, v2)
                if v2:
                    J = _mm(J, _rotmat(ct, st, True))      # v2 (B09): 入口の回転を行列に掛ける
            else:
                J = _rotmat(ct, st, True)
        else:
            c[2] -= dy
            J = np.eye(6)
        if dxr != 0.0:
            c[0] -= dxr
        _piece(J, T, D, B, want_jac, want_db, symp_fix)
    else:
        if th != 0.0:
            J = _rotmat(ct, st, False)
        else:
            J = np.eye(6)
        if dxr != 0.0:
            c[0] += dxr
        if th != 0.0:
            c[2] += st * dx0 + ct * dy
            _rot4(c, ct, -st)
        else:
            c[2] += dy
        _piece(J, T, D, B, want_jac, want_db, symp_fix)


@njit(**NJ)
def k_drift(cod, T, D, B, S, PH, L, use_dv, symp_fix, want_jac, want_db):
    """tdrife0。戻り 0: 成功、1: 解析ヤコビアン不可 (Python に戻す)"""
    if L == 0.0:
        return 0
    if want_jac:
        dvzL = _dvz(PH, cod[5]) * L if use_dv else 0.0
        ok, J = _jac_drift(cod, L, dvzL)
        if not ok:
            return 1
        _piece(J, T, D, B, want_jac, want_db, symp_fix)
    _drift_map(cod, L, PH, use_dv)
    return 0


@njit(**NJ)
def _fringe_raw(x, y, p, q, aki, n):
    """多極 (nord = 2n、n ≥ 3) の非線形フリンジ ttfrin (track6._ttfrin_gen・_ttfrin6 と同じ式) を実数の演算だけで書いたもの。
    (x, y, p = px/(1+δ), q = py/(1+δ)) → (x1, y1, p1, q1, Δz)。+・−・×・÷ だけなので複素ステップ微分 (虚部の 1e-30) で厳密な導関数が取れる。
    aki = ak/al/n!/4/(1+δ) (n = 3 で 6 極: ak/al/24)"""
    u1 = x * 0.0 + 1.0
    v1 = x * 0.0
    for _i in range(n - 1):
        u1, v1 = u1 * x - v1 * y, u1 * y + v1 * x          # (x + i y)^(n−1)
    u = u1 * x - v1 * y                                     # cz = cz1·cx
    v = u1 * y + v1 * x
    a = 2.0 * aki * v
    ar = u / (n + 1.0) - u                                  # cz/(n+1) − conj(cz)
    ai = v / (n + 1.0) + v
    car = -aki * (x * ar - y * ai)                          # ca = −aki·cx·(…)
    cai = -aki * (x * ai + y * ar)
    d = 1.0 + a * a - (aki * n)**2 * (u * u + v * v)
    m = x * p + y * q                                       # cx·cp (cp = p − i q)
    nn = y * p - x * q
    cp1r = ((p + a * q) - aki * n * (u1 * m + v1 * nn)) / d
    cp1i = ((a * p - q) - aki * n * (v1 * m - u1 * nn)) / d
    return x + car, y + cai, cp1r, -cp1i, -(car * cp1r - cai * cp1i)


@njit(**NJ)
def _fringe_step(c, ak, al, n, want_jac):
    """多極の非線形フリンジ ttfrin (入口は +ak、出口は −ak) の写像 (c を in-place で更新) と解析ヤコビアン (複素ステップ微分 + 連鎖)。
    c は CALC・Emittance[] の正準運動量 (px は p0 で規格化)、追跡は px/(1+δ) なので p = px/(1+δ)・aki ∝ 1/(1+δ) に直して同じ式を使う。
    戻り: J (want_jac のとき。そうでなければ単位行列)"""
    J = np.eye(6)
    if ak == 0.0 or al == 0.0:
        return J
    pr = 1.0 + c[5]
    x = c[0]
    y = c[2]
    p = c[1] / pr
    q = c[3] / pr
    aki = ak / al / _factorial(n) / 4.0 / pr
    r = _fringe_raw(x, y, p, q, aki, n)
    if want_jac:
        h = 1e-30
        P = np.zeros((5, 5))
        for k in range(5):
            v = np.zeros(5, dtype=np.complex128)
            v[0] = x
            v[1] = y
            v[2] = p
            v[3] = q
            v[4] = aki
            v[k] += 1j * h
            rc = _fringe_raw(v[0], v[1], v[2], v[3], v[4], n)
            P[0, k] = rc[0].imag / h
            P[1, k] = rc[1].imag / h
            P[2, k] = rc[2].imag / h
            P[3, k] = rc[3].imag / h
            P[4, k] = rc[4].imag / h
        dpd = -p / pr                    # ∂p/∂δ (px 一定)
        dqd = -q / pr
        dad = -aki / pr
        for o in range(3):               # 出力 x1 (行 0)、y1 (行 2)、z (行 4): P の行 0、1、4
            po = (0, 1, 4)[o]
            ro = (0, 2, 4)[o]
            J[ro, 0] = P[po, 0]
            J[ro, 2] = P[po, 1]
            J[ro, 1] = P[po, 2] / pr
            J[ro, 3] = P[po, 3] / pr
            J[ro, 5] = P[po, 2] * dpd + P[po, 3] * dqd + P[po, 4] * dad
        J[4, 4] = 1.0
        J[1, 0] = pr * P[2, 0]
        J[1, 2] = pr * P[2, 1]
        J[1, 1] = P[2, 2]
        J[1, 3] = P[2, 3]
        J[1, 5] = r[2] + pr * (P[2, 2] * dpd + P[2, 3] * dqd + P[2, 4] * dad)
        J[3, 0] = pr * P[3, 0]
        J[3, 2] = pr * P[3, 1]
        J[3, 1] = P[3, 2]
        J[3, 3] = P[3, 3]
        J[3, 5] = r[3] + pr * (P[3, 2] * dpd + P[3, 3] * dqd + P[3, 4] * dad)
    c[0] = r[0]
    c[2] = r[1]
    c[1] = r[2] * pr
    c[3] = r[3] * pr
    c[4] += r[4]
    return J


@njit(**NJ)
def k_thin(cod, T, D, B, S, PH, L, ak, kord, dx, dy, th, use_dv, symp_fix, want_jac, want_db, v2=False, fringe=False):
    """tthine (emit6d._tthine_gen): tchge(入口) → [フリンジ] → ドリフト/キック + dv → [フリンジ] → tchge(出口)。状態の複製で計算し成功時だけ書き戻す。
    fringe=True (v2 の L≠0 の 6 極以上。追跡の ttfrin と同じ非線形フリンジをヤコビアンと COD の写像に入れる。K64 の tthine は入れない): n = kord + 1"""
    c = cod.copy()
    Tl = T.copy()
    Dl = D.copy()
    Bl = B.copy()
    shift = ak != 0.0 and (dx != 0.0 or dy != 0.0 or th != 0.0)
    fr = fringe and ak != 0.0 and L != 0.0
    if shift:
        _tchge_simple(c, Tl, Dl, Bl, dx, dy, th, True, want_jac, want_db, symp_fix, v2)
    if fr:
        Jf = _fringe_step(c, ak, L, kord + 1, want_jac)
        _piece(Jf, Tl, Dl, Bl, want_jac, want_db, symp_fix)
    if want_jac:
        dvzL = _dvz(PH, c[5]) * L if use_dv else 0.0
        ok, J = _jac_thin(c, L, ak, kord, dvzL, v2)
        if not ok:
            return 1
        _piece(J, Tl, Dl, Bl, want_jac, want_db, symp_fix)
    _tthine_map(c, L, ak, kord, v2)
    if use_dv:
        c[4] -= _dv(PH, c[5]) * L
    if fr:
        Jf = _fringe_step(c, -ak, L, kord + 1, want_jac)
        _piece(Jf, Tl, Dl, Bl, want_jac, want_db, symp_fix)
    if shift:
        _tchge_simple(c, Tl, Dl, Bl, dx, dy, th, False, want_jac, want_db, symp_fix, v2)
    cod[:] = c
    T[:, :] = Tl
    D[:, :] = Dl
    B[:, :] = Bl
    return 0


# ---------------------------------------------------------------------------
# tchge (一般: BEND の据え付け誤差。lattice._tchge_cod / _tbrote_jac / _tbgrote_mat / _tdrife1)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _tdrife1(c, dv, al):
    """toplvl.f tdrife1: (trans1, cod 更新 in-place)"""
    pr = 1.0 + c[5]
    pxi = c[1]
    pyi = c[3]
    a = pxi**2 + pyi**2
    pzi = pr * math.sqrt(1.0 - a / pr**2)
    ale = al / pzi
    alz = ale / pzi**2
    t = np.eye(6)
    t[0, 1] = ale + pxi**2 * alz
    t[0, 3] = pxi * pyi * alz
    t[0, 5] = -pxi * pr * alz
    t[2, 1] = t[0, 3]
    t[2, 3] = ale + pyi**2 * alz
    t[2, 5] = -pyi * pr * alz
    t[4, 1] = t[0, 5]
    t[4, 3] = t[2, 5]
    t[4, 5] = a * alz
    c[0] += pxi / pzi * al
    c[2] += pyi / pzi * al
    c[4] -= (a / (pr + pzi) / pzi + dv) * al
    return t


@njit(**NJ)
def _tbgrote(c, dr, dx, dz, v2=False):
    """tchge.f tbgrote: cod 更新 (in-place) と trans"""
    xi = c[0] + dx
    yi = c[2]
    pxi = c[1]
    pyi = c[3]
    pr = 1.0 + c[5]
    a = (pxi / pr)**2 + (pyi / pr)**2
    dpzi = -a / (1.0 + math.sqrt(max(1e-200, 1.0 - a))) * pr
    pzi = pr + dpzi
    zi = dz
    xf = dr[0, 0] * xi + dr[0, 1] * yi + dr[0, 2] * zi + c[0]
    yf = dr[1, 0] * xi + dr[1, 1] * yi + dr[1, 2] * zi + yi
    zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * zi
    r11 = 1.0 + dr[0, 0]
    r22 = 1.0 + dr[1, 1]
    r33 = 1.0 + dr[2, 2]
    pxf = r11 * pxi + dr[0, 1] * pyi + dr[0, 2] * pzi
    pyf = dr[1, 0] * pxi + r22 * pyi + dr[1, 2] * pzi
    dpzf = dr[2, 0] * pxi + dr[2, 1] * pyi + r33 * dpzi + dr[2, 2] * pr
    pzf = pr + dpzf
    pxfpzf = pxf / pzf
    pyfpzf = pyf / pzf
    zfpzf = zf / pzf
    prpzf = pr / pzf
    c[1] = pxf
    c[3] = pyf
    c[0] = xf - pxf * zfpzf
    c[2] = yf - pyf * zfpzf
    c[4] = c[4] + pr * zfpzf
    dpzidpx = -pxi / pzi
    dpzidpy = -pyi / pzi
    dpzidp = pr / pzi
    t = np.zeros((6, 6))
    t[1, 1] = r11 + dr[0, 2] * dpzidpx
    t[1, 3] = dr[0, 1] + dr[0, 2] * dpzidpy
    t[1, 5] = dr[0, 2] * dpzidp
    t[3, 1] = dr[1, 0] + dr[1, 2] * dpzidpx
    t[3, 3] = r22 + dr[1, 2] * dpzidpy
    t[3, 5] = dr[1, 2] * dpzidp
    t[0, 0] = r11 - dr[2, 0] * pxfpzf
    t[0, 1] = -zfpzf * (t[1, 1] - pxfpzf * (dr[2, 0] + r33 * dpzidpx))
    t[0, 2] = dr[0, 1] - dr[2, 1] * pxfpzf
    t[0, 3] = -zfpzf * (t[1, 3] - pxfpzf * (dr[2, 1] + r33 * dpzidpy))
    t[0, 5] = -zfpzf * (dr[0, 2] - r33 * pxfpzf) * dpzidp
    t[2, 0] = dr[1, 0] - dr[2, 0] * pyfpzf
    t[2, 1] = -zfpzf * (t[3, 1] - pyfpzf * (dr[2, 0] + r33 * dpzidpx))
    t[2, 2] = r22 - dr[2, 1] * pyfpzf
    t[2, 3] = -zfpzf * (t[3, 3] - pyfpzf * (dr[2, 1] + r33 * dpzidpy))
    t[2, 5] = -zfpzf * (dr[1, 2] - r33 * pyfpzf) * dpzidp
    t[4, 0] = dr[2, 0] * prpzf
    t[4, 1] = -zfpzf * prpzf * (dr[2, 0] + r33 * dpzidpx)
    t[4, 2] = dr[2, 1] * prpzf
    t[4, 3] = -zfpzf * prpzf * (dr[2, 1] + r33 * dpzidpy)
    t[4, 4] = 1.0
    if v2:
        t[4, 5] = zfpzf * (1.0 - prpzf * r33 * dpzidp)      # v2 (B109): 写像の厳密な微分 (K64 のコメントアウトされた式)
    else:
        t[4, 5] = r33 * (dpzf * prpzf + dpzi / pzi) - dr[2, 2]
    t[5, 5] = 1.0
    return t


@njit(**NJ)
def _tbrote(c, alg, phig, dx, dz, dtheta, ent, v2=False):
    """tchge.f tbrote (dchi2=0): cod 更新 (in-place) と trans (lattice._tbrote_jac)"""
    if phig == 0.0:
        c0 = 1.0
        s0 = 0.0
        dl = alg
    else:
        c0 = math.cos(phig)
        s0 = math.sin(phig)
        dl = alg * s0 / phig
    s2 = 0.0
    c2 = 1.0
    dc2 = 0.0
    st, _xst, ct, dct = _xsincos(dtheta)
    rr = np.zeros((3, 3))
    rr[0, 0] = dct * c0**2 + dc2 * s0**2 + c0 * s0 * s2 * st
    rr[1, 1] = c2 * dct + dc2
    rr[2, 2] = dc2 * c0**2 + dct * s0**2 - c0 * s0 * s2 * st
    if ent:
        dx1 = dx * c0 - dz * s0
        dz1 = dz * c0 + dx * s0
        rr[0, 1] = s0 * s2 - c0 * c2 * st
        rr[0, 2] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[1, 0] = -ct * s0 * s2 + c0 * st
        rr[1, 2] = c0 * ct * s2 + s0 * st
        rr[2, 0] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[2, 1] = -c0 * s2 - c2 * s0 * st
        c[0] -= dx1
        if dz1 != 0.0:
            J1 = _tdrife1(c, 1.0, dz1)
            J2 = _tbgrote(c, rr, dl * s0, -dl * c0, v2)
            return _mm(J2, J1)
        return _tbgrote(c, rr, dl * s0, -dl * c0, v2)
    rr[1, 0] = s0 * s2 - c0 * c2 * st
    rr[2, 0] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
    rr[0, 1] = -ct * s0 * s2 + c0 * st
    rr[2, 1] = c0 * ct * s2 + s0 * st
    rr[0, 2] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
    rr[1, 2] = -c0 * s2 - c2 * s0 * st
    J = _tbgrote(c, rr, dl * s0, -dl * c0, v2)
    if v2:
        dx1 = dx * c0 - dz * s0         # v2 (A2): lattice._tbrote_cod と同じ (出口の (dx, dz) の回し方の符号)
        dz1 = dz * c0 + dx * s0
    else:
        dx1 = dx * c0 + dz * s0
        dz1 = dz * c0 - dx * s0
    if dz1 != 0.0:
        J2 = _tdrife1(c, 1.0, -dz1)
        J = _mm(J2, J)
    c[0] += dx1
    return J


@njit(**NJ)
def _tchge(c, T, D, B, dx0, dy, dz, theta, dtheta, alg, phig, ent, want_jac, want_db, symp_fix, v2=False):
    """lattice._tchge_cod + _tchge_jac を 1 つの piece として適用 (in-place)"""
    th = theta if phig != 0.0 else theta + dtheta
    if th != 0.0:
        ct = math.cos(th)
        st = math.sin(th)
        dxr = ct * dx0 - st * dy
    else:
        ct = 1.0
        st = 0.0
        dxr = dx0
    brote = (dxr != 0.0 or dz != 0.0 or dtheta != 0.0)
    simple = (dz == 0.0 and dtheta == 0.0 and (phig == 0.0 or dxr == 0.0))
    J = np.eye(6)
    if ent:
        if th != 0.0:
            _rot4(c, ct, st)
            c[2] -= st * dx0 + ct * dy
            if brote:
                if simple:
                    J = _quirk56(c, v2)
                    _tbrote(c, alg, phig, dxr, dz, dtheta, True, v2)
                else:
                    J = _tbrote(c, alg, phig, dxr, dz, dtheta, True, v2)
                if v2:
                    J = _mm(J, _rotmat(ct, st, True))      # v2 (B09): 入口の回転を行列に掛ける (K64 は tbrote が上書きして落とす)
            else:
                J = _rotmat(ct, st, True)
        else:
            c[2] -= dy
            if brote:
                if simple:
                    _tbrote(c, alg, phig, dxr, dz, dtheta, True, v2)
                else:
                    J = _tbrote(c, alg, phig, dxr, dz, dtheta, True, v2)
        _piece(J, T, D, B, want_jac, want_db, symp_fix)
        return
    if brote:
        if simple:
            _tbrote(c, alg, phig, dxr, dz, dtheta, False, v2)
        else:
            J = _tbrote(c, alg, phig, dxr, dz, dtheta, False, v2)
    if th != 0.0:
        c[2] += st * dx0 + ct * dy
        _rot4(c, ct, -st)
        J = _mm(_rotmat(ct, st, False), J)
    else:
        c[2] += dy
    _piece(J, T, D, B, want_jac, want_db, symp_fix)


# ---------------------------------------------------------------------------
# BEND のエッジ (tbfrie.f) と本体 (tbendebody0, tsteee)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _tbedgedrift(c, psi):
    """lattice._tbedgedrift_mat (更新前の cod で trans) と _tbedgedrift_cod (in-place)"""
    pr = 1.0 + c[5]
    sinp = math.sin(psi)
    cosp = math.cos(psi)
    pxi = min(pr, max(-pr, c[1]))
    pyi = min(pr, max(-pr, c[3]))
    a = pxi**2 + pyi**2
    dpzi = pr * (math.sqrt(max(0.0, 1.0 - a / pr**2)) - 1.0)
    s = a
    pzi = pr + dpzi
    f = cosp - pxi / pzi * sinp
    x1 = c[0] / f
    sx = x1 * sinp
    sxa = sx / f / pzi**3
    phsq = pzi**2 + pxi**2
    t = np.zeros((6, 6))
    t[0, 0] = 1.0 / f
    t[0, 1] = phsq * sxa
    t[0, 3] = pxi * pyi * sxa
    t[0, 5] = -pxi * pr * sxa
    t[1, 1] = f
    t[1, 3] = -pyi / pzi * sinp
    t[1, 5] = pr / pzi * sinp
    t[2, 0] = pyi / f / pzi * sinp
    t[2, 1] = pyi * sxa * (pxi * cosp + pzi * sinp)
    t[2, 2] = 1.0
    t[2, 3] = sxa * (f * pzi**2 + pyi**2 * cosp)
    t[2, 5] = -pr * pyi * sxa * (cosp / pzi)
    t[3, 3] = 1.0
    t[4, 0] = -t[0, 0] * t[1, 5]
    t[4, 1] = t[1, 1] * t[0, 5] - t[0, 1] * t[1, 5]
    t[4, 3] = t[1, 3] * t[0, 5] - t[0, 3] * t[1, 5] + t[2, 5]
    t[4, 4] = 1.0
    t[4, 5] = sxa * (f * pzi * s + pr**2 * pxi * sinp) / pzi
    t[5, 5] = 1.0
    # cod (lattice._tbedgedrift_cod: dpzi は _pxy2dpz(pxi, pyi, pr) = sqrt(max(0, pr²−a)) − pr)
    dpzi2 = _pxy2dpz3(pxi, pyi, pr)
    pzi2 = pr + dpzi2
    f2 = cosp - pxi / pzi2 * sinp
    c[0] = c[0] / f2
    sx2 = c[0] * sinp
    c[1] = pxi * cosp + pzi2 * sinp
    c[2] = c[2] + pyi / pzi2 * sx2
    c[4] = c[4] - pr / pzi2 * sx2
    return t


@njit(**NJ)
def _tbedgemaxwell(c, rhob, p2=False):
    """lattice._tbedgemaxwell_mat / _cod (in-place)。p2 = v2 (B324): COD の更新の Δx に bf·p²、Δz に −bf·px·p (k64 は p¹・p⁰)"""
    pr = 1.0 + c[5]
    pxi = min(pr, max(-pr, c[1]))
    pvi2 = (pr - pxi) * (pr + pxi)
    pvi = pr * math.sqrt(pvi2 / pr**2)
    yi = c[2]
    yr = (yi / rhob)**2
    a = 1.0 - yr / 6.0
    af = a * yi / rhob / pvi
    b = 0.5 - yr / 24.0
    bf = b * yi**2 / rhob / pvi2 / pvi
    cc = 1.0 - yr / 2.0
    t = np.eye(6)
    t[0, 1] = 3.0 * bf * pr**2 * pxi / pvi2
    t[0, 2] = af * pr**2 / pvi2
    t[0, 5] = -bf * pr * (pr**2 + 2.0 * pxi**2) / pvi2
    t[3, 1] = -t[0, 2]
    t[3, 2] = -cc * pxi / rhob / pvi
    t[3, 5] = af * pr * pxi / pvi2
    t[4, 1] = t[0, 5]
    t[4, 2] = -t[3, 5]
    t[4, 5] = bf * pxi * (2.0 * pr**2 + pxi**2) / pvi2
    # cod (lattice._tbedgemaxwell_cod)
    pvi2c = max(1e-30, pvi2)
    pvic = math.sqrt(pvi2c)
    afc = a * yi / rhob / pvic
    bfc = b * yi**2 / rhob / pvi2c / pvic
    pp = pr if p2 else 1.0
    c[0] = c[0] + bfc * pr * pp
    c[3] = min(pr, max(-pr, c[3] - afc * pxi))
    c[4] = c[4] - bfc * pxi * pp
    return t


@njit(**NJ)
def _tbedgebody(c, rhob, psi):
    """lattice._tbedgebody_mat / _cod (in-place)"""
    ampmin = 1e-6
    ampmax = 1.0 - ampmin
    pr = 1.0 + c[5]
    cosp = math.cos(psi)
    sinp = math.sin(psi)
    sinsq = 2.0 * math.sin(0.5 * psi)**2
    pxi = min(pr, max(-pr, c[1]))
    pyi = min(pr, max(-pr, c[3]))
    s = min(ampmax, pxi**2 + pyi**2)
    sqrs = pr * math.sqrt(max(ampmin, 1.0 - s / pr**2))
    dpzi = -s / (pr + sqrs)
    pzi = pr + dpzi
    phsq = pzi**2 + pxi**2
    px0 = pxi * cosp - pzi * sinp
    dpz0 = pxi * sinp + dpzi * cosp - pr * sinsq
    pz0 = pr + dpz0
    xi = c[0]
    pxf = px0 + xi * sinp / rhob
    s1 = min(ampmax, pyi**2 + pxf**2)
    sqrs1 = pr * math.sqrt(max(ampmin, 1.0 - s1 / pr**2))
    dpzf = -s1 / (pr + sqrs1)
    pzf = pr + dpzf
    ddpz = (s - s1) / (sqrs + sqrs1) + (dpzi + pr) * sinsq - pxi * sinp
    a = pr * xi * sinp / rhob + dpz0 * pxf - dpzf * px0
    aa = a * cosp - (px0 * pxf + pz0 * pzf) * sinp
    cc = (aa / phsq + sinp) / pzi / pzf
    b = dpzi + pr * sinsq - dpzf * cosp + pxf * sinp
    theta = math.asin(min(ampmax, max(-ampmax, a / phsq)))
    t = np.zeros((6, 6))
    t[0, 0] = cosp - pxf / pzf * sinp
    t[0, 1] = -rhob / pzi / pzf * a
    t[0, 3] = -pyi * rhob / pzi / pzf * b
    t[0, 5] = pr * rhob / pzi / pzf * b
    t[1, 0] = sinp / rhob
    t[1, 1] = pz0 / pzi
    t[1, 3] = pyi / pzi * sinp
    t[1, 5] = -pr / pzi * sinp
    t[2, 0] = -pyi / pzf * sinp
    t[2, 1] = pyi * rhob / pzi / pzf * (dpzf - dpz0)
    t[2, 2] = 1.0
    t[2, 3] = -rhob * (pyi**2 * cc + theta)
    t[2, 5] = pr * pyi * rhob * cc
    t[3, 3] = 1.0
    t[4, 0] = t[1, 0] * t[0, 5] - t[0, 0] * t[1, 5]
    t[4, 1] = t[1, 1] * t[0, 5] - t[0, 1] * t[1, 5]
    t[4, 3] = t[1, 3] * t[0, 5] - t[0, 3] * t[1, 5] + t[2, 5]
    t[4, 4] = 1.0
    t[4, 5] = rhob * (theta - pr**2 * cc)
    t[5, 5] = 1.0
    # cod (lattice._tbedgebody_cod)
    c[1] = min(pr, max(-pr, pxf))
    c[0] = rhob * ddpz + xi * cosp
    c[2] = c[2] - rhob * theta * pyi
    c[4] = c[4] + rhob * theta * pr
    return t


@njit(**NJ)
def _edge_chain(c, rhob, psi, entrance, p2=False):
    """lattice._bend_edge_chain: (cod in-place, M = 各段の trans の積)"""
    if psi == 0.0:
        return _tbedgemaxwell(c, rhob if entrance else -rhob, p2)
    if entrance:
        M = _tbedgedrift(c, psi)
        M = _mm(_tbedgemaxwell(c, rhob, p2), M)
        M = _mm(_tbedgebody(c, rhob, psi), M)
        return M
    M = _tbedgebody(c, rhob, psi)
    M = _mm(_tbedgemaxwell(c, -rhob, p2), M)
    M = _mm(_tbedgedrift(c, psi), M)
    return M


@njit(**NJ)
def _tbendebody0(c, phi0, rhob, drhob, dvzL, exact_dw=False):
    """lattice._tbendebody0_mat (+ (5,6) に dvz·aln) と _tbendebody0_cod (in-place)"""
    xi = c[0]
    dp = c[5]
    p = 1.0 + dp
    pxi = min(p, max(-p, c[1]))
    yi = c[2]
    pyi = min(p, max(-p, c[3]))
    sn = math.sin(phi0)
    cn = math.cos(phi0)
    dcn = cn - 1.0
    dpzi = _pxy2dpz3(pxi, pyi, p)
    pzi = p + dpzi
    dpza = dpzi + (drhob - xi) / rhob
    pza = dpza + dp
    dpxf = dcn * pxi + sn * pza
    pxf = pxi + dpxf
    dpzf = _pxy2dpz3(pxf, pyi, p)
    pzf = p + dpzf
    s1 = pxf / pzf
    phsq = pzf**2 + pxf**2
    da = math.asin((-p * dpxf - pxf * dpzi + pxi * dpzf) / phsq)
    omega = phi0 + da
    dpzadx = -1.0 / rhob
    dpzadpx = -pxi / pzi
    dpzadpy = -pyi / pzi
    ddpzadp = dpzi / pzi
    dpzadp = p / pzi
    u = (s1 + dpzadpx) / phsq
    t = np.eye(6)
    t[1, 0] = sn * dpzadx
    dm22 = dcn + sn * dpzadpx
    t[1, 1] = 1.0 + dm22
    t[1, 3] = sn * dpzadpy
    t[1, 5] = sn * dpzadp
    dodx = -t[1, 0] / pzf
    dodpx = ((dpzf - dpzi) / pzi - dm22) / pzf
    dodpy = -u * pyi - t[1, 3] / pzf
    if exact_dw:
        # v2 (B164 = B139): dω/dp = ∂θ_i/∂p − ∂θ_f/∂p (lattice._tbendebody0_mat と同じ)
        dpxfdp = sn * dpzadp
        dpzfdp = (p - pxf * dpxfdp) / pzf
        dodp = (-pxi * dpzadp - (pzf * dpxfdp - pxf * dpzfdp)) / phsq
        ddodp = dodp + sn / pzf
    else:
        ddodp = u * p - sn * ddpzadp / pzf
        dodp = ddodp - sn / pzf
    t[0, 0] = cn + sn * s1
    t[0, 1] = rhob * (sn - cn * dpzadpx - s1 * t[1, 1])
    t[0, 3] = rhob * (-cn * dpzadpy - s1 * t[1, 3] - pyi / pzf)
    t[0, 5] = rhob * (-cn * dpzadp - s1 * t[1, 5] + p / pzf)
    t[2, 0] = rhob * pyi * dodx
    t[2, 1] = rhob * pyi * dodpx
    t[2, 3] = rhob * (omega + pyi * dodpy)
    t[2, 5] = rhob * pyi * dodp
    t[4, 0] = -t[0, 0] * t[1, 5] + t[1, 0] * t[0, 5]
    t[4, 1] = -t[0, 1] * t[1, 5] + t[1, 1] * t[0, 5]
    t[4, 3] = -t[0, 3] * t[1, 5] + t[1, 3] * t[0, 5] + t[2, 5]
    xsn = phi0 - sn
    t[4, 5] = -rhob * (da + xsn + sn * dpzf / pzf + ddodp * p) + dvzL
    # cod (lattice._tbendebody0_cod: cod[1] は clip なしの値で、da の guard つき)
    pxi0 = c[1]
    pyi0 = c[3]
    dpzi0 = _pxy2dpz3(pxi0, pyi0, p)
    dpza0 = dpzi0 + (drhob - xi) / rhob
    pza0 = dpza0 + dp
    dpxf0 = dcn * pxi0 + sn * pza0
    pxf0 = pxi0 + dpxf0
    dpzf0 = _pxy2dpz3(pxf0, pyi0, p)
    pzf0 = p + dpzf0
    xf = rhob * (sn * pxi0 - cn * dpzi0 - dcn * dp + dpzf0) + cn * xi - dcn * drhob
    phsq0 = pzf0**2 + pxf0**2
    if abs(phsq0) > 1e-30:
        val = (-p * dpxf0 - pxf0 * dpzi0 + pxi0 * dpzf0) / phsq0
        da0 = math.asin(min(1.0, max(-1.0, val)))
    else:
        da0 = 0.0
    omega0 = phi0 + da0
    c[0] = xf
    c[1] = min(p, max(-p, pxf0))
    c[2] = yi + rhob * omega0 * pyi0
    c[3] = min(p, max(-p, pyi0))
    c[4] = c[4] - da0 * p * rhob - phi0 * (drhob + dp * rhob)
    return t


@njit(**NJ)
def _tsteee_body(c, alx, rhob, dvemit, dvzL):
    """lattice._tsteee_body: (trans, cod in-place)。(5,6) に dvz·alx を足す"""
    pr = 1.0 + c[5]
    pxi = c[1] / pr
    pyi = c[3] / pr
    rhoe = rhob * pr
    s = pxi**2 + pyi**2
    dpz1 = -s / (1.0 + math.sqrt(1.0 - s))
    pz1 = 1.0 + dpz1
    dpx = alx / rhoe
    pxf = pxi + dpx
    s = pxf**2 + pyi**2
    dpz2 = -s / (1.0 + math.sqrt(1.0 - s))
    pz2 = 1.0 + dpz2
    phsq = pxi**2 + pz1**2
    d = pxf * pz1 + pxi * pz2
    if d == 0.0 or pxi * pxf < 0.0:
        u = (dpx + pxf * dpz1 - pxi * dpz2) / phsq
        w = (pxf * dpz1 - pxi * dpz2 + dpx * pyi**2) / phsq
    else:
        u = dpx * (pxf + pxi) / d
        w = -dpx * (pxf * dpz1 + pxi * dpz2) / d
    dl = rhoe * (_asinx(u) + w)
    spz = pz2 + pz1
    spx = pxf + pxi
    t = np.eye(6)
    t[0, 1] = (2.0 * spz + spx * (pxf / pz2 + pxi / pz1)) / spz**2 * alx / pr
    t[0, 5] = -spx / pz1 / pz2 / spz * alx / pr
    t[0, 3] = -pyi * t[0, 5]
    t[4, 1] = t[0, 5]
    t[4, 5] = rhob * u / pz1 / pz2
    t[4, 3] = -pyi * t[4, 5]
    t[2, 1] = t[0, 3]
    t[2, 3] = (alx + dl) / pr - pyi * t[4, 3]
    t[2, 5] = t[4, 3]
    t[4, 5] = t[4, 5] - (alx + dl) / pr + dvzL
    c[0] = c[0] + alx * spx / spz
    c[1] = pxf * pr
    c[2] = c[2] + pyi * (alx + dl)
    c[4] = c[4] - (dl + dvemit * alx)
    return t


# ---------------------------------------------------------------------------
# 四極 (tquade: tqfrie, tqlfre, tqente, tsolque 本体) の写像と解析ヤコビアン (emit6d._jac_*)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _tqente(c, al):
    """lattice._tqente_cod (in-place)"""
    pr = 1.0 + c[5]
    pxi = c[1]
    pyi = c[3]
    a = pxi**2 + pyi**2
    if a > 0.0 and abs(pr) > 1e-15:
        arg = max(0.0, 1.0 - a / pr**2)
        dpz = -a / pr / (1.0 + math.sqrt(arg))
    else:
        dpz = 0.0
    pz = pr + dpz
    if abs(pz) < 1e-15:
        pz = pr
    r = -dpz / pz / pr * al
    rpz = pz / pr
    c[0] += r * pxi
    c[2] += r * pyi
    c[4] -= (2.0 + rpz) * a / 2.0 / (pr + pz) * r


@njit(**NJ)
def _jac_tqente(c, al):
    pr = 1.0 + c[5]
    px = c[1]
    py = c[3]
    a = px * px + py * py
    s2 = pr * pr - a
    J = np.eye(6)
    if s2 <= 1e-12:
        return False, J
    s = math.sqrt(s2)
    s3 = s2 * s
    r = al * (1.0 / s - 1.0 / pr)
    r_px = al * px / s3
    r_py = al * py / s3
    r_pr = al * (-pr / s3 + 1.0 / (pr * pr))
    J[0, 1] = r + px * r_px
    J[0, 3] = px * r_py
    J[0, 5] = px * r_pr
    J[2, 1] = py * r_px
    J[2, 3] = r + py * r_py
    J[2, 5] = py * r_pr
    g1 = 2.0 + s / pr
    g2 = a / (2.0 * (pr + s))
    G = g1 * g2
    g1_pr = 1.0 / s - s / (pr * pr)
    g2_pr = -a / (2.0 * (pr + s) ** 2) * (1.0 + pr / s)
    G_pr = g1_pr * g2 + g1 * g2_pr
    J[4, 5] = -(G_pr * r + G * r_pr)
    return True, J


@njit(**NJ)
def _tsolque_body(c, aln1, akkp):
    """lattice._tsolque_body_cod (dx0=dy0=0, in-place)"""
    pr = 1.0 + c[5]
    w1 = math.sqrt(akkp)
    phi1 = min(100.0, max(-100.0, aln1 * w1))
    c1 = math.cos(phi1)
    s1 = math.sin(phi1)
    dc1 = c1 - 1.0
    xs1 = phi1 - s1
    ch2 = math.cosh(phi1)
    sh2 = math.sinh(phi1)
    dch2 = ch2 - 1.0
    xsh2 = phi1 - sh2
    xi = c[0]
    yi = c[2]
    pxi = c[1] / pr
    pyi = c[3] / pr
    u1 = xi * dc1 + pxi * s1 / w1
    u2 = -xi * w1 * s1 + pxi * dc1
    v1 = yi * dch2 + pyi * sh2 / w1
    v2 = yi * w1 * sh2 + pyi * dch2
    pxf = pxi + u2
    pyf = pyi + v2
    c[0] = c[0] + u1
    c[1] = pxf * pr
    c[2] = c[2] + v1
    c[3] = pyf * pr
    c[4] -= 0.25 * (
        w1 * (xi**2 * xs1 - yi**2 * xsh2)
        + (pxi**2 + pyi**2) * aln1
        + u1 * pxf + xi * pxi * dc1
        + v1 * pyf + yi * pyi * dch2)


@njit(**NJ)
def _jac_tsolque_body(c, aln1, kappa, dvz):
    pr = 1.0 + c[5]
    x = c[0]
    y = c[2]
    p = c[1] / pr
    q = c[3] / pr
    w = math.sqrt(kappa / pr)
    phi = aln1 * w
    J = np.eye(6)
    if phi > 100.0:
        return False, J
    cs = math.cos(phi)
    sn = math.sin(phi)
    ch = math.cosh(phi)
    sh = math.sinh(phi)
    dc = cs - 1.0
    dch = ch - 1.0
    xs = phi - sn
    xsh = phi - sh
    dw = -w / (2.0 * pr)
    dphi = aln1 * dw
    p_r = -p / pr
    q_r = -q / pr
    u1 = x * dc + p * sn / w
    pxf = p * cs - x * w * sn
    v1 = y * dch + q * sh / w
    pyf = q * ch + y * w * sh
    J[0, 0] = cs
    J[0, 1] = sn / (w * pr)
    J[0, 5] = -x * sn * dphi + p_r * sn / w + p * (cs * dphi / w - sn * dw / (w * w))
    J[1, 0] = -pr * w * sn
    J[1, 1] = cs
    pxf_r = p_r * cs - p * sn * dphi - x * (dw * sn + w * cs * dphi)
    J[1, 5] = pxf + pr * pxf_r
    J[2, 2] = ch
    J[2, 3] = sh / (w * pr)
    J[2, 5] = y * sh * dphi + q_r * sh / w + q * (ch * dphi / w - sh * dw / (w * w))
    J[3, 2] = pr * w * sh
    J[3, 3] = ch
    pyf_r = q_r * ch + q * sh * dphi + y * (dw * sh + w * ch * dphi)
    J[3, 5] = pyf + pr * pyf_r
    u1_r = -x * sn * dphi + p_r * sn / w + p * (cs * dphi / w - sn * dw / (w * w))
    v1_r = y * sh * dphi + q_r * sh / w + q * (ch * dphi / w - sh * dw / (w * w))
    Z_r = (dw * (x * x * xs - y * y * xsh) + w * dphi * (-x * x * dc + y * y * dch)
           + 2.0 * (p * p_r + q * q_r) * aln1
           + u1_r * pxf + u1 * pxf_r
           + x * (p_r * dc - p * sn * dphi)
           + v1_r * pyf + v1 * pyf_r
           + y * (q_r * dch + q * sh * dphi))
    J[4, 5] = -0.25 * Z_r + dvz * aln1
    return True, J


@njit(**NJ)
def _tqfrie(c, ak, al):
    """lattice._tqfrie_cod (in-place)"""
    xmax = 1e10
    akk = ak / al / 4.0
    pr = 1.0 + c[5]
    xi = min(xmax, max(-xmax, c[0]))
    yi = min(xmax, max(-xmax, c[2]))
    aki = akk / pr
    a = aki * xi**2
    b = aki * yi**2
    d = a + b
    ab = aki * (xi - yi) * (xi + yi)
    t = ab**2 / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    c[0] = xi + dx1
    c[2] = yi + dy1
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h**2
    if f == 0.0:
        return
    pxi = c[1]
    pyi = c[3]
    px1 = (pxi * yy + pyi * h) / f
    py1 = (pyi * xx - pxi * h) / f
    c[1] = px1
    c[3] = py1
    dz1 = -(c[1] * (dx1 + t * xi) + c[3] * (dy1 + t * yi)) / pr
    c[4] = c[4] + dz1


@njit(**NJ)
def _jac_tqfrie(c, ak, al):
    pr = 1.0 + c[5]
    akk = ak / al / 4.0
    k = akk / pr
    x = c[0]
    y = c[2]
    px = c[1]
    py = c[3]
    x2 = x * x
    y2 = y * y
    a = k * x2
    b = k * y2
    ab = a - b
    t = ab * ab / 6.0
    X1 = 1.0 + a / 3.0 + b + t
    Y1 = 1.0 - a - b / 3.0 + t
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    q = 1.0 - ab / 3.0
    h = 2.0 * k * x * y * q
    f = xx * yy + h * h
    J = np.eye(6)
    if f == 0.0:
        return False, J
    a_x = 2.0 * k * x
    a_y = 0.0
    a_k = x2
    b_x = 0.0
    b_y = 2.0 * k * y
    b_k = y2
    ab_x = a_x - b_x
    ab_y = a_y - b_y
    ab_k = a_k - b_k
    t_x = ab * ab_x / 3.0
    t_y = ab * ab_y / 3.0
    t_k = ab * ab_k / 3.0
    X1_x = a_x / 3.0 + b_x + t_x
    X1_y = a_y / 3.0 + b_y + t_y
    X1_k = a_k / 3.0 + b_k + t_k
    Y1_x = -a_x - b_x / 3.0 + t_x
    Y1_y = -a_y - b_y / 3.0 + t_y
    Y1_k = -a_k - b_k / 3.0 + t_k
    xx_x = a_x + b_x + (ab_x * (5.0 * a - b) + ab * (5.0 * a_x - b_x)) / 6.0
    xx_y = a_y + b_y + (ab_y * (5.0 * a - b) + ab * (5.0 * a_y - b_y)) / 6.0
    xx_k = a_k + b_k + (ab_k * (5.0 * a - b) + ab * (5.0 * a_k - b_k)) / 6.0
    yy_x = -(a_x + b_x) + (ab_x * (a - 5.0 * b) + ab * (a_x - 5.0 * b_x)) / 6.0
    yy_y = -(a_y + b_y) + (ab_y * (a - 5.0 * b) + ab * (a_y - 5.0 * b_y)) / 6.0
    yy_k = -(a_k + b_k) + (ab_k * (a - 5.0 * b) + ab * (a_k - 5.0 * b_k)) / 6.0
    h_x = 2.0 * k * y * q - 2.0 * k * x * y * ab_x / 3.0
    h_y = 2.0 * k * x * q - 2.0 * k * x * y * ab_y / 3.0
    h_k = 2.0 * x * y * q - 2.0 * k * x * y * ab_k / 3.0
    f_x = xx_x * yy + xx * yy_x + 2.0 * h * h_x
    f_y = xx_y * yy + xx * yy_y + 2.0 * h * h_y
    f_k = xx_k * yy + xx * yy_k + 2.0 * h * h_k
    pxn = (px * yy + py * h) / f
    pyn = (py * xx - px * h) / f
    xn_x = X1 + x * X1_x
    xn_y = x * X1_y
    xn_k = x * X1_k
    yn_x = y * Y1_x
    yn_y = Y1 + y * Y1_y
    yn_k = y * Y1_k
    pxn_x = (px * yy_x + py * h_x - pxn * f_x) / f
    pxn_y = (px * yy_y + py * h_y - pxn * f_y) / f
    pxn_k = (px * yy_k + py * h_k - pxn * f_k) / f
    pyn_x = (py * xx_x - px * h_x - pyn * f_x) / f
    pyn_y = (py * xx_y - px * h_y - pyn * f_y) / f
    pyn_k = (py * xx_k - px * h_k - pyn * f_k) / f
    dkdp = -k / pr
    J[0, 0] = xn_x
    J[0, 2] = xn_y
    J[0, 5] = xn_k * dkdp
    J[2, 0] = yn_x
    J[2, 2] = yn_y
    J[2, 5] = yn_k * dkdp
    J[1, 0] = pxn_x
    J[1, 2] = pxn_y
    J[1, 1] = yy / f
    J[1, 3] = h / f
    J[1, 5] = pxn_k * dkdp
    J[3, 0] = pyn_x
    J[3, 2] = pyn_y
    J[3, 1] = -h / f
    J[3, 3] = xx / f
    J[3, 5] = pyn_k * dkdp
    A = x * (a / 3.0 + b + 2.0 * t)
    Bq = y * (-a - b / 3.0 + 2.0 * t)
    A_k = x * (a_k / 3.0 + b_k + 2.0 * t_k)
    B_k = y * (-a_k - b_k / 3.0 + 2.0 * t_k)
    Sg = pxn * A + pyn * Bq
    S_pr = (pxn_k * A + pxn * A_k + pyn_k * Bq + pyn * B_k) * dkdp
    J[4, 5] = Sg / (pr * pr) - S_pr / pr
    return True, J


@njit(**NJ)
def _tqlfre(c, al, ak, f1, f2):
    """lattice._tqlfre_cod (in-place)"""
    af1 = -ak / al * f1 * abs(f1) / 24.0
    af2 = ak / al * f2
    p = 1.0 + c[5]
    a = af1 / p
    b = af2 / p
    ea = math.exp(a)
    bp = b / p
    xf = ea * c[0] + bp * c[1]
    yf = c[2] / ea - bp * c[3]
    pxf = c[1] / ea
    pyf = c[3] * ea
    c[4] = c[4] - ((a * ea * c[0] + bp * ea * (1.0 + 0.5 * a) * pxf) * pxf
                   - (a * c[2] / ea + bp / ea * (1.0 - 0.5 * a) * pyf) * pyf) / p
    c[1] = pxf
    c[3] = pyf
    c[0] = xf
    c[2] = yf


@njit(**NJ)
def _jac_tqlfre(c, al, ak, f1, f2):
    af1 = -ak / al * f1 * abs(f1) / 24.0
    af2 = ak / al * f2
    p = 1.0 + c[5]
    a = af1 / p
    b = af2 / p
    ea = math.exp(a)
    bp = b / p
    x = c[0]
    px = c[1]
    y = c[2]
    py = c[3]
    pxf = px / ea
    pyf = py * ea
    t11 = ea
    t12 = bp
    t16 = (-a * ea * x - bp * 2.0 * px) / p
    t33 = 1.0 / ea
    t34 = -bp
    t36 = (a / ea * y + bp * 2.0 * py) / p
    t22 = 1.0 / ea
    t26 = a / p * pxf
    t44 = ea
    t46 = -a / p * pyf
    J = np.eye(6)
    J[0, 0] = t11; J[0, 1] = t12; J[0, 5] = t16
    J[1, 1] = t22; J[1, 5] = t26
    J[2, 2] = t33; J[2, 3] = t34; J[2, 5] = t36
    J[3, 3] = t44; J[3, 5] = t46
    J[4, 0] = -t11 * t26
    J[4, 1] = -t12 * t26 + t22 * t16
    J[4, 2] = -t33 * t46
    J[4, 3] = -t34 * t46 + t44 * t36
    J[4, 5] = (((2.0 + a) * a * ea * x - bp * ea * (2.0 + 0.5 * a) * pxf) / p * pxf
               - ((2.0 - a) * a / ea * y - bp / ea * (2.0 - 0.5 * a) * pyf) / p * pyf
               - (a * ea * x + bp * ea * (2.0 + a) * pxf) * t26
               + (a / ea * y + bp / ea * (2.0 - a) * pyf) * t46) / p
    return J


@njit(**NJ)
def _rollmat(th):
    """emit6d._roll_matrices の R (入口)。出口は転置"""
    ct = math.cos(th)
    st = math.sin(th)
    R = np.eye(6)
    R[0, 0] = ct; R[1, 1] = ct; R[2, 2] = ct; R[3, 3] = ct
    R[0, 2] = -st; R[1, 3] = -st
    R[2, 0] = st; R[3, 1] = st
    return R


@njit(**NJ)
def _rotin_mat():
    Ri = np.eye(6)
    Ri[0, 0] = 0.0; Ri[1, 1] = 0.0; Ri[2, 2] = 0.0; Ri[3, 3] = 0.0
    Ri[0, 2] = 1.0
    Ri[1, 3] = 1.0
    Ri[2, 0] = -1.0
    Ri[3, 1] = -1.0
    return Ri


@njit(**NJ)
def _quad_r_in(c, dx, dy, th, rot):
    """emit6d._tquade r_in (in-place): ずらす → ロール θ → K1<0 なら −π/2 回転"""
    c[0] -= dx
    c[2] -= dy
    if th != 0.0:
        _rot4(c, math.cos(th), math.sin(th))
    if rot:
        x = c[0]
        px = c[1]
        c[0] = c[2]
        c[1] = c[3]
        c[2] = -x
        c[3] = -px


@njit(**NJ)
def _quad_r_out(c, dx, dy, th, rot):
    if rot:
        x = c[0]
        px = c[1]
        c[0] = -c[2]
        c[1] = -c[3]
        c[2] = x
        c[3] = px
    if th != 0.0:
        _rot4(c, math.cos(th), -math.sin(th))
    c[0] += dx
    c[2] += dy


@njit(**NJ)
def _quad_full(c, L, K1, fringe, f1in, f2in, f1out, f2out, mfring, dx, dy, th, use_dv, PH, want_jac):
    """emit6d._jac_quad_full + lattice._quad_offset_map (+ dv 項): (ok, J)。c は in-place で出口の cod に更新"""
    J = np.eye(6)
    if abs(K1) < 1e-10 or abs(L) < 1e-15:
        return False, J
    ak1 = abs(K1) * L
    rot = K1 * L < 0.0
    dp = c[5]
    plain = (dx == 0.0 and dy == 0.0 and th == 0.0)
    cc = c.copy()
    _quad_r_in(cc, dx, dy, th, rot)
    if want_jac:
        if th != 0.0:
            J = _rollmat(th)
        if rot:
            J = _mm(_rotin_mat(), J)
    if fringe and mfring >= 0 and mfring != 2:
        if want_jac:
            ok, Jf = _jac_tqfrie(cc, ak1, L)
            if not ok:
                return False, J
            J = _mm(Jf, J)
        _tqfrie(cc, ak1, L)
    if mfring == 1 or mfring == 3:
        if want_jac:
            J = _mm(_jac_tqlfre(cc, L, ak1, f1in, f2in), J)
        _tqlfre(cc, L, ak1, f1in, f2in)
    kappa = abs(K1)
    sch = int(PH[QBSCH])
    ndiv = min(1000, _qb.ndiv_calc_nb(L * ak1, sch))
    aln1 = L / ndiv
    pr = 1.0 + dp
    for g in range(_qb.QB_M[sch] * ndiv):
        dl, bl, jb, rad, fin = _qb.group_nb(g, aln1, ndiv, sch)
        if want_jac:
            ok, J1 = _jac_tqente(cc, dl)
            if not ok:
                return False, J
        _tqente(cc, dl)
        if want_jac:
            ok, J2 = _jac_tsolque_body(cc, bl, kappa, 0.0)
            if not ok:
                return False, J
            J = _mm(J2, _mm(J1, J))
        _tsolque_body(cc, bl, kappa / pr)
    dl, bl, jb, rad, fin = _qb.group_nb(_qb.QB_M[sch] * ndiv, aln1, ndiv, sch)
    if want_jac:
        ok, J1 = _jac_tqente(cc, dl)
        if not ok:
            return False, J
        J = _mm(J1, J)
    _tqente(cc, dl)
    if mfring == 2 or mfring == 3:
        if want_jac:
            J = _mm(_jac_tqlfre(cc, L, ak1, -f1out, f2out), J)
        _tqlfre(cc, L, ak1, -f1out, f2out)
    if fringe and mfring >= 0 and mfring != 1:
        if want_jac:
            ok, Jf = _jac_tqfrie(cc, -ak1, L)
            if not ok:
                return False, J
            J = _mm(Jf, J)
        _tqfrie(cc, -ak1, L)
    if want_jac:
        if rot:
            J = _mm(_rotin_mat().T.copy(), J)
        if th != 0.0:
            J = _mm(_rollmat(th).T.copy(), J)
        if use_dv:
            J[4, 5] += _dvz(PH, dp) * L
    _quad_r_out(cc, dx, dy, th, rot)
    if use_dv:
        cc[4] -= _dv(PH, cc[5]) * L
    c[:] = cc
    return True, J


@njit(**NJ)
def _wrapj(Jc, Rin, Rout, plain):
    if plain:
        return Jc
    return _mm(Rout, _mm(Jc, Rin))


@njit(**NJ)
def k_quad(cod, T, D, B, codr0, transr, S, PH, L, K1, dx, dy, th, fringe, f1in, f2in, f1out, f2out, mfring,
           krad, use_dv, symp_fix, want_jac, want_db, radcod):
    """emit6d._tquade (L≠0, K1≠0)。戻り 0: 成功、1: Python に戻す"""
    c = cod.copy()
    Tl = T.copy()
    Dl = D.copy()
    Bl = B.copy()
    r0 = codr0.copy()
    trl = transr.copy()
    Sl = S.copy()
    if (not krad) or (c[0] == 0.0 and c[1] == 0.0 and c[2] == 0.0 and c[3] == 0.0
                      and dx == 0.0 and dy == 0.0 and th == 0.0):
        ok, J = _quad_full(c, L, K1, fringe, f1in, f2in, f1out, f2out, mfring, dx, dy, th, use_dv, PH, want_jac)
        if not ok:
            return 1
        _piece(J, Tl, Dl, Bl, want_jac, want_db, symp_fix)
        if krad:
            _setr0(c, Tl, r0, trl)
        cod[:] = c
        T[:, :] = Tl
        D[:, :] = Dl
        B[:, :] = Bl
        codr0[:] = r0
        transr[:, :] = trl
        return 0
    ak1 = abs(K1) * L
    rot = K1 * L < 0.0
    plain = (not rot) and th == 0.0
    Rin = np.eye(6)
    if th != 0.0:
        Rin = _rollmat(th)
    if rot:
        Rin = _mm(_rotin_mat(), Rin)
    Rout = Rin.T.copy()
    _setr0(c, Tl, r0, trl)
    if fringe and mfring >= 0 and mfring != 2:
        if want_jac:
            cr = c.copy()
            _quad_r_in(cr, dx, dy, th, rot)
            ok, Jc = _jac_tqfrie(cr, ak1, L)
            if not ok:
                return 1
            _piece(_wrapj(Jc, Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
        _quad_r_in(c, dx, dy, th, rot)
        _tqfrie(c, ak1, L)
        _quad_r_out(c, dx, dy, th, rot)
    if mfring == 1 or mfring == 3:
        if want_jac:
            cr = c.copy()
            _quad_r_in(cr, dx, dy, th, rot)
            _piece(_wrapj(_jac_tqlfre(cr, L, ak1, f1in, f2in), Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
        _quad_r_in(c, dx, dy, th, rot)
        _tqlfre(c, L, ak1, f1in, f2in)
        _quad_r_out(c, dx, dy, th, rot)
    if f1in != 0.0:
        _tradke(c, Tl, Dl, Bl, r0, trl, Sl, PH, f1in, 0.0, want_jac, want_db, radcod)
    else:
        _setr0(c, Tl, r0, trl)
    # itgetqraddiv
    cr = c.copy()
    _quad_r_in(cr, dx, dy, th, rot)
    xk = ak1 * max(1e-6, abs(cr[0]) + abs(cr[2]))
    b = PH[BRHOZ] * (xk / abs(L))
    nrad = int(abs(L * PH[CRAD] / EPSRAD * (PH[H0] * b)**2))
    nq = max(int(PH[EMIDIV] * PH[EMIDIQ] * nrad),
             int(abs(PH[H0] * b / PH[BRHOZ] * PH[ANRAD]) / EPSRAD / 1e6 * PH[EMIDIV] * PH[EMIDIB]))
    nq = min(1000, max(1, nq))
    aln = L / nq
    akn = ak1 / nq
    sch = int(PH[QBSCH])
    for _q in range(nq):
        ndiv = min(1000, _qb.ndiv_calc_nb(aln * akn, sch))
        aln1 = aln / ndiv
        for g in range(_qb.QB_M[sch] * ndiv):
            dl, bl, jb, rad, fin = _qb.group_nb(g, aln1, ndiv, sch)
            cr = c.copy()
            _quad_r_in(cr, dx, dy, th, rot)
            if want_jac:
                ok, J1 = _jac_tqente(cr, dl)
                if not ok:
                    return 1
                c1 = cr.copy()
                _tqente(c1, dl)
                ok, J2 = _jac_tsolque_body(c1, bl, akn / aln, _dvz(PH, c1[5]) if use_dv else 0.0)
                if not ok:
                    return 1
                _piece(_wrapj(_mm(J2, J1), Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
            pr = 1.0 + cr[5]
            _tqente(cr, dl)
            _tsolque_body(cr, bl, akn / aln / pr)
            if use_dv:
                cr[4] -= _dv(PH, cr[5]) * bl
            _quad_r_out(cr, dx, dy, th, rot)
            c[:] = cr
            if rad:
                _tradke(c, Tl, Dl, Bl, r0, trl, Sl, PH, aln1, 0.0, want_jac, want_db, radcod)
        dl, bl, jb, rad, fin = _qb.group_nb(_qb.QB_M[sch] * ndiv, aln1, ndiv, sch)
        cr = c.copy()
        _quad_r_in(cr, dx, dy, th, rot)
        if want_jac:
            ok, J1 = _jac_tqente(cr, dl)
            if not ok:
                return 1
            _piece(_wrapj(J1, Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
        _tqente(cr, dl)
        _quad_r_out(cr, dx, dy, th, rot)
        c[:] = cr
        _tradke(c, Tl, Dl, Bl, r0, trl, Sl, PH, aln1, 0.0, want_jac, want_db, radcod)
    if mfring == 2 or mfring == 3:
        if want_jac:
            cr = c.copy()
            _quad_r_in(cr, dx, dy, th, rot)
            _piece(_wrapj(_jac_tqlfre(cr, L, ak1, -f1out, f2out), Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
        _quad_r_in(c, dx, dy, th, rot)
        _tqlfre(c, L, ak1, -f1out, f2out)
        _quad_r_out(c, dx, dy, th, rot)
    if fringe and mfring >= 0 and mfring != 1:
        if want_jac:
            cr = c.copy()
            _quad_r_in(cr, dx, dy, th, rot)
            ok, Jc = _jac_tqfrie(cr, -ak1, L)
            if not ok:
                return 1
            _piece(_wrapj(Jc, Rin, Rout, plain), Tl, Dl, Bl, want_jac, want_db, symp_fix)
        _quad_r_in(c, dx, dy, th, rot)
        _tqfrie(c, -ak1, L)
        _quad_r_out(c, dx, dy, th, rot)
    if f1out != 0.0:
        _tradke(c, Tl, Dl, Bl, r0, trl, Sl, PH, f1out, 0.0, want_jac, want_db, radcod)
    cod[:] = c
    T[:, :] = Tl
    D[:, :] = Dl
    B[:, :] = Bl
    codr0[:] = r0
    transr[:, :] = trl
    S[:] = Sl
    return 0


# ---------------------------------------------------------------------------
# BEND (tbende.f: tchge → [tsteee | tquade | 本体 (K1=0 / K1≠0)] → tchge)
# ---------------------------------------------------------------------------

@njit(**NJ)
def _tsteee(c, T, D, B, codr0, transr, S, PH, al, k0, e1, e2, enarad, use_dv, symp_fix, want_jac, want_db, radcod, edgefix=False):
    """emit6d._tsteee (ずれは呼び出し側の tchge。F1 なし)。エッジ・本体の行列は SAD の解析式 (lattice._bend_edge_chain, _tsteee_body)"""
    if al == 0.0:
        if k0 != 0.0:
            if want_jac:
                _piece(np.eye(6), T, D, B, want_jac, want_db, symp_fix)
            c[1] -= k0
        return 0
    if k0 == 0.0:
        return k_drift(c, T, D, B, S, PH, al, use_dv, symp_fix, want_jac, want_db)
    krad = enarad
    if krad:
        _setr0(c, T, codr0, transr)
    rhoe = al / k0
    rhob = -rhoe
    ndiv = 1
    if krad:
        b = -PH[BRHOZ] / rhob
        nrad = 1 + int(abs(al / EPSRAD * PH[CRAD] * (PH[H0] * b)**2))
        ndiv = max(int(nrad * PH[EMIDIV] * PH[EMIDIB]),
                   1 + int(abs(k0 * PH[H0] * PH[ANRAD]) / EPSRAD / 1e6 * PH[EMIDIV] * PH[EMIDIB]))
    aln = al / ndiv
    M = _edge_chain(c, rhoe, e1, True, edgefix)
    _piece(M, T, D, B, want_jac, want_db, symp_fix)
    # v2 (B108): 端の角 (AE1・AE2) のある面では場の中の行路が x·tan(e) だけ短い (_bend_core の edgefix と同じ)
    tin = math.tan(e1) if edgefix else 0.0
    tout = math.tan(e2) if edgefix else 0.0
    for n in range(1, ndiv + 1):
        dvemit = _dv(PH, c[5]) if use_dv else 0.0
        dvzL = _dvz(PH, c[5]) * aln if use_dv else 0.0
        t = _tsteee_body(c, aln, rhob, dvemit, dvzL)
        _piece(t, T, D, B, want_jac, want_db, symp_fix)
        if n != ndiv and krad:
            _tradke(c, T, D, B, codr0, transr, S, PH, aln, 0.0, want_jac, want_db, radcod, tin, 0.0)
            tin = 0.0
    M = _edge_chain(c, rhoe, e2, False, edgefix)
    _piece(M, T, D, B, want_jac, want_db, symp_fix)
    if krad:
        _tradke(c, T, D, B, codr0, transr, S, PH, aln, 0.0, want_jac, want_db, radcod, tin, tout)
    return 0


@njit(**NJ)
def _bend_core(c, T, D, B, codr0, transr, S, PH, al0, phi0, k0, k1, e1, e2, enarad, use_dv, symp_fix,
               want_jac, want_db, radcod, edgefix=False, klin=0):
    """emit6d._tbende_core (F1=0, mfring=0, ずれ・ロールなし)。戻り 0: 成功、1: Python に戻す、2: tmulte (未対応)"""
    phib = phi0 + k0
    ak = k1 * al0
    if phi0 == 0.0:
        if ak == 0.0:
            return _tsteee(c, T, D, B, codr0, transr, S, PH, al0, k0, e1, e2, enarad, use_dv, symp_fix,
                           want_jac, want_db, radcod, edgefix)
        elif phib == 0.0:
            # tquade (ずれなし): emit6d は Quadrupole(L, K1=ak/al0) を作って _tquade(krad=enarad)
            return k_quad(c, T, D, B, codr0, transr, S, PH, al0, ak / al0, 0.0, 0.0, 0.0, True, 0.0, 0.0, 0.0, 0.0, 0,
                          enarad, use_dv, symp_fix, want_jac, want_db, radcod)
        return 2
    al = al0
    if enarad:
        _setr0(c, T, codr0, transr)
    phibl = phib / al
    rhob = 1.0 / phibl
    rho0 = al / phi0
    rbc = 1.0
    phic = phi0 * rbc
    b = PH[BRHOZ] / rhob
    ndiv = 1 if ak == 0.0 else 1 + int(abs(phic / EPSBEND))
    if enarad:
        nrad = 1 + int(abs(al0 * rbc / EPSRAD * PH[CRAD] * (PH[H0] * b)**2))
        ndiv = max(ndiv, int(nrad * PH[EMIDIV] * PH[EMIDIB]),
                   1 + int(abs(phib * PH[H0] * PH[ANRAD]) / EPSRAD / 1e6 * PH[EMIDIV] * PH[EMIDIB]))
    M = _edge_chain(c, rhob, e1, True, edgefix)
    _piece(M, T, D, B, want_jac, want_db, symp_fix)
    # v2 (edgefix): 端の角 e のある面では、場の中の行路が x*tan(e) だけ短い。最初の tradke に入口、最後の tradke に出口の分を渡す
    tin = math.tan(e1) if edgefix else 0.0
    tout = math.tan(e2) if edgefix else 0.0
    alc = al * rbc
    aln = alc / ndiv
    phin = phic / ndiv
    if ak == 0.0:
        for n in range(1, ndiv + 1):
            dvzL = _dvz(PH, c[5]) * aln if use_dv else 0.0
            t = _tbendebody0(c, phin, rhob, rhob - rho0, dvzL, edgefix)
            if use_dv:
                c[4] -= _dv(PH, c[5]) * aln
            _piece(t, T, D, B, want_jac, want_db, symp_fix)
            if enarad and n < ndiv:
                _tradke(c, T, D, B, codr0, transr, S, PH, aln, phin, want_jac, want_db, radcod, tin, 0.0)
                tin = 0.0
        alr = aln
        phi1 = phin
    else:
        alx0 = 0.0
        akx0 = 0.0
        alr = aln
        phi1 = phin
        for n in range(1, ndiv + 1):
            alx = aln
            alr = aln
            akx = ak * alx / al0
            phi1 = phi0 * alx / al0
            t, cn = _bendecorr(c, (akx + akx0) * 0.5, (alx + alx0) * 0.5, rho0, PSQMAX, klin)
            _apply_tdb(t, T, D, B, want_db)
            c[:] = cn
            akx0 = akx
            alx0 = alx
            h1emit = S[H1EMIT]
            dvemit = S[DVEMIT] if use_dv else 0.0
            t, cn = _bendebody(c, alx, phi1, akx, rhob, rho0, PH[H0], h1emit, dvemit, PSQMAX, klin)
            _apply_tdb(t, T, D, B, want_db)
            c[:] = cn
            if enarad and n != ndiv:
                _tradke(c, T, D, B, codr0, transr, S, PH, alr, phi1, want_jac, want_db, radcod, tin, 0.0)
                tin = 0.0
        t, cn = _bendecorr(c, akx0 * 0.5, alx0 * 0.5, rho0, PSQMAX, klin)
        _apply_tdb(t, T, D, B, want_db)
        c[:] = cn
    M = _edge_chain(c, rhob, e2, False, edgefix)
    _piece(M, T, D, B, want_jac, want_db, symp_fix)
    if enarad:
        _tradke(c, T, D, B, codr0, transr, S, PH, alr, phi1, want_jac, want_db, radcod, tin, tout)
    return 0


@njit(**NJ)
def k_bend(cod, T, D, B, codr0, transr, S, PH, al0, phi0, k0, k1, e1, e2,
           has_chg, dx, dy, dz, th, dth, alg, phig,
           enarad, use_dv, symp_fix, want_jac, want_db, radcod, edgefix=False, klin=0):
    """emit6d._tbende: (has_chg なら tchge 入口) → _bend_core → (tchge 出口)。状態の複製で計算し成功時だけ書き戻す。
    edgefix = compat.v2(): 端の角の放射の行路 (JX1) と tchge の入口の回転行列 (B09) を直した動作にする。
    klin = compat.bend_ecorr_flags(): bit0 = K1≠0 の本体の K1 の定義を B' 一様にする (v2 のオプション、B160)、bit1 = B161 の dk の修正 (v2)。
    戻り 0: 成功、それ以外: Python に戻す"""
    c = cod.copy()
    Tl = T.copy()
    Dl = D.copy()
    Bl = B.copy()
    r0 = codr0.copy()
    trl = transr.copy()
    Sl = S.copy()
    if has_chg:
        _tchge(c, Tl, Dl, Bl, dx, dy, dz, th, dth, alg, phig, True, want_jac, want_db, symp_fix, edgefix)
    st = _bend_core(c, Tl, Dl, Bl, r0, trl, Sl, PH, al0, phi0, k0, k1, e1, e2, enarad, use_dv, symp_fix,
                    want_jac, want_db, radcod, edgefix, klin)
    if st != 0:
        return st
    if has_chg:
        if phig != 0.0:
            a0 = al0
            p0 = phi0
        else:
            a0 = 0.0
            p0 = 0.0
        _tchge(c, Tl, Dl, Bl, dx, dy, dz, th, dth, alg - a0, phig - p0, False, want_jac, want_db, symp_fix, edgefix)
    cod[:] = c
    T[:, :] = Tl
    D[:, :] = Dl
    B[:, :] = Bl
    codr0[:] = r0
    transr[:, :] = trl
    S[:] = Sl
    return 0
