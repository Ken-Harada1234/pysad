"""
touschek.py — SAD の INTRA (バンチ内散乱 + Touschek 寿命表) の直訳
====================================================================

SAD での経路 (INTRA フラグ ON で emit / Emittance[] を実行したとき):
  temit (temit.f)
    converge ループ:
      tturne (tturne.f) : 放射込み 6D 一周追跡。calint=T なら各要素の入口で tintrb を呼ぶ
        tintrb (tintrb.f): スライス中心でのビーム行列から
            (a) IBS (Bjorken-Mtingwa 型) の拡散を beam(1:21) に加える
            (b) caltouck=T なら Touschek の損失率表 touckl / touckm / toucke を積算する
      tecalc (temit.f) : 固有値分解 + 平衡ビーム行列 → eemx, eemy, eemz, sige, sigz
      tintraconv (temit.f): IBS 反復の収束判定・エミッタンスの更新 (MINCOUP の下限も)、
            収束時に TouschekTable (SAD リスト) を組み立てる
  Packages/emit.n Touschek クラス: TouschekTable を Spline で補間して TouschekLifetime[nx,ny,nz]

このモジュールは emit6d.py (Emittance[] の直訳) の _State / 要素写像 / tecalc をそのまま使い、
beam(22:42) (前回反復の平衡ビーム行列を伝搬したもの) を _IntraState.B2 として追加する。
他の既存ファイルは変更しない。ffscmd.py への接続案は docs/history/work_touschek.md を参照。

主な公開関数:
  touschek_table(lattice, momentum_eV, em, pbunch, ...) -> dict
      SAD の emit (INTRA) 相当。戻り値 'table' が SAD の TouschekTable と同じ構造:
        [id, [dp(120), 1/τ(120)], [[nx(34), ny(34), nz(34)], mx(34×25), my(34×25)], toucke(nlat×120)]
  touschek_lifetime(table, nx, ny, nz, sige) -> float
      emit.n の TouschekLifetime[nx, ny, nz] (INF は None または math.inf)
  touschek_f(x): emit.n / tintrb.f の TouschekF (C(x)/x の近似)
"""

import math
from typing import List, Optional

import numpy as np

from . import emit6d as _e
from .emit6d import (_State, _Phys, _p2h, _sqrt1, _tinv6, _iaidx, _pack, _unpack,
                     tecalc, CVELOC, ELRADI)
from .elements import Drift, Quadrupole, Bend
from . import lattice as _lat

# ---------------------------------------------------------------------------
# touschek_table モジュール (temit.f 先頭)
# ---------------------------------------------------------------------------
NTOUCKX = 34      # 横方向振幅 (nx, ny) の表の長さ
NTOUCKZ = 25      # 縦方向振幅 (nz) の表の長さ
NTOUCKL = 120     # 運動量受容 dp の表の長さ (dp = 0.004 .. 0.242, 0.002 刻み)

# tintraconv (temit.f) の定数
_ITMAX = 100
_RESIB = 3e-6
_DCMIN = 1e-6


def initialize_tampl() -> np.ndarray:
    """temit.f initialize_tampl: 振幅表 tampl(34, 3) (列 0:x, 1:y, 2:z)。0 始まりの配列で返す。"""
    t = np.zeros((NTOUCKX, 3))
    t[0:4, 0] = [1.0, 2.0, 3.0, 4.0]
    t[4:9, 0] = [5.0, 6.0, 7.0, 8.0, 9.0]
    t[9:27, 0] = [10.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 24.0, 28.0,
                  32.0, 36.0, 41.0, 46.0, 53.0, 60.0, 68.0, 77.0, 88.0]
    # tampl(28:ntouckx,1) = tampl(10:ntouckx-18,1)*10
    t[27:NTOUCKX, 0] = t[9:NTOUCKX - 18, 0] * 10.0
    t[:, 1] = t[:, 0]
    # tampl(7:ntouckx-3,3) = tampl(10:ntouckx,1)
    t[6:NTOUCKX - 3, 2] = t[9:NTOUCKX, 0]
    t[0:6, 2] = [4.6, 5.3, 6.0, 6.8, 7.7, 8.8]
    return t


# ---------------------------------------------------------------------------
# tintrb.f: Touschek 積分の近似 touckf、IBS 積分 fintrbu と bint
# ---------------------------------------------------------------------------
_EEULER = 7.982212789187259
_TA = 5.62966
_TB = 0.75159


def touschek_f(x):
    """tintrb.f touckf / emit.n TouschekF: C(x)/x ∝ 1/τ の近似 (x=0 では 1e200)。
    numpy 配列も受け付ける。"""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = (np.log(1.0 / x / _EEULER + 1.0) * np.exp(-x)
             * (_TB + _EEULER * x) / (_TB + x * (_TA + 2.0 * x))) / x
    v = np.where(x == 0.0, 1e200, v)
    return float(v) if v.ndim == 0 else v


touckf = touschek_f


def fintrbu(u, a1, a2):
    """tintrb.f fintrbu: IBS 積分の被積分関数 (u 変数版)"""
    u21 = 1.0 - u**2
    u214 = u21**4
    u2141 = 1.0 - u214
    return 4.0 * u * u21 * math.sqrt(abs(u214 * u2141 / (
        (u2141 + a1 * u214) * (u2141 + a2 * u214))))


def bint(f, a, b, epslon, epsabs):
    """bint.f: 台形則の逐次細分 + Simpson 外挿による数値積分 (直訳)。
    収束判定が |s - s20| < max(epslon|s|, epsabs) なので、SAD と同じ値を得るには
    同じ反復をたどる必要がある。"""
    itmax = 30
    ba = b - a
    s0 = ba * (f(a) + f(b)) * 0.5
    if math.isnan(s0):
        s0 = 0.0
    xstep = ba
    s20 = s0
    n = 1
    s = 0.0
    for it in range(1, itmax + 1):
        s = 0.0
        x = a - xstep * 0.5
        for _ in range(n):
            x = x + xstep
            s = s + f(x)
            if math.isnan(s):
                s = 0.0
        s = s * xstep
        s1 = (s + s0) * 0.5
        s = (2.0 * s + s0) / 3.0
        if abs(s - s20) < max(epslon * abs(s), epsabs):
            return s
        elif it >= itmax:
            print('BINT Convergence failed. ', s, s20)
            return s
        n = n * 2
        xstep = xstep * 0.5
        s0 = s1
        s20 = s
    return s


# ---------------------------------------------------------------------------
# 3×3 対称行列の固有値分解と連立方程式 (teigens33 / sols33 の代わり)
# ---------------------------------------------------------------------------
# 位置の共分散 xx は静止系で z が γ 倍されるため、対角成分が 1e+2 (z) 〜 1e-12 (y) と
# 14 桁の幅を持つ。LAPACK の dsyevd (numpy.linalg.eigh) は絶対誤差 ~eps·‖A‖ なので
# 最小固有値 (σy²) が 1% 狂い、SAD と要素ごとに ±5e-3 の差が出た (work_touschek.md 参照)。
# SAD の teigens33 (thess3+tqr3) は vx による釣り合い (balancing) 付きの QR 法で相対精度が
# 出る。ここでは正定値対称行列に対して相対精度が保証される Jacobi 法 (巡回, 2 側) を使う。

def _eigs33(a):
    """対称 3×3 行列の固有値 (昇順) と固有ベクトル (列) を Jacobi 法で求める。"""
    A = np.array(a, dtype=float)
    V = np.eye(3)
    for _ in range(60):
        off = 0.0
        for p in range(3):
            for q in range(p + 1, 3):
                dpq = math.sqrt(abs(A[p, p] * A[q, q]))
                if dpq > 0.0:
                    off = max(off, abs(A[p, q]) / dpq)
        if off < 1e-15:
            break
        for p in range(3):
            for q in range(p + 1, 3):
                if A[p, q] == 0.0:
                    continue
                theta = (A[q, q] - A[p, p]) / (2.0 * A[p, q])
                t = (1.0 if theta >= 0.0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                sn = t * c
                tau = sn / (1.0 + c)
                apq = A[p, q]
                A[p, p] -= t * apq
                A[q, q] += t * apq
                A[p, q] = A[q, p] = 0.0
                for r in range(3):
                    if r != p and r != q:
                        arp = A[r, p]
                        arq = A[r, q]
                        A[r, p] = A[p, r] = arp - sn * (arq + tau * arp)
                        A[r, q] = A[q, r] = arq + sn * (arp - tau * arq)
                for r in range(3):
                    vrp = V[r, p]
                    vrq = V[r, q]
                    V[r, p] = vrp - sn * (vrq + tau * vrp)
                    V[r, q] = vrq + sn * (vrp - tau * vrq)
    eig = np.diag(A).copy()
    order = np.argsort(eig)
    return eig[order], V[:, order]


def _sols33(a, b):
    """sols33: a x = b (a: 対称正定値 3×3, b: 3×3)。対角で釣り合わせてから解く (桁幅対策)。"""
    d = np.sqrt(np.abs(np.diag(a)))
    d = np.where(d > 0.0, d, 1.0)
    As = a / d[:, None] / d[None, :]
    bs = b / d[:, None]
    return np.linalg.solve(As, bs) / d[:, None]


# ---------------------------------------------------------------------------
# 追跡状態: beam(22:42) を B2 として追加
# ---------------------------------------------------------------------------

class _IntraState(_State):
    """emit6d._State に beam(22:42) (= B2) を加えたもの。
    SAD では tmulbs(beam, trans, .true.) が calint=T のとき beam(22:42) も同じ行列で変換し、
    tradke は量子励起 dbeam を beam(1:21) と beam(22:42) の両方に加える (photons.f:1518)。"""

    def __init__(self, cod, ph: _Phys, B2):
        super().__init__(cod, ph)
        self.B2 = np.array(B2, dtype=float).copy()

    def apply(self, J):
        super().apply(J)
        self.B2 = J @ self.B2 @ J.T

    def tradke(self, al, phir0):
        # emit6d の tradke は B(1:21) に対して B ← J B J^T + dbeam(B) (J = I + tr1) を行う。
        # SAD は beam(22:42) ← J beam(22:42) J^T + dbeam(beam(1:21)) なので
        # B2_new = B_new + J (B2 - B) J^T。tr1 は D の更新 D_new = D0 + tr1 (T + D0) から復元する。
        D0 = self.D.copy()
        T0 = self.T.copy()
        B0 = self.B.copy()
        super().tradke(al, phir0)
        if not np.array_equal(self.D, D0):
            tr1 = np.linalg.solve((T0 + D0).T, (self.D - D0).T).T
            J = np.eye(6) + tr1
            self.B2 = self.B + J @ (self.B2 - B0) @ J.T


# ---------------------------------------------------------------------------
# temw / touschek_table / tmacro の変数をまとめた文脈
# ---------------------------------------------------------------------------

class _IntraContext:
    def __init__(self, nlat: int, pbunch: float, rclassic: float, coumin: float):
        self.nlat = nlat
        self.pbunch = pbunch
        self.rclassic = rclassic
        self.coumin = coumin        # MINCOUP
        self.intra = True
        self.calint = False
        self.caltouck = False
        # temw: 正規化行列 r (diagr), ri (diagri)、エミッタンス、減衰時間
        self.r = np.eye(6)
        self.ri = np.eye(6)
        self.eemx = self.eemy = self.eemz = 0.0
        self.taurdx = self.taurdy = self.taurdz = 0.0
        self.cintrb = 0.0
        # touschek_table モジュール
        self.tampl = initialize_tampl()
        self.touckl = np.zeros(NTOUCKL)
        self.touckm = np.zeros((NTOUCKZ, NTOUCKX, 2))
        self.toucke = np.zeros((NTOUCKL, nlat))
        self.id = 0

    def reset_tables(self):
        """tturne.f 先頭: calint のとき touckl, touckm, toucke を 0 にする"""
        self.touckl[:] = 0.0
        self.touckm[:] = 0.0
        self.toucke = np.zeros((NTOUCKL, self.nlat))


# ---------------------------------------------------------------------------
# tintrb.f (wspac なし、optics=.false.)
# ---------------------------------------------------------------------------

def _tintrb(st: _IntraState, al: float, al1: float, l: int, ctx: _IntraContext):
    """tintrb の直訳。
    al  : このスライスの長さ (前の要素の中点 〜 この要素の中点 など)
    al1 : 現在位置からスライス始点までの距離 (スライス中心は現在位置から .5*al-al1 だけ先)
    l   : 要素番号 (1 始まり、SAD の l)。toucke の列。
    st.B に IBS の拡散を加え、ctx.caltouck なら Touschek 表を積算する。戻り値は bmi (6×6, 入口座標)。"""
    if al == 0.0:
        return np.zeros((6, 6))
    ph = st.ph
    p0 = ph.p0
    h0 = ph.h0
    cod = st.cod
    pxi = cod[1]
    pyi = cod[3]
    s = pxi**2 + pyi**2
    pr = 1.0 + cod[5]
    pzi = pr * (1.0 + _sqrt1(-s / pr**2))
    trans1 = np.eye(6)
    alx = 0.5 * al - al1
    ale = alx / pzi
    alz = ale / pzi**2
    trans1[0, 1] = ale + pxi**2 * alz
    trans1[0, 3] = pxi * pyi * alz
    trans1[0, 5] = -pxi * alz * pr
    trans1[2, 1] = trans1[0, 3]
    trans1[2, 3] = ale + pyi**2 * alz
    trans1[2, 5] = -pyi * alz * pr
    trans1[4, 1] = trans1[0, 5]
    trans1[4, 3] = trans1[2, 5]
    hi = _p2h(p0 * pr)
    trans1[4, 5] = h0 / hi**3 * alx + s * alz
    # wspac (空間電荷) は未対応
    if not ctx.intra:
        return np.zeros((6, 6))
    bmi = st.B2 + st.B
    if ctx.caltouck:
        transw = st.T @ ctx.r
    a = p0**2 / (hi + 1.0)
    b = a / hi
    d = 1.0 / hi
    trans2 = np.eye(6)
    trans2[0, 0] = 1.0 + a * pxi**2
    trans2[0, 2] = a * pxi * pyi
    trans2[0, 4] = a * pxi * pzi
    trans2[2, 0] = trans2[0, 2]
    trans2[2, 2] = 1.0 + a * pyi**2
    trans2[2, 4] = a * pyi * pzi
    trans2[4, 0] = trans2[0, 4]
    trans2[4, 2] = trans2[2, 4]
    trans2[4, 4] = 1.0 + a * pzi**2
    trans2[1, 1] = d + b * (pzi**2 + pyi**2)
    trans2[1, 3] = -b * pxi * pyi
    trans2[1, 5] = -b * pxi * pzi
    trans2[3, 1] = trans2[1, 3]
    trans2[3, 3] = d + b * (pzi**2 + pxi**2)
    trans2[3, 5] = -b * pyi * pzi
    trans2[5, 1] = trans2[1, 5]
    trans2[5, 3] = trans2[3, 5]
    trans2[5, 5] = d + b * s
    trans2[:, 1] = trans2[:, 1] - pxi / pzi * trans2[:, 5]
    trans2[:, 3] = trans2[:, 3] - pyi / pzi * trans2[:, 5]
    trans2[:, 4] = (pxi * trans2[:, 0] + pyi * trans2[:, 2] + pzi * trans2[:, 4]) / pr
    trans2[:, 5] = pr / pzi * trans2[:, 5]
    trans1 = trans2 @ trans1
    bmi = trans1 @ bmi @ trans1.T
    ix = [0, 2, 4]
    ip = [1, 3, 5]
    # 位置の 3×3 共分散 xx (SAD は下三角のみ使う)
    xx = bmi[np.ix_(ix, ix)]
    xx = np.tril(xx) + np.tril(xx, -1).T
    eig, _ = _eigs33(xx)
    vol1 = math.sqrt(max(1e-80, eig[0] * eig[1] * eig[2]))
    vol = math.sqrt((4.0 * math.pi)**3) * vol1
    bm = math.sqrt(min(eig))
    # xp = xx^{-1} <x p>  (sols33),  pl = <p p> - <p x> xp  (位置を固定したときの運動量分散)
    xp = _sols33(xx, bmi[np.ix_(ix, ip)])
    pl = bmi[np.ix_(ip, ip)] - bmi[np.ix_(ip, ix)] @ xp
    pl = np.tril(pl) + np.tril(pl, -1).T
    eig, r3 = _eigs33(pl)
    ptrans = math.sqrt(eig[0] + eig[1] + eig[2])
    pvol = math.sqrt(max(1e-80, eig[0] * eig[1] * eig[2]))
    # ---- Touschek 損失率の積算 ----
    if vol != 0.0 and ctx.caltouck:
        if ptrans != 0.0:
            t2 = _tinv6(trans1 @ transw)
            extrans = (t2[0, 5]**2 + t2[1, 5]**2) * ptrans**2
            eytrans = (t2[2, 5]**2 + t2[3, 5]**2) * ptrans**2
            eztrans = (t2[4, 5]**2 + t2[5, 5]**2) * ptrans**2
            tf = al / vol / (ptrans * p0)**3
            if getattr(ctx, 'debug', None) is not None:
                ctx.debug.append((l, al, vol, ptrans, extrans, eytrans, eztrans, bm, pvol))
            if eztrans != 0.0:
                ri = ctx.ri
                i = np.arange(1, NTOUCKL + 1)
                aez = (0.002 + i * 0.002)**2 * (ri[4, 5]**2 + ri[5, 5]**2)
                v = touschek_f(aez / eztrans) * tf
                ctx.touckl += v
                ctx.toucke[:, l - 1] += v
                tampl = ctx.tampl
                exy = abs(ctx.eemx) + abs(ctx.eemy)
                aex0 = tampl[:, 0]**2 * exy                      # (34,)
                aey0 = tampl[:, 1]**2 * exy
                aez0 = tampl[:NTOUCKZ, 2]**2 * abs(ctx.eemz)     # (25,)
                aexz = 1.0 / (extrans / aex0[None, :] + eztrans / aez0[:, None])   # (25,34)
                ctx.touckm[:, :, 0] += touschek_f(aexz) * tf
                aeyz = 1.0 / (eytrans / aey0[None, :] + eztrans / aez0[:, None])
                ctx.touckm[:, :, 1] += touschek_f(aeyz) * tf
    else:
        ctx.touckl[:] = 1e20
        ctx.touckm[:] = 1e20
    # ---- IBS (バンチ内散乱) の拡散 ----
    a1 = eig[0] / eig[1]
    a2 = eig[0] / eig[2]
    f1 = 2.0 * bint(lambda u: fintrbu(u, a1, a2), 0.0, 1.0, 1e-3, 3e-17) * eig[0]
    a1 = eig[1] / eig[2]
    a2 = eig[1] / eig[0]
    f2 = 2.0 * bint(lambda u: fintrbu(u, a1, a2), 0.0, 1.0, 1e-3, 3e-17) * eig[1]
    a1 = eig[2] / eig[0]
    a2 = eig[2] / eig[1]
    f3 = 2.0 * bint(lambda u: fintrbu(u, a1, a2), 0.0, 1.0, 1e-3, 3e-17) * eig[2]
    e1 = f2 + f3 - 2.0 * f1
    e2 = f3 + f1 - 2.0 * f2
    e3 = f1 + f2 - 2.0 * f3
    pl = r3 @ np.diag([e1, e2, e3]) @ r3.T
    pbunch = ctx.pbunch
    if pbunch != 0.0:
        bn = abs(vol / pbunch)**(1.0 / 3.0)
    else:
        bn = math.inf
    bmax = max(1e-80, min(bm, bn))
    taur = max(ctx.taurdx, ctx.taurdy, ctx.taurdz)
    if pbunch != 0.0 and taur != 0.0:
        bmin = max(ctx.rclassic / (ptrans * p0)**2,
                   math.sqrt(abs(vol / math.pi / (ptrans * p0 / h0 * CVELOC) / pbunch / taur)))
    else:
        bmin = math.inf
    bmi6 = np.zeros((6, 6))
    if bmax > bmin:
        ci = ctx.cintrb * al * math.log(max(1e-300, 2.0 * bmax / bmin)) / vol1 / pvol / h0**4
        # bmi(ia(2,2)) = ci pl(1,1), bmi(ia(4,2)) = ci pl(2,1), ... (運動量ブロックのみ)
        bmi6[np.ix_(ip, ip)] = ci * pl
        t2 = _tinv6(trans1)
        bmi6 = t2 @ bmi6 @ t2.T
        st.B = st.B + bmi6
        bmi6 = _tinv6(st.T) @ bmi6 @ _tinv6(st.T).T
    return bmi6


# ---------------------------------------------------------------------------
# tturne1 (calint=T): 要素ループに tintrb を挟む
# ---------------------------------------------------------------------------

def _tturne_intra(els, cod0, ph: _Phys, trf0: float, rfsw: bool, radcod: bool,
                  ctx: _IntraContext, B2) -> _IntraState:
    """emit6d.tturne と同じ要素分岐に、tturne.f の calint 処理 (tintrb の呼び出し、
    ドリフトの 2 分割) を加えたもの。"""
    st = _IntraState(cod0, ph, B2)
    st.trf0 = trf0
    st.rfsw = rfsw
    st.radcod = radcod
    st.setr0()
    if ctx.calint:
        ctx.reset_tables()
    alid = 0.0
    nels = len(els)
    for idx, el in enumerate(els):
        l = idx + 1
        al = float(el.L)
        if ctx.calint:
            als = al
            if als != 0.0:
                if isinstance(el, Drift):
                    alib = als * 0.25 + alid
                    alid1 = als * 0.25
                else:
                    alib = als * 0.5 + alid
                    alid1 = als * 0.5
                _tintrb(st, alib, alid, l, ctx)
                alid = alid1
        if isinstance(el, Drift):
            if ctx.calint and al != 0.0:
                al = al * 0.5
                _e._tdrife0(st, al)
                _tintrb(st, al, al * 0.5, l, ctx)
            _e._tdrife0(st, al)
        elif isinstance(el, Bend):
            _e._tbende(st, el, enarad=True)
        elif isinstance(el, Quadrupole):
            _e._tquade(st, el, krad=True)
        elif type(el).__name__ in _lat._THIN_KORD:
            _e._tthine(st, el)          # tthine (放射なし)。emit6d.tturne と同じ
        elif type(el).__name__ == 'Cavity':
            _e._tcave(st, el)
        elif type(el).__name__ == 'Multipole':
            _e._tmulte(st, el, enarad=True)
        elif el.L == 0.0:
            continue
        else:
            raise NotImplementedError(f"touschek: 要素 {type(el).__name__} は未対応")
    if ctx.calint:
        if alid != 0.0:
            _tintrb(st, alid, alid, nels + 1, ctx)
    return st


# ---------------------------------------------------------------------------
# eintrb (temit.f): IBS 反復の加速
# ---------------------------------------------------------------------------

def _eintrb(em0, y, emr):
    eps = 1e-10
    itmax = 30
    em = em0
    a = (y - emr) * em0**2
    y1 = y
    for _ in range(itmax):
        em = em + (y1 - em) / (em**3 + 2.0 * a) * em**3
        y1 = emr + a / em**2
        if abs(y1 - em) < eps * em:
            return min(100.0 * em0, max(emr, y1, 0.01 * em0))
    return min(100.0 * em0, max(0.01 * em0, em))


# ---------------------------------------------------------------------------
# temit の converge ループ + tintraconv
# ---------------------------------------------------------------------------

def touschek_table(lattice, momentum_eV: float, em: dict, pbunch: float,
                   nbunch: float = 1.0, mincoup: float = 1e-2,
                   rfsw: bool = True, radcod: bool = False,
                   charge: float = 1.0, amass: float = _e.ELMASS,
                   verbose: bool = True) -> dict:
    """SAD の `emit` (INTRA ON) 相当: IBS 反復を収束させ、TouschekTable を作る。

    Parameters
    ----------
    lattice     : pysad.lattice.Lattice
    momentum_eV : 運動量 [eV/c]
    em          : emit6d.emittance6d(lattice, momentum_eV, rfsw=rfsw, radcod=radcod) の戻り値。
                  closed_orbit と trf0 だけを使う (SAD の tcod / 同期位相の決定に相当)。
    pbunch      : バンチあたり粒子数 (SAD の PBUNCH)。0 なら IBS も Touschek も 0。
    nbunch      : SAD の NBUNCH。Touschek/IBS の計算には使われない (互換のため受け取るだけ)。
    mincoup     : SAD の MINCOUP (tintraconv の emmin=(eemx+eemy)*coumin; SAD 既定 1e-2)
    rfsw, radcod: emittance6d と同じフラグ

    Returns dict:
      'table'        : SAD の TouschekTable と同じ構造 (Python の入れ子リスト)
      'EMITXC','EMITYC','EMITZC' : 収束した IBS 込みエミッタンス (SAD のグローバル)
      'emittance_x', ... : 最終パスの tecalc の値 (SAD の EMITX 等)、'tecalc' に全体
      'iterations'   : 反復回数 (SAD の it)
      'converged'    : 収束したか
      'touckl','touckm','toucke','tampl' : tf を掛ける前の生の表 (numpy)
    """
    els = lattice.elements
    circ = sum(el.L for el in els)
    ph = _Phys(momentum_eV, circ, amass=amass, charge=charge)
    rclassic = charge**2 * ELRADI
    nlat = len(els) + 1
    ctx = _IntraContext(nlat, float(pbunch), rclassic, float(mincoup))
    cod0 = np.array(em['closed_orbit'], dtype=float)
    trf0 = float(em.get('trf0', 0.0))
    p0 = ph.p0
    h0 = ph.h0
    omega0 = ph.omega0

    # temit の初期化
    it = 0
    emx0 = emy0 = emz0 = 0.0
    emxmin = emymin = emzmin = 1e-200
    emxmax = emymax = emzmax = 1e200
    emxr = emyr = emzr = 0.0
    demin = 1e100
    ctx.calint = False
    ctx.caltouck = False
    ceig0 = np.zeros(6, dtype=complex)
    B2 = np.zeros((6, 6))
    table = None
    converged = False
    tec = None
    eemx = eemy = eemz = 0.0
    while True:
        # ---- tturne + tecalc ----
        st = _tturne_intra(els, cod0, ph, trf0, rfsw, radcod, ctx, B2)
        tec = tecalc(st, circ, rfsw=rfsw)
        r = np.array(tec['normal_r'])
        ri = np.array(tec['normal_ri'])
        emitp = np.array(tec['beam_matrix'])
        emitn = ri @ emitp @ ri.T          # 正規モード座標での平衡ビーム行列 (tecalc の emitn)
        eemx = tec['emittance_x']
        eemy = tec['emittance_y']
        eemz = tec['emittance_z']
        sige = tec['momentum_spread']
        sigz = tec['bunch_length']
        tunes = np.array(tec['tunes'])
        synchm = bool(rfsw) and tunes[2] != 0.0
        damp = np.array(tec['damping_rate'])
        ctx.taurdx = -2.0 * math.pi / omega0 / damp[0] if damp[0] != 0.0 else 0.0
        ctx.taurdy = -2.0 * math.pi / omega0 / damp[1] if damp[1] != 0.0 else 0.0
        ctx.taurdz = -2.0 * math.pi / omega0 / damp[2] if damp[2] != 0.0 else 0.0
        # dc = Σ|ceig - ceig0| (tecalc)。T は反復で変わらないので 2 回目以降は 0
        ceig = np.zeros(6, dtype=complex)
        ceig[0::2] = np.exp(2j * math.pi * tunes)
        ceig[1::2] = np.conj(ceig[0::2])
        dc = float(np.sum(np.abs(ceig - ceig0)))
        ceig0 = ceig
        ctx.r = r.copy()
        ctx.ri = ri.copy()
        # ---- tintraconv ----
        emx1, emy1, emz1 = eemx, eemy, eemz
        emmin = (eemx + eemy) * ctx.coumin
        eemx = max(emmin, eemx)
        eemy = max(emmin, eemy)
        eemz = max(emz0 * 0.1, eemz)
        de = 0.0
        if eemx <= 0.0 or eemy <= 0.0 or eemz <= 0.0:
            print(' Negative emittance, No intrabeam/space charge calculation:'
                  f' eem(x,y,z) = {eemx} {eemy} {eemz}, iter = {it}')
            it = _ITMAX + 1
        else:
            eemx = min(emxmax, max(emxmin, eemx))
            eemy = min(emymax, max(emymin, eemy))
            eemz = min(emzmax, max(emzmin, eemz))
            de = ((1.0 - emx0 / eemx)**2 + (1.0 - emy0 / eemy)**2
                  + (1.0 - emz0 / eemz)**2)
            demin = min(de, demin)
            if it >= 20 and verbose:
                if it == 20:
                    print(' Poor convergence... ')
                    print('     EMITX          EMITY          EMITZ           conv')
                print(f'{eemx:15.7g}{eemy:15.7g}{eemz:15.7g}{de:15.7g}')
        if verbose:
            print(f'touschek: it={it:3d} eemx={eemx:.10e} eemy={eemy:.10e} '
                  f'eemz={eemz:.10e} de={de:.3e} dc={dc:.1e} caltouck={ctx.caltouck}')
        # label 7301
        while True:
            if (it > 1 and dc < _DCMIN and de < _RESIB) or it > _ITMAX:
                if de >= _RESIB or dc >= _DCMIN:
                    print(' Intrabeam/space charge convergence failed.')
                elif ctx.intra and not ctx.caltouck:
                    de = _RESIB * 1.01
                    continue          # go to 7301 → もう 1 周 (caltouck=T) 回す
                converged = de < _RESIB and dc < _DCMIN
                tf = rclassic**2 * pbunch * math.sqrt(math.pi) / h0 * omega0 / 2.0 / math.pi / p0 * h0
                if ctx.caltouck:
                    ctx.id += 1
                    dp = np.array([(i + 1) * 2e-3 for i in range(1, NTOUCKL + 1)])
                    rate = ctx.touckl * tf
                    toucke = ctx.toucke * tf
                    mx = (ctx.touckm[:, :, 0] * tf).T      # (34, 25): nx ごとの nz 列
                    my = (ctx.touckm[:, :, 1] * tf).T
                    table = [float(ctx.id),
                             [dp.tolist(), rate.tolist()],
                             [ctx.tampl.T.tolist(), mx.tolist(), my.tolist()],
                             toucke.T.tolist()]
                finished = True
                break
            else:
                ctx.caltouck = ctx.intra and de < _RESIB * 10.0
                if ctx.calint:
                    if ctx.intra:
                        rx = _eintrb(emx0, eemx, emxr) / eemx
                        ry = _eintrb(emy0, eemy, emyr) / eemy
                        rz = _eintrb(emz0, eemz, emzr) / eemz
                        rr = min(100.0, max(0.01, (rx * ry * rz)**(1.0 / 3.0)))
                        eemx = eemx * rr
                        eemy = eemy * rr
                        eemz = eemz * rr
                    elif emx0 != 0.0:
                        # wspac のみの経路 (intra=F)。ここでは到達しない
                        if it >= 20:
                            emxmax = min(max(eemx, emx0), emxmax)
                            emxmin = max(min(eemx, emx0), emxmin)
                            eemx = math.sqrt(emxmax * emxmin)
                        else:
                            eemx = math.sqrt(eemx * emx0)
                        if it > 30:
                            emymax = min(max(eemy, emy0), emymax)
                            emymin = max(min(eemy, emy0), emymin)
                            eemy = math.sqrt(emymax * emymin)
                        else:
                            eemy = math.sqrt(eemy * emy0)
                        emzmax = min(max(eemz, emz0), emzmax)
                        emzmin = max(min(eemz, emz0), emzmin)
                        eemz = math.sqrt(emzmax * emzmin)
                else:
                    emxr = eemx
                    emyr = eemy
                    emzr = eemz
                emx0 = eemx
                emy0 = eemy
                emz0 = eemz
                ctx.calint = True
                ctx.cintrb = rclassic**2 / 4.0 / math.pi * pbunch
                en = _pack(emitn)
                if emx1 > 0.01 * eemx:
                    rx = math.sqrt(eemx / emx1)
                else:
                    en[_iaidx(0, 0)] = eemx
                    en[_iaidx(1, 1)] = eemx
                    rx = 1.0
                if emy1 > 0.01 * eemy:
                    ry = math.sqrt(eemy / emy1)
                else:
                    en[_iaidx(2, 2)] = eemy
                    en[_iaidx(3, 3)] = eemy
                    ry = 1.0
                if emz1 > 0.01 * eemz:
                    rz = math.sqrt(eemz / emz1)
                else:
                    en[_iaidx(4, 4)] = eemz
                    en[_iaidx(5, 5)] = eemz
                    rz = 1.0
                trans1 = np.diag([rx, rx, ry, ry, rz, rz])
                en6 = trans1 @ _unpack(en) @ trans1.T
                rr6 = r.copy()
                if not synchm:
                    # RF なし: z 面は sigz, sige の対角、r の z ブロックは単位行列にする
                    # (SAD は temw の r 自体を書き換えるので、次の tturne の tintrb にもこの r が使われる)
                    en6[4, 0:4] = 0.0
                    en6[0:4, 4] = 0.0
                    en6[4, 4] = sigz**2
                    en6[4, 5] = en6[5, 4] = 0.0
                    en6[5, 5] = sige**2
                    rr6[0:4, 4] = 0.0
                    rr6[4, 4] = 1.0
                    rr6[5, 4] = 0.0
                    rr6[5, 0:5] = 0.0
                    rr6[5, 5] = 1.0
                    rr6[4, 0:4] = 0.0
                    rr6[4, 5] = 0.0
                    ctx.r = rr6
                B2 = rr6 @ en6 @ rr6.T
                it += 1
                finished = False
                break
        ctx.eemx, ctx.eemy, ctx.eemz = eemx, eemy, eemz
        if finished:
            break
    out = dict(tec)
    out['table'] = table
    out['EMITXC'] = eemx
    out['EMITYC'] = eemy
    out['EMITZC'] = eemz
    out['iterations'] = it
    out['converged'] = converged
    out['tecalc'] = tec
    out['touckl'] = ctx.touckl.copy()
    out['touckm'] = ctx.touckm.copy()
    out['toucke'] = ctx.toucke.copy()
    out['tampl'] = ctx.tampl.copy()
    out['pbunch'] = pbunch
    out['B2'] = B2            # 最終パスで使った beam(22:42) の初期値 (診断用)
    out['_ctx'] = ctx         # 最終パスの文脈 (診断用)
    return out


# ---------------------------------------------------------------------------
# Spline (spline.f mode=0, Packages/D.n SplineData) と emit.n の TouschekLifetime
# ---------------------------------------------------------------------------

def _spline_ddy(x, y):
    """spline.f spline(np,x,y,ddy,work,mode=0,dy): 自然スプライン (両端 y''=0)。
    ddy は y''/6。y は (n,) または (n, m) (成分ごとのスプライン = D.n の {_Real,_List} 版)。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    ddy = np.zeros_like(y)
    if n <= 2:
        return ddy
    work = np.zeros(n)
    ddy[1] = (y[2] - y[1]) / (x[2] - x[1]) - (y[1] - y[0]) / (x[1] - x[0])
    work[1] = (x[2] - x[0]) * 2.0
    for i in range(1, n - 2):
        f = x[i + 1] - x[i]
        work[i + 1] = (x[i + 2] - x[i]) * 2.0 - f**2 / work[i]
        ddy[i + 1] = ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1])
                      - (y[i + 1] - y[i]) / f - f / work[i] * ddy[i])
    ddy[n - 1] = 0.0
    ddy[n - 2] = ddy[n - 2] / work[n - 2]
    for i in range(n - 3, 0, -1):
        ddy[i] = (ddy[i] - ddy[i + 1] * (x[i + 1] - x[i])) / work[i]
    return ddy


def _find_index(xl, v):
    """spline.f tffindindex (FindIndex[xl, x]) の実数版: xl[i] <= v < xl[i+1] の i (0 始まり)。
    範囲外は 0 または n-2 (端の区間で外挿される)。"""
    i1 = 0
    i2 = len(xl) - 1
    while i2 > i1 + 1:
        ih = i1 + (i2 - i1) // 2
        if xl[ih] == v:
            return ih
        elif xl[ih] > v:
            i2 = ih
        else:
            i1 = ih
    return i1


class SadSpline:
    """Spline[{{x,y},...}] (D.n) の直訳。y はスカラー列でもベクトル列でもよい。"""

    def __init__(self, x, y):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        if len(self.x) == 1:
            self.ddy = np.zeros_like(self.y)
        else:
            self.ddy = _spline_ddy(self.x, self.y)

    def __call__(self, v):
        """SplineData[{xl,y,ddy}][x]"""
        if len(self.x) == 1:
            return self.y[0]
        i = _find_index(self.x, v)
        x0 = self.x[i]
        dx = self.x[i + 1] - x0
        b = (v - x0) / dx
        a = 1.0 - b
        g = -a * b * dx**2
        c = (a + 1.0) * g
        d = (b + 1.0) * g
        return a * self.y[i] + b * self.y[i + 1] + c * self.ddy[i] + d * self.ddy[i + 1]


def _is_inf(v) -> bool:
    return v is None or (isinstance(v, float) and math.isinf(v))


def _restrict(v, lo, hi):
    return min(max(v, lo), hi)


class TouschekLifetimeCalculator:
    """emit.n の Touschek クラス。Setup[TouschekTable] でスプラインを準備する。"""

    def __init__(self, table=None):
        self.id = None
        if table is not None:
            self.setup(table)

    def setup(self, d):
        if d is None:
            raise ValueError('TouschekTable Undefined; Do Emittance[] with INTRA;')
        if self.id is not None and d[0] == self.id:
            return
        self.id = d[0]
        self.spz = SadSpline(d[1][0], d[1][1])
        self.spxz = SadSpline(d[2][0][0], d[2][1])
        self.spyz = SadSpline(d[2][0][1], d[2][2])
        nzl = len(self.spxz.y[0])
        self.zl = np.asarray(d[2][0][2][:nzl], dtype=float)
        self.dpmin, self.dpmax = min(d[1][0]), max(d[1][0])
        self.nxmin, self.nxmax = min(d[2][0][0]), max(d[2][0][0])
        self.nymin, self.nymax = min(d[2][0][1]), max(d[2][0][1])
        self.nzmin, self.nzmax = float(self.zl.min()), float(self.zl.max())

    def lifetime(self, nx, ny, nz, sige):
        """TouschekLifetime[nx, ny, nz] (INF は None か math.inf)。SIGE は引数で渡す。"""
        if _is_inf(nx) and _is_inf(ny):
            sz = nz * sige
            if sz < self.dpmin:
                return (sz / self.dpmin)**2 / self.spz(self.dpmin)
            return 1.0 / self.spz(sz)
        if _is_inf(ny):
            if nz < self.nzmin:
                return (nz / self.nzmin)**2 * self.lifetime(nx, ny, self.nzmin, sige)
            if nx < self.nxmin:
                return (nx / self.nxmin)**2 * self.lifetime(self.nxmin, ny, nz, sige)
            v = self.spxz(min(nx, self.nxmax))
            return 1.0 / SadSpline(self.zl, v)(min(nz, self.nzmax))
        if _is_inf(nx):
            v = self.spyz(_restrict(ny, self.nymin, self.nymax))
            return 1.0 / SadSpline(self.zl, v)(_restrict(nz, self.nzmin, self.nzmax))
        return math.nan


def touschek_lifetime(table, nx, ny, nz, sige) -> float:
    """emit.n TouschekLifetime[nx, ny, nz] の直訳。
    table: touschek_table()['table'] (SAD の TouschekTable と同じ構造)
    nx, ny: None または math.inf で INF。両方 INF なら nz*sige の運動量受容に対する寿命 [s]、
    片方だけ INF なら横振幅 nx (または ny) と縦振幅 nz に対する寿命。
    sige: SAD の SIGE (相対運動量広がり)。"""
    return TouschekLifetimeCalculator(table).lifetime(nx, ny, nz, sige)
