"""
match6d.py — SAD の解析微分 (FFS$NumericalDerivative=False の既定経路) の直訳

SAD 側の対応:
  src/qdcell.f  tffsqu      … 変数ごと・条件ごとに qu(条件, 変数) = Σ_k coup·dtwiss(kf) を組み立てる (本ファイル tffsqu)
  src/qdcell.f  qdcell      … 要素 k の変数 iv による位置 l の Twiss 微分 dtwiss。CELL では周期解の補正を加える
  src/qdcell.f  dgammag     … GMX/GMY (と BMAG) の微分
  src/qdcell.f  tffscoupmatrix … CouplingMatrix フック (未移植: 常に無し)
  src/qdtwis.f  qdtwis      … 輸送 (INS 的) の微分: 要素の微分行列 dtrans を l まで運び、β/α/ψ/η/R の微分を作る
  src/qdtwis.f  qddtwiss    … 閉軌道の数値伝搬 (ctrans)。CELL では DZ/DDP 条件以外で呼ばれないので未移植
  src/qdtwis.f  qdtrans     … 要素写像の数値微分 (MULT/SOL 等)。未移植 (数値微分列で代用)
  src/qdquad.f  qdquad      … QUAD K1 (iv=2) の微分行列、qdrotate … QUAD ROTATE (iv=4)
  src/qdbend.f  qdbend      … BEND ANGLE (iv=2) / K1 (iv=8) の微分行列
  src/qdmdia.f  qdmdia      … 4×4 一周行列の対角化行列 R の微分
  src/qgettr.f  qgettru     … Twiss から作る k1→k2 の 4×5 転送行列 (正規モード → 物理座標へ R で戻す)
  src/tftmat.f  tftmatu     … qgettru に入口側の R 変換を掛けた物理座標の転送行列
  src/qtent.f   qtentu      … 位置の R 変換行列 (4×4) と分散列・閉軌道
  src/qchg.f    qchg        … 要素のずれ・回転 (DX, DY, ROT) の座標変換
  src/toplvl.f  tmultr45    … 4×5 行列の合成 c = b∘a

qu の定義 (tffsmatch.f 402 行目): qu[条件, 変数] = (tffsqu の値)·wiq[条件]/wvar[変数]。
match6.py の数値微分 qu = (ddf1−ddf2)/2/dvkc/wvar = −∂(df·wiq)/∂v/wvar = +∂f/∂v·wiq/wvar と同じ符号・正規化
(df = 目標 − f なので)。BX/BY は log の微分 dβ/β、AX/AY は atan の微分 dα/(1+α²) を dtwiss が直接返す。

使い方:
  from pysad.match6d import analytic_qu, AnalyticMatcher
  qu, numeric_cols = analytic_qu(matcher, ids, rows, wvar)   # match6.Matcher の go() 内の数値微分の代わり
  AnalyticMatcher は Matcher の go() を解析微分 (FitFunction 行は SAD と同じ片側差分) に置き換えたサブクラス。
検証 (docs/history/work_analytic_deriv.md): match6 の数値微分と Level 13B で最大相対差 2e-8、ESRF ARCA (BEND K1 含む) で 4e-7。
AnalyticMatcher で Level 13B と esrf/matching (matching1) の SAD の反復列・残差・K1 (12 桁) を再現。
注意: SAD の mathfun sinc(x) は x·cos x − sin x (sin x / x ではない)。qdbend はこれを使う。
制限: CELL のみ検証済み。変数は QUAD K1/ROTATE, BEND ANGLE/K1 (他は数値微分列にフォールバック)。
     条件は AX..EPY, R1..R4, DETR, DX..DPY, PEX..PEPY, GMX/GMY, TRX/TRY, LENG。
     L=0 の QUAD (qdthin)、MARK 先頭の変数 (qdini)、DX≠0 で geocal の QUAD は未移植 (微分 0)。
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from .match6 import Matcher, _sad_get, EPS, EPS1, _fuzz, _log10, _fmincube, _Res, PI2, M_EULER, \
    FLIM1, FLIM2, AIMP1, AIMP2, BADC1, BADC2, AMEDC1, AMEDC2, ALIT, WLMIN, RTOL1, AMTOL
from .sadscript import fmt_g

# ---------------------------------------------------------------------------
# tffitcode の添字 (tffs.f 104 行目〜)
# ---------------------------------------------------------------------------
AX, BX, NX, AY, BY, NY = 1, 2, 3, 4, 5, 6
EX, EPX, EY, EPY = 7, 8, 9, 10
R1, R2, R3, R4, DETR = 11, 12, 13, 14, 15
DX, DPX, DY, DPY, DZ, DDP = 16, 17, 18, 19, 20, 21
AZ, BZ, NZ = 22, 23, 24
ZX, ZPX, ZY, ZPY = 25, 26, 27, 28
PEX, PEPX, PEY, PEPY = 29, 30, 31, 32
PZX, PZPX, PZY, PZPY = 33, 34, 35, 36
GMX, GMY, GMZ = 37, 38, 39
BMAGX, BMAGY, BMAGZ = 40, 41, 42
TRX, TRY, TRZ = 43, 44, 45
LENG = 46
NTWISSFUN = ZPY            # ntwissfun
MFITO = BMAGZ              # mfito
MFITTRY = TRY
NDT = LENG + 1             # dtwiss の長さ (添字 1..46 を使う)

_KF_INDEX = {'AX': AX, 'BX': BX, 'NX': NX, 'AY': AY, 'BY': BY, 'NY': NY, 'EX': EX, 'EPX': EPX, 'EY': EY, 'EPY': EPY,
             'R1': R1, 'R2': R2, 'R3': R3, 'R4': R4, 'DETR': DETR, 'DX': DX, 'DPX': DPX, 'DY': DY, 'DPY': DPY,
             'DZ': DZ, 'DDP': DDP, 'AZ': AZ, 'BZ': BZ, 'NZ': NZ, 'ZX': ZX, 'ZPX': ZPX, 'ZY': ZY, 'ZPY': ZPY,
             'PEX': PEX, 'PEPX': PEPX, 'PEY': PEY, 'PEPY': PEPY, 'PZX': PZX, 'PZPX': PZPX, 'PZY': PZY, 'PZPY': PZPY,
             'GMX': GMX, 'GMY': GMY, 'GMZ': GMZ, 'BMAGX': BMAGX, 'BMAGY': BMAGY, 'BMAGZ': BMAGZ,
             'TRX': TRX, 'TRY': TRY, 'TRZ': TRZ, 'LENG': LENG}

# 要素型ごとのキーワード → SAD の iv (initb1.f: ky_K1_QUAD=2, ky_ROT_QUAD=4, ky_ANGL_BEND=2, ky_K1_BEND=8)
KY_K1_QUAD, KY_ROT_QUAD = 2, 4
KY_ANGL_BEND, KY_K1_BEND = 2, 8
_IV_TABLE = {('Quadrupole', 'K1'): KY_K1_QUAD, ('Quadrupole', 'ROTATE'): KY_ROT_QUAD,
             ('Quadrupole', 'ROTATION'): KY_ROT_QUAD,
             ('Bend', 'ANGLE'): KY_ANGL_BEND, ('Bend', 'K1'): KY_K1_BEND}


# ---------------------------------------------------------------------------
# 小道具
# ---------------------------------------------------------------------------
def _sinc(x: float) -> float:
    """SAD mathfun (tfloor.f) の sinc: x·cos(x) − sin(x) (sin(x)/x ではない)。|x| ≤ 0.1 は級数"""
    if abs(x) > 0.1:
        return x * math.cos(x) - math.sin(x)
    x2 = -x**2
    return x * x2 * (1.0 / 3.0 + x2 * (1.0 / 30.0 + x2 * (1.0 / 840.0 + x2 * (1.0 / 45360.0 + x2 / 3991680.0))))


def _sinhc(x: float) -> float:
    """SAD mathfun の sinhc: x·cosh(x) − sinh(x)"""
    if abs(x) > 0.1:
        return x * math.cosh(x) - math.sinh(x)
    x2 = x**2
    return x * x2 * (1.0 / 3.0 + x2 * (1.0 / 30.0 + x2 * (1.0 / 840.0 + x2 * (1.0 / 45360.0 + x2 / 3991680.0))))


def _sign(a: float, b: float) -> float:
    """Fortran の sign(a, b): |a| に b の符号"""
    return abs(a) if b >= 0.0 else -abs(a)


def tmultr45(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """toplvl.f tmultr45: c = b∘a (先に a、次に b)。c(:,1:4) = b(:,1:4)·a, c(:,5) = b(:,1:4)·a(:,5) + b(:,5)"""
    c = b[:, 0:4] @ a
    c[:, 4] += b[:, 4]
    return c


def _resetnan(a: np.ndarray) -> np.ndarray:
    a[~np.isfinite(a)] = 0.0
    return a


def build_utwiss(rows) -> np.ndarray:
    """twiss_table の行 (rows[l-1] = SAD の位置 l、1 始まり) → utwiss(l, 1..ntwissfun) 配列 (行 0 は未使用)"""
    n = len(rows)
    u = np.zeros((n + 1, NTWISSFUN + 1))
    for l, r in enumerate(rows, start=1):
        u[l, AX], u[l, BX], u[l, NX] = r.alfx, r.betx, r.mux
        u[l, AY], u[l, BY], u[l, NY] = r.alfy, r.bety, r.muy
        u[l, EX], u[l, EPX], u[l, EY], u[l, EPY] = r.etx, r.etpx, r.ety, r.etpy
        u[l, R1], u[l, R2], u[l, R3], u[l, R4], u[l, DETR] = r.r1, r.r2, r.r3, r.r4, r.detr
        u[l, DX], u[l, DPX], u[l, DY], u[l, DPY] = r.x_cod, r.px_cod, r.y_cod, r.py_cod
        u[l, DZ], u[l, DDP] = r.z_cod, r.delta_cod
        u[l, AZ], u[l, BZ], u[l, NZ] = r.alfz, r.betz, r.muz
        u[l, ZX], u[l, ZPX], u[l, ZY], u[l, ZPY] = r.zx, r.zpx, r.zy, r.zpy
    return u


# ---------------------------------------------------------------------------
# qtent.f qtentu / qchg.f qchg
# ---------------------------------------------------------------------------
def qtentu(u1: np.ndarray, disp: bool) -> Tuple[np.ndarray, np.ndarray]:
    """位置の R 変換 (物理座標 → 正規モード) 4×4 と、disp なら分散列 (5 列目)、閉軌道 cod(6)"""
    r1, r2, r3, r4 = u1[R1], u1[R2], u1[R3], u1[R4]
    detr = r1 * r4 - r2 * r3
    cc = math.sqrt(1.0 - detr)
    normal = u1[DETR] < 1.0
    t = np.zeros((4, 5))
    if normal:
        t[0, 0], t[0, 2], t[0, 3] = cc, r4, -r2
        t[1, 1], t[1, 2], t[1, 3] = cc, -r3, r1
        t[2, 0], t[2, 1], t[2, 2] = -r1, -r2, cc
        t[3, 0], t[3, 1], t[3, 3] = -r3, -r4, cc
    else:
        t[0, 0], t[0, 1], t[0, 2] = r4, -r2, cc
        t[1, 0], t[1, 1], t[1, 3] = -r3, r1, cc
        t[2, 0], t[2, 2], t[2, 3] = cc, -r1, -r2
        t[3, 1], t[3, 2], t[3, 3] = cc, -r3, -r4
    if disp:
        ex, epx, ey, epy = u1[EX], u1[EPX], u1[EY], u1[EPY]
        if normal:
            t[0, 4] = cc * ex + r4 * ey - r2 * epy
            t[1, 4] = cc * epx - r3 * ey + r1 * epy
            t[2, 4] = cc * ey - r1 * ex - r2 * epx
            t[3, 4] = cc * epy - r3 * ex - r4 * epx
        else:
            t[0, 4] = -r1 * ex - r2 * epx + cc * ey
            t[1, 4] = -r3 * ex - r4 * epx + cc * epy
            t[2, 4] = cc * ex + r4 * ey - r2 * epy
            t[3, 4] = cc * epx - r3 * ey + r1 * epy
    cod = np.array([u1[DX], u1[DPX], u1[DY], u1[DPY], u1[DZ], u1[DDP]])
    return t, cod


def qchg(trans: np.ndarray, cod: np.ndarray, dx: float, dy: float, theta: float, enter: bool) -> None:
    """qchg.f: 要素のずれ (dx, dy) と回転 theta。in place"""
    if enter:
        cod[0] += dx
        cod[2] += dy
    if theta != 0.0:
        cost, sint = math.cos(theta), math.sin(theta)
        x = trans[0, :].copy()
        trans[0, :] = cost * x - sint * trans[2, :]
        trans[2, :] = sint * x + cost * trans[2, :]
        y = trans[1, :].copy()
        trans[1, :] = cost * y - sint * trans[3, :]
        trans[3, :] = sint * y + cost * trans[3, :]
        xi = cod[0]
        cod[0] = cost * xi - sint * cod[2]
        cod[2] = sint * xi + cost * cod[2]
        pxi = cod[1]
        cod[1] = cost * pxi - sint * cod[3]
        cod[3] = sint * pxi + cost * cod[3]
    if not enter:
        cod[0] += dx
        cod[2] += dy


# ---------------------------------------------------------------------------
# qgettr.f qgettru / tftmat.f tftmatu
# ---------------------------------------------------------------------------
def qgettru(u1: np.ndarray, u2: np.ndarray, amux: float, amuy: float, k1: int, k2: int,
            disp: bool, fold: bool, trpt: bool) -> np.ndarray:
    """Twiss (u1 at k1, u2 at k2) から k1→k2 の 4×5 転送行列。k2<k1 (または k2==k1 で fold) は一周をまたぐ。
    disp=True なら 5 列目は k2 の分散そのもの、False なら分散の非斉次項 η2 − M·η1 (gammab は一定 → gr=1)。
    出口側 (u2) の R で正規モード → 物理座標に戻す。"""
    trans = np.zeros((4, 5))
    sqrtb1, sqrtb2 = math.sqrt(u1[BX]), math.sqrt(u2[BX])
    dpsi1 = u2[NX] - u1[NX]
    if not trpt and (k2 < k1 or (k2 == k1 and fold)):
        dpsi1 = amux + dpsi1
    th = math.tan(0.5 * dpsi1)
    sinpsi = 2.0 * th / (1.0 + th**2)
    cospsi = 1.0 - th * sinpsi
    trans[0, 0] = sqrtb2 / sqrtb1 * (cospsi + u1[AX] * sinpsi)
    trans[0, 1] = sqrtb2 * sqrtb1 * sinpsi
    trans[1, 0] = -((u2[AX] - u1[AX]) * cospsi + (u1[AX] * u2[AX] + 1.0) * sinpsi) / sqrtb1 / sqrtb2
    trans[1, 1] = sqrtb1 / sqrtb2 * (cospsi - u2[AX] * sinpsi)
    gr = 1.0
    if disp:
        trans[0, 4] = u2[EX]
        trans[1, 4] = u2[EPX]
    else:
        trans[0, 4] = u2[EX] - gr * (trans[0, 0] * u1[EX] + trans[0, 1] * u1[EPX])
        trans[1, 4] = u2[EPX] - gr * (trans[1, 0] * u1[EX] + trans[1, 1] * u1[EPX])
    sqrtb1, sqrtb2 = math.sqrt(u1[BY]), math.sqrt(u2[BY])
    dpsi1 = u2[NY] - u1[NY]
    if not trpt and (k2 < k1 or (k2 == k1 and fold)):
        dpsi1 = amuy + dpsi1
    th = math.tan(0.5 * dpsi1)
    sinpsi = 2.0 * th / (1.0 + th**2)
    cospsi = 1.0 - th * sinpsi
    trans[2, 2] = sqrtb2 / sqrtb1 * (cospsi + u1[AY] * sinpsi)
    trans[2, 3] = sqrtb2 * sqrtb1 * sinpsi
    trans[3, 2] = -((u2[AY] - u1[AY]) * cospsi + (u1[AY] * u2[AY] + 1.0) * sinpsi) / sqrtb1 / sqrtb2
    trans[3, 3] = sqrtb1 / sqrtb2 * (cospsi - u2[AY] * sinpsi)
    if disp:
        trans[2, 4] = u2[EY]
        trans[3, 4] = u2[EPY]
    else:
        trans[2, 4] = u2[EY] - gr * (trans[2, 2] * u1[EY] + trans[2, 3] * u1[EPY])
        trans[3, 4] = u2[EPY] - gr * (trans[3, 2] * u1[EY] + trans[3, 3] * u1[EPY])
    r1, r2, r3, r4 = u2[R1], u2[R2], u2[R3], u2[R4]
    detr = r1 * r4 - r2 * r3
    cc = math.sqrt(1.0 - detr)
    t1 = np.zeros((4, 5))
    if u2[DETR] < 1.0:
        t1[0, :] = cc * trans[0, :] + r4 * trans[2, :] - r2 * trans[3, :]
        t1[1, :] = cc * trans[1, :] - r3 * trans[2, :] + r1 * trans[3, :]
        t1[2, :] = cc * trans[2, :] - r1 * trans[0, :] - r2 * trans[1, :]
        t1[3, :] = cc * trans[3, :] - r3 * trans[0, :] - r4 * trans[1, :]
    else:
        t1[0, :] = cc * trans[2, :] - r1 * trans[0, :] - r2 * trans[1, :]
        t1[1, :] = cc * trans[3, :] - r3 * trans[0, :] - r4 * trans[1, :]
        t1[2, :] = cc * trans[0, :] + r4 * trans[2, :] - r2 * trans[3, :]
        t1[3, :] = cc * trans[1, :] - r3 * trans[2, :] + r1 * trans[3, :]
    return t1


def tftmatu(u1: np.ndarray, u2: np.ndarray, amux: float, amuy: float, l1: int, l2: int,
            fold: bool, trpt: bool) -> np.ndarray:
    """tftmat.f tftmatu: qgettru (disp=.false.) に入口側 (u1) の R 変換を右から掛け、物理座標 l1 → 物理座標 l2 にする"""
    trans = qgettru(u1, u2, amux, amuy, l1, l2, False, fold, trpt)
    r1, r2, r3, r4 = u1[R1], u1[R2], u1[R3], u1[R4]
    detr = r1 * r4 - r2 * r3
    cc = math.sqrt(1.0 - detr)
    x, px, y, py = trans[:, 0].copy(), trans[:, 1].copy(), trans[:, 2].copy(), trans[:, 3].copy()
    if u1[DETR] < 1.0:
        trans[:, 0] = cc * x + r1 * y + r3 * py
        trans[:, 1] = cc * px + r2 * y + r4 * py
        trans[:, 2] = cc * y - r4 * x + r3 * px
        trans[:, 3] = cc * py + r2 * x - r1 * px
    else:
        trans[:, 0] = cc * y + (-r4 * x + r3 * px)
        trans[:, 1] = cc * py + (r2 * x - r1 * px)
        trans[:, 2] = cc * x + (r1 * y + r3 * py)
        trans[:, 3] = cc * px + (r2 * y + r4 * py)
    return trans


# ---------------------------------------------------------------------------
# qdquad.f qdquad / qdrotate
# ---------------------------------------------------------------------------
def qdquad(al: float, ak: float, u_k: np.ndarray, u_k1: np.ndarray, amux: float, amuy: float, k1: int,
           dx: float, dy: float, theta: float, iv: int, trpt: bool) -> Tuple[np.ndarray, np.ndarray]:
    """QUAD の微分行列 dtrans(4,5) と閉軌道の微分 dcod(6)。iv=2: K1 (ak は SAD の積分値)、iv=4: ROTATE。
    dtrans = (∂M/∂v)·S_k (S_k = 入口の R 変換、5 列目は ∂M/∂v·η_k + 閉軌道由来の項)"""
    dtrans = np.zeros((4, 5))
    dcod = np.zeros(6)
    if iv == KY_ROT_QUAD:
        return qdrotate(u_k, u_k1, amux, amuy, k1, dx, dy, trpt)
    if iv != KY_K1_QUAD:
        return dtrans, dcod
    if al == 0.0:
        # qdthin (薄い四極の微分) は未移植
        return dtrans, dcod
    dtrans, cod = qtentu(u_k, True)
    qchg(dtrans, cod, -dx, -dy, theta, True)
    pr = 1.0 + cod[5]
    xi, pxi, yi, pyi = cod[0], cod[1], cod[2], cod[3]
    akk = math.sqrt(abs(ak / pr / al))
    phi = akk * al
    if ak * al * pr > 0.0:
        dphi = 0.5 * phi / ak
        sinphi, cosphi = math.sin(phi), math.cos(phi)
        a1 = cosphi + sinphi / phi
        a2 = (1.0 - phi) * (1.0 + phi) / phi * sinphi
        a11 = -sinphi * dphi
        a12 = (cosphi - sinphi / phi) * dphi / akk
        a21 = -a1 * dphi * akk
        shphi, chphi = math.sinh(phi), math.cosh(phi)
        b1 = chphi + shphi / phi
        b2 = (1.0 + phi**2) * shphi / phi
        b11 = shphi * dphi
        b12 = (chphi - shphi / phi) * dphi / akk
        b21 = b1 * dphi * akk
        a16 = 0.25 * al / pr * (a1 * xi / pr + (-a2 + cosphi) / ak * pxi)
        a26 = 0.25 / pr * ((a2 + 3.0 * cosphi) * xi / pr + a1 * al * pxi)
        b16 = -0.25 * al / pr * (b1 * yi / pr - (-b2 + chphi) / ak * pyi)
        b26 = -0.25 / pr * ((b2 + 3.0 * chphi) * yi / pr + b1 * al * pyi)
    elif ak * al * pr < 0.0:
        dphi = 0.5 * phi / ak
        sinphi, cosphi = math.sin(phi), math.cos(phi)
        a1 = cosphi + sinphi / phi
        a2 = (1.0 - phi) * (1.0 + phi) / phi * sinphi
        b11 = -sinphi * dphi
        b12 = (cosphi - sinphi / phi) * dphi / akk
        b21 = -a1 * dphi * akk
        shphi, chphi = math.sinh(phi), math.cosh(phi)
        b1 = chphi + shphi / phi
        b2 = (1.0 + phi**2) * shphi / phi
        a11 = shphi * dphi
        a12 = (chphi - shphi / phi) * dphi / akk
        a21 = b1 * dphi * akk
        b16 = -0.25 * al / pr * (a1 * yi / pr - (-a2 + cosphi) / ak * pyi)
        b26 = -0.25 / pr * ((a2 + 3.0 * cosphi) * yi / pr + a1 * al * pyi)
        a16 = 0.25 * al * (b1 * xi / pr + (-b2 + chphi) / ak * pxi)        # SAD のまま (/pr 無し)
        a26 = 0.25 / pr * ((b2 + 3.0 * chphi) * xi / pr + b1 * al * pxi)
    else:
        a11 = al * 0.5
        a12 = -al**2 / 6.0
        a21 = -1.0
        b11 = -al * 0.5
        b12 = al**2 / 6.0
        b21 = -1.0
        a16 = -(al * a21 * xi + (-a12 + al * a11) * pxi) * 0.5 / pr
        a26 = -((a21 - ak * a11 - 1.0) * xi + al * a21 * pxi) * 0.5 / pr
        b16 = -(al * b21 * yi + (-b12 + al * b11) * pyi) * 0.5 / pr
        b26 = -((b21 + ak * b11 + 1.0) * yi + al * b21 * pyi) * 0.5 / pr
    for i in range(5):
        x = dtrans[0, i]
        dtrans[0, i] = a11 * x + a12 * dtrans[1, i]
        dtrans[1, i] = a21 * x + a11 * dtrans[1, i]
        y = dtrans[2, i]
        dtrans[2, i] = b11 * y + b12 * dtrans[3, i]
        dtrans[3, i] = b21 * y + b11 * dtrans[3, i]
    dtrans[0, 4] += a16
    dtrans[1, 4] += a26
    dtrans[2, 4] += b16
    dtrans[3, 4] += b26
    dcod[0] = a11 * xi + a12 * pxi
    dcod[1] = a21 * xi + a11 * pxi
    dcod[2] = b11 * yi + b12 * pyi
    dcod[3] = b21 * yi + b11 * pyi
    dcod[4] = 0.0
    qchg(dtrans, dcod, 0.0, 0.0, -theta, False)
    return dtrans, dcod


def qdrotate(u_k: np.ndarray, u_k1: np.ndarray, amux: float, amuy: float, k1: int,
             dx: float, dy: float, trpt: bool) -> Tuple[np.ndarray, np.ndarray]:
    """qdquad.f qdrotate: QUAD の ROTATE による微分 (要素 k1→k1+1 の転送行列を回転の生成子で挟む)"""
    trans2, cod = qtentu(u_k, True)
    _, cod2 = qtentu(u_k1, True)
    cod[0] -= dx
    cod[2] -= dy
    trans1 = tftmatu(u_k, u_k1, amux, amuy, k1, k1 + 1, False, True)
    dtrans = np.zeros((4, 5))
    dcod1 = np.zeros(4)
    for i in range(4):
        dtrans[i, 2] = -trans1[i, 0]
        dtrans[i, 0] = trans1[i, 2]
        dtrans[i, 3] = -trans1[i, 1]
        dtrans[i, 1] = trans1[i, 3]
        dtrans[i, 4] = 0.0
        dcod1[i] = cod2[i] - (trans1[i, 0] * cod[0] + trans1[i, 1] * cod[1]
                              + trans1[i, 2] * cod[2] + trans1[i, 3] * cod[3])
    dcod1[0] -= dx
    dcod1[2] -= dy
    for i in range(5):
        dtrans[0, i] += trans1[2, i]
        dtrans[2, i] -= trans1[0, i]
        dtrans[1, i] += trans1[3, i]
        dtrans[3, i] -= trans1[1, i]
    dcod = np.zeros(6)
    for i in range(4):
        dcod[i] = dtrans[i, 0] * cod[0] + dtrans[i, 1] * cod[1] + dtrans[i, 2] * cod[2] + dtrans[i, 3] * cod[3]
    dcod[0] += dcod1[2]
    dcod[1] += dcod1[3]
    dcod[2] -= dcod1[0]
    dcod[3] -= dcod1[1]
    for i in range(4):
        x, px, y, py = dtrans[i, 0], dtrans[i, 1], dtrans[i, 2], dtrans[i, 3]
        dtrans[i, 0] = x * trans2[0, 0] + px * trans2[1, 0] + y * trans2[2, 0] + py * trans2[3, 0]
        dtrans[i, 1] = x * trans2[0, 1] + px * trans2[1, 1] + y * trans2[2, 1] + py * trans2[3, 1]
        dtrans[i, 2] = x * trans2[0, 2] + px * trans2[1, 2] + y * trans2[2, 2] + py * trans2[3, 2]
        dtrans[i, 3] = x * trans2[0, 3] + px * trans2[1, 3] + y * trans2[2, 3] + py * trans2[3, 3]
        dtrans[i, 4] = x * trans2[0, 4] + px * trans2[1, 4] + y * trans2[2, 4] + py * trans2[3, 4] + dtrans[i, 4]
    return dtrans, dcod


# ---------------------------------------------------------------------------
# qdbend.f qdbend
# ---------------------------------------------------------------------------
def qdbend(al: float, phib0: float, phi0: float, psi1: float, psi2: float, ak: float, u: np.ndarray,
           dx: float, dy: float, theta: float, iv: int) -> Tuple[np.ndarray, np.ndarray]:
    """BEND の微分行列。iv=2: ANGLE、iv=8: K1。psi1, psi2 は SAD の E1, E2 (ANGLE に対する比)。
    phib0 = ANGLE + K0 (実際の曲げ角)、phi0 = ANGLE (設計角)、ak = K1 (積分値)"""
    dtrans = np.zeros((4, 5))
    dcod = np.zeros(6)
    if al == 0.0:
        return dtrans, dcod
    if phi0 == 0.0 or phib0 == 0.0:
        phi = 1e-7
        phib = 1e-7
    else:
        phi = phi0
        phib = phib0
    dx1 = dpx1 = 0.0
    a12 = a16 = a22 = a26 = a34 = a44 = a46 = a36 = 0.0
    dtrans, cod = qtentu(u, True)
    qchg(dtrans, cod, -dx, -dy, theta, True)
    rhob = al / phib
    rho0 = al / phi
    aind = ak * rho0 / phi
    if aind == -1.0:
        aind = -1.0 + 1e-8
    ax = math.sqrt(abs(1.0 + aind))
    ay = math.sqrt(abs(aind))
    pr = 1.0 + u[DDP]
    rho = rhob * pr
    tanp1 = math.tan(psi1 * phi0)
    tanp2 = math.tan(psi2 * phi0)
    ak1 = -tanp1 / rho
    ak2 = -tanp2 / rho
    xi = cod[0]
    pxi = cod[1] - ak1 * xi
    yi = cod[2]
    pyi = cod[3] + ak1 * yi
    rhos = math.sqrt(rho * rho0)
    if ax != 0.0:
        rhox = rhos / ax
        phix = al / rhox
        rhoa = rho / (1.0 + aind)
    else:
        phix = 0.0
        rhox = 1e19
        rhoa = 1e19
    if ay != 0.0:
        rhoy = rhos / ay
        phiy = al / rhoy
    else:
        phiy = 0.0
        rhoy = 1e19
    if iv == KY_ANGL_BEND:
        bpr = rho / rho0
        if aind > -1.0:
            dphix = phib / pr / phix
            csphix, snphix = math.cos(phix), math.sin(phix)
            sinsq = snphix**2 / (1.0 + csphix) if csphix >= 0.0 else 1.0 - csphix
            scphix = _sinc(phix)
            xcs = phix * csphix + snphix
            drho = (rho - rho0) / (1.0 + aind)
            dx1 = drho * ((aind - 1.0) * sinsq / phix + snphix) * dphix
            dpx1 = drho * (xcs + (aind - 1.0) * snphix) * dphix / al
            xd = xi - drho
            t11 = csphix
            t12 = rhox * snphix
            t21 = -snphix / rhox
            t16 = ((phix * snphix * xd - rhox * scphix * pxi) * 0.5 + rhoa * sinsq) / pr
            if aind > 0.0:
                t33 = math.cosh(phiy)
                t34 = rhoy * math.sinh(phiy)
            elif aind == 0.0:
                t33 = 1.0
                t34 = al
            else:
                t33 = math.cos(phiy)
                t34 = rhoy * math.sin(phiy)
            b11 = -snphix * dphix
            b12 = rhox / phix * scphix * dphix
            b21 = -xcs * dphix / al
            b16 = (drho * (1.0 - aind) * snphix
                   + 2.0 * rhoa * ((aind - 1.0) / phix * sinsq + snphix)
                   + xd * xcs
                   + pxi * rhox / phix * (scphix + phix**2 * snphix)) * dphix * 0.5 / pr
            b26 = (pr * rhoa * (aind - 1.0) * snphix
                   + (pr * (rhoa + pxi * al * 0.5) + drho * (1.0 - aind) * 0.5) * xcs
                   + xd * (3.0 * phix * csphix + (1.0 - phix**2) * snphix) * 0.5) * dphix / al / pr
        elif aind == -1.0:
            snphix = phib / pr
            dx1 = (bpr - 1.0) / bpr * al
            dpx1 = dx1 / al
            t11 = 1.0
            t12 = al
            t21 = 0.0
            t16 = 0.0
            t33 = math.cos(phiy)
            t34 = rhoy * math.sin(phiy)
            b11 = -snphix
            b12 = 0.0
            b21 = 0.0
            b16 = -2.0 * bpr / pr * phi * al
            b26 = b16 / al
        else:
            dphix = phib / pr / phix
            chphix, shphix = math.cosh(phix), math.sinh(phix)
            sinsq = shphix**2 / (1.0 + chphix)
            scphix = _sinhc(phix)
            xcs = phix * chphix + shphix
            drho = (rho - rho0) / (1.0 + aind)
            dx1 = drho * ((aind - 1.0) * sinsq / phix + shphix) * dphix
            dpx1 = drho * (xcs + (aind - 1.0) * shphix) * dphix / al
            xd = xi - drho
            t11 = chphix
            t12 = rhox * shphix
            t21 = shphix / rhox
            t16 = ((-phix * shphix * xd - rhox * scphix * pxi) * 0.5 - rhoa * sinsq) / pr
            t33 = math.cos(phiy)
            t34 = rhoy * math.sin(phiy)
            b11 = -shphix * dphix
            b12 = -rhox / phix * scphix * dphix
            b21 = -xcs * dphix / al
            b16 = (drho * (1.0 - aind) * shphix
                   + 2.0 * rhoa * ((aind - 1.0) / phix * sinsq + shphix)
                   + xd * xcs
                   - pxi * rhox / phix * (scphix - phix**2 * shphix)) * dphix * 0.5 / pr
            b26 = (pr * rhoa * (aind - 1.0) * shphix
                   + (pr * (rhoa + pxi * al * 0.5) + drho * (1.0 - aind) * 0.5) * xcs
                   + xd * (3.0 * phix * chphix + (1.0 + phix**2) * shphix) * 0.5) * dphix / al / pr
        dk1 = 2.0 * psi1 / rho
        dk2 = 2.0 * psi2 / rho
        a11 = b11 + dk1 * t12
        a12 = b12
        a21 = b21 + (dk1 + dk2) * t11
        a22 = b11 + dk2 * t12
        a33 = -dk1 * t34
        a34 = 0.0
        a43 = -t33 * (dk1 + dk2)
        a44 = -dk2 * t34
        a16 = b16
        a26 = b26 + dk2 * t16
        a36 = 0.0
        a46 = 0.0
    elif iv == KY_K1_BEND:
        bpr = rho / rho0
        if aind > 0.0:
            dphix = al / 2.0 / phix / bpr
            dphiy = -al / 2.0 / phiy / bpr
            csphix, snphix = math.cos(phix), math.sin(phix)
            sinsq = snphix**2 / (1.0 + csphix) if csphix >= 0.0 else 1.0 - csphix
            chphiy, shphiy = math.cosh(phiy), math.sinh(phiy)
            scphix = _sinc(phix)
            scphiy = _sinhc(phiy)
            xcs = phix * csphix + snphix
            ycs = phiy * chphiy + shphiy
            drho = (rho - rho0) / (1.0 + aind)
            ddrho = -rho0 / phi * drho / (1.0 + aind)
            dx1 = ddrho * sinsq + drho * snphix * dphix
            dpx1 = ddrho * snphix / rhox + drho / al * xcs * dphix
            a11 = -snphix * dphix
            a12 = rhox / phix * scphix * dphix
            a21 = -xcs / al * dphix
            a16 = (drho * snphix
                   + rhoa * (-2.0 * sinsq / phix + snphix)
                   - (drho - xi) * xcs * 0.5
                   + pxi * rhox / phix * (scphix + phix**2 * snphix) * 0.5) * dphix / pr
            a26 = (-rhoa * 2.0 * snphix
                   + (rhoa + pxi * al * 0.5 + drho) * xcs
                   - (drho - xi) * 0.5 * (3.0 * phix * csphix + (1.0 - phix**2) * snphix)) / al / pr * dphix
            a33 = -shphiy * dphiy
            a34 = -rhoy / phiy * scphiy * dphiy
            a43 = -(shphiy + phiy * chphiy) / al * dphiy
            a36 = (yi * ycs - pyi * rhoy / phiy * (scphiy - phiy**2 * shphiy)) * 0.5 / pr * dphiy
            a46 = (pyi * ycs - yi * rhoy / phiy * (3.0 * phiy * chphiy + (1.0 + phiy**2) * shphiy)) * 0.5 / pr * dphiy
        elif aind == 0.0:
            uu = phi**2 / al
            dphix = phix / 2.0 / uu
            csphix, snphix = math.cos(phix), math.sin(phix)
            sinsq = snphix**2 / (1.0 + csphix) if csphix >= 0.0 else 1.0 - csphix
            scphix = _sinc(phix)
            drho = (rho - rho0)
            ddrho = -rho0 / phi * drho
            dx1 = ddrho * sinsq + drho * snphix * dphix
            dpx1 = ddrho * snphix / rhox + drho / al * (snphix + phix * csphix) * dphix
            a11 = -snphix * dphix
            a12 = rhox / phix * scphix * dphix
            a21 = -(snphix + phix * csphix) / al * dphix
            a16 = ((-2.0 * sinsq / phix + phix * (rho0 - rho) / rho * csphix / 2.0
                    + (3.0 * rho - rho0) / rho * snphix / 2.0) * rhox**2 / rho0
                   + xi / 2.0 * (phix * csphix + snphix)
                   + pxi * rhox / 2.0 * (scphix / phix + phix * snphix)) * dphix / pr
            a26 = ((scphix / phix * (1.0 + rho0 / rho) + (1.0 - rho0 / rho) * phix * snphix) * rhox / rho0
                   + xi * (3.0 * csphix + (1.0 / phix - phix) * snphix) / rhox
                   + pxi * (phix * csphix + snphix)) * dphix / 2.0 / pr
            a33 = -rho0 * al / rho / 2.0
            a34 = rho0 / rho / 6.0 * al**2
            a43 = -rho0 / rho
            a36 = (yi * 2.0 + pyi * al / 1.5) * rho0 / rho / 4.0 / pr
            a46 = -(yi * 4.0 / al + pyi * 2.0) * rho0 / rho / 4.0 / pr
        elif aind > -1.0:
            dphix = al / 2.0 / phix / bpr
            dphiy = -al / 2.0 / phiy / bpr
            csphix, snphix = math.cos(phix), math.sin(phix)
            sinsq = snphix**2 / (1.0 + csphix) if csphix >= 0.0 else 1.0 - csphix
            csphiy, snphiy = math.cos(phiy), math.sin(phiy)
            scphix = _sinc(phix)
            scphiy = _sinc(phiy)
            xcs = phix * csphix + snphix
            ycs = phiy * csphiy + snphiy
            drho = (rho - rho0) / (1.0 + aind)
            ddrho = -rho0 / phi * drho / (1.0 + aind)
            dx1 = ddrho * sinsq + drho * snphix * dphix
            dpx1 = ddrho * snphix / rhox + drho / al * xcs * dphix
            a11 = -snphix * dphix
            a12 = rhox / phix * scphix * dphix
            a21 = -xcs / al * dphix
            a16 = (drho * snphix
                   + rhoa * (-2.0 * sinsq / phix + snphix)
                   - (drho - xi) * xcs * 0.5
                   + pxi * rhox / phix * (scphix + phix**2 * snphix) * 0.5) * dphix / pr
            a26 = (-rhoa * 2.0 * snphix
                   + (rhoa + pxi * al * 0.5 + drho) * xcs
                   - (drho - xi) * 0.5 * (3.0 * phix * csphix + (1.0 - phix**2) * snphix)) / al / pr * dphix
            a33 = -snphiy * dphiy
            a34 = rhoy / phiy * scphiy * dphiy
            a43 = -(snphiy + phiy * csphiy) / al * dphiy
            a36 = (yi * ycs + pyi * rhoy / phiy * (scphiy + phiy**2 * snphiy)) * 0.5 / pr * dphiy
            a46 = (pyi * ycs + yi * rhoy / phiy * (3.0 * phiy * csphiy + (1.0 - phiy**2) * snphiy)) * 0.5 / pr * dphiy
        else:
            dphix = al / 2.0 / phix / bpr
            dphiy = -al / 2.0 / phiy / bpr
            chphix, shphix = math.cosh(phix), math.sinh(phix)
            sinsq = shphix**2 / (1.0 + chphix)
            csphiy, snphiy = math.cos(phiy), math.sin(phiy)
            scphix = _sinhc(phix)
            scphiy = _sinc(phiy)
            xcs = phix * chphix + shphix
            ycs = phiy * csphiy + snphiy
            drho = (rho - rho0) / (1.0 + aind)
            ddrho = -rho0 / phi * drho / (1.0 + aind)
            dx1 = -ddrho * sinsq + drho * shphix * dphix
            dpx1 = -ddrho * shphix / rhox + drho / al * xcs * dphix
            a11 = -shphix * dphix
            a12 = -rhox / phix * scphix * dphix
            a21 = -xcs / al * dphix
            a16 = (drho * shphix
                   + rhoa * (-2.0 * sinsq / phix + shphix)
                   - (drho - xi) * xcs * 0.5
                   - pxi * rhox / phix * (scphix - phix**2 * shphix) * 0.5) * dphix / pr
            a26 = (-rhoa * 2.0 * shphix
                   + (rhoa + pxi * al * 0.5 + drho) * xcs
                   - (drho - xi) * 0.5 * (3.0 * phix * chphix + (1.0 + phix**2) * shphix)) / al / pr * dphix
            a33 = -snphiy * dphiy
            a34 = rhoy / phiy * scphiy * dphiy
            a43 = -(snphiy + phiy * csphiy) / al * dphiy
            a36 = (yi * ycs + pyi * rhoy / phiy * (scphiy + phiy**2 * snphiy)) * 0.5 / pr * dphiy
            a46 = (pyi * ycs + yi * rhoy / phiy * (3.0 * phiy * csphiy + (1.0 - phiy**2) * snphiy)) * 0.5 / pr * dphiy
        a22 = a11
        a44 = a33
    else:
        raise ValueError(f'qdbend: invalid iv={iv}')
    dtrans[1, :] -= ak1 * dtrans[0, :]
    dtrans[3, :] += ak1 * dtrans[2, :]
    x = dtrans[0, :].copy()
    dtrans[0, :] = a11 * x + a12 * dtrans[1, :]
    dtrans[1, :] = a21 * x + a22 * dtrans[1, :]
    y = dtrans[2, :].copy()
    dtrans[2, :] = a33 * y + a34 * dtrans[3, :]
    dtrans[3, :] = a43 * y + a44 * dtrans[3, :]
    dtrans[0, 4] += a16
    dtrans[1, 4] += a26
    dtrans[2, 4] += a36
    dtrans[3, 4] += a46
    dcod[0] = a11 * xi + a12 * pxi + dx1
    dcod[1] = a21 * xi + a22 * pxi + dpx1
    dcod[2] = a33 * yi + a34 * pyi
    dcod[3] = a43 * yi + a44 * pyi
    dcod[5] = pr
    dtrans[1, :] -= ak2 * dtrans[0, :]
    dtrans[3, :] += ak2 * dtrans[2, :]
    dcod[1] -= ak2 * dcod[0]
    dcod[3] += ak2 * dcod[2]
    qchg(dtrans, dcod, 0.0, 0.0, -theta, False)
    return dtrans, dcod


# ---------------------------------------------------------------------------
# qdmdia.f qdmdia
# ---------------------------------------------------------------------------
def qdmdia(tm: np.ndarray, dtm: np.ndarray, r1: float) -> Tuple[float, float, float, float]:
    """4×4 一周行列 tm とその微分 dtm から、対角化行列 S = ((c, J·rᵗ·J), (r, c)) の r の微分 dr1..dr4 を返す。
    tm = S⁻¹·diag(A, B)·S。通常は (TrA−TrB)² > 0 の分岐。縮退分岐も SAD のまま (1300 の dr4 未設定の癖も含む)"""
    eps = 1e-7
    (tm11, tm12, tm13, tm14), (tm21, tm22, tm23, tm24), (tm31, tm32, tm33, tm34), (tm41, tm42, tm43, tm44) = tm
    (dtm11, dtm12, dtm13, dtm14), (dtm21, dtm22, dtm23, dtm24), (dtm31, dtm32, dtm33, dtm34), \
        (dtm41, dtm42, dtm43, dtm44) = dtm
    trpq = (tm11 + tm22 - tm33 - tm44)
    dtrpq = (dtm11 + dtm22 - dtm33 - dtm44)
    trab2 = 4.0 * ((tm31 + tm24) * (tm42 + tm13) + (tm14 - tm32) * (tm41 - tm23)) + trpq**2
    dtrab2 = 4.0 * ((dtm31 + dtm24) * (tm42 + tm13) + (dtm14 - dtm32) * (tm41 - tm23)
                    + (tm31 + tm24) * (dtm42 + dtm13) + (tm14 - tm32) * (dtm41 - dtm23)) + 2.0 * trpq * dtrpq
    degenerate = False
    if trab2 < eps:
        if abs(dtrab2) < eps:
            degenerate = True
        else:
            trab2 = eps
    if not degenerate:
        c2 = 0.5 * abs(trpq) / math.sqrt(trab2)
        c = math.sqrt(0.5 + c2)
        dc = 0.5 * c2 * (dtrpq / trpq - 0.5 * dtrab2 / trab2) / c
        if r1 * (tm31 + tm24) <= 0.0:
            aa = -1.0 / c / math.sqrt(trab2)
        else:
            aa = 1.0 / c / math.sqrt(trab2)
        da = aa * (-dc / c - 0.5 * dtrab2 / trab2)
        dr1 = da * (tm31 + tm24) + aa * (dtm31 + dtm24)
        dr2 = da * (tm32 - tm14) + aa * (dtm32 - dtm14)
        dr3 = da * (tm41 - tm23) + aa * (dtm41 - dtm23)
        dr4 = da * (tm42 + tm13) + aa * (dtm42 + dtm13)
        return dr1, dr2, dr3, dr4
    # ---- (TrA−TrB)² = 0 (label 100)
    detn = tm13 * tm24 - tm14 * tm23
    ddetn = dtm13 * tm24 + tm13 * dtm24 - dtm14 * tm23 - tm14 * dtm23
    detn_zero = False
    if abs(detn) < eps:
        if abs(ddetn) < 1.0:
            detn_zero = True
        else:
            detn = _sign(eps, detn)
    if not detn_zero:
        trqk = tm33 - tm44
        dtrqk = dtm33 - dtm44
        goto400 = False
        if abs(trqk) < eps:
            if abs(dtrqk) < eps:
                goto400 = True
            else:
                trqk = _sign(eps, trqk)
        if not goto400:
            trqj = tm43 - tm34
            dtrqj = dtm43 - dtm34
            trqkj2 = (trqk + trqj) * (trqk - trqj)
            dtrqkj = (dtrqk + dtrqj) * (trqk - trqj) + (trqk + trqj) * (dtrqk - dtrqj)
            goto200 = False
            if abs(trqkj2) < eps:
                if abs(dtrqkj) < eps:
                    goto200 = True
                else:
                    trqkj2 = eps
            if not goto200:
                bbb = -detn / (trqj + trqk)
                dbb = -((dtrqj + dtrqk) * bbb + ddetn) / (trqj + trqk)
            else:
                if abs(trqj) < eps:
                    trqj = _sign(eps, trqj)
                bbb = -0.5 * detn / trqj
                dbb = -(dtrqj * bbb + 0.5 * ddetn) / trqj
            aaa = (-detn - bbb * trqj) / trqk
            daa = (-dtrqk * aaa - ddetn - dbb * trqj - bbb * dtrqj) / trqk
            ccc = 0.0
            dcc = 0.0
        else:
            trqj = tm43 - tm34
            dtrqj = dtm43 - dtm34
            goto500 = False
            if abs(trqj) < eps:
                if abs(dtrqj) < eps:
                    goto500 = True
                else:
                    trqj = _sign(eps, trqj)
            if not goto500:
                aaa = detn / trqj
                bbb = -detn / trqj
                ccc = 0.0
                daa = (-dtrqj * aaa + ddetn) / trqj
                dbb = -(dtrqj * bbb + ddetn) / trqj
                dcc = 0.0
            else:
                trql = tm34 + tm43
                dtrql = dtm34 + dtm43
                goto600 = False
                if abs(trql) < eps:
                    if abs(dtrql) < eps:
                        goto600 = True
                    else:
                        trql = _sign(eps, trql)
                if not goto600:
                    aaa = 0.0
                    bbb = detn / trql
                    ccc = -bbb
                    daa = 0.0
                    dbb = (-dtrql * bbb + ddetn) / trql
                    dcc = -dbb
                else:
                    if detn > 0.0:
                        aaa = 0.0
                        bbb = math.sqrt(detn / 2.0)
                        ccc = 0.0
                        daa = 0.0
                        dbb = 0.5 * ddetn / math.sqrt(2.0 * detn)
                        dcc = 0.0
                    else:
                        aaa = math.sqrt(-detn / 2.0)
                        bbb = 0.0
                        ccc = 0.0
                        daa = -0.5 * ddetn / math.sqrt(-2.0 * detn)
                        dbb = 0.0
                        dcc = 0.0
        # label 700
        xx11, xx12, xx21, xx22 = aaa, bbb + ccc, -bbb + ccc, -aaa
        dx11, dx12, dx21, dx22 = daa, dbb + dcc, -dbb + dcc, -daa
        dni11 = -ddetn * tm24 / detn + dtm24
        dni12 = ddetn * tm14 / detn - dtm14
        dni21 = ddetn * tm23 / detn - dtm23
        dni22 = -ddetn * tm13 / detn + dtm13
        dr1 = (dx11 * tm24 - dx12 * tm23 + xx11 * dni11 + xx12 * dni21) / detn
        dr2 = (-dx11 * tm14 + dx12 * tm13 + xx11 * dni12 + xx12 * dni22) / detn
        dr3 = (dx21 * tm24 - dx22 * tm23 + xx21 * dni11 + xx22 * dni21) / detn
        dr4 = (-dx21 * tm14 + dx22 * tm13 + xx21 * dni12 + xx22 * dni22) / detn
        return dr1, dr2, dr3, dr4
    # ---- det(N) = 0 (label 800)
    trqk = tm33 - tm44
    trqj = tm43 - tm34
    trql = tm43 + tm34
    dtrqk = dtm33 - dtm44
    dtrqj = dtm43 - dtm34
    dtrql = dtm43 + dtm34
    sel = None
    if abs(trqk) >= eps:
        sel = 'k'
    elif abs(trqj) >= eps:
        sel = 'j'
    elif abs(trql) >= eps:
        sel = 'l'
    elif abs(dtrqk) >= eps:
        trqk = _sign(eps, trqk)
        sel = 'k'
    elif abs(dtrqj) >= eps:
        trqj = _sign(eps, trqj)
        sel = 'j'
    elif abs(dtrql) >= eps:
        trql = _sign(eps, trql)
        sel = 'l'
    if sel is not None:
        if sel == 'k':
            aaa, bbb, ccc = 0.0, 0.0, 1.0 / trqk
            daa, dbb, dcc = 0.0, 0.0, -dtrqk * ccc / trqk
        elif sel == 'j':
            aaa, bbb, ccc = -1.0 / trqj, 0.0, 0.0
            daa, dbb, dcc = -dtrqj * aaa / trqj, 0.0, 0.0
        else:
            aaa, bbb, ccc = 0.0, -1.0 / trql, 0.0
            daa, dbb, dcc = 0.0, -dtrql * bbb / trql, 0.0
        yy11, yy12, yy21, yy22 = aaa + bbb, ccc, ccc, aaa - bbb
        dy11, dy12, dy21, dy22 = daa + dbb, dcc, dcc, daa - dbb
        dr1 = dy11 * (-tm23) + dy12 * (-tm24) + yy11 * (-dtm23) + yy12 * (-dtm24)
        dr2 = dy11 * tm13 + dy12 * tm14 + yy11 * dtm13 + yy12 * dtm14
        # SAD (qdmdia.f label 1300) は dr3 を 2 度代入し dr4 を設定しない。ここでは dr4=0 とする
        dr3 = dy21 * tm13 + dy22 * tm14 + yy21 * dtm13 + yy22 * dtm14
        dr4 = 0.0
        return dr1, dr2, dr3, dr4
    # ---- Tr(QJ) = Tr(QL) = Tr(QK) = 0 (label 1400)
    trp = tm11 + tm22
    p1 = tm11 - 0.5 * trp
    p3 = tm21
    dtrp = dtm11 + dtm22
    dp1 = dtm11 - 0.5 * dtrp
    dp3 = dtm21
    if abs(p3) < eps:
        if abs(dp3) < eps:
            if abs(p1) < eps:
                if abs(dp1) < eps:
                    return 0.0, 0.0, 0.0, 0.0
                p1 = _sign(eps, p1)
            dr1 = -tm24 * dp1 / p1 / p1 + dtm24 / p1
            dr3 = tm23 * dp1 / p1 / p1 - dtm23 / p1
            return dr1, 0.0, dr3, 0.0
        p3 = _sign(eps, p3)
    dr2 = -dp3 * tm24 / p3 / p3 + dtm24 / p3
    dr4 = dp3 * tm23 / p3 / p3 - dtm23 / p3
    return 0.0, dr2, 0.0, dr4


# ---------------------------------------------------------------------------
# qdcell.f dgammag
# ---------------------------------------------------------------------------
def dgammag(da: float, db: float, a: float, b: float, a1: float, b1: float) -> Tuple[float, float]:
    """(dγ/γ 相当, dBMAG): da = dα/(1+α²), db = dβ/β。a1, b1 は参照 Twiss (BMAG 用)"""
    dg = 2.0 * a * da - db
    g = (1.0 + a**2) / b
    g1 = (1.0 + a1**2) / b1
    dbmag = 0.5 * (b1 * g * dg + g1 * b * db) - a1 * g * b * da
    return dg, dbmag


# ---------------------------------------------------------------------------
# 要素情報 → 微分行列 (qdtwis の select case 部分)
# ---------------------------------------------------------------------------
class DerivContext:
    """qdcell/qdtwis が参照する SAD のモジュール変数の束: utwiss (位置 l → Twiss)、nlat、nut、cell、trpt と要素列"""

    def __init__(self, elements, u: np.ndarray, cell: bool = True):
        self.elements = elements           # pysad の要素列 (SAD の位置 k = elements[k-1])
        self.u = u                         # utwiss[l, mfit]  (l = 1..nlat)
        self.nlat = len(elements) + 1
        self.nut = self.nlat               # itwissp(nlat)
        self.cell = cell
        self.trpt = not cell

    def element_dtrans(self, k: int, iv: int) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """要素 k (1 始まり) の変数 iv に対する微分行列 (dtrans, dcod)。未対応なら None (数値微分へ)"""
        el = self.elements[k - 1]
        tn = type(el).__name__
        u = self.u
        if tn == 'Quadrupole':
            if iv in (KY_K1_QUAD, KY_ROT_QUAD):
                al = float(el.L)
                ak = _sad_get(el, 'K1')
                theta = float(getattr(el, 'ROTATION', 0.0) or 0.0)
                return qdquad(al, ak, u[k], u[k + 1], u[self.nut, NX], u[self.nut, NY], k,
                              float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0),
                              theta, iv, self.trpt)
            return None
        if tn == 'Bend':
            if iv in (KY_ANGL_BEND, KY_K1_BEND):
                al = float(el.L)
                angle = float(el.ANGLE)
                k0 = float(getattr(el, 'K0', 0.0) or 0.0)
                ak = _sad_get(el, 'K1')
                # pysad は E1/E2 を rad で持つ (parser.py で変換) → SAD の E1 (ANGLE 比) に戻す
                e1 = float(getattr(el, 'E1', 0.0) or 0.0)
                e2 = float(getattr(el, 'E2', 0.0) or 0.0)
                psi1 = e1 / angle if angle != 0.0 else 0.0
                psi2 = e2 / angle if angle != 0.0 else 0.0
                theta = float(getattr(el, 'ROTATION', 0.0) or 0.0) + float(getattr(el, 'ROTATE', 0.0) or 0.0)   # SAD の ROT (設計ロール ROTATE + ROTATION)
                return qdbend(al, angle + k0, angle, psi1, psi2, ak, u[k],
                              float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0),
                              theta, iv)
            return None
        if tn == 'Drift' and getattr(el, '_sad_type', '') == 'BEND':
            # ANGLE=0 の BEND (VD/HD): qdbend は al=0 または L≠0 でも角度 1e-7 として扱う
            if iv in (KY_ANGL_BEND, KY_K1_BEND):
                return qdbend(float(el.L), 0.0, 0.0, 0.0, 0.0, 0.0, u[k], 0.0, 0.0, 0.0, iv)
            return None
        return None


# ---------------------------------------------------------------------------
# qdtwis.f qdtwis (輸送の微分)
# ---------------------------------------------------------------------------
def qdtwis(ctx: DerivContext, k0: int, l: int, iv: int) -> np.ndarray:
    """要素 k0 の変数 iv による位置 l の Twiss の微分 (k0 の入口 Twiss は固定とみなす = INS 的)。
    戻り dtwiss[1..mfittry]: AX は dα/(1+α²)、BX は dβ/β、NX は dψ (rad)、EX.. は dη、R1.. は dR、DX.. は d(閉軌道)"""
    u = ctx.u
    nut = ctx.nut
    trpt = ctx.trpt
    dtwiss = np.zeros(NDT)
    k = k0
    res = ctx.element_dtrans(k, iv)
    if res is None:
        return dtwiss
    dtrans, dcod = res
    dcod[4] = 0.0 if len(dcod) > 4 and not np.isfinite(dcod[4]) else dcod[4]
    k1 = k + 1
    gr = 1.0                                   # sqrt(gammab(k1)/gammab(l)): 加速なし
    dcod2 = np.zeros(6)
    dcod2[0:4] = dcod[0:4] * gr
    dcod2[4] = dcod[4]
    dcod2[5] = dcod[5] * gr**2
    dtrans[:, 4] = dtrans[:, 4] / gr
    trans = tftmatu(u[k1], u[l], u[nut, NX], u[nut, NY], k1, l, False, trpt)
    trans1 = trans[:, 0:4] @ dtrans
    dcod1 = np.zeros(6)
    dcod1[0:4] = trans[:, 0:4] @ dcod2[0:4]
    dcod1[4] = dcod2[4]
    # qddtwiss (閉軌道の数値伝搬) は CELL では DZ/DDP 条件のときだけ → 未移植
    trans = qgettru(u[k], u[l], u[nut, NX], u[nut, NY], k, l, True, True, trpt)
    r1, r2, r3, r4 = u[l, R1], u[l, R2], u[l, R3], u[l, R4]
    detp = r1 * r4 - r2 * r3
    sqrdet = math.sqrt(1.0 - detp)
    ddetp = trans1[0, 0] * trans[1, 1] + trans[0, 0] * trans1[1, 1] - trans1[0, 1] * trans[1, 0] - trans[0, 1] * trans1[1, 0]
    normal = u[l, DETR] < 1.0
    if normal:
        dsqr = 0.5 * ddetp / sqrdet
        x11, x12, x21, x22 = trans[0, 0] / sqrdet, trans[0, 1] / sqrdet, trans[1, 0] / sqrdet, trans[1, 1] / sqrdet
        dx11 = (trans1[0, 0] - x11 * dsqr) / sqrdet
        dx12 = (trans1[0, 1] - x12 * dsqr) / sqrdet
        dx21 = (trans1[1, 0] - x21 * dsqr) / sqrdet
        dx22 = (trans1[1, 1] - x22 * dsqr) / sqrdet
        dr1 = -trans1[2, 0] * x22 - trans[2, 0] * dx22 + trans1[2, 1] * x21 + trans[2, 1] * dx21
        dr2 = trans1[2, 0] * x12 + trans[2, 0] * dx12 - trans1[2, 1] * x11 - trans[2, 1] * dx11
        dr3 = -trans1[3, 0] * x22 - trans[3, 0] * dx22 + trans1[3, 1] * x21 + trans[3, 1] * dx21
        dr4 = trans1[3, 0] * x12 + trans[3, 0] * dx12 - trans1[3, 1] * x11 - trans[3, 1] * dx11
        dtwiss[EX] = (trans[0, 4] * dsqr + trans1[0, 4] * sqrdet - dr4 * trans[2, 4] - r4 * trans1[2, 4]
                      + dr2 * trans[3, 4] + r2 * trans1[3, 4])
        dtwiss[EPX] = (trans[1, 4] * dsqr + trans1[1, 4] * sqrdet + dr3 * trans[2, 4] + r3 * trans1[2, 4]
                       - dr1 * trans[3, 4] - r1 * trans1[3, 4])
        dtwiss[DETR] = dr1 * r4 + r1 * dr4 - dr2 * r3 - r2 * dr3
    else:
        dsqr = -0.5 * ddetp / sqrdet
        x11, x12, x21, x22 = trans[2, 0] / sqrdet, trans[2, 1] / sqrdet, trans[3, 0] / sqrdet, trans[3, 1] / sqrdet
        dx11 = (trans1[2, 0] - x11 * dsqr) / sqrdet
        dx12 = (trans1[2, 1] - x12 * dsqr) / sqrdet
        dx21 = (trans1[3, 0] - x21 * dsqr) / sqrdet
        dx22 = (trans1[3, 1] - x22 * dsqr) / sqrdet
        dr1 = (-trans1[0, 0] * x22 - trans[0, 0] * dx22 + trans1[0, 1] * x21 + trans[0, 1] * dx21)
        dr2 = (trans1[0, 0] * x12 + trans[0, 0] * dx12 - trans1[0, 1] * x11 - trans[0, 1] * dx11)
        dr3 = (-trans1[1, 0] * x22 - trans[1, 0] * dx22 + trans1[1, 1] * x21 + trans[1, 1] * dx21)
        dr4 = (trans1[1, 0] * x12 + trans[1, 0] * dx12 - trans1[1, 1] * x11 - trans[1, 1] * dx11)
        dtwiss[EX] = (trans[2, 4] * dsqr + trans1[2, 4] * sqrdet
                      + (-dr4 * trans[0, 4] - r4 * trans1[0, 4] + dr2 * trans[1, 4] + r2 * trans1[1, 4]))
        dtwiss[EPX] = (trans[3, 4] * dsqr + trans1[3, 4] * sqrdet
                       + (dr3 * trans[0, 4] + r3 * trans1[0, 4] - dr1 * trans[1, 4] - r1 * trans1[1, 4]))
        dtwiss[DETR] = -dr1 * r4 - r1 * dr4 + dr2 * r3 + r2 * dr3
    dtwiss[R1], dtwiss[R2], dtwiss[R3], dtwiss[R4] = dr1, dr2, dr3, dr4
    dtwiss[PEX] = trans1[0, 4]
    dtwiss[PEPX] = trans1[1, 4]
    ax0, bx0 = u[k, AX], u[k, BX]
    gx0 = (1.0 + ax0**2) / bx0
    dax = ((dx11 * x22 + x11 * dx22 + dx12 * x21 + x12 * dx21) * ax0
           - (dx11 * x21 + x11 * dx21) * bx0 - (dx12 * x22 + x12 * dx22) * gx0)
    dtwiss[AX] = dax / (1.0 + u[l, AX]**2)
    dbx = 2.0 * (-(dx11 * x12 + x11 * dx12) * ax0 + dx11 * x11 * bx0 + dx12 * x12 * gx0)
    dtwiss[BX] = dbx / u[l, BX]
    dtwiss[NX] = ((dx12 * (bx0 * x11 - ax0 * x12) - x12 * (bx0 * dx11 - ax0 * dx12))
                  / (x12**2 + (bx0 * x11 - ax0 * x12)**2))
    dtwiss[DX] = dcod1[0]
    dtwiss[DPX] = dcod1[1]
    if normal:
        y11, y12, y21, y22 = trans[2, 2] / sqrdet, trans[2, 3] / sqrdet, trans[3, 2] / sqrdet, trans[3, 3] / sqrdet
        dy11 = (trans1[2, 2] - y11 * dsqr) / sqrdet
        dy12 = (trans1[2, 3] - y12 * dsqr) / sqrdet
        dy21 = (trans1[3, 2] - y21 * dsqr) / sqrdet
        dy22 = (trans1[3, 3] - y22 * dsqr) / sqrdet
        dtwiss[EY] = (trans[2, 4] * dsqr + trans1[2, 4] * sqrdet + dr1 * trans[0, 4] + r1 * trans1[0, 4]
                      + dr2 * trans[1, 4] + r2 * trans1[1, 4])
        dtwiss[EPY] = (trans[3, 4] * dsqr + trans1[3, 4] * sqrdet + dr3 * trans[0, 4] + r3 * trans1[0, 4]
                       + dr4 * trans[1, 4] + r4 * trans1[1, 4])
    else:
        y11, y12, y21, y22 = trans[0, 2] / sqrdet, trans[0, 3] / sqrdet, trans[1, 2] / sqrdet, trans[1, 3] / sqrdet
        dy11 = (trans1[0, 2] - y11 * dsqr) / sqrdet
        dy12 = (trans1[0, 3] - y12 * dsqr) / sqrdet
        dy21 = (trans1[1, 2] - y21 * dsqr) / sqrdet
        dy22 = (trans1[1, 3] - y22 * dsqr) / sqrdet
        dtwiss[EY] = (trans[0, 4] * dsqr + trans1[0, 4] * sqrdet
                      + (dr1 * trans[2, 4] + r1 * trans1[2, 4] + dr2 * trans[3, 4] + r2 * trans1[3, 4]))
        dtwiss[EPY] = (trans[1, 4] * dsqr + trans1[1, 4] * sqrdet
                       + (dr3 * trans[2, 4] + r3 * trans1[2, 4] + dr4 * trans[3, 4] + r4 * trans1[3, 4]))
    dtwiss[PEY] = trans1[2, 4]
    dtwiss[PEPY] = trans1[3, 4]
    ay0, by0 = u[k, AY], u[k, BY]
    gy0 = (1.0 + ay0**2) / by0
    day = ((dy11 * y22 + y11 * dy22 + dy12 * y21 + y12 * dy21) * ay0
           - (dy11 * y21 + y11 * dy21) * by0 - (dy12 * y22 + y12 * dy22) * gy0)
    dtwiss[AY] = day / (1.0 + u[l, AY]**2)
    dby = 2.0 * (-(dy11 * y12 + y11 * dy12) * ay0 + dy11 * y11 * by0 + dy12 * y12 * gy0)
    dtwiss[BY] = dby / u[l, BY]
    dtwiss[NY] = ((dy12 * (by0 * y11 - ay0 * y12) - y12 * (by0 * dy11 - ay0 * dy12))
                  / (y12**2 + (by0 * y11 - ay0 * y12)**2))
    dtwiss[DY] = dcod1[2]
    dtwiss[DPY] = dcod1[3]
    dtwiss[DZ] = dcod1[4]
    az0 = u[k, AZ]
    dtwiss[GMX] = 2.0 * dtwiss[AX] * ax0 / (1.0 + ax0**2) - dtwiss[BX]
    dtwiss[GMY] = 2.0 * dtwiss[AY] * ay0 / (1.0 + ay0**2) - dtwiss[BY]
    dtwiss[GMZ] = 2.0 * dtwiss[AZ] * az0 / (1.0 + az0**2) - dtwiss[BZ]
    if l == ctx.nlat:
        bxx = math.sqrt(u[l, BX] / u[1, BX])
        cosmx, sinmx = math.cos(u[l, NX]), math.sin(u[l, NX])
        dtwiss[TRX] = (bxx * (0.5 * (cosmx + u[1, AX] * sinmx) * dtwiss[BX]
                              + (-sinmx + u[1, AX] * cosmx) * dtwiss[NX])
                       + (-0.5 * (cosmx - u[l, AX] * sinmx) * dtwiss[BX]
                          + (-sinmx - u[l, AX] * cosmx) * dtwiss[NX] - sinmx * dax) / bxx)
        byy = math.sqrt(u[l, BY] / u[1, BY])
        cosmy, sinmy = math.cos(u[l, NY]), math.sin(u[l, NY])
        dtwiss[TRY] = (byy * (0.5 * (cosmy + u[1, AY] * sinmy) * dtwiss[BY]
                              + (-sinmy + u[1, AY] * cosmy) * dtwiss[NY])
                       + (-0.5 * (cosmy - u[l, AY] * sinmy) * dtwiss[BY]
                          + (-sinmy - u[l, AY] * cosmy) * dtwiss[NY] - sinmy * day) / byy)
    return _resetnan(dtwiss)


# ---------------------------------------------------------------------------
# qdcell.f qdcell (周期解の補正)
# ---------------------------------------------------------------------------
def qdcell(ctx: DerivContext, kk: int, ll: int, iv: int) -> np.ndarray:
    """要素 kk の変数 iv による位置 ll の Twiss 微分 dtwiss (CELL なら周期解として)。
    流れ: qdtwis で「kk の入口を固定した輸送」の微分 dtwiss を得る → それで ll の一周行列の微分 dM を作る
    (dM = d(M_{kk→ll})·M_{ll→kk}) → R の微分 (qdmdia) と正規モードブロック A, B の微分 da, db →
    tr A = 2cos μ, a12 = β sin μ, a11−a22 = 2α sin μ から dμ, dβ/β, dα → 位相 ψ_ll の微分 (dμ と ll→終端 の
    位相進みの変化) → 閉軌道・分散の周期解の微分 (I−M)⁻¹"""
    u = ctx.u
    nut = ctx.nut
    k, l = kk, ll
    if ctx.cell:
        dpsix = u[nut, NX] - u[1, NX]
        dpsiy = u[nut, NY] - u[1, NY]
        thx = math.tan(0.5 * dpsix)
        sinmux = 2.0 * thx / (1.0 + thx**2)
        dcosx = thx * sinmux
        cosmux = 1.0 - dcosx
        detimx = 2.0 * dcosx
        thy = math.tan(0.5 * dpsiy)
        sinmuy = 2.0 * thy / (1.0 + thy**2)
        dcosy = thy * sinmuy
        cosmuy = 1.0 - dcosy
        detimy = 2.0 * dcosy
        # (trx/try の安定判定は qddtwiss の分岐にしか効かない → 省略)
    else:
        sinmux, cosmux, sinmuy, cosmuy = 0.0, 1.0, 0.0, 1.0
        dpsix = dpsiy = 0.0
        detimx = detimy = 0.0
    dtwiss = qdtwis(ctx, kk, ll, iv)
    if not ctx.cell:
        return _finish_gamma(dtwiss, u, l)
    # ====== 一周行列 tm (位置 l、物理座標) ======
    axl, bxl, ayl, byl = u[l, AX], u[l, BX], u[l, AY], u[l, BY]
    x11 = cosmux + axl * sinmux
    x22 = cosmux - axl * sinmux
    x12 = bxl * sinmux
    x21 = -(1.0 + axl**2) * sinmux / bxl
    y11 = cosmuy + ayl * sinmuy
    y22 = cosmuy - ayl * sinmuy
    y12 = byl * sinmuy
    y21 = -(1.0 + ayl**2) * sinmuy / byl
    dalphx = (1.0 + axl**2) * dtwiss[AX]
    dbetax = 0.5 * dtwiss[BX]
    dalphy = (1.0 + ayl**2) * dtwiss[AY]
    dbetay = 0.5 * dtwiss[BY]
    r11, r12, r21, r22 = u[l, R1], u[l, R2], u[l, R3], u[l, R4]
    detr = r11 * r22 - r12 * r21
    c1 = math.sqrt(1.0 - detr)
    # s: 物理 → 正規モード、sr: 正規モード → 物理 (qtentu と同じ形)
    s = np.array([[c1, 0.0, -r22, r12],
                  [0.0, c1, r21, -r11],
                  [r11, r12, c1, 0.0],
                  [r21, r22, 0.0, c1]])
    sr = np.array([[c1, 0.0, r22, -r12],
                   [0.0, c1, -r21, r11],
                   [-r11, -r12, c1, 0.0],
                   [-r21, -r22, 0.0, c1]])
    xm = np.array([[x11, x12], [x21, x22]])
    ym = np.array([[y11, y12], [y21, y22]])
    dm = np.zeros((4, 4))
    dm[0:2, 0:2] = xm
    dm[2:4, 2:4] = ym
    tm = sr @ dm @ s                                     # M = S⁻¹·Md·S
    # ====== dM = dS⁻¹·Md·S + S⁻¹·dMd·S ======
    dx11 = dbetax * x11 + (-sinmux + axl * cosmux) * dtwiss[NX]
    dx12 = dbetax * x12 + bxl * cosmux * dtwiss[NX]
    dx21 = -x21 * dbetax + (-(axl * sinmux + cosmux) * dalphx - (1.0 + axl**2) * cosmux * dtwiss[NX]) / bxl
    dx22 = -dbetax * x22 - sinmux * dalphx - (sinmux + axl * cosmux) * dtwiss[NX]
    dy11 = dbetay * y11 + (-sinmuy + ayl * cosmuy) * dtwiss[NY]
    dy12 = dbetay * y12 + byl * cosmuy * dtwiss[NY]
    dy21 = -y21 * dbetay + (-(ayl * sinmuy + cosmuy) * dalphy - (1.0 + ayl**2) * cosmuy * dtwiss[NY]) / byl
    dy22 = -dbetay * y22 - sinmuy * dalphy - (sinmuy + ayl * cosmuy) * dtwiss[NY]
    dr11, dr12, dr21, dr22 = dtwiss[R1], dtwiss[R2], dtwiss[R3], dtwiss[R4]
    ddetr = dr11 * r22 + dr22 * r11 - dr12 * r21 - dr21 * r12
    dc1 = -0.5 * ddetr / c1
    dsr = np.array([[dc1, 0.0, dr22, -dr12],
                    [0.0, dc1, -dr21, dr11],
                    [-dr11, -dr12, dc1, 0.0],
                    [-dr21, -dr22, 0.0, dc1]])
    ddm = np.zeros((4, 4))
    ddm[0:2, 0:2] = [[dx11, dx12], [dx21, dx22]]
    ddm[2:4, 2:4] = [[dy11, dy12], [dy21, dy22]]
    dtm = dsr @ dm @ s + sr @ ddm @ s
    # ====== R の微分 (周期境界条件) ======
    dr1, dr2, dr3, dr4 = qdmdia(tm, dtm, r11)
    dtwiss[R1], dtwiss[R2], dtwiss[R3], dtwiss[R4] = dr1, dr2, dr3, dr4
    ddetr0 = dr1 * r22 + r11 * dr4 - dr2 * r21 - r12 * dr3
    normal = u[l, DETR] < 1.0
    dtwiss[DETR] = ddetr0 if normal else -ddetr0
    dc = -0.5 * ddetr0 / c1
    ds = np.array([[dc, 0.0, -dr4, dr2],
                   [0.0, dc, dr3, -dr1],
                   [dr1, dr2, dc, 0.0],
                   [dr3, dr4, 0.0, dc]])
    dsr2 = np.array([[dc, 0.0, dr4, -dr2],
                     [0.0, dc, -dr3, dr1],
                     [-dr1, -dr2, dc, 0.0],
                     [-dr3, -dr4, 0.0, dc]])
    # ====== 対角化された一周行列の微分 da, db  (d(S·M·S⁻¹) の対角ブロック) ======
    ts = tm @ sr
    dts = dtm @ sr
    tds = tm @ dsr2
    dmd = ds @ ts + s @ dts + s @ tds
    da11, da12, da21, da22 = dmd[0, 0], dmd[0, 1], dmd[1, 0], dmd[1, 1]
    db11, db12, db21, db22 = dmd[2, 2], dmd[2, 3], dmd[3, 2], dmd[3, 3]
    # ====== Twiss の微分 ======
    if sinmux != 0.0:
        dmux = -0.5 * (da11 + da22) / sinmux
        dtwiss[AX] = (0.5 * (da11 - da22) - axl * cosmux * dmux) / sinmux
        dtwiss[BX] = (da12 - bxl * cosmux * dmux) / x12
    else:
        dmux = 0.0
    if k < l:
        sinl = math.sin(dpsix - u[l, NX])
        cosl = math.cos(dpsix - u[l, NX])
        dtwiss[NX] = dmux + sinl * ((cosl + axl * sinl) * dtwiss[BX] - sinl * dtwiss[AX])
    else:
        sinl = math.sin(u[l, NX])
        cosl = math.cos(u[l, NX])
        dtwiss[NX] = -sinl * ((cosl - axl * sinl) * dtwiss[BX] + sinl * dtwiss[AX])
    dtwiss[AX] = dtwiss[AX] / (1.0 + axl**2)
    if sinmuy != 0.0:
        dmuy = -0.5 * (db11 + db22) / sinmuy
        dtwiss[AY] = (0.5 * (db11 - db22) - ayl * cosmuy * dmuy) / sinmuy
        dtwiss[BY] = (db12 - byl * cosmuy * dmuy) / y12
    else:
        dmuy = 0.0
    if k < l:
        sinl = math.sin(dpsiy - u[l, NY])
        cosl = math.cos(dpsiy - u[l, NY])
        dtwiss[NY] = dmuy + sinl * ((cosl + ayl * sinl) * dtwiss[BY] - sinl * dtwiss[AY])
    else:
        sinl = math.sin(u[l, NY])
        cosl = math.cos(u[l, NY])
        dtwiss[NY] = -sinl * ((cosl - ayl * sinl) * dtwiss[BY] + sinl * dtwiss[AY])
    dtwiss[AY] = dtwiss[AY] / (1.0 + ayl**2)
    # ====== 閉軌道の微分: 2×2 の正規モードで dx = (I−M)⁻¹·(dM·x + δx) ======
    codv = np.array([u[l, DX], u[l, DPX], u[l, DY], u[l, DPY]])
    dcodv = np.array([dtwiss[DX], dtwiss[DPX], dtwiss[DY], dtwiss[DPY]])
    ddcod = s @ dcodv
    if normal:
        ddcod1, ddcod2, ddcod3, ddcod4 = ddcod
    else:
        ddcod3, ddcod4, ddcod1, ddcod2 = ddcod
    if detimx != 0.0:
        dcod1 = ((1.0 - x22) * ddcod1 + x12 * ddcod2) / detimx
        dcod2 = (x21 * ddcod1 + (1.0 - x11) * ddcod2) / detimx
    else:
        dcod1 = dcod2 = 0.0
    if detimy != 0.0:
        dcod3 = ((1.0 - y22) * ddcod3 + y12 * ddcod4) / detimy
        dcod4 = (y21 * ddcod3 + (1.0 - y11) * ddcod4) / detimy
    else:
        dcod3 = dcod4 = 0.0
    back = sr @ np.array([dcod1, dcod2, dcod3, dcod4])
    if normal:
        dtwiss[DX], dtwiss[DPX], dtwiss[DY], dtwiss[DPY] = back
    else:
        dtwiss[DY], dtwiss[DPY], dtwiss[DX], dtwiss[DPX] = back
    # ====== 分散の微分 ======
    dpe = s @ np.array([dtwiss[PEX], dtwiss[PEPX], dtwiss[PEY], dtwiss[PEPY]])
    if normal:
        da13, da23, db13, db23 = dpe
    else:
        db13, db23, da13, da23 = dpe
    if detimx != 0.0:
        dea = (da13 + x12 * da23 - x22 * da13) / detimx
        depa = (da23 + x21 * da13 - x11 * da23) / detimx
    else:
        dea = depa = 0.0
    if detimy != 0.0:
        deb = (db13 + y12 * db23 - y22 * db13) / detimy
        depb = (db23 + y21 * db13 - y11 * db23) / detimy
    else:
        deb = depb = 0.0
    ex, epx, ey, epy = u[l, EX], u[l, EPX], u[l, EY], u[l, EPY]
    pe = sr @ np.array([ex, epx, ey, epy])                       # 物理分散 pex..
    de = dsr2 @ pe                                                # ds·pe (R の変化による正規モード分散の変化)
    dtwiss[EX] = dea + de[0]
    dtwiss[EPX] = depa + de[1]
    dtwiss[EY] = deb + de[2]
    dtwiss[EPY] = depb + de[3]
    dpev = sr @ np.array([dea, depa, deb, depb])
    if normal:
        dtwiss[PEX], dtwiss[PEPX], dtwiss[PEY], dtwiss[PEPY] = dpev
    else:
        dtwiss[PEY], dtwiss[PEPY], dtwiss[PEX], dtwiss[PEPX] = dpev
    return _finish_gamma(dtwiss, u, l)


def _finish_gamma(dtwiss: np.ndarray, u: np.ndarray, l: int) -> np.ndarray:
    """qdcell 末尾の dgammag。BMAG の参照 Twiss (twiss(ll,-1,·)) は持たないので同じ点を渡す (BMAG は未対応)"""
    dtwiss[GMX], dtwiss[BMAGX] = dgammag(dtwiss[AX], dtwiss[BX], u[l, AX], u[l, BX], u[l, AX], u[l, BX])
    dtwiss[GMY], dtwiss[BMAGY] = dgammag(dtwiss[AY], dtwiss[BY], u[l, AY], u[l, BY], u[l, AY], u[l, BY])
    return dtwiss


# ---------------------------------------------------------------------------
# qdcell.f tffsqu (qu の組み立て)
# ---------------------------------------------------------------------------
def _variable_iv(v: dict) -> Optional[int]:
    return _IV_TABLE.get((v['type'], v['key']))


def _positions(m: Matcher, v: dict) -> List[int]:
    """変数 v の要素インスタンスのラティス上の位置 k (1 始まり) の一覧 (同じオブジェクトが複数回現れても全部)"""
    ids = {id(e) for e in v['els']}
    pos = [k for k, el in enumerate(m.lat.elements, start=1) if id(el) in ids]
    if not pos:                                    # 念のため名前で
        pos = [k for k, el in enumerate(m.lat.elements, start=1) if el.name == v['family']]
    return pos


def _row_wiq(m: Matcher, cid: int, idp: int, rows) -> float:
    """tffscalc の wiq(i): wfit (twfit) × FitWeight フック × オフモメンタム重み wi^(1/wexponent)"""
    if cid < 0:
        return 1.0
    c = m._cond(cid)
    wiq = float(c.get('wfit', 1.0))
    if m.fit_weight is not None:
        dps = m._dp_list(m._nfr())
        w2 = m.fit_weight(m._elname(c['kp']), c['kf'], m.dp0 + dps.get(idp, 0.0), wiq)
        if isinstance(w2, (int, float)):
            wiq = float(w2)
    if idp != 0:
        wi = (m.offmw / 2.0 / math.sqrt(max(1, abs(c['mfitp']))))**2
        wiq = wiq * wi**(1.0 / m.wexponent)
    return wiq


def tffsqu(m: Matcher, ids: List[tuple], rows_by: Dict[int, list]) -> Tuple[np.ndarray, List[int]]:
    """tffsqu の直訳。ids は match6._tffscalc が返す dfs のキー [(cid, idp), ...] (行の順序)。
    rows_by[idp] は各オフモメンタム点の twiss_table。戻り: (qu0 [条件, 変数] (wiq, wvar を掛ける前), 数値微分が必要な列)。
    SAD: rlist(kqu) = Σ_k coup(k) · (dtwiss(kf) at lfp(1)) − Σ_k coup(k) · (dtwiss(kf) at lfp(2))"""
    nvar = len(m.variables)
    qu = np.zeros((len(ids), nvar))
    numeric_cols: List[int] = []
    ctxs: Dict[int, DerivContext] = {}
    for idp, rws in rows_by.items():
        if rws is not None:
            ctxs[idp] = DerivContext(m.lat.elements, build_utwiss(rws), m.cell)
    # 条件行の情報
    rowinfo = []
    for i, (cid, idp) in enumerate(ids):
        if cid < 0:
            rowinfo.append(None)                  # FitFunction の行: 数値微分 (tffsmatch の EVDeriv)
            continue
        c = m._cond(cid)
        kf = _KF_INDEX.get(c['kf'])
        kpb, kpe = min(c['kp'], c['kp1']), max(c['kp'], c['kp1'])
        lfp1 = kpe
        lfp2 = kpb if c['kp'] != c['kp1'] else 0
        rowinfo.append((kf, lfp1, lfp2, idp))
    for ii, v in enumerate(m.variables):
        iv = _variable_iv(v)
        if iv is None:
            numeric_cols.append(ii)
            continue
        for k in _positions(m, v):
            coup = 1.0                                     # couple(k)*wk (COUPLE 未使用、境界の分数なし)
            cache: Dict[tuple, np.ndarray] = {}
            for i, info in enumerate(rowinfo):
                if info is None:
                    continue
                kf, lfp1, lfp2, idp = info
                ctx = ctxs.get(idp)
                if ctx is None or kf is None:
                    continue
                if kf == LENG:
                    if v['key'] == 'L':
                        posk = m._pos(rows_by[idp], k)
                        for lp, sgn in ((lfp1, 1.0), (lfp2, -1.0)):
                            if lp <= 0:
                                continue
                            if posk < m._pos(rows_by[idp], lp) and posk >= 0.0:
                                qu[i, ii] += sgn * coup
                            elif posk >= m._pos(rows_by[idp], lp) and posk < 0.0:
                                qu[i, ii] -= sgn * coup
                    continue
                if kf > MFITTRY:
                    continue                           # BMAG / 幾何 (tdgeo) は未対応
                for lp, sgn in ((lfp1, 1.0), (lfp2, -1.0)):
                    if lp <= 0:
                        continue
                    if not (k < lp or m.cell or m._pos(rows_by[idp], k) < m._pos(rows_by[idp], lp)):
                        continue
                    key = (lp, idp)
                    if key not in cache:
                        cache[key] = qdcell(ctx, k, lp, iv)
                    qu[i, ii] += sgn * coup * cache[key][kf]
    return qu, numeric_cols


def analytic_qu(m: Matcher, ids: List[tuple], rows, wvar: np.ndarray,
                rows_by: Optional[Dict[int, list]] = None) -> Tuple[np.ndarray, List[int]]:
    """match6.Matcher.go の数値微分ブロックの代替: qu[i, kc] = tffsqu·wiq[i]/wvar[kc] (tffsmatch.f 400-403 行目)。
    rows は idp=0 の twiss_table (直前の _tffscalc の戻り)。オフモメンタム条件があれば rows_by を渡す
    (無ければ m._optics_all で計算する)。戻り (qu, 数値微分で埋めるべき列の番号)。FitFunction 行 (cid<0) は 0 のまま"""
    if rows_by is None:
        rows_by = {0: rows}
        if any(idp != 0 for (_, idp) in ids):
            rb, _, _ = m._optics_all(m._nfr())
            rows_by = dict(rb)
            rows_by[0] = rows
    qu, numeric_cols = tffsqu(m, ids, rows_by)
    for i, (cid, idp) in enumerate(ids):
        if cid < 0:
            continue
        qu[i, :] *= _row_wiq(m, cid, idp, rows)
    qu /= wvar[None, :]
    return qu, numeric_cols


# ---------------------------------------------------------------------------
# match6.Matcher の go() を解析微分に置き換えたサブクラス
# ---------------------------------------------------------------------------
class AnalyticMatcher(Matcher):
    """互換のための別名。解析微分は match6.Matcher の既定になった (Matcher._derivative が analytic_qu を呼ぶ)。
    Twiss 条件は tffsqu の解析微分、FitFunction 行は片側差分、未対応の変数列・条件は中心差分。
    numerical_derivative=True で中心差分のみ (FFS$NumericalDerivative=True)。"""
    # go() / _derivative() / _numeric_columns() は match6.Matcher に統合済み (既定が解析微分)
    pass


def numeric_qu(m: Matcher, ids: List[tuple], dfs: dict, wvar: np.ndarray, step_scale: float = 1.0) -> np.ndarray:
    """検証用: match6.Matcher.go と同じ中心差分で qu を作る (dfs は現在点の _tffscalc の戻り)。
    step_scale で刻み dvkc を何倍かにできる (雑音判定用)"""
    nvar = len(m.variables)
    qu = np.zeros((len(ids), nvar))
    for kc, v in enumerate(m.variables):
        dvkc = max(abs(v['val']) * EPS, abs(EPS1 / wvar[kc])) * step_scale
        valbest = v['val']
        v['val'] = valbest + dvkc
        dfs1, _, _ = m._tffscalc(False)
        v['val'] = valbest - dvkc
        dfs2, _, _ = m._tffscalc(False)
        v['val'] = valbest
        for i, cid in enumerate(ids):
            ddf1 = dfs[cid] - dfs1.get(cid, dfs[cid])
            ddf2 = dfs[cid] - dfs2.get(cid, dfs[cid])
            qu[i, kc] = (ddf1 - ddf2) / 2.0 / dvkc / wvar[kc]
    m._apply_vars()
    return qu


def compare_derivatives(m: Matcher, verbose: bool = True, tol: float = 1e-6):
    """現在の変数・条件で解析微分と数値微分を比べる。戻り (qu_analytic, qu_numeric, ids, 最大相対差)。
    数値微分は刻み dvkc と 2·dvkc の 2 通りで計算し、両者が 1% 以上違う成分は数値微分の雑音 (閉軌道探索の
    収束誤差や対称点でゼロになる微分の丸め) とみなして比較から外す ("noise" と表示)。"""
    m.lat_brho = getattr(m.lat, 'brho', 1.0)
    dfs, r, rows = m._tffscalc(True)
    ids = list(dfs.keys())
    wvar = np.array([m._tweigh(v, v['val']) for v in m.variables])
    qa, ncols = analytic_qu(m, ids, rows, wvar)
    qn = numeric_qu(m, ids, dfs, wvar)
    qn2 = numeric_qu(m, ids, dfs, wvar, 2.0)
    names = [f"{v['family']}.{v['key']}" for v in m.variables]
    worst = 0.0
    nnoise = 0
    for i, (cid, idp) in enumerate(ids):
        c = m._cond(cid)
        tag = f"{c['kf']}@{m._elname(c['kp'])}" + (f"-{m._elname(c['kp1'])}" if c['kp'] != c['kp1'] else '')
        for kc in range(len(names)):
            a, n, n2 = qa[i, kc], qn[i, kc], qn2[i, kc]
            noise = abs(n - n2) > 0.01 * max(abs(n), abs(n2)) or (n == 0.0 and a == 0.0)
            rel = abs(a - n) / max(abs(n), 1e-300)
            if noise:
                flag = '  (noise)'
                nnoise += 1
            else:
                worst = max(worst, rel)
                flag = '' if rel < tol else '  <-- 不一致'
            if verbose:
                print(f'{tag:16s} {names[kc]:14s} analytic {a: .10e} numeric {n: .10e} rel {rel: .2e}{flag}')
    if verbose:
        print(f'最大相対差 {worst:.3e} (雑音として除外 {nnoise}/{qa.size} 成分)')
    return qa, qn, ids, worst
