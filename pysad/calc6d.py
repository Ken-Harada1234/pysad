"""
calc6d.py — SAD の FFS フラグ CALC6D (6 次元 CALC) の直訳
========================================================

SAD での呼び出し経路 (CELL; CALC6D; CALC;):
  tffscalc (tffscalc.f) → qcell/qcell1 (qcell.f) → `if(calc6d) call qcell6d`
  qcell6d (qcell.f:404):
    codplt=.true. にして temit (temit.f) を plot=.true. で呼ぶ
    temit:
      1. tcod (tcod.f)        6D 閉軌道 (RFSW ON なら RFSW OFF で求めてから z を同期位相に置き直す)
      2. tturne (plot=.false.) 一周写像 trans(6,12) = [T | D] (T: シンプレクティック部, D: 放射減衰部)
      3. tecalc (temit.f)      rx = T (RFSW OFF なら rx(6,1:5)=0, rx(6,6)=1) → teigenc/tnorm → r, ri=inv(r)
                               stab = (3 モードとも |Re log λ| < 1e-6) かつ閉軌道が見つかった
      4. plot=.true. でもう一周 tturne: 各要素入口 l で tsetetwiss (tturne.f:775)
             ti = ri · inv(T_{1→l});  twiss(l,:) = tfetwiss(ti, cod_l, normali)
             位相: l=1 は 0。l>1 は twi(nx) (<-toln なら +2π)。最後に tfadjustn で NX, NY を単調化 (NZ はしない)
         要素 nlat ($$$) は tfsetplot → tsetetwiss(一周の T, 一周後の cod)
    不安定 (stab=.false.) なら '*****qcell6d---> Unstable optics' を出し cell=.false. にして
      ri0 = etwiss2ri(直前の twiss(1)) から transport 線として伝搬 (tturne0, plot=.true.)

  拡張 Twiss 28 成分の定義は tfetwiss (tffs.f:2413)、逆変換は etwiss2ri (tffs.f:2569)。
  要素内の分数位置 (Twiss[key, 200.5]) は qtwissfrac (qtwiss.f:1008 の calc6d 分岐):
      trans = inv(etwiss2ri(twiss(l))) を部分要素で伝搬し ftwiss = tfetwiss(inv(trans), cod) に twiss(l) の位相を足す

RAD フラグ: temit の経路 (tturne/tbende/tradke) は ffs_flag の rad を参照しない (tbende の enarad は
  要素の RAD キーワード)。SAD で 18a (RAD) と 18c (NORAD) の出力が完全一致することを確認した。
  → calc6d は RAD フラグに依らず放射減衰込みで計算する (radiation=False は実験用)。

pysad での対応:
  SAD                      | ここ
  -------------------------|---------------------------------------------
  qcell6d                  | calc6d()
  tcod                     | _closed_orbit() (emit6d.emittance6d と同じ手順: tcod の直訳ではない)
  tturne/tturne1 (要素列)  | _turn()  (emit6d.tturne と同じ要素分岐 + 入口ごとの T, cod の記録)
  tecalc (固有値, r, ri)   | emit6d._normal_form / _tsymp / _tinv6 (+ emit6d.tecalc で Emittance 相当も返す)
  tsetetwiss               | _rows_from_turn()
  tfadjustn                | _adjust_phase()
  tfetwiss                 | tfetwiss()
  etwiss2ri                | etwiss2ri()
  qtwissfrac (calc6d)      | twiss_frac()
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from . import lattice as _lat
from .lattice import TwissRow
from .elements import Element, Drift, Quadrupole, Bend
from . import emit6d as _em
from .emit6d import (_Phys, _State, _tdrife0, _tbende, _tquade, _tcave, _tmulte, _tthine,
                     _rf_effective, _trf0_update, _normal_form, _tsymp, _tinv6,
                     tecalc, ELMASS)

XYTH = 0.375      # tffs.f:210  xyth (x-y モードの入れ替え判定)
TOLN = 0.1        # tffs.f:2316 toln (位相の 2π 補正の許容)
PI2 = 2.0 * math.pi


# ---------------------------------------------------------------------------
# tfetwiss / etwiss2ri の直訳
# ---------------------------------------------------------------------------

def tfetwiss(r: np.ndarray, cod, normi: bool = True) -> np.ndarray:
    """tffs.f tfetwiss の直訳。
    r : 6×6 (物理座標 → 正規モード座標 の行列。temw の ri に相当)
    cod : 閉軌道 6 成分
    normi : etwiss2ri が返す normal (x-y モードが入れ替わっていない)
    返り値: 長さ 28 の配列 (SAD の twiss(l,idp,1:28)、index 0 = mfitax ... 27 = mfitzpy)"""
    r = np.asarray(r, dtype=float)
    az = math.sqrt(r[4, 4] * r[5, 5] - r[5, 4] * r[4, 5])
    uz11 = r[4, 4] / az
    uz12 = r[4, 5] / az
    uz21 = r[5, 4] / az
    uz22 = r[5, 5] / az
    hx11 = r[5, 1] * uz11 - r[4, 1] * uz21
    hx12 = r[5, 1] * uz12 - r[4, 1] * uz22
    hx21 = -r[5, 0] * uz11 + r[4, 0] * uz21
    hx22 = -r[5, 0] * uz12 + r[4, 0] * uz22
    hy11 = r[5, 3] * uz11 - r[4, 3] * uz21
    hy12 = r[5, 3] * uz12 - r[4, 3] * uz22
    hy21 = -r[5, 2] * uz11 + r[4, 2] * uz21
    hy22 = -r[5, 2] * uz12 + r[4, 2] * uz22
    f = 1.0 / (1.0 + az)
    ax = 1.0 - (hx11 * hx22 - hx21 * hx12) * f
    ay = 1.0 - (hy11 * hy22 - hy21 * hy12) * f
    hi = np.zeros((6, 6))
    # Fortran の hi(i,j) → hi[i-1, j-1]
    hi[1, 1] = ax
    hi[0, 1] = 0.0
    hi[3, 1] = (hx12 * hy21 - hx11 * hy22) * f
    hi[2, 1] = -(-hx12 * hy11 + hx11 * hy12) * f
    hi[5, 1] = -hx11
    hi[4, 1] = hx12
    hi[1, 0] = 0.0
    hi[0, 0] = ax
    hi[3, 0] = -(hx22 * hy21 - hx21 * hy22) * f
    hi[2, 0] = (-hx22 * hy11 + hx21 * hy12) * f
    hi[5, 0] = hx21
    hi[4, 0] = -hx22
    hi[1, 3] = hi[2, 0]
    hi[0, 3] = -hi[2, 1]
    hi[3, 3] = ay
    hi[2, 3] = 0.0
    hi[5, 3] = -hy11
    hi[4, 3] = hy12
    hi[1, 2] = -hi[3, 0]
    hi[0, 2] = hi[3, 1]
    hi[3, 2] = 0.0
    hi[2, 2] = ay
    hi[5, 2] = hy21
    hi[4, 2] = -hy22
    hi[1, 5] = hx22
    hi[0, 5] = hx21
    hi[3, 5] = hy22
    hi[2, 5] = hy21
    hi[5, 5] = az
    hi[4, 5] = 0.0
    hi[1, 4] = hx12
    hi[0, 4] = hx11
    hi[3, 4] = hy12
    hi[2, 4] = hy11
    hi[5, 4] = 0.0
    hi[4, 4] = az
    hi = r @ hi
    detm = (hi[0, 0] * hi[1, 1] - hi[1, 0] * hi[0, 1]
            + hi[2, 2] * hi[3, 3] - hi[3, 2] * hi[2, 3]) * 0.5
    normal = detm > XYTH
    if not normal:
        his = hi[0, 0:4].copy()
        hi[0, 0:4] = hi[2, 0:4]
        hi[2, 0:4] = his
        his = hi[1, 0:4].copy()
        hi[1, 0:4] = hi[3, 0:4]
        hi[3, 0:4] = his
        detm = (hi[0, 0] * hi[1, 1] - hi[1, 0] * hi[0, 1]
                + hi[2, 2] * hi[3, 3] - hi[3, 2] * hi[2, 3]) * 0.5
    axy = math.sqrt(detm)
    r11 = (hi[3, 3] * hi[2, 0] - hi[2, 3] * hi[3, 0]) / axy
    r12 = (hi[3, 3] * hi[2, 1] - hi[2, 3] * hi[3, 1]) / axy
    r21 = (-hi[3, 2] * hi[2, 0] + hi[2, 2] * hi[3, 0]) / axy
    r22 = (-hi[3, 2] * hi[2, 1] + hi[2, 2] * hi[3, 1]) / axy
    crx = math.hypot(hi[0, 1], hi[1, 1])
    cx = hi[1, 1] / crx
    sx = -hi[0, 1] / crx
    bx21 = (-sx * hi[0, 0] + cx * hi[1, 0]) / axy
    bx22 = crx / axy
    cry = math.hypot(hi[2, 3], hi[3, 3])
    cy = hi[3, 3] / cry
    sy = -hi[2, 3] / cry
    by21 = (-sy * hi[2, 2] + cy * hi[3, 2]) / axy
    by22 = cry / axy
    crz = math.hypot(uz12, uz22)
    cz = uz22 / crz
    sz = -uz12 / crz
    bz21 = -sz * uz11 + cz * uz21
    bz22 = crz
    tw = np.zeros(28)
    tw[10] = r11
    tw[11] = r12
    tw[12] = r21
    tw[13] = r22
    if not normi:
        normal = not normal
    if normal:
        tw[0] = bx21 * bx22
        tw[1] = bx22 ** 2
        tw[2] = math.atan2(sx, cx)
        tw[3] = by21 * by22
        tw[4] = by22 ** 2
        tw[5] = math.atan2(sy, cy)
        tw[6] = axy * hx12 - r22 * hy12 + r12 * hy22
        tw[7] = axy * hx22 + r21 * hy12 - r11 * hy22
        tw[8] = axy * hy12 + r11 * hx12 + r12 * hx22
        tw[9] = axy * hy22 + r21 * hx12 + r22 * hx22
        tw[14] = r11 * r22 - r12 * r21
        tw[24] = axy * hx11 - r22 * hy11 + r12 * hy21
        tw[25] = axy * hx21 + r21 * hy11 - r11 * hy21
        tw[26] = axy * hy11 + r11 * hx11 + r12 * hx21
        tw[27] = axy * hy21 + r21 * hx11 + r22 * hx21
    else:
        tw[3] = bx21 * bx22
        tw[4] = bx22 ** 2
        tw[5] = math.atan2(sx, cx)
        tw[0] = by21 * by22
        tw[1] = by22 ** 2
        tw[2] = math.atan2(sy, cy)
        tw[8] = axy * hx12 - r22 * hy12 + r12 * hy22
        tw[9] = axy * hx22 + r21 * hy12 - r11 * hy22
        tw[6] = axy * hy12 + r11 * hx12 + r12 * hx22
        tw[7] = axy * hy22 + r21 * hx12 + r22 * hx22
        tw[14] = 1.0 + XYTH - r11 * r22 + r12 * r21
        tw[26] = axy * hx11 - r22 * hy11 + r12 * hy21
        tw[27] = axy * hx21 + r21 * hy11 - r11 * hy21
        tw[24] = axy * hy11 + r11 * hx11 + r12 * hx21
        tw[25] = axy * hy21 + r21 * hx11 + r22 * hx21
    tw[15:21] = np.asarray(cod, dtype=float)
    tw[21] = bz21 * bz22
    tw[22] = bz22 ** 2
    tw[23] = math.atan2(sz, cz)
    return tw


def etwiss2ri(tw) -> tuple:
    """tffs.f etwiss2ri の直訳。拡張 Twiss 28 成分 → ri (物理 → 正規モード座標) と normal (detr < 1)。"""
    tw = np.asarray(tw, dtype=float)
    r1, r2, r3, r4, detr = tw[10], tw[11], tw[12], tw[13], tw[14]
    normal = detr < 1.0
    if normal:
        amu = math.sqrt(1.0 - detr)
        ex, epx, ey, epy = tw[6], tw[7], tw[8], tw[9]
        zx, zpx, zy, zpy = tw[24], tw[25], tw[26], tw[27]
        sqrbx = math.sqrt(tw[1])
        ax = tw[0]
        sqrby = math.sqrt(tw[4])
        ay = tw[3]
    else:
        amu = math.sqrt(1.0 - r1 * r4 + r2 * r3)
        ey, epy, ex, epx = tw[6], tw[7], tw[8], tw[9]
        zy, zpy, zx, zpx = tw[24], tw[25], tw[26], tw[27]
        sqrby = math.sqrt(tw[1])
        ay = tw[0]
        sqrbx = math.sqrt(tw[4])
        ax = tw[3]
    sqrbz = math.sqrt(tw[22])
    az = tw[21]
    hx11 = amu * zx - r2 * zpy + r4 * zy
    hx12 = amu * ex - r2 * epy + r4 * ey
    hx21 = amu * zpx + r1 * zpy - r3 * zy
    hx22 = amu * epx + r1 * epy - r3 * ey
    hy11 = amu * zy - r2 * zpx - r1 * zx
    hy12 = amu * ey - r2 * epx - r1 * ex
    hy21 = amu * zpy - r4 * zpx - r3 * zx
    hy22 = amu * epy - r4 * epx - r3 * ex
    dethx = hx11 * hx22 - hx21 * hx12
    dethy = hy11 * hy22 - hy21 * hy12
    azz = math.sqrt(1.0 - dethx - dethy)
    aa = 1.0 / (1.0 + azz)
    amux = 1.0 - dethx * aa
    amuy = 1.0 - dethy * aa
    h = np.zeros((4, 6))
    h[0, 0] = amux
    h[0, 1] = 0.0
    h[0, 2] = (hx12 * hy21 - hx11 * hy22) * aa
    h[0, 3] = (-hx12 * hy11 + hx11 * hy12) * aa
    h[0, 4] = -hx11
    h[0, 5] = -hx12
    h[1, 0] = 0.0
    h[1, 1] = amux
    h[1, 2] = (hx22 * hy21 - hx21 * hy22) * aa
    h[1, 3] = (-hx22 * hy11 + hx21 * hy12) * aa
    h[1, 4] = -hx21
    h[1, 5] = -hx22
    h[2, 0] = h[1, 3]
    h[2, 1] = -h[0, 3]
    h[2, 2] = amuy
    h[2, 3] = 0.0
    h[2, 4] = -hy11
    h[2, 5] = -hy12
    h[3, 0] = -h[1, 2]
    h[3, 1] = h[0, 2]
    h[3, 2] = 0.0
    h[3, 3] = amuy
    h[3, 4] = -hy21
    h[3, 5] = -hy22
    ria = np.zeros((6, 6))
    ria[4, :] = [hx22, -hx12, hy22, -hy12, azz, 0.0]
    ria[5, :] = [-hx21, hx11, -hy21, hy11, 0.0, azz]
    br = np.zeros((4, 4))
    br[0, :] = [amu / sqrbx, 0.0, -r4 / sqrbx, r2 / sqrbx]
    br[1, :] = [amu * ax / sqrbx, amu * sqrbx, r3 * sqrbx - r4 * ax / sqrbx,
                r2 * ax / sqrbx - r1 * sqrbx]
    br[2, :] = [r1 / sqrby, r2 / sqrby, amu / sqrby, 0.0]
    br[3, :] = [r1 * ay / sqrby + r3 * sqrby, r2 * ay / sqrby + r4 * sqrby,
                amu * ay / sqrby, amu * sqrby]
    ria[0:4, :] = br @ h
    ria[4, :] = ria[4, :] / sqrbz
    ria[5, :] = ria[5, :] * sqrbz + ria[4, :] * az
    return ria, normal


# ---------------------------------------------------------------------------
# TwissRow との相互変換
# ---------------------------------------------------------------------------

_TW_FIELDS = ('alfx', 'betx', 'mux', 'alfy', 'bety', 'muy', 'etx', 'etpx', 'ety', 'etpy',
              'r1', 'r2', 'r3', 'r4', 'detr', 'x_cod', 'px_cod', 'y_cod', 'py_cod', 'z_cod',
              'delta_cod', 'alfz', 'betz', 'muz', 'zx', 'zpx', 'zy', 'zpy')


def row_from_tw(tw, s: float, name: str) -> TwissRow:
    """長さ 28 の twiss 配列 → TwissRow (lattice.py の列構成)"""
    return TwissRow(s=s, name=name, **{k: float(v) for k, v in zip(_TW_FIELDS, tw)})


def tw_from_row(row: TwissRow) -> np.ndarray:
    """TwissRow → 長さ 28 の twiss 配列"""
    return np.array([getattr(row, k) for k in _TW_FIELDS], dtype=float)


# ---------------------------------------------------------------------------
# 一周追跡 (emit6d.tturne と同じ要素分岐 + 要素入口ごとの記録)
# ---------------------------------------------------------------------------

class _State6d(_State):
    """emit6d._State に放射の ON/OFF を付けたもの (SAD は CALC6D で RAD を見ないので既定は ON)"""

    def __init__(self, cod, ph, radiation=True):
        super().__init__(cod, ph)
        self.radiation = radiation

    def tradke(self, al, phir0):
        if self.radiation:
            super().tradke(al, phir0)
        else:
            self.setr0()        # tradke の末尾 (codr0/transr の更新) だけ


def _turn(elements: List[Element], cod0, ph: _Phys, trf0: float = 0.0,
          rfsw: bool = False, radcod: bool = False, radiation: bool = True,
          record: bool = False):
    """tturne1 相当 (emit6d.tturne の要素分岐をそのまま使う)。
    record=True なら各要素入口 l (l=1..nlat) で (T_{1→l}, cod_l) を記録して st.rec に入れる
    (SAD の tsetetwiss が呼ばれる時点。要素 nlat ($$$) は一周後の T と cod)。"""
    st = _State6d(cod0, ph, radiation)
    st.trf0 = trf0
    st.rfsw = rfsw
    st.radcod = radcod
    st.setr0()
    rec = [] if record else None
    for el in elements:
        if record:
            rec.append((st.T.copy(), st.cod.copy()))
        st.setr0()                       # tturne1: 要素ごとに tsetr0
        if isinstance(el, Drift):
            _tdrife0(st, el.L)
        elif isinstance(el, Bend):
            _tbende(st, el, enarad=True)
        elif isinstance(el, Quadrupole):
            _tquade(st, el, krad=True)
        elif type(el).__name__ in _lat._THIN_KORD:
            _tthine(st, el)
        elif type(el).__name__ == 'Cavity':
            _tcave(st, el)
        elif type(el).__name__ == 'Multipole':
            _tmulte(st, el, enarad=True)
        elif el.L == 0.0:
            continue          # MARK, MONI など
        else:
            raise NotImplementedError(f"calc6d: 要素 {type(el).__name__} は未対応")
    if record:
        rec.append((st.T.copy(), st.cod.copy()))
    st.rec = rec
    return st


# ---------------------------------------------------------------------------
# 閉軌道 (emit6d.emittance6d と同じ手順)
# ---------------------------------------------------------------------------

def _closed_orbit(els, ph, rfsw, radcod, radiation):
    """temit の tcod 部分。emit6d.emittance6d の閉軌道の決め方をそのまま使う (tcod の直訳ではない):
      RFSW OFF            : 4D 閉軌道 (δ 固定) で一周して z が閉じる δ を割線法で求める
      RFSW ON, RADCOD OFF : 同上で z=0、trf0 は 1 周目の U0 と実効電圧から決める (tturne update)
      RFSW ON, RADCOD ON  : trf0 を固定し放射込みの一周写像の 6D 固定点を Newton 法で解く (radiation=False の経路)
    返り値: (cod, trf0, found)

    radiation=True (Emittance[] と CALC6D の既定) は emit6d.closed_orbit6d (tcod の直訳・cod だけ追う z_slip・解析ヤコビアン) に任せる。
    EMIT コマンドと同じ cod / trf0 になる (RADCOD ON の Equilibrium dz が SAD と一致)。"""
    if radiation:
        from .emit6d import closed_orbit6d
        cod, trf0, found, _ = closed_orbit6d(els, ph, rfsw, radcod, True, norf_trf0=True)
        return cod, trf0, found
    def z_slip(dp):
        c, _, _ = _lat._qcod_periodic(els, dp)
        c = np.array(c, dtype=float)
        c[4] = 0.0
        c[5] = dp
        st0 = _turn(els, c, ph, radiation=radiation)
        return st0.cod[4] - c[4], c
    if rfsw and radcod:
        c, _, _ = _lat._qcod_periodic(els, 0.0)
        cod0 = np.array(c, dtype=float)
        cod0[4] = 0.0
        cod0[5] = 0.0
        st1 = _turn(els, cod0, ph, trf0=0.0, rfsw=True, radiation=radiation)
        wrf, vce = _rf_effective(st1)
        trf0 = _trf0_update(st1, wrf, vce) if wrf != 0.0 else 0.0
        cod = cod0.copy()
        found = False
        for _ in range(30):
            st = _turn(els, cod, ph, trf0=trf0, rfsw=True, radcod=True, radiation=radiation)
            res = st.cod - cod
            if np.max(np.abs(res)) < 1e-15:
                found = True
                break
            M = st.T + st.D
            cod = cod - np.linalg.solve(M - np.eye(6), res)
        if not found:
            found = np.max(np.abs(res)) < 1e-9
        return cod, trf0, found
    dp0 = 0.0
    z0, cod0 = z_slip(dp0)
    if z0 != 0.0:
        dp1 = 1e-7
        z1, cod1 = z_slip(dp1)
        for _ in range(30):      # 4 回では K0 の総和が 0 でないとき (δ が 1e-3 近く動く) 収束しない
            if z1 == z0:
                break
            dp2 = dp1 - z1 * (dp1 - dp0) / (z1 - z0)
            dp0, z0 = dp1, z1
            dp1 = dp2
            z1, cod1 = z_slip(dp1)
            if abs(z1) < 1e-15:
                break
        cod0 = cod1
    cod0 = np.array(cod0, dtype=float)
    cod0[4] = 0.0
    if not rfsw:
        return cod0, 0.0, True
    st1 = _turn(els, cod0, ph, trf0=0.0, rfsw=True, radiation=radiation)
    wrf, vce = _rf_effective(st1)
    trf0 = _trf0_update(st1, wrf, vce) if wrf != 0.0 else 0.0
    return cod0, trf0, True


# ---------------------------------------------------------------------------
# 結果
# ---------------------------------------------------------------------------

@dataclass
class Calc6DResult:
    """calc6d の結果。rows[l-1] = 要素 l の入口 (SAD の twiss(l,0,1:28))、rows[-1] = $$$"""
    rows: List[TwissRow]
    stable: bool
    closed_orbit: np.ndarray
    trf0: float
    r: Optional[np.ndarray] = None          # temw の r  (正規モード → 物理)
    ri: Optional[np.ndarray] = None         # temw の ri (物理 → 正規モード)
    tunes: Optional[np.ndarray] = None      # 固有値から (νx, νy, νz)
    emit: Optional[dict] = None             # emit6d.tecalc の結果 (Emittance[] 相当)
    elements: List[Element] = field(default_factory=list)
    ph: Optional[_Phys] = None
    rfsw: bool = False
    radcod: bool = False
    radiation: bool = True
    messages: List[str] = field(default_factory=list)


def _adjust_phase(rows: List[TwissRow], key: str):
    """tfadjustn: l=2..nlat について前の要素より toln 以上減っていたら 2π を足す (累積)"""
    phi0 = 0.0
    for l in range(1, len(rows)):
        v = getattr(rows[l], key) + phi0
        if v + TOLN < getattr(rows[l - 1], key):
            phi0 += PI2
            v += PI2
        setattr(rows[l], key, v)


def _rows_from_turn(st, ri, normali, elements, cod_start=None, stable=True, row0=None):
    """tsetetwiss (lorg=1) の直訳: 各要素入口で ti = ri·inv(T_{1→l}), twiss = tfetwiss(ti, cod)。
    位相は l=1 で 0、l>1 は twi (<-toln なら +2π) を twiss(1) の位相 (=0) に足す。NZ は補正なし。"""
    rows = []
    s = 0.0
    for l, (T, cod) in enumerate(st.rec):
        ti = ri @ _tinv6(T)
        tw = tfetwiss(ti, cod, normali) if (l > 0 or row0 is None) else np.array(row0, dtype=float)
        if l == 0:
            tw[2] = 0.0
            tw[5] = 0.0
            tw[23] = 0.0
        else:
            if tw[2] < -TOLN:
                tw[2] += PI2
            if tw[5] < -TOLN:
                tw[5] += PI2
        name = '*' if l == 0 else elements[l - 1].name
        rows.append(row_from_tw(tw, s, name))
        if l < len(elements):
            s += elements[l].L
    _adjust_phase(rows, 'mux')
    _adjust_phase(rows, 'muy')
    return rows


def calc6d(lattice, momentum_eV: float, rfsw: bool = False, radcod: bool = False,
           rad: bool = True, charge: float = 1.0, amass: float = ELMASS,
           twiss0: Optional[TwissRow] = None, radiation: Optional[bool] = None,
           out=None, with_emit: bool = True, fshift: float = 0.0) -> Calc6DResult:
    """
    SAD の CELL; CALC6D; CALC; (qcell6d) の直訳。

    Parameters
    ----------
    lattice     : pysad.lattice.Lattice (emit6d と同じ要素列)
    momentum_eV : 運動量 [eV/c]
    rfsw, radcod: FFS フラグ。rad は受け取るが SAD と同じく結果に影響しない
                  (radiation=False を明示したときだけ放射減衰を切る: 実験用)
    twiss0      : 不安定なときに使う入口の Twiss (SAD の「直前の twiss(1)」)。None なら βx=βy=βz=1
    out         : メッセージの出力先 (write を持つもの)。None なら print しない
    with_emit   : True なら emit6d.tecalc も呼んで Emittance[] 相当 (emit) を結果に入れる

    Returns
    -------
    Calc6DResult: rows (長さ nelem+1 の TwissRow、rows[0] は入口 '*'、rows[i] は要素 i の出口 = 要素 i+1 の入口)
    """
    els = list(lattice.elements)
    circ = sum(el.L for el in els)
    ph = _Phys(momentum_eV, circ, amass=amass, charge=charge)
    ph.dvfs = fshift if rfsw else 0.0          # tsetdvfs (FSHIFT)
    if radiation is None:
        radiation = True
    msgs = []

    def say(m):
        msgs.append(m)
        if out is not None:
            out.write(m + chr(10))

    rfsw0 = _lat.RFSW
    _lat.RFSW = rfsw
    try:
        cod, trf0, found = _closed_orbit(els, ph, rfsw, radcod, radiation)
        # temit: 一周して trans(6,12)、tecalc で r, ri、安定判定
        st = _turn(els, cod, ph, trf0=trf0, rfsw=rfsw, radcod=radcod,
                   radiation=radiation, record=True)
    finally:
        _lat.RFSW = rfsw0
    rx = st.T.copy()
    if not rfsw:
        rx[5, :5] = 0.0
        rx[5, 5] = 1.0
    stable = found
    r = ri = None
    tunes = None
    emit = None
    if stable:
        try:
            r, ceig = _normal_form(rx)
            r = _tsymp(r)
            ri = _tinv6(r)
            cd = np.log(ceig[0::2])
            stable = bool(np.all(np.abs(cd.real) < 1e-6)) and np.all(np.isfinite(r))
            tunes = cd.imag / PI2
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            stable = False
    if stable:
        rows = _rows_from_turn(st, ri, True, els)
        if with_emit:
            try:
                emit = tecalc(st, circ, rfsw=rfsw)
                emit['closed_orbit'] = cod.copy()
            except Exception:       # 平衡ビーム行列が解けなくても Twiss は返す
                emit = None
    else:
        # qcell6d: '*****qcell6d---> Unstable optics' → cell=.false. で transport 線として伝搬
        say(' *****qcell6d---> Unstable optics')
        if twiss0 is None:
            twiss0 = TwissRow(betx=1.0, bety=1.0, betz=1.0)
        tw0 = tw_from_row(twiss0)
        if tw0[22] <= 0.0:
            tw0[22] = 1.0
        ri0, normal0 = etwiss2ri(tw0)
        cod0 = tw0[15:21].copy()
        st = _turn(els, cod0, ph, trf0=trf0, rfsw=rfsw, radcod=radcod,
                   radiation=radiation, record=True)
        rows = _rows_from_turn(st, ri0, normal0, els)
        # twiss(1) = tfetwiss(ri0, cod0, .true.) の位相 (qcell.f:446) に合わせて先頭の位相は 0 のまま
        cod = cod0
    return Calc6DResult(rows=rows, stable=stable, closed_orbit=np.array(cod, dtype=float),
                        trf0=trf0, r=r, ri=ri, tunes=tunes, emit=emit, elements=els, ph=ph,
                        rfsw=rfsw, radcod=radcod, radiation=radiation, messages=msgs)


def calc6d_ins(lattice, momentum_eV: float, tw0, rfsw: bool = False, radcod: bool = False,
               trf0: float = 0.0, rad: bool = True, charge: float = 1.0, amass: float = ELMASS,
               radiation: Optional[bool] = None, fshift: float = 0.0) -> Calc6DResult:
    """
    SAD の INS; CALC6D; CALC; (qcell.f:404 qcell6d の cell=.false. の枝) の直訳。docs/history/work_emit_unstable.md B 節。

      qcell: `if(calc6d) call qcell6d` — CELL か INS かに関わらず 6D 経路に入る (qcell.f:62)。INS (cell=.false.) では
        ri0  = etwiss2ri(twiss(1,idp,1:28), normali)         先頭 MARK の 28 成分 (AX..ZPY: tffscalc の twmov(right) で MARK から読む)
        cod0 = twiss(1,idp,mfitdx:mfitddp)                   (DX DPX DY DPY DZ DDP)
        twiss(1) = tfetwiss(ri0, cod0, .true.)               (28 成分を ri0 から作り直す)
        beam = tfinibeam(1)  → temw の ri = etwiss2ri(twiss(1))  (tsetetwiss の lorg=1 が使う ri; normal は局所)
        tinitr(trans); tturne0(trans, cod0, beam, srot, fbound, iaez, idp, .true., rt, .true.)
      tturne0 は閉軌道を求めず (tcod も temit も呼ばない) cod0 から 6D で伝搬する。各要素の出口で tsetetwiss:
        ti = ri · inv(T_{1→l});  twiss(l+1) = tfetwiss(ti, cod_l, normali)   (normali は 1 回目の etwiss2ri の normal)
      位相は l=1 が 0 (tfetwiss(ri0) の位相は 0)、以降は +tfetwiss の位相 (<-toln なら +2π)、最後に tfadjustn。
      空洞の位相基準 trf0 は直前の temit/tcod が残した tmacro の trf0 (呼び出し側が trf0 に渡す)。
      放射は CELL の CALC6D と同じく RAD フラグに依らず tradke の経路 (_turn(radiation=True))。

    tw0: 長さ 28 の先頭 MARK の値 (index は TWISS_FUNS: AX BX PSIX AY BY PSIY EX EPX EY EPY R1..R4 DETR DX DPX DY DPY DZ DDP AZ BZ PSIZ ZX ZPX ZY ZPY)
    """
    els = list(lattice.elements)
    circ = sum(el.L for el in els)
    ph = _Phys(momentum_eV, circ, amass=amass, charge=charge)
    ph.dvfs = fshift if rfsw else 0.0
    if radiation is None:
        radiation = True
    tw0 = np.array(tw0, dtype=float)
    if tw0[22] <= 0.0:
        tw0[22] = 1.0
    ri0, normal0 = etwiss2ri(tw0)
    cod0 = tw0[15:21].copy()
    tw1 = tfetwiss(ri0, cod0, True)                  # twiss(1) = tfetwiss(ri0, cod, .true.)
    ri1, _ = etwiss2ri(tw1)                          # tfinibeam(1): ri = etwiss2ri(twiss(1))
    rfsw0 = _lat.RFSW
    _lat.RFSW = rfsw
    try:
        st = _turn(els, cod0, ph, trf0=trf0, rfsw=rfsw, radcod=radcod, radiation=radiation, record=True)
    finally:
        _lat.RFSW = rfsw0
    rows = _rows_from_turn(st, ri1, normal0, els, row0=tw1)
    return Calc6DResult(rows=rows, stable=True, closed_orbit=cod0, trf0=trf0, r=None, ri=ri1, tunes=None,
                        emit=None, elements=els, ph=ph, rfsw=rfsw, radcod=radcod, radiation=radiation, messages=[])


# ---------------------------------------------------------------------------
# 要素内の分数位置 (qtwissfrac の calc6d 分岐)
# ---------------------------------------------------------------------------

def twiss_frac(res: Calc6DResult, l: int, fr: float) -> TwissRow:
    """qtwiss.f qtwissfrac (calc6d=.true.) の直訳:
    要素 l (1 始まり) の入口から割合 fr の位置の拡張 Twiss。
    trans = inv(etwiss2ri(twiss(l))) を部分要素 (qfracseg → ffscmd.FFS.frac_element) で伝搬し、
    ftwiss = tfetwiss(inv(trans), cod, normal)、位相は twiss(l) の NX/NY/NZ を足す。
    irad=6 なので放射の行列部は積まないが T (シンプレクティック部) は同じ。"""
    rows = res.rows
    if fr == 0.0 or l >= len(rows):
        return rows[l - 1]
    from .ffscmd import FFS          # 循環 import を避けるため遅延 import
    el = res.elements[l - 1]
    part = FFS.frac_element(el, fr)
    row = rows[l - 1]
    tw_l = tw_from_row(row)
    ri_l, normal = etwiss2ri(tw_l)
    cod_l = tw_l[15:21].copy()
    rfsw0 = _lat.RFSW
    _lat.RFSW = res.rfsw
    try:
        st = _turn([part], cod_l, res.ph, trf0=res.trf0, rfsw=res.rfsw, radcod=res.radcod,
                   radiation=res.radiation)
    finally:
        _lat.RFSW = rfsw0
    # trans = T_part · inv(ri_l) → tfetwiss(inv(trans)) = tfetwiss(ri_l · inv(T_part))
    tw = tfetwiss(ri_l @ _tinv6(st.T), st.cod, normal)
    tw[2] += tw_l[2]
    tw[5] += tw_l[5]
    tw[23] += tw_l[23]
    return row_from_tw(tw, row.s + part.L, el.name)


__all__ = ['calc6d', 'Calc6DResult', 'twiss_frac', 'tfetwiss', 'etwiss2ri',
           'row_from_tw', 'tw_from_row']
