"""
pysad.lattice
=============
加速器格子クラス。SAD の LINE[] / Twiss[] / CalculateEmittance[] に相当。

使い方:
    from pysad.elements import make_element
    from pysad.lattice import Lattice

    elems = [
        make_element('QUAD', 'QF', L=0.3, K1=1.2),
        make_element('DRIFT', 'D1', L=2.5),
        make_element('BEND', 'B1', L=1.0, ANGLE=0.1),
        make_element('DRIFT', 'D2', L=2.5),
        make_element('QUAD', 'QD', L=0.3, K1=-1.2),
        make_element('DRIFT', 'D3', L=2.5),
    ]
    lat = Lattice(elems)
    twiss = lat.twiss()
    emit = lat.emittance(energy_eV=3e9)
"""

from __future__ import annotations

import math
import os
import numpy as np
from numpy.linalg import eig, inv
from typing import List, Optional, Dict
from dataclasses import dataclass, field

from .elements import Element, Bend, Quadrupole, Drift
from .optics import periodic_optics, periodic_r, symplectic_j
from .constants import (SpeedOfLight, ElectronRadius, ElectronMass,
                        Cgamma, Cq)


# ---------------------------------------------------------------------------
# SAD mfit* 列定数 (SAD の qtwiss.f で定義される整数インデックス)
# ---------------------------------------------------------------------------
# twiss(nlat, -1:1, 28) の第 3 次元インデックス (1-based → 0-based はアクセス時に -1)
MFITAX  = 1;  MFITBX  = 2;  MFITNX  = 3    # αx, βx, ψx
MFITAY  = 4;  MFITBY  = 5;  MFITNY  = 6    # αy, βy, ψy
MFITEX  = 7;  MFITEPX = 8                  # ηx, η'x
MFITEY  = 9;  MFITEPY = 10                 # ηy, η'y
MFITR1  = 11; MFITR2  = 12; MFITR3  = 13   # R 行列成分
MFITR4  = 14; MFITDETR = 15                # R 行列 det
MFITDX  = 16; MFITDPX = 17                 # x_COD, px_COD
MFITDY  = 18; MFITDPY = 19                 # y_COD, py_COD
MFITDZ  = 20; MFITDDP = 21                 # z_COD, δ_COD
MFITALZ = 22; MFITBZ  = 23; MFITMUZ = 24  # αz, βz, ψz
MFITZX  = 25; MFITZPX = 26                 # z-x 結合
MFITZY  = 27; MFITZPY = 28                 # z-y 結合


@dataclass
class TwissRow:
    """
    SAD の twiss(i, j, 28) 配列の 1 行分。

    SAD の mfit* 定数と同じ列順序 (1-based) に対応。
    LINEのDX/DYは要素の据え付けオフセット (Quadrupole.DX/DY)。
    TwissのDX/DYはCOD/バンプ軌道オフセット (x_cod, y_cod)。

    Attributes
    ----------
    s    : s 座標 [m]  (要素出口位置)
    name : 要素名
    alfx, betx, mux   : α_x, β_x [m], ψ_x [rad]   (mfitax/bx/nx = 1/2/3)
    alfy, bety, muy   : α_y, β_y [m], ψ_y [rad]   (mfitay/by/ny = 4/5/6)
    etx, etpx         : η_x [m], η'_x              (mfitex/epx  = 7/8)
    ety, etpy         : η_y [m], η'_y              (mfitey/epy  = 9/10)
    r1..r4, detr      : 結合行列成分               (mfitr1..detr = 11..15)
    x_cod, px_cod     : 閉軌道変位 x,px            (mfitdx/dpx  = 16/17)
    y_cod, py_cod     : 閉軌道変位 y,py            (mfitdy/dpy  = 18/19)
    z_cod, delta_cod  : 閉軌道変位 z,δ             (mfitdz/ddp  = 20/21)
    alfz, betz, muz   : 縦方向 Twiss              (mfitalz/bz/muz = 22/23/24)
    zx, zpx, zy, zpy  : z-横方向結合              (mfitzx..zpy = 25..28)
    """
    s:         float = 0.0
    name:      str   = ''
    # col 1-3: x Twiss
    alfx:  float = 0.0
    betx:  float = 1.0
    mux:   float = 0.0
    # col 4-6: y Twiss
    alfy:  float = 0.0
    bety:  float = 1.0
    muy:   float = 0.0
    # col 7-10: 分散
    etx:   float = 0.0
    etpx:  float = 0.0
    ety:   float = 0.0
    etpy:  float = 0.0
    # col 11-15: 結合行列
    r1:    float = 0.0
    r2:    float = 0.0
    r3:    float = 0.0
    r4:    float = 0.0
    detr:  float = 0.0
    # col 16-21: 閉軌道 (Twiss DX/DY は x_cod/y_cod に相当)
    x_cod:     float = 0.0
    px_cod:    float = 0.0
    y_cod:     float = 0.0
    py_cod:    float = 0.0
    z_cod:     float = 0.0
    delta_cod: float = 0.0
    # col 22-24: 縦方向 Twiss
    alfz:  float = 0.0
    betz:  float = 1.0      # SAD の 4 次元 CALC/INS では BZ は常に 1 (Twiss["*"] で確認)
    muz:   float = 0.0
    # col 25-28: z-横結合
    zx:    float = 0.0
    zpx:   float = 0.0
    zy:    float = 0.0
    zpy:   float = 0.0

    def to_dict(self) -> dict:
        """辞書形式に変換 (後方互換性のため)"""
        return {
            's': self.s, 'name': self.name,
            'ax': self.alfx, 'bx': self.betx, 'nx': self.mux,
            'ay': self.alfy, 'by': self.bety, 'ny': self.muy,
            'ex': self.etx,  'epx': self.etpx,
            'ey': self.ety,  'epy': self.etpy,
            'r1': self.r1, 'r2': self.r2, 'r3': self.r3, 'r4': self.r4,
            'detr': self.detr,
            'x_cod': self.x_cod, 'px_cod': self.px_cod,
            'y_cod': self.y_cod, 'py_cod': self.py_cod,
            'z_cod': self.z_cod, 'delta_cod': self.delta_cod,
            # SAD alias names
            'alfx': self.alfx, 'betx': self.betx, 'mux': self.mux,
            'alfy': self.alfy, 'bety': self.bety, 'muy': self.muy,
            'etx': self.etx,   'etpx': self.etpx,
            'ety': self.ety,   'etpy': self.etpy,
        }


def nan_rows(lattice):
    """qcod が溢れて直前の twiss も無いときの表: 全要素 + 末尾 ($$$) の行を NaN で埋める (name と s だけ入れる)。
    SAD は壊れた値 (1e200 や前回の値) で続行するので、pysad は NaN で続行して Twiss[] を NaN にする"""
    import dataclasses
    nan = float('nan')
    rows = []
    s = 0.0
    names = [el.name for el in lattice.elements] + ['$$$']
    for k, nm in enumerate(names):
        r = TwissRow(name=nm, s=s)
        for f in dataclasses.fields(TwissRow):
            if f.type in (float, 'float') and f.name != 's':
                setattr(r, f.name, nan)
        rows.append(r)
        if k < len(lattice.elements):
            s += float(getattr(lattice.elements[k], 'L', 0.0) or 0.0)
    return rows


# ---------------------------------------------------------------------------
# 補助関数
# ---------------------------------------------------------------------------

def _transport_6x6(elements: List[Element],
                   delta: float = 0.0) -> np.ndarray:
    """要素リストを順に転送行列を掛け合わせて全体の 6×6 行列を返す。"""
    M = np.eye(6)
    for el in elements:
        M = el.transfer_matrix(delta) @ M
    return M


def _twiss_at_elements(M_one_turn: np.ndarray,
                       elements: List[Element],
                       twiss0: dict) -> List[dict]:
    """
    閉軌道の Twiss パラメータを各要素出口で計算する。

    Parameters
    ----------
    M_one_turn : ndarray (6,6) – 一周転送行列
    elements   : list of Element
    twiss0     : dict – 入口 (=出口) の Twiss パラメータ

    Returns
    -------
    twiss_list : list of dict, len = len(elements)+1
        twiss_list[0]   = twiss0 (入口)
        twiss_list[i]   = i 番目要素の出口
    """
    twiss_list = [twiss0]
    tw_cur = twiss0
    for el in elements:
        Me = el.transfer_matrix()
        # 前の要素出口の Twiss を初期値として要素行列で伝搬する。
        # （累積行列に arctan2 を使うと ±π で折り返すため誤った位相累積になる）
        tw_out = _transport_twiss(tw_cur, Me)
        twiss_list.append(tw_out)
        tw_cur = tw_out
    return twiss_list


_XYTH = 0.375      # tffs.f:210 xyth (x-y モードの入れ替え判定)


from .perfkern import twiss_chain as _twiss_chain_nb, TW_KEYS as _TW_KEYS


def _chain(tw_cur: dict, el_mats, clamp: bool, check: bool):
    """_transport_twiss を el_mats の全要素に順に適用する (numba)。戻り (out (n,15), nbad)"""
    vec = np.array([float(tw_cur.get(k, 0.0)) for k in _TW_KEYS])
    if not el_mats:
        return np.empty((0, 15)), -1
    mats = np.array([m[0] for m in el_mats])
    return _twiss_chain_nb(vec, mats, clamp, check)


def _transport_twiss(twiss0: dict, M: np.ndarray, clamp: bool = True, gr: float = 1.0) -> dict:
    """
    Twiss パラメータを転送行列 M で伝搬する。SAD qtwiss.f の qmat2twiss の直訳 (結合 R と分散を含む)。

    M は 6×6 (行 0-3 と列 0-3, 5 を SAD の trans(4,5) = x px y py と δ 列として使う)。
    要素行列に x-y 結合ブロック (M[0:2,2:4], M[2:4,0:2]) があれば結合の式 (coup=.true.)、なければ
    非結合の式 (R だけ t1..t4 で更新) を使う。twiss0['detr'] < 1 なら normal (x モードが x 面)、
    そうでなければ x-y が入れ替わったモード (SAD の `normal=twiss(ipa,mfitdetr) < 1`)。

    Parameters
    ----------
    twiss0 : dict  入口 Twiss パラメータ (ax bx ay by nx ny ex epx ey epy r1..r4 detr)
    M      : ndarray (6,6) 転送行列

    gr     : float  qmat2twiss の trpt の gr = gammab(l)/gammab(l+1) (加速要素)。非結合の枝で入口の分散に掛ける
                    (結合の枝では SAD も掛けない)

    Returns
    -------
    twiss_out : dict  出口 Twiss パラメータ
    """
    ax0, bx0 = twiss0['ax'], twiss0['bx']
    ay0, by0 = twiss0['ay'], twiss0['by']
    nx0, ny0 = twiss0.get('nx', 0.0), twiss0.get('ny', 0.0)
    ex0, epx0 = twiss0.get('ex', 0.0), twiss0.get('epx', 0.0)
    ey0, epy0 = twiss0.get('ey', 0.0), twiss0.get('epy', 0.0)
    r1, r2 = twiss0.get('r1', 0.0), twiss0.get('r2', 0.0)
    r3, r4 = twiss0.get('r3', 0.0), twiss0.get('r4', 0.0)
    detr0 = twiss0.get('detr', 0.0)
    normal = detr0 < 1.0

    t = np.empty((4, 5))
    t[:, :4] = M[:4, :4]
    t[:, 4] = M[:4, 5]
    coup = bool(np.any(t[0:2, 2:4]) or np.any(t[2:4, 0:2]))

    if coup:
        sq = math.sqrt(max(0.0, 1.0 - r1 * r4 + r2 * r3))
        u = np.empty((4, 5))
        if normal:
            r15 = sq * ex0 + r4 * ey0 - r2 * epy0
            r25 = sq * epx0 - r3 * ey0 + r1 * epy0
            r35 = sq * ey0 - r1 * ex0 - r2 * epx0
            r45 = sq * epy0 - r3 * ex0 - r4 * epx0
            for i in range(4):
                a1, a2, a3, a4, a5 = t[i]
                u[i, 0] = a1 * sq - a3 * r1 - a4 * r3
                u[i, 1] = a2 * sq - a3 * r2 - a4 * r4
                u[i, 2] = a1 * r4 - a2 * r3 + a3 * sq
                u[i, 3] = -a1 * r2 + a2 * r1 + a4 * sq
                u[i, 4] = a1 * r15 + a2 * r25 + a3 * r35 + a4 * r45 + a5
        else:
            r15 = -r1 * ex0 - r2 * epx0 + sq * ey0
            r25 = -r3 * ex0 - r4 * epx0 + sq * epy0
            r35 = sq * ex0 + r4 * ey0 - r2 * epy0
            r45 = sq * epx0 - r3 * ey0 + r1 * epy0
            for i in range(4):
                a1, a2, a3, a4, a5 = t[i]
                u[i, 0] = -a1 * r1 - a2 * r3 + a3 * sq
                u[i, 1] = -a1 * r2 - a2 * r4 + a4 * sq
                u[i, 2] = a1 * sq + a3 * r4 - a4 * r3
                u[i, 3] = a2 * sq - a3 * r2 + a4 * r1
                u[i, 4] = a1 * r15 + a2 * r25 + a3 * r35 + a4 * r45 + a5
        (u11, u12, u13, u14, u15), (u21, u22, u23, u24, u25), \
            (u31, u32, u33, u34, u35), (u41, u42, u43, u44, u45) = u
        detp = (u11 * u22 - u21 * u12 + u33 * u44 - u43 * u34) * 0.5
        if detp > _XYTH:
            sd = math.sqrt(detp)
            u11 /= sd; u12 /= sd; u21 /= sd; u22 /= sd
            u33 /= sd; u34 /= sd; u43 /= sd; u44 /= sd
            r1n = -u31 * u22 + u32 * u21
            r2n = u31 * u12 - u32 * u11
            r3n = -u41 * u22 + u42 * u21
            r4n = u41 * u12 - u42 * u11
            detr = r1n * r4n - r2n * r3n
            sd = math.sqrt(1.0 - detr)
            ex1 = u15 * sd - r4n * u35 + r2n * u45
            epx1 = u25 * sd + r3n * u35 - r1n * u45
            ey1 = u35 * sd + r1n * u15 + r2n * u25
            epy1 = u45 * sd + r3n * u15 + r4n * u25
            aa = bx0 * u11 - ax0 * u12
            cc = (u12 - ax0 * aa) / bx0
            ax1 = -(u21 * aa + u22 * cc)
            bx1 = u11 * aa + u12 * cc
            dpsix = math.atan2(u12, aa)
            bb = by0 * u33 - ay0 * u34
            dd = (u34 - ay0 * bb) / by0
            ay1 = -(u43 * bb + u44 * dd)
            by1 = u33 * bb + u34 * dd
            dpsiy = math.atan2(u34, bb)
        else:
            detm = (u31 * u42 - u32 * u41 + u13 * u24 - u14 * u23) * 0.5
            sd = math.sqrt(abs(detm))
            u31 /= sd; u32 /= sd; u41 /= sd; u42 /= sd
            u13 /= sd; u14 /= sd; u23 /= sd; u24 /= sd
            r1n = -u11 * u42 + u12 * u41
            r2n = u11 * u32 - u12 * u31
            r3n = -u21 * u42 + u22 * u41
            r4n = u21 * u32 - u22 * u31
            sd = math.sqrt(1.0 - r1n * r4n + r2n * r3n)
            ex1 = u35 * sd - r4n * u15 + r2n * u25
            epx1 = u45 * sd + r3n * u15 - r1n * u25
            ey1 = u15 * sd + r1n * u35 + r2n * u45
            epy1 = u25 * sd + r3n * u35 + r4n * u45
            detr = 1.0 + _XYTH - r1n * r4n + r2n * r3n
            aa = bx0 * u31 - ax0 * u32
            cc = (u32 - ax0 * aa) / bx0
            ax1 = -(u41 * aa + u42 * cc)
            bx1 = u31 * aa + u32 * cc
            dpsix = math.atan2(u32, aa)
            bb = by0 * u13 - ay0 * u14
            dd = (u14 - ay0 * bb) / by0
            ay1 = -(u23 * bb + u24 * dd)
            by1 = u13 * bb + u14 * dd
            dpsiy = math.atan2(u14, bb)
    else:
        (a11, a12, a13, a14, a15), (a21, a22, a23, a24, a25), \
            (a31, a32, a33, a34, a35), (a41, a42, a43, a44, a45) = t
        if gr != 1.0:
            ex0, epx0, ey0, epy0 = ex0 * gr, epx0 * gr, ey0 * gr, epy0 * gr
        if normal:
            t1 = r1 * a22 - r2 * a21
            t2 = -r1 * a12 + r2 * a11
            t3 = r3 * a22 - r4 * a21
            t4 = -r3 * a12 + r4 * a11
            r1n = a33 * t1 + a34 * t3
            r2n = a33 * t2 + a34 * t4
            r3n = a43 * t1 + a44 * t3
            r4n = a43 * t2 + a44 * t4
            detr = r1n * r4n - r2n * r3n
            sd = math.sqrt(1.0 - detr)
            ex1 = a11 * ex0 + a12 * epx0 + sd * a15 - r4n * a35 + r2n * a45
            epx1 = a21 * ex0 + a22 * epx0 + sd * a25 + r3n * a35 - r1n * a45
            ey1 = a33 * ey0 + a34 * epy0 + sd * a35 + r1n * a15 + r2n * a25
            epy1 = a43 * ey0 + a44 * epy0 + sd * a45 + r3n * a15 + r4n * a25
            aa = bx0 * a11 - ax0 * a12
            cc = (a12 - ax0 * aa) / bx0
            ax1 = -(a21 * aa + a22 * cc)
            bx1 = a11 * aa + a12 * cc
            dpsix = math.atan2(a12, aa)
            bb = by0 * a33 - ay0 * a34
            dd = (a34 - ay0 * bb) / by0
            ay1 = -(a43 * bb + a44 * dd)
            by1 = a33 * bb + a34 * dd
            dpsiy = math.atan2(a34, bb)
        else:
            sd = math.sqrt(1.0 - r1 * r4 + r2 * r3)
            t1 = r1 * a44 - r2 * a43
            t2 = -r1 * a34 + r2 * a33
            t3 = r3 * a44 - r4 * a43
            t4 = -r3 * a34 + r4 * a33
            r1n = a11 * t1 + a12 * t3
            r2n = a11 * t2 + a12 * t4
            r3n = a21 * t1 + a22 * t3
            r4n = a21 * t2 + a22 * t4
            detr = 1.0 + _XYTH - r1n * r4n + r2n * r3n
            ex1 = a33 * ex0 + a34 * epx0 + sd * a35 + (-r4n * a15 + r2n * a25)
            epx1 = a43 * ex0 + a44 * epx0 + sd * a45 + (r3n * a15 - r1n * a25)
            ey1 = a11 * ey0 + a12 * epy0 + sd * a15 + (r1n * a35 + r2n * a45)
            epy1 = a21 * ey0 + a22 * epy0 + sd * a25 + (r3n * a35 + r4n * a45)
            aa = bx0 * a33 - ax0 * a34
            cc = (a34 - ax0 * aa) / bx0
            ax1 = -(a43 * aa + a44 * cc)
            bx1 = a33 * aa + a34 * cc
            dpsix = math.atan2(a34, aa)
            bb = by0 * a11 - ay0 * a12
            dd = (a12 - ay0 * bb) / by0
            ay1 = -(a21 * bb + a22 * dd)
            by1 = a11 * bb + a12 * dd
            dpsiy = math.atan2(a12, bb)

    return {
        'ax': ax1, 'bx': max(bx1, 1e-12) if clamp else bx1,
        'ay': ay1, 'by': max(by1, 1e-12) if clamp else by1,
        'nx': nx0 + dpsix, 'ny': ny0 + dpsiy,
        'ex': ex1, 'epx': epx1,
        'ey': ey1, 'epy': epy1,
        'r1': r1n, 'r2': r2n, 'r3': r3n, 'r4': r4n,
        'detr': detr,
    }


def _one_turn_with_distortion(elements: List[Element],
                              delta: float = 0.0):
    """
    全要素の一周転送行列 M (6×6) と定常軌道歪みベクトル d (6,) を計算する。

    各要素を x_{i+1} = M_i @ x_i + d_i と表したとき、
    全体では x_N = M_total @ x_0 + d_total。

    d_i は要素の DX/DY フィードダウンおよび K0/ROTATE 補正コイルキック
    (element.constant_kick(delta)) から得る。

    Returns
    -------
    M_total : ndarray (6, 6)   一周転送行列
    d_total : ndarray (6,)     一周の定常キック累積
    """
    M_total = np.eye(6)
    d_total = np.zeros(6)
    for el in elements:
        M_el = el.transfer_matrix(delta)
        d_el = el.constant_kick(delta)
        d_total = M_el @ d_total + d_el
        M_total = M_el @ M_total
    return M_total, d_total


# ---------------------------------------------------------------------------
# SADのqcod相当: COD依存DRIFT Jacobianと周期閉軌道の計算
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# SADのtbendebody0・tbedge COD更新 (Fortranソースからの直接翻訳)
# ---------------------------------------------------------------------------

def _pxy2dpz(pxi: float, pyi: float, p: float) -> float:
    """SADのpxy2dpz: dpz = sqrt(p^2 - pxi^2 - pyi^2) - p"""
    arg = p**2 - pxi**2 - pyi**2
    return np.sqrt(max(0.0, arg)) - p


def _tbedgedrift_cod(cod: np.ndarray, psi: float) -> np.ndarray:
    """SADのtbedgedrift COD更新のみ (Fortranソース直訳)。
    tbfrie.f subroutine tbedgedrift"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    cosp = np.cos(psi)
    sinp = np.sin(psi)
    pxi = np.clip(cod[1], -pr, pr)
    pyi = np.clip(cod[3], -pr, pr)
    dpzi = _pxy2dpz(pxi, pyi, pr)
    pzi = pr + dpzi
    f = cosp - pxi / pzi * sinp
    cod[0] = cod[0] / f         # cod(1)=cod(1)/f
    sx = cod[0] * sinp          # sx=cod(1)*sinp  (updated cod(1))
    cod[1] = pxi * cosp + pzi * sinp  # cod(2)
    cod[2] = cod[2] + pyi / pzi * sx  # cod(3)
    cod[4] = cod[4] - pr / pzi * sx   # cod(5)
    return cod


def _tbedgebody_cod(cod: np.ndarray, rhob: float, psi: float) -> np.ndarray:
    """SADのtbedgebody COD更新のみ (Fortranソース直訳)。
    tbfrie.f subroutine tbedgebody"""
    import math
    cod = cod.copy()
    pr = 1.0 + cod[5]
    cosp = np.cos(psi)
    sinp = np.sin(psi)
    sinsq = 2.0 * np.sin(0.5 * psi)**2   # = 1 - cos(psi)
    ampmin = 1e-6
    ampmax = 1.0 - ampmin
    pxi = np.clip(cod[1], -pr, pr)
    pyi = np.clip(cod[3], -pr, pr)
    s = min(ampmax, pxi**2 + pyi**2)
    sqrs = pr * np.sqrt(max(ampmin, 1.0 - s / pr**2))
    dpzi = -s / (pr + sqrs)
    pzi = pr + dpzi
    phsq = pzi**2 + pxi**2
    px0 = pxi * cosp - pzi * sinp
    dpz0 = pxi * sinp + dpzi * cosp - pr * sinsq
    xi = cod[0]
    pxf = px0 + xi * sinp / rhob
    cod[1] = np.clip(pxf, -pr, pr)
    s1 = min(ampmax, pyi**2 + pxf**2)
    sqrs1 = pr * np.sqrt(max(ampmin, 1.0 - s1 / pr**2))
    dpzf = -s1 / (pr + sqrs1)
    ddpz = (s - s1) / (sqrs + sqrs1) + (dpzi + pr) * sinsq - pxi * sinp
    cod[0] = rhob * ddpz + xi * cosp
    a = pr * xi * sinp / rhob + dpz0 * pxf - dpzf * px0
    val = np.clip(a / phsq, -ampmax, ampmax)
    theta = math.asin(val)
    cod[2] = cod[2] - rhob * theta * pyi
    cod[4] = cod[4] + rhob * theta * pr
    return cod


def _tbedgemaxwell_cod(cod: np.ndarray, rhob: float) -> np.ndarray:
    """SADのtbedgemaxwell COD更新のみ (Fortranソース直訳)。
    tbfrie.f subroutine tbedgemaxwell"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    pxi = np.clip(cod[1], -pr, pr)
    pvi2 = (pr - pxi) * (pr + pxi)   # pr^2 - pxi^2
    pvi2 = max(1e-30, pvi2)
    pvi = np.sqrt(pvi2)
    yi = cod[2]
    yr = (yi / rhob)**2
    a = 1.0 - yr / 6.0
    af = a * yi / rhob / pvi
    b = 0.5 - yr / 24.0
    bf = b * yi**2 / rhob / pvi2 / pvi
    cod[0] = cod[0] + bf * pr
    cod[3] = np.clip(cod[3] - af * pxi, -pr, pr)
    cod[4] = cod[4] - bf * pxi
    return cod


def _tbedge_cod(cod: np.ndarray, rhob: float, psi: float,
                is_entrance: bool) -> np.ndarray:
    """SADのtbedge COD更新 (Fortranソース直訳)。
    tbfrie.f subroutine tbedge
    入口 (is_entrance=True):  tbedgedrift → tbedgemaxwell → tbedgebody
    出口 (is_entrance=False): tbedgebody → tbedgemaxwell(-rhob) → tbedgedrift"""
    if psi == 0.0:
        return _tbedgemaxwell_cod(cod, rhob if is_entrance else -rhob)
    if is_entrance:
        cod = _tbedgedrift_cod(cod, psi)
        cod = _tbedgemaxwell_cod(cod, rhob)
        cod = _tbedgebody_cod(cod, rhob, psi)
    else:
        cod = _tbedgebody_cod(cod, rhob, psi)
        cod = _tbedgemaxwell_cod(cod, -rhob)
        cod = _tbedgedrift_cod(cod, psi)
    return cod


def _tbendebody0_cod(cod: np.ndarray, phi0: float, rhob: float,
                     drhob: float = 0.0) -> np.ndarray:
    """SADのtbendebody0 COD更新のみ (K1=0, drhob=0) (Fortranソース直訳)。
    tbende.f subroutine tbendebody0
    phi0  : 偏向角 [rad]  (= phin = phic/ndiv)
    rhob  : 設計曲率半径 = L/phi0
    drhob = 0 (理想セクターベンド: rhob==rho0)"""
    cod = cod.copy()
    xi  = cod[0]; pxi = cod[1]
    yi  = cod[2]; pyi = cod[3]
    dp  = cod[5]
    p   = 1.0 + dp
    # SADのdcn = cos(phi0) - 1  (負の値)
    cn  = np.cos(phi0);  sn  = np.sin(phi0);  dcn = cn - 1.0
    # dpzi = p * pxy2dpz(pxi/p, pyi/p)
    dpzi = _pxy2dpz(pxi, pyi, p)  # = sqrt(p^2-pxi^2-pyi^2) - p
    pzi  = p + dpzi
    # drhob=0 の場合: dpza = dpzi + (drhob - xi)/rhob = dpzi - xi/rhob
    dpza = dpzi + (drhob - xi) / rhob     # tbende.f: dpza=dpzi+(drhob-xi)/rhob
    pza  = dpza + dp
    dpxf = dcn * pxi + sn * pza   # dcn = cos-1 (SADのdcn)
    pxf  = pxi + dpxf             # = cos*pxi + sin*pza
    dpzf = _pxy2dpz(pxf, pyi, p)  # = sqrt(p^2-pxf^2-pyi^2) - p
    pzf  = p + dpzf
    # xf = rhob*(sn*pxi - cn*dpzi - dcn*dp + dpzf) + cn*xi  (drhob=0)
    xf   = rhob * (sn*pxi - cn*dpzi - dcn*dp + dpzf) + cn*xi - dcn*drhob
    phsq = pzf**2 + pxf**2
    if abs(phsq) > 1e-30:
        val = (-p*dpxf - pxf*dpzi + pxi*dpzf) / phsq
        da  = np.arcsin(np.clip(val, -1.0, 1.0))
    else:
        da = 0.0
    omega = phi0 + da
    yf    = yi + rhob * omega * pyi
    # zf = cod(5) - da*p*rhob - phi0*(drhob + dp*rhob)  drhob=0:
    zf    = cod[4] - da * p * rhob - phi0 * (drhob + dp * rhob)
    cod[0] = xf
    cod[1] = np.clip(pxf, -p, p)
    cod[2] = yf
    cod[3] = np.clip(pyi, -p, p)
    cod[4] = zf
    return cod


def _bend_fringe_params(el) -> tuple:
    """BEND の (fb1, fb2, mfring, fringe): fb1 = F1+FB1, fb2 = F1+FB2, mfring = FRINGE, fringe = (DISFRIN==0)。
    逆向き要素 (dir<0) の入れ替えは未対応。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    return (g('F1') + g('FB1'), g('F1') + g('FB2'), int(round(g('FRINGE'))), g('DISFRIN') == 0.0)


def _tblfre_cod(cod: np.ndarray, dxfr: float, dyfr: float, dyfra: float) -> np.ndarray:
    """SAD tbende.f tblfre の COD 更新のみ (BEND の F1 フリンジの線形部)。"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    ysq = cod[2]**2
    dyfraysq = ysq * dyfra
    x_new = cod[0] + dxfr * cod[5] / pr
    py_new = cod[3] + (dyfr - dyfraysq) * cod[2] / pr
    cod[4] = cod[4] + (dxfr * cod[1] + (0.5 * dyfr - 0.25 * dyfraysq) * ysq) / pr**2
    cod[0] = x_new
    cod[3] = py_new
    return cod


_RBH = 0.5 + 1.0 / math.sqrt(12.0)
_RBL = 1.0 / 6.0 / _RBH


def _bend_cod_update(cod: np.ndarray, L: float, phi0: float,
                     E1: float, E2: float, phib: float = None,
                     fb1: float = 0.0, fb2: float = 0.0, mfring: int = 0,
                     fringe: bool = True) -> np.ndarray:
    """SADのtbende COD更新 (K1=0 の分岐、入口エッジ→ボディ→出口エッジ)。
    Fortranソース: tbende.f subroutine tbende (K1=0ブランチ), tblfre, tbendef1
    E1, E2 : 入口・出口エッジ角 [rad] (= pe1, pe2 = psi1*phi0, psi2*phi0)
    fb1, fb2: F1 フリンジの長さ (F1+FB1, F1+FB2)、mfring = FRINGE、fringe = (DISFRIN==0)"""
    if abs(L) < 1e-30 or abs(phi0) < 1e-30:
        return cod.copy()
    # SAD tbende: phib = ANGLE + K0 (実際の曲げ), phi0 = ANGLE (基準)。
    if phib is None:
        phib = phi0
    al0 = L
    # alcorr (qbend/tturne とも .true.): F1 フリンジがあるときの長さ補正 (tbende.f:381-386)
    if mfring != 0 and phi0 != 0.0:
        psi1, psi2 = E1 / phi0, E2 / phi0
        al = al0 - ((phi0 * fb1)**2 + (phi0 * fb2)**2) / 48.0 / al0 \
            * math.sin(0.5 * (phi0 * (1.0 - psi1 - psi2))) / math.sin(0.5 * phi0)
    else:
        al = al0
    phibl = phib / al
    rhob = 1.0 / phibl
    rho0 = al / phi0
    drhob = rhob - rho0
    f1r = 0.0
    f2r = 0.0
    if fb1 != 0.0 and (mfring > 0 or mfring == -1):
        dxfr1 = fb1**2 * phibl / 24.0
        dyfr1 = fb1 * phibl**2 / 6.0
        dyfra1 = 4.0 * dyfr1 / fb1**2 if fringe else 0.0
        cod = _tblfre_cod(cod, dxfr1, dyfr1, dyfra1)
        f1r = 0.5 * fb1
    if (fb2 != 0.0 and mfring > 0) or mfring == -2:
        f2r = 0.5 * fb2
    fs = f1r + f2r
    rbc = 1.0 - fs / al0
    phic = phi0 * rbc
    # 入口エッジ (ent=.true.)
    cod = _tbedge_cod(cod, rhob, E1, is_entrance=True)

    def ef1(c, fb, rb):
        # tbendef1: 長さ fb/2 を比 rb, 1/(6 rb) の 2 区間に分けて tbendebody0
        for _ in range(2):
            fl = 0.5 * fb * rb
            c = _tbendebody0_cod(c, phi0 * fl / al0, rhob, drhob)
            rb = 1.0 / 6.0 / rb
        return c
    if f1r != 0.0:
        cod = ef1(cod, fb1, _RBL)
    # ボディ (K1=0, ndiv=1: phin=phic)
    cod = _tbendebody0_cod(cod, phic, rhob, drhob)
    if f2r != 0.0:
        cod = ef1(cod, fb2, _RBH)
    # 出口エッジ (ent=.false.)
    cod = _tbedge_cod(cod, rhob, E2, is_entrance=False)
    if f2r != 0.0:
        dxfr2 = -fb2**2 / rhob / 24.0
        dyfr2 = fb2 / rhob**2 / 6.0
        dyfra2 = 4.0 * dyfr2 / fb2**2 if fringe else 0.0
        cod = _tblfre_cod(cod, dxfr2, dyfr2, dyfra2)
    return cod



# ---------------------------------------------------------------------------
# BEND (K1=0) の解析転送行列: SAD tbfrie.f / tbende.f の trans1 の直訳
# (COD 写像の厳密な微分ではなく SAD が使う解析式。軌道が大きいとき数値微分と 1e-7 程度ずれる。
#  Fortran の trans(i,j) は 1 始まりなので Python では [i-1, j-1]。z 行は 4 次元の CALC では使われない)
# ---------------------------------------------------------------------------

def _tbedgemaxwell_mat(cod: np.ndarray, rhob: float) -> np.ndarray:
    """tbfrie.f tbedgemaxwell の trans (cod は更新前)。"""
    pr = 1.0 + cod[5]
    pxi = min(pr, max(-pr, cod[1]))
    pvi2 = (pr - pxi) * (pr + pxi)
    pvi = pr * math.sqrt(pvi2 / pr**2)
    yi = cod[2]
    yr = (yi / rhob)**2
    a = 1.0 - yr / 6.0
    af = a * yi / rhob / pvi
    b = 0.5 - yr / 24.0
    bf = b * yi**2 / rhob / pvi2 / pvi
    c = 1.0 - yr / 2.0
    t = np.eye(6)
    t[0, 1] = 3.0 * bf * pr**2 * pxi / pvi2
    t[0, 2] = af * pr**2 / pvi2
    t[0, 5] = -bf * pr * (pr**2 + 2.0 * pxi**2) / pvi2
    t[3, 1] = -t[0, 2]
    t[3, 2] = -c * pxi / rhob / pvi
    t[3, 5] = af * pr * pxi / pvi2
    t[4, 1] = t[0, 5]
    t[4, 2] = -t[3, 5]
    t[4, 5] = bf * pxi * (2.0 * pr**2 + pxi**2) / pvi2
    return t


def _tbedgedrift_mat(cod: np.ndarray, psi: float) -> np.ndarray:
    """tbfrie.f tbedgedrift の trans (cod は更新前)。"""
    pr = 1.0 + cod[5]
    sinp, cosp = math.sin(psi), math.cos(psi)
    pxi = min(pr, max(-pr, cod[1]))
    pyi = min(pr, max(-pr, cod[3]))
    a = pxi**2 + pyi**2
    dpzi = pr * (math.sqrt(max(0.0, 1.0 - a / pr**2)) - 1.0)
    s = a
    pzi = pr + dpzi
    f = cosp - pxi / pzi * sinp
    x1 = cod[0] / f
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
    return t


def _tbedgebody_mat(cod: np.ndarray, rhob: float, psi: float) -> np.ndarray:
    """tbfrie.f tbedgebody の trans (cod は更新前)。"""
    ampmin = 1e-6
    ampmax = 1.0 - ampmin
    pr = 1.0 + cod[5]
    cosp, sinp = math.cos(psi), math.sin(psi)
    sinsq = 2.0 * math.sin(0.5 * psi)**2
    pxi = min(pr, max(-pr, cod[1]))
    pyi = min(pr, max(-pr, cod[3]))
    s = min(ampmax, pxi**2 + pyi**2)
    sqrs = pr * math.sqrt(max(ampmin, 1.0 - s / pr**2))
    dpzi = -s / (pr + sqrs)
    pzi = pr + dpzi
    phsq = pzi**2 + pxi**2
    px0 = pxi * cosp - pzi * sinp
    dpz0 = pxi * sinp + dpzi * cosp - pr * sinsq
    pz0 = pr + dpz0
    xi = cod[0]
    pxf = px0 + xi * sinp / rhob
    s1 = min(ampmax, pyi**2 + pxf**2)
    sqrs1 = pr * math.sqrt(max(ampmin, 1.0 - s1 / pr**2))
    dpzf = -s1 / (pr + sqrs1)
    pzf = pr + dpzf
    ddpz = (s - s1) / (sqrs + sqrs1) + (dpzi + pr) * sinsq - pxi * sinp
    a = pr * xi * sinp / rhob + dpz0 * pxf - dpzf * px0
    aa = a * cosp - (px0 * pxf + pz0 * pzf) * sinp
    c = (aa / phsq + sinp) / pzi / pzf
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
    t[2, 3] = -rhob * (pyi**2 * c + theta)
    t[2, 5] = pr * pyi * rhob * c
    t[3, 3] = 1.0
    t[4, 0] = t[1, 0] * t[0, 5] - t[0, 0] * t[1, 5]
    t[4, 1] = t[1, 1] * t[0, 5] - t[0, 1] * t[1, 5]
    t[4, 3] = t[1, 3] * t[0, 5] - t[0, 3] * t[1, 5] + t[2, 5]
    t[4, 4] = 1.0
    t[4, 5] = rhob * (theta - pr**2 * c)
    t[5, 5] = 1.0
    return t


def _tbendebody0_mat(cod: np.ndarray, phi0: float, rhob: float, drhob: float = 0.0) -> np.ndarray:
    """tbende.f tbendebody0 の trans1 (cod は更新前)。z 行の h0/h1emit**3*aln 項 (dvemit 由来) は省く。"""
    xi = cod[0]
    dp = cod[5]
    p = 1.0 + dp
    pxi = min(p, max(-p, cod[1]))
    yi = cod[2]
    pyi = min(p, max(-p, cod[3]))
    sn, cn = math.sin(phi0), math.cos(phi0)
    dcn = cn - 1.0
    dpzi = _pxy2dpz(pxi, pyi, p)
    pzi = p + dpzi
    dpza = dpzi + (drhob - xi) / rhob
    pza = dpza + dp
    dpxf = dcn * pxi + sn * pza
    pxf = pxi + dpxf
    dpzf = _pxy2dpz(pxf, pyi, p)
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
    xsn = phi0 - sn          # mathfun xsincos の xs = x − sin x (以前は符号が逆で、z 行 (5,6) が rhob·2(φ−sinφ) だけずれていた。4D の CALC は z 行を使わない)
    t[4, 5] = -rhob * (da + xsn + sn * dpzf / pzf + ddodp * p)
    return t


def _bend_edge_chain(c: np.ndarray, M: np.ndarray, rhob: float, psi: float, entrance: bool):
    """tbfrie.f tbedge の (COD, 解析行列): 入口は (psi≠0 なら drift → ) maxwell(rhob) (→ body)、
    出口は (psi≠0 なら body → ) maxwell(-rhob) (→ drift)。各段の行列はその段の入口の COD で評価し、M に左から掛ける。"""
    if psi == 0.0:
        stages = [(_tbedgemaxwell_mat, (rhob if entrance else -rhob,), _tbedgemaxwell_cod)]
    elif entrance:
        stages = [(_tbedgedrift_mat, (psi,), _tbedgedrift_cod),
                  (_tbedgemaxwell_mat, (rhob,), _tbedgemaxwell_cod),
                  (_tbedgebody_mat, (rhob, psi), _tbedgebody_cod)]
    else:
        stages = [(_tbedgebody_mat, (rhob, psi), _tbedgebody_cod),
                  (_tbedgemaxwell_mat, (-rhob,), _tbedgemaxwell_cod),
                  (_tbedgedrift_mat, (psi,), _tbedgedrift_cod)]
    for mat_f, args, cod_f in stages:
        M = mat_f(c, *args) @ M
        c = cod_f(c, *args)
    return c, M


def _bend_cod_matrix(cod: np.ndarray, L: float, phi0: float, E1: float, E2: float, phib: float = None):
    """K1=0・フリンジ F1=0・据え付け誤差なしの BEND の (出口の COD, SAD の解析転送行列)。
    tbende.f の K1=0 の分岐 (ndiv=1): tbedge(入口) → tbendebody0 → tbedge(出口) の各 trans を掛ける。
    COD の更新は _bend_cod_update と同じ関数を使い、行列は各段の入口の COD で評価する。"""
    if abs(L) < 1e-30 or abs(phi0) < 1e-30:
        return cod.copy(), None
    if phib is None:
        phib = phi0
    rhob = L / phib
    rho0 = L / phi0
    drhob = rhob - rho0
    c, M = _bend_edge_chain(cod.copy(), np.eye(6), rhob, E1, True)
    M = _tbendebody0_mat(c, phi0, rhob, drhob) @ M
    c = _tbendebody0_cod(c, phi0, rhob, drhob)
    c, M = _bend_edge_chain(c, M, rhob, E2, False)
    return c, M


# ---------------------------------------------------------------------------
# ANGLE=0 の BEND (ステアリング: tbende.f の phi0==0 の分岐 → tsteee.f / tthine nord=2)
#   K1=0: L=0 → 薄い双極キック px −= K0 (tthine(nord=2, ak=−phib_arg), phib_arg = −K0)、
#         K0=0 → ドリフト、他 → tsteee (入口エッジ → 本体 (rhob = −L/K0) → 出口エッジ。エッジは rhob = +L/K0)
#   K1≠0: K0=0 → tquade、他 → tmulte (akm(0)=K0, akm(1)=K1·L)
#   符号: K0>0 は px を減らす (ANGLE≠0 の phib = ANGLE+K0 が曲げを強める向きと同じ)。
# ---------------------------------------------------------------------------

_XTHASX = (2.0**-52 * 13312.0 / 231.0 / 6.0)**0.1     # tfloor.f xthasx (epsieee = 2^-52)


def _asinx(x: float) -> float:
    """mathfun asinx = asin(x) − x (小さい x は級数)"""
    if abs(x) < _XTHASX:
        x2 = x * x
        return x * x2 * (1.0 / 6.0 + x2 * (3.0 / 40.0 + x2 * (5.0 / 112.0 + x2 * (35.0 / 1152.0 + x2 * 63.0 / 2816.0))))
    return math.asin(min(1.0, max(-1.0, x))) - x


def _tsteee_body(cod: np.ndarray, alx: float, rhob: float, dvemit: float = 0.0):
    """tsteee.f の本体 1 スライス (COD 更新と trans1)。alx = スライス長、rhob = L/phib_arg (= −L/K0)。
    trans1 の (5,6) の h0/h1emit³·alx 項は CALC (4 次元) では使わないので省く (emit6d は数値微分の piece で組む)。"""
    pr = 1.0 + cod[5]
    pxi = cod[1] / pr
    pyi = cod[3] / pr
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
    t[4, 5] = t[4, 5] - (alx + dl) / pr
    c = cod.copy()
    c[0] = cod[0] + alx * spx / spz
    c[1] = pxf * pr
    c[2] = cod[2] + pyi * (alx + dl)
    c[4] = cod[4] - (dl + dvemit * alx)
    return c, t


def _tsteee_cod_matrix(el, cod: np.ndarray, dvemit: float = 0.0):
    """ANGLE=0, K1=0 の BEND (tsteee) の (出口 COD, 転送行列)。放射なし (ndiv=1)、F1 フリンジなし。
    DX/DY は tchge (ずらすだけ)。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    al = el.L
    k0 = g('K0')
    c = cod.copy()
    dx, dy = g('DX'), g('DY')
    c[0] -= dx
    c[2] -= dy
    M = np.eye(6)
    if al == 0.0:
        c[1] -= k0                      # tthine(nord=2, ak=K0): cod(2) −= aki (1/(1+δ) なし)
    elif k0 == 0.0:
        pr = 1.0 + c[5]
        M = _drift_matrix_at_cod(al, c[5], c[1], c[3])
        c = _exact_drift_cod(c, al)
    else:
        rhoe = al / k0                  # tbedge(al, −phib_arg = +K0)
        c, M = _bend_edge_chain(c, M, rhoe, g('E1'), True)
        M = _tsteee_body_mat(c, al, -rhoe) @ M
        c, _ = _tsteee_body(c, al, -rhoe, dvemit)
        c, M = _bend_edge_chain(c, M, rhoe, g('E2'), False)
    c = c.copy()
    c[0] += dx
    c[2] += dy
    return c, M


def _bend_as_mult(el):
    """ANGLE=0, K1≠0, K0≠0 の BEND を tbende の tmulte 呼び出し (akm(0)=K0, akm(1)=K1·L) と同じ MULT 相当の代理に変える。
    設計ロール ROTATE・据え付けロール ROTATION・DX/DY は呼び出し側 (_bend_element_cod_matrix / emit6d._tbende) が tchge で外から付けるので 0 にする。"""
    from types import SimpleNamespace
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    return SimpleNamespace(L=el.L, KnL={0: g('K0'), 1: g('K1') * el.L}, ROTATION=0.0, DROTATE=g('DROTATE'),
                           EPS=g('EPS'), DX=g('DX'), DY=g('DY'), DISFRIN=g('DISFRIN'), FRINGE=g('FRINGE'), K0FR=0.0)


def _tsteee_body_mat(cod, alx, rhob):
    return _tsteee_body(cod, alx, rhob)[1]


def _tqente_cod(cod: np.ndarray, al: float) -> np.ndarray:
    """SADのtqente COD更新 (bz=0, tqente.f直訳).
    非線形ドリフト: r=al*(-dpz)/(pz*pr), x+=r*px, y+=r*py。"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    pxi = cod[1]
    pyi = cod[3]
    a = pxi**2 + pyi**2
    if a > 0.0 and abs(pr) > 1e-15:
        arg = max(0.0, 1.0 - a / pr**2)
        dpz = -a / pr / (1.0 + np.sqrt(arg))
    else:
        dpz = 0.0
    pz = pr + dpz
    if abs(pz) < 1e-15:
        pz = pr
    r = -dpz / pz / pr * al
    rpz = pz / pr
    cod[0] += r * pxi
    cod[2] += r * pyi
    cod[4] -= (2.0 + rpz) * a / 2.0 / (pr + pz) * r
    return cod


def _tsolque_body_cod(cod: np.ndarray, aln1: float, akkp: float,
                      dx0: float = 0.0, dy0: float = 0.0) -> np.ndarray:
    """SADのtsolque ボディ1スライス COD更新 (bz=0, tsolque.f直訳).
    akkp = |K1|/pr > 0 必須。dx0 = ak0x/ak, dy0 = ak0y/ak (MULT の K0 による軸のずれ)。
    x方向: cos/sin (集束), y方向: cosh/sinh (発散)。"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    w1 = np.sqrt(akkp)
    phi1 = np.clip(aln1 * w1, -100.0, 100.0)
    c1 = np.cos(phi1); s1 = np.sin(phi1)
    # SAD の xsincos/xsincosh: xs = x − sin x, xsh = x − sinh x (tfloor.f)。符号を逆にすると z の x²/y² 項が半分になる
    dc1 = c1 - 1.0; xs1 = phi1 - s1
    ch2 = np.cosh(phi1); sh2 = np.sinh(phi1)
    dch2 = ch2 - 1.0; xsh2 = phi1 - sh2
    xi = cod[0] + dx0; yi = cod[2] + dy0
    pxi = cod[1] / pr; pyi = cod[3] / pr
    u1 = xi * dc1 + pxi * s1 / w1
    u2 = -xi * w1 * s1 + pxi * dc1
    v1 = yi * dch2 + pyi * sh2 / w1
    v2 = yi * w1 * sh2 + pyi * dch2
    pxf = pxi + u2
    pyf = pyi + v2
    cod[0] = cod[0] + u1    # = xi0 + u1 (SAD: cod(1)=xi0+u1)
    cod[1] = pxf * pr
    cod[2] = cod[2] + v1
    cod[3] = pyf * pr
    # z補正 (tsolque.f lines 420-425)
    cod[4] -= 0.25 * (
        w1 * (xi**2 * xs1 - yi**2 * xsh2)
        + (pxi**2 + pyi**2) * aln1
        + u1 * pxf + xi * pxi * dc1
        + v1 * pyf + yi * pyi * dch2)
    return cod


def _tsolque_cod(cod: np.ndarray, L: float, akkp: float,
                 ak: float = None) -> np.ndarray:
    """SADのtsolque COD更新ループ (bz=0, 放射なし, tsolque.f直訳).
    akkp = |K1|/pr > 0。ak = 積分K1 (SADの引数ak; 分割数の決定に使用)。
    構造: tqente(L/2N) + [body+tqente(L/N)]*(N-1) + body + tqente(L/2N)
    ndiv = min(1000, 1+int(|al*hypot(ak,bz)|/eps)), eps=0.1 (tsolque.f)"""
    if ak is None:
        ak = akkp * (1.0 + cod[5]) * L
    ndiv = min(1000, 1 + int(abs(L * ak) / 0.1))
    aln1 = L / ndiv
    al1 = aln1 * 0.5   # 最初は半幅
    for _n in range(ndiv):
        cod = _tqente_cod(cod, al1)
        cod = _tsolque_body_cod(cod, aln1, akkp)
        al1 = aln1      # 2回目以降は全幅
    cod = _tqente_cod(cod, aln1 * 0.5)   # 最後の半ドリフト
    return cod


def _akang(ak_re: float, ak_im: float, al: float):
    """tfloor.f akang: 複素 K1 (K1 + i SK1) の主軸の角 theta1 と cr1 = exp(−i theta1)。
    実数で K1·L < 0 なら theta1 = π/2, cr1 = −i (四極を 90° 回して |K1| にする)。"""
    a = complex(ak_re, ak_im) if al == 0.0 else complex(ak_re, ak_im) * al
    if a.imag == 0.0:
        if a.real < 0.0:
            return 0.5 * math.pi, complex(0.0, -1.0)
        return 0.0, complex(1.0, 0.0)
    th = math.atan2(a.imag, a.real) * 0.5
    return th, complex(math.cos(th), -math.sin(th))


def _mult_params(el) -> dict:
    """MULT (ANGLE=0, BZ=0, SK なし) の tpara (tffs.f) + tmulte/tmulti 共通の前処理。
    akr[n] = K_n·cr1^(n+1) (回した座標系での複素係数)、theta2 = ROTATE + DROTATE + akang。
    ndiv = max over n≥2 of int(sqrt(ampmax^(n−1)/6/(n−1)!/eps·|K_n·L|))+1 (eps = 0.005·EPS)。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    al = el.L
    knl = {int(k): float(v) for k, v in (getattr(el, 'KnL', {}) or {}).items() if v}
    nmmax = max(knl) if knl else -1
    th1, cr1 = _akang(knl.get(1, 0.0), 0.0, al)
    theta2 = g('ROTATION') + g('DROTATE') + th1
    akr = []
    cr = cr1
    for n in range(nmmax + 1):
        akr.append(knl.get(n, 0.0) * cr)
        cr = cr * cr1
    eps = 0.005 * (g('EPS') or 1.0)
    ndiv = 1
    if al != 0.0:
        for n in range(2, nmmax + 1):
            ndiv = max(ndiv, int(math.sqrt(0.05**(n - 1) / 6.0 / math.factorial(n - 1) / eps
                                           * abs(knl.get(n, 0.0) * al))) + 1)
    return dict(al=al, knl=knl, nmmax=nmmax, cr1=cr1, theta2=theta2, akr=akr, ndiv=ndiv, eps=eps,
                dx=g('DX'), dy=g('DY'), fringe=g('DISFRIN') == 0.0, mfring=int(round(g('FRINGE'))),
                bfrm=g('K0FR') == 0.0)


def _mult_rot_in(c, dx, dy, ct, st):
    """tsolrote (bz=0, dz=chi=0) の入口: ずらす → theta2 で回す"""
    c = c.copy()
    c[0] -= dx
    c[2] -= dy
    x, px, y, py = c[0], c[1], c[2], c[3]
    c[0] = ct * x - st * y; c[2] = st * x + ct * y
    c[1] = ct * px - st * py; c[3] = st * px + ct * py
    return c


def _mult_rot_out(c, dx, dy, ct, st):
    c = c.copy()
    x, px, y, py = c[0], c[1], c[2], c[3]
    c[0] = ct * x + st * y; c[2] = -st * x + ct * y
    c[1] = ct * px + st * py; c[3] = -st * px + ct * py
    c[0] += dx
    c[2] += dy
    return c


def _tbfrme_cod(cod: np.ndarray, ak: complex, ent: bool) -> np.ndarray:
    """tbfrie.f tbfrme (fb1=0) の COD 更新: MULT の K0 (双極) のハードエッジフリンジ。
    ak = ak0n/al1 (長さあたりの複素 K0)。akx = Re, aky = −Im。"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    akx = ak.real
    aky = -ak.imag
    xi, yi = cod[0], cod[2]
    pxi = min(pr, max(-pr, cod[1]))
    pyi = min(pr, max(-pr, cod[3]))
    y1 = -aky * xi + akx * yi
    px1 = akx * pxi + aky * pyi
    a = 1.0 / (akx**2 + aky**2) / pr
    dx = 0.5 * a * akx * y1**2
    dpx = a * aky * px1 * y1
    dy = 0.5 * a * aky * y1**2
    dpy = -a * akx * px1 * y1
    dz = -0.5 * a * px1 * y1**2 / pr
    cod[0] = xi + dx
    cod[1] = pxi + dpx
    cod[2] = yi + dy
    cod[3] = pyi + dpy
    cod[4] += dz
    return cod


def _asinz(x: float) -> float:
    """mathfun asinz"""
    if abs(x) < 0.0032:
        x2 = x * x
        return x + x * x2 * (1.0 / 6.0 + x2 * 3.0 / 40.0)
    return math.asin(max(-1.0, min(1.0, x)))


def _tsolconv(pz0, plz, pbz, bpr, phi):
    """tdrift.f tsolconv: phi の Newton 反復。戻りの sin 等は最後に評価した phi (更新前) の値 (SAD のまま)"""
    if plz == 0.0:
        sn, xs, cs, dc = _xsincos(phi)
        return phi, sn, xs, cs, dc
    for _ in range(15):
        sn, xs, cs, dc = _xsincos(phi)
        sv = plz * xs + pz0 * sn - pbz * dc
        u = -plz * dc + pz0 * cs + pbz * sn
        dphi = (bpr - sv) / u if u != 0.0 else (bpr - sv) / pz0
        phi0 = phi
        phi = phi + dphi
        if phi0 == phi or abs(dphi) <= 1e-15 * abs(phi):
            break
    return phi, sn, xs, cs, dc


def _tsoldz_cod(cod: np.ndarray, al: float, bxs: float, bys: float, dv: float) -> np.ndarray:
    """tdrift.f tsoldz (bzs=0, drift=.true.) の COD 更新: 一定の双極キック付きドリフト。
    tdrife からは bxs = ak0y/al, bys = ak0x/al で呼ばれる。cod の運動量は px/p0。"""
    cod = cod.copy()
    babs = math.hypot(bxs, bys)
    pr = 1.0 + cod[5]
    pxi, pyi = cod[1], cod[3]
    a = pxi**2 + pyi**2
    dpz0 = -a / pr / (1.0 + math.sqrt(max(0.0, 1.0 - a / pr**2)))
    pz0 = pr + dpz0
    if babs < 1e-20:
        r = al / pz0
        cod[0] += r * pxi
        cod[2] += r * pyi
        cod[4] += al * (dpz0 / pz0 - dv)
        return cod
    phix = bxs / babs
    phiy = bys / babs
    phiz = 0.0
    alb = 1.0 / babs
    albabs = al * babs
    dphizsq = phix**2 + phiy**2
    dpl = pxi * phix + pyi * phiy + dpz0 * phiz
    pl = pr * phiz + dpl
    plx, ply, plz = pl * phix, pl * phiy, pl * phiz
    ptx = pxi - plx
    pty = pyi - ply
    ptz = dpz0 - dpl * phiz + pr * dphizsq
    pbx = pty * phiz - ptz * phiy
    pby = ptz * phix - ptx * phiz
    pbz = ptx * phiy - pty * phix
    db = albabs - pbz
    ap = math.hypot(pz0, pbz)
    dphi = -math.atan2(pbz, pz0)
    phi = _asinz(db / ap) - dphi
    if al != 0.0:
        phi, sinphi, xsinphi, cosphi, dcosphi = _tsolconv(pz0, plz, pbz, albabs, phi)
    else:
        phi, sinphi, xsinphi, cosphi, dcosphi = 0.0, 0.0, 0.0, 1.0, 0.0
    cod[0] += (plx * phi + ptx * sinphi - pbx * dcosphi) * alb
    cod[2] += (ply * phi + pty * sinphi - pby * dcosphi) * alb
    cod[1] = pxi + ptx * dcosphi + pbx * sinphi
    cod[3] = pyi + pty * dcosphi + pby * sinphi
    cod[4] += ((dpl * phiz - dphizsq * pr) * xsinphi + dpz0 * sinphi - pbz * dcosphi) * alb - dv * al
    return cod


def _mult_sequence(el, dvf, ndiv_override: int = None) -> list:
    """tmulte.f tmulte (ANGLE=0, BZ=0, RF なし, F1 なし, SK なし) を写像の列に分解する。
    戻り値: [('map', f), ('rad', al), ...]。'rad' は放射ありのときに tradke(al, 0) を入れる位置
    (CALC では無視)。dvf(δ) は tgetdv (CALC では 0)。ndiv_override は放射ありの分割数 (itgetqraddiv)。
    順序 (tmulte.f): tsolrote → [tqfrie(ak1)] [tbfrme(ak0)] → [tsolque(半) → tmulke キック] ×ndiv
    (2 回目以降は全幅) → tsolque(半) → [tbfrme(−ak0)] [tqfrie(−ak1)] → tsolrote 戻し。"""
    pm = _mult_params(el)
    al, nmmax, akr, ndiv = pm['al'], pm['nmmax'], pm['akr'], pm['ndiv']
    if ndiv_override is not None:
        ndiv = min(300, max(ndiv, ndiv_override))
    ndiv = min(ndiv, 300)
    ct, st = math.cos(pm['theta2']), math.sin(pm['theta2'])
    dx, dy = pm['dx'], pm['dy']
    r_in = lambda c: _mult_rot_in(c, dx, dy, ct, st)
    r_out = lambda c: _mult_rot_out(c, dx, dy, ct, st)
    seq = []
    if nmmax < 0:
        # tmulte.f:130-135: K が全部 0 (VOLT も 0) なら tdrife0(al) = 長さ al の厳密ドリフト (tqente は tsolque の中の非線形部分だけ)
        seq.append(('map', lambda c: _exact_drift_cod(c, al) if dvf(0.0) == 0.0 else _dv_drift(c, al, dvf)))
        return seq
    aln = al / ndiv
    akn = [a / ndiv for a in akr]
    akn1 = akn[1].real if nmmax >= 1 else 0.0
    akn0 = akn[0]
    nz = al != 0.0
    fringe, mfring, bfrm = pm['fringe'], pm['mfring'], pm['bfrm']
    ak1h, al1h, ak0h = akn1 * 0.5, aln * 0.5, akn0 * 0.5

    def wrap(f):
        return lambda c, f=f: r_out(f(r_in(c)))
    if nz:
        if fringe and mfring != 2 and ak1h != 0.0:
            seq.append(('map', wrap(lambda c: _tqfrie_cod(c, ak1h, al1h))))
        if bfrm and ak0h != 0.0:
            seq.append(('map', wrap(lambda c: _tbfrme_cod(c, ak0h / al1h, True))))

    def tsolque_pieces(al1, ak1, ak0):
        """tsolque(al1, ak1, 0, Re ak0, Im ak0): 部分スライスごとの写像と放射点"""
        out = []
        if ak1 == 0.0:
            out.append(('map', wrap(lambda c: _tsoldz_cod(c, al1, ak0.imag / al1, ak0.real / al1, dvf(c[5])))))
            out.append(('rad', al1))
            return out
        ndiv_i = min(1000, 1 + int(abs(al1 * ak1) / 0.1))
        aln1 = al1 / ndiv_i
        dx0 = ak0.real / ak1
        dy0 = ak0.imag / ak1
        a1 = aln1 * 0.5
        for n in range(ndiv_i):
            def body(c, a1=a1):
                dv = dvf(c[5])
                c = _tqente_cod(c, a1)
                c = _tsolque_body_cod(c, aln1, ak1 / al1 / (1.0 + c[5]), dx0, dy0)
                c = c.copy()
                c[4] -= dv * aln1
                return c
            out.append(('map', wrap(body)))
            if n != ndiv_i - 1:
                out.append(('rad', aln1))
            a1 = aln1
        out.append(('map', wrap(lambda c: _tqente_cod(c, aln1 * 0.5))))
        out.append(('rad', aln1))
        return out

    def kick(c):
        c = c.copy()
        cx0 = complex(c[0], c[2])
        cx = 0j
        if nz:
            for kord in range(nmmax, 1, -1):
                cx = (cx + akn[kord]) * cx0 / kord
            cx = cx * cx0
        else:
            for kord in range(nmmax, 0, -1):
                cx = (cx + akn[kord]) * cx0 / kord
            cx = cx + akn[0]
        re = max(-0.9999, min(0.9999, cx.real))
        im = max(-0.9999, min(0.9999, cx.imag))
        c[1] -= re
        c[3] += im
        return c
    al1, ak1, ak0 = al1h, ak1h, ak0h
    for m in range(ndiv):
        if nz:
            seq += tsolque_pieces(al1, ak1, ak0)
            al1, ak1, ak0 = aln, akn1, akn0
        seq.append(('map', wrap(kick)))
    if nz:
        seq += tsolque_pieces(al1 * 0.5, ak1 * 0.5, ak0 * 0.5)
        if bfrm and ak0 != 0.0:
            seq.append(('map', wrap(lambda c: _tbfrme_cod(c, -ak0 / al1, False))))
        if fringe and mfring != 1 and ak1 != 0.0:
            seq.append(('map', wrap(lambda c: _tqfrie_cod(c, -ak1, al1))))
    return seq


def _dv_drift(c, al, dvf):
    c = _exact_drift_cod(c, al)
    c = c.copy()
    c[4] -= dvf(c[5]) * al
    return c


def _mult_cod_map(el):
    """CALC (qtwiss) 用: MULT の非線形 COD 写像 (放射なし, dv なし)"""
    seq = [f for kind, f in _mult_sequence(el, lambda d: 0.0) if kind == 'map']

    def f(c):
        for g in seq:
            c = g(c)
        return c
    return f


def _tqfrie_cod(cod: np.ndarray, ak: float, al: float) -> np.ndarray:
    """SADのtqfrie COD更新のみ (bz=0, tqfrie.f直訳)。
    四極のハードエッジフリンジ。ak = 積分K1 (入口 +ak1, 出口 -ak1), al = 長さ。"""
    cod = cod.copy()
    xmax = 1e10
    akk = ak / al / 4.0
    pr = 1.0 + cod[5]
    xi = min(xmax, max(-xmax, cod[0]))
    yi = min(xmax, max(-xmax, cod[2]))
    aki = akk / pr
    a = aki * xi**2
    b = aki * yi**2
    d = a + b
    ab = aki * (xi - yi) * (xi + yi)
    t = ab**2 / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    cod[0] = xi + dx1
    cod[2] = yi + dy1
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h**2
    if f == 0.0:
        return cod
    pxi = cod[1]
    pyi = cod[3]
    px1 = (pxi * yy + pyi * h) / f
    py1 = (pyi * xx - pxi * h) / f
    cod[1] = px1
    cod[3] = py1
    dz1 = -(cod[1] * (dx1 + t * xi) + cod[3] * (dy1 + t * yi)) / pr
    cod[4] = cod[4] + dz1
    return cod


def _tqlfre_cod(cod: np.ndarray, al: float, ak: float, f1: float, f2: float) -> np.ndarray:
    """SAD tqlfre.f の COD 更新のみ (bz=0)。四極の線形フリンジ。
    af1 = −ak/al·f1·|f1|/24, af2 = ak/al·f2, a = af1/p, ea = exp(a), bp = af2/p²
    x → ea·x + bp·px, px → px/ea, y → y/ea − bp·py, py → py·ea。z の式は SAD のまま (コード転記)。"""
    cod = cod.copy()
    af1 = -ak / al * f1 * abs(f1) / 24.0
    af2 = ak / al * f2
    p = 1.0 + cod[5]
    a = af1 / p
    b = af2 / p
    ea = math.exp(a)
    bp = b / p
    xf = ea * cod[0] + bp * cod[1]
    yf = cod[2] / ea - bp * cod[3]
    pxf = cod[1] / ea
    pyf = cod[3] * ea
    cod[4] = cod[4] - ((a * ea * cod[0] + bp * ea * (1.0 + 0.5 * a) * pxf) * pxf
                       - (a * cod[2] / ea + bp / ea * (1.0 - 0.5 * a) * pyf) * pyf) / p
    cod[1] = pxf
    cod[3] = pyf
    cod[0] = xf
    cod[2] = yf
    return cod


def _tquade_cod(cod: np.ndarray, L: float, K1: float,
                fringe: bool = True, f1in: float = 0.0, f2in: float = 0.0,
                f1out: float = 0.0, f2out: float = 0.0, mfring: int = 0) -> np.ndarray:
    """SADのtquade COD更新 (bz=0, dx=dy=0, 放射なし, tquade.f直訳).
    K1<0の場合: 座標を-90°回転 → tsolque(|K1|) → +90°回転 (tquade.f lines 30-36)。
    fringe=True (SADデフォルト DISFRIN=0, FRINGE=0): 両端で tqfrie を適用。"""
    if abs(K1) < 1e-10 or abs(L) < 1e-15:
        # K1=0: 純ドリフト (tquade.f line 26-28: tdrife0(al))。tqente は tsolque の中の非線形部分だけで、ドリフトの線形項を含まない
        return _exact_drift_cod(cod, L)
    pr = 1.0 + cod[5]
    ak1 = abs(K1) * L   # SAD の ak1 (積分K1, 回転後は正)

    def _body(c, akkp):
        # tquade.f: tqfrie (fringe かつ mfring>=0, 入口 /=2, 出口 /=1)、tqlfre (入口 mfring=1,3 / 出口 2,3)
        if fringe and mfring >= 0 and mfring != 2:
            c = _tqfrie_cod(c, ak1, L)
        if mfring in (1, 3):
            c = _tqlfre_cod(c, L, ak1, f1in, f2in)
        c = _tsolque_cod(c, L, akkp, ak1)
        if mfring in (2, 3):
            c = _tqlfre_cod(c, L, ak1, -f1out, f2out)
        if fringe and mfring >= 0 and mfring != 1:
            c = _tqfrie_cod(c, -ak1, L)
        return c

    if K1 * L < 0.0:
        # -90°回転: (x,px,y,py) → (y,py,-x,-px)  (tquade.f lines 30-36)
        cod = cod.copy()
        x, px, y, py = cod[0], cod[1], cod[2], cod[3]
        cod[0] = y; cod[1] = py; cod[2] = -x; cod[3] = -px
        akkp = (-K1) / pr
        cod = _body(cod, akkp)
        # +90°回転: (x,px,y,py) → (-y,-py,x,px)
        x, px, y, py = cod[0], cod[1], cod[2], cod[3]
        cod = cod.copy()
        cod[0] = -y; cod[1] = -py; cod[2] = x; cod[3] = px
    else:
        akkp = K1 / pr
        cod = _body(cod, akkp)
    return cod

def _drift_matrix_at_cod(L: float, delta: float,
                         pxi: float, pyi: float) -> np.ndarray:
    """
    SADの qtwiss.f DRIFT case の忠実実装。

    COD値 (pxi, pyi) 周りで線形化したDRIFT転送行列(6×6)を返す。
    a15 = -pr*pxi*alz, a35 = -pr*pyi*alz の項がCOD依存の核心。

    SADソース (qtwiss.f lines ~205-240):
        pr=1+cod(6), pxi=cod(2), pyi=cod(4)
        dpz=-a/pr/(1+sqrt(max(0,1-a/pr^2)))  where a=pxi^2+pyi^2
        pzi=pr+dpz, ale=L/pzi, alz=ale/pzi^2
        a12=ale+pxi^2*alz, a14=pxi*pyi*alz, a15=-pr*pxi*alz
        a32=a14, a34=ale+pyi^2*alz, a35=-pr*pyi*alz
    """
    pr = 1.0 + delta
    a = pxi**2 + pyi**2
    if a > 0.0 and abs(pr) > 1e-15:
        arg = max(0.0, 1.0 - a / pr**2)
        dpz = -a / pr / (1.0 + np.sqrt(arg))
    else:
        dpz = 0.0
    pzi = pr + dpz
    if abs(pzi) < 1e-15:
        pzi = pr
    ale = L / pzi
    alz = ale / pzi**2

    M = np.eye(6)
    M[0, 1] = ale + pxi**2 * alz   # a12
    M[0, 3] = pxi * pyi * alz       # a14
    M[0, 5] = -pr * pxi * alz       # a15  COD依存キー項
    M[2, 1] = pxi * pyi * alz       # a32 = a14
    M[2, 3] = ale + pyi**2 * alz    # a34
    M[2, 5] = -pr * pyi * alz       # a35  COD依存キー項
    return M


# 薄レンズ多極のクラス名 → (強度属性名, kord)。SAD tthine の kord = nord/2-1。
# RFSW フラグ (SAD の FFS フラグ。既定 OFF)。OFF のとき CALC の空洞はドリフト。ffs 層で切り替える
RFSW = False

_THIN_KORD = {'Sextupole': ('K2', 2), 'Octupole': ('K3', 3),
              'Decapole': ('K4', 4), 'Dodecapole': ('K5', 5)}


def _exact_drift_cod(c: np.ndarray, L: float) -> np.ndarray:
    """厳密ドリフトの COD 更新 (SAD tdrife1 / qtwiss DRIFT と同じ式, dv 項なし)。"""
    c = c.copy()
    pr = 1.0 + c[5]
    a = c[1]**2 + c[3]**2
    pzi = pr * np.sqrt(max(0.0, 1.0 - a / pr**2))
    c[0] += c[1] / pzi * L
    c[2] += c[3] / pzi * L
    c[4] -= a / (pr + pzi) / pzi * L
    return c


def _tthine_cod(c: np.ndarray, L: float, ak: float, kord: int) -> np.ndarray:
    """
    SAD tthine.f の COD 更新のみ (放射なし, dx=dy=theta=0)。
    薄レンズ多極 (SEXT kord=2, OCT 3, DECA 4, DODECA 5)。
    ak = 積分強度 (SAD の K)。aki = ak/kord!。
    長さありのとき: ドリフト L/6 → キック(aki/2) → ドリフト 2L/3 → キック(aki/2) → ドリフト L/6
    キック: cx = (x − i y)^kord,  px −= aki·Re(cx),  py −= aki·Im(cx)  (正準運動量、1/(1+δ) なし)
    """
    if ak == 0.0:
        return _exact_drift_cod(c, L)
    aki = ak / math.factorial(kord)

    def kick(c, a):
        c = c.copy()
        cx = complex(c[0], -c[2]) ** kord
        c[1] -= a * cx.real
        c[3] -= a * cx.imag
        return c

    if L == 0.0:
        return kick(c, aki)
    c = _exact_drift_cod(c, L / 6.0)
    c = kick(c, aki * 0.5)
    c = _exact_drift_cod(c, L / 1.5)
    c = kick(c, aki * 0.5)
    return _exact_drift_cod(c, L / 6.0)


def _thin_offset_map(L: float, ak: float, kord: int,
                     dx: float = 0.0, dy: float = 0.0, theta: float = 0.0):
    """
    tthine の COD 写像に SAD tchge のずれ・回転を付けたもの (alg=phig=dz=dtheta=dchi2=0)。
    入口: (x,px),(y,py) を θ 回転 → x −= cosθ·dx − sinθ·dy, y −= sinθ·dx + cosθ·dy
    出口: 逆 (tchge.f, tbrote)。θ の向きの SAD との照合は回転ありのベンチで行う。
    """
    if dx == 0.0 and dy == 0.0 and theta == 0.0:
        return lambda c: _tthine_cod(c, L, ak, kord)
    ct, st = math.cos(theta), math.sin(theta)
    sx = ct * dx - st * dy
    sy = st * dx + ct * dy

    def f(c):
        c = c.copy()
        x, y, px, py = c[0], c[2], c[1], c[3]
        c[0] = ct * x - st * y - sx
        c[2] = st * x + ct * y - sy
        c[1] = ct * px - st * py
        c[3] = st * px + ct * py
        c = _tthine_cod(c, L, ak, kord)
        x, y, px, py = c[0] + sx, c[2] + sy, c[1], c[3]
        c = c.copy()
        c[0] = ct * x + st * y
        c[2] = -st * x + ct * y
        c[1] = ct * px + st * py
        c[3] = -st * px + ct * py
        return c
    return f


def _thin_cod_matrix(L: float, ak: float, kord: int, dx: float, dy: float, th: float, cod: np.ndarray):
    """薄い多極 (tthine) の (出口 COD, 転送行列)。ずれ・ロールは tthine の先頭と末尾の tchge (入口: 回す→ずらす、出口: 逆)。
    行列は J_出口·M·J_入口 (入口は SAD の癖で th≠0 かつずれ≠0 のとき回転が入らない: _tchge_jac)。ak=0 は tdrife0 だけ (tchge を呼ばない)。"""
    _f = _thin_offset_map(L, ak, kord)

    def _jm(c):
        # tthine の解析行列 (ドリフト L/6 → キック → 2L/3 → キック → L/6 の各段の微分の積)。pz が潰れるときだけ数値微分
        from .emit6d import _jac_thin
        J = _jac_thin(c, L, ak, kord, 0.0)
        return J if J is not None else _num_jacobian(_f, c)
    if ak == 0.0 or (dx == 0.0 and dy == 0.0 and th == 0.0):
        return _f(cod), _jm(cod)
    c1 = _tchge_cod(cod, dx, dy, 0.0, th, 0.0, 0.0, 0.0, True)
    J1 = _tchge_jac(cod, dx, dy, 0.0, th, 0.0, 0.0, 0.0, True)
    c2 = _f(c1)
    M = _jm(c1)
    c3 = _tchge_cod(c2, dx, dy, 0.0, th, 0.0, 0.0, 0.0, False)
    return c3, _tchge_jac(c2, dx, dy, 0.0, th, 0.0, 0.0, 0.0, False) @ M @ J1


def _xsincos(x):
    """mathfun xsincos: (sin, x−sin, cos, cos−1)"""
    sn = math.sin(x)
    return sn, x - sn, math.cos(x), -2.0 * math.sin(0.5 * x)**2


def _tbgrote_cod(cod, dr, dx, dz):
    """SAD tchge.f tbgrote の COD 更新のみ (dr = 回転行列 − I)。"""
    cod = cod.copy()
    xi = cod[0] + dx
    yi = cod[2]
    pxi = cod[1]
    pyi = cod[3]
    pr = 1.0 + cod[5]
    a = (pxi / pr)**2 + (pyi / pr)**2
    dpzi = -a / (1.0 + math.sqrt(max(1e-200, 1.0 - a))) * pr
    pzi = pr + dpzi
    zi = dz
    xf = dr[0][0] * xi + dr[0][1] * yi + dr[0][2] * zi + cod[0]
    yf = dr[1][0] * xi + dr[1][1] * yi + dr[1][2] * zi + yi
    zf = dr[2][0] * xi + dr[2][1] * yi + dr[2][2] * zi
    r11 = 1.0 + dr[0][0]
    r22 = 1.0 + dr[1][1]
    r33 = 1.0 + dr[2][2]
    pxf = r11 * pxi + dr[0][1] * pyi + dr[0][2] * pzi
    pyf = dr[1][0] * pxi + r22 * pyi + dr[1][2] * pzi
    dpzf = dr[2][0] * pxi + dr[2][1] * pyi + r33 * dpzi + dr[2][2] * pr
    pzf = pr + dpzf
    cod[1] = pxf
    cod[3] = pyf
    cod[0] = xf - pxf * zf / pzf
    cod[2] = yf - pyf * zf / pzf
    cod[4] = cod[4] + pr * zf / pzf
    return cod


def _tbrote_cod(cod, alg, phig, dx, dz, dtheta, dchi2, ent):
    """SAD tchge.f tbrote の COD 更新のみ (dz≠0 のドリフトは未対応)。"""
    if phig == 0.0:
        c0, s0, dl = 1.0, 0.0, alg
    else:
        c0, s0 = math.cos(phig), math.sin(phig)
        dl = alg * s0 / phig
    s2, _, c2, dc2 = _xsincos(dchi2)
    st, _, ct, dct = _xsincos(dtheta)
    rr = [[0.0] * 3 for _ in range(3)]
    rr[0][0] = dct * c0**2 + dc2 * s0**2 + c0 * s0 * s2 * st
    rr[1][1] = c2 * dct + dc2
    rr[2][2] = dc2 * c0**2 + dct * s0**2 - c0 * s0 * s2 * st
    cod = cod.copy()
    if ent:
        dx1 = dx * c0 - dz * s0
        dz1 = dz * c0 + dx * s0
        rr[0][1] = s0 * s2 - c0 * c2 * st
        rr[0][2] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[1][0] = -ct * s0 * s2 + c0 * st
        rr[1][2] = c0 * ct * s2 + s0 * st
        rr[2][0] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[2][1] = -c0 * s2 - c2 * s0 * st
        cod[0] -= dx1
        if dz1 != 0.0:
            cod = _tdrife1_cod(cod, 1.0, dz1)      # tchge.f: call tdrife1(trans,cod,1.d0,0.d0,dz1)
        return _tbgrote_cod(cod, rr, dl * s0, -dl * c0)
    rr[1][0] = s0 * s2 - c0 * c2 * st
    rr[2][0] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
    rr[0][1] = -ct * s0 * s2 + c0 * st
    rr[2][1] = c0 * ct * s2 + s0 * st
    rr[0][2] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
    rr[1][2] = -c0 * s2 - c2 * s0 * st
    cod = _tbgrote_cod(cod, rr, dl * s0, -dl * c0)
    dx1 = dx * c0 + dz * s0
    dz1 = dz * c0 - dx * s0
    if dz1 != 0.0:
        cod = _tdrife1_cod(cod, 1.0, -dz1)         # tchge.f: call tdrife1(trans2,cod,1.d0,0.d0,-dz1)
    cod = cod.copy()
    cod[0] += dx1
    return cod


def _tdrife1_cod(cod, dv, al):
    """toplvl.f tdrife1 の COD 更新 (引数 dv はそのまま。tbrote からは dv=1 で呼ばれる)"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    pxi, pyi = cod[1], cod[3]
    a = pxi**2 + pyi**2
    pzi = pr * math.sqrt(1.0 - a / pr**2)
    cod[0] += pxi / pzi * al
    cod[2] += pyi / pzi * al
    cod[4] -= (a / (pr + pzi) / pzi + dv) * al
    return cod


def _bend_theta(el) -> float:
    """BEND の ROT (tbende の theta) = 設計ロール ROTATE + 据え付けロール ROTATION (旧 Python API 用)。SAD の ROT は 1 つだけ。"""
    return float(getattr(el, 'ROTATE', 0.0) or 0.0) + float(getattr(el, 'ROTATION', 0.0) or 0.0)


def _bend_chg_args(el):
    """BEND の tchge の引数 (dx, dy, dz, theta, dtheta, alg, phig)。要らなければ None (tchge を呼ばない/恒等)。
    tbende.f の分岐に合わせる:
      ANGLE≠0 (一般・phib=0・al=0 すべて) : tchge(dx,dy,dz,theta,dtheta,dchi2, alg=L/2, phig=ANGLE/2)  (th=theta。dtheta は tbrote の BEND 中点まわり)
      ANGLE=0, K1=0 : tsteee。K0≠0 のときだけ tchge(dx,dy,0,theta,0,...)。K0=0 はドリフト (ずれ・回転とも効かない)。DROTATE は渡らない
      ANGLE=0, K1≠0, K0=0 : tquade(dx,dy,theta) (DROTATE は渡らない)
      ANGLE=0, K1≠0, K0≠0 : tmulte (theta2 = theta + dtheta + akang → ここでは theta+dtheta をまとめて回す)
    CHI2 は未対応 (parser が読まない)。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    dx, dy, dz, dth, th = g('DX'), g('DY'), g('DZ'), g('DROTATE'), _bend_theta(el)
    if g('ANGLE') != 0.0:
        if dx == 0.0 and dy == 0.0 and dz == 0.0 and dth == 0.0 and th == 0.0:
            return None
        return (dx, dy, dz, th, dth, 0.5 * el.L, 0.5 * g('ANGLE'))
    k0, k1 = g('K0'), g('K1')
    if k1 == 0.0:
        if k0 == 0.0:
            return None
        dth = 0.0
    elif k0 == 0.0:
        dth = 0.0
    else:
        th, dth = th + dth, 0.0
    if dx == 0.0 and dy == 0.0 and th == 0.0:
        return None
    return (dx, dy, 0.0, th, 0.0, 0.0, 0.0)


def _bend_noerr(el):
    """ずれ・ロールを 0 にした BEND のコピー (tchge の外側で扱った後に中身へ渡す)。"""
    import copy
    e = copy.copy(el)
    for k in ('DX', 'DY', 'DZ', 'DROTATE', 'ROTATE', 'ROTATION'):
        if hasattr(e, k):
            setattr(e, k, 0.0)
    return e


def _tchge_cod(cod, dx0, dy, dz, theta, dtheta, alg, phig, ent):
    """SAD tchge.f の COD 更新 (dchi2=0)。ent=True 入口、False 出口 (出口は alg, phig に alg−al0, phig−phi0 を渡す)。
    th = theta (phig≠0) / theta+dtheta (phig=0)。入口: theta で回す (x'=c x − s y, y'=s x + c y。px, py も) → 回した座標系でずらす
    (x −= c dx0 − s dy, y −= s dx0 + c dy) → tbrote。出口はその逆。"""
    th = theta if phig != 0.0 else theta + dtheta
    c = np.array(cod, dtype=float)
    if ent:
        if th != 0.0:
            ct, st = math.cos(th), math.sin(th)
            dx = ct * dx0 - st * dy
            xi, pxi = c[0], c[1]
            c[0] = ct * xi - st * c[2]
            c[2] = st * xi + ct * c[2]
            c[1] = ct * pxi - st * c[3]
            c[3] = st * pxi + ct * c[3]
            c[2] -= st * dx0 + ct * dy
        else:
            dx = dx0
            c[2] -= dy
        if dx != 0.0 or dz != 0.0 or dtheta != 0.0:
            c = _tbrote_cod(c, alg, phig, dx, dz, dtheta, 0.0, True)
        return c
    if th != 0.0:
        ct, st = math.cos(th), math.sin(th)
        dx = ct * dx0 - st * dy
    else:
        ct, st = 1.0, 0.0
        dx = dx0
    if dx != 0.0 or dz != 0.0 or dtheta != 0.0:
        c = _tbrote_cod(c, alg, phig, dx, dz, dtheta, 0.0, False)
    if th != 0.0:
        c[2] += st * dx0 + ct * dy
        xi, pxi = c[0], c[1]
        c[0] = ct * xi + st * c[2]
        c[2] = -st * xi + ct * c[2]
        c[1] = ct * pxi + st * c[3]
        c[3] = -st * pxi + ct * c[3]
    else:
        c[2] += dy
    return c


def _tdrife1_mat(cod, dv, dvdp, al):
    """toplvl.f tdrife1 の trans1 (cod の更新は _tdrife1_cod)"""
    pr = 1.0 + cod[5]
    pxi, pyi = cod[1], cod[3]
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
    t[4, 5] = dvdp * al + a * alz
    return t


def _tbgrote_mat(cod, dr, dx, dz):
    """tchge.f tbgrote の (出口 cod, trans)。cod の更新は _tbgrote_cod と同じ式。
    trans は x,px,y,py,z の行が写像の厳密な微分で、(5,6) だけ SAD の式 r33·(dpzf·pr/pzf + dpzi/pzi) − dr33
    (写像の微分 zf/pzf·(1 − pr/pzf·r33·dpzidp) ではない。dr=0 でも dpzi·(pr+1)/pzi が残る) をそのまま使う。"""
    xi = cod[0] + dx
    yi = cod[2]
    pxi = cod[1]
    pyi = cod[3]
    pr = 1.0 + cod[5]
    a = (pxi / pr)**2 + (pyi / pr)**2
    dpzi = -a / (1.0 + math.sqrt(max(1e-200, 1.0 - a))) * pr
    pzi = pr + dpzi
    zi = dz
    xf = dr[0][0] * xi + dr[0][1] * yi + dr[0][2] * zi + cod[0]
    yf = dr[1][0] * xi + dr[1][1] * yi + dr[1][2] * zi + yi
    zf = dr[2][0] * xi + dr[2][1] * yi + dr[2][2] * zi
    r11 = 1.0 + dr[0][0]
    r22 = 1.0 + dr[1][1]
    r33 = 1.0 + dr[2][2]
    pxf = r11 * pxi + dr[0][1] * pyi + dr[0][2] * pzi
    pyf = dr[1][0] * pxi + r22 * pyi + dr[1][2] * pzi
    dpzf = dr[2][0] * pxi + dr[2][1] * pyi + r33 * dpzi + dr[2][2] * pr
    pzf = pr + dpzf
    pxfpzf = pxf / pzf
    pyfpzf = pyf / pzf
    zfpzf = zf / pzf
    prpzf = pr / pzf
    c = cod.copy()
    c[1] = pxf
    c[3] = pyf
    c[0] = xf - pxf * zfpzf
    c[2] = yf - pyf * zfpzf
    c[4] = cod[4] + pr * zfpzf
    dpzidpx = -pxi / pzi
    dpzidpy = -pyi / pzi
    dpzidp = pr / pzi
    t = np.zeros((6, 6))
    t[1, 1] = r11 + dr[0][2] * dpzidpx
    t[1, 3] = dr[0][1] + dr[0][2] * dpzidpy
    t[1, 5] = dr[0][2] * dpzidp
    t[3, 1] = dr[1][0] + dr[1][2] * dpzidpx
    t[3, 3] = r22 + dr[1][2] * dpzidpy
    t[3, 5] = dr[1][2] * dpzidp
    t[0, 0] = r11 - dr[2][0] * pxfpzf
    t[0, 1] = -zfpzf * (t[1, 1] - pxfpzf * (dr[2][0] + r33 * dpzidpx))
    t[0, 2] = dr[0][1] - dr[2][1] * pxfpzf
    t[0, 3] = -zfpzf * (t[1, 3] - pxfpzf * (dr[2][1] + r33 * dpzidpy))
    t[0, 5] = -zfpzf * (dr[0][2] - r33 * pxfpzf) * dpzidp
    t[2, 0] = dr[1][0] - dr[2][0] * pyfpzf
    t[2, 1] = -zfpzf * (t[3, 1] - pyfpzf * (dr[2][0] + r33 * dpzidpx))
    t[2, 2] = r22 - dr[2][1] * pyfpzf
    t[2, 3] = -zfpzf * (t[3, 3] - pyfpzf * (dr[2][1] + r33 * dpzidpy))
    t[2, 5] = -zfpzf * (dr[1][2] - r33 * pyfpzf) * dpzidp
    t[4, 0] = dr[2][0] * prpzf
    t[4, 1] = -zfpzf * prpzf * (dr[2][0] + r33 * dpzidpx)
    t[4, 2] = dr[2][1] * prpzf
    t[4, 3] = -zfpzf * prpzf * (dr[2][1] + r33 * dpzidpy)
    t[4, 4] = 1.0
    t[4, 5] = r33 * (dpzf * prpzf + dpzi / pzi) - dr[2][2]
    t[5, 5] = 1.0
    return c, t


def _tbrote_jac(cod, alg, phig, dx, dz, dtheta, dchi2, ent):
    """tchge.f tbrote の (出口 cod, trans)。cod の更新は _tbrote_cod と同じ。
    入口: x −= dx1 → (dz1≠0 なら tdrife1(dz1)) → tbgrote。出口: tbgrote → (dz1≠0 なら tdrife1(−dz1)) → x += dx1。
    行列は後の段を左から掛ける (tmultr5)。"""
    if phig == 0.0:
        c0, s0, dl = 1.0, 0.0, alg
    else:
        c0, s0 = math.cos(phig), math.sin(phig)
        dl = alg * s0 / phig
    s2, _, c2, dc2 = _xsincos(dchi2)
    st, _, ct, dct = _xsincos(dtheta)
    rr = [[0.0] * 3 for _ in range(3)]
    rr[0][0] = dct * c0**2 + dc2 * s0**2 + c0 * s0 * s2 * st
    rr[1][1] = c2 * dct + dc2
    rr[2][2] = dc2 * c0**2 + dct * s0**2 - c0 * s0 * s2 * st
    cod = cod.copy()
    if ent:
        dx1 = dx * c0 - dz * s0
        dz1 = dz * c0 + dx * s0
        rr[0][1] = s0 * s2 - c0 * c2 * st
        rr[0][2] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[1][0] = -ct * s0 * s2 + c0 * st
        rr[1][2] = c0 * ct * s2 + s0 * st
        rr[2][0] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[2][1] = -c0 * s2 - c2 * s0 * st
        cod[0] -= dx1
        if dz1 != 0.0:
            J1 = _tdrife1_mat(cod, 1.0, 0.0, dz1)
            cod = _tdrife1_cod(cod, 1.0, dz1)
            cod, J2 = _tbgrote_mat(cod, rr, dl * s0, -dl * c0)
            return cod, J2 @ J1
        return _tbgrote_mat(cod, rr, dl * s0, -dl * c0)
    rr[1][0] = s0 * s2 - c0 * c2 * st
    rr[2][0] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
    rr[0][1] = -ct * s0 * s2 + c0 * st
    rr[2][1] = c0 * ct * s2 + s0 * st
    rr[0][2] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
    rr[1][2] = -c0 * s2 - c2 * s0 * st
    cod, J = _tbgrote_mat(cod, rr, dl * s0, -dl * c0)
    dx1 = dx * c0 + dz * s0
    dz1 = dz * c0 - dx * s0
    if dz1 != 0.0:
        J2 = _tdrife1_mat(cod, 1.0, 0.0, -dz1)
        cod = _tdrife1_cod(cod, 1.0, -dz1)
        J = J2 @ J
    cod[0] += dx1
    return cod, J


def _tchge_jac(cod, dx0, dy, dz, theta, dtheta, alg, phig, ent):
    """_tchge_cod の cod での 6x6 Jacobian (SAD tchge.f の trans の作り方。11/05 から解析式: _tbrote_jac/_tbgrote_mat)。
    tbrote を呼ぶのは「回した座標系のずれ dxr = cosθ·dx0−sinθ·dy」「dz」「dtheta」のどれかが非 0 のとき。
    tbrote の trans は、dz1 = dz·cosφ + dxr·sinφ ≠ 0 なら先に tdrife1 (長さ dz1 のドリフト) が入り、dz=dtheta=0 かつ φ=0 (または dxr=0) なら単位行列。
    SAD の癖 (再現する): 入口で th≠0 かつ tbrote を呼ぶとき、回転行列 trans1 を渡された tbrote が行列を「上書き」するので、
    行列には回転が入らない (COD は回る)。出口は tmultr5 で掛けるので回転が入る → 一周行列は R^T·M (非対称) になり、
    結合 (R1..R4)・チューンが物理とずれる。(5,6) は tbgrote の式のまま (_tbgrote_mat 参照)。"""
    th = theta if phig != 0.0 else theta + dtheta
    if th != 0.0:
        ct, st = math.cos(th), math.sin(th)
        dxr = ct * dx0 - st * dy
    else:
        ct, st = 1.0, 0.0
        dxr = dx0
    brote = (dxr != 0.0 or dz != 0.0 or dtheta != 0.0)
    simple = (dz == 0.0 and dtheta == 0.0 and (phig == 0.0 or dxr == 0.0))   # tbrote の trans が単位行列 (tbgrote の (5,6) を除く)
    c = np.array(cod, dtype=float)
    if ent and th != 0.0 and brote:
        # 回して y をずらした後の COD で tbrote の行列だけ (回転は入らない)
        xi, pxi = c[0], c[1]
        c[0] = ct * xi - st * c[2]
        c[2] = st * xi + ct * c[2]
        c[1] = ct * pxi - st * c[3]
        c[3] = st * pxi + ct * c[3]
        c[2] -= st * dx0 + ct * dy
        if simple:
            # tbgrote (dr=0: 回転なし): trans は単位行列で、(5,6) だけ dpzi·(pr+1)/pzi が残る (SAD のまま)
            J = np.eye(6)
            pr = 1.0 + c[5]
            a = (c[1] / pr)**2 + (c[3] / pr)**2
            dpzi = (math.sqrt(1.0 - a) - 1.0) * pr
            J[4, 5] = dpzi * (pr + 1.0) / (pr + dpzi)
            return J
        return _tbrote_jac(c, alg, phig, dxr, dz, dtheta, 0.0, True)[1]
    if simple:
        J = np.eye(6)
        if th != 0.0:
            if ent:
                J[:4, :4] = [[ct, 0., -st, 0.], [0., ct, 0., -st], [st, 0., ct, 0.], [0., st, 0., ct]]
            else:
                J[:4, :4] = [[ct, 0., st, 0.], [0., ct, 0., st], [-st, 0., ct, 0.], [0., -st, 0., ct]]
        return J
    if ent:
        # th == 0 (th≠0 は上で返した): y をずらして tbrote
        c[2] -= dy
        return _tbrote_jac(c, alg, phig, dxr, dz, dtheta, 0.0, True)[1]
    # 出口: tbrote → 回転 (tmultr5 で掛かるので回転が入る)
    _, J = _tbrote_jac(c, alg, phig, dxr, dz, dtheta, 0.0, False)
    if th != 0.0:
        R = np.eye(6)
        R[:4, :4] = [[ct, 0., st, 0.], [0., ct, 0., st], [-st, 0., ct, 0.], [0., -st, 0., ct]]
        J = R @ J
    return J


def _quad_fringe_params(el) -> tuple:
    """tsetfringepe (tffs.f:982): f1in=F1+F1K1F, f1out=F1+F1K1B, f2in=F2+F2K1F, f2out=F2+F2K1B、
    mfring = FRINGE (FRMD)。逆向き要素 (ori=.false.) の入れ替えは未対応。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if el.L == 0.0 or g('K1') == 0.0:
        return (0.0, 0.0, 0.0, 0.0, 0)
    return (g('F1') + g('F1K1F'), g('F2') + g('F2K1F'),
            g('F1') + g('F1K1B'), g('F2') + g('F2K1B'), int(round(g('FRINGE'))))


def _quad_offset_map(L: float, K1: float, dx: float = 0.0, dy: float = 0.0,
                     theta: float = 0.0, fringe: bool = True, fr: tuple = (0.0, 0.0, 0.0, 0.0, 0)):
    """
    四極の COD 写像 (tquade)。SAD tsolrote (bz=0, dz=chi1=chi2=0) と同じく
    入口: x −= dx, y −= dy → ロール θ で回す (x' = cosθ·x − sinθ·y, y' = sinθ·x + cosθ·y)
    出口: 逆回転 → x += dx, y += dy。K1<0 の −π/2 回転は _tquade_cod の中 (θ−π/2 と同じ)。
    """
    if dx == 0.0 and dy == 0.0 and theta == 0.0:
        return lambda c: _tquade_cod(c, L, K1, fringe, *fr)
    ct, st = math.cos(theta), math.sin(theta)

    def f(c):
        c = c.copy()
        x, y = c[0] - dx, c[2] - dy
        px, py = c[1], c[3]
        c[0] = ct * x - st * y
        c[2] = st * x + ct * y
        c[1] = ct * px - st * py
        c[3] = st * px + ct * py
        c = _tquade_cod(c, L, K1, fringe, *fr)
        x, y, px, py = c[0], c[2], c[1], c[3]
        c = c.copy()
        c[0] = ct * x + st * y + dx
        c[2] = -st * x + ct * y + dy
        c[1] = ct * px + st * py
        c[3] = -st * px + ct * py
        return c
    return f

def _bend_k1_cod_matrix(el, cod: np.ndarray):
    """K1≠0 BEND の COD 更新と COD 位置での転送行列 (emit6d の tbende 直訳を使用)。
    運動量 (h0) は z 列にしか効かないので仮に 1 GeV/c を使う。
    据え付け誤差なしの BEND は、tbende の入口・出口エッジも SAD の解析行列 (tbfrie.f の trans) で組む
    (emit6d の _tbende はエッジを数値微分にしており、軌道が大きいと SAD と 1e-8 程度ずれる)。"""
    from .emit6d import _Phys, _State, _tbende, _tbendecorr, _tbendebody, EPSBEND
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    st = _State(cod, _Phys(1e9, 1.0))
    if (g('DX') != 0.0 or g('DY') != 0.0 or g('DROTATE') != 0.0 or g('ANGLE') == 0.0 or g('K1') == 0.0
            or g('F1') != 0.0 or g('FB1') != 0.0 or g('FB2') != 0.0):
        _tbende(st, el, enarad=False)
        return st.cod.copy(), st.T.copy()
    # tbende.f の K1≠0 の分岐 (emit6d._tbende と同じ順序。エッジだけ解析行列)
    al = el.L
    phi0 = el.ANGLE
    phib = phi0 + g('K0')
    ak = g('K1') * al
    rhob = al / phib
    rho0 = al / phi0
    ndiv = 1 + int(abs(phi0 / EPSBEND))
    aln = al / ndiv
    c, M = _bend_edge_chain(st.cod.copy(), np.eye(6), rhob, g('E1'), True)
    st.cod = c
    st.apply(M)
    alx0 = 0.0
    akx0 = 0.0
    for n in range(1, ndiv + 1):
        akx = ak * aln / al
        phi1 = phi0 * aln / al
        _tbendecorr(st, (akx + akx0) * 0.5, (aln + alx0) * 0.5, rho0)
        akx0 = akx
        alx0 = aln
        _tbendebody(st, aln, phi1, akx, aln, rhob, rho0, False)
    _tbendecorr(st, akx0 * 0.5, alx0 * 0.5, rho0)
    c, M = _bend_edge_chain(st.cod.copy(), np.eye(6), rhob, g('E2'), False)
    st.cod = c
    st.apply(M)
    return st.cod.copy(), st.T.copy()


def _bend_core_cod_matrix(el, cod: np.ndarray):
    """ずれ・ロールのない BEND の (出口 COD, COD での 6x6 転送行列)。tbende.f の分岐ごと:
    ANGLE=0 → tsteee (K1=0) / tmulte (K1≠0 かつ K0≠0) / tquade (K1≠0 かつ K0=0)、ANGLE≠0 → K1=0 は tbendebody0、K1≠0 は tbendebody。"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if float(el.ANGLE) == 0.0 and (abs(g('K1')) < 1e-10 or g('K0') != 0.0):
        # ANGLE=0 の BEND (ステアリング): tbende.f の phi0==0 の分岐。K1=0 は tsteee / 薄いキック / ドリフト、
        # K1≠0 かつ K0≠0 は tmulte (K0=0 の K1≠0 は下の K1≠0 の分岐で tquade)
        if abs(g('K1')) < 1e-10:
            return _tsteee_cod_matrix(el, cod)
        _f = _mult_cod_map(_bend_as_mult(el))
        return _f(cod), _num_jacobian(_f, cod)
    if abs(g('K1')) < 1e-10:
        # K1=0のBEND: SADのtbende相当の非線形COD追跡
        # tbedge(入口) → tbendebody0(ボディ) → tbedge(出口) (Fortranソース直訳: tbende.f, tbfrie.f)
        E1 = getattr(el, 'E1', 0.0)
        E2 = getattr(el, 'E2', 0.0)
        _pb = el.ANGLE + g('K0')
        _fp = _bend_fringe_params(el)
        _fb = lambda c: _bend_cod_update(c, el.L, el.ANGLE, E1, E2, _pb, *_fp)
        if _fp[0] == 0.0 and _fp[1] == 0.0:
            # F1 フリンジなし: SAD の解析転送行列 (tbedge/tbendebody0 の trans1)。
            # 写像の厳密な微分 (数値微分) とは 1e-7 程度ずれ、軌道が大きい INS で Twiss が SAD と 1e-6 違ってくる。
            cod_out, M_el = _bend_cod_matrix(cod, el.L, el.ANGLE, E1, E2, _pb)
            if M_el is None:
                M_el = _num_jacobian(_fb, cod)
            return cod_out, M_el
        # SADのqbend: tbende(transe,cod,...,mat=.true.) がCOD位置でJacobianを計算するのと等価 (数値微分)
        return _fb(cod), _num_jacobian(_fb, cod)
    # K1≠0 の BEND: SAD qbend → tbende(enarad=.false.) の直訳
    # (tbedge → [tbendecorr + tbendebody]×ndiv + tbendecorr → tbedge)
    return _bend_k1_cod_matrix(el, cod)


def _bend_element_cod_matrix(el, cod: np.ndarray, delta: float = 0.0):
    """BEND の (出口 COD, 転送行列)。ずれ DX/DY・ロール ROTATE/ROTATION・DROTATE は SAD tbende と同じく
    tchge (入口: 回す→ずらす→tbrote、出口はその逆) で外から付け、中身 (_bend_core_cod_matrix) は誤差なしで計算する。
    行列は J_出口 · M_中身 · J_入口 (J は tchge の解析 Jacobian。DROTATE≠0/DZ≠0 のときだけ数値微分)。"""
    chg = _bend_chg_args(el)
    if chg is None:
        return _bend_core_cod_matrix(el, cod)
    dx, dy, dz, th, dth, alg, phig = chg
    core = _bend_noerr(el)
    c1 = _tchge_cod(cod, dx, dy, dz, th, dth, alg, phig, True)
    J1 = _tchge_jac(cod, dx, dy, dz, th, dth, alg, phig, True)
    c2, M = _bend_core_cod_matrix(core, c1)
    # 出口は alg−al0, phig−phi0 (tbende.f)。ANGLE=0 では (0, 0)
    al0, phi0 = (el.L, float(el.ANGLE)) if phig != 0.0 else (0.0, 0.0)
    c3 = _tchge_cod(c2, dx, dy, dz, th, dth, alg - al0, phig - phi0, False)
    J2 = _tchge_jac(c2, dx, dy, dz, th, dth, alg - al0, phig - phi0, False)
    return c3, J2 @ M @ J1


def _num_jacobian(f, cod: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    """COD 位置での要素写像 f の 6×6 Jacobian (中心差分)。"""
    M = np.zeros((6, 6))
    for j in range(6):
        cp = cod.copy(); cp[j] += eps
        cm = cod.copy(); cm[j] -= eps
        M[:, j] = (f(cp) - f(cm)) / (2.0 * eps)
    return M


# 4D CALC の要素カーネル (elkern: 放射なし・dv なし・z 行の組み直しなし)。使えない要素は None を返し従来の Python 版で計算する
_USE_ELKERN = os.environ.get('PYSAD_NO_ELKERN') != '1'
_EK_PH4 = None
_EK_Z66 = np.zeros((6, 6))
_EK_Z6 = np.zeros(6)
_EK_S3 = np.zeros(3)


def _p2h_ek(p):
    return math.hypot(1.0, p)


def _ek_ph4():
    global _EK_PH4
    if _EK_PH4 is None:
        from .emit6d import _Phys
        from .elkern import phys_array
        ph = _Phys(1e9, 1.0)
        _EK_PH4 = phys_array(ph)
        _EK_S3[1] = _p2h_ek(ph.p0)       # h1emit (K1≠0 BEND の z 行 (5,6) にだけ入る。4D の CALC は z 行を使わないが有限にしておく)
    return _EK_PH4


def _bend_step_fast(el, cod):
    """BEND の (出口 COD, 転送行列) を elkern.k_bend で (tbende: tchge → 本体 → tchge)。F1 フリンジ・tmulte は None"""
    from . import elkern as _ek
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if g('F1') != 0.0 or g('FB1') != 0.0 or g('FB2') != 0.0:
        return None
    ang = float(el.ANGLE)
    k0 = g('K0')
    k1 = g('K1')
    if ang == 0.0 and k1 != 0.0 and k0 != 0.0:
        return None
    chg = _bend_chg_args(el)
    if chg is None:
        has_chg = False
        dx = dy = dz = th = dth = alg = phig = 0.0
    else:
        has_chg = True
        dx, dy, dz, th, dth, alg, phig = chg
    c = np.array(cod, dtype=float)
    T = np.eye(6)
    r = _ek.k_bend(c, T, _EK_Z66, _EK_Z66, _EK_Z6, _EK_Z66, _EK_S3, _ek_ph4(), float(el.L), ang, k0, k1,
                   g('E1'), g('E2'), has_chg, dx, dy, dz, th, dth, alg, phig,
                   False, False, False, True, False, False)
    if r != 0:
        return None
    return c, T


def _quad_step_fast(el, cod, K1):
    """QUAD (L≠0) の (出口 COD, 転送行列) を elkern._quad_full で (tquade の解析行列)"""
    from . import elkern as _ek
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    f1in, f2in, f1out, f2out, mfring = _quad_fringe_params(el)
    c = np.array(cod, dtype=float)
    ok, J = _ek._quad_full(c, float(el.L), K1, g('DISFRIN') == 0.0, f1in, f2in, f1out, f2out, int(mfring),
                           g('DX'), g('DY'), g('ROTATION'), False, _ek_ph4(), True)
    if not ok:
        return None
    return c, J


def _thin_step_fast(L, ak, kord, dx, dy, th, cod):
    """薄い多極 (tthine + tchge) の (出口 COD, 転送行列) を elkern.k_thin で"""
    from . import elkern as _ek
    c = np.array(cod, dtype=float)
    T = np.eye(6)
    r = _ek.k_thin(c, T, _EK_Z66, _EK_Z66, _EK_S3, _ek_ph4(), float(L), float(ak), int(kord), dx, dy, th,
                   False, False, True, False)
    if r != 0:
        return None
    return c, T


def _el_has_offset(el, delta) -> bool:
    """設計軌道 (delta=0, cod=0) でも線形転送行列で済ませられない要素か (ずれ・ロール・キック・フリンジ・薄い四極・MULT の双極/四極)"""
    return (bool(np.any(el.constant_kick(delta)))
            or float(getattr(el, 'DX', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'DY', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'ROTATION', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'ROTATE', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'DROTATE', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'F1', 0.0) or 0.0) != 0.0
            or float(getattr(el, 'FRINGE', 0.0) or 0.0) != 0.0
            or (type(el).__name__ == 'Quadrupole' and el.L == 0.0
                and float(getattr(el, 'K1', 0.0) or 0.0) != 0.0)      # 薄い四極 (SQ): transfer_matrix は恒等
            or (type(el).__name__ == 'Multipole'
                and any(getattr(el, 'KnL', {}).get(k, 0.0) for k in (0, 1))))


def _mult_is_drift(el) -> bool:
    """K_n が全部 0 (VOLT, ANGLE, BZ も 0) の MULT: tmulte は tdrife0 (厳密ドリフト) だけ"""
    if type(el).__name__ != 'Multipole':
        return False
    if any(v for v in (getattr(el, 'KnL', {}) or {}).values()) or any(v for v in (getattr(el, 'KsL', {}) or {}).values()):
        return False
    return not any(float(getattr(el, k, 0.0) or 0.0) for k in ('VOLT', 'ANGLE', 'BZ', 'K0'))


def _el_cod_step(el, cod, delta):
    """1 要素分の (COD 依存 Jacobian M_el, 出口 COD)。_one_turn_cod_matrix の要素ごとの本体 (キャッシュしない版)"""
    from .elements import Drift as _Drift, Quadrupole as Quad
    # DX/DY オフセットや補正コイルキック (K0/ROTATE) があれば近道しない (設計軌道上のときだけ評価する)
    if delta == 0.0 and not np.any(cod[:4]) and not _el_has_offset(el, delta):
        # 設計軌道上 (dp=0, cod=0): COD 依存 Jacobian = 線形転送行列。
        # 数値微分の丸め誤差を避けるため解析行列を使う。
        if type(el).__name__ == 'Cavity' and not RFSW:
            M_el = _drift_matrix_at_cod(el.L, delta, 0.0, 0.0)   # RFSW OFF: 空洞はドリフト (tcave rfsw=.false.)
        else:
            M_el = el.transfer_matrix(delta)
        return M_el, M_el @ cod

    if isinstance(el, _Drift) or (type(el).__name__ == 'Cavity' and not RFSW) or _mult_is_drift(el):
        # SAD の CALC は RFSW OFF が既定: 空洞 (tcave, rfsw=.false.) は cod を変えず、行列は長さ分のドリフト。
        # K が全部 0 の MULT も tdrife0 (tmulte.f:130-135) = 厳密ドリフト (数値微分を使わない)
        pxi = cod[1]
        pyi = cod[3]
        M_el = _drift_matrix_at_cod(el.L, delta, pxi, pyi)
        # COD位置更新 (qtwiss.f: cod(1)=cod(1)+pxi*ale)
        pr = 1.0 + delta
        a = pxi**2 + pyi**2
        if a > 0.0 and abs(pr) > 1e-15:
            arg = max(0.0, 1.0 - a / pr**2)
            dpz = -a / pr / (1.0 + np.sqrt(arg))
        else:
            dpz = 0.0
        pzi = pr + dpz
        if abs(pzi) < 1e-15:
            pzi = pr
        ale = el.L / pzi
        cod_out = cod.copy()
        cod_out[0] += pxi * ale
        cod_out[2] += pyi * ale
        cod_out[4] += dpz * ale          # qtwiss.f:224 cod(5)=cod(5)+dpz*ale (dvemit=0): 経路長の変化 -(px²+py²)/2·L
    elif isinstance(el, Bend):
        r = _bend_step_fast(el, cod) if _USE_ELKERN else None
        if r is not None:
            cod_out, M_el = r
        else:
            cod_out, M_el = _bend_element_cod_matrix(el, cod, delta)

    elif isinstance(el, Quad):
        # QUAD: SADのtquade相当の非線形COD追跡
        # bz=0, dx=dy=0, 放射なし (tquade.f直訳)
        # K1<0の場合: -90°回転 → tsolque(|K1|) → +90°回転
        # Jacobian は COD 位置での数値微分 (SADのqquad: tquade(mat=.true.) 相当)
        # DX/DY: tsolrote (入口で x-=dx, y-=dy → 回転, 出口で逆)
        K1 = getattr(el, 'K1', 0.0)
        if el.L == 0.0:
            # tquade.f: al=0 → tthine(nord=4, ak=K1 (積分値), dx, dy, theta) の薄い四極 (pf2 の SQ など。K1<0 の −90° 回転はない)
            r = _thin_step_fast(0.0, float(K1), 1, float(getattr(el, 'DX', 0.0) or 0.0),
                                float(getattr(el, 'DY', 0.0) or 0.0),
                                float(getattr(el, 'ROTATION', 0.0) or 0.0), cod) if _USE_ELKERN else None
            if r is not None:
                cod_out, M_el = r
            else:
                cod_out, M_el = _thin_cod_matrix(0.0, float(K1), 1, float(getattr(el, 'DX', 0.0) or 0.0),
                                                 float(getattr(el, 'DY', 0.0) or 0.0),
                                                 float(getattr(el, 'ROTATION', 0.0) or 0.0), cod)
        elif _USE_ELKERN and (r := _quad_step_fast(el, cod, float(K1))) is not None:
            cod_out, M_el = r
        else:
            # DISFRIN≠0 なら tqfrie (ハードエッジのフリンジ) を切る (SAD: fringe = (DISFRIN==0))
            _f = _quad_offset_map(el.L, K1, getattr(el, 'DX', 0.0),
                                  getattr(el, 'DY', 0.0),
                                  float(getattr(el, 'ROTATION', 0.0) or 0.0),
                                  float(getattr(el, 'DISFRIN', 0.0) or 0.0) == 0.0,
                                  _quad_fringe_params(el))
            cod_out = _f(cod)
            # tquade の解析行列 (tsolrote のずれ・ロール、tqfrie、tqlfre、tqente + tsolque 本体の積)。
            # SAD は CALC で数値微分をしない。解析式が作れない (pz が潰れる等) ときだけ数値微分
            from .emit6d import _jac_quad_full
            M_el = _jac_quad_full(cod, el.L, K1, float(getattr(el, 'DISFRIN', 0.0) or 0.0) == 0.0,
                                  _quad_fringe_params(el), None, float(getattr(el, 'DX', 0.0) or 0.0),
                                  float(getattr(el, 'DY', 0.0) or 0.0), float(getattr(el, 'ROTATION', 0.0) or 0.0))
            if M_el is None:
                M_el = _num_jacobian(_f, cod)

    elif type(el).__name__ in _THIN_KORD:
        # SEXT/OCT/DECA/DODECA: SAD qthin (qtwiss.f) → tthine.f の COD 写像。
        # pysad の K は長さあたり → SAD の積分値 K·L を渡す (L=0 なら K をそのまま積分値とみなす)。
        # DX/DY/ROT は SAD tchge (alg=phig=dz=dtheta=dchi2=0) と同じく座標をずらして通す。
        _kn, _kord = _THIN_KORD[type(el).__name__]
        _K = float(getattr(el, _kn, 0.0) or 0.0)
        _ak = _K * el.L if el.L != 0.0 else _K
        _dx, _dy, _th = (float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0),
                         float(getattr(el, 'ROTATION', 0.0) or 0.0))
        r = _thin_step_fast(el.L, _ak, _kord, _dx, _dy, _th, cod) if _USE_ELKERN else None
        if r is not None:
            cod_out, M_el = r
        else:
            cod_out, M_el = _thin_cod_matrix(el.L, _ak, _kord, _dx, _dy, _th, cod)

    elif type(el).__name__ == 'Multipole':
        # MULT (ANGLE=0, BZ=0): SAD qtwiss → tmulte の非線形 COD 写像 (tmulte.f 直訳の列)
        _f = _mult_cod_map(el)
        cod_out = _f(cod)
        M_el = _num_jacobian(_f, cod)

    else:
        M_el = el.transfer_matrix(delta)
        cod_out = M_el @ cod + el.constant_kick(delta)
    return M_el, cod_out


# 要素ごとの結果キャッシュ: 要素の属性・入口 COD・delta・RFSW が前回と同じなら前回の (M_el, cod_out) を返す (結果は同一)。
# key = id(el) → (属性の浅いコピー, cod.tobytes(), delta, RFSW, M_el, cod_out)。M_el/cod_out は読み出し専用
_EL_SCALARS = (float, int, str, bool)
_EL_CACHE: dict = {}
_EL_CACHE_ON = True


def _el_cod_step_cached(el, cod, delta):
    d = el.__dict__
    ent = _EL_CACHE.get(id(el))
    cb = cod.tobytes()
    if ent is not None and ent[1] == cb and ent[2] == delta and ent[3] == RFSW and ent[0] is el and ent[4] == d:
        return ent[5], ent[6].copy()
    M_el, cod_out = _el_cod_step(el, cod, delta)
    for k, v in d.items():
        if k != '_raw' and v.__class__ not in _EL_SCALARS:
            _EL_CACHE.pop(id(el), None)
            return M_el, cod_out
    if len(_EL_CACHE) > 60000:
        _EL_CACHE.clear()               # USE のたびに要素が作り直される場合の肥大化を防ぐ
    M_ro = np.array(M_el, dtype=float)
    M_ro.flags.writeable = False
    c_ro = np.array(cod_out, dtype=float)
    c_ro.flags.writeable = False
    _EL_CACHE[id(el)] = (el, cb, delta, RFSW, dict(d), M_ro, c_ro)
    return M_ro, c_ro.copy()


# ---------------------------------------------------------------------------
# 4D 1 周カーネル (pysad/turnkern.py): 要素のパラメータを配列に詰めて (_TurnPack) 1 周を numba 1 回で回す。
# pack は (要素数, delta, RFSW) ごとに保持し、要素ごとの __dict__ の写しと比べて変わった要素だけ詰め直す。
# 設計軌道上の線形行列 T0 は要るときにカーネルが Python に戻して作る (lazy)。
# ---------------------------------------------------------------------------
_USE_TURNKERN = os.environ.get('PYSAD_NO_TURNKERN') != '1'
_TURN_PACKS: dict = {}


class _TurnPack:
    __slots__ = ('els', 'snaps', 'kind', 'has_off', 'P', 'Mlin', 'Klin', 'T0', 'hasT0', 'delta')

    def __init__(self, n, delta):
        from .turnkern import NP
        self.els = [None] * n
        self.snaps = [None] * n
        self.kind = np.zeros(n, dtype=np.int64)
        self.has_off = np.zeros(n, dtype=np.bool_)
        self.P = np.zeros((n, NP))
        self.Mlin = np.zeros((n, 6, 6))
        self.Klin = np.zeros((n, 6))
        self.T0 = np.zeros((n, 6, 6))
        self.hasT0 = np.zeros(n, dtype=np.bool_)
        self.delta = delta


def _pack_el(pk: '_TurnPack', i: int, el, delta: float):
    """要素 el を pack の i 番目に詰める (_el_cod_step の分岐と同じ種別分け)"""
    from .elements import Drift as _Drift, Quadrupole as Quad
    from . import turnkern as _tk
    d = el.__dict__
    pk.els[i] = el
    snap = dict(d)
    for k, v in d.items():
        if k != '_raw' and v.__class__ not in _EL_SCALARS:
            snap = None
            break
    pk.snaps[i] = snap
    pk.hasT0[i] = False
    pk.has_off[i] = _el_has_offset(el, delta) if delta == 0.0 else False
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    tn = type(el).__name__
    P = pk.P[i]
    P[:] = 0.0
    if isinstance(el, _Drift) or (tn == 'Cavity' and not RFSW):
        pk.kind[i] = _tk.K_DRIFT
        P[0] = float(el.L)
        P[1] = 1.0 if tn == 'Cavity' else 0.0
    elif isinstance(el, Bend):
        ang, k0, k1 = float(el.ANGLE), g('K0'), g('K1')
        if g('F1') != 0.0 or g('FB1') != 0.0 or g('FB2') != 0.0 or (ang == 0.0 and k1 != 0.0 and k0 != 0.0):
            pk.kind[i] = _tk.K_PY
            return
        pk.kind[i] = _tk.K_BEND
        chg = _bend_chg_args(el)
        P[0:6] = (float(el.L), ang, k0, k1, g('E1'), g('E2'))
        if chg is not None:
            P[6] = 1.0
            P[7:14] = chg
    elif isinstance(el, Quad):
        K1 = float(getattr(el, 'K1', 0.0))
        if el.L == 0.0:
            pk.kind[i] = _tk.K_THIN
            P[0:6] = (0.0, K1, 1, g('DX'), g('DY'), g('ROTATION'))
        else:
            pk.kind[i] = _tk.K_QUAD
            f1in, f2in, f1out, f2out, mfring = _quad_fringe_params(el)
            P[0:11] = (float(el.L), K1, 1.0 if g('DISFRIN') == 0.0 else 0.0, f1in, f2in, f1out, f2out, int(mfring),
                       g('DX'), g('DY'), g('ROTATION'))
    elif tn in _THIN_KORD:
        _kn, _kord = _THIN_KORD[tn]
        _K = g(_kn)
        pk.kind[i] = _tk.K_THIN
        P[0:6] = (float(el.L), _K * el.L if el.L != 0.0 else _K, _kord, g('DX'), g('DY'), g('ROTATION'))
    elif tn == 'Multipole':
        pk.kind[i] = _tk.K_PY
    else:
        pk.kind[i] = _tk.K_LIN
        pk.Mlin[i] = el.transfer_matrix(delta)
        pk.Klin[i] = el.constant_kick(delta)


def _get_turn_pack(elements, delta: float) -> '_TurnPack':
    """elements の pack を (作って/検証して) 返す。変わった要素 (別のオブジェクト、属性が違う) だけ詰め直す"""
    n = len(elements)
    key = (n, delta, RFSW)
    pk = _TURN_PACKS.get(key)
    if pk is None:
        if len(_TURN_PACKS) >= 8:
            _TURN_PACKS.clear()
        pk = _TurnPack(n, delta)
        _TURN_PACKS[key] = pk
        for i, el in enumerate(elements):
            _pack_el(pk, i, el, delta)
        return pk
    els, snaps = pk.els, pk.snaps
    for i, el in enumerate(elements):
        sn = snaps[i]
        if els[i] is not el or sn is None:
            _pack_el(pk, i, el, delta)
            continue
        try:
            same = el.__dict__ == sn
        except Exception:
            same = False
        if not same:
            _pack_el(pk, i, el, delta)
    return pk


def _design_matrix(el) -> np.ndarray:
    """設計軌道上 (delta=0, cod=0, ずれなし) の要素行列 (_el_cod_step の近道と同じ)"""
    if type(el).__name__ == 'Cavity' and not RFSW:
        return _drift_matrix_at_cod(el.L, 0.0, 0.0, 0.0)
    return el.transfer_matrix(0.0)


def _one_turn_cod_matrix(elements, delta: float,
                         cod_init: np.ndarray, pack: '_TurnPack' = None):
    """
    SADのqtwiss1(mat=.true.)相当。

    初期COD cod_init (6,) でリングを1周し、
    COD依存Jacobian (6×6) と出口COD (6,) を返す。
    戻り値の el_mats は [(要素行列 M_el, 入口 COD), ...]。
    pack: _get_turn_pack の結果 (呼び出し側が反復の間使い回す。None なら作る)
    """
    cod = np.array(cod_init, dtype=float)
    M_total = np.eye(6)
    step = _el_cod_step_cached if _EL_CACHE_ON else _el_cod_step
    if not (_USE_ELKERN and _USE_TURNKERN):
        el_mats = []
        for el in elements:
            M_el, cod_out = step(el, cod, delta)
            el_mats.append((M_el, cod.copy()))
            M_total = M_el @ M_total
            cod = cod_out
        return M_total, cod, el_mats

    from .turnkern import turn4
    if pack is None:
        pack = _get_turn_pack(elements, delta)
    n = len(elements)
    mats = np.empty((n, 6, 6))
    cods = np.empty((n, 6))
    i = 0
    while i < n:
        i, why = turn4(i, pack.kind, pack.has_off, pack.P, pack.Mlin, pack.Klin, pack.T0, pack.hasT0,
                       cod, M_total, mats, cods, delta, _ek_ph4(), _EK_S3)
        if i >= n:
            break
        el = elements[i]
        if why == 1:
            pack.T0[i] = _design_matrix(el)
            pack.hasT0[i] = True
            continue
        M_el, cod_out = step(el, cod, delta)
        cods[i] = cod
        mats[i] = M_el
        M_total[:, :] = M_el @ M_total
        cod[:] = cod_out
        i += 1
    return M_total, cod, list(zip(mats, cods))


def _tt_scalar_only(el) -> bool:
    for k, v in el.__dict__.items():
        if k != '_raw' and v.__class__ not in _EL_SCALARS:
            return False
    return True


def _tt_copy(v):
    """memo に出し入れする値のコピー (配列・辞書・リストは中身まで 1 段コピー。読み出し専用配列は共有)"""
    if isinstance(v, np.ndarray):
        return v if not v.flags.writeable else v.copy()
    if isinstance(v, dict):
        return dict(v)
    if isinstance(v, list):
        return list(v)
    return v


_ROW_NEW = object.__new__


def _clone_rows(rows):
    """TwissRow のリストの浅いコピー (行ごとに独立。呼び出し側が属性を書き換えても memo に影響しない)"""
    out = []
    ap = out.append
    for r in rows:
        c = _ROW_NEW(TwissRow)
        c.__dict__.update(r.__dict__)
        ap(c)
    return out



def _qcod_periodic(elements, delta: float,
                   cod_init: np.ndarray = None, fixed: bool = False) -> tuple:
    """
    SADのqcod (qtwiss.f) 直訳: 周期閉軌道を反復で求める。

    - 収束判定: r = (Courant-Snyder 不変量で測った |cod_exit - cod0|^2) <= conv
        conv = min(1e-10, 1e-18*(1+(dp/0.001)^2+nlat))
    - r が増えたらダンピング (fact) して後退、減れば Newton 更新
    - 返す M_cod は「最後に評価した入口 cod0」での一周 Jacobian
      (SAD も収束点ではなく最後の評価点の trans で Twiss を計算する)
    - cod_init: 初期値 (SAD は直前の CALC の閉軌道を引き継ぐ)。None なら 0。
    - fixed: qcod を codfnd=.true. で呼んだとき (fam 点): 閉軌道を探さず、cod_init のまわりの一周行列 (1 周だけ) を返す

    注意: SAD の収束判定が緩いため、off-momentum の Twiss は初期値に依存して
    相対 1e-6 程度揺らぐ (SAD 自身の数値誤差)。

    Returns: (cod_periodic(入口), M_cod, el_mats)
    """
    nlat = len(elements) + 1
    itmax = 10
    conv0 = 1e-18
    conv1 = 1e-10
    factmin = 1e-3
    cod0 = np.zeros(6) if cod_init is None else np.array(cod_init, dtype=float)
    cod0[4] = 0.0
    cod0[5] = delta
    it = 0
    r0 = 1e100
    cod00 = cod0.copy()
    fact = 0.5
    conv = min(conv1, conv0 * (1.0 + (cod0[5] / 0.001)**2 + nlat))
    M_cod = None
    el_mats = []
    cod_in = cod0.copy()
    found = True
    over = False
    _mmax = 0.0
    pack = _get_turn_pack(elements, delta) if (_USE_ELKERN and _USE_TURNKERN) else None
    while it <= itmax:
        cod_in = cod0.copy()
        with np.errstate(all='ignore'):
            M_cod, cod, el_mats = _one_turn_cod_matrix(elements, delta, cod0, pack)
        # qtwiss1 (mat=.true.): 各要素の入口で一周行列の跡 (対角和) が NaN なら over=.true. で打ち切り、qcod は codfnd=.false. で戻る
        # (Newton の途中で試した軌道が大振幅になり NaN になると、SAD は後退せずに 'Overflow & closed orbit not found' で終える)
        if np.isnan(M_cod[0, 0] + M_cod[1, 1] + M_cod[2, 2] + M_cod[3, 3]):
            over = True
            M_cod = M_cod.copy()
            for _blk in (M_cod[:4, :4], M_cod[:4, 5:6]):
                _blk[np.isnan(_blk)] = _TRMAX
                np.clip(_blk, -_TRMAX, _TRMAX, out=_blk)
            found = False
            break
        # qcod: call limitnan(tr1v,trmax): 一周の行列 (4x4 と分散の列) の成分を ±1e10 に抑え、NaN は +1e10 にする。
        # 不安定でこの値を超えるとき (M の成分が 1e10〜) の Twiss は、抑えた行列から作られる
        _mmax = float(np.max(np.abs(M_cod[:4, :4]), initial=0.0))
        if _mmax > _TRMAX or np.any(np.isnan(M_cod[:4, :4])):
            M_cod = M_cod.copy()
            for _blk in (M_cod[:4, :4], M_cod[:4, 5:6]):
                _blk[np.isnan(_blk)] = _TRMAX
                np.clip(_blk, -_TRMAX, _TRMAX, out=_blk)
        if fixed:
            found = True                  # qcod: codfnd が真 (fam 点) なら cod0=cod で即 return (与えた軌道のまわりの一周行列)
            break
        trans = M_cod
        cx = 0.5 * (trans[0, 0] + trans[1, 1])
        cy = 0.5 * (trans[2, 2] + trans[3, 3])
        if abs(cx) > 1.0:
            cx = 1.0 / cx
        if abs(cy) > 1.0:
            cy = 1.0 / cy
        sx = np.copysign(np.sqrt(max(1e-6, 1.0 - cx**2)), trans[0, 1])
        ax = 0.5 * (trans[0, 0] - trans[1, 1]) / sx
        bx = trans[0, 1] / sx
        sy = np.copysign(np.sqrt(max(1e-6, 1.0 - cy**2)), trans[2, 3])
        ay = 0.5 * (trans[2, 2] - trans[3, 3]) / sy
        by = trans[2, 3] / sy
        dcod = cod - cod0
        r = (dcod[0]**2 / bx + bx * (dcod[1] + ax / bx * dcod[0])**2
             + dcod[2]**2 / by + by * (dcod[3] + ay / by * dcod[2])**2)
        if not np.isfinite(r):
            r = 1e300
        if r <= conv:
            break
        it += 1
        if r > r0:
            if fact < factmin:
                fact = fact * 16.0
                cod0 = (1.0 + fact) * cod00 - fact * cod0
            else:
                cod0 = (1.0 - fact) * cod00 + fact * cod0
                fact = abs(fact) * 0.5
        else:
            fact = min(0.5, fact * 2.0)
            r0 = r
            cod00 = cod0.copy()
            rhs = -cod[:4] + trans[:4, :4] @ cod0[:4]
            cod0 = cod0.copy()
            try:
                cod0[:4] = _solve_cod(trans[:4, :4] - np.eye(4), rhs)
            except np.linalg.LinAlgError:
                found = False
                break
        cod0[5] = delta
    else:
        found = False                 # itmax を使い切った (SAD: cod=cod0; codfnd=.false.)

    global _LAST_COD_EXIT, _LAST_QCOD_FOUND, _LAST_MMAX, _LAST_QCOD_OVER
    _LAST_QCOD_OVER = over
    _LAST_MMAX = _mmax if np.all(np.isfinite(M_cod)) else float('inf')
    _LAST_COD_EXIT = np.array(cod, dtype=float)     # 最後に評価した一周の出口 COD (SAD の twiss(nlat) の COD。入口と conv 分だけずれる)
    _LAST_QCOD_FOUND = found
    return cod_in, M_cod, el_mats


_LAST_COD_EXIT = None
_TRMAX = 1e10                  # qcod の trmax
_LAST_MMAX = 0.0               # 直近の _qcod_periodic の一周行列の最大成分 (クリップ前)
_MOVERFLOW = 1e18              # これを超える不安定な一周行列は SAD では COD の丸め雑音 (1e-18) が非線形領域まで増幅されて NaN/Overflow になる
_LAST_QCOD_FOUND = True        # 直近の _qcod_periodic が収束 (r <= conv) したか (SAD の codfnd)
_LAST_QCOD_OVER = False        # 直近の _qcod_periodic が qtwiss1 の over (一周行列の跡が NaN) で打ち切られたか


def _solve_cod(a: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """qcod の tsolvg(trans1, cod, cod0) (= tsvd, 閾値 1e-8)。条件の良い行列は厳密解 (numpy の solve) と同じなのでそれを使い、
    特異に近い (整数・半整数共鳴の近く: 不安定ラインの閉軌道探索) ときは SAD と同じ tsvd (暗黙 QR の滑らかな減衰) で解く"""
    if np.all(np.isfinite(a)):
        sv = np.linalg.svd(a, compute_uv=False)
        if sv[0] > 0.0 and sv[-1] > 1e-8 * sv[0]:
            return np.linalg.solve(a, rhs)
    from .tsvd import tsvd
    return np.nan_to_num(tsvd(a, rhs, 1e-8))


def _bend_slice_matrix(h, k, dl):
    """
    INTGRL.f BMAG に相当 (TETA=0, 非結合) 。
    h   = ANG/L  (全体曲率 1/ρ)
    k   = K1     (単位長あたり四極子強度 [m^-2]; el.K1 は既に per-length)
    dl  = スライス長

    Returns: (m11,m12,m21,m22, vx1,vx2,  m33,m34,m43,m44)
      [x-plane 2×2, x-dispersion orbit 2, y-plane 2×2]
    """
    import math as _m
    kx = h*h + k          # x 方向フォーカス強度
    ky = -k               # y 方向

    # --- x-plane ---
    if kx > 1e-30:
        skx = _m.sqrt(kx)
        phi = skx * dl
        c = _m.cos(phi); s = _m.sin(phi)
        m11, m12, m21, m22 = c, s/skx, -skx*s, c
        vx1 = h/kx * (1.0 - c)
        vx2 = h/skx * s
    elif kx < -1e-30:
        skx = _m.sqrt(-kx)
        phi = skx * dl
        ch = _m.cosh(phi); sh = _m.sinh(phi)
        m11, m12, m21, m22 = ch, sh/skx, skx*sh, ch
        vx1 = h/kx * (1.0 - ch)
        vx2 = h/skx * sh
    else:
        m11, m12, m21, m22 = 1.0, dl, 0.0, 1.0
        vx1 = 0.5 * h * dl*dl
        vx2 = h * dl

    # --- y-plane ---
    if ky > 1e-30:
        sky = _m.sqrt(ky)
        psi = sky * dl
        c = _m.cos(psi); s = _m.sin(psi)
        m33, m34, m43, m44 = c, s/sky, -sky*s, c
    elif ky < -1e-30:
        sky = _m.sqrt(-ky)
        psi = sky * dl
        ch = _m.cosh(psi); sh = _m.sinh(psi)
        m33, m34, m43, m44 = ch, sh/sky, sky*sh, ch
    else:
        m33, m34, m43, m44 = 1.0, dl, 0.0, 1.0

    return (m11, m12, m21, m22, vx1, vx2, m33, m34, m43, m44)


def _symps(z, n):
    """
    INTGRL.f SYMPS(Z, N) の Python 版。
    z : 長さ 2n+1 のリスト (0-indexed)。Fortran の Z(1)..Z(2N+1) に対応。
    n : IR (パネル数)
    composite Simpson 1/3 則:
      (z[0] + z[2n] + 2*Σ_even_interior + 4*Σ_odd) / 3
    """
    if n == 0:
        return z[0]
    s2 = sum(z[2*i]     for i in range(1, n))   # z[2],z[4],...,z[2n-2]
    s3 = sum(z[2*i - 1] for i in range(1, n+1)) # z[1],z[3],...,z[2n-1]
    return (z[0] + z[2*n] + 2.0*s2 + 4.0*s3) / 3.0


def _radiation_integrals_bend(el: Bend, tw_in: dict) -> dict:
    """
    偏向電磁石一本の放射積分寄与 I1..I5 を計算。
    SAD INTGRL.f の BEND 処理を忠実に再現:
      * IR = min(14,max(1,int(|ANG|/0.002))) スライス, IRP=2*IR サブ区間
      * composite Simpson 1/3 則 (SYMPS)
      * スライス毎 Twiss 伝播 (TWTRANS 相当、非結合 TETA=0 前提)
      * BEDGE(E1/E2) + TEDGE → I4 端部補正

    Parameters
    ----------
    el     : Bend 要素
    tw_in  : 要素入口での Twiss dict  {'bx','ax','ex','epx', ...}

    Notes
    -----
    非結合 (TETA=0) を前提とするため R0I = 単位行列。
    SAD の B' 積分 K1: el.K1 は単位長あたり強度 [m^-2] (pysad パーサが L で割る)。
    tbend.f tbendiinit: akk=ak1/al → ak1 (ky_K1_BEND) は integrated [m^-1]、
    INTGRL.f: DKL1=RKL1/IRP → DKL1 = K1_per_length * DL (per-length × slice length)。
    """
    import math as _m
    L     = el.L
    angle = el.ANGLE
    K1    = el.K1   # el.K1 は単位長あたり [m^-2] (pysad パーサが integrated/L を格納)
    E1    = getattr(el, 'E1', 0.0)
    E2    = getattr(el, 'E2', 0.0)

    if abs(L) < 1e-12 or abs(angle) < 1e-12:
        return {f'I{i}': 0.0 for i in range(1, 6)}

    h  = angle / L
    k  = K1          # K1 は既に per-unit-length [m^-2]; INTGRL.f akk=ak1/al に対応

    # SAD INTGRL.f と同一のスライス数決定
    IR  = min(14, max(1, int(abs(angle) / 0.002)))
    IRP = 2 * IR
    DL  = L   / IRP
    DKL0 = angle / IRP
    DKL1 = K1   * DL  # INTGRL.f: DKL1=RKL1/IRP = K1_integrated/IRP = K1_per_length*DL

    CURV  = DKL0 / DL       # = h
    CURV2 = CURV * CURV
    CURV3 = CURV2 * abs(CURV)

    # --- Twiss 初期化 (tw_in から) ---
    bx   = tw_in['bx']
    ax   = tw_in['ax']
    gx   = (1.0 + ax*ax) / bx
    eta  = tw_in['ex']
    etap = tw_in['epx']

    # --- BEDGE(E1) + TEDGE: 入口端部補正 ---
    # DIEG(4) -= h^2 * tan(E1) * eta_in   (TEDGE 相当)
    # TWTRANS: eta 不変, etap += h*tan(E1)*eta, ax -= h*tan(E1)*bx
    dieg4 = 0.0
    if abs(E1) > 1e-15:
        tanE1 = _m.tan(E1)
        dieg4 -= h * h * tanE1 * eta
        etap  += h * tanE1 * eta
        ax    -= h * tanE1 * bx
        gx     = (1.0 + ax*ax) / bx

    # --- スライス行列 (全スライス同一: 均一磁場) ---
    m11, m12, m21, m22, vx1, vx2, _, _, _, _ =         _bend_slice_matrix(h, k, DL)

    # --- FUNI 配列 (IRP+1 点) ---
    funi1 = [0.0] * (IRP + 1)
    funi2 = [0.0] * (IRP + 1)
    funi3 = [0.0] * (IRP + 1)
    funi4 = [0.0] * (IRP + 1)
    funi5 = [0.0] * (IRP + 1)

    for I in range(1, IRP + 2):      # I=1..IRP+1 (Fortran の I=1..IRP+1 に対応)
        # TCAL: 現在 Twiss で FUNI 計算
        Hx = gx*eta*eta + 2.0*ax*eta*etap + bx*etap*etap
        idx = I - 1
        funi1[idx] = DKL0 * eta
        funi2[idx] = CURV2 * DL
        funi3[idx] = CURV3 * DL
        if CURV > 0.0:
            funi4[idx] = (DKL0 + 2.0*DKL1/CURV) * CURV2 * eta
        else:
            funi4[idx] = 0.0
        funi5[idx] = CURV3 * Hx * DL

        if I == IRP + 1:
            break   # IFLG=1: 最後の点は TWTRANS しない

        # TWTRANS: Twiss を次スライス位置へ伝播
        bx_new  = m11*m11*bx - 2.0*m11*m12*ax + m12*m12*gx
        ax_new  = -m11*m21*bx + (1.0 + 2.0*m12*m21)*ax - m12*m22*gx
        eta_new  = m11*eta + m12*etap + vx1
        etap_new = m21*eta + m22*etap + vx2
        bx, ax = bx_new, ax_new
        gx = (1.0 + ax*ax) / bx
        eta, etap = eta_new, etap_new

    # --- BEDGE(E2) + TEDGE: 出口端部補正 ---
    # この時点の eta = ボディ出口 (エッジ前) の分散
    # 注: BEDGE では eta 変化なし (VM=0)
    if abs(E2) > 1e-15:
        dieg4 -= h * h * _m.tan(E2) * eta

    # --- Simpson 積分 ---
    DI1 = _symps(funi1, IR)
    DI2 = _symps(funi2, IR)
    DI3 = _symps(funi3, IR)
    DI4 = _symps(funi4, IR) + dieg4
    DI5 = _symps(funi5, IR)

    return {'I1': DI1, 'I2': DI2, 'I3': DI3, 'I4': DI4, 'I5': DI5}



# ---------------------------------------------------------------------------
# メイン Lattice クラス
# ---------------------------------------------------------------------------

def _qmdiag(tm: np.ndarray):
    """SAD qmdiag.f の直訳: 4x4 一周行列 tm を対角化する変換 (r1..r4, c)。
    戻り: (r1, r2, r3, r4, c, stab, nanq, msg)。stab=False は和共鳴 ((TrA-TrB)^2 < 0。msg に SAD のメッセージ)、
    nanq=True は行列が NaN/Inf。戻りの r, c は stab=False のとき意味を持たない"""
    (tm11, tm12, tm13, tm14), (tm21, tm22, tm23, tm24), \
        (tm31, tm32, tm33, tm34), (tm41, tm42, tm43, tm44) = [[float(x) for x in row] for row in tm[:4, :4]]
    r1 = r2 = r3 = r4 = c = 0.0
    trpq = tm11 + tm22 - tm33 - tm44
    trab2 = 4.0 * ((tm31 + tm24) * (tm42 + tm13) + (tm14 - tm32) * (tm41 - tm23)) + trpq**2
    if not math.isfinite(trab2):
        return r1, r2, r3, r4, c, False, True, ' ***qmdiag---> Matrix not a number: '
    if trab2 >= 1e-14:
        c2 = 0.5 * (1.0 + abs(trpq) / math.sqrt(trab2))
        c = math.sqrt(c2)
        aa = -1.0 / c / math.sqrt(trab2) if trpq >= 0.0 else 1.0 / c / math.sqrt(trab2)
        r1 = aa * (tm31 + tm24)
        r2 = aa * (tm32 - tm14)
        r3 = aa * (tm41 - tm23)
        r4 = aa * (tm42 + tm13)
        return r1, r2, r3, r4, c, True, False, None
    if abs(trab2) <= 1e-14:
        # (TrA-TrB)^2 = 0: 必然的に TrP = TrQ
        detn = tm13 * tm24 - tm14 * tm23
        if abs(detn) >= 1e-14:
            trqk = tm33 - tm44
            if abs(trqk) >= 1e-7:
                trqj = tm43 - tm34
                trqkj2 = (trqk + trqj) * (trqk - trqj)
                if abs(trqkj2) >= 1e-14:
                    c = 1.0
                    bbb = -detn / (trqj + trqk)
                    ccc = 0.0
                else:
                    c = 1.0
                    bbb = -0.5 * detn / trqj if abs(trqj) >= 1e-7 else 0.0
                    ccc = 0.0
                aaa = (detn * (1.0 - 2.0 * c * c) / c - bbb * trqj) / trqk
            else:
                trqj = tm43 - tm34
                if abs(trqj) >= 1e-7:
                    c = 1.0
                    aaa = detn / trqj
                    bbb = -detn / trqj
                    ccc = 0.0
                else:
                    trql = tm34 + tm43
                    if abs(trql) >= 1e-7:
                        c = 1.0
                        aaa = 0.0
                        bbb = detn / trql
                        ccc = -detn / trql
                    else:
                        c = 1.0 / math.sqrt(2.0)
                        if detn > 0.0:
                            aaa, bbb, ccc = 0.0, math.sqrt(detn / 2.0), 0.0
                        else:
                            aaa, bbb, ccc = math.sqrt(-detn / 2.0), 0.0, 0.0
            xx11 = aaa
            xx12 = bbb + ccc
            xx21 = -bbb + ccc
            xx22 = -aaa
            r1 = (xx11 * tm24 - xx12 * tm23) / detn
            r2 = (-xx11 * tm14 + xx12 * tm13) / detn
            r3 = (xx21 * tm24 - xx22 * tm23) / detn
            r4 = (-xx21 * tm14 + xx22 * tm13) / detn
        else:
            trqk = tm33 - tm44
            trqj = tm43 - tm34
            trql = tm43 + tm34
            done = True
            if abs(trqk) >= 1e-7:
                aaa, bbb, ccc = 0.0, 0.0, 1.0 / trqk
            elif abs(trqj) >= 1e-7:
                aaa, bbb, ccc = -1.0 / trqj, 0.0, 0.0
            elif abs(trql) >= 1e-7:
                aaa, bbb, ccc = 0.0, -1.0 / trql, 0.0
            else:
                done = False
            if done:
                c = 1.0
                yy11 = aaa + bbb
                yy12 = ccc
                yy21 = ccc
                yy22 = aaa - bbb
                r1 = yy11 * (-tm23) + yy12 * (-tm24)
                r2 = yy11 * tm13 + yy12 * tm14
                r3 = yy21 * (-tm23) + yy22 * (-tm24)
                r4 = yy21 * tm13 + yy22 * tm14
            else:
                trp = tm11 + tm22
                p1 = tm11 - 0.5 * trp
                p3 = tm21
                if abs(p3) >= 1e-7:
                    c, r1, r2, r3, r4 = 1.0, 0.0, tm24 / p3, 0.0, -tm23 / p3
                elif abs(p1) >= 1e-7:
                    c, r1, r2, r3, r4 = 1.0, tm24 / p1, 0.0, -tm23 / p1, 0.0
                else:
                    c = 1.0 / math.sqrt(2.0)
                    r1, r2, r3, r4 = 1.0, 0.0, 0.0, c
        return r1, r2, r3, r4, c, True, False, None
    return r1, r2, r3, r4, c, False, False, f' ***qmdiag---> Sum resonance: (TrA-TrB)^2 ={trab2:15.7e}'


def _qcell_periodic_sad(M_cod: np.ndarray, intres: bool, halfres: bool, init: Optional[dict] = None):
    """SAD qcell1 の「周期 Twiss」の部分 (qmdiag の後〜分散) の直訳。|cos μ| > 1 のとき cos μ → 1/cos μ に置き換えて
    チューンらしい値を作り、INTRES (cos>1) / HALFRES (cos<-1) が真ならその面を不安定とする。
    戻り (twiss0, info)。twiss0: ax bx ay by ex epx ey epy r1..r4 detr mux muy (置換後の μ。不安定でも値は入る)。
    info: stab (qmdiag の stab)、stabx/staby (ここまでの判定)、cosx/cosy (置換前)、
          msg (qmdiag のメッセージ)。和共鳴 (stab=False) のとき twiss0 は init (SAD の twiss(lb) の以前の値)
    init: 以前の Twiss (SAD の twiss(lb,idp,·) = 先頭 MARK の値)。不安定な面の ex/epx, ey/epy と和共鳴時の全値に使う"""
    base = dict(ax=0.0, bx=1.0, ay=0.0, by=1.0, ex=0.0, epx=0.0, ey=0.0, epy=0.0,
                r1=0.0, r2=0.0, r3=0.0, r4=0.0, detr=0.0)
    if init:
        base.update({k: float(v) for k, v in init.items() if k in base})
    r1, r2, r3, r4, c1, stab, nanq, msg = _qmdiag(M_cod)
    info = dict(stab=stab, nanq=nanq, msg=msg, stabx=True, staby=True, cosx=0.0, cosy=0.0)
    if nanq:
        info['stabx'] = info['staby'] = False
        return base, info
    if not stab:
        info['stabx'] = info['staby'] = False
        return base, info
    tm = np.asarray(M_cod, dtype=float)
    (tm11, tm12, tm13, tm14), (tm21, tm22, tm23, tm24), \
        (tm31, tm32, tm33, tm34), (tm41, tm42, tm43, tm44) = tm[:4, :4]
    tm15, tm25, tm35, tm45 = tm[0, 5], tm[1, 5], tm[2, 5], tm[3, 5]
    s11 = c1; s13 = -r4; s14 = r2
    s22 = c1; s23 = r3; s24 = -r1
    s31 = r1; s32 = r2; s33 = c1
    s41 = r3; s42 = r4; s44 = c1
    a11 = (s11 * (tm11 * s11 - tm13 * s31 - tm14 * s41)
           + s13 * (tm31 * s11 - tm33 * s31 - tm34 * s41)
           + s14 * (tm41 * s11 - tm43 * s31 - tm44 * s41))
    a12 = (s11 * (tm12 * s22 - tm13 * s32 - tm14 * s42)
           + s13 * (tm32 * s22 - tm33 * s32 - tm34 * s42)
           + s14 * (tm42 * s22 - tm43 * s32 - tm44 * s42))
    a21 = (s22 * (tm21 * s11 - tm23 * s31 - tm24 * s41)
           + s23 * (tm31 * s11 - tm33 * s31 - tm34 * s41)
           + s24 * (tm41 * s11 - tm43 * s31 - tm44 * s41))
    a22 = (s22 * (tm22 * s22 - tm23 * s32 - tm24 * s42)
           + s23 * (tm32 * s22 - tm33 * s32 - tm34 * s42)
           + s24 * (tm42 * s22 - tm43 * s32 - tm44 * s42))
    b11 = (s31 * (-tm11 * s13 - tm12 * s23 + tm13 * s33)
           + s32 * (-tm21 * s13 - tm22 * s23 + tm23 * s33)
           + s33 * (-tm31 * s13 - tm32 * s23 + tm33 * s33))
    b12 = (s31 * (-tm11 * s14 - tm12 * s24 + tm14 * s44)
           + s32 * (-tm21 * s14 - tm22 * s24 + tm24 * s44)
           + s33 * (-tm31 * s14 - tm32 * s24 + tm34 * s44))
    b21 = (s41 * (-tm11 * s13 - tm12 * s23 + tm13 * s33)
           + s42 * (-tm21 * s13 - tm22 * s23 + tm23 * s33)
           + s44 * (-tm41 * s13 - tm42 * s23 + tm43 * s33))
    b22 = (s41 * (-tm11 * s14 - tm12 * s24 + tm14 * s44)
           + s42 * (-tm21 * s14 - tm22 * s24 + tm24 * s44)
           + s44 * (-tm41 * s14 - tm42 * s24 + tm44 * s44))
    out = dict(base)
    out.update(r1=r1, r2=r2, r3=r3, r4=r4, detr=r1 * r4 - r2 * r3)
    a13 = s11 * tm15 + s13 * tm35 + s14 * tm45
    a23 = s22 * tm25 + s23 * tm35 + s24 * tm45
    b13 = s31 * tm15 + s32 * tm25 + s33 * tm35
    b23 = s41 * tm15 + s42 * tm25 + s44 * tm45
    with np.errstate(all='ignore'):
        for pl, (p11, p12, p21, p22, p13, p23, kx, kpx, ka, kb, kmu) in enumerate((
                (a11, a12, a21, a22, a13, a23, 'ex', 'epx', 'ax', 'bx', 'mux'),
                (b11, b12, b21, b22, b13, b23, 'ey', 'epy', 'ay', 'by', 'muy'))):
            cosmu = np.float64(0.5 * (p11 + p22))
            info['cosx' if pl == 0 else 'cosy'] = float(cosmu)
            stabp = True
            if cosmu > 1.0:
                cosmu = 1.0 / cosmu
                if intres:
                    stabp = False
            if cosmu < -1.0:
                cosmu = 1.0 / cosmu
                if halfres:
                    stabp = False
            sinmu = np.copysign(np.sqrt(1.0 - cosmu**2), p12)
            out[ka] = float((p11 - cosmu) / sinmu)
            out[kb] = float(p12 / sinmu)
            amu = np.arctan2(sinmu, cosmu)
            dcosmu = 2.0 * np.sin(0.5 * amu)**2
            out[kmu] = float(amu)
            if stabp and dcosmu != 0.0:
                out[kx] = float(0.5 * (p13 + p12 * p23 - p22 * p13) / dcosmu)
                out[kpx] = float(0.5 * (p23 + p21 * p13 - p11 * p23) / dcosmu)
            info['stabx' if pl == 0 else 'staby'] = stabp
    return out, info


def ins_twiss(elements, twiss0: dict, cod0=None, delta: float = 0.0, acc: dict = None) -> List[dict]:
    """
    輸送ライン (SAD の INS モード) の Twiss と軌道。
    入口の Twiss (twiss0: ax, bx, ay, by, ex, epx, (ey, epy)) と軌道 cod0 (x,px,y,py) から、
    各要素の非線形 COD 写像 (CELL と同じ _one_turn_cod_matrix の要素処理) で軌道を追い、
    その軌道上の行列で Twiss と分散を伝搬する (結合なしの伝搬式)。
    戻り値: 長さ nelem+1 のリスト。[0] = 入口、[i] = i 番目の要素の出口 (= i+1 番目の入口)。
    各要素は Twiss の辞書に 'cod' (6 成分) を加えたもの。
    acc: trptacc.make_acc の結果 (TRPT で加速要素があるとき)。加速要素は acc['step'] の写像で進め、行列は qcopymatg の
         sqrt(gammab(l+1)/gammab(l)) 倍、入口の分散は gr = gammab(l)/gammab(l+1) 倍 (qmat2twiss)。軌道は要素順に逐次追う
    """
    c0 = np.zeros(6)
    if cod0 is not None:
        c0[:len(cod0)] = cod0
    c0[5] = delta
    if acc is not None:
        return _ins_twiss_acc(elements, twiss0, c0, delta, acc)
    _, cod_end, el_mats = _one_turn_cod_matrix(elements, delta, c0)
    tw = dict(twiss0)
    tw.setdefault('nx', 0.0)
    tw.setdefault('ny', 0.0)
    for k in ('ey', 'epy', 'r1', 'r2', 'r3', 'r4'):
        tw.setdefault(k, 0.0)
    # qtwiss1 の入口: ex.. 以降の |値| < 1e-30 は 0 に丸め、normal = (DETR < 1) から DETR を作り直す
    for k in ('ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4'):
        if abs(tw[k]) < 1e-30:
            tw[k] = 0.0
    if tw.get('detr', 0.0) < 1.0:
        tw['detr'] = tw['r1'] * tw['r4'] - tw['r2'] * tw['r3']
    else:
        tw['detr'] = 1.0 + _XYTH - tw['r1'] * tw['r4'] + tw['r2'] * tw['r3']
    out = [dict(tw, cod=c0.copy())]
    cods = [m[1] for m in el_mats[1:]] + [cod_end]
    for (M_el, _), c in zip(el_mats, cods):
        tw = _transport_twiss(tw, M_el)
        out.append(dict(tw, cod=np.array(c, dtype=float)))
    return out


def _ins_twiss_init(twiss0: dict) -> dict:
    tw = dict(twiss0)
    tw.setdefault('nx', 0.0)
    tw.setdefault('ny', 0.0)
    for k in ('ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4'):
        tw.setdefault(k, 0.0)
    for k in ('ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4'):
        if abs(tw[k]) < 1e-30:
            tw[k] = 0.0
    if tw.get('detr', 0.0) < 1.0:
        tw['detr'] = tw['r1'] * tw['r4'] - tw['r2'] * tw['r3']
    else:
        tw['detr'] = 1.0 + _XYTH - tw['r1'] * tw['r4'] + tw['r2'] * tw['r3']
    return tw


def _ins_twiss_acc(elements, twiss0: dict, c0: np.ndarray, delta: float, acc: dict) -> List[dict]:
    """ins_twiss の TRPT 加速版: 要素を順に進める (加速要素は acc['step']、それ以外は _el_cod_step)。
    加速要素の行列は qcopymatg (sqrt(gb1/gb0) 倍)、分散は gr = gb0/gb1 (qmat2twiss の trpt)"""
    gb = acc['gb']
    idx = acc['idx']
    step_acc = acc['step']
    step = _el_cod_step_cached if _EL_CACHE_ON else _el_cod_step
    tw = _ins_twiss_init(twiss0)
    cod = np.array(c0, dtype=float)
    out = [dict(tw, cod=cod.copy())]
    for i, el in enumerate(elements):
        if i in idx:
            M, cod_out = step_acc(el, i, cod)
            g0, g1 = gb[i], gb[i + 1]
            if g1 != g0:
                M = M.copy()
                M[:4, :] *= math.sqrt(g1 / g0)
                gr = g0 / g1
            else:
                gr = 1.0
        else:
            # qtwiss の各要素は軌道の cod(6) を運動量のずれに使う (加速の後は cod(6) が DP0 からずれる)
            M, cod_out = step(el, cod, float(cod[5]))
            gr = 1.0
        tw = _transport_twiss(tw, M, gr=gr)
        cod = np.array(cod_out, dtype=float)
        out.append(dict(tw, cod=cod.copy()))
    return out


class Lattice:
    """
    SAD の LINE[] に相当する加速器格子クラス。

    Parameters
    ----------
    elements : list of Element
        格子を構成する要素のリスト (順番に伝搬)
    name : str
        格子名

    Examples
    --------
    >>> from pysad.elements import make_element
    >>> from pysad.lattice import Lattice
    >>> elems = [make_element('QUAD','QF',L=0.3,K1=1.2),
    ...          make_element('DRIFT','D',L=2.0),
    ...          make_element('QUAD','QD',L=0.3,K1=-1.2),
    ...          make_element('DRIFT','D',L=2.0)]
    >>> lat = Lattice(elems)
    >>> tw = lat.twiss()
    >>> print(tw[0]['bx'])
    """

    def __init__(self, elements: List[Element], name: str = 'LATTICE'):
        self.elements = list(elements)
        self.name = name
        self._twiss_cache: Optional[List[dict]] = None
        # 直前の閉軌道 (SAD の CALC 間で引き継がれる初期値)
        self._cod_prev: Optional[np.ndarray] = None
        # 直前の CELL CALC で先頭 MARK に書き戻した周期解 (SAD twmov)。SAD キーワード名 → 値
        self.cell_initial: Optional[Dict[str, float]] = None
        # SAD の先頭 MARK に残っている Twiss (twmov)。不安定な CELL CALC の続行 (sad_stab) で、不安定な面の分散や
        # 和共鳴のときの Twiss の「以前の値」として使う (SAD は twiss(lb,idp,·) を上書きしないので古い値が残る)
        self._mark_init: Optional[dict] = None
        self._qstat: Optional[dict] = None        # 直近の twiss(sad_stab=..) の安定判定
        self._last_cod_periodic: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # 格子情報
    # ------------------------------------------------------------------

    @property
    def circumference(self) -> float:
        """周長 [m]"""
        return sum(el.L for el in self.elements)

    @property
    def nelem(self) -> int:
        return len(self.elements)

    def positions(self) -> np.ndarray:
        """各要素出口の s 座標 [m]、長さ nelem+1 (先頭は 0)"""
        pos = np.zeros(len(self.elements) + 1)
        for i, el in enumerate(self.elements):
            pos[i+1] = pos[i] + el.L
        return pos

    # ------------------------------------------------------------------
    # 転送行列
    # ------------------------------------------------------------------

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """全格子の 6×6 一周転送行列を返す。"""
        return _transport_6x6(self.elements, delta)

    def partial_transfer_matrix(self, i_from: int, i_to: int,
                                delta: float = 0.0) -> np.ndarray:
        """i_from 番目の要素出口から i_to 番目の要素出口までの転送行列。"""
        return _transport_6x6(self.elements[i_from:i_to], delta)

    # ------------------------------------------------------------------
    # Twiss パラメータ
    # ------------------------------------------------------------------

    def twiss(self, delta: float = 0.0,
              full: bool = True, sad_stab: Optional[dict] = None, fixed_cod: Optional[np.ndarray] = None) -> List[dict]:
        """
        各要素出口の Twiss パラメータを計算する。
        SAD の Twiss[] コマンドに相当。

        Parameters
        ----------
        delta : float  運動量偏差 δ = Δp/p
        full  : bool   全要素位置で計算 (True) か、入口のみ (False)
        sad_stab : dict or None
            None (既定): 従来どおり。不安定 (|cos μ|>1 等) なら RuntimeError。
            dict: SAD の qcell1 と同じく不安定でも続行する (|cos μ|>1 は cos μ → 1/cos μ に置き換えて Twiss を作り、
            INTRES/HALFRES が真ならその面を不安定とする)。キー intres, halfres (既定 True)、init (以前の Twiss = 先頭 MARK の値)。
            判定結果は self._qstat (stabx, staby, over, codfnd, tracex, tracey, msgs) に残す

        Returns
        -------
        twiss_list : list of dict
            長さ nelem+1 のリスト。各要素は Twiss パラメータ辞書。
            辞書キー: ax, bx, nx, ay, by, ny, ex, epx, ey, epy,
                      r1, r2, r3, r4, detr

        Raises
        ------
        RuntimeError : 格子が不安定な場合
        """
        # SAD qcod → qtwiss の流れを再現
        # _qcod_periodic() で周期閉軌道を求め、COD依存一回転行列を構築
        # fixed_cod: 閉軌道を探さず与えた軌道 (fam 点: 親の閉軌道 + 振幅) のまわりの一周行列で周期解を作る (qcod の codfnd=.true.)。
        # 閉軌道の保存や先頭 MARK への書き戻しはしない
        cod_periodic, M_cod, el_mats = _qcod_periodic(
            self.elements, delta, cod_init=self._cod_prev if fixed_cod is None else fixed_cod, fixed=fixed_cod is not None)
        codfnd = _LAST_QCOD_FOUND
        if _LAST_QCOD_OVER and sad_stab is not None:
            raise RuntimeError("格子が不安定です: 一周行列の跡が NaN になりました (SAD: qcod Overflow & closed orbit not found)。")
        if not (np.all(np.isfinite(M_cod)) and np.all(np.isfinite(cod_periodic))):
            raise RuntimeError("格子が不安定です: 閉軌道が発散しました (SAD: Overflow)。")
        sad = sad_stab is not None
        if sad and _LAST_MMAX > _MOVERFLOW:
            # SAD の COD には 1e-18 程度の丸め雑音があり、一周行列の成分が 1e18 を超えると一周で非線形領域に入って NaN になり
            # qcod の over (`*****qcod---> Overflow & closed orbit not found`) になる。pysad の COD は厳密に 0 で収束するので、
            # 同じ結果 (Twiss は前回の値のまま) にするためここで overflow とする。丸めの再現ではなく目安 (bench_level24d の D27 など)
            raise RuntimeError("格子が不安定です: 一周行列の成分が 1e18 を超えました (SAD: qcod Overflow)。")
        qinfo = None                        # SAD 流の不安定続行 (qcell1) の判定。None = 従来の経路 (周期解が安定)
        if sad:
            twiss0 = None
            try:
                with np.errstate(all='ignore'):
                    _R, _A, _B = periodic_r(M_cod)
                    _cx = 0.5 * float(np.trace(_A))
                    _cy = 0.5 * float(np.trace(_B))
                if abs(_cx) <= 1.0 and abs(_cy) <= 1.0 and _qmdiag(M_cod)[5]:
                    twiss0 = periodic_optics(M_cod)         # 安定: 従来の式 (ビット単位で同じ結果)
            except (ValueError, FloatingPointError, ZeroDivisionError):
                twiss0 = None
            if twiss0 is None:
                init = sad_stab.get('init')
                if init is None:
                    init = self._mark_init
                twiss0, qinfo = _qcell_periodic_sad(M_cod, bool(sad_stab.get('intres', True)),
                                                    bool(sad_stab.get('halfres', True)), init)
                if qinfo['nanq']:
                    raise RuntimeError("格子が不安定です: 一周行列が数ではありません (SAD: qmdiag Matrix not a number)。")
        else:
            twiss0 = periodic_optics(M_cod)
            if twiss0 is None:
                raise RuntimeError(
                    "格子が不安定です: COD依存一回転行列の固有値が虚数単位上にありません。"
                )
        # SAD は直前の CALC の閉軌道を次の初期値に使う (tffscalc: twmov)。
        # 従来は安定なときだけ保存 (ここまで来れば安定)。SAD 流の続行 (sad_stab) では、安定と判定されたときだけ
        # 保存する (SAD は全 idp が安定のときだけ twiss(1,0,·) を先頭 MARK に書き戻す。不安定なら次の初期値は前の安定解のまま)
        if not sad and fixed_cod is None:
            self._cod_prev = cod_periodic.copy()
        self._last_el_mats = el_mats
        self._last_cod_exit = None if _LAST_COD_EXIT is None else _LAST_COD_EXIT.copy()
        self._last_cod_periodic = cod_periodic.copy()
        self._cell_twiss0 = twiss0          # 書き戻し用 (分散は下で入れる)

        if qinfo is None:
            # 分散: qcell.f qcell1 の直訳 ("Dispersion is defined in 2*2 world")。
            #   結合行列 R (r1..r4) と c = sqrt(1-det R) から s 行列を作り、一周行列を
            #   対角化した A, B ブロックと、δ 列を s で変換した (a13,a23), (b13,b23) から
            #   ex  = 0.5*(a13 + a12*a23 - a22*a13)/dcosμx,  epx = 0.5*(a23 + a21*a13 - a11*a23)/dcosμx
            #   (ey, epy も B ブロックで同様。結合がある場合はモードの分散)。x_cod/δ ではない。
            tm = M_cod
            r1 = float(twiss0.get('r1', 0.0)); r2 = float(twiss0.get('r2', 0.0))
            r3 = float(twiss0.get('r3', 0.0)); r4 = float(twiss0.get('r4', 0.0))
            c1 = math.sqrt(max(0.0, 1.0 - (r1 * r4 - r2 * r3)))
            s11 = c1; s13 = -r4; s14 = r2
            s22 = c1; s23 = r3; s24 = -r1
            s31 = r1; s32 = r2; s33 = c1
            s41 = r3; s42 = r4; s44 = c1
            (tm11, tm12, tm13, tm14), (tm21, tm22, tm23, tm24), \
                (tm31, tm32, tm33, tm34), (tm41, tm42, tm43, tm44) = tm[:4, :4]
            tm15, tm25, tm35, tm45 = tm[0, 5], tm[1, 5], tm[2, 5], tm[3, 5]
            A11 = (s11 * (tm11 * s11 - tm13 * s31 - tm14 * s41)
                   + s13 * (tm31 * s11 - tm33 * s31 - tm34 * s41)
                   + s14 * (tm41 * s11 - tm43 * s31 - tm44 * s41))
            A12 = (s11 * (tm12 * s22 - tm13 * s32 - tm14 * s42)
                   + s13 * (tm32 * s22 - tm33 * s32 - tm34 * s42)
                   + s14 * (tm42 * s22 - tm43 * s32 - tm44 * s42))
            A21 = (s22 * (tm21 * s11 - tm23 * s31 - tm24 * s41)
                   + s23 * (tm31 * s11 - tm33 * s31 - tm34 * s41)
                   + s24 * (tm41 * s11 - tm43 * s31 - tm44 * s41))
            A22 = (s22 * (tm22 * s22 - tm23 * s32 - tm24 * s42)
                   + s23 * (tm32 * s22 - tm33 * s32 - tm34 * s42)
                   + s24 * (tm42 * s22 - tm43 * s32 - tm44 * s42))
            B11 = (s31 * (-tm11 * s13 - tm12 * s23 + tm13 * s33)
                   + s32 * (-tm21 * s13 - tm22 * s23 + tm23 * s33)
                   + s33 * (-tm31 * s13 - tm32 * s23 + tm33 * s33))
            B12 = (s31 * (-tm11 * s14 - tm12 * s24 + tm14 * s44)
                   + s32 * (-tm21 * s14 - tm22 * s24 + tm24 * s44)
                   + s33 * (-tm31 * s14 - tm32 * s24 + tm34 * s44))
            B21 = (s41 * (-tm11 * s13 - tm12 * s23 + tm13 * s33)
                   + s42 * (-tm21 * s13 - tm22 * s23 + tm23 * s33)
                   + s44 * (-tm41 * s13 - tm42 * s23 + tm43 * s33))
            B22 = (s41 * (-tm11 * s14 - tm12 * s24 + tm14 * s44)
                   + s42 * (-tm21 * s14 - tm22 * s24 + tm24 * s44)
                   + s44 * (-tm41 * s14 - tm42 * s24 + tm44 * s44))
            a13 = s11 * tm15 + s13 * tm35 + s14 * tm45
            a23 = s22 * tm25 + s23 * tm35 + s24 * tm45
            b13 = s31 * tm15 + s32 * tm25 + s33 * tm35
            b23 = s41 * tm15 + s42 * tm25 + s44 * tm45
            for (P11, P12, P21, P22, p13, p23, kx, kpx) in (
                    (A11, A12, A21, A22, a13, a23, 'ex', 'epx'),
                    (B11, B12, B21, B22, b13, b23, 'ey', 'epy')):
                cosmu = 0.5 * (P11 + P22)
                if abs(cosmu) > 1.0:
                    cosmu = 1.0 / cosmu
                sinmu = math.copysign(math.sqrt(max(0.0, 1.0 - cosmu**2)), P12)
                dcosmu = 2.0 * math.sin(0.5 * math.atan2(sinmu, cosmu))**2
                if dcosmu != 0.0:
                    twiss0[kx] = 0.5 * (p13 + P12 * p23 - P22 * p13) / dcosmu
                    twiss0[kpx] = 0.5 * (p23 + P21 * p13 - P11 * p23) / dcosmu
                else:
                    twiss0[kx] = 0.0
                    twiss0[kpx] = 0.0

        # SAD tffscalc: CELL で安定なら周期解 (入口の Twiss と閉軌道) を先頭要素 (MARK) の値に書き戻す
        # (twmov(1,twiss,nlat,ndim,.false.))。以後の INS CALC はこの値を入口条件に使う。
        if not sad and fixed_cod is None:
            self._twmov_to_marker(twiss0, cod_periodic)

        if not full and not sad:
            return [twiss0]

        # 全要素位置で Twiss を伝搬 (el_mats は _qcod_periodic が返した各要素行列)
        tw_cur = dict(twiss0)
        tw_cur['nx'] = 0.0
        tw_cur['ny'] = 0.0
        twiss_list = []
        over = False
        if sad:
            # qtwiss1 の直訳: 次の要素に入る前に βx, βy > 0 でなければ over。以降の行は直前の正常な行のコピーにする
            twiss_list = [dict(tw_cur)]               # [0] = 入口 (位相 0)
            # numba カーネル (perfkern.twiss_chain): 入口の β>0 判定と出力の有限判定も同じ。異常が出た要素 nbad で打ち切る
            _out, _nbad = _chain(tw_cur, el_mats, False, True)
            _n_ok = len(el_mats) if _nbad < 0 else _nbad
            twiss_list.extend([dict(zip(_TW_KEYS, row)) for row in _out[:_n_ok].tolist()])
            if _nbad >= 0:
                if True:
                    over = True
                    # 現在の行 (twiss_list[-1]) が悪い。それより前の最後の正常な行 (なければ入口) のコピーで以降を埋める
                    p = len(twiss_list) - 1
                    j = 0
                    for jj in range(p - 1, -1, -1):
                        if twiss_list[jj]['bx'] > 0.0 and twiss_list[jj]['by'] > 0.0:
                            j = jj
                            break
                    good = twiss_list[j]
                    twiss_list = twiss_list[:j + 1] + [dict(good) for _ in range(len(el_mats) - j)]
            if twiss_list:
                tw_cur = twiss_list[-1]
            if False:
                pass
            tw_entrance = twiss_list[0]
            twiss_list = twiss_list[1:]
            if not over:
                # qcell1 の末尾: 入口と出口の Twiss から一周の「トレース」x11+x22 (= 2cos μ ではなく伝搬した Twiss の不一致込み)。
                # > 2 で INTRES (< -2 で HALFRES) なら不安定、そうでなければ 4/trace に置き換える
                end = twiss_list[-1]
                tr = []
                with np.errstate(all='ignore'):
                    for kb, ka, kn in (('bx', 'ax', 'nx'), ('by', 'ay', 'ny')):
                        d = end[kn] - tw_entrance[kn]
                        cd, sd = math.cos(d), math.sin(d)
                        x11 = math.sqrt(end[kb] / tw_entrance[kb]) * (cd + tw_entrance[ka] * sd)
                        x22 = math.sqrt(tw_entrance[kb] / end[kb]) * (cd - end[ka] * sd)
                        tr.append(x11 + x22)
            else:
                tr = [0.0, 0.0]
            intres = bool(sad_stab.get('intres', True))
            halfres = bool(sad_stab.get('halfres', True))
            if qinfo is None:
                stabx = staby = True
                msgs = []
                stabq = True
            else:
                stabx, staby, stabq = qinfo['stabx'], qinfo['staby'], qinfo['stab']
                msgs = [qinfo['msg']] if qinfo.get('msg') else []
            if not over:
                for i in (0, 1):
                    if tr[i] > 2.0:
                        if intres:
                            if i == 0:
                                stabx = False
                            else:
                                staby = False
                        else:
                            tr[i] = 4.0 / tr[i]
                    if tr[i] < -2.0:
                        if halfres:
                            if i == 0:
                                stabx = False
                            else:
                                staby = False
                        else:
                            tr[i] = 4.0 / tr[i]
                stabx = stabx and stabq and codfnd
                staby = staby and stabq and codfnd
            else:
                stabx = staby = False
            self._qstat = dict(stabx=bool(stabx), staby=bool(staby), over=over, codfnd=bool(codfnd),
                               tracex=float(tr[0]), tracey=float(tr[1]), msgs=msgs, stab=bool(stabq))
            if stabx and staby and not over and fixed_cod is None:
                self._cod_prev = cod_periodic.copy()
                self._mark_init = {k: float(twiss0[k]) for k in ('ax', 'bx', 'ay', 'by', 'ex', 'epx', 'ey', 'epy',
                                                                  'r1', 'r2', 'r3', 'r4', 'detr') if k in twiss0}
                self._twmov_to_marker(twiss0, cod_periodic)
            if qinfo is None:
                # 安定 (従来の式): 一周の位相の小数部を一周行列のモードに揃える (下の従来の処理と同じ)
                tw_entrance = dict(twiss0)
                tw_cur = twiss_list[-1] if twiss_list else tw_cur
                if not over:
                    for _k, _m in (('nx', 'mux'), ('ny', 'muy')):
                        _mode = float(twiss0.get(_m, tw_cur[_k])) % (2.0 * np.pi)
                        if np.isfinite(tw_cur[_k]) and np.isfinite(_mode):
                            _n = round((tw_cur[_k] - _mode) / (2.0 * np.pi))
                            _tot = 2.0 * np.pi * _n + _mode
                        else:
                            _tot = tw_cur[_k]
                        tw_entrance[_k] = _tot
                        twiss_list[-1][_k] = _tot
            twiss_list = [tw_entrance] + twiss_list
            self._twiss_cache = twiss_list
            return twiss_list
        _out, _nbad = _chain(tw_cur, el_mats, True, False)
        if _nbad < 0:
            twiss_list = [dict(zip(_TW_KEYS, row)) for row in _out.tolist()]
            if twiss_list:
                tw_cur = twiss_list[-1]
        else:
            # 非有限が出た (不安定・縮退): 例外の出方まで同じにするため pure Python の元の経路で計算し直す
            for M_el, _cod_el in el_mats:
                tw_out = _transport_twiss(tw_cur, M_el)
                twiss_list.append(tw_out)
                tw_cur = tw_out
        # 全周分の位相前進を入口 Twiss に反映
        tw_entrance = dict(twiss0)
        # 一周の位相: 小数部は一周行列のモード (periodic_optics の mux/muy, SAD qcell と同じ)、
        # 整数部だけ伝搬から取る。途中の要素の位相は _transport_twiss (SAD qmat2twiss の直訳、結合込み) の積算。
        for _k, _m in (('nx', 'mux'), ('ny', 'muy')):
            _mode = float(twiss0.get(_m, tw_cur[_k])) % (2.0 * np.pi)
            if np.isfinite(tw_cur[_k]) and np.isfinite(_mode):
                _n = round((tw_cur[_k] - _mode) / (2.0 * np.pi))
                _tot = 2.0 * np.pi * _n + _mode
            else:
                _tot = tw_cur[_k]   # 縮退した格子などで NaN のときは従来どおり
            tw_entrance[_k] = _tot
            twiss_list[-1][_k] = _tot
        twiss_list = [tw_entrance] + twiss_list

        self._twiss_cache = twiss_list
        return twiss_list

    # ------------------------------------------------------------------
    # SAD twmov: CELL の周期解を先頭 MARK に書き戻す
    # ------------------------------------------------------------------

    # SAD の MARK キーワード名 ↔ twiss0 のキー (tffscalc twmov: cmp%value(1:ntwissfun)=twiss(1,0,:))。
    # PSIX/PSIY は入口なので 0。DZ/DDP は 4D CALC では z=0, δ=DP0。
    _TWMOV_KEYS = (('AX', 'ax'), ('BX', 'bx'), ('AY', 'ay'), ('BY', 'by'), ('EX', 'ex'), ('EPX', 'epx'),
                   ('EY', 'ey'), ('EPY', 'epy'), ('R1', 'r1'), ('R2', 'r2'), ('R3', 'r3'), ('R4', 'r4'),
                   ('DETR', 'detr'))

    def _twmov_to_marker(self, twiss0: dict, cod: np.ndarray) -> None:
        """SAD tffscalc の twmov(1,twiss,nlat,ndim,.false.) の値 (入口の周期 Twiss と閉軌道) を self.cell_initial
        (SAD キーワード名の辞書) に残す。先頭が MARK でなければ何もしない (SAD は常に先頭 MARK を仮定)。
        MARK への書き戻し自体は ffscmd 層の FFSState._twmov_to_mark が行う (下のコメント参照)。"""
        if not self.elements or type(self.elements[0]).__name__ != 'Marker':
            return
        m0 = self.elements[0]
        d = {}
        for key, tk in self._TWMOV_KEYS:
            d[key] = float(twiss0.get(tk, 0.0))
        d['PSIX'] = 0.0
        d['PSIY'] = 0.0
        d['DX'] = float(cod[0]); d['DPX'] = float(cod[1])
        d['DY'] = float(cod[2]); d['DPY'] = float(cod[3])
        d['DZ'] = float(cod[4]); d['DDP'] = float(cod[5])
        # Marker オブジェクトの属性は書き換えない: Lattice.twiss() はオフモメンタムや微分用にも呼ばれるが、SAD の twmov は
        # tffscalc の本計算 (idp=0) でしか呼ばれない。MARK の値への書き戻しは ffscmd (FFSState._twmov_to_mark) が行う。
        self.cell_initial = d

    # ------------------------------------------------------------------
    # 閉軌道 (COD) 計算 — SAD qcod 方式
    # ------------------------------------------------------------------

    def tune(self, delta: float = 0.0) -> tuple:
        """
        チューン (Qx, Qy) を返す。SAD の Twiss["NX","$$$"]/2/Pi に相当。
        (10/11 系から移植。twiss() の最後の行の位相 = 一周分の位相)
        """
        tw_list = self.twiss(delta, full=True)
        return (tw_list[-1]['nx'] / (2.0 * np.pi),
                tw_list[-1]['ny'] / (2.0 * np.pi))

    def chromaticity(self, ddelta: float = 1e-4) -> tuple:
        """
        色収差 (ξx, ξy)。10/11 系から移植した pysad 独自の計算法
        (SAD の直訳ではない。SAD との照合は未実施 = 次のベンチ Level 3 の対象)。

        自然色収差: 線形一周行列のトレースの数値微分
            ξ = (1/2π)·(−1/sin 2πQ)·d(cos 2πQ)/dδ
        六極の寄与: 解析式 ±Σ K2 η β L/(4π) を加算 (線形行列に K2 が現れないため)
        """
        from .elements import Sextupole
        M_p = np.eye(6)
        M_m = np.eye(6)
        M_0 = np.eye(6)
        for el in self.elements:
            M_p = el.transfer_matrix(+ddelta) @ M_p
            M_m = el.transfer_matrix(-ddelta) @ M_m
            M_0 = el.transfer_matrix(0.0) @ M_0
        cos_x0 = (M_0[0, 0] + M_0[1, 1]) / 2.0
        cos_y0 = (M_0[2, 2] + M_0[3, 3]) / 2.0
        sign_sin_x = np.sign(M_0[0, 1]) if abs(M_0[0, 1]) > 1e-15 else 1.0
        sign_sin_y = np.sign(M_0[2, 3]) if abs(M_0[2, 3]) > 1e-15 else 1.0
        sin_x0 = np.sqrt(max(1e-30, 1.0 - cos_x0**2)) * sign_sin_x
        sin_y0 = np.sqrt(max(1e-30, 1.0 - cos_y0**2)) * sign_sin_y
        d_cos_x = ((M_p[0, 0] + M_p[1, 1]) - (M_m[0, 0] + M_m[1, 1])) / (2.0 * ddelta) / 2.0
        d_cos_y = ((M_p[2, 2] + M_p[3, 3]) - (M_m[2, 2] + M_m[3, 3])) / (2.0 * ddelta) / 2.0
        xi_x = (1.0 / (2.0 * np.pi)) * (-1.0 / sin_x0) * d_cos_x
        xi_y = (1.0 / (2.0 * np.pi)) * (-1.0 / sin_y0) * d_cos_y

        tw_list = self.twiss(delta=0.0, full=True)
        for i, el in enumerate(self.elements):
            if isinstance(el, Sextupole) and abs(el.K2) > 1e-15:
                tw_in  = tw_list[i]
                tw_out = tw_list[i + 1]
                bx = 0.5 * (tw_in['bx'] + tw_out['bx'])
                by = 0.5 * (tw_in['by'] + tw_out['by'])
                ex = 0.5 * (tw_in['ex'] + tw_out['ex'])  # η_x
                xi_x += el.K2 * ex * bx * el.L / (4.0 * np.pi)
                xi_y -= el.K2 * ex * by * el.L / (4.0 * np.pi)

        return xi_x, xi_y

    # ------------------------------------------------------------------
    # 放射積分 / エミッタンス
    # ------------------------------------------------------------------

    def find_cod(self, delta: float = 0.0,
                 tol: float = 1e-10,
                 max_iter: int = 20) -> np.ndarray:
        """
        閉軌道 (COD) を SAD の qcod 直訳 (_qcod_periodic) で求める。
        四極の DX/DY は tsolrote と同じく座標をずらして非線形 tquade を通す。
        初期値は直前の計算の閉軌道 (SAD と同じ)。tol, max_iter は互換用 (未使用;
        収束判定は SAD の conv 式)。

        Returns
        -------
        cod : ndarray (6,)  リング入口の閉軌道 (x, px, y, py, z, delta)
        """
        cod, _, _ = _qcod_periodic(self.elements, delta,
                                   cod_init=self._cod_prev)
        self._cod_prev = cod.copy()
        return cod

    # ------------------------------------------------------------------
    # パラメータ設定: 要素 (ファミリー) と誤差 (インスタンス) を区別する
    # ------------------------------------------------------------------

    def _instances(self, family: str) -> List[Element]:
        return [el for el in self.elements if el.name == family]

    def set_element(self, family: str, key: str, value: float) -> int:
        """
        要素 (ファミリー) のパラメータ設定。同名の全インスタンスに入る。
        SAD の Element["DX","QD"]=v、または LINE["DX","QD.1"]=v の後の
        CELL; CALC (要素値が全インスタンスに展開される) に相当。
        Returns: 設定したインスタンス数
        """
        els = self._instances(family)
        if not els:
            raise KeyError(f"要素 {family} がありません")
        for el in els:
            setattr(el, key, value)
        self._twiss_cache = None
        return len(els)

    def set_error(self, instance: str, key: str, value: float) -> Element:
        """
        誤差として 1 台だけに設定する。instance は SAD 形式 "QD.1"
        (同名要素の出現順, 1 始まり)。SAD の FFS["QD.1 DX 1E-4"] に相当。
        """
        fam, _, num = instance.rpartition('.')
        if not fam or not num.isdigit():
            raise ValueError(f"インスタンス名は 'NAME.n' 形式: {instance}")
        els = self._instances(fam)
        n = int(num)
        if not 1 <= n <= len(els):
            raise KeyError(f"{instance} がありません ({fam} は {len(els)} 台)")
        setattr(els[n - 1], key, value)
        self._twiss_cache = None
        return els[n - 1]

    # ------------------------------------------------------------------

    def twiss_table(self, delta: float = 0.0,
                    cod: Optional[np.ndarray] = None,
                    full: bool = True, sad_stab: Optional[dict] = None,
                    fixed_cod: Optional[np.ndarray] = None) -> List[TwissRow]:
        """
        SAD の twiss(nlat, -1:1, 28) に相当する Twiss テーブルを計算する。

        twiss() (SAD qcod + qtwiss) と同じ計算で、各要素出口の Twiss と
        非線形追跡による閉軌道を TwissRow として返す。
        四極の DX/DY は閉軌道・Twiss の両方に反映される。

        Parameters
        ----------
        delta : float          エネルギー偏差 δ = Δp/p
        cod   : ndarray (6,)   閉軌道探索の初期値 (None なら直前の閉軌道)
        full  : bool           全要素位置で計算
        sad_stab : dict or None  twiss() 参照。dict なら不安定でも SAD 流に続行し、判定を self._qstat に残す

        Returns
        -------
        table : List[TwissRow]  長さ nelem+1 のリスト
            table[0]   = 入口 (*) での Twiss
            table[i]   = i 番目要素出口での Twiss
        """
        if fixed_cod is not None:
            # fam 点 (MatchingAmplitude): 与えた軌道のまわり。閉軌道の初期値や memo には触らない
            return self._twiss_table_impl(delta, None, full, sad_stab, np.array(fixed_cod, dtype=float))
        if cod is not None:
            self._cod_prev = np.array(cod, dtype=float)
        # 前回と入力 (要素の属性・delta・閉軌道の初期値・安定判定の設定・MARK の値) が同じなら前回の結果を返す (性能)
        key = self._tt_key(delta, full, sad_stab)
        if key is not None:
            hit = self._tt_lookup(key)
            if hit is not None:
                return hit
            before = dict(self.__dict__)
        table = self._twiss_table_impl(delta, cod, full, sad_stab)
        if key is not None:
            self._tt_store(key, before, table)
        return table

    # ---- twiss_table の直前結果キャッシュ (1 スロット)。副作用 (self の属性、モジュールの _LAST_*) も復元する
    _TT_GLOBALS = ('_LAST_COD_EXIT', '_LAST_QCOD_FOUND', '_LAST_MMAX', '_LAST_QCOD_OVER')

    def _tt_key(self, delta, full, sad_stab):
        els = self.elements
        for el in els:
            if id(el) not in _EL_CACHE:         # 要素が非スカラーの属性を持つ (キャッシュ不可) か、まだ一度も計算されていない
                if not _tt_scalar_only(el):
                    return None
        cp = self._cod_prev
        mi = self._mark_init
        try:
            ss = None if sad_stab is None else tuple(sorted((k, repr(v)) for k, v in sad_stab.items()))
        except Exception:
            return None
        return (float(delta), bool(full), ss, None if cp is None else cp.tobytes(),
                None if mi is None else tuple(sorted(mi.items())), RFSW, _EL_CACHE_ON, len(els))

    def _tt_lookup(self, key):
        m = getattr(self, '_tt_memo', None)
        if m is None or m['key'] != key:
            return None
        els = self.elements
        if len(m['ids']) != len(els) or any(a is not b for a, b in zip(m['els'], els)):
            return None
        if m['snaps'] != [el.__dict__ for el in els]:
            return None
        for k, v in m['attrs'].items():
            setattr(self, k, _tt_copy(v))
        g = globals()
        for k, v in m['globs'].items():
            g[k] = _tt_copy(v)
        return _clone_rows(m['table'])

    def _tt_store(self, key, before, table):
        els = self.elements
        if not all(id(el) in _EL_CACHE for el in els):
            self._tt_memo = None
            return
        now = self.__dict__
        attrs = {}
        for k, v in now.items():
            if k == '_tt_memo':
                continue
            if k not in before or before[k] is not v:
                attrs[k] = _tt_copy(v)
        g = globals()
        self._tt_memo = dict(key=key, els=list(els), ids=[id(e) for e in els], snaps=[dict(el.__dict__) for el in els],
                             attrs=attrs, globs={k: _tt_copy(g[k]) for k in self._TT_GLOBALS},
                             table=_clone_rows(table))

    def _twiss_table_impl(self, delta, cod, full, sad_stab, fixed_cod=None):
        tws = self.twiss(delta, full=True, sad_stab=sad_stab, fixed_cod=fixed_cod)
        el_mats = self._last_el_mats
        cod0 = self._cod_prev if (sad_stab is None and fixed_cod is None) else self._last_cod_periodic
        # 各要素出口の COD: 次要素の入口 COD、最後は一周後 (= 周期解)
        # 最後の行 ($$$) は SAD と同じく最後に評価した一周の出口 COD (収束判定が緩いので入口 cod0 と ~sqrt(conv·β) だけずれる。
        # 例: bench_level23 の K0 の COD で SAD は $$$ の DX が ^^^ と 1e-4 (相対) 違う)。x, px, y, py だけ入れ替える
        cod_end = np.array(cod0, dtype=float)
        ce = getattr(self, '_last_cod_exit', None)
        if ce is not None:
            cod_end[:4] = ce[:4]
        cods = [m[1] for m in el_mats[1:]] + [cod_end]

        def _row(s, name, tw, c):
            return TwissRow(
                s=s, name=name,
                alfx=tw['ax'], betx=tw['bx'], mux=tw['nx'],
                alfy=tw['ay'], bety=tw['by'], muy=tw['ny'],
                etx=tw['ex'], etpx=tw['epx'],
                ety=tw.get('ey', 0.0), etpy=tw.get('epy', 0.0),
                r1=tw.get('r1', 0.0), r2=tw.get('r2', 0.0),
                r3=tw.get('r3', 0.0), r4=tw.get('r4', 0.0),
                detr=tw.get('detr', 0.0),
                x_cod=float(c[0]), px_cod=float(c[1]),
                y_cod=float(c[2]), py_cod=float(c[3]),
                z_cod=float(c[4]), delta_cod=float(c[5]),
            )

        tw0 = dict(tws[0])
        tw0['nx'] = 0.0
        tw0['ny'] = 0.0
        table = [_row(0.0, '*', tw0, cod0)]
        if not full:
            return table
        pos = 0.0
        for el, tw, c in zip(self.elements, tws[1:], cods):
            pos += el.L
            table.append(_row(pos, el.name, tw, c))
        return table

    def summary(self, delta: float = 0.0) -> dict:
        """
        格子の主要パラメータをまとめた辞書を返す。

        Returns
        -------
        summary : dict
            circumference, tune_x, tune_y, beta_max_x, beta_max_y,
            beta_min_x, beta_min_y, alpha_x, alpha_y (入口),
            chrom_x, chrom_y (chromaticity, 数値微分)
        """
        tw = self.twiss(delta)
        tw0 = tw[0]

        bx_vals = np.array([t['bx'] for t in tw])
        by_vals = np.array([t['by'] for t in tw])

        # tune
        tune_x = tw[-1]['nx'] / (2 * np.pi)
        tune_y = tw[-1]['ny'] / (2 * np.pi)

        # chromaticity (数値微分)
        ddelta = 1e-4
        try:
            tw_p = self.twiss(+ddelta, full=False)
            tw_m = self.twiss(-ddelta, full=False)
            tune_xp = periodic_optics(self.transfer_matrix(+ddelta))
            tune_xm = periodic_optics(self.transfer_matrix(-ddelta))
            chrom_x = ((tune_xp['nx'] - tune_xm['nx']) / (2*np.pi)
                       / (2 * ddelta)) if tune_xp and tune_xm else 0.0
            chrom_y = ((tune_xp['ny'] - tune_xm['ny']) / (2*np.pi)
                       / (2 * ddelta)) if tune_xp and tune_xm else 0.0
        except Exception:
            chrom_x = chrom_y = 0.0

        return {
            'name': self.name,
            'circumference': self.circumference,
            'tune_x': tune_x,
            'tune_y': tune_y,
            'beta_max_x': float(bx_vals.max()),
            'beta_min_x': float(bx_vals.min()),
            'beta_max_y': float(by_vals.max()),
            'beta_min_y': float(by_vals.min()),
            'alpha_x': tw0['ax'],
            'alpha_y': tw0['ay'],
            'dispersion_max_x': float(np.array([t['ex'] for t in tw]).max()),
            'chrom_x': chrom_x,
            'chrom_y': chrom_y,
        }

    # ------------------------------------------------------------------
    # 放射積分 / エミッタンス
    # ------------------------------------------------------------------

    def radiation_integrals(self, delta: float = 0.0) -> dict:
        """
        放射積分 I1..I5 を計算する。
        SAD の CalculateEmittance[] の計算部分に相当。

        計算手順:
        1. Twiss パラメータを全要素出口で計算
        2. 各偏向電磁石 (Bend) の放射積分寄与を積算

        Notes
        -----
        I1..I5 は設計軌道上の曲率 h = angle/L に基づく積分。
        分離関数格子 (separated-function lattice) では QUAD/SEXT の
        設計軌道上の曲率は 0 (h=0) なので I1..I5 への寄与はゼロ。

        QUAD/SEXT の粒子放射 (オフアクシス位置での実効曲率による放射) は
        トラッキング中の _radiation_kick_quad() で扱われる。
        これらの寄与を含む平衡エミッタンスは 6D 拡散行列計算 (temit) が必要。

        Returns
        -------
        ri : dict
            {'I1': .., 'I2': .., 'I3': .., 'I4': .., 'I5': ..}
            単位: I1, I2 [m^-1], I3 [m^-2], I4 [m^-1], I5 [m^-1]
        """
        tw = self.twiss(delta)  # 全要素位置での Twiss

        ri = {'I1': 0.0, 'I2': 0.0, 'I3': 0.0, 'I4': 0.0, 'I5': 0.0}

        for i, el in enumerate(self.elements):
            if not isinstance(el, Bend):
                continue
            # SAD INTGRL.f 準拠: tw_in (入口 Twiss) をそのまま渡す
            # 内部で Simpson 則 + スライス毎 Twiss 伝播を行う
            contrib = _radiation_integrals_bend(el, tw[i])
            for k in ri:
                ri[k] += contrib[k]

        return ri

    def emittance(self, energy_eV: float,
                  delta: float = 0.0) -> dict:
        """
        平衡エミッタンス・エネルギースプレッド・減衰定数を計算する。
        SAD の CalculateEmittance[] に相当。

        Parameters
        ----------
        energy_eV : float  ビームエネルギー [eV]
        delta     : float  運動量偏差 (閉軌道計算用)

        Returns
        -------
        result : dict
            emittance_x  [m·rad]   – 自然エミッタンス
            emittance_y  [m·rad]   – (垂直 = 0 as lower bound)
            energy_spread           – 相対エネルギースプレッド σ_E/E
            energy_loss  [eV]      – 1 周あたり放射エネルギー損失
            momentum_compaction     – 位相 slip 係数 α_c = I1/C
            damping_x               – 水平減衰定数 (1/τ)
            damping_y               – 垂直減衰定数
            damping_z               – 縦減衰定数
            I1..I5                 – 放射積分

        Notes
        -----
        電子 / 陽電子向け。陽子・イオンでは別途スケールが必要。
        """
        ri = self.radiation_integrals(delta)
        I1, I2, I3, I4, I5 = (ri['I1'], ri['I2'], ri['I3'],
                               ri['I4'], ri['I5'])
        C = self.circumference

        gamma_rel = energy_eV / ElectronMass   # Lorentz factor

        # 放射エネルギー損失 U0 [eV/turn]
        # U0 = C_gamma * E^4 * I2 / (2π)  (E in GeV)
        E_GeV = energy_eV * 1e-9
        U0 = Cgamma * E_GeV**4 * I2 / (2 * np.pi)  # m·GeV³ × GeV⁴ × m⁻¹ ÷ (2π)
        # Cgamma の単位: m/GeV³ → U0 は GeV → eV に変換
        U0_eV = U0 * 1e9   # eV

        # 位相 slip 係数
        alpha_c = I1 / C

        # 減衰パーティション数
        Jx = 1.0 - I4 / I2
        Jy = 1.0
        Jz = 2.0 + I4 / I2

        # 自然エミッタンス
        # ε_x = Cq * γ² * I5 / (Jx * I2)
        emit_x = Cq * gamma_rel**2 * I5 / (Jx * I2) if Jx * I2 != 0 else 0.0

        # 自然エネルギースプレッド
        # σ_E²/E² = Cq * γ² * I3 / (Jz * I2)
        sigma_E2 = Cq * gamma_rel**2 * I3 / (Jz * I2) if Jz * I2 != 0 else 0.0
        sigma_E = np.sqrt(max(sigma_E2, 0.0))

        # 減衰時間 [1/turn] (正規化)
        # 1/τ_x = (Jx * U0) / (E * T0)  → [s^-1]
        # 簡易版: 減衰定数をターン数で
        T0 = C / SpeedOfLight  # 1 周時間 [s]
        damp_x = Jx * U0_eV / (energy_eV * 2)   # 1/turn (approximation)
        damp_y = Jy * U0_eV / (energy_eV * 2)
        damp_z = Jz * U0_eV / (energy_eV * 2)

        return {
            'emittance_x':        emit_x,
            'emittance_y':        0.0,
            'energy_spread':      sigma_E,
            'energy_loss_eV':     U0_eV,
            'momentum_compaction': alpha_c,
            'damping_x':          damp_x,
            'damping_y':          damp_y,
            'damping_z':          damp_z,
            **ri,
        }

    # ------------------------------------------------------------------
    # 粒子トラッキング (線形)
    # ------------------------------------------------------------------

    def track_linear(self, coords_in: np.ndarray,
                     nturns: int = 1) -> np.ndarray:
        """
        線形転送行列による粒子トラッキング。

        Parameters
        ----------
        coords_in : ndarray, shape (6,) or (6, N)
            初期座標 (x, px, y, py, z, δ)
        nturns    : int  周回数

        Returns
        -------
        coords_out : ndarray, shape (6, nturns+1) or (6, N, nturns+1)
            各ターンの座標 (ターン 0 = 初期値を含む)
        """
        one_dim = (coords_in.ndim == 1)
        x = coords_in.copy() if not one_dim else coords_in[:, np.newaxis].copy()
        M = self.transfer_matrix()
        traj = np.zeros((6, x.shape[1], nturns + 1))
        traj[:, :, 0] = x
        for t in range(nturns):
            x = M @ x
            traj[:, :, t+1] = x
        if one_dim:
            return traj[:, 0, :]
        return traj

    def track_through(self, coords_in: np.ndarray) -> np.ndarray:
        """
        格子を 1 周して各要素出口の座標を返す (線形)。

        Parameters
        ----------
        coords_in : ndarray, shape (6,)

        Returns
        -------
        coords : ndarray, shape (nelem+1, 6)
            coords[0] = 入口、coords[i] = i 番目要素出口
        """
        x = coords_in.copy()
        traj = np.zeros((len(self.elements) + 1, 6))
        traj[0] = x
        for i, el in enumerate(self.elements):
            x = el.transfer_matrix() @ x
            traj[i+1] = x
        return traj

    # ------------------------------------------------------------------
    # 文字列表示
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (f"Lattice('{self.name}', "
                f"nelem={self.nelem}, C={self.circumference:.3f} m)")

    def __len__(self) -> int:
        return self.nelem

    def __getitem__(self, idx):
        return self.elements[idx]

    def print_twiss(self, delta: float = 0.0,
                    fmt: str = '{:>6} {:>8.4f} {:>8.4f} {:>8.4f} '
                               '{:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f}'):
        """
        Twiss パラメータを SAD 風に表示する。
        """
        tw = self.twiss(delta)
        pos = self.positions()
        header = (f"{'s[m]':>6} {'βx[m]':>8} {'αx':>8} {'ψx/2π':>8} "
                  f"{'βy[m]':>8} {'αy':>8} {'ψy/2π':>8} {'ηx[m]':>8}")
        print(header)
        print('-' * len(header))
        names = ['*'] + [el.name for el in self.elements]
        for i, (t, s, nm) in enumerate(zip(tw, pos, names)):
            print(fmt.format(
                s,
                t['bx'], t['ax'], t['nx'] / (2*np.pi),
                t['by'], t['ay'], t['ny'] / (2*np.pi),
                t['ex']
            ))


# ---------------------------------------------------------------------------
# 格子構築ユーティリティ
# ---------------------------------------------------------------------------

def fodo_cell(L_half: float = 5.0,
              K1F: float = 0.5,
              K1D: float = -0.5,
              L_bend: float = 0.0,
              angle: float = 0.0) -> Lattice:
    """
    FODO セルを生成する。

    Parameters
    ----------
    L_half  : 四極磁石間のドリフト長 (各半セル) [m]
    K1F     : 集束四極の K1 [1/m²]
    K1D     : 発散四極の K1 [1/m²]
    L_bend  : 偏向電磁石長 [m]  (0 の場合は偏向なし)
    angle   : 偏向電磁石の偏向角 [rad]

    Returns
    -------
    Lattice – FODO セル
    """
    from .elements import make_element

    L_q = 0.3  # 四極磁石の半長
    elems = []
    elems.append(make_element('QUAD', 'QF', L=L_q, K1=K1F))
    if L_bend > 0:
        L_d = (L_half - L_bend) / 2.0
        if L_d > 0:
            elems.append(make_element('DRIFT', 'D', L=L_d))
        elems.append(make_element('BEND', 'B', L=L_bend, ANGLE=angle))
        if L_d > 0:
            elems.append(make_element('DRIFT', 'D', L=L_d))
    else:
        elems.append(make_element('DRIFT', 'D', L=L_half))
    elems.append(make_element('QUAD', 'QD', L=L_q, K1=K1D))
    if L_bend > 0:
        L_d = (L_half - L_bend) / 2.0
        if L_d > 0:
            elems.append(make_element('DRIFT', 'D', L=L_d))
        elems.append(make_element('BEND', 'B', L=L_bend, ANGLE=angle))
        if L_d > 0:
            elems.append(make_element('DRIFT', 'D', L=L_d))
    else:
        elems.append(make_element('DRIFT', 'D', L=L_half))

    return Lattice(elems, name='FODO')


__all__ = ['Lattice', 'TwissRow', 'fodo_cell',
           'MFITAX', 'MFITBX', 'MFITNX', 'MFITAY', 'MFITBY', 'MFITNY',
           'MFITEX', 'MFITEPX', 'MFITEY', 'MFITEPY',
           'MFITR1', 'MFITR2', 'MFITR3', 'MFITR4', 'MFITDETR',
           'MFITDX', 'MFITDPX', 'MFITDY', 'MFITDPY', 'MFITDZ', 'MFITDDP',
           'MFITALZ', 'MFITBZ', 'MFITMUZ',
           'MFITZX', 'MFITZPX', 'MFITZY', 'MFITZPY']
