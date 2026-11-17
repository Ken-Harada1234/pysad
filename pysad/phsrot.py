# -*- coding: utf-8 -*-
"""
PHSROT 要素 (6 次元 Twiss による位相空間の回転、K64 の phsrot.f、K. Ohmi 1994)。v2 だけで動く。

K64 では動作しない (docs/21 B38): phsinit の行列がどんな入力でも 0 になり粒子の座標が 0 になる。位相 PX PY PZ は名前表に
未登録で指定できず、dv の更新は未定義の添字 i_p0 を使う。v2 は phsrot.f と matrix.f の式から意図された写像を実装する。

  M0 = N · Rot(2π PX, 2π PY, 2π PZ) · N⁻¹           (シンプレクティック)
  M  = N · Rot · (1 - D) · N⁻¹                      (RAD ON のとき。D = [Dij] は正規座標での減衰行列)
  x ← M x + N U (√λ ⊙ ξ)                            (RAD かつ FLUC のとき。B = [Bij] = U diag(λ) Uᵀ は正規座標での 1 回あたりの拡散行列、ξ は正規乱数)

  N = H R B (matrix.f m_NtoP_set_6d): B は (AX, BX), (AY, BY), (AZ, BZ) の正規化、R は結合 R1..R4、H は分散 EX EPX EY EPY ZX ZPX ZY ZPY。
  Rot は (x, px) と (y, py) が [[cos, sin], [-sin, cos]]、(z, δ) は [[cos, -sin], [sin, cos]] (matrix.f m_prot_set_6d。z は向きが逆)。
  座標は (x, px, y, py, z, δ)、px・py は正準運動量 (p0 単位)。BX・BY・BZ が 0 以下なら 1 にする (K64 と同じ)。

K64 との違い (v2 の修正):
  - 拡散: K64 は v = U·√λ (ベクトル) を作り、ξ に成分ごとに掛けてから N を掛ける (U が対角でないと共分散が B にならない)。
    v2 は N U (√λ ⊙ ξ) (共分散が N B Nᵀ になる)。
  - dv: K64 は blist(i_p0+1)、blist(i_p0+2) (i_p0 は未定義) を p0、h0 として使う。v2 は p0、h0。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import beambeam as _bb
from .track6n import lcg

NJ = dict(cache=True, error_model='numpy')
K_PHSROT = 41


def _g(p, k, d=0.0):
    return float(p.get(k, d) or d)


def matrices(params, rad=False, fluc=False):
    """戻り値: dict(N, Ninv, rot, M0, M, dJ, G)。M0: シンプレクティック、M: 減衰つき (rad のとき)、dJ = M - M0、
    G: 拡散の平方根 (6×6、x ← x + G ξ。fluc でないとき None)"""
    p = params
    bx = _g(p, 'BX')
    by = _g(p, 'BY')
    bz = _g(p, 'BZ')
    bx = bx if bx > 0.0 else 1.0
    by = by if by > 0.0 else 1.0
    bz = bz if bz > 0.0 else 1.0
    B = _bb.m_B(_g(p, 'AX'), bx, _g(p, 'AY'), by, _g(p, 'AZ'), bz)
    R = _bb.m_R(_g(p, 'R1'), _g(p, 'R2'), _g(p, 'R3'), _g(p, 'R4'))
    H = _bb.m_H(_g(p, 'EX'), _g(p, 'EPX'), _g(p, 'EY'), _g(p, 'EPY'), _g(p, 'ZX'), _g(p, 'ZPX'), _g(p, 'ZY'), _g(p, 'ZPY'))
    N = H @ R @ B
    S = np.zeros((6, 6))
    for i in range(3):
        S[2 * i, 2 * i + 1] = 1.0
        S[2 * i + 1, 2 * i] = -1.0
    Ninv = -S @ N.T @ S                       # m_sym_inv_6d (シンプレクティック行列の逆)
    rot = np.zeros((6, 6))
    for i, k in enumerate(('PX', 'PY', 'PZ')):
        a = 2.0 * math.pi * _g(p, k)
        c, s = math.cos(a), math.sin(a)
        sg = 1.0 if i < 2 else -1.0
        rot[2 * i, 2 * i] = rot[2 * i + 1, 2 * i + 1] = c
        rot[2 * i, 2 * i + 1] = sg * s
        rot[2 * i + 1, 2 * i] = -sg * s
    M0 = N @ rot @ Ninv
    M = M0
    if rad:
        D = np.array([[_g(p, 'D%d%d' % (i + 1, j + 1)) for j in range(6)] for i in range(6)])
        M = N @ rot @ (np.eye(6) - D) @ Ninv
    G = None
    if rad and fluc:
        Bn = np.zeros((6, 6))
        for i in range(6):
            for j in range(i, 6):
                Bn[i, j] = Bn[j, i] = _g(p, 'B%d%d' % (i + 1, j + 1))
        lam, U = np.linalg.eigh(Bn)
        if lam.min() < -1e-12 * max(abs(lam).max(), 1e-300):
            raise ValueError('PHSROT: 拡散行列 B に負の固有値がある')
        G = N @ U @ np.diag(np.sqrt(np.maximum(lam, 0.0)))
    return dict(N=N, Ninv=Ninv, rot=rot, M0=M0, M=M, dJ=M - M0, G=G)


def compile_phsrot(el):
    """追跡用: (M0 と減衰つき M、拡散 G) を rad・fluc の組ごとに作っておく"""
    from . import compat
    if not compat.v2():
        raise NotImplementedError('track6n: PHSROT は K64 では行列 (PtoN など) がすべて 0 になり、粒子の座標が 0 になる (docs/21 B38)。写さない')
    m0 = matrices(el.params)
    m1 = matrices(el.params, rad=True, fluc=True)
    return (np.ascontiguousarray(m0['M0']), np.ascontiguousarray(m1['M']), np.ascontiguousarray(m1['G']),
            bool(np.any(m1['G'] != 0.0)))


@njit(**NJ)
def k_phsrot(S, n, M0, Md, G, hasg, rad, fluc, rs, p0, h0):
    """S[:, (x, px/(1+δ), y, py/(1+δ), z, δ, dv, ..)]。乱数は SAD の正規乱数 (Box-Muller、粒子ごとに 6 個 = 3 組)"""
    M = Md if rad else M0
    v = np.empty(6)
    xi = np.empty(6)
    for i in range(n):
        pn = 1.0 + S[i, 5]
        v[0] = S[i, 0]
        v[1] = S[i, 1] * pn
        v[2] = S[i, 2]
        v[3] = S[i, 3] * pn
        v[4] = S[i, 4]
        v[5] = S[i, 5]
        o = M @ v
        if rad and fluc and hasg:
            for k in range(3):
                rs[0] = lcg(rs[0])
                a = rs[0] * (math.pi / 2147483648.0)
                rs[0] = lcg(rs[0])
                r = math.sqrt(-2.0 * math.log(rs[0] * (1.0 / 4294967296.0) + 0.5))
                xi[2 * k] = r * math.cos(a)
                xi[2 * k + 1] = r * math.sin(a)
            o = o + G @ xi
        pn = 1.0 + o[5]
        S[i, 0] = o[0]
        S[i, 1] = o[1] / pn
        S[i, 2] = o[2]
        S[i, 3] = o[3] / pn
        S[i, 4] = o[4]
        S[i, 5] = o[5]
        h1 = math.sqrt(1.0 + (pn * p0)**2)
        S[i, 6] = -o[5] * (1.0 + pn) / h1 / (h1 + pn * h0) + S[i, 10]


# ---------------------------------------------------------------------------- 光学・エミッタンス
def annotate(els, momentum_eV):
    """v2: PHSROT に el._solmap = (K_PHSROT, params, momentum) と transfer_matrix を付ける (CALC の Twiss は M0、
    Emittance[] は emit())。k64 モードでは付けない (恒等のまま)"""
    from . import compat
    for el in els:
        if type(el).__name__ != 'PhsRot':
            continue
        if hasattr(el, '_solmap'):
            del el._solmap
            el.__dict__.pop('transfer_matrix', None)
            el.__dict__.pop('constant_kick', None)
        if not compat.v2():
            continue
        m = matrices(el.params)
        el._solmap = (K_PHSROT, (el.params,), momentum_eV)
        M0 = m['M0']
        el.transfer_matrix = (lambda delta=0.0, M0=M0: M0.copy())
        el.constant_kick = (lambda delta=0.0: np.zeros(6))


def step(el, cod):
    """solmap.step から: (行列, 出口 COD)。CALC (4 次元の光学) はシンプレクティックな M0"""
    M0 = matrices(el.params)['M0']
    return M0.copy(), M0 @ np.asarray(cod, dtype=float)


def emit(st, params):
    """Emittance[]・EMIT: T ← M0 T、放射ありなら D ← M0 D + dJ (T + D)、B ← M0 B M0ᵀ + G Gᵀ (emit6d._State の約束: 全体の行列は T + D)"""
    rad = bool(getattr(st, 'radiation', True))
    m = matrices(params, rad=rad, fluc=rad)
    M0, dJ = m['M0'], m['dJ']
    full = st.T + st.D
    st.cod = (M0 + dJ if st.radcod else M0) @ st.cod
    if st.want_jac:
        st.apply(M0)
        if rad and st.want_db:
            st.D = st.D + dJ @ full
            if m['G'] is not None:
                st.B = st.B + m['G'] @ m['G'].T
    st.setr0()
