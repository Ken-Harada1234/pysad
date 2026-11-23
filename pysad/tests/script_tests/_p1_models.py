# -*- coding: utf-8 -*-
"""担当 P1 項目 1 の補助: Polarization の脱偏極率 (線形の式) と追跡の差 (4.2 GeV・12 MV の 12 セル FODO で 2.3 倍) の非線形の寄与を、
FLUC なし・放射の損失なしの決定論的な追跡で分解する (E の `_e_spinkick.py` の続き)。

蹴り (放射の光子 1 個 = 座標の小さな蹴り h) に対するスピン軸の応答 G = (a_B+ − a_B−)/(2h) (スピンチューンの線の複素振幅、NAFF) を、次の模型で比べる:
  M1  線形の軌道 + 要素ごとの線形のスピン軌道結合 (厳密な回転。B3 の MC と同じ。Polarization の式の導出の前提)
  M2  追跡 (非線形の軌道 + 厳密なスピン)
  Mx  非線形の軌道 (追跡) + 線形のスピン軌道結合 s' = exp(−[D (u − c)]×) R s   (軌道の非線形だけを残す)
  My  線形の軌道 + 厳密なスピン (要素ごとに追跡のスピンの写像を、線形の軌道が与える座標で評価)   (スピンの非線形だけを残す)
  My[sel]  My で、厳密なスピンの写像に渡す座標を選ぶ (sel = 'd': δ だけ、'dx': (x, px, δ)、'xy': ベータトロンだけ、…。残りは閉軌道の値)
各模型の応答の比 <|G|²>_M / <|G|²>_M1 を、同じ初期条件の組で出す。"""
import numpy as np

from pysad import emit6d, track6
import _e_spinkick as EK


def build_models(ring):
    """始点 i0 の輪の要素ごとの (el, T, R, D, c_in) (放射の損失なし、閉軌道 ring.cod の上)"""
    ctx = ring.ctx
    ph = emit6d._Phys(ctx['p'], ctx['circ'])
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    secs = []
    c = ring.cod.copy()
    try:
        for el in ctx['els']:
            st = emit6d.tturne([el], c, ph, trf0=ctx['trf0'], rfsw=True, radcod=False)
            secs.append(dict(el=el, T=st.T.copy(), R=st.srot[:, :3].copy(), D=st.srot[:, 3:9].copy(), c_in=c.copy()))
            c = st.cod.copy()
    finally:
        emit6d.POL_MODE = old
    return secs


def _el_track(ring, el, c, S):
    ctx = ring.ctx
    out, s1 = track6.track_particles([el], np.array(c), 1, ctx['p'], ctx['trf0'], True, True, False,
                                     spins=[np.array(x) for x in S], radpol=False, circ=ctx['circ'])
    return np.array(out), np.array(s1)


_SEL = {'d': [5], 'x': [0, 1], 'y': [2, 3], 'z': [4], 'dx': [0, 1, 5], 'dxy': [0, 1, 2, 3, 5], 'xy': [0, 1, 2, 3], 'all': [0, 1, 2, 3, 4, 5],
        'dz': [4, 5], 'dxz': [0, 1, 4, 5]}


def track_pair_model(ring, models, u, K, hs, nu_spin, kind, sel='all'):
    """u (6) と u + h e_δ の粒子を、模型 kind ('Mx'・'My') で K 周追い、スピンチューンの線の複素振幅 [(周波数, 振幅), …] を返す"""
    c0 = ring.cod
    pts = np.array([u] + [u + np.array([0, 0, 0, 0, 0, h]) for h in hs])
    n = len(pts)
    S = np.tile(ring.n0, (n, 1))
    zs = np.zeros((K, n), complex)
    idx = _SEL[sel]
    if kind == 'My':
        uu = pts - c0
    for k in range(K):
        if kind == 'Mx':
            c = pts.copy()
            for m in models:
                S = S @ m['R'].T
                w = (c - m['c_in']) @ m['D'].T
                S = EK._rot(-w, S)
                c, _ = _el_track(ring, m['el'], c, S)
            pts = c
        else:
            for m in models:
                c_in = m['c_in']
                cc = np.tile(c_in, (n, 1))
                cc[:, idx] = c_in[idx] + uu[:, idx]
                if type(m['el']).__name__ == 'Cavity':
                    pass
                else:
                    _o, S = _el_track(ring, m['el'], cc, S)
                uu = uu @ m['T'].T
        S = S / np.linalg.norm(S, axis=1)[:, None]
        zs[k] = ring.zeta(S)
    return [K_line(zs[:, j], nu_spin, K) for j in range(n)]


def K_line(z, nu_spin, Kt):
    return EK.line_amp(z, nu_spin, Kt)


def quad_map(ring, eps=None):
    """追跡の 1 周写像の 2 次までのテイラー展開 F(u) = M u + ½ Q(u, u) (u は閉軌道からの差)。M = ヤコビアン、Q_ijk = ∂J_ij/∂u_k を、
    閉軌道を ±eps e_k だけずらした点のヤコビアンの中心差分で作る (K64 の厳密な追跡の写像)"""
    eps = np.array([1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]) if eps is None else eps
    ctx = ring.ctx
    c0 = ring.cod
    M = EK.jac(ctx, c0)
    Q = np.zeros((6, 6, 6))
    for k in range(6):
        d = np.zeros(6)
        d[k] = eps[k]
        Jp = EK.jac(ctx, c0 + d)
        Jm = EK.jac(ctx, c0 - d)
        Q[:, :, k] = (Jp - Jm) / (2 * eps[k])
    return M, Q


def track_pair_quadmap(ring, secs, u, K, hs, nu_spin, M, Q, mask=None):
    """線形のスピン (M1 と同じ要素ごとの厳密な回転) + 軌道は 1 周ごとに 2 次の写像 F(u) = M u + ½ Q(u, u) (mask で Q の成分を選ぶ。None で全部)。
    要素の途中の軌道は線形の部分写像 (u_i = Φ_i u) で、1 周の終わりで u ← F(u_始点) に置き換える"""
    c0 = ring.cod
    uu = np.array([u - c0] + [u - c0 + np.array([0, 0, 0, 0, 0, h]) for h in hs])
    n = len(uu)
    S = np.tile(ring.n0, (n, 1))
    zs = np.zeros((K, n), complex)
    Qm = Q if mask is None else Q * mask
    for k in range(K):
        u0 = uu.copy()
        for T, R, D in secs:
            S = S @ R.T
            w = uu @ D.T
            S = EK._rot(-w, S)
            uu = uu @ T.T
        S /= np.linalg.norm(S, axis=1)[:, None]
        # 非線形の 1 周写像
        uu = u0 @ M.T + 0.5 * np.einsum('ijk,nj,nk->ni', Qm, u0, u0)
        zs[k] = ring.zeta(S)
    return [K_line(zs[:, j], nu_spin, K) for j in range(n)]
