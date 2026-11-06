"""turnkern.py — 4D 閉軌道 (CALC: SAD qtwiss1 mat=.true.) の 1 周を numba 1 回の呼び出しで回す。

lattice._one_turn_cod_matrix が要素ごとのパラメータを配列に詰めて (lattice._TurnPack) `turn4` を呼ぶ。
要素ごとの処理は lattice._el_cod_step と同じ (DRIFT の式、薄い多極 = elkern.k_thin、QUAD = elkern._quad_full、
BEND = elkern.k_bend、それ以外は転送行列 + 定数キック)。設計軌道上 (delta=0, cod=0, ずれなし) は
Python 版と同じく要素の線形転送行列 (T0) を使い、未計算なら Python に戻して作ってもらう (戻り値 why=1)。
カーネルで扱えない要素 (MULT、F1 付き BEND、解析ヤコビアンが作れない状態など) は why=2 で Python に戻す。
Python 側は要素 i を処理したあと i+1 から再開する。
"""
import numpy as np
from numba import njit

from .elkern import NJ, k_thin, k_bend, _quad_full, _mm

K_DRIFT, K_THIN, K_QUAD, K_BEND, K_LIN, K_PY = range(6)
NP = 16            # P の列数

# P の列 (種別ごと)
#   K_DRIFT: 0 L, 1 is_cavity (RFSW OFF の空洞: 設計軌道上でも _drift_matrix_at_cod を使う)
#   K_THIN : 0 L, 1 ak, 2 kord, 3 dx, 4 dy, 5 th
#   K_QUAD : 0 L, 1 K1, 2 fringe, 3 f1in, 4 f2in, 5 f1out, 6 f2out, 7 mfring, 8 dx, 9 dy, 10 th
#   K_BEND : 0 L, 1 ANGLE, 2 K0, 3 K1, 4 E1, 5 E2, 6 has_chg, 7 dx, 8 dy, 9 dz, 10 th, 11 dth, 12 alg, 13 phig


@njit(**NJ)
def _matvec(M, v):
    out = np.empty(6)
    for i in range(6):
        s = 0.0
        for k in range(6):
            s += M[i, k] * v[k]
        out[i] = s
    return out


@njit(**NJ)
def _eye(M):
    for i in range(6):
        for j in range(6):
            M[i, j] = 1.0 if i == j else 0.0


@njit(**NJ)
def _drift_at_cod(L, delta, pxi, pyi, M, cod, move):
    """lattice._drift_matrix_at_cod と _el_cod_step の DRIFT の cod 更新 (move=True)"""
    pr = 1.0 + delta
    a = pxi ** 2 + pyi ** 2
    if a > 0.0 and abs(pr) > 1e-15:
        arg = max(0.0, 1.0 - a / pr ** 2)
        dpz = -a / pr / (1.0 + np.sqrt(arg))
    else:
        dpz = 0.0
    pzi = pr + dpz
    if abs(pzi) < 1e-15:
        pzi = pr
    ale = L / pzi
    alz = ale / pzi ** 2
    _eye(M)
    M[0, 1] = ale + pxi ** 2 * alz
    M[0, 3] = pxi * pyi * alz
    M[0, 5] = -pr * pxi * alz
    M[2, 1] = pxi * pyi * alz
    M[2, 3] = ale + pyi ** 2 * alz
    M[2, 5] = -pr * pyi * alz
    if move:
        cod[0] += pxi * ale
        cod[2] += pyi * ale
        cod[4] += dpz * ale


@njit(**NJ)
def turn4(i0, kind, has_off, P, Mlin, Klin, T0, hasT0, cod, Mtot, mats, cods, delta, PH, S3):
    """要素 i0 から最後まで回す。戻り (i, why): i == n なら完了。why=1: 要素 i の T0 (設計軌道の線形行列) が要る、
    why=2: 要素 i は Python で処理してから i+1 で再開。cod・Mtot は in-place、mats[i]・cods[i] (入口 cod) を埋める"""
    n = kind.shape[0]
    Z66 = np.zeros((6, 6))
    Z6 = np.zeros(6)
    M = np.empty((6, 6))
    i = i0
    while i < n:
        k = kind[i]
        for j in range(6):
            cods[i, j] = cod[j]
        design = (delta == 0.0 and cod[0] == 0.0 and cod[1] == 0.0 and cod[2] == 0.0 and cod[3] == 0.0
                  and not has_off[i])
        if design:
            if k == K_LIN:
                M[:, :] = Mlin[i]
            elif k == K_DRIFT and P[i, 1] == 0.0:
                _eye(M)                              # Drift.transfer_matrix(0)
                M[0, 1] = P[i, 0]
                M[2, 3] = P[i, 0]
            elif k == K_DRIFT:
                _drift_at_cod(P[i, 0], 0.0, 0.0, 0.0, M, cod, False)   # RFSW OFF の空洞
            else:
                if not hasT0[i]:
                    return i, 1
                M[:, :] = T0[i]
            cod[:] = _matvec(M, cod)
        elif k == K_DRIFT:
            _drift_at_cod(P[i, 0], delta, cod[1], cod[3], M, cod, True)
        elif k == K_THIN:
            c = cod.copy()
            T = np.eye(6)
            rc = k_thin(c, T, Z66, Z66, S3, PH, P[i, 0], P[i, 1], int(P[i, 2]), P[i, 3], P[i, 4], P[i, 5],
                        False, False, True, False)
            if rc != 0:
                return i, 2
            cod[:] = c
            M[:, :] = T
        elif k == K_QUAD:
            c = cod.copy()
            ok, J = _quad_full(c, P[i, 0], P[i, 1], P[i, 2] != 0.0, P[i, 3], P[i, 4], P[i, 5], P[i, 6],
                               int(P[i, 7]), P[i, 8], P[i, 9], P[i, 10], False, PH, True)
            if not ok:
                return i, 2
            cod[:] = c
            M[:, :] = J
        elif k == K_BEND:
            c = cod.copy()
            T = np.eye(6)
            rc = k_bend(c, T, Z66, Z66, Z6, Z66, S3, PH, P[i, 0], P[i, 1], P[i, 2], P[i, 3], P[i, 4], P[i, 5],
                        P[i, 6] != 0.0, P[i, 7], P[i, 8], P[i, 9], P[i, 10], P[i, 11], P[i, 12], P[i, 13],
                        False, False, False, True, False, False)
            if rc != 0:
                return i, 2
            cod[:] = c
            M[:, :] = T
        elif k == K_LIN:
            M[:, :] = Mlin[i]
            cod[:] = _matvec(M, cod) + Klin[i]
        else:
            return i, 2
        mats[i] = M
        Mtot[:, :] = _mm(M, Mtot)
        i += 1
    return n, 0
