"""
pysad.optics
============
SAD の Optics.n / emit.n の計算部分を Python/NumPy に翻訳。

対応する SAD 関数:
  SymplecticJ          -> symplectic_j(n)
  SymplecticInverse    -> symplectic_inverse(M)
  NormalCoordinates    -> normal_coordinates(twiss)
  TransferMatrix       -> transfer_matrix(twiss1, twiss2, ...)
  BeamMatrix           -> beam_matrix(twiss, emitx, emity, dp, ...)
  PeriodicR            -> periodic_r(M)
  PeriodicOptics       -> periodic_optics(M)

Twiss パラメータは辞書で保持:
  {
    'ax': αx,  'bx': βx,  'nx': ψx,     # x 面
    'ay': αy,  'by': βy,  'ny': ψy,     # y 面
    'ex': ηx,  'epx': η'x,              # x 分散
    'ey': ηy,  'epy': η'y,              # y 分散
    'r1': R11, 'r2': R12, 'r3': R21, 'r4': R22,  # 結合行列要素
    'detr': det(R),                      # 結合行列行列式
    'nz': ψz,                            # 縦 (6D計算時)
  }
"""

import numpy as np
from numpy.linalg import det


# ---------------------------------------------------------------------------
# シンプレクティック行列演算 (emit.n より)
# ---------------------------------------------------------------------------

def symplectic_j(n: int) -> np.ndarray:
    """
    2n × 2n のシンプレクティック単位行列 J を返す。
    SAD の SymplecticJ[n] に相当。

    J = block_diag([[0,1],[-1,0]], [[0,1],[-1,0]], ...)

    Parameters
    ----------
    n : int  (行列の次元は 2n)

    Returns
    -------
    J : ndarray, shape (2n, 2n)
    """
    m = np.zeros((2 * n, 2 * n))
    for i in range(n):
        m[2*i,   2*i+1] =  1.0
        m[2*i+1, 2*i  ] = -1.0
    return m


def symplectic_inverse(M: np.ndarray) -> np.ndarray:
    """
    シンプレクティック行列 M の逆行列 M^{-1} = -J M^T J を返す。
    SAD の SymplecticInverse[M] に相当。

    Parameters
    ----------
    M : ndarray, shape (2n, 2n)

    Returns
    -------
    M_inv : ndarray, shape (2n, 2n)
    """
    n = M.shape[0]
    J = symplectic_j(n // 2)
    return -J @ M.T @ J


# ---------------------------------------------------------------------------
# 正規座標変換行列 (Optics.n: NormalCoordinates)
# ---------------------------------------------------------------------------

def normal_coordinates(twiss: dict) -> np.ndarray:
    """
    Twiss パラメータから正規座標変換行列 N を返す。
    座標変換: u_normal = N . u_phys
    SAD の NormalCoordinatesPrimitive に相当 (4D, 結合あり)。

    Parameters
    ----------
    twiss : dict  (ax, bx, nx, ay, by, ny, ex, epx, ey, epy,
                   r1, r2, r3, r4, detr)

    Returns
    -------
    N : ndarray, shape (6, 6)
        正規座標変換行列 (B . H . R)
    """
    ax, bx = twiss['ax'], twiss['bx']
    ay, by = twiss['ay'], twiss['by']
    ex, epx = twiss['ex'], twiss['epx']
    ey, epy = twiss['ey'], twiss['epy']
    r1, r2 = twiss.get('r1', 0.0), twiss.get('r2', 0.0)
    r3, r4 = twiss.get('r3', 0.0), twiss.get('r4', 0.0)
    detr   = twiss.get('detr', 0.0)

    # --- 結合行列 R (SAD の rmatrix) ---
    a1 = np.sqrt(max(0.0, 1.0 - detr))
    if detr < 1:  # x-dominated
        R = np.array([
            [ a1,  0,  -r4,  r2,  0,  0],
            [  0, a1,   r3, -r1,  0,  0],
            [ r1, r2,   a1,   0,  0,  0],
            [ r3, r4,    0,  a1,  0,  0],
            [  0,  0,    0,   0,  1,  0],
            [  0,  0,    0,   0,  0,  1],
        ])
    else:           # y-dominated
        R = np.array([
            [-r4, r2,  a1,  0,   0,  0],
            [ r3,-r1,   0, a1,   0,  0],
            [ a1,  0,  r1, r2,   0,  0],
            [  0, a1,  r3, r4,   0,  0],
            [  0,  0,   0,  0,   1,  0],
            [  0,  0,   0,  0,   0,  1],
        ])

    # --- Courant-Snyder 行列 B ---
    B = np.array([
        [1/np.sqrt(bx),          0,           0,           0,   0,  0],
        [ax/np.sqrt(bx), np.sqrt(bx),          0,           0,   0,  0],
        [         0,          0,  1/np.sqrt(by),           0,   0,  0],
        [         0,          0,  ay/np.sqrt(by),  np.sqrt(by),  0,  0],
        [         0,          0,           0,           0,   1,  0],
        [         0,          0,           0,           0,   0,  1],
    ])

    # --- 分散補正行列 H ---
    H = np.array([
        [1, 0, 0, 0, 0, -ex ],
        [0, 1, 0, 0, 0, -epx],
        [0, 0, 1, 0, 0, -ey ],
        [0, 0, 0, 1, 0, -epy],
        [epx, -ex, epy, -ey, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ])

    return B @ H @ R


# ---------------------------------------------------------------------------
# 転送行列 (Optics.n: TransferMatrix)
# ---------------------------------------------------------------------------

def transfer_matrix(twiss1: dict, twiss2: dict,
                    gamma1: float = 1.0, gamma2: float = 1.0,
                    symplectic: bool = False,
                    calc6d: bool = False) -> np.ndarray:
    """
    2点間の転送行列を Twiss パラメータから計算する。
    SAD の TransferMatrix[from, to] に相当。

    4D (結合なし) の場合は 4×4 行列、6D の場合は 6×6 行列。

    Parameters
    ----------
    twiss1 : dict  起点の Twiss パラメータ
    twiss2 : dict  終点の Twiss パラメータ
    gamma1, gamma2 : float  各点の γβ (エネルギー、减衰計算用)
    symplectic : bool  True なら γ の比を無視 (シンプレクティック補正なし)
    calc6d : bool  True なら 6×6 行列

    Returns
    -------
    M : ndarray, shape (4, 4) or (6, 6)
    """
    gr = 1.0 if symplectic else np.sqrt(gamma1 / gamma2)

    ax1, bx1, nx1 = twiss1['ax'], twiss1['bx'], twiss1['nx']
    ax2, bx2, nx2 = twiss2['ax'], twiss2['bx'], twiss2['nx']
    ay1, by1, ny1 = twiss1['ay'], twiss1['by'], twiss1['ny']
    ay2, by2, ny2 = twiss2['ay'], twiss2['by'], twiss2['ny']

    dpsx = nx2 - nx1
    dpsy = ny2 - ny1
    cpx, spx = np.cos(dpsx), np.sin(dpsx)
    cpy, spy = np.cos(dpsy), np.sin(dpsy)

    # --- 4×4 Courant-Snyder 転送行列 ---
    m4 = np.zeros((4, 4))
    # x-x ブロック
    m4[0, 0] =  np.sqrt(bx2 / bx1) * (cpx + ax1 * spx)
    m4[0, 1] =  np.sqrt(bx2 * bx1) * spx
    m4[1, 0] = -((ax2 - ax1) * cpx + (1 + ax2 * ax1) * spx) / np.sqrt(bx2 * bx1)
    m4[1, 1] =  np.sqrt(bx1 / bx2) * (cpx - ax2 * spx)
    # y-y ブロック
    m4[2, 2] =  np.sqrt(by2 / by1) * (cpy + ay1 * spy)
    m4[2, 3] =  np.sqrt(by2 * by1) * spy
    m4[3, 2] = -((ay2 - ay1) * cpy + (1 + ay2 * ay1) * spy) / np.sqrt(by2 * by1)
    m4[3, 3] =  np.sqrt(by1 / by2) * (cpy - ay2 * spy)

    if not calc6d:
        # 結合行列 R を考慮した補正
        N1 = _rmatrix(twiss1)
        N2 = _rmatrix(twiss2)
        M4 = symplectic_inverse(N2) @ m4 @ N1
        detr1 = twiss1.get('detr', 0.0)
        detr2 = twiss2.get('detr', 0.0)
        # 水平/垂直 dominant が異なる場合は行/列入れ替え
        if (1 - detr1) * (1 - detr2) < 0:
            M4 = M4[[2, 3, 0, 1], :][:, [2, 3, 0, 1]]
        return gr * M4

    # 6D の場合 (縦方向含む)
    nz1 = twiss1.get('nz', 0.0)
    nz2 = twiss2.get('nz', 0.0)
    dpsz = nz2 - nz1
    cpz, spz = np.cos(dpsz), np.sin(dpsz)
    N1 = normal_coordinates(twiss1)
    N2 = normal_coordinates(twiss2)
    # 正規座標での位相回転行列 (6×6)
    rot = np.zeros((6, 6))
    rot[0, 0] =  cpx; rot[0, 1] = spx
    rot[1, 0] = -spx; rot[1, 1] = cpx
    rot[2, 2] =  cpy; rot[2, 3] = spy
    rot[3, 2] = -spy; rot[3, 3] = cpy
    rot[4, 4] =  cpz; rot[4, 5] = spz
    rot[5, 4] = -spz; rot[5, 5] = cpz
    row = N1[0]   # 正規座標の最初の行ベクトル（未使用）
    M6 = symplectic_inverse(N2) @ rot @ N1
    detr1 = twiss1.get('detr', 0.0)
    detr2 = twiss2.get('detr', 0.0)
    if (1 - detr1) * (1 - detr2) < 0:
        idx = [2, 3, 0, 1, 4, 5]
        M6 = M6[idx, :][:, idx]
    return gr * M6


def _rmatrix(twiss: dict) -> np.ndarray:
    """
    Twiss から 4×4 結合補正行列 (rmatrix) を作る内部関数。
    SAD の NormalCoordinatesPrimitive 内の rmatrix に相当。
    """
    r1, r2 = twiss.get('r1', 0.0), twiss.get('r2', 0.0)
    r3, r4 = twiss.get('r3', 0.0), twiss.get('r4', 0.0)
    detr   = twiss.get('detr', 0.0)
    a1 = np.sqrt(max(0.0, 1.0 - detr))
    if detr < 1:
        return np.array([
            [ a1,  0, -r4,  r2],
            [  0, a1,  r3, -r1],
            [ r1, r2,  a1,   0],
            [ r3, r4,   0,  a1],
        ])
    else:
        return np.array([
            [-r4, r2, a1,  0],
            [ r3,-r1,  0, a1],
            [ a1,  0, r1, r2],
            [  0, a1, r3, r4],
        ])


# ---------------------------------------------------------------------------
# ビーム行列 (Optics.n: BeamMatrix)
# ---------------------------------------------------------------------------

def beam_matrix(twiss: dict,
                emitx: float,
                emity: float,
                dp: float,
                emitz: float = 0.0,
                gamma_beta: float = 1.0,
                gauss: bool = True) -> np.ndarray:
    """
    Twiss パラメータとエミッタンスからビーム行列 Σ (6×6) を返す。
    SAD の BeamMatrix[i] に相当。

    Σ = N^{-1} . diag(εx,εx,εy,εy,0,dp²/γ²) . (N^{-1})^T

    Parameters
    ----------
    twiss : dict
    emitx, emity : float  水平・垂直エミッタンス [m·rad]
    dp : float  エネルギー広がり δp/p
    emitz : float  縦方向エミッタンス (6D 計算時)
    gamma_beta : float  この点での γβ (エネルギー比、不均一マシンの場合)
    gauss : bool  True なら dp^2, False なら (dp/√3)^2 (一様分布)

    Returns
    -------
    Sigma : ndarray, shape (6, 6)
    """
    N_inv = symplectic_inverse(normal_coordinates(twiss))
    # γ の変化分 (エネルギー変化のある場合)
    d = 1.0  # gamma_beta_ref / gamma_beta ← 必要なら引数で渡す
    dp2 = dp**2 if gauss else (dp / np.sqrt(3))**2
    diag = np.array([emitx, emitx, emity, emity, 0.0, d * dp2]) * d
    D = np.diag(diag)
    return N_inv @ D @ N_inv.T


# ---------------------------------------------------------------------------
# 一周転送行列からの Twiss パラメータ計算 (Optics.n: PeriodicR, PeriodicOptics)
# ---------------------------------------------------------------------------

def periodic_r(M: np.ndarray) -> tuple:
    """
    4×4 一周転送行列 M から結合行列 R、および
    x 面行列 A・y 面行列 B を求める。
    SAD の PeriodicR[m] に相当。

    参考: Willeke & Ripken (1989) の結合分解。

    Parameters
    ----------
    M : ndarray, shape (4, 4)
        一周転送行列

    Returns
    -------
    R : ndarray (2, 2)  x-y 結合行列
    A : ndarray (2, 2)  x 面の Courant-Snyder 行列
    B : ndarray (2, 2)  y 面の Courant-Snyder 行列
    """
    P = M[0:2, 0:2]
    T = M[0:2, 2:4]
    S = M[2:4, 0:2]
    Q = M[2:4, 2:4]
    J = np.array([[0.0, 1.0], [-1.0, 0.0]])

    ds = det(S)
    st = S @ T
    pq = np.trace(P - Q)
    d  = 8.0 * ds + 4.0 * np.trace(st) + pq**2

    a  = np.sqrt((1.0 + abs(pq) / np.sqrt(d)) / 2.0)
    x  = (2.0 * a**2 - 1.0) / (a * pq)

    R  = x * (J @ T.T @ J - S)
    A  = P + (x / a) * (ds * np.eye(2) + T @ S)
    B  = Q - (x / a) * (ds * np.eye(2) + st)

    return R, A, B


def periodic_optics(M: np.ndarray,
                    closed_orbit: np.ndarray | None = None) -> dict:
    """
    4×4 一周転送行列から Twiss パラメータを計算する。
    SAD の PeriodicOptics[m] に相当。

    Parameters
    ----------
    M : ndarray, shape (4, 4)
        一周転送行列
    closed_orbit : ndarray, shape (4,), optional
        閉軌道ベクトル [x, x', y, y'] (デフォルト: 零ベクトル)

    Returns
    -------
    twiss : dict
        {ax, bx, mux, ay, by, muy, r1, r2, r3, r4, detr,
         dx, dpx, dy, dpy}
    """
    R, A, B = periodic_r(M)

    # x 面
    c = np.trace(A) / 2.0
    c = np.clip(c, -1.0, 1.0)
    s = np.sqrt((1.0 - c) * (1.0 + c)) * np.sign(A[0, 1])
    if abs(s) < 1e-15:
        raise ValueError("不安定ラティス: x 面の sin(μx) ≒ 0")
    bx  =  A[0, 1] / s
    ax  = (A[0, 0] - c) / s
    mux =  np.arctan2(s, c)

    # y 面
    c = np.trace(B) / 2.0
    c = np.clip(c, -1.0, 1.0)
    s = np.sqrt((1.0 - c) * (1.0 + c)) * np.sign(B[0, 1])
    if abs(s) < 1e-15:
        raise ValueError("不安定ラティス: y 面の sin(μy) ≒ 0")
    by  =  B[0, 1] / s
    ay  = (B[0, 0] - c) / s
    muy =  np.arctan2(s, c)

    detr = det(R)

    twiss = {
        'ax': ax, 'bx': bx, 'mux': mux,
        'ay': ay, 'by': by, 'muy': muy,
        'nx': mux,          # SAD では nx=ψx (累積位相進み)
        'ny': muy,
        'r1': R[0, 0], 'r2': R[0, 1],
        'r3': R[1, 0], 'r4': R[1, 1],
        'detr': detr,
        'ex': 0.0, 'epx': 0.0,
        'ey': 0.0, 'epy': 0.0,
    }

    # 閉軌道 (分散は別途計算)
    if closed_orbit is not None:
        twiss['cod'] = closed_orbit
    else:
        # (I - M) x = 0 から閉軌道を求める (ゼロ以外ならエラー)
        try:
            from numpy.linalg import solve
            cod = solve(np.eye(4) - M, np.zeros(4))
            twiss['cod'] = cod
        except Exception:
            twiss['cod'] = np.zeros(4)

    return twiss


# ---------------------------------------------------------------------------
# 分散関数の計算 (6×6 行列から)
# ---------------------------------------------------------------------------

def dispersion_from_6d(M6: np.ndarray) -> tuple:
    """
    6×6 一周転送行列から分散関数 (ηx, η'x, ηy, η'y) を求める。

    Parameters
    ----------
    M6 : ndarray, shape (6, 6)

    Returns
    -------
    eta : ndarray, shape (4,)  [ηx, η'x, ηy, η'y]
    """
    from numpy.linalg import solve
    M4 = M6[:4, :4]
    d4 = M6[:4, 5]
    # (I - M4) η = d4
    eta = solve(np.eye(4) - M4, d4)
    return eta


# ---------------------------------------------------------------------------
# 軌道ジオメトリ (Optics.n: GeoBase, OrbitGeo)
# ---------------------------------------------------------------------------

def geo_base(chi: np.ndarray) -> np.ndarray:
    """
    3 つの Euler 角から 3×3 回転行列を返す。
    SAD の GeoBase[{chi1, chi2, chi3}] に相当。

    Parameters
    ----------
    chi : array-like, shape (3,)  [χ1, χ2, χ3] [rad]

    Returns
    -------
    R : ndarray, shape (3, 3)
    """
    chi1, chi2, chi3 = chi
    c1, c2, c3 = np.cos(chi1), np.cos(chi2), np.cos(chi3)
    s1, s2, s3 = np.sin(chi1), np.sin(chi2), np.sin(chi3)
    return np.array([
        [c3*s1 - c1*s2*s3, -c1*c3 - s1*s2*s3,  c2*s3],
        [c1*c3*s2 + s1*s3,  c3*s1*s2 - c1*s3, -c2*c3],
        [c1*c2,             c2*s1,              s2   ],
    ])


__all__ = [
    'symplectic_j', 'symplectic_inverse',
    'normal_coordinates', 'transfer_matrix',
    'beam_matrix',
    'periodic_r', 'periodic_optics',
    'dispersion_from_6d',
    'geo_base',
]
