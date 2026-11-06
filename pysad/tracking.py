"""
pysad.tracking
==============
SAD の tturne / TrackParticles[] に相当する粒子トラッキング。

座標系 (SAD 規約):
  (x, px, y, py, z, δ)
  x, y  : 横方向変位 [m]
  px, py: 規格化正準運動量 p_x/p0, p_y/p0 (傾き角ではない!)
  z     : 縦方向座標 z = -c·Δt [m]
  δ     : 運動量偏差 Δp/p

重要: SAD のトラッキングは厳密非線形 (近軸近似なし)。
  dpz = sqrt(1 - px² - py²) - 1  (exact)
  drift: x → x + px·L/(1+dpz)

参照: src/tdrift.f (tdrift_free), src/tquad.f (tthin), src/tbend.f
"""

from __future__ import annotations

import numpy as np
from typing import List, Optional, Union
from dataclasses import dataclass

from .elements import (Element, Drift, Quadrupole, Bend,
                       Sextupole, Cavity, TCavity, Monitor, Marker, Multipole,
                       Octupole, Decapole, Dodecapole, Solenoid)
from .constants import SpeedOfLight, Cgamma, Cq, ElectronMass, TwoPi


# ---------------------------------------------------------------------------
# シミュレーション制御フラグ (SAD の type flagset / module ffs_flag に対応)
# ---------------------------------------------------------------------------

@dataclass
class SimulationFlags:
    """
    SAD の src/tffs.f で定義される type flagset に相当するフラグ集合。
    デフォルト値は SAD の src/initb1.f による起動時初期値と一致。

    主要フラグ
    ----------
    RFSW     : RF 有効化。OFF 時は毎ターン末 z=0 リセット＋空洞をドリフト扱い。
    RAD      : 放射制動マスタースイッチ。OFF=放射なし。
    FLUC     : 量子揺動 (デフォルト **ON**)。RAD=False の場合は無効。
    RADCOD   : 閉軌道計算での放射モード (通常トラッキングでは未使用)。
    RADTAPER : エネルギーテーパ補正 (将来実装予定)。
    TRPT     : LINE/RING モード (True=LINE モード、False=RING モード)。
    POL      : スピン偏極計算 (将来実装予定)。
    CALC6D   : 6 次元計算 (将来実装予定)。

    SAD の SADScript コマンド対応
    --------------------------------
    On[RFSW]  → flags.RFSW = True
    Off[RAD]  → flags.RAD  = False
    On[FLUC]  → flags.FLUC = True   (デフォルトから変化なし)
    On[TRPT]  → flags.TRPT = True   (LINE モード)
    """
    RFSW:     bool = False   # RF 有効化。OFF: z=0 リセット + 空洞→ドリフト扱い
    RAD:      bool = False   # 放射制動マスタースイッチ
    RADCOD:   bool = False   # 閉軌道計算での放射モード (通常トラッキングでは未使用)
    RADTAPER: bool = False   # エネルギーテーパ補正 (将来実装予定)
    FLUC:     bool = True    # 量子揺動 — デフォルト ON に注意！
    TRPT:     bool = False   # LINE/RING モード切替 (True=LINE モード)
    POL:      bool = False   # スピン偏極計算 (将来実装予定)
    CALC6D:   bool = False   # 6 次元計算 (将来実装予定)


# ---------------------------------------------------------------------------
# SAD の座標系補助関数
# ---------------------------------------------------------------------------

def pxy2dpz(px: np.ndarray, py: np.ndarray) -> np.ndarray:
    """
    正準運動量 (px, py) から縦成分偏差 dpz を計算。
    SAD の pxy2dpz() に相当。

    dpz = sqrt(1 - px² - py²) - 1

    Parameters
    ----------
    px, py : ndarray  規格化横方向運動量

    Returns
    -------
    dpz : ndarray  縦方向運動量偏差
    """
    s = np.clip(px**2 + py**2, 0.0, 0.9999)
    return np.sqrt(1.0 - s) - 1.0


# ---------------------------------------------------------------------------
# 要素別の厳密シンプレクティックトラッキング関数
# ---------------------------------------------------------------------------

def track_drift_exact(x, px, y, py, z, g, L: float) -> tuple:
    """
    ドリフト空間の厳密トラッキング (SAD tdrift_free に相当)。

    x_out  = x  + px · L / (1 + dpz)
    y_out  = y  + py · L / (1 + dpz)
    z_out  = z  + dpz · L / (1 + dpz)

    L は符号付き長さ。
    """
    dpz = pxy2dpz(px, py)
    pz1 = 1.0 + dpz
    al1 = L / pz1
    x_out  = x  + px * al1
    y_out  = y  + py * al1
    z_out  = z  + dpz * al1
    return x_out, px, y_out, py, z_out, g


def track_thin_quad(x, px, y, py, z, g, K1L: float) -> tuple:
    """
    薄レンズ四極キック。
    Δpx = -K1L/(1+δ) · x
    Δpy = +K1L/(1+δ) · y
    """
    p = 1.0 + g
    kl = K1L / p
    px_out = px - kl * x
    py_out = py + kl * y
    # 縦方向シンプレクティック補正
    z_out  = z - 0.5 * kl * (x * (px_out + px) - y * (py_out + py))
    return x, px_out, y, py_out, z_out, g


def track_thin_sext(x, px, y, py, z, g, K2L: float) -> tuple:
    """
    薄レンズ六極キック。
    Δpx = -K2L/(2(1+δ)) · (x² - y²)
    Δpy = +K2L/(1+δ) · x·y
    """
    p = 1.0 + g
    k2l = K2L / p
    dpx = -0.5 * k2l * (x**2 - y**2)
    dpy =        k2l * x * y
    px_out = px + dpx
    py_out = py + dpy
    z_out  = z - x * dpx - y * dpy  # symplectic
    return x, px_out, y, py_out, z_out, g


# ---------------------------------------------------------------------------
# RF 空洞の厳密トラッキング
# ---------------------------------------------------------------------------

def track_cavity(el: 'Cavity',
                 x, px, y, py, z, g,
                 rfsw: bool = True) -> tuple:
    """
    RF 空洞の非線形シンプレクティックトラッキング。

    縦方向キック (非線形):
        δ_out = δ + (VOLT/E0) * sin(PHI - k_rf * z)

    SAD 座標系では z = -c·Δt なので、
    基準粒子からの位相ずれ = -k_rf * z (z > 0 → 先行 → 位相が遅れ側)

    横方向: 空洞長 L がある場合はドリフトとして扱う (TM モードの横方向成分は無視)

    Parameters
    ----------
    el   : Cavity  空洞要素
    x..g : ndarray  粒子座標
    rfsw : bool  SAD の RFSW フラグ (デフォルト True)。
                 False の場合は空洞をドリフトとして扱う
                 (SAD tcav.f の ``if (rfsw)`` 分岐に対応)。

    Returns
    -------
    (x, px, y, py, z, g) : 更新後の座標
    """
    # RFSW=OFF → 空洞をドリフトとして扱う (SAD の tcav.f の動作)
    if not rfsw:
        if el.L > 0.0:
            return track_drift_exact(x, px, y, py, z, g, el.L)
        return x, px, y, py, z, g

    if el.VOLT == 0.0 or el.E0 == 0.0:
        # 電圧なし → ドリフトのみ
        if el.L > 0.0:
            return track_drift_exact(x, px, y, py, z, g, el.L)
        return x, px, y, py, z, g

    k_rf = el.k_rf
    v0   = el.v0   # VOLT / E0

    # 縦方向エネルギーキック
    #   δ += v0 * sin(PHI - k_rf * z)
    g_out = g + v0 * np.sin(el.PHI - k_rf * z)

    # 横方向: ドリフト (空洞ギャップを通過する分)
    if el.L > 0.0:
        return track_drift_exact(x, px, y, py, z, g_out, el.L)

    return x, px, y, py, z, g_out


# ---------------------------------------------------------------------------
# 進行波横偏向空洞 (TCAV) の厳密トラッキング  — SAD ttcav.f の忠実な翻訳
# ---------------------------------------------------------------------------

def track_tcav(el: 'TCavity',
               x, px, y, py, z, g) -> tuple:
    """
    進行波横偏向空洞 (TCAV) の非線形シンプレクティックトラッキング。

    SAD の ttcav.f を Python/NumPy に忠実に翻訳したもの。

    構造: 半ドリフト(L/2) → 空洞キック → 半ドリフト(L/2)

    主要物理:
      * 横方向キック: Δ(px/p₀) = −AK·sin(ph) / p₁r
          (Panofsky-Wenzel 定理に基づく x-z 結合効果の一部)
      * エネルギー変化: Δh = v·x·cos(ph),  v = AK·P0·w
      * 相対論的運動量更新: Newton 反復 2 回で 2 次精度
      * py は p1r/p2r で再スケール (エネルギー変化の影響)

    参照: SAD src/ttcav.f

    Parameters
    ----------
    el : TCavity
        TCAV 要素 (AK, FREQ, HARM, PHI, P0, H0, CHARGE, OMEGA0)
    x..g : ndarray
        粒子座標 (SAD 規約):
          x, y   : 横方向変位 [m]
          px, py : 規格化正準運動量 p_t/p₀
          z      : 縦方向座標 z = −c·Δt [m]
          g      : 運動量偏差 δ = (p − p₀)/p₀

    Returns
    -------
    (x, px, y, py, z, g) : 更新後の座標タプル

    Notes
    -----
    SAD のグローバル変数との対応:
      al    ← el.L      (要素長)
      ak    ← el.AK     (偏向パラメータ)
      harm  ← el.HARM
      phi   ← el.PHI
      freq  ← el.FREQ
      p0    ← el.P0     (基準運動量 [自然単位])
      h0    ← el.H0     (基準エネルギー [自然単位])
      charge← el.CHARGE
      omega0← el.OMEGA0
      dvfs  ≈ 0          (小さい補正項; 近似)
    """
    al     = el.L
    ak     = el.AK
    p0     = el.P0
    h0     = el.H0
    charge = el.CHARGE
    dvfs   = 0.0          # SAD のグローバル補正項; 近似 = 0

    # --- 入口ハーフドリフト ---
    if al != 0.0:
        x, px, y, py, z, g = track_drift_exact(x, px, y, py, z, g, al * 0.5)

    # --- RF 波数の計算 (SAD: w = pi2*freq/c または omega0*harm/c) ---
    if el.HARM == 0.0:
        w = TwoPi * el.FREQ / SpeedOfLight
    else:
        w = el.OMEGA0 * el.HARM / SpeedOfLight

    # --- 電荷依存位相補正 (SAD: phic = phi*sign(1.,charge)) ---
    phic = el.PHI * float(np.sign(charge))

    # --- 偏向電圧パラメータ ---
    v = ak * p0 * w    # [自然単位: 長さ^−1]

    # --- 粒子ループ (SAD do 10 i=1,np をベクトル化) ---
    dp1r = np.asarray(g, dtype=float)          # = g(i) in Fortran
    p1r  = 1.0 + dp1r                          # p/p₀
    p1   = p0 * p1r                            # 全運動量 [自然単位]

    # h1: 全エネルギー [自然単位] = sqrt(1 + p1²)
    # SAD の分岐: dv>0.1 → p2h(p1)=sqrt(1+p1²), else 近似式
    # pysad では dv を追跡しないため常に厳密式を使用
    h1 = np.sqrt(1.0 + p1**2)

    # 通過時刻 t (SAD: t=min(1e4,max(-1e4,-z*h1/p1)))
    t = np.clip(-z * h1 / p1, -1.0e4, 1.0e4)

    # RF 位相
    ph = w * t + phic

    # 横方向キック (SAD: px(i) = px(i) - ak*sin(ph)/p1r)
    px = px - ak * np.sin(ph) / p1r

    # エネルギー変化 (SAD: dh=v*x(i)*cos(ph))
    dh = v * x * np.cos(ph)
    h2 = h1 + dh

    # --- 相対論的運動量更新 (SAD Newton 反復 2 回) ---
    # a = max(-1, dh*(h1+h2)/p1²)
    a = np.maximum(-1.0, dh * (h1 + h2) / p1**2)

    # 1st 近似: 3次展開
    dpr = a * (0.5 - a * (0.125 - a * 0.0625))
    # Newton 1回目
    dpr = (dpr**2 + a) / (2.0 + 2.0 * dpr)
    # Newton 2回目 (SAD と同じ 2 回)
    dpr = (dpr**2 + a) / (2.0 + 2.0 * dpr)

    dp2r = dp1r + p1r * dpr    # 新しい δ = (p₂ − p₀)/p₀
    p2r  = 1.0 + dp2r

    # g を更新 (SAD: g(i) = dp2r)
    g_out = dp2r

    # dv の更新 (SAD: dv(i)=-(1+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs)
    # pysad では dv は独立変数として保持しないが、次のドリフトで使う
    # ため、この dv は内部計算のみ (次の drift_exact では近似を使う)

    # 横方向運動量の再スケール (エネルギー変化による)
    scale = p1r / p2r
    px_out = px * scale      # SAD: px(i)*p1r/p2r
    py_out = py * scale      # SAD: py(i)*p1r/p2r

    # z 座標の更新 (SAD: z(i) = -t*p2r*p0/h2)
    z_out = -t * p2r * p0 / h2

    # --- 出口ハーフドリフト ---
    if al != 0.0:
        x, px_out, y, py_out, z_out, g_out = track_drift_exact(
            x, px_out, y, py_out, z_out, g_out, al * 0.5)

    return x, px_out, y, py_out, z_out, g_out


# ---------------------------------------------------------------------------
# 放射光効果 (電子ビーム用)
# ---------------------------------------------------------------------------

# 電子静止質量 [GeV]
_Me_GeV = ElectronMass * 1e-9   # 0.51099895e-3 GeV


def _radiation_kick(x, px, y, py, z, g,
                    h: float, ds: float,
                    E0_GeV: float,
                    rng=None) -> tuple:
    """
    偏向電磁石内の放射光エネルギー損失キック (電子ビーム用)。

    SAD の synradcl.f の簡略版。

    確定的損失 (放射ダンピング):
        Δδ = -Cgamma * E0³ [GeV³] * (1+δ)⁴ * h² * ds

    量子励起 (確率的):
        Var(Δδ) ≈ Cgamma * Cq * E0⁵ * (1+δ)⁷ * |h|³ * ds / me²
        σ(Δδ)  = sqrt(Var(Δδ))
        Δδ     ~ N(0, σ)

    Parameters
    ----------
    h      : float  曲率 1/ρ [1/m]  (符号付き)
    ds     : float  ステップ長 [m]
    E0_GeV : float  基準エネルギー [GeV]
    rng    : numpy.random.Generator or None
             None の場合は量子励起なし (確定的損失のみ)

    Returns
    -------
    (x, px, y, py, z, g) : 更新後の座標
    """
    if E0_GeV <= 0.0 or h == 0.0:
        return x, px, y, py, z, g

    h2  = h * h
    ah3 = abs(h) ** 3

    # --- 確定的エネルギー損失 ---
    # Δδ = -Cgamma * E0³ * (1+δ)⁴ * h² * ds
    damping = Cgamma * E0_GeV**3 * (1.0 + g)**4 * h2 * ds
    g_out   = g - damping

    # --- 量子励起 ---
    if rng is not None and ah3 > 0.0:
        # Var(Δδ) = Cgamma * Cq * E0⁵ * (1+δ)⁷ * |h|³ * ds / me²  [無次元]
        # me を [GeV] 単位で:  me = 0.51099895e-3 GeV
        var   = (Cgamma * Cq / _Me_GeV**2) * E0_GeV**5 * (1.0 + g)**7 * ah3 * ds
        sigma = np.sqrt(np.maximum(var, 0.0))
        kick  = rng.normal(0.0, 1.0, size=np.shape(g)) * sigma
        g_out = g_out + kick

    return x, px, y, py, z, g_out


def radiation_loss_per_turn(E0_GeV: float, rho: float) -> float:
    """
    等磁場リングでの 1 周あたり放射エネルギー損失 U0 [GeV]。

    U0 = Cgamma * E0⁴ * I2,  I2 = 2π/ρ (等磁場)

    Parameters
    ----------
    E0_GeV : float  ビームエネルギー [GeV]
    rho    : float  曲率半径 [m]

    Returns
    -------
    U0 : float  エネルギー損失 [GeV]
    """
    I2 = 2.0 * np.pi / rho
    return Cgamma * E0_GeV**4 * I2


def damping_times(E0_GeV: float, circumference: float,
                  I2: float, I4: float) -> dict:
    """
    放射ダンピング時間 τ_x, τ_y, τ_s [turns] を計算。

    α_i = Cgamma * E0³ * I2 * (1 - I4/I2) / T0   (Sands の式)

    Parameters
    ----------
    E0_GeV        : float  ビームエネルギー [GeV]
    circumference : float  リング周長 [m]
    I2            : float  放射積分 I2 [1/m]
    I4            : float  放射積分 I4 [1/m]

    Returns
    -------
    {'tau_x', 'tau_y', 'tau_s'} : ダンピング時間 [周回数]
    """
    U0   = Cgamma * E0_GeV**4 * I2            # GeV/turn
    Jx   = 1.0 - I4 / I2 if abs(I2) > 0 else 1.0
    Jy   = 1.0
    Js   = 4.0 - Jx - Jy

    alpha_x = U0 * Jx / E0_GeV   # 1/turn
    alpha_y = U0 * Jy / E0_GeV
    alpha_s = U0 * Js / E0_GeV

    eps = 1e-30
    return {
        'tau_x': 1.0 / (alpha_x + eps),
        'tau_y': 1.0 / (alpha_y + eps),
        'tau_s': 1.0 / (alpha_s + eps),
        'U0_GeV': U0,
        'Jx': Jx, 'Jy': Jy, 'Js': Js,
    }


def track_multipole_kick(x, px, y, py, z, g, KnL: dict) -> tuple:
    """
    多極キック (Multipole 要素)。
    複素数形式: w = x - i·y,  F = -sum_n KnL[n]/(n!) · w^n
    px + i·py += F
    """
    p = 1.0 + g
    w = x + 1j * y
    F = np.zeros_like(w, dtype=complex)
    for n, knl in KnL.items():
        kn = int(n)
        if kn >= 1:
            F -= (knl / p) * w**(kn - 1)  # 規格化 (factorial is in knl)
    dpx = F.real
    dpy = -F.imag
    px_out = px + dpx
    py_out = py + dpy
    z_out  = z - x * dpx - y * dpy
    return x, px_out, y, py_out, z_out, g


def _z_update(z, px, py, ds):
    """
    z 座標のみの更新 (x,y は変更しない)。
    SAD tsolqu の Strang 分割で使用。

    dz = dpz * ds / (1 + dpz),  dpz = sqrt(1 - px² - py²) - 1
    """
    dpz = np.sqrt(np.clip(1.0 - px**2 - py**2, 0.0, 1.0)) - 1.0
    return z + dpz * ds / (1.0 + dpz)


def _quad_matrix_step(x, px, y, py, k1_eff, aln0):
    """
    四極電磁石の厳密 4×4 線形行列を 1 ステップ適用。
    SAD tsolqu (bz=0) の内部ループに相当。

    k1_eff = K1 / (1+g)  (スカラーまたは 1-D 配列)
    aln0   : ステップ長 [m]

    戻り値: x, px, y, py (z, g は変更しない)
    """
    scalar = np.ndim(k1_eff) == 0
    if scalar:
        k1_eff_arr = np.atleast_1d(np.asarray(k1_eff, dtype=float))
        x_a  = np.atleast_1d(np.asarray(x,  dtype=float))
        px_a = np.atleast_1d(np.asarray(px, dtype=float))
        y_a  = np.atleast_1d(np.asarray(y,  dtype=float))
        py_a = np.atleast_1d(np.asarray(py, dtype=float))
    else:
        k1_eff_arr = np.asarray(k1_eff, dtype=float)
        x_a  = np.asarray(x,  dtype=float)
        px_a = np.asarray(px, dtype=float)
        y_a  = np.asarray(y,  dtype=float)
        py_a = np.asarray(py, dtype=float)

    x_new  = np.empty_like(x_a)
    px_new = np.empty_like(x_a)
    y_new  = np.empty_like(x_a)
    py_new = np.empty_like(x_a)

    pos = k1_eff_arr > 0
    neg = ~pos

    # K1_eff > 0: x 収束 (cos/sin), y 発散補正 (cosh/sinh)
    if np.any(pos):
        w   = np.sqrt(k1_eff_arr[pos])
        phi = w * aln0
        cp  = np.cos(phi);   sp  = np.sin(phi)
        chp = np.cosh(phi);  shp = np.sinh(phi)
        x_new[pos]  =  cp  * x_a[pos] + (sp  / w) * px_a[pos]
        px_new[pos] = -(w * sp)        * x_a[pos] + cp           * px_a[pos]
        y_new[pos]  =  chp * y_a[pos] + (shp / w) * py_a[pos]
        py_new[pos] =  (w * shp)       * y_a[pos] + chp          * py_a[pos]

    # K1_eff <= 0: x 発散 (cosh/sinh), y 収束 (cos/sin)
    if np.any(neg):
        w   = np.sqrt(-k1_eff_arr[neg])
        phi = w * aln0
        cp  = np.cos(phi);   sp  = np.sin(phi)
        chp = np.cosh(phi);  shp = np.sinh(phi)
        x_new[neg]  =  chp * x_a[neg] + (shp / w) * px_a[neg]
        px_new[neg] =  (w * shp)       * x_a[neg] + chp          * px_a[neg]
        y_new[neg]  =  cp  * y_a[neg] + (sp  / w) * py_a[neg]
        py_new[neg] = -(w * sp)        * y_a[neg] + cp            * py_a[neg]

    if scalar:
        return (float(x_new[0]), float(px_new[0]),
                float(y_new[0]), float(py_new[0]))
    return x_new, px_new, y_new, py_new


def _radiation_kick_quad(x, px, y, py, z, g,
                         K1: float, L: float,
                         E0_GeV: float,
                         rng=None) -> tuple:
    """
    四極電磁石体内の放射エネルギー損失キック (SAD tquad krad に相当)。

    SAD の tquad.f では体内通過後に tradk(..., al, 0.d0) を呼ぶ。
    phi0=0 は設計軌道上での平均偏向角がゼロを意味し、曲率は粒子の
    横方向変位 (x, y) と四極強度 K1 から決まる実効曲率:

        h_eff_x = K1 * x / (1+g),  h_eff_y = K1 * y / (1+g)
        h_eff²  = K1² * (x² + y²) / (1+g)²

    確定的エネルギー損失:
        Δδ = −C_γ * E0³[GeV] * (1+g)⁴ * h_eff² * L
           = −C_γ * E0³ * (1+g)² * K1² * (x² + y²) * L

    量子励起 (|h_eff|³ ∝ (x²+y²)^(3/2)):
        σ(Δδ)² = C_γ * C_q / m_e² * E0⁵ * (1+g)⁷ * |h_eff|³ * L

    Parameters
    ----------
    K1     : float  四極強度 [m⁻²] (局所値)
    L      : float  四極長 [m]
    E0_GeV : float  基準エネルギー [GeV]  (0 なら放射なし)
    rng    : Generator or None  量子励起用乱数生成器

    Returns
    -------
    (x, px, y, py, z, g) : 更新後の座標 (x, px, y, py, z は変化なし)
    """
    if E0_GeV <= 0.0:
        return x, px, y, py, z, g

    p    = 1.0 + g
    r2   = x * x + y * y               # 横方向変位² [m²]

    # 実効曲率²: h² = K1² * r² / p²
    h2   = (K1 / p) ** 2 * r2

    # 確定的放射損失
    damping = Cgamma * E0_GeV**3 * p**4 * h2 * L
    g_out   = g - damping

    # 量子励起
    if rng is not None:
        ah3  = np.abs(K1 / p) ** 3 * np.maximum(r2, 0.0) ** 1.5  # |h|³ * L
        var  = (Cgamma * Cq / _Me_GeV**2) * E0_GeV**5 * p**7 * ah3 * L
        sigma = np.sqrt(np.maximum(var, 0.0))
        kick  = rng.normal(0.0, 1.0, size=np.shape(g)) * sigma
        g_out = g_out + kick

    return x, px, y, py, z, g_out


def track_thick_quad(x, px, y, py, z, g,
                     L: float, K1: float,
                     nsteps: int = 4,
                     E0_GeV: float = 0.0,
                     rng=None) -> tuple:
    """
    厚肉四極電磁石の厳密行列トラッキング。
    SAD の tsolqu() (bz=0 の場合) を忠実に再現する。

    アルゴリズム (SAD tsolqu, bz=0):
      ndiv = max(1, 1 + int(|K1*L| / 0.2))   [SADのデフォルト eps=0.2]
      aln0 = L / ndiv

      Strang 分割 (2次):
        z-only half_drift(aln0/2)
        [quad_matrix(aln0) → z-only full_drift(aln0)] × (ndiv-1)
        quad_matrix(aln0)
        z-only half_drift(aln0/2)

      x, y, px, py は厳密な cos/sin/cosh/sinh 行列で更新。
      z は行列の前後で Strang 分割による z-only 更新。
      K1_eff = K1 / (1+g) は粒子ごとに計算。

    放射光 (E0_GeV > 0 の場合):
      body 通過後に _radiation_kick_quad() を呼ぶ (SAD tquad krad に対応)。
      phi0=0 で tradk() を呼ぶ SAD の実装と同等。

    Parameters
    ----------
    E0_GeV : float  基準エネルギー [GeV]  (0 なら放射なし)
    rng    : Generator or None  量子励起用乱数生成器 (None なら確定的損失のみ)
    """
    if abs(K1) < 1e-15:
        return track_drift_exact(x, px, y, py, z, g, L)

    # SAD: ndiv = 1 + int(|K1*L|/eps), eps=0.2
    ndiv = max(1, 1 + int(abs(K1 * L) / 0.2))
    aln0 = L / ndiv

    # Strang 分割: 入口 z-only 半ドリフト
    z = _z_update(z, px, py, aln0 * 0.5)

    for i in range(ndiv):
        # 粒子ごとの有効 K1
        k1_eff = K1 / (1.0 + g)
        # 厳密行列適用 (x, px, y, py のみ更新; z, g は変えない)
        x, px, y, py = _quad_matrix_step(x, px, y, py, k1_eff, aln0)
        # 中間ステップ: 行列と行列の間の全 z ドリフト (隣り合う半ドリフトの統合)
        if i < ndiv - 1:
            z = _z_update(z, px, py, aln0)

    # Strang 分割: 出口 z-only 半ドリフト (更新後の px, py を使用)
    z = _z_update(z, px, py, aln0 * 0.5)

    # 放射光キック (SAD tquad krad: tradk(..., al, 0.d0) body 後に 1 回)
    if E0_GeV > 0.0:
        x, px, y, py, z, g = _radiation_kick_quad(
            x, px, y, py, z, g, K1, L, E0_GeV, rng)

    return x, px, y, py, z, g


def track_thick_bend(x, px, y, py, z, g,
                     L: float, angle: float, K1: float = 0.0,
                     E1: float = 0.0, E2: float = 0.0,
                     HGAP: float = 0.0, FINT: float = 0.5,
                     nsteps: int = 10,
                     E0_GeV: float = 0.0,
                     rng=None) -> tuple:
    """
    偏向電磁石 (セクターベンド) の厳密トラッキング。
    SAD の tbend.f の近似実装。

    ボディ: Leap-Frog (Verlet) シンプレクティック分割
    フリンジ: 標準端部集束
    放射光:  E0_GeV > 0 かつ rng が与えられると量子励起を含む確率的キック

    Parameters
    ----------
    E0_GeV : float  基準エネルギー [GeV]  (0 で放射光なし)
    rng    : Generator or None  量子励起用乱数生成器
    """
    if abs(L) < 1e-15:
        # SAD tbthin() の薄レンズ BEND キック
        # phi0 = angle (= 全偏向角), phib/(1+g) は規格化磁場強度
        # Δpx = phi0 - phib/(1+g)  ≈  angle * (1 - 1/(1+g))  = angle * g/(1+g)
        # ただし phib = angle (等価な integrated field), SAD: phib = phi0 と見なす
        # → Δpx = phi0 - phi0/(1+g) = phi0 * g/(1+g)
        # Δz   = -x * phi0
        p    = 1.0 + g
        phi0 = angle          # 薄レンズの全偏向角
        # SAD tbthin: px += phi0 - phib/(1+g); z -= x*phi0
        # phib は磁場から決まる量で、近似として phib = phi0 と設定
        phib = phi0
        px = px + phi0 - phib / p
        z  = z  - x * phi0
        return x, px, y, py, z, g

    h = angle / L if abs(L) > 1e-15 else 0.0  # 曲率

    # 入口端部集束 (E1 > 0 の場合)
    if abs(E1) > 1e-15:
        psi1 = HGAP * FINT * h * (1 + np.sin(E1)**2) / np.cos(E1)
        px = px + h * np.tan(E1) * x
        py = py - h * np.tan(E1 - psi1) * y

    # ボディ: Leap-Frog シンプレクティック分割
    ds = L / nsteps

    for _ in range(nsteps):
        dpz = pxy2dpz(px, py)
        pz1 = 1.0 + dpz
        p   = 1.0 + g
        h_eff = h / p

        # step 1: half drift
        al1 = 0.5 * ds / pz1
        x   = x  + px * al1
        y   = y  + py * al1
        z   = z  + dpz * al1

        # step 2: bend kick (thin element approximation per step)
        dpx_bend = h_eff * pz1 - h_eff**2 * x - K1 / p * x
        dpy_bend = K1 / p * y
        px = px + dpx_bend * ds
        py = py + dpy_bend * ds

        # step 3: half drift
        dpz = pxy2dpz(px, py)
        pz1 = 1.0 + dpz
        al1 = 0.5 * ds / pz1
        x   = x  + px * al1
        y   = y  + py * al1
        z   = z  + dpz * al1

        # step 4: 放射光キック (各ステップ中央で印加)
        if E0_GeV > 0.0:
            x, px, y, py, z, g = _radiation_kick(
                x, px, y, py, z, g, h, ds, E0_GeV, rng)

    # 出口端部集束 (E2 > 0 の場合)
    if abs(E2) > 1e-15:
        psi2 = HGAP * FINT * h * (1 + np.sin(E2)**2) / np.cos(E2)
        px = px + h * np.tan(E2) * x
        py = py - h * np.tan(E2 - psi2) * y

    return x, px, y, py, z, g



# ---------------------------------------------------------------------------
# 一般薄レンズ多極キック (SAD tthin() thin 限界に対応)
# ---------------------------------------------------------------------------

def track_thin_multipole(x, px, y, py, z, g, kord: int, akL: float) -> tuple:
    """
    SAD の tthin() 薄レンズ多極キック。

    複素表現:
        cx = (x - i·y)^kord
        akf = akL / kord!
        Δpx = -akf / (1+δ) · Re(cx)
        Δpy = -akf / (1+δ) · Im(cx)

    Parameters
    ----------
    kord : int    多極次数インデックス (1=QUAD, 2=SEXT, 3=OCTU, 4=DECA, 5=DODECA)
    akL  : float  積分強度 K_n·L (規格化は kord! で行う)

    Note
    ----
    SAD の tthin() は薄キック内で z を更新しない。
    DKD スキームの中間キックとして使うため z は変更しない。
    """
    from math import factorial
    p   = 1.0 + g
    akf = akL / factorial(kord)
    w   = x - 1j * y
    cx  = w ** kord
    dpx = -(akf / p) * cx.real
    dpy = -(akf / p) * cx.imag
    return x, px + dpx, y, py + dpy, z, g


# ---------------------------------------------------------------------------
# 厚肉多極電磁石 (Yoshida 4 次 DKD) — SAD tthin() thick 版に対応
# ---------------------------------------------------------------------------

def track_tthin(x, px, y, py, z, g, nord: int, L: float, ak: float) -> tuple:
    """
    SAD の tthin() 厚肉モード (Yoshida 4 次 Drift-Kick-Drift) に相当。

    DKD 構造 (7 ステップ):
        drift(c1·L) – kick(d1·ak·L) – drift(c2·L) – kick(d2·ak·L)
        – drift(c2·L) – kick(d1·ak·L) – drift(c1·L)

    Parameters
    ----------
    nord : int    多極次数 (4=QUAD, 6=SEXT, 8=OCTU, 10=DECA, 12=DODECA)
    L    : float  有効長 [m]
    ak   : float  多極強度 K_n [m^-(n+1)]  (per unit length)
    """
    kord = nord // 2 - 1

    if abs(L) < 1e-15:
        return track_thin_multipole(x, px, y, py, z, g, kord, ak)

    # Yoshida 4 次係数 (track_thick_quad と共通)
    c1 =  0.6756035959798288
    c2 = -0.17560359597982883
    d1 =  1.3512071919596576
    d2 = -1.7024143839193152

    x, px, y, py, z, g = track_drift_exact(x, px, y, py, z, g, c1 * L)
    x, px, y, py, z, g = track_thin_multipole(x, px, y, py, z, g, kord, d1 * ak * L)
    x, px, y, py, z, g = track_drift_exact(x, px, y, py, z, g, c2 * L)
    x, px, y, py, z, g = track_thin_multipole(x, px, y, py, z, g, kord, d2 * ak * L)
    x, px, y, py, z, g = track_drift_exact(x, px, y, py, z, g, c2 * L)
    x, px, y, py, z, g = track_thin_multipole(x, px, y, py, z, g, kord, d1 * ak * L)
    x, px, y, py, z, g = track_drift_exact(x, px, y, py, z, g, c1 * L)

    return x, px, y, py, z, g


# ---------------------------------------------------------------------------
# ソレノイドトラッキング (SAD tsol.f の tsolrot に相当)
# ---------------------------------------------------------------------------

def track_solenoid(x, px, y, py, z, g, L: float, KS: float) -> tuple:
    """
    ソレノイド磁石のトラッキング (SAD tsol.f: tsolrot の簡略版)。

    ソレノイドは (x, px, y, py) の線形回転変換であるため、
    運動量偏差 g = δ を含む有効 Ks で粒子ごとに厳密計算できる。

    転送行列:
        phi = Ks·L,  C = cos(phi),  S = sin(phi)
        [x']   [ C²     CS/Ks   CS    S²/Ks ] [x ]
        [px'] = [-Ks·CS  C²    -Ks·S²  CS   ] [px]
        [y']   [-CS    -S²/Ks   C²   -CS/Ks ] [y ]
        [py']  [ Ks·S² -CS    -Ks·CS   C²   ] [py]

    Parameters
    ----------
    KS : float  ソレノイド強度 [1/m] (= Bz / (2·Bρ))
    """
    if abs(L) < 1e-15:
        return x, px, y, py, z, g

    Ks = KS / (1.0 + g)           # delta 依存; ndarray if g is array
    ks_tiny = np.abs(Ks) < 1e-15

    phi  = np.where(ks_tiny, 0.0, Ks * L)
    C    = np.cos(phi)
    S    = np.sin(phi)
    Ks_s = np.where(ks_tiny, 1.0, Ks)   # safe divisor

    x_out  =  C**2*x  + (C*S/Ks_s)*px - C*S*y        - (S**2/Ks_s)*py
    px_out = -Ks_s*C*S*x + C**2*px    + Ks_s*S**2*y  - C*S*py
    y_out  =  C*S*x  + (S**2/Ks_s)*px + C**2*y       + (C*S/Ks_s)*py
    py_out = -Ks_s*S**2*x + C*S*px    - Ks_s*C*S*y   + C**2*py

    dpz_out = pxy2dpz(px_out, py_out)
    z_out   = z + dpz_out * L / (1.0 + dpz_out)

    # Ks ≈ 0 → ドリフト代替
    dpz_in = pxy2dpz(px, py)
    Leff   = L / (1.0 + dpz_in)
    x_out  = np.where(ks_tiny, x + px*Leff,     x_out)
    px_out = np.where(ks_tiny, px,               px_out)
    y_out  = np.where(ks_tiny, y + py*Leff,     y_out)
    py_out = np.where(ks_tiny, py,               py_out)
    z_out  = np.where(ks_tiny, z + dpz_in*Leff, z_out)

    return x_out, px_out, y_out, py_out, z_out, g


# ---------------------------------------------------------------------------
# 要素ディスパッチャ
# ---------------------------------------------------------------------------

def track_element(el: Element,
                  x, px, y, py, z, g,
                  exact: bool = True,
                  nsteps: int = 4,
                  E0_GeV: float = 0.0,
                  rng=None,
                  rfsw: bool = True) -> tuple:
    """
    要素 el を通して粒子座標を更新する。

    Parameters
    ----------
    el     : Element  磁気要素
    x..g   : ndarray  座標ベクトル (全て同じ shape)
    exact  : bool     True=厳密トラッキング, False=線形行列
    nsteps : int      分割ステップ数 (厳密モード)
    E0_GeV : float    基準エネルギー [GeV] (0 で放射光なし)
    rng    : Generator or None  量子励起乱数生成器
    rfsw   : bool     SAD の RFSW フラグ。False の場合は Cavity をドリフト扱い。

    Returns
    -------
    (x, px, y, py, z, g) : tuple of ndarray
    """
    if not exact:
        # 線形転送行列を使用
        coords = np.stack([x, px, y, py, z, g], axis=0)
        M = el.transfer_matrix(float(np.mean(g)))
        coords_out = M @ coords
        return tuple(coords_out[i] for i in range(6))

    # 厳密非線形トラッキング
    if isinstance(el, Drift):
        return track_drift_exact(x, px, y, py, z, g, el.L)

    elif isinstance(el, Quadrupole):
        if el.L == 0.0:
            return track_thin_quad(x, px, y, py, z, g, el.K1)
        else:
            return track_thick_quad(x, px, y, py, z, g,
                                    el.L, el.K1, nsteps)

    elif isinstance(el, Bend):
        return track_thick_bend(x, px, y, py, z, g,
                                el.L, el.ANGLE, el.K1,
                                el.E1, el.E2, el.HGAP, el.FINT,
                                nsteps,
                                E0_GeV=E0_GeV, rng=rng)

    elif isinstance(el, Sextupole):
        # SAD tthin() に対応: Yoshida 4 次 DKD
        return track_tthin(x, px, y, py, z, g, nord=6, L=el.L, ak=el.K2)

    elif isinstance(el, Multipole):
        return track_multipole_kick(x, px, y, py, z, g, el.KnL)

    elif isinstance(el, Octupole):
        return track_tthin(x, px, y, py, z, g, nord=8,  L=el.L, ak=el.K3)

    elif isinstance(el, Decapole):
        return track_tthin(x, px, y, py, z, g, nord=10, L=el.L, ak=el.K4)

    elif isinstance(el, Dodecapole):
        return track_tthin(x, px, y, py, z, g, nord=12, L=el.L, ak=el.K5)

    elif isinstance(el, Solenoid):
        return track_solenoid(x, px, y, py, z, g, el.L, el.KS)

    elif isinstance(el, (Monitor, Marker)):
        if el.L > 0:
            return track_drift_exact(x, px, y, py, z, g, el.L)
        return x, px, y, py, z, g

    elif isinstance(el, Cavity):
        # RF 空洞 (CAVI): 縦方向加速
        return track_cavity(el, x, px, y, py, z, g, rfsw=rfsw)

    elif isinstance(el, TCavity):
        # 横偏向空洞 (TCAV): 横方向キック + Panofsky-Wenzel エネルギー変化
        # RFSW=OFF の場合はドリフトとして扱う (SAD ttcave.f の動作に準拠)
        if not rfsw:
            if el.L > 0.0:
                return track_drift_exact(x, px, y, py, z, g, el.L)
            return x, px, y, py, z, g
        return track_tcav(el, x, px, y, py, z, g)

    else:
        # 未知要素: 転送行列で代用
        coords = np.stack([x, px, y, py, z, g], axis=0)
        M = el.transfer_matrix()
        coords_out = M @ coords
        return tuple(coords_out[i] for i in range(6))


# ---------------------------------------------------------------------------
# Lattice トラッキング
# ---------------------------------------------------------------------------

@dataclass
class TrackResult:
    """トラッキング結果コンテナ"""
    coords: np.ndarray      # shape (nturns+1, 6, nparticles)
    survived: np.ndarray    # shape (nparticles,) bool

    @property
    def x(self):   return self.coords[:, 0, :]
    @property
    def px(self):  return self.coords[:, 1, :]
    @property
    def y(self):   return self.coords[:, 2, :]
    @property
    def py(self):  return self.coords[:, 3, :]
    @property
    def z(self):   return self.coords[:, 4, :]
    @property
    def delta(self): return self.coords[:, 5, :]


def track_particles(elements: List[Element],
                    coords_in: np.ndarray,
                    nturns: int = 1024,
                    aperture: float = 0.1,
                    exact: bool = True,
                    nsteps_bend: int = 20,
                    nsteps_quad: int = 4,
                    radiation: bool = False,
                    E0_GeV: float = 0.0,
                    quantum: bool = True,
                    seed: Optional[int] = None,
                    flags: Optional['SimulationFlags'] = None) -> TrackResult:
    """
    SAD の TrackParticles[] に相当する多粒子長期トラッキング。

    Parameters
    ----------
    elements   : list of Element
    coords_in  : ndarray, shape (6,) or (6, N)
                 初期座標 (x, px, y, py, z, δ)
    nturns     : int   周回数
    aperture   : float アパーチャ半径 [m] (越えた粒子は消失)
    exact      : bool  True=厳密非線形, False=線形転送行列
    nsteps_bend: int   偏向電磁石の分割数
    nsteps_quad: int   四極の分割数
    radiation  : bool  True=放射光効果を含む (電子ビーム用)。
                       ``flags`` が与えられた場合は ``flags.RAD`` が優先。
    E0_GeV     : float 基準エネルギー [GeV]  (radiation=True の場合に必要)
    quantum    : bool  True=量子励起 (確率的) を含む, False=ダンピングのみ。
                       ``flags`` が与えられた場合は ``flags.FLUC`` が優先。
    seed       : int or None  量子励起の乱数シード
    flags      : SimulationFlags or None
                 SAD の type flagset に対応するフラグ集合。与えられた場合は
                 以下の個別パラメータを上書きする:

                 フラグ    対応パラメータ  動作
                 -------- -------------- -----------------------------------
                 RFSW     (新規)          False: 毎ターン末 z=0 リセット＋
                                                 空洞をドリフト扱い
                 RAD      radiation       False: 放射なし / True: 放射あり
                 FLUC     quantum         False: 確定的ダンピングのみ
                                         True:  量子励起あり (デフォルト)
                 TRPT     (新規)          True: LINE モード (z リセットなし)

    Returns
    -------
    result : TrackResult
        .coords  : ndarray, shape (nturns+1, 6, N)
        .survived: ndarray, shape (N,) bool

    Examples
    --------
    SAD デフォルト設定 (RFSW=OFF, RAD=OFF, FLUC=ON):

    >>> flags = SimulationFlags()          # デフォルト値で作成
    >>> result = track_particles(lat, x0, nturns=1024, flags=flags)

    RF 入 + 放射あり + 量子励起あり:

    >>> flags = SimulationFlags(RFSW=True, RAD=True, FLUC=True)
    >>> result = track_particles(lat, x0, nturns=2048, E0_GeV=6.0, flags=flags)
    """
    # ------------------------------------------------------------------
    # フラグ解決: flags が与えられた場合はそちらを優先
    # ------------------------------------------------------------------
    if flags is not None:
        _rad     = flags.RAD
        _quantum = flags.FLUC
        _rfsw    = flags.RFSW
        _trpt    = flags.TRPT
    else:
        _rad     = radiation
        _quantum = quantum
        # flags なし → 後方互換のため RFSW=ON (z リセットなし) として扱う
        _rfsw    = True
        _trpt    = False   # RING モード

    # 入力の正規化
    if coords_in.ndim == 1:
        coords_in = coords_in[:, np.newaxis]
    nparticles = coords_in.shape[1]

    coords = np.zeros((nturns + 1, 6, nparticles))
    coords[0] = coords_in
    survived = np.ones(nparticles, dtype=bool)

    x  = coords_in[0].copy()
    px = coords_in[1].copy()
    y  = coords_in[2].copy()
    py = coords_in[3].copy()
    z  = coords_in[4].copy()
    g  = coords_in[5].copy()

    # 放射光用乱数生成器
    # FLUC=False (quantum=False) → rng=None → 確定的ダンピングのみ
    # RAD=False  (_rad=False)    → 放射計算スキップ
    rng = None
    if _rad and _quantum and E0_GeV > 0.0:
        rng = np.random.default_rng(seed)

    _E0 = E0_GeV if _rad else 0.0

    for turn in range(nturns):
        for el in elements:
            # アパーチャチェック
            lost = (np.abs(x) > aperture) | (np.abs(y) > aperture)
            lost |= np.isnan(x) | np.isnan(y)
            survived[lost] = False
            # 消失粒子をゼロ固定
            x[~survived]  = 0.0
            px[~survived] = 0.0
            y[~survived]  = 0.0
            py[~survived] = 0.0
            z[~survived]  = 0.0
            g[~survived]  = 0.0

            nsteps = (nsteps_bend if isinstance(el, Bend)
                      else nsteps_quad)
            x, px, y, py, z, g = track_element(
                el, x, px, y, py, z, g,
                exact=exact, nsteps=nsteps,
                E0_GeV=_E0, rng=rng,
                rfsw=_rfsw)

        # ------------------------------------------------------------------
        # ターン末処理
        # ------------------------------------------------------------------
        # RFSW=OFF かつ RING モード (TRPT=OFF): 毎ターン末に z=0 リセット。
        # 参照: SAD src/tturn.f  ``if (.not. rfsw) z(:) = 0.d0``
        # LINE モード (TRPT=ON) では一周しないため z リセットは不要。
        if not _rfsw and not _trpt:
            z[survived] = 0.0

        coords[turn + 1, 0] = x
        coords[turn + 1, 1] = px
        coords[turn + 1, 2] = y
        coords[turn + 1, 3] = py
        coords[turn + 1, 4] = z
        coords[turn + 1, 5] = g

    return TrackResult(coords=coords, survived=survived)


def dynamic_aperture(elements: List[Element],
                     nturns: int = 1024,
                     nangles: int = 9,
                     x_max: float = 0.05,
                     dp: float = 0.0,
                     exact: bool = True) -> dict:
    """
    動的アパーチャ計算 (SAD の DynamicAperture[] に相当)。

    各方位角ごとに振幅を二分法で探索し、
    nturns 周回生き残れる最大振幅を求める。

    Parameters
    ----------
    elements : list of Element
    nturns   : int   探索周回数
    nangles  : int   方位角分割数
    x_max    : float 最大探索振幅 [m]
    dp       : float 運動量偏差
    exact    : bool  True=厳密トラッキング

    Returns
    -------
    result : dict
        'angles'     : ndarray, shape (nangles,)  方位角 [rad]
        'amplitudes' : ndarray, shape (nangles,)  動的アパーチャ [m]
        'x'          : ndarray  (x 成分)
        'y'          : ndarray  (y 成分)
    """
    angles = np.linspace(0, np.pi / 2, nangles)
    amplitudes = np.zeros(nangles)

    for ia, theta in enumerate(angles):
        lo, hi = 0.0, x_max
        for _ in range(20):  # bisection iterations
            mid = 0.5 * (lo + hi)
            x0 = mid * np.cos(theta)
            y0 = mid * np.sin(theta)
            coords = np.array([[x0], [0.0], [y0], [0.0], [0.0], [dp]])
            result = track_particles(elements, coords, nturns=nturns,
                                     aperture=x_max * 1.5, exact=exact)
            if result.survived[0]:
                lo = mid
            else:
                hi = mid
        amplitudes[ia] = lo

    return {
        'angles':     angles,
        'amplitudes': amplitudes,
        'x':          amplitudes * np.cos(angles),
        'y':          amplitudes * np.sin(angles),
    }


def tune_by_fft(track_result: TrackResult,
                particle: int = 0) -> dict:
    """
    FFT によるチューン測定 (SAD の FrequencyAnalysis に相当)。

    Parameters
    ----------
    track_result : TrackResult
    particle     : int  粒子インデックス

    Returns
    -------
    {'tune_x': Qx, 'tune_y': Qy}  (0 < Q < 0.5)
    """
    # coords shape: (nturns+1, 6, N)
    x  = track_result.coords[:, 0, particle]   # x coordinate each turn
    y  = track_result.coords[:, 2, particle]   # y coordinate each turn
    nturns = len(x) - 1

    # Hanning 窓
    win = np.hanning(nturns)
    Fx = np.abs(np.fft.rfft(x[1:] * win))
    Fy = np.abs(np.fft.rfft(y[1:] * win))
    freqs = np.fft.rfftfreq(nturns)

    tune_x = freqs[np.argmax(Fx)]
    tune_y = freqs[np.argmax(Fy)]

    return {'tune_x': tune_x, 'tune_y': tune_y}


__all__ = [
    'SimulationFlags',
    'pxy2dpz',
    'track_drift_exact', 'track_thin_quad', 'track_thin_sext',
    'track_thin_multipole', 'track_tthin', 'track_solenoid',
    'track_thick_quad', 'track_thick_bend', 'track_multipole_kick',
    'track_cavity', 'track_tcav',
    '_radiation_kick', '_radiation_kick_quad', 'radiation_loss_per_turn', 'damping_times',
    'track_element', 'track_particles', 'TrackResult',
    'dynamic_aperture', 'tune_by_fft',
]
