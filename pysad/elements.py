"""
pysad.elements
==============
SAD の各磁気要素の転送行列実装。
C/Fortran コア (libsad.a) の要素計算部分に相当。

SAD 要素タイプ:
  DRIFT  - ドリフト空間
  QUAD   - 四極電磁石
  BEND   - 偏向電磁石 (二極子)
  SEXT   - 六極電磁石
  MONI   - ビーム位置モニタ
  MARK   - マーカー
  MULT   - 多極子
  CAVI   - RF 空洞
  CORS   - 軌道補正コイル (水平/垂直)

転送行列の規則:
  相座標 (x, x', y, y', l, δ) の 6×6 行列。
  ここで δ = Δp/p (運動量偏差)、l = -c Δt (縦方向座標)。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# 基底クラス
# ---------------------------------------------------------------------------

@dataclass
class Element:
    """
    SAD 要素の基底クラス。
    全要素は transfer_matrix() メソッドで 6×6 転送行列を返す。
    """
    name: str = "MARK"
    L: float = 0.0            # 要素長 [m]
    K0: float = 0.0           # 水平薄レンズキック [rad] — COD 生成 (補正コイル等)
    ROTATE: float = 0.0       # 垂直薄レンズキック [rad] — 垂直ステアリング補正
    DX: float = 0.0           # 水平設置誤差オフセット [m] (四極→等価キック)
    DY: float = 0.0           # 垂直設置誤差オフセット [m]
    ROTATION: float = 0.0     # ロール誤差 [rad] — 据え付け回転ミスアライメント
                              # 主要要素の transfer_matrix() で適用 (モジュール末尾の _wrap_rotation)。
                              # 向きは SAD の ROT (tsolrote/tchge) と同じ: x' = cosθ·x − sinθ·y。

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """6×6 一次転送行列 (delta: 運動量偏差 δ = Δp/p)"""
        return np.eye(6)

    def _apply_rotation(self, M: np.ndarray) -> np.ndarray:
        """
        ROTATION ロール誤差を転送行列に適用する: M_rot = R(-θ) · M · R(θ)。
        R(θ) は横断面内のロール回転 (x, px, y, py)。
        """
        rotation = getattr(self, 'ROTATION', 0.0)
        if abs(rotation) < 1e-10:
            return M
        # SAD (tsolrote/tchge) の向き: 素子座標 x' = cosθ·x − sinθ·y, y' = sinθ·x + cosθ·y
        c, s = np.cos(rotation), -np.sin(rotation)
        R = np.array([[c,  0., s,  0.],
                      [0., c,  0., s ],
                      [-s, 0., c,  0.],
                      [0., -s, 0., c ]])
        R6 = np.eye(6); R6[:4, :4] = R
        Rt6 = np.eye(6); Rt6[:4, :4] = R.T
        return Rt6 @ M @ R6

    def radiation_integrals(self) -> dict:
        """放射積分への寄与 I1..I5 を返す (デフォルト: 0)"""
        return {f'I{i}': 0.0 for i in range(1, 6)}

    def constant_kick(self, delta: float = 0.0) -> np.ndarray:
        """
        定常軌道キック 6-vector (閉軌道計算 find_cod に使用)。

        K0     : 水平ステアリング補正コイルキック [rad]
        ROTATE : 垂直ステアリング補正コイルキック [rad]

        サブクラス (Quadrupole 等) は DX/DY フィードダウンもここに含める。

        Returns
        -------
        d : ndarray (6,)  [Δx, Δpx, Δy, Δpy, Δz, Δg]
        """
        d = np.zeros(6)
        d[1] += self.K0       # 水平補正コイルキック → Δpx
        d[3] += self.ROTATE   # 垂直補正コイルキック → Δpy
        return d


# ---------------------------------------------------------------------------
# ドリフト空間
# ---------------------------------------------------------------------------

@dataclass
class Drift(Element):
    """
    DRIFT 要素。
    SAD コマンド: DRIFT L=<長さ>;
    """
    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """
        qtwiss.f の DRFT 線形行列 (COD = 0 の場合)。
        delta ≠ 0 の場合は有効長 L/(1+delta) を使用 (色収差補正)。
        """
        L = self.L / (1.0 + delta)   # 有効ドリフト長 (SAD: ale = al/pzi = al/(1+delta))
        M = np.eye(6)
        M[0, 1] = L
        M[2, 3] = L
        return M


# ---------------------------------------------------------------------------
# 四極電磁石
# ---------------------------------------------------------------------------

@dataclass
class Quadrupole(Element):
    """
    QUAD 要素。
    SAD コマンド: QUAD L=<長さ> K1=<四極強度>;
    K1 > 0: x-focusing (水平集束), K1 < 0: y-focusing

    K1 は局所勾配 [m⁻²] = B'/(Bρ)。
    SAD .n ファイルからロードする場合、parser.py が SAD 積分強度 K1_SAD [m⁻¹] を
    K1_SAD / L [m⁻²] に変換してからこのクラスに渡す。
    """
    K1: float = 0.0        # 局所四極強度 [m⁻²]  (= B'/(Bρ))
    K1S: float = 0.0       # スキュー四極成分

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        L = self.L
        M = np.eye(6)
        if abs(L) < 1e-15:
            return M
        # K1 は局所勾配 [m⁻²] として格納済み
        K1 = self.K1 / (1.0 + delta)   # per-unit-length [m⁻²]
        if abs(K1) < 1e-15:
            # 四極強度 ≈ 0 → ドリフト (色収差考慮)
            Leff = L / (1.0 + delta)
            M[0, 1] = Leff
            M[2, 3] = Leff
            return M
        if K1 > 0:
            # x: 集束, y: 発散
            kx = np.sqrt(K1)
            ky = np.sqrt(K1)
            # x-x block (集束)
            M[0, 0] =  np.cos(kx * L)
            M[0, 1] =  np.sin(kx * L) / kx
            M[1, 0] = -kx * np.sin(kx * L)
            M[1, 1] =  np.cos(kx * L)
            # y-y block (発散)
            M[2, 2] =  np.cosh(ky * L)
            M[2, 3] =  np.sinh(ky * L) / ky
            M[3, 2] =  ky * np.sinh(ky * L)
            M[3, 3] =  np.cosh(ky * L)
        else:
            # x: 発散, y: 集束
            kx = np.sqrt(-K1)
            ky = np.sqrt(-K1)
            # x-x block (発散)
            M[0, 0] =  np.cosh(kx * L)
            M[0, 1] =  np.sinh(kx * L) / kx
            M[1, 0] =  kx * np.sinh(kx * L)
            M[1, 1] =  np.cosh(kx * L)
            # y-y block (集束)
            M[2, 2] =  np.cos(ky * L)
            M[2, 3] =  np.sin(ky * L) / ky
            M[3, 2] = -ky * np.sin(ky * L)
            M[3, 3] =  np.cos(ky * L)
        # --- 正準座標補正 (SAD canonical: px = p_x/p_0 = (1+delta)*x') ---
        pr = 1.0 + delta
        M[0, 1] /= pr
        M[1, 0] *= pr
        M[2, 3] /= pr
        M[3, 2] *= pr
        return M

    def constant_kick(self, delta: float = 0.0) -> np.ndarray:
        """
        四極電磁石 DX/DY 据え付け誤差オフセットによる定常フィードダウンキック。

        磁石中心が (DX, DY) だけ設計軌道からずれている場合、
        設計軌道 (x=0, y=0) を通る粒子は磁石フレームでは (-DX, -DY) の
        位置を通過するため、等価な定常四極キックを受ける。

        厚肉四極の場合、転送行列 M を用いた等価定常キック:
            d[0] = (1 − M[0,0]) * DX    (x 変位成分)
            d[1] =    −M[1,0]   * DX    (px キック成分)
            d[2] = (1 − M[2,2]) * DY    (y 変位成分)
            d[3] =    −M[3,2]   * DY    (py キック成分)

        また基底クラスの K0/ROTATE 補正コイルキックも含める。

        SAD の qcod / qtwiss に相当する閉軌道計算で使用。
        """
        d = super().constant_kick(delta)  # K0, ROTATE から
        if abs(self.DX) < 1e-15 and abs(self.DY) < 1e-15:
            return d
        M = self.transfer_matrix(delta)
        d[0] += (1.0 - M[0, 0]) * self.DX
        d[1] +=       -M[1, 0]  * self.DX
        d[2] += (1.0 - M[2, 2]) * self.DY
        d[3] +=       -M[3, 2]  * self.DY
        return d


# ---------------------------------------------------------------------------
# 偏向電磁石 (セクターベンド + 端部集束)
# ---------------------------------------------------------------------------

@dataclass
class Bend(Element):
    """
    BEND 要素 (SBEND 相当)。
    SAD コマンド: BEND L=<長さ> ANGLE=<偏向角> K1=<四極> E1=<入口角> E2=<出口角>;

    全偏向角 ANGLE = L / rho  [rad]
    K1: 四極成分 [1/m²]
    E1, E2: 入口/出口の端部角 [rad]
    HGAP: ギャップ半幅 [m]
    FINT: フリンジ積分
    """
    ANGLE: float = 0.0     # 全偏向角 [rad]
    K1: float = 0.0        # 四極成分 [1/m²]
    E1: float = 0.0        # 入口端部角 [rad]
    E2: float = 0.0        # 出口端部角 [rad]
    HGAP: float = 0.0      # ポールギャップ半幅 [m]
    FINT: float = 0.5      # フリンジ積分係数

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        L, angle = self.L, self.ANGLE
        if abs(L) < 1e-15:
            return np.eye(6)

        rho = L / angle if abs(angle) > 1e-15 else 1e15
        h   = 1.0 / rho        # 曲率

        # K1 は局所勾配 [m⁻²] = B'/(Bρ) として格納済み
        # (dp スケーリングは _sbend_body 内部で行う)

        # --- ボディ転送行列 (機械的形式) ---
        M = _sbend_body(L, h, self.K1, delta)

        # --- ボディのみ先に正準座標補正 ---
        # SAD tbendebody の trans1(1,2)=sxkx/p, trans1(2,1)=aksx*p は既に正準形。
        # 端部行列 _edge_matrix の M[1,0]=h*tan(E1) は SAD tbedge
        # (tbedgedrift+tbedgebody の合成) から導かれた正準形式で pr 係数なし。
        # tbedgedrift が x を 1/cos 伸縮し、tbedgebody が sin/rho で蹴るため
        # 合成すると h*tan (pr 不含) が得られる。
        # よってボディのみに pr スケーリングを適用し、端部には適用しない。
        pr = 1.0 + delta
        M[0, 1] /= pr
        M[1, 0] *= pr
        M[2, 3] /= pr
        M[3, 2] *= pr

        # --- 入口端部 (E1): 正準形式で右から掛ける ---
        if abs(self.E1) > 1e-15:
            M = M @ _edge_matrix(h, self.E1, self.HGAP, self.FINT)

        # --- 出口端部 (E2): 正準形式で左から掛ける ---
        if abs(self.E2) > 1e-15:
            M = _edge_matrix(h, self.E2, self.HGAP, self.FINT) @ M

        # --- ROTATE (ロール角): ビーム軸まわりの素子回転 ---
        rotate = self.ROTATE  # [rad]
        if abs(rotate) > 1e-10:
            # SAD (tchge) の向き: 素子座標 x' = cosθ·x − sinθ·y
            c, s = np.cos(rotate), -np.sin(rotate)
            R = np.array([[c, 0.,  s, 0.],
                          [0., c, 0.,  s],
                          [-s, 0.,  c, 0.],
                          [0., -s, 0.,  c]])
            Rt = R.T
            R6  = np.eye(6);  R6[:4, :4] = R
            Rt6 = np.eye(6);  Rt6[:4, :4] = Rt
            M = Rt6 @ M @ R6

        return M

    def constant_kick(self, delta: float = 0.0) -> np.ndarray:
        """Bend の ROTATE は設計ロール角 (transfer_matrix 内で処理) であってキックではない。
        基底クラスの d[3] += ROTATE (垂直キック) を継がない。K0 だけ水平キックとして返す (旧 ffs.py の線形 COD 用)。"""
        d = np.zeros(6)
        d[1] += self.K0
        return d

    def radiation_integrals(self) -> dict:
        """放射積分 I1..I5 への寄与 (横方向分散 η は外部から渡す必要あり)"""
        # η が未知のため完全計算は lattice.py 側で行う
        return {f'I{i}': 0.0 for i in range(1, 6)}


def _sbend_body(L: float, h: float, K1: float, delta: float = 0.0) -> np.ndarray:
    """
    セクターベンドのボディ転送行列 (6×6)。

    Parameters
    ----------
    L     : 要素長 [m]
    h     : 曲率 1/ρ [m⁻¹] (設計値)
    K1    : 四極強度 [m⁻²] (per-unit-length, SAD の ak1/al に相当; dp スケーリング前)
    delta : 運動量偏差 dp = (p-p0)/p0

    SAD 対応:
      akxsq = -(h^2 + K1) / (1+dp)   (tbend.f: akxsq=-1/rhosq-akk/p)
      akysq =         K1  / (1+dp)   (tbend.f: akysq=akk/p)
    """
    pr = 1.0 + delta
    M = np.eye(6)
    k2_0 = h**2 + K1      # スケーリング前 (分散の分母に使用)
    k2   = k2_0 / pr      # SAD: akxsq = -(h^2+K1)/p  → kx^2 = (h^2+K1)/p
    # --- x (偏向面) ---
    if abs(k2) < 1e-15:
        M[0, 0] = 1.0;  M[0, 1] = L
        M[1, 0] = 0.0;  M[1, 1] = 1.0
        M[0, 5] = L * h / 2   # 分散成分 (k2→0 の極限)
        M[1, 5] = h * L
    elif k2 > 0:
        kx = np.sqrt(k2)
        phix = kx * L
        sx, cx = np.sin(phix), np.cos(phix)
        M[0, 0] =  cx
        M[0, 1] =  sx / kx
        M[1, 0] = -kx * sx
        M[1, 1] =  cx
        # SAD tbendebody trans1(1,6) and trans1(2,6) at zero COD:
        # drhopp=h/p^2, drhop=delta*h/p, dcxkx=p*(1-cx)/k2_0, dcxkxp=[(1-cx)-0.5*sx*phix]/k2_0
        # sxkxp=0.5*(sx/kx - L*cx)/p
        # M[0,5] = drhopp*dcxkx + drhop*dcxkxp
        #        = h/k2_0 * [(1-cx) - 0.5*delta*sx*phix/pr]
        # M[1,5] = p*(drhopp*sxkx + drhop*sxkxp) + drhop*sxkx
        #        = h*sx/kx + 0.5*delta*h*(sx/kx - L*cx)/pr
        M[0, 5] = h/k2_0 * ((1-cx) - 0.5*delta*sx*phix/pr)
        M[1, 5] = h*sx/kx + 0.5*delta*h*(sx/kx - L*cx)/pr
    else:
        kx = np.sqrt(-k2)
        M[0, 0] =  np.cosh(kx * L)
        M[0, 1] =  np.sinh(kx * L) / kx
        M[1, 0] =  kx * np.sinh(kx * L)
        M[1, 1] =  np.cosh(kx * L)
        # For k2<0 (akxsq>0 in SAD): k2_0<0, kx=sqrt(-k2_0/p)
        # dcxkxp = [(chx-1)-0.5*shx*phix]/|k2_0|, drhop=delta*h/p
        # M[0,5] = h/|k2_0| * [(chx-1) - 0.5*delta*shx*phix/pr] = -h/k2_0 * [...]
        # M[1,5] = h*shx/kx + 0.5*delta*h*(shx/kx - L*chx)/pr
        M[0, 5] = -h/k2_0 * ((np.cosh(kx*L)-1) - 0.5*delta*np.sinh(kx*L)*kx*L/pr)
        M[1, 5] = h*np.sinh(kx*L)/kx + 0.5*delta*h*(np.sinh(kx*L)/kx - L*np.cosh(kx*L))/pr

    # --- y (非偏向面) ---
    # SAD: akysq = K1/p  → ky^2 = K1/pr
    K1_eff = K1 / pr
    if abs(K1_eff) < 1e-15:
        M[2, 2] = 1.0;  M[2, 3] = L
    elif K1_eff > 0:
        # K1 > 0: x-focusing → y-defocusing
        ky = np.sqrt(K1_eff)
        M[2, 2] =  np.cosh(ky * L)
        M[2, 3] =  np.sinh(ky * L) / ky
        M[3, 2] =  ky * np.sinh(ky * L)
        M[3, 3] =  np.cosh(ky * L)
    else:
        # K1 < 0: x-defocusing → y-focusing
        ky = np.sqrt(-K1_eff)
        M[2, 2] =  np.cos(ky * L)
        M[2, 3] =  np.sin(ky * L) / ky
        M[3, 2] = -ky * np.sin(ky * L)
        M[3, 3] =  np.cos(ky * L)

    # --- 縦方向 (path length): symplectic 1 次近似 ---
    M[4, 0] = -M[1, 5]
    M[4, 1] = -M[0, 5]
    M[4, 5] = h**2 * L / (abs(k2_0) + 1e-30) * pr - M[0, 5] * h
    return M


def _edge_matrix(h: float, e: float,
                 hgap: float = 0.0, fint: float = 0.5) -> np.ndarray:
    """偏向電磁石の端部集束行列 (6×6)"""
    M = np.eye(6)
    psi = hgap * fint * h * (1 + np.sin(e)**2) / np.cos(e)
    M[1, 0] =  h * np.tan(e)
    M[3, 2] = -h * np.tan(e - psi)
    return M


# ---------------------------------------------------------------------------
# 六極電磁石 (薄レンズ近似)
# ---------------------------------------------------------------------------

@dataclass
class Sextupole(Element):
    """
    SEXT 要素。薄レンズ近似では転送行列は非線形項のみ。
    通常の線形転送行列はドリフトと同じ。
    SAD コマンド: SEXT L=<長さ> K2=<六極強度>;
    """
    K2: float = 0.0        # 六極強度 [1/m³]

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """六極の 1 次転送行列はドリフト (qthin() 相当)。色収差補正あり。"""
        M = np.eye(6)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        return M

    def thin_kick(self, x: float, y: float, delta: float = 0.0) -> tuple:
        """
        薄レンズ近似での六極キック (非線形)。
        Δx' = -K2*L/2 * (x² - y²)
        Δy' = +K2*L * x*y
        """
        k2l = self.K2 * self.L / (1 + delta)
        dpx = -k2l / 2.0 * (x**2 - y**2)
        dpy =  k2l * x * y
        return dpx, dpy


# ---------------------------------------------------------------------------
# RF 空洞
# ---------------------------------------------------------------------------

@dataclass
class Cavity(Element):
    """
    CAVI 要素 (RF 空洞)。
    SAD コマンド: CAVI FREQ=<周波数[Hz]> VOLT=<電圧[eV]> PHI=<位相[rad]> E0=<ビームエネルギー[eV]>;

    パラメータ
    ----------
    FREQ : float  RF 周波数 [Hz]
    VOLT : float  ピーク電圧 [eV]  (= e * V_peak, V_peak はボルト単位)
    PHI  : float  同期位相 [rad]
                  PHI = π/2 → on-crest (最大加速), PHI = 0 → 最大縦方向集束
    E0   : float  基準ビームエネルギー [eV]  (空洞の相対電圧 VOLT/E0 に使用)

    物理
    ----
    縦方向正準キック (非線形):
        δ_out = δ + (VOLT/E0) * sin(PHI - k_rf * z)
    ここで k_rf = 2π * FREQ / c, z = -c * Δt (SAD 規約)

    線形転送行列 (同期粒子まわりの線形化):
        M[5,4] = -(VOLT/E0) * k_rf * cos(PHI)   [縦方向集束]
    """
    FREQ: float = 0.0      # RF 周波数 [Hz]
    VOLT: float = 0.0      # ピーク電圧 [eV]
    PHI: float  = 0.0       # SAD の既定値 (以前は π/2)# 同期位相 [rad]  (デフォルト: on-crest)
    E0: float   = 1.0e9    # 基準エネルギー [eV]  (デフォルト: 1 GeV)

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """
        空洞の 6×6 線形転送行列。

        横方向: 空洞長 L のドリフト (通常 L ≈ 0 の薄レンズ近似)
        縦方向: 同期位相 PHI まわりの線形化フォーカス
        """
        from .constants import SpeedOfLight
        M = np.eye(6)

        # 横方向ドリフト (空洞長がある場合)
        if self.L > 0:
            M[0, 1] = self.L
            M[2, 3] = self.L

        # 縦方向集束
        if self.VOLT != 0.0 and self.E0 != 0.0 and self.FREQ != 0.0:
            k_rf = 2.0 * np.pi * self.FREQ / SpeedOfLight
            v0   = self.VOLT / self.E0
            # ∂(δ_kick)/∂z = -v0 * k_rf * cos(PHI)
            M[5, 4] = -v0 * k_rf * np.cos(self.PHI)

        return M

    @property
    def k_rf(self) -> float:
        """RF 波数 [1/m]"""
        from .constants import SpeedOfLight
        return 2.0 * np.pi * self.FREQ / SpeedOfLight if self.FREQ > 0 else 0.0

    @property
    def v0(self) -> float:
        """規格化電圧 VOLT/E0 (無次元)"""
        return self.VOLT / self.E0 if self.E0 != 0 else 0.0

    def synchrotron_tune(self, alpha_c: float, circumference: float,
                         n_cav: int = 1) -> float:
        """
        シンクロトロン振動チューン Qs の近似推定。

        Qs ≈ sqrt(-α_c * h_rf * v0 * cos(PHI) / (2π))

        Parameters
        ----------
        alpha_c       : float  モメンタムコンパクション α_c
        circumference : float  リング周長 C [m]
        n_cav         : int    空洞個数 (合計電圧 = n_cav * VOLT)

        Returns
        -------
        Qs : float  シンクロトロンチューン (負の場合は不安定)
        """
        from .constants import SpeedOfLight
        h_rf = self.FREQ * circumference / SpeedOfLight  # ハーモニック数
        arg  = -alpha_c * h_rf * n_cav * self.v0 * np.cos(self.PHI) / (2.0 * np.pi)
        return np.sqrt(max(arg, 0.0))


# ---------------------------------------------------------------------------
# 進行波横偏向空洞 (TCAV)
# SAD コマンド: TCAV L=<長さ> AK=<偏向パラメータ> FREQ=<周波数[Hz]> PHI=<位相[rad]>;
# 参照: src/ttcav.f (非線形トラッキング), src/ttcave.f (線形転送行列)
# ---------------------------------------------------------------------------

@dataclass
class TCavity(Element):
    """
    TCAV 要素 (進行波横偏向空洞 / Transverse deflecting cavity)。

    SAD コマンド: TCAV L=<長さ> AK=<偏向パラメータ> FREQ=<周波数[Hz]>
                       PHI=<位相[rad]> HARM=<高調波数>;

    物理
    ----
    偏向空洞は横方向 (x) に粒子を蹴り、Panofsky-Wenzel の定理により
    縦方向 (z, δ) にも x に比例したエネルギー変化をもたらす。

    基本作用 (非線形; ttcav.f より):
      波数:         w = 2π·FREQ/c    (HARM=0 時), または w = OMEGA0·HARM/c
      位相:         phic = PHI·sign(CHARGE)
      通過時刻:     t = −z·h₁/p₁    (SAD 規約: z = −c·Δt)
      RF 位相:      ph = w·t + phic
      横方向キック: Δ(px/p₀) = −AK·sin(ph) / p₁r
      エネルギー変化: Δh = v·x·cos(ph),  v = AK·P0·w
      相対論的運動量更新: Newton 反復 (2回) で p₂r を求める

    転送行列 (線形化; ttcave.f, 基準軌道 x=z=0 での近似):
      主要結合要素 (0-indexed):
        R[1,4] = AK·w/β₀·cos(phic)    (px ← z 結合; SAD trans1(2,5))
        R[5,0] = AK·w/β₀·cos(phic)    (δ ← x 結合; SAD trans1(6,1))

    パラメータ
    ----------
    AK      : float  偏向パラメータ [rad] (SAD の ak; ≈ 最大偏向角)
                     AK = e·V_peak / (p₀·c)  (p₀ は基準運動量)
    FREQ    : float  RF 周波数 [Hz]
    HARM    : float  高調波数 (0 なら FREQ を使用)
    PHI     : float  RF 位相 [rad]  (同期位相)
    P0      : float  基準運動量 [自然単位: p₀/(m₀c), 無次元]
                     = βγ = Lorentz factor × β
                     例: 電子 6 GeV → P0 ≈ 11741
    H0      : float  基準全エネルギー [自然単位: E₀/(m₀c²), 無次元]
                     = γ = sqrt(1 + P0²).  指定しなければ自動計算。
    CHARGE  : float  粒子電荷符号 (+1 or −1; 電子は −1)
    OMEGA0  : float  周回角周波数 [rad/s]  (HARM ≠ 0 の場合に必要)

    Notes
    -----
    * P0, H0 は「自然単位」(無次元): p₀/(m₀c) および E₀/(m₀c²)。
      SI 単位からの変換:
        P0 = p₀[eV/c] / (m₀c²[eV])   (例: 電子 m₀c² = 0.511 MeV)
        H0 = E₀[eV] / (m₀c²[eV])

    * SAD のグローバル変数 ``dvfs`` は小さい補正項のため dvfs=0 と近似。
    """
    AK:     float = 0.0    # 偏向パラメータ [rad] (SAD の ak)
    FREQ:   float = 0.0    # RF 周波数 [Hz]
    HARM:   float = 0.0    # 高調波数 (0 なら FREQ を直接使用)
    PHI:    float = 0.0    # RF 位相 [rad]
    P0:     float = 1.0    # 基準運動量 [自然単位] (βγ)
    H0:     float = 0.0    # 基準エネルギー [自然単位] (γ). 0 なら自動計算
    CHARGE: float = -1.0   # 電荷符号 (電子: −1, 陽電子: +1, 陽子: +1)
    OMEGA0: float = 0.0    # 周回角周波数 [rad/s] (HARM≠0 時)

    def __post_init__(self):
        # H0 が未指定 (=0) なら P0 から自動計算: h0 = sqrt(1 + P0^2)
        if self.H0 == 0.0:
            self.H0 = float(np.sqrt(1.0 + self.P0**2))

    def _wave_number(self) -> float:
        """RF 波数 w [1/m] を返す。"""
        from .constants import SpeedOfLight, TwoPi
        if self.HARM == 0.0:
            return TwoPi * self.FREQ / SpeedOfLight
        else:
            return self.OMEGA0 * self.HARM / SpeedOfLight

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """
        TCAV の 6×6 線形転送行列 (SAD ttcave.f の近似実装)。

        基準軌道 x=0, z=0 (on-axis, 同期粒子) での線形化。
        空洞長 L がある場合は 半ドリフト–空洞キック–半ドリフト の構造。

        算出される非自明な行列要素 (0-indexed):
          R[1,4] = AK·w/β₀·cos(phic)   (px ← z)
          R[5,0] = AK·w/β₀·cos(phic)   (δ ← x)

        ここで β₀ = P0/H0 (基準粒子の速度/光速)。
        """
        w    = self._wave_number()
        phic = self.PHI * float(np.sign(self.CHARGE))
        cosp = np.cos(phic)
        # sinp = np.sin(phic)  # COD 更新に使うが transfer_matrix には不要

        # β₀ = v₁ = p₁/h₁ = P0/H0  (基準粒子速度)
        beta0 = self.P0 / self.H0

        # --- 空洞中央の薄レンズ転送行列 ---
        # (基準軌道 x=z=0, delta=0 での近似)
        Mcav = np.eye(6)
        if self.AK != 0.0 and w != 0.0:
            coupling = self.AK * w / beta0 * cosp
            Mcav[1, 4] = coupling   # Δ(px/p₀) ← z  [SAD trans1(2,5)]
            Mcav[5, 0] = coupling   # Δδ ← x         [SAD trans1(6,1)]
            # delta ≠ 0 の場合の高次補正 (SAD trans1(2,6), (6,6) 等) は
            # 基準軌道 t=0 の近似のため省略

        if self.L == 0.0:
            return Mcav

        # --- 半ドリフト–キック–半ドリフト ---
        Lhalf = self.L * 0.5 / (1.0 + delta)
        Md = np.eye(6)
        Md[0, 1] = Lhalf
        Md[2, 3] = Lhalf
        # z 方向ドリフト (線形近似: δ が大きくない範囲で)
        # SAD tdrift_free: z += dpz * L; 線形化 dpz ≈ -0.5*(px^2+py^2) - δ/(1+δ)
        # 一次近似: M[4,5] = -Lhalf (色収差ドリフト)
        Md[4, 5] = -Lhalf

        return Md @ Mcav @ Md


# ---------------------------------------------------------------------------
# モニタ / マーカー
# ---------------------------------------------------------------------------

@dataclass
class Monitor(Element):
    """MONI: ビーム位置モニタ (転送行列はドリフト)"""
    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        M = np.eye(6)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        return M


@dataclass
class Marker(Element):
    """
    MARK: マーカー要素 (単位行列)。

    ビームライン先頭の MARK の DX/DY/DPX/DPY は、アライメント誤差ではなく
    ins calc の初期軌道座標として使われる (ffs.py ins_calc)。
    なお SAD の CELL 計算では、求めた入口の閉軌道が先頭 MARK に書き戻され、
    DUMP では ORG.1 の DELX として表示される (SAD 実行で確認)。
    """
    L: float = 0.0
    DPX: float = 0.0   # 初期ビーム角 x' [rad]
    DPY: float = 0.0   # 初期ビーム角 y' [rad]

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        return np.eye(6)

    def constant_kick(self, delta: float = 0.0) -> np.ndarray:
        """MARK はキックを与えない (常にゼロ)。"""
        return np.zeros(6)


# ---------------------------------------------------------------------------
# 多極子 (薄レンズ)
# ---------------------------------------------------------------------------

@dataclass
class Multipole(Element):
    """
    MULT 要素 (集中多極キック)。
    係数は {KnL} の辞書で与える (n: 次数, 値: 強度×長さ)。
    """
    KnL: dict = field(default_factory=dict)   # {次数: 強度}

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """多極子の 1 次転送行列 (ドリフト + 四極成分キック)"""
        M = np.eye(6)
        # ドリフト成分 (L > 0 の場合, 色収差補正あり)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        # 四極成分の線形キック (薄レンズ近似)
        k1l = self.KnL.get(1, 0.0) / (1 + delta)
        if abs(k1l) > 0:
            M[1, 0] = -k1l
            M[3, 2] =  k1l
        return M


# ---------------------------------------------------------------------------
# 高次多極子 (OCTU / DECA / DODECA)
# SAD では qtwiss.f が qthin(nord=...) を呼び出し、
# 転送行列はドリフトと同じ (非線形項のみ)。
# tturn.f では tthin(nord=...) を呼び出す。
# ---------------------------------------------------------------------------

@dataclass
class Octupole(Element):
    """
    OCTU 要素 (八極電磁石)。
    SAD コマンド: OCTU L=<長さ> K3=<強度>;
    転送行列はドリフトと同じ。非線形トラッキングは tthin(nord=8) に相当。
    """
    K3: float = 0.0   # 八極強度 [m⁻⁴]

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        M = np.eye(6)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        return M


@dataclass
class Decapole(Element):
    """
    DECA 要素 (十極電磁石)。
    SAD コマンド: DECA L=<長さ> K4=<強度>;
    転送行列はドリフトと同じ。非線形トラッキングは tthin(nord=10) に相当。
    """
    K4: float = 0.0   # 十極強度 [m⁻⁵]

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        M = np.eye(6)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        return M


@dataclass
class Dodecapole(Element):
    """
    DODECA 要素 (十二極電磁石)。
    SAD コマンド: DODECA L=<長さ> K5=<強度>;
    転送行列はドリフトと同じ。非線形トラッキングは tthin(nord=12) に相当。
    """
    K5: float = 0.0   # 十二極強度 [m⁻⁶]

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        M = np.eye(6)
        Leff = self.L / (1.0 + delta)
        M[0, 1] = Leff
        M[2, 3] = Leff
        return M


# ---------------------------------------------------------------------------
# ソレノイド (SOL)
# SAD qtwiss.f: tsolrote() を呼び出し (tbende.f/tsol.f に実装)。
# SAD tturn.f:  tsol() を呼び出し (tsol.f に実装)。
# 線形行列: 標準ソレノイド行列 (Wiedemann, Sec.8)
# ---------------------------------------------------------------------------

@dataclass
class Solenoid(Element):
    """
    SOL 要素 (ソレノイド電磁石)。
    SAD コマンド: SOL L=<長さ> BZ=<軸方向磁場[T]> P0=<基準運動量[eV/c]>;

    BzL = BZ * L  [T·m]  (積分磁場)
    Ks  = BzL / (2 * Brho)  [1/m]  (ソレノイド強度パラメータ)
    Brho = P0 / (e*c) [T·m]

    線形行列 (SAD tsolrote() 相当、回転座標系):
    SADのqtwiss.fはソレノイドの線形行列を「implementation error」で
    スキップしているため、ここでは標準的なソレノイド行列を実装する。
    """
    BZ: float  = 0.0    # 軸方向磁場 [T]
    P0: float  = 3.0e9  # 基準運動量 [eV/c]  (デフォルト: 3 GeV/c)
    KS: float  = 0.0    # ソレノイド強度 Ks = BZ*L/(2*Brho) [1/m]
                        # BZ と P0 が与えられた場合は自動計算

    def __post_init__(self):
        if self.BZ != 0.0 and self.P0 > 0.0 and self.KS == 0.0:
            # Brho = P0 [eV/c] / (c [m/s]) / (e [C]) * c = P0 / (e*c) [T·m]
            # P0 は eV/c 単位: Brho = P0[eV/c] / (c/c) / e = P0 / (3e8) T·m
            # 簡略化: Brho ≈ P0[GeV/c] / 0.3 T·m
            from .constants import SpeedOfLight
            Brho = self.P0 / SpeedOfLight      # [T·m]  (P0 in eV/c)
            self.KS = self.BZ * self.L / (2.0 * Brho)

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        """
        ソレノイドの 6×6 線形転送行列。
        参照: Wolski, Beam Dynamics in High-Energy Particle Accelerators, 2014.
        """
        Ks = self.KS / (1.0 + delta)
        L  = self.L
        M  = np.eye(6)
        if abs(Ks) < 1e-15 or abs(L) < 1e-15:
            # 弱場 → ドリフト近似
            Leff = L / (1.0 + delta)
            M[0, 1] = Leff
            M[2, 3] = Leff
            return M
        phi = Ks * L    # ソレノイド位相
        C = np.cos(phi)
        S = np.sin(phi)
        # ソレノイド転送行列 (Courant-Snyder 座標, 2x4 結合)
        M[0, 0] =  C**2
        M[0, 1] =  C * S / Ks
        M[0, 2] =  C * S
        M[0, 3] =  S**2 / Ks
        M[1, 0] = -Ks * C * S
        M[1, 1] =  C**2
        M[1, 2] = -Ks * S**2
        M[1, 3] =  C * S
        M[2, 0] = -C * S
        M[2, 1] = -S**2 / Ks
        M[2, 2] =  C**2
        M[2, 3] = -C * S / Ks
        M[3, 0] =  Ks * S**2
        M[3, 1] = -C * S
        M[3, 2] = -Ks * C * S
        M[3, 3] =  C**2
        # 縦方向: ドリフト近似 (ソレノイドの縦結合は微小)
        M[4, 5] = 0.0  # no longitudinal focusing in thin-lens approx
        return M


# ---------------------------------------------------------------------------
# 要素ファクトリ (SAD コマンド風)
# ---------------------------------------------------------------------------

def _wrap_rotation(cls):
    """cls.transfer_matrix の結果に ROTATION ロール誤差を適用するよう包む。"""
    orig = cls.transfer_matrix

    def transfer_matrix(self, delta: float = 0.0) -> np.ndarray:
        return self._apply_rotation(orig(self, delta))
    transfer_matrix.__doc__ = orig.__doc__
    cls.transfer_matrix = transfer_matrix
    return cls


for _cls in (Drift, Quadrupole, Bend, Sextupole, Monitor, Multipole,
             Octupole, Decapole, Dodecapole, Solenoid):
    _wrap_rotation(_cls)


def _bend_constant_kick(self, delta: float = 0.0) -> np.ndarray:
    """
    偏向電磁石のキックベクトル。ROTATE は設計ロール角 (transfer_matrix 内で処理)
    なのでキックに含めない。K0 のみ水平薄レンズキックとして返す。
    """
    return np.array([0.0, float(self.K0), 0.0, 0.0, 0.0, 0.0])


Bend.constant_kick = _bend_constant_kick


ELEMENT_TYPES = {
    'DRIFT':   Drift,
    'QUAD':    Quadrupole,
    'BEND':    Bend,
    'SEXT':    Sextupole,
    'OCTU':    Octupole,
    'DECA':    Decapole,
    'DODECA':  Dodecapole,
    'SOL':     Solenoid,
    'CAVI':    Cavity,
    'TCAV':    TCavity,
    'MONI':    Monitor,
    'MARK':    Marker,
    'MULT':    Multipole,
}


def make_element(etype: str, name: str = '', **kwargs) -> Element:
    """
    要素タイプ名と引数から要素を生成する。

    Parameters
    ----------
    etype : str  要素タイプ ('DRIFT', 'QUAD', ...)
    name : str   要素名
    **kwargs     要素パラメータ

    Returns
    -------
    elem : Element のサブクラスのインスタンス
    """
    etype_upper = etype.upper()
    cls = ELEMENT_TYPES.get(etype_upper)
    if cls is None:
        raise ValueError(f"未知の要素タイプ: {etype}")
    return cls(name=name, **kwargs)


__all__ = [
    'Element', 'Drift', 'Quadrupole', 'Bend', 'Sextupole',
    'Octupole', 'Decapole', 'Dodecapole', 'Solenoid',
    'Cavity', 'TCavity', 'Monitor', 'Marker', 'Multipole',
    'make_element', 'ELEMENT_TYPES',
]
