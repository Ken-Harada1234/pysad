"""
pysad.ffs
=========
SAD FFS (Final Focus System) コマンドレベルの Python 実装。

SAD の FFS コマンドファイルで使われる主要コマンドを
Python オブジェクトとして提供する。ループ等は Python 側で実装すること。

対応コマンド:
  USE name           アクティブなビームラインを設定
  cell calc          周期 Twiss (閉軌道解) の計算
  ins calc           指定初期条件からの Twiss 伝播
  emit               放射積分・エミッタンスの計算
  free pattern param フィッティングのフリーパラメータ設定
  fit element param target nsteps   フィット制約の追加
  go [n]             フィット実行

クエリ関数 (Python API):
  ffs.Twiss("BX", "QF")     → Twiss パラメータを返す
  ffs.LINE("K1", "QF")      → 要素パラメータを返す
  ffs.LINE("K1", "QF", val) → 要素パラメータをセット

使い方:
    from pysad import load_sad_file, FFSState

    interp = load_sad_file("pf_ring.n")
    ffs = FFSState(interp.env)
    ffs.use("RING")
    ffs.cell_calc()

    print(ffs.Twiss("BX", "QF"))
    print(ffs.LINE("K1", "Q*"))

    # fit example
    ffs.free("Q*", "K1")
    ffs.fit_add(element="$$$", param="NX", target=9.5975, nsteps=3)
    ffs.fit_add(element="$$$", param="NY", target=5.2790, nsteps=3)
    ffs.go(niter=6)
    ffs.cell_calc()
"""

from __future__ import annotations

import re
import fnmatch
import math
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union


# ---------------------------------------------------------------------------
# SADList: SAD スタイルの 1-based インデックスリスト
# ---------------------------------------------------------------------------

class SADList:
    """
    SAD スタイルの 1-based インデックスリスト。

    SAD では配列インデックスが 1 から始まる。
    Python での利用時は通常の 0-based も使える。

    Examples
    --------
    >>> sl = SADList([1.1, 2.2, 3.3])
    >>> sl[1]   # SAD 1-based → 1.1
    1.1
    >>> list(sl)  # Python iteration (0-based) → [1.1, 2.2, 3.3]
    [1.1, 2.2, 3.3]
    """

    def __init__(self, items):
        self._items = list(items)

    def __getitem__(self, i: int):
        """1-based インデックスアクセス (SAD 互換)。"""
        if isinstance(i, int):
            if 1 <= i <= len(self._items):
                return self._items[i - 1]
            raise IndexError(f"SADList index {i} out of range [1, {len(self._items)}]")
        # slice はそのまま委譲 (Python 流)
        return self._items[i]

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __repr__(self) -> str:
        return f"SADList({self._items})"

    def __add__(self, other):
        if isinstance(other, SADList):
            return SADList(self._items + other._items)
        return NotImplemented

    def to_list(self) -> list:
        return list(self._items)

    def to_array(self) -> np.ndarray:
        return np.array(self._items)


# ---------------------------------------------------------------------------
# パターンマッチング
# ---------------------------------------------------------------------------

_SPECIAL_ELEMENTS = ('$$$', '***')


def _match_pattern(name: str, pattern: str) -> bool:
    """
    SAD ワイルドカードパターンでマッチング。

    SAD では * がワイルドカード。 .* も * と同じ意味。
    {FIO} のようなキャラクタクラスも対応 (例: PM{FIO}* → PMF*, PMI*, PMO* の和)。
    大文字小文字を区別しない。
    """
    import re as _re
    name_up = name.upper()
    pat = pattern.upper().replace('.*', '*')

    # {ABC} 形式を正規表現 [ABC] に展開してからマッチ
    if '{' in pat:
        # fnmatch の * → .* に変換して正規表現に
        rx = ''
        i = 0
        while i < len(pat):
            ch = pat[i]
            if ch == '{':
                j = pat.index('}', i)
                rx += '[' + pat[i+1:j] + ']'
                i = j + 1
            elif ch == '*':
                rx += '.*'
            elif ch == '?':
                rx += '.'
            else:
                rx += _re.escape(ch)
            if i >= len(pat):
                break
            i = i if ch == '{' else i + 1  # already advanced past {..}
        # fix: above loop has off-by-one for normal chars
        # rebuild cleanly
        rx = ''
        i = 0
        while i < len(pat):
            ch = pat[i]
            if ch == '{':
                j = pat.index('}', i + 1)
                rx += '[' + pat[i+1:j] + ']'
                i = j + 1
            elif ch == '*':
                rx += '.*'
                i += 1
            elif ch == '?':
                rx += '.'
                i += 1
            else:
                rx += _re.escape(ch)
                i += 1
        return bool(_re.fullmatch(rx, name_up))

    return fnmatch.fnmatch(name_up, pat)


# ---------------------------------------------------------------------------
# メインの FFS 状態クラス
# ---------------------------------------------------------------------------

class FFSState:
    """
    SAD FFS コマンドレベルの Python 実装。

    SADEnvironment を受け取り、その上で FFS 計算を行う。
    parser.py の SADScript.ffs プロパティからアクセスできる。

    Parameters
    ----------
    env : SADEnvironment
        pysad パーサーの環境オブジェクト
    momentum_eV : float, optional
        ビームの運動量 [eV]。省略時は env の MOMENTUM 変数から取得。
    """

    def __init__(self, env, momentum_eV: Optional[float] = None):
        self._env = env
        self._momentum_eV = momentum_eV  # 上書き可

        # 現在のビームライン
        self._beamline_name: Optional[str] = None
        self._elements: list = []   # Element オブジェクト
        self._names: list = []      # 要素名 (parallel to _elements)
        self._s: list = []          # s 座標 (各要素出口) [m]

        # 計算済み Twiss
        # _twiss[0] = 入口, _twiss[i] = i 番目要素の出口
        self._twiss: list = []
        self._twiss0: Optional[dict] = None   # 周期解

        # エミッタンス
        self._emitx: Optional[float] = None
        self._emity: Optional[float] = None
        self._sige: Optional[float] = None   # エネルギー広がり

        # フィッティング
        self._free_params: list = []       # [(idx, name, param)]
        self._constraints: list = []       # list of dict

        # 多点マッチング (chromaticity-aware fitting)
        self._dp: float = 0.0             # DP: 運動量範囲 ±DP (FITP 用)
        self._fit_npoints: int = 1        # FITP N: 運動量サンプル点数

        # SAVE / SAVE ALL — マッチング後のデフォルト値スナップショット
        self._saved_params: dict = {}     # {element_name: {attr: value}}

    # ------------------------------------------------------------------
    # USE / ビームライン設定
    # ------------------------------------------------------------------

    def use(self, beamline_name: str) -> 'FFSState':
        """
        USE beamline_name — アクティブなビームラインを設定する。

        Parameters
        ----------
        beamline_name : str
            ビームライン名 (大文字・小文字を区別しない)

        Returns
        -------
        self : FFSState (チェーン呼び出し可)
        """
        name = beamline_name.upper()
        lat = self._env.build_lattice(name)
        self._beamline_name = name
        self._elements = lat.elements
        self._names = [e.name for e in self._elements]
        # s 座標を計算
        s = 0.0
        self._s = []
        for e in self._elements:
            s += e.L
            self._s.append(s)
        # 計算済みキャッシュをクリア
        self._twiss = []
        self._twiss0 = None
        return self

    # ------------------------------------------------------------------
    # cell calc / ins calc / emit
    # ------------------------------------------------------------------

    def cell_calc(self) -> 'FFSState':
        """
        cell calc — 周期 Twiss (閉軌道解) を計算する。

        SAD の "FFS[\"cell calc\"]" に相当。
        """
        from .optics import periodic_optics
        from .lattice import _twiss_at_elements

        if not self._elements:
            raise RuntimeError("cell_calc: USE でビームラインを設定してください")

        # 6×6 一周転送行列
        M6 = np.eye(6)
        for el in self._elements:
            M6 = el.transfer_matrix() @ M6

        M4 = M6[:4, :4]
        try:
            twiss0 = periodic_optics(M4)
        except Exception as e:
            raise RuntimeError(f"cell_calc: 周期 Twiss の計算に失敗 — {e}")

        # 分散 (5×5 行列から閉軌道分散を求める)
        try:
            # η = (I - M4)^{-1} · d  where d = M6[:4, 5]
            I_M4 = np.eye(4) - M4
            d = M6[:4, 5]
            eta = np.linalg.solve(I_M4, d)
            twiss0['ex'] = float(eta[0])
            twiss0['epx'] = float(eta[1])
            twiss0['ey'] = float(eta[2])
            twiss0['epy'] = float(eta[3])
        except np.linalg.LinAlgError:
            twiss0.setdefault('ex', 0.0)
            twiss0.setdefault('epx', 0.0)
            twiss0.setdefault('ey', 0.0)
            twiss0.setdefault('epy', 0.0)

        twiss0['nx'] = 0.0   # 累積位相は 0 スタート (setdefault では上書きされない)
        twiss0['ny'] = 0.0
        twiss0.setdefault('dx', 0.0)
        twiss0.setdefault('dy', 0.0)

        self._twiss0 = twiss0
        self._twiss = _twiss_at_elements(M6, self._elements, twiss0)

        # dx/dy (COD) を各 Twiss 辞書に初期化してから COD 計算
        for tw in self._twiss:
            tw.setdefault('dx', 0.0)
            tw.setdefault('dy', 0.0)
        self.calc_cod()
        return self

    def ins_calc(self, twiss_init: Optional[dict] = None) -> 'FFSState':
        """
        ins calc — 指定した初期条件から Twiss を伝播する。

        SAD の "FFS[\"ins calc\"]" に相当。
        twiss_init を省略した場合は cell_calc() の結果を使う。

        SAD の規則:
          ビームライン先頭の MARK 要素が存在する場合、その DX/DPX/DY/DPY を
          初期軌道オフセット (x, x', y, y') として使用する。
          これは設置誤差ではなく、ビームの出発点座標を意味する。
          他の要素の DX/DY は通常のアライメント誤差として扱われる。
        """
        from .lattice import _twiss_at_elements
        from .elements import Marker

        tw0 = twiss_init
        if tw0 is None:
            if self._twiss0 is None:
                self.cell_calc()
            tw0 = self._twiss0

        M6 = np.eye(6)
        self._twiss = _twiss_at_elements(M6, self._elements, tw0)

        # 先頭 MARK の DX/DPX/DY/DPY を初期軌道オフセットとして設定
        if self._elements and isinstance(self._elements[0], Marker):
            mk = self._elements[0]
            x0 = np.array([
                float(mk.DX),  float(mk.DPX),
                float(mk.DY),  float(mk.DPY),
            ])
            # _twiss[0] = ビームライン入口 (MARK 直後)
            self._twiss[0]['dx'] = float(x0[0])
            self._twiss[0]['dy'] = float(x0[2])
            # 各要素の転送行列 + キックで軌道を逐次伝播
            x = x0.copy()
            for i, el in enumerate(self._elements):
                M4 = el.transfer_matrix()[:4, :4]
                d4 = el.constant_kick()[:4]
                x = M4 @ x + d4
                if i + 1 < len(self._twiss):
                    self._twiss[i + 1]['dx'] = float(x[0])
                    self._twiss[i + 1]['dy'] = float(x[2])
        return self

    def emit_calc(self, momentum_eV: Optional[float] = None) -> dict:
        """
        emit — 放射積分・エミッタンスを計算する。

        SAD の "FFS[\"cell calc emit\"]" の emit 部分に相当。
        """
        from .lattice import Lattice

        lat = Lattice(self._elements)

        # 運動量の決定
        p = momentum_eV or self._momentum_eV
        if p is None:
            raw = self._env.get('MOMENTUM')
            if raw is not None:
                p = float(raw)
        if p is None:
            raise RuntimeError(
                "emit_calc: 運動量が不明です。"
                "FFSState(..., momentum_eV=...) か"
                " MOMENTUM 変数を設定してください"
            )

        ri = lat.radiation_integrals()
        em = lat.emittance(float(p))
        self._emitx = em.get('emitx')
        self._emity = em.get('emity')
        self._sige  = em.get('sige')
        return {**ri, **em}

    def calc(self, mode: str = 'cell') -> 'FFSState':
        """
        calc — FFS の calc コマンド。

        Parameters
        ----------
        mode : 'cell' or 'ins'
        """
        if mode == 'cell':
            self.cell_calc()
        else:
            self.ins_calc()
        return self

    def calc_cod(self) -> 'FFSState':
        """
        K0 / ROTATE / DX / DY から閉軌道偏差 (COD) を計算して
        各要素位置の _twiss[i]['dx'], _twiss[i]['dy'] に格納する。

        SAD の定義:
          K0     — 偏向電磁石等の追加水平キック [rad]  (中心軌道と別)
          ROTATE — 垂直ステアリング用キック [rad]
          DX     — 水平設置誤差 [m] → 四極では等価キック K1*L*DX が発生
          DY     — 垂直設置誤差 [m]

        計算アルゴリズム (1 周期解法):
          各要素 i での薄レンズキック d_i を収集し、
            (I - M_ring) · x_0 = Σ_i M_{i→end} · d_i
          を解いて入口 COD x_0 を求め、逐次伝播する。

        Returns
        -------
        self : FFSState
        """
        n = len(self._elements)
        if n == 0 or not self._twiss:
            return self

        # 各要素の 4×4 転送行列と薄レンズキックベクトルを収集
        Ms = []    # 4×4 転送行列リスト
        ds = []    # キックベクトル [Δx, Δx', Δy, Δy'] リスト

        for el in self._elements:
            M6 = el.transfer_matrix()    # ROTATION ロール誤差を含む 6×6 行列
            Ms.append(M6[:4, :4])

            # constant_kick() は各要素タイプの慣例を正確に反映する:
            #   - Element (基底): [0, K0, 0, ROTATE, 0, 0]
            #   - Quadrupole: K0/ROTATE + DX/DY フィードダウンキック
            #   - Bend: [0, K0, 0, 0, 0, 0]  (ROTATE は設計ロール角でキックではない)
            d6 = el.constant_kick()
            ds.append(d6[:4])

        # 1 周転送行列 (4×4)
        M_ring = np.eye(4)
        for M in Ms:
            M_ring = M @ M_ring

        # 合計蹴り D = Σ M_{i+1→n} · d_i  (要素 i のキックがリング終端に伝播)
        D = np.zeros(4)
        for i in range(n):
            if np.any(np.abs(ds[i]) > 1e-15):
                # M_{i+1 → n}: 要素 i のキック後からリング終端までの転送行列
                M_tail = np.eye(4)
                for j in range(i + 1, n):
                    M_tail = Ms[j] @ M_tail
                D += M_tail @ ds[i]

        # (I - M_ring) · x_0 = D を解いて入口 COD を求める
        try:
            A = np.eye(4) - M_ring
            if abs(np.linalg.det(A)) < 1e-12:
                # 縮退 (整数チューン等): COD は無限大に発散 → 0 で代替
                x_cod = np.zeros(4)
            else:
                x_cod = np.linalg.solve(A, D)
        except np.linalg.LinAlgError:
            x_cod = np.zeros(4)

        # _twiss[0] (入口) に COD を格納
        if self._twiss:
            self._twiss[0]['dx'] = float(x_cod[0])
            self._twiss[0]['dy'] = float(x_cod[2])

        # 逐次伝播: 各要素出口の COD を計算して _twiss[i+1] に格納
        x_cur = x_cod.copy()
        for i, (M, d) in enumerate(zip(Ms, ds)):
            x_after = M @ x_cur + d
            if i + 1 < len(self._twiss):
                self._twiss[i + 1]['dx'] = float(x_after[0])
                self._twiss[i + 1]['dy'] = float(x_after[2])
            x_cur = x_after

        return self

    # ------------------------------------------------------------------
    # COD 応答行列 (Response Matrix)
    # ------------------------------------------------------------------

    def cod_response(
        self,
        corr_names,
        bpm_names,
        K0_kick: float = 0.5e-3,
        plane: str = 'both',
    ) -> dict:
        """
        COD 応答行列を数値計算する。

        SAD の write_response_COD スクリプト相当:
          LINE["K0", corr] = K0_kick → FFS["cell calc"] → Twiss["DX"/"DY", BPM]

        Parameters
        ----------
        corr_names : list[str]
            水平補正磁石名リスト (B* パターン)、または
            垂直補正磁石名リスト (V* パターン)。
            ``plane`` 引数で自動振り分けされる。
        bpm_names : list[str]
            BPM 名リスト (SAD の PM{FIO}* パターン等)。
        K0_kick : float
            K0 摂動量 [rad]。SAD デフォルト = 0.5/1000 = 5e-4 rad。
        plane : str
            'x'   — 水平 COD 応答 (corr は B* 等の水平補正磁石)
            'y'   — 垂直 COD 応答 (corr は V* 等の垂直補正磁石)
            'both'— 両者を計算 (corr_names が水平・垂直両方を含む場合は 'x'/'y' を明示)

        Returns
        -------
        dict with keys:
            'Rx'  : np.ndarray, shape (n_bpm, n_corrx), 水平応答 [mm / (K0 kick)]
            'Ry'  : np.ndarray, shape (n_bpm, n_corry), 垂直応答 [mm / (K0 kick)]
            'Rx_dict': dict[(j,i)] = value  (j=1-based corr, i=1-based BPM)
            'Ry_dict': dict[(j,i)] = value

        Notes
        -----
        * 水平補正磁石 (B*): SAD の K0 は素子中心での薄レンズキック。
          正の K0 → より内向きに偏向 → x_co < 0。実効キック = -K0_kick。
          素子が厚い (L > 0) 場合は BEND の前半行列 (E1 含む) で中心位置を求める。
        * 垂直補正磁石 (V*): ROTATE=-90 DEG で定義。局所 x' キック → 実験室 y' へ -K0。
          pysad のパーサが V* を Drift(L=0) と解析するため、
          K0 の方向は ROTATE 変換 (R(-90)) に基づき -K0_kick とする。
        """
        els = self._elements
        nms = self._names
        N = len(els)

        # --- 累積転送行列 T[i]: 素子 i の入口まで ---
        Mel = [e.transfer_matrix() for e in els]
        T = [np.eye(6)]
        for M in Mel:
            T.append(M @ T[-1])
        M_ot = T[N]
        M4 = M_ot[:4, :4]
        invIminusM4 = np.linalg.inv(np.eye(4) - M4)

        # --- 名前 → 最初の出現インデックス ---
        name_to_idx: Dict[str, int] = {}
        for i, nm in enumerate(nms):
            if nm not in name_to_idx:
                name_to_idx[nm] = i

        bpm_indices = [name_to_idx[n] for n in bpm_names if n in name_to_idx]

        from pysad.elements import Bend as _Bend

        def _half_bend_mat(el):
            """素子中心までの転送行列 (E1 エッジ込み)。"""
            if el.L < 1e-12:
                return np.eye(6)
            # 前半: BEND(L/2, ANGLE/2) に E1 を適用, E2=0
            return _Bend(
                name='_half',
                L=el.L / 2,
                ANGLE=getattr(el, 'ANGLE', 0.) / 2,
                K1=getattr(el, 'K1', 0.),
                E1=getattr(el, 'E1', 0.),
                E2=0.,
                HGAP=getattr(el, 'HGAP', 0.),
                FINT=getattr(el, 'FINT', 0.5),
            ).transfer_matrix()

        def _orbit_at_bpms(corr_idx: int, kick_x: float, kick_y: float):
            """キック摂動に対する BPM での 4 ベクトル [x, x', y, y'] を返す。"""
            el = els[corr_idx]
            if el.L < 1e-12:
                T_kick = T[corr_idx]
            else:
                T_kick = _half_bend_mat(el) @ T[corr_idx]

            M_after = T[N] @ np.linalg.inv(T_kick)
            d4 = np.array([0., kick_x, 0., kick_y])
            x0 = invIminusM4 @ (M_after[:4, :4] @ d4)

            out = []
            for bi in bpm_indices:
                xa = T[bi][:4, :4] @ x0
                if bi > corr_idx:
                    xa = xa + (T[bi] @ np.linalg.inv(T_kick))[:4, :4] @ d4
                out.append(xa)
            return out

        # --- 水平応答 ---
        Rx_dict: Dict[tuple, float] = {}
        if plane in ('x', 'both') and corr_names:
            for j, cn in enumerate(corr_names):
                ci = name_to_idx.get(cn)
                if ci is None:
                    continue
                orbs = _orbit_at_bpms(ci, -K0_kick, 0.)
                for i, orb in enumerate(orbs):
                    Rx_dict[(j + 1, i + 1)] = orb[0] * 1000   # [mm]

        # --- 垂直応答 ---
        Ry_dict: Dict[tuple, float] = {}
        if plane in ('y', 'both') and corr_names:
            for j, cn in enumerate(corr_names):
                ci = name_to_idx.get(cn)
                if ci is None:
                    continue
                # V* (ROTATE=-90): 局所 x' キック → 実験室 y' = -K0
                orbs = _orbit_at_bpms(ci, 0., -K0_kick)
                for i, orb in enumerate(orbs):
                    Ry_dict[(j + 1, i + 1)] = orb[2] * 1000   # [mm]

        nb = len(bpm_indices)
        nc = len(corr_names)

        if Rx_dict:
            Rx = np.array([[Rx_dict.get((j + 1, i + 1), 0.)
                            for j in range(nc)] for i in range(nb)])
        else:
            Rx = np.empty((nb, 0))

        if Ry_dict:
            Ry = np.array([[Ry_dict.get((j + 1, i + 1), 0.)
                            for j in range(nc)] for i in range(nb)])
        else:
            Ry = np.empty((nb, 0))

        return {'Rx': Rx, 'Ry': Ry, 'Rx_dict': Rx_dict, 'Ry_dict': Ry_dict}

    # ------------------------------------------------------------------
    # Twiss クエリ
    # ------------------------------------------------------------------

    # SAD Twiss キー → pysad 内部キー
    _TWISS_KEY_MAP = {
        'BX': 'bx', 'BY': 'by', 'AX': 'ax', 'AY': 'ay',
        'NX': 'nx', 'NY': 'ny',
        'EX': 'ex', 'EPX': 'epx',
        'EY': 'ey', 'EPY': 'epy',
        'DX': 'dx', 'DY': 'dy',   # COD (Closed Orbit Distortion)
        'GX': None, 'GY': None,    # ガンマ (計算)
    }

    def Twiss(self, param: str, element_arg) -> Any:
        """
        Twiss["BX", element] — 要素位置での Twiss パラメータを取得する。

        Parameters
        ----------
        param : str
            Twiss パラメータ名 (BX, BY, AX, AY, NX, NY, EX, EPX, ...)
        element_arg : str, int, or float
            要素名 / ワイルドカードパターン / 要素番号 (1-based, 端数可)
            端数の場合: Twiss["BX", 10.5] → 10 番目要素の中央点で線形補間

        Returns
        -------
        float or SADList
            1 要素の場合はスカラー、複数の場合は SADList。

        Examples
        --------
        >>> ffs.Twiss("BX", "QF")
        11.32
        >>> ffs.Twiss("NX", "$$$")    # リング末端 (= 2π × 全チューン)
        60.28...
        >>> ffs.Twiss("BX", "Q*")     # 全 Q 要素
        SADList([...])
        >>> ffs.Twiss("BX", 10.5)     # 10 番目要素の中央 (線形補間)
        9.87
        """
        if not self._twiss:
            self.cell_calc()

        key = self._TWISS_KEY_MAP.get(param.upper())

        # --- 端数インデックス (線形補間) ---
        if isinstance(element_arg, float) and element_arg != int(element_arg):
            return self._twiss_interpolate(key, param, element_arg)

        if key is None:
            # ガンマ: γ = (1+α²)/β (入口 Twiss を使用)
            if param.upper() in ('GX',):
                vals = []
                for idx in self._resolve_indices(element_arg):
                    tw = self._twiss[idx]
                    ax, bx = tw.get('ax', 0), tw.get('bx', 1)
                    vals.append((1 + ax**2) / bx)
                return vals[0] if len(vals) == 1 else SADList(vals)
            return None

        # *** / $$$ は特別: リング末端の全チューン (_twiss[N]) を返す
        if isinstance(element_arg, str) and element_arg.upper() in _SPECIAL_ELEMENTS:
            N = len(self._elements)
            val = self._twiss[N].get(key, 0.0)
            return val

        indices = self._resolve_indices(element_arg)
        if not indices:
            return None

        # SAD 互換の「入口 Twiss」慣習:
        # SAD の Twiss["BX","B02"] は B02 の入口 (= 直前要素の出口) の Twiss を返す。
        # 零長要素 (MONI, MARK, 零角 BEND) では入口 = 出口なので従来と同じ結果。
        # 有長要素 (B*, Q* など) では入口 Twiss が必要。
        vals = [self._twiss[i].get(key, 0.0) for i in indices]
        return vals[0] if len(vals) == 1 else SADList(vals)

    def _twiss_interpolate(self, key: Optional[str], param: str, frac_idx: float) -> float:
        """
        端数インデックスでの Twiss 値を線形補間する。

        SAD では Twiss["BX", 10.5] は 10 番目要素の入口と出口の
        中点 (s 方向で 50%) での値を意味する。

        _twiss のインデックス対応:
          _twiss[0]   = ビームライン入口 (要素 1 の入口)
          _twiss[i]   = i 番目要素の出口 (SAD 1-based: 要素 i の出口)
          _twiss[n]   = 最終要素の出口 (= $$$ の位置)

        frac_idx = 10.5 → 10 番目要素の中点
          入口 = _twiss[9]  (10 番目要素の入口 = 9 番目要素の出口)
          出口 = _twiss[10] (10 番目要素の出口)
          frac = 0.5
        """
        if key is None:
            return 0.0

        n = len(self._elements)
        floor_i = int(math.floor(frac_idx))    # SAD 1-based 整数部
        frac    = frac_idx - floor_i

        # _twiss インデックス: 要素 floor_i の入口 = _twiss[floor_i - 1]
        idx_lo = floor_i - 1    # 入口側 (_twiss インデックス)
        idx_hi = floor_i        # 出口側 (_twiss インデックス)

        idx_lo = max(0, min(idx_lo, n))
        idx_hi = max(0, min(idx_hi, n))

        tw_lo = self._twiss[idx_lo].get(key, 0.0)
        tw_hi = self._twiss[idx_hi].get(key, 0.0)
        return float(tw_lo + frac * (tw_hi - tw_lo))

    # ------------------------------------------------------------------
    # LINE クエリ / セット
    # ------------------------------------------------------------------

    # 要素属性マップ (SAD パラメータ名 → Python 属性名)
    _EL_ATTR_MAP = {
        'K1': 'K1', 'K2': 'K2', 'K3': 'K3',
        'L': 'L', 'ANGLE': 'ANGLE',
        'E1': 'E1', 'E2': 'E2',
        'K0': 'K0',              # 水平薄レンズキック [rad] (COD 用)
        'ROTATE': 'ROTATE',      # 垂直薄レンズキック [rad] (垂直 COD 用)
        'DX': 'DX',              # 水平設置誤差 [m]
        'DY': 'DY',              # 垂直設置誤差 [m]
        'ROTATION': 'ROTATION',  # 据え付けロール誤差 [rad] (全磁石共通)
        'BY': 'K1',              # 旧 SAD 表記 (四極: K1 = BY)
        'BX': 'K0',              # 水平キック (旧表記)
    }

    def LINE(self, param: str,
             element_arg=None,
             value=None) -> Any:
        """
        LINE["param", element] / LINE["param", element] = value

        Parameters
        ----------
        param : str
            要素パラメータ / クエリ種別 (K1, K2, L, NAME, POSITION, ...)
        element_arg : str or int or None
            要素名・パターン・番号。省略時は全要素。
        value : float or None
            None → get モード, 値あり → set モード

        Returns
        -------
        get モード: スカラー or SADList
        set モード: None

        Examples
        --------
        >>> ffs.LINE("K1", "QF")        # QF の K1 を取得
        >>> ffs.LINE("K1", "QF", 5.2)   # QF の K1 を 5.2 に設定
        >>> ffs.LINE("NAME", "Q*")      # Q で始まる全要素名を取得
        >>> ffs.LINE("POSITION", "B*")  # B で始まる全要素の s 座標
        """
        param_up = param.upper()

        # ---- NAME クエリ ----
        if param_up == 'NAME':
            indices = self._resolve_indices(element_arg)
            names = [self._names[i] for i in indices]
            return SADList(names) if names else SADList([])

        # ---- POSITION クエリ ----
        # SAD 互換: 1-based 整数インデックスを返す (s 座標ではない)
        # SAD の CalculateBeta[] 等は Twiss["AX", pos+1] のように使う
        if param_up == 'POSITION':
            indices = self._resolve_indices(element_arg)
            positions = [i + 1 for i in indices]   # 0-based → 1-based 整数
            return SADList(positions) if positions else SADList([])

        # ---- L クエリ ----
        if param_up == 'L':
            indices = self._resolve_indices(element_arg)
            lengths = [self._elements[i].L for i in indices]
            if value is not None:
                for i in indices:
                    self._elements[i].L = float(value)
                return None
            return lengths[0] if len(lengths) == 1 else SADList(lengths)

        # ---- K0 / K1 / K2 / ANGLE etc. ----
        attr = self._EL_ATTR_MAP.get(param_up, param_up)
        indices = self._resolve_indices(element_arg)

        # COD 関連パラメータ (Twiss 再計算は不要、COD のみ再計算)
        # ROTATION は転送行列を変更するため Twiss を完全無効化する (K1 等と同扱い)
        _COD_PARAMS = {'K0', 'ROTATE', 'DX', 'DY'}

        if value is not None:
            # set モード
            fval = float(value)
            for i in indices:
                el = self._elements[i]
                if hasattr(el, attr):
                    setattr(el, attr, fval)
                # 環境辞書も更新 (再ビルド時に反映される)
                raw = self._env._elements.get(el.name)
                if isinstance(raw, dict):
                    raw[param_up] = fval
            if param_up in _COD_PARAMS:
                # K0/ROTATE/DX/DY 変更時: Twiss は維持、COD のみ再計算
                if self._twiss:
                    self.calc_cod()
            else:
                # K1 等の変更: Twiss キャッシュを完全無効化
                self._twiss = []
                self._twiss0 = None
            return None

        # get モード
        results = []
        for i in indices:
            el = self._elements[i]
            val = getattr(el, attr, None)
            if val is None and hasattr(el, 'params'):
                val = el.params.get(param_up)
            results.append(val)

        if not results:
            return None
        return results[0] if len(results) == 1 else SADList(results)

    # ------------------------------------------------------------------
    # フィッティング
    # ------------------------------------------------------------------

    def free(self, pattern: str, param: str = 'K1') -> 'FFSState':
        """
        free Q* K1 — フィッティングのフリーパラメータを設定する。

        Parameters
        ----------
        pattern : str
            要素名のワイルドカードパターン
        param : str
            フリーにするパラメータ名 (デフォルト: K1)
        """
        self._free_params = []  # 毎回リセット (SAD の動作に倣う)
        indices = self._resolve_indices(pattern)
        for i in indices:
            self._free_params.append((i, self._names[i], param.upper()))
        return self

    def fit_add(self, element: str, param: str,
                target: float, nsteps: int = 3,
                npoints: Optional[int] = None) -> 'FFSState':
        """
        fit target を追加する。

        Parameters
        ----------
        element : str
            制約をかける要素名。 '$$$' で一周後 (チューン等)。
        param : str
            制約パラメータ (NX, NY, BX, BY, EX, ...)
        target : float
            目標値
        nsteps : int
            フィット手順ステップ数 (SAD の nsteps に相当)
        npoints : int or None
            運動量サンプル点数。None のとき fitp() で設定した値 (_fit_npoints) を使う。
            1 = 設計運動量のみ (通常)、2 = ±DP の 2 点、N = N 点均等配置。
            SAD の「FIT NX ncell 2」は npoints=2 に相当。

        Notes
        -----
        npoints > 1 かつ self._dp > 0 のとき、多点マッチングが有効になる。
        各運動量点で目標チューンが実現されるようにフィットする。
        """
        self._constraints.append({
            'element': element,
            'param': param.upper(),
            'target': float(target),
            'nsteps': nsteps,
            'npoints': npoints,   # None → go() 時に _fit_npoints を参照
        })
        return self

    def fit_clear(self) -> 'FFSState':
        """フィット制約をすべてクリア。"""
        self._constraints = []
        self._free_params = []
        return self

    def go(self, niter: int = 10) -> bool:
        """
        go — フィッターを実行する。

        scipy.optimize.minimize (Nelder-Mead) でフリーパラメータを
        最適化して制約を満たす。

        多点マッチング (fitp / fit_add の npoints > 1):
          制約ごとに指定された npoints 個の運動量点 (-DP ～ +DP) で
          チューン等の目標値を実現するようにフィットする。
          npoints=2 の場合: delta = -DP, +DP の 2 点で評価。
          これにより色収差を暗黙的に制御できる。

        Parameters
        ----------
        niter : int
            最大反復回数の係数 (内部の maxiter = niter * 200)

        Returns
        -------
        bool
            収束したら True
        """
        if not self._free_params:
            raise RuntimeError("go: free() でフリーパラメータを設定してください")
        if not self._constraints:
            raise RuntimeError("go: fit_add() でフィット制約を設定してください")

        from scipy.optimize import minimize

        # 初期値
        x0 = []
        for (i, name, param) in self._free_params:
            el = self._elements[i]
            x0.append(float(getattr(el, param, 0.0) or 0.0))

        def _set_params(x):
            for j, (i, name, param) in enumerate(self._free_params):
                setattr(self._elements[i], param, float(x[j]))
                raw = getattr(self._env, '_elements', {}).get(name)
                if isinstance(raw, dict):
                    raw[param] = float(x[j])

        def objective(x):
            _set_params(x)
            try:
                self.cell_calc()
            except Exception:
                return 1e12

            total = 0.0
            for c in self._constraints:
                # 制約ごとの npoints を解決 (None → グローバル設定)
                npts = c.get('npoints') or self._fit_npoints
                npts = max(1, int(npts))

                if npts > 1 and self._dp != 0.0:
                    # 多点マッチング: -DP ～ +DP を npts 点で評価
                    deltas = np.linspace(-self._dp, self._dp, npts)
                    for delta in deltas:
                        tw = self._cell_calc_at_delta(delta)
                        if tw is None:
                            total += 1e6
                            continue
                        val = self._eval_constraint_from_twiss(c, tw)
                        if val is None:
                            continue
                        total += (val - c['target']) ** 2
                else:
                    # 通常の 1 点マッチング (設計運動量)
                    val = self._eval_constraint(c)
                    if val is None:
                        continue
                    total += (val - c['target']) ** 2
            return total

        result = minimize(
            objective, x0,
            method='Nelder-Mead',
            options={
                'maxiter': niter * 200,
                'xatol': 1e-10, 'fatol': 1e-12,
                'adaptive': True,
            }
        )
        _set_params(result.x)
        self.cell_calc()
        return result.success

    def _cell_calc_at_delta(self, delta: float) -> Optional[dict]:
        """
        指定運動量偏差 delta = Δp/p での周期 Twiss を計算する。

        transfer_matrix(delta=delta) を各要素で呼び出し、
        periodic_optics() でチューン等を取得する。
        多点マッチングの内部ヘルパー。

        Parameters
        ----------
        delta : float
            運動量偏差 δ = Δp/p

        Returns
        -------
        dict or None
            Twiss パラメータ辞書 (bx, by, ax, ay, nx, ny, ...) or None (失敗時)
        """
        from .optics import periodic_optics

        M6 = np.eye(6)
        for el in self._elements:
            try:
                M6 = el.transfer_matrix(delta=delta) @ M6
            except TypeError:
                # delta 引数を受け付けない古い要素
                M6 = el.transfer_matrix() @ M6

        M4 = M6[:4, :4]
        try:
            tw = periodic_optics(M4)
            tw.setdefault('nx', 0.0)
            tw.setdefault('ny', 0.0)
            return tw
        except Exception:
            return None

    def _eval_constraint_from_twiss(self, c: dict, tw: Optional[dict]) -> Optional[float]:
        """
        Twiss 辞書から制約の現在値を評価する (多点マッチング用)。

        Parameters
        ----------
        c : dict   制約辞書
        tw : dict  Twiss パラメータ (bx, by, nx, ny, ...)

        Returns
        -------
        float or None
        """
        if tw is None:
            return None
        param = c['param']
        if param == 'NX':
            return tw.get('nx', 0.0) / (2 * math.pi)
        if param == 'NY':
            return tw.get('ny', 0.0) / (2 * math.pi)
        # Twiss パラメータの直接参照
        key = self._TWISS_KEY_MAP.get(param)
        if key:
            return tw.get(key)
        return None

    def _eval_constraint(self, c: dict) -> Optional[float]:
        """制約の現在値を評価する。"""
        param = c['param']
        elem  = c['element']

        # チューン系
        if elem in ('$$$', '***', None) or param in ('NX', 'NY'):
            tw = self._twiss[-1] if self._twiss else {}
            if param == 'NX':
                return tw.get('nx', 0.0) / (2 * math.pi)
            if param == 'NY':
                return tw.get('ny', 0.0) / (2 * math.pi)
        # Twiss パラメータ
        val = self.Twiss(param, elem)
        if isinstance(val, SADList):
            val = val[1]  # 最初の一致を使う
        return val

    # ------------------------------------------------------------------
    # SAVE / SAVE ALL — マッチング後の値をデフォルトとして格納
    # ------------------------------------------------------------------

    # 各要素タイプの「主変数」 (SAD initb1.f の kyparam 定義に基づく)
    # BEND の主変数: ANGLE (ky=2); K0 も保存 (実質的に ANGLE と対応)
    # QUAD の主変数: K1   (ky=2)
    # SEXT の主変数: K2
    # DRIFT の主変数: L
    _PRIMARY_PARAMS: dict = {}   # クラス変数として初期化は後で行う

    def save(self, all_params: bool = False) -> 'FFSState':
        """
        SAVE / SAVE ALL — マッチング計算で求めた値を新しいデフォルトとして格納する。

        SAD の FFS["SAVE"] / FFS["SAVE ALL"] に相当。
        フィッティング後に呼び出すことで、最適化後の要素パラメータ値を
        次回のエラーばらまき計算の基準値として確定させる。

        Parameters
        ----------
        all_params : bool
            False (SAVE):     各要素タイプの主変数のみ保存
            True  (SAVE ALL): 全パラメータを保存

        主変数の定義 (SAD kyparam モジュールに基づく):
            Drift      : L
            Bend       : ANGLE, K0 (両方)
            Quadrupole : K1
            Sextupole  : K2
            Octupole   : K3
            Cavity     : VOLT
            その他     : L のみ

        Returns
        -------
        self : FFSState (メソッドチェーン可能)

        Examples
        --------
        >>> ffs.free("SF*", "K2")
        >>> ffs.fit_add("$$$", "NX", 9.5)
        >>> ffs.go()
        >>> ffs.save()       # SF* の K2 を新しいデフォルトとして保存
        >>> ffs.save(all_params=True)   # 全要素の全パラメータを保存
        """
        from .elements import (Drift, Bend, Quadrupole, Sextupole,
                                Octupole, Cavity, Marker)

        # 主変数マッピング (クラス → パラメータ名リスト)
        primary_map = {
            Drift: ['L'],
            Bend: ['ANGLE', 'K0'],      # ANGLE が主; K0 は別名
            Quadrupole: ['K1'],
            Sextupole: ['K2'],
        }
        # Octupole, Cavity は importできない場合があるので try
        try:
            from .elements import Octupole
            primary_map[Octupole] = ['K3']
        except ImportError:
            pass
        try:
            from .elements import Cavity
            primary_map[Cavity] = ['VOLT']
        except ImportError:
            pass

        # 全パラメータリスト (SAVE ALL 用)
        all_attrs = [
            'L', 'K0', 'K1', 'K2', 'K3',
            'ANGLE', 'E1', 'E2',
            'DX', 'DY', 'DPX', 'DPY',
            'ROTATION', 'ROTATE', 'VOLT',
            'HGAP', 'FINT',
        ]

        self._saved_params = {}
        for el in self._elements:
            saved = {}
            if all_params:
                for attr in all_attrs:
                    if hasattr(el, attr):
                        saved[attr] = float(getattr(el, attr))
            else:
                # 主変数のみ
                primary = ['L']   # デフォルト: L のみ
                for cls, params in primary_map.items():
                    if isinstance(el, cls):
                        primary = params
                        break
                for attr in primary:
                    if hasattr(el, attr):
                        saved[attr] = float(getattr(el, attr))
            self._saved_params[el.name] = saved

        return self

    def restore_saved(self) -> 'FFSState':
        """
        save() / save(all_params=True) で保存した値を要素に復元する。

        エラーばらまき後に名目値に戻す場合などに使う。
        save() を呼んでいない場合は何もしない。

        Returns
        -------
        self : FFSState
        """
        if not self._saved_params:
            return self

        changed_twiss = False
        changed_cod = False
        _COD_ONLY = {'K0', 'ROTATE', 'DX', 'DY'}

        for el in self._elements:
            saved = self._saved_params.get(el.name, {})
            for attr, val in saved.items():
                if hasattr(el, attr):
                    setattr(el, attr, val)
                if attr.upper() in _COD_ONLY:
                    changed_cod = True
                else:
                    changed_twiss = True

        if changed_twiss:
            self._twiss = []
            self._twiss0 = None
        elif changed_cod and self._twiss:
            self.calc_cod()

        return self

    # ------------------------------------------------------------------
    # FITP — 多点マッチング設定
    # ------------------------------------------------------------------

    @property
    def dp(self) -> float:
        """DP — 多点マッチング用の運動量範囲 ±DP (Δp/p)。"""
        return self._dp

    @dp.setter
    def dp(self, value: float):
        """DP を設定する。"""
        self._dp = float(value)

    def fitp(self, n: int) -> 'FFSState':
        """
        FITP N — 以降のフィット計算で使う運動量サンプル点数を設定する。

        SAD の FFS["FITP N"] に相当。

        Parameters
        ----------
        n : int
            運動量点数:
              1 → 通常マッチング (設計運動量のみ)  ← FITP 1 で元に戻す
              2 → ±DP の 2 点 (chromaticity-aware)
              N → -DP ～ +DP を N 点均等配置

        Notes
        -----
        多点マッチングが有効になるには self.dp > 0 であることも必要。
        個々の fit_add() の npoints 引数が優先される (None のとき本設定を使う)。

        Examples
        --------
        >>> ffs.dp = 0.005          # ±0.5% の運動量範囲を設定
        >>> ffs.fitp(2)             # ±DP の 2 点でマッチング
        >>> ffs.fit_add("$$$", "NX", 9.5)
        >>> ffs.fit_add("$$$", "NY", 5.3)
        >>> ffs.go()
        >>> ffs.fitp(1)             # 1 点マッチングに戻す
        """
        self._fit_npoints = max(1, int(n))
        return self

    # ------------------------------------------------------------------
    # disp コマンド
    # ------------------------------------------------------------------

    def disp(self, element_arg=None, verbose: bool = True) -> Optional[dict]:
        """
        disp element — 要素位置での Twiss を表示する。

        Parameters
        ----------
        element_arg : str or None
            要素名・パターン。 None / '$$$' で現在のビームライン終端。
        verbose : bool
            True で標準出力にも表示。

        Returns
        -------
        dict : Twiss パラメータ辞書 (最初のマッチ)
        """
        if not self._twiss:
            self.cell_calc()

        if element_arg in (None, '$$$', '***'):
            tw = self._twiss0 or self._twiss[-1]
            label = self._beamline_name or 'RING'
        else:
            indices = self._resolve_indices(str(element_arg))
            if not indices:
                if verbose:
                    print(f"  [disp] 要素 '{element_arg}' が見つかりません")
                return None
            idx = indices[0]
            tw = self._twiss[idx + 1]
            label = self._names[idx]

        if verbose:
            bx, by = tw.get('bx', 0), tw.get('by', 0)
            ax, ay = tw.get('ax', 0), tw.get('ay', 0)
            ex, epx = tw.get('ex', 0), tw.get('epx', 0)
            nx = tw.get('nx', 0) / (2 * math.pi)
            ny = tw.get('ny', 0) / (2 * math.pi)
            print(f"  [{label}]")
            print(f"    BX={bx:12.6f}  BY={by:12.6f}")
            print(f"    AX={ax:12.6f}  AY={ay:12.6f}")
            print(f"    EX={ex:12.6f}  EPX={epx:12.6f}")
            print(f"    NX/2pi={nx:10.6f}  NY/2pi={ny:10.6f}")
        return tw

    # ------------------------------------------------------------------
    # 内部ユーティリティ
    # ------------------------------------------------------------------

    def _resolve_indices(self, element_arg) -> List[int]:
        """要素引数を元素インデックスのリストに変換する。"""
        if element_arg is None:
            return list(range(len(self._elements)))

        if isinstance(element_arg, (int, float)):
            idx = int(element_arg) - 1   # SAD 1-based
            if 0 <= idx < len(self._elements):
                return [idx]
            return []

        if isinstance(element_arg, str):
            if element_arg in _SPECIAL_ELEMENTS:
                return [len(self._elements) - 1]
            # ワイルドカードパターン
            return [i for i, n in enumerate(self._names)
                    if _match_pattern(n, element_arg)]

        # SADList または名前の list/tuple — 各要素を個別解決
        if isinstance(element_arg, (list, tuple, SADList)):
            result = []
            name_to_idx = {n: i for i, n in enumerate(self._names)}
            for item in element_arg:
                if isinstance(item, str):
                    idx = name_to_idx.get(item.upper())
                    if idx is None:
                        idx = name_to_idx.get(item)
                    if idx is not None:
                        result.append(idx)
                elif isinstance(item, (int, float)):
                    idx = int(item) - 1
                    if 0 <= idx < len(self._elements):
                        result.append(idx)
            return result

        return []

    # ------------------------------------------------------------------
    # プロパティ / 便利メソッド
    # ------------------------------------------------------------------

    @property
    def EMITX(self) -> Optional[float]:
        """水平エミッタンス [m·rad] (emit_calc() 後に有効)。"""
        return self._emitx

    @property
    def EMITY(self) -> Optional[float]:
        """垂直エミッタンス [m·rad]。"""
        return self._emity

    @property
    def SIGE(self) -> Optional[float]:
        """エネルギー広がり (σ_E/E)。"""
        return self._sige

    def summary(self) -> str:
        """現在の状態のサマリーを返す。"""
        lines = [f"FFSState — beamline: {self._beamline_name}"]
        if self._twiss0:
            tw = self._twiss0
            lines.append(
                f"  BX={tw['bx']:.4f}  BY={tw['by']:.4f}"
                f"  AX={tw['ax']:.4f}  AY={tw['ay']:.4f}"
            )
            nx = tw.get('nx', 0) / (2 * math.pi)
            ny = tw.get('ny', 0) / (2 * math.pi)
            lines.append(f"  Qx_frac={nx:.4f}  Qy_frac={ny:.4f}")
        if self._emitx is not None:
            lines.append(
                f"  EMITX={self._emitx:.4e} m·rad"
                f"  EMITY={self._emity:.4e} m·rad"
            )
        return '\n'.join(lines)

    def __repr__(self) -> str:
        return self.summary()

    # ------------------------------------------------------------------
    # 誤差ばらまき (Error Spreading)
    # SAD の FFS["DELX ae Q*;"] 等に相当
    # ------------------------------------------------------------------

    def apply_misalignment(self,
                           pattern: str,
                           dx_sigma: float = 0.0,
                           dy_sigma: float = 0.0,
                           gcut: float = 3.0,
                           seed: Optional[int] = None) -> 'FFSState':
        """
        アライメント誤差 (DX/DY) を指定パターンの要素にばらまく。

        SAD の FFS["DELX ae Q*;"] / FFS["DELY ae Q*;"] に相当。
        各要素に個別の独立したランダム値を設定する (モンテカルロ試行 1 回分)。

        Parameters
        ----------
        pattern   : str    要素名のワイルドカードパターン (例: "Q*", "B*")
        dx_sigma  : float  水平オフセット DX の標準偏差 [m]
        dy_sigma  : float  垂直オフセット DY の標準偏差 [m]
        gcut      : float  打ち切りガウス分布のカット (デフォルト 3σ)
                           SAD の GCUT パラメータに相当
        seed      : int    乱数シード (None なら自動)

        Returns
        -------
        self : FFSState (メソッドチェーン可能)

        Examples
        --------
        >>> ffs.apply_misalignment("Q*", dx_sigma=1e-4, dy_sigma=1e-4, gcut=3)
        >>> ffs.apply_misalignment("B*", dx_sigma=5e-5, dy_sigma=5e-5)
        """
        from scipy.stats import truncnorm

        rng = np.random.default_rng(seed)
        indices = self._resolve_indices(pattern)

        for i in indices:
            el = self._elements[i]
            if dx_sigma > 0.0:
                val = float(truncnorm.rvs(-gcut, gcut,
                                         scale=dx_sigma,
                                         random_state=rng))
                el.DX = val
            if dy_sigma > 0.0:
                val = float(truncnorm.rvs(-gcut, gcut,
                                         scale=dy_sigma,
                                         random_state=rng))
                el.DY = val

        # DX/DY 変更 → COD のみ再計算 (Twiss は維持)
        if self._twiss:
            self.calc_cod()
        return self

    def apply_field_error(self,
                          pattern: str,
                          param: str = 'K1',
                          sigma: float = 0.0,
                          relative: bool = True,
                          gcut: float = 3.0,
                          seed: Optional[int] = None) -> 'FFSState':
        """
        強さ誤差 (K0/K1/K2 の相対/絶対誤差) を指定パターンの要素にばらまく。

        SAD の FFS["DK fe B*;"] (二極偏向角誤差) に相当。
        各要素に個別の独立したランダム値を設定する。

        Parameters
        ----------
        pattern  : str    要素名のワイルドカードパターン
        param    : str    対象パラメータ名 ('K0', 'K1', 'K2')
        sigma    : float  誤差の標準偏差
                          relative=True  なら σ は相対値 (例: 1e-3 → 0.1%)
                          relative=False なら σ は絶対値 [m⁻ⁿ]
        relative : bool   True: 相対誤差 K_new = K * (1 + σ·ε),  False: 絶対誤差
        gcut     : float  打ち切りガウス分布のカット (デフォルト 3σ)
        seed     : int    乱数シード

        Returns
        -------
        self : FFSState

        Examples
        --------
        >>> ffs.apply_field_error("B*", param='K0', sigma=1e-3, relative=True)
        >>> ffs.apply_field_error("Q*", param='K1', sigma=5e-4, relative=True)
        """
        from scipy.stats import truncnorm

        rng = np.random.default_rng(seed)
        attr = self._EL_ATTR_MAP.get(param.upper(), param.upper())
        indices = self._resolve_indices(pattern)

        for i in indices:
            el = self._elements[i]
            if not hasattr(el, attr):
                continue
            k_nom = float(getattr(el, attr) or 0.0)
            eps = float(truncnorm.rvs(-gcut, gcut, scale=sigma,
                                      random_state=rng))
            if relative:
                k_new = k_nom * (1.0 + eps)
            else:
                k_new = k_nom + eps
            setattr(el, attr, k_new)

        # K0/K1/K2 変更 → Twiss キャッシュを完全無効化
        self._twiss = []
        self._twiss0 = None
        return self

    def apply_rotation(self,
                       pattern: str,
                       sigma: float = 0.0,
                       gcut: float = 3.0,
                       seed: Optional[int] = None) -> 'FFSState':
        """
        ロール誤差 (ROTATION) を指定パターンの要素にばらまく。

        SAD の FFS["DTHETA te B*;"] に相当。
        各要素の ROTATION 属性に個別の独立したランダム値を設定する。
        ROTATION は据え付けロール誤差で、全要素の transfer_matrix() で
        R(-θ)·M·R(θ) として適用される。

        Parameters
        ----------
        pattern : str    要素名のワイルドカードパターン
        sigma   : float  ロール誤差の標準偏差 [rad]
        gcut    : float  打ち切りガウス分布のカット (デフォルト 3σ)
        seed    : int    乱数シード

        Returns
        -------
        self : FFSState

        Examples
        --------
        >>> ffs.apply_rotation("B*", sigma=1e-3)   # 偏向電磁石に 1 mrad ロール誤差
        >>> ffs.apply_rotation("Q*", sigma=5e-4)   # 四極に 0.5 mrad ロール誤差
        """
        from scipy.stats import truncnorm

        rng = np.random.default_rng(seed)
        indices = self._resolve_indices(pattern)

        for i in indices:
            el = self._elements[i]
            if hasattr(el, 'ROTATION'):
                val = float(truncnorm.rvs(-gcut, gcut, scale=sigma,
                                          random_state=rng))
                el.ROTATION = val

        # ROTATION は転送行列を変更するため Twiss キャッシュを完全無効化
        self._twiss = []
        self._twiss0 = None
        return self

    def reset_errors(self,
                     pattern: Optional[str] = None,
                     dx: bool = True,
                     dy: bool = True,
                     rotation: bool = True,
                     field: Optional[List[str]] = None) -> 'FFSState':
        """
        誤差パラメータをゼロにリセットする。

        Parameters
        ----------
        pattern  : str or None  対象要素パターン (None なら全要素)
        dx       : bool  DX (水平設置誤差) をリセット
        dy       : bool  DY (垂直設置誤差) をリセット
        rotation : bool  ROTATION (ロール誤差) をリセット
        field    : list or None  強さ誤差をリセットするパラメータ名リスト
                   例: ['K0', 'K1'] — ただし nominal 値が不明なため 0 にリセット

        Returns
        -------
        self : FFSState
        """
        indices = self._resolve_indices(pattern)
        for i in indices:
            el = self._elements[i]
            if dx and hasattr(el, 'DX'):
                el.DX = 0.0
            if dy and hasattr(el, 'DY'):
                el.DY = 0.0
            if rotation and hasattr(el, 'ROTATION'):
                el.ROTATION = 0.0
            if field:
                # field は ['K0', 'K1', 'K2'] 等のリスト
                for param in field:
                    if hasattr(el, param.upper()):
                        setattr(el, param.upper(), 0.0)

        # Twiss / COD キャッシュを適切に無効化
        if rotation or field:
            self._twiss = []
            self._twiss0 = None
        elif (dx or dy) and self._twiss:
            self.calc_cod()
        return self
