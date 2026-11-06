"""
pysad.matching
==============
SAD の Match[] / Fit[] コマンドに相当する光学マッチング。

使い方
------
>>> from pysad import Lattice, Match
>>> lat = fodo_cell()
>>> m = Match(lat)
>>> m.vary('QF', 'K1', step=0.01)
>>> m.vary('QD', 'K1', step=0.01)
>>> m.constraint_tune(Qx=0.25, Qy=0.25)
>>> result = m.minimize()
>>> print(result)

対応する SAD コマンド
---------------------
  Match[LINE, {Vary["QF", "K1", 0.1, 5.0], ...},
              {Constraint[Tune[LINE] == {0.25, 0.25}], ...}]
"""

from __future__ import annotations

import copy
import numpy as np
from typing import Callable, Dict, List, Optional, Tuple, Any

try:
    from scipy.optimize import minimize as _scipy_minimize
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


# ---------------------------------------------------------------------------
# Vary 記述子
# ---------------------------------------------------------------------------

class Vary:
    """
    マッチングで変化させる変数の定義。

    Parameters
    ----------
    element_name : str   要素名
    param        : str   パラメータ名 ('K1', 'ANGLE', 'VOLT', ...)
    lower        : float 下限 (None で無制限)
    upper        : float 上限 (None で無制限)
    step         : float 数値微分のステップ幅
    """
    def __init__(self, element_name: str, param: str,
                 lower: float = None, upper: float = None,
                 step: float = 1e-4):
        self.element_name = element_name
        self.param = param.upper()
        self.lower = lower
        self.upper = upper
        self.step = step

    def __repr__(self):
        return (f"Vary('{self.element_name}', '{self.param}', "
                f"lower={self.lower}, upper={self.upper})")


# ---------------------------------------------------------------------------
# Constraint 記述子
# ---------------------------------------------------------------------------

class Constraint:
    """
    マッチングの制約条件。

    Parameters
    ----------
    func   : callable(lattice) -> float or array
             制約値を計算する関数
    target : float or array
             目標値
    weight : float
             ペナルティ重み (デフォルト 1.0)
    label  : str
             表示用ラベル
    """
    def __init__(self, func: Callable, target, weight: float = 1.0,
                 label: str = ''):
        self.func = func
        self.target = np.atleast_1d(np.asarray(target, dtype=float))
        self.weight = weight
        self.label = label

    def residual(self, lat) -> np.ndarray:
        val = np.atleast_1d(np.asarray(self.func(lat), dtype=float))
        return self.weight * (val - self.target)

    def __repr__(self):
        return f"Constraint(target={self.target}, weight={self.weight}, label={self.label!r})"


# ---------------------------------------------------------------------------
# マッチング結果
# ---------------------------------------------------------------------------

class MatchResult:
    """
    Match.minimize() の結果コンテナ。

    Attributes
    ----------
    success      : bool
    lattice      : Lattice  マッチ後の格子
    residuals    : dict     各制約の残差
    niter        : int      反復回数
    fun          : float    最終目的関数値
    message      : str
    """
    def __init__(self, success, lattice, residuals, niter, fun, message):
        self.success = success
        self.lattice = lattice
        self.residuals = residuals
        self.niter = niter
        self.fun = fun
        self.message = message

    def __repr__(self):
        return (f"MatchResult(success={self.success}, fun={self.fun:.3e}, "
                f"niter={self.niter})")

    def print_summary(self):
        print(f"{'Match result':=^50}")
        print(f"  Status : {'OK' if self.success else 'FAILED'}")
        print(f"  Merit  : {self.fun:.4e}")
        print(f"  Iters  : {self.niter}")
        print(f"  Msg    : {self.message}")
        if self.residuals:
            print("  Residuals:")
            for lbl, r in self.residuals.items():
                print(f"    {lbl:20s}: {r}")


# ---------------------------------------------------------------------------
# Match クラス (メイン)
# ---------------------------------------------------------------------------

class Match:
    """
    SAD の Match[] コマンドに相当する光学マッチング。

    使い方
    ------
    >>> m = Match(lattice)
    >>> m.vary('QF', 'K1', lower=0.1, upper=5.0)
    >>> m.vary('QD', 'K1', lower=-5.0, upper=-0.1)
    >>> m.constraint_tune(Qx=0.25, Qy=0.25)
    >>> m.constraint_beta(element_idx=0, bx=5.0, by=5.0)
    >>> result = m.minimize()
    """

    def __init__(self, lattice):
        self._lat0 = lattice          # 元格子 (変更しない)
        self._vary_list: List[Vary] = []
        self._constr_list: List[Constraint] = []

    # ------------------------------------------------------------------
    # 変数の登録
    # ------------------------------------------------------------------

    def vary(self, element_name: str, param: str,
             lower: float = None, upper: float = None,
             step: float = 1e-4) -> 'Match':
        """変化させる要素パラメータを登録する。"""
        self._vary_list.append(
            Vary(element_name, param, lower, upper, step))
        return self

    # ------------------------------------------------------------------
    # 制約の登録 (プリセット)
    # ------------------------------------------------------------------

    def constraint_tune(self, Qx: float = None, Qy: float = None,
                        weight: float = 1.0) -> 'Match':
        """チューン Qx, Qy の制約。"""
        targets = []
        components = []
        if Qx is not None:
            targets.append(Qx)
            components.append('x')
        if Qy is not None:
            targets.append(Qy)
            components.append('y')

        def _tune_func(lat):
            tw = lat.twiss(full=False)
            tw_end = lat.twiss()[-1]
            vals = []
            if 'x' in components:
                vals.append(tw_end['nx'] / (2 * np.pi))
            if 'y' in components:
                vals.append(tw_end['ny'] / (2 * np.pi))
            return vals

        self._constr_list.append(
            Constraint(_tune_func, targets, weight,
                       label=f"Tune({'Qx' if Qx else ''}"
                             f"{'Qy' if Qy else ''})"))
        return self

    def constraint_beta(self, element_idx: int,
                        bx: float = None, by: float = None,
                        weight: float = 1.0) -> 'Match':
        """特定要素位置での β 関数の制約。"""
        targets = []
        components = []
        if bx is not None:
            targets.append(bx)
            components.append('bx')
        if by is not None:
            targets.append(by)
            components.append('by')

        def _beta_func(lat):
            tw = lat.twiss()
            tw_at = tw[element_idx]
            return [tw_at[c] for c in components]

        self._constr_list.append(
            Constraint(_beta_func, targets, weight,
                       label=f"Beta[{element_idx}]"))
        return self

    def constraint_dispersion(self, element_idx: int,
                               eta: float = None, etap: float = None,
                               weight: float = 1.0) -> 'Match':
        """特定位置での分散 η, η' の制約。"""
        targets, comps = [], []
        if eta is not None:
            targets.append(eta); comps.append('ex')
        if etap is not None:
            targets.append(etap); comps.append('epx')

        def _disp_func(lat):
            tw = lat.twiss()
            return [tw[element_idx][c] for c in comps]

        self._constr_list.append(
            Constraint(_disp_func, targets, weight,
                       label=f"Dispersion[{element_idx}]"))
        return self

    def constraint_chromaticity(self, Cx: float = None, Cy: float = None,
                                 weight: float = 1.0) -> 'Match':
        """色収差 ξx, ξy の制約。"""
        targets, comps = [], []
        if Cx is not None:
            targets.append(Cx); comps.append('chrom_x')
        if Cy is not None:
            targets.append(Cy); comps.append('chrom_y')

        def _chrom_func(lat):
            s = lat.summary()
            return [s[c] for c in comps]

        self._constr_list.append(
            Constraint(_chrom_func, targets, weight,
                       label=f"Chromaticity"))
        return self

    def constraint_custom(self, func: Callable, target,
                           weight: float = 1.0, label: str = '') -> 'Match':
        """任意の制約関数を登録する。"""
        self._constr_list.append(
            Constraint(func, target, weight, label=label))
        return self

    # ------------------------------------------------------------------
    # 格子 <-> パラメータベクトル の変換
    # ------------------------------------------------------------------

    def _get_params(self, lat) -> np.ndarray:
        """現在の格子から vary パラメータを読み出す。"""
        vals = []
        for v in self._vary_list:
            el = self._find_element(lat, v.element_name)
            if el is not None:
                vals.append(float(getattr(el, v.param, 0.0)))
            else:
                vals.append(0.0)
        return np.array(vals)

    def _set_params(self, lat, params: np.ndarray):
        """vary パラメータをまとめて格子に設定する (格子をコピー済みの前提)。"""
        for v, p in zip(self._vary_list, params):
            el = self._find_element(lat, v.element_name)
            if el is not None:
                setattr(el, v.param, float(p))

    @staticmethod
    def _find_element(lat, name: str):
        for el in lat.elements:
            if el.name == name:
                return el
        return None

    def _bounds(self) -> list:
        return [(v.lower, v.upper) for v in self._vary_list]

    # ------------------------------------------------------------------
    # 目的関数
    # ------------------------------------------------------------------

    def _merit(self, params: np.ndarray) -> float:
        """制約残差の二乗和。"""
        lat = copy.deepcopy(self._lat0)
        self._set_params(lat, params)
        total = 0.0
        try:
            for c in self._constr_list:
                r = c.residual(lat)
                total += float(np.sum(r**2))
        except Exception:
            total = 1e30
        return total

    # ------------------------------------------------------------------
    # 最適化
    # ------------------------------------------------------------------

    def minimize(self, method: str = 'auto',
                 tol: float = 1e-8,
                 maxiter: int = 2000,
                 verbose: bool = False) -> MatchResult:
        """
        マッチングを実行する。

        Parameters
        ----------
        method  : 'auto' | 'nelder-mead' | 'cobyla' | 'lbfgsb' | 'simplex'
                  'auto' → scipy があれば L-BFGS-B (境界あり), なければ独自 Simplex
        tol     : 収束許容誤差
        maxiter : 最大反復回数
        verbose : True で各反復を表示

        Returns
        -------
        MatchResult
        """
        x0 = self._get_params(self._lat0)
        bounds = self._bounds()

        if method == 'auto':
            method = 'lbfgsb' if _HAS_SCIPY else 'simplex'

        if _HAS_SCIPY and method in ('lbfgsb', 'L-BFGS-B'):
            scipy_bounds = [(lo if lo is not None else -np.inf,
                             hi if hi is not None else  np.inf)
                            for lo, hi in bounds]
            res = _scipy_minimize(
                self._merit, x0, method='L-BFGS-B',
                bounds=scipy_bounds,
                options={'maxiter': maxiter, 'ftol': tol, 'gtol': tol*1e-3,
                         'disp': verbose})
            x_final = res.x
            success  = res.success
            niter    = res.nit
            message  = res.message

        elif _HAS_SCIPY and method in ('nelder-mead', 'nelder'):
            res = _scipy_minimize(
                self._merit, x0, method='Nelder-Mead',
                options={'maxiter': maxiter, 'xatol': tol, 'fatol': tol,
                         'disp': verbose})
            x_final = res.x
            success  = res.success
            niter    = res.nit
            message  = res.message

        elif _HAS_SCIPY and method == 'cobyla':
            constrs = [
                {'type': 'ineq', 'fun': lambda x, lo=lo: x[i] - lo}
                for i, (lo, hi) in enumerate(bounds) if lo is not None
            ] + [
                {'type': 'ineq', 'fun': lambda x, hi=hi: hi - x[i]}
                for i, (lo, hi) in enumerate(bounds) if hi is not None
            ]
            res = _scipy_minimize(
                self._merit, x0, method='COBYLA',
                constraints=constrs,
                options={'maxiter': maxiter, 'rhobeg': 0.1, 'disp': verbose})
            x_final = res.x
            success  = res.success
            niter    = res.nfev
            message  = res.message

        else:
            # 独自 Nelder-Mead (scipy 不要)
            x_final, success, niter, message = _nelder_mead(
                self._merit, x0, bounds, tol=tol, maxiter=maxiter,
                verbose=verbose)

        # ---- 結果格子を構築 ----
        lat_final = copy.deepcopy(self._lat0)
        self._set_params(lat_final, x_final)

        residuals = {}
        for c in self._constr_list:
            try:
                residuals[c.label or str(c)] = c.residual(lat_final)
            except Exception:
                pass

        fun = self._merit(x_final)

        return MatchResult(
            success=success,
            lattice=lat_final,
            residuals=residuals,
            niter=niter,
            fun=fun,
            message=message,
        )


# ---------------------------------------------------------------------------
# 独自 Nelder-Mead (scipy 不要バックアップ)
# ---------------------------------------------------------------------------

def _nelder_mead(func, x0, bounds=None, tol=1e-8, maxiter=2000,
                 verbose=False):
    """
    シンプルな Nelder-Mead 実装。scipy がない環境用。
    """
    n = len(x0)
    alpha = 1.0; beta = 0.5; gamma = 2.0; sigma = 0.5

    # 初期シンプレックス
    simplex = [x0.copy()]
    for i in range(n):
        xi = x0.copy()
        xi[i] += max(0.05 * abs(x0[i]), 1e-4)
        simplex.append(xi)

    def _clip(x):
        if bounds is None:
            return x
        y = x.copy()
        for i, (lo, hi) in enumerate(bounds):
            if lo is not None: y[i] = max(y[i], lo)
            if hi is not None: y[i] = min(y[i], hi)
        return y

    fvals = [func(_clip(s)) for s in simplex]
    niter = 0

    for niter in range(maxiter):
        # ソート
        order = np.argsort(fvals)
        simplex = [simplex[i] for i in order]
        fvals   = [fvals[i]   for i in order]

        if verbose and niter % 100 == 0:
            print(f"[Match] iter={niter} f={fvals[0]:.4e}")

        # 収束チェック
        spread = max(abs(f - fvals[0]) for f in fvals[1:])
        if spread < tol:
            break

        centroid = np.mean(simplex[:-1], axis=0)
        worst    = simplex[-1]

        # Reflection
        xr = _clip(centroid + alpha * (centroid - worst))
        fr = func(xr)

        if fr < fvals[0]:
            # Expansion
            xe = _clip(centroid + gamma * (xr - centroid))
            fe = func(xe)
            if fe < fr:
                simplex[-1] = xe; fvals[-1] = fe
            else:
                simplex[-1] = xr; fvals[-1] = fr
        elif fr < fvals[-2]:
            simplex[-1] = xr; fvals[-1] = fr
        else:
            # Contraction
            xc = _clip(centroid + beta * (worst - centroid))
            fc = func(xc)
            if fc < fvals[-1]:
                simplex[-1] = xc; fvals[-1] = fc
            else:
                # Shrink
                best = simplex[0]
                for i in range(1, n + 1):
                    simplex[i] = _clip(best + sigma * (simplex[i] - best))
                    fvals[i]   = func(simplex[i])

    success = fvals[0] < tol * 1e3
    message = "Converged" if success else f"Max iterations ({maxiter}) reached"
    return _clip(simplex[0]), success, niter, message


# ---------------------------------------------------------------------------
# DynamicApertureSurvey クラス
# ---------------------------------------------------------------------------

class DynamicApertureSurvey:
    """
    SAD の DynamicAperture[] の強化版。
    運動量偏差 dp をスキャンして動的アパーチャの dp 依存性を調べる。

    使い方
    ------
    >>> das = DynamicApertureSurvey(lat, nturns=256, dp_range=(-0.03, 0.03), ndp=7)
    >>> das.run()
    >>> das.print_summary()
    >>> data = das.result()  # dict: dp -> {'x': [...], 'y': [...]}
    """

    def __init__(self, lattice,
                 nturns: int = 1024,
                 nangles: int = 9,
                 x_max: float = 0.05,
                 dp_range: Tuple[float, float] = (0.0, 0.0),
                 ndp: int = 1,
                 aperture: float = 0.1,
                 exact: bool = True):
        """
        Parameters
        ----------
        lattice  : Lattice
        nturns   : int   各振幅探索の周回数
        nangles  : int   方位角分割数 (0 〜 π/2)
        x_max    : float 最大探索振幅 [m]
        dp_range : (dp_min, dp_max)  運動量偏差のスキャン範囲
        ndp      : int   dp スキャン点数
        aperture : float アパーチャ半径 [m]
        exact    : bool  厳密トラッキング使用
        """
        self._lat = lattice
        self.nturns  = nturns
        self.nangles = nangles
        self.x_max   = x_max
        self.aperture = aperture
        self.exact   = exact
        self.dp_values = (np.linspace(dp_range[0], dp_range[1], ndp)
                          if ndp > 1 else np.array([dp_range[0]]))
        self._results: Dict[float, dict] = {}

    def run(self, verbose: bool = True) -> 'DynamicApertureSurvey':
        """動的アパーチャサーベイを実行する。"""
        from .tracking import dynamic_aperture as _da
        for dp in self.dp_values:
            if verbose:
                print(f"[DAS] dp={dp:+.4f} ... ", end='', flush=True)
            res = _da(self._lat.elements, nturns=self.nturns,
                      nangles=self.nangles, x_max=self.x_max,
                      dp=dp, exact=self.exact)
            self._results[float(dp)] = res
            if verbose:
                a_max = res['amplitudes'].max()
                print(f"A_max={a_max*1e3:.2f} mm")
        return self

    def result(self) -> Dict[float, dict]:
        """
        結果辞書を返す。
        {dp_value: {'angles': [...], 'amplitudes': [...], 'x': [...], 'y': [...]}}
        """
        return dict(self._results)

    def max_aperture(self) -> Dict[float, float]:
        """各 dp での最大振幅 [m] を返す。"""
        return {dp: r['amplitudes'].max()
                for dp, r in self._results.items()}

    def print_summary(self):
        print(f"{'Dynamic Aperture Survey':=^50}")
        print(f"  nturns={self.nturns}, nangles={self.nangles}")
        print(f"  {'dp':>8}  {'Amax[mm]':>10}  {'Ax[mm]':>9}  {'Ay[mm]':>9}")
        for dp, res in sorted(self._results.items()):
            A  = res['amplitudes']
            ax = res['x']
            ay = res['y']
            print(f"  {dp:>+8.4f}  {A.max()*1e3:>10.3f}"
                  f"  {ax.max()*1e3:>9.3f}  {ay.max()*1e3:>9.3f}")

    def to_array(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        (dp_array, angles, amplitudes_2d) を返す。
        amplitudes_2d.shape = (ndp, nangles)
        """
        dps = sorted(self._results.keys())
        angles = next(iter(self._results.values()))['angles']
        amps = np.array([self._results[dp]['amplitudes'] for dp in dps])
        return np.array(dps), angles, amps


__all__ = [
    'Vary', 'Constraint', 'Match', 'MatchResult',
    'DynamicApertureSurvey',
]
