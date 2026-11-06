"""
match6.py — SAD の FIT / FREE / FIX / GO (tffsmatch.f, tffscalc.f, qdcell.f tdfun, tffsmatch.f tweigh/tfsolv/tgrad,
tgetfv.f, tffsfreefix.f) の直訳によるマッチング

SAD の流れ (tffsmatch):
  1. tffscalc: 現在の変数で光学 (qcell1) を計算し、FIT 条件ごとの残差 df (tdfun: 目標 − 現在値、BX/BY は log、
     AX/AY は atan) に重み wiq (twfit + FitWeight フック) を掛け、r = Σ|df|^w (w = ExponentOfResidual = 2) と
     不安定数 nstab を出す。
  2. 各変数の微分: 数値微分 (FFS$NumericalDerivative の経路): val ± dvkc で tffscalc を 2 回。
     qu[:,kc] = (df(−) − df(+))/2/dvkc/wvar (wvar = tweigh の重み)。
  3. Newton: tfsolv (SVD, 打ち切り tstol/fact) で quw·dval = df を解き、dg<0 なら更新
     val = bestval + dval·fact/wvar·wlimit。改善しなければ fmincube で fact を縮める。
     収束が悪いと DESCEND (tgrad の勾配) に切り替える (fuzz 論理)。
  4. 収束: res2r(r) ≤ CONVERGENCE·max(nfcol,1) → "Matched."。itmax (MAXITERATION=40)。

SAD の解析微分 (tffsqu / qdcell の dtwiss) は使わず数値微分で代用する (SAD でも FFS$NumericalDerivative=True で
同じ経路になる)。結果は収束解として一致するが、反復の経路は厳密には違い得る。
対応: CELL のリング (INS は未対応)、FITP 1 (オンモメンタムのみ)、条件 AX/BX/NX/AY/BY/NY/EX/EPX/EY/EPY/LENG、
MAX フィット ('M' 付き)、範囲フィット {min,max}、2 要素の相対フィット、FitValue/FitWeight フック。
変数: 要素ファミリーのキーワード (既定: DRIFT L, BEND ANGLE, QUAD K1, SEXT K2, OCT K3, CAVI VOLT)。
値は SAD の単位 (K1, K2, K3 は積分値)。
"""
from __future__ import annotations

import fnmatch
import math
from typing import Callable, Dict, List, Optional

import numpy as np

from . import lattice as _lat
from .lattice import TwissRow
from .sadscript import fmt_g

import os as _os
_DEBUG = bool(_os.environ.get('PYSAD_MATCH_DEBUG'))
PI2 = 2.0 * math.pi
# tffsmatch.f module match のパラメータ
FLIM1, FLIM2 = -4.0, -3.0
AIMP1, AIMP2 = -1.8, -0.8
BADC1, BADC2 = -3.5, -2.5
AMEDC1, AMEDC2 = -2.3, -1.3
ALIT = 0.75
WLMIN = 0.009
EPS = 1e-5
EPS1 = 1e-8
RTOL = 1.05
RTOL1 = 1.02
TSTOL = 1e-6
AMTOL = 1e-9
M_EULER = 0.57721566490153286061   # tfstk.f の m_euler はオイラー定数 γ (ネイピア数 e ではない)。fmincube の上限 min(m_euler·f1, …)。以前は e を入れていた

# tffitcode の関数名 (nlist)
_FUNS = ['AX', 'BX', 'NX', 'AY', 'BY', 'NY', 'EX', 'EPX', 'EY', 'EPY', 'R1', 'R2', 'R3', 'R4', 'DETR',
         'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AZ', 'BZ', 'NZ', 'ZX', 'ZPX', 'ZY', 'ZPY',
         'PEX', 'PEPX', 'PEY', 'PEPY', 'PZX', 'PZPX', 'PZY', 'PZPY', 'GMX', 'GMY', 'GMZ',
         'BMAGX', 'BMAGY', 'BMAGZ', 'TRX', 'TRY', 'TRZ', 'LENG', 'GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3']
_LOGFIT = ('BX', 'BY', 'BZ', 'GMX', 'GMY', 'GMZ')
_GEO_KF = ('GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3')      # 幾何の条件 (tgfun は geo から、微分は数値微分)


def _normal_to_physical(tw: dict) -> np.ndarray:
    """tfnormaltophysical (tffs.f): 正規モード座標 → 物理座標の 4×4 (R1-R4、DETR<1 で分岐)"""
    r1, r2, r3, r4 = (float(tw.get(k, 0.0) or 0.0) for k in ('r1', 'r2', 'r3', 'r4'))
    detr = r1 * r4 - r2 * r3
    cc = math.sqrt(max(0.0, 1.0 - detr))
    if float(tw.get('detr', 0.0) or 0.0) < 1.0:
        return np.array([[cc, 0.0, r4, -r2], [0.0, cc, -r3, r1], [-r1, -r2, cc, 0.0], [-r3, -r4, 0.0, cc]])
    return np.array([[-r1, -r2, cc, 0.0], [-r3, -r4, 0.0, cc], [cc, 0.0, r4, -r2], [0.0, cc, -r3, r1]])
_ATANFIT = ('AX', 'AY', 'AZ')
_ROW_ATTR = {'AX': 'alfx', 'BX': 'betx', 'NX': 'mux', 'AY': 'alfy', 'BY': 'bety', 'NY': 'muy',
             'EX': 'etx', 'EPX': 'etpx', 'EY': 'ety', 'EPY': 'etpy', 'R1': 'r1', 'R2': 'r2', 'R3': 'r3', 'R4': 'r4',
             'DETR': 'detr', 'DX': 'x_cod', 'DPX': 'px_cod', 'DY': 'y_cod', 'DPY': 'py_cod', 'DZ': 'z_cod',
             'DDP': 'delta_cod', 'LENG': 's'}

# 既定の変数キーワード (tfinit.f: ival=2 は各要素型の 2 番目のキーワード)
_DEFAULT_KEY = {'Drift': 'L', 'Bend': 'ANGLE', 'Quadrupole': 'K1', 'Sextupole': 'K2', 'Octupole': 'K3',
                'Decapole': 'K4', 'Cavity': 'VOLT', 'Solenoid': 'BZ', 'Multipole': 'K1'}
# SAD の値 (積分値) ↔ pysad の属性 (長さあたり) の変換が要るもの
_PER_LENGTH = {('Quadrupole', 'K1'), ('Bend', 'K1'), ('Sextupole', 'K2'), ('Octupole', 'K3'), ('Decapole', 'K4')}


def _sad_get(el, key):
    tn = type(el).__name__
    if tn == 'Multipole' and key[0] == 'K' and key[1:].isdigit():
        return float(getattr(el, 'KnL', {}).get(int(key[1:]), 0.0))
    v = float(getattr(el, key, 0.0) or 0.0)
    if (tn, key) in _PER_LENGTH and el.L != 0.0:
        v *= el.L
    return v


def _sad_set(el, key, v):
    tn = type(el).__name__
    if tn == 'Multipole' and key[0] == 'K' and key[1:].isdigit():
        el.KnL[int(key[1:])] = v
        return
    if (tn, key) in _PER_LENGTH and el.L != 0.0:
        v = v / el.L
    setattr(el, key, v)


class _Res:
    """ffs_res: (r, nstab)"""
    __slots__ = ('r', 'nstab')

    def __init__(self, r=0.0, nstab=0):
        self.r = r
        self.nstab = nstab

    def le(self, other, tol=None):
        """resle"""
        if tol is None:
            return self.nstab < other.nstab or (self.nstab == other.nstab and self.r <= other.r)
        return self.nstab < other.nstab or (self.nstab == other.nstab and self.r <= other.r / tol)

    def res2r(self):
        rb = 10.0
        while rb < self.r:
            rb *= 10.0
        return rb * self.nstab + self.r


def _fuzz(x, a, b):
    """tffsmatch の fuzz。x が NaN (log10 の負の引数) のときは gfortran の MIN/MAX が NaN でない方を返すので 1 になる"""
    if x != x:
        return 1.0
    return max(0.0, min(1.0, (x - a) / (b - a)))


def _log10(x):
    """Fortran の log10: 負 (と NaN) なら NaN、0 なら -Infinity (残差が全く変わらない crate=0 のとき badcnv=1 → モード切替になる。
    以前は 0 も NaN にしていたので、不安定ラインで勾配が 0 になる点 (crate=0) で無限に反復した)"""
    if x != x:
        return float('nan')
    if x > 0.0:
        return math.log10(x)
    return float('-inf') if x == 0.0 else float('nan')


def _fmincube(f1, f2, g1, g2, g0, dg):
    """findr fmincube"""
    if f2 == 0.0:
        f = -0.5 * f1 * dg / ((g1 - g0) / f1 - dg)
    else:
        a = ((g1 - g0) / f1**2 - (g2 - g0) / f2**2) / (f1 - f2) + dg / f1 / f2
        b = (-f2 * (g1 - g0) / f1**2 + f1 * (g2 - g0) / f2**2) - dg * (f1 + f2) / f1 / f2
        s = math.sqrt(max(0.0, b**2 - 3.0 * a * dg))
        if b > 0.0:
            f = -dg / (s + b)
        else:
            f = (s - b) / 3.0 / a
    return min(M_EULER * f1, max(f1 / 16.0, f))


def _tdlfun(fmin, fmax, f, typ=0):
    """qdcell.f tdlfun: 目標 (fmin, fmax) と現在値 f の差。typ 1: log、2: atan"""
    if typ == 1:
        return _tdlfun(math.log(fmin), math.log(fmax), math.log(f))
    if typ == 2:
        return _tdlfun(math.atan(fmin), math.atan(fmax), math.atan(f))
    if fmin == fmax:
        return fmin - f
    if fmin < fmax:
        if f < fmin:
            return fmin - f
        if f > fmax:
            return fmax - f
        return 0.0
    fmid = (fmin + fmax) * 0.5
    if f > fmid and f < fmin:
        return fmin - f
    if f < fmid and f < fmax:
        return fmax - f
    return 0.0


def _tdfun1(fmin, fmax, v, kf, maxfit):
    """qdcell.f tdfun1 (ttrans=.false.)"""
    if kf in _LOGFIT:
        if maxfit:
            return math.log(fmax / v) if v > fmax else 0.0
        return _tdlfun(fmin, fmax, v, 1)
    if kf in _ATANFIT:
        if maxfit:
            vfa = abs(fmin)
            if v > vfa:
                return math.atan(vfa) - math.atan(v)
            if v < -vfa:
                return -math.atan(vfa) - math.atan(v)
            return 0.0
        return _tdlfun(fmin, fmax, v, 2)
    if kf in ('EX', 'EY'):
        if maxfit:
            vfa = abs(fmax)
            return vfa - v if v > vfa else max(-vfa - v, 0.0)
        return _tdlfun(fmin, fmax, v)
    if kf in ('NX', 'NY', 'NZ'):
        if maxfit:
            return min(0.0, fmin - v)
        return _tdlfun(fmin, fmax, v)
    if maxfit:
        vfa = abs(fmin)
        return vfa - v if v > vfa else max(-vfa - v, 0.0)
    return _tdlfun(fmin, fmax, v)


def _trdfun(kf, vb, ve):
    return ve / vb if kf in _LOGFIT else ve - vb


EMMINV = 1e-15     # toplvl.f の emminv


def _fmax(a, b):
    """gfortran の MAX(a,b): NaN は無視する (片方だけ NaN なら他方)"""
    if a != a:
        return b
    if b != b:
        return a
    return a if a >= b else b


def _fmin(a, b):
    """gfortran の MIN(a,b): NaN は無視する"""
    if a != a:
        return b
    if b != b:
        return a
    return a if a <= b else b


class Matcher:
    """SAD の FIT/FREE/FIX/GO。lattice は pysad.lattice.Lattice (CELL のリング)。

    emx, emy, dpmax は SAD と同じく先頭 MARK の EMITX/EMITY/DP (なければ既定 EMITX=EMITY=0, DP=0.01) を使う。
    fit_value(loc, fun, dp, vgoal, vnow) → None なら条件を外す、数値なら目標を置き換える (SAD の FitValue)。
    fit_weight(loc, fun, dp, w) → 重み (SAD の FitWeight)。"""

    def __init__(self, lattice, cell: bool = True, dp0: float = 0.0, dpmax: Optional[float] = None,
                 emx: Optional[float] = None, emy: Optional[float] = None, itmax: int = 40,
                 convergence: float = 1e-9, exponent: float = 2.0, offmw: float = 1.0,
                 out: Optional[Callable[[str], None]] = print):
        self.lat = lattice
        self.cell = cell
        self.dp0 = dp0
        m0 = lattice.elements[0]
        self.dpmax = dpmax if dpmax is not None else (max(0.0, float(getattr(m0, 'DP', 0.0) or 0.0)) or 0.01)
        self.emx = emx if emx is not None else float(getattr(m0, 'EMITX', 0.0) or 0.0)
        self.emy = emy if emy is not None else float(getattr(m0, 'EMITY', 0.0) or 0.0)
        self.itmax = itmax
        self.coumin = 0.01          # MINCOUP (tfsetparam)
        self.rsconv = convergence
        self.wexponent = max(0.1, min(4.0, exponent))
        self.offmw = offmw
        self.out = out or (lambda s: None)
        self.variables: List[dict] = []
        self.conditions: List[dict] = []
        self.fit_value: Optional[Callable] = None
        self.fit_weight: Optional[Callable] = None
        self.fit_function: Optional[Callable] = None    # SAD の FitFunction: () → 実数 or 実数リスト
        self.var_range: Optional[Callable] = None       # SAD の VariableRange[element, keyword, v_]: (family, key, v) → None / 実数 / (min, max)
        self.var_weight: Optional[Callable] = None      # SAD の VariableWeight[element, keyword, v_]: (family, key, gw) → None / 実数
        # FITP の fam 点 (tffamsetup: MatchingAmplitude {{dp, xamp, yamp},…} と NFAMP の有限振幅ファミリー)
        self.fam_setup: Optional[Callable] = None       # () → None / dict(amp=[[dp,x,y],…], nfp, em, tw=dict(ax bx ay by r1..r4 detr))
        self.optics_prolog: Optional[Callable] = None   # tffscalc.f:64 OpticsProlog (光学計算の前)、398 OpticsEpilog (後) の評価
        self.optics_epilog: Optional[Callable] = None
        self._res1: Dict[int, float] = {}               # 直近の tffscalc の idp ごとの残差 (FFS$SHOW 用)
        self._dps_show: Dict[int, float] = {}
        self.calc_optics_fixed: Optional[Callable] = None    # (delta, cod0) → rows: 与えた軌道のまわり (閉軌道を探さない)。None なら lat
        self._fam: Dict[int, dict] = {}                 # ix → dict(jfam=親の idp, kfam=番号 (x: 1..nfp, y: -1..-nfp), dfam=軌道のずれ 4)
        self._fam_def_done = False                      # GO (tffsa) ごとに 1 回だけ tffamsetup する
        self._cur_fam_id = 0                            # FitValue の id (kfam)
        self.dlim = False                               # VariableRange が値を返したら真 (tffsvlimit の dlim)
        self.nlat = len(lattice.elements) + 1
        self._mfpnt = self.nlat
        self._mfpnt1 = self.nlat
        self.nfp = 1
        self.residual = _Res()
        self.matching_residual = 0.0
        self.unstable_level = 0
        self._rows = None
        self.bipol = True
        self.absweit = True
        # SAD の不安定ラインの扱い (qcell1 / tffscalc): idp ごとの stabx/staby (optstat)、cellstab、qcod 等のメッセージ回数 nmes
        self.cellstab = True
        self.last_stab: Dict[int, tuple] = {}
        self.qstat_by: Dict[int, Optional[dict]] = {}
        self._last_qstat: Optional[dict] = None
        self._rows_idp: Dict[int, list] = {}   # idp ごとの最後に計算できた Twiss (SAD の utwiss の残り。overflow の点の残差はこれで評価される)
        self.nmes = 0                 # cellm の nmes (メッセージは nmes < nmmax=4 のときだけ。反復表を 1 行出すたびに 0 に戻る)
        self.lout_on = True           # qmdiag の lfno=0 (和共鳴のメッセージを 1 回出したら次の表示まで黙る)

    # ------------------------------------------------------------------ 要素の参照
    def _families(self):
        seen = {}
        for i, el in enumerate(self.lat.elements):
            seen.setdefault(el.name, []).append(i)
        return seen

    def locate(self, name: str) -> int:
        """ielmf: 要素名 → SAD の位置 (1 始まり、入口)。'NAME' はファミリーの最初、'NAME.n' は n 台目、$$$ は末尾"""
        if name in ('$$$', '***'):
            return self.nlat
        if name == '^^^':
            return 1
        fam, _, num = name.partition('.')
        idx = self._families().get(fam)
        if not idx:
            raise KeyError(f'?Undefined location {name}')
        n = int(num) if num else 1
        return idx[n - 1] + 1

    # ------------------------------------------------------------------ FREE / FIX
    def free(self, pattern: str, key: Optional[str] = None):
        """FREE pattern [keyword]: pattern は SAD のワイルドカード (* ?)。複数はカンマ区切り"""
        fams = self._families()
        for pat in pattern.split(','):
            pat = pat.strip()
            if not pat:
                continue
            hit = False
            for fam, idx in fams.items():
                if not fnmatch.fnmatchcase(fam, pat):
                    continue
                el0 = self.lat.elements[idx[0]]
                tn = type(el0).__name__
                k = (key or _DEFAULT_KEY.get(tn, '')).upper()
                if not k:
                    continue
                if tn in ('Marker', 'Monitor'):
                    continue
                hit = True
                if any(v['family'] == fam and v['key'] == k for v in self.variables):
                    continue
                els = [self.lat.elements[i] for i in idx]
                self.variables.append(dict(family=fam, key=k, els=els, type=tn, val=_sad_get(el0, k),
                                           vlim=[0.1, 1e10] if tn == 'Drift' and k == 'L' else [-1e10, 1e10]))
            if not hit:
                self.out(f'?Undefined element {pat}')
        self.variables.sort(key=lambda v: self.lat.elements.index(v['els'][0]))

    def fix(self, pattern: str = '*', key: Optional[str] = None):
        pats = [p.strip() for p in pattern.split(',') if p.strip()]
        self.variables = [v for v in self.variables
                          if not any(fnmatch.fnmatchcase(v['family'], p) and (key is None or v['key'] == key.upper())
                                     for p in pats)]

    # ------------------------------------------------------------------ FIT
    def fit_at(self, loc: Optional[str] = None, loc1: Optional[str] = None):
        """FIT [loc [loc1]]: 以後の条件の位置 (省略時は $$$ / 直前の位置を保持)"""
        if loc is None:
            self._mfpnt = self.nlat
            self._mfpnt1 = self.nlat
        else:
            self._mfpnt = self.locate(loc)
            self._mfpnt1 = self.locate(loc1) if loc1 else self._mfpnt

    def fit(self, fun: str, value, npoints: Optional[int] = None, loc: Optional[str] = None,
            loc1: Optional[str] = None):
        """条件を登録: fun は 'BX', 'NX', ... ('M' を付けると最大フィット)。value は数値または (min, max)。
        loc を与えると fit_at(loc[, loc1]) も行う (SAD: fit loc fun value)。"""
        if loc is not None:
            self.fit_at(loc, loc1)
        fun = fun.upper()
        maxfit = fun.endswith('M') and fun[:-1] in _FUNS
        kf = fun[:-1] if maxfit else fun
        if kf not in _FUNS:
            raise ValueError(f'unknown fit function {fun}')
        sc = PI2 if kf in ('NX', 'NY', 'NZ') else (math.pi / 180.0 if kf.startswith('CHI') else 1.0)
        if isinstance(value, (tuple, list)):
            fmin, fmax = float(value[0]) * sc, float(value[1]) * sc
        else:
            fmin = fmax = float(value) * sc
        if kf in ('BX', 'BY') and fmin <= 0.0:
            raise ValueError(f'Non-positive value for {kf}')
        kp, kp1 = min(self._mfpnt, self._mfpnt1), max(self._mfpnt, self._mfpnt1)
        for c in self.conditions:
            if c['kp'] == kp and c['kp1'] == kp1 and c['kf'] == kf:
                c.update(fmin=fmin, fmax=fmax)
                if npoints is not None:
                    c['mfitp'] = npoints + 1
                c['mfitp'] = -abs(c['mfitp']) if maxfit else abs(c['mfitp'])
                return
        mfitp = (npoints + 1) if npoints is not None else 2
        self.conditions.append(dict(id=len(self.conditions) + 1, kp=kp, kp1=kp1, kf=kf, fmin=fmin, fmax=fmax,
                                    mfitp=-abs(mfitp) if maxfit else abs(mfitp)))

    def fit_clear(self):
        """FIT; (条件を消す: SAD では mfitp=0 にする)"""
        self.conditions = []
        self._mfpnt = self._mfpnt1 = self.nlat

    def fitp(self, n: int):
        self.nfp = n
        for c in self.conditions:
            if abs(c['mfitp']) > 2:
                c['mfitp'] = int(math.copysign(abs(n) + 1, c['mfitp']))

    # ------------------------------------------------------------------ 光学の計算 (tffscalc の qcell1 部分)
    def _apply_vars(self):
        for v in self.variables:
            for el in v['els']:
                _sad_set(el, v['key'], v['val'])
        self.lat._twiss_cache = None

    def _msg(self, text: str, count: bool = False):
        """qcell1 のメッセージ (lfno>0 かつ nmes<nmmax=4 のとき)。count=True なら nmes を行末に出す (list-directed の i12)"""
        if not self.lout_on or self.nmes >= 4:
            return
        self.nmes += 1
        self.out(text + (f'{self.nmes:12d}' if count else ''))

    def _calc_optics(self, delta=None, cod0=None):
        """qcell1: Twiss を各要素の入口 (SAD の twiss(l)) で返す。失敗 (NaN) なら None。
        delta: 運動量偏差 (None なら dp0)。cod0: 閉軌道探索の初期値 (オフモメンタムの連鎖用)。
        |cos μ|>1 などの不安定は SAD と同じく続行して Twiss を返し、判定 (stabx, staby, over, codfnd) を self._last_qstat に残す"""
        self._last_qstat = None
        flags = getattr(self, 'resflags', (True, True, True, False))
        try:
            rows = self.lat.twiss_table(delta=self.dp0 if delta is None else delta, cod=cod0,
                                        sad_stab={'intres': bool(flags[0]), 'halfres': bool(flags[1])})
        except Exception:
            self._msg(' *****qcod---> Overflow & closed orbit not found ', True)
            return None
        q = dict(self.lat._qstat) if self.lat._qstat is not None else None
        self._last_qstat = q
        if q is not None:
            if not q['codfnd'] and not q['over']:
                self._msg(' *****qcod---> Closed orbit not found ', True)
            for m in q['msgs']:
                if self.lout_on and self.nmes < 4:
                    self.nmes += 1
                    self.out(m)
                    if 'Sum resonance' in m:
                        self.lout_on = False
            if q['over']:
                if self.lout_on and self.nmes < 4:
                    self.nmes += 1
                    self.out(f'***qtwiss---> Overflow (fb) {1:8d}{self.nlat:8d}{fmt_g(0.0, 15, 7)}')
        if any(not np.isfinite(r.betx) or not np.isfinite(r.bety) for r in (rows[0], rows[-1])):
            return None
        return rows

    geo_fn = None                     # 幾何 (GeoResult) を現在の要素から作る関数 (ffscmd が ORG 込みで差し替える)
    _geo = None                       # 直近の _tffscalc の幾何 (条件に GX..CHI3 があるときだけ作る)

    def _compute_geo(self):
        if self.geo_fn is not None:
            return self.geo_fn()
        from . import geo as _geo
        return _geo.compute_geo(self.lat.elements)

    def _tgfun(self, rows, kf, kp):
        """tgfun: 条件 kf の位置 kp の値。GX/GY/GZ は geo(:,4,kp)、CHI1-3 は tfchi(geo(:,:,kp)) (tgfun.f)"""
        if kf in _GEO_KF:
            g = self._geo
            if g is None:
                g = self._geo = self._compute_geo()
            if kf[0] == 'G':
                return float(g.gxyz(kp)['XYZ'.index(kf[1])])
            return float(g.chi(kp)[int(kf[3]) - 1])
        r = rows[kp - 1]
        return float(getattr(r, _ROW_ATTR[kf]))

    def _pos(self, rows, kp):
        return float(rows[kp - 1].s)

    # ------------------------------------------------------------------ 残差 (tffscalc + tdfun + twfit)
    def _nfr(self):
        """オフモメンタム点の数 (tffsa: nfr = max((|mfitp|-1)/2))"""
        nfr = 0
        for c in self.conditions:
            if c['mfitp'] != 0:
                nfr = max(nfr, (abs(c['mfitp']) - 1) // 2)
        return nfr

    def _dp_list(self, nfr):
        """dp(idp) = dpmax*idp/nfr (DPM 指定時は ±1 を DPM に)"""
        dps = {0: 0.0}
        for i in range(1, nfr + 1):
            dps[i] = self.dpmax * i / max(nfr, 1)
            dps[-i] = -self.dpmax * i / max(nfr, 1)
        dpm2 = getattr(self, 'dpm2', 0.0)
        if dpm2 > 0.0 and nfr >= 1:
            dps[1], dps[-1] = dpm2, -dpm2
        return dps

    def _optics_all(self, nfr):
        """tffscalc: idp=0 (δ=dp0) を計算し、±1, ±2, … を前の点の閉軌道を初期値に連鎖で計算する。
        戻り: (rows_by_idp, stab_by_idp, over_by_idp)。SAD の INTRES/HALFRES/SUMRES/DIFFRES の安定判定も行う。
        その後 fam 点 (MatchingAmplitude) を親の点の軌道 + 振幅で計算する (tffscalc の fam=.true. のループ)"""
        rows_by, stab, over = self._optics_idp(nfr)
        self._fam_points(nfr, rows_by, stab, over)
        return rows_by, stab, over

    # ---- fam 点 (tffamsetup.f / tffscalc.f)
    def _fam_define(self, nfr):
        """tffamsetup: MatchingAmplitude = {{dp, xamp, yamp}, …} と NFAMP から、FITP の各点 k のまわりの有限振幅の軌道を作る。
        先頭 MARK の Twiss (twmov(1,·)) で正規モード振幅 (x0 = sqrt(βx em), px0 = x0/βx, em = |EMITX|+|EMITY|) の NFAMP 本の位相の軌道
        dx(i) = (cos θ x0, sin θ px0 − α/β·x, …) を tfnormaltophysical (R1-R4) で物理座標に移し、xamp/yamp 倍する。
        添字 ix は nfr+1, −(nfr+1), nfr+2, … と交互 (nfa=merge(−nfa, −nfa+1, nfa≥0))。dp は親と同じ"""
        self._fam = {}
        if self.fam_setup is None or nfr <= 0:
            return
        d = self.fam_setup()
        if not d:
            return
        amp = d.get('amp')
        if not isinstance(amp, list) or not amp:
            return
        nfp = d.get('nfp', 0)
        if not isinstance(nfp, (int, float)) or nfp <= 0:
            self.out(f'Non-positive NFAMP - FAM is skipped:{nfp}')
            return
        nfp = min(12, int(nfp))
        tw = d['tw']
        em = float(d['em'])
        bx, by, ax, ay = tw['bx'], tw['by'], tw['ax'], tw['ay']
        x0 = math.sqrt(max(0.0, bx * em))
        px0 = x0 / bx
        y0 = math.sqrt(max(0.0, by * em))
        py0 = y0 / by
        T = _normal_to_physical(tw)
        dx = np.zeros((4, nfp * 2))
        for i in range(nfp):
            c = math.cos(PI2 / nfp * i)
            sn = math.sin(PI2 / nfp * i)
            dx[0, i] = c * x0
            dx[1, i] = sn * px0 - ax / bx * dx[0, i]
            dx[2, i + nfp] = c * y0
            dx[3, i + nfp] = sn * py0 - ay / by * dx[2, i + nfp]
        dxp = T @ dx
        dps = self._dp_list(nfr)
        dpm = {}
        for i in range(-nfr, nfr):
            dpm[i + 1] = (dps[i] + dps[i + 1]) * 0.5
        dpm[-nfr] = dps[-nfr]
        dpm[nfr + 1] = dps[nfr]
        nfa = nfr + 1
        fam = {}
        for e in amp:
            try:
                dpi, x, y = (float(v) for v in e)
            except (TypeError, ValueError):
                self.out('MatchingAmplitude must be {{dp,xamp,yamp},...}.')
                break
            for k in range(-nfr, nfr + 1):
                if dpm[k] <= dpi < dpm[k + 1]:
                    if x != 0.0:
                        for j in range(nfp):
                            fam[nfa] = dict(jfam=k, kfam=j + 1, dfam=dxp[:, j] * x)
                            nfa = -nfa if nfa >= 0 else -nfa + 1
                    if y != 0.0:
                        for j in range(nfp, nfp * 2):
                            fam[nfa] = dict(jfam=k, kfam=nfp - (j + 1), dfam=dxp[:, j] * y)
                            nfa = -nfa if nfa >= 0 else -nfa + 1
                    break
        self._fam = fam

    def _fam_dps(self, dps):
        """fam 点の dp (親と同じ) を dps に足した辞書"""
        if not self._fam:
            return dps
        d = dict(dps)
        for ix, f in self._fam.items():
            d[ix] = dps[f['jfam']]
        return d

    def _eval_idps(self, mp):
        """tdfun: 条件を評価する点。通常の点 −mp..mp と、親 jfam がその範囲にある fam 点"""
        out = list(range(-mp, mp + 1))
        out.extend(ix for ix, f in self._fam.items() if -mp <= f['jfam'] <= mp)
        return out

    def _calc_optics_fixed_default(self, delta, cod0):
        """fam 点の光学: 閉軌道を探さず cod0 のまわりの一周行列で周期解 (qcell1 の fam=.true.: INTRES/HALFRES で不安定にしない)"""
        self._last_qstat = None
        try:
            rows = self.lat.twiss_table(delta=delta, sad_stab={'intres': False, 'halfres': False}, fixed_cod=cod0)
        except Exception:
            return None
        self._last_qstat = dict(self.lat._qstat) if self.lat._qstat is not None else None
        if any(not np.isfinite(r.betx) or not np.isfinite(r.bety) for r in (rows[0], rows[-1])):
            return None
        return rows

    def _fam_points(self, nfr, rows_by, stab, over):
        """tffscalc の fam ループ: 親 i2=jfam(ix) が over なら utwiss をコピー (over)、そうでなければ親の入口の軌道 + dfam を初期軌道に
        qcell1(fam=.true.)。fam 点の stabx/staby は常に真 (optstat(ix)%stabx = fam .or. …)"""
        if not self._fam_def_done:
            self._fam_define(nfr)
            self._fam_def_done = True
        if not self._fam:
            return
        dps = self._dp_list(nfr)
        calc = self.calc_optics_fixed or self._calc_optics_fixed_default
        for ix, f in self._fam.items():
            k = f['jfam']
            prow = rows_by.get(k)
            if prow is None or over.get(k):
                rows_by[ix] = self._rows_idp.get(ix)
                over[ix] = True
                stab[ix] = (True, True)
                self.qstat_by[ix] = None
                continue
            r0 = prow[0]
            dfam = f['dfam']
            c0 = np.array([r0.x_cod + dfam[0], r0.px_cod + dfam[1], r0.y_cod + dfam[2], r0.py_cod + dfam[3],
                           r0.z_cod, self.dp0 + dps[k]])
            rk = calc(self.dp0 + dps[k], c0)
            qk = self._last_qstat
            self.qstat_by[ix] = qk
            rows_by[ix] = rk if rk is not None else self._rows_idp.get(ix)
            over[ix] = rk is None or bool(qk and qk['over'])
            if rk is not None:
                self._rows_idp[ix] = rk
            stab[ix] = (True, True)

    def _optics_idp(self, nfr):
        """tffscalc の idp=0, ±1, ±2, … の部分 (_optics_all の本体)"""
        dps = self._dp_list(nfr)
        rows_by = {}
        stab = {}
        over = {}
        self.qstat_by = {}
        r0 = self._calc_optics(self.dp0, None)
        q0 = self._last_qstat
        self.qstat_by[0] = q0
        rows_by[0] = r0
        over[0] = r0 is None or bool(q0 and q0['over'])
        if r0 is not None:
            self._rows_idp[0] = r0
        if r0 is None:
            stab[0] = (False, False)
            for i in dps:
                rows_by[i] = None if i == 0 else self._rows_idp.get(i)
                over[i] = True
                stab[i] = (False, False)
                self.qstat_by[i] = None
            return rows_by, stab, over
        # idp=0 の判定は qcell1 の結果そのもの (|cos μ|>1 での INTRES/HALFRES、tracex、codfnd、qmdiag の stab)
        stab[0] = (bool(q0['stabx']), bool(q0['staby'])) if (q0 is not None and self.cell) else (True, True)
        pi = math.pi
        nux0, nuy0 = r0[-1].mux, r0[-1].muy
        ref = (math.floor(nux0 / PI2), math.floor(nuy0 / PI2), math.floor(nux0 / pi), math.floor(nuy0 / pi),
               math.floor((nux0 + nuy0) / PI2), math.floor((nux0 - nuy0) / PI2))
        flags = getattr(self, 'resflags', (True, True, True, False))
        for sign in (1, -1):
            prev = r0
            prev_over = over[0]
            for k in range(1, nfr + 1):
                idp = sign * k
                if prev is None or prev_over:
                    # tffscalc: 前の点が over なら utwiss をコピーして over (optstat は両面不安定)。条件はこの idp の前回の値で評価される
                    rows_by[idp] = self._rows_idp.get(idp)
                    over[idp] = True
                    stab[idp] = (False, False)
                    self.qstat_by[idp] = None
                    continue
                c0 = np.array([prev[0].x_cod, prev[0].px_cod, prev[0].y_cod, prev[0].py_cod, 0.0, self.dp0 + dps[idp]])
                rk = self._calc_optics(self.dp0 + dps[idp], c0)
                qk = self._last_qstat
                self.qstat_by[idp] = qk
                rows_by[idp] = rk if rk is not None else self._rows_idp.get(idp)
                over[idp] = rk is None or bool(qk and qk['over'])
                if rk is None:
                    stab[idp] = (False, False)
                    prev = None
                    continue
                self._rows_idp[idp] = rk
                prev_over = over[idp]
                nux, nuy = rk[-1].mux, rk[-1].muy
                cur = (math.floor(nux / PI2), math.floor(nuy / PI2), math.floor(nux / pi), math.floor(nuy / pi),
                       math.floor((nux + nuy) / PI2), math.floor((nux - nuy) / PI2))
                intres, halfres, sumres, diffres = flags
                sx = ((not intres or cur[0] == ref[0]) and (not halfres or cur[2] == ref[2])
                      and (not sumres or cur[4] == ref[4]) and (not diffres or cur[5] == ref[5]))
                sy = ((not intres or cur[1] == ref[1]) and (not halfres or cur[3] == ref[3])
                      and (not sumres or cur[4] == ref[4]) and (not diffres or cur[5] == ref[5]))
                if qk is not None:
                    sx = sx and bool(qk['stabx'])
                    sy = sy and bool(qk['staby'])
                stab[idp] = (bool(sx), bool(sy)) if self.cell else (True, True)
                prev = rk
        return rows_by, stab, over

    def _tffscalc(self, wcal: bool, lout: bool = True):
        """tffscalc + tdfun + twfit。戻り: (dfs: dict (id, idp)→df·wiq, r: _Res, rows (idp=0))
        wcal: 重み wfit を計算し直す"""
        self._apply_vars()
        if self.optics_prolog is not None:
            self.optics_prolog()          # tffscalc.f:64: OpticsProlog を評価してから光学計算
        self._geo = None                  # 変数を反映した要素で幾何を作り直す (必要なときに _tgfun が作る)
        nfr = self._nfr()
        rows_by, stab, over = self._optics_all(nfr)
        rows = rows_by[0]
        if rows is None:
            # qcod のオーバーフロー: SAD の qcell1 は twiss を書かずに戻る → twiss 配列には前回の値が残り、残差・条件はその値で評価される
            # (先頭行だけは tffscalc の twmov で先頭 MARK の値 = 最後に安定だった周期解に差し替わる)
            rows = self._stale_rows()
            rows_by[0] = rows
        else:
            self._rows = rows
        nlat = self.nlat
        if wcal:
            maxf = nlat if self.cell else max([nlat] + [c['kp1'] for c in self.conditions])
            self.etamax = max([1e-2] + [max(abs(rows[l - 1].etx), abs(rows[l - 1].ety)) for l in range(1, maxf + 1)])
            self.avebeta = (self._pos(rows, maxf) - self._pos(rows, 1)) / max(rows[maxf - 1].mux, rows[maxf - 1].muy)
            # tffscalc.f 442-452: emminv=1e-15, coumin=MINCOUP (既定 0.01), coum0=1e-3。
            # gfortran の MAX/MIN は NaN を無視する (EMIT 失敗で EMITX=EMITY=NaN のとき emminv 由来の値になる)
            emx, emy = self.emx, self.emy
            em = _fmax(EMMINV, abs(emx) + abs(emy))
            rr = abs(emx / _fmax(EMMINV, emy)) + abs(emy / _fmax(EMMINV, emx))
            inv = (1.0 / rr) if rr != 0.0 else math.inf          # Fortran は 1/0 = Inf
            coum = _fmin(1.0, _fmax(_fmax(self.coumin, 1e-3), inv))
            coum = coum / (1.0 + coum)
            emxx = _fmax(emx, coum * em)
            emyy = _fmax(emy, coum * em)
            dpm = max(0.001, self.dpmax)
            coup = math.sqrt(abs(emxx / emyy) + abs(emyy / emxx))
            for c in self.conditions:
                k, j = c['kf'], c['kp']
                bx, ax, by, ay = rows[j - 1].betx, rows[j - 1].alfx, rows[j - 1].bety, rows[j - 1].alfy
                if k in ('EX', 'PEX'):
                    w = dpm / math.sqrt(emxx * bx)
                elif k in ('EPX', 'PEPX'):
                    w = dpm * math.sqrt(bx / (1.0 + ax**2) / emxx)
                elif k in ('EY', 'PEY'):
                    w = dpm / math.sqrt(emyy * by)
                elif k in ('EPY', 'PEPY'):
                    w = dpm * math.sqrt(by / (1.0 + ay**2) / emyy)
                elif k == 'R1':
                    w = coup * math.sqrt(bx / by)
                elif k == 'R2':
                    w = coup / math.sqrt(bx * by)
                elif k == 'R3':
                    w = coup * math.sqrt(bx * by)
                elif k == 'R4':
                    w = coup * math.sqrt(by / bx)
                elif k == 'DX':
                    w = 1.0 / math.sqrt(bx * emxx)
                elif k == 'DPX':
                    w = math.sqrt(bx / (1.0 + ax**2) / emxx)
                elif k == 'DY':
                    w = 1.0 / math.sqrt(by * emyy)
                elif k == 'DPY':
                    w = math.sqrt(by / (1.0 + ay**2) / emyy)
                elif k == 'LENG':
                    w = 1.0 / (self._pos(rows, maxf) - self._pos(rows, 1)) * max(rows[maxf - 1].mux, rows[maxf - 1].muy)
                elif k == 'DZ':
                    w = 0.01 * math.sqrt(max(rows[maxf - 1].mux, rows[maxf - 1].muy) / em
                                         / max(1.0, self._pos(rows, maxf) - self._pos(rows, 1)))
                elif k in ('BMAGX', 'BMAGY', 'BMAGZ'):
                    w = 3.0
                elif k in ('GX', 'GY', 'GZ'):
                    w = 1e-3 * math.sqrt(max(rows[maxf - 1].mux, rows[maxf - 1].muy) / em
                                         / max(1.0, self._pos(rows, maxf) - self._pos(rows, 1)))
                elif k.startswith('CHI'):
                    w = 1e3
                else:
                    w = 1.0
                c['wfit'] = w
        # tdfun: 条件ごと・idp ごとの残差 (fam 点 (|idp| > nfr) も親 jfam が範囲内なら評価する)
        dps = self._fam_dps(self._dp_list(nfr))
        xix = getattr(self, 'xix', 0.0)
        xiy = getattr(self, 'xiy', 0.0)
        dfs = {}
        for c in self.conditions:
            if c['mfitp'] == 0:
                continue
            kf, kp, kp1 = c['kf'], c['kp'], c['kp1']
            maxfit = c['mfitp'] < 0
            mp = (abs(c['mfitp']) - 1) // 2
            name = self._elname(kp)
            for idp in self._eval_idps(mp):
                rws = rows_by.get(idp)
                if rws is None:
                    continue                      # overflow: 条件なし (SAD は前の family の値を使う)
                fmin, fmax = c['fmin'], c['fmax']
                if idp == 0 and mp * 2 + 1 == abs(c['mfitp']):
                    continue                      # tdfun: |mfitp| が奇数 (FITP 偶数) の条件は中心 (idp=0) では評価しない
                if kf in ('DX', 'DPX', 'DY', 'DPY') and abs(idp) > nfr:
                    continue                      # tdfun: 軌道の条件は fam 点 (|idp| > nfr) では評価しない
                if idp != 0 and kf in ('LENG', 'GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3', 'TRZ'):
                    continue                      # tdfun: kf > mfittry (TRX, TRY は含まない) は idp≠0 では条件にしない
                self._cur_fam_id = self._fam[idp]['kfam'] if idp in self._fam else 0     # FitValue の id (kfam)
                if (idp != 0 or mp * 2 + 1 != abs(c['mfitp'])) and kp == self.nlat and -1 <= idp <= 1:
                    if kf == 'NX':
                        fmin, fmax = fmin + dps[idp] * xix, fmax + dps[idp] * xix
                    elif kf == 'NY':
                        fmin, fmax = fmin + dps[idp] * xiy, fmax + dps[idp] * xiy
                dp = self.dp0 + dps[idp]
                if kp != kp1:
                    vb = self._tgfun(rws, kf, kp)
                    ve = self._tgfun(rws, kf, kp1)
                    if self.fit_value is not None:
                        g = self.fit_value(name, kf, dp, fmin if fmin == fmax else (fmin, fmax), ve, self._elname(kp1), vb)
                        if g is None:
                            continue
                        if isinstance(g, (int, float)):
                            fmin = fmax = fmin + _trdfun(kf, vb, float(g))
                        elif isinstance(g, tuple):
                            fmin, fmax = g
                    df = _tdfun1(fmin, fmax, _trdfun(kf, vb, ve), kf, maxfit)
                    if c.get('cellcond') and abs(df) < 1e-8:
                        continue                  # tdfun: cell かつ tfsetupcell の条件で |df| < abmax なら捨てる
                else:
                    if kf in ('TRX', 'TRY'):
                        # tgfun: optstat(idp)%tracex / tracey (qcell1 の末尾の「トレース」。不安定なら 4/trace に置き換え済み)
                        qq = self.qstat_by.get(idp)
                        v = float(qq['tracex' if kf == 'TRX' else 'tracey']) if qq else 0.0
                    else:
                        v = self._tgfun(rws, kf, kp)
                    if self.fit_value is not None:
                        g = self.fit_value(name, kf, dp, fmin if fmin == fmax else (fmin, fmax), v)
                        if g is None:
                            continue
                        if isinstance(g, (int, float)):
                            fmin = fmax = float(g)
                        elif isinstance(g, tuple):
                            fmin, fmax = g
                    df = _tdfun1(fmin, fmax, v, kf, maxfit)
                    if maxfit and df == 0.0:
                        continue
                wiq = c['wfit']
                if self.fit_weight is not None:
                    w2 = self.fit_weight(name, kf, dp, wiq)
                    if isinstance(w2, (int, float)):
                        wiq = float(w2)
                dfs[(c['id'], idp)] = (df, wiq)
        # FitFunction (tffsfitfun): 実数か実数リストを条件として追加 (重み 1、idp=0)
        if self.fit_function is not None:
            ffv = self.fit_function()
            if ffv is not None:
                vals = ffv if isinstance(ffv, (list, tuple)) else [ffv]
                for k, x in enumerate(vals):
                    dfs[(-(k + 1), 0)] = (float(x), 1.0)
        # 重みと残差 (tffscalc)
        rw = 0.0
        wsum = 0.0
        res1 = {}
        out = {}
        for (cid, idp), (df, wiq) in dfs.items():
            if idp != 0:
                wi = (self.offmw / 2.0 / math.sqrt(max(1, abs(self._cond(cid)['mfitp']))))**2
                wsum += wi
                wiq = wiq * wi**(1.0 / self.wexponent)
            else:
                wsum += 1.0
            dfw = df * wiq
            out[(cid, idp)] = dfw
            if dfw != 0.0:
                drw = min(1e50, max(1e-50, abs(dfw)))**self.wexponent
                rw += drw
                res1[idp] = res1.get(idp, 0.0) + drw
        r = _Res(wsum * (max(rw, 1e-50) / wsum)**(2.0 / self.wexponent) if rw > 0.0 and wsum > 0 else 0.0, 0)
        self.last_stab = dict(stab)
        if self.cell:
            # tffscalc.f: 残差が 0 でない idp について stabx/staby の偽を 1 ずつ数える。どの idp も残差 0 (条件なし) なら
            # nstab=0 で、cellstab だけ全 idp の stabx・staby の積にする
            cond_idps = self._cond_idps(nfr)
            cellstab = True
            zerores = True
            for idp in sorted(stab.keys()):
                has = res1.get(idp, 0.0) != 0.0
                if not has and over.get(idp) and rows_by.get(idp) is None and idp in cond_idps:
                    # オーバーフローした点は条件を評価できない (pysad) が、SAD は前の Twiss で評価した残差が非 0 → 両面 unstable を数える
                    has = True
                    if idp != 0:
                        r.nstab += 2
                        zerores = False
                        cellstab = False
                        continue
                if has:
                    zerores = False
                    sx, sy = stab.get(idp, (True, True))
                    if not sx:
                        r.nstab += 1
                        cellstab = False
                    if not sy:
                        r.nstab += 1
                        cellstab = False
            if over.get(0) and rows_by.get(0) is None and rw > 0.0:
                r.nstab += 2
            if zerores:
                cellstab = all(sx and sy for sx, sy in stab.values())
            self.cellstab = cellstab
        else:
            self.cellstab = True
            for idp, (sx, sy) in stab.items():
                if not sx:
                    r.nstab += 2
        self.wsum = wsum
        self._res1 = dict(res1)
        self._dps_show = dict(dps)
        if self.optics_epilog is not None:
            self.optics_epilog()          # tffscalc.f:398: 残差を出した後に OpticsEpilog
        return out, r, rows

    def _stale_rows(self):
        """SAD の twiss(·,0,·) に残っている前回の値 (qcell1 が overflow で戻ったとき)。先頭行は MARK の値に差し替える"""
        prev = self._rows
        if prev is None:
            from .lattice import nan_rows
            return nan_rows(self.lat)          # SAD は壊れた値で続行する (pysad は NaN の表)
        rows = list(prev)
        init = getattr(self.lat, '_mark_init', None)
        if init:
            r0 = TwissRow(**{**prev[0].__dict__})
            r0.alfx, r0.betx, r0.alfy, r0.bety = init['ax'], init['bx'], init['ay'], init['by']
            r0.etx, r0.etpx, r0.ety, r0.etpy = init['ex'], init['epx'], init['ey'], init['epy']
            r0.r1, r0.r2, r0.r3, r0.r4, r0.detr = init['r1'], init['r2'], init['r3'], init['r4'], init['detr']
            r0.mux = r0.muy = 0.0
            rows[0] = r0
        return rows

    def _cond_idps(self, nfr):
        """条件 (tdfun が対象にする) が 1 つでもある idp の集合"""
        out = set()
        for c in self.conditions:
            if c['mfitp'] == 0:
                continue
            mp = (abs(c['mfitp']) - 1) // 2
            for idp in range(-mp, mp + 1):
                if idp == 0 and mp * 2 + 1 == abs(c['mfitp']):
                    continue
                if idp != 0 and c['kf'] in ('LENG', 'GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3', 'TRZ'):
                    continue
                out.add(idp)
        if self.fit_function is not None:
            out.add(0)
        return out

    def _elname(self, kp):
        if kp >= self.nlat:
            return '$$$'
        return self.lat.elements[kp - 1].name

    # ------------------------------------------------------------------ 重み・制限
    def _tweigh(self, v, val0):
        """tffsmatch.f tweigh (absweit=.true.)"""
        key, tn = v['key'], v['type']
        gw = 1.0
        if key == 'K1':
            gw = 1.0
        elif key == 'L':
            gw = 1.0 / self.avebeta**2
        elif key in ('ANGLE', 'K0', 'SK0'):
            gw = 1.0 / self.avebeta
        elif key == 'K2':
            gw = max(1e-3, self.dpmax) * max(1e-2, self.etamax)
        elif key == 'K3':
            gw = (max(1e-3, self.dpmax) * max(1e-2, self.etamax))**2
        elif key == 'K4':
            gw = (max(1e-3, self.dpmax) * max(1e-2, self.etamax))**3
        elif tn == 'Multipole' and key[0] == 'K' and key[1:].isdigit() and int(key[1:]) >= 1:
            gw = (max(1e-3, self.dpmax) * max(1e-2, self.etamax))**int(key[1:])
        elif key in ('VOLT', 'FREQ'):
            gw = 100.0 / max(1.0, abs(val0))
        elif key in ('DX', 'DY'):
            if tn in ('Quadrupole', 'Bend', 'Multipole'):
                gw = 0.1 / self.avebeta
        elif key == 'BZ':
            gw = 1.0 / self.lat_brho
        if not self.absweit:
            vmin = {'K1': 1e-5, 'K2': 1e-3, 'K3': 1e-1, 'K4': 1e1}.get(key, 1e-6)
            gw = math.sqrt(max(vmin, abs(val0)) / gw)
        if self.var_weight is not None:
            # tffsmatch.f:924 tffsvarfun(2,..): VariableWeight[element, keyword, gw] が実数を返せばそれを重みにする
            w = self.var_weight(v['family'], key, gw)
            if w is not None:
                gw = w
        return gw

    def _vlimit(self, v, val, val0):
        """tffsvlimit: (vl, limited)"""
        vl1, vl2 = v['vlim']
        vl = val
        limited = False
        if vl < vl1:
            vl, limited = vl1, True
        elif vl > vl2:
            vl, limited = vl2, True
        vl0 = max(vl1, min(vl2, val0))
        if not self.bipol and vl * vl0 < 0.0:
            vl, limited = 0.0, True
        if self.var_range is not None:
            # tffsvlimit の VariableRange: 実数なら 0 (偽) で現在値に戻す (真なら dlim だけ)、{min,max} なら範囲を狭める
            r = self.var_range(v['family'], v['key'], val)
            if isinstance(r, float):
                self.dlim = True
                if r == 0.0:
                    vl, limited = vl0, True
            elif isinstance(r, tuple):
                self.dlim = True
                vl1 = max(vl1, r[0])
                vl2 = min(vl2, r[1])
                if vl < vl1:
                    vl, limited = vl1, True
                elif vl > vl2:
                    vl, limited = vl2, True
        return vl, vl1, vl2, limited

    # ------------------------------------------------------------------ ソルバ
    def _tfsolv(self, qu, dfv, wlimit, ids, fact):
        """tfsolv: 最大フィット条件は最初は外し、解が条件を越えるなら入れて解き直す。戻り (dval, dg)"""
        nq, nv = qu.shape
        fitsel = np.array([self._cond(i)['mfitp'] > 0 for i in ids])
        if not fitsel.any():
            fitsel[:] = True
        eps = TSTOL / fact
        nagain = 0
        while True:
            quw = qu[fitsel] * wlimit[None, :]
            b = dfv[fitsel]
            dval = self._tsolva(quw, b, eps)
            s = qu @ (wlimit * dval)
            dg = 0.0
            again = False
            for i in range(nq):
                if dfv[i] != 0.0:
                    if self.wexponent == 2.0:
                        dg -= dfv[i] * s[i]
                    else:
                        dg -= abs(dfv[i])**self.wexponent / dfv[i] * s[i]
                if not fitsel[i] and (s[i] - dfv[i]) * dfv[i] < 0.0:
                    fitsel[i] = True
                    again = True
            if again:
                nagain += 1
                if nagain <= nq:
                    continue
            break
        return dval, self.wexponent * dg

    @staticmethod
    def _tsolva(a, b, eps):
        """tsolva → tsvd (src/tsolvg.f)。SAD の 2 重対角化 + quick path / 暗黙 QR の直訳 (pysad/tsvd.py)。
        numpy の硬い打ち切りの SVD とは、条件数が 1e6 以上の J (INS のフィットなど) で解が大きく違う。
        環境変数 PYSAD_TSVD_NUMPY=1 で従来の numpy 版 (特異値 ≤ eps·最大 を落とす最小ノルム解) に戻る"""
        if a.size == 0:
            return np.zeros(a.shape[1])
        if not _os.environ.get('PYSAD_TSVD_NUMPY'):
            from .tsvd import tsvd
            return np.nan_to_num(tsvd(a, b, eps))
        u, s, vt = np.linalg.svd(a, full_matrices=False)
        smax = s[0] if len(s) else 0.0
        inv = np.array([1.0 / x if x > eps * smax and x > 0 else 0.0 for x in s])
        x = vt.T @ (inv * (u.T @ b))
        return np.nan_to_num(x)

    def _tgrad(self, qu, dfv, wlimit):
        quw = qu * wlimit[None, :]
        r = 0.0
        dfw = np.zeros_like(dfv)
        for i, d in enumerate(dfv):
            if d != 0.0:
                dfwi = abs(d)**self.wexponent
                r += dfwi
                dfw[i] = dfwi / d
        grad = quw.T @ dfw
        sg = float(np.sum(grad**2))
        if sg != 0.0:
            grad = grad * r / sg
        return grad

    def _cond(self, cid):
        if isinstance(cid, tuple):
            cid = cid[0]
        if cid < 0:
            return {'mfitp': 2, 'kf': 'FITFUNCTION', 'kp': self.nlat, 'kp1': self.nlat}
        return self.conditions[cid - 1]

    # ------------------------------------------------------------------ GO (tffsmatch)
    def go(self, calc_only: bool = False) -> _Res:
        """GO: マッチング。calc_only=True は CALC (残差の評価だけ)。戻りは残差 _Res。
        CELL では tffsa.f の tfsetupcell と同じく、条件の最後に「周期性」の 2 条件 AX(1→$$$)=0, AY(1→$$$)=0 を加える
        (点数は既存条件の最大 |mfitp|)。周期解が安定なら AX, AY が 1 周で閉じて |df|<1e-8 で捨てられる (tdfun: cell かつ ka が
        この 2 条件で |df|<abmax なら cycle) ので効かないが、不安定で cosμ→1/cosμ に置き換えた「周期解」では 1 周で閉じず、
        残差と微分に効く (SAD の不安定ラインでの反復の経路を決める)"""
        n0 = len(self.conditions)
        if self.cell:
            mfc = max([1] + [abs(c['mfitp']) for c in self.conditions])
            for kf in ('AX', 'AY'):
                self.conditions.append(dict(id=len(self.conditions) + 1, kp=1, kp1=self.nlat, kf=kf, fmin=0.0, fmax=0.0,
                                            mfitp=mfc, cellcond=True))
        try:
            return self._go_impl(calc_only)
        finally:
            del self.conditions[n0:]

    def _go_impl(self, calc_only: bool = False) -> _Res:
        self.lat_brho = getattr(self.lat, 'brho', 1.0)
        self._fam_def_done = False        # tffsa: GO/CALC のたびに tffamsetup (MatchingAmplitude) をやり直す
        nvar = len(self.variables)
        fitflg = not calc_only
        self.nmes = 0
        self.lout_on = True
        rl = abs(self.rsconv) * max(len(self.conditions), 1)
        out = self.out
        bestval = np.array([v['val'] for v in self.variables])
        sexp = '  CALEXP'
        wcal = True
        chgmod = True
        iter_ = 0
        newton = True
        imprv = True
        fact = 1.0
        rbest = _Res(1e300, 10000) if fitflg else _Res()
        rprint = rbest
        ra = 0.0
        aimprv = 0.0
        crate = 1.0
        outt = True
        nretry = 1 if self.cell else 0
        aitm1 = self.itmax * ALIT
        aitm2 = self.itmax
        nvara = nvar
        wlimit = np.ones(nvar)
        wvar = np.ones(nvar)
        f1 = f2 = g1 = g2 = 0.0
        rp0 = 0.0
        qu = None
        df0 = None
        ids0 = None
        r = _Res()
        self.no_variable = False
        if fitflg and nvar <= 0:
            out('?No variable.')
            self.no_variable = True
            fitflg = False
        # 変数の初期制限 (tffssetlimit)
        for i, v in enumerate(self.variables):
            vl, _, _, lim = self._vlimit(v, v['val'], v['val'])
            if lim:
                v['val'] = vl
        bestval = np.array([v['val'] for v in self.variables])
        while True:                                  # do9000
            exit9000 = False
            while True:                              # do200
                if not fitflg:
                    for i, v in enumerate(self.variables):
                        v['val'] = bestval[i]
                beforecal = np.array([v['val'] for v in self.variables])
                # tffsmatch.f: wcal=chgmod。モード切替・最初の反復の直後は重み wfit・etamax・avebeta を現在の光学から計算し直す
                wcal = chgmod
                dfs, r, rows = self._tffscalc(wcal)
                if _DEBUG:
                    print(f'[match6] iter={iter_} r={r.r:.16e} nstab={r.nstab} rbest={rbest.r:.16e} fact={fact:.4g} newton={newton} chgmod={chgmod} fitflg={fitflg} vals={[round(v["val"], 8) for v in self.variables]}')
                ids = list(dfs.keys())
                dfv = np.array([dfs[i] for i in ids])
                convgo = r.res2r() <= rl
                if convgo:
                    out(f' Matched. ({fmt_g(r.res2r())}) DP ={self.dpmax:8.5f}  DP0 ={self.dp0:8.5f}'
                        f'  ExponentOfResidual ={self.wexponent:4.1f}  OffMomentumWeight ={self.offmw:8.3f}{sexp}  ')
                    fitflg = False
                    exit9000 = True
                    break
                if not fitflg:
                    if r.nstab == 0:
                        out(f' Residual ={fmt_g(r.r)}  DP ={self.dpmax:8.5f}  DP0 ={self.dp0:8.5f}'
                            f'  ExponentOfResidual ={self.wexponent:4.1f}  OffMomentumWeight ={self.offmw:8.3f}{sexp}  ')
                    else:
                        out(f'Unstable = {r.nstab:5d} Residual ={fmt_g(r.r)}  DP ={self.dpmax:8.5f}  DP0 ={self.dp0:8.5f}'
                            f'  ExponentOfResidual ={self.wexponent:4.1f}  OffMomentumWeight ={self.offmw:8.3f}{sexp}  ')
                    exit9000 = True
                    break
                # ---- do1082
                exit200 = False
                while True:                          # do1082 (1 回だけ回る。exit で抜ける)
                    iter_ += 1
                    if chgmod:
                        fact = min(1.0, fact * 2.0)
                        f1 = 0.0
                        g1 = r.r
                        chgmod = False
                        aimprv = 0.0
                        crate = 1.0
                        if r.le(rbest):
                            rbest = r
                            rp0 = r.r
                            rprint = rbest
                            bestval = beforecal.copy()
                            ra = rbest.r * (1.0 + AMTOL)
                    else:
                        imprv = r.le(rbest)
                        if imprv:
                            if r.le(rprint, RTOL1):
                                if outt:
                                    outt = False
                                self.lout_on = True
                                out(f'{iter_:8d}   {r.nstab:6d}  {fmt_g(r.r)}' + ('  (NEWTON)  ' if newton else '  (DESCEND) ')
                                    + f'{fmt_g(fact)}{nvara:7d}')
                                self.nmes = 0
                                rprint = r
                            aimprv = _fuzz(_log10((ra - r.r) / ra), AIMP1, AIMP2)
                            crate = (rbest.r - r.r) / rbest.r if rbest.r else 0.0
                            fact = min(1.0, fact * (4.0 if newton else 2.0))
                            f1 = 0.0
                            g1 = r.r
                            bestval = beforecal.copy()
                            rp0 = r.r
                            rbest = r
                        else:
                            aimprv = 0.0
                            crate = (rbest.r - r.r) / rbest.r if rbest.r else 0.0
                        alate = _fuzz(float(iter_), aitm1, aitm2)
                        smallf = 1.0 - _fuzz(_log10(fact), FLIM1, FLIM2)
                        badcnv = 1.0 - _fuzz(_log10(crate), BADC1, BADC2)
                        amedcv = 1.0 - _fuzz(_log10(crate), AMEDC1, AMEDC2)
                        if _DEBUG:
                            print(f'[match6]   imprv={imprv} aimprv={aimprv:.3g} crate={crate:.3g} alate={alate:.3g} smallf={smallf:.3g} badcnv={badcnv:.3g} amedcv={amedcv:.3g}')
                        if newton:
                            chgmod = max(smallf, badcnv, min(alate, amedcv)) > 0.5
                        else:
                            if iter_ > self.itmax * 10:
                                fitflg = False
                                chgmod = False
                            elif aimprv > 0.5:
                                chgmod = True
                            elif max(smallf, badcnv, min(alate, amedcv)) > 0.5:
                                chgmod = True
                                if nretry > 0:
                                    nretry -= 1
                                else:
                                    fitflg = False
                                    chgmod = False
                            else:
                                chgmod = False
                        if chgmod or not fitflg:
                            r = rbest
                            newton = not newton
                            fact = 1.0
                            for i, v in enumerate(self.variables):
                                v['val'] = bestval[i]
                            exit200 = True
                            break
                        elif not imprv:
                            if newton:
                                f2 = f1
                                g2 = g1
                                f1 = fact
                                g1 = r.r
                                dg_ls = self._last_dg
                                fact = _fmincube(f1, f2, g1, g2, rp0, dg_ls)
                            else:
                                f1 = fact
                                fact = fact * 0.5
                            if nvara == nvar:
                                a = fact / f1
                                for i, v in enumerate(self.variables):
                                    v['val'] = v['val'] * a + bestval[i] * (1.0 - a)
                                exit200 = True
                                break
                            else:
                                dfv = df0.copy()
                                ids = list(ids0)
                                nvara = nvar
                                wlimit = np.maximum(WLMIN, wlimit)
                                break        # exit do1082 → 1082 (solve again with same qu)
                    # 微分 (tffsmatch.f の tffsqu / 数値微分): 既定は解析微分、FFS$NumericalDerivative で中心差分
                    for kc, v in enumerate(self.variables):
                        wvar[kc] = self._tweigh(v, bestval[kc])
                    qu = self._derivative(ids, dfs, rows, wvar)
                    self.last_qu = (list(ids), qu.copy(), wvar.copy())
                    df0 = dfv.copy()
                    ids0 = list(ids)
                    wlimit = np.ones(nvar)
                    nvara = nvar
                    break
                if exit200:
                    break
                # ---- 1082: 解く
                while True:
                    if newton:
                        dval, dg = self._tfsolv(qu, dfv, wlimit, ids, fact)
                        if self.wexponent != 2.0:
                            dg = dg * (rp0 / self.wsum)**(1.0 - self.wexponent / 2.0)
                        self._last_dg = dg
                        if dg > 0.0:
                            newton = False
                            dfv = df0.copy()
                            wlimit = np.ones(nvar)
                            nvara = nvar
                            continue
                    else:
                        dval = self._tgrad(qu, dfv, wlimit)
                        self._last_dg = 0.0
                    if _DEBUG:
                        print(f'[match6] DBGC iter={iter_} newton={newton} fact={fact!r} dg={self._last_dg!r} rp0={rp0!r} nvara={nvara}')
                        for _nm, _a in (('DBGD', dval), ('DBGW', wvar), ('DBGV', bestval), ('DBGL', wlimit)):
                            print(f'[match6] {_nm} ' + ' '.join(repr(float(x)) for x in _a))
                        if _os.environ.get('PYSAD_MATCH_DEBUG') == '2':
                            print(f'[match6] DBGN {iter_} {qu.shape[0]} {qu.shape[1]}')
                            print('[match6] DBGF ' + ' '.join(repr(float(x)) for x in dfv))
                            for _kc in range(qu.shape[1]):
                                print(f'[match6] DBGQ {_kc + 1} ' + ' '.join(repr(float(x)) for x in qu[:, _kc]))
                    limited = False
                    for ii, v in enumerate(self.variables):
                        dv = dval[ii] * fact / wvar[ii] * wlimit[ii]
                        v['val'] = bestval[ii] + dv
                        vl, vl1, vl2, lim1 = self._vlimit(v, v['val'], bestval[ii])
                        if lim1:
                            v['val'] = min(vl2, max(vl1, bestval[ii]))
                            if dv != 0.0:
                                wlimit[ii] = wlimit[ii] * min(abs((vl - v['val']) / dv), 0.3)
                                limited = True
                                if wlimit[ii] < WLMIN:
                                    wlimit[ii] = 0.0
                                    nvara -= 1
                    if limited and nvara > 0:
                        dfv = df0.copy()
                        continue
                    break
            self._apply_vars()
            if self.dlim:
                # tffsmatch.f:519: VariableRange が効いているときは、反復の終わりごとに tffssetlimit で現在値を範囲に収める
                lim_any = False
                for v in self.variables:
                    vl, _, _, lim = self._vlimit(v, v['val'], v['val'])
                    if lim:
                        v['val'] = vl
                        lim_any = True
                if lim_any:
                    self._apply_vars()
            if exit9000:
                break
        self.residual = r
        self.matching_residual = r.res2r()
        self.unstable_level = r.nstab
        self.iterations = iter_
        self.nmes = 0
        return r

    # ------------------------------------------------------------------ 微分 (qu)
    numerical_derivative = False      # True: FFS$NumericalDerivative=True と同じ中心差分のみ
    couples_active = False            # ffscmd が COUPLE の有無を設定 (解析微分は coup=1 のみ対応)
    last_derivative_mode = ''         # 直近の微分の種類 ('analytic' / 'numeric' / 'analytic+numeric') とフォールバック理由
    last_fallback_reason = ''

    def _analytic_unsupported(self, ids) -> Optional[str]:
        """解析微分 (match6d) が対応しない条件・状況なら理由の文字列、対応なら None (その場合は全体を中心差分にする)"""
        from .match6d import _KF_INDEX, MFITTRY, LENG
        if not self.cell and _os.environ.get('PYSAD_INS_NUMERIC'):
            return 'INS (CELL でない) は PYSAD_INS_NUMERIC=1 で数値微分に切替'
        if self.couples_active:
            return 'COUPLE 使用中 (coup≠1 は未対応)'
        for (cid, idp) in ids:
            if idp != 0:
                return 'オフモメンタム条件 (FITP>1) は未検証'
            if cid < 0:
                continue                                   # FitFunction は片側差分で別処理
            kf = _KF_INDEX.get(self._cond(cid)['kf'])
            if kf is None or (kf > MFITTRY and kf != LENG):
                return f"条件 {self._cond(cid)['kf']} は未対応"
        for v in self.variables:
            if any(getattr(e, '_errk', 1.0) != 1.0 for e in v['els']):
                return '誤差 (DK) 付きの変数は未対応'
            for e in v['els']:
                if v['type'] == 'Quadrupole' and getattr(e, 'L', 1.0) == 0.0:
                    return 'L=0 の QUAD (qdthin) は未対応'
        return None

    def _derivative(self, ids, dfs, rows, wvar):
        """qu[条件, 変数] = ∂(df)/∂v の符号・重み付き行列 (tffsmatch.f 395-460)。
        既定は match6d.analytic_qu (解析微分)。FFS$NumericalDerivative=True、未対応の条件・状況、未対応の変数列は中心差分。
        FitFunction の行は SAD と同じ v+eps1/wvar の片側差分。"""
        nvar = len(self.variables)
        if self.numerical_derivative:
            self.last_derivative_mode, self.last_fallback_reason = 'numeric', ''
            return self._numeric_columns(ids, dfs, wvar, range(nvar))
        reason = self._analytic_unsupported(ids)
        if reason is not None:
            self.last_derivative_mode, self.last_fallback_reason = 'numeric', reason
            return self._numeric_columns(ids, dfs, wvar, range(nvar))
        from .match6d import analytic_qu
        qu, numeric_cols = analytic_qu(self, ids, rows, wvar)
        self.last_derivative_mode, self.last_fallback_reason = 'analytic', ''
        if numeric_cols:
            qn = self._numeric_columns(ids, dfs, wvar, numeric_cols)
            for kc in numeric_cols:
                qu[:, kc] = qn[:, kc]
            self.last_derivative_mode = 'analytic+numeric'
        ff_rows = [i for i, (cid, _) in enumerate(ids) if cid < 0]
        if ff_rows:
            # tffsmatch.f 420-460 行目: FitFunction の行は v + eps1/wvar での片側差分
            for kc, v in enumerate(self.variables):
                if kc in numeric_cols:
                    continue
                valbest = v['val']
                v['val'] = valbest + EPS1 / wvar[kc]
                dfs1, _, _ = self._tffscalc(False)
                v['val'] = valbest
                for i in ff_rows:
                    cid = ids[i]
                    qu[i, kc] = (dfs[cid] - dfs1.get(cid, dfs[cid])) / EPS1
            self._apply_vars()
        return qu

    def _numeric_columns(self, ids, dfs, wvar, cols):
        """数値微分 (中心差分、tffsddf: df − df1。無い条件は 0) を指定の列だけ。他の列は 0"""
        nvar = len(self.variables)
        qu = np.zeros((len(ids), nvar))
        for kc in cols:
            v = self.variables[kc]
            dvkc = max(abs(v['val']) * EPS, abs(EPS1 / wvar[kc]))
            valbest = v['val']
            v['val'] = valbest + dvkc
            dfs1, _, _ = self._tffscalc(False)
            v['val'] = valbest - dvkc
            dfs2, _, _ = self._tffscalc(False)
            v['val'] = valbest
            for i, cid in enumerate(ids):
                ddf1 = dfs[cid] - dfs1.get(cid, dfs[cid])
                ddf2 = dfs[cid] - dfs2.get(cid, dfs[cid])
                qu[i, kc] = (ddf1 - ddf2) / 2.0 / dvkc / wvar[kc]
        self._apply_vars()
        return qu

    def values(self) -> Dict[str, float]:
        return {f"{v['family']} {v['key']}": v['val'] for v in self.variables}
