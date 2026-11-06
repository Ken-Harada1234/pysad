"""
emit6d.py — SAD の Emittance[] (6D 平衡エミッタンス) の直訳
==========================================================

SAD での呼び出し経路:
  Emittance[] (Packages/emit.n) → CalculateEmittance → temit (temit.f)
    → tturne (tturne.f): 放射込み 6D 一周追跡
         各要素: tdrife0 / tquade / tbende ...
         放射点: tsetr0 で基準を保存 → tradke (photons.f) で減衰と量子励起を加算
    → tecalc (temit.f): 固有値分解 + 平衡ビーム行列を SVD で解く

trans(6,12) の構造 (SAD と同じ):
  T = trans(:,1:6)   放射なしのシンプレクティック部分
  D = trans(:,7:12)  放射減衰による摂動 (一周行列 = T + D)
  B = beam(1:21)     量子励起による 2 次モーメント (ここでは 6x6 対称行列で保持)

対応範囲 (bench_level2 で使う範囲):
  DRIFT, QUAD (bz=0, 誤差なし), BEND (K1=0 / K1≠0, F1=0, FRINGE=0, 誤差なし), MARK
  RFSW=F (空洞なし), RADCOD=F, INTRA/WSPAC/偏極計算なし
要素写像は lattice.py の COD 写像 (SAD 直訳) を使い、Jacobian は数値微分。
K1≠0 の BEND 本体 (tbendecorr/tbendebody) は SAD の解析行列もそのまま直訳。
"""

import math
import os
from typing import List, Optional

import numpy as np

from . import lattice as _lat
from .elements import Element, Drift, Quadrupole, Bend
from .perfkern import tradke_db as _tradke_db
from .perfkern_bend import apply_tdb as _apply_tdb, bendecorr as _k_bendecorr, bendebody as _k_bendebody
from . import elkern as _ek

_USE_ELKERN = os.environ.get('PYSAD_NO_ELKERN') != '1'    # False: 要素カーネルを使わず従来の Python の piece で追う (比較用)

# ---------------------------------------------------------------------------
# 物理定数 (tfstk.f)
# ---------------------------------------------------------------------------
CVELOC = 299792458.0
PLANK = 6.62607015e-34
PLANKR = PLANK / (2.0 * math.pi)
ELEMCH = 1.602176634e-19
FINEST = 1.0 / 137.035999084
ELMASS = 0.51099895000e6
ELRADI = FINEST * PLANKR * CVELOC / (ELEMCH * ELMASS)

# photons.f / tmacro / bendib
CAVE = 8.0 / 15.0 / math.sqrt(3.0)
CUU = 11.0 / 27.0
GMIN = -0.999999
EPSRAD = 1e-6
EPSBEND = 1e-3
PSQMAX = 0.9999


def _p2h(p):
    """mathfun p2h: sqrt(1+p^2)"""
    if p > 1e3:
        p2 = 1.0 / p**2
        return p * (1.0 + p2 * (0.5 - p2 * 0.125))
    return math.hypot(1.0, p)


def _pxy2dpz(px, py):
    """mathfun pxy2dpz: sqrt(1-px^2-py^2)-1"""
    x = px**2 + py**2
    if x < 8e-6:
        return -x * (0.5 + x * (0.125 + x * 0.0625))
    return -x / (1.0 + math.sqrt(max(1e-200, 1.0 - x)))


def _sqrt1(x):
    """mathfun sqrt1: sqrt(1+x)-1"""
    if abs(x) < 8e-6:
        return x * (0.5 + x * (-0.125 + x * 0.0625))
    if abs(x) < 3.4e-3:
        return x * (0.5 + x * (-0.125 + x * (0.0625 + x * (
            -.0390625 + x * (.02734375 - .0205078125 * x)))))
    return x / (1.0 + math.sqrt(max(1e-200, 1.0 + x)))


def _asinz(x):
    return math.asin(max(-1.0, min(1.0, x)))


def _tinv6(ra):
    """toplvl.f tinv6: シンプレクティック逆行列 -J r^T J"""
    J = np.zeros((6, 6))
    for i in range(0, 6, 2):
        J[i, i + 1] = 1.0
        J[i + 1, i] = -1.0
    return -J @ ra.T @ J


def _tsymp(t):
    """toplvl.f tsymp: シンプレクティック化"""
    rx = t @ _tinv6(t)
    return -0.5 * (rx @ t) + 1.5 * t


def _iaidx(i, j):
    """SAD iaidx (0-based 入力, 0-based 出力): 対称行列の下三角パック番号"""
    if i < j:
        i, j = j, i
    return i * (i + 1) // 2 + j


def _pack(B):
    v = np.zeros(21)
    for i in range(6):
        for j in range(i + 1):
            v[_iaidx(i, j)] = B[i, j]
    return v


def _unpack(v):
    B = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            B[i, j] = B[j, i] = v[_iaidx(i, j)]
    return B


# ---------------------------------------------------------------------------
# tphyzp: エネルギーから決まる定数
# ---------------------------------------------------------------------------

class _Phys:
    def __init__(self, pgev: float, circ: float, amass: float = ELMASS,
                 charge: float = 1.0, emidiv: float = 1.0,
                 emidib: float = 1.0, emidiq: float = 1.0):
        hp = PLANKR
        e = ELEMCH
        c = CVELOC
        self.pgev = pgev
        self.amass = amass
        self.charge = charge
        self.brhoz = pgev / c
        self.p0 = pgev / amass
        self.h0 = _p2h(self.p0)
        re0 = ELRADI
        rclassic = charge**2 * re0
        sgn = 1.0 if charge >= 0 else -1.0
        self.crad = sgn * abs(rclassic * (c / amass)**2 / self.p0 / 1.5)
        self.urad = sgn * abs(1.5 * hp * c / self.p0 / amass / e)
        rcratio = rclassic / (hp * c / amass / e)
        self.cuc = 1.5 * rclassic / rcratio
        self.anrad = 5.0 / 2.0 / math.sqrt(3.0) * rcratio
        self.omega0 = (2.0 * math.pi * c * self.p0 / self.h0 / circ
                       if circ != 0 else 0.0)
        self.emidiv = emidiv
        self.emidib = emidib
        self.emidiq = emidiq
        self.dvfs = 0.0          # tsetdvfs: RFSW ON のとき FSHIFT (RF 周波数のずらし = Orbit dilation)。emittance6d が設定する

    def dv(self, dp):
        """tesetdv / tgetdv: dvemit = -(p1+p0)/h1/(p0 h1+p1 h0) dp + dvfs"""
        p1 = self.p0 * (1.0 + dp)
        h1 = _p2h(p1)
        return -(p1 + self.p0) / h1 / (self.p0 * h1 + p1 * self.h0) * dp + self.dvfs


# ---------------------------------------------------------------------------
# 追跡状態 (trans(6,12), cod, beam, 放射基準点)
# ---------------------------------------------------------------------------

class _State:
    def __init__(self, cod, ph: _Phys):
        self.T = np.eye(6)
        self.D = np.zeros((6, 6))
        self.B = np.zeros((6, 6))
        self.cod = np.array(cod, dtype=float)
        self.ph = ph
        # RF (tcave / tturneg の積算量。tmacro の vc0, hvc0, vcacc, dvcacc, ddvcacc)
        self.trf0 = 0.0
        self.rfsw = False
        self.radcod = False
        self.vc0 = 0.0
        self.hvc0 = 0.0
        self.vcacc = 0.0
        self.dvcacc = 0.0
        self.ddvcacc = 0.0
        self.codr0 = self.cod.copy()
        self.transr = self.T.copy()
        self.want_jac = True    # False: cod だけ追う (ヤコビアン T も作らない)。閉軌道の収束判定だけに使う反復用
        self.want_db = True     # False: SAD の irad=6 (tcod の反復用)。放射による cod の変化と u0 だけ追い、D・B 行列は作らない
        p1 = ph.p0 * (1.0 + self.cod[5])
        self._S = np.array([0.0, _p2h(p1), ph.dv(self.cod[5])])   # [u0, h1emit, dvemit] (elkern のカーネルと共有)
        self.PH = _ek.phys_array(ph)

    @property
    def u0(self):
        return self._S[0]

    @u0.setter
    def u0(self, v):
        self._S[0] = v

    @property
    def h1emit(self):
        return self._S[1]

    @h1emit.setter
    def h1emit(self, v):
        self._S[1] = v

    @property
    def dvemit(self):
        return self._S[2]

    @dvemit.setter
    def dvemit(self, v):
        self._S[2] = v

    # tmultr5 + tmulbs
    def apply(self, J):
        # in-place (perfkern_bend.apply_tdb): T ← J·T、want_db なら D ← J·D、B ← J·B·J^T
        _apply_tdb(J, self.T, self.D, self.B, self.want_db)

    def piece(self, f, symp_z=True, jac=None):
        """COD 写像 f を適用し、COD 位置での Jacobian (jac(cod) の解析式、なければ数値微分) で行列を伝搬。
        z の行 (列 x,px,y,py) は SAD と同じくシンプレクティック条件から組み立てる:
          trans1(5,j) = -t(1,j)t(2,6) + t(2,j)t(1,6) - t(3,j)t(4,6) + t(4,j)t(3,6)
        (tqfrie.f, tbende.f tbendebody0, tsolque.f 等)。SAD の z 写像は近似式を含むので、
        数値微分のままだと一周行列がシンプレクティックからずれ、軌道がずれた場合の
        εy が 2% 程度狂う (Level 3C で確認)。(5,6) 成分は数値微分のまま。"""
        if not self.want_jac:
            self.cod = f(self.cod)
            return
        J = jac(self.cod) if (jac is not None and _USE_ANALYTIC_JAC) else None
        if J is None:
            J = _lat._num_jacobian(f, self.cod)
        if symp_z:
            for j in range(4):
                J[4, j] = (-J[0, j] * J[1, 5] + J[1, j] * J[0, 5]
                           - J[2, j] * J[3, 5] + J[3, j] * J[2, 5])
        self.cod = f(self.cod)
        self.apply(J)

    # temw tsetr0
    def setr0(self):
        self.codr0 = self.cod.copy()
        self.transr = self.T.copy()

    # photons.f tradke (bzh=0, radcod=F, calpol=F, irad=12)
    def tradke(self, al, phir0):
        if not self.want_jac and not self.radcod:
            return          # cod も変えず、u0・D・B も使わない一周 (z_slip 等)
        ph = self.ph
        cod = self.cod
        codr0 = self.codr0
        gi = codr0[5]
        pr = 1.0 + gi
        th = math.tan(0.5 * phir0)
        sp = 2.0 * th / (1.0 + th**2)
        cp = 1.0 - th * sp
        pxi = codr0[1]
        pyi = codr0[3]
        pz00 = pr * (1.0 + _pxy2dpz(pxi / pr, pyi / pr))
        pxr0 = cp * pxi + sp * pz00
        pz0 = -sp * pxi + cp * pz00
        px = cod[1]
        py = cod[3]
        pz = pr * (1.0 + _pxy2dpz(px / pr, py / pr))
        dpx = px - pxr0
        dpy = py - pyi
        xpx = (py * pz0 - pz * pyi)
        xpy = (pz * pxr0 - px * pz0)
        xpz = (px * pyi - py * pxr0)
        xpa = math.sqrt(xpx**2 + xpy**2 + xpz**2) / pr**2
        theta = _asinz(xpa)
        p = ph.p0 * pr
        h1 = _p2h(p)
        al1 = al - cod[4] + codr0[4]
        anp = ph.anrad * h1 * theta
        uct = ph.cuc * h1**3 / ph.p0 / al1
        uc = uct * theta
        dg = -CAVE * anp * uc
        self.u0 -= dg
        ddpx = 0.5 * dpx * dg
        ddpy = 0.5 * dpy * dg
        c1 = al1 / pr / 3.0
        if self.radcod:
            # photons.f tradke (radcod): 放射で cod を変える (bzh=0)
            g = max(GMIN, gi + dg)
            pr1 = 1.0 + g
            cod = self.cod.copy()
            cod[0] += ddpx * c1
            cod[2] += ddpy * c1
            cod[1] = px * pr1 / pr + ddpx
            cod[3] = py * pr1 / pr + ddpy
            cod[5] = g
            p2 = ph.p0 * pr1
            h2 = _p2h(p2)
            cod[4] = cod[4] * p2 / h2 * h1 / p
            self.cod = cod
            # tesetdv(g)
            self.h1emit = _p2h(p2)
            self.dvemit = ph.dv(g)
        # radcod=F: cod は変えない
        if not self.want_db:
            self.codr0 = self.cod.copy()
            self.transr = self.T.copy()
            return
        # D・B の更新は numba カーネル (perfkern.tradke_db: 下の numpy 版と同じ式。numpy の小行列演算の連続が 0.5 ms/回かかるため)
        ok, Dn, Bn = _tradke_db(self.T, self.transr, self.D, self.B, pr, px, py, pxi, pyi, pz, pz00, pz0, pxr0, cp, sp,
                                xpx, xpy, xpz, xpa, anp, uc, uct, al1, dg, ddpx, ddpy, c1, h1, p, ph.p0, ph.anrad, CAVE, CUU)
        if ok:
            self.D = Dn
            self.B = Bn
        self.codr0 = self.cod.copy()
        self.transr = self.T.copy()


# ---------------------------------------------------------------------------
# 解析ヤコビアン: st.piece(f, jac=...) に渡す jac(c) = 写像 f の COD 位置 c での 6x6 微分
# (数値微分 _num_jacobian の置き換え。12 回の写像評価が要らないので 3〜5 倍速い)。
# 写像 f は変えない。z の行 (列 x,px,y,py) は piece() が SAD と同じくシンプレクティック条件から組み直すので、
# ここで要るのは x,px,y,py の行と z の (5,6) 成分だけ (z 列は z に依存しない写像なので単位)。
# 検証: pysad/tests/test_emit_jac.py が各ヤコビアンを数値微分と比べる (相対 1e-9 以内)。
# 解析式が使えない状況 (運動量が pr 以上でクリップされる等) は None を返し、piece() が数値微分に戻る。
# ---------------------------------------------------------------------------

_EYE6 = np.eye(6)
_USE_ANALYTIC_JAC = True      # False にすると解析ヤコビアンを使わず従来どおり数値微分 (テスト・デバッグ用)


def _dvz(ph, dp):
    """z -= dv(δ)·L の δ 微分 ∂z/∂δ = +(h0/h1³)·L の係数。h1 = √(1+(p0(1+δ))²)。
    SAD の h0/h1emit**3 (tbende.f の trans(5,6)) で、-dv'(δ) と代数的に一致する (数値で 1e-12 まで確認)。"""
    h1 = _p2h(ph.p0 * (1.0 + dp))
    return ph.h0 / (h1 * h1 * h1)


def _jac_drift(c, L, dvz=0.0):
    """厳密ドリフト (lattice._exact_drift_cod / emit6d._drift_map の dv 項を除いた部分) の微分。
    dvz は _dvz(...) (dv 項を足すときだけ指定): ∂z/∂δ に dvz·L を加える。
      x += px L/pz, y += py L/pz, z += L(1 - pr/pz) (pz = √(pr² - px² - py²))"""
    pr = 1.0 + c[5]
    px = c[1]
    py = c[3]
    pz2 = pr * pr - px * px - py * py
    if pz2 <= 1e-12:
        return None
    k = L / (pz2 * math.sqrt(pz2))
    J = _EYE6.copy()
    J[0, 1] = k * (pz2 + px * px)
    J[0, 3] = k * px * py
    J[0, 5] = -k * px * pr
    J[2, 1] = J[0, 3]
    J[2, 3] = k * (pz2 + py * py)
    J[2, 5] = -k * py * pr
    J[4, 1] = -k * pr * px
    J[4, 3] = -k * pr * py
    J[4, 5] = k * (px * px + py * py) + dvz * L
    return J


def _jac_kick(c, aki, kord):
    """薄い多極のキック px -= aki Re((x-iy)^k), py -= aki Im((x-iy)^k) の微分"""
    g = kord * complex(c[0], -c[2]) ** (kord - 1)
    J = _EYE6.copy()
    J[1, 0] = -aki * g.real
    J[1, 2] = -aki * g.imag
    J[3, 0] = -aki * g.imag
    J[3, 2] = aki * g.real
    return J


def _jac_thin(c, L, ak, kord, dvz):
    """lattice._tthine_cod (ドリフト L/6 → キック → ドリフト 2L/3 → キック → ドリフト L/6) + 全長 L の dv 項"""
    if ak == 0.0:
        return _jac_drift(c, L, dvz)
    aki = ak / math.factorial(kord)
    if L == 0.0:
        return _jac_kick(c, aki, kord)
    J = _jac_drift(c, L / 6.0)
    if J is None:
        return None
    c1 = _lat._exact_drift_cod(c, L / 6.0)
    K = _jac_kick(c1, aki * 0.5, kord)
    J = K @ J
    c1 = c1.copy()
    cx = complex(c1[0], -c1[2]) ** kord
    c1[1] -= aki * 0.5 * cx.real
    c1[3] -= aki * 0.5 * cx.imag
    J2 = _jac_drift(c1, L / 1.5)
    if J2 is None:
        return None
    J = J2 @ J
    c1 = _lat._exact_drift_cod(c1, L / 1.5)
    K = _jac_kick(c1, aki * 0.5, kord)
    J = K @ J
    c1 = c1.copy()
    cx = complex(c1[0], -c1[2]) ** kord
    c1[1] -= aki * 0.5 * cx.real
    c1[3] -= aki * 0.5 * cx.imag
    J3 = _jac_drift(c1, L / 6.0)
    if J3 is None:
        return None
    J = J3 @ J
    J[4, 5] += dvz * L
    return J


def _jac_tqente(c, al):
    """lattice._tqente_cod: x += r px, y += r py, z -= (2+s/pr) a/(2(pr+s)) r  (s=pz, r = al(1/s - 1/pr), a=px²+py²)"""
    pr = 1.0 + c[5]
    px = c[1]
    py = c[3]
    a = px * px + py * py
    s2 = pr * pr - a
    if s2 <= 1e-12:
        return None
    s = math.sqrt(s2)
    s3 = s2 * s
    r = al * (1.0 / s - 1.0 / pr)
    r_px = al * px / s3
    r_py = al * py / s3
    r_pr = al * (-pr / s3 + 1.0 / (pr * pr))
    J = _EYE6.copy()
    J[0, 1] = r + px * r_px
    J[0, 3] = px * r_py
    J[0, 5] = px * r_pr
    J[2, 1] = py * r_px
    J[2, 3] = r + py * r_py
    J[2, 5] = py * r_pr
    # z' = z - G r,  G = g1 g2,  g1 = 2 + s/pr,  g2 = a/(2(pr+s))
    g1 = 2.0 + s / pr
    g2 = a / (2.0 * (pr + s))
    G = g1 * g2
    # ds/dpx = -px/s, ds/dpy = -py/s, ds/dpr = pr/s
    g1_pr = 1.0 / s - s / (pr * pr)
    g2_pr = -a / (2.0 * (pr + s) ** 2) * (1.0 + pr / s)
    G_pr = g1_pr * g2 + g1 * g2_pr
    J[4, 5] = -(G_pr * r + G * r_pr)
    return J


def _jac_tqfrie(c, ak, al):
    """lattice._tqfrie_cod (四極のハードエッジフリンジ) の微分。ak = 積分 K1 (入口 +, 出口 -)。"""
    pr = 1.0 + c[5]
    akk = ak / al / 4.0
    k = akk / pr
    x = c[0]
    y = c[2]
    px = c[1]
    py = c[3]
    x2 = x * x
    y2 = y * y
    a = k * x2
    b = k * y2
    ab = a - b
    t = ab * ab / 6.0
    X1 = 1.0 + a / 3.0 + b + t
    Y1 = 1.0 - a - b / 3.0 + t
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    q = 1.0 - ab / 3.0
    h = 2.0 * k * x * y * q
    f = xx * yy + h * h
    if f == 0.0:
        return None
    # 変数 v = x, y, k (k = akk/pr → ∂k/∂pr = -k/pr) での偏微分
    # a_v, b_v, ab_v, t_v
    a_x, a_y, a_k = 2.0 * k * x, 0.0, x2
    b_x, b_y, b_k = 0.0, 2.0 * k * y, y2
    ab_x, ab_y, ab_k = a_x - b_x, a_y - b_y, a_k - b_k
    t_x, t_y, t_k = ab * ab_x / 3.0, ab * ab_y / 3.0, ab * ab_k / 3.0
    X1_x, X1_y, X1_k = a_x / 3.0 + b_x + t_x, a_y / 3.0 + b_y + t_y, a_k / 3.0 + b_k + t_k
    Y1_x, Y1_y, Y1_k = -a_x - b_x / 3.0 + t_x, -a_y - b_y / 3.0 + t_y, -a_k - b_k / 3.0 + t_k
    xx_x = a_x + b_x + (ab_x * (5.0 * a - b) + ab * (5.0 * a_x - b_x)) / 6.0
    xx_y = a_y + b_y + (ab_y * (5.0 * a - b) + ab * (5.0 * a_y - b_y)) / 6.0
    xx_k = a_k + b_k + (ab_k * (5.0 * a - b) + ab * (5.0 * a_k - b_k)) / 6.0
    yy_x = -(a_x + b_x) + (ab_x * (a - 5.0 * b) + ab * (a_x - 5.0 * b_x)) / 6.0
    yy_y = -(a_y + b_y) + (ab_y * (a - 5.0 * b) + ab * (a_y - 5.0 * b_y)) / 6.0
    yy_k = -(a_k + b_k) + (ab_k * (a - 5.0 * b) + ab * (a_k - 5.0 * b_k)) / 6.0
    h_x = 2.0 * k * y * q - 2.0 * k * x * y * ab_x / 3.0
    h_y = 2.0 * k * x * q - 2.0 * k * x * y * ab_y / 3.0
    h_k = 2.0 * x * y * q - 2.0 * k * x * y * ab_k / 3.0
    f_x = xx_x * yy + xx * yy_x + 2.0 * h * h_x
    f_y = xx_y * yy + xx * yy_y + 2.0 * h * h_y
    f_k = xx_k * yy + xx * yy_k + 2.0 * h * h_k
    pxn = (px * yy + py * h) / f
    pyn = (py * xx - px * h) / f
    # 位置
    xn_x = X1 + x * X1_x
    xn_y = x * X1_y
    xn_k = x * X1_k
    yn_x = y * Y1_x
    yn_y = Y1 + y * Y1_y
    yn_k = y * Y1_k
    # 運動量
    pxn_x = (px * yy_x + py * h_x - pxn * f_x) / f
    pxn_y = (px * yy_y + py * h_y - pxn * f_y) / f
    pxn_k = (px * yy_k + py * h_k - pxn * f_k) / f
    pyn_x = (py * xx_x - px * h_x - pyn * f_x) / f
    pyn_y = (py * xx_y - px * h_y - pyn * f_y) / f
    pyn_k = (py * xx_k - px * h_k - pyn * f_k) / f
    dkdp = -k / pr
    J = _EYE6.copy()
    J[0, 0] = xn_x
    J[0, 2] = xn_y
    J[0, 5] = xn_k * dkdp
    J[2, 0] = yn_x
    J[2, 2] = yn_y
    J[2, 5] = yn_k * dkdp
    J[1, 0] = pxn_x
    J[1, 2] = pxn_y
    J[1, 1] = yy / f
    J[1, 3] = h / f
    J[1, 5] = pxn_k * dkdp
    J[3, 0] = pyn_x
    J[3, 2] = pyn_y
    J[3, 1] = -h / f
    J[3, 3] = xx / f
    J[3, 5] = pyn_k * dkdp
    # z' = z - (pxn A + pyn B)/pr,  A = x(a/3 + b + 2t),  B = y(-a - b/3 + 2t)
    A = x * (a / 3.0 + b + 2.0 * t)
    B = y * (-a - b / 3.0 + 2.0 * t)
    A_k = x * (a_k / 3.0 + b_k + 2.0 * t_k)
    B_k = y * (-a_k - b_k / 3.0 + 2.0 * t_k)
    S = pxn * A + pyn * B
    S_pr = (pxn_k * A + pxn * A_k + pyn_k * B + pyn * B_k) * dkdp
    J[4, 5] = S / (pr * pr) - S_pr / pr
    return J


def _jac_tsolque_body(c, aln1, kappa, dvz):
    """lattice._tsolque_body_cod(c, aln1, kappa/pr) + z -= dv(δ)·aln1 の微分 (kappa = akn/aln = ak/L·…, w² = kappa/pr)。
    x: 集束 (cos/sin)、y: 発散 (cosh/sinh)。変数は x, p=px/pr, y, q=py/pr, w=√(kappa/pr), φ=aln1·w。"""
    pr = 1.0 + c[5]
    x = c[0]
    y = c[2]
    p = c[1] / pr
    q = c[3] / pr
    w = math.sqrt(kappa / pr)
    phi = aln1 * w
    if phi > 100.0:
        return None
    cs = math.cos(phi)
    sn = math.sin(phi)
    ch = math.cosh(phi)
    sh = math.sinh(phi)
    dc = cs - 1.0
    dch = ch - 1.0
    xs = phi - sn
    xsh = phi - sh
    dw = -w / (2.0 * pr)            # dw/dpr
    dphi = aln1 * dw                # dφ/dpr
    p_r = -p / pr
    q_r = -q / pr
    u1 = x * dc + p * sn / w
    pxf = p * cs - x * w * sn
    v1 = y * dch + q * sh / w
    pyf = q * ch + y * w * sh
    J = _EYE6.copy()
    # x 面
    J[0, 0] = cs
    J[0, 1] = sn / (w * pr)
    J[0, 5] = -x * sn * dphi + p_r * sn / w + p * (cs * dphi / w - sn * dw / (w * w))
    J[1, 0] = -pr * w * sn
    J[1, 1] = cs
    pxf_r = p_r * cs - p * sn * dphi - x * (dw * sn + w * cs * dphi)
    J[1, 5] = pxf + pr * pxf_r
    # y 面
    J[2, 2] = ch
    J[2, 3] = sh / (w * pr)
    J[2, 5] = y * sh * dphi + q_r * sh / w + q * (ch * dphi / w - sh * dw / (w * w))
    J[3, 2] = pr * w * sh
    J[3, 3] = ch
    pyf_r = q_r * ch + q * sh * dphi + y * (dw * sh + w * ch * dphi)
    J[3, 5] = pyf + pr * pyf_r
    # z' = z - Z/4,
    #   Z = w(x² xs - y² xsh) + (p²+q²)aln1 + u1 pxf + x p dc + v1 pyf + y q dch
    u1_r = -x * sn * dphi + p_r * sn / w + p * (cs * dphi / w - sn * dw / (w * w))
    v1_r = y * sh * dphi + q_r * sh / w + q * (ch * dphi / w - sh * dw / (w * w))
    Z_r = (dw * (x * x * xs - y * y * xsh) + w * dphi * (-x * x * dc + y * y * dch)
           + 2.0 * (p * p_r + q * q_r) * aln1
           + u1_r * pxf + u1 * pxf_r
           + x * (p_r * dc - p * sn * dphi)
           + v1_r * pyf + v1 * pyf_r
           + y * (q_r * dch + q * sh * dphi))
    J[4, 5] = -0.25 * Z_r + dvz * aln1
    return J


def _jac_tqlfre(c, al, ak, f1, f2):
    """tqlfre.f (bz=0) の trans1 の直訳: 四極の線形フリンジ (FRINGE=1,2,3)。ak = 積分 K1 (回転後は正)、f1/f2 は
    入口 (f1in, f2in)、出口は (-f1out, f2out) を渡す。行 1〜4 は cod 写像の厳密な微分、z の行は SAD のシンプレクティック条件の式のまま"""
    af1 = -ak / al * f1 * abs(f1) / 24.0
    af2 = ak / al * f2
    p = 1.0 + c[5]
    a = af1 / p
    b = af2 / p
    ea = math.exp(a)
    bp = b / p
    x, px, y, py = c[0], c[1], c[2], c[3]
    pxf = px / ea
    pyf = py * ea
    t11 = ea
    t12 = bp
    t16 = (-a * ea * x - bp * 2.0 * px) / p
    t33 = 1.0 / ea
    t34 = -bp
    t36 = (a / ea * y + bp * 2.0 * py) / p
    t22 = 1.0 / ea
    t26 = a / p * pxf
    t44 = ea
    t46 = -a / p * pyf
    J = _EYE6.copy()
    J[0, 0], J[0, 1], J[0, 5] = t11, t12, t16
    J[1, 1], J[1, 5] = t22, t26
    J[2, 2], J[2, 3], J[2, 5] = t33, t34, t36
    J[3, 3], J[3, 5] = t44, t46
    J[4, 0] = -t11 * t26
    J[4, 1] = -t12 * t26 + t22 * t16
    J[4, 2] = -t33 * t46
    J[4, 3] = -t34 * t46 + t44 * t36
    J[4, 5] = (((2.0 + a) * a * ea * x - bp * ea * (2.0 + 0.5 * a) * pxf) / p * pxf
               - ((2.0 - a) * a / ea * y - bp / ea * (2.0 - 0.5 * a) * pyf) / p * pyf
               - (a * ea * x + bp * ea * (2.0 + a) * pxf) * t26
               + (a / ea * y + bp / ea * (2.0 - a) * pyf) * t46) / p
    return J


def _roll_matrices(th):
    """tsolrote (bz=0) のロール θ の入口/出口の 6x6: x' = c x − s y, y' = s x + c y (px, py も同じ)。出口は転置"""
    ct, st = math.cos(th), math.sin(th)
    R = np.eye(6)
    R[0, 0] = R[1, 1] = R[2, 2] = R[3, 3] = ct
    R[0, 2] = R[1, 3] = -st
    R[2, 0] = R[3, 1] = st
    return R, R.T.copy()


def _jac_quad_full(c, L, K1, fringe, fr, ph, dx=0.0, dy=0.0, th=0.0):
    """lattice._tquade_cod (+ _quad_map の dv 項 z -= dv(δ)·L) の微分: tsolrote (ずれ・ロール) → [K1<0 は −π/2 回転] →
    tqfrie(+) → tqlfre(入口) → tsolque (tqente + 本体 ×ndiv) → tqlfre(出口) → tqfrie(−) → 逆回転。
    ph=None は dv 項なし (lattice の CALC 経路)。ずれ dx,dy は微分に効かない。ロール θ は R_出口 · J(R_入口 (c − ずれ)) · R_入口。"""
    f1in, f2in, f1out, f2out, mfring = fr
    if abs(K1) < 1e-10 or abs(L) < 1e-15:
        return None
    ak1 = abs(K1) * L
    rot = K1 * L < 0.0
    Rin = Rout = None
    cc = c
    if dx != 0.0 or dy != 0.0 or th != 0.0:
        cc = c.copy()
        cc[0] -= dx
        cc[2] -= dy
        if th != 0.0:
            Rin, Rout = _roll_matrices(th)
            cc = Rin @ cc
    J = _EYE6 if Rin is None else Rin
    if rot:
        cc = cc.copy()
        cc[0], cc[1], cc[2], cc[3] = cc[2], cc[3], -cc[0], -cc[1]
        J = _ROT_IN @ J
    dp = c[5]
    if fringe and mfring >= 0 and mfring != 2:
        Jf = _jac_tqfrie(cc, ak1, L)
        if Jf is None:
            return None
        J = Jf @ J
        cc = _lat._tqfrie_cod(cc, ak1, L)
    if mfring in (1, 3):
        J = _jac_tqlfre(cc, L, ak1, f1in, f2in) @ J
        cc = _lat._tqlfre_cod(cc, L, ak1, f1in, f2in)
    kappa = abs(K1)
    ndiv = min(1000, 1 + int(abs(L * ak1) / 0.1))
    aln1 = L / ndiv
    al1 = aln1 * 0.5
    for _n in range(ndiv):
        J1 = _jac_tqente(cc, al1)
        if J1 is None:
            return None
        c1 = _lat._tqente_cod(cc, al1)
        pr = 1.0 + dp
        J2 = _jac_tsolque_body(c1, aln1, kappa, 0.0)
        if J2 is None:
            return None
        J = J2 @ (J1 @ J)
        cc = _lat._tsolque_body_cod(c1, aln1, kappa / pr)
        al1 = aln1
    J1 = _jac_tqente(cc, aln1 * 0.5)
    if J1 is None:
        return None
    J = J1 @ J
    cc = _lat._tqente_cod(cc, aln1 * 0.5)
    if mfring in (2, 3):
        J = _jac_tqlfre(cc, L, ak1, -f1out, f2out) @ J
        cc = _lat._tqlfre_cod(cc, L, ak1, -f1out, f2out)
    if fringe and mfring >= 0 and mfring != 1:
        Jf = _jac_tqfrie(cc, -ak1, L)
        if Jf is None:
            return None
        J = Jf @ J
    if rot:
        J = _ROT_OUT @ J
    if Rout is not None:
        J = Rout @ J
    J = J.copy() if J is _EYE6 else J
    if ph is not None:
        J[4, 5] += _dvz(ph, dp) * L
    return J


def _rot_in_matrix(rot):
    """tquade の K1<0 の -π/2 回転 (x,px,y,py) → (y,py,-x,-px) の 6x6 と、その逆 (出口)"""
    Ri = np.eye(6)
    Ri[0, 0] = Ri[1, 1] = Ri[2, 2] = Ri[3, 3] = 0.0
    Ri[0, 2] = 1.0
    Ri[1, 3] = 1.0
    Ri[2, 0] = -1.0
    Ri[3, 1] = -1.0
    return Ri, Ri.T.copy()


_ROT_IN, _ROT_OUT = _rot_in_matrix(True)


# ---------------------------------------------------------------------------
# 要素ごとの 6D 追跡
# ---------------------------------------------------------------------------

def _drift_map(L, ph):
    """tdrife0/tdrife1 (bz=0) COD 写像: 厳密ドリフト + dv 項"""
    def f(c):
        c = c.copy()
        pr = 1.0 + c[5]
        a = c[1]**2 + c[3]**2
        dpz = -a / pr / (1.0 + math.sqrt(max(0.0, 1.0 - a / pr**2)))
        pzi = pr + dpz
        ale = L / pzi
        c[0] += c[1] * ale
        c[2] += c[3] * ale
        c[4] += dpz * ale - ph.dv(c[5]) * L
        return c
    return f


def _tdrife0(st: _State, L):
    if L == 0.0:
        return
    ph = st.ph
    st.piece(_drift_map(L, ph), jac=lambda c: _jac_drift(c, L, _dvz(ph, c[5])))


def _quad_map(L, K1, ph, dx=0.0, dy=0.0, theta=0.0, fringe=True, fr=(0.0, 0.0, 0.0, 0.0, 0)):
    g = _lat._quad_offset_map(L, K1, dx, dy, theta, fringe, fr)

    def f(c):
        c = g(c)
        c = c.copy()
        c[4] -= ph.dv(c[5]) * L
        return c
    return f


def _tquade(st: _State, el: Quadrupole, krad: bool = True):
    """tquade (bz=0, 誤差なし, F1=0, FRINGE=0, DISFRIN=0)。
    放射: tsetr0 → tqfrie → tsetr0 → tsolque スライス毎に tradke → tqfrie
    (軌道ゼロでは tradke の寄与は厳密に 0)"""
    L = el.L
    K1 = getattr(el, 'K1', 0.0)
    ph = st.ph
    if L == 0.0 and K1 != 0.0:
        # tquade.f: al=0 → tthine(nord=4, ak=K1 (積分値), dx, dy, theta) (放射なし。K1<0 の −90° 回転はない)
        _tthine_gen(st, 0.0, float(K1), 1, float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0),
                    float(getattr(el, 'ROTATION', 0.0) or 0.0))
        return
    if L == 0.0 or K1 == 0.0:
        _tdrife0(st, L)
        return
    dx = getattr(el, 'DX', 0.0)
    dy = getattr(el, 'DY', 0.0)
    th = float(getattr(el, 'ROTATION', 0.0) or 0.0)   # SAD の ROT (tsolrote の chi3)
    # tquade.f: DISFRIN (tqfrie の有無)、F1/F2/FRINGE (tqlfre, tsetfringepe の長さ)
    fringe = float(getattr(el, 'DISFRIN', 0.0) or 0.0) == 0.0
    fr = _lat._quad_fringe_params(el)
    f1in, f2in, f1out, f2out, mfring = fr
    if not krad or (not np.any(st.cod[:4]) and dx == 0.0 and dy == 0.0 and th == 0.0):
        jq = lambda c: _jac_quad_full(c, L, K1, fringe, fr, ph, dx, dy, th)
        st.piece(_quad_map(L, K1, ph, dx, dy, th, fringe, fr), jac=jq)
        if krad:
            st.setr0()
        return
    # 軌道があるとき: SAD と同じ分割で放射を入れる
    ak1 = abs(K1) * L
    rot = K1 * L < 0.0

    # tsolrote: 入口で (x,y) -= (dx,dy) → ロール θ で回転 (→ K1<0 なら更に −π/2)、出口で逆
    ct, stt = math.cos(th), math.sin(th)

    def r_in(c):
        c = c.copy()
        c[0] -= dx
        c[2] -= dy
        if th != 0.0:
            x, px, y, py = c[0], c[1], c[2], c[3]
            c[0] = ct * x - stt * y; c[2] = stt * x + ct * y
            c[1] = ct * px - stt * py; c[3] = stt * px + ct * py
        if rot:
            x, px, y, py = c[0], c[1], c[2], c[3]
            c[0] = y; c[1] = py; c[2] = -x; c[3] = -px
        return c

    def r_out(c):
        c = c.copy()
        if rot:
            x, px, y, py = c[0], c[1], c[2], c[3]
            c[0] = -y; c[1] = -py; c[2] = x; c[3] = px
        if th != 0.0:
            x, px, y, py = c[0], c[1], c[2], c[3]
            c[0] = ct * x + stt * y; c[2] = -stt * x + ct * y
            c[1] = ct * px + stt * py; c[3] = -stt * px + ct * py
        c[0] += dx
        c[2] += dy
        return c

    # 解析ヤコビアン: ずれは微分に効かない。ロール θ と K1<0 の −π/2 回転は定数行列 (R_出口 · J(r_in(c)) · R_入口) で挟む
    Rin_t = np.eye(6)
    if th != 0.0:
        Rin_t, _ = _roll_matrices(th)
    if rot:
        Rin_t = _ROT_IN @ Rin_t
    Rout_t = Rin_t.T.copy()
    plain = (not rot) and th == 0.0

    def wrapj(Jc):
        if Jc is None:
            return None
        return Rout_t @ Jc @ Rin_t

    def jfun(core):
        """core(c_in) → 核の 6x6 (回転後の座標)。ずれ・ロール・−π/2 回転の定数行列で挟んで返す"""
        if plain:
            return core
        return lambda c: wrapj(core(r_in(c)))

    st.setr0()
    if fringe and mfring >= 0 and mfring != 2:
        st.piece(lambda c: r_out(_lat._tqfrie_cod(r_in(c), ak1, L)), jac=jfun(lambda c: _jac_tqfrie(c, ak1, L)))
    if mfring in (1, 3):
        st.piece(lambda c: r_out(_lat._tqlfre_cod(r_in(c), L, ak1, f1in, f2in)),
                 jac=jfun(lambda c: _jac_tqlfre(c, L, ak1, f1in, f2in)))
    if f1in != 0.0:
        st.tradke(f1in, 0.0)         # SAD のまま: 長さの引数に F1 (長さ) を渡す
    else:
        st.setr0()
    # itgetqraddiv
    c = r_in(st.cod)
    xk = ak1 * max(1e-6, abs(c[0]) + abs(c[2]))
    xpd = max(1e-6, abs(c[1]) + abs(c[3]))
    a = min(1e-2, xk + xpd)
    b = ph.brhoz * (xk / abs(L))
    nrad = int(abs(L * ph.crad / EPSRAD * (ph.h0 * b)**2))
    nq = max(int(ph.emidiv * ph.emidiq * nrad),
             int(abs(ph.h0 * b / ph.brhoz * ph.anrad) / EPSRAD / 1e6
                 * ph.emidiv * ph.emidib))
    nq = min(1000, max(1, nq))
    aln = L / nq
    akn = ak1 / nq
    for _ in range(nq):
        # tsolque(aln, akn, krad): 内部分割 + スライス毎 tradke
        ndiv = min(1000, 1 + int(abs(aln * akn) / 0.1))
        aln1 = aln / ndiv
        al1 = aln1 * 0.5
        for n in range(ndiv):
            def f_body(cc, al1=al1):
                cc = r_in(cc)
                pr = 1.0 + cc[5]
                cc = _lat._tqente_cod(cc, al1)
                cc = _lat._tsolque_body_cod(cc, aln1, akn / aln / pr)
                cc[4] -= ph.dv(cc[5]) * aln1
                return r_out(cc)

            def core_body(c, al1=al1):
                # tqente(al1) → tsolque 本体 (+ dv 項)。ヤコビアンは各段の入口の cod で作って掛ける
                J1 = _jac_tqente(c, al1)
                if J1 is None:
                    return None
                c1 = _lat._tqente_cod(c, al1)
                J2 = _jac_tsolque_body(c1, aln1, akn / aln, _dvz(ph, c1[5]))
                return None if J2 is None else J2 @ J1
            st.piece(f_body, jac=jfun(core_body))
            if n != ndiv - 1:
                st.tradke(aln1, 0.0)
            al1 = aln1
        st.piece(lambda cc: r_out(_lat._tqente_cod(r_in(cc), aln1 * 0.5)),
                 jac=jfun(lambda c: _jac_tqente(c, aln1 * 0.5)))
        st.tradke(aln1, 0.0)
    if mfring in (2, 3):
        st.piece(lambda c: r_out(_lat._tqlfre_cod(r_in(c), L, ak1, -f1out, f2out)),
                 jac=jfun(lambda c: _jac_tqlfre(c, L, ak1, -f1out, f2out)))
    if fringe and mfring >= 0 and mfring != 1:
        st.piece(lambda c: r_out(_lat._tqfrie_cod(r_in(c), -ak1, L)), jac=jfun(lambda c: _jac_tqfrie(c, -ak1, L)))
    if f1out != 0.0:
        st.tradke(f1out, 0.0)


# --- K1≠0 BEND: tbendeinit/tbendiinit, tbendecorr, tbendebody (tbende.f, tbend.f)

def _sxsin(x):
    s = math.sin(x)
    return s, x - s


def _sxsinh(x):
    sh = math.sinh(x)
    return sh, x - sh


class _BendIBody:
    """bendib + bendeb モジュール変数 (tbendiinit / tbendeinit)"""

    def __init__(self, ak1, al, dp, rhob, rho0):
        self.rhob = rhob
        self.rho0 = rho0
        p = 1.0 + dp
        self.p = p
        rhoe = rhob * p
        self.rhoe = rhoe
        # tbendiinit
        rhosq = rho0 * rhoe
        drhop = (rhoe - rho0) / rhosq
        akk = ak1 / al
        akxsq = -1.0 / rhosq - akk / p
        akysq = akk / p
        if akxsq > 0.0:
            akx = math.sqrt(akxsq)
            phix = akx * al
            dcx = 2.0 * math.sinh(0.5 * phix)**2
            spx, xspx = _sxsinh(phix)
            aksx = akx * spx
            dcxkx = dcx / akxsq
            sxkx = spx / akx
            xsxkx = xspx / akx / akxsq
        elif akxsq < 0.0:
            akx = math.sqrt(-akxsq)
            phix = akx * al
            dcx = -2.0 * math.sin(0.5 * phix)**2
            spx, xspx = _sxsin(phix)
            aksx = -akx * spx
            dcxkx = dcx / akxsq
            sxkx = spx / akx
            xsxkx = -xspx / akx / akxsq
        else:
            akx = 0.0; phix = 0.0; dcx = 0.0; spx = 0.0; aksx = 0.0
            dcxkx = 0.5 * al**2
            sxkx = al
            xsxkx = -1.0 / 6.0 * al**3
        if akysq > 0.0:
            aky = math.sqrt(akysq)
            phiy = aky * al
            dcy = 2.0 * math.sinh(0.5 * phiy)**2
            spy = math.sinh(phiy)
            aksy = aky * spy
            syky = spy / aky
        elif akysq < 0.0:
            aky = math.sqrt(-akysq)
            phiy = aky * al
            dcy = -2.0 * math.sin(0.5 * phiy)**2
            spy = math.sin(phiy)
            aksy = -aky * spy
            syky = spy / aky
        else:
            aky = 0.0; phiy = 0.0; dcy = 0.0; spy = 0.0; aksy = 0.0
            syky = al
        # tbendeinit
        cx = 1.0 + dcx
        cy = 1.0 + dcy
        drhopp = 1.0 / rhoe / p
        akxsqp = -akxsq / p
        phixp = -0.5 * phix / p
        if akxsq > 0.0:
            dcxp = spx * phixp
            dcxkxp = dcxp / akxsq + dcxkx / p
            sxkxp = 0.5 * sxkx / p + cx * phixp / akx
            xsxkxp = (-dcx * phixp / akxsq - 1.5 * xsxkx / p) / akx
            aksxp = -0.5 * aksx / p + akx * cx * phixp
        elif akxsq < 0.0:
            dcxp = -spx * phixp
            dcxkxp = dcxp / akxsq + dcxkx / p
            sxkxp = 0.5 * sxkx / p + cx * phixp / akx
            xsxkxp = (-dcx * phixp / akxsq - 1.5 * xsxkx / p) / akx
            aksxp = -0.5 * aksx / p - akx * cx * phixp
        else:
            dcxp = 0.0; dcxkxp = 0.0; sxkxp = 0.0; xsxkxp = 0.0; aksxp = 0.0
        akysqp = -akysq / p
        phiyp = -0.5 * phiy / p
        if akysq > 0.0:
            dcyp = phiyp * spy
            sykyp = 0.5 * syky / p + cy * phiyp / aky
            aksyp = -0.5 * aksy / p + aky * cy * phiyp
        elif akysq < 0.0:
            dcyp = -phiyp * spy
            sykyp = 0.5 * syky / p + cy * phiyp / aky
            aksyp = -0.5 * aksy / p - aky * cy * phiyp
        else:
            dcyp = 0.0; sykyp = 0.0
            aksyp = 0.0   # SAD では未設定 (akysq=0 は K1≠0 では起きない)
        self.__dict__.update(dict(
            drhop=drhop, akxsq=akxsq, akysq=akysq, dcx=dcx, aksx=aksx,
            dcxkx=dcxkx, sxkx=sxkx, xsxkx=xsxkx, dcy=dcy, aksy=aksy,
            syky=syky, cx=cx, cy=cy, drhopp=drhopp, akxsqp=akxsqp,
            dcxp=dcxp, dcxkxp=dcxkxp, sxkxp=sxkxp, xsxkxp=xsxkxp,
            aksxp=aksxp, akysqp=akysqp, dcyp=dcyp, sykyp=sykyp,
            aksyp=aksyp))


def _tbendecorr(st: _State, ak1, al, rho0):
    """tbende.f tbendecorr (解析行列ごと直訳。本体は perfkern_bend.bendecorr)"""
    t, cod = _k_bendecorr(st.cod, ak1, al, rho0, PSQMAX)
    st.apply(t)
    st.cod = cod


def _tbendebody(st: _State, al, phin, ak1, alr, rhob, rho0, enarad):
    """tbende.f tbendebody (tbinit=.true. で毎回 tbendeinit。本体は perfkern_bend.bendebody)"""
    t, cod = _k_bendebody(st.cod, al, phin, ak1, rhob, rho0, st.ph.h0, st.h1emit, st.dvemit, PSQMAX)
    st.apply(t)
    st.cod = cod
    if enarad:
        st.tradke(alr, phin)


def _tsteee(st: _State, el: Bend, enarad: bool = True):
    """tsteee.f (ANGLE=0, K1=0 の BEND = ステアリング)。K0>0 は px を減らす。
    L=0 → 薄い双極キック (tthine nord=2, 放射なし)、K0=0 → ドリフト。L>0 は入口エッジ (rhob_e = L/K0, psi=AE1)
    → 本体 (rhob = −L/K0) ×ndiv (放射ありはスライスごとに tradke) → 出口エッジ。F1 フリンジ (mfring≠0) は未対応。
    ずれ DX/DY は tchge (ずらすだけ。DROTATE は tsteee に渡らない)。"""
    ph = st.ph
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    al = el.L
    k0 = g('K0')
    dx, dy = g('DX'), g('DY')
    if al == 0.0:
        if k0 != 0.0:
            def kick(c):
                c = c.copy()
                c[1] -= k0
                return c
            st.piece(kick)
        return
    if k0 == 0.0:
        _tdrife0(st, al)
        return
    if g('F1') != 0.0 or g('FB1') != 0.0 or g('FB2') != 0.0:
        raise NotImplementedError('emit6d: ANGLE=0 の BEND (tsteee) の F1/FB1/FB2 フリンジは未対応')
    if dx != 0.0 or dy != 0.0:
        def shift_in(c):
            c = c.copy()
            c[0] -= dx
            c[2] -= dy
            return c
        st.piece(shift_in)
    krad = enarad
    if krad:
        st.setr0()
    rhoe = al / k0
    rhob = -rhoe
    ndiv = 1
    if krad:
        b = -ph.brhoz / rhob
        nrad = 1 + int(abs(al / EPSRAD * ph.crad * (ph.h0 * b)**2))
        ndiv = max(int(nrad * ph.emidiv * ph.emidib),
                   1 + int(abs(k0 * ph.h0 * ph.anrad) / EPSRAD / 1e6 * ph.emidiv * ph.emidib))
    aln = al / ndiv
    st.piece(lambda c: _lat._tbedge_cod(c, rhoe, g('E1'), True))
    for n in range(1, ndiv + 1):
        st.piece(lambda c: _lat._tsteee_body(c, aln, rhob, ph.dv(c[5]))[0])
        if n != ndiv and krad:
            st.tradke(aln, 0.0)
    st.piece(lambda c: _lat._tbedge_cod(c, rhoe, g('E2'), False))
    if krad:
        st.tradke(aln, 0.0)
    if dx != 0.0 or dy != 0.0:
        def shift_out(c):
            c = c.copy()
            c[0] += dx
            c[2] += dy
            return c
        st.piece(shift_out)


def _tbende(st: _State, el: Bend, enarad: bool = True):
    """tbende.f tbende。ずれ DX/DY・ロール ROTATE/ROTATION・DROTATE は SAD と同じく tchge (入口: 回す→ずらす→tbrote、出口はその逆) を
    外側に付け、中身 (_tbende_core) は誤差なしの BEND として動かす (放射の setr0/tradke は tchge の内側)。
    tchge の引数は lattice._bend_chg_args (tbende.f の分岐ごとの theta/dtheta の渡し方)。"""
    chg = _lat._bend_chg_args(el)
    if chg is None:
        _tbende_core(st, el, enarad)
        return
    dx, dy, dz, th, dth, alg, phig = chg
    core = _lat._bend_noerr(el)
    st.piece(lambda c: _lat._tchge_cod(c, dx, dy, dz, th, dth, alg, phig, True),
             jac=lambda c: _lat._tchge_jac(c, dx, dy, dz, th, dth, alg, phig, True))
    _tbende_core(st, core, enarad)
    al0, phi0 = (el.L, float(el.ANGLE)) if phig != 0.0 else (0.0, 0.0)
    st.piece(lambda c: _lat._tchge_cod(c, dx, dy, dz, th, dth, alg - al0, phig - phi0, False),
             jac=lambda c: _lat._tchge_jac(c, dx, dy, dz, th, dth, alg - al0, phig - phi0, False))


def _tbende_core(st: _State, el: Bend, enarad: bool = True):
    """tbende.f tbende (phib=phi0=ANGLE, F1=0, mfring=0, ずれ・ロールなし)"""
    ph = st.ph
    al0 = el.L
    phi0 = el.ANGLE
    phib = phi0 + float(getattr(el, 'K0', 0.0) or 0.0)   # SAD: phib = ANGLE + K0
    pe1 = getattr(el, 'E1', 0.0)    # pysad の E1 は rad (= psi1*phi0)
    pe2 = getattr(el, 'E2', 0.0)
    ak = getattr(el, 'K1', 0.0) * al0   # SAD の K1 (積分値)
    if phi0 == 0.0:
        # tbende.f の phi0==0 の分岐: K1=0 → tsteee、K1≠0 で K0=0 → tquade、K1≠0 で K0≠0 → tmulte (akm(0)=K0, akm(1)=K1·L)
        if ak == 0.0:
            _tsteee(st, el, enarad)
        elif phib == 0.0:
            _tquade(st, Quadrupole(name=el.name, L=al0, K1=ak / al0), enarad)
        else:
            _tmulte(st, _lat._bend_as_mult(el), enarad)
        return
    al = al0
    if enarad:
        st.setr0()
    phibl = phib / al
    rhob = 1.0 / phibl
    rho0 = al / phi0
    rbc = 1.0
    phic = phi0 * rbc
    eps = EPSBEND
    b = ph.brhoz / rhob
    ndiv = 1 if ak == 0.0 else 1 + int(abs(phic / eps))
    if enarad:
        epsr1 = EPSRAD
        nrad = 1 + int(abs(al0 * rbc / epsr1 * ph.crad * (ph.h0 * b)**2))
        ndiv = max(ndiv, int(nrad * ph.emidiv * ph.emidib),
                   1 + int(abs(phib * ph.h0 * ph.anrad) / epsr1 / 1e6
                           * ph.emidiv * ph.emidib))
    # 入口エッジ (解析ヤコビアン: tbfrie.f の trans を段ごとに掛けたもの)
    st.piece(lambda c: _lat._tbedge_cod(c, rhob, pe1, True),
             jac=lambda c: _lat._bend_edge_chain(c, _EYE6, rhob, pe1, True)[1])
    alc = al * rbc
    aln = alc / ndiv
    phin = phic / ndiv
    if ak == 0.0:
        for n in range(1, ndiv + 1):
            def f(c):
                c = _lat._tbendebody0_cod(c, phin, rhob, rhob - rho0)
                c = c.copy()
                c[4] -= ph.dv(c[5]) * aln
                return c

            def jbody(c, phin=phin, aln=aln):
                J = _lat._tbendebody0_mat(c, phin, rhob, rhob - rho0)
                J[4, 5] += _dvz(ph, c[5]) * aln          # z -= dv(δ)·aln の δ 微分 (SAD の h0/h1emit**3·aln)
                return J
            st.piece(f, jac=jbody)
            if enarad and n < ndiv:
                st.tradke(aln, phin)
        alr = aln
        phi1 = phin
    else:
        alx0 = 0.0
        akx0 = 0.0
        akx = 0.0
        alr = aln
        phi1 = phin
        for n in range(1, ndiv + 1):
            alx = aln
            alr = aln
            akx = ak * alx / al0
            phi1 = phi0 * alx / al0
            _tbendecorr(st, (akx + akx0) * 0.5, (alx + alx0) * 0.5, rho0)
            akx0 = akx
            alx0 = alx
            _tbendebody(st, alx, phi1, akx, alr, rhob, rho0,
                        enarad and n != ndiv)
        _tbendecorr(st, akx0 * 0.5, alx0 * 0.5, rho0)
    # 出口エッジ
    st.piece(lambda c: _lat._tbedge_cod(c, rhob, pe2, False),
             jac=lambda c: _lat._bend_edge_chain(c, _EYE6, rhob, pe2, False)[1])
    if enarad:
        st.tradke(alr, phi1)


def _h2p(h):
    return math.sqrt(max(0.0, (h - 1.0) * (h + 1.0)))


def _tcavfrie_map(al, v, w, dphis, s0, p0):
    """tcave.f tcavfrie の COD 写像 (autophi=F)。出口は v → −v で呼ぶ。"""
    vf = v / al / p0 * 0.5

    def f(cc):
        cc = cc.copy()
        dp1r = cc[5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        h1 = _p2h(p1)
        v1 = p1 / h1
        t = -cc[4] / v1
        ph = w * t - dphis
        sph, cph = math.sin(ph), math.cos(ph)
        dpt = vf * (sph + s0)
        cc[1] += cc[0] * dpt
        cc[3] += cc[2] * dpt
        wc = -w * vf
        dha = wc * (cc[0]**2 + cc[2]**2) * 0.5
        dh = dha * cph * p0
        h2 = h1 + dh
        a = dh * (h1 + h2) / p1**2
        dpr = a / (1.0 + math.sqrt(1.0 + a))
        dp2r = dp1r + p1r * dpr
        p2r = 1.0 + dp2r
        v2 = p2r * p0 / h2
        cc[5] = dp2r
        cc[4] = cc[4] * v2 / v1
        return cc
    return f


def _tcave(st: _State, el):
    """
    SAD tcave.f の直訳 (V1/V20/V11/V02 キーワード=0, autophi=F, 加速なし gammab 一定)。DX/DY/ROTATE は入口・出口の tchge。
    w = omega0*HARM/c (HARM≠0) または 2π FREQ/c。phis = trf0*w, phic = PHI。
    ndiv = 1+int(min(|w L|, sqrt((v(1/h0+1/h1))²/3/eps))), eps=1e-2。
    RFSW ON: L≠0 なら入口 tcavfrie → [ドリフト (最初は L/2ndiv) → キック]×ndiv → ドリフト L/2ndiv → 出口 tcavfrie。
      キック: dh = -va (sin φ + offset1), φ = w t + phic − phis, t = −z/β。
      offset1 = sin(phis) (RADCOD OFF) / 0 (RADCOD ON)。
    RFSW OFF: 積算量だけ足して長さ L のドリフト。
    行列は COD 写像の数値微分 (SAD は解析式 trans1)。空洞はエネルギーを変えるので
    z 行のシンプレクティック組み直しはしない。
    """
    ph = st.ph
    al = el.L
    vc = float(getattr(el, 'VOLT', 0.0) or 0.0)
    harm = float(getattr(el, 'HARM', 0.0) or getattr(el, '_harm', 0.0) or 0.0)
    freq = float(getattr(el, 'FREQ', 0.0) or 0.0)
    phi = float(getattr(el, 'PHI', 0.0) or 0.0)
    fringe = float(getattr(el, 'DISFRIN', 0.0) or 0.0) == 0.0
    mfring = int(round(float(getattr(el, 'FRINGE', 0.0) or 0.0)))
    cdx, cdy, cth = (float(getattr(el, k, 0.0) or 0.0) for k in ('DX', 'DY', 'ROTATION'))
    if cdx != 0.0 or cdy != 0.0 or cth != 0.0:
        # tcave.f: 入口 tchge(dx,dy,0,theta,0,0,0,0,.true.) (回してからずらす。th≠0 かつずれ≠0 の行列の癖 A25 は _tchge_jac)
        st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True),
                 jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True))
    c = CVELOC
    w = ph.omega0 * harm / c if harm != 0.0 else 2.0 * math.pi * freq / c
    wi = 1.0 / w if w != 0.0 else 0.0
    p0 = ph.p0
    h0 = ph.h0
    v = vc / ph.amass * abs(ph.charge)
    ndiv = 1 + int(min(abs(w * al), math.sqrt((v * (1.0 / h0 + 1.0 / h0))**2 / 3.0 / 1e-2)))
    vn = v / ndiv
    vcn = vc / ndiv
    aln = al / ndiv
    phis = st.trf0 * w
    phic = phi * (1.0 if ph.charge >= 0 else -1.0)
    v20a = vn * (w * (0.5 / p0 + 0.5 / p0))**2 / 4.0
    v02a = v20a
    st.vc0 += vc
    if ph.omega0 != 0.0:
        st.hvc0 += (c * w) / ph.omega0 * vc
    oneev = 1.0 + 3.83e-12
    if st.rfsw:
        if st.radcod:
            offset1 = 0.0
            s0 = 0.0
        else:
            s0 = math.sin(phis)
            offset1 = math.sin(phis)
        if al != 0.0 and fringe and mfring >= 0 and mfring != 2:
            st.piece(_tcavfrie_map(al, v, w, phis - phic, s0, p0), symp_z=False)

        def f(cc):
            cc = cc.copy()
            p1 = p0 * (1.0 + cc[5])
            h1 = _p2h(p1)
            v1 = p1 / h1
            t = -cc[4] / v1
            va = vn + v20a * cc[0]**2 + v02a * cc[2]**2
            phii = w * t + phic - phis
            sp, cp = math.sin(phii), math.cos(phii)
            dh = max(oneev - h1, -va * (sp + offset1))
            h2 = h1 + dh
            p2 = _h2p(h2)
            pf = (h2 + h1) / (p2 + p1) * dh
            p2 = p1 + pf
            v2 = p2 / h2
            av = -(cp * wi - offset1 * t) / p0
            cc[1] += 2.0 * v20a * cc[0] * av
            cc[3] += 2.0 * v02a * cc[2] * av
            cc[5] += pf / p0
            cc[4] = -t * v2
            return cc
        for n in range(ndiv):
            if al != 0.0:
                _tdrife0(st, aln * 0.5 if n == 0 else aln)    # tgetdvh(dgb=0) → dv=0
            # 積算量 (実際の cod で)
            p1 = p0 * (1.0 + st.cod[5])
            v1 = p1 / _p2h(p1)
            t = -st.cod[4] / v1
            phii = w * t + phic - phis
            sp, cp = math.sin(phii), math.cos(phii)
            st.dvcacc += vcn * cp * w
            st.ddvcacc += vcn * sp * w**2
            st.vcacc -= vcn * sp
            st.piece(f, symp_z=False)
        if al != 0.0:
            _tdrife0(st, aln * 0.5)
            if fringe and mfring >= 0 and mfring != 1:
                st.piece(_tcavfrie_map(al, -v, w, phis - phic, s0, p0), symp_z=False)
    else:
        sp, cp = math.sin(phic), math.cos(phic)
        st.dvcacc += vc * cp * w
        st.ddvcacc += vc * sp * w**2
        st.vcacc -= vc * sp
        if al != 0.0:
            _tdrife0(st, al)
    if cdx != 0.0 or cdy != 0.0 or cth != 0.0:
        st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False),
                 jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False))     # 出口 tchge


def _itgetqraddiv(ph, c, ak1, L):
    """tsolque.f itgetqraddiv 相当 (bz=0): 放射の分割数。c は要素座標の cod。_tquade と同じ式"""
    xk = abs(ak1) * max(1e-6, abs(c[0]) + abs(c[2]))
    b = ph.brhoz * (xk / abs(L))
    nrad = int(abs(L * ph.crad / EPSRAD * (ph.h0 * b)**2))
    nq = max(int(ph.emidiv * ph.emidiq * nrad),
             int(abs(ph.h0 * b / ph.brhoz * ph.anrad) / EPSRAD / 1e6
                 * ph.emidiv * ph.emidib))
    return min(1000, max(1, nq))


def _tmulte(st: _State, el, enarad: bool = True):
    """
    SAD tmulte.f の直訳 (ANGLE=0, BZ=0, SK なし, RF なし, F1 なし)。K0/K1/K2… とずれ・回転に対応。
    写像の列は lattice._mult_sequence (tsolrote → tqfrie/tbfrme → [tsolque → tmulke]×ndiv → tsolque →
    tbfrme/tqfrie → 戻し)。放射あり: tsolrote の後に tsetr0、tsolque の部分スライスごとに tradke。
    ndiv (放射あり) = min(300, max(ndiv_kick, itgetqraddiv(cod, K1, L)))。
    追跡側 (tmulti) と違い、六極などのフリンジ ttfrin は入らず、キックは (1+δ) で割らない。
    """
    al = el.L
    if any(float(getattr(el, k, 0.0) or 0.0) for k in ('ANGLE', 'F1', 'FRINGE', 'VOLT')):
        raise NotImplementedError('emit6d: MULT の曲げ・F1・RF は未対応')
    pm = _lat._mult_params(el)
    krad = enarad and al != 0.0
    if pm['nmmax'] < 0:
        _tdrife0(st, al)
        return
    nq = None
    if krad:
        ct, stt = math.cos(pm['theta2']), math.sin(pm['theta2'])
        c = _lat._mult_rot_in(st.cod, pm['dx'], pm['dy'], ct, stt)
        nq = _itgetqraddiv(st.ph, c, pm['knl'].get(1, 0.0), al)
        st.setr0()
    seq = _lat._mult_sequence(el, st.ph.dv, nq)
    first = True
    for kind, f in seq:
        if kind == 'map':
            st.piece(f)
            if first and krad:
                st.setr0()      # tmulte: フリンジの後 (f1in=0) に tsetr0
            first = False
        elif krad:
            st.tradke(f, 0.0)

def _rf_effective(st: _State):
    """tturneg: wrfeff, vceff (実効角周波数 [1/m] と実効電圧)。"""
    if st.vcacc != 0.0:
        wrfeff = math.sqrt(abs(st.ddvcacc / st.vcacc))
    elif st.vc0 != 0.0:
        wrfeff = abs(st.dvcacc / st.vc0)
    else:
        wrfeff = 0.0
    if wrfeff == 0.0 and st.vc0 != 0.0:
        wrfeff = st.hvc0 / st.vc0 * st.ph.omega0 / CVELOC
    vceff = math.hypot(st.vcacc, st.dvcacc / wrfeff) if wrfeff != 0.0 else 0.0
    if vceff == 0.0:
        vceff = st.vc0
    return wrfeff, vceff


def _trf0_update(st: _State, wrfeff, vceff):
    """tturne (update=.true., RADCOD OFF): 同期位置 trf0。"""
    u0e = st.u0 * st.ph.pgev
    alambdarf = 2.0 * math.pi / wrfeff
    if vceff > u0e:
        a = math.asin(max(-1.0, min(1.0, u0e / vceff)))
        trf0 = (a / wrfeff) if st.T[4, 5] < 0.0 else ((math.pi - a) / wrfeff)
    else:
        trf0 = (0.5 * math.pi) / wrfeff
    # Fortran の mod は被除数の符号 (Python の % とは違う) → math.fmod
    if trf0 < 0.0:
        trf0 = -math.fmod(-trf0 + 0.5 * alambdarf, alambdarf)
    else:
        trf0 = math.fmod(trf0 - 0.5 * alambdarf, alambdarf)
    return trf0 + 0.5 * alambdarf


def _tturne_update_trf0(st: _State, trf0):
    """tturne.f の update=.true. (RADCOD OFF): tturneg の後で trf0 を更新する。
    wrfeff≠0 なら vceff=hypot(vcacc, dvcacc/wrfeff) (vc0 へのフォールバックなし)、vceff≠0 のときだけ trf0 を作り直す
    (trans(5,6)<0 かで asin / π-asin)。vceff=0 (空洞なし等) なら trf0 はそのまま"""
    wrf, _ = _rf_effective(st)
    if wrf == 0.0:
        return trf0
    vce = math.hypot(st.vcacc, st.dvcacc / wrf)
    if vce == 0.0:
        return trf0
    return _trf0_update(st, wrf, vce)


def _tthine(st: _State, el):
    """SEXT/OCT/DECA/DODECA: tturne1 → tthine (enarad=.false., tturne.f:592)。放射なし。
    ずれ DX/DY・ロール ROTATE(ROTATION) は tthine の先頭と末尾の tchge (入口: 回す→ずらす、出口: 逆。dz=dtheta=0, alg=phig=0)。
    tchge の入口は SAD の癖 (th≠0 かつずれ≠0 のとき行列に回転が入らない) を _tchge_jac が再現する。"""
    kn, kord = _lat._THIN_KORD[type(el).__name__]
    K = float(getattr(el, kn, 0.0) or 0.0)
    ak = K * el.L if el.L != 0.0 else K
    _tthine_gen(st, el.L, ak, kord, float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0),
                float(getattr(el, 'ROTATION', 0.0) or 0.0))


def _tthine_gen(st: _State, L, ak, kord, dx, dy, th):
    """tthine の本体 (ak = 積分強度, kord = nord/2−1。薄い四極 (L=0 の QUAD, tquade → tthine(nord=4)) も kord=1 で使う)"""
    ph = st.ph
    f0 = _lat._thin_offset_map(L, ak, kord)

    def f_core(c):
        # tthine の中のドリフトは tdrife0(dvon=.true.) なので dv 項 (z -= L·dv(δ)) が付く。
        # 薄いキックは δ を変えないので全長まとめて引ける (Level 16: ESRF の R56 が 5.6e-6 ずれていた原因)
        c2 = f0(c)
        c2[4] -= ph.dv(c2[5]) * L
        return c2

    jcore = lambda c: _jac_thin(c, L, ak, kord, _dvz(ph, c[5]))
    if ak == 0.0 or (dx == 0.0 and dy == 0.0 and th == 0.0):
        # ak=0 は tdrife0 だけ (tchge を呼ばない)
        st.piece(f_core, jac=jcore)
        return
    st.piece(lambda c: _lat._tchge_cod(c, dx, dy, 0.0, th, 0.0, 0.0, 0.0, True),
             jac=lambda c: _lat._tchge_jac(c, dx, dy, 0.0, th, 0.0, 0.0, 0.0, True))
    st.piece(f_core, jac=jcore)
    st.piece(lambda c: _lat._tchge_cod(c, dx, dy, 0.0, th, 0.0, 0.0, 0.0, False),
             jac=lambda c: _lat._tchge_jac(c, dx, dy, 0.0, th, 0.0, 0.0, 0.0, False))


def _g0(el, k):
    return float(getattr(el, k, 0.0) or 0.0)


def _bende_fast(st, el, enarad):
    """_tbende を elkern.k_bend で。戻り 0: 成功、それ以外: Python 版で計算し直す (tmulte、tsteee の F1 など)"""
    ang = float(el.ANGLE)
    k0 = _g0(el, 'K0')
    k1 = _g0(el, 'K1')
    if ang == 0.0 and k1 == 0.0 and (_g0(el, 'F1') != 0.0 or _g0(el, 'FB1') != 0.0 or _g0(el, 'FB2') != 0.0):
        return 1
    chg = _lat._bend_chg_args(el)
    if chg is None:
        has_chg = False
        dx = dy = dz = th = dth = alg = phig = 0.0
    else:
        has_chg = True
        dx, dy, dz, th, dth, alg, phig = chg
    return _ek.k_bend(st.cod, st.T, st.D, st.B, st.codr0, st.transr, st._S, st.PH, float(el.L), ang, k0, k1,
                      _g0(el, 'E1'), _g0(el, 'E2'), has_chg, dx, dy, dz, th, dth, alg, phig,
                      enarad, True, True, st.want_jac, st.want_db, st.radcod)


def _quade_fast(st, el, krad):
    """_tquade を elkern.k_quad (L=0 は k_thin、K1=0 は k_drift) で"""
    L = float(el.L)
    K1 = _g0(el, 'K1')
    if L == 0.0 and K1 != 0.0:
        return _ek.k_thin(st.cod, st.T, st.D, st.B, st._S, st.PH, 0.0, K1, 1, _g0(el, 'DX'), _g0(el, 'DY'),
                          _g0(el, 'ROTATION'), True, True, st.want_jac, st.want_db)
    if L == 0.0 or K1 == 0.0:
        return _ek.k_drift(st.cod, st.T, st.D, st.B, st._S, st.PH, L, True, True, st.want_jac, st.want_db)
    f1in, f2in, f1out, f2out, mfring = _lat._quad_fringe_params(el)
    return _ek.k_quad(st.cod, st.T, st.D, st.B, st.codr0, st.transr, st._S, st.PH, L, K1,
                      _g0(el, 'DX'), _g0(el, 'DY'), _g0(el, 'ROTATION'), _g0(el, 'DISFRIN') == 0.0,
                      f1in, f2in, f1out, f2out, int(mfring), krad, True, True, st.want_jac, st.want_db, st.radcod)


def _thine_fast(st, el):
    """_tthine を elkern.k_thin で"""
    kn, kord = _lat._THIN_KORD[type(el).__name__]
    K = _g0(el, kn)
    L = float(el.L)
    ak = K * L if L != 0.0 else K
    return _ek.k_thin(st.cod, st.T, st.D, st.B, st._S, st.PH, L, ak, int(kord), _g0(el, 'DX'), _g0(el, 'DY'),
                      _g0(el, 'ROTATION'), True, True, st.want_jac, st.want_db)


def _nonfinite_state(st) -> bool:
    """cod が有限でない (tcod の反復が発散して NaN/Inf になった) とき True にして、残りの要素を NaN で埋める。
    Python 版の写像に NaN を通すと要素ごとに数値微分 (QUAD 1 個 0.09 s) を繰り返すだけで結果はどのみち NaN"""
    if np.all(np.isfinite(st.cod)):
        return False
    st.cod[:] = np.nan
    st.T[:, :] = np.nan
    st.D[:, :] = np.nan
    st.B[:, :] = np.nan
    return True


def tturne(elements: List[Element], cod0, ph: _Phys, trf0: float = 0.0,
           rfsw: bool = False, radcod: bool = False, want_db: bool = True, want_jac: bool = True) -> _State:
    """tturne1 相当: 放射込み 6D 一周追跡。want_db=False は SAD の irad=6 (tcod の反復用: D・B 行列を作らない)"""
    st = _State(cod0, ph)
    st.want_jac = want_jac
    st.want_db = want_db and want_jac
    st.trf0 = trf0
    st.rfsw = rfsw
    st.radcod = radcod
    st.setr0()
    fast = _USE_ELKERN
    for el in elements:
        if isinstance(el, Drift):
            if not fast or _ek.k_drift(st.cod, st.T, st.D, st.B, st._S, st.PH, float(el.L), True, True,
                                       st.want_jac, st.want_db) != 0:
                if fast and _nonfinite_state(st):
                    break
                _tdrife0(st, el.L)
        elif isinstance(el, Bend):
            if not fast or _bende_fast(st, el, True) != 0:
                if fast and _nonfinite_state(st):
                    break
                _tbende(st, el, enarad=True)
        elif isinstance(el, Quadrupole):
            if not fast or _quade_fast(st, el, True) != 0:
                if fast and _nonfinite_state(st):
                    break
                _tquade(st, el, krad=True)
        elif type(el).__name__ in _lat._THIN_KORD:
            if not fast or _thine_fast(st, el) != 0:
                if fast and _nonfinite_state(st):
                    break
                _tthine(st, el)
        elif type(el).__name__ == 'Cavity':
            _tcave(st, el)
        elif type(el).__name__ == 'Multipole':
            _tmulte(st, el, enarad=True)
        elif el.L == 0.0:
            continue          # MARK など
        else:
            raise NotImplementedError(
                f"emit6d: 要素 {type(el).__name__} は未対応")
    return st


# ---------------------------------------------------------------------------
# tecalc: 固有値分解と平衡ビーム行列
# ---------------------------------------------------------------------------

def _tnorm(r, ceig, warn=None):
    """tnorm.f (toplvl.f) の直訳。teigenc の固有ベクトル r (列 i,i+1 = 実部,虚部) を正規化し、モードを並べ替える。
    実固有値の対 (不安定) もそのまま扱う。r, ceig を in-place で更新。
    警告 (s==0) のとき warn() を呼ぶ (SAD: '???-tnorm-Unstable transfer matrix.')"""
    for i in (0, 2, 4):
        if ceig[i].imag == 0.0:
            sa = float(np.sum(r[:, i] ** 2))
            sb = float(np.sum(r[:, i + 1] ** 2))
            a = math.sqrt(math.sqrt(sa / sb))
            r[:, i] = r[:, i] / a
            r[:, i + 1] = r[:, i + 1] * a
        s = 0.0
        for j in (0, 2, 4):
            s = s + r[j, i] * r[j + 1, i + 1] - r[j, i + 1] * r[j + 1, i]
        sa = math.sqrt(abs(s))
        if s > 0.0:
            sb = sa
        elif s == 0.0:
            if warn is not None:
                warn()
            sa = 1.0
            sb = 1.0
        else:
            sb = -sa
            ceig[i] = np.conj(ceig[i])
            ceig[i + 1] = np.conj(ceig[i + 1])
        r[:, i] = r[:, i] / sa
        r[:, i + 1] = r[:, i + 1] / sb
    s1 = (r[4, 0] * r[5, 1] - r[4, 1] * r[5, 0]) ** 2
    s2 = (r[4, 2] * r[5, 3] - r[4, 3] * r[5, 2]) ** 2
    s3 = (r[4, 4] * r[5, 5] - r[4, 5] * r[5, 4]) ** 2
    j = int(np.argmax([s1, s2, s3])) * 2           # maxloc は最初の最大
    if j != 4:
        v = r[:, 4:6].copy()
        r[:, 4:6] = r[:, j:j + 2]
        r[:, j:j + 2] = v
        cc = ceig[4]; ceig[4] = ceig[j]; ceig[j] = cc
        cc = ceig[5]; ceig[5] = ceig[j + 1]; ceig[j + 1] = cc
    s1 = (r[0, 0] * r[1, 1] - r[0, 1] * r[1, 0]) ** 2
    s2 = (r[0, 2] * r[1, 3] - r[0, 3] * r[1, 2]) ** 2
    if s2 > s1:
        v = r[:, 0:2].copy()
        r[:, 0:2] = r[:, 2:4]
        r[:, 2:4] = v
        cc = ceig[0]; ceig[0] = ceig[2]; ceig[2] = cc
        cc = ceig[1]; ceig[1] = ceig[3]; ceig[3] = cc
    for i in (0, 2, 4):
        if ceig[i].imag != 0.0:
            a = math.hypot(r[i, i], r[i, i + 1])
            cost = r[i, i] / a
            sint = r[i, i + 1] / a
            v1 = r[:, i].copy()
            r[:, i] = v1 * cost + r[:, i + 1] * sint
            r[:, i + 1] = -v1 * sint + r[:, i + 1] * cost
        else:
            a = r[i, i]
            if a != 0.0:
                b = r[i, i + 1]
                r[:, i + 1] = r[:, i + 1] * a - b * r[:, i]
                r[:, i] = r[:, i] / a
    return r, ceig


def _normal_form(rx, warn=None):
    """teigenc + tnorm の直訳 (teigen.py と _tnorm)。RF なしの z の Jordan 形 (λ=1 の重根) や
    不安定 (実固有値の対) も SAD と同じ並び・規格化になる。"""
    import os
    if os.environ.get('PYSAD_OLD_NF') == '1':       # 回帰比較用 (旧: numpy.linalg.eig 版)
        return _normal_form_old(rx)
    from .teigen import teigen
    ceig, r = teigen(rx)
    ceig = ceig.astype(complex)
    _order_real_pairs(rx, r, ceig)
    return _tnorm(r, ceig, warn)


def _order_real_pairs(rx, r, ceig):
    """不安定な面 (実固有値の対 λ, 1/λ が teigen の並びで隣り合う) の並び順。
    teigen は 2×2 ブロックの対角の差 (a-d) の符号で λ と 1/λ のどちらを先にするかを決める。始点が鏡映対称 (α=0) のリングでは
    a-d は厳密には 0 で、SAD でも pysad でも丸め雑音の符号で決まる (再現できない)。SAD の実測 (bench_level28c/28d) は
    増大する方 (λ>1) が先だったので、a-d が雑音の大きさ (相対 1e-8 未満) のときだけ増大側を先にする。
    a-d がはっきりあるときは teigen の並びのまま (SAD と同じ)。ceig, r を in-place で更新"""
    for i in (0, 2, 4):
        a, b = ceig[i], ceig[i + 1]
        if a.imag != 0.0 or b.imag != 0.0:
            continue
        if not (a.real < b.real and b.real > 1.0 and abs(a.real * b.real - 1.0) < 1e-6):
            continue
        u = r[:, i]
        wx = abs(u[0]) + abs(u[1])
        wy = abs(u[2]) + abs(u[3])
        if max(wx, wy) == 0.0 or wx == wy:
            continue
        k = 0 if wx > wy else 2
        da = rx[k, k] - rx[k + 1, k + 1]
        if abs(da) <= 1e-8 * max(abs(rx[k, k]), abs(rx[k + 1, k + 1])):
            r[:, [i, i + 1]] = r[:, [i + 1, i]]
            ceig[i], ceig[i + 1] = b, a


def _normal_form_old(rx):
    """teigenc + tnorm 相当。
    複素固有値のモードは tnorm と同じ正規化 (シンプレクティック積=1、
    r(i,i+1)=0, r(i,i)>0)。RF なしで z の固有値が 1 (Jordan) の場合は
    tnorm の実固有値分岐の結果と同じ形 r(:,5)=e5, r(:,6)=(分散,0,1) を直接作る。"""
    w, V = np.linalg.eig(rx)
    r = np.zeros((6, 6))
    ceig = np.zeros(6, dtype=complex)
    used = np.zeros(6, dtype=bool)
    cols = []
    for k in range(6):
        if used[k] or abs(w[k].imag) < 1e-12:
            continue
        # 共役のペアを探す
        kc = [j for j in range(6) if not used[j] and j != k
              and abs(w[j] - np.conj(w[k])) < 1e-9]
        used[k] = True
        if kc:
            used[kc[0]] = True
        cols.append((w[k], V[:, k]))
    # 振動モードの並び: x 成分が大きいものを 1 番目 (tnorm の並べ替え)
    def xweight(v):
        a = np.real(v); b = np.imag(v)
        return (a[0] * b[1] - a[1] * b[0])**2
    def zweight(v):
        a = np.real(v); b = np.imag(v)
        return (a[4] * b[5] - a[5] * b[4])**2
    if len(cols) >= 3:
        kz = max(range(len(cols)), key=lambda k: zweight(cols[k][1]))
        zc = cols.pop(kz)
        cols.sort(key=lambda t: -xweight(t[1]))
        cols.append(zc)
    else:
        cols.sort(key=lambda t: -xweight(t[1]))
    for m, (lam, v) in enumerate(cols[:2]):
        i = 2 * m
        a = np.real(v).copy()
        b = np.imag(v).copy()
        s = sum(a[j] * b[j + 1] - b[j] * a[j + 1] for j in (0, 2, 4))
        sa = math.sqrt(abs(s))
        if s > 0:
            sb = sa
            ceig[i] = lam
        else:
            sb = -sa
            ceig[i] = np.conj(lam)
        ceig[i + 1] = np.conj(ceig[i])
        a /= sa
        b /= sb
        hyp = math.hypot(a[i], b[i])
        cost = a[i] / hyp
        sint = b[i] / hyp
        r[:, i] = a * cost + b * sint
        r[:, i + 1] = -a * sint + b * cost
    # z モード
    if len(cols) >= 3:
        lam, v = cols[2]
        a = np.real(v).copy(); b = np.imag(v).copy()
        s = sum(a[j] * b[j + 1] - b[j] * a[j + 1] for j in (0, 2, 4))
        sa = math.sqrt(abs(s))
        if s > 0:
            sb = sa; ceig[4] = lam
        else:
            sb = -sa; ceig[4] = np.conj(lam)
        ceig[5] = np.conj(ceig[4])
        a /= sa; b /= sb
        hyp = math.hypot(a[4], b[4])
        cost = a[4] / hyp; sint = b[4] / hyp
        r[:, 4] = a * cost + b * sint
        r[:, 5] = -a * sint + b * cost
    else:
        M4 = rx[:4, :4]
        eta = np.linalg.solve(np.eye(4) - M4, rx[:4, 5])
        r[:, 4] = [0, 0, 0, 0, 1.0, 0]
        r[:, 5] = [eta[0], eta[1], eta[2], eta[3], 0.0, 1.0]
        ceig[4] = 1.0
        ceig[5] = 1.0
    return r, ceig


def _tsvd_tsolva(A, b, eps=1e-8):
    from .tsvd import tsvd
    return tsvd(A, b, eps)


def _tsolva(A, b, eps=1e-8):
    """tsvd: 特異値 < max*eps を捨てた SVD 解"""
    U, s, Vt = np.linalg.svd(A)
    thr = s.max() * eps
    sinv = np.where(s > thr, 1.0 / s, 0.0)
    return Vt.T @ (sinv * (U.T @ b))


def tecalc(st: _State, circ: float, rfsw: bool = False, emit_floor=(0.0, 0.0, 0.0)) -> dict:
    ph = st.ph
    T = st.T
    rx = T.copy()
    if not rfsw:
        rx[5, :5] = 0.0
        rx[5, 5] = 1.0
    r, ceig = _normal_form(rx)
    r = _tsymp(r)
    ri = _tinv6(r)
    cc = np.zeros(6, dtype=complex)
    cc[0::2] = ceig[0::2]
    cc[1::2] = np.conj(cc[0::2])
    cd = np.zeros(6, dtype=complex)
    with np.errstate(divide='ignore', invalid='ignore'):      # 不安定で固有値が 0 のときは NaN のまま (警告は出さない)
        cd[3:6] = np.log(cc[0::2])
    # tturneg は RFSW の ON/OFF に関係なく空洞の積算量 (vc0, vcacc, dvcacc …) から実効電圧 vceff・角周波数 wrfeff を作る。
    # RFSW OFF でも Vc・h・バケット高さ・縦方向エミッタンス・バンチ長が SAD では 0 にならない (以前は 0 を返していた)
    wrfeff, vceff = _rf_effective(st)
    u0e = st.u0 * ph.pgev
    if vceff != 0.0:
        phirf = math.asin(max(-1.0, min(1.0, u0e / vceff)))
        heff = wrfeff * CVELOC / ph.omega0
    else:
        phirf = 0.0
        heff = 0.0
    synchm = rfsw and cd[5].imag != 0.0
    if synchm:
        if wrfeff != 0.0:
            alphap = (-cd[5].imag * abs(cd[5].imag)
                      / (CVELOC * 2 * math.pi / ph.omega0) / (st.dvcacc / ph.pgev))
        else:
            alphap = 0.0
        omegaz = abs(cd[5].imag) * ph.omega0 / (2 * math.pi)
    else:
        alphap = -rx[4, 5] / (2 * math.pi) / CVELOC / ph.p0 * ph.h0 * ph.omega0
        # temit.f tecalc: omegaz=sqrt(|alphap 2π heff vceff/pgev cos(phirf)|) ω0/2π (RF の効果を RFSW OFF でも見積もる)
        omegaz = (math.sqrt(abs(alphap * 2.0 * math.pi * heff * vceff / ph.pgev * math.cos(phirf)))
                  * ph.omega0 / (2.0 * math.pi)) if vceff != 0.0 else 0.0
    if vceff != 0.0 and alphap != 0.0 and heff != 0.0:
        bh = math.sqrt(abs(vceff / math.pi / abs(alphap) / heff / ph.pgev
                           * (2.0 * math.cos(phirf) - (math.pi - 2.0 * phirf) * u0e / vceff)))
    elif vceff != 0.0:
        bh = float('inf')           # α_c=0 (等時性): SAD の Fortran は Inf
    else:
        bh = 0.0
    params = {}
    params['tunes'] = cd[3:6].imag / (2 * math.pi)
    rd = ri @ st.D @ r
    with np.errstate(divide='ignore', invalid='ignore'):
        for k, i in enumerate((0, 2, 4)):
            cd[k] = complex((rd[i, i] + rd[i + 1, i + 1]) * 0.5,
                            (rd[i, i + 1] - rd[i + 1, i]) * 0.5) / cc[i]
    params['damping'] = cd[0:3].real
    beam = _pack(ri @ st.B @ ri.T)
    rdb = np.zeros((6, 12))
    rdb[:, :6] = rd
    if not synchm:
        for jj in range(6):
            beam[_iaidx(4, jj)] = 0.0
        rdb[0:6, 4] = 0.0
        rdb[4, 0:6] = 0.0
    for i in (0, 2, 4):
        tune = cd[i // 2 + 3].imag
        rdb[i, i + 6] = math.cos(tune)
        rdb[i, i + 7] = math.sin(tune)
        rdb[i + 1, i + 6] = -math.sin(tune)
        rdb[i + 1, i + 7] = math.cos(tune)
    btr = np.zeros((21, 21))
    for i in range(6):
        for j in range(i + 1):
            k = _iaidx(i, j)
            for mm in range(6):
                for n in range(6):
                    l = _iaidx(mm, n)
                    btr[k, l] -= ((rdb[i, mm] + rdb[i, mm + 6])
                                  * (rdb[j, n] + rdb[j, n + 6]))
    sqr2 = math.sqrt(0.5)
    for i in (0, 2, 4):
        k1 = _iaidx(i, i)
        k2 = _iaidx(i + 1, i + 1)
        k3 = _iaidx(i, i + 1)
        bbv = btr[k1, :].copy()
        btr[k1, :] = (bbv + btr[k2, :]) * sqr2
        btr[k2, :] = (-bbv + btr[k2, :]) * sqr2
        btr[k3, :] = btr[k3, :] / sqr2
        bb = beam[k1]
        beam[k1] = (bb + beam[k2]) * sqr2
        beam[k2] = (-bb + beam[k2]) * sqr2
        beam[k3] = beam[k3] / sqr2
        bbv = btr[:, k1].copy()
        btr[:, k1] = (bbv + btr[:, k2]) * sqr2
        btr[:, k2] = (-bbv + btr[:, k2]) * sqr2
        btr[:, k3] = btr[:, k3] * sqr2
    for i in range(21):
        btr[i, i] += 1.0
    for i in (0, 2, 4):
        k = _iaidx(i, i)
        btr[k, k] = (-(rdb[i, i]**2 + rdb[i, i + 1]**2
                       + rdb[i + 1, i]**2 + rdb[i + 1, i + 1]**2) * 0.5
                     - rdb[i, i + 6] * (rdb[i, i] + rdb[i + 1, i + 1])
                     - rdb[i, i + 7] * (rdb[i, i + 1] - rdb[i + 1, i]))
    if not synchm:
        btr[14, 20] = 0.0
        btr[14, 14] *= 2.0
        btr[20, 14] = -btr[20, 20]
        beam[20] = 0.0
    for i in (0, 2, 4):
        k1 = _iaidx(i, i)
        k2 = _iaidx(i + 1, i + 1)
        if btr[k2, k2] != 0.0 and btr[k1, k1] != 0.0:
            ab = math.sqrt(abs(btr[k1, k1] / btr[k2, k2]))
            btr[k1, :] /= ab
            btr[k2, :] *= ab
            beam[k1] /= ab
            beam[k2] *= ab
    # tsolva (tsolvg.f) = tsvd: numpy の SVD (特異値の硬い打ち切り) ではなく SAD の 2 重対角化 + 暗黙 QR の直訳 (pysad/tsvd.py)。
    # y 方向が厳密に 0 のとき numpy 版は εy ~1e-21 の雑音を作る (SAD は ~1e-38)。他の量は 1e-13 以内で同じ
    emitn = _tsvd_tsolva(btr, beam, 1e-8)
    for i in (0, 2, 4):
        k1 = _iaidx(i, i)
        k2 = _iaidx(i + 1, i + 1)
        bb = emitn[k1]
        emitn[k1] = (bb - emitn[k2]) * sqr2
        emitn[k2] = (bb + emitn[k2]) * sqr2
        emitn[_iaidx(i, i + 1)] *= sqr2
    # temit.f:456-467: EMITXE/EMITYE/EMITZE (正規モードのエミッタンスの下限。符号は元のまま)
    for i, fl in zip((0, 2, 4), emit_floor):
        if fl:
            for k in (_iaidx(i, i), _iaidx(i + 1, i + 1)):
                emitn[k] = math.copysign(max(abs(emitn[k]), float(fl)), emitn[k])
    ch = ph.charge

    def _sgn(a, b):
        return abs(a) if b >= 0 else -abs(a)
    e = emitn
    eemx = _sgn(math.sqrt(abs(e[_iaidx(0, 0)] * e[_iaidx(1, 1)]
                              - e[_iaidx(0, 1)]**2)), e[_iaidx(1, 1)] * ch)
    eemy = _sgn(math.sqrt(abs(e[_iaidx(2, 2)] * e[_iaidx(3, 3)]
                              - e[_iaidx(2, 3)]**2)), e[_iaidx(3, 3)] * ch)
    emitp = r @ _unpack(emitn) @ r.T
    sige = math.sqrt(abs(emitp[5, 5]))
    if synchm:
        sigz = math.sqrt(abs(emitp[4, 4]))
        eemz = _sgn(math.sqrt(abs(e[_iaidx(4, 4)] * e[_iaidx(5, 5)]
                                  - e[_iaidx(4, 5)]**2)), e[_iaidx(5, 5)] * ch)
    else:
        sigz = (abs(alphap) * sige * CVELOC * ph.p0 / ph.h0 / omegaz
                if omegaz != 0.0 else 0.0)
        eemz = sigz * sige
    return {
        'emittance_x': eemx, 'emittance_y': eemy, 'emittance_z': eemz,
        'momentum_spread': sige, 'bunch_length': sigz,
        'momentum_compaction': alphap,
        'rf_voltage_eff': vceff, 'wrfeff': wrfeff, 'trf0': st.trf0, 'bucket_height': bh,
        'closed_orbit': st.cod.copy(),
        'energy_loss_eV': st.u0 * ph.pgev,
        'tunes': params['tunes'], 'damping_rate': params['damping'],
        'beam_matrix': emitp, 'one_turn_matrix': T + st.D,
        'rx': rx,                                 # temit の rx (RFSW OFF なら z 行を 0,1 に置き換えた T)。temitprint の symplectic 偏差用
        'normal_r': r, 'normal_ri': ri,          # temw の r / ri (正規モード座標 → 物理座標とその逆)
        # temitprint の報告用 (Extended Twiss, h, fs, Imag.tune, Tune shift, Damping time/partition)
        'cd': cd.copy(), 'omegaz': omegaz, 'heff': heff, 'omega0': ph.omega0, 'h0': ph.h0,
        'ext_twiss': _ext_twiss(ri, st.cod),
    }


def _ext_twiss(ri, cod):
    """tecalc: params(iptwiss:..) = tfetwiss(ri, cod, .true.) (temitprint の Extended Twiss Parameters)。cod は最後の一周の出口"""
    from .calc6d import tfetwiss
    try:
        return np.asarray(tfetwiss(ri, cod, True), dtype=float)
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return np.full(28, float('nan'))


# ---------------------------------------------------------------------------
# tcod (tcod.f): RADCOD ON の 6D 閉軌道。反復の履歴 (z を trf0 に移す) まで直訳する
# ---------------------------------------------------------------------------

_CODW0 = np.array([1e-6, 1e-7, 1e-7, 1e-8, 1e-5, 1e-6])      # codm: codw0
_CODM = dict(conv0=1e-4, fmin=1e-4, dpthre=3e-4)


def _tsvd_solve(a, b, eps=1e-8):
    """tsvd (tsolvg.f): a は (n, m) で n<=m でもよい。特異値 < max*eps を捨てた最小ノルム解"""
    U, sv, Vt = np.linalg.svd(a, full_matrices=False)
    thr = sv.max() * eps if sv.size else 0.0
    sinv = np.where(sv > thr, 1.0 / np.where(sv > thr, sv, 1.0), 0.0)
    return Vt.T @ (sinv * (U.T @ b))


# 前回の Newton の動きがこの値未満ならヤコビアンを作り直さず前のものを使う (chord 法) こともできるが、
# 軌道が 1.7e-14 変わるだけで εx 等が 2e-8 変わる (量子励起の分割数 int() の境界) ので既定では使わない (0 = 無効)
_REUSE_STEP = 0.0
_R_PRED = 1e-1         # 前回の残差 r がこれ未満なら次の反復は収束しそう → まず cod だけ追い、収束していなければヤコビアン付きでやり直す
_R_FULL = 100.0        # 前回の残差 r がこれ未満なら次の反復が収束しそう → D・B 行列も含めた完全な一周にして最終の追跡と兼ねる


def _tcod_radcod(els, ph: _Phys, cod_in, rfsw: bool, trace: bool = False, final_full: bool = False, radcod: bool = True):
    """tcod.f の直訳 (RADTAPER OFF)。radcod=True は RADCOD ON、False は RADCOD OFF (各反復の tturne が update=.true. で trf0 を更新し、
    収束判定に (dtrf0/trw)^2 が入る。z を trf0 に移さない。docs/history/work_emit_unstable.md C 節)。
    戻り値 (cod, fndcod, trf0, st, st_full)。st は最後の tturne の状態。
    rfsw のときは先に RFSW OFF で再帰呼び出しし、cod(5)=asin((u0 pgev - vcacc)/vceff)/wrfeff - trf0 とおいてから本体に入る。
    本体は SAD と同じく Newton 法 (ヤコビアン T-1 を特異値分解で解く) で、一周ごとに z の平均 dz を
    codi(5) から引いて trf0 に足す (RADCOD のとき RF の位相の基準を動かす)。収束判定は r = Σ(dcod1/codw)² < 1e-4。
    SAD の収束判定は緩いので、z と trf0 の最終値は反復の履歴に依存する (Equilibrium dz の 2e-5 の差はここから)。

    速度のための工夫 (反復の経路と数値は SAD の直訳と同じ。ビット単位で一致):
      - 前回の r < _R_PRED (収束しそう) の反復は、まず cod だけ追う一周 (want_jac=False) で判定し、収束していなければ
        同じ cod でヤコビアン付きの一周をやり直す。ヤコビアンは収束判定の後の Newton 法でしか使わない
      - final_full=True のとき、前回の r < _R_FULL の反復は D・B も作る完全な一周にし、収束したらその状態を st_full として返す
        (呼び出し側は temit の最後の tturne をやり直さなくてよい)"""
    trf0 = 0.0
    cod = np.array(cod_in, dtype=float)
    if rfsw:
        cod, fnd, trf0, st, _ = _tcod_radcod(els, ph, cod, False, trace, radcod=radcod)
        if fnd:
            wrf, vce = _rf_effective(st)
            if vce != 0.0:
                cod[4] = _asinz((st.u0 * ph.pgev - st.vcacc) / vce) / wrf - trf0
        im = 6
    else:
        im = 5
    v0 = ph.p0 / ph.h0
    conv = _CODM['conv0']
    fmin = _CODM['fmin']
    dpthre = _CODM['dpthre']
    codw = _CODW0.copy()
    if not rfsw:
        codw = codw * 0.01
    dcod = np.zeros(6)
    dcod0 = np.zeros(6)
    r0 = 1e200
    fact = 1.0
    loop = 100
    codi = cod.copy()
    st = None
    trs_cache = None
    prev_step = np.inf
    r_prev = None
    trf0s = trf0
    trw = codw[4]
    while True:
        loop -= 1
        if loop <= 0:
            return cod, False, trf0, st, None
        full = radcod and final_full and rfsw and r_prev is not None and r_prev < _R_FULL
        need_j = trs_cache is None or prev_step >= _REUSE_STEP
        cod_first = radcod and (not full) and r_prev is not None and r_prev < _R_PRED
        trf00 = trf0
        if full:
            st = tturne(els, codi, ph, trf0=trf0, rfsw=rfsw, radcod=True)
        elif cod_first or (radcod and not need_j):
            st = tturne(els, codi, ph, trf0=trf0, rfsw=rfsw, radcod=radcod, want_jac=False)
        else:
            st = tturne(els, codi, ph, trf0=trf0, rfsw=rfsw, radcod=radcod, want_db=False)
        if not radcod:
            trf0 = _tturne_update_trf0(st, trf0)        # tturne(update=.true.): RADCOD OFF では毎回 trf0 を更新する
        codf = st.cod
        dz = (codi[4] + codf[4]) * 0.5
        dcod1 = codi - codf
        dtrf0 = trf0 - trf00
        if not rfsw:
            dcod1[4] = 0.0
        r = float(np.sum((dcod1[:im] / codw[:im]) ** 2))
        if not radcod:
            r = r + (dtrf0 / trw) ** 2
        if trace:
            kind = 'full' if full else ('cod' if (cod_first or (radcod and not need_j)) else 'jac')
            print(f'  tcod rfsw={rfsw} loop={100 - loop} [{kind}] r={r:.3e} fact={fact:.3g} z={codi[4]:.4e} trf0={trf0:.9e}')
        if r < conv:
            return codi, True, trf0, st, (st if full else None)
        if r != r:
            trf0 = trf0s
            return cod, False, trf0, st, None
        r_prev = r
        if not st.want_jac and need_j:
            # 収束していなかった: 同じ cod でヤコビアン付きの一周をやり直す
            st = tturne(els, codi, ph, trf0=trf0, rfsw=rfsw, radcod=radcod, want_db=False)
        if st.want_jac:
            trs_cache = st.T - np.eye(6)
        if r >= r0:
            fact *= 0.5
            codi = cod + fact * dcod
            prev_step = float(np.max(np.abs(fact * dcod)))
            if fact < fmin:
                cod = codi
                return cod, False, trf0, st, None
            continue
        red = r / r0
        trf0s = trf0
        r0 = r
        sdot = float(np.sum(dcod * dcod0 / codw ** 2))
        if red < 1.0:
            fact = min(fact * (2.0 - max(red - 0.7, 0.0) / 0.3), 1.0)
        elif sdot < 0.0:
            fact = max(fact * 0.5, fmin)
        trs = trs_cache.copy()
        if not rfsw:
            trs[:, 4] = 0.0
            trs[:, 5] = 0.0
            dcod1[5] = 0.0
            dcod1[4] = 0.0
        dcod0 = dcod
        ddp = dcod1[5]
        dcod = _tsvd_solve(trs[:im, :], dcod1[:im])
        if radcod:
            wrf, _ = _rf_effective(st)
            alambdarf = 2.0 * math.pi / wrf
            codi = codi.copy()
            codi[4] = codi[4] - dz                       # dtrf0 = 0 (RADCOD では tturne が trf0 を更新しない)
            trf0 = trf0 + dz / v0
            sg = 1.0 if trf0 >= 0.0 else -1.0            # Fortran sign(1.d0, x)
            trf0 = sg * math.fmod(sg * (trf0 - 0.5 * alambdarf), alambdarf) + 0.5 * alambdarf
            if abs(dz * st.dvcacc / ph.pgev) + abs(ddp) > dpthre:
                dcod = dcod.copy()
                dcod[:4] = 0.0
        cod = codi.copy()
        prev_step = float(np.max(np.abs(fact * dcod)))
        codi = codi + fact * dcod


def _qcod_accepts_start(els, dp, dcod4):
    """lattice._qcod_periodic が初期値 (横方向 0) を最初の反復で受理するか。受理条件 r <= conv (qtwiss.f qcod) の r を、
    線形の一周行列 (解析式, 1.6e-2 s) の Courant-Snyder パラメータで評価する。誤差のない理想ラティスでは受理される
    (qcod は受理すると初期値をそのまま返す)。境界 (r が conv の 10 分の 1 以上) は受理しない = 本物の qcod に任せる。"""
    if np.max(np.abs(dcod4)) < 1e-13:
        return True
    if np.max(np.abs(dcod4)) > 1e-6:
        return False
    M = np.eye(6)
    for el in els:
        M = el.transfer_matrix(dp) @ M
    tr = M
    cx = 0.5 * (tr[0, 0] + tr[1, 1])
    cy = 0.5 * (tr[2, 2] + tr[3, 3])
    if abs(cx) > 1.0:
        cx = 1.0 / cx
    if abs(cy) > 1.0:
        cy = 1.0 / cy
    sx = math.copysign(math.sqrt(max(1e-6, 1.0 - cx * cx)), tr[0, 1])
    ax = 0.5 * (tr[0, 0] - tr[1, 1]) / sx
    bx = tr[0, 1] / sx
    sy = math.copysign(math.sqrt(max(1e-6, 1.0 - cy * cy)), tr[2, 3])
    ay = 0.5 * (tr[2, 2] - tr[3, 3]) / sy
    by = tr[2, 3] / sy
    d = dcod4
    r = (d[0] ** 2 / bx + bx * (d[1] + ax / bx * d[0]) ** 2
         + d[2] ** 2 / by + by * (d[3] + ay / by * d[2]) ** 2)
    conv = min(1e-10, 1e-18 * (1.0 + (dp / 0.001) ** 2 + len(els) + 1))
    return bool(np.isfinite(r) and r <= 0.1 * conv)


def beam_tilt_sizes(emitp):
    """temitprint (temit.f:1126-1141) の Beam tilt / Beam size xi, eta。emitp は入口の平衡ビーム行列 (6x6, tecalc の 'beam_matrix')。
    btilt = atan2(-2 σxy, σxx-σyy)/2 (どちらかが 0 なら 0)、xi/eta は主軸方向の rms ビームサイズ [m]。
    戻り値 (btilt, sigx, sigy)。"""
    xxs = emitp[0, 0] - emitp[2, 2]
    yys = -2.0 * emitp[0, 2]
    btilt = 0.5 * math.atan2(yys, xxs) if (xxs != 0.0 and yys != 0.0) else 0.0
    sig1 = abs(emitp[0, 0] + emitp[2, 2]) / 2.0
    sig2 = 0.5 * math.hypot(emitp[0, 0] - emitp[2, 2], 2.0 * emitp[0, 2])
    sigx = max(math.sqrt(sig1 + sig2), math.sqrt(abs(sig1 - sig2)))
    sigy = min(math.sqrt(sig1 + sig2), math.sqrt(abs(sig1 - sig2)))
    return btilt, sigx, sigy


def closed_orbit6d(els, ph: _Phys, rfsw: bool, radcod: bool, close_z: bool = True, norf_trf0: bool = False):
    """temit の tcod 部分 (放射あり)。戻り値 (cod, trf0, found, st_full)。
    st_full は RFSW+RADCOD の tcod の最後の反復が D・B 付きの完全な一周のときその状態 (temit の最後の tturne と同じ)、なければ None。

      RFSW ON, RADCOD ON  : tcod の直訳 (_tcod_radcod)。CALC (4D の閉軌道, z=0, dp=0) を初期値に、z のずれを trf0 に移しながら反復する。
                            収束判定が緩い (r<1e-4) ので cod の z と trf0 (Equilibrium dz) の最終値は反復の履歴で決まる。
                            空洞がない・収束しないときは trf0 を固定して放射込みの一周写像の 6D 固定点を Newton 法で解く (z と trf0 の分け方は SAD と違う)
      RFSW ON, RADCOD OFF : 4D 閉軌道 (δ 固定) で一周して z が閉じる δ を割線法で求める (z=0)。誤差で周長が変わると δ≠0
                            (SAD Level 7: DELX/DTHETA 入りで δ=-3.21e-6, z=0 を確認)。trf0 は 1 周目 (trf0=0) の U0・実効電圧から (tturne の update)
      RFSW OFF            : 同じ z 閉じ。trf0 は norf_trf0=True のときだけここで求める (SAD は RFSW OFF でも tturne の update で更新する)
    """
    if rfsw and radcod:
        c, _, _ = _lat._qcod_periodic(els, 0.0)
        cod0 = np.array(c, dtype=float)
        cod0[4] = 0.0
        cod0[5] = 0.0
        if any(type(el).__name__ == 'Cavity' for el in els):
            cod_t, fnd, trf0_t, _, st_full = _tcod_radcod(els, ph, cod0, True, final_full=True)
            if fnd:
                return cod_t, trf0_t, True, st_full
        st1 = tturne(els, cod0, ph, trf0=0.0, rfsw=True)
        wrf, vce = _rf_effective(st1)
        trf0 = _trf0_update(st1, wrf, vce) if wrf != 0.0 else 0.0
        cod = cod0.copy()
        res = np.ones(6)
        found = False
        for _ in range(30):
            st = tturne(els, cod, ph, trf0=trf0, rfsw=True, radcod=True)
            res = st.cod - cod
            if np.max(np.abs(res)) < 1e-15:
                found = True
                break
            M = st.T + st.D
            cod = cod - np.linalg.solve(M - np.eye(6), res)
        if not found:
            found = bool(np.max(np.abs(res)) < 1e-9)
        return cod, trf0, found, None

    if rfsw and not radcod and any(type(el).__name__ == 'Cavity' for el in els):
        # RFSW ON, RADCOD OFF: tcod の直訳 (_tcod_radcod(radcod=False))。RFSW OFF の 4D 閉軌道 (trf0 は各反復の tturne の update で更新) →
        # cod(5)=asin((u0 pgev - vcacc)/vceff)/wrfeff - trf0 → 6D の Newton 法 (FSHIFT≠0 では δ が動き z を閉じる)。
        # 収束判定は r=Σ(dcod1/codw)²+(dtrf0/trw)² < 1e-4 で緩く、最終値は反復の履歴に依存する (以前は z を閉じる δ の割線法で 5e-5 ずれた)
        c, _, _ = _lat._qcod_periodic(els, 0.0)
        cod0 = np.array(c, dtype=float)
        cod0[4] = 0.0
        cod0[5] = 0.0
        cod_t, fnd, trf0_t, _, _ = _tcod_radcod(els, ph, cod0, True, radcod=False)
        if fnd:
            return cod_t, trf0_t, True, None
    if radcod and not rfsw:
        # RFSW OFF, RADCOD ON: SAD の tcod は rfsw=.false. の枝 (im=5) をそのまま実行する。一周の放射損失 (δ が U0/P0 だけ下がる) を
        # 含めた閉軌道 (x,px,y,py を閉じ、z のずれは trf0 に移す)。_tcod_radcod の第 1 段 (RFSW ON のときの前半) と同じ
        c, _, _ = _lat._qcod_periodic(els, 0.0)
        cod0 = np.array(c, dtype=float)
        cod0[4] = 0.0
        cod0[5] = 0.0
        cod_t, fnd, trf0_t, _, _ = _tcod_radcod(els, ph, cod0, False)
        if fnd:
            return cod_t, trf0_t, True, None
    # SAD の Emittance[] は z も閉じる 6D 閉軌道を使う (RF なしでは δ が決まる)。tcod の直訳ではなく、
    # 「4D 閉軌道 (δ 固定) を求めて一周し、z の戻りが 0 になる δ」を割線法で求める近似。
    # z 写像は放射なしの COD 写像 (RADCOD=F なので cod は放射で変わらない)。
    def z_slip(dp):
        # まず軌道 0 から cod だけ追う一周を試す。横方向が閉じていれば (誤差のない理想ラティス)、
        # qcod (_qcod_periodic) も最初の反復で受理して初期値の 0 を返すので、ヤコビアン付きの反復を省く
        c = np.zeros(6)
        c[5] = dp
        st0 = tturne(els, c, ph, want_jac=False)      # z の戻りだけ見るので cod だけ追う
        if _qcod_accepts_start(els, dp, st0.cod[:4]):
            return st0.cod[4] - c[4], c
        c, _, _ = _lat._qcod_periodic(els, dp)
        c = np.array(c, dtype=float)
        c[4] = 0.0
        c[5] = dp
        st0 = tturne(els, c, ph, want_jac=False)
        return st0.cod[4] - c[4], c
    dp0 = 0.0
    z0, cod0 = z_slip(dp0)
    if close_z and z0 != 0.0:
        dp1 = 1e-7
        z1, cod1 = z_slip(dp1)
        for _ in range(30):      # 4 回では K0 の総和が 0 でないとき (δ が 1e-3 近く動く) 収束しない
            if z1 == z0:
                break
            dp2 = dp1 - z1 * (dp1 - dp0) / (z1 - z0)
            dp0, z0 = dp1, z1
            dp1 = dp2
            z1, cod1 = z_slip(dp1)
            if abs(z1) < 1e-15:
                break
        cod0 = cod1
    cod0 = np.array(cod0, dtype=float)
    cod0[4] = 0.0
    if rfsw and not close_z:
        c, _, _ = _lat._qcod_periodic(els, 0.0)
        cod0 = np.array(c, dtype=float)
        cod0[4] = 0.0
        cod0[5] = 0.0
    trf0 = 0.0
    if rfsw or norf_trf0:
        st1 = tturne(els, cod0, ph, trf0=0.0, rfsw=rfsw, want_db=False)      # U0・RF の積算量と T(5,6) の符号だけ必要
        wrf, vce = _rf_effective(st1)
        trf0 = _trf0_update(st1, wrf, vce) if wrf != 0.0 else 0.0
    return cod0, trf0, True, None


def emittance6d(lattice, momentum_eV: float, charge: float = 1.0,
                amass: float = ELMASS, close_z: bool = True,
                rfsw: bool = False, radcod: bool = False, fshift: float = 0.0, emit_floor=(0.0, 0.0, 0.0)) -> dict:
    """
    SAD の Emittance[] 相当 (temit: tcod → tturne (irad=12) → tecalc)。emit_floor = (EMITXE, EMITYE, EMITZE)。

    Parameters
    ----------
    lattice     : pysad.lattice.Lattice
    momentum_eV : 運動量 [eV/c] (SAD の MOMENTUM)
    rfsw, radcod: FFS の RFSW / RADCOD フラグ (閉軌道の求め方は closed_orbit6d)

    Returns dict: emittance_x, emittance_y, momentum_spread,
                  momentum_compaction, energy_loss_eV, ...
    """
    els = lattice.elements
    circ = sum(el.L for el in els)
    ph = _Phys(momentum_eV, circ, amass=amass, charge=charge)
    ph.dvfs = fshift if rfsw else 0.0          # tsetdvfs
    cod, trf0, _found, st_full = closed_orbit6d(els, ph, rfsw, radcod, close_z)
    # temit の最後の tturne (irad=12, update なし): tcod の最後の反復が D・B 付きの完全な一周ならそれと同じ。
    # RFSW OFF でも RADCOD ON なら放射込みで追う (SAD の tturne は rfsw に関係なく radcod を見る。以前は RADCOD OFF として追っていた)
    st = st_full if st_full is not None else tturne(els, cod, ph, trf0=trf0, rfsw=rfsw, radcod=radcod)
    out = tecalc(st, circ, rfsw=rfsw, emit_floor=emit_floor)
    # temitprint の Closed orbit の Entrance (tcod の cod) と Exit (最後の一周の出口)
    out['orbit_entrance'] = np.array(cod, dtype=float)
    out['orbit_exit'] = np.array(st.cod, dtype=float)
    if radcod:
        out['closed_orbit'] = cod.copy()
    elif not rfsw and not radcod:
        # RFSW OFF でも SAD の tcod の tturne (update=.true.) は trf0 を更新するので、Equilibrium dz が報告される
        wrf, vce = _rf_effective(st)
        if wrf != 0.0:
            out['trf0'] = _trf0_update(st, wrf, vce)
    return out
