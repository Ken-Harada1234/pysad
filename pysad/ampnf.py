# -*- coding: utf-8 -*-
"""v2.4.0 A1: 振幅・δ に依る (J²、δ 依存、交差項まで) チューンを、ハミルトニアン模型の 1 周の写像のべき級数 (TPSA) の Birkhoff 正規形で求める。

【何を求めるか】 ν_x(Jx, Jy, δ)、ν_y(Jx, Jy, δ) を J・δ のべき級数 Σ c_abp Jx^a Jy^b δ^p として、全次数 a・b・p の係数。
  J¹ の係数 (a+b = 1、p = 0) = ADTS (resdrive.adts と同じ量)、J の係数 (a=b=0、p) = クロマティシティ ξ_p (∂ᵖν/∂δᵖ = p! c_00p)、
  J^a δ^p (a ≥ 1、p ≥ 1) = ADTS の δ 依存 (振幅と δ の交差項)、J² (a+b = 2、p = 0) = ADTS の振幅 4 次の項 (チューンの振幅の 4 次)。
【摂動論の次数】 模型の 1 周の写像 (x, px, y, py, δ の 5 変数) を全次数 D (既定 5) までべき級数で厳密に作る (近似は次数の切り捨てだけ)。次数 D の写像は、ハミルトニアンの
  (x, px, y, py, δ の合計) D+1 次までの項を厳密に含む。ν の係数 c_abp のうち 2a+2b+p ≤ D−1 が正しい。強さのべきでは、J¹ が 6 極の 2 次 (K2²)・8 極の 1 次 (K3)、
  J² が 6 極の 4 次 (K2⁴)・K2² K3・K3²、ξ_p が 6 極の p 次まで。正規形は写像の近恒等変換を次数ごとに非共鳴の項を消して求める (Birkhoff、共鳴の分母 λ^α − λ_i)。
  δ は動かない 5 番目の変数 (固有値 1) で、n = 1 の δ の単項式 (分散による閉軌道のずれ) もこの変換で消える。追跡の値に合わせる項・経験的な重みは入れていない。
【模型 (入れている効果)】 線形の素子 (DRIFT・QUAD・BEND のセクター型 E1・E2 つき) は H = (px²+py²)/(2(1+δ)) + ((K1+h²)x² − K1 y²)/2 − hδx の厳密な転送行列
  (δ のべきに展開した係数。複素 δ の円周上の値からの Cauchy の公式)。薄い・厚い (nslice 分割) 6 極・8 極の kick (強さは δ に依らない)。オプション:
  bend_cubic=True: 偏向磁石の曲率の 3 次の項 (1+δ)(h x/2) u (u = (px²+py²)/(1+δ)²) と K1 の 3 次の項 ('log' が既定、compat.set_bend_k1_def('lin') で 'lin')、
    E1・E2≠0 の端の 3 次の生成関数 (bendnl.bend_edge_generators)。
  kin=True: 運動学的な項 H_kin = (1+δ) Σ_{k≥2} c_k u^k (−√(1−u) の展開。ドリフト・四極・偏向磁石・6 極のドリフト部。偏向磁石は (1+hx) 倍)。
  fringe=True: 四極と K1≠0 の偏向磁石の端の Lee–Whiting のフリンジ H_f = σ (K1/12)[(x³+3xy²)px − (y³+3x²y)py]/(1+δ)^FRINGE_POWER (FRINGE_POWER = 1)。
  非線形の項はハミルトン流を Lie 級数 exp(ds {·, H}) で次数 D まで厳密に求め、Strang 分割で線形の部分と組み合わせる (nslice・nslice_bend・nslice_quad、誤差 O(1/n²))。
【v2.4.0 F1】 ROTATE・ROTATION した BEND・QUAD・SEXT・OCT は R⁻¹ ∘ (回していない素子) ∘ R で入る (全次数の θ)。結合した線形光学は 4×4 の固有モード (linear_modes) の正規形で扱う。
  rdt_series は 1 周の写像の 2 次の部分から共鳴駆動項 h を出す (結合した光学でも動く)。docs/49 §14。
【v2.4.0 G】 K0≠0 の偏向磁石 (K1 = 0。挿入光源・逆曲げ・MPW・傾けた MPW) は、K64 の追跡の閉じた式を磁石内の軌道のまわりで TPSA 評価した写像 (ampnf_k0.py、op 'tmap') で入る。
  v2.4.1: その非線形を他の素子のオプションと独立に入れる・落とす `k0_nonlinear` (True・False。None = 従来: 3 つのオプションが全部 False のときだけ落とす)。
  閉軌道が原点にない 1 周の写像は 2 段 (不動点 c を求め、c を定数項に入れて通す。res['orbit'])。docs/49 §14.4。
【入れていない効果】 K1≠0 かつ K0≠0 の偏向磁石 (NotImplementedError)、F1・DX・DY・DROTATE、誤差、K1 の偏向磁石の 4 次以上、K64 の四極の積分器の台形則、10 極以上、放射・RF。
  差の原因は docs/49 §13.3。使い方・検証は docs/49 §13、記録は docs/history/v24_A.md。
"""
import math
import itertools
import numpy as np

_TAB = {}


class _Tab:
    """nvar 変数・次数 D 以下の単項式の表と、積の添字表"""

    def __init__(self, nvar, D):
        self.nvar, self.D = nvar, D
        mons = []
        for d in range(D + 1):
            for e in _comps(nvar, d):
                mons.append(e)
        self.mons = np.array(mons, dtype=np.int64)
        self.N = len(mons)
        self.deg = self.mons.sum(axis=1)
        self.idx = {tuple(e): i for i, e in enumerate(mons)}
        I, J, K = [], [], []
        for i, a in enumerate(mons):
            da = self.deg[i]
            for j, b in enumerate(mons):
                if da + self.deg[j] <= D:
                    I.append(i)
                    J.append(j)
                    K.append(self.idx[tuple(x + y for x, y in zip(a, b))])
        self.I, self.J, self.K = np.array(I), np.array(J), np.array(K)

    def mul(self, a, b):
        nz_a = np.nonzero(a)[0]
        if len(nz_a) == 0:
            return np.zeros(self.N, dtype=complex)
        out = np.zeros(self.N, dtype=complex)
        # a の非零の項だけ: 各 i について (j, k) の表
        sel = np.isin(self.I, nz_a)
        I, J, K = self.I[sel], self.J[sel], self.K[sel]
        w = a[I] * b[J]
        out.real = np.bincount(K, weights=w.real, minlength=self.N)
        out.imag = np.bincount(K, weights=w.imag, minlength=self.N)
        return out

    def var(self, j, c=1.0):
        v = np.zeros(self.N, dtype=complex)
        e = [0] * self.nvar
        e[j] = 1
        v[self.idx[tuple(e)]] = c
        return v

    def const(self, c):
        v = np.zeros(self.N, dtype=complex)
        v[0] = c
        return v

    def shift_delta(self, p):
        """δ^p をかける添字の表 (δ = 最後の変数)"""
        key = ('sh', p)
        if key in self.__dict__:
            return self.__dict__[key]
        src, dst = [], []
        for i, e in enumerate(self.mons):
            if self.deg[i] + p <= self.D:
                e2 = list(e)
                e2[-1] += p
                src.append(i)
                dst.append(self.idx[tuple(e2)])
        self.__dict__[key] = (np.array(src), np.array(dst))
        return self.__dict__[key]

    def deriv(self, f, j):
        out = np.zeros(self.N, dtype=complex)
        for i in np.nonzero(f)[0]:
            e = self.mons[i]
            if e[j] > 0:
                e2 = list(e)
                e2[j] -= 1
                out[self.idx[tuple(e2)]] += e[j] * f[i]
        return out

    def poisson(self, f, g):
        """{f, g} = sum (df/dq dg/dp - df/dp dg/dq)、(q, p) = (0,1), (2,3)"""
        out = np.zeros(self.N, dtype=complex)
        for q, pp in ((0, 1), (2, 3)):
            out += self.mul(self.deriv(f, q), self.deriv(g, pp)) - self.mul(self.deriv(f, pp), self.deriv(g, q))
        return out

    def flow(self, H, ds):
        """H (物理の 5 変数の TPSA、δ は変数 4) のハミルトン流 (dq/ds = dH/dp) の ds の写像 Phi_j = exp(ds L) q_j、L f = {f, H} (j = x, px, y, py)。
        H の次数が 3 以上なので k 回の括弧は次数が上がり、k <= D-1 で厳密 (次数 D まで)"""
        out = []
        for j in range(4):
            term = self.var(j)
            acc = term.copy()
            for k in range(1, self.D + 1):
                term = self.poisson(term, H) * (ds / k)
                if not np.any(term):
                    break
                acc = acc + term
            out.append(acc)
        return out

    def mul_dser(self, c, a):
        """δ の級数 c[p] (p=0..D) × a"""
        out = np.zeros(self.N, dtype=complex)
        for p, cp in enumerate(c):
            if cp == 0:
                continue
            s, d = self.shift_delta(p)
            out[d] += cp * a[s]
        return out


def _comps(nvar, d):
    if nvar == 1:
        yield (d,)
        return
    for k in range(d, -1, -1):
        for rest in _comps(nvar - 1, d - k):
            yield (k,) + rest


def get_tab(nvar, D):
    if (nvar, D) not in _TAB:
        _TAB[(nvar, D)] = _Tab(nvar, D)
    return _TAB[(nvar, D)]


# ---------------------------------------------------------------------------------------------
# 線形の素子の δ の級数 (Cauchy の公式)
# ---------------------------------------------------------------------------------------------
NTH = 64
RAD = 0.25


def _ent(k, L):
    """c = cos(√k L)、S = sin(√k L)/√k、C1 = (1−cos(√k L))/k (k の整関数)。k は複素の配列"""
    z = k * L * L
    small = np.abs(z) < 0.5
    c = np.empty_like(z, dtype=complex)
    S = np.empty_like(z, dtype=complex)
    C1 = np.empty_like(z, dtype=complex)
    zs = z[small]
    # 級数: cos = Σ (−z)^n/(2n)!、S/L = Σ (−z)^n/(2n+1)!、C1/L² = Σ (−1)^n z^n/(2n+2)!
    cs = np.zeros_like(zs)
    Ss = np.zeros_like(zs)
    Cs = np.zeros_like(zs)
    t = np.ones_like(zs)
    for n in range(14):
        cs += t / math.factorial(2 * n)
        Ss += t / math.factorial(2 * n + 1)
        Cs += t / math.factorial(2 * n + 2)
        t = -t * zs
    c[small], S[small], C1[small] = cs, Ss * L, Cs * L * L
    big = ~small
    if big.any():
        kb = k[big]
        w = np.sqrt(kb + 0j)
        c[big] = np.cos(w * L)
        S[big] = np.sin(w * L) / w
        C1[big] = (1 - np.cos(w * L)) / kb
    return c, S, C1


def _plane(Kc, L, h, dl):
    """1 平面 (x, px) の 2×2 行列と強制項。H = p²/(2(1+δ)) + Kc x²/2 − h δ x。dl: 複素 δ の配列。戻り M (n,2,2)、d (n,2)"""
    k = Kc / (1 + dl)
    c, S, C1 = _ent(np.asarray(k, dtype=complex), L)
    n = len(dl)
    M = np.zeros((n, 2, 2), dtype=complex)
    M[:, 0, 0] = c
    M[:, 0, 1] = S / (1 + dl)
    M[:, 1, 0] = -Kc * S
    M[:, 1, 1] = c
    d = np.zeros((n, 2), dtype=complex)
    if h:
        d[:, 0] = h * dl * C1 / (1 + dl)
        d[:, 1] = h * dl * S
    return M, d


def _lin4(Mx, dx, My, dy):
    n = Mx.shape[0]
    M = np.zeros((n, 4, 4), dtype=complex)
    M[:, :2, :2] = Mx
    M[:, 2:, 2:] = My
    d = np.concatenate([dx, dy], axis=1)
    return M, d


def rot_matrix(theta):
    """ROTATE (ロール角 θ) の座標変換 R (4×4): 素子座標 = R · 実験室座標 (x' = cosθ x − sinθ y、y' = sinθ x + cosθ y、px・py も同じ。K64 の tchge と同じ向き。
    ROTATE した素子の写像は R⁻¹ ∘ (回していない素子の写像) ∘ R (pysad の追跡のヤコビ行列と 2e-15 で一致、v2.4.0 F1)"""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, 0.0, -s, 0.0], [0.0, c, 0.0, -s], [s, 0.0, c, 0.0], [0.0, s, 0.0, c]])


def roll_angle(el):
    """素子のロール角 ROTATE + ROTATION [rad] (BEND・QUAD・SEXT・OCT)。DROTATE≠0 は未対応"""
    if float(getattr(el, 'DROTATE', 0.0) or 0.0):
        raise NotImplementedError('ampnf: DROTATE≠0: %s' % el.name)
    return float(getattr(el, 'ROTATE', 0.0) or 0.0) + float(getattr(el, 'ROTATION', 0.0) or 0.0)


def _expand_rolls(elements):
    """ロール角のある BEND・QUAD・SEXT・OCT を ('rot', θ)、ロール角 0 の素子、('rot', −θ) に展開する (v2.4.0 F1: 結合した光学・垂直分散を通す)。
    v2.4.0 G: K0≠0 の BEND は、ANGLE≠0 ならそのまま (ampnf_k0 の写像が DX・DY・ROTATE・DROTATE・CHI2 を tbshift・tbrot の式で内部で扱う)、
    ANGLE=0 (ステアリング・MPW 型。K64 は tent・texit) は ('shift', −DX, −DY)・('rot', θ + DROTATE)・素子・('rot', −θ)・('shift', DX, DY) に展開する"""
    import dataclasses
    import copy
    out = []
    for el in elements:
        tn = type(el).__name__
        if tn == 'Bend' and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0:
            if float(getattr(el, 'ANGLE', 0.0) or 0.0) != 0.0 and float(getattr(el, 'L', 0.0) or 0.0) > 0.0:
                out.append(el)
                continue
            g = lambda k: float(getattr(el, k, 0.0) or 0.0)
            th = g('ROTATE') + g('ROTATION') + g('DROTATE')
            dx, dy = g('DX'), g('DY')
            if th == 0.0 and dx == 0.0 and dy == 0.0:
                out.append(el)
                continue
            clean = copy.copy(el)
            for k in ('ROTATE', 'ROTATION', 'DROTATE', 'DX', 'DY'):
                if hasattr(clean, k):
                    setattr(clean, k, 0.0)
            if dx or dy:
                out.append(('shift', -dx, -dy))
            if abs(th) > 1e-15:
                out.append(('rot', th))
            out.append(clean)
            if abs(th) > 1e-15:
                out.append(('rot', -th))
            if dx or dy:
                out.append(('shift', dx, dy))
            continue
        if tn in ('Bend', 'Quadrupole', 'Sextupole', 'Octupole'):
            th = roll_angle(el)
            if abs(th) > 1e-15:
                out += [('rot', th), dataclasses.replace(el, ROTATE=0.0, ROTATION=0.0), ('rot', -th)]
                continue
        out.append(el)
    return out


def element_linear(el, dl):
    """素子の線形の写像 (M (n,4,4), d (n,4)) を複素 δ の配列 dl で返す。None = 非線形の kick だけの素子 (薄い 6 極・8 極)。
    戻り (M, d, kind)。kind: 'drift'|'quad'|'bend'|'other'"""
    tn = type(el).__name__
    L = float(getattr(el, 'L', 0.0) or 0.0)
    n = len(dl)
    if tn == 'Drift':
        Mx, dx = _plane(0.0, L, 0.0, dl)
        My, dy = _plane(0.0, L, 0.0, dl)
        return (*_lin4(Mx, dx, My, dy), 'drift')
    if tn == 'Quadrupole':
        if float(getattr(el, 'K1S', 0.0) or 0.0) or float(getattr(el, 'ROTATE', 0.0) or 0.0):
            raise NotImplementedError('ampnf: スキュー・回転つきの QUAD: %s' % el.name)
        K1 = float(el.K1)
        Mx, dx = _plane(K1, L, 0.0, dl)
        My, dy = _plane(-K1, L, 0.0, dl)
        return (*_lin4(Mx, dx, My, dy), 'quad')
    if tn == 'Bend' and L <= 0:
        A = np.asarray(el.transfer_matrix(0.0), dtype=float)[:4, :4]
        d = np.zeros((n, 4), dtype=complex)
        d[:, 1] = -float(getattr(el, 'K0', 0.0) or 0.0)       # v2.4.0 G: 薄い K0 の kick (K64 の tsteer0: px −= K0/(1+δ) は内部の px/(1+δ)。正準運動量では δ に依らない −K0)
        return np.tile(A.astype(complex), (n, 1, 1)), d, 'other'
    if tn == 'Bend':
        if float(getattr(el, 'K0', 0.0) or 0.0) or float(getattr(el, 'ROTATE', 0.0) or 0.0):
            raise NotImplementedError('ampnf: K0≠0・ROTATE の BEND: %s' % el.name)
        if L <= 0:
            return np.tile(np.eye(4, dtype=complex), (n, 1, 1)), np.zeros((n, 4), dtype=complex), 'other'
        h = float(el.ANGLE) / L
        K1 = float(el.K1)
        Mx, dx = _plane(K1 + h * h, L, h, dl)
        My, dy = _plane(-K1, L, 0.0, dl)
        M, d = _lin4(Mx, dx, My, dy)
        from .elements import _edge_matrix
        if el.E1:
            E = _edge_matrix(h, el.E1, el.HGAP, el.FINT)[:4, :4]
            M = M @ E
        if el.E2:
            E = _edge_matrix(h, el.E2, el.HGAP, el.FINT)[:4, :4]
            M = E @ M
            d = d @ E.T
        return M, d, 'bend'
    # その他: δ に依らない線形の行列 (4×4)
    A = np.asarray(el.transfer_matrix(0.0), dtype=float)[:4, :4]
    return np.tile(A.astype(complex), (n, 1, 1)), np.zeros((n, 4), dtype=complex), 'other'


def _series_from_values(vals, D, dl):
    """vals (n, ...) の複素 δ の円周上の値から Taylor 係数 (D+1, ...)"""
    n = len(dl)
    out = []
    for p in range(D + 1):
        out.append((vals * (dl ** (-p)).reshape((n,) + (1,) * (vals.ndim - 1))).mean(axis=0))
    return np.array(out)


class Seg:
    """線形の区間: x_out = M(δ) x_in + d(δ)。M・d は δ の級数 (D+1, 4, 4)・(D+1, 4)。δ̃ = δ/sd で規格化した係数"""

    def __init__(self, D, sd):
        self.D, self.sd = D, sd
        self.M = np.zeros((D + 1, 4, 4), dtype=complex)
        self.M[0] = np.eye(4)
        self.d = np.zeros((D + 1, 4), dtype=complex)

    @staticmethod
    def from_series(M, d, sd):
        D = M.shape[0] - 1
        s = Seg(D, sd)
        sc = (sd ** np.arange(D + 1))
        s.M = M * sc[:, None, None]
        s.d = d * sc[:, None]
        return s

    def then(self, other):
        """self の後に other"""
        D = self.D
        M = np.zeros_like(self.M)
        d = np.zeros_like(self.d)
        for p in range(D + 1):
            for q in range(D + 1 - p):
                M[p + q] += other.M[q] @ self.M[p]
                d[p + q] += other.M[q] @ self.d[p]
        d += other.d
        out = Seg(D, self.sd)
        out.M, out.d = M, d
        return out


_FLOW = {}


def _bend_h3(tab, h, K1, kin=False):
    """偏向磁石の本体の 3 次のハミルトニアン (物理の 5 変数 x, px, y, py, δ): (h/2) x (px^2+py^2)/(1+δ) + K1 の項 ('log': (K1 h/6) x^3、'lin': (K1 h/3) x^3 - (K1 h/2) x y^2)"""
    from . import compat
    x, px, y, py = (tab.var(j) for j in range(4))
    p2 = tab.mul(px, px) + tab.mul(py, py)
    inv = np.array([(-1.0) ** k for k in range(tab.D + 1)])
    H = 0.5 * h * tab.mul_dser(inv, tab.mul(x, p2))
    x3 = tab.mul(tab.mul(x, x), x)
    if K1:
        if compat.bend_k1_def() == 'lin' and compat.v2():
            H = H + (K1 * h / 3.0) * x3 - 0.5 * K1 * h * tab.mul(x, tab.mul(y, y))
        else:
            H = H + (K1 * h / 6.0) * x3
    if kin:
        hx = h * tab.var(0)
        H = H + _kin_h(tab, hx)
    return H


_KIN_C = [0.5, 1.0 / 8, 1.0 / 16, 5.0 / 128, 7.0 / 256, 21.0 / 1024, 33.0 / 2048]      # −√(1−u) = −1 + Σ c_k u^k


def _kin_h(tab, hx=None):
    """運動学的な項のハミルトニアン (u = (px²+py²)/(1+δ)²): H_kin = (1+δ) Σ_{k≥2} c_k u^k。hx (TPSA、偏向磁石の曲率 h x) があれば (1+hx) 倍 (偏向磁石の座標系の因子)"""
    px, py = tab.var(1), tab.var(3)
    inv2 = np.array([(-1.0) ** k * (k + 1) for k in range(tab.D + 1)])         # (1+δ)^−2
    u = tab.mul_dser(inv2, tab.mul(px, px) + tab.mul(py, py))
    uk = u.copy()
    acc = np.zeros(tab.N, dtype=complex)
    for k in range(2, len(_KIN_C) + 1):
        uk = tab.mul(uk, u)
        if not np.any(uk):
            break
        acc = acc + _KIN_C[k - 1] * uk
    acc = tab.mul_dser(np.array([1.0, 1.0] + [0.0] * (tab.D - 1)), acc)         # (1+δ) 倍
    if hx is not None:
        acc = acc + tab.mul(hx, acc)
    return acc


def _binom_gen(a, k):
    r = 1.0
    for j in range(k):
        r *= (a - j) / (j + 1.0)
    return r


FRINGE_POWER = 1      # 端部の Lee-Whiting のフリンジの強さは K1/(1+δ)^FRINGE_POWER (既定 1: 生成関数はスロープ x' = px/(1+δ) で書かれる。0 は δ に依らない場合の診断。FODO の ξ3 が CALC と 6e-4 で一致するのは 1 だけ)


def _fringe_h(tab, K1, sigma):
    """四極・K1≠0 の偏向磁石の端の Lee-Whiting の非線形フリンジ (K64 の ttfrin4 の 3 次の項、docs/49 §2): H_f = σ (K1/12) [(x³ + 3 x y²) px − (y³ + 3 x² y) py] (σ = +1 入口、−1 出口。δ に依らない)"""
    x, px, y, py = (tab.var(j) for j in range(4))
    x2, y2 = tab.mul(x, x), tab.mul(y, y)
    a = tab.mul(tab.mul(x, x2) + 3.0 * tab.mul(x, y2), px)
    b = tab.mul(tab.mul(y, y2) + 3.0 * tab.mul(x2, y), py)
    H = sigma * K1 / 12.0 * (a - b)
    if FRINGE_POWER:
        c = np.array([_binom_gen(-FRINGE_POWER, k) for k in range(tab.D + 1)])
        H = tab.mul_dser(c, H)
    return H


def _flow_cached(key, tab, Hfunc, ds):
    k = key + (tab.D, float(ds))
    if k not in _FLOW:
        _FLOW[k] = tab.flow(Hfunc(), ds)
    return _FLOW[k]


def _poly_from_gen(tab, g):
    """bendnl の端の生成関数 {(ex, epx, ey, epy): 係数} を TPSA に"""
    out = np.zeros(tab.N, dtype=complex)
    for e, c in g.items():
        out[tab.idx[tuple(e) + (0,)]] += c
    return out


def _edge_seg(h, E, hgap, fint):
    from .elements import _edge_matrix
    if not E:
        return np.eye(4), np.zeros(4)
    return np.asarray(_edge_matrix(h, E, hgap, fint), dtype=float)[:4, :4], np.zeros(4)


def _edge_series(h, E, hgap, fint, D):
    M, d = _edge_seg(h, E, hgap, fint)
    Ms = np.zeros((D + 1, 4, 4), dtype=complex)
    Ms[0] = M
    return Ms, np.zeros((D + 1, 4), dtype=complex)


class NSeg:
    """固定した数値の δ での線形の区間 (M (4,4)、d (4,))"""

    def __init__(self, M=None, d=None):
        self.M = np.eye(4) if M is None else M
        self.d = np.zeros(4) if d is None else d

    def then(self, o):
        return NSeg(o.M @ self.M, o.M @ self.d + o.d)


def build_ops(lat, D, sd, nslice=8, delta=None, bend_cubic=False, nslice_bend=8, Dnl=5, kin=False, fringe=False, nslice_quad=4, k0_nonlinear=None):
    """素子列を [('lin', Seg) | ('kick', k2l, k3l)] に。厚い 6 極・8 極は nslice 分割 (中点則 = 厳密な極限は nslice→∞、誤差 O(1/nslice²))"""
    dl = RAD * np.exp(2j * np.pi * np.arange(NTH) / NTH)

    def seg_of(el):
        if delta is not None:
            M, d, kind = element_linear(el, np.array([delta + 0j]))
            return NSeg(M[0].real, d[0].real), kind
        M, d, kind = element_linear(el, dl)
        return Seg.from_series(_series_from_values(M, D, dl), _series_from_values(d, D, dl), sd), kind

    ops = []
    new0 = (lambda: NSeg()) if delta is not None else (lambda: Seg(D, sd))
    cur = new0()
    warn = []
    for el in _expand_rolls(lat.elements):
        if isinstance(el, tuple):
            if el[0] == 'shift':          # v2.4.0 G: x += dx、y += dy (ANGLE=0 の K0≠0 の磁石の据え付けのずれ)
                Rm = np.eye(4)
                dd = np.array([el[1], 0.0, el[2], 0.0])
            else:
                Rm = rot_matrix(el[1])
                dd = np.zeros(4)
            if delta is not None:
                cur = cur.then(NSeg(Rm, dd))
            else:
                Ms = np.zeros((D + 1, 4, 4), dtype=complex)
                Ms[0] = Rm
                ds_ = np.zeros((D + 1, 4), dtype=complex)
                ds_[0] = dd
                cur = cur.then(Seg.from_series(Ms, ds_, sd))
            continue
        tn = type(el).__name__
        L = float(getattr(el, 'L', 0.0) or 0.0)
        k2l = k3l = 0.0
        if tn == 'Bend' and L > 0 and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0:
            # v2.4.0 G: K0≠0 の偏向磁石は K64 の追跡の閉じた式を軌道のまわりで TPSA 評価した写像 (ampnf_k0)。kin・bend_cubic・fringe のオプションは効かない (全部入り)
            from . import ampnf_k0 as _k0
            ops.append(('lin', cur))
            if float(getattr(el, 'ANGLE', 0.0) or 0.0) == 0.0 and float(getattr(el, 'CHI2', 0.0) or 0.0) != 0.0:
                warn.append('ANGLE=0 の K0≠0 の BEND %s の CHI2 は K64 の追跡 (tsteer・tmulti) が使わないので入れていない (追跡に合わせる)' % el.name)
            # v2.4.1: k0_nonlinear=None は従来どおり (bend_cubic・kin・fringe が全部 False のときだけ 1 次)。True・False は K0≠0 の磁石の自身の非線形を、
            # ほかの素子のオプションと独立に入れる・落とす (曲率・運動学・楔・端は閉じた式に分けられず、個別には切れない。docs/46)
            lin = (not (bend_cubic or kin or fringe)) if k0_nonlinear is None else (not k0_nonlinear)
            ops.append(('tmap', _k0.bend_k0_map(el, linear=lin), el, lin))
            cur = new0()
            continue
        tabn = get_tab(5, max(Dnl, 3)) if (kin or fringe or bend_cubic) else None
        if kin and tn == 'Drift' and L > 0:
            # ドリフトの運動学的な項: 線形のドリフトと可換なので厳密に 1 回の流れ
            sg, _ = seg_of(el)
            cur = cur.then(sg)
            ops.append(('lin', cur))
            ops.append(('nlphi', _flow_cached(('kdrift',), tabn, lambda: _kin_h(tabn), L)))
            cur = new0()
            continue
        if (kin or fringe) and tn == 'Quadrupole' and L > 0 and not float(getattr(el, 'K1S', 0.0) or 0.0) and not float(getattr(el, 'ROTATE', 0.0) or 0.0):
            import dataclasses
            K1q = float(el.K1)
            if fringe and K1q != 0.0:
                ops.append(('lin', cur))
                ops.append(('nlphi', _flow_cached(('qfr', K1q, 1), tabn, lambda: _fringe_h(tabn, K1q, +1.0), 1.0)))
                cur = new0()
            if kin:
                ds = L / nslice_quad
                half, _ = seg_of(dataclasses.replace(el, L=ds / 2))
                full, _ = seg_of(dataclasses.replace(el, L=ds))
                Phi = _flow_cached(('kdrift',), tabn, lambda: _kin_h(tabn), ds)
                cur = cur.then(half)
                for i in range(nslice_quad):
                    ops.append(('lin', cur))
                    ops.append(('nlphi', Phi))
                    cur = full if i < nslice_quad - 1 else half
            else:
                sg, _ = seg_of(el)
                cur = cur.then(sg)
            if fringe and K1q != 0.0:
                ops.append(('lin', cur))
                ops.append(('nlphi', _flow_cached(('qfr', K1q, -1), tabn, lambda: _fringe_h(tabn, K1q, -1.0), 1.0)))
                cur = new0()
            continue
        if tn == 'Bend' and (bend_cubic or kin) and L > 0 and float(el.ANGLE) != 0.0 and not float(getattr(el, 'K0', 0.0) or 0.0) and not float(getattr(el, 'ROTATE', 0.0) or 0.0):
            import dataclasses
            h = float(el.ANGLE) / L
            ds = L / nslice_bend
            tab = get_tab(5, Dnl)
            Phi = _flow_cached(('bend', h, float(el.K1), kin), tab, lambda: _bend_h3(tab, h, float(el.K1), kin), ds)
            from .bendnl import bend_edge_generators
            g1, g2 = bend_edge_generators(el)
            fe1 = tab.flow(_poly_from_gen(tab, g1), 1.0)
            fe2 = tab.flow(_poly_from_gen(tab, g2), 1.0)
            body = dataclasses.replace(el, E1=0.0, E2=0.0)
            half = dataclasses.replace(body, L=ds / 2, ANGLE=float(el.ANGLE) * ds / 2 / L)
            sh, _ = seg_of(half)
            full = dataclasses.replace(body, L=ds, ANGLE=float(el.ANGLE) * ds / L)
            sf, _ = seg_of(full)
            ops.append(('lin', cur))
            if fringe and float(el.K1) != 0.0:
                K1b = float(el.K1)
                ops.append(('nlphi', _flow_cached(('qfr', K1b, 1), tab, lambda: _fringe_h(tab, K1b, +1.0), 1.0)))
            ops.append(('nlphi', fe1))
            if delta is not None:
                e1 = NSeg(*_edge_seg(h, el.E1, el.HGAP, el.FINT))
                e2 = NSeg(*_edge_seg(h, el.E2, el.HGAP, el.FINT))
            else:
                e1 = Seg.from_series(*_edge_series(h, el.E1, el.HGAP, el.FINT, D), sd)
                e2 = Seg.from_series(*_edge_series(h, el.E2, el.HGAP, el.FINT, D), sd)
            cur = e1.then(sh)
            for i in range(nslice_bend):
                ops.append(('lin', cur))
                ops.append(('nlphi', Phi))
                cur = sf if i < nslice_bend - 1 else sh
            ops.append(('lin', cur))
            ops.append(('nlphi', fe2))
            if fringe and float(el.K1) != 0.0:
                K1b = float(el.K1)
                ops.append(('nlphi', _flow_cached(('qfr', K1b, -1), tab, lambda: _fringe_h(tab, K1b, -1.0), 1.0)))
            cur = e2
            continue
        if tn == 'Sextupole':
            k2l = float(el.K2) * (L if L > 0 else 1.0)
        elif tn == 'Octupole':
            k3l = float(el.K3) * (L if L > 0 else 1.0)
        elif tn == 'Multipole':
            knl = getattr(el, 'KnL', {}) or {}
            ksl = getattr(el, 'KsL', {}) or {}
            if any(v for v in ksl.values() if v) or float(getattr(el, 'ROTATION', 0.0) or 0.0):
                raise NotImplementedError('ampnf: スキュー・回転つきの MULT: %s' % el.name)
            k2l, k3l = float(knl.get(2, 0.0)), float(knl.get(3, 0.0))
            if (k2l or k3l) and L > 0 and knl.get(1, 0.0):
                raise NotImplementedError('ampnf: K1≠0 の厚い MULT: %s' % el.name)
            if not (k2l or k3l):
                sg, kind = seg_of(el)
                cur = cur.then(sg)
                continue
        if tn in ('Sextupole', 'Octupole', 'Multipole') and (k2l or k3l):
            from .elements import Drift
            if L > 0:
                h = L / nslice
                sg, _ = seg_of(Drift(name='_d', L=h / 2))
                sg2, _ = seg_of(Drift(name='_d', L=h))
                cur = cur.then(sg)
                for i in range(nslice):
                    ops.append(('lin', cur))
                    if kin:
                        ops.append(('nlphi', _flow_cached(('kdrift',), tabn, lambda: _kin_h(tabn), h)))
                    ops.append(('kick', k2l / nslice, k3l / nslice))
                    cur = sg2 if i < nslice - 1 else sg
            else:
                if tn == 'Multipole' and knl.get(1, 0.0):
                    sg, _ = seg_of(el)
                    cur = cur.then(sg)
                ops.append(('lin', cur))
                ops.append(('kick', k2l, k3l))
                cur = new0()
            continue
        if tn in ('Sextupole', 'Octupole'):       # 強さ 0 の 6 極・8 極はドリフト
            from .elements import Drift
            sg, _ = seg_of(Drift(name='_d', L=L))
            cur = cur.then(sg)
            if kin and L > 0:
                ops.append(('lin', cur))
                ops.append(('nlphi', _flow_cached(('kdrift',), tabn, lambda: _kin_h(tabn), L)))
                cur = new0()
            continue
        sg, kind = seg_of(el)
        if kind == 'other' and L > 0 and tn not in ('Cavity', 'Marker', 'Monitor', 'Mark', 'Moni'):
            warn.append('素子 %s (%s) の線形の写像を δ に依らない行列で入れた' % (el.name, tn))
        cur = cur.then(sg)
    ops.append(('lin', cur))
    return ops, warn


def linear_modes(A):
    """1 周の線形写像 A (4×4、実、シンプレクティック) の正規化: 物理座標 = Am · zn、zn = (x̂, p̂, ŷ, q̂)、モード ζ = x̂ − i p̂ (固有値 e^{iμ})。
    横方向の結合がなければ (A の非対角ブロック < 1e-12) 各面の Twiss (α, β)、結合があれば 4×4 の固有ベクトル (v2.4.0 F1): Am = [Re e1, Im e1, Re e2, Im e2]、
    e は ē^T S e = 2i (S = (0 1; −1 0) の直和。正のクレイン符号) に規格化し、e1 は x 成分、e2 は y 成分が実で正になる位相。
    戻り (Am, Ai, (μ1, μ2) ∈ [0, 2π)、tw = (αx, βx, αy, βy) (結合があれば None)、結合の大きさ)"""
    cpl = float(max(np.abs(A[:2, 2:]).max(), np.abs(A[2:, :2]).max()))
    if cpl < 1e-12:
        tw, mus = [], []
        for b in (0, 2):
            m = A[b:b + 2, b:b + 2]
            cm = 0.5 * (m[0, 0] + m[1, 1])
            sm = math.copysign(math.sqrt(max(1 - cm * cm, 0.0)), m[0, 1])
            mus.append(math.atan2(sm, cm) % (2 * math.pi))
            tw += [(m[0, 0] - m[1, 1]) / (2 * sm), m[0, 1] / sm]
        (ax, bx, ay, by) = tw
        Am = np.array([[math.sqrt(bx), 0, 0, 0], [-ax / math.sqrt(bx), 1 / math.sqrt(bx), 0, 0],
                       [0, 0, math.sqrt(by), 0], [0, 0, -ay / math.sqrt(by), 1 / math.sqrt(by)]])
        return Am, np.linalg.inv(Am), tuple(mus), tuple(tw), cpl
    S4 = np.zeros((4, 4))
    S4[0, 1], S4[1, 0], S4[2, 3], S4[3, 2] = 1.0, -1.0, 1.0, -1.0
    w, V = np.linalg.eig(A)
    if np.abs(np.abs(w) - 1.0).max() > 1e-6:
        raise ValueError('ampnf: 結合した線形の光学が不安定 (固有値の絶対値 %s)' % np.abs(w))
    modes = []
    for k in range(4):
        e = V[:, k]
        kr = (np.conj(e) @ S4 @ e) / 1j          # 実数。正なら正のクレイン符号
        if kr.real > 0:
            modes.append((np.angle(w[k]) % (2 * math.pi), e * math.sqrt(2.0 / kr.real)))
    if len(modes) != 2:
        raise ValueError('ampnf: 結合した光学の正のクレイン符号のモードが 2 つでない (共鳴の上か不安定)')
    (m1, e1), (m2, e2) = modes
    if abs(e1[0]) * abs(e2[2]) < abs(e2[0]) * abs(e1[2]):       # 横 (x) 成分が大きい方を x のモードに
        (m1, e1), (m2, e2) = (m2, e2), (m1, e1)
    e1 = e1 * np.exp(-1j * np.angle(e1[0]))
    e2 = e2 * np.exp(-1j * np.angle(e2[2]))
    Am = np.array([e1.real, e1.imag, e2.real, e2.imag]).T
    return Am, np.linalg.inv(Am), (m1, m2), None, cpl


# ---------------------------------------------------------------------------------------------
# 1 周のべき級数の写像
# ---------------------------------------------------------------------------------------------
def _run_ops(tab, ops, S, D, sd):
    """素子列 ops を TPSA の状態 S (4 成分。定数項を持ってよい) に順に作用させる"""
    for o in ops:
        if o[0] == 'lin':
            sg = o[1]
            new = []
            for i in range(4):
                acc = tab.const(0.0)
                for j in range(4):
                    ser = sg.M[:, i, j]
                    if np.any(ser):
                        acc = acc + tab.mul_dser(ser, S[j])
                dser = sg.d[:, i]
                if np.any(dser):
                    acc = acc + tab.mul_dser(dser, tab.const(1.0))
                new.append(acc)
            S = new
        elif o[0] == 'tmap':
            S = [np.asarray(a) for a in o[1](S, tab, sd)]
        elif o[0] == 'nlphi':
            Q5 = S + [sd * tab.var(4)]
            S = _compose(tab, o[1], Q5, D)
        else:
            _, k2l, k3l = o
            x, y = S[0], S[2]
            if k2l:
                x2 = tab.mul(x, x)
                y2 = tab.mul(y, y)
                xy = tab.mul(x, y)
                S[1] = S[1] - 0.5 * k2l * (x2 - y2)
                S[3] = S[3] + k2l * xy
            if k3l:
                x2 = tab.mul(x, x)
                y2 = tab.mul(y, y)
                x3 = tab.mul(x2, x)
                y3 = tab.mul(y2, y)
                xy2 = tab.mul(x, y2)
                x2y = tab.mul(x2, y)
                S[1] = S[1] - k3l / 6.0 * (x3 - 3 * xy2)
                S[3] = S[3] - k3l / 6.0 * (y3 - 3 * x2y)
    return S


def _eval_poly(tab, P, u):
    mono = np.prod(u[None, :] ** tab.mons, axis=1)
    return np.dot(P, mono)


def _fixed_point(tab, Fp, s, maxit=40, tol=1e-15):
    """物理座標の 1 周の写像 Fp (u = x/s の級数。δ̃ = 0) の不動点 c (m, rad) とそのまわりのヤコビ行列 A (4×4 実)。不動点の反復 (Newton 法、解析的な微分)"""
    mons = tab.mons
    c = np.zeros(4)

    def jac(u):
        A = np.zeros((4, 4), dtype=complex)
        for k in range(4):
            e = mons.copy()
            ek = e[:, k].copy()
            e[:, k] = np.maximum(ek - 1, 0)
            mono = np.prod(u[None, :] ** e, axis=1) * ek
            for j in range(4):
                A[j, k] = np.dot(Fp[j], mono) / s
        return A.real
    for _ in range(maxit):
        u = np.concatenate([c / s, [0.0]])
        f = np.array([_eval_poly(tab, Fp[j], u).real for j in range(4)])
        r = f - c
        if np.abs(r).max() < tol:
            break
        c = c - np.linalg.solve(jac(u) - np.eye(4), r)
    u = np.concatenate([c / s, [0.0]])
    return c, jac(u)


def one_turn_series(lat, D=5, nslice=8, sd=1e-2, s=3e-4, nsplit=None, ops=None, bend_cubic=False, nslice_bend=8, kin=False, fringe=False, nslice_quad=4, k0_nonlinear=None):
    """1 周の写像を 5 変数 (ζx, ζ̄x, ζy, ζ̄y, δ̃) のべき級数で返す (ζ = s ζ̃、δ = sd δ̃ の規格化)。戻り dict(G (4 成分の TPSA)、tab、mux、muy、
    twiss0 = (αx, βx, αy, βy)、coupling (線形の結合の最大値)、warnings)"""
    tab = get_tab(5, D)
    if ops is None:
        ops, warn = build_ops(lat, D, sd, nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=D, kin=kin, fringe=fringe, nslice_quad=nslice_quad, k0_nonlinear=k0_nonlinear)
    else:
        warn = []
    try:
        from . import ampnf_coupled as _ac
        if _ac.has_errors(lat):
            warn.append('DX・DY・DROTATE のある素子 (%s ...) の据え付けのずれは `one_turn_series` では入れていない (閉軌道 ≠ 0 の展開が要る。`amplitude_tune` は ampnf_coupled に回す。H1)' % _ac.has_errors(lat)[0])
    except Exception:
        pass
    has_t = any(o[0] == 'tmap' for o in ops)
    if has_t and D < 4:
        warn.append('K0≠0 の磁石があり D < 4: 軌道の傾きによる運動学 (4 次のハミルトニアン) のフィードダウンが 2 次の写像に入らない (D ≥ 4 を勧める。docs/49 §14.4)')
    if has_t:
        # v2.4.0 G: K0≠0 の偏向磁石がある (閉軌道が原点にない・写像に定数項)。第 1 段: 物理座標 (原点が入口) の 1 周の写像 F(u) を作り、不動点 c (閉軌道) と
        # そのまわりのヤコビ行列 A を求める。第 2 段: ζ = (モード) の正規化座標で、c を定数項に入れた状態から出発して全素子を通す (磁石の内部の軌道は閉じた式の展開で厳密)
        u = [tab.var(j) for j in range(5)]
        Fp = _run_ops(tab, ops, [s * u[j] for j in range(4)], D, sd)
        orbit, A = _fixed_point(tab, Fp, s)
        Am, Ai, mus, tw, cpl = linear_modes(A)
        (ax, bx, ay, by) = tw if tw is not None else (None,) * 4
        zz = [tab.var(j) for j in range(5)]
        zn = [0.5 * s * (zz[0] + zz[1]), 0.5j * s * (zz[0] - zz[1]), 0.5 * s * (zz[2] + zz[3]), 0.5j * s * (zz[2] - zz[3])]
        Sin = [tab.const(orbit[i]) + sum(Am[i, j] * zn[j] for j in range(4) if Am[i, j] != 0.0) for i in range(4)]
        if float(np.abs(orbit).max()) > 1e-12:
            # 閉軌道が原点にない: 閉軌道を定数項に入れた状態で全素子をもう一度通す (原点のまわりの展開を組み替えない = 次数 D の打ち切りの誤差が出ない)
            S = _run_ops(tab, ops, Sin, D, sd)
        else:
            S = _compose(tab, Fp, [Sin[i] / s for i in range(4)] + [tab.var(4)], D)
        S = [S[i] - tab.const(orbit[i]) for i in range(4)]
    else:
        orbit = np.zeros(4)
        # 線形の 1 周から Twiss (模型自身の)
        Mtot = Seg(D, sd)
        for o in ops:
            if o[0] == 'lin':
                Mtot = Mtot.then(o[1])
        # δ=0 の 1 周の行列には kick は寄与しない
        A = Mtot.M[0].real
        Am, Ai, mus, tw, cpl = linear_modes(A)
        (ax, bx, ay, by) = tw if tw is not None else (None,) * 4
        # 初期の状態: 物理座標を ζ の線形式で (正規化座標 zn = (x̂, p̂, ŷ, q̂)、x̂ = s Re ζx、p̂ = −s Im ζx。物理座標 = Am zn)
        z = [tab.var(j) for j in range(5)]
        def xhat(zp, zm):
            return 0.5 * s * (zp + zm)
        def phat(zp, zm):
            return 0.5j * s * (zp - zm)
        zn = [xhat(z[0], z[1]), phat(z[0], z[1]), xhat(z[2], z[3]), phat(z[2], z[3])]
        S = [sum(Am[i, j] * zn[j] for j in range(4) if Am[i, j] != 0.0) for i in range(4)]
        S = _run_ops(tab, ops, S, D, sd)
    # 出力を正規化座標に
    xo, po, yo, qo = [sum(Ai[i, j] * S[j] for j in range(4) if Ai[i, j] != 0.0) for i in range(4)]
    G = [(xo - 1j * po) / s, (xo + 1j * po) / s, (yo - 1j * qo) / s, (yo + 1j * qo) / s]
    return dict(G=G, tab=tab, mux=mus[0], muy=mus[1], twiss0=(ax, bx, ay, by), Amat=Am, Ainv=Ai, coupling=float(cpl), warnings=warn, s=s, sd=sd, D=D, orbit=orbit)


# ---------------------------------------------------------------------------------------------
# Birkhoff の正規形
# ---------------------------------------------------------------------------------------------
def _compose(tab, P_list, Q, nmax):
    """P_list の各 TPSA に Q (5 成分) を代入。次数 nmax 以下の単項式だけ使い、結果は次数 nmax まで有効"""
    N = tab.N
    pw = {}
    one = tab.const(1.0)
    pw[0] = one
    order = [i for i in range(N) if tab.deg[i] <= nmax]
    for i in order:
        if i == 0:
            continue
        e = tab.mons[i]
        j = int(np.nonzero(e)[0][0])
        e2 = list(e)
        e2[j] -= 1
        parent = tab.idx[tuple(e2)]
        pw[i] = tab.mul(Q[j], pw[parent])
    outs = []
    for P in P_list:
        acc = np.zeros(N, dtype=complex)
        for i in order:
            if P[i] != 0:
                acc += P[i] * pw[i]
        outs.append(acc)
    return outs


def normal_form(res, warn_den=1e-4, skip_den=0.0, keep=None):
    """res = one_turn_series の結果。戻り dict: Gn (正規形の写像 4 成分)、phi (変換)、min_den (共鳴の分母の最小値)、near (近い共鳴の項)。
    keep (v2.4.1): 消さずに残す単項式の集合 {(i, (e0, e1, e2, e3, e4))} (i = 方程式 0〜3、e = ζx, ζ̄x, ζy, ζ̄y, δ̃ のべき)。
    1 本の共鳴を残した正規形 (島の幅。`ampnf_ext.resonance_normal_form`)。残した係数 Gn[i][m] はその共鳴の駆動項 (低い次数の非共鳴の項の変換を含む)"""
    keep = set(keep) if keep else None
    tab = res['tab']
    D = res['D']
    mux, muy = res['mux'], res['muy']
    lam = np.array([np.exp(1j * mux), np.exp(-1j * mux), np.exp(1j * muy), np.exp(-1j * muy)])
    e = tab.mons
    phase = (e[:, 0] - e[:, 1]) * mux + (e[:, 2] - e[:, 3]) * muy
    lam_mon = np.exp(1j * phase)
    G = [g.copy() for g in res['G']]
    N = tab.N
    # 線形の部分は対角のはず
    zvars = [tab.var(j) for j in range(5)]
    Tm = [zvars[j].copy() for j in range(4)]      # 変換 T = w + φ
    Gp = [lam[i] * zvars[i] for i in range(4)]    # 正規形の写像 G' = Λw + ...
    min_den = np.inf
    near = []
    for n in range(1, D + 1):
        Q5 = Tm + [zvars[4]]
        GT = _compose(tab, G, Q5, n)
        G5 = Gp + [zvars[4]]
        TG = _compose(tab, Tm, G5, n)
        for i in range(4):
            R = GT[i] - TG[i]
            for m in range(N):
                if tab.deg[m] != n:
                    continue
                r = R[m]
                if r == 0:
                    continue
                d = lam_mon[m] - lam[i]
                if abs(d) < 1e-9 or abs(d) < skip_den or (keep is not None and (i, tuple(int(v) for v in tab.mons[m])) in keep):      # skip_den > 0: 分母が小さい (共鳴に近い) 項は消さずに残す (診断。ν の級数にはその項は入らない)
                    Gp[i][m] += r                     # 共鳴の項は残す (G' = R)
                else:
                    Tm[i][m] += r / d                 # φ = R/(λ^α − λ_i)
                    min_den = min(min_den, abs(d))
                    if abs(d) < warn_den:
                        near.append((i, tuple(int(v) for v in tab.mons[m]), abs(d)))
    return dict(Gn=Gp, phi=Tm, min_den=float(min_den), near=near, lam=lam)


def tune_series(res, nf, plane=0):
    """正規形から ν_plane の級数。戻り dict {(a, b, p): c_abp} (物理の単位: Jx^a Jy^b δ^p。J の単位 m·rad)、ν0 (端数)"""
    tab = res['tab']
    D = res['D']
    s, sd = res['s'], res['sd']
    Gn = nf['Gn'][2 * plane]
    lam0 = nf['lam'][2 * plane]
    # c(X, Y, δ̃) = G'_ζ/(λ ζ) の単項式 ζx^(a+1) ζ̄x^a ζy^b ζ̄y^b δ̃^p (plane 0)、ζx^a ζ̄x^a ζy^(b+1) ζ̄y^b δ̃^p (plane 1)
    c = {}
    for m in range(tab.N):
        if Gn[m] == 0:
            continue
        ex = tab.mons[m]
        if plane == 0 and ex[0] == ex[1] + 1 and ex[2] == ex[3]:
            c[(int(ex[1]), int(ex[2]), int(ex[4]))] = Gn[m] / lam0
        elif plane == 1 and ex[2] == ex[3] + 1 and ex[0] == ex[1]:
            c[(int(ex[0]), int(ex[3]), int(ex[4]))] = Gn[m] / lam0
    # c の係数 (c[(0,0,0)] − 1 = 0 を除く) を 1 + c として log の級数
    cc = {k: v for k, v in c.items()}
    cc[(0, 0, 0)] = cc.get((0, 0, 0), 0) - 1.0       # 線形の λζ の項 → 0
    maxw = D - 1                                      # 重み 2a+2b+p ≤ D−1

    def mulc(A, B):
        out = {}
        for ka, va in A.items():
            for kb, vb in B.items():
                k = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
                if 2 * k[0] + 2 * k[1] + k[2] <= maxw:
                    out[k] = out.get(k, 0) + va * vb
        return out
    cc = {k: v for k, v in cc.items() if 2 * k[0] + 2 * k[1] + k[2] <= maxw and abs(v) > 0}
    logs = {}
    powk = dict(cc)
    for kk in range(1, maxw + 1):
        for k, v in powk.items():
            logs[k] = logs.get(k, 0) + ((-1) ** (kk + 1)) * v / kk
        powk = mulc(powk, cc)
        if not powk:
            break
    out = {}
    for (a, b, p), v in logs.items():
        # X̃ = |ζ̃|² = 2J/s²、δ̃ = δ/sd
        ph = v.imag / (2 * math.pi)
        out[(a, b, p)] = ph * (2.0 / s ** 2) ** (a + b) / sd ** p
    return out


class AmpTuneResult:
    def __init__(self, nux, nuy, nu0, nf, res):
        self.nux, self.nuy, self.nu0, self.nf, self.res = nux, nuy, nu0, nf, res
        self.min_den = nf['min_den']
        self.near = nf['near']
        self.warnings = list(res['warnings'])

    def nu(self, plane, jx=0.0, jy=0.0, delta=0.0, jorder=None):
        """ν_plane(Jx, Jy, δ)。jorder (v2.4.1): J のべきの和 a+b ≤ jorder だけの級数 (1 = J¹ のみ = ADTS、2 = J² まで。比較・次数の違いの表示用)"""
        ser = self.nux if plane == 0 else self.nuy
        return self.nu0[plane] + sum(v * jx ** a * jy ** b * delta ** p for (a, b, p), v in ser.items() if jorder is None or a + b <= jorder)

    def d2nu(self, plane):
        """J の 2 次の ADTS (v2.4.1): ∂²ν_plane/∂Ji∂Jj (δ = 0) の 2×2 行列 [1/(m rad)²] (J² の係数 c_20・c_11・c_02 の ∂ 二階微分: [[2c_20, c_11], [c_11, 2c_02]])。D ≥ 5 で入る"""
        c20, c11, c02 = self.coef(plane, 2, 0, 0), self.coef(plane, 1, 1, 0), self.coef(plane, 0, 2, 0)
        return np.array([[2 * c20, c11], [c11, 2 * c02]])

    def coef(self, plane, a, b, p=0):
        return (self.nux if plane == 0 else self.nuy).get((a, b, p), 0.0)

    @property
    def adts(self):
        """J¹ の ADTS {'xx','xy','yx','yy'} [1/(m rad)] (∂νx/∂Jx = c_100 …)"""
        return {'xx': self.coef(0, 1, 0), 'xy': self.coef(0, 0, 1), 'yx': self.coef(1, 1, 0), 'yy': self.coef(1, 0, 1)}

    @property
    def adts2(self):
        """J² の項の係数 c_ab (∂²ν/∂J²/2 型。ν = … + c_20 Jx² + c_11 JxJy + c_02 Jy²) [1/(m rad)²]"""
        return {pl + k: self.coef(i, *ab) for i, pl in enumerate(('x', 'y')) for k, ab in (('20', (2, 0, 0)), ('11', (1, 1, 0)), ('02', (0, 2, 0)))}

    def chrom(self, plane, p=1):
        """∂ᵖν/∂δᵖ (p! c_00p)"""
        return math.factorial(p) * self.coef(plane, 0, 0, p)

    def adts_delta(self, plane, a, b, p=1):
        """∂ᵖ/∂δᵖ の ADTS (J^aJ^b の係数) p! c_abp"""
        return math.factorial(p) * self.coef(plane, a, b, p)


def amplitude_tune(lat, D=5, nslice=8, sd=1e-2, s=3e-4, ops=None, bend_cubic=False, nslice_bend=8, kin=False, fringe=False, nslice_quad=4, skip_den=0.0, k0_nonlinear=None):
    """lat の ν_x, ν_y(Jx, Jy, δ) の級数 (D 次の 1 周の写像の Birkhoff 正規形)。戻り AmpTuneResult"""
    if ops is None:
        # v2.4.0 G: DX・DY (・DZ・DROTATE) のある素子 (BEND・QUAD・SEXT・OCT) は、以前は黙って無視されていた (pysad 自身の不具合)。H1 の `ampnf_coupled` (tchge の TPSA 版、
        # 閉軌道 ≠ 0 の展開) に渡す。D は 4 以上、skip_den は未対応
        from . import ampnf_coupled as _ac
        if _ac.has_errors(lat):
            if k0_nonlinear is not None:
                raise NotImplementedError('ampnf: k0_nonlinear は DX・DY・DROTATE のある lat (ampnf_coupled の経路) では未対応')
            return _ac.amplitude_tune_ext(lat, D=D, nslice=nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad, sd=sd, s=s)
    res = one_turn_series(lat, D=D, nslice=nslice, sd=sd, s=s, ops=ops, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad, k0_nonlinear=k0_nonlinear)
    nf = normal_form(res, skip_den=skip_den)
    nux = tune_series(res, nf, 0)
    nuy = tune_series(res, nf, 1)
    nu0 = (res['mux'] / (2 * math.pi), res['muy'] / (2 * math.pi))
    return AmpTuneResult(nux, nuy, nu0, nf, res)


# ---------------------------------------------------------------------------------------------
# 共鳴駆動項 h_jklmp を 1 周の写像のべき級数 (2 次の部分) から (v2.4.0 F1: 結合した光学・ROTATE した偏向磁石)
# ---------------------------------------------------------------------------------------------
_PZ = np.array([[1, -1j, 0, 0], [1, 1j, 0, 0], [0, 0, 1, -1j], [0, 0, 1, 1j]], dtype=complex)       # ζ = P zn (zn = (x̂, p̂, ŷ, q̂)、ζ = x̂ − i p̂)


def _quad_mats(res):
    """1 周の写像の 2 次の部分 G_i^(2) (p=0) の対称行列 Q_i (ζ̃ の 4×4)、δ の 1 次の定数 g1、ζ̃δ の係数 d (4×4: d[i,j])、δ² の定数 g2 (δ は物理の値)"""
    tab, sd = res['tab'], res['sd']
    Q = np.zeros((4, 4, 4), dtype=complex)
    g1 = np.zeros(4, dtype=complex)
    g2 = np.zeros(4, dtype=complex)
    d = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        G = res['G'][i]
        for m in np.nonzero(G)[0]:
            e = tab.mons[m]
            deg = int(e.sum())
            ed = int(e[4])
            ez = [int(v) for v in e[:4]]
            if ed == 0 and deg == 2:
                idx = [a for a in range(4) for _ in range(ez[a])]
                a, b = idx
                if a == b:
                    Q[i, a, a] += G[m]
                else:
                    Q[i, a, b] += 0.5 * G[m]
                    Q[i, b, a] += 0.5 * G[m]
            elif ed == 1 and deg == 1:
                g1[i] += G[m] / sd
            elif ed == 2 and deg == 2:
                g2[i] += G[m] / sd ** 2
            elif ed == 1 and deg == 2:
                d[i, ez.index(1)] += G[m] / sd
    return Q, g1, g2, d


def rdt_series(res, with_delta=True):
    """1 周の写像の 2 次の部分 (D ≥ 2 の one_turn_series の結果) から、共鳴駆動項 h_jklmp (3 次: p=0 の 20 項、δ: p=1 の 8 項と h10002)。
    結合した線形光学 (4×4 の固有モード、linear_modes) でも同じ式 (モード変数 ζ1・ζ2 が基準)。方法は `rdttrack` と同じ (写像の 2 次の部分 → kick の生成関数
    → 角度・作用の調和)。**模型の写像の 2 次の部分は、6 極と偏向磁石の 3 次の項 (bend_cubic=True の曲率・端) の kick の和そのもの**で、追跡との違いは入れた効果だけ。
    δ の項は閉軌道のまわりのヤコビ行列の δ 微分と 2 次の閉軌道 (rdttrack.delta_quadratic_total・delta2_orbit_h10002 と同じ定義)。戻り {名前: 複素数}"""
    from . import rdttrack as T
    from . import bendnl as B
    s = res['s']
    P = _PZ
    Pi = np.linalg.inv(P)
    lam = np.array([np.exp(1j * res['mux']), np.exp(-1j * res['mux']), np.exp(1j * res['muy']), np.exp(-1j * res['muy'])])
    Q, g1, g2, d = _quad_mats(res)
    # 実の正規化座標での 2 次の部分: zn' = R zn + T(zn, zn)、T[j] = Σ_i Pi[j,i] P^T Q_i P / s
    Tn = np.einsum('ji,iab->jab', Pi, np.einsum('ca,icd,db->iab', P, Q, P)) / s
    Rr = (Pi @ np.diag(lam) @ P)
    rm = T.RingMap()
    rm.R, rm.T = Rr.real, Tn.real
    out = {}
    h3 = T.rdt_from_map(rm)
    out.update(h3)
    out['_imag'] = float(max(np.abs(Rr.imag).max(), np.abs(Tn.imag).max()))
    if with_delta:
        I4 = np.eye(4)
        L = np.diag(lam)
        z1 = np.linalg.solve(I4 - L, g1)                                  # ζ̃ の 1 次の閉軌道 (δ あたり)
        J1 = d + 2.0 * np.einsum('ijk,k->ij', Q, z1)                       # δ の 1 次のヤコビ行列 (ζ̃ 空間)
        Xr = (Pi @ J1 @ P).real
        Qm = -T.S4 @ np.linalg.solve(rm.R, Xr)
        out.update(T.quad_to_h(Qm, T.DELTA1_MONOS))
        z2 = np.linalg.solve(I4 - L, g2 + np.einsum('ijk,j,k->i', Q, z1, z1) + d @ z1)
        c2 = (Pi @ z2).real * s                                            # 実の正規化座標の δ² の係数
        gradg = -T.S4 @ np.linalg.solve(rm.R, (I4 - rm.R) @ c2)
        amp = B.harmonic_amplitudes({(1, 0, 0, 0): gradg[0], (0, 1, 0, 0): gradg[1], (0, 0, 1, 0): gradg[2], (0, 0, 0, 1): gradg[3]}, T.W_NORM, [(1, 0, 0, 0)])[0, 0]
        out['h10002'] = complex(amp)
    return out


# ---------------------------------------------------------------------------------------------
# 同じ模型の数値の追跡 (検証用): 固定した δ で 1 周の写像 (線形 + kick) を数値で作る
# ---------------------------------------------------------------------------------------------
def model_tracker(lat, delta=0.0, nslice=8, ops=None, bend_cubic=False, nslice_bend=8, Dnl=5, kin=False, fringe=False, nslice_quad=4):
    """模型の 1 周の写像 v(n,4) → v (δ 固定、厳密な線形の行列 + 薄い 6 極・8 極の kick)。ops を省略すると build_ops(delta=δ)"""
    if ops is None:
        ops, _ = build_ops(lat, 0, 1.0, nslice, delta=delta, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=Dnl, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    tabn = get_tab(5, Dnl)

    def tr(v):
        v = np.array(v, dtype=float)
        for o in ops:
            if o[0] == 'lin':
                v = v @ o[1].M.T + o[1].d
            elif o[0] == 'tmap':
                from . import ampnf_k0 as _k0
                if len(o) > 3 and o[3]:
                    raise NotImplementedError('model_tracker: K0≠0 の磁石の線形だけの写像 (bend_cubic・kin・fringe を全部切った模型) の数値の追跡は未対応')
                v = _k0.numeric_map(o[2], delta)(v)
            elif o[0] == 'nlphi':
                full = np.column_stack([v, np.full(len(v), delta)])
                mono = np.prod(full[:, None, :] ** tabn.mons[None, :, :], axis=2)
                v = np.column_stack([(mono * o[1][j].real[None, :]).sum(axis=1) for j in range(4)])
            else:
                _, k2l, k3l = o
                x, y = v[:, 0], v[:, 2]
                if k2l:
                    v[:, 1] -= 0.5 * k2l * (x * x - y * y)
                    v[:, 3] += k2l * x * y
                if k3l:
                    v[:, 1] -= k3l / 6.0 * (x ** 3 - 3 * x * y * y)
                    v[:, 3] -= k3l / 6.0 * (y ** 3 - 3 * x * x * y)
        return v
    return tr


def model_closed_orbit(tr, tol=1e-14):
    z = np.zeros((1, 4))
    for _ in range(30):
        f = tr(z) - z
        J = np.zeros((4, 4))
        for j in range(4):
            e = np.zeros((1, 4))
            e[0, j] = 1e-7
            J[:, j] = ((tr(z + e) - tr(z - e)) / 2e-7)[0] - 0 * z[0]
        J -= np.eye(4) * 0
        A = (J - np.eye(4))
        dz = np.linalg.solve(A, -f[0])
        z = z + dz
        if np.abs(dz).max() < tol:
            break
    return z[0]


def model_tune(lat, jx, jy, delta=0.0, nturn=1024, nslice=8, window=1, phase=(0.3, 0.9), ops=None, nfres=None, bend_cubic=False, nslice_bend=8, Dnl=5, kin=False, fringe=False, nslice_quad=4):
    """模型の追跡 + NAFF の (νx, νy) の端数。リング先頭の Twiss と閉軌道から、Jx, Jy (m rad) の粒子。ν0 に近い側 (線形の端数) に揃えて返す"""
    from . import naff
    tr = model_tracker(lat, delta, nslice, ops=ops, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=Dnl, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    # Twiss (δ=0 の模型の線形の 1 周): ops の線形の積から
    ops0, _ = build_ops(lat, 0, 1.0, nslice, delta=0.0, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=Dnl, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    Mt = NSeg()
    for o in ops0:
        if o[0] == 'lin':
            Mt = Mt.then(o[1])
    A = Mt.M
    orb = model_closed_orbit(tr)
    if any(o[0] == 'tmap' for o in ops0):
        from .bendnl import jac_at
        A = jac_at(lambda Z: tr(Z), orb)         # v2.4.0 G: K0≠0 の磁石があるときは模型の追跡の閉軌道まわりのヤコビ行列
    Am, Ai, mus, tw, cpl = linear_modes(A)
    ph, pv = phase
    if nfres is None:
        # アクション・位相 → 物理座標 (正規化座標 x̂ = √(2J) cos φ、p̂ = −√(2J) sin φ)
        zn = np.array([math.sqrt(2 * jx) * math.cos(ph), -math.sqrt(2 * jx) * math.sin(ph), math.sqrt(2 * jy) * math.cos(pv), -math.sqrt(2 * jy) * math.sin(pv)])
        v = np.array([Am @ zn]) + orb
    else:
        v = np.array([initial_from_nf(nfres, jx, jy, ph, pv, delta)])
    rec = np.zeros((nturn, 4))
    for t in range(nturn):
        rec[t] = v[0]
        v = tr(v)
    zrec = (Ai @ (rec - orb).T).T
    out = []
    for p in range(2):
        xh, ph_ = zrec[:, 2 * p], zrec[:, 2 * p + 1]
        zz = (xh - xh.mean()) - 1j * (ph_ - ph_.mean())
        f0 = (mus[p] / (2 * math.pi) + 0.5) % 1.0 - 0.5
        f, _a = naff.naff(zz, 1, window=window, fmin=f0 - 0.05, fmax=f0 + 0.05)
        out.append(f[0] % 1.0)
    return tuple(out)


def initial_from_nf(r, jx, jy, phx, phy, delta):
    """正規形の座標 (J, φ, δ) に対応する物理の初期座標 (x, px, y, py)。変換 T (正規形の計算で得た φ) で戻す (T は次数 D まで)。
    r = AmpTuneResult。J・φ から解析式の ν(J, δ) に対応する粒子を作る (追跡の初期条件の J が正規形の J と O(√J) 違う影響を除く)"""
    res, nf = r.res, r.nf
    tab, s, sd = res['tab'], res['s'], res['sd']
    Am = res['Amat']
    w = np.array([math.sqrt(2 * jx) / s * np.exp(1j * phx), math.sqrt(2 * jx) / s * np.exp(-1j * phx),
                  math.sqrt(2 * jy) / s * np.exp(1j * phy), math.sqrt(2 * jy) / s * np.exp(-1j * phy), delta / sd])
    mono = np.prod(w[None, :] ** tab.mons, axis=1)
    z = [s * np.dot(nf['phi'][i], mono) for i in range(4)]
    zn = np.array([0.5 * (z[0] + z[1]).real, (0.5j * (z[0] - z[1])).real, 0.5 * (z[2] + z[3]).real, (0.5j * (z[2] - z[3])).real])
    return list(Am @ zn + res.get('orbit', 0.0))


def nf_action(r, coords, delta=0.0, niter=60):
    """物理の座標 (x, px, y, py) (リング先頭) に対応する正規形のアクション (Jx, Jy) [m rad] と位相。T の逆を不動点の反復 w = z − φ(w) で求める (T は次数 D まで)。
    初期条件の J と正規形の J が O(√J) 違うことを補正して、追跡のチューンと ν(J, δ) の級数を比べるために使う"""
    res, nf = r.res, r.nf
    tab, s, sd = res['tab'], res['s'], res['sd']
    zn = res['Ainv'] @ (np.asarray(coords, dtype=float) - res.get('orbit', 0.0))
    xh, ph, yh, qh = zn
    z = np.array([(xh - 1j * ph) / s, (xh + 1j * ph) / s, (yh - 1j * qh) / s, (yh + 1j * qh) / s])
    w = z.copy()
    for _ in range(niter):
        full = np.concatenate([w, [delta / sd]])
        mono = np.prod(full[None, :] ** tab.mons, axis=1)
        zz = np.array([np.dot(nf['phi'][i], mono) for i in range(4)])
        w = w - (zz - z)
    jx = abs(w[0]) ** 2 * s * s / 2
    jy = abs(w[2]) ** 2 * s * s / 2
    return jx, jy, np.angle(w[0]), np.angle(w[2]), float(np.abs(zz - z).max())


# ---------------------------------------------------------------------------------------------
# 検証用: 実際の追跡 (pysad の追跡、全非線形・全素子) + NAFF の ADTS (結合した光学でも使える。v2.4.0 F1)
# ---------------------------------------------------------------------------------------------
def tracked_tunes(lat, r, jx, jy, phase=(0.3, 0.9), delta=0.0, nturn=1024, momentum_eV=3e9):
    """AmpTuneResult r の正規形の座標 (Jx, Jy, φ) から出発した、pysad の追跡 (k64・v2 の track6n。RF・放射なし) + NAFF の (νx, νy)。
    ν は r.nu0 に近い側の端数。結合した光学では、モード ζ1・ζ2 の正規化座標 (r.res['Ainv']) の NAFF"""
    from . import adtsdirect as D, naff
    c0 = initial_from_nf(r, jx, jy, phase[0], phase[1], delta)
    c = np.zeros((1, 6))
    c[0, :4] = c0
    c[0, 5] = delta
    rec = D.track_record(lat.elements, momentum_eV, c, nturn)
    from . import rdttrack as T
    orb = r.res.get('orbit', np.zeros(4)) if delta == 0.0 else T.closed_orbit(lat.elements, momentum_eV, delta)[0]
    zr = (r.res['Ainv'] @ (rec[:, 0, :4] - orb).T).T
    out = []
    for p in range(2):
        xh, ph = zr[:, 2 * p], zr[:, 2 * p + 1]
        zz = (xh - xh.mean()) - 1j * (ph - ph.mean())
        f0 = (r.nu0[p] + 0.5) % 1.0 - 0.5
        f, _ = naff.naff(zz, 1, window=1, fmin=f0 - 0.05, fmax=f0 + 0.05)
        out.append(f[0] % 1.0)
    return tuple(out)


def tracked_adts(lat, r, j=3e-7, nturn=1024, momentum_eV=3e9):
    """追跡 + NAFF の J¹ の ADTS {'xx','xy','yx','yy'}: 正規形の J の格子 (j, 2j, 3j の組 8 点。両方の面に振幅) の ν を、ν0 + c10 Jx + c01 Jy + c20 Jx² + c11 JxJy + c02 Jy² で最小二乗フィット
    (初期条件を正規形の J・φ から作るので、O(√J) の J のずれの影響を受けない)"""
    pts = [(j, j), (2 * j, j), (j, 2 * j), (2 * j, 2 * j), (3 * j, j), (j, 3 * j), (3 * j, 2 * j), (2 * j, 3 * j)]       # 両方の面に振幅 (振幅 0 の面の NAFF は雑音を拾う)
    nus = np.array([tracked_tunes(lat, r, a, b, nturn=nturn, momentum_eV=momentum_eV) for a, b in pts])
    Amat = np.array([[a, b, a * a, a * b, b * b] for a, b in pts])
    out = {}
    for pl, nm in ((0, 'x'), (1, 'y')):
        y = nus[:, pl] - r.nu0[pl] % 1.0
        y = (y + 0.5) % 1.0 - 0.5
        co, *_ = np.linalg.lstsq(Amat, y, rcond=None)
        out[nm + 'x'] = co[0]
        out[nm + 'y'] = co[1]
    return out
