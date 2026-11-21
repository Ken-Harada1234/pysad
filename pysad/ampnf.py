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
【入れていない効果】 K0≠0・ROTATE の偏向磁石 (NotImplementedError)、結合、誤差、K1 の偏向磁石の 4 次以上、K64 の四極の積分器の台形則、10 極以上、放射・RF。
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
        return np.tile(A.astype(complex), (n, 1, 1)), np.zeros((n, 4), dtype=complex), 'other'
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


def build_ops(lat, D, sd, nslice=8, delta=None, bend_cubic=False, nslice_bend=8, Dnl=5, kin=False, fringe=False, nslice_quad=4):
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
    for el in lat.elements:
        tn = type(el).__name__
        L = float(getattr(el, 'L', 0.0) or 0.0)
        k2l = k3l = 0.0
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


# ---------------------------------------------------------------------------------------------
# 1 周のべき級数の写像
# ---------------------------------------------------------------------------------------------
def one_turn_series(lat, D=5, nslice=8, sd=1e-2, s=3e-4, nsplit=None, ops=None, bend_cubic=False, nslice_bend=8, kin=False, fringe=False, nslice_quad=4):
    """1 周の写像を 5 変数 (ζx, ζ̄x, ζy, ζ̄y, δ̃) のべき級数で返す (ζ = s ζ̃、δ = sd δ̃ の規格化)。戻り dict(G (4 成分の TPSA)、tab、mux、muy、
    twiss0 = (αx, βx, αy, βy)、coupling (線形の結合の最大値)、warnings)"""
    tab = get_tab(5, D)
    if ops is None:
        ops, warn = build_ops(lat, D, sd, nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=D, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    else:
        warn = []
    # 線形の 1 周から Twiss (模型自身の)
    Mtot = Seg(D, sd)
    for o in ops:
        if o[0] == 'lin':
            Mtot = Mtot.then(o[1])
    # δ=0 の 1 周の行列には kick は寄与しない
    A = Mtot.M[0].real
    cpl = max(np.abs(A[:2, 2:]).max(), np.abs(A[2:, :2]).max())
    tw = []
    mus = []
    for b in (0, 2):
        m = A[b:b + 2, b:b + 2]
        cm = 0.5 * (m[0, 0] + m[1, 1])
        sm = math.copysign(math.sqrt(max(1 - cm * cm, 0.0)), m[0, 1])
        mu = math.atan2(sm, cm)
        beta = m[0, 1] / sm
        alpha = (m[0, 0] - m[1, 1]) / (2 * sm)
        tw += [alpha, beta]
        mus.append(mu % (2 * math.pi))
    (ax, bx, ay, by) = tw
    # 初期の状態: 物理座標を ζ の線形式で
    z = [tab.var(j) for j in range(5)]
    def xhat(zp, zm):
        return 0.5 * s * (zp + zm)
    def phat(zp, zm):
        return 0.5j * s * (zp - zm)
    xh, ph = xhat(z[0], z[1]), phat(z[0], z[1])
    yh, qh = xhat(z[2], z[3]), phat(z[2], z[3])
    S = [math.sqrt(bx) * xh, (ph - ax * xh) / math.sqrt(bx), math.sqrt(by) * yh, (qh - ay * yh) / math.sqrt(by)]
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
    # 出力を正規化座標に
    xo = S[0] * (1 / math.sqrt(bx))
    po = S[1] * math.sqrt(bx) + ax * xo
    yo = S[2] * (1 / math.sqrt(by))
    qo = S[3] * math.sqrt(by) + ay * yo
    G = [(xo - 1j * po) / s, (xo + 1j * po) / s, (yo - 1j * qo) / s, (yo + 1j * qo) / s]
    return dict(G=G, tab=tab, mux=mus[0], muy=mus[1], twiss0=(ax, bx, ay, by), coupling=float(cpl), warnings=warn, s=s, sd=sd, D=D)


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


def normal_form(res, warn_den=1e-4, skip_den=0.0):
    """res = one_turn_series の結果。戻り dict: Gn (正規形の写像 4 成分)、phi (変換)、min_den (共鳴の分母の最小値)、near (近い共鳴の項)"""
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
                if abs(d) < 1e-9 or abs(d) < skip_den:      # skip_den > 0: 分母が小さい (共鳴に近い) 項は消さずに残す (診断。ν の級数にはその項は入らない)
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

    def nu(self, plane, jx=0.0, jy=0.0, delta=0.0):
        ser = self.nux if plane == 0 else self.nuy
        return self.nu0[plane] + sum(v * jx ** a * jy ** b * delta ** p for (a, b, p), v in ser.items())

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


def amplitude_tune(lat, D=5, nslice=8, sd=1e-2, s=3e-4, ops=None, bend_cubic=False, nslice_bend=8, kin=False, fringe=False, nslice_quad=4, skip_den=0.0):
    """lat の ν_x, ν_y(Jx, Jy, δ) の級数 (D 次の 1 周の写像の Birkhoff 正規形)。戻り AmpTuneResult"""
    res = one_turn_series(lat, D=D, nslice=nslice, sd=sd, s=s, ops=ops, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    nf = normal_form(res, skip_den=skip_den)
    nux = tune_series(res, nf, 0)
    nuy = tune_series(res, nf, 1)
    nu0 = (res['mux'] / (2 * math.pi), res['muy'] / (2 * math.pi))
    return AmpTuneResult(nux, nuy, nu0, nf, res)


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
    tw = []
    for b in (0, 2):
        m = A[b:b + 2, b:b + 2]
        cm = 0.5 * (m[0, 0] + m[1, 1])
        sm = math.copysign(math.sqrt(1 - cm * cm), m[0, 1])
        tw.append(((m[0, 0] - m[1, 1]) / (2 * sm), m[0, 1] / sm))
    (ax, bx), (ay, by) = tw
    orb = model_closed_orbit(tr)
    ph, pv = phase
    if nfres is None:
        v = np.array([[math.sqrt(2 * jx * bx) * math.cos(ph), -math.sqrt(2 * jx / bx) * (math.sin(ph) + ax * math.cos(ph)),
                       math.sqrt(2 * jy * by) * math.cos(pv), -math.sqrt(2 * jy / by) * (math.sin(pv) + ay * math.cos(pv))]]) + orb
    else:
        v = np.array([initial_from_nf(nfres, jx, jy, ph, pv, delta)])
    rec = np.zeros((nturn, 4))
    for t in range(nturn):
        rec[t] = v[0]
        v = tr(v)
    out = []
    for p, (al, be) in enumerate(((ax, bx), (ay, by))):
        x = rec[:, 2 * p] - orb[2 * p]
        px = rec[:, 2 * p + 1] - orb[2 * p + 1]
        xh, ph_ = naff.twiss_normalize(x, px, al, be)
        zz = (xh - xh.mean()) - 1j * (ph_ - ph_.mean())
        cm = 0.5 * (A[2 * p, 2 * p] + A[2 * p + 1, 2 * p + 1])
        sm = math.copysign(math.sqrt(max(1 - cm * cm, 0.0)), A[2 * p, 2 * p + 1])
        f0 = (math.atan2(sm, cm) / (2 * math.pi) + 0.5) % 1.0 - 0.5
        f, _a = naff.naff(zz, 1, window=window, fmin=f0 - 0.05, fmax=f0 + 0.05)
        out.append(f[0] % 1.0)
    return tuple(out)


def initial_from_nf(r, jx, jy, phx, phy, delta):
    """正規形の座標 (J, φ, δ) に対応する物理の初期座標 (x, px, y, py)。変換 T (正規形の計算で得た φ) で戻す (T は次数 D まで)。
    r = AmpTuneResult。J・φ から解析式の ν(J, δ) に対応する粒子を作る (追跡の初期条件の J が正規形の J と O(√J) 違う影響を除く)"""
    res, nf = r.res, r.nf
    tab, s, sd = res['tab'], res['s'], res['sd']
    ax, bx, ay, by = res['twiss0']
    w = np.array([math.sqrt(2 * jx) / s * np.exp(1j * phx), math.sqrt(2 * jx) / s * np.exp(-1j * phx),
                  math.sqrt(2 * jy) / s * np.exp(1j * phy), math.sqrt(2 * jy) / s * np.exp(-1j * phy), delta / sd])
    mono = np.prod(w[None, :] ** tab.mons, axis=1)
    z = [s * np.dot(nf['phi'][i], mono) for i in range(4)]
    xh = 0.5 * (z[0] + z[1]).real
    ph = (0.5j * (z[0] - z[1])).real
    yh = 0.5 * (z[2] + z[3]).real
    qh = (0.5j * (z[2] - z[3])).real
    return [math.sqrt(bx) * xh, (ph - ax * xh) / math.sqrt(bx), math.sqrt(by) * yh, (qh - ay * yh) / math.sqrt(by)]


def nf_action(r, coords, delta=0.0, niter=60):
    """物理の座標 (x, px, y, py) (リング先頭) に対応する正規形のアクション (Jx, Jy) [m rad] と位相。T の逆を不動点の反復 w = z − φ(w) で求める (T は次数 D まで)。
    初期条件の J と正規形の J が O(√J) 違うことを補正して、追跡のチューンと ν(J, δ) の級数を比べるために使う"""
    res, nf = r.res, r.nf
    tab, s, sd = res['tab'], res['s'], res['sd']
    ax, bx, ay, by = res['twiss0']
    x, px, y, py = coords
    xh = x / math.sqrt(bx)
    ph = math.sqrt(bx) * px + ax * xh
    yh = y / math.sqrt(by)
    qh = math.sqrt(by) * py + ay * yh
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
