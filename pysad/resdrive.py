# -*- coding: utf-8 -*-
"""6 極 (と 8 極) による共鳴駆動項 h_jklmp、振幅依存チューンシフト (ADTS)、クロマティシティの高次項 (v2.3.0 Q2、K64 にない機能)。

解説と式は docs/49_resonance_terms_adts.md。ここは実装の要約。

ハミルトニアン (s 方向、パラキシャル): H = ½(px²+py²) + Σ_n (1/(n+1)!) Re[(K_n + i S_n)(x+iy)^(n+1)] (n=2 が 6 極、n=3 が 8 極。
K_n は SAD の KnL (積分値) を長さで割った密度。SAD の符号の向きは pysad の追跡の kick で確認済み: Δpx = −(K2L/2)(x²−y²)、Δpy = +K2L·x·y)。
アクション・角度 x = √(2Jx) Re[w e^{iφx}]、w(s) = √β e^{i(μ(s)−μ_ref)}、δ に比例する分散 x → x + η δ。

  h_jklmp = Σ_素子 ∫ds (H の係数)  :  H = Σ h_jklmp (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) δ^p exp{i[(j−k)φx + (l−m)φy]}
  1 周のハミルトニアン (1 次) は h の全周の和。厚い素子の積分 ∫ w^j w̄^k w_y^l w̄_y^m η^r ds は
    method='analytic' : 素子内の w(s)・η(s) を指数・多項式の和で表して厳密に積分 (ドリフト = 多項式、四極成分つきの MULT = 指数)。
    method='numeric'  : 素子を nslice 等分して中点で評価 (= nslice 個の薄レンズ列)。
  長さ 0 の素子 (薄レンズ) では積分を素子の位置での値に置き換える (厚い式の L→0 の極限)。

  ADTS (6 極の 2 次): K2(J) = −½ Σ_k Σ_(q,q') D·S,  S = (1/(4π sin πx_k)) ∬ ds1 ds2 F_q(s1) F̄_q'(s2) [e^{−iπx_k} (s1>s2) + e^{+iπx_k} (s1<s2)],
  x_k = k·Q。∂Q/∂J = ∂²K/∂J∂J。8 極は 1 次: K1 = (1/2π)∮⟨H4⟩ds。
  v2.3.0 R2: 運動学的な項 p⁴/8 と四極の端部のフリンジの 1 次 (kinematic_adts・quad_fringe_adts)。
  v2.3.0 R3: 偏向磁石の 3 次の項 (曲率 (h/2)x(px²+py²)・硬い端 ±(h/2)px y²・複合機能の本体・E1/E2≠0 の端) の 2 次の対の項と 6 極との交差項
  (adts(bend_cubic=True)、係数は pysad/bendnl.py が K64 の追跡の写像から取り出す。bend_cubic_adts)。
"""
import math
from math import comb
import numpy as np

TWO_PI = 2.0 * math.pi
_KAPPA_L_MIN = 1e-3          # |κL| がこれ未満の素子内の四極成分はドリフトとして扱う (誤差 ~ (κL)²)


# ---------------------------------------------------------------------------------------------
# 指数・多項式の和 Σ c s^n exp(λ s)。キー (n, cx, cy): λ = cx·λx + cy·λy (λx, λy は素子ごと)。ドリフトでは λ=0 (キーは (n,0,0))
# ---------------------------------------------------------------------------------------------
class _EP:
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = d if d is not None else {}

    @staticmethod
    def const(c):
        return _EP({(0, 0, 0): complex(c)})

    def __mul__(self, o):
        if not isinstance(o, _EP):
            return _EP({k: v * o for k, v in self.d.items()})
        r = {}
        for (n1, a1, b1), c1 in self.d.items():
            for (n2, a2, b2), c2 in o.d.items():
                k = (n1 + n2, a1 + a2, b1 + b2)
                r[k] = r.get(k, 0j) + c1 * c2
        return _EP(r)

    def __add__(self, o):
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = r.get(k, 0j) + v
        return _EP(r)

    def conj(self, flipx, flipy):
        return _EP({(n, -a if flipx else a, -b if flipy else b): np.conj(v) for (n, a, b), v in self.d.items()})

    def integral(self, L, rate):
        """∫0^L"""
        s = 0j
        for (n, a, b), c in self.d.items():
            lam = rate(a, b)
            if lam == 0:
                s += c * L ** (n + 1) / (n + 1)
            else:
                if n != 0:
                    raise RuntimeError('resdrive: s^n exp(λs) with λ≠0, n>0 は未対応')
                s += c * np.expm1(lam * L) / lam
        return s

    def anti(self, rate):
        """原始関数 F(s) = ∫0^s"""
        r = {}
        for (n, a, b), c in self.d.items():
            lam = rate(a, b)
            if lam == 0:
                k = (n + 1, a, b)
                r[k] = r.get(k, 0j) + c / (n + 1)
            else:
                if n != 0:
                    raise RuntimeError('resdrive: s^n exp(λs) with λ≠0, n>0 は未対応')
                r[(0, a, b)] = r.get((0, a, b), 0j) + c / lam
                r[(0, 0, 0)] = r.get((0, 0, 0), 0j) - c / lam
        return _EP(r)

    def value(self, t, rate):
        t = np.asarray(t, dtype=float)
        out = np.zeros(t.shape, dtype=complex)
        for (n, a, b), c in self.d.items():
            out = out + c * t ** n * np.exp(rate(a, b) * t)
        return out


def _term_names(t):
    return 'h%d%d%d%d%d' % t if t[4] < 10 else 'h%d%d%d%d_%d' % t


def parse_term(name):
    """'h21000' → (2,1,0,0,0)、'h11001' → (1,1,0,0,1)。5 文字目まで (j,k,l,m,p)。タプルはそのまま通す"""
    if isinstance(name, (tuple, list)):
        t = tuple(int(v) for v in name)
        return t if len(t) == 5 else t + (0,) * (5 - len(t))
    s = str(name).lower().lstrip('h')
    if len(s) != 5 or not s.isdigit():
        raise ValueError('resdrive: 項の名前は h + 5 桁 (jklmp) 例 h21000: %r' % (name,))
    return tuple(int(c) for c in s)


DEFAULT_TERMS = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200', 'h11001', 'h00111', 'h20001', 'h00201', 'h10002')


# ---------------------------------------------------------------------------------------------
# 素子の抽出
# ---------------------------------------------------------------------------------------------
class _Elem:
    pass


def _tn(el):
    return type(el).__name__


def collect_nonlinear(elements, rows, delta=0.0, ref_phase=(0.0, 0.0), warnings=None):
    """多極成分をもつ素子の一覧。rows は Lattice.twiss_table() の行 (rows[i] = i 番目の素子の入口)。
    強さは積分値 (薄レンズでも厚くても同じ意味)。戻り: list[_Elem]"""
    warnings = warnings if warnings is not None else []
    sc = 1.0 / (1.0 + delta)         # 追跡の素子の kick は正準運動量 px (p0 で規格化) に対する場の強さそのもので、δ に依らない (v2.3.0 R6: gs と pysad の追跡で確認。多極の強さは scs=1)。MULT の K1 だけ内部の光学のために 1/(1+δ)
    scs = 1.0
    out = []
    skipped = {}
    for i, el in enumerate(elements):
        tn = _tn(el)
        L = float(getattr(el, 'L', 0.0) or 0.0)
        strengths = {}
        k1 = 0.0
        if tn == 'Sextupole':
            k = float(el.K2)
            strengths[2] = complex(k * L if L > 0 else k)
            if L > 0 and k != 0.0:
                # v2.3.0 R6: 軌道がずれた (δ≠0 の分散軌道・閉軌道) 6 極は、内部で 4 極 K1 = K2·x0 (フィードダウン) として働く。内部の光学に入れる (x0 = 素子内の軌道の平均)
                x0 = 0.5 * (float(getattr(rows[i], 'x_cod', 0.0) or 0.0) + float(getattr(rows[i + 1], 'x_cod', 0.0) or 0.0))
                k1 = k * x0
        elif tn == 'Octupole':
            k = float(el.K3)
            strengths[3] = complex(k * L if L > 0 else k)
            if k != 0.0:
                # v2.3.0 R6: 軌道がずれた 8 極のフィードダウン: H4 = (K3/24)(x⁴ − 6x²y² + y⁴) を x → x0 + x で展開した 3 次の項は 6 極 K2 = K3·x0
                x0 = 0.5 * (float(getattr(rows[i], 'x_cod', 0.0) or 0.0) + float(getattr(rows[i + 1], 'x_cod', 0.0) or 0.0))
                if x0 != 0.0:
                    strengths[2] = complex(k * x0 * (L if L > 0 else 1.0))
        elif tn == 'Multipole':
            knl = {int(n): float(v) for n, v in (getattr(el, 'KnL', {}) or {}).items() if v}
            ksl = {int(n): float(v) for n, v in (getattr(el, 'KsL', {}) or {}).items() if v}
            for n in (2, 3):
                if n in knl or n in ksl:
                    strengths[n] = complex(knl.get(n, 0.0), ksl.get(n, 0.0))
            if 1 in knl and L > 0:
                k1 = knl[1] / L
            if ksl.get(1, 0.0) or float(getattr(el, 'ANGLE', 0.0) or 0.0) or knl.get(0, 0.0) or ksl.get(0, 0.0):
                if strengths:
                    warnings.append('MULT %s: K0/SK0/SK1/ANGLE は内部の光学に入れていない (rows の光学を入口の値に使うが、内部は K1 のみ)' % el.name)
            if any(n >= 4 for n in list(knl) + list(ksl)):
                warnings.append('MULT %s: 4 次以上 (10 極以上) の成分は無視' % el.name)
        elif tn in ('Decapole', 'Dodecapole'):
            if getattr(el, 'K4', 0.0) or getattr(el, 'K5', 0.0):
                skipped[tn] = skipped.get(tn, 0) + 1
            continue
        else:
            continue
        if not strengths:
            continue
        rot = float(getattr(el, 'ROTATION', 0.0) or 0.0) + float(getattr(el, 'DROTATE', 0.0) or 0.0)
        if rot != 0.0:
            for n in list(strengths):
                strengths[n] *= complex(math.cos((n + 1) * rot), math.sin((n + 1) * rot))      # pysad の ROTATION の向き (追跡の kick と照合済み): c → c e^{+iNθ}
            if k1 != 0.0:
                warnings.append('MULT %s: 回転した素子の K1 は光学に入れていない' % el.name)
                k1 = 0.0
        if float(getattr(el, 'DX', 0.0) or 0.0) or float(getattr(el, 'DY', 0.0) or 0.0):
            warnings.append('素子 %s: DX/DY (据え付けのずれ) は未対応 (誤差なしの設計光学を想定)' % el.name)
        e = _Elem()
        e.idx, e.name, e.L = i, el.name, L
        e.strengths = {n: c * scs for n, c in strengths.items()}
        e.k1 = k1 * sc
        e.fl = sc                    # v2.3.0 R6: 素子内のドリフトの有効長 L/(1+δ) (追跡は x += px·s/(1+δ)。K64 の6 極の厚い本体)
        r0, r1 = rows[i], rows[i + 1]
        e.opt0 = (float(r0.betx), float(r0.alfx), float(r0.mux) - ref_phase[0], float(r0.bety), float(r0.alfy), float(r0.muy) - ref_phase[1],
                  float(r0.etx), float(r0.etpx))
        e.opt1 = (float(r1.betx), float(r1.alfx), float(r1.mux) - ref_phase[0], float(r1.bety), float(r1.alfy), float(r1.muy) - ref_phase[1],
                  float(r1.etx), float(r1.etpx))
        e.s0 = float(rows[i].s)
        kap = math.sqrt(abs(e.k1)) if L > 0 else 0.0
        if L > 0 and kap * L >= _KAPPA_L_MIN:
            e.kind = 'xfoc' if e.k1 > 0 else 'xdef'
            e.kappa = kap
        else:
            e.kind = 'drift'
            e.kappa = 0.0
        out.append(e)
    for k, v in skipped.items():
        warnings.append('%s %d 台: 10 極以上は未対応 (無視)' % (k, v))
    return out


def _rates(e):
    """(rate(cx,cy), flipx, flipy)"""
    if e.kind == 'drift':
        return (lambda a, b: 0.0), False, False
    k = e.kappa
    if e.kind == 'xfoc':
        lx, ly = 1j * k, k + 0j
        return (lambda a, b: a * lx + b * ly), True, False
    lx, ly = k + 0j, 1j * k
    return (lambda a, b: a * lx + b * ly), False, True


def _m_ep(e, plane):
    """(m11, m12) を _EP で。plane 0 = x (kx = k1)、1 = y (ky = −k1)"""
    if e.kind == 'drift':
        return _EP({(0, 0, 0): 1.0 + 0j}), _EP({(1, 0, 0): float(getattr(e, 'fl', 1.0)) + 0j})
    k = e.kappa
    u = (1, 0) if plane == 0 else (0, 1)
    kp = e.k1 if plane == 0 else -e.k1
    p, m = (0, u[0], u[1]), (0, -u[0], -u[1])
    m11 = _EP({p: 0.5 + 0j, m: 0.5 + 0j})
    if kp > 0:
        m12 = _EP({p: 1.0 / (2j * k), m: -1.0 / (2j * k)})
    else:
        m12 = _EP({p: 1.0 / (2 * k) + 0j, m: -1.0 / (2 * k) + 0j})
    return m11, m12


def _m_num(e, plane, t):
    t = np.asarray(t, dtype=float)
    if e.kind == 'drift':
        return np.ones_like(t), t * float(getattr(e, 'fl', 1.0))
    kp = e.k1 if plane == 0 else -e.k1
    k = e.kappa
    if kp > 0:
        return np.cos(k * t), np.sin(k * t) / k
    return np.cosh(k * t), np.sinh(k * t) / k


def _funcs_ep(e):
    """素子内の w_x, w̄_x, w_y, w̄_y, η を _EP で (s は素子の入口からの距離)"""
    bx, ax, mx, by, ay, my, ex, epx = e.opt0
    rate, fx, fy = _rates(e)
    w0x = math.sqrt(bx) * np.exp(1j * mx)
    w0px = (1j - ax) / math.sqrt(bx) * np.exp(1j * mx)
    w0y = math.sqrt(by) * np.exp(1j * my)
    w0py = (1j - ay) / math.sqrt(by) * np.exp(1j * my)
    x11, x12 = _m_ep(e, 0)
    y11, y12 = _m_ep(e, 1)
    wx = x11 * w0x + x12 * w0px
    wy = y11 * w0y + y12 * w0py
    eta = x11 * ex + x12 * epx
    return dict(wx=wx, wxb=wx.conj(fx, fy), wy=wy, wyb=wy.conj(fx, fy), eta=eta.conj(fx, fy) if False else eta), rate, fx, fy


def _funcs_num(e, t):
    bx, ax, mx, by, ay, my, ex, epx = e.opt0
    w0x = math.sqrt(bx) * np.exp(1j * mx)
    w0px = (1j - ax) / math.sqrt(bx) * np.exp(1j * mx)
    w0y = math.sqrt(by) * np.exp(1j * my)
    w0py = (1j - ay) / math.sqrt(by) * np.exp(1j * my)
    x11, x12 = _m_num(e, 0, t)
    y11, y12 = _m_num(e, 1, t)
    wx = x11 * w0x + x12 * w0px
    wy = y11 * w0y + y12 * w0py
    return dict(wx=wx, wxb=np.conj(wx), wy=wy, wyb=np.conj(wy), eta=x11 * ex + x12 * epx)


def _thin_funcs(e):
    bx, ax, mx, by, ay, my, ex, epx = e.opt0
    wx = math.sqrt(bx) * np.exp(1j * mx)
    wy = math.sqrt(by) * np.exp(1j * my)
    return dict(wx=wx, wxb=np.conj(wx), wy=wy, wyb=np.conj(wy), eta=ex)


def _mono_thin(f, q):
    j, k, l, m, r = q
    return f['wx'] ** j * f['wxb'] ** k * f['wy'] ** l * f['wyb'] ** m * f['eta'] ** r


def _mono_ep(F, q):
    j, k, l, m, r = q
    res = _EP.const(1.0)
    for key, n in (('wx', j), ('wxb', k), ('wy', l), ('wyb', m), ('eta', r)):
        for _ in range(n):
            res = res * F[key]
    return res


def _mono_num(f, q):
    j, k, l, m, r = q
    return f['wx'] ** j * f['wxb'] ** k * f['wy'] ** l * f['wyb'] ** m * f['eta'] ** r


def _nodes(e, nslice):
    t = (np.arange(nslice) + 0.5) * e.L / nslice
    return t, np.full(nslice, e.L / nslice)


# ---------------------------------------------------------------------------------------------
# H の単項式の係数
# ---------------------------------------------------------------------------------------------
_RE_I = (1.0, 0.0, -1.0, 0.0)
_IM_I = (0.0, 1.0, 0.0, -1.0)


def poly_coeffs(n, c):
    """H_n = (1/(n+1)!) Re[c (x+iy)^(n+1)] の x^a y^b の係数 (a+b = n+1)。c = K_n + i S_n (積分値、または密度)"""
    N = n + 1
    P = {}
    for b in range(N + 1):
        a = N - b
        P[(a, b)] = comb(N, a) / math.factorial(N) * (c.real * _RE_I[b % 4] - c.imag * _IM_I[b % 4])
    return P


def mono_coef(strengths, q):
    """素子の強さ (積分値の辞書) から、単項式 q = (j,k,l,m,p) の H の係数 (∫ds の前の因子、(2J) のべきを付ける前の形)"""
    j, k, l, m, p = q
    a = j + k + p
    b = l + m
    N = a + b
    n = N - 1
    if n not in strengths:
        return 0.0
    P = poly_coeffs(n, strengths[n])
    return P[(a, b)] * comb(a, p) * comb(j + k, j) * comb(l + m, l) / 2.0 ** (j + k + l + m)


def elem_integrals(e, monos, method='analytic', nslice=8):
    """素子 e について、各単項式 q の ∫ w^j w̄^k w_y^l w̄_y^m η^r ds (薄レンズは強さを掛けた「値」の形ではなく、
    ここでは「厚さの積分」の形: 厚い素子は ∫ds、薄い素子は 1 × 値)。強さ (密度または積分値) は呼び出し側で掛ける"""
    if e.L == 0.0:
        f = _thin_funcs(e)
        return np.array([_mono_thin(f, q) for q in monos], dtype=complex)
    if method == 'analytic':
        F, rate, fx, fy = _funcs_ep(e)
        return np.array([_mono_ep(F, q).integral(e.L, rate) for q in monos], dtype=complex)
    t, wt = _nodes(e, nslice)
    f = _funcs_num(e, t)
    return np.array([np.sum(wt * _mono_num(f, q)) for q in monos], dtype=complex)


def elem_pair_integrals(e, monos, method='analytic', nslice=8):
    """素子内の二重積分 Igt[q, q'] = ∬_{s1>s2} mono_q(s1) conj(mono_q'(s2)) ds1 ds2 (薄レンズは ½ a_q conj(a_q'))"""
    n = len(monos)
    out = np.zeros((n, n), dtype=complex)
    if e.L == 0.0:
        f = _thin_funcs(e)
        v = np.array([_mono_thin(f, q) for q in monos])
        return 0.5 * np.outer(v, np.conj(v))
    if method == 'analytic':
        F, rate, fx, fy = _funcs_ep(e)
        eps = [_mono_ep(F, q) for q in monos]
        anti = [ep.conj(fx, fy).anti(rate) for ep in eps]
        for a in range(n):
            for b in range(n):
                out[a, b] = (eps[a] * anti[b]).integral(e.L, rate)
        return out
    t, wt = _nodes(e, nslice)
    f = _funcs_num(e, t)
    v = np.array([wt * _mono_num(f, q) for q in monos])         # (n, nslice) 重み込み
    for a in range(n):
        for b in range(n):
            M = np.outer(v[a], np.conj(v[b]))
            out[a, b] = np.sum(np.tril(M, -1)) + 0.5 * np.sum(np.diag(M))
    return out


def _strength_density(e, n):
    c = e.strengths.get(n, 0j)
    return c / e.L if e.L > 0 else c       # 厚い素子: 密度 (積分は ∫ds)、薄い: 積分値そのもの


def _dens_strengths(e):
    return {n: (c / e.L if e.L > 0 else c) for n, c in e.strengths.items()}


# ---------------------------------------------------------------------------------------------
# 作動点の確認
# ---------------------------------------------------------------------------------------------
def working_point_check(qx, qy, tol=0.01, order=3):
    """整数・半整数共鳴と 3 次までの共鳴線 (m Qx + n Qy = p、|m|+|n| ≤ order) への距離。戻り: dict"""
    qx, qy = float(qx), float(qy)
    d = lambda v: abs(v - round(v))
    info = {'dist_x_int': d(qx), 'dist_y_int': d(qy), 'dist_x_half': d(2 * qx) / 2, 'dist_y_half': d(2 * qy) / 2}
    lines = []
    for m in range(-order, order + 1):
        for n in range(-order, order + 1):
            if (m, n) == (0, 0) or abs(m) + abs(n) > order:
                continue
            if m < 0 or (m == 0 and n < 0):
                continue
            v = m * qx + n * qy
            lines.append((m, n, d(v)))
    lines.sort(key=lambda t: t[2])
    info['nearest_lines'] = lines[:5]
    info['near'] = [t for t in lines if t[2] < tol] + ([('x int', 0, info['dist_x_int'])] if info['dist_x_int'] < tol else []) \
        + ([('y int', 0, info['dist_y_int'])] if info['dist_y_int'] < tol else []) \
        + ([('x half', 0, info['dist_x_half'])] if info['dist_x_half'] < tol else []) \
        + ([('y half', 0, info['dist_y_half'])] if info['dist_y_half'] < tol else [])
    return info


def _check_uncoupled(rows, warnings):
    r = max(max(abs(float(r_.r1)), abs(float(r_.r2)), abs(float(r_.r3)), abs(float(r_.r4))) for r_ in rows)
    if r > 1e-6:
        warnings.append('線形の結合がある (|R| 最大 %.2e): 結合なしの式なので結果は近似 (横方向の結合は未対応)' % r)
    ey = max(abs(float(r_.ety)) for r_ in rows)
    if ey > 1e-6:
        warnings.append('垂直分散がある (最大 %.2e m): 垂直分散の項は入れていない' % ey)


def _tables(lat_or_elements, momentum=None, delta=0.0, rows=None):
    if hasattr(lat_or_elements, 'twiss_table'):
        lat = lat_or_elements
        elements = lat.elements
        if rows is None:
            rows = lat.twiss_table(delta=delta)
    else:
        elements = lat_or_elements
        if rows is None:
            raise ValueError('resdrive: 素子のリストを渡すときは rows (Twiss の行) も渡す')
    return elements, rows


# ---------------------------------------------------------------------------------------------
# 共鳴駆動項
# ---------------------------------------------------------------------------------------------
class RDTResult:
    def __init__(self):
        self.warnings = []

    def table(self):
        lines = ['共鳴駆動項 h_jklmp  (H の係数、(2J)^(1/2 次数) で規格化; 条件: %s、δ=%g、基準位置 %s)' % (self.method, self.delta, self.ref_name)]
        lines.append('  作動点 νx=%.6f νy=%.6f; 整数・半整数共鳴までの距離: x %.4f/%.4f  y %.4f/%.4f' % (
            self.qx, self.qy, self.wp['dist_x_int'], self.wp['dist_x_half'], self.wp['dist_y_int'], self.wp['dist_y_half']))
        for nm in self.names:
            v = self.total[nm]
            lines.append('  %-8s % .6e %+.6e i   |h| = %.6e' % (nm, v.real, v.imag, abs(v)))
        for w in self.warnings:
            lines.append('  (注意) ' + w)
        return '\n'.join(lines)

    __str__ = table

    def element_table(self, name, top=10):
        k = self.names.index(name)
        v = self.per_element[:, k]
        order = np.argsort(-np.abs(v))[:top]
        return [(self.elem_names[i], self.elem_s[i], v[i]) for i in order]


# 偏向磁石の 3 次の項を共鳴駆動項 h に入れるか (v2.3.0 R4)。検証 (docs/49 §9) の結果で既定を決めた
BEND_RDT_DEFAULT = True


def bend_rdt(elements, rows, monos, nslice=16, momentum=3e9, mu_ref=(0.0, 0.0), warnings=None, delta=0.0):
    """偏向磁石の 3 次の項 (曲率 (h/2)x(px²+py²)・硬い端 ±(h/2)px y²・複合機能の本体・E1/E2≠0 の端。`bendnl.py` が K64 の追跡の写像から取り出す) の
    共鳴駆動項 h_(j,k,l,m,p) への寄与 (p ≥ 1 は分散による軌道のずれ x → x + ηδ の展開)。j+k+l+m+p = 3 の項だけが 0 でない。
    戻り: list[(素子の添字, 名前, s 入口, ndarray(len(monos)) 複素)] (対象の偏向磁石だけ)"""
    from . import bendnl
    warnings = warnings if warnings is not None else []
    skipped = {}
    out = []
    for i, el in enumerate(elements):
        if _tn(el) != 'Bend':
            if bendnl.needs_other(el, rows[i], delta):
                # v2.3.0 R6: δ≠0・閉軌道がずれたときの四極・ドリフトの 3 次の項 (運動学的な項・フリンジの軌道のずれによる展開)。p=0 だけ
                p0 = [q for q in monos if q[4] == 0]
                try:
                    a, info = bendnl.other_atoms(el, rows[i], rows[i + 1], [q[:4] for q in p0], delta=delta, momentum=momentum)
                except Exception as ex:
                    warnings.append('%s: δ・軌道のずれの 3 次の項の取り出しに失敗 (%s: %s)' % (el.name, type(ex).__name__, ex))
                    continue
                tot = a.sum(axis=1) / 2.0 ** 1.5
                vec = np.zeros(len(monos), dtype=complex)
                ip = 0
                for iq, q in enumerate(monos):
                    if q[4] == 0:
                        vec[iq] = tot[ip] * np.exp(-1j * ((q[0] - q[1]) * mu_ref[0] + (q[2] - q[3]) * mu_ref[1]))
                        ip += 1
                out.append((i, el.name, float(rows[i].s), vec))
            continue
        L = float(getattr(el, 'L', 0.0) or 0.0)
        ang = float(getattr(el, 'ANGLE', 0.0) or 0.0)
        k0 = float(getattr(el, 'K0', 0.0) or 0.0)
        if L <= 0 or (ang == 0.0 and k0 == 0.0):
            continue
        if float(getattr(el, 'ROTATE', 0.0) or 0.0) or float(getattr(el, 'ROTATION', 0.0) or 0.0):
            skipped.setdefault('ROTATE・ROTATION≠0', []).append(el.name)
            continue
        if any(float(getattr(el, k, 0.0) or 0.0) for k in ('F1', 'F2', 'FRINGE', 'DISFRIN')):
            skipped.setdefault('F1・F2・FRINGE・DISFRIN の指定 (入れるが線形のフリンジ・端の有無は反映しない)', []).append(el.name)
        try:
            if delta != 0.0:
                # v2.3.0 R6: δ≠0 は、δ の閉軌道・光学 (rows は δ の Twiss) のまわりで取り出した原子の p=0 の項だけ (p≥1 の展開は δ=0 の分散のまわりの量なので 0)
                p0 = [q for q in monos if q[4] == 0]
                ph = np.array([np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[1]), np.exp(-1j * mu_ref[1])])
                a, info = bendnl.bend_atoms_auto(el, rows[i], rows[i + 1], [q[:4] for q in p0], nslice=nslice, momentum=momentum, delta=delta, check=True)
                tot = a.sum(axis=1) / 2.0 ** 1.5
                vec = np.zeros(len(monos), dtype=complex)
                ip = 0
                for iq, q in enumerate(monos):
                    if q[4] == 0:
                        vec[iq] = tot[ip] * np.exp(-1j * ((q[0] - q[1]) * mu_ref[0] + (q[2] - q[3]) * mu_ref[1]))
                        ip += 1
                out.append((i, el.name, float(rows[i].s), vec))
                continue
            if k0 != 0.0:
                # v2.3.0 R6: K0≠0 (挿入光源) は、p=0 は閉軌道のまわりの原子 (bend_atoms_orbit)、p=1・2 は要素の追跡から直接 (bendnl.element_delta_terms。自然な項を含む全体)
                p0 = [q for q in monos if q[4] == 0]
                ph = np.array([np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[1]), np.exp(-1j * mu_ref[1])])
                d_ = bendnl.orbit_decomposition(el, bendnl.row_orbit(rows[i]), 0.0, nslice, momentum)
                w = bendnl.w_from_row(rows[i]) * ph
                am = []
                for kind, Rpre, f, Rpost in d_['units']:
                    wc = Rpre @ w
                    am.append(bendnl.harmonic_amplitudes(f, wc, [q[:4] for q in p0]).sum(axis=1))
                    w = Rpost @ wc
                tot = np.sum(am, axis=0)
                dt = bendnl.element_delta_terms(el, rows[i], rows[i + 1], momentum, mu_ref=mu_ref)
                vec = np.zeros(len(monos), dtype=complex)
                ip = 0
                for iq, q in enumerate(monos):
                    if q[4] == 0:
                        vec[iq] = tot[ip]
                        ip += 1
                    else:
                        vec[iq] = dt.get(tuple(q), 0.0)
                out.append((i, el.name, float(rows[i].s), vec))
                continue
            a, info = bendnl.bend_atoms_rdt(el, rows[i], rows[i + 1], monos, nslice=nslice, momentum=momentum, mu_ref=mu_ref)
        except Exception as ex:
            warnings.append('偏向磁石 %s: 3 次の項の取り出しに失敗 (%s: %s)。この磁石の項は入れていない' % (el.name, type(ex).__name__, ex))
            continue
        if info['w_err'] > 1e-6 or info['e_err'] > 1e-4:
            warnings.append('偏向磁石 %s: スライスの線形写像で伝播した w・分散が Twiss と違う (w %.1e、分散 %.1e): 項は不正確' % (el.name, info['w_err'], info['e_err']))
        out.append((i, el.name, float(rows[i].s), a.sum(axis=1)))
    for why, names in skipped.items():
        warnings.append('偏向磁石の 3 次の項: %s の偏向磁石 %d 台 (%s ...) は %s' % (why, len(names), names[0], '入れていない' if why[:2] in ('K0', 'RO') else '注意'))
    return out


def _resonance_terms_cubic(lat, terms=DEFAULT_TERMS, method='analytic', nslice=8, delta=0.0, ref=None, rows=None,
                           bend_cubic=None, bend_nslice=None, momentum_eV=3e9):
    """6 極 (と 8 極) の共鳴駆動項 h_jklmp (名前 'h21000' など、(j,k,l,m,p) のタプルも可)。
    lat: Lattice (または素子のリストと rows)。ref: 位相の基準 (素子の添字、既定 None = リングの先頭)。
    v2.3.0 R4: bend_cubic=True (既定) で、偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2。bendnl.py が追跡の写像から取り出す。p ≥ 1 は分散による
    軌道のずれ) の寄与を足す (.total = .total_sext + .total_bend)。False で 6 極 (と 8 極) だけ (R3 までの値)。δ≠0 では入れない (警告)。
    戻り: RDTResult (.total[名前] 複素数、.per_element[素子, 項] (6 極と偏向磁石。s の順)、.elem_names・.elem_idx・.elem_s、.warnings、
    .total_sext・.total_bend)"""
    elements, rows = _tables(lat, delta=delta, rows=rows)
    res = RDTResult()
    _check_uncoupled(rows, res.warnings)
    qx, qy = float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI
    mu_ref = (0.0, 0.0)
    ref_name = 'リング先頭'
    if ref is not None:
        mu_ref = (float(rows[ref].mux), float(rows[ref].muy))
        ref_name = '%s の入口' % elements[ref].name
    nl = collect_nonlinear(elements, rows, delta, mu_ref, res.warnings)
    tl = [parse_term(t) for t in terms]
    monos = [(j, k, l, m, p) for (j, k, l, m, p) in tl]
    pe = np.zeros((len(nl), len(tl)), dtype=complex)
    for ie, e in enumerate(nl):
        integ = elem_integrals(e, monos, method, nslice)
        dens = _dens_strengths(e)
        for it, q in enumerate(monos):
            c = mono_coef(dens if e.L > 0 else e.strengths, q)
            pe[ie, it] = c * integ[it]
    res.names = [_term_names(t) for t in tl]
    res.terms = tl
    res.total_sext = {nm: complex(pe[:, i].sum()) for i, nm in enumerate(res.names)}
    res.total_bend = {nm: 0j for nm in res.names}
    res.bend_included = False
    names_e = [e.name for e in nl]
    idx_e = [e.idx for e in nl]
    s_e = [e.s0 for e in nl]
    if bend_cubic is None:
        bend_cubic = BEND_RDT_DEFAULT
    if bend_cubic:
        if True:
            bl = bend_rdt(elements, rows, monos, nslice=bend_nslice or BEND_NSLICE_DEFAULT, momentum=momentum_eV, mu_ref=mu_ref, warnings=res.warnings, delta=delta)
            if bl:
                rows_b = np.array([b[3] for b in bl])
                pe = np.vstack([pe, rows_b])
                names_e += [b[1] for b in bl]
                idx_e += [b[0] for b in bl]
                s_e += [b[2] for b in bl]
                order = np.argsort(idx_e, kind='stable')
                pe = pe[order]
                names_e = [names_e[k] for k in order]
                idx_e = [idx_e[k] for k in order]
                s_e = [s_e[k] for k in order]
                res.total_bend = {nm: complex(rows_b[:, i].sum()) for i, nm in enumerate(res.names)}
            res.bend_included = True
    res.per_element = pe
    res.total = {nm: complex(pe[:, i].sum()) for i, nm in enumerate(res.names)}
    res.elem_names = names_e
    res.elem_idx = idx_e
    res.elem_s = s_e
    res.method, res.delta, res.ref_name = method + ('(nslice=%d)' % nslice if method == 'numeric' else ''), delta, ref_name
    res.qx, res.qy = qx, qy
    res.wp = working_point_check(qx, qy)
    res.momentum_note = '運動量は光学の計算に使うだけ (強さは SAD の K_nL = 運動量に対する規格化済み)'
    return res


# v2.3.0 R5: 4 次の項 (j+k+l+m = 4、p = 0) と 6 極の 2 次 (rdt2nd.py)。既定は order2=True
ORDER2_DEFAULT = True


def resonance_terms(lat, terms=DEFAULT_TERMS, method='analytic', nslice=8, delta=0.0, ref=None, rows=None,
                    bend_cubic=None, bend_nslice=None, momentum_eV=3e9, order2=None, rdt4_options=None):
    """共鳴駆動項 h_jklmp。3 次 (j+k+l+m = 3) と δ の項 (p ≥ 1) は 6 極の 1 次 (+ 偏向磁石の 3 次の項、R4)。
    **v2.3.0 R5: j+k+l+m = 4 (p=0) の名前 ('h22000'、'h31000'…) は 4 次の項 (写像の 4 次の生成関数の係数: 8 極 1 次 + 運動学 + 端のフリンジ + 偏向磁石の 4 次の項 + 3 次 (6 極・偏向磁石) の 2 次、
    `rdt2nd.rdt4`、docs/49 §10)**。order2=True (既定) で 3 次の 2 次 (Magnus の第 2 項) を入れる、False で 4 次のハミルトニアンの 1 次だけ。
    4 次の項は ref=None・δ=0 のリング先頭が位相の基準。rdt4_options: rdt4 への追加のオプション (dict)。戻り: RDTResult
    (4 次の項は .total に入る。.per_element は 3 次・δ の項だけ。.rdt4 に RDT4Result (内訳 .parts))"""
    tl = [parse_term(t) for t in terms]
    four = [t for t in tl if sum(t[:4]) == 4 and t[4] == 0]
    if not four:
        return _resonance_terms_cubic(lat, terms, method, nslice, delta, ref, rows, bend_cubic, bend_nslice, momentum_eV)
    rest = [t for t in tl if t not in four]
    if delta != 0.0 or ref is not None:
        raise ValueError('resonance_terms: 4 次の項は δ=0、ref=None (リング先頭が基準) だけ')
    from . import rdt2nd
    if rest:
        res = _resonance_terms_cubic(lat, rest, method, nslice, delta, ref, rows, bend_cubic, bend_nslice, momentum_eV)
    else:
        elements, rows_ = _tables(lat, delta=delta, rows=rows)
        res = RDTResult()
        res.names, res.terms = [], []
        res.total, res.total_sext, res.total_bend = {}, {}, {}
        res.per_element = np.zeros((0, 0), dtype=complex)
        res.elem_names, res.elem_idx, res.elem_s = [], [], []
        res.method, res.delta, res.ref_name = method, delta, 'リング先頭'
        res.qx, res.qy = float(rows_[-1].mux) / TWO_PI, float(rows_[-1].muy) / TWO_PI
        res.wp = working_point_check(res.qx, res.qy)
    opts = dict(rdt4_options or {})
    o2 = ORDER2_DEFAULT if order2 is None else order2
    r4 = rdt2nd.rdt4(lat, order2=o2, method=method, nslice=nslice, bend_nslice=bend_nslice or BEND_NSLICE_DEFAULT, rows=rows, momentum_eV=momentum_eV, **opts)
    res.rdt4 = r4
    for t in four:
        nm = _term_names(t)
        res.names.append(nm)
        res.terms.append(t)
        v = r4.total['h%d%d%d%d0' % t[:4]]
        res.total[nm] = v
        res.total_sext[nm] = 0j
        res.total_bend[nm] = 0j
    res.warnings += [w for w in r4.warnings if w not in res.warnings]
    return res


# ---------------------------------------------------------------------------------------------
# ADTS
# ---------------------------------------------------------------------------------------------
def _monos_of_degree(N):
    out = []
    for j in range(N + 1):
        for k in range(N + 1 - j):
            for l in range(N + 1 - j - k):
                m = N - j - k - l
                out.append((j, k, l, m, 0))
    return out


class ADTSResult:
    def __init__(self):
        self.warnings = []

    def table(self):
        lines = ['振幅依存チューンシフト  J = アクション (ε 相当、m·rad)。条件: %s、δ=%g、6 極は 2 次・8 極は 1 次、δ の項と結合なし' % (self.method, self.delta)]
        lines.append('  作動点 νx=%.6f νy=%.6f; 整数・半整数共鳴までの距離: x %.4f/%.4f  y %.4f/%.4f' % (
            self.qx, self.qy, self.wp['dist_x_int'], self.wp['dist_x_half'], self.wp['dist_y_int'], self.wp['dist_y_half']))
        lines.append('                       6 極 (2 次)      8 極 (1 次)       合計   [1/(m·rad)]')
        for nm, a, b in (('dνx/dJx', self.sext['xx'], self.octu['xx']), ('dνx/dJy', self.sext['xy'], self.octu['xy']), ('dνy/dJy', self.sext['yy'], self.octu['yy'])):
            lines.append('  %-8s % .6e  % .6e  % .6e' % (nm, a, b, a + b))
        lines.append('  (6 極の虚部の残り %.1e)' % self.imag_residual)
        lines.append('  運動学的な項 (xx, xy, yy) = %.6e %.6e %.6e; 端部のフリンジ = %.6e %.6e %.6e; 全部入りの合計 = %.6e %.6e %.6e' % (
            self.kin['xx'], self.kin['xy'], self.kin['yy'], self.fringe['xx'], self.fringe['xy'], self.fringe['yy'],
            self.total['xx'], self.total['xy'], self.total['yy']))
        if getattr(self, 'bend_included', False):
            lines.append('  偏向磁石の 3 次の項 (R3): 3 次どうしの対 (xx, xy, yy) = %.6e %.6e %.6e; 6 極との交差項 = %.6e %.6e %.6e  (原子 %d 個)' % (
                self.bend['xx'], self.bend['xy'], self.bend['yy'], self.cross['xx'], self.cross['xy'], self.cross['yy'], self.bend_natoms))
        for w in self.warnings:
            lines.append('  (注意) ' + w)
        return '\n'.join(lines)

    __str__ = table

    @property
    def dnux_dJx(self): return self.sext['xx'] + self.octu['xx']

    @property
    def dnux_dJy(self): return self.sext['xy'] + self.octu['xy']

    @property
    def dnuy_dJy(self): return self.sext['yy'] + self.octu['yy']

    @property
    def total(self):
        """6 極 2 次 + 8 極 1 次 + 運動学的な項 + 端部のフリンジ + (v2.3.0 R3、bend_cubic=True の既定) 偏向磁石の 3 次どうし + 6 極との交差項
        (追跡の ADTS に対応する全部入りの値): {'xx','xy','yy'}"""
        b4 = getattr(self, 'bend4', None) or {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
        return {k: self.sext[k] + self.octu[k] + self.kin[k] + self.fringe[k] + b4[k] + self.bend[k] + self.cross[k] for k in ('xx', 'xy', 'yy')}

    total_xx = property(lambda self: self.total['xx'])
    total_xy = property(lambda self: self.total['xy'])
    total_yy = property(lambda self: self.total['yy'])


def _deriv_from_K(Kc):
    """K = Σ c[(ax,ay)] Jx^ax Jy^ay (2 次) → (∂²K/∂Jx², ∂²K/∂Jx∂Jy, ∂²K/∂Jy²)"""
    return 2.0 * Kc.get((2, 0), 0.0), Kc.get((1, 1), 0.0), 2.0 * Kc.get((0, 2), 0.0)


def _second_order(byk, a_all, igt, Ne, qx, qy, warnings):
    """3 次の項の 2 次の摂動 (対の項) の係数 Kc (キー (ax, ay) = Jx^ax Jy^ay、2 次)、要素ごとの寄与 (3, Ne)、調和ごとの寄与。
    a_all[q][e] = e 番目の要素の ∫F_q (強さ・係数込み)、igt[(q, q')][e] = 素子内の対の積分 (薄い要素は ½ a_q conj(a_q'))。要素は s の順"""
    Kc = {}
    per_e3 = np.zeros((3, Ne), dtype=complex)
    harm = {}
    for kv, qs in (byk.items() if Ne > 0 else ()):
        xk = kv[0] * qx + kv[1] * qy
        sn = math.sin(math.pi * xk)
        if abs(sn) < 1e-9:
            warnings.append('調和 k=(%d,%d): k·Q = %.6f が整数 (共鳴) で 2 次の式が発散' % (kv[0], kv[1], xk))
            continue
        em, ep = np.exp(-1j * math.pi * xk), np.exp(1j * math.pi * xk)
        pref = -0.5 / (4.0 * math.pi * sn)
        for q in qs:
            a = a_all[q]
            ca = np.cumsum(a)
            for q2 in qs:
                b = np.conj(a_all[q2])
                cb = np.cumsum(b)
                Tb, Ta = cb[-1], ca[-1]
                cb_before = cb - b
                ca_before = ca - a
                ig = igt[(q, q2)]
                ab = a * b
                intra = em * ig + ep * (ab - ig)
                # 行の和 (s1 側がこの要素) と列の和 (s2 側がこの要素): 要素ごとの寄与は両者の平均
                row = a * (em * cb_before + ep * (Tb - cb)) + intra
                col = b * (em * (Ta - ca) + ep * ca_before) + intra
                contrib = 0.5 * (row + col)
                px, py = (q[0] + q[1] + q2[0] + q2[1]) / 2.0, (q[2] + q[3] + q2[2] + q2[3]) / 2.0
                for (ax, ay), co in (((px - 1, py), kv[0] * px), ((px, py - 1), kv[1] * py)):
                    if co == 0:
                        continue
                    key = (int(round(ax)), int(round(ay)))
                    Kc[key] = Kc.get(key, 0j) + pref * co * contrib.sum()
                    idx = {(2, 0): 0, (1, 1): 1, (0, 2): 2}[key]
                    d = {(2, 0): 2.0, (1, 1): 1.0, (0, 2): 2.0}[key]
                    per_e3[idx] += pref * co * contrib * d
                    h = harm.setdefault(kv, np.zeros(3, dtype=complex))
                    h[idx] += pref * co * contrib.sum() * d
    return Kc, per_e3, harm


# 偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2≠0 の端) を入れるか。検証 (docs/49 §8) の結果で既定を決めた
BEND_CUBIC_DEFAULT = True
BEND_NSLICE_DEFAULT = 16
# 4 次の項 (運動学・フリンジ・偏向磁石の 4 次) の既定は 'r3' (物理的な項だけ: 四極は厳密な積分、偏向磁石は運動学的な項だけ)。'r6' = 精密版のオプション (四極の K64 積分器の台形則 + 偏向磁石の 4 次の項。0.1〜0.3% の端数を含む)
FOURTH_DEFAULT = 'r3'


def bend_cubic_adts(elements, rows, nl3, a_sext, igt_sext, qx, qy, nslice=BEND_NSLICE_DEFAULT, momentum=3e9, warnings=None, delta=0.0, orbit_policy='auto'):
    """偏向磁石の 3 次の項 (R3): K64 の追跡の写像から取り出した偏向磁石のリー生成関数 (`bendnl.py`) を、6 極 (nl3 の厚い・薄い要素。解析式) と s の順に並べて、
    2 次の摂動 (対の項) にかける。戻り dict: 'bend' (偏向磁石の 3 次の項どうしの対)、'cross' (6 極 × 偏向磁石)、'sext' (6 極だけ。adts の .sext と一致するはず)、
    'all' (全部)。各 {'xx','xy','yy'} [1/(m rad)]。nslice: 偏向磁石の本体の分割数 (誤差 ~ 1/nslice²)"""
    from . import bendnl
    warnings = warnings if warnings is not None else []
    monos = _monos_of_degree(3)
    byk = {}
    for q in monos:
        byk.setdefault((q[0] - q[1], q[2] - q[3]), []).append(q)
    # 単位の列 (s の順): 6 極 (1 要素 = 1 単位)、偏向磁石 (nslice + 2 個の薄い kick)
    sext_at = {e.idx: ie for ie, e in enumerate(nl3)}
    units = []                       # (種類, 添字 or 振幅 (20, n))
    skipped = {}
    for i, el in enumerate(elements):
        tn = _tn(el)
        if i in sext_at:
            units.append(('s', sext_at[i]))
        if tn != 'Bend':
            if bendnl.needs_other(el, rows[i], delta):
                try:
                    a, info = bendnl.other_atoms(el, rows[i], rows[i + 1], monos, delta=delta, momentum=momentum)
                    units.append(('b', a))
                except Exception as ex:
                    warnings.append('%s (%s): δ・軌道のずれの 3 次の項の取り出しに失敗 (%s: %s)' % (el.name, tn, type(ex).__name__, ex))
            continue
        L = float(getattr(el, 'L', 0.0) or 0.0)
        ang = float(getattr(el, 'ANGLE', 0.0) or 0.0)
        k0 = float(getattr(el, 'K0', 0.0) or 0.0)
        if L <= 0 or (ang == 0.0 and k0 == 0.0):
            continue
        if float(getattr(el, 'ROTATE', 0.0) or 0.0) or float(getattr(el, 'ROTATION', 0.0) or 0.0):
            skipped.setdefault('ROTATE・ROTATION≠0', []).append(el.name)
            continue
        if float(getattr(el, 'DX', 0.0) or 0.0) or float(getattr(el, 'DY', 0.0) or 0.0):
            skipped.setdefault('DX・DY≠0 (入れるが据え付けのずれは無視)', []).append(el.name)
        if any(float(getattr(el, k, 0.0) or 0.0) for k in ('F1', 'F2', 'FRINGE', 'DISFRIN')):
            skipped.setdefault('F1・F2・FRINGE・DISFRIN の指定 (入れるが線形のフリンジ・端の有無は反映しない)', []).append(el.name)
        try:
            # v2.3.0 R6: K0≠0 (挿入光源)・δ≠0・閉軌道がずれた偏向磁石は、閉軌道・δ のまわりで取り出す (bend_atoms_orbit。それ以外は同じ値)
            a, info = bendnl.bend_atoms_auto(el, rows[i], rows[i + 1], monos, nslice=nslice, momentum=momentum, delta=delta, policy=orbit_policy)
        except Exception as ex:          # 追跡からの取り出しに失敗した偏向磁石は入れない (従来どおり)
            warnings.append('偏向磁石 %s: 3 次の項の取り出しに失敗 (%s: %s)。この磁石の 3 次の項は入れていない' % (el.name, type(ex).__name__, ex))
            continue
        if info['w_err'] > 1e-6:
            warnings.append('偏向磁石 %s: スライスの線形写像で伝播した w が Twiss と %.1e 違う (3 次の項は不正確)' % (el.name, info['w_err']))
        units.append(('b', a))
    for why, names in skipped.items():
        warnings.append('偏向磁石の 3 次の項: %s の偏向磁石 %d 台 (%s ...) は %s' % (why, len(names), names[0], '入れていない' if why[:2] in ('K0', 'RO') else '注意'))
    # 並べて配列にする
    n_tot = sum(1 if k == 's' else v.shape[1] for k, v in units)
    A = {q: np.zeros(n_tot, dtype=complex) for q in monos}
    IG = {}
    kind = np.zeros(n_tot, dtype=int)          # 0: 6 極、1: 偏向磁石
    pos = 0
    pairs = [(q, q2) for qs in byk.values() for q in qs for q2 in qs]
    for pr_ in pairs:
        IG[pr_] = np.zeros(n_tot, dtype=complex)
    for k, v in units:
        if k == 's':
            for iq, q in enumerate(monos):
                A[q][pos] = a_sext[q][v]
            for pr_ in pairs:
                IG[pr_][pos] = igt_sext[pr_][v]
            kind[pos] = 0
            pos += 1
        else:
            n = v.shape[1]
            for iq, q in enumerate(monos):
                A[q][pos:pos + n] = v[iq]
            for (q, q2) in pairs:
                IG[(q, q2)][pos:pos + n] = 0.5 * v[monos.index(q)] * np.conj(v[monos.index(q2)])
            kind[pos:pos + n] = 1
            pos += n

    def run(mask):
        idx = np.nonzero(mask)[0]
        if len(idx) == 0:
            return {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
        a_ = {q: A[q][idx] for q in monos}
        ig_ = {p: IG[p][idx] for p in pairs}
        Kc, _, _ = _second_order(byk, a_, ig_, len(idx), qx, qy, [])
        xx, xy, yy = _deriv_from_K({k_: v_.real for k_, v_ in Kc.items()})
        return {'xx': xx, 'xy': xy, 'yy': yy}
    r_all = run(np.ones(n_tot, dtype=bool))
    r_s = run(kind == 0)
    r_b = run(kind == 1)
    cross = {k: r_all[k] - r_s[k] - r_b[k] for k in ('xx', 'xy', 'yy')}
    return {'all': r_all, 'sext': r_s, 'bend': r_b, 'cross': cross, 'natoms': int((kind == 1).sum())}


def adts(lat, method='analytic', nslice=8, delta=0.0, include_octupole=True, rows=None, detail=False, include_kinematic=True,
         bend_cubic=None, bend_nslice=None, momentum_eV=3e9, fourth=None, optics=None):
    """6 極の 2 次 + 8 極の 1 次の ADTS ∂νx/∂Jx、∂νx/∂Jy (= ∂νy/∂Jx)、∂νy/∂Jy。J はアクション (x = √(2Jβ) cos φ)。
    delta≠0 ならそのエネルギーでの Twiss・閉軌道の光学を使う (6 極の強さは 1/(1+δ) 倍、閉軌道のずれによる 8 極→6 極のフィードダウンは入れない)。
    戻り: ADTSResult (.sext・.octu の辞書 'xx','xy','yy'、.dnux_dJx など (6 極 + 8 極)、.per_element・.by_harmonic)。
    v2.3.0 R2: 追跡の全部入りの値に近づけるため、磁石の多極以外の 1 次の項を別に出す: .kin (運動学的な項 p⁴/8、kinematic_adts)、
    .fringe (四極・K1 つきの偏向磁石の端部の非線形のフリンジ、K64 の既定 ON、quad_fringe_adts)。
    v2.3.0 R3: bend_cubic=True (既定) で、偏向磁石の 3 次の項 (曲率 (h/2)x(px²+py²)・硬い端 ±(h/2)px y²・複合機能の本体・E1/E2≠0 の端。`bendnl.py` が K64 の追跡の写像から取り出す)
    の 2 次の対の項 .bend と、6 極との交差項 .cross (6 極の強さに比例) を足す。bend_nslice = 偏向磁石の本体の分割数 (既定 16、誤差 ~ 1/n²)。
    v2.3.0 V2: 4 次の項の既定は fourth='r3' (物理的な項だけ)。fourth='r6' は K64 の四極の積分器の台形則と偏向磁石の 4 次の項 (0.1〜0.3% の端数) を含む精密版。
    δ≠0 の光学の既定は CALC (optics=None・'calc')。optics='tracked' は追跡の閉軌道・ヤコビ行列から作る光学。
    .total_xx などが 6 つの和 (include_kinematic=False で .kin・.fringe は 0、bend_cubic=False で .bend・.cross は 0 = R2 の値)。δ≠0 では 3 次の項は入れない (警告)"""
    elements, rows = _tables(lat, delta=delta, rows=rows) if (rows is not None or delta == 0.0 or optics in (None, 'calc')) else (lat.elements, None)
    if rows is None:
        # optics='tracked' (オプション): δ≠0 の光学を追跡の閉軌道 (2 次の分散を含む) とヤコビ行列から作る (rdttrack.tracked_rows)。既定 (None・'calc') は CALC の光学 (lat.twiss_table(delta))
        from . import rdttrack
        rows = rdttrack.tracked_rows(elements, momentum_eV, delta, lat.twiss_table(delta=delta))
    res = ADTSResult()
    _check_uncoupled(rows, res.warnings)
    if delta != 0.0:
        res.warnings.append('δ≠0 の解析式 (既定: fourth="r3"・CALC の光学): ADTS の値は追跡と 4e-3 以内 (ASTRID2・FODO、δ = ±0.002)、δ 微分は yy などで外れる。精度が要るときは nform.adts_at_delta / adts_delta_derivative(method="nform") (追跡の写像の数値の正規形)、または fourth="r6"・optics="tracked" (精密版)')
    qx, qy = float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI
    nl = collect_nonlinear(elements, rows, delta, (0.0, 0.0), res.warnings)
    nl3 = [e for e in nl if 2 in e.strengths]
    nl4 = [e for e in nl if 3 in e.strengths]
    # ---- 6 極 2 次
    monos = _monos_of_degree(3)
    byk = {}
    for q in monos:
        byk.setdefault((q[0] - q[1], q[2] - q[3]), []).append(q)
    Ne = len(nl3)
    a_all = {q: np.zeros(Ne, dtype=complex) for q in monos}          # ∫ F_q (強さ・係数込み、(2J) → J の 2^{p/2} 込み)
    igt = {}                                                          # igt[(q,q')][e]
    for ie, e in enumerate(nl3):
        dens = _dens_strengths(e)
        integ = elem_integrals(e, monos, method, nslice)
        pi = elem_pair_integrals(e, monos, method, nslice)
        cf = np.array([mono_coef(dens, q) * 2.0 ** ((q[0] + q[1] + q[2] + q[3]) / 2.0) for q in monos])
        for iq, q in enumerate(monos):
            a_all[q][ie] = cf[iq] * integ[iq]
        for iq, q in enumerate(monos):
            for jq, q2 in enumerate(monos):
                igt.setdefault((q, q2), np.zeros(Ne, dtype=complex))[ie] = cf[iq] * cf[jq] * pi[iq, jq]
    Kc, per_e3, harm = _second_order(byk, a_all, igt, Ne, qx, qy, res.warnings)
    xx, xy, yy = _deriv_from_K({k: v.real for k, v in Kc.items()})
    res.imag_residual = max([abs(v.imag) for v in Kc.values()] + [0.0])
    res.sext = {'xx': xx, 'xy': xy, 'yy': yy}
    res.per_element = per_e3.real                      # (3, Ne): xx, xy, yy
    res.elem_names = [e.name for e in nl3]
    res.elem_idx = [e.idx for e in nl3]
    res.elem_s = [e.s0 for e in nl3]
    res.by_harmonic = {kv: tuple(v.real) for kv, v in harm.items()}
    # ---- v2.3.0 R3: 偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2) の 2 次の対の項と 6 極との交差項
    if fourth is None:
        fourth = FOURTH_DEFAULT
    res.bend = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    res.cross = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    res.bend_included = False
    if bend_cubic is None:
        bend_cubic = BEND_CUBIC_DEFAULT
    if bend_cubic:
        # v2.3.0 R6: δ≠0 でも、δ の閉軌道・光学のまわりで取り出した偏向磁石の項を入れる (bendnl.bend_atoms_orbit)。R5 までは δ≠0 では入れなかった
        bc = bend_cubic_adts(elements, rows, nl3, a_all, igt, qx, qy, nslice=bend_nslice or BEND_NSLICE_DEFAULT, momentum=momentum_eV,
                             warnings=res.warnings, delta=delta, orbit_policy=('auto' if fourth == 'r6' else 'k0'))
        res.bend, res.cross, res.bend_included = bc['bend'], bc['cross'], True
        res.bend_check = dict(bc['sext'])
        res.bend_natoms = bc['natoms']
    # ---- 8 極 1 次 (k = 0 の単項式)
    ok = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    if include_octupole and nl4:
        m4 = [q for q in _monos_of_degree(4) if q[0] == q[1] and q[2] == q[3]]
        Kc4 = {}
        for e in nl4:
            dens = _dens_strengths(e)
            integ = elem_integrals(e, m4, method, nslice)
            for iq, q in enumerate(m4):
                c = mono_coef(dens, q) * 2.0 ** ((q[0] + q[1] + q[2] + q[3]) / 2.0) * integ[iq] / TWO_PI
                key = (q[0] + q[1], q[2] + q[3])
                key = (key[0] // 2, key[1] // 2)
                Kc4[key] = Kc4.get(key, 0j) + c
        a, b, c_ = _deriv_from_K({k: v.real for k, v in Kc4.items()})
        ok = {'xx': a, 'xy': b, 'yy': c_}
    res.octu = ok
    res.kin = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    res.fringe = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    res.bend4 = {'xx': 0.0, 'xy': 0.0, 'yy': 0.0}
    res.fourth = None
    if include_kinematic and fourth == 'r6':
        # 精密版のオプション (fourth='r6'): 四極は K64 の積分器の台形則 (kinematic_h4(quad_rule='k64'))、偏向磁石 (K0≠0 を含む) は要素の追跡の写像から取り出した 4 次の項の全部
        # (運動学・端・K1 の端のフリンジ・積分器の誤差。数値積分器の 0.1〜0.3% の端数を含む)。既定は物理的な項だけの 'r3' (docs/49 §8.1、§11)
        from . import rdt2nd
        bu = rdt2nd.bend_units(elements, rows, bend_nslice or BEND_NSLICE_DEFAULT, momentum_eV, res.warnings, delta=delta)
        f4 = rdt2nd.fourth_order_adts(elements, rows, bu, momentum_eV, delta=delta, warnings=res.warnings)
        res.kin, res.fringe, res.bend4 = f4['kin'], f4['fringe'], f4['bend4']
        res.fourth = 'r6'
    elif include_kinematic:
        res.kin = kinematic_adts(elements, rows, warnings=res.warnings)
        res.fringe = quad_fringe_adts(elements, rows, warnings=res.warnings)
        if delta != 0.0:           # p は p0 で規格化した正準運動量: p⁴/(8(1+δ)³)。四極端部のフリンジは δ に依らない (R6: 追跡の kick は p0 で規格化した場の強さ)
            res.kin = {k: v / (1.0 + delta) ** 3 for k, v in res.kin.items()}
        res.fourth = 'r3'
    res.method, res.delta = method + ('(nslice=%d)' % nslice if method == 'numeric' else ''), delta
    res.qx, res.qy = qx, qy
    res.wp = working_point_check(qx, qy)
    return res


# ---------------------------------------------------------------------------------------------
# 運動学的な項 (px²+py²)²/8 の 1 次の ADTS (v2.3.0 R2)
# ---------------------------------------------------------------------------------------------
def _slice_for_optics(el, l):
    """素子 el を長さ l に縮めた複製 (内部の光学を線形の転送行列で求めるため)。四極・偏向磁石 (端部なし、角度は長さに比例)・多極 (K1 は長さに比例)・ドリフト類"""
    import dataclasses
    tn = _tn(el)
    L = float(el.L)
    f = l / L
    if tn == 'Quadrupole':
        return dataclasses.replace(el, L=l)
    if tn == 'Bend':
        return dataclasses.replace(el, L=l, ANGLE=float(el.ANGLE) * f, E1=0.0, E2=0.0)
    if tn == 'Multipole':
        return dataclasses.replace(el, L=l, KnL={n: v * f for n, v in (el.KnL or {}).items()})
    return None          # ドリフトと同じ扱い (SEXT・OCT・DRIFT・MARK・CAVI ...)


def kinematic_adts(lat_or_elements, rows=None, nsub=8, warnings=None):
    """運動学的な項 H4 = (px² + py²)²/8 (厳密なドリフトの展開 −√(1−p²) の 4 次) の 1 次の ADTS。px = x' = √(2J) Re[w' e^{iφ}]、|w'|² = γ = (1+α²)/β として
      ∂νx/∂Jx = (3/16π) ∮ γx² ds、∂νx/∂Jy = (1/8π) ∮ γx γy ds、∂νy/∂Jy = (3/16π) ∮ γy² ds。
    γ(s) は素子の内部を線形の転送行列で伝播する (四極・偏向磁石は nsub 区間の Simpson 則、ドリフト・6 極・8 極は γ 一定)。
    戻り dict 'xx'、'xy'、'yy' [1/(m rad)]。偏向磁石の曲率の項 (h x p²/2) と端部の非線形は入っていない"""
    import math
    from .elements import Drift
    if hasattr(lat_or_elements, 'twiss_table'):
        elements = lat_or_elements.elements
        rows = rows if rows is not None else lat_or_elements.twiss_table()
    else:
        elements = lat_or_elements
    if nsub % 2:
        nsub += 1
    sxx = sxy = syy = 0.0
    warned = set()
    k0b = [el.name for el in elements if _tn(el) == 'Bend' and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0]
    if k0b and warnings is not None:
        warnings.append('K0 のある BEND %d 台 (挿入光源 (ウィグラー) など: %s ...): 内部の軌道の振動・端の非線形は入れていない (線形の光学は CALC の値)' % (len(k0b), k0b[0]))
    for i, el in enumerate(elements):
        L = float(getattr(el, 'L', 0.0) or 0.0)
        if L <= 0:
            continue
        r0 = rows[i]
        bx, ax, by, ay = float(r0.betx), float(r0.alfx), float(r0.bety), float(r0.alfy)
        gx0, gy0 = (1 + ax * ax) / bx, (1 + ay * ay) / by
        sl = _slice_for_optics(el, L / nsub)
        if sl is None:
            if _tn(el) in ('Solenoid', 'Undulator', 'TCavity') and warnings is not None and _tn(el) not in warned:
                warned.add(_tn(el))
                warnings.append('運動学的な項: %s は内部の光学を求めずドリフトとして扱う' % _tn(el))
            sxx += gx0 * gx0 * L
            sxy += gx0 * gy0 * L
            syy += gy0 * gy0 * L
            continue
        w = np.ones(nsub + 1)
        w[1:-1:2] = 4.0
        w[2:-1:2] = 2.0
        w *= L / nsub / 3.0
        gx = np.zeros(nsub + 1)
        gy = np.zeros(nsub + 1)
        for k in range(nsub + 1):
            if k == 0:
                gx[k], gy[k] = gx0, gy0
                continue
            M = np.asarray(_slice_for_optics(el, L * k / nsub).transfer_matrix(0.0), dtype=float)
            for (m11, m12, m21, m22, a0, b0), arr in (((M[0, 0], M[0, 1], M[1, 0], M[1, 1], ax, bx), gx), ((M[2, 2], M[2, 3], M[3, 2], M[3, 3], ay, by), gy)):
                g0 = (1 + a0 * a0) / b0
                beta = m11 * m11 * b0 - 2 * m11 * m12 * a0 + m12 * m12 * g0
                alf = -m11 * m21 * b0 + (m11 * m22 + m12 * m21) * a0 - m12 * m22 * g0
                arr[k] = (1 + alf * alf) / beta
        sxx += float(np.sum(w * gx * gx))
        sxy += float(np.sum(w * gx * gy))
        syy += float(np.sum(w * gy * gy))
    return {'xx': 3.0 / (16 * math.pi) * sxx, 'xy': 1.0 / (8 * math.pi) * sxy, 'yy': 3.0 / (16 * math.pi) * syy}


def quad_fringe_adts(elements, rows, warnings=None):
    """四極磁石と、K1≠0 の偏向磁石 (複合機能) の端部 (K64 の ttfrin4: Lee-Whiting の非線形のフリンジ、SAD の既定は ON = DISFRIN=0) の 1 次の ADTS。
    端部は H_f = σ (K1/12)[(x³ + 3xy²) px − (y³ + 3x²y) py] δ(s − s_edge) の時間 1 の流れ (σ = +1 入口、−1 出口。K1 は符号つきの勾配 [1/m²]。
    pysad の ttfrin4_1 の写像の 3 次の項と一致)。平均すると
      ⟨H_f⟩ = σ (K1/12)[−(3/2) αx βx Jx² + (3/2) αy βy Jy² + 3 (αy βx − αx βy) Jx Jy]
    から ∂νx/∂Jx = −(1/8π) Σ σ K1 αx βx、∂νx/∂Jy = (1/8π) Σ σ K1 (αy βx − αx βy)、∂νy/∂Jy = (1/8π) Σ σ K1 αy βy (α, β は端の値)。
    FRINGE (mfring) の指定 (入口だけ・出口だけ・なし) と DISFRIN は追跡と同じに読む。F1・F2 (線形のフリンジ) と K1 のある BEND の端は入っていない"""
    import math
    g = lambda el, k: float(getattr(el, k, 0.0) or 0.0)
    sxx = sxy = syy = 0.0
    for i, el in enumerate(elements):
        if _tn(el) not in ('Quadrupole', 'Bend'):
            continue
        L = float(el.L)
        k1 = float(el.K1)
        if L <= 0 or k1 == 0.0:
            continue
        if g(el, 'DISFRIN') != 0.0:
            continue
        mf = int(round(g(el, 'FRINGE')))
        if g(el, 'F1') or g(el, 'F2') or g(el, 'F1K1F') or g(el, 'F2K1F') or g(el, 'F1K1B') or g(el, 'F2K1B'):
            if warnings is not None:
                warnings.append('四極 %s: F1・F2 (線形のフリンジ) は入れていない' % el.name)
        for sig, r, on in ((1.0, rows[i], mf > -4 and mf != 2), (-1.0, rows[i + 1], mf > -4 and mf != 1)):
            if not on:
                continue
            ax, bx, ay, by = float(r.alfx), float(r.betx), float(r.alfy), float(r.bety)
            sxx += -sig * k1 * ax * bx
            sxy += sig * k1 * (ay * bx - ax * by)
            syy += sig * k1 * ay * by
    c = 1.0 / (8.0 * math.pi)
    return {'xx': c * sxx, 'xy': c * sxy, 'yy': c * syy}


# ---------------------------------------------------------------------------------------------
# クロマティシティ
# ---------------------------------------------------------------------------------------------
def chromaticity_sextupole(lat, method='analytic', nslice=8, rows=None):
    """6 極の 1 次のクロマティシティ ξx = h11001/π、ξy = h00111/π (= (1/4π)Σ K2L η β。η は分散、1 次の分散だけ)。戻り (ξx, ξy)"""
    r = resonance_terms(lat, ('h11001', 'h00111'), method=method, nslice=nslice, rows=rows, bend_cubic=False)
    return r.total['h11001'].real / math.pi, r.total['h00111'].real / math.pi


def _beta_integrals(el, row, nsub=8):
    """素子 el (入口の Twiss の行 row) の ∫βx ds、∫βy ds (線形の転送行列で伝播した β の Simpson 則)"""
    L = float(getattr(el, 'L', 0.0) or 0.0)
    if L <= 0:
        return 0.0, 0.0
    bx, ax, by, ay = float(row.betx), float(row.alfx), float(row.bety), float(row.alfy)
    w = np.ones(nsub + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    w *= L / nsub / 3.0
    bxs, bys = [bx], [by]
    for k in range(1, nsub + 1):
        sl = _slice_for_optics(el, L * k / nsub)
        if sl is None:
            s_ = L * k / nsub
            M = np.array([[1, s_, 0, 0], [0, 1, 0, 0], [0, 0, 1, s_], [0, 0, 0, 1]], dtype=float)
        else:
            M = np.asarray(sl.transfer_matrix(0.0), dtype=float)
        for (m11, m12, a0, b0), arr in (((M[0, 0], M[0, 1], ax, bx), bxs), ((M[2, 2], M[2, 3], ay, by), bys)):
            arr.append(m11 * m11 * b0 - 2 * m11 * m12 * a0 + m12 * m12 * (1 + a0 * a0) / b0)
    return float(np.sum(w * np.array(bxs))), float(np.sum(w * np.array(bys)))


def chromaticity_natural(lat, rows=None, nsub=8):
    """自然クロマティシティの解析式 (v2.3.0 R4): 四極・K1 つきの偏向磁石・MULT の K1 と偏向磁石の弱収束 h² と端の収束 (角度 E1・E2 [rad]、h tan E)
      ξx = −(1/4π) [∮ (K1 + h²) βx ds − Σ_端 h tan E βx]、  ξy = +(1/4π) [∮ K1 βy ds − Σ_端 h tan E βy]。
    (δ に対する 1/(1+δ) の依存は K1 と端の収束、h² は δ に依らない。**偏向磁石の 3 次の項を分散の軌道のずれで展開した分** (`chromaticity_bend`、h11001/π など)
    を足すと、CALC の δ 依存の 1 次 (6 極なし) と一致する)。K0≠0 の偏向磁石は飛ばす (v2.4.0 A2: `chromaticity_bend` の全体の値が自然な項を含むので二重に数えない)。入っていない: ROTATE した素子・F1/F2・FRINGE の指定・ソレノイド。戻り (ξx, ξy)"""
    elements, rows = _tables(lat, rows=rows)
    sx = sy = hx = ex = ey = 0.0
    for i, el in enumerate(elements):
        tn = _tn(el)
        L = float(getattr(el, 'L', 0.0) or 0.0)
        if L <= 0 or tn not in ('Quadrupole', 'Bend', 'Multipole'):
            continue
        if tn == 'Bend' and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0:
            continue      # v2.4.0 A2: K0≠0 の偏向磁石の自然な項 (K1・h²・端) は chromaticity_bend (element_delta_terms) が全体で持つ (二重に数えない)
        if tn == 'Multipole':
            knl = getattr(el, 'KnL', {}) or {}
            k1 = float(knl.get(1, 0.0) or 0.0) / L
        else:
            k1 = float(getattr(el, 'K1', 0.0) or 0.0)
        ang = float(getattr(el, 'ANGLE', 0.0) or 0.0) if tn == 'Bend' else 0.0
        if k1 == 0.0 and ang == 0.0:
            continue
        ix, iy = _beta_integrals(el, rows[i])
        sx += k1 * ix
        sy += k1 * iy
        if ang != 0.0:
            h = ang / L
            hx += h * h * ix
            for r, e in ((rows[i], float(getattr(el, 'E1', 0.0) or 0.0)), (rows[i + 1], float(getattr(el, 'E2', 0.0) or 0.0))):
                if e != 0.0:
                    t = h * math.tan(e)
                    ex += t * float(r.betx)
                    ey += t * float(r.bety)
    c = 1.0 / (4.0 * math.pi)
    return -c * (sx + hx - ex), c * (sy - ey)


def chromaticity_bend(lat, rows=None, nslice=None, momentum_eV=3e9):
    """偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2) の、分散による軌道のずれの寄与 (h11001/π、h00111/π の偏向磁石の分。追跡の写像から取り出す)。戻り (ξx, ξy)"""
    r = resonance_terms(lat, ('h11001', 'h00111'), rows=rows, bend_cubic=True, bend_nslice=nslice, momentum_eV=momentum_eV)
    return r.total_bend['h11001'].real / math.pi, r.total_bend['h00111'].real / math.pi


def chromaticity_numeric(lat, order=3, dmax=2e-3, ncheck=True):
    """CALC (4 次元の閉軌道 + 光学) の δ 依存から ∂^n ν/∂δ^n (n = 1..order) を数値微分で求める。δ = ±dmax·(Chebyshev 節点) の
    (2·order+1) 点を多項式でフィット。戻り dict: {'nu0': (νx,νy), 'd1': (..), 'd2': (..), ..., 'dmax': dmax}
    (チューンは小数部の連続値、閉軌道が見つからない δ では例外)。微分は ∂^n ν/∂δ^n (テイラー係数ではなく n! 倍した値)"""
    npts = 2 * order + 3
    nodes = dmax * np.cos(np.pi * (np.arange(npts) + 0.5) / npts)
    nus = []
    for d in nodes:
        t = lat.twiss_table(delta=float(d))[-1]
        nus.append((float(t.mux) / TWO_PI, float(t.muy) / TWO_PI))
    nus = np.array(nus)
    out = {'dmax': dmax, 'order': order}
    deg = order + 2
    co = [np.polynomial.polynomial.polyfit(nodes, nus[:, i], deg) for i in range(2)]
    out['nu0'] = tuple(float(c[0]) for c in co)
    for n in range(1, order + 1):
        out['d%d' % n] = tuple(float(c[n]) * math.factorial(n) for c in co)
    return out


# ---------------------------------------------------------------------------------------------
# 追跡による基準値
# ---------------------------------------------------------------------------------------------
def tracking_tunes(elements, momentum_eV, twiss0, jx, jy, nturn=1024, phix=0.0, phiy=0.0, delta=0.0, tracker=None):
    """アクション (jx, jy) の粒子を nturn 周追跡 (放射・RF なし) して NAFF で (νx, νy) (符号つき、端数) を返す。
    twiss0 = ((αx, βx), (αy, βy)) はリング先頭 (追跡の始点) の Twiss。位相 (phix, phiy) は x = √(2Jβ) cos φ、px = −√(2J/β)(sin φ + α cos φ)"""
    from .track6n import track_particles_fast
    from . import naff
    (ax, bx), (ay, by) = twiss0
    jx = np.atleast_1d(np.asarray(jx, dtype=float))
    jy = np.atleast_1d(np.asarray(jy, dtype=float))
    ph = np.atleast_1d(np.asarray(phix, dtype=float)) * np.ones_like(jx)
    pv = np.atleast_1d(np.asarray(phiy, dtype=float)) * np.ones_like(jx)
    n = len(jx)
    c = np.zeros((n, 6))
    c[:, 0] = np.sqrt(2 * jx * bx) * np.cos(ph)
    c[:, 1] = -np.sqrt(2 * jx / bx) * (np.sin(ph) + ax * np.cos(ph))
    c[:, 2] = np.sqrt(2 * jy * by) * np.cos(pv)
    c[:, 3] = -np.sqrt(2 * jy / by) * (np.sin(pv) + ay * np.cos(pv))
    c[:, 5] = delta
    rec = np.zeros((nturn, n, 6))
    tr = tracker or (lambda cc: track_particles_fast(elements, cc, 1, momentum_eV, 0.0, False, False, False))
    for t in range(nturn):
        rec[t] = c
        c = tr(c)
    out = []
    for i in range(n):
        (fx, fy, _fz), _a = naff.tunes_from_coords(rec[:, i, :], twiss=((ax, bx), (ay, by), (0.0, 1.0)))
        out.append((fx, fy))
    return np.array(out)


def adts_tracking(elements, momentum_eV, twiss0, jref=(1e-8, 1e-8), nturn=1024, delta=0.0, tracker=None):
    """追跡の NAFF のチューンから ADTS を求める (小振幅の極限)。基準点 (b, b)、b = s·jref (s = 1, 2) のまわりで、
    Jx を b ± b/2 (Jy は b)、Jy を b ± b/2 (Jx は b) と動かした 8 粒子の中心差分 f(s) = ∂ν/∂J |_(b,b) = q + c·s を作り、
    q = 2 f(1) − f(2) (J の 1 次の補正を消す)。どちらの振幅も同程度にするのは、片方を極小にすると 6 極の強制振動 (2νy など) が
    もう一方の面の NAFF のピークを乱すため。戻り: dict 'xx' = ∂νx/∂Jx、'xy' = ∂νx/∂Jy、'yx' = ∂νy/∂Jx、'yy' = ∂νy/∂Jy
    (jref = (jx0, jy0): 基準の振幅。2·jref まで使う)"""
    jx0, jy0 = jref
    f = []
    fr = lambda a: (a + 0.5) % 1.0 - 0.5
    for sc_ in (1.0, 2.0):
        bx, by = sc_ * jx0, sc_ * jy0
        hx, hy = bx / 2, by / 2
        jx = [bx + hx, bx - hx, bx, bx]
        jy = [by, by, by + hy, by - hy]
        nu = tracking_tunes(elements, momentum_eV, twiss0, jx, jy, nturn=nturn, delta=delta, tracker=tracker)
        f.append({'xx': fr(nu[0][0] - nu[1][0]) / (2 * hx), 'yx': fr(nu[0][1] - nu[1][1]) / (2 * hx),
                  'xy': fr(nu[2][0] - nu[3][0]) / (2 * hy), 'yy': fr(nu[2][1] - nu[3][1]) / (2 * hy)})
    out = {k: 2 * f[0][k] - f[1][k] for k in f[0]}
    out['jref'] = jref
    return out


# ---------------------------------------------------------------------------------------------
# 薄レンズ列への分割 (検証用)
# ---------------------------------------------------------------------------------------------
def slice_thin(lat, n, names=None):
    """厚い SEXT (と K2≠0 の厚い MULT) を n 個の薄い MULT (KnL[2] = K2L/n) とドリフトの列に置き換えた新しい Lattice (元は変えない)。
    線形光学は元と同じ (SEXT の内部は線形ではドリフト)。K1≠0 の厚い MULT は未対応 (例外)。n 個の薄レンズは等分した各区間の中央"""
    from .elements import Drift, Multipole
    from .lattice import Lattice
    out = []
    for el in lat.elements:
        tn = _tn(el)
        L = float(getattr(el, 'L', 0.0) or 0.0)
        if tn == 'Sextupole' and L > 0 and (names is None or el.name in names):
            k2l = float(el.K2) * L
            h = L / n
            for i in range(n):
                out.append(Drift(name=el.name + '_d%d' % i, L=h / 2 if True else h))
                m = Multipole(name=el.name + '_m%d' % i, L=0.0)
                m.KnL = {2: k2l / n}
                out.append(m)
                out.append(Drift(name=el.name + '_e%d' % i, L=h / 2))
        else:
            out.append(el)
    return Lattice(out, name=getattr(lat, 'name', 'LATTICE'))


# ---------------------------------------------------------------------------------------------
# SAD の台本から: FFS の関数
# ---------------------------------------------------------------------------------------------
def ffs_momentum(ffs):
    """FFS の運動量 [eV] (偏向磁石の写像を取り出す追跡に渡すだけ。4 次元・δ=0 の写像は運動量によらない)"""
    for k in ('_momentum_eV', 'momentum', 'p0'):
        v = getattr(ffs, k, None)
        if isinstance(v, (int, float)) and v > 0:
            return float(v)
    return 3e9


def install(ffs):
    """FFS の組み込み関数: ResonanceDrivingTerms[opts]、AmplitudeTuneShift[opts]、ChromaticityTerms[opts]"""
    from .sadscript import Rule, Sym, is_num
    from .sadscript import SADError
    B = ffs.I.builtins

    def _opts(args):
        o = {}
        for a in args:
            for r in (a if isinstance(a, list) else [a]):
                if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                    o[r.lhs.name] = r.rhs
        return o

    def _rows_lat():
        if ffs.lat is None:
            raise SADError('No beam line (USE first)')
        rows = ffs.rows if ffs.rows is not None else ffs.calc(quiet=True)
        return ffs.lat, list(rows)

    def _method(o):
        m = o.get('Method', 'Analytic')
        m = str(getattr(m, 'name', m)).strip('"').lower()
        return 'numeric' if m.startswith('n') else 'analytic'

    def _nslice(o):
        v = o.get('Slices', 8)
        return int(v) if is_num(v) else 8

    def _truth(v, default):
        if v is None:
            return default
        if isinstance(v, Sym):
            return v.name.upper() not in ('FALSE', 'OFF', 'NO')
        if is_num(v):
            return bool(v)
        return str(v).strip('"').lower() not in ('false', 'off', 'no', '0')

    def fn_rdt(args):
        o = _opts(args)
        lat, rows = _rows_lat()
        terms = o.get('Terms', None)
        if isinstance(terms, list) and terms:
            terms = [str(getattr(t, 'name', t)).strip('"') for t in terms]
        else:
            terms = DEFAULT_TERMS
        bc = _truth(o.get('BendCubic'), BEND_RDT_DEFAULT)
        bs = int(o['BendSlices']) if is_num(o.get('BendSlices')) else BEND_NSLICE_DEFAULT
        o2 = _truth(o.get('Order2'), ORDER2_DEFAULT)
        r = resonance_terms(lat, terms, method=_method(o), nslice=_nslice(o), rows=rows, bend_cubic=bc, bend_nslice=bs, momentum_eV=ffs_momentum(ffs), order2=o2)
        out = []
        for nm in r.names:
            v = r.total[nm]
            out.append(Rule(Sym(nm.upper()), [float(v.real), float(v.imag)]))
        return out

    def fn_adts(args):
        o = _opts(args)
        lat, rows = _rows_lat()
        bc = _truth(o.get('BendCubic'), BEND_CUBIC_DEFAULT)
        bs = int(o['BendSlices']) if is_num(o.get('BendSlices')) else BEND_NSLICE_DEFAULT
        fo = o.get('Fourth')
        fo = fo.lower() if isinstance(fo, str) else None      # v2.3.0 R6: Fourth->"r6" で精密版 (既定 "r3")
        r = adts(lat, method=_method(o), nslice=_nslice(o), rows=rows, bend_cubic=bc, bend_nslice=bs, momentum_eV=ffs_momentum(ffs), fourth=fo)
        return [Rule(Sym('DNuxDJx'), float(r.dnux_dJx)), Rule(Sym('DNuxDJy'), float(r.dnux_dJy)), Rule(Sym('DNuyDJy'), float(r.dnuy_dJy)),
                Rule(Sym('Sextupole'), [float(r.sext['xx']), float(r.sext['xy']), float(r.sext['yy'])]),
                Rule(Sym('Octupole'), [float(r.octu['xx']), float(r.octu['xy']), float(r.octu['yy'])]),
                Rule(Sym('Kinematic'), [float(r.kin['xx']), float(r.kin['xy']), float(r.kin['yy'])]),
                Rule(Sym('QuadFringe'), [float(r.fringe['xx']), float(r.fringe['xy']), float(r.fringe['yy'])]),
                Rule(Sym('BendFourth'), [float(r.bend4['xx']), float(r.bend4['xy']), float(r.bend4['yy'])]),
                Rule(Sym('BendCubic'), [float(r.bend['xx']), float(r.bend['xy']), float(r.bend['yy'])]),
                Rule(Sym('SextBendCross'), [float(r.cross['xx']), float(r.cross['xy']), float(r.cross['yy'])]),
                Rule(Sym('Total'), [float(r.total['xx']), float(r.total['xy']), float(r.total['yy'])])]

    def fn_chrom(args):
        o = _opts(args)
        lat, rows = _rows_lat()
        order = int(o['Order']) if is_num(o.get('Order')) else 3
        dmax = float(o['DeltaMax']) if is_num(o.get('DeltaMax')) else 2e-3
        c = chromaticity_numeric(lat, order=order, dmax=dmax)
        sx, sy = chromaticity_sextupole(lat, method=_method(o), nslice=_nslice(o), rows=rows)
        nx, ny = chromaticity_natural(lat, rows=rows)
        bx, by = chromaticity_bend(lat, rows=rows, momentum_eV=ffs_momentum(ffs))
        out = [Rule(Sym('Tunes'), [c['nu0'][0], c['nu0'][1]])]
        for n in range(1, order + 1):
            out.append(Rule(Sym('D%d' % n), [c['d%d' % n][0], c['d%d' % n][1]]))
        out.append(Rule(Sym('SextupolePart'), [sx, sy]))
        out.append(Rule(Sym('NaturalPart'), [nx, ny]))
        out.append(Rule(Sym('BendPart'), [bx, by]))
        return out

    B['ResonanceDrivingTerms'] = fn_rdt
    B['AmplitudeTuneShift'] = fn_adts
    B['ChromaticityTerms'] = fn_chrom
