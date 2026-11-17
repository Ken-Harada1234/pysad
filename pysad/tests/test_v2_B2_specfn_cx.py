# -*- coding: utf-8 -*-
"""複素引数の Zeta・HurwitzZeta・PolyLog・HypergeometricU (pysad/specfn_cx.py。v2 の追加で、k64 モードでも同じ値)。

根拠は mpmath ではない独立な計算 (この環境に mpmath は無い):
  - Zeta: Borwein のアルゴリズム 2 (交代級数) と関数等式。HurwitzZeta: 別に書いた Euler–Maclaurin (固定の 120 項)
  - PolyLog: |z| < 0.4 は冪級数、s = 2 は scipy.special.spence の複素版、s > 0 は積分表示 z/Γ(s) ∫ t^{s-1}/(e^t - z) dt (数値積分)、
    整数次の |z| 大は Bernoulli 多項式の反転公式
  - HypergeometricU: 積分表示 (被積分関数が減衰する向きに回した経路。Re a > 0)、接続公式 Γ(1-b)/Γ(a-b+1) M(a,b,z) + Γ(b-1)/Γ(a) z^{1-b} M(a-b+1,2-b,z)
    (scipy の hyp1f1)、a が負の整数の多項式、漸近展開 (|z| 大)
gs (SAD K64、2026-10-19 ビルド) の値は、gs が正しい点 (相対 1e-14 で一致) と gs が誤る点 (下の `test_gs_defects_*`) を分けて書いた。
使った範囲: s・z・a・b は |s| ≤ 100、|z| ≤ 1e8 (PolyLog)、|z| ≤ 30 (HypergeometricU)。それ以外は確かめていない。
"""
import cmath
import io
import math
import re
import warnings

import pytest

scipy = pytest.importorskip('scipy')
from scipy import integrate, special as sp                       # noqa: E402

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _quad_c(f, a, b):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        re_ = integrate.quad(lambda t: f(t).real, a, b, epsabs=0, epsrel=1e-13, limit=500)[0]
        im_ = integrate.quad(lambda t: f(t).imag, a, b, epsabs=0, epsrel=1e-13, limit=500)[0]
    return complex(re_, im_)


def ref_zeta(s, n=60):
    s = complex(s)
    if s.real < 0.5:
        return 2**s * math.pi**(s - 1) * cmath.sin(math.pi * s / 2) * complex(sp.gamma(1 - s)) * ref_zeta(1 - s, n)
    d = []
    for i in range(n + 1):
        t = 0.0
        for k in range(i + 1):
            t += math.factorial(n + k - 1) * 4**k / (math.factorial(n - k) * math.factorial(2 * k))
        d.append(n * t)
    tot = sum((-1)**k * (d[k] - d[n]) / (k + 1)**s for k in range(n))
    return -tot / (d[n] * (1 - 2**(1 - s)))


def ref_hurwitz(s, a, N=120):
    s, a = complex(s), complex(a)
    head = 0j
    while a.real < 0.5:
        head += cmath.exp(-s * cmath.log(a))
        a += 1
    tot = head + sum((a + k)**(-s) for k in range(N))
    x = a + N
    tot += x**(1 - s) / (s - 1) + 0.5 * x**(-s)
    B = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6, -3617 / 510, 43867 / 798, -174611 / 330]
    p, xp = s, x**(-s - 1)
    for j in range(1, len(B) + 1):
        tot += B[j - 1] / math.factorial(2 * j) * p * xp
        p *= (s + 2 * j - 1) * (s + 2 * j)
        xp /= x * x
    return tot


def ref_polylog(s, z):
    z = complex(z)
    if abs(z) < 0.4:
        tot, t = 0j, z
        for k in range(1, 4000):
            tot += t / k**s
            t *= z
            if abs(t) < 1e-18:
                break
        return tot
    if s == 2.0:
        return complex(sp.spence(1 - z))
    if s == 1.0:
        return -cmath.log(1 - z)
    lg = math.lgamma(s)

    def f(t):
        return cmath.exp((s - 1) * math.log(t) - lg) / (cmath.exp(t) - z) if 0 < t < 700 else 0j
    lz = math.log(abs(z))
    tot = _quad_c(f, 0.0, 50.0 + 10 * s) + _quad_c(f, 50.0 + 10 * s, 700.0)
    return z * tot


def ref_hyperu_int(a, b, z):
    z = complex(z)
    phi = -cmath.phase(z)
    e = cmath.exp(1j * phi)

    def f(r):
        if r <= 0:
            return 0j
        t = r * e
        return cmath.exp(-z * t + (a - 1) * cmath.log(t) + (b - a - 1) * cmath.log(1 + t)) * e
    tot = _quad_c(f, 0.0, 1.0) + _quad_c(f, 1.0, 40.0) + _quad_c(f, 40.0, 800.0 / abs(z))
    return tot / complex(sp.gamma(a))


def ref_hyperu_conn(a, b, z):
    z = complex(z)
    return complex(sp.gamma(1 - b) / sp.gamma(a - b + 1) * sp.hyp1f1(a, b, z)
                   + sp.gamma(b - 1) / sp.gamma(a) * cmath.exp((1 - b) * cmath.log(z)) * sp.hyp1f1(a - b + 1, 2 - b, z))


def _rel(a, b):
    return abs(complex(a) - complex(b)) / max(abs(complex(b)), 1e-300)


# ---------------------------------------------------------------- 直接の呼び出し
def test_zeta_complex_vs_borwein():
    from pysad import specfn_cx as X
    for s in [2 + 3j, 0.5j, -1 + 2j, 1 + 1j, 2j, 3 + 0.001j, -5 + 1j, -20.5 + 3j, 10 + 50j, 50 + 1j, 1.5 - 7j, -2.5 + 30j]:
        assert _rel(X.zeta_c(s), ref_zeta(s)) < 1e-11, s
    assert abs(X.zeta_c(0.5 + 14.134725141734693j)) < 1e-13                   # 臨界線上の最初の零点
    assert _rel(X.zeta_c(0.2 + 100j), ref_hurwitz(0.2 + 100j, 1, N=300)) < 1e-11    # Borwein (n=60) が |Im s| = 100 で足りない所は Euler–Maclaurin (N = 300)


def test_hurwitz_complex_vs_euler_maclaurin():
    from pysad import specfn_cx as X
    for s, a in [(2 + 1j, .5), (2 + 1j, 3.5), (3 + 2j, 10), (.5 + 1j, 1), (2 + 1j, 1 + 1j), (2, .5 + 1j), (-1 + 1j, 1.5), (1.5 + 5j, 0.2 + 0.1j), (3, 2 - 2j)]:
        assert _rel(X.hurwitz_c(s, a), ref_hurwitz(s, a)) < 1e-10, (s, a)
    assert _rel(X.hurwitz_c(2 + 1j, .5), (2**(2 + 1j) - 1) * ref_zeta(2 + 1j)) < 1e-11       # ζ(s, 1/2) = (2^s - 1) ζ(s)
    # a が負 (主値の冪)
    v = X.hurwitz_c(2 + 1j, -.5 + .3j)
    assert _rel(v, ref_hurwitz(2 + 1j, -.5 + .3j)) < 1e-10


def test_polylog_complex_vs_independent():
    from pysad import specfn_cx as X
    zs = [.5 + .5j, 2 + 1j, -3 + 2j, .3 + .4j, -.5 - .5j, .5 - 1j, 2 - .1j, 2 + .1j, .9j, 1j, .99 + .05j, 1.01 + .01j, -1 + 1e-9j, .7 * cmath.exp(2.5j),
          1.2 * cmath.exp(-2.9j), 5 + 3j, -1.2 - .2j]
    for s in [2, 3, 2.5, .5, 1.5, 1, 4, 7.3]:
        for z in zs:
            if abs(z) > 3 and s not in (1, 2):
                continue                                                          # 参照の積分は |z| が大きいと不安定
            assert _rel(X.polylog_c(s, z), ref_polylog(s, z)) < 2e-11, (s, z)
    # |z| 大: 整数次は反転公式 Li_n(z) = Li_n(1/z) - (2πi)^n/n! B_n(1/2 + ln(-z)/(2πi)) (参照は 1/z の冪級数)
    from pysad.specfn import _bernoulli
    for n in (2, 3, 4, 5):
        for z in (5 + 3j, -20 + 1j, 100j, -1000 + 500j, 1e4 + 1e3j, -1e8 + 1j):
            x = .5 + cmath.log(-z) / (2j * math.pi)
            Bn = sum(math.comb(n, k) * float(_bernoulli(k)) * x**(n - k) for k in range(n + 1))
            ref = ref_polylog(n, 1 / z) - (-1)**n * 0 - (2j * math.pi)**n / math.factorial(n) * Bn * (1 if n % 2 == 0 else 1)
            ref = (-(-1)**n * ref_polylog(n, 1 / z)) - (2j * math.pi)**n / math.factorial(n) * Bn
            assert _rel(X.polylog_c(n, z), ref) < 1e-11, (n, z)
    # s が 0 以下の整数: 有理関数
    z = 2 + 1j
    assert _rel(X.polylog_c(0, z), z / (1 - z)) < 1e-15
    assert _rel(X.polylog_c(-1, z), z / (1 - z)**2) < 1e-15
    z = .3 + .4j
    assert _rel(X.polylog_c(-3, z), sum(k**3 * z**k for k in range(1, 400))) < 1e-12        # Li_{-3}(z) = Σ k³ z^k


def test_hyperu_complex_vs_independent():
    from pysad import specfn_cx as X
    for a, b, z in [(1, 2, .5 + .5j), (.5, 1.5, 1 + 1j), (.5, 1.5, -1 + .1j), (.5, 1.5, -1 - .1j), (2, 3, -.5 + .5j), (.5, 1.5, 3j), (1.5, 3, 20 + 5j), (2.5, 1.2, .3 + .2j),
                    (.2, 3.3, -2 + 2j), (3, 1, 2 + 1j), (.05, 1.5, 1 + 2j), (1, 1, 5 + 0j), (4.5, .5, -3 + 4j), (.5, .5, 30 + 30j), (1.5, 2.5, -4 + 1e-3j)]:
        assert _rel(X.hyperu_c(a, b, z), ref_hyperu_int(a, b, z)) < 1e-10, (a, b, z)
    # a が 0 以下 (非整数): 接続公式 (b が整数でない)
    for a, b, z in [(-1.3, 2.5, .5 + .5j), (-.7, .4, 2 + 1j), (-2.5, -.6, 1 - 1j), (-3.4, 2.2, -1.5 + .5j), (-.5, 1.7, 1 + 1j), (-.5, 1.7, -1 + 1j), (-2.2, 3.3, 4 + 2j)]:
        assert _rel(X.hyperu_c(a, b, z), ref_hyperu_conn(a, b, z)) < 1e-10, (a, b, z)
    # a - b + 1 が 0 以下の整数 (接続公式の Γ が極): b を ±1e-3 振った値の平均 (U は b について解析的。誤差 O(δ²))
    for a, b, z in [(-.5, 1.5, 1 + 1j), (-.7, .3, 2 + 1j), (-2.5, -.5, 1 - 1j), (-.5, 1.5, -1 + 1j)]:
        d = 1e-3
        ref = 0.5 * (ref_hyperu_conn(a, b + d, z) + ref_hyperu_conn(a, b - d, z))
        assert _rel(X.hyperu_c(a, b, z), ref) < 1e-5, (a, b, z)
    # a が負の整数: U(-2, 3/2, z) = z² - 5z + 15/4
    z = 2 + 1j
    assert abs(X.hyperu_c(-2, 1.5, z) - (z * z - 5 * z + 3.75)) < 1e-13
    # U(a, a+1, z) = z^{-a}、U(1, 2, z) = 1/z (b が整数)
    assert abs(X.hyperu_c(1, 2, .5 + .5j) - 1 / (.5 + .5j)) < 1e-12
    assert abs(X.hyperu_c(2, 3, -.5 + .5j) - (-.5 + .5j)**-2) < 1e-12
    # 負の実数の z は上側 (arg z = +π): U(1/2, 3/2, -2) = (-2)^{-1/2} = -i/√2
    assert abs(X.hyperu_c(.5, 1.5, -2.0) + 1j / math.sqrt(2.0)) < 1e-14


def test_hyperu_large_argument_matches_asymptotic_series():
    from pysad import specfn_cx as X
    a, b, z = 1.5, 3.0, 20 + 5j                 # 漸近展開 z^{-a} Σ (a)_k (a-b+1)_k / k! (-z)^-k (打ち切りの誤差 ~1e-9 の範囲で)
    tot, t = 0, 1
    for k in range(40):
        tot += t
        t *= -(a + k) * (a - b + 1 + k) / ((k + 1) * z)
        if abs(t) < 1e-20:
            break
    asym = cmath.exp(-a * cmath.log(z)) * tot
    v = X.hyperu_c(a, b, z)
    assert abs(v - asym) / abs(asym) < 1e-7
    assert _rel(v, ref_hyperu_int(a, b, z)) < 1e-11


# ---------------------------------------------------------------- SADScript から (gs の値)
def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=300;\n' + body + '\nExit;\n', out=buf)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _cx(line, tag):
    for l in line:
        if l.startswith(tag + ' '):
            t = l[len(tag) + 1:].strip()
            if t.startswith('('):
                t = t[1:-1].strip()
                assert t.endswith('I')
                t = t[:-1].strip()
                i = max(t.rfind('+'), t.rfind('-'))
                while i > 0 and t[i - 1] in 'eE':
                    i = max(t.rfind('+', 0, i - 1), t.rfind('-', 0, i - 1))
                im = t[i:]
                if im in ('+', '-'):
                    im += '1'
                return complex(float(t[:i]) if t[:i] else 0.0, float(im))
            return complex(float(t))
    raise AssertionError('no %r in %r' % (tag, line))


# gs で正しい点 (gs の出力。pysad は相対 1e-13 以内で一致する)
_GS_OK = {
    'P1': ('PolyLog[2,0.5+0.5*I]', .4539852691502954 + .6437673328892692j),
    'P2': ('PolyLog[2,2+I]', 1.1866885370000568 + 2.4077407693457724j),
    'P3': ('PolyLog[3,-3+2*I]', -2.466410148606688 + 1.271214825519137j),
    'P4': ('PolyLog[2.5,0.3+0.4*I]', .2786247499211948 + .4435630363020921j),
    'P5': ('PolyLog[2,-0.5-0.5*I]', -.4834728063610784 - .4001561537420078j),
    'P7': ('PolyLog[2,2-0.1*I]', 2.312941798841537 - 2.181342177469578j),
    'P8': ('PolyLog[2,2+0.1*I]', 2.312941798841535 + 2.1813421774695785j),
    'P10': ('PolyLog[0.5,0.9*I]', -.3730455882234203 + .6371956283180455j),
    'Z1': ('Zeta[2+3*I]', .7980219851462758 - .1137443080529386j),
    'Z3': ('Zeta[-1+2*I]', .16891566977083436 - .07051598890825445j),
    'Z4': ('Zeta[0.5*I]', -.2991978789224702 - .35937117926480217j),
    'Z7': ('Zeta[2*I]', .314725764042099 - .2316796487505208j),
    'H2': ('HurwitzZeta[2,0.5+I]', .03672455194101715 - 1.1170686578295996j),
    'U2': ('HypergeometricU[0.5,1.5,1+I]', .7768869870150186 - .3217971264527913j),
    'U3': ('HypergeometricU[0.5,1.5,-1+0.1*I]', .04968994017403841 - .9962771201104087j),
    'U5': ('HypergeometricU[0.5,1.5,-2]', -.7071067811865475j),
    'U6': ('HypergeometricU[-2,1.5,2+I]', -3.25 - 1j),
    'U8': ('HypergeometricU[0.5,1.5,3*I]', .4082482904638631 - .408248290463863j),
}


def test_complex_values_agree_with_gs_where_gs_is_right():
    r = _run('\n'.join('Print["%s ",%s];' % (k, e) for k, (e, _v) in _GS_OK.items()))
    for k, (_e, v) in _GS_OK.items():
        assert _rel(_cx(r, k), v) < 2e-13, (k, _cx(r, k), v)


def test_gs_defects_are_not_reproduced():
    """gs (K64) が誤る点: pysad は独立な参照と一致する値を返す (gs の値は下の通りで、相対誤差は大きい)"""
    r = _run('Print["H1 ",HurwitzZeta[2+I,0.5]];Print["H3 ",HurwitzZeta[2+I,3.5]];Print["U1 ",HypergeometricU[1,2,0.5+0.5*I]];'
             'Print["U9 ",HypergeometricU[2,3,-0.5+0.5*I]];Print["L1 ",PolyLog[2.5,-1.2-0.2*I]];Print["L2 ",PolyLog[1.5,-3.2843023372255384-3.5397879004724295*I]];'
             'Print["U10 ",HypergeometricU[1.5,3,20+5*I]];')
    assert _rel(_cx(r, 'H1'), (2**(2 + 1j) - 1) * ref_zeta(2 + 1j)) < 1e-11            # gs: 3.9548810631401707+2.179111683655945 I (相対 12% の誤り)
    assert _rel(_cx(r, 'H3'), ref_hurwitz(2 + 1j, 3.5)) < 1e-10                         # gs: .3721142639290967-.07448815908916451 I
    assert _rel(_cx(r, 'U1'), 1 / (.5 + .5j)) < 1e-12                                   # gs: NaN (b が整数)
    assert _rel(_cx(r, 'U9'), 2j) < 1e-12                                               # gs: NaN
    assert _rel(_cx(r, 'L1'), ref_polylog(2.5, -1.2 - .2j)) < 1e-10                     # gs: 17.81-15.59 I (非整数次の第 3 象限・|z| > 1 の誤り)
    assert _rel(_cx(r, 'L2'), ref_polylog(1.5, -3.2843023372255384 - 3.5397879004724295j)) < 1e-10   # gs: 3.466+9.430 I
    assert _rel(_cx(r, 'U10'), ref_hyperu_int(1.5, 3, 20 + 5j)) < 1e-11                 # gs: .010280882117903572 (相対 1e-8 の誤り)


def test_real_arguments_are_unchanged():
    r = _run('Print["a ",Zeta[3]];Print["b ",PolyLog[2,0.5]];Print["c ",HurwitzZeta[2,0.5]];Print["d ",HypergeometricU[1,2,0.5]];Print["e ",PolyLog[2,2]];')
    v = {l[0]: l[2:] for l in r if l and l[0] in 'abcde'}
    assert abs(float(v['a']) - 1.2020569031595942) < 1e-14
    assert abs(float(v['b']) - (math.pi**2 / 12 - math.log(2)**2 / 2)) < 1e-14
    assert abs(float(v['c']) - math.pi**2 / 2) < 1e-13
    assert abs(float(v['d']) - 2.0) < 1e-12
    assert 'I' in v['e']                                                                # 実軸上の z > 1 は従来どおり (切断線の下側)
