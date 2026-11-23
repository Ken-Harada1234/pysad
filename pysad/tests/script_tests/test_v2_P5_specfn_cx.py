# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 3: 複素引数の特殊関数の k64 の挙動を gs と比べ、足りなかった関数を足した (pysad/specfn_cx2.py)。

走査: 12 個の複素数 (|z| 0.05〜22) × 60 の関数・形 = 732 式を gs と pysad (k64) で評価した (example/10_math_data/fn_SpecialFunctions_complex_scan.*)。
結果: 追加前は gs と 1e-9 で一致するのが 436 式。追加後は **641 式** (相対 1e-9)。残り 91 式は gs の値が誤る・不正確で、独立な基準 (scipy.special、級数) と比べて pysad が正しい:
  - gs の誤り (大きな差): LogGamma (Im z < 0 の Re z < 0 で虚部の分枝が 4π 違う)、Erf・Erfc (Re z < 0 で 2 − 真値など)、ExpIntegralEi (複素引数は全 12 点で誤り: D.n の式、
    実引数の B44 と同じ)、ArcCosh (複素数で符号が逆: 実数 0.5 の K64 の動作 −1.047 I と同じ向き)、FresnelS・FresnelC (|z| ≥ 12 で 1e38 など)、BesselY の複素引数 (gs は NaN)、
    HypergeometricU の整数 b (B120)・HurwitzZeta の Re s < −8 の点 (B118 の周辺)・PolyLog (B119)。
  - gs の精度の低さ: BesselK (|z| 大で 1e-9〜1e-1)、AiryAi・AiryBi (|z| 大)、ExpIntegralE[1,-20-10I] (4e-6)、Factorial (1e-9〜3e-8)。
  k64 モードは gs が誤る点でも真値を返す (docs/30 B118〜B120 と同じ扱い。gs の誤りの再現には K64 の gamma.f・D.n の直訳が要るので、要判断のまま)。
追加した関数 (実数引数は従来どおり。複素数が 1 つでもあれば複素の実装): ExpIntegralE[n,z] (整数 n)・ExpIntegralEi・PolyGamma・HarmonicNumber・ProductLog・SinIntegral・CosIntegral・
FresnelS・FresnelC・Beta・Factorial・GammaRegularized (P・Q)・Gamma[a,z]・LegendreP[ν,z]・WhittakerM・WhittakerW・BesselJ[ν,x] (複素の次数、|x| ≤ 25)。
従来の pysad は虚部を黙って捨てた実数を返す (ExpIntegralEi・SinIntegral・FresnelS/C・HarmonicNumber・LegendreP) か、何も返さなかった。
"""
import cmath
import math
import re

import numpy as np
import pytest

from . import _p5_gscx as G
from . import _p5_util as U

scipy = pytest.importorskip('scipy')
from scipy import integrate, special as sp                        # noqa: E402


def _evaluate(exprs, mode):
    src = 'FFS;\n$FORM="22.17";\n' + ''.join('Print["R%d= ",%s];\n' % (i, e) for i, e in enumerate(exprs)) + 'Exit;\n'
    out = U.run(src, mode)
    got = {}
    for l in out.splitlines():
        m = re.match(r'R(\d+)= (.*)', l)
        if m:
            got[int(m.group(1))] = m.group(2).strip()
    return got


def _val(s):
    if s is None:
        return None
    s = s.replace('INF', '1e999')
    if 'NaN' in s:
        return complex('nan')
    m = re.match(r'\(\s*(%s)\s*([-+])\s*(%s)\s*I\)' % (U.NUM, U.NUM), s)
    if m:
        return complex(float(m.group(1)), float(m.group(3)) * (1 if m.group(2) == '+' else -1))
    m = re.match(r'^\s*(%s)\s*$' % U.NUM, s)
    return complex(float(m.group(1)), 0.0) if m else None


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_added_complex_functions_match_gs(mode):
    """gs が正しく返す点 (追加した 15 関数) が gs と相対 1e-9 で一致 (k64・v2 とも)"""
    exprs = [e for e, _, _ in G.MATCH]
    got = _evaluate(exprs, mode)
    bad = []
    for i, (e, a, b) in enumerate(G.MATCH):
        v = _val(got.get(i))
        if v is None or abs(v - complex(a, b)) > 1e-9 * max(abs(complex(a, b)), 1e-300) + 1e-300:
            bad.append((e, (a, b), v))
    assert not bad, bad[:5]


def _ref(expr):
    """独立な基準 (scipy.special の複素対応)。gs が誤る式"""
    m = re.match(r'(\w+)\[(.*)\]$', expr)
    f, arg = m.group(1), m.group(2)
    z = complex(eval(arg.split(',')[-1].replace('I', 'j')))
    if f == 'LogGamma':
        return sp.loggamma(z)
    if f == 'Erf':
        return sp.erf(z)
    if f == 'Erfc':
        return sp.erfc(z)
    if f == 'ExpIntegralEi':
        return sp.expi(z)
    if f == 'ArcCosh':
        return cmath.acosh(z)
    if f == 'FresnelS':
        return sp.fresnel(z)[0]
    if f == 'FresnelC':
        return sp.fresnel(z)[1]
    if f == 'BesselY':
        return sp.yv(float(arg.split(',')[0]), z)
    return None


def test_gs_defects_pysad_k64_returns_the_true_value():
    """gs が誤る点 (DEFECT の LogGamma・Erf・Erfc・ExpIntegralEi・ArcCosh・Fresnel・BesselY): k64 の pysad は独立な基準の値 (gs の値とは違う)"""
    gsd = dict(G.DEFECT)
    exprs = [e for e, _ in G.DEFECT if _ref(e) is not None]
    got = _evaluate(exprs, 'k64')
    n_checked = 0
    for i, e in enumerate(exprs):
        r = _ref(e)
        v = _val(got.get(i))
        assert v is not None and abs(v - r) <= 1e-9 * abs(r) + 1e-300, (e, v, r)
        g = gsd[e]
        assert g is None or abs(complex(*g) - r) > 1e-6 * abs(r), (e, g, r)
        n_checked += 1
    assert n_checked >= 40


# ---------------------------------------------------------------- 独立な恒等式・微分方程式 (gs にも scipy の同じ関数にも頼らない)
Z = [0.3 + 0.4j, -1.2 + 0.7j, 2.1 - 1.3j, 0.9 - 2.2j, 5.0 + 1.0j]


def _lit(a):
    return ('(%r)' % a).replace('j', 'I') if isinstance(a, complex) else repr(a)


def _f(name, *args):
    expr = name + '[' + ','.join(_lit(a) for a in args) + ']'
    return _val(_evaluate([expr], 'v2')[0])


def _d12(fn, z, h=2e-4):
    f0, fp, fm = fn(z), fn(z + h), fn(z - h)
    return (fp - fm) / (2 * h), (fp - 2 * f0 + fm) / h ** 2


@pytest.mark.parametrize('z', Z)
def test_expintegral_ei_series_and_en_recurrence(z):
    ei = 0.5772156649015329 + cmath.log(z)
    t, term = 0j, 1 + 0j
    for k in range(1, 100):
        term *= z / k
        t += term / k
    assert _f('ExpIntegralEi', z) == pytest.approx(ei + t, rel=1e-11)
    e1, e2, e3 = _f('ExpIntegralE', 1, z), _f('ExpIntegralE', 2, z), _f('ExpIntegralE', 3, z)
    assert e2 == pytest.approx(cmath.exp(-z) - z * e1, rel=1e-11)
    assert 2 * e3 == pytest.approx(cmath.exp(-z) - z * e2, rel=1e-10)
    d = (_f('ExpIntegralE', 2, z + 1e-4) - _f('ExpIntegralE', 2, z - 1e-4)) / 2e-4
    assert d == pytest.approx(-e1, rel=1e-6)                                  # dE_n/dz = -E_{n-1}


@pytest.mark.parametrize('z', Z)
def test_polygamma_harmonic_beta_factorial_identities(z):
    psi = lambda w: _f('PolyGamma', w)
    assert psi(z + 1) == pytest.approx(psi(z) + 1 / z, rel=1e-11)
    n = 200000                                                                # ψ'(z) = Σ 1/(z+k)^2 (n 項 + 尾の補正)
    k = np.arange(n)
    s1 = np.sum(1.0 / (z + k) ** 2) + 1.0 / (z + n) + 0.5 / (z + n) ** 2
    assert _f('PolyGamma', 1, z) == pytest.approx(s1, rel=1e-9)
    assert _f('HarmonicNumber', z) == pytest.approx(psi(z + 1) + 0.5772156649015329, rel=1e-12)
    assert _f('HarmonicNumber', z) - _f('HarmonicNumber', z - 1) == pytest.approx(1 / z, rel=1e-10)
    f = lambda w: _f('Factorial', w)
    assert f(z + 1) == pytest.approx((z + 1) * f(z), rel=1e-11)
    assert _f('Beta', z, 2.5) == pytest.approx(_f('Gamma', z) * math.gamma(2.5) / _f('Gamma', z + 2.5), rel=1e-11)
    assert _f('Beta', z, 1.0) == pytest.approx(1 / z, rel=1e-11)


def _line_integral(g, z):
    """∫_0^1 z g(s z) ds (原点から z への直線の複素の線積分)。実部・虚部を別に求積する"""
    re_ = integrate.quad(lambda s: (z * g(s * z)).real, 0, 1, epsabs=0, epsrel=1e-13, limit=200)[0]
    im_ = integrate.quad(lambda s: (z * g(s * z)).imag, 0, 1, epsabs=0, epsrel=1e-13, limit=200)[0]
    return complex(re_, im_)


@pytest.mark.parametrize('z', Z)
def test_productlog_sici_fresnel_incomplete_gamma(z):
    w = _f('ProductLog', z)
    assert w * cmath.exp(w) == pytest.approx(z, rel=1e-12)
    si = sum((-1) ** k * z ** (2 * k + 1) / ((2 * k + 1) * math.factorial(2 * k + 1)) for k in range(40))
    ci = 0.5772156649015329 + cmath.log(z) + sum((-1) ** k * z ** (2 * k) / (2 * k * math.factorial(2 * k)) for k in range(1, 40))
    assert _f('SinIntegral', z) == pytest.approx(si, rel=1e-10)
    assert _f('CosIntegral', z) == pytest.approx(ci, rel=1e-10)
    assert _f('FresnelS', z) == pytest.approx(_line_integral(lambda t: cmath.sin(math.pi * t * t / 2), z), rel=1e-9)
    assert _f('FresnelC', z) == pytest.approx(_line_integral(lambda t: cmath.cos(math.pi * t * t / 2), z), rel=1e-9)
    a = 2.5                                                                   # γ(a,z) = ∫_0^z t^{a-1} e^{-t} dt
    low = _line_integral(lambda t: t ** (a - 1) * cmath.exp(-t), z)
    assert _f('GammaRegularizedP', a, z) == pytest.approx(low / math.gamma(a), rel=1e-9)
    assert _f('GammaRegularized', a, z) + _f('GammaRegularizedP', a, z) == pytest.approx(1.0, abs=1e-12)
    assert _f('Gamma', a, z) == pytest.approx(math.gamma(a) - low, rel=1e-9)


@pytest.mark.parametrize('z', [0.3 + 0.4j, -1.2 + 0.7j, 2.1 - 1.3j])
def test_legendre_whittaker_bessel_satisfy_their_odes(z):
    nu = 2.5                                                                  # (1 - z²) P'' - 2 z P' + ν(ν+1) P = 0
    d1, d2 = _d12(lambda w: _f('LegendreP', nu, w), z)
    p = _f('LegendreP', nu, z)
    assert (1 - z * z) * d2 - 2 * z * d1 + nu * (nu + 1) * p == pytest.approx(0, abs=3e-6 * abs(p))
    k, m = 0.5, 1.5                                                           # W'' + (-1/4 + κ/z + (1/4 - μ²)/z²) W = 0
    for name in ('WhittakerM', 'WhittakerW'):
        d1, d2 = _d12(lambda w: _f(name, k, m, w), z)
        w = _f(name, k, m, z)
        assert d2 + (-0.25 + k / z + (0.25 - m * m) / z ** 2) * w == pytest.approx(0, abs=3e-6 * abs(w))
    nu = 0.3 + 0.4j                                                           # 複素の次数: z² J'' + z J' + (z² - ν²) J = 0
    d1, d2 = _d12(lambda w: _f('BesselJ', nu, w), z)
    j = _f('BesselJ', nu, z)
    assert z * z * d2 + z * d1 + (z * z - nu * nu) * j == pytest.approx(0, abs=3e-6 * abs(j))


def test_complex_order_bessel_reduces_to_the_real_order_value():
    """複素の次数の実装 (級数) は、次数の虚部を 0 に近づけた極限で実数の次数 (scipy) に一致"""
    expr = 'BesselJ[(1.5+1e-12 I),2.5]'
    v = _val(_evaluate([expr], 'v2')[0])
    assert v.real == pytest.approx(float(sp.jv(1.5, 2.5)), rel=1e-9)


def test_real_arguments_are_unchanged():
    """実数の引数は従来の実装のまま (ExpIntegralEi[1.5] = 3.3013、実数の Beta・Gamma[a,z])"""
    got = _evaluate(['ExpIntegralEi[1.5]', 'Beta[2.5,3.5]', 'Gamma[2.5,1.5]', 'SinIntegral[2]', 'FresnelS[1.2]', 'LegendreP[2.5,.3]'], 'k64')
    ref = [3.301285449130, 0.036815538909, math.gamma(2.5) * float(sp.gammaincc(2.5, 1.5)), 1.605412976803, 0.623400918546, -0.463661013457]
    for i, r in enumerate(ref):
        assert _val(got[i]).real == pytest.approx(r, rel=1e-9)
        assert _val(got[i]).imag == 0.0
