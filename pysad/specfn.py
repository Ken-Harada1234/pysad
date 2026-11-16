"""
特殊関数・数学関数。K64 の組み込み (src/tfinitn.f の itfunaloc、tfunc*.f) と Packages/D.n の SADScript 定義の置き換え。

SAD の実数引数の結果を scipy.special で計算する (K64 は Fortran の自前実装)。値は標準の数学的な値で、
K64 が誤った値を返す箇所 (docs/19 の B43〜) はそちらに合わせず、正しい値を返す:
  - ExpIntegralEi: D.n の式が実引数で複素数の誤値になる (Ei[1.5] → .164-2.985 I; 正しくは 3.3012)
  - HermiteH[3,.7] (整数 n、Gamma[(1-n)/2] の極で NaN)、LegendreP[n,m,z] (整数 m で NaN)、Binomial[-3,2] (NaN)
  - HypergeometricU[1,2,.5] (INF; 正しくは 2)、BesselY[0,.1] (NaN)、BesselK[1.5,.7] (INF)、BesselJ[-1,1.5]/BesselY[1,2.5] (gs が応答しない)
  - BesselJZero[0,1] (gs が応答しない)
対応: Arc 双曲線関数、LogGamma、Gamma (複素)、GammaRegularized/P/Q、PolyGamma、Beta、Zeta、HurwitzZeta、PolyLog、LerchPhi、BernoulliB、ProductLog、
  Pochhammer、超幾何関数 (0F1 1F1 2F1 U PFQ とその Regularized)、Bessel J/Y/I/K、D.n の特殊関数 (ExpIntegralE/Ei, LogIntegral, Sin/Cos/Sinh/CoshIntegral,
  Whittaker M/W, HermiteH, Fresnel S/C, Airy, LegendreP, BesselJZero, Binomial, Multinomial, InverseErfc, HarmonicNumber, Factorial2, XSin, Gamma0)。
未対応: 複素引数の多くの関数 (Bessel、超幾何、PolyLog)、PolyLog の |z|>=150 (非整数次)、Derivative の記号微分 (D.n の D)。
"""
from __future__ import annotations

import os
import cmath
import math
from fractions import Fraction

import numpy as np

from .sadscript import SADError, is_num, INF, NAN, Sym, Expr, S

try:
    from scipy import special as sp
except ImportError:  # pragma: no cover
    sp = None

EULER = 0.5772156649015328606


def _thread(f, args):
    """リストの引数があれば要素ごとに適用 (長さが違えば SADError)"""
    ls = [a for a in args if isinstance(a, list)]
    if not ls:
        return f(*args)
    n = len(ls[0])
    for l in ls:
        if len(l) != n:
            raise SADError('General::wrongleng: list length mismatch')
    return [_thread(f, [a[k] if isinstance(a, list) else a for a in args]) for k in range(n)]


def _numeric(x):
    if isinstance(x, list):
        return all(_numeric(t) for t in x)
    return is_num(x)


def _real(z):
    """虚部が 0 の複素数は実数にする"""
    if isinstance(z, complex) and z.imag == 0.0:
        return z.real
    return z


def _zeta_real(s):
    if s == 1.0:
        return -INF
    if s > 1.0:
        return float(sp.zeta(s, 1.0))
    return float(sp.zetac(s)) + 1.0


def _hurwitz(s, q):
    if q <= 0.0:
        if q == math.floor(q):
            return NAN
        # ζ(s, q) = Σ_{n<m} (n+q)^{-s} + ζ(s, q+m)。負の (n+q)^{-s} は主値 (gs: ζ(2.5, -0.5) = 6.2471 - 5.6569 i)
        m = int(math.ceil(-q)) + 1
        tot = 0j
        for n in range(m):
            tot += complex(n + q) ** (-s)
        tot += _hurwitz(s, q + m)
        return tot if abs(tot.imag) > 1e-12 * max(abs(tot.real), 1e-300) else tot.real
    if s > 1.0:
        return float(sp.zeta(s, q))
    # Euler-Maclaurin
    n = 24
    tot = sum((k + q) ** (-s) for k in range(n))
    x = n + q
    tot += x ** (1.0 - s) / (s - 1.0) + 0.5 * x ** (-s)
    fact = 1.0
    poch = s
    for j in range(1, 13):
        b2j = float(_bernoulli(2 * j))
        fact = math.factorial(2 * j)
        # (s)_{2j-1}
        p = 1.0
        for k in range(2 * j - 1):
            p *= (s + k)
        tot += b2j / fact * p * x ** (-s - 2 * j + 1)
    return tot


_BCACHE = {}


def _bernoulli(n):
    """ベルヌーイ数 B_n (B_1 = -1/2) を Fraction で"""
    if n in _BCACHE:
        return _BCACHE[n]
    a = [Fraction(0)] * (n + 1)
    for m in range(n + 1):
        a[m] = Fraction(1, m + 1)
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
        _BCACHE[m] = -a[0] if m == 1 else a[0]
    return _BCACHE[n]


def _bernoulli_poly(n, x):
    return sum(math.comb(n, k) * float(_bernoulli(k)) * x ** (n - k) for k in range(n + 1))


def _loggamma(x):
    if isinstance(x, complex):
        return complex(sp.loggamma(x))
    if x > 0.0:
        return math.lgamma(x)
    if x == math.floor(x):
        return INF
    n = int(math.ceil(-x))
    r = complex(math.lgamma(x + n))
    for k in range(n):
        r -= cmath.log(complex(x + k))
    return r


def _gamma(x):
    if isinstance(x, complex):
        return _real(complex(sp.gamma(x)))
    if x == math.floor(x) and x <= 0.0:
        return INF if x == 0.0 else NAN
    try:
        return math.gamma(x)
    except OverflowError:
        return INF


def _polygamma(n, z):
    n = int(n)
    z = float(z)
    if z <= 0.0 and z == math.floor(z):
        return NAN
    if n == 0:
        return float(sp.digamma(z))
    # ψ^(n)(z) = ψ^(n)(z+1) + (-1)^(n+1) n! / z^(n+1)
    tot = 0.0
    while z < 1.0:
        tot += (-1) ** (n + 1) * math.factorial(n) / z ** (n + 1)
        z += 1.0
    return float(sp.polygamma(n, z)) + tot


def _polylog(s, z):
    if isinstance(z, complex):
        z = _real(z)
        if isinstance(z, complex):
            raise SADError('PolyLog: complex z is not supported')
    if z == 0.0:
        return 0.0
    if z == 1.0:
        return _zeta_real(s)
    if z == -1.0:
        return -(1.0 - 2.0 ** (1.0 - s)) * _zeta_real(s)
    if s == 0.0:
        return z / (1.0 - z)
    if s == 1.0:
        return -math.log(1.0 - z) if z < 1.0 else NAN
    if s == 2.0 and z <= 1.0:
        return float(sp.spence(1.0 - z))
    if abs(z) <= 0.8:
        tot = 0.0
        t = z
        k = 1
        while True:
            term = t / k ** s
            tot += term
            if abs(term) < 1e-17 * max(abs(tot), 1e-300) or k > 5000:
                break
            k += 1
            t *= z
        return tot
    if (z > 1.0 or z < -1.0) and s != int(s) and abs(z) < 150.0:
        # Li_s(z) = Γ(1-s) (-μ)^{s-1} + Σ ζ(s-k) μ^k / k! (μ = ln z、|μ| < 2π)。z > 1 は μ = ln z で (-μ) = μ e^{iπ} (gs の符号)、z < -1 は μ = ln|z| + iπ
        if z > 1.0:
            mu = complex(math.log(z), 0.0)
            neg = cmath.exp(complex(math.log(math.log(z)), math.pi))            # (-μ) = μ e^{iπ} (gs の符号)
        else:
            mu = complex(math.log(-z), math.pi)
            neg = -mu
        tot = math.gamma(1.0 - s) * neg ** (s - 1.0)
        for k in range(0, 60):
            tot += _zeta_real(s - k) * mu ** k / math.factorial(k)
        return tot if abs(tot.imag) > 1e-12 * max(abs(tot.real), 1e-300) else tot.real
    if z > 1.0:
        if s == int(s) and s >= 2:
            n = int(s)
            # z > 1 は切断線の下側 (gs: Li2(2) = 2.4674 - 2.1776 i)。Li_n(z) + (-1)^n Li_n(1/z) = -(2πi)^n/n! B_n(1/2 + ln(-z)/(2πi)) の ln(-z) に ln z + iπ を入れる (符号は gs に合わせた)
            lz = complex(math.log(z), math.pi)
            val = -((2j * math.pi) ** n) / math.factorial(n) * sum(
                math.comb(n, k) * float(_bernoulli(k)) * (0.5 + lz / (2j * math.pi)) ** (n - k) for k in range(n + 1))
            return val - (-1) ** n * _polylog(s, 1.0 / z)
        return NAN
    if z < -1.0:
        if s == int(s) and s >= 1:
            n = int(s)
            # Li_n(z) + (-1)^n Li_n(1/z) = -(2πi)^n/n! B_n(1/2 + ln(z)/(2πi))
            lz = cmath.log(complex(-z))                 # ln(-z) (z<-1 なので実数)
            val = -((2j * math.pi) ** n) / math.factorial(n) * sum(
                math.comb(n, k) * float(_bernoulli(k)) * (0.5 + lz / (2j * math.pi)) ** (n - k) for k in range(n + 1))
            return _real_part(val - (-1) ** n * _polylog(s, 1.0 / z))
        return NAN
    # 0.8 < |z| < 1: μ = ln z の展開
    mu = cmath.log(complex(z))
    if s == int(s) and s >= 1:
        n = int(s)
        hn = sum(1.0 / k for k in range(1, n))
        tot = mu ** (n - 1) / math.factorial(n - 1) * (hn - cmath.log(-mu))
        for k in range(0, 40):
            if k == n - 1:
                continue
            tot += _zeta_real(n - k) * mu ** k / math.factorial(k)
        return _real_part(tot)
    tot = math.gamma(1.0 - s) * (-mu) ** (s - 1.0)
    for k in range(0, 40):
        tot += _zeta_real(s - k) * mu ** k / math.factorial(k)
    return _real_part(tot)


def _real_part(v):
    return v.real if isinstance(v, complex) else v


def _lerch(z, s, a):
    if abs(z) < 0.99:
        tot = 0.0
        t = 1.0
        for k in range(20000):
            term = t / (k + a) ** s
            tot += term
            if abs(term) < 1e-17 * abs(tot):
                break
            t *= z
        return tot
    if z == 1.0:
        return _hurwitz(s, a)
    if z == -1.0:
        return 2.0 ** (-s) * (_hurwitz(s, a / 2.0) - _hurwitz(s, (a + 1.0) / 2.0))
    return NAN


def _pfq(a, b, z, regularized=False):
    if isinstance(z, complex):
        raise SADError('HypergeometricPFQ: complex argument is not supported')
    tot = 1.0
    t = 1.0
    for k in range(0, 5000):
        num = 1.0
        for ai in a:
            num *= (ai + k)
        den = float(k + 1)
        for bi in b:
            den *= (bi + k)
        if den == 0.0:
            return NAN
        t *= num / den * z
        tot += t
        if t == 0.0 or abs(t) < 1e-17 * abs(tot):
            break
        if abs(t) > 1e300:
            return NAN
    if len(a) > len(b) + 1 or (len(a) == len(b) + 1 and abs(z) >= 1.0):
        return NAN
    if regularized:
        g = 1.0
        for bi in b:
            g *= _gamma(bi)
        return tot / g
    return tot


def _pochhammer(a, n):
    if n == int(n) and n >= 0:
        r = 1.0
        for k in range(int(n)):
            r *= (a + k)
        return r
    return _gamma(a + n) / _gamma(a)


def _hermite(n, z):
    if n == int(n) and n >= 0:
        return float(sp.eval_hermite(int(n), z))
    z2 = z * z
    sq = abs(z)
    g = _gamma((1.0 - n) / 2.0)
    return (2.0 ** n / (sq * g)) * (math.sqrt(math.pi) * (sq - z) * float(sp.hyp1f1(-n / 2.0, 0.5, z2))
                                   + z * g * float(sp.hyperu(-n / 2.0, 0.5, z2)))


def _legendre_p(n, m, z):
    """D.n: ((1+z)/(1-z))^(m/2) 2F1reg[-n, n+1, 1-m, (1-z)/2]。整数 m は lpmv (Condon-Shortley の位相を除く)"""
    if m == int(m) and n == int(n) and n >= 0:
        m = int(m)
        n = int(n)
        if m >= 0:
            return (-1) ** m * float(sp.lpmv(m, n, z)) if m <= n else 0.0
    if m == 0:
        return float(sp.hyp2f1(-n, n + 1, 1.0, (1.0 - z) / 2.0))
    g = _gamma(1.0 - m)
    if g in (INF, NAN) or g != g:
        return NAN
    return ((1.0 + z) / (1.0 - z)) ** (m / 2.0) * float(sp.hyp2f1(-n, n + 1, 1.0 - m, (1.0 - z) / 2.0)) / g


def _besseljzero(n, k):
    from scipy.optimize import brentq
    k = int(k)
    if k == 0:
        return NAN if n == 0 else 0.0
    if n == int(n) and n >= 0:
        return float(sp.jn_zeros(int(n), k)[-1])
    x = (k + n / 2.0 - 0.25) * math.pi
    x0 = x - (2 * n + 1) * (2 * n - 1) / x * (0.125 + (28 * n * n - 31) / 384.0 / x ** 2)
    lo, hi = x0 - 1.0, x0 + 1.0
    f = lambda t: float(sp.jv(n, t))
    for _ in range(40):
        if f(lo) * f(hi) < 0:
            break
        lo -= 0.5
        hi += 0.5
    return float(brentq(f, lo, hi, xtol=1e-14))


def inverseerf(x):
    """gamma.f inverseerf の直訳。K64 は Newton の係数に sqrt(pi/2) を使い (正しくは sqrt(pi)/2)、反復を df^2 <= no f1^2 eps で止めるので、
    真値から相対 1e-8 ほどずれる (InverseErf[.3] = .272462711811; 真値は .272462714727)。Fit.n の FitGaussian がこれを使う"""
    nolog = 8
    eps = 5e-17
    c = math.sqrt(math.pi / 2.0)
    if x != x:
        return NAN
    if x < 0.0:
        return -inverseerf(-x)
    if x == 1.0:
        return INF
    if x > 1.0:
        return NAN
    if x > 0.8:
        f1 = math.log(math.sqrt(2.0 / math.pi) / (1.0 - x))
        f1 = math.sqrt(f1 - 0.5 * math.log(2.0 * f1))
        no = nolog * 2 + 10
    else:
        x2 = math.pi * x ** 2
        f1 = 0.5 * math.sqrt(math.pi) * x * (1.0 + x2 * (1.0 / 12.0 + x2 * (7.0 / 480.0 + x2 * (127.0 / 40320.0 + x2 * 4369.0 / 5806080.0))))
        no = 11
    for _ in range(1000):
        df = (x - math.erf(f1)) * c * math.exp(f1 ** 2)
        f1 += df
        if df ** 2 <= no * f1 ** 2 * eps:
            break
        no += nolog
    return f1


def _ei(x):
    if x == 0.0:
        return -INF
    return float(sp.expi(x))


def _cosint(x):
    if x > 0.0:
        return float(sp.sici(x)[1])
    if x == 0.0:
        return -INF
    return complex(float(sp.sici(-x)[1]), math.pi)


def _xsin(x):
    if abs(x) < 0.1:
        t = x ** 3 / 6.0
        s = t
        k = 1
        while abs(t) > 1e-18 * abs(s):
            t *= -x * x / ((2 * k + 2) * (2 * k + 3))
            s += t
            k += 1
        return s
    return x - math.sin(x)


def _fourier(x, inverse=False):
    a = np.array(x, dtype=complex)
    n = len(a)
    r = (np.fft.fft(a) if inverse else np.conj(np.fft.fft(np.conj(a)))) / math.sqrt(n)
    return [_real(complex(v)) for v in r]


def install(I):
    B = I.builtins
    if sp is None:
        return

    def reg(name, f, nargs=None):
        def w(a, f=f, name=name):
            if nargs is not None:
                if isinstance(nargs, tuple):
                    if len(a) not in nargs:
                        raise SADError(f'General::narg: Number of arguments is expected to be {" or ".join(map(str, nargs))} in {name}')
                elif len(a) != nargs:
                    raise SADError(f'General::narg: Number of arguments is expected to be {nargs} in {name}')
            if not all(_numeric(x) for x in a):
                return Expr(S(name), list(a))          # 記号を含むときは評価しない (D.n の XSin[x] など)
            return _thread(f, a)
        B[name] = w

    # ---- 基本関数
    def arcsinh(x):
        return cmath.asinh(x) if isinstance(x, complex) else math.asinh(x)

    def arccosh(x):
        if isinstance(x, complex):
            return cmath.acosh(x)
        if x >= 1.0:
            return math.acosh(x)
        if x >= -1.0:
            return complex(0.0, -math.acos(x))      # K64: ArcCosh[.5] = -1.047 I
        return cmath.acosh(x)

    def arctanh(x):
        if isinstance(x, complex):
            return cmath.atanh(x)
        if abs(x) < 1.0:
            return math.atanh(x)
        if abs(x) == 1.0:
            return math.copysign(INF, x)
        return NAN                                    # K64: ArcTanh[2.] = NaN
    reg('ArcSinh', arcsinh, 1)
    reg('ArcCosh', arccosh, 1)
    reg('ArcTanh', arctanh, 1)
    reg('LogGamma', _loggamma, 1)
    reg('LogGamma1', lambda x: _loggamma(x + 1), 1)         # LogGamma1[x] = Log[x!] = LogGamma[x+1] (gs と 10.5 などで一致)
    reg('Gamma', lambda *a: _gamma(a[0]) if len(a) == 1 else float(sp.gammaincc(a[0], a[1]) * sp.gamma(a[0])))
    reg('Gamma0', lambda x: float(sp.exp1(x)), 1)
    reg('XSin', _xsin, 1)
    reg('GammaRegularized', lambda a, z: float(sp.gammaincc(a, z)), 2)
    reg('GammaRegularizedQ', lambda a, z: float(sp.gammaincc(a, z)), 2)
    reg('GammaRegularizedP', lambda a, z: float(sp.gammainc(a, z)), 2)
    reg('PolyGamma', lambda *a: _polygamma(0, a[0]) if len(a) == 1 else _polygamma(a[0], a[1]), (1, 2))
    reg('Beta', lambda a, b: float(sp.beta(a, b)), 2)
    reg('Zeta', _zeta_real, 1)
    reg('HurwitzZeta', _hurwitz, 2)
    reg('PolyLog', _polylog, 2)
    reg('LerchPhi', _lerch, 3)
    reg('HurwitzLerchPhi', _lerch, 3)
    reg('BernoulliB', lambda n, *x: float(_bernoulli(int(n))) if not x else _bernoulli_poly(int(n), x[0]), (1, 2))
    reg('ProductLog', lambda z: float(sp.lambertw(z).real) if z >= -1.0 / math.e else NAN, 1)
    reg('Pochhammer', _pochhammer, 2)
    reg('Hypergeometric0F1', lambda b, z: float(sp.hyp0f1(b, z)), 2)
    reg('Hypergeometric0F1Regularized', lambda b, z: float(sp.hyp0f1(b, z)) / _gamma(b), 2)
    reg('Hypergeometric1F1', lambda a, b, z: float(sp.hyp1f1(a, b, z)), 3)
    reg('Hypergeometric1F1Regularized', lambda a, b, z: float(sp.hyp1f1(a, b, z)) / _gamma(b), 3)
    reg('Hypergeometric2F1', lambda a, b, c, z: float(sp.hyp2f1(a, b, c, z)), 4)
    reg('Hypergeometric2F1Regularized', lambda a, b, c, z: float(sp.hyp2f1(a, b, c, z)) / _gamma(c), 4)
    reg('HypergeometricU', lambda a, b, z: float(sp.hyperu(a, b, z)), 3)
    B['HypergeometricPFQ'] = lambda a: _pfq(a[0], a[1], a[2])
    B['HypergeometricPFQRegularized'] = lambda a: _pfq(a[0], a[1], a[2], True)
    reg('BesselJ', lambda n, x: float(sp.jv(n, x)), 2)
    reg('BesselY', lambda n, x: float(sp.yv(n, x)), 2)
    reg('BesselI', lambda n, x: float(sp.iv(n, x)), 2)
    reg('BesselK', lambda n, x: float(sp.kv(n, x)), 2)
    B['Fourier'] = lambda a: _fourier(a[0])
    B['InverseFourier'] = lambda a: _fourier(a[0], True)

    # ---- D.n の特殊関数
    reg('InverseErf', inverseerf, 1)
    reg('InverseErfc', lambda z: inverseerf(1.0 - z), 1)
    reg('HarmonicNumber', lambda z, *r: float(EULER + sp.digamma(z + 1.0)) if not r else _zeta_real(r[0]) - _hurwitz(r[0], z + 1.0), (1, 2))
    reg('Factorial2', lambda n: (2.0 / math.pi) ** (math.sin(n * math.pi / 2.0) ** 2 / 2.0) * 2.0 ** (n / 2.0) * _gamma(1.0 + n / 2.0), 1)
    reg('ExpIntegralE', lambda n, z: float(sp.expn(int(n), z)) if n == int(n) and n >= 0 else z ** (n - 1) * math.exp(-z) * float(sp.hyperu(n, n, z)), 2)
    reg('ExpIntegralEi', _ei, 1)
    reg('LogIntegral', lambda z: _ei(math.log(z)) if z > 0 else NAN, 1)
    reg('SinIntegral', lambda z: float(sp.sici(z)[0]), 1)
    reg('CosIntegral', _cosint, 1)
    reg('SinhIntegral', lambda z: float(sp.shichi(z)[0]), 1)
    reg('CoshIntegral', lambda z: float(sp.shichi(z)[1]), 1)
    reg('WhittakerM', lambda n, m, z: z ** (m + 0.5) * float(sp.hyp1f1(m - n + 0.5, 2 * m + 1, z)) * math.exp(-z / 2.0), 3)
    reg('WhittakerW', lambda n, m, z: z ** (m + 0.5) * float(sp.hyperu(m - n + 0.5, 2 * m + 1, z)) * math.exp(-z / 2.0), 3)
    reg('HermiteH', _hermite, 2)
    reg('FresnelS', lambda z: float(sp.fresnel(z)[0]), 1)
    reg('FresnelC', lambda z: float(sp.fresnel(z)[1]), 1)
    reg('AiryAi', lambda z: float(sp.airy(z)[0]), 1)
    reg('AiryBi', lambda z: float(sp.airy(z)[2]), 1)
    reg('AiryAiPrime', lambda z: float(sp.airy(z)[1]), 1)
    reg('AiryBiPrime', lambda z: float(sp.airy(z)[3]), 1)
    B['LegendreP'] = lambda a: _thread(lambda *x: _legendre_p(x[0], 0.0, x[1]) if len(x) == 2 else _legendre_p(x[0], x[1], x[-1]), a) if len(a) in (2, 3, 4) else (_ for _ in ()).throw(SADError('General::narg: LegendreP'))
    reg('BesselJZero', _besseljzero, 2)

    def binom(n, k):
        if k == int(k) and k >= 0:
            r = 1.0
            for j in range(int(k)):
                r *= (n - j) / (j + 1.0)
            return r
        return _pochhammer(k + 1.0, n - k) / _gamma(n - k + 1.0)
    reg('Binomial', binom, 2)

    def multinomial(a):
        s = sum(a)
        return math.exp(math.lgamma(s + 1.0) - sum(math.lgamma(x + 1.0) for x in a))
    B['Multinomial'] = lambda a: float(round(multinomial(a))) if all(is_num(x) and x == int(x) for x in a) else multinomial(a)

    # ---- 数値・リスト・文字列の小物 (K64 の組み込み)
    def restrict(x, lo, hi):
        return min(max(x, lo), hi)
    reg('Restrict', restrict, 3)
    def _minmax(a):
        v = [y for x in a for y in _flat(x)]        # MinMax[ref, ll]: すべての引数の最小と最大 (Plots.n の ColumnPlot)
        return [min(v), max(v)]
    B['MinMax'] = _minmax
    B['NaNQ'] = lambda a: 1.0 if (is_num(a[0]) and not isinstance(a[0], complex) and a[0] != a[0]) else 0.0
    B['BitwiseAnd'] = lambda a: _thread(lambda x, y: float(int(x) & int(y)), a)
    B['BitwiseOr'] = lambda a: _thread(lambda x, y: float(int(x) | int(y)), a)
    B['BitwiseXor'] = lambda a: _thread(lambda x, y: float(int(x) ^ int(y)), a)
    B['DigitQ'] = lambda a: 1.0 if isinstance(a[0], str) and a[0] != '' and a[0].isdigit() else 0.0
    B['LetterQ'] = lambda a: 1.0 if isinstance(a[0], str) and a[0] != '' and a[0].isalpha() else 0.0
    B['RealListQ'] = lambda a: 1.0 if isinstance(a[0], list) and all(is_num(x) and not isinstance(x, complex) for x in a[0]) else 0.0
    B['SymbolName'] = lambda a: a[0].name if isinstance(a[0], Sym) else a[0]
    B['Complex'] = lambda a: complex(a[0], a[1])
    B['Duplicate'] = lambda a: _copy(a[0])
    B['FlattenAt'] = lambda a: _flattenat(a[0], int(a[1]))

    def inner(a):
        f, x, y, g = a
        r = I.apply_value(f, [x, y], f)
        if isinstance(r, list):
            return I.apply_value(g, list(r), g)
        return r
    B['Inner'] = inner
    reg('InversePower', lambda a, b: float(b) ** a, 2)
    B['Union1'] = lambda a: sorted(set(a[0]))
    B['RotateRight1'] = lambda a: a[0][-1:] + a[0][:-1]
    B['Range$'] = lambda a: I.builtins['Range'](a)


def _flat(x):
    out = []
    for t in (x if isinstance(x, list) else [x]):
        out.extend(_flat(t) if isinstance(t, list) else [t])
    return out


def _copy(x):
    return [_copy(t) for t in x] if isinstance(x, list) else x


def _flattenat(lst, n):
    k = n - 1 if n > 0 else len(lst) + n
    out = []
    for i, t in enumerate(lst):
        if i == k and isinstance(t, list):
            out.extend(t)
        else:
            out.append(t)
    return out
