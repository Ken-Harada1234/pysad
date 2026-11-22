"""複素引数の特殊関数の追加 (v2next P5 項目 3)。両モード (K64 の動作の追加。gs が正しく返すものは gs と一致し、gs が誤るものは真値)。

追加するもの (どれも scipy.special の複素対応と、級数・漸化式・既存の specfn_cx で書く。mpmath は使わない):
  ExpIntegralE[n,z] (整数 n >= 0)、ExpIntegralEi[z]、PolyGamma[z]・PolyGamma[n,z]、HarmonicNumber[z]、ProductLog[z]、SinIntegral[z]・CosIntegral[z]、
  FresnelS[z]・FresnelC[z]、Beta[a,b]、Factorial[z]、GammaRegularized[a,z]、Gamma[a,z] (z が複素数)、LegendreP[nu,z]、WhittakerM・WhittakerW[k,m,z]、BesselJ[nu,x] (複素の次数)。
従来の pysad はこれらに複素数を渡すと虚部を捨てた実数を返したり (ExpIntegralEi・SinIntegral・FresnelS/C・HarmonicNumber・LegendreP)、何も返さなかったりした
(ExpIntegralE・PolyGamma・WhittakerM/W・ProductLog・CosIntegral・GammaRegularized・Beta・Factorial・BesselJ の複素の次数)。
独立な根拠は pysad/tests/test_v2_P5_specfn_cx.py (微分方程式の残差、漸化式、恒等式、複素の線積分)。使える範囲: |z| <~ 30 (WhittakerW・Legendre は |z| <~ 20)。
"""
from __future__ import annotations

import cmath
import math

import numpy as np

from .sadscript import SADError

try:
    from scipy import special as sp
except ImportError:  # pragma: no cover
    sp = None

EULER = 0.5772156649015329


def _c(z):
    return complex(z)


def expn_c(n, z):
    """E_n(z) = ∫_1^∞ e^{-zt} t^{-n} dt (整数 n >= 0)。E_0 = e^{-z}/z、E_1 = exp1(z)、E_{k+1} = (e^{-z} - z E_k)/k (|z| が小さい / n が大きい時は下向きに安定でないので、
    n > |z| のときは級数 E_n = (-z)^{n-1}/(n-1)! (-ln z + ψ(n)) - Σ_{k != n-1} (-z)^k / ((k - n + 1) k!) を使う)"""
    z = _c(z)
    n = int(round(n.real if isinstance(n, complex) else n))
    if n < 0:
        raise SADError('General::wrongtype: ExpIntegralE order must be a non-negative integer for complex z')
    if n == 0:
        return cmath.exp(-z) / z
    if n > abs(z) + 1:
        tot = 0j
        term = 1.0 + 0j
        for k in range(0, 80):
            if k > 0:
                term *= -z / k
            if k != n - 1:
                tot -= term / (k - n + 1)
        tot += ((-z) ** (n - 1) / math.factorial(n - 1)) * (-cmath.log(z) + float(sp.digamma(n)))
        return tot
    e = complex(sp.exp1(z))
    ez = cmath.exp(-z)
    for k in range(1, n):
        e = (ez - z * e) / k
    return e


def ei_c(z):
    return complex(sp.expi(_c(z)))


def digamma_c(z):
    return complex(sp.digamma(_c(z)))


def polygamma_c(n, z):
    n = int(round(n))
    if n == 0:
        return digamma_c(z)
    from . import specfn_cx
    return (-1) ** (n + 1) * math.factorial(n) * complex(specfn_cx.hurwitz_c(float(n + 1), _c(z)))


def harmonic_c(z):
    return digamma_c(_c(z) + 1.0) + EULER


def lambertw_c(z):
    return complex(sp.lambertw(_c(z)))


def sici_c(z):
    """(Si(z), Ci(z))。scipy の sici は複素数に対応する"""
    si, ci = sp.sici(_c(z))
    return complex(si), complex(ci)


def fresnel_c(z):
    s, c = sp.fresnel(_c(z))
    return complex(s), complex(c)


def beta_c(a, b):
    a, b = _c(a), _c(b)
    return complex(sp.gamma(a)) * complex(sp.gamma(b)) / complex(sp.gamma(a + b))


def factorial_c(z):
    return complex(sp.gamma(_c(z) + 1.0))


def hyp1f1_c(a, b, z):
    """1F1(a; b; z) (実数の a・b、複素の z)。|z| < 6 は scipy。それより大きい |z| の scipy は複素数で不正確 (|z| = 25〜30 で相対 1e-4〜1) なので、
    べき級数を十進の多倍長 (桁数 = 25 + 0.45|z|。項の最大値 e^{|z|} の桁落ちを吸収) で足す (v2.4.0 E。独立な基準は mpmath との比較)"""
    z = _c(z)
    a = float(a)
    b = float(b)
    if abs(z) < 6.0:
        return complex(sp.hyp1f1(a, b, z))
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 30 + int(0.45 * abs(z))
        da, db = Decimal(a), Decimal(b)
        zr, zi = Decimal(z.real), Decimal(z.imag)
        tr, ti = Decimal(1), Decimal(0)
        sr, si = Decimal(1), Decimal(0)
        eps = Decimal(10) ** (-25)
        for k in range(0, 20000):
            den = (db + k) * (k + 1)
            if den == 0:
                raise SADError('General::wrongval: Hypergeometric1F1: b is a non-positive integer')
            f = (da + k) / den
            tr, ti = (tr * zr - ti * zi) * f, (tr * zi + ti * zr) * f
            sr += tr
            si += ti
            if (tr == 0 and ti == 0) or (k > abs(z) and abs(tr) + abs(ti) < eps * (abs(sr) + abs(si))):
                break
        return complex(float(sr), float(si))


def gamma_p_c(a, z):
    """P(a, z) = γ(a, z)/Γ(a) = z^a e^{-z}/Γ(a+1) 1F1(1; a+1; z) (複素の z、実数の a > 0)"""
    z = _c(z)
    a = float(a)
    if z == 0:
        return 0j
    return cmath.exp(a * cmath.log(z) - z - float(sp.gammaln(a + 1.0))) * hyp1f1_c(1.0, a + 1.0, z)


def _gamma_q_cf(a, z):
    """Q(a, z) の連分数 (Lentz。Re z > 0 かつ |z| ≥ 10。1 - P の桁落ちを避ける)"""
    tiny = 1e-300
    b = z + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        dl = d * c
        h *= dl
        if abs(dl - 1.0) < 1e-16:
            break
    return cmath.exp(-z + a * cmath.log(z) - float(sp.gammaln(a))) * h


def gamma_reg_c(a, z):
    """Q(a, z) = Γ(a, z)/Γ(a)"""
    zc = _c(z)
    if zc.real > 0.0 and abs(zc) >= 10.0 and float(a) > 0.0:
        return _gamma_q_cf(float(a), zc)
    return 1.0 - gamma_p_c(a, z)


def gamma_upper_c(a, z):
    return float(sp.gamma(a)) * gamma_reg_c(a, z)


def legendre_c(nu, z):
    """P_ν(z) = 2F1(-ν, ν+1; 1; (1-z)/2) (複素の z。|1-z|/2 < 1 の外は scipy の解析接続)"""
    return complex(sp.hyp2f1(-float(nu), float(nu) + 1.0, 1.0, (1.0 - _c(z)) / 2.0))


def whittaker_m_c(k, m, z):
    z = _c(z)
    return cmath.exp(-z / 2.0 + (m + 0.5) * cmath.log(z)) * hyp1f1_c(m - k + 0.5, 2.0 * m + 1.0, z)


def whittaker_w_c(k, m, z):
    from . import specfn_cx
    z = _c(z)
    return cmath.exp(-z / 2.0 + (m + 0.5) * cmath.log(z)) * complex(specfn_cx.hyperu_c(m - k + 0.5, 2.0 * m + 1.0, z))


def besselj_cord(nu, x):
    """J_ν(x) の複素の次数 ν (実数 x > 0 か複素数 x): (x/2)^ν Σ (-x²/4)^k / (k! Γ(ν+k+1))。|x| <~ 20 で使える (それ以上は SADError)"""
    nu, x = _c(nu), _c(x)
    if abs(x) > 25.0:
        raise SADError('General::wrongval: BesselJ with a complex order is limited to |x| <= 25')
    q = -x * x / 4.0
    term = 1.0 / complex(sp.gamma(nu + 1.0))
    tot = term
    for k in range(1, 200):
        term *= q / (k * (nu + k))
        tot += term
        if abs(term) < 1e-17 * abs(tot):
            break
    return cmath.exp(nu * cmath.log(x / 2.0)) * tot


# ---------------------------------------------------------------- v2.4.0 E: 複素の次数 (PolyLog・ExpIntegralE)
def polylog_cs(s, z):
    """Li_s(z) (複素の次数 s、複素の z)。|z| ≤ 0.5 は級数、それ以外で |ln z| ≤ 5 は μ = ln z の展開
    Li_s(z) = Γ(1-s)(-μ)^{s-1} + Σ ζ(s-k) μ^k/k! (s が正の整数でないとき)。|Log z| > 5 は反転公式 (|z| ≥ 2 なので |1/z| ≤ 0.5 の級数)。独立な基準は mpmath"""
    from . import specfn_cx
    s = _c(s)
    z = _c(z)
    if z == 0:
        return 0j
    if s.imag == 0.0:
        return specfn_cx.polylog_c(s.real, z)
    if abs(z) <= 0.5:
        tot = 0j
        t = z
        for k in range(1, 400):
            term = t * cmath.exp(-s * math.log(k))
            tot += term
            if abs(term) < 1e-18 * max(abs(tot), 1e-300):
                break
            t *= z
        return tot
    mu = cmath.log(z)
    if z == 1.0:
        return specfn_cx.zeta_c(s)
    if abs(mu) > 5.0:
        # |Log z| > 5 (|z| > 148): 反転公式 Li_s(z) = -e^{iπs} Li_s(1/z) + (2π)^s e^{iπs/2}/Γ(s) ζ(1-s, 1/2 + ln(-z)/(2πi)) (|1/z| < 0.007 は級数)。
        # 実軸上の z > 1 は ln(-z) = ln z + iπ (mpmath の値と同じ側。独立な基準は mpmath)
        lnmz = cmath.log(-z)
        if z.imag == 0.0 and z.real > 0.0:
            lnmz = complex(math.log(z.real), math.pi)
        a = 0.5 + lnmz / (2j * math.pi)
        head = cmath.exp(s * math.log(2.0 * math.pi) + 0.5j * math.pi * s)
        term = head / complex(sp.gamma(s)) * specfn_cx.hurwitz_c(1.0 - s, a)
        return -cmath.exp(1j * math.pi * s) * polylog_cs(s, 1.0 / z) + term
    lm = complex(math.log(mu.real), math.pi) if (mu.imag == 0.0 and mu.real > 0.0) else cmath.log(-mu)       # z > 1 は (-μ) = μ e^{iπ} (-0.0 の虚部で符号が逆にならないように)
    tot = complex(sp.gamma(1.0 - s)) * cmath.exp((s - 1.0) * lm)
    for k in range(0, 120):
        tot += specfn_cx.zeta_c(s - k) * mu ** k / math.factorial(k)
    return tot


def _gamma_upper_ca(a, z):
    """Γ(a, z) (複素の a・z、a が非正の整数でない)。|z| ≤ 8 は Γ(a) - γ(a, z) の級数、それ以上は連分数 (|arg z| < 2.6)"""
    a, z = _c(a), _c(z)
    if abs(z) > 8.0 and abs(cmath.phase(z)) >= 2.6:
        if abs(z) <= 60.0:
            # Re z < 0 の近く: γ(a,z) = z^a Σ (-z)^k/(k!(a+k)) (項が同符号でそろい桁落ちが小さい)。Γ(a,z) = Γ(a) - γ(a,z)
            t = 1.0
            tot = 1.0 / a
            for k in range(1, 800):
                t *= -z / k
                u = t / (a + k)
                tot += u
                if abs(u) < 1e-18 * abs(tot) and k > abs(z):
                    break
            return complex(sp.gamma(a)) - cmath.exp(a * cmath.log(z)) * tot
        # |z| > 60: 漸近展開 Γ(a,z) ~ z^{a-1} e^{-z} Σ (a-1)(a-2)...(a-k)/z^k (|arg z| < 3π/2。最小項 ~ e^{-|z|})
        t = 1.0 + 0j
        tot = 1.0 + 0j
        best = 1e300
        for k in range(1, 400):
            t *= (a - k) / z
            if abs(t) > best:
                break
            best = abs(t)
            tot += t
            if abs(t) < 1e-18:
                break
        return cmath.exp((a - 1.0) * cmath.log(z) - z) * tot
    if abs(z) <= 8.0:
        # γ(a,z) = z^a e^{-z} Σ z^k / (a)_{k+1}
        t = 1.0 / a
        tot = t
        for k in range(1, 500):
            t *= z / (a + k)
            tot += t
            if abs(t) < 1e-18 * abs(tot):
                break
        g = cmath.exp(a * cmath.log(z) - z) * tot
        return complex(sp.gamma(a)) - g
    tiny = 1e-300
    b = z + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 3000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        dl = d * c
        h *= dl
        if abs(dl - 1.0) < 1e-16:
            break
    return cmath.exp(-z + a * cmath.log(z)) * h


def expn_cn(n, z):
    """E_ν(z) = z^{ν-1} Γ(1-ν, z) (複素の次数 ν・複素の z。ν が正の整数なら expn_c)"""
    n, z = _c(n), _c(z)
    if n.imag == 0.0 and n.real == int(n.real) and n.real >= 0:
        return expn_c(int(n.real), z)
    return cmath.exp((n - 1.0) * cmath.log(z)) * _gamma_upper_ca(1.0 - n, z)
