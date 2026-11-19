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


def gamma_p_c(a, z):
    """P(a, z) = γ(a, z)/Γ(a) = z^a e^{-z}/Γ(a+1) 1F1(1; a+1; z) (複素の z、実数の a > 0)"""
    z = _c(z)
    a = float(a)
    if z == 0:
        return 0j
    return cmath.exp(a * cmath.log(z) - z - float(sp.gammaln(a + 1.0))) * complex(sp.hyp1f1(1.0, a + 1.0, z))


def gamma_reg_c(a, z):
    """Q(a, z) = Γ(a, z)/Γ(a)"""
    return 1.0 - gamma_p_c(a, z)


def gamma_upper_c(a, z):
    return float(sp.gamma(a)) * gamma_reg_c(a, z)


def legendre_c(nu, z):
    """P_ν(z) = 2F1(-ν, ν+1; 1; (1-z)/2) (複素の z。|1-z|/2 < 1 の外は scipy の解析接続)"""
    return complex(sp.hyp2f1(-float(nu), float(nu) + 1.0, 1.0, (1.0 - _c(z)) / 2.0))


def whittaker_m_c(k, m, z):
    z = _c(z)
    return cmath.exp(-z / 2.0 + (m + 0.5) * cmath.log(z)) * complex(sp.hyp1f1(m - k + 0.5, 2.0 * m + 1.0, z))


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
