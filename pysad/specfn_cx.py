"""複素引数の Zeta・HurwitzZeta・PolyLog・HypergeometricU (v2 の追加。K64 の gamma.f の複素版は誤る範囲があるので、真値を返す)。

mpmath が無い環境でも動くよう、級数・Euler–Maclaurin・解析接続の公式・積分表示だけで書いた。数値の根拠 (独立な参照) は
`pysad/tests/script_tests/test_v2_B2_specfn_cx.py` (別の方法: Borwein のアルゴリズム、関数等式、scipy.special.spence の複素版、積分表示、漸近展開) で確かめる。

- Zeta・HurwitzZeta: 複素の s と a (a は非正の整数でない)。Euler–Maclaurin (N 項 + 補正 16 項、N は |s| に応じて増やす)。Re s < 0 の Zeta は関数等式。
  使える範囲: |s| ≲ 500 (それ以上は SADError)。
- PolyLog: 実数の次数 s、複素数の z。主値 (切断線は z ∈ [1, ∞)。実軸上の z > 1 は specfn の実数版が下側の値を返す)。
  |z| ≤ 0.5 は冪級数、|ln z| ≤ 4 は μ = ln z の展開 (Γ(1-s)(-μ)^{s-1} + Σ ζ(s-k) μ^k/k!、整数 s は調和数の項つき)、それより大きい |z| は
  倍角の公式 Li_s(z) = 2^{s-1} (Li_s(√z) + Li_s(-√z)) を再帰。s が 0 以下の整数は有理関数の閉じた式。
- HypergeometricU: 実数または複素数の z (主値、負の実軸は上側の値 arg z = +π)、実数の a・b (複素の a・b は Re a > 0 だけ)。
  Re a > 0 は積分表示 U = 1/Γ(a) ∫ e^{-zt} t^{a-1} (1+t)^{b-a-1} dt を、被積分関数が減衰する向き (t = r e^{iφ}) に回した経路で数値積分する。
  a が 0 以下: a - b + 1 > 0 なら Kummer の変換 U(a,b,z) = z^{1-b} U(a-b+1, 2-b, z)、a が負の整数なら多項式、それ以外は a の漸化式。b が整数でもよい。
"""
from __future__ import annotations

import cmath
import math

import numpy as np

from .sadscript import SADError, INF, NAN

try:
    from scipy import special as sp
    from scipy import integrate
except ImportError:  # pragma: no cover
    sp = None
    integrate = None

_B2J = None


def _bern(n):
    from .specfn import _bernoulli
    return float(_bernoulli(n))


def _zeta_real(s):
    from .specfn import _zeta_real as z
    return z(s)


# ---------------------------------------------------------------- Hurwitz ζ (複素の s・a)
def hurwitz_c(s, a):
    """ζ(s, a) (s ≠ 1、a は非正の整数でない)。複素の s・a。Re a < 1/2 なら Σ_{k<m} (a+k)^-s の項を足して Re a ≥ 1/2 に移す"""
    s = complex(s)
    a = complex(a)
    if s == 1.0:
        return complex(INF, 0.0)
    if abs(s) > 500.0:
        raise SADError('HurwitzZeta: |s| > 500 is not supported for complex arguments')
    if a.imag == 0.0 and a.real <= 0.0 and a.real == math.floor(a.real):
        return complex(NAN, 0.0)
    head = 0j
    if a.real < 0.5:
        m = int(math.ceil(0.5 - a.real))
        for k in range(m):
            head += cmath.exp(-s * cmath.log(a + k))             # 主値 (a+k が負の実数なら ln = ln|.| + iπ。gs の実数版と同じ)
        a = a + m
    N = int(20 + abs(s) / 1.5 + abs(a.imag) / 2.0)
    M = 16
    re, im = [head.real], [head.imag]
    for k in range(N):
        t = cmath.exp(-s * cmath.log(a + k))
        re.append(t.real)
        im.append(t.imag)
    x = a + N
    lx = cmath.log(x)
    for t in (cmath.exp((1.0 - s) * lx) / (s - 1.0), 0.5 * cmath.exp(-s * lx)):
        re.append(t.real)
        im.append(t.imag)
    p = s
    xp = cmath.exp((-s - 1.0) * lx)
    for j in range(1, M + 1):
        t = _bern(2 * j) / math.factorial(2 * j) * p * xp
        re.append(t.real)
        im.append(t.imag)
        p *= (s + 2 * j - 1) * (s + 2 * j)
        xp /= x * x
    return complex(math.fsum(re), math.fsum(im))


def zeta_c(s):
    """ζ(s) (複素の s ≠ 1)。Re s ≥ 0 は Hurwitz (a = 1)、Re s < 0 は関数等式 ζ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s) ζ(1-s)"""
    s = complex(s)
    if s == 1.0:
        return complex(INF, 0.0)
    if s.real >= 0.0:
        return hurwitz_c(s, 1.0)
    return (cmath.exp(s * math.log(2.0) + (s - 1.0) * math.log(math.pi)) * cmath.sin(0.5 * math.pi * s)
            * complex(sp.gamma(1.0 - s)) * hurwitz_c(1.0 - s, 1.0))


# ---------------------------------------------------------------- PolyLog (実数の s、複素数の z)
def _li_negint(n, z):
    """Li_{-n}(z) (n ≥ 0) = Σ_{k=0}^{n} k! S(n+1, k+1) (z/(1-z))^{k+1} (S は第 2 種のスターリング数)"""
    # S(n+1, k+1) の表
    S = [[0] * (n + 2) for _ in range(n + 2)]
    S[0][0] = 1
    for i in range(1, n + 2):
        for j in range(1, i + 1):
            S[i][j] = j * S[i - 1][j] + S[i - 1][j - 1]
    w = z / (1.0 - z)
    tot = 0j
    wk = w
    for k in range(n + 1):
        tot += math.factorial(k) * S[n + 1][k + 1] * wk
        wk *= w
    return tot


def _li_series(s, z):
    tot = 0j
    t = z
    for k in range(1, 400):
        term = t / (float(k) ** s)
        tot += term
        if abs(term) < 1e-18 * max(abs(tot), 1e-300):
            break
        t *= z
    return tot


def _li_mu(s, z):
    """μ = ln z (|μ| ≤ 4) の展開。z ≠ 1"""
    mu = cmath.log(z)
    if s == int(s) and s >= 1:
        n = int(s)
        hn = sum(1.0 / k for k in range(1, n))
        tot = mu ** (n - 1) / math.factorial(n - 1) * (hn - cmath.log(-mu))
        for k in range(0, 110):
            if k == n - 1:
                continue
            tot += _zeta_real(n - k) * mu ** k / math.factorial(k)
        return tot
    tot = complex(sp.gamma(1.0 - s)) * cmath.exp((s - 1.0) * cmath.log(-mu))
    for k in range(0, 110):
        tot += _zeta_real(s - k) * mu ** k / math.factorial(k)
    return tot


def polylog_c(s, z, _depth=0):
    """Li_s(z) (実数の s、複素数の z、z ∉ [1, ∞))。主値"""
    z = complex(z)
    s = float(s)
    if z == 0.0:
        return 0j
    if z == 1.0:
        return complex(_zeta_real(s), 0.0)
    if s == int(s) and s <= 0:
        return _li_negint(int(-s), z)
    if s == 1.0:
        return -cmath.log(1.0 - z)
    if abs(z) <= 0.5:
        return _li_series(s, z)
    mu = cmath.log(z)
    if abs(mu) <= 4.0:
        return _li_mu(s, z)
    if _depth > 40:
        raise SADError('PolyLog: argument is too large')
    w = cmath.sqrt(z)
    return 2.0 ** (s - 1.0) * (polylog_c(s, w, _depth + 1) + polylog_c(s, -w, _depth + 1))


# ---------------------------------------------------------------- HypergeometricU
def _quad_c(f, a, b):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        re = integrate.quad(lambda t: f(t).real, a, b, epsabs=0, epsrel=1e-13, limit=400)[0]
        im = integrate.quad(lambda t: f(t).imag, a, b, epsabs=0, epsrel=1e-13, limit=400)[0]
    return complex(re, im)


def _u_integral(a, b, z):
    """Re a > 0。U = 1/Γ(a) ∫_0^∞ e^{-zt} t^{a-1} (1+t)^{b-a-1} dt を t = r e^{iφ} の経路で。φ = -arg z (負の実軸の近くは 0.4 だけ手前)"""
    z = complex(z)
    th = cmath.phase(z)
    if th > math.pi - 0.4:
        phi = -th + 0.4
    elif th < -math.pi + 0.4:
        phi = -th - 0.4
    else:
        phi = -th
    e = cmath.exp(1j * phi)
    ra = a.real if isinstance(a, complex) else a
    p = 1.0 / ra if ra < 1.0 else 1.0               # t = u^p で t^{a-1} の端点の特異性を除く (a の実部 < 1)
    absz = abs(z)
    rmax = 745.0 / (absz * max(math.cos(phi + th), 0.05)) + 1.0     # e^{-|z| r cos} が 1e-323 を下回る所

    def f(u):
        if u <= 0.0:
            return 0j
        r = u ** p
        t = r * e
        dr = p * u ** (p - 1.0)
        return cmath.exp(-z * t + (a - 1.0) * cmath.log(t) + (b - a - 1.0) * cmath.log(1.0 + t)) * e * dr
    umax = rmax ** (1.0 / p)
    # 区間: 0 〜 1/|z| 付近 (被積分関数の山) 〜 十分遠く
    pts = sorted({0.0, min(umax, (0.3 / absz) ** (1.0 / p)), min(umax, (3.0 / absz) ** (1.0 / p)), min(umax, (30.0 / absz) ** (1.0 / p)), umax})
    tot = 0j
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi > lo:
            tot += _quad_c(f, lo, hi)
    return tot / complex(sp.gamma(a))


def _u_poly(n, b, z):
    """U(-n, b, z) = (-1)^n Σ_k (-n)_k (b+k)_{n-k} z^k / k!"""
    tot = 0j
    for k in range(n + 1):
        c = 1.0
        for j in range(k):
            c *= (-n + j)
        c /= math.factorial(k)
        for j in range(n - k):
            c *= (b + k + j)
        tot += c * z ** k
    return (-1) ** n * tot


def hyperu_c(a, b, z):
    """U(a, b, z)。z は複素数 (実数の負の値は上側 arg z = +π)。a・b は実数 (複素の a・b は Re a > 0 だけ)"""
    z = complex(z)
    if z == 0.0:
        raise SADError('HypergeometricU: z = 0 is not supported')
    if isinstance(a, complex) or isinstance(b, complex):
        a = complex(a)
        b = complex(b)
        if a.real <= 0.0:
            raise SADError('HypergeometricU: complex a or b needs Re a > 0')
        return _u_integral(a, b, z)
    a = float(a)
    b = float(b)
    if z.imag == 0.0 and z.real < 0.0:
        z = complex(z.real, 0.0)                # 主値の切断線の上側 (cmath.phase = +π)
    if a > 0.0:
        return _u_integral(a, b, z)
    if a == int(a):
        return _u_poly(int(-a), b, z)
    if a - b + 1.0 > 0.0:
        return cmath.exp((1.0 - b) * cmath.log(z)) * hyperu_c(a - b + 1.0, 2.0 - b, z)
    # a の漸化式: U(a-1) = (2a - b + z) U(a) - a (a - b + 1) U(a+1)。a0 = a + m > 0 から下ろす
    m = int(math.ceil(-a)) + 1
    a0 = a + m
    um1 = hyperu_c(a0 + 1.0, b, z)
    u0 = hyperu_c(a0, b, z)
    for j in range(m):
        ak = a0 - j                              # U(ak - 1) を作る
        un = (2.0 * ak - b + z) * u0 - ak * (ak - b + 1.0) * um1
        um1, u0 = u0, un
    return u0

