"""v2 の特殊関数 (担当 B: B43〜B47、B53、B54)。

根拠は gs (K64 は誤値・NaN・INF・無応答) ではなく、scipy.special と独立な計算:
閉じた式、級数、AGM、積分表示の数値積分 (台形則・scipy.integrate.quad)、恒等式 (Wronskian、erf(InverseErf(x)) = x)。
mpmath はこの環境に無いので使っていない。"""
import io
import math
import re

import numpy as np
import pytest
from scipy import integrate

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
EULER = 0.5772156649015328606


def _vals(exprs, mode='v2'):
    """式のリストを pysad で評価して数値 (複素は (re, im)) のリストを返す"""
    from pysad import compat
    from pysad.ffscmd import run_sad_string
    old = compat._mode
    compat.set_mode(mode)
    try:
        src = 'FFS;\n$FORM="25.17";\nPageWidth=1999;\n' + ''.join('Print["R%d= ",%s];\n' % (i, e) for i, e in enumerate(exprs)) + 'Exit;\n'
        buf = io.StringIO()
        run_sad_string(src, out=buf)
    finally:
        compat.set_mode(old)
    got = {}
    for l in buf.getvalue().splitlines():
        m = re.match(r'R(\d+)= (.*)', l)
        if m:
            t = m.group(2).replace('-INF', '-1e999').replace('INF', '1e999')
            if 'NaN' in t:
                got[int(m.group(1))] = [float('nan')]
            else:
                got[int(m.group(1))] = [float(x) for x in re.findall(_NUM, t)]
    return [got.get(i) for i in range(len(exprs))]


def _close(a, b, rel=1e-12, ab=0.0):
    return abs(a - b) <= rel * abs(b) + ab


# ---------------------------------------------------------------- 独立な計算
def _agm_K(m):
    """第 1 種完全楕円積分 K(m) = π / (2 AGM(1, sqrt(1-m)))"""
    a, b = 1.0, math.sqrt(1.0 - m)
    for _ in range(60):
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return math.pi / (2.0 * a)


def _jn_int(n, x, N=4000):
    """整数次の J_n(x) = (1/π) ∫_0^π cos(nθ - x sinθ) dθ (周期関数の台形則は指数的に収束)"""
    th = (np.arange(N) + 0.5) * math.pi / N
    return float(np.mean(np.cos(n * th - x * np.sin(th))))


def _knu_int(nu, x):
    """K_ν(x) = ∫_0^∞ exp(-x cosh t) cosh(νt) dt"""
    return integrate.quad(lambda t: math.exp(-x * math.cosh(t)) * math.cosh(nu * t), 0.0, 20.0, epsabs=0, epsrel=1e-13, limit=400)[0]


def _y0_series(x):
    """Y_0(x) = (2/π)(ln(x/2)+γ) J_0(x) + (2/π) Σ_{k≥1} (-1)^{k+1} H_k (x²/4)^k/(k!)²"""
    j0 = sum((-1) ** k * (x * x / 4.0) ** k / math.factorial(k) ** 2 for k in range(60))
    s = 0.0
    h = 0.0
    for k in range(1, 60):
        h += 1.0 / k
        s += (-1) ** (k + 1) * h * (x * x / 4.0) ** k / math.factorial(k) ** 2
    return 2.0 / math.pi * (math.log(x / 2.0) + EULER) * j0 + 2.0 / math.pi * s


def _ei_series(x):
    """Ei(x) = γ + ln|x| + Σ x^k/(k k!)"""
    return EULER + math.log(abs(x)) + sum(x ** k / (k * math.factorial(k)) for k in range(1, 120))


def _fresnel_quad(x):
    """S(x) = ∫_0^x sin(πt²/2) dt、C(x) = ∫_0^x cos(πt²/2) dt (Simpson、刻み 2.5e-5)"""
    n = 2 * int(x / 5e-5)
    t = np.linspace(0.0, x, n + 1)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    h = x / n
    return (float(np.sum(w * np.sin(0.5 * math.pi * t * t)) * h / 3.0), float(np.sum(w * np.cos(0.5 * math.pi * t * t)) * h / 3.0))


# ---------------------------------------------------------------- 試験
def test_b43_hyp2f1_besselj_precision():
    v = _vals(['Hypergeometric2F1[.5,.5,1,.99]', 'BesselJ[.5,30.]', 'Hypergeometric2F1[.5,.5,1,.999999]'])
    assert _close(v[0][0], 2.0 / math.pi * _agm_K(0.99), 1e-13)                 # K64: 2.348414 (相対 1.8e-3 の誤差)
    assert _close(v[1][0], math.sqrt(2.0 / (math.pi * 30.0)) * math.sin(30.0), 1e-13)      # K64: -.143924709693
    assert _close(v[2][0], 2.0 / math.pi * _agm_K(0.999999), 1e-12)


def test_b44_b54_expintegralei():
    xs = [1.5, -1.5, 0.3, -0.1, 10.0, 1e-3]                      # 級数は x ≲ -5 で桁落ちするので、負の大きい側は下の積分で確かめる
    v = _vals(['ExpIntegralEi[%r]' % x for x in xs])
    for x, r in zip(xs, v):
        assert len(r) == 1                                                       # 実数 (K64 は複素数を返す)
        assert _close(r[0], _ei_series(x), 2e-12), (x, r)
    # 大きい引数: Ei(x) ~ e^x/x (1 + 1/x + 2/x² + 6/x³ + …)
    big = _vals(['ExpIntegralEi[50.]', 'ExpIntegralEi[700.]'])
    for x, r in zip((50.0, 700.0), big):
        asym = math.exp(x) / x * sum(math.factorial(k) / x ** k for k in range(12))
        assert _close(r[0], asym, 1e-10), (x, r)
    # 小さい負の引数: Ei(-x) = -E1(x) = -∫_x^∞ e^-t/t dt
    e1 = integrate.quad(lambda t: math.exp(-t) / t, 0.1, 60.0, epsabs=0, epsrel=1e-13)[0]
    assert _close(_vals(['ExpIntegralEi[-.1]'])[0][0], -e1, 1e-11)
    e1 = integrate.quad(lambda t: math.exp(-t) / t, 10.0, 80.0, epsabs=0, epsrel=1e-13)[0]
    assert _close(_vals(['ExpIntegralEi[-10.]'])[0][0], -e1, 1e-10)
    # LogIntegral[z] = Ei(ln z)
    assert _close(_vals(['LogIntegral[2.]'])[0][0], _ei_series(math.log(2.0)), 1e-12)


def test_b45_hermite_legendre_binomial():
    v = _vals(['HermiteH[3,.7]', 'HermiteH[4,-1.3]', 'HermiteH[0,2.]', 'LegendreP[3,1,.4]', 'LegendreP[2,2,.3]', 'LegendreP[3,.4]',
               'LegendreP[4,3,-.3]', 'Binomial[-3,2]', 'Binomial[5,2]', 'Binomial[.5,3]'])
    x = 0.7
    assert _close(v[0][0], 8 * x ** 3 - 12 * x)                                  # H3 = 8x³ - 12x (K64: NaN)
    x = -1.3
    assert _close(v[1][0], 16 * x ** 4 - 48 * x ** 2 + 12)
    assert v[2][0] == 1.0
    # D.n の定義 P_n^m(z) = ((1+z)/(1-z))^(m/2) 2F1reg(-n, n+1; 1-m; (1-z)/2) (Condon-Shortley の位相つき)
    z = 0.4
    assert _close(v[3][0], -math.sqrt(1 - z * z) * (15 * z * z - 3) / 2.0, 1e-13)       # P_3^1 = -(1-z²)^(1/2) (15z²-3)/2 = +.2749545…
    z = 0.3
    assert _close(v[4][0], 3 * (1 - z * z), 1e-13)                                # P_2^2 = 3(1-z²)
    z = 0.4
    assert _close(v[5][0], (5 * z ** 3 - 3 * z) / 2.0, 1e-13)
    z = -0.3
    assert _close(v[6][0], -105 * z * (1 - z * z) ** 1.5, 1e-13)                  # P_4^3 = -105 z (1-z²)^(3/2)
    assert v[7][0] == 6.0 and v[8][0] == 10.0                                     # (-3)(-4)/2 = 6 (K64: NaN)
    assert _close(v[9][0], 0.5 * -0.5 * -1.5 / 6.0)


def test_b45_legendre_integer_m_is_limit_of_formula():
    """整数 m の値は、D.n の式の m → 整数 の極限と同じ符号・値 (奇数 m の符号の確認)"""
    from pysad.specfn import _legendre_p
    for n, m, z in ((3, 1, 0.4), (4, 3, -0.3), (5, 1, 0.2), (4, 2, 0.6)):
        lo = _legendre_p(n, m - 1e-6, z)
        hi = _legendre_p(n, m + 1e-6, z)
        assert abs(_legendre_p(n, m, z) - 0.5 * (lo + hi)) < 1e-8 * max(1.0, abs(lo))


def test_b46_integer_order_u_bessely_besselk():
    v = _vals(['HypergeometricU[1,2,.5]', 'BesselY[0,.1]', 'BesselK[1.5,.7]', 'BesselJ[-1,1.5]', 'BesselY[1,2.5]', 'BesselJZero[0,1]',
               'BesselJZero[1,2]', 'BesselK[0,.3]', 'BesselY[0,2.5]'])
    assert _close(v[0][0], 2.0)                                                   # U(1,2,z) = 1/z (K64: INF)
    assert _close(v[1][0], _y0_series(0.1), 1e-13)                                # K64: NaN
    assert _close(v[2][0], math.sqrt(math.pi / 1.4) * math.exp(-0.7) * (1 + 1 / 0.7), 1e-13)     # K_{3/2} の閉じた式 (K64: INF)
    assert _close(v[3][0], -_jn_int(1, 1.5), 1e-13)                               # J_{-1} = -J_1 (K64: 応答なし)
    # Wronskian: J_1 Y_0 - J_0 Y_1 = 2/(πx)
    x = 2.5
    y1 = (_jn_int(1, x) * _y0_series(x) - 2.0 / (math.pi * x)) / _jn_int(0, x)
    assert _close(v[4][0], y1, 1e-12)                                             # K64: 応答なし
    assert _close(v[8][0], _y0_series(x), 1e-12)
    assert abs(_jn_int(0, v[5][0])) < 1e-14 and 2.40 < v[5][0] < 2.41            # J_0 の最初の零点
    assert abs(_jn_int(1, v[6][0])) < 1e-14 and 7.0 < v[6][0] < 7.03             # J_1 の 2 番目の零点
    assert _close(v[7][0], _knu_int(0.0, 0.3), 1e-11)


def test_b53_large_argument_bessel():
    xs = [33.0, 100.0, 1000.0]
    v = _vals(['BesselJ[0,%r]' % x for x in xs] + ['BesselJ[2,%r]' % x for x in xs] + ['BesselK[0,12.5]', 'BesselK[1,40.]', 'BesselY[0,33.]',
                                                                                    'BesselI[0,30.]'])
    for k, x in enumerate(xs):
        assert _close(v[k][0], _jn_int(0, x, 20000), 0, 2e-14)                    # K64: BesselJ[0,33.] は相対 3e-3 の誤差
        assert _close(v[3 + k][0], _jn_int(2, x, 20000), 0, 2e-14)
    assert _close(v[6][0], _knu_int(0.0, 12.5), 1e-11)                            # K64: 相対 1e-5 の誤差
    assert _close(v[7][0], _knu_int(1.0, 40.0), 1e-11)
    # Y_0(x) = (2/π) ∫… の代わりに Wronskian を J と Y の漸近形でなく、J_1 Y_0 - J_0 Y_1 = 2/(πx) を満たすことだけ確認
    y1 = _vals(['BesselY[1,33.]'])[0][0]
    assert abs(_jn_int(1, 33.0, 20000) * v[8][0] - _jn_int(0, 33.0, 20000) * y1 - 2.0 / (math.pi * 33.0)) < 1e-14
    # I_0(x) = (1/π) ∫_0^π exp(x cosθ) dθ
    th = (np.arange(20000) + 0.5) * math.pi / 20000
    assert _close(v[9][0], float(np.mean(np.exp(30.0 * np.cos(th)))), 1e-12)


def test_b53_large_argument_fresnel_airy():
    v = _vals(['FresnelS[6.1]', 'FresnelC[6.1]', 'FresnelS[20.]', 'FresnelC[20.]', 'AiryAi[15.5]', 'AiryAi[3.]', 'AiryBi[3.]', 'AiryAi[-12.]'])
    s, c = _fresnel_quad(6.1)
    assert _close(v[0][0], s, 0, 1e-11) and _close(v[1][0], c, 0, 1e-11)          # K64: 1.3e7、-2.1e7 (桁落ち)
    s, c = _fresnel_quad(20.0)
    assert _close(v[2][0], s, 0, 1e-10) and _close(v[3][0], c, 0, 1e-10)
    assert abs(v[0][0] - 0.5) < 0.06 and abs(v[1][0] - 0.5) < 0.06                # 大きい x で 1/2 に近づく

    def ai(x):                                                                    # Ai(x) = (1/π) sqrt(x/3) K_{1/3}(2/3 x^{3/2})
        return math.sqrt(x / 3.0) / math.pi * _knu_int(1.0 / 3.0, 2.0 / 3.0 * x ** 1.5)
    assert _close(v[4][0], ai(15.5), 1e-10)                                       # K64: -16
    assert _close(v[5][0], ai(3.0), 1e-10)
    # Ai Bi' - Ai' Bi = 1/π
    w = _vals(['AiryAi[3.]*AiryBiPrime[3.]-AiryAiPrime[3.]*AiryBi[3.]', 'AiryAi[-12.]*AiryBiPrime[-12.]-AiryAiPrime[-12.]*AiryBi[-12.]'])
    assert _close(w[0][0], 1.0 / math.pi, 1e-12) and _close(w[1][0], 1.0 / math.pi, 1e-11)
    # Ai(-x) = (1/π) ∫_0^∞ cos(t³/3 - x t) dt は収束が遅いので、Airy の微分方程式 y'' = x y を差分で確認
    h = 1e-3
    y = _vals(['AiryAi[%r]' % (-12.0 + d) for d in (-h, 0.0, h)])
    assert abs((y[0][0] - 2 * y[1][0] + y[2][0]) / h ** 2 - (-12.0) * y[1][0]) < 1e-4


def test_b54_polylog_real():
    v = _vals(['PolyLog[2,.5]', 'PolyLog[2,-1.]', 'PolyLog[3,.5]', 'PolyLog[2,.9]', 'PolyLog[2,2.]', 'PolyLog[1.5,.3]'])
    assert len(v[0]) == 1 and _close(v[0][0], math.pi ** 2 / 12 - math.log(2) ** 2 / 2, 1e-13)      # 実数 (K64 は虚部 1e-15 の雑音つきの複素数)
    assert _close(v[1][0], -math.pi ** 2 / 12, 1e-13)
    assert _close(v[2][0], 7.0 / 8.0 * 1.2020569031595942 - math.pi ** 2 * math.log(2) / 12 + math.log(2) ** 3 / 6, 1e-13)
    # Li_2(z) = -∫_0^z ln(1-t)/t dt
    li2 = -integrate.quad(lambda t: math.log1p(-t) / t if t else -1.0, 0.0, 0.9, epsabs=0, epsrel=1e-13)[0]
    assert len(v[3]) == 1 and _close(v[3][0], li2, 1e-11)
    # z > 1: Li_2(2) = π²/4 - iπ ln 2 (切断線の下側)
    assert _close(v[4][0], math.pi ** 2 / 4, 1e-12) and _close(v[4][1], -math.pi * math.log(2), 1e-12)
    assert _close(v[5][0], sum(0.3 ** k / k ** 1.5 for k in range(1, 200)), 1e-13)


def test_b47_inverseerf_v2_true_value():
    xs = [0.3, 0.999, -0.7, 1e-8, 0.8, 0.9999999]
    v = _vals(['InverseErf[%r]' % x for x in xs] + ['InverseErfc[.2]', 'InverseErfc[1e-30]', 'InverseErf[1]', 'InverseErf[1.5]'])
    for x, r in zip(xs, v):
        y = r[0]
        # 独立な Newton 法 (正しい係数 sqrt(π)/2 e^{y²})
        t = y
        for _ in range(4):
            t -= (math.erf(t) - x) * math.sqrt(math.pi) / 2.0 * math.exp(t * t)
        assert _close(y, t, 1e-13, 1e-300), (x, y, t)
        assert abs(math.erf(y) - x) < 4e-16
    assert _close(v[0][0], 0.2724627147267543, 1e-13)                            # K64: .272462711811 (相対 1e-8 の誤差)
    assert abs(math.erfc(v[6][0]) - 0.2) < 1e-15
    assert abs(math.erfc(v[7][0]) / 1e-30 - 1.0) < 1e-12                          # 小さい z でも桁落ちしない
    assert v[8][0] == float('inf') and v[9][0] != v[9][0]


def test_b47_inverseerf_k64_unchanged():
    """k64 モードは K64 と同じ値 (係数 sqrt(π/2) と早い打ち切り)"""
    v = _vals(['InverseErf[.3]', 'InverseErf[.999]', 'InverseErfc[.2]'], mode='k64')
    assert abs(v[0][0] - .272462711811) < 1e-12
    assert abs(v[1][0] - 2.326753739975) < 1e-11
    from pysad import compat
    from pysad.specfn import inverseerf
    old = compat._mode
    compat.set_mode('k64')
    try:
        assert v[2][0] == pytest.approx(inverseerf(0.8), rel=1e-15)
    finally:
        compat.set_mode(old)


def test_fitgaussian_uses_true_inverseerf_in_v2(v2_mode):
    """FitGaussian (Fit.n) は InverseErf を使う。v2 では正規分布の分位点から σ が 1e-9 以内で戻る (K64 は 8 桁目がずれる)"""
    from pysad.ffscmd import run_sad_string
    from scipy.special import erfinv
    n = 400
    xs = [math.sqrt(2.0) * 1.7 * float(erfinv(2.0 * (k + 0.5) / n - 1.0)) + 0.3 for k in range(n)]
    buf = io.StringIO()
    run_sad_string('FFS;\n$FORM="25.17";\nr=FitGaussian[{%s}];\nPrint["S= ",r[[1]]," ",r[[2]]];\nExit;\n' % ','.join(repr(x) for x in xs), out=buf)
    m = re.search(r'S= (.*)', buf.getvalue())
    s, mu = [float(t) for t in re.findall(_NUM, m.group(1))][:2]
    assert abs(s - 1.7) < 0.05 and abs(mu - 0.3) < 1e-6


def _cplx(r):
    return complex(r[0], r[1] if len(r) > 1 else 0.0)


def test_b54_polylog_large_argument_noninteger_order():
    """非整数次の PolyLog は |z| の上限なし (以前は |z| ≥ 150 で NaN、z < -12・z > 60 で展開の精度が落ちていた)。
    根拠: Fermi–Dirac 積分 -Li_s(-y) = 1/Γ(s) ∫ t^{s-1}/(e^t/y + 1) dt (s > 0)、Im Li_s(x > 1) = -π ln^{s-1} x / Γ(s)、Li_{s-1} = z dLi_s/dz"""
    def fd(s, y):
        w = math.log(y)
        f = lambda t: t ** (s - 1) / (math.exp(t - w) + 1.0)
        a = integrate.quad(f, 0, w, epsabs=0, epsrel=1e-13, limit=400)[0] + integrate.quad(f, w, w + 80, epsabs=0, epsrel=1e-13, limit=400)[0]
        return -a / math.gamma(s)
    cases = [(s, z) for s in (0.5, 1.5, 2.5, 3.7, 7.9) for z in (-20.0, -140.0, -200.0, -1e6)]
    v = _vals(['PolyLog[%r,%r]' % c for c in cases])
    for (s, z), r in zip(cases, v):
        assert len(r) == 1, (s, z, r)                                            # z < -1 は実数
        assert _close(r[0], fd(s, -z), 1e-11), (s, z, r)
    cases = [(s, z) for s in (0.5, 1.5, 2.5, 3.7, 7.9) for z in (80.0, 200.0, 1e6, 1e12)]
    v = _vals(['PolyLog[%r,%r]' % c for c in cases])
    for (s, z), r in zip(cases, v):
        assert _close(r[1], -math.pi * math.log(z) ** (s - 1) / math.gamma(s), 1e-11), (s, z, r)
    # 実部: 従来の展開が正確な z = 80 と一致することは上の run で使われる経路の切り替え (z > 60) の両側で連続
    a = _vals(['PolyLog[2.5,59.99]', 'PolyLog[2.5,60.01]', 'PolyLog[.5,-11.99]', 'PolyLog[.5,-12.01]', 'PolyLog[3.7,59.99]', 'PolyLog[3.7,60.01]'])
    for k in (0, 2, 4):
        d = abs(_cplx(a[k + 1]) - _cplx(a[k])) / abs(_cplx(a[k]))
        assert d < 2e-3 and d > 0                                                 # 0.02 だけ離れた点: 滑らか (段差なし)
    # s < 0: Li_{s-1}(z) = z d/dz Li_s(z) (中心差分)
    for s, z in ((0.5, -200.0), (-0.5, 300.0), (-1.5, -1e4)):
        h = 1e-4 * abs(z)
        w = _vals(['PolyLog[%r,%r]' % (s, z + h), 'PolyLog[%r,%r]' % (s, z - h), 'PolyLog[%r,%r]' % (s - 1, z)])
        d = (_cplx(w[0]) - _cplx(w[1])) / (2 * h) * z
        assert abs(d - _cplx(w[2])) < 1e-7 * abs(d), (s, z)


def test_complex_argument_bessel_hypergeometric():
    """複素引数の Bessel・超幾何・Airy (以前は TypeError)。根拠: 積分表示・級数・Wronskian"""
    v = _vals(['BesselJ[0,1.+2.*I]', 'BesselJ[3,2.-1.5*I]', 'BesselK[0.5,2.-1.*I]', 'BesselK[0,1.+1.*I]', 'BesselI[0,1.5+.5*I]',
               'BesselY[0,.5+.5*I]', 'BesselY[1,.5+.5*I]', 'BesselJ[1,.5+.5*I]', 'BesselJ[0,.5+.5*I]',
               'Hypergeometric2F1[.5,1.5,2,.3+.4*I]', 'Hypergeometric1F1[1,2,1.+1.*I]', 'Hypergeometric0F1[1.5,-2.+1.*I]',
               'AiryAi[1.+1.*I]', 'AiryAiPrime[1.+1.*I]', 'AiryBi[1.+1.*I]', 'AiryBiPrime[1.+1.*I]'])
    th = (np.arange(4000) + 0.5) * math.pi / 4000

    def jn(n, z):
        return complex(np.mean(np.cos(n * th - z * np.sin(th))))
    assert abs(_cplx(v[0]) - jn(0, 1 + 2j)) < 1e-13
    assert abs(_cplx(v[1]) - jn(3, 2 - 1.5j)) < 1e-13
    z = 2 - 1j
    assert abs(_cplx(v[2]) - cmath_sqrt(math.pi / (2 * z)) * np.exp(-z)) < 1e-13         # K_{1/2}(z) = sqrt(π/2z) e^-z
    z = 1 + 1j
    k0 = complex(integrate.quad(lambda t: (np.exp(-z * math.cosh(t))).real, 0, 30, epsabs=0, epsrel=1e-13)[0],
                 integrate.quad(lambda t: (np.exp(-z * math.cosh(t))).imag, 0, 30, epsabs=0, epsrel=1e-13)[0])
    assert abs(_cplx(v[3]) - k0) < 1e-12
    assert abs(_cplx(v[4]) - complex(np.mean(np.exp((1.5 + 0.5j) * np.cos(th))))) < 1e-13   # I_0(z) = (1/π) ∫ exp(z cosθ) dθ
    z = 0.5 + 0.5j
    assert abs(_cplx(v[7]) * _cplx(v[5]) - _cplx(v[8]) * _cplx(v[6]) - 2.0 / (math.pi * z)) < 1e-13      # J_1 Y_0 - J_0 Y_1 = 2/(πz)
    assert abs(_cplx(v[7]) - jn(1, z)) < 1e-13 and abs(_cplx(v[8]) - jn(0, z)) < 1e-13

    def series(a, b, z):
        t, tot = 1.0 + 0j, 1.0 + 0j
        for k in range(400):
            num = 1.0
            for x in a:
                num *= x + k
            den = k + 1.0
            for x in b:
                den *= x + k
            t *= num / den * z
            tot += t
        return tot
    assert abs(_cplx(v[9]) - series([.5, 1.5], [2], .3 + .4j)) < 1e-13
    assert abs(_cplx(v[10]) - series([1], [2], 1 + 1j)) < 1e-13
    assert abs(_cplx(v[11]) - series([], [1.5], -2 + 1j)) < 1e-13
    # Airy: Ai Bi' - Ai' Bi = 1/π、y'' = z y
    assert abs(_cplx(v[12]) * _cplx(v[15]) - _cplx(v[13]) * _cplx(v[14]) - 1.0 / math.pi) < 1e-13
    h = 1e-4
    y = _vals(['AiryAi[%r+1.*I]' % (1.0 + d) for d in (-h, 0.0, h)])
    assert abs((_cplx(y[0]) - 2 * _cplx(y[1]) + _cplx(y[2])) / h ** 2 - (1 + 1j) * _cplx(y[1])) < 1e-6
    # 実数の引数の値は変わらない
    r = _vals(['BesselJ[0,2.]', 'Hypergeometric2F1[.5,1.5,2,.3]'])
    assert _close(r[0][0], jn(0, 2.0).real, 1e-13) and _close(r[1][0], series([.5, 1.5], [2], .3).real, 1e-13)


def cmath_sqrt(z):
    import cmath
    return cmath.sqrt(z)
