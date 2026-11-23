# -*- coding: utf-8 -*-
r"""17_collective.md の途中式の確認 (本文は紙の上で追える形に書いてある。ここは検算)。

(a) 自己項 1/2: 幅 ε の対称な分布どうしの W(0+) P(s>0) → W(0+)/2 (§5.2)
(b) ガウスの畳み込みの平方完成 (sympy) と 1/r の積分表示 (§7.3)
(c) 円筒壁の Green 関数: K0 のフーリエ表示、Graf の加法定理、像の項が r=R で全体を 0 にする (§7.1)
(d) Rutherford の小角の関係 θ = 2k/(μ v² b) と dσ/dΩ の整合 (§9.2)
(e) CSR の出発点: 円軌道の角度分布 (K_{2/3}, K_{1/3}) を ψ で積分して (17.37) の F(y) に一致 (§10.1)
(f) -ik の 1/3 乗の逆変換の有限部分 (Γ の解析接続) と部分積分の形 (§10.3)
"""
import math
import _common as C

C.start(__file__)
import mpmath as mp
import sympy as sp
mp.mp.dps = 25

print('== (a) 自己項: 幅 ε の対称な分布の差 s ~ N(0, 2ε²)、W(s) = W0 exp(-s/ℓ) (s>0) ==')
for eps in (1e-1, 1e-2, 1e-3):
    ell = 1.0
    # <W> = W0 ∫_0^∞ e^{-s/ℓ} N(s;2ε²) ds
    val = mp.quad(lambda s: mp.e**(-s/ell)*mp.e**(-s**2/(4*eps**2))/(2*eps*mp.sqrt(mp.pi)), [0, mp.inf])
    print('  ε/ℓ = %g: <W>/W0 = %.8f  (1/2 との差 %.2e)' % (eps/ell, val, val-0.5))

print('\n== (b) ∫ N(x\';σ²) exp(-(x-x\')²τ²) dx\' = (1+2σ²τ²)^{-1/2} exp(-x²τ²/(1+2σ²τ²)) ==')
x, xp, s, t = sp.symbols('x xp sigma tau', positive=True)
x = sp.Symbol('x', real=True)
xp = sp.Symbol('xp', real=True)
a = 1/(2*s**2) + t**2
# 平方完成: -[a xp² - 2 τ² x xp + τ² x²] = -a (xp - τ² x/a)² - τ² x² (1 - τ²/a)
expo = -xp**2/(2*s**2) - (x-xp)**2*t**2
chk1 = sp.simplify(expo - (-a*(xp - t**2*x/a)**2 - t**2*x**2*(1 - t**2/a)))
print('  平方完成の恒等式の差:', chk1)
print('  1 - τ²/a - 1/(1+2σ²τ²) =', sp.simplify(1 - t**2/a - 1/(1+2*s**2*t**2)))
print('  ガウス積分の前因子 sqrt(π/a)/sqrt(2πσ²) - (1+2σ²τ²)^{-1/2} =', sp.simplify(sp.sqrt(sp.pi/a)/sp.sqrt(2*sp.pi*s**2) - (1+2*s**2*t**2)**sp.Rational(-1, 2)))
v = mp.quad(lambda u: mp.e**(-u**2/(2*0.7**2))/(0.7*mp.sqrt(2*mp.pi))*mp.e**(-(1.3-u)**2*0.9**2), [-mp.inf, mp.inf])
w = (1+2*0.7**2*0.9**2)**-0.5*mp.e**(-1.3**2*0.9**2/(1+2*0.7**2*0.9**2))
print('  数値 (σ=0.7, τ=0.9, x=1.3): 積分 %.12f  閉じた式 %.12f' % (v, w))
print('  ∫_0^∞ exp(-r²τ²)dτ = √π/(2r):', sp.integrate(sp.exp(-sp.Symbol('r', positive=True)**2*t**2), (t, 0, sp.oo)))

print('\n== (c) 円筒壁の Green 関数 ==')
d = mp.mpf('0.8'); k = mp.mpf('1.7')
ft = 2*mp.quad(lambda z: mp.cos(k*z)/mp.sqrt(d**2+z**2), [0, mp.inf])  # 振動積分
ft = 2*mp.quadosc(lambda z: mp.cos(k*z)/mp.sqrt(d**2+z**2), [0, mp.inf], omega=k)
print('  ∫dζ e^{-ikζ}/sqrt(d²+ζ²) = %.10f   2K0(kd) = %.10f' % (ft, 2*mp.besselk(0, k*d)))
r, rp, R, dth = mp.mpf('0.6'), mp.mpf('0.9'), mp.mpf('1.5'), mp.mpf('0.7')
dist = mp.sqrt(r**2+rp**2-2*r*rp*mp.cos(dth))
lhs = mp.besselk(0, k*dist)
rhs = mp.nsum(lambda m: (1 if m == 0 else 2)*mp.besselk(int(m), k*max(r, rp))*mp.besseli(int(m), k*min(r, rp))*mp.cos(m*dth), [0, 80])
print('  Graf: K0(k|r-r\'|) = %.12f   Σ ε_m K_m I_m cos = %.12f' % (lhs, rhs))
def G(r, rp, dth, k, R, M=60):
    tot = 0
    for m in range(M):
        em = 1 if m == 0 else 2
        free = 2*em*mp.besselk(m, k*max(r, rp))*mp.besseli(m, k*min(r, rp))
        img = -2*em*mp.besselk(m, k*R)/mp.besseli(m, k*R)*mp.besseli(m, k*r)*mp.besseli(m, k*rp)
        tot += (free+img)*mp.cos(m*dth)
    return tot
print('  r = R での全体 G(R, r\', Δθ): %.3e   (r\'=0.9, Δθ=0.7)' % G(R, rp, dth, k, R))
print('  項ごと (m=0..3) の free + image の和:', [mp.nstr(2*(1 if m == 0 else 2)*mp.besselk(m, k*R)*mp.besseli(m, k*rp) - 2*(1 if m == 0 else 2)*mp.besselk(m, k*R)/mp.besseli(m, k*R)*mp.besseli(m, k*R)*mp.besseli(m, k*rp), 3) for m in range(4)])
print('  r < R の内側 (r=0.6) の値 G = %.8f、自由空間 2K0 = %.8f (差は像の項: 正則な調和関数)' % (G(r, rp, dth, k, R), 2*lhs))

print('\n== (d) Rutherford の小角: θ = 2k/(μ v² b)、dσ/dΩ ==')
b_, th, kk, mu, vv, be, qq = sp.symbols('b theta k mu v beta q', positive=True)
b_exact = kk/(mu*vv**2)/sp.tan(th/2)   # b = (k/μv²) cot(θ/2) (Rutherford、斥力 k/r)
print('  b θ (θ→0) =', sp.limit(b_exact*th, th, 0), ' (= 2k/(μ v²))')
sub = {mu: sp.Symbol('m')/2, vv: 2*be*sp.Symbol('c'), kk: sp.Symbol('r_e')*sp.Symbol('m')*sp.Symbol('c')**2}
print('  2k/(μv²) =', sp.simplify((2*kk/(mu*vv**2)).subs(sub)), '  (= r_e/β*²)')
# dσ = 2π b |db|、b = r_e/(β*² θ) → dσ/dΩ
bb = sp.Symbol('r_e')/(be**2*th)
re_p = sp.Symbol('r_e', positive=True)
bb = re_p/(be**2*th)
dsig = 2*sp.pi*bb*(-sp.diff(bb, th))   # |db| = -db
dsdo = sp.simplify(dsig/(2*sp.pi*th))
print('  dσ/dΩ = 2π b|db/dθ|/(2π θ) =', dsdo, ' ; 16 r_e²/(|q|⁴θ⁴) with β*=|q|/2:', sp.simplify(dsdo.subs(be, qq/2) - 16*re_p**2/(qq**4*th**4)))

print('\n== (e) CSR: 角度分布 → (17.37) ==')
# Jackson の円軌道の d²I/dωdΩ (SI、ω>0)。ξ = (ωR/3c)(1/γ²+ψ²)^{3/2}
# 1 周のエネルギー = 2π ∫dψ d²I/dωdΩ、パワー = (c/2πR) × それ → P_ω = (c/R) ∫dψ d²I/dωdΩ
gam = mp.mpf(50); Rr = mp.mpf(1); cc = mp.mpf(1); e2 = mp.mpf(1)   # e²/(ε0) を 1 とする単位
def d2I(w, psi):
    g2 = 1/gam**2 + psi**2
    xi = (w*Rr/(3*cc))*g2**1.5
    return e2/(12*mp.pi**3*cc)*(w*Rr/cc)**2*g2**2*(mp.besselk(mp.mpf(2)/3, xi)**2 + psi**2/g2*mp.besselk(mp.mpf(1)/3, xi)**2)
wc = 1.5*gam**3*cc/Rr
for y in (0.05, 0.3, 1.0, 3.0):
    w = y*wc
    Pnum = (cc/Rr)*2*mp.quad(lambda p: d2I(w, p), [0, 1/gam, 5/gam, 40/gam])
    Fy = y*mp.quad(lambda t: mp.besselk(mp.mpf(5)/3, t), [y, y+1, mp.inf])
    P37 = mp.sqrt(3)*e2*gam/(8*mp.pi**2*Rr)*Fy
    print('  y = %.2f: ψ 積分 %.10e   (17.37) %.10e   比 %.8f' % (y, Pnum, P37, Pnum/P37))
print('  (SI の e²/ε0 を 1 とした単位。Jackson の (14.79) の e²→e²/(4πε0) を入れた形)')
xi = mp.mpf('0.7')
z_ = (1.5*xi)**(mp.mpf(2)/3)
I1 = (1.5*xi)**(-mp.mpf(1)/3)*2*mp.pi*mp.airyai(z_)
I2 = -(1.5*xi)**(-mp.mpf(2)/3)*2*mp.pi*mp.airyai(z_, derivative=1)
print("  Airy 経由: int cos[(3/2)xi(x+x^3/3)]dx = (3xi/2)^(-1/3) 2pi Ai(z) = %.10f   (2/sqrt3)K_{1/3}(xi) = %.10f" % (I1, 2/mp.sqrt(3)*mp.besselk(mp.mpf(1)/3, xi)))
print("  Airy 経由: int x sin[(3/2)xi(x+x^3/3)]dx = -(3xi/2)^(-2/3) 2pi Ai'(z) = %.10f   (2/sqrt3)K_{2/3}(xi) = %.10f" % (I2, 2/mp.sqrt(3)*mp.besselk(mp.mpf(2)/3, xi)))
ep = mp.mpf('1e-3')
J1 = 2*mp.quad(lambda u: mp.cos(1.5*xi*(u+u**3/3))*mp.e**(-ep*u**2), mp.linspace(0, 60, 400))
print("  直接の数値積分 (収束因子 e^{-1e-3 x^2}): int cos = %.5f (Airy 経由 %.5f)" % (J1, I1))


print('\n== (f) (-ik)^{1/3} の逆変換: 有限部分 ==')
for p in (mp.mpf('0.5'), mp.mpc('0.3', -1.2)):
    lhs = mp.quad(lambda s_: s_**(-mp.mpf(4)/3)*(mp.e**(-p*s_) - 1), [0, 1, mp.inf])
    rhs = mp.gamma(-mp.mpf(1)/3)*p**(mp.mpf(1)/3)
    print('  p = %s: ∫s^{-4/3}(e^{-ps}-1)ds = %s,  Γ(-1/3) p^{1/3} = %s' % (mp.nstr(p, 4), mp.nstr(lhs, 12), mp.nstr(rhs, 12)))
print('  Γ(-1/3) + 3Γ(2/3) =', mp.nstr(mp.gamma(-mp.mpf(1)/3) + 3*mp.gamma(mp.mpf(2)/3), 3))
sig = mp.mpf(1)
lam = lambda z: mp.e**(-z**2/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
dlam = lambda z: -z/sig**2*lam(z)
for z0 in (mp.mpf(2), mp.mpf(-2), mp.mpf(0)):
    fp = mp.quad(lambda s_: s_**(-mp.mpf(4)/3)*(lam(z0+s_) - lam(z0)), [0, 1, 10, mp.inf]) * (-1/(3*mp.gamma(mp.mpf(2)/3)))
    der = -1/mp.gamma(mp.mpf(2)/3)*mp.quad(lambda s_: s_**(-mp.mpf(1)/3)*dlam(z0+s_), [0, 1, 10, mp.inf])
    print('  z=%s: ∫w(s)λ(z+s)ds (有限部分) = %s   部分積分の形 = %s' % (z0, mp.nstr(fp, 12), mp.nstr(der, 12)))
