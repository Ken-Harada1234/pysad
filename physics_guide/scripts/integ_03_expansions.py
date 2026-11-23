# -*- coding: utf-8 -*-
r"""03 章の途中式の確認 (本文に書き足した展開)。

(a) §5: Yoshida の係数 (3.5)(3.6) の代数 (数値)。一般化 (p 次から p+2 次へ) の係数 2 w1^{p+1} + w0^{p+1} = 0。
(b) §8.1: A = -:T:, B = -:V: (T = p^2/2、1 自由度) での [A,B]、W = [B,[A,B]]、W5 = [W,[A,B]] の関数表現。
    ([f,g] = ∂f/∂x ∂g/∂p - ∂f/∂p ∂g/∂x。A, B, W, W5 は :f: の形に書いた f の関数)
(c) §8.2: (3.11) の係数の代数 (a = 1/2 - t を代入)。
(d) §8.3: 力勾配型の 1 刻みの誤差 (V = x^4/4、50 桁) / h^5 を、(3.13) が作るベクトル場と比べる。
"""
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)

# (a) ------------------------------------------------------------------
print('(a) Yoshida の係数')
w1 = 1 / (2 - sp.cbrt(2))
w0 = 1 - 2 * w1
print('  2^(1/3) = %.10f、2 - 2^(1/3) = %.10f、w1 = %.10f、w0 = %.10f'
      % (float(sp.cbrt(2)), float(2 - sp.cbrt(2)), float(w1), float(w0)))
print('  2 w1 + w0 = %.12f、2 w1^3 + w0^3 = %.3e' % (float(2 * w1 + w0), float(2 * w1**3 + w0**3)))
print('  (1 - 2 w1)^3 = -2 w1^3 から 1 - 2 w1 = -2^(1/3) w1:  %.10f = %.10f'
      % (float(1 - 2 * w1), float(-sp.cbrt(2) * w1)))
w1_6 = 1 / (2 - 2**(sp.Rational(1, 5)))
print('  p = 4 の次 (6 次): w1 = %.10f、2 w1^5 + w0^5 = %.2e' % (float(w1_6), float(2 * w1_6**5 + (1 - 2 * w1_6)**5)))

# (b) ------------------------------------------------------------------
print('(b) 関数表現 (1 自由度。[f,g] = f_x g_p - f_p g_x)')
x, p = sp.symbols('x p')
V = sp.Function('V')(x)
T = p**2 / 2


def br(f, g):
    return sp.simplify(sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x))


# A = -:T:, B = -:V:  ->  [A,B] = [:T:, :V:] = :[T,V]:
AB = br(T, V)
print('  [A,B] = :f_AB:、f_AB = [T,V] =', AB)
# W = [B,[A,B]] = [-:V:, :f_AB:] = -:[V, f_AB]:
W = -br(V, AB)
print('  W = [B,[A,B]] = :f_W:、f_W = -[V, f_AB] =', sp.simplify(W))
# W5 = [W,[A,B]] = [:f_W:, :f_AB:] = :[f_W, f_AB]:
W5 = br(W, AB)
print('  W5 = [W,[A,B]] = :f_W5:、f_W5 = [f_W, f_AB] =', sp.simplify(W5))

# (c) ------------------------------------------------------------------
print('(c) (3.11) の係数')
t, c3 = sp.symbols('t c3')
a = sp.Rational(1, 2) - t
cAAB = t**2 / 3 - a**2 / 6 - a * t / 3
cBBA = -t / 12 - c3 + a / 6
print('  [A,[A,B]] の係数 =', sp.factor(sp.simplify(cAAB)), ' (= (12t^2-1)/24:', sp.simplify(cAAB - (12 * t**2 - 1) / 24) == 0, ')')
print('  [B,[B,A]] の係数 =', sp.simplify(cBBA), ' (= 1/12 - t/4 - c3:', sp.simplify(cBBA - (sp.Rational(1, 12) - t / 4 - c3)) == 0, ')')
tt = 1 / sp.sqrt(12)
print('  t = 1/sqrt(12) で c3 = 1/12 - t/4 = %.12f、beta/2 = %.12f' % (float(sp.Rational(1, 12) - tt / 4), float((sp.Rational(1, 6) - tt / 2) / 2)))
print('  t = 1/3、c3 = 0 (CALC の tthine): [A,[A,B]] 係数 = %s、[B,[B,A]] 係数 = %s'
      % (sp.simplify(cAAB.subs(t, sp.Rational(1, 3))), sp.simplify(cBBA.subs({t: sp.Rational(1, 3), c3: 0}))))

# (d) ------------------------------------------------------------------
print('(d) 力勾配型の 1 刻みの誤差 / h^5 (H = p^2/2 + x^4/4、mpmath 50 桁)')
mp.mp.dps = 50
tt_ = 1 / mp.sqrt(12)
a_ = mp.mpf(1) / 2 - tt_
c3_ = mp.mpf(1) / 12 - tt_ / 4
c5_ = (9 - 5 * mp.sqrt(3)) / 720


def Vp(x): return x**3          # V' (V = x^4/4)
def Vpp(x): return 3 * x**2
def Vppp(x): return 6 * x
def Vpppp(x): return mp.mpf(6)


def step(x0, p0, h, c5):
    x_, p_ = x0, p0
    x_ += a_ * h * p_                       # exp(a h A): x += a h p
    p_ -= h / 2 * Vp(x_)                    # exp(h B/2): p -= (h/2) V'
    x_ += tt_ * h * p_
    # 補正キック: exp(:f:)、f = c3 h^3 V'^2 + c5 h^5 (-2 V'^2 V'')、Δp = ∂f/∂x
    dfdx = c3_ * h**3 * 2 * Vp(x_) * Vpp(x_) + c5 * h**5 * (-2) * (2 * Vp(x_) * Vpp(x_)**2 + Vp(x_)**2 * Vppp(x_))
    p_ += dfdx
    x_ += tt_ * h * p_
    p_ -= h / 2 * Vp(x_)
    x_ += a_ * h * p_
    return x_, p_


def exact(x0, p0, h):
    f = lambda s, y: [y[1], -y[0]**3]
    sol = mp.odefun(f, 0, [x0, p0])
    y = sol(h)
    return y[0], y[1]


def H5(x, pp, c5):
    s3 = mp.sqrt(3)
    Vp_, Vpp_, Vppp_, Vpppp_ = Vp(x), Vpp(x), Vppp(x), Vpppp(x)
    return ((-720 * c5 + 9 - 5 * s3) / 360 * Vp_**2 * Vpp_
            + ((s3 / 288 - mp.mpf(1) / 120) * Vp_ * Vppp_ + (mp.mpf(1) / 240 - s3 / 288) * Vpp_**2) * pp**2
            + Vpppp_ / 4320 * pp**4)


x0, p0 = mp.mpf('0.5'), mp.mpf('0.3')
for label, c5 in (('c5 あり (γ/2)', c5_), ('c5 なし', mp.mpf(0))):
    print('  --- %s' % label)
    for h in (mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01')):
        xs, ps = step(x0, p0, h, c5)
        xe, pe = exact(x0, p0, h)
        ex, ep = (xs - xe) / h**5, (ps - pe) / h**5
        # H5 の流れ (x' = ∂H5/∂p, p' = -∂H5/∂x)
        dHdp = mp.diff(lambda q: H5(x0, q, c5), p0)
        dHdx = mp.diff(lambda q: H5(q, p0, c5), x0)
        print('   h=%s: (S-exact)/h^5 = (%.6f, %.6f)  H5 のベクトル場 (∂H5/∂p, -∂H5/∂x) = (%.6f, %.6f)'
              % (mp.nstr(h, 3), float(ex), float(ep), float(dHdp), float(-dHdx)))
