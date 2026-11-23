# -*- coding: utf-8 -*-
r"""BEND の章 (05_bend.md) で、途中式を足した箇所の確認 (sympy・mpmath。pysad は使わない)。

(a) hi = (H2 + p)/p と、保存量 C_x・C_y との関係 (§4.2)
(b) z の更新 (5.11) を、厳密解の積分 (mpmath の数値積分、30 桁) と比べる (§4.2)
(c) 面までの直進 (5.20)(5.21)、出口の面の内積 E·N2 と E·e_f2 (§5.2、§5.3)、Theta の式
(d) 設計軌道のまわりの 1 次の展開: 入口・円弧・出口の各段の 2x2 行列と、その積 (5.25) (§5.4)
(e) tbendcore の t2・t3・t2t3・px3・dpz32・x3・da の意味 (§7): Fortran の式を転記して、幾何の式と 40 桁で比べる
(f) 端の磁場 (5.30)(5.31): ラプラス方程式と回転、力の式 (5.32) (§6.1)
(g) 硬い端の Maxwell フリンジ (5.33)〜(5.36): スパイクの 1 次の摂動による Delta x・Delta z、C 付録の括弧での exp(:f:) (§6.2)
(h) F1 の 3 つの積分 (直線の立ち上がり・tanh)、部分積分 (§6.3)、(5.42)
(i) L' の補正 (5.43) を、面の変位 ΔX = -F1^2/(24 rho) の幾何 (弦の長さの条件) から出す (§6.4)
(j) rbh・rbl の 2 点ガウス・ルジャンドルの階段近似としての意味 (§6.4)
(k) tbendicorr の各行 (残りの運動学の流れ)・K1 の非線形の蹴りの係数 (§8)
(l) 分散 R16・R26・R56 の導出 (§4.1)
"""
import math
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)
mp.mp.dps = 40


def rnd(seed, lo, hi):
    return mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * mp.mpf((seed * 7919 % 10007)) / 10007


print('== (a) hi = (H2 + p)/p、C_x、C_y との関係 ==')
x, px, y, py, p, h, hb, K = sp.symbols('x p_x y p_y p h h_b K', real=True)
H2 = -p + (hb - h * p) * x + h * hb * x**2 / 2 + K * (x**2 - y**2) / 2 + (px**2 + py**2) / (2 * p)
qx, qy = sp.symbols('q_x q_y', real=True)
akxsq = -(h * hb + K) / p
akysq = K / p
drhop = h - hb / p
hi = sp.Rational(1, 2) * (qx**2 + qy**2 - akxsq * x**2 - akysq * y**2) - drhop * x
d = sp.simplify(((H2 + p).subs({px: p * qx, py: p * qy}) / p) - hi)
print('  (H2 + p)/p - hi  (px = p qx, py = p qy を代入) =', d)
w2, xe = sp.symbols('omega2 x_e', positive=True)
Cx = qx**2 + w2 * (x - xe)**2
Cy = qy**2 - akysq * y**2
hi_alt = Cx / 2 + Cy / 2 - w2 * xe**2 / 2
sub = {w2: -akxsq, xe: drhop / (-akxsq)}
print('  hi - (C_x/2 + C_y/2 - omega_x^2 x_e^2/2) =', sp.simplify((hi - hi_alt).subs(sub)))

print('\n== (b) z の更新 (5.11) と、厳密解の数値積分 (30 桁) ==')
rho0, rhob, pp, Kq, L = mp.mpf('10'), mp.mpf('9.7'), mp.mpf('1.01'), mp.mpf('0.31'), mp.mpf('1.3')
for Kq in (mp.mpf('0.0'), mp.mpf('0.31'), mp.mpf('-0.4')):
    hh, hbb = 1 / rho0, 1 / rhob
    rhoe = rhob * pp
    ak_x = -(hh * hbb + Kq) / pp
    ak_y = Kq / pp
    dr = hh - hbb / pp
    dv = mp.mpf('-0.0001')
    x0, q0, y0, qy0 = mp.mpf('0.002'), mp.mpf('0.001'), mp.mpf('0.0015'), mp.mpf('-0.0007')
    wx = mp.sqrt(-ak_x)
    xe_ = dr / (-ak_x)
    if ak_y < 0:
        wy = mp.sqrt(-ak_y)
        yf = lambda s: y0 * mp.cos(wy * s) + qy0 / wy * mp.sin(wy * s)
        qyf = lambda s: -wy * y0 * mp.sin(wy * s) + qy0 * mp.cos(wy * s)
    elif ak_y > 0:
        wy = mp.sqrt(ak_y)
        yf = lambda s: y0 * mp.cosh(wy * s) + qy0 / wy * mp.sinh(wy * s)
        qyf = lambda s: wy * y0 * mp.sinh(wy * s) + qy0 * mp.cosh(wy * s)
    else:
        yf = lambda s: y0 + qy0 * s
        qyf = lambda s: qy0
    xf = lambda s: xe_ + (x0 - xe_) * mp.cos(wx * s) + q0 / wx * mp.sin(wx * s)
    qf = lambda s: -wx * (x0 - xe_) * mp.sin(wx * s) + q0 * mp.cos(wx * s)
    zint = mp.quad(lambda s: -dv - hh * xf(s) - (qf(s)**2 + qyf(s)**2) / 2, [0, L])
    # (5.11)
    phix = wx * L
    sxkx = mp.sin(phix) / wx
    dcxkx = (1 - mp.cos(phix)) / wx**2
    xsxkx = (L - sxkx) / ak_x
    T = dr * xsxkx - sxkx * x0 - dcxkx * q0
    dxf = xf(L) - x0
    dqf = qf(L) - q0
    dyf = yf(L) - y0
    dqyf = qyf(L) - qy0
    hi_ = (q0**2 + qy0**2 - ak_x * x0**2 - ak_y * y0**2) / 2 - dr * x0
    dz = -(dxf * q0 + x0 * dqf + dxf * dqf + dyf * qy0 + y0 * dqyf + dyf * dqyf) / 4 \
        + (5 / rho0 - 1 / rhoe) * T / 4 - hi_ * L / 2 - dv * L
    print('  K = %-6s  数値積分 %s  (5.11) %s  差 %s' % (mp.nstr(Kq, 3), mp.nstr(zint, 20), mp.nstr(dz, 20), mp.nstr(zint - dz, 3)))
    # 途中: -∫x = T、∫q^2 = (C l + Δ(uq))/2
    ix = mp.quad(xf, [0, L])
    print('     -∫x dl - T = %s' % mp.nstr(-ix - T, 3))
    Cxv = q0**2 + wx**2 * (x0 - xe_)**2
    uq = lambda s: (xf(s) - xe_) * qf(s)
    iq2 = mp.quad(lambda s: qf(s)**2, [0, L])
    print('     ∫q_x^2 dl - (C_x L + Δ(u q_x))/2 = %s' % mp.nstr(iq2 - (Cxv * L + uq(L) - uq(0)) / 2, 3))

print('\n== (c) 面までの直進・出口の面の内積・Θ (sympy) ==')
psi1, psi2, phi, rho = sp.symbols('psi1 psi2 phi rho', real=True)
omega = phi - psi1 - psi2
E = sp.Matrix([-rho * (1 - sp.cos(phi)), rho * sp.sin(phi)])
ef2 = sp.Matrix([sp.cos(phi - psi2), sp.sin(phi - psi2)])
N2 = sp.Matrix([-sp.sin(phi - psi2), sp.cos(phi - psi2)])
ex1 = sp.Matrix([sp.cos(psi1), sp.sin(psi1)])
ez1 = sp.Matrix([-sp.sin(psi1), sp.cos(psi1)])
print('  E·e_f2 - rho[cos psi2 - cos(omega+psi1)] =', sp.simplify((E.dot(ef2) - rho * (sp.cos(psi2) - sp.cos(omega + psi1))).rewrite(sp.cos)))
print('  E·N2 - rho[sin psi2 + sin(omega+psi1)]   =', sp.simplify(E.dot(N2) - rho * (sp.sin(psi2) + sp.sin(omega + psi1))))
print('  e_x1·e_f2 - cos(omega) =', sp.simplify(ex1.dot(ef2) - sp.cos(omega)), ' e_z1·e_f2 - sin(omega) =', sp.simplify(ez1.dot(ef2) - sp.sin(omega)))
print('  e_x1·N2 + sin(omega)   =', sp.simplify(ex1.dot(N2) + sp.sin(omega)), ' e_z1·N2 - cos(omega) =', sp.simplify(ez1.dot(N2) - sp.cos(omega)))
# 基準面 (法線 N_ref = (-sin phi, cos phi)、接線 e_ref = (cos phi, sin phi)) との内積
Nref = sp.Matrix([-sp.sin(phi), sp.cos(phi)])
eref = sp.Matrix([sp.cos(phi), sp.sin(phi)])
print('  e_f2·e_ref - cos psi2 =', sp.simplify(ef2.dot(eref) - sp.cos(psi2)), '  e_f2·N_ref + sin psi2 =', sp.simplify(ef2.dot(Nref) + sp.sin(psi2)),
      '  N2·e_ref - sin psi2 =', sp.simplify(N2.dot(eref) - sp.sin(psi2)))
# (5.20)(5.21): 入口の直進
PX, PZ, X1, Pp, py_ = sp.symbols('P_X P_Z x1 p p_y', real=True)
t = X1 * sp.sin(psi1) / (PZ * sp.cos(psi1) - PX * sp.sin(psi1))
cond = sp.simplify((X1 + PX * t) * sp.sin(psi1) - PZ * t * sp.cos(psi1))
x2 = X1 / (sp.cos(psi1) - PX / PZ * sp.sin(psi1))
print('  面の条件の残差 (t を代入) =', cond, '   x2 - (x1 + P_X t)/cos... : 面上の x2 = (x1 + P_X t)cos psi1 + P_Z t sin psi1 の差 =',
      sp.simplify((X1 + PX * t) * sp.cos(psi1) + PZ * t * sp.sin(psi1) - x2))
print('  P_Z t - x2 sin psi1 =', sp.simplify(PZ * t - x2 * sp.sin(psi1)), ' (これが z の更新の p/P_Z x2 sin psi1 の根拠: 経路長 = p t)')
# Θ: 角度の足し算
th1, th2 = sp.symbols('theta1 theta2', real=True)
print('  Θ = (phi - psi2 - theta2) - (psi1 - theta1) = ', sp.simplify((phi - psi2 - th2) - (psi1 - th1)), ' = omega + theta1 - theta2')

print('\n== (d) 設計軌道のまわりの 1 次の展開 (§5.4): 入口・円弧・出口の各段の 2x2 行列 ==')
# 運動量を p_x = sin(theta) (p_z = cos(theta)、|p| = 1) で表し、θ で微分してから dp_x = cos(theta) dθ で戻す
x1s, th1s = sp.symbols('x1 theta1', real=True)
rb = rho
px2s = -x1s / rb * sp.sin(omega) + sp.cos(th1s) * sp.sin(omega) + sp.sin(th1s) * sp.cos(omega) - (rho / rb) * (sp.sin(psi2) + sp.sin(omega + psi1))
sqrt_term = sp.sqrt(1 - px2s**2)
x2s = x1s * sp.cos(omega) + rb * (sqrt_term - sp.cos(th1s) * sp.cos(omega) + sp.sin(th1s) * sp.sin(omega)) + rho * (sp.cos(omega + psi1) - sp.cos(psi2))
dsn = {x1s: 0, th1s: psi1}
px2d = sp.simplify(sp.expand_trig(px2s.subs(dsn)).subs(phi, omega + psi1 + psi2))
print('  設計点 (x1 = 0、p_x1 = sin psi1): p_x2 =', px2d, ' (= -sin psi2)')
sq = sp.cos(psi2)   # 設計点で sqrt(1 - px2^2) = cos(psi2)
# d sqrt(1 - px2^2) = -px2 d(px2)/sqrt = tan(psi2) d(px2) (設計点)
dx2_dx1 = sp.cos(omega) + rb * (-px2s / sq) * sp.diff(px2s, x1s)
dx2_dth = rb * ((-px2s / sq) * sp.diff(px2s, th1s) + sp.sin(th1s) * sp.cos(omega) + sp.cos(th1s) * sp.sin(omega))
def sub_(e):
    e = sp.expand_trig(e.subs(dsn)).subs(phi, omega + psi1 + psi2)
    return sp.simplify(e)
J_x = [sub_(dx2_dx1), sub_(dx2_dth / sp.cos(th1s))]
J_p = [sub_(sp.diff(px2s, x1s)), sub_(sp.diff(px2s, th1s) / sp.cos(th1s))]
print('  円弧の 2x2 [[∂x2/∂x1, ∂x2/∂p_x1], [∂p_x2/∂x1, ∂p_x2/∂p_x1]] (omega = phi - psi1 - psi2、rho_b = rho_0 = rho。ここで -px2/cos(psi2) -> tan(psi2) を使用):')
print('    ', [str(e) for e in J_x])
print('    ', [str(e) for e in J_p])
Jarc = sp.Matrix([J_x, J_p])
th2s, x2f = sp.symbols('theta2 x2f', real=True)
pxr = sp.sin(th2s + psi2)
xout = x2f * (sp.cos(psi2) + sp.tan(th2s + psi2) * sp.sin(psi2))
d2 = {x2f: 0, th2s: -psi2}
Jout = sp.Matrix([[sp.diff(xout, x2f).subs(d2), 0], [0, sp.diff(pxr, th2s).subs(d2) / sp.cos(th2s).subs(d2)]]).applyfunc(sp.simplify)
print('  出口の直進 (面の座標系 -> 基準面) の 2x2 =', Jout.tolist(), '  (設計点で p_xr = sin(theta2 + psi2) =', sp.simplify(pxr.subs(d2)), ')')
Jin = sp.Matrix([[1 / sp.cos(psi1), 0], [0, sp.cos(psi1)]])
print('  入口の直進 (基準面 -> 面の座標系) の 2x2 =', Jin.tolist())
Mtot = (Jout * Jarc * Jin).applyfunc(lambda e: sp.simplify(e))
Ex = lambda ps: sp.Matrix([[1, 0], [sp.tan(ps) / rho, 1]])
Sx = sp.Matrix([[sp.cos(phi), rho * sp.sin(phi)], [-sp.sin(phi) / rho, sp.cos(phi)]])
Mref = Ex(psi2) * Sx * Ex(psi1)
num = {psi1: 0.03, psi2: 0.07, phi: 0.1, rho: 10.0}
diff = (Mtot.subs(omega, phi - psi1 - psi2) - Mref)
print('  積 Jout·Jarc·Jin - E_x(psi2) S_x E_x(psi1) の数値 (psi1 = 0.03、psi2 = 0.07、phi = 0.1、rho = 10) =',
      [[float(diff[i, j].subs(num)) for j in range(2)] for i in range(2)])
Jarc_n = Jarc.subs(omega, phi - psi1 - psi2)
print('  数値: Jin =', [[round(float(Jin[i, j].subs(num)), 9) for j in range(2)] for i in range(2)],
      ' Jarc =', [[round(float(Jarc_n[i, j].subs(num)), 9) for j in range(2)] for i in range(2)],
      ' Jout =', [[round(float(Jout[i, j].subs(num)), 9) for j in range(2)] for i in range(2)])
print('  積の (2,1) 成分を symbolic に簡約 (omega を phi - psi1 - psi2 に戻す):', sp.simplify(Mtot[1, 0].subs(omega, phi - psi1 - psi2)))
print('  積の (1,1)・(1,2)・(2,2) と (5.25) の差を symbolic に:', [sp.simplify((Mtot - Mref.subs(phi, omega + psi1 + psi2))[i, j]) for (i, j) in ((0, 0), (0, 1), (1, 1))])

print('\n== (e) tbendcore の t2・t3・t2t3・px3・dpz32・x3・da (mpmath、40 桁) ==')
mp.mp.dps = 40
def tbendcore_rows(p, rhob, rho0, psi1, psi2, phi0, x2, px1, py2):
    """K64 tbendcore の 465-500 行を転記 (q の単位)。入力: p (全運動量)、x2 (F1・Maxwell 後の x)、px1 (面の座標系の q_x)、py2 (Maxwell 後の q_y)"""
    sinp1, cosp1, sinp2, cosp2 = mp.sin(psi1), mp.cos(psi1), mp.sin(psi2), mp.cos(psi2)
    w = phi0 - psi1 - psi2
    sinw, cosw = mp.sin(w), mp.cos(w)
    sqwh = 1 - cosw
    sinwp1 = mp.sin(phi0 - psi2)
    rhoe = rhob * p
    drhob = rhob - rho0
    dp = p - 1
    ph2 = mp.sqrt(1 - py2**2)
    dph2 = ph2 - 1
    pz2 = mp.sqrt(ph2**2 - px1**2)
    # pxi cosp1 + (dpzi - dph2) sinp1 = px1 - ph2 sinp1
    dpx2 = px1 - ph2 * sinp1
    drho1 = drhob + rhob * dp
    drho = drho1 + rhoe * dph2
    sinp2wp1 = sinp2 + sinwp1
    dcosp = cosp1 - cosp2
    t2 = (px1 + ph2 * sinp1) / (pz2 + ph2 * cosp1)
    dpx3a = (x2 * sinw - drho * sinp2wp1) / rhoe
    dpx3 = dpx3a - dpx2 * (cosw - sinw * t2)
    px3 = ph2 * sinp2 + dpx3
    sqwht2 = sqwh + sinw * t2
    px1px3 = ph2 * (sinp2 + sinp1) + dpx3a + dpx2 * sqwht2
    pz3 = mp.sqrt(ph2**2 - px3**2)
    dpz32 = px1px3 * (px1 - px3) / (pz2 + pz3)
    pz3ph2 = pz3 + ph2 * cosp2
    pz2px1 = pz2 * cosp1 + px1 * sinp1
    t3 = (px3 + ph2 * sinp2) / pz3ph2
    t2t3 = (ph2 * sinp2 + px1px3) / pz3ph2 + (ph2 * sinp1 + px1 * (dpz32 - ph2 * dcosp) / pz3ph2) / (pz2 + ph2 * cosp1)
    t4 = (cosp2 + t3 * sinp2) * pz2px1
    x3 = x2 * (cosw - rho0 / rhoe * t3 * sinw / ph2) + (rho0 * (cosw * t2t3 + sinw * (1 - t2 * t3)) * dpx2 - drho * (-sinp2wp1 * rho0 / rhoe * t3 - dpz32 - sqwh * pz2 - sinw * px1)) / ph2
    sda = (dpx2 * (sinp1 * (t2 * (pz3 * cosp2 + px3 * sinp2) - px1 * cosp2) - sinp2 * (t3 * pz2px1 - px3 * cosp1) + cosp1 * cosp2 * dpz32 + t4 * sqwht2) + dpx3a * t4) / ph2**2
    return dict(dict_=None, t2=t2, t3=t3, t2t3=t2t3, px3=px3, px1px3=px1px3, dpz32=dpz32, pz2=pz2, pz3=pz3, x3=x3, sda=sda, ph2=ph2, dpx2=dpx2, t4=t4, dpx3a=dpx3a, w=w)

cases = [dict(p=mp.mpf('1.03'), rhob=mp.mpf('10'), rho0=mp.mpf('10'), psi1=mp.mpf('0.03'), psi2=mp.mpf('0.07'), phi0=mp.mpf('0.1'), x2=mp.mpf('0.004'), px1=mp.mpf('0.04'), py2=mp.mpf('0.01')),
         dict(p=mp.mpf('0.98'), rhob=mp.mpf('9.5'), rho0=mp.mpf('10'), psi1=mp.mpf('-0.02'), psi2=mp.mpf('0.05'), phi0=mp.mpf('0.15'), x2=mp.mpf('-0.003'), px1=mp.mpf('-0.01'), py2=mp.mpf('-0.02'))]
for cs in cases:
    r = tbendcore_rows(**cs)
    p, rhob, rho0, psi1, psi2, phi0, x2, px1, py2 = [cs[k] for k in 'p rhob rho0 psi1 psi2 phi0 x2 px1 py2'.split()]
    ph2, pz2, pz3, px3, w = r['ph2'], r['pz2'], r['pz3'], r['px3'], r['w']
    rhoe = rhob * p
    th1 = mp.asin(px1 / ph2)
    th3 = mp.asin(px3 / ph2)
    # 幾何の式 (q の単位): (5.22) と (5.23)
    qx2 = -x2 * mp.sin(w) / rhoe + pz2 * mp.sin(w) + px1 * mp.cos(w) - (rho0 / rhoe) * (mp.sin(psi2) + mp.sin(w + psi1))
    xg = x2 * mp.cos(w) + rhoe * (mp.sqrt(ph2**2 - qx2**2) - pz2 * mp.cos(w) + px1 * mp.sin(w)) + rho0 * (mp.cos(w + psi1) - mp.cos(psi2))
    a1 = th1 - psi1
    a3 = th3 - psi2
    print('  ケース p = %s、psi1 = %s、psi2 = %s' % (mp.nstr(p, 4), mp.nstr(psi1, 3), mp.nstr(psi2, 3)))
    print('    px3 + q_x2 (5.22 を p で割った値) = %s     (px3 = -q_x2)' % mp.nstr(px3 + qx2, 3))
    print('    t2 - tan((theta1 + psi1)/2) = %s、t3 - tan((theta3 + psi2)/2) = %s   (px1 = ph2 sin theta1、px3 = ph2 sin theta3)' % (
        mp.nstr(r['t2'] - mp.tan((th1 + psi1) / 2), 3), mp.nstr(r['t3'] - mp.tan((th3 + psi2) / 2), 3)))
    print('    t2t3 - (t2 + t3) = %s   (t2t3 は積ではなく和 t2 + t3 を、桁落ちなく書いたもの)' % mp.nstr(r['t2t3'] - (r['t2'] + r['t3']), 3))
    A_, B_ = (th1 + psi1) / 2, (th3 + psi2) / 2
    print('    (t2 + t3) cos w + sin w (1 - t2 t3) - sin(A + B + w)/(cos A cos B) = %s   (A = (theta1 + psi1)/2、B = (theta3 + psi2)/2)' % mp.nstr((r['t2'] + r['t3']) * mp.cos(w) + mp.sin(w) * (1 - r['t2'] * r['t3']) - mp.sin(A_ + B_ + w) / (mp.cos(A_) * mp.cos(B_)), 3))
    print('    pz2 - (ph2 cos psi1 - dpx2 t2) = %s   (dpx2 = px1 - ph2 sin psi1)' % mp.nstr(pz2 - (ph2 * mp.cos(psi1) - r['dpx2'] * r['t2']), 3))
    print('    px1px3 - (px1 + px3) = %s、dpz32 - (pz3 - pz2) = %s' % (mp.nstr(r['px1px3'] - (px1 + px3), 3), mp.nstr(r['dpz32'] - (pz3 - pz2), 3)))
    print('    t4 - ph2 cos(a1)(cos psi2 + t3 sin psi2) = %s   (pz2px1 = ph2 cos a1、a1 = theta1 - psi1)' % mp.nstr(r['t4'] - ph2 * mp.cos(a1) * (mp.cos(psi2) + r['t3'] * mp.sin(psi2)), 3))
    print('    x3 - (5.23) = %s' % mp.nstr(r['x3'] - xg, 3))
    print('    sin(da) の引数 - sin(a1 + a3) = %s   (da = Theta - phi0 = a1 + a3、a1 = theta1 - psi1、a3 = theta3 - psi2)' % mp.nstr(r['sda'] - mp.sin(a1 + a3), 3))
    print('    Theta = theta1 + theta3 + omega: asin(px1/ph2) - asin(q_x2/ph2) + omega - (phi0 + da) = %s' % mp.nstr(th1 - mp.asin(qx2 / ph2) + w - (phi0 + mp.asin(r['sda'])), 3))
    # x3 の分解: rhoe K (K = dpz32 + sqwh pz2 + sinw px1)
    K_ = r['dpz32'] + (1 - mp.cos(w)) * pz2 + mp.sin(w) * px1
    print('    pz3 - pz2 cos(w) + px1 sin(w) - K = %s  (K = dpz32 + sqwh pz2 + sinw px1)' % mp.nstr(pz3 - pz2 * mp.cos(w) + px1 * mp.sin(w) - K_, 3))

print('\n== (f) 端の磁場 (5.30)(5.31)(5.32) (sympy) ==')
s, yy, rb_ = sp.symbols('s y rho_b', real=True)
g = sp.Function('g')
V = (yy * g(s) - yy**3 / 6 * g(s).diff(s, 2) + yy**5 / 120 * g(s).diff(s, 4)) / rb_
lap = sp.simplify(sp.diff(V, yy, 2) + sp.diff(V, s, 2))
print('  ∇²V (y と s だけ) =', sp.expand(lap), '  (残りは y^5 g^(6) の 1 項だけ。V の打ち切りの次数)')
by = sp.diff(V, yy)
bs = sp.diff(V, s)
print('  b_y =', sp.expand(by))
print('  b_s =', sp.expand(bs))
print('  ∇·b = ∂_y b_y + ∂_s b_s =', sp.expand(sp.diff(by, yy) + sp.diff(bs, s)))
bx_, by_, bs_ = sp.symbols('b_x b_y b_s')
Pv = sp.Matrix([sp.Symbol('p_x'), sp.Symbol('p_y'), sp.Symbol('p_z')])
bv = sp.Matrix([0, by_, bs_])
print('  P × b =', list(Pv.cross(bv)), '  -> dp_x/ds = (p_y b_s - p_z b_y)/p_z、dp_y/ds = -p_x b_s/p_z')

print('\n== (g) 硬い端の Maxwell フリンジ (5.33)〜(5.36) ==')
pxs, ps, kap, ysym, rr = sp.symbols('p_x p kappa y rho', positive=True)
pz = sp.sqrt(ps**2 - pxs**2)
F = pxs / pz
print('  dx/ds = p_x/p_z の p_x に、スパイク kappa g\' (kappa = y^2/(2 rho)) を足す。∫ g\' ds = 1 の 1 次: Delta x = kappa dF/dp_x =',
      sp.simplify(kap * sp.diff(F, pxs)), '  (p^2/p_z^3 に一致:', sp.simplify(sp.diff(F, pxs) - ps**2 / pz**3), ')')
print('  dz/ds = -p/p_z: Delta z = -p kappa d(1/p_z)/dp_x =', sp.simplify(-ps * kap * sp.diff(1 / pz, pxs)), '  (-(y^2/(2rho)) p_x p/p_z^3 =',
      sp.simplify(-ps * kap * sp.diff(1 / pz, pxs) + kap * pxs * ps / pz**3), ' の差)')
fnew = -ysym**2 / (2 * rr) * pxs / pz
print('  f = -(y^2/(2 rho)) p_x/p_z (C 付録の約束)。 -∂f/∂p_x =', sp.simplify(-sp.diff(fnew, pxs)), '、 +∂f/∂y =', sp.simplify(sp.diff(fnew, ysym)), '、 -∂f/∂p =', sp.simplify(-sp.diff(fnew, ps)))
# C 付録のポアソン括弧での exp(:f:): 変数 (x, px), (y, py), (z, p) (p は 1+delta)
xx, pyy, zz = sp.symbols('x p_y z')
def PB(f1, f2):
    pairs = [(xx, pxs), (ysym, pyy), (zz, ps)]
    return sum(sp.diff(f1, q) * sp.diff(f2, pp_) - sp.diff(f1, pp_) * sp.diff(f2, q) for q, pp_ in pairs)
print('  C 付録の括弧で :f: x =', sp.simplify(PB(fnew, xx)), '、:f: p_y =', sp.simplify(PB(fnew, pyy)), '、:f: z =', sp.simplify(PB(fnew, zz)))
print('  :f:^2 x =', sp.simplify(PB(fnew, PB(fnew, xx))), '(級数は 1 項で止まる)')
# y^4 の係数 a, b と f の関係
a_ = 1 - ysym**2 / (6 * rr**2)
b_ = sp.Rational(1, 2) - ysym**2 / (24 * rr**2)
f4 = -(1 / (2 * rr)) * (ysym**2 - ysym**4 / (12 * rr**2)) * pxs / pz
print('  f_4 = -(1/(2 rho))(y^2 - y^4/(12 rho^2)) p_x/p_z:  ∂f/∂y + a y/rho p_x/p_z =', sp.simplify(sp.diff(f4, ysym) + a_ * ysym / rr * pxs / pz),
      '、 -∂f/∂p_x - b y^2/rho p^2/p_z^3 =', sp.simplify(-sp.diff(f4, pxs) - b_ * ysym**2 / rr * ps**2 / pz**3))

print('\n== (h) F1 の積分 (6.3) ==')
ell = sp.Symbol('ell', positive=True)
sg = sp.Symbol('s', real=True)
gl = sp.Piecewise((0, sg < -ell / 2), (sg / ell + sp.Rational(1, 2), sg < ell / 2), (1, True))
theta = sp.Piecewise((0, sg < 0), (1, True))
Gl_int = sp.integrate(sg / ell + sp.Rational(1, 2) - sp.Rational(1, 2) * 0, (sg, -ell / 2, ell / 2))
print('  直線の立ち上がり: ∫g(1-g) ds =', sp.simplify(sp.integrate((sg / ell + sp.Rational(1, 2)) * (sp.Rational(1, 2) - sg / ell), (sg, -ell / 2, ell / 2))), '(= ell/6 -> F1 = ell)')
print('                   ∫g\'^2 ds =', sp.simplify(sp.integrate(1 / ell**2, (sg, -ell / 2, ell / 2))))
I0 = sp.integrate(sg / ell + sp.Rational(1, 2), (sg, -ell / 2, 0)) + sp.integrate(sg / ell - sp.Rational(1, 2), (sg, 0, ell / 2))
I1 = sp.integrate(sg * (sg / ell + sp.Rational(1, 2)), (sg, -ell / 2, 0)) + sp.integrate(sg * (sg / ell - sp.Rational(1, 2)), (sg, 0, ell / 2))
print('                   I_0 = ∫(g - Θ) ds =', sp.simplify(I0), '、 I_1 = ∫ s (g - Θ) ds =', sp.simplify(I1), '(= -ell^2/24)')
dd = mp.mpf('0.002')
gt = lambda u: (1 + mp.tanh(u / dd)) / 2
print('  tanh (d = 2 mm): F1 = 6∫g(1-g) =', mp.nstr(6 * mp.quad(lambda u: gt(u) * (1 - gt(u)), [-mp.inf, 0, mp.inf]), 12), ' (3d = 0.006)')
gp = lambda u: mp.diff(gt, u)
print('                   ∫g\'^2 ds · F1 =', mp.nstr(mp.quad(lambda u: (1 / (4 * dd**2 * mp.cosh(u / dd)**4)), [-mp.inf, 0, mp.inf]) * 3 * dd, 12))
I1t = 2 * mp.quad(lambda u: u * (gt(u) - 1), [0, mp.inf])
print('                   I_1 = ∫ s (g - Θ) ds =', mp.nstr(I1t, 12), '  -pi^2 d^2/24 =', mp.nstr(-mp.pi**2 * dd**2 / 24, 12), '  比 I_1/(-F1^2/24) =', mp.nstr(I1t / (-(3 * dd)**2 / 24), 10), ' pi^2/9 =', mp.nstr(mp.pi**2 / 9, 10))
# 誤差関数の端
sig = mp.mpf('0.002')
ge = lambda u: (1 + mp.erf(u / (sig * mp.sqrt(2)))) / 2
F1e = 6 * mp.quad(lambda u: ge(u) * (1 - ge(u)), [-mp.inf, 0, mp.inf])
gpe = lambda u: mp.exp(-u**2 / (2 * sig**2)) / (sig * mp.sqrt(2 * mp.pi))
print('  誤差関数 (sigma = 2 mm): F1 =', mp.nstr(F1e, 10), ' ∫g\'^2 F1 =', mp.nstr(mp.quad(lambda u: gpe(u)**2, [-mp.inf, 0, mp.inf]) * F1e, 8), ' (3/pi 型の 0.955)')
# 部分積分
print('  部分積分: ∫G g\'\'\' ds = [G g\'\'] - ∫ g g\'\' ds = -[g g\'] + ∫ g\'^2 ds = ∫ g\'^2 ds (境界項は g\', g\'\' -> 0 で消える)。tanh で数値確認:')
Gt = lambda u: dd * (mp.log(2 * mp.cosh(u / dd)) / 2 + u / (2 * dd))  # ∫ g
g2 = lambda u: mp.diff(gt, u, 2)
g3 = lambda u: mp.diff(gt, u, 3)
lhs = mp.quad(lambda u: Gt(u) * g3(u), [-20 * dd, 0, 20 * dd])
rhs = mp.quad(lambda u: gp(u)**2, [-20 * dd, 0, 20 * dd])
print('    ∫G g\'\'\' ds = %s、 ∫g\'^2 ds = %s' % (mp.nstr(lhs, 10), mp.nstr(rhs, 10)))
# (5.42): 符号を反転した f_lin
pp, yv, pxv, F1s, rhos = sp.symbols('p y p_x F1 rho', positive=True)
dxfr = F1s**2 / (24 * rhos)
dyfr = F1s / (6 * rhos**2)
dyfra = 4 * dyfr / F1s**2
flin = -pxv * dxfr * (1 - 1 / pp) + (dyfr * yv**2 / 2 - dyfra * yv**4 / 4) / pp
dl = pp - 1
print('  f_lin (C 付録の約束): -∂f/∂p_x =', sp.simplify(-sp.diff(flin, pxv)), '= dxfr δ/p:', sp.simplify(-sp.diff(flin, pxv) - dxfr * dl / pp))
print('                        +∂f/∂y  =', sp.simplify(sp.diff(flin, yv)), ' vs (dyfr - dyfra y^2) y / p:', sp.simplify(sp.diff(flin, yv) - (dyfr - dyfra * yv**2) * yv / pp))
print('                        -∂f/∂p  =', sp.simplify(-sp.diff(flin, pp)), ' vs [dxfr p_x + (dyfr/2 - dyfra y^2/4) y^2/p]/p:',
      sp.simplify(-sp.diff(flin, pp) - (dxfr * pxv + (dyfr / 2 - dyfra * yv**2 / 4) * yv**2 / pp) / pp))
Vhelp = -(F1s**2 / (24 * rhos) * pxv + F1s / (12 * rhos**2) * yv**2 - yv**4 / (6 * F1s * rhos**2)) / pp
print('  ヘルプ (平文版) の V - f_old (f_old = -f_lin) =', sp.simplify(Vhelp - (-flin)), '  (= -dxfr p_x: p に依らない x のずれの項。L\' で吸収される)')

print('\n== (i) L\' の補正 (5.43): 面の変位 ΔX = -F1^2/(24 rho) からの幾何 ==')
s0, rho_, phi_, ps1, ps2 = sp.symbols('s0 rho phi psi1 psi2', positive=True)
Ev = sp.Matrix([-rho_ * (1 - sp.cos(phi_)), rho_ * sp.sin(phi_)])
ex1v = sp.Matrix([sp.cos(ps1), sp.sin(ps1)])
ef2v = sp.Matrix([sp.cos(phi_ - ps2), sp.sin(phi_ - ps2)])
Ehat = sp.Matrix([-sp.sin(phi_ / 2), sp.cos(phi_ / 2)])
print('  E の方向 E/|E| - (-sin(phi/2), cos(phi/2)) =', sp.simplify(Ev / (2 * rho_ * sp.sin(phi_ / 2)) - Ehat).T)
# 実際の内側の円は、入口の面の点 -s0 e_x1 と、出口の面の点 E - s0 e_f2 を結ぶ (硬い端の円弧を、各端で DX = -s0 (s0 = F1^2/(24 rho)) だけずらしたもの)
chord_vec = Ev - s0 * (ef2v - ex1v)
dchord = sp.simplify((chord_vec.dot(Ehat) - 2 * rho_ * sp.sin(phi_ / 2)))
print('  弦の長さの変化 (s0 の 1 次) =', sp.simplify(sp.expand_trig(dchord)))
om2 = (phi_ - ps1 - ps2) / 2
print('  (ef2 - ex1)·E_hat - 2 sin(omega/2) cos((psi1 - psi2)/2) =', sp.simplify(sp.expand_trig(((ef2v - ex1v).dot(Ehat)) - 2 * sp.sin(om2) * sp.cos((ps1 - ps2) / 2))))
# L = rho phi, 弦 = 2 rho sin(phi/2): 角 phi を保ったまま弦が dchord だけ長くなる -> 半径 rho の変化 -> L' の変化
Lp = sp.Symbol('Lp', positive=True)
# 弦 C = 2 rho sin(phi/2) + dC を保つ: rho' = rho + dC/(2 sin(phi/2)),  L' = rho' phi
dC = sp.Symbol('dC')
dLp = sp.simplify(dC * phi_ / (2 * sp.sin(phi_ / 2)))
print('  L\' の変化 = dC phi/(2 sin(phi/2)) =', dLp)
# dC = -s0 [sin(phi/2 - psi1) + sin(phi/2 - psi2)]。対称 psi1 = psi2 では dC = -2 s0 sin(omega/2) -> dL' = -s0 phi sin(omega/2)/sin(phi/2)
dL_sym = sp.simplify(dLp.subs(dC, -2 * s0 * sp.sin(om2)))
F1sym, Lsym = sp.symbols('F1 L', positive=True)
dL_expr = dL_sym.subs(s0, F1sym**2 / (24 * rho_)).subs(rho_, Lsym / phi_)
target = -(phi_ * F1sym)**2 / (24 * Lsym) * sp.sin(om2) / sp.sin(phi_ / 2)
print("  dC = -2 s0 sin(omega/2) とすると dL' - [-(phi F1)^2/(24 L) sin(omega/2)/sin(phi/2)] =", sp.simplify(dL_expr - target), "  (K64・ヘルプの式の補正 L' - L に一致)")

print('\n== (j) rbh・rbl: 階段近似の 2 点ガウス・ルジャンドル ==')
rbh = mp.mpf(1) / 2 + 1 / mp.sqrt(12)
rbl = 1 / (6 * rbh)
print('  rbh = %s、rbl = %s、rbl - (1/2 - 1/sqrt(12)) = %s、rbh + rbl = %s、rbh rbl = %s' % (mp.nstr(rbh, 12), mp.nstr(rbl, 12), mp.nstr(rbl - (mp.mpf(1) / 2 - 1 / mp.sqrt(12)), 3), mp.nstr(rbh + rbl, 12), mp.nstr(rbh * rbl, 12)))
print('  線形の立ち上がり g(u) = u (0 <= u <= 1) の ∫ g^k du = 1/(k+1) と、等しい幅 1/2 の 2 段の階段 (高さ rbl、rbh) の (1/2)(rbl^k + rbh^k):')
for k in range(0, 5):
    print('    k = %d: 厳密 %s、階段 %s、差 %s' % (k, mp.nstr(mp.mpf(1) / (k + 1), 12), mp.nstr((rbl**k + rbh**k) / 2, 12), mp.nstr((rbl**k + rbh**k) / 2 - mp.mpf(1) / (k + 1), 3)))

print('\n== (k) tbendicorr の各行と K1 の非線形の蹴り (sympy) ==')
u = sp.Symbol('u')
ser = sp.series((1 + u) * sp.log(1 + u) - u, u, 0, 6).removeO()
code = u**2 * (sp.Rational(1, 2) - u * (2 - u) / 12)
print('  (1+u) ln(1+u) - u の級数 =', sp.expand(ser))
print('  K64 の x_r^2 (1/2 - x_r (2 - x_r)/12) =', sp.expand(code), ' 差 (u^4 まで) =', sp.expand(ser - code))
# 運動学の残り
qxs, qys, hs, xs, ell_ = sp.symbols('q_x q_y h x ell')
pzq = sp.sqrt(1 - qxs**2 - qys**2)
zrest = qxs**2 / 2 + qys**2 / 2 + (1 + hs * xs) * (1 - 1 / pzq)
zexact = 1 - (1 + hs * xs) / pzq
zh2 = -hs * xs - (qxs**2 + qys**2) / 2
print('  z\' の残り = [1 - (1+hx)/pz] - [-hx - q^2/2] と q^2/2 + (1+hx)(1 - 1/pz) の差 =', sp.simplify(zexact - zh2 - zrest))
print('  (1 - 1/pz) = dpzi/pzi (dpzi = pz - 1):', sp.simplify((1 - 1 / pzq) - (pzq - 1) / pzq))
a_x = -qxs * (pzq - 1) / pzq * 0  # 記号のみ
dpz, pzi_ = sp.symbols('dpzi pzi')
a_c = -qxs * dpz / pzi_   # x' の定数項: q_x (1/pz - 1) = -q_x dpzi/pzi
b_c = hs * qxs / pzi_    # x の係数
x_exact = xs * sp.exp(b_c * ell_) + a_c * (sp.exp(b_c * ell_) - 1) / b_c
x_code = (xs + a_c * ell_) / (1 - b_c * ell_)
print('  x\' = a + b x (a = -q_x dpzi/pzi、b = h q_x/pzi) の厳密解と tbendicorr の (x + a l)/(1 - b l) の級数の差:',
      sp.simplify(sp.series(x_exact - x_code, ell_, 0, 3).removeO()))
print('  tbendicorr の x = (-al dpzi pxi + pzi x)/(pzi - phi pxi) は上の式 (a l = -l dpzi q_x/pzi、b l = phi q_x/pzi の分母分子に pzi をかけたもの):',
      sp.simplify(x_code - (-ell_ * dpz * qxs + pzi_ * xs) / (pzi_ - hs * ell_ * qxs)))

print('\n== (l) 分散 R16・R26・R56 (sympy、§4.1) ==')
dl_, L_, rho_l, gam = sp.symbols('delta L rho gamma0', positive=True)
hh_ = 1 / rho_l
pl = 1 + dl_
xe_l = (hh_ * pl - hh_) / (hh_ * hh_)
w2_l = hh_ * hh_ / pl
wl = sp.sqrt(w2_l)
xL = xe_l * (1 - sp.cos(wl * L_))
qL = xe_l * wl * sp.sin(wl * L_)
print('  x_e =', sp.simplify(xe_l), ' (rho delta)')
print('  R16 = ∂x/∂δ|0 =', sp.simplify(sp.diff(xL, dl_).subs(dl_, 0)), '、 R26 = ∂q/∂δ|0 =', sp.simplify(sp.diff(qL, dl_).subs(dl_, 0)))
# z: z' = -dv - h x - q^2/2、dv = -delta/gamma0^2 (1 次)
zz1 = sp.integrate(sp.expand((-(-dl_ / gam**2) - hh_ * xe_l * (1 - sp.cos(sp.Symbol('t') / rho_l))).subs(sp.Symbol('t'), sp.Symbol('t'))), (sp.Symbol('t'), 0, L_))
print('  R56 = ∂z/∂δ|0 (1 次) =', sp.simplify(sp.diff(zz1, dl_).subs(dl_, 0)), ' (= -(L - rho sin(L/rho)) + L/gamma0^2)')
