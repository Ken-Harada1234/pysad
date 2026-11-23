# -*- coding: utf-8 -*-
r"""BEND の章 (05_bend.md) で、紙の上で 1 段ずつ追う導出の各行を sympy・mpmath (40 桁) で確かめる (pysad は使わない)。

(a) §5.2: 出口の面から基準面への直進 (内積・時間 t・x_out)
(b) §5.4: 各段の 1 次の微分と通分、E_x(psi2) S_x E_x(psi1) の 4 成分、J_in・J_out が対角になる理由
(c) §7.1: tbendcore の t2t3・x3・da の導出の中間式 (d3 - d2、t2t3 = t2 + t3、加法定理、(5.23) との一致、da の引数 = sin(Δ1 + Δ3))
(d) §8: K1 の b_y と力、xi の厳密解との差
(e) §6.2: p_z → p_v の近似の大きさ
"""
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)
mp.mp.dps = 40

print('== (a) 出口の面から基準面への直進 (§5.2) ==')
phi, ps2, x2, pxf, pzf, py, pp = sp.symbols('phi psi2 x2 p_x2 p_z2 p_y p', real=True)
ef2 = sp.Matrix([sp.cos(phi - ps2), sp.sin(phi - ps2)])
N2 = sp.Matrix([-sp.sin(phi - ps2), sp.cos(phi - ps2)])
eref = sp.Matrix([sp.cos(phi), sp.sin(phi)])
Nref = sp.Matrix([-sp.sin(phi), sp.cos(phi)])
chk = lambda e: sp.simplify(sp.expand_trig(e))
print('  e_ref·e_f2 - cos ψ2 =', chk(eref.dot(ef2) - sp.cos(ps2)), '、e_ref·N2 - sin ψ2 =', chk(eref.dot(N2) - sp.sin(ps2)))
print('  N_ref·e_f2 + sin ψ2 =', chk(Nref.dot(ef2) + sp.sin(ps2)), '、N_ref·N2 - cos ψ2 =', chk(Nref.dot(N2) - sp.cos(ps2)))
P = pxf * ef2 + pzf * N2
pxr = chk(P.dot(eref))
pzr = chk(P.dot(Nref))
print('  p_xr - (p_x2 cos ψ2 + p_z2 sin ψ2) =', chk(pxr - (pxf * sp.cos(ps2) + pzf * sp.sin(ps2))))
print('  p_zr - (-p_x2 sin ψ2 + p_z2 cos ψ2) =', chk(pzr - (-pxf * sp.sin(ps2) + pzf * sp.cos(ps2))))
pzr_s, pxr_s = sp.symbols('p_zr p_xr')
t = x2 * sp.sin(ps2) / pzr_s       # (r + t P - E)·N_ref = 0: -x2 sin ψ2 + t p_zr = 0
r_out = x2 * (eref.dot(ef2)) + t * pxr_s
print('  x_out - x2 (cos ψ2 + (p_xr/p_zr) sin ψ2) =', sp.simplify(r_out.subs({}) - x2 * (sp.cos(ps2) + pxr_s / pzr_s * sp.sin(ps2)) - 0).subs(sp.cos(ps2), sp.cos(ps2)) if False else sp.simplify(x2 * sp.cos(ps2) + t * pxr_s - x2 * (sp.cos(ps2) + pxr_s / pzr_s * sp.sin(ps2))))
# 数値で、出口の基準面に実際に当たること
vals = {phi: sp.Rational(1, 10), ps2: sp.Rational(7, 100), x2: sp.Rational(3, 1000), pxf: sp.Rational(-7, 100), pzf: sp.sqrt(1 - sp.Rational(49, 10000))}
rr = sp.Matrix([x2, 0]).subs(vals)
pt = (ef2 * x2 + N2 * 0).subs(vals)
tt = (x2 * sp.sin(ps2) / pzr).subs(vals)
hit = (pt + tt * P.subs(vals))
print('  面の点 + t P の基準面からの距離 (N_ref·(r - 0)) = 0 の確認 (E を原点とした座標):', sp.N(Nref.subs(vals).dot(hit), 20))

print('\n== (b) 各段の 1 次の微分と通分 (§5.4) ==')
ps1, om = sp.symbols('psi1 omega', real=True)
rho = sp.Symbol('rho', positive=True)
phi_ = ps1 + ps2 + om
t1_, t2_ = sp.tan(ps1), sp.tan(ps2)
dpx2_dpx1 = -sp.sin(om) * t1_ + sp.cos(om)
print('  ∂p_x2/∂p_x1 - cos(ω+ψ1)/cosψ1 =', chk(dpx2_dpx1 - sp.cos(om + ps1) / sp.cos(ps1)))
dx2_dx1 = sp.cos(om) + rho * t2_ * (-sp.sin(om) / rho)
print('  ∂x2/∂x1 - cos(ω+ψ2)/cosψ2 =', chk(dx2_dx1 - sp.cos(om + ps2) / sp.cos(ps2)), '  (ω + ψ2 = φ - ψ1)')
dx2_dpx1 = rho * (t2_ * dpx2_dpx1 + t1_ * sp.cos(om) + sp.sin(om))
num = sp.sin(ps2) * sp.cos(om + ps1) + sp.cos(ps2) * (sp.sin(ps1) * sp.cos(om) + sp.cos(ps1) * sp.sin(om))
print('  通分後の分子 sinψ2 cos(ω+ψ1) + cosψ2 sin(ω+ψ1) - sin(ω+ψ1+ψ2) =', chk(num - sp.sin(om + ps1 + ps2)))
print('  ∂x2/∂p_x1 - ρ sin φ/(cosψ1 cosψ2) =', chk(dx2_dpx1 - rho * sp.sin(om + ps1 + ps2) / (sp.cos(ps1) * sp.cos(ps2))))
Jin = sp.diag(1 / sp.cos(ps1), sp.cos(ps1))
Jout = sp.diag(sp.cos(ps2), 1 / sp.cos(ps2))
Jarc = sp.Matrix([[dx2_dx1, dx2_dpx1], [-sp.sin(om) / rho, dpx2_dpx1]])
M = (Jout * Jarc * Jin)
c, s_ = sp.cos(phi_), sp.sin(phi_)
Ex = lambda ps: sp.Matrix([[1, 0], [sp.tan(ps) / rho, 1]])
Sx = sp.Matrix([[c, rho * s_], [-s_ / rho, c]])
SE1 = Sx * Ex(ps1)
print('  S_x E_x(ψ1) =', [sp.simplify(e) for e in SE1])
full = Ex(ps2) * SE1
print('  E_x(ψ2) S_x E_x(ψ1) の 4 成分 - J_out J_arc J_in =', [chk(e) for e in (full - M)])
m21 = (sp.tan(ps2) * (c + s_ * t1_) - s_ + c * t1_) * 1 / rho
print('  M21 の括弧: c (t1 + t2) - s (1 - t1 t2) - (-sin ω/(cosψ1 cosψ2)) =', chk(c * (t1_ + t2_) - s_ * (1 - t1_ * t2_) + sp.sin(om) / (sp.cos(ps1) * sp.cos(ps2))))
print('  tanψ1 + tanψ2 - sin(ψ1+ψ2)/(cosψ1 cosψ2) =', chk(t1_ + t2_ - sp.sin(ps1 + ps2) / (sp.cos(ps1) * sp.cos(ps2))), '、1 - tanψ1 tanψ2 - cos(ψ1+ψ2)/(cosψ1 cosψ2) =', chk(1 - t1_ * t2_ - sp.cos(ps1 + ps2) / (sp.cos(ps1) * sp.cos(ps2))))
# J_in が対角になる理由: x1 = x/(cosψ1 - (p_x/p_z) sinψ1) で x = 0 なら ∂x1/∂p_x = x (..) = 0
xs, pxs = sp.symbols('x p_x', real=True)
x1e = xs / (sp.cos(ps1) - pxs / sp.sqrt(1 - pxs**2) * sp.sin(ps1))
print('  ∂x1/∂p_x |_{x=0} =', sp.simplify(sp.diff(x1e, pxs).subs(xs, 0)), '、∂x1/∂x |_{x=0, p_x=0} =', sp.simplify(sp.diff(x1e, xs).subs({xs: 0, pxs: 0})), '(= 1/cosψ1)')

print('\n== (c) tbendcore の t2t3・x3・da の導出の中間式 (§7.1。40 桁) ==')
mp.mp.dps = 40
mpf = mp.mpf


def build(seed, psi1, psi2, phi0, h, th1_off, x2v, rho0, rhoe):
    """(5.22) から出口の運動量を作り、すべての量を返す"""
    sinp1, cosp1, sinp2, cosp2 = mp.sin(psi1), mp.cos(psi1), mp.sin(psi2), mp.cos(psi2)
    omv = phi0 - psi1 - psi2
    cosw, sinw, sqwh = mp.cos(omv), mp.sin(omv), 1 - mp.cos(omv)
    th1 = psi1 + th1_off
    qx1, pz2 = h * mp.sin(th1), h * mp.cos(th1)
    drho = rhoe * h - rho0
    sinp2wp1 = sinp2 + mp.sin(omv + psi1)
    dpx2 = qx1 - h * sinp1
    t2 = (qx1 + h * sinp1) / (pz2 + h * cosp1)
    dpx3a = (x2v * sinw - drho * sinp2wp1) / rhoe
    dpx3 = dpx3a - dpx2 * (cosw - sinw * t2)
    px3 = h * sinp2 + dpx3
    pz3 = mp.sqrt(h * h - px3 * px3)
    th3 = mp.asin(px3 / h)
    dcosp = cosp1 - cosp2
    dpz32 = (pz3 - pz2)
    return dict(psi1=psi1, psi2=psi2, phi0=phi0, om=omv, h=h, th1=th1, th3=th3, qx1=qx1, pz2=pz2, px3=px3, pz3=pz3, drho=drho,
                rho0=rho0, rhoe=rhoe, x2=x2v, dpx2=dpx2, dpx3a=dpx3a, dpx3=dpx3, t2=t2, dcosp=dcosp, dpz32=dpz32,
                sinp1=sinp1, cosp1=cosp1, sinp2=sinp2, cosp2=cosp2, cosw=cosw, sinw=sinw, sqwh=sqwh, sinp2wp1=sinp2wp1)


worst = {}


def upd(key, v):
    worst[key] = max(worst.get(key, mpf(0)), abs(v))


cases = [(mpf('0.03'), mpf('0.07'), mpf('0.1'), mpf('1.02'), mpf('0.004'), mpf('0.0031'), mpf(10), mpf('10.2')),
         (mpf('0.3'), mpf('-0.1'), mpf('0.1'), mpf('0.99'), mpf('-0.006'), mpf('-0.002'), mpf(10), mpf('9.9')),
         (mpf('-0.05'), mpf('0.25'), mpf('0.12'), mpf('1.0'), mpf('0.0'), mpf('0.001'), mpf(12), mpf('12.0'))]
for seed, args in enumerate(cases):
    v = build(seed, *args)
    h, om_, th1, th3 = v['h'], v['om'], v['th1'], v['th3']
    psi1, psi2 = v['psi1'], v['psi2']
    A, B = (th1 + psi1) / 2, (th3 + psi2) / 2
    a1, a3 = (th1 - psi1) / 2, (th3 - psi2) / 2
    D1_, D3_ = th1 - psi1, th3 - psi2
    t2, t3 = v['t2'], (v['px3'] + h * v['sinp2']) / (v['pz3'] + h * v['cosp2'])
    d2 = v['pz2'] + h * v['cosp1']
    d3 = v['pz3'] + h * v['cosp2']
    upd('t2 = tan A', t2 - mp.tan(A))
    upd('t3 = tan B', t3 - mp.tan(B))
    upd('D3 - D2 = dpz32 - h dcosp', (d3 - d2) - (v['dpz32'] - h * v['dcosp']))
    px1px3 = v['qx1'] + v['px3']
    upd('px1px3 (Fortran) = qx1 + px3', h * (v['sinp2'] + v['sinp1']) + v['dpx3a'] + v['dpx2'] * (v['sqwh'] + v['sinw'] * t2) - px1px3)
    t2t3 = ((h * v['sinp2'] + px1px3) / d3 + (h * v['sinp1'] + v['qx1'] * (v['dpz32'] - h * v['dcosp']) / d3) / d2)
    upd('t2t3 (Fortran) = t2 + t3', t2t3 - (t2 + t3))
    upd('cosω (t2+t3) + sinω (1 - t2 t3) = sin((Θ+φ)/2)/(cosA cosB)',
        v['cosw'] * (t2 + t3) + v['sinw'] * (1 - t2 * t3) - mp.sin(A + B + om_) / (mp.cos(A) * mp.cos(B)))
    upd('A + B + ω = (Θ + φ0)/2', (A + B + om_) - ((th1 + th3 + om_) + v['phi0']) / 2)
    upd('dpx2 = 2 h cosA sin a1', v['dpx2'] - 2 * h * mp.cos(A) * mp.sin(a1))
    upd('dpx3a = 2h [cosB sin a3 + sin a1 cos(ω+A)]', v['dpx3a'] - 2 * h * (mp.cos(B) * mp.sin(a3) + mp.sin(a1) * mp.cos(om_ + A)))
    # x3
    Bq = v['dpz32'] + v['sqwh'] * v['pz2'] + v['sinw'] * v['qx1']
    Kq = mp.cos(om_ + psi1) - mp.cos(psi2)
    upd('B = h [cosθ3 - cos(θ1+ω)]', Bq - h * (mp.cos(th3) - mp.cos(th1 + om_)))
    lhs = (v['cosw'] * (t2 + t3) + v['sinw'] * (1 - t2 * t3)) * v['dpx2'] - t3 * v['dpx3a']
    upd('(I1) (cosω T + sinω (1 - t2 t3)) dpx2 - t3 dpx3a = B + h K', lhs - (Bq + h * Kq))
    upd('(I1) の左辺/h = 2 sin a1 sin(ω+A) - 2 sinB sin a3', lhs / h - (2 * mp.sin(a1) * mp.sin(om_ + A) - 2 * mp.sin(B) * mp.sin(a3)))
    sinp2wp1 = v['sinp2wp1']
    rho0, rhoe, drho, x2v = v['rho0'], v['rhoe'], v['drho'], v['x2']
    t2t3_ = t2 + t3
    x3code = (x2v * (v['cosw'] - rho0 / rhoe * t3 * v['sinw'] / h)
              + (rho0 * (v['cosw'] * t2t3_ + v['sinw'] * (1 - t2 * t3)) * v['dpx2']
                 - drho * (-sinp2wp1 * rho0 / rhoe * t3 - v['dpz32'] - v['sqwh'] * v['pz2'] - v['sinw'] * v['qx1'])) / h)
    x3_523 = x2v * v['cosw'] + rhoe * (v['pz3'] - v['pz2'] * v['cosw'] + v['qx1'] * v['sinw']) + rho0 * Kq
    upd('x3 (Fortran) = (5.23)', x3code - x3_523)
    # da
    sinp1, cosp1, sinp2, cosp2 = v['sinp1'], v['cosp1'], v['sinp2'], v['cosp2']
    pz2px1 = v['pz2'] * cosp1 + v['qx1'] * sinp1
    t4 = (cosp2 + t3 * sinp2) * pz2px1
    sqwht2 = v['sqwh'] + v['sinw'] * t2
    upd('pz2px1 = h cos(θ1 - ψ1)', pz2px1 - h * mp.cos(D1_))
    upd('t4 = h cos(Δ1) cos(a3)/cos B', t4 - h * mp.cos(D1_) * mp.cos(a3) / mp.cos(B))
    upd('sqwht2 = (cosA - cos(ω+A))/cosA', sqwht2 - (mp.cos(A) - mp.cos(om_ + A)) / mp.cos(A))
    Nn = (v['dpx2'] * (sinp1 * (t2 * (v['pz3'] * cosp2 + v['px3'] * sinp2) - v['qx1'] * cosp2)
                       - sinp2 * (t3 * pz2px1 - v['px3'] * cosp1) + cosp1 * cosp2 * v['dpz32'] + t4 * sqwht2) + v['dpx3a'] * t4)
    arg = Nn / h**2
    upd('da の引数 = sin(Δ1 + Δ3)', arg - mp.sin(D1_ + D3_))
    # 紙の上の途中: ω の項が消える、残りが 0
    w_terms = v['dpx2'] * t4 * sqwht2 + v['dpx3a'] * t4
    upd('ω の項: dpx2 t4 sqwht2 + (dpx3a の ω の部分) t4 = 2 h sin a1 cosA t4', v['dpx2'] * t4 * sqwht2 + 2 * h * mp.sin(a1) * mp.cos(om_ + A) * t4 - 2 * h * mp.sin(a1) * mp.cos(A) * t4)
    upd('dpx3a の残り 2h cosB sin a3 × t4 = h^2 cosΔ1 sinΔ3', 2 * h * mp.cos(B) * mp.sin(a3) * t4 - h**2 * mp.cos(D1_) * mp.sin(D3_))
    Bk0 = (sinp1 * (t2 * h * mp.cos(D3_) - v['qx1'] * cosp2) - sinp2 * (t3 * h * mp.cos(D1_) - v['px3'] * cosp1) + cosp1 * cosp2 * v['dpz32'])
    upd('Bk0/h + cos a3 cosΔ1/cosB - cos a1 cosΔ3/cosA = 0', (Bk0 / h + mp.cos(a3) * mp.cos(D1_) / mp.cos(B) - mp.cos(a1) * mp.cos(D3_) / mp.cos(A)))
    upd('Θ - φ0 = Δ1 + Δ3', (th1 + th3 + om_ - v['phi0']) - (D1_ + D3_))
for k, w in worst.items():
    print('  %-62s 最大 %s' % (k, mp.nstr(w, 3)))

print('\n== (d) §8: K1 の b_y と xi の厳密解との差 ==')
x, h_, K, hb = sp.symbols('x h K h_b', positive=True)
by = hb + sp.integrate(K / (1 + h_ * x), x)
print('  ∫ K dx/(1 + h x) =', by - hb, '  (b_y(0) = h_b を境界条件として h_b + (K/h) ln(1 + h x))')
print('  (1 + h x) ∂b_y/∂x - K =', sp.simplify((1 + h_ * x) * sp.diff(hb + K / h_ * sp.log(1 + h_ * x), x) - K))
u = sp.Symbol('u')
print('  p_x\' の K の分 -(1 + u)(K/h) ln(1 + u) の u 級数:', sp.series(-(1 + u) * sp.log(1 + u), u, 0, 6))
ell, a, b_, xx = sp.symbols('ell a b x', real=True)
exact = xx * sp.exp(b_ * ell) + a * (sp.exp(b_ * ell) - 1) / b_
k64 = (xx + a * ell) / (1 - b_ * ell)
diff2 = sp.simplify(sp.series(exact - k64, ell, 0, 3).removeO())
print('  厳密 - K64 の ℓ の 2 次まで =', diff2, ' 期待 -(b ℓ^2/2)(a + b x):', sp.simplify(diff2 + b_ * ell**2 / 2 * (a + b_ * xx)))
q, dpzi, pzi, hh = sp.symbols('q dpzi pzi h', real=True)
sub = {a: -q * dpzi / pzi, b_: hh * q / pzi}
print('  a, b を代入: (厳密 - K64) - ℓ^2 h q^2/(2 pzi^2) (dpzi - h x) =', sp.simplify(diff2.subs(sub) - ell**2 * hh * q**2 / (2 * pzi**2) * (dpzi - hh * xx)))
# K64 の式の通分
print('  (x + aℓ)/(1 - bℓ) - (-ℓ dpzi q + pzi x)/(pzi - h ℓ q) =', sp.simplify(k64.subs(sub) - (-ell * dpzi * q + pzi * xx) / (pzi - hh * ell * q)))

print('\n== (e) §6.2 (5.33): p_z → p_v の近似の大きさ ==')
pv, pyy = mp.mpf('1.0'), mp.mpf('0.002')
for pyv in (mp.mpf('0.001'), mp.mpf('0.002'), mp.mpf('0.01')):
    pz = mp.sqrt(pv**2 - pyv**2)
    print('  p_v = 1、p_y = %s: (1/p_z)/(1/p_v) - 1 = %s、p_y^2/(2 p_v^2) = %s' % (mp.nstr(pyv, 3), mp.nstr(pv / pz - 1, 6), mp.nstr(pyv**2 / (2 * pv**2), 6)))
