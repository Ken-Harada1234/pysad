# -*- coding: utf-8 -*-
"""11 章 (11_misalignment_and_coord.md) の途中式の確認。本文の各中間式を sympy で 1 つずつ確かめる。

(A) R3 R2 R1 と R1 R2(-) R3 の 3×3 の積: 中間の R2 R1、R1 R2(-) と、最後の積の各行 (COORD、(11.6)(11.7))。
(B) D = diag(1,-1,-1) による共役: D R1 D、D R2 D、D R3 D の成分。
(C) ロドリゲスの式による R_t(Θ) R_n(-χ2) の積の 9 成分 (中間の 3 項の和と、Fortran の rr との差)。
(D) tbgrote の 6×6: 写像を記号のまま微分した 36 成分と、Fortran の式 (tchge.f:249-289) の差。(5,6) 成分の導出。
(E) スナップショットの座標 (11.3) が dr∧dP を dx∧dpx + dy∧dpy + dz∧dδ に写すこと (シンプレクティック性の本文の議論)。
(F) 回転で剛体運動の x 成分 dx1 = dx c0 - dz s0 など。
"""
import sympy as sp
import _common as C

C.start(__file__)

c1, s1, c2, s2, c3, s3 = sp.symbols('c1 s1 c2 s2 c3 s3', real=True)
R3 = sp.Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
R2 = sp.Matrix([[1, 0, 0], [0, c2, s2], [0, -s2, c2]])
R1 = sp.Matrix([[c1, 0, s1], [0, 1, 0], [-s1, 0, c1]])
R2m = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])


def show(name, M):
    print('  %s = %s' % (name, M.tolist()))


print('== (A) 3×3 の積の中間式 ==')
R2R1 = R2 * R1
show('R2 R1', R2R1)
show('R2 R1 の期待 [[c1,0,s1],[-s1 s2,c2,s2 c1],[-s1 c2,-s2,c1 c2]] との差', R2R1 - sp.Matrix([[c1, 0, s1], [-s1 * s2, c2, s2 * c1], [-s1 * c2, -s2, c1 * c2]]))
for i in range(3):
    row = R3[i, :] * R2R1
    print('  R3 の第 %d 行 × (R2 R1) = %s' % (i + 1, [sp.expand(e) for e in row]))
M0 = sp.expand(R3 * R2R1)
exp0 = sp.Matrix([[c1 * c3 + s1 * s2 * s3, -c2 * s3, s1 * c3 - c1 * s2 * s3],
                  [-s1 * s2 * c3 + c1 * s3, c2 * c3, c1 * s2 * c3 + s1 * s3],
                  [-s1 * c2, -s2, c1 * c2]])
show('R3 R2 R1 - (11.6 の本文の行列)', sp.simplify(M0 - exp0))
R1R2m = R1 * R2m
show('R1 R2(-)', R1R2m)
show('R1 R2(-) の期待 [[c1,s1 s2,s1 c2],[0,c2,-s2],[-s1,c1 s2,c1 c2]] との差', R1R2m - sp.Matrix([[c1, s1 * s2, s1 * c2], [0, c2, -s2], [-s1, c1 * s2, c1 * c2]]))
for j in range(3):
    col = R1R2m * R3[:, j]
    print('  (R1 R2(-)) × R3 の第 %d 列 = %s' % (j + 1, [sp.expand(e) for e in col]))
M7 = sp.expand(R1R2m * R3)
exp7 = sp.Matrix([[c1 * c3 + s1 * s2 * s3, s1 * s2 * c3 - c1 * s3, s1 * c2],
                  [c2 * s3, c2 * c3, -s2],
                  [-s1 * c3 + c1 * s2 * s3, c1 * s2 * c3 + s1 * s3, c1 * c2]])
show('R1 R2(-) R3 - (11.7)', sp.simplify(M7 - exp7))
show('R3 R2(-) R1 - (11.7)  (非零 = 一致しない)', sp.simplify(R3 * R2m * R1 - exp7))
print('  R3 R2(-) R1 の (1,1) 成分 = %s' % sp.expand((R3 * R2m * R1)[0, 0]))
# 付録 A の MULT の式: R3 R2(-) R1(-)
R1m = sp.Matrix([[c1, 0, -s1], [0, 1, 0], [s1, 0, c1]])
R21 = R2m * R1m
show('R2(-) R1(-)', R21)
show('R3 R2(-) R1(-) = R3 R2 R1 の s1,s2 を -s1,-s2 に替えたもの',
     sp.simplify(R3 * R21 - M0.subs({s1: -s1, s2: -s2}, simultaneous=True)))
show('R3 R2(-) R1(-)', sp.expand(R3 * R21))

print('\n== (B) D = diag(1,-1,-1) による共役 ==')
D = sp.diag(1, -1, -1)
show('D R1 D', D * R1 * D)
show('D R1 D - R1(s1 -> -s1)', D * R1 * D - R1.subs(s1, -s1))
show('D R2 D', D * R2 * D)
show('D R2 D - R2', D * R2 * D - R2)
show('D R3 D', D * R3 * D)
show('D R3 D - R3(s3 -> -s3)', D * R3 * D - R3.subs(s3, -s3))
Mt = sp.expand(M0.T)
show('R1(c1,s1) R2(-) R3(c3,s3) - D M0^T D', sp.simplify(M7 - D * Mt * D))
show('R1^T R2^T R3^T - M0^T', sp.simplify(R1.T * R2.T * R3.T - Mt))
show('R2^T - R2(-)', R2.T - R2m)

print('\n== (C) ロドリゲスの式 R_t(Θ) R_n(-χ2): 9 成分 ==')
c0, s0, cT, sT, cc2, ss2 = sp.symbols('c0 s0 cT sT cc2 ss2', real=True)
t = sp.Matrix([-s0, 0, c0])
n = sp.Matrix([c0, 0, s0])


def cross(a):
    return sp.Matrix([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])


ru = {c0**2: 1 - s0**2}
Kt, Kn = cross(t), cross(n)
show('[t]x', Kt)
show('[n]x', Kn)
show('[t]x^2', sp.expand(Kt * Kt))
show('[n]x^2', sp.expand(Kn * Kn))
At = sp.eye(3) + sT * Kt + (1 - cT) * Kt * Kt
Bn = sp.eye(3) - ss2 * Kn + (1 - cc2) * Kn * Kn
At = At.subs(s0**2, 1 - c0**2)
show('R_t(Θ)', sp.simplify(At))
show('R_n(-χ2)', sp.simplify(Bn))
Rfull = sp.expand(At * Bn)
# Fortran (tchge.f:161-172) rr の入口
dct, dc2 = cT - 1, cc2 - 1
rr = sp.Matrix([
    [1 + dct * c0**2 + dc2 * s0**2 + c0 * s0 * ss2 * sT, s0 * ss2 - c0 * cc2 * sT, c0 * ((dct - dc2) * s0 - c0 * ss2 * sT)],
    [-cT * s0 * ss2 + c0 * sT, 1 + cc2 * dct + dc2, c0 * cT * ss2 + s0 * sT],
    [s0 * (c0 * (dct - dc2) + s0 * ss2 * sT), -c0 * ss2 - cc2 * s0 * sT, 1 + dc2 * c0**2 + dct * s0**2 - c0 * s0 * ss2 * sT]])
for i in range(3):
    for j in range(3):
        terms = [sp.expand(At[i, k] * Bn[k, j]) for k in range(3)]
        dif = sp.simplify((Rfull[i, j] - rr[i, j]).subs(s0, sp.sqrt(1 - c0**2)))
        print('  R_%d%d = Σ_k A_%dk B_k%d = %s   [Fortran との差 %s]' % (i + 1, j + 1, i + 1, j + 1, ' + '.join('(%s)' % sp.factor(x) for x in terms), dif))

print('\n== (D) tbgrote の 6×6 ==')
x, px, y, py, z, dl_ = sp.symbols('x px y py z delta', real=True)
dr = sp.Matrix(3, 3, lambda i, j: sp.Symbol('dr%d%d' % (i + 1, j + 1), real=True))
xl0, dzz = sp.symbols('xl0 dzz', real=True)
pr = 1 + dl_
pzi = sp.sqrt(pr**2 - px**2 - py**2)
Rm = dr + sp.eye(3)
r_vec = sp.Matrix([x + xl0, y, dzz])        # (xi, yi, zi): 支点からの位置
P_vec = sp.Matrix([px, py, pzi])
# 本文: r' - r = (R - 1)(r - r_p) = dr (r - r_p)、 (xf, yf, zf) = r + dr (r - r_p) の面内成分
xf = dr[0, :].dot(r_vec) + x
yf = dr[1, :].dot(r_vec) + y
zf = dr[2, :].dot(r_vec)
Pf = Rm * P_vec
pxf, pyf, pzf = Pf
xn = xf - pxf * zf / pzf
yn = yf - pyf * zf / pzf
zn = z + pr * zf / pzf
fmap = sp.Matrix([xn, pxf, yn, pyf, zn, dl_])
vars6 = [x, px, y, py, z, dl_]
Jsym = fmap.jacobian(vars6)
# Fortran の式
r11, r22, r33 = 1 + dr[0, 0], 1 + dr[1, 1], 1 + dr[2, 2]
dpzidpx, dpzidpy, dpzidp = -px / pzi, -py / pzi, pr / pzi
pxfpzf, pyfpzf, zfpzf, prpzf = pxf / pzf, pyf / pzf, zf / pzf, pr / pzf
dpzf = (Pf[2] - pr)                                    # (R P)_z - pr
dpzi = pzi - pr
T = sp.zeros(6, 6)
T[1, 1] = r11 + dr[0, 2] * dpzidpx
T[1, 3] = dr[0, 1] + dr[0, 2] * dpzidpy
T[1, 5] = dr[0, 2] * dpzidp
T[3, 1] = dr[1, 0] + dr[1, 2] * dpzidpx
T[3, 3] = r22 + dr[1, 2] * dpzidpy
T[3, 5] = dr[1, 2] * dpzidp
T[0, 0] = r11 - dr[2, 0] * pxfpzf
T[0, 1] = -zfpzf * (T[1, 1] - pxfpzf * (dr[2, 0] + r33 * dpzidpx))
T[0, 2] = dr[0, 1] - dr[2, 1] * pxfpzf
T[0, 3] = -zfpzf * (T[1, 3] - pxfpzf * (dr[2, 1] + r33 * dpzidpy))
T[0, 5] = -zfpzf * (dr[0, 2] - r33 * pxfpzf) * dpzidp
T[2, 0] = dr[1, 0] - dr[2, 0] * pyfpzf
T[2, 1] = -zfpzf * (T[3, 1] - pyfpzf * (dr[2, 0] + r33 * dpzidpx))
T[2, 2] = r22 - dr[2, 1] * pyfpzf
T[2, 3] = -zfpzf * (T[3, 3] - pyfpzf * (dr[2, 1] + r33 * dpzidpy))
T[2, 5] = -zfpzf * (dr[1, 2] - r33 * pyfpzf) * dpzidp
T[4, 0] = dr[2, 0] * prpzf
T[4, 1] = -zfpzf * prpzf * (dr[2, 0] + r33 * dpzidpx)
T[4, 2] = dr[2, 1] * prpzf
T[4, 3] = -zfpzf * prpzf * (dr[2, 1] + r33 * dpzidpy)
T[4, 4] = 1
T[5, 5] = 1
T56_live = r33 * (dpzf * prpzf + dpzi / pzi) - dr[2, 2]
T56_comm = zfpzf * (1 - pr / pzf * r33 * dpzidp)
import random
random.seed(3)
subs = {s: sp.Rational(random.randint(-9, 9), 997) for s in list(dr) + [x, px, y, py, z, xl0, dzz]}
subs[dl_] = sp.Rational(2, 1000)
worst = 0
names = ['x', 'px', 'y', 'py', 'z', 'δ']
for i in range(6):
    for j in range(6):
        if (i, j) == (4, 5):
            continue
        d = sp.N((Jsym[i, j] - T[i, j]).subs(subs), 20)
        worst = max(worst, abs(d))
print('  記号のまま微分した 6×6 と Fortran の式の差 (乱数の dr、5 成分除く 35 成分) の最大 = %s' % sp.N(worst, 5))
print('  (5,6): 記号微分 - コメントアウトされた式 zf/pzf (1 - pr/pzf r33 dpzidp) = %s' % sp.N((Jsym[4, 5] - T56_comm).subs(subs), 5))
print('  (5,6): 記号微分 - 現在の式 r33 (dpzf prpzf + dpzi/pzi) - dr33 = %s (非零 = 一致しない)' % sp.N((Jsym[4, 5] - T56_live).subs(subs), 5))
subs0 = dict(subs)
for s in dr:
    subs0[s] = 0
print('  dr = 0: 記号微分 (5,6) = %s、現在の式 = %s (= dpzi (pr + 1)/pzi。dpzi = %s)' % (
    sp.N(Jsym[4, 5].subs(subs0), 6), sp.N(T56_live.subs(subs0), 6), sp.N(dpzi.subs(subs0), 6)))
print('  dr = 0 の現在の式の閉じた形: dpzi (pr/pzi + 1/pzi) - 差 = %s' % sp.N((T56_live - dpzi * (pr + 1) / pzi).subs(subs0), 5))
# 一般の関係式 (本文): 位置型 v∈{x,y}、運動量型 w∈{px,py,δ}
print('  公式の確認: ∂x'' /∂x = ∂xf/∂x - (pxf/pzf) ∂zf/∂x、∂z'' /∂x = (pr/pzf) ∂zf/∂x の差 = %s、%s' % (
    sp.N((Jsym[0, 0] - (sp.diff(xf, x) - pxfpzf * sp.diff(zf, x))).subs(subs), 5),
    sp.N((Jsym[4, 0] - prpzf * sp.diff(zf, x)).subs(subs), 5)))

print('\n== (E) スナップショット座標 (11.3) の pullback ==')
pz_ = sp.sqrt(pr**2 - px**2 - py**2)
u = sp.Matrix([px, py, pz_]) / pr
rr3 = sp.Matrix([x, y, 0]) + z * u
PP = pr * u
J0 = sp.Matrix.vstack(rr3, PP).jacobian([x, px, y, py, z, dl_])      # (r, P) の x,px,... による微分
Om3 = sp.zeros(6, 6)
for k in range(3):
    Om3[k, 3 + k] = 1
    Om3[3 + k, k] = -1
Pull = sp.simplify(J0.T * Om3 * J0)
print('  J0^T Ω3 J0 = %s' % Pull.tolist())
Om6 = sp.zeros(6, 6)
for a, b in ((0, 1), (2, 3), (4, 5)):
    Om6[a, b] = 1
    Om6[b, a] = -1
print('  dx∧dpx + dy∧dpy + dz∧dδ の行列との差 = %s' % sp.simplify(Pull - Om6).tolist())
print('  同じ pullback を数値で: 差の最大 = %s' % max(abs(sp.N(e.subs({x: 0.3, y: -0.2, z: 0.5, px: 0.11, py: -0.07, dl_: 0.02}), 6)) for e in (Pull - Om6)))

print('\n== (F) dx1・dz1: (dx, dz) を φg だけ回した成分 ==')
phi, dx_, dz_ = sp.symbols('phi dx dz', real=True)
Rphi = sp.Matrix([[sp.cos(phi), -sp.sin(phi)], [sp.sin(phi), sp.cos(phi)]])
v = Rphi * sp.Matrix([dx_, dz_])
print('  入口: dx1 = dx c0 - dz s0 = (R(φg)(dx,dz))_1: 差 = %s、dz1 = dz c0 + dx s0 = (R(φg)(dx,dz))_2: 差 = %s' % (
    sp.simplify(v[0] - (dx_ * sp.cos(phi) - dz_ * sp.sin(phi))), sp.simplify(v[1] - (dz_ * sp.cos(phi) + dx_ * sp.sin(phi)))))
