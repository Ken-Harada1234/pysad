# -*- coding: utf-8 -*-
"""BEND の転送行列 (CALC・Emittance[]) と追跡の写像の関係 (05_bend.md §9)。

K1 = 0 の BEND 1 台 (L = 1 m、ANGLE = 0.1 rad、電子 3 GeV)。
(a) 追跡の写像の中心差分ヤコビアン J_track (6x6、刻み 1e-6) と、CALC の行列 (lattice._el_cod_step) の 4x4 と δ の列の差。
    場合: 扇形、矩形 (E1 = E2 = 0.5)、非対称、F1 = 0.1 (FRINGE = 1)。設計軌道と、ずれた軌道。
(b) CALC の COD の更新と追跡の写像の差 (δ と y が大きいとき)。
(c) CALC の解析行列 (lattice._bend_cod_matrix) と、CALC 自身の COD の更新の数値微分の差。
(d) 追跡のヤコビアンのシンプレクティック性 |J^T S J - S|、z の行の条件。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
from pysad import lattice as lat

S6 = np.zeros((6, 6))
for i in range(3):
    S6[2 * i, 2 * i + 1] = 1
    S6[2 * i + 1, 2 * i] = -1


def jac_track(ffs, c, h=1e-6):
    c = np.array(c, float)
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (C.track_one(ffs, cp) - C.track_one(ffs, cm)) / (2 * h)
    return M


cases = [('扇形', ''), ('矩形 (E1 = E2 = 0.5)', 'E1=0.5 E2=0.5'), ('非対称 (E1 = 0.3、E2 = 0.7)', 'E1=0.3 E2=0.7'),
         ('F1 = 0.1、FRINGE = 1、E1 = E2 = 0.5', 'E1=0.5 E2=0.5 F1=0.1 FRINGE=1')]
orbits = [('設計軌道 (δ = 0)', [0, 0, 0, 0, 0, 0.0]),
          ('小さいずれ (δ = 0.002)', [0.001, 0.0005, 0.0003, 0.0002, 0, 0.002]),
          ('大きいずれ (y = 10 mm、δ = 0.05)', [0.001, 0.0005, 0.01, 0.0002, 0, 0.05])]
print('== (a)(b) CALC の 4x4・δ の列と COD の更新を、追跡と比べる (k64 と v2 で同じ値) ==')
print('  %-36s %-34s %-10s %-10s %-10s' % ('場合', '軌道', '4x4', 'δ の列', 'COD (x px y py z)'))
for name, extra in cases:
    ffs, el = C.make_bend(1.0, 0.1, extra)
    for oname, c in orbits:
        c = np.array(c, float)
        Jt = jac_track(ffs, c)
        Mc, cod = lat._el_cod_step(el, c, c[5])
        tr = C.track_one(ffs, c)
        print('  %-36s %-34s %-10.1e %-10.1e %-10.1e' % (name, oname, np.abs(Mc[:4, :4] - Jt[:4, :4]).max(), np.abs(Mc[:4, 5] - Jt[:4, 5]).max(), np.abs(cod[:5] - tr[:5]).max()))

print('\n== (c) CALC の解析行列 (tbedge・tbendebody0 の trans) と、CALC 自身の COD の更新の数値微分 (F1 = 0) ==')
print('  %-36s %-34s %-12s' % ('場合', '軌道', '4x4・δ の列の差'))
for name, extra in cases[:3]:
    e1 = 0.5 * 0.1 if 'E1=0.5' in extra else (0.3 * 0.1 if 'E1=0.3' in extra else 0.0)
    e2 = 0.5 * 0.1 if 'E2=0.5' in extra else (0.7 * 0.1 if 'E2=0.7' in extra else 0.0)
    for oname, c in orbits:
        c = np.array(c, float)
        cod_out, M = lat._bend_cod_matrix(c, 1.0, 0.1, e1, e2)
        Mn = np.zeros((6, 6))
        for j in range(6):
            cp, cm = c.copy(), c.copy()
            cp[j] += 1e-6
            cm[j] -= 1e-6
            Mn[:, j] = (lat._bend_cod_update(cp, 1.0, 0.1, e1, e2) - lat._bend_cod_update(cm, 1.0, 0.1, e1, e2)) / 2e-6
        sel = [0, 1, 2, 3, 5]
        print('  %-36s %-34s %-12.1e' % (name, oname, np.abs(M[np.ix_(sel, sel)] - Mn[np.ix_(sel, sel)]).max()))

print('\n== (d) 追跡のヤコビアンのシンプレクティック性と z の行 ==')
c = np.array([0.001, 0.0005, 0.01, 0.0002, 0.0, 0.05])
for name, extra in cases:
    ffs, el = C.make_bend(1.0, 0.1, extra)
    Jt = jac_track(ffs, c)
    Mc, _ = lat._el_cod_step(el, c, c[5])
    zf = np.array([-Jt[0, j] * Jt[1, 5] + Jt[1, j] * Jt[0, 5] - Jt[2, j] * Jt[3, 5] + Jt[3, j] * Jt[2, 5] for j in range(4)])
    print('  %-36s |J_track^T S J_track - S| = %.1e、z の行 T(5,1..4) と条件で組んだ値の差 %.1e、CALC の 6x6 の |M^T S M - S| = %.1e' % (
        name, np.abs(Jt.T @ S6 @ Jt - S6).max(), np.abs(zf - Jt[4, :4]).max(), np.abs(Mc.T @ S6 @ Mc - S6).max()))

print('\n== (c2) (c) の差の出所: 段ごとに、解析行列 (trans) と同じ段の COD 更新の数値微分を比べる (y = 10 mm、δ = 0.05、ρb = 10 m、ψ = 0.05 rad) ==')
rhob_ = 10.0
c = np.array([0.001, 0.0005, 0.01, 0.0002, 0.0, 0.05])


def numj(f, c, h=1e-6):
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (f(cp) - f(cm)) / (2 * h)
    return M


sel = [0, 1, 2, 3, 5]
stages = [('tbedgemaxwell (Maxwell フリンジ)', lat._tbedgemaxwell_mat(c, rhob_), lambda cc: lat._tbedgemaxwell_cod(cc, rhob_)),
          ('tbendebody0 (円弧の本体)', lat._tbendebody0_mat(c, 0.1, rhob_, 0.0), lambda cc: lat._tbendebody0_cod(cc, 0.1, rhob_, 0.0)),
          ('tbedgedrift (面までの直進)', lat._tbedgedrift_mat(c, 0.05), lambda cc: lat._tbedgedrift_cod(cc, 0.05)),
          ('tbedgebody (面から扇形の面まで)', lat._tbedgebody_mat(c, rhob_, 0.05), lambda cc: lat._tbedgebody_cod(cc, rhob_, 0.05))]
for name, Ma, f in stages:
    d = np.abs(Ma - numj(f, c))[np.ix_(sel, sel)]
    i, j = np.unravel_index(np.argmax(d), d.shape)
    print('  %-34s 最大の差 %.1e (行 %s、列 %s)' % (name, d.max(), 'x px y py δ'.split()[i], 'x px y py δ'.split()[j]))
Mm = lat._tbedgemaxwell_mat(c, rhob_)
pr = 1 + c[5]
pv = math.sqrt(pr**2 - c[1]**2)
yr = (c[2] / rhob_)**2
a_ = 1 - yr / 6
print('  tbedgemaxwell の ∂x/∂y: 解析行列 trans(1,3) = %.6e、COD 更新 Δx = bf pr の数値微分 = %.6e、比 = %.5f (= 1 + δ = %.2f)' % (
    Mm[0, 2], numj(lambda cc: lat._tbedgemaxwell_cod(cc, rhob_), c)[0, 2], Mm[0, 2] / numj(lambda cc: lat._tbedgemaxwell_cod(cc, rhob_), c)[0, 2], pr))

print('\n== (e) 設計軌道の 2x2 行列: 端面の収束 (x: +tan ψ/ρ、y: -tan ψ/ρ) と扇形の本体の積 ==')
rho_ = 10.0
th = 0.1


def Ex(psi):
    return np.array([[1, 0], [math.tan(psi) / rho_, 1]])


def Ey(psi):
    return np.array([[1, 0], [-math.tan(psi) / rho_, 1]])


Sx_ = np.array([[math.cos(th), rho_ * math.sin(th)], [-math.sin(th) / rho_, math.cos(th)]])
Sy_ = np.array([[1, 1.0], [0, 1]])
for name, extra, e1, e2 in (('扇形', 'E1=0 E2=0', 0.0, 0.0), ('矩形 (E1 = E2 = 0.5)', 'E1=0.5 E2=0.5', 0.5 * th, 0.5 * th), ('非対称 (E1 = 0.3、E2 = 0.7)', 'E1=0.3 E2=0.7', 0.3 * th, 0.7 * th)):
    ffs, el = C.make_bend(1.0, th, extra)
    Jt = jac_track(ffs, np.zeros(6))
    Mx = Ex(e2) @ Sx_ @ Ex(e1)
    My = Ey(e2) @ Sy_ @ Ey(e1)
    dx_ = np.abs(Jt[np.ix_([0, 1], [0, 1])] - Mx).max()
    dy_ = np.abs(Jt[np.ix_([2, 3], [2, 3])] - My).max()
    print('  %-28s 追跡のヤコビアンと 2x2 の積の差: x %.1e、y %.1e   (x の R21 = %+.6f、y の R43 = %+.6f)' % (name, dx_, dy_, Jt[1, 0], Jt[3, 2]))
print('  (ψ の符号は、E > 0 が面を磁石の外側へ開く向き。x の端の収束は +tan ψ/ρ、y は -tan ψ/ρ)')

print('\n== (f) F1 ≠ 0 (E = 0、FRINGE = 1、F1 = 0.1) の CALC と追跡のヤコビアンの差 (小さいずれ: x = 1 mm、px = 0.5 mrad、y = 0.3 mm、py = 0.2 mrad、δ = 0.002) ==')
np.set_printoptions(precision=2, linewidth=200)
ffs, el = C.make_bend(1.0, 0.1, 'FRINGE=1 F1=0.1')
c = np.array([0.001, 0.0005, 0.0003, 0.0002, 0, 0.002])
Jt = jac_track(ffs, c)
Mc, _ = lat._el_cod_step(el, c, c[5])
fb = lat._bend_fringe_params(el)


def calc_map(cc):
    return lat._bend_cod_update(cc, 1.0, 0.1, 0.0, 0.0, 0.1, fb[0], fb[1], fb[2], fb[3])


Jc = numj(calc_map, c)
print('  CALC の解析行列 - 追跡のヤコビアン (行 x px y py、列 x px y py z δ):')
print(str((Mc - Jt)[:4, :]).replace('\n', '\n  '))
print('  CALC の COD 更新の数値微分 - 追跡のヤコビアン (同じ形): 最大 %.1e (x の行の y・py の列が主)、解析行列 - 自身の数値微分: 最大 %.1e' % (np.abs((Jc - Jt)[:4, :]).max(), np.abs((Mc - Jc)[:4, :]).max()))
print('  F1 = 0 では同じ軌道で 4x4 の差は 1.4e-12 (上の (a))。F1 ≠ 0 では、CALC の写像そのもの (COD の更新の微分) が追跡の写像と違う。原因は特定していない。')
