# -*- coding: utf-8 -*-
r"""付録 B-1 (CALC の tbedgemaxwell の COD の更新が因子 p を落とす疑い) の大きさの見積りと、差し替えの試験 (05_bend.md §9、付録 B-1)。

(a) 1 つの面での、(5.36) の p^2 型と K64 の p 型 (bf*pr) の Delta x、Delta z の差の見積り:
      Delta x の差 = (b y^2/rho_b) (p^2 - p)/p_v^3、Delta z の差 = -(b y^2/rho_b) p_x (p - 1)/p_v^3。
    入口と出口で符号が逆 (f_out = -f_in) なので、要素全体では差の差しか残らない。
(b) lattice._tbedgemaxwell_cod を p^2 型 (Delta x に pr^2、Delta z に pxi*pr) に差し替えたとき、lattice._bend_cod_update (CALC の COD の更新、
    Python 版) と追跡の写像の差がどうなるか。差し替えは試験のためだけで、pysad のファイルは変えない。
(c) F1 ≠ 0 のときの CALC の 4x4 と追跡のヤコビアンの差 (4.4e-5 の成分と、y・δ・F1・E への依存)。これは (b) の B-1 とは別の差。
"""
import math
import numpy as np
import _common as C

C.start_modes(__file__)
from pysad import lattice as lat

orig = lat._tbedgemaxwell_cod


def p2_version(cod, rhob):
    """(5.36) の p^2 型 (K64 の tbedgemaxwell の解析行列と整合する形)。Delta x に pr^2、Delta z に pxi*pr。"""
    cod = cod.copy()
    pr = 1.0 + cod[5]
    pxi = np.clip(cod[1], -pr, pr)
    pvi2 = max(1e-30, (pr - pxi) * (pr + pxi))
    pvi = math.sqrt(pvi2)
    yi = cod[2]
    yr = (yi / rhob)**2
    a = 1.0 - yr / 6.0
    af = a * yi / rhob / pvi
    b = 0.5 - yr / 24.0
    bf = b * yi**2 / rhob / pvi2 / pvi
    cod[0] = cod[0] + bf * pr * pr
    cod[3] = np.clip(cod[3] - af * pxi, -pr, pr)
    cod[4] = cod[4] - bf * pxi * pr
    return cod


print('== (a) 1 つの面での p^2 型と p 型の差の見積り (rho_b = 10 m、y = 10 mm) ==')
rhob = 10.0
y = 0.01
for dl in (0.05, 0.02, -0.05):
    pr = 1 + dl
    for qx in (0.0, 0.05, 0.1):
        px = qx * pr
        pv = math.sqrt(pr**2 - px**2)
        b = 0.5 - (y / rhob)**2 / 24
        ddx = b * y**2 / rhob * (pr**2 - pr) / pv**3
        ddz = -b * y**2 / rhob * px * (pr - 1) / pv**3
        print('  δ = %+.2f、q_x = %.2f:  Δx の差 %.2e、Δz の差 %.2e   (δ y^2/(2 rho_b) = %.2e)' % (dl, qx, ddx, ddz, dl * y**2 / (2 * rhob)))
print('  (δ y^2/(2 rho_b) = 2.5e-7 は 1 つの面の Δx の差の大きさ。入口と出口の差が打ち消し合い、要素全体の差はずっと小さい: (b)。)')

print('\n== (b) CALC の COD の更新 (Python 版 _bend_cod_update) と追跡の写像の差。tbedgemaxwell を p^2 型に差し替える試験 ==')
print('  初期値 (x, px, y, py, z, δ) = (1 mm, 0.5 mrad, y, 0.2 mrad, 0, δ)。差 = CALC - 追跡 (x px y py z)')
for extra, e1, e2, name in (('', 0.0, 0.0, '扇形'), ('E1=0.5 E2=0.5', 0.05, 0.05, '矩形 (E1 = E2 = 0.5)'), ('E1=0.3 E2=0.7', 0.03, 0.07, '非対称 (E1 = 0.3、E2 = 0.7)')):
    ffs, el = C.make_bend(1.0, 0.1, extra)
    for yv, dl in ((0.01, 0.05), (0.005, 0.05), (0.01, 0.02), (0.01, -0.05)):
        c = np.array([0.001, 0.0005, yv, 0.0002, 0.0, dl])
        tr = C.track_one(ffs, c)
        res = []
        for f in (orig, p2_version):
            lat._tbedgemaxwell_cod = f
            res.append(lat._bend_cod_update(c, 1.0, 0.1, e1, e2)[:5] - tr[:5])
        lat._tbedgemaxwell_cod = orig
        print('  %-26s y = %4.1f mm、δ = %+.2f:' % (name, yv * 1e3, dl))
        print('      K64 (p 型)  :', np.array2string(res[0], precision=2))
        print('      p^2 型に差替:', np.array2string(res[1], precision=2))
        if name == '扇形' and yv == 0.01 and dl == 0.05:
            # 入口と出口の差の差し引き: y_out は追跡の y
            pr = 1 + dl
            b = 0.5 - (yv / rhob)**2 / 24
            eps_in = b * yv**2 / rhob * (pr**2 - pr) / (pr**3)
            yo = tr[2]
            eps_out = b * yo**2 / rhob * (pr**2 - pr) / (pr**3)
            print('      見積り: 入口の Δx の差 %.2e、出口 (y_out = %.4f m) %.2e、符号が逆なので、CALC - 追跡 = -入口 + 出口 = %+.2e (計算値の x の差 +9.7e-9 と同じ桁。残りは円弧を通る間の伝わり方の違い)' % (eps_in, yo, eps_out, eps_out - eps_in))
            print('      見積り: 入口の Δx の差が本体で dz/ds = -h x を通して効く分 eps_in L/rho_0 = %.2e (z の差 +2.0e-8 と同じ桁)' % (eps_in * 1.0 / 10.0))
print('  p^2 型に差し替えると、x と px の差は 1e-16 台に下がり、z の差は 2e-8 から 1.4e-9 に下がる。残る z の差 (1.4e-9) は y によらず δ に比例 (δ = 0.05、0.02、-0.05 で')
print('  -1.35e-9、-5.6e-10、+1.6e-9)。B-1 とは別の効果で、原因は特定していない。')

print('\n== (c) F1 ≠ 0 のときの CALC の 4x4 と追跡のヤコビアンの差 ==')


def jac(ffs, c, h=1e-6):
    c = np.array(c, float)
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (C.track_one(ffs, cp) - C.track_one(ffs, cm)) / (2 * h)
    return M


print('  %-34s %-12s %-8s %-12s %-20s' % ('場合', '軌道 (y, δ)', '4x4 最大', '最大の成分', '最大の差'))
labels = 'x px y py'.split()
for extra in ('', 'F1=0.1 FRINGE=1', 'E1=0.5 E2=0.5 F1=0.1 FRINGE=1', 'E1=0.5 E2=0.5 F1=0.05 FRINGE=1', 'E1=0.5 E2=0.5 F1=0.2 FRINGE=1'):
    ffs, el = C.make_bend(1.0, 0.1, extra)
    for yv, dl in ((0.005, 0.05), (0.01, 0.05), (0.01, 0.02), (0.01, -0.05), (0.01, 0.0), (0.0, 0.05)):
        c = np.array([0.0, 0.0, yv, 0.0, 0.0, dl])
        Jt = jac(ffs, c)
        Mc, _ = lat._el_cod_step(el, c, c[5])
        d = (Mc - Jt)[:4, :4]
        i, j = np.unravel_index(np.argmax(np.abs(d)), d.shape)
        print('  %-34s y = %4.1f mm、δ = %+.2f  %.1e  ∂%s_out/∂%s   %+.2e' % (extra or 'F1 = 0、扇形', yv * 1e3, dl, np.abs(d).max(), labels[i], labels[j], d[i, j]))
print('  F1 = 0 では、最大でも 3e-8 (y = 10 mm、|δ| = 0.05)。F1 ≠ 0 では ∂x_out/∂p_y の成分が 2e-5〜6e-5 になり、')
print('  y に比例し、ほぼ δ に比例し、F1 (0.05〜0.2) にも E1・E2 (0、0.5) にもよらない。したがって、この差は (b) の B-1 (F1 = 0 で 3e-10〜2e-8) とは別の、')
print('  F1 ≠ 0 の経路の違いによる。F1 つきの 4x4 の差 (4.4e-5) と、段ごとの解析行列の不整合 (c2) の 4.5e-5 は、数値が近いだけで同じ原因とは言えない。')
