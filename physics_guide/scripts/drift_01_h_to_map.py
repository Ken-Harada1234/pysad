# -*- coding: utf-8 -*-
"""ドリフトのハミルトニアンから写像まで (02_drift.md §3〜§6、式 (2.1)〜(2.12)) を確かめる。

H = -sqrt(p^2 - px^2 - py^2) + E/(p0 v0)   (p = 1 + δ、E/(p0 v0) = h1 h0 / p0^2)
(a) sympy: ハミルトンの方程式 x' = px/pz、z' = -p/pz + v/v0 を偏微分から作り、dv の閉じた式 (01 章 (1.9)) と z' = -(px^2+py^2)/((p+pz) pz) - dv が一致すること。
(b) 閉じた式 (x += px L/pz など) と、正準方程式の数値積分 (DOP853) と、pysad の追跡 (DRIFT 1 台) の比較。
(c) 厳密なドリフト = 近軸のドリフト (H2 で K = 0) と運動学的ドリフト D (ΔH) の積。両者は運動量だけの関数なので交換する。
(d) DISKIN (近軸のドリフト): Fortran の式 (tturn.f:296 の近軸式) と、pysad が DISKIN を無視する (厳密な式を使う) こと。
"""
import math
import io
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

# 電子、3 GeV (台本の MOMENTUM = 3GEV)。p0 は mc 単位の設計運動量 β0γ0、h0 = γ0
MC2 = 0.51099895e6
P_EV = 3.0e9
p0 = P_EV / MC2
h0 = math.sqrt(1 + p0**2)
print('電子、運動量 3 GeV/c: p0 = β0γ0 = %.4f、h0 = γ0 = %.4f' % (p0, h0))


def dv_closed(dp):
    """01 章 (1.9): dv = -δ(2+δ)/(h1 (h1 + (1+δ) h0))"""
    h1 = math.sqrt(1 + (p0 * (1 + dp))**2)
    return -dp * (2 + dp) / (h1 * (h1 + (1 + dp) * h0))


print('\n== (a) sympy: H の偏微分と dv ==')
x, y, px, py, p = sp.symbols('x y p_x p_y p', real=True)
P0, H0 = sp.symbols('p0 h0', positive=True)
h1 = sp.sqrt(1 + (P0 * p)**2)
Epart = h1 * H0 / P0**2                                   # E/(p0 v0) (E = h1 mc^2、v0 = β0 c = p0/h0 c、p0 の運動量は p0 mc)
pz = sp.sqrt(p**2 - px**2 - py**2)
Hd = -pz + Epart
xp = sp.diff(Hd, px)
yp = sp.diff(Hd, py)
zp = sp.diff(Hd, p)
pxp = -sp.diff(Hd, x)
print('  x\' = ∂H/∂px = %s、y\' = ∂H/∂py = %s' % (sp.simplify(xp), sp.simplify(yp)))
print('  px\' = -∂H/∂x = %s、py\' = %s (H は x, y を含まない)' % (pxp, -sp.diff(Hd, y)))
print('  z\' = ∂H/∂p = %s' % sp.simplify(zp))
vv0 = sp.simplify(sp.diff(Epart, p))
print('  ∂[E/(p0 v0)]/∂p = %s  (= v/v0 = p h0/h1 と一致するか: 差 %s)' % (vv0, sp.simplify(vv0 - p * H0 / h1)))
# z' = -(px^2+py^2)/((p+pz) pz) - dv   (dv = 1 - v/v0)
zp_claim = -(px**2 + py**2) / ((p + pz) * pz) - (1 - p * H0 / h1)
d = sp.simplify(zp - zp_claim)
vals = {px: 0.013, py: -0.021, p: 1.004, P0: p0, H0: h0}
print('  z\' と -(px^2+py^2)/((p+pz)pz) - (1 - p h0/h1) の差 (記号: %s、数値 %.1e)' % (d, abs(float((zp - zp_claim).subs(vals)))))
# dv の閉じた式 (1.9) と 1 - p h0/h1
for dp in (1e-6, 1e-3, -0.02, 0.05):
    ref = 1 - (1 + dp) * h0 / math.sqrt(1 + (p0 * (1 + dp))**2)
    print('  δ = %+.0e: dv = %+.12e (閉じた式 (1.9))、1 - v/v0 の直接の式 %+.12e、差 %.1e' % (dp, dv_closed(dp), ref, abs(dv_closed(dp) - ref)))
print('  δ → 0: dv ≈ -δ/γ0^2 = %+.6e (δ = 1e-3 のとき)、閉じた式 %+.6e' % (-1e-3 / h0**2, dv_closed(1e-3)))
dvdp_sym = sp.simplify(sp.diff(1 - p * H0 / h1, p))
print('  ∂dv/∂δ = -h0/h1^3 の確認 (差 %.1e。Fortran の dvdp = h0/h1^3 は -∂dv/∂δ = ∂(v/v0)/∂δ)' %
      abs(float((dvdp_sym + H0 / h1**3).subs({p: 1.003, P0: p0, H0: h0}))))

print('\n== (b) 閉じた式・DOP853・pysad の追跡の比較 (L = 1 m、3 GeV) ==')
Ld = 1.0


def exact_map(c, L):
    """x += px L/pz、y += py L/pz、z += -(px^2+py^2) L/((p+pz) pz) - dv L"""
    x_, px_, y_, py_, z_, dp = c
    pr = 1 + dp
    pz_ = math.sqrt(pr * pr - px_ * px_ - py_ * py_)
    a = px_ * px_ + py_ * py_
    return np.array([x_ + px_ * L / pz_, px_, y_ + py_ * L / pz_, py_, z_ - a / ((pr + pz_) * pz_) * L - dv_val(dp) * L, dp])


def dv_val(dp):
    return dv_closed(dp)


def ode(c, L):
    def rhs(s, v):
        pr = 1 + c[5]
        pz_ = math.sqrt(pr * pr - c[1]**2 - c[3]**2)
        return [c[1] / pz_, c[3] / pz_, -pr / pz_ + pr * h0 / math.sqrt(1 + (p0 * pr)**2)]      # z' = ∂H/∂p = -p/pz + v/v0
    sol = solve_ivp(rhs, (0, L), [c[0], c[2], c[4]], method='DOP853', rtol=1e-13, atol=1e-16)
    return np.array([sol.y[0, -1], c[1], sol.y[1, -1], c[3], sol.y[2, -1], c[5]])


ffs, el = C.make_drift(Ld)
cases = (np.array([0.001, 0.01, -0.002, 0.02, 0.0, 0.003]),
         np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.01]),
         np.array([0.0, 0.1, 0.0, -0.05, 0.0, -0.02]),
         np.array([0.002, 0.3, 0.0, 0.2, 0.001, 0.0]))
for c in cases:
    f = exact_map(c, Ld)
    o = ode(c, Ld)
    tr = C.track_one(ffs, c)
    print('  x, px, y, py, z, δ = %s' % np.array2string(c, precision=4))
    print('    閉じた式 - DOP853: %.1e、閉じた式 - pysad 追跡: %.1e    (Δz = %+.6e)' % (np.abs(f - o).max(), np.abs(f - tr).max(), f[4] - c[4]))

print('\n== (c) 厳密なドリフト = 近軸 (H2、K = 0) × 運動学的ドリフト D (ΔH) ==')
print('  H2 (K=0) の流れ: x += q L、z += -(q^2/2 + dv) L (q = px/p)。D の流れ (04 章 (4.16)): x += q (1/ζ - 1) L、z += -(2+ζ)(1-ζ)^2/(2ζ) L (ζ = sqrt(1-q^2))')
worst = 0.0
for c in cases:
    x_, px_, y_, py_, z_, dp = c
    pr = 1 + dp
    q, r = px_ / pr, py_ / pr
    zeta = math.sqrt(1 - q * q - r * r)
    para = np.array([x_ + q * Ld, px_, y_ + r * Ld, py_, z_ - (0.5 * (q * q + r * r) + dv_val(dp)) * Ld, dp])
    dmap = np.array([para[0] + q * (1 / zeta - 1) * Ld, px_, para[2] + r * (1 / zeta - 1) * Ld, py_,
                     para[4] - (2 + zeta) * (1 - zeta)**2 / (2 * zeta) * Ld, dp])
    # 逆の順序 (D → H2) でも同じ (どちらも運動量だけの関数)
    dmap2 = np.array([x_ + q * (1 / zeta - 1) * Ld + q * Ld, px_, y_ + r * (1 / zeta - 1) * Ld + r * Ld, py_,
                      z_ - (2 + zeta) * (1 - zeta)**2 / (2 * zeta) * Ld - (0.5 * (q * q + r * r) + dv_val(dp)) * Ld, dp])
    f = exact_map(c, Ld)
    worst = max(worst, np.abs(dmap - f).max(), np.abs(dmap2 - dmap).max())
    print('  q = %.3f, r = %.3f: H2 → D と厳密な式の差 %.1e、D → H2 との差 %.1e、D の寄与 Δx = %+.3e、Δz = %+.3e' %
          (q, r, np.abs(dmap - f).max(), np.abs(dmap2 - dmap).max(), dmap[0] - para[0], dmap[4] - para[4]))
print('  最大差 %.2e' % worst)

print('\n== (d) DISKIN (近軸のドリフト) ==')
print('  Fortran (tturn.f:296): DISKIN ≠ 0 の DRIFT は x += q L、y += r L、z -= ((q^2 + r^2)/2 + dv) L (q, r は規格化した運動量)。H2 (K = 0) の流れと同じ式')
ffs_k, el_k = C.make_drift(Ld, 'DISKIN=1')
for c in cases[2:]:
    x_, px_, y_, py_, z_, dp = c
    pr = 1 + dp
    q, r = px_ / pr, py_ / pr
    para = np.array([x_ + q * Ld, px_, y_ + r * Ld, py_, z_ - (0.5 * (q * q + r * r) + dv_val(dp)) * Ld, dp])
    f = exact_map(c, Ld)
    tr0 = C.track_one(ffs, c)
    tr1 = C.track_one(ffs_k, c)
    print('  q = %.3f: pysad の DISKIN=1 の追跡 - 厳密な式 %.1e、- 近軸の式 %.1e、(DISKIN なし - DISKIN=1) %.1e。近軸と厳密の差 %.1e' %
          (q, np.abs(tr1 - f).max(), np.abs(tr1 - para).max(), np.abs(tr1 - tr0).max(), np.abs(para - f).max()))
print('  → pysad の追跡は DISKIN を読むが使わず、DISKIN=1 でも厳密な式のまま (未実装)。Fortran は上の近軸式')
