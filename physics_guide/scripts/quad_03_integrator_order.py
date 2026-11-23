# -*- coding: utf-8 -*-
"""四極の本体の積分器: K64 の Strang 分割 (2 次) と、v2 の 4 次の Yoshida 合成 (04_quad.md 式 (4.18)〜(4.23)) の誤差の次数を確かめる。

(a) Yoshida の係数 w1 = 1/(2 - 2^(1/3))、w0 = 1 - 2 w1 が 2 w1 + w0 = 1、2 w1^3 + w0^3 = 0 を満たすこと。
(b) 運動学的な運動方程式 (x' = px/pz, px' = -k x, ...) を DOP853 (rtol 1e-13) で積分した厳密な写像を基準に、
    刻みの数 nd = 1, 2, 4, ... の誤差。2 倍にしたときの誤差の比は Strang が 4、Yoshida が 16。
(c) 同じ計算量 (本体 B の評価の数) あたりの誤差。刻みの数の式 (ndiv) の表。
図: physics_guide/figures/quad_integrator_order.png
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6, quadbody as qb

W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
W0 = 1.0 - 2.0 * W1
print('== (a) Yoshida の係数 ==')
print('  w1 = 1/(2 - 2^(1/3)) = %.15f、w0 = 1 - 2 w1 = %.15f' % (W1, W0))
print('  2 w1 + w0 - 1 = %.2e (= 0)、2 w1^3 + w0^3 = %.2e (= 0)' % (2 * W1 + W0 - 1.0, 2 * W1**3 + W0**3))
print('  pysad の係数 quadbody.QB_BW[1] = (%.15f, %.15f, %.15f)' % qb.QB_BW[1])


def ref_map(c, L, k):
    """運動学的な運動方程式の厳密な写像 (δ = 0、dv = 0)。k = K1/L (x に収束で正)。戻り: x, px, y, py, z"""
    def rhs(s, v):
        x, px, y, py, z = v
        pz = math.sqrt(1.0 - px * px - py * py)
        return [px / pz, -k * x, py / pz, k * y, -(1.0 / pz - 1.0)]
    return solve_ivp(rhs, (0.0, L), list(c) + [0.0], method='DOP853', rtol=1e-13, atol=1e-16).y[:, -1]


def integrate(c, L, K1L, nd, sch):
    """quadbody の D と B を、刻みの数 nd・方式 sch (0: Strang、1: Yoshida) で組んで 1 台の四極 (K1L > 0、δ = 0、dv = 0) を通す"""
    x, px, y, py = c
    z = 0.0
    akk = K1L / L
    aln = L / nd
    pa = t6._tzsetparam0(0.0, aln * qb.QB_BW[sch][0], akk)
    pb = t6._tzsetparam0(0.0, aln * qb.QB_BW[sch][1], akk) if sch == 1 else pa
    m = qb.QB_M[sch]
    for g in range(m * nd):
        dl, bl, jb, _r, _f = qb.group(g, aln, nd, sch)
        x, y, z = t6._qb_drift(x, px, y, py, z, dl)
        x, px, y, py, z = t6._qb_body(x, px, y, py, z, 0.0, pb if jb == 1 else pa, bl, 0.0, 0.0, False)
    dl, _b, _j, _r, _f = qb.group(m * nd, aln, nd, sch)
    x, y, z = t6._qb_drift(x, px, y, py, z, dl)
    return np.array([x, px, y, py, z])


L, K1L = 1.0, 1.5
cases = {'中振幅 (10 mm、10 mrad)': (0.01, 0.01, -0.008, 0.008), '大振幅 (30 mm、30 mrad)': (0.03, 0.03, -0.025, 0.025)}
nds = [1, 2, 4, 8, 16, 32, 64]
res = {}
print('\n== (b) 刻みの数 nd と誤差 (四極 L = %g m、K1L = %g、最大差は x, px, y, py, z の中の最大) ==' % (L, K1L))
for name, c in cases.items():
    ref = ref_map(c, L, K1L / L)
    print(' ', name)
    print('    nd    Strang 誤差    比      Yoshida 誤差   比')
    errs = {0: [], 1: []}
    for nd in nds:
        for sch in (0, 1):
            errs[sch].append(np.abs(integrate(c, L, K1L, nd, sch) - ref).max())
    for i, nd in enumerate(nds):
        r0 = ('%6.2f' % (errs[0][i - 1] / errs[0][i])) if i else '     -'
        r1 = ('%6.2f' % (errs[1][i - 1] / errs[1][i])) if i else '     -'
        print('   %3d   %.3e   %s   %.3e   %s' % (nd, errs[0][i], r0, errs[1][i], r1))
    res[name] = errs

print('\n== (c) 計算量あたりの誤差 (本体 B の評価の数: Strang は nd、Yoshida は 3 nd) と、刻みの数の式 ==')
c = cases['中振幅 (10 mm、10 mrad)']
print('  追跡の刻みの数: K64 は 1 + int(L K1/0.2)、v2 は 1 + int(L K1/0.1) (K1 は SAD の積分強度、EPS = 1)')
for sch, nd_ in ((0, qb.ndiv(L * K1L, 0.2, 0)), (1, qb.ndiv(L * K1L, 0.2, 1))):
    e = np.abs(integrate(c, L, K1L, nd_, sch) - ref_map(c, L, K1L / L)).max()
    print('   %s: nd = %2d、B の評価 %3d 回、誤差 %.2e' % ('Strang  (K64)' if sch == 0 else 'Yoshida (v2) ', nd_, nd_ * qb.QB_M[sch], e))
print('\n  刻みの数の表 (EPS = 1): ヘルプの式 1 + floor(10 L K1 / EPS)、追跡・CALC の K64 と v2')
print('    L [m]   K1 (積分強度)   L K1     ヘルプ   追跡 K64   追跡 v2   CALC K64   CALC v2')
for Lq, K1q in ((0.3, 0.6), (0.5, 1.5), (1.0, 1.5), (0.5, 4.0), (2.0, 3.0)):
    a = Lq * K1q
    print('   %5.2f   %9.2f     %6.3f   %5d   %7d   %7d   %8d   %7d' % (
        Lq, K1q, a, 1 + int(10 * a), qb.ndiv(a, 0.2, 0), qb.ndiv(a, 0.2, 1), qb.ndiv_calc(a, 0), qb.ndiv_calc(a, 1)))

# ---------------------------------------------------------------- 図
matplotlib = C.font()
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.4))
cols = {'Strang': '#1f77b4', 'Yoshida': '#d62728'}
name = '中振幅 (10 mm、10 mrad)'
e0, e1 = res[name][0], res[name][1]
ax[0].loglog(nds, e0, 'o-', color=cols['Strang'], label='Strang (K64、2 次)')
ax[0].loglog(nds, e1, 's-', color=cols['Yoshida'], label='Yoshida (v2、4 次)')
n_ = np.array([1, 64.0])
ax[0].loglog(n_, e0[0] * n_**-2.0, ':', color=cols['Strang'], lw=1, label='傾き -2')
ax[0].loglog(n_, e1[0] * n_**-4.0, ':', color=cols['Yoshida'], lw=1, label='傾き -4')
ax[0].set_xlabel('刻みの数 nd'); ax[0].set_ylabel('厳密な写像との最大差 [m または rad]')
ax[0].set_title('(a) 刻みの数と誤差\n(L = 1 m、K1L = 1.5、10 mm・10 mrad)'); ax[0].legend(fontsize=8); ax[0].grid(True, which='both', alpha=0.3)
ax[1].loglog([n * qb.QB_M[0] for n in nds], e0, 'o-', color=cols['Strang'], label='Strang (K64)')
ax[1].loglog([n * qb.QB_M[1] for n in nds], e1, 's-', color=cols['Yoshida'], label='Yoshida (v2)')
ax[1].set_xlabel('本体 B の評価の数 (計算量)'); ax[1].set_ylabel('最大差')
ax[1].set_title('(b) 計算量あたりの誤差'); ax[1].legend(fontsize=8); ax[1].grid(True, which='both', alpha=0.3)
# (c) 1 刻み h の構成: D (ΔH の流れ) と B (H2 の流れ) を、それぞれの「時間」(適用した長さの累計、負も数える) の軸に、適用の順に 1 行ずつ矢印で描く
row = 0.0
for sch, lab in ((0, 'Strang (K64)'), (1, 'Yoshida (v2)')):
    bw = qb.QB_BW[sch]
    m = qb.QB_M[sch]
    d = [bw[0] / 2] + [(bw[j] + bw[j - 1]) / 2 for j in range(1, m)] + [bw[m - 1] / 2]
    ax[2].text(-0.55, row, lab, fontsize=9, va='center', ha='left', fontweight='bold')
    row -= 0.8
    tD = tB = 0.0
    k = 1
    for j in range(m + 1):
        for kind, ln in (('D', d[j]), ('B', bw[j] if j < m else None)):
            if ln is None:
                continue
            t0 = tD if kind == 'D' else tB
            col = '#1f77b4' if kind == 'D' else '#d62728'
            ax[2].annotate('', xy=(t0 + ln, row), xytext=(t0, row), arrowprops=dict(arrowstyle='->', lw=2.2, color=col))
            ax[2].text(-0.55, row, '%d %s' % (k, kind), fontsize=8, va='center', color=col)
            ax[2].text(max(t0, t0 + ln) + 0.04, row, '%+.3f h' % ln, fontsize=7.5, va='center')
            if kind == 'D':
                tD += ln
            else:
                tB += ln
            k += 1
            row -= 0.6
    row -= 0.8
ax[2].set_xlim(-0.6, 1.9)
ax[2].set_ylim(row + 0.3, 0.5)
ax[2].set_yticks([])
ax[2].axvline(0, color='k', lw=0.5)
ax[2].axvline(1, color='k', lw=0.5)
ax[2].set_xlabel('各流れの「時間」(適用した長さの累計 / h)')
ax[2].set_title('(c) 1 刻みの構成 (青 D: ΔH の流れ、赤 B: H2 の厳密な流れ。\n左の数字は適用の順。Yoshida は逆向きの矢印を含む)', fontsize=9)
plt.tight_layout()
os_path = C.FIG_DIR + '/quad_integrator_order.png'
plt.savefig(os_path, dpi=110)
print('\n図を保存: physics_guide/figures/quad_integrator_order.png')
