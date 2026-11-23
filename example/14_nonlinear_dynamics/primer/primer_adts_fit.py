# -*- coding: utf-8 -*-
"""入門 §6(a): 追跡 + NAFF + 多項式フィットで ADTS (∂ν/∂J) を求める手順の中身。FODO リング (3 GeV、105.6 m、6 極 K2L = 4, −6、リング先頭から)。
J の格子 (Jx, Jy ∈ {0, 1/5, ..., 1}·Jmax の 35 点 × (φ, φ+π) の対 = 70 粒子) を、強さを変えずに 2048 周追跡し (pysad k64、RF なし、放射なし、δ = 0)、
各粒子の NAFF (Hann p=1、線形のチューン ± 0.02 を探す) のチューンを ±対で平均し、ν(Jx, Jy) を J の全次数 ≤ order の多項式でフィットする。
  (1) 図左: Jy = 0 の点の νx − ν0 と、フィット (order = 1, 2) の曲線。Jmax = 1e-5 (J が大きいと 1 次のフィットは外れる)
  (2) 表: Jmax = 1e-7, 1e-6, 1e-5 と order = 1, 2, 3 の ADTS を解析の合計 (adts().total) と比べる
  (3) 表: 周回数 N = 128 … 2048 の ADTS (Jmax = 1e-7、order = 2)
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/primer/primer_adts_fit.py   出力: primer_adts_fit.png、標準出力  (所要 約 1 分)
"""
import math
import numpy as np
import primer_common as C
plt = C.setup_plot()
from pysad import adtsdirect as D, resdrive as R

lat, rows, tw, nu, circ = C.fodo()
tw_, orb, g, nu_ = D.lattice_start(lat, 0.0, rows)
A = R.adts(lat, rows=rows)
ana = np.array([A.total['xx'], A.total['xy'], A.total['yy']])
print('FODO リング (3 GeV、105.6 m)、νx = %.6f、νy = %.6f。解析の合計 (adts().total) xx, xy, yy = %.3f, %.3f, %.3f [1/(m rad)]' % (nu[0], nu[1], *ana))
NT = 2048
data = {}


def track_grid(jmax):
    jx, jy = D.grid_actions(jmax, 5)
    JX, JY = np.r_[jx, jx], np.r_[jy, jy]
    ph = np.r_[np.zeros(len(jx)), np.full(len(jx), math.pi)]
    c = D.initial_coords(tw, JX, JY, ph, ph)
    return jx, jy, D.track_record(lat.elements, C.P_FODO, c, NT)


def fit(jx, jy, rec, order, nturn):
    nuh = D.record_tunes(rec, tw, guess=g, halfwidth=0.02, nturn=nturn)
    nuh = D.pair_average(nuh)
    r = D.adts_from_tunes(jx, jy, nuh, order=order, nu_lin=g)
    return r, nuh


for jmax in (1e-7, 1e-6, 1e-5):
    data[jmax] = track_grid(jmax)
print()
print('(2) J の範囲とフィットの次数 (1024 周)。相対差 = (フィット − 解析)/|解析|')
print('%-9s %-6s %12s %12s %12s    %10s %10s %10s' % ('Jmax', 'order', 'xx', 'xy', 'yy', 'dxx', 'dxy', 'dyy'))
for jmax in (1e-7, 1e-6, 1e-5):
    jx, jy, rec = data[jmax]
    for order in (1, 2, 3):
        r, _ = fit(jx, jy, rec, order, 1024)
        v = np.array([r['xx'], 0.5 * (r['xy'] + r['yx']), r['yy']])
        d = (v - ana) / np.abs(ana)
        print('%-9.0e %-6d %12.3f %12.3f %12.3f    %+9.2f%% %+9.2f%% %+9.2f%%' % (jmax, order, *v, *(100 * d)))
print()
print('(3) 周回数 N (Jmax = 1e-7、order = 2)')
jx, jy, rec = data[1e-7]
for N in (128, 256, 512, 1024, 2048):
    r, _ = fit(jx, jy, rec, 2, N)
    v = np.array([r['xx'], 0.5 * (r['xy'] + r['yx']), r['yy']])
    print('N = %4d  xx, xy, yy = %.3f, %.3f, %.3f   相対差 %+.3f%%, %+.3f%%, %+.3f%%   フィットの残差 rms = %.1e' % (N, *v, *(100 * (v - ana) / np.abs(ana)), r['fitx']['rms']))

# 図
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
jx, jy, rec = data[1e-5]
r1, nuh = fit(jx, jy, rec, 1, 1024)
r2, _ = fit(jx, jy, rec, 2, 1024)
m = (jy == 0)
nux = D.fr(nuh[m, 0] - g[0])
a = axs[0]
a.plot(jx[m] * 1e6, nux / jx[m], 'ko', label='追跡 + NAFF (Jy = 0) の (νx − ν0)/Jx')
jj = np.linspace(jx[m].min(), jx.max(), 100)
a.plot(jj * 1e6, np.full_like(jj, r1['fitx']['c'][(1, 0)]), 'C0--', label='1 次のフィット: ∂νx/∂Jx = %.1f' % r1['xx'])
a.plot(jj * 1e6, r2['fitx']['c'][(1, 0)] + r2['fitx']['c'][(2, 0)] * jj, 'C3-', label='2 次のフィット: ∂νx/∂Jx = %.1f' % r2['xx'])
a.set_xlabel(r'$J_x$ [$10^{-6}$ m rad]'); a.set_ylabel(r'$(\nu_x-\nu_0)/J_x$ [1/(m rad)]')
a.set_title('FODO、Jmax = 1e-5: J の 2 次の項で 1 次のフィットは外れる', fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)
a = axs[1]
Ns = [128, 256, 512, 1024, 2048]
jx, jy, rec = data[1e-7]
dd = []
for N in Ns:
    r, _ = fit(jx, jy, rec, 2, N)
    v = np.array([r['xx'], 0.5 * (r['xy'] + r['yx']), r['yy']])
    dd.append(np.abs((v - ana) / ana))
dd = np.array(dd)
for k, nm in enumerate(('xx', 'xy', 'yy')):
    a.semilogy(Ns, dd[:, k], 'o-', label=nm)
a.set_xscale('log', base=2)
a.set_xlabel('周回数 N'); a.set_ylabel('|フィット − 解析| / |解析|')
a.set_title('Jmax = 1e-7、order = 2: 周回数への依存', fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(C.os.path.join(C.HERE, 'primer_adts_fit.png'), dpi=110)
