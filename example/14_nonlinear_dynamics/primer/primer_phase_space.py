# -*- coding: utf-8 -*-
"""入門 §2: FODO リングの 1 粒子の位相空間 (x, px) と正規化座標。
  左: 6 極なし、小振幅 (J = 1e-6 m rad)。リング先頭の (x, px) の 1 周ごとの点は楕円、正規化座標 (x/√β, (αx+βpx)/√β) では円。
  中・右: 6 極あり (K2L = 4, −6)、振幅を大きくすると円が三角形にゆがみ、さらに大きいと外れる (J = 5e-6, 1.5e-5, 3e-5)。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/primer/primer_phase_space.py   出力: primer_phase_space.png、標準出力
"""
import math
import numpy as np
import primer_common as C
plt = C.setup_plot()
from pysad import adtsdirect as D

NT = 2000
lat0, rows0, tw0, nu0, circ = C.fodo(k2=(0.0, 0.0))
lat1, rows1, tw1, nu1, _ = C.fodo()
(ax, bx), (ay, by) = tw1
print('FODO リング: p = %.1f GeV/c、周長 %.3f m、24 セル。リング先頭 βx = %.4f m、αx = %.4f、νx = %.6f、νy = %.6f' % (C.P_FODO / 1e9, circ, bx, ax, nu1[0], nu1[1]))
print('J と振幅: x_max = √(2 J βx)。J = 1e-6 m rad → %.3f mm' % (1e3 * math.sqrt(2e-6 * bx)))
fig, axs = plt.subplots(2, 4, figsize=(15, 7))


def run(lat, tw, J, label, col):
    c = D.initial_coords(tw, [J], [1e-12], [0.0], [0.0])
    rec = D.track_record(lat.elements, C.P_FODO, c, NT)
    x, px = rec[:, 0, 0], rec[:, 0, 1]
    lost = (not np.all(np.isfinite(x))) or np.abs(x).max() > 0.5
    if lost:
        print('  %s: J = %.1e  失われた (|x| > 0.5 m または NaN)' % (label, J))
        return None
    xh, ph = x / math.sqrt(tw[0][1]), (tw[0][0] * x + tw[0][1] * px) / math.sqrt(tw[0][1])
    nux = D.record_tunes(rec, tw, guess=(float(D.fr(nu1[0])), float(D.fr(nu1[1]))), halfwidth=0.1)[0, 0]
    print('  %s: J = %.1e  x_max = %.3f mm  px_max = %.3f mrad  NAFF のチューン端数 νx = %.6f  (線形 %.6f)  正規化の半径の最大/最小 = %.4f' % (
        label, J, 1e3 * np.abs(x).max(), 1e3 * np.abs(px).max(), nux, float(D.fr(nu1[0])), np.hypot(xh, ph).max() / np.hypot(xh, ph).min()))
    return x, px, xh, ph


cases = [(lat0, tw0, 1e-6, '6 極なし', 'C0'), (lat1, tw1, 5e-6, '6 極あり', 'C1'), (lat1, tw1, 1.5e-5, '6 極あり', 'C2'), (lat1, tw1, 3e-5, '6 極あり', 'C3')]
for k, (lat, tw, J, lab, col) in enumerate(cases):
    r = run(lat, tw, J, lab, col)
    for row in (0, 1):
        a = axs[row, k]
        if r is not None:
            x, px, xh, ph = r
            if row == 0:
                a.plot(x * 1e3, px * 1e3, '.', ms=1.5, color=col)
                a.set_xlabel('x [mm]'); a.set_ylabel('px [mrad]')
            else:
                a.plot(xh * 1e3 / math.sqrt(1.0), ph * 1e3, '.', ms=1.5, color=col)
                a.set_xlabel(r'$\hat x = x/\sqrt{\beta}$ [$10^{-3}\sqrt{\mathrm{m}}$]'); a.set_ylabel(r'$\hat p=(\alpha x+\beta p_x)/\sqrt{\beta}$ [$10^{-3}\sqrt{\mathrm{m}}$]')
                a.set_aspect('equal', adjustable='datalim')
        else:
            a.text(0.5, 0.5, 'lost', ha='center', transform=a.transAxes)
        a.set_title('%s, J=%.1e m rad' % (lab, J), fontsize=9)
        a.grid(alpha=0.3)
fig.suptitle(u'FODO リング (3 GeV、105.6 m) 先頭の位相空間: 上 = (x, px)、下 = 正規化座標 (%d 周)' % NT)
fig.tight_layout()
fig.savefig(C.os.path.join(C.HERE, 'primer_phase_space.png'), dpi=110)
