# -*- coding: utf-8 -*-
"""入門 §6(a): NAFF と単純な FFT の精度。FODO リング (3 GeV、105.6 m) の微小振幅 (J = 1e-14 m rad。ADTS による ν のずれ ~ 4e-12 は無視できる) の
1 粒子を追跡し、周回数 N を変えて、(1) FFT のピークのビン、(2) Hann 窓つきの FFT のピークのビン、(3) NAFF (窓なし)、(4) NAFF (Hann p=1)、(5) NAFF (Hann p=2) の
チューンの誤差 |ν − ν(CALC)| を比べる。基準の ν は CALC の線形のチューン (6.696329...)。信号は z = x̂ − i p̂ (複素)。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/primer/primer_naff.py   出力: primer_naff.png、標準出力
"""
import math
import numpy as np
import primer_common as C
plt = C.setup_plot()
from pysad import adtsdirect as D, naff

lat, rows, tw, nu, circ = C.fodo()
(ax, bx), (ay, by) = tw
NMAX = 4096
J = 1e-14
c = D.initial_coords(tw, [J], [J], [0.3], [0.9])
rec = D.track_record(lat.elements, C.P_FODO, c, NMAX)
x, px = rec[:, 0, 0], rec[:, 0, 1]
xh, ph = naff.twiss_normalize(x, px, ax, bx)
z_all = (xh - xh.mean()) - 1j * (ph - ph.mean())
nu_true = float(D.fr(nu[0]))          # 符号つき端数 (−0.303671...。z = x̂ − i p̂ の向きで ν = 0.696329 − 1)
print('FODO リング (3 GeV、105.6 m)、J = %.0e m rad、基準 ν = CALC の線形のチューンの端数 %.9f (= 0.696329 − 1)' % (J, nu_true))
Ns = [32, 64, 128, 256, 512, 1024, 2048, 4096]
meth = ['FFT のピークのビン (窓なし)', 'FFT のピークのビン (Hann)', 'NAFF (窓なし)', 'NAFF (Hann p=1)', 'NAFF (Hann p=2)']
err = np.zeros((len(meth), len(Ns)))
for k, N in enumerate(Ns):
    z = z_all[:N] - 0   # 先頭から N 周
    n = np.arange(N)
    f = np.fft.fftfreq(N)
    w1 = 0.5 * (1 - np.cos(2 * np.pi * (n + 0.5) / N))
    for i, w in enumerate((np.ones(N), w1)):
        sp = np.abs(np.fft.fft(w * z))
        err[i, k] = abs(f[int(np.argmax(sp))] - nu_true)
    for i, p in zip((2, 3, 4), (0, 1, 2)):
        fr_, _ = naff.naff(z, 1, window=p, fmin=nu_true - 0.1, fmax=nu_true + 0.1)
        err[i, k] = abs(fr_[0] - nu_true)
print('%-30s' % '誤差 |ν − ν_CALC|' + ''.join('%11d' % N for N in Ns) + '   (周回数 N)')
for i, m in enumerate(meth):
    print('%-30s' % m + ''.join('%11.2e' % e for e in err[i]))
fig, a = plt.subplots(figsize=(7, 4.5))
for i, m in enumerate(meth):
    a.loglog(Ns, np.maximum(err[i], 1e-17), 'o-', label=m, ms=4)
a.loglog(Ns, 0.5 / np.array(Ns), 'k:', lw=0.8)
a.text(Ns[1], 0.5 / Ns[1] * 1.3, '1/(2N)', fontsize=8)
a.set_xlabel('周回数 N'); a.set_ylabel('チューンの誤差 |ν − ν(CALC)|')
a.set_title('FODO リング、J = 1e-14 m rad: チューンの精度', fontsize=10)
a.legend(fontsize=7.5); a.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(C.os.path.join(C.HERE, 'primer_naff.png'), dpi=110)
