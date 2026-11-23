# -*- coding: utf-8 -*-
"""docs/18: 1 GeV の偏向磁石 (L=1.2 m、ANGLE=0.2618) の RADLIGHT 軌道から RadiationField と RadiationSpectrum で求めた
軸上の d²I/dωdΩ と、Jackson の式 (K_{2/3}) の比較。計算の手順は pysad/tests/script_tests/test_v2_C3_radfield.py と同じ。
実行: py -3.12 tools/figures/make_18_radlight.py"""
import math, io
from _common import *
sys.path.insert(0, str(ROOT / "pysad" / "tests"))
import numpy as np
from scipy.special import kv
import _c3_indep as C
from pysad import ffscmd
from pysad.radfield import radiation_field, radiation_spectrum
src = '''MOMENTUM = 1e9;
DRIFT D1 = (L=0.5);
BEND B = (L=1.2 ANGLE=0.2617993877991494);
MARK ST = ();
LINE RING = (ST D1 B D1);
FFS;
USE RING; INS; CALC;
RADLIGHT;
np=1;
x=Table[0,{7},{np}];
Do[x[[7,i]]=1,{i,np}];
d={1,x};
r=TrackParticles[d,4];
traj=r[[2,1]];
Exit;
'''
ffs, I = ffscmd.run_sad_string(src, out=io.StringIO())
traj = np.array(I.globals['traj'], float)
gamma = math.hypot(1.0, ffs.momentum() / 0.51099895e6)
rho = 1.2 / 0.2617993877991494
D = 5000.0; phi = 0.13089969389957471
f = radiation_field(traj.T, [D * math.cos(phi), D * math.sin(phi), 0.0])
tau = f[:, 0]
wc = 3 * gamma ** 3 * C.C / (2 * rho); lc = 2 * math.pi * C.C / wc
tot = 0.0
for comp in (1, 2, 3):
    k, c, s = radiation_spectrum(tau, f[:, comp], [0.2 * lc, 6 * lc, 0.2 * lc])
    tot = tot + np.array(c) ** 2 + np.array(s) ** 2
k = np.array(k); w = k * C.C
i_num = D ** 2 * C.EPS0 * C.C * tot / math.pi
q2 = C.E0 ** 2 / (4 * math.pi * C.EPS0)
def jackson(w):
    xi = w * rho / (3 * C.C * gamma ** 3)
    return q2 / (3 * math.pi ** 2 * C.C) * (w * rho / C.C) ** 2 * gamma ** -4 * kv(2.0 / 3.0, xi) ** 2
ww = np.linspace(0.05, 6, 400) * wc
print("lambda_c [nm]", lc * 1e9, "max rel diff", np.abs(i_num / jackson(w) - 1).max(), "n", traj.shape)
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.semilogy(ww / wc, jackson(ww), color=BLUE, label="Jackson の式 (軸上、K₂/₃)" )
ax.semilogy(w / wc, i_num, "o", ms=3.5, color=ORANGE, label="RadiationField + RadiationSpectrum")
ax.set_xlabel("ω / ωc (ωc = 3γ³c/2ρ、λc = %.2f nm)" % (lc * 1e9)); ax.set_ylabel("d²I/dωdΩ [J s/sr]")
ax.set_title("1 GeV 偏向磁石 (ρ = %.2f m): 軸上のスペクトル" % rho, fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); save(fig, "18_radlight_bend_spectrum.png")
