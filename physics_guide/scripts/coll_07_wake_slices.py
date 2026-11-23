# -*- coding: utf-8 -*-
r"""ウェイクの切片の方法 (txwake の直訳 wake.WakeKick) の粗視化の誤差の粒子数依存 (17_collective.md §6.3)。

coll_01 (c) と同じ設定 (ガウス分布 σz = 5 mm、W(s) = W0 exp(-s/l)、l = 5 mm、表は -0.05 m から 0.2 m、1.7 GeV) で、粒子数 np を変える。
同じ切片の内部の対は、表が負の z から始まるとき (この表) 全重み W(0) で数えられるので、解析式 (17.12) に対する過大評価は 1/(2 n_s k)
(n_s = 切片の数 ~ √np、k = (1/2) e^{σ²/l²} erfc(σ/l)) 程度のはず。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
from pysad import wake as W

mc2 = 0.51099895e6
ELEMCH = 1.602176634e-19
p0w = 1.7e9 / mc2
g0w = math.sqrt(1 + p0w**2)
sigz, W0, Nb = 5e-3, 1e12, 1e10
for ell in (2e-3, 5e-3, 20e-3):
    k = 0.5 * math.exp(sigz**2 / ell**2) * math.erfc(sigz / ell)
    exp_d = -(ELEMCH * Nb * W0 * k / mc2) / p0w
    print('\nl = %.0f mm: k = %.5f' % (ell * 1e3, k))
    print('  np      切片数 n_s   <Δδ>/解析式 - 1    1/(2 n_s k)    比')
    for npart in (2500, 10000, 40000, 160000):
        rat = []
        for seed in (1, 2, 3):
            rng = np.random.default_rng(seed)
            zz = rng.normal(0.0, sigz, npart)
            zt = np.concatenate([[-0.05, -1e-5], np.linspace(0.0, 0.2, 4001), [0.21]])
            wt = np.where(zt >= 0.0, W0 * np.exp(-zt / ell), 0.0)
            tab = np.stack([zt, wt], axis=1)
            S = np.zeros((npart, 16))
            S[:, 4] = zz
            wk = W.WakeKick(p0w, g0w, W.kick_factor(-1, Nb, 1, npart), 1)
            wk.apply(S, npart, tab, None, 0.0, 0.0, True)
            wk.apply(S, npart, tab, None, 0.0, 0.0, False)
            rat.append(S[:, 5].mean() / exp_d - 1)
        ns = int(math.sqrt(npart))
        pred = 1 / (2 * ns * k)
        print('  %-7d %-10d %+.4f (3 種の平均; 振れ ±%.4f)  %.4f    %.2f' % (npart, ns, np.mean(rat), np.std(rat), pred, np.mean(rat) / pred))
