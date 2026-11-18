# -*- coding: utf-8 -*-
"""v2next E 項目 3: Polarization の脱偏極率と FLUC 追跡の差 (B3: 4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け。
放射の光子 1 個 = δ の蹴り h。蹴りへのスピンの不変軸の応答 G = ∂n/∂δ を、追跡 (M2: 非線形な軌道 + 厳密なスピン) と線形の軌道モデル (M1: B3 の MC と同じ) で、
同じ初期条件 (ベータトロン振幅だけ / シンクロトロン振幅だけ / 両方、振幅は平衡ビームの f 倍) の組で測り、<|G2|²>/<|G1|²> を出す。
脱偏極率 = Σ_i ½ B_δδ,i |G_i|² なので、この比が追跡と式の比に当たる (補助: pysad/tests/_e_spinkick.py)。
使い方: py -3.12 example/v2_next/E_pol_kick_response.py E[GeV] V[MV] 試料数 K[周] 乱数の種 out.json [modes] [f の列]   (1 試料 約 1〜2 分)
要約: py -3.12 example/v2_next/E_pol_kick_response.py --sum out1.json out2.json ..."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'pysad', 'tests'))
import _e_spinkick as K

if sys.argv[1] == '--sum':
    for fn in sys.argv[2:]:
        d = json.load(open(fn))
        print(fn, 'E', d['E'], 'V', d['V'], 'N', d['N'], '線形の |N⊥δ|', d.get('Nperp_delta'))
        S = d['samples']
        for k in sorted(set((s['mode'], s['f']) for s in S)):
            g1 = np.array([abs(complex(*s['G1'])) ** 2 for s in S if (s['mode'], s['f']) == k])
            g2 = np.array([abs(complex(*s['G2'])) ** 2 for s in S if (s['mode'], s['f']) == k])
            rng = np.random.default_rng(0)
            bs = [g2[i].mean() / g1[i].mean() for i in (rng.integers(0, len(g1), len(g1)) for _ in range(300))]
            print('  %-5s f=%.1f n=%2d  <G1^2>=%.3f <G2^2>=%.3f  比 %.3f +- %.3f' % (k[0], k[1], len(g1), g1.mean(), g2.mean(), g2.mean() / g1.mean(), np.std(bs)))
    sys.exit(0)
E, V, N, Kt, seed, out = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
modes_ = sys.argv[7].split(',') if len(sys.argv) > 7 else ['beta', 'synch', 'both']
fs = [float(x) for x in sys.argv[8].split(',')] if len(sys.argv) > 8 else [1.0]
ring = K.Ring(E, V, 0)
nu = abs(ring.sp['nup'])
secs = K.build_secs(ring)
N0 = np.array(ring.sp['spin_orbit'])
h = 3e-5
res = dict(E=E, V=V, K=Kt, N=N, seed=seed, h=h, nu_spin=nu, Nperp_delta=float(np.linalg.norm(N0[:, 5] - ring.n0 * (N0[:, 5] @ ring.n0))), samples=[])
rng = np.random.default_rng(seed)
for mode in modes_:
    for f in fs:
        for ip in range(N):
            u = ring.sample(mode, f, rng)
            r1 = K.track_pair_linear(ring, secs, u, Kt, [h, -h], nu)
            r2 = K.track_pair(ring, u, Kt, [h, -h], nu)
            G1 = (r1[1][1] - r1[2][1]) / (2 * h)
            G2 = (r2[1][1] - r2[2][1]) / (2 * h)
            res['samples'].append(dict(mode=mode, f=f, G1=[G1.real, G1.imag], G2=[G2.real, G2.imag], u=(u - ring.cod).tolist()))
            print(mode, f, ip, abs(G1), abs(G2), flush=True)
            json.dump(res, open(out, 'w'))
K.restore_cave()
