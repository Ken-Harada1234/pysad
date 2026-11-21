# -*- coding: utf-8 -*-
"""v07_14 (5 GeV) の六極の整理 (R7、検討の 1): ファミリー・台数・配置・K2L・ηx と β、八極・制約、鏡面対称の軸、`cell calc emit` (リング全周、RF・放射 ON)、クロマティシティ。
使い方: py -3.12 sextupole_ring_summary.py [--ring RING]   出力: results/ring_summary.out.txt
"""
import argparse
import io
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import resdrive as R
from pysad import sxtopt as S
from pysad.trksetup import TrackSetup

ap = argparse.ArgumentParser()
ap.add_argument('--ring', default='RING')
args = ap.parse_args()
LOG = io.StringIO()


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')


src, lat = C.load_lattice(args.ring)
els = lat.elements
tw = lat.twiss_table()
cnt = Counter(type(e).__name__ for e in els)
L = sum(float(getattr(e, 'L', 0) or 0) for e in els)
P('リング %s: %d 素子、周長 %.4f m、%s' % (args.ring, len(els), L, dict(cnt)))
oct_k3 = sorted({float(e.K3) for e in els if type(e).__name__ == 'Octupole'})
P('八極 OC1 %d 台の K3 = %s (0 なので 4 次の項への寄与なし。プレースホルダ)' % (cnt['Octupole'], oct_k3))
names, fam, k = S.sextupole_families(lat)
P('\n六極 %d 台、ファミリー %d (セルの種類 LH/MH/LHA/MHA × S1〜S5)。K2L は SAD の K2 (積分強さ、長さ 0.1 m)' % (sum(len(v) for v in fam.values()), len(names)))
P('| ファミリー | 台数 | K2L [1/m²] | ηx [m] (最初の 1 台) | βx [m] | βy [m] | s の最初の位置 [m] |')
P('|---|---|---|---|---|---|---|')
for n, kk in zip(names, k):
    i = fam[n][0]
    r = tw[i]
    P('| %s | %d | %+.3f | %.4f | %.2f | %.2f | %.1f |' % (n, len(fam[n]), kk, r.etx, r.betx, r.bety, r.s))
P('\n分散のない場所の六極 (|ηx| < 1e-4 m): %s (クロマティシティに寄与しない調和の六極)' % [n for n in names if abs(tw[fam[n][0]].etx) < 1e-4])
axes = C.mirror_axes(lat)
P('鏡面対称の軸 (素子の添字): %s、s = %s m、軸での α = %s' % (axes[:2], [round(float(tw[i].s), 3) for i in axes[:2]], [(round(float(tw[i].alfx), 12), round(float(tw[i].alfy), 12)) for i in axes[:2]]))
P('K2L の自由度: ファミリー 20。同じセルの種類の六極は同じファミリーなので、鏡面の像の六極は同じファミリー (鏡面対称の拘束で自由度は減らない)。')
P('等式の条件 (線形): ξx・ξy の 2 本。1 次の RDT は、鏡面の軸から始めると h e^{-iπΔν} が実数なので、5 項 = 実数 5 本、δ の項 3 本 → 線形の条件は計 10 本、残りの自由度 10 が 2 次以上の量 (ADTS 3、4 次の調和) に使える。')
P('\nチューン (CALC): (%.5f, %.5f)' % (tw[-1].mux / 2 / np.pi, tw[-1].muy / 2 / np.pi))
num = R.chromaticity_numeric(lat, order=3, dmax=2e-3)
m = S.SxtModel(lat, C.P0, quad_forms=False)
P('クロマティシティ ξ: CALC (δ 依存) = (%.4f, %.4f)、∂²ν/∂δ² = (%.0f, %.0f)、解析式 (自然 %.3f, %.3f + 偏向磁石の δ の項 + 六極) = (%.4f, %.4f)' % (
    *num['d1'], *num['d2'], *R.chromaticity_natural(lat), *m.xi(m.k_ref)))
st, ffs = TrackSetup.from_source(src, args.ring, rad=True, rfsw=True, alost=0.02)
e = ffs.emit
P('\ncell calc emit (リング全周、RF 6.5 MV・放射 ON): エネルギー損失 U0 = %.3f MeV、運動量圧縮 α = %.3e、RF バケット高さ %.2f%%、シンクロトロンチューン |νs| = %.5f、' % (
    e['energy_loss_eV'] / 1e6, e['momentum_compaction'], 100 * e['bucket_height'], abs(e['tunes'][2])))
P('  emittance_x = %s、emittance_y = %s、momentum_spread = %s' % (e.get('emittance_x'), e.get('emittance_y'), e.get('momentum_spread')))
open(os.path.join(C.OUT, 'ring_summary.out.txt'), 'w', encoding='utf-8').write(LOG.getvalue())
