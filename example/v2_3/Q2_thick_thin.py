# -*- coding: utf-8 -*-
"""v2.3.0 Q2: 厚い 6 極 (SEXT L=0.2 m) と、それを n 分割した薄レンズ列の ADTS・共鳴駆動項の差 (n を増やすと 1/n² で 0 に収束)。
条件: 小さな FODO リング (example/00_lattice/FODO_ring.sad の SF・SD に K2 = 4, −6 [SAD の積分値] を入れたもの。3 GeV、24 セル、105.6 m、
νx = 6.696329、νy = 4.906056。整数・半整数共鳴までの距離 0.09 以上)。6 極の 2 次 + 共鳴駆動項、δ = 0 (分散の項は h11001 など)。
使い方 (リポジトリ直下から): py -3.12 example/v2_3/Q2_thick_thin.py
出力: 厚い式 (method='analytic')、数値積分 (method='numeric'、nslice = n の中点)、分割した薄レンズ列のラティス (slice_thin(n) を解析式で) の 3 通り。"""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, ROOT)
from pysad import resdrive as R, compat
from pysad.parser import load_sad_string
compat.set_mode('v2')
src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
lat = load_sad_string(src).build_lattice('RING')
A = R.adts(lat)
T = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200', 'h11001', 'h00111', 'h20001', 'h00201', 'h10002')
H = R.resonance_terms(lat, T, bend_cubic=False)         # 6 極だけ (R4 以後の既定は偏向磁石の項つき)
print('条件: 3 GeV、FODO 24 セル 105.6 m、SF L=0.2 K2L=4、SD L=0.2 K2L=-6、δ=0、結合なし、6 極の 2 次 (ADTS)・1 次 (h)')
print('厚い式 (解析) ADTS [1/(m rad)]: xx = %.6f  xy = %.6f  yy = %.6f' % (A.sext['xx'], A.sext['xy'], A.sext['yy']))
print()
print('ADTS の差 (分割 − 厚い解析)。slice = 分割した薄レンズ列のラティスを解析式で、numeric = 同じ素子を nslice 中点で数値積分')
print('  n    slice xx       slice xy       slice yy       numeric xx     numeric xy     numeric yy     比 (slice xx, 前の n / この n)')
prev = None
for n in (1, 2, 4, 8, 16, 32, 64):
    S = R.adts(R.slice_thin(lat, n))
    N = R.adts(lat, method='numeric', nslice=n)
    d = [S.sext[k] - A.sext[k] for k in ('xx', 'xy', 'yy')] + [N.sext[k] - A.sext[k] for k in ('xx', 'xy', 'yy')]
    print('%3d  %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e   %s' % ((n,) + tuple(d) + ('%.3f' % (prev / d[0]) if prev else '-',)))
    prev = d[0]
print()
print('共鳴駆動項 |h(薄レンズ n 分割) − h(厚い解析)| (全周の和)')
print('  n   ' + ' '.join('%-9s' % t for t in T))
for n in (1, 2, 4, 8, 16, 32, 64):
    S = R.resonance_terms(R.slice_thin(lat, n), T, bend_cubic=False)
    print('%3d  ' % n + ' '.join('%.2e ' % abs(S.total[nm] - H.total[nm]) for nm in H.names))
print()
print('長さ → 0 の極限: SEXT の長さ L を小さくして K2L を一定にした厚い式の ADTS xx と、長さ 0 の MULT (薄レンズ) の ADTS xx')
thin_src = src.replace('SEXT  SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)', 'MULT  SF = (L=0 K2=0.8) SD = (L=0 K2=-1.2)')
t = R.adts(load_sad_string(thin_src).build_lattice('RING'))
print('  L = 0 (MULT)  xx = %.8f  xy = %.8f  yy = %.8f   ← 薄レンズ (積分を素子の位置での値に置き換えた式)' % (t.sext['xx'], t.sext['xy'], t.sext['yy']))
for L in (0.2, 0.02, 0.002, 2e-4, 2e-5):
    s = src.replace('SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)', 'SF = (L=%g K2=0.8) SD = (L=%g K2=-1.2)' % (L, L))    # K2 は SAD の積分値 (L によらず 0.8, −1.2)
    a = R.adts(load_sad_string(s).build_lattice('RING'))
    print('  L = %-8g xx = %.8f  xy = %.8f  yy = %.8f   (厚い − 薄)xx = %+.2e' % (L, a.sext['xx'], a.sext['xy'], a.sext['yy'], a.sext['xx'] - t.sext['xx']))
print('  (差は L に比例: 薄レンズの値との差は素子が長さ L だけリングを伸ばすことと、素子内の位相の積分による。L → 0 で 0 に収束)')
