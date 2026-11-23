# -*- coding: utf-8 -*-
"""相対運動量の縦成分 q_s の効果を入れた 3 次元の積分 (touschek_bunch.qs_ratio) と、実験室系の対のモンテカルロ (touschek_bunch_mc.py の出力) の比較。
3 次元の絶対値 / kick の損失率の比 (δ_c ~ N(0, σδ²/2) の平均と q_s を含む) を、スライスと受容 a ごとに出し、MC の「絶対値 / kick − 1」(touschek_bunch_mc_<名前>.out.txt) と並べる。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_bunch_qs_3d.py astrid2|kekls   (出力は標準出力。保存は touschek_bunch_qs_3d_<名前>.out.txt)"""
import os, re, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
name = sys.argv[1]
from touschek_formula_scan import gather
from pysad.touschek import piwinski_slice_rate
from pysad.touschek_bunch import qs_ratio, piwinski_slice_rate_bunch
tt, tt2, ref = gather(name, 1e9, 0.01)
info = tt2['rest_info']; rest = tt2['rest']
mc = {}
lines = open(os.path.join(HERE, 'touschek_bunch_mc_%s.out.txt' % name), encoding='utf-8').read().split('\n')
for ln in lines:
    m = re.match(r'\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+\S+\s+\S+\s+\S+\s+([+-][0-9.]+) \(([0-9.]+)\)\s+([+-][0-9.]+)', ln)
    if m:
        mc[(int(m.group(1)), float(m.group(3)))] = (float(m.group(5)), float(m.group(6)))
print('%s: σδ = %.3e' % (name, ref['SIGE']))
print(' slice    a    a/σδ   MC abs/kick−1 (±)    3D (δ_c + q_s)−1 (±)   δ_c だけ (1D)−1   時間')
for (j, a), (inc, err) in sorted(mc.items()):
    s = rest[j]
    t = time.time()
    r, e = qs_ratio(s, a, a, nu=2048)
    r0 = piwinski_slice_rate(s, a, a, info['circ'], 1e9)
    rc = piwinski_slice_rate_bunch(s, a, a, info['circ'], 1e9) / r0
    print(' %5d %.4f %5.1f   %+.4f (%.4f)      %+.4f (%.4f)        %+.4f        %.0f s' % (j, a, a / ref['SIGE'], inc, err, r - 1, e, rc - 1, time.time() - t), flush=True)
