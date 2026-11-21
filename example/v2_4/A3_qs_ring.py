# -*- coding: utf-8 -*-
"""v2.4.0 A10: ASTRID2 全周の Touschek の損失率 (一様な受容 a) の、δ_c = 0 (elegant・K64・Piwinski と同じ基準) に対する増え: δ_c だけ (1 次元) と δ_c + q_s (3 次元を 6 / 12 / 24 スライスで補間)。
使い方: PYTHONUTF8=1 py -3.12 example/v2_4/A3_qs_ring.py   (出力 A3_qs_ring.out.txt)"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'example', 'v2_3'))
from T2_scan import gather
from pysad.touschek_bunch import touschek_piwinski_bunch
tt, tt2, ref = gather('astrid2', 1e9, 0.01)
print('ASTRID2、IBS 込み (PBUNCH 1e9、MINCOUP 0.01)、σδ = %.3e' % ref['SIGE'])
for a in (0.004, 0.006, 0.010):
    r1 = touschek_piwinski_bunch(tt2, a=a)
    print('a = %.3f (a/σδ = %.1f): δ_c = 0 の寿命 %.1f s、δ_c だけ +%.2f%% (寿命 %.1f s)' % (a, a / ref['SIGE'], 1 / r1['rate0'], 100 * r1['increase'], r1['lifetime']), flush=True)
    for ns in (6, 12, 24):
        t = time.time()
        r = touschek_piwinski_bunch(tt2, a=a, qs=True, qs_nsample=ns)
        print('   qs=True (%d スライス): +%.2f%% (寿命 %.1f s)、補正の係数 %.4f〜%.4f、%.0f s' % (ns, 100 * r['increase'], r['lifetime'], min(r['qs']['factor']), max(r['qs']['factor']), time.time() - t), flush=True)
