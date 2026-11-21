# -*- coding: utf-8 -*-
"""v2.3.0 S2 (1): KEK-LS (kek-ls_v03_12.lat) の誤差つき (補正なし) の線形の性質: 強さ ×0.03〜×3、種 1〜5 の閉軌道・線形の安定性・結合・フィードダウンの見積もり。
周波数マップを走らせず、fmapcheck.errors_summary だけ (1 本 数 s)。使い方: py -3.12 example/v2_3/S2_kekls_orbit.py"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from pysad import compat, fmapres as FR, fmapcheck as FC
from pysad.trksetup import TrackSetup
from S2_errors_scan import kind_sigmas
compat.set_mode('v2')
st0 = TrackSetup.from_source(open(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), encoding='utf-8').read(), 'RING', rad=False, rfsw=False, alost=0.1)[0]
print('KEK-LS: 誤差 (all: 四極・6 極・偏向磁石・8 極の DX・DY 50 μm、四極のロール 0.2 mrad、偏向磁石のロール 0.1 mrad) × k、補正なし')
print('%6s %4s | %-9s %-7s | %-18s | %-8s %-8s %-8s | %s' % ('k', '種', '閉軌道', '安定', 'rms x,y [mm]', '|c-|予想', 'split予', 'split数', 'ΔQ予/数'))
for k in (0.03, 0.1, 0.3, 1.0, 3.0):
    for seed in range(1, 6):
        st, _ = FR.apply_alignment_errors(st0, seed=seed, streams='type', **kind_sigmas('kekls', 'all', k))
        e = FC.errors_summary(st)
        co = e.get('closed_orbit', True)
        print('%6g %4d | %-9s %-7s | %7.3f %7.3f   | %8.4f %8.4f %8.4f | %.4f / %.4f' % (
            k, seed, 'あり' if co else '求まらない', e['stable'], e['orbit_rms'][0] * 1e3, e['orbit_rms'][1] * 1e3, e['cminus_pred'], e['split_pred'], e['split_num'],
            np.hypot(*e['dnu_pred']), np.hypot(*e['dnu_lin'])))
