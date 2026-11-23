# -*- coding: utf-8 -*-
"""fmap_check_errors_scan.py の結果 (fmap_check_data/errors_<ring>.jsonl) を (誤差の種類, 強さ) ごとに集計する。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_errors_summary.py <astrid2|kekls>     標準出力 (fmap_check_errors_summary_<ring>.out.txt に保存して使う)
列: 種の数、D > −5 の割合 (残った点のうち) の平均 ± 標準偏差 [最小〜最大]、失われた点 (全 336 点のうち) の平均 ± 標準偏差、D の中央値、最小振幅の点の D (D0)、
閉軌道の rms (x, y) [mm] の平均、|c−| の見積もり (フィードダウン + ロール)、差の共鳴までの距離の見積もり (split_pred) と数値 (split_num)。
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ring = sys.argv[1]
rows = [json.loads(l) for l in open(os.path.join(HERE, 'fmap_check_data', 'errors_%s.jsonl' % ring), encoding='utf-8')]
g = defaultdict(list)
for d in rows:
    g[(d['kind'], d['k'])].append(d)
base = g.get(('none', 0.0))
print('%s: 21 x 16 = 336 点、1024 周、%d 組' % (ring, len(rows)))
print('%-22s %4s | %-26s | %-18s | %6s | %6s | %-11s | %7s | %7s %7s %7s' % (
    '種類 k', 'n', 'D>-5 の割合 平均±σ [最小,最大]', '失われた点 平均±σ (D>-5 か失われた割合)', 'D中央', 'D0', '軌道rms mm', '|c-|予想', 'split予', 'split数', 'ΔQ予/数'))
order = ['none', 'all', 'all_cs', 'quad_dxy', 'sext_dxy', 'oct_dxy', 'bend_dxy', 'bend_rot', 'quad_rot']
for kind in order:
    for k in sorted(set(kk for (kn, kk) in g if kn == kind)):
        L = g[(kind, k)]
        fh = np.array([d['stats']['frac_hi'] for d in L])
        nl = np.array([d['stats']['nlost'] for d in L])
        fa = np.array([d['stats']['frac_hi_all'] for d in L])
        dm = np.array([d['stats']['Dmed'] for d in L])
        d0 = np.array([d['D0'] if d['D0'] is not None else np.nan for d in L])
        orb = np.mean([d['errors']['orbit_rms'] for d in L], axis=0) * 1e3
        cm = np.mean([d['errors']['cminus_pred'] for d in L])
        sp = np.mean([d['errors']['split_pred'] for d in L])
        sn = np.nanmean([d['errors']['split_num'] for d in L])
        dq = np.mean([np.hypot(*np.array(d['errors']['dnu_lin'])) for d in L])
        dp = np.mean([np.hypot(*np.array(d['errors']['dnu_pred'])) for d in L])
        print('%-14s ×%-5g %4d | %.3f ± %.3f [%.3f, %.3f]   | %6.1f ± %5.1f (%.2f) | %6.2f | %6.2f | (%.2f, %.2f) | %7.4f | %7.4f %7.4f  %.4f/%.4f' % (
            kind, k, len(L), fh.mean(), fh.std(), fh.min(), fh.max(), nl.mean(), nl.std(), fa.mean(), np.mean(dm), np.nanmean(d0), orb[0], orb[1], cm, sp, sn, dp, dq))
# 相関: D > -5 の割合 と 各量 (all の全部)
A = [d for d in rows if d['kind'] == 'all' and d['stats']['frac_hi'] == d['stats']['frac_hi'] and d['errors'].get('stable', True)]
if len(A) > 4:
    fh = np.array([d['stats']['frac_hi'] for d in A])
    for lab, v in (('軌道 rms x', [d['errors']['orbit_rms'][0] for d in A]), ('軌道 rms y', [d['errors']['orbit_rms'][1] for d in A]),
                   ('|c-| 予想', [d['errors']['cminus_pred'] for d in A]), ('|ΔQ| 数値', [np.hypot(*d['errors']['dnu_lin']) for d in A]),
                   ('6 極の中心の軌道 rms y', [d['errors']['rms_orbit_sext'][1] for d in A]), ('失われた点', [d['stats']['nlost'] for d in A])):
        print('相関 (all の %d 組): D>-5 の割合 と %-18s r = %+.2f' % (len(A), lab, np.corrcoef(fh, v)[0, 1]))
