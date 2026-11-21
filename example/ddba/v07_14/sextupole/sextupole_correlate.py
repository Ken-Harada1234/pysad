# -*- coding: utf-8 -*-
"""解析量 (rdt_table.json) と追跡の判定 (results/track_<名前>.json) の順位相関 (Spearman の ρ と p 値、R7)。
使い方: py -3.12 sextupole_correlate.py   (sextupole_rdt_table.py と sextupole_tracking_check.py の後)
群: all (追跡できた全部)、near (目的関数が最終の 30 倍未満の組 = 最終の近く。単項の最適化 3 つを除く)、same_xi (ξ が最終と同じ組だけ = rev9・v07_11 を除く)、random (最終の周りの乱数 16 個)、historic (歴史的な 5 つのうち追跡できたもの)、
optimized (最適化の結果 8 個)。追跡の指標は「大きいほど良い」(DA 面積・x-δ の幅・LMA) と「小さいほど良い」(周波数マップのカオス割合・D の平均) が混ざるので、
符号をそろえて「良いほど大きい」にして相関を出す (解析量は「小さいほど良い」。ρ < 0 が「解析量が小さい組ほど追跡が良い」)。
出力: results/correlation.md、results/correlation.out.txt、fig_corr.png、fig_track_summary.png
"""
import glob
import io
import re
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']

LOG = io.StringIO()


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')


rdt = json.load(open(os.path.join(C.OUT, 'rdt_table.json'), encoding='utf-8'))
tr = {}
for f in glob.glob(os.path.join(C.OUT, 'track_*.json')):
    if f.endswith('_quick.json') or re.search(r'_(err\d+|n\d+)(_|\.json)', os.path.basename(f)):
        continue
    d = json.load(open(f, encoding='utf-8'))
    tr[d['name']] = d

TRACK = {   # 指標: (取り出し関数, 「良いほど大きい」への符号)
    'DA面積': (lambda d: d['da']['area_mm2'], +1),
    'x_DA': (lambda d: d['da']['x_da_mm'], +1),
    'y_DA': (lambda d: d['da']['y_da_mm'], +1),
    'xδ幅': (lambda d: d['xdp']['mean_width_mm'], +1),
    'FMカオス割合': (lambda d: d['fmap']['frac_chaotic'], -1),
    'FM失われた割合': (lambda d: d['fmap']['frac_lost'], -1),
    'LMA': (lambda d: d['lma']['mean_min'], +1),
}
ANA = {     # 解析量 (小さいほど良い)
    '目的関数(全項)': lambda r: r['obj'],
    'W3 (1次RDT)': lambda r: r['parts']['w3'],
    'W4 (4次非対角)': lambda r: r['parts']['w4'],
    'ADTS項': lambda r: r['parts']['adts'],
    'δの1次(chrom)': lambda r: r['parts']['chrom'],
    'Σ|h3| (5項)': lambda r: r['h3_sum'],
    '|ADTS|': lambda r: r['adts_norm'],
    '|ADTS_xx|': lambda r: abs(r['adts'][0]),
    '|ADTS_xy|': lambda r: abs(r['adts'][1]),
    '|ADTS_yy|': lambda r: abs(r['adts'][2]),
    '4次非対角Σ|h|': lambda r: r['h4_nondiag_sum'],
    '|ξ2| (CALC)': lambda r: float(np.hypot(*r['xi2_calc'])),
    '|∂ADTS/∂δ|': lambda r: r['dadts_norm'],
}
GROUPS = {
    'all': lambda n: True,
    'same_xi': lambda n: n not in ('rev9', 'v07_11'),
    'near': lambda n: rdt[n]['obj_rel'] < 30 and n not in ('opt_rdt1', 'opt_adts', 'opt_rdt4'),
    'random': lambda n: n.startswith('rnd'),
    'historic': lambda n: rdt[n]['group'] == 'historic',
    'optimized': lambda n: rdt[n]['group'] == 'optimized',
}
names_all = [n for n in tr if n in rdt]
P('追跡できた組 %d 個: %s' % (len(names_all), ' '.join(names_all)))
if any(n not in tr for n in rdt):
    P('追跡がない組:', ' '.join(n for n in rdt if n not in tr))

# 追跡の表
P('\n追跡の判定の表')
P('| 名前 | DA面積 mm² | x_DA mm | y_DA mm | x-δ幅 mm | FM カオス割合 | FM 失われた | meanD | LMA 平均 % | LMA δ+ % | LMA |δ−| % |')
P('|' + '---|' * 11)
for n in names_all:
    d = tr[n]
    g = lambda k, f: (f % k) if k is not None else '-'
    P('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        n, '%.1f' % d['da']['area_mm2'] if 'da' in d else '-', '%.2f' % d['da']['x_da_mm'] if 'da' in d else '-', '%.2f' % d['da']['y_da_mm'] if 'da' in d else '-',
        '%.1f' % d['xdp']['mean_width_mm'] if 'xdp' in d else '-', '%.2f' % d['fmap']['frac_chaotic'] if 'fmap' in d else '-',
        '%.2f' % d['fmap']['frac_lost'] if 'fmap' in d else '-', '%.2f' % d['fmap']['meanD'] if 'fmap' in d else '-',
        '%.2f' % (100 * d['lma']['mean_min']) if 'lma' in d else '-', '%.2f' % (100 * d['lma']['mean_plus']) if 'lma' in d else '-',
        '%.2f' % (100 * d['lma']['mean_minus']) if 'lma' in d else '-'))

res = {}
for gname, sel in GROUPS.items():
    ns = [n for n in names_all if sel(n)]
    P('\n[群 %s: %d 組] Spearman ρ (行 = 解析量 (小さいほど良い)、列 = 追跡 (符号をそろえて「良いほど大きい」)。ρ < 0 なら解析量が小さい組ほど追跡が良い。括弧内は p 値)' % (gname, len(ns)))
    if len(ns) < 5:
        P('  組が少ない (<5) ので出さない')
        continue
    P('| 解析量 | ' + ' | '.join(TRACK) + ' |')
    P('|' + '---|' * (len(TRACK) + 1))
    for an, fa in ANA.items():
        cells = []
        for tn, (ft, sg) in TRACK.items():
            xs, ys = [], []
            for n in ns:
                try:
                    y = sg * ft(tr[n])
                    x = fa(rdt[n])
                except (KeyError, TypeError):
                    continue
                if y is None or not np.isfinite(y):
                    continue
                xs.append(x)
                ys.append(y)
            if len(xs) < 5 or len(set(ys)) < 2:
                cells.append('-')
                continue
            rho, p = spearmanr(xs, ys)
            res[(gname, an, tn)] = (float(rho), float(p), len(xs))
            cells.append('%+.2f (%.2g)' % (rho, p))
        P('| %s | ' % an + ' | '.join(cells) + ' |')

# 図: 目的関数 / ADTS 項 / W3 と DA 面積・LMA の散布
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    col = {'historic': 'C0', 'optimized': 'C1', 'random': 'C2'}
    fig, axs = plt.subplots(2, 4, figsize=(17, 8))
    pairs = [('目的関数(全項)', 'DA面積'), ('ADTS項', 'DA面積'), ('W3 (1次RDT)', 'DA面積'), ('W4 (4次非対角)', 'DA面積'),
             ('目的関数(全項)', 'LMA'), ('ADTS項', 'FMカオス割合'), ('W3 (1次RDT)', 'FMカオス割合'), ('目的関数(全項)', 'xδ幅')]
    for ax, (an, tn) in zip(axs.ravel(), pairs):
        for n in names_all:
            try:
                x, y = ANA[an](rdt[n]), TRACK[tn][0](tr[n])
            except (KeyError, TypeError):
                continue
            ax.scatter(x, y, c=col[rdt[n]['group']], s=25)
            if n in ('final', 'v07_12', 'retune1', 'opt_full_1', 'opt_adts', 'opt_rdt1', 'opt_rdt4'):
                ax.annotate(n, (x, y), fontsize=6)
        ax.set_xscale('log')
        ax.set_xlabel(an + ' (小さいほど良い)', fontsize=8)
        ax.set_ylabel(tn, fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle('解析量と追跡の判定 (青: 歴史的な組、橙: 最適化、緑: 最終の周りの乱数)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_corr.png'), dpi=120)
    # DA の形 (光線) の比較
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for n in ('v07_12', 'retune1', 'final', 'opt_full_1', 'opt_adts', 'opt_rdt1', 'rnd10_0', 'rnd30_5'):
        if n in tr and 'da' in tr[n]:
            d = tr[n]['da']
            xs = [0.0] + d['x_mm']
            ys = [d['y_mm'][0] * 0] + d['y_mm']
            ax[0].plot(d['x_mm'], d['y_mm'], 'o-', ms=3, lw=1, label=n)
    ax[0].set_xlabel('x [mm] (リング先頭、βx=5.4 m)')
    ax[0].set_ylabel('y [mm] (βy=5.1 m)')
    ax[0].set_title('x-y の DA (512 周、δ=0、RF・放射なし、LOSSAMPL ±50 mm)', fontsize=9)
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)
    for n in ('v07_12', 'retune1', 'final', 'opt_full_1', 'opt_adts', 'opt_rdt1'):
        if n in tr and 'xdp' in tr[n]:
            d = tr[n]['xdp']
            ax[1].plot(np.array(d['dp']) * 100, d['x_plus_mm'], '-o', ms=3, lw=1, label=n)
            ax[1].plot(np.array(d['dp']) * 100, -np.array(d['x_minus_mm']), '-o', ms=3, lw=1, color=ax[1].lines[-1].get_color())
    ax[1].set_xlabel('δ [%]')
    ax[1].set_ylabel('x の DA (±) [mm]')
    ax[1].set_title('x-δ の DA (y=0.1 mm、512 周)', fontsize=9)
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_track_summary.png'), dpi=120)
except Exception as ex:
    P('図の作成に失敗:', ex)

# DA の比較表と図: 最終 (final)・rev9・最適化の結果。512 周 (誤差なし) / 1024 周 / 六極の位置ずれ 50 μm (3 種の乱数)
def load_json(f):
    return json.load(open(f, encoding='utf-8'))


cmp_names = [n for n in ('v07_12', 'retune1', 'final', 'rev9', 'opt_full_1', 'opt_full_2', 'opt_full_k125', 'opt_J3e-07', 'opt_J3e-06', 'opt_J1e-05', 'opt_dadts', 'opt_adts') if n in tr]
P('\nDA の比較 (x-y の DA 面積 mm²、リング先頭、δ=0。誤差なし 512 周 / 1024 周 / 六極の位置ずれ rms 50 μm の 3 種 (補正なし) の平均と範囲) と x-δ の幅 (δ = −3〜+3% の平均、mm)、LMA (4 点の min(δ+,|δ−|) の平均 %)')
P('| 名前 | DA 512 周 | DA 1024 周 | 位置ずれ DA (3 種) 平均 [最小, 最大] | x-δ 幅 512 周 | x-δ 幅 位置ずれ 平均 | LMA % | 最終に対する DA の比 (512 周) |')
P('|---|---|---|---|---|---|---|---|')
da_rows = {}
for n in cmp_names:
    d0 = tr[n]
    f1 = os.path.join(C.OUT, 'track_%s_n1024.json' % n)
    d1 = load_json(f1) if os.path.exists(f1) else None
    es = [load_json(os.path.join(C.OUT, 'track_%s_err%d.json' % (n, sd))) for sd in (1, 2, 3) if os.path.exists(os.path.join(C.OUT, 'track_%s_err%d.json' % (n, sd)))]
    ea = [e['da']['area_mm2'] for e in es]
    ex = [e['xdp']['mean_width_mm'] for e in es]
    da_rows[n] = dict(da=d0['da']['area_mm2'], da1024=d1['da']['area_mm2'] if d1 else None, err=ea, xdp=d0['xdp']['mean_width_mm'], xdp_err=ex)
    P('| %s | %.1f | %s | %s | %.1f | %s | %s | %.2f |' % (
        n, d0['da']['area_mm2'], '%.1f' % d1['da']['area_mm2'] if d1 else '-', '%.1f [%.1f, %.1f]' % (np.mean(ea), min(ea), max(ea)) if ea else '-',
        d0['xdp']['mean_width_mm'], '%.1f' % np.mean(ex) if ex else '-', '%.2f' % (100 * d0['lma']['mean_min']) if 'lma' in d0 else '-',
        d0['da']['area_mm2'] / tr['final']['da']['area_mm2']))
try:
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    xs_ = np.arange(len(cmp_names))
    ax[0].bar(xs_ - 0.2, [da_rows[n]['da'] for n in cmp_names], 0.4, label='誤差なし 512 周')
    for i, n in enumerate(cmp_names):
        e = da_rows[n]['err']
        if e:
            ax[0].plot([i + 0.2] * len(e), e, 'o', color='C1', ms=4, label='六極の位置ずれ 50 μm (3 種)' if n == 'retune1' else None)
        if da_rows[n]['da1024']:
            ax[0].plot(i - 0.2, da_rows[n]['da1024'], '_', color='k', ms=14, mew=2, label='誤差なし 1024 周' if n == 'retune1' else None)
    ax[0].set_xticks(xs_)
    ax[0].set_xticklabels(cmp_names, rotation=60, fontsize=8)
    ax[0].set_ylabel('x-y の DA 面積 [mm²] (リング先頭)')
    ax[0].legend(fontsize=8)
    ax[0].set_title('DA: 最終の組・rev9・pysad の解析式の最適化の結果', fontsize=9)
    ax[1].bar(xs_, [100 * tr[n]['lma']['mean_min'] if 'lma' in tr[n] else 0 for n in cmp_names])
    ax[1].set_xticks(xs_)
    ax[1].set_xticklabels(cmp_names, rotation=60, fontsize=8)
    ax[1].set_ylabel('LMA の平均 [%] (4 点、RF 6.5 MV・放射あり、±20 mm)')
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_da_compare.png'), dpi=120)
except Exception as ex_:
    P('図に失敗:', ex_)

open(os.path.join(C.OUT, 'correlation.out.txt'), 'w', encoding='utf-8').write(LOG.getvalue())
open(os.path.join(C.OUT, 'correlation.md'), 'w', encoding='utf-8').write(LOG.getvalue())
