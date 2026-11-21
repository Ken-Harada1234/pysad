# -*- coding: utf-8 -*-
"""入門 §5: ADTS の内訳 (大きなリング A = kek-ls_v03_12.lat、3 GeV、603.12 m、6 極の強さ ×1) の棒グラフと表。
数値は保存済みの出力 R3_scan_kekls.pysad.txt の「6 極の強さ ×1」の節から読む (走らせ直さない)。単位 1/(m rad)。
使い方 (リポジトリ直下から): py -3.12 example/v2_3/primer_breakdown.py   出力: primer_breakdown.png、標準出力
"""
import os
import re
import numpy as np
import primer_common as C
plt = C.setup_plot()

txt = open(os.path.join(C.HERE, 'R3_scan_kekls.pysad.txt'), encoding='utf-8').read()
blk = txt.split('--- 6 極の強さ ×1 (')[1].split('--- 6 極の強さ ×2')[0]
rows = {}
for ln in blk.splitlines():
    m = re.match(r'\s+(\S.*?)\s+xx=\s*(\S+)\s+xy=\s*(\S+)\s+yy=\s*(\S+)', ln)
    if m:
        rows[m.group(1).strip()] = [float(m.group(i)) for i in (2, 3, 4)]
rows = {k: v for k, v in rows.items() if not k.startswith('解析の合計 (R2') and not k.startswith('追跡')}
print(txt.splitlines()[0])
print('6 極の強さ ×1 の内訳 [1/(m rad)]:')
for k, v in rows.items():
    print('  %-34s xx = %12.3f  xy = %12.3f  yy = %12.3f' % (k, *v))
names = ['6 極 2 次', '8 極 1 次', '運動学的 p^4/8', 'K1 の端のフリンジ', '偏向磁石の 3 次どうし (R3)', '6 極 × 偏向磁石の交差項 (R3)']
vals = np.array([rows[n] for n in names])
tot = np.array(rows['解析の合計 (R3、偏向磁石の 3 次あり)'])
print('和 (6 項) = %s、解析の合計 (R3) = %s' % (np.array2string(vals.sum(axis=0), precision=3), np.array2string(tot, precision=3)))
fig, a = plt.subplots(figsize=(9, 4.8))
w = 0.13
for i, n in enumerate(names):
    a.bar(np.arange(3) + (i - 2.5) * w, vals[i], w, label=n)
a.bar(np.arange(3) + 3.5 * w, tot, w, color='k', label='合計 (Total)')
a.set_yscale('symlog', linthresh=1)
a.set_xticks(range(3)); a.set_xticklabels(['xx = ∂νx/∂Jx', 'xy = ∂νx/∂Jy', 'yy = ∂νy/∂Jy'])
a.axhline(0, color='gray', lw=0.5)
a.set_ylabel('1/(m rad)  (symlog)')
a.set_title('大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m) の ADTS の内訳 (6 極 ×1)', fontsize=10)
a.set_ylim(-3e5, 5e7)
a.legend(fontsize=7.5, ncol=2, loc='upper right'); a.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(C.HERE, 'primer_breakdown.png'), dpi=110)
