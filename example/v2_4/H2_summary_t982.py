# -*- coding: utf-8 -*-
"""H2: v03_12 (RAD ON、982 周、9 位置の 8 つ = s = 0, 75.4, …, 527.8 m。最後の位置は最初と同じ) の H2_kekls12_rad_t982_p9_*.out.txt をまとめる。
使い方: py -3.12 H2_summary_t982.py > H2_kekls12_t982_summary.out.txt"""
import re, glob, numpy as np
def parse(f):
    t = open(f, encoding='utf-8', errors='replace').read()
    g = lambda k: np.array([float(x) for x in re.search(k + r' \{([^}]*)\}', t).group(1).replace(' ', '').split(',')]) * 100
    cod = re.findall(r'COD rms \[m\] x y:\s+(\S+)\s+(\S+)', t)
    return g('LMA DPPLUS')[:-1], g('LMA DPMINUS')[:-1], cod, float(re.search(r'LMA WALL \[s\]\s+(\S+)', t).group(1))
print('状態 | COD rms 補正前→後 (x, y) [μm] | δ+ [%] 平均 (最小) | δ− [%] 平均 (最悪 = 絶対値が最小) | LMA 所要 [s]')
rows = {}
for f in sorted(glob.glob('H2_kekls12_rad_t982_p9_*.out.txt')):
    tag = f.split('p9_')[1][:-8]
    p, m, cod, w = parse(f)
    rows[tag] = (p, m)
    c = '%.0f・%.0f → %.0f・%.0f' % (float(cod[0][0]) * 1e6, float(cod[0][1]) * 1e6, float(cod[-1][0]) * 1e6, float(cod[-1][1]) * 1e6)
    print('%s | %s | %.3f (%.3f) | %.3f (%.3f) | %.0f' % (tag, c, p.mean(), p.min(), m.mean(), m.max(), w))
cs = [t for t in rows if t.startswith('corr')]
P = np.concatenate([rows[t][0] for t in cs]); M = np.concatenate([rows[t][1] for t in cs])
print('補正あり %d 種 × %d 位置: δ+ 平均 %.3f 最小 %.3f、δ− 平均 %.3f 中央値 %.3f 最悪 %.3f。種ごとの δ− の平均 %s (標準偏差 %.2f)' % (
    len(cs), len(rows[cs[0]][0]), P.mean(), P.min(), M.mean(), np.median(M), M.max(), ' '.join('%.2f' % rows[t][1].mean() for t in cs), np.std([rows[t][1].mean() for t in cs])))
