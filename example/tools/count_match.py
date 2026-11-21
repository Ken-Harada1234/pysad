"""example/*/README.md の表の最後の列の判定 (○ 一致・△ 差あり・× 未対応・失敗) を数える。
使い方: py -3.12 example/tools/count_match.py   (リポジトリ直下から)。README の一致度の表の数値の元。"""
import glob, os, re, sys, collections
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
tot = collections.Counter()
rows = []
for f in sorted(glob.glob(os.path.join(root, '*', 'README.md'))):
    d = os.path.basename(os.path.dirname(f))
    c = collections.Counter()
    for line in open(f, encoding='utf-8'):
        line = line.rstrip('\r\n')
        if not line.startswith('|') or set(line) <= set('|-: '):
            continue
        cells = [x.strip() for x in line.strip().strip('|').split('|')]
        if len(cells) < 2:
            continue
        last = cells[-1]
        m = re.match(r'^(○|△|×)', last)
        if m:
            c[m.group(1)] += 1
    if sum(c.values()):
        rows.append((d, c)); tot.update(c)
print('| フォルダ | ○ | △ | × | 計 |\n|---|---|---|---|---|')
for d, c in rows:
    print('| %s | %d | %d | %d | %d |' % (d, c['○'], c['△'], c['×'], sum(c.values())))
print('| **合計** | %d | %d | %d | %d |' % (tot['○'], tot['△'], tot['×'], sum(tot.values())))
