# -*- coding: utf-8 -*-
"""gs と pysad の出力の数値を比べる (各 .sad の Print の出力行だけ)。"""
import re, sys, os, glob
num = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?')
def out_lines(path, sad):
    txt = open(path, encoding='utf-8', errors='replace').read().split('\n')
    return [l for l in txt if l.strip() and not l.startswith(' ') and 'pysad' not in l and not l.startswith('Warning') and 'RuntimeWarning' not in l]
def cmp(base):
    a = out_lines(base + '.gs.txt', base); b = out_lines(base + '.pysad.txt', base)
    # 出力部だけ: Print の出力は最後の part。両方で共通の最後の n 行を使う
    ga = [l for l in a if num.search(l)]; pb = [l for l in b if num.search(l)]
    n = min(len(ga), len(pb)); ga = ga[-n:]; pb = pb[-n:]
    mx = 0; worst = None; cnt = 0
    for x, y in zip(ga, pb):
        nx = [float(t) for t in num.findall(x)]; ny = [float(t) for t in num.findall(y)]
        if len(nx) != len(ny): cnt += 1; continue
        for u, v in zip(nx, ny):
            d = abs(u - v) / max(abs(u), abs(v), 1e-6)
            if d > mx: mx = d; worst = (x[:80], y[:80])
    return mx, cnt, worst
for p in sys.argv[1:]:
    mx, cnt, w = cmp(p)
    print('%-60s maxrel %.2e  lines_with_different_count %d' % (os.path.basename(p), mx, cnt))
    if mx > 1e-6 and w: print('    ', w)
