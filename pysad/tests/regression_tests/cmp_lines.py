# -*- coding: utf-8 -*-
"""cmp_lines.py SAD.result PYSAD.out — Print[] のラベル付き行 ("LABEL: 数 数 ..." / "LABEL 数 数") を順に突き合わせ、
行ごとの最大 相対差 (max(|SAD|,1) で割る) を出す。ラベルは行頭の英数字+記号の語 (同じ語の行は出現順に対応づける)。
ECHO 行 (先頭 1 空白で始まる SAD のエコー) は除く。run_regression.py からは rows()/compare() を使う。"""
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
LAB = re.compile(r'^([A-Za-z][A-Za-z0-9_.\[\]]*):?\s')


def rows(lines):
    out = []
    for l in lines:
        l = l.rstrip('\n')
        if l.startswith(' ') or not l.strip():
            continue
        m = LAB.match(l)
        if not m:
            continue
        nums = [float(x) for x in NUM.findall(l[m.end():])]
        if nums:
            out.append((m.group(1), nums, l))
    return out


def compare(a, b):
    """a (SAD), b (pysad) の rows → [(ラベル, 最大相対差 or None (個数違い))]"""
    dic = {}
    for lab, nums, l in b:
        dic.setdefault(lab, []).append((nums, l))
    cnt = {}
    res = []
    for lab, nums, l in a:
        k = cnt.get(lab, 0)
        cnt[lab] = k + 1
        if lab not in dic or k >= len(dic[lab]):
            continue
        pn, pl = dic[lab][k]
        if len(pn) != len(nums):
            res.append((lab, None))
            continue
        res.append((lab, max(abs(x - y) / max(abs(x), 1.0) for x, y in zip(nums, pn))))
    return res


if __name__ == '__main__':
    a = rows(open(sys.argv[1], encoding='utf-8', errors='replace'))
    b = rows(open(sys.argv[2], encoding='utf-8', errors='replace'))
    worst = 0.0
    for lab, d in compare(a, b):
        if d is None:
            print(f'{lab}: 個数が違う')
            continue
        worst = max(worst, d)
        print(f'{lab:28s} {d:9.2e}')
    print('worst', worst)
