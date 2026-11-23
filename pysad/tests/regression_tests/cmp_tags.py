#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""cmp_tags.py — SAD の出力 (.result) と pysad の出力を、行頭のタグ (T001 など 英大文字 1 字 + 数字) 単位で比較する。

使い方: python cmp_tags.py [--num ABS REL] SAD.result PYSAD.out [許容する差のタグ ...]
  - 既定: 同じタグの行群が文字列として完全一致すれば OK ($FORM の桁まで含めて比較。末尾の空白は無視)。
  - --num ABS REL: 行を「数値とそれ以外」に分け、非数値部分が一致し、各数値が |a-b| <= ABS + REL*|b| なら OK
    (R1 などの丸め誤差 1e-16 級の値や DETR の雑音を許す。数値以外の部分の空白差は無視)。
  - ブロック: スクリプトが `Print["L001 begin"]; ...; Print["L001 end"];` で挟んだ範囲の出力行 (複数行) を 1 つのタグ L001 として比較する。
    SAD の ECHO (文の先頭が 1 空白 + 既知の文頭語の行) は除いて比較する。
  - 引数のタグ (SAD が未定義で pysad が拡張として評価するもの等) は「既知の拡張」として数えるだけで失敗にしない。
終了コード: 既知以外の差があれば 1。
"""
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
# SAD は読んだ文をそのまま (先頭に空白 1 つ付けて) エコーする
ECHO = re.compile(r'^ (Print|PrintBeamLine|WriteBeamLine|USE|VISIT|BYE|CELL|CALC|FFS|aaa|stop|Element\[|SetElement|LINE|DP0|DP|CHRO|CHROMATICITY)\b')


def load(path):
    return parse(open(path, encoding='utf-8', errors='replace'))


def parse(lines):
    """行の列 → {タグ: 行のリスト}"""
    d = {}
    block = None
    for l in lines:
        l = l.rstrip('\r\n').rstrip()
        m = re.match(r'([A-Z]\d+) begin$', l)
        if m:
            block = m.group(1)
            d[block] = []
            continue
        if block is not None:
            if l == block + ' end':
                block = None
            elif l and not ECHO.match(l):
                d[block].append(l)
            continue
        m = re.match(r'([A-Z]\d+)(?: |$)', l)
        if m:
            d.setdefault(m.group(1), []).append(l)
    return d


def same_num(a, b, tabs, trel):
    if a == b:
        return True
    pa = NUM.findall(a)
    pb = NUM.findall(b)
    if len(pa) != len(pb):
        return False
    ra = re.sub(r'\s+', '', NUM.sub('#', a))
    rb = re.sub(r'\s+', '', NUM.sub('#', b))
    if ra != rb:
        return False
    return all(abs(float(x) - float(y)) <= tabs + trel * abs(float(y)) for x, y in zip(pa, pb))


def compare(da, db, known=(), tabs=None, trel=None, verbose=True):
    """(不一致のタグのリスト, 既知の拡張の数, タグ総数)"""
    known = set(known)
    bad = []
    ext = 0
    for k in sorted(set(da) | set(db)):
        la, lb = da.get(k), db.get(k)
        ok = la == lb
        if not ok and tabs is not None and la is not None and lb is not None and len(la) == len(lb):
            ok = all(same_num(x, y, tabs, trel) for x, y in zip(la, lb))
        if not ok:
            if k in known:
                ext += 1
                continue
            bad.append(k)
            if verbose:
                print('DIFF', k)
                print('  SAD  :', str(la)[:600])
                print('  pysad:', str(lb)[:600])
    return bad, ext, len(set(da) | set(db))


def main():
    argv = sys.argv[1:]
    tabs = trel = None
    if argv and argv[0] == '--num':
        tabs, trel = float(argv[1]), float(argv[2])
        argv = argv[3:]
    da, db = load(argv[0]), load(argv[1])
    bad, ext, n = compare(da, db, argv[2:], tabs, trel)
    print(f'tags {n}  diff {len(bad)}  known-extension {ext}')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
