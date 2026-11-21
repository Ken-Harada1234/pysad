# -*- coding: utf-8 -*-
"""example/ddba/v07_14 の .sad を作る。 ddba_v07_14_src/<フォルダ>/NAME.body の前に、ラティスを貼って
example/ddba/v07_14/<フォルダ>/NAME.sad にする。
  - 先頭の行が `!@LAT 2.5GeV` なら lattice/v07_14_2.5GeV.sad、無ければ lattice/v07_14_5GeV.sad を貼る (その行は残す)。
  - 本体に `!@PARAM name = 値` の行があり、環境変数 V0714_FULL=1 のときは、`!@FULL name = 値` の行で置き換える
    (完全な設定の版は NAME_full.sad として別に作る: `--full`)。
使い方: py -3.12 example/tools/mk_ddba_v07_14.py [NAME ...]   (名前なしは全部)
         py -3.12 example/tools/mk_ddba_v07_14.py --full [NAME ...]  (完全な設定の版 NAME_full.sad)"""
import os, sys, glob, re
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.join(os.path.dirname(HERE), 'ddba', 'v07_14')
SRC = os.path.join(HERE, 'ddba_v07_14_src')


def build(only=(), full=False):
    out = []
    for b in sorted(glob.glob(os.path.join(SRC, '*', '*.body'))):
        d = os.path.basename(os.path.dirname(b)); n = os.path.basename(b)[:-5]
        if only and n not in only: continue
        body = open(b, encoding='utf-8').read().strip('\n')
        m = re.match(r'!@LAT (\S+)', body)
        lat = m.group(1) if m else '5GeV'
        has_full = '!@FULL' in body
        if full and not has_full: continue
        lines = []
        full_over = {}
        for l in body.split('\n'):
            mm = re.match(r'!@FULL\s+(\w+)\s*=\s*(.*)$', l)
            if mm: full_over[mm.group(1)] = mm.group(2)
        for l in body.split('\n'):
            if l.startswith('!@FULL'): continue
            mm = re.match(r'!@PARAM\s+(\w+)\s*=\s*(.*)$', l)
            if mm:
                name, val = mm.group(1), mm.group(2)
                if full and name in full_over:
                    lines.append('%s = %s;    ! 完全な設定' % (name, full_over[name]))
                else:
                    lines.append('%s = %s;' % (name, val))
            else:
                lines.append(l)
        body2 = '\n'.join(lines)
        latt = open(os.path.join(EX, 'lattice', 'v07_14_%s.sad' % lat), encoding='utf-8').read().rstrip('\n')
        p = os.path.join(EX, d, n + ('_full' if full else '') + '.sad'); os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, 'w', encoding='utf-8', newline='\n').write(latt + '\n' + body2 + '\n')
        out.append(p)
    return out


if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    for p in build(a, '--full' in sys.argv): print(p)
