# -*- coding: utf-8 -*-
"""example/kek_ls の .sad を作る。 kek_ls_src/<フォルダ>/NAME.body (FFS; 以降も含む台本の本体) の前に、
lattice/kek_ls_v03_68.sad を貼って example/kek_ls/<フォルダ>/NAME.sad にする。
使い方: py -3.12 example/tools/mk_kekls.py [NAME ...]   (名前なしは全部)"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.join(os.path.dirname(HERE), 'kek_ls')
SRC = os.path.join(HERE, 'kek_ls_src')
def build(only=()):
    lat = open(os.path.join(EX, 'lattice', 'kek_ls_v03_68.sad'), encoding='utf-8').read().rstrip('\n')
    out = []
    for b in sorted(glob.glob(os.path.join(SRC, '*', '*.body'))):
        d = os.path.basename(os.path.dirname(b)); n = os.path.basename(b)[:-5]
        if only and n not in only: continue
        body = open(b, encoding='utf-8').read().strip('\n')
        p = os.path.join(EX, d, n + '.sad'); os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, 'w', encoding='utf-8', newline='\n').write(lat + '\n' + body + '\n')
        out.append(p)
    return out
if __name__ == '__main__':
    for p in build(sys.argv[1:]): print(p)
