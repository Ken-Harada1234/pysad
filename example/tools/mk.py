# -*- coding: utf-8 -*-
"""example の .sad を作る: mk.py <フォルダ/名前> <本体ファイル> [ラティス名]  (本体は FFS; 以降)。ラティスを先頭に貼って .sad にする。"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(HERE)
def mk(rel, body, lat='FODO_ring'):
    L = open(os.path.join(EX, '00_lattice', lat + '.sad'), encoding='utf-8').read()
    p = os.path.join(EX, rel + '.sad'); os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w', encoding='utf-8', newline='\n').write(L.rstrip('\n') + '\nFFS;\n' + body.strip('\n') + '\nExit;\n')
    return p
if __name__ == '__main__':
    mk(sys.argv[1], open(sys.argv[2], encoding='utf-8').read(), *(sys.argv[3:4]))
