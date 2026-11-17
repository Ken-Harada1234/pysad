# -*- coding: utf-8 -*-
"""example_PF_ASTRID の .sad を作る。mkpf(rel, body, lat) は lattice/<lat>.sad を先頭に貼り、本体 (FFS; 以降) と Exit; をつなぐ。
走らせる: py -3.12 tools/runpf.py T2_samples/NAME.sad [timeout秒]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(HERE)
def mk(rel, body, lat='v9_16_asym', pre=''):
    L = open(os.path.join(EX, 'lattice', lat + '.sad'), encoding='utf-8').read()
    p = os.path.join(EX, rel + '.sad'); os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w', encoding='utf-8', newline='\n').write(L.rstrip('\n') + '\n' + pre.strip('\n') + '\nFFS;\n' + body.strip('\n') + '\nExit;\n')
    return p
