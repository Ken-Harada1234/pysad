# -*- coding: utf-8 -*-
"""保存済みの結果 (pysad_*.ref.json、elegant_*.touschek.txt、elegant_*.lma.txt) から、寿命と LMA の比較の表を作る (計算はしない)。

使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_summary.py   出力: 標準出力 (elegant_summary.out.txt に保存した)
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = {'astrid2': 'ASTRID2 (0.58 GeV、45.7 m、RF 50 kV)', 'kekls': 'kek-ls_v03_12.lat (3 GeV、603 m、RF 1.4 MV)'}


def load_ele(name):
    out = {}
    p = os.path.join(HERE, 'elegant_%s.touschek.txt' % name)
    if not os.path.exists(p):
        return out
    for ln in open(p, encoding='utf-8'):
        if ln.startswith('#') or ':' not in ln:
            continue
        lb, rest = ln.split(':', 1)
        m = re.search(r'tLifetime \(hour\) =\s*([-0-9.e+]+)', rest)
        if m:
            out[lb.strip()] = float(m.group(1)) * 3600.0
    return out


def main():
    for name in NAMES:
        ref = json.load(open(os.path.join(HERE, 'pysad_%s.ref.json' % name)))
        ele = load_ele(name)
        print('## %s' % NAMES[name])
        print('PBUNCH = %.3g 粒子/バンチ、IBS 込み: εx = %.4e、εy = %.4e (εy/εx = %.4f)、σδ = %.4e、σz = %.3f mm' % (
            ref['PBUNCH'], ref['EMITXC'], ref['EMITYC'], ref['EMITYC'] / ref['EMITXC'], ref['SIGE'], ref['SIGZ'] * 1e3))
        print('| 受容 δ | pysad [s] | elegant [s] | elegant / pysad − 1 |')
        print('|---|---|---|---|')
        for a in ref['uniform_a']:
            key = '%.10g' % a
            tp = ref['pysad_uniform'][key]
            k = [x for x in ele if x.startswith('uniform') and abs(float(x.split()[1]) - a) < 6e-5]
            te = ele[k[0]] if k else None
            tab = ' (表の点)' if any(abs(a - 0.002 * i) < 1e-9 for i in range(2, 122)) else ' (表の点の間: 対数補間)'
            if te:
                print('| 一様 %.4f%%%s | %.1f | %.1f | %+.1f%% |' % (a * 100, tab, tp, te, 100 * (te / tp - 1)))
        if 'pysad LMA' in ele:
            print('| 局所の LMA (pysad の δ±(s) と pysad の式 = TouschekLifetimeLMA) | %.1f | (elegant の式 + pysad の LMA: %.1f) | %+.1f%% |' % (
                ref.get('pysad_local', float('nan')), ele['pysad LMA'], 100 * (ele['pysad LMA'] / ref.get('pysad_local', float('nan')) - 1)))
        if 'elegant LMA' in ele:
            print('| 局所の LMA (pysad の式 ← pysad の LMA / elegant の式 ← elegant の LMA) | %.1f | %.1f | %+.1f%% |' % (
                ref.get('pysad_local', float('nan')), ele['elegant LMA'], 100 * (ele['elegant LMA'] / ref.get('pysad_local', float('nan')) - 1)))
        print()


if __name__ == '__main__':
    main()
