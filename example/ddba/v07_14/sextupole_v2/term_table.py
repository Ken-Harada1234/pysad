# -*- coding: utf-8 -*-
"""候補の解析式の値の表 (RDT・ADTS・ξ・誤差感度)。完成版との比較用 (README の表)。
  py -3.12 term_table.py NAME [NAME ...] [--md results/term_table.md]"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

S = V.S
ap = argparse.ArgumentParser()
ap.add_argument('names', nargs='+')
ap.add_argument('--md', default=None)
ap.add_argument('--err', type=float, default=0.03)
args = ap.parse_args()
src, lat, m, ks, xt = V.load()
cands = json.load(open(os.path.join(V.OUT, 'cands.json'), encoding='utf-8'))
terms = ['h21000', 'h30000', 'h10110', 'h10020', 'h10200', 'h20001', 'h00201', 'h10002', 'h20200', 'h20020', 'h31000', 'h11200', 'h00310', 'h40000', 'h20110']
so = S.SxtSpecObjective(m, ', '.join('%s=objective' % t for t in terms) + ', adts=objective, xi=ignore', warn=False)
o = S.SxtObjective(m)
head = ['候補'] + terms + ['ADTS xx', 'ADTS xy', 'ADTS yy', 'ξx', 'ξy', '誤差感度 h21000/h10110/h10020 (rms)', 'max abs K2L', 'w3', 'w4', 'adts']
rows = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
for n in args.names:
    k = ks['final'] if n == 'final' and n not in cands else np.array(cands[n]['k'])
    tb = {r['name']: r for r in so.table(k, xt)}
    er = m.err_rms(k, args.err)
    p = o.parts(k)
    cells = [n] + ['%.3g' % abs(tb[t]['value']) for t in terms] + ['%.0f' % tb['adts_xx']['value'], '%.0f' % tb['adts_xy']['value'], '%.0f' % tb['adts_yy']['value'],
             '%.2f' % m.xi(k)[0], '%.2f' % m.xi(k)[1], '%.1f / %.1f / %.1f' % (er['h21000'], er['h10110'], er['h10020']), '%.1f' % np.max(np.abs(k)),
             '%.2e' % p['w3'], '%.2e' % p['w4'], '%.2e' % p['adts']]
    rows.append('| ' + ' | '.join(cells) + ' |')
txt = '\n'.join(rows)
sys.stdout.buffer.write(txt.encode('utf-8') + bytes([10]))
if args.md:
    open(args.md, 'w', encoding='utf-8').write(txt + '\n')
