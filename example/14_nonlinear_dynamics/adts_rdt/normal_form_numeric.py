# -*- coding: utf-8 -*-
"""数値の正規形 (pysad/nform.py): 追跡の 1 周の写像の 3 次のテイラー展開から ∂ν/∂J を取り出し、解析式 (rdt4 + 正規形の補正)・直接の追跡 (NAFF、rdt4_adts_*.pysad.json) と比べる。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/adts_rdt/normal_form_numeric.py fodo astrid astrid_mpw kekls   (rdt4_tracking.py adts <ring> を先に走らせておく)
条件: pysad k64 の追跡、RF・放射なし、閉軌道のまわり、δ=0。ring_map3 の amp3 = FODO・ASTRID2 1e-4、大きなリング A 1e-5 √m (正規化座標)。出力 normal_form_numeric.out.txt"""
import json
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import nform, rdttrack as T, rdt2nd as Q, compat
compat.set_mode('k64')
import adts_bend as M

A3 = {'fodo': 1e-4, 'astrid': 1e-4, 'astrid_mpw': 1e-4, 'kekls': 1e-5}
out = ['R5 数値の正規形 (J の 1 次まで): [1/(m rad)]。nform = 追跡の写像の h3・h4 + 正規形の補正、解析 = rdt4 + 補正、追跡 = 直接 (追跡 + NAFF + フィット)']
for key in (sys.argv[1:] or list(A3)):
    S = M.setup(key)
    rm = T.ring_map3(S['lat'].elements, S['p'], S['rows'][0], amp3=A3[key])
    a, info = nform.adts_from_map(rm)
    a4, _ = Q.adts_from_rdt4(S['lat'], rows=S['rows'])
    d = json.load(open(os.path.join(HERE, 'rdt4_adts_%s.pysad.json' % key), encoding='utf-8'))['track']
    out.append('%s (%.2f GeV、周長 %.1f m、ν = %.6f, %.6f; f4 のシンプレクティック残差 %.1e)' % (key, S['p'] / 1e9, S['circ'], info['tunes'][0], info['tunes'][1], info['symp']))
    for k in ('xx', 'xy', 'yy'):
        out.append('   %s: nform %+.4f  解析 %+.4f  追跡 %+.4f   (nform−追跡)/追跡 %+.1e   (解析−追跡)/追跡 %+.1e' % (k, a[k], a4[k], d[k], (a[k] - d[k]) / d[k], (a4[k] - d[k]) / d[k]))
print('\n'.join(out))
open(os.path.join(HERE, 'normal_form_numeric.out.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
