# -*- coding: utf-8 -*-
"""元のラティス定義 (読み取り専用) から、6 極 20 ファミリーの K2L を読み出して sextupole_sets.json に書く (一度だけ走らせる)。
使い方: 環境変数 V07_14_SRC_DIR に元のラティス定義のあるフォルダを指定して  py -3.12 sxt_extract_sets.py
読み出す版 (ファイル名の一部は個人名を含むので、パターンで探し、名前は JSON に書かない):
  v07_11 (初期の 3 個 (S1=S2=0、S3=S5))、v07_12、6 極の再調整 第 1 版 (v07_12 系の `*_SX`、v07_14_2.5GeV と同じ組)、v07_14 の 20231205 再調整版 (最終)、v07_16 の sx_rev9 (さらに後の再調整)。
最終の組は example/ddba/v07_14/lattice/ のラティスの値と一致することを確かめる。
"""
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import sxtopt as S
from pysad.parser import load_sad_file

D = os.environ.get('V07_14_SRC_DIR')
if not D or not os.path.isdir(D):
    sys.exit('V07_14_SRC_DIR に元のラティス定義のフォルダを指定する')
files = sorted(os.listdir(D))


def find(pat):
    m = [f for f in files if re.fullmatch(pat, f)]
    if len(m) != 1:
        sys.exit('パターン %s の一致が %d 個' % (pat, len(m)))
    return os.path.join(D, m[0])


SPEC = [
    ('v07_11', r'v07_11', 'v07_11 の6 極 (初期の組。S1=S2=0、S3=S5 の 3 個)'),
    ('v07_12', r'v07_12', 'v07_12 の6 極 (v07_11 とほぼ同じ 3 個。M-A セルだけ少し違う)'),
    ('retune1', r'v07_12_.*_SX', '6 極の再調整 第 1 版 (v07_12 系のファイル。v07_14_2.5GeV と同じ組。2023-10)'),
    ('final', r'v07_14_[A-Za-z_]*20231205', '6 極の再調整版 (20231205。v07_14 の最終の組)'),
    ('rev9', r'v07_16_sx_rev9_5GeV', 'v07_16 の6 極の再調整 (sx_rev9。2024-05。v07_14 より後)'),
]
src0 = C.load_source()
_, lat0 = C.load_lattice('RING', src=src0)
names, fam, kref = S.sextupole_families(lat0)
out = {'family_names': names, 'note': 'K2L [1/m^2] (SAD の SEXT の K2 = 積分強さ、長さ 0.1 m)。ラティスの配置 (RING) は含まない', 'sets': {}}
for key, pat, label in SPEC:
    p = load_sad_file(find(pat))
    lat = p.build_lattice('RING')
    k = []
    for n in names:
        es = [e for e in lat.elements if type(e).__name__ == 'Sextupole' and e.name == n]
        ks = {round(float(e.K2) * float(e.L), 9) for e in es}
        assert len(ks) == 1, (key, n, ks)
        k.append(ks.pop())
    out['sets'][key] = dict(k=k, label=label)
    print(key, np.round(k, 3).tolist())
kf = np.array(out['sets']['final']['k'])
print('最終の組とラティスの K2L の最大差:', np.max(np.abs(kf - kref)))
assert np.max(np.abs(kf - kref)) < 1e-6
json.dump(out, open(C.SETS_JSON, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
