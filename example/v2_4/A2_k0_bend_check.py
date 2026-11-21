# -*- coding: utf-8 -*-
"""v2.4.0 A2: K0≠0 の偏向磁石 (挿入光源型) で E1・E2≠0 (ANGLE≠0 のとき SAD の E1・E2 は ANGLE の倍率)、ANGLE≠0 かつ K0≠0、ROTATE がある場合の、
resdrive の解析式 (3 次の共鳴駆動項 h・クロマティシティ・ADTS) を、独立な数値 (追跡の 1 周の写像のテイラー展開) と比べる。
リング: 小さな FODO (3 GeV、105.6 m、24 セル、SF・SD の K2L = 4, -6) の 1 セルに、K0 の逆符号の 2 台の偏向磁石 (ブロック) を挿入。
比較: h (h21000, h30000, h10110, h10020, h10200) = rdttrack.ring_map + rdt_from_map (追跡の 2 次のテイラー展開)、
ADTS = adts(fourth='r6') と nform.adts_from_map (3 次のテイラー展開 + 数値の正規形)、クロマティシティ = 自然 + 偏向磁石の分と CALC。
使い方: py -3.12 example/v2_4/A2_k0_bend_check.py [variant ...]    (variant: base, ang0, ang, edge, rot)"""
import math
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import compat, nform, rdt2nd as Q, rdttrack as T, resdrive as R
from pysad.parser import load_sad_string
compat.set_mode('k64')

FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9
NAMES = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')

# ブロック: 5 台の偏向磁石 (長さ 0.1 m)。K0 の列 (k/2, -k, k, -k, k/2) で閉軌道は局所のバンプ (ブロックの外に漏れない)。ANGLE の列も同じ符号の列 (a/2, -a, a, -a, a/2)。
# E1・E2 は SAD の規則 (ANGLE の倍率)、ROTATE [rad]、K1。
VARIANTS = {
    'ang0': dict(k=0.002, a=0.0, e=(0, 0), rot=0.0, k1=0.0, doc='ANGLE=0、K0 の列 (±0.002)、E なし (MPW 型)'),
    'ang': dict(k=0.002, a=0.004, e=(0, 0), rot=0.0, k1=0.0, doc='ANGLE の列 (±0.004) かつ K0 の列 (±0.002)、E なし'),
    'edge': dict(k=0.002, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.0, doc='ANGLE=±0.004 かつ K0=±0.002、E1=E2=ANGLE/2'),
    'edgeasym': dict(k=0.002, a=0.004, e=(1.0, 0.0), rot=0.0, k1=0.0, doc='ANGLE=±0.004 かつ K0=±0.002、E1=ANGLE (入口だけ)'),
    'k1': dict(k=0.002, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.1, doc='上 (E1=E2=ANGLE/2) + K1=0.1'),
    'k1_k00': dict(k=0.0, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.1, doc='K0 = 0 (普通の偏向磁石) + K1=0.1、E1=E2=ANGLE/2 (K0 の効果の切り分け)'),
    'k1_noedge': dict(k=0.002, a=0.004, e=(0, 0), rot=0.0, k1=0.1, doc='K0 の列 + K1=0.1、E なし (端の効果の切り分け)'),
    'k1_big': dict(k=0.0, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.3, doc='K0 = 0 + K1=0.3 (K1 を大きく)'),
    'k1_k00_leak': dict(k=0.0, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.1, steer=4e-5, doc='K0 = 0 の K1 = 0.1 のブロック + 別の ANGLE=0 の K0 の磁石 (閉軌道を 1e-4 m 漏らす。K0 と K1 を同じ磁石に持たせず、閉軌道の漏れだけを再現)'),
    'k1_kq': dict(k=0.0005, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.1, doc='K0 を 1/4 (±0.0005)、K1=0.1、E1=E2=ANGLE/2'),
    'k1_k2': dict(k=0.004, a=0.004, e=(0.5, 0.5), rot=0.0, k1=0.1, doc='K0 を 2 倍 (±0.004)、K1=0.1、E1=E2=ANGLE/2'),
    'rot': dict(k=0.002, a=0.004, e=(0, 0), rot=0.02, k1=0.0, doc='ANGLE の列 (±0.004) かつ K0 の列 (±0.002)、ROTATE=0.02 rad'),
    'rot2': dict(k=0.002, a=0.004, e=(0, 0), rot=0.2, k1=0.0, doc='ANGLE の列 (±0.004) かつ K0 の列 (±0.002)、ROTATE=0.2 rad (大きな傾き)'),
    'bigk0': dict(k=0.01, a=0.02, e=(0.5, 0.5), rot=0.0, k1=0.0, doc='K0 の列 ±0.01、ANGLE の列 ±0.02、E1=E2=ANGLE/2 (大きな K0)'),
    'rot0': dict(k=0.002, a=0.0, e=(0, 0), rot=0.02, k1=0.0, doc='ANGLE=0、K0 の列 (±0.002)、ROTATE=0.02 rad'),
}


def build(var):
    v = VARIANTS[var]
    pat = v.get('pat', (0.5, -1.0, 1.0, -1.0, 0.5))
    defs = []
    for i, f in enumerate(pat):
        x = 'L=0.1 K0=%r' % (v['k'] * f)
        if v['a']:
            x += ' ANGLE=%r E1=%r E2=%r' % (v['a'] * f, v['e'][0], v['e'][1])
        if v['k1']:
            x += ' K1=%r' % v['k1']
        if v['rot']:
            x += ' ROTATE=%r RAD' % v['rot']
        defs.append('W%d = (%s)' % (i + 1, x))
    s = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
    if v.get('steer'):
        defs.append('WS = (L=0.05 K0=%r)' % v['steer'])
    nl = chr(10)
    s = s.replace('LINE CELL =', 'BEND ' + ' '.join(defs) + ';' + nl + 'DRIFT D1X = (L=0.25);' + nl
                  + 'LINE CELLX = (SF DS QF D2 B D2 QD DS SD D1X W1 W2 W3 W4 W5' + (' WS' if v.get('steer') else '') + ' D1X);' + nl + 'LINE CELL =')
    s = s.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF 23*CELL CELLX);')
    return load_sad_string(s).build_lattice('RING')


def P(*a):
    print(*a, flush=True)


def run(var):
    lat = build(var)
    rows = lat.twiss_table()
    P('## %s: %s' % (var, VARIANTS[var]['doc']))
    orb = max(abs(r.x_cod) for r in rows)
    P('  閉軌道の最大 |x| = %.3e m、νx νy = %.5f %.5f' % (orb, rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi))
    r = R.resonance_terms(lat, NAMES, rows=rows)
    rm = T.ring_map(lat.elements, P0, rows[0], amp=1e-5)
    h = T.rdt_from_map(rm)
    rel = {n: abs(r.total[n] - h[n]) / abs(h[n]) for n in NAMES}
    P('  h (解析 対 追跡の 2 次のテイラー展開) の相対差: ' + ' '.join('%s %.1e' % (n, rel[n]) for n in NAMES))
    P('  警告: %s' % [w for w in r.warnings if '偏向' in w or 'K0' in w or 'ROTATE' in w])
    ref = R.chromaticity_numeric(lat, order=1)['d1']
    nat = R.chromaticity_natural(lat, rows=rows)
    bend = R.chromaticity_bend(lat, rows=rows)
    P('  ξ (自然 + 偏向磁石の分) − CALC (6 極なしでないので 6 極の 1 次を足す): ' + str(
        [float(nat[i] + bend[i] + R.chromaticity_sextupole(lat, rows=rows)[i] - ref[i]) for i in (0, 1)]))
    A = R.adts(lat, rows=rows, fourth='r6')
    rm3 = T.ring_map3(lat.elements, P0, rows[0], amp3=1e-4)
    an, _ = nform.adts_from_map(rm3)
    P('  ADTS (fourth=r6) 対 数値の正規形 (追跡の 3 次のテイラー展開): ' + ' '.join(
        '%s %.4g/%.4g (%.1e)' % (k, A.total[k], an[k], abs(A.total[k] - an[k]) / abs(an[k])) for k in ('xx', 'xy', 'yy')))
    A3 = R.adts(lat, rows=rows)
    P('  ADTS (既定 r3) 対 数値の正規形: ' + ' '.join('%s %.1e' % (k, abs(A3.total[k] - an[k]) / abs(an[k])) for k in ('xx', 'xy', 'yy')))


if __name__ == '__main__':
    for v in (sys.argv[1:] or list(VARIANTS)):
        try:
            run(v)
        except Exception as e:
            import traceback
            P('## %s: 例外 %r' % (v, e))
            traceback.print_exc()
