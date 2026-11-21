# -*- coding: utf-8 -*-
"""v2.3.0 S2 (2): 挿入光源 (ASTRID2_v01_02_id.sad の RINGMPW) で多くの点が失われる原因の切り分け。
MPW の 4νy = 9 の共鳴 (νy = 2.25) が原因かを、(A) MPW の強さ s = 0〜1 (偏向磁石の K0 を s 倍) で νy が動いたときの失われる点の割合、(B) s = 1 のまま四極 (QF3・QF2 の 2 家族)
の微調整で νy を ±0.01〜0.05 動かしたときの失われる点の割合、(C) s = 1・0.5 で νx・νy を MPW なしの作動点に戻した (MPW の線形のずれを四極で打ち消した) ときの割合、で調べる。
(C) は MPW の内部の非線形 (K0 ≠ 0 の偏向磁石) だけの寄与。
使い方: py -3.12 example/v2_3/S2_mpw.py [--jobs 3] [--only A|B|C]    出力: 標準出力、S2_data/S2_mpw.jsonl (ケースごとに 1 行、走った組は飛ばす)
条件: v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、位置 2、x = 0.5〜18 mm・y = 0.5〜14 mm の 21 × 16 (336 点)、1024 周、Hann 窓、規格化 = 標本の楕円。
"""
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, 'S2_data')
PATH = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
K0 = '0.0640435140143602'
Q3, Q2 = '5.939500', '6.070000'


def source(s=1.0, a=0.0, b=0.0, edge=1.0):
    src = open(PATH, encoding='utf-8').read()
    n0 = src.count(K0)
    src = src.replace(K0, repr(float(K0) * s))
    if edge != 1.0:         # MPW の磁極の端 (E1・E2 = 端の角の ANGLE に対する比) を edge 倍 (0 で端面なし)
        lines = src.split(chr(10))
        for i, l in enumerate(lines):
            if 'BMPW' in l and 'K0' in l:
                lines[i] = re.sub(r'(E[12]\s*=\s*)([0-9.]+)', lambda m: m.group(1) + repr(float(m.group(2)) * edge), l)
        src = chr(10).join(lines)
    assert n0 >= 4
    src = src.replace('K1 =%s*0.1800' % Q3, 'K1 =%r*0.1800' % (float(Q3) * (1 + float(a))))
    src = src.replace('K1 =%s*0.1800' % Q2, 'K1 =%r*0.1800' % (float(Q2) * (1 + float(b))))
    return src


def setup_of(src, line='RINGMPW'):
    from pysad.trksetup import TrackSetup
    st, _ = TrackSetup.from_source(src, line, rad=False, rfsw=False, alost=0.1)
    return st


def tunes(s, a, b, edge=1.0):
    st = setup_of(source(s, a, b, edge))
    return np.array(st.tunes_full())


def tune_match(s, target, a=0.0, b=0.0, it=6, edge=1.0):
    """QF3・QF2 の K1 の倍率 (1 + a)・(1 + b) を、全周のチューンが target (整数部つき) になるまで Newton 法で求める"""
    h = 2e-4
    for _ in range(it):
        t0 = tunes(s, a, b, edge)
        r = np.array(target) - t0
        if np.max(np.abs(r)) < 2e-6:
            break
        J = np.column_stack([(tunes(s, a + h, b, edge) - t0) / h, (tunes(s, a, b + h, edge) - t0) / h])
        da, db = np.linalg.solve(J, r)
        a, b = a + da, b + db
    return a, b, tunes(s, a, b, edge)


def cases():
    out = []
    for s in (0.0, 0.25, 0.5, 0.75, 1.0):
        out.append(dict(group='A', name='A_s%g' % s, s=s, target=None))
    for d in (-0.05, -0.03, -0.01, 0.01, 0.03, 0.05):
        out.append(dict(group='B', name='B_dnuy%+.2f' % d, s=1.0, target='mpw', dnuy=d))
    for s in (1.0, 0.5):
        out.append(dict(group='C', name='C_s%g_retuned' % s, s=s, target='base'))
    return out


def _init():
    from pysad import compat
    from pysad.trksetup import lower_priority
    compat.set_mode('v2')
    lower_priority()


def _task(c):
    from pysad import fmapcheck as FC
    from pysad.freqmap import frequency_map, resonance_check
    t0 = time.time()
    s = c['s']
    base = tunes(0.0, 0.0, 0.0)
    a = b = 0.0
    if c['target'] == 'base':
        a, b, nu = tune_match(s, base, edge=c.get('edge', 1.0))
    elif c['target'] == 'mpw':
        nu1 = tunes(1.0, 0.0, 0.0)
        a, b, nu = tune_match(1.0, (nu1[0], nu1[1] + c['dnuy']))
    else:
        nu = tunes(s, 0.0, 0.0)
    st = setup_of(source(s, a, b, c.get('edge', 1.0)))
    xs, ys = np.linspace(0.0005, 0.018, 21), np.linspace(0.0005, 0.014, 16)
    r = frequency_map(st, 2, x=xs, y=ys, nturn=1024, workers=1)
    stt = FC.map_stats(r, -5.0)
    lost = (r.lost_turn > 0).reshape(21, 16)
    near = [(m, n, p, round(d, 4)) for (m, n, p, d, o, nm) in resonance_check(nu[0], nu[1], 5, maxn=6)]
    out = dict(case=c, a=float(a), b=float(b), nu=[float(v) for v in nu], seconds=round(time.time() - t0, 1), stats=stt,
               lost_by_y=[int(v) for v in lost.sum(axis=0)], lost_by_x=[int(v) for v in lost.sum(axis=1)],
               d4nuy9=float(abs(4 * nu[1] - round(4 * nu[1]))) / 4, near=near,
               D=[round(float(v), 3) if np.isfinite(v) else None for v in r.D], lost=[int(v) for v in r.lost_turn])
    return out


def main():
    jobs = int(sys.argv[sys.argv.index('--jobs') + 1]) if '--jobs' in sys.argv else 3
    only = sys.argv[sys.argv.index('--only') + 1] if '--only' in sys.argv else None
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, 'S2_mpw.jsonl')
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding='utf-8'):
            done.add(json.loads(l)['case']['name'])
    cs = [c for c in cases() if c['name'] not in done and (only is None or c['group'] == only)]
    print('cases:', [c['name'] for c in cs], flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init) as ex:
        for n, d in enumerate(ex.map(_task, cs)):
            with open(out, 'a', encoding='utf-8') as f:
                f.write(json.dumps(d) + '\n')
            s = d['stats']
            print('[%d/%d %.0fs] %s: nu = (%.5f, %.5f)、|4νy − 9|/4 = %.4f、lost %d / %d、D>-5 %.3f、y 列の落ちた数 %s' % (
                n + 1, len(cs), time.time() - t0, d['case']['name'], d['nu'][0], d['nu'][1], d['d4nuy9'], s['nlost'], s['n'], s['frac_hi'], d['lost_by_y']), flush=True)


if __name__ == '__main__':
    main()
