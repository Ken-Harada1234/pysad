# -*- coding: utf-8 -*-
"""v2.3.0 S2 (1): 誤差つき周波数マップの統計 (複数の乱数の種 × 複数の強さ × 誤差の種類の切り分け)。
各タスク = (リング, 種, 強さ k, 誤差の種類) の周波数マップ (小さな格子 21 × 16、1024 周) を走らせ、統計を JSON の 1 行で S2_data/ の .jsonl に追記する (走った組は飛ばす)。
誤差の標準偏差 (k = 1): 四極・6 極・偏向磁石の DX・DY 50 μm (KEK-LS は 8 極も)、四極・6 極のロール 0.2 mrad (6 極は 'all' のときだけ入れず、'sext_rot' の切り分けだけ)、偏向磁石のロール 0.1 mrad。
種類: all_cs (all のあと 6 極 (KEK-LS は 8 極も) のずれを閉軌道に合わせてフィードダウンを切る)、all (quad_dxy + quad_rot + sext_dxy + bend_dxy + bend_rot [+ oct_dxy])、quad_dxy、quad_rot、sext_dxy、bend_dxy、bend_rot、oct_dxy、none (k = 0)。
乱数は pysad.fmapres.apply_alignment_errors(streams='type') (種類ごとに別の列。種類を 1 つにしても、全部のときのその種類の誤差と同じ)。
使い方 (リポジトリ直下から): py -3.12 example/v2_3/S2_errors_scan.py astrid2 --seeds 1-8 --levels 0.3,1,3 --kinds all --jobs 3
                              py -3.12 example/v2_3/S2_errors_scan.py astrid2 --seeds 1-5 --levels 1,3 --kinds quad_dxy,sext_dxy,bend_rot,quad_rot,bend_dxy --jobs 3
条件: v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、px = py = z = 0 から出発、NAFF の窓 Hann、規格化 = 標本の楕円。ASTRID2 は位置 2、x = 0.5〜18 mm・y = 0.5〜14 mm、
KEK-LS (pysad/tests/kek-ls_v03_12.lat) は位置 3、x = 0.2〜6 mm・y = 0.2〜8 mm。格子は 21 × 16 (S1 の full の 1 つおき、336 点)。
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, 'S2_data')

RINGS = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), line='RING', loc=2, x=(0.0005, 0.018), y=(0.0005, 0.014)),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), line='RING', loc=3, x=(0.0002, 0.006), y=(0.0002, 0.008)),
}
SIG1 = dict(quad_dxy=50e-6, quad_rot=0.2e-3, sext_dxy=50e-6, bend_dxy=50e-6, bend_rot=0.1e-3, sext_rot=0.2e-3, oct_dxy=50e-6)
ALL = ('quad_dxy', 'quad_rot', 'sext_dxy', 'bend_dxy', 'bend_rot')
_W = {}


def kind_sigmas(ring, kind, k):
    if kind == 'none':
        return {}
    if kind.endswith('_cs'):
        kind = kind[:-3]
    if kind == 'all':
        ks = list(ALL) + (['oct_dxy'] if ring == 'kekls' else [])
    else:
        ks = kind.split('+')
    return {a: SIG1[a] * k for a in ks}


def _init(ring):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    c = RINGS[ring]
    src = open(c['path'], encoding='utf-8').read()
    st, _ = TrackSetup.from_source(src, c['line'], rad=False, rfsw=False, alost=0.1)
    _W['st'] = st
    _W['ring'] = ring
    from pysad.trksetup import lower_priority
    lower_priority()


def _task(a):
    ring, seed, k, kind = a
    from pysad import fmapres as FR, fmapcheck as FC
    from pysad.freqmap import frequency_map
    c = RINGS[ring]
    st0 = _W['st']
    sig = kind_sigmas(ring, kind, k)
    t0 = time.time()
    st, cnt = FR.apply_alignment_errors(st0, seed=seed, streams='type', **sig)
    if kind.endswith('_cs'):            # 6 極 (KEK-LS は 8 極も) を閉軌道に合わせる = フィードダウンを切る切り分け
        st = FC.center_sextupoles(st, kinds=('Sextupole', 'Octupole') if ring == 'kekls' else ('Sextupole',))
    es = FC.errors_summary(st, loc=c['loc'])
    xs = np.linspace(*c['x'], 21)
    ys = np.linspace(*c['y'], 16)
    r = frequency_map(st, c['loc'], x=xs, y=ys, nturn=1024, workers=1)
    s = FC.map_stats(r, -5.0)
    s4 = FC.map_stats(r, -4.0)
    out = dict(ring=ring, seed=seed, k=k, kind=kind, sigma=sig, count=cnt, seconds=round(time.time() - t0, 1),
               stats=s, frac_hi4=s4['frac_hi'], D=[round(float(v), 3) if np.isfinite(v) else None for v in r.D], lost=[int(v) for v in r.lost_turn],
               nux1=[round(float(v), 6) if np.isfinite(v) else None for v in r.nux1], nuy1=[round(float(v), 6) if np.isfinite(v) else None for v in r.nuy1],
               D0=float(r.D[0]) if np.isfinite(r.D[0]) else None, errors=_clean(es))
    return out


def _clean(d):
    o = {}
    for k, v in d.items():
        if isinstance(v, (tuple, list)):
            o[k] = [float(t) for t in v]
        elif isinstance(v, (bool, np.bool_)):
            o[k] = bool(v)
        else:
            o[k] = float(v)
    return o


def parse_range(s):
    out = []
    for p in s.split(','):
        if '-' in p:
            a, b = p.split('-')
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(p))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ring')
    ap.add_argument('--seeds', default='1-5')
    ap.add_argument('--levels', default='1')
    ap.add_argument('--kinds', default='all')
    ap.add_argument('--jobs', type=int, default=3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    out = a.out or os.path.join(DATA, 'S2_errors_%s.jsonl' % a.ring)
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding='utf-8'):
            d = json.loads(l)
            done.add((d['seed'], d['k'], d['kind']))
    tasks = []
    for kind in a.kinds.split(','):
        for k in [float(v) for v in a.levels.split(',')]:
            for seed in (parse_range(a.seeds) if kind != 'none' else [0]):
                if (seed, k, kind) not in done:
                    tasks.append((a.ring, seed, k, kind))
    print('tasks:', len(tasks), flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init, initargs=(a.ring,)) as ex:
        for n, d in enumerate(ex.map(_task, tasks)):
            with open(out, 'a', encoding='utf-8') as f:
                f.write(json.dumps(d) + '\n')
            s = d['stats']
            print('[%d/%d %.0fs] %s seed %d k %g %s: lost %d, D>-5 %.3f, Dmed %.2f, D0 %s, split_pred %.4f split_num %.4f cminus %.4f' % (
                n + 1, len(tasks), time.time() - t0, a.ring, d['seed'], d['k'], d['kind'], s['nlost'], s['frac_hi'], s['Dmed'], d['D0'],
                d['errors']['split_pred'], d['errors']['split_num'], d['errors']['cminus_pred']), flush=True)


if __name__ == '__main__':
    main()
