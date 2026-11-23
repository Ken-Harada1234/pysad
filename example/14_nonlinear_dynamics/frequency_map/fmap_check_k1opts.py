# -*- coding: utf-8 -*-
"""v2 と k64 で動作が違う設定での周波数マップ (ASTRID2、21 × 16、1024 周)。
ケース: default (v2、radial・log)、face (compat.set_bend_k1_edge('face'): K1 ≠ 0 の BEND の四極の磁場も端面から始まる)、lin (compat.set_bend_k1_def('lin'): B' 一様、
elegant・MAD と同じ定義)、face_lin (両方)、k64 (互換モード k64)。設定は TrackSetup を作る前に切り替える (setup が edge・def の値を持ち、frequency_map が合わせる)。
default と比べる: 作動点、D の相関、|ΔD|、D > −5 の割合、失われた点、チューンの差 (最小振幅・全点)。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_k1opts.py [--jobs 2]    出力: 標準出力、k1opts.jsonl
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, 'fmap_check_data')
CASES = {'default': ('v2', 'radial', 'log'), 'face': ('v2', 'face', 'log'), 'lin': ('v2', 'radial', 'lin'), 'face_lin': ('v2', 'face', 'lin'), 'k64': ('k64', 'radial', 'log'),
         # rect: 偏向磁石 BH の端の角を E1 = E2 = 0.02 (ANGLE に対する比。端の角 0.0105 rad = 0.6°。0.1 以上は作動点が動いて全部失われる) にした変形。ASTRID2 の BH は K1 ≠ 0 だが E1 = E2 = 0 なので 'face' は効かない (端の角 ≠ 0 のときだけ効く)
         'rect_default': ('v2', 'radial', 'log'), 'rect_face': ('v2', 'face', 'log'), 'rect_lin': ('v2', 'radial', 'lin'), 'rect_face_lin': ('v2', 'face', 'lin'), 'rect_k64': ('k64', 'radial', 'log')}


def _task(name):
    from pysad import compat, fmapcheck as FC
    from pysad.freqmap import frequency_map
    from pysad.trksetup import TrackSetup, lower_priority
    lower_priority()
    mode, edge, kdef = CASES[name]
    compat.set_mode(mode)
    compat.set_bend_k1_edge(edge)
    compat.set_bend_k1_def(kdef)
    t0 = time.time()
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    if name.startswith('rect'):
        assert 'K1 = -1.678200*0.8500  E1=0 E2=0' in src
        src = src.replace('K1 = -1.678200*0.8500  E1=0 E2=0', 'K1 = -1.678200*0.8500  E1=0.02 E2=0.02')
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    xs, ys = np.linspace(0.0005, 0.018, 21), np.linspace(0.0005, 0.014, 16)
    r = frequency_map(st, 2, x=xs, y=ys, nturn=1024, workers=1)
    return dict(case=name, cfg=CASES[name], nu_lin=[float(v) for v in st.tunes_full()], seconds=round(time.time() - t0, 1), stats=FC.map_stats(r, -5.0),
                setup_k1=[st.bend_k1_edge, st.bend_k1_def], D=[round(float(v), 4) if np.isfinite(v) else None for v in r.D], lost=[int(v) for v in r.lost_turn],
                nux1=[float(v) for v in r.nux1], nuy1=[float(v) for v in r.nuy1])


def main():
    jobs = int(sys.argv[sys.argv.index('--jobs') + 1]) if '--jobs' in sys.argv else 2
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, 'k1opts.jsonl')
    done = set()
    if os.path.exists(out):
        done = set(json.loads(l)['case'] for l in open(out, encoding='utf-8'))
    names = [n for n in CASES if n not in done]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for d in ex.map(_task, names):
            with open(out, 'a', encoding='utf-8') as f:
                f.write(json.dumps(d) + '\n')
            s = d['stats']
            print('%s: 作動点 (%.6f, %.6f)、lost %d、D>-5 %.3f、D 中央値 %.2f、%s' % (d['case'], d['nu_lin'][0], d['nu_lin'][1], s['nlost'], s['frac_hi'], s['Dmed'], d['setup_k1']), flush=True)
    rows = {json.loads(l)['case']: json.loads(l) for l in open(out, encoding='utf-8')}
    for ref in ('default', 'rect_default'):
        if ref not in rows:
            continue
        b = rows[ref]
        Db = np.array([np.nan if v is None else v for v in b['D']])
        for n, d in rows.items():
            if n == ref or (n.startswith('rect') != ref.startswith('rect')):
                continue
            Dn = np.array([np.nan if v is None else v for v in d['D']])
            ok = np.isfinite(Db) & np.isfinite(Dn)
            dnu = np.hypot(np.array(d['nux1']) - np.array(b['nux1']), np.array(d['nuy1']) - np.array(b['nuy1']))
            print('%s と %s: 作動点の差 (%.2e, %.2e)、D の相関 %.3f、|ΔD| 中央値 %.3f・最大 %.2f、失われた点 %d 対 %d、最小振幅のチューン差 %.2e、チューン差の中央値 %.2e' % (
                n, ref, d['nu_lin'][0] - b['nu_lin'][0], d['nu_lin'][1] - b['nu_lin'][1], np.corrcoef(Db[ok], Dn[ok])[0, 1], np.median(np.abs(Db[ok] - Dn[ok])),
                np.max(np.abs(Db[ok] - Dn[ok])), sum(1 for v in d['lost'] if v), sum(1 for v in b['lost'] if v), dnu[0], np.nanmedian(dnu)))


if __name__ == '__main__':
    main()
