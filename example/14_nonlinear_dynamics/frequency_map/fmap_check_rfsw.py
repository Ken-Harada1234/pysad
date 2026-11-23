# -*- coding: utf-8 -*-
"""RFSW ON・RAD ON・FLUC ON の周波数マップ (ASTRID2、位置 2)。既定 (RFSW・RAD OFF) と比べる。
ケース: base (NORFSW・NORAD)、rf (RFSW ON・RAD OFF)、rfrad (RFSW ON・RAD ON)、rfradfluc (RFSW・RAD・FLUC ON、乱数の種 17 と 29)、radfluc (RFSW OFF・RAD ON・FLUC ON)。
各ケースの 21 × 16 の周波数マップ (1024 周、Hann 窓、規格化 = 標本の楕円、LOSSAMPL = 0.1 m、LOSSDZ は SAD の既定 100 m ではなく RF バケットの半分 (RF ON のケース)) を走らせ、
base との差 (D の相関・差・失われた点・チューンのずれ) と、代表の 4 点の z・δ の振幅を出す。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_rfsw.py [case ...] [--jobs 3]     出力: 標準出力、rfsw.jsonl
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
CASES = {
    'base': dict(rfsw=False, rad=False, fluc=False, seed=None),
    'rf': dict(rfsw=True, rad=False, fluc=False, seed=None),
    'rfrad': dict(rfsw=True, rad=True, fluc=False, seed=None),
    'rfradfluc17': dict(rfsw=True, rad=True, fluc=True, seed=17),
    'rfradfluc29': dict(rfsw=True, rad=True, fluc=True, seed=29),
    'radfluc17': dict(rfsw=False, rad=True, fluc=True, seed=17),
}
PTS = [(1e-3, 1e-3), (6e-3, 3e-3), (10e-3, 6e-3), (14e-3, 9e-3)]


def _init():
    from pysad import compat
    from pysad.trksetup import lower_priority
    compat.set_mode('v2')
    lower_priority()


def _task(name):
    from pysad import fmapcheck as FC
    from pysad.freqmap import frequency_map
    from pysad.trksetup import TrackSetup
    c = CASES[name]
    t0 = time.time()
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    st, _ = TrackSetup.from_source(src, 'RING', rad=c['rad'], rfsw=c['rfsw'], fluc=c['fluc'], alost=0.1)
    if c['seed'] is not None:
        st.seed = c['seed']
    if st.rfsw:
        zl = st.bucket_zlimit()
        st.zlost = zl if zl else st.zlost
    xs, ys = np.linspace(0.0005, 0.018, 21), np.linspace(0.0005, 0.014, 16)
    import warnings
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        r = frequency_map(st, 2, x=xs, y=ys, nturn=1024, workers=1)
    # 代表の 4 点の z・δ
    base = st.closed_orbit([2])[0] if st.rfsw else np.zeros(6)
    co = np.array([base + np.array([x, 0, y, 0, 0, 0]) for (x, y) in PTS])
    tr = st.track(2, co, 1024, record=True)
    rec = tr['rec']
    zs = [[float(np.nanmax(np.abs(rec[:, i, 4]))), float(np.nanmax(np.abs(rec[:, i, 5]))), float(np.nanmean(rec[-64:, i, 5]))] for i in range(len(PTS))]
    out = dict(case=name, cfg=c, seconds=round(time.time() - t0, 1), stats=FC.map_stats(r, -5.0), cond=dict(
        rfsw=st.rfsw, rad=st.rad, fluc=st.fluc, bucket=st.bucket, zlost=st.zlost), zdp=zs, warnings=[str(w.message)[:200] for w in wl],
        D=[round(float(v), 3) if np.isfinite(v) else None for v in r.D], lost=[int(v) for v in r.lost_turn],
        nux1=[round(float(v), 7) if np.isfinite(v) else None for v in r.nux1], nuy1=[round(float(v), 7) if np.isfinite(v) else None for v in r.nuy1],
        nux2=[round(float(v), 7) if np.isfinite(v) else None for v in r.nux2], nuy2=[round(float(v), 7) if np.isfinite(v) else None for v in r.nuy2])
    return out


def main():
    jobs = int(sys.argv[sys.argv.index('--jobs') + 1]) if '--jobs' in sys.argv else 3
    names = [a for a in sys.argv[1:] if a in CASES] or list(CASES)
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, 'rfsw.jsonl')
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding='utf-8'):
            done.add(json.loads(l)['case'])
    names = [n for n in names if n not in done]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init) as ex:
        for n, d in enumerate(ex.map(_task, names)):
            with open(out, 'a', encoding='utf-8') as f:
                f.write(json.dumps(d) + '\n')
            s = d['stats']
            print('[%d/%d %.0fs] %s: lost %d / %d、D>-5 %.3f、D 中央値 %.2f、(max|z|, max|δ|, 末尾の <δ>) %s、警告 %d' % (
                n + 1, len(names), time.time() - t0, d['case'], s['nlost'], s['n'], s['frac_hi'], s['Dmed'],
                [[round(v, 5) for v in z] for z in d['zdp']], len(d['warnings'])), flush=True)


if __name__ == '__main__':
    main()
