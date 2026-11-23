# -*- coding: utf-8 -*-
"""結合 (線形のスキュー四極) だけで D > −5 の割合が上がるかの確かめ。
誤差なしの ASTRID2 に、**1 台の四極 (最初の QF3) のロール θ だけ**を入れる (他の誤差なし、閉軌道のずれなし)。θ は、フィードダウンの見積もり (fmapcheck.feeddown_estimate) の
結合の強さ |c−| が目標値 (0.0003、0.001、0.003、0.01) になるように決める。21 × 16 の周波数マップ (1024 周) の D > −5 の割合・失われた点・D の中央値を、
誤差なし (0.030) と比べ、差の共鳴までの距離 (数値) も出す。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_coupling.py [--jobs 2]     出力: 標準出力、coupling.jsonl
"""
import copy
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
TARGETS = (0.0003, 0.001, 0.003, 0.01)


def _task(target):
    from pysad import compat, fmapcheck as FC
    from pysad.freqmap import frequency_map
    from pysad.trksetup import TrackSetup, lower_priority
    lower_priority()
    compat.set_mode('v2')
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    st0 = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    iq = [i for i, e in enumerate(st0.els) if type(e).__name__ == 'Quadrupole'][0]

    def make(th):
        st = st0.copy()
        st.els = list(st0.els)
        e = copy.copy(st.els[iq])
        e.ROTATION = th
        st.els[iq] = e
        return st
    th = 1e-3
    c1 = FC.feeddown_estimate(make(th))['cminus_roll']
    th = th * target / abs(c1)
    st = make(th)
    es = FC.errors_summary(st)
    t0 = time.time()
    r = frequency_map(st, 2, x=np.linspace(0.0005, 0.018, 21), y=np.linspace(0.0005, 0.014, 16), nturn=1024, workers=1, warn=False)
    return dict(target=target, theta=th, errors=dict(cminus=float(es['cminus_pred']), split_num=float(es['split_num']), split_pred=float(es['split_pred'])),
                stats=FC.map_stats(r, -5.0), seconds=round(time.time() - t0, 1), D=[round(float(v), 3) if np.isfinite(v) else None for v in r.D])


def main():
    jobs = int(sys.argv[sys.argv.index('--jobs') + 1]) if '--jobs' in sys.argv else 2
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, 'coupling.jsonl')
    done = set()
    if os.path.exists(out):
        done = set(json.loads(l)['target'] for l in open(out, encoding='utf-8'))
    ts = [t for t in TARGETS if t not in done]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for d in ex.map(_task, ts):
            with open(out, 'a', encoding='utf-8') as f:
                f.write(json.dumps(d) + '\n')
    for l in open(out, encoding='utf-8'):
        d = json.loads(l)
        s = d['stats']
        print('|c−| の目標 %.4f: θ = %.3e rad (1 台の QF3 のロール)、|c−| の見積もり %.4f、差の共鳴までの距離 (見積もり %.4f、数値 %.4f)、D > −5: %.3f、失われた点 %d、D の中央値 %.2f' % (
            d['target'], d['theta'], d['errors']['cminus'], d['errors']['split_pred'], d['errors']['split_num'], s['frac_hi'], s['nlost'], s['Dmed']))


if __name__ == '__main__':
    main()
