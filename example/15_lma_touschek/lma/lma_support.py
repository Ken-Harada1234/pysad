# -*- coding: utf-8 -*-
"""LMA (局所運動量アクセプタンス) の FLUC・誤差・結合・k64 の対応範囲の確認。
  k64:   同じ ASTRID2・同じ位置で v2 と k64 の LMA (check_every = 1、同じ引数) を比べる。
  coupl: ASTRID2 に位置のずれ・ロールの誤差を level 倍 (既定 5 = ロール 1 mrad) で与え、結合 (差の共鳴からの距離と結合の見積もり) と LMA を誤差なしと比べる。
  kekls: kek-ls_v03_12.lat (大きなリング A) に誤差 (level を指定、seed 1) を与えた LMA (4 位置)。level 1 は閉軌道が求まらない (補正なし)。
使い方: PYSAD_SAD_PACKAGES=<...> PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_support.py k64|coupl|kekls [workers] [level] [nloc]
条件: RFSW・RAD ON、FLUC なし、LOSSAMPL = ±10 mm (例の値)、バケット判定あり (既定)。"""
import os
import sys
import time
import warnings

import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from lma_common import CASES


def load(name, mode):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode(mode)
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    return TrackSetup.from_source(src, cs['ring'], rad=True, rfsw=True, alost=cs['aper'])


def sel(st, n):
    from pysad.trksetup import select_locations
    return select_locations(st, 1, None, 1, ds=st.circ / n)


def P(*a):
    print(*a, flush=True)


def main(kind, workers=2, level=5.0, nloc=8):
    from pysad.lma import momentum_acceptance, periodicity_report
    from pysad import fmapres as FR
    if kind == 'k64':
        res = {}
        for mode in ('v2', 'k64'):
            st, ffs = load('astrid2', mode)
            locs = sel(st, nloc)
            t0 = time.time()
            r = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0)
            res[mode] = r
            P('%s: %d 位置 %.0f s、警告 %s' % (mode, len(locs), time.time() - t0, r.warnings))
        a, b = res['v2'], res['k64']
        P('v2 と k64 の δ± の差: δ+ 最大 %.2e、δ− 最大 %.2e (δ+ の範囲 %.5f〜%.5f)' % (np.abs(a.dp_plus - b.dp_plus).max(), np.abs(a.dp_minus - b.dp_minus).max(), a.dp_plus.min(), a.dp_plus.max()))
        for i in range(len(a)):
            P('  loc %4d s=%7.3f  δ+ v2 %.5f k64 %.5f | δ− v2 %.5f k64 %.5f' % (a.loc[i], a.s[i], a.dp_plus[i], b.dp_plus[i], a.dp_minus[i], b.dp_minus[i]))
        return
    name = 'astrid2' if kind == 'coupl' else 'kekls'
    lv = level
    st0, ffs = load(name, 'v2')
    locs = sel(st0, nloc)
    P('%s: %d 位置、誤差 level %.2f' % (CASES[name]['label'], len(locs), lv))
    t0 = time.time()
    r0 = momentum_acceptance(st0, locs, workers=workers, quiet=True, turncheck=0)
    P('誤差なし: %.0f s、警告 %s' % (time.time() - t0, r0.warnings))
    sig = dict(quad_dxy=50e-6 * lv, quad_rot=0.2e-3 * lv, sext_dxy=50e-6 * lv, bend_dxy=50e-6 * lv, bend_rot=0.1e-3 * lv)
    if kind == 'kekls':
        sig['oct_dxy'] = 50e-6 * lv
    st, cnt = FR.apply_alignment_errors(st0, seed=1, **{k: v for k, v in sig.items() if k in ('quad_dxy', 'quad_rot', 'sext_dxy', 'bend_dxy', 'bend_rot')})
    st.cod1 = st.refine_orbit(st.cod1)
    P('誤差を与えた要素 %s、閉軌道 cod1 (x, y) = (%.3e, %.3e) m、残差 %.1e' % (cnt, st.cod1[0], st.cod1[2], st.info.get('orbit_resid', float('nan'))))
    try:
        from pysad import fmapcheck as FC
        w = FC.check_setup(st) if hasattr(FC, 'check_setup') else None
        P('fmapcheck:', w if w is None else [str(x)[:100] for x in w][:5])
    except Exception as e:
        P('fmapcheck 不可: %r' % (e,))
    t0 = time.time()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        r = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0)
    P('誤差あり: %.0f s、警告 %s' % (time.time() - t0, [str(x.message)[:100] for x in w]))
    P('δ+ の差 (あり − なし): 平均 %+.2e 最小 %+.2e 最大 %+.2e、δ− 平均 %+.2e。δ+ の範囲 %.5f〜%.5f' % (
        (r.dp_plus - r0.dp_plus).mean(), (r.dp_plus - r0.dp_plus).min(), (r.dp_plus - r0.dp_plus).max(), (r.dp_minus - r0.dp_minus).mean(), r.dp_plus.min(), r.dp_plus.max()))
    for i in range(len(r)):
        P('  loc %4d s=%7.3f δ+ %.5f → %.5f | δ− %.5f → %.5f' % (r.loc[i], r.s[i], r0.dp_plus[i], r.dp_plus[i], r0.dp_minus[i], r.dp_minus[i]))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0] if a else 'k64', int(a[1]) if len(a) > 1 else 2, float(a[2]) if len(a) > 2 else 5.0, int(a[3]) if len(a) > 3 else 8)
