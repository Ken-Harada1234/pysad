# -*- coding: utf-8 -*-
"""誤差つきのリングと MPW つきのリングの LMA。
  errors: ASTRID2 に位置のずれ・ロール (fmapres.apply_alignment_errors: 四極 50 µm・0.2 mrad、6 極 50 µm、偏向 50 µm・0.1 mrad の rms ×level) を、
          複数の種で与える。periodicity_report が周期性の崩れを検出して警告することを確かめ、全周の等間隔 25 位置で LMA を求め、誤差なしの全周と比べる。
          1 セルの結果を使った場合の誤差 (同じ位置の別のセルとの差) も見る。
  mpw:    ASTRID2_v01_02_id.sad の RINGMPW (MPW 入り)。周期性、全周 (等間隔 ~40 位置) の LMA、バケット・アパーチャの限界の種類。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_errors.py errors|mpw [workers] [level] [seeds 例 1,2,3]"""
import sys
import time
import warnings
import numpy as np
from lma_common import *


def sel(st, n):
    from pysad.trksetup import select_locations
    return select_locations(st, 1, None, 1, ds=st.circ / n)


def main(kind, workers=2, level=1.0, seeds=(1, 2, 3)):
    from pysad.lma import momentum_acceptance, periodicity_report, cell_locations
    from pysad import fmapres as FR
    from pysad import lmabucket as lb
    if kind == 'mpw':
        st, ffs = load('mpw')
        rep = periodicity_report(st, warn=False)
        print('MPW 版 ASTRID2: nper = %d、周期的 %s、メッセージ %s' % (rep['nper'], rep['periodic'], rep['messages']))
        b = st.bucket
        print('U0 = %.1f eV、V = %.0f V、δ_rf = %.5f、νs = %.5f、nturn(2 周期) = %d、周長 %.3f m' % (b['U0'], b['Vrf'], b['height'], b['nus'], st.nturn_for_synchrotron(2.0), st.circ))
        locs = sel(st, 40)
        t0 = time.time()
        r = momentum_acceptance(st, locs, workers=workers, quiet=True)
        print('全周 %d 位置: %.0f s、粒子 %d、警告 %s' % (len(locs), time.time() - t0, r.cond['tracked_particles'], r.warnings))
        lt = lb.limit_types(r.setup, r, workers=workers)
        print(lt.table())
        print('δ+: 最小 %.5f 最大 %.5f、δ−: 最小 %.5f 最大 %.5f。種類 δ+ %s' % (r.dp_plus.min(), r.dp_plus.max(), r.dp_minus.max(), r.dp_minus.min(),
                                                                             {k: lt.type_plus.count(k) for k in set(lt.type_plus)}))
        print('LMA − min(rf, aper) の最大: δ+ %.2e、δ− %.2e' % (np.abs(lt.dev_plus).max(), np.abs(lt.dev_minus).max()))
        r.save(os.path.join(HERE, 'lma_errors_mpw.lma.txt'))
        return
    st0, ffs = load('astrid2')
    locs = sel(st0, 25)
    print('ASTRID2、全周 %d 位置 (等間隔)、誤差 level %.2f' % (len(locs), level))
    t0 = time.time()
    r0 = momentum_acceptance(st0, locs, workers=workers, quiet=True, turncheck=0)
    print('誤差なし (全周): %.0f s' % (time.time() - t0))
    sig = dict(quad_dxy=50e-6 * level, quad_rot=0.2e-3 * level, sext_dxy=50e-6 * level, bend_dxy=50e-6 * level, bend_rot=0.1e-3 * level)
    for seed in seeds:
        st, cnt = FR.apply_alignment_errors(st0, seed=seed, **sig)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            rep = periodicity_report(st)
        print('\nseed %d: 誤差を与えた要素 %s。周期性: nper = %d、周期的 %s。警告 %d 件: %s' % (seed, cnt, rep['nper'], rep['periodic'], len(w), [str(x.message)[:90] for x in w]))
        locs_c, rep2 = cell_locations(st, warn=False)
        print('  cell_locations: %d 位置 (周期的でないので全周)' % len(locs_c))
        st.cod1 = st.refine_orbit(st.cod1)
        print('  閉軌道の残差 %.1e、cod1 (x, y) = (%.3e, %.3e) m' % (st.info.get('orbit_resid', float('nan')), st.cod1[0], st.cod1[2]))
        t0 = time.time()
        r = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0)
        d = r.dp_plus - r0.dp_plus
        dm = r.dp_minus - r0.dp_minus
        print('  LMA (全周): %.0f s、警告 %s' % (time.time() - t0, r.warnings))
        print('  誤差あり − なし: δ+ 平均 %+.2e 最小 %+.2e 最大 %+.2e、δ− 平均 %+.2e 最小 %+.2e 最大 %+.2e。δ+ の範囲 %.5f〜%.5f' % (
            d.mean(), d.min(), d.max(), dm.mean(), dm.min(), dm.max(), r.dp_plus.min(), r.dp_plus.max()))
        # 1 セルの結果を使うと: 誤差ありの位置 (最初のセル内) と、別のセルの同じ位置の差
        nper = 6
        starts = rep['starts'] if rep['nper'] > 1 else None
        r.save(os.path.join(HERE, 'lma_errors_seed%d.lma.txt' % seed))
    # 別のセルの同じ位置 (誤差なしは周期的、誤差ありは崩れる) の差: 誤差 seed 1 で、3 位置 × 6 セル
    from pysad.trksetup import find_periodicity
    per = find_periodicity(st0.els, st0.rows)
    st, _c = FR.apply_alignment_errors(st0, seed=seeds[0], **sig)
    st.cod1 = st.refine_orbit(st.cod1)
    first = [l for l in locs if per['starts'][0] <= l <= per['starts'][1] - 1]
    picks = [first[1], first[len(first) // 2], first[-2]]
    ll, tag = [], []
    for k in range(per['nper']):
        for l in picks:
            ll.append(per['starts'][k] + (l - per['starts'][0]))
            tag.append((k + 1, l))
    r = momentum_acceptance(st, ll, workers=workers, quiet=True, turncheck=0)
    print('\n誤差あり (seed %d): 同じ位置の別のセルとの差 (1 セルの結果を全周に延ばしたときの誤差の目安)' % seeds[0])
    for l in picks:
        v = [(r.dp_plus[i], r.dp_minus[i]) for i, t in enumerate(tag) if t[1] == l]
        print('  loc %d: δ+ %s  (幅 %.2e)  δ− (幅 %.2e)' % (l, ' '.join('%.5f' % x[0] for x in v), max(x[0] for x in v) - min(x[0] for x in v), max(x[1] for x in v) - min(x[1] for x in v)))


if __name__ == '__main__':
    kind = sys.argv[1] if len(sys.argv) > 1 else 'errors'
    main(kind, int(sys.argv[2]) if len(sys.argv) > 2 else 2, float(sys.argv[3]) if len(sys.argv) > 3 else 1.0,
         tuple(int(x) for x in sys.argv[4].split(',')) if len(sys.argv) > 4 else (1, 2, 3))
