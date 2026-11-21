# -*- coding: utf-8 -*-
"""v2.3.0 L2 項目 2: RF バケットの限界とアパーチャの限界の区別 (pysad/lmabucket.py)。
使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/v2_3/L2_limit_types.py astrid2|kekls [workers]
  1. 縦方向の写像のパラメーター (α1・α2、放射の δ 依存、空洞の力の確認) と、解析式のバケットの高さ δ_rf、Emittance[] の値
  2. 1 次元の写像・6 次元の追跡 (アパーチャ無限大) のバケットの限界 (1 つの位置、周回数を変えて)
  3. 位置ごとに LMA (既定の check_every = 1)・バケットの限界 (6 次元)・アパーチャの限界 (4 次元) を比べ、種類を決める
"""
import sys
import time
import numpy as np
from L2_common import *


def main(name, workers=2):
    from pysad.lma import momentum_acceptance
    from pysad import lmabucket as lb
    from pysad.trksetup import find_periodicity, select_locations
    st, ffs = load(name)
    cs = CASES[name]
    print(cs['label'], 'LOSSAMPL = %.3f m' % cs['aper'])
    an = lb.rf_bucket_analytic(st)
    print('解析式 δ_rf = %.5f (Emittance[] の bucket_height %.5f)、h = %.2f、バケットの |z| の限界 %.4f m' % (an['delta_rf'], an['emit_bucket_height'], an['h'], an['zlimit']))
    par = lb.longitudinal_params(st)
    print('縦の写像 (空洞の後ろの loc %d): a1 = %.6f (αc の %.4f 倍)、a2 = %.5f (a2/a1 = %.3f)、放射 r1 = %.4e r2 = %.3e (理想 −2U0/E0 = %.3e、−U0/E0 = %.3e)、'
          '空洞の力 (V/E0) sin(φs+kz) − U0/E0 の追跡との最大の差 %.2e (V/E0 = %.3e)' % (
              par['loc'], par['a1'], par['a1_over_alpha'], par['a2'], par['a2'] / par['a1'], par['r1'], par['r2'],
              -2 * st.bucket['U0'] / st.p0, -st.bucket['U0'] / st.p0, par['kick_max_diff'], par['kick_scale']))
    nturn = st.nturn_for_synchrotron(2.0)
    print('周回数 (2 周期) %d、νs = %.5f' % (nturn, st.bucket['nus']))
    # 1 次元の写像の項の寄与と周回数への依存
    print('\n--- 1 次元の写像のバケットの限界 (δ+, δ−) ---')
    for N in sorted(set([nturn // 4, nturn // 2, nturn, 2 * nturn, 4 * nturn])):
        m = lb.bucket_limits_1d(st, N, par=par)
        print('nturn = %6d (%.2f 周期): %s %s' % (N, N * abs(st.bucket['nus']), m['dp_plus'], m['dp_minus']))
    for nm, kw in (('a2 = 0', dict(a2=0.0)), ('放射の δ 依存なし (r1 = r2 = 0)', dict(r_on=False)), ('a2 = 0 かつ r なし', dict(a2=0.0, r_on=False))):
        m = lb.bucket_limits_1d(st, nturn, par=par, **kw)
        print('nturn = %6d、%s: %.5f %.5f' % (nturn, nm, m['dp_plus'], m['dp_minus']))
    # 位置ごと
    per = find_periodicity(st.els, st.rows)
    locs = select_locations(st, per['starts'][0], per['starts'][1] - 1, 1, cs['ds'])
    if name == 'kekls':
        locs = locs[::4]
    print('\n位置の数', len(locs))
    t0 = time.time()
    r = momentum_acceptance(st, locs, workers=workers, quiet=True)
    print('LMA (既定 check_every = %d): %.0f s、追跡した粒子 %d' % (r.cond['check_every'], time.time() - t0, r.cond['tracked_particles']))
    t0 = time.time()
    lt = lb.limit_types(st.copy(check_every=r.cond['check_every']), r, workers=workers, nturn_aper=256)
    print('バケット (6D)・アパーチャ (4D、256 周): %.0f s' % (time.time() - t0))
    print(lt.table())
    for sg in ('plus', 'minus'):
        ty = getattr(lt, 'type_' + sg)
        dev = np.abs(getattr(lt, 'dev_' + sg))
        print('δ%s: 種類 rf %d、aper %d、both %d。LMA − min(rf, aper) の最大 %.2e、rms %.2e' % (
            '+' if sg == 'plus' else '−', ty.count('rf'), ty.count('aper'), ty.count('both'), dev.max(), np.sqrt((dev ** 2).mean())))
    print('1 次元の写像 (nturn = %d): δ+ = %.5f、δ− = %.5f。6 次元 (アパーチャ無限大) の δ+ の範囲 %.5f〜%.5f、|δ−| の範囲 %.5f〜%.5f' % (
        lt.nturn, lt.map1d['dp_plus'], lt.map1d['dp_minus'], lt.rf_plus.min(), lt.rf_plus.max(), -lt.rf_minus.max(), -lt.rf_minus.min()))
    print('LMA の損失の直接の理由 (δ+): ' + ', '.join('%s %d' % (k, lt.reason_plus.count(k)) for k in ('rf', 'aper', 'sat')))
    r.save(os.path.join(HERE, 'L2_limit_types_%s.lma.txt' % name))
    open(os.path.join(HERE, 'L2_limit_types_%s.table.txt' % name), 'w', encoding='utf-8').write(lt.table())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 2)
