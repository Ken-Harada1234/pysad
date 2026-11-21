# -*- coding: utf-8 -*-
"""v2.4.0 A5: 周波数マップの予想 (S3) が ASTRID2・KEK-LS で改善しない原因の切り分け: 保存済みの周波数マップ (Q1_fmap_*.table.txt、41×31 の格子、1024 周) の各点のチューン
(追跡の前半・後半 NAFF の平均) を、ampnf の ν(Jx, Jy) (正規形の級数、J¹ だけ / J² まで) と比べる。
  J の決め方: 格子の初期条件 (x, 0, y, 0) を ampnf.nf_action で正規形のアクションにする (追跡の初期条件の J と正規形の J の違いを補正)。
  ampnf の模型: 線形 + 6 極・8 極 + 偏向磁石の曲率・端 + 運動学的な項 + 端部のフリンジ (bend_cubic・kin・fringe)。位置 loc から始まるリング (要素の並びを回す)。
出力: 振幅の層ごとの |Δν| の中央値 (J¹ のみ、J² まで)、D < −5 の規則的な点だけの値、J² を入れた効果。
使い方: PYSAD_COMPAT=k64 PYTHONUTF8=1 py -3.12 example/v2_4/A5_nu_of_J_grid.py astrid2|kekls [D]"""
import math
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'v2')
from pysad import ampnf as A, compat
from pysad.parser import load_sad_file

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), ring='RING', loc=2, table='Q1_fmap_astrid2.table.txt'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), ring='RING', loc=3, table='Q1_fmap_kekls.table.txt'),
}


class Rot:
    def __init__(self, elements, k):
        self.elements = list(elements[k:]) + list(elements[:k])


def fr(a):
    return (np.asarray(a) + 0.5) % 1.0 - 0.5


def main(name, D=5, skip_den=0.0):
    cs = CASES[name]
    compat.set_mode('v2')
    lat = load_sad_file(cs['path']).build_lattice(cs['ring'])
    rot = Rot(lat.elements, cs['loc'] - 1)
    t0 = time.time()
    r = A.amplitude_tune(rot, D=D, nslice=8, nslice_bend=8, bend_cubic=True, kin=True, fringe=True, skip_den=skip_den)
    print('%s: ampnf D=%d、%.0f s、ν0 = (%.5f, %.5f)、最小の共鳴の分母 %.3g' % (name, D, time.time() - t0, r.nu0[0], r.nu0[1], r.min_den), flush=True)
    d = np.loadtxt(os.path.join(HERE, '..', 'v2_3', cs['table']))
    x, y, nx1, ny1, nx2, ny2, Dd, lost = d[:, 0], d[:, 1], d[:, 3], d[:, 4], d[:, 5], d[:, 6], d[:, 7], d[:, 8]
    ok = lost == 0
    nxt = 0.5 * (nx1 + nx2)
    nyt = 0.5 * (ny1 + ny2)
    rows = []
    pred = []
    for i in np.nonzero(ok)[0]:
        jx, jy, _, _, err = A.nf_action(r, (x[i], 0.0, y[i], 0.0), 0.0)
        if err > 1e-6:
            rows.append((np.nan,) * 6)
            pred.append((np.nan,) * 4)
            continue
        out = []
        for pl, tr in ((0, nxt[i]), (1, nyt[i])):
            j1 = r.nu0[pl] + sum(r.coef(pl, a, b, 0) * jx ** a * jy ** b for a, b in ((1, 0), (0, 1)))
            j2 = r.nu(pl, jx, jy, 0.0)
            out += [abs(fr(tr - j1)), abs(fr(tr - j2))]
        rows.append((jx, jy) + tuple(out))
        pred.append((r.nu(0, jx, jy, 0.0), r.nu(1, jx, jy, 0.0), r.nu0[0] + r.coef(0, 1, 0) * jx + r.coef(0, 0, 1) * jy, r.nu0[1] + r.coef(1, 1, 0) * jx + r.coef(1, 0, 1) * jy))
    arr = np.array(rows)
    idx = np.nonzero(ok)[0]
    jt = arr[:, 0] + arr[:, 1]
    good = np.isfinite(jt) & (Dd[idx] < -5.0)
    print('点: 失われない %d / %d、nf_action が収束した点 %d、規則的 (D < −5) %d' % (ok.sum(), len(d), np.isfinite(jt).sum(), good.sum()))
    qs = np.nanquantile(jt, [0.0, 1 / 3, 2 / 3, 1.0])
    print('振幅の層 (J = Jx + Jy [m rad]): 小 %.1e〜%.1e、中 〜%.1e、大 〜%.1e' % (qs[0], qs[1], qs[2], qs[3]))
    print('|Δν| の中央値 (追跡 − 級数): 列 = νx の J¹ のみ / J² まで、νy の J¹ のみ / J² まで   (規則的な点 D < −5 だけ | 全点)')
    for k in range(3):
        sel = (jt >= qs[k]) & (jt <= qs[k + 1]) & np.isfinite(jt)
        s2 = sel & good
        print('  層 %s: 規則的 %4d 点  %.2e %.2e | %.2e %.2e    全点 %4d 点  %.2e %.2e | %.2e %.2e' % (
            ('小', '中', '大')[k], s2.sum(), *[np.nanmedian(arr[s2, c]) for c in (2, 3, 4, 5)], sel.sum(), *[np.nanmedian(arr[sel, c]) for c in (2, 3, 4, 5)]))
    print('J² を入れた効果 (規則的な点、中央値の比 J² まで / J¹ のみ): νx %.2f、νy %.2f' % (np.nanmedian(arr[good, 3]) / np.nanmedian(arr[good, 2]), np.nanmedian(arr[good, 5]) / np.nanmedian(arr[good, 4])))


    # v2.4.0 A9: 最小の共鳴の分母の項の特定。ASTRID2 は ζ̄x³ ζ̄y² など (4νx + 2νy = 25 の 6 次の和の共鳴、距離 3.2e-4) の項が最小 (|λ^α − λ_i| = 0.002)。
    # 追跡のチューンから共鳴までの距離 |4νx + 2νy − 25| と |Δν| の関係を見る
    if name == 'astrid2':
        from scipy.stats import spearmanr
        dist = np.abs(4 * (5 + nxt[idx]) + 2 * (2 + nyt[idx]) - 25.0)
        ok2 = np.isfinite(jt) & good
        print('共鳴 4νx + 2νy = 25 までの距離 (追跡のチューン) と |Δν| (J² まで/級数全体) の Spearman: νx %.2f、νy %.2f' % (
            spearmanr(dist[ok2], arr[ok2, 3])[0], spearmanr(dist[ok2], arr[ok2, 5])[0]))
        for lo, hi in ((0, 0.005), (0.005, 0.02), (0.02, 0.1), (0.1, 9)):
            s = ok2 & (dist >= lo) & (dist < hi)
            if s.sum():
                print('  距離 %.3f〜%.3f: %4d 点、|Δν| の中央値 νx %.2e νy %.2e' % (lo, hi, s.sum(), np.nanmedian(arr[s, 3]), np.nanmedian(arr[s, 5])))

        # 差の共鳴 νx − νy = 3 (線形の作動点から 0.04) にロックされた点 (追跡のチューンが |νx − νy − 3| < 0.002) と、そうでない点
        d2 = np.abs((5 + nxt[idx]) - (2 + nyt[idx]) - 3.0)
        for lab, s in (('ロック (|νx−νy−3| < 0.002)', ok2 & (d2 < 0.002)), ('近い (0.002〜0.01)', ok2 & (d2 >= 0.002) & (d2 < 0.01)), ('それ以外', ok2 & (d2 >= 0.01))):
            if s.sum():
                print('  差の共鳴 %s: %4d 点、|Δν| の中央値 νx %.2e νy %.2e (J の中央値 %.2e)' % (lab, s.sum(), np.nanmedian(arr[s, 3]), np.nanmedian(arr[s, 5]), np.nanmedian(jt[s])))

    # v2.4.0 A13: ロックした点の識別を共鳴の島の位置の予想に使えるか。線形の作動点に最も近い 2 次の和・差の共鳴 (mνx + nνy = p、|m|=|n|=1 または 2) を選び、
    # 級数のチューン (J¹ のみ / J² まで) からその共鳴までの距離で、追跡のチューンが共鳴にロックした点 (距離 < 0.002) を当てられるか (AUC) を見る
    from scipy.stats import rankdata
    tf = [float(v) for v in open(os.path.join(HERE, '..', 'v2_3', cs['table']), encoding='utf-8').read().split('# tunes_full = (')[1].split(')')[0].split(',')]
    best = None
    for m in (1, 2):
        for n in (-2, -1, 1, 2):
            if abs(n) > m: continue
            d0 = m * tf[0] + n * tf[1]
            p0 = round(d0)
            if best is None or abs(d0 - p0) < best[0]:
                best = (abs(d0 - p0), m, n, p0)
    _, m, n, p0 = best
    pr = np.array(pred)
    ip = np.arange(len(idx))
    trk = np.abs(m * (math.floor(tf[0]) + nxt[idx]) + n * (math.floor(tf[1]) + nyt[idx]) - p0)
    lock = (trk < 0.002) & np.isfinite(arr[:, 0])
    print('共鳴 %d νx %+d νy = %d (線形の作動点から %.4f)。ロックした点 (追跡の距離 < 0.002) %d / %d' % (m, n, p0, best[0], lock.sum(), np.isfinite(arr[:, 0]).sum()))
    def auc(score, lab):
        ok = np.isfinite(score)
        r = rankdata(score[ok]); l = lab[ok]
        return (r[l].sum() - l.sum() * (l.sum() + 1) / 2) / (l.sum() * (~l).sum()) if l.sum() and (~l).sum() else float('nan')
    for lab, (px, py) in (('線形の ADTS (J¹ のみ)', (pr[:, 2], pr[:, 3])), ('J² まで (級数全体)', (pr[:, 0], pr[:, 1]))):
        dist = np.abs(m * (math.floor(tf[0]) + px) + n * (math.floor(tf[1]) + py) - p0)
        # ロックした点では追跡のチューンが線上にある。予想の距離が小さいほど「ロック」と予想する: スコア = −距離
        print('  %s: 予想の距離の中央値 ロック %.4f、それ以外 %.4f、AUC (−距離でロックを識別) %.3f' % (lab, np.nanmedian(dist[lock]), np.nanmedian(dist[~lock & np.isfinite(dist)]), auc(-dist, lock)))
    jt2 = arr[:, 0] + arr[:, 1]
    print('  参考: 振幅 J だけのスコア (J が大きいほどロック) の AUC %.3f' % auc(jt2, lock))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 5, float(sys.argv[3]) if len(sys.argv) > 3 else 0.0)
