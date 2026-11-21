# -*- coding: utf-8 -*-
"""v2.3.0 S2 (5): D の周回数 (512〜4096)・NAFF の窓・格子の細かさへの依存。
1 回の追跡 (4096 周、周ごとに記録) から、最初の N 周 (N = 256, 512, 1024, 2048, 4096。前半・後半 N/2 周) の D_N を作り、点ごとの D_N と N の関係 (傾き) を調べる。
窓 (なし = 0、Hann の次数 1・2・3) も同じ追跡から作る。格子は S1 の full の格子 (41 × 31、Q1 の表) から 1 つおき・2 つおきなどの部分格子で、D の分布を比べる。
使い方: py -3.12 example/v2_3/S2_nturn.py <astrid2|kekls> [--nturn 4096] [--errors 種,強さ]     出力: 標準出力、S2_data/S2_nturn_<ring>.npz (D_N の配列)
条件: v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、px = py = z = 0 から出発、規格化 = 標本の楕円。格子 21 × 16 (336 点、S2_errors_scan.py と同じ)。
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from S2_errors_scan import RINGS, DATA  # noqa


def dvals(rec, lost, tw, N, window=1, norm='data'):
    """rec: (nturn, n, 4+)。最初の N 周から D (失われた点・nan は nan)"""
    from pysad.freqmap import _tunes_of
    n = rec.shape[1]
    out = np.full(n, np.nan)
    for i in range(n):
        if lost[i] > 0 and lost[i] <= N:
            continue
        sr = rec[:N, i, :4]
        if not np.all(np.isfinite(sr)):
            continue
        try:
            f, _a = _tunes_of(sr, tw, window, norm)
        except Exception:
            continue
        d = np.hypot(((f[2] - f[0] + 0.5) % 1.0) - 0.5, ((f[3] - f[1] + 0.5) % 1.0) - 0.5)
        out[i] = np.log10(d) if d > 0 else -16.0
    return out


def main():
    ring = sys.argv[1]
    nturn = int(sys.argv[sys.argv.index('--nturn') + 1]) if '--nturn' in sys.argv else 4096
    from pysad import compat
    from pysad.trksetup import TrackSetup, lower_priority
    compat.set_mode('v2')
    lower_priority()
    c = RINGS[ring]
    st, _ = TrackSetup.from_source(open(c['path'], encoding='utf-8').read(), c['line'], rad=False, rfsw=False, alost=0.1)
    tag = ring
    if '--errors' in sys.argv:          # --errors <種>,<強さ>: S2_errors_scan.py の 'all' と同じ誤差
        from S2_errors_scan import kind_sigmas
        from pysad import fmapres as FR
        sd, lv = sys.argv[sys.argv.index('--errors') + 1].split(',')
        st, _cnt = FR.apply_alignment_errors(st, seed=int(sd), streams='type', **kind_sigmas(ring, 'all', float(lv)))
        tag = '%s_err_s%s_k%s' % (ring, sd, lv)
    tw0 = st.twiss_at(c['loc'])
    tw = ((tw0['alfx'], tw0['betx']), (tw0['alfy'], tw0['bety']))
    xs, ys = np.linspace(*c['x'], 21), np.linspace(*c['y'], 16)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    X, Y = X.ravel(), Y.ravel()
    coords = np.zeros((len(X), 6))
    coords[:, 0], coords[:, 2] = X, Y
    t0 = time.time()
    r = st.track(c['loc'], coords, nturn, record=True)
    print('%s: 追跡 %d 点 × %d 周 %.0f s、失われた点 %d' % (ring, len(X), nturn, time.time() - t0, int((r['lost_turn'] > 0).sum())), flush=True)
    rec, lost = r['rec'], r['lost_turn']
    Ns = [n for n in (256, 512, 1024, 2048, 4096) if n <= nturn]
    D = {}
    for N in Ns:
        t0 = time.time()
        D[N] = dvals(rec, lost, tw, N)
        d = D[N][np.isfinite(D[N])]
        print('N = %4d: 残った点 %d (落ちた %d)、D の百分位 5/25/50/75/95 = %s、D > -5: %.3f、D > -4: %.3f、D の最小 %.2f (%.1f s)' % (
            N, len(d), int(np.isnan(D[N]).sum()), ' '.join('%.2f' % v for v in np.percentile(d, [5, 25, 50, 75, 95])), np.mean(d > -5), np.mean(d > -4), d.min(), time.time() - t0), flush=True)
    # 窓 (N = 1024)
    W = {}
    for w in (0, 1, 2, 3):
        W[w] = dvals(rec, lost, tw, 1024, window=w)
        d = W[w][np.isfinite(W[w])]
        print('N = 1024、窓 p = %d: D の中央値 %.2f、5%% %.2f、D > -5: %.3f、D の最小 %.2f' % (w, np.median(d), np.percentile(d, 5), np.mean(d > -5), d.min()), flush=True)
    # 点ごとの傾き: D_N と log10 N の最小二乗の傾き (N = 512〜4096 のうち、全部の N で残った点)
    use = [N for N in Ns if N >= 512]
    ok = np.all([np.isfinite(D[N]) for N in use], axis=0)
    ln = np.log10(use)
    sl = np.full(len(X), np.nan)
    for i in np.nonzero(ok)[0]:
        sl[i] = np.polyfit(ln, [D[N][i] for N in use], 1)[0]
    print('傾き d D / d log10 N (N = %s、%d 点): 百分位 5/25/50/75/95 = %s' % (use, ok.sum(), ' '.join('%.2f' % v for v in np.nanpercentile(sl, [5, 25, 50, 75, 95]))))
    # D_1024 の層ごとの傾き
    D1 = D[1024]
    for lo, hi in ((-20, -8), (-8, -6), (-6, -5), (-5, -4), (-4, 0)):
        m = ok & (D1 >= lo) & (D1 < hi)
        if m.any():
            print('  D_1024 が [%g, %g) の点 %d: 傾きの中央値 %.2f (25%%〜75%% = %.2f〜%.2f)' % (lo, hi, m.sum(), np.nanmedian(sl[m]), np.nanpercentile(sl[m], 25), np.nanpercentile(sl[m], 75)))
    os.makedirs(DATA, exist_ok=True)
    np.savez_compressed(os.path.join(DATA, 'S2_nturn_%s.npz' % tag), x=X, y=Y, lost=lost, slope=sl, **{'D%d' % N: D[N] for N in Ns}, **{'Dw%d' % w: W[w] for w in W})


if __name__ == '__main__':
    main()
