# -*- coding: utf-8 -*-
"""六極の組ごとの追跡の判定 (簡略化した設定): x-y の DA、x-δ の DA、周波数マップ、LMA。リング全周 (周期性なし)。
使い方 (どこからでも。1 つの組ごとに 1 プロセス。結果は results/track_<名前>.json):
  py -3.12 sextupole_tracking_check.py NAME [NAME ...] [--cand results/candidates.json] [--parts da,xdp,fmap,lma] [--quick]
候補 (K2L の組) は results/candidates.json (sextupole_optimize.py が作る)。ラティスは sxt_common の LATTICE の RING (= RINGSP)。六極の K2L は
ソース文字列の SEXT の定義を置き換えて追跡に渡す (pysad.sxtopt.replace_k2_in_source)。

簡略化した設定 (完全な再現ではない。README の表に元の設定と並べた):
  da    : RF・放射なし、δ = 0、リング先頭 (RING = RINGSP の先頭 = A 型の長直線部 (約 10 m) の中心、βx = 5.4 m、βy = 5.1 m) から px = py = 0 で出発。x-y の規格化した振幅 u (x = u cosθ √βx、y = u sinθ √βy) の
          9 本の光線 (θ = 0〜90°、y の下限 1e-5 m)、u = 0〜4e-3 (x ≈ 18 mm、y ≈ 9 mm まで) の粗い走査 24 点 + 二分 2 回 (各 6 点)。512 周。LOSSAMPL = ±50 mm (実機の真空槽ではない)。
          DA = 最初に失われた振幅の手前。面積は 4 象限に拡げた x-y の多角形 [mm²]
  xdp   : RF・放射なし。δ = −3〜+3% (1% 刻み 7 点)、y = 1e-4 m、x を ± (最大 24 mm) に走査 (δ の閉軌道から) して x-δ の DA の幅 (±の合計) [mm]、δ について平均。512 周
  fmap  : FrequencyMap (RF・放射なし)。位置 1、x = 0.3〜7 mm (14 点)、y = 0.3〜5 mm (8 点)、1024 周 (前半・後半 512 周ずつ)、NAFF の窓 Hann。D = log10 の拡散。
          指標 = 失われた点の割合、D > −3 の点 (カオス) の割合、生き残った点の D の平均
  lma   : MomentumAcceptance (RF 6.5 MV・放射あり、LOSSAMPL = ±20 mm、シンクロトロン振動 2 周期、tol = 2e-4、二分は 3 点)。リング全周に約 190 m おき 4 点 (周期性がないので全周の数点)
          指標 = 4 点の min(δ+, |δ−|) の平均と、δ+・|δ−| の平均 [%]
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import sxtopt as S

ap = argparse.ArgumentParser()
ap.add_argument('names', nargs='+')
ap.add_argument('--cand', default=os.path.join(C.OUT, 'candidates.json'))
ap.add_argument('--parts', default='da,xdp,fmap,lma')
ap.add_argument('--quick', action='store_true', help='点の数を減らす (試験・確認用)')
ap.add_argument('--nturn', type=int, default=512)
ap.add_argument('--out', default=C.OUT)
ap.add_argument('--errors', type=int, default=0, help='位置ずれの乱数の種 (0 = 誤差なし)。四極・六極に DX・DY を rms --err-um で与える (閉軌道の補正なし。閉軌道からの DA)')
ap.add_argument('--err-um', type=float, default=50.0)
ap.add_argument('--err-magnets', default='sext', help='sext (既定。六極だけ。軌道がほとんど動かず補正なしで意味がある) または quad+sext (補正なしだと軌道が 1 mm 以上ずれて全部の DA が 0 近くになる)')
args = ap.parse_args()
C.lower_priority()
CAND = json.load(open(args.cand, encoding='utf-8'))
names_fam = CAND['family_names']
src0 = C.load_source()
Q = 3 if args.quick else 1


from sxt_track import da_rays, da_xdp, fmap, lma


os.makedirs(args.out, exist_ok=True)
for name in args.names:
    c = CAND['candidates'][name]
    k = np.array(c['k'])
    src = S.replace_k2_in_source(src0, names_fam, k)
    t0 = time.time()
    out = dict(name=name, desc=c['desc'], group=c['group'], k=c['k'], cond=dict(nturn=args.nturn, quick=args.quick), time={})
    st = C.make_track_setup(src, rad=False, rfsw=False, alost=0.05)
    if args.errors:
        rng = np.random.default_rng(args.errors)
        for e in st.els:
            if type(e).__name__ in (('Sextupole',) if args.err_magnets == 'sext' else ('Quadrupole', 'Sextupole')):
                e.DX = float(args.err_um * 1e-6 * rng.standard_normal())
                e.DY = float(args.err_um * 1e-6 * rng.standard_normal())
        out['cond']['errors'] = dict(seed=args.errors, rms_um=args.err_um, magnets=args.err_magnets, correction='なし', cod0=[float(v) for v in st.refine_orbit(np.zeros(6))[:4]])
    parts = args.parts.split(',')
    for p in parts:
        t1 = time.time()
        if p == 'da':
            out['da'] = da_rays(st, 1, args.nturn, args.quick)
        elif p == 'xdp':
            out['xdp'] = da_xdp(st, 1, args.nturn, args.quick)
        elif p == 'fmap':
            out['fmap'] = fmap(st, 1, args.quick)
        elif p == 'lma':
            out['lma'] = lma(src, args.quick)
        out['time'][p] = time.time() - t1
        print('%s %s: %.0f s' % (name, p, out['time'][p]), flush=True)
    suffix = ('_err%d%s' % (args.errors, '' if args.err_um == 50.0 else 'u%g' % args.err_um) if args.errors else '') + ('_n%d' % args.nturn if args.nturn != 512 else '') + ('_quick' if args.quick else '')
    path = os.path.join(args.out, 'track_%s%s.json' % (name, suffix))
    json.dump(out, open(path, 'w', encoding='utf-8'), ensure_ascii=False)
    d = out.get('da')
    print('%s 完了 (%.0f s): DA面積 %s mm², x_DA %s mm, y_DA %s mm | xdp 平均幅 %s mm | fmap カオス割合 %s | LMA 平均 %s' % (
        name, time.time() - t0, '%.1f' % d['area_mm2'] if d else '-', '%.1f' % d['x_da_mm'] if d else '-', '%.1f' % d['y_da_mm'] if d else '-',
        '%.1f' % out['xdp']['mean_width_mm'] if 'xdp' in out else '-', '%.2f' % out['fmap']['frac_chaotic'] if 'fmap' in out else '-',
        '%.2f%%' % (100 * out['lma']['mean_min']) if 'lma' in out else '-'), flush=True)
