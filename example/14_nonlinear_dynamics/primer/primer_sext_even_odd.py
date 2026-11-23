# -*- coding: utf-8 -*-
"""入門 §5: 6 極の「1 次ではチューンシフトが出ず、2 次 (強さの 2 乗) で出る」ことを数値で見る。FODO リング (3 GeV、105.6 m、24 セル) の SF・SD の強さを λ 倍 (K2L = λ(4, −6))。
  (1) 解析: adts(lat).sext (6 極の 2 次) が λ² に比例する (λ = 0.5, 1, 2、−1 で同じ値)
  (2) 追跡 + NAFF + フィット (強さを変えない直接の方法。adtsdirect.direct_adts。Jmax = 1e-7、1024 周、order = 2) を λ = +1, −1, 0 で行い、
      偶数部 C = (f(1) + f(−1) − 2 f(0))/2 が解析の 6 極 2 次、奇数部 O = (f(1) − f(−1))/2 が 6 極 × 偏向磁石 (曲率・硬い端) の交差項 (解析の `.cross`)、
      f(0) が 6 極に依らない項 (運動学 + 四極端のフリンジ + 偏向磁石の 3 次どうし + 8 極) に当たることを見る。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/primer/primer_sext_even_odd.py   (所要 約 30 s)
"""
import numpy as np
import primer_common as C
from pysad import adtsdirect as D, resdrive as R

JMAX = 1e-7
lat, rows, tw, nu, circ = C.fodo()
print('FODO リング (3 GeV、105.6 m)、基準の強さ K2L = (4, −6) [1/m²]。単位 1/(m rad)、(xx, xy, yy)')
print('(1) 解析の 6 極 2 次 (adts().sext) の λ 依存:')
base = None
for lam in (0.5, 1.0, 2.0, -1.0):
    l2, r2, *_ = C.fodo(k2=(4.0 * lam, -6.0 * lam))
    a = R.adts(l2, rows=r2)
    v = np.array([a.sext['xx'], a.sext['xy'], a.sext['yy']])
    if lam == 1.0:
        base = v
    print('  λ = %+4.1f  sext = %9.3f %9.3f %9.3f   cross = %8.3f %8.3f %8.3f' % (lam, *v, a.cross['xx'], a.cross['xy'], a.cross['yy']))
print('  (λ = 1 の値との比 λ = 0.5: %s、2: %s、−1: %s)' % tuple(
    np.array2string(np.array([R.adts(C.fodo(k2=(4.0 * s, -6.0 * s))[0], rows=C.fodo(k2=(4.0 * s, -6.0 * s))[1]).sext[k] for k in ('xx', 'xy', 'yy')]) / base, precision=4) for s in (0.5, 2.0, -1.0)))
print('(2) 追跡 + NAFF + フィットの直接の ADTS (λ = +1, −1, 0):')
f = {}
for lam in (1.0, -1.0, 0.0):
    l2, r2, tw2, nu2, _ = C.fodo(k2=(4.0 * lam, -6.0 * lam))
    r = D.direct_adts(l2, C.P_FODO, JMAX, nturn=1024, order=2, rows=r2)
    f[lam] = np.array([r['xx'], r['xy_mean'], r['yy']])
    print('  λ = %+4.1f  追跡 f = %9.3f %9.3f %9.3f' % (lam, *f[lam]))
A = R.adts(lat, rows=rows)
C_even = (f[1.0] + f[-1.0] - 2 * f[0.0]) / 2
O_odd = (f[1.0] - f[-1.0]) / 2
sx = np.array([A.sext['xx'], A.sext['xy'], A.sext['yy']])
cr = np.array([A.cross['xx'], A.cross['xy'], A.cross['yy']])
rest = np.array([A.kin[k] + A.fringe[k] + A.bend[k] + A.octu[k] for k in ('xx', 'xy', 'yy')])
print('  偶数部 C = (f(1)+f(−1)−2f(0))/2 = %9.3f %9.3f %9.3f   解析の 6 極 2 次 = %9.3f %9.3f %9.3f   差/|解析| = %.1e %.1e %.1e' % (*C_even, *sx, *((C_even - sx) / np.abs(sx))))
print('  奇数部 O = (f(1)−f(−1))/2      = %9.3f %9.3f %9.3f   解析の交差項 .cross  = %9.3f %9.3f %9.3f   差 = %.3f %.3f %.3f' % (*O_odd, *cr, *(O_odd - cr)))
print('  f(0) (6 極なし)               = %9.3f %9.3f %9.3f   解析の 6 極に依らない項   = %9.3f %9.3f %9.3f   差 = %.3f %.3f %.3f' % (*f[0.0], *rest, *(f[0.0] - rest)))
