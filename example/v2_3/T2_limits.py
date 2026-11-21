# -*- coding: utf-8 -*-
"""v2.3.0 T2: 解析式との比較 (極限) と、近似の大きさの見積もり。
 A. 平坦ビームの極限 (λ2 → 0): touschek_piwinski の 1 次元の式は、Bruck・Piwinski の関数 D(ε) (記憶に基づく式: Handbook of Accelerator Physics の Touschek の節、
    D(ε) = √ε [ −3/2 e^{−ε} + ε/2 ∫_ε^∞ (ln u/u) e^{−u} du + 1/2 (3ε − ε ln ε + 2) ∫_ε^∞ (e^{−u}/u) du ]) に一致するか。 ε = x = (2βa)²/(4λ1)。
    比 J √(π α2) c² / (π D(x)/√x) が 1 なら一致。
 B. K64 の当てはめ touschek_f (tintrb.f の touckf、C(x)/x ∝ 1/τ) と厳密な関数 g(x) = π D(x)/√x の比: f(x) π x / g(x) − 1。
 C. 受容 a と σδ の比による補正の見積もり (バンチ内の δ の分布 a − δ_old を入れたときの損失率の増え。K64・elegant・touschek_piwinski は入れていない):
    δ_old ~ N(0, σδ²) と仮定し、受容を (a − δ_old, a + δ_old) に変えた損失率の平均 / 変えない損失率 − 1 (13 点の Gauss-Hermite。δ_old と q_s の相関は入れない概算)。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/T2_limits.py [astrid2 kekls]   (T2_limits.out.txt に保存した)
"""
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import exp1

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)


def D_bruck(e):
    i1 = quad(lambda u: math.log(u) / u * math.exp(-u), e, np.inf)[0]
    return math.sqrt(e) * (-1.5 * math.exp(-e) + e / 2 * i1 + 0.5 * (3 * e - e * math.log(e) + 2) * exp1(e))


def main():
    from pysad.touschek import _piw_J, touschek_f
    print('A. 平坦ビームの極限: J(c) √(π α2) c² / (π D(x)/√x)  (α2 = 1e9、α1 = x/c²、c = 1)')
    for x in (1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 3.0):
        J = _piw_J(1.0, x, 1e9, npan=60)
        print('  x = %-8g  J√(πα2)c² = %.6g   π D(x)/√x = %.6g   比 %.6f' % (x, J * math.sqrt(math.pi * 1e9), math.pi * D_bruck(x) / math.sqrt(x),
                                                                         J * math.sqrt(math.pi * 1e9) / (math.pi * D_bruck(x) / math.sqrt(x))))
    print('\nB. K64 の当てはめ touschek_f(x) と厳密な関数 π D(x)/(π x √x)… の比: f(x) π x / g(x) − 1,  g(x) = π D(x)/√x')
    worst = 0.0
    for x in np.logspace(-6, 1, 29):
        g = math.pi * D_bruck(x) / math.sqrt(x)
        r = touschek_f(x) * math.pi * x / g - 1
        worst = max(worst, abs(r))
        if abs(math.log10(x) * 2 - round(math.log10(x) * 2)) < 1e-9 or True:
            pass
    for x in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0):
        g = math.pi * D_bruck(x) / math.sqrt(x)
        print('  x = %-8g  f π x/g − 1 = %+.5f' % (x, touschek_f(x) * math.pi * x / g - 1))
    print('  x = 1e-6 … 10 (29 点) の最大の |差| = %.4f' % worst)
    names = sys.argv[1:]
    if names:
        from T2_scan import gather
        from pysad.touschek import piwinski_slice_rate
        print('\nC. a − δ_old の効き (概算)')
        for name in names:
            tt, tt2, ref = gather(name, 1e9, 0.01)
            info = tt2['rest_info']
            rest = tt2['rest']
            idx = np.arange(0, len(rest), max(1, len(rest) // 150))
            t, w = np.polynomial.hermite_e.hermegauss(13)
            w = w / w.sum()
            sd = ref['SIGE']
            for a in ([0.004, 0.006, 0.010, 0.016] if name == 'astrid2' else [0.004, 0.006, 0.010, 0.020]):
                r0 = sum(piwinski_slice_rate(rest[i], a, a, info['circ'], 1e9, 20) for i in idx)
                r1 = sum(wi * sum(piwinski_slice_rate(rest[i], max(a - sd * ti, 0.2 * a), max(a + sd * ti, 0.2 * a), info['circ'], 1e9, 20) for i in idx)
                         for ti, wi in zip(t, w))
                print('  %-8s a = %.4f (a/σδ = %5.1f)  損失率の増え %+.4f' % (name, a, a / sd, r1 / r0 - 1), flush=True)


if __name__ == '__main__':
    main()
