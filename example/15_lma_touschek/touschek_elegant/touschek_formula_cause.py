# -*- coding: utf-8 -*-
"""K64 の Touschek の式と独立な式の差の原因を、近似を 1 つずつ置き換えて切り分ける。
スライスごとの損失率を次の式で作り、リング全体の寿命を比べる (一様な受容 a):
  E0 : touschek_piwinski (厳密、x-y 面の 2 固有値 λ1 ≥ λ2)
  E1 : 平坦ビームの極限 (λ2 → 0) の厳密な関数 F_flat(x) を、λ1 だけで評価
  E2 : E1 と同じ関数を、K64 と同じ運動量の広がり ptrans² = tr(P) (3 固有値の和) で評価
  E3 : K64 の当てはめ touschek_f を ptrans² = tr(P) で評価 (K64 の式の再現。TouschekTable と比べて確かめる)
  E4 : K64 の当てはめ touschek_f を λ1 だけで評価
  F_flat(x): R = 2 · r0² c N/(γ V* √λ1) (2/√π) c⁻² g(x)、g(x) = ∫_1^∞ dκ κ⁻³ 2π(κ² − 1 − ln κ) e^{−xκ²}、x = c²/(4λ1)、c = 2βa (2 = 正負の受容)
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_cause.py astrid2|kekls [pbunch:mincoup ...]
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
R0 = 2.8179403262e-15
CL = 299792458.0


def g_flat(x, n=2000):
    """g(x) = ∫_1^∞ dκ κ⁻³ 2π (κ² − 1 − ln κ) e^{−x κ²}  (c = 1 の単位。F_flat = c⁻² g(x))。κ の対数で積分"""
    from numpy.polynomial.legendre import leggauss
    kmax = math.sqrt(45.0 / x) if x > 0 else 1e3
    kmax = max(kmax, 1.001)
    gx, gw = leggauss(24)
    edges = np.linspace(0.0, math.log(kmax), 60)
    tot = 0.0
    for a_, b_ in zip(edges[:-1], edges[1:]):
        u = 0.5 * (b_ - a_) * gx + 0.5 * (a_ + b_)
        k = np.exp(u)
        tot += np.sum(2 * np.pi * (k**2 - 1 - np.log(k)) / k**3 * np.exp(-x * k**2) * 0.5 * (b_ - a_) * gw * k)
    return tot


def run(name, pb, mc, amps):
    from touschek_formula_scan import gather
    from pysad.touschek import touschek_lifetime, touschek_piwinski, touschek_f
    tt, tt2, ref = gather(name, pb, mc)
    tab = tt['table']
    info = tt2['rest_info']
    circ, h0, p0 = info['circ'], info['h0'], info['p0']
    beta = p0 / h0
    f0 = CL / circ
    print('\nPBUNCH=%g MINCOUP=%g  εy/εx=%.4f' % (pb, mc, ref['EMITYC'] / ref['EMITXC']))
    print('   a     K64(table)  E3(replica)   E0(独立)  E1(flat,λ1)  E2(flat,tr)   E4(K64fit,λ1)  [寿命 s]   E3/E0-1  E1/E0-1  E2/E0-1  E4/E0-1')
    for a in amps:
        c = 2 * beta * a
        r = dict(E1=0.0, E2=0.0, E3=0.0, E4=0.0)
        for sl in tt2['rest']:
            vol = (4 * math.pi)**1.5 * math.sqrt(abs(np.linalg.det(sl['xx'])))
            P = np.asarray(sl['pl']) * p0**2
            l1 = np.sort(np.linalg.eigvalsh(P[:2, :2]))[1]
            tr = np.trace(P)
            w = sl['al'] / circ
            for key, lam, fn in (('E1', l1, 'flat'), ('E2', tr, 'flat'), ('E3', tr, 'k64'), ('E4', l1, 'k64')):
                x = c**2 / (4 * lam) if fn == 'flat' else (beta * a)**2 / lam
                if fn == 'flat':
                    # R = r0² c N/(γ V* √λ) (2/√π) c⁻² g(x)  [per ring fraction w]
                    r[key] += w * R0**2 * CL * pb / (h0 * vol * math.sqrt(lam)) * 2.0 * (2 / math.sqrt(math.pi)) * g_flat(x) / c**2          # 2.0 = 正負の受容の 2 項
                else:
                    # K64: r0² N √π f0 /p0 · al/vol · ptrans_mc⁻³ · f(x)  (al/C を w に)
                    r[key] += w * R0**2 * pb * math.sqrt(math.pi) * CL / p0 * (1.0 / vol) * lam**-1.5 * touschek_f(x) * circ / circ * (1.0)
        k64 = touschek_lifetime(tab, None, None, a / ref['SIGE'], ref['SIGE'])
        e0 = touschek_piwinski(tt2, a=a)['lifetime']
        L = {k: 1.0 / v for k, v in r.items()}
        print('  %.4f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f      %+.4f  %+.4f  %+.4f  %+.4f' % (
            a, k64, L['E3'], e0, L['E1'], L['E2'], L['E4'], L['E3'] / e0 - 1, L['E1'] / e0 - 1, L['E2'] / e0 - 1, L['E4'] / e0 - 1), flush=True)


if __name__ == '__main__':
    name = sys.argv[1]
    sets = [tuple(float(v) for v in s.split(':')) for s in sys.argv[2:]] or [(1e9, 0.01)]
    amps = [0.004, 0.006, 0.010, 0.016] if name == 'astrid2' else [0.004, 0.006, 0.010, 0.020]
    for pb, mc in sets:
        run(name, pb, mc, amps)
