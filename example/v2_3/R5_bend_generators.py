# -*- coding: utf-8 -*-
"""v2.3.0 R5: 偏向磁石の 3 次の生成関数 (bendnl.bend_decomposition が K64 の追跡の写像の 2 次のテイラー展開から取り出す) と解析式の比較 (docs/49 §8・§10.3)。
  (式 B1) スライス (長さ ds、E=0): f_b = ds [ (h/2) x (px² + py²) + (K1 h/6) x³ ]       (h = ANGLE/L)
  (式 B2) 硬い端 (Maxwell、E=0): f_e1 = +(h/2) px y²、f_e2 = −(h/2) px y²
  (式 B3) E1・E2 ≠ 0 の端: 解析式は入れていない (取り出した係数だけを示す。楔の幾何 + Maxwell の端の 2 次のテイラー展開)
使い方: py -3.12 example/v2_3/R5_bend_generators.py   (出力 R5_bend_generators.out.txt)
条件: pysad k64 の追跡、E=0 の端は 1 要素の写像、取り出しの有限差分は振幅 5e-3 (extract_map の既定)・hlin=1e-4、スライス数 16 (ds = L/16)、運動量 3 GeV (4 次元・δ=0 の写像は運動量によらない)。
四極のフリンジと運動学的な項の比較: 四極 1 個 (L=0.3 m、K1=2 1/m²) の写像の 4 次の生成関数 (3 次のテイラー展開から) と解析 (rdt2nd.kinematic_h4・fringe_h4)。"""
import dataclasses
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import bendnl as B, rdttrack as T, rdt2nd as Q, compat
from pysad.elements import Quadrupole
compat.set_mode('k64')
import R3_adts_bend as M

OUT = []


def P(s=''):
    print(s)
    OUT.append(s)


def first_bend(key):
    S = M.setup(key)
    return next(e for e in S['lat'].elements if type(e).__name__ == 'Bend' and e.L > 0 and e.ANGLE != 0), S


def show(label, el, nslice=16):
    d = B.bend_decomposition(el, nslice, 3e9)
    L, ang = float(el.L), float(el.ANGLE)
    h = ang / L
    ds = L / nslice
    K1 = float(el.K1)
    P('%s: L = %.4f m、ANGLE = %.5f rad (h = %.5f 1/m)、K1 = %.4f 1/m²、E1 = %g、E2 = %g、スライス %d (ds = %.5f m)' % (label, L, ang, h, K1, float(el.E1), float(el.E2), nslice, ds))
    fb = d['f_b']
    ana = {(1, 2, 0, 0): 0.5 * h * ds, (1, 0, 0, 2): 0.5 * h * ds, (3, 0, 0, 0): K1 * h / 6 * ds}
    P('  (B1) スライスの生成関数 f_b の係数 [x px²、x py²、x³]: 取り出し | 解析 | 相対差')
    for mon, nm in (((1, 2, 0, 0), 'x px^2'), ((1, 0, 0, 2), 'x py^2'), ((3, 0, 0, 0), 'x^3')):
        v, a = fb.get(mon, 0.0), ana[mon]
        P('     %-7s %+.6e | %+.6e | %s' % (nm, v, a, ('%.1e' % (abs(v - a) / abs(a))) if a != 0 else '(解析 0、取り出し %.1e)' % v))
    others = {k: v for k, v in fb.items() if k not in ana}
    big = max([abs(v) for v in others.values()] + [0.0])
    P('     そのほかの係数の最大 |c| = %.1e (解析 0、(h/2) ds = %.3e に対し %.1e)' % (big, 0.5 * h * ds, big / (0.5 * h * ds)))
    for nm, key, sg in (('入口の端 f_e1', 'f_e1', +1), ('出口の端 f_e2', 'f_e2', -1)):
        f = d[key]
        a = sg * 0.5 * h
        v = f.get((0, 1, 2, 0), 0.0)
        oth = max([abs(c) for k, c in f.items() if k != (0, 1, 2, 0)] + [0.0])
        if float(el.E1 if key == 'f_e1' else el.E2) == 0.0:
            P('  (B2) %s: px y² の係数 取り出し %+.6e | 解析 %+.6e | 相対差 %.1e、そのほかの最大 %.1e' % (nm, v, a, abs(v - a) / abs(a), oth))
        else:
            P('  (B3) %s (E ≠ 0): 取り出した係数 %s (解析式なし)' % (nm, ', '.join('%s: %+.3e' % (''.join(map(str, k)), c) for k, c in sorted(f.items(), key=lambda kv: -abs(kv[1]))[:5])))
    P('')


def main():
    P('R5 偏向磁石の 3 次の生成関数: 取り出し (K64 の追跡の写像) 対 解析式')
    P('')
    el, S = first_bend('fodo_rect')
    el0 = dataclasses.replace(el, E1=0.0, E2=0.0)
    show('FODO の偏向磁石 (セクター型、E=0)', el0)
    el, S = first_bend('astrid')
    show('ASTRID2 の複合機能の偏向磁石 (K1≠0。E は設計の値)', el)
    el, S = first_bend('fodo_asym')
    show('FODO の偏向磁石 (E1 = 0.8·ANGLE、E2 = −0.3·ANGLE)', el)
    # 四極 1 個
    P('四極 1 個の 4 次の生成関数: 追跡の写像の 3 次のテイラー展開 対 解析 (運動学的な項の台形則 + 端のフリンジ)  (L=0.3 m、K1=2、3 GeV、amp3=2e-4)')
    q = Quadrupole(name='q', L=0.3, K1=2.0)
    row0 = SimpleNamespace(betx=1.0, alfx=0.0, mux=0.0, bety=1.0, alfy=0.0, muy=0.0, etx=0.0, etpx=0.0, s=0.0)
    Mx = np.asarray(q.transfer_matrix(0.0))[:4, :4]
    row1 = SimpleNamespace(betx=Mx[0, 0] ** 2 + Mx[0, 1] ** 2, alfx=-Mx[0, 0] * Mx[1, 0] - Mx[0, 1] * Mx[1, 1], mux=math.atan2(Mx[0, 1], Mx[0, 0]),
                           bety=Mx[2, 2] ** 2 + Mx[2, 3] ** 2, alfy=-Mx[2, 2] * Mx[3, 2] - Mx[2, 3] * Mx[3, 3], muy=math.atan2(Mx[2, 3], Mx[2, 2]), etx=0.0, etpx=0.0, s=0.0)
    for dis, lab in ((1.0, 'DISFRIN=1 (フリンジ切る: 運動学的な項だけ)'), (0.0, 'DISFRIN=0 (K64 の既定: フリンジあり)')):
        q.DISFRIN = dis
        rm = T.ring_map3([q], 3e9, row0, amp3=2e-4, orbit=np.zeros(4))
        trk, _ = T.rdt4_from_map(rm)
        kin, fr = Q.kinematic_h4([q], [row0, row1]), Q.fringe_h4([q], [row0, row1])
        kin_ex = Q.kinematic_h4([q], [row0, row1], quad_rule='exact')
        P('  %s' % lab)
        for nm in ('h22000', 'h11110', 'h00220', 'h40000', 'h20020', 'h31000'):
            i = Q._IDX4[tuple(int(c) for c in nm[1:5])]
            tot = kin[i] + fr[i]
            P('     %-7s 追跡 %+.6e%+.6ej | 解析 運動学 %+.6e (厳密な積分 %+.6e) + フリンジ %+.6e%+.6ej = %+.6e | 相対差 %.1e' % (
                nm, trk[nm].real, trk[nm].imag, kin[i].real, kin_ex[i].real, fr[i].real, fr[i].imag, tot.real, abs(tot - trk[nm]) / abs(trk[nm])))
    open(os.path.join(HERE, 'R5_bend_generators.out.txt'), 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
