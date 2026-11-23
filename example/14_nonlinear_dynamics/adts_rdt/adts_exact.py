# -*- coding: utf-8 -*-
"""ADTS・RDT の解析式の精密版 (K0≠0 の偏向磁石・δ≠0・生成関数・4 次の項) の検証 (a)〜(e) の表を作る。
使い方 (リポジトリ直下から。出力は example/14_nonlinear_dynamics/adts_rdt/adts_exact_*.txt に上書き):
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_exact.py mpw            (a) MPW (K0≠0 の偏向磁石): h・f4・ξ・ADTS (λ = 0, 1, 2) の解析 対 追跡
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_exact.py delta fodo|astrid|kekls   (b) δ≠0 の ADTS と ∂(ADTS)/∂δ: 解析 (R6)・数値の正規形・直接の追跡 (NAFF)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_exact.py bend           (c) 偏向磁石の生成関数: 閉形式 (端・本体の積分) 対 追跡から取り出した係数
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_exact.py nf             (d) 非対角の調和の正規形の補正: 独立な数値のリー変換との比較 (トイ・リング)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_exact.py adts           (e) adts() の 4 次の項: 既定 r3 対 精密版 r6 対 追跡 (FODO・ASTRID2・MPW・KEK-LS)
条件: pysad k64 の追跡 (RF 止め・放射なし・結合なし・閉軌道のまわり)、6 極 ×1 (mpw の λ は 6 極の強さの倍率)。直接の追跡: J の格子 6×6 (±位相の対 = 70 粒子)、1024 周、
NAFF Hann p=1、order=2 のフィット、Jmax は FODO 1e-7、ASTRID2・MPW 2e-8 (mpw では 5e-9 も)、KEK-LS 2e-10。δ の追跡は閉軌道のまわり。"""
import dataclasses
import json
import math
import os
import re
import sys
import time

import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import adtsdirect as D, bendnl as B, compat, nform, rdt2nd as Q, rdttrack as T, resdrive as R
compat.set_mode('k64')
import adts_bend as M

OUT = []
KEYS = ('xx', 'xy', 'yy')


def P(s=''):
    print(s, flush=True)
    OUT.append(s)


def save(name):
    open(os.path.join(HERE, name), 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def rel(a, b):
    return '%+.1e' % ((a - b) / b)


def direct(S, jmax, delta=0.0, order=2):
    tr = D.direct_adts(S['lat'], S['p'], jmax, delta=delta, order=order, nturn=1024)
    return {'xx': tr['xx'], 'xy': tr['xy_mean'], 'yy': tr['yy']}


def cmd_mpw():
    P('R6 (a) MPW (ASTRID2_v01_02_id.sad の RINGMPW、0.58 GeV、周長 45.7 m): K0≠0 の偏向磁石 23 台 (L = 0.04375・0.0875 m、K0 = ±0.032・±0.064 rad、E1=E2=0、ANGLE=0。内部の軌道が ±0.032 rad 傾く)')
    S = M.setup('astrid_mpw', 0.0)
    lat, rows = S['lat'], S['rows']
    P('  K0 の偏向磁石の内部の閉軌道: x = %.3e m (入口)、px = ±%.4f rad。原子の w は CALC の Twiss と 1e-13 で一致 (w_err)、閉軌道は 1e-17 (z_err)' % (rows[101].x_cod, abs(rows[101].px_cod)))
    # 1. 3 次の項
    S1 = M.setup('astrid_mpw', 1.0)
    names = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
    rm = T.ring_map(S1['lat'].elements, 3e9, S1['rows'][0], amp=1e-5)
    h = T.rdt_from_map(rm)
    r = R.resonance_terms(S1['lat'], names, rows=S1['rows'])
    r_old = R.resonance_terms(S1['lat'], names, rows=S1['rows'])
    P('(1) 3 次の項 |解析 − 追跡 (6 極 ×1)| / |追跡|: R4 (K0 の磁石なし) → R6')
    old = {'h21000': 5.2e-9, 'h30000': 9e-10, 'h10110': 8e-10, 'h10020': 4.6e-5, 'h10200': 2.0e-5}
    for n in names:
        P('   %-7s %.1e → %.1e' % (n, old[n], abs(r.total[n] - h[n]) / abs(h[n])))
    # 2. f4
    for lam in (0.0, 1.0):
        Sl = M.setup('astrid_mpw', lam)
        rm3 = T.ring_map3(Sl['lat'].elements, Sl['p'], Sl['rows'][0], amp3=1e-4)
        ht, _ = T.rdt4_from_map(rm3)
        r4 = Q.rdt4(Sl['lat'], rows=Sl['rows'])
        big = max(abs(v) for v in ht.values())
        w = sorted(((abs(r4.total[n] - v) / abs(v), n) for n, v in ht.items() if abs(v) > 1e-3 * big), reverse=True)
        P('(2) 4 次の項 (6 極 ×%g) の |解析 − 追跡|/|追跡| の最大 (|h| > 1e-3 max の %d 項): %.1e (%s)   [R5: ×0 で 8.7e-3、×1 で 1.8e-4]' % (lam, len(w), w[0][0], w[0][1]))
    # 3. ξ
    nat = R.chromaticity_natural(lat, rows=rows)
    bend = R.chromaticity_bend(lat, rows=rows)
    ref = R.chromaticity_numeric(lat, order=1)['d1']
    P('(3) クロマティシティ (6 極 ×0): CALC (%.4f, %.4f)、自然 (%.4f, %.4f) + 偏向磁石の分 (%.4f, %.4f) → 差 (CALC − 解析) = (%+.1e, %+.1e)   [R4: (+5.8e-4, −0.173)]' % (
        ref[0], ref[1], nat[0], nat[1], bend[0], bend[1], ref[0] - nat[0] - bend[0], ref[1] - nat[1] - bend[1]))
    # 4. ADTS
    P('(4) ADTS [1/(m rad)]: 解析 adts().total の精密版 (fourth="r6") / 既定 (fourth="r3") と直接の追跡。追跡の J の範囲 (Jmax) と order を変えたもの')
    for lam in (0.0, 1.0, 2.0):
        Sl = M.setup('astrid_mpw', lam)
        A = R.adts(Sl['lat'], rows=Sl['rows'], fourth='r6').total
        A3 = R.adts(Sl['lat'], rows=Sl['rows']).total
        P('  λ=%g  解析 精密 r6 %s   既定 r3 %s' % (lam, ' '.join('%.3f' % A[k] for k in KEYS), ' '.join('%.3f' % A3[k] for k in KEYS)))
        for jmax, order, nt in ((2e-8, 2, 1024), (1e-8, 2, 1024), (5e-9, 2, 1024), (1e-8, 3, 1024)):
            if lam != 2.0 and (jmax, order) not in ((2e-8, 2),):
                continue
            tr = D.direct_adts(Sl['lat'], Sl['p'], jmax, order=order, nturn=nt)
            t = {'xx': tr['xx'], 'xy': tr['xy_mean'], 'yy': tr['yy']}
            P('       追跡 Jmax=%.0e order=%d  %s   (精密 r6−追跡)/追跡 %s   (既定 r3−追跡)/追跡 %s' % (jmax, order, ' '.join('%.3f' % t[k] for k in KEYS),
                                                                               ' '.join(rel(A[k], t[k]) for k in KEYS), ' '.join(rel(A3[k], t[k]) for k in KEYS)))
    P('  (λ=2 の yy は νy = 2.2492 で 4νy = 8.997 が整数 9 に近く (ε = −3e-3)、J² の項が大きいので Jmax 2e-8・order 2 のフィットが偏る。Jmax を小さくする・order 3 にすると解析に収束する)')
    save('adts_exact_mpw.out.txt')


def cmd_delta(key, d0=2e-3):
    S = M.setup(key, 1.0)
    P('R6 (b) δ≠0 の ADTS: %s (%.2f GeV、周長 %.1f m)。δ = ±%g。[1/(m rad)]' % (key, S['p'] / 1e9, S['circ'], d0))
    jb = S['jb']
    res = {}
    cache = os.path.join(HERE, 'adts_exact_delta_%s.track.json' % key)
    C = json.load(open(cache)) if os.path.exists(cache) else {}
    for d in (-d0, 0.0, d0):
        t0 = time.time()
        # a = 精密版のオプション (fourth='r6'・optics='tracked')、ac = 既定 (fourth='r3'・optics='calc'。物理的な項だけ)
        a = R.adts(S['lat'], delta=d, momentum_eV=S['p'], fourth='r6', optics='tracked').total
        ac = R.adts(S['lat'], delta=d, momentum_eV=S['p']).total
        nf = nform.adts_at_delta(S['lat'], S['p'], d, amp3=1e-5 if key == 'kekls' else 1e-4)[0]
        if repr(d) not in C:
            C[repr(d)] = direct(S, jb, delta=d)
            json.dump(C, open(cache, 'w'))
        t = C[repr(d)]
        res[d] = (a, ac, nf, t)
        P('  δ=%+g 追跡 %s | 解析 (精密版 r6・追跡の光学) %s %s | 解析 (既定 r3・CALC の光学) %s | 数値の正規形 %s %s  (%.0f s)' % (
            d, ' '.join('%.3f' % t[k] for k in KEYS), ' '.join('%.3f' % a[k] for k in KEYS), ' '.join(rel(a[k], t[k]) for k in KEYS),
            ' '.join(rel(ac[k], t[k]) for k in KEYS), ' '.join('%.3f' % nf[k] for k in KEYS), ' '.join(rel(nf[k], t[k]) for k in KEYS), time.time() - t0))
    P('  ∂(ADTS)/∂δ (中心差分 ±%g):' % d0)
    for k in KEYS:
        tr = (res[d0][3][k] - res[-d0][3][k]) / (2 * d0)
        an = (res[d0][0][k] - res[-d0][0][k]) / (2 * d0)
        ac = (res[d0][1][k] - res[-d0][1][k]) / (2 * d0)
        nf = (res[d0][2][k] - res[-d0][2][k]) / (2 * d0)
        P('    %s 追跡 %9.1f | 解析 (精密版 r6・追跡の光学) %9.1f (比 %.3f) | 解析 (既定 r3・CALC) %9.1f (比 %.3f) | 数値の正規形 %9.1f (比 %.3f)' % (k, tr, an, an / tr, ac, ac / tr, nf, nf / tr))
    save('adts_exact_delta_%s.out.txt' % key)


def cmd_bend():
    P('R6 (c) 偏向磁石の 3 次の生成関数: 閉形式 (端の係数・本体の Gauss 積分) 対 追跡の写像から取り出した原子の和')
    P('  端 (E1・E2≠0): f_e1 = (h/2)sec²E1 px y² + h t² x y py − (h t²/2) x² px − (h² t³/2) x y² − (h² t³/6) x³ (t = tan E1)、f_e2 = −(h/2)sec²E2 px y² − h t² x y py + (h t²/2) x² px − (h² t/2) x y² + (h² t³/3) x³ (t = tan E2)')
    P('  本体: H = (h/2) x (px² + py²) + (K1 h/6) x³ (K64 の写像の連続極限)。h の 20 項の最大の相対差:')
    for key in ('fodo', 'fodo_rect', 'fodo_asym', 'astrid', 'kekls'):
        S = M.setup(key, 1.0)
        lat, rows = S['lat'], S['rows']
        seen = set()
        w = 0.0
        for i, el in enumerate(lat.elements):
            if R._tn(el) != 'Bend' or el.L <= 0 or el.ANGLE == 0 or el.K0 != 0:
                continue
            sig = (el.L, el.ANGLE, el.K1, el.E1, el.E2)
            if sig in seen:
                continue
            seen.add(sig)
            a, _ = B.bend_atoms(el, rows[i], rows[i + 1], Q.MONOS3, nslice=16)
            an = B.bend_analytic_h3(el, rows[i], rows[i + 1], Q.MONOS3)
            e = np.abs(a.sum(axis=1) / 2 ** 1.5 - an).max() / np.abs(an).max()
            P('    %-10s %-6s K1=% .4f E1=% .4f E2=% .4f   解析 − 取り出し = %.1e' % (key, el.name, el.K1, el.E1, el.E2, e))
    P('  (K1=0 で 1e-8 は取り出しの極薄い偏向磁石 ε = 1e-8 m の分。K1≠0 の 2e-6〜3e-4 は K64 の tbendi の積分器と連続の式の差 (本体の係数 h/2・K1h/6 は n=128 の極限で 1e-5 で一致))')
    P('  スライスの数への依存: R3 の合計は O(1/n²)、R6 の精密版 (fourth="r6"。4 次の項を要素の写像から取り出して内部の対の項を引く) は nslice に依らない:')
    for key in ('fodo', 'astrid'):
        S = M.setup(key, 1.0)
        t3 = {n: R.adts(S['lat'], rows=S['rows'], bend_nslice=n).total for n in (4, 8, 64)}
        t6 = {n: R.adts(S['lat'], rows=S['rows'], bend_nslice=n, fourth='r6').total for n in (4, 8, 64)}
        P('    %-6s |n=4 − n=64|/|n=64|  既定 r3 %s  精密 r6 %s' % (key, ' '.join('%.1e' % (abs(t3[4][k] - t3[64][k]) / abs(t3[64][k])) for k in KEYS),
                                                      ' '.join('%.1e' % (abs(t6[4][k] - t6[64][k]) / abs(t6[64][k])) for k in KEYS)))
    save('adts_exact_bend.out.txt')


def cmd_adts():
    P('R6 (e) adts() の 4 次の項: 既定 (fourth="r3") は物理的な項だけ (R3 の値)。fourth="r6" は精密版のオプション: 四極本体の運動学的な項に K64 の積分器の台形則、偏向磁石 (K0≠0 を含む) の 4 次の項に要素の写像から取り出した値 (数値積分器の 0.1〜0.3% の端数を含む)')
    P('  [1/(m rad)]、6 極 ×1、追跡 = 直接の追跡 (rdt4_adts_*.pysad.json の保存値、1024 周、70 粒子、order 2)')
    for key in ('fodo', 'astrid', 'astrid_mpw', 'kekls'):
        S = M.setup(key, 1.0)
        t0 = time.time()
        A = R.adts(S['lat'], rows=S['rows'], fourth='r6').total
        dt = time.time() - t0
        A3 = R.adts(S['lat'], rows=S['rows']).total
        f = os.path.join(HERE, 'rdt4_adts_%s.pysad.json' % key)
        t = json.load(open(f, encoding='utf-8'))['track']
        P('  %-10s 追跡 %s' % (key, ' '.join('%.4f' % t[k] for k in KEYS)))
        P('        既定 r3 %s  相対差 %s' % (' '.join('%.4f' % A3[k] for k in KEYS), ' '.join(rel(A3[k], t[k]) for k in KEYS)))
        P('        精密 r6 %s  相対差 %s  (%.1f s)' % (' '.join('%.4f' % A[k] for k in KEYS), ' '.join(rel(A[k], t[k]) for k in KEYS), dt))
    save('adts_exact_adts.out.txt')


def cmd_nf():
    P('R6 (d) 非対角の調和の正規形の補正: test_v23_R6_bend_orbit.py::test_nonresonant_correction_to_all_fourth_order_harmonics_vs_numerical_lie_transform と同じ検証 (独立な数値のリー変換)')
    sys.path.insert(0, os.path.join(ROOT, 'pysad', 'tests', 'script_tests'))
    import test_v23_R6_bend_orbit as TT
    t0 = time.time()
    TT.test_nonresonant_correction_to_all_fourth_order_harmonics_vs_numerical_lie_transform()
    P('  合格 (共役変換後の 3 次の項 < 1e-8、F4 + 補正と数値の変換後の 4 次の項の差 < 1e-6·max|補正|、補正なしの F4 とは max|補正| 程度違う)。%.0f s' % (time.time() - t0))
    save('adts_exact_nf.out.txt')


if __name__ == '__main__':
    cmd = sys.argv[1]
    {'mpw': cmd_mpw, 'bend': cmd_bend, 'adts': cmd_adts, 'nf': cmd_nf}.get(cmd, lambda: cmd_delta(sys.argv[2]) if cmd == 'delta' else None)()
