# -*- coding: utf-8 -*-
"""v2.3.0 U2: K1 の定義のオプション ('log'|'lin'、B160) と B161 (CALC の dk の修正) の検証の表を作る。すべて v2 モード。
使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/v2_3/U2_verify.py   (標準出力を U2_verify.out.txt に保存した)
 1. n = 0.4 の扇形磁石 16 台のリング (3 GeV 相当): dν/dδ (CALC・追跡・独立な積分・解析式・elegant (T1 の保存値))。追跡のチューンと積分の差 (5 つの δ)
 2. ASTRID2 の自然な色収差 (6 極なし) を、独立な積分 (ドリフト・四極・複合機能の扇形磁石 + 端の縦方向の標準の集束。pysad・SAD のコードを使わない) と比べる
 3. B161: 同じリングの CALC のチューンの δ 依存 (dk 修正の有無) と ASTRID2・大きなリング A の ξ・∂²ν/∂δ²・Emittance[] への影響
独立な積分は pysad/tests/_u1_indep.py・_u2_indep.py。"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'pysad', 'tests'))
sys.path.insert(0, HERE)
import _u1_indep as U1  # noqa: E402
import _u2_indep as U2  # noqa: E402

from pysad import compat, resdrive  # noqa: E402
from pysad.trksetup import TrackSetup  # noqa: E402

N, LB, NIDX = 16, 1.0, 0.4
ANG = 2 * math.pi / N
RHO = LB / ANG
K1 = -NIDX / RHO**2
H = 1.0 / RHO
nux, nuy = math.sqrt(1 - NIDX), math.sqrt(NIDX)


def frac(v):
    return min(v % 1, 1 - v % 1)


def ring_src(extra=''):
    return ('MOMENTUM = 3 GEV;\nMARK M0 = ();\nBEND B = (L = %g ANGLE = %.12f K1 = %.12f*%g %s);\nLINE RING = (M0 %d*B);\n' % (LB, ANG, K1, LB, extra, N))


def part1():
    print('## 1. n = 0.4 の扇形磁石 16 台 (3 GeV 相当、K1 = %.6f /m^2、ρ0 = %.4f m、h = %.4f /m)' % (K1, RHO, H))
    print('解析式 (B\' 一様): dνx/dδ = %.5f、dνy/dδ = %.5f   [ν = (%.6f, %.6f)]' % (
        -NIDX * (1 + NIDX) / (2 * nux * (1 - NIDX)), nuy * (1 + NIDX) / (2 * (1 - NIDX)), nux, nuy))
    print('elegant CSBEND (T1 の保存値 T1_bend_k1_elegant.out.txt): K2 = 0 → (-0.60246, +0.73786)、K2 = -hK1 → (-0.17213, +0.21082)')
    print('\n%-6s %-24s %-24s %-24s' % ('K1 定義', 'CALC (解析行列)', '独立な積分', '比 (CALC/積分)'))
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        _, ffs = TrackSetup.from_source(ring_src(), 'RING', rad=False, rfsw=False, alost=1.0)
        r = resdrive.chromaticity_numeric(ffs.lat)['d1']
        cx, cy, _, _ = U1.dnu_ddelta(d, RHO, K1, N)
        print('%-6s (%9.5f, %9.5f)   (%9.5f, %9.5f)   (%.5f, %.5f)' % (d, r[0], r[1], cx, cy, r[0] / cx, r[1] / cy))
    print('\n追跡 (閉軌道の周りの 1 周の写像の線形化、DISFRIN=1) と CALC のチューン − 独立な積分のチューン (端数):')
    print('%-6s %8s %12s %12s %12s %12s' % ('K1 定義', 'δ', 'track x', 'track y', 'CALC x', 'CALC y'))
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st, ffs = TrackSetup.from_source(ring_src('DISFRIN=1'), 'RING', rad=False, rfsw=False, alost=1.0)
        for delta in (-0.004, -0.002, 0.0, 0.002, 0.004):
            p = 1 + delta
            xc = U1.closed_orbit_radius(U1.make_field(d, RHO, K1), RHO, p) - RHO
            s0 = np.array([xc, 0, 0, 0, 0, delta], float)

            def f(c):
                cc = np.array(c, float)
                cc[1] *= p
                cc[3] *= p
                r = st.track(1, cc.reshape(1, 6), 1, record=True)['rec'][0, 0]
                return np.array([r[0], r[1] / p, r[2], r[3] / p])
            m = np.zeros((4, 4))
            for j in range(4):
                dd = np.zeros(6)
                dd[j] = 1e-7
                m[:, j] = (f(s0 + dd) - f(s0 - dd)) / 2e-7
            tx = math.acos(0.5 * (m[0, 0] + m[1, 1])) / 2 / math.pi
            ty = math.acos(0.5 * (m[2, 2] + m[3, 3])) / 2 / math.pi
            ix, iy, _ = U1.tunes(d, RHO, K1, N, delta)
            t = ffs.lat.twiss_table(delta=delta)[-1]
            print('%-6s %8.3f %12.2e %12.2e %12.2e %12.2e' % (d, delta, tx - frac(ix), ty - frac(iy),
                                                              frac(t.mux / 2 / math.pi) - frac(ix), frac(t.muy / 2 / math.pi) - frac(iy)))
    print('\n追跡の 1 周 (16 台) の写像 − 独立な積分の最大差 (DISFRIN=1、x と y の振幅):')
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st = TrackSetup.from_source(ring_src('DISFRIN=1'), 'RING', rad=False, rfsw=False, alost=1.0)[0]
        f = U1.make_field(d, RHO, K1)
        for delta in (0.0, 0.004, -0.004):
            p = 1 + delta
            xc = U1.closed_orbit_radius(f, RHO, p) - RHO
            for lab, s0 in (('x 1 mm, y 0.15 mm', [xc + 1e-3, 0.5e-3, 1.5e-4, -0.3e-3]), ('x 5 mm, y 7.5 mm', [xc + 5e-3, 2.5e-3, 7.5e-3, -1.5e-3])):
                s0 = np.array(s0)
                c = np.array([s0[0], s0[1] * p, s0[2], s0[3] * p, 0.0, delta])
                r = st.track(1, c.reshape(1, 6), 1, record=True)['rec'][0, 0]
                s = s0.copy()
                for _ in range(N):
                    s = U1.sector_map(f, RHO, ANG, p, s)
                print('%-6s δ=%+.3f %-18s %.2e' % (d, delta, lab, np.max(np.abs(s - np.array([r[0], r[1] / p, r[2], r[3] / p])))))
    compat.set_bend_k1_def(None)


def astrid_src(nat=True):
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    if nat:
        src = src.replace('K2 =-11.729486011842393', 'K2 = 0').replace('K2 =8.681418961367605', 'K2 = 0')
    return src


def part2():
    print('\n## 2. ASTRID2 の自然な色収差 (6 極なし。DISFRIN=1 で四極・偏向磁石の端のフリンジなし) を独立な積分と比べる')
    src = astrid_src()
    nof = (src.replace('QF3 = ( L =0.1800', 'QF3 = ( DISFRIN=1 L =0.1800').replace('QF2 = ( L =0.1800', 'QF2 = ( DISFRIN=1 L =0.1800')
           .replace('BH = (L = 0.8500', 'BH = (DISFRIN=1 L = 0.8500'))
    st, ffs = TrackSetup.from_source(nof, 'RING', rad=False, rfsw=False, alost=0.01)
    ring = []
    for e in st.els:
        t = type(e).__name__
        L = float(e.L)
        if t in ('Drift', 'Sextupole'):
            ring.append(('D', L))
        elif t == 'Quadrupole':
            ring.append(('Q', L, float(e.K1)))
        elif t == 'Bend' and L != 0:
            ring.append(('B', '%s', L / float(e.ANGLE), float(e.ANGLE), float(e.K1)))
    DM = 3e-3
    ds = DM * np.cos(np.pi * (np.arange(9) + 0.5) / 9)
    print('%-6s %-44s %-44s' % ('K1 定義', 'pysad CALC: ξ (x, y) | ∂²ν/∂δ² (x, y)', '独立な積分 (端の縦集束あり)'))
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        r = resdrive.chromaticity_numeric(ffs.lat, dmax=DM)
        rr = [x if x[0] != 'B' else ('B', d) + x[2:] for x in ring]
        fl = {}
        t = np.array([U2.tunes(rr, x, fields=fl, edge=True) for x in ds])
        co = [np.polynomial.polynomial.polyfit(ds, t[:, i], 6) for i in range(2)]
        sg = [np.sign(r['d1'][i]) * np.sign(co[i][1]) for i in range(2)]      # 端数 (整数 − ν など) の符号を pysad の ξ の符号にそろえる
        ind1 = [sg[i] * co[i][1] for i in range(2)]
        ind2 = [sg[i] * 2 * co[i][2] for i in range(2)]
        print('%-6s (%8.4f, %8.4f) | (%7.2f, %7.2f)   (%8.4f, %8.4f) | (%7.2f, %7.2f)' % (d, r['d1'][0], r['d1'][1], r['d2'][0], r['d2'][1], ind1[0], ind1[1], ind2[0], ind2[1]))
    compat.set_bend_k1_def(None)


def part3():
    print('\n## 3. B161: CALC の dk の修正 (compat.set_bend_dk_fix) の影響')
    compat.set_mode('v2')
    st, ffs = TrackSetup.from_source(ring_src(), 'RING', rad=False, rfsw=False, alost=1.0)
    print('n = 0.4 のリング: CALC の ν(δ) − 独立な積分 (K1 の定義 log)')
    print('%8s %14s %14s %14s %14s' % ('δ', 'K64 の dk x', 'K64 の dk y', '修正 x', '修正 y'))
    for d in (0.0, 0.01, 0.02, 0.04, 0.06):
        ix, iy, _ = U1.tunes('log', RHO, K1, N, d)
        row = []
        for fix in (False, True):
            compat.set_bend_dk_fix(fix)
            t = ffs.lat.twiss_table(delta=d)[-1]
            row += [t.mux / 2 / math.pi - ix, t.muy / 2 / math.pi - iy]
        print('%8.3f %14.2e %14.2e %14.2e %14.2e' % (d, row[0], row[1], row[2], row[3]))
    for name, path in (('ASTRID2', os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')),
                       ('大きなリング A (kek-ls_v03_12.lat)', os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))):
        src = open(path, encoding='utf-8').read()
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.01)
        print('\n%s:' % name)
        for fix in (False, True):
            compat.set_bend_dk_fix(fix)
            r = resdrive.chromaticity_numeric(ffs.lat, dmax=3e-3)
            ffs.run_emit(quiet=True)
            g = ffs.I.globals
            tt = ffs.lat.twiss_table(delta=0.02)[-1]
            print('  dk 修正 %-5s ξ = (%.6f, %.6f)  ∂²ν/∂δ² = (%.4f, %.4f)  ν(δ=2%%) = (%.7f, %.7f)  EMITX = %.8e  SIGE = %.8e' % (
                fix, r['d1'][0], r['d1'][1], r['d2'][0], r['d2'][1], tt.mux / 2 / math.pi, tt.muy / 2 / math.pi, float(g['EMITX']), float(g['SIGE'])))
    compat.set_bend_dk_fix(None)


if __name__ == '__main__':
    compat.set_mode('v2')
    part1()
    part2()
    part3()
