# -*- coding: utf-8 -*-
"""v2.3.0 T2: Touschek 寿命の式 (touschek_piwinski) の検証。

(1) 平坦ビームの極限が Bruck・Piwinski の関数 D(ε) (記憶に基づく式) に一致、K64 の当てはめ touschek_f の誤差 (≤ 1.1%)、
(2) 1 次元の式が、3 次元の一般の積分・実験室系の対のモンテカルロ (厳密な相対論の運動学) と一致 (ASTRID2 のスライス)、
(3) ASTRID2 のリング: K64 の表・Piwinski 型の式・保存済みの elegant (example/v2_3/T2_scan_*.json) の比、
(4) 引数 (一様・要素ごとの配列、片側ずつ)、局所の受容 (lma.touschek_lifetime_local_piwinski)、台本の TouschekLifetime[INF, INF, nz] の既定 (v2 = Piwinski、k64 = 表)・Formula->。
"""
import json
import math
import os
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import exp1

from pysad import compat
from pysad.sadscript import Rule, S
from pysad.touschek import (_piw_J, piwinski_slice_rate, touschek_f, touschek_lifetime, touschek_piwinski, touschek_table)
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EX = os.path.join(ROOT, 'example', 'v2_3')
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
sys.path.insert(0, EX)


def D_bruck(e):
    i1 = quad(lambda u: math.log(u) / u * math.exp(-u), e, np.inf)[0]
    return math.sqrt(e) * (-1.5 * math.exp(-e) + e / 2 * i1 + 0.5 * (3 * e - e * math.log(e) + 2) * exp1(e))


def test_flat_limit_is_bruck_piwinski_D():
    """λ2 → 0 で J(c) √(π α2) c² = π D(x)/√x (x = c² α1)。厳密な関数の解析式 (記憶に基づく) との一致"""
    for x in (1e-4, 1e-2, 0.3, 1.0, 3.0):
        J = _piw_J(1.0, x, 1e9, npan=60)
        assert abs(J * math.sqrt(math.pi * 1e9) / (math.pi * D_bruck(x) / math.sqrt(x)) - 1) < 1e-5


def test_k64_fit_error_vs_exact_function():
    """K64 の touschek_f (C(x)/x の当てはめ) と厳密な関数 g(x) = π D(x)/√x の比 f π x/g: x = 1e-6 … 10 で 1.1% 以内"""
    for x in np.logspace(-6, 1, 15):
        r = touschek_f(x) * math.pi * x / (math.pi * D_bruck(x) / math.sqrt(x))
        assert abs(r - 1) < 0.011


@pytest.fixture(scope='module')
def ring():
    old = compat._mode
    compat.set_mode('v2')
    try:
        src = open(ASTRID, encoding='utf-8').read()
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)
        g = ffs.I.globals
        ffs.flags['INTRA'] = True
        g['PBUNCH'] = 1e9
        ffs.flags['RFSW'] = True
        ffs.flags['RAD'] = True
        ffs.run_emit(quiet=True)
        yield st, ffs
    finally:
        compat.set_mode(old)


def test_ring_values_k64_and_piwinski(ring):
    """ASTRID2 (0.58 GeV)、IBS 込み (PBUNCH 1e9)、一様な a = 0.4% / 1.6%: K64 の表 453.45 s、Piwinski 型 427.07 s / 18912.4 s (T2 の値)"""
    st, ffs = ring
    tt = ffs.touschek
    sige = float(ffs.I.globals['SIGE'])
    k64 = touschek_lifetime(tt['table'], None, None, 0.004 / sige, sige)
    piw = touschek_piwinski(tt, a=0.004)['lifetime']
    assert abs(k64 / 453.4547 - 1) < 1e-4 and abs(piw / 427.0707 - 1) < 1e-4
    assert abs(touschek_piwinski(tt, a=0.016)['lifetime'] / 18912.42 - 1) < 1e-4


def test_stored_elegant_scan_agrees_with_independent_formula():
    """保存済みの比較 (example/v2_3/T2_scan_*.json、elegant は全要素を 8 分割して標本化): elegant は独立な式と −1.6% … +0.3% (elegant の損失率が少し大きい)。
    K64 の式は丸いビーム (MINCOUP 1) などで −14%〜+21%"""
    ks = []
    for name in ('astrid2', 'kekls'):
        for r in json.load(open(os.path.join(EX, 'T2_scan_%s.json' % name))):
            e = r['elegant'] / r['piwinski'] - 1
            assert -0.016 < e < 0.003
            ks.append(r['k64'] / r['piwinski'] - 1)
    assert min(ks) < -0.12 and max(ks) > 0.15


def test_argument_forms(ring):
    st, ffs = ring
    tt = ffs.touschek
    nl = tt['rest_info']['nlat']
    r0 = touschek_piwinski(tt, a=0.006)
    r1 = touschek_piwinski(tt, ap=0.006, am=0.006)
    r2 = touschek_piwinski(tt, ap=np.full(nl, 0.006), am=np.full(nl, -0.006))      # 負の値も大きさとして扱う
    assert r0['rate'] == r1['rate'] == r2['rate']
    ra = touschek_piwinski(tt, ap=0.006, am=0.010)['rate']
    rb = touschek_piwinski(tt, ap=0.010, am=0.006)['rate']
    assert abs(ra / rb - 1) < 1e-12 and abs(ra - 0.5 * (r0['rate'] + touschek_piwinski(tt, a=0.010)['rate'])) < 1e-12 * ra       # 片側ずつ半分
    assert abs(r0['rate_per_slice'].sum() / r0['rate'] - 1) < 1e-14


def test_local_acceptance_uniform_equals_uniform(ring):
    from pysad import lma
    st, ffs = ring
    tt = ffs.touschek
    s = np.array([0.0, 10.0, 20.0])
    uni = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [0.006] * 3, [0.006] * 3, [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    res = lma.touschek_lifetime_local_piwinski(tt, uni, setup=st)
    assert abs(res['rate'] / touschek_piwinski(tt, a=0.006)['rate'] - 1) < 1e-12
    half = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [0.006, 0.003, 0.006], [0.006, 0.003, 0.006], [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    assert lma.touschek_lifetime_local_piwinski(tt, half, setup=st)['rate'] > res['rate']


def test_script_default_and_formula_option(ring):
    """TouschekLifetime[INF, INF, nz]: v2 の既定は Piwinski 型、k64 モードと Formula->"K64" は表。nx か ny が有限のときは表 (Formula->"Piwinski" を明示するとエラー)"""
    from pysad.sadscript import SADError
    st, ffs = ring
    tt = ffs.touschek
    sige = float(ffs.I.globals['SIGE'])
    nz = 0.004 / sige
    tab = touschek_lifetime(tt['table'], None, None, nz, sige)
    piw = touschek_piwinski(tt, a=nz * sige)['lifetime']
    old = compat._mode
    try:
        compat.set_mode('v2')
        assert abs(ffs.fn_touschek([None, None, nz]) / piw - 1) < 1e-12
        assert abs(ffs.fn_touschek([None, None, nz, Rule(S('Formula'), 'K64')]) / tab - 1) < 1e-12
        assert abs(ffs.fn_touschek([None, None, nz, Rule(S('Formula'), 'Piwinski')]) / piw - 1) < 1e-12
        assert abs(ffs.fn_touschek([30.0, None, nz]) / touschek_lifetime(tt['table'], 30.0, None, nz, sige) - 1) < 1e-12
        with pytest.raises(SADError):
            ffs.fn_touschek([30.0, None, nz, Rule(S('Formula'), 'Piwinski')])
        compat.set_mode('k64')
        assert abs(ffs.fn_touschek([None, None, nz]) / tab - 1) < 1e-12
    finally:
        compat.set_mode(old)


def test_one_dimensional_formula_vs_3d_and_relativistic_mc(ring):
    """ASTRID2 の 1 スライスで、1 次元の式を 3 次元の一般の積分 (T2_check3d) と実験室系の対のモンテカルロ (T2_mc、厳密な相対論) と比べる (小さな標本。本番は T2_validate: 3D は 1e-3、MC は 0.3% 以内)"""
    from T2_check3d import slice_rate3d
    from T2_mc import mc_rate
    st, ffs = ring
    old = compat._mode
    compat.set_mode('v2')
    try:
        tt = touschek_table(ffs.lat, ffs.momentum(), ffs.emit, 1e9, mincoup=0.01, rfsw=True, radcod=ffs.flags['RADCOD'], verbose=False, record_rest=True)
    finally:
        compat.set_mode(old)
    s = tt['rest'][76]
    info = tt['rest_info']
    for a in (0.004, 0.016):
        r1 = piwinski_slice_rate(s, a, a, info['circ'], 1e9)
        r3, e3 = slice_rate3d(s['al'], None, s['pl'], s['p0'], s['h0'], s['xx'], info['circ'], 1e9, a, a, nu=1024)
        rm, em = mc_rate(s['lab'], s['al'], info['circ'], 1e9, s['p0'], a, a, npairs=2048)
        assert abs(r3 / r1 - 1) < 0.02
        assert abs(rm / r1 - 1) < 0.02
