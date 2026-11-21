# -*- coding: utf-8 -*-
"""v2.4.0 A3: Touschek の損失率にバンチ内の δ の分布 (対の重心の δ_c ~ N(0, σδ²/2)、touschek_bunch.py) を入れるオプション。
リング: ASTRID2 (0.58 GeV、IBS 込み PBUNCH 1e9)。検証: (1) σδ → 0 で touschek_piwinski と一致、(2) a ≫ σδ の展開 R(a)(1 + n(n+1)/4 σδ²/a²) に合う、
(3) 打ち切り floor に依らない、(4) 保存した実験室系の対のモンテカルロ (example/v2_4/A3_mc_bunch_astrid2.out.txt) の δ_c の効果の比と 5% 以内、
(5) 局所の LMA の受け渡し lma.touschek_lifetime_local_piwinski(bunch_delta=True)・台本 TouschekLifetimeLMA[BunchDelta->True]、既定は変わらない。"""
import math
import os
import re

import numpy as np
import pytest

from pysad import compat
from pysad.touschek import piwinski_slice_rate, touschek_piwinski
from pysad.touschek_bunch import piwinski_slice_rate_bunch, slice_sigma_delta, touschek_piwinski_bunch
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


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


def test_zero_spread_is_piwinski_and_sigma_matches_emit(ring):
    st, ffs = ring
    tt = ffs.touschek
    r0 = touschek_piwinski(tt, a=0.006)
    rb = touschek_piwinski_bunch(tt, a=0.006, sigd=0.0)
    assert abs(rb['rate'] / r0['rate'] - 1) < 1e-13
    sl = tt['rest'][76]
    assert abs(slice_sigma_delta(sl) / float(ffs.I.globals['SIGE']) - 1) < 5e-3     # スライスの lab 共分散の σδ = SIGE


def test_increase_follows_expansion_and_is_independent_of_floor(ring):
    """損失率の増え = R(a − δ_c) の平均 ≈ 1 + n(n+1)/4 σδ²/a² (n = −dlnR/dlna は a の差分から)。a = 0.006 (a/σδ = 13.7) で 2 次の展開と 10% 以内、floor 0.1 と 0.5 で 1e-6 以内"""
    st, ffs = ring
    tt = ffs.touschek
    sl = tt['rest'][76]
    info = tt['rest_info']
    sd = slice_sigma_delta(sl)
    a = 0.006
    r = lambda x: piwinski_slice_rate(sl, x, x, info['circ'], 1e9)
    h = 0.02 * a
    n = -(math.log(r(a + h)) - math.log(r(a - h))) / (math.log(a + h) - math.log(a - h))
    inc = piwinski_slice_rate_bunch(sl, a, a, info['circ'], 1e9) / r(a) - 1
    expd = n * (n + 1) / 4.0 * sd ** 2 / a ** 2
    assert inc > 0 and abs(inc / expd - 1) < 0.1
    i1 = piwinski_slice_rate_bunch(sl, a, a, info['circ'], 1e9, floor=0.1)
    i5 = piwinski_slice_rate_bunch(sl, a, a, info['circ'], 1e9, floor=0.5)
    assert abs(i1 / i5 - 1) < 1e-6


def test_matches_saved_monte_carlo_delta_c_effect(ring):
    """保存した対のモンテカルロ (A3_mc_bunch_astrid2.out.txt) の δ_c の効果 (abs/c − 1) の列 (スライス 76・152、a = 0.004・0.006) と、理論 (σδ²/2) が 5% 以内。
    (理論はスライス 76 でこの出力の「σδ²/2: δ_c のみ」の列と同じ)"""
    st, ffs = ring
    tt = ffs.touschek
    info = tt['rest_info']
    txt = open(os.path.join(ROOT, 'example', 'v2_4', 'A3_mc_bunch_astrid2.out.txt'), encoding='utf-8').read().split('\n')
    rows = []
    for i, ln in enumerate(txt):
        m = re.match(r'\s+\[δ_c 基準の kick .*abs/c = ([+-][0-9.]+)\]', ln)
        if m:
            m2 = re.match(r'\s+(\d+)\s+(\d+)\s+([0-9.]+)\s', txt[i + 1])
            rows.append((int(m2.group(1)), float(m2.group(3)), float(m.group(1))))
    assert len(rows) >= 4
    for j, a, mc_inc in rows:
        sl = tt['rest'][j]
        r0 = piwinski_slice_rate(sl, a, a, info['circ'], 1e9)
        th = piwinski_slice_rate_bunch(sl, a, a, info['circ'], 1e9) / r0 - 1
        assert abs(th / mc_inc - 1) < 0.07, (j, a, th, mc_inc)


def test_local_lma_and_script_option(ring):
    from pysad import lma
    st, ffs = ring
    tt = ffs.touschek
    s = np.array([0.0, 10.0, 20.0])
    uni = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [0.006] * 3, [0.006] * 3, [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    r0 = lma.touschek_lifetime_local_piwinski(tt, uni, setup=st)
    rb = lma.touschek_lifetime_local_piwinski(tt, uni, setup=st, bunch_delta=True)
    assert abs(r0['rate'] / touschek_piwinski(tt, a=0.006)['rate'] - 1) < 1e-12        # 既定は変わらない
    assert 0.005 < rb['rate'] / r0['rate'] - 1 < 0.05 and rb['lifetime'] < r0['lifetime']
    assert abs(rb['rate'] / touschek_piwinski_bunch(tt, a=0.006)['rate'] - 1) < 1e-12
