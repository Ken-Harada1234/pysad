# -*- coding: utf-8 -*-
"""v2.4.0 A1: ampnf (1 周の写像のべき級数の Birkhoff 正規形) による ν(Jx, Jy, δ) の J²・δ 依存・交差項。
リング: 小さな FODO (3 GeV、105.6 m、24 セル、SF・SD の K2L = 4, -6)。条件: k64、RF・放射なし、結合なし。
検証の基準: (1) J¹ は resdrive.adts の 6 極 2 次 + 8 極 1 次 (解析式) に nslice → ∞ で 1/nslice² で収束、(2) 同じ模型の数値の追跡 + NAFF (初期条件は正規形の J から作る) と
J²・δ・交差項で 1e-8 以内、(3) 偏向磁石の曲率の 3 次の項 (bend_cubic) を入れると 6 極なしの ξ が CALC と 5e-4 以内、入れないと h11001/π の分ずれる。"""
import math
import os

import numpy as np
import pytest

from pysad import ampnf as A
from pysad import resdrive as R
from pysad.parser import load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'


def fodo(k_sf=4.0, k_sd=-6.0):
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))
    return load_sad_string(src).build_lattice('RING')


def test_j1_converges_to_analytic_sextupole_adts():
    """ampnf の J¹ (ADTS) は resdrive.adts の 6 極 2 次の解析式に、厚い 6 極の分割数 n で 1/n² (n を 2 倍で 1/4) で収束する (n=32 で 2e-4 以内)"""
    lat = fodo()
    a = R.adts(lat).sext
    errs = []
    for n in (8, 16, 32):
        r = A.amplitude_tune(lat, D=5, nslice=n)
        errs.append([abs(r.adts[k] - a[k]) / abs(a[k]) for k in ('xx', 'xy', 'yy')])
    for j in range(3):
        assert errs[2][j] < 2e-4
        assert 3.5 < errs[0][j] / errs[1][j] < 4.5 and 3.5 < errs[1][j] / errs[2][j] < 4.5


def test_series_matches_model_tracking_j2_delta_cross():
    """級数 (D=7、bend_cubic なし) と同じ模型の追跡 + NAFF (初期条件は正規形の J・φ から) の差が 1e-8 以内 (J の 2 次・δ・Jδ の項を含む点)。
    J¹ と ξ だけの予想は 1e-6 以上外れる (J²・Jδ・δ² の項が要る)"""
    lat = fodo()
    r = A.amplitude_tune(lat, D=7, nslice=8)
    assert r.min_den > 0.05
    for (jx, jy, d) in [(2e-6, 2e-6, 0.0), (2e-6, 2e-6, 1e-3), (1e-6, 3e-6, -5e-4)]:
        f = A.model_tune(lat, jx, jy, d, nturn=1024, nslice=8, nfres=r)
        s = (r.nu(0, jx, jy, d), r.nu(1, jx, jy, d))
        assert max(abs(f[0] % 1 - s[0] % 1), abs(f[1] % 1 - s[1] % 1)) < 1e-8
        s1 = (r.nu0[0] + r.coef(0, 1, 0) * jx + r.coef(0, 0, 1) * jy + sum(r.coef(0, 0, 0, p) * d ** p for p in range(1, 7)),
              r.nu0[1] + r.coef(1, 1, 0) * jx + r.coef(1, 0, 1) * jy + sum(r.coef(1, 0, 0, p) * d ** p for p in range(1, 7)))
        assert max(abs(f[0] % 1 - s1[0]), abs(f[1] % 1 - s1[1])) > 1e-6


def test_series_is_converged_in_order_and_normal_form_symmetry():
    """次数 D = 5, 6, 7 で ν(δ=1e-3) が 1e-8 以内 (級数が収束している)。J² の係数は D=5 と 7 で同じ (J² は D=5 で入る)"""
    lat = fodo()
    rs = [A.amplitude_tune(lat, D=D, nslice=8) for D in (5, 6, 7)]
    for r in rs[1:]:
        assert abs(r.nu(0, 0, 0, 1e-3) - rs[0].nu(0, 0, 0, 1e-3)) < 1e-7
    for k in ('x20', 'x11', 'x02', 'y20', 'y11', 'y02'):
        assert abs(rs[0].adts2[k] - rs[2].adts2[k]) < 1e-9 * abs(rs[2].adts2[k]) + 1e-6


def test_bend_cubic_makes_chromaticity_match_calc():
    """6 極なしの ξ: 偏向磁石の曲率の 3 次の項 (bend_cubic=True) を入れた模型が CALC の ξ と 5e-4 以内。入れない模型との差が h11001/π・h00111/π (偏向磁石の分、R4)"""
    lat = fodo(0.0, 0.0)
    ref = R.chromaticity_numeric(lat, order=1)['d1']
    r0 = A.amplitude_tune(lat, D=5, nslice=4)
    r1 = A.amplitude_tune(lat, D=5, nslice=4, bend_cubic=True, nslice_bend=8)
    bend = R.chromaticity_bend(lat)
    for i in (0, 1):
        assert abs(r1.chrom(i, 1) - ref[i]) < 5e-4
        assert abs(r0.chrom(i, 1) + bend[i] - ref[i]) < 5e-4
        assert abs(r0.chrom(i, 1) - ref[i]) > 0.1


def test_bend_cubic_series_matches_model_tracking():
    """bend_cubic あり (D=5) の級数と、同じ模型 (曲率の 3 次の項の流れつき) の追跡 + NAFF の差が 5e-7 以内 (D=5 の打ち切り誤差。D=7 では 1e-9)"""
    lat = fodo()
    r = A.amplitude_tune(lat, D=5, nslice=4, bend_cubic=True, nslice_bend=4)
    jx, jy, d = 1e-6, 1e-6, 1e-3
    f = A.model_tune(lat, jx, jy, d, nturn=1024, nslice=4, nfres=r, bend_cubic=True, nslice_bend=4, Dnl=5)
    assert max(abs(f[0] % 1 - r.nu(0, jx, jy, d) % 1), abs(f[1] % 1 - r.nu(1, jx, jy, d) % 1)) < 5e-7


def test_near_resonance_reported():
    """作動点が 3νx = 整数の近くだと、共鳴の分母の最小値 (min_den) が小さく、警告のリスト near に項が出る"""
    lat = fodo()
    r = A.amplitude_tune(lat, D=7, nslice=4)
    assert r.min_den < 0.1 and isinstance(r.near, list)


def test_full_model_matches_adts_total_and_calc_chromaticity():
    """運動学的な項・端部のフリンジ・偏向磁石の曲率の 3 次の項を全部入れた模型 (bend_cubic・kin・fringe、D = 5) で、6 極なしの FODO の
    (1) J¹ = resdrive.adts().total (運動学 + フリンジ + 偏向磁石) と 3e-3 以内、(2) ξ1・ξ2 が CALC と 1e-3・3% 以内 (x・y)、ξ3 が CALC と 3% (x)。
    端部のフリンジの強さは K1/(1+δ) (FRINGE_POWER = 1)。δ に依らない (0) と ξ3 が 22% ずれる"""
    lat = fodo(0.0, 0.0)
    an = R.adts(lat)
    ref = R.chromaticity_numeric(lat, order=3)
    r = A.amplitude_tune(lat, D=5, nslice=8, nslice_bend=8, bend_cubic=True, kin=True, fringe=True)
    for k in ('xx', 'xy', 'yy'):
        assert abs(r.adts[k] - an.total[k]) / abs(an.total[k]) < 3e-3, k
    for i in (0, 1):
        assert abs(r.chrom(i, 1) - ref['d1'][i]) < 1e-3
        assert abs(r.chrom(i, 2) / ref['d2'][i] - 1) < 0.03
    assert abs(r.chrom(0, 3) / ref['d3'][0] - 1) < 0.03
    A.FRINGE_POWER = 0
    A._FLOW.clear()
    try:
        r0 = A.amplitude_tune(lat, D=5, nslice=8, nslice_bend=8, bend_cubic=True, kin=True, fringe=True)
    finally:
        A.FRINGE_POWER = 1
        A._FLOW.clear()
    assert abs(r0.chrom(0, 3) / ref['d3'][0] - 1) > 0.15
