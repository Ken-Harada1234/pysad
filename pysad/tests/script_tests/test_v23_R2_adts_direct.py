# -*- coding: utf-8 -*-
"""v2.3.0 R2: ADTS の直接の検証 (実際の強さの追跡 + NAFF + J の多項式フィット、pysad/adtsdirect.py) と、
運動学的な項・四極の端部のフリンジの 1 次の ADTS (resdrive.kinematic_adts・quad_fringe_adts)。
追跡は k64 (conftest の既定)。小さな FODO (3 GeV、105.6 m、24 セル、SF・SD の K2L = 4, -6)、ASTRID2 (0.58 GeV、45.7 m)。
条件: RF 空洞は止める (rfsw=False)、放射なし、δ=0、結合なし、NAFF は Hann p=1、1024 周 (試験では 512 周のものもある)。"""
import json
import math
import os

import numpy as np
import pytest

from pysad import adtsdirect as D
from pysad import resdrive as R
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9


def fodo(k_sf=4.0, k_sd=-6.0, nobend=False):
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))
    if nobend:
        src = src.replace('ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747', 'ANGLE=0 E1=0 E2=0')
    return load_sad_string(src).build_lattice('RING')


def rel(a, b):
    return abs(a - b) / abs(b)


def test_fit_tune_surface_recovers_coefficients():
    """合成した ν(Jx, Jy) (3 次の多項式 + 小さな雑音) から、1 次・2 次・3 次の係数を戻す。±位相の対の平均、νx は Jx>0、νy は Jy>0 の点だけ"""
    jx, jy = D.grid_actions(1e-6, 5)
    cx = {(1, 0): 300.0, (0, 1): -450.0, (2, 0): 2e7, (1, 1): -5e7, (0, 2): 1e7, (3, 0): -1e12}
    cy = {(1, 0): -450.0, (0, 1): 500.0, (2, 0): 3e6, (1, 1): 2e7, (0, 2): -4e7, (0, 3): 5e11}

    def f(c):
        return -0.3 + sum(v * jx ** a * jy ** b for (a, b), v in c.items())
    nu = np.column_stack([f(cx), f(cy) * 0 + (0.14 + sum(v * jx ** a * jy ** b for (a, b), v in cy.items()))])
    nu[:, 0] += 0.0
    r = D.adts_from_tunes(jx, jy, nu, order=3)
    assert abs(r['xx'] - 300.0) < 1e-6 * 300 and abs(r['xy'] + 450.0) < 1e-6 * 450
    assert abs(r['yx'] + 450.0) < 1e-6 * 450 and abs(r['yy'] - 500.0) < 1e-6 * 500
    assert abs(r['fitx']['c'][(2, 0)] - 2e7) < 1e-4 * 2e7 and abs(r['fity']['c'][(0, 2)] + 4e7) < 1e-4 * 4e7
    # 次数が足りないと (order=1) 大きな J の曲率で ADTS がずれる (J の範囲を変える意味)
    r1 = D.adts_from_tunes(jx, jy, nu, order=1)
    assert abs(r1['xx'] - 300.0) > 1e-3 * 300


def test_fit_delta_curve():
    d = np.linspace(-2e-3, 2e-3, 9)
    nu = 0.3 + 5.0 * d - 2.0e3 * d ** 2 + 4e5 * d ** 3
    c = D.fit_delta_curve(d, nu, 3)
    assert abs(c['d1'] - 5.0) < 1e-8 and abs(c['d2'] + 4e3) < 1e-5 and abs(c['d3'] - 2.4e6) < 1


def test_paraxial_tracker_matches_analytic_sextupole_adts():
    """線形行列 + 薄い 6 極の追跡 (解析式と同じ模型) の直接のフィット = 解析の 6 極 2 次 (FODO、2e-4 以内。薄レンズ 30 分割の誤差を含む)"""
    lat = fodo()
    A = R.adts(lat)
    tr = D.paraxial_tracker(lat, 30)
    r = D.direct_adts(lat, P0, 1e-7, order=2, nturn=512, nstep=3, tracker=tr)
    for k, key in (('xx', 'xx'), ('xy', 'xy_mean'), ('yy', 'yy')):
        assert rel(r[key], A.sext[k]) < 2e-4


def test_kinematic_and_fringe_without_sextupoles():
    """6 極なし (K2=0) の FODO の ADTS = 運動学的な項 (p⁴/8) + 四極の端部のフリンジ (ttfrin4、K64 の既定 ON): 解析と追跡が 1% 以内 (J = 1e-7)。
    端部のフリンジを切る (DISFRIN=1) と運動学的な項だけになる (1 次の式と 3% 以内、絶対値で 0.3。yy は 6% 残る: SAD の四極本体の写像の 4 次の項と未確認の差)"""
    for nobend in (True, False):
        lat = fodo(0.0, 0.0, nobend=nobend)
        A = R.adts(lat)
        assert abs(A.sext['xx']) < 1e-12 and abs(A.octu['xx']) < 1e-12
        r = D.direct_adts(lat, P0, 1e-7, order=2, nturn=512, nstep=3)
        for k, key in (('xx', 'xx'), ('xy', 'xy_mean'), ('yy', 'yy')):
            assert rel(r[key], A.total[k]) < 1e-2, (nobend, k, r[key], A.total[k])
    lat = fodo(0.0, 0.0, nobend=True)
    for e in lat.elements:
        if type(e).__name__ in ('Quadrupole', 'Bend'):
            e.DISFRIN = 1.0
    A = R.adts(lat)
    assert A.fringe['xx'] == 0.0 and A.kin['xx'] > 0
    r = D.direct_adts(lat, P0, 1e-7, order=2, nturn=512, nstep=3)
    for k, key in (('xx', 'xx'), ('xy', 'xy_mean'), ('yy', 'yy')):
        assert abs(r[key] - A.kin[k]) < 0.3, (k, r[key], A.kin[k])


def test_quad_fringe_formula_matches_map():
    """端部のフリンジの 1 次の式 (σ K1/12 [(x³+3xy²)px − (y³+3x²y)py]) が、pysad の ttfrin4 の写像の 3 次の項と一致する (入口・出口、K1 の符号)"""
    from pysad.track6n import track_particles_fast
    from pysad.elements import Quadrupole
    # 短い四極 (L=0.3)。入口の端だけ・出口の端だけを FRINGE で切り替え、レンズの本体は K1=0 に近づけて端の写像だけを見る
    for k1 in (0.7, -0.6):
        q = Quadrupole(name='Q', L=1e-6, K1=k1)
        q.FRINGE = 1.0
        a = 3e-3
        c = np.array([[a, 0.0, 0.7 * a, 0.0, 0.0, 0.0]])
        out = track_particles_fast([q], c, 1, P0, 0.0, False, False, False)
        # 入口のフリンジ (px=py=0 のとき): Δx = (K1/12)(x³+3xy²)、Δy = −(K1/12)(y³+3x²y)
        x, y = c[0, 0], c[0, 2]
        # 薄い本体の寄与 (L = 1e-6、K1L = 1e-6 K1 ほぼ 0) と出口 (FRINGE=1 → 入口のみ) を除く
        assert abs((out[0, 0] - x) - k1 / 12 * (x ** 3 + 3 * x * y * y)) < 1e-3 * abs(k1 / 12 * (x ** 3 + 3 * x * y * y)) + 1e-14
        assert abs((out[0, 2] - y) + k1 / 12 * (y ** 3 + 3 * x * x * y)) < 1e-3 * abs(k1 / 12 * (y ** 3 + 3 * x * x * y)) + 1e-14


def test_adts_result_total_and_ffs_function():
    lat = fodo()
    A = R.adts(lat)
    assert abs(A.total['xx'] - (A.sext['xx'] + A.octu['xx'] + A.kin['xx'] + A.fringe['xx'] + A.bend['xx'] + A.cross['xx'])) < 1e-12      # v2.3.0 R3: 偏向磁石の 3 次の項 (既定 ON) も入る
    A2 = R.adts(lat, bend_cubic=False)
    assert abs(A2.total['xx'] - (A.sext['xx'] + A.octu['xx'] + A.kin['xx'] + A.fringe['xx'])) < 1e-12 and A2.bend['xx'] == 0.0
    assert A.dnux_dJx == pytest.approx(A.sext['xx'] + A.octu['xx'])          # 従来の合計 (6 極 + 8 極) は変えない
    A0 = R.adts(lat, include_kinematic=False)
    assert A0.kin['xx'] == 0.0 and A0.fringe['xx'] == 0.0
    assert 'Total' in A.table() or '全部入り' in A.table()


def test_direct_fodo_matches_stored_gs_if_present():
    """保存した gs (NORFSW、70 粒子) の直接のフィットが、pysad の追跡と 1e-4 以内 (gs が走る環境で作った出力 example/14_nonlinear_dynamics/adts_rdt/adts_direct_fodo.gs.json)"""
    f = os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'adts_rdt', 'adts_direct_fodo.gs.json')
    if not os.path.exists(f):
        pytest.skip('gs の保存出力がない')
    d = json.load(open(f))
    for k in ('xx', 'xy', 'yy'):
        assert rel(d['gs'][k], d['pysad'][k]) < 2e-4


def test_tune_vs_delta_matches_calc_fodo():
    """微小振幅の NAFF のチューンの δ 依存 (追跡) が CALC の線形のチューンの δ 依存と合う (FODO、∂ν/∂δ・∂²ν/∂δ² を 1e-3 の相対で。6 極の 2 次の ADTS は J=1e-12 で無視できる)"""
    lat = fodo()
    d0 = 2e-3
    nodes = d0 * np.cos(np.pi * (np.arange(9) + 0.5) / 9)
    nut, ref = D.tune_vs_delta(lat, P0, nodes, j=1e-12, nturn=512)
    for p in (0, 1):
        val = D.fr(nut[:, p] - D.fr(ref[:, p])) + D.fr(ref[:, p])
        fit = D.fit_delta_curve(nodes, val, 4)
        cal = D.fit_delta_curve(nodes, D.fr(ref[:, p]), 4)
        assert abs(fit['d1'] - cal['d1']) < 2e-3 * abs(cal['d1']) + 1e-4
        assert abs(fit['d2'] - cal['d2']) < 2e-2 * abs(cal['d2']) + 1.0
