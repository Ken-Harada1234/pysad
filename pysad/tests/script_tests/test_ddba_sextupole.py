# -*- coding: utf-8 -*-
"""R7 (試作 pysad/sxtopt.py): 6 極ファミリーの K2L に対する RDT・ADTS の多項式モデルと、解析式を目的関数にした最適化。
小さな条件 (FODO リング、6 極 2 ファミリー SF・SD、3 GeV) で、(1) 多項式モデルが解析式 (resonance_terms・rdt4・adts) と一致すること、
(2) 1 次の RDT の6 極の分 (強さの変化分) が追跡の 1 周の写像の 2 次のテイラー展開から取り出した値と合うこと、(3) 二次形式のフィットが厳密なこと、
(4) 勾配が有限差分と合うこと、(5) 最適化が走って制約を満たし目的関数を下げること、(6) ソース文字列の K2 の置き換え。
v07_14 のリング (example/ddba/v07_14/lattice/v07_14_5GeV.sad) があれば、20 ファミリーのモデルも解析式と照合する。
条件: pysad k64 の追跡 (conftest の既定)、RF・放射なし、リング先頭が位相の基準。"""
import os

import numpy as np
import pytest

from pysad import rdttrack as T
from pysad import resdrive as R
from pysad import sxtopt as S
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9
V0714 = os.path.join(ROOT, 'example', 'ddba', 'v07_14', 'lattice', 'v07_14_5GeV.sad')


def fodo_src(k_sf=4.0, k_sd=-6.0):
    return FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))


@pytest.fixture(scope='module')
def fodo_model():
    lat = load_sad_string(fodo_src()).build_lattice('RING')
    return lat, S.SxtModel(lat, P0)


def test_families_and_k2l(fodo_model):
    lat, m = fodo_model
    assert m.names == ['SF', 'SD'] and m.F == 2
    # K2L = K2 × L (SAD の SEXT の K2 は積分強さ、pysad は K2/L を持つ)
    assert np.allclose(m.k_ref, [4.0, -6.0])
    lat2 = S.with_k2l(lat, m.names, [1.5, -2.5])
    assert np.allclose(S.sextupole_families(lat2)[2], [1.5, -2.5])
    assert np.allclose(S.sextupole_families(lat)[2], [4.0, -6.0])        # 元は変わらない


def test_model_equals_direct_analytic_formulas(fodo_model):
    """多項式モデルの 1 次の項・3 次の和・4 次の項 (正規形の前)・ADTS が、解析式の直接の呼び出しと 1e-9 で一致 (組み立てたときの K2L と別の K2L の両方)"""
    lat, m = fodo_model
    for k in (m.k_ref, np.array([2.3, -7.1]), np.array([-3.0, 5.0]), np.zeros(2)):
        c = m.check(k)
        assert c['h1'] < 1e-9 and c['h3'] < 1e-9 and c['h4'] < 1e-9, c
        if np.any(k != 0):
            assert c['adts'] < 1e-9, c
        assert c['forms'] < 1e-9, c


def test_first_order_rdt_change_matches_tracked_one_turn_map(fodo_model):
    """6 極の強さの変化に対する 1 次の RDT の変化 (モデル) = 追跡の 1 周の写像の 2 次のテイラー展開から取り出した h の変化。2 次の部分は kick の和なので 1 次の摂動で厳密 (1e-6)"""
    lat, m = fodo_model
    k1 = np.array([4.0, -6.0])
    k0 = np.zeros(2)
    t1 = T.tracked_rdt_table(S.with_k2l(lat, m.names, k1).elements, P0, lat.twiss_table()[0], with_delta=False)['cubic']
    t0 = T.tracked_rdt_table(S.with_k2l(lat, m.names, k0).elements, P0, lat.twiss_table()[0], with_delta=False)['cubic']
    h1, h0 = m.h1(k1), m.h1(k0)
    for nm in S.CUBIC5:
        d_model = h1[nm] - h0[nm]
        d_track = t1[nm] - t0[nm]
        assert abs(d_model - d_track) / abs(d_track) < 1e-6, nm


def test_quad_fit_is_exact():
    rng = np.random.default_rng(0)
    n = 4
    c = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    b = rng.standard_normal((n, 3))
    A = rng.standard_normal((n, n, 3))
    A = A + A.transpose(1, 0, 2)

    def f(k):
        return c + k @ b + np.einsum('i,ijm,j->m', k, A, k)

    cbS = S.quad_fit(f, n)
    k = rng.standard_normal(n)
    assert np.allclose(S.quad_eval(cbS, k), f(k), atol=1e-12)
    # 勾配
    J = S.quad_jac(cbS, k)
    e = np.zeros(n)
    for i in range(n):
        e[:] = 0
        e[i] = 1e-6
        assert np.allclose(J[:, i], (f(k + e) - f(k - e)) / 2e-6, atol=1e-6)


def test_objective_gradient_matches_finite_difference(fodo_model):
    lat, m = fodo_model
    obj = S.SxtObjective(m)
    k = np.array([3.1, -5.2])
    g = obj.grad(k)
    for i in range(2):
        e = np.zeros(2)
        e[i] = 1e-5
        fd = (obj.value(k + e) - obj.value(k - e)) / 2e-5
        assert abs(g[i] - fd) < 1e-5 * abs(fd) + 1e-30


def test_optimize_runs_respects_constraints_and_lowers_objective():
    """2 ファミリーでは ξx・ξy の 2 本の等式で k が決まってしまう (自由度 0) ので、ファミリーを 4 つに増やした FODO (SF・SD を 2 つずつ) で走らせる"""
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=3) SD = (L=0.2 K2=-5) SF2 = (L=0.2 K2=2) SD2 = (L=0.2 K2=-4)')
    src = src.replace('LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);', 'LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);\nLINE CELLB = (SF2 DS QF D2 B D2 QD DS SD2 D1);').replace(
        'LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF 12*CELL 12*CELLB);')
    lat = load_sad_string(src).build_lattice('RING')
    m = S.SxtModel(lat, P0)
    assert m.F == 4
    obj = S.SxtObjective(m, J_ref=1e-6)
    k0 = m.k_ref.copy()
    xi0, M = m.xi_matrix()
    tgt = xi0 + M @ k0
    start = k0 * np.array([1.3, 0.7, 1.2, 0.8])
    start = S.project_to_xi(m, start, tgt)
    assert np.allclose(xi0 + M @ start, tgt, atol=1e-9)
    k, res = S.optimize(obj, start, tgt, kmax=12.0)
    assert np.allclose(xi0 + M @ k, tgt, atol=1e-7)                 # 等式制約
    assert np.max(np.abs(k)) <= 12.0 + 1e-9                         # 上限
    assert obj.value(k) <= obj.value(start) * (1 + 1e-9)            # 目的関数が下がる
    assert obj.value(k) < obj.value(start)


def test_replace_k2_in_source():
    src = fodo_src(4.0, -6.0)
    out = S.replace_k2_in_source(src, ['SF', 'SD'], [1.25, -2.5])
    lat = load_sad_string(out).build_lattice('RING')
    assert np.allclose(S.sextupole_families(lat)[2], [1.25, -2.5])
    with pytest.raises(ValueError):
        S.replace_k2_in_source(src, ['NOPE'], [1.0])


@pytest.mark.skipif(not os.path.exists(V0714), reason='v07_14 のラティスがない')
def test_v07_14_twenty_families_model_matches_analytic():
    """v07_14 (5 GeV、749.5 m、6 極 20 ファミリー) の多項式モデルが解析式と 1e-9 で一致し、ξ の行列と K2L の鏡面対称 (ファミリー数 20) を確かめる"""
    lat = load_sad_file(V0714).build_lattice('RING')
    m = S.SxtModel(lat, 5e9)
    assert m.F == 20
    k = m.k_ref * (1 + 0.1 * np.random.default_rng(1).standard_normal(20))
    c = m.check(k)
    assert c['h1'] < 1e-9 and c['h4'] < 1e-9 and c['adts'] < 1e-9 and c['forms'] < 1e-9, c
    xi0, M = m.xi_matrix()
    assert M.shape == (2, 20)
    assert np.allclose(xi0 + M @ k, m.xi(k))
