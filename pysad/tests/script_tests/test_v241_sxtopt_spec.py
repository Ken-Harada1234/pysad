# -*- coding: utf-8 -*-
"""v2.4.1: pysad/sxtopt.py の項の指定 (SxtSpec・SxtSpecObjective・design_k2l)・作動点の周りの共鳴 (near_resonances・resonance_warnings) の試験。
(1) 既定の指定 (legacy_spec) が従来の SxtObjective と同じ値・勾配、(2) 指定の書き方 (辞書・文字列・JSON・規則の列) が同じ結果、(3) 勾配が有限差分と合う (重み・目標値・4 次の項・ADTS・ξ を入れた指定)、
(4) 制約 (線形の等式・非線形の等式・tol つきの不等式) が満たされる、(5) 'ignore' の項が目的に入らない、(6) 共鳴の一覧と警告、(7) 不正な名前の拒否。
条件: FODO リング (6 極 4 ファミリー、3 GeV)。pysad k64 の追跡の設定 (conftest の既定)。"""
import os
import warnings

import numpy as np
import pytest

from pysad import sxtopt as S
from pysad.parser import load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9


@pytest.fixture(scope='module')
def model4():
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=3) SD = (L=0.2 K2=-5) SF2 = (L=0.2 K2=2) SD2 = (L=0.2 K2=-4)')
    src = src.replace('LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);', 'LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);\nLINE CELLB = (SF2 DS QF D2 B D2 QD DS SD2 D1);').replace(
        'LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF 12*CELL 12*CELLB);')
    lat = load_sad_string(src).build_lattice('RING')
    return S.SxtModel(lat, P0)


def _fd(obj, k, xt, i, h=1e-5):
    e = np.zeros(len(k))
    e[i] = h
    return (obj.value(k + e, xt) - obj.value(k - e, xt)) / (2 * h)


def test_default_spec_equals_legacy_objective(model4):
    m = model4
    k = np.array([3.1, -5.2, 2.2, -3.7])
    old = S.SxtObjective(m)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        new = S.SxtSpecObjective(m)                       # spec=None → legacy_spec()
    xt = np.array(m.xi(k))
    assert abs(new.value(k, xt) - old.value(k)) <= 1e-12 * old.value(k)
    assert np.allclose(new.grad(k, xt), old.grad(k), rtol=1e-10, atol=0)
    # 重みを変えた旧指定も同じ
    old2 = S.SxtObjective(m, weights=dict(adts=5.0, w4=0.3), terms=('w3', 'w4', 'adts'))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        new2 = S.SxtSpecObjective(m, S.legacy_spec(weights=dict(adts=5.0, w4=0.3), terms=('w3', 'w4', 'adts')))
    assert abs(new2.value(k, xt) - old2.value(k)) <= 1e-12 * old2.value(k)
    assert np.allclose(new2.grad(k, xt), old2.grad(k), rtol=1e-10, atol=0)


def test_spec_forms_are_equivalent(model4):
    d = {'adts': 'objective', 'h10020': {'mode': 'objective', 'weight': 4, 'target': 0.0}, 'xi': ('constraint', 1, None)}
    s1 = S.SxtSpec.parse(d)
    s2 = S.SxtSpec.parse('adts=objective, h10020=objective:4:0, xi=constraint')
    s3 = S.SxtSpec.parse(s1.to_json())
    s4 = S.SxtSpec.parse([('adts', 'objective'), ('h10020', ('objective', 4, 0.0)), ('xi', ('constraint', 1, None))])
    assert s1.to_dict() == s2.to_dict() == s3.to_dict() == s4.to_dict()
    assert S.SxtSpec.parse('h21000=objective:2:1+2j').terms['h21000'].target == complex(1, 2)
    k = np.array([3.1, -5.2, 2.2, -3.7])
    xt = np.array(model4.xi(k))
    vals = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for s in (s1, s2, s3, s4):
            vals.append(S.SxtSpecObjective(model4, s).value(k, xt))
    assert np.allclose(vals, vals[0], rtol=1e-13)


def test_gradient_with_weights_targets_and_terms(model4):
    m = model4
    spec = {'h21000': ('objective', 2.0, 1e-9 + 2e-9j), 'h10020': ('objective', 0.5), 'h31000': ('objective', 3.0, 1e-3),
            'adts_xx': ('objective', 1.0, 1.0e3), 'adts_xy': ('objective', 2.0), 'h20001': ('objective', 1.0),
            'xi_x': ('objective', 1e-3, 1.0), 'xi_y': ('ignore',), 'w4': ('ignore',)}
    k = np.array([3.1, -5.2, 2.2, -3.7])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obj = S.SxtSpecObjective(m, spec)
    g = obj.grad(k)
    assert np.all(np.isfinite(g))
    for i in range(4):
        fd = _fd(obj, k, None, i)
        assert abs(g[i] - fd) < 1e-5 * abs(fd) + 1e-30, (i, g[i], fd)
    names = {r['name'] for r in obj.table(k)}
    assert 'xi_y' not in names and 'h31000' in names and 'adts_xx' in names


def test_ignore_terms_do_not_enter(model4):
    k = np.array([3.1, -5.2, 2.2, -3.7])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        only_adts = S.SxtSpecObjective(model4, 'adts=objective')
        both = S.SxtSpecObjective(model4, 'adts=objective, w3=objective')
        w0 = S.SxtSpecObjective(model4, 'adts=objective:0, w3=ignore')
    assert both.value(k, [0, 0]) > only_adts.value(k, [0, 0]) > 0
    assert w0.value(k, [0, 0]) == 0.0


def test_constraints_are_satisfied(model4):
    m = model4
    k0 = np.array([3.0, -5.0, 2.0, -4.0])
    xt = np.array(m.xi(k0))
    # 目標値は ξ = xt の多様体の別の点 k1 の ADTS xy (実現できる値)
    k1 = S.project_to_xi(m, k0 * np.array([1.2, 0.9, 1.1, 0.85]), xt)
    tgt = m.adts(k1)['xy']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # ξ の等式 (線形、自動) + ADTS xy = tgt (非線形の等式)。目的は ADTS xx・yy
        r = S.design_k2l(m, {'adts_xx': ('objective',), 'adts_yy': ('objective',), 'adts_xy': ('constraint', 1.0, tgt)}, k0=k0, xi_target=xt, kmax=40.0)
    k = r['k']
    assert r['success']
    assert np.allclose(m.xi(k), xt, atol=1e-7)
    assert abs(m.adts(k)['xy'] - tgt) < 1e-6 * max(abs(tgt), 1.0)
    # tol つきの不等式: |h21000| ≤ tol
    h0 = abs(m.h1(k0)['h21000'])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r2 = S.design_k2l(m, {'adts': ('objective',), 'h21000': ('constraint', 1.0, 0.0, 0.5 * h0)}, k0=k0, xi_target=xt, kmax=40.0)
    assert abs(m.h1(r2['k'])['h21000']) <= 0.5 * h0 * (1 + 1e-6)
    assert np.allclose(m.xi(r2['k']), xt, atol=1e-7)
    assert np.max(np.abs(r2['k'])) <= 40.0 + 1e-9


def test_design_returns_table_and_lowers_objective(model4):
    m = model4
    k0 = np.array([3.0, -5.0, 2.0, -4.0])
    start = S.project_to_xi(m, k0 * np.array([1.3, 0.8, 1.2, 0.7]), m.xi(k0))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = S.design_k2l(m, 'adts=objective, h10020=objective:4', k0=start, xi_target=m.xi(k0), kmax=12.0)
        v0 = S.SxtSpecObjective(m, 'adts=objective, h10020=objective:4').value(start, m.xi(k0))
    assert r['objective'] < v0
    assert {x['name'] for x in r['table']} >= {'adts_xx', 'adts_xy', 'adts_yy', 'h10020', 'xi_x', 'xi_y'}
    assert r['constraint_violation'] < 1e-6


def test_near_resonances_and_warnings():
    res = S.near_resonances(10.505, 5.20, max_order=4, radius=0.05)
    first = res[0]
    assert (first['a'], first['b'], first['p']) == (2, 0, 21) and first['order'] == 2
    assert abs(first["dist"] - 0.005) < 1e-9
    names = [t[0] for t in first['terms']]
    assert 'h20001' in names and 'h31000' in names            # 半整数 2νx は δ の項と 4 次の項が駆動
    # midplane 対称: y の奇数次の項 (h00030 など) は駆動項に出ない
    allnames = {t[0] for d in S.near_resonances(10.3, 5.1, max_order=4, radius=0.5) for t in d['terms']}
    assert 'h00030' not in allnames and 'h30000' in allnames
    w = S.resonance_warnings(10.505, 5.20)
    assert w and '2 次共鳴' in w[0]
    assert all('2 次共鳴' not in x for x in S.resonance_warnings(10.31, 5.17))
    txt = S.format_near_resonances(res)
    assert 'h31000' in txt


def test_bad_names_rejected(model4):
    with pytest.raises(ValueError):
        S.SxtSpecObjective(model4, 'h99999=objective')
    with pytest.raises(ValueError):
        S.SxtSpecObjective(model4, 'h22000=objective')        # 対角は adts_* で指定
    with pytest.raises(ValueError):
        S.SxtSpec.parse('adts=maximize')
    with pytest.raises(ValueError):
        S.TermSpec('objective', -1.0)


def test_error_sensitivity_terms(model4):
    """'err_*' (K2L のランダム誤差による 1 次の項の rms の期待値): 台ごとの寄与のモンテカルロ (400 回) と 15% 以内で一致し、勾配が有限差分と合い、誤差の大きさに比例する。"""
    m = model4
    k = np.array([3.1, -5.2, 2.2, -3.7])
    eps = 0.03
    pred = m.err_rms(k, eps)
    nm_idx = {}
    for ir, nm in enumerate(m.r1.elem_names):
        nm_idx.setdefault(nm, []).append(ir)
    rows = []
    for f, n in enumerate(m.names):
        for ir in nm_idx[n]:
            rows.append(m.r1.per_element[ir] / m.k_ref[f] * k[f])
    rows = np.array(rows)
    rng = np.random.default_rng(3)
    dh = np.array([eps * rng.standard_normal(len(rows)) @ rows for _ in range(400)])
    for t in S.CUBIC5:
        i = m.terms1.index(t)
        mc = np.sqrt(np.mean(np.abs(dh[:, i]) ** 2))
        assert abs(mc - pred[t]) < 0.15 * pred[t], (t, mc, pred[t])
    assert np.isclose(m.err_rms(k, 2 * eps)['h21000'], 2 * pred['h21000'])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obj = S.SxtSpecObjective(m, {'err3': ('objective', 1.0), 'adts': ('objective', 0.5), 'xi_x': ('ignore',)}, err_rms=eps)
    g = obj.grad(k, None)
    for i in range(4):
        fd = _fd(obj, k, None, i)
        assert abs(g[i] - fd) < 1e-5 * abs(fd) + 1e-30
    assert {r['name'] for r in obj.table(k)} >= {'err_h21000', 'err_h10110', 'adts_xy'}


def test_trust_region_dk_limits_change(model4):
    m = model4
    k0 = np.array([3.0, -5.0, 2.0, -4.0])
    xt = np.array(m.xi(k0))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = S.design_k2l(m, 'adts=objective', k0=k0, xi_target=xt, kmax=40.0, dk=0.5)
    assert np.max(np.abs(r['k'] - k0)) <= 0.5 + 1e-9
