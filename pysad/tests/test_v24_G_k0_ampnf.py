# -*- coding: utf-8 -*-
"""v2.4.0 G: K0≠0 の偏向磁石 (MPW・挿入光源・逆曲げ・傾けた MPW) を ampnf (ν(J, δ)・ADTS・共鳴駆動項の TPSA の正規形) に入れる。
磁石の写像は K64 の追跡の閉じた式 (tbendcore1・tsteer) を、磁石内の軌道のまわりで TPSA の演算で評価したもの (pysad/ampnf_k0.py)。独立な基準: (1) 追跡の素子の写像そのもの
(軌道を持つ入口・δ≠0・E1/E2・ANGLE=0/≠0)、(2) pysad の実際の追跡 + NAFF、(3) 追跡の 1 周の写像の 2 次のテイラー展開から取り出した h、(4) K0 → 0 で従来の K0=0 の模型、
(5) ROTATE = π と K0 の符号反転、(6) 摂動の次数 (K0²、傾き θ)。リング: 小さな FODO (3 GeV、24 セル、6 極 K2 = 4・−6) の 1 セルの直線部に 3 極の挿入光源 (K0 = a・(1, −2, 1) rad)。"""
import math
import os
import sys

import numpy as np
import pytest

from pysad import ampnf as A
from pysad import ampnf_k0 as K
from pysad import compat
from pysad import rdttrack as T
from pysad.elements import Bend

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'normal_form'))
import k0_ampnf as G  # noqa: E402

OPT = G.OPT
CASES = {'steer': dict(L=0.3, K0=0.06), 'steer_edge': dict(L=0.3, K0=0.06, E1=0.1, E2=-0.05),
         'bend': dict(L=0.5, ANGLE=0.1, K0=0.03), 'bend_edge': dict(L=0.5, ANGLE=0.1, K0=-0.02, E1=0.04, E2=0.03),
         'reverse': dict(L=0.5, ANGLE=-0.1, K0=0.05, E1=0.04, E2=0.03), 'weak_field': dict(L=0.5, ANGLE=0.1, K0=-0.08),
         # v2.4.0 G の続き: K1≠0 (tbendi の積分器)・F1 のフリンジ・据え付け (DX・DY・ROTATE・DROTATE・CHI2)
         'k1': dict(L=0.5, ANGLE=0.1, K0=0.03, K1=0.4), 'k1_edge': dict(L=0.5, ANGLE=0.1, K0=-0.02, K1=-0.3, E1=0.04, E2=0.03),
         'f1': dict(L=0.5, ANGLE=0.1, K0=0.03, F1=0.2, FRINGE=1.0), 'f1_steer': dict(L=0.3, K0=0.05, F1=0.2, FRINGE=1.0, E1=0.1, E2=0.05),
         'k1_f1': dict(L=0.5, ANGLE=0.1, K0=0.03, K1=0.4, F1=0.2, FRINGE=1.0),
         'k1_f1_edge': dict(L=0.5, ANGLE=0.1, K0=-0.02, K1=-0.3, E1=0.04, E2=0.03, F1=0.15, FB1=0.05, FRINGE=1.0),
         'k1_fb2': dict(L=0.5, ANGLE=0.1, K0=0.01, K1=0.2, FB2=0.2, FRINGE=-2.0),
         'shift': dict(L=0.5, ANGLE=0.1, K0=0.03, DX=2e-3, DY=-1e-3),
         'roll_shift': dict(L=0.5, ANGLE=0.1, K0=0.03, ROTATE=0.3, DX=1e-3, DY=1e-3, DROTATE=0.01, CHI2=0.005, E1=0.02)}


def _mk(kw):
    """Bend (dataclass のフィールド以外 (F1・FRINGE・DX・DY・DROTATE・CHI2・ROTATE) は属性として足す)"""
    base = {k: v for k, v in kw.items() if k in ('L', 'ANGLE', 'K0', 'K1', 'E1', 'E2')}
    el = Bend(name='B', **base)
    for k, v in kw.items():
        if k not in base:
            setattr(el, k, v)
    return el


def test_tpsa_functions_are_exact_series():
    """TP の sqrt・1/x・asin が、定数部分 a0 のまわりの厳密な級数 (数値の関数と 1e-5。次数 5 の打ち切り)"""
    tab = A.get_tab(5, 5)
    u = np.array([0.7, 0.0, -0.4, 0.0, 0.0])
    mono = np.prod(u[None, :] ** tab.mons, axis=1)
    for a0 in (0.3, 0.45, 0.9):
        arr = 0.05 * (tab.var(0) + 0.5 * tab.var(2))
        arr[0] = a0
        x = K.TP(tab, arr)
        z = a0 + 0.05 * (0.7 - 0.5 * 0.4)
        assert abs(np.dot(x.sqrt().a, mono) - np.sqrt(z)) < 1e-5
        assert abs(np.dot(x.inv().a, mono) - 1.0 / z) < 1e-5
        assert abs(np.dot(x.asin().a, mono) - np.arcsin(z)) < 1e-5


@pytest.mark.parametrize('name', list(CASES))
@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('orbit,delta', [((0, 0, 0, 0), 0.0), ((2e-3, -0.02, 1e-3, 0.01), 0.0), ((0, 0, 0, 0), 2e-3), ((-3e-3, 0.03, 0.0, 0.0), -1e-3)])
def test_element_series_equals_tracker_map(name, mode, orbit, delta):
    """K0≠0 の磁石 1 台の写像 (TPSA 評価) が、pysad の追跡の素子の写像と、入口の軌道のまわりの振幅 1e-3 の点で 1e-13 以内 (δ≠0・E1/E2・軌道の傾き・ANGLE=0/≠0。k64・v2 とも。
    ANGLE=0 の出口の端 E2 は v2 だけ常に入る (A15) ので、模型も追跡と同じ互換モードに従う)"""
    old = compat.mode()
    compat.set_mode(mode)
    try:
        el = _mk(CASES[name])
        e, _ = K.check_against_tracker(el, amp=1e-3, orbit=orbit, delta=delta)
    finally:
        compat.set_mode(old)
    assert e < 1e-13


def test_unsupported_k0_bends_raise():
    """未対応の組み合わせは NotImplementedError: ANGLE=0 で K1≠0 かつ K0≠0 の F1・FB1・FB2 (K64 の tmulti。gs が SIGSEGV で基準がなく、pysad の tmulti は積分器の精度が 4e-7〜1.6e-5)。
    ANGLE≠0 の K1≠0 かつ K0≠0 の F1・FB は tbendi の fb の項を入れた (element 試験の k1_f1 など)。ANGLE=0 の CHI2 は K64 の追跡が使わないので無視する (`ampnf.build_ops` が警告)"""
    with pytest.raises(NotImplementedError):
        K.bend_k0_map(_mk(dict(L=0.5, K0=0.01, K1=0.1, F1=0.1, FRINGE=1.0)))
    K.bend_k0_map(_mk(dict(L=0.5, ANGLE=0.1, K0=0.01, K1=0.1, F1=0.1, FRINGE=1.0)))      # tbendi の fb: 例外なし
    K.bend_k0_map(_mk(dict(L=0.5, K0=0.01, CHI2=0.01)))          # 無視 (例外なし)
    ops, warn = A.build_ops(type('L', (), {'elements': [_mk(dict(L=0.5, K0=0.01, CHI2=0.01))]})(), 3, 1e-2, 8)
    assert any('CHI2' in w for w in warn)


@pytest.mark.parametrize('kw', [dict(L=0.3, K0=0.06, K1=0.1), dict(L=0.3, K0=-0.12, K1=0.1), dict(L=0.5, K0=0.03, K1=0.4)])
def test_straight_k0_k1_magnet_equals_independent_flow(kw):
    """ANGLE = 0 で K1≠0 かつ K0≠0 の直線の磁石 (K64 の tmulti。gs は SIGSEGV で基準なし): 厳密な Hamiltonian の対称な 4 次の分割の TPSA (`_straight_k1_core`、Lee–Whiting の四極のフリンジと
    双極子の Maxwell の端つき) が、独立な運動方程式の積分 (四極 + 双極子の Maxwell を満たす場、端の幅 w、`k0_independent_flow.lorentz_straight_k1`) と w に比例以上で収束する
    (w = 4e-3 → 1e-3 で 1.6e-7 → 1e-8 以下)。pysad の tmulti の追跡は 4e-7〜1.6e-5 外れる (積分器の精度。K0 = −0.12 の極で最大)"""
    import k0_independent_flow as F
    tab = A.get_tab(5, 5)
    sd, s = 1e-2, 1e-3
    rng = np.random.default_rng(3)
    pts = np.column_stack([rng.uniform(-1e-3, 1e-3, 3), rng.uniform(-1e-3, 1e-3, 3), rng.uniform(-1.5e-4, 1.5e-4, 3), rng.uniform(-1e-3, 1e-3, 3)])
    el = _mk(kw)
    S1 = K.bend_k0_map(el)([s * tab.var(j) for j in range(4)], tab, sd)
    u = np.column_stack([pts[:, j] / s for j in range(4)] + [np.zeros(3)])
    mono = np.prod(u[:, None, :] ** tab.mons[None, :, :], axis=2)
    got = np.column_stack([(mono * S1[j][None, :]).sum(axis=1).real for j in range(4)])
    e = [float(np.abs(got - F.lorentz_straight_k1(el, pts, 0.0, w)).max()) for w in (4e-3, 1e-3)]
    assert e[1] < 3e-8 and e[0] / e[1] > 3.5
    trk = float(np.abs(K.numeric_map(el)(pts) - F.lorentz_straight_k1(el, pts, 0.0, 1e-3)).max())
    assert trk > 3e-7


@pytest.mark.parametrize('kd', ['log', 'lin'])
def test_curved_k1_body_equals_independent_planar_flow(kd):
    """K1≠0 の曲がった磁石 (tbendi の積分器。K0≠0 を含む) の本体 (水平面 y = 0、E = 0、DISFRIN = 1 でフリンジなし) が、独立な運動方程式の積分 (曲率中心のまわりの軸対称な場
    B_Y(r) = b0 + K1 ρ0 ln(r/ρ0) ('log' = K64 の暗黙の K1 の定義) または b0 + K1 (r − ρ0) ('lin' = v2 のオプション)) と 1e-9 以内。垂直方向 (Maxwell の端・K1 の場の y の依存) は未確認"""
    import k0_independent_flow as F
    old = compat.mode()
    compat.set_mode('v2')
    compat.set_bend_k1_def(kd)
    try:
        rng = np.random.default_rng(2)
        pts = np.column_stack([rng.uniform(-3e-3, 3e-3, 4), rng.uniform(-3e-3, 3e-3, 4), np.zeros(4), np.zeros(4)])
        for kw in (dict(L=0.5, ANGLE=0.1, K0=0.0, K1=0.4), dict(L=0.5, ANGLE=0.1, K0=0.03, K1=-0.3), dict(L=0.85, ANGLE=0.5236, K0=0.0, K1=-1.6782)):
            el = _mk(dict(kw, DISFRIN=1.0))
            e, _ = F.compare_body(el, pts, 2e-3, kd)
            assert e < 1e-9, (kw, e)
        # 定義を取り違えると見えること (log の追跡を lin の場と比べる)
        if kd == 'log':
            el = _mk(dict(L=0.85, ANGLE=0.5236, K0=0.0, K1=-1.6782, DISFRIN=1.0))
            assert F.compare_body(el, pts, 0.0, 'lin')[0] > 1e-7
    finally:
        compat.set_bend_k1_def(None)
        compat.set_mode(old)


@pytest.fixture(scope='module')
def mpw():
    lat = G.build(a=0.06)
    r = A.amplitude_tune(lat, D=5, nslice=8, **OPT)
    return lat, r


def test_mpw_orbit_is_closed_bump_and_series_matches_tracking(mpw):
    """3 極の挿入光源 (K0 = 0.06・−0.12・0.06 rad、軌道の振幅 1.1e-2 m): 軌道は局所のバンプ (リング先頭で 0)、ν(J, δ) の級数が pysad の実際の追跡 + NAFF と
    4e-6 以内 (K0=0 の A1 の FODO と同じ水準 1.2e-6。δ = 0・±1e-3)。挿入光源の効果が見える (K0 = 0 のリングと ν0・ADTS が違う)"""
    lat, r = mpw
    assert np.abs(r.res['orbit']).max() < 1e-12
    for (jx, jy, d) in ((1e-6, 1e-6, 0.0), (2e-6, 5e-7, 0.0), (1e-6, 1e-6, 1e-3), (1e-6, 1e-6, -1e-3)):
        t = A.tracked_tunes(lat, r, jx, jy, delta=d)
        assert abs(t[0] - r.nu(0, jx, jy, d)) < 4e-6 and abs(t[1] - r.nu(1, jx, jy, d)) < 4e-6
    r0 = A.amplitude_tune(G.build(a=0.0), D=3, nslice=8, **OPT)
    assert abs(r.nu0[1] - r0.nu0[1]) > 1e-3 and abs(r.adts['yy'] - r0.adts['yy']) > 10.0


def test_mpw_model_tracking_converges_to_series(mpw):
    """同じ模型 (磁石は追跡の素子の写像、そのほかは模型) の追跡 + NAFF と級数が 2e-8 以内 (D=5 の打ち切り。J = 1e-6、δ=0)"""
    lat, r = mpw
    m = A.model_tune(lat, 1e-6, 1e-6, 0.0, nturn=512, nfres=r, **OPT)
    assert abs(m[0] - r.nu(0, 1e-6, 1e-6)) < 2e-8 and abs(m[1] - r.nu(1, 1e-6, 1e-6)) < 2e-8


@pytest.mark.parametrize('th,tol', [(0.0, 1.5e-3), (0.1, 1.5e-3)])
def test_mpw_rdt_matches_tracked_one_turn_map(th, tol):
    """3 極の挿入光源 (ロールなし、0.1 rad の ROTATE = 結合 0.025) の共鳴駆動項 (3 次 20 項・δ 9 項、ampnf.rdt_series) が、追跡の 1 周の写像の 2 次のテイラー展開と
    1.5e-3 以内 (最大は δ² の閉軌道の h10002 の 1.3e-3、ほかは 1e-3 以下。K0=0 の F1 の水準は 6e-4)。スキューの項は θ=0 で 0"""
    lat = G.build(a=0.06, rot=th)
    res = A.one_turn_series(lat, D=3, nslice=16, bend_cubic=True, nslice_bend=8)
    h = A.rdt_series(res)
    assert h.pop('_imag') < 1e-5
    Am = res['Amat']
    rm = T.ring_map(lat.elements, 3e9, Am, amp=1e-5)
    ht = T.rdt_from_map(rm)
    tr = T.tracked_rdt_table(lat.elements, 3e9, Am)
    ht.update(tr['delta1'])
    ht['h10002'] = tr['h10002']
    worst = max(abs(h[k] - ht[k]) / max(abs(ht[k]), 1e-3) for k in h if k in ht)
    assert worst < tol, worst
    if th == 0.0:
        assert abs(h['h10101']) < 1e-9 and abs(h['h10011']) < 1e-5
    else:
        assert abs(h['h10101']) > 1e-3


def test_rotate_pi_equals_negated_k0():
    """ROTATE = π の K0≠0 の磁石 = K0 と ANGLE の符号を反転した磁石 (SAD の E1・E2 は曲げの向きに対する相対で符号そのまま)。ANGLE=0 の MPW 型と ANGLE≠0 の両方で
    模型の ν0・J¹・ξ が一致 (x → −x の鏡映)"""
    for ang, e in ((0.0, (0.0, 0.0)), (0.01, (0.5, 0.4))):
        la = G.build(a=-0.06, ang=-ang, e=e)
        lb = G.build(a=0.06, ang=ang, e=e, rot=math.pi)
        ra = A.amplitude_tune(la, D=3, nslice=8, **OPT)
        rb = A.amplitude_tune(lb, D=3, nslice=8, **OPT)
        assert abs(ra.nu0[0] - rb.nu0[0]) < 1e-8 and abs(ra.nu0[1] - rb.nu0[1]) < 1e-8
        for k in ('xx', 'xy', 'yy'):
            assert abs(ra.adts[k] - rb.adts[k]) < 1e-6 * abs(ra.adts[k]), (ang, k)
        for pl in (0, 1):
            assert abs(ra.chrom(pl, 1) - rb.chrom(pl, 1)) < 1e-5


def test_k0_to_zero_limit_matches_k0_zero_model():
    """K0 → 0 (1e-9) の磁石は K0 = 0 の模型 (曲率・端の 3 次の生成関数) と、ν0 が 2e-7、J¹ が 2e-4、ξ1 が 1e-3 (相対) 以内 (D=5。D=3 では K0=0 の模型の 4 次の項が切れて J¹ が 0.4% 違う)。差は K0=0 の模型が入れていない
    4 次以上の曲率・端の項 (K0≠0 の写像は閉じた式で全次数) と線形の端の式の違い"""
    from pysad.parser import load_sad_string
    src = G.FODO.replace('E2=0.1308996938995747)', 'E2=0.1308996938995747 K0=%s)')
    out = []
    for k0 in ('0', '1e-9'):
        lat = load_sad_string(src % k0).build_lattice('RING')
        out.append(A.amplitude_tune(lat, D=5, nslice=8, **OPT))
    a, b = out
    assert abs(a.nu0[0] - b.nu0[0]) < 2e-7 and abs(a.nu0[1] - b.nu0[1]) < 2e-7
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.adts[k] - b.adts[k]) < 2e-4 * abs(a.adts[k]), k
    for pl in (0, 1):
        assert abs(a.chrom(pl) - b.chrom(pl)) < 1e-3 * abs(a.chrom(pl))


def test_perturbation_orders_in_k0_and_tilt():
    """摂動の次数: 挿入光源による ν0・ADTS の変化は K0² (a = 0.0075 と 0.015 の比 4。a = 0.03 と 0.06 では K0⁴ の項で 10% ずれる)、傾き θ のスキューの h は θ¹ (0.05 と 0.1 の比 2)、ADTS の傾きの分は θ² (0.1 と 0.2 の比 4)"""
    r0 = A.amplitude_tune(G.build(a=0.0), D=3, nslice=8, **OPT)
    d = []
    for a in (0.0075, 0.015):
        r = A.amplitude_tune(G.build(a=a), D=3, nslice=8, **OPT)
        d.append((r.nu0[1] - r0.nu0[1], r.adts['yy'] - r0.adts['yy']))
    assert abs(d[0][0] / d[1][0] - 0.25) < 0.005 and abs(d[0][1] / d[1][1] - 0.25) < 0.01
    h = [A.rdt_series(A.one_turn_series(G.build(a=0.06, rot=th), D=3, nslice=8, bend_cubic=True, nslice_bend=8), with_delta=True) for th in (0.05, 0.1)]
    for k in ('h10101', 'h10011'):
        assert abs(h[1][k] / h[0][k] - 2.0) < 0.05, k
    base = A.amplitude_tune(G.build(a=0.06), D=3, nslice=8, **OPT).adts
    sh = []
    for th in (0.1, 0.2):
        rr = A.amplitude_tune(G.build(a=0.06, rot=th), D=3, nslice=8, **OPT).adts
        sh.append(rr['yy'] - base['yy'])
    assert 3.5 < sh[1] / sh[0] < 4.5


def test_resdrive_adts_with_tilted_mpw_includes_rotate_term():
    """傾けた (0.2 rad の ROTATE) 3 極の挿入光源の adts().total (ロールを外した lat の解析式 + ampnf の差 .rotate。以前は K0≠0 の傾いた磁石を飛ばして警告) が、
    実際の追跡 + NAFF のフィットと 3e-3 以内 (.rotate = yy −8.9 [1/(m rad)] = 1.5% を入れないと yy が 1.6% 外れる)。ロールなしの値は 1e-3 以内"""
    from pysad import resdrive as R
    lat = G.build(a=0.06, rot=0.2)
    a = R.adts(lat)
    assert a.rolled and abs(a.rotate['yy']) > 5.0 and not any('未対応' in w for w in a.warnings)
    r = A.amplitude_tune(lat, D=5, nslice=16, **OPT)
    tr = A.tracked_adts(lat, r)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.total[k] - tr[k]) / abs(tr[k]) < 3e-3, k
    no = a.total['yy'] - a.rotate['yy']
    assert abs(no - tr['yy']) / abs(tr['yy']) > 1e-2
    a0 = R.adts(G.build(a=0.06))
    r0 = A.amplitude_tune(G.build(a=0.06), D=5, nslice=16, **OPT)
    tr0 = A.tracked_adts(G.build(a=0.06), r0)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a0.total[k] - tr0[k]) / abs(tr0[k]) < 1.5e-3, k


# ---------------------------------------------------------------------------------------------
# G の続き (2026-11-22): K1≠0・F1・据え付け・場が 0・オプション・独立な運動方程式の積分
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize('name,kw', [
    ('k1_edge', dict(a=0.06, ang=0.01, e=(0.5, 0.5), k1=0.03)),
    ('steer_roll_shift', dict(a=0.06, extra='DX=1e-3 DY=-5e-4 DROTATE=0.05 ROTATE=0.1')),
    ('bend_shift_f1', dict(a=0.06, ang=0.01, extra='DX=1e-3 DY=5e-4 DROTATE=0.01 F1=0.05 FRINGE=1')),
])
def test_ring_with_k1_shift_f1_matches_tracking(name, kw):
    """K1≠0 かつ K0≠0 (tbendi の積分器の TPSA 評価)、ANGLE=0 の磁石の DX・DY・ROTATE・DROTATE (tent・texit)、ANGLE≠0 の DX・DY・DROTATE・F1 のフリンジの 3 極の挿入光源つき FODO
    (D=4、J = 3e-7): ν(J, δ) の級数が実際の追跡 + NAFF と 3e-6 (δ = 0)・1.5e-5 (δ = 1e-3) 以内。閉軌道が漏れる列 (軌道 1e-4〜8e-4 m) を含む"""
    lat = G.build(**kw)
    r = A.amplitude_tune(lat, D=4, nslice=8, **OPT)
    for d in (0.0, 1e-3):
        t = A.tracked_tunes(lat, r, 3e-7, 3e-7, delta=d)
        tol = 3e-6 if d == 0.0 else 1.5e-5            # δ ≠ 0 は νy = 0.94 が整数に近く D = 4 の δ の次数の打ち切りが効く (D = 7 で 1e-6)
        assert abs(t[0] - r.nu(0, 3e-7, 3e-7, d)) < tol and abs(t[1] - r.nu(1, 3e-7, 3e-7, d)) < tol, (name, d)


def test_zero_field_bend_is_free_flight_limit():
    """ANGLE + K0 = 0 (場が 0 の曲がった幾何の磁石): K64 の式は phib = 0 で割る (追跡は NaN) ので phib → 0 の極限 (1e-9 ANGLE) で入れる。
    独立な基準 = 場が 0 の自由な直進 (実験室の直線を、入口の面から出口の座標系へ。`k0_independent_flow.lorentz_map` の b = 0) と 3e-8 以内 (E1・E2 あり・なし、δ = 0・2e-3)"""
    import k0_independent_flow as F
    tab = A.get_tab(5, 5)
    sd, s = 1e-2, 1e-3
    rng = np.random.default_rng(1)
    pts = rng.uniform(-1e-3, 1e-3, (6, 4))
    for kw in (dict(L=0.5, ANGLE=0.1, K0=-0.1, E1=0.03, E2=0.02), dict(L=0.5, ANGLE=0.1, K0=-0.1)):
        el = _mk(kw)
        S1 = K.bend_k0_map(el)([s * tab.var(j) for j in range(4)], tab, sd)
        for delta in (0.0, 2e-3):
            u = np.column_stack([pts[:, j] / s for j in range(4)] + [np.full(6, delta / sd)])
            mono = np.prod(u[:, None, :] ** tab.mons[None, :, :], axis=2)
            got = np.column_stack([(mono * S1[j][None, :]).sum(axis=1).real for j in range(4)])
            assert np.abs(got - F.lorentz_map(el, pts, delta, w=1e-3)).max() < 3e-8


@pytest.mark.parametrize('name,kw', [('bend', dict(L=0.5, ANGLE=0.1, K0=0.03)), ('bend_edge', dict(L=0.5, ANGLE=0.1, K0=-0.02, E1=0.04, E2=0.03)),
                                     ('steer', dict(L=0.3, K0=0.06)), ('reverse', dict(L=0.5, ANGLE=-0.1, K0=0.05, E1=0.04, E2=0.03))])
def test_k64_formula_equals_independent_equations_of_motion(name, kw):
    """K64 の追跡の式 (tbendcore1・k_steer。v2 モード) が、K64 の式を使わない独立な基準 (実験室の直交座標での荷電粒子の運動方程式の数値積分。端の幅 w の
    Maxwell を満たす場、`example/14_nonlinear_dynamics/normal_form/k0_independent_flow.py`) と: 水平面 (y = 0、w = 2e-5) で 1e-11 以内 (楔・曲がる軌道の厳密な幾何)、4 次元 (Maxwell の項を含む) の差が w = 2e-3 → 1e-3 でほぼ半分になり
    (w に比例して 0 に向かう = 硬い端の極限で一致)、w = 1e-3 で 1e-7 以下"""
    import k0_independent_flow as F
    old = compat.mode()
    compat.set_mode('v2')
    try:
        el = _mk(kw)
        rng = np.random.default_rng(5)
        best = None
        for sig in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            e, _ = F.compare(el, [[1e-3, 2e-3, 0.0, 0.0]], 0.0, ws=(2e-5,), sig=sig)
            if best is None or e[0] < best[0]:
                best = (e[0], sig)
        sig = best[1]
        ph = np.column_stack([rng.uniform(-3e-3, 3e-3, 4), rng.uniform(-3e-3, 3e-3, 4), np.zeros(4), np.zeros(4)])
        eh, _ = F.compare(el, ph, 0.0, ws=(2e-5,), sig=sig)
        assert eh[0] < 1e-11
        pv = np.column_stack([rng.uniform(-1e-3, 1e-3, 3), rng.uniform(-1e-3, 1e-3, 3), rng.uniform(-1.5e-4, 1.5e-4, 3), rng.uniform(-1e-3, 1e-3, 3)])
        ev, _ = F.compare(el, pv, 0.0, ws=(2e-3, 1e-3), sig=sig)
        assert ev[1] < 1e-7 and 1.6 < ev[0] / ev[1] < 2.4
    finally:
        compat.set_mode(old)


def test_independent_flow_on_steer_wedge_k64_vs_v2():
    """独立な運動方程式の積分が見つけた K64 の ANGLE=0 の磁石 (tsteer) の端の近似: (1) k64 モードは出口の楔 E2 を入れない (B140、fb2 = 0 のとき)。入れない写像は E2 = −0.05 で
    9e-5 外れる (面の符号は v2 の線形の項で決めた sig = (−1, 1))。(2) v2 の `tsteer` は、以前は楔を薄い蹴り x tan(E)/ρ (線形) だけで入れ 3e-7 外れていたが、B221 (放射担当) が端の面を
    CALC と同じ厳密な硬い端の鎖にしたので、独立な基準と 1e-11 以内で一致する (K0 = 0.06 rad・E1 = 0.1・E2 = −0.05・L = 0.3 m。水平面、w = 2e-5)"""
    import k0_independent_flow as F
    el = _mk(dict(L=0.3, K0=0.06, E1=0.1, E2=-0.05))
    out = {}
    for mode in ('k64', 'v2'):
        old = compat.mode()
        compat.set_mode(mode)
        try:
            out[mode] = F.compare(el, [[0.0, 0.0, 0.0, 0.0], [3e-3, -2e-3, 0.0, 0.0]], 0.0, ws=(2e-5,), sig=(-1, 1))[0][0]
        finally:
            compat.set_mode(old)
    assert 5e-5 < out['k64'] < 2e-4 and out['v2'] < 1e-11


def test_options_linear_k0_magnet_and_nonlinear_difference():
    """bend_cubic・kin・fringe を全部切ると K0≠0 の磁石は入口の軌道のまわりの 1 次 (軌道 + ヤコビ行列。δ の依存は全次数) だけ (`ampnf_k0.linearize`)。
    線形の光学 (ν0) は全部入りと同じ (1e-12)、磁石自身の非線形の分が J¹ の差に出る (yy で −58 [1/(m rad)] = 10%)。RDT の `.total_sext`・`.total_bend` の分離がこれで成り立つ"""
    lat = G.build(a=0.06)
    rf = A.amplitude_tune(lat, D=4, nslice=8, bend_cubic=True, nslice_bend=4)
    rl = A.amplitude_tune(lat, D=4, nslice=8, bend_cubic=False, nslice_bend=4)
    assert abs(rf.nu0[0] - rl.nu0[0]) < 1e-12 and abs(rf.nu0[1] - rl.nu0[1]) < 1e-12
    assert 40.0 < rl.adts['yy'] - rf.adts['yy'] < 80.0
    with pytest.raises(NotImplementedError):
        A.model_tune(lat, 1e-6, 1e-6, 0.0, nturn=64, nslice=8, bend_cubic=False, nslice_bend=4)


def test_amplitude_tune_no_longer_ignores_dx_dy():
    """ampnf.amplitude_tune は素子の DX・DY を黙って無視していた (pysad 自身の不具合。H の報告)。DX・DY・DZ・DROTATE のある lat は H1 の `ampnf_coupled.amplitude_tune_ext`
    (tchge の TPSA 版、閉軌道 ≠ 0 の展開) に回る。四極 QF の DX = 2 mm: 閉軌道 3.7e-3 m、ν(J) が実際の追跡 + NAFF と 1e-5 以内 (D = 4、J = 3e-7)。無視していたら ν0 が 1e-2 動くことを見落とす"""
    from pysad.parser import load_sad_string
    src = G.FODO.replace('QUAD  QF = (L=0.3 K1=0.7)', 'QUAD  QF = (L=0.3 K1=0.7 DX=2e-3)')
    lat = load_sad_string(src).build_lattice('RING')
    r = A.amplitude_tune(lat, D=4, nslice=8, **OPT)
    assert np.abs(r.res['orbit']).max() > 1e-3
    t = A.tracked_tunes(lat, r, 3e-7, 3e-7)
    assert abs(t[0] - r.nu(0, 3e-7, 3e-7)) < 1e-5 and abs(t[1] - r.nu(1, 3e-7, 3e-7)) < 1e-5
    r0 = A.amplitude_tune(load_sad_string(G.FODO).build_lattice('RING'), D=4, nslice=8, **OPT)
    assert abs(r.nu0[0] - r0.nu0[0]) > 1e-5 or abs(r.nu0[1] - r0.nu0[1]) > 1e-5
    res = A.one_turn_series(lat, D=4, nslice=8, **OPT)
    assert any('DX' in w for w in res['warnings'])


def test_k0_magnets_with_errors_go_through_coupled_path():
    """K0≠0 の磁石に DX・DY・DROTATE がある lat (閉軌道が漏れる) は、`amplitude_tune` が H1 の結合した経路 (`ampnf_coupled`) に回り、K0≠0 の磁石は G の `tmap` (追跡と同じ tbshift・tbrot の式) で入る。
    G の内部の経路 (one_turn_series 直接) と ν0・J¹ が 1e-8 (相対) 以内で一致し、実際の追跡と ν が 2e-6 (δ = 0) 以内"""
    lat = G.build(a=0.06, ang=0.01, extra='DX=1e-3 DY=5e-4 DROTATE=0.01')
    r = A.amplitude_tune(lat, D=4, nslice=8, **OPT)
    res = A.one_turn_series(lat, D=4, nslice=8, **OPT)
    nf = A.normal_form(res)
    r2 = A.AmpTuneResult(A.tune_series(res, nf, 0), A.tune_series(res, nf, 1), (res['mux'] / (2 * math.pi), res['muy'] / (2 * math.pi)), nf, res)
    assert abs(r.nu0[0] - r2.nu0[0]) < 1e-7 and abs(r.nu0[1] - r2.nu0[1]) < 1e-7
    for k in ('xx', 'xy', 'yy'):
        assert abs(r.adts[k] - r2.adts[k]) < 1e-6 * abs(r2.adts[k]), k
    t = A.tracked_tunes(lat, r, 3e-7, 3e-7)
    assert abs(t[0] - r.nu(0, 3e-7, 3e-7)) < 2e-6 and abs(t[1] - r.nu(1, 3e-7, 3e-7)) < 2e-6


def test_tbendi_xsxkx_sign():
    """B250: tbendi の xsxkx の符号 (akxsq < 0 の枝が極限と逆) は行路 z だけに出る。θ = 5e-4 の 1 スライスの K1 = 0.3 の磁石 (L = 1 m) の z が、同じ磁石を 64 個に分けた z と、
    k64 では 8e-10 (δ = 0.01) 違い (K64 の誤り)、v2 では 1e-15 以内で一致 (連続な式)。K1 = −2 (akxsq > 0 の枝) は両モードとも 1e-15 以内"""
    from pysad.track6n import track_particles_fast

    def z_of(els, delta):
        c = np.zeros((1, 6))
        c[0, 5] = delta
        return track_particles_fast(els, c, 1, 3e9, 0.0, True, False, False)[0, 4]
    old = compat.mode()
    out = {}
    try:
        for mode in ('k64', 'v2'):
            compat.set_mode(mode)
            for K1 in (0.3, -2.0):
                one = [Bend(name='B', L=1.0, ANGLE=5e-4, K1=K1)]
                many = [Bend(name='B', L=1.0 / 64, ANGLE=5e-4 / 64, K1=K1) for _ in range(64)]
                out[(mode, K1)] = z_of(one, 1e-2) - z_of(many, 1e-2)
    finally:
        compat.set_mode(old)
    assert 5e-10 < abs(out[('k64', 0.3)]) < 2e-9 and abs(out[('v2', 0.3)]) < 1e-15
    assert abs(out[('k64', -2.0)]) < 1e-15 and abs(out[('v2', -2.0)]) < 1e-15
    # CALC (perfkern_bend.bendebody) の旗 bit2 (4 = v2 の B250): 旗なし (K64) と旗つき (v2) の出口の z の差が追跡の k64 − v2 の差と同じ (8.1e-10)
    from pysad.perfkern_bend import bendebody
    L, phi = 1.0, 5e-4
    cod = np.zeros(6)
    cod[5] = 1e-2

    def zc(fl):
        return bendebody(cod.copy(), L, phi, 0.3 * L, L / phi, L / phi, 1.0, 0.0, 0.0, 0.9999, fl)[1][4]
    assert abs((zc(0) - zc(4)) - (out[('k64', 0.3)] - out[('v2', 0.3)])) < 2e-11
    assert compat.bend_ecorr_flags() & 4 == 0 or compat.v2()


def test_delta_terms_of_h_with_leaking_closed_orbit_match_tracked_table():
    """閉軌道がリング先頭で 0 でない (磁石の外に漏れる、1.6e-4 m) リングの δ の項の h (h10002 を含む 9 項) の検証。以前は基準の `rdttrack.tracked_rdt_table` が δ=0 の閉軌道ではなく原点を
    引いていて h10002 が 6.6e3 と出ていた (展開の中心を δ=0 の閉軌道にして直した)。D = 5 の ampnf の h (3 次 20 項 + δ 9 項) が追跡の表と 1e-3 以内 (実測 1.2e-4 以下)"""
    lat = G.build(**G.KINDS['mpw_ang'])
    res = A.one_turn_series(lat, D=5, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True)
    assert np.abs(res['orbit']).max() > 1e-4
    h = A.rdt_series(res)
    h.pop('_imag')
    Am = res['Amat']
    ht = T.rdt_from_map(T.ring_map(lat.elements, 3e9, Am, amp=1e-5))
    tt = T.tracked_rdt_table(lat.elements, 3e9, Am)
    ht.update(tt['delta1'])
    ht['h10002'] = tt['h10002']
    assert abs(ht['h10002']) < 10.0
    worst = max(abs(h[k] - ht[k]) / max(abs(ht[k]), 1e-3) for k in h if k in ht)
    assert worst < 1e-3, worst


def test_astrid2_delta_difference_is_bend_slicing():
    """ASTRID2 (MPW なし、0.58 GeV) の δ = ±1e-3 の ν の級数と追跡の差 1e-5 の原因は偏向磁石 BH (0.52 rad・K1≠0) の分割数 nslice_bend (Strang 分割 O(1/n²)): 4 で 1e-5、16 で 3e-7 (νx)・3e-6 (νy)。
    級数と同じ模型の追跡の差は 1e-9 以下 (分割に依らず、追跡との差は模型の分割の誤差)"""
    from pysad.parser import load_sad_file
    lat = load_sad_file(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')).build_lattice('RING')
    err = {}
    for nb in (4, 16):
        r = A.amplitude_tune(lat, D=5, nslice=8, bend_cubic=True, kin=True, fringe=True, nslice_bend=nb)
        t = A.tracked_tunes(lat, r, 5e-9, 5e-9, delta=1e-3, momentum_eV=0.58e9)
        err[nb] = (abs(t[0] - r.nu(0, 5e-9, 5e-9, 1e-3)), abs(t[1] - r.nu(1, 5e-9, 5e-9, 1e-3)))
    assert err[4][0] > 5e-6 and err[16][0] < 1e-6 and err[16][1] < 6e-6 and err[16][1] < 0.6 * err[4][1]


@pytest.mark.parametrize('name,kw,edge,tight', [
    ('k1_noedge', dict(L=0.5, ANGLE=0.1, K0=0.0, K1=0.4), 'radial', True),
    ('k1_k0', dict(L=0.5, ANGLE=0.1, K0=0.03, K1=-0.3), 'radial', True),
    ('k1_edge_face', dict(L=0.5, ANGLE=0.1, K0=-0.02, K1=0.4, E1=0.04, E2=0.03), 'face', True),
    ('k1_edge_radial', dict(L=0.5, ANGLE=0.1, K0=-0.02, K1=0.4, E1=0.04, E2=0.03), 'radial', False)])
def test_curved_k1_with_edges_equals_independent_flow(name, kw, edge, tight):
    """K1≠0 の曲がった磁石 (K0≠0・E1・E2 を含む) の全 4 次元・端つきの独立な基準 (`k0_independent_flow.lorentz_curved_k1`: 曲率中心のまわりの軸対称な中心面の場 b(r) ('log') に面の関数 G を掛け、
    多極 (調和) 展開の端 = Lee–Whiting の前提の x³ ΔG の項を足し、Maxwell を満たす Y³ までの場。端の幅 w)。追跡 (tbendi の積分器 + Lee–Whiting + Maxwell の端 + 楔) が
    w に対して w² で 0 に向かう (w = 2e-3 → 1e-3 で減り (E なしは 4 倍、E ありの face は 1.5 倍で 2.6e-8 に近づく)、1e-3 で 5e-8 以下): E なし (K0 = 0・K0 = 0.03)、E1 = 0.04・E2 = 0.03 は v2 の 'face' (四極の場も端面から始まる模型 = 独立な基準と同じ)。
    **K64 の既定の 'radial' (四極の場は動径方向の面から始まる) は、E≠0・K1≠0・K0≠0 で 3.5e-7 に止まる** (面の取り方の違い。K64 の模型の選択で、誤りとは言えない)"""
    import warnings
    import k0_independent_flow as F
    warnings.filterwarnings('ignore')
    old = compat.mode()
    compat.set_mode('v2')
    compat.set_bend_k1_edge(edge)
    try:
        rng = np.random.default_rng(5)
        pts = np.column_stack([rng.uniform(-1e-3, 1e-3, 2), rng.uniform(-1e-3, 1e-3, 2), rng.uniform(-1.2e-4, 1.2e-4, 2), rng.uniform(-1e-3, 1e-3, 2)])
        el = _mk(kw)
        ref = K.numeric_map(el)(pts)
        e = [float(np.abs(ref - F.lorentz_curved_k1(el, pts, 0.0, w)).max()) for w in (2e-3, 1e-3)]
    finally:
        compat.set_bend_k1_edge(None)
        compat.set_mode(old)
    if tight:
        assert e[1] < 5e-8 and e[0] / e[1] > 1.4, e
    else:
        assert e[1] > 2e-7, e
