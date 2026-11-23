# -*- coding: utf-8 -*-
"""v2.4.0 H1: 結合した光学 (ROTATE した磁石) の共鳴駆動項の拡張 (`ampnf_coupled.py`): 4 次の h、δ≠0、ref≠None、DX・DY・DROTATE。
独立な基準: (1) 要素 1 個の写像 (tchge の入口・出口を TPSA にしたもの) と pysad の追跡 (v2 モード)、(2) 1 周の写像の 3 次までのテイラー展開 (`rdttrack.ring_map3`、結合したモードの正規化) から
取り出した h3・h4 (閉軌道は追跡で求める)、(3) 模型の h3・h4 から作る正規形の ADTS が ampnf の正規形 (D=5) と一致 (同じ摂動論の量)、(4) 傾けない lat で resdrive の従来の解析式 (rdt2nd)、
(5) 摂動の次数 (スキューの h4 は θ の 1 次、ADTS の変化は θ²、閉軌道は DX の 1 次)。リング: 小さな FODO (3 GeV、24 セル、6 極 K2 = 4・−6)。"""
import copy
import math
import os

import numpy as np
import pytest

from pysad import ampnf as A
from pysad import ampnf_coupled as C
from pysad import compat
from pysad import nform
from pysad import rdttrack as T
from pysad import resdrive as R
from pysad.lattice import Lattice

from pysad.tests import test_v24_F1_rotated_bend as F

OPT = dict(nslice=16, bend_cubic=True, kin=True, fringe=True, nslice_bend=4)


def _rel(h, ht, frac=1e-3):
    m = max(abs(v) for v in ht.values())
    return max(abs(h[k] - ht[k]) for k in ht if abs(ht[k]) > frac * m) / m


def _compare(lat, delta=0.0, D=None, start=None):
    rm = C.taylor_map(lat, delta=delta, D=D, start=start, **OPT)
    els = list(lat.elements)
    if start is not None:
        els = els[start:] + els[:start]
    rmt = T.ring_map3(els, 3e9, rm.A, delta=delta, amp=1e-5, amp3=2e-4)
    h3, (h4, _) = C.rdt3(rm), C.rdt4(rm)
    t3 = T.rdt_from_map(rmt)
    t4, _ = T.rdt4_from_map(rmt)
    return rm, rmt, _rel(h3, t3), _rel(h4, t4)


def _with(lat, **spec):
    """lat の素子を複製して、{種類: [(k 番目, {属性: 値})]} の属性を与える"""
    els = [copy.copy(e) for e in lat.elements]
    for tn, lst in spec.items():
        idx = [i for i, e in enumerate(els) if type(e).__name__ == tn]
        for k, attrs in lst:
            for a, v in attrs.items():
                setattr(els[idx[k]], a, v)
    return Lattice(els, name='X')


@pytest.fixture(scope='module')
def base():
    old = compat._mode
    compat.set_mode('v2')
    lat = F.rolled(0.05, nrot=1)
    out = (lat,) + _compare(lat)
    yield out
    compat.set_mode(old)


def test_element_maps_with_tchge_match_tracking_both_modes():
    """DX・DY・DROTATE・ROTATE のある BEND、DX・DY・ROTATION のある QUAD・SEXT の写像 (入口・出口の tchge を TPSA にした ops) が pysad の追跡と 2e-10 以内 (振幅 3e-4、δ = 0・1e-3)。
    v2 と k64 の両方のモードで (以前は k64 の BEND の DX が 4e-8 違った: K64 の CALC・Emittance の tchge の出口の符号は追跡 (tbshift) と逆で、模型は追跡と同じ幾何として正しい符号に固定した)"""
    old = compat._mode
    try:
        for mode in ('v2', 'k64'):
            compat.set_mode(mode)
            _element_maps_check()
    finally:
        compat.set_mode(old)


def _element_maps_check():
    from pysad.elements import Bend, Quadrupole, Sextupole
    from pysad.track6n import track_particles_fast
    D, sd, s = 5, 1e-2, 1e-3
    tab = A.get_tab(5, D)

    def build(el):
        ops = C.build_ops_ext([el], D, sd, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True)
        Phi = C.propagate(ops, tab, [s * tab.var(j) for j in range(4)], sd, D)
        return lambda x4, d: np.array([np.dot(Phi[i], np.prod(np.concatenate([np.asarray(x4) / s, [d / sd]])[None, :] ** tab.mons, axis=1)) for i in range(4)]).real

    def trk(el, x4, d):
        c = np.zeros((1, 6))
        c[0, :4] = x4
        c[0, 5] = d
        return np.array(track_particles_fast([el], c, 1, 3e9, 0.0, False, False, False))[0, :4]

    def mk(el, **kw):
        for k, v in kw.items():
            setattr(el, k, v)
        return el
    bnd = dict(name='B', L=1.2, ANGLE=0.2618, E1=0.13, E2=0.13)
    els = [mk(Bend(**bnd), DX=2e-4, DY=-1e-4), mk(Bend(**bnd, ROTATE=0.05), DROTATE=3e-3, DX=2e-4, DY=-1e-4), mk(Quadrupole(name='Q', L=0.3, K1=0.7, ROTATION=0.03), DX=2e-4, DY=-3e-4),
           mk(Sextupole(name='S', L=0.2, K2=4.0), DX=2e-4, DY=-3e-4)]
    rng = np.random.default_rng(3)
    for el in els:
        ev = build(el)
        for d in (0.0, 1e-3):
            for x in rng.normal(size=(3, 4)) * 3e-4:
                assert np.abs(ev(x, d) - trk(el, x, d)).max() < 2e-10


def test_h3_h4_match_old_path_and_tracked_map(base):
    """ROTATE を 1 台 (θ = 50 mrad、結合 9e-3、垂直分散あり): 3 次の h は従来の `ampnf.rdt_series` と 1e-12、追跡の写像のテイラー展開から取り出した h3 と 2e-5、h4 と 5e-4 (最大の項の相対。h4 の差は模型に入れていない
    K64 の四極の積分器の台形則など、docs/25 §10.5 の FODO の 1.3e-4 と同じ水準)。ADTS の対角が `resdrive.rdt2nd` の解析式と 4e-4 以内 (傾けない lat)"""
    lat, rm, rmt, e3, e4 = base
    assert rm.coupling > 5e-3 and rm.const < 1e-12 and np.abs(rm.orbit).max() < 1e-15
    r1 = A.one_turn_series(lat, D=3, nslice=16, bend_cubic=True, nslice_bend=4)
    ho = A.rdt_series(r1)
    h3 = C.rdt3(rm)
    assert max(abs(h3[k] - ho[k]) for k in h3) < 1e-10
    assert e3 < 2e-5 and e4 < 1e-3      # h4: 5e-4 → 1e-3 (v2 の追跡の四極が 4 次の合成 (B226) になり、K64 の台形則の誤差との相殺が無くなった)
    # 傾けない lat: 従来の解析式 (rdt2nd) と
    lat0 = F.rolled(0.0, nrot=1)
    rm0 = C.taylor_map(lat0, **OPT)
    h4, _ = C.rdt4(rm0)
    rl = R.resonance_terms(lat0, ['h22000', 'h11110', 'h00220', 'h31000'], nslice=16)
    for t in ('h22000', 'h11110', 'h00220', 'h31000'):
        assert abs(h4[t] - rl.total[t]) < 1e-3 * abs(rl.total[t])


def test_adts_from_h3_h4_equals_ampnf_normal_form(base):
    """模型の写像の h3・h4 から作る ADTS (h の対角 + 3 次の非共鳴の補正) が ampnf の正規形 (D=5) の J¹ と 1e-9 以内で一致 (同じ摂動論の量。結合したモードで 4 次の h が正しい確認)"""
    lat, rm, *_ = base
    a_map, _ = nform.adts_from_map(rm)
    r = A.amplitude_tune(lat, D=5, **OPT)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a_map[k] - r.adts[k]) < 1e-9 * abs(r.adts[k])


def test_h4_perturbation_order_in_roll(v2_mode):
    """摂動の次数: スキューの h4 (結合する項、θ=0 で 0) は θ の 1 次 (θ を 2 倍にすると 2.0)、ADTS の傾きの分は θ² (4 倍)"""
    out = {}
    for th in (0.025, 0.05):
        lat = F.rolled(th, nrot=1)
        rm = C.taylor_map(lat, **OPT)
        h4, _ = C.rdt4(rm)
        a, _ = nform.adts_from_map(rm)
        out[th] = (h4, a)
    rm0 = C.taylor_map(F.rolled(0.0, nrot=1), **OPT)
    a0, _ = nform.adts_from_map(rm0)
    h4_0, _ = C.rdt4(rm0)
    skew = [k for k in out[0.05][0] if abs(out[0.05][0][k]) > 0.5 and abs(h4_0[k]) < 1e-9]
    assert len(skew) > 4
    for k in skew:
        assert abs(out[0.05][0][k]) / abs(out[0.025][0][k]) == pytest.approx(2.0, rel=0.02)
    d1 = out[0.025][1]['xx'] - a0['xx']
    d2 = out[0.05][1]['xx'] - a0['xx']
    assert d2 / d1 == pytest.approx(4.0, rel=0.1)


def test_delta_nonzero_matches_tracked_map_at_delta(v2_mode):
    """δ = 1e-3 の閉軌道 (分散による 0.57 mm) のまわりの h3・h4 が、追跡の写像 (ring_map3(delta=δ)) と h3 3e-5・h4 6e-4 以内。閉軌道は 1e-9 以内。D = 5 (T3 の δ の 1 次まで)。
    δ→0 の連続性: δ = 0 の値との差が δ の 1 次 (比例)"""
    lat = F.rolled(0.05, nrot=1)
    rm, rmt, e3, e4 = _compare(lat, delta=1e-3, D=5)
    assert np.abs(rm.orbit - rmt.orbit).max() < 1e-9 and np.abs(rm.orbit).max() > 3e-4
    assert e3 < 3e-5 and e4 < 6e-4
    rm0 = C.taylor_map(lat, delta=0.0, D=5, **OPT)
    h0, (h0_4, _) = C.rdt3(rm0), C.rdt4(rm0)
    h1, (h1_4, _) = C.rdt3(rm), C.rdt4(rm)
    assert max(abs(h1_4[k] - h0_4[k]) for k in h0_4) > 1.0           # δ で h4 が有意に変わる (クロマティックな効果)


def test_ref_is_phase_reference_and_start_is_map_origin(v2_mode):
    """ref = resdrive の ref と同じ定義 (位相の基準): 傾けない lat で従来の `resonance_terms(ref=)` と 1e-4 以内 (複素数。位相込み)。|h| は ref に依らない。
    start (1 周の写像の始点) は別の量で、追跡の写像 (始点を start にした ring_map3) と一致する (h3 2e-5)"""
    lat0 = F.rolled(0.0, nrot=1)
    terms = ['h30000', 'h21000', 'h10110', 'h10200']
    for ref in (30, 57):
        rr = R.resonance_terms(lat0, terms, nslice=16, ref=ref, coupled=False)
        rc = R.resonance_terms(lat0, terms, nslice=16, ref=ref, coupled=True)
        for t in terms:
            assert abs(rr.total[t] - rc.total[t]) < 1e-4 * max(abs(rr.total[t]), 1.0)
        r0 = R.resonance_terms(lat0, terms, nslice=16, coupled=True)
        for t in terms:
            assert abs(abs(rc.total[t]) - abs(r0.total[t])) < 1e-9
    lat = F.rolled(0.05, nrot=1)
    rm, rmt, e3, e4 = _compare(lat, start=30)
    assert e3 < 2e-5 and e4 < 5e-4
    h_start = C.rdt3(rm)
    h_ring = C.rdt3(C.taylor_map(lat, **OPT))
    assert max(abs(abs(h_start[k]) - abs(h_ring[k])) for k in h_ring) > 0.1      # 始点が違えば h (非共鳴の部分) は変わる


def test_errors_dx_dy_drotate_match_tracked_map(v2_mode):
    """据え付けのずれ (四極 3 台・6 極 2 台・BEND 1 台の DX・DY) と BEND の DROTATE (3 台、閉軌道 y 2 mm): 閉軌道 1e-8 以内、h3 2e-4 以内・h4 7e-4 以内 (追跡の写像、閉軌道のまわり)。
    閉軌道は DX の 1 次 (DX を 2 倍にして比 2.0)"""
    lat = F.rolled(0.05, nrot=1)
    spec = dict(Quadrupole=[(2, dict(DX=1.5e-4, DY=-1e-4)), (7, dict(DX=-2e-4, DY=1.2e-4)), (12, dict(DX=1e-4, DY=2e-4))],
                Sextupole=[(3, dict(DX=2e-4, DY=-1e-4)), (10, dict(DX=-1e-4, DY=3e-4))], Bend=[(5, dict(DX=2e-4, DY=1e-4))])
    late = _with(lat, **spec)
    rm, rmt, e3, e4 = _compare(late)
    assert np.abs(rm.orbit).max() > 3e-4 and np.abs(rm.orbit - rmt.orbit).max() < 5e-8
    assert e3 < 1e-4 and e4 < 1e-3
    spec2 = {k: [(i, {a: 2 * v for a, v in at.items()}) for i, at in lst] for k, lst in spec.items()}
    rm2 = C.taylor_map(_with(lat, **spec2), **OPT)
    assert np.abs(rm2.orbit).max() / np.abs(rm.orbit).max() == pytest.approx(2.0, rel=0.05)
    lat2 = _with(lat, Bend=[(1, dict(DROTATE=3e-3)), (7, dict(DROTATE=-2e-3)), (13, dict(DROTATE=4e-3))])
    rmd, rmtd, e3d, e4d = _compare(lat2)
    assert np.abs(rmd.orbit).max() > 1e-3 and np.abs(rmd.orbit - rmtd.orbit).max() < 5e-8
    assert e3d < 3e-4 and e4d < 1e-3


def test_resonance_terms_api_coupled_extensions(v2_mode):
    """`resonance_terms`: ロールのある lat に 4 次の名前・δ・ref・DX/DY/DROTATE を渡すと結合した経路 (以前は NotImplementedError → 従来の式でロールを飛ばす)。δ の項 (p ≥ 1) はずれがあると NotImplementedError → 従来に戻って警告"""
    lat = F.rolled(0.05, nrot=1)
    r = R.resonance_terms(lat, ['h30000', 'h22000', 'h31000'], bend_cubic=False, bend_nslice=8)
    assert r.coupled and 'ampnf_coupled' in r.method
    rm = C.taylor_map(lat, **dict(OPT, bend_cubic=False, nslice_bend=8))
    assert abs(r.total['h22000'] - C.rdt4(rm)[0]['h22000']) < 1e-9 * abs(r.total['h22000'])
    late = _with(lat, Quadrupole=[(2, dict(DX=1.5e-4, DY=-1e-4))])
    re_ = R.resonance_terms(late, ['h30000', 'h22000'], bend_cubic=False)
    assert re_.coupled and any('閉軌道' in w for w in re_.warnings) and np.abs(re_.orbit).max() > 1e-5
    re2 = R.resonance_terms(late, ['h11001'], bend_cubic=False)         # p ≥ 1 は結合した経路で扱えない → 従来の式 + 警告
    assert not getattr(re2, 'coupled', False) and any('未対応' in w for w in re2.warnings)


def test_adts_rotate_term_follows_delta_series(v2_mode):
    """`adts(delta≠0)` の傾きの項 .rotate は、模型の ν の級数 (D=5) の J¹ の δ 依存を δ で評価 (以前は δ=0 の値で警告)。全 24 台を 0.1 rad 傾けた FODO、δ = 1e-3:
    δ=0 の値と 8% 違い (入れる意味がある)、δ の閉軌道のまわりの模型の写像 (h3・h4 の正規形、D=5) の傾きの差と 2e-3 (.rotate の相対) 以内"""
    lat = F.rolled(0.1, nrot=24)
    lat0 = R.unrolled_lattice(lat)
    r0 = R.adts(lat)
    r1 = R.adts(lat, delta=1e-3)
    a = {}
    for nm, ll in (('r', lat), ('u', lat0)):
        rm = C.taylor_map(ll, delta=1e-3, D=5, **OPT)
        a[nm], _ = nform.adts_from_map(rm)
    for k in ('xx', 'xy', 'yy'):
        want = a['r'][k] - a['u'][k]
        assert abs(r1.rotate[k] - want) < 3e-3 * abs(want)
        assert abs(r1.rotate[k] - r0.rotate[k]) > 0.04 * abs(r0.rotate[k])
    assert not any('δ=0 の値' in w for w in r1.warnings)


def test_h4_difference_to_tracking_is_the_k64_quadrature(k64_mode):
    """h4 の追跡との差 (3e-4〜6e-4) の原因の切り分け: 厚い 6 極を連続に近く (nslice=64) して四極の運動学を K64 の積分器の台形則 (quad_rule='k64') にすると差が 1.5e-4 以下 (7e-5) に下がり、
    厳密な積分 (quad_rule='exact') のままだと 5e-4 以上残る。K64 の厚い 6 極の積分器の分 (docs/25 §10.5) が残りの 7e-5。傾けない FODO。k64 の追跡で見る (v2 の追跡の四極は 4 次の合成で台形則ではない。B226)"""
    lat = F.rolled(0.0, nrot=1)
    kw = dict(nslice=64, bend_cubic=True, kin=True, fringe=True, nslice_bend=8)
    rm_k = C.taylor_map(lat, D=4, quad_rule='k64', **kw)
    rm_e = C.taylor_map(lat, D=4, quad_rule='exact', **kw)
    rmt = T.ring_map3(list(lat.elements), 3e9, rm_k.A, amp=1e-5, amp3=2e-4)
    t4, _ = T.rdt4_from_map(rmt)
    r_k = _rel(C.rdt4(rm_k)[0], t4)
    r_e = _rel(C.rdt4(rm_e)[0], t4)
    assert r_k < 1.5e-4 and r_e > 4e-4


def test_adts_with_closed_orbit_matches_tracking_nafffit_with_constant(v2_mode):
    """閉軌道 ≠ 0 (四極 3・6 極 2 の DX・DY、BEND 2 台の DROTATE、閉軌道 1.6 mm・結合 0.06) の ADTS: `amplitude_tune_ext` の J¹ が、追跡 + NAFF の定数項つきフィット (`tracked_adts_fit`、
    J = 1e-7・3e-8) と 3e-3 以内。定数項なし (`ampnf.tracked_adts`) は J = 3e-8 で 2〜4% 外れる (模型と追跡の線形のチューンの差 5e-7〜1e-6 が入る)"""
    lat = F.rolled(0.05, nrot=1)
    late = _with(lat, Quadrupole=[(2, dict(DX=1.5e-4, DY=-1e-4)), (7, dict(DX=-2e-4, DY=1.2e-4)), (12, dict(DX=1e-4, DY=2e-4))],
                 Sextupole=[(3, dict(DX=2e-4, DY=-1e-4)), (10, dict(DX=-1e-4, DY=3e-4))], Bend=[(1, dict(DROTATE=3e-3)), (7, dict(DROTATE=-2e-3))])
    r = C.amplitude_tune_ext(late, D=4, nslice=16, nslice_bend=4)
    assert np.abs(r.res['orbit']).max() > 1e-3
    for j in (1e-7, 3e-8):
        t, off = C.tracked_adts_fit(late, r, j=j)
        for k in ('xx', 'xy', 'yy'):
            assert abs(t[k] / r.adts[k] - 1) < 3e-3
        assert abs(off['x']) < 3e-6 and abs(off['y']) < 3e-6
    t0 = A.tracked_adts(late, r, j=3e-8)
    assert max(abs(t0[k] / r.adts[k] - 1) for k in ('xx', 'xy', 'yy')) > 1e-2
