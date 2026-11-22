# -*- coding: utf-8 -*-
"""v2.4.0 F1: ROTATE (ロール角) した偏向磁石・四極を、結合した光学 (4×4 の固有モード、垂直分散) で RDT・ADTS (resdrive) と ν(J, δ) の級数 (ampnf) に入れる。
傾いた磁石の写像は R⁻¹ ∘ (傾けない磁石の写像) ∘ R (全次数の θ)。独立な基準: (1) 同じ模型の追跡 + NAFF、(2) pysad の実際の追跡 + NAFF (全素子・全非線形)、
(3) 追跡の 1 周の写像の 2 次のテイラー展開から取り出した h (rdttrack、結合したモードの正規化)、(4) ROTATE = π の偏向磁石 = 逆曲げ (ANGLE を負にし、SAD の E1・E2 は曲げの向きに対する相対なので符号そのまま)、
(5) 摂動の次数 (ADTS は θ の偶関数なので傾き効果は θ²、スキューの h は θ の 1 次)。リング: 小さな FODO (3 GeV、24 セル、6 極 K2 = 4・−6)、k64 (模型は K1 = 0 の BEND で
k64 と v2 に差がない)。"""
import math
import os

import numpy as np
import pytest

from pysad import ampnf as A
from pysad import rdttrack as T
from pysad import resdrive as R
from pysad.parser import load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
ANG, EDG = 0.2617993877991494, 0.1308996938995747


def ring(bend_extra, nrot=1, quad_extra=''):
    """FODO の先頭の nrot セルを、BR (bend_extra の BEND) を使うセルにする"""
    src = FODO.replace('BEND  B  = (', 'BEND BR = (%s)\nBEND  B  = (' % bend_extra)
    src = src.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE CELLR = (SF DS QF D2 BR D2 QD DS SD D1);\nLINE RING = (ORIGIN RF %d*CELLR %d*CELL);' % (nrot, 24 - nrot))
    if quad_extra:
        src = src.replace('QUAD  QF = (L=0.3 K1=0.7)', 'QUAD  QF = (L=0.3 K1=0.7 %s)' % quad_extra)
    return load_sad_string(src).build_lattice('RING')


def rolled(theta, nrot=1):
    return ring('ROTATE=%.16g L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g' % (theta, ANG, EDG, EDG), nrot)


OPT = dict(bend_cubic=True, kin=True, fringe=True, nslice_bend=4)


def test_linear_modes_uncoupled_is_twiss_and_coupled_is_symplectic():
    """結合がなければ Twiss の正規化と同じ (Am の逆が正規化行列)。結合があれば Ainv A Am = 回転 (固有値 e^{±iμ})・ē^T S e = 2i の規格化で Am は シンプレクティック"""
    lat = rolled(0.0)
    r0 = A.one_turn_series(lat, D=2)
    assert r0['coupling'] == 0.0 and r0['twiss0'][1] > 0
    lat = rolled(0.1, nrot=24)
    r = A.one_turn_series(lat, D=3, bend_cubic=True)
    assert r['coupling'] > 1e-2 and r['twiss0'][0] is None
    Am = r['Amat']
    S4 = np.zeros((4, 4))
    S4[0, 1], S4[1, 0], S4[2, 3], S4[3, 2] = 1, -1, 1, -1
    assert np.abs(Am.T @ S4 @ Am - S4).max() < 1e-10             # 正規化座標は正準 (シンプレクティック)
    # 1 周の線形写像が正規化座標で 2 つの回転 (μx, μy) になる
    ops, _ = A.build_ops(lat, 0, 1.0, 8, delta=0.0, bend_cubic=True, nslice_bend=4, Dnl=3)
    Mt = A.NSeg()
    for o in ops:
        if o[0] == 'lin':
            Mt = Mt.then(o[1])
    Rn = r['Ainv'] @ Mt.M @ Am
    for k, mu in enumerate((r['mux'], r['muy'])):
        blk = Rn[2 * k:2 * k + 2, 2 * k:2 * k + 2]
        want = np.array([[math.cos(mu), math.sin(mu)], [-math.sin(mu), math.cos(mu)]])
        assert np.abs(blk - want).max() < 1e-10
    assert np.abs(Rn[:2, 2:]).max() < 1e-10 and np.abs(Rn[2:, :2]).max() < 1e-10


@pytest.fixture(scope='module')
def single():
    lat = rolled(0.05)
    r = A.amplitude_tune(lat, D=5, nslice=16, **OPT)
    return lat, r


def test_series_with_rolled_bend_matches_model_and_real_tracking(single):
    """1 台を 50 mrad 傾けた FODO (結合 |A12| = 9e-3、垂直分散あり): 級数と、同じ模型の追跡 + NAFF が 3e-8 以内 (D=5 の打ち切り)、
    pysad の実際の追跡 + NAFF が 3e-6 以内 (ロールなしの A1 の FODO と同じ水準 1.2e-6)"""
    lat, r = single
    assert r.res['coupling'] > 5e-3
    (jx, jy) = (1e-6, 1e-6)
    f = A.model_tune(lat, jx, jy, 0.0, nturn=1024, nslice=16, nfres=r, **OPT)
    s = (r.nu(0, jx, jy, 0.0), r.nu(1, jx, jy, 0.0))
    assert max(abs(f[0] - s[0]), abs(f[1] - s[1])) < 3e-8
    t = A.tracked_tunes(lat, r, jx, jy)
    assert max(abs(t[0] - s[0]), abs(t[1] - s[1])) < 3e-6


def _tracked_h(lat, res):
    Am = res['Amat']
    rm = T.ring_map(lat.elements, 3e9, Am, amp=1e-5)
    h = T.rdt_from_map(rm)
    tr = T.tracked_rdt_table(lat.elements, 3e9, Am)
    h.update(tr['delta1'])
    h['h10002'] = tr['h10002']
    return h


def test_rdt_with_rolled_bend_matches_tracked_one_turn_map():
    """結合した光学の共鳴駆動項 (3 次 p=0 の 20 項・δ の 9 項、ampnf.rdt_series) が、追跡の 1 周の写像の 2 次のテイラー展開 (結合したモードの正規化で取り出す) と 6e-4 以内
    (ロールなしと同じ水準 3.3e-4 = 厚い 6 極の分割と写像の取り出し)。スキューの項 (h10101・h10011 ほか) は 1e-5 以内"""
    lat = rolled(0.05)
    res = A.one_turn_series(lat, D=3, nslice=16, bend_cubic=True, nslice_bend=16)
    h = A.rdt_series(res)
    assert h.pop('_imag') < 1e-9
    ht = _tracked_h(lat, res)
    skew = ('h10101', 'h10011', 'h02100', 'h20010', 'h11010', 'h00030', 'h00120')
    for k in h:
        if k in ht:
            assert abs(h[k] - ht[k]) < 6e-4, k
    for k in skew:
        assert abs(h[k]) > 1e-3 or k in ('h02100', 'h20010', 'h11010', 'h00030', 'h00120')
        assert abs(h[k] - ht[k]) < 1e-5, k


def test_skew_rdt_is_first_order_in_theta_and_adts_shift_is_second_order():
    """スキューの h は θ の 1 次 (θ = 0.025 と 0.05 の比が 2)、ADTS の傾き効果 (.rotate) は θ の 2 次 (θ = 0.05 と 0.1 の比が 4、全 24 台を傾ける)。
    ADTS は θ → −θ・y → −y の鏡映で不変なので 1 次の項がない"""
    h = []
    for th in (0.025, 0.05):
        res = A.one_turn_series(rolled(th), D=3, nslice=8, bend_cubic=True, nslice_bend=8)
        h.append(A.rdt_series(res, with_delta=True))
    for k in ('h10101', 'h10011', 'h20010'):
        assert abs(h[1][k] / h[0][k] - 2.0) < 0.02, k
    # 偶関数の項 (h30000) は θ² の変化
    d1 = abs(h[0]['h30000'] - A.rdt_series(A.one_turn_series(rolled(0.0), D=3, nslice=8, bend_cubic=True, nslice_bend=8), with_delta=False)['h30000'])
    d2 = abs(h[1]['h30000'] - A.rdt_series(A.one_turn_series(rolled(0.0), D=3, nslice=8, bend_cubic=True, nslice_bend=8), with_delta=False)['h30000'])
    assert 3.5 < d2 / d1 < 4.5
    rot = []
    for th in (0.05, 0.1):
        a = R.adts(rolled(th, nrot=24))
        rot.append(a.rotate)
    for k in ('xx', 'xy', 'yy'):
        assert 3.7 < rot[1][k] / rot[0][k] < 4.3, k


@pytest.fixture(scope='module')
def allrot():
    lat = rolled(0.1, nrot=24)
    r = A.amplitude_tune(lat, D=5, nslice=16, **OPT)
    return lat, r


def test_adts_with_rolled_bends_matches_ampnf_and_tracking(allrot):
    """全 24 台を 0.1 rad 傾けた FODO (結合 0.03): adts().total (= ロールなしの解析式 + .rotate) が、結合した光学の模型の J¹ (ampnf D=5 の全部入り) と 5e-4、
    pysad の実際の追跡 + NAFF のフィット (正規形の J・φ から出発、J² まで) と 1.5e-3 以内。.rotate を足さない値 (従来) は追跡と 0.1〜0.3% 違う (xx で 2.7e-3)"""
    lat, r = allrot
    a = R.adts(lat)
    assert a.rolled and 0.5 < abs(a.rotate['xx']) < 3.0 and any('rotate' in w or '結合した光学' in w for w in a.warnings)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.total[k] - r.adts[k]) / abs(r.adts[k]) < 5e-4, k
    tr = A.tracked_adts(lat, r)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.total[k] - tr[k]) / abs(tr[k]) < 1.5e-3, k
    no = a.total['xx'] - a.rotate['xx']
    assert abs(no - tr['xx']) / abs(tr['xx']) > 2e-3 and abs(a.total['xx'] - tr['xx']) < abs(no - tr['xx'])


def test_reverse_bend_equals_rotation_by_pi():
    """逆曲げ: ROTATE = π の偏向磁石 (ANGLE = +a) は ANGLE = −a の偏向磁石 (SAD の E1・E2 は曲げの向きに対する相対で、符号は同じ) と同じ写像 (x → −x の鏡映)。
    模型の級数 (ν0・J¹・ξ) は 1e-9 で一致、adts() の合計 (ロールの項つき) は ANGLE < 0 の解析式 (ロールなし) と 1e-3 以内 (これが θ = π の全次数の検証。ロールを入れないと yy が 1.6% 違う)"""
    la = ring('L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g' % (-ANG, EDG, EDG), nrot=2)
    lb = ring('L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g ROTATE=%.16g' % (ANG, EDG, EDG, math.pi), nrot=2)
    ra = A.amplitude_tune(la, D=5, nslice=8, **OPT)
    rb = A.amplitude_tune(lb, D=5, nslice=8, **OPT)
    assert abs(ra.nu0[0] - rb.nu0[0]) < 1e-9 and abs(ra.nu0[1] - rb.nu0[1]) < 1e-9
    for k in ('xx', 'xy', 'yy'):
        assert abs(ra.adts[k] - rb.adts[k]) < 1e-8 * abs(ra.adts[k]), k
    for pl in (0, 1):
        assert abs(ra.chrom(pl, 1) - rb.chrom(pl, 1)) < 1e-9
    aa = R.adts(la)
    ab = R.adts(lb)
    assert not hasattr(aa, 'rotate') or not getattr(aa, 'rolled', None)
    for k in ('xx', 'xy', 'yy'):
        assert abs(ab.total[k] - aa.total[k]) / abs(aa.total[k]) < 1e-3, k
    assert abs(ab.rotate['yy']) / abs(aa.total['yy']) > 0.01
    # RDT も同じ
    h1 = R.resonance_terms(la, ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')).total
    h2 = R.resonance_terms(lb, ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')).total
    for k in h1:
        assert abs(h1[k] - h2[k]) < 3e-4 * max(abs(h1[k]), 1.0), k


def test_rotated_sextupole_and_quad_in_ampnf():
    """ロール角つきの 6 極 (パーサは ROTATE を ROTATION に入れる) が ampnf の級数に入る (以前は黙って無視)。J¹ が adts().total と 2e-3 以内。傾けた四極 (結合) は追跡と 3e-6 以内"""
    src = FODO.replace('SF = (L=0.2 K2=4) SD', 'SF = (L=0.2 K2=4 ROTATE=0.3) SD')
    assert 'ROTATE=0.3' in src
    lat = load_sad_string(src).build_lattice('RING')
    r = A.amplitude_tune(lat, D=5, nslice=16, **OPT)
    a = R.adts(lat)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.total[k] - r.adts[k]) / abs(r.adts[k]) < 2e-3, k
    base = R.adts(load_sad_string(FODO).build_lattice('RING'))
    assert abs(a.total['xx'] / base.total['xx'] - 1) > 0.02              # 傾けると ADTS が変わる (無視していたら変わらない)
    latq = ring('L=1.2 ANGLE=%.16g E1=%.16g E2=%.16g' % (ANG, EDG, EDG), nrot=1, quad_extra='ROTATE=0.05')
    assert R.rolled_elements(latq) == ['QF'] or 'QF' in R.rolled_elements(latq)
    rq = A.amplitude_tune(latq, D=5, nslice=8, **OPT)
    assert rq.res['coupling'] > 1e-3
    t = A.tracked_tunes(latq, rq, 1e-6, 1e-6)
    s = (rq.nu(0, 1e-6, 1e-6, 0.0), rq.nu(1, 1e-6, 1e-6, 0.0))
    assert max(abs(t[0] - s[0]), abs(t[1] - s[1])) < 3e-6
