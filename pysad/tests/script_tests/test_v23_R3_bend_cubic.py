# -*- coding: utf-8 -*-
"""v2.3.0 R3: 偏向磁石の 3 次の項 (曲率 (h/2) x (px²+py²)・硬い端 ±(h/2) px y²・複合機能の本体・E1/E2≠0 の端) を入れた解析の ADTS
(resdrive.adts(bend_cubic=True)、pysad/bendnl.py) を、追跡 (k64、conftest の既定) の直接のフィット (adtsdirect.direct_adts) と比べる。
リング: 小さな FODO (3 GeV、105.6 m、24 セル、SF・SD の K2L = 4, −6、E1=E2≠0 の偏向磁石)、ASTRID2 (0.58 GeV、45.7 m、30° の複合機能偏向磁石)、
大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)。条件: RF 空洞止め、放射なし、δ=0、結合なし、NAFF Hann p=1、512 周 (保存した追跡は 1024 周)。"""
import copy
import dataclasses
import json
import os

import numpy as np
import pytest

from pysad import adtsdirect as D
from pysad import bendnl as B
from pysad import resdrive as R
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
KEKLS = os.path.join(HERE, '..', 'kek-ls_v03_12.lat')
EX = os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'adts_rdt')


def fodo(k_sf=4.0, k_sd=-6.0):
    return load_sad_string(FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))).build_lattice('RING')


def find(lat, name):
    return [e for e in lat.elements if e.name == name][0]


def rel(a, b):
    return abs(a - b) / abs(b)


def test_taylor2_compose_and_inverse():
    rng = np.random.default_rng(3)
    R_ = np.eye(4) + 0.1 * rng.normal(size=(4, 4))
    T = rng.normal(size=(4, 4, 4))
    T = 0.5 * (T + T.transpose(0, 2, 1))
    m = B.Taylor2(R_, T)
    mi = m.inverse()
    z = 1e-4 * rng.normal(size=4)
    f = lambda mm, z_: mm.R @ z_ + np.einsum('iab,a,b->i', mm.T, z_, z_)
    # 2 次まで: m⁻¹(m(z)) = z + O(z³)
    assert np.abs(f(mi, f(m, z)) - z).max() < 1e-9
    c = m.then(mi)
    assert np.abs(c.R - np.eye(4)).max() < 1e-12 and np.abs(c.T).max() < 1e-12


def test_sector_bend_generator_is_curvature_and_hard_edge_terms():
    """E1=E2=0 の偏向磁石 (FODO の B): スライスの生成関数 = (h/2)·ds·(x px² + x py²)、端 = ±(h/2) px y² (R2 が追跡の奇数次から推定した 2 項)。他の係数は 1e-5·h 未満"""
    b = dataclasses.replace(find(fodo(), 'B'), E1=0.0, E2=0.0)
    h = b.ANGLE / b.L
    d = B.bend_decomposition(b, 16)
    ds = d['ds']
    fb = d['f_b']
    assert fb[(1, 2, 0, 0)] / (0.5 * h * ds) == pytest.approx(1.0, abs=2e-4)
    assert fb[(1, 0, 0, 2)] / (0.5 * h * ds) == pytest.approx(1.0, abs=2e-4)
    for k, v in fb.items():
        if k not in ((1, 2, 0, 0), (1, 0, 0, 2)):
            assert abs(v) < 1e-5 * h * ds, (k, v)
    assert d['f_e1'][(0, 1, 2, 0)] / (0.5 * h) == pytest.approx(1.0, abs=1e-9)
    assert d['f_e2'][(0, 1, 2, 0)] / (0.5 * h) == pytest.approx(-1.0, abs=1e-9)
    assert all(abs(v) < 1e-9 * h for k, v in d['f_e1'].items() if k != (0, 1, 2, 0))
    assert d['R_b_check'] < 1e-12


def test_combined_function_body_cubic_is_k1_h_x3_over_6():
    """K1≠0 の偏向磁石 (ASTRID2 の BH、30°、K1=-1.68): 本体の 3 次の項 = (h/2) x (px²+py²) + (K1 h/6) x³ (K64 の bendecorr の kick Δpx = −K1 h x²/2 と一致)。
    Maxwell の方程式を満たす項 (K1 h x³/3 − K1 h x y²/2) とは違う (K64 の定義。docs/25 §8 の注)"""
    b = find(load_sad_file(ASTRID).build_lattice('RING'), 'BH')
    h = b.ANGLE / b.L
    d = B.bend_decomposition(b, 16)
    ds = d['ds']
    assert d['f_b'][(3, 0, 0, 0)] / (b.K1 * h / 6 * ds) == pytest.approx(1.0, abs=2e-3)
    assert d['f_b'][(1, 2, 0, 0)] / (0.5 * h * ds) == pytest.approx(1.0, abs=2e-3)
    assert abs(d['f_b'].get((1, 0, 2, 0), 0.0)) < 2e-3 * abs(b.K1 * h / 6 * ds)


def test_sub_bend_chain_equals_whole_bend():
    """偏向磁石を n 個の部分偏向磁石に分けた連続は 1 つの偏向磁石と一致する (端の硬いフリンジが隣り合うスライスで打ち消す。分解の前提)"""
    from pysad.track6n import track_particles_fast
    rng = np.random.default_rng(1)
    Z = np.zeros((6, 6))
    Z[:, :4] = rng.normal(size=(6, 4)) * 1e-3
    for lat, name, p in ((fodo(), 'B', 3e9), (load_sad_file(ASTRID).build_lattice('RING'), 'BH', 0.58e9)):
        b = find(lat, name)
        ref = track_particles_fast([b], Z.copy(), 1, p, 0.0, False, False, False)
        n = 5
        subs = [dataclasses.replace(b, L=b.L / n, ANGLE=b.ANGLE / n, E1=b.E1 if k == 0 else 0.0, E2=b.E2 if k == n - 1 else 0.0) for k in range(n)]
        out = track_particles_fast(subs, Z.copy(), 1, p, 0.0, False, False, False)
        assert np.abs(np.array(out) - np.array(ref))[:, :4].max() < 1e-10


def test_w_propagation_matches_twiss_at_exit():
    """分解した線形写像 (E1 の端・スライス・E2 の端) で入口の w を出口まで伝播すると、Twiss の出口の値と一致する"""
    for lat, name in ((fodo(), 'B'), (load_sad_file(ASTRID).build_lattice('RING'), 'BH')):
        rows = lat.twiss_table()
        i = [k for k, e in enumerate(lat.elements) if e.name == name][0]
        a, info = B.bend_atoms(lat.elements[i], rows[i], rows[i + 1], R._monos_of_degree(3), nslice=8)
        assert info['w_err'] < 1e-8 and a.shape[1] == 10


def test_fodo_adts_with_bend_cubic_matches_tracking():
    """FODO (K2L = 4, −6): 偏向磁石の 3 次の項を入れた合計が追跡の直接のフィットと 2e-3 以内 (R2 の合計は xx 2.7%・xy 1.6%・yy 9.5% 外れていた)"""
    lat = fodo()
    A0 = R.adts(lat, bend_cubic=False)
    A1 = R.adts(lat, bend_cubic=True)
    r = D.direct_adts(lat, P0, 1e-7, order=2, nturn=512, nstep=3)
    for k, key in (('xx', 'xx'), ('xy', 'xy_mean'), ('yy', 'yy')):
        assert rel(A1.total[k], r[key]) < 2e-3, (k, A1.total[k], r[key])
    assert rel(A0.total['yy'], r['yy']) > 5e-2 and rel(A0.total['xx'], r['xx']) > 1e-2
    # 保存した 1024 周・70 粒子の追跡 (adts_bend_scan_fodo.pysad.json)
    f = os.path.join(EX, 'adts_bend_scan_fodo.pysad.json')
    if os.path.exists(f):
        d = [x for x in json.load(open(f))['rows'] if x['lam'] == 1.0][0]
        for k in ('xx', 'xy', 'yy'):
            assert rel(A1.total[k], d['track'][k]) < 1e-3


def test_cross_term_is_linear_in_sextupole_strength_and_bend_part_is_independent():
    A1, A2, A0 = R.adts(fodo(4.0, -6.0)), R.adts(fodo(8.0, -12.0)), R.adts(fodo(0.0, 0.0))
    for k in ('xx', 'xy', 'yy'):
        assert A2.cross[k] == pytest.approx(2.0 * A1.cross[k], rel=1e-9, abs=1e-9)
        assert A2.bend[k] == pytest.approx(A1.bend[k], rel=1e-9, abs=1e-12)
        assert abs(A0.cross[k]) < 1e-9
        assert A2.sext[k] == pytest.approx(4.0 * A1.sext[k], rel=1e-9)        # 6 極 2 次は強さの 2 乗
    assert A1.bend_included and A1.bend_natoms == 24 * 18


def test_sextupole_subset_of_merged_sum_equals_adts_sext():
    """偏向磁石の原子と並べて 2 次の和を作る経路で、6 極だけの部分が従来の adts().sext と一致する (並べ替えの検算)"""
    for lat in (fodo(), load_sad_file(ASTRID).build_lattice('RING')):
        A = R.adts(lat)
        for k in ('xx', 'xy', 'yy'):
            assert A.bend_check[k] == pytest.approx(A.sext[k], rel=1e-10)


def test_nslice_convergence_second_order():
    """偏向磁石の本体の分割数 nslice への収束 (O(1/nslice²))。FODO・ASTRID2 で nslice=4 と 64 の合計の差が 1e-4 以内、8 と 64 は 4 以内の約 1/4"""
    for lat in (fodo(), load_sad_file(ASTRID).build_lattice('RING')):
        rows = lat.twiss_table()
        t = {n: R.adts(lat, rows=rows, bend_cubic=True, bend_nslice=n).total for n in (4, 8, 64)}
        for k in ('xx', 'xy', 'yy'):
            e4, e8 = abs(t[4][k] - t[64][k]), abs(t[8][k] - t[64][k])
            assert e4 < 1e-4 * abs(t[64][k])
            assert e8 < 0.35 * e4 + 1e-12


def test_astrid2_and_large_ring_a_vs_saved_tracking():
    """ASTRID2 (30° の複合機能偏向磁石) と大きなリング A (kek-ls_v03_12.lat) の、保存した追跡 (adts_bend_scan_*.pysad.json、1024 周・70 粒子・6 極 ×1) との差:
    ASTRID2 は 4e-3 以内 (R2: 11%・15%・6.5%)、大きなリング A は xx・xy が 1e-3、yy が 5e-3 以内 (R2: 0.2%・1.2%・3.6%)"""
    for path, line, name, tol in ((ASTRID, 'RING', 'astrid', (4e-3, 4e-3, 4e-3)), (KEKLS, 'RING', 'kekls', (1e-3, 1e-3, 5e-3))):
        f = os.path.join(EX, 'adts_bend_scan_%s.pysad.json' % name)
        if not os.path.exists(f):
            pytest.skip('保存した追跡がない')
        d = [x for x in json.load(open(f))['rows'] if x['lam'] == 1.0][0]
        lat = load_sad_file(path).build_lattice(line)
        A = R.adts(lat)
        for k, t in zip(('xx', 'xy', 'yy'), tol):
            assert rel(A.total[k], d['track'][k]) < t, (name, k, A.total[k], d['track'][k])


def test_disfrin_changes_hard_edge_of_k1_bend():
    """K1≠0 の偏向磁石は DISFRIN=1 で硬い端のフリンジ (Maxwell の項) も切れる (K64 の tbendi)。分解がそれを追う: ASTRID2 の端の生成関数が 0 になり、
    DISFRIN=1 の追跡の直接の合計と 6 極 ×1 で 1e-2 以内 (xx・xy・yy。J=2e-8)"""
    lat = load_sad_file(ASTRID).build_lattice('RING')
    l1 = copy.deepcopy(lat)
    for e in l1.elements:
        if type(e).__name__ in ('Quadrupole', 'Bend'):
            e.DISFRIN = 1.0
    d = B.bend_decomposition(find(l1, 'BH'), 8)
    assert max([abs(v) for v in d['f_e1'].values()] + [0.0]) < 1e-9
    A = R.adts(l1, bend_cubic=True)
    assert A.fringe['xx'] == 0.0
    r = D.direct_adts(l1, 0.58e9, 2e-8, order=2, nturn=512, nstep=3)
    for k, key in (('xx', 'xx'), ('xy', 'xy_mean'), ('yy', 'yy')):
        assert rel(A.total[k], r[key]) < 1e-2, (k, A.total[k], r[key])


def test_ffs_amplitude_tune_shift_bend_option():
    import io
    from pysad import ffscmd, compat
    old = compat._mode
    compat.set_mode('v2')
    try:
        src = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)') + '''
FFS;
USE RING; CELL; CALC;
Print["ON ", AmplitudeTuneShift[BendCubic->True]];
Print["OFF ", AmplitudeTuneShift[BendCubic->False]];
Exit;
'''
        buf = io.StringIO()
        ffscmd.run_sad_string(src, out=buf)
    finally:
        compat.set_mode(old)
    txt = buf.getvalue().replace(chr(10), '')
    on = txt[txt.index('ON {'):txt.index('OFF {')]
    off = txt[txt.index('OFF {'):]
    assert 'SextBendCross->{10.9345' in on and 'BendCubic->{-.0291' in on and 'Total->{408.843' in on
    assert 'BendCubic->{0,0,0}' in off and 'SextBendCross->{0,0,0}' in off and 'Total->{397.937' in off
