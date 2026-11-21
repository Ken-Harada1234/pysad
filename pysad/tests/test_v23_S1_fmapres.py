# -*- coding: utf-8 -*-
"""v2.3.0 S1: 周波数マップと RDT・ADTS の接続 (pysad/fmapres.py)、誤差つき・δ つき・k64・MPW つきの周波数マップ、gs との比較 (KEK-LS)。
追跡の重い計算は小さな格子にする。保存済みの出力 (example/v2_3/S1_*) を読む試験は、その出力が無ければ skip する。"""
import math
import os

import numpy as np
import pytest

from pysad import compat, fmapres as FR
from pysad.freqmap import FMapResult, frequency_map
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
EX = os.path.join(ROOT, 'example', 'v2_3')
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
ASTRID_ID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
FODO = os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad')


@pytest.fixture(autouse=True)
def _v2():
    old = compat.mode()
    compat.set_mode('v2')
    yield
    compat.set_mode(old)


def _src(path):
    return open(path, encoding='utf-8').read()


def _fodo_src():
    return _src(FODO).replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)').replace('QF = (L=0.3 K1=0.7)', 'QF = (L=0.3 K1=0.697)')


# ---------------------------------------------------------------------------------------------------------------------
def test_pendulum_criterion_vs_map():
    """共鳴の島に捕らえられる条件 (振り子の式、位相つき) を、同じ生成関数の 1 周の写像 (J, φ) の追跡と比べる: 3νx の項 V = 2 Re[h (2J)^(3/2) e^{3iφ}]、
    ν(J) = ν0 + aJ。予想と追跡 (最後の 3φ − 2πpn の位相の変化が 2π 未満 = 位相が巡らない = 捕らえられる) が 90% 以上一致する"""
    h = 0.1 * np.exp(0.7j)
    a = 2000.0
    nu0 = 20.0 / 3.0 - 0.012
    pr = FR.ResonancePredictor({(3, 0): [(3, 0, 0, 0, h)]}, [[a, 0.0], [0.0, 0.0]], (nu0, 4.9), (0.0, 1.0, 0.0, 1.0, 0.0, 0.0))
    jres = 0.012 / a                    # Δ0 = 0 になる J (3(ε + aJ) = 0)
    J0 = np.linspace(0.3 * jres, 2.2 * jres, 160)
    x = np.sqrt(2 * J0)
    P = pr.predict(x, np.full_like(x, 1e-9))
    pred = P['by_harmonic'][(3, 0)]['trapped']
    assert pred.any() and (~pred).any()
    # 写像 (px = 0 の初期条件 = φ0 = 0): J' = J − ∂V/∂φ、φ' = φ + 2π(ν0 + aJ') + ∂V/∂J
    N = 4000
    J = J0.copy()
    phi = np.zeros_like(J)
    th0 = 3 * phi
    winding = np.zeros_like(J)
    for n in range(N):
        e = h * (2 * J) ** 1.5 * np.exp(3j * phi)
        dVdphi = 2 * np.real(3j * e)
        J = J - dVdphi
        e2 = h * (2 * J) ** 1.5 * np.exp(3j * phi)
        dVdJ = 2 * np.real(1.5 * e2 / J)
        phi = phi + 2 * math.pi * (nu0 + a * J) + dVdJ
        th = 3 * phi - 2 * math.pi * 20 * (n + 1)
        winding = np.abs(th - th0)
    traj_trapped = winding < 2 * math.pi
    agree = np.mean(traj_trapped == pred)
    assert agree > 0.9, agree
    # 捕らえられる J0 の範囲の端が追跡と 0.15 jres 以内で一致
    assert abs(J0[pred].min() - J0[traj_trapped].min()) < 0.15 * jres and abs(J0[pred].max() - J0[traj_trapped].max()) < 0.15 * jres


def test_first_order_terms_fodo():
    from pysad import resdrive as R
    from pysad.parser import load_sad_string
    lat = load_sad_string(_fodo_src()).build_lattice('RING')
    rows = lat.twiss_table()
    harms, res = FR.first_order_terms(lat, rows)
    assert set(harms) == {(3, 0), (1, 0), (1, 2), (1, -2)}
    assert abs(abs(harms[(3, 0)][0][4]) - abs(res.total['h30000'])) < 1e-12
    r2 = R.resonance_terms(lat, terms=['h30000'], rows=rows)
    assert abs(r2.total['h30000'] - res.total['h30000']) < 1e-12
    # 作動点は 20/3 のすぐ下、3νx の共鳴線までの距離
    nux = rows[-1].mux / 2 / math.pi
    assert abs(3 * nux - 20) < 0.03


def test_coverage_stats_synthetic():
    nx, ny = 30, 20
    rng = np.random.default_rng(0)
    D = np.full((nx, ny), -7.0)
    M = np.zeros((nx, ny), dtype=bool)
    M[10:14, 3:15] = True
    D[M] = np.where(rng.random(M.sum()) < 0.5, -4.0, -7.0)
    s = FR.coverage_stats(D.ravel(), M.ravel(), (nx, ny), thr=-5.0, nshift=300)
    assert s['cov'] == 1.0 and 0.05 < s['area'] < 0.12 and s['enrich'] > 8
    assert s['pval'] < 0.05 and s['shift_mean'] < 0.3
    assert s['p_out'] == 0.0 and 0.3 < s['p_in'] < 0.7


def test_lines_distance_and_other_mn():
    d = FR.lines_distance([0.3334, 0.5], [0.2, 0.25], [(3, 0), (1, 0)])
    assert d[0] < 1e-3 and abs(d[1] - 0.5 / 3) < 1e-3 or d[1] < 0.2
    oth = FR.other_mn(4)
    assert (1, -1) in oth and (3, 0) not in oth and (1, 2) not in oth


@pytest.mark.skipif(not os.path.exists(os.path.join(EX, 'S1_fmap_rdt_fodo.table.txt')), reason='S1 の出力がない')
def test_fodo_fmap_ridge_matches_3nux_prediction():
    """FODO (QF の K1 = 0.697、νx = 6.6578): D > -5 の点がすべて 3νx の予想の領域に入り、共鳴の中心の予想 (ADTS) が D の尾根の x と 1 mm 以内"""
    from pysad.parser import load_sad_string
    lat = load_sad_string(_fodo_src()).build_lattice('RING')
    rows = lat.twiss_table()
    pr = FR.ResonancePredictor.from_lattice(lat, 1, rows)
    r = FMapResult.load_table(os.path.join(EX, 'S1_fmap_rdt_fodo.table.txt'))
    nu = (r.nux1 + 6, r.nuy1 + 4)
    P = pr.predict(r.x, r.y, nu=nu)
    s = FR.coverage_stats(r.D, P['any'], r.grid, thr=-5.0)
    assert s['n_hi'] >= 30 and s['cov'] > 0.95 and s['enrich'] > 2.5 and s['p_out'] < 0.01
    nx, ny = r.grid[0], r.grid[1]
    X, D = r.x.reshape(nx, ny), r.D.reshape(nx, ny)
    d0 = pr.predict(r.x, r.y)['by_harmonic'][(3, 0)]['delta0'].reshape(nx, ny)
    diffs = []
    for j in range(0, ny, 3):
        col = D[:, j]
        if not (col > -5).any():
            continue
        idx = np.nonzero(np.sign(d0[:-1, j]) * np.sign(d0[1:, j]) < 0)[0]
        if len(idx) == 0:
            continue
        i = idx[0]
        zc = X[i, j] + (X[i + 1, j] - X[i, j]) * d0[i, j] / (d0[i, j] - d0[i + 1, j])
        diffs.append(abs(zc - X[np.argmax(col), j]))
    assert len(diffs) >= 5 and max(diffs) < 1.0e-3 and np.mean(diffs) < 0.6e-3


def test_adts_tune_prediction_small_amplitude_astrid2():
    """線形の ADTS ν0 + A·J は、追跡 (NAFF) のチューンと、小振幅 (νの動きが 0.003 以下) で 5e-4 以内 (保存済みの Q1 の表)"""
    from pysad.parser import load_sad_string
    p = os.path.join(EX, 'Q1_fmap_astrid2.table.txt')
    if not os.path.exists(p):
        pytest.skip('Q1 の出力がない')
    lat = load_sad_string(_src(ASTRID)).build_lattice('RING')
    pr = FR.ResonancePredictor.from_lattice(lat, 2)
    r = FMapResult.load_table(p)
    Jx, Jy, _a, _b = pr.actions(r.x, r.y)
    nl = pr.tunes_lin(Jx, Jy)
    small = (Jx + Jy) < 3e-6
    assert small.sum() > 30
    ex = np.abs(nl[0] - 5 - r.nux1)[small]
    ey = np.abs(nl[1] - 2 - r.nuy1)[small]
    assert np.median(ex) < 5e-4 and np.median(ey) < 5e-4


# ---------------------------------------------------------------------------------------------------------------------
def _astrid_setup(src=None, line='RING'):
    return TrackSetup.from_source(src or _src(ASTRID), line, rad=False, rfsw=False, alost=0.1)[0]


def test_alignment_errors_run_and_zero_is_identity():
    st0 = _astrid_setup()
    xs, ys = [2e-3, 8e-3], [2e-3, 6e-3]
    r0 = frequency_map(st0, 2, x=xs, y=ys, nturn=256)
    st1, cnt = FR.apply_alignment_errors(st0, seed=3, quad_dxy=0.0, quad_rot=0.0)
    assert cnt == {'quad': 0, 'sext': 0, 'bend': 0}
    r1 = frequency_map(st1, 2, x=xs, y=ys, nturn=256)
    assert np.allclose(r0.D, r1.D, atol=0, rtol=0)
    st2, cnt = FR.apply_alignment_errors(st0, seed=3, quad_dxy=50e-6, quad_rot=2e-4, sext_dxy=50e-6, bend_dxy=50e-6, bend_rot=1e-4)
    assert cnt['quad'] == 24 and cnt['sext'] == 24 and cnt['bend'] > 0
    # 元の setup は変わらない
    assert np.allclose(frequency_map(st0, 2, x=xs, y=ys, nturn=256).D, r0.D, atol=0, rtol=0)
    c0 = st0.closed_orbit_dp(2, 0.0)
    c2 = st2.closed_orbit_dp(2, 0.0)
    assert np.abs(c0).max() < 1e-12 and np.abs(c2[:4]).max() > 1e-6          # 閉軌道がずれる (全周の追跡)
    assert st2.info['orbit_resid'] < 1e-10
    r2 = frequency_map(st2, 2, x=xs, y=ys, nturn=256)
    assert np.isfinite(r2.D).all() and r2.cond['nlost'] == 0
    assert not np.allclose(r2.D, r0.D)


def test_delta_map_runs_and_tune_follows_chromaticity():
    from pysad.parser import load_sad_string
    src = _src(ASTRID)
    st = _astrid_setup(src)
    dps = np.array([-0.01, 0.0, 0.01])
    r = frequency_map(st, 2, x=[0.5e-3], y=[0.5e-3], dp=dps, nturn=256)
    assert r.grid == (1, 1, 3) and np.isfinite(r.D).all() and r.cond['nlost'] == 0
    lat = load_sad_string(src).build_lattice('RING')
    for k, d in enumerate(dps):
        rows = lat.twiss_table(delta=float(d))
        nux, nuy = rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi
        assert abs(r.nux1[k] - (nux - 5)) < 2e-3 and abs(r.nuy1[k] - (nuy - 2)) < 2e-3
    # x–δ の格子
    r2 = frequency_map(st, 2, x=np.linspace(1e-3, 6e-3, 3), y=[2e-3], dp=dps, nturn=256)
    assert r2.grid == (3, 1, 3) and np.isfinite(r2.D).all()


def test_k64_mode_runs_same_map():
    st = _astrid_setup()
    xs, ys = [2e-3, 10e-3], [2e-3, 6e-3]
    rv = frequency_map(st, 2, x=xs, y=ys, nturn=256)
    compat.set_mode('k64')
    st64 = _astrid_setup()
    rk = frequency_map(st64, 2, x=xs, y=ys, nturn=256)
    assert np.isfinite(rk.D).all()
    assert np.allclose(rk.nux1, rv.nux1, atol=1e-6) and np.allclose(rk.nuy1, rv.nuy1, atol=1e-6)


def test_mpw_ring_runs():
    st = _astrid_setup(_src(ASTRID_ID), 'RINGMPW')
    nux, nuy = st.tunes_full()
    assert abs(nux - 5.18001) < 1e-4 and abs(nuy - 2.24920) < 1e-4
    r = frequency_map(st, 2, x=[2e-3, 8e-3], y=[2e-3], nturn=256)
    assert np.isfinite(r.D).all()


def test_naff_window_and_turns_regular_point():
    """規則的な点の D は周回数・窓を変えても同じ向きに下がる (周回数を 2 倍にすると D の下限が下がる)"""
    from pysad.freqmap import _tunes_of
    st = _astrid_setup()
    tw = st.twiss_at(2)
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))
    c = np.zeros((1, 6))
    c[0, 0], c[0, 2] = 3e-3, 2e-3
    rec = st.track(2, c, 2048, record=True)['rec'][:, 0, :4]

    def D(N, w):
        f, _ = _tunes_of(rec[:N], twp, w)
        return math.log10(math.hypot(f[2] - f[0], f[3] - f[1]))
    assert D(2048, 1) < D(1024, 1) < D(512, 1) < -6
    assert abs(D(1024, 1) - D(1024, 2)) < 3.0


# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not os.path.exists(os.path.join(EX, 'S1_gs_fmap_kekls.out.txt')), reason='S1 の gs の出力がない')
def test_kekls_gs_fmap_tracks_match():
    """KEK-LS の周波数マップの追跡を gs と比べた保存済みの出力: 規則的な点 (D < -3.5) の座標が 1e-9 m 以内、D の差が小さい、失われる点 (1 周目) の落ちた周が一致、
    カオス的な点の食い違いは初期値 1e-13 m の摂動の増幅と同じ大きさ、gs の終了後に残りなし"""
    import re
    txt = open(os.path.join(EX, 'S1_gs_fmap_kekls.out.txt'), encoding='utf-8').read()
    m = re.search(r'規則的な点 \(D < -3.5\) の座標列の最大の差: ([0-9.e+-]+)、\|D\(gs の座標\) − D\(frequency_map\)\| の最大 = ([0-9.e+-]+)', txt)
    assert m and float(m.group(1)) < 1e-9 and float(m.group(2)) < 1e-5
    assert re.search(r'落ちた周: API 1 / 台本 1 / gs 1', txt)
    m = re.search(r'座標列の最大の差 \(gs と pysad の台本、残った点.*?: ([0-9.e+-]+)', txt)
    gs_diff = float(m.group(1))
    m = re.search(r'dx = 1e-13: 最大の差 (.*?)、落ちた周', txt)
    pert = max(float(v) for v in m.group(1).split() if v != 'nan')
    assert gs_diff < 10 * pert                         # gs との差は 1e-13 m の摂動の増幅の 10 倍以内
    assert 'sad1.exe の残り: なし' in txt
