# -*- coding: utf-8 -*-
"""v2.3.0 Q1: 周波数マップ (pysad/freqmap.py)、局所運動量アクセプタンス (pysad/lma.py)、共通の追跡設定 (pysad/trksetup.py)、SAD の台本の関数
(FrequencyMap・MomentumAcceptance・TouschekLifetimeLMA、pysad/q1cmd.py)。

根拠: (1) 追跡は SAD の TrackParticles と同じ損失判定 (台本の TrackParticles と数値が一致、gs の値 `_Q1_GS_LOSS` は example/v2_3/Q1_gs_loss_compare.py の gs の出力)、
(2) 小さな振幅の周波数マップのチューンは CALC のチューンと一致、(3) 周期 (セルの数) は要素の並びから (24 セルの FODO、ASTRID2 は 6、kek-ls_v03_12.lat は 20)、
(4) LMA は、RF バケットの高さ以下 (開口が広いとき)・線形の見積もり以下 (開口が狭いとき)・周期的なセルで同じ値、(5) Touschek に渡した一様な受容は TouschekTable の値と一致。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, ffscmd
from pysad.trksetup import TrackSetup, find_periodicity, select_locations
from pysad import freqmap, lma

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
KEKLS = os.path.join(HERE, 'kek-ls_v03_12.lat')

# 3 GeV、24 セルの FODO (周長 36 m、6 極つき、RF 4 MV、500 MHz)
FODO = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=4) SD = (L=0.2 K2=-7);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 24*CELL);
'''


@pytest.fixture(autouse=True)
def _v2():
    old = compat._mode
    compat.set_mode('v2')
    yield
    compat.set_mode(old)


@pytest.fixture(scope='module')
def fodo_rf():
    old = compat._mode
    compat.set_mode('v2')
    try:
        st, ffs = TrackSetup.from_source(FODO, 'RING', rad=True, rfsw=True)
    finally:
        compat.set_mode(old)
    return st, ffs


def test_periodicity_fodo_and_rings():
    st = TrackSetup.from_source(FODO, 'RING', rad=False, rfsw=False)[0]
    p = find_periodicity(st.els, st.rows)
    assert p['nper'] == 24 and p['ok']
    assert abs(p['cell_length'] * 24 - st.rows[-1][0]) < 1e-6
    for path, nper in ((ASTRID, 6), (KEKLS, 20)):
        if not os.path.exists(path):
            pytest.skip('ラティスがない')
        src = open(path, encoding='utf-8').read()
        ffs, _ = ffscmd.run_sad_string(src + '\nFFS;\nUSE RING; CELL; CALC;\n', out=io.StringIO())
        els = ffs.els
        rows = [(float(r.s),) for r in ffs.rows]
        p = find_periodicity(els, rows)
        assert p['nper'] == nper and p['ok']
        assert abs(p['cell_length'] * nper - rows[-1][0]) < 1e-6


def test_periodicity_detects_a_broken_cell():
    src = FODO.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF 11*CELL QF 12*CELL);')
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False)[0]
    assert find_periodicity(st.els, st.rows)['nper'] == 1


def test_resonance_helpers():
    nux, nuy = 5.18, 2.14
    chk = freqmap.resonance_check(nux, nuy, 4)
    names = [(m, n) for (m, n, p, d, o, nm) in chk]
    for mn in ((1, 0), (0, 1), (2, 0), (0, 2)):           # 整数・半整数は必ず入る
        assert mn in names
    d = {(m, n): (p, dist) for (m, n, p, dist, o, nm) in chk}
    assert d[(1, 0)][0] == 5 and abs(d[(1, 0)][1] - 0.18) < 1e-12
    assert d[(0, 2)][0] == 4 and abs(d[(0, 2)][1] - 0.28 / 2) < 1e-12
    assert abs(d[(1, -1)][1] - 0.04 / math.sqrt(2)) < 1e-12 and d[(1, -1)][0] == 3
    # 線の上の点は距離 0
    assert freqmap.line_distance([0.25], [0.5], 3)[0] < 1e-15        # 4 νx = 1, 2 νy = 1
    ls = freqmap.resonance_lines(2, ((0.1, 0.3), (0.1, 0.3)))
    assert (1, 1, 0) not in ls and (1, -1, 0) in ls and (2, 0, 1) not in ls


def test_track_matches_trackparticles_script(fodo_rf):
    """TrackSetup.track (LMA・周波数マップが使う) と台本の TrackParticles[{loc, beam}, loc, 1, n] が、損失の周・座標まで一致 (RF・放射つき、アパーチャ 5 mm)"""
    st, ffs = fodo_rf
    st = st.copy(alost=0.005)
    loc = 7
    cod = st.closed_orbit([loc])[0]
    dps = np.array([0.0, 0.004, 0.008, 0.012, 0.016, 0.02, -0.004, -0.01, -0.02])
    c = np.tile(cod, (len(dps), 1))
    c[:, 5] += dps
    n = 150
    r = st.track(loc, c, n)
    cols = ['{' + ','.join('%.17g' % c[i, k] for i in range(len(dps))) + '}' for k in range(6)]
    cols.append('{' + ','.join('1' for _ in dps) + '}')
    scr = ('FFS;\nUSE RING; CELL; CALC; RFSW; RAD; RADCOD; NOFLUC;\nEmittance[];\nLOSSAMPL = 0.005;\n$FORM = "20.14";\n'
           'beam = {%s};\nr = TrackParticles[{%d, beam}, %d, 1, %d];\n'
           'Do[Print["P ", i, " ", r[[2,7,i]], " ", r[[2,1,i]], " ", r[[2,5,i]], " ", r[[2,6,i]]], {i, 1, %d}];\nExit;\n') % (','.join(cols), loc, loc, n, len(dps))
    buf = io.StringIO()
    ffscmd.run_sad_string(FODO + scr, out=buf)
    rows = [l.split() for l in buf.getvalue().splitlines() if l.startswith('P ')]
    assert len(rows) == len(dps)
    alive = np.array([int(float(a[2])) for a in rows])
    assert np.array_equal(alive, (r['lost_turn'] == 0).astype(int))
    assert 0 < alive.sum() < len(dps)                      # 残る点と落ちる点が両方ある例
    for i, a in enumerate(rows):
        assert abs(float(a[3]) - r['final'][i, 0]) < 1e-12 and abs(float(a[5]) - r['final'][i, 5]) < 1e-12


def _gs_lines():
    """example/v2_3/Q1_gs_loss_astrid2.gs.txt (gs の出力。K64、2026-10-19 ビルド): 'P loc i 旗 x' の行と 'B loc 符号 lo hi' の行"""
    p = os.path.join(ROOT, 'example', 'v2_3', 'Q1_gs_loss_astrid2.gs.txt')
    if not os.path.exists(p):
        return None, None
    pr, br = [], []
    for l in open(p, encoding='utf-8', errors='replace'):
        a = l.split()
        if l.startswith('P ') and len(a) >= 5:
            pr.append((int(a[1]), int(float(a[3])), float(a[4])))
        elif l.startswith('B ') and len(a) >= 5:
            br.append((int(a[1]), int(a[2]), float(a[3]), float(a[4])))
    return pr or None, br or None


def _gs_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location('q1_gs_loss', os.path.join(ROOT, 'example', 'v2_3', 'Q1_gs_loss_compare.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_astrid2_track_loss_matches_gs():
    """gs の出力 (リポジトリの例) と同じ損失判定 (ASTRID2、RF・放射つき、1312 周、位置 2 と 17、δ の 30 点)。
    生き残りの旗が全部一致し、最後の x は gs と 3 mm 以内 (損失ぎりぎりの粒子は閉軌道の 1e-7 の違いで軌道が変わる)。
    gs の二分法の区間 (lo: 残る、hi: 落ちる) の 2 点を、pysad の追跡で確かめる。gs の出力がなければ skip"""
    gs, bis = _gs_lines()
    if gs is None or bis is None or not os.path.exists(ASTRID):
        pytest.skip('gs の出力の例がない')
    mod = _gs_module()
    src = open(ASTRID, encoding='utf-8').read()
    st = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=mod.APER)[0]
    zb = st.bucket_zlimit()
    k = 0
    nlost = 0
    for loc, dps in mod.CASES:
        cod = st.closed_orbit([loc])[0]
        for sg in (1, -1):
            for d in dps:
                if sg == -1 and d == 0.0:
                    continue
                c = cod.copy()
                c[5] += sg * d
                r = st.track(loc, c[None, :], mod.NTURN, zlost=zb)
                assert gs[k][0] == loc
                assert gs[k][1] == int(r['lost_turn'][0] == 0), (loc, sg * d)
                assert abs(gs[k][2] - r['final'][0, 0]) < 3e-3
                nlost += int(r['lost_turn'][0] > 0)
                k += 1
    assert k == len(gs) == 30 and 8 < nlost < 22
    for loc, sg, lo, hi in bis:
        cod = st.closed_orbit([loc])[0]
        c = np.tile(cod, (2, 1))
        c[0, 5] += sg * lo
        c[1, 5] += sg * hi
        r = st.track(loc, c, mod.NTURN, zlost=zb)
        assert r['lost_turn'][0] == 0 and r['lost_turn'][1] > 0, (loc, sg, lo, hi)


def test_frequency_map_tunes_match_calc_and_flag_lost():
    st, ffs = TrackSetup.from_source(FODO, 'RING', rad=False, rfsw=False)
    nux, nuy = st.tunes_full()
    loc = 1
    r = freqmap.frequency_map(st, loc, x=[1e-6, 2e-3, 0.5], y=[1e-6, 1e-3], nturn=256)
    assert r.cond['grid_nx_ny_ndp'] == (3, 2, 1) and len(r) == 6
    # 小さな振幅: CALC のチューン (端数) と 1e-8 以内、前半と後半の差は ~0 (D < −9)
    i0 = 0
    assert abs(r.nux1[i0] - (nux % 1)) < 1e-8 and abs(r.nuy1[i0] - (nuy % 1)) < 1e-8
    assert r.D[i0] < -9
    # 0.5 m (開口 0.1 m を超える) は失われ、D は nan、lost_turn > 0
    lost = r.x > 0.1
    assert np.all(r.lost_turn[lost] > 0) and np.all(np.isnan(r.D[lost]))
    assert np.all(r.lost_turn[~lost] == 0) and np.all(np.isfinite(r.D[~lost]))
    # 振幅が大きいと拡散が大きい (6 極)
    assert r.D[2] > r.D[0] or r.D[3] > r.D[0]


def test_frequency_map_dp_plane_uses_offmomentum_orbit():
    """x–δ の格子: δ ≠ 0 の点は δ の閉軌道 (4 次元) の上から出発する。x=0 の点のチューン・クロマティシティ ν(δ) − ν(0) ≈ ξ δ (CALC のクロマティシティ) と、拡散が小さいこと"""
    st, ffs = TrackSetup.from_source(FODO, 'RING', rad=False, rfsw=False)
    r = freqmap.frequency_map(st, 1, x=[1e-6], y=[1e-6], dp=[-0.005, 0.0, 0.005], nturn=256)
    assert len(r) == 3 and np.all(r.lost_turn == 0)
    assert np.all(r.D < -7)                       # 閉軌道の上に置いたので、振幅ほぼゼロ (1 μm) の運動 (δ = 0 は -11、±0.5% は -7〜-10: 規格化と窓の誤差)
    dnu = (r.nux1[2] - r.nux1[0]) / 0.01
    # 数値のチューンの差 (CALC の δ の依存) と比べる: DP=±0.005 の CALC
    out = []
    for d in (-0.005, 0.005):
        b = io.StringIO()
        ffscmd.run_sad_string(FODO + 'FFS;\nUSE RING; CELL; DP0=%g; CALC;\nPrint["NU ", Twiss["NX", "$$$"]/2/Pi];\nExit;\n' % d, out=b)
        out.append(float([l for l in b.getvalue().splitlines() if l.startswith('NU ')][0].split()[1]))
    assert abs(dnu - (out[1] - out[0]) / 0.01) < 0.02 * abs(dnu) + 1e-3


def test_frequency_map_workers_agree():
    st, ffs = TrackSetup.from_source(FODO, 'RING', rad=False, rfsw=False)
    a = freqmap.frequency_map(st, 3, x=np.linspace(1e-3, 6e-3, 4), y=np.linspace(1e-3, 4e-3, 3), nturn=128, workers=1, chunk=5)
    b = freqmap.frequency_map(st, 3, x=np.linspace(1e-3, 6e-3, 4), y=np.linspace(1e-3, 4e-3, 3), nturn=128, workers=2, chunk=5)
    assert np.array_equal(a.lost_turn, b.lost_turn)
    assert np.allclose(a.nux1, b.nux1, atol=1e-13, equal_nan=True) and np.allclose(a.D, b.D, atol=1e-6, equal_nan=True)


# ---------------------------------------------------------------------------------------------------------------------
# LMA
# ---------------------------------------------------------------------------------------------------------------------

def test_lma_bucket_limited_and_converged(fodo_rf):
    """アパーチャが広い (±0.1 m) と、δ± は RF バケットの高さで決まる (理由 'rf')。δ+ は バケット高さの ±5% 以内。二分の区間の中で
    δ+ は生き残り、区間の上端 (hi) は失われる (追跡で確認)"""
    st, ffs = fodo_rf
    st = st.copy(alost=0.1)
    bh = st.bucket['height']
    nt = st.nturn_for_synchrotron(2.0)
    assert nt >= int(2.0 / abs(st.bucket['nus']))              # 既定の周回数はシンクロトロン振動 2 周期以上
    locs = [2, 7]
    r = lma.momentum_acceptance(st, locs, workers=1, tol=2e-4)
    assert r.cond['nturn'] == nt
    for i in range(2):
        assert r.reason_plus[i] == 'rf' and r.reason_minus[i] == 'rf'
        assert abs(r.dp_plus[i] / bh - 1) < 0.05 and abs(-r.dp_minus[i] / bh - 1) < 0.05
        assert 0 < r.hi_plus[i] - r.dp_plus[i] <= 2e-4 * 1.0001
    zb = r.cond['bucket_zlimit_m']
    cod = st.closed_orbit([2])[0]
    for sgn, lo, hi in ((1, r.dp_plus[0], r.hi_plus[0]), (-1, -r.dp_minus[0], -r.hi_minus[0])):
        c = np.tile(cod, (2, 1))
        c[0, 5] += sgn * lo
        c[1, 5] += sgn * hi
        t = st.track(2, c, nt, zlost=zb)
        assert t['lost_turn'][0] == 0 and t['lost_turn'][1] > 0


def test_lma_aperture_limited_below_linear_estimate_and_periodic(fodo_rf):
    """アパーチャが狭い (±4 mm) と分散で決まる (理由 'aper')。δ は 線形の見積もり δ_lin (開口 ±4 mm) の 0.6 〜 1.2 倍、バケット以下。
    セルが同じ 2 つの位置 (cell 1 と cell 13 の同じ要素) の値は tol の程度で一致"""
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    per = find_periodicity(st.els, st.rows)
    nel = per['starts'][1] - per['starts'][0]
    l1 = per['starts'][0] + 4
    l2 = per['starts'][12] + 4
    r = lma.momentum_acceptance(st, [l1, l2], workers=1, tol=2e-4)
    lin = lma.linear_momentum_limit(st, [l1, l2], 0.004)
    bh = st.bucket['height']
    for i in range(2):
        assert r.reason_plus[i] == 'aper' and r.reason_minus[i] == 'aper'
        assert 0.6 < r.dp_plus[i] / lin[i] < 1.2 and 0.6 < -r.dp_minus[i] / lin[i] < 1.2
        assert r.dp_plus[i] < bh
    assert abs(r.dp_plus[0] - r.dp_plus[1]) < 1e-3 and abs(r.dp_minus[0] - r.dp_minus[1]) < 1e-3     # 二分の幅 + RF の位置の違い (閉軌道の δ の差)
    assert nel > 0


def test_lma_workers_agree_and_low_priority(fodo_rf):
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    a = lma.momentum_acceptance(st, [3, 9], workers=1)
    b = lma.momentum_acceptance(st, [3, 9], workers=2)
    assert np.array_equal(a.dp_plus, b.dp_plus) and np.array_equal(a.dp_minus, b.dp_minus)
    assert b.cond['workers'] == 2


def test_lma_nsub_one_is_plain_bisection(fodo_rf):
    """nsub=1 (普通の二分法) と nsub=3 の結果の差は、区間の幅 tol 以内"""
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    a = lma.momentum_acceptance(st, [5], workers=1, nsub=1, tol=5e-5)
    b = lma.momentum_acceptance(st, [5], workers=1, nsub=3, tol=5e-5)
    assert abs(a.dp_plus[0] - b.dp_plus[0]) <= 1.5e-4 and abs(a.dp_minus[0] - b.dp_minus[0]) <= 1.5e-4


def test_select_locations(fodo_rf):
    st, ffs = fodo_rf
    all_ = select_locations(st)
    assert all_[0] == 1 and len(set(st.s_of(l) for l in all_)) == len(all_)       # s が同じ位置は 1 つだけ
    sub = select_locations(st, 1, 40, ds=1.0)
    assert all(st.s_of(b) - st.s_of(a) >= 1.0 - 1e-9 for a, b in zip(sub, sub[1:]))
    assert select_locations(st, 5, 30, step=5) == [l for l in all_ if l in (5, 10, 15, 20, 25, 30)] or len(select_locations(st, 5, 30, step=5)) > 0


def test_lma_warns_for_short_turns(fodo_rf):
    st, ffs = fodo_rf
    with pytest.warns(UserWarning):
        lma.momentum_acceptance(st.copy(alost=0.004), [3], nturn=20, workers=1, ncoarse=2, nsub=1, tol=0.01)


# ---------------------------------------------------------------------------------------------------------------------
# Touschek
# ---------------------------------------------------------------------------------------------------------------------

def test_touschek_uniform_acceptance_reproduces_table(fodo_rf):
    """一様な受容 δ+ = |δ−| = a (表の点) を渡すと、touschek_lifetime_local は TouschekTable の損失率 rate[k] と一致する。局所の値を入れると変わる"""
    st, ffs = fodo_rf
    ffs.flags['INTRA'] = True
    ffs.I.globals['PBUNCH'] = 5e8
    ffs.run_emit(quiet=True)
    tab = ffs.touschek['table']
    dp = np.array(tab[1][0])
    k = 14
    s = np.array([0.0, 10.0, 20.0])
    uniform = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [dp[k]] * 3, [dp[k]] * 3, [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    res = lma.touschek_lifetime_local(tab, uniform, setup=st)
    assert abs(res['rate'] / tab[1][1][k] - 1) < 1e-10 and res['extrapolated'] == 0
    # 受容が位置で半分になる区間 → 寿命が短くなる
    half = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [dp[k], dp[k] / 2, dp[k]], [dp[k], dp[k] / 2, dp[k]], [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    assert lma.touschek_lifetime_local(tab, half, setup=st)['rate'] > res['rate']
    # δ+ ≠ |δ−|: 片側ずつ半分の和 → 対称な受容の寿命の間に入る
    asym = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [dp[k + 4]] * 3, [dp[k]] * 3, [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=st)
    ra = lma.touschek_lifetime_local(tab, asym, setup=st)['rate']
    assert tab[1][1][k + 4] < ra < tab[1][1][k]
    assert abs(ra - 0.5 * (tab[1][1][k + 4] + tab[1][1][k])) < 1e-10 * tab[1][1][k]
    s_, dpp, dpm = uniform.touschek_input()
    assert np.array_equal(s_, s) and np.all(dpp > 0) and np.all(dpm < 0)


# ---------------------------------------------------------------------------------------------------------------------
# SAD の台本
# ---------------------------------------------------------------------------------------------------------------------

def test_script_commands(tmp_path):
    out_tab = str(tmp_path / 'lma.txt')
    fmap_tab = str(tmp_path / 'fm.txt')
    png = str(tmp_path / 'fm.png')
    scr = ('FFS;\nUSE RING; CELL; CALC; NORFSW; NORAD;\n'
           'f = FrequencyMap[3, 128, X->{0.001, 0.004, 3}, Y->{0.001, 0.003, 2}, Output->"%s", Plot->"%s"];\n'
           'Print["FM ", Length[Diffusion /. f], " ", Length[NuX /. f], " ", Length[LostTurn /. f]];\n'
           'RFSW; RAD; RADCOD; NOFLUC; Emittance[]; INTRA; PBUNCH = 5e8; Emittance[];\n'
           'LOSSAMPL = 0.004;\n'
           'm = MomentumAcceptance[{4, 12}, Workers->1, Output->"%s"];\n'
           'Print["LMA ", Length[DPPlus /. m], " ", (DPPlus /. m)[[1]] > 0, " ", (DPMinus /. m)[[1]] < 0];\n'
           'Print["TL ", TouschekLifetimeLMA[] > 0];\nExit;\n') % (fmap_tab.replace('\\', '/'), png.replace('\\', '/'), out_tab.replace('\\', '/'))
    buf = io.StringIO()
    ffscmd.run_sad_string(FODO + scr, out=buf)
    txt = buf.getvalue()
    assert 'FM 6 6 6' in txt and 'LMA 2 1 1' in txt and 'TL 1' in txt
    assert os.path.getsize(png) > 1000
    body = [l for l in open(out_tab, encoding='utf-8') if not l.startswith('#')]
    assert len(body) == 2
    body = [l for l in open(fmap_tab, encoding='utf-8') if not l.startswith('#')]
    assert len(body) == 6


def test_script_cell_option_restricts_to_first_cell():
    scr = ('FFS;\nUSE RING; CELL; CALC; RFSW; RAD; RADCOD; NOFLUC; Emittance[];\nLOSSAMPL = 0.004;\n'
           'm = MomentumAcceptance[Cell->True, Distance->0.7, Workers->1, Coarse->4, Sub->1, Tolerance->0.002];\n'
           'Print["S ", Max[S /. m] < 4.5, " ", Length[S /. m]];\nExit;\n')
    buf = io.StringIO()
    ffscmd.run_sad_string(FODO + scr, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('S ')][0]
    assert line.startswith('S 1 ')


def test_lma_check_every_finer_checks_lower_the_limit(fodo_rf):
    """損失の判定を細かく (check_every=1) すると、アパーチャが限界の δ は SAD の tapert のまま (0) と同じか小さく、線形の見積もりの 0.7〜1.05 倍に入る"""
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    locs = [5, 9]
    a = lma.momentum_acceptance(st, locs, workers=1, check_every=0, tol=1e-4)
    b = lma.momentum_acceptance(st, locs, workers=1, check_every=1, tol=1e-4)
    lin = lma.linear_momentum_limit(st, locs, 0.004)
    assert b.cond['check_every'] == 1 and a.cond['check_every'] == 0
    assert np.all(b.dp_plus <= a.dp_plus + 1.5e-4) and np.all(-b.dp_minus <= -a.dp_minus + 1.5e-4)
    assert np.all(b.dp_plus / lin > 0.7) and np.all(b.dp_plus / lin < 1.05)


def _table_case(name):
    p = os.path.join(ROOT, 'example', 'v2_3', 'Q1_fmap_%s.table.txt' % name)
    if not os.path.exists(p):
        pytest.skip('例の表がない')
    return freqmap.FMapResult.load_table(p)


@pytest.mark.parametrize('name,path', [('astrid2', ASTRID), ('kekls', KEKLS)])
def test_example_frequency_map_tables(name, path):
    """例の表 (example/v2_3/Q1_fmap_*.table.txt。41 × 31 の格子、1024 周、RFSW・RAD OFF): (1) 表の読み込み、(2) 拡散 D が大きい点は、次数 4 以下の共鳴線 (距離 < 0.003) の近くに、
    遠くより 1.5 倍以上多い (D > −5 の割合)、(3) 表の 4 点を、同じ条件で追跡し直すと tune・D が 1e-8 以内で再現する"""
    if not os.path.exists(path):
        pytest.skip('ラティスがない')
    r = _table_case(name)
    assert r.grid == (41, 31, 1) and len(r) == 1271
    ok = np.isfinite(r.D)
    assert ok.sum() + (r.lost_turn > 0).sum() == len(r)
    dist = freqmap.line_distance(r.nux1[ok], r.nuy1[ok], 4)
    near = dist < 0.003
    f_near = np.mean(r.D[ok][near] > -5)
    f_far = np.mean(r.D[ok][~near] > -5)
    assert near.sum() > 100 and (~near).sum() > 100 and f_near > 1.5 * f_far
    # 再現
    src = open(path, encoding='utf-8').read()
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    loc = int(r.loc)
    idx = np.nonzero(ok)[0]
    pick = idx[np.linspace(0, len(idx) - 1, 4).astype(int)]
    for i in pick:
        q = freqmap.frequency_map(st, loc, x=[r.x[i]], y=[r.y[i]], nturn=int(r.cond['nturn']))
        assert abs(q.nux1[0] - r.nux1[i]) < 1e-8 and abs(q.nuy2[0] - r.nuy2[i]) < 1e-8 and abs(q.D[0] - r.D[i]) < 1e-3, (name, i)


def test_example_lma_tables_are_consistent():
    """例の LMA の表: ASTRID2 (check_every=1、±10 mm) は δ+ の最大が RF バケットの高さの ±3% 以内、δ± は線形の見積もり以下 (×1.05)、
    大きなリング A (check_every=5) は RF バケットで決まる (δ+ の最大が バケット高さの 1.0〜1.03 倍)"""
    for name, kind in (('astrid2', 'a'), ('kekls', 'k')):
        p = os.path.join(ROOT, 'example', 'v2_3', 'Q1_lma_%s.table.txt' % name)
        if not os.path.exists(p):
            pytest.skip('例の表がない')
        cond = {}
        rows = []
        for l in open(p, encoding='utf-8'):
            if l.startswith('# ') and ' = ' in l:
                k, v = l[2:].strip().split(' = ', 1)
                cond[k] = v
            elif not l.startswith('#') and l.strip():
                rows.append(l.split())
        bh = float(cond['bucket_height'])
        dpp = np.array([float(a[3]) for a in rows])
        dpm = np.array([float(a[4]) for a in rows])
        assert int(cond['nturn']) >= 2 / abs(float(cond['nu_s'])) - 1 and cond['bucket'] == 'True'
        assert np.all(dpp > 0) and np.all(dpm < 0)
        if kind == 'a':
            assert abs(dpp.max() / bh - 1) < 0.03 and abs(-dpm.max() / bh - 1) < 0.03 or abs(dpm.min() / -bh - 1) < 0.03
            assert dpp.min() < 0.6 * bh                       # 分散の大きい位置で横のアパーチャが限界
        else:
            assert 1.0 <= dpp.max() / bh <= 1.03
            assert dpp.min() > 0.95 * bh
