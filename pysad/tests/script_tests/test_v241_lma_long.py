# -*- coding: utf-8 -*-
"""v2.4.1: 誤差つき LMA と 10^4 周級の長時間追跡 (放射・FLUC) の gs との比較 (保存した gs の出力 + 短い live の追跡)。
保存物は scripts/v241_lma_* (台本 .sad、gs の出力 .gs.txt、pysad の出力 .pysad.txt (v2)・.pysadk64.txt)。生成は `_v241_lma_errors_run.py`・`_v241_lma_long_run.py`・`_v241_lma_k64_run.py`。
誤差つき LMA: ASTRID2、4 つの種 × 6 位置 × ±、2624 周、gs の二分法 (9 回) と pysad k64 の同じ台本が 48 / 48 で同じ区間。
長時間追跡: KEK-LS 全周・ASTRID2 の 10 粒子 (NOFLUC)・40 粒子 (FLUC) の 10^4 周で、gs と pysad k64 の座標が 1e-13 (NOFLUC)・1e-14 (FLUC、K64 の乱数列まで同じ)。
v2 は K64 の誤り (放射・四極の積分器) を直すので gs と違う (差の大きさを確認する)。短い追跡 (250 周) は live、長い比較は保存物の読み込み。
"""
import os
import sys

import numpy as np
import pytest

from pysad import compat

HERE = os.path.dirname(os.path.abspath(__file__))
SP = os.path.join(HERE, 'scripts')
sys.path.insert(0, HERE)
import _v241_lma_errors_summary as ES       # noqa: E402
import _v241_lma_long_analyze as LA         # noqa: E402
import _v241_lma_tools as T                 # noqa: E402

SLOW = pytest.mark.skipif(not os.environ.get('PYSAD_SLOW'), reason='誤差つき LMA の live の二分法 (k64、約 1 分)。PYSAD_SLOW=1 で走らせる')


def _need(*names):
    for n in names:
        if not os.path.exists(os.path.join(SP, n)):
            pytest.skip('保存物がない: ' + n)


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat._mode = old


# ---------------------------------------------------------------------------------------------
# 誤差つき LMA
# ---------------------------------------------------------------------------------------------
def test_err_lma_saved_gs_vs_pysad_k64_identical():
    """4 つの種 × 6 位置 × ± = 48 区間が gs と pysad k64 で完全に同じ (二分法 9 回の幅 3.9e-5 の区間)"""
    n = 0
    for s in (1, 2, 3, 4):
        base = 'v241_lma_err_s%d_t2624' % s
        _need(base + '.gs.txt', base + '.pysadk64.txt')
        g = ES.parse_b(os.path.join(SP, base + '.gs.txt'))
        k = ES.parse_b(os.path.join(SP, base + '.pysadk64.txt'))
        keys, d = ES.compare(g, k)
        assert len(keys) == 12
        assert np.all(d == 0), (s, d)
        assert all(g[x] == k[x] for x in keys)
        n += len(keys)
    assert n == 48


def test_err_lma_saved_gs_vs_pysad_v2_statistics():
    """v2 (四極の 4 次の積分器・放射の修正) は gs の区間と 38 / 48 が同じ、44 / 48 が 1 刻み (3.9e-5) 以内、最大の差 7e-4 (0.07%)。δ± の平均・標準偏差は同じ"""
    same = adj = 0
    worst = 0.0
    vg = {1: [], -1: []}
    vp = {1: [], -1: []}
    for s in (1, 2, 3, 4):
        base = 'v241_lma_err_s%d_t2624' % s
        _need(base + '.gs.txt', base + '.pysad.txt')
        g = ES.parse_b(os.path.join(SP, base + '.gs.txt'))
        p = ES.parse_b(os.path.join(SP, base + '.pysad.txt'))
        for k in g:
            d = abs(p[k][0] - g[k][0])
            same += d < 1e-9
            adj += 1e-9 < d < 1.5 * ES.BIN
            worst = max(worst, d)
            vg[k[1]].append(0.5 * sum(g[k]))
            vp[k[1]].append(0.5 * sum(p[k]))
    assert same >= 36 and same + adj >= 44
    assert worst < 1e-3
    for sg in (1, -1):
        a, b = np.abs(vg[sg]), np.abs(vp[sg])
        assert abs(a.mean() - b.mean()) < 0.01 * a.mean()
        assert abs(a.std(ddof=1) - b.std(ddof=1)) < 0.05 * a.std(ddof=1)


@SLOW
def test_err_lma_live_k64_matches_saved_gs():
    """seed 1・位置 17 の ± を、台本の二分法 (誤差の値を要素定義として書いた台本) で pysad k64 が live で再現する (保存した gs の区間と同じ)"""
    _need('v241_lma_err_s1_t2624.gs.txt')
    import _v241_lma_errors_run as ER
    script, st, info, zb = ER.make_script(1, 2624, [17])
    txt, _t = T.run_pysad_text(script, mode='k64')
    b = ER.parse_b(txt)
    g = ES.parse_b(os.path.join(SP, 'v241_lma_err_s1_t2624.gs.txt'))
    for sg in (1, -1):
        assert b[(17, sg)] == g[(17, sg)]


def test_err_lma_turn_dependence_saved():
    """誤差つきの LMA は周回数に依存する (1312 周 = 2 周期 → 2624 周 = 4 周期): 位置 17 (δ−) と 150 の限界が 1e-4〜6e-4 動く。長直線部・バケットの位置は動かない"""
    for s in (1, 2):
        a, b = 'v241_lma_err_s%d_t2624.gs.txt' % s, 'v241_lma_err_s%d_t1312.gs.txt' % s
        _need(a, b)
        g4 = ES.parse_b(os.path.join(SP, a))
        g2 = ES.parse_b(os.path.join(SP, b))
        diff = [abs(g4[k][0] - g2[k][0]) for k in g4]
        assert max(diff) > 5e-5
        assert max(diff) < 1e-3
        assert g4[(2, 1)] == g2[(2, 1)]


# ---------------------------------------------------------------------------------------------
# 長時間追跡
# ---------------------------------------------------------------------------------------------
def _stored(ring, mode, n, fn):
    nm = 'v241_lma_long_%s_%s_n%d_t10000%s' % (ring, mode, n, fn)
    _need(nm)
    return LA.load(os.path.join(SP, nm))


@pytest.mark.parametrize('ring,mode,n', [('kekls', 'nofluc', 10), ('astrid2', 'nofluc', 10), ('kekls', 'fluc', 40), ('astrid2', 'fluc', 40)])
def test_long_saved_gs_vs_pysad_k64(ring, mode, n):
    """10^4 周・全周の追跡 (放射。FLUC は量子励起つき): gs と pysad k64 の座標が全 40 点 (250 周ごと) で 1e-12 以内 (FLUC でも K64 の乱数列が同じ)"""
    g = _stored(ring, mode, n, '.gs.txt')
    k = _stored(ring, mode, n, '.pysadk64.txt')
    assert g['X'].shape == (40, n, 6) and np.all(g['flag'] == 1)
    assert np.abs(g['X'] - k['X']).max() < 1e-12
    if mode == 'fluc':
        assert np.abs(g['X'][-1, :, 5]).max() > 1e-4        # 量子励起で実際に δ が広がっている (自明な 0 の一致でない)


@pytest.mark.parametrize('ring,tol_250,tol_end', [('kekls', 1e-5, 1e-3), ('astrid2', 1e-4, 1e-3)])
def test_long_saved_v2_differs_from_gs_by_v2_fixes(ring, tol_250, tol_end):
    """v2 は K64 の誤り (四極の積分器 B226・B243 と放射の修正) を直すので gs と同じではない: 250 周で 4e-6 (KEK-LS)・1.7e-5 (ASTRID2)、10^4 周で 3e-4・2e-4。0 ではない"""
    g = _stored(ring, 'nofluc', 10, '.gs.txt')
    v = _stored(ring, 'nofluc', 10, '.pysad.txt')
    d = np.abs(g['X'] - v['X'])
    assert 1e-8 < d[0].max() < tol_250
    assert 1e-5 < d[-1].max() < tol_end


def test_long_saved_damping_rates_vs_emittance():
    """追跡の x の減衰率 (10 粒子の J の対数の傾き) と Emittance[] の α: gs (= k64) は KEK-LS 0.945 倍・ASTRID2 0.762 倍 (K64 の追跡の x の分配が Emittance[] の放射積分より小さい。B76 (K1≠0 の BEND の放射) の候補、切り分けは未実施)、v2 は 0.99・1.04 倍"""
    exp = {'kekls': ((0.93, 0.96), (0.97, 1.01)), 'astrid2': ((0.70, 0.82), (0.95, 1.12))}
    for ring, (rng_gs, rng_v2) in exp.items():
        tw, cod = LA.twiss(ring)
        em = _stored(ring, 'nofluc', 10, '.gs.txt')['emit']
        for fn, rng in (('.gs.txt', rng_gs), ('.pysad.txt', rng_v2)):
            d = _stored(ring, 'nofluc', 10, fn)
            X = d['X'] - cod[None, None, :]
            x, px, dp = X[:, :, 0], X[:, :, 1], X[:, :, 5]
            xb, pb = x - tw['etx'] * dp, px - tw['etpx'] * dp
            Jx = (1 + tw['alfx'] ** 2) / tw['betx'] * xb ** 2 + 2 * tw['alfx'] * xb * pb + tw['betx'] * pb ** 2
            sl = np.array([-np.polyfit(d['n'], np.log(Jx[:, i]), 1)[0] / 2 for i in range(Jx.shape[1])])
            r = sl.mean() / em['alpha'][0]
            assert rng[0] < r < rng[1], (ring, fn, r)


def _live(nm, k64=True):
    _need(nm + '.sad', nm + '.gs.txt')
    txt = open(os.path.join(SP, nm + '.sad'), encoding='utf-8').read().replace('{k, 1, 40}', '{k, 1, 1}')
    g = LA.load(os.path.join(SP, nm + '.gs.txt'))
    out, _ = T.run_pysad_text(txt, mode='k64' if k64 else 'v2')
    tmp = os.path.join(SP, '_tmp_%s_live.txt' % nm)
    try:
        open(tmp, 'w', encoding='utf-8').write(out)
        k = LA.load(tmp)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return g['X'][0], k['X'][0]


def test_long_live_nofluc_250_turns_k64_matches_gs():
    """KEK-LS 全周 (3563 要素) の 10 粒子を 250 周 (放射・RF) live で追跡して、保存した gs の 250 周の座標と 1e-12 以内 (k64)"""
    g, k = _live('v241_lma_long_kekls_nofluc_n10_t10000')
    assert np.abs(k - g).max() < 1e-12


def test_long_live_fluc_250_turns_k64_matches_gs():
    """KEK-LS 全周の FLUC つき 40 粒子を 250 周 live (k64) で追跡して、保存した gs の座標と 1e-12 以内 (量子励起の乱数列が gs と同じ)"""
    g, k = _live('v241_lma_long_kekls_fluc_n40_t10000')
    assert np.abs(k - g).max() < 1e-12
    assert np.abs(k[:, 5]).max() > 1e-5


def test_long_fluc_equilibrium_vs_emittance_saved():
    """FLUC の σδ²・εx が Emittance[] の平衡値の指数的な立ち上がり (1 − exp(−2αn)) に従う (KEK-LS、10^4 周 = 0.8 τx)。N = 40 の gs: σδ² 0.86、εx 0.90。許容は標本の統計 (σδ² 22%/点、εx 16%/点) の数倍"""
    g = _stored('kekls', 'fluc', 40, '.gs.txt')
    tw, cod = LA.twiss('kekls')
    s = LA.stats(g, tw, cod)
    em = g['emit']
    t = s['n']
    m = t >= 1500
    fz = 1 - np.exp(-2 * em['alpha'][2] * t)
    fx = 1 - np.exp(-2 * em['alpha'][0] * t)
    assert 0.6 < (s['sd'] ** 2 / (em['sigd'] ** 2 * fz))[m].mean() < 1.4
    assert 0.6 < (s['ex'] / (em['eps'][0] * fx))[m].mean() < 1.4
    assert s['ey'].max() < 1e-30            # 垂直方向は結合なし・光子の角度なし (K64 も同じ)


def _fluc_ratios(ring, n, fn):
    """t ≥ 1500 周の平均: (σδ²/予想、εx/予想、σz²/予想)。予想 = Emittance[] の平衡値 × (1 − exp(−2αn))"""
    d = _stored(ring, 'fluc', n, fn)
    tw, cod = LA.twiss(ring)
    s = LA.stats(d, tw, cod)
    # 予想は k64 (gs があれば gs) の Emittance[] の値に統一 (v2 の Emittance[] ではなく。比較表 _v241_lma_long_compare.py と同じ)
    em = _stored(ring, 'fluc', n, '.gs.txt' if n == 40 else '.pysadk64.txt')['emit']
    t = s['n']
    m = t >= 1500
    fz = 1 - np.exp(-2 * em['alpha'][2] * t)
    fx = 1 - np.exp(-2 * em['alpha'][0] * t)
    return ((s['sd'] ** 2 / (em['sigd'] ** 2 * fz))[m].mean(), (s['ex'] / (em['eps'][0] * fx))[m].mean(),
            (s['sz'] ** 2 / (em['sigz'] ** 2 * fz))[m].mean(), s['ey'].max(), s['nalive'][-1])


def test_long_fluc_equilibrium_astrid2_saved_gs_k64_v2():
    """ASTRID2 FLUC N = 40 (10^4 周 = 0.07 τx。平衡の比較でなく立ち上がりの傾き): gs = pysad k64 は σδ² 0.950・εx 1.195・σz² 0.953、v2 は 0.803・0.943・0.831。
    1 点の統計誤差は σδ² で 22%・εx で 16% (N = 40) なので、gs の 1.195 も v2 の 0.80 も 1 と区別できない (gs と v2 の違いも統計の範囲)。
    時間の点は同じ 40 粒子の追跡なので互いに独立でなく、文書の 6 点の標準誤差 (0.10・0.04) は過小"""
    exp = {'.gs.txt': (0.950, 1.195, 0.953), '.pysadk64.txt': (0.950, 1.195, 0.953), '.pysad.txt': (0.803, 0.943, 0.831)}
    for fn, e in exp.items():
        r = _fluc_ratios('astrid2', 40, fn)
        assert r[:3] == pytest.approx(e, abs=0.005), (fn, r)
        assert r[3] < 1e-60 and r[4] == 40
        assert 0.6 < r[0] < 1.4 and 0.6 < r[1] < 1.4        # 平衡の予想の 40% 以内 (統計 N=40 の範囲)
    g, v = _stored('astrid2', 'fluc', 40, '.gs.txt'), _stored('astrid2', 'fluc', 40, '.pysad.txt')
    d = np.abs(g['X'] - v['X'])
    assert d[0].max() > 1e-3 and d[-1].max() > 1e-2          # v2 は gs と別の軌道 (放射・積分器の修正。0 の一致でない)


def test_long_fluc_kekls_n120_saved_k64_v2():
    """KEK-LS FLUC N = 120 (gs なし。pysad k64 と v2): k64 σδ² 0.990・εx 0.951、v2 1.035・1.016 (N = 40 の gs 0.86・0.90、v2 1.18・0.92 より 1 に近い = N = 40 の揺らぎ)。
    1 点の統計は N = 120 で σδ² 13%・εx 9%"""
    exp = {'.pysadk64.txt': (0.990, 0.951, 1.028), '.pysad.txt': (1.035, 1.016, 0.897)}
    for fn, e in exp.items():
        r = _fluc_ratios('kekls', 120, fn)
        assert r[:3] == pytest.approx(e, abs=0.005), (fn, r)
        assert r[4] == 120 and r[3] < 1e-30
        assert 0.8 < r[0] < 1.2 and 0.8 < r[1] < 1.2
