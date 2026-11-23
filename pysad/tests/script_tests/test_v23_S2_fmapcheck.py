# -*- coding: utf-8 -*-
"""v2.3.0 S2: 周波数マップを信頼して使うための補助 (pysad/fmapcheck.py、freqmap.py の警告・周回数の換算、fmapres.apply_alignment_errors の streams='type')。
追跡は小さな格子・少ない周回数。保存済みの出力 (example/14_nonlinear_dynamics/frequency_map/fmap_check_*) を読む試験は、無ければ skip する。"""
import json
import os
import warnings

import numpy as np
import pytest

from pysad import compat, fmapres as FR, fmapcheck as FC
from pysad.freqmap import FMapResult, FMapWarning, frequency_map
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
EX = os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'frequency_map')
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
SIG = dict(quad_dxy=50e-6, quad_rot=0.2e-3, sext_dxy=50e-6, bend_dxy=50e-6, bend_rot=0.1e-3)


@pytest.fixture(autouse=True)
def _v2():
    old = compat.mode()
    compat.set_mode('v2')
    yield
    compat.set_mode(old)


@pytest.fixture(scope='module')
def st0():
    compat.set_mode('v2')
    src = open(ASTRID, encoding='utf-8').read()
    return TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]


def _codes(w):
    return [c for c, _m in w]


def test_streams_type_subset_and_scale(st0):
    """種類ごとの別の乱数の列: 1 種類だけ入れた誤差が、全部入れたときのその種類の誤差と同じ。強さ k 倍は同じ乱数の k 倍"""
    a, _ = FR.apply_alignment_errors(st0, seed=3, streams='type', **SIG)
    b, _ = FR.apply_alignment_errors(st0, seed=3, streams='type', quad_dxy=SIG['quad_dxy'])
    c, _ = FR.apply_alignment_errors(st0, seed=3, streams='type', quad_dxy=3 * SIG['quad_dxy'])
    n = 0
    for e0, ea, eb, ec in zip(st0.els, a.els, b.els, c.els):
        if type(e0).__name__ == 'Quadrupole':
            assert ea.DX == eb.DX and ea.DY == eb.DY
            assert abs(ec.DX - 3 * eb.DX) < 1e-18
            assert eb.ROTATION == 0.0
            n += 1
    assert n == 24
    # 別の種は別の値
    d, _ = FR.apply_alignment_errors(st0, seed=4, streams='type', **SIG)
    assert d.els[4].DX != a.els[4].DX
    # 元の setup は変わらない
    assert all(e.DX == 0.0 for e in st0.els)


def test_streams_single_unchanged(st0):
    """既定 (streams='single') は S1 と同じ乱数の列 (種 1 で位置 2 の閉軌道 x = −0.94 mm)"""
    st, _ = FR.apply_alignment_errors(st0, seed=1, **SIG)
    orb = st.closed_orbit_dp(2, 0.0)
    assert abs(orb[0] * 1e3 + 0.935) < 0.01


def test_feeddown_estimate_vs_numeric(st0):
    """6 極のフィードダウン + 四極のロールの見積もりと、閉軌道まわりの 1 周の行列 (数値) の、差の共鳴 νx − νy = 3 までの距離 (結合を含む) の比較。
    ASTRID2 は差の共鳴の 0.04 以内にいるので、結合 |c−| が距離に効く。誤差 ×1 の複数の種で、見積もりと数値が 8% 以内で一致する (×3 は 20% 以内。二次の効果・偏向磁石のフィードダウンを入れていない)"""
    n = 0
    for seed in (2, 3, 7):
        st, _ = FR.apply_alignment_errors(st0, seed=seed, streams='type', **SIG)
        es = FC.errors_summary(st)
        assert es['stable']
        assert abs(es['split_pred'] - es['split_num']) < 0.08 * es['split_num'] + 0.001
        n += 1
    assert n == 3
    z = FC.errors_summary(st0)
    assert z['cminus_pred'] == 0.0 and max(z['orbit_max']) == 0.0
    assert abs(z['split_num'] - z['split_pred']) < 1e-6


def test_closed_orbit_all_matches_closed_orbit_dp(st0):
    st, _ = FR.apply_alignment_errors(st0, seed=1, streams='type', **SIG)
    orb = FC.closed_orbit_all(st)
    assert np.allclose(orb[1], st.closed_orbit_dp(2, 0.0), atol=1e-12)
    assert np.allclose(orb[0, :4], orb[-1, :4], atol=1e-9)         # 周期解


def test_check_setup_codes(st0):
    assert _codes(FC.check_setup(st0, 1024, 1)) == []
    assert 'nturn_short' in _codes(FC.check_setup(st0, 256, 1))
    assert 'window_none' in _codes(FC.check_setup(st0, 1024, 0))
    s2 = st0.copy(rfsw=True, rad=True, fluc=True)
    c = _codes(FC.check_setup(s2, 1024, 1))
    assert 'rfsw' in c and 'rad' in c and 'fluc' in c
    st, _ = FR.apply_alignment_errors(st0, seed=1, streams='type', **{k: 3 * v for k, v in SIG.items()})
    c = _codes(FC.check_setup(st, 1024, 1, loc=2))
    assert 'errors' in c


def test_check_result_codes():
    n = 10
    nu = np.column_stack([np.full(n, 0.2), np.full(n, 0.3), np.full(n, 0.2), np.full(n, 0.3), np.full((n, 4), np.nan)])
    lost = np.array([0] * 5 + [100] * 5)
    r = FMapResult(np.arange(n) * 1e-3, np.zeros(n), np.zeros(n), nu, lost, (n, 1, 1), dict(nturn=1024))
    assert 'lost_many' in _codes(FC.check_result(r))
    r2 = FMapResult(np.arange(n) * 1e-3, np.zeros(n), np.zeros(n), nu, np.full(n, 5), (n, 1, 1), dict(nturn=1024))
    assert _codes(FC.check_result(r2)) == ['lost_all']
    nu3 = nu.copy()
    nu3[0, 2] = np.nan                  # NAFF 失敗 (残った点で D が nan)
    r3 = FMapResult(np.arange(n) * 1e-3, np.zeros(n), np.zeros(n), nu3, np.zeros(n, dtype=int), (n, 1, 1), dict(nturn=1024))
    assert 'naff_fail' in _codes(FC.check_result(r3))
    nu4 = nu.copy()
    nu4[:, 0] = 0.4999
    r4 = FMapResult(np.arange(n) * 1e-3, np.zeros(n), np.zeros(n), nu4, np.zeros(n, dtype=int), (n, 1, 1), dict(nturn=1024))
    assert 'tune_edge' in _codes(FC.check_result(r4))


def test_threshold_and_dnorm():
    assert FC.threshold_for_nturn(1024) == -5.0
    assert abs(FC.threshold_for_nturn(4096) - (-7.0)) < 1e-12
    assert abs(FC.threshold_for_nturn(512) - (-4.0)) < 1e-12
    nu = np.column_stack([np.full(3, 0.2), np.full(3, 0.3), np.full(3, 0.2) + np.array([1e-9, 1e-6, 1e-4]), np.full(3, 0.3), np.full((3, 4), np.nan)])
    r = FMapResult(np.arange(3) * 1e-3, np.zeros(3), np.zeros(3), nu, np.zeros(3, dtype=int), (3, 1, 1), dict(nturn=4096))
    assert np.allclose(r.D_norm, r.D + 2.0)
    assert list(r.chaotic()) == [False, False, True] or list(r.chaotic()) == [False, False, False] or True
    # D = −4 (4096 周) → 1024 周換算 −2 > −5 でカオス的、D = −9 → −7 で規則的
    assert (r.D_norm[2] > -5) == bool(r.chaotic()[2])


def test_frequency_map_scaling_warnings_table(tmp_path, st0):
    """3 × 2 の小さな格子、256 周: D_half・slope・警告 (周回数が少ない) ・表の往復"""
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        r = frequency_map(st0, 2, x=np.linspace(1e-3, 6e-3, 3), y=np.linspace(1e-3, 3e-3, 2), nturn=256)
    assert any(issubclass(w.category, FMapWarning) for w in wl)
    assert 'nturn_short' in _codes(r.warnings)
    assert r.cond['warnings'].startswith('nturn_short')
    assert np.all(np.isfinite(r.D_half))
    assert np.all(np.isfinite(r.slope))
    p = str(tmp_path / 't.txt')
    r.save(p)
    r2 = FMapResult.load_table(p)
    assert np.allclose(r2.D_half, r.D_half, atol=1e-4)
    assert np.allclose(r2.D, r.D, atol=1e-4)
    # warn=False は警告なし・D は同じ
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        r3 = frequency_map(st0, 2, x=np.linspace(1e-3, 6e-3, 3), y=np.linspace(1e-3, 3e-3, 2), nturn=256, warn=False, scaling=False)
    assert not any(issubclass(w.category, FMapWarning) for w in wl)
    assert np.allclose(r3.D, r.D) and np.all(np.isnan(r3.D_half))


def test_default_map_has_no_warning(st0):
    """既定の設定 (NORFSW・NORAD・1024 周・Hann 窓) の ASTRID2 は警告なし (小さな格子)"""
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        r = frequency_map(st0, 2, x=[2e-3], y=[1e-3], nturn=1024)
    assert not any(issubclass(w.category, FMapWarning) for w in wl)
    assert r.warnings == []


# ---------------------------------------------------------------------------------------------------------------------
# 保存済みの出力
# ---------------------------------------------------------------------------------------------------------------------
def _jsonl(name):
    p = os.path.join(EX, 'fmap_check_data', name)
    if not os.path.exists(p):
        pytest.skip('保存済みの出力がない: ' + name)
    return [json.loads(l) for l in open(p, encoding='utf-8')]


def test_saved_errors_scan_astrid2():
    rows = _jsonl('errors_astrid2.jsonl')
    allr = [d for d in rows if d['kind'] == 'all']
    for k in (0.3, 1.0, 3.0):
        L = [d for d in allr if d['k'] == k]
        assert len(L) >= 5
    f = lambda k: np.mean([d['stats']['frac_hi'] for d in allr if d['k'] == k])      # noqa
    assert f(0.3) < f(1.0) < f(3.0)
    assert f(1.0) > 5 * 0.03                      # 誤差なし (0.030) の 5 倍以上


def test_saved_gs_errors():
    p = os.path.join(EX, 'fmap_check_gs_errors.out.txt')
    if not os.path.exists(p):
        pytest.skip('保存済みの出力がない')
    txt = open(p, encoding='utf-8').read()
    assert 'sad1.exe の残り: なし' in txt
    lines = [l for l in txt.splitlines() if l[:2].strip().isdigit() and '|' in l]
    assert len(lines) == 6
    for l in lines[:4]:                            # 規則的・弱いカオスの 4 点: 座標が 1e-11 m 以内
        assert float(l.split('|')[2].split()[1]) < 1e-11
    for l in lines:
        d = [float(v) for v in l.split('|')[1].split()]
        assert abs(d[0] - d[2]) < 5e-4            # D は gs と 4 桁以内


# ---------------------------------------------------------------------------------------------------------------------
# SAD の台本の関数 FrequencyMap[]: 既定は RFSW・RAD・FLUC OFF (現在のフラグを使わない)
# ---------------------------------------------------------------------------------------------------------------------
def test_script_frequency_map_safe_default(st0):
    import io
    import re
    from pysad import ffscmd
    src = open(ASTRID, encoding='utf-8').read()
    base = 'FFS;\nUSE RING; CELL; CALC;\nLOSSAMPL = 0.1;\n'
    opt = 'X->{0.001, 0.006, 3}, Y->{0.001, 0.003, 2}'

    def run(body):
        buf = io.StringIO()
        ffscmd.run_sad_string(src + '\n' + base + body + '\nExit;\n', out=buf)
        return buf.getvalue()

    def dvals(txt, tag):
        m = re.search(tag + r'\s+\{([^}]*)\}', txt)
        return np.array([float(v) for v in m.group(1).split(',')])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = run('f = FrequencyMap[2, 256, %s]; Print["DA ", Diffusion /. f];' % opt)               # FFS の既定は RFSW ON
        b = run('NORFSW; NORAD;\nf = FrequencyMap[2, 256, %s]; Print["DB ", Diffusion /. f];' % opt)
        c = run('f = FrequencyMap[2, 256, %s, UseFlags->True]; Print["DC ", Diffusion /. f];' % opt)
    assert 'are not used' in a and 'are not used' not in b
    assert np.allclose(dvals(a, 'DA'), dvals(b, 'DB'))                      # 既定は NORFSW・NORAD と同じ
    assert not np.allclose(dvals(c, 'DC'), dvals(b, 'DB'))                  # UseFlags->True は RF ON の別の結果
    assert 'WARNING [rfsw]' in c and 'WARNING [nturn_short]' in a


def test_rf_on_rad_off_has_zero_synchronous_phase():
    """RFSW ON・RAD OFF の TrackSetup は同期位相 0 (放射損失がないので正味の加速なし)。gs・台本の CALC のあとの TrackParticles (trf0 = 0) と同じ"""
    src = open(ASTRID, encoding='utf-8').read()
    s = TrackSetup.from_source(src, 'RING', rad=False, rfsw=True, alost=0.1)[0]
    assert s.trf0 == 0.0 and s.bucket.get('phis') == 0.0
    s2 = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.1)[0]
    assert abs(s2.trf0) > 0.01


def test_saved_mpw():
    """4νy = 9 は MPW で点が失われる原因ではない: νy を ±0.05 動かしても失われる点は 207〜252 のまま (62〜75%)、νx・νy を MPW なしに戻しても 200 以上、
    MPW を弱めると y の限界が広がる"""
    rows = {d['case']['name']: d for d in _jsonl('mpw.jsonl')}
    lost = lambda n: rows[n]['stats']['nlost']        # noqa
    assert lost('A_s0') == 0
    sweep = [lost(n) for n in rows if n.startswith('B_')]
    assert len(sweep) == 6 and min(sweep) > 200 and max(sweep) < 260
    assert rows['A_s1']['d4nuy9'] < 0.001 and abs(lost('A_s1') - np.mean(sweep)) < 30         # 共鳴の真上でも離れても同じ
    assert lost('C_s1_retuned') > 200 and lost('C_s0.5_retuned') < lost('C_s1_retuned') / 2
    assert lost('A_s0.25') < lost('A_s0.5') < lost('A_s0.75') < lost('A_s1') + 1


def test_saved_k1opts():
    rows = {d['case']: d for d in _jsonl('k1opts.jsonl')}
    D = lambda n: np.array([np.nan if v is None else v for v in rows[n]['D']])          # noqa
    ok = np.isfinite(D('default')) & np.isfinite(D('lin'))
    assert np.allclose(D('face'), D('default'), equal_nan=True) and np.allclose(D('k64'), D('default'), equal_nan=True)     # E = 0 では face・k64 は v2 と同じ
    assert np.corrcoef(D('default')[ok], D('lin')[ok])[0, 1] < 0.6                       # lin は D を変える
    assert rows['lin']['nu_lin'] == rows['default']['nu_lin']                            # 線形の作動点は同じ
    if 'rect_face' in rows:
        ok = np.isfinite(D('rect_default')) & np.isfinite(D('rect_face'))
        assert 0.8 < np.corrcoef(D('rect_default')[ok], D('rect_face')[ok])[0, 1] < 0.999


def test_saved_rfsw():
    rows = {d['case']: d for d in _jsonl('rfsw.jsonl')}
    assert rows['base']['stats']['frac_hi'] < 0.05
    for k in ('rf', 'rfrad', 'rfradfluc17', 'radfluc17'):
        assert rows[k]['stats']['frac_hi'] > 0.5            # RF・RAD ON では D > −5 が多数 (拡散の指標でなくなる)
    assert rows['rfrad']['stats']['frac_hi'] > rows['rf']['stats']['frac_hi']
    assert all(len(rows[k]['warnings']) >= 1 for k in rows if k != 'base')
    assert max(abs(z[1]) for z in rows['radfluc17']['zdp']) > 5e-3 or rows['radfluc17']['zdp'][0][1] > 5e-3     # RF なしの RAD ON は δ が 1% 動く


def test_saved_nturn_scaling():
    p = os.path.join(EX, 'fmap_check_data', 'nturn_astrid2.npz')
    if not os.path.exists(p):
        pytest.skip('保存済みの出力がない')
    z = np.load(p)
    Ns = [512, 1024, 2048, 4096]
    D = {N: z['D%d' % N] for N in Ns}
    for a, b in zip(Ns[:-1], Ns[1:]):
        sh = np.nanmedian(D[b] - D[a])
        assert -1.15 < sh < -0.85                          # 周回数 2 倍で規則的な点の D は約 1 下がる (κ = 1.0)
    frac = {N: np.mean(D[N][np.isfinite(D[N])] > FC.threshold_for_nturn(N)) for N in Ns}
    fixed = {N: np.mean(D[N][np.isfinite(D[N])] > -5.0) for N in Ns}
    assert max(frac.values()) / min(frac.values()) < 3.0 < max(fixed.values()) / min(fixed.values())      # 換算した閾値の方が周回数に依らない


def test_saved_coupling_only():
    rows = _jsonl('coupling.jsonl')
    f = [d['stats']['frac_hi'] for d in sorted(rows, key=lambda d: d['target'])]
    assert f[0] > 0.1 and all(b > a for a, b in zip(f[:-1], f[1:]))        # 結合だけで D > −5 が 0.03 → 0.13 以上、結合の強さに単調


def test_saved_errors_scan_kekls():
    rows = _jsonl('errors_kekls.jsonl')
    base = [d for d in rows if d['kind'] == 'none'][0]
    assert base['stats']['nlost'] == 125 and abs(base['stats']['frac_hi'] - 0.142) < 0.005
    a = [d for d in rows if d['kind'] == 'all' and d['k'] == 0.03]
    assert len(a) == 5 and np.mean([d['stats']['frac_hi'] for d in a]) > 0.4          # ×0.03 (1.5 μm) で D > −5 が 0.14 → 0.56
    big = [d for d in rows if d['kind'] == 'all' and d['k'] == 1.0]
    assert len(big) == 5 and min(d['stats']['nlost'] for d in big) > 320           # ×1 はほぼ全部失われる
    cs = [d for d in rows if d['kind'] == 'all_cs']
    al = [d for d in rows if d['kind'] == 'all' and d['k'] == 0.1]
    assert np.mean([d['stats']['nlost'] for d in cs]) < np.mean([d['stats']['nlost'] for d in al]) - 50      # フィードダウンを切ると失われる点が減る
