# -*- coding: utf-8 -*-
"""v2.3.0 S3: 共鳴線の予想の作り直し (pysad/fmappred.py)。保存済みの周波数マップの表 (example/14_nonlinear_dynamics/frequency_map/) を読むので速い。"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '..', '..', 'example', '14_nonlinear_dynamics', 'frequency_map'))
sys.path.insert(0, EX)

from pysad import compat  # noqa: E402
from pysad import fmapres as FR  # noqa: E402
from pysad import fmappred as FP  # noqa: E402
from pysad.freqmap import FMapResult  # noqa: E402
from pysad.parser import load_sad_string  # noqa: E402


@pytest.fixture(scope='module')
def v2():
    old = compat.get_mode() if hasattr(compat, 'get_mode') else None
    compat.set_mode('v2')
    yield
    if old:
        compat.set_mode(old)


def _case(name):
    import fmap_rdt as S1
    cs = S1.CASES[name]
    lat = load_sad_string(S1.source(name)).build_lattice(cs['line'])
    rows = lat.twiss_table()
    tab = cs['table'] or os.path.join(EX, 'fmap_rdt_%s.table.txt' % name)
    return lat, rows, cs['loc'], FMapResult.load_table(tab)


@pytest.fixture(scope='module')
def fodo(v2):
    lat, rows, loc, r = _case('fodo')
    return FP.ResonancePredictorS3.from_lattice(lat, loc, rows), r, lat, rows, loc


@pytest.fixture(scope='module')
def astrid(v2):
    lat, rows, loc, r = _case('astrid2')
    return FP.ResonancePredictorS3.from_lattice(lat, loc, rows), r


def test_auc_and_spearman_synthetic():
    s = np.array([1.0, 2, 3, 4, 5, 6])
    pos = np.array([0, 0, 0, 1, 1, 1], dtype=bool)
    assert FP.auc_score(s, pos) == 1.0
    assert FP.auc_score(-s, pos) == 0.0
    assert FP.auc_score(np.ones(6), pos) == 0.5
    assert abs(FP.spearman(s, s ** 3) - 1.0) < 1e-12
    assert abs(FP.spearman(s, -s) + 1.0) < 1e-12
    rng = np.random.default_rng(1)
    a = rng.normal(size=200)
    b = a + 0.3 * rng.normal(size=200)
    c = a.copy()
    # c が a と同じなら、c を除いた偏相関は定義できるが、a だけの効果は消える
    assert abs(FP.partial_spearman(a, b, c + 1e-9 * rng.normal(size=200))) < 0.2


def test_score_stats_top_decile():
    D = np.linspace(-10, -2, 100)
    st = FP.score_stats(D, D, thr=-5.0, top=0.1)
    assert st['auc'] == 1.0 and st['top_Dmed'] > -3 and st['top_frac_pos'] == 1.0


def test_harmonics_and_parity(fodo):
    pr = fodo[0]
    assert (3, 0) in pr.harms and (2, -2) in pr.harms
    # 3 次 (奇数) と 4 次 (偶数) は同じ調和に同居しない
    for a, ts in pr.harms.items():
        assert len({j + k + l + m for (j, k, l, m, _h) in ts}) == 1
    # 3 次の h は resonance_terms と同じ
    from pysad import resdrive as R
    r = R.resonance_terms(fodo[2], ['h30000'], rows=fodo[3])
    h30 = [t for t in pr.harms[(3, 0)] if (t[0], t[1], t[2], t[3]) == (3, 0, 0, 0)][0][4]
    assert abs(h30 - r.total['h30000']) < 1e-9 * abs(h30)


def test_no_tracking_fit_in_predictor(fodo):
    """方針: 追跡に合わせた当てはめの項を持たない (ν(J) は解析の線形の ADTS)"""
    pr = fodo[0]
    assert not hasattr(pr, 'fit_tunes')
    Jx = np.array([1e-7, 2e-7])
    Jy = np.array([3e-8, 1e-8])
    nx, ny = pr.tunes_model(Jx, Jy)
    assert np.allclose(nx, pr.nu0[0] + pr.A[0, 0] * Jx + pr.A[0, 1] * Jy)


def test_fodo_3nux_island_discriminates(fodo):
    pr, r = fodo[0], fodo[1]
    m = pr.predict(r.x, r.y)
    e = FP.evaluate_prediction(r, m, nshift=100)
    assert e['cov']['cov'] > 0.99                     # D > -5 の 44 点がすべて 3νx の包絡に入る (S1 と同じ)
    assert e['score']['auc'] > 0.9
    assert e['score']['rho_partial'] > 0.4             # 振幅の効果を除いても順位相関がある
    assert e['score']['top_frac_pos'] > 5 * e['score']['base_frac_pos']
    # 3νx が作動点の 0.0089 以内 → 警告
    assert any('1 次の RDT の共鳴線の近く' in w for w in pr.warnings_for(m))


def test_predict_matches_s1_for_first_order_only(fodo):
    """FODO は 3νx だけが包絡に入る。S1 の予想 (1 次だけ、同じ ADTS 行列) と同じ包絡 (A が同じなら)"""
    pr, r, lat, rows, loc = fodo
    s1 = FR.ResonancePredictor.from_lattice(lat, loc, rows)
    s1.A = pr.A.copy()
    a = s1.predict(r.x, r.y)['envelope']
    b = pr.predict(r.x, r.y)['envelope']
    assert (a & b).sum() >= 0.95 * a.sum()


def test_model_hook(fodo):
    pr, r, lat, rows, loc = fodo
    p2 = FR.ResonancePredictor.from_lattice(lat, loc, rows, model='s3')
    assert isinstance(p2, FP.ResonancePredictorS3)


def test_normal_form_correction_changes_off_diagonal(astrid):
    pr, r = astrid
    m0 = pr.predict(r.x, r.y, nf=False)
    m1 = pr.predict(r.x, r.y, nf=True)
    W0 = m0['by_harmonic'][(2, 0)]['W']
    W1 = m1['by_harmonic'][(2, 0)]['W']
    assert np.nanmedian(W1 / W0) > 1.5                # h20110 の補正 (4.6 → 91) で (2,0) の幅が大きく変わる
    # 補正のない調和 (3 次) は同じ
    assert np.allclose(m0['by_harmonic'][(3, 0)]['W'], m1['by_harmonic'][(3, 0)]['W'], equal_nan=True)


def test_astrid2_island_locking_2nux_2nuy(astrid):
    """D を使わない島の予想の検証: ASTRID2 の 2νx−2νy (νx−νy = 3) の予想の島の中の点は、測ったチューンが線に固着している"""
    pr, r = astrid
    m = pr.predict(r.x, r.y)
    v = m['by_harmonic'][(2, -2)]
    ip = (np.floor(pr.nu0[0]), np.floor(pr.nu0[1]))
    nx, ny = r.nux1 + ip[0], r.nuy1 + ip[1]
    d0m = ((2 * nx - 2 * ny) - np.round(2 * nx - 2 * ny)) / 2
    d0 = v['delta0']
    inn = v['trapped']
    ring = ~v['trapped'] & (np.abs(d0) >= v['W']) & (np.abs(d0) < 3 * v['W'])
    lock = (np.abs(d0m) < 0.003) & (np.abs(d0m) < 0.25 * np.abs(d0))
    assert inn.sum() > 50 and lock[inn].mean() > 0.8 and lock[ring].mean() < 0.4


def test_reliability_flags(astrid):
    pr, r = astrid
    m = pr.predict(r.x, r.y)
    rel = m['reliability']
    assert set(np.unique(rel)) <= {0, 1, 2}
    # 重なり (2 つ以上の包絡) は信頼できない
    assert (rel[m['n_env'] >= 2] == 0).all()
    # 線形の移動が 0.05 を超える点は信頼できない
    Jx, Jy = m['Jx'], m['Jy']
    lin = pr.tunes_lin(Jx, Jy)
    far = np.hypot(lin[0] - pr.nu0[0], lin[1] - pr.nu0[1]) > 0.05
    assert (rel[far] == 0).all()
    # 最小振幅の点は信頼できる側 (0 でない)
    assert rel.ravel()[0] >= 1
    assert any('線形の ADTS' in w for w in pr.warnings_for(m))


def test_nearest_order_stats_synthetic():
    nux = np.array([0.5, 0.3333, 0.25, 0.123])
    nuy = np.full(4, 0.1234)
    D = np.array([-4.0, -4.0, -8.0, -8.0])
    o = FP.nearest_order_stats(D, nux, nuy, tol=0.003, orders=(1, 2, 3, 4))
    # 0.5: 2νx = 1 (次数 2)、0.3333: 3νx = 1 (次数 3)
    assert o['frac_hi'][1] == 0.5 and o['frac_hi'][2] == 0.5


def test_delta_unsupported(fodo):
    pr, r, lat, rows, loc = fodo
    with pytest.raises(NotImplementedError):
        FP.ResonancePredictorS3.from_lattice(lat, loc, rows, delta=0.001)
