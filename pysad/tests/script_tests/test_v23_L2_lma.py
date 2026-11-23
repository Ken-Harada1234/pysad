# -*- coding: utf-8 -*-
"""v2.3.0 L2: 局所運動量アクセプタンス (LMA) を信頼して使えるようにする変更の試験。
(1) check_every の既定 = 1 (要素ごとの判定。0 は SAD の tapert のままで、アパーチャが限界の位置があると警告)、(2) RF バケットの限界の独立な確認 (解析式・1 次元の写像・
6 次元の追跡) と、バケット・アパーチャの区別 (limit_types)、(3) verify (狭い損失の窓)、(4) turncheck (周回数への依存)、(8) 警告 (アパーチャ未設定・FLUC・周期性の崩れ・走査の上限・周回数)、
periodicity_report (誤差つきで周期性が崩れることを検出する)。保存した例の表の整合 (ASTRID2・kek-ls)。
"""
import os
import warnings

import numpy as np
import pytest

from pysad import compat, lma, lmabucket as lb
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
EX = os.path.join(ROOT, 'example', '15_lma_touschek', 'lma')

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


# ---------------------------------------------------------------------------------------------------------------------
# (1) check_every の既定
# ---------------------------------------------------------------------------------------------------------------------

def test_default_check_every_is_one_and_zero_warns(fodo_rf):
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    kw = dict(workers=1, tol=2e-4, turncheck=0)
    a = lma.momentum_acceptance(st, [5], quiet=True, **kw)
    b = lma.momentum_acceptance(st, [5], check_every=1, quiet=True, **kw)
    assert a.cond['check_every'] == 1 and np.array_equal(a.dp_plus, b.dp_plus) and np.array_equal(a.dp_minus, b.dp_minus)
    with pytest.warns(UserWarning, match='check_every = 0'):
        z = lma.momentum_acceptance(st, [5], check_every=0, **kw)
    assert z.cond['check_every'] == 0 and z.dp_plus[0] >= a.dp_plus[0] - 1.5e-4
    # setup 側で check_every を決めておけばそれに従う (None のとき)
    c = lma.momentum_acceptance(st.copy(check_every=5), [5], quiet=True, **kw)
    assert c.cond['check_every'] == 5


# ---------------------------------------------------------------------------------------------------------------------
# (2) RF バケット
# ---------------------------------------------------------------------------------------------------------------------

def test_rf_bucket_analytic_and_longitudinal_params(fodo_rf):
    st, ffs = fodo_rf
    an = lb.rf_bucket_analytic(st)
    assert abs(an['delta_rf'] / an['emit_bucket_height'] - 1) < 2e-3          # 解析式 = Emittance[] の bucket_height (2 cos φs − (π − 2φs) sin φs の式)
    assert abs(an['h'] - an['k_rf'] * st.circ / (2 * np.pi)) < 1e-9
    par = lb.longitudinal_params(st)
    assert par['kick_max_diff'] < 1e-3 * par['kick_scale']                    # 空洞の力 = (V/E0) sin(φs + k z) − U0/E0 (追跡と 1e-3 以内)
    assert abs(par['a1_over_alpha'] - 1) < 0.1                                # 測った a1 は αc に近い (α2 の分だけ違う)
    assert par['r1'] < 0                                                      # 放射の減衰 (δ が大きいほど失う)


def test_bucket_1d_map_matches_6d_tracking(fodo_rf):
    """縦方向の 1 次元の写像 (解析式の RF + 測った α1・α2・放射) の RF バケットの限界は、6 次元の追跡 (アパーチャ無限大) と 3e-4 以内 (または 3%) で一致する"""
    st, ffs = fodo_rf
    nt = st.nturn_for_synchrotron(2.0)
    m = lb.bucket_limits_1d(st, nt)
    r = lb.bucket_limits_6d(st, [2], nturn=nt, tol=5e-5)
    bh = st.bucket['height']
    assert abs(m['dp_plus'] - r.dp_plus[0]) < 0.03 * bh
    assert abs(m['dp_minus'] - r.dp_minus[0]) < 0.03 * bh
    # 周回数が短いと限界は上がる (バケットを出るのに時間がかかる)
    m1 = lb.bucket_limits_1d(st, int(0.5 / abs(st.bucket['nus'])))
    assert m1['dp_plus'] is None or m1['dp_plus'] > m['dp_plus']


def test_limit_types_separates_rf_and_aperture(fodo_rf):
    st, ffs = fodo_rf
    wide = st.copy(alost=0.1)
    r = lma.momentum_acceptance(wide, [3], workers=1, tol=1e-4, quiet=True, turncheck=0)
    lt = lb.limit_types(wide, r, nturn_aper=64)
    assert lt.type_plus[0] in ('rf', 'both') and lt.type_minus[0] in ('rf', 'both')
    assert abs(lt.dev_plus[0]) < 3e-4 and abs(lt.dev_minus[0]) < 3e-4
    narrow = st.copy(alost=0.004)
    r = lma.momentum_acceptance(narrow, [5], workers=1, tol=1e-4, quiet=True, turncheck=0)
    lt = lb.limit_types(narrow, r, nturn_aper=64)
    assert lt.type_plus[0] == 'aper' and lt.type_minus[0] == 'aper'
    assert abs(lt.dev_plus[0]) < 4e-4 and abs(lt.dev_minus[0]) < 4e-4         # LMA = min(rf, aper)
    assert 'aper' in lt.table()


# ---------------------------------------------------------------------------------------------------------------------
# (3) verify・(4) turncheck: 模擬の追跡 (損失の判定を差し替える) で窓・周回数の依存を確かめる
# ---------------------------------------------------------------------------------------------------------------------

def _fake_track(window=None, slow_above=None, top=0.010):
    def f(setup, loc, cod, pts, nturn, zb, x0, y0, alost):
        out = {}
        for sg, d in pts:
            lost = d > top
            if window and window[0] < d < window[1]:
                lost = True
            if slow_above is not None and nturn > 600 and d > slow_above:
                lost = True
            out[(sg, d)] = (7 if lost else 0, 'aper')
        return out
    return f


def test_verify_finds_a_narrow_window(monkeypatch):
    monkeypatch.setattr(lma, '_track_points', _fake_track(window=(0.0033, 0.0040)))
    args = (None, 1, np.zeros(6), 500, 0.0125, 12, 3, 1e-4, None, 0.0, 0.0)
    plain = lma._lma_location(*args, verify=0)
    assert abs(plain['dp_plus'] - 0.010) < 2e-4 and plain['nwin_plus'] == 0                       # 窓を見逃す (粗い走査の刻み 1.0e-3 のあいだ)
    ver = lma._lma_location(*args, verify=100)
    assert ver['nwin_plus'] > 0 and ver['nwin_minus'] > 0
    assert 0.0033 < ver['win_plus'] < 0.0041 and ver['dp_plus'] < 0.0034                             # 限界が窓の手前に下がる


def test_turncheck_flags_turn_dependence(monkeypatch):
    monkeypatch.setattr(lma, '_track_points', _fake_track(slow_above=0.0085))
    args = (None, 1, np.zeros(6), 500, 0.0125, 12, 3, 1e-4, None, 0.0, 0.0)
    r = lma._lma_location(*args, turncheck=2.0)
    assert r['tdep_plus'] and r['tdep_minus']
    r = lma._lma_location(*args, turncheck=0)
    assert not r['tdep_plus']


def test_flags_and_table_columns(fodo_rf):
    st, ffs = fodo_rf
    st = st.copy(alost=0.004)
    r = lma.momentum_acceptance(st, [5], workers=1, tol=2e-4, quiet=True)
    t = r.table().splitlines()
    assert t[-1].split()[-1] == r.flags(0)
    assert 'flags' in [ln for ln in t if ln.startswith('# loc')][0]
    # 保存した表の読み込み (T1 の比較の読み方: 先頭の 9 列) と互換
    cols = t[-1].split()
    assert int(cols[0]) == 5 and abs(float(cols[3]) - r.dp_plus[0]) < 1e-8


# ---------------------------------------------------------------------------------------------------------------------
# (8) 警告・周期性
# ---------------------------------------------------------------------------------------------------------------------

def test_warnings_for_unusable_conditions(fodo_rf):
    st, ffs = fodo_rf
    # アパーチャ未設定 (LOSSAMPL = 1 m) と、シンクロトロン振動 1 周期未満
    with pytest.warns(UserWarning) as rec:
        r = lma.momentum_acceptance(st.copy(alost=1.0), [3], nturn=20, workers=1, ncoarse=2, nsub=1, tol=0.01, turncheck=0)
    msg = ' '.join(str(w.message) for w in rec)
    assert 'LOSSAMPL' in msg and 'shorter than one synchrotron period' in msg
    assert r.cond['warnings'] and r.warnings == r.cond['warnings']
    # FLUC
    with pytest.warns(UserWarning, match='FLUC'):
        lma.momentum_acceptance(st.copy(alost=0.004, fluc=True), [3], nturn=300, workers=1, ncoarse=2, nsub=1, tol=0.01, turncheck=0)
    # 走査の上限まで残る
    r = lma.momentum_acceptance(st.copy(alost=0.1), [3], nturn=300, dpmax=0.001, workers=1, ncoarse=2, nsub=1, tol=0.01, turncheck=0, quiet=True)
    assert any('scan top' in w for w in r.warnings) and r.reason_plus[0] == 'sat'


def test_periodicity_report_detects_errors_and_mpw_like_breaks(fodo_rf):
    from pysad import fmapres
    st, ffs = fodo_rf
    rep = lma.periodicity_report(st, warn=False)
    assert rep['periodic'] and rep['nper'] == 24 and rep['cod_spread'] < rep['cod_tol']
    assert abs(rep['delta_sawtooth'] - st.bucket['U0'] / st.p0) < 1e-15
    locs, rep = lma.cell_locations(st)
    assert len(locs) < 20 and rep['periodic']
    st2, cnt = fmapres.apply_alignment_errors(st, seed=1, quad_dxy=50e-6, quad_rot=0.2e-3)
    with pytest.warns(UserWarning, match='not periodic'):
        locs2, rep2 = lma.cell_locations(st2)
    assert not rep2['periodic'] and len(locs2) > 100
    # 要素の並びが同じでも、閉軌道がセルで違えば (周期性の崩れ) 検出する: 偏向磁石だけの誤差のように、find_periodicity が見ない場合の想定
    rep3 = lma.periodicity_report(st, cod_tol=0.0, warn=False)             # 許容 0: 数値の揺れでも「周期的でない」にする
    assert (not rep3['periodic']) or rep3['cod_spread'] == 0.0


# ---------------------------------------------------------------------------------------------------------------------
# 保存した例の表 (ASTRID2・kek-ls)
# ---------------------------------------------------------------------------------------------------------------------

def _read(path):
    rows = []
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        rows.append(ln.split())
    return rows


@pytest.mark.parametrize('name,n,lo', [('astrid2', 25, ('aper',)), ('kekls', 11, ('rf', 'both'))])
def test_saved_limit_types(name, n, lo):
    """保存した limit_types の表 (ASTRID2 25 位置、kek-ls 11 位置): LMA = min(RF バケットの限界 (6 次元)、実効のアパーチャの限界 (4 次元)) が 4e-4 以内で成り立つ。
    ASTRID2 (±10 mm) は全位置がアパーチャ、kek-ls は全位置が RF バケット"""
    p = os.path.join(EX, 'lma_limit_types_%s.table.txt' % name)
    if not os.path.exists(p):
        pytest.skip('例の表がない')
    rows = _read(p)
    assert len(rows) == n
    # 列: loc s name LMA+ rf+ aper+ type+ reason+ | LMA- rf- aper- type- reason- | aper_eff
    for t in rows:
        lma_p, rf_p, eff = float(t[3]), float(t[4]), float(t[15])
        lma_m, rf_m = float(t[9]), float(t[10])
        assert abs(lma_p - min(rf_p, eff)) < 4e-4 and abs(-lma_m - min(-rf_m, eff)) < 4e-4
        assert t[6] in lo and t[12] in lo
