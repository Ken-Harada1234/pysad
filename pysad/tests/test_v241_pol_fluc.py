# -*- coding: utf-8 -*-
"""v2.4.1: 偏光 (POL) の FLUC つきスピン追跡の gs との比較 (4 つのリング: 12 セル FODO の 6 MV・24 MV、PF、ASTRID2)。

台本は `scripts/v241_pol_*.sad`、gs の参照出力は同名の `.gs.txt`、pysad (k64) の出力は `.pysad.txt`。生成・解析は `_v241_pol_tools.py`。
- ident (同一性): gs と pysad (k64) に同じ台本 (粒子ごとの追跡、乱数 SeedRandom・NPARA=1) を渡し、数値が一致する。
- tune (スピンチューン): 閉軌道上の 1 粒子の毎周のスピンの回転角から求めたスピンチューンと、Emittance[] の SpinTune。
- stat (統計): gs の大粒子数の追跡 (台本・出力を保存) の脱偏極率 (1 − <S·n0> の周あたりの増え) とその統計誤差。
  pysad の式 (v2 の Polarization の `depol_rate`) との比は docs/46 と docs/43 の表。長い計算は `PYSAD_SLOW=1`。
記録は docs/history/v241_pol_fluc.md。
"""
import io
import os
import sys

import numpy as np
import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _v241_pol_tools as T  # noqa: E402
SC = os.path.join(HERE, 'scripts')
SLOW = pytest.mark.skipif(not os.environ.get('PYSAD_SLOW'), reason='pysad の Python 追跡 (数分) か v2 の Emittance[] (POL ON、FODO で約 100 秒)')


def _p(name):
    return os.path.join(SC, name)


def _have(*names):
    for n in names:
        if not os.path.exists(_p(n)):
            pytest.skip(n + ' がない')


def _run_pysad(sad):
    out = io.StringIO()
    ffscmd.run_sad_string(open(_p(sad), encoding='utf-8').read(), out=out)
    tmp = _p(sad.replace('.sad', '._tmp_out.txt'))
    with open(tmp, 'w', encoding='utf-8') as fh:
        fh.write(out.getvalue())
    try:
        return T.read_stat(tmp)
    finally:
        os.remove(tmp)


def _same(a, b, rtol):
    """追跡の量 (STAGE・DAMP・DV・AV・LV・EQ) は rtol。Emittance[] の量は、閉軌道が収束しにくい ASTRID2 (y 12 mm) で gs と 1e-6 ほど違う
    (SpinTune 4e-7、n0 の成分 2e-5。追跡は全桁一致) ので、時間・チューンは 1e-5、n0 は絶対 1e-4 で比べる"""
    assert a['EM'] is not None and len(a['EM']) == len(b['EM'])
    for i in (0, 1, 2):
        assert abs(a['EM'][i] - b['EM'][i]) <= 1e-5 * abs(a['EM'][i]), ('EM', i)
    for x, y in zip(a['EM'][4:], b['EM'][4:]):
        assert abs(x - y) <= 1e-4, ('EM n0', x, y)
    assert len(a['EQ']) == len(b['EQ'])
    for x, y in zip(a['EQ'], b['EQ']):
        assert abs(x - y) <= rtol * max(abs(x), 1e-30), ('EQ', x, y)
    assert len(a['stages']) == len(b['stages']) > 0
    for sa, sb in zip(a['stages'], b['stages']):
        for k in ('A', 'L', 'rms_x', 'rms_y', 'rms_dp', 'rms_z'):
            assert abs(sa[k] - sb[k]) <= rtol * max(abs(sa[k]), 1e-30), (k, sa[k], sb[k])
    for sa, sb in zip(a['damp'], b['damp']):
        assert abs(sa['D'] - sb['D']) <= rtol * abs(sa['D']), (sa, sb)
    for tag in ('DV', 'AV', 'LV'):
        for k in a[tag]:
            assert np.allclose(a[tag][k], b[tag][k], rtol=rtol, atol=0), (tag, k)


# ---------------------------------------------------------------------------------------------------------------------
# 同一性 (小さな台本を pysad (k64) で走らせて、保存してある gs の出力と比べる)

def test_ident_fodo6_pysad_k64_matches_gs_per_particle(k64_mode):
    """4.2 GeV・6 MV の 12 セル FODO (縦ずれ ±3 mm)、3 粒子、RAD FLUC + POL、スピンの 1 − <S·n0> が粒子ごとに gs と全桁近く一致 (相対 1e-6 以下)"""
    _have('v241_pol_identq_fodo6.sad', 'v241_pol_identq_fodo6.gs.txt')
    g = T.read_stat(_p('v241_pol_identq_fodo6.gs.txt'))
    b = _run_pysad('v241_pol_identq_fodo6.sad')
    _same(g, b, 1e-6)


@pytest.mark.parametrize('ring', ['fodo24', 'pf', pytest.param('astrid2', marks=SLOW)])
def test_ident_other_rings_pysad_k64_matches_gs(ring, k64_mode):
    """24 MV の FODO・PF・ASTRID2 (k64 の閉軌道が収束する縦ずれ)"""
    _have('v241_pol_identq_%s.sad' % ring, 'v241_pol_identq_%s.gs.txt' % ring)
    g = T.read_stat(_p('v241_pol_identq_%s.gs.txt' % ring))
    b = _run_pysad('v241_pol_identq_%s.sad' % ring)
    _same(g, b, 1e-6)


@pytest.mark.parametrize('ring', ['fodo6', 'fodo24', 'pf', 'astrid2'])
def test_ident_stored_runs_gs_vs_pysad_agree(ring):
    """保存した大きめの同一性の実行 (gs と pysad k64、10 粒子・数千周ほか) の gs と pysad の出力が一致する"""
    _have('v241_pol_ident_%s.gs.txt' % ring, 'v241_pol_ident_%s.pysad.txt' % ring)
    g = T.read_stat(_p('v241_pol_ident_%s.gs.txt' % ring))
    b = T.read_stat(_p('v241_pol_ident_%s.pysad.txt' % ring))
    _same(g, b, 1e-6)


# ---------------------------------------------------------------------------------------------------------------------
# スピンチューン (追跡) と Emittance[] の SpinTune

# gs の SpinTune (Emittance[]、k64)、gs の追跡 (毎周のスピンの回転角の直線当てはめ、n0 まわり)、v2 の Emittance[] の nup
TUNE = {
    #            gs EM        gs 追跡 |ν|   v2 EM |ν|
    'fodo6': (0.23625181, 0.23587558, 0.23589590),
    'fodo24': (0.23625181, 0.23587557, 0.23589590),
    'pf': (0.32765303, 0.32673071, 0.32662037),
    'astrid2': (0.31275848, 0.31752023, 0.31490465),
}


@pytest.mark.parametrize('ring', ['fodo6', 'fodo24', 'pf', 'astrid2'])
def test_spin_tune_from_gs_tracking_vs_emittance(ring):
    """gs の追跡のスピンチューン (毎周のスピン) と Emittance[] の SpinTune を比べる (保存した gs の出力)。
    k64 の Emittance[] は追跡と 3.8e-4 (FODO、閉軌道 y 6 mm) 〜 5e-3 (ASTRID2、y 12 mm) 違う。原因は B58 (Emittance[] が偏向磁石の出口の端の縦磁場 ψ を捨てる)。"""
    _have('v241_pol_tune_%s.gs.txt' % ring, 'v241_pol_prep_%s.json' % ring)
    import json
    n0 = json.load(open(_p('v241_pol_prep_%s.json' % ring)))['n0']
    em, S = T.read_tune(_p('v241_pol_tune_%s.gs.txt' % ring))
    nu, res, amp, par = T.tune_from_spin(S, n0)
    ref = TUNE[ring]
    assert abs(abs(em[1]) - ref[0]) < 2e-8
    assert abs(abs(nu) - ref[1]) < 2e-8
    assert amp == pytest.approx(1.0, abs=1e-6) and abs(par) < 1e-3
    # k64 の Emittance[] の SpinTune は追跡より大きく外れる (v2 の Emittance[] は追跡に近い: FODO 2e-6、PF 1e-4、ASTRID2 は 2.6e-3 (別の原因: 追跡側の刻み、docs/46))
    assert abs(abs(em[1]) - abs(nu)) > 3e-4
    if ring in ('fodo6', 'fodo24'):
        assert abs(ref[2] - abs(nu)) < 5e-5


def test_spin_tune_tracking_pysad_k64_equals_gs_fodo6(k64_mode):
    """pysad (k64) の追跡の毎周のスピンが gs と 1e-8 以下で一致 (台本は同じ、閉軌道上の粒子)"""
    _have('v241_pol_tune_fodo6.sad', 'v241_pol_tune_fodo6.gs.txt')
    out = io.StringIO()
    ffscmd.run_sad_string(open(_p('v241_pol_tune_fodo6.sad'), encoding='utf-8').read(), out=out)
    tmp = _p('v241_pol_tune_fodo6._tmp_out.txt')
    with open(tmp, 'w', encoding='utf-8') as fh:
        fh.write(out.getvalue())
    try:
        em_b, S_b = T.read_tune(tmp)
    finally:
        os.remove(tmp)
    em_g, S_g = T.read_tune(_p('v241_pol_tune_fodo6.gs.txt'))
    assert np.allclose(em_b, em_g, rtol=1e-8)
    assert np.abs(S_b - S_g).max() < 1e-8


# ---------------------------------------------------------------------------------------------------------------------
# 統計 (gs の保存した出力から脱偏極率と統計誤差)

def _stat_rate(fn, tag):
    o = T.read_stat(_p(fn))
    return o, T.rate(o, tag)


# 台本ごとの測定値 (gs の追跡 k64、D = 減衰してから測った 1 − <S·n0>。粒子ごとの直線当てはめの傾きの平均と標準誤差 [/周])
STAT = {
    'fodo6_s101': ('v241_pol_stat_fodo6_s101.gs.txt', 'DV', 1.128112e-07, 1.30848e-08),
    'fodo6_s102': ('v241_pol_stat_fodo6_s102.gs.txt', 'DV', 1.843274e-07, 7.58306e-08),
    'fodo24_s201': ('v241_pol_stat_fodo24_s201.gs.txt', 'DV', 7.121243e-08, 2.99986e-09),
    'pf_s301': ('v241_pol_stat_pf_s301.gs.txt', 'DV', 2.379943e-09, 2.41411e-10),
    # ASTRID2 (40 粒子 × 18000 周): 率は 0 と区別できない (1.5σ)。式 5.4e-14 の 1 σ の感度の 1/1900。下の「検出限界」の試験
    'astrid2_s402': ('v241_pol_stat_astrid2_s402.gs.txt', 'DV', 1.036588e-10, 6.82535e-11),
}
# v2 の `Polarization` の脱偏極率 (側帯・クロマティックな変調込み。prep.json に保存)
RING_OF = {'fodo6_s101': 'fodo6', 'fodo6_s102': 'fodo6', 'fodo24_s201': 'fodo24', 'pf_s301': 'pf', 'astrid2_s402': 'astrid2'}


@pytest.mark.parametrize('name', sorted(STAT))
def test_stat_gs_depolarization_rate_regression(name):
    fn, tag, val, err = STAT[name]
    _have(fn)
    o, r = _stat_rate(fn, tag)
    assert r[0] == pytest.approx(val, rel=1e-3)
    assert r[1] == pytest.approx(err, rel=1e-2)
    # D と L (線形の n(u)) は 5% 以内で、6 MV の FODO の 2 つ目の種だけ重い裾のため誤差が大きい
    if name != 'astrid2_s402':      # ASTRID2 は率が 0 と区別できず D と L が一致しない (下の検出限界の試験)
        assert T.rate(o, 'LV')[0] == pytest.approx(val, rel=0.1)


def test_stat_formula_over_tracking_ratio():
    """v2 の式 / gs の追跡の脱偏極率。24 MV の FODO と PF は統計 (4%・10%) の範囲で 1 (1.05 ± 0.04、0.98 ± 0.10)。
    6 MV の FODO は式が 0.4〜0.65 倍だが、追跡の平均が少数の粒子に支配される (下の裾の試験)。"""
    import json
    for name, lo, hi in (('fodo24_s201', 0.95, 1.15), ('pf_s301', 0.8, 1.2), ('fodo6_s101', 0.45, 0.85), ('fodo6_s102', 0.2, 0.7)):
        _have(STAT[name][0], 'v241_pol_prep_%s.json' % RING_OF[name])
        f = json.load(open(_p('v241_pol_prep_%s.json' % RING_OF[name])))['v2']['depol_rate']
        r = f / T.rate(T.read_stat(_p(STAT[name][0])), 'DV')[0]
        assert lo < r < hi, (name, r)


def test_stat_astrid2_rate_is_below_detection_limit():
    """ASTRID2 (0.58 GeV、40 粒子 × 18000 周、gs の追跡 k64): 1 − <S·n0> は (2〜5)e-6 で時間に対して平らで、率は D (1.04 ± 0.68)e-10、
    A・L は負 (−1.2e-10 ± 0.4e-10)。v2 の式 5.4e-14 /周 (18000 周で 1e-9 の増え) は 1 − <S·n0> の平均の標準誤差 (4e-7〜8e-7) の 1/500 で、**式の検証にならない**。
    比 (式 / 追跡) は意味を持たないので出さない。追跡の率の符号が D と A・L で逆であることが、統計の揺らぎの大きさの証拠。"""
    import json
    fn = STAT['astrid2_s402'][0]
    _have(fn, 'v241_pol_prep_astrid2.json')
    o = T.read_stat(_p(fn))
    assert len(o['DV'][1]) == 40 and [s['turns'] for s in o['stages']] == [6000.0, 12000.0, 18000.0]
    f = json.load(open(_p('v241_pol_prep_astrid2.json')))['v2']['depol_rate']
    assert f == pytest.approx(5.44e-14, rel=0.01)
    rd, ed = T.rate(o, 'DV')[:2]
    ra, ea = T.rate(o, 'AV')[:2]
    assert abs(rd) < 2.5 * ed                       # 0 と区別できない
    assert ra < 0 < rd                              # 符号が揃わない = 雑音
    assert f * 18000 < np.mean([s['sdA'] for s in o['stages']]) / 150   # 式の 18000 周の増えは平均の標準誤差の 1/150 以下
    assert min(ed, ea) > 500 * f                    # 1σ の感度が式の 500 倍以上悪い (実測 690〜1250 倍)
    for s in o['stages']:                           # 平衡の分布は保たれている (σδ は Emittance[] の 4.24e-4 の ±15%)
        assert 3.6e-4 < s['rms_dp'] < 4.9e-4


def test_stat_fodo6_rate_is_dominated_by_rare_particles():
    """6 MV の FODO: 粒子ごとの傾きの中央値 5e-8 に対し平均 1.5e-7、2000 粒子のうち 1 粒子が合計の 25%、上位 1% を除くと 8.8e-8。
    24 MV・PF では最大の 1 粒子が合計の 1〜5% しか占めない。式 (線形の軌道のガウス分布) は裾の共鳴を含まない。"""
    def slopes(fn):
        o = T.read_stat(_p(fn))
        stg = {s['k']: s['turns'] for s in o['stages']}
        ks = sorted(o['DV'])
        tt = np.array([stg[k] for k in ks], float)
        Y = np.array([o['DV'][k] for k in ks])
        tm = tt.mean()
        return ((tt - tm)[:, None] * (Y - Y.mean(0))).sum(0) / ((tt - tm) ** 2).sum()
    _have('v241_pol_stat_fodo6_s101.gs.txt', 'v241_pol_stat_fodo6_s102.gs.txt', 'v241_pol_stat_fodo24_s201.gs.txt')
    s6 = np.concatenate([slopes('v241_pol_stat_fodo6_s101.gs.txt'), slopes('v241_pol_stat_fodo6_s102.gs.txt')])
    s24 = slopes('v241_pol_stat_fodo24_s201.gs.txt')
    assert np.median(s6) == pytest.approx(4.96e-8, rel=0.02)
    assert s6.max() / s6.sum() > 0.15
    assert s24.max() / s24.sum() < 0.03
    q = np.sort(s6)
    assert q[:-len(s6) // 100].mean() == pytest.approx(8.8e-8, rel=0.02)


# ---------------------------------------------------------------------------------------------------------------------
# B270 (pysad 自身の誤り): 台本の Emittance[] (emit_extra の経路) が、L ≠ 0 の 6 極の B240 (放射・スピンの歳差) を渡さなかった

_B270_SRC = '''MOMENTUM = 1GEV;
DRIFT D1 = (L=0.6) D2 = (L=0.3);
QUAD QF = (L=0.2 K1=1.2 DY=2e-3) QD = (L=0.2 K1=-1.0 DY=-2e-3);
SEXT SF = (L=0.1 K2=30);
BEND B = (L=0.8 ANGLE=0.7853981633974483 E1=0.39269908169872414 E2=0.39269908169872414);
CAVI RFC = (VOLT=2e5 HARM=16);
LINE CELL = (QF D2 SF D2 B D2 QD D1);
LINE RING = (4*CELL RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
em = Emittance[];
Print["SPT ", SpinTune/.em];
Exit;
'''


def test_b270_script_emittance_spin_tune_includes_sextupole_body_v2(v2_mode):
    """縦ずれの四極 (閉軌道 y 2.5 mm) と L ≠ 0 の 6 極のリング: 台本の `SpinTune/.em` が `emittance6d` の nup と一致 (修正前は 3.2e-4 違った)。
    ASTRID2 (y 12 mm) で 3.6e-3、PF で 1e-3 違った原因。k64 は 6 極の放射なしのまま (gs と同じ)"""
    from pysad import emit6d
    o = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(_B270_SRC, out=o)
    spt = [float(ln.split()[1]) for ln in o.getvalue().splitlines() if ln.startswith('SPT')][0]
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=False, pol=True)
    assert spt == pytest.approx(e['spin']['nup'], abs=1e-9)
