# -*- coding: utf-8 -*-
"""v2.4.0 F2: v2 の既定の見直し (K64 の誤りと言えるもの)。k64 は変えない。
(a) B200: Touschek の損失率の δ_old = 0 の基準 (K64・Piwinski・elegant) はバンチ内の δ の分布を無視する近似。v2 の既定は BunchDelta・BunchQs つき (->False で切る、k64 は False)。
(b) B201: DynamicApertureSurvey の ExpandElementValues の既定 True は補正磁石の成分ごとの値 (K0・位置ずれ) を代表の値に揃えて消す。v2 の既定は False (k64 は True)。
(c) pysad 自身: FrequencyMap・MomentumAcceptance・NormalFormTune の TrackSetup への受け渡しが NOCALEXP でも無条件に揃えていた (補正後の状態を渡す経路)。
(d) WSPACSC の既定は 0 のまま (理由は docs/28)。"""
import contextlib
import io
import os

import pytest

from pysad import compat
from pysad.ffscmd import run_sad_string
from pysad.sadscript import Rule, S

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


def _run(body, cc='cell calc noexpand'):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        run_sad_string(FODO + '\nFFS;\nUSE RING NOEXPAND;\ncc = "%s";\nFFS[cc];\n' % cc + body + '\nstop;\n')
    return buf.getvalue().split('\n')


def _after(lines):
    return float([l for l in lines if l.startswith('AFTER')][0].split()[1])


DAS = 'LINE["DX","QF.3"] = 1.0e-4;\nFFS[cc];\nres = DynamicApertureSurvey[{{1,10},{1,1},{0}}, 50, Output->0%s];\nPrint["AFTER ", LINE["DX","QF.3"]];'


def test_das_expand_element_values_default_v2_false_k64_true(v2_mode):
    """1 台だけ変えた QF.3 の DX: v2 の既定 (False) は 1e-4 のまま、->True で代表 (QF.1) の 0 に揃う。k64 の既定は True (gs と同じ)"""
    assert abs(_after(_run(DAS % '')) - 1e-4) < 1e-12
    assert _after(_run(DAS % ',ExpandElementValues->True')) == 0.0
    assert abs(_after(_run(DAS % ',ExpandElementValues->False')) - 1e-4) < 1e-12
    compat.set_mode('k64')
    assert _after(_run(DAS % '')) == 0.0
    assert abs(_after(_run(DAS % ',ExpandElementValues->False')) - 1e-4) < 1e-12


def test_das_v2_default_keeps_corrected_state_dynamic_aperture(v2_mode):
    """ずれのある状態 (1 台の QF の DX = 5 mm) の DA: v2 の既定 (False) は ->False と同じ (ずれた状態の値。実測 15)、->True はずれが消えた状態の値 (51 = ずれなしと同じ)"""
    def da(opt):
        body = ('LINE["DX","QF.3"] = 5.0e-3;' + chr(10) + 'FFS[cc];' + chr(10) + 'res = DynamicApertureSurvey[{{1,10},{1,1},{0}}, 300, Output->0, DAPWIDTH->51' + opt + '];' + chr(10) + 'Print["SCORE ", res[[1]]];')
        lines = _run(body)
        return float([l for l in lines if l.startswith('SCORE')][0].split()[1])
    d_def, d_false, d_true = da(''), da(',ExpandElementValues->False'), da(',ExpandElementValues->True')
    assert d_def == d_false
    assert d_false < 30 and d_true == 51


FM = 'LINE["DX","QF.3"] = 1.0e-4;\nFFS[cc];\nr = FrequencyMap[1, 32, X->{1e-4, 2e-4, 2}, Y->1e-4, DP->0];\nPrint["AFTER ", LINE["DX","QF.3"]];'


def test_frequency_map_keeps_component_values_with_nocalexp():
    """補正後の状態を TrackSetup に渡す経路 (FrequencyMap・MomentumAcceptance・NormalFormTune): NOCALEXP の台本は成分ごとの値を保つ (以前は無条件に揃えて 0 に戻った)。
    既定の CALEXP の台本は従来どおり代表の値に揃える (CALC と同じ規則)"""
    assert abs(_after(_run(FM)) - 1e-4) < 1e-12
    lines = _run(FM, cc='cell calc')
    assert _after(lines) == 0.0


def test_touschek_bunch_default_v2_on_k64_off(ring_touschek):
    """TouschekLifetime[INF, INF, nz]: v2 の既定はバンチ内の δ の分布つき (BunchDelta・BunchQs)、->False で δ_c = 0 の Piwinski 型 (T2)、k64 は既定で従来のまま (Formula->"Piwinski" を明示しても)"""
    from pysad.touschek import touschek_piwinski
    from pysad.touschek_bunch import touschek_piwinski_bunch
    ffs, tt, nz, a = ring_touschek
    zero = Rule(S('BunchDelta'), 0.0)
    nq = Rule(S('BunchQs'), 0.0)
    piw = touschek_piwinski(tt, a=a)['lifetime']
    bun = touschek_piwinski_bunch(tt, a=a)['lifetime']
    old = compat._mode
    try:
        compat.set_mode('v2')
        assert abs(ffs.fn_touschek([None, None, nz, zero]) / piw - 1) < 1e-12
        assert abs(ffs.fn_touschek([None, None, nz, nq]) / bun - 1) < 1e-12          # 既定 BunchDelta + BunchQs->False = δ_c だけ
        assert bun < piw and 0.005 < piw / bun - 1 < 0.06
        full = ffs.fn_touschek([None, None, nz])                                      # 既定 = BunchDelta + BunchQs (q_s の効果も入る)
        assert full < bun and abs(full / bun - 1) < 0.05
        compat.set_mode('k64')
        tab = ffs.fn_touschek([None, None, nz])
        assert tab > 0
        assert abs(ffs.fn_touschek([None, None, nz, Rule(S('Formula'), 'Piwinski')]) / piw - 1) < 1e-12      # k64 の Piwinski 型は従来 (δ_c = 0)
    finally:
        compat.set_mode(old)


@pytest.fixture(scope='module')
def ring_touschek():
    from pysad.trksetup import TrackSetup
    old = compat._mode
    compat.set_mode('v2')
    try:
        src = open(ASTRID, encoding='utf-8').read()
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)
        g = ffs.I.globals
        ffs.flags['INTRA'] = True
        g['PBUNCH'] = 1e9
        ffs.flags['RFSW'] = True
        ffs.flags['RAD'] = True
        ffs.run_emit(quiet=True)
        sige = float(g['SIGE'])
        a = 0.006
        yield ffs, ffs.touschek, a / sige, a
    finally:
        compat.set_mode(old)


def test_local_lma_default_follows_compat(ring_touschek):
    """lma.touschek_lifetime_local_piwinski の既定 (bunch_delta=None): v2 = True・k64 = False。明示すればその値"""
    import numpy as np
    from pysad import lma
    from pysad.touschek import touschek_piwinski
    ffs, tt, nz, a = ring_touschek
    st = None
    from pysad.trksetup import TrackSetup
    s = np.array([0.0, 10.0, 20.0])
    uni = lma.LMAResult([1, 2, 3], s, ['a', 'b', 'c'], [a] * 3, [a] * 3, [0.0] * 3, [0.0] * 3, [''] * 3, [''] * 3, {}, setup=None)
    # setup は s_el を渡すので不要
    nl = int(tt['rest_info']['nlat'])
    s_el = np.linspace(0.0, 20.0, nl)
    old = compat._mode
    try:
        compat.set_mode('k64')
        r_k = lma.touschek_lifetime_local_piwinski(tt, uni, s_el=s_el)
        assert abs(r_k['rate'] / touschek_piwinski(tt, a=a)['rate'] - 1) < 1e-12
        compat.set_mode('v2')
        r_n = lma.touschek_lifetime_local_piwinski(tt, uni, s_el=s_el, bunch_delta=False)
        assert abs(r_n['rate'] / r_k['rate'] - 1) < 1e-12
        r_b = lma.touschek_lifetime_local_piwinski(tt, uni, s_el=s_el, bunch_qs=False)
        assert r_b['rate'] > r_n['rate']
    finally:
        compat.set_mode(old)
