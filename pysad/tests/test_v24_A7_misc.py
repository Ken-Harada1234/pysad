# -*- coding: utf-8 -*-
"""v2.4.0 A7〜A10: (1) Touschek の q_s の効果 (3 次元の積分 touschek_bunch.qs_ratio) が対のモンテカルロ (touschek_bunch_mc_astrid2.out.txt) と一致、
(2) fmapcheck.errors_summary が RFSW・RAD ON でも誤判定しない (stable=False の誤判定を直した)、(3) 台本の関数 NormalFormTune[] (pysad のみ)。"""
import contextlib
import io
import os
import sys

import pytest

from pysad import compat
from pysad.ffscmd import run_sad_string
from pysad.trksetup import TrackSetup

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


@pytest.fixture(scope='module')
def ring():
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
        yield st, ffs
    finally:
        compat.set_mode(old)


def test_qs_ratio_matches_monte_carlo(ring):
    """3 次元の積分の (絶対値 / kick の損失率 − 1) (δ_c と q_s) が、対のモンテカルロ (ASTRID2 スライス 76、統計 ±0.3%) と a = 0.004・0.006 で +4.84%・+2.22% に 0.6% 以内。
    δ_c だけの 1 次元の値 (+2.6%・+1.3%) とは違う (q_s の効果が +2.0%・+0.9%)"""
    from pysad.touschek import piwinski_slice_rate
    from pysad.touschek_bunch import piwinski_slice_rate_bunch, qs_ratio
    st, ffs = ring
    tt = ffs.touschek
    sl = tt['rest'][76]
    info = tt['rest_info']
    for a, mc in ((0.004, 0.0484), (0.006, 0.0222)):
        r, e = qs_ratio(sl, a, a, nu=1024)
        assert abs((r - 1) - mc) < 0.006
        dc = piwinski_slice_rate_bunch(sl, a, a, info['circ'], 1e9) / piwinski_slice_rate(sl, a, a, info['circ'], 1e9) - 1
        assert (r - 1) - dc > 0.005


def test_bunch_rate_with_qs_option(ring):
    from pysad.touschek_bunch import touschek_piwinski_bunch
    st, ffs = ring
    tt = ffs.touschek
    r1 = touschek_piwinski_bunch(tt, a=0.006)
    r2 = touschek_piwinski_bunch(tt, a=0.006, qs=True, qs_nsample=3, qs_nu=512)
    assert r2['rate'] > r1['rate'] > r2['rate0'] and r2['qs'] is not None
    assert 0.01 < r2['rate'] / r2['rate0'] - 1 < 0.05


def test_errors_summary_is_stable_with_rad_and_rfsw(ring):
    """RAD・RFSW ON の追跡の 1 周の行列は固有値の絶対値が 1 より少し小さい (減衰)。これを不安定としていた。誤差 level 1 (ASTRID2) で stable = True、結合の警告が出る (level 5)"""
    from pysad import fmapcheck as FC
    from pysad import fmapres as FR
    st0, ffs = ring
    sig = dict(quad_dxy=250e-6, quad_rot=1e-3, sext_dxy=250e-6, bend_dxy=250e-6, bend_rot=0.5e-3)
    st, _ = FR.apply_alignment_errors(st0, seed=1, **sig)
    st.cod1 = st.refine_orbit(st.cod1)
    es = FC.errors_summary(st)
    assert es['stable'] and es['closed_orbit']
    codes = [c for c, _ in FC.check_setup(st)]
    assert 'unstable' not in codes and 'coupling' in codes


def test_normal_form_tune_script_function():
    """NormalFormTune[] (pysad のみ): FODO の J¹ が ampnf.amplitude_tune と同じ、チューン・ξ が返る"""
    from pysad import ampnf as A
    from pysad.parser import load_sad_string
    FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
        'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
    body = FODO + '\nFFS;\nUSE RING;\nr = NormalFormTune[Order->5, BendCubic->False, Kinematic->False, Fringe->False, Slices->8];\n' \
        'Print["ADTS ", ADTS /. r];\nPrint["TUNES ", Tunes /. r];\nPrint["CHR ", (Chromaticity /. r)[1]];\nstop;\n'
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        run_sad_string(body)
    out = buf.getvalue()
    lat = load_sad_string(FODO).build_lattice('RING')
    r = A.amplitude_tune(lat, D=5, nslice=8)
    line = [l for l in out.split('\n') if l.startswith('ADTS')][0]
    v = [float(x) for x in line[line.index('{') + 1:line.index('}')].split(',')]
    assert abs(v[0] - r.adts['xx']) < 1e-6 * abs(v[0]) and abs(v[1] - r.adts['xy']) < 1e-6 * abs(v[1]) and abs(v[2] - r.adts['yy']) < 1e-6 * abs(v[2])
    t = [l for l in out.split('\n') if l.startswith('TUNES')][0]
    assert abs(float(t[t.index('{') + 1:t.index(',')]) - r.nu0[0]) < 1e-9
    c = [l for l in out.split('\n') if l.startswith('CHR')][0]
    assert abs(float(c[c.index('{') + 1:c.index(',')]) - r.chrom(0, 1)) < 1e-9
