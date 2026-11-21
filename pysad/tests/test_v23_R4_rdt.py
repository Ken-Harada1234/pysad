# -*- coding: utf-8 -*-
"""v2.3.0 R4: 共鳴駆動項 (RDT) h_jklmp の解析式 (resdrive.resonance_terms) を、追跡から取り出した値 (pysad/rdttrack.py: 1 周の写像の 2 次のテイラー展開、
閉軌道のヤコビ行列の δ 微分、周期点) と比べる。偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2、δ の項は分散による軌道のずれ) の寄与
(resonance_terms(bend_cubic=True)、既定) の検証。リング: 小さな FODO (3 GeV、105.6 m、SF・SD の K2L = 4, −6)、ASTRID2 (0.58 GeV、45.7 m、30° の複合機能偏向磁石)。
条件: pysad k64 の追跡 (conftest の既定)、RF・放射なし、閉軌道のまわり、リング先頭が位相の基準。"""
import copy
import dataclasses
import io
import math
import os

import numpy as np
import pytest

from pysad import bendnl as B
from pysad import rdttrack as T
from pysad import resdrive as R
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
P0 = 3e9
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
CUB5 = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
DEL5 = ('h11001', 'h00111', 'h20001', 'h00201', 'h10002')


def fodo(k_sf=4.0, k_sd=-6.0):
    return load_sad_string(FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))).build_lattice('RING')


def astrid(lam=1.0):
    lat = load_sad_file(ASTRID).build_lattice('RING')
    for e in lat.elements:
        if type(e).__name__ == 'Sextupole':
            e.K2 = e.K2 * lam
    return lat


def relc(a, b):
    return abs(a - b) / abs(b)


def tracked(lat, p, **kw):
    return T.tracked_rdt_table(lat.elements, p, lat.twiss_table()[0], **kw)


@pytest.fixture(scope='module')
def fodo_pair():
    l1, l0 = fodo(), fodo(0.0, 0.0)
    return l1, tracked(l1, P0), tracked(l0, P0)


@pytest.fixture(scope='module')
def astrid_pair():
    l1, l0 = astrid(1.0), astrid(0.0)
    return l1, tracked(l1, 0.58e9), tracked(l0, 0.58e9)


def test_fodo_sextupole_part_matches_one_turn_map_to_1e_6(fodo_pair):
    """6 極の強さ ×1 − ×0 の 1 周の写像の 2 次のテイラー展開 (追跡) から取り出した h は、解析 (6 極のみ) と 1e-6 以内で一致 (2 次の部分は kick の和なので 1 次の摂動で厳密)。
    位相・符号・厚いレンズの積分・規格化 2^{-(j+k+l+m)} を含む"""
    lat, t1, t0 = fodo_pair
    r = R.resonance_terms(lat, CUB5, bend_cubic=False)
    for nm in CUB5:
        assert relc(r.total[nm], t1['cubic'][nm] - t0['cubic'][nm]) < 1e-6, nm


def test_fodo_total_with_bend_matches_one_turn_map_and_old_does_not(fodo_pair):
    """偏向磁石の項を入れた解析 (既定) は ×1 の追跡と 1e-6 以内、6 極のみは 1.4〜3.4% 違う (R3 までの値)"""
    lat, t1, t0 = fodo_pair
    rb = R.resonance_terms(lat, CUB5)
    ra = R.resonance_terms(lat, CUB5, bend_cubic=False)
    assert rb.bend_included and not ra.bend_included
    errs_old = []
    for nm in CUB5:
        assert relc(rb.total[nm], t1['cubic'][nm]) < 1e-6, nm
        assert abs(rb.total_bend[nm] - t0['cubic'][nm]) < 1e-6 * abs(rb.total[nm]), nm        # 偏向磁石の分 (追跡の ×0)
        errs_old.append(relc(ra.total[nm], t1['cubic'][nm]))
        assert abs(rb.total[nm] - (rb.total_sext[nm] + rb.total_bend[nm])) < 1e-12
    assert 0.010 < min(errs_old) and max(errs_old) < 0.04


def test_other_cubic_terms_are_zero(fodo_pair):
    lat, t1, t0 = fodo_pair
    for q in T.CUBIC_MONOS:
        nm = T.mono_name(q)
        if nm in CUB5 or nm in [x[::-1] for x in CUB5] or abs(t1['cubic'][nm]) > 1e-6:
            continue
        assert abs(t1['cubic'][nm]) < 1e-8


def test_astrid2_bend_terms_matter_and_match(astrid_pair):
    """ASTRID2 (30° の複合機能偏向磁石): 偏向磁石の項は h の 0.8〜11% (6 極のみは 0.8〜11% の差)、入れると追跡 (×1) と 1e-5 以内。6 極のみの部分は 1e-6"""
    lat, t1, t0 = astrid_pair
    rb = R.resonance_terms(lat, CUB5)
    ra = R.resonance_terms(lat, CUB5, bend_cubic=False)
    old = []
    for nm in CUB5:
        assert relc(rb.total[nm], t1['cubic'][nm]) < 1e-5, nm
        assert relc(ra.total[nm], t1['cubic'][nm] - t0['cubic'][nm]) < 1e-6, nm
        old.append(relc(ra.total[nm], t1['cubic'][nm]))
    assert max(old) > 0.08 and min(old) > 0.005


def test_delta_terms_sextupole_part_matches_tracking(fodo_pair, astrid_pair):
    """h11001・h00111・h20001・h00201 (閉軌道のヤコビ行列の δ 微分、強さ ×1 − ×0) と h10002 (2 次の閉軌道、×1 − ×0) が解析 (6 極のみ) と 1e-5 以内"""
    for lat, t1, t0, tol in (fodo_pair + (1e-5,), astrid_pair + (1e-5,)):
        r = R.resonance_terms(lat, DEL5, bend_cubic=False)
        for nm in DEL5:
            if nm == 'h10002':
                d = t1['h10002'] - t0['h10002']
            else:
                d = t1['delta1'][nm] - t0['delta1'][nm]
            assert relc(r.total[nm], d) < tol, (nm, r.total[nm], d)


def test_thin_sextupole_closed_forms_of_delta_terms():
    """薄い 6 極 (FODO の SF・SD を slice_thin(4) した薄レンズ列) の δ の項の閉じた式 (H = (K2L/6)(x³ − 3xy²)、x → x + ηδ から導出):
    h11001 = ¼ΣK2L η βx、h00111 = −¼ΣK2L η βy、h20001 = ⅛ΣK2L η βx e^{2iμx}、h00201 = −⅛ΣK2L η βy e^{2iμy}、h10002 = ¼ΣK2L η² √βx e^{iμx}"""
    lat = R.slice_thin(fodo(), 4)
    rows = lat.twiss_table()
    s = dict.fromkeys(DEL5, 0j)
    for i, e in enumerate(lat.elements):
        if type(e).__name__ != 'Multipole' or not (e.KnL or {}).get(2):
            continue
        k, r = float(e.KnL[2]), rows[i]
        bx, by, eta, mx, my = float(r.betx), float(r.bety), float(r.etx), float(r.mux), float(r.muy)
        s['h11001'] += k * eta * bx / 4
        s['h00111'] += -k * eta * by / 4
        s['h20001'] += k * eta * bx * np.exp(2j * mx) / 8
        s['h00201'] += -k * eta * by * np.exp(2j * my) / 8
        s['h10002'] += k * eta * eta * math.sqrt(bx) * np.exp(1j * mx) / 4
    for method in ('analytic', 'numeric'):
        r = R.resonance_terms(lat, DEL5, rows=rows, method=method, nslice=4, bend_cubic=False)
        for nm in DEL5:
            assert abs(r.total[nm] - s[nm]) < 1e-12 * max(1.0, abs(s[nm])), (method, nm)


def test_sector_bend_atoms_equal_independent_integral_of_curvature_and_hard_edge():
    """E1=E2=0・K1=0 の偏向磁石 (FODO の B): 偏向磁石の分 (追跡の写像から取り出した原子の和) = 曲率の項 (h/2)x(px²+py²) の ds の積分 + 硬い端 ±(h/2)px y²
    (x → x+ηδ、px → px+η'δ の展開も) の、扇型の偏向磁石の解析的な w・η・η' (cos/sin) の中点則 (400 分割) の積分と、p=0 は 2e-3、p≥1 は 5e-3 以内
    (差は 16 分割の O(ds²))"""
    lat = fodo()
    for e in lat.elements:
        if type(e).__name__ == 'Bend':
            e.E1 = 0.0
            e.E2 = 0.0
    rows = lat.twiss_table()
    terms = [(2, 1, 0, 0, 0), (3, 0, 0, 0, 0), (1, 0, 1, 1, 0), (1, 0, 2, 0, 0), (1, 0, 0, 2, 0), (1, 1, 0, 0, 1), (0, 0, 1, 1, 1), (2, 0, 0, 0, 1), (0, 0, 2, 0, 1), (1, 0, 0, 0, 2)]
    rb = R.resonance_terms(lat, terms, bend_cubic=True)
    ref = np.zeros(len(terms), dtype=complex)
    nsl = 400
    for i, e in enumerate(lat.elements):
        if type(e).__name__ != 'Bend':
            continue
        h = e.ANGLE / e.L
        rho = 1.0 / h
        r0 = rows[i]
        w0 = np.array([math.sqrt(r0.betx) * np.exp(1j * r0.mux), (1j - r0.alfx) / math.sqrt(r0.betx) * np.exp(1j * r0.mux),
                       math.sqrt(r0.bety) * np.exp(1j * r0.muy), (1j - r0.alfy) / math.sqrt(r0.bety) * np.exp(1j * r0.muy)])
        e0 = np.array([r0.etx, r0.etpx, 0.0, 0.0])

        def prop(s):
            c, sn = math.cos(s * h), math.sin(s * h)
            Mx = np.array([[c, rho * sn], [-sn / rho, c]])
            w = np.concatenate([Mx @ w0[:2], np.array([w0[2] + s * w0[3], w0[3]])])
            ee = np.concatenate([Mx @ e0[:2] + np.array([rho * (1 - c), sn]), [0.0, 0.0]])
            return w, ee
        ds = e.L / nsl
        f_b = {(1, 2, 0, 0): 0.5 * h * ds, (1, 0, 0, 2): 0.5 * h * ds}      # (x, px, y, py) の指数: x px² と x py²
        f_e1, f_e2 = {(0, 1, 2, 0): 0.5 * h}, {(0, 1, 2, 0): -0.5 * h}
        for k in range(nsl):
            w, ee = prop((k + 0.5) * ds)
            ref += B.harmonic_amplitudes_shift(f_b, w, ee, terms)[:, 0]
        w, ee = prop(0.0)
        ref += B.harmonic_amplitudes_shift(f_e1, w, ee, terms)[:, 0]
        w, ee = prop(e.L)
        ref += B.harmonic_amplitudes_shift(f_e2, w, ee, terms)[:, 0]
    for q, v in zip(terms, ref):
        nm = T.mono_name(q)
        got = rb.total_bend[nm]
        tol = 2e-3 if q[4] == 0 else 5e-3
        assert abs(got - v) < tol * max(abs(v), 1e-3 * abs(rb.total[nm])), (nm, got, v)


def test_natural_plus_bend_chromaticity_matches_calc_without_sextupoles():
    """6 極なし: CALC の ξ = 自然クロマティシティの式 (−(1/4π)∮(K1+h²)β + 端) + 偏向磁石の 3 次の項の分散による軌道のずれ (h11001/π の偏向磁石の分)。
    FODO は 3e-3 (E1=E2≠0 の端の 3 次以上)、ASTRID2 は 1e-3 以内 (ξ = −6.2, −11.8)。偏向磁石の分を入れないと FODO で 0.34・0.11、ASTRID2 で 0.08・3.5 ずれる"""
    for lat, tol in ((fodo(0.0, 0.0), 4e-3), (astrid(0.0), 1.2e-3)):
        c = R.chromaticity_numeric(lat, order=2)['d1']
        n = R.chromaticity_natural(lat)
        b = R.chromaticity_bend(lat)
        for i in (0, 1):
            assert abs(c[i] - n[i] - b[i]) < tol, (i, c[i], n[i], b[i])
            assert abs(b[i]) > 0.05


def test_ring_map_independent_of_amplitude_and_phase_reference():
    lat = fodo()
    rows = lat.twiss_table()
    p = P0
    h1 = T.rdt_from_map(T.ring_map(lat.elements, p, rows[0], amp=1e-5))
    h2 = T.rdt_from_map(T.ring_map(lat.elements, p, rows[0], amp=2e-5))
    for nm in CUB5:
        assert relc(h2[nm], h1[nm]) < 1e-5
    i0 = 100
    els = lat.elements[i0:] + lat.elements[:i0]
    hr = T.rdt_from_map(T.ring_map(els, p, rows[i0], amp=1e-5))
    ra = R.resonance_terms(lat, CUB5, ref=i0, rows=rows)
    qx, qy = float(rows[-1].mux) / (2 * math.pi), float(rows[-1].muy) / (2 * math.pi)
    # ref=i は位相の基準を変えるだけ (全要素に同じ位相)。ref から始めた 1 周の写像の h は、ref より前の要素の位相に 1 周分の e^{2πi((j-k)νx+(l-m)νy)} がつく
    before = np.array([idx < i0 for idx in ra.elem_idx])
    for nm, q in zip(ra.names, ra.terms):
        if nm not in CUB5:
            continue
        k = ra.names.index(nm)
        wrap = np.exp(2j * math.pi * ((q[0] - q[1]) * qx + (q[2] - q[3]) * qy))
        h_map = ra.per_element[~before, k].sum() + wrap * ra.per_element[before, k].sum()
        assert relc(h_map, hr[nm]) < 1e-6, nm


def test_resonance_terms_bend_option_delta_warning_and_per_element():
    lat = fodo()
    rb = R.resonance_terms(lat, CUB5)
    assert abs(rb.per_element.sum(axis=0) - np.array([rb.total[n] for n in rb.names])).max() < 1e-9
    assert rb.elem_idx == sorted(rb.elem_idx) and 'B' in rb.elem_names
    rd = R.resonance_terms(lat, CUB5, delta=1e-3)
    assert rd.bend_included          # v2.3.0 R6: δ≠0 でも偏向磁石 (δ の閉軌道のまわり) の p=0 の項を入れる
    ra = R.resonance_terms(lat, CUB5, bend_cubic=False)
    assert max(abs(rb.total_sext[n] - ra.total[n]) for n in rb.names) < 1e-8          # twiss_table の再計算の丸め
    assert all(abs(v) == 0 for v in ra.total_bend.values())


def test_fixed_points_near_3nux_resonance_give_h30000_fodo():
    """FODO の νx を 3νx = 20 − 3e-3 の近くに動かし (QF の K1 を調整)、3 周写像の不動点 (y=0) の J(φ) の 3φ 成分から |h30000| を求める。
    24 セルの対称リングで共鳴は駆動されず (h30000 ∝ ε)、不動点は曲線になる。正規形では J(φ) の 3φ 成分の振幅は 2|h|u0/(π a)、位相は arg h30000 (+3π ε の程度)。
    解析 (6 極 + 偏向磁石) と 0.5% 以内 (6 極のみは 1.6% 違う)"""
    from scipy.optimize import brentq
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')

    def build(k1):
        return load_sad_string(src.replace('QF = (L=0.3 K1=0.7)', 'QF = (L=0.3 K1=%r)' % k1)).build_lattice('RING')
    target = 20.0 / 3.0 - 1e-3
    k1 = brentq(lambda k: float(build(k).twiss_table()[-1].mux) / (2 * math.pi) - target, 0.68, 0.72, xtol=1e-13)
    lat = build(k1)
    rows = lat.twiss_table()
    eps = float(rows[-1].mux) / (2 * math.pi) - 20.0 / 3.0
    hb = R.resonance_terms(lat, ('h30000',), rows=rows).total['h30000']
    hs = R.resonance_terms(lat, ('h30000',), rows=rows, bend_cubic=False).total['h30000']
    A = R.adts(lat, rows=rows)
    a_nr = A.total['xx'] + 12.0 * abs(hb) ** 2 / (math.pi ** 2 * eps)
    pred = T.predict_periodic_J(eps, a_nr, abs(hb))
    pts = []
    for jg in pred:
        pts += T.periodic_points(lat.elements, P0, rows[0], 3, jg, nphase=16)
    assert len(pts) >= 20
    ph = np.array([p[1] for p in pts])
    J = np.array([p[0] for p in pts])
    cols = [np.ones_like(ph)] + [f(k * ph) for k in range(1, 7) for f in (np.cos, np.sin)]
    co = np.linalg.lstsq(np.array(cols).T, J, rcond=None)[0]
    A3 = math.hypot(co[5], co[6])
    h_inf = math.pi * a_nr * A3 / (2 * math.sqrt(2 * co[0]))
    assert abs(h_inf / abs(hb) - 1) < 5e-3
    assert abs(h_inf / abs(hs) - 1) > 0.008
    assert abs(math.remainder(math.atan2(-co[6], co[5]) - math.atan2(hb.imag, hb.real), 2 * math.pi)) < 0.03


def test_ffs_rdt_bend_option_and_chromaticity_parts():
    from pysad import ffscmd, compat
    old = compat._mode
    compat.set_mode('v2')
    try:
        src = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)') + '''
FFS;
USE RING; CELL; CALC;
Print["OFF ", ResonanceDrivingTerms[Terms->{"h30000","h21000"}, BendCubic->False]];
Print["ON ", ResonanceDrivingTerms[Terms->{"h30000","h21000"}]];
Print["CHR ", ChromaticityTerms[Order->2]];
Exit;
'''
        buf = io.StringIO()
        ffscmd.run_sad_string(src, out=buf)
    finally:
        compat.set_mode(old)
    txt = buf.getvalue()
    assert 'OFF {H30000->{-.37313' in txt.replace('\n', ' ').replace('  ', ' ') or 'H30000->{-.37313' in txt
    assert 'H30000->{-.38503' in txt or 'H30000->{-.38504' in txt
    assert 'NaturalPart->{-7.8537' in txt and 'BendPart->{.3446' in txt


def test_v2_mode_bend_terms_follow_the_current_tracker(v2_mode):
    """v2 モードでも、偏向磁石の項は現在のモードの追跡の写像から取り出すので、同じモードの追跡の 1 周の写像と 1e-6 以内で一致する (FODO・ASTRID2。R3 は v2 の追跡との照合が未実施だった)"""
    for lat, p in ((fodo(), P0), (astrid(1.0), 0.58e9)):
        t1 = tracked(lat, p, with_delta=False)
        r = R.resonance_terms(lat, CUB5)
        for nm in CUB5:
            assert relc(r.total[nm], t1['cubic'][nm]) < 1e-6, nm
