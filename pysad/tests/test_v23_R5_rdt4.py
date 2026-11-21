# -*- coding: utf-8 -*-
"""v2.3.0 R5: 4 次の共鳴駆動項 h_jklm (6 極の 2 次を含む、pysad/rdt2nd.py) の検証。
追跡の 1 周の写像の 3 次のテイラー展開 (rdttrack.ring_map3) から Dragt-Finn の f4 を取り出して比べる。
リング: 小さな FODO (3 GeV、105.6 m、SF・SD の K2L = 4, −6)、ASTRID2 (0.58 GeV、45.7 m)。条件: pysad k64 の追跡、RF・放射なし、閉軌道のまわり、リング先頭が位相の基準、amp3 = 1e-4 √m。"""
import io
import math
import os
from types import SimpleNamespace

import numpy as np
import pytest

from pysad import rdt2nd as Q
from pysad import rdttrack as T
from pysad import resdrive as R
from pysad.elements import Quadrupole, Bend
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
P0 = 3e9


def fodo(k_sf=4.0, k_sd=-6.0):
    return load_sad_string(FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd))).build_lattice('RING')


def astrid():
    return load_sad_file(ASTRID).build_lattice('RING')


def tracked(lat, p=P0, amp3=1e-4):
    rows = lat.twiss_table()
    rm = T.ring_map3(lat.elements, p, rows[0], amp3=amp3)
    return T.rdt4_from_map(rm)


def maxrel(r, trk, thr=1e-3):
    mx = max(abs(trk[n]) for n in r.names)
    return max(abs(r.total[n] - trk[n]) / abs(trk[n]) for n in r.names if abs(trk[n]) > thr * mx)


def test_poisson_bracket_table_and_hermiticity():
    # {ζ, ζ̄} = 2i: {ζ^2 ζ̄, ζ} = 2i (2·0 − 1·1) ζ^2 → −2i ζ²。表に (2,1,0,0)×(1,0,0,0)… は 3 次どうしなので、4 次の h は h_kjml = conj(h_jklm)
    lat = fodo()
    r = Q.rdt4(lat)
    for q in Q.MONOS4:
        qb = (q[1], q[0], q[3], q[2])
        assert abs(r.total_arr[Q._IDX4[qb]] - np.conj(r.total_arr[Q._IDX4[q]])) < 1e-9 * (1 + abs(r.total_arr[Q._IDX4[q]]))


def test_thin_sextupoles_second_order_matches_tracking():
    """厚い 6 極を薄いレンズの列にした FODO: 6 極 2 次 + 偏向磁石の 3 次 + 交差項 + 8 極なし + 運動学 + フリンジ + 偏向磁石の 4 次の全体が追跡の f4 と一致 (K64 の厚い 6 極の積分器の誤差が無い場合)"""
    lat = R.slice_thin(fodo(), 4)
    trk, info = tracked(lat)
    r = Q.rdt4(lat, rows=lat.twiss_table())
    assert info['symp'] < 1e-6
    assert maxrel(r, trk) < 2e-6


def test_fodo_lambda_scan_vs_tracking_and_quadratic_in_lambda():
    h = {}
    for lam in (0.0, 1.0, 2.0, 0.5):
        lat = fodo(4.0 * lam, -6.0 * lam)
        trk, info = tracked(lat)
        r = Q.rdt4(lat, rows=lat.twiss_table())
        err = maxrel(r, trk)
        assert err < (1e-6 if lam == 0.0 else 3e-4), (lam, err)
        h[lam] = r.total
    # 解析は λ の厳密な 2 次式
    for n in ('h22000', 'h31000', 'h20020', 'h00310'):
        b = (h[2.0][n] - 2 * h[1.0][n] + h[0.0][n]) / 2
        a = h[1.0][n] - h[0.0][n] - b
        assert abs(h[0.0][n] + a * 0.5 + b * 0.25 - h[0.5][n]) < 1e-9 * abs(h[1.0][n])


def test_astrid2_vs_tracking():
    lat = astrid()
    trk, info = tracked(lat)
    r = Q.rdt4(lat, rows=lat.twiss_table())
    assert maxrel(r, trk) < 2e-3
    assert abs(r.total['h22000'] - trk['h22000']) / abs(trk['h22000']) < 2e-4


def test_adts_identity_map_plus_normal_form_correction():
    """写像の h22000・h11110・h00220 (6 極 2 次 + 偏向磁石 + 交差項) + 3 次の項の正規形の補正 (nf_correction) = adts() の同じ内訳 (R3 の厚いレンズの式)。"""
    for lat in (fodo(), astrid()):
        rows = lat.twiss_table()
        A = R.adts(lat, rows=rows)
        r0 = Q.rdt4(lat, rows=rows, kinematic=False, fringe=False, octupole=False, bend4=False)
        m0, c0 = r0.adts(), Q.nf_correction(r0.h3, r0.qx, r0.qy)
        for k in ('xx', 'xy', 'yy'):
            ex = A.sext[k] + A.bend[k] + A.cross[k]
            assert abs(m0[k] + c0[k] - ex) < 1e-9 * abs(ex)


def test_full_adts_from_rdt4_closer_to_direct_tracking_than_r3():
    """rdt4 の全部入り (+正規形の補正) と直接の追跡 (NAFF) との差 < adts().total との差 (FODO、保存値 追跡 408.990、-412.737、519.946 [1/(m rad)])"""
    lat = fodo()
    a4, _ = Q.adts_from_rdt4(lat)
    A = R.adts(lat)
    trk = {'xx': 408.9901, 'xy': -412.7365, 'yy': 519.9457}
    for k in trk:
        assert abs(a4[k] - trk[k]) / abs(trk[k]) < 2e-4
    assert sum(abs(a4[k] - trk[k]) for k in trk) < 0.5 * sum(abs(A.total[k] - trk[k]) for k in trk)


def _row_after(el, row):
    M = np.asarray(el.transfer_matrix(0.0))[:4, :4]
    b0, a0 = row.betx, row.alfx
    bx = M[0, 0] ** 2 * b0 - 2 * M[0, 0] * M[0, 1] * a0 + M[0, 1] ** 2 * (1 + a0 ** 2) / b0
    ax = -M[0, 0] * M[1, 0] * b0 + (M[0, 0] * M[1, 1] + M[0, 1] * M[1, 0]) * a0 - M[0, 1] * M[1, 1] * (1 + a0 ** 2) / b0
    mux = row.mux + math.atan2(M[0, 1], M[0, 0] * b0 - M[0, 1] * a0)
    by = M[2, 2] ** 2 * 1 + M[2, 3] ** 2
    ay = -M[2, 2] * M[3, 2] - M[2, 3] * M[3, 3]
    muy = math.atan2(M[2, 3], M[2, 2])
    return SimpleNamespace(betx=bx, alfx=ax, mux=mux, bety=by, alfy=ay, muy=muy, etx=0.0, etpx=0.0, s=0.0)


@pytest.mark.parametrize('disfrin', [1.0, 0.0])
def test_quadrupole_kinematic_trapezoid_and_fringe_vs_one_element_map(disfrin):
    """四極 1 個の写像の f4 = 運動学的な項 (K64 の積分器の台形則) + 端のフリンジ (Lee-Whiting、DISFRIN=0 のとき)。K64 の四極の本体は分割数 1 + int(L·K1L/0.2) の (線形写像 + 補正のドリフト)"""
    q = Quadrupole(name='q', L=0.3, K1=2.0)
    q.DISFRIN = disfrin
    row0 = SimpleNamespace(betx=1.0, alfx=0.0, mux=0.0, bety=1.0, alfy=0.0, muy=0.0, etx=0.0, etpx=0.0, s=0.0)
    rows = [row0, _row_after(q, row0)]
    rm = T.ring_map3([q], P0, row0, amp3=2e-4, orbit=np.zeros(4))
    trk, _ = T.rdt4_from_map(rm)
    tot = Q.kinematic_h4([q], rows) + Q.fringe_h4([q], rows)
    for nm in ('h22000', 'h11110', 'h00220', 'h40000', 'h20020', 'h31000'):
        i = Q._IDX4[tuple(int(c) for c in nm[1:5])]
        assert abs(tot[i] - trk[nm]) < 1e-7 * max(abs(trk[nm]), 1e-2), nm
    if disfrin == 0.0:
        assert abs(Q.fringe_h4([q], rows)[Q._IDX4[(1, 1, 1, 1)]]) > 0.05          # フリンジは効いている


def test_quadrupole_kinematic_exact_vs_k64_rule_differ_by_integrator_error():
    q = Quadrupole(name='q', L=0.6, K1=2.0)           # ndiv = 1 + int(0.6·1.2/0.2) = 4
    assert Q._quad_ndiv(q) == 4
    row0 = SimpleNamespace(betx=1.0, alfx=0.0, mux=0.0, bety=1.0, alfy=0.0, muy=0.0)
    rows = [row0, row0]
    a = Q.kinematic_h4([q], rows, quad_rule='k64')
    b = Q.kinematic_h4([q], rows, quad_rule='exact')
    i = Q._IDX4[(2, 2, 0, 0)]
    assert 1e-4 < abs(a[i] - b[i]) / abs(b[i]) < 0.1


def test_sector_bend_fourth_order_remainder_is_kinematic_and_edges():
    """扇型偏向磁石 (E=0、K1=0) の要素の写像の f4 から、内部の 3 次どうしの対 (薄い kick の列) を引いた残りは、運動学的な項 p⁴/8 の厳密な積分 (誤差 < 1%)"""
    lat = fodo()
    rows = lat.twiss_table()
    i = next(k for k, e in enumerate(lat.elements) if type(e).__name__ == 'Bend')
    el = lat.elements[i]
    from pysad import bendnl
    if float(el.E1) != 0.0 or float(el.E2) != 0.0:
        import dataclasses
        el = dataclasses.replace(el, E1=0.0, E2=0.0)
    elements = [el]
    r0 = rows[i]
    row1 = _row_after(el, r0)
    bu = Q.bend_units(elements, [r0, row1], 16, P0)
    h, done = Q.bend_h4_extracted(elements, [r0, row1], bu, P0)
    kin = Q.kinematic_h4(elements, [r0, row1], quad_rule='exact')
    for q in ((2, 2, 0, 0), (1, 1, 1, 1), (0, 0, 2, 2), (4, 0, 0, 0)):
        j = Q._IDX4[q]
        assert abs(h[j] - kin[j]) < 0.02 * abs(kin[j]) + 1e-3 * abs(kin[Q._IDX4[(2, 2, 0, 0)]]), q


def test_resonance_terms_accepts_fourth_order_names_and_order2_switch():
    lat = fodo()
    r = R.resonance_terms(lat, ['h21000', 'h22000', 'h31000', 'h11110'])
    r2 = Q.rdt4(lat)
    assert abs(r.total['h31000'] - r2.total['h31000']) < 1e-12
    assert abs(r.total['h21000'] - R.resonance_terms(lat, ['h21000']).total['h21000']) < 1e-12
    r1 = R.resonance_terms(lat, ['h22000'], order2=False)
    assert abs(r1.total['h22000'] - r.total['h22000']) > 100.0          # 2 次を入れないと h22000 は 313 → 59
    r_no = Q.rdt4(lat, order2=False)
    assert abs(r1.total['h22000'] - r_no.total['h22000']) < 1e-12


def test_ffs_rdt_fourth_order_names():
    from pysad import ffscmd
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)') + '''
FFS;
USE RING; CELL; CALC;
Print["R5 ", ResonanceDrivingTerms[Terms->{"h22000","h31000"}]];
Exit;
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    txt = buf.getvalue()
    assert 'H22000' in txt.upper() and 'H31000' in txt.upper()


def test_numerical_normal_form_adts_matches_direct_tracking_fodo():
    """数値の正規形 (nform.adts_from_map): 追跡の写像の h3・h4 + 正規形の補正 = 直接の追跡 (NAFF、保存値 408.9901, -412.7365, 519.9457) と 1e-5 以内"""
    from pysad import nform
    lat = fodo()
    rm = T.ring_map3(lat.elements, P0, lat.twiss_table()[0], amp3=1e-4)
    a, info = nform.adts_from_map(rm)
    for k, v in (('xx', 408.9901), ('xy', -412.7365), ('yy', 519.9457)):
        assert abs(a[k] - v) / abs(v) < 1e-5


def test_fmapres_second_order_harmonics():
    from pysad import fmapres as FR
    lat = astrid()
    rows = lat.twiss_table()
    h2, r4 = FR.second_order_terms(lat, rows)
    # 4 次の調和: 2νx−2νy (a = (2, −2))、4νx・4νy・2νx・2νy が駆動される。νx−νy (1, −1) は上下対称なので駆動されない
    assert (2, -2) in h2 and (4, 0) in h2 and (0, 4) in h2
    assert (1, -1) not in h2
    pr1 = FR.ResonancePredictor.from_lattice(lat, 2, rows)
    pr2 = FR.ResonancePredictor.from_lattice(lat, 2, rows, order2=True)
    assert (2, -2) not in pr1.harms and (2, -2) in pr2.harms
    assert set(pr1.harms) <= set(pr2.harms)
