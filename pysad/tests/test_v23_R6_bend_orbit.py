# -*- coding: utf-8 -*-
"""v2.3.0 R6: 偏向磁石の閉軌道・δ のまわりの分解 (bendnl.bend_atoms_orbit)、K0≠0 の偏向磁石 (挿入光源 = ASTRID2 MPW)、δ≠0 の ADTS、
偏向磁石の生成関数の閉形式、adts() の既定 (四極の台形則 + 偏向磁石の 4 次の項)、'log'/'lin' (U2: BEND の K1 の定義)。
リング: 小さな FODO (3 GeV、105.6 m)、ASTRID2 (0.58 GeV、45.7 m)、ASTRID2 MPW (`ASTRID2_v01_02_id.sad` の RINGMPW)。条件: pysad の追跡 (RF・放射なし・閉軌道のまわり)。"""
import math
import os
import re

import numpy as np
import pytest

from pysad import bendnl as B
from pysad import compat
from pysad import nform
from pysad import rdt2nd as Q
from pysad import rdttrack as T
from pysad import resdrive as R
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
MPW = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
KEKLS = os.path.join(HERE, 'kek-ls_v03_12.lat')
M3 = Q.MONOS3


def fodo(e1=None, e2=None):
    s = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
    if e1 is not None:
        s = s.replace('E1=0.1308996938995747 E2=0.1308996938995747', 'E1=%r E2=%r' % (e1, e2))
    return load_sad_string(s).build_lattice('RING')


def astrid():
    return load_sad_file(ASTRID).build_lattice('RING')


def mpw(lam=1.0):
    s = open(MPW, encoding='utf-8').read()
    if lam != 1.0:
        s = re.sub(r'(\bK2\s*=\s*)([-+]?[0-9.]+(?:[eE][-+]?\d+)?)', lambda m: m.group(1) + repr(float(m.group(2)) * lam), s)
    return load_sad_string(s).build_lattice('RINGMPW')


def bends(lat, k0=False):
    return [i for i, e in enumerate(lat.elements) if type(e).__name__ == 'Bend' and e.L > 0 and ((float(e.K0) != 0.0) == k0) and (k0 or e.ANGLE != 0)]


def test_orbit_atoms_equal_legacy_atoms_on_axis():
    """閉軌道 0・δ=0・K0=0 では、閉軌道のまわりの取り出し (R6、スライスは E=0 の部分偏向磁石の写像そのもの) が R3 の分解 (A_in の解析式つき) と原子の和で 1e-12 一致する"""
    for lat in (fodo(), astrid(), fodo(0.2, -0.07)):
        rows = lat.twiss_table()
        i = bends(lat)[0]
        el = lat.elements[i]
        a0, _ = B.bend_atoms(el, rows[i], rows[i + 1], M3, nslice=8)
        a1, i1 = B.bend_atoms_orbit(el, rows[i], rows[i + 1], M3, nslice=8)
        assert np.abs(a0.sum(axis=1) - a1.sum(axis=1)).max() < 1e-12 * np.abs(a0).max()
        assert i1['w_err'] < 1e-7 and i1['z_err'] < 1e-12


def test_mpw_k0_bends_cubic_terms_match_tracked_map():
    """MPW の K0≠0 の偏向磁石 23 台 (内部の軌道が ±0.032 rad 傾く) の 3 次の項を入れると、h の全部が追跡の 1 周の写像のテイラー展開と 1e-7 以内
    (R4: h10020・h10200 で 4.6e-5・2.0e-5 ずれていた)。w と閉軌道は CALC の Twiss と 1e-12 以内で一致"""
    lat = mpw()
    rows = lat.twiss_table()
    i = bends(lat, k0=True)[1]
    _, info = B.bend_atoms_orbit(lat.elements[i], rows[i], rows[i + 1], M3, nslice=16)
    assert info['w_err'] < 1e-10 and info['z_err'] < 1e-12
    names = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
    r = R.resonance_terms(lat, names, rows=rows)
    rm = T.ring_map(lat.elements, 3e9, rows[0], amp=1e-5)
    h = T.rdt_from_map(rm)
    for n in names:
        assert abs(r.total[n] - h[n]) / abs(h[n]) < 1e-7, n
    assert not any('K0≠0' in w for w in r.warnings)


def test_mpw_fourth_order_without_sextupoles_matches_tracked_map():
    """MPW の 6 極 ×0 の 4 次の項 (偏向磁石の 4 次の項 = 軌道が傾いた K0 の磁石の運動学・端を含む) が、追跡の 3 次のテイラー展開 (Dragt-Finn の f4) と 1e-4 以内
    (R5: 8.7e-3)"""
    lat = mpw(0.0)
    rows = lat.twiss_table()
    rm = T.ring_map3(lat.elements, 0.58e9, rows[0], amp3=1e-4)
    ht, _ = T.rdt4_from_map(rm)
    r = Q.rdt4(lat, rows=rows)
    big = max(abs(v) for v in ht.values())
    for nm, v in ht.items():
        if abs(v) > 1e-3 * big:
            assert abs(r.total[nm] - v) / abs(v) < 1e-4, nm


def test_mpw_chromaticity_with_k0_bends():
    """MPW の自然クロマティシティ + 偏向磁石の分 (K0≠0 の磁石は要素の追跡から直接) が CALC の ξ (6 極なし) と 2e-3 以内 (R4: ξy が 0.17 ずれた)。
    素子ごとの δ の項 (スロープ座標) の全素子の和 (ドリフト・四極を含む) も CALC と 1e-4 (4e-7・3e-5。自然な項の分け方に依らない)"""
    lat = mpw(0.0)
    rows = lat.twiss_table()
    nat = R.chromaticity_natural(lat, rows=rows)
    bend = R.chromaticity_bend(lat, rows=rows)
    ref = R.chromaticity_numeric(lat, order=1)['d1']
    assert abs(nat[0] + bend[0] - ref[0]) < 2e-3 and abs(nat[1] + bend[1] - ref[1]) < 2e-3
    tot = np.zeros(2)
    for i, el in enumerate(lat.elements):
        if float(getattr(el, 'L', 0.0) or 0.0) > 0:
            d = B.element_delta_terms(el, rows[i], rows[i + 1], 3e9)
            tot += [d[(1, 1, 0, 0, 1)].real / math.pi, d[(0, 0, 1, 1, 1)].real / math.pi]
    assert abs(tot[0] - ref[0]) < 1e-4 and abs(tot[1] - ref[1]) < 1e-4


def test_mpw_adts_r6_is_rdt4_route_and_close_to_converged_tracking():
    """MPW λ=2 の adts(fourth='r6').total = rdt4 + 正規形の補正 (1e-10)。R3 の値 (fourth='r3') との差は 4 次の項 (四極の台形則 + 偏向磁石)。
    直接の追跡 (Jmax 2e-8、order 2) の 6 極 ×2 の yy は 4νy ≈ 9 の近くで J² の項が大きく (νy = 2.2492)、フィットが −0.5% 偏る。Jmax 5e-9 では 5e-4 以内
    (原因が物理の欠落でなくフィットの偏りであることの確認。追跡は 1 回)"""
    from pysad import adtsdirect as D
    lat = mpw(2.0)
    rows = lat.twiss_table()
    A = R.adts(lat, rows=rows, fourth='r6')
    a4, _ = Q.adts_from_rdt4(lat, rows=rows)
    for k in ('xx', 'xy', 'yy'):
        assert abs(A.total[k] - a4[k]) < 1e-9 * abs(a4[k])
    tr = D.direct_adts(lat, 0.58e9, 5e-9, order=2, nturn=1024)
    t = {'xx': tr['xx'], 'xy': tr['xy_mean'], 'yy': tr['yy']}
    for k in t:
        assert abs(A.total[k] - t[k]) / abs(t[k]) < 1e-3, k
    A3 = R.adts(lat, rows=rows, fourth='r3')
    assert abs(A3.total['yy'] - t['yy']) > abs(A.total['yy'] - t['yy'])


def test_adts_r6_option_equals_rdt4_route_and_default_is_r3():
    """adts(fourth='r6') (精密版のオプション。既定は 'r3') は adts_from_rdt4 と 1e-12 で一致 (FODO・ASTRID2)、fourth='r3' は R3 の値 (kin は Simpson・fringe は偏向磁石を含む、bend4=0)。
    追跡 (保存値) との差は R3 (既定) の 1e-4〜2.4e-3 から、'r6' で 2e-4 以内 (四極の積分器の台形則など数値積分器の端数を含む)"""
    for lat in (fodo(), astrid()):
        rows = lat.twiss_table()
        A = R.adts(lat, rows=rows, fourth='r6')
        a4, _ = Q.adts_from_rdt4(lat, rows=rows)
        A3 = R.adts(lat, rows=rows)
        for k in ('xx', 'xy', 'yy'):
            assert abs(A.total[k] - a4[k]) < 1e-12 * abs(a4[k]) + 1e-12
            assert A3.bend4[k] == 0.0 and A3.fourth == 'r3'
    lat = fodo()
    A = R.adts(lat, fourth='r6')
    trk = {'xx': 408.9901, 'xy': -412.7365, 'yy': 519.9457}
    for k in trk:
        assert abs(A.total[k] - trk[k]) / abs(trk[k]) < 2e-4


def test_bend_adts_independent_of_nslice_with_extracted_fourth_order():
    """fourth='r6' では偏向磁石の 4 次の項が要素の写像から取り出した値 − 内部の対の項 (原子の数に依る) なので、合計は nslice (4・8・64) に依らない
    (R3 の値は O(1/n²): 1e-6〜6e-5)。FODO は 1e-9 以内、ASTRID2 は 1e-6 以内"""
    for lat, tol in ((fodo(), 1e-9), (astrid(), 1e-6)):
        rows = lat.twiss_table()
        t = {n: R.adts(lat, rows=rows, bend_nslice=n, fourth='r6').total for n in (4, 8, 64)}
        for k in ('xx', 'xy', 'yy'):
            assert abs(t[4][k] - t[64][k]) < tol * abs(t[64][k]) and abs(t[8][k] - t[64][k]) < tol * abs(t[64][k])


def test_bend_edge_closed_form_and_analytic_h3():
    """偏向磁石の端の生成関数の閉形式 (E1・E2≠0。h = ANGLE/L、t = tan E: f_e1 = (h/2)sec²E px y² + h t² x y py − (h t²/2) x² px − (h² t³/2) x y² − (h² t³/6) x³ など) が、
    追跡の写像から取り出した係数 (bend_decomposition) と 1e-6 以内。本体 H = (h/2) x (px² + py²) + (K1 h/6) x³ の Gauss 積分 + 端の閉形式 = 取り出した原子の和 (K1=0 の
    FODO 3 種 (E=0.034・0.131・非対称) で 1e-7、ASTRID2 (K1≠0) で 1e-5)"""
    import dataclasses
    base = fodo().elements[[i for i, e in enumerate(fodo().elements) if type(e).__name__ == 'Bend'][0]]
    for ang, L, E1, E2 in ((0.2618, 1.2, 0.2094, -0.0785), (0.5, 0.8, 0.3, 0.1), (0.2, 1.5, -0.25, 0.4)):
        el = dataclasses.replace(base, name='x', L=L, ANGLE=ang, E1=E1, E2=E2)
        d = B.bend_decomposition(el, 4, 3e9)
        f1, f2 = B.bend_edge_generators(el)
        for ex, fa in ((d['f_e1'], f1), (d['f_e2'], f2)):
            keys = set(ex) | set(fa)
            assert max(abs(ex.get(k, 0.0) - fa.get(k, 0.0)) for k in keys) < 1e-6 * max(abs(v) for v in fa.values())
    for lat, tol in ((fodo(), 1e-7), (fodo(0.5, 0.5), 1e-7), (fodo(0.8, -0.3), 1e-7), (astrid(), 1e-5)):
        rows = lat.twiss_table()
        i = bends(lat)[0]
        el = lat.elements[i]
        a, _ = B.bend_atoms(el, rows[i], rows[i + 1], M3, nslice=16)
        an = B.bend_analytic_h3(el, rows[i], rows[i + 1], M3)
        assert np.abs(a.sum(axis=1) / 2 ** 1.5 - an).max() < tol * np.abs(an).max()


def test_delta_adts_uses_unscaled_multipole_strengths_tracked_optics_option():
    """δ≠0: 追跡の kick (6 極・四極・端部のフリンジ) は p0 で規格化した場の強さで δ に依らない (gs と pysad の追跡で確認)ので、解析も強さを 1/(1+δ) しない。
    FODO・ASTRID2 の δ = ±0.002 の adts(delta, optics='tracked', fourth='r6') (R6 の精密版のオプション。既定は CALC の光学・fourth='r3') と、追跡の写像の数値の正規形 (nform.adts_at_delta) が 1.5e-3 以内。δ の微分 (中心差分) は数値の正規形と
    FODO で 5%、ASTRID2 の xx・xy で 1% 以内 (R2 の解析は 0.3〜7 倍)"""
    lat = fodo()
    nl = R.collect_nonlinear(lat.elements, lat.twiss_table(), 0.01, (0, 0), [])
    nl0 = R.collect_nonlinear(lat.elements, lat.twiss_table(), 0.0, (0, 0), [])
    assert abs(nl[0].strengths[2]) == abs(nl0[0].strengths[2])
    for lat, p, tol in ((fodo(), 3e9, 1.5e-3), (astrid(), 0.58e9, 1.5e-3)):
        for d in (0.002, -0.002):
            a = R.adts(lat, delta=d, momentum_eV=p, optics='tracked', fourth='r6').total
            n, _ = nform.adts_at_delta(lat, p, d)
            for k in ('xx', 'xy', 'yy'):
                assert abs(a[k] - n[k]) / abs(n[k]) < tol, (d, k)
    dn = nform.adts_delta_derivative(astrid(), 0.58e9, 0.002, method='nform')
    da = nform.adts_delta_derivative(astrid(), 0.58e9, 0.002, optics='tracked', fourth='r6')
    assert abs(da['xx'] - dn['xx']) < 1e-2 * abs(dn['xx']) and abs(da['xy'] - dn['xy']) < 1e-2 * abs(dn['xy'])


def test_tracked_rows_reproduce_calc_on_axis():
    """δ=0 の追跡の光学 (tracked_rows) が CALC の Twiss と一致 (β 1e-6、位相進み 1e-7、閉軌道 1e-12)"""
    lat = fodo()
    rows = lat.twiss_table()
    tr = T.tracked_rows(lat.elements, 3e9, 0.0, rows)
    for a, b in zip(tr[::7], rows[::7]):
        assert abs(a.betx - float(b.betx)) < 1e-6 * float(b.betx) and abs(a.bety - float(b.bety)) < 1e-6 * float(b.bety)
        assert abs(a.mux - float(b.mux)) < 1e-6 and abs(a.muy - float(b.muy)) < 1e-6
        assert abs(a.x_cod) < 1e-12


def test_nonresonant_correction_to_all_fourth_order_harmonics_vs_numerical_lie_transform():
    """3 次の項の非共鳴の摂動が 4 次の全 35 調和に足す補正 nf_correction_h4 (F4' = F4 + ½{R⁻¹W, W}、W_q = h_q/(1 − e^{iφ_q})) の独立な検証:
    2 つの薄い 6 極 + 回転の 1 周写像で、W のハミルトン流 (RK4) で写像を共役変換して (3 次の項が 1e-9 以下に消える) 4 次の項を取り出し、F4 + 補正と 1e-6 以内
    (補正なしの F4 とは最大 |補正| ≈ 1 違う)。対角 (k=0) は nf_correction と 1e-12 で一致"""
    from types import SimpleNamespace
    qx, qy = 0.2137, 0.3319

    def rot(mu):
        return np.array([[math.cos(mu), math.sin(mu)], [-math.sin(mu), math.cos(mu)]])

    def Rm(ax, ay):
        R4 = np.zeros((4, 4))
        R4[:2, :2], R4[2:, 2:] = rot(2 * math.pi * ax), rot(2 * math.pi * ay)
        return R4

    def kick(Z, K):
        Z = Z.copy()
        x, y = Z[:, 0], Z[:, 2]
        Z[:, 1] -= 0.5 * K * (x * x - y * y)
        Z[:, 3] += K * x * y
        return Z
    Ra = Rm(0.37 * qx, 0.52 * qy)
    Rb = Rm(qx, qy) @ np.linalg.inv(Ra)

    def ring(Z):
        return kick(kick(Z, 3.0) @ Ra.T, -2.2) @ Rb.T
    row = SimpleNamespace(betx=1.0, alfx=0.0, mux=0.0, bety=1.0, alfy=0.0, muy=0.0)

    def hh(f):
        rm = T.ring_map3(None, 3e9, row, amp=1e-5, amp3=2e-3, orbit=np.zeros(4), tracker=f)
        return T.rdt_from_map(rm), T.rdt4_from_map(rm)[0]
    h3, h4 = hh(ring)
    h3v = np.array([h3['h%d%d%d%d0' % q] for q in M3])
    h4v = np.array([h4['h%d%d%d%d0' % q] for q in Q.MONOS4])
    c = Q.nf_correction_h4(h3v, qx, qy)
    ix = Q._IDX4
    ref = Q.nf_correction(h3v, qx, qy)
    assert abs(4 * c[ix[(2, 2, 0, 0)]].real / math.pi - ref['xx']) < 1e-9 and abs(2 * c[ix[(1, 1, 1, 1)]].real / math.pi - ref['xy']) < 1e-9
    ph = np.array([2 * math.pi * ((q[0] - q[1]) * qx + (q[2] - q[3]) * qy) for q in M3])
    W = h3v / (1 - np.exp(1j * ph))

    def flow(Z, sgn, nst=16):
        Z = Z.copy()
        dt = sgn / nst

        def grad(Zc):
            zx, zy = Zc[:, 0] - 1j * Zc[:, 1], Zc[:, 2] - 1j * Zc[:, 3]
            zxb, zyb = np.conj(zx), np.conj(zy)
            g = np.zeros((len(Zc), 4), dtype=complex)
            for i, (j, k, l, m) in enumerate(M3):
                def pw(z, n):
                    return z ** n
                dzx = j * pw(zx, j - 1) * pw(zxb, k) * pw(zy, l) * pw(zyb, m) * W[i] if j else 0
                dzxb = k * pw(zx, j) * pw(zxb, k - 1) * pw(zy, l) * pw(zyb, m) * W[i] if k else 0
                dzy = l * pw(zx, j) * pw(zxb, k) * pw(zy, l - 1) * pw(zyb, m) * W[i] if l else 0
                dzyb = m * pw(zx, j) * pw(zxb, k) * pw(zy, l) * pw(zyb, m - 1) * W[i] if m else 0
                g[:, 0] += dzx + dzxb
                g[:, 1] += -1j * dzx + 1j * dzxb
                g[:, 2] += dzy + dzyb
                g[:, 3] += -1j * dzy + 1j * dzyb
            return g.real

        def vf(Zc):
            g = grad(Zc)
            return np.stack([g[:, 1], -g[:, 0], g[:, 3], -g[:, 2]], axis=1)
        for _ in range(nst):
            k1 = vf(Z)
            k2 = vf(Z + 0.5 * dt * k1)
            k3 = vf(Z + 0.5 * dt * k2)
            k4 = vf(Z + dt * k3)
            Z = Z + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return Z
    h3p, h4p = hh(lambda Z: flow(ring(flow(Z, -1)), +1))
    assert max(abs(v) for v in h3p.values()) < 1e-8
    h4pv = np.array([h4p['h%d%d%d%d0' % q] for q in Q.MONOS4])
    assert np.abs(h4pv - (h4v + c)).max() < 1e-6 * np.abs(c).max()
    assert np.abs(h4pv - h4v).max() > 0.5 * np.abs(c).max()
    # 低レベルの共鳴の近く (k·ν が整数) の調和は補正を 0 にする (発散を避ける)
    c2 = Q.nf_correction_h4(h3v, 1.0 / 3.0, qy)
    assert np.all(np.isfinite(c2))


@pytest.mark.parametrize('mode', ['log', 'lin'])
def test_k1_definition_log_and_lin(mode):
    """U2 の BEND の K1 の定義 ('log' 既定・'lin'、v2 のみ) の両方で、adts()・resonance_terms()・4 次の項が追跡の写像のテイラー展開と合う
    (ASTRID2: h の 3 次 1e-8、4 次の対角 5e-4)。キャッシュのキーに定義を含むので、切り替えで値が変わる ('lin' の xx は 'log' と 1% 以上違う)"""
    old = compat._mode
    compat.set_mode('v2')
    old_def = compat.bend_k1_def()
    try:
        compat.set_bend_k1_def(mode)
        lat = astrid()
        rows = lat.twiss_table()
        rm = T.ring_map3(lat.elements, 0.58e9, rows[0], amp3=1e-4)
        h3 = T.rdt_from_map(rm)
        ht, _ = T.rdt4_from_map(rm)
        names = ('h21000', 'h10110', 'h10020')
        r3 = R.resonance_terms(lat, names, rows=rows)
        for n in names:
            assert abs(r3.total[n] - h3[n]) / abs(h3[n]) < 1e-8
        r4 = Q.rdt4(lat, rows=rows)
        for n in ('h22000', 'h11110', 'h00220'):
            assert abs(r4.total[n] - ht[n]) / abs(ht[n]) < 5e-4
        xx = R.adts(lat, rows=rows).total['xx']
        compat.set_bend_k1_def('lin' if mode == 'log' else 'log')
        lat2 = astrid()
        xx2 = R.adts(lat2, rows=lat2.twiss_table()).total['xx']
        assert abs(xx - xx2) > 0.01 * abs(xx)
    finally:
        compat.set_bend_k1_def(old_def)
        compat.set_mode(old)


def test_octupole_feed_down_to_sextupole_on_displaced_orbit():
    """軌道が x0 ずれた 8 極は 6 極 K2 = K3·x0 のフィードダウンを持つ (H4 = (K3/24)(x⁴ − 6x²y² + y⁴) を x → x0 + x で展開)。collect_nonlinear が
    強さ K3 L x0 の 6 極成分を足す (R6: KEK-LS の δ = −0.002 の ADTS の追跡との差 1.4e-2・2.0e-2・7.8e-2 → 1.8e-5・6.6e-5・3.7e-4 の原因)。
    δ=0 (軌道 0) では足さない"""
    from types import SimpleNamespace
    from pysad.elements import Octupole
    el = Octupole(name='O', L=0.1, K3=500.0)
    rows = [SimpleNamespace(betx=1.0, alfx=0, mux=0, bety=1.0, alfy=0, muy=0, etx=0, etpx=0, s=0.0, x_cod=1e-3),
            SimpleNamespace(betx=1.0, alfx=0, mux=0, bety=1.0, alfy=0, muy=0, etx=0, etpx=0, s=0.1, x_cod=3e-3)]
    e = R.collect_nonlinear([el], rows, 0.002, (0, 0), [])[0]
    assert abs(e.strengths[3] - 500.0 * 0.1) < 1e-12
    assert abs(e.strengths[2] - 500.0 * 0.1 * 2e-3) < 1e-12           # K3 L x0 (x0 = 軌道の平均 2 mm)
    rows0 = [SimpleNamespace(**{**vars(r), 'x_cod': 0.0}) for r in rows]
    e0 = R.collect_nonlinear([el], rows0, 0.0, (0, 0), [])[0]
    assert 2 not in e0.strengths


def test_default_adts_is_r3_calc_and_independent_of_calc_history():
    """V2: adts() の既定は fourth='r3'・optics='calc' (物理的な項だけ)。δ≠0 の CALC を呼んだあと (閉軌道に 1e-10 の雑音が残る) でも δ=0 の値が変わらない
    (orbit 版と R3 版の偏向磁石の分解は K1≠0 で内部の対の項の分け方が違うので、雑音で経路が切り替わらない: bendnl.ORBIT_TOL)。
    δ≠0 の既定の値は追跡の数値の正規形と 4e-3 以内 (ASTRID2・FODO、δ = ±0.002)、fourth='r6'・optics='tracked' の精密版は 1e-3 以内"""
    lat = astrid()
    a0 = R.adts(lat)
    assert a0.fourth == 'r3'
    lat.twiss_table(delta=-0.002)
    R.adts(lat, delta=-0.002, momentum_eV=0.58e9)
    a1 = R.adts(lat)
    for k in ('xx', 'xy', 'yy'):
        assert a1.total[k] == pytest.approx(a0.total[k], rel=1e-6)
    for d in (0.002, -0.002):
        n, _ = nform.adts_at_delta(lat, 0.58e9, d)
        a = R.adts(lat, delta=d, momentum_eV=0.58e9).total
        p = R.adts(lat, delta=d, momentum_eV=0.58e9, fourth='r6', optics='tracked').total
        for k in ('xx', 'xy', 'yy'):
            assert abs(a[k] - n[k]) / abs(n[k]) < 4e-3, (d, k)
            assert abs(p[k] - n[k]) / abs(n[k]) < 1e-3, (d, k)
