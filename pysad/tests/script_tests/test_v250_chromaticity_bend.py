# -*- coding: utf-8 -*-
"""v2.5.0: `Lattice.chromaticity()` が偏向磁石のあるリングで追跡・gs と一致する (docs/30 B325)。

旧い実装は、δ ごとの要素の線形行列 (閉軌道 = 0) の 1 周のトレースの微分 + 6 極の解析式だった。
分散で動いた閉軌道が偏向磁石の端面・本体 (K1 の有無) を通るときの δ 依存 (実効の端面角 e + η'δ による垂直・水平の収束など) を
落とし、偏向磁石のあるリングで水平 −1.456 対 −0.780、垂直 −1.180 対 −0.879 と追跡・gs から大きくずれた。
新しい実装は twiss_table (SAD の CALC と同じ閉軌道 + 光学) の位相進みの δ = ±1e-4 の中心差分。

基準 (1) 追跡: 閉軌道 (Newton) のまわりの 1 周写像のヤコビアンから出したチューンの差分 (pysad の光学を使わない)。
基準 (2) gs (K64, 2026-10-22 ビルド) の `Twiss["NX",-1]` の DP0 = ±1e-4 の中心差分 (CELL; NORFSW)。値は下の定数に保存した。
ラティス: 3 GeV、周長 41.6 m、8 セル (QF D1 B D2 QD D1 B D2 [+ SF/SD])、偏向磁石 16 台 (ANGLE = 2π/16、L = 1.2 m)。
厚い 6 極 (B326): k64 の CALC の厚い 6 極 (tthine、2 次の分割) は追跡 (tthin、4 次) と dν/dδ が 4e-4 (相対) ずれる (K64 の再現。gs と 1e-7 で一致)。
v2 は CALC も追跡と同じ 4 次の分割 (力勾配型の補正キック) にしたので、追跡と 1e-5 以内 (本リングで 1.5e-9) で一致し、gs (K64) とは 8e-5〜4e-4 離れる。
限界: 6 極なしの 4 リング・薄い 6 極は 1e-6〜1e-5 で一致。
1 つのエネルギー・1 つのリング構成だけ。RF・放射なし。δ = ±1e-4 の 1 次の中心差分 (2 次以上の δ 依存は見ない)。
FFS の `CHROMATICITY` コマンドは QUAD・SEXT・MULT だけの閉形式 (ヘルプの simplest formula) で BEND を数えず、値は 4π·ξ の単位。
ここでは、四極だけのリングで TOTAL = 4π·ξ、偏向磁石があると TOTAL/4π が dν/dδ と違うことだけを確かめる (コマンドは変えていない)。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat, lattice as lat_mod, rdttrack
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')

ANG = 2 * math.pi / 16
H = 1e-4

# 名前: (BEND の属性、セルの中身)。gs の参照 (dνx/dδ, dνy/dδ) = (NX, NY の差分 / 2h / 2π)
CASES = {
    'drift': ('DRIFT B=(L=1.2)', 'QF D1 B D2 QD D1 B D2', (-1.2143680927, -1.1798992371)),
    'bend': ('BEND B=(L=1.2 ANGLE=%.15g)' % ANG, 'QF D1 B D2 QD D1 B D2', (-0.7797363826, -0.8786795941)),
    'bend_edge': ('BEND B=(L=1.2 ANGLE=%.15g E1=0.1 E2=0.1)' % ANG, 'QF D1 B D2 QD D1 B D2', (-0.6367595612, -0.8285175391)),
    'bend_k1': ('BEND B=(L=1.2 ANGLE=%.15g K1=0.02)' % ANG, 'QF D1 B D2 QD D1 B D2', (-0.8275066962, -0.9456916018)),
    'bend_k1_edge_sext': ('BEND B=(L=1.2 ANGLE=%.15g E1=0.05 E2=0.05 K1=0.02)' % ANG, 'QF D1 B D2 SF QD D1 B D2 SD',
                          (-21.1586913997, -1.4480861266)),
}


def _ring(name):
    bend, cell, _ = CASES[name]
    src = ('MOMENTUM = 3GEV;\nDRIFT D1=(L=0.8) D2=(L=0.3);\n' + bend + ';\n'
           'QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);\nSEXT SF=(L=0.2 K2=2) SD=(L=0.2 K2=-3);\nMARK ST=();\n'
           'LINE CELL=(' + cell + ');\nLINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);\nFFS;\nUSE RING;\n')
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def _tracking(ffs, h=H):
    """追跡の 1 周写像 (閉軌道 + ヤコビアン) の小数部のチューンの δ 差分"""
    nus = []
    for dd in (h, -h):
        f = rdttrack.one_turn(ffs.lat.elements, ffs.momentum(), dd)
        z = np.zeros(4)
        for _ in range(12):
            r = f(z)[0] - z
            z = z - np.linalg.solve(rdttrack.jacobian(f, z) - np.eye(4), r)
        J = rdttrack.jacobian(f, z)
        cx, cy = np.trace(J[:2, :2]) / 2, np.trace(J[2:, 2:]) / 2
        mx = math.atan2(math.sqrt(1 - cx ** 2) * np.sign(J[0, 1]), cx)
        my = math.atan2(math.sqrt(1 - cy ** 2) * np.sign(J[2, 3]), cy)
        nus.append((mx / (2 * math.pi), my / (2 * math.pi)))
    return ((nus[0][0] - nus[1][0]) / (2 * h), (nus[0][1] - nus[1][1]) / (2 * h))


@pytest.fixture(params=['k64', 'v2'])
def mode(request):
    old = compat._mode
    compat.set_mode(request.param)
    try:
        yield request.param
    finally:
        compat.set_mode(old)


@pytest.mark.parametrize('name', list(CASES))
def test_chromaticity_matches_tracking_and_gs(name, mode):
    ffs = _ring(name)
    xi = np.array(ffs.lat.chromaticity())
    trk = np.array(_tracking(ffs))
    gs = np.array(CASES[name][2])
    # 追跡のヤコビアンの差分と 5e-6 (相対。ヤコビアンの数値微分の誤差 ~6e-7)
    # 6 極つき: k64 は gs の CALC の線形化 (2 次の分割。pysad も同じ) と追跡の厳密なヤコビアンが 4e-4 (相対) ずれる (h = 2e-5〜3e-4 で同じ)。
    # v2 (B326) は CALC も 4 次の分割にしたので追跡と 1e-5 以内 (ここでは 1.5e-9)、gs (K64) とは離れる
    sext = 'sext' in name
    rt = 1e-3 if (sext and mode == 'k64') else (1e-5 if sext else 5e-6)
    assert np.allclose(xi, trk, rtol=rt, atol=1e-6), (name, xi, trk)
    if sext and mode == 'v2':
        assert np.abs(xi / gs - 1).max() > 5e-5, (name, xi, gs)         # gs の 2 次の分割から動いた (水平 7.7e-5、垂直 4.2e-4)
        assert np.abs(xi / gs - 1).max() < 2e-3
        return
    # gs と (gs の出力は 10 桁)
    assert np.allclose(xi, gs, rtol=1e-7, atol=1e-8), (name, xi, gs)


def test_old_formula_missed_the_bend_part():
    """旧い実装 (閉軌道 = 0 の線形行列のトレースの微分 + 6 極の解析式) は、偏向磁石のあるリングで追跡と大きく違った (再現)。"""
    ffs = _ring('bend')
    els = ffs.lat.elements
    mp, mm = np.eye(6), np.eye(6)
    for el in els:
        mp = el.transfer_matrix(+H) @ mp
        mm = el.transfer_matrix(-H) @ mm
    m0 = np.eye(6)
    for el in els:
        m0 = el.transfer_matrix(0.0) @ m0
    cx0 = (m0[0, 0] + m0[1, 1]) / 2
    sx0 = math.sqrt(1 - cx0 ** 2) * np.sign(m0[0, 1])
    old_x = (-1.0 / sx0) * ((mp[0, 0] + mp[1, 1]) - (mm[0, 0] + mm[1, 1])) / (2 * H) / 2 / (2 * math.pi)
    new_x = ffs.lat.chromaticity()[0]
    assert abs(old_x - new_x) > 0.5          # −1.456 対 −0.780
    assert abs(new_x - CASES['bend'][2][0]) < 1e-6


def test_chromaticity_restores_state():
    ffs = _ring('bend_k1_edge_sext')
    lat = ffs.lat
    lat.twiss_table(delta=0.0)
    cod0 = lat._cod_prev.copy()
    rf0 = lat_mod.RFSW
    a = lat.chromaticity()
    b = lat.chromaticity()
    assert a == b
    assert np.array_equal(lat._cod_prev, cod0)
    assert lat_mod.RFSW == rf0


def _chro_total(name):
    bend, cell, _ = CASES[name]
    src = ('MOMENTUM = 3GEV;\nDRIFT D1=(L=0.8) D2=(L=0.3);\n' + bend + ';\n'
           'QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);\nSEXT SF=(L=0.2 K2=2) SD=(L=0.2 K2=-3);\nMARK ST=();\n'
           'LINE CELL=(' + cell + ');\nLINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);\nFFS;\nUSE RING;\nCELL;\nCALC;\nCHROMATICITY;\n')
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    tot = [ln for ln in buf.getvalue().splitlines() if ln.strip().startswith('TOTAL')][-1]
    return np.array([float(v) for v in tot.split()[1:3]])


def test_chromaticity_command_is_quad_sext_formula_in_4pi_units(mode):
    """四極だけのリングで TOTAL/(4π) = dν/dδ (薄肉近似の積分 = 厚い四極の厳密な ∫Kβ ds の形)。偏向磁石があると BEND の分だけ違う。
    gs (CELL; CALC; CHROMATICITY) の TOTAL: drift (−15.260, −14.827)、bend (−10.613, −14.827)。"""
    t0 = _chro_total('drift')
    xi0 = np.array(_ring('drift').lat.chromaticity())
    assert np.allclose(t0 / (4 * math.pi), xi0, rtol=2e-3)
    assert np.allclose(t0, [-15.260, -14.827], atol=2e-3)
    t1 = _chro_total('bend')
    assert np.allclose(t1, [-10.613, -14.827], atol=2e-3)
    xi1 = np.array(_ring('bend').lat.chromaticity())
    assert abs(t1[1] / (4 * math.pi) - xi1[1]) > 0.2        # 垂直: コマンドは −1.180、実際は −0.879 (BEND の分)


# ---------------------------------------------------------------------------------------------
# 6 極つきのずれの原因 (docs/30 B326): CALC の厚い 6 極 (tthine、2 次の分割) と追跡 (tthin、4 次) の差。v2 は CALC も 4 次にして直した
# ---------------------------------------------------------------------------------------------
def _sext_ring(K2, L, nslice=0):
    """3 GeV、8 セル、QF D1 B D2 SF QD D1 B D2。SF は厚い 6 極 (L)、nslice > 0 なら L=0 の薄い 6 極 nslice 個 (中点則) に分ける"""
    if nslice:
        d = L / nslice
        sf = ('DRIFT DS=(L=%.15g);\nDRIFT DH=(L=%.15g);\nSEXT KS=(L=0 K2=%.15g);\n' % (d, d / 2, K2 / nslice)
              + 'LINE SF=(DH ' + ' '.join(['KS'] + ['DS KS'] * (nslice - 1)) + ' DH);\n')
    else:
        sf = 'SEXT SF=(L=%.15g K2=%.15g);\n' % (L, K2)
    src = ('MOMENTUM = 3GEV;\nDRIFT D1=(L=0.8) D2=(L=0.3);\nBEND B=(L=1.2 ANGLE=%.15g);\n'
           'QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);\n' % ANG + sf +
           'MARK ST=();\nLINE CELL=(QF D1 B D2 SF QD D1 B D2);\nLINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);\nFFS;\nUSE RING;\n')
    return run_sad_string(src, out=io.StringIO())[0]


def _trk(ffs):
    return np.array(_tracking(ffs))


def test_sextupole_thin_calc_equals_tracking(mode):
    """薄い (L=0) 6 極では、CALC の線形化と追跡のヤコビアンが 2e-5 (相対) 以内 (6 極の強さ 2 と 4)"""
    for k2 in (2.0, 4.0):
        ffs = _sext_ring(k2, 0.0)
        assert np.allclose(ffs.lat.chromaticity(), _trk(ffs), rtol=2e-5, atol=1e-5), k2


def test_sextupole_thick_calc_error_scales_L2_and_tracking_is_exact():
    """k64: 厚い 6 極の差 (CALC − 追跡) は L² に比例 (L = 0.2 と 0.1 で 4 倍) し、薄い 6 極に分けた極限 (厳密) は追跡に一致、CALC は離れる。
    原因: CALC の tthine は L/6・キック/2・2L/3・キック/2・L/6 (2 次)。追跡の tthin は 4 次"""
    old = compat._mode
    compat.set_mode('k64')
    try:
        d = {}
        for L in (0.2, 0.1):
            ffs = _sext_ring(2.0, L)
            d[L] = np.array(ffs.lat.chromaticity()) - _trk(ffs)
        assert abs(d[0.2][0]) > 1e-3 and abs(d[0.2][0] / d[0.1][0] - 4.0) < 0.2
        ex = np.array(_sext_ring(2.0, 0.2, nslice=40).lat.chromaticity())      # 厳密 (40 分割の中点則。80 分割との差は 2e-6)
        trk = _trk(_sext_ring(2.0, 0.2))
        calc = np.array(_sext_ring(2.0, 0.2).lat.chromaticity())
        assert np.allclose(trk, ex, rtol=0, atol=1e-4)
        assert abs(calc[0] - ex[0]) > 1e-3
    finally:
        compat.set_mode(old)


def test_sextupole_thick_calc_equals_tracking_v2():
    """v2 (B326): 厚い 6 極の CALC も追跡と同じ 4 次の分割なので、dν/dδ が追跡と 1e-5 (相対) 以内で一致する (L = 0.2・0.1、K2 = 2。L = 0.4・K2 = 4)。
    薄い 6 極 40 個の極限 (厳密) との差は 4 次の誤差 (L = 0.2 で水平 3e-6) だけ"""
    old = compat._mode
    compat.set_mode('v2')
    try:
        for K2, L in ((2.0, 0.2), (2.0, 0.1), (4.0, 0.4)):
            ffs = _sext_ring(K2, L)
            assert np.allclose(ffs.lat.chromaticity(), _trk(ffs), rtol=1e-5, atol=1e-5), (K2, L)
        ex = np.array(_sext_ring(2.0, 0.2, nslice=40).lat.chromaticity())
        assert np.allclose(_sext_ring(2.0, 0.2).lat.chromaticity(), ex, rtol=0, atol=1e-4)
    finally:
        compat.set_mode(old)


def _nums(t, tag):
    import re
    m = re.search(r'^' + tag + r' \{([^}]*)\}', t, re.M | re.S)
    return np.array([float(v) for v in m.group(1).replace('\n', ' ').split(',')])


def test_k64_tracking_with_sextupole_matches_gs_trackparticles():
    """6 極つきのリング (厚い SF) の TrackParticles 1 周 (5 粒子、δ = 1e-3 を含む) が gs と一致 (k64: 1e-16)。gs の出力は scripts/ に保存"""
    import os
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts')
    src = open(os.path.join(here, 'v250_chromaticity_sext__ring.sad'), encoding='utf-8').read()
    gs = open(os.path.join(here, 'v250_chromaticity_sext__ring.gs.txt'), encoding='utf-8', errors='replace').read()
    old = compat._mode
    compat.set_mode('k64')
    try:
        buf = io.StringIO()
        run_sad_string(src, out=buf)
    finally:
        compat.set_mode(old)
    for tag in 'XPYQ':
        assert np.abs(_nums(buf.getvalue(), tag) - _nums(gs, tag)).max() < 1e-15, tag
