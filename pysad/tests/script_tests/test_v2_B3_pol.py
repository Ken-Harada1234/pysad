# -*- coding: utf-8 -*-
"""担当 B3 (1・6): v2 の Polarization のシンクロトロン側帯 (スピンの位相変調) の項 (spin.sideband_depolarization)。

導出 (docs/history/v2next_B3.md §2): スピンのずれ α の 1 周の写像 α' = e^{iμ}(e^{iP·u} α + G·u)。位相の変調 P·u は主にエネルギー δ (aγ δ)。
不変な場 α*(u) = −Σ_k e^{−iμk} e^{−iS_k} G·M^k u と、その勾配から脱偏極率 ½Σ_i<G_i B_i G_i^†> (u はガウス分布)。変調を切ると線形の SLIM の式に一致する。
根拠: (1) 変調 OFF で線形の式と 1e-5、(2) 振幅を小さくした極限 (放射の跳びと分布を f 倍、f=0.01) で線形の式 (rate/f が 2% 以内)、
(3) 線形の軌道 + 厳密なスピン回転の Monte Carlo (example/16_spin_polarization/pol_sideband_spin_mc.py、3000 粒子 × 18000 周、傾きの統計誤差 ~2%): 実際の振幅 (f=1) で
    傾き 4.97e-8 /周 (式 4.46e-8: −10%)、f=0.01 で 2.37e-8 (線形 2.32e-8、式 2.33e-8)。
条件: 12 セル FODO (周長 45.6 m、偏向磁石 15 度 × 12 = 180 度、四極の縦ずれ ±3 mm、RF 調和数 64)、4.2 GeV・RF 12 MV (νs = 0.030、スピンチューンの端数 −0.236、σδ = 1.84e-3)。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, spin


@pytest.fixture(autouse=True, scope='module')
def _emitdiv_default_1():
    """参照値は EMITDIV = 1 で作った。v2.4.0 から v2 の EMITDIV の既定は 4 だが、v2 の Polarization の式 (線形の脱偏極率) は EMITDIV に強く動く
    (4.2 GeV・12 MV の FODO: 1 → 2.318e-8、4 → −3%、16 → −16%。B137 の刻み依存の疑い) ので、スピンの試験は 1 に固定する"""
    old = os.environ.get('PYSAD_EMITDIV_DEFAULT')
    os.environ['PYSAD_EMITDIV_DEFAULT'] = '1'
    try:
        yield
    finally:
        if old is None:
            os.environ.pop('PYSAD_EMITDIV_DEFAULT', None)
        else:
            os.environ['PYSAD_EMITDIV_DEFAULT'] = old

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
PST = 8.0 / (5.0 * math.sqrt(3.0))


def _lat(E, volt, dy=3e-3):
    return f'''MOMENTUM = {E}GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY={dy}) QDE = (L=0.3 K1=-0.6 DY=-{dy});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={volt}e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''


def _setup(E, volt, dy=3e-3):
    compat.set_mode('v2')
    try:
        ffs, _ = ffscmd.run_sad_string(_lat(E, volt, dy), out=io.StringIO())
        p = ffs.momentum()
        e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
        els = ffs.lat.elements
        cod = np.array(e['orbit_entrance'])
        trf0 = e.get('trf0', 0.0)
        ph = emit6d._Phys(p, emit6d._emit_circ(els))
        old = emit6d.POL_MODE
        emit6d.POL_MODE = True
        secs = []
        c = cod.copy()
        try:
            for el in els:
                st = emit6d.tturne([el], c, ph, trf0=trf0, rfsw=True, radcod=True)
                secs.append((st.T + st.D, st.B.copy(), st.srot[:, :3].copy(), st.srot[:, 3:9].copy()))
                c = st.cod.copy()
        finally:
            emit6d.POL_MODE = old
        sl = spin.slim_depolarization(els, ph, cod, trf0=trf0, rfsw=True, radcod=True)
    finally:
        compat.set_mode('k64')
    return dict(e=e, secs=secs, sl=sl)


@pytest.fixture(scope='module')
def ring():
    return _setup(4.2, 12)


def test_sideband_without_modulation_equals_linear_slim(ring):
    sb = spin.sideband_depolarization(ring['secs'], ring['sl']['n0'], ring['e']['beam_matrix'], modulate=False)
    assert abs(sb['rate'] / ring['sl']['rate_el'] - 1) < 1e-5            # 要素の出口の N での線形の率 (P1: 'rate' は放射の刻みごとの N の積分で 0.2% 違う)
    assert abs(sb['rate'] / ring['sl']['rate'] - 1) < 5e-3
    assert abs(sb['spread'] / ring['e']['spin']['spin_spread'] - 1) < 5e-3       # ½<|α*|²> = 入口の n の広がり


def test_sideband_small_amplitude_limit_is_linear(ring):
    """放射の跳びとビーム分布を f=0.01 倍 (変調指数が 1/10) にすると、率/f は線形の式に 2% で近づく。f=1 では変調で 1.9 倍"""
    s0, secs, n0 = ring['e']['beam_matrix'], ring['secs'], ring['sl']['n0']
    lin = ring['sl']['rate']
    f = 0.01
    r = spin.sideband_depolarization([(T, B * f, R, D) for (T, B, R, D) in secs], n0, s0 * f)['rate'] / f
    assert abs(r / lin - 1) < 0.02
    r1 = spin.sideband_depolarization(secs, n0, s0)['rate']
    assert 1.8 < r1 / lin < 2.05


def test_sideband_matches_linear_orbit_exact_spin_monte_carlo(ring):
    """線形の軌道 + 厳密なスピン回転の MC (3000 粒子、傾きの統計誤差 ~2%): f=1 で 4.97e-8、f=0.01 で 2.37e-8 /周。式は 15% 以内 (結合の u の 2 次の項は入れていない)"""
    s0, secs, n0 = ring['e']['beam_matrix'], ring['secs'], ring['sl']['n0']
    r1 = spin.sideband_depolarization(secs, n0, s0)['rate']
    assert abs(r1 / 4.97e-8 - 1) < 0.15, r1
    f = 0.01
    r = spin.sideband_depolarization([(T, B * f, R, D) for (T, B, R, D) in secs], n0, s0 * f)['rate'] / f
    assert abs(r / 2.37e-8 - 1) < 0.05, r


def test_polarization_includes_sideband_v2(v2_mode):
    """Emittance[] の spin 結果: 平衡偏極度は側帯の項を入れた率を使い (線形より低い)、depol_rate_linear が線形の値を残す"""
    ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12), out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    sp = e['spin']
    assert sp['depol_sideband'] is True
    assert abs(sp['depol_rate_linear'] / 2.40756e-8 - 1) < 0.01      # v2.4.0 D (B165): 2.318e-8 → 2.342e-8 (+1.0%。D の y 列の本体の B∥)、D3 (B230: 四極の本体の B∥): → 2.4076e-8 (+2.8%)
    assert abs(sp['depol_rate'] / 1.09e-7 - 1) < 0.05                    # P1: ベータトロンのクロマティックな変調を入れた値 (B3 は 4.46e-8)
    assert sp['depol_chromatic'] is True
    sdamp, srel = sp['sdamp'], sp['sdamp_rel']
    assert abs(sp['equpol'][0] - PST * (sdamp - sp['kinetic']) / (srel + sp['depol_rate'])) < 1e-12     # P1: 運動論的偏極の項
    old_c = spin.CHROMATIC
    spin.CHROMATIC = False
    try:
        e1 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    finally:
        spin.CHROMATIC = old_c
    assert abs(e1['spin']['depol_rate'] / 4.46e-8 - 1) < 0.05            # B3 の値 (軌道は線形)
    old = spin.SIDEBAND
    spin.SIDEBAND = False
    try:
        e0 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    finally:
        spin.SIDEBAND = old
    assert e0['spin']['depol_rate'] == e0['spin']['depol_rate_linear']
    assert e0['spin']['equpol'][0] > sp['equpol'][0]


def test_sideband_is_small_at_low_energy_and_planar_ring_is_unchanged(v2_mode):
    """平面リング (ずれなし): 脱偏極率 0、P = ST。3 GeV・1.5 MV (νs = 0.012) のずれつきリングは側帯の効果が +13% (線形 5.52e-11 → 6.26e-11)"""
    ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12, 0.0), out=io.StringIO())
    sp = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)['spin']
    assert sp['depol_rate'] == pytest.approx(0.0, abs=1e-20)
    assert abs(sp['equpol'][0] - PST) < 1e-7
    ffs, _ = ffscmd.run_sad_string(_lat(3.0, 1.5), out=io.StringIO())
    sp = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)['spin']
    assert abs(sp['depol_rate_linear'] / 5.8504e-11 - 1) < 0.02      # v2.4.0 D (B165): 5.52e-11 → 5.850e-11 (+6%。垂直の D への本体の B∥。以前はこの値が「平面リングで変わらない」基準だった。B∥ なしとの差は spin_bpar=False で戻る)
    assert 1.05 < sp['depol_rate'] / sp['depol_rate_linear'] < 1.25


def test_sideband_k64_unchanged(k64_mode):
    ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12, 0.0), out=io.StringIO())
    sp = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)['spin']
    assert math.isnan(sp['equpol'][0])                # K64 は NaN (B55)
    assert 'depol_rate' not in sp


def test_astrid2_mpw_polarization_after_sideband_v2(v2_mode):
    """項目 6: ASTRID2 + MPW (逆向きの極、0.58 GeV) の P・τ を独立な式 (符号つき・絶対値の ∮ρ^-3) と比べる。側帯の項を入れた後も、
    平面のリング (ずれなし、n は鉛直) は脱偏極率 0 で P・τ は A2 の値 (P = 0.7019、τ = 21915 s)。独立な式 P = 0.70163、τ = 21886 s と 2e-3"""
    path = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
    if not os.path.exists(path):
        pytest.skip('ASTRID2_v01_02_id.sad がない')
    src = open(path, encoding='utf-8').read() + '\nFFS;\nUSE RINGMPW;\nCELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    sp = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)['spin']
    RE, LC, C, ME = 2.8179403262e-15, 3.8615926796e-13, 299792458.0, 0.51099895e6
    i3s = i3a = ctot = 0.0
    for el in ffs.lat.elements:
        ctot += el.L
        if type(el).__name__ == 'Bend' and el.L > 0:
            h = (el.ANGLE + getattr(el, 'K0', 0.0)) / el.L
            i3s += h ** 3 * el.L
            i3a += abs(h) ** 3 * el.L
    gam = math.sqrt(1 + (ffs.momentum() / ME) ** 2)
    tau = 1.0 / (5 * math.sqrt(3) / 8 * RE * LC * C * gam ** 5 * i3a / ctot)
    pst = PST * i3s / i3a
    assert abs(sp['equpol'][0] / pst - 1) < 2e-3
    assert abs(sp['taup'] / tau - 1) < 2e-3
    assert sp['depol_rate'] == pytest.approx(0.0, abs=1e-20)
