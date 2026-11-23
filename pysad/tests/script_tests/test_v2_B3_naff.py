# -*- coding: utf-8 -*-
"""担当 B3 (5): pysad/naff.py (NAFF 風のチューン抽出)。

K64 の MAIN の TRACK の Fourier 解析 (wjfit・wfres・wfbin・wiord・wjdraw、FOURIER フラグ) は FFS から使えず (gs の MAIN の TRACK は SIGSEGV で落ちる、
docs/history/work_main_track.md)、pysad の maintrack は FOURIER を無視する。ここでは追跡の座標列からチューンを取り出す Python の関数として
NAFF (Laskar) を実装し、(1) 既知の周波数の信号、(2) 小さなリング・大きなリングの追跡のチューンを CALC の値と比べる。
根拠: 信号は自分で作った Σ a_k exp(2πi ν_k n) (周波数・振幅・位相が既知)、リングは CALC のチューン (別の計算: 要素の行列の積)。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, ffscmd, naff

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


def _frac(a):
    return (a + 0.5) % 1.0 - 0.5


@pytest.mark.parametrize('N', [256, 1024, 4096])
def test_single_frequency_is_exact(N):
    n = np.arange(N)
    nu = 0.2347182913
    z = 0.7 * np.exp(2j * np.pi * (nu * n + 0.3))
    f, a = naff.naff(z, 1)
    assert abs(f[0] - nu) < 1e-11
    assert abs(abs(a[0]) - 0.7) < 1e-10
    assert abs(np.angle(a[0] * np.exp(-2j * np.pi * 0.3))) < 1e-7


def test_negative_and_aliased_frequencies():
    n = np.arange(1024)
    for nu in (-0.4821, 0.0123, 0.4999, -0.25):
        f, a = naff.naff(np.exp(2j * np.pi * nu * n), 1)
        assert abs(_frac(f[0] - nu)) < 1e-10


def test_two_close_components_and_amplitudes():
    """振幅の差が 14 倍、間隔 0.0041 (= 4.2/N) の 2 成分。周波数 1e-6、振幅 1e-4 (相対)"""
    n = np.arange(1024)
    nu1, nu2 = 0.2347182913, 0.2347182913 + 0.0041
    z = 0.7 * np.exp(2j * np.pi * (nu1 * n + 0.3)) + 0.05 * np.exp(2j * np.pi * (nu2 * n + 1.1))
    f, a = naff.naff(z, 2)
    assert abs(f[0] - nu1) < 1e-6 and abs(f[1] - nu2) < 1e-6
    assert abs(abs(a[0]) / 0.7 - 1) < 1e-6 and abs(abs(a[1]) / 0.05 - 1) < 1e-4


def test_real_signal_with_two_tones():
    n = np.arange(1024)
    x = 0.5 * np.cos(2 * np.pi * (0.31415926 * n) + 0.4) + 0.02 * np.cos(2 * np.pi * 0.1234 * n)
    f, amp, ph = naff.naff_real(x, 2)
    assert abs(f[0] - 0.31415926) < 1e-9 and abs(f[1] - 0.1234) < 1e-7
    assert abs(amp[0] - 0.5) < 1e-8 and abs(amp[1] - 0.02) < 1e-6
    assert abs(ph[0] - 0.4) < 1e-7


def test_noise_error_scales_as_expected():
    """SN 比 14 (振幅 0.7、ノイズ rms 0.05√2) でも周波数は 1e-5 以内 (FFT の分解能 1/N = 1e-3 の 1/100)"""
    n = np.arange(1024)
    nu = 0.2347182913
    z = 0.7 * np.exp(2j * np.pi * nu * n)
    rng = np.random.default_rng(1)
    err = []
    for _ in range(20):
        zn = z + 0.05 * (rng.standard_normal(1024) + 1j * rng.standard_normal(1024))
        err.append(naff.naff(zn, 1)[0][0] - nu)
    assert np.sqrt(np.mean(np.square(err))) < 1e-5


def test_window_is_better_than_fft_peak():
    """窓なしの FFT の放物線補間 (maintrack.wtunem 相当) より精度が 2 桁以上よい"""
    from pysad.maintrack import wtunem
    n = np.arange(1024)
    nu = 0.2347182913
    z = np.zeros((1024, 3), dtype=complex)
    z[:, 0] = np.exp(-2j * np.pi * nu * n)              # wtunem は z = x + i px (正のチューンで exp(−iφ))
    f1, _fa = wtunem(z)
    e_fft = abs(f1[0][0] - nu)
    e_naff = abs(naff.naff(np.exp(2j * np.pi * nu * n), 1)[0][0] - nu)
    assert e_naff < e_fft / 100


_TW = '''FFS;
USE %s; CELL; CALC;
Print["TW ",{Twiss["NX","$$$"],Twiss["NY","$$$"],Twiss["AX",1],Twiss["BX",1],Twiss["AY",1],Twiss["BY",1]}];
Exit;
'''


def _ring_tunes(src, ring, amp=2e-7, N=1024):
    """CALC のチューンと、小さな振幅 (線形) の追跡の座標列から NAFF で出したチューン (x, y)。追跡は放射・RF なし (高速追跡)"""
    from pysad.track6n import track_particles_fast
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(src + '\n' + _TW % ring, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('TW ')][0]
    nx, ny, ax, bx, ay, by = [float(v) for v in line[3:].strip('{}').split(',')]
    els = ffs.lat.elements
    p = ffs.momentum()
    cc = np.array([[amp * math.sqrt(bx), 0.0, 0.7 * amp * math.sqrt(by), 0.0, 0.0, 0.0]])
    rec = []
    for _ in range(N):
        rec.append(cc[0].copy())
        cc = track_particles_fast(els, cc, 1, p, 0.0, False, False, False)
    rec = np.array(rec)
    (fx, fy, _fz), amps = naff.tunes_from_coords(rec, twiss=((ax, bx), (ay, by), (0.0, 1.0)))
    return nx / (2 * math.pi), ny / (2 * math.pi), fx, fy, amps, rec


def test_small_fodo_ring_tune_matches_calc(k64_mode):
    src = '''MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.2 K2=20) SD = (L=0.2 K2=-30);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 SF D2 B D2 SD D2 QD D1);
LINE RING = (ST 12*CELL1);
'''
    nx, ny, fx, fy, amps, rec = _ring_tunes(src, 'RING', amp=1e-7)
    assert abs(_frac(fx - nx)) < 2e-9 and abs(_frac(fy - ny)) < 2e-9
    # 振幅が大きいと6 極でチューンが振幅の 2 乗で動く (1e-5 → 1e-4 で 100 倍)
    d1 = _frac(_ring_tunes(src, 'RING', amp=1e-5)[2] - nx)
    d2 = _frac(_ring_tunes(src, 'RING', amp=1e-4)[2] - nx)
    assert 80 < d2 / d1 < 120
    # Twiss で規格化しなくても (z = x − i px) 1e-6 で合う (振幅は歪む)
    f_raw = naff.tune_xp(rec[:, 0], rec[:, 1])[0]
    assert abs(_frac(f_raw - nx)) < 1e-6


@pytest.mark.parametrize('name,path,ring,tol', [
    ('大きなリング A (3 GeV、603 m、3563 要素)', os.path.join(HERE, '..', 'kek-ls_v03_12.lat'), 'RING', 2e-9),
    ('PF リング (2.5 GeV、187 m)', os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), 'RING', 2e-9),
    ('ASTRID2 (0.58 GeV、45.7 m)', os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING', 2e-9),
])
def test_ring_tune_by_tracking_matches_calc(k64_mode, name, path, ring, tol):
    """1024 周の追跡 (x = 2e-7 m √β) から NAFF で出したチューンが CALC と 2e-9 以内 (実測: A 8e-10、PF 1.4e-10、ASTRID2 5e-11。FFT の放物線補間は 1e-4 の桁)"""
    if not os.path.exists(path):
        pytest.skip(name + ' がない')
    src = open(path, encoding='utf-8').read()
    nx, ny, fx, fy, _a, _r = _ring_tunes(src, ring)
    assert abs(_frac(fx - nx)) < tol, (name, fx, nx)
    assert abs(_frac(fy - ny)) < tol, (name, fy, ny)
