# -*- coding: utf-8 -*-
"""v2 (C2): DELCOR の corlength≠0 の規格化。

K64 (tvcorr.f): 1/sqrt(4Σw²+2) = 複素数 c の絶対値の 2 乗平均を 1 にする。DX は Re(c)·delx、DY は Im(c)·dely なので、DX の rms は delx/√2。
v2: 1/sqrt(2Σw²+1) = 実部・虚部それぞれの分散を 1 にする。DX の rms = delx、DY の rms = dely。
根拠:
  - 実部と虚部は独立な 2 つの実数の場 (DX 用と DY 用) を入れる器で、|c| という量はどこにも使われない。delx・dely は別々に指定する。
  - corlength=0 の枝 (ランダムウォーク) は成分ごとの規格化 (実部の 1 歩の分散 = al/n、つまり 1 m あたり 1) で、1/2 は入らない。
  - 相関長 ≫ 周長の極限では全要素が同じ値 delx·g (g: 標準ガウス乱数) だけずれる。その rms は delx であるべき。
  - 相関関数の形 exp(−Δ²/(2 x0²)) は規格化の定数によらず保たれる。
"""
import io
import math
import re

import numpy as np

from pysad.sadrandom import SadRandom


class _Stub:
    def __init__(self, seed):
        self.rng = SadRandom()
        self.rng.restore(seed)

    def gcut(self):
        return 1.0e300


def _samples(x0, al, n, m, seed=17):
    from pysad.ffscmd import FFS
    st = _Stub(seed)
    return np.array([FFS._tvcorr(st, x0, al, n) for _ in range(m)])


def test_rms_of_each_component_is_one(v2_mode):
    # 相関長が格子間隔程度・中間・周長より長い、の 3 通り。標本 400 × 512 点 (独立な自由度は al/x0 × 400 程度)
    for x0, tol in ((3.0, 0.02), (15.0, 0.04), (5000.0, 0.08)):
        c = _samples(x0, 500.0, 512, 400)
        assert abs(math.sqrt(np.mean(c.real**2)) - 1.0) < tol, x0
        assert abs(math.sqrt(np.mean(c.imag**2)) - 1.0) < tol, x0
        assert abs(np.mean(c.real * c.imag)) < tol                 # DX と DY は無相関


def test_long_correlation_limit_is_a_common_gaussian_shift(v2_mode):
    c = _samples(1.0e6, 500.0, 128, 2000)
    assert np.max(np.abs(c - c[:, :1])) < 1e-9                     # リング全体が同じ値
    g = c[:, 0].real
    assert abs(g.std() - 1.0) < 0.05 and abs(g.mean()) < 0.07      # その値は標準ガウス乱数 (rms 1)
    assert abs(np.mean(np.abs(g) < 1.0) - 0.6827) < 0.035


def test_correlation_function_shape_is_kept(v2_mode):
    n, al, x0, m = 512, 500.0, 15.0, 300
    c = _samples(x0, al, n, m)
    acc = np.zeros(n, dtype=complex)
    for cv in c:
        acc += np.fft.ifft(np.abs(np.fft.fft(cv))**2) / n          # 巡回自己相関 <c_j c*_(j+k)>
    r = acc / m
    assert abs(r[0].real - 2.0) < 0.1                              # <|c|²> = 2 (実部・虚部それぞれ 1)
    for k in (1, 5, 10, 20, 30):
        d = k * al / n
        assert abs(r[k].real / r[0].real - math.exp(-d * d / (2 * x0 * x0))) < 0.05
        assert abs(r[k].imag / r[0].real) < 0.02


def test_k64_mode_keeps_k64_normalization(k64_mode):
    c = _samples(15.0, 500.0, 512, 400)
    assert abs(math.sqrt(np.mean(np.abs(c)**2)) - 1.0) < 0.04      # K64: <|c|²> = 1
    assert abs(math.sqrt(np.mean(c.real**2)) - math.sqrt(0.5)) < 0.03


def test_random_walk_branch_is_per_component_and_same_in_both_modes():
    from pysad import compat
    old = compat._mode
    try:
        compat.set_mode('v2')
        a = _samples(0.0, 100.0, 128, 1500)
        compat.set_mode('k64')
        b = _samples(0.0, 100.0, 128, 1500)
    finally:
        compat.set_mode(old)
    assert np.array_equal(a, b)                                    # corlength = 0 の枝はモードによらずビット単位で同じ
    # 両端を引き戻したランダムウォーク: 実部の分散 = s (1 − s/al) [m] (成分ごと。1/2 は入らない)。中点で al/4 = 25
    v = np.mean(a[:, 64].real**2)
    assert abs(v / 25.0 - 1.0) < 0.1


SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
$FORM = "18.10";
USE RING;
CELL;
Do[
  FFS["DELCOR 5 1e-4 2e-4 Q*"];
  Print["DX ", LINE["DX","Q*"]];
  Print["DY ", LINE["DY","Q*"]], {150}];
Exit;
'''


def _run_rms():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    o = {'DX': [], 'DY': []}
    for key, body in re.findall(r'(DX|DY) (\{[^}]*\})', buf.getvalue()):      # 長いリストは行を折り返して表示される
        o[key] += [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', body)]
    return o


def test_delcor_command_rms_of_dx_is_delx(v2_mode):
    o = _run_rms()
    assert len(o['DX']) == len(o['DY']) == 150 * 24
    rx = math.sqrt(np.mean(np.square(o['DX'])))
    ry = math.sqrt(np.mean(np.square(o['DY'])))
    # 周長 62.4 m、相関長 5 m、150 回: 独立な自由度は約 150 × 62.4/(5√(2π)) ≈ 750 → rms の統計誤差は約 2.6 %
    assert abs(rx / 1e-4 - 1.0) < 0.08
    assert abs(ry / 2e-4 - 1.0) < 0.08


def test_delcor_command_k64_rms_is_delx_over_sqrt2(k64_mode):
    o = _run_rms()
    rx = math.sqrt(np.mean(np.square(o['DX'])))
    ry = math.sqrt(np.mean(np.square(o['DY'])))
    assert abs(rx / 1e-4 - math.sqrt(0.5)) < 0.06
    assert abs(ry / 2e-4 - math.sqrt(0.5)) < 0.06
