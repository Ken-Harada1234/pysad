# -*- coding: utf-8 -*-
"""v2: FFS の MEASURE / TRCOD / LOGMEASURE と MAIN の TRACK TURNS=-1 (pysad/measure.py)、B41。

K64 (gs) は MEASURE が無出力、TRCOD が SIGSEGV (B42) で参照が取れないので、根拠は解析式と pysad の独立な経路:
  - 初期分布 (ttinit): <x²> = εβ + (ηδ)²、<x px> = -αε + ηη'<δ²>、<px²> = εγ + η'²<δ²>、<z δ> = -σz·DP·AZ/√(1+AZ²)、
    一様分布の dp は rms = DP/√3、TRUNI (単位円内の一様) の横方向は rms = σ/2。
  - tbdecoup: 非結合の emittance = 4×4 行列 Σ·J の固有値 (固有エミッタンス)。入れた R1..R4・β・α・η が戻る。
  - MEASURE: 出口の分散行列 = M Σ0 Mᵀ。M は CALC の Twiss (β, α, ψ, η) から作る 2×2 行列と分散 (追跡とは独立)。
    リングでは Emittance (EMIT) の EMITX・SIGE と Twiss から σx = √(εβ + (η σδ)²)。
  - TRCOD: 重心の軌道 = CALC (INS) の軌道 DX DPX DY DPY、相関 <xδ>/<δ²> = CALC の分散 EX EPX。
  - B41: 軌道の揺らぎは JDX JDPX JDY JDPY × 正規乱数 (v2)。K64 の読み方 (ZPX ZY ZPY DP) は measure.ttinit(b41=True) で再現 (ソースの読みのみ)。
  - JITTER OFF の重心合わせ: v2 は分布の作り方と同じ基準 (K64 の式は px に余分な -α/β·η·dp0 があり、結合も入っていない)。
k64 モードは MEASURE / TRCOD を従来どおり「未対応 (無視)」にする (K64 の動作を変えない)。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat, measure as ms
from pysad.ffscmd import run_sad_string
from pysad.sadrandom import SadRandom

# ---------------------------------------------------------------------------- 格子
TL_HEAD = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 %s) MM=();
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
BEND B1 = (L=1 ANGLE=0.05) ST=(L=0.1 ANGLE=0 K0=1E-3);
LINE TL = (M0 QF D1 B1 D1 QD ST D1 MM B1 D1 QF D1 QD D1);
'''
FFS_TL = r'''
FFS;
USE TL;
INS; CALC;
'''
RING = r'''
MOMENTUM = 2 GEV;
MARK M0 = () MID=();
DRIFT D1 = (L=0.5);
QUAD QF = (L=0.4 K1=0.32) QD=(L=0.4 K1=-0.32);
BEND B1 = (L=1.5 ANGLE=0.1308996938995747);
CAVI CA = (VOLT=5E5 HARM=40);
LINE CELL1 = (QF D1 B1 D1 QD D1 B1 D1);
LINE RING = (M0 12*CELL1 MID 12*CELL1 CA);
FFS;
USE RING;
CELL; CALC;
EMIT;
EMITY=EMITX*0.01; DP=SIGE;
NP=5000; GAUSS;
'''


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def _tl(mark='', body=''):
    return _run(TL_HEAD % mark + FFS_TL + body)


MK = dict(AX=0.5, BX=8.0, AY=-0.3, BY=3.0, EX=0.2, EPX=0.01, EY=0.05, EPY=-0.02, EMITX=1e-8, EMITY=2e-9, DP=1e-3, SIGZ=3e-3,
          AZ=0.4, DX=1e-3, DPX=2e-4, DY=-5e-4, DPY=1e-4)


# ---------------------------------------------------------------------------- ttinit (初期分布)
def test_ttinit_moments_match_analytic(v2_mode):
    n = 100000
    x = ms.ttinit(MK, n, SadRandom(17), gauss=True, jitter=True, dp0=0.0)
    sa, ss = ms.ttstat(x)
    ex, bx, ax, eta, etap, sd = MK['EMITX'], MK['BX'], MK['AX'], MK['EX'], MK['EPX'], MK['DP']
    ey, by, ay, etay, etapy = MK['EMITY'], MK['BY'], MK['AY'], MK['EY'], MK['EPY']
    tol = 4.0 / math.sqrt(n)                       # 分散の推定の統計誤差 √(2/n) の約 3 倍
    assert ss[0, 0] == pytest.approx(ex * bx + (eta * sd) ** 2, rel=tol)
    assert ss[1, 1] == pytest.approx(ex * (1 + ax * ax) / bx + (etap * sd) ** 2, rel=tol)
    assert ss[0, 1] == pytest.approx(-ax * ex + eta * etap * sd * sd, rel=3 * tol)
    assert ss[2, 2] == pytest.approx(ey * by + (etay * sd) ** 2, rel=tol)
    assert ss[3, 3] == pytest.approx(ey * (1 + ay * ay) / by + (etapy * sd) ** 2, rel=tol)
    assert ss[5, 5] == pytest.approx(sd * sd, rel=tol)
    assert ss[4, 4] == pytest.approx(MK['SIGZ'] ** 2, rel=tol)
    assert ss[4, 5] == pytest.approx(-MK['SIGZ'] * sd * MK['AZ'] / math.sqrt(1 + MK['AZ'] ** 2), rel=4 * tol)
    assert ss[0, 5] / ss[5, 5] == pytest.approx(eta, rel=0.02)
    # 重心: 設計値 ± 3σ/√n
    for j, d in enumerate((MK['DX'], MK['DPX'], MK['DY'], MK['DPY'])):
        assert abs(sa[j] - d) < 4.0 * math.sqrt(ss[j, j] / n)


def test_ttinit_uniform_dp_and_truni(v2_mode):
    n = 60000
    mk = dict(MK, EX=0.0, EPX=0.0, EY=0.0, EPY=0.0, AZ=0.0)
    x = ms.ttinit(mk, n, SadRandom(17), gauss=False, trgauss=False, jitter=True)
    sa, ss = ms.ttstat(x)
    assert abs(x[:, 5]).max() <= mk['DP']                                   # 一様分布: |dp| <= DP
    assert ss[5, 5] == pytest.approx(mk['DP'] ** 2 / 3.0, rel=0.02)
    assert ss[0, 0] == pytest.approx(mk['EMITX'] * mk['BX'] / 4.0, rel=0.02)   # 単位円内の一様: <r1²> = 1/4
    assert ss[2, 2] == pytest.approx(mk['EMITY'] * mk['BY'] / 4.0, rel=0.02)


def test_ttinit_quiet_centroid_v2_is_consistent_k64_is_not(v2_mode):
    """JITTER OFF: v2 の重心は分布の作り方と同じ基準 (大きな n の JITTER ON の平均と一致)。K64 の式は px が -ax/bx·ex·dp0 ずれる"""
    dp0 = 2e-3
    v2c = ms.ttinit(MK, 500, SadRandom(17), gauss=True, jitter=False, dp0=dp0).mean(axis=0)
    k64c = ms.ttinit(MK, 500, SadRandom(17), gauss=True, jitter=False, dp0=dp0, center_fix=False).mean(axis=0)
    want = [MK['DX'] + MK['EX'] * dp0, MK['DPX'] + MK['EPX'] * dp0, MK['DY'] + MK['EY'] * dp0, MK['DPY'] + MK['EPY'] * dp0, 0.0, dp0]
    assert v2c == pytest.approx(want, abs=1e-17)
    assert k64c[0] == pytest.approx(want[0], abs=1e-17)
    assert k64c[1] == pytest.approx(want[1] - MK['AX'] / MK['BX'] * MK['EX'] * dp0, abs=1e-17)
    assert k64c[3] == pytest.approx(want[3] - MK['AY'] / MK['BY'] * MK['EY'] * dp0, abs=1e-17)
    big = ms.ttinit(MK, 200000, SadRandom(5), gauss=True, jitter=True, dp0=dp0)
    sa, ss = ms.ttstat(big)
    err = math.sqrt(ss[1, 1] / 200000)
    assert abs(sa[1] - v2c[1]) < 4 * err                     # 真の平均は v2 の基準
    assert abs(sa[1] - k64c[1]) > 20 * err                   # K64 の基準は統計誤差の 20 倍以上ずれる


def test_ttinit_quiet_centroid_with_coupling(v2_mode):
    mk = dict(MK, R1=0.1, R2=0.3, R3=-0.05, R4=0.2, DETR=0.1 * 0.2 + 0.3 * 0.05)
    dp0 = 2e-3
    v2c = ms.ttinit(mk, 400, SadRandom(17), gauss=True, jitter=False, dp0=dp0).mean(axis=0)
    big = ms.ttinit(mk, 200000, SadRandom(5), gauss=True, jitter=True, dp0=dp0)
    sa, ss = ms.ttstat(big)
    for j in range(4):
        assert abs(sa[j] - v2c[j]) < 4 * math.sqrt(ss[j, j] / 200000)


def test_b41_orbit_jitter_reads_jdx_in_v2_and_zpx_in_k64():
    mk = dict(MK, JDX=1e-3, JDPX=2e-4, JDY=3e-4, JDPY=4e-5, ZPX=7e-3, ZY=6e-3, ZPY=5e-3)
    g = SadRandom(17).gauss(4, 1e300)                       # 最初の 4 つの正規乱数 (dx, dpx, dy, dpy の順)
    compat.set_mode('v2')
    try:
        c = ms.ttinit(mk, 100, SadRandom(17), gauss=True, jitter=False).mean(axis=0)
        assert c[:4] == pytest.approx([MK['DX'] + 1e-3 * g[0], MK['DPX'] + 2e-4 * g[1], MK['DY'] + 3e-4 * g[2],
                                       MK['DPY'] + 4e-5 * g[3]], abs=1e-16)
        compat.set_mode('k64')                              # K64 (ソースの読みのみ): rlist(kx+26..29) = ZPX ZY ZPY DP
        c = ms.ttinit(mk, 100, SadRandom(17), gauss=True, jitter=False).mean(axis=0)
        assert c[0] == pytest.approx(MK['DX'] + 7e-3 * g[0], abs=1e-16)
        assert c[1] == pytest.approx(MK['DPX'] + 6e-3 * g[1], abs=1e-16)
        assert c[2] == pytest.approx(MK['DY'] + 5e-3 * g[2], abs=1e-16)
        assert c[3] == pytest.approx(MK['DPY'] + MK['DP'] * g[3], abs=1e-16)      # py の揺らぎが DP·tgauss になる (B41)
    finally:
        compat.set_mode(None)


# ---------------------------------------------------------------------------- tbdecoup
def test_tbdecoup_eigen_emittance_and_recovers_coupling(v2_mode):
    mk = dict(MK, R1=0.1, R2=0.3, R3=-0.05, R4=0.2, DETR=0.1 * 0.2 + 0.3 * 0.05, DX=0.0, DPX=0.0, DY=0.0, DPY=0.0)
    n = 200000
    x = ms.ttinit(mk, n, SadRandom(17), gauss=True)
    sa, ss = ms.ttstat(x)
    p, p1, emx, emy, emu, emv, sigp = ms.tbdecoup(ss)
    eta = ss[5, :4] / ss[5, 5]
    b4 = ss[:4, :4] - np.outer(eta, eta) * ss[5, 5]
    J = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], dtype=float)
    eig = sorted(abs(np.linalg.eigvals(b4 @ J).imag))[::2]
    assert sorted([emu, emv]) == pytest.approx(eig, rel=1e-10)               # 非結合の emittance = 固有エミッタンス
    assert emu == pytest.approx(mk['EMITX'], rel=0.01) and emv == pytest.approx(mk['EMITY'], rel=0.01)
    assert emx < emu and emy > emv                                            # 射影は結合で混ざる
    assert p1[11:15] == pytest.approx([0.1, 0.3, -0.05, 0.2], abs=5e-3)       # R1..R4
    assert (p1[1], p1[2], p1[4], p1[5]) == pytest.approx((0.5, 8.0, -0.3, 3.0), rel=0.02)   # au bu av bv
    assert p1[7:11] == pytest.approx([0.2, 0.01, 0.05, -0.02], abs=2e-3)      # eu epu ev epv
    assert sigp == pytest.approx(1e-6, rel=0.02)


# ---------------------------------------------------------------------------- 書式
def test_ttstat_lines_format(v2_mode):
    assert ms._e103(1.23456e-5) == ' 1.235E-05' and ms._e103(-0.5) == '-5.000E-01' and ms._e103(0.0) == ' 0.000E+00'
    x = ms.ttinit(dict(MK, AZ=0.0), 2000, SadRandom(17), gauss=True)
    L, sa, ss, es = ms.ttstat_lines(x, 'QF.1', {'RAD': False, 'RFSW': True, 'GAUSS': True, 'DP': 1e-3, 'DP0': 0.0, 'GCUT': 1e300})
    assert L[0] == ' Statistics at QF.1:   particles = 2000'
    assert L[1].startswith(' RAD: F, RFSW: T, GAUSS: T, DP =') and ', DP0 =' in L[1] and ', GCUT =' in L[1]
    assert L[2] == ' ' * 12 + '        x     px/p0         y     py/p0         z     dp/p0 '
    assert L[3].startswith('   C of M : ') and len(L[3]) == 12 + 60
    assert [len(l) for l in L[4:10]] == [12 + 10 * k for k in range(1, 7)]    # 下三角
    assert L[10] == ' x-y projected(coupled) parameters:' and L[13] == ' x-y decoupled parameters:'
    assert L[11].startswith(' emitx:') and ' bx:' in L[11] and ' epx:' in L[11]
    assert L[-1].startswith(' sigp/p:') and L[-1].endswith('/sigz')
    assert es is None and len(L) == 20
    assert ms.ttstat_lines(x, ' ')[0] == []                                   # 題名が空白なら表示なし


def test_ttstat_lum_effective_area_and_vertical_moments(v2_mode):
    """lum=.true. の枝 (K64 の MEASURE は呼ばない): 実効面積 ≈ σx·σy (同じ正規分布 2 つの重なり)、垂直の 4・6・8 次のモーメント ≈ 3, 15, 105"""
    rng = np.random.default_rng(11)
    n = 40000
    x = np.zeros((n, 6))
    x[:, 0] = rng.normal(0.0, 1e-3, n)
    x[:, 2] = rng.normal(0.0, 2e-4, n)
    x[:, 1] = rng.normal(0.0, 1e-4, n)
    x[:, 3] = rng.normal(0.0, 1e-4, n)
    L, sa, ss, es = ms.ttstat_lines(x, 'IP', lum=True)
    assert es == pytest.approx(1e-3 * 2e-4, rel=0.1)
    assert L[-2].startswith('   Effective beam area:') and L[-2].endswith('(141*141 bin)')
    mom = [float(v) for v in L[-1].split(':')[1].split()]
    assert mom[0] == pytest.approx(3.0, rel=0.05) and mom[1] == pytest.approx(15.0, rel=0.15) and mom[2] == pytest.approx(105.0, rel=0.4)


# ---------------------------------------------------------------------------- MEASURE (FFS)
def _m2(b1, a1, b2, a2, dpsi):
    c, s = math.cos(dpsi), math.sin(dpsi)
    return np.array([[math.sqrt(b2 / b1) * (c + a1 * s), math.sqrt(b1 * b2) * s],
                     [-((1 + a1 * a2) * s + (a2 - a1) * c) / math.sqrt(b1 * b2), math.sqrt(b1 / b2) * (c - a2 * s)]])


def test_measure_covariance_equals_twiss_transfer_of_initial_sample(v2_mode):
    """MEASURE の出口の分散行列 = M Σ0 Mᵀ。M は CALC の Twiss (追跡と独立) から。エミッタンスと DP を小さくして非線形を消す"""
    ffs, I, out = _tl('EMITX=1E-12 EMITY=2E-13 DP=1E-6', 'NP=1500; GAUSS; r=FFS["MEASURE MM"];')
    r = I.globals['r']
    assert r[0] == 1500.0
    ss = np.array(r[2])
    x0 = ms.Tracker(ffs).x0                                 # FIXSEED: 同じ種なので MEASURE と同じ初期分布
    k = ffs.positions('MM')[0]
    T = lambda key, kk: ffs.twiss_value(key, kk)
    M = np.zeros((5, 5))
    M[:2, :2] = _m2(T('BX', 1), T('AX', 1), T('BX', k), T('AX', k), T('NX', k) - T('NX', 1))
    M[2:4, 2:4] = _m2(T('BY', 1), T('AY', 1), T('BY', k), T('AY', k), T('NY', k) - T('NY', 1))
    M[4, 4] = 1.0
    e1 = np.array([T('EX', 1), T('EPX', 1), T('EY', 1), T('EPY', 1)])
    e2 = np.array([T('EX', k), T('EPX', k), T('EY', k), T('EPY', k)])
    M[:4, 4] = e2 - M[:4, :4] @ e1
    c0 = x0[:, [0, 1, 2, 3, 5]]
    c0 = c0 - c0.mean(axis=0)
    s0 = c0.T @ c0 / 1500
    want = M @ s0 @ M.T
    got = ss[np.ix_([0, 1, 2, 3, 5], [0, 1, 2, 3, 5])]
    scale = np.sqrt(np.outer(np.diag(want), np.diag(want)))
    assert np.abs((got - want) / scale).max() < 2e-4
    # 表示の値: emittance は保存 (線形・非結合)
    p, p1, emx0, emy0, *_ = ms.tbdecoup(ms.ttstat(x0)[1])
    p, p1, emx1, emy1, *_ = ms.tbdecoup(ss)
    assert emx1 == pytest.approx(emx0, rel=1e-4) and emy1 == pytest.approx(emy0, rel=1e-4)


def test_measure_output_and_beam_size_vs_twiss(v2_mode):
    ffs, I, out = _tl('EX=0.2 EPX=0.01 EMITX=1E-8 EMITY=2E-9 DP=1E-3 SIGZ=3E-3',
                      'NP=8000; GAUSS; MEASURE; LOGMEASURE; r=FFS["MEA"]; s=SeedRandom[];')
    assert ' Statistics at $$$:   particles = 8000' in out
    assert 'emitx:' in out and 'sigp/p:' in out
    r = I.globals['r']
    ss = np.array(r[2])
    n = ffs.nlat
    sx = math.sqrt(ffs.twiss_value('BX', n) * 1e-8 + (ffs.twiss_value('EX', n) * 1e-3) ** 2)
    sy = math.sqrt(ffs.twiss_value('BY', n) * 2e-9)
    assert math.sqrt(ss[0, 0]) == pytest.approx(sx, rel=0.03)
    assert math.sqrt(ss[2, 2]) == pytest.approx(sy, rel=0.03)
    p, p1, emx, emy, *_ = ms.tbdecoup(ss)
    assert emx == pytest.approx(1e-8, rel=0.03) and emy == pytest.approx(2e-9, rel=0.03)
    assert p[2] == pytest.approx(ffs.twiss_value('BX', n), rel=0.06) and p[1] == pytest.approx(ffs.twiss_value('AX', n), rel=0.06)
    assert p[7] == pytest.approx(ffs.twiss_value('EX', n), abs=0.03)             # 統計誤差 σxβ/(σδ√N) = 0.006
    # LOGMEASURE: xa ya xxa xya yya (1p,5g15.7)
    ln = [l for l in out.split('\n') if len(l) == 75 and 'E-' in l and ':' not in l]
    v = [float(t) for t in ln[-1].split()]
    assert v == pytest.approx([r[1][0], r[1][2], ss[0, 0], ss[0, 2], ss[2, 2]], rel=1e-6)
    # FIXSEED (既定): MEASURE は乱数の種を進めない。TRPT も元に戻る
    assert I.globals['s'][1] == 17.0
    assert ffs.flags['TRPT'] is False


def test_measure_moveseed_advances_seed(v2_mode):
    ffs, I, out = _tl('EMITX=1E-8 EMITY=2E-9', 'NP=100; MOVESEED; MEASURE; s=SeedRandom[];')
    assert I.globals['s'][1] != 17.0


def test_measure_ring_matched_beam_vs_emittance(v2_mode):
    """リング: CELL CALC の周期解が先頭 MARK に入り、EMIT の EMITX・SIGE の整合した分布を半周追う。
    σx = √(EMITX·BX + (EX·SIGE)²)、σy = √(EMITY·BY) (Twiss と Emittance: 追跡と独立)"""
    ffs, I, out = _run(RING + 'r=FFS["MEASURE MID"];')
    g = I.globals
    k = ffs.positions('MID')[0]
    ss = np.array(g['r'][2])
    emx, sige = float(g['EMITX']), float(g['SIGE'])
    assert emx > 1e-8 and sige > 1e-4
    sx = math.sqrt(emx * ffs.twiss_value('BX', k) + (ffs.twiss_value('EX', k) * sige) ** 2)
    sy = math.sqrt(0.01 * emx * ffs.twiss_value('BY', k))
    assert math.sqrt(ss[0, 0]) == pytest.approx(sx, rel=0.03)
    assert math.sqrt(ss[2, 2]) == pytest.approx(sy, rel=0.03)
    p, p1, e1, e2, *_ = ms.tbdecoup(ss)
    assert e1 == pytest.approx(emx, rel=0.03) and e2 == pytest.approx(0.01 * emx, rel=0.03)
    assert ffs.flags['TRPT'] is False


def test_measure_particle_loss_and_errors(v2_mode):
    src = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3);
DRIFT D1 = (L=1.0);
APERT AP=(DX1=-1E-4 DX2=1E-4 DY1=-1 DY2=1);
LINE NOMARK = (D1 M0 D1);
LINE T2 = (M0 D1 AP D1);
FFS;
USE NOMARK; INS; CALC;
MEASURE;
TRCOD;
USE T2; INS; CALC;
NP=1; MEASURE;
NP=300; EMITX=1E-6; EMITY=1E-9; GAUSS;
r=FFS["MEASURE",6];
TRCOD FIXCOD;
'''
    ffs, I, out = _run(src)
    assert out.count(' Missing entrance condition of beam.') == 2          # 先頭が MARK でない (ttinit.f)
    assert '?MEASURE: NP (number of particles) must be 2 or more' in out
    n = int(I.globals['r'][0])
    assert 2 <= n < 60                                                      # APERT ±0.1 mm、σx = 2.8 mm: 3% ほどが生き残る
    assert f' Statistics at $$$:   particles ={n:5d}' in out
    assert '  trcodf' in out                                                   # FIXCOD は C2 で実装した (v2_changes_C2_measure.md)


# ---------------------------------------------------------------------------- TRCOD
def test_trcod_orbit_and_dispersion_match_calc(v2_mode):
    body = r'''
NP=2000; GAUSS;
n=LINE["LENGTH"];
c0=Table[Twiss[{"DX","DPX","DY","DPY","EX","EPX","EY","EPY"},i],{i,n}];
bx0=Twiss["BX","*"];
TRCOD;
c1=Table[Twiss[{"DX","DPX","DY","DPY","EX","EPX","EY","EPY"},i],{i,n}];
bx1=Twiss["BX","*"];
CALC;
c2=Table[Twiss[{"DX","DPX","DY","DPY","EX","EPX","EY","EPY"},i],{i,n}];
'''
    ffs, I, out = _tl('DX=1E-3 DPX=-2E-4 DY=5E-4 DPY=1E-4 EMITX=1E-16 EMITY=1E-16 DP=1E-5', body)
    g = I.globals
    c0, c1, c2 = np.array(g['c0']), np.array(g['c1']), np.array(g['c2'])
    assert abs(c0[:, 0]).max() > 5e-3                                        # 軌道は mm 程度 (ステアリング ST と入口のずれ)
    assert np.abs(c1[:, :4] - c0[:, :4]).max() < 1e-9                        # 重心の軌道 = CALC の軌道
    assert np.abs(c1[:, 4:6] - c0[:, 4:6]).max() < 2e-4                      # <xδ>/<δ²> = 分散 (EX は 0.6 m 程度)
    assert abs(c0[:, 4]).max() > 0.5
    assert g['bx1'] == g['bx0']                                              # 他の Twiss は変えない
    assert np.array_equal(c2, c0)                                            # 次の CALC で元に戻る
    assert ' Statistics at $$$:   particles = 2000' in out


def test_trcod_centroid_equals_single_particle_tracking(v2_mode):
    """TRCOD の終端の重心 = 入口の重心の 1 粒子の TrackParticles (独立な呼び出し)。位置を指定すると、その後ろの Twiss は変えない"""
    body = r'''
NP=500; GAUSS; DP=0;
n=LINE["LENGTH"];
TRCOD MM;
k=LINE["POSITION","MM"];
a=Twiss[{"DX","DPX","DY","DPY"},k];
b=Twiss[{"DX","DPX","DY","DPY"},n];
p=TrackParticles[{1,{{1E-3},{-2E-4},{5E-4},{1E-4},{0},{0},{1}}},k];
q=Flatten[p[[2,{1,2,3,4}]]];
'''
    ffs, I, out = _tl('DX=1E-3 DPX=-2E-4 DY=5E-4 DPY=1E-4 EMITX=1E-16 EMITY=1E-16', body + 'CALC; c=Twiss[{"DX","DPX","DY","DPY"},n];')
    g = I.globals
    assert np.array(g['a']) == pytest.approx(np.array(g['q']), abs=1e-10)
    assert ' Statistics at MM:' in out
    ffs2, I2, out2 = _tl('DX=1E-3 DPX=-2E-4 DY=5E-4 DPY=1E-4 EMITX=1E-16 EMITY=1E-16', body)
    assert I2.globals['b'] == pytest.approx(g['c'], abs=1e-15)               # MM より後ろは CALC の値のまま


def test_trcod_dp_zero_keeps_dispersion(v2_mode):
    body = r'''
NP=100; GAUSS; DP=0;
n=LINE["LENGTH"];
e0=Twiss["EX","*"]; TRCOD; e1=Twiss["EX","*"];
'''
    ffs, I, out = _tl('EMITX=1E-10 EMITY=1E-10', body)
    assert I.globals['e1'] == I.globals['e0']                                # <δ²>=0 なら EX.. は書かない (trcoda.f: sv(5) /= 0)


# ---------------------------------------------------------------------------- k64 モード
def test_k64_mode_measure_and_trcod_are_ignored_as_before(k64_mode):
    ffs, I, out = _tl('EMITX=1E-8 EMITY=2E-9', 'NP=100; d0=Twiss["DX","*"]; MEASURE; TRCOD; d1=Twiss["DX","*"]; r=FFS["MEASURE"];')
    assert 'Statistics at' not in out
    assert out.count('未対応 (無視)') == 2
    assert I.globals['d1'] == I.globals['d0']
    assert not isinstance(I.globals.get('r'), list)


# ---------------------------------------------------------------------------- MAIN の TRACK TURNS=-1
MAIN_TL = TL_HEAD % 'EX=0.2 EMITX=1E-10 EMITY=2E-11 DP=1E-5 SIGMAZ=3E-3' + r'''
ON GAUSS;
TRACK USE=TL NPART=3000 TURNS=-1;
FFS;
USE TL;
INS; CALC;
'''


def test_main_track_trpt_entrance_exit(v2_mode):
    ffs, I, out = _run(MAIN_TL)
    assert ' Statistics at Entrance:   particles = 3000' in out and ' Statistics at Exit:   particles = 3000' in out
    assert 'Tracking             end:  CPU time' in out
    res = ffs.track_results[0]
    sa0, ss0 = res['entrance']
    sa1, ss1 = res['exit']
    assert math.sqrt(ss0[0, 0]) == pytest.approx(math.sqrt(8 * 1e-10 + (0.2 * 1e-5) ** 2), rel=0.03)
    assert math.sqrt(ss0[4, 4]) == pytest.approx(3e-3, rel=0.03)             # MARK の SIGMAZ
    n = ffs.nlat
    p, p1, emx, emy, *_ = ms.tbdecoup(ss1)
    p0, p10, emx0, emy0, *_ = ms.tbdecoup(ss0)
    assert emx == pytest.approx(emx0, rel=1e-3) and emy == pytest.approx(emy0, rel=1e-3)      # emittance の保存
    # 出口の β: 入口の標本の (β, α) を Twiss の転送行列で運んだ値 (追跡と独立)
    M = _m2(8.0, 0.5, ffs.twiss_value('BX', n), ffs.twiss_value('AX', n), ffs.twiss_value('NX', n))
    b0, a0 = p0[2], p0[1]
    g0 = (1 + a0 * a0) / b0
    b1 = M[0, 0] ** 2 * b0 - 2 * M[0, 0] * M[0, 1] * a0 + M[0, 1] ** 2 * g0
    assert p[2] == pytest.approx(b1, rel=2e-3)
    assert res['alive'] == 3000


def test_main_track_turns_minus1_k64_mode_unchanged(k64_mode):
    """k64 モードは従来どおりリングの経路 (EMIT の正規モード) に入る: 輸送路では周期解が無く、Statistics は出ない"""
    ffs, I, out = _run(MAIN_TL)
    assert 'Statistics at' not in out and 'Tracking             end:' in out
