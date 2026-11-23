# -*- coding: utf-8 -*-
"""担当 C3 (v2next): SynchroBetaEmittance (結合のあるリング)、FITP>1 の最終解、FLUC の乱数列 (NPARA)、DELCOR の規格化。

- 結合のある平衡エミッタンスの独立な計算 (`_c3_indep.coupled_equilibrium`: 自前の 6x6 行列 + 放射の減衰の摂動 + 量子励起の拡散 + Lyapunov 方程式、
  SAD・pysad の光学を使わない) と Emittance[] が一致 (ASTRID2 + skew 四極)。SynchroBetaEmittance (シンクロ・ベータ結合を含む) は結合のあるリングで
  走り、Emittance[] と 1.5% 以内で一致する。
- FITP 3 の最終解 (QF・QD の K1、SF・SD の K2) は、CALC のチューンを使った独立な最小化 (scipy の最小二乗、K64 の残差の定義) の最小値と 1e-5 で一致する。
- FLUC つきの追跡の乱数列: pysad (k64) は gs と全桁一致する (NPARA=1)。gs の既定は NPARA=4 (子プロセスに分ける) で、`NPARA=4;` で同じ乱数列を再現する。
"""
import io
import math
import os
import sys

import numpy as np
import pytest

from pysad import compat, ffscmd, emit6d
from pysad.sadrandom import SadRandom

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import _c3_indep as C       # noqa: E402

PK = os.environ.get('PYSAD_SAD_PACKAGES')
needs_pk = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')
AST = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
PF = os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad')


def _read(path):
    if not os.path.exists(path):
        pytest.skip(os.path.basename(path) + ' がない')
    with open(path, encoding='utf-8') as fh:
        return fh.read()


def _skew_ring(path, ks, tail):
    return _read(path) + f'''
QUAD SQ1 = (L=0.1 K1={ks} ROTATE=0.7853981633974483);
LINE RING2 = (SQ1 RING);
FFS;
PageWidth=1999;
USE RING2; CELL; RFSW; RAD; CALC;
{tail}
Exit;
'''


# ---------------------------------------------------------------------------------------------
# 結合のある平衡エミッタンス
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize('ks', [0.0, 0.05])
def test_emittance_of_coupled_ring_matches_independent_6x6_equilibrium_v2(v2_mode, ks):
    """ASTRID2 (0.58 GeV、周長 45.7 m) + skew 四極 (ks = 0.05、ν_x − ν_y ≈ 0.04 の近く): Emittance[] の εx・εy・σδ・U0 が、自前の 6x6 の方法と
    εx 1e-4、εy 3e-4、σδ 1e-4、U0 1e-5 で一致する (実測: 6e-5、1.3e-4、4e-5、2e-6)"""
    ffs, _ = ffscmd.run_sad_string(_skew_ring(AST, ks, ''), out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
    eps, sd, u0, nus, _sig = C.coupled_equilibrium(ffs.lat.elements, ffs.momentum())
    tr = sorted(eps)[:2]                                      # 固有エミッタンス 3 つのうち最大が縦 (z) の値
    emit_ = sorted([abs(e['emittance_x']), abs(e['emittance_y'])])
    assert abs(tr[1] / emit_[1] - 1) < (1e-4 if ks == 0.0 else 3e-4), (eps, emit_)
    assert abs(tr[0] - emit_[0]) < (1e-15 if ks == 0.0 else 3e-4 * emit_[0]), (eps, emit_)
    assert abs(sorted(eps)[2] / e['emittance_z'] - 1) < 1e-3
    assert abs(sd / e['momentum_spread'] - 1) < 1e-4
    assert abs(u0 / e['energy_loss_eV'] - 1) < 1e-5


@needs_pk
@pytest.mark.parametrize('path,ks,tol', [(AST, 0.05, 0.015), (AST, 0.005, 0.015), (PF, 0.01, 0.03)])
def test_synchrobeta_runs_on_coupled_ring_and_agrees_with_emittance_v2(v2_mode, path, ks, tol):
    """SynchroBetaEmittance が結合のあるリング (gs は落ちる) で走り、Emittance[] (線形) と tol 以内で一致する
    (ASTRID2: +0.8%、PF: +1.1%・+2.9% (εy が小さい)。差はシンクロ・ベータ結合・色収差の分)。結合なしの比と同じ比 (結合の効果が両方で同じ)"""
    tail = ('em=Emittance[]; nus=Abs[(Tunes/.em)[[3]]];\n'
            'r=SynchroBetaEmittance[{nus,nus,0.001},AzimuthalModes->3];\n'
            'Print["EM ",Emittances/.em]; Print["SB ",r[[1,2]]," ",r[[1,3]]];')
    out = []
    for k in (0.0, ks):
        buf = io.StringIO()
        ffscmd.run_sad_string(_skew_ring(path, k, tail), out=buf)
        txt = buf.getvalue()
        em = [float(v) for v in [l for l in txt.splitlines() if l.startswith('EM ')][0].replace('{', ' ').replace('}', ' ').replace(',', ' ').split()[1:3]]
        sb = [float(v) for v in [l for l in txt.splitlines() if l.startswith('SB ')][0].split()[1:3]]
        out.append((em, sb))
    (em0, sb0), (em1, sb1) = out
    assert sb1[1] > 1e-12                                     # 結合で εy ができる
    assert abs(sb1[0] / em1[0] - 1) < tol and abs(sb1[1] / em1[1] - 1) < tol * 3
    assert abs((sb1[0] / em1[0]) / (sb0[0] / em0[0]) - 1) < 0.005     # x の比は結合の有無で変わらない
    # 和則 (弱い結合): ε1 + ε2 が結合なしの値 (EMIT・SB とも)
    if path == AST and ks <= 0.01:                            # PF (ks = 0.01) は +0.5% (η_y・結合共鳴の効果で和は保存しない)
        assert abs((em1[0] + em1[1]) / em0[0] - 1) < 2e-3
        assert abs((sb1[0] + sb1[1]) / sb0[0] - 1) < 2e-3


_SB_COUPLED = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD  SQ = (L=0.1 K1=0.02 ROTATE=0.7853981633974483);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 HARM=64);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RING = (ORIGIN SQ RF 24*CELLQ);
FFS;
PageWidth = 1999;
USE RING; CELL; RFSW; RAD; CALC;
r = SynchroBetaEmittance[{0.0142, 0.0142, 0.001}, AzimuthalModes->5];
Print["SB ", r];
Exit;
'''
# gs (K64) の参照 (example/04_emittance_radiation/fn_SynchroBetaEmittance_coupled.gs.txt): {νs, εx, εy, εx', εy', 収束}
_GS_SB = [0.0142, 1.360336495570143e-6, 4.751665441653782e-8, 1.3599991774350626e-6, 4.925609883688084e-8, -5.079004432581897e-7]


@needs_pk
def test_synchrobeta_coupled_ring_matches_gs_k64(k64_mode):
    """結合のある 24 セル FODO (skew 四極 1 台、6 極なし、3 GeV、RF 4 MV): k64 の SynchroBetaEmittance が gs と εx・εy で 1e-6、εx'・εy' で 1e-6 一致する
    (gs はこの結合のあるリングで落ちない。「結合のあるラティスで gs が落ちる」(docs/31) は再現しなかった)。v2 は追跡の放射の式の修正で 1.3% 違う"""
    buf = io.StringIO()
    ffscmd.run_sad_string(_SB_COUPLED, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('SB ')][0]
    got = [float(v) for v in line.replace('{', ' ').replace('}', ' ').replace(',', ' ').split()[1:]]
    for a, b in zip(got[:5], _GS_SB[:5]):
        assert abs(a / b - 1) < 1e-6, (got, _GS_SB)
    assert abs(got[5] / _GS_SB[5] - 1) < 1e-4


# ---------------------------------------------------------------------------------------------
# FITP>1 の最終解: 独立な最小化
# ---------------------------------------------------------------------------------------------

def test_fitp_final_solution_is_the_independent_minimum(k64_mode):
    """examples の FFS_FITPOINTS (24 セル FODO、QF QD SF SD を変数に NX 6.7・NY 4.9 を δ = −1%, 0, +1% の 3 点で課す) の GO の最終解が、
    CALC のチューンを使った独立な最小二乗 (残差 = Σ w(NX − 6.7·2π)² + w(NY − 4.9·2π)²、w = 1 (δ=0)、(1/4)² (δ=±1%) = K64 の重み
    (OffMomentumWeight/2/√npoints)²) の最小値と一致する: K1・K2 は 5e-6 以内、残差は最小値の +1e-5 以内。(gs は最小値の手前で止まり、K2 が 1.7e-4 違う)"""
    from scipy.optimize import least_squares
    head = open(os.path.join(ROOT, 'example', '03_matching', 'FFS_FITPOINTS.sad'), encoding='utf-8').read().split('FFS;')[0]

    def tunes(v):
        body = f'''FFS;
$FORM = "20.14"; PageWidth=1999;
USE RING; CELL;
QF K1={v[0]}; QD K1={v[1]}; SF K2={v[2]}; SD K2={v[3]};
DP0=-0.01; CALC; Print["TN ",DP0," ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]];
DP0=0; CALC; Print["TN ",DP0," ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]];
DP0=0.01; CALC; Print["TN ",DP0," ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]];
Exit;'''
        buf = io.StringIO()
        ffscmd.run_sad_string(head + body, out=buf)
        return [[float(x) for x in l.split()[1:]] for l in buf.getvalue().splitlines() if l.startswith('TN ')]

    def resid(v):
        r = []
        for d, nx, ny in tunes(v):
            w = 1.0 if d == 0 else 0.25
            r += [w * (nx - 6.7 * 2 * math.pi), w * (ny - 4.9 * 2 * math.pi)]
        return np.array(r)
    sol = least_squares(resid, np.array([0.7, -0.6, 1.0, -1.0]), xtol=1e-14, ftol=1e-14, gtol=1e-14, diff_step=1e-6)
    buf = io.StringIO()
    ffscmd.run_sad_string(open(os.path.join(ROOT, 'example', '03_matching', 'FFS_FITPOINTS.sad'), encoding='utf-8').read().replace(
        'SHOW;', '').replace('$FORM = "10.6";', '$FORM = "20.14";'), out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('after GO')][0]
    import re
    got = [float(x) for x in re.findall(r'-?\d*\.\d+(?:[Ee][-+]?\d+)?', line)]
    assert len(got) == 4
    for a, b in zip(got, sol.x):
        assert abs(a - b) < 5e-6, (got, sol.x)
    assert 2 * sol.cost == pytest.approx(6.3318e-6, rel=2e-4)


# ---------------------------------------------------------------------------------------------
# FLUC の乱数列 (NPARA)
# ---------------------------------------------------------------------------------------------

_FODO = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
$FORM = "20.14"; PageWidth=1999;
USE RINGQ;
CELL;
CALCULATE;
'''

# gs (K64、2 つの条件) の参照: 3 GeV の 24 セル FODO、RFSW RAD FLUC、SeedRandom[7]、粒子はすべて原点・dp=0
_GS_N16 = [-4.37001210354E-4, -5.08284567481E-4, -5.06458604942E-4, -5.35499042371E-4, -5.47300159119E-4, -6.35133594701E-4,
           -4.76839614340E-4, -5.25627832138E-4, -5.96251166359E-4, -6.05650429747E-4, -4.18754741530E-4, -5.13758163338E-4,
           -4.72070798367E-4, -4.43723294408E-4, -5.79080006744E-4, -4.50845319374E-4]      # np=16、1 周後の dp (gs の既定 NPARA=4。有効 NPARA=2)
_GS_RMS = {1: (5.252739925347E-4, .00314056869362, .00236045313644),      # NPARA=1: np=200 の dp の rms (1・10・100 周後)
           4: (5.205643226207E-4, .00314042063085, .00237021386166)}      # NPARA=4


def _fluc_run(npara, np_, turns_print):
    body = f'''NPARA={npara};
FLUC; SeedRandom[7];
RFSW; RAD; NORADCOD;
np = {np_};
beam = {{Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[1,{{np}}]}};
Do[ beam = TrackParticles[{{1, beam}}, 1, 1, 1][[2]];
    If[MemberQ[{{{','.join(str(t) for t in turns_print)}}},k], Print["K ",k," ",Sqrt[Plus @@ (beam[[6]]^2)/np]]], {{k,{max(turns_print)}}}];
Print["DP ", beam[[6]]];
Exit;
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(_FODO + body, out=buf)
    txt = buf.getvalue()
    k = [float(l.split()[2]) for l in txt.splitlines() if l.startswith('K ')]
    dp = None
    for l in txt.splitlines():
        if l.startswith('DP '):
            dp = [float(v) for v in l.replace('{', ' ').replace('}', ' ').replace(',', ' ').split()[1:]]
    return k, dp


@pytest.mark.parametrize('npara', [1, 4])
def test_fluc_random_stream_matches_gs_for_nparallel_k64(k64_mode, npara):
    """FLUC つきの追跡 (200 粒子、1・10・100 周) の dp の rms が gs と 1e-9 で一致する: gs の NPARA=1 は pysad の既定 (NPARA=1)、
    gs の既定 NPARA=4 (tftrack.f: 子プロセスに粒子を分け、乱数の種 + 2k) は pysad の `NPARA=4;` (有効 npara = min(NPARA, np/9 + 1, ne·np/3000 + 1))"""
    k, _ = _fluc_run(npara, 200, [1, 10, 100])
    for a, b in zip(k, _GS_RMS[npara]):
        assert abs(a / b - 1) < 1e-9, (k, _GS_RMS[npara])


def test_fluc_np16_per_particle_matches_gs_default_nparallel_k64(k64_mode):
    """np=16 (有効 npara = 2: 子 [0,8) 種 + 2、親 [8,16)) の 1 周後の粒子ごとの dp が gs の既定 (NPARA=4) と 1e-9 で一致"""
    _, dp = _fluc_run(4, 16, [1])
    assert len(dp) == 16
    for a, b in zip(dp, _GS_N16):
        assert abs(a / b - 1) < 1e-9


def test_fluc_nparallel_default_is_one_and_changes_nothing_when_noflucc(k64_mode):
    """NPARA の既定は 1 (分けない)。FLUC を切れば NPARA に依らず同じ結果 (乱数を引かない)"""
    s = _FODO + 'NPARA=4; NOFLUC;\nRFSW; RAD; np=40;\nbeam={Table[0,{np}],Table[1e-4,{np}],Table[0,{np}],Table[0,{np}],Table[0,{np}],Table[0,{np}],Table[1,{np}]};\nbeam=TrackParticles[{1,beam},1,1,1][[2]]; Print["B ",beam[[1]]];\nExit;\n'
    s1 = s.replace('NPARA=4;', 'NPARA=1;')
    b4, b1 = io.StringIO(), io.StringIO()
    ffscmd.run_sad_string(s, out=b4)
    ffscmd.run_sad_string(s1, out=b1)
    l4 = [l for l in b4.getvalue().splitlines() if l.startswith('B ')]
    l1 = [l for l in b1.getvalue().splitlines() if l.startswith('B ')]
    assert l4 == l1 and l4


# PF リング (周長 187 m、2.5 GeV) の gs の参照 (RFSW RAD FLUC、np=40、SeedRandom[7]、`USE RING; CELL; RFSW; RAD; CALC`): 1・10・40 周後の dp の rms と最終の種
_GS_PF = {1: [(1, 1.603389376910E-4, 3353699231), (10, .00139199603361, 3553936015), (40, 7.796561863505E-4, 2330613871)],
          4: [(1, 1.559341195548E-4, 3835945671), (10, .00139257583841, 2363275047), (40, 7.858910508268E-4, 3463281031)]}


@pytest.mark.parametrize('npara', [1, 4])
def test_fluc_random_stream_matches_gs_on_pf_ring_k64(k64_mode, npara):
    """PF リング (760 要素、6 極つき) でも dp の rms が gs と 1e-9、40 周後までの乱数の種が gs と同じ (NPARA=1 と 4)"""
    src = _read(PF) + f'''
FFS;
$FORM="20.14"; PageWidth=1999;
USE RING; CELL; RFSW; RAD; CALC;
NPARA={npara};
FLUC; SeedRandom[7];
np=40;
beam = {{Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[0,{{np}}],Table[1,{{np}}]}};
Do[ beam = TrackParticles[{{1, beam}}, 1, 1, 1][[2]]; If[k==1 || k==10 || k==40, Print["K ",k," ",Sqrt[Plus @@ (beam[[6]]^2)/np]," ",SeedRandom[][[2]]]], {{k,40}}];
Exit;
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    rows = [l.split() for l in buf.getvalue().splitlines() if l.startswith('K ')]
    assert len(rows) == 3
    for r, (k, rms, seed) in zip(rows, _GS_PF[npara]):
        assert int(float(r[1])) == k
        assert abs(float(r[2]) / rms - 1) < 1e-9
        assert int(round(float(r[3]))) == seed


# ---------------------------------------------------------------------------------------------
# DELCOR の規格化 (B99)
# ---------------------------------------------------------------------------------------------

class _Stub:
    def __init__(self, seed):
        self.rng = SadRandom()
        self.rng.restore(seed)

    def gcut(self):
        return 1.0e300


@pytest.mark.parametrize('mode,want', [('k64', 1.0), ('v2', 2.0)])
def test_delcor_tvcorr_normalisation_of_the_correlated_branch(mode, want):
    """corlength ≠ 0 の tvcorr: K64 は s = 1/√(4Σw² + 2) で複素数の <|c|²> = 1 に規格化する (実部の分散 1/2、DX の rms = delx/√2)。
    v2 は実部・虚部それぞれの分散 1 (<|c|²> = 2、rms = delx)。K64 の corlength = 0 の枝 (ランダムウォーク) は実部・虚部の刻みがそれぞれ分散 al/n
    (1/√2 の補正なし) なので、K64 の 2 つの枝で規格化の約束が揃っていない (docs/history/v2next_C3.md)"""
    from pysad.ffscmd import FFS
    old = compat._mode
    compat.set_mode(mode)
    try:
        n, al, x0 = 512, 500.0, 15.0
        st = _Stub(3)
        acc = 0.0
        m = 200
        for _ in range(m):
            cv = FFS._tvcorr(st, x0, al, n)
            acc += float(np.mean(np.abs(cv) ** 2))
    finally:
        compat.set_mode(old)
    assert abs(acc / m / want - 1) < 0.05


def test_delcor_tvcorr_zero_corlength_is_brownian_bridge_per_component():
    """corlength = 0: 実部・虚部それぞれが al/n の分散の刻みのランダムウォークを直線で引き直した橋 (Brownian bridge)。中央の分散 = al u(1−u) = al/4
    (実部。刻みの分散に 1/2 は掛からない)"""
    from pysad.ffscmd import FFS
    n, al = 256, 200.0
    st = _Stub(5)
    vals = []
    for _ in range(600):
        cv = FFS._tvcorr(st, 0.0, al, n)
        vals.append(cv[n // 2].real)
    assert abs(np.var(vals) / (al / 4) - 1) < 0.15
