# -*- coding: utf-8 -*-
"""担当 C の追跡まわり: B10 (tdusrn の NaN)、B33 (最初の tapert の l_track=0)、B35 (SOL 領域の未対応要素)、LOSSMAP (戻りの 2 行)。"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _nums(line):
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', line.split('=', 1)[-1])]


# ---------------------------------------------------------------- B10
def test_b10_tdusrn_nan_draws_one_random_number(v2_mode):
    """v2: anp が NaN・Inf なら光子を出さず、乱数は 1 個 (K64 は NaN で 10000 回まわり 20001 個)。
    根拠: 光子数の期待値が定義できない粒子 (発散した粒子) は放射させない。有限の anp では K64 と同じ乱数列"""
    from pysad.track6n import make_rs, tdusrn, lcg
    for bad in (float('nan'), float('inf'), float('-inf')):
        rs = make_rs(12345)
        assert rs[1] == 1
        r = tdusrn(rs, bad)
        assert r == (0.0, 0.0, 0.0, 0.0)
        assert rs[0] == lcg(12345)
    a = make_rs(777)
    compat.set_mode('k64')
    b = make_rs(777)
    compat.set_mode('v2')
    assert b[1] == 0
    assert tdusrn(a, 0.7) == tdusrn(b, 0.7) and a[0] == b[0]


def test_b10_k64_mode_loops(k64_mode):
    from pysad.track6n import make_rs, tdusrn, lcg
    rs = make_rs(12345)
    tdusrn(rs, float('nan'))
    s = 12345
    for _ in range(20001):
        s = lcg(s)
    assert rs[0] == s
    from pysad.sadrandom import SadRandom, tdusrn as td2
    r = SadRandom()
    r.restore(12345)
    td2(r, float('nan'))
    assert r.seed == s


def test_b10_python_path(v2_mode):
    from pysad.sadrandom import SadRandom, tdusrn as td2, _lcg
    r = SadRandom()
    r.restore(12345)
    assert td2(r, float('nan')) == (0.0, 0.0, 0.0, 0.0)
    assert r.seed == _lcg(12345)


# ---------------------------------------------------------------- B33
B33 = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1);
FFS;
$FORM="18.14";
PageWidth=1999;
USE RING;
LOSSAMPL=1e-3;
NORAD; LOSSMAP;
x0={1e-4,1e-6,2e-2};px0={0,0,1e-4};y0={0,1e-6,0};py0={0,0,0};z0={0,0,0};dp0={0,0,0};fl0={1,1,1};
r=TrackParticles[{1,{x0,px0,y0,py0,z0,dp0,fl0}},1,1,2];
Print["S =",r[[2,1]]," ",r[[2,7]]," ",r[[2,8]]," ",r[[2,9]]];
r=TrackParticles[{1,{x0,px0,y0,py0,z0,dp0,fl0}},1,1,1];
Print["S2 =",r[[2,1]]," ",r[[2,7]]," ",r[[2,8]]," ",r[[2,9]]];
Exit;
'''


def test_b33_first_tapert_drops_particle_outside_lossampl(v2_mode):
    """v2: 最初から LOSSAMPL の外にいる粒子は、セッションの最初の TrackParticles でも最初の tapert で落ちる
    (座標は入射のまま、落ちた周は 1、落ちた要素は始点の 1 つ前 = 最後の要素)。K64 は l_track=0 で印が付かず追跡を続ける"""
    out = _out(B33)
    s1 = _nums([l for l in out.splitlines() if l.startswith('S =')][0])
    s2 = _nums([l for l in out.splitlines() if l.startswith('S2 =')][0])
    # x (3), flag (3), 要素 (3), 周 (3)
    assert abs(s1[2] - 3e-3) < 1e-15 and s1[5] == 0.0 and s1[8] == 37.0 and s1[11] == 1.0      # 座標は 3·LOSSAMPL に頭打ち (limittrack)
    assert abs(s2[2] - 3e-3) < 1e-15 and s2[5] == 0.0 and s2[11] == 1.0
    assert s1[3:5] == [1.0, 1.0] and s1[6:8] == [0.0, 0.0] and s1[9:11] == [3.0, 3.0]      # 生きている粒子: 0 と nend+1


def test_b33_k64_mode_keeps_tracking_on_first_call(k64_mode):
    out = _out(B33)
    s1 = _nums([l for l in out.splitlines() if l.startswith('S =')][0])
    s2 = _nums([l for l in out.splitlines() if l.startswith('S2 =')][0])
    assert abs(s1[2] - 3e-3) > 1e-6      # 最初の呼び出し: 落とす印が 0 で、追跡が続く
    assert abs(s2[2] - 3e-3) < 1e-15     # 2 回目からは最初の tapert で落ちる


# ---------------------------------------------------------------- LOSSMAP (gs の参照: example/05_tracking/flag_LOSSMAP.gs.txt)
LOSSMAP = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 24*CELL);
FFS;
USE RING; CELL; CALC;
NORFSW; NORAD; NOFLUC;
LOSSAMPL = 0.05;
$FORM = "12.6";
x0 = {1e-3, 0.02, 0.025, 0.03, 0.04};
n = Length[x0];
zero = Table[0, {i, n}];
beam = {x0, zero, zero, zero, zero, zero, Table[1, {i, n}]};
r = TrackParticles[{1, beam}, 1, 1, 30];
Print["dim0 =", Dimensions[r[[2]]]];
LOSSMAP;
r = TrackParticles[{1, beam}, 1, 1, 30];
Print["dim1 =", Dimensions[r[[2]]]];
Print["flag =", r[[2,7]]];
Print["elem =", r[[2,8]]];
Print["turn =", r[[2,9]]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_lossmap_rows_match_gs(mode):
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(LOSSMAP)
    finally:
        compat.set_mode(old)
    g = lambda t: _nums([l for l in out.splitlines() if l.startswith(t + ' =')][0])
    assert g('dim0') == [7.0, 5.0] and g('dim1') == [9.0, 5.0]
    assert g('flag') == [1.0, 1.0, 0.0, 0.0, 0.0]
    assert g('elem') == [0.0, 0.0, 242.0, 242.0, 242.0]          # gs: D1.24 (周の終わりの tapert)
    assert g('turn') == [31.0, 31.0, 4.0, 1.0, 1.0]


# ---------------------------------------------------------------- B35
B35 = '''MOMENTUM = 1GEV;
DRIFT D1 = (L=0.5);
SOL SB = (BZ=0.5 BOUND=1) SE = (BOUND=1);
SEXT SX = (L=0.1 K2=1);
MARK ST = ();
LINE BL = (ST D1 SB D1 SX D1 SE D1);
FFS;
USE BL;
NORAD;
r=TrackParticles[{1,{{1e-3},{0},{0},{0},{0},{0},{1}}},8,1,1];
Print["X =",r[[2,1]]];
Exit;
'''


def test_b35_unsupported_element_in_sol_region_is_a_clear_error(v2_mode):
    """K64 は SIGSEGV。v2 は tsol.f が出すはずの文で止める"""
    from pysad.track6 import SolRegionError
    with pytest.raises(SolRegionError) as e:
        _out(B35)
    assert 'are supported within SOL' in str(e.value) and 'SX' in str(e.value)


# ---------------------------------------------------------------- RADTAPER (参照は gs: example/04_emittance_radiation/flag_RADTAPER.gs.txt)
RADTAPER = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
$FORM = "18.12";
PageWidth=1999;
USE RINGQ;
CELL;
CALCULATE;
RFSW; RAD; RADCOD;
e0 = Emittance[];
Print["A = ", (OrbitAtExit/.e0)[[1]], " ", EnergyLossU0/.e0];
RADTAPER;
e1 = Emittance[];
Print["B = ", (OrbitAtExit/.e1)[[1]], " ", (OrbitAtExit/.e1)[[6]], " ", EnergyLossU0/.e1, " ", (Emittances/.e1)[[1]], " ", (Emittances/.e1)[[3]], " ", (Tunes/.e1)[[1]], " ", (Tunes/.e1)[[2]], " ", (Tunes/.e1)[[3]]];
Print["DDP = ", Twiss["DDP", {1, 60, 121, 243}]];
b = {{0},{0},{0},{0},{0},{0},{1}};
Do[b = TrackParticles[{1, b}, 1, 1, 1][[2]]; ,{2}];
Print["T = ", b[[1,1]], " ", b[[2,1]], " ", b[[5,1]], " ", b[[6,1]]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_radtaper_matches_gs(mode):
    """RADTAPER: 磁石の強さを閉軌道の局所の運動量に合わせる (radtaper.py)。gs の値:
    テーパーなし x=-7.887E-5、U0=1563263.2。テーパーあり x=-2.198E-7、δ=3.5E-9、U0=1564904.0、εx=1.7372E-7、εz=1.9077E-5、
    チューン -.303633 -.094058 -.027619、DDP {0, 3.9113E-4, 2.6068E-4, 3.5E-9}、2 周の追跡 {4.8092E-5, 1.4477E-5, .., -3.693E-5, 8.0067E-5}。
    gs の tcod は収束判定が緩い (radtaper では 100 倍) ので、軌道の残り (0.2 μm) と追跡は 1〜2% の一致"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(RADTAPER)
    finally:
        compat.set_mode(old)
    g = lambda t: _nums([l for l in out.splitlines() if l.startswith(t + ' =')][0])
    a, b, d, t = g('A'), g('B'), g('DDP'), g('T')
    # v2 (A1 v2next 項目 7): tradke の ddpx の符号が追跡 (tradkn) と同じになり、U0 などが gs と 1e-6 程度 (U0 は相対 2e-6) 違う。k64 は gs と同じ
    tu = 1e-7 if mode == 'k64' else 5e-6
    assert abs(a[0] / -7.887e-5 - 1) < 1e-3 and abs(a[1] / 1563263.2 - 1) < tu
    assert abs(b[0]) < 3e-7 and abs(b[1]) < 1e-8                 # 軌道が設計軌道に戻る
    assert abs(b[2] / 1564904.0 - 1) < tu
    if mode == 'k64':             # エミッタンスは v2 で放射の式が直っている (担当 A) ので、gs との比較は k64 だけ
        assert abs(b[3] / 1.7372e-7 - 1) < 1e-4 and abs(b[4] / 1.9077e-5 - 1) < 1e-4
    for u, v in zip(b[5:8], [-.303633, -.094058, -.027619]):
        assert abs(u - v) < (2e-6 if mode == 'k64' else 5e-5)       # v2: tradke の ddpx の符号 (A1 項目 7)・tradke の微分と RADTAPER の係数 (P2 項目 3・4)
    assert d[0] == 0.0 and abs(d[1] / 3.9113e-4 - 1) < 1e-4 and abs(d[2] / 2.6068e-4 - 1) < 1e-4 and abs(d[3]) < 1e-8
    if mode == 'k64':             # v2 は追跡の閉軌道が閉じる (A1 項目 7) ので、閉じない gs の追跡の値 (1 周 1e-5 m 動く) とは比べない (閉じ方は test_v2_A1_next.py)
        for u, v in zip(t, [4.8092e-5, 1.4477e-5, -3.693e-5, 8.0067e-5]):
            assert abs(u / v - 1) < 0.06


# ---------------------------------------------------------------- ResponseMatrix の TRPT (参照は gs)
RM = '''MOMENTUM = 1GEV;
DRIFT D1 = (L=1);
QUAD QF = (L=0.3 K1=0.5) QD = (L=0.3 K1=-0.5);
BEND ST1 = (L=0.1 K0=0) ST2=(L=0.1 K0=0);
CAVI AC = (L=1 VOLT=1e8 FREQ=2.856e9 PHI=-1.5707963267948966);
MONI M1=() M2=();
MARK ST = (BX=5 BY=5);
LINE BL = (ST ST1 D1 QF D1 M1 AC QD D1 ST2 AC D1 QF D1 M2);
FFS;
USE BL;
TRPT; INS; CALC;
$FORM="18.12";
PageWidth=1999;
n=LINE["LENGTH"];opt={Table[Twiss["*",k],{k,1,n}]};
st={LINE["POSITION","ST1"],LINE["POSITION","ST2"]};
Print["RM =",ResponseMatrix[LINE["POSITION","M2"],st,st+1,opt,0,"X"]];
Print["R1 =",ResponseMatrix[LINE["POSITION","M1"],st,st+1,opt,0,"X"]];
NOTRPT;
Print["RN =",ResponseMatrix[LINE["POSITION","M2"],st,st+1,opt,0,"X"]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_response_matrix_trpt_matches_gs(mode):
    """pgrmat.f: TRPT では応答に sqrt(gammab(stt+1)/gammab(mon)) が掛かる (加速による断熱減衰)。1 → 1.2 GeV の加速
    v2 (B330): 空洞の端の 1/β (1 GeV で 1.3e-7) のため応答が gs と相対 4e-9 (RM[0] で 1.9e-8) 違う。k64 は gs と 1e-9 以内、v2 は 1e-7 以内"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(RM)
    finally:
        compat.set_mode(old)
    tol = 1e-9 if mode == 'k64' else 1e-7
    g = lambda t: _nums([l for l in out.splitlines() if l.startswith(t + ' =')][0])
    for a, b in zip(g('RM'), [-5.180683831097, -1.896985774280]):
        assert abs(a - b) < tol
    for a, b in zip(g('R1'), [-1.459575549598, 0.0]):
        assert abs(a - b) < tol
    for a, b in zip(g('RN'), [-5.675154846643, -1.896985774280]):
        assert abs(a - b) < tol
    assert abs(g('RM')[0] / g('RN')[0] - math.sqrt(1.0 / 1.2)) < 2e-5        # 運動量の比の平方根



# ---------------------------------------------------------------- スピン追跡 (POL + RAD) の損失 (参照は gs)
SPIN_LOSS = 'MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="16.12";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nLOSSAMPL=0.02; LOSSMAP;\nx0={1e-4,0.012,0.03};px0={0,0,0};y0={0,1e-4,0};py0={0,0,0};z0={0,0,0};dp0={0,1e-3,0};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0,0},{0.1,0.2,0.3},{1,1,1}}}, 1, 1, 6];\nPrint["D ",Dimensions[r[[2]]]];\nPrint["F ",r[[2,9]]," ",r[[2,10]]," ",r[[2,11]]];\nPrint["X ",r[[2,1]]];\nPrint["S ",r[[2,7]]," ",r[[2,8]]];\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0,0},{0.1,0.2,0.3},{1,1,1}}}, 1, 1, 6];\nPrint["F2 ",r[[2,9]]," ",r[[2,10]]," ",r[[2,11]]];\nPrint["X2 ",r[[2,1]]];\nExit;\n'


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_spin_tracking_loss_matches_gs(mode):
    """POL + RAD の TrackParticles (track6 の経路) でも LOSSAMPL・LOSSDZ の tapert で粒子が落ち、LOSSMAP の 2 行が付く。
    gs: 粒子 2 (x0=12 mm) は 1 周目に要素 74 で落ち x=.022037636775、粒子 1 は生き残る。2 回目の呼び出しでは粒子 3 (x0=30 mm > LOSSAMPL) が
    最初の tapert (要素 74 = 前の追跡の最後、1 周目) で落ちる。1 回目の粒子 3 は K64 の癖 (B33: l_track=0 で落ちない) で gs は要素 60、
    pysad の k64 は要素 15 (track6 の経路は癖を近似)、v2 は最初の tapert で落とす"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(SPIN_LOSS)
    finally:
        compat.set_mode(old)
    g = lambda t: _nums([l for l in out.splitlines() if l.startswith(t + ' ')][0].replace(t + ' ', '=', 1))
    f, x, s, f2, x2 = g('F'), g('X'), g('S'), g('F2'), g('X2')
    # v2 は放射の式の修正 (担当 A) で 1e-7 ほど違い、さらに四極の本体の 4 次の合成 (B226) で、12 mm の粒子 (7 周) の x が gs と 6.5e-7 違う
    # (K64 の Strang 分割 (2 次) の誤り。下の方式 0 の確認で 2e-7 以内に戻る)
    tol = 1e-12 if mode == 'k64' else 2e-6
    assert f[:2] == [1.0, 0.0] and f[3:5] == [0.0, 74.0] and f[6:8] == [7.0, 1.0]
    assert abs(x[0] - -2.715851550e-4) < tol and abs(x[1] - .022037636775) < tol
    # スピンの角: v2 (C2) は刻みごとの歳差の角を厳密な θ にした (K64 は 2 sin(θ/2): 7 周で 1.6e-3 rad 小さい)
    stol = 10 * tol + 1e-11 if mode == 'k64' else 3e-3
    assert abs(s[0] - 2.650424039172) < stol and abs(s[1] - 2.456033665635) < stol
    assert f2 == [1.0, 0.0, 0.0, 0.0, 74.0, 74.0, 7.0, 1.0, 1.0]
    assert abs(x2[2] - 0.03) < 1e-15
    assert f[2] == 0.0 and f[8] == 1.0
    if mode == 'v2':
        assert f[5] == 74.0 and abs(x[2] - 0.03) < 1e-15
    if mode == 'v2':
        # 四極を K64 の方式 (方式 0) にすると、従来の許容 (2e-7) で gs と一致する = 差は四極の積分器の次数の違い
        from pysad import quadbody
        old_s = quadbody.set_scheme(0)
        compat.set_mode(mode)
        try:
            out0 = _out(SPIN_LOSS)
        finally:
            compat.set_mode(old)
            quadbody.set_scheme(old_s)
        x0 = _nums([l for l in out0.splitlines() if l.startswith('X ')][0].replace('X ', '=', 1))
        assert abs(x0[0] - -2.715851550e-4) < 2e-7 and abs(x0[1] - .022037636775) < 2e-7
        assert abs(x[1] - x0[1]) > 1e-7                    # 4 次の合成が効いている
