"""v2.4.0 担当 A: パーサ・評価器・FIT まわりの gs (K64) との照合。参照値は gs で取った。
(1) ユーザー定義の関数名 At がメンバー参照 `obj@name` と衝突しない、(2) NaN の FitFunction は条件として効かない、
(3) 末尾の数 `FIT .. 3.` が語として評価される、(4) USE のたびに DP が MARK (または 0.01) に戻る、
(5) FitFunction の行は gs では最大フィット型 (不等式) の扱い (k64)、(6) 繰り返しの USE → マッチングが同じ結果になる。"""
import io
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

RING = '''MOMENTUM=3GEV;
MARK M0=();
DRIFT D1=(L=1);
BEND B=(L=1 ANGLE=0.3926990817);
QUAD QF=(L=0.3 K1=0.5) QD=(L=0.3 K1=-0.5);
LINE C=(M0 QF D1 B D1 QD D1 B D1 QF);
LINE R=(C C C C);
'''


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _line(txt, tag):
    for l in txt.splitlines():
        if l.startswith(tag + ' '):
            return l[len(tag) + 1:].split()
    raise AssertionError(tag + ' not printed:\n' + txt[-1500:])


def test_user_function_At():
    r = _run('FFS;\nAt[x_] := x + 100; Print["A ", At[1]];\nExit;\n')
    assert _line(r, 'A') == ['101']


@pytest.mark.parametrize('ff', ['FitFunction := 0./0.;', 'FitFunction := 0.;', ''])
def test_nan_fitfunction_is_ignored(ff, k64_mode):
    # gs: NaN の FitFunction は、FitFunction を置かないのと反復の履歴まで同じ (QF の K1 の gs の値 1.3827321232325072)
    src = RING + 'FFS;\nUSE R;\nCELL; INS;\nFREE QF QD;\nfit NX 1.9;\nfit NY 1.1;\n' + ff + \
        '\nGO;\nPrint["K ", LINE["K1","QF"], " ", Twiss["NX","$$$"]/2/Pi];\nExit;\n'
    k = _line(_run(src), 'K')
    assert float(k[0]) == pytest.approx(1.3827321232325072, abs=1e-8)
    assert float(k[1]) == pytest.approx(1.9, abs=1e-6)


def test_trailing_number_word_in_fit():
    r = _run(RING + 'FFS;\nUSE R;\nCELL;\nFIT M0 BX 3.;\nFIT M0 BY 3.\n;\nFIT M0 AX 3.0;\nExit;\n')
    assert 'Bad value' not in r


def test_use_resets_dp():
    # gs (tffsinitparam): USE のたびに、MARK の DP が正ならそれが DP になる。MARK に無ければ DP のまま (0 以下なら 0.01)
    src = ('MOMENTUM=3GEV;\nMARK M0=(DP=0.02) M1=() EMITX=1;\nDRIFT D1=(L=1);\nQUAD QF=(L=0.3 K1=0.5) QD=(L=0.3 K1=-0.5);\n'
           'LINE R0=(M0 QF D1 QD D1 QF QD);\nLINE R1=(M1 QF D1 QD D1 QF QD);\nFFS;\nUSE R0;\n'
           'DP = 0.005; USE R0; Print["a ", DP];\nDP = 0; USE R0; Print["b ", DP];\n'
           'DP = 0.005; USE R1; Print["c ", DP];\nDP = 0; USE R1; Print["d ", DP];\nExit;\n')
    r = _run(src)
    assert _line(r, 'a') == ['.02'] and _line(r, 'b') == ['.02']
    assert _line(r, 'c') == ['.005'] and _line(r, 'd') == ['.01']


ARCA_TOP = '''FFS;
USE ARCA;
LINE["K1","BPI7H"]  = -.006398424930727633;
LINE["K1","BPI6H"]  = -1.0100957067176413;
LINE["K1","QFA8H"]  =  .9417515241499581;
LINE["K1","QFA6H"]  =  .8782619605554324;
LINE["K1","QDA5"]   = -.4052202821210436;
LINE["K1","QFA4"]   =  .5166819921440575;
LINE["K1","QDA3"]   = -.49958031983962237;
LINE["K1","QDA2"]   = -.725123012102831;
LINE["K1","QFA1H"]  =  .46236106609239286;
LINE["K2","SFH"]    = 15.13439810682049;
LINE["K2","SDH"]    = -12.38248640793169;
nxcell = 48.58/20; nycell = 17.62/20; tay = 1.19; tby = 6.4; tex = 0.14; tbyc = 1.3317;
CELL; CALC;
FIT MKS BY 3.6;
FIT NX nxcell; FIT NY nycell;
FIT MKC AX 0; FIT MKC AY 0; FIT MKC EPX 0; FIT MKC BY tbyc;
FIT MKSFA1 AY tay; FIT MKSFA1 BY tby; FIT MKSFA1 EX tex; FIT MKSFA1 BX 3.;
FitValue["MKSFA1","BX",dp_,vgoal_,vnow_]:=If[vnow<5,5,If[vnow>15,15,Null]];
FitValue["MKS","BY",dp_,vgoal_,vnow_]:=If[vnow<0,0,If[vnow>3.6,3.6,Null]];
FitFunction:={(Abs[Twiss["BX","MKSFA1.1"]*Sin[Twiss["NX","MKSFA1.2"]-Twiss["NX","MKSFA1.1"]] - 0.3] + Abs[Sin[Twiss["NY","MKSFA1.2"]-Twiss["NY","MKSFA1.1"]]])*1};
FREE Q*;
FREE BPI6H K1;
FREE BPI7H K1;
GO;
Print["A ", LINE["K1","QFA1H"], " ", LINE["K1","QDA2"]];
Exit;
'''


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, 'example', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad')), reason='no example lattice')
def test_fitfunction_row_is_not_fit_in_k64(k64_mode):
    # gs の値 (QFA1H 0.5115137397134801、QDA2 -0.7640460317262181)。FitFunction の行を等式の条件にすると、反復が gs と別の谷に落ちる
    # (pysad v1.x は 0.5002)。k64 の扱いに固定する (conftest は k64)
    lat = open(os.path.join(ROOT, 'example', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad'), encoding='utf-8').read()
    r = _run(lat + chr(10) + ARCA_TOP)
    a = _line(r, 'A')
    assert float(a[0]) == pytest.approx(0.5115137397, abs=1e-6)
    assert float(a[1]) == pytest.approx(-0.7640460317, abs=1e-6)


PM1 = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
MULT MM = (L=0.1 K1=0.05 K2=2.0 DY=1e-4);
MULT MQ = (L=0.2 K1=0.4 K3=300 DY=-1e-4 DX=1e-4);
LINE CELL1 = (QF D2 B D2 MM QD D1 MQ);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW; RAD; NOFLUC;
POL;
em = Emittance[];
Print["OUT ",{PolarizationTime, SpinTune, EquilibriumPosition, PolarizationVector}/.em];
Exit;
'''


def test_emittance_starts_tcod_from_zero_when_qcod_fails(k64_mode):
    # K3 が強くずれの大きい MULT リング。CALC の閉軌道は gs でも見つからない (Overflow)。gs の qcod は見つからないとき初期値 0 を返すので、
    # Emittance[] の tcod は 0 から始まり小さい閉軌道 (1e-4) に収束する。gs の値 (k64): PolarizationTime 133.6984911450866、
    # SpinTune 0.40269221047113、EquilibriumPosition 0.06700392801591、PolarizationVector {4.465895419431E-5, 0.99999846587407, -0.00175107255363}。
    # (このリングは不安定 (Stable=0) で、DampingRate・TuneShiftByRadiation などは gs でも意味のない値なので比べない)
    r = _run(PM1)
    nums = [float(x) for x in re.findall(r'[-+]?\d*\.\d+(?:[Ee][-+]?\d+)?', _line_raw(r, 'OUT'))]
    assert nums[0] == pytest.approx(133.6984911450866, rel=1e-9)
    assert nums[1] == pytest.approx(0.40269221047113, rel=1e-9)
    assert nums[2] == pytest.approx(0.06700392801591, rel=1e-9)
    assert nums[3] == pytest.approx(4.465895419431e-5, rel=1e-7)
    assert nums[4] == pytest.approx(0.99999846587407, rel=1e-9)


def _line_raw(txt, tag):
    for l in txt.splitlines():
        if l.startswith(tag + ' '):
            return l
    raise AssertionError(tag + ' not printed:\n' + txt[-1500:])


SIMPLE = '''MOMENTUM=3GEV;
MARK M0=();
DRIFT D1=(L=1);
BEND B=(L=1 ANGLE=0.3926990817);
QUAD QF=(L=0.3 K1=0.5) QD=(L=0.3 K1=-0.5);
LINE C=(M0 QF D1 B D1 QD D1 B D1 QF);
FFS;
USE C;
'''


def test_division_is_a_times_inverse_in_k64(k64_mode):
    # gs: a/b は a*(1/b) (最後の桁が違う)。gs の値
    r = _run('FFS;\nPrint["g ", Log[100]/Log[10], " ", 7/3, " ", 0.3/0.1, " ", 100/3., " ", 5./7];\nExit;\n')
    assert _line(r, 'g') == ['1.9999999999999998', '2.333333333333333', '3', '33.33333333333333', '.7142857142857142']


def test_line_position_list_is_resolved_per_element():
    # gs: 位置の入れ子のリストは入れ子のまま、一致が 1 個ならスカラー・複数 (0 個も) ならリスト
    r = _run(SIMPLE + 'Print["a ", LINE["POSITION",{{1,2},{3,4}}]];\nPrint["b ", LINE["POSITION",{"QF","B*"}]];\n'
             'Print["c ", LINE["POSITION",{"B*"}]];\nPrint["d ", LINE["POSITION",{"Z*",1}]];\nPrint["e ", LINE["POSITION","B*"]];\nExit;\n')
    assert _line(r, 'a') == ['{{1,2},{3,4}}'] and _line(r, 'b') == ['{2,{4,8}}'] and _line(r, 'c') == ['{{4,8}}']
    assert _line(r, 'd') == ['{{},1}'] and _line(r, 'e') == ['{4,8}']


def test_too_many_indices_message():
    r = _run('FFS;\nPrint["a ", {f,g}[1,2]];\nPrint["b ", {{1,2},{3,4}}[1,2]];\nExit;\n')
    assert '???General::toomany:  Too many indices in {f,g}[1,2]' in r
    assert _line(r, 'b') == ['2']


def test_k1_readback_is_exact():
    r = _run(SIMPLE + 'Element["K1","QF"] = 0.93; Print["k ", Element["K1","QF"] - 0.93, " ", LINE["K1","QF"] - 0.93];\nExit;\n')
    assert _line(r, 'k') == ['0', '0']


def test_main_assignment_before_syntax_error():
    # gs: `A = 3 B = 4;` は A = 3 を代入する (B = 4 は構文エラー)
    r = _run('A = 3 B = 4;\nC = 5;\nDRIFT D1=(L=1);\nMARK M0=();\nLINE X=(M0 D1);\nFFS;\nUSE X;\nPrint["m ", A, " ", C];\nExit;\n')
    assert _line(r, 'm') == ['3', '5']


def test_line_emit_of_first_mark_uses_global_at_use():
    # gs: MARK に EMITY が無いとき LINE["EMITY",1]・Element["EMITY","M0"] は USE の時点のグローバル (USE のあとの代入は次の USE まで効かない)
    src = ('EMITY = 2e-11;\nEMITX = 1e-9;\nDRIFT D1=(L=1);\nMARK M0=(EMITX=5e-9);\nLINE X=(M0 D1);\nFFS;\nUSE X;\n'
           'Print["m2 ", LINE[{"EMITX","EMITY"},1], " ", EMITX, " ", LINE["EMITY","M0"], " ", Element["EMITY","M0"]];\n'
           'EMITY = 3e-11;\nPrint["m3 ", LINE["EMITY",1], " ", Element["EMITY","M0"]];\nUSE X;\nPrint["m4 ", LINE["EMITY",1], " ", EMITY];\nExit;\n')
    r = _run(src)
    assert _line(r, 'm2') == ['{5e-9,2e-11}', '5e-9', '2e-11', '2e-11']
    assert _line(r, 'm3') == ['2e-11', '2e-11']
    assert _line(r, 'm4') == ['3e-11', '3e-11']


COUPLE_FIT = '''MOMENTUM = 3GEV;
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
ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]};
CALC;
FREE QF;
FIT NX 7.0;
GO;
Print["R: ", LINE["K1","QF.1"], " ", LINE["K1","QD.1"], " ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
FIT NY 4.4;
GO;
Print["S: ", LINE["K1","QF.1"], " ", LINE["K1","QD.1"], " ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
Exit;
'''


def test_fit_with_element_values_coupling(k64_mode):
    # gs: ElementValues で QD を QF に結んだままフィットする (QF の変化に QD が追随する)。gs の値:
    # R: .7194073833 -.5755259067 7.0000013217 4.4701423728  S: .7163478092 -.5730782474 6.9661048595 4.444188379
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('PYSAD_SAD_PACKAGES が要る (Couple.n)')
    r = _run(COUPLE_FIT)
    rr = [float(x) for x in _line(r, 'R:')]
    ss = [float(x) for x in _line(r, 'S:')]
    for a, b in zip(rr, [.7194073833, -.5755259067, 7.0000013217, 4.4701423728]):
        assert a == pytest.approx(b, abs=2e-9)
    for a, b in zip(ss, [.7163478092, -.5730782474, 6.9661048595, 4.444188379]):
        assert a == pytest.approx(b, abs=2e-9)


def test_standalone_noexpand_is_silent():
    # gs: 単独の `NOEXPAND;`・`EXPAND;` は何も出さない (pysad は ?Undefined command を出していた。結果には影響なし)
    r = _run(SIMPLE + 'NOEXPAND;\nEXPAND;\nUSE C NOEXPAND;\nPrint["ok"];\nExit;\n')
    assert 'Undefined command' not in r and 'ok' in r


def test_couple_command_with_fit_matches_gs(k64_mode):
    # gs: COUPLE QD QF -0.8; FREE QF; FIT NX 7.0; GO → 2 段 (NX → NY を足す)。gs の値 R: .8158475149 -.6526780119 7.0000000426 4.4167559261、S: .8150246162 -.652019693 6.9919193435 4.4105963378
    src = COUPLE_FIT.replace('SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);\n', '').replace('LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);', 'LINE CELL = (DS QF D2 B D2 QD DS D1);')
    src = src.replace('ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]};\nCALC;\n', 'COUPLE QD QF -0.8;\n')
    r = _run(src)
    for tag, ref in (('R:', [.8158475149, -.6526780119, 7.0000000426, 4.4167559261]), ('S:', [.8150246162, -.652019693, 6.9919193435, 4.4105963378])):
        for a, b in zip([float(x) for x in _line(r, tag)], ref):
            assert a == pytest.approx(b, abs=2e-7)


def test_use_resets_dp0_to_zero(k64_mode):
    # gs: USE (文・FFS["USE .."]) のたびに MARK が作り直され DP0 (= LINE["DDP",1]) は 0。定義の DDP は無視される。CELL・CALC では戻らない。
    # gs の値: a .001, c 0, d 0, e 0, f 0, g 0, h .001 (DDP=0.004 の MARK でも USE のあとは 0)
    src = ('MOMENTUM=3GEV;\nMARK Ma=() Mb=(DDP=0.004);\nDRIFT D1=(L=1);\nQUAD QF=(L=0.3 K1=0.5) QD=(L=0.3 K1=-0.5);\n'
           'LINE La=(Ma QF D1 QD D1);\nLINE Lb=(Mb QF D1 QD D1);\nFFS;\nUSE La;\n'
           'DP0=0.001; Print["a ", DP0];\nDP0=0.001; USE La; Print["c ", DP0];\nDP0=0.001; FFS["USE Lb"]; Print["d ", DP0, " ", LINE["DDP",1]];\n'
           'DP0=0.001; CELL; Print["h ", DP0];\nExit;\n')
    r = _run(src)
    assert _line(r, 'a') == ['.001'] and _line(r, 'c') == ['0'] and _line(r, 'd') == ['0', '0'] and _line(r, 'h') == ['.001']
