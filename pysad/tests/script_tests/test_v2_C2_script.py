# -*- coding: utf-8 -*-
"""C2: 台本の読み方の不具合 2 件 (pysad 自身の不具合。互換モードによらない)。参照値は gs (K64、2026-10-19 ビルド) の出力。

(3) `FFS` の文の読み方: `FFS;NPARA=1;` (同じ行に続き)、`FFS USE RING;`、行の途中の `FFS;`、小文字。
(4) TRACK の無い台本でも MAIN 部の一般の代入 `NAME = 式;` が実行され、FFS に入ってからも見える。
"""
import io
import re

import pytest

LAT = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ST = ();
LINE RING = (ST QF D2 QD D1);'''
P = 'Print["R ", LINE["L","D1"], " ", Twiss["BX",1]];'
GS_BX = 3.6391435504154614            # gs: Twiss["BX",1] (CELL; CALC の後)


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _line(out, tag):
    rows = [l for l in out.splitlines() if l.startswith(tag + ' ')]
    assert rows, out[-600:]
    return rows[-1][len(tag) + 1:]


def _nums(s):
    return [float(v) for v in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', s)]


# ---------------------------------------------------------------- (3) FFS の文
@pytest.mark.parametrize('head, body', [
    ('FFS;NPARA=3;', 'USE RING; CELL; CALC;'),            # 同じ行の続きは FFS 部では実行されない (gs: NPARA は変わらない)
    ('FFS USE RING;', 'CELL; CALC;'),
    ('FFS USE=RING;', 'CELL; CALC;'),
    ('ffs;', 'USE RING; CELL; CALC;'),
    ('  FFS  ;  ', 'USE RING; CELL; CALC;'),
    ('FFS;   ! comment', 'USE RING; CELL; CALC;'),
    ('FFS', 'USE RING;\nCELL; CALC;'),                    # gs: FFS USE RING; と同じ
    ('FFS;Print["Z same line"]; NPARA=3;', 'USE RING; CELL; CALC;'),
])
def test_ffs_statement_forms(head, body):
    out = _run('%s\n%s\n%s\n%s\nExit;\n' % (LAT, head, body, P))
    v = _nums(_line(out, 'R'))
    assert v[0] == 1.0 and abs(v[1] / GS_BX - 1.0) < 1e-12
    assert not [l for l in out.splitlines() if l.startswith('Z same line')]


def test_ffs_in_the_middle_of_a_line():
    # gs: `LINE …; FFS; Print[…];` は FFS に入り、同じ行の Print は実行されず、次の行から FFS 部
    src = LAT.replace('LINE RING = (ST QF D2 QD D1);', 'LINE RING = (ST QF D2 QD D1); FFS; Print["Z dropped"];')
    out = _run('%s\nUSE RING; CELL; CALC;\n%s\nExit;\n' % (src, P))
    v = _nums(_line(out, 'R'))
    assert v[0] == 1.0 and abs(v[1] / GS_BX - 1.0) < 1e-12
    assert not [l for l in out.splitlines() if l.startswith('Z dropped')]
    assert 'FFS の文と同じ行の残り' in out                 # 捨てたことを知らせる


def test_rest_of_ffs_line_is_not_executed():
    # gs: `FFS USE RING; NPARA=3; CELL; CALC; Print[…]` の後ろは実行されない → 次の行の Print は CALC 前の値 (BX = 1)
    out = _run('%s\nFFS USE RING; NPARA=3; CELL; CALC; %s\n%s\nExit;\n' % (LAT, P.replace('"R "', '"Q "'), P))
    assert _nums(_line(out, 'R')) == [1.0, 1.0]
    assert not [l for l in out.splitlines() if l.startswith('Q ')]


def test_ffs_with_unknown_argument_is_reported():
    # gs: `FFS NPARA=3;` は「$FFS syntax error: undefined variable $NPARA」を出して FFS に入る
    out = _run('%s\nFFS NPARA=3;\nUSE RING; CELL; CALC;\n%s\nExit;\n' % (LAT, P))
    assert '$FFS' in out and 'NPARA' in out
    assert abs(_nums(_line(out, 'R'))[1] / GS_BX - 1.0) < 1e-12


def test_ffs_word_inside_other_text_is_not_the_command():
    from pysad.ffscmd import _locate_ffs, split_main_ffs, main_part
    assert _locate_ffs('A = 1;\nLINE X = (FFS1 FFS2);\n! FFS;\n') is None
    assert _locate_ffs('A = 1;\nB = FFS;\n') is None
    src = 'A = 1;  ! FFS; in a comment\nON RAD;\nFFS;\nCALC;\n'
    main, rest, flags = split_main_ffs(src)
    assert rest == 'CALC;\n' and flags == [('RAD', True)] and 'CALC' not in main and 'A = 1' in main
    assert main_part(src).rstrip().endswith('ON RAD;')
    # 前の文が ; で終わっていなくても、FFS だけの行は従来どおり FFS の文
    assert _locate_ffs('LINE X = (A B)\nFFS;\nCALC;\n')['main_end'] == 15


# ---------------------------------------------------------------- (4) MAIN 部の代入
SRC4 = '''MOMENTUM = 3 GEV;
A = 1.5;
B = A*2+0.25;
C = 10 DEG;
E1 = 2 GEV;
aa = 7;
F9 = SQRT(A+0.75)*2;
DRIFT D1 = (L=A) D2 = (L=B/5) D3 = (L=F9) D4 = (L=aa);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND BB = (L=1 ANGLE=C);
MARK ST = ();
A = 9;
LINE RING = (ST QF D2 QD D1 BB D3 D4);
FFS;
USE RING;
$FORM="20.14";
Print["R ", LINE["L","D1"], " ", LINE["L","D2"], " ", LINE["ANGLE","BB"], " ", LINE["L","D3"], " ", LINE["L","D4"]];
Print["V ", A, " ", B, " ", C, " ", E1, " ", F9];
Print["W ", aa, " ", AA, " ", a];
Print["M ", MOMENTUM];
A = 4;
Print["N ", A, " ", LINE["L","D1"]];
Exit;
'''


def test_main_assignments_reach_elements_and_ffs():
    out = _run(SRC4)
    # gs: R 1.5 .65 .17453292519943 (3 7)、V 9 3.25 .17453292519943 2000000000 3、W 7 AA a、M 3000000000、N 4 1.5
    r = _nums(_line(out, 'R'))
    for a, b in zip(r, [1.5, 0.65, 0.17453292519943, 3.0, 7.0]):      # 要素定義は定義の時点の値 (後の A = 9 は効かない)
        assert abs(a - b) < 1e-13
    v = _nums(_line(out, 'V'))
    for a, b in zip(v, [9.0, 3.25, 0.17453292519943, 2e9, 3.0]):      # FFS では最後の値
        assert abs(a - b) <= 1e-13 * max(1.0, abs(b))
    assert _line(out, 'W').split() == ['7.00000000000000', 'AA', 'a']  # 名前は大文字小文字を区別する (gs と同じ)
    assert _nums(_line(out, 'M')) == [3e9]
    assert _nums(_line(out, 'N')) == [4.0, 1.5]                        # FFS での代入は要素を変えない


SRC4E = '''MOMENTUM = 3 GEV;
A = 1.5;
D = 2^3;
F = 1E-3 M;
G = 5 MM;
P = PI;
X = 2*F;
Y = A/4-1;
T = -A*2+D;
U = 2*(A+1);
W1 = 1.5E0+A;
E2 = -A;
F1 = COS(0);
F2 = EXP(1);
F3 = LOG(100);
F4 = LN(10);
F5 = ATAN(1);
F6 = TAN(0.5);
F8 = sin(1);
F7 = COSH(0.5);
K2 = 2 deg;
Q1 = MOMENTUM/1E9;
Q2 = A^2^2;
Q3 = -A^2;
Q5 = 6/2/3;
Q6 = 1.5 M + 20 CM;
H1 = 2^3^2;
H3 = 8/2*4;
H4 = 2-3-4;
H5 = 2^(1+1)*3;
H6 = 2E0^.5;
H7 = EXP(-1);
H8 = 2*(-3);
E = -A;
MARK ST = ();
DRIFT D1 = (L=1);
LINE RING = (ST D1);
FFS;
USE RING;
$FORM="20.14";
Print["V ", {A, D, F, G, P, X, Y, T, U, W1, E2}];
Print["W ", {F1, F2, F3, F4, F5, F6, F8, F7, K2}];
Print["N ", {Q1, Q2, Q3, Q5, Q6, H1, H3, H4, H5, H6, H7, H8}];
Print["M ", E];
Exit;
'''
# gs の出力 ($FORM="20.14")
GS4E = {
    'V': [1.5, 8.0, 0.001, 0.005, 3.14159265358979, 0.002, -0.625, 5.0, 5.0, 3.0, -1.5],
    'W': [1.0, 2.71828182845905, 2.0, 2.30258509299405, 0.78539816339745, 0.54630248984379, 0.84147098480790,
          1.12762596520638, 0.03490658503989],
    'N': [3.0, 5.0625, -2.25, 1.0, 1.7, 512.0, 16.0, -5.0, 12.0, 1.41421356237310, 0.36787944117144, -6.0],
}


def test_main_expression_grammar_matches_gs():
    out = _run(SRC4E)
    flat = ' '.join(out.split())                           # 長いリストは行を折り返して表示される
    for tag, ref in GS4E.items():
        m = re.search(tag + r' (\{[^}]*\})', flat)
        got = _nums(m.group(1))
        assert len(got) == len(ref), (tag, got)
        for a, b in zip(got, ref):
            assert abs(a - b) <= 2e-14 * max(1.0, abs(b)), (tag, a, b)
    assert abs(_nums(_line(out, 'M'))[0] - 2.71828182845905) < 1e-13     # gs: MAIN の `E = -A` は FFS の定数 E を変えない


def test_main_expressions_gs_rejects_are_not_assigned():
    # gs でも構文エラーになる形 (gs は変数に無意味な値を入れることがある。pysad は代入しない)
    src = '''MOMENTUM = 3 GEV;
A = 1.5;
Q = UNDEF*2;
Z = 2 A;
G = Sqrt[16];
V1 = A**2;
Q4 = 2*-A;
H2 = 2^-1;
P2 = pi;
F4 = ABS(-2);
MARK ST = ();
DRIFT D1 = (L=1);
LINE RING = (ST D1);
FFS;
USE RING;
Print["V ", {A, Q, Z, G, V1, Q4, H2, P2, F4}];
Exit;
'''
    out = _run(src)
    assert _line(out, 'V').replace(' ', '') == '{1.5,Q,Z,G,V1,Q4,H2,P2,F4}'


def test_main_special_variables_reach_ffs():
    src = '''MOMENTUM = 3 GEV;
NPARA=3; PBUNCH=2E10; MINCOUP=0.02; DP=0.002;
EMITY=2E-11; EMITZ=3E-6; SIGE=1.1E-3; SIGZ=4E-3; LOSSAMPL=0.07; LOSSDZ=55; FSHIFT=1E-5;
TDXI=1E-4; TDYI=2E-4; TDZI=3E-4; DP0=0.003; MAXORDER=5; ADDTERMS=2; CHARGE=-1; MASS=1E8;
EMITDIV=7; EMITDIVB=8; EMITDIVQ=9; EMITDIVS=11; FRINGDIV=3; NBUNCH=6; PHSPLOTS=12;
EMITX=1E-9;
MARK ST = (EMITX=5E-9);
DRIFT D1 = (L=1);
LINE RING = (ST D1);
FFS;
USE RING;
Print["V ", {EMITY,EMITZ,SIGE,SIGZ,LOSSAMPL,LOSSDZ,FSHIFT}];
Print["W ", {TDXI,TDYI,TDZI,DP0,MAXORDER,ADDTERMS,CHARGE,MASS}];
Print["M ", {EMITDIV,EMITDIVB,EMITDIVQ,EMITDIVS,FRINGDIV,NBUNCH,PHSPLOTS}];
Print["N ", {PBUNCH, MINCOUP, EMITX}];
Print["Q ", DP];
Exit;
'''
    out = _run(src)
    # gs の出力そのまま。DP0・CHARGE は FFS の入口で既定値に戻る (0、1)。DP は 0.01 (単精度の .00999999977648 を gs は表示)。
    # MARK に書いた EMITX は MAIN の EMITX より優先、MARK に無いものは MAIN の値が残る
    assert _line(out, 'V') == '{2e-11,3e-6,.0011,.004,.07,55,1e-5}'
    assert _line(out, 'W') == '{.0001,.0002,.0003,0,5,2,1,1e8}'
    assert _line(out, 'M') == '{7,8,9,11,3,6,12}'
    assert _nums(_line(out, 'N')) == [2e10, 0.02, 5e-9]
    assert abs(_nums(_line(out, 'Q'))[0] - 0.01) < 1e-9


def test_main_assignments_with_track_statement():
    # MAIN に TRACK がある台本 (命令を順に実行する経路) でも一般の代入は FFS に届く
    src = LAT + '''
A = 1.5; B = A+1;
TRACK USE=RING NPART=1 TURNS=1 NX=(1) NY=(0) NZ=(0);
FFS;
Print["V ", A, " ", B];
Exit;
'''
    out = _run(src)
    assert _nums(_line(out, 'V')) == [1.5, 2.5]


def test_parse_main_assignments_units():
    from pysad.parser import parse_main_assignments, _eval_main_expr
    a = parse_main_assignments('A = 1.5; ! c\nDRIFT D = (L = A);\nf[x_] := (y = 2; x);\nB = A + 1 (* c *) ;\nMOMENTUM = 2 GEV; C = MOMENTUM/1E9')
    assert [(x['name'], x['kind'], x['value']) for x in a] == [('A', 'VAR', 1.5), ('B', 'VAR', 2.5), ('MOMENTUM', 'MOMENTUM', 2e9),
                                                               ('C', 'VAR', 2.0)]
    # gs (2026-11-18 に実測): 0 での割り算は Infinity、SQRT(-1) は NaN (エラーにならない)
    assert _eval_main_expr('1/0', {}) == float('inf')
    assert _eval_main_expr('SQRT(-1)', {}) != _eval_main_expr('SQRT(-1)', {})
