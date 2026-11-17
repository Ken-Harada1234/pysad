# -*- coding: utf-8 -*-
"""v2 の続き (担当 B2): SADScript・FFS の小項目。参照は gs (SAD K64、2026-10-19 ビルド) の出力 (2026-11-18)。

(4) `$Context`・`$ContextPath` への代入、`WithContext`、入れ子の文脈、`Begin["A`"]` (文字列) の拒否。
(5) ElementValues など 5 つの記号の VISIT 退避・BYE 復元・USE の消去。
(6) `Symbol["x"]` の値、`Message[f::tag]` の表示。
(7) `Ffs;`・`FFs;` (FFS と ffs 以外の綴り) の拒否、`USE ring` と `USE RING` の区別、`ffs use ring;`。
(3) `Parallelize` (POSIX は fork、Windows は逐次)。
(8) Packages (functions.n・Misc.n・element.n・beamline.n・emit.n) の関数の gs との差の修正。
"""
import io
import os
import re
import sys

import pytest

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _run(body, cwd=None, main=''):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(main + 'FFS;\nPageWidth=1999;\n' + body + '\nExit;\n', out=buf, cwd=cwd)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _line(lines, tag):
    for l in lines:
        if l.startswith(tag):
            return l[len(tag):]
    raise AssertionError('no line %r in %r' % (tag, lines))


def _nums(s):
    return [float(t) for t in re.findall(_NUM, s)]


_LAT = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ST = ();
LINE RING = (ST QF D2 QD D1);
LINE R2 = (ST QD D2 QF);
'''


# ---------------------------------------------------------------- (6) Symbol・Message
def test_symbol_returns_value_and_message_prints(k64_mode):
    r = _run('Print["A ",Symbol["zq"]];zq=3;Print["B ",Symbol["zq"]+1];Symbol["zq"]=9;Print["C ",zq," ",Head[Symbol["zq"]]];'
             'f::x="a `1` b `2`";Message[f::x,10,"q"];Message[f::x];Message[f::none];g::usage="usage of g";Message[g::usage];'
             'Print["D ",Message[g::usage]];Message["abc"];Print["E ",Symbol["Sin"][0.5]]')
    assert _line(r, 'A ') == 'zq'                                               # gs: 未定義なら記号のまま、値があれば値
    assert _line(r, 'B ') == '4'
    assert _line(r, 'C ') == '9 Real'
    assert '???f::x:  a 10 b q' in r                                            # gs の表示
    assert '???f::x:  a `1` b `2`' in r
    assert '???f::none:  Undefined MessageName General::none' in r
    assert r.count('???g::usage:  usage of g') == 2
    assert _line(r, 'D ') == 'Null'
    assert abs(_nums(_line(r, 'E '))[0] - 0.479425538604203) < 1e-15


# ---------------------------------------------------------------- (7) Ffs・USE
def _ffs_main(head, tail='USE RING; CELL; CALC;'):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_LAT + head + '\n' + tail + '\nPrint["R ", LINE["L","D1"], " ", Twiss["BX",1]];\nExit;\n', out=buf)
    return buf.getvalue()


@pytest.mark.parametrize('head', ['Ffs;', 'FFs;'])
def test_ffs_misspelled_is_a_syntax_error(head):
    out = _ffs_main(head)
    assert 'LALR(1) parser: syntax error' in out                                # gs: 読めない。FFS 部に入らず何も実行されない
    assert 'R ' not in out


def test_ffs_lowercase_and_use_names_are_case_sensitive():
    out = _ffs_main('ffs;')
    assert abs(float(out.split('R 1 ')[1]) / 3.6391435504154614 - 1.0) < 1e-12   # ffs は読める
    out = _ffs_main('ffs use ring;', 'CELL; CALC;')                              # gs: ring という LINE は無い (LINE 名は書かれたとおり)
    assert '???-FFS-Error-Missing beamline in USE.' in out
    out = _ffs_main('FFS;', 'USE ring; CELL; CALC;')
    assert '???-FFS-Error-Missing beamline in USE.' in out
    out = _ffs_main('FFS;', 'USE RING; CELL; CALC;')
    assert 'Missing beamline' not in out


# ---------------------------------------------------------------- (5) VISIT / BYE
def test_visit_bye_save_and_restore_blocked_symbols(k64_mode):
    r = _run('USE RING; CELL; CALC;\nElementValues={K1[QF]:>0.5}; FitFunction:=1+2; FitValue=7; OpticsProlog=11; OpticsEpilog=12;\n'
             'Print["A ",ElementValues," ",FitFunction," ",FitValue," ",OpticsProlog," ",OpticsEpilog];\nVISIT R2;\n'
             'Print["B ",ElementValues," ",FitFunction," ",FitValue," ",OpticsProlog," ",OpticsEpilog];\n'
             'ElementValues={K1[QD]:>0.1}; FitValue=99;\nBYE;\n'
             'Print["D ",ElementValues," ",FitFunction," ",FitValue," ",OpticsProlog," ",OpticsEpilog];\n'
             'VISIT R2; ElementValues={K1[QD]:>0.1}; VISIT RING; Print["E ",ElementValues]; BYE; Print["F ",ElementValues]; BYE; Print["G ",ElementValues];\n'
             'USE R2; Print["H ",ElementValues," ",FitFunction," ",FitValue," ",OpticsProlog," ",OpticsEpilog];\nBYE;',
             main=_LAT.replace('FFS;', ''))
    assert _line(r, 'A ') == '{K1[QF]:>.5} 3 7 11 12'
    assert _line(r, 'B ') == 'ElementValues FitFunction FitValue OpticsProlog OpticsEpilog'      # gs: 訪問先では 5 つとも未定義
    assert _line(r, 'D ') == '{K1[QF]:>.5} 3 7 11 12'                                           # BYE で戻る
    assert _line(r, 'E ') == 'ElementValues'
    assert _line(r, 'F ') == '{K1[QD]:>.1}'
    assert _line(r, 'G ') == '{K1[QF]:>.5}'
    assert _line(r, 'H ') == 'ElementValues FitFunction FitValue OpticsProlog OpticsEpilog'      # USE は 5 つとも消す
    assert any('BYE without VISIT' in l for l in r)


# ---------------------------------------------------------------- (4) 文脈
_CTX_BODY = '''Begin[A`]; Print["1 ", $Context]; x = 1; Print["2 ", Context[x], " ", x];
Begin[A`B`]; Print["3 ", $Context]; w6 = 6; Print["4 ", Context[w6], " ", x]; End[]; Print["5 ", $Context]; End[];
Print["6 ", $Context, " ", x, " ", A`B`w6];
Begin[`B`]; Print["7 ", $Context]; End[];
$Context = A`; Print["8 ", $Context]; $Context = `; Print["9 ", $Context];
$Context = "A`B`"; Print["10 ", $Context]; w7 = 6; Print["11 ", Context[w7]]; $Context = `;
$ContextPath = {A`, System`}; Print["12 ", $ContextPath, " ", x];
$ContextPath = {`, System`};
WithContext[{Zed`}, zz = 4; Print["13 ", $Context, " ", $ContextPath]];
Print["14 ", $Context, " ", $ContextPath];
WithContext[{Zed`, A`}, Print["15 ", $Context, " ", $ContextPath]];
Print["16 ", Symbol["A`x"], " ", Symbol["x"]];'''


def test_context_assignments_nested_and_withcontext_match_gs(k64_mode):
    r = _run(_CTX_BODY)
    # gs の出力 (Get ではなく同じ文を FFS に与えた実測)
    assert _line(r, '1 ') == 'A`'
    assert _line(r, '2 ') == '` 1'                          # x は init.n が大域に作っているので、文脈 A` に新しい記号は作られない
    assert _line(r, '3 ') == 'A`B`'                         # 入れ子の文脈
    assert _line(r, '4 ') == 'A`B` 1'                       # 大域に無い名前 (w6) は現在の文脈に作られる
    assert _line(r, '5 ') == 'A`'
    assert _line(r, '6 ') == '` 1 6'                         # A`B`w6 は文脈 A`B` に作った変数 (gs の c2 実測は大域の y を使い A`B`y が未評価だった)
    assert _line(r, '7 ') == 'B`'
    assert _line(r, '8 ') == 'A`' and _line(r, '9 ') == '`'
    assert _line(r, '10 ') == 'A`B`'                        # $Context = "A`B`" (文字列) は $SetContext を呼ぶ
    assert _line(r, '11 ') == 'A`B`'                        # 読み取りの文脈が変わる
    assert _line(r, '12 ') == '{`A`,`System`} 1'            # gs は ` が $ContextPath に無いと、文脈の記号を `A` と表示する
    assert _line(r, '13 ') == '`Zed` {`Zed`,`System`}'
    assert _line(r, '14 ') == '` {`,System`}'
    assert _line(r, '15 ') == '`Zed` {`Zed`,`A`,`System`}'
    assert _line(r, '16 ') == 'A`x 1'


def test_begin_with_string_is_rejected_in_k64_and_accepted_in_v2():
    from pysad import compat
    old = compat._mode
    try:
        compat.set_mode('k64')
        r = _run('Begin["A`"]; Print["1 ", $Context]; End[]; Begin["`B`"]; Print["2 ", $Context]; Print["3 ", Context[Sin]];')
        assert '???General::wrongtype:  Argument must be Context (Symbol ending `) in $SetContext["A`"]' in r      # gs
        assert _line(r, '1 ') == '`' and _line(r, '2 ') == '`'
        compat.set_mode('v2')
        r = _run('Begin["A`"]; Print["1 ", $Context]; End[]; Print["2 ", Context["A`x"]];')
        assert _line(r, '1 ') == 'A`'                                  # v2: Mathematica と同じく文字列も受ける
        assert _line(r, '2 ') == 'A`'
    finally:
        compat.set_mode(old)


def test_set_context_with_invalid_value_keeps_the_context(k64_mode):
    r = _run('$Context = "abc"; Print["1 ", $Context]; $ContextPath = {1, 2}; Print["2 ", $ContextPath];')
    assert any('Argument must be Context (Symbol ending `) in $SetContext["abc"]' in l for l in r)
    assert _line(r, '1 ') == '`'
    assert any('List of Contexts' in l for l in r)
    assert _line(r, '2 ') == '{`,System`}'


# ---------------------------------------------------------------- (3) Parallelize
_PAR = ('Print["A ",Parallelize[(#^2)&,5,200,2]];Print["B ",Parallelize[Function[k,Table[k*i,{i,3}]],4,500,4]];'
        'Print["C ",Parallelize[Function[k,If[k==2,1/0,k]],3,200,3]];x=0;Print["E ",Parallelize[Function[k,x=k;k],3,100,3]," ",x];'
        'Print["D ",Parallelize[Function[k,Table[k,{i,100}]],2,100,2]];Print["S ",Parallelize[Function[k,k+1],3,100,1]];')


def _check_par(r):
    assert _line(r, 'A ') == '{1,4,9,16,25}'                       # gs と同じ
    assert _line(r, 'B ') == '{{1,2,3},{2,4,6},{3,6,9},{4,8,12}}'
    assert _line(r, 'C ') == '{1,Infinity,3}'
    assert _line(r, 'S ') == '{2,3,4}'


@pytest.mark.skipif(not hasattr(os, 'fork'), reason='fork の無い環境 (Windows) では逐次')
def test_parallelize_forks_on_posix():
    r = _run(_PAR)
    _check_par(r)
    assert _line(r, 'E ') == '{1,2,3} 0'                           # gs と同じく、子の代入は親に戻らない
    assert _line(r, 'D ') == '{$Failed,$Failed}'                   # lshare を超える結果は $Failed (子が Shared::toolarge を出す)


@pytest.mark.skipif(hasattr(os, 'fork'), reason='Windows (fork が無い環境) の動作')
def test_parallelize_is_serial_on_windows():
    r = _run(_PAR)
    _check_par(r)
    assert _line(r, 'E ') == '{1,2,3} 3'                           # 逐次なので代入が残る (gs と違う。docs/history/v2next_B2.md)


def test_parallelize_serial_by_environment(monkeypatch):
    monkeypatch.setenv('PYSAD_PARALLELIZE', 'serial')
    r = _run('Print["A ",Parallelize[(#^2)&,5,200,2]];')
    assert _line(r, 'A ') == '{1,4,9,16,25}'


# ---------------------------------------------------------------- (8) Packages
def test_functions_n_misc_n_against_gs(k64_mode):
    r = _run('Print["a1 ",Log2$[10,100]," ",Log2$[2,8]];'
             'Print["a5 ",Union[{3,1,2,1,3}]," ",Union[{1,1.0000001,2},SameTest->(Abs[#-#2]<1e-3&)]," ",Union[{5,1,3,1,5},SameTest:>(#==#2&)]];'
             'sq[x_Real]:=x^2; sq[x_String]:=x//x;'
             'Print["c7 ",FunctionDefinition[Subtract[1,2]]];Print["c8 ",FunctionDefinition[sq[3.]]];Print["c9 ",FunctionDefinition[sq[{1}]]];'
             'zq1=1; f1[x_]:=x; Print["n1 ",Names["Sin"]," ",Names["Para*"]," ",Names["zq*"]," ",Names["f1"]];'
             'Print["r ",First[a->1]," ",Last[a->1]," ",Length[a->1]," ",Map[First,{a->1,b->2}]," ",First[a:>b]];'
             'Print["s ",Fill$String["abc",9],"|",Remove$Paren["(abc)"]," ",Remove$Paren["abc"]];')
    assert _nums(_line(r, 'a1 ')) == [1.9999999999999998, 2.9999999999999996]   # gs: K64 の a/b は a*(1/b) (最後の桁)
    a5 = _line(r, 'a5 ')
    assert a5 == '{1,2,3} {1} {5}'                                              # gs: SameTest をつけると結合した先頭の要素だけが残る (K64 の誤り)
    assert _line(r, 'c7 ') == 'Null'
    assert _line(r, 'c8 ') == '{Hold[sq[x_Real]],Hold[x^2],Hold[3^2]}'
    assert _line(r, 'c9 ') == 'Null'                                            # sq[{1}] に合う定義は無い
    assert _line(r, 'n1 ') == '{} {} {"zq1"} {"f1"}'                            # gs: 大域に作った記号だけ
    assert _line(r, 'r ') == 'a 1 2 {a,b} a'
    assert _line(r, 's ') == 'abc      |abc abc'


def test_union_samplest_true_value_in_v2(v2_mode):
    r = _run('Print["a ",Union[{3,1,2,1,3},SameTest->(#===#2&)]," ",Union[{1,2,3,4,5,6},SameTest->(Abs[#-#2]<1.5&)]];')
    assert _line(r, 'a ') == '{1,2,3} {1,3,5}'


def test_element_op_assign_goes_through_element_and_setelement_returns_main_values(k64_mode):
    r = _run('USE RING; CELL; CALC;\nElement["K1","QF"]=0.9; Element["K1","QF"]+=0.1; Print["a ",Element["K1","QF"]," ",LINE["K1","QF"]];'
             'Element["VALUE","QF"]=0.5; Print["b ",Element["K1","QF"]];Element["K1","QD"]*=2; Print["c ",Element["K1","QD"]," ",LINE["K1","QD"]];'
             'LINE["K1","QF"]+=0.25; Print["d ",LINE["K1","QF"]];'
             'Print["e ",SetElement["QF"][[3,2]]];Print["f ",SetElement["QF","QUAD",{"L"->0.35}][[3,{1,2}]]," ",Element["L","QF"]];', main=_LAT.replace('FFS;', ''))
    assert abs(_nums(_line(r, 'a '))[0] - 1.0) < 1e-12 and abs(_nums(_line(r, 'a '))[1] - 1.0) < 1e-12
    assert abs(_nums(_line(r, 'b '))[0] - 0.5) < 1e-12                       # 以前は += のあと Element の代入が効かなくなっていた
    assert abs(_nums(_line(r, 'c '))[0] + 1.2) < 1e-12 and abs(_nums(_line(r, 'c '))[1] + 1.2) < 1e-12
    assert abs(_nums(_line(r, 'd '))[0] - 0.75) < 1e-12
    assert _line(r, 'e ') == '"K1"->.7'                                     # gs: MAIN の定義の値 (FFS で変えた値ではない)
    assert _line(r, 'f ').startswith('{"L"->.35,"K1"->.7} .35')


def test_writebeamline_name_string_k64_vs_v2():
    from pysad import compat
    old = compat._mode
    body = 'PrintBeamLine[BeamLine[QF,D1],Format->"MAIN",Name->"L2"];'
    try:
        compat.set_mode('k64')
        r = _run('USE RING; CELL; CALC;\n' + body, main=_LAT.replace('FFS;', ''))
        assert 'WriteBeamLine - Insufficient number of Names' in r         # gs (K64 の誤り: Name は {"L2"} と書く)
        compat.set_mode('v2')
        r = _run('USE RING; CELL; CALC;\n' + body, main=_LAT.replace('FFS;', ''))
        assert any(l.startswith('LINE L2') for l in r)
    finally:
        compat.set_mode(old)


def test_emittance_false_omits_radiation_rules(k64_mode):
    main = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ORIGIN = ();
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (QF D2 B D2 QD DS D1);
LINE RING = (ORIGIN RF 8*CELL);
'''
    r = _run('USE RING; CELL; CALC; RFSW; RAD; NOFLUC;\nr5=Emittance[Emittance->False]; Print["e0 ",Map[First,r5]," ",Stable/.r5];'
             'r=Emittance[]; Print["e1 ",Length[r]];', main=main)
    assert _line(r, 'e0 ') == ('{Stable,Region,Tunes,EnergyLossU0,RfVoltageVc,EquilibriumPosition,MomentumCompaction,OrbitDilation,BucketHeight,'
                               'HarmonicNumber,OrbitAtExit} 1')                    # gs と同じ並び
    assert _line(r, 'e1 ') == '22'
