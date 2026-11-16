"""example/10_math_data〜13_packages と PF の T5 の gs との差を直した箇所の試験。参照は SAD K64 の gs の出力 (2026-11-17)。
SAD の Packages が見つからなければ、パッケージを使う試験は skip。"""
import io
import os
import re

import pytest

_HERE = os.path.dirname(__file__)


def _need_pkg():
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('SAD の Packages がない')


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _tags(out):
    d = {}
    for l in out.splitlines():
        m = re.match(r'([A-Za-z]\w*) (.*)', l)
        if m:
            d[m.group(1)] = m.group(2)
    return d


def test_makebump_upvalue_set_loaded_lazily():
    """b=Bump[…] は bump.n の `(a_=Bump[…])^:=…` (Set の UpValue)。最初の使用が Bump でも (MakeBump より前でも) 効く。gs: ybump=.0010000141"""
    _need_pkg()
    src = open(os.path.join(_HERE, 'mkbump_ring.sad'), encoding='utf-8').read()
    out = _run(src)
    assert 'list assignment' not in out
    assert 'steerers used=12' in out
    v = {}
    for l in out.splitlines():
        m = re.match(r'R:(\w+)=(\S+)', l)
        if m:
            v[m.group(1)] = float(m.group(2))
    assert abs(v['ybump'] - 1.0000141e-3) < 1e-9
    assert abs(v['ybefore'] - (-7.80747637435e-6)) < 1e-10
    assert abs(v['yafter'] - 9.706571350986e-6) < 1e-10


_RING4 = '''MARK ORIGIN=();
DRIFT D=(L=1);
QUAD QF=(L=.3 K1=.7) QD=(L=.3 K1=-.6);
LINE R=(ORIGIN QF D QD D QF D QD D);
FFS;
USE R; CELL; CALC;
'''


def test_string_head_expression_stays_symbolic_and_with_hold():
    """gs: "K1"["QD"] は評価しない式 (Head は K1)、With[{y},Hold[y]] は Hold[値]、(Hold[f[#]]&)[3] は Hold[f[#]]"""
    out = _run('FFS;\nPrint["A ",Head["K1"["QD"]]];\nPrint["B ",{"K1"["QD"]->-0.8*"K1"["QF"]}];\n'
               'y={1,2};\nPrint["C ",With[{y},Hold[y]]];\nPrint["D ",With[{y=5},Hold[y]]];\nExit;\n')
    d = _tags(out)
    assert d['A'] == 'K1'
    assert d['B'] == '{"K1"["QD"]->(-.8 "K1"["QF"])}'
    assert d['C'] == 'Hold[{1,2}]'
    assert d['D'] == 'Hold[5]'


def test_element_values_applied_at_calc_like_gs():
    """ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]} は次の CALC で QD の K1 = -0.56 にする (Couple.n)。USE で消える。
    gs: b -.56、USE 後は ElementValues が未定義に戻り QD の K1 は -.6"""
    _need_pkg()
    out = _run(_RING4 + '''ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]};
CALC;
Print["B ",LINE["K1","QD"]];
Print["N ",Twiss["NX","$$$"]/(2Pi)];
USE R;
Print["C ",ElementValues];
CALC;
Print["D ",LINE["K1","QD"]];
Exit;
''')
    d = _tags(out)
    assert abs(float(d['B'].split(',')[0].strip('{}')) - (-0.56)) < 1e-12
    assert d['C'] == 'ElementValues'
    assert abs(float(d['D'].strip('{}')) - (-0.6)) < 1e-12


def test_rational_arithmetic_uses_package_upvalues_like_gs():
    """rational.n の `Plus[Rational[…],…]^:=…` などの UpValue が組み込みの Plus・Times・Minus・Divide・Subtract に効く (gs と同じ)。
    ReducedRational は p^2>b で打ち切るので、1/3+1/6 は Rational[3,6] (gs も同じ)"""
    _need_pkg()
    out = _run('''FFS;
a=Rational[1,3];b=Rational[2,3];c=Rational[1,6];
Print["A ",a+b];
Print["B ",a+c+b];
Print["C ",a*b];
Print["D ",-a];
Print["E ",1/a];
Print["F ",0.5+a];
Print["H ",a-c];
Print["K ",a/b];
Print["N ",Rational[1,3]+Rational[1,6]];
Print["L ",1+2];
Exit;
''')
    d = _tags(out)
    assert d == {'A': '1', 'B': 'Rational[7,6]', 'C': 'Rational[2,9]', 'D': 'Rational[-1,3]', 'E': '3', 'F': 'Rational[2.5,3]',
                 'H': 'Rational[1,6]', 'K': 'Rational[3,6]', 'N': 'Rational[3,6]', 'L': '3'}


def test_nintegrate_complex_integrand_like_gs():
    """gs: NIntegrate[Exp[I*x],{x,0,Pi}] = (8.6e-16+1.999999999999996 I)。実数の場合は従来どおり"""
    out = _run('FFS;\nPrint["A ",NIntegrate[Exp[I*x],{x,0,Pi}]];\nPrint["B ",NIntegrate[Exp[I*x],{x,0,Pi},Spline->False]];\n'
               'Print["C ",NIntegrate[Sin[x],{x,0,Pi}]];\nExit;\n')
    d = _tags(out)
    m = re.match(r'\(([-+.\deE]+)\+([-+.\deE]+) I\)', d['A'])
    assert abs(float(m.group(1))) < 1e-13 and abs(float(m.group(2)) - 1.999999999999996) < 1e-12
    m = re.match(r'\(([-+.\deE]+)\+([-+.\deE]+) I\)', d['B'])
    assert abs(float(m.group(2)) - 2.0) < 1e-3
    assert abs(float(d['C']) - 1.999999999999996) < 1e-12


def test_loggamma1_is_loggamma_of_x_plus_1_like_gs():
    """gs: LogGamma1[10.5]=16.292000476567242、[0.5]=-.12078223763524522、[-1.5]=(1.2655121234846454-3.141592653589793 I)"""
    out = _run('FFS;\nPrint["A ",LogGamma1[10.5]];\nPrint["B ",LogGamma1[0.5]];\nPrint["C ",LogGamma1[{1.5,2.5}]];\n'
               'Print["D ",LogGamma1[-1.5]];\nPrint["E ",LogGamma1[x]];\nExit;\n')
    d = _tags(out)
    assert abs(float(d['A']) - 16.292000476567242) < 1e-12
    assert abs(float(d['B']) - (-0.12078223763524522)) < 1e-12
    v = [float(t) for t in d['C'].strip('{}').split(',')]
    assert abs(v[0] - 0.2846828704729192) < 1e-12 and abs(v[1] - 1.2009736023470743) < 1e-12
    assert d['D'].startswith('(1.2655121234846') and d['D'].endswith('-3.141592653589793 I)')
    assert d['E'] == 'LogGamma1[x]'


def test_parametricplot_pi_in_range_and_curve_list_like_gs():
    """gs: {t,0,2*Pi} (Pi を含む範囲) は x_Real に照合され Graphics、曲線のリストは List (曲線数 + 1 要素)。
    Evaluate[…] は保持引数の位置で評価される"""
    _need_pkg()
    out = _run('''FFS;
$DisplayFunction=Identity;
g1=ParametricPlot[{Sin[3*t],Cos[2*t]},{t,0,2*Pi}];
Print["A ",Head[g1]," ",Length[g1]];
g2=ParametricPlot[{{2*Cos[t],2*Sin[t]},{Cos[t],Sin[t]}},{t,0,6.28}];
Print["B ",Head[g2]," ",Length[g2]];
g3=ParametricPlot[{{2*Cos[t],2*Sin[t]},{Cos[t],Sin[t]}},{t,0,2*Pi},Legend->"Expressions"];
Print["C ",Head[g3]," ",Length[g3]];
a=Hold[{Cos[t],Sin[t]}];
g4=ParametricPlot[Evaluate[ReleaseHold[a]],{t,0,Pi}];
Print["D ",Head[g4]," ",Length[g4]];
Exit;
''')
    d = _tags(out)
    assert d == {'A': 'Graphics 2', 'B': 'List 3', 'C': 'List 3', 'D': 'Graphics 2'}


def test_held_typed_blank_matches_constant_expression_like_gs():
    """gs: HoldAll の f[x_Real] に 2*Pi・Sin[1.]・1+2. は定数式として照合する。変数 y (値あり) は照合しない"""
    out = _run('''FFS;
SetAttributes[f,HoldAll];
f[x_Real]:=1;
f[x_]:=2;
y=3.;
Print["A ",f[2*Pi]];
Print["B ",f[Sin[1.]]];
Print["C ",f[1+2.]];
Print["D ",f[y]];
Print["E ",f[Pi]];
Exit;
''')
    d = _tags(out)
    assert d == {'A': '1', 'B': '1', 'C': '1', 'D': '2', 'E': '1'}


def test_pattern_test_error_is_mismatch_and_listcontourplot_sadhelp_data():
    """gs: _?(Second[#]>3&) の検査の中のエラー (Second[List]) は不一致 (メッセージなし)。Position の level 1 は頭部 List も検査する。
    SADHelp の ListContourPlot の例 (MeshRange つき) が Graphics を返す"""
    _need_pkg()
    out = _run('''FFS;
$DisplayFunction=Identity;
Print["A ",Position[{{1,2},{3,4}},_?(Second[#]>3&),1,1]];
Print["B ",Cases[{{1,2},{3,4}},_?(Second[#]>3&),1]];
xr0=0;xr=Table[xr0+=Sin[i*Pi/21],{i,20}];xr-=(xr[[1]]+xr[[-1]])/2;
yr=xr;
table=Outer[Cos[Sqrt[(#^2+#2^2)/2]]&,xr,yr];
g=ListContourPlot[table,MeshRange->{xr,yr},AspectRatio->1,Contours->8,FrameLabel->{"x","y"}];
Print["C ",Head[g]," ",Length[g]];
Print["D ",DataRange/.g[[2]]];
Exit;
''')
    assert '???' not in out
    d = _tags(out)
    assert d['A'] == '{{2}}' and d['B'] == '{{3,4}}' and d['C'] == 'Graphics 2'
    assert d['D'].startswith('{{-6.5975151867')


def test_optimizeoptics_output0_prints_no_iterations_like_gs():
    """OptimizeOptics[Output->0, MaxIteration->60] は Join[{opt},{VariableRange->limit}] (規則のリスト) を DownhillSimplex に 1 引数で渡す。
    gs は Output->0 で反復の行を出さない。MaxIteration も効く"""
    _need_pkg()
    src = open(os.path.join(_HERE, 'oo1_optimize.sad'), encoding='utf-8').read()
    out = _run(src)
    assert 'DownhillSimplex' not in out
    assert not re.search(r'^\s+\d+\.\d{6}\s+\d+\.\d{6}', out, re.M)
    d = _tags(out)
    assert d['K'].split() == ['.7059715383', '-.5989199800']
    # Output->2 は従来どおり途中経過を出す
    out2 = _run(src.replace('Output->0', 'Output->6'))
    assert 'DownhillSimplex' in out2


def test_mapthread_broadcasts_non_list_elements_like_gs():
    """gs: MapThread[List,{{1,2},3}] = {{1,3},{2,3}}、False は 0。リストが 1 つも無ければエラー"""
    out = _run('''FFS;
Print["A ",MapThread[List,{{1,2},3}]];
Print["B ",MapThread[List,{{1,2},False}]];
Print["C ",MapThread[List,{{1,2},{3,4},"x"}]];
Print["D ",MapThread[#+#2&,{{1,2},{10,20}}]];
Print["E ",MapThread[List,{{1,2},{3}}]];
Exit;
''')
    d = _tags(out)
    assert d['A'] == '{{1,3},{2,3}}' and d['B'] == '{{1,0},{2,0}}' and d['C'] == '{{1,3,"x"},{2,4,"x"}}' and d['D'] == '{11,22}'
    assert 'equalleng' in out


def test_structured_pattern_test_rule_is_tried_before_general_rule():
    """OpticsPlot.n: OpticsList1$[pr_,{f_}?(Length[Dimensions[#]]>=3&),opt___] は OpticsList1$[{pr1_,pr2_},f_,opt___] より先 (gs は定義の順)"""
    out = _run('''FFS;
l1={{1,2},{3,4},{5,6}};
g[pr_,{f_}?(Length[Dimensions[#]]>=3&),opt___]:=(Print["U"];g[pr,f,opt]);
g[{p1_,p2_},f_,opt___]:=Print["G ",Length[f]];
g[{1,2},{l1},PlotColor->"red"];
Exit;
''')
    assert re.search(r'^U\nG 3$', out, re.M)


def test_opticsplot_with_own_data_lists_like_gs():
    """gs: OpticsPlot[{{{line1},{line2},{line3},PlotJoined->True,Dashing->{1,1,1,1,1},PlotColor->{…},Legend->False}},Lattice->True,Names->"B*"] が
    Graphics を返す (データの x は要素番号。Legend->False の MapThread、{line}?(Dimensions…) の規則の順)"""
    _need_pkg()
    out = _run('''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (DS QF D2 B D2 QD DS D1);
LINE RING = (ORIGIN RF 24*CELL);
FFS;
USE RING; CELL; CALC;
$DisplayFunction=Identity;
L=LINE["LENGTH"];
pos=Range[1,L];
line1=Thread[{pos,Table[Twiss["BX",i],{i,1,L}]}];
line2=Thread[{pos,Table[Twiss["BY",i],{i,1,L}]}];
line3=Thread[{pos,Table[Twiss["EX",i],{i,1,L}]}];
p2=OpticsPlot[{{{line1},{line2},{line3},PlotJoined->True,Dashing->{1,1,1,1,1},PlotColor->{"black","red","blue","forest green"},FrameLabel->"Dx(mm)",TextSize->1.8,Legend->False}},Lattice->True,Names->"B*"];
Print["P ",Head[p2]," ",Length[p2]];
Exit;
''')
    assert '???' not in out
    assert _tags(out)['P'] == 'Graphics 2'


def test_stringjoin_does_not_concatenate_lists_like_gs():
    """gs: StringJoin[{"a","b"}] と {"a","b"}//"" は {"a","b"} の形の文字列 (リストを連結しない)。PointSymbol の名前の文字の元"""
    out = _run('''FFS;
sym={{"CI","CL","1O"},{Bar,"BA","BR"}};
Print["A ",(sym[[1]]//"")[2,-2]];
Print["B ",(sym[[2]]//"")[2,-2]];
Print["C ",StringJoin[{"a","b"}]];
Print["D ",StringJoin["x",{"b",1}]];
Print["E ",StringJoin["a","b",3]];
Print["F ","ab"//"cd"];
Print["G ",{1.5,"b"}//"-"];
Exit;
''')
    d = _tags(out)
    assert d['A'] == '"CI","CL","1O"' and d['B'] == 'Bar,"BA","BR"'
    assert d['C'] == '{"a","b"}' and d['D'] == 'x{"b",1}' and d['E'] == 'ab3' and d['F'] == 'abcd' and d['G'] == '{1.5,"b"}-'
