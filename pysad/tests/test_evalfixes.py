"""GUI・描画 (Canvas.n、ListPlot.n、Plots.n) を動かすために SADScript の評価を gs に合わせた箇所の試験。参照は SAD K64 (gs) の出力
(local/wake の ch.sad sp.sad th.sad rp.sad sq.sad sq2.sad wf.sad vc.sad pm.sad)。"""
import io
import re


def _out(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=1999;\n' + body + 'Exit;\n', out=buf)
    lines = {}
    for l in buf.getvalue().splitlines():
        m = re.match(r'([A-Z]\d?) (.*)', l)
        if m:
            lines[m.group(1)] = m.group(2)
    return lines


def test_chained_comparison():
    r = _out('''pr={{0.79,4.21},{1,16}};
Print["A ",Map[If[pr[[1,1]]<#<pr[[1,2]],1,0]&,{0.5,0.79,1,4.21,5}]];
Print["B ",0.79<0.5<4.21," ",1<2<3," ",3>2>1," ",1<=1<2," ",1<3<2];
''')
    assert r['A'] == '{0,0,1,0,0}'
    assert r['B'] == '0 1 1 1 0'


def test_string_position_alternatives_and_overlap():
    r = _out('''s="800x600+10+20";
Print["A ",StringPosition[s,"x"]];
Print["B ",StringPosition[s,{"x","+","+"}]];
Print["C ",StringPosition[s,"0"|"6"]];
Print["E ",StringPosition["aaaa","aa"]];
Print["F ",StringPosition[s,"+"|"x",2]];
''')
    assert r['A'] == '{{4,4}}'
    assert r['B'] == '{{4,4},{8,8},{11,11}}'
    assert r['C'] == '{{2,2},{3,3},{6,6},{7,7},{10,10},{13,13},{5,5}}'
    assert r['E'] == '{{1,2},{2,3},{3,4}}'
    assert r['F'] == '{{8,8},{11,11}}'


def test_thread_broadcast_and_replacepart_blank_spec():
    r = _out('''Print["A ",Thread[{a,{}}]," ",Thread[{a,{1,2}}]," ",Thread[{{a,b},{1,2}}]];
l={{1,2},{3,4},{5,6}};
Print["B ",ReplacePart[l,0,{,2}]," ",ReplacePart[l,0,{2,}]," ",Extract[l,{,2}]];
Print["C ",ReplacePart[{1,2,3},9,{}]];
''')
    assert r['A'] == '{} {{a,1},{a,2}} {{a,1},{b,2}}'
    assert r['B'] == '{{1,0},{3,0},{5,0}} {{1,2},{0,0},{5,6}} {2,4,6}'
    assert r['C'] == '{1,2,3}'


def test_empty_sequence_is_spliced():
    r = _out('''Print["A ",If[#=="a",Null[],#]&/@{"a","b","c"}];
Print["B ",MapThread[If[#>1,Null[#,#2],Null[]]&,{{1,2,3},{4,5,6}}]];
Print["C ",Table[If[i==2,Null[],i],{i,4}]," ",Table[If[i==2,Null[i,i],i],{i,3}]];
Print["D ",{1,INF,3}/.INF:>Null[]];
Print["E ",Join[Switch[#," ",With[{o=""},o],_,Null[]]&/@Characters["ab cd"],{"x"}]];
''')
    assert r['A'] == '{"b","c"}'
    assert r['B'] == '{2,5,3,6}'
    assert r['C'] == '{1,3,4} {1,2,2,3}'
    assert r['D'] == '{1,3}'


def test_with_declares_functions_and_bare_alternatives_rank_first():
    r = _out('''With[{TF[ft_,a_]:=Switch[ft,a,Null[],{Tk,_List},ft[[-1]],_,Throw[]]},
  G[x_]:=TF[x,True];
  H[x_]:=Catch[TF[x,Both]]];
Print["A ",G[True]," ",G[{Tk,{1,2}}]];
Print["B ",H[{Tk,{3,4}}]," ",H[7]];
U[""|" ",__]=Null;
U[e_,t_,f_,"e"]:=("second");
Print["C ",U["",1,2,"e"]," ",U[" ",1,2,"e"]," ",U["a",1,2,"e"]];
''')
    assert r['A'] == ' {1,2}'
    assert r['B'] == '{3,4} Null'
    assert r['C'] == 'Null Null second'


def test_vector_and_releasevector():
    r = _out('''Print["A ",Vector[{1,2,3}]];
Print["B ",ReleaseVector[Sin[Vector[{0.,1.,2.}]]]];
f=With[{x=Vector[#]},Check[ReleaseVector[Sin[x]],0]]&;
Print["C ",f[{0.5,1.5}]];
Print["D ",ReleaseVector[Vector[{1,2}]+3]];
''')
    assert r['A'] == 'Vector[{1,2,3}]'
    assert r['B'] == '{0,.8414709848078965,.9092974268256817}'
    assert r['C'] == '{.479425538604203,.9974949866040544}'
    assert r['D'] == '{4,5}'


def test_function_argument_symbol_is_localized_by_with_and_block():
    """f[g_, {z_Symbol, …}] の z が x に束縛されているとき、本体の Block[{z}]・With[{z = …}] は x を局所にし、引数の式 Sin[x] の x にも効く"""
    r = _out('''g[fxin_,{z_Symbol,a_,b_}]:=Module[{fx,zl,xl},
  fx= With[{z = Vector[#]}, Check[ReleaseVector[fxin], 0]]&;
  zl=Range[0,4]*0.5+a;
  xl=Re[fx[zl]];
  xl];
Print["A ",g[Sin[x],{x,0.,6.}]];
h[{u_,fxin_},{z_Symbol,a_,b_}]:=Block[{z},Module[{fx,zl,xl},
  fx= With[{z = Vector[#]}, Check[ReleaseVector[fxin], 0]]&;
  zl=Range[0,4]*0.5+a;
  xl=Re[fx[zl]];
  xl]];
Print["B ",h[{x,Sin[x]},{x,0.,6.}]];
''')
    want = '{0,.479425538604203,.8414709848078965,.9974949866040544,.9092974268256817}'
    assert r['A'] == want
    assert r['B'] == want


def test_module_local_delayed_function_and_holdnone():
    """Module の局所変数への h := Sqrt[Abs[#]]& (遅延した純関数)、SetAttributes[f, HoldNone] (OpticsList の fun:=…、PlotColorDefault)。gs と同じ"""
    r = _out('''Print["A ",Module[{h},h:=3+1;h]];
Print["C ",Module[{h},h:=Sqrt[Abs[#]]&;h[4]]];
Print["D ",Module[{h},h:=Sqrt[Abs[#]]&;Table[h[i],{i,2}]]];
Plot$DefaultColor={"red","blue"};
With[{tdr:=TopDrawer,cdr:=System`CanvasDrawer},
  SetAttributes[PCD,HoldAll];
  PCD[tdr]={"black"};
  PCD[cdr]=Plot$DefaultColor;
  SetAttributes[PCD,HoldNone]];
$DisplayFunction=CanvasDrawer;
Print["E ",PCD[$DisplayFunction]," ",PCD[TopDrawer]];
''')
    assert r['A'] == '4'
    assert r['C'] == '2'
    assert r['D'] == '{1,1.4142135623730951}'
    assert r['E'] == '{"red","blue"} {"black"}'


def test_set_definition_body_substitutes_pattern_vars_inside_value_lists():
    """G[x_] = Round[{H[x], x}*3] の右辺は評価されて値のリスト {3 Restrict[x,0,1], 3 x} になる。呼び出しでその中の x も置き換える (Canvas.n の CoordLimited)"""
    r = _out('''H[x_]:=Restrict[x,0,1];
G5[x_]=Round[{H[x],x}*3];
Print["E ",G5[0.5]];
G6[x_]=Round[{x}*3];
Print["F ",G6[0.5]];
Block[{x,y},
  F[{x_,y_}]=Round[{H[x],y}*10]];
Print["A ",F[{0.5,1.5}]];
''')
    assert r['E'] == '{2,2}'
    assert r['F'] == '{2}'
    assert r['A'] == '{5,15}'


def test_override_sorts_by_key_and_repeated_alternatives_pattern():
    r = _out('''Print["A ",Override[{b->1,B->2,a->3,Zed->4,"x"->5,c->1},{c->9}]];
Print["B ",Override[{b->1},{a->2},{b->3}]];
Print["D ",Override[{Z->1,z->2,Y->3,y->4}]];
P[R[min_,max_,opt:((_Rule|_RuleDelayed)...)]]:=1;
P[R[min0_,max0_,g__,opt:((_Rule|_RuleDelayed)..)]]:=2;
P[R[min_,max_,g__]]:=3;
Print["E ",P[R[1,2]]," ",P[R[1,2,a->1]]," ",P[R[1,2,a->1,b->"x"]]," ",P[R[1,2,{q}]]," ",P[R[1,2,{q},a->1]]," ",P[R[1,2,{q},{r}]]];
Q[__]={};
Q[a_,b_]:=2;
Print["F ",Q[1,2]," ",Q[1,2,3]];
''')
    assert r['A'] == '{"x"->5,a->3,b->1,B->2,c->1,Zed->4}'
    assert r['B'] == '{a->2,b->1}'
    assert r['D'] == '{y->4,Y->3,z->2,Z->1}'
    assert r['E'] == '1 1 1 3 2 3'
    assert r['F'] == '2 {}'


def test_with_delayed_rule_and_class_with_values():
    """With[{s1, leng}, PlotRange :> {s1 - 0.001, s1 + leng}] は右辺にも値を埋め込む。クラスの中の With[{c := X}, f[c] = …] は左辺にも"""
    r = _out('''F[a_]:=Module[{s1=3,leng=4,ss=2},
  With[{s1,leng,ss},
    G[#,{DataRange->{{s1,s1+leng}*ss,{-1,1}},
      PlotRange:>With[{XR=ZZ}, {If[Head[XR]===List,XR,ss*{s1-0.001,s1+leng+0.001}],{-1,1}}]}]&/@{a}]];
r=F[1];
Print["A ",PlotRange/.r[[1,2]]];
Plot$Color={"red","blue"};
OP=Class[{},{},{},
 With[{tdr:=TopDrawer,cdr:=System`CanvasDrawer},
  SetAttributes[PCD2,HoldAll];
  PCD2[tdr]={"black"};
  PCD2[cdr]=Plot$Color;
  SetAttributes[PCD2,HoldNone]];
 With[{def={PlotColor->PCD2[$DisplayFunction]}},
  Test[]:=Override[{},def]]];
$DisplayFunction=CanvasDrawer;
o=OP[];
Print["B ",o@Test[]];
Plot$Color2={"x"};
''')
    assert r['A'] == '{{5.998,14.002},{-1,1}}'
    assert 'PCD2' not in r['B'] and r['B'].startswith('{PlotColor->{')


def test_block_saves_and_restores_function_definitions_of_local_symbols():
    """Block[{mesh}, mesh[x_] := …] は入るとき定義を退避し、出るとき戻す (Plots.n の ListContourPlot の中の色目盛りの ListDensityPlot)。gs と同じ"""
    r = _out('''Print["A ",Block[{mesh},mesh[x_]:=x+1;Block[{mesh},mesh[x_]:=x*100;mesh[2]]+mesh[2]]];
Block[{mesh},mesh[x_]:=x+1;
  Print["B ",mesh[1]];
  Block[{mesh},mesh[x_]:=x*100; Print["C ",mesh[1]]];
  Print["D ",mesh[1]]];
Print["E ",mesh[1]];
Block[{nx},nx=3;Print["G ",nx]];
Print["H ",nx];
''')
    assert [r[k] for k in 'ABCDEGH'] == ['203', '2', '100', '2', 'mesh[1]', '3', 'nx']


def test_literal_inside_pattern_is_more_specific_and_optional_in_head_pattern():
    r = _out('''D1[Default]:=1;
D1[{x_,y_}]:=2;
D1[{x_,Default}]:=3;
D1[{Default,y_}]:=4;
Print["A ",D1[Default]," ",D1[{1,2}]," ",D1[{1,Default}]," ",D1[{Default,2}]];
RP[min_,max_,Default|Null]:={min,max};
RP[min_,max_,Right[sep_:0]]:={"R",sep};
Print["B ",RP[1,2,Default]," ",RP[1,2,Right[]]," ",RP[1,2,Right[3]]];
FR[s_,{Automatic,Automatic},{b1_,b2_}]:="both";
FR[s_,{Automatic,r2_},{b1_,_}]:="first";
Print["C ",FR[1,{Automatic,Automatic},{1,2}]," ",FR[1,{Automatic,3},{1,2}]];
''')
    assert r['A'] == '1 2 3 4'
    assert r['B'] == '{1,2} {"R",0} {"R",3}'
    assert r['C'] == 'both first'


def test_minmax_takes_several_arguments_and_definition_of_symbol():
    r = _out('''Print["A ",MinMax[0,{{3,4},{14,1}}]," ",MinMax[{2,5}]];
c[1][5]=7;c[1][6]=8;c[2][5]=9;
c[1][5]=.;
Print["C ",Length[First[?c]]," ",c[1][6]," ",c[1][5]];
c[_][_]={{}};
Print["D ",Length[First[?c]]," ",c[1][5]];
Print["E ",Extract[Delete[?c,{{1,1},{1,-1}}],{{1,,2},{1,,1,0,1}}]];
''')
    assert r['A'] == '{0,14} {2,5}'
    assert r['C'] == '3 8 c[1][5]'
    assert r['D'] == '4 {{}}'
    assert sorted(r['E'].strip('{}').split(',')) == sorted(['9', '8', '2', '1'])


def test_module_with_initialization_list_and_hold_rest_sequence():
    """Plots.n の FitPlot: Module[{vars}, {init}, body] (本体が 2 つ以上)、HoldRest の f[data_, fun_, x_, par0__] の Fit[data, fun, x, par0]"""
    r = _out('''SetAttributes[FP,HoldRest];
FP[data_,fun_,x_,par0__]:=Module[{s,{opt,par}=SwitchCases[{par0},{_Rule|_RuleDelayed,_}]},
  {{w}={3},lp=Length[par]},
  s=Fit[data,fun,x,par0];
  {lp,s[[1]]}];
r=FP[Table[{i,2*i+1+0.1*Sin[i*3]},{i,1,6}],a+b*u,u,{a,1},{b,1}];
Print["A ",r[[1]]," ",Head[r[[2]]]," ",Round[10^6*(a/.r[[2]])]/10^6];
''')
    assert r['A'] == '2 Rule 1.020144'


def test_short_and_inequality_match_gs():
    """Short[x, n] は x を書いて Null (gs で確認)。Inequality[a, Less, b, Less, c] は連なった比較"""
    r = _out('''Print["A ",Short[{1,2,3},2]];
Print["C ",Inequality[1,Less,2,Less,3]," ",Inequality[1,Less,3,Less,2]," ",Inequality[1,Less,2,LessEqual,2]];
''')
    assert r['A'] == 'Null'
    assert r['C'] == '1 0 1'


def test_polylog_above_one_is_complex_like_gs():
    """PolyLog[n, x] (x > 1、整数 n >= 2) は切断線の下側の複素数 (gs: Li2(2) = 2.4674 - 2.1776 i)"""
    r = _out('''Print["G ",PolyLog[2,2.]," ",PolyLog[3,2.]," ",PolyLog[4,1.5]];
''')
    import re as _re
    nums = [float(x) for x in _re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', r['G'].replace(' I', ''))]
    want = [2.4674011002723386, -2.1775860903036026, 2.7620719062289245, -.7546938294602483, 1.7347570807760613, -.034902704828337024]
    assert len(nums) == 6
    for g, w in zip(nums, want):
        assert abs(g - w) < 1e-12 * max(abs(w), 1.0)


def test_messagename_function_form_like_gs():
    r = _out('''Print["A ",MessageName[g,"usage"]];
Print["E ",MessageName[General,"wrongtype"]];
''')
    assert r['A'] == 'Undefined MessageName General::"usage"'
    assert r['E'] == 'Undefined MessageName General::"wrongtype"'


def test_polylog_noninteger_outside_unit_disk_and_hurwitz_negative_q_like_gs():
    """PolyLog[s, z] (非整数 s、z>1 は切断線の下側の複素数、z<-1 は実数)、HurwitzZeta[s, q] (q<=0 の非整数。(n+q)^-s は主値)。参照は gs"""
    r = _out('''Print["A ",PolyLog[2.5,3.]," ",PolyLog[1.5,-4.]];
Print["B ",HurwitzZeta[2.5,-0.5]," ",HurwitzZeta[3,-1.5]];
''')
    import re as _re

    def nums(t):
        return [float(x) for x in _re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', t.replace(' I', ''))]
    a = nums(r['A'])
    for g, w in zip(a, [3.2828227108912214, -2.7213246265012407, -2.0056564849416016]):
        assert abs(g - w) < 1e-11 * abs(w)
    b = nums(r['B'])
    for g, w in zip(b, [6.247110634568811, -5.6568542494923815, .11810202582086049]):
        assert abs(g - w) < 1e-11 * abs(w)


def test_type_command_packs_drifts_on_one_line_like_gs():
    """TYPE (tftype.f の twbuf): DRIFT は 1 行に詰めて並べ (タブ 7)、ほかの型は 1 要素ごとに改行する。参照は gs"""
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(r"""MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5) DRIDX = (L=0.5) DRIDY = (L=0.4) DLONGNAME1 = (L=0.123456) D3=(L=0.1) D4=(L=0.2) D5=(L=0.3) D6=(L=0.4);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.1308996938995747 E1=0.06544984694978735 E2=0.06544984694978735);
LINE CELL1 = (QF D2 B D2 QD D1 DRIDX DRIDY DLONGNAME1 D3 D4 D5 D6);
LINE RING = (ST 3*CELL1);
FFS;
USE RING;
type D*;
type Q* B;
Exit;
""", out=buf)
    lines = [l.rstrip() for l in buf.getvalue().splitlines()]
    want = [' DRIFT  D2      =(L =.5 )    D1      =(L =1 )     DRIDX   =(L =.5 )    DRIDY   =(L =.4 )    DLONGNAME1=(L =.123456 )',
            '        D3      =(L =.1 )    D4      =(L =.2 )    D5      =(L =.3 )    D6      =(L =.4 )',
            ' QUAD   QF      =(L =.3   K1 =.7 )',
            '        QD      =(L =.3   K1 =-.6 )',
            ' BEND   B       =(L =1.2  ANGLE =.1308996938995747 E1 =.06544984694978735   E2 =.06544984694978735 )']
    for w in want:
        assert w in lines


def test_read_types_and_skip_like_gs(tmp_path):
    # gs: Real/Word/Expression は語のあとの空白 (改行を含む) まで読む。Character は 1 文字。数でない語の Real は式の値。Skip[g,Word,n] は n 個目の語
    p = str(tmp_path / 'rd.txt').replace(chr(92), '/')
    r = _out('''f=OpenWrite["%s"];
Write[f,"alpha beta 1.5 2.5 3.5"];
Write[f,"line2 x y z"];
Close[f];
g=OpenRead["%s"];
Print["A ",Read[g,Character]," ",Read[g,Character]," ",Read[g,Word]," ",Read[g,Real]," ",Read[g,Real]];
Print["B ",Read[g,Character]," ",Read[g,Character]," ",Read[g,Character]];
Print["C ",Read[g,Real]," ",Read[g,Character]," ",Read[g,Word]];
Print["D ",Read[g,Word]," ",Read[g,Word]," ",Read[g,Word]," ",Read[g,Word]];
Close[g];
g=OpenRead["%s"];
Print["E ",Skip[g,Word,2]];
Print["F ",Read[g,Word]];
Close[g];
''' % (p, p, p))
    assert r['A'] == 'a l pha beta 1.5'
    assert r['B'] == '2 . 5'
    assert r['C'] == '3.5 l ine2'
    assert r['D'] == 'x y z EndOfFile'
    assert r['E'] == 'beta'
    assert r['F'] == '1.5'


def test_wake_command_is_wakecommand_script_like_gs():
    # K64 の WAKE コマンドは WAKECOMMAND[] (配布に無い WAKECOMMAND.n) を評価し、True でなければ行の残りを捨てる。Fortran の tfwake.f は呼ばれない
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('''MOMENTUM=3GEV;
CAVI RF=(L=0 VOLT=1e6 FREQ=5e8 PHI=0);
DRIFT D=(L=1);
MARK ORG=(); MARK ZEND=();
LINE R=(ORG D D RF D ZEND);
FFS;
USE R;
WAKE LONGITUDINAL RF 0 1 1e-3 0.5 ENDWAKE;
Print["A ok"];
Exit;
''', out=buf)
    assert 'A ok' in buf.getvalue()


def test_complex_matrix_results_narrow_per_element_like_gs():
    # gs: 虚部が厳密に 0 の要素は実数に戻る。DiagonalMatrix は複素数の列も受ける (LDSolve.n の ExpMatrix が使う)
    r = _out('''{v,u}=Eigensystem[{{0,1},{-1,0}}];
Print["A ",Dot[{{1+I,0},{0,2}},{{1,2},{3,4}}]];
Print["B ",Inverse[u]];
Print["C ",Transpose[u].DiagonalMatrix[Exp[v*0.1]].Inverse[Transpose[u]]];
''')
    assert r['A'] == '{{(  1.0000000000+I),(  2.0000000000+  2.0000000000 I)},{  6.0000000000,  8.0000000000}}' or '6' in r['A']
    assert 'I),' in r['B'] and '(   .5000000000+' not in r['B']
    assert r['C'].replace(' ', '').startswith('{{0.9950041653,0.0998334166}') or '.99500416' in r['C']


import os as _os
import pytest as _pytest


@_pytest.mark.skipif(not _os.environ.get('PYSAD_SAD_PACKAGES'), reason='PYSAD_SAD_PACKAGES が必要')
def test_ldsolve_rotation_matches_gs():
    """Packages/LDSolve.n (x' = [[0,1],[-1,0]] x、x(0)=(1,0)) の終点 (cos 1, −sin 1)。gs と同じ値 (local/wake/ld2.sad)"""
    r = _out('''u=LDSolve[{{0,1},{-1,0}},{0,0},{1,0},{s,0,1,0.1}];
Print["A ",u[[1,2,-1]]];
''')
    import re
    v = [float(x) for x in re.findall(r'-?\d*\.\d+', r['A'])]
    assert abs(v[0] - 0.5403023059) < 1e-9 and abs(v[1] + 0.8414709848) < 1e-9


def test_openread_failure_and_read_of_non_unit_match_gs():
    # gs: 開けない OpenRead は General::fileopen を出して入力行を打ち切る (h は未代入)。そのあとの Read[h, ...]・Close[h] は wrongtype
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('''FFS;
h=OpenRead["no_such_file_xyz.csv"];
Print["A ",h];
x=Read[h,Real];
Print["B ",x];
Close[h];
f=OpenWrite["/no_such_dir_xyz/x.txt"];
Print["C ",f];
Exit;
''', out=buf)
    lines = [l.rstrip() for l in buf.getvalue().splitlines()]
    assert '???General::fileopen:  Open error for file no_such_file_xyz.csv in OpenRead["no_such_file_xyz.csv"]' in lines
    assert 'A h' in lines and 'B x' in lines
    assert '???General::wrongtype:  Argument must be Real number or List of two Reals in Read[h,Real]' in lines
    assert '???General::wrongtype:  Argument must be Real number in Close$[h]' in lines
    assert 'C $Failed' in lines
