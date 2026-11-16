"""DownhillSimplex (optimize.n)、StringToStream など。参照は SAD K64 (gs)。"""
import io
import re

SRC = r'''FFS;
f=((#[[1]]-1)^2+(#[[2]]+2)^2)&;
p0=Map[{f[#],#}&,{{0,0},{1,0},{0,1}}];
Print["DS ",DownhillSimplex[Sort[p0],f,Output->0]];
g=(100*(#[[2]]-#[[1]]^2)^2+(1-#[[1]])^2)&;
p1=Map[{g[#],#}&,{{-1.2,1},{-1,1},{-1.2,.8}}];
r=DownhillSimplex[Sort[p1],g,Output->0,Tolerance->1e-12,MaxIteration->500];
Print["RB ",r[[1]]];
r=DownhillSimplex[Sort[p1],g,Output->0,ReflectRange->{2,-1},VariableRange->{{-2,-2},{2,.5}},MaxIteration->3];
Print["RE ",r];
f=StringToStream["abc def 12 34\nline2"];
Print["R ",Read[f,Word]," ",Read[f,Word]," ",Read[f,Real]];
Print["LR ",ListRandom[]];
Exit;
'''


def test_misc_matches_gs():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    d = {}
    for l in buf.getvalue().splitlines():
        k, _, v = l.partition(' ')
        d[k] = v
    assert d['DS'] == '{{0,{1,-2}},{0,{1,-2}},{0,{1,-2}}}'
    nums = [float(t) for t in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', d['RB'])]
    ref = [3.8499258148115327, -.961728395066932, .9209876543306752]
    for a, b in zip(nums, ref):
        assert abs(a - b) < 1e-9
    assert d['RE'] == '{{4,{-1,1}},{4.3569140625,{-1.025,1}},{4.672539062499999,{-1.025,.975}}}'
    assert d['R'] == 'abc def 12'
    assert d['LR'] == '{"SAD"}'


def test_fromdate_without_list_is_current_time_and_mapat_on_rule():
    """FromDate[] / FromDate[Tick->dt] は現在時刻 (Tkinter.n・KBMainFrame.n が使う)。MapAt は a->b の右辺 {2} も置き換える (gs と同じ)"""
    import io
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPrint["A ",FromDate[]>3e9," ",FromDate[Tick->0.01]>3e9," ",FromDate[{2020,1,2,3,4,5}]];\n'
                   'Print["B ",MapAt[f,a->b,{2}]];\nExit;\n', out=buf)
    out = buf.getvalue()
    assert 'A  1.000000' in out or 'A 1' in out
    assert '3786923045' in out
    assert 'a->f[b]' in out
