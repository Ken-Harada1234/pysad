"""記号微分 D (Packages/D.n を pysad が SADScript のまま実行) と、そのためのパターン照合 (Orderless の Plus・Times、Flat の積、f/@(a+b)、
f[Unevaluated$[g[y,x]]]、ReplaceRepeated など)。参照は gs の出力 (空白を詰めた文字列。symbolic_ref.json)。SAD の Packages が見つからなければ skip。"""
import io
import json
import os
import re

import pytest


def _norm(s):
    return re.sub(r'\s+', '', s)


def test_symbolic_derivative_matches_gs():
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('SAD の Packages がない')
    from pysad.ffscmd import run_sad_string
    ref = json.load(open(os.path.join(os.path.dirname(__file__), 'symbolic_ref.json'), encoding='utf-8'))
    src = 'FFS;\n$FORM="18.12";\n' + ''.join('Print["R%d= ",%s];\n' % (i, e) for i, (e, _) in enumerate(ref)) + 'Exit;\n'
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    got = {}
    for l in buf.getvalue().splitlines():
        m = re.match(r'R(\d+)= (.*)', l)
        if m:
            got[int(m.group(1))] = _norm(m.group(2))
    bad = [(e, g, got.get(i)) for i, (e, g) in enumerate(ref) if got.get(i) != g]
    assert not bad, bad[:5]


def test_orderless_matching_and_replace_repeated():
    from pysad.ffscmd import run_sad_string
    src = r'''FFS;
f1[x_*y_]:={x,y};
f9[b___*Sin[x_]*Cos[x_]]:={b,x};
f10[Sin[x_]^2+Cos[x_]^2]:=one[x];
Print["A ",f1[2 x]," ",f1[x*2]," ",f9[w*Cos[z]*Sin[z]]," ",f10[Cos[t]^2+Sin[t]^2]];
Print["B ",(f[f[f[a]]])//.f[z_]:>z," ",f/@(a+b+c)];
Print["C ",(2 x)^2," ",-x//Length];
Exit;
'''
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    lines = {l[0]: l[2:] for l in buf.getvalue().splitlines() if l[:2] in ('A ', 'B ', 'C ')}
    assert lines['A'] == '{2,x} {2,x} {w,z} one[t]'
    assert lines['B'] == 'a f[a]+f[b]+f[c]'
    assert lines['C'].startswith('(4 x^2)')


def test_spline_matches_gs():
    """Spline・SplineData (D.n) と Spline$ (自然スプライン y''/6)。参照は gs"""
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('SAD の Packages がない')
    from pysad.ffscmd import run_sad_string
    src = r'''FFS;
$FORM="14.10";
s=Spline[{{0,0},{1,1},{2,4},{3,9},{4,16}}];
Print["1 ",s[1.5]," ",s[0.2]," ",s[3.9]," ",s[{0.5,2.5}]];
Print["2 ",Head[s]," ",Derivative[1][s][1.5]," ",Derivative[2][s][1.5]];
Print["3 ",NIntegrate[s[x],{x,0,4}]," ",Integrate[s[x],{x,0.5,3.2}]];
Print["5 ",Spline[{{0,{0,1}},{1,{1,2}},{2,{4,3}}}][0.5]];
Exit;
'''
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    lines = {l[0]: _norm(l[2:]) for l in buf.getvalue().splitlines() if l[:2] in ('1 ', '2 ', '3 ', '5 ')}
    assert lines['1'] == '2.2321428571.117714285715.2575714286{.3392857143,6.2321428571}'
    assert lines['2'] == 'SplineData3.03571428572.1428571429'
    assert lines['3'] == '21.428571428610.8837035714'
    assert lines['5'] == '{.3125000000,1.5000000000}'
