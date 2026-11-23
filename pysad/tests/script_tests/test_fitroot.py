"""FindRoot / Fit / FitGaussian / NIntegrate / FixedPointList / SelectCases。参照は SAD K64 (gs)。
K64 は記号微分、pysad は数値微分なので、Fit・FindRoot は相対 1e-6 程度、NIntegrate は 1e-12 で比べる。"""
import io

SRC = r'''FFS;
Print["NI ",NIntegrate[x^2,{x,0,1}]];
Print["NI2 ",NIntegrate[Sin[x],{x,0,Pi},AccuracyGoal->1e-8]];
r=Fit[{{0,1},{1,3},{2,7},{3,13.2}},a+b x+c x^2,x,{a,0},{b,0},{c,1}];
Print["Fa ",{a,b,c}/.r];
Print["Fc ",ChiSquare/.r];
Print["Fg ",GoodnessOfFit/.r];
Print["Fi ",ConfidenceInterval/.r];
r=Fit[{{0,1,.1},{1,3,.1},{2,7,.2},{3,13.2,.1}},a+b x+c x^2,x,{a,0},{b,0},{c,1}];
Print["Fb ",{a,b,c}/.r];
Print["Fd ",ChiSquare/.r];
r=FindRoot[x^2==2,{x,1}];
Print["FR ",x/.r];
r=FindRoot[{x+y==3,x y==2},{x,0.5},{y,2.5}];
Print["FS ",{x,y}/.r];
Print["FG ",FitGaussian[Table[Sin[1.7 k]+Cos[3.1 k^2]*.3,{k,100}]]];
r=Fit[{{0,1.05,.1,.1},{1,2.9,.1,.1},{2,5.2,.2,.1},{3,6.8,.1,.2},{4,9.3,.1,.1}},a+b*x,x,{a,0},{b,1}];
Print["FX ",{a,b}/.r," ",ChiSquare/.r];
Print["SC ",SelectCases[{1,2,3,4,5},{EvenQ,True&}]];
Exit;
'''


def _run():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    d = {}
    for l in buf.getvalue().splitlines():
        if l[:2] in ('NI', 'Fa', 'Fc', 'Fg', 'Fi', 'Fb', 'Fd', 'FR', 'FS', 'FG', 'SC', 'FX'):
            k, _, v = l.partition(' ')
            d[k] = v
    return d


def _nums(s):
    import re
    return [float(t) for t in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s)]


def test_fitroot_matches_gs():
    d = _run()
    ref = {'NI': [1 / 3], 'NI2': [1.9999999995870694],
           'Fa': [1.01, .91, 1.05], 'Fc': [.002], 'Fg': [.3173105078629139],
           'Fi': [.08185845541235924, .13145732335355942, .04199250159410424],
           'Fb': [1.004255319148937, .9425531914893619, 1.0404255319148934], 'Fd': [.08510638297872336],
           'FR': [2 ** .5], 'FS': [1.0, 2.0], 'FX': [.9489011872, 2.0510404949, 2.0388543524],
           'FG': [.8275559516159802, .012067673076785912, .011240680298617984, .010540347287583861, 403.26814032862774]}
    for k, rv in ref.items():
        if k == 'NI2':
            continue
        got = _nums(d[k])
        assert len(got) == len(rv), (k, d[k])
        tol = 1e-12 if k == 'NI' else 1e-6
        for a, b in zip(got, rv):
            assert abs(a - b) <= tol * max(abs(b), 1e-3), (k, a, b)
    assert abs(_nums(d['NI2'])[0] - 2.0) < 1e-8


def test_selectcases():
    assert _run()['SC'] == '{{2,4},{1,3,5}}'
