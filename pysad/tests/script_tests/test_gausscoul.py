"""GaussianCoulomb / GaussianCoulombU (tspac.f bbkick1、tintrb.f twspu)。参照は SAD K64 (gs) の出力 (12 桁、15 点)。"""
import io
import re

SRC = r'''FFS;
$FORM="18.12";
l = {{0.3,0.2,1,0.5},{1,2,1,0.5},{0,0,1,0.5},{0.5,0,1,1},{-1,2,1,0.5},{0.7,-0.4,1,1},{3,4,2,0.1},{10,3,1,0.2},{0.001,0.002,1,0.3},{2,0,0.5,1},{0.3,0.4,1,0.99999},{6,5,1,0.5},{20,10,1,0.5},{0.0,2,1,0.5},{-3,-1,1,0.7}};
Do[Print["GC ", GaussianCoulomb@@l[[i]]], {i, Length[l]}];
Do[Print["GU ", GaussianCoulombU@@l[[i]]], {i, Length[l]}];
Exit;
'''
REF = ['GC (    -.380080864473    -.499796970274 I)', 'GC (    -.304429918079    -.759588115982 I)', 'GC      .000000000000', 'GC     -.470012389662', 'GC (     .304429918079    -.759588115982 I)', 'GC (    -.597633392156+     .341504795517 I)', 'GC (    -.188354751434    -.322614357177 I)', 'GC (    -.184572576793    -.056386275890 I)', 'GC (    -.001538453194    -.010256312798 I)', 'GC     -.871339083353', 'GC (    -.282010822633    -.376014430178 I)', 'GC (    -.195066832416    -.166632076142 I)', 'GC (    -.080023585072    -.040132442289 I)', 'GC (    -.871339083353 I)', 'GC (     .614292292936+     .230293800534 I)', 'GU     -.177561319578', 'GU     1.872275132772', 'GU     -.287682072452', 'GU      .121199763806', 'GU     1.872275132772', 'GU      .300390433929', 'GU     2.424292039632', 'GU     5.078843955729', 'GU     -.430771890502', 'GU     1.710805699712', 'GU      .121196099916', 'GU     4.280623603702', 'GU     6.385459616123', 'GU     1.710805699712', 'GU     2.308482621335']
_NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def test_gaussian_coulomb_matches_gs():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    got = [l.rstrip() for l in buf.getvalue().splitlines() if l.startswith(('GC', 'GU'))]
    assert len(got) == len(REF)
    for a, b in zip(got, REF):
        assert a[:2] == b[:2]
        va = [float(v) for v in re.findall(_NUM, a[2:])]
        vb = [float(v) for v in re.findall(_NUM, b[2:])]
        assert len(va) == len(vb), (a, b)
        for x, y in zip(va, vb):
            assert abs(x - y) <= 1e-11 * max(abs(y), 1.0), (a, b)


def test_gaussian_coulomb_fitted_matches_gs():
    """GaussianCoulombFitted[x, y, sigx, sigy] = {fx, fy, fu, fxx, fyy, fxy} (tintrb.f twspfu: 21×61×61 の表と 3 次のスプライン補間)。参照は gs (相対 1e-9)"""
    import io
    import re
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=1999;\nPrint["B ",GaussianCoulombFitted[1e-3,2e-3,1e-3,0.5e-3]];\n'
                   'Print["C ",GaussianCoulombFitted[5e-3,1e-3,1e-3,0.5e-3]];\nPrint["F ",GaussianCoulombFitted[2e-4,5e-3,2e-3,0.7e-3]];\nExit;\n', out=buf)
    ref = {'B': [-304.42102299048554, -759.5940473238176, 1.8722712437633469, -235268.35375096588, 233980.23432509732, 200962.37753776336],
           'C': [-394.55546054410104, -84.23024764220699, 3.4023010813079275, 75939.26902299032, -76085.17452104327, 35450.98977903649],
           'F': [-11.795681831580751, -358.21439052516445, 2.2301426243649614, -58782.377051045194, 58866.935355580856, 3616.5705459520077]}
    got = {}
    for l in buf.getvalue().splitlines():
        if l[:2] in ('B ', 'C ', 'F '):
            got[l[0]] = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', l[2:])]
    for k, want in ref.items():
        assert len(got[k]) == 6
        for g, w in zip(got[k], want):
            assert abs(g - w) < 1e-9 * abs(w)
