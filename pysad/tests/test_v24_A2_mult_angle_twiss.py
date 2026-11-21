"""MULT ANGLE≠0 の CALC (Twiss) の gs (K64) との照合 (v2.4.0)。24 セルでなく 6 周期のリング (νx = 2.75037…。βx が大きく鋭いので、
一致は 10 桁で確かめている)。gs の値: F1・FB1・FB2・FRINGE・DISFRIN は Twiss に効かない、ROTATE・E1・K1 は効く。"""
import io

import pytest

SRC = '''MOMENTUM=3GEV;
MARK M0=();
DRIFT D1=(L=1);
QUAD QF=(L=0.3 K1=0.7) QD=(L=0.3 K1=-0.6);
MULT MA=(L=1 ANGLE=0.2617993877991494 K1=0.05);
LINE C=(M0 QF D1 MA D1 QD D1 MA D1);
LINE R=(C C C C C C);
FFS;
$FORM="14.10";
USE R; CELL; CALC;
h[k_,v_]:=(LINE[k,"MA"]=v; FFS["CALC"]; Print["T ",k," ",Twiss["NX","$$$"]/2/Pi," ",Twiss["NY","$$$"]/2/Pi," ",Twiss["BX","M0"]]; LINE[k,"MA"]=0);
h["F1",0.3]; h["FB1",0.5]; h["FB2",0.5]; h["FRINGE",1]; h["DISFRIN",1]; h["ROTATE",0.2]; h["E1",0.2]; h["K1",0.3];
Exit;
'''

GS = {'F1': (2.7503712185, 2.1263061470, 2897.59453799), 'FB1': (2.7503712185, 2.1263061470, 2897.59453799),
      'FB2': (2.7503712185, 2.1263061470, 2897.59453799), 'FRINGE': (2.7503712185, 2.1263061470, 2897.59453799),
      'DISFRIN': (2.7503712185, 2.1263061470, 2897.59453799), 'ROTATE': (2.7504116381, 2.1597735608, 2659.77552252),
      'E1': (2.7504131924, 2.2081555903, 2741.76881808), 'K1': (2.7501876945, .4255465250, 2634.90966087)}


def test_mult_angle_twiss_matches_gs(k64_mode):
    from pysad.ffscmd import run_sad_string
    b = io.StringIO()
    run_sad_string(SRC, out=b)
    got = {}
    for l in b.getvalue().splitlines():
        if l.startswith('T '):
            p = l.split()
            got[p[1]] = tuple(float(x) for x in p[2:5])
    assert set(got) == set(GS)
    for k, v in GS.items():
        for a, c in zip(got[k], v):
            assert a == pytest.approx(c, rel=1e-9)
