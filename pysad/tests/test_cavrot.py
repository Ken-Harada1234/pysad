"""CAVI の DX/DY/ROTATE (tcave.f の tchge) と VOLT 付き MULT の DX/DY/ROTATE (tmulte.f の tsolrote) (docs/history/work_trpt_acc.md 9 節)。
参照値は SAD K64 (<work-perf>/cavrot.sad.out, linac_rot.sad.out)。"""
import io

import numpy as np

from pysad.ffscmd import run_sad_string

RING = '''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5) DS = (L=0.2);
BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
SEXT SF = (L=0.1 K2=2) SD = (L=0.1 K2=-3);
CAVI RFC = (L=0.5 VOLT=2E6 HARM=24 DX=2E-3 DY=-1E-3 ROTATE=30 DEG);
MARK MK = ();
LINE CELL = (MK QF DS SF D B D QD DS SD D B D);
LINE RING = (7*CELL MK QF DS SF D B D QD DS SD D B RFC D);
FFS;
USE RING;
CELL;
CALC6D;
RFSW;
CALC;
emit;
'''

LINAC = '''
MOMENTUM = 20MEV;
DRIFT D1 = (L=0.5) D2 = (L=0.3);
CAVI C1 = (L=0.6 VOLT=8E6 PHI=-1.5707963267948966 FREQ=2.856E9 DX=1.5E-3 DY=-0.5E-3 ROTATE=20 DEG);
MULT AC1 = (L=1.5 VOLT=12E6 PHI=-1.5707963267948966 FREQ=2.856E9 K1=0 DX=-1E-3 DY=2E-3 ROTATE=-15 DEG) AC0 = (L=0.25 PHI=-1.5707963267948966 FREQ=2.856E9 K1=0);
QUAD QF = (L=0.3 K1=1.2) QD = (L=0.3 K1=-1.1) Q0 = (L=0.2 K1=0);
BEND B1 = (L=1.0 ANGLE=0.3 F1=0.06 FRINGE=1 E1=0.5 E2=0.5);
MARK M1 = (AX=-1.2 BX=8 AY=0.7 BY=6 DX=1E-4 DY=-2E-4) MM = (OFFSET=-0.5) MC = (OFFSET=1.25) M2 = ();
LINE LinacT = (M1 D1 QF D2 C1 D1 AC0 AC1 Q0 D1 QD MM D2 B1 D1 MC QF D2 M2);
FFS;
USE LinacT;
TRPT;
INS;
CALC;
'''
KEYS = ('DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AX', 'BX', 'AY', 'BY')
SAD_ROWS = {
    6: [2.36047128680E-4, -3.6043638782E-5, -4.7358615375E-4, -2.2000413543E-4, -6.3682868462E-8, 2.26147172100E-7,
        -.546517440612, .132839527258, -12.555382796739, 25.282373925743],
    8: [2.09014405036E-4, -3.6043638782E-5, -6.3858922211E-4, -2.2000413543E-4, -8.2264257763E-8, 2.26147172100E-7,
        -7.878753614344, 6.451791536462, -17.261357416891, 47.644925731388],
    12: [-9.9088939174E-5, -1.0230517822E-4, -5.3456904082E-4, 1.62230965350E-4, -1.2414363913E-7, 6.48235935818E-7,
         -84.149948705414, 102.277924774895, 52.235692061594, 143.634009268578],
    14: [-1.6730387208E-4, -1.6180564563E-4, -3.5499632185E-4, 4.45571883770E-4, -1.6664181286E-7, 6.48235935818E-7,
         -230.18320191935, 258.899204550447, 81.650824164406, 58.440290298892],
}


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def test_rotated_offset_cavity_emit_matches_sad():
    ffs, I, text = _run(RING)
    g = I.globals
    # SAD: {EMITX, EMITY, EMITZ, SIGE, SIGZ} (<work-perf>/cavrot.sad.out)。回した・ずらした空洞が x-y 結合を作り EMITY≠0、Beam tilt .2927
    ref = [2.43939537896E-7, 1.73971983835E-7, 1.55573994265E-5, 6.53617781988E-4, .024117497099]
    np.testing.assert_allclose([g['EMITX'], g['EMITY'], g['EMITZ'], g['SIGE'], g['SIGZ']], ref, rtol=1e-9)
    assert 'Beam tilt              =  .29272601 rad' in text
    # DX/DY/ROTATE を外すと平らなリング (EMITY=0)
    ffs2, I2, text2 = _run(RING.replace(' DX=2E-3 DY=-1E-3 ROTATE=30 DEG', ''))
    assert abs(I2.globals['EMITY']) < 1e-20 and 'Beam tilt              =  .00000000 rad' in text2


def test_rotated_offset_cavity_and_mult_trpt_match_sad():
    ffs, I, text = _run(LINAC)
    # DZ だけ 1.5e-10 の差 (δ≠0 の軌道のドリフトの dv 項: test_trpt_acc と同じ、追わない)
    for k, ref in SAD_ROWS.items():
        got = [ffs.twiss_value(key, k) for key in KEYS]
        np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-9, err_msg=f'row {k}')
