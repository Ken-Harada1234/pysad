"""ANGLE=0 で K1≠0 の BEND の追跡 (tbend.f:566-590 → tquad。SAD K64 の癖: K1 の符号を失う、11_sad_quirks A38) と
TrackParticles の粒子損失 (tapert: LOSSAMPL/LOSSDZ、落ちた粒子は alive=0 で座標は limittrack で抑えた値)。
参照値は SAD K64 (<work-perf>/bendk1_{plain,dxy,drot,quad}.sad.out)。"""
import io

import numpy as np

from pysad.ffscmd import run_sad_string


def _src(bq, bq_is_quad=False):
    if bq_is_quad:
        els = f'''BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5) BQ2 = (L=0.3 ANGLE=0 K1=0.45);
QUAD QF = (L=0.3 K1=0.6) {bq};'''
    else:
        els = f'''BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5) {bq} BQ2 = (L=0.3 ANGLE=0 K1=0.45);
QUAD QF = (L=0.3 K1=0.6);'''
    return f'''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5) DS = (L=0.2);
{els}
SEXT SF = (L=0.1 K2=2) SD = (L=0.1 K2=-3);
MARK MK = ();
LINE CELL = (MK QF DS SF D B D BQ DS SD D B D);
LINE RING = (7*CELL MK BQ2 DS SF D B D BQ DS SD D B D);
FFS;
USE RING;
CELL;
CALC;
x0  = {{ 1E-3, 3E-3, 0,     1E-3,  0    }};
px0 = {{ 0,    0,    0,     1E-4,  0    }};
y0  = {{ 0,    0,    1E-3,  1E-3,  0    }};
py0 = {{ 0,    0,    0,     0,     0    }};
z0  = {{ 0,    0,    0,     0,     1E-3 }};
dp0 = {{ 0,    0,    0,     1E-3,  1E-3 }};
fl0 = {{ 1,    1,    1,     1,     1    }};
beam0 = {{1, {{x0, px0, y0, py0, z0, dp0, fl0}}}};
NORAD; NOFLUC;
b1 = TrackParticles[beam0, 1, 1, 3];
emit;
RAD; NOFLUC;
b2 = TrackParticles[beam0, 1, 1, 1];
'''


LOST3 = [3.0, .9999, 3.0, .9999, 1e10, 0.0]
LOST4 = [3.0, 1.0008999, 3.0, 1.0008999, 1e10, .001]
LOSTR = [3.0, 9.99900000029E-7, 3.0, 9.99900000029E-7, 1e10, -.999999]
SAD = {
    'plain': ('BQ = (L=0.3 ANGLE=0 K1=-0.5)', False, {
        1: [-8.7963388069E-4, 9.19593687075E-5, 0.0, 0.0, 2.21804024814E-4, 0.0],
        2: [-.002612011476, 2.83988644831E-4, 0.0, 0.0, 5.97540057079E-4, 0.0],
        3: LOST3, 4: LOST4,
        5: [.001249186719, 1.45556390325E-4, 0.0, 0.0, -.009934011618, .001]},
        {1: [-2.9389329757E-4, -6.1825704761E-4, 0.0, 0.0, -2.1177562040E-4, -3.4737866883E-5],
         2: [-8.1245878712E-4, -.001826357922, 0.0, 0.0, -7.7135413310E-4, -3.4738485011E-5],
         3: LOSTR, 4: LOSTR,
         5: [5.80539503511E-4, 3.66263373223E-4, 0.0, 0.0, -.002387727124, 9.65200643616E-4]}),
    'dxy': ('BQ = (L=0.3 ANGLE=0 K1=-0.5 DX=1E-3 DY=5E-4)', False, {
        1: [2.661130223577, -.001099556230, .640077922955, .004366811936, -.145569730751, 0.0],
        2: LOST3, 3: LOST3, 4: LOST4, 5: LOST4},
        {1: [2.548362931280, -.007847542754, .662655648905, .002709641207, -.119686772859, -2.8143062414E-5],
         2: LOSTR, 3: LOSTR, 4: LOSTR, 5: LOSTR}),
    'drot': ('BQ = (L=0.3 ANGLE=0 K1=-0.5 DROTATE=0.01)', False, {
        1: LOST3, 2: LOST3, 3: LOST3, 4: LOST4,
        5: [-3.0, 1.0008999, -3.0, .323691050445, 84.398999080554, .001]},
        {1: LOSTR, 2: LOSTR, 3: LOSTR, 4: LOSTR,
         5: [-3.0, -.146649012415, 3.0, .196422124199, 4.706010223089, .001311245008]}),
    'quad': ('BQ = (L=0.3 K1=-0.5 ROTATE=20 DEG DX=1E-3 DY=5E-4)', True, {
        1: [.003108690192, 1.06481899024E-4, -.001965218918, -6.3416188467E-4, .030141576455, 0.0],
        2: [.006132774227, 2.49944680990E-4, -.003779486512, -5.9507066115E-4, .032709692632, 0.0],
        3: [.002741739483, -1.8227712599E-4, -8.1887681215E-4, -.001126467162, .027819884975, 0.0],
        4: [7.01588812993E-4, -8.9151134765E-5, .002192002848, -4.1601019831E-4, -.012565750047, .001],
        5: [-.001168871755, 2.08747744949E-4, .002173688715, 2.32231549763E-4, -.009897683756, .001]},
        None),
}


def _check(got, ref, name):
    for i, r in ref.items():
        row = [float(got[1][j][i - 1]) for j in range(6)]
        np.testing.assert_allclose(row, r, rtol=1e-8, atol=1e-12, err_msg=f'{name} P{i}')
        alive = float(got[1][6][i - 1])
        # tapert: |x|,|y| < LOSSAMPL (=1) かつ |z| ≤ LOSSDZ (=100) でなければ落ちる (座標は落ちた時点の値)
        lost = not (abs(r[0]) < 1.0 and abs(r[2]) < 1.0 and abs(r[4]) <= 100.0)
        assert alive == (0.0 if lost else 1.0), f'{name} P{i} alive={alive} row={r}'


def test_angle0_k1_bend_tracking_matches_sad_quirk():
    for name in ('plain', 'dxy', 'drot', 'quad'):
        bq, isq, ref3, refr = SAD[name]
        out = io.StringIO()
        ffs, I = run_sad_string(_src(bq, isq), out=out)
        _check(I.get('b1'), ref3, name + ' NORAD T3')
        if refr is not None:
            _check(I.get('b2'), refr, name + ' RAD T1')


def test_bend_k1_sign_is_lost_only_in_tracking():
    # CALC は BEND(ANGLE=0, K1=-0.5) と QUAD(K1=-0.5) で同じ (符号あり)、追跡は BEND が QUAD(K1=+0.5) と同じ (SAD K64 の癖)
    src_b = _src('BQ = (L=0.3 ANGLE=0 K1=-0.5)')
    src_q = _src('BQ = (L=0.3 K1=-0.5)', True)
    src_qp = _src('BQ = (L=0.3 K1=0.5)', True)
    ffs_b, I_b = run_sad_string(src_b, out=io.StringIO())
    ffs_q, I_q = run_sad_string(src_q, out=io.StringIO())
    ffs_qp, I_qp = run_sad_string(src_qp, out=io.StringIO())
    assert abs(ffs_b.twiss_value('NX', ffs_b.nlat) - ffs_q.twiss_value('NX', ffs_q.nlat)) < 1e-12
    assert abs(ffs_b.twiss_value('NX', ffs_b.nlat) - ffs_qp.twiss_value('NX', ffs_qp.nlat)) > 1e-3
    b, qp = I_b.get('b1'), I_qp.get('b1')
    assert [float(v) for v in b[1][0]] == [float(v) for v in qp[1][0]]
