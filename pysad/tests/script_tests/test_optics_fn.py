"""優先度 高・中の関数・変数 (docs/history/work_trpt_acc.md 11 節): TransferMatrix, BeamMatrix, NormalCoordinates, SurvivedParticles, OrbitGeo,
GeoBase, SymplecticJ, FFS$SHOW, EMIT の DTSYNCH/PHICAV/EFFVC/EFFRFFREQ/OMEGA0, EMITYE, OpticsProlog/Epilog, ParabolaRandom,
AccelerateParticles, SADScript の単位 (1GEV)。参照値は SAD K64 (<work-perf>/hiprio.sad)。"""
import io
import re

import numpy as np

from pysad.ffscmd import run_sad_string

RING = '''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5) DS = (L=0.2);
BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5 ROTATE=3 DEG);
SEXT SF = (L=0.1 K2=2) SD = (L=0.1 K2=-3);
CAVI RFC = (L=0.5 VOLT=2E6 HARM=24);
MARK MK = ();
LINE CELL = (MK QF DS SF D B D QD DS SD D B D);
LINE RING = (7*CELL MK QF DS SF D B D QD DS SD D B RFC D);
FFS;
USE RING;
CELL;
CALC;
'''


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def _f(v):
    return np.array(v, dtype=float)


def test_transfer_matrix_4d_matches_sad():
    ffs, I, text = _run(RING + '''
tm1 = TransferMatrix["QF.2", "QD.5"];
tm2 = TransferMatrix[10, 3];
tm3 = TransferMatrix["QD.5", "QF.2", Fold->True];
tm4 = TransferMatrix[3, 40, Symplectic->True];
''')
    ref1 = [[.675109609983, -.933800671911, .177019726207, .403150304502], [-.810606020024, 2.518831959390, -.212222311402, -.164392441736],
            [-.587338431876, 2.039788920212, 2.569057994712, 1.169883364084], [-.143521684707, .402318569693, .628148006898, .653314608884]]
    np.testing.assert_allclose(_f(I.get('tm1')), ref1, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(_f(I.get('tm2'))[0], [2.235193054906, -3.200844962721, -.126737166168, .043963745002], rtol=1e-9)
    np.testing.assert_allclose(_f(I.get('tm3'))[0], [-.1093486197, 2.697500459952, .280025252204, -2.83326805549], rtol=1e-9)
    np.testing.assert_allclose(_f(I.get('tm4'))[0], [-1.654773043616, -6.938850833964, .256230818135, .620165617897], rtol=1e-9)
    # 4D の転送行列はシンプレクティック (det 1)
    assert abs(np.linalg.det(_f(I.get('tm1'))) - 1.0) < 1e-9


def test_beam_matrix_and_normal_coordinates_match_sad():
    ffs, I, text = _run(RING + '''
EMITX = 2E-9; EMITY = 1E-11; DP = 0.004;
bm1 = BeamMatrix[5];
nc = NormalCoordinates[7];
GAUSS;
bm3 = BeamMatrix[5];
''')
    ref = [1.10309720316E-5, -3.9894528689E-6, -1.8618723567E-6, -3.3577826833E-7, 2.0094495095E-10, 7.66578520900E-6]
    np.testing.assert_allclose(_f(I.get('bm1'))[0], ref, rtol=1e-9, atol=1e-20)
    bm1, bm3 = _f(I.get('bm1')), _f(I.get('bm3'))
    assert abs(bm3[5, 5] / bm1[5, 5] - 3.0) < 1e-12                      # GAUSS: DP² と DP²/3
    nc = I.get('nc')
    np.testing.assert_allclose(_f(nc[0])[0], [.640862562345, 0.0, -.294981797152, .549546289878, 0.0, -.638614978551], rtol=1e-9, atol=1e-12)


def test_beam_matrix_6d_and_transfer_matrix_6d():
    ffs, I, text = _run(RING + '''
CALC6D;
RFSW;
CALC;
EMITX = 2E-9; EMITY = 1E-11; DP = 0.004;
bm6 = BeamMatrix[10];
tm6 = TransferMatrix[10, 40];
''')
    # SAD (<work-perf>/hiprio11.sad.out: UNIFORM): BM6 = r.Diag(.., (DP/√3)²).rᵀ
    np.testing.assert_allclose(_f(I.get('bm6'))[0], [2.41777682627E-6, 2.25985486356E-6, -2.9328841803E-8, -8.5393911507E-8, -3.6958190902E-9, 3.57302329304E-6],
                               rtol=1e-9, atol=1e-20)
    np.testing.assert_allclose(_f(I.get('tm6'))[0], [-.164462773127, -2.872721859337, .112522464626, -.107339172665, 0.0, 3.335573545045], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(_f(I.get('tm6'))[4], [1.037392013566, -1.713080651512, .095451032013, .563898078828, 1.0, -3.789856673856], rtol=1e-9, atol=1e-12)


def test_geo_functions_and_symplectic_j():
    ffs, I, text = _run(RING + '''
gb = GeoBase[{0.1, -0.2, 0.3}];
og = OrbitGeo["QD.3"];
sj = SymplecticJ[6];
''')
    ref = [[.153791997989, -.944702485995, .289629477626], [-.159345079308, -.312991825785, -.936293363584], [.975170327202, .097843395007, -.198669330795]]
    np.testing.assert_allclose(_f(I.get('gb')), ref, rtol=1e-9)
    np.testing.assert_allclose(_f(I.get('og')), [6.483169788248, 8.804785807831, 0.0], rtol=1e-9, atol=1e-12)
    sj = _f(I.get('sj'))
    assert sj.shape == (6, 6) and sj[0, 1] == 1 and sj[1, 0] == -1 and sj[4, 5] == 1 and sj[5, 4] == -1 and abs(sj).sum() == 6


def test_emit_variables_and_emittance_floor():
    ffs, I, text = _run(RING + '''
om = OMEGA0;
emit;
ev = {DTSYNCH, PHICAV, EFFVCRATIO, EFFVC, EFFRFFREQ, EMITX, EMITY, SIGE};
g1 = ?GAUSS;
EMITYE = 1E-7;
emit;
eye = {EMITX, EMITY, EMITZ};
''')
    assert abs(float(I.get('om')) - 44742311.6717202) < 1e-6
    ref = [.004849571625, .017370500636, 1.0, 2000000.0, 170903041.629899, 2.55730893071E-7, 4.39521978048E-8, 6.49875311752E-4]
    np.testing.assert_allclose([float(v) for v in I.get('ev')], ref, rtol=1e-8)
    assert float(I.get('g1')) == 1.0                                      # EMIT コマンドは GAUSS を ON にする (SAD で確認)
    np.testing.assert_allclose([float(v) for v in I.get('eye')], [2.55730893071E-7, 1.74202227029E-7, 1.52790364513E-5], rtol=1e-8)


def test_optics_prolog_epilog_and_ffs_show():
    ffs, I, text = _run(RING + '''
OpticsProlog := Print["PROLOG"];
OpticsEpilog := Print["EPILOG ", Twiss["NX","$$$"]];
calc;
OpticsProlog =.; OpticsEpilog =.;
fit $$$ NX 2.3 NY 1.4;
fit QF.3 BX 7;
fit QD.3 BY 3;
free QF K1 QD K1;
calc;
sh = FFS$SHOW[];
sh2 = FFS["SHOW"];
''')
    assert 'PROLOG' in text
    m = re.search(r'EPILOG\s+([-0-9.Ee+]+)', text)              # Print の数値書式は SAD と異なるので値で照合
    assert m and abs(float(m.group(1)) - 13.857817572298) < 1e-9
    sh = I.get('sh')
    assert [list(c[:3]) for c in sh] == [['$$$', '', 'NX'], ['$$$', '', 'NY'], ['QF.3', '', 'BX'], ['QD.3', '', 'BY']]
    np.testing.assert_allclose([[float(c[3]), float(c[4]), float(c[5])] for c in sh],
                               [[2.3, 1, 6.283185307180], [1.4, 1, 6.283185307180], [7, 1, 1], [3, 1, 1]], rtol=1e-12)
    assert I.get('sh2') == sh


def test_survived_particles_parabola_random_and_accelerate():
    ffs, I, text = _run(RING + '''
x0  = { 1E-3, 0.5, 0,     1E-3,  0    };
px0 = { 0,    0,    0,     1E-4,  0    };
y0  = { 0,    0,    1E-3,  0.9,  0    };
py0 = { 0,    0,    0,     0,     0    };
z0  = { 0,    0,    0,     0,     1E-3 };
dp0 = { 0,    0,    0,     1E-3,  1E-3 };
fl0 = { 1,    1,    1,     1,     1    };
beam0 = {1, {x0, px0, y0, py0, z0, dp0, fl0}};
emit;
NORAD; NOFLUC;
b1 = TrackParticles[beam0, 1, 1, 2];
sv = SurvivedParticles[b1[2]];
SeedRandom[5];
pr = ParabolaRandom[6];
b2 = AccelerateParticles[beam0, 1E9*(1+0.01*n), {n, 3}];
mom = MOMENTUM;
''')
    sv = I.get('sv')
    assert len(sv) == 7 and len(sv[0]) == 3 and all(float(f) == 1.0 for f in sv[6])
    np.testing.assert_allclose([float(v) for v in sv[0]], [-.001846166488, .001744433096, .003043457122], rtol=1e-9)
    np.testing.assert_allclose([float(v) for v in sv[5]], [8.19053000146E-6, 4.26603634674E-6, 8.17338126852E-4], rtol=1e-9)
    np.testing.assert_allclose([float(v) for v in I.get('pr')], [.113686837722, .108501863200, -.341520651709, -.413843683433, .091772363055, .207556743640], rtol=1e-11)
    assert abs(float(I.get('mom')) - 1.03e9) < 1.0                        # AccelerateParticles は MOMENTUM を変えたまま
    b2 = I.get('b2')
    # 生き残った粒子 1 は SAD と 1e-8 で一致 (落ちた粒子の座標の扱いは照合しない)
    assert abs(float(b2[1][0][0]) - (-.039110705781)) < 1e-9 * 4
    assert float(b2[1][6][1]) == 0.0 and float(b2[1][6][0]) == 1.0


def test_sadscript_unit_suffix():
    ffs, I, text = _run(RING + 'a = 2.5GEV; b = 30DEG; c = 1GEV*(1+0.01*2); d = 3MM;\n')
    assert float(I.get('a')) == 2.5e9 and abs(float(I.get('b')) - 30 * np.pi / 180) < 1e-15
    assert float(I.get('c')) == 1.02e9 and float(I.get('d')) == 3e-3
