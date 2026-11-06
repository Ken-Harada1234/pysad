"""TRPT の加速 (gammab、tcave / tmulte の acc、qcopymatg)、MARK OFFSET、USE の大文字小文字、K=0 の QUAD/MULT のドリフト
(docs/history/work_trpt_acc.md)。参照値は SAD K64 (<work-perf>/linac_tiny.sad.out, linac_tiny2.sad.out)。"""
import io
import math

import numpy as np
import pytest

from pysad.ffscmd import run_sad_string
from pysad.sadscript import SADError

LINAC = '''
MOMENTUM = 20MEV;
DRIFT D1 = (L=0.5) D2 = (L=0.3);
CAVI C1 = (L=0.6 VOLT=8E6 PHI=-1.5707963267948966 FREQ=2.856E9);
MULT AC1 = (L=1.5 VOLT=12E6 PHI=-1.5707963267948966 FREQ=2.856E9 K1=0) AC0 = (L=0.25 PHI=-1.5707963267948966 FREQ=2.856E9 K1=0);
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

# SAD: Twiss[{"DX","DPX","DY","DPY","DZ","DDP","AX","BX","AY","BY"}, k] (<work-perf>/linac_tiny2.sad.out)
SAD_ROWS = {
    5: [4.86550036061E-5, -1.1292851659E-4, -3.1349146988E-4, -2.5466140660E-4, -1.5405310945E-8, .000000000000,
        5.229735523561, 2.979064168042, -9.336191167182, 12.070387549647],
    6: [-1.6602911843E-5, -9.6194418644E-5, -3.8936947830E-4, -1.9995387875E-4, -2.0980465491E-8, 1.54866190982E-8,
        -.546517508742, .132839521365, -12.555384878086, 25.282377221801],
    8: [-8.8748726485E-5, -9.6194418644E-5, -5.3933488873E-4, -1.9995387875E-4, -3.9439693837E-8, 1.54866190982E-8,
        -7.878756050547, 6.451794807000, -17.261361372194, 47.644938007790],
    12: [-2.6961267063E-4, -2.1900450971E-4, -7.9172630417E-4, 2.85881216632E-4, -6.2885953693E-8, 3.72900765822E-8,
         -84.149974391683, 102.277951926786, 52.235793440617, 143.634187084369],
    14: [-4.2764567010E-4, -3.7843456453E-4, -5.0589926573E-4, 7.03656705028E-4, -1.8510448005E-7, 3.72900765822E-8,
         -230.18310698879, 258.899433065078, 81.650886555468, 58.440181536793],
}
SAD_END = [-6.8208887506E-4, 8.05171557357E-4, .001389030550, .001585325192, 1.84210120292E-4, 3.72900765822E-8,
           781.371489235087, 662.826191422139, -509.42380263396, 446.687507432669]
KEYS = ('DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AX', 'BX', 'AY', 'BY')


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def _tw(ffs, k):
    return [ffs.twiss_value(key, k) for key in KEYS]


def test_use_mixed_case_line_name():
    ffs, I, text = _run(LINAC)
    assert 'Missing beamline' not in text
    assert ffs.name == 'LinacT' and ffs.nlat == 20


def test_linac_acceleration_matches_sad():
    ffs, I, text = _run(LINAC)
    # 加速前 (空洞の入口まで) は 1e-10、空洞 (CAVI C1, MULT AC1) の後も 1e-10 で一致 (写像も解析ヤコビアンも tcave/tmulte の直訳)
    # (DZ だけ 1e-11 の差: SAD の tdrife0 (dvon=.true.) の dv 項 (δ≠0 の軌道の速度差) を pysad の 4D 経路のドリフトは入れない。追わない)
    for k, ref in SAD_ROWS.items():
        np.testing.assert_allclose(_tw(ffs, k), ref, rtol=1e-9, atol=1e-10, err_msg=f'row {k}')
    # F1 付き BEND を δ≠0 の軌道で通った後は 1e-9 (BEND のフリンジの経路の差、追わない)
    np.testing.assert_allclose(_tw(ffs, ffs.nlat), SAD_END, rtol=3e-8, atol=1e-10)
    # 基準運動量: 20 MeV → C1 (+8 MV) → AC1 (+12 MV)
    from pysad import trptacc
    gb = trptacc.gammab_profile(ffs.els, 20e6 / 0.51099895e6)
    assert abs(gb[0] * 0.51099895e6 - 20e6) < 1e-3
    assert gb[5] > gb[4] and gb[8] > gb[7] and gb[-1] == gb[8]
    e_end = np.sqrt(1 + gb[-1]**2) * 0.51099895e6
    assert abs(e_end - (np.sqrt(20e6**2 + 0.51099895e6**2) + 20e6)) < 1.0      # エネルギーは 20 MV 分増える


def test_mark_offset_positions():
    ffs, I, text = _run(LINAC)
    # MM (OFFSET=-.5): 直前の QD の中央、MC (OFFSET=1.25): 直後の QF の 1/4 (SAD の LINE["S"] と Twiss)
    assert abs(ffs.line_value('S', 12) - 4.8) < 1e-12
    assert abs(ffs.twiss_value('S', 12) - 4.8) < 1e-12
    assert abs(ffs.twiss_value('S', 16) - 6.825) < 1e-12
    assert abs(ffs.twiss_value('S', ffs.nlat) - 7.35) < 1e-12
    ref_mc = [-113.53339183081, 1442.59862803150, -132.73127296013, 94.623634397591]
    np.testing.assert_allclose([ffs.twiss_value(k, 16) for k in ('AX', 'BX', 'AY', 'BY')], ref_mc, rtol=3e-8)


def test_no_acceleration_without_trpt():
    # RING (TRPT OFF) の CALC では空洞は RFSW に依らずドリフト (qcav: rfsw = rfsw .and. trpt)
    ffs, I, text = _run(LINAC.replace('TRPT;', 'RING;'))
    bx = ffs.twiss_value('BX', 6)
    r5 = ffs.row(5)
    drift = r5.betx - 2 * r5.alfx * 0.6 + (1 + r5.alfx**2) / r5.betx * 0.36
    assert abs(bx - drift) < 1e-6 * drift        # 軌道 (DX=1e-4) 上の厳密ドリフトなので px² 程度の差は残る


def test_strengthless_quad_and_mult_drift_off_axis():
    """K1=0 の QUAD と K が全部 0 の MULT は軌道がずれていても厳密ドリフト (以前は tqente の非線形項だけで直進しなかった)"""
    from pysad import lattice as L
    from pysad.elements import Quadrupole, Multipole
    c = np.array([1e-4, 1e-3, -2e-4, 5e-4, 0.0, 0.0])
    for el in (Quadrupole(name='Q', L=1.0, K1=0.0), Multipole(name='M', L=1.0)):
        M, co = L._el_cod_step(el, c, 0.0)
        pz = np.sqrt(1 - 1e-6 - 2.5e-7)
        assert abs(co[0] - (1e-4 + 1e-3 / pz)) < 1e-15 and abs(co[2] - (-2e-4 + 5e-4 / pz)) < 1e-15
        assert abs(M[0, 1] - 1.0) < 1e-5 and abs(M[2, 3] - 1.0) < 1e-5


def test_disp_value_column_for_zero_length_drift():
    ffs, I, text = _run('''
MOMENTUM = 1GEV;
DRIFT D0 = (L=0) D1 = (L=1);
QUAD Q0 = (L=0.2 K1=0);
MARK M1 = (BX=5 BY=5);
LINE L1 = (M1 D0 D1 Q0);
FFS;
USE L1;
INS;
CALC;
disp;
''')
    # tfdisp: 既定キーのある要素は値が 0 でも ' .0'、MARK は ' 0'
    lines = [ln for ln in text.splitlines() if ' D0 ' in ln or ' Q0 ' in ln or ' M1 ' in ln]
    assert any(' D0        .0      .0  ' in ln for ln in lines), lines
    assert any(' Q0        .20000  .0  ' in ln for ln in lines), lines
    assert any(' M1        .0     0   ' in ln for ln in lines), lines


def test_line_gammabeta_follows_acceleration():
    # tftwiss.f:837: LINE["GAMMABETA", k] = gammab(k) (p/mc)。TRPT の加速があれば要素ごとに増え、RING では一定
    ffs, I, text = _run(LINAC)
    gb1 = ffs.line_value('GAMMABETA', 1)
    gbe = ffs.line_value('GAMMABETA', ffs.nlat)
    assert abs(gb1 - 20e6 / 0.51099895e6) < 1e-9
    assert gbe > gb1 * 1.9                                   # 20 MeV → 40 MeV
    assert abs(ffs.line_value('GAMMA', 1) - math.sqrt(1 + gb1**2)) < 1e-12
    # 分数位置は補間 (C1 の中)
    assert ffs.line_value('GAMMABETA', 5) < ffs.line_value('GAMMABETA', 5, 0.5) < ffs.line_value('GAMMABETA', 6)
    ffs, I, text = _run(LINAC.replace('TRPT;', 'RING;'))
    assert ffs.line_value('GAMMABETA', ffs.nlat) == ffs.line_value('GAMMABETA', 1)
