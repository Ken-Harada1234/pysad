"""POL (スピン): Emittance[] の PolarizationTime / SpinTune / PolarizationVector / SpinPrecessionMatrix と EMIT の報告。
参照は SAD K64 (gs) の出力 (local/wake/pol1.sad, pol2.sad, pol5.sad)。

Polarization (平衡偏極度) は K64 では NaN (spdepol の cerfc(cb) が |cb|>18 の漸近式で exp(-cb²)=0、exp(cb²)=∞ の積が NaN: docs/19 B55)。
pysad も K64 の式をそのまま写すので NaN になる (どちらも NaN) ことを確かめる。"""
import io
import math
import re

import numpy as np

CELLS = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY=1e-4) QDE = (L=0.3 K1=-0.6 DY=-1.5e-4 ROTATE=0.001);
QUAD QT = (L=0.2 K1=0.02 ROTATE=0.5);
SEXT SF = (L=0.1 K2=2.0 DY=1e-4);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND BR = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747 ROTATE=0.0015);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D3 SF D3 B D2 QDE D1);
LINE CELL3 = (QF D2 BR D2 QT D1 QD D1);
'''
PLANAR = CELLS + '''LINE RING = (ST 12*CELL1 RFC);
'''
COUPLED = CELLS + '''LINE RING = (ST 8*CELL1 CELL2 CELL3 CELL2 CELL3 RFC);
'''
TAIL = '''FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW; RAD; NOFLUC;
POL;
'''


def _run(lat, body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(lat + TAIL + body + 'Exit;\n', out=buf)
    return buf.getvalue()


def _nums(s):
    return [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|NaN', s.replace('NaN', 'nan'))] if False else \
        [float(x) for x in re.findall(r'nan|[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', s.replace('NaN', 'nan'))]


def _out(lat):
    out = _run(lat, 'em = Emittance[];\nPrint["OUT ", {PolarizationTime, SpinTune, NominalSpinTune, Polarization, PolarizationVector}/.em];\n')
    line = [l for l in out.splitlines() if l.startswith('OUT ')][0]
    v = _nums(line[4:])
    return v[0], v[1], v[2], v[3], np.array(v[4:7])


def test_planar_ring_matches_gs():
    tau, nu, nnu, pol, vec = _out(PLANAR)
    assert abs(tau - 123.8154010596169) < 1e-9
    assert abs(nu - 0.40407272052341) < 1e-12
    assert abs(nnu - 6.80814822478217) < 1e-12
    assert math.isnan(pol)
    assert np.allclose(vec, [0.0, 1.0, 0.0], atol=1e-12)


def test_coupled_ring_matches_gs():
    """縦ずれ・傾いた四極と偏向・六極のあるリング: 偏極時間・スピンチューン・偏極ベクトルが gs と一致する"""
    tau, nu, nnu, pol, vec = _out(COUPLED)
    assert abs(tau - 130.3526952378181) < 1e-8
    assert abs(nu - 0.40403473710144) < 1e-11
    assert math.isnan(pol)
    assert np.allclose(vec, [-1.35258874105e-4, 0.99953579967395, 0.03046583127374], atol=1e-10)


def test_emit_report_and_spin_precession_matrix_match_gs():
    out = _run(COUPLED, 'EMIOUT;\nEMIT;\n')
    lines = out.splitlines()
    assert 'Nominal spin tune      =  6.8081482     Polarization time      =  2.1725449 min' in lines
    assert '   Polarization vector at entrance:' in lines
    assert '          x : -1.352589E-04   y :  0.999536       z :  3.046583E-02' in lines
    assert 'Spin tune              =  .40403474     Equil. polarization    =        NaN %' in lines
    assert '        sx -0.8236564      1.7030292E-02 -0.5668335    ' in lines
    assert '        sy -1.7523395E-02  0.9983073      5.5456713E-02' in lines
    assert '        sz  0.5668184      5.5610121E-02 -0.8219637    ' in lines
    assert '        gx -0.3963091      0.1842614      -13.26426       8.737059       0.000000       1.151410    ' in lines


def test_pol_off_leaves_spin_outputs_zero():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(PLANAR + TAIL.replace('POL;', 'NOPOL;') + 'em = Emittance[];\nPrint["OUT ", {PolarizationTime, SpinTune}/.em];\nExit;\n', out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('OUT ')][0]
    assert _nums(line[4:]) == [0.0, 0.0]


def test_quad_offset_jacobian_in_python_path():
    """スロー経路 _tquade: ずれ (DX, DY) だけで回転なしの四極のヤコビアンはずらした座標で評価する (以前は未評価で T が 3e-6 ずれた)"""
    from pysad import emit6d as E
    from pysad.elements import Quadrupole
    ph = E._Phys(3e9, 48.0)
    for kw in (dict(K1=0.7, DY=1e-4), dict(K1=0.7, DX=2e-4, DY=1e-4)):
        el = Quadrupole(name='Q', L=0.3, **kw)
        res = []
        for fast in (True, False):
            st = E._State(np.array([5e-4, 3e-4, 3.5e-3, -9e-4, 0, -4e-5]), ph)
            st.setr0()
            if fast:
                assert E._quade_fast(st, el, True) == 0
            else:
                E._tquade(st, el, True)
            res.append(st)
        assert np.allclose(res[0].T, res[1].T, atol=1e-12)
        assert np.allclose(res[0].D, res[1].D, rtol=1e-9, atol=1e-18)


def test_mult_ring_matches_gs():
    """MULT (K1・K2・K3 とずれ) のあるリング: gs (local/wake/pm2.sad, pm3.sad) と偏極時間・スピンチューン・偏極ベクトルが一致する。
    K3 が強くずれも大きい場合 (pm1.sad) は pysad の閉軌道が gs と違う解に行って NaN になる (未解決、docs/history/work_unimpl.md)"""
    lat = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
MULT MM = (L=0.1 K1=0.05 K2=2.0 DY=1e-4);
MULT MQ = (L=0.2 K1=0.4 K3=300 DY=-3e-5 DX=2e-5);
LINE CELL1 = (QF D2 B D2 MM QD D1 MQ);
LINE RING = (ST 12*CELL1 RFC);
'''
    tau, nu, nnu, pol, vec = _out(lat)
    assert abs(tau - 133.6117332301700) < 1e-8
    assert abs(nu - 0.40379908026905) < 1e-11
    assert np.allclose(vec, [4.722099251248e-6, 0.99999998116143, -1.94048561191e-4], atol=1e-10)


def test_solenoid_ring_matches_gs():
    """SOL (GEO は 1 つだけ) を含むリング: 偏極時間・スピンチューン・偏極ベクトルが gs (local/wake/ps8.sad) と一致する。
    K64 は SOL (BOUND≠0) の長さを周長に入れない (周回周波数・偏極時間・MomentumCompaction に効く) ので pysad も入れない。
    GEO がない SOL では K64 の S 座標が途切れて周長がさらに狂う (偏極時間が半分になる) が、pysad は再現しない"""
    lat = CELLS.replace('MARK ST = ();', 'SOL SO = (L=0.3 BZ=0.05 BOUND=1 GEO=1);\nSOL SO2 = (L=0 BZ=0 BOUND=1);\nMARK ST = ();') + \
        'LINE SC = (SO QT SO2);\nLINE RING = (ST 6*CELL1 SC 6*CELL1 RFC);\n'
    tau, nu, nnu, pol, vec = _out(lat)
    assert abs(tau - 124.3584681268823) < 1e-8
    assert abs(nu - 0.40407273290258) < 1e-11
    assert np.allclose(vec, [0.0, 0.99999986279488, -5.23841797019e-4], atol=1e-10)


def test_thin_mult_ring_spin_tune_matches_gs():
    """L=0 の MULT (薄いキック) のあるリング。tmulke の bsir0 更新は al≠0 だけ (al=0 で 0 割りになっていた)。参照は gs (local/wake/pt7b.sad)"""
    lat = CELLS.replace('MARK ST = ();', 'MARK ST = ();\nMULT TH = (L=0 K2=0.8 DY=5e-5);') + 'LINE CELL1T = (QF D2 B D2 TH QD D1);\nLINE RING = (ST 12*CELL1T RFC);\n'
    out = _run(lat, 'em = Emittance[];\nPrint["OUT ", {SpinTune, NominalSpinTune}/.em];\n')
    v = _nums([l for l in out.splitlines() if l.startswith('OUT ')][0][4:])
    assert abs(v[0] - 0.40407273544612) < 1e-11
