"""APERT 要素 (tapert.f / tapert1) の追跡。楕円 AX AY、矩形 DX1 DX2 DY1 DY2、中心のずれ DX DY、回転 ROTATE、既定 0 (全粒子が落ちる)、
JDPX JDPY DP (落とさず境界へ再入射)、複数周、APERT が 2 個ある並び。TrackParticles の戻りの flag と座標 (落ちた粒子は落ちた時点の値)
を SAD K64 (gs) と比べる。参照値は gs (NPARA=1; $FORM="14.9"、粒子 16 個、3 GeV の FODO リング 24 セル、RFSW/RAD なし)。
numba 版 (track6n) と Python 版 (track6) の両方を確かめる。
未確認: JDPX/JDPY の再入射は粒子が 11 個以下のとき gs と同じ乱数列で一致するが、12 個以上では gs が別の乱数列を使う (原因は未調査)。
SOL 領域の中の APERT、TRPT での損失の報告は gs と比べていない。"""
import io
import re

import numpy as np
import pytest

from pysad.ffscmd import run_sad_string

LAT = """
MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
APERT APE = (AX=0.004 AY=0.003);
APERT APR = (DX1=0.003 DX2=-0.002 DY1=-0.004 DY2=0.0025);
APERT APO = (AX=0.004 AY=0.003 DX=0.001 DY=-0.001);
APERT APB = (AX=0.005 AY=0.004 DX1=-0.003 DX2=0.0045 DY1=-0.003 DY2=0.002);
APERT APT = (AX=0.006 AY=0.002 ROTATE=0.6 DX=0.0005);
APERT APZ = ();
APERT APJ = (DX1=-0.002 DX2=0.002 DY1=-0.0015 DY2=0.0015 JDPX=1e-4 JDPY=2e-4 DP=0.001);
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE R0 = (ORIGIN RF 24*CELL);
LINE RE = (ORIGIN APE RF 24*CELL);
LINE RR = (ORIGIN RF 3*CELL APR 21*CELL);
LINE RO = (ORIGIN RF 5*CELL APO 19*CELL);
LINE RB = (ORIGIN RF 24*CELL APB);
LINE RT = (ORIGIN APT RF 24*CELL);
LINE RZ = (ORIGIN RF 2*CELL APZ 22*CELL);
LINE R2 = (ORIGIN APE RF 10*CELL APR 14*CELL);
LINE RJ = (ORIGIN RF 24*CELL APJ);
NPARA = 1;
FFS;
NORFSW; NORAD; NOFLUC;
$FORM = "14.9";
"""
XS = [0, 1e-3, -1e-3, 2e-3, -2e-3, 3.5e-3, -3.5e-3, 0, 0, 0, 0, 2e-3, -2e-3, 4e-3, 2.5e-3, 1e-3]
YS = [0, 0, 0, 0, 0, 0, 0, 1e-3, -1e-3, 2e-3, -2.5e-3, 2e-3, -2e-3, 3e-3, -1.5e-3, 4e-3]

GS_FLAGS = {
    'R0': '1111111111111111',
    'RE': '1111111111111010',
    'RR': '1100000111000000',
    'RO': '1110100110101000',
    'RB': '1111000111110011',
    'RT': '1111110111000010',
    'RZ': '0000000000000000',
    'R2': '1100000110100000',
    'RE3': '1110000111100000',
}
GS_X = {  # (x, y) �S���q (���������q�͗��������_�̍��W)
    'RE': ([0.0, 0.001713765, -0.001713994, 0.003430061, -0.003430962, 0.0060153, -0.006017935, 5.77745611e-07, 5.77745611e-07, 2.31239932e-06, 3.61478068e-06, 0.003435024, -0.003431276, 0.004, 0.004293187, 0.001], [0.0, 2.5456729e-19, -2.541052e-19, 5.0876671e-19, -5.069188e-19, 8.8755744e-19, -8.819022e-19, 0.000159000901, -0.0001590009, 0.000320135476, -0.00040216947, 0.000323476978, -0.00032294953, 0.003, -0.00024302145, 0.004]),
    'RR': ([0.0, 0.001713765, -0.002369553, 0.004743835, -0.004737406, 0.00830548, -0.00828579, 5.77745611e-07, 5.77745611e-07, 2.31239932e-06, -1.209007e-06, 0.004742916, -0.004738036, 0.009490952, 0.005930173, 0.00236777], [0.0, 2.5456729e-19, 2.0325111e-19, -4.090089e-19, 4.055486e-19, -7.176094e-19, 7.0701365e-19, 0.000159000901, -0.0001590009, 0.000320135476, 0.003864927, -0.003093709, 0.00309045, -0.004643691, 0.002320644, -0.006185653]),
    'RO': ([0.0, 0.001713765, -0.001713994, -0.004236274, -0.003430962, -0.007414221, 0.007402489, 5.77745611e-07, 5.77745611e-07, 3.83734977e-07, 3.61478068e-06, -0.004235437, -0.003431276, -0.00847053, -0.005295061, -0.002115406], [0.0, 2.5456729e-19, -2.541052e-19, 9.5700838e-19, -5.069188e-19, 1.6730415e-18, -1.676438e-18, 0.000159000901, -0.0001590009, 0.002315362, -0.00040216947, 0.002322257, -0.00032294953, 0.003496083, -0.00174299, 0.004639605]),
    'RB': ([0.0, 0.001713765, -0.001713994, 0.003430061, -0.003430962, 0.0060153, -0.006017935, 5.77745611e-07, 5.77745611e-07, 2.31239932e-06, 3.61478068e-06, 0.003435024, -0.003431276, 0.00689843, 0.004293187, 0.001728322], [0.0, 2.5456729e-19, -2.541052e-19, 5.0876671e-19, -5.069188e-19, 8.8755744e-19, -8.819022e-19, 0.000159000901, -0.0001590009, 0.000320135476, -0.00040216947, 0.000323476978, -0.00032294953, 0.000504835358, -0.00024302145, 0.000659153424]),
    'RT': ([0.0, 0.001713765, -0.001713994, 0.003430061, -0.003430962, 0.0060153, -0.0035, 5.77745611e-07, 5.77745611e-07, 2.31239932e-06, 0.0, 0.002, -0.002, 0.004, 0.004293187, 0.001], [0.0, 2.5456729e-19, -2.541052e-19, 5.0876671e-19, -5.069188e-19, 8.8755744e-19, 0.0, 0.000159000901, -0.0001590009, 0.000320135476, -0.0025, 0.002, -0.002, 0.003, -0.00024302145, 0.004]),
    'R2': ([0.0, -0.00030377281, -0.002346697, 0.004697562, -0.004692346, 0.008225403, -0.008209402, -1.0272603e-06, -1.0272603e-06, 1.20935756e-06, -6.4210084e-06, 0.004698966, -0.004691332, 0.004, 0.005873827, 0.001], [0.0, 3.4847422e-19, 2.1891178e-19, -4.362011e-19, 4.3891292e-19, -7.635565e-19, 7.7185311e-19, -0.001381731, 0.001381731, 0.002586054, 0.003447732, 0.002588172, -0.0025856, 0.003, -0.001941515, 0.004]),
    'RE3': ([0.0, -0.00030377281, 0.000302791118, -0.004265, 0.004261251, 0.0060153, -0.006017935, -1.0272603e-06, -1.0272603e-06, -4.1092855e-06, -6.4210084e-06, -0.004262562, 0.004255992, 0.004, 0.004293187, 0.001], [0.0, 3.4847422e-19, -3.490848e-19, 1.6452056e-18, -1.644208e-18, 8.8755744e-19, -8.819022e-19, -0.001381731, 0.001381731, -0.002760451, 0.003447732, -0.001458573, 0.001465893, 0.003, -0.00024302145, 0.004]),
}
GS_JDP = {
    'x': [-0.00069791797, 0.001495669, 0.001749059, 0.001793683, 0.001857716, 0.000236831728],
    'y': [-0.001405707, -0.00099582784, -0.00072745802, 0.000663528372, -0.00055565444, -0.001489446],
    'px': [2.53116308e-05, 4.05750355e-05, -1.7543083e-05, 5.37170539e-05, 9.70186035e-05, 7.10144738e-05],
    'd': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
}


def _fmt(v):
    return '{' + ', '.join(repr(float(a)) for a in v) + '}'


def _track(line, nturn, xs=XS, ys=YS, fallback=False):
    """USE line; CELL; CALC; TrackParticles → (flag 文字列, x, y, px, dp)。fallback=True は track6 (Python 版) を強制"""
    n = len(xs)
    z = _fmt([0] * n)
    src = LAT + f'beam = {{{_fmt(xs)}, {z}, {_fmt(ys)}, {z}, {z}, {z}, {_fmt([1] * n)}}};\n'
    src += f'USE {line}; CELL; CALC;\n'
    src += f'r = TrackParticles[{{1, beam}}, 1, 1, {nturn}];\n'
    src += 'Print["F ", r[[2,7]]]; Print["X ", r[[2,1]]]; Print["Y ", r[[2,3]]]; Print["P ", r[[2,2]]]; Print["D ", r[[2,6]]];\n'
    if fallback:
        from pysad import track6n
        old = track6n.TrackRun

        def boom(*a, **k):
            raise NotImplementedError('forced')
        track6n.TrackRun = boom
    try:
        out = io.StringIO()
        run_sad_string(src, out=out)
    finally:
        if fallback:
            track6n.TrackRun = old
    res = {}
    for m in re.finditer(r'^([FXYPD]) \{([^}]*)\}', out.getvalue(), re.M | re.S):
        res[m.group(1)] = [float(v) for v in m.group(2).replace('\n', ' ').split(',')]
    return (''.join(str(int(v)) for v in res['F']), np.array(res['X']), np.array(res['Y']),
            np.array(res['P']), np.array(res['D']))


CASES = [('RE', 1), ('RR', 1), ('RO', 1), ('RB', 1), ('RT', 1), ('RZ', 1), ('R2', 3), ('RE', 3)]


def _key(line, nt):
    return 'RE3' if (line, nt) == ('RE', 3) else line


@pytest.mark.parametrize('line,nt', CASES)
def test_apert_matches_gs(line, nt):
    k = _key(line, nt)
    flag, x, y, _, _ = _track(line, nt)
    assert flag == GS_FLAGS[k]
    if k in GS_X:
        assert np.allclose(x, GS_X[k][0], atol=2e-9, rtol=0)
        assert np.allclose(y, GS_X[k][1], atol=2e-9, rtol=0)


def test_no_apert_all_alive():
    flag, _, _, _, _ = _track('R0', 1)
    assert flag == GS_FLAGS['R0']


def test_apert_default_zero_drops_everything():
    """全キーワード既定 0 の APERT: xl<xa<xh が成り立たないので全粒子が落ちる (gs と同じ。RZ)"""
    flag, _, _, _, _ = _track('RZ', 1)
    assert flag == '0' * 16


def test_apert_python_tracker_matches_numba():
    """track6 (Python 版) の APERT も numba 版・gs と同じ (落ちた粒子の並べ替え、複数の APERT、3 周)"""
    for line, nt in (('R2', 3), ('RT', 1), ('RO', 1)):
        f0, x0, y0, _, _ = _track(line, nt)
        f1, x1, y1, _, _ = _track(line, nt, fallback=True)
        assert f1 == f0 == GS_FLAGS[_key(line, nt)]
        assert np.allclose(x1, x0, atol=1e-12, rtol=0) and np.allclose(y1, y0, atol=1e-12, rtol=0)


def test_apert_jdp_reinjection_matches_gs():
    """JDPX/JDPY ≠ 0: 落とさず、外にいる粒子を境界上 (x=xa cosφ1, y=ya sinφ1, px=JDPX cosφ2, py=JDPY sinφ2) へ乱数で入れ直し、DP を足す。
    参照は gs (粒子 6 個 x {2e-3,-2e-3,3.5e-3,-3.5e-3,2e-3,-2e-3}、y {0,0,0,0,2e-3,-2e-3}、APERT は最後の要素)"""
    xs = [2e-3, -2e-3, 3.5e-3, -3.5e-3, 2e-3, -2e-3]
    ys = [0, 0, 0, 0, 2e-3, -2e-3]
    for fb in (False, True):
        flag, x, y, px, d = _track('RJ', 1, xs, ys, fallback=fb)
        assert flag == '1' * 6
        assert np.allclose(x, GS_JDP['x'], atol=2e-9, rtol=0)
        assert np.allclose(y, GS_JDP['y'], atol=2e-9, rtol=0)
        assert np.allclose(px, GS_JDP['px'], atol=2e-9, rtol=0)
        assert np.allclose(d, GS_JDP['d'], atol=1e-12, rtol=0)
