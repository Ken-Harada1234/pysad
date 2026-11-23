# -*- coding: utf-8 -*-
"""v2next 担当 E 項目 2: ANGLE=0 の BEND (補正磁石、tsteee) と MULT の F1 (四極の線形フリンジ tqlfre・FRINGE・FB1・FB2) の Emittance[]・CALC。

(1) k64 モード: gs (K64 の 10-19 ビルド) と全桁一致 (gs の表示 10 桁)。条件: 3 GeV、周長 45.6 m の 12 セル FODO (QF D2 B D2 QD D1)。
    ステアリング (L=0.2 m、K0 = ±8e-4 を交互に、F1・FRINGE=−2..3・FB1・FB2・AE1・AE2 の組み合わせ) と、
    QF・QD を MULT (K1 = 0.7・−0.6、F1 などを付ける) に替えたリング。台本は example/01_lattice_elements/elem_BEND_steer_F1.sad・elem_MULT_F1.sad。
(2) B108 (ANGLE=0 の BEND の端の角 AE1・AE2): K64 は端の蹴りを偏向角に入れるが場の中の行路を変えないので、端の項が 2 倍になる
    (ANGLE≠0 の BEND の B75 と同じ構造)。v2 は行路を x·tan(e) だけ短くする。独立な値 = 硬い端の一様磁場の極 (円弧と直線の厳密な幾何) の
    ロビンソンの定理 Je = 2 + dln(∮B² dl)/dδ。極 3 個の組 (K0・−2K0・K0。軌道が極の外で 0 に戻る) を 24 セルの FODO に入れる。
"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat

STEER_CASES = {
    'v1': ('K0=0.0008 F1=0.06 FRINGE=3', 'K0=-0.0008 F1=0.06 FRINGE=3', {'calc': [3.1335807272, 2.2665026666], 'tunes': [0.1335415957, 0.2665186605, -0.0156476727], 'u0': 781688.042159, 'alpha': 0.0248981165, 'eqpos': 0.0330388023, 'damp': [-0.00012536927, -0.00013028143, -0.00026539875], 'emit': [1.96530549e-07, 1.80194107e-05]}),
    'v2': ('K0=0.0008 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01', 'K0=-0.0008 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01', {'calc': [3.1335808642, 2.2665024001], 'tunes': [0.1335420095, 0.2665186693, -0.0156476735], 'u0': 781687.534306, 'alpha': 0.0248981179, 'eqpos': 0.0330387803, 'damp': [-0.0001253694, -0.00013028136, -0.00026539836], 'emit': [1.96531005e-07, 1.80194239e-05]}),
    'v3': ('K0=0.0008 F1=0.06 FRINGE=1', 'K0=-0.0008 F1=0.06 FRINGE=2', {'calc': [3.1335807272, 2.2665026666], 'tunes': [0.1335415957, 0.2665186605, -0.0156476727], 'u0': 781688.042159, 'alpha': 0.0248981165, 'eqpos': 0.0330388023, 'damp': [-0.00012536927, -0.00013028143, -0.00026539875], 'emit': [1.96530549e-07, 1.80194107e-05]}),
    'v4': ('K0=0.0008 F1=0.06 FRINGE=-1', 'K0=-0.0008 F1=0.06 FRINGE=-2', {'calc': [3.1335807272, 2.2664941086], 'tunes': [0.1335416033, 0.2665100869, -0.0156476533], 'u0': 781690.228063, 'alpha': 0.0248980598, 'eqpos': 0.0330388968, 'damp': [-0.00012536965, -0.00013028179, -0.00026539946], 'emit': [1.96529561e-07, 1.80193449e-05]}),
    'v7': ('K0=0.0008 FRINGE=-2', 'K0=-0.0008 FRINGE=-1', {'calc': [3.1335807272, 2.2664949867], 'tunes': [0.133541611, 0.2665110169, -0.0156476693], 'u0': 781692.414007, 'alpha': 0.0248981158, 'eqpos': 0.0330389914, 'damp': [-0.00012537002, -0.00013028216, -0.00026540019], 'emit': [1.96529413e-07, 1.8019319e-05]}),
    'v8': ('K0=0.003 F1=0.08 FRINGE=3', 'K0=-0.003 F1=0.08 FRINGE=3', {'calc': [3.1337043264, 2.2668052098], 'tunes': [0.1338957163, 0.2670302115, -0.0156422917], 'u0': 782388.4751, 'alpha': 0.02488262, 'eqpos': 0.033069113, 'damp': [-0.00012606038, -0.00013040776, -0.00026504736], 'emit': [1.95375322e-07, 1.80342196e-05]}),
}
MULT_CASES = {
    'm1': ('F1=0.05 FRINGE=3', 'F1=0.05 FRINGE=3', {'calc': [2.7046079095, 1.8643910287], 'tunes': [-0.2954682507, -0.1356089713, -0.0179324755], 'u0': 781632.856818, 'alpha': 0.0326998112, 'eqpos': 0.0298900891, 'damp': [-0.00012428006, -0.00013026822, -0.00026645667], 'emit': [2.68456802e-07, 1.86272495e-05]}),
    'm2': ('F1=0.05 F2=0.02 FRINGE=3', 'F1=0.04 F2=0.03 FRINGE=3', {'calc': [2.6791309296, 1.8839998954], 'tunes': [-0.3209439361, -0.1160001046, -0.017783914], 'u0': 781632.856818, 'alpha': 0.0321602527, 'eqpos': 0.0298900891, 'damp': [-0.00012438316, -0.00013026786, -0.00026635308], 'emit': [2.66007585e-07, 1.84793035e-05]}),
    'm3': ('F1=0.05 FRINGE=1', 'F1=0.05 FRINGE=2', {'calc': [2.7053477284, 1.8652245546], 'tunes': [-0.2947283686, -0.1347754454, -0.0179262099], 'u0': 781632.856818, 'alpha': 0.0326769649, 'eqpos': 0.0298900891, 'damp': [-0.00012428448, -0.00013026821, -0.00026645225], 'emit': [2.68190098e-07, 1.86210881e-05]}),
    'm5': ('F1=0.05 FRINGE=3 K0=0.0004', 'F1=0.05 FRINGE=3 K0=-0.0004', {'calc': [2.7045594712, 1.8643438989], 'tunes': [-0.2978777303, -0.1379282973, -0.0179621035], 'u0': 783085.550645, 'alpha': 0.0328123994, 'eqpos': 0.0299469681, 'damp': [-0.00011848344, -0.00013039093, -0.00027309632], 'emit': [2.83730389e-07, 1.82733654e-05]}),
    'm6': ('F1=0.05 FRINGE=3 K0=0.0004 FB1=0.03 FB2=0.02', 'F1=0.05 FRINGE=3 K0=-0.0004 FB1=0.03 FB2=0.02', {'calc': [2.7045594712, 1.8643437672], 'tunes': [-0.2978777314, -0.13792843, -0.0179621001], 'u0': 783085.551323, 'alpha': 0.0328123869, 'eqpos': 0.0299469681, 'damp': [-0.00011848344, -0.00013039093, -0.00027309632], 'emit': [2.83730292e-07, 1.82733619e-05]}),
    'm8': ('F1=0.05 FRINGE=3 K2=20', 'F1=0.05 FRINGE=3 K2=-20', {'calc': [2.7046079095, 1.8643910287], 'tunes': [-0.2954682507, -0.1356089713, -0.0179324755], 'u0': 781632.856818, 'alpha': 0.0326998112, 'eqpos': 0.0298900891, 'damp': [-0.00012428006, -0.00013026822, -0.00026645667], 'emit': [2.68456802e-07, 1.86272495e-05]}),
    'm10': ('F1K1F=0.02 F1K1B=0.04 F1=0.03 FRINGE=3', 'F1K1F=0.02 F1K1B=0.04 F1=0.03 FRINGE=3', {'calc': [2.7035775109, 1.8632145023], 'tunes': [-0.2964987418, -0.1367854977, -0.0179416558], 'u0': 781632.856818, 'alpha': 0.0327333004, 'eqpos': 0.0298900891, 'damp': [-0.00012427353, -0.00013026823, -0.0002664632], 'emit': [2.68842541e-07, 1.86362752e-05]}),
    'm11': ('K0=0.0004 FRINGE=3', 'K0=-0.0004 FRINGE=3', {'calc': [2.706709397, 1.8668251839], 'tunes': [-0.2957345968, -0.1354527349, -0.017946361], 'u0': 783087.591686, 'alpha': 0.0327549153, 'eqpos': 0.029947048, 'damp': [-0.0001185008, -0.00013039108, -0.00027308015], 'emit': [2.82982651e-07, 1.82586342e-05]}),
}
BASE_STEER = "{'calc': [3.1335807272, 2.26650439], 'tunes': [0.133541611, 0.2665203978, -0.0156476693], 'u0': 781692.414007, 'alpha': 0.0248981158, 'eqpos': 0.0330389914, 'damp': [-0.00012537002, -0.00013028216, -0.00026540019], 'emit': [1.96529413e-07, 1.8019319e-05]}"
BASE_MULT = "{'calc': [2.7067578519, 1.8668722975], 'tunes': [-0.2933181479, -0.1331277025, -0.0179166854], 'u0': 781632.856818, 'alpha': 0.0326422504, 'eqpos': 0.0298900891, 'damp': [-0.0001242913, -0.00013026821, -0.00026644543], 'emit': [2.67760673e-07, 1.86117406e-05]}"


def _floats(s):
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', s)]


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _ = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


_RING_TAIL = '''FFS;
$FORM="14.10";
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
RFSW; RAD; NOFLUC;
e = Emittance[];
Print["Tunes ", Tunes/.e];
Print["U0 ", EnergyLossU0/.e, " alpha ", MomentumCompaction/.e, " EqPos ", EquilibriumPosition/.e];
Print["Emittances ", Emittances/.e];
Print["Damp ", DampingRate/.e];
Exit;
'''

_STEER = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) DS = (L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND SA = (L=0.2 ANGLE=0 %s);
BEND SB = (L=0.2 ANGLE=0 %s);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELLA = (QF D2 B D2 QD DS SA D1);
LINE CELLB = (QF D2 B D2 QD DS SB D1);
LINE CC = (CELLA CELLB);
LINE RING = (ST 6*CC RFC);
''' + _RING_TAIL

_MULT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
MULT QF = (L=0.3 K1=0.7 %s);
MULT QD = (L=0.3 K1=-0.6 %s);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
''' + _RING_TAIL


def _values(out):
    d = {}
    for l in out.splitlines():
        if l.startswith('CALC tunes'):
            d['calc'] = _floats(l[10:])
        elif l.startswith('Tunes'):
            d['tunes'] = _floats(l[5:])
        elif l.startswith('U0'):
            x = _floats(l[2:])
            d['u0'], d['alpha'], d['eqpos'] = x[:3]
        elif l.startswith('Emittances'):
            x = _floats(l[10:])
            d['emit'] = [x[0], x[2]]
        elif l.startswith('Damp'):
            d['damp'] = _floats(l[4:])
    return d


def _check(d, ref):
    for k in ('calc', 'tunes', 'damp'):
        assert np.allclose(d[k], ref[k], rtol=3e-9, atol=2e-10), (k, d[k], ref[k])
    assert np.allclose(d['emit'], ref['emit'], rtol=3e-8), (d['emit'], ref['emit'])
    for k in ('u0', 'alpha', 'eqpos'):
        assert abs(d[k] / ref[k] - 1.0) < 3e-9, (k, d[k], ref[k])


@pytest.mark.parametrize('case', list(STEER_CASES))
def test_steering_bend_f1_matches_gs_k64(k64_mode, case):
    """ANGLE=0 の BEND (tsteee) の F1・FB1・FB2・FRINGE・AE1・AE2 の Emittance[] と CALC が gs と全桁一致 (K0 = ±8e-4 を交互に 12 台)"""
    a, b, ref = STEER_CASES[case]
    _f, out = _run(_STEER % (a, b))
    _check(_values(out), ref)


@pytest.mark.parametrize('case', list(MULT_CASES))
def test_mult_f1_matches_gs_k64(k64_mode, case):
    """MULT の F1・F2・F1K1F・F1K1B・FRINGE (四極の線形フリンジ tqlfre)・K0 つきの FB1・FB2 (tbfrme の線形フリンジ) の Emittance[] と CALC が gs と全桁一致"""
    a, b, ref = MULT_CASES[case]
    _f, out = _run(_MULT % (a, b))
    _check(_values(out), ref)


def test_steering_and_mult_f1_change_the_result_v2(v2_mode):
    """v2 でも走り、F1 の有無で結果が変わる (ステアリングは U0 が約 4.4 eV 小さくなる、MULT は運動量圧縮が 2e-3 変わる)"""
    _f, out0 = _run(_STEER % ('K0=0.0008', 'K0=-0.0008'))
    _f, out1 = _run(_STEER % STEER_CASES['v1'][:2])
    d0, d1 = _values(out0), _values(out1)
    assert d1['u0'] < d0['u0'] - 1.0
    _f, out0 = _run(_MULT % ('', ''))
    _f, out1 = _run(_MULT % MULT_CASES['m1'][:2])
    assert abs(_values(out1)['alpha'] / _values(out0)['alpha'] - 1.0) > 1e-4


def test_steering_calc_tunes_equal_emittance_tunes_without_rf_and_rad(v2_mode):
    """検算 (gs を使わない): NORFSW NORAD の Emittance[] のチューン (6D の一周写像) = CALC のチューン (4D)。F1・FB1・FB2・AE1・AE2 つきのステアリング"""
    src = (_STEER % STEER_CASES['v2'][:2]).replace('RFSW; RAD; NOFLUC;', 'NORFSW; NORAD;')
    _f, out = _run(src)
    d = _values(out)
    fc = d['calc'][1] % 1.0
    fe = (-d['tunes'][1]) % 1.0
    assert min(abs(fc - fe), abs(fc - (1 - fe))) < 2e-8, (d['calc'], d['tunes'])


# ---------------------------------------------------------------- B108: 独立な厳密な幾何 (円弧と直線)
_TH = 2 * np.pi / 24
_LB = 1.2
_RHO = _LB / _TH


def _pole(x, phi, d, L, K0, a1, a2):
    """硬い端の一様磁場の極 (参照軌道は直線 +Z、面は参照軌道に垂直な面から a1・a2 だけ傾く)。K0>0 は −X へ曲がる。K0<0 は鏡映 (x → −x、a → −a)。
    戻り (出口の参照面での x、向き φ、場の中の路程)"""
    if K0 < 0:
        x, phi, lf = _pole(-x, -phi, d, L, -K0, -a1, -a2)
        return -x, -phi, lf
    p = np.array([x, 0.0])
    dr = np.array([math.sin(phi), math.cos(phi)])
    f1 = np.array([math.cos(a1), math.sin(a1)])
    t, _s = np.linalg.solve(np.array([[dr[0], -f1[0]], [dr[1], -f1[1]]]), -p)
    p = p + t * dr
    P2 = np.array([0.0, L])
    f2 = np.array([math.cos(a2), -math.sin(a2)])
    n2 = np.array([-f2[1], f2[0]])
    rho = (L / K0) * (1 + d)
    c = p + rho * np.array([-dr[1], dr[0]])

    def arc(a):
        ca, sa = math.cos(a), math.sin(a)
        r = p - c
        return (c + np.array([ca * r[0] - sa * r[1], sa * r[0] + ca * r[1]]), np.array([ca * dr[0] - sa * dr[1], sa * dr[0] + ca * dr[1]]))
    lo, hi = 0.0, 1.2
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if np.dot(arc(mid)[0] - P2, n2) >= 0.0:
            hi = mid
        else:
            lo = mid
    a = 0.5 * (lo + hi)
    p, dr = arc(a)
    tt = (L - p[1]) / dr[1]
    p = p + tt * dr
    return p[0], math.atan2(dr[0], dr[1]), rho * a


def _bend(x, phi, d):
    """扇形の偏向磁石 (E1=E2=0) の厳密な幾何 (test_v2_A_damping.jx_geometric と同じ)"""
    p = np.array([x, 0.0])
    dr = np.array([math.sin(phi), math.cos(phi)])
    f1 = np.array([1.0, 0.0])
    t, _s = np.linalg.solve(np.array([[dr[0], -f1[0]], [dr[1], -f1[1]]]), -p)
    p = p + t * dr
    P2 = np.array([-_RHO + _RHO * math.cos(_TH), _RHO * math.sin(_TH)])
    f2 = np.array([math.cos(_TH), math.sin(_TH)])
    n2 = np.array([-f2[1], f2[0]])
    rho = _RHO * (1 + d)
    c = p + rho * np.array([-dr[1], dr[0]])

    def arc(a):
        ca, sa = math.cos(a), math.sin(a)
        r = p - c
        return (c + np.array([ca * r[0] - sa * r[1], sa * r[0] + ca * r[1]]), np.array([ca * dr[0] - sa * dr[1], sa * dr[0] + ca * dr[1]]))
    lo, hi = 0.0, 2 * _TH
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if np.dot(arc(mid)[0] - P2, n2) >= 0.0:
            hi = mid
        else:
            lo = mid
    a = 0.5 * (lo + hi)
    p, dr = arc(a)
    nout = np.array([math.cos(_TH), math.sin(_TH)])
    h = np.array([-math.sin(_TH), math.cos(_TH)])
    p = p - np.dot(p - P2, h) / np.dot(dr, h) * dr
    return np.dot(p - P2, nout), math.atan2(np.dot(dr, nout), np.dot(dr, h)), rho * a


def _quad(x, phi, d, L, K1):
    k = K1 / L / (1 + d)
    xp = math.tan(phi)
    if k > 0:
        w = math.sqrt(k)
        x, xp = math.cos(w * L) * x + math.sin(w * L) / w * xp, -w * math.sin(w * L) * x + math.cos(w * L) * xp
    else:
        w = math.sqrt(-k)
        x, xp = math.cosh(w * L) * x + math.sinh(w * L) / w * xp, w * math.sinh(w * L) * x + math.cosh(w * L) * xp
    return x, math.atan(xp)


def _jx_exact(K0, a1, a2, LP=0.2, DW=0.05):
    """24 セル (QF D2 B D2 QD DS [極 K0 DW 極 −2K0 DW 極 K0] DS D1) の Jx = 3 − Je。Je = 2 + dln(Σ B_i² l_i(δ))/dδ (U ∝ E² ∮B² dl。
    l_i は運動量 δ の閉軌道が磁場の中を通る路程。磁場の強さ B_i は固定)。閉軌道は 1 セルの写像の不動点"""
    wB = 1 / _RHO**2

    def cell(w, d):
        x, phi = _quad(w[0], w[1], d, 0.3, 0.7)
        x += 0.5 * math.tan(phi)
        x, phi, lb = _bend(x, phi, d)
        wsum = wB * lb
        x += 0.5 * math.tan(phi)
        x, phi = _quad(x, phi, d, 0.3, -0.6)
        x += 0.1 * math.tan(phi)
        for i, k in enumerate((K0, -2 * K0, K0)):
            x, phi, lp = _pole(x, phi, d, LP, k, a1, a2)
            wsum += (k / LP)**2 * lp
            if i < 2:
                x += DW * math.tan(phi)
        x += 1.1 * math.tan(phi)
        return np.array([x, phi]), wsum

    def wtot(d):
        w = np.zeros(2)
        for _ in range(30):
            r = cell(w, d)[0] - w
            J = np.zeros((2, 2))
            for i in range(2):
                dw = np.zeros(2)
                dw[i] = 1e-7
                J[:, i] = (cell(w + dw, d)[0] - cell(w - dw, d)[0]) / 2e-7
            w = w - np.linalg.solve(J - np.eye(2), r)
        return cell(w, d)[1]
    h = 1e-5
    w0 = wtot(0.0)
    return 3.0 - (2.0 + (wtot(h) - wtot(-h)) / (2 * h) / w0)


def _pole_ring(K0, ae, extra=''):
    return f'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1) DW = (L=0.05);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494);
BEND  P1 = (L=0.2 ANGLE=0 K0={K0} AE1={ae[0]} AE2={ae[1]});
BEND  P2 = (L=0.2 ANGLE=0 K0={-2 * K0} AE1={ae[0]} AE2={ae[1]});
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4E6 FREQ=5.0E8 PHI=0);
LINE CELL = (QF D2 B D2 QD DS P1 DW P2 DW P1 DS D1);
LINE RING = (ORIGIN RF 24*CELL);
FFS;
USE RING; CELL; CALCULATE; RFSW; RAD; NOFLUC; RADCOD;
{extra}
e=Emittance[]; Print[DampingRate/.e]; Print[EnergyLossU0/.e];
'''


def _emit_jx(K0, ae):
    _f, out = _run(_pole_ring(K0, ae, 'EMITDIV=20;'))
    o = [l for l in out.split(chr(10)) if l.strip()]
    d = [float(x) for x in o[-2].strip('{}').split(',')]
    k = float(o[-1]) / 6e9
    return -d[0] / k


@pytest.mark.parametrize('ae', [(0.3, 0.3), (-0.2, 0.4)])
def test_steering_edge_angle_term_v2_matches_exact_geometry_and_k64_doubles_it(ae):
    """B108。極の端の角 (AE1・AE2) による Jx の変化 ΔJx (AE なしとの差) が、厳密な幾何 (円弧と直線のロビンソン) と一致するのは v2 (行路を x·tan(e) だけ短くする)。
    K64 は端の項が約 2 倍 (比 2.0〜2.2)。K0 = 0.005 の極 (軌道の傾き 5 mrad) で比べる。残る差は 3〜4% (有限の分割 EMITDIV=20 と 2 次の項)"""
    K0 = 0.005
    ex = _jx_exact(K0, *ae) - _jx_exact(K0, 0.0, 0.0)
    out = {}
    for mode in ('v2', 'k64'):
        old = compat._mode
        compat.set_mode(mode)
        try:
            out[mode] = _emit_jx(K0, ae) - _emit_jx(K0, (0.0, 0.0))
        finally:
            compat.set_mode(old)
    assert abs(ex) > 3e-3
    assert abs(out['v2'] / ex - 1.0) < 0.06, (ex, out)
    assert 1.8 < out['k64'] / ex < 2.4, (ex, out)


def test_steering_edge_angle_wedge_cancels_to_first_order(v2_mode):
    """AE1 = −AE2 (くさび): 行路の 1 次の変化 (x·(tan e1 + tan e2)) が打ち消す。v2 の ΔJx は小さく (K0 = 0.005 で 1e-3 以下)、厳密な幾何と 1e-3 で一致"""
    K0 = 0.005
    ex = _jx_exact(K0, 0.3, -0.3) - _jx_exact(K0, 0.0, 0.0)
    v2 = _emit_jx(K0, (0.3, -0.3)) - _emit_jx(K0, (0.0, 0.0))
    assert abs(ex) < 2e-3 and abs(v2 - ex) < 1e-3, (ex, v2)
