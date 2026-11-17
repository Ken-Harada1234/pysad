# -*- coding: utf-8 -*-
"""担当 A1 (v2next): N12 (INS の Twiss["DETR"])、TrackParticles の MASS・CHARGE (陽子)、ほか。
v2next 項目の試験。根拠は解析式 (相対論の運動学・ラーマー運動・古典放射・シンクロトロン振動)。"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat

E_MASS = 0.51099895e6
P_MASS = 938.272088e6          # 陽子の質量 [eV/c^2]
C_LIGHT = 299792458.0


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _floats(s):
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', s)]


def _tag(text, tag):
    """text の中の `tag ` で始まる出力 (複数行にまたがる) の数値列"""
    i = text.index('\n' + tag + ' ') + 1 if ('\n' + tag + ' ') in text else text.index(tag + ' ')
    j = text.find('\n', i)
    seg = text[i + len(tag):]
    # 次の行頭が英字で始まる別のタグ (または終端) まで
    m = re.search(r'\n(?=[A-Za-z?])', seg)
    seg = seg[:m.start()] if m else seg
    return _floats(seg.replace('{', ' ').replace('}', ' ').replace(',', ' '))


# ---------------------------------------------------------------- N12: INS の Twiss["DETR"]
_RING_C = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QS = (L=0.2 K1=0.05 ROTATE=0.7853981633974483);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = (); MARK M2 = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 11*CELL1 QS M2 CELL1);
FFS;
USE RING; CELL; CALC;
"""


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_n12_ins_calc_keeps_detr(mode):
    """INS の CALC の後の Twiss["DETR"] が、R 行列の行列式 r1 r4 - r2 r3 (DETR < 1 のとき) になる (pysad の不具合: 0 になっていた。両モード同じ)。
    gs (K64) の値 (結合したリング、CELL CALC の後の INS CALC の M2 での値) は 0.0021332033027 と -0.0800626854892 (R1)"""
    compat.set_mode(mode)
    try:
        src = _RING_C.replace('K1=0.05', 'K1=0.8') + 'INS; CALC;\n' + \
            'Print["RES ", Twiss["DETR","M2"], " ", Twiss["R1","M2"], " ", Twiss["R2","M2"], " ", Twiss["R3","M2"], " ", Twiss["R4","M2"]];\n'
        v = _tag(_out(src), 'RES')
        detr, r1, r2, r3, r4 = v
        assert abs(detr - (r1 * r4 - r2 * r3)) < 1e-12
        assert detr != 0.0
        assert abs(detr - 0.0021332033027) < 1e-10 and abs(r1 + 0.0800626854892) < 1e-10      # gs と一致
    finally:
        compat.set_mode(None)


# ---------------------------------------------------------------- TrackParticles の MASS・CHARGE
def _beta(p, m):
    return p / math.hypot(p, m)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_track_proton_drift_z_exact(mode):
    """陽子 (MASS=938.272088e6) の 10 m のドリフト: z = L (β1/β0 - 1) (z は δ と正準共役な β c (t0 - t) 型の座標)。
    px だけの粒子は z = -L px^2/2。MASS を渡さない (電子の質量のまま) 旧動作では 1e-3 の δ で 3e-10 しか出ず、一致しない"""
    compat.set_mode(mode)
    try:
        src = f"""MASS = {P_MASS!r}; CHARGE = 1; MOMENTUM = 3.0e9;
DRIFT D1 = (L=10); MARK M0 = (BX=1 BY=1); LINE TL = (M0 D1);
FFS;
USE TL; INS; CALC; NORAD;
bm={{{{1e-3,0,1e-3}},{{2e-3,0,0}},{{1e-3,0,-1e-3}},{{0,0,0}},{{0,0,0}},{{0,1e-3,-1e-2}},{{1,1,1}}}};
$FORM="16.12";
r=TrackParticles[{{1,bm}},1,1,1];
Print["RZ ", r[[2,5]]];
"""
        z = _tag(_out(src), 'RZ')
        p0 = 3.0e9
        for dz, d in ((z[1], 1e-3), (z[2], -1e-2)):
            ex = 10.0 * (_beta(p0 * (1 + d), P_MASS) / _beta(p0, P_MASS) - 1.0)
            assert abs(dz - ex) < 1e-12 * max(1.0, abs(ex) / 1e-3), (dz, ex)
        assert abs(z[0] + 10.0 * (2e-3) ** 2 / 2) < 1e-9        # 2 次 (4 次の項 6e-11 は除く)
    finally:
        compat.set_mode(None)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('charge', [1.0, -1.0, 2.0])
def test_track_solenoid_larmor_charge(mode, charge):
    """ソレノイド (BOUND=1 の SOL 領域、L=2 m、BZ=1 T、p=3 GeV/c) の通過: 厚いソレノイドの行列 (k = charge BZ c / (2 p)、C = cos kL、S = sin kL)
    x_out = C^2 x0、y_out = -S C x0 (入口で x0、運動量 0)。電荷の符号で回転の向きが変わり、電荷 2 は k が 2 倍。磁気剛性は p/(charge c)"""
    compat.set_mode(mode)
    try:
        src = f"""MASS = {P_MASS!r}; MOMENTUM = 3.0e9;
SOL SA = (L=0 BZ=1.0 BOUND=1) SB = (L=0 BZ=0 BOUND=1); DRIFT DD = (L=2); MARK M0 = (BX=1 BY=1);
LINE TL = (M0 SA DD SB);
FFS;
USE TL; INS; NORAD; CHARGE = {charge};
bm={{{{1e-3}},{{0}},{{0}},{{0}},{{0}},{{0}},{{1}}}};
$FORM="16.12";
r=TrackParticles[{{1,bm}},1,1,1];
Print["RX ", r[[2,1]], " ", r[[2,3]]];
"""
        x, y = _tag(_out(src), 'RX')
        k = charge * 1.0 * C_LIGHT / (2 * 3.0e9)
        C, S = math.cos(2 * k), math.sin(2 * k)
        assert abs(x - C * C * 1e-3) < 5e-12
        assert abs(y + S * C * 1e-3) < 5e-12
    finally:
        compat.set_mode(None)


_PRING = """MASS = {mass!r}; MOMENTUM = 30.0e9;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT={volt!r} HARM=24 PHI=0);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 {rf});
FFS;
USE RING; CELL; CALC; RFSW; NORAD;
"""


def _ring(volt, rf='RFC', mass=P_MASS):
    return _PRING.format(mass=mass, volt=volt, rf=rf)


def _oneturn_zd(src_head, eps_z=1e-7, eps_d=1e-9):
    """1 周写像の (z, δ) 部分 (2x2)。z=±eps_z、δ=±eps_d (分散軌道上: x = EX δ、px = EPX δ) の 4 粒子を 1 周追跡して中心差分"""
    src = src_head + f"""ex=Twiss["EX",1]; epx=Twiss["EPX",1];
bm={{{{0,0,ex*{eps_d},-ex*{eps_d}}},{{0,0,epx*{eps_d},-epx*{eps_d}}},{{0,0,0,0}},{{0,0,0,0}},{{{eps_z},-{eps_z},0,0}},{{0,0,{eps_d},-{eps_d}}},{{1,1,1,1}}}};
$FORM="20.14";
r=TrackParticles[{{1,bm}},1,1,1];
Print["OT ", r[[2,5]], " ", r[[2,6]]];
"""
    v = _tag(_out(src), 'OT')
    z, d = v[:4], v[4:]
    return np.array([[(z[0] - z[1]) / (2 * eps_z), (z[2] - z[3]) / (2 * eps_d)],
                     [(d[0] - d[1]) / (2 * eps_z), (d[2] - d[3]) / (2 * eps_d)]])


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_proton_rf_kick_and_synchrotron_tune(mode):
    """陽子 (30 GeV/c、γ = 32、周長 45.6 m、h = 24、V = 2 MV、PHI = 0) の RF と縦方向の 1 周写像 (TrackParticles の 1 周を中心差分):
    - z → δ の傾き = 2π h V / (C β^2 E)  (E = sqrt(p^2 + m^2)、ΔE = eV sin φ に対し Δδ = ΔE / (β p c))。電子の質量のままだと E が同じでも β が違う
    - スリップ因子 η は RF なしの 1 周の z 変化から (z → z - C η δ)、Qs = sqrt(h η V / (2π β^2 E)) が 1 周写像のトレース (2 cos 2πQs) と 1e-5 で一致"""
    compat.set_mode(mode)
    try:
        C = 45.6
        p = 30.0e9
        E = math.hypot(p, P_MASS)
        b2 = (p / E) ** 2
        h, V = 24.0, 2.0e6
        M0 = _oneturn_zd(_ring(0.0, rf='ST'))
        M = _oneturn_zd(_ring(V))
        eta = -M0[0, 1] / C
        assert abs(M[1, 0] / (2 * math.pi * h * V / (C * b2 * E)) - 1.0) < 1e-3
        qs_an = math.sqrt(h * eta * V / (2 * math.pi * b2 * E))
        qs_tr = math.acos(0.5 * np.trace(M)) / (2 * math.pi)
        assert abs(qs_tr / qs_an - 1.0) < 1e-4, (qs_tr, qs_an)
        assert abs(np.linalg.det(M) - 1.0) < 1e-9
    finally:
        compat.set_mode(None)


def _u0_analytic(mass, charge=1.0):
    """24 偏向 (12 個、θ = 0.2617993877991494、L = 1.2 m) のリングの 1 周の放射損失 [eV]: U0 = (2/3) r_c mc^2 β^3 γ^4 Σ θ^2/L (r_c = 電荷^2 e^2/(4π ε0 m c^2))"""
    from pysad.emit6d import ELRADI, ELMASS
    rc = charge ** 2 * ELRADI * ELMASS / mass
    p = 30.0e9
    E = math.hypot(p, mass)
    g = E / mass
    b = p / E
    return 2.0 / 3.0 * rc * mass * b ** 3 * g ** 4 * 12 * 0.2617993877991494 ** 2 / 1.2


def _turn_delta(mass, volt=0.0):
    src = _ring(volt, rf='ST', mass=mass) + """RAD; NOFLUC;
bm={{0},{0},{0},{0},{0},{0},{1}};
$FORM="16.8";
r=TrackParticles[{1,bm},1,1,1];
Print["DD ", r[[2,6,1]]];
"""
    return _tag(_out(src), 'DD')[0]


def test_proton_radiation_loss_v2(v2_mode):
    """v2: 陽子の放射損失が古典式 U0 = (2/3) r_p mc^2 β^3 γ^4 Σ θ^2/L (r_p = r_e m_e / m_p) と一致する。
    TrackParticles の 1 周の δ の変化 (RF なし) は -U0 / (β^3 p c) (相対 1e-6): K64 の tradkn は光子数 (5/(2√3)) α γ θ と臨界エネルギーから
    δ を Δδ = -(エネルギー) / (p c) だけ下げる超相対論的な式で、β の因子 (β^3) を持たない。β = 1 - 1/(2γ^2) なので γ ≫ 1 (放射が意味を持つ範囲) では
    U0 との差は 1/γ^2 (γ = 32 で 1e-3)。Emittance[] の EnergyLossU0 も同じ式で U0/β^3 と 1e-6 で一致"""
    for mass in (P_MASS, 2.0 * P_MASS):
        u0 = _u0_analytic(mass)
        p = 30.0e9
        b = p / math.hypot(p, mass)
        dd = _turn_delta(mass)
        assert abs(dd * p * b ** 3 / (-u0) - 1.0) < 1e-6, (mass, dd, u0)
    src = _ring(2e6) + "e=Emittance[];" + chr(10) + 'Print["U0 ", EnergyLossU0/.e];' + chr(10)
    b = 30.0e9 / math.hypot(30.0e9, P_MASS)
    assert abs(_tag(_out(src), 'U0')[0] * b ** 3 / _u0_analytic(P_MASS) - 1.0) < 1e-6      # Emittance[] も同じ式 (β^3 の因子なし)


def test_proton_radiation_loss_k64(k64_mode):
    """k64: K64 (tfsetparam.f: re0 = elradi) は質量によらず電子の古典半径を使うので、陽子の放射損失は mp/me 倍 (1836 倍) 大きい (gs で確認済み)"""
    u0 = _u0_analytic(P_MASS)
    p = 30.0e9
    b = p / math.hypot(p, P_MASS)
    dd = _turn_delta(P_MASS)
    assert abs(dd * p * b ** 3 / (-u0) / (P_MASS / E_MASS) - 1.0) < 1e-6


def test_electron_radiation_same_in_both_modes():
    """電子 (MASS 既定) の放射は v2 でも k64 と同じ (ビット単位)"""
    out = {}
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        try:
            src = _ring(2e6, mass=E_MASS) + "RAD; NOFLUC;\nbm={{1e-3},{0},{0},{0},{0},{1e-3},{1}};\n$FORM=\"20.14\";\nr=TrackParticles[{1,bm},1,1,5];\nPrint[\"EL \", r[[2]]];\n"
            out[m] = _out(src)
        finally:
            compat.set_mode(None)
    assert out['k64'] == out['v2']


# ---------------------------------------------------------------- emit6d._tquade の解析ヤコビアン (軌道がずれたとき)
def _numjac(f, c, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        J[:, j] = (f(cp) - f(cm)) / (2 * h)
    return J


def test_quad_jacobian_full_v2(v2_mode):
    """v2: _jac_quad_full (四極 1 個の解析ヤコビアン) が、x・y・δ がずれた軌道で数値微分 (刻み 1e-7) と 6x6 全成分で 1e-8 以内に一致する。
    K64 の動作 (k64 モード) では z の行 (列 x px y py) が 0 のままで、emit6d.State.piece がシンプレクティック条件から組み直していた (結果は同じ)"""
    from pysad import emit6d as e
    ph = e._Phys(3e9, 100.0)
    rng = np.random.default_rng(5)
    for K in (0.7, -0.6, 2.5):
        for _ in range(4):
            c = rng.normal(size=6) * 1e-3
            c[4] = 3e-4
            c[5] = rng.normal() * 2e-3
            Jn = _numjac(e._quad_map(0.3, K, ph), c)
            Ja = e._jac_quad_full(c, 0.3, K, True, (0, 0, 0, 0, 0), ph)
            assert np.abs(Ja - Jn).max() < 2e-7        # 数値微分の丸め (強い K=2.5 で 3e-8)


def test_quad_jacobian_k64_keeps_zero_z_row(k64_mode):
    from pysad import emit6d as e
    ph = e._Phys(3e9, 100.0)
    c = np.array([1e-3, 2e-4, -5e-4, 1e-4, 3e-4, 1e-3])
    Ja = e._jac_quad_full(c, 0.3, 0.7, True, (0, 0, 0, 0, 0), ph)
    assert np.all(Ja[4, :4] == 0.0)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_tquade_T_matches_numeric_with_orbit(mode):
    """Python 版 _tquade と numba 版 (elkern) の 1 周の T が、DX・DY・ROTATION のある四極と軌道のずれで 1e-9 以内に一致し、
    単独の四極の T は写像の数値微分と 2e-7 以内 (z 行は組み直しのシンプレクティック条件)"""
    from pysad import emit6d as e
    from pysad.elements import Quadrupole
    compat.set_mode(mode)
    try:
        ph = e._Phys(3e9, 100.0)
        rng = np.random.default_rng(1)
        for _ in range(24):
            K = rng.choice([0.7, -0.6, 2.5, -3.0])
            L = rng.choice([0.3, 0.6, 1.5])
            dx, dy = rng.normal(size=2) * 2e-3
            th = rng.choice([0.0, 0.3, -0.7])
            cod = rng.normal(size=6) * rng.choice([1e-4, 1e-3, 5e-3])
            cod[4] = rng.normal() * 1e-3
            cod[5] = rng.normal() * 2e-3
            q = Quadrupole(name='Q', L=float(L), K1=float(K))
            q.DX, q.DY, q.ROTATION = float(dx), float(dy), float(th)
            s_py = e._State(cod.copy(), ph)
            e._tquade(s_py, q, True)
            s_fast = e._State(cod.copy(), ph)
            assert e._quade_fast(s_fast, q, True) == 0
            assert np.abs(s_py.T - s_fast.T).max() < 1e-9
            Jn = _numjac(e._quad_map(float(L), float(K), ph, float(dx), float(dy), float(th)), cod)
            assert np.abs(s_py.T - Jn).max() < 2e-7
    finally:
        compat.set_mode(None)


# ---------------------------------------------------------------- BEND の F1・FRINGE + RAD の追跡 (tfgeo の iprev・inext)
_BEND_F1_LINE = """MOMENTUM = 3GEV;
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747 F1=0.1 FRINGE=3);
MARK M0 = (BX=1 BY=1);
LINE TL = (M0 B);
FFS;
USE TL; INS; CALC; NORAD;
$FORM="18.12";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
RAD; NOFLUC;
r=TrackParticles[{1,bm},1,1,1]; Print["T1 ",r[[2]]];
"""
# gs (K64) の出力 (local/wake/gsrun.sh)。1 個しかない BEND の輸送路は、先頭の MARK が長さ 0 なので tfgeo の inext・iprev が自分自身に繋がり
# (iprev = inext = 自分)、F1 の端の分割 (n = -1, 0, ndiv+1, ndiv+2) に放射を入れない
_GS_LINE_X = [1.092318870e-3, 1.548953264374e-4, -6.04233811313e-4]
_GS_LINE_DP = [-2.17142054679e-5, 9.782447092820e-4, -1.0021281980e-2]

_BEND_F1_RING = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=0.6 ANGLE=0.1308996938995747 E1=0.0654498469497874 E2=0 F1=0.1 FRINGE=3);
BEND B2 = (L=0.6 ANGLE=0.1308996938995747 E1=0 E2=0.0654498469497874 F1=0.1 FRINGE=3);
MARK MK = ();
CAVI RFC = (L=0 VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B B2 D2 QD D1);
LINE CELL2 = (QF D2 B MK B2 D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
LINE RINGM = (ST 12*CELL2 RFC);
FFS;
$FORM="14.8";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; r=TrackParticles[{1,bm},1,1,1]; Print["T1 ",r[[2]]];
USE RINGM; CELL; CALC; RFSW; RAD; NOFLUC; r=TrackParticles[{1,bm},1,1,1]; Print["T2 ",r[[2]]];
"""
# gs の出力 (隣り合う BEND の対 B B2 の 12 セル): 行 x, px, y, py, z, dp (粒子 3 個ずつ)
_GS_RING = [[0.00138678, 2.98238682e-4, -0.00143511], [8.36993148e-4, -5.4351648e-5, 0.00138246],
            [-9.9254076e-4, -1.10806e-20, 7.39808960e-4], [0.00105006, -6.686232e-21, -9.7279926e-4],
            [5.78360966e-4, -0.00138341, 0.01585603], [-2.5073905e-4, 7.40115642e-4, -0.0101792]]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_bend_f1_single_bend_line_matches_gs(mode):
    """F1=0.1 FRINGE=3 の BEND 1 個の輸送路 (INS) の RAD の追跡が gs と一致 (1e-7)。1 個しかない BEND は自分自身に繋がる (iprev = inext = 自分) ので
    F1 の端の断片に放射を入れない (旧 pysad は入れて、δ の放射損失が 3% 違った)。両モード同じ (K64 の動作)"""
    compat.set_mode(mode)
    try:
        tol = 1e-7 if mode == "k64" else 1e-4       # v2 は BEND の端の角の場の中の行路 (tan ψ x の補正) が入り、gs と 1e-7 ほど違う
        v = _tag(_out(_BEND_F1_LINE), 'T1')
        x, dp = v[0:3], v[15:18]
        for a, b in zip(x + dp, _GS_LINE_X + _GS_LINE_DP):
            assert abs(a - b) <= tol * max(abs(b), 1e-5), (a, b)
    finally:
        compat.set_mode(None)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_bend_f1_adjacent_bends_match_gs(mode):
    """F1=0.1 FRINGE=3 の隣り合う BEND (B と B2 の間が距離 0 か長さ 0 の MARK だけ) の RAD の追跡が gs と一致 (3e-6。8 桁表示の丸め)。
    隣り合う側では F1 の断片 (放射) を入れない (tfgeo の inext・iprev)。旧 pysad は入れて、px で相対 3.6% 違った"""
    compat.set_mode(mode)
    try:
        out = _out(_BEND_F1_RING)
        for tag in ('T1', 'T2'):
            v = _tag(out, tag)
            rows = [v[3 * k:3 * k + 3] for k in range(6)]
            for r, g in zip(rows, _GS_RING):
                for a, b in zip(r, g):
                    assert abs(a - b) <= (3e-6 if mode == 'k64' else 1e-4) * max(abs(b), 1e-5) + 2e-11, (tag, a, b)
    finally:
        compat.set_mode(None)


def test_bend_links_function():
    """tfgeo.f の inext・iprev: 隣り合う BEND (間は長さ 0 の要素だけ) は繋がり、F1 の断片は繋がっていない側にだけ入る"""
    from pysad.track6n import bend_links
    from pysad import parser
    src = """BEND B = (L=0.5 ANGLE=0.1); MARK MK = (); DRIFT D = (L=1);
LINE R1 = (D B B D); LINE R2 = (D B MK B D); LINE R3 = (D B D B D); LINE R4 = (MK B); LINE R5 = (D B);
"""
    sc = parser.load_sad_string(src)
    assert bend_links(sc.build_lattice('R1').elements) == {1: (True, False), 2: (False, True)}
    assert bend_links(sc.build_lattice('R2').elements) == {1: (True, False), 3: (False, True)}
    assert bend_links(sc.build_lattice('R3').elements) == {1: (True, True), 3: (True, True)}
    assert bend_links(sc.build_lattice('R4').elements) == {1: (False, False)}      # 最後の BEND が終わりから先頭へ回り、長さ 0 の MARK を越えて自分に繋がる
    assert bend_links(sc.build_lattice('R5').elements) == {1: (True, True)}      # 先頭が長さのある DRIFT なら繋がらない


# ---------------------------------------------------------------- CAVI の DX・DY・ROTATE と V1・V20・V11・V02 (tcav.f の TENT・TEXIT と v1a…v02a)
_CAV_RING = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (@CAV@);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; NORAD;
$FORM="14.8";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
RAD; NOFLUC;
r=TrackParticles[{1,bm},1,1,1]; Print["T1 ",r[[2]]];
"""
# gs (K64) の 1 周後 (RAD): 行 x, px, y, py, z, dp (粒子 3 個ずつ)
_GS_CAV_DXDY = [0.00152394, 0.0004501705, -0.00269744, 0.00078021027, 3.25916411e-06, 0.000855858945, 4.89197643e-05, 6.12812188e-09,
                -0.00026872608, 0.000823699236, -1.64157e-13, -0.00069224483, 0.00061989206, -0.00140615, 0.01626493,
                -0.00025781297, 0.000732851343, -0.01018546]
_GS_CAV_V = [0.001523, 0.000449237006, -0.00269834, 0.000776436144, -4.7754473e-07, 0.000852305637, 4.91167645e-05, 2.05327786e-07,
             -0.00026852043, 0.000824487035, 7.97381208e-07, -0.00069143058, 0.000619892633, -0.00140615, 0.01626493,
             -0.00025781298, 0.000732851734, -0.01018548]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('cav,gs', [
    ('L=0.5 VOLT=1.5e6 HARM=64 DX=1e-3 DY=-2e-3', _GS_CAV_DXDY),
    ('L=0.5 VOLT=1.5e6 HARM=64 V1=1e5 V11=0.5e6 V20=1e6 V02=-3e5 DX=1e-3 DY=-2e-3 ROTATE=0.2', _GS_CAV_V)])
def test_cavity_offset_and_transverse_voltage_match_gs(mode, cav, gs):
    """L=0.5 の CAVI に DX・DY (・ROTATE) と V1・V20・V11・V02 を付けた RAD の追跡が gs と一致する (表示 9 桁。v2 は BEND の端の角の補正で 1e-4 まで)。
    旧 pysad は CAVI の DX・DY・ROTATE を無視し (x=y=0 の粒子に入る縦の蹴り y 6e-9 が無く、y が 1e-3 違った)、
    V1・V20・V11・V02 は読まず (パーサが捨てていた)、k_cav に式が無かった"""
    compat.set_mode(mode)
    try:
        v = _tag(_out(_CAV_RING.replace('@CAV@', cav)), 'T1')[:18]
        tol = 3e-8 if mode == 'k64' else 2e-4
        for a, b in zip(v, gs):
            assert abs(a - b) <= tol * max(abs(b), 1e-5) + (1e-13 if mode == "k64" else 5e-9), (a, b)
    finally:
        compat.set_mode(None)


def test_cavity_offset_changes_orbit():
    """V11 のある CAVI の DX=1e-3 は、y = 0 の粒子に入る縦の蹴り (V11 × (x − DX)) を変える (旧 pysad は DX を無視して y が同じだった)"""
    base = _out(_CAV_RING.replace('@CAV@', 'L=0.5 VOLT=1.5e6 HARM=64 V11=0.5e6'))
    off = _out(_CAV_RING.replace('@CAV@', 'L=0.5 VOLT=1.5e6 HARM=64 V11=0.5e6 DX=1e-3'))
    y_base = _tag(base, 'T1')[6:9]
    y_off = _tag(off, 'T1')[6:9]
    assert abs(y_off[1] - y_base[1]) > 1e-9      # 4.8e-9 (DX = 0: -2.1e-9、DX = 1e-3: +2.6e-9)


# ---------------------------------------------------------------- EMIT のあとの USE (trf0 が残る)
_USE_AFTER_EMIT = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=1.5e6 HARM=64 PHI=0.1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="14.8";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
USE RING; CELL; CALC; RFSW; RAD; NOFLUC;
r=TrackParticles[{1,bm},1,1,1]; Print["A1 ",r[[2]]];
@EM@
r=TrackParticles[{1,bm},1,1,1]; Print["B1 ",r[[2]]];
USE RING;
r=TrackParticles[{1,bm},1,1,1]; Print["C1 ",r[[2]]];
CELL; CALC;
r=TrackParticles[{1,bm},1,1,1]; Print["D1 ",r[[2]]];
"""
_GS_USE_B1_DP = [-4.1421124e-05, 9.50018237e-4, -0.00997733]       # EMIT 後 (USE の後も同じ)
_GS_USE_A1_DP = [-3.077161e-4, 6.82902209e-4, -0.01023444]          # EMIT の前 (trf0 = 0)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('em', ['EMIT;', 'e=Emittance[];'])
def test_use_after_emit_keeps_trf0(mode, em):
    """EMIT (Emittance[]) のあとに USE しても、追跡の RF の位相 (trf0 = 同期位置の時間) は EMIT の値のまま (K64 の trf0 は大域変数で、USE では
    0 に戻らず、temit と CALC6D だけが更新する。gs で EMIT・Emittance[]・CALC6D のどれでも確認)。旧 pysad は USE で trf0 = 0 に戻した。
    EMIT 前 (trf0 = 0) と後の δ は、gs の値 (表示 9 桁) と一致する (k64 は 1e-8。v2 は BEND の端の角の補正で 1e-3 まで)"""
    compat.set_mode(mode)
    try:
        out = _out(_USE_AFTER_EMIT.replace('@EM@', em))
        a1, b1, c1, d1 = (_tag(out, t) for t in ('A1', 'B1', 'C1', 'D1'))
        tol = 1e-7 if mode == 'k64' else 2e-3
        for a, b in zip(b1[15:18], _GS_USE_B1_DP):
            assert abs(a - b) <= tol * max(abs(b), 1e-4)
        for a, b in zip(a1[15:18], _GS_USE_A1_DP):
            assert abs(a - b) <= tol * max(abs(b), 1e-4)
        assert max(abs(x - y) for x, y in zip(c1, b1)) == 0.0 and max(abs(x - y) for x, y in zip(d1, b1)) == 0.0
        assert max(abs(x - y) for x, y in zip(a1, b1)) > 1e-4         # EMIT で trf0 が変わった
    finally:
        compat.set_mode(None)


# ---------------------------------------------------------------- RADTAPER・RADCOD の閉軌道と追跡 (tradke の ddpx の符号)
_TAPER_RING = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SD = (L=0.2 K2=0) SF = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
$FORM = "20.14";
USE RINGQ; CELL; CALC;
RFSW; RAD; RADCOD; NOFLUC;
@TAPER@
e0 = Emittance[];
c0 = OrbitAtExit/.e0; c0[[5]]=0;
bm = Append[Table[{c0[[i]]},{i,1,6}],{1}];
r = TrackParticles[{1,bm},1,1,1];
Print["D1 ", Flatten[Take[r[[2]],6]]-c0];
Print["OB ", OrbitAtExit/.e0, " ", EnergyLossU0/.e0];
"""


@pytest.mark.parametrize('taper', ['', 'RADTAPER;'])
def test_radcod_closed_orbit_closes_in_tracking_v2(taper, v2_mode):
    """v2: Emittance[] (RFSW・RAD・RADCOD、RADTAPER あり・なし) の閉軌道から追跡 (RAD、NOFLUC) を 1 周始めると、閉じる (x・px・z の動きが 1e-8 以下)。
    K64 は tradke (Emittance[]・閉軌道の放射) の ddpx = +0.5 dpx dg と、tradkn (追跡) の ddpx = −0.5 dpx dg で符号が逆で、閉軌道が追跡で閉じない
    (3 GeV、周長 105.6 m、4 MV のリングで x 4e-7 m、z 2e-6 m (RADTAPER なし)。RADTAPER あり: 2e-6 m・1e-5 m)。
    エネルギーを失うと同じ磁場で偏向が強くなる (平均の運動量が違う) ので追跡の符号が正しい。RADTAPER の追跡の BEND の係数も tturne と同じ (3δ_l + δ_l+1)/4 にする"""
    out = _out(_TAPER_RING.replace('@TAPER@', taper))
    d1 = _tag(out, 'D1')
    assert abs(d1[0]) < 1e-8 and abs(d1[1]) < 1e-8 and abs(d1[4]) < 1e-7 and abs(d1[5]) < 1e-8, d1


@pytest.mark.parametrize('taper', ['', 'RADTAPER;'])
def test_radcod_closed_orbit_does_not_close_k64(taper, k64_mode):
    """k64: K64 のとおり、Emittance[] の閉軌道が追跡で 4e-7 m 以上動く (gs も同じ。x −4.349e-7 m、z 2.375e-6 m、RADTAPER なし)"""
    d1 = _tag(_out(_TAPER_RING.replace('@TAPER@', taper)), 'D1')
    assert abs(d1[0]) > 3e-7 and abs(d1[4]) > 1e-6, d1


def test_tradke_python_and_numba_agree_v2(v2_mode):
    """v2: tradke の ddpx の符号反転が、numba カーネル (elkern) と Python 版 (emit6d._State.tradke) の両方に入り、閉軌道が 1e-12 で一致する"""
    from pysad import emit6d as e, parser
    src = _TAPER_RING.split('FFS;')[0]
    lat = parser.load_sad_string(src).build_lattice('RINGQ')
    ph = e._Phys(3e9, e._emit_circ(lat.elements))
    ph.dvfs = 0.0
    cods = []
    for fast in (True, False):
        old = e._USE_ELKERN
        e._USE_ELKERN = fast
        try:
            cod, trf0, fnd, _ = e.closed_orbit6d(lat.elements, ph, True, True)
        finally:
            e._USE_ELKERN = old
        assert fnd
        cods.append(np.array(cod))
    assert np.abs(cods[0] - cods[1]).max() < 1e-12


# ---------------------------------------------------------------- PSPAC のメッシュの分解能 (gs との差の原因: メッシュ)
_PSPAC_LINE = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.5);
SPCH SP1 = ();
MARK ORIGIN = ();
LINE TL = (ORIGIN D1 SP1 D1 SP1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
@MESH@
a = 1e-3; lb = @LB@; n = 4096;
u = Table[Mod[i*0.6180339887, 1], {i,1,n}];
v = Table[Mod[i*0.7548776662, 1], {i,1,n}];
w = Table[Mod[i*0.5698402910, 1], {i,1,n}];
bm = {a*Sqrt[u]*Cos[2*Pi*v], Table[0,{n}], a*Sqrt[u]*Sin[2*Pi*v], Table[0,{n}], lb*(w-0.5), Table[0,{n}], Table[1,{n}]};
rms[l_] := Sqrt[Plus@@(l^2)/Length[l]];
PSPAC;
r = TrackParticles[{1,bm},1,1,1];
Print["OUT ", rms[r[[2,1]]], " ", rms[r[[2,3]]]];
"""


def _pspac_rms(mesh, lb):
    out = _out(_PSPAC_LINE.replace('@MESH@', mesh).replace('@LB@', repr(lb)))
    return _tag(out, 'OUT'), out


def _pspac_analytic(lb, N=1e11, a=1e-3, p=5e6):
    """薄いキックの解析式: 一様な円のビームのキック px = K L x/a^2 (K = 2 r_e N/(lb β^2 γ^3))、SPCH の間 0.5 m のドリフトで x が (1 + s/2) 倍"""
    RE = 2.8179403205e-15
    pp = p / E_MASS
    g = math.sqrt(1 + pp * pp)
    b = pp / g
    s = 2 * RE * N / lb / (b * b * g ** 3) * 0.5 / a ** 2
    return 0.5e-3 * (1 + 0.5 * s)


def test_pspac_resolution_vs_analytic_and_gs(v2_mode):
    """PSPAC の gs との差 (v2next 項目 8) はメッシュの分解能。KV 風の一様な円 (半径 1 mm)・一様な長さ lb のバンチ、5 MeV、1e11 個、D(0.5) SPCH D(0.5) SPCH。
    分解能が足りているとき (a = 10 セル、lb = 20 セル) は v2 が解析式と 1.7%、gs (K64) と 0.6% で合い、x と y が 0.1% で等しい。
    (gs: x 1.2338、y 1.2356 mm、v2: 1.2270、1.2282 mm、解析式 1.2481 mm、4096 個。1/γ (5〜50 MeV)・粒子数 (500〜4096 個) には 0.3% しか依らない)
    縦のメッシュが粗いと v2 の λ(z) (2 次の B スプラインでぼかす) が gs (K64 の rhocal はほぼ最近接) より小さく出て、lb が 4 セルで v2 は解析式の −12% (gs −8%)、
    2 セルで −27% (gs −17%)、1 セルで −48% (gs −20%)。警告 (???-PSPAC-Warning) を出す"""
    mesh = "PSPACNX = 128; PSPACNY = 128; PSPACDX = 1e-4; PSPACDY = 1e-4; PSPACNZ = 32; PSPACDZ = 0.005;"
    (x, y), out = _pspac_rms(mesh, 0.1)
    assert abs(x / _pspac_analytic(0.1) - 1) < 0.02 and abs(y / x - 1) < 1e-2
    assert abs(x - 1.2338e-3) < 0.01 * 1.2338e-3                    # gs
    assert 'Warning' not in out
    for lb, lo, hi in ((0.02, -0.15, -0.08), (0.01, -0.32, -0.22)):
        (x, y), out = _pspac_rms(mesh, lb)
        r = x / _pspac_analytic(lb) - 1
        assert lo < r < hi, (lb, r)
    assert 'longitudinal mesh is too coarse' in out
    # 縦のメッシュを自動 (PSPACNZ<4) にすると同じバンチが 1% で合う
    (x, y), out = _pspac_rms("PSPACNX = 128; PSPACNY = 128; PSPACDX = 1e-4; PSPACDY = 1e-4; PSPACNZ = 1;", 0.02)
    assert abs(x / _pspac_analytic(0.02) - 1) < 0.015
    # 横のメッシュが粗い (半径が 1 セル): −13%、警告
    (x, y), out = _pspac_rms("PSPACNX = 128; PSPACNY = 128; PSPACDX = 1e-3; PSPACDY = 1e-3; PSPACNZ = 1;", 0.1)
    assert -0.2 < x / _pspac_analytic(0.1) - 1 < -0.08
    assert 'transverse mesh is too coarse' in out


def test_emittance_no_longer_warns_about_bend_f1(v2_mode):
    """(v2next B1 で変更) Emittance[]・EMIT・CALC6D は BEND の F1・FRINGE を入れるようになった (tbende の tblfre・tbendef1)。旧 A1 の警告は出さない"""
    src = _BEND_F1_RING.split('USE RING')[0] + "USE RING; CELL; CALC; RFSW; RAD; NOFLUC;\ne = Emittance[];\nUSE RING;\ne = Emittance[];\n"
    assert 'Emittance[]-Warning' not in _out(src)
    compat.set_mode('k64')
    try:
        assert 'Emittance[]-Warning' not in _out(src)
    finally:
        compat.set_mode('v2')
