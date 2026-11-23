# -*- coding: utf-8 -*-
"""担当 C2: SOL 領域 (bz≠0) の中の CAVI・RF つき MULT の TRPT の加速、wake、AUTOPHI、RF つき MULT の Emittance[]。
担当 C の実装では、どれも明確なエラーだった。

根拠 (gs を使わない):
  - 加速: 一様な縦電場 Ez(t) = −(V/L) sin(ω (t − t_ref) + φ) と一様な Bz の中の運動方程式を、基準粒子と一緒に数値積分した結果と一致。
  - CALC (TRPT) の軌道 = 追跡。加速の後の領域の要素は bz をその場所の基準運動量で規格化する。
  - wake: 領域内の空洞の中で分割ごとに入れた蹴りの合計 = 領域の外の同じ空洞。
  - AUTOPHI: 位相の基準を設計軌道の到着時刻 (Twiss の DZ) にする。設計軌道の粒子は公称の加速を受け、ずれた粒子は
    −V [sin(φ − ω (z − DZ)/v) − sin φ]。
  - RF つき MULT の Emittance[]: 多極なしなら CAVI と同じ。多極つきは追跡の 1 周の行列 (数値微分) の固有値と同じチューン。

新しく見つけた K64 の誤り: tcav.f・tmulti.f の wake の重み fw0·(ws(n) + ws(n−1)·0.5) (分割が 2 以上だと合計が 1.5 倍)。v2 は台形の重み。
"""
import io
import math
import re

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, emit6d as _e, track6n as tn, track6n_ext as x

ME = 0.51099895e6
CL = 299792458.0
P = 1.0e9


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, I = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _vals(out, tag):
    l = [s for s in out.splitlines() if s.startswith(tag + ' =')][0]
    return [float(v) for v in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', l.split('=', 1)[1])]


def _elem(defs):
    import pysad.ffscmd as F
    holder = {}
    orig = F.FFS.use

    def use(self, *a, **k):
        r = orig(self, *a, **k)
        holder['f'] = self
        return r
    F.FFS.use = use
    try:
        F.run_sad_string('MOMENTUM = 1GEV;\n' + defs + '\nMARK ST=();\nLINE BL=(ST E1);\nFFS;\nUSE BL;\nExit;\n', out=io.StringIO())
    finally:
        F.FFS.use = orig
    return [e for e in holder['f'].els if e.name == 'E1'][0]


class _Col:
    def __init__(self):
        self.items = []

    def _add(self, k, p, name, cnt, sm):
        self.items.append((k, p))


def _ode(c, L, bz, volt, w, phic, p0):
    """s を独立変数に、試験粒子と基準粒子 (軸上、δ = 0) の時刻 (ct) と運動量 (mc 単位) を積分。
    dPz/ds = Ez γ/Pz、dP⊥/ds = (P × b ẑ)⊥/Pz、Ez = −(V/(L·mc²)) sin(ω (t − t_ref) + φ)。出口は基準粒子の運動量 p2 で規格化"""
    eps = volt / ME / L
    b = bz * p0

    def rhs(s, u):
        xx, yy, t, px, py, pz, tr, pzr = u
        g = math.sqrt(1 + px * px + py * py + pz * pz)
        gr = math.sqrt(1 + pzr * pzr)
        ez = -eps * math.sin(w * (t - tr) + phic)
        return [px / pz, py / pz, g / pz, py * b / pz, -px * b / pz, ez * g / pz, gr / pzr, -eps * math.sin(phic) * gr / pzr]
    pt = p0 * (1 + c[5])
    pxk = (c[1] + 0.5 * bz * c[2]) * p0
    pyk = (c[3] - 0.5 * bz * c[0]) * p0
    pz = math.sqrt(pt * pt - pxk * pxk - pyk * pyk)
    v = pt / math.sqrt(1 + pt * pt)
    sol = solve_ivp(rhs, (0, L), [c[0], c[2], -c[4] / v, pxk, pyk, pz, 0.0, p0], rtol=1e-13, atol=1e-16, method='DOP853')
    xx, yy, t, px, py, pz, tr, p2 = sol.y[:, -1]
    pp = math.sqrt(px * px + py * py + pz * pz)
    bz2 = bz * p0 / p2
    return np.array([xx, px / p2 - 0.5 * bz2 * yy, yy, py / p2 + 0.5 * bz2 * xx, -(t - tr) * pp / math.sqrt(1 + pp * pp), pp / p2 - 1]), p2


@pytest.mark.parametrize('volt,phi,tol', [(2e7, -1.2, 5e-10), (2e8, 1.0, 3e-8)])
def test_acceleration_in_solenoid_matches_equation_of_motion(v2_mode, volt, phi, tol):
    """領域内の CAVI の TRPT の加速 (1 GeV から +19 MV、−168 MV)。bz = 0.3 と 3 1/m。出口の基準運動量と、出口の基準運動量で
    規格化した座標 (正準運動量は出口の bz·p0/p2) が、運動方程式の積分と合う。差は K64 の空洞の分割の精度 (200 MV で 7e-9)"""
    el = _elem('CAVI E1 = (L=0.5 VOLT=%g FREQ=1.3E9 PHI=%g DISFRIN=1);' % (volt, phi))
    w = 2 * math.pi * 1.3e9 / CL
    for bz in (0.3, 3.0):
        for c in ([0, 0, 0, 0, 0, 0], [1e-4, 2e-5, -1e-4, 1e-5, 1e-3, 1e-3], [2e-4, -3e-5, 5e-5, 4e-5, -2e-3, -5e-3]):
            ph = _e._Phys(P, 100.0)
            col = _Col()
            x._compile_sol_v2(col, el, 'Cavity', bz, ph, 2)
            k, p = col.items[0]
            S = tn.make_state(np.array([c], float), ph)
            ph2 = x.step_kind(k, p, S, 1, ph, 0.0, True, False, False, tn.make_rs(1), tm=10.0)
            o = tn.state_out(S, 1)[0]
            r, p2 = _ode(c, 0.5, bz, volt, w, phi, ph.p0)           # (この ph は charge = +1: phic = PHI)
            assert abs(ph2.p0 / p2 - 1) < 1e-14
            assert np.abs(o - r).max() < tol, (bz, c, o - r)


_TL = '''MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 EX=0.1 DX=1E-4 DPX=2E-5 DY=-1E-4 DPY=1E-5) M1=();
DRIFT D1 = (L=1.0) DS=(L=0.3);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25) QS=(L=0.2 K1=0.1);
CAVI CA = (L=0.5 VOLT=2E8 FREQ=1.3E9 PHI=-1.2);
MULT MA = (L=0.4 VOLT=1E8 FREQ=1.3E9 PHI=-1.0 K1=0.05);
SOL SA=(L=0 BZ=1 BOUND=1 GEO=1) SB=(L=0 BZ=0 BOUND=1);
LINE TL = (M0 QF D1 SA DS CA DS QS DS MA DS SB D1 QD D1 M1);
FFS;
USE TL;
TRPT; RFSW; INS; CALC;
$FORM="20.14";
PageWidth=1999;
Print["TW =",Twiss[{"DX","DPX","DY","DPY","DZ","DDP"},"M1"]];
Print["GB =",LINE["GAMMABETA",{1,"M1"}]];
b={{1E-4},{2E-5},{-1E-4},{1E-5},{0},{0},{1}};
r=TrackParticles[{1,b},LINE["LENGTH"]][[2]];
Print["TR =",Table[r[[k,1]],{k,6}]];
Exit;
'''


def test_trpt_calc_orbit_through_accelerating_region_equals_tracking(v2_mode):
    """輸送路 (1 GeV → 1.46 GeV): 領域 (BZ = 1 T) の中に CAVI (+186 MV)・四極・RF つき MULT (+84 MV、K1 つき)。
    CALC (TRPT) の軌道 (入口のずれ 0.1 mm から) = 追跡。運動量は 1 + (ΔE の和)/E"""
    ffs, out = _out(_TL)
    tw, gb, tr = _vals(out, 'TW'), _vals(out, 'GB'), _vals(out, 'TR')
    for a, b in zip(tw, tr):
        assert abs(a - b) < 2e-9 + 1e-6 * abs(b), (tw, tr)
    de = (2e8 * math.sin(1.2) + 1e8 * math.sin(1.0)) / ME
    h1 = math.sqrt(1 + gb[0] ** 2) + de
    assert abs(gb[1] / math.sqrt(h1 * h1 - 1) - 1) < 1e-12
    assert abs(tw[0]) > 1e-5                                   # 軌道は 0 でない (結合つきで運ばれている)


_WK = '''MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3) M1=();
DRIFT D1 = (L=1.0) DS=(L=0.3);
CAVI CA = (L=%s VOLT=2E7 FREQ=1.3E9 PHI=-1.2);
MULT MA = (L=0.5 VOLT=2E7 FREQ=1.3E9 PHI=-1.2 K1=0.1);
SOL SA=(L=0 BZ=%s BOUND=1 GEO=1) SB=(L=0 BZ=0 BOUND=1);
LINE TL = (M0 D1 SA DS %s DS SB D1 M1) T1=(M0 D1 DS %s DS D1 M1);
FFS;
$FORM="20.14";
PageWidth=1999;
PBUNCH=1e11;
WakeFunction[Longitudinal,"%s"]={{-1,0},{0,0},{0,1e13},{0.001,1e13},{0.01,5e12},{1,0}};
WakeFunction[Transverse,"%s"]={{-1,0},{0,0},{0.001,1e15},{0.01,5e14},{1,0}};
b={{1e-4,2e-4,-1e-4,3e-4},{0,2e-5,0,1e-5},{1e-4,-1e-4,1e-4,2e-4},{0,1e-5,0,0},{0,1e-3,-1e-3,2e-3},{0,1e-3,0,-1e-3},{1,1,1,1}};
USE %s;
TRPT; RFSW; INS; CALC; LWAKE; TWAKE;
r=TrackParticles[{1,b},LINE["LENGTH"]][[2]];
Print["W =",r[[{1,2,3,4,5,6}]]];
NOLWAKE; NOTWAKE;
r=TrackParticles[{1,b},LINE["LENGTH"]][[2]];
Print["N =",r[[{1,2,3,4,5,6}]]];
Exit;
'''


def _wk(line, bz, el, L='0.5', volt='2E7'):
    ffs, out = _out((_WK % (L, bz, el, el, el, el, line)).replace('CAVI CA = (L=%s VOLT=2E7' % L, 'CAVI CA = (L=%s VOLT=%s' % (L, volt)))
    return np.array(_vals(out, 'W')).reshape(6, 4) - np.array(_vals(out, 'N')).reshape(6, 4)


@pytest.mark.parametrize('el', ['CA', 'MA'])
def test_wake_inside_region_cavity(v2_mode, el):
    """領域内の CAVI・RF つき MULT の wake (tsol.f の ewak: tmultiacc の中で分割ごとに txwake)。bz がごく小さい領域と領域の外で、
    エネルギーの変化 (縦の wake) の合計が同じ。横の蹴りは分割して入るので位置の変化が少し違う (蹴りの合計は 2% 以内)。
    BZ = 1 T でも縦の wake の効果は同じ (位置によらない表)"""
    a = _wk('TL', '1E-9', el)
    c = _wk('T1', '0', el)
    d = _wk('TL', '1', el)
    assert np.abs(c[5]).max() > 5e-5
    assert np.abs(a[5] - c[5]).max() < 1e-9 * 1e3 and np.abs(d[5] - c[5]).max() < 1e-8
    assert np.abs(a[1] - c[1]).max() < 0.02 * np.abs(c[1]).max()
    assert np.abs(a[0] - c[0]).max() < 0.1 * np.abs(c[0]).max()


def test_wake_weight_sum_is_one_in_v2_and_one_and_half_in_k64():
    """空洞の中の wake の重み (K64 の tcav.f: fw0·(ws(n) + ws(n−1)·0.5)) は、分割が 2 以上だと合計が (1.5 n − 0.5)/n 倍になる。
    v2 は台形の重みで、薄い空洞 (分割 1) と同じエネルギーの変化。k64 モードは K64 のまま"""
    old = compat._mode
    try:
        compat.set_mode('v2')
        # 200 MV の空洞 (加速が大きいので 4 分割以上になる) と 20 MV (分割 1) を比べる:
        # wake の表は位置によらないので、エネルギーの変化 (eV) は電圧によらず同じはず
        def dE(volt):
            ffs, out = _out((_WK % ('0.5', '0', 'CA', 'CA', 'CA', 'CA', 'T1')).replace('CAVI CA = (L=0.5 VOLT=2E7', 'CAVI CA = (L=0.5 VOLT=%s' % volt)
                            .replace('Print["W =",r[[{1,2,3,4,5,6}]]];', 'Print["W =",r[[{1,2,3,4,5,6}]]]; Print["G =",LINE["GAMMABETA","M1"]];', 1))
            d = np.array(_vals(out, 'W')).reshape(6, 4)[5] - np.array(_vals(out, 'N')).reshape(6, 4)[5]
            return d[:3] * _vals(out, 'G')[0]                     # δ × 出口の運動量 = エネルギーの変化 (mc² 単位)
        e20, e200 = dE('2E7'), dE('2E8')
        assert np.abs(e200 / e20 - 1).max() < 5e-3
        compat.set_mode('k64')
        r = dE('2E8') / dE('2E7')
        assert np.all(r > 1.25) and np.all(r < 1.5)
    finally:
        compat.set_mode(old)


_AP = '''MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 DZ=0.003) M1=() M2=();
DRIFT D1 = (L=1.0) DS=(L=0.3);
CAVI CA = (L=0.5 VOLT=2E7 FREQ=1.3E9 PHI=-1.2 AUTOPHI=%s);
SOL SA=(L=0 BZ=1 BOUND=1 GEO=1) SB=(L=0 BZ=0 BOUND=1);
LINE TL = (M0 D1 CA M1 D1 SA DS CA DS SB D1 M2);
FFS;
USE TL;
TRPT; RFSW; INS; CALC;
$FORM="20.14";
PageWidth=1999;
Print["TW =",Twiss[{"DZ","DDP"},{"M1","M2"}]];
Print["GB =",LINE["GAMMABETA",{1,"M1","M2"}]];
b={{0,0},{0,0},{0,0},{0,0},{0.003,0},{0,0},{1,1}};
r=TrackParticles[{1,b},LINE["POSITION","M1"]][[2]];
Print["T1 =",r[[{5,6}]]];
r=TrackParticles[{1,b},LINE["LENGTH"]][[2]];
Print["T2 =",r[[{5,6}]]];
Exit;
'''


def test_autophi_locks_cavity_phase_to_design_orbit_arrival(v2_mode):
    """AUTOPHI: 設計軌道が z = DZ = 3 mm (位相 ω·DZ = 0.082 rad) だけ早く着く輸送路。AUTOPHI なしでは設計軌道の粒子が
    公称と違う加速を受ける (CALC の DDP = 5.2e-4、追跡も同じ)。AUTOPHI ありでは設計軌道の粒子 (z = DZ) の δ が 0 のままで、
    z = 0 の粒子は δ = −V [sin(φ − ω(0 − DZ)/v) − sin φ]/(β p c) (φ = PHI·charge)。領域の外の空洞と領域 (BZ = 1 T) の中の空洞の両方"""
    ffs, out0 = _out(_AP % '0')
    tw0, t1_0 = _vals(out0, 'TW'), _vals(out0, 'T1')
    assert abs(tw0[2] - 5.1978e-4) < 1e-7 and abs(t1_0[2] / tw0[2] - 1) < 1e-5 and abs(t1_0[3]) < 1e-8
    ffs, out1 = _out(_AP % '1')
    tw, gb, t1, t2 = _vals(out1, 'TW'), _vals(out1, 'GB'), _vals(out1, 'T1'), _vals(out1, 'T2')
    assert abs(tw[2]) < 1e-12 and abs(tw[3]) < 1e-9                    # CALC: 設計軌道は公称の加速
    assert abs(t1[2]) < 1e-9 and abs(t2[2]) < 1e-9                     # 追跡: z = DZ の粒子
    w = 2 * math.pi * 1.3e9 / CL
    h0 = math.sqrt(1 + gb[0] ** 2)
    v = gb[0] / h0
    dh = -(2e7 / ME) * (math.sin(-1.2 + w * 0.003 / v) - math.sin(-1.2))     # z = 0 の粒子の、公称との差
    h1 = math.sqrt(1 + gb[1] ** 2)
    want = (math.sqrt((h1 + dh) ** 2 - 1) - gb[1]) / gb[1]
    assert abs(t1[3] / want - 1) < 2e-4, (t1[3], want)
    # gs (TRPT、領域の外の空洞): AUTOPHI=1 の追跡 T1 = {{.002999998831, -1.22409844074E-9}, {-1.5352159821E-10, -6.41977395147E-4}} と全桁
    assert abs(t1[3] / -6.41977395147e-4 - 1) < 1e-9 and abs(t1[0] - .002999998831) < 1e-11
    assert abs(t2[3]) > 1.9 * abs(t1[3]) * (gb[1] / gb[2]) * 0.9       # 領域の中の空洞でも同じ向きにずれる
    compat.set_mode('k64')
    try:
        with pytest.raises(Exception):
            _out(_AP % '1')                                            # k64 モードは従来どおり未対応
    finally:
        compat.set_mode('v2')


_RING = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.25);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI CA = (L=0.3 VOLT=4.0E6 FREQ=5.0E8 PHI=0 %s);
MULT MA = (L=0.3 VOLT=4.0E6 FREQ=5.0E8 PHI=0 %s);
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA DS %s DS SB DS);
LINE RINGQ = (ORIGIN 12*CELLQ CELLS 11*CELLQ);
FFS;
USE RINGQ;
CELL; CALC;
$FORM="20.14";
PageWidth=1999;
RFSW; %s; NOFLUC;
e=Emittance[];
Print["EM =",Emittances/.e," ",EnergyLossU0/.e," ",Tunes/.e];
Exit;
'''


def _em(cav, k, bz, el, flags='RAD; RADCOD'):
    ffs, out = _out(_RING % (cav, k, bz, el, flags))
    return ffs, np.array(_vals(out, 'EM'))


def test_emittance_with_rf_mult_in_region(v2_mode):
    """Emittance[]: 領域 (BZ = 2 T) の中の RF つき MULT。多極なしなら CAVI と同じ (1e-9)。bz → 0 の領域の中の CAVI は領域の外の
    空洞の写像と同じ (端の RF のフリンジも入る: チューン 1e-8)。K1・K2 つきの RF MULT は、追跡 (RFSW、放射なし) の 1 周の行列の
    固有値と同じチューン"""
    ffs, a = _em('', '', '2', 'CA')
    ffs, b = _em('', '', '2', 'MA')
    assert np.abs(a[[0, 1, 2, 3]] / b[[0, 1, 2, 3]] - 1).max() < 1e-8 and np.abs(a[4:7] - b[4:7]).max() < 1e-9
    ffs, z0 = _em('', '', '0', 'CA')
    ffs, z1 = _em('', '', '1E-6', 'CA')
    assert np.abs(z0[4:7] - z1[4:7]).max() < 1e-8 and abs(z0[0] / z1[0] - 1) < 1e-7
    ffs, c = _em('', 'K1=0.05 K2=1', '2', 'MA', flags='NORAD; NORADCOD')
    ffs.flags['RFSW'] = True
    ffs.flags['RAD'] = False
    n = len(ffs.els) + 1
    h = [1e-6, 1e-7, 1e-6, 1e-7, 1e-5, 1e-6]
    pts = [[0.0] * 6]
    for j in range(6):
        for s in (h[j], -h[j]):
            q = [0.0] * 6
            q[j] = s
            pts.append(q)
    pts = np.array(pts)
    out = ffs.fn_track([[1.0, [list(pts[:, j]) for j in range(6)] + [[1.0] * len(pts)]], float(n)])
    o = np.array(out[1][:6]).T
    M = np.array([(o[1 + 2 * j] - o[2 + 2 * j]) / (2 * h[j]) for j in range(6)]).T
    fr = sorted(abs(math.atan2(e.imag, e.real)) / (2 * math.pi) for e in np.linalg.eigvals(M))[::2]
    want = sorted(min(abs(v) % 1, 1 - abs(v) % 1) for v in c[4:7])
    for u, v in zip(fr, want):
        assert abs(u - v) < 3e-6, (fr, want)


def test_emittance_autophi_uses_the_cavity_phase_at_the_orbit(v2_mode):
    """AUTOPHI (Emittance[]): 軌道は位相 PHI (= 0) で空洞を通る。AUTOPHI なしのリングは、放射損失を補う同期位相 φs
    (sin φs = U0/V) で通るので、RF の傾きが cos φs 倍: シンクロトロンチューンの比 = 1/sqrt(cos φs)。横のチューンは同じ (1e-4)"""
    ffs, a = _em('', '', '2', 'CA', flags='NORAD; NORADCOD')
    ffs, b = _em('AUTOPHI=1', '', '2', 'CA', flags='NORAD; NORADCOD')
    phis = math.asin(a[3] / 4.0e6)
    assert abs(b[6] / a[6] - 1.0 / math.sqrt(math.cos(phis))) < 1e-3
    assert np.abs(a[4:6] - b[4:6]).max() < 1e-4
