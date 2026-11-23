# -*- coding: utf-8 -*-
"""担当 C2: SOL のずれ・回転 (DX DY DZ CHI1..3 = tsol.f の trots)、端の放射 F1 (tradks)、領域の幾何。

K64 では、GEO≠0 の SOL に DX DY (位置) と DPX DPY (向き) を書くと、tsgeo がその軌道を領域の中で追って幾何を作り、
両端の SOL の CHI1..3 と、もう一方の SOL の DX DY DZ を決める。追跡 (tsol.f) は入口と出口でその値の座標変換 trots をする。
pysad の v1 は幾何だけ対応で、追跡・Twiss・Emittance[] は未対応だった。v2 で実装した。

根拠 (gs を使わない、独立な計算):
  - 実験室系の合成: 幾何 (GX GY GZ と軸の向き) から、ソレノイドの軸と両端の面を実験室系に置く。粒子を入口の基準面から
    ソレノイドの入口の面まで直進させ、端の蹴り (力学的運動量 = 正準運動量 + (bz/2)(y, −x))、一様磁場 (と四極) の中の
    ローレンツ力 dp/dl = (p × B)/|p| の数値積分 (solve_ivp)、出口の蹴り、出口の基準面までの直進をつなぐ。
    z は z = −v (t − t0): Δz = (v/v0)·(設計の行路長) − (粒子の行路長)。
  - 設計軌道の粒子 (0,0,0,0,0,0) は領域の出口で (0,0,0,0,0,0) に戻る。
  - F1 の放射: 端の蹴り Δθ = bz r/2 を長さ F1 の一様な偏向とした放射損失 (2/3) r_e γ³ Δθ²/F1。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, geo as G

RE = 2.8179403262e-15
ME = 0.51099895e6


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, I = run_sad_string(src, out=buf)
    return ffs, I, buf.getvalue()


def _track(ffs, coords, k0, k1, rad=False):
    """要素 k0 の入口から k1 の入口まで (1 始まり)。coords (n, 6)"""
    c = np.asarray(coords, dtype=float)
    cols = [list(c[:, j]) for j in range(6)] + [[1.0] * len(c)]
    ffs.flags['RAD'] = bool(rad)
    ffs.flags['FLUC'] = False
    out = ffs.fn_track([[float(k0), cols], float(k1)])
    return np.array(out[1][:6]).T


def _lab_reference(ffs, c, k1, ke, quads=(), frin=None):
    """実験室系の合成で、SOL (位置 k1) の入口から終端の SOL (位置 ke) の出口までを追う。c は 6 成分 (px py は p0 で規格化した運動量)。
    quads: [(s0, s1, k)] か [(s0, s1, k, by, bx)] ソレノイドの系の s (入口の面から) の区間と K1/L (と一様な双極 K0/L)。frin: 端の非線形フリンジ (pysad の tsfrin1) を使うとき (dbz_in, dbz_out)"""
    from pysad import emit6d as _e, track6 as _t6, track6n as tn, track6n_ext as x
    mom = ffs.momentum()
    res = G.compute_geo(ffs.els, warn=lambda s: None, momentum_eV=mom)
    g0, g1, g2, g3 = res.geo[k1 - 1], res.geo[k1], res.geo[ke - 1], res.geo[ke]
    ds_design = res.pos[ke] - res.pos[k1 - 1]
    ph = _e._Phys(mom, 100.0)
    bz = _t6._tfbzs(ffs.els[k1 - 1], ph)
    x0, px0, y0, py0, z0, d0 = c
    pt = 1.0 + d0
    r = g0[:, 3] + g0[:, 0] * x0 + g0[:, 1] * y0
    p = g0[:, 0] * px0 + g0[:, 1] * py0 + g0[:, 2] * math.sqrt(pt * pt - px0 * px0 - py0 * py0)
    ex, ey, ez = g1[:, 0], g1[:, 1], g1[:, 2]
    path = 0.0

    def to_plane(r, p, o, n):
        t = -np.dot(r - o, n) / np.dot(p, n) * pt          # 行路長 (p/|p| の向き)
        return r + p / pt * t, t

    def edge(r, p, o, sgn, dbz):
        xs, ys = np.dot(r - o, ex), np.dot(r - o, ey)
        a, b, cz = np.dot(p, ex), np.dot(p, ey), np.dot(p, ez)
        if frin is not None and dbz != 0.0:
            S = tn.make_state(np.array([[xs, a if sgn > 0 else a - 0.5 * bz * ys, ys, b if sgn > 0 else b + 0.5 * bz * xs, 0.0, d0]]), ph)
            x.tsfrin1(0, S, dbz)
            o_ = tn.state_out(S, 1)[0]
            xs2, a2, ys2, b2 = o_[0], o_[1], o_[2], o_[3]
            if sgn > 0:
                a2, b2 = a2 + 0.5 * bz * ys2, b2 - 0.5 * bz * xs2
            r = o + ex * xs2 + ey * ys2
            dz_f = o_[4]
        else:
            a2, b2 = a + sgn * 0.5 * bz * ys, b - sgn * 0.5 * bz * xs
            dz_f = 0.0
        return r, ex * a2 + ey * b2 + ez * math.sqrt(pt * pt - a2 * a2 - b2 * b2), dz_f

    r, t = to_plane(r, p, g1[:, 3], ez)
    path += t
    r, p, dzf1 = edge(r, p, g1[:, 3], +1.0, frin[0] if frin else 0.0)

    stot = np.dot(g2[:, 3] - g1[:, 3], ez)
    bounds = sorted({0.0, stot} | {q[0] for q in quads} | {q[1] for q in quads})
    u = np.concatenate([r, p])
    for sa, sb in zip(bounds[:-1], bounds[1:]):
        ev = lambda l, u, sb=sb: np.dot(u[:3] - g1[:, 3], ez) - sb
        ev.terminal = True
        ev.direction = 1.0
        # 区間の途中で場が切り替わらないよう、区間の中点の s で場を決める
        sm = 0.5 * (sa + sb)

        def rhs_seg(l, u, sm=sm):
            rr, pp = u[:3], u[3:]
            xs, ys = np.dot(rr - g1[:, 3], ex), np.dot(rr - g1[:, 3], ey)
            B = bz * ez
            for q in quads:
                if q[0] <= sm < q[1]:
                    B = B + q[2] * (ys * ex + xs * ey)
                    if len(q) > 3:
                        B = B + q[3] * ey + q[4] * ex          # 一様な双極 (By, Bx)
            return np.concatenate([pp / pt, np.cross(pp, B) / pt])
        sol = solve_ivp(rhs_seg, (0.0, 2.0 * (sb - sa) + 1.0), u, events=ev, rtol=1e-13, atol=1e-16, method='DOP853')
        assert sol.t_events[0].size == 1
        u = sol.y_events[0][0]
        path += float(sol.t_events[0][0])
    r, p = u[:3], u[3:]
    r, p, dzf2 = edge(r, p, g2[:, 3], -1.0, frin[1] if frin else 0.0)
    r, t = to_plane(r, p, g3[:, 3], g3[:, 2])
    path += t
    p0 = mom / ME
    v0 = p0 / math.sqrt(1.0 + p0 * p0)
    pp_ = p0 * pt
    v = pp_ / math.sqrt(1.0 + pp_ * pp_)
    z = z0 + (v / v0) * ds_design - path + dzf1 + dzf2
    return np.array([np.dot(r - g3[:, 3], g3[:, 0]), np.dot(p, g3[:, 0]), np.dot(r - g3[:, 3], g3[:, 1]), np.dot(p, g3[:, 1]), z, d0])


_SRC = '''MOMENTUM = %s;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.3 DISFRIN=1) QD = (L=0.3 K1=-0.25 DISFRIN=1);
SOL SA = (L=0 BZ=%s BOUND=1 %s) SB = (L=0 BZ=0 BOUND=1 %s);
MARK ST = (BX=5 BY=5);
LINE BL = (ST DR SA %s SB DR);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''

_PARTS = np.array([[0, 0, 0, 0, 0, 0],
                   [1e-3, 2e-4, -1e-3, 1e-4, 1e-3, 1e-3],
                   [-2e-3, -3e-4, 5e-4, 4e-4, -2e-3, -5e-3],
                   [3e-3, 1e-3, 2e-3, -1e-3, 0, 2e-2]], dtype=float)


@pytest.mark.parametrize('case', ['geo_in', 'geo_out', 'geo_in_quads', 'geo_out_quads', 'lowE'])
def test_tilted_offset_solenoid_matches_lab_frame_lorentz_integration(v2_mode, case):
    """ずれて傾いたソレノイド (端のフリンジは DISFRIN=1: 正準運動量が連続な硬い端) を、実験室系のローレンツ力の積分と比べる。
    GEO が入口の SOL / 出口の SOL (入口の SOL に DZ が入る)、領域内に四極、低エネルギー (30 MeV: v ≠ c の z の項)"""
    off = 'GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002 DISFRIN=1'
    mom, bzv = '3 GEV', '3'
    inner = 'DR2 DR2 DR2' if 'quads' not in case else 'DR2 QF DR2 QD DR2'
    if case == 'lowE':
        mom, bzv = '30 MEV', '0.05'
    a, b = (off, 'DISFRIN=1') if 'out' not in case else ('DISFRIN=1', off)
    ffs, I, out = _run(_SRC % (mom, bzv, a, b, inner))
    k1 = [i for i, e in enumerate(ffs.els) if e.name == 'SA'][0] + 1
    ke = [i for i, e in enumerate(ffs.els) if e.name == 'SB'][0] + 1
    got = _track(ffs, _PARTS, k1, ke + 1)
    quads = [(0.5, 0.8, 0.3 / 0.3), (1.3, 1.6, -0.25 / 0.3)] if 'quads' in case else []
    assert np.max(np.abs(got[0])) < 1e-13                     # 設計軌道は設計軌道に戻る
    for c, g in zip(_PARTS, got):
        ref = _lab_reference(ffs, c, k1, ke, quads)
        # ドリフトだけなら丸めまで一致。領域内の四極はカーネル (tsolqu の bz≠0 の枝) の精度 (3 mm・1 mrad・δ=2% の粒子で 4e-9 m: 運動学的な高次の項)
        tol = 1e-8 if quads else 2e-12
        assert np.max(np.abs(g[:4] - ref[:4])) < tol, (case, c, g - ref)
        assert abs(g[4] - ref[4]) < tol, (case, c, g[4] - ref[4])


def test_region_frames_are_inverse_pair_in_the_field_free_limit(v2_mode):
    """BZ=0 の領域 (ずれ・傾きだけ): 入口と出口の座標変換は互いに打ち消し、領域全体は長さ ΣL の直進 (ドリフト) と同じ"""
    off = 'GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002'
    ffs, I, out = _run(_SRC % ('3 GEV', '0', off, '', 'DR2 DR2 DR2'))
    k1 = [i for i, e in enumerate(ffs.els) if e.name == 'SA'][0] + 1
    ke = [i for i, e in enumerate(ffs.els) if e.name == 'SB'][0] + 1
    got = _track(ffs, _PARTS, k1, ke + 1)
    res = G.compute_geo(ffs.els, warn=lambda s: None, momentum_eV=ffs.momentum())
    g0, g3 = res.geo[k1 - 1], res.geo[ke]
    for c, g in zip(_PARTS, got):
        pt = 1.0 + c[5]
        r = g0[:, 3] + g0[:, 0] * c[0] + g0[:, 1] * c[2]
        p = g0[:, 0] * c[1] + g0[:, 1] * c[3] + g0[:, 2] * math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)
        t = -np.dot(r - g3[:, 3], g3[:, 2]) / np.dot(p, g3[:, 2])
        r2 = r + p * t
        ref = [np.dot(r2 - g3[:, 3], g3[:, 0]), np.dot(p, g3[:, 0]), np.dot(r2 - g3[:, 3], g3[:, 1]), np.dot(p, g3[:, 1])]
        assert np.max(np.abs(g[:4] - ref)) < 1e-14
    # 軸を外れて傾いていても、磁場がなければ幾何は直線: 出口の軸は入口の軸と同じ向きで、位置は 1.5 m 先
    # (領域の長さ 1.5 m はソレノイドの軸に沿った長さなので、設計軌道に沿っては 1.5/cos(傾き))
    cs = math.cos(0.003) * math.cos(0.002)
    assert np.max(np.abs(g3[:, :3] - g0[:, :3])) < 1e-12 and np.max(np.abs(g3[:, 3] - g0[:, 3] - 1.5 / cs * g0[:, 2])) < 1e-12


def test_tilted_solenoid_with_bend_and_mult_inside(v2_mode):
    """ずれて傾いた領域の中の BEND (直線の座標の一様な双極)・MULT (K0 + K1): 幾何 (tsgeo1 の icBEND・icMULT) が設計軌道を
    その中でも追うので、設計軌道の粒子は出口で 0 に戻る。ほかの粒子は実験室系のローレンツ力の積分と合う。GEO が出口側でも同じ"""
    for a, b in (('GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002 DISFRIN=1', 'DISFRIN=1'),
                 ('DISFRIN=1', 'GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002 DISFRIN=1')):
        src = '''MOMENTUM = 3 GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
BEND BK = (L=0.4 K0=0.002 ROTATE=0.3);
MULT MK = (L=0.3 K0=-0.001 K1=0.2 DISFRIN=1);
SOL SA = (L=0 BZ=3 BOUND=1 %s) SB = (L=0 BZ=0 BOUND=1 %s);
MARK ST = (BX=5 BY=5);
LINE BL = (ST DR SA DR2 BK DR2 MK DR2 SB DR);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
''' % (a, b)
        ffs, I, out = _run(src)
        k1 = [i for i, e in enumerate(ffs.els) if e.name == 'SA'][0] + 1
        ke = [i for i, e in enumerate(ffs.els) if e.name == 'SB'][0] + 1
        got = _track(ffs, _PARTS, k1, ke + 1)
        assert np.max(np.abs(got[0])) < 1e-12
        fld = [(0.5, 0.9, 0.0, 0.002 / 0.4 * math.cos(0.3), 0.002 / 0.4 * math.sin(0.3)),
               (1.4, 1.7, 0.2 / 0.3, -0.001 / 0.3, 0.0)]
        for c, g in zip(_PARTS, got):
            ref = _lab_reference(ffs, c, k1, ke, fld)
            assert np.max(np.abs(g[:5] - ref[:5])) < 1e-8, (c, g - ref)


# ---------------------------------------------------------------- F1 (端の放射)
_F1 = '''MOMENTUM = 3 GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
SOL SA = (L=0 BZ=3 BOUND=1 GEO=1 F1=%s) SM=(L=0 BZ=1 F1=%s) SB = (L=0 BZ=0 BOUND=1 F1=%s);
MARK ST = (BX=5 BY=5);
LINE BL = (ST DR SA DR2 SM DR2 SB DR);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''


def test_sol_f1_edge_radiation_matches_uniform_kick_formula(v2_mode):
    """SOL の F1 (tradks): 端の力学的運動量の向きの変化 Δθ = Δbz·r/2 を長さ F1 の一様な偏向とみた放射損失
    Δδ = −(2/3) r_e γ³ Δθ²/F1。軸に平行に半径 r0 で入った粒子の半径は、ラーモア系で r'' = −(bz/2)² r に従う:
    入口 r0、領域内の SOL (0.5 m 先、bz → bz/3) で r0 cos(bz/4)、出口 (さらに 0.5 m) は下の式。F1 を 2 倍にすると半分"""
    c = np.array([[5e-3, 0, 0, 0, 0, 0], [0, 0, 3e-3, 0, 0, 0]])
    res = {}
    for k, f in (('0', (0, 0, 0)), ('a', (0.1, 0, 0)), ('m', (0, 0.1, 0)), ('b', (0, 0, 0.1)), ('a2', (0.2, 0, 0))):
        ffs, I, out = _run(_F1 % f)
        res[k] = _track(ffs, c, 3, 8, rad=True)[:, 5]
    p0 = 3e9 / ME
    gam = math.sqrt(1 + p0 * p0)
    bz = 3 * 299792458.0 / 3e9
    k1, k2 = bz / 2, bz / 6
    for j, r0 in enumerate((5e-3, 3e-3)):
        a = -(2 / 3) * RE * gam ** 3 * (bz * r0 / 2) ** 2 / 0.1
        rm = r0 * math.cos(k1 * 0.5)
        rpm = -r0 * k1 * math.sin(k1 * 0.5)
        rb = rm * math.cos(k2 * 0.5) + rpm / k2 * math.sin(k2 * 0.5)
        am = -(2 / 3) * RE * gam ** 3 * ((2 * bz / 3) * rm / 2) ** 2 / 0.1
        ab = -(2 / 3) * RE * gam ** 3 * ((bz / 3) * rb / 2) ** 2 / 0.1
        assert abs((res['a'][j] - res['0'][j]) / a - 1) < 1e-5
        assert abs((res['m'][j] - res['0'][j]) / am - 1) < 1e-4
        assert abs((res['b'][j] - res['0'][j]) / ab - 1) < 1e-4
        assert abs((res['a2'][j] - res['0'][j]) / (a / 2) - 1) < 1e-5


# ---------------------------------------------------------------- リング (CALC・Emittance[])
_RING = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.25);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QS1=(L=0.2 K1=0.1) QS2=(L=0.2 K1=-0.1);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1 %s F1=%s) SB = (L=0 BZ=0 BOUND=1 F1=%s);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA DS QS1 DS QS2 DS SB);
LINE RINGQ = (ORIGIN RF 12*CELLQ CELLS 11*CELLQ);
FFS;
USE RINGQ;
CELL; CALC;
$FORM="20.14";
PageWidth=1999;
Print["NU =",Twiss["NX","$$$"]/2/Pi," ",Twiss["NY","$$$"]/2/Pi];
Print["COD =",Twiss[{"DX","DPX","DY","DPY"},{1,"SA+1","SB","$$$"}]];
RFSW; RAD; RADCOD; NOFLUC;
e=Emittance[];
Print["EM =",Emittances/.e," ",EnergyLossU0/.e," ",DampingRate/.e];
Exit;
'''
_OFF = 'DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002'


def _vals(out, tag):
    import re
    l = [x for x in out.splitlines() if x.startswith(tag + ' =')][0]
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', l.split('=', 1)[1])]


def test_ring_with_tilted_solenoid_calc_orbit_tunes_and_emittance(v2_mode):
    """リング (3 GeV、周長 91.35 m、BZ = 2 T の領域 1.15 m に四極 2 台) にずれて傾いたソレノイド:
    CALC の閉軌道は設計軌道 (入口・出口で 0、領域の中はソレノイドの系で (DX, DY) から)、CALC のチューン = 追跡の 1 周の行列
    (数値微分) の固有値、追跡は CALC の軌道の上を 1 周して戻る。Emittance[]: 減衰率の和 = 2 U0/E·f0 の関係 (Robinson)、
    F1 による U0 の増分 = E·(2/3) r_e γ³ (bz/2)² (r_in² + r_out²)/F1 (r は閉軌道の、ソレノイドの軸からの距離)"""
    ffs, I, out = _run(_RING % ('2', _OFF, '0', '0'))
    nu = _vals(out, 'NU')
    cod = np.array(_vals(out, 'COD')).reshape(4, 4)
    assert np.max(np.abs(cod[:, 0])) < 1e-12 and np.max(np.abs(cod[:, 3])) < 1e-12
    assert abs(cod[0, 1] - 0.002) < 1e-8 and abs(cod[2, 1] + 0.001) < 1e-8          # 領域の入口: ソレノイドの系で (DX, DY)
    # 追跡 (NORAD、空洞は RFSW ON でも δ=0・z=0 の粒子を蹴らない位相) の 1 周の行列
    ffs.flags['RFSW'] = False
    n = len(ffs.els) + 1
    h = 1e-6
    base = [[0.0] * 6]
    pts = []
    for j in range(4):
        for s in (h, -h):
            c = [0.0] * 6
            c[j] = s
            pts.append(c)
    o = _track(ffs, np.array(base + pts), 1, n)
    assert np.max(np.abs(o[0, :4])) < 1e-12
    M = np.array([(o[1 + 2 * j, :4] - o[2 + 2 * j, :4]) / (2 * h) for j in range(4)]).T
    ev = np.linalg.eigvals(M)
    fr = sorted(abs(math.atan2(e.imag, e.real)) / (2 * math.pi) for e in ev)[::2]
    want = sorted(min(v % 1, 1 - v % 1) for v in nu)
    assert abs(fr[0] - want[0]) < 2e-7 and abs(fr[1] - want[1]) < 2e-7
    em0 = _vals(out, 'EM')
    u0, damp = em0[3], em0[4:7]
    E = 3e9
    f0 = 299792458.0 / 91.35
    # 減衰率 (1 周あたりの振幅の減衰率、負) の和 = −2 U0/E (Robinson の和則。pysad の DampingRate は 1 周あたり)
    assert abs(sum(damp) / (-2 * u0 / E) - 1) < 2e-3
    # F1 の放射
    # (RADCOD OFF で比べる: 閉軌道が設計軌道のままなので r が幾何から決まる。RADCOD ON では閉軌道が 2e-5 m 動いて 3% 変わる)
    nrc = lambda s: s.replace('RFSW; RAD; RADCOD; NOFLUC;', 'RFSW; RAD; NORADCOD; NOFLUC;')
    u0 = _vals(_run(nrc(_RING % ('2', _OFF, '0', '0')))[2], 'EM')[3]
    u1 = _vals(_run(nrc(_RING % ('2', _OFF, '0.1', '0.1')))[2], 'EM')[3]
    from pysad import sgeo, emit6d as _e
    i = [k for k, e in enumerate(ffs.els) if e.name == 'SA'][0]
    ent, ext = sgeo.region_frames(ffs.els, i, _e._Phys(E, 91.35))
    p0 = E / ME
    gam = math.sqrt(1 + p0 * p0)
    bz = 2 * 299792458.0 / E
    du = E * (2 / 3) * RE * gam ** 3 * (bz / 2) ** 2 * (ent[1] ** 2 + ent[2] ** 2 + ext[1] ** 2 + ext[2] ** 2) / 0.1
    assert du > 1e-3 and abs((u1 - u0) / du - 1) < 2e-3


def test_ring_tilted_solenoid_zero_field_limit(v2_mode):
    """BZ → 0 のずれて傾いた領域は直進と同じ: Emittance[] と CALC のチューンが、SOL に何も書かないリングと合う
    (傾きの分だけ設計軌道の長さが 1.15·(1/cos − 1) = 7e-6 m 伸びる程度の差)"""
    k0 = lambda s: s.replace('K1=0.1', 'K1=0').replace('K1=-0.1', 'K1=0')       # 領域の中の四極は切る (軸がソレノイドの軸なので、傾けると四極もずれる)
    a = _run(k0(_RING % ('0', _OFF, '0', '0')))[2]
    b = _run(k0(_RING % ('0', '', '0', '0')))[2]
    for u, v in zip(_vals(a, 'NU'), _vals(b, 'NU')):
        assert abs(u - v) < 1e-6
    ea, eb = _vals(a, 'EM'), _vals(b, 'EM')
    for j in (0, 2, 3):                                        # εx、εz、U0 (εy は両方 0)
        assert abs(ea[j] / eb[j] - 1) < 1e-5
    assert abs(ea[1]) < 1e-20 and abs(eb[1]) < 1e-20


def test_geometry_fit_through_tilted_solenoid_region(v2_mode):
    """幾何の fit (GY・GZ) が、ずれて傾いた領域の中の変数 (四極の K1、ドリフトの L) で収束する。領域の中の変数の微分は中心差分
    (K64 の tdgeo は軸方向の伸びだけの式で、軌道の回転の変化が入らない)。fit の後も設計軌道の粒子は出口で 0 に戻る"""
    src = '''MOMENTUM = 3 GEV;
DRIFT DR = (L=1) DR2 = (L=0.5) DV=(L=0.5);
QUAD QF = (L=0.3 K1=0.3) QD = (L=0.3 K1=-0.25);
SOL SA = (L=0 BZ=3 BOUND=1 GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002) SB = (L=0 BZ=0 BOUND=1);
MARK ST = (BX=5 BY=5) M2=();
LINE BL = (ST DR SA DR2 QF DV QD DR2 SB DR M2);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
$FORM="18.14";
PageWidth=1999;
FIT M2 GY 0.0008; FIT M2 GZ -0.0035;
FREE QF QD;
GO;
Print["G1 =",LINE[{"GY","GZ"},"M2"]," ",LINE["K1","QF"]," ",LINE["K1","QD"]];
FIX QF QD; FREE DV;
FIT M2 GZ -0.0035 0;
FIT M2 GY 0.0009;
GO;
Print["G2 =",LINE[{"GY","GZ"},"M2"]," ",LINE["L","DV"]];
Exit;
'''
    ffs, I, out = _run(src)
    assert out.count(' Matched.') >= 3
    g1, g2 = _vals(out, 'G1'), _vals(out, 'G2')
    assert abs(g1[0] - 0.0008) < 1e-8 and abs(g1[1] + 0.0035) < 1e-8 and abs(g1[2] - 0.3) > 1e-3
    assert abs(g2[0] - 0.0009) < 1e-8 and abs(g2[2] - 0.5) > 1e-2
    got = _track(ffs, np.zeros((1, 6)), 3, len(ffs.els))
    assert np.max(np.abs(got[0])) < 1e-12
