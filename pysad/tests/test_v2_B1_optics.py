# -*- coding: utf-8 -*-
"""担当 B1 (v2next): 光学・Emittance の未対応の追加。
(1) Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE (K64 の tbende の tblfre・tbendef1・alcorr、tbend の周長の補正)。gs と全桁一致 (k64)。
(2) K1≠0 の BEND の CALC (tbende) と追跡 (tbendi) の 1e-5 の差: 独立な「厳密な写像」(ローレンツ力の数値積分) との比較で、誤りは追跡の薄い端の写像。v2 で直す。
(3) WSPAC (空間電荷) の Emittance[]・CALC6D・CALC: 一次のラスレットの式との一致、gs との比較 (Emittance[] だけ。gs の CALC・CALC6D は NaN)。
(4) CALC6D・RADCOD つきの CALC と SOL 領域、INTRA と SOL、RADCOD の CALC6D のあとの DP0 (gs と一致)。
(5) RADTAPER: 4 次元の CALC は RADTAPER を見ない (K64 も同じ)、CALC6D は k64 でも RADTAPER を入れる (gs と一致)。
"""
import io
import math
import os
import re

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _ = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _floats(s):
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', s)]


def _tag(out, tag, n=None):
    i = out.index('\n' + tag + ' ') + 1 if ('\n' + tag + ' ') in out else out.index(tag + ' ')
    seg = out[i + len(tag):]
    m = re.search(r'\n(?=[^\s,{}0-9.+-])', seg)
    seg = seg[:m.start()] if m else seg
    v = _floats(seg.replace('{', ' ').replace('}', ' ').replace(',', ' '))
    return v if n is None else v[:n]


def _emit_rule(out, key):
    """Print[e] の出力 ({Stable->..., Tunes->{..}, ...}) から key の数値列"""
    i = out.rindex('{Stable')
    seg = out[i:].replace('\n', ' ')
    m = re.search(key + r'->\s*(\{[^}]*\}|[^,]*)', seg)
    return _floats(m.group(1))


# ====================================================================== (1) BEND の F1・FRINGE
_RING_F1 = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747 %BEND%);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="14.10";
USE RING; CELL; CALC;
%FLAGS%
e = Emittance[];
Print[e];
Exit;
"""

# gs (10-19 ビルド、2026-11-18) の Emittance[] (RFSW RAD NOFLUC): (Tunes, U0 [eV]、EquilibriumPosition、MomentumCompaction)
_GS_F1 = {
    'F1=0.1 FRINGE=3': ([-0.2933271700, -0.1395704757, -0.0179539450], 739673.742271, 0.0282503369, 0.0326541595),
    'F1=0.1 FRINGE=1': ([-0.2933271700, -0.1395704757, -0.0179539450], 739673.742271, 0.0282503369, 0.0326541595),
    'F1=0.1 FRINGE=-1': ([-0.2933269738, -0.1354189066, -0.0179290223], 770787.825581, 0.0294655665, 0.0326544681),
    'F1=0.1 FRINGE=-2': ([-0.2933270813, -0.1372812762, -0.0179420505], 750530.231592, 0.0286739774, 0.0326421595),
    'F1=0.05 FB1=0.07 FB2=0.02 FRINGE=3': ([-0.2933268315, -0.1387818190, -0.0179497834], 746836.300110, 0.0285297923, 0.0326595983),
    'K1=-0.05 F1=0.1 FRINGE=3': ([-0.4307280113, 0.0499557121, -0.0188202545], 759910.329928, 0.0290403357, 0.0359460096),
    'K1=0.04 F1=0.1 FRINGE=1': ([-0.1863156899, -0.3035190683, -0.0173132656], 759910.329928, 0.0290403357, 0.0304198935),
    'F1=0.1 FRINGE=3 K0=1e-3': ([-0.3032263221, -0.1484899393, -0.0179785839], 751039.044037, 0.0286938422, 0.0327767101),
}


def _f1_emit(bend, flags='RFSW; RAD; NOFLUC;'):
    ffs, out = _run(_RING_F1.replace('%BEND%', bend).replace('%FLAGS%', flags))
    return out


@pytest.mark.parametrize('bend', list(_GS_F1))
def test_emittance_bend_f1_matches_gs_k64(k64_mode, bend):
    """k64: Emittance[] の BEND の F1・FRINGE (tblfre・tbendef1・alcorr・周長の補正) が gs と全桁 (表示の 10 桁) で一致する。
    FRINGE = 1・2・3 は同じ (mfring > 0)、−1 は入口側だけ、−2 は出口側だけ、FB1・FB2、K1≠0 (tbendebody の分割)、K0 つき"""
    tun, u0, eq, al = _GS_F1[bend]
    out = _f1_emit(bend)
    t = _emit_rule(out, 'Tunes')
    for a, b in zip(t, tun):
        assert abs(a - b) < 2e-9, (bend, a, b)
    assert abs(_emit_rule(out, 'EnergyLossU0')[0] / u0 - 1) < 1e-9
    assert abs(_emit_rule(out, 'EquilibriumPosition')[0] - eq) < 2e-9
    assert abs(_emit_rule(out, 'MomentumCompaction')[0] / al - 1) < 1e-8


def test_emittance_bend_f1_radcod_and_pol_match_gs_k64(k64_mode):
    """k64: RADCOD つき・POL つきの Emittance[] も gs と一致 (RADCOD では閉軌道が放射で動く)"""
    out = _f1_emit('F1=0.1 FRINGE=3', 'RFSW; RAD; RADCOD; NOFLUC;')
    for a, b in zip(_emit_rule(out, 'Tunes'), [-0.2933113150, -0.1395570999, -0.0179539039]):
        assert abs(a - b) < 2e-9
    assert abs(_emit_rule(out, 'EnergyLossU0')[0] / 739672.942309 - 1) < 1e-9
    out = _f1_emit('K1=-0.05 F1=0.1 FRINGE=3', 'RFSW; RAD; RADCOD; NOFLUC;')
    for a, b in zip(_emit_rule(out, 'Tunes'), [-0.4307138905, 0.0499670036, -0.0188202161]):
        assert abs(a - b) < 2e-9
    out = _f1_emit('F1=0.1 FRINGE=3', 'RFSW; RAD; NOFLUC; POL;')
    assert abs(_emit_rule(out, 'DampingRate')[0] + 1.1182420e-4) < 2e-10        # POL ON は放射の分割が細かい (gs の値)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('bend', ['', 'F1=0.1 FRINGE=3', 'K1=-0.05 F1=0.1 FRINGE=3', 'F1=0.1 FRINGE=-1', 'K1=-0.05 F1=0.05 FB1=0.07 FB2=0.02 FRINGE=3'])
def test_emittance_tunes_equal_calc_tunes_with_bend_f1(mode, bend):
    """両モード: RF・放射なしの Emittance[] のチューン (6D の一周写像) が CALC (4D) のチューンと 1e-9 で一致する (F1 を含む。
    解析的な条件: 放射・RF なしの 6D の x, y のチューンは 4D の光学と同じ)。F1 を無視していた旧版は 0.006 ずれた"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        src = _RING_F1.replace('%BEND%', bend).replace('%FLAGS%', 'NORFSW; NORAD;').replace(
            'e = Emittance[];\nPrint[e];', 'e = Emittance[];\nPrint["T ", Tunes/.e, " ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];')
        _, out = _run(src)
        v = _tag(out, 'T')
        nx, ny = v[3] % 1.0, v[4] % 1.0
        assert abs(abs(v[0]) - min(nx, 1 - nx)) < 1e-9 and abs(abs(v[1]) - min(ny, 1 - ny)) < 1e-9, (bend, v)
    finally:
        compat.set_mode(old)


def test_bend_f1_changes_u0_in_v2(v2_mode):
    """v2: F1 の BEND の U0・チューンが F1 なしと違い、gs の比 (F1 あり/なし = 739673.74/781633.07 = 0.94632) と 2e-3 で合う
    (v2 は端の角の場の中の行路の補正で gs と絶対値が少し違う)"""
    u1 = _emit_rule(_f1_emit('F1=0.1 FRINGE=3'), 'EnergyLossU0')[0]
    u0 = _emit_rule(_f1_emit(''), 'EnergyLossU0')[0]
    assert abs(u1 / u0 / 0.94632 - 1) < 2e-3, (u1, u0)


def _ring_a_f1(flags):
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kek-ls_v03_12.lat')
    if not os.path.exists(path):
        pytest.skip('kek-ls_v03_12.lat がない')
    with open(path, encoding='utf-8') as fh:
        src = fh.read()
    src = src.replace('ANGLE =', 'F1=0.03 FRINGE=3 ANGLE =')
    src += '\nFFS;\n$FORM="14.10";\nUSE RING; CELL; CALC;\n' + flags + '\ne = Emittance[];\nPrint[e];\nExit;\n'
    return _run(src)[1]


def test_large_ring_a_bend_f1_emittance_matches_gs_k64(k64_mode):
    """大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m、全偏向磁石 (K1≠0 の BPI6・BPI7 を含む) に F1=0.03 FRINGE=3 を付けた):
    Emittance[] のチューン・U0・運動量圧縮・減衰率が gs と一致 (F1 なしとは U0 が 3.6%、νy が 5e-3 違う)"""
    out = _ring_a_f1('RFSW; RAD; NOFLUC;')
    for a, b in zip(_emit_rule(out, 'Tunes'), [-0.4199675211, -0.3846720252, -0.0040789829]):
        assert abs(a - b) < 2e-9
    assert abs(_emit_rule(out, 'EnergyLossU0')[0] / 285476.566791 - 1) < 1e-9
    assert abs(_emit_rule(out, 'MomentumCompaction')[0] / 2.30667956e-4 - 1) < 1e-8
    assert abs(_emit_rule(out, 'DampingRate')[0] / -7.4678385e-5 - 1) < 1e-7
    assert abs(_emit_rule(out, 'Emittances')[0] / 1.5180365e-10 - 1) < 1e-7


def test_large_ring_a_bend_f1_v2_consistent(v2_mode):
    """v2: 大きなリング A の F1 つき BEND の Emittance[] (RF・放射なし) のチューンが CALC と 1e-9"""
    out = _ring_a_f1('NORFSW; NORAD;').replace('Print[e];', '')
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kek-ls_v03_12.lat')
    src = open(path, encoding='utf-8').read().replace('ANGLE =', 'F1=0.03 FRINGE=3 ANGLE =')
    src += ('\nFFS;\n$FORM="14.10";\nUSE RING; CELL; CALC;\nNORFSW; NORAD;\ne = Emittance[];\n'
            'Print["T ", Tunes/.e, " ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nExit;\n')
    v = _tag(_run(src)[1], 'T')
    assert abs(abs(v[0]) - min(v[3] % 1, 1 - v[3] % 1)) < 1e-9 and abs(abs(v[1]) - min(v[4] % 1, 1 - v[4] % 1)) < 1e-9, v


def _inject_f1(path, ins, ring, flags):
    src = open(path, encoding='utf-8', errors='replace').read()
    src = re.sub(r'(?<![A-Za-z0-9_])ANGLE\s*=', ins + 'ANGLE =', src)
    return src + '\nFFS;\n$FORM="14.10";\nUSE %s; CELL; CALC;\n%s\ne = Emittance[];\nPrint[e];\nExit;\n' % (ring, flags)


_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.mark.parametrize('name,path,ring,ins,tun,u0', [
    ('pf', 'example/pf/lattice/v9_16_asym.sad', 'RING', 'F1=0.05 FRINGE=3 ', [-0.4024695092, 0.2721266455, -0.0146840676], 392531.84006),
    ('astrid2', 'example/astrid/lattice/ASTRID2_v01_01.sad', 'RING', 'F1=0.04 FB2=0.02 FRINGE=3 ', [0.1799349388, 0.0824834394, -0.0015263854], 6045.54954484),
    ('astrid2_mpw', 'example/astrid/lattice/ASTRID2_v01_02_id.sad', 'RINGMPW', 'F1=0.04 FRINGE=2 ', [0.1799932932, 0.2006305242, -0.001521812], 7704.31936504),
])
def test_other_rings_bend_f1_emittance_matches_gs_k64(k64_mode, name, path, ring, ins, tun, u0):
    """PF リング (2.5 GeV、187 m)、ASTRID2 (0.58 GeV)、ASTRID2 + MPW の全偏向磁石に F1 を付けた Emittance[] が gs と一致"""
    p = os.path.join(_ROOT, path)
    if not os.path.exists(p):
        pytest.skip(path)
    out = _run(_inject_f1(p, ins, ring, 'RFSW; RAD; NOFLUC;'))[1]
    for a, b in zip(_emit_rule(out, 'Tunes'), tun):
        assert abs(a - b) < 2e-9, (name, a, b)
    assert abs(_emit_rule(out, 'EnergyLossU0')[0] / u0 - 1) < 1e-9


def test_emittance_f1_of_steering_bend_is_implemented_by_e(v2_mode):
    """ANGLE=0 の BEND (tsteee) の F1 は v2next E で実装した (以前は NotImplementedError。test_v2_E_steer.py に gs との比較)。F1 を付けると U0 などが変わる"""
    from pysad.emit6d import _tsteee, _State, _Phys
    from pysad.elements import Bend
    res = []
    for f1 in (0.0, 0.1):
        st = _State(np.zeros(6), _Phys(3e9, 10.0))
        el = Bend(name='S', L=0.3, ANGLE=0.0)
        el.K0 = 1e-3
        el.F1 = f1
        el.FRINGE = 3
        _tsteee(st, el, True)
        res.append((st.T.copy(), st.u0))
    assert np.abs(res[0][0] - res[1][0]).max() > 1e-9


# ====================================================================== (2) K1≠0 の BEND の CALC と追跡の差
_L = 1.2
_TH = 0.2617993877991494


def _exact_map(c, k, e1, e2, phi=_TH, L=_L, cq=-0.5):
    """硬い端の厳密な写像 (ローレンツ力の数値積分、DOP853 rtol 1e-13)。磁石の形: 双極磁場 B_y = h (入口・出口の端面の間、端の角 e1, e2 [rad])、
    勾配の磁場 B_y = k x + cq k h x² (x = r − ρ0。中面。cq = −1/2 は K64 の tbendicorr の蹴り −K1 x² h/2 から)、B_r = k y (半径方向、y の 1 次)
    が、半径方向の入口・出口の面の間。b = 軸上の曲率中心側を +r 向きにとる (ANGLE>0 は SAD の座標と同じ x を外向きにとる)。
    戻り (x, px) = 半径方向の出口の面での外向きの座標と運動量。端面の外はドリフト"""
    rho0 = L / phi
    h = 1 / rho0
    x, px, y, py, z, dl = c
    p = 1 + dl
    pz = math.sqrt(p * p - px * px - py * py)
    Pin = np.array([rho0, 0.0, 0.0])
    Pout = np.array([rho0 * math.cos(phi), 0.0, rho0 * math.sin(phi)])
    tin = np.array([0, 0, 1.0])
    rin = np.array([1.0, 0, 0])
    tout = np.array([-math.sin(phi), 0, math.cos(phi)])
    rout = np.array([math.cos(phi), 0, math.sin(phi)])
    n1 = math.cos(e1) * tin - math.sin(e1) * rin
    n2 = math.cos(e2) * tout + math.sin(e2) * rout
    s1 = lambda w: (w[0] - Pin[0]) * n1[0] + (w[2] - Pin[2]) * n1[2]
    s2 = lambda w: (w[0] - Pout[0]) * n2[0] + (w[2] - Pout[2]) * n2[2]
    r1 = lambda w: w[2]
    r2 = lambda w: -w[0] * math.sin(phi) + w[2] * math.cos(phi)
    P = np.array([rho0 + x, y, 0.0])
    pv = np.array([px, py, pz])
    u = pv / np.linalg.norm(pv)
    w = np.concatenate([P - u * 0.4, pv])
    planes = [s1, r1, s2, r2]
    inside_dip = lambda q: s1(q) > 0 and s2(q) < 0
    inside_quad = lambda q: r1(q) > 0 and r2(q) < 0
    for _ in range(12):
        uu = w[3:] / np.linalg.norm(w[3:])
        q = w.copy()
        q[:3] = q[:3] + uu * 1e-9
        dip, quad = inside_dip(q), inside_quad(q)

        def f(l, ww):
            X, Y, Z, pX, pY, pZ = ww
            pp = math.sqrt(pX * pX + pY * pY + pZ * pZ)
            ux, uy, uz = pX / pp, pY / pp, pZ / pp
            rr = math.hypot(X, Z)
            xx = rr - rho0
            By = (h if dip else 0.0) + (k * xx + cq * k * h * xx * xx if quad else 0.0)
            Br = k * Y if quad else 0.0
            Bx, Bz = Br * X / rr, Br * Z / rr
            return [ux, uy, uz, uy * Bz - uz * By, uz * Bx - ux * Bz, ux * By - uy * Bx]
        evs = []
        for pl in planes:
            sg = np.sign(pl(q))

            def mk(pl=pl, sg=sg):
                def e(l, ww):
                    return pl(ww) * sg
                e.terminal = True
                e.direction = -1
                return e
            evs.append(mk())
        sol = solve_ivp(f, [0, 3.0], w, method='DOP853', rtol=1e-13, atol=1e-15, events=evs)
        hit = [j for j in range(4) if len(sol.t_events[j])]
        if not hit:
            w = sol.y[:, -1]
            break
        j = min(hit, key=lambda j: sol.t_events[j][0])
        w = sol.y_events[j][0]
        if s2(w) > 0 and r2(w) > 0 and j in (2, 3):
            break
    u = w[3:] / np.linalg.norm(w[3:])
    lstr = -np.dot(w[:3], tout) / np.dot(u, tout)
    Pf = w[:3] + u * lstr
    rr = math.hypot(Pf[0], Pf[2])
    return np.array([rr - rho0, np.dot(w[3:], rout)])


def _k1_bend(k1, e):
    from pysad import ffscmd
    src = f'''MOMENTUM = 3GEV;
BEND BE = (L={_L} ANGLE={_TH} E1={e} E2={e} K1={k1});
LINE X = (BE);
FFS;
USE X; INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return [x for x in ffs.lat.elements if x.name == 'BE'][0], ffs


def _calc_and_track(el, ffs, c):
    from pysad import lattice, track6n
    cc, M = lattice._bend_element_cod_matrix(el, np.array(c, float))
    t = np.asarray(track6n.track_particles_fast([el], np.array([c], float), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * _L)[0])
    return cc, M, t


_ORBITS = [(1e-3, 0, 0), (3e-3, 0, 0), (0, 5e-4, 0), (2e-3, 3e-4, 0), (1e-3, 0, 5e-3)]


@pytest.mark.parametrize('e', [0.5, 0.25])
def test_k1_bend_calc_map_equals_exact_hard_edge_map(v2_mode, e):
    """CALC (tbende、厳密な tbedge) の写像 (x, px) が硬い端の厳密な写像 (ローレンツ力の積分) と 1e-10 で一致する (軌道 1〜3 mm、δ = 5e-3、
    端の角 e·ANGLE = 0.13・0.065 rad)。v2 でも k64 でも同じ (CALC は変えていない)"""
    el, ffs = _k1_bend(-0.05, e)
    for (x0, px0, dl) in _ORBITS:
        c = [x0, px0, 0, 0, 0, dl]
        cc, _, _ = _calc_and_track(el, ffs, c)
        ex = _exact_map(c, el.K1, e * _TH, e * _TH)
        assert np.abs(cc[:2] - ex).max() < 1e-10, (e, c, cc[:2] - ex)


@pytest.mark.parametrize('e', [0.5, 0.25])
def test_k1_bend_tracking_map_equals_exact_in_v2_not_in_k64(e):
    """追跡 (tbendi) の写像: v2 は厳密な写像と 1e-10 で一致 (CALC と同じ tbedge)。K64 の追跡は薄い端の写像 (入口と出口で px の蹴りと
    x の補正だけ) で、軌道 3 mm で 1.6e-9、δ = 5e-3 + 1 mm で 9e-9 ずれる (CALC は 1e-12 で合うので、誤りは追跡の側)"""
    old = compat._mode
    try:
        compat.set_mode('v2')
        el, ffs = _k1_bend(-0.05, e)
        err2 = []
        for (x0, px0, dl) in _ORBITS:
            c = [x0, px0, 0, 0, 0, dl]
            _, _, t = _calc_and_track(el, ffs, c)
            err2.append(np.abs(t[:2] - _exact_map(c, el.K1, e * _TH, e * _TH)).max())
        assert max(err2) < 1e-10, err2
        compat.set_mode('k64')
        el, ffs = _k1_bend(-0.05, e)
        err64 = []
        for (x0, px0, dl) in _ORBITS:
            c = [x0, px0, 0, 0, 0, dl]
            _, _, t = _calc_and_track(el, ffs, c)
            err64.append(np.abs(t[:2] - _exact_map(c, el.K1, e * _TH, e * _TH)).max())
        if e == 0.5:
            assert err64[1] > 1e-9 and err64[4] > 5e-9 and err64[3] > 3e-9, err64
        assert max(err64) > 5 * max(err2)
    finally:
        compat.set_mode(old)


def test_k1_bend_calc_matrix_vs_tracking_derivative(v2_mode):
    """x0 = 1 mm の行列 (4x4): CALC (tbende の解析行列) と追跡の写像の数値微分の差が、v2 で 5e-8 (K64 の追跡の端は 1.1e-5: 矩形 E=0.5)。
    e=0 (端の角なし) は両モードで 4e-8"""
    from pysad import lattice, track6n
    for e, lim in ((0.5, 6e-8), (0.25, 6e-8), (0.0, 6e-8)):
        el, ffs = _k1_bend(-0.05, e)
        c0 = np.array([1e-3, 0, 0, 0, 0, 3e-3])
        Mc = lattice._bend_element_cod_matrix(el, c0)[1]
        J = np.zeros((4, 4))
        for j in range(4):
            d = np.zeros(6)
            d[j] = 1e-6
            tp = np.asarray(track6n.track_particles_fast([el], np.array([c0 + d]), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * _L)[0])
            tm = np.asarray(track6n.track_particles_fast([el], np.array([c0 - d]), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * _L)[0])
            J[:, j] = (tp[:4] - tm[:4]) / 2e-6
        assert np.abs(Mc[:4, :4] - J).max() < lim, (e, np.abs(Mc[:4, :4] - J).max())


def test_k1_bend_calc_matrix_vs_tracking_derivative_k64(k64_mode):
    """k64: 同じ行列の差が 1.1e-5 (矩形 E=0.5)・2.7e-6 (E=0.25)、E=0 は 4e-8 (記録)"""
    from pysad import lattice, track6n
    res = {}
    for e in (0.5, 0.25, 0.0):
        el, ffs = _k1_bend(-0.05, e)
        c0 = np.array([1e-3, 0, 0, 0, 0, 3e-3])
        Mc = lattice._bend_element_cod_matrix(el, c0)[1]
        J = np.zeros((4, 4))
        for j in range(4):
            d = np.zeros(6)
            d[j] = 1e-6
            tp = np.asarray(track6n.track_particles_fast([el], np.array([c0 + d]), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * _L)[0])
            tm = np.asarray(track6n.track_particles_fast([el], np.array([c0 - d]), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * _L)[0])
            J[:, j] = (tp[:4] - tm[:4]) / 2e-6
        res[e] = np.abs(Mc[:4, :4] - J).max()
    assert 8e-6 < res[0.5] < 1.5e-5 and 2e-6 < res[0.25] < 4e-6 and res[0.0] < 1e-7, res


# ====================================================================== (3) WSPAC (空間電荷)
_WS_SRC = """MOMENTUM = 0.2GEV;
PBUNCH = %(pb)s;
MINCOUP = 0.5;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=0.2e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="14.10";
USE RING; CELL; CALC;
RFSW; RAD; NOFLUC;
e0 = Emittance[];
Print["ZERO ", Tunes/.e0, Emittances/.e0];
WSPAC;
e = Emittance[];
Print["SC ", Tunes/.e, Emittances/.e];
Print["EMITC ", EMITXC, EMITYC, EMITZC];
Exit;
"""


def _ws(pb, mode='k64'):
    old = compat._mode
    compat.set_mode(mode)
    try:
        _, out = _run(_WS_SRC % dict(pb=pb))
    finally:
        compat.set_mode(old)
    return _tag(out, 'ZERO'), _tag(out, 'SC'), _tag(out, 'EMITC')


def test_wspac_emittance_matches_gs_k64(k64_mode):
    """k64: WSPAC ON の Emittance[] (PBUNCH = 1e8、MINCOUP = 0.5、0.2 GeV、12 セルの FODO): 空間電荷のチューンシフトが出る。
    gs の値 (Tunes -.2976763273、-.1419755207、-.0182639185、EMITXC 1.27780330E-9) と全桁一致 (D1: 以前の 4% (x)・3% (y) の差は、K64 の最後の
    スライスのレンズが gammab(nlat+1) = 0 を読んで 8 倍強いことだった。k64 はそれを再現する)。MINCOUP の下限 εy = 0.5 (εx + εy) は gs と同じ"""
    z, s, c = _ws('1E8')
    assert abs(z[0] + 0.2933208697) < 1e-9 and abs(z[1] + 0.1331276964) < 1e-9
    gs = np.array([-0.2976763273, -0.1419755207, -0.0182639185])
    assert np.max(np.abs(np.array(s[:3]) - gs)) < 2e-9, (s[:3], gs)
    assert abs(c[0] / 1.27780330e-9 - 1) < 3e-3
    assert abs(c[1] / c[0] - 0.5) < 1e-9                            # MINCOUP = 0.5: εy = 0.5 εx (下限)


def test_wspac_tune_shift_equals_first_order_laslett_formula(v2_mode):
    """v2: WSPAC の Emittance[] のチューンシフトが、一次のラスレットの式 Δν = −(1/4π) ∮ β_u 2 N r0 /(β²γ³ √(2π) σz σ_u (σx+σy)) ds
    (ガウス分布、ビームの大きさは β ε + η² σδ²。ラスレットの公称の式 Δν = −N r0 C/(4π β²γ³ ε √(2π) σz) の非円形・分散つきへの一般化)
    と 0.15% で一致する。要素を 30 分割したリング (Twiss は分割の端)。PBUNCH = 1e5 (低電流)。シフトは PBUNCH に比例 (1e5・1e7 で比 100 ± 0.3)"""
    from pysad import emit6d, touschek, ffscmd
    n = 30
    src = f"""MOMENTUM = 0.2GEV;
DRIFT D1 = (L=1/{n}) D2 = (L=0.5/{n});
QUAD QF = (L=0.3/{n} K1=0.7/{n}) QD = (L=0.3/{n} K1=-0.6/{n});
BEND B = (L=1.2/{n} ANGLE=0.2617993877991494/{n});
CAVI RFC = (L=0 VOLT=0.2e6 HARM=64);
MARK ST = ();
LINE CELL1 = ({n}*QF {n}*D2 {n}*B {n}*D2 {n}*QD {n}*D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
USE RING; CELL; CALC;
Exit;
"""
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    els = ffs.lat.elements
    p = 0.2e9
    em = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=False)
    shifts = {}
    for pb in (1e5, 1e7):
        t = touschek.touschek_table(ffs.lat, p, em, pb, rfsw=True, intra=False, wspac=True, verbose=False, mincoup=0.5)
        shifts[pb] = (np.array(t['tunes']) - np.array(em['tunes']))[:2]
        if pb == 1e5:
            ex, ey, sige, sigz = t['EMITXC'], t['EMITYC'], t['momentum_spread'], t['bunch_length']
    ph = emit6d._Phys(p, sum(e.L for e in els))
    fac = 2 * 1e5 * emit6d.ELRADI / (ph.p0**2 * ph.h0 * math.sqrt(2 * math.pi) * sigz)
    rows = [ffs.row(k) for k in range(1, len(els) + 2)]
    dx = dy = 0.0
    for i, el in enumerate(els):
        if el.L == 0:
            continue
        r1, r2 = rows[i], rows[i + 1]
        bx, by, eta = (r1.betx + r2.betx) / 2, (r1.bety + r2.bety) / 2, (r1.etx + r2.etx) / 2
        sx, sy = math.sqrt(bx * ex + eta**2 * sige**2), math.sqrt(by * ey)
        w = fac * el.L / (sx + sy)
        dx += -(1 / (4 * math.pi)) * bx * w / sx
        dy += -(1 / (4 * math.pi)) * by * w / sy
    assert abs(shifts[1e5][0] / dx - 1) < 1.5e-3 and abs(shifts[1e5][1] / dy - 1) < 1.5e-3, (shifts, dx, dy)
    assert abs(shifts[1e7][0] / shifts[1e5][0] / 100 - 1) < 3e-3


def test_wspac_calc6d_and_calc_rows_have_the_space_charge_tune_shift_v2(v2_mode):
    """v2: WSPAC ON の CELL CALC6D (と 4 次元の CALC。RFSW ON のとき同じ経路) の Twiss が空間電荷のチューンシフトを持ち、周期的 (BX($$$) = BX(1))、
    Emittance[] のチューンと一致する。K64 (gs) は WSPAC つきの CALC・CALC6D の入口以外の Twiss が NaN になる (k64 モードは空間電荷を入れない)"""
    base = _WS_SRC % dict(pb='1E8')
    src = base.split('e0 = Emittance[];')[0] + '''CALC6D;
CALC;
Print["N0 ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",1], Twiss["BX","$$$"]];
WSPAC;
CALC;
Print["N1 ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",1], Twiss["BX","$$$"]];
e = Emittance[];
Print["TUNES ", Tunes/.e];
NOCALC6D;
CALC;
Print["N2 ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",1], Twiss["BX","$$$"]];
Exit;
'''
    _, out = _run(src)
    n0, n1, tu, n2 = _tag(out, 'N0', 4), _tag(out, 'N1', 4), _tag(out, 'TUNES', 3), _tag(out, 'N2', 4)
    assert abs(n1[3] / n1[2] - 1) < 1e-9 and abs(n2[3] / n2[2] - 1) < 1e-9          # 周期的
    assert n1[0] < n0[0] - 3e-3 and n1[1] < n0[1] - 7e-3                              # 空間電荷でチューンが下がる (PBUNCH = 1e8: −4.1e-3、−8.5e-3)
    assert abs(n1[0] - (3 - abs(tu[0]))) < 2e-6 and abs(n1[1] - (2 - abs(tu[1]))) < 2e-6          # Emittance[] のチューンと一致
    assert abs(n2[0] - n1[0]) < 2e-3 and abs(n2[1] - n1[1]) < 2e-3                    # 4 次元の経路も空間電荷つき (RF の焦点の差だけ)


def test_wspac_calc6d_k64_keeps_ignoring_space_charge(k64_mode):
    """k64: WSPAC ありの CALC6D は今までどおり空間電荷を入れない (gs は NaN を返す K64 の誤り)"""
    src = (_WS_SRC % dict(pb='1E8')).split('e0 = Emittance[];')[0] + 'CALC6D; CALC;\nWSPAC; CALC;\nPrint["N ", Twiss["NX","$$$"]/(2Pi)];\nExit;\n'
    _, out = _run(src)
    assert abs(_tag(out, 'N')[0] - 2.7066791303) < 1e-8


# ====================================================================== (4) SOL 領域: CALC6D・RADCOD の CALC・INTRA
_SOL_RING = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.25);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QFE = (L=0.3 K1=0.7 DX=3e-4 DY=-2e-4) QS = (L=0.2 K1=0.05 ROTATE=0.3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
SOL SA = (L=0 BZ=2 BOUND=1 GEO=1) SM = (L=0 BZ=0.8) SB = (L=0 BZ=0 BOUND=1);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELLS = (QFE D2 B D2 QD D3 SA D3 QS D3 SB D3);
LINE RING = (ST 11*CELL1 CELLS RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
'''


def test_sol_region_calc6d_radcod_calc_match_gs_k64(k64_mode):
    """k64: SOL 領域のあるリング (3 GeV、BZ = 2 T の領域に回転四極) の CALC (4D)、CALC6D、RADCOD つき CALC6D、RADCOD の CALC6D の後の 4D CALC
    (DP0 は LINE["DDP",1] の別名なので閉軌道の δ で光学を計算する) のチューン・β が gs と 1e-9 (SOL の数値微分の精度) で一致する。
    以前は SOL があると SADError だった"""
    src = _SOL_RING + '''Print["C4D ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",{1,30,70}]];
RFSW; NORAD; NOFLUC; CALC6D;
CALC;
Print["C6D ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["NZ","$$$"]/(2Pi), Twiss["BX",{1,30,70}]];
RAD; RADCOD;
CALC;
Print["C6DRC ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["NZ","$$$"]/(2Pi), Twiss["BX",{1,30,70}], Twiss["DX",1], Twiss["DDP",1]];
NOCALC6D;
CALC;
Print["C4DRC ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",{1,30,70}]];
Exit;
'''
    _, out = _run(src)
    ref = {
        'C4D': [2.73475252978800, 1.87531224259842, 5.18402239147695, 1.50681636461491, 4.03654103677484],
        'C6D': [2.73495437251036, 1.87554377922481, -0.01173531080122, 5.18451031535258, 1.50628681331677, 4.03528365708710],
        'C6DRC': [2.73497629185364, 1.87556637674447, -0.01173482495696, 5.18470402895559, 1.50673394160417, 4.03628835938675, 4.219482453095e-4, 3.560338125976e-5],
        'C4DRC': [2.73465818826401, 1.87522109782155, 5.18340468531468, 1.50681743651318, 4.03674570416948],
    }
    for tag, v in ref.items():
        got = _tag(out, tag, len(v))
        assert len(got) == len(v)
        for a, b in zip(got, v):
            assert abs(a - b) < 3e-9 * max(1.0, abs(b)) + 1.2e-9, (tag, a, b)


def test_calc6d_radcod_sets_dp0_for_the_following_4d_calc_k64(k64_mode):
    """k64 (SOL なしのリング、gs と同じ): CELL CALC6D (RADCOD) のあと DP0 は先頭の DDP (閉軌道の δ = 1.2868e-4) になり、次の 4 次元の CALC の
    チューンがそのエネルギーで変わる (gs: 2.70642542976106)。DP0 = 0.001 を与えればそのエネルギー、DP0 = 0 に戻せば元に戻る"""
    src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
USE RING; CELL; CALC;
RFSW; RAD; RADCOD; CALC6D; CALC;
Print["C6DRC ", Twiss["NX","$$$"]/(2Pi), Twiss["DDP",1], DP0];
NOCALC6D; NORADCOD; CALC;
Print["A ", Twiss["NX","$$$"]/(2Pi), Twiss["DDP",{1,5}]];
RADCOD; DP0=0.001; CALC;
Print["B ", Twiss["NX","$$$"]/(2Pi), Twiss["DDP",{1,5}]];
DP0=0; NOCALC6D; CALC;
Print["C ", Twiss["NX","$$$"]/(2Pi), Twiss["DDP",{1,5}]];
Exit;
'''
    _, out = _run(src)
    assert abs(_tag(out, 'C6DRC', 3)[0] - 2.70669909844398) < 1e-9 and abs(_tag(out, 'C6DRC', 3)[2] - 1.286841123935e-4) < 1e-12
    assert abs(_tag(out, 'A')[0] - 2.70642542976106) < 1e-9
    assert abs(_tag(out, 'B')[0] - 2.70418147848253) < 1e-9
    assert abs(_tag(out, 'C')[0] - 2.70675785170376) < 1e-9


def test_sol_region_intra_emittance_matches_gs_k64(k64_mode):
    """k64: SOL 領域のあるリングの INTRA ON の Emittance[] (IBS の反復。領域内の要素は solemit) が gs と一致 (εx 3e-9、εz 1e-7)。
    以前は NotImplementedError。TouschekLifetime は gs・pysad とも NaN (このリングの条件)"""
    src = _SOL_RING.replace('MOMENTUM = 3GEV;', 'MOMENTUM = 3GEV;\nPBUNCH = 1E10;') + '''RFSW; RAD; NOFLUC; INTRA;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["EMC ", EMITXC, EMITYC, EMITZC];
Exit;
'''
    _, out = _run(src)
    ref = [2.706114556595e-7, 2.245979910482e-8, 2.778701736555e-5]
    for a, b in zip(_tag(out, 'EMT'), ref):
        assert abs(a / b - 1) < 2e-7, (a, b)
    for a, b in zip(_tag(out, 'EMC'), ref):
        assert abs(a / b - 1) < 2e-7


# ====================================================================== (5) RADTAPER
_TAPER = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
USE RING; CELL; CALC;
Print["C4D ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",{1,30,70}]];
RFSW; RAD; RADCOD; RADTAPER; CALC;
Print["C4DTAPER ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["BX",{1,30,70}], Twiss["DDP",{1,30,70}]];
CALC6D; CALC;
Print["C6DTAPER ", Twiss["NX","$$$"]/(2Pi), Twiss["NY","$$$"]/(2Pi), Twiss["DDP",{1,30,70}]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_calc_4d_ignores_radtaper_like_k64(mode):
    """両モード (gs も同じ): 4 次元の CALC は RADTAPER を見ない (K64 の qtwiss には rtaper がなく、テーパーは 6D の tturne・qtwissfrac だけ)。
    RFSW・RAD・RADCOD・RADTAPER ON でも、チューン・β・DDP が RADTAPER なしの CALC と同じ。実装しない理由: 4 次元の光学には放射による
    エネルギー変化がないのでテーパー (局所の運動量に磁石の強さを合わせる) は何も変えない。テーパーの光学は CALC6D"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        _, out = _run(_TAPER)
        a, b = _tag(out, 'C4D', 5), _tag(out, 'C4DTAPER', 8)
        assert abs(a[0] - 2.70675785189057) < 1e-9
        assert np.abs(np.array(a) - np.array(b[:len(a)])).max() < 1e-9
        assert np.abs(b[len(a):]).max() < 1e-15
    finally:
        compat.set_mode(old)


def test_calc6d_radtaper_matches_gs_k64(k64_mode):
    """k64: CALC6D + RADTAPER が RADTAPER を入れる (Emittance[] と同じ閉軌道・係数。gs で確かめた: チューンが Emittance[] の値 2.7067335)。
    gs 2.70673351546438・DDP {0, −1.08537004e-4, −2.38719232e-4} と、gs の Emittance[] との差と同じ 1e-7 で一致 (以前はテーパーなし 2.70669910)"""
    _, out = _run(_TAPER)
    v = _tag(out, 'C6DTAPER', 5)
    assert abs(v[0] - 2.70673351546438) < 2e-7 and abs(v[1] - 1.86681334434236) < 2e-7
    assert abs(v[3] + 1.08537004e-4) < 2e-9 and abs(v[4] + 2.38719232e-4) < 2e-9
