# -*- coding: utf-8 -*-
"""v2next 担当 E 項目 1: SOL 領域 (ずれて傾いたソレノイドを含む) の Emittance[] の POL の量 (スピン回転・n0・SpinTune・スピン軌道結合 D)。

K64 の tradke は、運動量の向きの磁場の積分 (縦磁場 bz·al、端の ψ) を z 軸の向きの磁場として足す。運動量が z 軸から傾くと (傾いて・ずれたソレノイドの
軌道、軸を外れた粒子)、(1) 向きの変化から作る B⊥ に垂直成分が二重に入り、(2) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない。傾き 3〜10 mrad・
bz = 2.4 1/m で SpinTune が 1.5e-4〜2.6e-4、n0 が 4e-3 rad 違った。v2 は追跡 (track6._sprot の bpar) と同じく B⊥ と B∥ を分けて足す
(solemit._trk・spin.tradke_spin の hint)。あわせて直したもの (gate なし、pysad 自身の不具合): (a) 軌道ゼロの四極の Emittance[] (POL) が tradke を通さずスピン軌道結合
D が 0 だった、(b) 領域内の四極の ROTATE・K1<0 でスピンの座標系の回転 (tsolrote の srot) が抜けていた (Emittance[] と track6 の両方)、(c) 領域内の CAVI を
追跡 (v2) が放射・スピンの刻みに通さなかった。v2 だけ: 領域内の DRIFT・QUAD の放射・スピンの刻み (bz·L/n ≤ 0.02)、向きの変化の大きさを θ にする
(spin.SPIN_ANGLE_EXACT、N8)。

独立な根拠: (1) 回転座標系の定常な歳差の閉じた式 (analytic_region)、(2) Thomas-BMT の数値積分 (test_v2_C2_spin._tbmt)、(3) 追跡のスピン (track6)、(4) gs (k64)。
"""
import io
import math
import re

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, emit6d, sgeo, spin, track6, track6n as _n, track6n_ext as _x

GS = 0.001159652181280002
ME = 0.51099895e6
CL = 299792458.0


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _ = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _track(ffs, c, k0, k1, spin_=None):
    """POL + RAD (FLUC なし) のスピン追跡。戻り (座標 (n, 6)、スピン (n, 3))"""
    c = np.asarray(c, float)
    n = len(c)
    cols = [list(c[:, j]) for j in range(6)]
    ffs.flags['RAD'] = True
    ffs.flags['FLUC'] = False
    ffs.flags['POL'] = True
    try:
        sp = np.asarray(spin_, float)
        out = ffs.fn_track([[float(k0), cols + [[math.atan2(s[2], s[0]) for s in sp], [s[1] for s in sp], [1.0] * n]], float(k1)])
    finally:
        ffs.flags['POL'] = False
    r = out[1]
    so = [[math.sqrt(max(0.0, 1 - y * y)) * math.cos(a), y, math.sqrt(max(0.0, 1 - y * y)) * math.sin(a)] for a, y in zip(r[6], r[7])]
    return np.array(r[:6]).T, np.array(so)


def _rot(w):
    th = np.linalg.norm(w)
    if th < 1e-300:
        return np.eye(3)
    k = w / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * K @ K


# ---------------------------------------------------------------- (1) 傾いて・ずれたソレノイドの 1 領域: 閉じた式
def _region(bzT, p, geo, L=0.5):
    src = f'''MOMENTUM = {p} GEV;
DRIFT E1 = (L={L});
SOL SA = (L=0 BZ={bzT} BOUND=1 GEO=1 DISFRIN=1 {geo}) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA E1 SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _out = _run(src)
    els = list(ffs.lat.elements)
    i = [k for k, e in enumerate(els) if e.name == 'SA'][0]
    ent, ext = sgeo.region_frames(els, i, emit6d._Phys(ffs.momentum(), 10.0))
    return ffs, els, i, ent[0], ext[0]


def _after_entry(ffs, els, i):
    """設計軌道の粒子の、SOL の入口の端 (trots + tsfrin) の後の座標 c (ソレノイド系、正準運動量) と領域の bz"""
    ph = emit6d._Phys(ffs.momentum(), 100.0)

    class Col:
        def __init__(s):
            s.items = []

        def _add(s, k, p, name, cnt, sm):
            s.items.append((k, p))
    col = Col()
    _x.compile_sol_region(col, els, i, ph)
    kind, pr = col.items[0]
    assert kind == _x.K_SOL_EDGE
    (_rot_, r, dx, dy, dz, dbz, bz0, bz1, f1, _radok) = pr
    S = _n.make_state(np.zeros((1, 6)), ph)
    _x.trots1(0, S, r, dx, dy, dz, True)
    _x.tsfrin1(0, S, dbz)
    return _n.state_out(S, 1)[0], bz1


def analytic_region(ffs, els, i, rent, rext, bzT, p, L):
    """設計軌道の粒子のスピンの 1 領域の写像 (出口の座標系の成分) = r_出口 · M · r_入口。M はソレノイド系の一様な縦磁場 B = bz ẑ の中の
    歳差: 運動量 p̂ は B のまわりを角速度 −B で回り (路程 L/p̂_z)、スピンは Ω = (1 + aγ) B − (aγ − a)(B·p̂) p̂ (Thomas-BMT)。
    p̂ と一緒に回る系では Ω_rot = aγ B − (aγ − a)(B·p̂0) p̂0 が定常なので、M = R(−B l) R(−Ω_rot l)。p̂0 は入口の端の後の力学的運動量の向き"""
    gam = math.sqrt(1 + (p * 1e9 / ME) ** 2)
    c, bz = _after_entry(ffs, els, i)
    pm = np.array([c[1] + 0.5 * bz * c[2], c[3] - 0.5 * bz * c[0], 0.0])
    pm[2] = math.sqrt(1 - pm[0] ** 2 - pm[1] ** 2)
    Bn = bz * np.array([0, 0, 1.0])
    Orot = GS * gam * Bn - (GS * gam - GS) * (Bn @ pm) * pm
    lp = L / pm[2]
    return rext @ _rot(-Bn * lp) @ _rot(-Orot * lp) @ rent


_GEOS = [(2.0, 0.5, 'DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002'), (2.0, 0.5, 'DPX=0.02'), (1.0, 1.0, 'DPY=0.05'),
         (3.0, 0.3, 'DPX=0.03 DPY=-0.02'), (2.0, 0.5, 'DX=0.004 DY=0.003 DPX=0.001')]


@pytest.mark.parametrize('bzT,p,geo', _GEOS)
def test_tilted_solenoid_spin_matrix_equals_closed_form_in_tracking_and_emittance_v2(v2_mode, bzT, p, geo):
    """ずれて傾いたソレノイド 1 台 (長さ 0.5 m、傾き 1〜50 mrad、1〜3 T) を通る設計軌道のスピンの回転行列が、閉じた式 (analytic_region。回転した縦磁場の中の
    歳差で、傾き θ のとき磁場の横成分 bz sin θ が (1 + aγ) の歳差を与える) と、追跡 (track6、5e-7 以下 = 放射の刻みの粗さ) で一致し、Emittance[] の srot (1e-9) は追跡と一致する"""
    ffs, els, i, rent, rext = _region(bzT, p, geo)
    n = len(ffs.els) + 1
    Mtr = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _b, s = _track(ffs, np.zeros((1, 6)), 1, n, spin_=[s0])
        Mtr[:, j] = s[0]
    Man = analytic_region(ffs, els, i, rent, rext, bzT, p, 0.5)
    ph = emit6d._Phys(ffs.momentum(), 100.0)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne(els[i:], np.zeros(6), ph, rfsw=False, radcod=False)
    finally:
        emit6d.POL_MODE = old
    Rem = st.srot[:, :3]
    assert np.abs(Man - np.eye(3)).max() > 0.1                 # 歳差は大きい (bz L = 0.1〜0.9 rad)
    assert np.abs(Mtr - Man).max() < 1e-6, np.abs(Mtr - Man).max()
    assert np.abs(Rem - Mtr).max() < 1e-8, np.abs(Rem - Mtr).max()


# ---------------------------------------------------------------- (2) リング: SpinTune・n0
def _ring(sol, region='DS DS', p=0.5):
    return f'''MOMENTUM = {p} GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
SOL SA = (L=0 BZ={sol['bz']} BOUND=1 GEO=1 {sol.get('geo', '')}) SB = (L=0 BZ=0 BOUND=1);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA {region} SB);
LINE RINGQ = (ORIGIN 12*CELLQ CELLS 11*CELLQ);
FFS;
USE RINGQ;
CELL; CALC;
NORFSW; NORAD; NOFLUC; POL;
Exit;
'''


def _ring_compare(sol, region='DS DS', p=0.5):
    ffs, _out = _run(_ring(sol, region, p))
    ffs.run_emit(quiet=True)
    sp = ffs.emit['spin']
    n = len(ffs.els) + 1
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _b, s = _track(ffs, np.zeros((1, 6)), 1, n, spin_=[s0])
        M[:, j] = s[0]
    cmu = (np.trace(M) - 1) / 2
    nu_tr = math.acos(max(-1.0, min(1.0, cmu))) / (2 * math.pi)
    w, v = np.linalg.eig(M)
    k = int(np.argmin(abs(w - 1)))
    n0 = np.real(v[:, k])
    n0 /= np.linalg.norm(n0)
    n_em = sp['sps'][:, 0]
    if n0 @ n_em < 0:
        n0 = -n0
    return dict(nu_em=abs(sp['nup']), nu_tr=nu_tr, ang=math.acos(min(1.0, n0 @ n_em)), srot1=sp['srot1'], M=M, sp=sp, ffs=ffs)


_SOLS = [('axis', dict(bz=2)),
         ('tilt', dict(bz=2, geo='DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002')),
         ('tilt-big', dict(bz=2, geo='DX=0.004 DY=-0.002 DPX=0.01 DPY=0.006'))]


@pytest.mark.parametrize('name,sol', _SOLS)
def test_ring_spin_tune_and_n0_with_tilted_solenoid_v2(v2_mode, name, sol):
    """平面リング (0.5 GeV、aγ = 1.135、24 台の偏向磁石) の 12 セル目に長さ 1 m・BZ = 2 T (bz = 1.2 1/m、スピンの歳差 2.4 rad) のソレノイドを置く。
    Emittance[] の SpinTune が追跡の 1 周のスピン行列の回転角と 1e-6 以内 (残りは追跡の放射のエネルギー損失 δ = −2.4e-6 による 5e-7)、n0 の向きが 1e-5 rad 以内。
    傾き・ずれがあっても軸上と同じ精度 (K64 の式では傾きがあると SpinTune が 1.5e-4〜2.6e-4、n0 が 4e-3 rad 違った。次の試験)"""
    r = _ring_compare(sol)
    assert abs(r['nu_em'] - r['nu_tr']) < 1e-6, (r['nu_em'], r['nu_tr'])
    assert r['ang'] < 1e-5, r['ang']
    assert np.abs(r['srot1'] - r['M']).max() < 1e-5


def test_ring_with_tilted_solenoid_k64_formula_is_off_in_spin_tune_and_n0(v2_mode):
    """旧来の式 (spin.SOL_BPAR = False: bsi を z 軸の向きに足す。K64 と同じ) では、傾いたソレノイドで SpinTune が 1.5e-4 以上・n0 が 1e-3 rad 以上ずれ、
    軸上では 5e-7 のまま (傾き・ずれが原因)。向きの変化の大きさを sin θ のままにする (SPIN_ANGLE_EXACT = False、N8) と、平面のリングで 1.3e-5 ずれる"""
    old = spin.SOL_BPAR
    spin.SOL_BPAR = False
    try:
        ra = _ring_compare(_SOLS[0][1])
        rt = _ring_compare(_SOLS[1][1])
    finally:
        spin.SOL_BPAR = old
    assert abs(ra['nu_em'] - ra['nu_tr']) < 1e-6
    assert abs(rt['nu_em'] - rt['nu_tr']) > 1.5e-4
    assert rt['ang'] > 1e-3
    old = spin.SPIN_ANGLE_EXACT
    spin.SPIN_ANGLE_EXACT = False
    try:
        r0 = _ring_compare(dict(bz=0))
    finally:
        spin.SPIN_ANGLE_EXACT = old
    assert 1e-5 < abs(r0['nu_em'] - r0['nu_tr']) < 5e-5


# ---------------------------------------------------------------- (3) 要素ごとのスピン軌道結合 D: Thomas-BMT の積分との比較
_BASE = '''MOMENTUM = %g GEV;
%s
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1 DISFRIN=1) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA E1 SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''


def _tbmt(c, s, L, field, p0mc):
    """運動方程式と Thomas-BMT の数値積分。c = (x, 力学的 px, y, 力学的 py, z, δ)、field(x, y, s) → B (1/m、p0 で規格化)"""
    pt = 1 + c[5]
    pv = np.array([c[1], c[3], math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)])
    gam = math.sqrt(1 + (p0mc * pt) ** 2)

    def rhs(l, u):
        r, p, S = u[:3], u[3:6], u[6:]
        B = field(r[0], r[1], r[2])
        ph = p / pt
        Bl = np.dot(B, ph) * ph
        Om = ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt
        return np.concatenate([ph, np.cross(p, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: u[2] - L
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2 * L + 1), np.concatenate([[c[0], c[2], 0.0], pv, s]), events=ev, rtol=1e-12, atol=1e-15, method='DOP853')
    return sol.y_events[0][0]


def _mult_field(kn, th):
    """MULT の横磁場 (p0 で規格化): By + i Bx = Σ kn[n] z^n/n! (z = 回した局所座標、kn = 積分強度/L)。実験室系のベクトルに戻す"""
    ct, st = math.cos(th), math.sin(th)

    def mk(bz):
        def f(x, y, s):
            z = complex(ct * x - st * y, st * x + ct * y)
            F = 0j
            zn = 1 + 0j
            fact = 1.0
            for n, k in enumerate(kn):
                F += k * zn / fact
                zn *= z
                fact *= (n + 1)
            return np.array([ct * F.imag + st * F.real, -st * F.imag + ct * F.real, bz])
        return f
    return mk


def _quad_field(k, th=0.0):
    ct, st = math.cos(th), -math.sin(th)

    def mk(bz):
        def f(x, y, s):
            xq, yq = ct * x + st * y, -st * x + ct * y
            bxq, byq = k * yq, k * xq
            return np.array([ct * bxq - st * byq, st * bxq + ct * byq, bz])
        return f
    return mk


_L = 0.4
_ELEMENTS = [
    ('drift', 'DRIFT E1 = (L=1.0);', 1.0, lambda bz: (lambda x, y, s: np.array([0.0, 0.0, bz])), 1.0),
    ('bend', 'BEND E1 = (L=0.4 K0=0.02 ROTATE=0.3);', _L, lambda bz: (lambda x, y, s: np.array([0.05 * math.sin(0.3), 0.05 * math.cos(0.3), bz])), _L),
    ('quad', 'QUAD E1 = (L=0.4 K1=0.3 DISFRIN=1);', _L, _quad_field(0.75), _L),
    ('quad_k1neg', 'QUAD E1 = (L=0.4 K1=-0.3 DISFRIN=1);', _L, _quad_field(-0.75), _L),
    ('quad_rot', 'QUAD E1 = (L=0.4 K1=0.3 ROTATE=0.4 DISFRIN=1);', _L, _quad_field(0.75, 0.4), _L),
    ('mult', 'MULT E1 = (L=0.4 K0=0.02 K1=0.3 DISFRIN=1);', _L, _mult_field([0.02 / _L, 0.3 / _L], 0.0), _L),
    ('mult_rot', 'MULT E1 = (L=0.4 K0=0.01 K1=0.3 ROTATE=0.5 DISFRIN=1);', _L, _mult_field([0.01 / _L, 0.3 / _L], 0.5), _L),
    ('cavi', 'CAVI E1 = (L=0.4 VOLT=0 FREQ=5e8);', _L, lambda bz: (lambda x, y, s: np.array([0.0, 0.0, bz])), _L),
    ('mult_rf', 'MULT E1 = (L=0.4 K0=0.01 K1=0.3 VOLT=1e3 FREQ=5e8 DISFRIN=1);', _L, _mult_field([0.01 / _L, 0.3 / _L], 0.0), _L),
]


def _study(el, fieldf, L, bzT, p, du=1e-4):
    ffs, _o = _run(_BASE % (p, el, bzT))
    els = list(ffs.lat.elements)
    i0 = [k for k, e in enumerate(els) if type(e).__name__ == 'Solenoid'][0]
    pm = ffs.momentum()
    ph = emit6d._Phys(pm, 100.0)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne(els[i0:], np.zeros(6), ph, rfsw=False, radcod=False)
    finally:
        emit6d.POL_MODE = old
    R = st.srot[:, :3].copy()
    D = st.srot[:, 3:9].copy()

    def tr(u, s):
        _o, sp = track6.track_particles(els[i0:], np.array([u]), 1, pm, 0.0, False, True, False, spins=[s], radpol=False, circ=100.0)
        return sp[0]

    def exact(c, s):
        bz = bzT * CL / (p * 1e9)
        ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
        return _tbmt(ck, s, L, fieldf(bz), p * 1e9 / ME)[6:]
    s0 = np.array([0.3, 0.8, 0.52])
    s0 /= np.linalg.norm(s0)
    rows = []
    for j in range(4):
        u = np.zeros(6)
        u[j] = du
        dx_ = (exact(u, s0) - exact(-u, s0)) / 2
        dt = (tr(u, s0) - tr(-u, s0)) / 2
        dm = ((np.eye(3) - spin._cross_mat(D @ u)) @ (R @ s0) - (np.eye(3) + spin._cross_mat(D @ u)) @ (R @ s0)) / 2
        rows.append((dx_, dt, dm))
    r0 = exact(np.zeros(6), np.array([1.0, 0, 0]))
    return rows, R, r0


@pytest.mark.parametrize('name,el,L,fieldf,_l', _ELEMENTS)
def test_solenoid_region_spin_orbit_coupling_emittance_equals_thomas_bmt_v2(v2_mode, name, el, L, fieldf, _l):
    """領域内の DRIFT・BEND・QUAD (K1<0・ROTATE つき)・MULT (K0 + K1、ROTATE つき)・CAVI (VOLT=0) の 1 要素を、強いソレノイド (BZ = 2 T、0.5 GeV、bz = 1.2 1/m) の
    中を通る粒子のスピンの軌道座標についての微分 (x, px, y, py の中心差分 1e-4): Thomas-BMT の数値積分 (厳密) に対し、追跡 (track6) は 1e-6 以内 (最大の成分の 2%。MULT の蹴りの分割)、
    Emittance[] の D (srot(:,4:9)) は 3e-6 + 最大の成分の 1% 以内 (放射の刻み bz·L/n ≤ 0.02 の弦の誤差と四極・MULT の分割)。軸上のスピン回転 R は厳密解と 3e-5 以内 (MULT の K0 の蹴りの分割が最大)"""
    rows, R, r0 = _study(el, fieldf, L, 2.0, 0.5)
    big = max(np.abs(r[0]).max() for r in rows)
    assert big > 1e-5
    for dx_, dt, dm in rows:
        assert np.abs(dt - dx_).max() < 1e-6 + 0.01 * big, (name, dt, dx_)
        assert np.abs(dm - dx_).max() < 3e-6 + 0.01 * big, (name, dm, dx_)
    assert np.abs(R @ np.array([1.0, 0, 0]) - r0).max() < 3e-5


def test_solenoid_region_quad_spin_frame_rotation_k1neg_and_roll_v2(v2_mode):
    """領域内の四極 (K1<0 は −90°、ROTATE は座標を回す) は、座標と一緒にスピンの座標系も回す (K64 の tsolrote の srot)。追跡 (track6) の抜けを直した (気づかずに
    2e-3 ずれていた): 軸を外れた粒子のスピンが Thomas-BMT の積分と 3e-6 以内。k64 モードの追跡は領域内の四極の bsi が K64 のままで (v2 の B80 の項) 比べられない。
    k64 の Emittance[] の srot の回転は次の gs との比較で確かめる"""
    for el, fieldf in (('QUAD E1 = (L=0.4 K1=-0.3 DISFRIN=1);', _quad_field(-0.75)), ('QUAD E1 = (L=0.4 K1=0.3 ROTATE=0.4 DISFRIN=1);', _quad_field(0.75, 0.4))):
        ffs, _o = _run(_BASE % (1.0, el, 1.0))
        bz = CL / 1e9
        worst = 0.0
        for s0 in ([0, 0, 1], [1, 0, 0], [0, 1, 0]):
            for c in ([1e-3, 2e-4, -1e-3, 1e-4, 0, 1e-3], [-2e-3, -3e-4, 5e-4, 4e-4, 0, -5e-3]):
                _b, s = _track(ffs, np.array([c]), 1, len(ffs.els) + 1, spin_=[s0])
                ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
                u = _tbmt(ck, s0, 0.4, fieldf(bz), 1e9 / ME)
                worst = max(worst, float(np.abs(s[0] - u[6:]).max()))
        assert worst < 3e-6, (el, worst)


# ---------------------------------------------------------------- (4) gs (k64): 軌道ゼロの四極のスピン軌道結合、SOL 内の四極のスピンの回転
_GS_PERFECT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QS = (L=0.2 K1=0.05 ROTATE=0.7853981633974483);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QF D2 B D2 QD QS D1);
LINE RING = (ST 11*CELL1 CELL2 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''

# gs の 'One-turn depolarization vectors' (x px y py z pz の欄、gx gy gz の行)。表示は 7 桁
_GS_PERFECT_D = np.array([[0.3243963, -1.059973, -2.498274, 16.24751, 0.0, 0.2222216],
                          [-8.480186, 26.12038, 0.4333444, 1.369173, 0.0, 20.83216],
                          [0.0, 0.0, 6.845606, 7.107342, 0.0, 0.0]])

_GS_SOL = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) DS = (L=0.3);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QS1 = (L=0.2 K1=-0.3) QS2 = (L=0.2 K1=0.3 ROTATE=0.5);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
SOL SA = (L=0 BZ=1 BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QF D2 B D2 QD SA DS QS1 DS QS2 DS SB D1);
LINE RING = (ST 11*CELL1 CELL2 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''

_GS_SOL_D = np.array([[0.06968153, -0.8052360, -0.9971148, 20.86659, 0.0, 3.710627],
                      [-8.714483, 27.47604, -1.080052, -3.592943, 0.0, 20.12429],
                      [0.0, 0.0, 6.845606, 7.107342, 0.0, 0.0]])
_GS_SOL_NU = 0.40428183
_GS_PERFECT_NU = 0.40407272


@pytest.mark.parametrize('src,ref,nu', [(_GS_PERFECT, _GS_PERFECT_D, _GS_PERFECT_NU), (_GS_SOL, _GS_SOL_D, _GS_SOL_NU)])
def test_one_turn_depolarization_vectors_match_gs_k64(k64_mode, src, ref, nu):
    """k64: Emittance[] (POL) の 'One-turn depolarization vectors' (スピン軌道結合 D = srot(:,4:9)) と SpinTune が gs と全桁 (表示 7 桁) 一致する。
    (a) 全部の四極の軌道が 0 の完全なリング (斜め四極つき): K64 の tquade は軌道ゼロでも tradke を通すので、四極の D が 0 でない
    (旧 pysad は軌道ゼロの四極を恒等にして D の x・y 列が違った)。(b) SOL 領域 (BZ = 1 T) の中に K1<0 の四極と ROTATE つきの四極: tsolrote がスピンの座標系を回す
    (旧 pysad は回さなかった。srot の x 列が 1.7 違った)"""
    ffs, _o = _run(src)
    ffs.run_emit(quiet=True)
    sp = ffs.emit['spin']
    D = sp['srot'][:, 3:9]
    assert np.abs(D - ref).max() < 3e-6 * max(1.0, np.abs(ref).max()), np.abs(D - ref).max()
    assert abs(abs(sp['nup']) - nu) < 2e-8
