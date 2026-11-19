# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): スピン (POL) と SOL 内の放射 — B39 (tradke の dbeam(8))、B58・B64 (端の縦磁場の積分 bsir0)、B55 (平衡偏極度)。

根拠:
- B39: 放射の拡散は力学的運動量 (px + b y, py − b x) の向きに入る。その共分散は A B Aᵀ (A = [[0,1,b,0],[−b,0,0,1]]) で、
  px・py の成分は B[py,px] − b (B[x,px] − B[y,py] + b B[x,y])。K64 は B[y,py] の代わりに B[py,py] (次元が違う) を使う。
- B58・B64: 磁石の端の縦磁場の積分は ∫Bs ds/(Bρ) = ψ(x, y) (横磁場のスカラーポテンシャル。偏向 y/ρ、四極 k x y、
  多極 Im Σ k_n z^(n+1)/(n+1)!)。入口で +ψ、出口で −ψ。Emittance[] (POL) の 1 要素のスピン回転が、別のコードである
  追跡 (track6 の sprot。tbend・tsolqu・tmulti はこの ψ を使う) と一致する。
- B55: exp(cb²)·erfc(cb) は Faddeeva 関数 w(i·cb)。平面リングの平衡偏極度は Sokolov–Ternov の 8/(5√3) = 92.38%。
  誤差のあるリングは、線形のスピン軌道結合 (Derbenev–Kondratenko の (11/18)|γ∂n/∂γ|² の項) の脱偏極率を使う。
  FLUC つきのスピン追跡 (32 粒子) で測った脱偏極と比べた記録は docs/history/v2_changes_A.md。

条件: 3 GeV、12 セルの FODO リング (周長 45.6 m)、偏向磁石 L=1.2 m・15 度、RF 1.5 MV (HARM=64)。
"""
import io
import math

import numpy as np
import pytest

from pysad import ffscmd, emit6d, spin, track6, compat


# ------------------------------------------------------------------ B39
def _tradke_b(bzh):
    ph = emit6d._Phys(3.0e10, 100.0)          # 30 GeV: 拡散 dee を大きくして丸めの上に出す
    cod = np.array([1.0e-3, 4.0e-3, -2.0e-3, 1.0e-3, 0.0, 1.0e-3])
    st = emit6d._State(cod, ph)
    st.calpol = False
    st.codr0 = np.array([1.0e-3, 1.0e-3, -2.0e-3, -1.0e-3, 0.0, 1.0e-3])
    st.bzhr0 = bzh
    rng = np.random.default_rng(1)
    a = rng.normal(size=(4, 4)) * 1e-3
    st.B = np.zeros((6, 6))
    st.B[:4, :4] = a @ a.T                    # δ・z の行と列は 0: 放射のあとの B[δ,δ] がほぼ dee になる
    b0 = st.B.copy()
    st.tradke(0.1, 0.0, bzh)
    return b0, st.B.copy()


def test_kinetic_momentum_covariance_formula():
    """<(px + b y)(py − b x)> = B[py,px] − b (B[x,px] − B[y,py] + b B[x,y]) (K64 は B[y,py] の所が B[py,py])"""
    rng = np.random.default_rng(2)
    a = rng.normal(size=(4, 4))
    B = a @ a.T
    b = 0.37
    A = np.array([[0.0, 1.0, b, 0.0], [-b, 0.0, 0.0, 1.0]])
    cov = (A @ B @ A.T)[1, 0]
    assert abs(cov - (B[3, 1] - b * (B[0, 1] - B[2, 3] + b * B[0, 2]))) < 1e-13
    assert abs(cov - (B[3, 1] - b * (B[0, 1] - B[3, 3] + b * B[0, 2]))) > 1e-2
    assert abs((A @ B @ A.T)[0, 0] - (B[1, 1] + b * (2 * B[1, 2] + b * B[2, 2]))) < 1e-13      # px・px と py・py の項は K64 も正しい
    assert abs((A @ B @ A.T)[1, 1] - (B[3, 3] + b * (-2 * B[0, 3] + b * B[0, 0]))) < 1e-13


def test_tradke_px_py_diffusion_uses_b_y_py_v2():
    """SOL 内 (bzh≠0) の tradke: v2 と k64 の B[py,px] の差 = b·(B[y,py] − B[py,py])·dee (b = bzh/2、dee = 放射後の B[δ,δ])"""
    bzh = 0.3
    out = {}
    for mode in ('k64', 'v2'):
        old = compat._mode
        compat.set_mode(mode)
        try:
            out[mode] = _tradke_b(bzh)
        finally:
            compat.set_mode(old)
    b0, bk = out['k64']
    _b0, bv = out['v2']
    dee = bk[5, 5]
    assert dee > 0.0
    bk_, bv_ = bzh * 0.5, bzh      # k64: b = bzh/2、v2: B106 で b = bzh (力学的運動量 px + bzh y)
    want = dee * (-bv_ * (b0[0, 1] - b0[2, 3] + bv_ * b0[0, 2]) + bk_ * (b0[0, 1] - b0[3, 3] + bk_ * b0[0, 2]))
    assert abs((bv[3, 1] - bk[3, 1]) / want - 1.0) < 2e-2
    # B106: px・px、py・py の拡散も b = bzh (K64 は bzh/2)
    w11 = dee * (bv_ * (2 * b0[1, 2] + bv_ * b0[2, 2]) - bk_ * (2 * b0[1, 2] + bk_ * b0[2, 2]))
    w33 = dee * (bv_ * (-2 * b0[0, 3] + bv_ * b0[0, 0]) - bk_ * (-2 * b0[0, 3] + bk_ * b0[0, 0]))
    assert abs((bv[1, 1] - bk[1, 1]) / w11 - 1.0) < 2e-2 and abs((bv[3, 3] - bk[3, 3]) / w33 - 1.0) < 2e-2
    for i in range(6):
        for j in range(6):
            if (i, j) not in ((3, 1), (1, 3), (1, 1), (3, 3)):
                assert abs(bv[i, j] - bk[i, j]) <= 1e-10 * max(abs(bk[i, j]), 1e-12)      # v2: tradke の ddpx の符号 (A1 項目 7) で丸めの差だけ


# ------------------------------------------------------------------ B58・B64: 1 要素のスピン回転 (Emittance[] の写像と追跡)
ELS = {
    'bend0': 'BEND E1 = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);',
    'bendk': 'BEND E1 = (L=1.2 ANGLE=0.2617993877991494 K1=0.3 E1=0.1308996938995747 E2=0.1308996938995747);',
    'bendt': 'BEND E1 = (L=1.2 ANGLE=0.2617993877991494 K1=1e-9 E1=0.1308996938995747 E2=0.1308996938995747);',
    'quad': 'QUAD E1 = (L=0.3 K1=0.7);',
    'quadn': 'QUAD E1 = (L=0.3 K1=-0.7);',
    'mult2': 'MULT E1 = (L=0.3 K2=30 EPS=0.1);',                  # EPS=0.1: 追跡の分割も細かくする
    'mult12': 'MULT E1 = (L=0.3 K1=0.5 K2=30 SK2=10 EPS=0.1);',
}
OFF = [8e-3, 1e-3, 6e-3, -1e-3, 0.0, 0.0]


def _one(kind, cod):
    src = 'MOMENTUM = 3GEV;\nDRIFT D1 = (L=1);\n%s\nMARK ST = ();\nLINE RING = (ST E1 D1);\nFFS;\nUSE RING;\n' % ELS[kind]
    ffs, _I = ffscmd.run_sad_string(src, out=io.StringIO())
    els = [e for e in ffs.lat.elements if e.name == 'E1']
    p = ffs.momentum()
    ph = emit6d._Phys(p, 100.0)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne(els, np.array(cod, float), ph, rfsw=False, radcod=True)
    finally:
        emit6d.POL_MODE = old
    s_emit = st.srot[:, :3].copy()
    s_trk = np.zeros((3, 3))
    for j in range(3):
        s = [0.0, 0.0, 0.0]
        s[j] = 1.0
        _out, sp = track6.track_particles(els, np.array([cod], float), 1, p, 0.0, False, True, False, spins=[s],
                                          radpol=False, circ=100.0)
        s_trk[:, j] = sp[0]
    return s_emit, s_trk


@pytest.mark.parametrize('kind,tol', [('quad', 8e-6), ('quadn', 1e-5), ('mult2', 1e-6), ('mult12', 1e-5),
                                       ('bendk', 5e-6), ('bend0', 1.5e-4)])
def test_element_spin_rotation_emit_equals_tracking_v2(v2_mode, kind, tol):
    """軌道のずれた 1 要素: Emittance[] のスピン回転行列 = 追跡のスピン回転行列。
    (bend0 の 1.5e-4: K1=0 の追跡は放射の分割が粗く、軌道 0 でも 2e-5 違う。QUAD の 8e-6: 追跡の四極は分割が 1〜2 個で、
    v2 の Emittance[] (POL) は四極を細かく分割するので、分割の差が残る。端の項は MULT と BEND と、下の微分の試験で確かめる)"""
    a, b = _one(kind, OFF)
    assert np.abs(a - b).max() < tol, np.abs(a - b).max()


@pytest.mark.parametrize('kind,low', [('mult2', 5e-5), ('mult12', 1e-4), ('bendk', 1e-3), ('bend0', 5e-4)])
def test_element_spin_rotation_k64_differs(k64_mode, kind, low):
    a, b = _one(kind, OFF)
    assert np.abs(a - b).max() > low


def _spin_orbit_d(kind, cod):
    """1 要素のスピン軌道結合 D (3x6): Emittance[] の解析式 (srot(:,4:9)) と、追跡のスピン回転の中心差分"""
    src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1);
%s
MARK ST = ();
LINE RING = (ST E1 D1);
FFS;
USE RING;
''' % ELS[kind]
    ffs, _I = ffscmd.run_sad_string(src, out=io.StringIO())
    el =[e for e in ffs.lat.elements if e.name == 'E1'][0]
    p = ffs.momentum()
    ph = emit6d._Phys(p, 100.0)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne([el], np.array(cod, float), ph, rfsw=False, radcod=True)
    finally:
        emit6d.POL_MODE = old
    de = st.srot[:, 3:9].copy()

    def trk(c0):
        _o, sp = track6.track_particles([el], np.tile(c0, (3, 1)), 1, p, 0.0, False, True, False, spins=np.eye(3),
                                        radpol=False, circ=100.0)
        return np.asarray(sp).T
    c = np.array(cod, float)
    r0 = trk(c)
    h = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-4, 1e-5])
    dt = np.zeros((3, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h[j]
        a = (trk(c + d) - trk(c - d)) @ r0.T / (2 * h[j])          # = −[D e_j]×
        dt[:, j] = -np.array([a[2, 1] - a[1, 2], a[0, 2] - a[2, 0], a[1, 0] - a[0, 1]]) * 0.5
    return de, dt


def test_bend_spin_orbit_coupling_includes_fringe_v2(v2_mode):
    """BEND のスピン軌道結合 (脱偏極ベクトルの元) の y・py 成分: v2 の解析式が追跡の差分と 3% で一致する。
    端の縦磁場の積分 y/ρ の y についての微分が入る (K64 は bsir0 を定数として扱う)"""
    de, dt = _spin_orbit_d('bend0', OFF)
    assert np.abs(de - dt).max() < 0.03 * np.abs(dt).max()
    assert abs(de[2, 2] / dt[2, 2] - 1) < 0.02            # 縦 (z) まわりの歳差の y 依存
    de, dt = _spin_orbit_d('bendk', OFF)
    assert np.abs(de - dt).max() < 0.03 * np.abs(dt).max()


def test_bend_spin_orbit_coupling_k64_lacks_fringe(k64_mode):
    de, dt = _spin_orbit_d('bend0', OFF)
    assert abs(de[2, 2] / dt[2, 2]) < 0.5                  # K64 の Emittance[] は 1/6 しかない


def test_bend_spin_is_continuous_in_k1_v2(v2_mode):
    """K1→0 の BEND (tbendebody・tbendi) のスピン回転は K1=0 (tbendebody0・tbend) に連続につながる"""
    a0, t0 = _one('bend0', OFF)
    a1, t1 = _one('bendt', OFF)
    assert np.abs(a0 - a1).max() < 1e-5
    assert np.abs(t0 - t1).max() < 1.5e-4


def test_bend_spin_tracking_k64_is_discontinuous_in_k1(k64_mode):
    _a0, t0 = _one('bend0', OFF)
    _a1, t1 = _one('bendt', OFF)
    assert np.abs(t0 - t1).max() > 1e-3


# ------------------------------------------------------------------ B55
def test_erfcx_is_exp_z2_erfc():
    from scipy.special import erfc
    for z in (0.3 + 0.2j, 2 - 3j, 1.5 + 0.0j, -1 + 2j, 0.5j, 3 + 4j):
        ref = np.exp(z * z) * erfc(z)
        assert abs(spin._erfcx(z) / ref - 1.0) < 1e-11
    # 大きい引数: 1/(z √π) (1 − 1/(2 z²) + …)
    z = 6.29 - 17.344j
    assert abs(spin._erfcx(z) * z * math.sqrt(math.pi) / (1 - 0.5 / z**2 + 0.75 / z**4) - 1.0) < 1e-6


def _lat(dy):
    return '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY=%s) QDE = (L=0.3 K1=-0.6 DY=-%s);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
em = Emittance[];
Print["OUT ", {Polarization, PolarizationTime}/.em];
''' % (dy, dy)


def _pol(dy):
    buf = io.StringIO()
    ffs, _I = ffscmd.run_sad_string(_lat(dy), out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('OUT ')][0]
    v = line[4:].strip('{}').split(',')
    return ffs, float(v[0].replace('NaN', 'nan')), float(v[1])


def test_planar_ring_polarization_is_sokolov_ternov_v2(v2_mode):
    _f, pol, tau = _pol('0')
    assert abs(pol - 8.0 / (5.0 * math.sqrt(3.0))) < 1e-9
    assert abs(tau - 123.8) < 0.1


def test_polarization_k64_is_nan(k64_mode):
    _f, pol, _tau = _pol('0')
    assert math.isnan(pol)


def test_polarization_with_orbit_errors_v2(v2_mode):
    """四極の縦ずれ 3 mm (閉軌道 y 6 mm、n の傾き 0.04 rad): v2 は 88.4%。脱偏極率は 5.5e-11 /周 (ずれの 2 乗に比例)、
    平衡分布での n の広がり (1/2)<|N u|²> = 6.4e-6 (FLUC つきの追跡で測った 1−<S·n0> の 3e-6〜1.2e-5 と合う)"""
    ffs, pol, _tau = _pol('3e-3')
    # P1 (v2.2.0): 運動論的偏極 (−0.7%)・厳密な緩和率 (−0.3%)・刻みごとの脱偏極率 (−0.03%) で 0.8798 → 0.8695 (線形の式は 0.8840)。
    # 以前は 0.8798 ± 3e-3 だった許容を緩めて 1.5e-2 にしていたが、内訳が分かったので 0.8695 ± 3e-3 に戻す (下の試験が内訳を確かめる)
    assert abs(pol - 0.8695) < 3e-3
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    sp = e['spin']
    assert abs(sp['depol_rate_linear'] / 5.52e-11 - 1) < 0.05        # 線形の式 (SLIM)。v2next B3 で側帯の項を足した depol_rate は +13% (test_v2_B3_pol.py)
    assert abs(sp['spin_spread'] / 6.36e-6 - 1) < 0.05
    _f1, pol1, _t = _pol('1e-3')
    e1 = emit6d.emittance6d(_f1.lat, _f1.momentum(), rfsw=True, radcod=True, pol=True)
    assert abs(sp['depol_rate_linear'] / e1['spin']['depol_rate_linear'] / 9.0 - 1) < 0.1     # ずれの 2 乗に比例


def test_p1_terms_account_for_the_change_of_polarization_v2(v2_mode, monkeypatch):
    """P1 の 3 つの項 (spin.KINETIC・STEP_RATE・ST_EXACT) を全て切ると v2.1.0 (A2・B3) の値 0.8798 に戻る (内訳: 運動論的偏極 −7.1e-3、
    刻みごとの脱偏極率 −2.5e-4、厳密な緩和率 −3.0e-3。合計 −1.0e-2)。緩和率の τ は ST_EXACT で 124.18 → 123.73 s (平面リングの 123.8 s に近づく)"""
    _f, pol, tau = _pol('3e-3')
    monkeypatch.setattr(spin, 'KINETIC', False)
    monkeypatch.setattr(spin, 'STEP_RATE', False)
    monkeypatch.setattr(spin, 'ST_EXACT', False)
    _f, pol0, tau0 = _pol('3e-3')
    assert abs(pol0 - 0.8798) < 3e-3
    assert abs(tau0 - 124.18) < 0.05 and abs(tau - 123.73) < 0.05
    assert -1.3e-2 < pol - pol0 < -0.7e-2


def test_spin_orbit_response_matches_tracking_v2(v2_mode):
    """slim_depolarization の要素ごとの写像 s' = (1 − [D u]×) R s が、追跡の 1 周のスピンの軌道依存 (中心差分) を再現する"""
    ffs, _pol_, _tau = _pol('3e-3')
    p = ffs.momentum()
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
    els = ffs.lat.elements
    cod = np.array(e['orbit_entrance'])
    trf0 = e.get('trf0', 0.0)
    circ = emit6d._emit_circ(els)
    ph = emit6d._Phys(p, circ)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    secs = []
    c = cod.copy()
    try:
        for el in els:
            st = emit6d.tturne([el], c, ph, trf0=trf0, rfsw=True, radcod=True)
            secs.append((st.T + st.D, st.srot[:, :3].copy(), st.srot[:, 3:9].copy()))
            c = st.cod.copy()
    finally:
        emit6d.POL_MODE = old

    def model(u, s):
        for T, R, D in secs:
            s = (np.eye(3) - spin._cross_mat(D @ u)) @ (R @ s)
            u = T @ u
        return s

    def track(u, s):
        _o, sp = track6.track_particles(els, np.array([cod + u]), 1, p, trf0, True, True, False, spins=[s],
                                        radpol=False, circ=circ)
        return sp[0]
    s0 = np.array([0.3, 0.8, 0.52])
    s0 /= np.linalg.norm(s0)
    for u in ([1e-4, 0, 0, 0, 0, 0], [0, 0, 1e-4, 0, 0, 0], [0, 0, 0, 1e-4, 0, 0], [0, 0, 0, 0, 0, 1e-4]):
        u = np.array(u, float)
        dt = (track(u, s0) - track(-u, s0)) / 2
        dm = (model(u, s0) - model(-u, s0)) / 2
        assert np.abs(dt - dm).max() < 0.08 * np.abs(dt).max()
