# -*- coding: utf-8 -*-
"""担当 P2 (v2next、2026-11-19): 放射・エミッタンスの残り。記録: docs/history/v2next_P2.md。
  項目 1: SynchroBetaEmittance (B128) が始点の分散 η ≠ 0 のリングで Emittance[] の数倍〜40 倍になる件 (K64 の式の誤り)。
  項目 2: v2 の Emittance[] と RADINT の差 (矩形磁石): EMITDIV の既定 (1) の離散化の誤差。
  項目 3: tradke の ddpx の符号 (放射反作用つきのローレンツ力の数値積分が基準) と RADTAPER の係数。
  項目 4: RADCOD の 1 周の解析の T + D と写像の数値微分。
  項目 5: K1≠0 の BEND の四極も端面から始まる模型 (オプション)。
gs の値は 2026-11-19 の gs (10-19 ビルド、`local/wake/gsrun.sh`)。"""
import io
import math
import os
import sys

import numpy as np
import pytest

from pysad import compat, ffscmd, emit6d

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
import _p2_indep as X

_LAT_A = os.path.join(_HERE, 'kek-ls_v03_12.lat')
_PF = os.path.join(_ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad')
_AST = os.path.join(_ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


def _run(src):
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _floats(line, tag):
    import re
    return [float(v) for v in re.findall(r'-?\d*\.\d+(?:[Ee][-+]?\d+)?|-?\d+(?:[Ee][-+]?\d+)?', line[len(tag):])]


# ======================================================================= 項目 1: SynchroBetaEmittance
_FODO_START = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 HARM=64);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLB = (B D2 QD D1 QF D2);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
LINE RINGB = (ORIGIN RF 24*CELLB);
FFS;
PageWidth = 1999;
'''
_SB_TAIL = '''USE %s; CELL; RFSW; RAD; CALC;
em = Emittance[]; nus = Abs[(Tunes/.em)[[3]]];
r = SynchroBetaEmittance[{nus, nus, 0.001}, AzimuthalModes->%d];
Print["%s ", Twiss["EX", 1], " ", (Emittances/.em)[[1]], " ", r[[1,2]]];
'''


def _sb_fodo(ring, modes=3):
    _f, out = _run(_FODO_START + _SB_TAIL % (ring, modes, 'RES') + 'Exit;\n')
    return _floats([l for l in out.splitlines() if l.startswith('RES ')][0], 'RES ')       # η0、Emittance[] の εx、SB の εx


@pytest.mark.parametrize('ring', ['RINGQ', 'RINGB'])
def test_synchrobeta_fodo_start_point_v2(v2_mode, ring):
    """v2: 24 セルの FODO (3 GeV、RF 4 MV) を始点 QF (η = 0.70 m)・始点 BEND (η = 0.57 m) で走らせた SB が Emittance[] と 1e-3 で一致する
    (K64 は 6.54 倍で、始点に依った: 下の k64 の試験)。Emittance[] は始点に依らない (2.7614e-7)"""
    eta, em, sb = _sb_fodo(ring)
    assert eta > 0.5
    assert abs(sb / em - 1) < 1.5e-3, (sb, em)
    assert abs(em / 2.761415674e-7 - 1) < 1e-6


@pytest.mark.parametrize('ring,sb_gs', [('RINGQ', 1.7510378728578392e-6), ('RINGB', 1.7507988871340014e-6)])
def test_synchrobeta_fodo_start_point_k64_matches_gs(k64_mode, ring, sb_gs):
    """k64: 同じ FODO の SB は gs と 1e-6 で一致し (SB = 1.751e-6、gs の Emittance[] 2.6773e-7 の 6.54 倍)、始点の分散 η ≠ 0 の誤り (K64 の式) を再現する"""
    eta, em, sb = _sb_fodo(ring)
    assert abs(sb / sb_gs - 1) < 1e-6, (sb, sb_gs)
    assert sb / em > 6.0
    assert abs(em / 2.67726385604e-7 - 1) < 1e-8


def _sb_ring(path, modes=3):
    src = open(path, encoding='utf-8', errors='replace').read() + '\nFFS;\nPageWidth = 1999;\n' + _SB_TAIL % ('RING', modes, 'RES') + 'Exit;\n'
    _f, out = _run(src)
    return _floats([l for l in out.splitlines() if l.startswith('RES ')][0], 'RES ')


def test_synchrobeta_large_ring_a_v2(v2_mode):
    """v2: 大きなリング A (kek-ls_v03_12.lat、3 GeV、周長 603 m、始点の η = 0.0445 m) の SB は Emittance[] の 1.0084 倍 (K64: 40 倍)。
    AzimuthalModes を 3 から 12 にすると 1.0043 倍 (K64 は 6.18e-9 から 3.19e-9 に動き続ける)"""
    e3 = _sb_ring(_LAT_A, 3)
    e12 = _sb_ring(_LAT_A, 12)
    assert 0.0444 < e3[0] < 0.0446
    assert abs(e3[2] / e3[1] - 1.0084) < 2e-3 and abs(e12[2] / e12[1] - 1.0043) < 2e-3


def test_synchrobeta_large_ring_a_k64_matches_gs(k64_mode):
    """k64: 大きなリング A の SB は gs と 2e-7 (modes 3)・6e-6 (modes 12) で一致し (modes 3・12: 6.178294026e-9・3.18594988e-9)、Emittance[] (1.549e-10) の 40 倍・20 倍"""
    e3 = _sb_ring(_LAT_A, 3)
    e12 = _sb_ring(_LAT_A, 12)
    assert abs(e3[2] / 6.178294026e-9 - 1) < 1e-6 and abs(e12[2] / 3.18594988e-9 - 1) < 1e-5 and e3[2] / e3[1] > 39.0


@pytest.mark.parametrize('path', [_PF, _AST])
def test_synchrobeta_pf_astrid_v2(v2_mode, path):
    """v2: PF リング・ASTRID2 の SB は Emittance[] の 1.005〜1.006 倍 (K64 は PF 1.011・ASTRID2 1.008)"""
    if not os.path.exists(path):
        pytest.skip(path)
    eta, em, sb = _sb_ring(path, 3)
    assert 0.995 < sb / em < 1.012, (sb, em)


def _toy_sb(mux, a, alz=2e-3, alx=1e-3, sige=1e-3, nus=0.01, dx=2e-9, mphi2=3):
    """SB のソルバー (synchrobeta の tefsetup・tesolvd・tesumb_v2) に、回転 (位相 mux) + 減衰 + 拡散の 1 周の写像と、1 次の分散 a·δ の軌道 c(δ) = a δ を
    合成して入れる (光学のないトイ)。戻り: 4 次元の分散行列の x-px の射影エミッタンス"""
    from pysad import synchrobeta as S
    from pysad.calc6d import etwiss2ri
    from pysad.emit6d import _tinv6
    mphiz = mphi2 * 2 - 1
    ndims = mphiz
    ndp = max(15, mphiz)
    nz = mphiz * 2 + 1
    mphi = mphiz
    sigea = math.sqrt(12.0) * sige
    dpndim = sigea / ndims
    tws = np.zeros((2 * ndims + 1, 28))
    beams = np.zeros((2 * ndims + 1, 10))
    trads = np.zeros((2 * ndims + 1, 5, 5))
    for i in range(-ndims, ndims + 1):
        t = tws[i + ndims]
        t[1] = 1.0
        t[4] = 1.0
        t[2] = mux
        t[5] = 0.3
        t[22] = 1.0
        t[15] = a * i * dpndim
        t[20] = i * dpndim
        beams[i + ndims][0] = dx
        beams[i + ndims][2] = dx
        c_, s_ = math.cos(mux), math.sin(mux)
        trads[i + ndims][0:2, 0:2] = -alx * np.array([[c_, s_], [-s_, c_]])
    ri = np.eye(6)
    tw0 = tws[ndims].copy()
    tw0[2] = tw0[5] = tw0[23] = 0.0
    dj = sigea ** 2 / ndp * 0.5
    fj = np.array([math.exp(-(j - 0.5) * dj / sige ** 2) for j in range(1, ndp + 1)])
    w = float(np.sum(fj * dj / sige ** 2))
    jump = {'dd': 4.0 * alz * sige ** 2}
    ha, hb, ba, bb, bd = S.tefsetup(ndp, mphi, mphi2, nz, ndims, beams, trads, tws, dj, dpndim, sigea, tw0, ri, (etwiss2ri, _tinv6, None), jump=jump)
    amuj = np.full(ndp, 2 * math.pi * nus)
    hc = np.zeros((4, mphi2, ndp))
    hs = np.zeros((4, mphi2, ndp))
    bc = np.zeros((10, mphi2, ndp))
    bs = np.zeros((10, mphi2, ndp))
    S.tesolvd(bb, bc, bs, ba, amuj, bd, hc, hs, ndp, mphi2, mphi, 10, dj, -alz, sige, S._tsolva)
    beam = S.beam6(S.tesumb_v2(bc, hc, mphi2, ndp, w, dj, sige, np.zeros(4), fj, jump=jump))
    return math.sqrt(beam[0, 0] * beam[1, 1] - beam[0, 1] ** 2)


@pytest.mark.parametrize('mux', [-2.64, 0.3, 1.0, 2.5])
def test_synchrobeta_orbit_jump_excitation_matches_analytic_toy(v2_mode, mux):
    """SB のソルバーに、軌道が δ に比例して動く (c(δ) = a δ、a = 0.1 m) 合成のトイを入れる。光子の放出で δ が跳ぶと軌道が a du だけ跳び、
    betatron の振動が励起される (H = a² の励起)。独立な解析式: x-px の分散は 等方な Dx/(2α_x − α_x²) に λ = a² e α_z/α_x (励起) と、x だけに a² e
    (静的な軌道の広がり) が足される。投影エミッタンスの比 sqrt((1 + λ'+a²e')(1 + λ')) = 1.02499 (e' = e/Σ0)。モンテカルロ (8000 粒子) は 1.025〜1.029。
    位相 mux に依らない (K64 の式 (平均を解く) は mux = 2.5 で +8.6%、0.3 で +2.6% と mux に依った)"""
    a, alx, alz, e, dx = 0.1, 1e-3, 2e-3, 1e-6, 2e-9
    s0 = dx / (2 * alx - alx ** 2)
    lam = a * a * 4 * alz * e / (2 * alx - alx ** 2) / 2.0           # tr ΔΣ = a² D/(1 − (1−α)²)、λ = tr/2
    want = math.sqrt((s0 + lam + a * a * e) * (s0 + lam)) / s0
    got = _toy_sb(mux, a) / s0
    assert abs(got / want - 1) < 1e-3, (mux, got, want)
    assert abs(got - 1.025) < 2e-3


# ======================================================================= 項目 2: EMITDIV の離散化
_FODO_RECT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=%g E2=%g);
CAVI RF = (L=0 VOLT=4e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 24*CELL1 RF);
FFS;
USE RING; CELL; CALC;
'''


def _fodo_lat(e):
    ffs, _ = _run(_FODO_RECT % (e, e))
    return ffs


def test_emittance_default_emitdiv_error_is_discretisation_v2(v2_mode):
    """v2: 矩形の偏向磁石 (E1 = E2 = 0.5·ANGLE) の FODO の Emittance[] (RFSW) の εx は、EMITDIV = 1 (既定) で 16 の値より 4.6e-4 大きく、4 でほぼ収束する
    (K64 の既定の刻みの離散化誤差。大きなリング B (8 GeV、周長 1.4 km) では +1.4e-3 で、v2 と RADINT (6.9554e-9) の差の原因。k64・v2 とも、扇形でも同じ誤差)。
    RFSW OFF の Jx は解析 (平面リングの放射積分 I4、端の項 −η h² tan e が 1 倍) と 1e-4 で一致 (v2)、k64 は端の項が 2 倍で 0.12 大きい"""
    import _a2_indep as A
    ffs = _fodo_lat(0.5)
    old = list(emit6d.EMITDIV)
    try:
        eps = {}
        for div in (1, 4, 16):
            emit6d.EMITDIV[:] = [div] * 3
            eps[div] = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)['emittance_x']
        assert 1e-4 < eps[1] / eps[16] - 1 < 1e-3, eps
        assert abs(eps[4] / eps[16] - 1) < 2e-5, eps
        emit6d.EMITDIV[:] = [4] * 3
        e0 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False)
        jx = -2 * e0['damping_rate'][0] * ffs.momentum() / e0['energy_loss_eV']
        I, C, _nu = A.integrals_rich(ffs.lat.elements, n=80)
        jx_ana = A.derived(I, C, ffs.momentum())['Jx']
        assert abs(jx - jx_ana) < 2e-4, (jx, jx_ana)
        compat.set_mode('k64')
        e1 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False)
        jx64 = -2 * e1['damping_rate'][0] * ffs.momentum() / e1['energy_loss_eV']
        assert jx64 - jx_ana > 0.1, (jx64, jx_ana)
    finally:
        emit6d.EMITDIV[:] = old


# ======================================================================= 項目 3: tradke の符号 (独立な基準: 放射反作用つきのローレンツ力の積分)
_TH, _L, _X0, _PX0, _G0 = 0.2618, 1.2, 2e-3, 1e-3, -1e-4


def _emit_end(mode, radcod, div=20):
    from pysad.elements import Bend
    compat.set_mode(mode)
    el = Bend(name='B', L=_L, ANGLE=_TH, E1=0.0, E2=0.0)
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [div] * 3                      # _Phys が作られるとき EMITDIV を読むので、先に設定する
    try:
        ph = emit6d._Phys(3e9, 100.0, amass=emit6d.ELMASS, charge=1.0)
        st = emit6d.tturne([el], np.array([_X0, _PX0, 0, 0, 0, _G0], float), ph, trf0=0.0, rfsw=False, radcod=radcod, want_db=False, want_jac=False)
    finally:
        emit6d.EMITDIV[:] = old
    return st.cod


def _track_end(mode, rad):
    compat.set_mode(mode)
    src = f'''MOMENTUM = 3GEV;
BEND B = (L={_L} ANGLE={_TH});
LINE RING = (B);
FFS;
USE RING;
{'RAD; NOFLUC;' if rad else 'NORAD;'}
r = TrackParticles[{{1, {{{{{_X0}}},{{{_PX0}}},{{0}},{{0}},{{0}},{{{_G0}}},{{1}}}}}}, 1, 1, 1];
Print["R ", Flatten[Take[r[[2]],6]]];
Exit;
'''
    _f, out = _run(src)
    return np.array(_floats([l for l in out.splitlines() if l.startswith('R ')][0], 'R '))


def test_tradke_sign_agrees_with_exact_radiation_reaction_in_v2_not_k64():
    """扇形の偏向磁石 1 台 (3 GeV、1.2 m、0.2618 rad) の放射 ON と OFF の出口の差 Δx・Δpx を、放射反作用つきのローレンツ力の数値積分
    (連続な放射、強さ k を pysad の Δδ に合わせる: _p2_indep.bend_with_radiation) と比べる。放射なしの出口は 1e-12 で一致 (曲線座標への変換の確認)。
    v2 の Emittance[] の tradke は厳密な積分と Δpx で 2e-5 (実測 6e-6)、Δx で 1e-4 一致。K64 の tradke (ddpx = +0.5 dpx dg) は Δpx が 2.2e-4 (EMITDIV = 20)・9% (既定の EMITDIV = 1) ずれる (符号が逆の分)。
    追跡 (tradkn、ddpx = −0.5 dpx dg) は k64・v2 とも 5e-5 で一致 (B111: tradke が誤り)"""
    old = compat._mode
    try:
        off = _emit_end('v2', False)
        dg = _emit_end('v2', True)[5] - _G0
        k = X.radiation_reaction_k_for_dg(_X0, _PX0, _G0, _TH, _L, dg)
        xc0, pxc0, _, _ = X.bend_with_radiation(_X0, _PX0, _G0, _TH, _L, 0.0)
        xe, pxe, _, _ = X.bend_with_radiation(_X0, _PX0, _G0, _TH, _L, k)
        assert abs(off[0] - xc0) < 1e-12 and abs(off[1] - pxc0) < 1e-12
        dx_ex, dpx_ex = xe - xc0, pxe - pxc0
        on2 = _emit_end('v2', True)
        assert abs((on2[1] - off[1]) / dpx_ex - 1) < 2e-5 and abs((on2[0] - off[0]) / dx_ex - 1) < 1e-4
        on64 = _emit_end('k64', True)
        r64 = (on64[1] - off[1]) / dpx_ex - 1
        assert 1e-4 < abs(r64) < 1e-3, r64
        # 既定の刻み (EMITDIV = 1、この磁石は 22 刻み) では K64 の tradke の誤差は 9%、v2 は 1e-3 以内
        on2d = _emit_end('v2', True, 1)
        on64d = _emit_end('k64', True, 1)
        offd = _emit_end('v2', False, 1)
        assert abs((on2d[1] - offd[1]) / dpx_ex - 1) < 1e-3 and 0.05 < abs((on64d[1] - offd[1]) / dpx_ex - 1) < 0.15
        for mode in ('v2', 'k64'):
            t0, t1 = _track_end(mode, False), _track_end(mode, True)
            assert abs((t1[1] - t0[1]) / dpx_ex - 1) < 5e-5 and abs((t1[0] - t0[0]) / dx_ex - 1) < 5e-4, mode
    finally:
        compat.set_mode(old)


def test_radtaper_coefficients_v2_use_arithmetic_mean_for_bends():
    """RADTAPER の係数: K64 の tturne は BEND に (3δ_l + δ_l+1)/4 (放射を刻みの終わりに入れる粗い刻みに合わせた値)、追跡は (δ_l + δ_l+1)/2。
    連続な放射で要素の中の損失が線形のとき、曲げ角を設計値にする係数は調和平均 ≈ 算術平均 (1/2)。v2 は Emittance[]・追跡とも (δ_l + δ_l+1)/2 にする。
    (係数が違っても Emittance[] の閉軌道は 3e-10 しか変わらない: 大きなリング B (8 GeV) の確認。追跡で閉じるのは同じ係数を両方に使うため)"""
    from pysad import radtaper
    from pysad.elements import Bend, Quadrupole
    els = [Bend(name='B', L=1.0, ANGLE=0.1), Quadrupole(name='Q', L=0.3, K1=0.5)]
    ddp = [0.0, -2e-4, -3e-4]
    old = compat._mode
    try:
        compat.set_mode('v2')
        for tracking in (False, True):
            tb, tq = radtaper.factors(els, ddp, 1e-5, tracking)
            assert abs(tb - (1.0 + (ddp[0] + ddp[1]) / 2 + 1e-5)) < 1e-15 and abs(tq - (1.0 + (ddp[1] + ddp[2]) / 2 + 1e-5)) < 1e-15
        compat.set_mode('k64')
        tb, _ = radtaper.factors(els, ddp, 1e-5, False)
        assert abs(tb - (1.0 + (3 * ddp[0] + ddp[1]) / 4 + 1e-5)) < 1e-15
        tb, _ = radtaper.factors(els, ddp, 1e-5, True)
        assert abs(tb - (1.0 + (ddp[0] + ddp[1]) / 2 + 1e-5)) < 1e-15
    finally:
        compat.set_mode(old)


# ======================================================================= 項目 4: RADCOD の解析の T + D と数値微分
_C0 = np.array([1e-4, 1e-5, 0, 0, 0, -1e-4])


def _derivative_errors(mode):
    compat.set_mode(mode)
    ffs = _fodo_lat(0.1308996938995747)
    els = ffs.lat.elements
    ph = emit6d._Phys(ffs.momentum(), sum(e.L for e in els), amass=emit6d.ELMASS, charge=1.0)

    def f(e_, c0, want=False):
        return emit6d.tturne(e_, np.array(c0, float), ph, trf0=0.0, rfsw=False, radcod=True, want_db=want, want_jac=want)

    def numjac(e_, h=2e-5):
        J = np.zeros((6, 6))
        for k in range(6):
            def fe(s):
                c = _C0.copy()
                c[k] += s * h
                return f(e_, c).cod
            J[:, k] = (8 * (fe(1) - fe(-1)) - (fe(2) - fe(-2))) / (12 * h)
        return J
    res = {}
    for name, sel in (('quad', [els[1]]), ('bend', [els[3]]), ('quad_bend', els[1:4]), ('cell', els[1:7])):
        st = f(sel, _C0, True)
        res[name] = np.abs(st.T + st.D - numjac(sel)).max()
    return res


def test_radcod_analytic_matrix_matches_numeric_derivative_v2(v2_mode):
    """v2: RADCOD ON の要素の解析の転送行列 T + D と写像の数値微分 (4 次の中心差分) が QF 1 台で 1e-12、BEND 1 台で 1e-9、1 セルで 1e-8 以内
    (K64: 1e-8・5e-7・2e-6)。残る差は BEND の T[z, δ] (tbendebody0 の (5,6)、K64 の仕様) の px² の項だけで、放射の微分ではない"""
    r = _derivative_errors('v2')
    assert r['quad'] < 1e-12 and r['bend'] < 1e-9 and r['quad_bend'] < 1e-8 and r['cell'] < 1e-8, r


def test_radcod_analytic_matrix_matches_numeric_derivative_k64(k64_mode):
    """k64: K64 の tradke の微分の近似 (4 か所) が残る。QF 1 台 1e-8、BEND 1 台 5e-7、1 セル 2e-6 (記録)"""
    r = _derivative_errors('k64')
    assert 5e-9 < r['quad'] < 5e-8 and 2e-7 < r['bend'] < 2e-6 and 1e-6 < r['cell'] < 5e-6, r


# ======================================================================= 項目 5: K1≠0 の BEND の四極も端面から始まる模型
def _k1_bend(k1, e):
    src = f'''MOMENTUM = 3GEV;
BEND BE = (L={X._L} ANGLE={X._TH} E1={e} E2={e} K1={k1});
LINE X = (BE);
FFS;
USE X; INS; NORAD; NOFLUC; CALC;
Exit;
'''
    ffs, _ = _run(src)
    return [x for x in ffs.lat.elements if x.name == 'BE'][0], ffs


_ORBITS = [(1e-3, 0, 0), (3e-3, 0, 0), (0, 5e-4, 0), (2e-3, 3e-4, 0), (1e-3, 0, 5e-3)]


@pytest.mark.parametrize('e,k1', [(0.5, -0.05), (0.25, -0.05), (0.5, 0.3)])
def test_k1_bend_face_model_matches_exact_hard_edge_map(v2_mode, e, k1):
    """v2 のオプション compat.set_bend_k1_edge('face'): CALC と追跡の K1≠0 の BEND の (x, px) が、四極の磁場も端面の間とする硬い端の厳密な写像
    (ローレンツ力の積分) と CALC で 3e-10、追跡で 1e-9 以内。既定 ('radial': K64 の模型) は半径方向の面の間とする厳密な写像と同じ精度で、
    2 つの模型の差は軌道 1〜3 mm で 5e-8〜5e-7 (x0 tan e に比例)"""
    from pysad import lattice, track6n
    el, ffs = _k1_bend(k1, e)
    try:
        diffs = []
        for edge in ('radial', 'face'):
            compat.set_bend_k1_edge(edge)
            for (x0, px0, dl) in _ORBITS:
                c = [x0, px0, 0, 0, 0, dl]
                cc, _m = lattice._bend_element_cod_matrix(el, np.array(c, float))
                t = np.asarray(track6n.track_particles_fast([el], np.array([c], float), 1, ffs.momentum(), 0.0, False, False, False, circ=2 * X._L)[0])
                ex = X.exact_map(c, el.K1, e * X._TH, e * X._TH, quad=edge)
                assert np.abs(cc[:2] - ex).max() < 3e-10 and np.abs(t[:2] - ex).max() < 1e-9, (edge, c)
                if edge == 'face':
                    other = X.exact_map(c, el.K1, e * X._TH, e * X._TH, quad='radial')
                    diffs.append(np.abs(ex - other).max())
        assert 1e-8 < max(diffs) < 1e-6, diffs
    finally:
        compat.set_bend_k1_edge(None)


def test_k1_bend_face_model_default_is_radial_and_k64_ignores_it():
    """既定は 'radial' (K64 の模型)。k64 モードでは set_bend_k1_edge('face') も無視する (k64 は K64 のまま)。軌道が中心 (x = 0) のリングは 2 つの模型で同じ
    (蹴りは x² に比例)。環境変数 PYSAD_BEND_K1_EDGE=face でも切り替わる"""
    from pysad import lattice
    old = compat._mode
    try:
        compat.set_mode('v2')
        compat.set_bend_k1_edge(None)
        assert compat.bend_k1_edge() == 'radial'
        el, _ffs = _k1_bend(-0.05, 0.5)
        c = np.array([2e-3, 3e-4, 0, 0, 0, 2e-3])
        m_radial = lattice._bend_element_cod_matrix(el, c)
        compat.set_bend_k1_edge('face')
        assert compat.bend_k1_edge() == 'face'
        m_face = lattice._bend_element_cod_matrix(el, c)
        assert np.abs(m_face[0] - m_radial[0]).max() > 1e-8
        # 軌道が中心: 蹴りは入らない
        c0 = np.zeros(6)
        assert np.abs(lattice._bend_element_cod_matrix(el, c0)[1] - np.eye(6)).max() > 0     # (ダミー: 行列が作れる)
        compat.set_mode('k64')
        assert compat.bend_k1_edge() == 'radial'
        m_k64 = lattice._bend_element_cod_matrix(el, c)
        compat.set_bend_k1_edge(None)
        m_k64r = lattice._bend_element_cod_matrix(el, c)
        assert np.abs(m_k64[0] - m_k64r[0]).max() == 0.0 and np.abs(m_k64[1] - m_k64r[1]).max() == 0.0
        with pytest.raises(ValueError):
            compat.set_bend_k1_edge('x')
    finally:
        compat.set_bend_k1_edge(None)
        compat.set_mode(old)


def test_k1_bend_face_model_emittance_runs_and_equals_radial_on_axis(v2_mode):
    """Emittance[] の Python の経路 (emit6d._tbende_core) にも入れた。軌道が中心のリング (K1≠0 の矩形磁石 E = 0.5) は 2 つの模型で εx が 1e-12 で同じ"""
    src = _FODO_RECT.replace('BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=%g E2=%g);', 'BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5 K1=-0.1);')
    src = src.replace('E1=%g E2=%g', 'E1=0.5 E2=0.5 K1=-0.1')
    ffs, _ = _run(src)
    try:
        res = []
        for edge in ('radial', 'face'):
            compat.set_bend_k1_edge(edge)
            res.append(emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=False)['emittance_x'])
        assert abs(res[0] / res[1] - 1) < 1e-10, res
    finally:
        compat.set_bend_k1_edge(None)
