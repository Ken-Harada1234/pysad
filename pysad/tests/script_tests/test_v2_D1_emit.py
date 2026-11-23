# -*- coding: utf-8 -*-
"""担当 D1 (v2next): Emittance[] の NORFSW + 閉軌道のずれ、WSPAC の最後のスライス、光子 (MULT・SEXT・OCT・BEAMBEAM + BEND・TRPT)、
SynchroBetaEmittance の警告、K64 の spdepol の性質、不安定なリングの POL。記録: docs/history/v2next_D1.md。
gs の値は 2026-11-18 の gs (10-19 ビルド、`local/wake/gsrun.sh`) の出力。"""
import io
import math
import os
import re
import sys

import numpy as np
import pytest

from pysad import compat, ffscmd, emit6d

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_HERE)))
_LAT_A = os.path.join(_HERE, '..', 'kek-ls_v03_12.lat')
_PF = os.path.join(_ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad')
_AST = os.path.join(_ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


def _run(src):
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _tag(out, tag):
    for l in out.splitlines():
        if l.startswith(tag):
            return [float(x) for x in re.findall(r'-?\d*\.\d+(?:[Ee][-+]?\d+)?|-?\d+(?:[Ee][-+]?\d+)?', l[len(tag):])]
    raise AssertionError(tag + ' がない: ' + out[-300:])


# ================================================================== 項目 1: NORFSW; RAD の Emittance[] と閉軌道のずれ
_FODO = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=%s) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=%s DX=%g) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="14.10";
USE RING; CELL; CALC;
%s
e = Emittance[];
Print["NR ", Emittances/.e, Tunes/.e, DampingRate/.e];
Exit;
'''


@pytest.mark.parametrize('flags,want', [
    ('NORFSW; RAD; NOFLUC;', (2.70000647e-7, 1.92080517e-5, -0.2932001422, -0.1330809406, -1.2328621e-4, -1.3027092e-4)),
    ('NORFSW; RAD; RADCOD; NOFLUC;', (2.68403700e-7, 1.92471880e-5, -0.2928504850, -0.1327463304, None, None)),
])
def test_norfsw_emittance_with_closed_orbit_offset_matches_gs_k64(k64_mode, flags, want):
    """k64: 12 セル FODO (3 GeV、QF の DX = 1e-4) の NORFSW + RAD の Emittance[] が gs と一致する。以前は εx が 2% 違い、チューンが 1e-3 違った
    (K64 の tcod は RFSW OFF で x・px・y・py だけ閉じ δ・z を動かさないが、pysad は z の戻りを δ で閉じていた: D1)"""
    _ffs, out = _run(_FODO % ('1', '0.7', 1e-4, flags))
    v = _tag(out, 'NR ')
    assert abs(v[0] / want[0] - 1) < 2e-8 and abs(v[2] / want[1] - 1) < 1e-8, v
    assert abs(v[3] - want[2]) < 2e-9 and abs(v[4] - want[3]) < 2e-9, v
    if want[4] is not None:
        assert abs(v[6] / want[4] - 1) < 1e-7 and abs(v[7] / want[5] - 1) < 1e-7, v


def test_norfsw_emittance_with_orbit_offset_large_rings_match_gs_k64(k64_mode):
    """k64: 大きなリング A (kek-ls_v03_12.lat、3 GeV、周長 603 m) と PF リング (2.5 GeV) の NORFSW + RAD、QF に DX = 1e-4・QD に DY = −8e-5:
    εx・εy が gs と一致する (A: εx 1.5803933e-10・εy 3.1222687e-11、PF: 3.62598585e-8・6.7367214e-12)"""
    for path, q1, q2, ex, ey in [(_LAT_A, 'QFA1H.1', 'QDA2.1', 1.5803933e-10, 3.1222687e-11), (_PF, 'Q1.1', 'Q2.1', 3.62598585e-8, 6.7367214e-12)]:
        if not os.path.exists(path):
            pytest.skip(path)
        src = open(path, encoding='utf-8', errors='replace').read() + f'''
FFS;
USE RING; CELL; CALC;
LINE["DX","{q1}"]=1e-4;
LINE["DY","{q2}"]=-8e-5;
NORFSW; RAD; NOFLUC;
e = Emittance[];
Print["NR ", Emittances/.e, Tunes/.e];
Exit;
'''
        _ffs, out = _run(src)
        v = _tag(out, 'NR ')
        assert abs(v[0] / ex - 1) < 1e-7 and abs(v[1] / ey - 1) < 1e-6, (path, v)


def _jx(dxf):
    ffs, _out = _run(_FODO % ('1.6', '0.7', dxf, 'NORFSW; RAD; NOFLUC;'))
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False)
    return -2.0 * e['damping_rate'][0] * ffs.momentum() / e['energy_loss_eV']


def test_norfsw_jx_change_from_orbit_offset_matches_planar_integrals_v2(v2_mode):
    """v2: NORFSW (δ = 0 の軌道、RF の周波数と無関係) の Emittance[] の Jx が、軌道が四極の中心からずれたとき (QF の DX = 1e-4) の変化について
    平面リングの放射積分 (_c3_indep.planar_quad_offset_jx: 3x3 行列、h = k x_o、I4 の 2 k² x_o η) の ΔJx = −0.0052 と 10% で一致する (実測 −0.0054)。
    リング: QF D2 B D2 QD D1 (1.6 m) × 12 (planar_quad_offset_jx と同じセル)"""
    sys.path.insert(0, _HERE)
    import _c3_indep as C
    d_emit = _jx(1e-4) - _jx(0.0)
    d_want = C.planar_quad_offset_jx(1e-4, 0.0) - C.planar_quad_offset_jx(0.0, 0.0)
    assert abs(d_emit / d_want - 1) < 0.10, (d_emit, d_want)


# ================================================================== 項目 2: WSPAC の最後のスライスの p1 (gammab(nlat+1) = 0)
_WS = '''MOMENTUM = 0.2GEV;
PBUNCH = 1E8;
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
Print["ZERO ", Tunes/.e0];
WSPAC;
e = Emittance[];
Print["SC ", Tunes/.e];
Exit;
'''


def test_wspac_last_slice_uses_gammab_zero_in_k64_and_physical_p1_in_v2(monkeypatch):
    """K64 の tintrb の最後のスライス (l = nlat) は twspace で gammab(nlat+1) = 0 を読み p1 が半分になり、そのレンズが 8 倍強い (gs で確認:
    GAMB 391.39 0.000 195.70)。k64 はこれを再現し (チューンシフト −4.35546e-3・−8.84782e-3 が gs と全桁一致)、v2 は物理的な p1 を使う
    (シフト −4.13862e-3・−8.46689e-3: ラスレットの式と 0.15%。v2 は A・C の修正でゼロ電流の光学も少し違う)"""
    # 参照値は EMITDIV = 1 で作った (v2.4.0 から v2 の既定は 4 で、この低エミッタンスの FODO ではシフトが 5% 動く: 1 → −4.129e-3、4 → −4.305e-3、32 → −4.339e-3)
    monkeypatch.setenv('PYSAD_EMITDIV_DEFAULT', '1')
    old = compat._mode
    try:
        res = {}
        for mode in ('k64', 'v2'):
            compat.set_mode(mode)
            _f, out = _run(_WS)
            z, s = _tag(out, 'ZERO ')[:2], _tag(out, 'SC ')[:2]
            res[mode] = [s[0] - z[0], s[1] - z[1]]
    finally:
        compat.set_mode(old)
    assert abs(res['k64'][0] + 4.3554576e-3) < 2e-9 and abs(res['k64'][1] + 8.8478243e-3) < 2e-9, res
    # v2next P2: tradke の微分を厳密にした (放射による D 行列が変わり、放射による tune shift の 1e-5 の量が変わる) ので v2 の値が 2e-3 (相対) 動いた (旧: −4.1386187e-3・−8.4668911e-3)
    assert abs(res['v2'][0] + 4.1294834e-3) < 2e-9 and abs(res['v2'][1] + 8.4361235e-3) < 2e-9, res
    assert 1.05 < res['k64'][0] / res['v2'][0] < 1.06 and 1.04 < res['k64'][1] / res['v2'][1] < 1.06


@pytest.mark.parametrize('path,ring,pb,want', [
    (_AST, 'RING', '1e9', (0.1799774472, 0.1401141997, 1.26803995e-8)),
    (_PF, 'RING', '1e10', (-0.4024743242, 0.2789345926, 3.45421908e-8)),
    (_LAT_A, 'RING', '1e10', (-0.4225238563, -0.401707958, 1.5497111e-10)),
])
def test_wspac_emittance_matches_gs_k64_other_rings(k64_mode, path, ring, pb, want):
    """k64: WSPAC の Emittance[] (MINCOUP = 0.1) が ASTRID2 (0.58 GeV)、PF リング (2.5 GeV)、大きなリング A (3 GeV、603 m) で gs と全桁一致"""
    if not os.path.exists(path):
        pytest.skip(path)
    src = open(path, encoding='utf-8', errors='replace').read() + f'''
PBUNCH={pb}; MINCOUP=0.1;
FFS;
$FORM="14.10";
USE {ring}; CELL; CALC;
RFSW; RAD; NOFLUC;
WSPAC;
e = Emittance[];
Print["SC ", Tunes/.e, Emittances/.e];
Exit;
'''
    _f, out = _run(src)
    v = _tag(out, 'SC ')
    assert abs(v[0] - want[0]) < 2e-9 and abs(v[1] - want[1]) < 2e-9 and abs(v[3] / want[2] - 1) < 1e-8, (v, want)


def test_wspac_v2_tune_shift_large_ring_a_matches_laslett_integral_v2(v2_mode):
    """v2: 大きなリング A (3 GeV、603 m、3563 要素) の WSPAC のチューンシフトが一次のラスレットの積分 (各要素を入口・中点・出口の Simpson。
    σx² = β εx + η² σδ²) と 0.5% (x)・1% (y) で一致する (実測 0.9982・1.0041)。k64 は最後のスライスのレンズが 8 倍強いので +0.9%・+1.1%"""
    from pysad import touschek
    if not os.path.exists(_LAT_A):
        pytest.skip(_LAT_A)
    src = open(_LAT_A, encoding='utf-8', errors='replace').read() + '\nFFS;\nUSE RING; CELL; CALC;\nRFSW; RAD; NOFLUC;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    els = ffs.lat.elements
    p = ffs.momentum()
    em = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=False)
    pb = 1e8
    t = touschek.touschek_table(ffs.lat, p, em, pb, rfsw=True, intra=False, wspac=True, verbose=False, mincoup=0.1)
    shift = (np.array(t['tunes']) - np.array(em['tunes']))[:2]
    ex, ey, sige, sigz = t['EMITXC'], t['EMITYC'], t['momentum_spread'], t['bunch_length']
    ph = emit6d._Phys(p, emit6d._emit_circ(els))
    fac = 2 * pb * emit6d.ELRADI / (ph.p0**2 * ph.h0 * math.sqrt(2 * math.pi) * sigz)
    rows = [ffs.row(k) for k in range(1, len(els) + 2)]

    def f(r):
        sx, sy = math.sqrt(r.betx * ex + r.etx**2 * sige**2), math.sqrt(r.bety * ey)
        return r.betx / (sx * (sx + sy)), r.bety / (sy * (sx + sy))
    dx = dy = 0.0
    for i, el in enumerate(els):
        if el.L == 0:
            continue
        fa, fb, fc = f(rows[i]), f(ffs.row_frac(i + 1, 0.5)), f(rows[i + 1])
        dx += -(1 / (4 * math.pi)) * fac * el.L * (fa[0] + 4 * fb[0] + fc[0]) / 6
        dy += -(1 / (4 * math.pi)) * fac * el.L * (fa[1] + 4 * fb[1] + fc[1]) / 6
    assert abs(shift[0] / dx - 1) < 5e-3 and abs(shift[1] / dy - 1) < 1e-2, (shift, dx, dy)


# ================================================================== 項目 9: 光子 (MULT・SEXT・OCT・BEAMBEAM + BEND・TRPT)
def _photons(src):
    _ffs, out = _run(src)
    A = [[float(x) for x in l.split()[1:]] for l in out.splitlines() if l.startswith('PA ')]
    B = [[float(x) for x in l.split()[1:]] for l in out.splitlines() if l.startswith('PB ')]
    return (np.hstack([np.array(A), np.array(B)]) if A else np.zeros((0, 12))), out


_PH_MULT_SRC = os.path.join(_ROOT, 'example', '07_spin_solenoid_special', 'flag_PHOTONS_mult.sad')
_PH_BB_SRC = os.path.join(_ROOT, 'example', '07_spin_solenoid_special', 'flag_PHOTONS_beambeam_bend.sad')


def test_photons_in_mult_match_gs_k64(k64_mode):
    """k64: BEND・QUAD・MULT (K0・K1・K2、ずれ・ROTATE つき) の光子 325 個の PhotonList が gs と一致する: 要素ごとの個数 (4: 66、10: 62、16: 67、24: 56、
    32: 69、34 (MULT): 3、36 (QUAD): 2)、E・位置・方向の列の和。E 列は K64 が γ でなく gammab = βγ を掛ける (相対 1/(2γ²) = 1.45e-8) ことまで再現する"""
    if not os.path.exists(_PH_MULT_SRC):
        pytest.skip('flag_PHOTONS_mult.sad')
    P, _out = _photons(open(_PH_MULT_SRC, encoding='utf-8').read())
    assert len(P) == 325
    cnt = {int(e): int((P[:, 11] == e).sum()) for e in set(P[:, 11])}
    assert cnt == {4: 66, 10: 62, 16: 67, 24: 56, 32: 69, 34: 3, 36: 2}, cnt
    gs = [1131316.469194416, 2624.2927423952, 1351.7752768154746, -0.003740676700538311, 237.3077858472, 184.05207215377834]
    for j, w in enumerate(gs):
        assert abs(P[:, j].sum() - w) <= 1e-9 * max(abs(w), 1.0) + 1e-9, (j, P[:, j].sum(), w)


def test_photons_beambeam_with_bend_and_trpt_match_gs_k64(k64_mode):
    """k64: BEAMBEAM と BEND (と QUAD) を同じ輸送路に置いた PHOTONS + TRPT が gs と一致する (以前は SADError)。光子 71 個 (BEAMBEAM 7、BEND 64)、
    追跡のあとの乱数の状態 (RND .5837396986 .9722504371) と座標の和 (TRK) まで一致する。TRPT の放射の分割 (ndivrad × 10) を track6 の経路でも効かせた"""
    if not os.path.exists(_PH_BB_SRC):
        pytest.skip('flag_PHOTONS_beambeam_bend.sad')
    P, out = _photons(open(_PH_BB_SRC, encoding='utf-8').read())
    assert len(P) == 71
    assert {int(e): int((P[:, 11] == e).sum()) for e in set(P[:, 11])} == {4: 7, 6: 64}
    rnd = _tag(out, 'RND ')
    assert abs(rnd[0] - 0.5837396986) < 2e-10 and abs(rnd[1] - 0.9722504371) < 2e-10, rnd
    trk = _tag(out, 'TRK ')
    assert abs(trk[0] + 0.0064064173) < 1e-9 and abs(trk[1] + 0.0041644990) < 1e-9 and abs(trk[2] + 3.015296742e-4) < 1e-12, trk
    for j, w in enumerate([109587.61024303322, 140.54692935429998, 0.519475588990291, 0.014846480557233598, 70.892933502, 3.284602775775321]):
        assert abs(P[:, j].sum() - w) <= 1e-7 * max(abs(w), 1.0), (j, P[:, j].sum(), w)


_GAM = 3e9 / 0.51099895e6
_HBARC = 1.973269804e-7           # [eV m]
_ALPHA = 1.0 / 137.035999084


def _thin_line(elem, n):
    return f'''MOMENTUM = 3 GEV;
DRIFT D = (L=1);
{elem}
MARK ST = (BX=5 BY=5);
LINE BL = (ST D XX D);
FFS;
USE BL; INS; NORAD; NOFLUC; CALC;
RAD; FLUC; PHOTONS;
SeedRandom[5];
n={n};
z0=Table[0.0,{{i,n}}]; fl0=Table[1,{{i,n}}];
r=TrackParticles[{{1,{{z0,z0,z0,z0,z0,z0,fl0}}}},1,1];
$FORM="16.10";
Do[Print["PA ",PhotonList[[k,1]]," ",PhotonList[[k,2]]," ",PhotonList[[k,3]]," ",PhotonList[[k,4]]," ",PhotonList[[k,5]]," ",PhotonList[[k,6]]],{{k,Length[PhotonList]}}];
Do[Print["PB ",PhotonList[[k,7]]," ",PhotonList[[k,8]]," ",PhotonList[[k,9]]," ",PhotonList[[k,10]]," ",PhotonList[[k,11]]," ",PhotonList[[k,12]]],{{k,Length[PhotonList]}}];
Exit;
'''


@pytest.mark.parametrize('elem,theta,l_el', [
    ('MULT XX = (L=0.4 K0=0.02);', 0.02, 0.4),
    ('SEXT XX = (L=0.2 K2=40 DX=1e-2);', 0.5 * 40 * 1e-4, 0.2),
    ('OCT XX = (L=0.1 K3=1.2e4 DX=1e-2);', 1.2e4 * 1e-6 / 6.0, 0.1),
])
def test_photons_from_mult_sext_oct_match_synchrotron_radiation_formulas_v2(v2_mode, elem, theta, l_el):
    """v2: 多極 (MULT の K0、DX だけずれた SEXT・OCT の曲げ角 θ = ½ K2 DX²、K3 DX³/6。SAD の K2・K3 は積分値) の光子の数・エネルギーが、偏向角 θ の一様な円弧 (ρ = L/θ) の
    放射の式と一致する: 光子数 N = (5/(2√3)) α γ θ per 粒子 (3 GeV。Poisson の 4σ)、平均 ⟨u⟩ = 8/(15√3) u_c (u_c = (3/2) ħ c γ³/ρ)、要素の番号 3。
    以前は記録されなかった"""
    n = 6000
    P, _out = _photons(_thin_line(elem, n))
    nexp = (5.0 / (2.0 * math.sqrt(3.0))) * _ALPHA * _GAM * theta * n
    assert abs(len(P) - nexp) < 4.0 * math.sqrt(nexp) + 0.02 * nexp, (len(P), nexp)
    assert set(P[:, 11].astype(int)) == {3}
    uc = 1.5 * _HBARC * _GAM**3 / (l_el / theta)
    mean = P[:, 0].mean() / uc
    assert abs(mean / (8.0 / (15.0 * math.sqrt(3.0))) - 1) < 5.0 * 0.76 / math.sqrt(len(P)) + 0.02, (mean, len(P))
    assert np.max(np.abs(P[:, 7])) <= 1.0 + 1e-9 and np.max(np.abs(P[:, 9])) <= 1.0 + 1e-9


# ================================================================== 項目 3: SynchroBetaEmittance の収束の警告
def test_synchrobeta_warns_when_azimuthal_modes_not_converged_v2(v2_mode):
    """(P2 で変更) v2: 大きなリング A (3 GeV、603 m、νs = 0.004) の AzimuthalModes->9 の emitx は、K64 の式の誤り (始点の分散 η ≠ 0。P2 項目 1) を直したので
    Emittance[] の 1.5% 以内 (1.0069 倍) で、14 の値と 1% 以内に収束し、警告を出さない (K64 では 5.1e-9 対 14 の 2.5e-9 で警告が出ていた)。PYSAD_SLOW=1 のときだけ (約 20 秒)"""
    if os.environ.get('PYSAD_SLOW') != '1':
        pytest.skip('PYSAD_SLOW=1 のときだけ (約 20 秒)')
    tail = '''
FFS; PageWidth=1999;
USE RING; CELL; RFSW; RAD; CALC;
em=Emittance[]; nus=Abs[(Tunes/.em)[[3]]];
r=SynchroBetaEmittance[{nus,nus,0.001},AzimuthalModes->9];
Print["SB ",r[[1,2]]];
Print["EM ",(Emittances/.em)[[1]]];
Exit;
'''
    _f, out = _run(open(_LAT_A, encoding='utf-8', errors='replace').read() + tail)
    assert 'SynchroBetaEmittance: AzimuthalModes' not in out
    assert abs(_tag(out, 'SB ')[0] / _tag(out, 'EM ')[0] - 1) < 0.015


def test_synchrobeta_does_not_warn_on_converged_ring_v2(v2_mode):
    """ASTRID2 (0.58 GeV) では AzimuthalModes->4 の SB が Emittance[] と 1% で、警告しない"""
    if not os.path.exists(_AST):
        pytest.skip(_AST)
    tail = '''
FFS; PageWidth=1999;
USE RING; CELL; RFSW; RAD; CALC;
em=Emittance[]; nus=Abs[(Tunes/.em)[[3]]];
r=SynchroBetaEmittance[{nus,nus,0.001},AzimuthalModes->4];
Print["SB ",r[[1,2]]];
Print["EM ",(Emittances/.em)[[1]]];
Exit;
'''
    _f, out = _run(open(_AST, encoding='utf-8', errors='replace').read() + tail)
    assert 'SynchroBetaEmittance: AzimuthalModes' not in out
    assert abs(_tag(out, 'SB ')[0] / _tag(out, 'EM ')[0] - 1) < 0.02


# ================================================================== 項目 10: K64 の spdepol の性質 (v2 の equpol_k64 = 桁あふれだけ直した K64 の式)
def test_k64_spdepol_rate_scales_as_emittance_over_damping_v2(v2_mode):
    """K64 の spdepol の脱偏極率 (4.2 GeV・12 MV の例で 4.3e-5 /周) は平衡エミッタンスに比例し **減衰率に反比例** する (減衰率を 0.5 倍・2 倍にすると
    1.95 倍・0.53 倍)。正しい脱偏極率は光子の跳びの拡散 D (減衰率に依らない) で決まる ½ Σ tr(N B Nᵀ) = 2.3e-8 /周 (線形)。K64 の式は
    2 次の摂動の零位相の項で 1/(1 − e^{−2λ}) (= 平衡エミッタンスが含む 1/(2λ) の二重数え) を掛けた形で、1000 倍以上大きい。docs/history/v2next_D1.md §6"""
    from pysad import spin
    src = open(os.path.join(_ROOT, 'example', '16_spin_polarization', 'pol_sideband.sad'), encoding='utf-8').read().split('FFS;')[0]
    ffs, _ = ffscmd.run_sad_string(src + 'FFS;\nUSE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n', out=io.StringIO())
    PST = 8.0 / (5.0 * math.sqrt(3.0))
    orig = spin.spdepol
    rates = {}
    lin = None
    try:
        for f in (0.5, 1.0, 2.0):
            def w(gxr, gxi, gyr, gyi, gzr, gzi, e1, e2, dx, ssprd, amx, ams, f=f):
                return orig(gxr, gxi, gyr, gyi, gzr, gzi, e1, e2, np.asarray(dx) * f, ssprd, amx, ams)
            spin.spdepol = w
            o = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
            sp = o['spin']
            rates[f] = sp['sdamp'] * (PST / sp['equpol_k64'][0] - 1.0)
            lin = sp['depol_rate_linear']
    finally:
        spin.spdepol = orig
    assert 4.0e-5 < rates[1.0] < 4.6e-5, rates
    assert abs(rates[0.5] * 0.5 / rates[1.0] - 1) < 0.1 and abs(rates[2.0] * 2.0 / rates[1.0] - 1) < 0.1, rates
    assert rates[1.0] / lin > 1000.0, (rates, lin)


# ================================================================== 項目 8: 不安定なリングの POL の Emittance[]
def test_pol_emittance_of_unstable_ring_returns_nan_without_exception_v2(v2_mode):
    """v2: 光学が不安定で閉軌道が見つからないリング (四極が強すぎる) の POL ON の Emittance[] は、30 回の Newton 法・NaN の 1 周・要素ごとの
    スピンの計算をやらずに NaN を返す (SP8 (大きなリング B、8 GeV、1998 要素) の四極 60 台に DX・DY 0.1 mm rms: 以前は 20 分以上、今は 数分以内)。
    ここでは小さなリングで例外にならないことを確かめる"""
    src = _FODO % ('1', '3.0', 0.0, 'RFSW; RAD; NOFLUC; POL;')
    _f, out = _run(src.replace('Print["NR ", Emittances/.e, Tunes/.e, DampingRate/.e];', 'Print["NR ", Polarization/.e];'))
    assert 'NR' in out


# ================================================================== 項目 5: POL + TRPT + 加速 (track6 の経路)
def _track_spin(ffs, c, k0, k1, spin):
    c = np.asarray(c, float)
    n = len(c)
    cols = [list(c[:, j]) for j in range(6)]
    ffs.flags['RAD'] = True
    ffs.flags['FLUC'] = False
    ffs.flags['POL'] = True
    try:
        sp = np.asarray(spin, float)
        out = ffs.fn_track([[float(k0), cols + [[math.atan2(s[2], s[0]) for s in sp], [s[1] for s in sp], [1.0] * n]], float(k1)])
    finally:
        ffs.flags['POL'] = False
    r = out[1]
    return np.array(r[:6]).T, np.array([[math.sqrt(max(0.0, 1 - y * y)) * math.cos(a), y, math.sqrt(max(0.0, 1 - y * y)) * math.sin(a)]
                                        for a, y in zip(r[6], r[7])])


_ACC = '''MOMENTUM = 1 GEV;
DRIFT D = (L=0.5);
BEND B1 = (L=1.0 ANGLE=0.3 E1=0 E2=0) B2 = (L=1.0 ANGLE=0.2 E1=0 E2=0);
CAVI RF = (L=0 VOLT=%g FREQ=5e8 PHI=-PI/2);
MARK ST = (BX=5 BY=5);
LINE BL = (ST D B1 D RF D B2 D);
FFS;
USE BL;
TRPT; INS; RFSW; RAD; NOFLUC; CALC;
Exit;
'''


@pytest.mark.parametrize('volt', [0.0, 2e8, 5e8])
def test_spin_tracking_with_trpt_acceleration_equals_adiabatic_precession_v2(v2_mode, volt):
    """v2: POL + RAD + TRPT + RFSW の加速のある CAVI (track6 の経路。以前は SADError) の後、スピンの歳差は断熱加速の解析式
    Φ = aγ1 θ1 + aγ2 θ2 (γ2 = γ1 + VOLT/(mc²)、a = 0.00115965) になる: 1 GeV、θ1 = 0.3・θ2 = 0.2 の扇形の偏向磁石の間に VOLT = 0〜0.5 GV の空洞。
    縦向きに始めたスピンの (sx, sz) = (sin Φ, cos Φ) と 5e-7 以内 (VOLT = 0.5 GV: Φ = 1.3616 rad、加速なしなら 1.1347 rad)"""
    ffs = ffscmd.run_sad_string(_ACC % volt, out=io.StringIO())[0]
    c, s = _track_spin(ffs, [[0, 0, 0, 0, 0, 0]], 1, len(ffs.els) + 1, [[0, 0, 1]])
    ph = emit6d._Phys(ffs.momentum(), 10.0)
    a = 0.00115965218
    h1 = ph.h0
    h2 = h1 + volt / ph.amass
    phi = a * (h1 * 0.3 + h2 * 0.2)
    assert abs(abs(s[0, 0]) - abs(math.sin(phi))) < 5e-7 and abs(s[0, 2] - math.cos(phi)) < 5e-7, (s, phi)
    assert abs(np.linalg.norm(s[0]) - 1.0) < 1e-9


def test_spin_tracking_with_trpt_acceleration_is_still_an_error_in_k64(k64_mode):
    """k64: POL + TRPT + 加速 (track6 の経路) は従来どおり明示のエラー (K64 の動作ではなく pysad の未対応。v2 は上のとおり対応)"""
    ffs = ffscmd.run_sad_string(_ACC % 2e8, out=io.StringIO())[0]
    with pytest.raises(Exception):
        _track_spin(ffs, [[0, 0, 0, 0, 0, 0]], 1, len(ffs.els) + 1, [[0, 0, 1]])
