# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 3: k64 モードの PSPAC (K64 の spch.c の形状関数まで直訳: pysad/spch_k64.py)。参照は gs (K64、2026-10-19 ビルド、2026-11-19)。

K64 の PSPAC: 電荷は横が CIC (双線形)、縦が「ほぼ最近接」(rhocal)、電位は 5 点差分のポアソン方程式 (FACR) の解を双 3 次補間 (makebcu・epotkick)、
縦の線密度は 3 次スプライン。k64 モードは gs と全桁 (rms の 12 桁以上) 一致する。v2 は別の形状関数 (2 次の B スプライン) で、メッシュが十分細かければ
解析式と 1〜2% で合うが、粗いメッシュ (縦・横が数セル) の値は gs と違う (v2 の値は下の表。`test_v2_A1_next.py` の警告の試験も参照)。
使った条件: 5 MeV の電子、1e11 個、半径 1 mm の一様な円 (準乱数 4096 個)、`D(0.5) SPCH D(0.5) SPCH`、PSPACNX = PSPACNY = 128、NPARA = 1
(gs の既定の NPARA は 1 より大きく、TrackParticles は粒子を子プロセスに分け、各プロセスは自分の粒子だけで電位を作るので、gs と比べるときは NPARA = 1)。
未確認: 実機のリング・多周回 (3・4 周は確かめた)、SPAC・WSPAC との併用。gs は 1 セッションの 2 回目の TrackParticles で PSPAC が働かない (SPC8) ので、
gs と比べるときは 1 セッション 1 回。
"""
import io
import re

import pytest

_TPL = '''MOMENTUM = 5MEV;
DRIFT D1 = (L=0.5);
SPCH SP1 = ();
MARK ORIGIN = ();
LINE TL = (ORIGIN D1 SP1 D1 SP1);
FFS;
NPARA = 1; RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
PSPACNX = %(nxy)d; PSPACNY = %(nxy)d; PSPACDX = %(dx)s; PSPACDY = %(dx)s; PSPACNZ = %(nz)d; PSPACDZ = %(dz)s;
PSPACNTURN = %(pn)d; PSPACNTURNCALC = %(pc)d;
a = 1e-3; n = 4096; lb = %(lb)s;
u = Table[Mod[i*0.6180339887, 1], {i,1,n}];
v = Table[Mod[i*0.7548776662, 1], {i,1,n}];
w = Table[Mod[i*0.5698402910, 1], {i,1,n}];
rms[l_] := Sqrt[Plus@@(l^2)/Length[l]];
$FORM = "20.14";
PSPAC;
bm = {a*Sqrt[u]*Cos[2*Pi*v], Table[0,{n}], a*Sqrt[u]*Sin[2*Pi*v], Table[0,{n}], lb*(w-0.5), Table[0,{n}], Table[1,{n}]};
r = TrackParticles[{1,bm},1,1,%(nturn)d];
Print["RMS ", rms[r[[2,1]]], " ", rms[r[[2,3]]], " ", rms[r[[2,2]]], " ", rms[r[[2,4]]], " ", Plus@@(r[[2,6]])/n, " ", rms[r[[2,6]]]];
Exit;
'''


def _run(**kw):
    from pysad.ffscmd import run_sad_string
    d = dict(nxy=128, dx='1e-4', nz=32, dz='0.005', pn=1, pc=0, lb='0.1', nturn=1)
    d.update(kw)
    buf = io.StringIO()
    run_sad_string(_TPL % d, out=buf)
    m = re.search(r'^RMS\s+(.*)$', buf.getvalue(), re.M)
    return [float(t) for t in m.group(1).split()]


# gs の値 (rms x, rms y, rms px, rms py, 平均 δ, rms δ)。ps1・ps2・ps3: lb = 0.1・0.02・0.01 m (縦 20・4・2 セル)、ps4: lb = 0.02、DZ = 1 mm (20 セル)、
# DX = 2.5e-4、ps5: 64×64、DX = 4e-4、DZ = 0.01、16 ノード、lb = 0.05 m
_GS = {
    'ps1': (dict(), [.00123400976450, .00123535491635, .00207032904309, .00207111530278, 1.818466684475E-5, 3.604206917924E-4]),
    'ps2': (dict(lb='0.02'), [.00391558304726, .00391696207641, .00723196098521, .00723305416468, 3.049772047593E-4, .00261404512982]),
    'ps3': (dict(lb='0.01'), [.00665638860781, .00665701817226, .01236208778810, .01236114013699, .00108651085631, .00721425190019]),
    'ps4': (dict(lb='0.02', dz='0.001', dx='2.5e-4'), [.00416592073137, .00416804115983, .00821562569900, .00821789421427, 4.844090086514E-4, .00974309896641]),
    'ps5': (dict(lb='0.05', nz=16, dz='0.01', dx='4e-4', nxy=64), [.00190864985013, .00191085239796, .00359783908812, .00360051409731, 5.480362947328E-5, .00118284299365]),
}


@pytest.mark.parametrize('name', sorted(_GS))
def test_k64_pspac_matches_gs_to_all_digits(k64_mode, name):
    kw, ref = _GS[name]
    got = _run(**kw)
    for g, r in zip(got, ref):
        assert g == pytest.approx(r, rel=3e-12, abs=1e-18), (name, got, ref)


# v2 の値 (同じ条件): 解析式の側の PSPAC (v2) は粗いメッシュで gs から離れる。lb = 0.1 (20 セル) は 0.6%、0.02 (4 セル) は 4.5%、0.01 (2 セル) は 12%
_V2 = {'ps1': [.00122696303226, .00122819881575], 'ps2': [.00373747598022, .00373842773943], 'ps3': [.00586189643299, .00586242099212]}


@pytest.mark.parametrize('name', sorted(_V2))
def test_v2_pspac_differs_from_gs_at_coarse_mesh_as_recorded(v2_mode, name):
    kw, ref = _GS[name]
    got = _run(**kw)
    assert got[:2] == pytest.approx(_V2[name], rel=2e-9)
    d = abs(got[0] / ref[0] - 1)
    assert {'ps1': 4e-3 < d < 8e-3, 'ps2': 0.03 < d < 0.06, 'ps3': 0.10 < d < 0.14}[name], d


# 多周回: PSPACNTURN = 3 で毎周電位を作り直す (PSPACNTURNCALC = 0)、PSPACNTURNCALC = 1 は 1 周目の電位を凍結 (λ は最後に計算した SPCH のものを全部で使う)、
# PSPACNTURN = 4・PSPACNTURNCALC = 2。PSPACNTURN = 1 のまま 2 周追うと gs は解放済みのメモリ (凍結した電位) を読み、凍結と同じ値
_GS_TURNS = {
    (3, 3, 0): [.00615192755698, .00615448379681, .00288664003471, .00288747505634, 4.536212853657E-5, 8.748833344156E-4],
    (3, 3, 1): [.00693066341280, .00693348158573, .00370934603264, .00370864745189, 5.464756049201E-5, 9.996393068090E-4],
    (4, 4, 2): [.00936903200871, .00937188603884, .00332556750639, .00332549014551, 5.952033105547E-5, .00107022648223],
    (2, 1, 0): [.00362801211377, .00362993127980, .00318843519373, .00318951594439, 4.384106347301E-5, 8.297625755486E-4],
}


@pytest.mark.parametrize('nturn,pn,pc', sorted(_GS_TURNS))
def test_k64_pspac_multi_turn_and_frozen_potential_match_gs(k64_mode, nturn, pn, pc):
    got = _run(nturn=nturn, pn=pn, pc=pc, dx='2e-4')
    for g, r in zip(got, _GS_TURNS[(nturn, pn, pc)]):
        assert g == pytest.approx(r, rel=3e-12, abs=1e-18)


def test_k64_pspac_nz_below_3_is_an_error(k64_mode):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_TPL % dict(nxy=128, dx='1e-4', nz=1, dz='0.005', pn=1, pc=0, lb='0.1', nturn=1), out=buf)
    assert 'PSPACNZ' in buf.getvalue()


def test_k64_shape_functions_against_the_c_code_building_blocks():
    """K64 の補助関数: 3 次スプライン (splined) は端の傾き 0 の自然な形で、cos の 2 階微分を与える。ポアソン方程式は 5 点差分を厳密に満たす。
    双 3 次補間の係数は 2 次式の電位を厳密に再現する (中心差分の微分が厳密になるので)"""
    import numpy as np
    from pysad.spch_k64 import K64Pspac, splined
    x = np.linspace(0.0, 1.0, 41)
    y = np.sin(2.0 * np.pi * x) ** 2
    y2 = splined(x, y)
    mid = slice(5, 36)
    assert np.abs(y2[mid] - (8 * np.pi ** 2 * np.cos(4 * np.pi * x))[mid]).max() < 5e-2 * 8 * np.pi ** 2
    K = K64Pspac(32, 32, 8, 0.1, 0.1, 0.1, 1e10, 0.511e6)
    rng = np.random.default_rng(3)
    rho = np.zeros((33, 33))
    rho[1:32, 1:32] = rng.random((31, 31))
    phi = K._poisson(rho)
    lap = (phi[2:, 1:-1] + phi[:-2, 1:-1] + phi[1:-1, 2:] + phi[1:-1, :-2] - 4 * phi[1:-1, 1:-1]) / 0.01
    assert np.abs(lap + 2 * np.pi * rho[1:-1, 1:-1] / 0.01).max() < 1e-9 * 2 * np.pi * 100
    I, J = np.meshgrid(np.arange(33), np.arange(33), indexing='ij')
    q = 3.0 + 0.7 * I - 0.2 * J + 0.05 * I * I + 0.03 * I * J - 0.04 * J * J
    B = K._bcu(q)
    c = B[10, 12]
    assert c[0] == pytest.approx(q[10, 12], abs=1e-12) and c[4] == pytest.approx(0.7 + 0.1 * 10 + 0.03 * 12, abs=1e-10)
