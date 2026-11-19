# -*- coding: utf-8 -*-
"""担当 C2 の光学まわり: N3 (CalculateOptics の TRPT の加速)、加速区間の分散 (qmat2twiss の gr の向き)。
参照は gs (SAD 1.4.4.2k64、2026-10-19 ビルド、WSL) の出力。どちらも pysad の不具合の修正で、モードによらない。"""
import io
import re

import pytest


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _nums(text, tag):
    m = re.search(tag + r'\s*(\{.*?\}\})', text.replace('\n', ''))
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', m.group(1))]


_LAT = '''MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 EX=0.1);
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
CAVI CA = (L=0.5 VOLT=2E8 FREQ=1.3E9 PHI=1.2) CB = (L=0 VOLT=2E8 FREQ=1.3E9 PHI=1.2);
LINE TL = (M0 QF D1 CA D1 QD D1 CA D1 QF D1 QD D1) T2=(M0 D1 CB D1) T3=(M0 D1 CA D1);
FFS;
TRPT; RFSW; INS;
'''

# gs: Twiss[{"BX","EX","EPX"},"*"] (M0 D1 CB D1 と M0 D1 CA D1。電子なので PHI=1.2 は減速: p は 1/1.2291 倍)
_GS_T2 = [8, 8, 7.15625, 5.822268844140862, 5.326814072797931, .1, .1, .1, .0733854778922643, .07338544900614298,
          0, 0, 0, -2.888612131518278e-8, -2.888612131518278e-8]
_GS_T3 = [8, 8, 7.15625, 6.783631774519875, 6.134412906367176, .1, .1, .1, .08096106142091347, .07922533591269021,
          0, 0, 0, -.0017357255082232632, -.0017357255082232632]
# gs: 輸送路 TL の終端の BX AX BY AY EX EPX NX NY
_GS_TL = [20.241775326529364, -4.547386967584856, 12.782565754236398, 1.2916215522778773, -.04217888436136674,
          -.014576829008643814, 2.42224243409817, .8747171317907109]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_dispersion_through_accelerating_cavity_matches_gs(mode):
    """加速 (減速) する空洞を通る分散。qmat2twiss の gr = gammab(出口)/gammab(入口)。
    薄い空洞では x は変わらず δ = Δp/p が p0/p1 倍になるので、分散は p1/p0 倍 (× SAD の座標の規格化 sqrt(p1/p0)):
    0.1 → 0.1·(p1/p0)^1.5 = 0.0733855"""
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(_LAT + 'USE T2; CALC; Print["T2 ",Twiss[{"BX","EX","EPX"},"*"]];\n'
                   'USE T3; CALC; Print["T3 ",Twiss[{"BX","EX","EPX"},"*"]];\nExit;\n')
    finally:
        compat.set_mode(old)
    for tag, ref in (('T2', _GS_T2), ('T3', _GS_T3)):
        v = _nums(out, tag)
        assert len(v) == len(ref)
        for a, b in zip(v, ref):
            assert abs(a - b) <= 1e-11 * max(1.0, abs(b)), (tag, a, b)
    g = 5.822268844140862 / 7.15625               # 薄い空洞の β の比 = p1/p0
    assert abs(_GS_T2[8] - 0.1 * g ** 1.5) < 1e-9


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_calculate_optics_includes_trpt_acceleration(mode):
    """N3: CalculateOptics (cell=0) は CALC と同じ結果 (TRPT かつ RFSW の加速を入れる)。区間の途中からでも同じ"""
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        out = _out(_LAT + '''USE TL; CALC;
n=LINE["LENGTH"];
k={"BX","AX","BY","AY","EX","EPX","NX","NY"};
t0=Twiss[k,n];
r=CalculateOptics[1,n,Twiss["*",1],0,2];
t1=Twiss[k,n];
r=CalculateOptics[5,n,Twiss["*",5],0,2];
t2=Twiss[k,n];
Print["RES ",{t0,{Max[Abs[t0-t1]],Max[Abs[t0-t2]]}}];
Exit;
''')
    finally:
        compat.set_mode(old)
    v = _nums(out, 'RES')
    for a, b in zip(v[:8], _GS_TL):
        assert abs(a - b) <= 1e-11 * max(1.0, abs(b))
    assert v[8] < 1e-13 and v[9] < 1e-13


# ---------------------------------------------------------------- RADTAPER (項目 4)
_RT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
$FORM = "18.12";
PageWidth=1999;
USE RINGQ;
CELL;
CALCULATE;
RFSW; RAD; RADCOD;
'''
_TRK = ('b = {{0},{0},{0},{0},{0},{0},{1}};\nDo[b = TrackParticles[{1, b}, 1, 1, 1][[2]]; ,{2}];\n'
        'Print["%s = ", b[[1,1]], " ", b[[2,1]], " ", b[[5,1]], " ", b[[6,1]]];\n')


def _line(out, tag):
    l = [x for x in out.splitlines() if x.startswith(tag + ' =')][0]
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', l.split('=', 1)[1])]


def test_radtaper_closed_orbit_is_the_converged_gs_solution(k64_mode):
    """RADTAPER の閉軌道。gs の既定の収束判定 (tcod: Σ(Δ/codw)² < 1e-4、codw_x = 1e-6 m → 軌道の残り 1e-8 m) では
    x = −2.19771e-7、δ = 3.5e-9 だが、`CODCONV = 1e-4` で判定を 1 万倍きつくすると x = −2.17319734889e-7、δ = 5.8e-13、
    U0 = 1564904.008828806、DDP(60) = 3.911239587615e-4 になる。pysad (radtaper.solve、判定 1e-10) はこの収束した解と合う。
    閉軌道の z は trf0 に移す (gs の収束した解は z = −4.8e-15、trf0 = 0.038358121964): これで追跡 (z は trf0 からの距離) も合う。
    gs (CODCONV=1e-4) の追跡 2 周: FLUC {4.808978204006e-5, 1.447636324091e-5, −3.69156359875e-5, 8.005891925389e-5}、
    NOFLUC {5.900361498408e-6, 1.929885522281e-6, −2.03579883126e-5, −1.25856551223e-7}、
    その後 NORADTAPER: FLUC x = 3.276317389820e-4、NOFLUC x = 2.894235513086e-4"""
    out = _out(_RT + 'RADTAPER;\ne1 = Emittance[];\n'
               'Print["B = ", (OrbitAtExit/.e1)[[1]], " ", (OrbitAtExit/.e1)[[6]], " ", EnergyLossU0/.e1, " ", (Emittances/.e1)[[1]], " ",'
               ' (Emittances/.e1)[[3]], " ", (Tunes/.e1)[[1]], " ", (Tunes/.e1)[[2]], " ", (Tunes/.e1)[[3]]];\n'
               'Print["Z = ", (OrbitAtExit/.e1)[[5]], " ", EquilibriumPosition/.e1];\n'
               'Print["DDP = ", Twiss["DDP", {1, 60, 121, 243}]];\n'
               + _TRK % 'T' + 'Exit;\n')
    out2 = _out(_RT + 'RADTAPER;\ne1 = Emittance[];\nNOFLUC;\n' + _TRK % 'TN' + 'FLUC; NORADTAPER;\n' + _TRK % 'TF0'
                + 'NOFLUC;\n' + _TRK % 'TN0' + 'Exit;\n')
    b, z, d = _line(out, 'B'), _line(out, 'Z'), _line(out, 'DDP')
    assert abs(b[0] / -2.17319734889e-7 - 1) < 1e-5 and abs(b[1]) < 1e-10
    assert abs(b[2] / 1564904.008828806 - 1) < 1e-9
    assert abs(b[3] / 1.737161355738e-7 - 1) < 1e-9 and abs(b[4] / 1.907727704002e-5 - 1) < 1e-8
    for u, v in zip(b[5:8], [-.303633024754, -.094058498417, -.027618670292]):
        assert abs(u - v) < 1e-10
    assert abs(z[0]) < 1e-12 and abs(z[1] - 0.038358121964) < 2e-9
    assert d[0] == 0.0 and abs(d[1] / 3.911239587615e-4 - 1) < 1e-8 and abs(d[2] / 2.606813180711e-4 - 1) < 1e-8
    for tag, ref, tol in (('T', [4.808978204006e-5, 1.447636324091e-5, -3.69156359875e-5, 8.005891925389e-5], 1e-5),
                          ('TN', [5.900361498408e-6, 1.929885522281e-6, -2.03579883126e-5, -1.25856551223e-7], 1e-4),
                          ('TF0', [3.276317389820e-4, 7.150983414142e-5, -.001231553142, 6.012899757899e-5], 1e-7),
                          ('TN0', [2.894235513086e-4, 6.085548893954e-5, -.001220333649, -7.89905013383e-6], 1e-6)):
        for u, v in zip(_line(out if tag == 'T' else out2, tag), ref):
            assert abs(u / v - 1) < tol, (tag, u, v)


def test_radtaper_solution_is_a_fixed_point_and_tolerance_independent(v2_mode):
    """radtaper.solve の解は不動点: 係数を作り直しても δ の分布が変わらず (1e-12)、判定 tol を 1e-10 → 1e-13 にしても
    δ の分布が 1e-11、一様な係数 u が 2e-10 以内で同じ (gs の既定の判定での差 2.5e-9 m は収束の残り)"""
    import numpy as np
    from pysad import radtaper as rt, emit6d as e
    from pysad.ffscmd import run_sad_string
    ffs, I = run_sad_string(_RT + 'Exit;\n', out=io.StringIO())
    r1 = rt.solve(ffs.els, ffs.momentum())
    r2 = rt.solve(ffs.els, ffs.momentum(), tol=1e-13, maxit=60)
    assert r1['found'] and r2['found']
    assert np.max(np.abs(r1['ddp'] - r2['ddp'])) < 1e-11 and abs(r1['u'] - r2['u']) < 2e-10
    ph = e._Phys(ffs.momentum(), e._emit_circ(ffs.els))
    elt = rt.tapered(ffs.els, r2['shape'], r2['u'])
    cod, trf0, fnd, _ = e.closed_orbit6d(elt, ph, True, True)
    prof = rt._profile(elt, cod, ph, trf0)
    assert fnd and np.max(np.abs(prof - prof[0] - r2['shape'])) < 1e-12 and abs(prof[0]) < 1e-12
    assert abs(prof[-1] - prof[0]) < 1e-9                     # 1 周で元の δ に戻る (閉軌道の探索の判定の範囲)


def test_radtaper_in_calc6d(v2_mode):
    """v2: CELL + CALC6D の CALC にも RADTAPER が入る (Emittance[] と同じ閉軌道)。K64 (gs) は 'qcell6d---> Unstable optics' を出し、
    入口の軌道がテーパーなしの値 (−7.887e-5) のままで動作しない。4 次元の CALC は放射を入れないのでテーパーは無い (K64 も v2 も)"""
    from pysad import compat
    src = (_RT + 'CALC6D;\nCALC;\nPrint["CA = ", Twiss["DX",1], " ", Twiss["DDP",60]];\nRADTAPER;\nCALC;\n'
           'Print["CB = ", Twiss["DX",1], " ", Twiss["DDP",1], " ", Twiss["DDP",60], " ", Twiss["DDP",121], " ", Twiss["DZ",1]];\n'
           'e1 = Emittance[];\nPrint["E = ", (OrbitAtExit/.e1)[[1]], " ", Twiss["DDP",60]];\nExit;\n')
    out = _out(src)
    ca, cb, em = _line(out, 'CA'), _line(out, 'CB'), _line(out, 'E')
    assert abs(ca[0] / -7.887e-5 - 1) < 1e-3
    assert abs(cb[0] / em[0] - 1) < 1e-6 and abs(cb[0] / -2.1732e-7 - 1) < 1.5e-2       # テーパーで軌道が設計軌道に戻る (v2: tradke の ddpx の符号 (A1 項目 7) で gs の値から 3e-3、BEND の係数を (δ_l + δ_l+1)/2 にした (P2 項目 3) ので 1.1e-2 ずれる)
    assert abs(cb[1]) < 1e-9 and abs(cb[2] / 3.911239e-4 - 1) < 1e-5 and abs(cb[3] / 2.606813e-4 - 1) < 1e-5 and abs(cb[4]) < 1e-12
    assert 'Unstable' not in out
    compat.set_mode('k64')
    try:
        out = _out(src)
    finally:
        compat.set_mode('v2')
    cbk = _line(out, 'CB')
    assert abs(cbk[0] / -2.1736e-7 - 1) < 1e-3                # k64 も CALC6D に RADTAPER を入れる (v2next B1: gs で確かめた。gs の DX(1) = -2.46e-7、チューン 6.696366849 と 2e-9 で一致)
