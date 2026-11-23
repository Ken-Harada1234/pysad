# -*- coding: utf-8 -*-
"""担当 C3 (v2next): MEASURE・TRCOD・ALIGN・FIXCOD の独立な検証 (大きなリング)。gs は MEASURE が無出力・TRCOD が SIGSEGV (B42) で比較できない。

独立な基準:
- 閉軌道の測定: 薄いキック θ が 1 か所にあるリングの閉軌道は、解析式 x(s) = θ √(β_s β_k) / (2 sin πν) · cos(πν − |φ_s − φ_k|) (Twiss のβ・位相から)。
  TRCOD が測った (粒子の追跡の重心の) 軌道が、CALC の閉軌道と、この式と一致するか。
- MEASURE の初期分布の統計: 粒子数を増やして、設定した EMITX・EMITY の固有エミッタンス (Σ·J の固有値) と、各要素の σx = √(εx βx + (ηx δ)²) に一致するか。
- ALIGN・FIXCOD: 揃えた磁石の DX・DY の残差 0、基準軌道との差が、キックの変化に対する解析式と一致するか。補正 (キックを戻す) 後の残差。
条件: PF リング (2.5 GeV、周長 187 m、`example/pf/lattice/v9_16_asym.sad`)、ASTRID2 (0.58 GeV、周長 45.7 m)、
大きなリング A (3 GeV、周長 603 m、`kek-ls_v03_12.lat`。TRCOD は 35 s かかるので `PYSAD_SLOW=1` のときだけ走る)。ラティスの最初に MARK と K0 のキック (L=0 の BEND) を足す。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, ffscmd, measure as ms
from pysad.sadrandom import SadRandom

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RINGS = {
    'PF': os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'),
    'AST': os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'),
    'A': os.path.join(HERE, '..', 'kek-ls_v03_12.lat'),
}


def _src(name):
    path = RINGS[name]
    if not os.path.exists(path):
        pytest.skip(os.path.basename(path) + ' がない')
    with open(path, encoding='utf-8') as fh:
        return fh.read()


def _kick_ring(name, th, extra=''):
    return _src(name) + f'''
MARK ORG0 = ();
BEND KICK = (L=0 ANGLE=0 K0={th});
LINE RINGK = (ORG0 KICK RING);
FFS;
USE RINGK; CELL; CALC;
{extra}
Exit;
'''


def _analytic_orbit(ffs, th):
    rows = ffs.ensure_rows()
    names = [e.name for e in ffs.lat.elements]
    ik = names.index('KICK')
    bx = np.array([r.betx for r in rows])
    mux = np.array([r.mux for r in rows])
    nu = mux[-1] / 2 / math.pi
    f = np.sqrt(bx * bx[ik]) / (2 * math.sin(math.pi * nu)) * np.cos(math.pi * nu - np.abs(mux - mux[ik]))
    return -th * f                       # SAD の K0 の符号: 正のキックは −x 側 (符号は 2 つの向きで当てて確かめた)


def _orbit_x(ffs):
    return np.array([r.x_cod for r in ffs.ensure_rows()])


@pytest.mark.parametrize('name,tol_an,tol_trcod', [('AST', 1e-3, 1e-5), ('PF', 2e-3, 1e-5),
                                                   pytest.param('A', 3e-3, 5e-5, marks=pytest.mark.skipif(not os.environ.get('PYSAD_SLOW'), reason='PYSAD_SLOW=1 (35 s)'))])
def test_trcod_measures_the_closed_orbit_of_a_kick_v2(v2_mode, name, tol_an, tol_trcod):
    """TRCOD (NP=2、ε=0、DP=0: 重心の追跡) で測った閉軌道 = CALC の閉軌道 (tol_trcod 以内、最大値で規格化) = 解析式 (tol_an 以内。6 極のフィードダウンで 1e-3 のずれ)。
    実測: ASTRID2 4.3e-4・1.8e-6、PF 6.4e-4・6e-7、大きなリング A 1.9e-3・1.3e-5"""
    th = 1e-5
    ffs, _ = ffscmd.run_sad_string(_kick_ring(name, th), out=io.StringIO())
    x_calc = _orbit_x(ffs)
    xa = _analytic_orbit(ffs, th)
    assert np.abs(x_calc - xa).max() / np.abs(xa).max() < tol_an
    ffs2, _ = ffscmd.run_sad_string(_kick_ring(name, th, 'NP=2; DP=0; EMITX=0; EMITY=0; SIGZ=0;\nTRCOD;'), out=io.StringIO())
    x_trcod = _orbit_x(ffs2)
    assert np.abs(x_trcod - x_calc).max() / np.abs(x_calc).max() < tol_trcod
    assert np.abs(x_trcod - xa).max() / np.abs(xa).max() < tol_an


def test_fixcod_ddx_equals_calc_difference_v2(v2_mode):
    """FIXCOD の基準 (キック 1e-5) と、キック 3e-5 の CALC の Twiss["DDX"] (現在 − 基準) が、2 つの CALC の DX の差と 1e-12 以内。
    元のキックに戻すと DDX・DDPX の最大値は 5e-8 以下 (補正後の残差。CALC の閉軌道の収束の精度で、軌道の最大値の 2e-4)"""
    th = 1e-5
    src = _src('PF') + f'''
MARK ORG0 = ();
BEND KICK = (L=0 ANGLE=0 K0={th});
LINE RINGK = (ORG0 KICK RING);
FFS;
$FORM="20.14"; PageWidth=1999;
USE RINGK; CELL; CALC;
TRCOD FIXCOD;
Print["B ", Twiss["DX", 100], " ", Twiss["DX", 300]];
KICK K0 {3 * th}; CALC;
Print["N ", Twiss["DX", 100], " ", Twiss["DX", 300], " ", Twiss["DDX", 100], " ", Twiss["DDX", 300]];
KICK K0 {th}; CALC;
Print["R ", Max[Abs[Table[Twiss["DDX", i], {{i, 1, 400}}]]], " ", Max[Abs[Table[Twiss["DDPX", i], {{i, 1, 400}}]]]];
Exit;
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    v = {l.split()[0]: [float(q) for q in l.split()[1:]] for l in buf.getvalue().splitlines() if l[:2] in ('B ', 'N ', 'R ')}
    for k in (0, 1):
        assert abs(v['N'][2 + k] - (v['N'][k] - v['B'][k])) < 1e-12
    assert v['R'][0] < 5e-8 and v['R'][1] < 5e-8                # CALC の閉軌道の収束の精度 (実測 1.6e-8 m。軌道の最大値 8e-5 m の 2e-4)
    assert abs(v['N'][2]) > 1e-6                              # 差は有限 (キックの変化が測れている)


def test_align_copies_dx_dy_exactly_v2(v2_mode):
    """TRCOD ALIGN 要素1 要素2: 区間の BEND QUAD SEXT OCT DECA DODECA MULT の DX・DY が要素1 の値に揃う (残差 0)。PF リング、QF を DX = 0.1 mm・DY = −0.2 mm にして
    次の SEXT・四極・偏向磁石 (要素2 の手前まで) に写す"""
    src = _src('PF')
    names = [w for w in __import__('re').findall(r'^\s*(?:QUAD|SEXT|BEND)\s+(\w+)\s*=', src, __import__('re').M)]
    assert len(names) > 5
    body = '''
FFS;
$FORM="20.14"; PageWidth=1999;
USE RING; CELL;
'''
    ffs, _ = ffscmd.run_sad_string(src + body + 'Exit;\n', out=io.StringIO())
    els = ffs.lat.elements
    idx = [i for i, e in enumerate(els) if type(e).__name__ in ('Quadrupole', 'Sextupole', 'Bend') and e.L > 0]
    def _pos(i):                       # 同じ名前の要素はファミリー: NAME.n (n 番目) で指定する
        n = sum(1 for e in els[:i + 1] if e.name == els[i].name)
        return f'{els[i].name}.{n}'
    ia, ib = idx[3], idx[12]
    script = (src + body + f'LINE["DX", "{_pos(ia)}"] = 1e-4; LINE["DY", "{_pos(ia)}"] = -2e-4;'
              f'TRCOD ALIGN {_pos(ia)} {_pos(ib)};'+chr(10))
    rng_idx = [i for i in idx if ia < i < ib]
    prints = ''.join(f'Print["V {i} ", {{LINE["DX", "{_pos(i)}"], LINE["DY", "{_pos(i)}"]}}];' + chr(10) for i in rng_idx[:6])
    buf = io.StringIO()
    ffscmd.run_sad_string(script + prints + 'Exit;'+chr(10), out=buf)
    got = [l for l in buf.getvalue().splitlines() if l.startswith('V ')]
    assert len(got) >= 3
    for l in got:
        nums = [float(q) for q in l.replace('{', ' ').replace('}', ' ').replace(',', ' ').split()[2:]]
        assert abs(nums[0] - 1e-4) < 1e-15 and abs(nums[1] + 2e-4) < 1e-15, l


# ---------------------------------------------------------------------------------------------
# MEASURE の初期分布の統計 (粒子数を増やす)
# ---------------------------------------------------------------------------------------------

def _mark_from_ring(name, ks, **kw):
    src = _src(name) + f'''
MARK ORGM = ();
QUAD SQ1 = (L=0.1 K1={ks} ROTATE=0.7853981633974483);
LINE RINGM = (ORGM SQ1 RING);
FFS;
USE RINGM; CELL; CALC;
Exit;
'''
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    r = ffs.ensure_rows()[1]            # ORGM
    mk = dict(AX=r.alfx, BX=r.betx, AY=r.alfy, BY=r.bety, EX=r.etx, EPX=r.etpx, EY=r.ety, EPY=r.etpy,
              R1=r.r1, R2=r.r2, R3=r.r3, R4=r.r4, DETR=r.detr, DX=r.x_cod, DPX=r.px_cod, DY=r.y_cod, DPY=r.py_cod)
    mk.update(kw)
    return mk


@pytest.mark.parametrize('np0,tol', [(20000, 0.03), (200000, 0.01)])
def test_ttinit_eigen_emittances_match_setting_for_a_coupled_large_ring_v2(v2_mode, np0, tol):
    """PF リングに skew 四極 (ks = 0.02) を入れた結合のある周期解の MARK (R1..R4 ≠ 0) から作った分布 (DP = 0) の 4x4 の固有エミッタンス (Σ·J の固有値) が
    設定した EMITX・EMITY と一致する (誤差 ∝ 1/√N: N = 2e4 で 3%、2e5 で 1% 以内)。第 1 ・第 2 モードの対応は EMITX・EMITY の大きさの順"""
    ex, ey = 1.5e-9, 4.0e-10
    mk = _mark_from_ring('PF', 0.02, EMITX=ex, EMITY=ey, DP=0.0, SIGZ=0.0)
    assert abs(mk['R1']) + abs(mk['R2']) + abs(mk['R3']) + abs(mk['R4']) > 0.05
    x = ms.ttinit(mk, np0, SadRandom(11), gauss=True, trgauss=True, jitter=True)
    c = np.cov(x[:, :4].T, bias=True)
    J = np.zeros((4, 4))
    J[0, 1] = J[2, 3] = 1.0
    J = J - J.T
    w = np.linalg.eigvals(c @ J)
    eig = sorted(abs(w.imag))[::2]
    want = sorted([ex, ey])
    for a, b in zip(sorted(eig), want):
        assert abs(a / b - 1) < tol, (eig, want)


def test_measure_beam_size_along_a_large_ring_matches_twiss_v2(v2_mode):
    """PF リング (skew なし) の周期解から作った分布 (EMITX = 1e-9、EMITY = 3e-10、DP = 1e-4 の正規分布、N = 20000) を 1 周の途中まで追った σx・σy が、
    Twiss の √(εx βx + (ηx δ)²)・√(εy βy + (ηy δ)²) と 4% 以内 (統計誤差 1/√(2N) = 0.5%、色収差・非線形の分)。MEASURE を 3 か所の要素で"""
    src = _src('PF') + '''
MARK ORGM = ();
LINE RINGM = (ORGM RING);
FFS;
$FORM="20.14"; PageWidth=1999;
USE RINGM; CELL; CALC;
EMITX=1e-9; EMITY=3e-10; DP=1e-4; SIGZ=1e-3;
NP=20000; GAUSS; TRGAUSS;
'''
    ffs, _ = ffscmd.run_sad_string(src + 'Exit;\n', out=io.StringIO())
    rows = ffs.ensure_rows()
    els = ffs.lat.elements
    sel = [i for i, e in enumerate(els) if type(e).__name__ == 'Quadrupole'][5::40][:3]
    def _pos(i):                       # 同じ名前の要素はファミリー: NAME.n (n 番目) で指定する
        n = sum(1 for e in els[:i + 1] if e.name == els[i].name)
        return f'{els[i].name}.{n}'
    prints = ''.join(f'r=FFS["MEASURE {_pos(i)}"]; Print["S {i} ", Sqrt[r[[3,1,1]]], " ", Sqrt[r[[3,3,3]]]];'+chr(10) for i in sel)
    buf = io.StringIO()
    ffscmd.run_sad_string(src + prints + 'Exit;\n', out=buf)
    got = {int(l.split()[1]): [float(q) for q in l.split()[2:]] for l in buf.getvalue().splitlines() if l.startswith('S ')}
    assert len(got) == 3
    for i, (sx, sy) in got.items():
        r = rows[i]            # MEASURE の位置は要素の入口 = Twiss の行 i (0 始まり) の要素の前
        wx = math.sqrt(1e-9 * r.betx + (r.etx * 1e-4) ** 2)
        wy = math.sqrt(3e-10 * r.bety + (r.ety * 1e-4) ** 2)
        assert abs(sx / wx - 1) < 0.04 and abs(sy / wy - 1) < 0.04, (i, sx, wx, sy, wy)
