# -*- coding: utf-8 -*-
"""D2 (4): CorrectOrbit の Micado (`Solver->{"Micado",n}`) の K64 の誤り。

K64 の pmicad.f は、選んだ nx 本の補正量をハウスホルダー変換したあと、新しい tsvd (配列の大きさを size() で取る) に替えたときに nx を渡し忘れ、
全 n 本の最小二乗解を解いて先頭の nx 個だけを使う (選ばなかった補正量の寄与が落ちる)。24 セルの FODO (BPM 48 台、垂直ステアラ 24 台) で、
Micado の補正後の垂直の rms が補正前 (9.41e-5 m) より悪化する (n=1: 1.46e-4)。v2 は上三角の nx x nx だけを解く (旧 `tsolva(a,b,x,nx,nx,mdim,eps)`)。
参照値は gs (K64、2026-10-19 ビルド)。
"""
import io
import os
import re

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, '..', '..', 'example', '13_packages', 'fn_CorrectOrbit_Monitor_Steer.sad')
PKG = os.environ.get('PYSAD_SAD_PACKAGES')


def _problem(seed=3, m=40, n=16):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(m, n)) * rng.uniform(0.3, 3.0, size=n)
    b = a @ rng.normal(size=n) * 0.5 + rng.normal(size=m) * 0.3
    return a, b


def _selected_ls(a, b, x):
    sel = np.nonzero(x)[0]
    xs = np.linalg.lstsq(a[:, sel], b, rcond=None)[0]
    return sel, xs


def test_micado_v2_is_least_squares_on_selected_correctors(v2_mode):
    from pysad.corr import linear_solve_conditioned
    a, b = _problem()
    prev = np.inf
    for nx in (1, 2, 3, 5, 8, 12, 16):
        x = linear_solve_conditioned(a, b, 1e-8, nx, True)
        sel, xs = _selected_ls(a, b, x)
        assert len(sel) == nx
        assert np.allclose(x[sel], xs, rtol=1e-9, atol=1e-12)          # 選んだ補正量の最小二乗解
        r = np.linalg.norm(b - a @ x)
        assert r <= prev * (1 + 1e-12)                                    # 本数を増やすと残差は減る
        prev = r
    assert np.allclose(x, np.linalg.lstsq(a, b, rcond=None)[0], rtol=1e-8, atol=1e-12)     # 全部使えば最小二乗


def test_micado_k64_keeps_truncated_full_least_squares(k64_mode):
    from pysad.corr import linear_solve_conditioned
    a, b = _problem()
    xf = np.linalg.lstsq(a, b, rcond=None)[0]
    x = linear_solve_conditioned(a, b, 1e-8, 3, True)
    sel = np.nonzero(x)[0]
    assert len(sel) == 3
    # K64: 全 16 本の最小二乗解の先頭 3 個 (ハウスホルダー変換後の基底での値) を選んだ補正量に割り当てる。選んだ補正量だけの最小二乗解とは違う
    _, xs = _selected_ls(a, b, x)
    assert not np.allclose(x[sel], xs, rtol=1e-3)
    assert np.linalg.norm(b - a @ x) > np.linalg.norm(b - a[:, sel] @ xs)


def _ring_script(n):
    src = open(EXAMPLE, encoding='utf-8').read()
    head = src[:src.index('\nFFS;')]
    return head + '''
FFS; PageWidth=1999;
USE RING; CELL; CALC;
mon=Monitor[Flatten[Table[{"BA"//ToString[k],"BB"//ToString[k]},{k,1,24}]]];
sv=Steer[Table["SV"//ToString[k],{k,1,24}]];
rmy:=Sqrt[Plus@@(Twiss["DY",mon[[1]]]^2)/Length[mon[[1]]]];
Print["R0 ",rmy];
opt=CalculateOptics[1,LINE["LENGTH"],Twiss["*",1],FLAG["CELL"],2];
r=%s;
CALC;
Print["R1 ",rmy];
Exit;
''' % ('CorrectOrbit["Y",mon,sv,opt]' if n == 0 else 'CorrectOrbit["Y",mon,sv,opt,Solver->{"Micado",%d}]' % n)


def _rms(n):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_ring_script(n), out=buf)
    out = buf.getvalue()
    v = {k: float(re.search(r'^%s ([-+0-9.eE]+)' % k, out, re.M).group(1)) for k in ('R0', 'R1')}
    return v['R0'], v['R1']


@pytest.mark.skipif(not PKG, reason='PYSAD_SAD_PACKAGES が無い')
def test_micado_24cell_ring_k64_matches_gs_and_gets_worse(k64_mode):
    # gs の補正後の rms (rmy は Sqrt[Σ DY^2 / 48]): n=0 (最小二乗) 2.0567e-5、Micado n=1 1.46285e-4、n=3 8.17755e-5、n=24 2.0567e-5。補正前 9.40603e-5
    r0, r1 = _rms(0)
    assert r0 == pytest.approx(9.406029839523e-5, rel=1e-9)
    assert r1 == pytest.approx(2.0566604e-5, rel=1e-5)
    assert _rms(1)[1] == pytest.approx(1.4628509e-4, rel=1e-4)
    assert _rms(3)[1] == pytest.approx(8.1775516e-5, rel=1e-4)
    assert _rms(1)[1] > r0                                   # K64 は補正前より悪化する


@pytest.mark.skipif(not PKG, reason='PYSAD_SAD_PACKAGES が無い')
def test_micado_24cell_ring_v2_decreases_monotonically(v2_mode):
    r0, _ = _rms(0)
    prev = r0
    for n in (1, 3, 6, 12, 24):
        _, r1 = _rms(n)
        assert r1 < prev * 1.0000001
        prev = r1
    assert _rms(1)[1] == pytest.approx(4.7399638e-5, rel=1e-4)        # pysad v2 (独立な貪欲法 + 最小二乗の線形モデルと整合: 4.1e-5 前後へ単調に下がる)
    assert _rms(6)[1] == pytest.approx(4.1117702e-5, rel=1e-4)
    assert _rms(24)[1] == pytest.approx(2.0566605e-5, rel=1e-5)


# ---------------------------------------------------------------- (5) CSROYMatrix (B60)
GS_OYMAT = [0.014472081731347494, -2.769057571207063e-07, 3.5612223931331085e-08, -1.4443896972909878e-07, 5.592772512398629e-08,
            -1.6286863687502863e-07, 1.861146310621842e-07, 0.05788834553896768, 1.0363496011494454e-07, 1.2162938318034985e-08,
            8.529427262506952e-08, -1.604931461594072e-08, -4.624243805276703e-08, -1.5244085983423994e-07, 0.014544362451563281,
            -9.314617692947756e-08, 1.695913851393607e-08, -9.108689368985506e-08, 1.3883194302495384e-07, -5.384864747518984e-08,
            7.721880506787548e-08, 0.05817744893760843, 6.386213879865394e-08, -2.5414344834357217e-08, -7.021513567630911e-08,
            -1.1619135840189766e-07, -1.8363655079360367e-08, -7.81951675158796e-08, 0.014616809377895017, -7.1913054564385e-08,
            1.300817498347352e-07, -1.4155669108334983e-08, 7.116609325715078e-08, 1.756417018040006e-08, 5.5542605522326526e-08,
            0.05846723466016907]
OY_SCRIPT = '''FFS; PageWidth=1999;
nr=64; nj=3; nphi=2; mphi=nphi*16;
wi=Table[0.01*Sin[0.37 i]/(1+i/10.),{i,1,nr}];
aj=Table[1e-8*i,{i,1,nj}];
om=Table[0.12+0.0003*i,{i,1,nj}];
zj=Table[2e-5*Sin[0.2 k+0.1]*(1+0.3 Mod[k,5]),{k,1,nj*mphi}];
dpj=Table[3e-4*Cos[0.17 k]*(1-0.1 Mod[k,7]),{k,1,nj*mphi}];
rj=Table[0.22/i,{i,1,nj}];
las={wi,aj,om,zj,dpj,rj};
m=CSROYMatrix[las,{2e-5,1e-8,0.01,3e-4,1e-3}];
Print["MD ",Dimensions[m]];
Print["MR ",Flatten[m]];
Exit;
'''


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_csroymatrix_hand_made_input_matches_gs(mode):
    """B60 の原因 (D2): gs の `CSROYMatrix` は、手で作った正しい形の実数リストなら動き、pysad と相対 1e-12 以内で一致する (gs: 36 要素、MD={6,6})。
    wrongtype になるのは `CSRSetupOY` の出力 (K64 は 4 番目の zj のリストの頭 (head) に zjn(1,1) を上書きし、zj が実数リストでなくなる)。"""
    from pysad import compat
    from pysad.ffscmd import run_sad_string
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        run_sad_string(OY_SCRIPT, out=buf)
    finally:
        compat.set_mode(old)
    out = buf.getvalue()
    assert re.search(r'^MD \{6,6\}', out, re.M)
    row = re.search(r'^MR \{(.*?)\}\s*$', out, re.M | re.S).group(1)
    v = np.array([float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', row)])
    assert v.shape == (36,)
    assert np.allclose(v, GS_OYMAT, rtol=1e-12, atol=1e-18)


# ---------------------------------------------------------------- (5) tffamsetup (fam 点): 独立な解析の基準
# gs は FITP>1 の GO で再起動するので fam 点の値は gs と比べられない (docs/history/v2next_C.md)。代わりに 2 つの解析式で確かめる:
#  (a) 初期条件の Courant-Snyder 不変量 = (xamp)^2 * em (x と y の面ごと)。tffamsetup.f の dx(1)=c x0、dx(2)=s px0 - (α/β) dx(1) の恒等式
#  (b) 厚さ 0 の八極 (K3L) が MARK の直後にある線形ラティスで、fam 点のチューン = 八極のキックのヤコビアンを一周行列に掛けた式
#      cos μ = cos μ0 + (J β / 2) sin μ0 (J = ∂Δpx/∂x = -K3L (x² - y²)/2、y 面は -J)
_FAM_SRC = '''
MOMENTUM = 1GEV;
DRIFT D = (L=0.7);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.65);
MARK MK = ();
MULT OC = (L=0 K3=%s);
LINE CELL = (QF D QD D);
LINE RING = (MK OC 8*CELL);
FFS;
USE RING;
CELL;
CALC;
EMITX = 1E-6; EMITY = 1E-7;
DP = 0.001;
fit $$$ NX 2.21 3 NY 1.31 3;
free QF K1 QD K1;
MatchingAmplitude = {{0, 3, 2}};
NFAMP = 4;
calc;
'''


def _fam_ffs(k3):
    from pysad.ffscmd import run_sad_string
    return run_sad_string(_FAM_SRC % k3, out=io.StringIO())[0]


def test_fam_initial_conditions_have_courant_snyder_invariant_of_emittance_times_amplitude():
    ffs = _fam_ffs(0.0)
    m = ffs.matcher_()
    m._fam_define(1)
    r = ffs.row(1)
    em = 1.1e-6                         # |EMITX| + |EMITY|
    n = 0
    for ix, f in m._fam.items():
        # 正規座標の物理座標への変換 (結合なし: R1..R4 = 0) なので dfam はそのまま (x, px, y, py)
        x, px, y, py = f['dfam']
        if f['kfam'] > 0:
            amp2 = 3.0 ** 2
            inv = (x * x + (r.alfx * x + r.betx * px) ** 2) / r.betx
            assert abs(inv - amp2 * em) < 1e-12 * em and abs(y) + abs(py) < 1e-15
        else:
            amp2 = 2.0 ** 2
            inv = (y * y + (r.alfy * y + r.bety * py) ** 2) / r.bety
            assert abs(inv - amp2 * em) < 1e-12 * em and abs(x) + abs(px) < 1e-15
        n += 1
    assert n == 8                       # NFAMP=4 本 x 2 面


def test_fam_point_tunes_match_thin_octupole_kick_jacobian():
    import math
    from pysad import lattice as L
    K3L = 6000.0
    ffs0 = _fam_ffs(0.0)
    rows = ffs0.lat.twiss_table(delta=0.0)
    r0 = ffs0.row(1)
    mu0 = (rows[-1].mux, rows[-1].muy)
    beta = (r0.betx, r0.bety)
    ffs = _fam_ffs(K3L)
    m = ffs.matcher_()
    m._fam_define(1)
    r1 = ffs.row(1)
    L.RFSW = False
    checked = 0
    for ix, f in m._fam.items():
        c0 = np.array([r1.x_cod, r1.px_cod, r1.y_cod, r1.py_cod, 0.0, 0.0])
        c0[:4] += f['dfam']
        rk = m._calc_optics_fixed_default(0.0, c0)
        x, y = c0[0], c0[2]
        jx = -K3L * (x * x - y * y) / 2.0
        for mu_fam, mu_0, b, j in ((rk[-1].mux, mu0[0], beta[0], jx), (rk[-1].muy, mu0[1], beta[1], -jx)):
            c = math.cos(mu_0) + j * b * math.sin(mu_0) / 2.0
            mu = math.acos(c)
            if math.sin(mu_0) < 0:
                mu = 2 * math.pi - mu
            # 線形の腕と八極 1 個の近似なので (二つ目の厚い四極の非線形項で) 1e-4 程度の差。八極のキックが効く点 (x または y が 3 mm 以上) では差の 1% 以内
            assert abs(mu_fam - mu) / (2 * math.pi) < 3e-4
        if abs(x) > 3e-3 or abs(y) > 3e-3:
            assert abs(rk[-1].mux - mu0[0]) / (2 * math.pi) > 1e-2 or abs(rk[-1].muy - mu0[1]) / (2 * math.pi) > 1e-2
            checked += 1
    assert checked == 4


# ---------------------------------------------------------------- (6) numba のキャッシュの依存の追跡 (pysad/nbcache.py)
def test_nbcache_purges_caller_when_callee_file_changes(tmp_path, monkeypatch):
    from pysad import nbcache
    pkg = tmp_path / 'pkg'
    cache = pkg / '__pycache__'
    cache.mkdir(parents=True)
    (pkg / 'consts.py').write_text('ELMASS = 1.0\nOTHER = 2.0\n', encoding='utf-8')
    (pkg / 'modb.py').write_text('from numba import njit\n@njit(cache=True)\ndef g(x):\n    return x + 1.0\n', encoding='utf-8')
    (pkg / 'moda.py').write_text('from numba import njit\nfrom .modb import g\nfrom .consts import ELMASS\n@njit(cache=True)\ndef f(x):\n    return g(x) * ELMASS\n',
                                 encoding='utf-8')
    for name in ('moda.f-3.py312.nbi', 'moda.f-3.py312.1.nbc', 'modb.g-2.py312.nbi', 'modb.g-2.py312.1.nbc'):
        (cache / name).write_bytes(b'x')
    monkeypatch.setattr(nbcache, '_PKG', str(pkg))
    monkeypatch.setattr(nbcache, '_CACHE', str(cache))
    monkeypatch.setattr(nbcache, '_RECORD', str(cache / 'nbdeps.json'))
    assert nbcache._closure('moda') == {'modb', 'consts:ELMASS'}
    assert nbcache.purge() == []                                   # 記録がない最初の 1 回は記録だけ
    assert nbcache.purge() == []                                   # 変化なし
    (pkg / 'consts.py').write_text('ELMASS = 1.0\nOTHER = 3.0\n', encoding='utf-8')     # 使っていない定数の変更では捨てない
    assert nbcache.purge() == []
    (pkg / 'modb.py').write_text('from numba import njit\n@njit(cache=True)\ndef g(x):\n    return x + 100.0\r\n', encoding='utf-8')
    gone = nbcache.purge()
    assert gone == ['moda.f-3.py312.1.nbc', 'moda.f-3.py312.nbi']  # 呼ぶ側 (moda) だけ捨てる。modb は自分のファイルのハッシュで numba が見る
    assert (cache / 'modb.g-2.py312.nbi').exists()
    (pkg / 'consts.py').write_text('ELMASS = 2.0\nOTHER = 3.0\n', encoding='utf-8')     # 使っている定数
    for name in ('moda.f-3.py312.nbi', 'moda.f-3.py312.1.nbc'):
        (cache / name).write_bytes(b'x')
    assert nbcache.purge() == ['moda.f-3.py312.1.nbc', 'moda.f-3.py312.nbi']
