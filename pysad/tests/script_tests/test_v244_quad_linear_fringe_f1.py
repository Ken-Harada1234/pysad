# -*- coding: utf-8 -*-
"""v2.4.4: 四極 (と K1≠0 の MULT) の線形フリンジの係数 F1 の符号 (docs/30 B323)。

端で K(s) が 0 → K0 に立ち上がるとき、K の 1 次の摂動から線形フリンジの係数は a = K0·I1/p (I1 = 硬い端との差の 1 次モーメント)。
ヘルプの定義 F1 = sgn(α)√|α|、α = 24(I0²/2 − I1) (I0 = 0 なら F1·|F1| = −24·I1) から a = −K0·F1·|F1|/(24 p)。F1 の符号が a の符号を決める。
K64 の CALC (tqlfre.f:12 の af1 = −ak/al·f1·|f1|/24) はこの形だが、追跡 (tffs.f tsetfringep:969 の −|akk·f1²|/24) は F1 の符号を捨てて常に a ≤ 0。
F1 < 0 (I1 > 0: 立ち上がりが硬い端より内側へ偏る・過剰に突き出る) で追跡と CALC が食い違う (K1 = 0.6・L = 0.3 m・FRINGE = 3・F1 = −0.1 で 4×4 が 2.1e-3)。
v2 は追跡の係数を CALC と同じ形 (符号を残す) に直した。k64 は K64 のまま (gs と同じ)。

独立な基準 (pysad のコードを使わない): 立ち上がり K(s) (tanh の端 + 内側のガウス状の盛り上がり) の線形運動 x'' = −K(s) x を solve_ivp で積分した 2×2 行列。
I0 = 0 になる s0 を求め、I1 > 0 (F1 < 0) を確かめる。
gs の保存した出力 (scripts/v244_quad_linear_fringe_f1__ring_*.gs.txt。FODO の 1 セル (L = 0.3 m、K1 = ±0.6、F1 = ±0.15、FRINGE = 3、DISFRIN = 1、3 GeV/c) の CALC の Twiss と
TrackParticles 1 周 (±1e-6 の 9 粒子) の保存) と比べる。
限界: 設計軌道まわりの線形部分だけ (軌道がずれたときの 2 次以上、δ の列、z の行は見ない)。放射の引数は K64 の癖 (長さの引数に係数 f1 を渡す) を従来のまま
−|f1| で渡すので、F1 < 0 の放射の結果は B323 の前と同じ (python 版と numba 版の一致と、引数の値だけ確かめる)。
"""
import io
import math
import os
import re

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from pysad import compat, lattice as lat, track6, track6n
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')


def _mode(m):
    old = compat._mode
    compat.set_mode(m)
    return old


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def make_elem(kind, L, K1, extra=''):
    """QUAD か MULT 1 台のラティス。K1 は積分強度 (K1/L が長さあたり)。戻り値: (ffs, 要素)"""
    el = ('QUAD' if kind == 'QUAD' else 'MULT') + ' SA=(L=%.15g K1=%.15g %s);' % (L, K1, extra)
    src = 'MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n' % el
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs, [e for e in ffs.lat.elements if e.name == 'SA'][0]


def track_jac4(ffs, h=1e-5):
    """追跡 (全要素 1 周、RFSW ON、放射なし) の設計軌道まわりの 4×4 ヤコビアン (中心差分)"""
    a = np.zeros((9, 6))
    for j in range(4):
        a[1 + 2 * j, j] = h
        a[2 + 2 * j, j] = -h
    o = np.asarray(track6.track_particles(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))
    M = np.zeros((4, 4))
    for j in range(4):
        M[:, j] = (o[1 + 2 * j, :4] - o[2 + 2 * j, :4]) / (2 * h)
    return M


def calc_mat4(el):
    M, _ = lat._el_cod_step(el, np.zeros(6), 0.0)
    return np.asarray(M)[:4, :4]


# ---------------------------------------------------------------- (i) v2: 追跡の 1 周写像の線形部分 = CALC の行列
CASES = [  # (extra, 説明)
    ('F1=0.1 FRINGE=3 DISFRIN=1', 'F1=+0.1 入口+出口'),
    ('F1=-0.1 FRINGE=3 DISFRIN=1', 'F1=-0.1 入口+出口'),
    ('F1=-0.25 FRINGE=3 DISFRIN=1', 'F1=-0.25 入口+出口'),
    ('F1=-0.1 FRINGE=1 DISFRIN=1', 'F1=-0.1 入口だけ'),
    ('F1=-0.1 FRINGE=2 DISFRIN=1', 'F1=-0.1 出口だけ'),
    ('F1=-0.1 FRINGE=0 DISFRIN=1', 'FRINGE=0 (F1 は効かない)'),
    ('F1=-0.1 FRINGE=3 DISFRIN=0', 'F1=-0.1 + 非線形フリンジ (DISFRIN=0)'),
    ('F1=0 F1K1F=-0.08 F1K1B=0.06 FRINGE=3 DISFRIN=1', 'F1K1F<0・F1K1B>0 (入口と出口で別の符号)'),
    ('F1=0 F1K1F=0.08 F1K1B=-0.06 FRINGE=3 DISFRIN=1', 'F1K1F>0・F1K1B<0'),
    ('F1=-0.1 F2=0.004 FRINGE=3 DISFRIN=1', 'F1=-0.1 + F2'),
    ('F1=0.1 F2=-0.004 FRINGE=3 DISFRIN=1', 'F1=+0.1 + F2<0'),
]


@pytest.mark.parametrize('kind', ['QUAD', 'MULT'])
@pytest.mark.parametrize('K1', [0.6, -0.6])
@pytest.mark.parametrize('extra,desc', CASES, ids=[c[1] for c in CASES])
def test_v2_tracking_matches_calc(kind, K1, extra, desc):
    """v2: F1 の符号・K1 の符号 (K1 < 0 は x・y の入れ替え)・入口と出口・FRINGE の各設定で、追跡のヤコビアンが CALC の行列と 1e-10 以下"""
    _mode('v2')
    ffs, el = make_elem(kind, 0.3, K1, extra)
    d = np.abs(track_jac4(ffs) - calc_mat4(el)).max()
    assert d < 1e-10, (kind, K1, desc, d)


def test_v2_f1_has_an_effect_with_its_sign():
    """F1 = −0.1 と +0.1 は別の行列を与える (符号が消えていない)。a = −K·F1·|F1|/24 の符号"""
    _mode('v2')
    ffs_p, el_p = make_elem('QUAD', 0.3, 0.6, 'F1=0.1 FRINGE=3 DISFRIN=1')
    ffs_m, el_m = make_elem('QUAD', 0.3, 0.6, 'F1=-0.1 FRINGE=3 DISFRIN=1')
    assert np.abs(track_jac4(ffs_p) - track_jac4(ffs_m)).max() > 1e-3
    # 係数 a = ∓K F1²/24 (K = K1/L = 2/m²: |a| = 8.3e-4)
    t_p = track6._tsetfringep(el_p, 0.6 / 0.3)
    t_m = track6._tsetfringep(el_m, 0.6 / 0.3)
    assert t_p[0] == pytest.approx(-2.0 * 0.01 / 24, rel=1e-14)
    assert t_m[0] == pytest.approx(+2.0 * 0.01 / 24, rel=1e-14)
    assert t_p[2] == t_p[0] and t_m[2] == t_m[0]
    # K の符号は係数に入らない (K1 < 0 は回転して |K| で通す)
    assert track6._tsetfringep(el_m, -2.0)[0] == t_m[0]


def test_v2_f1_nonnegative_coefficients_same_as_k64_bitwise():
    """F1 ≥ 0 (と F1K1F/B ≥ 0) の追跡の係数は、v2 が k64 とビット単位で同じ (copysign の形は F1 ≥ 0 で k64 の式と同じ値)。
    追跡の結果そのものは v2 が四極の積分を 4 次にしているので k64 と違う。係数以外を変えていないことは、修正前後の v2 の追跡の保存の比較で確認した"""
    for kind in ('QUAD', 'MULT'):
        for K1 in (0.6, -0.6):
            for ex in ('F1=0.1 F2=0.003 F1K1F=0.02', 'F1=0.37 F1K1B=0.01', 'F1=0 F2=-0.002'):
                ffs, el = make_elem(kind, 0.3, K1, ex + ' FRINGE=3 DISFRIN=0')
                coef = {}
                for m in ('k64', 'v2'):
                    _mode(m)
                    coef[m] = track6._tsetfringep(el, K1 / 0.3)
                assert coef['k64'] == coef['v2'], (kind, K1, ex)


# ---------------------------------------------------------------- (ii) k64: gs と同じ (F1 の符号を捨てる)
def test_k64_drops_the_sign_and_disagrees_with_calc():
    """k64: F1 = −0.1 と +0.1 の追跡は同じ。CALC とは F1 < 0 で 2e-3 級に食い違い、F1 > 0 では一致 (K1 = 0.6、L = 0.3 m、FRINGE = 3)"""
    _mode('k64')
    ffs_p, el_p = make_elem('QUAD', 0.3, 0.6, 'F1=0.1 FRINGE=3 DISFRIN=1')
    ffs_m, el_m = make_elem('QUAD', 0.3, 0.6, 'F1=-0.1 FRINGE=3 DISFRIN=1')
    Tp, Tm = track_jac4(ffs_p), track_jac4(ffs_m)
    assert np.abs(Tp - Tm).max() < 1e-9
    assert np.abs(Tp - calc_mat4(el_p)).max() < 1e-10
    d = np.abs(Tm - calc_mat4(el_m)).max()
    assert 1.5e-3 < d < 3e-3, d
    assert track6._tsetfringep(el_m, 2.0)[0] < 0.0          # 係数は常に ≤ 0


# ---------------------------------------------------------------- (iii) 独立な数値積分 (I1 > 0 の立ち上がり)
_D, _A, _C, _SG = 0.01, 0.3, 0.05, 0.02


def _g(s):
    """立ち上がり K(s)/K0: tanh の端 + 内側 (s = 0.05 m 付近) の過剰な盛り上がり (ガウス)"""
    return 0.5 * (1.0 + math.tanh(s / _D)) + _A * math.exp(-((s - _C) / _SG) ** 2)


def _moments():
    S, W = 0.3, 0.3

    def I(n, s0):
        f = lambda s: (s - s0) ** n * (_g(s) - (1.0 if s > s0 else 0.0))
        return quad(f, -S, W, points=[s0, 0.0, _C], limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    s0 = brentq(lambda s: I(0, s), -0.1, 0.2, xtol=1e-15)
    I0, I1, I2 = I(0, s0), I(1, s0), I(2, s0)
    al = 24.0 * (I0 ** 2 / 2.0 - I1)
    return s0, I0, I1, I2, math.copysign(math.sqrt(abs(al)), al), I2 - I0 ** 3 / 3.0


def _num_matrix(K, sign, s0, L=0.6, S=0.3):
    """x 平面の 2×2 (x, 傾き)。k(s) = sign·K·g(入口の端からの距離 + s0) (出口は鏡像)。両端 (硬い端の位置 0 と L) から S 外側まで積分して、外側のドリフトを戻す"""
    def rhs(s, v):
        sm = s if s < L / 2 else L - s
        k = sign * K * _g(sm + s0)
        M = v.reshape(2, 2)
        return np.array([[M[1, 0], M[1, 1]], [-k * M[0, 0], -k * M[0, 1]]]).ravel()
    sol = solve_ivp(rhs, (-S, L + S), np.eye(2).ravel(), method='DOP853', rtol=1e-13, atol=1e-17, max_step=0.005)
    M = sol.y[:, -1].reshape(2, 2)
    Di = np.array([[1.0, -S], [0.0, 1.0]])
    return Di @ M @ Di


def test_numeric_integration_reference():
    """I1 > 0 の立ち上がり: F1 < 0。数値積分の 2×2 と (i) 符号を残す式は K の 2 次の差 (K を 2 倍で 4 倍)、(ii) 符号を捨てる式は 2·K·I1 の差"""
    s0, I0, I1, I2, F1, F2 = _moments()
    assert abs(I0) < 1e-14 and I1 > 1e-4 and F1 < -0.05
    L = 0.6
    e_i, e_ii = [], []
    for K in (1e-3, 2e-3):
        Mn = _num_matrix(K, +1, s0)
        w = math.sqrt(K)
        body = np.array([[math.cos(w * L), math.sin(w * L) / w], [-w * math.sin(w * L), math.cos(w * L)]])
        fr = lambda a, b: np.array([[math.exp(a), b], [0.0, math.exp(-a)]])
        a_i, a_ii, b = -K * F1 * abs(F1) / 24.0, -abs(K * F1 * F1) / 24.0, K * F2
        assert a_i == pytest.approx(K * I1, rel=1e-12)
        e_i.append(np.abs(Mn - fr(-a_i, b) @ body @ fr(a_i, b)).max())
        e_ii.append(np.abs(Mn - fr(-a_ii, b) @ body @ fr(a_ii, b)).max())
    assert e_i[1] < 1e-9 and e_i[1] / e_i[0] == pytest.approx(4.0, rel=0.1)
    assert e_ii[0] > 1e-6 and e_ii[0] > 1000 * e_i[0]
    assert e_ii[1] / e_ii[0] == pytest.approx(2.0, rel=0.05)


@pytest.mark.parametrize('kind', ['QUAD', 'MULT'])
@pytest.mark.parametrize('sign', [+1, -1])
def test_pysad_vs_numeric_integration(kind, sign):
    """v2: 追跡と CALC の x ブロックが数値積分 (K の 1 次まで、差 K² ≲ 3e-10) と一致。k64 の追跡は 2.6e-6 食い違う (CALC は合う)。sign < 0 は x の発散"""
    s0, I0, I1, I2, F1, F2 = _moments()
    K = 2e-3
    Mn = _num_matrix(K, sign, s0)
    ex = 'F1=%.17g F2=%.17g FRINGE=3 DISFRIN=1' % (F1, F2)
    for m in ('v2', 'k64'):
        _mode(m)
        ffs, el = make_elem(kind, 0.6, sign * K * 0.6, ex)
        dc = np.abs(calc_mat4(el)[:2, :2] - Mn).max()
        dt = np.abs(track_jac4(ffs)[:2, :2] - Mn).max()
        assert dc < 5e-10, (m, dc)
        if m == 'v2':
            assert dt < 5e-10, dt
        else:
            assert dt > 2e-6, dt


# ---------------------------------------------------------------- gs の保存した出力との比較
def _read(name):
    return open(os.path.join(SCRIPTS, name), encoding='utf-8', errors='replace').read()


def _parse(t):
    d = {}
    m = re.search(r'^NXNY\s+([-0-9.]\S*)\s+(\S+)', t, re.M)
    d['nx'], d['ny'] = float(m.group(1)), float(m.group(2))
    cols = {}
    for k in 'XPYQ':
        m = re.search(r'^%s \{([^}]*)\}' % k, t, re.M)
        cols[k] = np.array([float(v) for v in m.group(1).replace('\n', ' ').split(',')])
    h = 1e-6
    M = np.array([[(cols[k][1 + 2 * j] - cols[k][2 + 2 * j]) / (2 * h) for j in range(4)] for k in 'XPYQ'])
    d['M'] = M
    d['tx'] = math.acos(0.5 * (M[0, 0] + M[1, 1]))
    d['ty'] = math.acos(0.5 * (M[2, 2] + M[3, 3]))
    return d


def _pysad_run(name, mode):
    _mode(mode)
    buf = io.StringIO()
    run_sad_string(_read('v244_quad_linear_fringe_f1__ring_%s.sad' % name), out=buf)
    return _parse(buf.getvalue())


# a: QF・QD とも F1 = +0.15、b: 両方 −0.15、c: QF だけ −0.15、d: QD だけ −0.15、e: F1 = 0
RINGS = {'a': 'F1=+0.15 (QF,QD)', 'b': 'F1=-0.15 (QF,QD)', 'c': 'F1=-0.15 (QF only)', 'd': 'F1=-0.15 (QD only)', 'e': 'F1=0'}


@pytest.mark.parametrize('name', list(RINGS), ids=list(RINGS.values()))
def test_k64_matches_gs_saved_output(name):
    """k64: CALC のチューン (Twiss NX・NY) と追跡の 1 周の 4×4 が gs の保存した出力と一致する (F1 の符号を捨てる動作の再現)"""
    g = _parse(_read('v244_quad_linear_fringe_f1__ring_%s.gs.txt' % name))
    k = _pysad_run(name, 'k64')
    assert abs(k['nx'] - g['nx']) < 1e-11 and abs(k['ny'] - g['ny']) < 1e-11
    assert np.abs(k['M'] - g['M']).max() < 1e-9


def test_gs_tune_calc_vs_tracking_disagree_only_for_negative_f1():
    """gs (保存した出力): 追跡の 1 周写像のチューン角が CALC の NX と、F1 > 0・F1 = 0 では 1e-9 以下で一致、F1 < 0 では 5e-3 以上食い違う"""
    for name in 'ae':
        g = _parse(_read('v244_quad_linear_fringe_f1__ring_%s.gs.txt' % name))
        assert abs(g['tx'] - g['nx']) < 1e-9 and abs(g['ty'] - g['ny']) < 1e-9
    for name in 'bcd':
        g = _parse(_read('v244_quad_linear_fringe_f1__ring_%s.gs.txt' % name))
        assert max(abs(g['tx'] - g['nx']), abs(g['ty'] - g['ny'])) > 5e-3
    # F1 = −0.15 の追跡は F1 = +0.15 と同じ (符号を捨てる)
    a, b = _parse(_read('v244_quad_linear_fringe_f1__ring_a.gs.txt')), _parse(_read('v244_quad_linear_fringe_f1__ring_b.gs.txt'))
    assert np.abs(a['M'] - b['M']).max() < 1e-9


@pytest.mark.parametrize('name', list(RINGS), ids=list(RINGS.values()))
def test_v2_tracking_tune_matches_calc_tune(name):
    """v2: 追跡の 1 周写像のチューン角が、F1 < 0 でも CALC の NX・NY (gs と同じ。F1 の符号を残す) と 1e-9 以下で一致する"""
    g = _parse(_read('v244_quad_linear_fringe_f1__ring_%s.gs.txt' % name))
    v = _pysad_run(name, 'v2')
    assert abs(v['nx'] - g['nx']) < 1e-11 and abs(v['ny'] - g['ny']) < 1e-11      # CALC は gs と同じ
    assert abs(v['tx'] - v['nx']) < 1e-9 and abs(v['ty'] - v['ny']) < 1e-9


# ---------------------------------------------------------------- python 版と numba 版の一致 (放射つきも)・放射の引数
@pytest.mark.parametrize('rad', [False, True])
@pytest.mark.parametrize('K1', [0.6, -0.6])
def test_v2_numba_equals_python_for_negative_f1(rad, K1):
    """v2・F1 < 0 (入口・出口で別の符号も): numba 版 (track6n) と python 版 (track6) の追跡が一致 (放射ありは FLUC なし)"""
    _mode('v2')
    ffs, el = make_elem('QUAD', 0.3, K1, 'F1=-0.1 F1K1B=0.25 F2=0.003 FRINGE=3 DISFRIN=0')
    rng = np.random.default_rng(5)
    P = rng.normal(0, 1, (6, 6)) * np.array([2e-3, 1e-3, 2e-3, 1e-3, 0.0, 2e-3])
    a = np.asarray(track6.track_particles(ffs.lat.elements, P.copy(), 3, ffs.momentum(), 0.0, False, rad, False))
    b = np.asarray(track6n.track_particles_fast(ffs.lat.elements, P.copy(), 3, ffs.momentum(), 0.0, False, rad, False))
    assert np.abs(a - b).max() < 1e-13


def test_radiation_argument_is_minus_abs_f1(monkeypatch):
    """放射の引数 (K64 の癖: tradk の長さの引数に係数 f1 を渡す) は B323 の前と同じ (−|f1|)。F1 = ±0.1 で同じ値を渡す"""
    _mode('v2')
    got = {}
    for f1 in (0.1, -0.1):
        args = []
        monkeypatch.setattr(track6, '_tradk', lambda p, al, *a, _args=args: _args.append(al))
        monkeypatch.setattr(track6, '_tsolqur_ens', lambda *a, **k: None)
        ffs, el = make_elem('QUAD', 0.3, 0.6, 'F1=%g FRINGE=3 DISFRIN=1' % f1)
        track6._tquad_ens([track6._Particle(np.array([1e-3, 0, 1e-3, 0, 0, 0.0]), track6._e._Phys(3e9, 100.0))], el, True)
        got[f1] = list(args)
    assert got[0.1] == got[-0.1] and len(got[0.1]) == 2
    assert all(v < 0.0 for v in got[0.1])
