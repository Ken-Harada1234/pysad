# -*- coding: utf-8 -*-
"""emit6d の解析ヤコビアン (st.piece(jac=...)) と tcod の直訳 (_tcod_radcod) の回帰テスト (docs/history/work_item6.md)。

(1) 各要素写像の解析ヤコビアンを数値微分 (刻み 1e-6) と比べる
(2) 要素ごとの piece の解析ヤコビアンを ESRF 1 セルのリング上で数値微分と比べる
(3) ESRF 1 セルのリング (bench_item6_tiny.sad) の Emittance 相当を SAD の結果 (bench_item6_tiny.result) と比べる
    RFSW OFF / RFSW / RFSW+RADCOD (RADCOD は tcod の z と trf0 の分け方まで SAD と一致する)
"""
import os

import numpy as np
import pytest

from pysad import lattice as _lat
from pysad import emit6d as E

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
LAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kek-ls_v03_12.lat')


def _fix(J):
    """piece() と同じ z 行の組み直し (列 x,px,y,py はシンプレクティック条件から)"""
    J = J.copy()
    for j in range(4):
        J[4, j] = (-J[0, j] * J[1, 5] + J[1, j] * J[0, 5] - J[2, j] * J[3, 5] + J[3, j] * J[2, 5])
    return J


def _check(f, jac, c, tol=1e-8, eps=1e-6):
    Jn = _fix(_lat._num_jacobian(f, c, eps))
    Ja = _fix(jac(c))
    d = np.abs(Ja - Jn) / np.maximum(np.abs(Jn), 1.0)
    assert d.max() < tol, (d.max(), np.unravel_index(np.argmax(d), d.shape))


PH = E._Phys(3e9, 603.12)


def _rand_cods(amp, n=4, seed=0):
    rng = np.random.default_rng(seed)
    for _ in range(n):
        c = rng.normal(size=6) * amp
        c[4] = 0.0
        c[5] = rng.normal() * 1e-3
        yield c


@pytest.mark.parametrize('amp', [0.0, 1e-5])
def test_drift_jacobian(amp):
    for c in _rand_cods(amp):
        f = E._drift_map(0.7, PH)
        _check(f, lambda c: E._jac_drift(c, 0.7, E._dvz(PH, c[5])), c)


@pytest.mark.parametrize('amp', [0.0, 1e-5])
def test_thin_jacobian(amp):
    for kord, ak in ((2, 0.4), (3, 12.0), (2, 0.0)):
        for L in (0.0, 0.05):
            f0 = lambda c, L=L, ak=ak, kord=kord: _lat._tthine_cod(c, L, ak, kord)

            def f(c, f0=f0, L=L):
                c2 = f0(c)
                c2[4] -= PH.dv(c2[5]) * L
                return c2
            for c in _rand_cods(amp):
                _check(f, lambda c, L=L, ak=ak, kord=kord: E._jac_thin(c, L, ak, kord, E._dvz(PH, c[5])), c)


@pytest.mark.parametrize('K1,L', [(1.03, 0.25), (-0.74, 0.2), (2.2, 0.5), (-1.7, 0.335)])
@pytest.mark.parametrize('amp', [0.0, 1e-5])
def test_quad_jacobians(K1, L, amp):
    ak1 = abs(K1) * L
    for c in _rand_cods(amp):
        _check(lambda c: _lat._tqfrie_cod(c, ak1, L), lambda c: E._jac_tqfrie(c, ak1, L), c)
        _check(lambda c: _lat._tqfrie_cod(c, -ak1, L), lambda c: E._jac_tqfrie(c, -ak1, L), c)
        _check(lambda c: _lat._tqente_cod(c, L / 3.0), lambda c: E._jac_tqente(c, L / 3.0), c)

        def fb(c):
            pr = 1.0 + c[5]
            c = _lat._tsolque_body_cod(c, L / 2.0, abs(K1) / pr)
            c[4] -= PH.dv(c[5]) * (L / 2.0)
            return c
        _check(fb, lambda c: E._jac_tsolque_body(c, L / 2.0, abs(K1), E._dvz(PH, c[5])), c)
        # 四極全体 (tqfrie + tsolque + tqfrie + dv 項)
        q = E._quad_map(L, K1, PH, 0, 0, 0, True, (0, 0, 0, 0, 0))
        _check(q, lambda c: E._jac_quad_full(c, L, K1, True, (0, 0, 0, 0, 0), PH), c)


def test_bend_jacobians():
    """K1=0 の BEND の入口・出口エッジと本体 (lattice の tbfrie/tbendebody0 の解析行列 + dv 項)"""
    rhob = 0.413 / 0.0199
    for amp in (0.0, 1e-5):
        for c in _rand_cods(amp):
            for ent in (True, False):
                _check(lambda c, ent=ent: _lat._tbedge_cod(c, rhob, 0.01, ent),
                       lambda c, ent=ent: _lat._bend_edge_chain(c, E._EYE6, rhob, 0.01, ent)[1], c, tol=1e-7, eps=1e-5)

            def f(c):
                c = _lat._tbendebody0_cod(c, 0.0199, rhob, 0.0)
                c[4] -= PH.dv(c[5]) * 0.413
                return c

            def j(c):
                J = _lat._tbendebody0_mat(c, 0.0199, rhob, 0.0)
                J[4, 5] += E._dvz(PH, c[5]) * 0.413
                return J
            _check(f, j, c, tol=1e-7, eps=1e-5)


def test_dv_derivative_is_h0_over_h1_cubed():
    for dp in (0.0, 1e-5, 1e-3, 1e-2):
        num = (PH.dv(dp + 1e-7) - PH.dv(dp - 1e-7)) / 2e-7
        assert abs(-num - E._dvz(PH, dp)) < 1e-10 * abs(num)


# ---------------------------------------------------------------------------
# ESRF 1 セルのリング: 解析ヤコビアンと数値微分の一致、SAD との照合
# ---------------------------------------------------------------------------

TINY = """
 CAVI RFT = (VOLT = 7e4 HARM = 50);
 LINE TINY = (MKS ARCA RFT MKE);
"""


@pytest.fixture(scope='module')
def tiny_lattice():
    if not os.path.exists(LAT):
        pytest.skip('pysad/tests/kek-ls_v03_12.lat がない')
    from pysad import parser
    with open(LAT, encoding='utf-8') as fh:
        src = fh.read()
    sc = parser.load_sad_string(src + TINY)
    return sc.build_lattice('TINY')


def test_piece_jacobians_on_cell(tiny_lattice):
    """リング 1 周の全 piece で、解析ヤコビアンが数値微分 (刻み 1e-5) と相対 1e-7 以内 (z 行は組み直し後)"""
    els = tiny_lattice.elements
    ph = E._Phys(3e9, sum(el.L for el in els))
    worst = [0.0]
    orig = E._State.piece

    def piece(self, f, symp_z=True, jac=None):
        if jac is not None and self.want_jac:
            Ja = jac(self.cod)
            if Ja is not None:
                Jn = _lat._num_jacobian(f, self.cod, 1e-5)
                if symp_z:
                    Ja, Jn = _fix(Ja), _fix(Jn)
                d = np.abs(Ja - Jn) / np.maximum(np.abs(Jn), 1.0)
                worst[0] = max(worst[0], d.max())
        return orig(self, f, symp_z, jac)
    E._State.piece = piece
    try:
        for c in (np.array([-1.8e-6, -2.4e-8, 0, 0, -1.1e-6, 9.8e-6]),
                  np.array([1.5e-4, 2e-5, -1e-4, 3e-5, 2e-4, 3e-4])):
            for rad in (False, True):
                E.tturne(els, c, ph, trf0=0.0206, rfsw=True, radcod=rad, want_db=False)
    finally:
        E._State.piece = orig
    assert worst[0] < 1e-7, worst[0]


def _read_sad_tiny():
    """bench_item6_tiny.result の K=1,2,3 (RFSW OFF / RFSW / RFSW+RADCOD) の値"""
    import re
    fn = os.path.join(ROOT, 'bench_item6_tiny.result')
    if not os.path.exists(fn):
        pytest.skip('bench_item6_tiny.result がない')
    d, k = {}, None
    for line in open(fn, encoding='utf-8', errors='replace'):
        m = re.match(r'^K\s+([\d.]+)', line)
        if m:
            k = int(float(m.group(1)))
            continue
        m = re.match(r'^(EMT|SIG|TUN|ALP|DMP|COD)\s*\{(.*)\}', line)
        if m and k:
            d[(k, m.group(1))] = [float(x) for x in m.group(2).split(',')]
    return d


@pytest.mark.parametrize('mode', [1, 2, 3])
def test_tiny_ring_vs_sad(tiny_lattice, mode):
    sad = _read_sad_tiny()
    rfsw = mode >= 2
    radcod = mode == 3
    e = E.emittance6d(tiny_lattice, 3e9, rfsw=rfsw, radcod=radcod)
    rel = lambda a, b: abs(a - b) / abs(b)
    # 横方向・チューン・U0 は全モードで 1e-9 以内 (εy は 0 同然なので見ない)
    assert rel(e['emittance_x'], sad[(mode, 'EMT')][0]) < 1e-9
    assert rel(e['tunes'][0], sad[(mode, 'TUN')][0]) < 1e-9
    assert rel(e['tunes'][1], sad[(mode, 'TUN')][1]) < 1e-9
    assert rel(e['energy_loss_eV'], sad[(mode, 'ALP')][1]) < 1e-9
    assert rel(e['momentum_compaction'], sad[(mode, 'ALP')][0]) < 1e-7
    for k in range(3):
        assert rel(e['damping_rate'][k], sad[(mode, 'DMP')][k]) < 1e-7
    # 縦方向: RFSW OFF でも SAD は空洞の積算量から Vc・h・バケット高さ・εz・σz・同期位置を作る (tecalc の omegaz 式)。
    # RFSW ON は νs=2e-4 と小さく、解析行列の 1e-13 の違いが 1e4 倍に拡大されるので 1e-6
    tz = 1e-6 if rfsw else 1e-9
    assert rel(e['emittance_z'], sad[(mode, 'EMT')][2]) < tz
    assert rel(e['momentum_spread'], sad[(mode, 'SIG')][0]) < tz
    assert rel(e['bunch_length'], sad[(mode, 'SIG')][1]) < tz
    assert rel(e['tunes'][2], sad[(mode, 'TUN')][2]) < tz if rfsw else True
    assert rel(e['bucket_height'], sad[(mode, 'ALP')][4]) < tz
    assert rel(e['rf_voltage_eff'], sad[(mode, 'ALP')][2]) < 1e-12
    # 同期位置 (Equilibrium dz): RFSW OFF/ON 非 RADCOD は全桁、RADCOD は tcod の直訳で SAD と同じ z/trf0 の分け方
    assert rel(e['trf0'], sad[(mode, 'ALP')][3]) < 1e-9
    if radcod:
        cod = e['closed_orbit']
        sc = sad[(mode, 'COD')]
        # SAD の OrbitAtExit は出口の cod、pysad の closed_orbit は入口。tcod の収束判定が緩いので両者は 1e-8 程度ずれる。
        # z (RF 位相の分け方: trf0 と z) が SAD と同じ分け方になっていること (旧実装は z≈0, trf0≈0.0205 だった) を見る
        assert abs(cod[4] - sc[4]) < 1e-6 * abs(sc[4])
        assert abs(cod[0] - sc[0]) < 1e-8
        assert abs(cod[5] - sc[5]) < 1e-9


def test_numeric_and_analytic_agree_on_tiny_ring(tiny_lattice):
    """解析ヤコビアンを切って (数値微分に戻して) も、非 RADCOD の結果は 1e-6 以内で一致する"""
    e1 = E.emittance6d(tiny_lattice, 3e9, rfsw=True)
    E._USE_ANALYTIC_JAC = False
    try:
        e0 = E.emittance6d(tiny_lattice, 3e9, rfsw=True)
    finally:
        E._USE_ANALYTIC_JAC = True
    for k in ('emittance_x', 'momentum_spread', 'bunch_length', 'momentum_compaction'):
        assert abs(e1[k] - e0[k]) < 1e-6 * abs(e0[k]), k


def test_beam_tilt_sizes():
    B = np.zeros((6, 6))
    B[0, 0], B[2, 2], B[0, 2] = 4e-8, 1e-8, 1e-8      # σxx, σyy, σxy
    tilt, sx, sy = E.beam_tilt_sizes(B)
    # 回転角: atan2(-2 σxy, σxx-σyy)/2 (temitprint)
    assert abs(tilt - 0.5 * np.arctan2(-2e-8, 3e-8)) < 1e-15
    assert sx > sy > 0.0
    # 結合なし → tilt = 0
    B[0, 2] = 0.0
    assert E.beam_tilt_sizes(B)[0] == 0.0


def test_squeeze_space_matches_tfsqueezespace():
    """tfdisp.f tfsqueezespace: 直前も空白の空白を先頭から n 個まで除き、末尾を空白で埋める"""
    from pysad.ffscmd import _squeeze_space
    assert _squeeze_space('ab  cd   ef', 2) == ('ab cd  ef  ', 2)
    assert _squeeze_space('ab  cd   ef', 0) == ('ab  cd   ef', 0)
    assert _squeeze_space('a b c', 3) == ('a b c', 0)       # 連続空白がなければ何も除かない


def test_disp_reversed_names_and_columns():
    """DISP: 反転配置の要素名に `-`、桁位置は SAD (bench_item6_dispnames.result) と同じ、ヘッダは 66 行ごと、`disp B*` は B* の行だけ"""
    import contextlib
    import io
    fn = os.path.join(ROOT, 'bench_item6_dispnames.sad')
    res = os.path.join(ROOT, 'bench_item6_dispnames.result')
    if not os.path.exists(fn) or not os.path.exists(res):
        pytest.skip('bench_item6_dispnames がない')
    from pysad.ffscmd import run_sad
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        run_sad([fn])
    txt = out.getvalue()
    sad = open(res, encoding='utf-8', errors='replace').read()

    def names_cols(t):
        rows = {}
        for line in t.splitlines():
            tok = line.split()
            if len(tok) >= 14 and tok[-1].isdigit() and tok[5][0] in '-ABCDEFGHIJKLMNOPQRSTUVWXYZ$':
                rows[int(tok[-1])] = line
        return rows
    a, b = names_cols(sad), names_cols(txt)
    assert len(a) == len(b) == 35
    # 名前から s(m) の桁までを比べる (雑音の EY 等は比べない)。squeeze の個数が雑音の桁に依るので、10 桁を超える名前の反転 2 行は除く
    for k in a:
        if len(a[k].split()[5]) > 14:
            continue
        assert a[k][40:len(a[k]) - 52] == b[k][40:len(b[k]) - 52], (k, a[k], b[k])
    assert sum(1 for l in txt.splitlines() if l.lstrip().startswith('AX')) == 1
    assert any(l.split()[5].startswith('-') for l in b.values())
