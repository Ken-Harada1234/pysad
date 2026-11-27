# -*- coding: utf-8 -*-
"""B137 (v2): BEND の端の縦磁場 ψ の歳差を、刻みの本体の歳差と別の回転として掛ける (compat.spin_edge_split)。

K64 の tradk の sprot (追跡)・tradke (Emittance[]) は入口の ψ を最初の刻み、出口の ψ を最後の刻みの回転ベクトルに足して 1 回で回す。回転は非可換なので
刻み 1 つ分の偏向角に比例する誤差が残り、BEND のスピン軌道結合 D の y の列が刻みの数に反比例してずれる (D の x 成分の誤差が −D/n)。
v2 は入口の ψ を最初の刻みの本体の前に、出口の ψ を最後の刻みの本体の後に掛ける。

独立な基準: Maxwell 方程式を満たす硬い端の扇形 BEND (15°、L = 1.2 m、0.5 GeV、K1 = 0、E1 = E2 = 0) の中の Thomas-BMT の数値積分
(`_spin_maxwell_indep.maxwell_R`、y0 = ±2e-4・py = ±2e-4 の差。出口の端のフリンジ幅 w = 1e-3・5e-4 m の値を w → 0 に外挿 (w に比例))。
回転ベクトル ω = rotvec(R(u) R(0)ᵀ) の (y, py) での微分 (出口の座標系)。基準値 (w → 0):
    ∂ω/∂y  = (+0.0639335, 0, +0.0095666)、 ∂ω/∂py = (−0.0437546, 0, +0.554514)。
(w = 4e-3 m の基準は y の列が 1% 小さく、py の列の x 成分が 0.9% 小さい (w に比例して収束する)。)"""
import io
import math
import os

import numpy as np
import pytest
from scipy.spatial.transform import Rotation as Rot

from pysad import compat, emit6d, ffscmd, track6

E, L, TH = 0.5, 1.2, 0.2617993877991494
H = 2e-4
REF_Y = np.array([0.0639335, 0.0, 0.0095666])
REF_PY = np.array([-0.0437546, 0.0, 0.554514])


def _els(nslice=1):
    names = [f'B{i}' for i in range(nslice)]
    src = f'MOMENTUM = {E} GEV;\nDRIFT D1=(L=1);\n'
    for n in names:
        src += f'BEND {n} = (L={L / nslice} ANGLE={TH / nslice} E1=0 E2=0);\n'
    src += f'MARK ST = ();\nLINE RING = (ST {" ".join(names)} D1);\nFFS;\nUSE RING;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    return [e for e in ffs.lat.elements if e.name in names], ffs.momentum()


def _R_track(els, p, y, py):
    cols = []
    for j in range(3):
        out, s1 = track6.track_particles(els, np.array([[0, 0, y, py, 0, 0.0]]), 1, p, 0.0, False, True, False, spins=[np.eye(3)[j]], radpol=False, circ=100.0)
        cols.append(s1[0])
    return np.array(cols).T, out[0]


def _R_emit(els, p, y, py, ndiv):
    emit6d.POL_MODE = True
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [ndiv, 1, 1]
    try:
        st = emit6d.tturne(els, np.array([0, 0, y, py, 0, 0.0]), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
        return st.srot[:, :3].copy()
    finally:
        emit6d.EMITDIV[:] = old
        emit6d.POL_MODE = False


def _D(Rf):
    """回転ベクトルの (y, py) 微分 (中心差分)"""
    R0 = Rf(0.0, 0.0)
    out = []
    for dy, dpy in ((H, 0.0), (0.0, H)):
        wp = Rot.from_matrix(Rf(dy, dpy) @ R0.T).as_rotvec()
        wm = Rot.from_matrix(Rf(-dy, -dpy) @ R0.T).as_rotvec()
        out.append((wp - wm) / (2 * H))
    return out


@pytest.fixture
def v2_split():
    compat.set_mode('v2')
    compat.set_spin_edge_split(None)
    yield
    compat.set_spin_edge_split(None)
    compat.set_mode(None)


def test_tracking_D_y_column_exact_for_any_step_count(v2_split, monkeypatch):
    """入口・出口の ψ を別の回転にすると、D の y の列は刻みの数に依らず基準 (w → 0) と 3e-5 で合う (刻み 1・6・既定の 53)。
    分けない (K64 のまま足す) と刻み 6 で −1.6e-2 (D の 1/6)"""
    els, p = _els()
    orig = track6._ndivrad
    for nd in (1, 6, 0):
        monkeypatch.setattr(track6, '_ndivrad', orig if nd == 0 else (lambda n: (lambda *a, **k: n))(nd))
        compat.set_spin_edge_split(True)
        dy, dpy = _D(lambda y, py: _R_track(els, p, y, py)[0])
        assert np.abs(dy - REF_Y).max() < 3e-5, (nd, dy)
        compat.set_spin_edge_split(False)
        dy0, _ = _D(lambda y, py: _R_track(els, p, y, py)[0])
        if nd == 6:
            assert abs(dy0[0] - REF_Y[0]) > 8e-3, dy0


def test_tracking_D_py_column_converges_with_the_step(v2_split, monkeypatch):
    """py の列は 1 刻みの偏向角の 1 次の誤差が残る (B235 の続き): 既定 (53 刻み = 偏向角 5e-3 rad) で x 成分 8e-5・z 成分 2e-5。刻みに比例して減る (刻み 6 の 8 倍近い)"""
    els, p = _els()
    orig = track6._ndivrad
    err = {}
    for nd in (6, 0):
        monkeypatch.setattr(track6, '_ndivrad', orig if nd == 0 else (lambda n: (lambda *a, **k: n))(nd))
        _, dpy = _D(lambda y, py: _R_track(els, p, y, py)[0])
        err[nd] = dpy - REF_PY
    assert abs(err[0][0]) < 2e-4 and abs(err[0][2]) < 1e-4, err[0]
    assert abs(err[6][0]) > 5 * abs(err[0][0])


def test_sliced_bend_equals_unsliced_with_edge_split(v2_split):
    """BEND を 2・8・24 本に分けても (各断片の入口・出口の ψ は隣と打ち消す)、y = 3 mm・py = 1e-3 のスピン回転 R は分けない BEND と 5e-8 以内。
    分けない足し方 (K64) は断片の数に応じて 4e-6〜7e-5 ずれる。独立な DK の積分 (`_p1_dk.py`、BEND・QUAD を 1 刻み ≤ 0.05 m に切る) の刻み依存はこれだった"""
    e1, p = _els(1)
    compat.set_spin_edge_split(True)
    r1, _ = _R_track(e1, p, 3e-3, 1e-3)
    d_on, d_off = [], []
    for ns in (2, 8, 24):
        en, p = _els(ns)
        compat.set_spin_edge_split(True)
        d_on.append(np.abs(_R_track(en, p, 3e-3, 1e-3)[0] - r1).max())
        compat.set_spin_edge_split(False)
        d_off.append(np.abs(_R_track(en, p, 3e-3, 1e-3)[0] - r1).max())
    assert max(d_on) < 5e-8, d_on
    assert d_off[0] > 2e-6 and d_off[2] > 3e-5, d_off


def test_emittance_D_y_column_exact_at_default_division(v2_split):
    """Emittance[] の R(y, py) の y の列 (EMITDIV = 16、POL ON の既定): 基準と 5e-5 (分けないと 1.5e-3)。py の列は刻みの 1 次の誤差が残る (EMITDIV = 16 で z 成分 6e-3)"""
    els, p = _els()
    compat.set_spin_edge_split(True)
    dy, dpy = _D(lambda y, py: _R_emit(els, p, y, py, 16))
    assert np.abs(dy - REF_Y).max() < 5e-5, dy
    assert abs(dpy[0] - REF_PY[0]) < 2e-3 and abs(dpy[2] - REF_PY[2]) < 8e-3, dpy
    compat.set_spin_edge_split(False)
    dy0, _ = _D(lambda y, py: _R_emit(els, p, y, py, 16))
    assert abs(dy0[0] - REF_Y[0]) > 1e-3, dy0


def test_emittance_D_columns_of_srot_match_the_reference(v2_split):
    """srot の D の列 (4〜9 列。符号は回転ベクトルと逆): y の列は基準と 5e-5 (EMITDIV = 16)"""
    els, p = _els()
    emit6d.POL_MODE = True
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [16, 1, 1]
    try:
        st = emit6d.tturne(els, np.zeros(6), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
    finally:
        emit6d.EMITDIV[:] = old
        emit6d.POL_MODE = False
    assert np.abs(-st.srot[:, 5] - REF_Y).max() < 5e-5, st.srot[:, 5]


def test_k64_is_unchanged(monkeypatch):
    """k64: 端の ψ は K64 のまま 1 回の回転に足す (spin_edge_split() は常に False)。D の y の列は刻み 6 で基準の 5/6 (x 成分 −1.6e-2 ずれる)"""
    els, p = _els()
    compat.set_mode('k64')
    try:
        assert compat.spin_edge_split() is False
        monkeypatch.setattr(track6, '_ndivrad', lambda *a, **k: 6)
        dy, _ = _D(lambda y, py: _R_track(els, p, y, py)[0])
    finally:
        compat.set_mode(None)
    assert abs(dy[0] - REF_Y[0] * 5.0 / 6.0) < 3e-3, dy


@pytest.mark.skipif(not os.environ.get('PYSAD_SLOW'), reason='PYSAD_SLOW=1 (独立な積分、約 1 分)')
def test_reference_regeneration():
    from . import _spin_maxwell_indep as M
    Ep = E * 1e9

    def D(w):
        def R(y, py):
            return M.maxwell_R(Ep, L, TH, py, w=w, zs=6 * w, y0=y)
        R0 = R(0.0, 0.0)
        out = []
        for dy, dpy in ((H, 0.0), (0.0, H)):
            out.append((Rot.from_matrix(R(dy, dpy) @ R0.T).as_rotvec() - Rot.from_matrix(R(-dy, -dpy) @ R0.T).as_rotvec()) / (2 * H))
        return out
    a, b = D(1e-3), D(5e-4)
    dy = 2 * b[0] - a[0]
    dpy = 2 * b[1] - a[1]
    assert np.abs(dy - REF_Y).max() < 2e-6 and np.abs(dpy - REF_PY).max() < 2e-6
