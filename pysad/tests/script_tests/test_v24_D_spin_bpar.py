# -*- coding: utf-8 -*-
"""v2.4.0 担当 D: 磁石の本体の運動量の向きの磁場 (1+a)B∥ ds (docs/30 B165。B136・B137)。

独立な基準 `_spin_maxwell_indep.maxwell_R`: Maxwell 方程式を満たす硬い端の扇形 BEND (15°、L = 1.2 m、0.5 GeV、K1 = 0、E1 = E2 = 0) の中の
Thomas-BMT の数値積分。スピン回転行列 R (入口のフレーム → 出口のフレーム) の py に関する 1 次の係数 dR/dpy (py = ±2e-4 の差)。
基準値 (w = 2e-3 m。w = 4e-3・1e-3・5e-4 で 1e-6 以内):
    dR01/dpy = −0.554514、dR10/dpy = +0.543079、dR12/dpy = −0.1203、dR21/dpy = −0.0439。
v2 (B165 を入れる) の追跡 (BEND 1 要素) と Emittance[] (EMITDIV = 256) の差は、R01 で追跡 −8e-4・Emittance[] −4e-4 (K64 のまま: −0.26)。
追跡の R12・R21 の差は 1e-2 (K64 のまま: 0.05・0.03。刻みの中点の向きで端の ψ を分解する K64 の近似)。"""
import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, track6

import io
import math

GS = 0.001159652181280002
E, L, TH, PY = 0.5, 1.2, 0.2617993877991494, 2e-4
REF = np.array([[0.0, -0.554514, 0.0], [0.543079, 0.0, -0.1203], [0.0, -0.0439, 0.0]])


def _setup():
    src = f'MOMENTUM = {E} GEV;\nDRIFT D1=(L=1);\nBEND E1 = (L={L} ANGLE={TH} E1=0 E2=0);\nMARK ST = ();\nLINE RING = (ST E1 D1);\nFFS;\nUSE RING;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    els = [e for e in ffs.lat.elements if e.name == 'E1']
    return els, ffs.momentum()


def _R_track(els, p, py):
    def one(s0):
        out, s1 = track6.track_particles(els, np.array([[0, 0, 0, py, 0, 0.0]]), 1, p, 0.0, False, True, False, spins=[s0], radpol=False, circ=100.0)
        return s1[0]
    return np.array([one(np.eye(3)[j]) for j in range(3)]).T


def _lin(f):
    return (f(PY) - f(-PY)) / (2 * PY)


def test_bend_body_bpar_tracking_v2():
    els, p = _setup()
    compat.set_mode('v2')
    try:
        compat.set_spin_bpar(False)
        d0 = _lin(lambda py: _R_track(els, p, py)) - REF
        compat.set_spin_bpar(True)
        d1 = _lin(lambda py: _R_track(els, p, py)) - REF
    finally:
        compat.set_spin_bpar(None)
        compat.set_mode(None)
    assert abs(d0[0, 1]) > 0.2 and abs(d0[1, 0]) > 0.2          # B∥ なし (K64・従来の v2): 厳密な値と 0.26 違う
    assert abs(d1[0, 1]) < 2e-3 and abs(d1[1, 0]) < 8e-3 and abs(d1[1, 2]) < 2e-2 and abs(d1[2, 1]) < 2e-2


def test_bend_body_bpar_emittance_v2():
    els, p = _setup()
    emit6d.POL_MODE = True
    compat.set_mode('v2')
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [256, 1, 1]
    try:
        def em(py):
            st = emit6d.tturne(els, np.array([0, 0, 0, py, 0, 0.0]), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
            return st.srot[:, :3].copy()
        compat.set_spin_bpar(True)
        d1 = _lin(em) - REF
        compat.set_spin_bpar(False)
        d0 = _lin(em) - REF
    finally:
        emit6d.EMITDIV[:] = old
        compat.set_spin_bpar(None)
        compat.set_mode(None)
        emit6d.POL_MODE = False
    assert abs(d0[0, 1]) > 0.2
    assert np.abs(d1).max() < 2e-3


def test_k64_spin_unchanged():
    """k64: 本体の B∥ は入らない (spin_bpar() は常に False。K64 の追跡の dR01/dpy = −0.8157 ... 厳密 −0.5545 との差 0.26)"""
    els, p = _setup()
    compat.set_mode('k64')
    try:
        assert compat.spin_bpar() is False
        d = _lin(lambda py: _R_track(els, p, py)) - REF
    finally:
        compat.set_mode(None)
    assert abs(d[0, 1]) > 0.2


@pytest.mark.skipif(not __import__('os').environ.get('PYSAD_SLOW'), reason='PYSAD_SLOW=1 (独立な積分、約 1 分)')
def test_reference_regeneration():
    from . import _spin_maxwell_indep as M
    R = (M.maxwell_R(E * 1e9, L, TH, PY, w=4e-3, zs=2.4e-2) - M.maxwell_R(E * 1e9, L, TH, -PY, w=4e-3, zs=2.4e-2)) / (2 * PY)
    assert np.abs(R - REF).max() < 2e-4
