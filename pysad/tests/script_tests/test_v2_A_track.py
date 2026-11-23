# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): 追跡の B14 (tsolqum の dv の項) と B08 (tsolqur が最後の粒子の bzp を全粒子に使う)。

根拠:
- B14: 同じ K1・L の QUAD と MULT (K1 だけ) は同じ物理なので、放射なしの追跡の z が一致しなければならない。K1→0 の MULT はドリフトに一致する。
  K64 の tsolqum は z の式の dv の項が -0.25·(ap - dv)·aln0 (= +0.25·dv·aln0) で、tsolqu の -dv·aln0 と 1.25·dv·L 違う。
- B08: 1 粒子の追跡結果は、一緒に追う他の粒子に依らない。K64 の tsolqur (RAD、BZ≠0 の SOL 領域内の四極) は最後の半ステップで
  最後の粒子の bzp = BZ/(1+δ) を全粒子に使うので、他の粒子の δ で結果が変わる。

条件: B14 は 20 MeV (dv = δ/γ² を大きくするため)、L=0.5 m、K1=0.3。B08 は 3 GeV、BZ=2 T、四極 L=0.3 m K1=-0.6。
"""
import io

import numpy as np
import pytest

from pysad import ffscmd, track6, track6n

V = [1e-3, 2e-4, -5e-4, 1e-4, 0.0, 5e-2]


def _through(elem, fn, v, momentum='0.02GEV'):
    src = 'MOMENTUM = %s;\nDRIFT D1 = (L=1);\n%s\nMARK ST = ();\nLINE RING = (ST E1);\nFFS;\nUSE RING;\n' % (momentum, elem)
    ffs, _I = ffscmd.run_sad_string(src, out=io.StringIO())
    els = [e for e in ffs.lat.elements if e.name == 'E1']
    out = fn(els, np.array([v]), 1, ffs.momentum(), 0.0, True, False, False, circ=100.0)
    return np.asarray(out)[0][:6]


@pytest.mark.parametrize('fn', [track6n.track_particles_fast, track6.track_particles])
def test_mult_z_equals_quad_and_drift_v2(v2_mode, fn):
    q = _through('QUAD E1 = (L=0.5 K1=0.3 DISFRIN=1);', fn, V)
    m = _through('MULT E1 = (L=0.5 K1=0.3 DISFRIN=1);', fn, V)
    assert np.abs(q - m).max() < 1e-12
    d = _through('DRIFT E1 = (L=0.5);', fn, V)
    m0 = _through('MULT E1 = (L=0.5 K1=1e-12 DISFRIN=1);', fn, V)
    assert abs(d[4] - m0[4]) < 1e-13
    assert abs(d[4]) > 1e-5                    # dv·L の項そのものは 1.5e-5 m ある


@pytest.mark.parametrize('fn', [track6n.track_particles_fast, track6.track_particles])
def test_mult_z_k64_differs_by_1p25_dv_l(k64_mode, fn):
    q = _through('QUAD E1 = (L=0.5 K1=0.3 DISFRIN=1);', fn, V)
    m = _through('MULT E1 = (L=0.5 K1=0.3 DISFRIN=1);', fn, V)
    d = _through('DRIFT E1 = (L=0.5);', fn, V)
    px, py = V[1], V[3]
    dvl = d[4] + 0.5 * (px * px + py * py) * 0.5      # ドリフトの z = -(px²+py²)/2·L - dv·L (内部座標の px, py は 1/(1+δ) 倍: 2 次の差は無視)
    assert np.abs(q[:4] - m[:4]).max() < 1e-12
    assert abs((q[4] - m[4]) / (1.25 * d[4]) - 1.0) < 2e-2


SOLRING = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=0.5);
QUAD QD = (L=0.3 K1=-0.6);
SOL SB = (BZ=2 BOUND=1 GEO=1) SE = (BOUND=1);
MARK ST = ();
LINE RING = (ST SB D1 QD D1 SE);
FFS;
USE RING;
'''
PA = [1e-3, 2e-4, -5e-4, 1e-4, 0.0, 1e-3]
PB = [2e-3, -1e-4, 5e-4, 2e-4, 0.0, -2e-2]


def _sol(fn, coords, rad):
    ffs, _I = ffscmd.run_sad_string(SOLRING, out=io.StringIO())
    return np.asarray(fn(ffs.lat.elements, np.array(coords), 1, ffs.momentum(), 0.0, True, rad, False, circ=100.0))[:, :6]


@pytest.mark.parametrize('fn', [track6n.track_particles_fast, track6.track_particles])
def test_solenoid_quad_rad_is_independent_of_other_particles_v2(v2_mode, fn):
    alone = _sol(fn, [PA], True)[0]
    pair = _sol(fn, [PA, PB], True)[0]
    assert np.abs(alone - pair).max() == 0.0
    # 放射ありと放射なしの差は放射の大きさ (1e-8 以下) だけ
    norad = _sol(fn, [PA, PB], False)[0]
    assert np.abs(pair[:4] - norad[:4]).max() < 1e-8


@pytest.mark.parametrize('fn', [track6n.track_particles_fast, track6.track_particles])
def test_solenoid_quad_rad_k64_depends_on_last_particle(k64_mode, fn):
    alone = _sol(fn, [PA], True)[0]
    pair = _sol(fn, [PA, PB], True)[0]
    assert np.abs(alone - pair).max() > 1e-6
