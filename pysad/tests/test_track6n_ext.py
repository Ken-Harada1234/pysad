"""
track6n (numba) の MULT / ソレノイド領域 (track6n_ext のカーネルを track6n.Prog/turn に統合したもの) の照合。
純 Python の track6.track_particles と 20 粒子 × 5 周で比較する:
  NORAD / RAD / FLUC とも座標差 1e-13 以下、FLUC は乱数の最終状態 (RNG.seed) も一致。
  L=0 の MULT (薄い多極) と SOL が黙って無視されないこと (K_NONE にならない)。
  tapert (損失判定) 付きの turn がソレノイド領域でも動くこと (損失なしで n 不変、x=5 m の粒子は落ちる)。
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

import numpy as np
import pytest

from pysad.elements import make_element
from pysad import track6 as t6
from pysad import track6n as tn
from pysad import emit6d as _e

A = 0.157079632679490


def _mult(name, L, knl, **kw):
    e = make_element('MULT', name, L=L, KnL=knl)
    for k, v in kw.items():
        setattr(e, k, v)
    return e


def _sol(name, bz, bound):
    e = make_element('SOL', name, L=0.0)
    e.BZ = bz
    e.BOUND = bound
    return e


def _lattice(case):
    ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
    ca.HARM = 12.0
    els = [make_element('MARK', 'ORG'), ca,
           _mult('MQF', 0.3, {1: 1.8 * 0.3}),
           make_element('DRIFT', 'DQS', L=0.1),
           _mult('MSF', 0.2, {0: 2e-4, 2: 0.6}),
           make_element('DRIFT', 'DS', L=0.9),
           make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
           _mult('MK0', 0.1, {0: 1e-4}),
           _mult('MTH', 0.0, {0: 5e-5, 2: 0.3}),            # 薄い多極 (L=0)
           _mult('MS2', 0.2, {2: 6.0}),                      # ndiv=4 (tpara の ndiv)
           make_element('DRIFT', 'D1S', L=1.0)]
    if case == 'sol':
        els += [_sol('SB', 2.0, 1.0),
                make_element('DRIFT', 'D1', L=0.1),
                make_element('QUAD', 'QD', L=0.3, K1=-1.7),
                _sol('SE', 2.0, 1.0)]
    else:
        els += [make_element('DRIFT', 'D1', L=0.1),
                _mult('MQD', 0.3, {1: -1.7 * 0.3, 3: 50.0}, DX=1e-4)]   # K3 付き (フリンジ dofr)
    els += [make_element('DRIFT', 'DQS', L=0.1),
            _mult('MSD', 0.2, {2: -1.2}, DY=1e-4, ROTATION=1e-3),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


def _coords(npart=20, seed=1):
    rng = np.random.default_rng(seed)
    c = np.zeros((npart, 6))
    c[:, 0] = rng.normal(0, 3e-4, npart)
    c[:, 1] = rng.normal(0, 1e-4, npart)
    c[:, 2] = rng.normal(0, 2e-4, npart)
    c[:, 3] = rng.normal(0, 1e-4, npart)
    c[:, 4] = rng.normal(0, 1e-3, npart)
    c[:, 5] = rng.normal(0, 5e-4, npart)
    return c


P0 = 2.5e9
NT = 5


@pytest.mark.parametrize('case', ['mult', 'sol'])
@pytest.mark.parametrize('mode', ['NORAD', 'RAD', 'FLUC'])
def test_track6n_vs_track6(case, mode):
    els = _lattice(case)
    coords = _coords()
    rad = mode != 'NORAD'
    fluc = mode == 'FLUC'
    # track_particles_fast が track6 にフォールバックしていないこと (Prog が compile できる)
    ph = _e._Phys(P0, sum(el.L for el in els))
    prog = tn.Prog(els, ph)
    assert not any(k == tn.K_NONE for k, el in zip(prog.kinds, els) if type(el).__name__ in ('Multipole', 'Solenoid'))
    ra = t6.track_particles(els, coords, NT, P0, trf0=1e-4, rad=rad, fluc=fluc, seed=12345)
    sa = t6.RNG.seed
    rb = tn.track_particles_fast(els, coords, NT, P0, trf0=1e-4, rad=rad, fluc=fluc, seed=12345)
    sb = t6.RNG.seed
    d = np.abs(ra - rb).max()
    assert d <= 1e-13, f'{case} {mode}: max|diff| = {d:.2e}'
    assert sa == sb, f'{case} {mode}: RNG.seed {sa} != {sb}'


def test_thin_mult_and_sol_not_ignored():
    """L=0 の MULT (薄い多極) と SOL を黙って無視しない: MULT/SOL を抜いた追跡と結果が違う"""
    els = _lattice('sol')
    coords = _coords()
    full = tn.track_particles_fast(els, coords, NT, P0, trf0=1e-4, rad=False, fluc=False)
    no_thin = [e for e in els if e.name != 'MTH']
    r1 = tn.track_particles_fast(no_thin, coords, NT, P0, trf0=1e-4, rad=False, fluc=False)
    assert np.abs(full - r1).max() > 1e-8
    no_sol = [e for e in els if type(e).__name__ != 'Solenoid']
    r2 = tn.track_particles_fast(no_sol, coords, NT, P0, trf0=1e-4, rad=False, fluc=False)
    assert np.abs(full - r2).max() > 1e-8


def test_apert_with_solenoid():
    """ソレノイド領域を含む turn(apert=...) : 損失なしで n 不変、x=5 m の粒子を 1 つ入れると n が 1 減って周番号が入る"""
    els = _lattice('sol')
    ph = _e._Phys(P0, sum(el.L for el in els))
    prog = tn.Prog(els, ph)
    assert any(sm == tn.SM_START for sm in prog.sol_mode) and any(sm == tn.SM_END for sm in prog.sol_mode)
    coords = _coords(10)
    rs = np.array([12345], dtype=np.int64)
    S = tn.make_state(coords, ph)
    kid = np.arange(10, dtype=np.int64)
    lost_turn = np.zeros(10, dtype=np.int64)
    n = tn.turn(prog, S, 10, ph, 1e-4, True, False, False, rs, apert=(kid, lost_turn, 1, 1.0, 100.0))
    assert n == 10
    coords2 = coords.copy()
    coords2[3, 0] = 5.0
    S = tn.make_state(coords2, ph)
    kid = np.arange(10, dtype=np.int64)
    lost_turn = np.zeros(10, dtype=np.int64)
    n = tn.turn(prog, S, 10, ph, 1e-4, True, False, False, rs, apert=(kid, lost_turn, 1, 1.0, 100.0))
    assert n == 9
    assert lost_turn[3] == 1 or 1 in lost_turn.tolist()


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q']))
