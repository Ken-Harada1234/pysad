# -*- coding: utf-8 -*-
"""TrackParticles の並列化 (NWORKERS・PYSAD_WORKERS、track6n.ParallelTrackRun) と TrackSetup.track(workers=) の試験。

粒子を連続した塊に分け、塊ごとの TrackRun をスレッドで同時に回す。粒子は互いに独立なので、FLUC なしの結果は逐次とビット単位で同じ。
FLUC の乱数は塊ごとに種を変えるので逐次とは別の列になる (統計は同じ)。粒子どうしの相互作用や Python のコールバックがあるものは並列にしない。
"""
import contextlib
import io
import os

import numpy as np
import pytest

from pysad.ffscmd import run_sad_string
from pysad import track6n

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')

NEED = max(1, int((os.cpu_count() or 1) * 0.9))


def _coords_script(nworkers, extra, np_=96, turns=30):
    """TrackParticles の戻りの座標 (6 行 + alive) を numpy 配列で返す"""
    body = ('FFS;\nUSE RING NOEXPAND;\nCELL; CALC;\n%s\nNWORKERS=%d;\nnp=%d;\nSeedRandom[1];\n'
            'beam={Table[2e-3*Random[],{np}],Table[0,{np}],Table[1e-3*Random[],{np}],Table[0,{np}],Table[0,{np}],Table[1e-3*Random[],{np}],Table[1,{np}]};\n'
            'out=TrackParticles[{1,beam},1,1,%d];\n'
            'Do[Print["ROW ",i," ",Table[ToString[out[[2,i,k]],FormatType->InputForm],{k,1,np}]], {i,1,7}];\nstop;\n') % (extra, nworkers, np_, turns)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        run_sad_string(FODO + '\n' + body)
    rows = {}
    for line in buf.getvalue().split('\n'):
        if line.startswith('ROW '):
            parts = line.split(' ', 2)
            rows[int(parts[1])] = parts[2]
    return rows


@pytest.mark.skipif(NEED < 2, reason='CPU が 1 コアだけ')
def test_ffs_parallel_equals_serial_without_fluc():
    a = _coords_script(1, 'NORFSW; NOFLUC;')
    b = _coords_script(2, 'NORFSW; NOFLUC;')
    assert set(a) == set(range(1, 8)) and a == b          # 17 桁まで同じ (ビット単位)


@pytest.mark.skipif(NEED < 2, reason='CPU が 1 コアだけ')
def test_ffs_parallel_with_rf_and_apert_loss_equals_serial():
    """RF ON・損失 (LOSSAMPL を小さくして粒子を落とす) でも逐次と同じ。落ちた粒子の座標・alive も"""
    ex = 'RFSW; NOFLUC; LOSSAMPL=3e-3;'
    a = _coords_script(1, ex, turns=40)
    b = _coords_script(2, ex, turns=40)
    assert a == b
    assert '0' in a[7] and '1' in a[7]            # 落ちた粒子と残った粒子の両方がある (alive 行)


def test_max_workers_caps():
    assert track6n.max_workers(1) == 1
    assert track6n.max_workers('x') == 1
    assert track6n.max_workers(10 ** 6) == NEED
    assert track6n.max_workers(8, nitems=20) == 1             # 1 本あたり 16 粒子未満
    assert track6n.max_workers(10 ** 6, nitems=64) <= 4


def test_split_blocks_cover_all_particles_in_order():
    for n, w in ((10, 3), (7, 7), (5, 9), (96, 4)):
        b = track6n.split_blocks(n, w)
        assert b[0][0] == 0 and b[-1][1] == n
        assert all(b[i][1] == b[i + 1][0] for i in range(len(b) - 1))
        assert max(y - x for x, y in b) - min(y - x for x, y in b) <= 1


@pytest.mark.skipif(NEED < 2, reason='CPU が 1 コアだけ')
def test_tracksetup_workers_equals_serial():
    from pysad.trksetup import TrackSetup
    ts, _ffs = TrackSetup.from_source(FODO, 'RING', rad=False, rfsw=False)
    rng = np.random.default_rng(3)
    c = np.zeros((64, 6))
    c[:, 0] = rng.uniform(-4e-3, 4e-3, 64)
    c[:, 2] = rng.uniform(-2e-3, 2e-3, 64)
    r1 = ts.track(1, c, 50, workers=1)
    r2 = ts.track(1, c, 50, workers=2)
    assert np.array_equal(r1['lost_turn'], r2['lost_turn'])
    assert np.array_equal(r1['final'], r2['final'])
    assert np.array_equal(r1['alive'], r2['alive'])


def test_ffs_workers_ignored_for_small_particle_count():
    """粒子が 16 個未満なら並列にしない (逐次と同じ)"""
    a = _coords_script(1, 'NORFSW; NOFLUC;', np_=8, turns=5)
    b = _coords_script(4, 'NORFSW; NOFLUC;', np_=8, turns=5)
    assert a == b


@pytest.mark.skipif(NEED < 2, reason='CPU が 1 コアだけ')
def test_ffs_parallel_path_is_really_used_and_not_for_interacting_flags(monkeypatch):
    used = []
    orig = track6n.ParallelTrackRun.segment

    def spy(self, a, b):
        used.append(len(self.runs))
        return orig(self, a, b)

    monkeypatch.setattr(track6n.ParallelTrackRun, 'segment', spy)
    _coords_script(2, 'NORFSW; NOFLUC;', np_=64, turns=3)
    assert used and max(used) == 2
    used.clear()
    _coords_script(1, 'NORFSW; NOFLUC;', np_=64, turns=3)          # NWORKERS=1 は逐次
    assert not used
    used.clear()
    _coords_script(2, 'NORFSW; NOFLUC; SPAC;', np_=64, turns=3)    # 空間電荷は並列にしない
    assert not used


@pytest.mark.skipif(NEED < 2, reason='CPU が 1 コアだけ')
def test_ffs_parallel_fluc_statistics_match_serial():
    """FLUC (RAD) の乱数は塊ごとに種が違い、逐次とは別の列。ただし分布 (δ の rms) は同じ程度"""
    ex = 'RFSW; RAD; FLUC;'
    def rms(nw):
        rows = _coords_script(nw, ex, np_=128, turns=200)
        v = np.array([float(s) for s in rows[6].replace('{', '').replace('}', '').replace('"', '').split(',') if s.strip()])
        return float(np.std(v))
    r1, r2 = rms(1), rms(2)
    assert r1 > 0 and r2 > 0 and abs(r2 / r1 - 1.0) < 0.5
