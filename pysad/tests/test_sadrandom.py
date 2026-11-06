"""SAD の既定乱数 (random_plugin_sad.c) と Random[]/GaussRandom[] の束縛の確認"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from pysad.sadrandom import SadRandom  # noqa: E402

A, M = 48828125, 4294967296


def _lcg(s):
    s = (A * s) % M
    return s - M if s >= 2147483648 else s


def test_uniform0_formula():
    # Random[] は RandomUniform0 (tfRandom_.c REG8("Random", RandomUniform0)): [0,1) = (seed + 2^31)/2^32
    r = SadRandom(187)
    s = r.seed
    for _ in range(5):
        s = _lcg(s)
        assert r.uniform0() == (s + 2147483648.0) / M
    # RANDOM_UNIFORM (0,1) は +0.5 (Random[] には使わない)
    s = _lcg(s)
    assert r.uniform() == (s + 2147483648.5) / M


def test_random_binding_is_uniform0():
    from pysad.sadscript import Interp
    I = Interp()
    I.rng = SadRandom(187)
    ref = SadRandom(187)
    v = I.run('Random[]', 'script')
    assert isinstance(v, float) and v == ref.uniform0()


def test_seed_list_matches_sad():
    # print_all_error.sad: SeedRandom[187]; Floor[Abs[Random[1000]*1e8]]*2-1 の先頭 (SAD の出力 esrf.print_all_error.sad.out)
    r = SadRandom(187)
    seeds = [int(abs(r.uniform0() * 1e8)) * 2 - 1 for _ in range(1000)]
    assert seeds[:3] == [125188771, 60490425, 145796967]
    # 以前 +0.5 の式だったときに 2 ずれていた seed (SAD: 19334753 ほか)
    assert 19334753 in seeds and 185519945 in seeds and 165373679 in seeds and 62119043 in seeds
