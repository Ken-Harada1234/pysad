# -*- coding: utf-8 -*-
"""F4: BOUND=0 の SOL (BZ≠0) が領域の外にあるときの追跡。K64 の tsol は BOUND≠0 でなくてもその SOL を領域の始点として、次の BOUND≠0 の SOL まで
ソレノイド領域として追う (幾何 tsgeo は BOUND≠0 だけを始点とみなして ` Missing BOUND of Solenoid` を出す)。gs (K64) は追跡では BOUND=1 の始点と全桁一致。
example/07_spin_solenoid_special/elem_SOL_bound0_track.sad (1.7 GeV、BZ=0.8 T、QUAD ×4、3 粒子 1 周 RAD OFF) の 6 座標 18 個が gs と 10 桁一致 (両モード)。
Twiss・CALC は gs が SIGSEGV で基準が無く未対応のまま。"""
import io
import os
import re

import numpy as np
import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), 'example', '07_spin_solenoid_special', 'elem_SOL_bound0_track.sad')


def _res(text):
    rows = []
    for ln in text.splitlines():
        if ln.startswith('RES'):
            rows.append([float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?', ln[3:].replace('{', ' ').replace('}', ' '))])
    return np.array(rows)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_bound0_region_matches_gs(mode):
    ref = _res(open(SAD.replace('.sad', '.gs.txt'), encoding='utf-8', errors='ignore').read())
    assert ref.shape == (6, 3)
    old = compat.mode() if hasattr(compat, 'mode') else 'k64'
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        ffscmd.run_sad_string(open(SAD, encoding='utf-8').read(), out=buf)
    finally:
        compat.set_mode(old)
    got = _res(buf.getvalue())
    assert got.shape == ref.shape
    assert np.max(np.abs(got - ref)) < 1e-9 * np.max(np.abs(ref)) + 1e-15
