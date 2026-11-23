# -*- coding: utf-8 -*-
"""F2: PHOTONS (PhotonList) と COORD・TCAVI を同じラティスに置く (track6 が COORD・TCAVI を通す)。
K64 は tcoord・ttcave で光子を記録しない (BEND・QUAD などだけ)。gs (K64) との照合: example/07_spin_solenoid_special/flag_PHOTONS_coord_tcavi.sad
(3 GeV、BEND ×4 + COORD の対 ×2 + TCAVI ×2、4 粒子 8 周、RAD・FLUC) の光子 278 個が、個数・全 12 列 (gs の出力桁 1e-10 まで) で一致。
MAP を含むラティスは未対応 (ExternalMap の呼び出しが TrackRun の状態に依る)。"""
import io
import os

import numpy as np

from pysad import ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SAD = os.path.join(ROOT, 'example', '07_spin_solenoid_special', 'flag_PHOTONS_coord_tcavi.sad')


def _parse_gs(path):
    a, b = [], []
    for ln in open(path, encoding='utf-8', errors='ignore'):
        if ln.startswith('PA'):
            a.append([float(v) for v in ln.split()[1:]])
        elif ln.startswith('PB'):
            b.append([float(v) for v in ln.split()[1:]])
    return np.hstack([np.array(a), np.array(b)])


def _run(src):
    _ffs, I = ffscmd.run_sad_string(src, out=io.StringIO())
    return np.array(I.globals['PhotonList'], dtype=float).reshape(-1, 12)


def test_matches_gs(k64_mode):
    ref = _parse_gs(SAD.replace('.sad', '.gs.txt'))
    ph = _run(open(SAD, encoding='utf-8').read())
    assert ph.shape == ref.shape == (278, 12)
    assert np.max(np.abs(ph - ref)) < 2e-9
    assert np.array_equal(ph[:, 10:], ref[:, 10:])        # 粒子番号・要素番号


def test_v2_records_photons():
    from pysad import compat
    compat.set_mode('v2')
    try:
        ph = _run(open(SAD, encoding='utf-8').read())
    finally:
        compat.set_mode('k64')
    assert 200 < ph.shape[0] < 360 and ph.shape[1] == 12        # v2 は乱数の引き方が違うので個数は統計の範囲 (278 ± 18%)
