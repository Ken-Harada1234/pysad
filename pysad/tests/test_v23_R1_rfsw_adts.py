# -*- coding: utf-8 -*-
"""v2.3.0 R1: gs (SAD K64) の TrackParticles で、RF つきのラティスの 6 極 ×0 の ADTS が pysad (k64、rfsw=False) と食い違った原因。
原因: FFS に入ると RFSW は ON (track.f:146-149) なので、gs の TrackParticles は RF つきのラティスで空洞を働かせる。1 周内の z (経路長の振幅依存) に応じて
空洞が δ を与え、(自然クロマティシティ) × δ(J) の偽の ADTS が入る。台本に NORFSW を入れると pysad (rfsw=False) と一致する。
基準値は gs の出力 (example/v2_3/R1_adts_gs_astrid_norfsw.txt・R1_adts_gs_astrid_rfsw.txt、1024 周、J = 2e-8、ASTRID2 0.58 GeV 45.7 m、6 極 ×0・8 極なし)。
(gs の値は 1 周ずつ z を 0 に戻して回した追跡の NAFF。pysad は同じ方法)"""
import os

import numpy as np
import pytest

from pysad import resdrive as R
from pysad.parser import load_sad_file
from pysad.track6n import track_particles_fast

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
ASTRID = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
P = 0.58e9
GS_NORFSW = np.array([83.20202, 235.512949, 238.574908])        # gs、6 極 ×0、NORFSW、1024 周
GS_RFSW = np.array([94.861447, 257.642491, 280.591272])         # gs、6 極 ×0、RFSW ON (FFS の既定)、1024 周


def _adts(rfsw, nturn=1024):
    lat = load_sad_file(ASTRID).build_lattice('RING')
    r0 = lat.twiss_table()[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    for e in lat.elements:
        if type(e).__name__ == 'Sextupole':
            e.K2 = 0.0

    def trk(c):
        c = c.copy()
        c[:, 4] = 0.0                  # gs の台本と同じ: 1 周ごとに z を 0 に戻す
        return track_particles_fast(lat.elements, c, 1, P, 0.0, rfsw, False, False)
    a = R.adts_tracking(lat.elements, P, tw, jref=(2e-8, 2e-8), nturn=nturn, tracker=trk)
    return np.array([a['xx'], a['xy'], a['yy']])


@pytest.mark.skipif(not os.path.exists(ASTRID), reason='ASTRID2 のラティスがない')
def test_astrid2_no_sextupole_adts_without_rf_matches_gs_norfsw():
    """RF なし (rfsw=False) の追跡は gs の NORFSW と 1e-6 で一致 (原因を除くと食い違いは消える)"""
    a = _adts(False)
    assert np.all(np.abs(a / GS_NORFSW - 1) < 1e-6), (a, GS_NORFSW)


@pytest.mark.skipif(not os.path.exists(ASTRID), reason='ASTRID2 のラティスがない')
def test_astrid2_rf_on_adts_reproduces_gs_default_rfsw_on():
    """RF ON (rfsw=True、z を毎周 0 に戻す) の追跡は gs の既定 (RFSW ON) の値 (94.86 ...) と 1e-3 で一致し、RF なしの値 (83.2) から 5% 以上ずれる。
    pysad の RF の追跡は gs と同じ。食い違いは台本の RFSW の扱いで、K64 の誤りではない (差は NAFF の窓の 1 周のずれ)"""
    a = _adts(True)
    assert np.all(np.abs(a / GS_RFSW - 1) < 1e-3), (a, GS_RFSW)
    assert np.all(np.abs(a / GS_NORFSW - 1) > 0.05)
