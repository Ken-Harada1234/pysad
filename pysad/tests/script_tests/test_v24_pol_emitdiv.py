# -*- coding: utf-8 -*-
"""v2.4.0: POL ON の Emittance[] (スピン) の EMITDIV への依存。

条件: 4.2 GeV・RF 12 MV の 12 セル FODO (`test_v2_B3_pol._lat`)、v2、四極の縦ずれ ±3 mm。
測定 (線形の脱偏極率 depol_rate_linear): EMITDIV = 1 で 2.3138e-8、2 で 2.3138e-8、4 で 2.2492e-8、16 で 1.958e-8、64 で 1.868e-8、200 で 1.886e-8。
四極の刻み (EMITDIVQ) には依らず、偏向磁石の刻みだけに依る。刻みごとのスピン・軌道の結合 (d ベクトル) を刻みの端の値で積算する誤差で、
EMITDIV = 1 は収束値より +24% 大きい。v2 の既定 (POL ON) は 16 (残りの誤差 4%)。EMITDIV = 64 以上では εx が壊れる (200 で 2 倍) ので上げすぎない。"""
import io
import os
import sys

import pytest

from pysad import compat, emit6d, ffscmd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_v2_B3_pol import _lat  # noqa: E402


def _lin(div):
    ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12), out=io.StringIO())
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [div, None, None] if div is not None else [None, None, None]
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    finally:
        emit6d.EMITDIV[:] = old
    return e['spin']['depol_rate_linear']


def test_v2_pol_default_emitdiv_is_16(v2_mode, monkeypatch):
    monkeypatch.delenv('PYSAD_EMITDIV_DEFAULT', raising=False)
    assert compat.emitdiv_default(True) == 16.0 and compat.emitdiv_default(False) == 4.0
    d = _lin(None)
    r16 = _lin(16.0)
    r1 = _lin(1.0)
    assert d == pytest.approx(r16, rel=1e-12)
    assert r1 == pytest.approx(2.40686e-8, rel=2e-4)      # v2.4.0 D (B165): 2.3138e-8 → 2.3423e-8 (+1.2%)、D3 (B230: 四極の本体の B∥): → 2.4076e-8 (+2.8%)
    assert 0.83 < r16 / r1 < 0.87      # 0.848


def test_k64_emitdiv_default_is_1(k64_mode):
    assert compat.emitdiv_default(True) == 1.0 and compat.emitdiv_default(False) == 1.0
