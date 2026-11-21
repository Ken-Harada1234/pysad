# -*- coding: utf-8 -*-
"""v2.4.0 A5: DynamicApertureSurvey の ExpandElementValues->False と ReferenceOrbit->{x0, px0, y0, py0, z0, dp0}。
gs (K64) は ExpandElementValues->True (既定) で成分ごとの値 (補正磁石の K0・位置ずれ) を代表の値に揃える (help: 'Machine errors may be reset')。False なら揃えない。
DDBA v07_14 の誤差の例 (example/ddba/v07_14/errors) で gs の DA は False で 34σ・29σ (True では 0 や 1σ: 補正した COD が消える)。pysad は以前 False を読まず、
常に揃えていた (gs と DA が違った原因。docs/history/v24_A.md)。ReferenceOrbit は dapert の追跡の初期軌道。
リング: 小さな FODO (3 GeV、24 セル、6 極 K2 = 4・−6)、k64。"""
import contextlib
import io
import os

from pysad.ffscmd import run_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')


def _run(body):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        run_sad_string(FODO + '\nFFS;\nUSE RING NOEXPAND;\ncc = "cell calc noexpand"; ccemit = "cell calc emit noexpand";\nFFS[ccemit];\n' + body + '\nstop;\n')
    return buf.getvalue().split('\n')


def test_expand_element_values_false_keeps_component_values():
    """1 台だけ変えた QF.3 の DX: 既定 (True) の DynamicApertureSurvey のあとは代表 (QF.1) の値 0 に揃う。False は 1e-4 のまま (gs と同じ)"""
    out = {}
    for opt in ('', ',ExpandElementValues->False', ',ExpandElementValues->True'):
        lines = _run('LINE["DX","QF.3"] = 1.0e-4;\nFFS[cc];\nres = DynamicApertureSurvey[{{1,10},{1,1},{0}}, 50, Output->0' + opt + '];\n'
                     'Print["AFTER ", LINE["DX","QF.3"]];')
        out[opt] = [l for l in lines if l.startswith('AFTER')][0].split()[1]
    assert float(out['']) == 0.0 and float(out[',ExpandElementValues->True']) == 0.0
    assert abs(float(out[',ExpandElementValues->False']) - 1e-4) < 1e-12


def test_reference_orbit_shifts_the_survey():
    """ReferenceOrbit->{0,0,0,0,0,0} は理想のリング (閉軌道 0) の既定と同じスコア。大きなずれ (x0 = 20 mm) の軌道のまわりでは違う (ずれた軌道のまわりの追跡)"""
    def score(opt):
        lines = _run('FFS[cc];\nLOSSAMPL = 0.02;\nres = DynamicApertureSurvey[{{1,60},{1,1},{0}}, 50, Output->0' + opt + '];\nPrint["SCORE ", res[1]];')
        return float([l for l in lines if l.startswith('SCORE')][0].split()[1])
    a = score('')
    b = score(',ReferenceOrbit->{0,0,0,0,0,0}')
    c = score(',ReferenceOrbit->{0.015,0,0,0,0,0}')
    assert a == b
    assert c != a
