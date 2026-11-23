"""tdgeo.f (幾何の量 GX GY GZ CHI1..3 の、DRIFT の L・BEND の ANGLE・L に対する解析微分) の直訳 match6d.tdgeo_term。
DRIFT の L と BEND の ANGLE (ROTATE 付きも) は、幾何を有限差分したものと一致する。BEND の L は K64 が出射の向き geo(:,3,k+1) で近似する式なので、
有限差分とは O(ANGLE^2) 違い、K64 の式そのものと一致することを確かめる。"""
import io

import numpy as np

SRC = r"""MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7);
BEND B1 = (L=1.2 ANGLE=0.1308996938995747) B2 = (L=1.0 ANGLE=0.09 ROTATE=0.3 DEG) B3 = (L=0.8 ANGLE=-0.05 ROTATE=20 DEG);
LINE RING = (ST QF D1 B1 D2 B2 D1 B3 D2 B1 QF);
FFS;
USE RING;
Exit;
"""


def _lat():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _I = run_sad_string(SRC, out=buf)
    return ffs.lat.elements


def _val(g, kname, ne):
    if kname[0] == 'G':
        return float(g.gxyz(ne)['XYZ'.index(kname[1])])
    return float(g.chi(ne)[int(kname[3]) - 1])


def _numeric(els, k, key, kname, ne, eps=1e-6):
    from pysad import geo as G
    el = els[k - 1]
    old = getattr(el, key)
    raw = getattr(el, '_raw', None)
    oldr = raw.get(key) if isinstance(raw, dict) else None
    vals = []
    for sg in (1.0, -1.0):
        setattr(el, key, old + sg * eps)
        if oldr is not None:
            raw[key] = oldr + sg * eps
        vals.append(_val(G.compute_geo(els), kname, ne))
    setattr(el, key, old)
    if oldr is not None:
        raw[key] = oldr
    return (vals[0] - vals[1]) / (2 * eps)


def test_tdgeo_matches_finite_differences_for_drift_l_and_bend_angle():
    from pysad import geo as G
    from pysad.match6d import tdgeo_term
    els = _lat()
    n = len(els)
    ne = n + 1
    geo = G.compute_geo(els).geo
    checked = 0
    for k in range(1, n + 1):
        el = els[k - 1]
        tn = type(el).__name__
        keys = ('L',) if tn == 'Drift' else ('ANGLE',) if tn == 'Bend' else ()
        for key in keys:
            for kname in ('GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3'):
                if key == 'L' and kname.startswith('CHI'):
                    continue
                an = tdgeo_term(geo, kname, ne, k, el, tn, key, 1.0)
                nd = _numeric(els, k, key, kname, ne)
                assert abs(an - nd) < 1e-5 * max(1.0, abs(nd)), (k, tn, key, kname, an, nd)
                checked += 1
    assert checked > 30


def test_tdgeo_bend_length_is_the_k64_exit_direction_formula():
    from pysad import geo as G
    from pysad.match6d import tdgeo_term
    els = _lat()
    n = len(els)
    geo = G.compute_geo(els).geo
    for k, el in enumerate(els, start=1):
        if type(el).__name__ != 'Bend':
            continue
        for j, kname in enumerate(('GX', 'GY', 'GZ')):
            an = tdgeo_term(geo, kname, n + 1, k, el, 'Bend', 'L', 1.0)
            assert abs(an - geo[k, j, 2]) < 1e-14          # geo(j, 3, k+1)
