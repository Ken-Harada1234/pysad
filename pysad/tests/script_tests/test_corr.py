"""軌道補正 (CorrectOrbit・Monitor・Steer = Packages/correction.n・corutility.n を pysad が SADScript のまま実行、ResponseMatrix・LinearSolveConditioned は corr.py)。
参照は SAD K64 (gs)。SAD の Packages が見つからなければ skip。"""
import io
import os
import re

import pytest

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _nums(s):
    return [float(t) for t in re.findall(_NUM, s)]


def test_correct_orbit_matches_gs():
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('SAD の Packages がない')
    from pysad.ffscmd import run_sad_string
    src = open(os.path.join(os.path.dirname(__file__), 'co5_corr.sad'), encoding='utf-8').read()
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    d = {}
    for l in buf.getvalue().splitlines():
        k, _, v = l.partition(' ')
        if k in ('TUN', 'RM0', 'A', 'R', 'RX', 'RM'):
            d[k] = v
    ref = {'TUN': [2.8173851815, 1.9637233095], 'RM0': [1.57207490e-5, 1.34249935e-4],
           'A': [24, 24, .0059489760, -10.977135647, -.0024236166, -4.8041058095],
           'R': [24, .1492707525, -2.1547311e-6, .1998077673, -3.2440806e-5, .0721791517, -5.8430972e-5],
           'RX': [24, -4.1636051e-4, .8831154312, 1.47288518e-4, -2.5537930504, 1.28814950e-4, 2.2358640256],
           'RM': [24, 0, 0, 0, 0, 0, 3.51893198e-5]}
    for k, rv in ref.items():
        got = _nums(d[k])[:len(rv)]
        assert len(got) == len(rv), (k, d[k])
        for a, b in zip(got, rv):
            assert abs(a - b) <= 1e-6 * max(abs(b), 1e-4) + 1e-9, (k, a, b)


def test_linear_solve_conditioned_and_response():
    import numpy as np
    from pysad.corr import linear_solve_conditioned
    a = np.array([[1., 2.], [3., 4.], [5., 6.5]])
    b = np.array([1., 2., 3.])
    x = linear_solve_conditioned(a, b)
    assert np.allclose(x, np.linalg.lstsq(a, b, rcond=None)[0], atol=1e-7)
    xc = linear_solve_conditioned(a, b, 1e-8, 0, True, [[1., 1.]], [1.])
    assert abs(xc[0] + xc[1] - 1.0) < 1e-9
    xm = linear_solve_conditioned(a, b, 1e-8, 1, True)
    assert (xm != 0).sum() == 1


def test_optimize_optics_matches_gs():
    """OptimizeOptics (Packages/optimize.n を SADScript のまま実行) と FFS["VAR"]。参照は gs"""
    from pysad import pkgload
    if pkgload.packages_dir() is None:
        pytest.skip('SAD の Packages がない')
    from pysad.ffscmd import run_sad_string
    src = open(os.path.join(os.path.dirname(__file__), 'oo1_optimize.sad'), encoding='utf-8').read()
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    d = {}
    for l in buf.getvalue().splitlines():
        k, _, v = l.partition(' ')
        d[k] = v
    ref = {'T0': [2.7067578519, 1.8668722975], 'R': [3, 1.0057934346],
           'K': [.7059715383, -.5989199800], 'T1': [2.7386574736, 1.8457460010]}
    assert abs(_nums(d['V'])[0] - 2) < 1e-9
    for k, rv in ref.items():
        got = _nums(d[k])
        for a, b in zip(got, rv):
            assert abs(a - b) <= 1e-7 * max(abs(b), 1e-3), (k, a, b)
