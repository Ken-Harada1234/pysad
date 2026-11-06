"""ROTATE (ROT/TILT) 引数への対応のテスト (docs/history/work_rotate.md)。

1. parser: ANGLE=0 の BEND (ROTATE=-90 DEG の垂直ステアリング VD/PV/VS) は Drift にせず Bend として保持する。
2. lattice._tchge_jac: SAD tchge.f の trans の作り方 (癖を含む) を、写像の数値微分と比べる。
3. SAD K64 の実行結果 (bench_level25a〜f.result, リポジトリ直下) と CALC の Twiss を突き合わせる
   (追跡・Emittance[] の手前までに切って走らせる。SAD は不要、保存済みの .result を参照する)。
"""
import math
import os
import re

import numpy as np
import pytest

from pysad import lattice as L
from pysad import parser as _parser
from pysad.elements import Bend

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?')
STRIP = re.compile(r'(?<![\d.])[A-Za-z_$^][\w.$^*]*')


def test_parser_keeps_zero_angle_rotated_bend():
    sc = _parser.load_sad_string(
        'MOMENTUM = 3GEV;\n BEND VD=(L=.15 ROTATE=-90 DEG) PV=(L=0 ROTATE=-90 DEG K0=1e-4) ;\n LINE T=(VD PV) ;\n')
    lat = sc.env.build_lattice('T')
    els = [e for e in lat.elements if e.name in ('VD', 'PV')]
    assert [type(e).__name__ for e in els] == ['Bend', 'Bend']
    assert abs(els[0].ROTATE + math.pi / 2) < 1e-12 and els[0].ANGLE == 0.0
    assert abs(els[1].K0 - 1e-4) < 1e-18


def _numjac(f, c, eps=1e-7):
    return L._num_jacobian(f, np.asarray(c, float), eps)


C0 = np.array([3e-4, 2e-5, -2e-4, 1e-5, 0.0, 1e-3])


@pytest.mark.parametrize('args', [
    (0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0),                       # 回転だけ
    (2e-4, -1e-4, 0.0, 0.0, 0.0, 0.0, 0.0),                    # ずれだけ (th=0)
    (2e-4, -1e-4, 0.0, 0.0, 0.0, 0.1, 1e-3),                   # BEND のずれ (phig≠0: tbrote のドリフト分)
    (2e-4, -1e-4, 0.0, 0.0, 1e-3, 0.1, 1e-3),                  # DROTATE
])
@pytest.mark.parametrize('ent', [True, False])
def test_tchge_jac_matches_map_when_no_quirk(args, ent):
    """th=0 か、ずれ・dtheta がないとき、_tchge_jac は _tchge_cod の微分に一致する (SAD の癖が出ない場合)"""
    dx, dy, dz, th, dth, alg, phig = args
    f = lambda c: L._tchge_cod(c, dx, dy, dz, th, dth, alg, phig, ent)
    J = L._tchge_jac(C0, dx, dy, dz, th, dth, alg, phig, ent)
    Jn = _numjac(f, C0)
    # z 行は tbgrote の近似式 (SAD の trans) なので x,px,y,py の行だけ比べる
    assert np.max(np.abs(J[:4, :] - Jn[:4, :])) < 5e-9


def test_tchge_jac_sad_quirk_loses_rotation_at_entrance():
    """SAD tchge の癖: 入口で th≠0 かつずれ≠0 のとき trans に回転が入らない (COD は回る)。出口は回転が入る"""
    th = 0.3
    Jin = L._tchge_jac(C0, 1e-4, 0.0, 0.0, th, 0.0, 0.0, 0.0, True)
    assert np.allclose(Jin[:4, :4], np.eye(4))
    Jex = L._tchge_jac(C0, 1e-4, 0.0, 0.0, th, 0.0, 0.0, 0.0, False)
    c, s = math.cos(th), math.sin(th)
    assert abs(Jex[0, 0] - c) < 1e-15 and abs(Jex[0, 2] - s) < 1e-15
    # ずれなしなら入口にも回転が入る
    Jin0 = L._tchge_jac(C0, 0.0, 0.0, 0.0, th, 0.0, 0.0, 0.0, True)
    assert abs(Jin0[0, 2] + s) < 1e-15


def test_tchge_exit_is_inverse_of_entrance_for_cod():
    for args in [(2e-4, -1e-4, 0.0, 0.4, 0.0, 0.1, 1e-3), (0.0, 0.0, 0.0, -1.5707963267948966, 0.0, 0.0, 0.0),
                 (1e-4, 2e-4, 0.0, 0.2, 0.0, 0.1, 1e-3)]:    # DROTATE≠0 は BEND 本体を挟まないと逆にならないので除く
        dx, dy, dz, th, dth, alg, phig = args
        c1 = L._tchge_cod(C0, dx, dy, dz, th, dth, alg, phig, True)
        al0, phi0 = (2 * alg, 2 * phig) if phig != 0.0 else (0.0, 0.0)
        c2 = L._tchge_cod(c1, dx, dy, dz, th, dth, alg - al0, phig - phi0, False)
        assert np.max(np.abs(c2 - C0)) < 1e-12


def test_vertical_steering_kick_is_rotated():
    """ANGLE=0, L=0, ROTATE=-90 DEG の BEND (PV): K0>0 の水平キック px −= K0 が垂直 (py) に出る"""
    el = Bend(name='PV', L=0.0, ANGLE=0.0, ROTATE=-math.pi / 2)
    el.K0 = 1e-4
    c, M = L._bend_element_cod_matrix(el, np.zeros(6))
    assert abs(c[1]) < 1e-12
    el0 = Bend(name='ST', L=0.0, ANGLE=0.0)
    el0.K0 = 1e-4
    c0, _ = L._bend_element_cod_matrix(el0, np.zeros(6))
    assert abs(c0[1] + 1e-4) < 1e-18 and abs(c0[3]) < 1e-18
    # 回転しても大きさは同じ、向きは R(-th)·(−K0, 0) = (−K0 cos th, +K0 sin th)
    th = -math.pi / 2
    assert abs(c[1] - (-1e-4) * math.cos(th)) < 1e-12 and abs(c[3] - 1e-4 * math.sin(th)) < 1e-12


# ----------------------------------------------------------------------------
# SAD の保存結果との照合 (CALC の Twiss)
# ----------------------------------------------------------------------------
def _rows(text, tags):
    out = []
    lines = text.replace('\r', '').split('\n')
    cur = None
    for ln in lines:
        w = ln.split(None, 1)
        if w and w[0] in tags:
            cur = [w[0], w[1] if len(w) > 1 else '']
            out.append(cur)
        elif cur is not None and ln.startswith(' ') and ln.strip() and ln.strip()[0] in '{}-.0123456789,' and \
                not ln.lstrip().startswith(('Matched', '$$$')):
            cur[1] += ' ' + ln
        else:
            cur = None
    res = []
    for t, r in out:
        r = r.split('}')[0]
        res.append((t, [float(x) for x in NUM.findall(STRIP.sub(' ', r.replace('{', '')))]))
    return res


def _run_calc_only(bench):
    from pysad.ffscmd import run_sad_string
    import io
    src = open(os.path.join(ROOT, 'bench_level%s.sad' % bench), encoding='utf-8').read()
    cut = len(src)
    for key in ('\nem = Emittance', '\nx0  ='):
        i = src.find(key)
        if i >= 0:
            cut = min(cut, i)
    src = src[:cut] + '\nstop;\n'
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


@pytest.mark.parametrize('bench,tags,tol', [
    ('25a', ('T0', 'T1', 'T2', 'GX', 'GY'), 5e-9),     # 垂直ステアリング (ANGLE=0, ROTATE=-90 DEG, K0)
    ('25b', ('T0', 'T1', 'T2'), 5e-9),                 # skew 四極
    ('25c', ('T0', 'GEO'), 5e-9),                      # ROTATE 付き BEND (5 DEG, 90 DEG, K1 付き)
    ('25d', ('T0',), 5e-8),                            # 六極・八極・十極・MULT の ROTATE
])
def test_calc_matches_sad_bench(bench, tags, tol):
    ref = os.path.join(ROOT, 'bench_level%s.result' % bench)
    if not os.path.exists(ref) or not os.path.exists(os.path.join(ROOT, 'bench_level%s.sad' % bench)):
        pytest.skip('ベンチファイルがない')
    a = _rows(open(ref, encoding='utf-8', errors='replace').read(), tags)
    b = _rows(_run_calc_only(bench), tags)
    assert len(a) == len(b) and len(a) > 5
    worst = 0.0
    for (ta, na), (tb, nb) in zip(a, b):
        assert ta == tb
        n = min(len(na), len(nb))
        assert n > 0
        worst = max(worst, max(abs(x - y) for x, y in zip(na[:n], nb[:n])))
    assert worst < tol
