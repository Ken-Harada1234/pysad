"""特殊関数・数学関数 (specfn.py と sadscript の基本関数)。参照は SAD K64 (gs) の出力 (specfn_ref.json: 式と数値のリスト)。
gs が誤値を返す式 (docs/21 B43〜B47) と、gs が応答しない式は入れていない。"""
import io
import json
import os
import re

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _nums(s):
    s = s.replace('-INF', '-1e999').replace('INF', '1e999')
    return [float(t) for t in re.findall(_NUM, s)]


def test_special_functions_match_gs():
    from pysad.ffscmd import run_sad_string
    ref = json.load(open(os.path.join(os.path.dirname(__file__), 'specfn_ref.json'), encoding='utf-8'))
    src = 'FFS;\n$FORM="18.12";\n' + ''.join('Print["R%d= ",%s];\n' % (i, e) for i, (e, _) in enumerate(ref)) + 'Exit;\n'
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    got = {}
    for l in buf.getvalue().splitlines():
        m = re.match(r'R(\d+)= (.*)', l)
        if m:
            got[int(m.group(1))] = m.group(2)
    bad = []
    for i, (e, gn) in enumerate(ref):
        pn = _nums(got.get(i, ''))
        if len(pn) != len(gn) or any(abs(a - b) > 1e-8 * max(abs(a), 1e-3) + 1e-12 for a, b in zip(gn, pn)):
            bad.append((e, gn, got.get(i)))
    assert not bad, bad[:5]


def test_k64_inverseerf_is_reproduced():
    from pysad.specfn import inverseerf
    assert abs(inverseerf(.3) - .272462711811) < 1e-12      # 真値は .272462714727 (K64 の反復は相対 1e-8 で止まる)
    assert abs(inverseerf(.999) - 2.326753739975) < 1e-11
