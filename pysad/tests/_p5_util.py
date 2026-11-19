# -*- coding: utf-8 -*-
"""v2next 担当 P5 の試験の補助: 例 (example/v2_next/P5_*) の台本・gs の出力を読み、pysad を指定の互換モードで走らせる。"""
import io
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EX = os.path.join(ROOT, 'example', 'v2_next')
LAT = {
    'A': os.path.join(HERE, 'kek-ls_v03_12.lat'),
    'PF': os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'),
    'AST': os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'),
}
NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def floats(s):
    return [float(x) for x in re.findall(NUM, s.replace('INF', '1e999'))]


def need(path):
    if not os.path.exists(path):
        pytest.skip(os.path.basename(path) + ' がない')
    return path


def example_src(name, lat=None):
    """example/v2_next/P5_<name>.sad (lat があれば、そのラティスを先頭に連結。gs に渡した台本と同じ)"""
    src = open(need(os.path.join(EX, 'P5_%s.sad' % name)), encoding='utf-8').read()
    if lat:
        src = open(need(LAT[lat]), encoding='utf-8', errors='replace').read() + '\n' + src
    return src


def gs_text(name):
    return open(need(os.path.join(EX, 'P5_%s.gs.txt' % name)), encoding='utf-8', errors='replace').read()


def run(src, mode):
    """pysad を互換モード mode ('k64' か 'v2') で走らせた画面表示"""
    from pysad import compat
    from pysad.ffscmd import run_sad_string
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        run_sad_string(src, out=buf)
        return buf.getvalue()
    finally:
        compat.set_mode(old)


def stat_blocks(text):
    """MEASURE・TRCOD の統計の表示 (' Statistics at' から ' sigp/p:' の行まで) のリスト。各ブロックは行のリスト"""
    lines = text.splitlines()
    out, cur = [], None
    for l in lines:
        if l.startswith(' Statistics at'):
            cur = [l.rstrip()]
        elif cur is not None:
            cur.append(l.rstrip())
            if l.startswith(' sigp/p:'):
                out.append(cur)
                cur = None
    return out


def blocks_close(a, b, rtol=2e-3, atol=1e-19):
    """2 つの統計ブロックが、表示の 4 桁で一致するか (丸め誤差だけの 1e-20 の値は atol で同じとみなす)"""
    import numpy as np
    if len(a) != len(b):
        return False
    for la, lb in zip(a, b):
        if la.startswith(' Statistics') or la.startswith(' RAD'):
            if la != lb:
                return False
            continue
        fa, fb = floats(la.split(':', 1)[-1]), floats(lb.split(':', 1)[-1])
        if len(fa) != len(fb) or not np.allclose(fa, fb, rtol=rtol, atol=atol):
            return False
    return True
