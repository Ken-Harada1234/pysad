# -*- coding: utf-8 -*-
"""MEASURE・DETR・FITP・Micado・OYStability・陽子の β・複素特殊関数の走査の例 (example の機能別フォルダ) の台本・gs の出力を読み、pysad を指定の互換モードで走らせる補助。"""
import io
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EX = os.path.join(ROOT, 'example')
# 例の名前 → example からの相対パス (拡張子なし。.sad と .gs.txt が並ぶ)
EXAMPLES = {
    'detr_astrid2': '02_optics_cod/flag_DETR_astrid2',
    'detr_pf': '02_optics_cod/flag_DETR_pf',
    'detr_ring_a': '02_optics_cod/flag_DETR_ring_a',
    'fitp_beta': '03_matching/flag_FITP_beta',
    'fitp_five': '03_matching/flag_FITP_five',
    'fitp_one': '03_matching/flag_FITP_one',
    'measure_astrid2': '02_optics_cod/fn_MEASURE_astrid2',
    'measure_default_np': '02_optics_cod/fn_MEASURE_default_np',
    'measure_fodo': '02_optics_cod/fn_MEASURE_fodo',
    'measure_pf': '02_optics_cod/fn_MEASURE_pf',
    'measure_ring_a': '02_optics_cod/fn_MEASURE_ring_a',
    'micado_bump': '13_packages/fn_CorrectOrbit_Bump',
    'micado_variants': '13_packages/fn_CorrectOrbit_Micado_variants',
    'oy_pipeline': '13_packages/fn_OYStability_pipeline',
    'oy_pipeline_unstable': '13_packages/fn_OYStability_pipeline_unstable',
    'proton_beta_p1p5': '04_emittance_radiation/var_MASS_proton_beta_p1p5',
    'proton_beta_p3': '04_emittance_radiation/var_MASS_proton_beta_p3',
    'proton_beta_p30': '04_emittance_radiation/var_MASS_proton_beta_p30',
    'specfn_cx_scan': '10_math_data/fn_SpecialFunctions_complex_scan',
}
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
    """EXAMPLES[name] の .sad (lat があれば、そのラティスを先頭に連結。gs に渡した台本と同じ)"""
    src = open(need(os.path.join(EX, EXAMPLES[name] + '.sad')), encoding='utf-8').read()
    if lat:
        src = open(need(LAT[lat]), encoding='utf-8', errors='replace').read() + '\n' + src
    return src


def gs_text(name):
    return open(need(os.path.join(EX, EXAMPLES[name] + '.gs.txt')), encoding='utf-8', errors='replace').read()


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
