# -*- coding: utf-8 -*-
"""KEK-LS v03_68 の例 (example/kek_ls): ラティスの `cell calc emit` の値 (周長・チューン・εx) と、速いスクリプトの数本が pysad で走ること。
gs の保存値 (example/kek_ls/*/NAME.gs.txt) とも k64 で一致する。スクリプトは `Get["../modules/..."]` を相対パスで読むので、
台本のあるフォルダをカレントにして走らせる (cwd 引数)。"""
import io
import os
import re

import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), 'example', 'kek_ls')


def _run(rel):
    path = os.path.join(EX, rel)
    src = open(path, encoding='utf-8').read()
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf, cwd=os.path.dirname(path))
    return buf.getvalue()


def _nums(txt, pat):
    # gs の出力には台本の復唱行も入るので、数値に読める最初の一致を使う
    for m in re.finditer(pat, txt):
        try:
            return [float(x) for x in m.groups()]
        except ValueError:
            continue
    raise AssertionError(pat)


def _saved(rel):
    p = os.path.join(EX, rel)
    if not os.path.exists(p):
        pytest.skip('gs output not saved')
    return open(p, encoding='utf-8', errors='replace').read()


def test_lattice_header_and_original_body():
    txt = open(os.path.join(EX, 'lattice', 'kek_ls_v03_68.sad'), encoding='utf-8').read()
    assert 'HMBA' in txt and 'Pantaleo Raimondi' in txt and 'Simone Liuzzo' in txt
    assert '\r' not in txt
    assert 'C:\\' not in txt and 'D:\\' not in txt


def test_cell_calc_emit_k64():
    compat.set_mode('k64')
    txt = _run('lattice/KEKLS_cell_calc_emit.sad')
    # RINGC (標準セル 20 個) の後に RING (入射セル込み)
    parts = txt.split('LINE RING\n')
    assert len(parts) == 2
    c, r = parts
    assert _nums(c, r'circumference \[m\]\s*=\s*(\S+)')[0] == pytest.approx(570.72, abs=1e-6)
    assert _nums(c, r'tune nx ny\s*=\s*(\S+)\s+(\S+)') == pytest.approx([48.5800234, 17.6200786], abs=1e-6)
    assert _nums(c, r'EMITX \[m rad\]\s*=\s*(\S+)')[0] == pytest.approx(1.343162e-10, rel=1e-5)
    assert _nums(c, r'SIGE\s*=\s*(\S+)')[0] == pytest.approx(6.405158e-4, rel=1e-5)
    assert _nums(r, r'circumference \[m\]\s*=\s*(\S+)')[0] == pytest.approx(570.7209201, abs=1e-6)
    assert _nums(r, r'EMITX \[m rad\]\s*=\s*(\S+)')[0] == pytest.approx(1.356809e-10, rel=1e-5)


def test_cell_calc_emit_matches_gs_saved():
    compat.set_mode('k64')
    g = _saved('lattice/KEKLS_cell_calc_emit.gs.txt')
    p = _run('lattice/KEKLS_cell_calc_emit.sad')
    for pat in (r'EMITX \[m rad\]\s*=\s*(\S+)', r'tune nx ny\s*=\s*(\S+)\s+(\S+)', r'circumference \[m\]\s*=\s*(\S+)'):
        assert _nums(p, pat) == pytest.approx(_nums(g, pat), rel=1e-6)


def test_match_ring_geo_returns_to_design_drift():
    compat.set_mode('k64')
    txt = _run('design/KEKLS_match_ring_geo.sad')
    # 標準の 2.94 m から始めて、設計の DRID1INJ = 2.76546005 m に戻る
    assert _nums(txt, r'start\s*:\s*DRID1INJ\s+(\S+)')[0] == pytest.approx(2.94)
    assert _nums(txt, r'result\s*:\s*DRID1INJ\s+(\S+)')[0] == pytest.approx(2.76546005, abs=1e-6)


def test_match_injcell_converges():
    compat.set_mode('k64')
    txt = _run('design/KEKLS_match_injcell.sad')
    last = [l for l in txt.split('\n') if l.startswith('iteration')][-1]
    nums = [float(x) for x in re.findall(r'[-+]?(?:\d+\.\d*|\.\d+)(?:E[-+]?\d+)?', last)]
    # 位相進み (1.2145, 0.4405) と BX(MKSI) = 30 m
    assert nums[1] == pytest.approx(1.2145007, abs=1e-5) and nums[2] == pytest.approx(0.4405045, abs=1e-5)
    assert nums[-1] == pytest.approx(30.0, abs=1e-3)


def test_calc_H_matches_emitx_within_few_percent():
    compat.set_mode('k64')
    txt = _run('display/KEKLS_calc_H.sad')
    ratio = _nums(txt, r'ratio\s*=\s*(\S+)')[0]
    assert ratio == pytest.approx(0.9767, abs=1e-3)
    assert _nums(txt, r'Jx\s*=\s*(\S+)')[0] == pytest.approx(1.3110389, abs=1e-6)


def test_chromaticity_has_two_tables():
    compat.set_mode('k64')
    txt = _run('display/KEKLS_chromaticity.sad')
    t1 = [l for l in txt.split('\n') if l.startswith('T1 ')]
    t0 = [l for l in txt.split('\n') if l.startswith('T0 ')]
    assert len(t1) == 13 and len(t0) == 13
    # δ = 0 のチューンは設計値
    z = [l for l in t1 if abs(float(l.split()[1])) < 1e-12][0].split()
    assert float(z[2]) == pytest.approx(48.5800234, abs=1e-6)


def test_cod_response_matches_sad_cod():
    compat.set_mode('k64')
    txt = _run('errors/KEKLS_cod_response.sad')
    errs = [float(m) for m in re.findall(r'max \|SAD \+ K0\*Response\| \[um\]\s*(\S+)', txt)]
    amps = [float(m) for m in re.findall(r'max \|COD\| \[um\]\s*(\S+)', txt)]
    assert len(errs) == 3
    # 解析式と SAD の COD が 1e-3 以内 (最大振幅に対して)
    assert all(e < 1e-3 * a + 1e-3 for e, a in zip(errs, amps))
