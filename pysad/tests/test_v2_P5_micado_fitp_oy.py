# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 4・5・6: OYStability の後半、Micado の変種、FITP>1 を gs と比べる。

(項目 5) Micado: 12 セルの FODO (3 GeV、周長 52.8 m、BPM 24 台、水平・垂直ステアラ各 12 台、四極の据え付け誤差 DX・DY)。水平 (`"X"`)・2 平面同時 (`{"X","Y"}`)・
  `SVDmode->1` の垂直を n = 0・1・2・3・6・12 で、k64 が gs と 5e-6 以内 (rms・最大のキック。D2 の 1e-6 と同じ程度)。v2 は B129 の修正で、選んだ本数 n を増やすと rms が単調に下がり、
  n = 12 は最小二乗 (n = 0) と同じ。`MakeBump` は gs と 1e-12。`CorrectOrbit[..., Bump->b]` は K64 の誤り (CorrectOrbitBump に引数を誤って渡す) で gs も pysad も
  Index エラーになり、補正は何も起きない (K64 の誤りの候補 P5-5)。`CorrectOrbitBump` を直接呼ぶ使い方は gs が 300 s で終わらず未確認。
(項目 6) FITP>1: 12 セルの FODO、δ = ±2% の 5 点・±1% の 3 点 (NX・NY・BX・BY)。k64 は gs と最終の残差が 2e-4 以内・K1 が 6 桁一致・K2 が 1e-3 以内 (K2 は sextupole の谷が浅い)。
  反復の経路は違う (gs は 2〜3 反復で `Overflow & closed orbit not found` の試行点で止まり、pysad は 1〜2 反復多い)。FITP 1 は反復の残差まで一致。
(項目 4) OYStability の後半 (`CSROYMatrix` → `Eigensystem` → 成長率): gs は `CSRSetupOY` の出力を `CSROYMatrix` に渡すと wrongtype (B60) で OYStability が完走しないので、
  pysad v2 の `CSRSetupOY` の 6 つのリストを台本にそのまま書き (nr=256、nj=6、nphi=3)、gs の `CSROYMatrix`・`Eigensystem` と比べた。行列の和、固有値 (18 個)、成長率が 1e-12 以内
  (安定 aw=1e-8: 成長率 −damp、不安定 aw=1e-3: 0.0272796226)。**未確認**: `OYStability`・`OYScan` の通し (CSRStability の構築・ZL の用意を含む) は gs が動かないので gs とは比べていない。
"""
import re

import numpy as np
import pytest

from . import _p5_util as U


def _lines(text, prefix):
    return [l for l in text.splitlines() if l.startswith(prefix)]


def _micado_rows(text, prefix):
    rows = []
    for l in _lines(text, prefix):
        m = re.match(re.escape(prefix) + r'\s*n=(\d+)\s+(.*)', l)
        v = U.floats(m.group(2))
        rows.append((int(m.group(1)), v))
    return rows


@pytest.fixture(scope='module')
def micado():
    src = U.example_src('micado_variants')
    return {'gs': U.gs_text('micado_variants'), 'k64': U.run(src, 'k64'), 'v2': U.run(src, 'v2')}


@pytest.mark.parametrize('prefix', ['H Micado', 'XY(sh only)', 'SVD1 Y'])
def test_micado_variants_match_gs_k64(micado, prefix):
    g, p = _micado_rows(micado['gs'], prefix), _micado_rows(micado['k64'], prefix)
    assert len(g) == len(p) == 6
    for (n1, a), (n2, b) in zip(g, p):
        assert n1 == n2
        a, b = np.array(a), np.array(b)
        assert np.allclose(a, b, rtol=5e-6, atol=1e-12), (prefix, n1, a, b)


def test_micado_horizontal_v2_decreases_monotonically_k64_does_not(micado):
    """B129 (D2) の水平面: k64 は n = 1 で補正前 (6.85e-4) より悪化する (gs と同じ)。v2 は n を増やすと rms が単調に下がり n = 12 は最小二乗"""
    k = {n: v[0] for n, v in _micado_rows(micado['k64'], 'H Micado')}
    v = {n: v[0] for n, v in _micado_rows(micado['v2'], 'H Micado')}
    # rmsX after は 'rmsX after X' の最初の数 (行は `rmsX after <x> max kick ...`)
    assert k[1] > 5.0e-4 and k[2] > 5.0e-4                       # 補正前 6.85e-4 に近く、n = 0 の 7.3e-5 より 1 桁悪い
    vs = [v[n] for n in (1, 2, 3, 6, 12)]
    assert all(a > b for a, b in zip(vs, vs[1:]))
    assert v[12] == pytest.approx(v[0], rel=1e-6)
    assert v[0] == pytest.approx(k[0], rel=1e-6)                 # n = 0 (最小二乗) は両モード同じ


def test_make_bump_matches_gs_and_corrector_option_is_broken_in_both():
    src = U.example_src('micado_bump')
    gs, k = U.gs_text('micado_bump'), U.run(src, 'k64')
    kg = U.floats([l for l in gs.splitlines() if l.startswith('MakeBump')][0].split('kicks', 1)[1])
    kp = U.floats([l for l in k.splitlines() if l.startswith('MakeBump')][0].split('kicks', 1)[1])
    assert np.allclose(kg, kp, rtol=1e-12, atol=1e-18)
    # CorrectOrbit の Bump-> は両方とも誤り (r は評価されず、rms は補正前のまま)
    for t in (gs, k):
        assert 'index' in t.lower() and 'out of range' in t.lower()
        assert re.search(r'^Bump-> result r\s*$', t, re.M)
    rb = U.floats([l for l in gs.splitlines() if l.startswith('rms before Y')][0].split('Y', 1)[1])[0]
    ra = U.floats([l for l in gs.splitlines() if l.startswith('after rms Y')][0].split('Y', 1)[1])[0]
    assert ra == pytest.approx(rb, rel=1e-9)


@pytest.mark.parametrize('name', ['five', 'beta', 'one'])
def test_fitp_final_solution_matches_gs_k64(name):
    src = U.example_src('fitp_' + name)
    gs = U.gs_text('fitp_' + name)
    out = {m: U.run(src, m) for m in ('k64', 'v2')}

    def final(t):
        r = [l for l in t.splitlines() if l.strip().startswith('Residual =')]
        res = U.floats(r[-1].split('DP')[0].split('=')[1])[0] if r else 0.0
        a = [l for l in t.splitlines() if l.startswith('after GO')][0]
        return res, U.floats(a.split('QD =')[1].replace('K2 SF SD =', ' '))
    rg, kg = final(gs)
    for m in ('k64', 'v2'):
        rp, kp = final(out[m])
        if name != 'one':
            assert rp == pytest.approx(rg, rel=2e-4), (m, rp, rg)
        assert kp[0] == pytest.approx(kg[0], abs=2e-6) and kp[1] == pytest.approx(kg[1], abs=2e-6)       # K1 QD
        if name != 'one':
            assert np.allclose(kp[2:], kg[2:], rtol=1e-3), (m, kp, kg)                                       # K2 SF SD
    assert final(out['k64']) == final(out['v2'])                       # FITP は両モードで同じ


def test_fitp_first_iteration_residual_matches_gs_for_fitp_1():
    src = U.example_src('fitp_one')
    gs, k = U.gs_text('fitp_one'), U.run(src, 'k64')
    g = U.floats([l for l in gs.splitlines() if re.match(r'^\s+2\s+0\s', l)][0].split('(')[0])[2]
    p = U.floats([l for l in k.splitlines() if re.match(r'^\s+2\s+0\s', l)][0].split('(')[0])[2]
    assert p == pytest.approx(g, rel=1e-4)


@pytest.mark.parametrize('name', ['oy_pipeline', 'oy_pipeline_unstable'])
def test_oy_pipeline_matches_gs(name):
    src = U.example_src(name)
    gs = U.gs_text(name)
    for mode in ('k64', 'v2'):
        out = U.run(src, mode)
        for tag in ('MSUM', 'GR'):
            g = U.floats([l for l in gs.splitlines() if l.startswith(tag + ' ')][0][len(tag):])
            p = U.floats([l for l in out.splitlines() if l.startswith(tag + ' ')][0][len(tag):])
            assert np.allclose(g, p, rtol=1e-12), (mode, tag, g, p)
        for tag in ('EVRE', 'EVIM'):
            g = U.floats([l for l in gs.splitlines() if l.startswith(tag + ' ')][0][len(tag):])
            p = U.floats([l for l in out.splitlines() if l.startswith(tag + ' ')][0][len(tag):])
            assert len(g) == len(p) == 18 and np.allclose(g, p, rtol=1e-11, atol=1e-13)


def test_oy_growth_rate_values():
    """安定 (aw = 1e-8) の成長率は −damp (減衰率 1e-4)、不安定 (aw = 1e-3) は +0.0272796226 (gs と同じ)"""
    for name, want in (('oy_pipeline', -1e-4), ('oy_pipeline_unstable', 0.027279622601971)):
        out = U.run(U.example_src(name), 'v2')
        g = U.floats([l for l in out.splitlines() if l.startswith('GR ')][0][3:])[0]
        assert g == pytest.approx(want, rel=1e-9)
