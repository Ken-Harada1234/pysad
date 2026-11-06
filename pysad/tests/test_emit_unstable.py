# -*- coding: utf-8 -*-
"""不安定な光学の EMIT / INS + CALC6D / RFSW ON + RADCOD OFF の tcod (docs/history/work_emit_unstable.md) の回帰テスト。

A. teigen.py (teigen.f の直訳) が固有値分解として正しい (A v = λ v)、tnorm が実固有値の対も扱える
   + 不安定な ESRF (bench_level28c) の EMIT 報告が SAD と行単位で一致 (丸め雑音の桁を除く)
B. INS + CALC6D (bench_level30a: 先頭 MARK の 28 成分から 6D 伝搬) が SAD と 1e-9
C. RFSW ON + RADCOD OFF + FSHIFT≠0 の Emittance[] (bench_level30b: tcod の 6D Newton の直訳) が SAD と 1e-10
"""
import io
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

from pysad import ffscmd
from pysad.teigen import teigen
from pysad import emit6d


def _run(name):
    src = open(os.path.join(ROOT, name), encoding='utf-8').read()
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _have(*names):
    return all(os.path.exists(os.path.join(ROOT, n)) for n in names)


def _symplectic(rng, scale=0.3):
    from scipy.linalg import expm
    J = np.zeros((6, 6))
    for i in range(0, 6, 2):
        J[i, i + 1] = 1.0
        J[i + 1, i] = -1.0
    M = np.eye(6)
    for _ in range(4):
        S = rng.normal(size=(6, 6))
        S = (S + S.T) * scale
        M = M @ expm(J @ S)
    return M


def test_teigen_is_an_eigen_decomposition():
    rng = np.random.default_rng(1)
    for _ in range(6):
        M = _symplectic(rng)
        w, V = teigen(M)
        i = 0
        while i < 6:
            if w[i].imag != 0.0:
                v = V[:, i] + 1j * V[:, i + 1]       # 複素固有値: 列 i = 実部, i+1 = 虚部
                assert np.abs(M @ v - w[i] * v).max() < 1e-9 * np.abs(v).max() * max(1.0, abs(w[i]))
                i += 2
            else:
                assert np.abs(M @ V[:, i] - w[i].real * V[:, i]).max() < 1e-9 * np.abs(V[:, i]).max() * max(1.0, abs(w[i]))
                i += 1


def test_tnorm_handles_a_real_pair():
    """y 面が不安定 (実固有値 λ, 1/λ) でも r が (ほぼ) シンプレクティックで、ceig の虚部は 0 の対になる"""
    M = np.eye(6)
    c, s = np.cos(1.3), np.sin(1.3)
    M[0:2, 0:2] = [[c, 5 * s], [-s / 5, c]]
    ch, sh = np.cosh(1.0), np.sinh(1.0)
    M[2:4, 2:4] = [[ch, 4 * sh], [sh / 4, ch]]
    M[4:6, 4:6] = [[np.cos(0.05), 80 * np.sin(0.05)], [-np.sin(0.05) / 80, np.cos(0.05)]]
    r, ceig = emit6d._normal_form(M)
    J = np.zeros((6, 6))
    for i in range(0, 6, 2):
        J[i, i + 1] = 1.0
        J[i + 1, i] = -1.0
    assert np.abs(r.T @ J @ r - J).max() < 1e-9 or np.abs(emit6d._tsymp(r).T @ J @ emit6d._tsymp(r) - J).max() < 1e-9
    sort = sorted(ceig.real[2:4])
    assert ceig[2].imag == 0.0 and ceig[3].imag == 0.0
    assert abs(sort[0] * sort[1] - 1.0) < 1e-9


@pytest.mark.skipif(not _have('bench_level28c.sad', 'bench_level28c.result'), reason='bench ファイルなし')
def test_unstable_emit_report_matches_sad():
    import cmp_emit_report as cer
    ffs, out = _run('bench_level28c.sad')
    a = cer.load_blocks(os.path.join(ROOT, 'bench_level28c.result'))
    b = cer.blocks(out.split('\n'))
    assert len(a) == len(b) == 1
    assert len(a[0]) == len(b[0])
    bad = [(x, y) for x, y in zip(a[0], b[0]) if not cer.same(x, y, 1e-10, 1e-9)[0]]
    assert bad == []


@pytest.mark.skipif(not _have('bench_level30a.sad', 'bench_level30a.result'), reason='bench ファイルなし')
def test_ins_calc6d_matches_sad():
    import cmp_lines
    ffs, out = _run('bench_level30a.sad')
    a = cmp_lines.rows(open(os.path.join(ROOT, 'bench_level30a.result'), encoding='utf-8', errors='replace'))
    b = cmp_lines.rows(out.split('\n'))
    res = cmp_lines.compare(a, b)
    assert len(res) >= 30
    assert all(r[1] is not None and r[1] <= 1e-9 for r in res)


@pytest.mark.skipif(not _have('bench_level30b.sad', 'bench_level30b.result'), reason='bench ファイルなし')
def test_rfsw_noradcod_fshift_tcod_matches_sad():
    import cmp_lines
    ffs, out = _run('bench_level30b.sad')
    a = cmp_lines.rows(open(os.path.join(ROOT, 'bench_level30b.result'), encoding='utf-8', errors='replace'))
    b = cmp_lines.rows(out.split('\n'))
    res = cmp_lines.compare(a, b)
    assert len(res) >= 30
    assert all(r[1] is not None and r[1] <= 1e-10 for r in res)
