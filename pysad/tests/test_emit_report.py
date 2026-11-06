# -*- coding: utf-8 -*-
"""EMIT 報告の全行 (temitprint) と、四極・薄い多極の解析行列 (tquade/tqlfre/tsolrote/tthine) のテスト (docs/history/work_emit_report.md)。

(1) 四極 (ずれ・ロール・線形フリンジ FRINGE=1,2,3・K1<0) の解析ヤコビアン `emit6d._jac_quad_full` を、COD 写像の数値微分と比べる
    (z の行 (列 x,px,y,py) は SAD がシンプレクティック条件で作り直すので比べない。(5,6) は SAD の式で数値微分とは 1e-7 違う)
(2) 薄い多極 (SEXT) のずれ・ロール付きの `lattice._thin_cod_matrix` が写像の数値微分と一致する
(3) SAD K64 の EMIT 報告 (bench_level26a.result: ESRF の RFSW × RADCOD の 4 通り) と pysad の報告を行単位で突き合わせる
    (0 同然の雑音の桁 (|値|<1e-12) の差だけは許す)
(4) SAD の FSHIFT (Orbit dilation) を Emittance の閉軌道に入れる (bench_level26a と同じ型の 1 ケース)
"""
import math
import os
import re
import subprocess
import sys

import numpy as np
import pytest

from pysad import lattice as L
from pysad import emit6d as E

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PH = E._Phys(3e9, 844.0)


def _num(f, c, eps=1e-6):
    return L._num_jacobian(f, np.asarray(c, float), eps)


def _cods(n=3, seed=3):
    rng = np.random.default_rng(seed)
    for _ in range(n):
        c = np.array([rng.normal() * 3e-3, rng.normal() * 3e-4, rng.normal() * 3e-3, rng.normal() * 3e-4, 0.0, 4e-3])
        yield c


@pytest.mark.parametrize('K1', [1.2, -0.9])
@pytest.mark.parametrize('mfring', [0, 1, 2, 3])
@pytest.mark.parametrize('off', [(0, 0, 0), (2e-4, -3e-4, 0), (0, 0, 8e-4), (1e-3, 5e-4, 2e-3)])
def test_quad_jacobian_offset_roll_fringe(K1, mfring, off):
    dx, dy, th = off
    fr = (0.05, 0.02, 0.06, 0.03, mfring) if mfring else (0.0, 0.0, 0.0, 0.0, 0)
    for fringe in (True, False):
        for c in _cods():
            f = E._quad_map(0.4, K1, PH, dx, dy, th, fringe, fr)
            Ja = E._jac_quad_full(c, 0.4, K1, fringe, fr, PH, dx, dy, th)
            Jn = _num(f, c)
            d = np.abs(Ja - Jn)
            d[4, :] = 0.0                      # z の行: SAD の式 (数値微分と違う)
            assert d.max() < 1e-8, (K1, mfring, off, fringe, d.max())


def test_tqlfre_jacobian():
    for f1, f2 in ((0.05, 0.01), (-0.05, 0.02)):
        for c in _cods():
            f = lambda cc, f1=f1, f2=f2: L._tqlfre_cod(cc, 0.3, 0.54, f1, f2)
            Ja = E._jac_tqlfre(c, 0.3, 0.54, f1, f2)
            Jn = _num(f, c)
            d = np.abs(Ja - Jn)
            d[4, :] = 0.0
            assert d.max() < 1e-9


@pytest.mark.parametrize('kord,ak', [(2, 0.6), (3, 40.0)])
def test_thin_multipole_jacobian_with_offset(kord, ak):
    for L_ in (0.0, 0.2):
        for off in ((0.0, 0.0, 0.0), (5e-4, -4e-4, 0.0), (0.0, 0.0, 2e-3)):
            dx, dy, th = off
            for c in _cods():
                cod_out, M = L._thin_cod_matrix(L_, ak, kord, dx, dy, th, c)
                f = lambda cc: L._tchge_cod(L._thin_offset_map(L_, ak, kord)(
                    L._tchge_cod(cc, dx, dy, 0.0, th, 0.0, 0.0, 0.0, True)), dx, dy, 0.0, th, 0.0, 0.0, 0.0, False)
                assert np.max(np.abs(f(c) - cod_out)) < 1e-14
                Jn = _num(f, c)
                d = np.abs(M - Jn)
                d[4, :] = 0.0
                if th != 0.0 and (dx != 0.0 or dy != 0.0):
                    continue                    # SAD の癖: 入口の行列に回転が入らない (lattice._tchge_jac)
                assert d.max() < 1e-8, (kord, L_, off, d.max())


# ---------------------------------------------------------------------------
# EMIT 報告 (bench_level26a.result との行単位の比較)
# ---------------------------------------------------------------------------
NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def _blocks(lines):
    out, cur = [], None
    for l in lines:
        l = l.rstrip('\n')
        if l.strip() == 'Closed orbit:':
            cur = [l]
            continue
        if cur is not None:
            cur.append(l)
            if l.startswith('Nominal spin tune'):
                out.append(cur)
                cur = None
    return out


def _same(lx, ly):
    if lx == ly:
        return True
    nx, ny = NUM.findall(lx), NUM.findall(ly)
    # 数値以外の並びは空白を全部除いて比べる ('tune:-0.0000000' と 'tune: 0.0000000' のように負の 0 で空白の位置が変わる)
    if len(nx) != len(ny) or ''.join(NUM.sub('#', lx).split()) != ''.join(NUM.sub('#', ly).split()):
        return False
    # 0 同然の量 (Beam tilt など、どちらも |v|<1e-9) は丸め雑音なので同じとみなす (方針: 1e-10 以下の差は追わない)
    return all(abs(float(p) - float(q)) <= 1e-11 or abs(float(p) - float(q)) <= 1e-9 * max(abs(float(p)), abs(float(q)))
               or (abs(float(p)) < 1e-9 and abs(float(q)) < 1e-9)
               for p, q in zip(nx, ny))


def _run(args, timeout=600):
    env = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1')
    r = subprocess.run([sys.executable, '-u', os.path.join(ROOT, 'run_sad.py')] + args, cwd=ROOT, env=env,
                       capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=timeout)
    return r.stdout.splitlines()


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, 'bench_level26a.result')), reason='SAD 参照結果なし')
def test_emit_report_lines_vs_sad_4_modes():
    ref = _blocks(open(os.path.join(ROOT, 'bench_level26a.result'), encoding='utf-8', errors='replace').readlines())
    got = _blocks(_run(['bench_level26a.sad']))
    assert len(ref) == len(got) == 4
    bad = []
    for k, (x, y) in enumerate(zip(ref, got)):
        assert len(x) == len(y), (k, len(x), len(y))
        for i, (lx, ly) in enumerate(zip(x, y)):
            if not _same(lx, ly):
                bad.append((k + 1, i, lx, ly))
    assert not bad, bad[:3]


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, 'bench_level26c.result')), reason='SAD 参照結果なし')
def test_emit_report_fshift():
    ref = _blocks(open(os.path.join(ROOT, 'bench_level26c.result'), encoding='utf-8', errors='replace').readlines())
    got = _blocks(_run(['bench_level26c.sad']))
    assert len(ref) == len(got) == 1
    bad = [(i, a, b) for i, (a, b) in enumerate(zip(ref[0], got[0])) if not _same(a, b)]
    assert not bad, bad[:3]
