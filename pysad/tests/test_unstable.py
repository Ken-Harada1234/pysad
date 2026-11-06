# -*- coding: utf-8 -*-
"""不安定ラインでの CELL CALC の続行と nstab の判定 (docs/history/work_unstable.md) のテスト。

(1) qmdiag (全分岐の直訳) の通常の場合が従来の periodic_r と同じ R を返す
(2) 安定な周期系では SAD 流の経路 (sad_stab=) が従来の Twiss と同じ値を返し、`_qstat` が安定を報告する
(3) |cosμ|>1 の周期系: cosμ → 1/cosμ に置き換えた値で続行し、INTRES で stabx/staby が偽、NOINTRES で安定扱い (tracex=4/trace)
(4) 一周行列の成分が 1e10 を超えると ±1e10 に抑える (qcod の limitnan)
(5) fmincube の上限はオイラー定数 γ、log10(0) は -Infinity (matching02_r / tune_sv の反復を決めた 2 つの定数・関数)
(6) SAD の bench_level28a.result (ESRF の QD ×1.19: 垂直だけ不安定) の Twiss と `Vertical Unstable.` が run_sad で再現できる
"""
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from pysad import lattice as L
from pysad import match6 as M
from pysad.elements import make_element
from pysad.optics import periodic_r

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))


def _fodo(kf=1.2, kd=-1.2, n=2):
    els = []
    for _ in range(n):
        els += [make_element('QUAD', 'QF', L=0.3, K1=kf), make_element('DRIFT', 'D', L=2.0),
                make_element('QUAD', 'QD', L=0.3, K1=kd), make_element('DRIFT', 'D', L=2.0)]
    return L.Lattice(els)


def _stable_sym(seed):
    rng = np.random.default_rng(seed)
    # 2 つの回転 (非結合) を結合の小さな相似変換で混ぜる: 安定で結合した 4x4 一周行列
    def rot(mu, beta, alpha):
        c, s = math.cos(mu), math.sin(mu)
        return np.array([[c + alpha * s, beta * s], [-(1 + alpha**2) / beta * s, c - alpha * s]])
    m = np.zeros((4, 4))
    m[:2, :2] = rot(1.1 + rng.random(), 5 * rng.random() + 1, rng.normal())
    m[2:, 2:] = rot(0.7 + rng.random(), 5 * rng.random() + 1, rng.normal())
    return m


def test_qmdiag_normal_matches_periodic_r_for_uncoupled():
    m = _stable_sym(1)
    r1, r2, r3, r4, c, stab, nanq, msg = L._qmdiag(m)
    assert stab and not nanq
    # 非結合の行列では R = 0、c=1 (変換は単位)
    assert max(abs(r1), abs(r2), abs(r3), abs(r4)) < 1e-12
    assert abs(c - 1.0) < 1e-12


def test_qmdiag_coupled_agrees_with_periodic_r():
    # 要素の積で作る結合した一周行列 (回転四極を含む)
    from pysad.elements import make_element as mk
    els = [mk('QUAD', 'QF', L=0.3, K1=1.1, ROTATE=0.05), mk('DRIFT', 'D', L=2.0), mk('QUAD', 'QD', L=0.3, K1=-1.0),
           mk('DRIFT', 'D', L=2.0)] * 2
    lat = L.Lattice(els)
    M6 = lat.transfer_matrix(0.0)
    R, A, B = periodic_r(M6[:4, :4])
    r1, r2, r3, r4, c, stab, nanq, msg = L._qmdiag(M6)
    assert stab
    assert np.allclose([r1, r2, r3, r4], [R[0, 0], R[0, 1], R[1, 0], R[1, 1]], atol=1e-10, rtol=1e-8)


def test_sad_stab_stable_ring_equals_legacy():
    lat1 = _fodo(1.1, -1.0, n=4)
    lat2 = _fodo(1.1, -1.0, n=4)
    ref = lat1.twiss_table()
    got = lat2.twiss_table(sad_stab={'intres': True, 'halfres': True})
    q = lat2._qstat
    assert q['stabx'] and q['staby'] and not q['over'] and q['codfnd']
    for a, b in zip(ref, got):
        assert a.betx == pytest.approx(b.betx, rel=1e-13)
        assert a.mux == pytest.approx(b.mux, rel=1e-13, abs=1e-13)
        assert a.etx == pytest.approx(b.etx, rel=1e-12, abs=1e-14)


def test_unstable_ring_continues_and_flags():
    # QF をきつくして 1 周のトレースが 2 を超える (cosμ>1)。従来は ValueError / RuntimeError
    lat = _fodo(2.4, -1.0, n=1)
    with pytest.raises((ValueError, RuntimeError)):
        lat.twiss_table()
    lat = _fodo(2.4, -1.0, n=1)
    rows = lat.twiss_table(sad_stab={'intres': True, 'halfres': True})
    q = lat._qstat
    assert not (q['stabx'] and q['staby'])
    assert all(math.isfinite(r.betx) and math.isfinite(r.bety) for r in rows)
    # NOINTRES / NOHALFRES では安定扱い (trace は 4/trace に置き換え)
    lat = _fodo(2.4, -1.0, n=1)
    lat.twiss_table(sad_stab={'intres': False, 'halfres': False})
    q2 = lat._qstat
    assert q2['stabx'] and q2['staby']
    assert abs(q2['tracex']) <= 2.0 and abs(q2['tracey']) <= 2.0


def test_trmax_clip_of_one_turn_matrix():
    lat = _fodo(2.4, -1.0, n=30)
    lat.twiss_table(sad_stab={'intres': True, 'halfres': True})
    assert L._LAST_MMAX > 1e10                       # クリップ前は 1e10 を超えている
    M6 = L._qcod_periodic(lat.elements, 0.0)[1]
    assert np.max(np.abs(M6[:4, :4])) <= 1e10 * (1 + 1e-12)


def test_fmincube_cap_is_euler_gamma_and_log10_zero():
    assert M.M_EULER == pytest.approx(0.5772156649015329, abs=1e-15)
    # f2=0 の枝: -0.5 f1 dg / ((g1-g0)/f1 - dg) = 0.0554 は γ·f1=0.0425 に抑えられる (matching02_r の 8 反復目)
    f = M._fmincube(7.3686067752573053e-02, 0.0, 4.068739537864065, 4.279836004676336, 4.279836004676336,
                    -8.559672009352672)
    assert f == pytest.approx(0.5772156649015329 * 7.3686067752573053e-02, rel=1e-14)
    assert M._log10(0.0) == -math.inf
    assert math.isnan(M._log10(-1.0))
    assert M._fuzz(M._log10(0.0), -3.5, -2.5) == 0.0     # badcnv = 1 - 0 = 1 → モード切替


def test_esrf_vertical_unstable_vs_sad(tmp_path):
    lat = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kek-ls_v03_12.lat'), encoding='utf-8').read()
    script = '''FFS;
$FORM="15.10"; PageWidth=1999;
USE RING;
CELL;
CALC;
k0={LINE["K1","QDA2"],LINE["K1","QDA3"],LINE["K1","QDA5"]};
LINE["K1","QDA2"]=k0[1]*1.19; LINE["K1","QDA3"]=k0[2]*1.19; LINE["K1","QDA5"]=k0[3]*1.19;
CALC;
Print["TWU_V: ", Twiss[{"NX","NY","BX","BY","AX","AY","EX","DX"},"$$$"]];
stop;
stop;
'''
    f = tmp_path / 'esrf_unstable.sad'
    f.write_text(lat + '\n' + script, encoding='utf-8')
    env = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1')
    out = subprocess.run([sys.executable, '-u', os.path.join(ROOT, 'run_sad.py'), str(f)], cwd=ROOT, env=env,
                         capture_output=True, text=True, encoding='utf-8', timeout=300).stdout
    assert ' Vertical Unstable.' in out
    assert 'MKS.1/$$$   BMAGY' in out
    sad = open(os.path.join(ROOT, 'bench_level28a.result'), encoding='utf-8', errors='replace').read()
    import re
    num = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
    pick = lambda t: [float(x) for x in num.findall(next(l for l in t.split('\n') if l.startswith('TWU_V:')))[:7]]
    a, b = pick(sad), pick(out)
    for x, y in zip(a, b):
        assert abs(x - y) <= 1e-9 * max(abs(x), 1.0)
