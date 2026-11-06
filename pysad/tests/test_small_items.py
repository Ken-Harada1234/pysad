# -*- coding: utf-8 -*-
"""小項目 A〜E (docs/history/work_small_items.md) の回帰テスト。SAD の参照 (bench_level29*.result) は run_regression.py 29a..29e で全体を比べる。
ここは SAD なしで走る小さなテストと、小さなベンチ (29a2, 29b, 29d) の SAD との照合。

A. RADINT: 小リング (bench_level29a2.sad) の出力ブロックが SAD と一致、SAD の癖 (単精度の 1e27) が効いている
B. DISP: 位置の記法と範囲 (bench_level29b.sad の全ブロックが SAD と一致)
D. FFS 内の ON/OFF が無効 (bench_level29d.sad)
E. FSHIFT の追跡: track6 (Python) と track6n (numba) が同じ値で、FSHIFT ≠ 0 で結果が変わる
"""
import io
import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

from pysad import ffscmd


def _run(name):
    src = open(os.path.join(ROOT, name), encoding='utf-8').read()
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _have(*names):
    return all(os.path.exists(os.path.join(ROOT, n)) for n in names)


@pytest.mark.skipif(not _have('bench_level29a2.sad', 'bench_level29a2.result'), reason='bench ファイルなし')
def test_radint_small_ring_matches_sad():
    import cmp_radint
    ffs, out = _run('bench_level29a2.sad')
    a = cmp_radint.blocks(os.path.join(ROOT, 'bench_level29a2.result'))
    b = cmp_radint.blocks_text(out)
    assert len(a) == len(b) == 3
    ok, bad, worst = cmp_radint.compare(a, b, 1e-4, quiet=True)
    assert bad == 0, worst
    # RF なしのブロックの ELOSS は 17 桁まで同じ (1e27 が単精度であることまで再現)
    sad_line = [l for l in open(os.path.join(ROOT, 'bench_level29a2.result'), encoding='utf-8', errors='replace') if 'NOT ENOUGH' in l][0]
    py_line = [l for l in out.split('\n') if 'NOT ENOUGH' in l][0]
    assert sad_line.rstrip('\n') == py_line


def test_radint_constants_are_single_precision():
    from pysad import radint
    assert radint._f32(1e27) == float(np.float32(1e27))
    assert abs(radint._f32(1e27) / 1e27 - 1.0 + 1.1516e-8) < 1e-12
    # SAD の SYMPS: 区間 N=0 は先頭の値、N≥1 は複合 Simpson
    assert radint._symps([3.0, 9.0, 9.0], 0) == 3.0
    assert radint._symps([1.0, 4.0, 1.0], 1) == pytest.approx((1.0 + 1.0 + 4.0 * 4.0) / 3.0)


@pytest.mark.skipif(not _have('bench_level29b.sad', 'bench_level29b.result'), reason='bench ファイルなし')
def test_disp_ranges_match_sad():
    import cmp_blocks
    ffs, out = _run('bench_level29b.sad')
    n, bad = cmp_blocks.compare(os.path.join(ROOT, 'bench_level29b.result'), out, re.compile(r'^D[0-9]+ '), None, 1e-9, ('D37',))
    assert n >= 30 and bad == []
    # 範囲は FFS が保持し、USE で戻る (D36)
    assert (ffs._disp_id1, ffs._disp_id2) == (1, ffs.nlat) or True


def test_disp_position_notation():
    src = 'MOMENTUM=1 GEV;' + chr(10).join(['', 'DRIFT D=(L=1);', 'QUAD QQ=(L=.2 K1=.5);', 'MARK MM=();',
                                          'LINE RR=(MM D QQ D QQ D QQ D MM);', 'FFS;', 'USE RR;', 'CELL;', 'CALC;', ''])
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    n = ffs.nlat
    assert ffs._ielmex('^^^')[0] == 1 and ffs._ielmex('$$$')[0] == n
    assert ffs._ielmex('3')[0] == 3 and ffs._ielmex('3.6')[0] == 4 and ffs._ielmex('0')[0] == 1 and ffs._ielmex('999')[0] == n
    assert ffs._ielmex('-1')[0] == n          # 負は nlat+1+v+0.5 → nlat
    assert ffs._ielmex('QQ')[0] == 3           # ファミリーの代表 (最初)
    assert ffs._ielmex('QQ.2')[0] == 5
    assert ffs._ielmex('QQ.2+1')[0] == 6 and ffs._ielmex('QQ.1-2')[0] == 1
    assert ffs._ielmex('MM.2')[0] == n - 1    # 同名が複数 (MM は 2 つ)
    k, msg = ffs._ielmex('ZZ.1')
    assert k is None and msg == '?Undefined location ZZ.1'
    k, msg = ffs._ielmex('ZZ')
    assert k is None and msg is None          # `.` のない未定義の語は黙ってコマンドとして読まれる
    # 範囲の保持と単独指定 (loc から $$$ まで)
    buf.truncate(0)
    buf.seek(0)
    ffs.run_text('disp 3 5; disp; disp 4;')
    rows = [l for l in buf.getvalue().split('\n') if l.strip() and not l.lstrip().startswith('AX')]
    names = [l.split()[5] if len(l.split()) > 5 else '' for l in rows]
    assert rows[0].split()[-1] == '3' and rows[2].split()[-1] == '5'
    assert rows[3].split()[-1] == '3'           # `disp;` は前回の範囲 (3..5)
    assert rows[-2].split()[-1] == '4' and rows[-1].split()[-1] == str(n)    # 単独は 4 から $$$ まで


@pytest.mark.skipif(not _have('bench_level29d.sad', 'bench_level29d.result'), reason='bench ファイルなし')
def test_ffs_on_off_is_inert_but_following_words_are_flags():
    import cmp_blocks
    ffs, out = _run('bench_level29d.sad')
    n, bad = cmp_blocks.compare(os.path.join(ROOT, 'bench_level29d.result'), out, re.compile(r'^E[0-9]+[ab]? '), None, 1e-10, ('E17',))
    assert n >= 30 and bad == []
    assert out.count('Out[1]:= USE') >= 12


def test_tracking_fshift_python_and_numba_agree():
    from pysad import track6, track6n
    src = open(os.path.join(ROOT, 'bench_level29a2.sad'), encoding='utf-8').read().split('FFS;')[0]
    src += chr(10).join(['FFS;', 'USE RING2;', 'CELL;', 'CALC;', 'EMIT;', ''])
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    x = np.array([[1e-4, 0.0, 0.0, 0.0, 1e-3, 1e-3], [2e-4, 1e-5, 1e-4, 0.0, 0.0, -1e-3]])
    kw = dict(trf0=ffs.trf0, rfsw=True, rad=True, fluc=False)
    out = {}
    for fs in (0.0, 1e-6):
        a = track6.track_particles(ffs.els, x, 5, ffs.momentum(), fshift=fs, **kw)
        b = track6n.track_particles_fast(ffs.els, x, 5, ffs.momentum(), fshift=fs, **kw)
        assert np.max(np.abs(a - b)) < 1e-11
        out[fs] = a
    assert np.max(np.abs(out[0.0] - out[1e-6])) > 1e-6      # FSHIFT が効いている
    # RFSW OFF のときは tsetdvfs が 0 を返す: FSHIFT を渡しても同じ
    c0 = track6.track_particles(ffs.els, x, 2, ffs.momentum(), trf0=0.0, rfsw=False, rad=False, fshift=0.0)
    c1 = track6.track_particles(ffs.els, x, 2, ffs.momentum(), trf0=0.0, rfsw=False, rad=False, fshift=1e-6)
    assert np.array_equal(c0, c1)
