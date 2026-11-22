# -*- coding: utf-8 -*-
"""F3: PUTMEASURE / GETMEASURE を Unix の FIFO (名前つきパイプ) で受け渡す (K64 の nfgetm.c は SADPOUT・SADPIN に FIFO を想定)。
書き手 (最初の PUTMEASURE の fopen) と読み手 (最初の GETMEASURE の fopen) は、相手が開くまで待つ (open が止まる)。
POSIX の mkfifo が無い環境 (Windows) では skip。gs の PUTMEASURE・GETMEASURE は FFS から呼ぶ MEASURE が先頭の MARK 以外では動かない (B151) ので
gs との比較は無い。基準は nfgetm.c の書式 (`%15.7lg` ×5 をコンマ区切り) とソースの読み。"""
import os
import threading

import pytest

from pysad import measure

pytestmark = pytest.mark.skipif(not hasattr(os, 'mkfifo'), reason='mkfifo (POSIX) が無い')


class _F:
    def __init__(self):
        self.out = []

    def write(self, s):
        self.out.append(s)


def test_put_get_through_fifo(tmp_path):
    path = str(tmp_path / 'm.fifo')
    os.mkfifo(path)
    old = (os.environ.get('SADPOUT'), os.environ.get('SADPIN'))
    os.environ['SADPOUT'] = path
    os.environ['SADPIN'] = path
    got = []
    rd = _F()

    def reader():                              # 読み手 (GETMEASURE 側): open が書き手を待つ
        for _ in range(3):
            got.append(measure.read_measure_record(rd))
        measure.close_measure_files(rd)

    th = threading.Thread(target=reader, daemon=True)
    th.start()
    wr = _F()
    vals = [[1.5e-4, -2.5e-5, 3.25e-9, -1.0e-10, 4.0e-9], [0.0, 1.0, 2.0, 3.0, 4.0], [-1e-3, 2e-3, 3e-6, 4e-7, 5e-8]]
    try:
        for v in vals:
            wr._measure_sv = [v[0], v[1], v[2], v[4], v[3]]       # cmd_putmeasure は (xa, ya, xxa, yya, xya) の順で持つ
            assert measure.cmd_putmeasure(wr, [], 0) == 1
        th.join(timeout=20)
        assert not th.is_alive()
        measure.close_measure_files(wr)
    finally:
        for k, v in zip(('SADPOUT', 'SADPIN'), old):
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    assert not wr.out and not rd.out
    for g, v in zip(got, vals):
        assert g is not None
        assert g == pytest.approx(v, rel=1e-6, abs=0.0)     # %15.7g の桁
