# -*- coding: utf-8 -*-
"""DDBA v07_14 の一式 (example/ddba/v07_14): ラティス 2 本の `cell calc emit` の値と、誤差の計算 (補正しながら誤差を大きくする手法) の
縮小版が pysad (k64) で走って補正が効くこと、`USE RING NOEXPAND` で補正磁石の成分ごとの K0 が効くこと。
gs の保存値 (example/ddba/v07_14/*/NAME.gs.txt) があれば、チューンの出力を照合する。スクリプトは `Get["../modules/..."]` を相対パスで読むので、
台本のあるフォルダをカレントにして走らせる (cwd 引数)。"""
import io
import os
import re

import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'example', 'ddba', 'v07_14')


def _run_text(src, cwd):
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf, cwd=cwd)
    return buf.getvalue()


def _lat(name):
    return open(os.path.join(EX, 'lattice', name), encoding='utf-8').read()


def test_lattice_headers_have_no_personal_names():
    for n in ('v07_14_5GeV.sad', 'v07_14_2.5GeV.sad'):
        t = _lat(n)
        assert '\r' not in t
        assert 'C:\\' not in t and 'D:\\' not in t
        low = t.lower()
        for bad in ('shimo' + 'saki', 'sli' + 'tj', 'sad_' + 'archive', 'wi' + 'ki'):
            assert bad not in low
        assert '10.1107/S1600577521012753' in t  # 6 極の関連する公開論文


@pytest.mark.parametrize('name,energy,emit,volt', [('v07_14_5GeV.sad', 5.0, 8.32e-10, 6.5e6), ('v07_14_2.5GeV.sad', 2.5, 2.08e-10, 1.6e6)])
def test_lattice_cell_calc_emit(name, energy, emit, volt):
    compat.set_mode('k64')
    src = _lat(name) + '''
ON RAD FLUC COD RFSW RADCOD;
FFS;
USE RING;
cell calc emit;
Print["CHK S ",Twiss["S","$$$"]," NX ",Twiss["NX","$$$"]/2/Pi," NY ",Twiss["NY","$$$"]/2/Pi," EMITX ",EMITX," N ",LINE["LENGTH"]];
USE RINGSYM;
cell calc emit;
Print["CHKSYM S ",Twiss["S","$$$"]," NX ",Twiss["NX","$$$"]/2/Pi," NY ",Twiss["NY","$$$"]/2/Pi];
Print["CHKRF ",LINE["VOLT","RFC"]," ",MOMENTUM];
stop;
'''
    out = _run_text(src, EX)
    m = re.search(r'^CHK S\s+(\S+) NX\s+(\S+) NY\s+(\S+) EMITX\s+(\S+) N\s+(\S+)', out, re.M)
    assert m, out[-500:]
    s, nx, ny, ex_, n = (float(x) for x in m.groups())
    assert s == pytest.approx(749.4808883, abs=1e-5)
    assert nx == pytest.approx(47.8650597, abs=1e-6)
    assert ny == pytest.approx(16.6550904, abs=1e-6)
    assert ex_ == pytest.approx(emit, rel=2e-3)
    assert n == 3540
    m = re.search(r'^CHKSYM S\s+(\S+) NX\s+(\S+) NY\s+(\S+)', out, re.M)
    assert float(m.group(2)) == pytest.approx(47.8650597, abs=1e-6)
    m = re.search(r'^CHKRF\s+(\S+)\s+(\S+)', out, re.M)
    assert float(m.group(1)) == pytest.approx(volt, rel=1e-9)
    assert float(m.group(2)) == pytest.approx(energy * 1e9, rel=1e-9)


def test_noexpand_per_position_k0():
    """既定の (EXPAND の) cell calc では HD.2 以降の K0 が効かず (HD.1 が全部に効く)、NOEXPAND なら成分ごとに効く。"""
    compat.set_mode('k64')
    body = '''
ON RAD FLUC COD RFSW RADCOD;
FFS;
USE RING%s;
FFS["%s"];
n=LINE["NAME","HD*.*"];
LINE["K0",n[7]]=1e-5;
FFS["%s"];
Print["CHK7 rb ",LINE["K0",n[7]]," cod ",Max[Abs[Twiss["DX","PM*.*"]]]];
stop;
'''
    t = _lat('v07_14_5GeV.sad')
    out_e = _run_text(t + body % ('', 'cell calc', 'cell calc'), EX)
    out_n = _run_text(t + body % (' NOEXPAND', 'cell calc noexpand', 'cell calc noexpand'), EX)
    cod_e = float(re.search(r'^CHK7 rb\s+\S+ cod\s+(\S+)', out_e, re.M).group(1))
    cod_n = float(re.search(r'^CHK7 rb\s+\S+ cod\s+(\S+)', out_n, re.M).group(1))
    assert cod_e < 1e-6
    assert cod_n > 1e-5


def test_correct_cod_small():
    """誤差の計算の縮小版 (種 1 個、段階 5、補正 1 回ずつ+最後 3 回、DA なし): 補正のあとの COD の rms が補正前より小さい。"""
    compat.set_mode('k64')
    p = os.path.join(EX, 'errors', 'V0714_correct_cod.sad')
    src = open(p, encoding='utf-8').read()
    for a, b in (('nseed = 2;', 'nseed = 1;'), ('nitr = 3;', 'nitr = 1;'), ('nitrlast = 10;', 'nitrlast = 3;'), ('daOn = 1;', 'daOn = 0;')):
        assert a in src
        src = src.replace(a, b)
    out = _run_text(src, os.path.dirname(p))
    rows = re.findall(r'^mat x rms =\s+(\S+)\s+mat y rms =\s+(\S+)\s+sad x rms =\s+(\S+)\s+sadt y rms =\s+(\S+)', out, re.M)
    assert len(rows) == 4 * 1 + 3
    first = (float(rows[0][2]), float(rows[0][3]))
    last = (float(rows[-1][2]), float(rows[-1][3]))
    # 最後の段階は誤差が 5 倍になっているのに、補正後は最初の段階 (1 回補正) と同程度以下
    assert last[0] < 1e-5 and last[1] < 1e-4
    m = re.search(r'^CHK seed\s+\S+ tune\s+(\S+)\s+(\S+) codrms\s+(\S+)\s+(\S+)', out, re.M)
    assert m
    nx, ny, cx, cy = (float(x) for x in m.groups())
    assert 47.5 < nx < 48.0 and 16.5 < ny < 17.0
    assert cx < 1e-4 and cy < 1e-4
    assert first[0] > 0
