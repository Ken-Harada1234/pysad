# -*- coding: utf-8 -*-
"""example/ddba/lattice/DDBA_NN_*.sad (設計メモの版 DDBA_01〜28 のラティス): `cell calc emit` が走り、周長・チューン・エミッタンスが
記録した値 (pysad k64 = gs と同じ式。example/ddba/results_k64.json) と一致する。内容が同じ台本 (DDBA_04・05 = 03、08・09 = 07、
26 = 25) は 1 本だけ走らせる (DDBA_06 の CAVI の `VOLT =1.7MV` も読める)。
追加したファイルに氏名・元のフォルダ名・内部のサイトの名前が入っていないことも調べる。"""
import io
import os
import re

import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
DD = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'example', 'ddba')
LAT = os.path.join(DD, 'lattice')
# (DDBA_NN, 周長 m, νx, νy, εx k64 [m rad])
CASES = [
    (1, 750.6000000, 52.8739311, 17.5347853, 5.20960e-10),
    (2, 731.4000000, 48.6000658, 12.7001616, 3.18170e-10),
    (3, 731.4000000, 48.6000682, 12.7001671, 1.27485e-09),
    (6, 731.4000000, 48.6000691, 12.7001725, 1.89100e-10),
    (7, 731.4000000, 48.6000691, 12.7001725, 7.56470e-10),
    (11, 731.4000000, 48.6000691, 12.7001725, 7.56770e-10),
    (12, 731.4000000, 47.1000664, 13.2001615, 8.75870e-10),
    (13, 731.4000000, 48.6000729, 12.7001569, 6.98430e-10),
    (14, 731.4000000, 48.6000691, 12.7001725, 7.56500e-10),
    (15, 731.4000000, 48.6000344, 12.7001425, 9.63190e-10),
    (16, 738.6898537, 48.6345358, 12.7446376, 9.21860e-10),
    (17, 738.6898537, 48.8650679, 12.6551782, 8.85630e-10),
    (18, 738.6898537, 48.6100651, 12.7201665, 8.62790e-10),
    (19, 744.8474014, 45.1301022, 12.2102246, 9.09150e-10),
    (20, 738.6898537, 48.6100651, 12.7201665, 8.62790e-10),
    (21, 738.6898537, 48.6100812, 12.7201742, 7.16960e-10),
    (22, 749.4811450, 47.8650520, 16.6650968, 8.30320e-10),
    (23, 757.0443181, 50.1350612, 16.3450645, 9.17670e-10),
    (24, 749.4808883, 47.8650597, 16.6550904, 8.31780e-10),
    (25, 749.4808883, 47.8650597, 16.6550904, 2.07940e-10),
    (27, 749.4808883, 47.8650597, 16.6550904, 8.32110e-10),
    (28, 749.4808883, 47.8650597, 16.6550904, 8.32190e-10),
]

CALC = '''
ON LOG;
ON RAD FLUC COD RFSW RADCOD;
FFS;
$FORM = "15.10"; PageWidth = 1999;
USE RING;
cell calc emit;
Print["circumference [m] = ", LINE["S","$$$"]];
Print["tune nx ny = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
Print["EMITX [m rad] = ", EMITX];
Exit;
'''


def _file(nn):
    names = [f for f in os.listdir(LAT) if f.startswith('DDBA_%02d_' % nn)]
    assert len(names) == 1
    return os.path.join(LAT, names[0])


def _run(nn):
    compat.set_mode('k64')
    src = open(_file(nn), encoding='utf-8').read()
    buf = io.StringIO()
    ffscmd.run_sad_string(src + CALC, out=buf)
    return buf.getvalue()


def _val(txt, pat):
    return [float(x) for x in re.search(pat, txt).groups()]


@pytest.mark.parametrize('nn,circ,nux,nuy,ex', CASES)
def test_cell_calc_emit_matches_recorded(nn, circ, nux, nuy, ex):
    txt = _run(nn)
    assert _val(txt, r'circumference \[m\] =\s*(\S+)')[0] == pytest.approx(circ, abs=1e-6)
    tx, ty = _val(txt, r'tune nx ny =\s*(\S+)\s+(\S+)')
    assert tx == pytest.approx(nux, abs=1e-6) and ty == pytest.approx(nuy, abs=1e-6)
    assert _val(txt, r'EMITX \[m rad\] =\s*(\S+)')[0] == pytest.approx(ex, rel=1e-4)


def test_all_28_lattice_files_exist_with_header_and_lf():
    for nn in range(1, 29):
        b = open(_file(nn), 'rb').read()
        assert b'\r' not in b
        txt = b.decode('utf-8')
        assert txt.startswith('! ====') and ('DDBA_%02d ' % nn) in txt.split('\n')[1]


def test_no_personal_names_or_site_names_in_ddba_files():
    parts = ['sli' + 'tj', 'sad_' + 'archive', 'wi' + 'ki', 'tr' + r'ac', 'ce' + 'rl', 'pfcon' + 'rg', 'v05_and_' + 'older', 'C:.Users', 'D:.work',
             'shimo' + 'saki', '下' + '崎']
    bad = re.compile('|'.join(parts), re.I)  # 名前は文字列を分けて書く (この試験ファイルに 1 つの語として載せない)
    for dp, dn, fn in os.walk(LAT):
        for f in fn:
            assert not bad.search(f), f
            assert not bad.search(open(os.path.join(dp, f), encoding='utf-8').read()), f
    for f in ('README.md', 'STAGES.md', 'ORIGINAL_MAP.md', 'FFS_cell_calc_emit_DDBA_28.md', 'FFS_cell_calc_emit_DDBA_28.sad'):
        assert not bad.search(open(os.path.join(DD, f), encoding='utf-8').read()), f
