"""MAIN レベルの TRACK (pysad/maintrack.py, parser.parse_main_commands) のテスト。
SAD との数値照合は bench_level27a_pysad.py (SADScript 版との比較)。ここは書式・規則・小さな端から端までの確認。"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, HERE)

from pysad.parser import parse_main_commands, parse_track_statement
from pysad import maintrack as mt


def test_parse_commands_order_and_filters():
    src = """
    ! comment
    ON RFSW RADCOD COD; OFF FLUC;
    EMIT USE=RING;
    PHSPLOTS=9;
    track[{nx_},n_]:=Module[{a}, a=1; a];     ! SADScript の関数定義は命令ではない
    emit = Emittance[];
    TRACK USE=RING NPART=3 TURNS=10 NX=(1 5) NY=(0,2) NZ=(1);
    STOP;
    TRACK USE=X NPART=1 TURNS=1;
    """
    c = parse_main_commands(src)
    assert [x['cmd'] for x in c] == ['ON', 'OFF', 'EMIT', 'SET', 'TRACK', 'STOP']
    assert c[0]['flags'] == ['RFSW', 'RADCOD', 'COD']
    assert c[2]['use'] == 'RING'
    assert c[3]['name'] == 'PHSPLOTS' and c[3]['value'] == 9.0
    t = c[4]
    assert (t['npart'], t['turns'], t['nx'], t['ny'], t['nz']) == (3, 10, [1.0, 5.0], [0.0, 2.0], [1.0])


def test_track_statement_defaults_and_override():
    t = parse_track_statement('TRACK USE=ring NPART=5 TURNS=1E3 NX=(1,2) NX=(7) PROTON')
    assert t['use'] == 'RING' and t['npart'] == 5 and t['turns'] == 1000
    assert t['nx'] == [7.0] and t['ny'] is None and t['ptype'] == 'PROTON'


def test_fortran_g_and_f_formats():
    # gfortran の 1P,G12.4 / G13.7 / F10.6 と照合した値
    assert mt.fmt_g(0.0, 12, 4) == '   0.000    '
    assert mt.fmt_g(1.0, 12, 4) == '   1.000    '
    assert mt.fmt_g(0.09999, 12, 4) == '  9.9990E-02'
    assert mt.fmt_g(1234.5678, 12, 4) == '   1235.    '
    assert mt.fmt_g(-0.5, 13, 7) == '-.5000000    '
    assert mt.fmt_g(3.0, 13, 7) == '3.000000    '.rjust(13)
    assert mt.fmt_f(-0.423015, 10, 6) == ' -0.423015'
    assert mt.ordinal(1) == 'st' and mt.ordinal(2) == 'nd' and mt.ordinal(11) == 'th' and mt.ordinal(22) == 'nd' and mt.ordinal(103) == 'rd'


def test_wtunem_pure_rotation():
    # z_n = x+i px = a exp(-2πi ν n) (正の tune で px = -x sinφ) → tcftr (exp(+i...)) のピークは f1 = +ν (tune 単位、[-.5,.5])。
    # 補間が効いて整数ビンの間の ν も 1e-3 以下で出る
    lsp = 256
    n = np.arange(lsp)
    nu = [0.2113, -0.3071, 0.0421]
    z = np.stack([np.exp(-2j * np.pi * nu[k] * n) * (k + 1) for k in range(3)], axis=1)
    f1, fa = mt.wtunem(z)
    for k in range(3):
        assert abs(f1[k, 0] - nu[k]) < 2e-3
        assert fa[k, 0] > fa[k, 1] * 1e3 or fa[k, 1] == 0.0


def test_echo_lists_wraps_at_70_columns():
    from pysad.sadscript import fmt_real
    ls = mt.echo_lists([float(i) for i in range(1, 12)], [0.0], None, fmt_real)
    assert ls[0] == ' NX=' + ''.join(str(i).ljust(7) for i in range(1, 11))     # 10 個で改行
    assert ls[1] == ' NX=' + '11'.ljust(7)
    assert ls[2] == ' NY=' + '0'.ljust(7)
    assert len(ls) == 3


def test_smear_report_format():
    sj = [np.array([2e-7, 0.0]), np.array([1e-9, 0.0]), np.array([3e-5, 0.0])]
    sjj = [sj[k] ** 2 * 1.01 for k in range(3)]
    lines, recs = mt.smear_report(100, 100, sj, sjj, np.array([0]), 2)
    assert lines[0] == ' Smear at    100th turn:'
    assert lines[1] == ' p. #   2Jx(m)      2Jy(m)      2Jz(m)      sigma2Jx  sigma2Jy  sigma2Jz'
    assert lines[2].startswith('     1  2.0000E-09  1.0000E-11  3.0000E-07')
    assert lines[3].startswith(' Avrg   2.0000E-09')
    assert recs[0][0] == 1


def test_initial_distribution_identity_map():
    c = mt.initial_distribution(4, [1.0, 2.0], None, [0.5], np.eye(6), np.zeros(6), 2.0, 3.0, tdx=0.1)
    assert np.allclose(c[:, 0], [2.0 * 1 + 0.1, 2.0 * 2 + 0.1, 2.0 * 1 + 0.1, 2.0 * 2 + 0.1])    # x = NX × emx (+TDXI)
    assert np.allclose(c[:, 2], 2.0)                                                          # NY 省略 → 1
    assert np.allclose(c[:, 5], 1.5) and np.allclose(c[:, 1], 0.0) and np.allclose(c[:, 4], 0.0)   # g = NZ × emz


def test_end_to_end_small_ring():
    from pysad.ffscmd import run_sad_string
    from pysad import track6 as t6
    lat = os.path.join(HERE, 'bench_level27a_main.sad')
    src = open(lat, encoding='utf-8').read().split('ON RFSW;')[0]
    src += """
ON RFSW;
OFF RAD;
EMIT USE=RING;
TRACK USE=RING NPART=3 TURNS=20 NX=(1,2) NY=(1) NZ=(0);
STOP;
"""
    t6.RNG.restore(17)
    buf = io.StringIO()
    ffs, I = run_sad_string(src, out=buf, cwd=HERE)
    out = buf.getvalue()
    assert ' RFSW RADCOD RAD  FLUC  INTRA POL   COD  DAPER EMIOU K64   FOURI SMEAR' in out
    assert ' Smear at     21th turn:' in out
    assert out.count(' nux       ampl        nuy       ampl        nuz       ampl') == 3
    assert ' Tracking             end:  CPU time =' in out
    assert ' NX=1      2      ' in out and ' NY=1      ' in out and ' NZ=0      ' in out
    res = ffs.track_results[0]
    assert res['alive'] == 3 and not res['losses']
    assert len(res['tunes']) == 3 and len(res['tunes'][0][1]) == mt.NFRE
