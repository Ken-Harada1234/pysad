"""v2 の FFS まわり (担当 B): B22 MON/MCOD/TCOD、B24 LINE["LENGTH",loc]、B26 ORG の負数、B27 SeedRandom、B30 MAIN の OFF RFSW、B50 FileByteCount、
B67 MOMENTUM 未設定の警告、NOKEEPEXP、RENUMBER と INDEX 属性。

根拠: NOKEEPEXP と INDEX は gs (K64) の実測 (2026-11-17)。B30 は「FFS の NORFSW と同じ追跡結果になる」、B26 は「K64 で曖昧さのない書き方 (0-2.5) と同じ結果」。
使った条件: 小さな試験用リング (下の _RING: 1 GeV、周長 6 m、または _FODO: 3 GeV の FODO 24 セル + 空洞)。"""
import io
import os
import re
import tempfile

import pytest

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

_RING = '''DRIFT D1=(L=1);
QUAD QF=(L=.2 K1=.1 %s) QD=(L=.2 K1=-.1) QS=(L=.1 K1=.01 INDEX=3);
BEND B3=(L=1 ANGLE=.1);
MARK M1=(BX=10 BY=10);
MONI PM1=();
LINE RING=(M1 QF D1 QD D1 QF D1 QD D1 QS D1 QF D1);
LINE RINGB=(M1 QF D1 B3 D1 QD D1 PM1 QF D1);
'''

_FODO = '''DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 24*CELL);
'''


def _run(main, body, mode):
    from pysad import compat
    from pysad.ffscmd import run_sad_string
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        ffs, I = run_sad_string(main + 'FFS;\nPageWidth=1999;\n' + body + '\nExit;\n', out=buf)
    finally:
        compat.set_mode(old)
    return [l.rstrip() for l in buf.getvalue().splitlines()], ffs


def _nums(s):
    return [float(t) for t in re.findall(_NUM, s)]


def _line(lines, tag):
    for l in lines:
        if l.startswith(tag):
            return l[len(tag):]
    raise AssertionError('no line %r in %r' % (tag, lines[-20:]))


# ---------------------------------------------------------------- B22 / B24 / B26
def test_b22_mon_mcod_tcod():
    main = 'MOMENTUM = 1 GEV;\n' + _RING % ''
    body = 'USE RINGB;\nINS; CALC;\nMON 0 0 0 0 PM*;\nPrint["after MON"];\nMCOD 1 2;\nTCOD;\nPrint["end"];'
    r, _ = _run(main, body, 'v2')
    assert any('MON is not a command of SAD K64' in l for l in r)
    assert any('MCOD is not a command' in l for l in r) and any('TCOD is not a command' in l for l in r)
    assert not any('Missing value' in l for l in r)                              # 別のエラーで止まらない
    assert 'after MON' in r and 'end' in r
    r, _ = _run(main, body, 'k64')
    assert any('Missing value for PM*' in l for l in r)                          # gs と同じ
    assert not any('is not a command' in l for l in r)


def test_b24_line_length_with_location():
    main = 'MOMENTUM = 1 GEV;\n' + _RING % ''
    body = 'USE RINGB;\nINS; CALC;\nPrint["L ",LINE["LENGTH"]," ",LINE["POSITION","$$$"]," ",LINE["S","$$$"]," ",LINE["LENG","$$$"]];\nPrint[LINE["LENGTH","$$$"]];'
    for mode in ('v2', 'k64'):
        r, _ = _run(main, body, mode)
        v = _nums(_line(r, 'L '))
        assert v[0] == 11.0 and v[1] == 11.0 and abs(v[2] - 5.6) < 1e-12 and abs(v[3] - 5.6) < 1e-12
        err = [l for l in r if 'General::narg' in l]
        assert err and 'LINE["LENGTH","$$$"]' in err[0]                          # K64 のソース (tftwiss.f) も引数 1 個だけ
        assert ('LINE["POSITION",loc]' in err[0]) == (mode == 'v2')              # v2 は位置の取り方を知らせる


def test_b26_org_negative_numbers():
    main = 'MOMENTUM = 1 GEV;\n' + _RING % ''
    pr = 'Print["G ",LINE["GX","B3"]," ",LINE["GY","B3"]," ",LINE["GZ","B3"]," ",LINE["GX","$$$"]," ",LINE["GY","$$$"]," ",LINE["GZ","$$$"]];'
    amb = 'USE RINGB;\nINS; CALC;\nORG B3 1.5 -2.5 0.5 10 20 30;\n' + pr
    una = 'USE RINGB;\nINS; CALC;\nORG B3 1.5 0-2.5 0.5 10 20 30;\n' + pr          # K64 でも 6 個の値として読まれる書き方
    r2, f2 = _run(main, amb, 'v2')
    rk, fk = _run(main, una, 'k64')
    assert not any('Missing parameters' in l for l in r2)
    assert f2.org_g == [1.5, -2.5, 0.5] and f2.org_g == fk.org_g and list(f2.org_chi) == list(fk.org_chi)
    assert _nums(_line(r2, 'G ')) == _nums(_line(rk, 'G '))
    # v2 でも式は書ける (演算子で終わる語・* / ^ で始まる語は続き)
    r3, f3 = _run(main, 'USE RINGB;\nINS; CALC;\nORG B3 3/2 -5 / 2 0.5 2 *5 20 30;\n' + pr, 'v2')
    assert f3.org_g == [1.5, -2.5, 0.5] and list(f3.org_chi) == list(f2.org_chi)
    # k64: 1.5 -2.5 は 1.5-2.5 = -1.0 と読まれ、値が 1 つ足りなくなる (K64 の getvad)
    r4, f4 = _run(main, amb, 'k64')
    assert f4.org_g == [-1.0, 0.5, 10.0] and any('Missing parameters' in l for l in r4)


# ---------------------------------------------------------------- B30 MAIN の OFF RFSW
_TRK = ('x=TrackParticles[{1,{{1e-4},{0},{0},{0},{1e-3},{1e-4},{1}}},LINE["LENGTH"]];\n'
        'Print["T ",x[[2,1,1]]," ",x[[2,5,1]]," ",x[[2,6,1]]];\nPrint["F ",FLAG[{"RFSW","RAD","FLUC"}]];')


def test_b30_main_off_rfsw_reaches_ffs():
    lat = 'MOMENTUM = 3 GEV;\n' + _FODO
    r_v2_off, _ = _run(lat + 'OFF RFSW RAD FLUC;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'v2')
    r_k_off, _ = _run(lat + 'OFF RFSW RAD FLUC;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'k64')
    r_k_no, _ = _run(lat + 'OFF RAD FLUC;\n', 'USE RING;\nNORFSW;\nCELL; CALC;\n' + _TRK, 'k64')
    r_v2_on, _ = _run(lat + 'OFF RAD FLUC;\nON RFSW;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'v2')
    assert _nums(_line(r_v2_off, 'F ')) == [0.0, 0.0, 0.0]                       # v2: 届く
    assert _nums(_line(r_k_off, 'F ')) == [1.0, 0.0, 0.0]                        # K64: 届かない (track.f:149)
    t_off = _nums(_line(r_v2_off, 'T '))
    r_v2_no, _ = _run(lat + 'OFF RAD FLUC;\n', 'USE RING;\nNORFSW;\nCELL; CALC;\n' + _TRK, 'v2')
    r_v2_rf, _ = _run(lat + 'OFF RFSW RAD FLUC;\n', 'USE RING;\nRFSW;\nCELL; CALC;\n' + _TRK, 'v2')
    r_k_on, _ = _run(lat + 'OFF RAD FLUC;\nON RFSW;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'k64')
    # 比べる相手は同じ互換モードの結果。k64 との比較は v2 の四極の 4 次の合成 (B226) の分 (1e-9 台) だけずれるので、
    # 比較は「同じ互換モード同士で全桁同じ」と「k64 との差が積分器の差の範囲」に分ける
    assert t_off == _nums(_line(r_v2_no, 'T '))                                  # FFS の NORFSW と全桁同じ追跡
    assert t_off[2] == 1e-4                                                      # 空洞なし・放射なし: δ は保存
    t_on = _nums(_line(r_v2_on, 'T '))
    assert t_on == _nums(_line(r_v2_rf, 'T ')) and abs(t_on[2] - 1e-4) > 1e-6    # 空洞ありは δ が変わる。FFS の RFSW (v2) も同じ
    # k64: MAIN の OFF RFSW は届かず (空洞あり)、ON RFSW と同じ追跡。v2 との差は四極の積分器だけ (1e-7 以内)
    assert _nums(_line(r_k_off, 'T ')) == _nums(_line(r_k_on, 'T '))
    for k, v in ((_nums(_line(r_k_no, 'T ')), t_off), (_nums(_line(r_k_on, 'T ')), t_on)):
        assert max(abs(a - b) for a, b in zip(k, v)) < 1e-7
    # 四極を K64 の方式 (方式 0) にすれば、v2 の追跡は k64 と全桁同じ
    from pysad import quadbody
    old = quadbody.set_scheme(0)
    try:
        r_v2_s0, _ = _run(lat + 'OFF RFSW RAD FLUC;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'v2')
        r_v2_s0on, _ = _run(lat + 'OFF RAD FLUC;\nON RFSW;\n', 'USE RING;\nCELL; CALC;\n' + _TRK, 'v2')
    finally:
        quadbody.set_scheme(old)
    assert _nums(_line(r_v2_s0, 'T ')) == _nums(_line(r_k_no, 'T '))
    assert _nums(_line(r_v2_s0on, 'T ')) == _nums(_line(r_k_on, 'T '))


# ---------------------------------------------------------------- B67 MOMENTUM
def test_b67_momentum_warning():
    w = 'MOMENTUM is not set'
    body = 'USE RING;\nNORAD;NOFLUC;\nCELL; CALC;\nPrint["before"];\n' + _TRK + '\n' + _TRK + '\nPrint["M ",MOMENTUM];'
    r, ffs = _run(_FODO, body, 'v2')
    assert sum(w in l for l in r) == 1                                           # 1 度だけ
    assert r.index('before') < [i for i, l in enumerate(r) if w in l][0]         # 光学 (CALC) では出さず、運動量を使う計算 (追跡) で出す
    assert _nums(_line(r, 'M ')) == [3e10] and ffs.momentum() == 30e9           # 値は K64 の既定のまま
    r, _ = _run(_FODO, body, 'k64')
    assert not any(w in l for l in r)                                            # K64: 黙って 30 GeV
    r, _ = _run('MOMENTUM = 3 GEV;\n' + _FODO, body, 'v2')
    assert not any(w in l for l in r)
    r, _ = _run(_FODO, 'USE RING;\nMOMENTUM=3e9;\nNORAD;NOFLUC;\nCELL; CALC;\n' + _TRK, 'v2')
    assert not any(w in l for l in r)                                            # FFS の中で設定しても出ない
    r, _ = _run(_FODO, 'USE RING;\nMOMENTUM=30e9;\nNORAD;NOFLUC;\nCELL; CALC;\n' + _TRK, 'v2')
    assert not any(w in l for l in r)                                            # 明示の 30 GeV は警告しない
    r, _ = _run(_FODO, 'USE RING;\nCELL; CALC;\nEMIT;', 'v2')
    assert sum(w in l for l in r) == 1                                           # EMIT でも出る
    r, _ = _run('', 'Print["x ",1+1];', 'v2')
    assert not any(w in l for l in r)                                            # ラティスを使わないスクリプトでは出ない


# ---------------------------------------------------------------- NOKEEPEXP / RENUMBER と INDEX
_KEXP = '''USE RING;
INS; CALC;
pr[t_]:=Print[t," ",LINE["K1","QF.*"]," BX=",Twiss["BX","$$$"]," ",FLAG[{"KEEPEXP","CALEXP"}]];
pr["a"];
LINE["K1","QF.2"]=0.3;
CALC NOEXPAND; pr["b"];
CALC; pr["c"];
NOKEEPEXP;
LINE["K1","QF.2"]=0.3;
CALC NOEXPAND; pr["d"];
CALC; pr["e"];
LINE["K1","QF.2"]=0.3;
GO NOEXPAND; pr["f"];
GO; pr["g"];
KEEPEXP;
LINE["K1","QF.2"]=0.3;
CALC NOEXPAND; pr["h"];
GO; pr["i"];
CALC EXPAND; pr["j"];'''

# gs (SAD 1.4.4.2k64、2026-11-17): 各行は LINE["K1","QF.*"] の 2 番目、BX($$$)、KEEPEXP、CALEXP
_KEXP_GS = {'a': (.1, 7.322092408212155, 1, 1), 'b': (.3, 1.682490579872757, 1, 0), 'c': (.3, 1.682490579872757, 1, 0),
            'd': (.3, 1.682490579872757, 0, 0), 'e': (.1, 7.322092408212155, 0, 1), 'f': (.3, 1.682490579872757, 0, 0),
            'g': (.1, 7.322092408212155, 0, 1), 'h': (.3, 1.682490579872757, 1, 0), 'i': (.3, 1.682490579872757, 1, 0),
            'j': (.1, 7.322092408212155, 1, 1)}


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_nokeepexp_matches_gs(mode):
    main = ('MOMENTUM = 1 GEV;\nDRIFT D1=(L=1);\nQUAD QF=(L=.2 K1=.1) QD=(L=.2 K1=-.1);\nMARK M1=(BX=10 BY=10);\n'
            'LINE RING=(M1 QF D1 QD D1 QF D1 QD D1 QF D1);\n')
    r, _ = _run(main, _KEXP, mode)
    for t, (k2, bx, ke, ce) in _KEXP_GS.items():
        v = _nums(_line(r, t + ' '))
        assert v[0] == .1 and v[2] == .1 and v[1] == k2, (t, v)
        assert abs(v[3] - bx) < 1e-12 * bx and v[4:] == [float(ke), float(ce)], (t, v)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_renumber_with_index_attribute_matches_gs(mode):
    main = 'MOMENTUM = 1 GEV;\n' + _RING % 'INDEX=5'
    body = ('USE RING;\nINS; CALC;\nPrint["N ",LINE["NAME"]];\nPrint["M ",LINE["MULT"]];\n'
            'Print["P ",LINE["POSITION","QF.6"]," ",LINE["POSITION","QF.1"]," ",LINE["POSITION","QF.7"]," ",LINE["POSITION","QS"]," ",LINE["INDEX","QF"]];\n'
            'Print["Q ",LINE["NAME","QF.*"]];\nRENUMBER QD.2;\nPrint["N2 ",LINE["NAME"]];\nPrint["M2 ",LINE["MULT"]];\nPrint["K ",LINE["K1","QF.7"]," ",LINE["POSITION","QF.5"]];')
    r, _ = _run(main, body, mode)
    # gs の出力 (INDEX=5 の QF は QF.5 から。1 台だけの QS は INDEX があっても番号なし)
    assert _line(r, 'N ') == '{"M1","QF.5","D1.1","QD.1","D1.2","QF.6","D1.3","QD.2","D1.4","QS","D1.5","QF.7","D1.6","$$$"}'
    assert _line(r, 'M ') == '{0,5,1,1,2,6,3,2,4,0,5,7,6,0}'
    assert _line(r, 'P ') == '6 {} 12 10 5'
    assert _line(r, 'Q ') == '{"QF.5","QF.6","QF.7"}'
    assert _line(r, 'N2 ') == '{"M1","QF.6","D1.4","QD.2","D1.5","QF.7","D1.6","QD.1","D1.1","QS","D1.2","QF.5","D1.3","$$$"}'
    assert _line(r, 'M2 ') == '{0,6,4,2,5,7,6,1,1,0,2,5,3,0}'
    assert _line(r, 'K ') == '.1 12'                                             # 名前と位置は LINE["NAME"] の表示と一致する
    # INDEX が無ければ従来どおり 1 から
    r, _ = _run('MOMENTUM = 1 GEV;\n' + _RING % '', 'USE RING;\nINS; CALC;\nPrint["M ",LINE["MULT"]];', mode)
    assert _line(r, 'M ') == '{0,1,1,1,2,2,3,2,4,0,5,3,6,0}'


# ---------------------------------------------------------------- B27 / B50
@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_b27_seedrandom_display_unsigned(mode):
    r, _ = _run('', 'SeedRandom[-5];Print["S ",SeedRandom[]];SeedRandom[17];Print["T ",SeedRandom[]];', mode)
    assert _line(r, 'S ').startswith('{"SAD",4294967291,')                       # 符号なし 32 bit (seed mod 2^32)
    assert _line(r, 'T ').startswith('{"SAD",17,')


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_b50_filebytecount_is_bytes(mode):
    d = tempfile.mkdtemp()
    p = os.path.join(d, 'two_lines.txt')
    with open(p, 'wb') as f:
        f.write(b'abcdef\n123456\n')                                             # 14 バイト、2 行 (K64 は wc -l で 2 を返す)
    try:
        r, _ = _run('', 'Print["B ",FileByteCount["%s"]," ",FileByteCount["%s"]];' % (p.replace(os.sep, '/'), (p + '.none').replace(os.sep, '/')), mode)
    finally:
        os.remove(p)
        os.rmdir(d)
    assert _line(r, 'B ') == '14 $Failed'


# ---------------------------------------------------------------- 実際の使い方: 文脈つきのパッケージ + DP0 の対称スキャン (B18)
_CHRO = '''USE RING;
NORAD; NOFLUC;
CELL; CALC;
BeginPackage[Chro`];
Begin[Chro`];
tune[dp_] := Module[{nx}, DP0 = dp; FFS["CALC"]; nx = Twiss["NX", "$$$"]/(2*Pi); nx];
`ChromScan[dps_List, n_] := Module[{d = Map[{#, tune[#]}&, dps], c}, c = PolynomialFit[d, n][[1]]; DP0 = 0; FFS["CALC"]; c];
End[];
EndPackage[];
c = ChromScan[{-0.002, -0.001, 0, 0.001, 0.002}, 2];
Print["C ", c];
DP0 = 0.001; CALC; np = Twiss["NX", "$$$"]/(2*Pi);
DP0 = -0.001; CALC; nm = Twiss["NX", "$$$"]/(2*Pi);
DP0 = 0; CALC; n0 = Twiss["NX", "$$$"]/(2*Pi);
Print["D ", n0, " ", (np - nm)/0.002, " ", (np - 2*n0 + nm)/0.001^2/2];
Print["E ", tune[0], " ", Chro`tune[0] - n0];'''


def test_chromaticity_scan_in_user_package():
    """利用者のパッケージ (文脈 Chro`) から FFS の CALC・Twiss と PolynomialFit を使い、DP0 を対称に振ってクロマティシティを求める。
    v2 は差分で求めた値と一致する。K64 (k64 モード) は x の平均が 0 なので係数が壊れる (B18)。条件: 3 GeV の FODO 24 セル (6 極 0)"""
    main = 'MOMENTUM = 3 GEV;\n' + _FODO
    r, _ = _run(main, _CHRO, 'v2')
    c = _nums(_line(r, 'C '))
    d = _nums(_line(r, 'D '))
    assert abs(c[0] - d[0]) < 1e-9                                               # ν(0)
    assert abs(c[1] - d[1]) < 2e-3 * abs(d[1]) and d[1] < -1.0                   # 自然クロマティシティ dν/dδ (負)。2 次の当てはめと中心差分
    assert abs(c[2] - d[2]) < 0.05 * abs(d[2])                                   # 2 次の係数
    e = _line(r, 'E ')
    assert e.startswith('tune[0] ') and abs(_nums(e)[-1]) < 1e-12                # 文脈の中の tune は外から Chro`tune でだけ呼べる
    rk, _ = _run(main, _CHRO, 'k64')
    ck = _nums(_line(rk, 'C '))
    assert abs(ck[1] - d[1]) > 0.5 * abs(d[1])                                   # K64: 1 次の係数が壊れる (c1 + 2 c2 になる)
