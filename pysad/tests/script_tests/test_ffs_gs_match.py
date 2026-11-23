"""FFS・光学・ラティス側の gs (K64) の実動作に合わせた箇所の試験。参照値は gs (2026-10-19 ビルド) で取った。
FFS 内の MOMENTUM、Twiss["LENGTH"/"GAMMA"]、USE/CALC の EXPAND・NOEXPAND、RENUMBER、LINE["MULT"]、TYPE の K5、
DISP の GA・Z・B・REFERENCE・DREFERENCE、REFERENCE と Twiss["R…"/"D…"]、DUMP (パターンは 20 文字まで)、DDK、OrbitGeo の DZ、`QF.3 K1 0.9;`。"""
import io
import re

import pytest

LAT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494);
MARK ORG=();
CAVI RF = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (QF D2 B D2 QD D1);
LINE RING = (ORG RF 24*CELL);
LINE RING8 = (ORG 8*CELL);
LINE RING4 = (ORG 4*CELL);
SEXT  SF = (L=0.2 K2=0);
LINE CELLS = (SF QF D2 B D2 QD D1);
LINE RINGS = (ORG RF 8*CELLS);
'''


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(LAT + 'FFS;\n$FORM="12.6";\n' + body + '\nExit;\n', out=buf)
    return buf.getvalue()


def _nums(s):
    return [float(x) for x in re.findall(r'-?\d+\.?\d*(?:[Ee][-+]?\d+)?|-?\.\d+(?:[Ee][-+]?\d+)?', s)]


def _lines(out, tag):
    return [l for l in out.splitlines() if l.startswith(tag + ' ')]


def test_momentum_set_in_ffs_survives_use():
    out = _run('''MOMENTUM = 1e9;
RAD; NOFLUC;
USE RING; CELL; CALC;
Print["M ", MOMENTUM];
r = TrackParticles[{1, {{3e-4},{0},{1e-4},{0},{0.01},{0},{1}}}, 1, 1, 50][[2]];
Print["T ", r[[1]], " ", r[[3]]];
EMIT;
Print["S ", SIGE];''')
    assert _nums(_lines(out, 'M')[0])[0] == 1e9
    t = _nums(_lines(out, 'T')[0])
    assert t[0] == pytest.approx(2.237707e-4, rel=1e-5)
    assert _nums(_lines(out, 'S')[0])[0] == pytest.approx(3.916263e-4, rel=1e-5)
    assert 'P0 = 1.0000000 GeV' in out


def test_twiss_length_gamma_s():
    out = _run('''USE RING4; CELL; CALC;
Print["L ", Twiss["LENGTH"]];
Print["G ", Twiss["GAMMA",1], " ", Twiss["GAMMABETA",2.5], " ", Twiss["S","B.1"]];''')
    assert _nums(_lines(out, 'L')[0])[0] == 26
    g = _nums(_lines(out, 'G')[0])
    assert g[0] == pytest.approx(5870.8536, rel=1e-7) and g[2] == pytest.approx(0.8)


def test_expand_noexpand():
    out = _run('''NORAD;
USE RING8; CELL; CALC;
t[s_]:=Module[{r}, r = TrackParticles[{1, {{3e-4},{0},{1e-4},{0},{0},{0},{1}}}, 1, 1, 20][[2]];
Print[s," ", LINE["K1","QF.2"], " ", r[[1]]]];
LINE["K1","QF.2"]=0.6; t["B"];
CALC; t["C"];
LINE["K1","QF.2"]=0.6; CALC NOEXPAND; t["D"];
CALC; t["E"];
LINE["K1","QF.2"]=0.5; CALC; t["F"];
LINE["K1","QF.2"]=0.6; CALC EXPAND; t["G"];
LINE["K1","QF.2"]=0.6; USE NOEXPAND RING8; t["H"];
USE RING8; t["I"];''')
    exp = {'B': (0.6, 8.4359e-4), 'C': (0.7, 3.3924e-4), 'D': (0.6, 8.4359e-4), 'E': (0.6, 8.4359e-4),
           'F': (0.5, -4.693e-4), 'G': (0.7, 3.3924e-4), 'H': (0.6, 8.4359e-4), 'I': (0.7, 3.3924e-4)}
    for k, (k1, x) in exp.items():
        v = _nums(_lines(out, k)[0])
        assert v[0] == pytest.approx(k1) and v[1] == pytest.approx(x, rel=2e-4), k
    assert out.count('Element values are not expanded') == 3
    assert 'NOCALEXP' in out


def test_family_dot_number_sets_nothing_and_save_spreads():
    out = _run('''USE RING4; CELL; CALC;
QF.3 K1 0.9;
Print["A ", LINE["K1",{"QF.1","QF.3"}]];
LINE["K1","QF.2"]=0.6; LINE["K1","QF.1"]=0.8; SAVE;
Print["B ", LINE["K1",{"QF.1","QF.2","QF.3"}]];''')
    assert _nums(_lines(out, 'A')[0]) == pytest.approx([0.7, 0.7])
    assert _nums(_lines(out, 'B')[0]) == pytest.approx([0.8, 0.8, 0.8])


def test_renumber_and_mult():
    out = _run('''USE RING4; CELL; CALC;
Print["m ", LINE["MULT",{"ORG","QF.1","QF.2","$$$","D2.2","D2.3","B.4"}]];
RENUMBER QF.3;
Print["d ", LINE["NAME",{1,2,3,4,5,6,7,8,9,10,11,12,13,14}]];
Print["e ", LINE["MULT",{2,3,14}]];
RENUMBER D2.5;
Print["g ", LINE["NAME",{1,2,3,4,5,6,7,8}]];
RENUMBER XXX;''')
    assert _nums(_lines(out, 'm')[0]) == [0, 1, 2, 0, 2, 3, 4]
    names = re.findall(r'"([^"]+)"', ''.join(l for l in out.splitlines() if l.startswith('d ') or l.startswith(' "')))
    assert names[:14] == ['ORG', 'QF.3', 'D2.5', 'B.3', 'D2.6', 'QD.3', 'D1.3', 'QF.4', 'D2.7', 'B.4', 'D2.8', 'QD.4', 'D1.4', 'QF.1']
    assert _nums(_lines(out, 'e')[0]) == [3, 5, 1]
    assert 'QF.4","D2.1","B.1","D2.2"' in out
    assert 'Missing origin component for RENUM_BER' in out


def test_type_k5_dodeca():
    out = _run('''DODECA D6 = (L=0.2 K5=10000);
Print[1];''')
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('DODECA D6 = (L=0.2 K5=10000);\nMARK M=();\nLINE T=(M D6);\nFFS;\nUSE T; INS; CALC;\nTYPE D6;\nPrint["K5 ", LINE["K5","D6"]];\nExit;\n', out=buf)
    o = buf.getvalue()
    assert 'K5 =10000' in o
    assert _nums(_lines(o, 'K5')[0])[-1] == 10000


def _table(out, start, n):
    ls = out.splitlines()
    i = [j for j, l in enumerate(ls) if l.strip().startswith(start)][0]
    return ls[i:i + n]


def test_disp_ga_z_reference():
    out = _run('''USE RINGS; CELL; CALC;
DISP GA 1 6;
DISP Z 1 4;
REFERENCE;
QF K1 0.75;
CALC;
Print["R ", Twiss["RBX",{5,10}], Twiss["DBX",{5,10}], Twiss["DNX",5]];
DISP REFERENCE 1 6;
DISP DREFERENCE 1 6;
DISP Z REF 1 4;''')
    ga = _table(out, 'AX      BX      NX      GMX', 7)
    assert ga[0].rstrip() == '   AX      BX      NX      GMX     GMX*L Element   Length   Value      s(m)       AY      BY      NY      GMY     GMY*L  DetR     #'
    assert ga[5][:95] == ' 1.59073 5.55659  .01495  .63536  .31768  D2.1      .50000  .5000000    .500000 -.50572 2.22702  .03313  .56387  .28194'[:95]
    z = _table(out, 'AZ      BZ      NZ', 5)
    assert z[0].rstrip() == '   AZ      BZ      NZ      DZ     DDP    Element   Length   Value      s(m)       ZX      ZPX     ZY      ZPY     GMZ    DetR     #'
    assert z[4].startswith('  .0     1.00000  .0      .0      .0      QF.1      .30000  .7000000    .200000  .0      .0      .0      .0     1.00000')
    assert _nums(_lines(out, 'R')[0])[:4] == pytest.approx([5.556585, 4.499806, -1.515e-4, -0.011829], rel=1e-3)
    r = _table(out, 'AXR', 5)
    assert r[0].rstrip() == '   AXR     BXR     NXR     EXR     EPXR  Element   Length   Value      s(m)       AYR     BYR     NYR    REY      EPYR   DetRR    #'
    assert r[4][:60] == ' -2.2343 5.34874  .00649  .64925  .21405  QF.1      .30000  .7500000    .200000'[:60]
    d = _table(out, 'dAX', 5)
    assert d[0].rstrip() == '  dAX     dBX/BX  dNX     dEX     dEPX   Element   Length   Value      s(m)      dAY     dBY/BY  dNY     dEY     dEPY   dDetR     #'
    assert d[4][:72].startswith(' -.14067 -3.4E-4 2.50E-4 -.08251 -.02012  QF.1      .30000  .7500000')
    zr = _table(out, 'AZR', 1)
    assert zr[0].rstrip() == '   AZR     BZR     NZR     DZR    DDP R  Element   Length   Value      s(m)       ZXR     ZPR     ZYR    RZPY     GMZ    DetRR    #'


def test_disp_b_beam():
    out = _run('''USE RINGS; CELL; CALC;
REFERENCE;
QF K1 0.75;
CALC;
EMIT;
DISP B 1 6;
DISP B REF 1 6;''')
    b = _table(out, 'AXp', 5)
    assert b[0].rstrip() == '   AXp     BXp    EMITXp   EXp     EPXp  Element   Sigx(mm)  Sigy(mm)  Rot(deg)   AYp     BYp    EMITYp   EYp     EPYp   sigp     #'
    assert b[4][:70] == ' -2.3750 5.34691 1.73E-7  .56675  .19393  QF.1     1.1711359 .06627814'[:70]
    assert 'REF and DREF not implemented for DISP B' in out


def test_dump_pattern_truncated_and_ddk():
    out = _run('''USE RINGS; CELL; CALC;
DELX UNIFORM 1E-4 QF.1|QF.2;
DL UNIFORM 1E-3 D1.1;
DK UNIFORM 1E-3 QF.3;
DDK UNIFORM 1E-5 B.1;
DUMP QF.1|QF.2|QF.3|QD.2|QD.3|D1.1|B.1;
DUMP QF.1|D1.1|B.1;
DDK 1E-5 QF.*;
DUMP QF.4;''')
    parts = out.split('  element ')
    first = [l.split()[0] for l in parts[1].splitlines()[1:] if l.strip() and l.split()[0][0] in 'QBD' and '.' in l.split()[0]]
    assert first == ['QF.1', 'QF.3']
    second = [l.split()[0] for l in parts[2].splitlines()[1:] if l.strip() and '.' in l.split()[0] and l.split()[0][0] in 'QBD']
    assert second == ['QF.1', 'B.1', 'D1.1']
    # QUAD に DDK は入らない (K0 のキーワードがない)。QF.4 は誤差なし: 全要素の行は出ない (パターンが合った)
    assert 'QF.4' not in parts[3]


def test_orbitgeo_end_uses_dz():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    # 閉軌道が 0 でないリング (BEND の K0 で軌道を曲げる): 周回の DZ が $$$ に入る (gs: Twiss["DZ","$$$"] = -5.338e-5 の例)
    run_sad_string(LAT.replace('B  = (L=1.2 ANGLE=0.2617993877991494)', 'B  = (L=1.2 ANGLE=0.2617993877991494 K0=1e-4)')
                   + 'FFS;\n$FORM="14.8";\nUSE RINGS; CELL; CALC;\nPrint["DZ ", Twiss["DZ",{5,"$$$"}], " ", Twiss["DX","$$$"]];\nExit;\n', out=buf)
    v = _nums(_lines(buf.getvalue(), 'DZ')[0])
    assert abs(v[1]) > 1e-9          # 一周の経路長のずれ (以前は $$$ だけ 0)
