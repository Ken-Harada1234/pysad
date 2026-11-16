"""FFS の小さなコマンド (ffsextra.py): VARIABLES ATTRIBUTE MATRIX SHOW VARY RECOVER REVERSE PRINT EXECUTE OUT/APP/CLOSE IN/READ。
参照は SAD K64 (gs) の出力 (表の書式、数値)。"""
import io

HEAD = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
USE RING;
CELL;
CALC;
'''


def _run(body, cwd=None):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(HEAD + body + 'Exit;\n', out=buf, cwd=cwd)
    return buf.getvalue()


def test_variables_show_attr_match_gs():
    out = _run('FREE QF QD;\nFIT NX 2.75 NY 1.85;\nVARIABLES;\nSHOW;\nATTR QF*;\nATTRIBUTE B;\n')
    lines = out.splitlines()
    assert '!Variable     Keyword    Now          !   Previous    Saved         Minimum      Maximum       Couple      Coefficient' in lines
    assert ' QF           K1        .700000000000 !  .700000000  .700000000  -1.00000E10   1.000000E10 <--           1.00000000' in lines
    assert ' QD           K1       -.600000000000 ! -.600000000 -.600000000  -1.00000E10   1.000000E10 <--           1.00000000' in lines
    assert '!    component1   component2   fun         goal-value     np     scale' in lines
    assert ' FIT $$$                       NX        2.750000000       1 ! * 6.283185307' in lines
    assert ' FIT $$$                       NY        1.850000000       1 ! * 6.283185307' in lines
    assert 'Element     Keyword    Value     Mimimum  Maximum  Couple         Coefficient' in lines
    assert 'QF.1        K1         .700000 -1.000E10 1.0000E10 <--           1.000000' in lines
    assert 'B.1         ANGLE      .261799 -1.000E10 1.0000E10 <--           1.000000' in lines


def test_matrix_matches_gs():
    out = _run('MATRIX;\nMAT QF.1 QD.1;\n')
    lines = out.splitlines()
    assert ' Symplectic transfer matrix from ST to $$$' in lines
    i = lines.index(' Symplectic transfer matrix from ST to $$$')
    assert lines[i + 1].split()[:2] == ['1.808320', '-4.825592'] and lines[i + 1].split()[-1] == '.633281'
    assert lines[i + 2].split()[-1] == '.071197'
    j = lines.index(' Symplectic transfer matrix from QF.1 to QD.1')
    assert lines[j + 1] == '    -.597902    2.204993     .000000     .000000     .286178 '
    assert lines[j + 3].split()[2:4] == ['2.673940', '2.726526']


def test_recover_and_print_execute():
    out = _run('FREE QF;\nFIT NX 2.75;\nGO;\nRECOVER;\nCALC;\nPRINT 1+2;\nEXECUTE "PRINT 5";\n')
    assert ' Residual = 7.3820E-02  DP = 0.01000  DP0 = 0.00000  ExponentOfResidual = 2.0  OffMomentumWeight =   1.000  CALEXP  ' in out.splitlines()
    assert 'Out[2]:= 3' in out and 'Out[3]:= 5' in out


def test_output_redirection(tmp_path):
    out = _run('Print["before"];\nOUT "o.txt";\nCALC;\nPrint["during"];\nCLOSE OUT;\nPrint["after"];\nAPP "o.txt";\nPrint["more"];\nCLOSE OUT;\n', cwd=str(tmp_path))
    assert 'before' in out and 'after' in out and 'during' not in out and 'more' not in out
    txt = (tmp_path / 'o.txt').read_text(encoding='utf-8')
    assert 'during' in txt and 'more' in txt and ' Matched.' in txt


def test_status_flag_lines_match_gs():
    out = _run('STATUS;\n')
    lines = out.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith(' CTIME='))
    assert lines[i + 1] == ' RING    CELL    RFSW    NORAD   FLUC    NORADCOD        NORADTAPER      CALC4D  KEEPEXP CALEXP  NOEMIOUT        NOCODPLOT'
    assert lines[i + 2].startswith(' COD     NOINTRA NOPOL   RADPOL  NOLWAKE NOTWAKE NOWSPAC NOSPAC  NOSELFCOD       ')
    assert lines[i + 4] == ' UNIFORM FIXSEED BIPOL   NOPSPAC ORBITCAL        CALOPT  NODAPERT        REAL    NOFOURIER       NOTRACKSIZ      SIMULATE'
    assert lines[i + 5] == ' ABSW    JITTER  TRGAUSS BARYCOD NOWAKEOPT'
    assert lines[i + 7].startswith(' EMITY=0 m       EMITZ=0 m       SIGZ=0 m        SIGE=0  XIX=0   XIY=0   NBUNCH=1        NP=1000 GCUT=1e+300')
    assert lines[i + 9] == ' FIT $$$ $$$     DISP_LAY ST $$$ ORG ST  MEA_SURE $$$'


def test_renumber_element():
    nl = chr(10)
    out = _run('RenumberElement["QF",{5,2}];' + nl + 'Print["R ",LINE["NAME",{2,6}]," ",LINE["POSITION","QF.5"]," ",LINE["POSITION","QF.2"]," ",LINE["POSITION","QF.6"]];' + nl)
    line = next(l for l in out.splitlines() if l.startswith('R '))
    assert line == 'R {"QF.5","QD.1"} 2 8 14'


def test_dump_format_matches_gs():
    """DUMP (terror.f の format 9001: f10.5 に 3p、g12.5 に 1p): DK・DDK の列は Fortran の G 形式 (0 は 0.0000、小さい値は 4.22829E-04)。参照は gs (local/wake/dm2.sad)"""
    import io
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nMARK ST = ();\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE RING = (ST 12*CELL1);\nFFS;\n$FORM="18.14";\nPageWidth=1999;\nUSE RING; CELL; CALC;\nDELX 1e-4 QF.2; DELY -2e-4 QD.3; DTHETA 3e-4 B.1; DK 1e-3 QF.5; DL 2e-4 D1.4;\n\nDUMP;\nExit;\n', out=buf)
    lines = [l.rstrip() for l in buf.getvalue().splitlines()]
    i = [k for k, l in enumerate(lines) if 'element     DL(mm)' in l][0]
    want = [' B.1          0.00000  0.0000      0.0000      -0.14289   0.00000   0.00000',
            ' QF.2         0.00000  0.0000      0.0000       0.00000   0.06172   0.00000',
            ' QD.3         0.00000  0.0000      0.0000       0.00000   0.00000  -0.33148',
            ' D1.4        -0.27985  0.0000      0.0000       0.00000   0.00000   0.00000',
            ' QF.5         0.00000 4.22829E-04  0.0000       0.00000   0.00000   0.00000']
    assert lines[i + 1:i + 6] == want


def test_matrix_accepts_bare_element_names_like_gs():
    """MATRIX [PHYSICAL] 名前 名前: 名前だけ (QF) は NAME.1。ラティスに無い名前は gs は読み飛ばし to=$$$ (esrf の disp_matrix で見つかった差)"""
    import io
    from pysad.ffscmd import run_sad_string
    src = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
PageWidth=1999;
USE RING; CELL; CALC;
MATRIX QF QD;
MATRIX PHYSICAL QF NOSUCH;
Exit;
"""
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    lines = buf.getvalue().splitlines()
    heads = [l.strip() for l in lines if 'transfer matrix from' in l]
    assert heads == ['Symplectic transfer matrix from QF.1 to QD.1', 'Physical transfer matrix from QF.1 to $$$']
    assert not any(l.startswith('?') for l in lines)
