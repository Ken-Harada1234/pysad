# -*- coding: utf-8 -*-
"""DDBA_01〜28 の共通部: 版の対応表 (MAP)、元の台本の読み込み、pysad での cell calc emit の実行。
元の台本のあるフォルダは環境変数 DDBA_SRC で与える (リポジトリには置かない)。"""
import os, sys, io, re
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE)); ROOT = os.path.dirname(EX)
sys.path.insert(0, ROOT)
LAT = os.path.join(EX, 'ddba', 'lattice')
# (NN, 元の台本のファイル名, 例のファイル名の中の台本名, 元の台本の相対位置のヒント(DDBA_SRC の下))
MAP = [
 (1,  'v03_12_30cells_DDBA', 'v03_12_30cells_DDBA'),
 (2,  'v05_08_30cell', 'v05_08_30cell'),
 (3,  'v05_12', 'v05_12'),
 (4,  'v05_12', 'v05_12'),
 (5,  'v05_12', 'v05_12'),
 (6,  'v05_14_2.5GeV', 'v05_14_2.5GeV'),
 (7,  'v05_14_spcod', 'v05_14_spcod'),
 (8,  'v05_14_spcod', 'v05_14_spcod'),
 (9,  'v05_14_spcod', 'v05_14_spcod'),
 (10, 'v05_14', 'v05_14'),
 (11, 'v05_14_[a-z]*_SX', 'v05_14_sx'),
 (12, 'v05_14_47.1_13.2', 'v05_14_47.1_13.2'),
 (13, 'v05_17', 'v05_17'),
 (14, 'v05_20', 'v05_20'),
 (15, 'v05_21', 'v05_21'),
 (16, 'v05_22', 'v05_22'),
 (17, 'v05_28', 'v05_28'),
 (18, 'v05_29', 'v05_29'),
 (19, 'v06_14_45.13_12.21', 'v06_14_45.13_12.21'),
 (20, 'v05_29_manyrings', 'v05_29_manyrings'),
 (21, 'v05_30', 'v05_30'),
 (22, 'v07_12', 'v07_12'),
 (23, 'v06_26', 'v06_26'),
 (24, 'v07_14_SP_transport', 'v07_14_SP_transport'),
 (25, 'v07_14_2.5GeV', 'v07_14_2.5GeV'),
 (26, 'v07_14_2.5GeV', 'v07_14_2.5GeV'),
 (27, 'v07_14_*_20231205', 'v07_14_sx_20231205'),
 (28, 'v07_16_sx_rev9_5GeV', 'v07_16_sx_rev9_5GeV'),
]
def find_src(name, base=None):
    base = base or os.environ['DDBA_SRC']
    import fnmatch
    for dp, dn, fn in os.walk(base):
        hit = [f for f in fn if fnmatch.fnmatch(f, name)]
        if hit and 'elegant' not in dp and 'csr' not in dp: return os.path.join(dp, hit[0])
    raise FileNotFoundError(name)
def lat_path(nn):
    for m in MAP:
        if m[0] == nn: return os.path.join(LAT, 'DDBA_%02d_%s.sad' % (nn, m[2]))
CALC = r'''
FFS;
$FORM = "15.10"; PageWidth = 1999;
USE RING;
cell calc emit;
Print["circumference [m] = ", LINE["S","$$$"]];
Print["number of elements = ", LINE["LENGTH"]];
Print["tune nx ny = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
Print["EMITX [m rad] = ", EMITX, "  SIGE = ", SIGE, "  SIGZ [m] = ", SIGZ];
Print["MOMENTUM [eV] = ", MOMENTUM];
nb = 0; nr = 0; nq = 0; ns = 0; no = 0; nc = 0; sa = 0; vol = 0; frq = 0; hm = 0;
Do[
  t = LINE["TYPENAME", i];
  If[t == "BEND" && LINE["L", i] > 0 && Abs[LINE["ANGLE", i]] > 0.02, nb = nb + 1; sa = sa + LINE["ANGLE", i]];
  If[t == "BEND" && LINE["L", i] > 0 && Abs[LINE["ANGLE", i]] <= 0.02 && LINE["ANGLE", i] <> 0, nr = nr + 1];
  If[t == "QUAD" && LINE["L", i] > 0 && LINE["K1", i] <> 0, nq = nq + 1];
  If[t == "SEXT" && LINE["L", i] > 0 && LINE["K2", i] <> 0, ns = ns + 1];
  If[t == "OCT", no = no + 1];
  If[t == "CAVI", nc = nc + 1; vol = vol + LINE["VOLT", i]; frq = LINE["FREQ", i]; hm = LINE["HARM", i]],
{i, 1, LINE["LENGTH"]}];
Print["bends = ", nb, "  sum of ANGLE [rad] = ", sa, "  small-angle bends = ", nr, "  quads = ", nq, "  sexts = ", ns, "  octs = ", no, "  cavities = ", nc];
Print["CAVI total VOLT [V] FREQ [Hz] HARM = ", vol, " ", frq, " ", hm];
Exit;
'''
def run(src, mode='k64', calc=CALC):
    os.environ['PYSAD_COMPAT'] = mode
    from pysad import compat
    compat.set_mode(mode)
    from pysad import ffscmd
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src + '\nON LOG;\nON RAD FLUC COD RFSW RADCOD;\n' + calc, out=buf)
    return buf.getvalue()
CHROM = r'''
FFS;
$FORM = "15.10"; PageWidth = 1999;
USE RING;
ChromFD[d_]:=(
  DP0 = d;  FFS["cell calc"];  x1 = Twiss["NX","$$$"]/2/Pi;  y1 = Twiss["NY","$$$"]/2/Pi;
  DP0 = -d; FFS["cell calc"];  x2 = Twiss["NX","$$$"]/2/Pi;  y2 = Twiss["NY","$$$"]/2/Pi;
  DP0 = 0;  FFS["cell calc"];
  {(x1 - x2)/(2 d), (y1 - y2)/(2 d)}
);
cell calc;
cs = ChromFD[0.001];
Print["chromaticity with sextupoles = ", cs[1], " ", cs[2]];
Do[If[LINE["TYPENAME", i] == "SEXT", LINE["K2", i] = 0], {i, 1, LINE["LENGTH"]}];
cell calc;
cn = ChromFD[0.001];
Print["natural chromaticity = ", cn[1], " ", cn[2]];
Exit;
'''
