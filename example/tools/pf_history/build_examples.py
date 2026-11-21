# -*- coding: utf-8 -*-
"""example/pf/optics/FFS_cell_calc_emit_PF_<年>_<名前>.sad を作る (4 本)。ラティス (example/pf/lattice/pf_*.sad) を先頭に貼り、
計算の台本をつなぐ。gs と pysad の画面表示は example/tools/run_example.py で作る (.gs.txt と .pysad.txt)。
使い方: py -3.12 example/tools/pf_history/build_examples.py [--run]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE))
OPT = os.path.join(EX, 'pf', 'optics')
# ラティスごとの違い: (ファイル名の末尾, ラティス, 直線部・標準セルの Twiss を出す MARK, 標準セルの偏向磁石)
CASES = {
 '1981_460nmrad': dict(marks=['IMARK', 'CAVI', 'MMARK', 'DMARK'], bends=['B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12']),
 '1986_128nmrad': dict(marks=['IMARK', 'CAVI', 'MMARK', 'DMARK'], bends=['B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12']),
 '1997_036nmrad': dict(marks=['U2', 'M16', 'U5', 'M13', 'V14', 'M26', 'M28', 'U19', 'CAVI'], bends=['B19', 'B20', 'B21', 'B22', 'B23', 'B24', 'B25', 'B26']),
 '2005_036nmrad': dict(marks=['U02', 'U05', 'U13', 'U14', 'U16', 'U19', 'U28', 'INJ', 'CAVI'], bends=['B20', 'B21', 'B22', 'B23', 'B24', 'B25', 'B26']),
}
BODY = r'''FFS;
$FORM = "15.10";
PageWidth = 1999;

(* 要素の数を型ごとに数える (偏向磁石は ANGLE <> 0、四極は L > 0 かつ K1 <> 0。ステアリング・スキュー四極は含めない) *)
CountEl[]:=(
  nb = 0; nq = 0; ns = 0; no = 0; nc = 0; sa = 0;
  Do[
    t = LINE["TYPENAME", i];
    If[t == "BEND" && LINE["ANGLE", i] <> 0, nb = nb + 1; sa = sa + LINE["ANGLE", i]];
    If[t == "QUAD" && LINE["L", i] > 0 && LINE["K1", i] <> 0, nq = nq + 1];
    If[t == "SEXT" && LINE["L", i] > 0, ns = ns + 1];
    If[t == "OCT", no = no + 1];
    If[t == "CAVI", nc = nc + 1],
  {i, 1, LINE["LENGTH"]}];
  Print["bends = ", nb, "  sum of ANGLE [rad] = ", sa, "  quads = ", nq, "  sexts = ", ns, "  octs = ", no, "  cavities = ", nc]
);
(* 色収差 = 運動量を ±d ずらしたチューンの差 *)
ChromFD[d_]:=(
  DP0 = d;  FFS["cell calc"];  x1 = Twiss["NX","$$$"]/2/Pi;  y1 = Twiss["NY","$$$"]/2/Pi;
  DP0 = -d; FFS["cell calc"];  x2 = Twiss["NX","$$$"]/2/Pi;  y2 = Twiss["NY","$$$"]/2/Pi;
  DP0 = 0;  FFS["cell calc"];
  {(x1 - x2)/(2 d), (y1 - y2)/(2 d)}
);
(* 場所での BX BY EX NX NY (NX NY は 2 pi で割ったチューンの位相) *)
ShowAt[nm_]:=Print[nm, "  #", LINE["POSITION", nm], "  S = ", LINE["S", nm], "  BX BY EX AX = ",
  Twiss[{"BX","BY","EX","AX"}, nm], "  NX NY/2pi = ", Twiss["NX", nm]/2/Pi, " ", Twiss["NY", nm]/2/Pi];

USE RING;
cell calc emit;
Print["circumference [m] = ", LINE["S","$$$"]];
Print["number of elements = ", LINE["LENGTH"]];
Print["tune nx ny = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
Print["EMITX [m rad] = ", EMITX, "  SIGE = ", SIGE, "  SIGZ [m] = ", SIGZ];
Print["MOMENTUM [eV] = ", MOMENTUM];
CountEl[];
Print["CAVI VOLT [V] FREQ [Hz] = ", LINE["VOLT","CAVI"], " ", LINE["FREQ","CAVI"]];
Print["--- Twiss at the marks (straight sections)"];
%MARKS%
Print["--- Twiss at the bends of the normal cells"];
%BENDS%
Print["--- chromaticity d(nu)/d(dp/p), dp/p = +-0.001"];
Print["with the design sextupoles: ", ChromFD[0.001]];
SF* 0; SD* 0;
FFS["cell calc"];
Print["natural (SF* SD* = 0):      ", ChromFD[0.001]];
Exit;
'''
def body(case):
    c = CASES[case]
    marks = '\n'.join('ShowAt["%s"];' % m for m in c['marks'])
    bends = '\n'.join('ShowAt["%s"];' % m for m in c['bends'])
    return BODY.replace('%MARKS%', marks).replace('%BENDS%', bends)
def make(case):
    lat = open(os.path.join(EX, 'pf', 'lattice', 'pf_' + case + '.sad'), encoding='utf-8').read()
    p = os.path.join(OPT, 'FFS_cell_calc_emit_PF_' + case + '.sad')
    open(p, 'w', encoding='utf-8', newline='\n').write(lat.rstrip('\n') + '\nON LOG;\nON RAD FLUC COD RFSW RADCOD;\n' + body(case))
    return p
if __name__ == '__main__':
    ps = [make(c) for c in CASES]
    print('\n'.join(ps))
    if '--run' in sys.argv:
        sys.path.insert(0, os.path.join(EX, 'tools'))
        import run_example
        for p in ps:
            run_example.run(p, 300)
