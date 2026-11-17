# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkpf import mk
FLAGS = 'ON LOG;\nON RAD FLUC COD RFSW RADCOD;'
DET = 'ON LOG;\nOFF RAD FLUC;\nON COD RFSW RADCOD;'   # tracking examples: deterministic
HDR = '$FORM = "15.10";\nPageWidth = 1999;\n'
S = {}

S['FFS_cell_calc_emit_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
SAVE ALL;
Print["circumference [m] = ", LINE["S","$$$"]];
Print["number of elements = ", LINE["LENGTH"]];
Print["tune nx ny = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
Print["EMITX [m rad] = ", EMITX, "  SIGE = ", SIGE, "  SIGZ [m] = ", SIGZ];
Print["MOMENTUM [eV] = ", MOMENTUM];
''')

S['FFS_disp_type_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
disp QF;
type QF;
type B02;
type SFH;
disp INJ;
''')

S['fn_Twiss_LINE_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
Print["Twiss BX at QF       : ", Twiss["BX","QF"]];
Print["Twiss BX,BY at INJ   : ", Twiss[{"BX","BY"},"INJ"]];
Print["number of QF*        : ", Length[Twiss["BX","QF*"]]];
Print["QF entrance/mid/exit : ", Twiss[{"BX","BY"},{LINE["POSITION","QF"],LINE["POSITION","QF"]+0.5,LINE["POSITION","QF"]+1}]];
Print["LINE items for QF    : ", LINE[{"POSITION","NAME","ELEMENT","S","L","K1"},"QF"]];
Print["GEO of last element  : ", LINE["GEO",LINE["LENGTH"]]];
Print["LINE[LENGTH]         : ", LINE["LENGTH"]];
Print["set K1 of QF to 0.75 : "];
LINE["K1","QF"] = 0.75;
Print[LINE["K1","QF"]];
''')

S['fn_PolynomialFit_PF'] = (FLAGS, HDR + '''
tunex = Table[{0,0},{41}];
tuney = Table[{0,0},{41}];
USE RING;
cell calc emit;
Do[
  DP0 = (i - 21) * 0.001;
  FFS["cell calc"];
  tunex[i] = {DP0, Twiss["NX","$$$"]/2/Pi};
  tuney[i] = {DP0, Twiss["NY","$$$"]/2/Pi};
  Print[i, " ", DP0, " ", tunex[i,2], " ", tuney[i,2]];
,{i,1,41,1}];
DP0 = 0;
Print["symmetric range  (-2%..+2%, mean 0)"];
Print[PolynomialFit[tunex,3]];
Print[PolynomialFit[tuney,3]];
Print["shifted range (-2%..+1.0%, mean not 0)"];
tx = Take[tunex,{1,31}];
ty = Take[tuney,{1,31}];
Print[PolynomialFit[tx,3]];
Print[PolynomialFit[ty,3]];
''')

S['FFS_FIT_chrom_PF'] = (FLAGS, HDR + '''
ChromFD[]:=(
  DP0 = 1E-4;  FFS["cell calc"]; np = Twiss[{"NX","NY"},"$$$"]/2/Pi;
  DP0 = -1E-4; FFS["cell calc"]; nm = Twiss[{"NX","NY"},"$$$"]/2/Pi;
  DP0 = 0;     FFS["cell calc"];
  (np - nm)/2E-4
);
USE RING;
cell calc emit;
{nx0, ny0} = Twiss[{"NX","NY"},"$$$"]/2/Pi;
Print["design tune ", nx0, " ", ny0];
Print["xi with design sextupoles : ", ChromFD[]];
type SFH; type SDH;
S* 0;
FFS["cell calc"];
Print["xi natural (S* 0)         : ", ChromFD[]];
fit nx nx0 2;
fit ny ny0 2;
DP=0.0001;
free S*;
go;
FITP 1;
type SFH; type SDH;
SAVE S*;
cell calc;
Print["tune after fit        : ", Twiss[{"NX","NY"},"$$$"]/2/Pi];
Print["xi after fit (1e-4)   : ", ChromFD[]];
''')

S['FFS_FIT_modules_PF'] = (FLAGS, HDR + '''
ChromCorrect[]:=(
  Print["Chromaticity Correction"];
  FFS["cell calc;"];
  {nx0, ny0} = Twiss[{"NX","NY"},"***"]/2/Pi;
  FFS["FIT;"];
  FFS["NX nx0 2;"];
  FFS["NY ny0 2;"];
  FFS["DP=0.0001;"];
  FFS["FREE SFH SDH;"];
  FFS["GO;"];
  FFS["FITP 1;"];
);
MatchingNcellTune[]:=(
  Print["Matching Tune to ", tunex, " ", tuney];
  FFS["cell calc;"];
  FFS["fit NX tunex;"];
  FFS["fit NY tuney;"];
  FFS["free QD QF;"];
  FFS["go;"];
  FFS["cell calc"];
);
USE NCELL;
cell calc;
Print["one normal cell: tune ", Twiss[{"NX","NY"},"$$$"]/2/Pi, " QF K1 ", LINE["K1","QF"], " QD K1 ", LINE["K1","QD"]];
Do[
  Do[
    FFS["USE NCELL;"];
    MatchingNcellTune[];
    ChromCorrect[];
    Print["target ", tunex, " ", tuney, " reached ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi,
          " K1 QF,QD ", LINE["K1","QF"], " ", LINE["K1","QD"],
          " K2 SFH,SDH ", LINE["K2","SFH"], " ", LINE["K2","SDH"]];
  , {tuney, 0.15, 0.20, 0.05}];
, {tunex, 0.30, 0.35, 0.05}];
''')

S['fn_DynamicApertureSurvey_PF'] = (DET, HDR + '''
NPARA = 2;
USE RING;
cell calc emit;
w = EMITX;
sigx = Sqrt[Twiss["BX",1] * w];
sigy = Sqrt[Twiss["BY",1] * w];
dp = 0.01 / SIGE;
Print["sigx sigy [mm] ", sigx*1000, " ", sigy*1000, "   1% of dp/p = ", dp, " sigma_E"];
DAPWIDTH = 51;
FFS["cell calc emit"];
zrange = Range[-1.5 * dp, 1.6 * dp, dp/2];
range = {{1,100},{1,1},zrange};
resx = DynamicApertureSurvey[range,1000,Output->6];
Do[Print[zrange[i]," ",resx[2,2,i,2]," ",(resx[2,2,i,2] - 1) * sigx * 1000 * 2],{i,1,Length[zrange]}];
range = {{1,1},{1,100},zrange};
resy = DynamicApertureSurvey[range,1000,Output->6];
Do[Print[zrange[i]," ",resy[2,2,i,2]," ",(resy[2,2,i,2] - 1) * sigy * 1000 * 2],{i,1,Length[zrange]}];
''')

S['fn_TrackParticles_PF'] = (DET, HDR + '''
USE RING;
cell calc;
iq = LINE["POSITION","QF"];
rqk = Sqrt[LINE["K1","QF"]/LINE["L","QF"]];
ql  = LINE["L","QF"];
Print["QF element number ", iq, "  sqrt(K) = ", rqk, "  L = ", ql];
qkmath = {{Cos[rqk*ql], Sin[rqk*ql]/rqk},{-rqk*Sin[rqk*ql], Cos[rqk*ql]}};
qkmatv = {{Cosh[rqk*ql], Sinh[rqk*ql]/rqk},{rqk*Sinh[rqk*ql], Cosh[rqk*ql]}};
Print["hand matrix H: ", qkmath];
Print["hand matrix V: ", qkmatv];
Print["det H, det V : ", Det[qkmath], " ", Det[qkmatv]];
a = 1E-6;
Do[
  c = {{0},{0},{0},{0},{0},{0},{1}};
  c[[k]] = {a};
  after = TrackParticles[{iq, c}, iq+1];
  Print["unit ", k, " -> x,px,y,py = ", after[2,1,1]/a, " ", after[2,2,1]/a, " ", after[2,3,1]/a, " ", after[2,4,1]/a];
,{k,1,4}];
before = {iq,{{0.001},{0},{0},{0},{0},{0},{1}}};
after = TrackParticles[before, iq+1];
Print["1 mm, 0 mrad -> ", after[2]*1000];
''')

S['FFS_draw_OpticsPlot_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
(* 1. the draw command (needs a display or TopDrawer; in this run it only produces a window/file where one exists) *)
draw bx & by & ex b*;
(* 2. OpticsPlot: with DisplayFunction->Identity it returns the Graphics expression and opens no window.
      On a desktop, w = KBMainFrame[...]; c1 = Canvas[...]; $DisplayFunction = CanvasDrawer; Canvas$Widget = c1; OpticsPlot[...]; TkWait[]; gives the window. *)
g = OpticsPlot[{{"BX","BY"},{"EX"}}, DisplayFunction->Identity];
Print["OpticsPlot returned: ", Head[g], " with ", Length[g], " parts"];
(* 3. the data behind such a plot, element by element: s, BX, BY, EX *)
ne = LINE["LENGTH"];
Print["PLOTDATA ", ne];
Do[ Print[Twiss["S",i], " ", Twiss["BX",i], " ", Twiss["BY",i], " ", Twiss["EX",i], " ", LINE["NAME",i]], {i,1,ne}];
''')

if __name__ == '__main__':
    names = sys.argv[1:] or list(S)
    for n in names:
        pre, body = S[n]
        mk('T1_intro/' + n, body, pre=pre)
