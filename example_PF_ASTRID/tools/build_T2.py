# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkpf import mk
FLAGS = 'ON LOG;\nON RAD FLUC COD RFSW RADCOD;'
DET = 'ON LOG;\nOFF RAD FLUC;\nON COD RFSW RADCOD;'
HDR = '$FORM = "15.10";\nPageWidth = 1999;\n'
S = {}   # name -> (pre, body, lattice)

S['calc_respx_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
kick = 1E-5;
dtnx = LINE["NAME","B*"];
npm  = Length[Twiss["BX","PM*"]];
Print["number of bends B* = ", Length[dtnx], "  number of PM* = ", npm];
tunex = Twiss["NX","$$$"]/2/Pi;
bpm = Twiss["BX","PM*"];
ppm = Twiss["NX","PM*"];
js = {1,10,20};
Do[
  If[mode == 0, FFS["OFF RFSW"], FFS["ON RFSW"]];
  Print["--- RFSW mode ", mode, "  (0: RF off = pure betatron response, 1: RF on = with path-length/energy term)"];
  Do[
    j = js[jj];
    pb  = LINE["POSITION", dtnx[j]];
    bb  = Twiss["BX", pb+0.5];
    pb2 = Twiss["NX", pb+0.5];
    LINE["K0", dtnx[j]] = kick;
    FFS["cell calc;"];
    sx = Twiss["DX","PM*"];
    LINE["K0", dtnx[j]] = 0;
    FFS["cell calc;"];
    rms = 0; mx = 0;
    Do[
      an = -kick*0.5*Sqrt[bpm[i]*bb]*Cos[tunex*Pi - Abs[ppm[i]-pb2]]/Sin[tunex*Pi];
      rms = rms + (sx[i]-an)^2; mx = Max[mx, Abs[sx[i]]];
      If[i<=3, Print[dtnx[j], " PM#", i, "  SAD DX [um] ", sx[i]*1E6, "  analytic [um] ", an*1E6]];
    ,{i,1,npm}];
    Print[dtnx[j], ": max|DX| [um] ", mx*1E6, "  rms(SAD-analytic) [um] ", Sqrt[rms/npm]*1E6];
  ,{jj,1,3}];
,{mode,0,1}];
''')

S['calc_dispersion_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
pmex = Twiss["EX","PM*"];
pml  = Twiss["S","PM*"];
Do[ Print[i, " ", pml[i], " ", pmex[i]], {i,1,Length[pml]}];
(* check 1: dispersion = d(closed orbit)/d(delta), RF off so that the energy stays fixed *)
FFS["OFF RFSW"];
dd = 1E-4;
DP0 =  dd; FFS["cell calc"]; dxp = Twiss["DX","PM*"];
DP0 = -dd; FFS["cell calc"]; dxm = Twiss["DX","PM*"];
DP0 = 0;   FFS["cell calc"];
mx = 0;
Do[ mx = Max[mx, Abs[(dxp[i]-dxm[i])/(2*dd) - pmex[i]]], {i,1,Length[pml]}];
Print["max |(DX(+d)-DX(-d))/2d - EX| over PM* [m] = ", mx];
(* check 2: momentum compaction alpha = (1/C) sum over bends of theta*<eta>  (Simpson) *)
ib = LINE["POSITION","B*"];
sum = 0;
Do[
  e = Twiss["EX", {ib[i], ib[i]+0.5, ib[i]+1}];
  sum = sum + LINE["ANGLE", ib[i]] * (e[1] + 4*e[2] + e[3])/6;
,{i,1,Length[ib]}];
Print["alpha_c from bends = ", sum / LINE["S","$$$"]];
''')

S['calc_chromaticity_PF'] = (FLAGS, HDR + '''
(* natural chromaticity from quadrupoles and sextupoles (thin-lens/Simpson sums, SAD K1 and K2 are integrated strengths) *)
ChromSum[]:=(
  xq = 0; yq = 0; xs = 0; ys = 0;
  qp = LINE["POSITION","Q*"];
  Do[
    i = qp[ii];
    If[LINE["L",i] > 0 && LINE["K1",i] <> 0,
      b = Twiss[{"BX","BY"}, {i, i+0.5, i+1}];
      xq = xq - LINE["K1",i]*(b[1,1]+4*b[1,2]+b[1,3])/6/4/Pi;
      yq = yq + LINE["K1",i]*(b[2,1]+4*b[2,2]+b[2,3])/6/4/Pi],
  {ii,1,Length[qp]}];
  sp = Join[LINE["POSITION","SFH*"],LINE["POSITION","SDH*"]];
  Do[
    i = sp[ii];
    b = Twiss[{"BX","BY","EX"}, i+0.5];
    xs = xs + LINE["K2",i]*b[3]*b[1]/4/Pi;
    ys = ys - LINE["K2",i]*b[3]*b[2]/4/Pi,
  {ii,1,Length[sp]}];
  {xq,yq,xs,ys}
);
Scan1[]:=(
  tunex = Table[{0,0},{31}];
  tuney = Table[{0,0},{31}];
  Do[
    DP0 = (i - 11) * 0.001;
    FFS["cell calc"];
    tunex[i] = {DP0, Twiss["NX","$$$"]/2/Pi};
    tuney[i] = {DP0, Twiss["NY","$$$"]/2/Pi};
  ,{i,1,31}];
  DP0 = 0;
  Print["PolynomialFit x : ", PolynomialFit[tunex,6]];
  Print["PolynomialFit y : ", PolynomialFit[tuney,6]];
);
USE RING;
cell calc emit;
Print["quad sums (x,y) and sextupole sums (x,y) with design sextupoles: ", ChromSum[]];
Print["--- with design sextupoles, dp/p = -1% .. +2%"];
Scan1[];
S* 0;
FFS["cell calc"];
Print["--- natural chromaticity (S* 0)"];
Scan1[];
''')

S['calc_RadInt_PF'] = (FLAGS, HDR + '''
(* Radiation integrals of the bends (PF bends have no gradient, K1 = 0, and no vertical bend; the edge term of I4 is included).
   Twiss functions are sampled inside each bend at 17 points (Twiss["EX", n + f], f = 0 .. 1) and integrated with the composite Simpson rule. *)
USE RING;
cell calc emit;
bposi = LINE["POSITION","B*"];
nsl = 16;
I1 = 0; I2 = 0; I3 = 0; I4 = 0; I5 = 0;
Do[
  n = bposi[i];
  L = LINE["L",n]; tha = LINE["ANGLE",n]; rho = L/tha;
  sa = 0; sh = 0; eta0 = Twiss["EX", n]; eta1 = Twiss["EX", n+1];
  e1 = LINE["E1",n]*tha; e2 = LINE["E2",n]*tha;
  Do[
    t = Twiss[{"BX","AX","EX","EPX"}, n + (k-1)/nsl];
    w = If[k==1 || k==nsl+1, 1, If[Mod[k,2]==0, 4, 2]];
    sa = sa + w*t[3];
    sh = sh + w*(((1 + t[2]^2)/t[1])*t[3]^2 + 2*t[2]*t[3]*t[4] + t[1]*t[4]^2);
  ,{k,1,nsl+1}];
  sa = sa*L/nsl/3;  sh = sh*L/nsl/3;          (* integral of eta ds and of H ds over the bend *)
  I1 = I1 + sa/rho;
  I2 = I2 + L/rho^2;
  I3 = I3 + L/Abs[rho]^3;
  I4 = I4 + sa/rho^3 - 2*(eta0*Tan[e1] + eta1*Tan[e2])/rho^2;   (* edge focusing of the rectangular bend = thin quad of strength -tan(e)/rho, entering 2*k*eta/rho *)
  I5 = I5 + sh/Abs[rho]^3;
,{i,1,Length[bposi]}];
Print["I1 [m]   = ", I1];
Print["I2 [1/m] = ", I2];
Print["I3 [1/m2]= ", I3];
Print["I4 [1/m] = ", I4];
Print["I5 [1/m] = ", I5];
DFac = I4/I2;
Jx = 1 - DFac; Jy = 1; Jz = 2 + DFac;
Print["Jx Jy Jz = ", Jx, " ", Jy, " ", Jz, "   (SAD damping partition: see emit table)"];
cq = 3.832E-13;
gamma = MOMENTUM/1E6/0.51099906;
Print["gamma = ", gamma];
Print["alpha_c = I1/C = ", I1/LINE["S","$$$"]];
Print["epsilon_x = Cq gamma^2 I5/(Jx I2)  [nm rad] = ", cq*gamma*gamma*I5/I2/Jx*1E9];
Print["epsilon_x from SAD emit            [nm rad] = ", EMITX*1E9];
Print["U0 = Cgamma E^4 I2/(2 pi)            [MeV] = ", 8.846E-5*(MOMENTUM/1E9)^4*I2/2/Pi];
Print["sigma_E = sqrt(Cq gamma^2 I3/(Jz I2)) = ", Sqrt[cq*gamma*gamma*I3/(Jz*I2)], "   SAD SIGE = ", SIGE];
''')

S['FFS_FIT_injbump_PF'] = (FLAGS, HDR + '''
USE RING;
K* 0;
KK* 0;
cell calc emit;
disp $$$;
SAVE ALL;
USE INJL;
ins calc;
Print["kickers at element numbers ", LINE["POSITION","KK*"], "  INJ at ", LINE["POSITION","INJ"], "  line length ", LINE["LENGTH"]];
fit INJ DX  15/1000;
fit INJ DPX 1/1000;
fit DX  0;
fit DPX 0;
free KK1 K0;
free KK2 K0;
free KK3 K0;
free KK4 K0;
go;
disp K*;
ins calc;
Print["DX, DPX at INJ [mm, mrad] = ", Twiss["DX","INJ"]*1000, " ", Twiss["DPX","INJ"]*1000];
Print["DX, DPX at the end of the line [mm, mrad] = ", Twiss["DX","$$$"]*1000, " ", Twiss["DPX","$$$"]*1000];
kn = {"KK1","KK2","KK3","KK4"};
Do[ Print[kn[k], " K0 [mrad] = ", LINE["K0", kn[k]]*1000], {k,1,4}];
(* analytic check: a kick theta at kicker k moves the orbit downstream by sqrt(beta_k beta_s) sin(psi_s - psi_k) theta *)
bi = Twiss["BX","INJ"]; pi0 = Twiss["NX","INJ"]; be = Twiss["BX","$$$"]; pe = Twiss["NX","$$$"];
x1 = 0; x2 = 0;
Do[
  th = LINE["K0", kn[k]]; kp = LINE["POSITION", kn[k]];
  bk = Twiss["BX", kp+0.5]; pk = Twiss["NX", kp+0.5];
  If[kp < LINE["POSITION","INJ"], x1 = x1 - th*Sqrt[bk*bi]*Sin[pi0-pk]];
  x2 = x2 - th*Sqrt[bk*be]*Sin[pe-pk];
, {k,1,4}];
Print["analytic DX at INJ [mm] (kickers upstream of INJ only) = ", x1*1000, "   at line end [mm] = ", x2*1000];
''', 'v9_16_inj_noSF')

S['FFS_mom_dep_ts_PF'] = (FLAGS, HDR + '''
USE RING;
cell calc emit;
disp $$$;
DP = 0;
Print["--- with the design sextupoles"];
Do[
  DP0 = (i - 11) * 0.002;
  FFS["cell calc"];
  Print[DP0, " ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
,{i,1,21}];
DP = 0;
S* 0;
Print["--- without sextupoles (S* 0)"];
Do[
  DP0 = (i - 11) * 0.002;
  FFS["cell calc"];
  Print[DP0, " ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
,{i,1,21}];
''')

S['FFS_amp_dep_ts_PF'] = (DET, HDR + '''
SetInitialParameters[]:=(
  FFS["cell calc emit fix"];
  alphax = Twiss["AX",1]; betax = Twiss["BX",1];
  alphay = Twiss["AY",1]; betay = Twiss["BY",1];
  w = EMITX;
);
CalcAmpDepTS[inidx_,inidpx_,inidy_,inidpy_]:=(
  ininame = LINE["NAME",1];
  LINE["DX",ininame] = inidx;  LINE["DPX",ininame] = inidpx;
  LINE["DY",ininame] = inidy;  LINE["DPY",ininame] = inidpy;
  FFS["ins calc"];
  {Twiss["NX","$$$"]/2/Pi, Twiss["NY","$$$"]/2/Pi}
);
(* nx, ny: amplitude in sigma (sigma = sqrt(EMITX beta)); average over thetadiv points of the Courant-Snyder ellipse *)
AmplitudeDependentTuneShift[nx_, ny_, thetadiv_]:=(
  FFS["cell calc emit"];
  dth = 2*Pi/thetadiv;
  deptune = {0,0};
  Do[
    theta = (ith-1)*dth;
    inidx  = nx*Sqrt[w]*Sqrt[betax]*Cos[theta];
    inidpx = nx*Sqrt[w]/Sqrt[betax]*(Sin[theta] - alphax*Cos[theta]);
    inidy  = ny*Sqrt[w]*Sqrt[betay]*Cos[theta];
    inidpy = ny*Sqrt[w]/Sqrt[betay]*(Sin[theta] - alphay*Cos[theta]);
    deptune = deptune + CalcAmpDepTS[inidx,inidpx,inidy,inidpy];
  ,{ith,1,thetadiv}];
  deptune/thetadiv
);
USE RING;
cell calc emit;
DP = 0;
Print["design tune ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
Print["--- horizontal amplitude nx sigma (y = tiny), 12 points on the ellipse"];
Do[
  SetInitialParameters[];
  res = AmplitudeDependentTuneShift[nn, 0.01, 12];
  Print[nn, " ", res[1], " ", res[2]];
,{nn,5,35,5}];
Print["--- vertical amplitude ny sigma (x = tiny)"];
Do[
  SetInitialParameters[];
  res = AmplitudeDependentTuneShift[0.01, nn, 12];
  Print[nn, " ", res[1], " ", res[2]];
,{nn,5,30,5}];
''')

S['phase_PF'] = (DET, HDR + '''
(* one-turn map at the start of the ring: tracking one particle for nturn turns at several initial x amplitudes *)
USE RING;
cell calc emit;
alphax = Twiss["AX",1]; betax = Twiss["BX",1]; gammax = (1 + alphax^2)/betax;
sigmax = Sqrt[EMITX/gammax];
Print["alpha beta gamma at start = ", alphax, " ", betax, " ", gammax, "  EMITX = ", EMITX];
nturn = 60;
ka = {1,10,20,30};
Do[
  k = ka[kk];
  x1 = {Sqrt[EMITX*betax]*k};
  res = {1,{x1,{1E-9},{1E-9},{1E-9},{1E-9},{1E-9},{1}}};
  Do[
    res = TrackParticles[res,1];
    x = res[2,1,1]; px = res[2,2,1];
    Print["amp ", k, " sigma  turn ", i, "  x [mm] ", x*1000, "  px [mrad] ", px*1000, "  y [mm] ", res[2,3,1]*1000, "  z [mm] ", res[2,5,1]*1000, "  CSinv/EMITX ", (gammax*x^2 + 2*alphax*x*px + betax*px^2)/EMITX];
  ,{i,1,nturn}];
,{kk,1,4}];
''')

S['fn_DynamicApertureSurvey_DAsv_PF'] = (DET, HDR + '''
NPARA = 2;
USE RING;
cell calc emit;
DAPWIDTH = 51;
FFS["cell calc emit"];
dp = 0.01 / SIGE;
w = EMITX;
sigx = Sqrt[Twiss["BX",1]*w];  sigy = Sqrt[Twiss["BY",1]*w];
zrange = Range[-2.0*dp, 2.01*dp, dp/2];
nturn = 2000;
Print["sigma_x, sigma_y at the start [mm] = ", sigx*1000, " ", sigy*1000, "   1% = ", dp, " sigma_E,  turns = ", nturn];
range = {{1,100},{1,1},zrange};
FFS["cell calc"];
res = DynamicApertureSurvey[range,nturn,Output->6];
Print["x survey, score per dp/p row:"];
Print[Table[res[2,2,i,2], {i,1,Length[zrange]}]];
range = {{1,1},{1,100},zrange};
FFS["cell calc"];
res = DynamicApertureSurvey[range,nturn,Output->6];
Print["y survey, score per dp/p row:"];
Print[Table[res[2,2,i,2], {i,1,Length[zrange]}]];
''')

S['fn_DynamicApertureSurvey_da_xy_PF'] = (DET, HDR + '''
NPARA = 2;
USE RING;
DAPWIDTH = 51;
cell calc emit;
rmin = 1; rmax = 100; nstep = 24; nturn = 1000;
dp = 0.01 / SIGE;
zrange = {-dp*2, 0, dp*2};
nzrange = 3;
score = Table[List[], {nzrange}];
z = Table[0, {nzrange}];
Do[
  ang = (2*Pi)*((i+0.5)/nstep);
  xmax = rmax*Cos[ang]; ymax = rmax*Sin[ang];
  xmin = rmin*Cos[ang]; ymin = rmin*Sin[ang];
  range = {{xmin,xmax},{ymin,ymax},zrange};
  res = DynamicApertureSurvey[range,nturn,Output->6];
  Do[
    z[j] = res[2,2,j,1];
    s = res[2,2,j,2];
    score[j] = Append[score[j], {(s-1)/50*xmax, (s-1)/50*ymax}];
  ,{j,1,nzrange}];
,{i,0,nstep-1}];
Do[
  Print["dp/p = ", z[j], " (in sigma_E), boundary (x,y) in sigma units, ", nstep, " rays"];
  Do[ Print[score[j,i,1], " ", score[j,i,2]], {i,1,nstep}];
,{j,1,nzrange}];
''')

S['FFS_draw_topdrawer_PF'] = (FLAGS, HDR + '''
(* A TopDrawer file written by hand, i.e. the kind of text that  FFS["draw bx & by & ex b*",12]  puts in fort.12.
   Print[...] sends it to the screen (unit 6); with  Write[12, ...]  it would go to the file fort.12. *)
USE RING;
cell calc emit;
ne = LINE["LENGTH"];
smax = LINE["S","$$$"];
Print["TDBEGIN"];
Print["NEWFRAME;SET CARD 80;SET FONT DUPLEX"];
Print["SET LIMIT X 0 ", smax];
Print["TITLE 7.3 9.3 CENTER SIZE 2.5 'PF 2.5 GeV  sqrt(beta) and eta'"];
Print["SET ORDER X Y"];
Do[ Print[Twiss["S",i], " ", Sqrt[Twiss["BX",i]]], {i,1,ne}];
Print["JOIN 1 TEXT 'sqrt(BX)'"];
Print["NEW PLOT"];
Do[ Print[Twiss["S",i], " ", Sqrt[Twiss["BY",i]]], {i,1,ne}];
Print["JOIN 1 TEXT 'sqrt(BY)'"];
Print["NEW PLOT"];
Do[ Print[Twiss["S",i], " ", Twiss["EX",i]], {i,1,ne}];
Print["JOIN 1 TEXT 'EX'"];
Print["TDEND"];
''')

if __name__ == '__main__':
    out = sys.argv[1]
    names = sys.argv[2:] or list(S)
    for n in names:
        v = S[n]; pre, body = v[0], v[1]; lat = v[2] if len(v) > 2 else 'v9_16_asym'
        mk(out + '/' + n, body, lat, pre)
