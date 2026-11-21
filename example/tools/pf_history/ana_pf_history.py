# -*- coding: utf-8 -*-
"""PF の履歴ラティス 4 本の解析 (pysad)。結果を JSON に書く: py -3.12 example/tools/pf_history/ana_pf_history.py OUT.json
各ラティス・各モード (k64 と v2) を別のプロセスで走らせる (互換モードは import 時に決まるため)。"""
import os, sys, json, subprocess, re, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pf_lib
RADINT = r'''
USE RING;
cell calc emit;
bposi = LINE["POSITION","B*"];
nsl = 16;
I1 = 0; I2 = 0; I3 = 0; I4 = 0; I5 = 0; Hs = 0; Ls = 0;
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
  sa = sa*L/nsl/3;  sh = sh*L/nsl/3;
  I1 = I1 + sa/rho; I2 = I2 + L/rho^2; I3 = I3 + L/Abs[rho]^3;
  I4 = I4 + sa/rho^3 - 2*(eta0*Tan[e1] + eta1*Tan[e2])/rho^2;
  I5 = I5 + sh/Abs[rho]^3;
  Hs = Hs + sh; Ls = Ls + L;
,{i,1,Length[bposi]}];
Print["RADINT ", I1, " ", I2, " ", I3, " ", I4, " ", I5, " ", Hs/Ls, " ", Length[bposi]];
'''
def one(lat, mode):
    import build_examples as be
    body = RADINT
    ffs, txt = pf_lib.run(lat, mode, body=body)
    d = pf_lib.emit_table(txt)
    m = re.search(r'RADINT\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', txt)
    d['I1'], d['I2'], d['I3'], d['I4'], d['I5'], d['Havg'], d['nb_rad'] = map(float, m.groups())
    rows, els = ffs.rows, ffs.els
    d['C'] = rows[-1].s; d['nelem_rows'] = len(rows)
    d['nux'] = rows[-1].mux / 2 / math.pi; d['nuy'] = rows[-1].muy / 2 / math.pi
    # 要素の数
    T = lambda e: type(e).__name__
    d['nb'] = sum(1 for e in els if T(e) == 'Bend' and getattr(e, 'ANGLE', 0) != 0)
    d['nq'] = sum(1 for e in els if T(e) == 'Quadrupole' and e.L > 0 and e.K1 != 0)
    d['ns'] = sum(1 for e in els if T(e) == 'Sextupole' and e.L > 0)
    d['no'] = sum(1 for e in els if T(e) == 'Octupole')
    d['nc'] = sum(1 for e in els if T(e) == 'Cavity')
    return d
def twiss_geometry(lat):
    """k64 で 1 回走らせて、リング 1 周の Twiss と要素の並びを返す (図・セル・直線部の解析用)"""
    ffs, txt = pf_lib.run(lat, 'k64')
    return ffs
if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == '--one':
        print(json.dumps(one(sys.argv[2], sys.argv[3])))
    else:
        out = {}
        for lat, yr, nm in pf_lib.LATS:
            out[lat] = {}
            for mode in ('k64', 'v2'):
                r = subprocess.run([sys.executable, os.path.abspath(__file__), '--one', lat, mode], capture_output=True, text=True, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
                out[lat][mode] = json.loads(r.stdout.strip().splitlines()[-1])
        json.dump(out, open(sys.argv[1], 'w'), indent=1)
        for lat in out:
            for mode in out[lat]:
                d = out[lat][mode]
                print(lat, mode, 'C=%.5f nux=%.5f nuy=%.5f ex=%.4e sige=%.4e alpha=%.5f U0=%.4f tau=%.3f/%.3f/%.3f J=%.4f/%.4f/%.4f h=%.3f fs=%.0f I2=%.4f I4=%.4f I5=%.4e Havg=%.4f' % (d['C'], d['nux'], d['nuy'], d['ex'], d['sige'], d['alpha'], d['U0'], d['taux']*1e3, d['tauy']*1e3, d['tauz']*1e3, d['Jx'], d['Jy'], d['Jz'], d['h'], d['fs'], d['I2'], d['I4'], d['I5'], d['Havg']))
