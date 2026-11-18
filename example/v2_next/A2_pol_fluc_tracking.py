# -*- coding: utf-8 -*-
"""v2next A2 項目 2(b): v2 の Polarization の脱偏極率を、FLUC つきのスピン追跡 (track6、粒子ごとの Python 追跡) で測って比べる。
使い方: py -3.12 example/v2_next/A2_pol_fluc_tracking.py E[GeV] V[MV] 粒子数 乱数の種 n1 ndamp  (例: 4.2 12 120 11 1500 4500。1 粒子 1 周 約 4.5 ms)
手順: 12 セルの FODO (四極の縦ずれ ±3 mm、閉軌道 y 約 6 mm、RF 調和数 64) の閉軌道・n0・v2 の depol_rate を Emittance[] から得る
  -> 高速追跡 (FLUC、スピンなし) で 8000 周かけて平衡にし、スピンを n(u) = n0 + N0 u で置く
  -> FLUC つきのスピン追跡を n1 周 -> FLUC を切って軌道を閉軌道へ減衰 (ndamp 周) -> 1 − <S·n0> (不可逆な脱偏極) を測る (これを 2 回繰り返し、差から率)。
結果は dep_<E>_<種>.json。記録は docs/history/v2next_A2.md §3.2。"""
import sys, io, json, time, math
import sys, io, math, time, numpy as np
from pysad import compat, emit6d, track6, ffscmd
compat.set_mode('v2')
def lat(E, volt, dy=3e-3):
    return f'''MOMENTUM = {E}GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY={dy}) QDE = (L=0.3 K1=-0.6 DY=-{dy});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={volt}e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''

# 使い方: py dep_run.py E volt N seed n1 ndamp  -> 結果を dep_<E>_<seed>.json に
import numpy as np
from pysad.track6n import track_particles_fast
E=float(sys.argv[1]); volt=float(sys.argv[2]); N=int(sys.argv[3]); seed=int(sys.argv[4]); n1=int(sys.argv[5]); nd=int(sys.argv[6])
ffs,_=ffscmd.run_sad_string(lat(E,volt),out=io.StringIO())
p=ffs.momentum()
e=emit6d.emittance6d(ffs.lat,p,rfsw=True,radcod=True,pol=True)
sp=e['spin']; els=ffs.lat.elements; circ=emit6d._emit_circ(els); trf0=e.get('trf0',0.0)
cod=np.array(e['orbit_entrance'])
M=np.zeros((3,3))
for j in range(3):
    sj=np.zeros(3); sj[j]=1
    _o,s_=track6.track_particles(els,np.array([cod]),1,p,trf0,True,True,False,spins=[sj],radpol=False,circ=circ)
    M[:,j]=s_[0]
w,v=np.linalg.eig(M); k=np.argmin(abs(w-1)); n0=np.real(v[:,k]); n0/=np.linalg.norm(n0)
if n0[1]<0: n0=-n0
N0=np.array(sp['spin_orbit']) if sp.get('spin_orbit') is not None else None
res=dict(E=E,volt=volt,N=N,seed=seed,n1=n1,nd=nd,rate_model=sp['depol_rate'],spread_model=sp['spin_spread'],n0=n0.tolist(),nup=sp['nup'])
t0=time.time()
c=np.tile(cod,(N,1))
c=track_particles_fast(els,c,8000,p,trf0,True,True,True,seed=seed,circ=circ)   # 平衡まで (スピンなし、高速)
S=np.tile(n0,(N,1))
if N0 is not None:
    S=S+(c-cod)@N0.T      # n(u) = n0 + N0 u (入口で)
    S/=np.linalg.norm(S,axis=1)[:,None]
res['t_equil']=time.time()-t0
out=[]
for stage in range(2):
    c,S1=track6.track_particles(els,c,n1,p,trf0,True,True,True,spins=list(S),radpol=False,circ=circ)
    S=np.array(S1)
    cd=c.copy()
    # 減衰 (FLUC なし)
    cd,Sd=track6.track_particles(els,cd,nd,p,trf0,True,True,False,spins=list(S),radpol=False,circ=circ)
    Sd=np.array(Sd)
    sn=Sd@n0
    out.append(dict(turns=(stage+1)*n1,D=float(1-sn.mean()),err=float(sn.std()/np.sqrt(N)),orb=float(np.abs(cd-cod).max())))
    res['stages']=out
    json.dump(res,open(f'dep_{E}_{seed}.json','w'))
    print(out[-1],time.time()-t0,flush=True)
