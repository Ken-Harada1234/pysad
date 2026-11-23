# -*- coding: utf-8 -*-
"""独立な検証 (b-1): 3 次元の一般の積分。touschek_piwinski の 1 次元の式 (q_s = 0、x-y ブロックの 2 固有値) の近似を外して評価する。
運動量の相対ベクトル q ~ N(0, 2P) (P は位置を固定したときの 3×3 の運動量の共分散。s 成分 σδ/γ と、すべての非対角を入れる) の方向を準モンテカルロ
(Sobol) で標本化し、動径方向の積分 ∫dk/k e^{−k²A} は k の対数の Gauss-Legendre、散乱角の半球の積分 (kick の基準 |Δδ| > a、最も近い粒子の対応)
は θ の対数の Gauss-Legendre と方位角 φ' の解析積分 (閾値を超える φ' の長さ = 2 arccos) で評価する。
  R_inst = (1/γ)(N²/V*) 2 r0² c ⟨k⁻³ Ŵ⟩,  Ŵ = ∫_半球 K(θ) [I(k|g| > 2β a+) + I(k|g| > 2β a−)] dΩ',  g = n'_s − û_s
関数 slice_rate3d(...) = 1 スライスの 1/τ への寄与 (touschek.piwinski_slice_rate と同じ規約)。
"""
import numpy as np, math
from scipy.stats import qmc
R0=2.8179403262e-15; CLIGHT=299792458.0
_gx,_gw=np.polynomial.legendre.leggauss(48)
def Kth(th):
    s2=np.sin(th)**2
    return 4/s2**2-3/s2
def What(t, ca):
    """hemisphere-kick W: int_{hemi} K [I(|g|>t)] dOmega' for arrays t, ca=cos(alpha)=u_s (same shape).  returns array"""
    t=np.asarray(t,float); ca=np.asarray(ca,float)
    sa=np.sqrt(np.maximum(1-ca**2,0))
    th_lo=2*np.arcsin(np.minimum(t/2,1.0))
    th_hi=np.full_like(t,math.pi/2)
    ok=th_lo<th_hi
    lo=np.where(ok,np.log(np.maximum(th_lo,1e-12)),0.0); hi=np.log(th_hi)
    # nodes
    x=_gx[None,:]; w=_gw[None,:]
    u=0.5*(hi-lo)[:,None]*x+0.5*(hi+lo)[:,None]
    th=np.exp(u); wt=0.5*(hi-lo)[:,None]*w*th
    st=np.sin(th); ct=np.cos(th)
    cs=ca[:,None]; ss=sa[:,None]; tt=t[:,None]
    den=ss*st
    with np.errstate(divide='ignore',invalid='ignore'):
        y1=(tt+cs*(1-ct))/den       # cos(phi') >= y1  -> g >= t
        y2=(-tt+cs*(1-ct))/den      # cos(phi') <= y2  -> g <= -t
    m1=np.where(y1>=1,0.0,np.where(y1<=-1,2*np.pi,2*np.arccos(np.clip(y1,-1,1))))
    m2=np.where(y2>=1,2*np.pi,np.where(y2<=-1,0.0,2*(np.pi-np.arccos(np.clip(y2,-1,1)))))
    # sin alpha ->0 : den=0
    small=den<1e-14
    if small.any():
        g=-cs*(1-ct)
        m1=np.where(small,np.where(g>=tt,2*np.pi,0.0),m1)
        m2=np.where(small,np.where(g<=-tt,2*np.pi,0.0),m2)
    f=st*Kth(th)*(m1+m2)*wt
    f=np.where(ok[:,None],f,0.0)
    return f.sum(axis=1)
_rx,_rw=np.polynomial.legendre.leggauss(40)
def slice_rate3d(al,P3,pl_unit_p0,p0,h0,xx,circ,N,ap,am,nu=2048,seed=1,nk=40,sigs=True):
    """P3 unused. pl_unit_p0 : 3x3 conditional cov in p0 units. generic 3D kick-criterion hemisphere integral."""
    beta=p0/h0
    P=np.asarray(pl_unit_p0)*p0**2
    if not sigs:
        P=P.copy(); P[2,:]=0; P[:,2]=0; P[2,2]=1e-30
    S=2*P
    Lc=np.linalg.cholesky(S)
    detS=np.linalg.det(S)
    sob=qmc.Sobol(3,scramble=True,seed=seed)
    z=qmc.MultivariateNormalQMC(mean=np.zeros(3),cov=np.eye(3),seed=seed).random(nu)
    q=z@Lc.T
    kk=np.linalg.norm(q,axis=1); u=q/kk[:,None]
    Sinv=np.linalg.inv(S)
    A=0.5*np.einsum('ni,ij,nj->n',u,Sinv,u)    # q-density exp(-k^2 A), A = u^T S^-1 u /2
    ca=u[:,2]
    tot=0.0
    out=np.zeros(nu)
    for c in (2*beta*ap,2*beta*am):
        # radial integral int dk/k exp(-k^2 A) What(c/k, ca)
        # k from klo=c/2 (W=0 below c/sqrt2) to kmax s.t. k^2 A=40
        kmax=np.sqrt(45.0/A)
        klo=np.full(nu,c/1.5)
        okk=kmax>klo
        lo=np.log(klo); hi=np.log(np.maximum(kmax,klo*1.0001))
        x=_rx[None,:]; w=_rw[None,:]
        lk=0.5*(hi-lo)[:,None]*x+0.5*(hi+lo)[:,None]
        k=np.exp(lk); wk=0.5*(hi-lo)[:,None]*w     # dk/k = dlnk
        t=(c/k).ravel(); cc=np.repeat(ca,nk)
        Wv=What(t,cc).reshape(nu,nk)
        out+= np.where(okk,(np.exp(-k**2*A[:,None])*Wv*wk).sum(axis=1),0.0)
    # <k^-3 What> = E_u[ 2 A^{3/2}/sqrt(pi) * int dk k^-1 e^{-k^2 A} What ]
    val=np.mean(4*A**1.5/math.sqrt(math.pi)*out)
    err=np.std(4*A**1.5/math.sqrt(math.pi)*out)/math.sqrt(nu)
    vol=(4*np.pi)**1.5*math.sqrt(abs(np.linalg.det(xx)))
    # R_inst = (1/gam)(N^2/vol)(1/2)*4 r0^2 c <k^-3 What>
    Rinst=(1/h0)*(N**2/vol)*2*R0**2*CLIGHT*val
    return al/circ*Rinst/N, al/circ*(1/h0)*(N/vol)*2*R0**2*CLIGHT*err
