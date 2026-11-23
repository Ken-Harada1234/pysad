# -*- coding: utf-8 -*-
"""独立な検証 (b-2): 実験室系の 6 次元ガウス分布からの 2 粒子の対のモンテカルロ (厳密な相対論の運動学)。
実験室系の 6 次元の共分散 Σ から、同じ位置の 2 粒子の運動量を作る (重心の運動量 ~ N(0, Σ_pp/2)、相対運動量 ~ N(0, 2P)、P は位置を固定したときの
実験室系の (px, py, δ) の共分散。準モンテカルロ)。各対を 4 元運動量にして対の重心系へ Lorentz 変換し、重心系で非相対論的な Møller の角分布
(4 r0²/k⁴ K(θ)、θ の対数の Gauss-Legendre × φ 一様) で散乱させ、実験室系へ戻して δ_new = |P|/p0 − 1 を厳密に求める。損失の判定は kick (δ_new − δ_old が a+ を超える、または −a− を下回る)。
流束は Møller の速度 v_M = √((P1·P2)² − 1)/(E1 E2)、密度は実験室系の位置の共分散 (4π)^{3/2} √det Σ_xx。
touschek_piwinski が使う静止系への写像 (z* = γz、ps* = δ/γ) や非相対論的な近似、q_s = 0 の近似を使わない。
"""
import numpy as np, math
from scipy.stats import qmc
R0=2.8179403262e-15; CLIGHT=299792458.0
def Kth(th):
    s2=np.sin(th)**2
    return 4/s2**2-3/s2
def lab_cond(S):
    ix=[0,2,4]; ip=[1,3,5]
    xx=S[np.ix_(ix,ix)]; xp=S[np.ix_(ix,ip)]; pp=S[np.ix_(ip,ip)]
    P=pp-xp.T@np.linalg.solve(xx,xp)
    return xx,pp,P
def mom4(px,py,dl,p0):
    P=p0*(1+dl); Px=px*p0; Py=py*p0
    Ps=np.sqrt(np.maximum(P**2-Px**2-Py**2,0))
    E=np.sqrt(1+P**2)
    return np.stack([E,Px,Py,Ps],axis=-1)
def boost_to(u, v):
    """boost 4-vectors v (...,4) into frame moving with 4-velocity u (..,4) (u0=gamma, u123=gamma*beta). returns v'."""
    g=u[...,0]; bvec=u[...,1:]/g[...,None]
    bp=np.sum(bvec*v[...,1:],axis=-1); b2=np.sum(bvec**2,axis=-1)
    gam=g
    f=(gam-1)/np.where(b2>0,b2,1)*bp - gam*v[...,0]
    E=gam*(v[...,0]-bp)
    p=v[...,1:]+f[...,None]*bvec
    return np.concatenate([E[...,None],p],axis=-1)
def boost_from(u, v):
    u2=u.copy(); u2[...,1:]*=-1
    return boost_to(u2,v)
def mc_rate(S, al, circ, N, p0, ap, am, npairs=8192, nth=40, nph=32, seed=1, chunk=1024, include_dold=False):
    h0=math.sqrt(1+p0**2)
    xx,pp,P=lab_cond(S)
    vol=(4*np.pi)**1.5*math.sqrt(np.linalg.det(xx))
    Lq=np.linalg.cholesky(2*P); Lc=np.linalg.cholesky(pp/2)
    z=qmc.MultivariateNormalQMC(mean=np.zeros(6),cov=np.eye(6),seed=seed).random(npairs)
    qrel=z[:,:3]@Lq.T; cmom=z[:,3:]@Lc.T
    p1=cmom+qrel/2; p2=cmom-qrel/2           # lab (px,py,delta)
    f1=mom4(p1[:,0],p1[:,1],p1[:,2],p0); f2=mom4(p2[:,0],p2[:,1],p2[:,2],p0)
    tot=f1+f2
    mtot=np.sqrt(tot[:,0]**2-np.sum(tot[:,1:]**2,axis=1))
    u=tot/mtot[:,None]                       # CM 4-velocity
    q1=boost_to(u,f1)                        # CM momentum of particle 1
    pm=np.linalg.norm(q1[:,1:],axis=1)       # |p*|
    k=2*pm
    uh=q1[:,1:]/pm[:,None]
    # v_Moller
    dot=f1[:,0]*f2[:,0]-np.sum(f1[:,1:]*f2[:,1:],axis=1)
    vM=np.sqrt(np.maximum(dot**2-1,0))/(f1[:,0]*f2[:,0])
    d1=np.linalg.norm(f1[:,1:],axis=1)/p0-1; d2=np.linalg.norm(f2[:,1:],axis=1)/p0-1
    # orthonormal frame about uh
    ref=np.where(np.abs(uh[:,0:1])<0.9,np.array([[1,0,0]]),np.array([[0,1,0]]))
    e1=np.cross(uh,ref); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(uh,e1)
    # theta nodes (log), phi uniform
    gx,gw=np.polynomial.legendre.leggauss(nth)
    phis=(np.arange(nph)+0.5)*2*np.pi/nph
    total=0.0; sq=0.0
    out=np.zeros(npairs)
    for c0 in range(0,npairs,chunk):
        sl=slice(c0,min(npairs,c0+chunk)); n=sl.stop-sl.start
        kk=k[sl]; 
        # theta lower bound: need possible |kick|>a ; kick <= 2 sin(th/2) * p*/... in delta: ~ 2 sin(th/2) pm*gamma/p0 ; use generous
        cmin=min(ap,am)
        tmin=cmin*p0/(h0*pm[sl])/ (1.0)   # = a beta/pm  (kick max ~ 2 sin(th/2)*pm/beta... ) -> sin(th/2)>= a*beta/(2 pm)*0.5
        sinhalf=np.minimum(0.5*cmin*(p0/h0)/ (2*pm[sl]) ,1.0)
        thlo=2*np.arcsin(np.minimum(sinhalf,1.0))
        thhi=np.full(n,np.pi/2)
        ok=thlo<thhi
        lo=np.log(np.maximum(thlo,1e-9)); hi=np.log(thhi)
        lt=0.5*(hi-lo)[:,None]*gx[None,:]+0.5*(hi+lo)[:,None]*np.ones((1,1)); th=np.exp(lt)
        wt=0.5*(hi-lo)[:,None]*gw[None,:]*th*np.sin(th)*Kth(th)*(2*np.pi/nph)   # dOmega weight K
        # n' directions in CM (n,nth,nph,3)
        st=np.sin(th)[:,:,None,None]; ct=np.cos(th)[:,:,None,None]
        nd=ct*uh[sl][:,None,None,:]+st*(np.cos(phis)[None,None,:,None]*e1[sl][:,None,None,:]+np.sin(phis)[None,None,:,None]*e2[sl][:,None,None,:])
        # new CM momenta
        q1n=np.empty(nd.shape[:-1]+(4,)); q1n[...,0]=q1[sl][:,0][:,None,None]; q1n[...,1:]=pm[sl][:,None,None,None]*nd
        q2n=q1n.copy(); q2n[...,1:]*=-1
        un=np.broadcast_to(u[sl][:,None,None,:],q1n.shape)
        l1=boost_from(un,q1n); l2=boost_from(un,q2n)
        dn1=np.linalg.norm(l1[...,1:],axis=-1)/p0-1; dn2=np.linalg.norm(l2[...,1:],axis=-1)/p0-1
        K1=dn1-d1[sl][:,None,None]; K2=dn2-d2[sl][:,None,None]
        L=((K1>ap)|(K1<-am)).astype(float)+((K2>ap)|(K2<-am)).astype(float)
        W=(L*wt[:,:,None]).sum(axis=(1,2))
        W=np.where(ok,W,0.0)
        out[sl]=vM[sl]*(4*R0**2/kk**4)*W
    val=out.mean(); err=out.std()/math.sqrt(npairs)
    # rate: Ndot = (N^2/vol)(1/2) c <v_M * sigma-weight>  ; hemisphere (L1+L2) form
    Ndot=(N**2/vol)*0.5*CLIGHT*val
    return al/circ*Ndot/N, al/circ*(N/vol)*0.5*CLIGHT*err
