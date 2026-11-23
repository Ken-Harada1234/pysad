# -*- coding: utf-8 -*-
"""担当 A2 の独立な計算 (SAD・pysad の光学・放射積分のコードを使わない)。試験 test_v2_A2_radint.py から使う。
  integrals / derived: ラティスの要素表から 3x3 行列 (x, x', δ) で I1..I5 と U0・Jx・Je・減衰時間・εx・σδ・α (平面リング)。
  coupled_partition: 自前の 6x6 行列と放射の 1 次の摂動論 (固有ベクトル) による各モードの減衰分配 (傾いた四極のある結合リング用)。
  coupled_i4: 結合リングの物理的な I4 = ∮ η_x h (h²+2k) ds + 端の項 (5x5 行列)。"""
"""自前の放射積分 (SAD・pysad の光学を使わない): ラティスの要素表から 3x3 行列 (x, x', δ) で I1..I5 を積分する。
平面リング (垂直の偏向・回転・結合なし) を仮定。6 極・八極・ステアリング・モニターはドリフト/無視。"""
import math, numpy as np
C_LIGHT=299792458.0; RE=2.8179403262e-15; MC2=0.51099895e6; LC=3.8615926796e-13
CG=4*math.pi/3*RE/(MC2*1e-9)**3
CQ=55/(32*math.sqrt(3))*LC
def drift(l): return np.array([[1,l,0],[0,1,0],[0,0,1.0]])
def quad(k,l):
    if k==0: return drift(l)
    w=math.sqrt(abs(k))
    if k>0: return np.array([[math.cos(w*l),math.sin(w*l)/w,0],[-w*math.sin(w*l),math.cos(w*l),0],[0,0,1.0]])
    return np.array([[math.cosh(w*l),math.sinh(w*l)/w,0],[w*math.sinh(w*l),math.cosh(w*l),0],[0,0,1.0]])
def bend(l,h,k):
    """曲率 h (=ANGLE/L)、勾配 k (=K1 per length、SAD では k>0 が収束) の複合機能の偏向 (3x3: x,x',δ)"""
    kx=k+h*h      # 水平の収束力 (SAD の K1 の符号の約束: 四極と同じ)。DX: x'' = -(k + h²) x + h δ
    if abs(kx)<1e-14: 
        return np.array([[1,l,h*l*l/2],[0,1,h*l],[0,0,1.0]])
    if kx>0:
        w=math.sqrt(kx); c=math.cos(w*l); s=math.sin(w*l)
        return np.array([[c,s/w,h*(1-c)/kx],[-w*s,c,h*s/w],[0,0,1.0]])
    w=math.sqrt(-kx); c=math.cosh(w*l); s=math.sinh(w*l)
    return np.array([[c,s/w,h*(1-c)/kx],[w*s,c,h*s/w],[0,0,1.0]])
def edge(h,e): return np.array([[1,0,0],[math.tan(e)*h,1,0],[0,0,1.0]])
def twiss_step(T,m):
    b,a,g=T; c,s,d,e=m[0,0],m[0,1],m[1,0],m[1,1]
    return np.array([c*c*b-2*c*s*a+s*s*g, -c*d*b+(c*e+s*d)*a-s*e*g, d*d*b-2*d*e*a+e*e*g])
def parts(els, n_per_m=40):
    """要素 -> [(種別, 行列/パラメータ)]"""
    out=[]
    for e in els:
        t=type(e).__name__; L=getattr(e,'L',0.0) or 0.0
        if t=='Bend' and L>0:
            th=getattr(e,'ANGLE',0.0)+getattr(e,'K0',0.0); h=th/L; k=getattr(e,'K1',0.0)
            geom_h=getattr(e,'ANGLE',0.0)/L          # 端の角は幾何の ANGLE に対する比 (pysad の E1 は rad)
            out.append(('bend',L,h,k,getattr(e,'E1',0.0),getattr(e,'E2',0.0)))
        elif t=='Quadrupole' and L>0: out.append(('quad',L,e.K1/1.0))
        elif L>0: out.append(('drift',L))
    return out
def integrals(els, n_per_m=40, edge_coeff=1.0):
    P=parts(els)
    C=sum((getattr(e,'L',0.0) or 0.0) for e in els)
    M=np.eye(3)
    for p in P:
        if p[0]=='drift': m=drift(p[1])
        elif p[0]=='quad': m=quad(p[2],p[1])
        else:
            _,L,h,k,e1,e2=p
            m=edge(h,e2)@bend(L,h,k)@edge(h,e1)
        M=m@M
    A=M[:2,:2]; cm=0.5*np.trace(A); mu=math.acos(cm)*(1 if A[0,1]>0 else -1)
    eta=np.linalg.solve(np.eye(2)-A,M[:2,2])
    b0=A[0,1]/math.sin(mu); a0=(A[0,0]-A[1,1])/(2*math.sin(mu)); g0=(1+a0*a0)/b0
    X=np.array([eta[0],eta[1],1.0]); T=np.array([b0,a0,g0]); I=np.zeros(5)
    H=lambda X,T: T[2]*X[0]**2+2*T[1]*X[0]*X[1]+T[0]*X[1]**2
    for p in P:
        if p[0]=='drift': m=drift(p[1]); X=m@X; T=twiss_step(T,m)
        elif p[0]=='quad': m=quad(p[2],p[1]); X=m@X; T=twiss_step(T,m)
        else:
            _,L,h,k,e1,e2=p
            if abs(h)==0 and k==0: m=drift(L); X=m@X; T=twiss_step(T,m); continue
            m=edge(h,e1); I[3]+=-edge_coeff*X[0]*math.tan(e1)*h*h; X=m@X; T=twiss_step(T,m)
            n=max(2,int(L*n_per_m)); ds=L/n; mb=bend(ds,h,k)
            for _ in range(n):
                X1=mb@X; T1=twiss_step(T,mb); xm=0.5*(X[0]+X1[0])
                I[0]+=xm*h*ds; I[1]+=h*h*ds; I[2]+=abs(h)**3*ds
                I[3]+=0.5*(X[0]*h*(h*h+2*k)+X1[0]*h*(h*h+2*k))*ds
                I[4]+=0.5*(H(X,T)+H(X1,T1))*abs(h)**3*ds
                X,T=X1,T1
            m=edge(h,e2); I[3]+=-edge_coeff*X[0]*math.tan(e2)*h*h; X=m@X; T=twiss_step(T,m)
    return I,C,mu/2/math.pi
def derived(I,C,E_eV,vrf=None,frf=None):
    E=E_eV/1e9; gam=E_eV/(MC2)
    I1,I2,I3,I4,I5=I
    U0=CG*E**4*I2/(2*math.pi)*1e9
    Jx=1-I4/I2; Je=2+I4/I2
    T0=C/C_LIGHT
    tau=[2*E*1e9*T0/(J*U0) for J in (Jx,1.0,Je)]
    ex=CQ*gam**2*I5/(Jx*I2); sd=math.sqrt(CQ*gam**2*I3/(Je*I2))
    return dict(U0=U0,Jx=Jx,Je=Je,taux=tau[0],tauy=tau[1],tauz=tau[2],ex=ex,sd=sd,alpha=I1/C)
def integrals_rich(els, n=80):
    """台形則の 2 次の誤差を Richardson で消す (n と 2n)"""
    a,C,nu=integrals(els,n_per_m=n); b,_,_=integrals(els,n_per_m=2*n)
    return b+(b-a)/3.0, C, nu


# ---------------------------------------------------------------------------------------------
# 結合リング (傾いた四極): 5x5 / 6x6 行列
# ---------------------------------------------------------------------------------------------
def _blk2(k, L):
    if k == 0:
        return np.array([[1, L], [0, 1.0]])
    w = math.sqrt(abs(k))
    if k > 0:
        return np.array([[math.cos(w * L), math.sin(w * L) / w], [-w * math.sin(w * L), math.cos(w * L)]])
    return np.array([[math.cosh(w * L), math.sinh(w * L) / w], [w * math.sinh(w * L), math.cosh(w * L)]])


def _roll(M, tilt):
    if tilt == 0.0:
        return M
    c, s = math.cos(tilt), math.sin(tilt)
    n = M.shape[0]
    R = np.eye(n)
    R[0, 0] = R[1, 1] = R[2, 2] = R[3, 3] = c
    R[0, 2] = R[1, 3] = s
    R[2, 0] = R[3, 1] = -s
    return np.linalg.inv(R) @ M @ R


def _elem(kind, L, par, tilt, n):
    """n=5: (x,x',y,y',δ)、n=6: (x,x',y,y',z,δ)。δ の列は最後"""
    M = np.eye(n)
    d = n - 1
    if kind == 'drift':
        M[0, 1] = L
        M[2, 3] = L
    elif kind == 'quad':
        M[:2, :2] = _blk2(par, L)
        M[2:4, 2:4] = _blk2(-par, L)
    else:
        h, k = par
        kx = k + h * h
        if kx == 0:
            bx = np.array([[1, L], [0, 1.0]])
            dv = np.array([h * L * L / 2, h * L])
        else:
            w = math.sqrt(abs(kx))
            if kx > 0:
                bx = np.array([[math.cos(w * L), math.sin(w * L) / w], [-w * math.sin(w * L), math.cos(w * L)]])
                dv = np.array([h * (1 - math.cos(w * L)) / kx, h * math.sin(w * L) / w])
            else:
                bx = np.array([[math.cosh(w * L), math.sinh(w * L) / w], [w * math.sinh(w * L), math.cosh(w * L)]])
                dv = np.array([h * (1 - math.cosh(w * L)) / kx, h * math.sinh(w * L) / w])
        M[:2, :2] = bx
        M[:2, d] = dv
        M[2:4, 2:4] = _blk2(-k, L)
    return _roll(M, tilt)


def _edge(h, e, n):
    M = np.eye(n)
    M[1, 0] = math.tan(e) * h
    M[3, 2] = -math.tan(e) * h
    return M


def _parts(els):
    out = []
    for x in els:
        t = type(x).__name__
        L = getattr(x, 'L', 0.0) or 0.0
        if L == 0:
            if t == 'Cavity':
                out.append(('rf', x))
            continue
        if t == 'Bend':
            out.append(('bend', L, x.ANGLE / L, x.K1, x.E1, x.E2))
        elif t == 'Quadrupole':
            out.append(('quad', L, x.K1, float((getattr(x, 'ROTATE', 0.0) or 0.0) + (getattr(x, 'ROTATION', 0.0) or 0.0))))
        else:
            out.append(('drift', L))
    return out


def coupled_i4(els, n_per_m=400):
    """結合リングの物理的な I4 = ∮ η_x h (h²+2k) ds - Σ η_x h² tan(e) (端) と I2。η_x は 4x4 の一周行列から"""
    P = _parts(els)
    Mt = np.eye(5)
    for p in P:
        if p[0] == 'rf':
            continue
        if p[0] == 'bend':
            _, L, h, k, e1, e2 = p
            m = _edge(h, e2, 5) @ _elem('bend', L, (h, k), 0.0, 5) @ _edge(h, e1, 5)
        elif p[0] == 'quad':
            m = _elem('quad', p[1], p[2], p[3], 5)
        else:
            m = _elem('drift', p[1], 0, 0.0, 5)
        Mt = m @ Mt
    D = np.linalg.solve(np.eye(4) - Mt[:4, :4], Mt[:4, 4])
    X = np.append(D, 1.0)
    I2 = I4 = 0.0
    for p in P:
        if p[0] == 'rf':
            continue
        if p[0] == 'bend':
            _, L, h, k, e1, e2 = p
            X = _edge(h, e1, 5) @ X
            I4 += -X[0] * math.tan(e1) * h * h
            n = max(2, int(L * n_per_m))
            ds = L / n
            mb = _elem('bend', ds, (h, k), 0.0, 5)
            for _ in range(n):
                X1 = mb @ X
                I2 += h * h * ds
                I4 += 0.5 * (X[0] + X1[0]) * h * (h * h + 2 * k) * ds
                X = X1
            X = _edge(h, e2, 5) @ X
            I4 += -X[0] * math.tan(e2) * h * h
        elif p[0] == 'quad':
            X = _elem('quad', p[1], p[2], p[3], 5) @ X
        else:
            X = _elem('drift', p[1], 0, 0.0, 5) @ X
    return I4, I2


def coupled_partition(els, E_eV, ndiv=40):
    """放射の 1 次の摂動論による各モードの減衰分配 (Jx 様、Jy 様、Je 様) と U0 [eV]。
    放射 (偏向磁石のスライスごと): Δδ = -ε(X)、Δ(x', y') = -ε0 (x', y')、ε = K (1+δ)² [(h+kx)² + (ky)²](1+hx) ds ≈ ε0 (1 + (2k/h + h) x + 2δ)、
    端の角: Δε = -K h² tan(e) x。RF は δ の線形キック (安定位相 sinφ = U0/V)。z 行は symplectic 条件 (R5k = -d^T S4 A_k) と R56 = -∫ h η_loc ds。
    M = M0 (1 - Σ T_i^-1 Γ_i T_i) の固有値の摂動 (左右の固有ベクトル) から |λ| の減衰を得る。戻り: [(ν, J)] を ν の大きい順、U0"""
    gam = E_eV / MC2
    Kc = (2.0 / 3.0) * RE * gam ** 3
    P = _parts(els)
    S4 = np.zeros((4, 4))
    S4[0, 1] = S4[2, 3] = 1
    S4 = S4 - S4.T
    steps = []
    I2 = 0.0

    def addz(M, hL):
        A, d = M[:4, :4], M[:4, 5]
        for k in range(4):
            M[4, k] = -float(d @ S4 @ A[:, k])
        M[4, 5] = -hL * float(d[0]) / 3.0
        return M
    for p in P:
        if p[0] == 'rf':
            steps.append(('rf', p[1]))
        elif p[0] == 'bend':
            _, L, h, k, e1, e2 = p
            for ee in (e1,):
                G = np.zeros((6, 6))
                G[5, 0] = -Kc * h * h * math.tan(ee)
                steps.append((_edge(h, ee, 6), G))
            ds = L / ndiv
            Mb = addz(_elem('bend', ds, (h, k), 0.0, 6), h * ds)
            eps0 = Kc * h * h * ds
            G = np.zeros((6, 6))
            G[5, 0] = eps0 * (2 * k / h + h)
            G[5, 5] = 2 * eps0
            G[1, 1] = eps0
            G[3, 3] = eps0
            for _ in range(ndiv):
                steps.append((Mb, G))
                I2 += h * h * ds
            G = np.zeros((6, 6))
            G[5, 0] = -Kc * h * h * math.tan(e2)
            steps.append((_edge(h, e2, 6), G))
        elif p[0] == 'quad':
            steps.append((addz(_elem('quad', p[1], p[2], p[3], 6), 0.0), np.zeros((6, 6))))
        else:
            steps.append((addz(_elem('drift', p[1], 0, 0.0, 6), 0.0), np.zeros((6, 6))))
    U0 = Kc * I2 * E_eV
    M = np.eye(6)
    Gs = np.zeros((6, 6))
    Ts = []
    for s in steps:
        if isinstance(s[0], str):
            cav = s[1]
            phis = math.asin(U0 / cav.VOLT)
            R = np.eye(6)
            R[5, 4] = -cav.VOLT * 2 * math.pi * cav.FREQ * math.cos(phis) / (E_eV * C_LIGHT)
            M = R @ M
        else:
            M = s[0] @ M
            if np.any(s[1]):
                Ts.append((M.copy(), s[1]))
    for T, G in Ts:
        Gs += np.linalg.inv(T) @ G @ T
    w, u = np.linalg.eig(M)
    wl, vl = np.linalg.eig(M.T)
    out = []
    for i in range(6):
        if w[i].imag <= 0:
            continue
        j = int(np.argmin(abs(wl - w[i])))
        vt = vl[:, j]
        dl = -(vt @ Gs @ u[:, i]) / (vt @ u[:, i])
        out.append((math.atan2(w[i].imag, w[i].real) / 2 / math.pi, 2 * E_eV * (-dl.real) / U0))
    return sorted(out, key=lambda t: -t[0]), U0
