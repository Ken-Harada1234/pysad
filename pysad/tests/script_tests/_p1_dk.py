# -*- coding: utf-8 -*-
"""担当 P1 の独立な計算: Derbenev-Kondratenko (DK) の式を、閉軌道に沿った数値積分で直接評価する (Emittance[] の POL の量 (tradke・gintd・SLIM の N) を使わない)。

    1/τ_ST = (5√3/8) r_e ƛ_C c γ⁵ ⟨ρ⁻³ (1 − 2/9 (n·v)²)⟩                      (Sokolov-Ternov の厳密式。反転の率)
    P_eq  = −(8/5√3) ⟨ρ⁻³ b·(n − ∂n/∂δ)⟩ / ⟨ρ⁻³ (1 − 2/9 (n·v)² + 11/18 |∂n/∂δ|²)⟩    (⟨⟩ は閉軌道の路程の平均。b = v×v̇/|v×v̇|)

使うもの: (1) 要素を細かく刻んだ配列 (BEND・QUAD は 1 刻み ≤ ds_max)、(2) 追跡 (track6.track_particles、RAD ON・損失は emit6d.CAVE の 1e-9 倍) の軌道とスピン
(Thomas-BMT を E・C2 が検証済み)、(3) 要素の磁場の定義 B/Bρ = (k y, h + k x, 0) (要素の系。ROTATE と DX・DY で ring の系に直す) から ρ⁻¹ = |B⊥|・b = B⊥/|B⊥|、
(4) スピン軸の場 n(u) = n0 + N u の N (3×6) を、要素 (刻み) ごとの軌道の写像 T_i・スピンの回転 R_i・スピン軌道結合 D_i (どれも追跡の有限差分) の鎖
N_{i+1} T_i = R_i N_i + [n_{i+1}]× D_i の周期解として求め (SLIM と同じ線形の式だが、要素ごとの行列は追跡から取る)、∂n/∂δ (位置・角度を固定した δ の偏微分
= 放射の跳びが n を動かす量) は N の δ 列 (刻みの位置で)。放射積分・K64 の tradke・gintd・Emittance[] の N は使わない。
n は ring の局所の系 (追跡のスピンの系) で、v̂ = (px, py, pz)/|p| も同じ系。
限界: SEXT・MULT の磁場は無視 (軌道が軸の近くの小さなリング)、DROTATE なし、BEND の F1・FB1・FB2 なし。SOL の領域 (BOUND の組) は 1 つの塊として追跡し、中の磁場は縦磁場 Bz だけ (B⊥ = Bz×軌道の傾き)。"""
import copy
import math

import numpy as np

from pysad import emit6d, track6

RE = 2.8179403262e-15
LC = 3.8615926796e-13
C = 299792458.0
ME = 0.51099895e6
PSTC = 8.0 / (5.0 * math.sqrt(3.0))
_CAVE0 = emit6d.CAVE


def restore_cave():
    emit6d.CAVE = _CAVE0


def slice_elements(els, ds_max=0.05):
    """BEND (L>0) を等分に刻む (QUAD は刻まない: 端のフリンジが刻みごとに入るので。四極の ρ⁻³ の寄与は偏向磁石の 1e-6 倍以下)。戻り: 塊の列 [(元の要素, [刻んだ要素の列])]。BEND は E1 を最初、E2 を最後の刻みだけに付ける"""
    out = []
    for e in els:
        t = type(e).__name__
        L = getattr(e, 'L', 0.0) or 0.0
        if t == 'Bend' and L > 0.0:
            n = max(1, int(math.ceil(L / ds_max - 1e-9)))
            pieces = []
            for i in range(n):
                s = copy.copy(e)
                s.L = L / n
                if t == 'Bend':
                    s.ANGLE = e.ANGLE / n
                    if getattr(e, 'K0', 0.0):
                        s.K0 = e.K0 / n
                    s.E1 = e.E1 if i == 0 else 0.0
                    s.E2 = e.E2 if i == n - 1 else 0.0
                pieces.append(s)
            out.append((e, pieces))
        else:
            out.append((e, [e]))
    return out


def _el_field(e, xr, yr):
    """要素の系の磁場 (Bx, By)/Bρ0 (δ=0)。xr, yr: その要素の系での軌道の位置"""
    t = type(e).__name__
    L = e.L
    if t == 'Bend' and L > 0.0:
        h = (e.ANGLE + getattr(e, 'K0', 0.0)) / L
        k = e.K1
        return k * yr, h + k * xr
    if t == 'Quadrupole' and L > 0.0:
        k = e.K1
        return k * yr, k * xr
    return 0.0, 0.0


def _to_el(e, x, y):
    """ring の系 (x, y) → 要素の系 (tchge の入口: 回して (θ = ROTATE)、ずらす)"""
    th = getattr(e, 'ROTATE', 0.0) or 0.0
    dx0 = getattr(e, 'DX', 0.0) or 0.0
    dy0 = getattr(e, 'DY', 0.0) or 0.0
    ct, st = math.cos(th), math.sin(th)
    xe = ct * x - st * y - (ct * dx0 - st * dy0)
    ye = st * x + ct * y - (st * dx0 + ct * dy0)
    return xe, ye, ct, st


def _axial(a):
    return 0.5 * np.array([a[2, 1] - a[1, 2], a[0, 2] - a[2, 0], a[1, 0] - a[0, 1]])


def _cross_mat(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


class DK:
    """els: RF つきの要素の列。p_eV: 運動量 [eV]。trf0: RF の同期位置 (放射の損失を消すので 0)"""
    H = np.array([1e-6, 1e-6, 1e-6, 1e-6, 1e-5, 1e-6])

    def __init__(self, els, p_eV, ds_max=0.05, trf0=0.0):
        self.els = list(els)
        self.p = p_eV
        self.trf0 = trf0
        self.circ = sum(getattr(e, 'L', 0.0) or 0.0 for e in els)
        self.groups = self._groups(ds_max)
        emit6d.CAVE = _CAVE0 * 1e-9

    def _groups(self, ds_max):
        sl = slice_elements(self.els, ds_max)
        groups = []
        i = 0
        while i < len(sl):
            e = sl[i][0]
            if type(e).__name__ == 'Solenoid' and getattr(e, 'BOUND', 0):
                j = i + 1
                while j < len(sl) and not (type(sl[j][0]).__name__ == 'Solenoid' and getattr(sl[j][0], 'BOUND', 0)):
                    j += 1
                grp = []
                for k in range(i, min(j + 1, len(sl))):
                    grp += sl[k][1]
                groups.append(('region', grp))
                i = j + 1
            else:
                for piece in sl[i][1]:
                    groups.append(('el', [piece]))
                i += 1
        return groups

    def _track(self, els, c, s):
        out, s1 = track6.track_particles(els, np.array(c, float), 1, self.p, self.trf0, True, True, False,
                                         spins=[np.array(x, float) for x in s], radpol=False, circ=self.circ)
        return np.array(out), np.array(s1)

    def closed_orbit(self):
        c = np.zeros(6)
        for it in range(12):
            cs = [c.copy()]
            for j in range(6):
                for sgn in (1, -1):
                    cc = c.copy()
                    cc[j] += sgn * self.H[j]
                    cs.append(cc)
            out, _ = self._track(self.els, cs, [[0, 1.0, 0]] * len(cs))
            f = out[0] - c
            J = np.zeros((6, 6))
            for j in range(6):
                J[:, j] = (out[1 + 2 * j] - out[2 + 2 * j]) / (2 * self.H[j])
            dc = np.linalg.solve(J - np.eye(6), -f)
            c = c + dc
            if np.abs(dc).max() < 1e-14:
                break
        self.J = J
        return c

    def maps(self, grp, c):
        """塊の写像 (T, R, D) と出口の軌道: 追跡の有限差分"""
        cs = [c.copy()]
        for j in range(6):
            for sgn in (1, -1):
                cc = c.copy()
                cc[j] += sgn * self.H[j]
                cs.append(cc)
        parts = []
        sp = []
        for cc in cs:
            for k in range(3):
                parts.append(cc)
                e = np.zeros(3)
                e[k] = 1.0
                sp.append(e)
        out, s1 = self._track(grp, parts, sp)
        T = np.zeros((6, 6))
        Rs = []
        for o in range(13):
            Rs.append(np.array([s1[3 * o + k] for k in range(3)]).T)
        for j in range(6):
            T[:, j] = (out[3 * (1 + 2 * j)] - out[3 * (2 + 2 * j)]) / (2 * self.H[j])
        R = Rs[0]
        D = np.zeros((3, 6))
        for j in range(6):
            A = Rs[1 + 2 * j] @ Rs[2 + 2 * j].T
            D[:, j] = -_axial(A - np.eye(3)) / (2 * self.H[j])
        return T, R, D, out[0]

    def run(self):
        c = self.closed_orbit()
        self.cod = c
        recs = []
        cc = c.copy()
        for kind, grp in self.groups:
            T, R, D, cout = self.maps(grp, cc)
            recs.append(dict(kind=kind, grp=grp, T=T, R=R, D=D, cin=cc.copy(), cout=cout.copy()))
            cc = cout
        self.closure = np.abs(cc - c).max()
        Rtot = np.eye(3)
        Ttot = np.eye(6)
        for r in recs:
            Rtot = r['R'] @ Rtot
            Ttot = r['T'] @ Ttot
        w, v = np.linalg.eig(Rtot)
        k = int(np.argmin(abs(w - 1)))
        n0 = np.real(v[:, k])
        n0 /= np.linalg.norm(n0)
        if n0[np.argmax(abs(n0))] < 0:
            n0 = -n0
        tinv = [np.linalg.inv(r['T']) for r in recs]

        def chain(N, store=False):
            n = n0
            out = []
            for r, ti in zip(recs, tinv):
                nin = n
                n = r['R'] @ n
                if store:
                    out.append((N.copy(), nin, n.copy()))
                N = (r['R'] @ N + _cross_mat(n) @ r['D']) @ ti
            return N, out
        Cn, _ = chain(np.zeros((3, 6)))
        Tinv_tot = np.linalg.inv(Ttot)
        A = np.kron(Rtot, Tinv_tot.T)
        N0 = np.linalg.solve(np.eye(18) - A, Cn.ravel()).reshape(3, 6)
        Nend, store = chain(N0, True)
        self.N0 = N0
        self.n0 = n0
        self.nchk = float(np.abs(Nend - N0).max())
        sum_w = sum_wn = sum_wm = sum_wnv = sum_wm2 = sum_wnb2 = 0.0
        for r, (Ni, nin, nout) in zip(recs, store):
            if r['kind'] != 'el':
                continue
            e = r['grp'][0]
            L = getattr(e, 'L', 0.0) or 0.0
            tn = type(e).__name__
            if L <= 0.0 or tn not in ('Bend', 'Quadrupole'):
                continue
            ci, co = r['cin'], r['cout']
            xm, ym = 0.5 * (ci[0] + co[0]), 0.5 * (ci[2] + co[2])
            pxm, pym = 0.5 * (ci[1] + co[1]), 0.5 * (ci[3] + co[3])
            xe, ye, ct, st = _to_el(e, xm, ym)
            bxe, bye = _el_field(e, xe, ye)
            B = np.array([ct * bxe + st * bye, -st * bxe + ct * bye, 0.0])
            pz = math.sqrt(max(0.0, 1.0 - pxm * pxm - pym * pym))
            v = np.array([pxm, pym, pz])
            v /= np.linalg.norm(v)
            Bp = B - (B @ v) * v
            rho_inv = float(np.linalg.norm(Bp))
            if rho_inv == 0.0:
                continue
            b = Bp / rho_inv
            nm = nin + nout
            nm /= np.linalg.norm(nm)
            Nm = 0.5 * (Ni + (r['R'] @ Ni + _cross_mat(nout) @ r['D']) @ np.linalg.inv(r['T']))
            m = Nm[:, 5] - nm * (nm @ Nm[:, 5])
            w_ = rho_inv ** 3 * L
            sum_w += w_
            sum_wn += w_ * float(b @ nm)
            sum_wm += w_ * float(b @ m)
            sum_wnv += w_ * float(nm @ v) ** 2
            sum_wm2 += w_ * float(m @ m)
            sum_wnb2 += w_ * float(b @ nm) ** 2
        # 全ての座標の放射の拡散 (K64 の 1 刻みの拡散行列 b_s = 刻みの出口の B。Emittance[] の POL の量・SLIM のスピン軸の場は使わない) と刻みの出口の N から、
        # 脱偏極率 ½ Σ tr(N⊥ b N⊥ᵀ) (Σ は刻み)
        ph = emit6d._Phys(self.p, self.circ)
        rate_full = 0.0
        for r, (Ni, nin, nout) in zip(recs, store):
            if r['kind'] != 'el':
                continue
            e = r['grp'][0]
            if (getattr(e, 'L', 0.0) or 0.0) <= 0.0 or type(e).__name__ not in ('Bend', 'Quadrupole'):
                continue
            st = emit6d.tturne([e], r['cin'], ph, trf0=self.trf0, rfsw=True, radcod=False)
            Nout = (r['R'] @ Ni + _cross_mat(nout) @ r['D']) @ np.linalg.inv(r['T'])
            Np = (np.eye(3) - np.outer(nout, nout)) @ Nout
            rate_full += 0.5 * float(np.trace(Np @ st.B @ Np.T))
        self.rate_full = rate_full
        gam = math.sqrt(1.0 + (self.p / ME) ** 2)
        k0 = 5.0 * math.sqrt(3.0) / 8.0 * RE * LC * C * gam ** 5
        c1 = sum_w - 2.0 / 9.0 * sum_wnv
        res = dict(n0=n0, sum_w=sum_w, tau_st=1.0 / (k0 * c1 / self.circ), tau_b2=1.0 / (k0 * sum_wnb2 / self.circ),
                   two_ninth=2.0 / 9.0 * sum_wnv / sum_w, m2_over_w=sum_wm2 / sum_w,
                   kin_over_num=(sum_wm / sum_wn) if sum_wn else float('nan'),
                   P_nokin=abs(PSTC * sum_wn / c1), P_dk=abs(PSTC * (sum_wn - sum_wm) / (c1 + 11.0 / 18.0 * sum_wm2)),
                   P_dk_nodep=abs(PSTC * (sum_wn - sum_wm) / c1), P_num_nokin_dep=abs(PSTC * sum_wn / (c1 + 11.0 / 18.0 * sum_wm2)),
                   k0=k0, c1=c1, wn=sum_wn, wm=sum_wm, wm2=sum_wm2, circ=self.circ)
        res['dep_over_c1'] = 11.0 / 18.0 * sum_wm2 / c1
        # rate_full は 1 周あたり。ST の緩和率の 1 周あたりは k0 c1 / c (k0 は秒あたり・c1 = Σ ρ⁻³ ds [m⁻²])
        res['rate_full'] = rate_full
        res['srel_turn'] = k0 * c1 / C
        res['P_dk_full'] = abs(PSTC * (sum_wn - sum_wm) / (c1 + rate_full * C / k0))
        res['dep_full_over_c1'] = rate_full * C / k0 / c1
        return res
