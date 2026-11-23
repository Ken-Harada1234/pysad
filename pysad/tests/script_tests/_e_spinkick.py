# -*- coding: utf-8 -*-
"""v2next E 項目 3 の補助: δ の蹴り (放射の光子 1 個) に対するスピンの不変軸の応答 G = ∂n/∂δ を、追跡 (M2) と線形の軌道モデル (M1) で測る。
方法: 軌道 u の粒子 A と u に δ の ±h を足した B± を、同じスピン (n0) から K 周 (放射の損失なし: emit6d.CAVE を 1e-9 倍) 追い、
ζ = s·(m + i l) (m, l は n0 に垂直な基底) の列の、スピンチューンの線の複素振幅 (NAFF) を a_A、a_B± とする。円錐の振幅の差 a_B+ − a_B− = 2h G なので
G = (a_B+ − a_B−)/(2h) (複素 2 成分 = n⊥ の 2 成分)。軌道の線 (νx, νy, νs とその組み合わせ) とスピンチューンの線は分離する。
M1 = 線形の軌道 + 要素ごとの線形のスピン軌道結合 s' = exp(−[D u]×) R s (B3 の MC と同じ)。放射の脱偏極率は Σ_i ½ B_δδ,i |G_i|² (δ の蹴りが支配的)。
Ring() は emit6d.CAVE を書き換える。呼び出し側が使い終えたら emit6d.CAVE を元に戻すこと (restore_cave)。"""
import io
import math

import numpy as np

from pysad import compat, emit6d, track6, ffscmd, naff as _naff

_CAVE0 = emit6d.CAVE


def restore_cave():
    emit6d.CAVE = _CAVE0


def lat(E, volt, dy=3e-3, harm=64):
    return f'''MOMENTUM = {E}GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY={dy}) QDE = (L=0.3 K1=-0.6 DY=-{dy});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={volt}e6 HARM={harm});
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''


def setup(ffs):
    p = ffs.momentum()
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
    els = ffs.lat.elements
    circ = emit6d._emit_circ(els)
    trf0 = e.get('trf0', 0.0)
    sp = e['spin']
    ctx = dict(p=p, e=e, els=els, circ=circ, trf0=trf0, sp=sp)
    return ctx


def oneturn(ctx, c, s=None):
    if s is None:
        out, _ = track6.track_particles(ctx['els'], np.array(c), 1, ctx['p'], ctx['trf0'], True, True, False,
                                        spins=[np.array([0, 1.0, 0])] * len(c), radpol=False, circ=ctx['circ'])
        return np.array(out)
    out, s1 = track6.track_particles(ctx['els'], np.array(c), 1, ctx['p'], ctx['trf0'], True, True, False,
                                    spins=[np.array(x) for x in s], radpol=False, circ=ctx['circ'])
    return np.array(out), np.array(s1)


def jac(ctx, c0, h=None):
    h = h if h is not None else np.array([1e-6, 1e-6, 1e-6, 1e-6, 1e-5, 1e-6])
    C = []
    for j in range(6):
        d = np.zeros(6)
        d[j] = h[j]
        C.append(c0 + d)
        C.append(c0 - d)
    out = oneturn(ctx, C)
    J = np.zeros((6, 6))
    for j in range(6):
        J[:, j] = (out[2 * j] - out[2 * j + 1]) / (2 * h[j])
    return J


def fixed_point(ctx, c0):
    c = np.array(c0, float)
    for it in range(8):
        f = oneturn(ctx, [c])[0] - c
        J = jac(ctx, c)
        dc = np.linalg.solve(J - np.eye(6), -f)
        c = c + dc
        if np.abs(dc).max() < 1e-13:
            break
    return c, J


def modes(J, tunes):
    """J の固有モード (x, y, z の順): 戻り V (6x6 複素、列 = [e_x, e_x*, e_y, e_y*, e_z, e_z*])"""
    w, v = np.linalg.eig(J)
    ph = np.angle(w) / (2 * np.pi)
    cols = []
    used = set()
    for t in tunes:
        # 固有値の位相が ±t に最も近いもの (正の側)
        idx = [i for i in range(6) if i not in used]
        i = min(idx, key=lambda i: min(abs(ph[i] - abs(t)), abs(ph[i] + abs(t))))
        # 共役のペアを探す
        j = min([k for k in idx if k != i], key=lambda k: abs(w[k] - np.conj(w[i])))
        used.update((i, j))
        cols += [v[:, i], v[:, j]]
    return np.array(cols).T




class Ring:
    def __init__(self, E, V, i0=0, lattice=None, post=None):
        ffs, _ = ffscmd.run_sad_string((lattice or lat)(E, V), out=io.StringIO())
        if post is not None:
            post(ffs)
        self.ffs = ffs
        self.ctx = ctx = setup(ffs)
        self.e = ctx['e']
        self.els_all = list(ctx['els'])
        self.i0 = i0
        self.sp = ctx['sp']
        self.cod_em = np.array(self.e['orbit_entrance'])
        emit6d.CAVE = _CAVE0 * 1e-9
        ctx['els'] = self.els_all[i0:] + self.els_all[:i0]
        # 始点を i0 に回した輪での閉軌道
        c = self.cod_em.copy()
        if i0 > 0:
            out, _ = track6.track_particles(self.els_all[:i0], np.array([c]), 1, ctx['p'], ctx['trf0'], True, True, False,
                                            spins=[np.array([0, 1.0, 0])], radpol=False, circ=ctx['circ'])
            c = out[0]
        self.cod, self.J = fixed_point(ctx, c)
        M = np.zeros((3, 3))
        for j in range(3):
            sj = np.zeros(3); sj[j] = 1
            _o, s1 = oneturn(ctx, [self.cod], [sj])
            M[:, j] = s1[0]
        w, v = np.linalg.eig(M)
        k = int(np.argmin(abs(w - 1)))
        n0 = np.real(v[:, k]); n0 /= np.linalg.norm(n0)
        if n0[1] < 0:
            n0 = -n0
        self.n0 = n0
        a = np.array([1.0, 0, 0]) if abs(n0[0]) < 0.9 else np.array([0, 1.0, 0])
        m = a - (a @ n0) * n0; m /= np.linalg.norm(m)
        self.mvec, self.lvec = m, np.cross(n0, m)
        # 共分散 Σ を始点 i0 へ運ぶ (線形の部分写像)
        Sig0 = np.array(self.e['beam_matrix'])
        if i0 > 0:
            Jp = self._partial_jac(self.cod_em, i0)
            self.Sig = Jp @ Sig0 @ Jp.T
        else:
            self.Sig = Sig0
        tunes = [abs(t) for t in self.e['tunes']]
        self.Vm = modes(self.J, tunes)
        Vi = np.linalg.inv(self.Vm)
        Ca = Vi @ self.Sig @ Vi.conj().T
        self.m = np.array([Ca[2 * i, 2 * i].real for i in range(3)])

    def _partial_jac(self, c0, i0):
        ctx = self.ctx
        h = np.array([1e-6, 1e-6, 1e-6, 1e-6, 1e-5, 1e-6])
        C = []
        for j in range(6):
            d = np.zeros(6); d[j] = h[j]
            C += [c0 + d, c0 - d]
        out, _ = track6.track_particles(self.els_all[:i0], np.array(C), 1, ctx['p'], ctx['trf0'], True, True, False,
                                        spins=[np.array([0, 1.0, 0])] * 12, radpol=False, circ=ctx['circ'])
        Jp = np.zeros((6, 6))
        for j in range(6):
            Jp[:, j] = (out[2 * j] - out[2 * j + 1]) / (2 * h[j])
        return Jp

    def sample(self, mode, f, rng):
        a = np.zeros(6, complex)
        for i in range(3):
            on = (mode == 'both') or (mode == 'beta' and i < 2) or (mode == 'synch' and i == 2)
            if on:
                amp = math.sqrt(self.m[i] * rng.exponential(1.0)) * f
                ph = rng.uniform(0, 2 * np.pi)
                a[2 * i] = amp * np.exp(1j * ph)
                a[2 * i + 1] = np.conj(a[2 * i])
        return self.cod + (self.Vm @ a).real

    def zeta(self, S):
        return S @ self.mvec + 1j * (S @ self.lvec)


def line_amp(z, nu_spin, K):
    """z の列のスピンチューンの線 (|ν| が nu_spin に最も近い周波数) の複素振幅 (k=0 に外挿) と周波数"""
    fr, am = _naff.naff(z, nfreq=8)
    i = int(np.argmin(np.abs(np.abs(fr) - nu_spin) + 0.0))
    return fr[i], am[i]


def track_pair(ring, u, K, hs, nu_spin):
    """u と u+h e_δ の粒子 (h は hs のリスト) を K 周追い、ζ の線の複素振幅 [a_A, a_B..] を返す"""
    ctx = ring.ctx
    c = np.array([u] + [u + np.array([0, 0, 0, 0, 0, h]) for h in hs])
    n = len(c)
    S = np.tile(ring.n0, (n, 1))
    zs = np.zeros((K, n), complex)
    for k in range(K):
        c, S = oneturn(ctx, c, S)
        zs[k] = ring.zeta(S)
    res = []
    for j in range(n):
        fr, am = line_amp(zs[:, j], nu_spin, K)
        res.append((fr, am))
    return res




def build_secs(ring):
    """始点 i0 の輪の要素ごとの (T, R, D) (放射の損失なし、閉軌道 ring.cod の上)"""
    ctx = ring.ctx
    ph = emit6d._Phys(ctx['p'], ctx['circ'])
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    secs = []
    c = ring.cod.copy()
    try:
        for el in ctx['els']:
            st = emit6d.tturne([el], c, ph, trf0=ctx['trf0'], rfsw=True, radcod=False)
            secs.append((st.T.copy(), st.srot[:, :3].copy(), st.srot[:, 3:9].copy()))
            c = st.cod.copy()
    finally:
        emit6d.POL_MODE = old
    return secs


def _rot(w, s):
    """回転ベクトル w (n,3) で s (n,3) を回す (Rodrigues)"""
    th = np.linalg.norm(w, axis=1)
    th = np.where(th < 1e-300, 1e-300, th)
    k = w / th[:, None]
    ct, st = np.cos(th)[:, None], np.sin(th)[:, None]
    return s * ct + np.cross(k, s) * st + k * np.einsum('ij,ij->i', k, s)[:, None] * (1 - ct)


def track_pair_linear(ring, secs, u, K, hs, nu_spin):
    c0 = ring.cod
    uu = np.array([u - c0] + [u - c0 + np.array([0, 0, 0, 0, 0, h]) for h in hs])
    n = len(uu)
    S = np.tile(ring.n0, (n, 1))
    zs = np.zeros((K, n), complex)
    for k in range(K):
        for T, R, D in secs:
            S = S @ R.T
            w = uu @ D.T            # 回転ベクトル (K64: ds = s × g。回転は −g)
            S = _rot(-w, S)
            uu = uu @ T.T
        S /= np.linalg.norm(S, axis=1)[:, None]
        zs[k] = ring.zeta(S)
    res = []
    for j in range(n):
        res.append(line_amp(zs[:, j], nu_spin, K))
    return res
