# -*- coding: utf-8 -*-
"""担当 (放射・エミッタンス) の独立な計算 (pysad の光学・放射積分のコードを使わない)。試験 test_v24_rad_pole.py から使う。
  integrals_pole: 要素表から 3x3 行列 (x, x', δ) で I1..I5。ANGLE=0・K0≠0 の極は「直線の基準軌道の中の一様な磁場」
                  (η'' + k η = h。扇形磁石の幾何学的な収束 h² も、曲線座標の行路の項 (I4 の h³) も無い) として積分する。
  Ring4: (x, px, δ, z) の放射つき 1 周写像 (運動方程式を RK4 で積分、薄い空洞つき)。固有値の絶対値から Jx = 2 E (1 − |λx|) / U0。
         曲線座標: dx/ds = (1+hx) px/pz、dpx/ds = h pz − (1+hx) F、F = hf + k1 x。直線の極は h = 0 で F = hp + k1 x。
         放射 dδ = −(2/3) re γ³ (1+δ)² F² dl、dl = (1+hx)(1+δ)/pz ds、px は放射で (1+δ) に比例して縮む。z は経路長の遅れ。"""
import math

import numpy as np

from . import _a2_indep as A

RE = 2.8179403262e-15
MC2 = 0.51099895e6


def _pole_matrix(l, h, k):
    if abs(k) < 1e-14:
        return np.array([[1, l, h * l * l / 2], [0, 1, h * l], [0, 0, 1.0]])
    if k > 0:
        w = math.sqrt(k)
        c, s = math.cos(w * l), math.sin(w * l)
        return np.array([[c, s / w, h * (1 - c) / k], [-w * s, c, h * s / w], [0, 0, 1.0]])
    w = math.sqrt(-k)
    c, s = math.cosh(w * l), math.sinh(w * l)
    return np.array([[c, s / w, h * (1 - c) / k], [w * s, c, h * s / w], [0, 0, 1.0]])


def _parts(els):
    out = []
    for e in els:
        t = type(e).__name__
        L = getattr(e, 'L', 0.0) or 0.0
        if t == 'Bend' and L > 0:
            ang = getattr(e, 'ANGLE', 0.0)
            k0 = getattr(e, 'K0', 0.0)
            if ang == 0.0 and k0 != 0.0:
                out.append(('pole', L, k0 / L, getattr(e, 'K1', 0.0), 0.0, 0.0))
            else:
                out.append(('bend', L, (ang + k0) / L, getattr(e, 'K1', 0.0), getattr(e, 'E1', 0.0), getattr(e, 'E2', 0.0)))
        elif t == 'Quadrupole' and L > 0:
            out.append(('quad', L, e.K1))
        elif L > 0:
            out.append(('drift', L))
    return out


def integrals_pole(els, n_per_m=40):
    """I1..I5。K1 は pysad の要素の値 (K1/L ではなく _a2_indep と同じ約束) をそのまま使う (極の K1 は 0 の場合だけ使う)"""
    P = _parts(els)
    C = sum((getattr(e, 'L', 0.0) or 0.0) for e in els)

    def mat(p):
        if p[0] == 'drift':
            return A.drift(p[1])
        if p[0] == 'quad':
            return A.quad(p[2], p[1])
        _, L, h, k, e1, e2 = p
        if p[0] == 'pole':
            return _pole_matrix(L, h, k)
        return A.edge(h, e2) @ A.bend(L, h, k) @ A.edge(h, e1)
    M = np.eye(3)
    for p in P:
        M = mat(p) @ M
    Am = M[:2, :2]
    mu = math.acos(0.5 * np.trace(Am)) * (1 if Am[0, 1] > 0 else -1)
    eta = np.linalg.solve(np.eye(2) - Am, M[:2, 2])
    b0 = Am[0, 1] / math.sin(mu)
    a0 = (Am[0, 0] - Am[1, 1]) / (2 * math.sin(mu))
    X = np.array([eta[0], eta[1], 1.0])
    T = np.array([b0, a0, (1 + a0 * a0) / b0])
    I = np.zeros(5)

    def H(X, T):
        return T[2] * X[0] ** 2 + 2 * T[1] * X[0] * X[1] + T[0] * X[1] ** 2
    for p in P:
        if p[0] in ('drift', 'quad'):
            m = mat(p)
            X = m @ X
            T = A.twiss_step(T, m)
            continue
        _, L, h, k, e1, e2 = p
        if p[0] == 'bend':
            m = A.edge(h, e1)
            I[3] += -X[0] * math.tan(e1) * h * h
            X = m @ X
            T = A.twiss_step(T, m)
        n = max(2, int(L * n_per_m))
        ds = L / n
        mb = _pole_matrix(ds, h, k) if p[0] == 'pole' else A.bend(ds, h, k)
        f4 = (2 * k * h) if p[0] == 'pole' else h * (h * h + 2 * k)
        for _ in range(n):
            X1 = mb @ X
            T1 = A.twiss_step(T, mb)
            I[0] += 0.5 * (X[0] + X1[0]) * h * ds
            I[1] += h * h * ds
            I[2] += abs(h) ** 3 * ds
            I[3] += 0.5 * (X[0] + X1[0]) * f4 * ds
            I[4] += 0.5 * (H(X, T) + H(X1, T1)) * abs(h) ** 3 * ds
            X, T = X1, T1
        if p[0] == 'bend':
            m = A.edge(h, e2)
            I[3] += -X[0] * math.tan(e2) * h * h
            X = m @ X
            T = A.twiss_step(T, m)
    return I, C


def integrals_pole_rich(els, n=80):
    a, C = integrals_pole(els, n)
    b, _ = integrals_pole(els, 2 * n)
    return b + (b - a) / 3.0, C


class Ring4:
    def __init__(self, E_eV, elems, V, harm, ds=0.025):
        self.g = E_eV / MC2
        self.E = E_eV
        self.V = V
        self.cr = (2 / 3) * RE * self.g ** 3
        self.C = sum(e['L'] for e in elems)
        self.k = 2 * math.pi * harm / self.C
        self.S = []
        for e in elems:
            L = e['L']
            if e['type'] == 'drift':
                self.S.append(e)
                continue
            n = max(1, int(math.ceil(L / ds)))
            for _ in range(n):
                self.S.append(dict(type=e['type'], L=L / n, h=e.get('h', 0.0) if e['type'] == 'bend' else 0.0,
                                   hf=e.get('hf', 0.0), k1=e.get('k1', 0.0)))
        self.phis = 0.0

    def _f(self, y, s):
        x, px, d, z = y
        h = s['h']
        F = s['hf'] + s['k1'] * x
        pz = math.sqrt((1 + d) ** 2 - px * px)
        dl = (1 + h * x) * (1 + d) / pz
        dd = -self.cr * (1 + d) ** 2 * F * F * dl
        return np.array([(1 + h * x) * px / pz, h * pz - (1 + h * x) * F + px * dd / (1 + d), dd, dl - 1.0])

    def slice(self, y, s):
        L = s['L']
        if s['type'] == 'drift':
            x, px, d, z = y
            pz = math.sqrt((1 + d) ** 2 - px * px)
            return np.array([x + L * px / pz, px, d, z + L * ((1 + d) / pz - 1)])
        if s['type'] == 'quad':
            s = dict(s, h=0.0, hf=0.0)
        k1 = self._f(y, s)
        k2 = self._f(y + L / 2 * k1, s)
        k3 = self._f(y + L / 2 * k2, s)
        k4 = self._f(y + L * k3, s)
        return y + L / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    def turn(self, y, rf=True):
        if rf:
            y = y.copy()
            y[2] += self.V / self.E * math.sin(self.phis + self.k * y[3])
        for s in self.S:
            y = self.slice(y, s)
        return y

    def jac(self, y, h=1e-6):
        J = np.zeros((4, 4))
        for i in range(4):
            d = np.zeros(4)
            d[i] = h
            J[:, i] = (self.turn(y + d) - self.turn(y - d)) / (2 * h)
        return J

    def solve(self):
        y = np.zeros(4)
        for _ in range(20):
            f = self.turn(y, rf=False) - y
            J = np.zeros((2, 2))
            for i in range(2):
                d = np.zeros(4)
                d[i] = 1e-7
                J[:, i] = ((self.turn(y + d, False) - self.turn(y - d, False)) / 2e-7)[:2]
            if np.max(np.abs(f[:2])) < 1e-16:
                break
            y[:2] -= np.linalg.solve(J - np.eye(2), f[:2])
        y1 = self.turn(y, rf=False)
        self.u0 = -(y1[2] - y[2])
        self.phis = math.pi - math.asin(self.u0 * self.E / self.V)
        for _ in range(30):
            f = self.turn(y) - y
            if np.max(np.abs(f)) < 1e-15:
                break
            y = y - np.linalg.solve(self.jac(y) - np.eye(4), f)
        self.orbit = y
        self.M = self.jac(y)
        return np.linalg.eigvals(self.M)


def toy_ring(nb=16, a0=0.06, kq=0.7):
    """16 セルの小リング (1 GeV)。セル: QF D B D QD D (扇形、E1=E2=0)。先頭の直線部に 3 極の挿入光源 (K0 = a0, −2 a0, a0、長さ 0.2・0.4・0.2 m)。
    戻り値: SAD の台本、Ring4 用の要素表"""
    Lb, Lq, Ld, Lw, Ldw = 1.0, 0.3, 0.6, 0.2, 0.3
    th = 2 * math.pi / nb
    sad = ['MOMENTUM = 1GEV;', 'DRIFT D1=(L=%g) DW=(L=%g);' % (Ld, Ldw), 'QUAD QF=(L=%g K1=%g) QD=(L=%g K1=%g);' % (Lq, kq, Lq, -kq),
           'BEND B=(L=%g ANGLE=%.12g E1=0 E2=0);' % (Lb, th)]
    elems = []
    if a0 != 0:
        sad.append('BEND P1=(L=%g ANGLE=0 K0=%.12g) P2=(L=%g ANGLE=0 K0=%.12g);' % (Lw, a0, 2 * Lw, -2 * a0))
        for a, l in [(a0, Lw), (-2 * a0, 2 * Lw), (a0, Lw)]:
            elems.append(dict(type='drift', L=Ldw))
            elems.append(dict(type='pole', L=l, hf=a / l))
        elems.append(dict(type='drift', L=Ldw))
    for _ in range(nb):
        elems += [dict(type='quad', L=Lq, k1=kq / Lq), dict(type='drift', L=Ld), dict(type='bend', L=Lb, h=th / Lb, hf=th / Lb),
                  dict(type='drift', L=Ld), dict(type='quad', L=Lq, k1=-kq / Lq), dict(type='drift', L=Ld)]
    sad += ['LINE CELL=(QF D1 B D1 QD D1);']
    if a0 != 0:
        sad.append('LINE WIG=(DW P1 DW P2 DW P1 DW);')
    sad += ['CAVI RF=(L=0 VOLT=2E6 HARM=40);', 'LINE RING=(RF %s %d*CELL);' % ('WIG' if a0 != 0 else '', nb)]
    return '\n'.join(sad), elems
