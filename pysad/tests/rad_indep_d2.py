# -*- coding: utf-8 -*-
"""v2.4.0 D2 の独立な基準: (x, px, δ, z) の放射つき 1 周写像 (RK4) + 薄い空洞 + 硬い端の極 (面を実験室の直線座標で厳密に積分)。SAD・pysad のコードを使わない。"""
import math
import numpy as np
RE = 2.8179403262e-15
MC2 = 0.51099895e6
LC = 3.8615926796e-13


class Ring4E:
    """(x,px,δ,z) 放射つき 1 周写像 + 薄い空洞 + 端の角 (1 次) + 量子励起の Lyapunov。
    elems: dict(type in drift|quad|bend|pole|cav, L, h, hf, k1, e1, e2)。cav は rf_at の位置 (elems の中の 'cav' 要素)。
    端: 入口 px += hf tan(e1) x、Δδ = +cr (1+δ)^2 hf^2 x tan(e1) (行路が x tan e 短い)。出口も同じ (x は出口の値)。"""

    def __init__(self, E_eV, elems, V, harm, ds=0.02, cu=1.0):
        self.g = E_eV / MC2
        self.E = E_eV
        self.V = V
        self.cr = (2 / 3) * RE * self.g ** 3
        self.cu = cu * (55 / (24 * math.sqrt(3))) * RE * LC * self.g ** 5
        self.C = sum(e['L'] for e in elems if e['type'] != 'cav')
        self.k = 2 * math.pi * harm / self.C
        self.S = []
        for e in elems:
            L = e['L']
            if e['type'] in ('drift', 'cav'):
                self.S.append(e)
                continue
            if e['type'] == 'bend' and e.get('e1', 0.0) != 0.0:
                self.S.append(dict(type='edge', h=e['hf'], e=e['e1']))
            n = max(1, int(math.ceil(L / ds)))
            for _ in range(n):
                self.S.append(dict(type=e['type'], L=L / n, h=e.get('h', 0.0) if e['type'] == 'bend' else 0.0,
                                   hf=e.get('hf', 0.0), k1=e.get('k1', 0.0)))
            if e['type'] == 'bend' and e.get('e2', 0.0) != 0.0:
                self.S.append(dict(type='edge', h=e['hf'], e=e['e2']))
        self.phis = 0.0

    def _f(self, y, s):
        x, px, d, z = y
        h = s['h']
        F = s['hf'] + s['k1'] * x
        pz = math.sqrt((1 + d) ** 2 - px * px)
        dl = (1 + h * x) * (1 + d) / pz
        dd = -self.cr * (1 + d) ** 2 * F * F * dl
        return np.array([(1 + h * x) * px / pz, h * pz - (1 + h * x) * F + px * dd / (1 + d), dd, dl - 1.0]), dl, F

    def slice(self, y, s):
        t = s['type']
        if t == 'cav':
            y = y.copy()
            y[2] += self.V / self.E * math.sin(self.phis + self.k * y[3])
            return y, 0.0
        if t == 'edge':
            y = y.copy()
            x, px, d, z = y
            y[1] = px + s['h'] * math.tan(s['e']) * x
            y[2] = d + self.cr * (1 + d) ** 2 * s['h'] ** 2 * x * math.tan(s['e'])
            return y, 0.0
        L = s['L']
        if t == 'drift':
            x, px, d, z = y
            pz = math.sqrt((1 + d) ** 2 - px * px)
            return np.array([x + L * px / pz, px, d, z + L * ((1 + d) / pz - 1)]), 0.0
        if t == 'quad':
            s = dict(s, h=0.0, hf=0.0)
        f = lambda v: self._f(v, s)[0]
        k1 = f(y)
        k2 = f(y + L / 2 * k1)
        k3 = f(y + L / 2 * k2)
        k4 = f(y + L * k3)
        yn = y + L / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        _, dl, F = self._f(0.5 * (y + yn), s)
        return yn, abs(F) ** 3 * dl * L

    def turn(self, y, rf=True):
        for s in self.S:
            if s['type'] == 'cav' and not rf:
                continue
            y, _ = self.slice(y, s)
        return y

    def jac_fn(self, fn, y, h=1e-6):
        J = np.zeros((4, 4))
        for i in range(4):
            d = np.zeros(4)
            d[i] = h
            J[:, i] = (fn(y + d) - fn(y - d)) / (2 * h)
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
            y = y - np.linalg.solve(self.jac_fn(self.turn, y) - np.eye(4), f)
        self.orbit = y
        self.M = self.jac_fn(self.turn, y)
        return np.linalg.eigvals(self.M)

    def equilibrium(self):
        y = self.orbit.copy()
        Qt = np.zeros((4, 4))
        for s in self.S:
            J = self.jac_fn(lambda v: self.slice(v, s)[0], y)
            yn, I = self.slice(y, s)
            w = np.array([0.0, y[1] / (1 + y[2]), 1.0, 0.0])
            Qt = J @ Qt @ J.T + self.cu * I * np.outer(w, w)
            y = yn
        A = np.eye(16) - np.kron(self.M, self.M)
        Sg = np.linalg.solve(A, Qt.reshape(16)).reshape(4, 4)
        return Sg

    def mode_emittances(self):
        """固有エミッタンス (x 様・z 様)。固有ベクトルを 4 次元のシンプレクティック形式で規格化"""
        Sg = self.equilibrium()
        w, V = np.linalg.eig(self.M)
        Jm = np.zeros((4, 4))
        Jm[0, 1] = 1
        Jm[1, 0] = -1
        Jm[2, 3] = 1
        Jm[3, 2] = -1
        # (x,px,δ,z): 正準対 (x,px)、(z,δ): z-δ の形式は J[3,2]=+1 (z が座標、δ が運動量) に直す: 順番を (x,px,z,δ) に並べ替える
        P = [0, 1, 3, 2]
        Mp = self.M[np.ix_(P, P)]
        Sp = Sg[np.ix_(P, P)]
        w, V = np.linalg.eig(Mp)
        Jf = np.zeros((4, 4))
        Jf[0, 1] = 1
        Jf[1, 0] = -1
        Jf[2, 3] = 1
        Jf[3, 2] = -1
        out = []
        cols = []
        for i in range(4):
            if w[i].imag > 1e-9:
                v = V[:, i]
                nrm = (v.conj() @ Jf @ v) / 1j
                if nrm.real < 0:
                    v = V[:, [j for j in range(4) if abs(w[j] - w[i].conjugate()) < 1e-9][0]]
                    nrm = (v.conj() @ Jf @ v) / 1j
                v = v / math.sqrt(abs(nrm.real))
                cols.append((abs(np.angle(w[i])), v, w[i]))
        cols.sort(key=lambda c: -c[0])      # 位相の大きい順: x 様 (大) が先、z 様 (小) が後ろ
        Vm = np.zeros((4, 4), dtype=complex)
        for k, (ph, v, wi) in enumerate(cols):
            Vm[:, 2 * k] = v
            Vm[:, 2 * k + 1] = v.conj()
        Vi = np.linalg.inv(Vm)
        Sc = Vi @ Sp @ Vi.conj().T
        return [float(np.real(Sc[0, 0])), float(np.real(Sc[2, 2]))], [c[2] for c in cols]


class _Ring4X(Ring4E):
    """硬い端の極 'xpole' (一様磁場、面 z = X tan e1・L − X tan e2) を含む"""

    def __init__(self, E_eV, elems, V, harm, ds=0.02, nz=40):
        self.nz = nz
        super().__init__(E_eV, [e for e in elems if e['type'] != 'xpole'], V, harm, ds)
        self.S = []
        for e in elems:
            L = e['L']
            if e['type'] in ('drift', 'cav', 'xpole'):
                self.S.append(e)
                continue
            n = max(1, int(math.ceil(L / ds)))
            for _ in range(n):
                self.S.append(dict(type=e['type'], L=L / n, h=e.get('h', 0.0) if e['type'] == 'bend' else 0.0, hf=e.get('hf', 0.0), k1=e.get('k1', 0.0)))
        self.C = sum(e['L'] for e in elems if e['type'] != 'cav')
        self.k = 2 * math.pi * harm / self.C

    def _fz(self, y, F):
        x, px, d, z = y
        pz = math.sqrt((1 + d) ** 2 - px * px)
        dl = (1 + d) / pz
        dd = -self.cr * (1 + d) ** 2 * F * F * dl
        return np.array([px / pz, -F + px * dd / (1 + d), dd, dl - 1.0])

    def xpole(self, y, s):
        L, F, e1, e2 = s['L'], s['hf'], s['e1'], s['e2']
        t1, t2 = math.tan(e1), math.tan(e2)
        x, px, d, z = y
        pz = math.sqrt((1 + d) ** 2 - px * px)
        zf = x * t1 / (1 - (px / pz) * t1)
        y = np.array([x + (px / pz) * zf, px, d, z + zf * ((1 + d) / pz - 1)])
        zend = L
        for _ in range(8):
            h = (zend - zf) / self.nz
            yy = y.copy()
            for _i in range(self.nz):
                k1 = self._fz(yy, F)
                k2 = self._fz(yy + h / 2 * k1, F)
                k3 = self._fz(yy + h / 2 * k2, F)
                k4 = self._fz(yy + h * k3, F)
                yy = yy + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            znew = L - yy[0] * t2
            if abs(znew - zend) < 1e-15:
                break
            zend = znew
        x, px, d, z = yy
        pz = math.sqrt((1 + d) ** 2 - px * px)
        dz = L - zend
        return np.array([x + (px / pz) * dz, px, d, z + dz * ((1 + d) / pz - 1)]), 0.0

    def slice(self, y, s):
        if s['type'] == 'xpole':
            return self.xpole(y, s)
        return super().slice(y, s)


def Ring4X(E_eV, ae, th, nb, a0=0.06, kq=0.7, V=2e6, harm=40):
    """test_v2_D2_edge の挿入光源リング: RF (薄い空洞) + 3 極 (K0 = a0・−2a0・a0、L = 0.2・0.4・0.2、AE1 = AE2 = ae) + nb セル (QF D B D QD D)"""
    Lb, Lq, Ld, Lw, Ldw = 1.0, 0.3, 0.6, 0.2, 0.3
    elems = [dict(type='cav', L=0.0)]
    for a, l in [(a0, Lw), (-2 * a0, 2 * Lw), (a0, Lw)]:
        elems.append(dict(type='drift', L=Ldw))
        elems.append(dict(type='xpole', L=l, hf=a / l, e1=ae, e2=ae))
    elems.append(dict(type='drift', L=Ldw))
    for _ in range(nb):
        elems += [dict(type='quad', L=Lq, k1=kq / Lq), dict(type='drift', L=Ld), dict(type='bend', L=Lb, h=th / Lb, hf=th / Lb),
                  dict(type='drift', L=Ld), dict(type='quad', L=Lq, k1=-kq / Lq), dict(type='drift', L=Ld)]
    return _Ring4X(E_eV, elems, V, harm)
