# -*- coding: utf-8 -*-
"""FODO_ring.sad (24 セル、3 GeV) の線形光学を、SAD を使わずに 2x2/3x3 行列で計算する検算用モジュール。
.md の「式と数値の照合」表の「式」の列はここで計算した (SAD の E1, E2 は ANGLE に対する比)。
使い方: py -3.12 example/tools/ana_fodo_ring.py  (主な値を表示)"""
import numpy as np
from math import cos, sin, cosh, sinh, sqrt, tan, pi, atan2, acos

THETA = 0.2617993877991494           # ANGLE (rad)
LB = 1.2                              # BEND の L (m)
E1 = E2 = 0.1308996938995747          # E1 E2 (ANGLE に対する比)
CELL = [('SF', 0.2, 0.0), ('D', 0.1, 0), ('QF', 0.3, 0.7), ('D', 0.5, 0), ('B', LB, THETA),
        ('D', 0.5, 0), ('QD', 0.3, -0.6), ('D', 0.1, 0), ('SD', 0.2, 0.0), ('D', 1.0, 0)]
NCELL = 24

def drift(L):
    return np.array([[1, L, 0], [0, 1, 0], [0, 0, 1.0]]), np.array([[1, L], [0, 1.0]])

def quad(L, K1, delta=0.0):
    k = K1 / L / (1 + delta)
    def foc(k):
        if abs(k) < 1e-14: return np.array([[1, L], [0, 1.0]])
        if k > 0:
            s = sqrt(k); return np.array([[cos(s*L), sin(s*L)/s], [-s*sin(s*L), cos(s*L)]])
        s = sqrt(-k); return np.array([[cosh(s*L), sinh(s*L)/s], [s*sinh(s*L), cosh(s*L)]])
    mx = np.eye(3); mx[:2, :2] = foc(k)
    return mx, foc(-k)

def bend(L, theta, e1, e2, delta=0.0):
    """セクター + 辺の収束 (e1, e2 は面の角度 rad)。出力: 水平 3x3 (x, x', delta)、垂直 2x2"""
    th = theta / (1 + delta); rho = L / th
    c, s = cos(th), sin(th)
    sec = np.array([[c, rho*s, rho*(1-c)], [-s/rho, c, s], [0, 0, 1.0]])
    ed = lambda e, sg: np.array([[1, 0], [sg*tan(e)/rho, 1.0]])
    ex1 = np.eye(3); ex1[:2, :2] = ed(e1, +1); ex2 = np.eye(3); ex2[:2, :2] = ed(e2, +1)
    mx = ex2 @ sec @ ex1
    my = ed(e2, -1) @ np.array([[1, L], [0, 1.0]]) @ ed(e1, -1)
    return mx, my

def element(name, L, v, delta=0.0, K1scale=(1.0, 1.0)):
    if name in ('D', 'SF', 'SD'): return drift(L)
    if name == 'QF': return quad(L, v*K1scale[0], delta)
    if name == 'QD': return quad(L, v*K1scale[1], delta)
    if name == 'B': return bend(L, v, E1*v, E2*v, delta)
    raise KeyError(name)

def seq_matrix(seq, delta=0.0, K1scale=(1.0, 1.0)):
    mx, my = np.eye(3), np.eye(2)
    for (n, L, v) in seq:
        a, b = element(n, L, v, delta, K1scale)
        mx, my = a @ mx, b @ my
    return mx, my

def twiss_from_matrix(m2):
    c = (m2[0, 0] + m2[1, 1]) / 2
    if abs(c) >= 1: return None
    mu = acos(c)
    if m2[0, 1] < 0: mu = 2*pi - mu
    beta = m2[0, 1] / sin(mu); alpha = (m2[0, 0] - m2[1, 1]) / (2*sin(mu))
    return beta, alpha, mu

def ring_tunes(delta=0.0, K1scale=(1.0, 1.0)):
    mx, my = seq_matrix(CELL, delta, K1scale)
    tx = twiss_from_matrix(mx[:2, :2]); ty = twiss_from_matrix(my)
    nx = NCELL * tx[2] / 2 / pi if tx else None
    ny = NCELL * ty[2] / 2 / pi if ty else None
    return nx, ny, tx, ty, mx, my

def cell_dispersion(mx):
    """1 セルの周期分散 (EX, EPX)"""
    A = mx[:2, :2]; d = mx[:2, 2]
    return np.linalg.solve(np.eye(2) - A, d)

def cell_integrals(nsl=200):
    """周期解 (1 セル) に沿って、BEND 以外の要素ごとに ∫βx ds, ∫βy ds, ∫ηβx ds, ∫ηβy ds を台形則 (nsl 分割) で求める。
    戻り値: {名前: [dict(Ibx, Iby, Iebx, Ieby, L, K)]}。クロマティシティ (FFS_CHROMATICITY 他) の式の値に使う。"""
    mx, my = seq_matrix(CELL)
    tx = twiss_from_matrix(mx[:2, :2]); ty = twiss_from_matrix(my)
    d0 = cell_dispersion(mx)
    bet = lambda m, b0, a0: m[0, 0]**2*b0 - 2*m[0, 0]*m[0, 1]*a0 + m[0, 1]**2*(1 + a0**2)/b0
    out = {}
    mxa = np.eye(3); mya = np.eye(2)
    for (n, L, v) in CELL:
        Ibx = Iby = Iebx = Ieby = 0.0
        if n == 'B':
            a, b = element(n, L, v); mxa = a @ mxa; mya = b @ mya
        else:
            h = L / nsl
            vals = []
            for i in range(nsl + 1):
                eta = (mxa @ np.array([d0[0], d0[1], 1.0]))[0]
                vals.append((bet(mxa[:2, :2], tx[0], tx[1]), bet(mya, ty[0], ty[1]), eta))
                if i < nsl:
                    a, b = element(n, h, v/nsl if n in ('QF', 'QD') else v)
                    mxa = a @ mxa; mya = b @ mya
            w = np.ones(nsl + 1); w[0] = w[-1] = 0.5
            bx, by, es = (np.array([x[k] for x in vals]) for k in range(3))
            Ibx, Iby, Iebx, Ieby = h*np.sum(w*bx), h*np.sum(w*by), h*np.sum(w*es*bx), h*np.sum(w*es*by)
        out.setdefault(n, []).append(dict(Ibx=Ibx, Iby=Iby, Iebx=Iebx, Ieby=Ieby, L=L, K=v))
    return out

if __name__ == '__main__':
    nx, ny, tx, ty, mx, my = ring_tunes()
    print('ring tune NX NY      = %.6f %.6f' % (nx, ny))
    print('cell BX AX (start SF) = %.6f %.6f' % (tx[0], tx[1]))
    print('cell BY AY           = %.6f %.6f' % (ty[0], ty[1]))
    print('cell trace x y       = %.6f %.6f' % (mx[0, 0] + mx[1, 1], my[0, 0] + my[1, 1]))
    print('EX EPX               = %s' % cell_dispersion(mx))
    dd = 1e-4
    a = ring_tunes(dd)[:2]; b = ring_tunes(-dd)[:2]
    print('xi from tunes of this 3x3 model (bend orbit terms not included; SAD gives -7.512 -6.565) = %.6f %.6f' % ((a[0]-b[0])/2/dd, (a[1]-b[1])/2/dd))


def errors_ring_cod():
    """00_lattice/errors_ring.sad の閉軌道 (線形、ずれのある四極は x_out = M (x_in - d) + d)。
    戻り値: 各要素入口の (x, x', y, y') のリスト (index 1 = ORIGIN, 2 = RF, 3.. = セル)"""
    cell = CELL
    ring = [('MARK', 0, 0, (0, 0)), ('MARK', 0, 0, (0, 0))]
    for icell in range(24):
        for (n, L, v) in cell:
            d = (0.0, 0.0)
            if icell == 0 and n == 'QF': d = (1.0e-4, -2.0e-4)
            if icell == 12 and n == 'QD': d = (0.0, 3.0e-4)
            ring.append((n, L, v, d))
    def step(el, vx, vy):
        n, L, v, d = el
        if n == 'MARK': return vx, vy
        a, b = element(n, L, v)
        a = a[:2, :2]
        dx = np.array([d[0], 0.0]); dy = np.array([d[1], 0.0])
        return a @ (vx - dx) + dx, b @ (vy - dy) + dy
    # 1 周の写像 (アフィン): 単位ベクトルで線形部と定数部を求める
    def one_turn(vx, vy):
        for el in ring: vx, vy = step(el, vx, vy)
        return vx, vy
    z = np.zeros(2)
    cx, cy = one_turn(z, z)
    Mx = np.column_stack([one_turn(np.array([1.0, 0]), z)[0] - cx, one_turn(np.array([0, 1.0]), z)[0] - cx])
    My = np.column_stack([one_turn(z, np.array([1.0, 0]))[1] - cy, one_turn(z, np.array([0, 1.0]))[1] - cy])
    x0 = np.linalg.solve(np.eye(2) - Mx, cx); y0 = np.linalg.solve(np.eye(2) - My, cy)
    out = []; vx, vy = x0, y0
    for el in ring:
        out.append((vx[0], vx[1], vy[0], vy[1]))
        vx, vy = step(el, vx, vy)
    return out
