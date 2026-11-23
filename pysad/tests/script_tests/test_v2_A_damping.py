# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): 放射減衰の分配 (Jx, Jy, Je) — 偏向磁石の端の角と、勾配つき偏向磁石。

K64 の動作と根拠 (docs/history/v2_changes_A.md の JX1・JX2・JX3):
- JX1: K64 の Emittance[] (tradke) と追跡 (tradk) は、端の角 e の蹴り x·tan(e)/ρ を偏向角に入れるが、場の中の行路を
  x·tan(e) だけ短くしないので、I4 の端の項 -η·tan(e)/ρ² が 2 倍になる (矩形磁石のリングで Jx が 9% 大きい)。
- JX2: K64 の追跡 (tbendi、K1≠0 の偏向磁石) は、放射で運動量が変わっても曲率半径 rhoe を入口の値のままにするので、
  減衰の和 Jx+Jy+Je が 4 にならない (K1=0.1 で 4.05)。
- 根拠は 2 つの独立な計算: (a) 解析式 Jx = 1 - I4/I2 (自前の 3x3 行列で積分)、(b) 硬い端の偏向磁石の厳密な幾何
  (円弧と直線) で運動量のずれた閉軌道を求め、場の中の行路 l(δ) から Je = 2 + dln(l)/dδ (Robinson)。(a) と (b) は 1e-6 で一致する。

条件: 3 GeV、24 セルの FODO リング (周長 105.6 m)、偏向磁石 L=1.2 m・15 度、RF 500 MHz 4 MV。
"""
import io
import math

import numpy as np
import pytest

from pysad import ffscmd, track6n

TH = 2 * np.pi / 24
LB = 1.2
RHO = LB / TH


# ------------------------------------------------------------------ (a) 解析式 Jx = 1 - I4/I2
def jx_analytic(er, k1b):
    e = er * TH
    kb = k1b / LB

    def drift(l):
        return np.array([[1, l, 0], [0, 1, 0], [0, 0, 1.0]])

    def quad(k, l):
        if k > 0:
            w = np.sqrt(k)
            return np.array([[np.cos(w * l), np.sin(w * l) / w, 0], [-w * np.sin(w * l), np.cos(w * l), 0], [0, 0, 1.0]])
        w = np.sqrt(-k)
        return np.array([[np.cosh(w * l), np.sinh(w * l) / w, 0], [w * np.sinh(w * l), np.cosh(w * l), 0], [0, 0, 1.0]])

    def bend(l):
        K = 1 / RHO**2 + kb
        m = quad(K, l)
        if K > 0:
            w = np.sqrt(K)
            m[0, 2] = (1 - np.cos(w * l)) / (RHO * K)
            m[1, 2] = np.sin(w * l) / (RHO * w)
        else:
            w = np.sqrt(-K)
            m[0, 2] = (1 - np.cosh(w * l)) / (RHO * K)
            m[1, 2] = np.sinh(w * l) / (RHO * w)
        return m

    def edge():
        return np.array([[1, 0, 0], [np.tan(e) / RHO, 1, 0], [0, 0, 1.0]])

    seq = [quad(0.7 / 0.3, 0.3), drift(0.5), 'B', drift(0.5), quad(-0.6 / 0.3, 0.3), drift(1.6)]
    M = np.eye(3)
    for m in seq:
        M = (edge() @ bend(LB) @ edge() @ M) if isinstance(m, str) else (m @ M)
    eta = np.linalg.solve(np.eye(2) - M[:2, :2], M[:2, 2])
    X = np.array([eta[0], eta[1], 1.0])
    I2 = I4 = 0.0
    n = 1000
    for m in seq:
        if isinstance(m, str):
            I4 -= X[0] * np.tan(e) / RHO**2
            X = edge() @ X
            ds = LB / n
            b = bend(ds)
            for _ in range(n):
                X1 = b @ X
                I4 += 0.5 * (X[0] + X1[0]) * (1 / RHO**3 + 2 * kb / RHO) * ds
                I2 += ds / RHO**2
                X = X1
            I4 -= X[0] * np.tan(e) / RHO**2
            X = edge() @ X
        else:
            X = m @ X
    return 1 - I4 / I2


# ------------------------------------------------------------------ (b) 厳密な幾何 (硬い端、K1=0) と Robinson の定理
def jx_geometric(er):
    e1 = er * TH

    def bend(x, phi, d):
        # 参照軌道は (0,0) から +Z へ進み -X 側へ曲がる。外側 = +X。磁場の中は半径 RHO(1+d) の円弧
        p = np.array([x, 0.0])
        dr = np.array([math.sin(phi), math.cos(phi)])
        f1 = np.array([math.cos(e1), math.sin(e1)])
        t, _s = np.linalg.solve(np.array([[dr[0], -f1[0]], [dr[1], -f1[1]]]), -p)
        p = p + t * dr
        P2 = np.array([-RHO + RHO * math.cos(TH), RHO * math.sin(TH)])
        f2 = np.array([math.cos(TH - e1), math.sin(TH - e1)])
        n2 = np.array([-f2[1], f2[0]])
        rho = RHO * (1 + d)
        c = p + rho * np.array([-dr[1], dr[0]])

        def arc(a):
            ca, sa = math.cos(a), math.sin(a)
            r = p - c
            return (c + np.array([ca * r[0] - sa * r[1], sa * r[0] + ca * r[1]]),
                    np.array([ca * dr[0] - sa * dr[1], sa * dr[0] + ca * dr[1]]))
        lo, hi = 0.0, 2 * TH
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if np.dot(arc(mid)[0] - P2, n2) >= 0.0:
                hi = mid
            else:
                lo = mid
        a = 0.5 * (lo + hi)
        p, dr = arc(a)
        nout = np.array([math.cos(TH), math.sin(TH)])
        h = np.array([-math.sin(TH), math.cos(TH)])
        p = p - np.dot(p - P2, h) / np.dot(dr, h) * dr
        return np.dot(p - P2, nout), math.atan2(np.dot(dr, nout), np.dot(dr, h)), rho * a

    def quad(x, phi, d, L, K1):
        k = K1 / L / (1 + d)
        xp = math.tan(phi)
        if k > 0:
            w = math.sqrt(k)
            x, xp = math.cos(w * L) * x + math.sin(w * L) / w * xp, -w * math.sin(w * L) * x + math.cos(w * L) * xp
        else:
            w = math.sqrt(-k)
            x, xp = math.cosh(w * L) * x + math.sinh(w * L) / w * xp, w * math.sinh(w * L) * x + math.cosh(w * L) * xp
        return x, math.atan(xp)

    def cell(w, d):
        x, phi = quad(w[0], w[1], d, 0.3, 0.7)
        x += 0.5 * math.tan(phi)
        x, phi, lf = bend(x, phi, d)
        x += 0.5 * math.tan(phi)
        x, phi = quad(x, phi, d, 0.3, -0.6)
        x += 1.6 * math.tan(phi)
        return np.array([x, phi]), lf

    def lfield(d):
        w = np.zeros(2)
        for _ in range(20):
            r = cell(w, d)[0] - w
            J = np.zeros((2, 2))
            for i in range(2):
                dw = np.zeros(2)
                dw[i] = 1e-7
                J[:, i] = (cell(w + dw, d)[0] - cell(w - dw, d)[0]) / 2e-7
            w = w - np.linalg.solve(J - np.eye(2), r)
        return cell(w, d)[1]
    h = 1e-5
    je = 2 + (lfield(h) - lfield(-h)) / (2 * h) / LB       # U ∝ (1+δ)² l(δ)
    return 3 - je


# ------------------------------------------------------------------ pysad
def _lat(er, k1b, extra=''):
    return '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=%s E2=%s K1=%s);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
USE RINGQ; CELL; CALCULATE; RFSW; RAD; NOFLUC; RADCOD;
%s
e=Emittance[]; Print[DampingRate/.e]; Print[EnergyLossU0/.e];
''' % (er, er, k1b, extra)


def _emit_j(er, k1b, extra=''):
    buf = io.StringIO()
    ffs, _I = ffscmd.run_sad_string(_lat(er, k1b, extra), out=buf)
    o = [l for l in buf.getvalue().split('\n') if l.strip()]
    d = [float(x) for x in o[-2].strip('{}').split(',')]
    k = float(o[-1]) / 6e9
    return ffs, [-x / k for x in d], k


def _track_j(ffs, k):
    """追跡 (RAD、FLUC なし) の一周写像を閉軌道のまわりで数値微分し、固有値の絶対値から減衰率を求める"""
    els = ffs.lat.elements
    p = ffs.momentum()

    def f(v):
        return np.asarray(track6n.track_particles_fast(els, np.array([v]), 1, p, 0.0, True, True, False))[0][:6]
    h = np.array([1e-6, 1e-6, 1e-6, 1e-6, 1e-5, 1e-6])

    def jac(v):
        J = np.zeros((6, 6))
        for i in range(6):
            d = np.zeros(6)
            d[i] = h[i]
            J[:, i] = (f(v + d) - f(v - d)) / (2 * h[i])
        return J
    v = np.zeros(6)
    for _ in range(6):
        v = v - np.linalg.solve(jac(v) - np.eye(6), f(v) - v)
    ev = np.linalg.eigvals(jac(v))
    # 対 (複素共役) ごとに 1 つ。x・y・z の順は虚部 (チューン) の大きさでなく、Jy≈1、Je≈2 で見分ける
    r = sorted(set(round(-math.log(abs(z)) / k, 9) for z in ev))
    return r


def test_references_agree():
    """解析式と厳密な幾何 (Robinson) が一致する: 根拠の 2 つが独立に同じ値"""
    for er in (0.0, 0.1308996938995747, 0.5):
        assert abs(jx_analytic(er, 0.0) - jx_geometric(er)) < 2e-6
    assert abs(jx_analytic(0.5, 0.0) - 1.0057510) < 1e-6


def test_emit_edge_term_v2(v2_mode):
    """矩形磁石 (E1=E2=0.5): v2 の Emittance[] の Jx は解析値 1.00575 に一致 (既定の分割で 1e-3、EMITDIV=20 で 3e-4)"""
    ref = jx_analytic(0.5, 0.0)
    _f, j, _k = _emit_j(0.5, 0.0)
    assert abs(j[0] - ref) < 1e-3
    assert abs(j[1] - 1.0) < 3e-4
    assert abs(sum(j) - 4.0) < 2e-3
    _f, j, _k = _emit_j(0.5, 0.0, 'EMITDIV=20;')
    assert abs(j[0] - ref) < 3e-4


def test_emit_edge_term_k64(k64_mode):
    """K64 は端の項がほぼ 2 倍: Jx = 1.0973 (扇形との差が解析値の差の 2 倍。分散が少し違うので 1e-2 で比べる)"""
    _f, j, _k = _emit_j(0.5, 0.0)
    assert abs(j[0] - 1.09735) < 2e-4
    assert abs(j[0] - (1.0 + 2 * (jx_analytic(0.5, 0.0) - jx_analytic(0.0, 0.0)) - (1 - jx_analytic(0.0, 0.0)))) < 1e-2


def test_emit_sector_unchanged(v2_mode):
    """端の角がなければ v2 と k64 は同じ式 (扇形磁石は解析値と 2e-3 で一致)"""
    _f, j, _k = _emit_j(0.0, 0.0)
    assert abs(j[0] - jx_analytic(0.0, 0.0)) < 2e-3


def test_track_edge_and_sum_rule_v2(v2_mode):
    """追跡の一周写像の減衰: 矩形磁石で Jx が解析値、勾配つき偏向磁石 (K1=0.1) で和が 4、Jx が解析値"""
    ffs, _j, k = _emit_j(0.5, 0.0)
    r = _track_j(ffs, k)
    assert abs(sum(r) - 4.0) < 1.5e-3
    jx = [x for x in r if abs(x - 1.0) > 1e-4 and x < 1.5][0]
    assert abs(jx - jx_analytic(0.5, 0.0)) < 4e-3
    ffs, _j, k = _emit_j(0.0, 0.1)
    r = _track_j(ffs, k)
    assert abs(sum(r) - 4.0) < 1.5e-3
    assert abs(min(r) - jx_analytic(0.0, 0.1)) < 2e-3


def test_track_k64_breaks_sum_rule(k64_mode):
    """K64 の追跡: 矩形磁石で Jx=1.096 (端の項 2 倍)、K1=0.1 の偏向磁石で和が 4.05"""
    ffs, _j, k = _emit_j(0.5, 0.0)
    r = _track_j(ffs, k)
    assert max(x for x in r if x < 1.5) > 1.09
    ffs, _j, k = _emit_j(0.0, 0.1)
    r = _track_j(ffs, k)
    assert abs(sum(r) - 4.05) < 5e-3


def test_emitdiv_is_used(k64_mode):
    """大域変数 EMITDIV (放射の分割の倍率) が Emittance[] に効く。gs: EMITDIV=20 で DampingRate の x が -.000286831320204 (E1=E2=0.5、
    セル QF D2 B D2 QD 3*DS D1 DS の 24 セル、周長 100.8 m)"""
    src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS DS DS D1 DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
USE RINGQ; CELL; CALCULATE; RFSW; RAD;
e=Emittance[]; Print[DampingRate/.e];
EMITDIV=20;
e=Emittance[]; Print[DampingRate/.e];
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    o = [l for l in buf.getvalue().split('\n') if l.startswith('{')]
    a = [float(x) for x in o[0].strip('{}').split(',')]
    b = [float(x) for x in o[1].strip('{}').split(',')]
    assert abs(a[0] / -.0002864491710759299 - 1) < 1e-9
    assert abs(b[0] / -.00028683132020400457 - 1) < 1e-9
    assert abs(b[1] / -.0002605558976010172 - 1) < 1e-9
    assert abs(b[2] / -.0004945592197001899 - 1) < 1e-9


def test_emitdiv_in_main_section(k64_mode):
    """MAIN 部 (FFS; の前) の `EMITDIV = 20;` も効く (gs: -.00028683132020400457)"""
    src = '''MOMENTUM = 3GEV;
EMITDIV = 20;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS DS DS D1 DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
USE RINGQ; CELL; CALCULATE; RFSW; RAD;
e=Emittance[]; Print[DampingRate/.e];
'''
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    o = [l for l in buf.getvalue().split('\n') if l.startswith('{')]
    a = [float(x) for x in o[0].strip('{}').split(',')]
    assert abs(a[0] / -.00028683132020400457 - 1) < 1e-9
