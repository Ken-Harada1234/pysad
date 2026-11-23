# -*- coding: utf-8 -*-
"""K1 = 0 の BEND の本体 (円弧) と端面角 (05_bend.md §4〜§6、式 (5.21)〜(5.40))。

(a) 端面までの直進の写像 (面の座標 x2 = x1/(cosψ - (px/pz) sinψ) など) がシンプレクティック (sympy、6 変数 (x, px, y, py, z, δ))。
(b) 独立な基準: 実験室系 (X, Z) の幾何で、1 粒子を (入口の基準面) → (入口の面) → (一様磁場の円) → (出口の面) → (出口の基準面) と通す。
    端面の Maxwell フリンジ (Δx = ff、Δq_y = -q f、Δz = -q ff) は面の座標系で、本文の式 (5.31) をそのまま入れる (面を横切る点での直進・円・直進は幾何だけで解く)。
    pysad の追跡 (track6._tbendcore: tbend.f tbendcore の直訳。k64・v2 とも) と 5 成分を比べる。
    場合: 扇形 (E1 = E2 = 0)、矩形 (E1 = E2 = 0.5)、非対称 (E1 = 0.3、E2 = 0.7、AE1・AE2)、K0 ≠ 0 (ρb ≠ ρ0)、δ = 0・+0.02・-0.03。
(c) 同じ場合の ρ0 ≠ ρb の幾何 (drhob)。
(d) 設計軌道が (0, 0, 0, 0, 0, 0) のまま出てくること。
"""
import math
import numpy as np
import sympy as sp
from scipy.optimize import brentq
import _common as C

C.start(__file__)

print('== (a) 端面までの直進の写像のシンプレクティック性 (6 変数 (x, px, y, py, z, δ)) ==')
X, PX, Y, PY, Z, DL, psi = sp.symbols('x p_x y p_y z delta psi', real=True)
pp = 1 + DL
pz = sp.sqrt(pp**2 - PX**2 - PY**2)
x2 = X / (sp.cos(psi) - PX / pz * sp.sin(psi))
mp_ = [x2, PX * sp.cos(psi) + pz * sp.sin(psi), Y + PY / pz * x2 * sp.sin(psi), PY, Z - pp / pz * x2 * sp.sin(psi), DL]
vars_ = [X, PX, Y, PY, Z, DL]
J = sp.Matrix(6, 6, lambda i, j: sp.diff(mp_[i], vars_[j]))
S = sp.zeros(6, 6)
for i in range(3):
    S[2 * i, 2 * i + 1] = 1
    S[2 * i + 1, 2 * i] = -1
sub = {X: sp.Rational(3, 1000), PX: sp.Rational(2, 1000), Y: sp.Rational(2, 1000), PY: sp.Rational(-1, 1000), Z: 0, DL: sp.Rational(2, 100), psi: sp.Rational(3, 10)}
Jn = J.subs(sub).evalf(40)
res = (Jn.T * S * Jn - S).applyfunc(lambda v: abs(v))
print('  x2 = x1/(cos ψ - (px/pz) sin ψ)、px2 = px cos ψ + pz sin ψ、y2 = y1 + (py/pz) x2 sin ψ、z2 = z1 - (p/pz) x2 sin ψ (p = 1 + δ、pz = sqrt(p^2 - px^2 - py^2))')
print('  |J^T S J - S| の最大 = %.1e (40 桁の評価。ψ = 0.3 rad、(x, px, y, py, δ) = (3 mm, 2 mrad, 2 mm, -1 mrad, 0.02))' % float(max(res)))
# 出口の面から出口の基準面への直進: px2 = px cosψ + pz sinψ、x2 = x1 (cosψ + (px2/pz2) sinψ)、... 同様に (入口の逆写像)
print('  出口の面から基準面への直進 (ヘルプの "drift from the exit face") は、この写像の ψ → -ψ の逆写像と同じ形: 式 (5.23)')

print('\n== (b) 実験室系の幾何 (X, Z) との比較 ==')
p0c = 3e9 / 0.51099895e6      # 電子 3 GeV の β0γ0 (mc 単位)
h0 = math.sqrt(1 + p0c**2)


def dv_of(dp):
    pr = 1 + dp
    h1 = math.sqrt(1 + (p0c * pr)**2)
    return -dp * (1 + pr) / h1 / (h1 + pr * h0)


def maxwell(xf, qx, y, qy, rho_b, p, enter):
    """端面の Maxwell フリンジ (本文 式 (5.31) の閉じた形)。入力: 面の座標 xf、面の座標系の q = px/p、y、qy。戻り値: (xf, qy, dz)"""
    rhoe = rho_b * p
    pv1sqi = 1.0 / max(1e-4, 1 - qx**2)
    fa = y / rhoe * math.sqrt(pv1sqi)
    f = (1 - (y / rho_b)**2 / 6) * fa
    ff = .25 * (f + fa) * y * pv1sqi
    if enter:
        return xf + ff, qy - qx * f, -qx * ff
    return xf - ff, qy + qx * f, +qx * ff


def cartesian(c, L, phi, K0, psi1, psi2, dp):
    """c = (x, px, y, py, z) 正準 (px, py は p0 で規格化)。戻り値: 出口の基準面の (x, px, y, py, z)"""
    x1, px1, y1, py1, z1 = c
    p = 1 + dp
    rho0 = L / phi
    rho_b = L / (phi + K0)
    hb = 1 / rho_b
    Px, Py = px1, py1
    Pz = math.sqrt(p * p - Px * Px - Py * Py)
    # --- 入口の基準面 Z = 0 から入口の面 (原点を通り、方向 (cosψ1, sinψ1)) まで直進
    ef1 = np.array([math.cos(psi1), math.sin(psi1)])
    nf1 = np.array([-math.sin(psi1), math.cos(psi1)])
    t1 = x1 * math.sin(psi1) / (Pz * math.cos(psi1) - Px * math.sin(psi1))
    r = np.array([x1 + Px * t1, Pz * t1])
    yy = y1 + Py * t1
    tau = t1
    # --- 面の座標系で Maxwell フリンジ
    xf = r @ ef1
    pxf = np.array([Px, Pz]) @ ef1
    xf2, qy2, dz = maxwell(xf, pxf / p, yy, Py / p, rho_b, p, True)
    zacc = z1 + dz
    Py = qy2 * p
    r = r + (xf2 - xf) * ef1
    pzf = math.sqrt(p * p - pxf**2 - Py**2)
    Pvec = pxf * ef1 + pzf * nf1
    Px, Pz = Pvec
    # --- 一様磁場の円 (向きは Θ だけ回る: Px' = Px cosΘ - Pz sinΘ、Pz' = Pz cosΘ + Px sinΘ。X' = X + ρb[Px sinΘ + Pz (cosΘ - 1)]、Z' = Z + ρb[Pz sinΘ + Px (1 - cosΘ)])
    ex2 = np.array([math.cos(phi), math.sin(phi)])
    n2 = np.array([-math.sin(phi), math.cos(phi)])
    E = np.array([-rho0 * (1 - math.cos(phi)), rho0 * math.sin(phi)])
    N2 = math.sin(psi2) * ex2 + math.cos(psi2) * n2         # 出口の面の法線
    ef2 = math.cos(psi2) * ex2 - math.sin(psi2) * n2         # 出口の面の方向

    def arc(Th):
        Xn = r[0] + rho_b * (Px * math.sin(Th) + Pz * (math.cos(Th) - 1))
        Zn = r[1] + rho_b * (Pz * math.sin(Th) + Px * (1 - math.cos(Th)))
        return np.array([Xn, Zn])

    g = lambda Th: (arc(Th) - E) @ N2
    # 面に当たる最初の Θ > 0 を探す (ほぼ φ)
    Ths = np.linspace(1e-9, 2.5 * phi + 0.2, 4000)
    gv = [g(T) for T in Ths]
    Th = None
    for i in range(len(Ths) - 1):
        if gv[i] * gv[i + 1] < 0:
            Th = brentq(g, Ths[i], Ths[i + 1], xtol=1e-15, rtol=1e-15)
            break
    if Th is None:
        raise RuntimeError('出口の面に当たらない')
    r = arc(Th)
    Pn = np.array([Px * math.cos(Th) - Pz * math.sin(Th), Pz * math.cos(Th) + Px * math.sin(Th)])
    yy = yy + Py * (Th / hb)
    tau += Th / hb
    # --- 出口の面で Maxwell フリンジ (面の座標系: x = (r - E)·ef2、px = P·ef2)
    xf = (r - E) @ ef2
    pxf = Pn @ ef2
    xf2, qy2, dz = maxwell(xf, pxf / p, yy, Py / p, rho_b, p, False)
    zacc += dz
    Py = qy2 * p
    r = r + (xf2 - xf) * ef2
    pzf = math.sqrt(p * p - pxf**2 - Py**2)
    Pn = pxf * ef2 + pzf * N2
    # --- 出口の面から出口の基準面 (n2·(r - E) = 0) まで直進 (向きは符号つき)
    t3 = -((r - E) @ n2) / (Pn @ n2)
    r = r + Pn * t3
    yy = yy + Py * t3
    tau += t3
    xo = (r - E) @ ex2
    pxo = Pn @ ex2
    zo = zacc - p * tau + (1 - dv_of(dp)) * L
    return np.array([xo, pxo, yy, Py, zo])


cases = [('扇形 (E1 = E2 = 0)', 'E1=0 E2=0', 0.0, 0.0, 0.0, 0.0),
         ('矩形 (E1 = E2 = 0.5)', 'E1=0.5 E2=0.5', 0.5, 0.5, 0.0, 0.0),
         ('非対称 (E1 = 0.3、E2 = 0.7)', 'E1=0.3 E2=0.7', 0.3, 0.7, 0.0, 0.0),
         ('AE1 = 0.02、AE2 = -0.01 (扇形)', 'E1=0 E2=0 AE1=0.02 AE2=-0.01', 0.0, 0.0, 0.02, -0.01),
         ('K0 = 0.005 (ρb ≠ ρ0)、E1 = E2 = 0.5', 'E1=0.5 E2=0.5 K0=0.005', 0.5, 0.5, 0.0, 0.0)]
L, phi = 1.0, 0.1
print('  (L = %.1f m、ANGLE = %.2f rad、電子 3 GeV、入射 (x, px, y, py) = (3 mm, 2 mrad, 2 mm, 1 mrad)。|差| は x, px, y, py, z の 5 成分の最大)' % (L, phi))
print('  %-40s %-8s %-12s %-12s' % ('場合', 'δ', '|差| k64', '|差| v2'))
worst = 0.0
for name, extra, e1, e2, ae1, ae2 in cases:
    K0 = 0.005 if 'K0' in extra else 0.0
    ffs, el = C.make_bend(L, phi, extra)
    for dp in (0.0, 0.02, -0.03):
        c = np.array([0.003, 0.002, 0.002, 0.001, 0.0, dp])
        out = []
        for mode in ('k64', 'v2'):
            C.set_mode(mode)
            tr = C.track_one(ffs, c)
            ref = cartesian(c[:5], L, phi, K0, e1 * phi + ae1, e2 * phi + ae2, dp)
            out.append(np.abs(tr[:5] - ref).max())
        C.set_mode('v2')
        worst = max(worst, max(out))
        print('  %-40s %+.2f   %-12.1e %-12.1e' % (name, dp, out[0], out[1]))
print('  最大 %.1e。幾何だけで解いた基準と、K64 の tbendcore (閉じた式) が一致する。' % worst)

print('\n== (c) 設計軌道と K0 ==')
for extra, name in (('E1=0.5 E2=0.5 F1=0.2 FRINGE=1', 'F1 = 0.2、FRINGE = 1、E1 = E2 = 0.5 (K0 = 0)'), ('E1=0.5 E2=0.5 K0=0.005', 'K0 = 0.005 (場の角 ANGLE + K0)')):
    ffs, el = C.make_bend(1.0, 0.1, extra)
    print('  %-46s 設計軌道 (0,0,0,0,0,0) の出口: %s' % (name, np.array2string(C.track_one(ffs, np.zeros(6)), precision=4)))
print('  K0 = 0 なら、F1 のフリンジつきでも設計軌道は 0 のまま (L\' の補正と ρb = L\'/ANGLE で保つ)。K0 ≠ 0 は、設計軌道が -K0 の偏向 (px = -K0) を受ける: 補正磁場として働く。')

print('\n== (d) ヘルプの本体の式 (入口の面 → 出口の面) と、実験室系の幾何 ==')


def arc_geo(x1, px1, y1, py1, p, L, phi, K0, psi1, psi2):
    """入口の面の座標 (x1, px1, y1, py1) → 出口の面の座標 (x2, px2, y2) と、円弧の 3 次元の行程パラメータ τ (経路長 = p τ)"""
    rho0 = L / phi
    rho_b = L / (phi + K0)
    hb = 1 / rho_b
    ef1 = np.array([math.cos(psi1), math.sin(psi1)])
    nf1 = np.array([-math.sin(psi1), math.cos(psi1)])
    pz1 = math.sqrt(p * p - px1**2 - py1**2)
    r = x1 * ef1
    Pv = px1 * ef1 + pz1 * nf1
    Px, Pz = Pv
    ex2 = np.array([math.cos(phi), math.sin(phi)])
    n2 = np.array([-math.sin(phi), math.cos(phi)])
    E = np.array([-rho0 * (1 - math.cos(phi)), rho0 * math.sin(phi)])
    N2 = math.sin(psi2) * ex2 + math.cos(psi2) * n2
    ef2 = math.cos(psi2) * ex2 - math.sin(psi2) * n2

    def arc(Th):
        return np.array([r[0] + rho_b * (Px * math.sin(Th) + Pz * (math.cos(Th) - 1)), r[1] + rho_b * (Pz * math.sin(Th) + Px * (1 - math.cos(Th)))])
    g = lambda Th: (arc(Th) - E) @ N2
    Ths = np.linspace(1e-9, 2.5 * phi + 0.2, 4000)
    gv = [g(T) for T in Ths]
    Th = next(brentq(g, Ths[i], Ths[i + 1], xtol=1e-15, rtol=1e-15) for i in range(len(Ths) - 1) if gv[i] * gv[i + 1] < 0)
    rr = arc(Th)
    Pn = np.array([Px * math.cos(Th) - Pz * math.sin(Th), Pz * math.cos(Th) + Px * math.sin(Th)])
    return (rr - E) @ ef2, Pn @ ef2, y1 + py1 * Th / hb, Th / hb


def arc_help(x1, px1, y1, py1, p, L, phi, K0, psi1, psi2):
    """ヘルプの本体の式 (ω = ANGLE - ψ1 - ψ2)。s は面内の経路長 ρb Θ ph (ph = sqrt(p^2 - py^2))"""
    rho0 = L / phi
    rho_b = L / (phi + K0)
    om = phi - psi1 - psi2
    ph = math.sqrt(p * p - py1**2)
    pz1 = math.sqrt(ph**2 - px1**2)
    px2 = -(rho0 / rho_b) * (math.sin(psi2) + math.sin(om + psi1)) + pz1 * math.sin(om) + px1 * math.cos(om) - x1 * math.sin(om) / rho_b
    pz2 = math.sqrt(ph**2 - px2**2)
    x2 = x1 * math.cos(om) + rho_b * (pz2 - pz1 * math.cos(om) + px1 * math.sin(om)) + rho0 * (math.cos(om + psi1) - math.cos(psi2))
    Th = math.asin(px1 / ph) - math.asin(px2 / ph) + om
    return x2, px2, y1 + py1 * rho_b * Th, rho_b * Th


print('  入口の面の座標 (x1, px1, y1, py1) = (3 mm, 2 mrad, 2 mm, 1 mrad)、p = 1 + δ。差は (x2, px2, y2, τ) の最大。τ = ρb Θ (面内の円弧の角 Θ。経路長 = p τ)')
worst = 0.0
for name, extra, e1, e2, ae1, ae2 in cases:
    K0 = 0.005 if 'K0' in extra else 0.0
    for dp in (0.0, 0.02, -0.03):
        a = arc_geo(0.003, 0.002, 0.002, 0.001, 1 + dp, L, phi, K0, e1 * phi + ae1, e2 * phi + ae2)
        b = arc_help(0.003, 0.002, 0.002, 0.001, 1 + dp, L, phi, K0, e1 * phi + ae1, e2 * phi + ae2)
        dd = np.abs(np.array(a) - np.array(b)).max()
        worst = max(worst, dd)
print('  5 つの場合 × δ = 0・0.02・-0.03 の最大の差 %.1e。ヘルプの px2・x2・Θ の式は、実験室系の幾何 (円の中心の座標) と一致する。' % worst)

print('\n== (e) DISFRIN = 1 でも、K1 = 0 の追跡は端面の Maxwell フリンジ (本文 式 (5.31)) を入れている ==')
ffs, el = C.make_bend(L, phi, 'E1=0.5 E2=0.5 DISFRIN=1')
c = np.array([0.003, 0.002, 0.002, 0.001, 0.0, 0.02])
tr = C.track_one(ffs, c)
ref = cartesian(c[:5], L, phi, 0.0, 0.5 * phi, 0.5 * phi, 0.02)
print('  DISFRIN = 1 の追跡と、Maxwell フリンジを入れた幾何の基準の差 %.1e (入れない基準との差は y = 2 mm・δ = 0.02 で下の値)' % np.abs(tr[:5] - ref).max())
import builtins
_mx = maxwell


def _no_maxwell(xf, qx, y, qy, rho_b, p, enter):
    return xf, qy, 0.0


maxwell = _no_maxwell
ref0 = cartesian(c[:5], L, phi, 0.0, 0.5 * phi, 0.5 * phi, 0.02)
maxwell = _mx
print('  Maxwell フリンジなしの基準との差 %.1e (フリンジの効果の大きさ。上の一致が偶然でないことの確認)' % np.abs(tr[:5] - ref0).max())
