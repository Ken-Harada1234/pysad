# -*- coding: utf-8 -*-
"""座標の変換 COORD (11_misalignment_and_coord.md §3〜§5): 剛体の運動 + 直線の面への戻し、で作った独立な写像と比べる。

(a) ヘルプの回転行列 (DIR = 0 と DIR != 0) と、K64 の tcoord.f の r_ij (pysad の coord._rot と同じ) を sympy で比べる。
    DIR = 0 は一致。DIR != 0 はヘルプの式と一致しない (K64 とヘルプの食い違い。疑い)。
(b) 独立な写像 (粒子の直線を 3 次元で剛体変換し、面 z' = 0 と交わらせる) と、pysad の追跡 (tcoord)・CALC の写像 (tcoorde) の比較。
(c) CALC の 6x6 行列 (tcoorde の解析式) と、独立な写像の数値微分 (中心差分)。シンプレクティック性。
(d) 小さな角度の 1 次。
(e) DIR = 0 の COORD のあとに DIR != 0 の COORD を置いても元に戻らない (ヘルプの "NOT the inverse")。真の逆は R^T。
"""
import io
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import coord


def mk(el, momentum='3GEV'):
    src = 'MOMENTUM = %s;\n%s\nMARK ST = ();\nLINE ONE = (ST CO);\nFFS;\nUSE ONE;\n' % (momentum, el)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


print('pysad:', C.pysad_src())

print('== (a) ヘルプの回転行列と tcoord.f の r_ij (sympy) ==')
c1, s1, c2, s2, c3, s3 = sp.symbols('c1 s1 c2 s2 c3 s3')
R3 = sp.Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
R2p = sp.Matrix([[1, 0, 0], [0, c2, s2], [0, -s2, c2]])      # ヘルプ DIR = 0 の 2 番目の行列
R2m = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])      # ヘルプ DIR != 0 の 2 番目の行列
R1 = sp.Matrix([[c1, 0, s1], [0, 1, 0], [-s1, 0, c1]])
M0_help = (R3 * R2p * R1).applyfunc(sp.expand)
M1_help = (R3 * R2m * R1).applyfunc(sp.expand)
# tcoord.f の r_ij (dir = .true. と .false.)
M0_code = sp.Matrix([[c1*c3 + s1*s2*s3, -c2*s3, s1*c3 - c1*s2*s3],
                     [-s1*s2*c3 + c1*s3, c2*c3, c1*s2*c3 + s1*s3],
                     [-s1*c2, -s2, c1*c2]])
M1_code = sp.Matrix([[c1*c3 + s1*s2*s3, s1*s2*c3 - c1*s3, s1*c2],
                     [c2*s3, c2*c3, -s2],
                     [-s1*c3 + c1*s2*s3, c1*s2*c3 + s1*s3, c1*c2]])
z33 = sp.zeros(3, 3)
print('DIR = 0: ヘルプの行列 - コードの r_ij = 0 か:', (M0_help - M0_code).applyfunc(sp.expand) == z33)
print('DIR != 0: ヘルプの行列 - コードの r_ij = 0 か:', (M1_help - M1_code).applyfunc(sp.expand) == z33)
Dg = sp.diag(1, -1, -1)
print('DIR != 0 のコードの r_ij = D (DIR = 0 の行列)^T D (D = diag(1,-1,-1)) か:',
      (M1_code - Dg * M0_code.T * Dg).applyfunc(sp.expand) == z33)
sub = {s1: math.sin(0.01), c1: math.cos(0.01), s2: math.sin(-0.02), c2: math.cos(-0.02), s3: math.sin(0.03), c3: math.cos(0.03)}
M0n = np.array(M0_code.subs(sub).evalf(), dtype=float)
M1n = np.array(M1_code.subs(sub).evalf(), dtype=float)
M1hn = np.array(M1_help.subs(sub).evalf(), dtype=float)
print('M0、M1 (コード)、M1 (ヘルプ) は直交行列か: %.1e %.1e %.1e' % (
    np.abs(M0n @ M0n.T - np.eye(3)).max(), np.abs(M1n @ M1n.T - np.eye(3)).max(), np.abs(M1hn @ M1hn.T - np.eye(3)).max()))
print('M1 M0 - 1 の最大 (CHI = 0.01, -0.02, 0.03): コードの M1 %.3e、ヘルプの M1 %.3e (DIR != 0 が DIR = 0 の逆なら 0)' % (
    np.abs(M1n @ M0n - np.eye(3)).max(), np.abs(M1hn @ M0n - np.eye(3)).max()))
print('M0^T M0 - 1 = %.1e (真の逆は転置)' % np.abs(M0n.T @ M0n - np.eye(3)).max())


def rot_help(chi1, chi2, chi3):
    """ヘルプ (DIR = 0) の R3 R2 R1"""
    a, b, c = np.cos(chi1), np.cos(chi2), np.cos(chi3)
    sa, sb, sc = np.sin(chi1), np.sin(chi2), np.sin(chi3)
    R1_ = np.array([[a, 0, sa], [0, 1, 0], [-sa, 0, a]])
    R2_ = np.array([[1, 0, 0], [0, b, sb], [0, -sb, b]])
    R3_ = np.array([[c, -sc, 0], [sc, c, 0], [0, 0, 1]])
    return R3_ @ R2_ @ R1_


def to_plane(r1, P1, p, dl):
    """3 次元の位置 r1 (基準粒子の時刻のもの)・運動量 P1 から、面 z' = 0 を基準にした (x, px, y, py, z, δ) を作る"""
    u1 = P1 / p
    z1 = r1[2] / u1[2]                                  # 面 z' = 0 から測った、進行方向の距離 (先が正)
    X1 = r1 - z1 * u1                                   # 面 z' = 0 を通る点 (X1[2] = 0)
    return np.array([X1[0], P1[0], X1[1], P1[1], z1, dl])


def snapshot(cod):
    x, px, y, py, z, dl = cod
    p = 1.0 + dl
    pz = math.sqrt(p * p - px * px - py * py)
    P = np.array([px, py, pz])
    return p, P, np.array([x, y, 0.0]) + z * P / p      # 面 s = 0 を通る点 + z だけ進行方向に先の点


def line_map(cod, d, chi):
    """独立な写像 (DIR = 0): r' = R (r - d)、P' = R P"""
    p, P, r = snapshot(cod)
    R = rot_help(*chi)
    return to_plane(R @ (r - np.array(d)), R @ P, p, cod[5])


print('\n== (b) 追跡 (tcoord) と CALC の写像 (tcoorde) を、独立な写像と比べる (DIR = 0) ==')
rng = np.random.default_rng(1)
cases = [((1e-3, -2e-3, 5e-3), (0.0, 0.0, 0.0)), ((0.0, 0.0, 0.0), (0.01, -0.02, 0.03)),
         ((1e-3, -2e-3, 5e-3), (0.01, -0.02, 0.03)), ((0.0, 0.0, 0.0), (0.2, 0.0, 0.0)), ((0.0, 0.0, 0.0), (0.0, 0.3, 0.0)),
         ((0.0, 0.0, 0.0), (0.0, 0.0, 0.5)), ((2e-2, 1e-2, -3e-2), (0.1, 0.15, -0.2))]
worst_t = worst_c = 0.0
for d, chi in cases:
    ffs = mk('COORD CO = (DX=%g DY=%g DZ=%g CHI1=%g CHI2=%g CHI3=%g);' % (d + chi))
    el = [e for e in ffs.lat.elements if e.name == 'CO'][0]
    prm = coord.params(el)
    for k in range(5):
        c0 = np.array([rng.normal(0, 1e-3), rng.normal(0, 5e-4), rng.normal(0, 1e-3), rng.normal(0, 5e-4), rng.normal(0, 0.01), rng.normal(0, 2e-3)])
        ref = line_map(c0, d, chi)
        tr = C.track_one(ffs, c0)
        cm = coord.coord_map(c0, prm)[1]
        worst_t = max(worst_t, np.abs(tr - ref).max())
        worst_c = max(worst_c, np.abs(cm - ref).max())
print('独立な写像との最大の差: 追跡 %.2e、CALC の写像 (coord_map) %.2e (7 組 x 5 粒子。pysad の COORD は K64 の tcoord/tcoorde の移植)' % (worst_t, worst_c))

print('\n== (b2) DIR != 0: ヘルプの式 (R3 R2(-s2) R1 (x,y,z) + (DX, -DY, 0)) との比較 ==')
d, chi = (1e-3, -2e-3, 5e-3), (0.01, -0.02, 0.03)
ffs = mk('COORD CO = (DX=%g DY=%g DZ=%g CHI1=%g CHI2=%g CHI3=%g DIR=1);' % (d + chi))
c0 = np.array([1e-3, 5e-4, -1e-3, 3e-4, 0.01, 2e-3])


def line_map_dir1(cod, d, chi, kind):
    p, P, r = snapshot(cod)
    sub_ = {s1: math.sin(chi[0]), c1: math.cos(chi[0]), s2: math.sin(chi[1]), c2: math.cos(chi[1]), s3: math.sin(chi[2]), c3: math.cos(chi[2])}
    if kind == 'help':
        R = np.array(M1_help.subs(sub_).evalf(), dtype=float)
        shift = np.array([d[0], -d[1], 0.0])         # ヘルプ: + (DX, -DY, 0)
    else:
        R = np.array(M1_code.subs(sub_).evalf(), dtype=float)
        shift = np.array([d[0], -d[1], -d[2]])       # tcoord.f: dx1, dy1 = -dy, dz1 = -dz
    return to_plane(R @ r + shift, R @ P, p, cod[5])


tr = C.track_one(ffs, c0)
mc = line_map_dir1(c0, d, chi, 'code')
mh = line_map_dir1(c0, d, chi, 'help')
print('追跡 (DIR = 1)             :', tr)
print('コードの r_ij で作った写像 : 差 %.1e' % np.abs(tr - mc).max())
print('ヘルプの行列で作った写像   : 差 %.1e (x の差 %.2e、px の差 %.2e)' % (np.abs(tr - mh).max(), abs(tr[0] - mh[0]), abs(tr[1] - mh[1])))

print('\n== (c) CALC の 6x6 (tcoorde の解析式) と独立な写像の数値微分、シンプレクティック性 ==')
S = np.zeros((6, 6))
for i in (0, 2, 4):
    S[i, i + 1] = 1
    S[i + 1, i] = -1


def jac_num(f, c0, h=1e-6):
    J = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = h
        J[:, j] = (f(c0 + e) - f(c0 - e)) / (2 * h)
    return J


d, chi = (1e-3, -2e-3, 5e-3), (0.01, -0.02, 0.03)
ffs = mk('COORD CO = (DX=%g DY=%g DZ=%g CHI1=%g CHI2=%g CHI3=%g);' % (d + chi))
el = [e for e in ffs.lat.elements if e.name == 'CO'][0]
prm = coord.params(el)
c0 = np.array([1e-3, 5e-4, -1e-3, 3e-4, 0.05, 2e-3])
Jn = jac_num(lambda c: line_map(c, d, chi), c0)
Ja = coord.coord_map(c0, prm)[0]
print('解析式 (tcoorde) と数値微分の最大の差 = %.2e' % np.abs(Ja - Jn).max())
print('シンプレクティック性 |J^T S J - S| = %.2e (数値微分)、%.2e (解析式)' % (np.abs(Jn.T @ S @ Jn - S).max(), np.abs(Ja.T @ S @ Ja - S).max()))
o = line_map(c0, d, chi)
print('z の更新: z = 0.05 -> %.6f (px, py, delta も変わるので、写像の z は 6 変数すべてに依存する)' % o[4])

print('\n== (d) 小さな角度の 1 次 ==')
eps = 1e-6
zero = np.zeros(6)
for nm, chi_ in (('CHI1', (eps, 0, 0)), ('CHI2', (0, eps, 0)), ('CHI3', (0, 0, eps))):
    o = line_map(zero, (0, 0, 0), chi_)
    print('%s = 1e-6 rad: 設計軌道 (0,0,0,0,0,0) は (x, px, y, py, z) = (%.2e, %.2e, %.2e, %.2e, %.2e) へ' % (nm, *o[:5]))
o = line_map(zero, (1e-3, 0, 0), (0, 0, 0))
print('DX = 1 mm: (x, px) = (%.3e, %.1e)' % (o[0], o[1]))
o = line_map(zero, (0, 0, 1e-2), (0, 0, 0))
print('DZ = 1 cm: 設計軌道は変わらない (x, px, y, py, z) = (%.1e, %.1e, %.1e, %.1e, %.1e)' % tuple(o[:5]))
o = line_map(np.array([0, 1e-3, 0, 0, 0.0, 0]), (0, 0, 1e-2), (0, 0, 0))
print('DZ = 1 cm で px = 1e-3: x = %.4e (= +px/pz * DZ = %.4e)' % (o[0], 1e-3 * 1e-2 / math.sqrt(1 - 1e-6)))

print('\n== (e) DIR = 0 のあとの DIR != 0 は元に戻らない ==')
d, chi = (1e-3, -2e-3, 5e-3), (0.01, -0.02, 0.03)
src = ('MOMENTUM = 3GEV;\nCOORD C0 = (DX=%g DY=%g DZ=%g CHI1=%g CHI2=%g CHI3=%g);\nCOORD C1 = (DX=%g DY=%g DZ=%g CHI1=%g CHI2=%g CHI3=%g DIR=1);\n'
       'MARK ST = ();\nLINE ONE = (ST C0 C1);\nFFS;\nUSE ONE;\n') % (d + chi + d + chi)
ffs, _ = run_sad_string(src, out=io.StringIO())
c0 = np.array([1e-3, 5e-4, -1e-3, 3e-4, 0.01, 2e-3])
tr = C.track_one(ffs, c0)
print('C0 (DIR = 0) の次に C1 (DIR = 1、同じパラメータ): 入力 - 出力の最大 = %.2e (元に戻れば 0)' % np.abs(tr - c0).max())


def inv_map(cod, d, chi):
    p, P, r = snapshot(cod)
    R = rot_help(*chi)
    return to_plane(R.T @ r + np.array(d), R.T @ P, p, cod[5])


print('独立な写像 (DIR = 0) の次に真の逆 (R^T、+d): 入力 - 出力の最大 = %.2e' % np.abs(inv_map(line_map(c0, d, chi), d, chi) - c0).max())
