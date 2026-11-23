# -*- coding: utf-8 -*-
"""転送行列 (CALC・Emittance 用) と追跡の写像の関係 (04_quad.md §8・§10、式 (4.34)〜(4.38)) を確かめる。

四極 1 台 (3 GeV、L = 0.3 m、K1 = -0.6 (積分強度)、DISFRIN = 1 でフリンジなし)。
(a) 追跡の写像を中心差分 (刻み 1e-6) で微分した 6x6 ヤコビアン J_track を作る (正準運動量の座標。z・δ の行と列を含む)。
(b) CALC の行列 (lattice._el_cod_step の 4x4) と比べる。基準は変分方程式の DOP853 (独立。δ = 0 だけ)。
(c) z の行がシンプレクティック条件 T(5,j) = -T(1,j)T(2,6) + T(2,j)T(1,6) - T(3,j)T(4,6) + T(4,j)T(3,6) で組めることを、J_track で確かめる。
(d) J_track のシンプレクティック性 |M^T J M - J|。
K64 の動作 (k64) と v2 の動作の両方で行う。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import lattice as lat

L, K1 = 0.3, -0.6
ffs, el = C.make_quad(L, K1, 'DISFRIN=1')
J6 = np.zeros((6, 6))
for i in range(3):
    J6[2 * i, 2 * i + 1] = 1.0
    J6[2 * i + 1, 2 * i] = -1.0


def jac_track(c, h=1e-6):
    """追跡の写像の中心差分ヤコビアン 6x6 (入出力とも正準運動量の (x, px, y, py, z, δ))"""
    c = np.array(c, float)
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (C.track_one(ffs, cp) - C.track_one(ffs, cm)) / (2 * h)
    return M


def ref_jac(c, Ln, KL):
    """独立な基準: 運動学的な運動方程式 (δ = 0) の変分方程式 4x4。KL < 0 は -90 度回して |K1| で通し、戻す"""
    rot = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float)
    k = abs(KL) / Ln
    c4 = np.array(c[:4], float)
    if KL < 0:
        c4 = rot @ c4

    def rhs(s, v):
        x, px, y, py = v[:4]
        M = v[4:].reshape(4, 4)
        pz = math.sqrt(1.0 - px * px - py * py)
        A = np.zeros((4, 4))
        A[0, 1] = 1 / pz + px * px / pz**3
        A[0, 3] = px * py / pz**3
        A[2, 3] = 1 / pz + py * py / pz**3
        A[2, 1] = px * py / pz**3
        A[1, 0] = -k
        A[3, 2] = k
        return np.concatenate([[px / pz, -k * x, py / pz, k * y], (A @ M).ravel()])
    sol = solve_ivp(rhs, (0.0, Ln), np.concatenate([c4, np.eye(4).ravel()]), method='DOP853', rtol=1e-13, atol=1e-16)
    M = sol.y[4:, -1].reshape(4, 4)
    return rot.T @ M @ rot if KL < 0 else M


cases = [('設計軌道 (δ = 0)', [0, 0, 0, 0, 0, 0.0]),
         ('設計軌道、δ = 0.01', [0, 0, 0, 0, 0, 0.01]),
         ('3 mm・2 mrad ずれた軌道 (δ = 0)', [0.003, 0.002, 0, 0, 0, 0.0]),
         ('x, y, px, py がずれた軌道 (δ = 0.003)', [0.003, 0.002, 0.001, -0.001, 0, 0.003])]
for mode in ('k64', 'v2'):
    C.set_mode(mode)
    print('\n########## %s ##########' % ('K64 の動作 (k64: 2 次の Strang)' if mode == 'k64' else 'v2 (4 次の Yoshida)'))
    for name, c in cases:
        c = np.array(c, float)
        Jt = jac_track(c)
        Mc, cod_out = lat._el_cod_step(el, c, c[5])
        d44 = np.abs(Mc[:4, :4] - Jt[:4, :4]).max()
        d6 = np.abs(Mc[:4, 5] - Jt[:4, 5]).max()
        s = '  %-38s CALC 4x4 と J_track の差 %.1e、δ の列の差 %.1e' % (name, d44, d6)
        if c[5] == 0.0:
            Mr = ref_jac(c, L, K1 * 1.0)
            s += '、CALC と変分方程式 %.1e、J_track と変分方程式 %.1e' % (np.abs(Mc[:4, :4] - Mr).max(), np.abs(Jt[:4, :4] - Mr).max())
        print(s)
    # (c) z の行
    c = np.array([0.003, 0.002, 0.001, -0.001, 0, 0.003])
    Jt = jac_track(c)
    z_form = np.array([-Jt[0, j] * Jt[1, 5] + Jt[1, j] * Jt[0, 5] - Jt[2, j] * Jt[3, 5] + Jt[3, j] * Jt[2, 5] for j in range(4)])
    print('  z の行 T(5,1..4): 追跡の数値微分 %s' % np.array2string(Jt[4, :4], precision=6))
    print('                    シンプレクティック条件で組んだ値 %s、差 %.1e' % (np.array2string(z_form, precision=6), np.abs(z_form - Jt[4, :4]).max()))
    # (d) シンプレクティック性
    print('  |J_track^T J J_track - J| = %.1e (中心差分の精度 ~1e-9 まで)、CALC の 4x4 は %.1e' % (
        np.abs(Jt.T @ J6 @ Jt - J6).max(),
        np.abs(lat._el_cod_step(el, c, c[5])[0][:4, :4].T @ J6[:4, :4] @ lat._el_cod_step(el, c, c[5])[0][:4, :4] - J6[:4, :4]).max()))
C.set_mode('v2')

print('\n########## (e) emit6d の 1 刻みの本体 B (解析ヤコビアン、v2) の z の行 ##########')
from pysad import emit6d
C.set_mode('v2')
c = np.array([0.003, 0.002, 0.001, -0.001, 0.0, 0.003])
JB = emit6d._jac_tsolque_body(c, 0.15, 2.0, 1e-6)
z_form = np.array([-JB[0, j] * JB[1, 5] + JB[1, j] * JB[0, 5] - JB[2, j] * JB[3, 5] + JB[3, j] * JB[2, 5] for j in range(4)])
print('  J_B の z の行 %s' % np.array2string(JB[4, :4], precision=6))
print('  シンプレクティック条件の値 %s、差 %.1e' % (np.array2string(z_form, precision=6), np.abs(z_form - JB[4, :4]).max()))
print('  |J_B^T J J_B - J| (6x6) = %.1e' % np.abs(JB.T @ J6 @ JB - J6).max())
JD = emit6d._jac_tqente(c, 0.15)
print('  D (tqente 相当) の 6x6: pysad の解析ヤコビアンは z の行の (5,2)・(5,4) を持たない (K64 の tqente の trans2(5,2) = trans2(1,6) に当たる、傾きの 3 乗の大きさ)。|J_D^T J J_D - J| = %.1e (その分)' % np.abs(JD.T @ J6 @ JD - J6).max())

print('\n########## (f) 線形フリンジの係数 F1 の符号: 追跡と CALC (K1 = 0.6、L = 0.3 m、FRINGE = 3、DISFRIN = 1、設計軌道) ##########')
for F1 in (0.1, -0.1):
    ffs_f, el_f = C.make_quad(0.3, 0.6, 'FRINGE=3 DISFRIN=1 F1=%g' % F1)
    for mode in ('k64', 'v2'):
        C.set_mode(mode)
        c = np.zeros(6)
        Jt = np.zeros((4, 4))
        for j in range(4):
            cp, cm = c.copy(), c.copy()
            cp[j] += 1e-6
            cm[j] -= 1e-6
            Jt[:, j] = (C.track_one(ffs_f, cp)[:4] - C.track_one(ffs_f, cm)[:4]) / 2e-6
        Mc, _ = lat._el_cod_step(el_f, c, 0.0)
        print('  F1 = %+.1f、%s: CALC 4x4 と追跡のヤコビアンの差 %.2e (追跡の x の行 %s)' % (F1, mode, np.abs(Mc[:4, :4] - Jt).max(), np.array2string(Jt[0, :2], precision=6)))
C.set_mode('v2')
