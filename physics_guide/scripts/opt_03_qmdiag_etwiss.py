# -*- coding: utf-8 -*-
"""12 章 (12_linear_optics.md) の補完: qmdiag の代数的な証明、特別な場合、拡張 Twiss (etwiss2ri・tfetwiss)。

(a) (12.6) が Riccati 式 (12.5) を満たすことの代数的な証明に使う恒等式を、ランダムなシンプレクティック行列で確かめる:
    X P - K X = tau Q、X N X = det(X) Q、det N = det Q、det P + det Q = 1。E = Q (1 + kappa tau - kappa² det X)。
    さらに c² = 1/(1 + kappa² det X)、tr A - tr B = tau - 2 kappa det X = sgn(tau) sqrt(trab2)。
(b) qmdiag の特別な場合 (trab2 ≈ 0): Fortran の分岐を numpy に写し、縮退した (等しいチューンで結合した) 行列で、
    c = 1 の別の式が R M R^-1 をブロック対角にすることを確かめる。通常の式では c² が発散する。
(c) 拡張 Twiss: K64 の etwiss2ri (twiss パラメータ → 6×6 行列) が、ヘルプの V = B R6 H (位相の P なし) と一致し、
    tfetwiss (行列 → パラメータ) がその逆であること。
"""
import math
import numpy as np
from scipy.linalg import expm
import _common as C

C.start(__file__)
rng = np.random.default_rng(11)
J2 = np.array([[0, 1], [-1, 0.0]])
J4 = np.kron(np.eye(2), J2)
bar = lambda A: J2 @ A.T @ J2            # 2×2: bar(A) = -adj(A)


def rand_sympl(scale=0.4):
    S = rng.normal(size=(4, 4))
    S = (S + S.T) / 2
    return expm(J4 @ S * scale)


print('== (a) Riccati 式の代数的な証明 (数値で恒等式を確認) ==')
mx = np.zeros(8)
for _ in range(200):
    M = rand_sympl()
    P, N, Q, K = M[:2, :2], M[:2, 2:], M[2:, :2], M[2:, 2:]
    X = Q - bar(N)
    tau = np.trace(P - K)
    dX = np.linalg.det(X)
    mx[0] = max(mx[0], np.abs(X @ P - K @ X - tau * Q).max())
    mx[1] = max(mx[1], np.abs(X @ N @ X - dX * Q).max())
    mx[2] = max(mx[2], abs(np.linalg.det(N) - np.linalg.det(Q)))
    mx[3] = max(mx[3], abs(np.linalg.det(P) + np.linalg.det(Q) - 1))
    trab2 = tau**2 + 4 * dX
    if trab2 > 1e-8:
        s = 1.0 if tau >= 0 else -1.0
        kap = -2.0 / (tau + s * math.sqrt(trab2))
        rh = kap * X
        E = Q + rh @ P - K @ rh - rh @ N @ rh
        mx[4] = max(mx[4], np.abs(E).max())
        c2 = 1.0 / (1.0 + kap**2 * dX)
        mx[5] = max(mx[5], abs(c2 - 0.5 * (1 + abs(tau) / math.sqrt(trab2))))
        A = P + bar(rh) @ Q
        B = K + rh @ N
        mx[6] = max(mx[6], abs(np.trace(A) - np.trace(B) - s * math.sqrt(trab2)))
        mx[7] = max(mx[7], abs(dX * kap**2 - tau * kap - 1))
print('  200 個のランダムなシンプレクティック行列 (4×4)、2×2 ブロック P N Q K、X = Q - bar(N)、tau = tr(P - K)、bar(A) = J A^T J:')
print('   max|X P - K X - tau Q|            = %.1e   (恒等式 (i))' % mx[0])
print('   max|X N X - det(X) Q|             = %.1e   (恒等式 (ii))' % mx[1])
print('   max|det N - det Q|                = %.1e' % mx[2])
print('   max|det P + det Q - 1|            = %.1e' % mx[3])
print('   max|Riccati の残差 E| (kappa は (12.6)) = %.1e' % mx[4])
print('   max|c² - (1 + |tau|/sqrt(trab2))/2|  = %.1e   (c² = 1/(1 + kappa² det X))' % mx[5])
print('   max|trA - trB - sgn(tau) sqrt(trab2)| = %.1e' % mx[6])
print('   max|det(X) kappa² - tau kappa - 1|   = %.1e   (kappa の 2 次方程式)' % mx[7])


# ------------------------------------------------------------------------------------------ (b)
def qmdiag(tm):
    """qmdiag.f の枝を写したもの。tm は 4×4 (tm[i-1][j-1] = tm(i,j))。戻り値 (r1, r2, r3, r4, c, stab, 分岐の名前)"""
    t = lambda i, j: tm[i - 1, j - 1]
    trpq = t(1, 1) + t(2, 2) - t(3, 3) - t(4, 4)
    trab2 = 4.0 * ((t(3, 1) + t(2, 4)) * (t(4, 2) + t(1, 3)) + (t(1, 4) - t(3, 2)) * (t(4, 1) - t(2, 3))) + trpq**2
    if trab2 >= 1e-14:
        c = math.sqrt(0.5 * (1.0 + abs(trpq) / math.sqrt(trab2)))
        aa = -1.0 / c / math.sqrt(trab2) if trpq >= 0 else 1.0 / c / math.sqrt(trab2)
        return (aa * (t(3, 1) + t(2, 4)), aa * (t(3, 2) - t(1, 4)), aa * (t(4, 1) - t(2, 3)), aa * (t(4, 2) + t(1, 3)), c, True, '通常')
    if abs(trab2) > 1e-14:
        return (0, 0, 0, 0, 0, False, '和共鳴 (trab2 < 0)')
    detn = t(1, 3) * t(2, 4) - t(1, 4) * t(2, 3)
    if abs(detn) >= 1e-14:
        trqk = t(3, 3) - t(4, 4)
        if abs(trqk) >= 1e-7:
            trqj = t(4, 3) - t(3, 4)
            trqkj2 = (trqk + trqj) * (trqk - trqj)
            if abs(trqkj2) >= 1e-14:
                c, bbb, ccc = 1.0, -detn / (trqj + trqk), 0.0
            elif abs(trqj) >= 1e-7:
                c, bbb, ccc = 1.0, -0.5 * detn / trqj, 0.0
            else:
                raise RuntimeError('logic broken')
            aaa = (detn * (1.0 - 2.0 * c * c) / c - bbb * trqj) / trqk
            name = 'trab2 = 0, det N ≠ 0, tr(QK) ≠ 0'
        else:
            trqj = t(4, 3) - t(3, 4)
            if abs(trqj) >= 1e-7:
                c, aaa, bbb, ccc = 1.0, detn / trqj, -detn / trqj, 0.0
            else:
                trql = t(3, 4) + t(4, 3)
                if abs(trql) >= 1e-7:
                    c, aaa, bbb, ccc = 1.0, 0.0, detn / trql, -detn / trql
                else:
                    c = 1.0 / math.sqrt(2.0)
                    if detn > 0:
                        aaa, bbb, ccc = 0.0, math.sqrt(detn / 2.0), 0.0
                    else:
                        aaa, bbb, ccc = math.sqrt(-detn / 2.0), 0.0, 0.0
            name = 'trab2 = 0, det N ≠ 0, tr(QK) = 0'
        xx11, xx12, xx21, xx22 = aaa, bbb + ccc, -bbb + ccc, -aaa
        r1 = (xx11 * t(2, 4) - xx12 * t(2, 3)) / detn
        r2 = (-xx11 * t(1, 4) + xx12 * t(1, 3)) / detn
        r3 = (xx21 * t(2, 4) - xx22 * t(2, 3)) / detn
        r4 = (-xx21 * t(1, 4) + xx22 * t(1, 3)) / detn
        return (r1, r2, r3, r4, c, True, name)
    trqk = t(3, 3) - t(4, 4)
    trqj = t(4, 3) - t(3, 4)
    trql = t(4, 3) + t(3, 4)
    if abs(trqk) >= 1e-7:
        aaa, bbb, ccc = 0.0, 0.0, 1.0 / trqk
    elif abs(trqj) >= 1e-7:
        aaa, bbb, ccc = -1.0 / trqj, 0.0, 0.0
    elif abs(trql) >= 1e-7:
        aaa, bbb, ccc = 0.0, -1.0 / trql, 0.0
    else:
        trp = t(1, 1) + t(2, 2)
        p1, p3 = t(1, 1) - 0.5 * trp, t(2, 1)
        if abs(p3) >= 1e-7:
            return (0.0, t(2, 4) / p3, 0.0, -t(2, 3) / p3, 1.0, True, 'trab2 = 0, det N = 0, 全部 0 (p3 ≠ 0)')
        if abs(p1) >= 1e-7:
            return (t(2, 4) / p1, 0.0, -t(2, 3) / p1, 0.0, 1.0, True, 'trab2 = 0, det N = 0, 全部 0 (p1 ≠ 0)')
        return (1.0, 0.0, 0.0, 1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), True, 'trab2 = 0, det N = 0, 全部 0')
    yy11, yy12, yy21, yy22 = aaa + bbb, ccc, ccc, aaa - bbb
    return (yy11 * (-t(2, 3)) + yy12 * (-t(2, 4)), yy11 * t(1, 3) + yy12 * t(1, 4),
            yy21 * (-t(2, 3)) + yy22 * (-t(2, 4)), yy21 * t(1, 3) + yy22 * t(1, 4), 1.0, True, 'trab2 = 0, det N = 0, tr(QK)/tr(QJ)/tr(QL) ≠ 0')


def build_R(r1, r2, r3, r4, c, third_form=False):
    r = np.array([[r1, r2], [r3, r4]])
    return np.block([[c * np.eye(2), bar(r).T if False else J2 @ r.T @ J2], [r, c * np.eye(2)]])


def offdiag(M, R):
    D = R @ M @ np.linalg.inv(R)
    return max(np.abs(D[:2, 2:]).max(), np.abs(D[2:, :2]).max())


print('\n== (b) qmdiag の特別な場合 (trab2 ≈ 0) ==')
# 縮退: 2 つのモードの回転角が等しい (tr A = tr B)。結合した行列 M = S^-1 diag(A, B) S
mu = 0.9
Rot = lambda m: np.array([[math.cos(m), math.sin(m)], [-math.sin(m), math.cos(m)]])
Dg = np.block([[Rot(mu), np.zeros((2, 2))], [np.zeros((2, 2)), Rot(mu)]])
# 結合させる変換 (正規形に近いもの): 小さいシンプレクティック行列
S = rand_sympl(0.12)
Mc = np.linalg.inv(S) @ Dg @ S
t = lambda i, j: Mc[i - 1, j - 1]
trpq = t(1, 1) + t(2, 2) - t(3, 3) - t(4, 4)
trab2 = 4 * ((t(3, 1) + t(2, 4)) * (t(4, 2) + t(1, 3)) + (t(1, 4) - t(3, 2)) * (t(4, 1) - t(2, 3))) + trpq**2
print('  縮退した結合行列: trab2 = %.1e (= 4 det X + tau²)、tau = %.1e' % (trab2, trpq))
res = qmdiag(Mc)
print('  qmdiag の枝: %s、c = %.6f、(r1..r4) = (%.5f, %.5f, %.5f, %.5f)' % (res[6], res[4], *res[:4]))
R = build_R(*res[:5])
print('  R M R^-1 の非対角ブロックの最大 = %.1e、|R^T J4 R - J4| = %.1e、c² + det r - 1 = %.1e' % (
    offdiag(Mc, R), np.abs(R.T @ J4 @ R - J4).max(), res[4]**2 + (res[0] * res[3] - res[1] * res[2]) - 1))
# 通常の式 (12.6) を無理に使うと
tau_ = trpq
Xn = np.array([[t(3, 1) + t(2, 4), t(3, 2) - t(1, 4)], [t(4, 1) - t(2, 3), t(4, 2) + t(1, 3)]])
print('  通常の式 (12.6): c² = (1 + |tau|/sqrt(trab2))/2 は trab2 → 0 で %s (|tau|/sqrt(trab2) = %.2e)' % ('発散' if True else '', abs(tau_) / math.sqrt(abs(trab2) + 1e-300)))
print('  (理由: kappa の 2 次方程式 det(X) kappa² - tau kappa - 1 = 0 の判別式 trab2 = tau² + 4 det X が 0 で重根 kappa = tau/(2 det X) = -2/tau、')
print('   c² = 1/(2 + tau kappa) = 1/0。実の r が c = 1 で取れる別の式に切り替える。)')
# detn = 0 の場合 (結合なし、チューンが等しい)
Mu = np.block([[Rot(mu), np.zeros((2, 2))], [np.zeros((2, 2)), Rot(mu)]])
resu = qmdiag(Mu)
print('  結合なし・等しいチューン (N = Q = 0): 枝 %s、c = %.3f、r = (%.1e, %.1e, %.1e, %.1e)、非対角ブロック = %.1e' % (
    resu[6], resu[4], *resu[:4], offdiag(Mu, build_R(*resu[:5]))))
# 和共鳴 (trab2 < 0): nu_x + nu_y が整数の近く (mu_y = -mu_x) で、x-y の結合 (H に k x y) を足す
Dsum = np.block([[Rot(0.9), np.zeros((2, 2))], [np.zeros((2, 2)), Rot(-0.9)]])
Sc = np.zeros((4, 4))
Sc[0, 2] = Sc[2, 0] = 1.0                      # H = (1/2) u^T Sc u = x y
for kk in (0.0, 0.1, 0.3):
    Ms = Dsum @ expm(J4 @ Sc * kk)
    rs = qmdiag(Ms)
    t_ = lambda i, j: Ms[i - 1, j - 1]
    tp = t_(1, 1) + t_(2, 2) - t_(3, 3) - t_(4, 4)
    tb = 4 * ((t_(3, 1) + t_(2, 4)) * (t_(4, 2) + t_(1, 3)) + (t_(1, 4) - t_(3, 2)) * (t_(4, 1) - t_(2, 3))) + tp**2
    print('  和共鳴付近 (mu_x = 0.9、mu_y = -0.9)、結合の強さ k = %.1f: trab2 = %+.3e、枝 = %s、固有値の絶対値の最大 = %.4f' % (kk, tb, rs[6], np.abs(np.linalg.eigvals(Ms)).max()))

# ------------------------------------------------------------------------------------------ (c)
print('\n== (c) 拡張 Twiss: etwiss2ri (K64) と V = B R6 H (ヘルプ)、tfetwiss (逆) ==')


def etwiss2ri(r1, r2, r3, r4, bx, ax, by, ay, bz, az_, hx, hy):
    """tffs.f etwiss2ri の normal = .true. の枝。hx, hy は 2×2 (hx11 hx12; hx21 hx22)"""
    detr = r1 * r4 - r2 * r3
    amu = math.sqrt(1 - detr)
    sqrbx, sqrby, sqrbz = math.sqrt(bx), math.sqrt(by), math.sqrt(bz)
    hx11, hx12, hx21, hx22 = hx[0, 0], hx[0, 1], hx[1, 0], hx[1, 1]
    hy11, hy12, hy21, hy22 = hy[0, 0], hy[0, 1], hy[1, 0], hy[1, 1]
    dethx = hx11 * hx22 - hx21 * hx12
    dethy = hy11 * hy22 - hy21 * hy12
    azz = math.sqrt(1 - dethx - dethy)
    aa = 1 / (1 + azz)
    amux = 1 - dethx * aa
    amuy = 1 - dethy * aa
    h = np.zeros((4, 6))
    h[0, 0] = amux
    h[0, 2] = (hx12 * hy21 - hx11 * hy22) * aa
    h[0, 3] = (-hx12 * hy11 + hx11 * hy12) * aa
    h[0, 4], h[0, 5] = -hx11, -hx12
    h[1, 1] = amux
    h[1, 2] = (hx22 * hy21 - hx21 * hy22) * aa
    h[1, 3] = (-hx22 * hy11 + hx21 * hy12) * aa
    h[1, 4], h[1, 5] = -hx21, -hx22
    h[2, 0] = h[1, 3]
    h[2, 1] = -h[0, 3]
    h[2, 2] = amuy
    h[2, 4], h[2, 5] = -hy11, -hy12
    h[3, 0] = -h[1, 2]
    h[3, 1] = h[0, 2]
    h[3, 3] = amuy
    h[3, 4], h[3, 5] = -hy21, -hy22
    ria = np.zeros((6, 6))
    ria[4, :] = [hx22, -hx12, hy22, -hy12, azz, 0]
    ria[5, :] = [-hx21, hx11, -hy21, hy11, 0, azz]
    br = np.array([
        [amu / sqrbx, 0, -r4 / sqrbx, r2 / sqrbx],
        [amu * ax / sqrbx, amu * sqrbx, r3 * sqrbx - r4 * ax / sqrbx, r2 * ax / sqrbx - r1 * sqrbx],
        [r1 / sqrby, r2 / sqrby, amu / sqrby, 0],
        [r1 * ay / sqrby + r3 * sqrby, r2 * ay / sqrby + r4 * sqrby, amu * ay / sqrby, amu * sqrby]])
    ria[:4, :] = br @ h
    ria[4, :] = ria[4, :] / sqrbz
    ria[5, :] = ria[5, :] * sqrbz + ria[4, :] * az_
    return ria


def help_V(r1, r2, r3, r4, bx, ax, by, ay, bz, az_, hx, hy):
    """ヘルプの V = B R6 H (位相の P なし)。(12.10)"""
    detr = r1 * r4 - r2 * r3
    mu_ = math.sqrt(1 - detr)
    r = np.array([[r1, r2], [r3, r4]])
    R = np.block([[mu_ * np.eye(2), bar(r).T if False else J2 @ r.T @ J2], [r, mu_ * np.eye(2)]])
    R6 = np.eye(6)
    R6[:4, :4] = R
    a = math.sqrt(1 - np.linalg.det(hx) - np.linalg.det(hy))
    I2 = np.eye(2)
    H = np.block([
        [(1 - np.linalg.det(hx) / (1 + a)) * I2, hx @ J2 @ hy.T @ J2 / (1 + a), -hx],
        [hy @ J2 @ hx.T @ J2 / (1 + a), (1 - np.linalg.det(hy) / (1 + a)) * I2, -hy],
        [-J2 @ hx.T @ J2, -J2 @ hy.T @ J2, a * I2]])
    Bm = lambda b, al: np.array([[1 / math.sqrt(b), 0], [al / math.sqrt(b), math.sqrt(b)]])
    B6 = np.zeros((6, 6))
    B6[:2, :2] = Bm(bx, ax)
    B6[2:4, 2:4] = Bm(by, ay)
    B6[4:, 4:] = Bm(bz, az_)
    return B6 @ R6 @ H


mxd = 0.0
for _ in range(5):
    r = rng.normal(0, 0.2, 4)
    bx, by, bz = rng.uniform(1, 8, 3)
    ax, ay, az_ = rng.normal(0, 0.5, 3)
    hx = rng.normal(0, 0.1, (2, 2))
    hy = rng.normal(0, 0.1, (2, 2))
    V1 = etwiss2ri(*r, bx, ax, by, ay, bz, az_, hx, hy)
    V2 = help_V(*r, bx, ax, by, ay, bz, az_, hx, hy)
    mxd = max(mxd, np.abs(V1 - V2).max())
print('  etwiss2ri (K64) - V = B R6 H (ヘルプ (12.10) の式から作ったもの) の最大 = %.1e  (5 組のランダムなパラメータ、P なし)' % mxd)
Vs = V1
print('  この V はシンプレクティック: |V^T J6 V - J6| = %.1e (J6 = diag(J, J, J))' % np.abs(Vs.T @ np.kron(np.eye(3), J2) @ Vs - np.kron(np.eye(3), J2)).max())


def tfetwiss(r):
    """tffs.f tfetwiss (normal の枝) の ex, epx, ey, epy, zx, zpx, zy, zpy, ax, bx, ay, by, r1..r4 を返す"""
    az = math.sqrt(r[4, 4] * r[5, 5] - r[5, 4] * r[4, 5])
    uz11, uz12, uz21, uz22 = r[4, 4] / az, r[4, 5] / az, r[5, 4] / az, r[5, 5] / az
    hx11 = r[5, 1] * uz11 - r[4, 1] * uz21
    hx12 = r[5, 1] * uz12 - r[4, 1] * uz22
    hx21 = -r[5, 0] * uz11 + r[4, 0] * uz21
    hx22 = -r[5, 0] * uz12 + r[4, 0] * uz22
    hy11 = r[5, 3] * uz11 - r[4, 3] * uz21
    hy12 = r[5, 3] * uz12 - r[4, 3] * uz22
    hy21 = -r[5, 2] * uz11 + r[4, 2] * uz21
    hy22 = -r[5, 2] * uz12 + r[4, 2] * uz22
    f = 1 / (1 + az)
    ax = 1 - (hx11 * hx22 - hx21 * hx12) * f
    ay = 1 - (hy11 * hy22 - hy21 * hy12) * f
    hi = np.zeros((6, 6))
    hi[1, 1] = ax
    hi[0, 1] = 0
    hi[3, 1] = (hx12 * hy21 - hx11 * hy22) * f
    hi[2, 1] = -(-hx12 * hy11 + hx11 * hy12) * f
    hi[5, 1] = -hx11
    hi[4, 1] = hx12
    hi[1, 0] = 0
    hi[0, 0] = ax
    hi[3, 0] = -(hx22 * hy21 - hx21 * hy22) * f
    hi[2, 0] = (-hx22 * hy11 + hx21 * hy12) * f
    hi[5, 0] = hx21
    hi[4, 0] = -hx22
    hi[1, 3] = hi[2, 0]
    hi[0, 3] = -hi[2, 1]
    hi[3, 3] = ay
    hi[2, 3] = 0
    hi[5, 3] = -hy11
    hi[4, 3] = hy12
    hi[1, 2] = -hi[3, 0]
    hi[0, 2] = hi[3, 1]
    hi[3, 2] = 0
    hi[2, 2] = ay
    hi[5, 2] = hy21
    hi[4, 2] = -hy22
    hi[1, 5] = hx22
    hi[0, 5] = hx21
    hi[3, 5] = hy22
    hi[2, 5] = hy21
    hi[5, 5] = az
    hi[4, 5] = 0
    hi[1, 4] = hx12
    hi[0, 4] = hx11
    hi[3, 4] = hy12
    hi[2, 4] = hy11
    hi[5, 4] = 0
    hi[4, 4] = az
    hi = r @ hi
    detm = (hi[0, 0] * hi[1, 1] - hi[1, 0] * hi[0, 1] + hi[2, 2] * hi[3, 3] - hi[3, 2] * hi[2, 3]) * 0.5
    axy = math.sqrt(detm)
    r11 = (hi[3, 3] * hi[2, 0] - hi[2, 3] * hi[3, 0]) / axy
    r12 = (hi[3, 3] * hi[2, 1] - hi[2, 3] * hi[3, 1]) / axy
    r21 = (-hi[3, 2] * hi[2, 0] + hi[2, 2] * hi[3, 0]) / axy
    r22 = (-hi[3, 2] * hi[2, 1] + hi[2, 2] * hi[3, 1]) / axy
    crx = math.hypot(hi[0, 1], hi[1, 1])
    cx, sx = hi[1, 1] / crx, -hi[0, 1] / crx
    bx21 = (-sx * hi[0, 0] + cx * hi[1, 0]) / axy
    bx22 = crx / axy
    cry = math.hypot(hi[2, 3], hi[3, 3])
    cy, sy = hi[3, 3] / cry, -hi[2, 3] / cry
    by21 = (-sy * hi[2, 2] + cy * hi[3, 2]) / axy
    by22 = cry / axy
    out = dict(r=(r11, r12, r21, r22), ax=bx21 * bx22, bx=bx22**2, ay=by21 * by22, by=by22**2,
               ex=axy * hx12 - r22 * hy12 + r12 * hy22, epx=axy * hx22 + r21 * hy12 - r11 * hy22,
               ey=axy * hy12 + r11 * hx12 + r12 * hx22, epy=axy * hy22 + r21 * hx12 + r22 * hx22,
               zx=axy * hx11 - r22 * hy11 + r12 * hy21, zpx=axy * hx21 + r21 * hy11 - r11 * hy21,
               zy=axy * hy11 + r11 * hx11 + r12 * hx21, zpy=axy * hy21 + r21 * hx11 + r22 * hx21)
    return out


r = np.array([0.12, -0.05, 0.08, 0.15])
bx, by, bz = 4.2, 2.7, 3.3
ax, ay, az_ = 0.6, -0.4, 0.2
hx = np.array([[0.05, 0.12], [-0.03, 0.07]])
hy = np.array([[0.02, -0.04], [0.06, 0.03]])
Vr = etwiss2ri(*r, bx, ax, by, ay, bz, az_, hx, hy)
o = tfetwiss(Vr)
# 期待: R [Hx; Hy] の列。ex, epx, ey, epy = R (Hx[:, 1]; Hy[:, 1]) 、zx.. = R (Hx[:, 0]; Hy[:, 0])
detr = r[0] * r[3] - r[1] * r[2]
mu_ = math.sqrt(1 - detr)
rr_ = np.array([[r[0], r[1]], [r[2], r[3]]])
Rm = np.block([[mu_ * np.eye(2), J2 @ rr_.T @ J2], [rr_, mu_ * np.eye(2)]])
HH = np.vstack([hx, hy])
expect = Rm @ HH
print('  tfetwiss(etwiss2ri(…)) の往復: β_x、α_x、β_y、α_y の差 = %.1e、r の差 = %.1e' % (
    max(abs(o['bx'] - bx), abs(o['ax'] - ax), abs(o['by'] - by), abs(o['ay'] - ay)), np.abs(np.array(o['r']) - r).max()))
print('  分散 (EX, EPX, EY, EPY) と R [H_x; H_y] の第 2 列の差 = %.1e、(ZX, ZPX, ZY, ZPY) と第 1 列の差 = %.1e' % (
    np.abs(np.array([o['ex'], o['epx'], o['ey'], o['epy']]) - expect[:, 1]).max(),
    np.abs(np.array([o['zx'], o['zpx'], o['zy'], o['zpy']]) - expect[:, 0]).max()))

# ------------------------------------------------------------------------------------------ (d)
import sympy as sp
print('\n== (d) Twiss の伝搬 (12.9)、位相の進み、薄いレンズの色収差 (sympy) ==')
mu_s, b1, a1 = sp.symbols('mu beta1 alpha1', real=True)
t11, t12, t21, t22 = sp.symbols('t11 t12 t21 t22', real=True)
g1 = (1 + a1**2) / b1
A1 = sp.Matrix([[sp.cos(mu_s) + a1 * sp.sin(mu_s), b1 * sp.sin(mu_s)], [-g1 * sp.sin(mu_s), sp.cos(mu_s) - a1 * sp.sin(mu_s)]])
print('  det A1 = %s' % sp.simplify(A1.det()))
T = sp.Matrix([[t11, t12], [t21, t22]])
detT = t11 * t22 - t12 * t21
A2 = (T * A1 * T.inv()).applyfunc(sp.simplify)
beta2 = sp.simplify(A2[0, 1] / sp.sin(mu_s))
alpha2 = sp.simplify((A2[0, 0] - A2[1, 1]) / (2 * sp.sin(mu_s)))
beta2_f = t11**2 * b1 - 2 * t11 * t12 * a1 + t12**2 * g1
alpha2_f = -t11 * t21 * b1 + (t11 * t22 + t12 * t21) * a1 - t12 * t22 * g1
print('  beta2 (A2 の (1,2) 成分 / sin mu) - (12.9) = %s' % sp.simplify((beta2 - beta2_f).subs(t22, (1 + t12 * t21) / t11)))
print('  alpha2 ((A2_11 - A2_22)/(2 sin mu)) - (12.9) = %s' % sp.simplify((alpha2 - alpha2_f).subs(t22, (1 + t12 * t21) / t11)))
Bm1 = sp.Matrix([[1 / sp.sqrt(b1), 0], [a1 / sp.sqrt(b1), sp.sqrt(b1)]])
b2s = sp.symbols('beta2s', positive=True)
a2s = sp.symbols('alpha2s', real=True)
Bm2 = sp.Matrix([[1 / sp.sqrt(b2s), 0], [a2s / sp.sqrt(b2s), sp.sqrt(b2s)]])
Rn = (Bm2 * T * Bm1.inv()).applyfunc(sp.simplify)
print('  tan(Δψ) = (B2 T B1^-1)_12 / (B2 T B1^-1)_11 = %s' % sp.simplify(Rn[0, 1] / Rn[0, 0]))
print('   (12.9 の下の式: t12/(beta1 t11 - alpha1 t12))')
# 薄いレンズ: 1 周行列 M0 (レンズの直後から)、M = M0 [[1,0],[-dk,1]]
m11, m12, m21, m22, dk = sp.symbols('m11 m12 m21 m22 dk')
M0 = sp.Matrix([[m11, m12], [m21, m22]])
Mn = M0 * sp.Matrix([[1, 0], [-dk, 1]])
print('  tr(M0 [[1,0],[-dk,1]]) - tr M0 = %s  → cos mu の変化 = -dk m12/2、m12 = beta sin mu より  d mu = dk beta/2' % sp.simplify(Mn.trace() - M0.trace()))
