# -*- coding: utf-8 -*-
r"""SLIM の脱偏極率 (15.15) のモンテカルロ確認 (15_spin.md §6.2)。

玩具の 1 周写像: 軌道 x_{n+1} = M x_n + ξ_n (M: 弱い減衰つきの 3 つの回転、ξ ~ N(0, B))、
スピン S_{n+1} = exp(-[D x_n]×) R0 S_n (K64 の約束)。(15.13) の N で n(x) = n0 + N x (規格化)。
偏極 P(n) = < S_n · n(x_n) > (n は x_n に対する軸) の減りの率 r を、予想 r = (1/2) tr(N B N^T) と比べる。
§6.1 (d) の代数の確認 (自己無撞着) と違い、これは確率過程の平均としての独立な検査である。
"""
import math
import numpy as np
import _common as C

C.start(__file__)


def cm(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def rodrigues_batch(V):
    """exp(-[v]×) を v の配列 (K×3) の各行について作る (K×3×3)"""
    th = np.linalg.norm(V, axis=1)
    th = np.where(th == 0, 1e-300, th)
    k = V / th[:, None]
    K = np.zeros((len(V), 3, 3))
    K[:, 0, 1], K[:, 0, 2] = -k[:, 2], k[:, 1]
    K[:, 1, 0], K[:, 1, 2] = k[:, 2], -k[:, 0]
    K[:, 2, 0], K[:, 2, 1] = -k[:, 1], k[:, 0]
    s = np.sin(th)[:, None, None]
    c1 = (1 - np.cos(th))[:, None, None]
    return np.eye(3)[None] - s * K + c1 * (K @ K)


def run(seed, lam, bscale, dscale, K=2000, T=None):
    rng = np.random.default_rng(seed)
    n0 = np.array([0.0, 1.0, 0.0])
    mus = 2 * math.pi * 0.403
    R0 = np.array([[math.cos(mus), 0, math.sin(mus)], [0, 1, 0], [-math.sin(mus), 0, math.cos(mus)]])   # n0 = ŷ のまわり
    M = np.zeros((6, 6))
    for k, nu in enumerate((0.27, 0.19, 0.03)):
        c, s_ = math.cos(2 * math.pi * nu), math.sin(2 * math.pi * nu)
        M[2 * k:2 * k + 2, 2 * k:2 * k + 2] = (1 - lam) * np.array([[c, s_], [-s_, c]])
    B = np.diag([1, 1, 0.5, 0.5, 0.2, 0.2]) * bscale
    # 平衡の Σ: Σ = M Σ M^T + B
    Sig = np.zeros((6, 6))
    for _ in range(200000):
        Snew = M @ Sig @ M.T + B
        if np.abs(Snew - Sig).max() < 1e-18:
            break
        Sig = Snew
    D = rng.normal(size=(3, 6)) * dscale
    Minv = np.linalg.inv(M)
    Kmat = np.eye(18) - np.kron(Minv.T, R0)
    rhs = cm(n0) @ D @ Minv
    N = np.linalg.lstsq(Kmat, rhs.reshape(-1, order='F'), rcond=None)[0].reshape((3, 6), order='F')
    rate = 0.5 * np.trace(N @ B @ N.T)
    if T is None:
        T = int(round(0.4 / rate / 10)) * 10      # 偏極が約 e^{-0.4} になる周回数
    # モンテカルロ
    L = np.linalg.cholesky(Sig)
    X = rng.normal(size=(K, 6)) @ L.T
    Lb = np.linalg.cholesky(B)

    def axis(X):
        n = n0[None] + X @ N.T
        return n / np.linalg.norm(n, axis=1)[:, None]
    S = axis(X)
    Ps = []
    marks = list(range(0, T + 1, T // 10))
    for t in range(T + 1):
        if t in marks:
            Ps.append((t, float(np.mean(np.sum(S * axis(X), axis=1)))))
        if t == T:
            break
        Rs = rodrigues_batch(X @ D.T)
        S = np.einsum('kij,kj->ki', Rs, S @ R0.T)
        S /= np.linalg.norm(S, axis=1)[:, None]
        X = X @ M.T + rng.normal(size=(K, 6)) @ Lb.T
    ts = np.array([p[0] for p in Ps], float)
    ps = np.array([p[1] for p in Ps])
    # 初期値を 1 に正規化して、ln P の傾き
    y = np.log(ps / ps[0])
    slope = np.sum(ts * y) / np.sum(ts * ts)
    return rate, -slope, Ps, np.sqrt(np.diag(Sig)), np.abs(N).max(), T


print('玩具の 1 周写像 (ν = 0.27, 0.19, 0.03、スピンチューン 0.403、n0 = ŷ)。K = 2000 個のスピン、周回数 T は偏極が約 e^{-0.4} になる長さ。')
print('列: 減衰 λ (1 周あたり)、B の規模、D の規模、予想 r = (1/2) tr(N B Nᵀ) [/周]、測った r (ln P の傾き) [/周]、比、x の rms、|N| の最大')
for seed, lam, bs, ds in ((1, 1e-2, 2e-5, 1.0), (2, 1e-2, 5e-5, 1.0), (3, 3e-3, 1e-5, 1.0), (4, 1e-2, 5e-5, 0.3)):
    rate, meas, Ps, xr, Nmax, T_ = run(seed, lam, bs, ds)
    print('λ = %.0e、B = %.0e、D = %.1f、T = %d: 予想 r = %.3e、測定 r = %.3e、比 = %.3f、x の rms 最大 = %.3f、|N| 最大 = %.2f' % (
        lam, bs, ds, T_, rate, meas, meas / rate, xr.max(), Nmax))
    print('   P(t)/P(0):', ' '.join('%.4f' % (p[1] / Ps[0][1]) for p in Ps))
