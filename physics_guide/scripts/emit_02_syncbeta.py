# -*- coding: utf-8 -*-
"""平衡エミッタンス (14_emittance.md) の補足の確認: 有限のシンクロトロンチューンの補正、tnorm・tsymp、btr、結合のあるリング。

リング: 14 章と同じ 8 セルの FODO (3 GeV、扇形偏向 16 台、QF = 0.25・QD = -0.30、周長 41.6 m)。RF は先頭 (ST のあと)。
tecalc に渡る状態 st (st.T = 対称な 1 周の行列、st.D = 放射の微分 (減衰)、st.B = 励起 Δ) を取り出して使う。

(A) 加速電圧を変えて (3 GeV 固定、νs だけが変わる) εx、Jx、Σa を放射積分の式と比べる。(ε_x/ε_std - 1)/νs² が一定になるか。
(B) 複素固有ベクトルで直接作ったモードの量 (a_k = -Re ℓDr/λ、Δ̄_k = ℓBℓ†/2、ε_k = Δ̄_k/(1-e^{-2a_k})) が tecalc と一致するか。
(C) κ (RF の kick の傾き) を 0 から動かす: T(κ) = A K(κ)。κ → 0 で放射積分の式に戻るか。1 次の摂動 (λ'、r') の式。κ と νs² の関係。
(D) Σa と 2U0/E0。Δ_δδ の観測点の依存 (RF が先頭か末尾か)。
(E) tsymp の収束 (E → -3E²/4)、tnorm の規格化、btr の連立方程式 (全体の解と tecalc の差、モード間の相関の見積り)。
(F) 結合のあるリング (スキュー四極を 1 台): εy の kx² 比例と、1 次の摂動 (混合振幅) の式。
"""
import io
import math
import numpy as np
import _common as C

C.start(__file__)
import pysad.calc6d as c6
from pysad.ffscmd import run_sad_string

print('pysad:', C.pysad_src())
np.set_printoptions(linewidth=200, precision=5)
mc2 = 0.51099895e6
re = 2.8179403262e-15
alpha = 1 / 137.035999084
lam_c = re / alpha
Cq = 55 / (32 * math.sqrt(3)) * lam_c
Cg = 4 * math.pi / 3 * re / (mc2 / 1e9) ** 3

_cap = {}
_orig = c6.tecalc


def _wrap(st, circ, **kw):
    _cap['st'] = st
    return _orig(st, circ, **kw)


c6.tecalc = _wrap


def run(E, V, rf_first=True, skew=None):
    """E [GeV]、V [V]。戻り値 (ffs, res, st)。skew が None でなければ ST のあとにスキュー四極 (L=0.1、K1=skew、ROTATE=π/4) を 1 台"""
    cells = ' '.join(['CELL'] * 8)
    sk = ' SK' if skew is not None else ''
    ring = '(ST%s RF %s)' % (sk, cells) if rf_first else '(ST%s %s RF)' % (sk, cells)
    src = '''MOMENTUM = %gGEV;
DRIFT D1=(L=0.8) D2=(L=0.3);
BEND B=(L=1.2 ANGLE=%.15g);
QUAD QF=(L=0.3 K1=0.25) QD=(L=0.3 K1=-0.30);
QUAD SK=(L=0.1 K1=%.15g ROTATE=0.7853981633974483);
CAVI RF=(L=0 VOLT=%g HARM=64);
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=%s;
FFS;
USE RING;
''' % (E, 2 * math.pi / 16, skew or 0.0, V, ring)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    res = c6.calc6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
    return ffs, res, _cap['st']


def standard(lat, E):
    """放射積分の式: U0 [eV]、Jx、εx、I 群"""
    ri = lat.radiation_integrals()
    I1, I2, I3, I4, I5 = [ri[k] for k in ('I1', 'I2', 'I3', 'I4', 'I5')]
    g = E * 1e9 / mc2
    Jx = 1 - I4 / I2
    U0 = Cg * E**4 * I2 / (2 * math.pi) * 1e9
    return dict(U0=U0, Jx=Jx, eps=Cq * g**2 * I5 / (Jx * I2), I=(I1, I2, I3, I4, I5), g=g, Jz=2 + I4 / I2,
                sd=math.sqrt(Cq * g**2 * I3 / ((2 + I4 / I2) * I2)))


S6 = np.zeros((6, 6))
for _i in (0, 2, 4):
    S6[_i, _i + 1] = 1.0
    S6[_i + 1, _i] = -1.0


def modes(T, D, B, want=(0, 1, 2)):
    """複素固有ベクトルによるモードの量。T は (対称な) 1 周の行列、D は放射の微分、B は励起 Δ。
    シンプレクティック規格化 i r†S r = 1/2 (r = (R_X - i R_P)/2、x = R (X, P))。ℓ = r の逆行列の行 (ℓ r = 1)。
    a = -Re(ℓ D r/λ)、Δ̄ = ℓ B ℓ†/2、ε = Δ̄/(1 - e^{-2a})。モードの割り当ては成分の重み (x: 0,1、y: 2,3、z: 4,5) の最大で行う"""
    w, V = np.linalg.eig(T)
    cand = []
    for k in range(6):
        if abs(w[k].imag) < 1e-9:
            continue
        r = V[:, k]
        s = (1j * (r.conj() @ S6 @ r)).real
        if s <= 0:
            continue
        r = r * math.sqrt(0.5 / s)
        cand.append((w[k], r))
    out = {}
    names = {}
    free = list(range(len(cand)))
    for nm, idx in (('x', (0, 1)), ('z', (4, 5)), ('y', (2, 3))):
        if not free:
            break
        j = max(free, key=lambda q: abs(cand[q][1][idx[0]])**2 + abs(cand[q][1][idx[1]])**2)
        free.remove(j)
        names[nm] = cand[j]
    for nm, (lam, r) in names.items():
        ci = None
        # ℓ は 6 本の固有ベクトル (r と共役) を並べた行列の逆から
        out[nm] = dict(lam=lam, r=r)
    cols = []
    order = [n for n in ('x', 'y', 'z') if n in names]
    for nm in order:
        cols += [names[nm][1], names[nm][1].conj()]
    # 足りないモード (z が Jordan など) は ℓ を作れない
    res = {}
    if len(order) == 3:
        Rc = np.array(cols).T
        L = np.linalg.inv(Rc)
        for q, nm in enumerate(order):
            lam, r = names[nm]
            ell = L[2 * q]
            a = -(ell @ D @ r / lam).real
            db = 0.5 * (ell @ B @ ell.conj()).real
            res[nm] = dict(lam=lam, mu=np.angle(lam), a=a, dbar=db, eps=db / (1 - math.exp(-2 * a)), r=r, ell=ell)
    else:
        # x だけが必要な場合 (κ = 0): 左固有ベクトルは固有値問題の転置から
        w2, V2 = np.linalg.eig(T.T)
        for nm in order:
            lam, r = names[nm]
            k2 = int(np.argmin(abs(w2 - lam)))
            ell = V2[:, k2]
            ell = ell / (ell @ r)
            a = -(ell @ D @ r / lam).real
            db = 0.5 * (ell @ B @ ell.conj()).real
            res[nm] = dict(lam=lam, mu=np.angle(lam), a=a, dbar=db, eps=db / (1 - math.exp(-2 * a)), r=r, ell=ell)
    return res


# CUT_A
print('\n== (A) 加速電圧を変える (3 GeV 固定。νs だけが変わる) ==')
print('V[MV]   νs      εx(calc6d)/εx(式)  (比-1)/νs²   Jx(calc6d)  Jx(式)   Jy      Σa/(2U0/E0)  a_x/a_x(式)  Δ̄_x/Δ̄_x(式)')
rowsA = []
for V in (3e6, 4e6, 6e6, 12e6, 24e6):
    ffs, res, st = run(3, V)
    em = res.emit
    sd = standard(ffs.lat, 3)
    dr = -np.asarray(em['damping_rate'])
    nus = abs(em['tunes'][2])
    Jx = dr[0] * 4 / dr.sum()
    Jy = dr[1] * 4 / dr.sum()
    ratio = em['emittance_x'] / sd['eps']
    ax_std = sd['Jx'] * sd['U0'] / (2 * 3e9)
    dbar_std = sd['eps'] * 2 * ax_std
    md = modes(st.T, st.D, st.B)
    rowsA.append((nus, ratio))
    print('%5.1f  %.4f   %.5f           %.2f        %.5f     %.5f  %.5f  %.5f      %.5f      %.5f' % (
        V / 1e6, nus, ratio, (ratio - 1) / nus**2, Jx, sd['Jx'], Jy, dr.sum() / (2 * sd['U0'] / 3e9), dr[0] / ax_std,
        md['x']['dbar'] / dbar_std))
print('(比-1)/νs² がほぼ一定: 補正は νs² に比例 (νs² ∝ V cos φs)')

print('\n== (A2) E と V を同時に変える (14 章の本文の 4 点。νs と U0/E0 が両方変わる) ==')
for E, V in ((3, 3e6), (1, 0.3e6), (1, 0.1e6), (0.5, 0.05e6)):
    ffs, res, st = run(E, V)
    em = res.emit
    sd = standard(ffs.lat, E)
    dr = -np.asarray(em['damping_rate'])
    nus = abs(em['tunes'][2])
    print('E=%.1f GeV V=%.2f MV: νs=%.4f εx比=%.4f (比-1)/νs²=%.2f Jx=%.5f (式 %.5f) Jy=%.5f Σa/(2U0/E0)=%.5f U0/E0=%.2e' % (
        E, V / 1e6, nus, em['emittance_x'] / sd['eps'], (em['emittance_x'] / sd['eps'] - 1) / nus**2,
        dr[0] * 4 / dr.sum(), sd['Jx'], dr[1] * 4 / dr.sum(), dr.sum() / (2 * sd['U0'] / (E * 1e9)), sd['U0'] / (E * 1e9)))

# 基準の 1 点 (3 GeV、3 MV)
ffs, res, st = run(3, 3e6)
em = res.emit
sd = standard(ffs.lat, 3)
T, D, B = st.T, st.D, st.B
print('\n== (B) 複素固有ベクトルのモードの量 と tecalc ==')
md = modes(T, D, B)
dr = -np.asarray(em['damping_rate'])
for k, nm in enumerate(('x', 'y', 'z')):
    m = md[nm]
    print('%s: a = %.6e (tecalc %.6e)、μ/2π = %.5f (tecalc %.5f)、Δ̄ = %.5e、ε = %.5e' % (
        nm, m['a'], dr[k], abs(m['mu']) / 2 / math.pi, abs(em['tunes'][k]), m['dbar'], m['eps']))
print('εx = %.5e (tecalc %.5e、比 %.5f)、εz = %.5e (tecalc %.5e)' % (
    md['x']['eps'], em['emittance_x'], md['x']['eps'] / em['emittance_x'], md['z']['eps'], em['emittance_z']))

print('\n== (C) κ を動かす: T(κ) = A K(κ)、K = 1 + κ e_δ e_z^T ==')
K = st.rec[2][0]
kap = K[5, 4]
A = T @ np.linalg.inv(K)
Dp = D @ np.linalg.inv(K)       # D = D' K: 放射の微分は 1 周の途中の要素ごとに T_{i→end} D_i T_{0→i} (T_{0→i} は RF の K を含む) の和なので、D は κ に依存する
E4 = np.outer(np.eye(6)[5], np.eye(6)[4])
print('D K^{-1} の z 列と、D の z 列の κ 依存: D[δ,z]/D[δ,δ] = %.6e (κ = %.6e)' % (D[5, 4] / D[5, 5], kap))
print('RF の kick の傾き κ = %.6e [1/m] (δ → δ + κ z)。K の他の成分は単位行列からの差 %.1e' % (kap, np.abs(K - np.eye(6) - kap * np.outer(np.eye(6)[5], np.eye(6)[4])).max()))
Reta = A[4, 5] + A[4, :4] @ np.linalg.solve(np.eye(4) - A[:4, :4], A[:4, 5])
print('縦の 2×2 (RF を含む) のトレース = %.6f、2 - 4 sin²(μs/2) = %.6f、κ R56eff = %.6f、-4 sin²(μs/2) = %.6f (R56eff = R56 + R5j η_j = %.5f m)' % (
    np.trace(T[4:, 4:]), 2 - 4 * math.sin(abs(md['z']['mu']) / 2)**2, kap * Reta, -4 * math.sin(abs(md['z']['mu']) / 2)**2, Reta))
ax_std = sd['Jx'] * sd['U0'] / (2 * 3e9)
dbar_std = sd['eps'] * 2 * ax_std
print('標準の式: Jx = %.6f、a_x = Jx U0/(2E0) = %.6e、Δ̄_x = ε Jx U0/E0 = %.6e' % (sd['Jx'], ax_std, dbar_std))
print(' κ/κ0   a_x/a_x(式)   Δ̄_x/Δ̄_x(式)   εx/εx(式)   (εx比-1)/κ        (a_x比-1)/κ       (Δ̄比-1)/κ')
base = {}
for f in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
    Kf = np.eye(6) + f * kap * E4
    Tk = A @ Kf
    m = modes(Tk, Dp @ Kf, B, want=(0,))
    x = m['x']
    ra, rd_, re_ = x['a'] / ax_std, x['dbar'] / dbar_std, x['eps'] / sd['eps']
    if f == 0:
        base = (ra, rd_, re_)
    print(' %.2f   %.6f     %.6f      %.6f    %.4e       %.4e       %.4e' % (
        f, ra, rd_, re_, (re_ - 1) / (f * kap) if f else 0, (ra - 1) / (f * kap) if f else 0, (rd_ - 1) / (f * kap) if f else 0))
print('κ = 0 (縦運動なし): εx の比 = %.6f、a_x の比 = %.6f (Jx に直すと %.6f、式 %.6f)' % (base[2], base[0], base[0] * sd['Jx'], sd['Jx']))
Kf1 = np.eye(6) + kap * E4
m1 = modes(A @ Kf1, Dp @ Kf1, B, want=(0,))['x']
sl = (m1['eps'] / sd['eps'] - 1) / kap
print('κ = κ0 での (εx比-1)/κ = %.4f m、(2π)²/|R56eff| = %.5f 1/m → (εx比-1)/νs² ≈ %.3f (νs → 0 の極限の係数は f = 0.25 の値 %.3f m で %.3f)' % (sl, (2 * math.pi)**2 / abs(Reta), sl * (2 * math.pi)**2 / abs(Reta), 3.4587, 3.4587 * (2 * math.pi)**2 / abs(Reta)))

# 1 次の摂動: λ' = ℓ V r、r' の δ 成分 (κ z0/(λ-1))
V1 = A @ np.outer(np.eye(6)[5], np.eye(6)[4])           # dT/dκ
T0 = A
w0, V0 = np.linalg.eig(T0)
k0 = int(np.argmax([abs(V0[0, k])**2 + abs(V0[1, k])**2 if abs(w0[k].imag) > 1e-9 else 0 for k in range(6)]))
r0 = V0[:, k0]
lam0 = w0[k0]
print('RF なし (κ = 0) の x チューン (小数部) = %.5f、RF あり (calc6d) = %.5f。twiss の位相 (整数部つき) = %.5f (RF なしは 1 + %.5f の符号を揃えたもの)' % (
    abs(np.angle(lam0)) / 2 / math.pi, abs(em['tunes'][0]), res.rows[-1].mux / 2 / math.pi, 1 - abs(np.angle(lam0)) / 2 / math.pi))
wy0, Vy0 = np.linalg.eig(T0[2:4, 2:4])
print('RF なし (κ = 0) の y チューン (小数部) = %.5f、RF あり (calc6d) = %.5f。twiss の位相 = %.5f' % (
    abs(np.angle(wy0[0])) / 2 / math.pi, abs(em['tunes'][1]), res.rows[-1].muy / 2 / math.pi))
w2, V2 = np.linalg.eig(T0.T)
l0 = V2[:, int(np.argmin(abs(w2 - lam0)))]
l0 = l0 / (l0 @ r0)
lp = l0 @ V1 @ r0
# r' = -(T0 - λ)^+ (V - λ') r。ℓ r' = 0 の条件つきで解く (縮退の無い固有値なので (T0 - λ) は特異: 境界つきの系)
Mb = np.zeros((7, 7), dtype=complex)
Mb[:6, :6] = T0 - lam0 * np.eye(6)
Mb[:6, 6] = -r0
Mb[6, :6] = l0
rhs = np.concatenate([-(V1 @ r0 - lp * r0), [0]])
sol = np.linalg.solve(Mb, rhs)
r1 = sol[:6]
Rr = r0 / (r0[0])
print("\n1 次の摂動: λ' = ℓ V r = %.4e%+.4ei (x チューンの κ あたりのずれ Δμ_x = Im(λ'/λ)κ、Im(λ'/λ) = %.4e [1/m])" % (lp.real, lp.imag, (lp / lam0).imag))
z0 = r0[4]
print('z 成分の式 z0 = (A[z, x 成分] r_x)/(λ - 1) との比 = %.8f%+.8fi (z 行は z(次) = z + c·x + R56 δ で δ = 0)' % (((A[4, :4] @ r0[:4]) / (lam0 - 1) / z0).real, ((A[4, :4] @ r0[:4]) / (lam0 - 1) / z0).imag))
print('δ 成分 (κ = 0) = %.2e (|r| の最大 %.3f)' % (abs(r0[5]), np.abs(r0).max()))
print('κ=0 の x モードの固有ベクトルの z 成分 / x 成分 = %.5f%+.5fi' % ((z0 / r0[0]).real, (z0 / r0[0]).imag))
print('1 次の r\' の δ 成分 / (z0/(λ-1)) = %.6f%+.6fi (理論値 1: δ1 = κ z0/(λ-1) が主要。δ の他の入り方は無い)' % ((r1[5] / (z0 / (lam0 - 1))).real, (r1[5] / (z0 / (lam0 - 1))).imag))
# 厳密な固有ベクトル (κ の値) との比較
Tk = A + 0.01 * kap * V1 / kap  # 小さな κ で確かめる
dk = 1e-4
Tk = T0 + dk * V1
wk, Vk = np.linalg.eig(Tk)
kk = int(np.argmin(abs(wk - lam0)))
rk = Vk[:, kk]
rk = rk / (l0 @ rk)
print('厳密 (κ = 1e-4) の δ 成分 / κ = %.6f%+.6fi、1 次の r\' の δ 成分 = %.6f%+.6fi' % (
    (rk[5] / dk).real, (rk[5] / dk).imag, r1[5].real, r1[5].imag))
# a_x の κ 微分: 1 次の摂動と有限差分
Ddamp = D
def ax_of(Tm):
    kk_ = (np.linalg.solve(A, Tm)[5, 4])          # T = A K(κ) から κ を取り出す
    Kk_ = np.eye(6) + kk_ * E4
    mm = modes(Tm, Dp @ Kk_, B, want=(0,))
    return mm['x']['a']
l1 = None
Mb2 = np.zeros((7, 7), dtype=complex)
Mb2[:6, :6] = (T0 - lam0 * np.eye(6)).T
Mb2[:6, 6] = -l0
Mb2[6, :6] = r0
rhs2 = np.concatenate([-(V1.T @ l0 - lp * l0), [0]])
l1 = np.linalg.solve(Mb2, rhs2)[:6]
# 規格化 (ℓ r = 1 を保つ) に r' の ℓ 成分 = 0 の条件を選んだので、ℓ' r0 + ℓ r' = 0 を満たすように ℓ' を直す
l1 = l1 - (l1 @ r0) * l0
D0 = Dp                  # κ = 0 の D
D1 = Dp @ E4             # dD/dκ
da_pt = -((l1 @ D0 @ r0 + l0 @ D0 @ r1 + l0 @ D1 @ r0) / lam0 - (l0 @ D0 @ r0) * lp / lam0**2).real
da_fd = (ax_of(T0 + dk * V1) - ax_of(T0)) / dk
print('da_x/dκ: 1 次の摂動 = %.5e、有限差分 (κ = 1e-4) = %.5e (a_x = %.4e なので、a_x の相対変化は κ = %.4e で %.2e)' % (
    da_pt, da_fd, ax_of(T0), kap, da_pt * kap / ax_of(T0)))

print('\n== (D) Σa と 2U0/E0、Δ_δδ の観測点 ==')
tr = -0.5 * np.trace(np.linalg.solve(T, D))
print('-½ Tr(T^{-1} D) = %.6e、2U0/E0 = %.6e (比 %.6f)、tecalc の Σa = %.6e (比 %.6f)' % (
    tr, 2 * sd['U0'] / 3e9, tr / (2 * sd['U0'] / 3e9), dr.sum(), dr.sum() / (2 * sd['U0'] / 3e9)))
print('-½ Tr(D) = %.6e (比 %.6f)' % (-0.5 * np.trace(D), -0.5 * np.trace(D) / (2 * sd['U0'] / 3e9)))
print('Σ a の定義の違い: tecalc の a_k は -Re(δλ_k/λ_k)。Σ_k Re(δλ_k/λ_k) = ½ Tr(T^{-1}D) (複素固有値の和)')
U0e = sd['U0']
phis_ = math.asin(U0e / 3e6)
apc = sd['I'][0] / 41.6
bh_std = math.sqrt(3e6 / (math.pi * apc * 64 * 3e9) * (2 * math.cos(phis_) - (math.pi - 2 * phis_) * U0e / 3e6))
print('バケット高さ: 標準の式 (α_p = I1/C) = %.5e、tecalc (bh、α_p は calc6d の値) = %.5e、α_p の違いだけの補正 √(α_p 比) = %.5e' % (
    bh_std, em['bucket_height'], bh_std * math.sqrt(apc / abs(em['momentum_compaction']))))
print('同期位相 φs = asin(U0/V) = %.5f rad、α_p(calc6d)/α_p(I1/C) = %.5f、κR56eff/(-4 sin²(μs/2)) = %.5f' % (
    phis_, abs(em['momentum_compaction']) / apc, kap * Reta / (-4 * math.sin(abs(md['z']['mu']) / 2)**2)))
th = 2 * math.pi / 16
rho = 1.2 / th
g = sd['g']
Nph = 5 / (2 * math.sqrt(3)) * alpha * g * th
uc = 1.5 * lam_c * g**3 / rho * mc2 / 3e9
dee = 16 * Nph * 11 / 27 * uc**2
print('Δ_δδ: RF 先頭 (観測点は RF の前) st.B[δδ] = %.5e、16 台の N<u²> の和 = %.5e (比 %.5f)' % (B[5, 5], dee, B[5, 5] / dee))
ffs2, res2, st2 = run(3, 3e6, rf_first=False)
print('Δ_δδ: RF 末尾 (観測点は RF のあと) st.B[δδ] = %.5e (比 %.5f)、εx = %.5e (RF 先頭 %.5e)、σδ = %.5e (%.5e)' % (
    st2.B[5, 5], st2.B[5, 5] / dee, res2.emit['emittance_x'], em['emittance_x'], res2.emit['momentum_spread'], em['momentum_spread']))
Bz = st.B
print('RF 末尾の B = K B_0 K^T: δδ = B0_δδ + 2κ B0_zδ + κ² B0_zz = %.5e (st.B の値 %.5e)' % (
    Bz[5, 5] + 2 * kap * Bz[4, 5] + kap**2 * Bz[4, 4], st2.B[5, 5]))
for E in (1, 2, 3):
    ffs3, res3, st3 = run(E, 3e6)
    s3 = standard(ffs3.lat, E)
    g3 = s3['g']
    uc3 = 1.5 * lam_c * g3**3 / rho * mc2 / (E * 1e9)
    d3 = 16 * 5 / (2 * math.sqrt(3)) * alpha * g3 * th * 11 / 27 * uc3**2
    print('E = %d GeV: B_δδ/(16 台の和) - 1 = %.3e、U0/E0 = %.3e' % (E, st3.B[5, 5] / d3 - 1, s3['U0'] / (E * 1e9)))

print('\n== (E1) tsymp: t → t(1 - E/2)、E = t^s t - 1。E → -3E²/4 (sympy は別) ==')
rng = np.random.default_rng(1)
Rn = res.emit['normal_r']
for eps_ in (1e-1, 1e-2, 1e-3):
    t = Rn + eps_ * rng.standard_normal((6, 6))
    def tinv6(a):
        return -S6 @ a.T @ S6
    E0 = tinv6(t) @ t - np.eye(6)
    t2 = t @ (1.5 * np.eye(6) - 0.5 * tinv6(t) @ t)
    E1 = tinv6(t2) @ t2 - np.eye(6)
    print('|E| = %.3e → |E\'| = %.3e (予想 3/4 |E²| ≈ %.3e)、行列の 2 乗 -3/4 E²との差 %.1e' % (
        np.abs(E0).max(), np.abs(E1).max(), 0.75 * np.abs(E0 @ E0).max(), np.abs(E1 + 0.75 * E0 @ E0).max()))
Rr_ = Rn
print('tecalc の R: R^T S R - S の最大 = %.2e、R^{-1} T R の x ブロック (回転): [[%.5f, %.5f],[%.5f, %.5f]] (cos μx = %.5f)' % (
    np.abs(Rr_.T @ S6 @ Rr_ - S6).max(), *(np.linalg.inv(Rr_) @ T @ Rr_)[:2, :2].ravel(), math.cos(2 * math.pi * abs(em['tunes'][0]))))
print('R の x モードの列 (R_X, R_P) の x・px・y・py 成分 (tnorm: r(i,i+1) の x 成分 = 0): R[0,1] = %.2e' % Rr_[0, 1])

print('\n== (E2) btr: 全体の連立方程式の解 と tecalc ==')
Rn_i = np.linalg.inv(Rn)
Lp = Rn_i @ (T + D) @ Rn            # Λ' = P + rd
rdm = Rn_i @ D @ Rn
Bn = Rn_i @ B @ Rn_i.T
# 全体の解: Σ = Λ' Σ Λ'^T + Bn → (1 - Λ'⊗Λ') vec Σ = vec Bn
Mfull = np.eye(36) - np.kron(Lp, Lp)
Sig_full = np.linalg.solve(Mfull, Bn.reshape(36)).reshape(6, 6)
Sn_tec = Rn_i @ em['beam_matrix'] @ Rn_i.T
print('εx: 全体の解 %.5e、tecalc %.5e (比 %.5f)、εz: 全体の解 %.5e、tecalc %.5e (比 %.5f)' % (
    0.5 * (Sig_full[0, 0] + Sig_full[1, 1]), 0.5 * (Sn_tec[0, 0] + Sn_tec[1, 1]), (Sig_full[0, 0] + Sig_full[1, 1]) / (Sn_tec[0, 0] + Sn_tec[1, 1]),
    0.5 * (Sig_full[4, 4] + Sig_full[5, 5]), 0.5 * (Sn_tec[4, 4] + Sn_tec[5, 5]), (Sig_full[4, 4] + Sig_full[5, 5]) / (Sn_tec[4, 4] + Sn_tec[5, 5])))
print('btr の対角ブロックの置き換え: 2a_k - a_k² = %s (tecalc の式)、1 - e^{-2a} = %s' % (
    ['%.6e' % (2 * x_ - x_**2) for x_ in dr], ['%.6e' % (1 - math.exp(-2 * x_)) for x_ in dr]))
# モード間の相関の見積り: Σ_kl = B_kl / (1 - λ_k λ_l) を回転系で
mu = [2 * math.pi * abs(t_) for t_ in em['tunes']]
Bn_c = Bn
def zeta(blk):
    return None
print('モード間の相関: 正規座標で複素 ζ_k = X_k + iP_k (ブロックの回転が e^{-iμ_k})。B_kl^{(±)} = (⟨ζ_k ζ_l^*⟩, ⟨ζ_k ζ_l⟩ の励起)')
Tc = np.zeros((3, 6), dtype=complex)
for k in range(3):
    Tc[k, 2 * k] = 1
    Tc[k, 2 * k + 1] = 1j
Bz_ = Tc @ Bn @ Tc.conj().T        # ⟨ζ_k ζ_l*⟩ の励起
Bzz = Tc @ Bn @ Tc.T               # ⟨ζ_k ζ_l⟩ の励起
mu_s = [np.angle(np.linalg.eigvals(Lp[2 * k:2 * k + 2, 2 * k:2 * k + 2])[0]) for k in range(3)]
# Λ' のブロックは e^{-a}P(μ): ζ → e^{-a}e^{∓iμ}ζ。符号は Λ' から直接取る
lam_blk = []
for k in range(3):
    blk = Lp[2 * k:2 * k + 2, 2 * k:2 * k + 2]
    lam_blk.append((blk[0, 0] + 1j * blk[0, 1] * 0 + 0))   # 後で厳密に
# ζ の 1 周の係数: ζ' = (X' + iP') = Σ blk を使い、ζ の固有値を求める
zl = []
for k in range(3):
    blk = Lp[2 * k:2 * k + 2, 2 * k:2 * k + 2]
    ev = np.linalg.eigvals(blk)
    # ζ = X + iP に対する固有値: (1, i) 左固有ベクトル
    row = np.array([1, 1j]) @ blk
    zl.append(row[0] / 1.0)       # ζ' = zl ζ (+ 共役の項は 0 なら)
pairs = [(0, 2, 'x-z'), (0, 1, 'x-y'), (1, 2, 'y-z')]
Smax = max(np.abs(Sn_tec[0:2, 4:]).max(), np.abs(Sn_tec[0:2, 2:4]).max(), np.abs(Sn_tec[2:4, 4:]).max())
print('tecalc の Σ_n のモード間 (x-z) の相関の最大 = %.3e' % np.abs(Sn_tec[0:2, 4:]).max())
k, l = 0, 2
lk, ll = zl[k], zl[l]
c1 = Bz_[k, l] / (1 - lk * ll.conjugate())
c2 = Bzz[k, l] / (1 - lk * ll)
print('見積り (x-z): ⟨ζ_x ζ_z*⟩ = B/(1 - λ_x λ_z*) = %.3e、⟨ζ_x ζ_z⟩ = B/(1 - λ_x λ_z) = %.3e。|1 - λ_x λ_z*| = %.4f、|1 - λ_x λ_z| = %.4f' % (
    abs(c1), abs(c2), abs(1 - lk * ll.conjugate()), abs(1 - lk * ll)))
Sn_cc = Tc @ Sn_tec @ Tc.conj().T
Sn_cn = Tc @ Sn_tec @ Tc.T
print('tecalc の ⟨ζ_x ζ_z*⟩ = %.3e、⟨ζ_x ζ_z⟩ = %.3e' % (abs(Sn_cc[0, 2]), abs(Sn_cn[0, 2])))
# モード間の相関の主な源: 放射の行列 rd の非対角ブロック (モード間の混合)。それを 0 にして解き直す
Lp0 = Lp.copy()
for kk_ in (0, 2, 4):
    for ll_ in (0, 2, 4):
        if kk_ != ll_:
            Lp0[kk_:kk_ + 2, ll_:ll_ + 2] = 0.0
Sig0 = np.linalg.solve(np.eye(36) - np.kron(Lp0, Lp0), Bn.reshape(36)).reshape(6, 6)
print('rd のモード間ブロックを 0 にして解くと、x-z の相関の最大 = %.3e (全体の解 %.3e)、x の ε = %.5e (全体 %.5e)' % (
    np.abs(Sig0[0:2, 4:]).max(), np.abs(Sig_full[0:2, 4:]).max(), 0.5 * (Sig0[0, 0] + Sig0[1, 1]), 0.5 * (Sig_full[0, 0] + Sig_full[1, 1])))
print('rd のモード間ブロックの大きさ: |rd_xz| の最大 = %.3e、|rd_zx| の最大 = %.3e (対角ブロックの減衰 a_x = %.3e)' % (
    np.abs(rdm[0:2, 4:]).max(), np.abs(rdm[4:, 0:2]).max(), dr[0]))

print('\n== (F) 結合のあるリング: スキュー四極 (L = 0.1 m、ROTATE = π/4) を 1 台 ==')
print('基準は K1 = 0 の同じリング (0.1 m のドリフトが入る)。T(K1) - T(0) を摂動 V とする')
ffs_b, res_b, st_b = run(3, 3e6, skew=0.0)
mp = modes(st_b.T, st_b.D, st_b.B)
print('基準 (K1 = 0): εy = %.2e、εx = %.5e' % (res_b.emit['emittance_y'], res_b.emit['emittance_x']))
print('K1      εx(calc6d)   εy(calc6d)    εy/εx      1 次の摂動 (εy/εx)    PT/厳密   Jx      Jy      Jz')
rowsF = []
for ks in (0.002, 0.004, 0.008):
    ffs4, res4, st4 = run(3, 3e6, skew=ks)
    e4 = res4.emit
    d4 = -np.asarray(e4['damping_rate'])
    Vp = st4.T - st_b.T
    # y モードの左固有ベクトルの 1 次の補正: ℓ_y' = ℓ_y + Σ_m ℓ_m (ℓ_y V r_m)/(λ_y - λ_m)  (m は y 自身を除く 5 本)
    ly, lyl = mp['y']['ell'], mp['y']['lam']
    allm = []
    for nm in ('x', 'z', 'y'):
        allm.append((mp[nm]['lam'], mp[nm]['r'], mp[nm]['ell']))
        allm.append((mp[nm]['lam'].conjugate(), mp[nm]['r'].conj(), mp[nm]['ell'].conj()))
    ly1 = np.zeros(6, dtype=complex)
    for lm, rm, ellm in allm:
        if abs(lm - lyl) < 1e-12:
            continue
        ly1 += ellm * (ly @ Vp @ rm) / (lyl - lm)
    lyp = ly + ly1
    dby = 0.5 * (lyp @ st4.B @ lyp.conj()).real
    epy = dby / (1 - math.exp(-2 * mp['y']['a']))
    rowsF.append((ks, e4['emittance_y'] / e4['emittance_x'], epy / e4['emittance_x']))
    print('%.3f   %.4e   %.4e   %.4e   %.4e        %.4f    %.5f %.5f %.5f' % (
        ks, e4['emittance_x'], e4['emittance_y'], e4['emittance_y'] / e4['emittance_x'], epy / e4['emittance_x'],
        epy / e4['emittance_y'], d4[0] * 4 / d4.sum(), d4[1] * 4 / d4.sum(), d4[2] * 4 / d4.sum()))
    contrib = []
    for lm, rm, ellm in allm:
        if abs(lm - lyl) < 1e-12:
            continue
        c = (ly @ Vp @ rm) / (lyl - lm)
        contrib.append(abs(c)**2 * 0.5 * (ellm @ st4.B @ ellm.conj()).real)
    print('        混合の内訳 Σ|c|²Δ̄_m: x %.2e、x* %.2e、z %.2e、z* %.2e、y* %.2e (y 自身の励起 Δ̄_y は 0)' % tuple(contrib))
print('εy/εx の比 (K1 = 0.004 / 0.002、0.008 / 0.004): %.3f、%.3f (弱い結合: K1² に比例なら 4)' % (rowsF[1][1] / rowsF[0][1], rowsF[2][1] / rowsF[1][1]))
