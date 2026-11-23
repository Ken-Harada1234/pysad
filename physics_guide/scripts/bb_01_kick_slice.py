# -*- coding: utf-8 -*-
"""BEAMBEAM (16_undulator_beambeam.md): Bassetti-Erskine のキック、スライスの写像、ビームストラールング。

(a) 1 スライス (SLICE = 1、z = 0 の粒子、px = py = 0) のキック Δpx、Δpy、Δδ を、ガウス分布の 2 次元ポテンシャル
    U = c ∫0^∞ [1 - exp(-x²/(2σx²+t) - y²/(2σy²+t))] / √((2σx²+t)(2σy²+t)) dt の数値積分 (scipy quad) の微分と比べる
    (丸い・σx > σy・σx < σy の 3 通り。Δδ は σ² の変化率の項: -gxy = -(σ12x ∂U/∂σx² + σ12y ∂U/∂σy²))。
    小振幅の線形の傾き 2 c/(σx(σx+σy))。
(b) 6 次元の写像の閉じ方: SLICE = 5、交差角、軌道のずれのある BEAMBEAM の写像の 6×6 ヤコビアンのシンプレクティック性
    (標準の J、(z, δ) が共役)。交差角だけ (NP = 0) の往復が恒等。
(c) ビームストラールング: 定数 cnbs = anrad γ、cubs = cuc γ²、光子エネルギー表の平均と 2 次モーメント (cave = 8/(15√3)、cuu = 11/27)。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
from pysad import beambeam as BB
from scipy.integrate import quad

mom = 1.7e9
me = 0.51099895e6
re = 2.8179403262e-15
gam = math.sqrt(1 + (mom / me) ** 2)


def run(params, rows, nrep=1):
    """rows: 各粒子の (x, px_canon, y, py_canon, z, δ)。BEAMBEAM 1 個を通す。戻り値: 同じ形 (canon)"""
    bb = BB.setup(params, mom)
    S = np.zeros((len(rows), 16))
    for i, r in enumerate(rows):
        pn = 1 + r[5]
        S[i, :6] = [r[0], r[1] / pn, r[2], r[3] / pn, r[4], r[5]]
    BB.apply(S, len(rows), bb)
    out = np.zeros((len(rows), 6))
    for i in range(len(rows)):
        pn = 1 + S[i, 5]
        out[i] = [S[i, 0], S[i, 1] * pn, S[i, 2], S[i, 3] * pn, S[i, 4], S[i, 5]]
    return out, bb


def U_pot(x, y, sx2, sy2):
    f = lambda t: (1 - math.exp(-x * x / (2 * sx2 + t) - y * y / (2 * sy2 + t))) / math.sqrt((2 * sx2 + t) * (2 * sy2 + t))
    # ∫0^∞: 区間を分けて積分
    v = 0.0
    a = 0.0
    for b in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, np.inf):
        v += quad(f, a, b, epsabs=0, epsrel=1e-12, limit=200)[0]
        a = b
    return v + 2 * math.log(math.sqrt(sx2) + math.sqrt(sy2))     # 物理のポテンシャル: 遠方で ln r² (σ に依らない定数) になるよう σ 依存の定数を足す


print('== (a) Bassetti-Erskine のキックと 2 次元ポテンシャル ==')
NP = 1e10
ce = NP * re / gam
print('γ = %.2f、NP = %.0e、ce = NP r_e/γ = %.5e m' % (gam, NP, ce))
cases = [('丸い (σ = 0.1 mm)', dict(BX=1.0, BY=1.0, EMIX=1e-8, EMIY=1e-8, AX=0.3, AY=0.3)),
         ('σx > σy (0.2 mm、0.05 mm)', dict(BX=4.0, BY=0.25, EMIX=1e-8, EMIY=1e-8, AX=0.5, AY=-0.4)),
         ('σx < σy (0.05 mm、0.2 mm)', dict(BX=0.25, BY=4.0, EMIX=1e-8, EMIY=1e-8, AX=-0.2, AY=0.6))]
pts = [(1.5e-4, 0.7e-4), (-0.6e-4, 2.5e-4), (3.0e-4, -1.0e-4), (0.2e-4, -0.1e-4)]
for name, kw in cases:
    params = dict(DP=1e-3, SIGZ=0.01, SLICE=1, NP=NP)
    params.update(kw)
    rows = [(px_, 0.0, py_, 0.0, 0.0, 0.0) for (px_, py_) in pts]
    out, bb = run(params, rows)
    b5 = bb['benv5']
    sx2, sy2 = kw['BX'] * kw['EMIX'], kw['BY'] * kw['EMIY']          # 強ビームの IP の 2 次モーメント (η = 0)
    s12x, s12y = -kw['AX'] * kw['EMIX'], -kw['AY'] * kw['EMIY']      # σ12 = -α ε (z = 0 のスライス)
    print('  (setup の傾き角 = %.4f rad。benv5 の σx = %.3e、σy = %.3e)' % (bb['xangle'][4], math.sqrt(b5[0, 0]), math.sqrt(b5[2, 2])))
    print('\n%s: σx = %.4e、σy = %.4e m、σ12x = %.3e、σ12y = %.3e (z = 0 のスライス)' % (name, math.sqrt(sx2), math.sqrt(sy2), s12x, s12y))
    for (xx, yy), o in zip(pts, out):
        h = 1e-4 * min(abs(xx), abs(yy)) + 1e-9
        dUdx = (U_pot(xx + h, yy, sx2, sy2) - U_pot(xx - h, yy, sx2, sy2)) / (2 * h)
        dUdy = (U_pot(xx, yy + h, sx2, sy2) - U_pot(xx, yy - h, sx2, sy2)) / (2 * h)
        hs = 1e-4 * sx2
        dUdsx = (U_pot(xx, yy, sx2 + hs, sy2) - U_pot(xx, yy, sx2 - hs, sy2)) / (2 * hs)
        hs = 1e-4 * sy2
        dUdsy = (U_pot(xx, yy, sx2, sy2 + hs) - U_pot(xx, yy, sx2, sy2 - hs)) / (2 * hs)
        ex = -ce * dUdx
        ey = -ce * dUdy
        eg = -ce * (s12x * dUdsx + s12y * dUdsy) + 0.25 * (ex * ex + ey * ey)     # 最後の項: exp(:D:) の戻しの px²/4 (kick の後の px)
        print('  (x, y) = (%+.1e, %+.1e): Δpx = %+.6e (積分 %+.6e、差 %.1e)、Δpy = %+.6e (積分 %+.6e、差 %.1e)、Δδ = %+.4e (積分 %+.4e、差 %.1e)' % (
            xx, yy, o[1], ex, o[1] - ex, o[3], ey, o[3] - ey, o[5], eg, o[5] - eg))
    # 線形の傾き
    o, _ = run(params, [(1e-9, 0.0, 0.0, 0.0, 0.0, 0.0)])
    sx, sy = math.sqrt(sx2), math.sqrt(sy2)
    print('  小振幅: Δpx/x = %.6e、標準の -2 ce/(σx(σx+σy)) = %.6e' % (o[0][1] / 1e-9, -2 * ce / (sx * (sx + sy))))


print('\n-- k64 モードの傾きの回転 (B34): σx < σy の強ビームは benv_norm_axis が 90 度回すので傾きが入る --')
old = C.set_mode('k64')
kw = cases[2][1]
params = dict(DP=1e-3, SIGZ=0.01, SLICE=1, NP=NP)
params.update(kw)
rows = [(px_, 0.0, py_, 0.0, 0.0, 0.0) for (px_, py_) in pts]
out_k, _ = run(params, rows)
C.set_mode(old)
out_v, _ = run(params, rows)
for (xx, yy), ok, ov in zip(pts, out_k, out_v):
    print('  (x, y) = (%+.1e, %+.1e): v2 Δ(px, py) = (%+.4e, %+.4e)、k64 = (%+.4e, %+.4e)' % (xx, yy, ov[1], ov[3], ok[1], ok[3]))

# ---------------------------------------------------------------- (b)
print('\n== (b) 6 次元写像のシンプレクティック性 ==')
params = dict(BX=0.5, BY=0.01, EMIX=1e-8, EMIY=1e-10, SIGZ=0.005, DP=1e-3, SLICE=5, NP=1e10, XANGLE=0.01,
              AX=0.1, AY=-0.05, DX=2e-5, DY=-1e-5, DPX=3e-7, DPY=1e-7)
z0 = np.array([8e-5, 2e-5, 2e-6, -3e-6, 0.004, 5e-4])
scale = np.array([1e-4, 1e-4, 1e-5, 1e-5, 1e-2, 1e-3])


def Fmap(v, prm):
    return run(prm, [tuple(v)])[0][0]


def jac(prm, v0):
    J = np.zeros((6, 6))
    for k in range(6):
        h = np.zeros(6)
        h[k] = 1e-5 * scale[k]
        J[:, k] = (Fmap(v0 + h, prm) - Fmap(v0 - h, prm)) / (2 * h[k])
    return J


Jst = np.zeros((6, 6))
for k in range(3):
    Jst[2 * k, 2 * k + 1] = 1.0
    Jst[2 * k + 1, 2 * k] = -1.0
for label, prm in (('SLICE = 5、交差角 10 mrad、軌道のずれ、AX・AY ≠ 0', params),
                   ('同じ、傾いた強ビーム (R1..R4 で x-y 結合)', dict(params, R1=0.3, R2=0.0, R3=0.0, R4=0.2)),
                   ('SLICE = 1、交差角なし', dict(params, SLICE=1, XANGLE=0.0, DX=0.0, DY=0.0, DPX=0.0, DPY=0.0))):
    Jm = jac(prm, z0)
    dev = Jm.T @ Jst @ Jm - Jst
    # 無次元化: 各成分を scale で割って (変数のスケールに依存しない) 偏差の最大値
    sc = np.diag(scale)
    dev_n = (sc @ dev @ sc) / (np.abs(sc @ Jst @ sc).max())
    print('%s: max|J^T S J - S| (規格化) = %.2e、ヤコビアンの行列式 = %.12f' % (label, np.abs(dev_n).max(), np.linalg.det(Jm)))

prm0 = dict(params, NP=0.0)
o = [Fmap(z0, prm0)]
print('NP = 0 (キックなし。交差角・軌道のずれ・スライスの drift の往復): |出力 - 入力| の最大 (成分ごと) =', np.abs(o[0] - z0))
# 交差角の線形の部分
v1 = np.array([0.0, 0.0, 0.0, 0.0, 0.01, 0.0])
o1 = Fmap(v1, dict(BX=0.5, BY=0.01, EMIX=1e-8, EMIY=1e-10, SIGZ=0.005, DP=1e-3, SLICE=1, NP=0.0, XANGLE=0.01))
print('(c1) 交差角 φ = 0.01 のみ・キックなし: z = 0.01 の粒子の往復の後 = ', o1, '(恒等)')

# ---------------------------------------------------------------- (c)
print('\n== (c) ビームストラールング ==')
alpha = 1 / 137.035999084
lamC = re / alpha
anrad = 5 * alpha / (2 * math.sqrt(3))
cuc = 1.5 * lamC
bbs = BB.setup(dict(BX=1.0, BY=1.0, EMIX=1e-8, EMIY=1e-8, DP=1e-3, SIGZ=0.01, SLICE=1, NP=1e10, BSTRL=1), mom)
print('cnbs = %.6e、anrad γ = %.6e (比 %.9f)  [N_γ = cnbs/ρ_i dsz = anrad γ θ]' % (bbs['cnbs'], anrad * gam, bbs['cnbs'] / (anrad * gam)))
print('cubs = %.6e、cuc γ² = %.6e (比 %.9f)  [u_c/E = cubs/ρ]' % (bbs['cubs'], cuc * gam**2, bbs['cubs'] / (cuc * gam**2)))
xit, swit = BB.besseltab()
print('光子エネルギーの表: 点数 %d、SwI の範囲 %.3f〜%.3f' % (len(xit), swit.min(), swit.max()))
# xx ∈ (0,1) の一様乱数 → ln(u/u_c) = xi を表の線形補間
uu = (np.arange(200000) + 0.5) / 200000
xi_s = np.interp(uu, swit[1:1001], xit[1:1001])      # SwI は x が大きい側から積算した光子数の分布 (1 始まりの表)
e1 = np.mean(np.exp(xi_s))
e2 = np.mean(np.exp(2 * xi_s))
print('<u>/u_c = %.6f (cave = 8/(15√3) = %.6f)、<u²>/u_c² = %.6f (cuu = 11/27 = %.6f)' % (e1, 8 / (15 * math.sqrt(3)), e2, 11 / 27))

# 表の精度: 同じ格子 (x_i = 1.1^40 / 1.02^i) で、S(x) = ∫_x^∞ K_{5/3} と累積 ∫_x^∞ S dx を scipy で厳密に作り直して同じ手順でモーメントを取る
from scipy.special import kv
from scipy.integrate import cumulative_trapezoid
xg = 1.1**40 / 1.02**np.arange(1, 1001)
xf = np.exp(np.linspace(np.log(1.1**40), np.log(xg[-1]), 400001))
Kf = kv(5 / 3, xf)
Sf = -cumulative_trapezoid(Kf[::-1][::-1], xf, initial=0.0)          # xf は減少: ∫_{xf}^{x_top} K を負の向きで
Sf = -Sf                                                                # S(x) = ∫_x^{x_top} K_{5/3} (x_top = 45 以上は無視)
Cf = -cumulative_trapezoid(Sf, xf, initial=0.0)
Cf = -Cf
lxf = np.log(xf)
Cn = Cf / Cf[-1]
xi_ex = np.interp(uu, Cn, lxf)
print('厳密な累積の逆引き (同じ手順): <u>/u_c = %.6f、<u²>/u_c² = %.6f (理論 %.6f、%.6f)' % (
    np.mean(np.exp(xi_ex)), np.mean(np.exp(2 * xi_ex)), 8 / (15 * math.sqrt(3)), 11 / 27))
