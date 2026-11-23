# -*- coding: utf-8 -*-
"""行列から光学へ (12_linear_optics.md)。8 セルの FODO リング (3 GeV、扇形偏向 16 台) で確かめる。

(a) 閉軌道の Newton 法: 1 周の追跡 f(c) の不動点。収束の速さ (残差の列)。1 回目の更新が線形解 (1 - M)^{-1} と一致するか。
(b) 分散: eta = (1 - M4)^{-1} m6 と、δ を変えた閉軌道の差分 (追跡) と、pysad の Twiss の ex, epx。
(c) Courant-Snyder の行列 A と (α, β, μ): B A B^{-1} = 回転 (sympy)。伝搬 (ドリフト) の式。
(d) x-y 結合の分解 (qmdiag): 結合した安定なシンプレクティック行列で、K64 の式 (qmdiag.f の通常の場合) が R M R^{-1} をブロック対角にするか。
    c^2 + det r = 1、R のシンプレクティック性、ヘルプの R^{-1} (r の符号を反転) の確認。pysad.optics.periodic_r との一致。
(e) 色収差 (偏向あり): 1 周の追跡の Jacobian から作ったチューンの δ 微分、lat.chromaticity()、SAD の CHROMATICITY (CELL; CALC) の比較。
"""
import io
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import rdttrack, optics

print('pysad:', C.pysad_src())
S2 = np.array([[0.0, 1.0], [-1.0, 0.0]])


def ring(extra_el='', extra_in_cell=''):
    src = '''MOMENTUM = 3GEV;
DRIFT D1=(L=0.8) D2=(L=0.3);
BEND B=(L=1.2 ANGLE=%.15g);
QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);
%s
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL %s);
FFS;
USE RING;
''' % (2 * math.pi / 16, extra_el, extra_in_cell)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


ffs = ring()
lat = ffs.lat
elements = ffs.lat.elements
mom = ffs.momentum()
M6 = lat.transfer_matrix()
M4 = M6[:4, :4]
tw = lat.twiss()
tw0 = tw[0]
print('リング: 8 セル、周長 %.3f m、チューン (νx, νy) = (%.4f, %.4f)' % (sum(e.L for e in elements), tw0['nx'] / (2 * math.pi), tw0['ny'] / (2 * math.pi)))

print('\n== (a) 閉軌道の Newton 法 (非線形: 6 極を 1 台、ずれた四極を 1 台) ==')
ffs2 = ring('QUAD KK=(L=0.2 K1=0.3 DX=3e-4);\nSEXT SX=(L=0.1 K2=30);', 'KK SX')
f = rdttrack.one_turn(ffs2.lat.elements, ffs2.momentum(), 0.0)
z = np.zeros(4)
res = []
for it in range(8):
    r = f(z)[0] - z
    res.append(np.abs(r).max())
    J = rdttrack.jacobian(f, z)
    z = z - np.linalg.solve(J - np.eye(4), r)
print('反復ごとの残差 max|f(c) - c|:', ' '.join('%.1e' % v for v in res))
print('閉軌道 (x, px, y, py) =', np.array2string(z, precision=6))
# 1 回目の更新: 非線形性を無視した線形解
f0 = rdttrack.one_turn(ffs2.lat.elements, ffs2.momentum(), 0.0)
z1 = np.zeros(4)
r0 = f0(z1)[0] - z1
J0 = rdttrack.jacobian(f0, z1)
dz = -np.linalg.solve(J0 - np.eye(4), r0)
print('1 回目の更新 (原点で線形化) の x = %.6e、収束した閉軌道の x = %.6e' % (dz[0], z[0]))

print('\n== (b) 分散 ==')
d4 = M6[:4, 5]
eta = np.linalg.solve(np.eye(4) - M4, d4)
print('eta = (1 - M4)^{-1} m6 : (ηx, ηx\', ηy, ηy\') =', np.array2string(eta, precision=6))
fa = rdttrack.one_turn(elements, mom, 0.0)
cods = []
for dd in (1e-5, -1e-5):
    fd = rdttrack.one_turn(elements, mom, dd)
    zz = np.zeros(4)
    for _ in range(8):
        r = fd(zz)[0] - zz
        zz = zz - np.linalg.solve(rdttrack.jacobian(fd, zz) - np.eye(4), r)
    cods.append(zz)
print('δ = ±1e-5 の閉軌道の差分 (追跡 + Newton):', np.array2string((cods[0] - cods[1]) / 2e-5, precision=6))
print('pysad の Twiss: (ex, epx) = (%.6f, %.6f)' % (tw0['ex'], tw0['epx']))

print('\n== (c) Courant-Snyder の行列と Twiss の伝搬 ==')
al, be, mu = sp.symbols('alpha beta mu', real=True)
ga = (1 + al**2) / be
A = sp.Matrix([[sp.cos(mu) + al * sp.sin(mu), be * sp.sin(mu)], [-ga * sp.sin(mu), sp.cos(mu) - al * sp.sin(mu)]])
Bm = sp.Matrix([[1 / sp.sqrt(be), 0], [al / sp.sqrt(be), sp.sqrt(be)]])
Rm = sp.Matrix([[sp.cos(mu), sp.sin(mu)], [-sp.sin(mu), sp.cos(mu)]])
print('B A B^{-1} - P(μ) = 0 か:', sp.simplify(Bm * A * Bm.inv() - Rm) == sp.zeros(2, 2))
print('det B = %s、B のシンプレクティック性 B^T J B - J =' % sp.simplify(Bm.det()), sp.simplify(Bm.T * sp.Matrix([[0, 1], [-1, 0]]) * Bm - sp.Matrix([[0, 1], [-1, 0]])))
# 一般の 2x2 T による伝搬: A2 = T A1 T^{-1} から (α2, β2)
t11, t12, t21, t22 = sp.symbols('t11 t12 t21 t22')
T = sp.Matrix([[t11, t12], [t21, t22]])
A2 = T * A * T.inv()
T_sym = sp.simplify(T.det())
b2 = sp.simplify(A2[0, 1] / sp.sin(mu))
b2 = sp.simplify(b2.subs(t22, (1 + t12 * t21) / t11))
print('β2 = (T A T^{-1})_{12} / sin μ =', sp.factor(b2), ' (det T = 1 を使用)')
print('  -> β2 = t11² β - 2 t11 t12 α + t12² γ (γ = (1 + α²)/β) との差:', sp.simplify(b2 - (t11**2 * be - 2 * t11 * t12 * al + t12**2 * ga)))
# pysad の行とドリフト
rows = tw
cnt = 0
worst = 0.0
for k, el in enumerate(elements):
    if type(el).__name__ == 'Drift' and k + 1 < len(rows) and k < 12:
        r1, r2 = rows[k], rows[k + 1]
        L = el.L
        g1 = (1 + r1['ax']**2) / r1['bx']
        b2p = r1['bx'] - 2 * r1['ax'] * L + g1 * L**2
        a2p = r1['ax'] - g1 * L
        worst = max(worst, abs(b2p - r2['bx']), abs(a2p - r2['ax']))
        cnt += 1
print('ドリフト %d 台: β2 = β - 2αL + γL²、α2 = α - γL と pysad の Twiss の差の最大 = %.2e' % (cnt, worst))
# 位相進み
r1, r2 = rows[2], rows[3]
print('(参考) 行 2、3 の nx の差 (SAD の Twiss の位相進み。2π で割る前): %.6f' % (r2['nx'] - r1['nx']))

print('\n== (d) x-y 結合の分解 (qmdiag) ==')
rng = np.random.default_rng(3)


def rand_symp4():
    """結合した安定なシンプレクティック 4x4 を作る: 安定な (x, y) 回転に、ランダムなシンプレクティック V を掛ける"""
    def rot(mu_, a_, b_):
        g_ = (1 + a_**2) / b_
        return np.array([[math.cos(mu_) + a_ * math.sin(mu_), b_ * math.sin(mu_)], [-g_ * math.sin(mu_), math.cos(mu_) - a_ * math.sin(mu_)]])
    D = np.zeros((4, 4))
    D[:2, :2] = rot(2 * math.pi * 0.31, 0.4, 5.0)
    D[2:, 2:] = rot(2 * math.pi * 0.17, -0.8, 3.0)
    H = rng.normal(0, 0.3, (4, 4))
    Hm = np.zeros((4, 4))
    J4 = np.zeros((4, 4))
    J4[:2, :2] = S2
    J4[2:, 2:] = S2
    Hs = (H + H.T) / 2
    from scipy.linalg import expm
    V = expm(J4 @ Hs * 0.5)
    return V @ D @ np.linalg.inv(V), J4


M, J4 = rand_symp4()
print('M のシンプレクティック性 |M^T J M - J| = %.1e' % np.abs(M.T @ J4 @ M - J4).max())
tm = lambda i, j: M[i - 1, j - 1]
trpq = tm(1, 1) + tm(2, 2) - tm(3, 3) - tm(4, 4)
trab2 = 4 * ((tm(3, 1) + tm(2, 4)) * (tm(4, 2) + tm(1, 3)) + (tm(1, 4) - tm(3, 2)) * (tm(4, 1) - tm(2, 3))) + trpq**2
c = math.sqrt(0.5 * (1 + abs(trpq) / math.sqrt(trab2)))
aa = -1 / c / math.sqrt(trab2) if trpq >= 0 else 1 / c / math.sqrt(trab2)
r = np.array([[aa * (tm(3, 1) + tm(2, 4)), aa * (tm(3, 2) - tm(1, 4))], [aa * (tm(4, 1) - tm(2, 3)), aa * (tm(4, 2) + tm(1, 3))]])
Rr = np.block([[c * np.eye(2), S2 @ r.T @ S2], [r, c * np.eye(2)]])
Rinv = np.block([[c * np.eye(2), -S2 @ r.T @ S2], [-r, c * np.eye(2)]])
Md = Rr @ M @ Rinv
print('c² + det r - 1 = %.1e (c = %.6f、det r = %.6f)' % (c**2 + np.linalg.det(r) - 1, c, np.linalg.det(r)))
print('R のシンプレクティック性 |R^T J R - J| = %.1e、ヘルプの R^{-1} (r の符号を逆) と R の積 |R R^{-1} - 1| = %.1e' % (np.abs(Rr.T @ J4 @ Rr - J4).max(), np.abs(Rr @ Rinv - np.eye(4)).max()))
print('R M R^{-1} の非対角ブロック (xy、yx) の最大 = %.1e、対角ブロックの最大 = %.2f' % (max(np.abs(Md[:2, 2:]).max(), np.abs(Md[2:, :2]).max()), max(np.abs(Md[:2, :2]).max(), np.abs(Md[2:, 2:]).max())))
Rp, Ap, Bp = optics.periodic_r(M)
print('pysad.optics.periodic_r の R と、上の r の差 = %.1e' % np.abs(Rp - r).max())
print('A = (R M R^{-1}) の左上ブロックと、pysad の A の差 = %.1e' % np.abs(Md[:2, :2] - Ap).max())
# (L) 式: Q + r̂ P - K r̂ - r̂ N r̂ = 0 (r̂ = r / c)
P_, N_, Q_, K_ = M[:2, :2], M[:2, 2:], M[2:, :2], M[2:, 2:]
rh = r / c
resL = Q_ + rh @ P_ - K_ @ rh - rh @ N_ @ rh
print('Riccati 式 Q + r̂P - K r̂ - r̂ N r̂ の残差 = %.1e' % np.abs(resL).max())
Nbar = S2 @ N_.T @ S2
kappa_num = -2 / (trpq + np.sign(trpq) * math.sqrt(trab2))
print('r̂ = κ (Q - N̄)、κ = -2/(tr(P-K) + sgn √(trab2)) との差 = %.1e' % np.abs(rh - kappa_num * (Q_ - Nbar)).max())

print('\n== (e) 色収差 ==')
chx, chy = lat.chromaticity()
src = ('MOMENTUM = 3GEV;\nDRIFT D1=(L=0.8) D2=(L=0.3);\nBEND B=(L=1.2 ANGLE=%.15g);\nQUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);\nMARK ST=();\n'
       'LINE CELL=(QF D1 B D2 QD D1 B D2);\nLINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);\nFFS;\nUSE RING;\nCELL;\nCALC;\nCHROMATICITY;\n') % (2 * math.pi / 16)
buf = io.StringIO()
run_sad_string(src, out=buf)
tot = [ln for ln in buf.getvalue().splitlines() if ln.strip().startswith('TOTAL')]
print('SAD の CHROMATICITY の TOTAL (CELL; CALC のあと。四極・6 極だけの閉形式、BEND は数えない。四極だけなら 4π dν/dδ の単位):', tot[-1].strip() if tot else '(出力なし)')
gx, gy = [float(v) for v in tot[-1].split()[1:3]] if tot else (float('nan'), float('nan'))
print('  -> 4π で割ると (BEND の寄与を含まない) (%.4f, %.4f)' % (gx / (4 * math.pi), gy / (4 * math.pi)))
nus = []
for dd in (1e-4, -1e-4):
    fd = rdttrack.one_turn(elements, mom, dd)
    zz = np.zeros(4)
    for _ in range(8):
        r = fd(zz)[0] - zz
        zz = zz - np.linalg.solve(rdttrack.jacobian(fd, zz) - np.eye(4), r)
    Jd = rdttrack.jacobian(fd, zz)
    cx_, cy_ = np.trace(Jd[:2, :2]) / 2, np.trace(Jd[2:, 2:]) / 2
    mux = np.arctan2(math.sqrt(1 - cx_ ** 2) * np.sign(Jd[0, 1]), cx_)
    muy = np.arctan2(math.sqrt(1 - cy_ ** 2) * np.sign(Jd[2, 3]), cy_)
    nus.append((mux / (2 * math.pi), muy / (2 * math.pi)))
print('追跡 + Newton + 1 周の Jacobian の μ (符号は M12 で決める) から: dν/dδ = (%.4f, %.4f) (小数部。整数部は 1 周の行列では決まらない)' % (
    (nus[0][0] - nus[1][0]) / 2e-4, (nus[0][1] - nus[1][1]) / 2e-4))
print('lat.chromaticity() (twiss_table(±δ) の位相進みの差分): (ξx, ξy) = (%.4f, %.4f)' % (chx, chy))
print('  gs (K64) の Twiss["NX",-1] の DP0 = ±1e-4 の差分 (CELL; NORFSW): (-0.7797363826, -0.8786795941) (保存した値。この環境で別に実行)')

print('\n== (e2) 四極だけのリング (偏向なし): 追跡の色収差と閉形式 -(1/4π) Σ K1L β の比較 ==')
srcq = '''MOMENTUM = 3GEV;
DRIFT D1=(L=0.8) D2=(L=0.3) B=(L=1.2);
QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);
FFS;
USE RING;
'''
ffq, _ = run_sad_string(srcq, out=io.StringIO())
elq = ffq.lat.elements
twq = ffq.lat.twiss(full=True) if False else ffq.lat.twiss()
nuq = (twq[0]['nx'] / (2 * math.pi), twq[0]['ny'] / (2 * math.pi))
sx = sy = 0.0
for k, el in enumerate(elq):
    if type(el).__name__ == 'Quadrupole':
        kl = el.K1 * el.L                  # 積分強度 (el.K1 は長さあたり)
        b0x, b1x = twq[k]['bx'], twq[k + 1]['bx']
        b0y, b1y = twq[k]['by'], twq[k + 1]['by']
        sx += kl * 0.5 * (b0x + b1x)
        sy -= kl * 0.5 * (b0y + b1y)
print('チューン (νx, νy) = (%.4f, %.4f)。薄いレンズの近似: dνx/dδ = %.4f、dνy/dδ = %.4f (= -(1/4π) Σ ±K1L β、β は要素の両端の平均)' % (
    nuq[0], nuq[1], -sx / (4 * math.pi), -sy / (4 * math.pi)))
nus = []
for dd in (1e-4, -1e-4):
    fd = rdttrack.one_turn(elq, ffq.momentum(), dd)
    Jd = rdttrack.jacobian(fd, np.zeros(4))
    cx_, cy_ = np.trace(Jd[:2, :2]) / 2, np.trace(Jd[2:, 2:]) / 2
    mux = np.arctan2(math.sqrt(1 - cx_ ** 2) * np.sign(Jd[0, 1]), cx_)
    muy = np.arctan2(math.sqrt(1 - cy_ ** 2) * np.sign(Jd[2, 3]), cy_)
    nus.append((mux / (2 * math.pi), muy / (2 * math.pi)))
print('追跡の Jacobian から dν/dδ = (%.4f, %.4f)' % ((nus[0][0] - nus[1][0]) / 2e-4, (nus[0][1] - nus[1][1]) / 2e-4))
lq = ffq.lat.chromaticity()
print('lat.chromaticity() = (%.4f, %.4f)' % lq)

print('\n== (f) ヘルプの V = P B R6 H の H (拡張 Twiss) ==')
J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
I2 = np.eye(2)
worst = 0.0
for trial in range(5):
    Hx = rng.normal(0, 0.3, (2, 2))
    Hy = rng.normal(0, 0.3, (2, 2))
    a2 = 1 - np.linalg.det(Hx) - np.linalg.det(Hy)
    a = math.sqrt(a2)
    H = np.block([[(1 - np.linalg.det(Hx) / (1 + a)) * I2, Hx @ J2 @ Hy.T @ J2 / (1 + a), -Hx],
                  [Hy @ J2 @ Hx.T @ J2 / (1 + a), (1 - np.linalg.det(Hy) / (1 + a)) * I2, -Hy],
                  [-J2 @ Hx.T @ J2, -J2 @ Hy.T @ J2, a * I2]])
    J6 = np.zeros((6, 6))
    for i in range(3):
        J6[2 * i:2 * i + 2, 2 * i:2 * i + 2] = J2
    worst = max(worst, np.abs(H.T @ J6 @ H - J6).max())
print('H (ヘルプの式、TeX 版) は 5 組のランダムな Hx, Hy (a^2 + det Hx + det Hy = 1) でシンプレクティック: 最大の |H^T J H - J| = %.1e' % worst)
