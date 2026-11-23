# -*- coding: utf-8 -*-
"""偏向磁石のあるリングの色収差 (12_linear_optics.md §7、付録 B)。

同じ 8 セルのリングで、次の 3 つの dν/dδ を比べる。
(1) 厳密なハミルトニアンの写像 (独立な計算): H = -(1 + h x) sqrt(p² - px² - py²) + h (x + h x²/2) + K (x² - y²)/2 を DOP853 で積分し、
    δ = ±1e-4 のそれぞれで閉軌道 (Newton 法) を求め、1 周のヤコビアンの位相から ν の差分。フリンジなし。
(2) pysad の追跡 (閉軌道 + ヤコビアンの位相の差分)。
(3) pysad の lat.chromaticity()。
さらに、単独の偏向磁石の追跡 (DISFRIN = 1) と厳密な H の写像の差 (y ≠ 0 での Δp_y が h² L δ y の程度)。
"""
import io
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import rdttrack

h = 2 * math.pi / 16 / 1.2           # 偏向磁石の曲率 [1/m]


def drift(s, L, d):
    x, px, y, py = s
    p = 1 + d
    pz = math.sqrt(p * p - px * px - py * py)
    return np.array([x + px / pz * L, px, y + py / pz * L, py])


def quad(s, L, k, d):
    def f(t, u):
        x, px, y, py = u
        p = 1 + d
        pz = math.sqrt(p * p - px * px - py * py)
        return [px / pz, -k * x, py / pz, k * y]
    return solve_ivp(f, [0, L], s, rtol=1e-13, atol=1e-16, method='DOP853').y[:, -1]


def bend(s, L, d):
    """H = -(1 + h x) pz + h (x + h x²/2) (曲率 h の扇形、端面・フリンジなし)"""
    def f(t, u):
        x, px, y, py = u
        p = 1 + d
        pz = math.sqrt(p * p - px * px - py * py)
        return [(1 + h * x) * px / pz, -(-h * pz + h + h * h * x), (1 + h * x) * py / pz, 0.0]
    return solve_ivp(f, [0, L], s, rtol=1e-13, atol=1e-16, method='DOP853').y[:, -1]


def turn(s, d):
    for _ in range(8):
        s = quad(s, 0.3, 0.35 / 0.3, d)
        s = drift(s, 0.8, d)
        s = bend(s, 1.2, d)
        s = drift(s, 0.3, d)
        s = quad(s, 0.3, -0.33 / 0.3, d)
        s = drift(s, 0.8, d)
        s = bend(s, 1.2, d)
        s = drift(s, 0.3, d)
    return s


def jac(f, c, e=1e-7):
    J = np.zeros((4, 4))
    for i in range(4):
        a = c.copy()
        b = c.copy()
        a[i] += e
        b[i] -= e
        J[:, i] = (f(a) - f(b)) / (2 * e)
    return J


def tunes_of(J):
    out = []
    for a in (0, 2):
        B = J[a:a + 2, a:a + 2]
        mu = math.acos(np.trace(B) / 2)
        if B[0, 1] < 0:
            mu = -mu
        out.append((mu % (2 * math.pi)) / (2 * math.pi))
    return out


print('== (1) 厳密なハミルトニアンの写像 (DOP853、フリンジなし) ==')
nus = []
for d in (1e-4, -1e-4):
    f = lambda s, d=d: turn(s, d)
    c = np.zeros(4)
    for _ in range(6):
        J = jac(f, c)
        c = c - np.linalg.solve(J - np.eye(4), f(c) - c)
    nus.append(tunes_of(jac(f, c)))
xi_exact = [(nus[0][i] - nus[1][i]) / 2e-4 for i in range(2)]
print('  dν/dδ = (%.5f, %.5f)   (δ = ±1e-4 の差分。閉軌道は分散で動かす)' % tuple(xi_exact))

print('\n== (2)(3) pysad: 追跡 (閉軌道 + ヤコビアン) と lat.chromaticity() ==')


def run(extraB='', extraQ=''):
    src = """MOMENTUM = 3GEV;
DRIFT D1=(L=0.8) D2=(L=0.3);
BEND B=(L=1.2 ANGLE=%.15g %s);
QUAD QF=(L=0.3 K1=0.35 %s) QD=(L=0.3 K1=-0.33 %s);
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);
FFS;
USE RING;
""" % (2 * math.pi / 16, extraB, extraQ, extraQ)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    el = ffs.lat.elements
    mom = ffs.momentum()
    nn = []
    for dd in (1e-4, -1e-4):
        fd = rdttrack.one_turn(el, mom, dd)
        zz = np.zeros(4)
        for _ in range(8):
            r = fd(zz)[0] - zz
            zz = zz - np.linalg.solve(rdttrack.jacobian(fd, zz) - np.eye(4), r)
        Jd = rdttrack.jacobian(fd, zz)
        cx, cy = np.trace(Jd[:2, :2]) / 2, np.trace(Jd[2:, 2:]) / 2
        mux = np.arctan2(math.sqrt(1 - cx**2) * np.sign(Jd[0, 1]), cx)
        muy = np.arctan2(math.sqrt(1 - cy**2) * np.sign(Jd[2, 3]), cy)
        nn.append((mux / 2 / math.pi, muy / 2 / math.pi))
    xi = [(nn[0][i] - nn[1][i]) / 2e-4 for i in range(2)]
    return xi, ffs.lat.chromaticity()


for lab, eb, eq in (('既定 (BEND・QUAD とも DISFRIN = 0)', '', ''), ('DISFRIN = 1 (BEND・QUAD)', 'DISFRIN=1', 'DISFRIN=1')):
    xi, lc = run(eb, eq)
    print('  %-40s: 追跡 dν/dδ = (%.5f, %.5f)、lat.chromaticity() = %s' % (lab, xi[0], xi[1], np.array2string(np.asarray(lc, float), precision=5)))
print('  (四極だけのリングでは追跡と lat.chromaticity() が一致する (opt_01 (e2))。偏向磁石があると lat.chromaticity() だけが合わない。)')

print('\n== 単独の偏向磁石 (L = 1.2 m、ANGLE = 2π/16、DISFRIN = 1) の追跡と厳密な H の写像の差 ==')
ffs, _ = C.make_bend(1.2, 2 * math.pi / 16, 'DISFRIN=1')


def cmp_(c):
    out = C.track_one(ffs, c)
    mine = bend(np.array(c[:4], float), 1.2, c[5])
    return out[:4] - mine


print('  x のみ (1 mm)          : 差 (x, px, y, py) = %s' % np.array2string(cmp_([1e-3, 0, 0, 0, 0, 0]), precision=2))
print('  δ のみ (1e-3)          : 差 = %s' % np.array2string(cmp_([0, 0, 0, 0, 0, 1e-3]), precision=2))
print('  y のみ (1 mm)          : 差 = %s' % np.array2string(cmp_([0, 0, 1e-3, 0, 0, 0]), precision=2))
print('  y = 1 mm、δ = 1e-3     : 差 = %s、h² L δ y = %.2e' % (np.array2string(cmp_([0, 0, 1e-3, 0, 0, 1e-3]), precision=2), h * h * 1.2 * 1e-3 * 1e-3))
print('  (x のみ・δ のみでは差は 0。y ≠ 0 で p_y に差が出る。h² L δ y と同じ大きさだが、原因 (硬い端の Maxwell フリンジか、本体の y 依存の項か) の内訳は未確認。)')
