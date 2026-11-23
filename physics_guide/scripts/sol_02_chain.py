# -*- coding: utf-8 -*-
"""SOL + QUAD + SOL (BOUND = 1) の追跡の全体 (09_solenoid.md §7): pysad の追跡と、独立な基準の比較。

基準 = [端のフリンジ exp(:f:) の厳密な流れ (tsfrin、edge_04 で f の流れと確認済み)] → [四極 + 一様な Bz の全 H の正準方程式
H = -pz(P) + K (x² - y²)/2、P = p_c - a、a = (-b y/2, b x/2) を DOP853 で積分。z' = 1 - p/pz、-dv L] → [出口のフリンジ]。
K64 の本体は H2 の厳密解 + ΔH の Strang 型の分割 (ndiv 個) なので、差は分割の誤差 (振幅と ndiv による)。
DISFRIN = 1 (フリンジなし) の追跡は、フリンジなしの基準と一致するはず。
"""
import io
import math
import warnings
warnings.filterwarnings("ignore")
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6 as t6

brho = 3e9 / 299792458.0           # Bρ [T m] (3 GeV の電子)
b, K, L = 1.3, 0.9, 0.4            # b_z = BZ/Bρ [1/m]、K = K1/L [1/m²]、長さ [m]
m_e, P0 = 0.51099895e6, 3e9


def run(c, disfrin):
    src = """MOMENTUM = 3GEV;
SOL S1=(BOUND=1 BZ=%.15g GEO=1 DISFRIN=%d);
QUAD Q=(L=%.15g K1=%.15g DISFRIN=1);
SOL S2=(BOUND=1 BZ=0 DISFRIN=%d);
MARK ST=();
LINE ONE=(ST S1 Q S2);
FFS;
USE ONE;
""" % (b * brho, disfrin, L, K * L, disfrin)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return C.track_one(ffs, c)


def tsf(c, dbz):
    q = SimpleNamespace(c=np.array(c, float))
    t6._tsfrin(q, dbz)
    return q.c


def dvf(dl):
    pr = (1 + dl) * P0 / m_e
    h = math.sqrt(1 + pr * pr)
    p0 = P0 / m_e
    h0 = math.sqrt(1 + p0 * p0)
    return 1 - (pr / h) / (p0 / h0)


def ref(c0, fr):
    x, px, y, py, z, dl = c0
    p = 1 + dl
    s = np.array([x, px / p, y, py / p, z, dl])
    if fr:
        s = tsf(s, b)
    x, qx, y, qy, z, _ = s

    def rhs(s_, v):
        x, px, y, py, z = v
        Px, Py = px + b * y / 2, py - b * x / 2
        pz = math.sqrt(p * p - Px * Px - Py * Py)
        return [Px / pz, (b / 2) * Py / pz - K * x, Py / pz, -(b / 2) * Px / pz + K * y, 1 - p / pz]
    v = solve_ivp(rhs, (0, L), [x, qx * p, y, qy * p, z], rtol=1e-13, atol=1e-18, method='DOP853').y[:, -1]
    s = np.array([v[0], v[1] / p, v[2], v[3] / p, v[4] - dvf(dl) * L, dl])
    if fr:
        s = tsf(s, -b)
    return np.array([s[0], s[1] * p, s[2], s[3] * p, s[4], dl])


print('条件: 3 GeV の電子、b_z = %.1f /m (BZ = %.4f T)、四極 K = %.1f /m² (K1 = %.2f、L = %.1f m)、入射は正準変数 (x, px, y, py, z, δ)' % (b, b * brho, K, K * L, L))
print('  (ndiv = 1 + int(L hypot(K L, b_z)/0.2) = %d)' % (1 + int(L * math.hypot(K * L, b) / 0.2)))
for c in ([0.003, 0.002, -0.002, 0.0015, 0.0, 0.01], [0.001, 0.0005, 0.001, -0.0005, 0.0, 0.0], [0.001, 0.0005, 0.001, -0.0005, 0.0, -0.02]):
    t0, t1 = run(c, 0), run(c, 1)
    r0, r1 = ref(c, True), ref(c, False)
    print('  入射 %s' % c)
    print('    フリンジあり:   |追跡 - 基準| の最大 = %.2e (x, px, y, py, z: %s)' % (np.abs(t0 - r0).max(), ' '.join('%.1e' % v for v in np.abs(t0 - r0)[:5])))
    print('    DISFRIN = 1:    |追跡 - 基準| の最大 = %.2e' % np.abs(t1 - r1).max())
    print('    フリンジの効果 (基準): %s   追跡 (DISFRIN 0 - 1): %s' % (' '.join('%.2e' % v for v in (r0 - r1)[:5]), ' '.join('%.2e' % v for v in (t0 - t1)[:5])))
