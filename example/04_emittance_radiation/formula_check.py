# -*- coding: utf-8 -*-
"""04_emittance_radiation の式の照合。
  1. FODO_ring_V4MV (RINGQ) の放射積分 I1..I5 を、SAD を使わずに 3x3 行列 (x, x', delta) で自前に積分する。
  2. 式 U0, tau_i, eps_x, sigma_delta, alpha, nu_s, sigma_z を計算する。
  3. FFS_RADINT.*.txt と FFS_EMITTANCE.*.txt から読んだ gs・pysad の値と並べて formula_check.txt に書く。
使い方: py -3.12 formula_check.py   (SAD は走らせない。先に run_example.py で .txt を作っておく)"""
import os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- ラティス (FODO_ring_V4MV.sad と同じ。K1 は積分強度なので k = K1/L)
E_GEV = 3.0
TH = 2 * np.pi / 24
L_B = 1.2
RHO = L_B / TH
EPS = 0.1308996938995747 * TH          # SAD の E1/E2 は ANGLE に対する比。端の角は E1*ANGLE
C_RING = 105.6
#            種別, 長さ, k
CELL = [('q', 0.3, 0.7 / 0.3), ('d', 0.5), ('b', 1.2), ('d', 0.5), ('q', 0.3, -0.6 / 0.3),
        ('d', 0.1), ('d', 0.2), ('d', 1.0), ('d', 0.2), ('d', 0.1)]   # QF D2 B D2 QD DS SD D1 SF DS
V_RF, F_RF = 4.0e6, 5.0e8


def drift(l): return np.array([[1, l, 0], [0, 1, 0], [0, 0, 1.0]])


def quad(k, l):
    if k > 0:
        w = np.sqrt(k)
        return np.array([[np.cos(w * l), np.sin(w * l) / w, 0], [-w * np.sin(w * l), np.cos(w * l), 0], [0, 0, 1.0]])
    w = np.sqrt(-k)
    return np.array([[np.cosh(w * l), np.sinh(w * l) / w, 0], [w * np.sinh(w * l), np.cosh(w * l), 0], [0, 0, 1.0]])


def bend(l):
    c, s = np.cos(l / RHO), np.sin(l / RHO)
    return np.array([[c, RHO * s, RHO * (1 - c)], [-s / RHO, c, s], [0, 0, 1.0]])


def edge():
    return np.array([[1, 0, 0], [np.tan(EPS) / RHO, 1, 0], [0, 0, 1.0]])


def mat(e):
    return drift(e[1]) if e[0] == 'd' else quad(e[2], e[1])


def one_cell():
    M = np.eye(3)
    for e in CELL:
        if e[0] == 'b':
            M = edge() @ M; M = bend(e[1]) @ M; M = edge() @ M
        else:
            M = mat(e) @ M
    return M


def twiss_step(T, m):
    b, a, g = T
    c, s, d, e = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    return np.array([c * c * b - 2 * c * s * a + s * s * g, -c * d * b + (c * e + s * d) * a - s * e * g,
                     d * d * b - 2 * d * e * a + e * e * g])


def integrals(n=400, edge_coeff=-2.0):
    """I1..I5。edge_coeff は I4 の端の項 edge_coeff*eta*tan(eps)/rho^2 (1 端あたり) の係数。"""
    M = one_cell()
    A = M[:2, :2]
    mu = np.arccos(0.5 * np.trace(A))
    eta = np.linalg.solve(np.eye(2) - A, M[:2, 2])
    b0 = A[0, 1] / np.sin(mu); a0 = (A[0, 0] - A[1, 1]) / (2 * np.sin(mu)); g0 = (1 + a0 * a0) / b0
    X = np.array([eta[0], eta[1], 1.0]); T = np.array([b0, a0, g0])
    I = np.zeros(5)

    def H(X, T): return T[2] * X[0] ** 2 + 2 * T[1] * X[0] * X[1] + T[0] * X[1] ** 2
    for e in CELL:
        if e[0] == 'b':
            for _ in range(2 if False else 1):
                I[3] += edge_coeff * X[0] * np.tan(EPS) / RHO ** 2
                m = edge(); X = m @ X; T = twiss_step(T, m)
            ds = e[1] / n
            for _ in range(n):
                m = bend(ds); X1 = m @ X; T1 = twiss_step(T, m)
                I[0] += 0.5 * (X[0] + X1[0]) / RHO * ds
                I[1] += ds / RHO ** 2
                I[2] += ds / RHO ** 3
                I[3] += 0.5 * (X[0] + X1[0]) / RHO ** 3 * ds
                I[4] += 0.5 * (H(X, T) + H(X1, T1)) / RHO ** 3 * ds
                X, T = X1, T1
            I[3] += edge_coeff * X[0] * np.tan(EPS) / RHO ** 2
            m = edge(); X = m @ X; T = twiss_step(T, m)
        else:
            m = mat(e); X = m @ X; T = twiss_step(T, m)
    return I * 24, mu / 2 / np.pi * 24


# ---- 物理定数
C_LIGHT = 299792458.0
RE = 2.8179403262e-15
MC2_GEV = 0.51099895e-3
CG = 4 * np.pi / 3 * RE / MC2_GEV ** 3          # m/GeV^3 (= 8.846e-5)
CQ = 55 / (32 * np.sqrt(3)) * 3.8615926796e-13  # m (= 3.832e-13)
GAMMA = E_GEV / MC2_GEV
T0 = C_RING / C_LIGHT


def formulas(I1, I2, I3, I4, I5, alpha=None):
    U0 = CG * E_GEV ** 4 * I2 / (2 * np.pi) * 1e9               # eV
    Jx = 1 - I4 / I2; Jy = 1.0; Je = 2 + I4 / I2
    tau = [2 * E_GEV * 1e9 * T0 / (J * U0) for J in (Jx, Jy, Je)]
    ex = CQ * GAMMA ** 2 * I5 / (Jx * I2)
    sd = np.sqrt(CQ * GAMMA ** 2 * I3 / (Je * I2))
    a = I1 / C_RING if alpha is None else alpha
    f0 = C_LIGHT / C_RING                                      # beta=1 近似 (差 1e-8)
    h = F_RF / f0
    cphi = np.sqrt(1 - (U0 / V_RF) ** 2)
    nus = np.sqrt(h * a * V_RF * cphi / (2 * np.pi * E_GEV * 1e9))
    sz = a * C_LIGHT * sd / (2 * np.pi * nus * f0)
    return dict(U0=U0, Jx=Jx, Je=Je, taux=tau[0], tauy=tau[1], tauz=tau[2], ex=ex, sd=sd, alpha=a, nus=nus, sz=sz)


def num(s): return float(s.replace('D', 'E'))


def parse_radint(path):
    t = open(path, encoding='utf-8', errors='replace').read()
    g = lambda k: num(re.search(k + r'\s*=\s*([-+0-9.Ee]+)', t).group(1))
    r = dict(I1=g('I1X'), I2=g('I2 '), I3=g('I3 '), I4=g('I4U'), I5=g('I5U'), U0=g(r'ENERGY LOSS PER TURN \(kV\)') * 1e3,
             Jx=g('DU=1-I4U/I2 '), taux=g(r'TAUU \(msec\)') * 1e-3, tauy=g(r'TAUV \(msec\)') * 1e-3,
             tauz=g(r'TAUE \(msec\)') * 1e-3, sd=g(r'ENERGY SPREAD \(SIG\) DE/E'), ex=g('EMITTANCE U '),
             alpha=g('MOMENTUM COMPACTION FAC X '), nus=g('NU-S '), sz=g(r'NATURAL BUNCH LENGTH \(mm\)') * 1e-3)
    return r


def parse_emit(path):
    t = open(path, encoding='utf-8', errors='replace').read()
    r = {}
    m = re.search(r'U0 =\s*([-+0-9.Ee]+) MV', t); r['U0'] = num(m.group(1)) * 1e6
    m = re.search(r'alpha =\s*([-+0-9.Ee]+)', t); r['alpha'] = num(m.group(1))
    m = re.search(r'Synchrotron frequency\s*=\s*([-+0-9.Ee]+) Hz', t); fs = num(m.group(1))
    m = re.search(r'f0 =\s*([-+0-9.Ee]+) Hz', t); r['nus'] = fs / num(m.group(1))
    m = re.search(r'Damping time \(sec\):\s*X :\s*([-+0-9.Ee]+)\s+Y :\s*([-+0-9.Ee]+)\s+Z :\s*([-+0-9.Ee]+)', t)
    r['taux'], r['tauy'], r['tauz'] = [num(x) for x in m.groups()]
    m = re.search(r'partition number:\s*X :\s*([-+0-9.Ee]+)\s+Y :\s*([-+0-9.Ee]+)\s+Z :\s*([-+0-9.Ee]+)', t)
    r['Jx'], r['Jy'], r['Je'] = [num(x) for x in m.groups()]
    m = re.search(r'Emittance X\s*=\s*([-+0-9.Ee]+) m', t); r['ex'] = num(m.group(1))
    m = re.search(r'Energy spread\s*=\s*([-+0-9.Ee]+)', t); r['sd'] = num(m.group(1))
    m = re.search(r'Bunch Length\s*=\s*([-+0-9.Ee]+) mm', t); r['sz'] = num(m.group(1)) * 1e-3
    return r


def main():
    I, qx = integrals()
    I_old, _ = integrals(edge_coeff=-1.0)
    out = []
    P = out.append
    P('自前の積分 (SAD なし): Qx = %.9f (SAD CALC 6.696329)' % qx)
    P('I1..I5 (I4 の端の項 係数 -2)  = ' + ' '.join('%.7g' % v for v in I))
    P('I1..I5 (I4 の端の項 係数 -1)  = ' + ' '.join('%.7g' % v for v in I_old))
    ra = {k: parse_radint(os.path.join(HERE, 'FFS_RADINT.%s.txt' % k)) for k in ('gs', 'pysad')}
    em = {k: parse_emit(os.path.join(HERE, 'FFS_EMITTANCE.%s.txt' % k)) for k in ('gs', 'pysad')}
    P('RADINT gs    I1..I5 = %.7g %.7g %.7g %.7g %.7g' % tuple(ra['gs'][k] for k in ('I1', 'I2', 'I3', 'I4', 'I5')))
    P('RADINT pysad I1..I5 = %.7g %.7g %.7g %.7g %.7g' % tuple(ra['pysad'][k] for k in ('I1', 'I2', 'I3', 'I4', 'I5')))
    fa = formulas(*I_old)            # 端の項 -1 (SAD RADINT と同じ I4)
    fb = formulas(*I)                # 端の項 -2
    fc = formulas(*I_old, alpha=em['gs']['alpha'])   # alpha だけ EMIT の値
    names = [('U0', 'U0 [eV]'), ('Jx', 'Jx'), ('taux', 'tau_x [s]'), ('tauy', 'tau_y [s]'), ('tauz', 'tau_z [s]'),
             ('ex', 'eps_x [m]'), ('sd', 'sigma_delta'), ('alpha', 'alpha_c'), ('nus', 'nu_s'), ('sz', 'sigma_z [m]')]
    P('')
    P('| 量 | 式 (I4 端 -1) | 式 (I4 端 -2) | RADINT gs | RADINT pysad | EMIT gs | EMIT pysad |')
    P('|---|---|---|---|---|---|---|')
    for k, nm in names:
        def f(d): return ('%.6g' % d[k]) if k in d else '-'
        P('| %s | %.6g | %.6g | %s | %s | %s | %s |' % (nm, fa[k], fb[k], f(ra['gs']), f(ra['pysad']), f(em['gs']), f(em['pysad'])))
    P('')
    P('相対差 (式 I4 端 -1 と RADINT gs):')
    for k, nm in names:
        if k in ra['gs']:
            P('  %-12s %+.2e' % (nm, fa[k] / ra['gs'][k] - 1))
    P('相対差 (式 I4 端 -2 と EMIT gs):')
    for k, nm in names:
        if k in em['gs']:
            P('  %-12s %+.2e' % (nm, fb[k] / em['gs'][k] - 1))
    P('nu_s を EMIT の alpha で計算した式の値: nu_s=%.6g (EMIT gs %.6g), sigma_z=%.6g m (EMIT gs %.6g)' % (
        fc['nus'], em['gs']['nus'], fc['sz'], em['gs']['sz']))
    g = em['gs']
    szz = g['alpha'] * C_LIGHT * g['sd'] / (2 * np.pi * g['nus'] * C_LIGHT / C_RING)
    P('sigma_z = alpha*c*sigma_delta/(2 pi f_s) に EMIT gs の alpha, sigma_delta, f_s を入れた値: %.6g m (EMIT gs の Bunch Length %.6g m)' % (szz, g['sz']))
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'formula_check.txt'), 'w', encoding='utf-8').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
