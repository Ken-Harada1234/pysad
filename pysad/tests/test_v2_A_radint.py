# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): RADINT (INTGRL.f) — B06・B07・B12・B13 と、量子寿命の単位 (R1)、色収差の積分 (R2)、MULT/SOL の警告 (R3)。

根拠は、SAD を使わない自前の放射積分 (3x3 行列 (x, x', δ) で分散と Twiss を進めて台形則で積分) と、定義式。
条件: 3 GeV、12 セルのリング (曲率の違う 2 種類の偏向磁石 B1・B2 を含む)、周長 62.4 m、RF 500 MHz。
"""
import io
import math
import re

import numpy as np
import pytest

from pysad import ffscmd, radint

TH1, L1 = 0.4, 1.2
TH2, L2 = 2 * math.pi / 12 - 0.4, 0.6
E1R = 0.3            # 端の角 (ANGLE に対する比)

LAT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.55);
SEXT  SF = (L=0.2 K2=%(k2)s) ST = (L=0 K2=%(kt)s);
BEND  B1 = (L=1.2 ANGLE=0.4 E1=0.3 E2=0.3) B2 = (L=0.6 ANGLE=%(th2)r E1=0.3 E2=0.3);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=2.5E6 FREQ=5.0E8 PHI=0) RF2 = (L=0 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (%(first)s QF D2 B1 D2 QD D2 B2 D1 SF ST DS);
LINE RING = (ORIGIN RF %(rf2)s 12*CELL);
FFS;
USE RING; CELL; CALCULATE;
%(extra)s
RADINT;
'''


def _run(k2='0', first='', rf2='RF2', extra='', kt='1.5'):
    src = LAT % dict(k2=k2, th2=TH2, first=first, rf2=rf2, extra=extra, kt=kt)
    buf = io.StringIO()
    ffs, _I = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _g(t, key):
    m = re.search(key + r'\s*=\s*([-+0-9.EeD]+|NaN)', t)
    return float(m.group(1).replace('D', 'E'))


# ------------------------------------------------------------------ 自前の放射積分
def _integrals(n=400):
    def drift(l):
        return np.array([[1, l, 0], [0, 1, 0], [0, 0, 1.0]])

    def quad(k, l):
        if k > 0:
            w = np.sqrt(k)
            return np.array([[np.cos(w * l), np.sin(w * l) / w, 0], [-w * np.sin(w * l), np.cos(w * l), 0], [0, 0, 1.0]])
        w = np.sqrt(-k)
        return np.array([[np.cosh(w * l), np.sinh(w * l) / w, 0], [w * np.sinh(w * l), np.cosh(w * l), 0], [0, 0, 1.0]])

    def bend(l, rho):
        c, s = np.cos(l / rho), np.sin(l / rho)
        return np.array([[c, rho * s, rho * (1 - c)], [-s / rho, c, s], [0, 0, 1.0]])

    def edge(e, rho):
        return np.array([[1, 0, 0], [np.tan(e) / rho, 1, 0], [0, 0, 1.0]])
    cell = [('q', 0.3, 0.6 / 0.3), ('d', 0.5), ('b', L1, TH1), ('d', 0.5), ('q', 0.3, -0.55 / 0.3), ('d', 0.5),
            ('b', L2, TH2), ('d', 1.0), ('d', 0.2), ('d', 0.1)]

    def mat(e):
        if e[0] == 'd':
            return drift(e[1])
        if e[0] == 'q':
            return quad(e[2], e[1])
        rho = e[1] / e[2]
        ed = edge(E1R * e[2], rho)
        return ed @ bend(e[1], rho) @ ed
    M = np.eye(3)
    for e in cell:
        M = mat(e) @ M
    A = M[:2, :2]
    mu = np.arccos(0.5 * np.trace(A))
    eta = np.linalg.solve(np.eye(2) - A, M[:2, 2])
    b0 = A[0, 1] / np.sin(mu)
    a0 = (A[0, 0] - A[1, 1]) / (2 * np.sin(mu))
    X = np.array([eta[0], eta[1], 1.0])
    T = np.array([b0, a0, (1 + a0 * a0) / b0])

    def tw(T, m):
        b, a, g = T
        c, s, d, e = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
        return np.array([c * c * b - 2 * c * s * a + s * s * g, -c * d * b + (c * e + s * d) * a - s * e * g,
                         d * d * b - 2 * d * e * a + e * e * g])

    def H(X, T):
        return T[2] * X[0]**2 + 2 * T[1] * X[0] * X[1] + T[0] * X[1]**2
    I = np.zeros(5)
    gam = 0.0
    for e in cell:
        if e[0] == 'b':
            rho = e[1] / e[2]
            ed = edge(E1R * e[2], rho)
            I[3] -= X[0] * np.tan(E1R * e[2]) / rho**2
            X = ed @ X
            T = tw(T, ed)
            ds = e[1] / n
            m = bend(ds, rho)
            for _ in range(n):
                X1 = m @ X
                T1 = tw(T, m)
                I[0] += 0.5 * (X[0] + X1[0]) / rho * ds
                I[1] += ds / rho**2
                I[2] += ds / abs(rho)**3
                I[3] += 0.5 * (X[0] + X1[0]) / rho**3 * ds
                I[4] += 0.5 * (H(X, T) + H(X1, T1)) / abs(rho)**3 * ds
                gam += 0.5 * (T[2] + T1[2]) * ds
                X, T = X1, T1
            I[3] -= X[0] * np.tan(E1R * e[2]) / rho**2
            X = ed @ X
            T = tw(T, ed)
        else:
            ns = 50
            m = mat((e[0], e[1] / ns) + tuple(e[2:]))
            for _ in range(ns):
                T1 = tw(T, m)
                gam += 0.5 * (T[2] + T1[2]) * e[1] / ns
                X, T = m @ X, T1
    return I * 12, -gam * 12 / (4 * np.pi)


def test_per_bend_radius_v2(v2_mode):
    """B06: 曲率の違う BEND があるリングで、v2 の I1..I5 は自前の積分と一致する"""
    I, _xi = _integrals()
    _f, out = _run()
    got = [_g(out, 'I1X'), _g(out, 'I2 '), _g(out, 'I3 '), _g(out, 'I4U'), _g(out, 'I5U')]
    for a, b in zip(got, I):
        assert abs(a / b - 1) < 2e-5, (got, I)


def test_per_bend_radius_k64(k64_mode):
    """K64: 全 BEND に最後の BEND (B2) の L/ANGLE を使うので、I1・I4・I5 (分散の伝搬を通る量) がずれる"""
    I, _xi = _integrals()
    _f, out = _run()
    got = [_g(out, 'I1X'), _g(out, 'I4U'), _g(out, 'I5U')]
    ref = [I[0], I[3], I[4]]
    assert max(abs(a / b - 1) for a, b in zip(got, ref)) > 1e-3


def test_sextupole_first_and_thin(v2_mode):
    """B07: 先頭が SEXT でも NaN にならず (K64 は DL=0 で NaN)、L=0 の六極は 0.2 m として足されない。
    軌道が 0 なら六極 (K2≠0) は I1..I5 を変えない"""
    I, _xi = _integrals()
    _f, out0 = _run()
    _f, out = _run(k2='3', first='SF')
    assert 'DL=0' not in out and 'SET L=0.2' not in out and 'NaN' not in out
    for key in ('I2 ', 'I3 '):
        assert abs(_g(out, key) / _g(out0, key) - 1) < 1e-12
    assert abs(_g(out, 'I2 ') / I[1] - 1) < 2e-5


def test_sextupole_first_k64_nan(k64_mode):
    _f, out = _run(k2='3', first='SF')
    radint._STATE.DL = 0.0          # プロセス内の前回の値を消す (K64 の初回と同じ状態)
    _f, out = _run(k2='3', first='SF')
    assert 'SET L=0.2' in out


def test_double_precision_constants(v2_mode):
    """B12: ELOSS = Cγ E⁴ I2/(2π) (Cγ = 4π/3·r_e/(mc²)³、倍精度)"""
    _f, out = _run()
    m = re.search(r'NOT ENOUGH RF VOLTS;  ELOSS,RFV =\s*([0-9.]+)', out)
    res, _msgs = radint.compute(_f, _f.ensure_rows())
    i2 = res['RI'][2]
    eloss_kv = 2.0e6 * radint._coeff() * i2 * res['energy']**4
    cg = 4 * math.pi / 3 * radint.ELRADI / (radint.ELMASS * 1e-9)**3
    ref_kv = cg * res['energy']**4 * i2 / (2 * math.pi) * 1e6
    assert abs(eloss_kv / ref_kv - 1) < 1e-13
    assert abs(radint._coeff() / radint.COEFF - 1 - 1.1516e-8) < 1e-11      # K64 との差は相対 1.15e-8


def test_cavity_voltage_is_summed(v2_mode):
    """B13: 空洞 2 台 (2.5 MV + 1.5 MV) の RF 電圧は 4000 kV (K64 は最後の 1500 kV)"""
    _f, out = _run()
    assert abs(_g(out, r'RF VOLTAGE \(kV\)') - 4000.0) < 1e-6
    # ν_s は h·α·V·cosφ/(2πE) の平方根
    I, _xi = _integrals()
    u0 = _g(out, r'ENERGY LOSS PER TURN \(kV\)')
    al = _g(out, 'MOMENTUM COMPACTION FAC X')
    h = _g(out, 'HARMONIC NUMBER')
    nus = math.sqrt(h * al * 4000.0e3 * math.sqrt(1 - (u0 / 4000.0)**2) / (2 * math.pi * 3.0e9))
    assert abs(_g(out, 'NU-S') / nus - 1) < 2e-3


def test_cavity_voltage_k64_last(k64_mode):
    _f, out = _run()
    assert abs(_g(out, r'RF VOLTAGE \(kV\)') - 1500.0) < 1e-6


def test_quantum_lifetime_hours(v2_mode):
    """R1: 量子寿命 [時間] = τ_E/2·exp(ξ)/ξ/3600、ξ = (δmax/σδ)²/2 (K64 は exp(2.3026·秒) を出す)"""
    _f, out = _run(extra='')
    if 'QUANTUM LIFETIME' not in out:
        pytest.skip('ξ >= 100')
    taue = _g(out, r'TAUE \(msec\)') * 1e-3
    xi = 0.5 * (_g(out, 'DEMAX/E') / _g(out, r'ENERGY SPREAD \(SIG\) DE/E'))**2
    ref = taue * math.exp(xi) / 2 / xi / 3600.0
    assert abs(_g(out, r'QUANTUM LIFETIME\(hours\)') / ref - 1) < 2e-2       # 表示桁 (δmax 6 桁、σδ 5 桁) の範囲


def test_chromaticity_is_tune_derivative(v2_mode):
    """R2: GX・GY は運動量をずらした CALC のチューンの微分 dν/dδ (六極 K2≠0 の寄与も入る)。
    K64 の -I9/4π は磁石の中だけの ∮γ ds で、この値にならない"""
    for k2, kt in (('0', '0'), ('3', '1.5')):
        src = LAT % dict(k2=k2, th2=TH2, first='', rf2='', kt=kt, extra='''
DP0=1e-4; CALC; a={Twiss["NX","$$$"],Twiss["NY","$$$"]}/2/Pi;
DP0=-1e-4; CALC; b={Twiss["NX","$$$"],Twiss["NY","$$$"]}/2/Pi;
Print[(a-b)/2e-4]; DP0=0; CALC;''')
        buf = io.StringIO()
        ffscmd.run_sad_string(src, out=buf)
        out = buf.getvalue()
        fd = [l for l in out.splitlines() if l.startswith('{')][0].strip('{}').split(',')
        assert abs(_g(out, 'LINEAR CHROMATICITY GX') - float(fd[0])) < 2e-5
        assert abs(_g(out, 'LINEAR CHROMATICITY GY') - float(fd[1])) < 2e-5
    _I, xi = _integrals()
    assert abs(float(fd[0]) - xi) > 5.0          # 六極を入れた色収差は、線形の ∮γ ds (-4.3) と別物


def test_chromaticity_k64_is_magnet_only_integral(k64_mode):
    _f, out = _run(kt='0')
    _I, xi = _integrals()
    assert abs(_g(out, 'LINEAR CHROMATICITY GX') / xi - 1) > 0.2


def test_mult_warning(v2_mode):
    """R3: 場のある MULT は RADINT に入らないので、その旨を表示する"""
    src = (LAT % dict(k2='0', th2=TH2, first='', rf2='', extra='', kt='0')).replace(
        'MARK  ORIGIN', 'MULT M1 = (L=0.2 K1=0.01);\nMARK  ORIGIN').replace('SF ST DS', 'M1 ST DS')
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    assert 'MULT/SOL with fields are NOT included' in buf.getvalue()


def _alpha(rfsw):
    src = LAT % dict(k2='0', th2=TH2, first='', rf2='RF2', kt='0', extra='RAD; %s;\ne=Emittance[];\nPrint["AL ",MomentumCompaction/.e];' % rfsw)
    buf = io.StringIO()
    ffscmd.run_sad_string(src.replace('\nRADINT;\n', '\n'), out=buf)
    return float([l for l in buf.getvalue().splitlines() if l.startswith('AL ')][0][3:])


def test_momentum_compaction_is_i1_over_c_v2(v2_mode):
    """ALPHA: Emittance[] の MomentumCompaction = I1/C (自前の積分)。RFSW OFF は −(R56 + Σ R5j η_j)/C、
    RFSW ON はシンクロトロンチューンから 4 sin²(μ/2) で求める (空洞や放射の分で 1e-3 以内)"""
    I, _xi = _integrals()
    ref = I[0] / 62.4
    assert abs(_alpha('NORFSW') / ref - 1) < 2e-4
    assert abs(_alpha('RFSW') / ref - 1) < 1.5e-3


def test_momentum_compaction_k64_uses_r56_only(k64_mode):
    """K64: RFSW OFF は −R56/C だけ (入口の分散の分が抜ける)"""
    I, _xi = _integrals()
    ref = I[0] / 62.4
    assert abs(_alpha("NORFSW") / ref - 1) > 1e-3          # このリングでは 0.14% (分散の大きい所から始まるリングでは 3.4%: 記録の ALPHA)
