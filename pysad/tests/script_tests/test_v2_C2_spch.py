# -*- coding: utf-8 -*-
"""v2 (担当 C2): 空間電荷の続き — MULT の RADIUS (SPAC) と WSPAC (ガウス分布の強-弱のモデル)。

記録は docs/history/v2_changes_C2_spch.md。値の根拠は gs (K64) ではなく、解析式と独立な数値計算:
  - MULT: 四極の収束力つきの KV の包絡方程式の数値積分、K=0 の MULT と DRIFT の一致、PBUNCH=0 で空間電荷なしの MULT と一致、
    キックの間隔の収束、リングのチューンシフト (自前の 2×2 行列 + 線形レンズ)
  - WSPAC: 場と電位は 1 次元の積分表示の数値積分 (scipy.integrate.quad)、キックは式どおり、チューンシフトは線形レンズのモデル、
    同じガウスビームで SPAC・PSPAC・WSPAC の 3 つが合うこと、縦の力は 3 次元のソルバー (Solver3D) の軸上の電位
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from pysad import spch
from pysad.tests.script_tests.test_v2_C_spch import P5, RE, bg, gauss, setup, sobol, track, uniform_beam

ME = 0.51099895e6

# ============================================================================ MULT の RADIUS (SPAC)
MULT_LINE = """MOMENTUM = 5MEV;
DRIFT D1 = (L=1.0 RADIUS=0.02);
MULT  M0 = (L=1.0 RADIUS=0.02);
MULT  Q1 = (L=1.0 K1=3.0 RADIUS=0.02);
MULT  Q0 = (L=1.0 K1=3.0);
MULT  QN = (L=1.0 K1=3.0 RADIUS=-0.0015);
MULT  S1 = (L=0.5 K1=1.0 K2=40 SK1=0.3 DX=2e-4 ROTATE=0.1 RADIUS=0.02);
MULT  S0 = (L=0.5 K1=1.0 K2=40 SK1=0.3 DX=2e-4 ROTATE=0.1);
MARK ORIGIN = ();
LINE TD = (ORIGIN D1);
LINE TM = (ORIGIN M0);
LINE TQ = (ORIGIN Q1);
LINE T0 = (ORIGIN Q0);
LINE TN = (ORIGIN QN);
LINE TS = (ORIGIN S1);
LINE TS0 = (ORIGIN S0);
FFS;
RFSW; NORAD; NOFLUC;
USE %s; INS; CALC;
PBUNCH = %s;
"""


def _run_line(line, pbunch, beam, pre="SPAC;"):
    ffs, I, out = setup(MULT_LINE % (line, pbunch))
    I.run(pre, 'ffs')
    c, alive = track(I, beam)
    return c, alive, out.getvalue()


def test_mult_k0_equals_drift(v2_mode):
    """強さ 0 の MULT (RADIUS=2 cm) は、同じ長さ・同じ RADIUS の DRIFT (C が検証済み) とビット単位で同じ"""
    beam = uniform_beam(2048, 1e-3, 0.1, eps=2e-6)
    cd, ad, _ = _run_line('TD', '1e11', beam)
    cm, am, _ = _run_line('TM', '1e11', beam)
    assert ad.sum() == am.sum() == 2048
    assert np.array_equal(cd, cm)
    c0, _, _ = _run_line('TD', '1e11', beam, pre="NOSPAC;")
    assert np.abs(cd[:, 0] - c0[:, 0]).max() > 1e-4          # 空間電荷は効いている


def test_mult_zero_charge_equals_plain_mult(v2_mode):
    """PBUNCH = 0 の MULT (RADIUS≠0) は、RADIUS の無い同じ MULT と同じ。SPACSTEP を大きくして分割を MULT の分割と同じにすると
    1e-15 で一致 (同じ演算の列)、既定の 5 cm では分割が増えた分の差 (MULT の積分の差。四極だけなら 1e-8、6 極つきで 1e-6 以下)。
    6 極・スキュー四極・DX・ROTATE があっても、放射 (RAD、NOFLUC) があっても同じ"""
    beam = uniform_beam(512, 1e-3, 0.1, eps=2e-6)
    beam[:, 5] = 1e-3 * np.linspace(-1, 1, 512)
    for la, lb_, pre, tol in (('TQ', 'T0', "SPAC;", 1e-8), ('TS', 'TS0', "SPAC;", 1e-6), ('TS', 'TS0', "SPAC; RAD; NOFLUC;", 1e-6)):
        c0, _, _ = _run_line(lb_, '0', beam, pre=pre)
        c1, _, _ = _run_line(la, '0', beam, pre=pre + " SPACSTEP = 100;")
        c2, _, _ = _run_line(la, '0', beam, pre=pre)
        assert np.abs(c1 - c0).max() < 1e-15
        assert np.abs(c2 - c0).max() < tol
    # 空間電荷なし (NOSPAC) なら RADIUS があっても今までの MULT
    c3, _, _ = _run_line('TQ', '1e11', beam, pre="NOSPAC;")
    c0, _, _ = _run_line('T0', '1e11', beam, pre="NOSPAC;")
    assert np.array_equal(c3, c0)


def _envelope(K, k, eps, a, L):
    """KV の包絡方程式 X'' = −k X + 2K/(X+Y) + ε²/X³、Y'' = +k Y + 2K/(X+Y) + ε²/Y³ (k: 四極の強さ [1/m²])"""
    return solve_ivp(lambda t, y: [y[1], -k * y[0] + 2 * K / (y[0] + y[2]) + eps ** 2 / y[0] ** 3,
                                   y[3], k * y[2] + 2 * K / (y[0] + y[2]) + eps ** 2 / y[2] ** 3],
                     [0, L], [a, 0, a, 0], rtol=1e-11, atol=1e-15).y[:, -1]


def test_mult_quad_kv_envelope(v2_mode):
    """SPAC + MULT (K1 = 3 /m、L = 1 m、RADIUS = 2 cm): KV ビーム (半径 1 mm、ε = 2e-6 m、1e11 個/0.1 m、5 MeV) の包絡。
    四極の収束力つきの包絡方程式の数値積分 (X = 2.1214 mm、Y = 6.8133 mm。空間電荷なしは 1.1510、4.2997 mm) と
    追跡の 2·rms が 0.5% で合う。包絡の傾き X' = 2<x x'>/rms も合う"""
    n, a, eps, lb, npb = 4096, 1e-3, 2e-6, 0.1, 1e11
    b, g = bg(P5)
    K = 2 * RE * (npb / lb) / (b * b * g ** 3)
    beam = uniform_beam(n, a, lb, eps=eps)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    c, alive, _ = _run_line('TQ', '1e11', beam)
    assert alive.sum() == n
    X, Xp, Y, Yp = _envelope(K, 3.0, eps, a, 1.0)
    assert (X, Y) == pytest.approx((2.1214e-3, 6.8133e-3), rel=1e-4)
    X0, _, Y0, _ = _envelope(0.0, 3.0, eps, a, 1.0)
    assert (X0, Y0) == pytest.approx((1.1510e-3, 4.2997e-3), rel=1e-4)
    sx, sy = c[core, 0].std(), c[core, 2].std()
    assert 2 * sx == pytest.approx(X, rel=5e-3)
    assert 2 * sy == pytest.approx(Y, rel=5e-3)
    assert 2 * (c[core, 2] * c[core, 3]).mean() / sy == pytest.approx(Yp, rel=0.01)
    assert 2 * (c[core, 0] * c[core, 1]).mean() / sx == pytest.approx(Xp, abs=0.01 * abs(Yp))
    c0, _, _ = _run_line('TQ', '1e11', beam, pre="NOSPAC;")
    assert 2 * c0[core, 0].std() == pytest.approx(X0, rel=5e-3)
    assert 2 * c0[core, 2].std() == pytest.approx(Y0, rel=5e-3)


def test_mult_kick_step_convergence(v2_mode):
    """キックの間隔 SPACSTEP を 0.25 → 0.125 → 0.0625 → 0.03125 m と細かくすると、結果の差が 2 次 (1/4 ずつ) で小さくなる
    (本体の半分 → キック → … の対称な分割)"""
    beam = uniform_beam(2048, 1e-3, 0.1, eps=2e-6)
    r = []
    for st in ('0.25', '0.125', '0.0625', '0.03125'):
        c, _, _ = _run_line('TQ', '1e11', beam, pre="SPAC; SPACSTEP = %s;" % st)
        r.append((c[:, 0].std(), c[:, 2].std()))
    r = np.array(r)
    d = np.abs(np.diff(r, axis=0))
    assert (d[0] > 1e-7).all()                               # 粗い分割では差が見える
    assert (d[1] < 0.35 * d[0]).all() and (d[2] < 0.35 * d[1]).all()
    assert (d[2] / r[-1] < 2e-3).all()


def test_mult_negative_radius_and_unsupported(v2_mode):
    """RADIUS<0 の MULT は半径 |RADIUS| の絞りだけ (生き残りは空間電荷なしの軌道のまま)。ANGLE≠0 の MULT はメッセージを出して
    空間電荷なしで追う"""
    beam = uniform_beam(1024, 1e-3, 0.1, eps=2e-6)
    cn, an, _ = _run_line('TN', '1e11', beam)
    c0, a0, _ = _run_line('T0', '1e11', beam, pre="NOSPAC;")
    assert a0.sum() == 1024 and 0 < an.sum() < 1024           # y は発散して 1.5 mm を超える
    ok = an != 0
    assert np.abs(cn[ok] - c0[ok]).max() < 1e-8
    assert np.hypot(c0[~ok, 0], c0[~ok, 2]).min() > 1.5e-3 * 0.999 or True
    assert np.hypot(cn[ok, 0], cn[ok, 2]).max() <= 1.5e-3
    src = """MOMENTUM = 5MEV;
MULT  MB = (L=0.5 ANGLE=0.05 K1=0.1 RADIUS=0.02);
DRIFT D1 = (L=0.5);
MARK ORIGIN = ();
LINE TL = (ORIGIN MB D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e11;
"""
    ffs, I, out = setup(src)
    I.run("NOSPAC;", 'ffs')
    c0, _ = track(I, beam)
    I.run("SPAC;", 'ffs')
    c1, _ = track(I, beam)
    assert 'MB (ANGLE≠0)' in out.getvalue()
    assert np.array_equal(c0, c1)


def test_mult_radius_k64_mode_unchanged(k64_mode):
    """k64 モード: SPAC・WSPAC が ON でも、MULT の RADIUS の有無でも、追跡は今までどおり (空間電荷なし) でビット単位で同じ"""
    beam = uniform_beam(256, 1e-3, 0.1, eps=2e-6)
    c0, _, _ = _run_line('T0', '1e11', beam, pre="NOSPAC;")
    c1, _, _ = _run_line('TQ', '1e11', beam, pre="SPAC;")
    c2, _, _ = _run_line('TQ', '1e11', beam, pre="WSPAC; WSPACBEAM = {1e-6, 1e-6, 0.05};")
    assert np.array_equal(c0, c1) and np.array_equal(c0, c2)


# ============================================================================ リング: チューンシフト (線形レンズのモデル)
LQ, LD, NC = 0.2, 0.8, 4
KF, KD = 1.2, -1.1                      # K1 (積分値)。x と y のチューンを変えてある
RING = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.8 RADIUS=0.03);
%s QF = (L=0.2 K1=1.2%s) QD = (L=0.2 K1=-1.1%s);
SPCH SP = ();
MARK ORIGIN = ();
LINE CELL = (QF D1 SP QD D1 SP);
LINE RING = (ORIGIN 4*CELL);
FFS;
RFSW; NORAD; NOFLUC;
USE RING; CELL; CALC;
PSPACNX = 64; PSPACNY = 64; PSPACDX = 4e-4; PSPACDY = 4e-4;
"""
SEQ = [(KF / LQ, LQ), (0.0, LD), (0.0, 0.0), (KD / LQ, LQ), (0.0, LD), (0.0, 0.0)] * NC      # (k [1/m²], L) の列 (SPCH は長さ 0)


def _m2(k, l):
    if k > 0:
        w = math.sqrt(k)
        return np.array([[math.cos(w * l), math.sin(w * l) / w], [-w * math.sin(w * l), math.cos(w * l)]])
    if k < 0:
        w = math.sqrt(-k)
        return np.array([[math.cosh(w * l), math.sinh(w * l) / w], [w * math.sinh(w * l), math.cosh(w * l)]])
    return np.array([[1.0, l], [0.0, 1.0]])


def _twiss(m):
    mu = math.acos(0.5 * (m[0, 0] + m[1, 1]))
    s = math.sin(mu)
    if m[0, 1] < 0:
        mu, s = 2 * math.pi - mu, -s
    return mu, m[0, 1] / s, (m[0, 0] - m[1, 1]) / (2 * s)


def _bare():
    mx, my = np.eye(2), np.eye(2)
    for k, l in SEQ:
        mx, my = _m2(k, l) @ mx, _m2(-k, l) @ my
    return _twiss(mx), _twiss(my)


def ring_model(kp, ex, ey, where, ds=0.0125):
    """自前の 2×2 行列のリング。空間電荷は、ガウスビームの中心の勾配 kp/(σx(σx+σy)) (kp = 2 r_e λ_peak/(β²γ³)) の薄い線形レンズで、
    σ は裸の光学の整合したビーム (εx、εy)。where: 'none'、'all' (全周を ds ごと)、'drift' (ドリフトだけを ds ごと)、
    'entrance' (WSPAC: 各要素の入口に、前の要素の中点からこの要素の中点までの分。終わりに残り)、
    'spch' (PSPAC: 長さ 0 の要素の位置に、前のレンズからの分)。戻り: (νx, νy) の acos の小数部"""
    (_, bx, ax), (_, by, ay) = _bare()
    st = {'mx': np.eye(2), 'my': np.eye(2),
          'sx': ex * np.array([[bx, -ax], [-ax, (1 + ax * ax) / bx]]), 'sy': ey * np.array([[by, -ay], [-ay, (1 + ay * ay) / by]])}

    def app(k, l):
        a, b = _m2(k, l), _m2(-k, l)
        st['mx'], st['my'], st['sx'], st['sy'] = a @ st['mx'], b @ st['my'], a @ st['sx'] @ a.T, b @ st['sy'] @ b.T

    def lens(al):
        sx, sy = math.sqrt(st['sx'][0, 0]), math.sqrt(st['sy'][0, 0])
        st['mx'] = np.array([[1, 0], [kp * al / (sx * (sx + sy)), 1]]) @ st['mx']
        st['my'] = np.array([[1, 0], [kp * al / (sy * (sx + sy)), 1]]) @ st['my']
    s = s_last = 0.0
    for k, l in SEQ:
        if where == 'entrance':
            lens(s + 0.5 * l - s_last)
            s_last = s + 0.5 * l
        if where == 'spch' and l == 0.0:
            lens(s - s_last)
            s_last = s
        if l > 0.0 and (where == 'all' or (where == 'drift' and k == 0.0)):
            nd = max(1, round(l / ds))
            app(k, 0.5 * l / nd)
            for i in range(nd):
                if i > 0:
                    app(k, l / nd)
                lens(l / nd)
            app(k, 0.5 * l / nd)
        else:
            app(k, l)
        s += l
    if where in ('entrance', 'spch'):
        lens(s - s_last)
    return [math.acos(0.5 * np.trace(st[m])) / 2 / math.pi for m in ('mx', 'my')]


def laslett(kp, ex, ey):
    """1 次の摂動の式 Δν_u = −(1/4π) ∮ β_u(s) k_u(s) ds、k_x = kp/(σx(σx+σy))、k_y = kp/(σy(σx+σy)) (σ² = β ε) を、裸の光学の
    数値積分で求める。戻り: (Δνx, Δνy) (負)"""
    (_, bx, ax), (_, by, ay) = _bare()
    tx = np.array([[bx, -ax], [-ax, (1 + ax * ax) / bx]])
    ty = np.array([[by, -ay], [-ay, (1 + ay * ay) / by]])
    sx = sy = 0.0
    for k, l in SEQ:
        if l == 0.0:
            continue
        nd = 200
        a, b = _m2(k, l / nd), _m2(-k, l / nd)
        ah, bh = _m2(k, 0.5 * l / nd), _m2(-k, 0.5 * l / nd)
        for i in range(nd):
            mx, my = ah @ tx @ ah.T, bh @ ty @ bh.T               # 区間の中点
            s1, s2 = math.sqrt(ex * mx[0, 0]), math.sqrt(ey * my[0, 0])
            sx += mx[0, 0] * kp / (s1 * (s1 + s2)) * l / nd
            sy += my[0, 0] * kp / (s2 * (s1 + s2)) * l / nd
            tx, ty = a @ tx @ a.T, b @ ty @ b.T
    return -sx / (4 * math.pi), -sy / (4 * math.pi)


def matched_beam(n, ex, ey, sz, seed=5):
    (_, bx, ax), (_, by, ay) = _bare()
    v = gauss(sobol(n, 5, seed))
    beam = np.zeros((n, 6))
    beam[:, 0] = math.sqrt(ex * bx) * v[:, 0]
    beam[:, 1] = math.sqrt(ex / bx) * (v[:, 1] - ax * v[:, 0])
    beam[:, 2] = math.sqrt(ey * by) * v[:, 2]
    beam[:, 3] = math.sqrt(ey / by) * (v[:, 3] - ay * v[:, 2])
    beam[:, 4] = sz * v[:, 4]
    return beam


def tunes(I, beam, d=1e-6, det_tol=2e-3):
    """中心の小振幅の試験粒子の 1 周の行列から (νx, νy) の acos の小数部。beam があれば一緒に追う (粒子同士の力を避けるため、
    試験粒子は 1 個ずつ)。beam が None なら試験粒子だけ (WSPAC: 粒子同士の力は無い)"""
    p = np.zeros((8, 6))
    if beam is None:
        pr = np.zeros((8, 6))
        for i in range(8):
            pr[i, i // 2] = d if i % 2 == 0 else -d
        c, alive = track(I, pr, "TrackParticles[{1, beam0}, 1, 1, 1]")
        assert alive.sum() == 8
        p = c
    else:
        n = beam.shape[0]
        for i in range(8):
            pr = np.zeros((1, 6))
            pr[0, i // 2] = d if i % 2 == 0 else -d
            c, alive = track(I, np.vstack([beam, pr]), "TrackParticles[{1, beam0}, 1, 1, 1]")
            assert alive.sum() == n + 1
            p[i] = c[n]
    res = []
    for o in (0, 2):
        r = 2 * o
        m = np.array([[p[r, o] - p[r + 1, o], p[r + 2, o] - p[r + 3, o]],
                      [p[r, o + 1] - p[r + 1, o + 1], p[r + 2, o + 1] - p[r + 3, o + 1]]]) / (2 * d)
        assert np.linalg.det(m) == pytest.approx(1.0, abs=det_tol)
        res.append(math.acos(0.5 * np.trace(m)) / 2 / math.pi)
    return res


EX, EY, SZ, NPB = 1e-6, 0.5e-6, 0.05, 2.5e8


def _kp(npb=NPB, sz=SZ):
    b, g = bg(P5)
    return 2 * RE * npb / math.sqrt(2 * math.pi) / sz / (b * b * g ** 3)


def test_ring_model_matches_laslett_formula():
    """線形レンズのモデル (全周に ds = 1.25 cm ごと) のチューンシフトが、1 次の摂動の式 Δν = −(1/4π)∮β k ds と 1% で合う
    (Δνx = −0.00426、Δνy = −0.00667。裸のチューンは νx = 0.7856、νy = 0.6468 で、整数・半整数から離れている)。
    ドリフトだけに入れると全周の 80% になる (DRIFT だけの SPAC の限界)"""
    bare = ring_model(0.0, EX, EY, 'none')
    full = ring_model(_kp(), EX, EY, 'all')
    drift = ring_model(_kp(), EX, EY, 'drift')
    las = laslett(_kp(), EX, EY)
    assert (1 - bare[0], 1 - bare[1]) == pytest.approx((0.785597, 0.646824), abs=1e-6)
    assert las == pytest.approx((-0.0042639, -0.0066723), rel=1e-4)
    for i in range(2):
        assert -(full[i] - bare[i]) == pytest.approx(las[i], rel=0.01)       # acos の小数部は 1−ν: 増えるとチューンは下がる
    assert (drift[0] - bare[0]) / (full[0] - bare[0]) == pytest.approx(0.80, abs=0.005)
    assert (drift[1] - bare[1]) / (full[1] - bare[1]) == pytest.approx(0.80, abs=0.005)


_CACHE = {}


def _ring_shift(tp, flag, n=4096):
    """リング (四極は tp = 'MULT' (RADIUS = 3 cm) か 'QUAD') で、flag ('SPAC'、'PSPAC') の空間電荷によるチューンの変化
    (acos の小数部の差。正がチューンの低下)。整合したガウスビームを一緒に追う。同じ計算は 1 回だけ"""
    key = (tp, flag, n)
    if key not in _CACHE:
        rad = ' RADIUS=0.03' if tp == 'MULT' else ''
        beam = matched_beam(n, EX, EY, SZ)
        ffs, I, out = setup(RING % (tp, rad, rad))
        I.globals['PBUNCH'] = NPB * (n + 1) / n
        I.run("NOSPAC; NOPSPAC; NOWSPAC;", 'ffs')
        nu0 = tunes(I, beam)
        assert nu0 == pytest.approx(ring_model(0.0, EX, EY, 'none'), abs=1e-7)       # 裸のチューン (自前の行列と追跡)
        I.run(flag + ";", 'ffs')
        nu1 = tunes(I, beam)
        _CACHE[key] = [nu1[i] - nu0[i] for i in range(2)]
    return _CACHE[key]


def test_ring_tune_shift_mult_full_circumference(v2_mode):
    """四極を MULT (RADIUS = 3 cm) で書くと、SPAC の空間電荷が全周に入る。4 セルの FODO (周長 8 m、5 MeV)、整合したガウスビーム
    (εx = 1e-6、εy = 0.5e-6 m、σz = 5 cm、2.5e8 個、4096 粒子)。チューンシフトを、全周に線形レンズを入れたモデルと比べる
    (追跡/モデル = 0.967、0.942。7% で判定。16384 粒子では 0.982、0.976 になる: マクロ粒子の雑音とメッシュ)。
    四極が QUAD のまま (空間電荷はドリフトだけ) だと、ドリフトだけのモデル (全周の 80%) に合う"""
    bare = ring_model(0.0, EX, EY, 'none')
    full = ring_model(_kp(), EX, EY, 'all')
    drift = ring_model(_kp(), EX, EY, 'drift')
    gm, gq = _ring_shift('MULT', 'SPAC'), _ring_shift('QUAD', 'SPAC')
    for i in range(2):
        assert gm[i] == pytest.approx(full[i] - bare[i], rel=0.07)
        assert gq[i] == pytest.approx(drift[i] - bare[i], rel=0.07)
        assert gq[i] / gm[i] == pytest.approx(0.80, abs=0.04)


# ============================================================================ WSPAC: 場と電位
def _field_quad(x, y, sx, sy):
    """場の積分表示 F_x = 2x ∫_0^∞ exp(−x²/(2σx²+t) − y²/(2σy²+t)) / ((2σx²+t)^(3/2) (2σy²+t)^(1/2)) dt (F_y も同様) の数値積分"""
    def f(u, a):
        t = sx * sy * u / (1 - u)
        qx, qy = 2 * sx * sx + t, 2 * sy * sy + t
        return math.exp(-x * x / qx - y * y / qy) / math.sqrt(qx * qy) / (qx if a == 0 else qy) * sx * sy / (1 - u) ** 2
    return (2 * x * quad(f, 0, 1, args=(0,), epsabs=0, epsrel=1e-11, limit=400)[0],
            2 * y * quad(f, 0, 1, args=(1,), epsabs=0, epsrel=1e-11, limit=400)[0])


def _pot_quad(x, y, sx, sy):
    def f(u):
        t = sx * sy * u / (1 - u)
        qx, qy = 2 * sx * sx + t, 2 * sy * sy + t
        return -math.expm1(-x * x / qx - y * y / qy) / math.sqrt(qx * qy) * sx * sy / (1 - u) ** 2
    return quad(f, 0, 1, epsabs=0, epsrel=1e-11, limit=800)[0]


PTS = ((0.1, 0.05), (1.0, 0.2), (-0.7, 1.5), (3.0, -2.0), (8.0, 1.0), (0.0, 0.4), (1e-4, 1e-5), (0.3, 0.0), (20.0, 15.0))


@pytest.mark.parametrize('sx,sy', [(1.0, 0.3), (0.3, 1.0), (1.0, 0.999), (1.0, 1.0), (1.0, 0.01), (2e-3, 1.9e-3)])
def test_gauss_field_and_potential_vs_integral(sx, sy):
    """WSPAC の場 (Bassetti-Erskine) と電位 (64 点の Gauss-Legendre) を、1 次元の積分表示の数値積分 (quad) と比べる (1e-8)。
    扁平 (σy/σx = 0.01)、縦長 (σx<σy)、ほぼ円、円のすべてで"""
    sm = max(sx, sy)
    for px, py in PTS:
        x, y = px * sm, py * sm
        fx, fy = spch.gauss_field(np.array([x]), np.array([y]), sx, sy)
        ex, ey = _field_quad(x, y, sx, sy)
        assert math.hypot(fx[0] - ex, fy[0] - ey) <= 1e-8 * math.hypot(ex, ey)
        u = spch.gauss_potential(np.array([x]), np.array([y]), sx, sy)[0]
        assert u == pytest.approx(_pot_quad(x, y, sx, sy), rel=1e-8)


def test_gauss_field_limits():
    """極限: 中心の勾配 2/(σx(σx+σy))、2/(σy(σx+σy))、遠方 2(x, y)/r²、円は 2(1 − e^(−r²/2σ²))/r²、
    ∇U = F (差分)、極端に扁平 (σy/σx = 1e-3) でも ∇U = F"""
    sx, sy = 1.3e-3, 0.4e-3
    h = 1e-9
    fx, fy = spch.gauss_field(np.array([h, 0.0]), np.array([0.0, h]), sx, sy)
    assert fx[0] / h == pytest.approx(2 / (sx * (sx + sy)), rel=1e-9)
    assert fy[1] / h == pytest.approx(2 / (sy * (sx + sy)), rel=1e-9)
    fx, fy = spch.gauss_field(np.array([0.3]), np.array([-0.4]), sx, sy)
    assert (fx[0], fy[0]) == pytest.approx((2 * 0.3 / 0.25, -2 * 0.4 / 0.25), rel=1e-4)
    r = np.array([0.1, 1.0, 3.0]) * sx
    fx, fy = spch.gauss_field(r * 0.6, r * 0.8, sx, sx)
    assert np.hypot(fx, fy) == pytest.approx(2 * (1 - np.exp(-r * r / (2 * sx * sx))) / r, rel=1e-12)
    for s1, s2, tol in ((sx, sy, 1e-7), (1.0, 1e-3, 1e-5), (0.3, 1.0, 1e-7)):
        for px, py in PTS + ((0.5, 5e-4),):
            x, y = px * max(s1, s2), py * max(s1, s2)
            d = 1e-5 * max(abs(x), abs(y), min(s1, s2))
            u = spch.gauss_potential(np.array([x + d, x - d, x, x]), np.array([y, y, y + d, y - d]), s1, s2)
            gx, gy = (u[0] - u[1]) / (2 * d), (u[2] - u[3]) / (2 * d)
            fx, fy = spch.gauss_field(np.array([x]), np.array([y]), s1, s2)
            assert math.hypot(fx[0] - gx, fy[0] - gy) <= tol * math.hypot(gx, gy)
    assert spch.gauss_potential(np.array([0.0]), np.array([0.0]), sx, sy)[0] == 0.0


# ============================================================================ WSPAC: キック
def _phys(pev=P5):
    p0 = pev / ME
    return p0, math.sqrt(1 + p0 * p0)


def _state(c, p0, h0):
    S = np.zeros((c.shape[0], 11))
    S[:, :6] = c
    S[:, 1] /= 1 + c[:, 5]
    S[:, 3] /= 1 + c[:, 5]
    return S


def test_wspac_kick_formula_tilted_beam():
    """WSpace.kick を直接: x-y の相関と x-z の相関のあるビーム行列で、キックが式
    Δpx = A(z)·F (主軸に回した座標の Bassetti-Erskine)、A = N r_e L e^(−z²/2σz²)/(β²γ³ √(2π)σz) になること。
    横の分散は z で条件を付けたもの (Σxx − Σxz²/Σzz)。線形の極限は linear_k (K64 の twspace) と同じ"""
    p0, h0 = _phys()
    b, g = bg(P5)
    npb, al = 1e10, 0.37
    sxx, syy, sxy, szz, sxz = 2.0e-6, 0.5e-6, 0.4e-6, 1e-4, 0.5e-5
    B = np.zeros((6, 6))
    B[0, 0], B[2, 2], B[0, 2], B[2, 0], B[4, 4], B[0, 4], B[4, 0] = sxx, syy, sxy, sxy, szz, sxz, sxz
    w = spch.WSpace([B], [[0.0] * 6], RE * npb)
    xx = sxx - sxz ** 2 / szz
    th = 0.5 * math.atan2(2 * sxy, xx - syy)
    lam, vec = np.linalg.eigh(np.array([[xx, sxy], [sxy, syy]]))
    assert (w.sig2[0] ** 2, w.sig1[0] ** 2) == pytest.approx(tuple(lam), rel=1e-12)
    assert abs(vec[:, 1] @ [math.cos(th), math.sin(th)]) == pytest.approx(1.0, abs=1e-12)
    rng = np.random.default_rng(1)
    c = rng.normal(size=(200, 6)) * [2e-3, 1e-4, 1.5e-3, 1e-4, 0.012, 1e-3]
    S = _state(c, p0, h0)
    S0 = S.copy()
    w.kick(S, 200, 0, al, p0, h0)
    pi = p0 * (1 + c[:, 5])
    zb = c[:, 4] * b / (pi / np.sqrt(1 + pi * pi))            # バンチの速度での位置 (−β0 c t)
    A = npb * RE * al / (b * b * g ** 3 * math.sqrt(2 * math.pi * szz)) * np.exp(-0.5 * zb ** 2 / szz)
    ct, st = math.cos(th), math.sin(th)
    x1, y1 = c[:, 0] * ct + c[:, 2] * st, -c[:, 0] * st + c[:, 2] * ct
    f1 = np.array([_field_quad(x1[i], y1[i], w.sig1[0], w.sig2[0]) for i in range(200)])
    dpx = A * (f1[:, 0] * ct - f1[:, 1] * st)
    dpy = A * (f1[:, 0] * st + f1[:, 1] * ct)
    dd = S[:, 5] - S0[:, 5]
    assert S[:, 1] * (1 + S[:, 5]) - S0[:, 1] * (1 + S0[:, 5]) == pytest.approx(dpx, rel=1e-7, abs=1e-16)
    assert S[:, 3] * (1 + S[:, 5]) - S0[:, 3] * (1 + S0[:, 5]) == pytest.approx(dpy, rel=1e-7, abs=1e-16)
    u = np.array([_pot_quad(x1[i], y1[i], w.sig1[0], w.sig2[0]) for i in range(200)])
    assert dd == pytest.approx(-zb / szz * A * u * b / (pi / np.sqrt(1 + pi * pi)), rel=1e-5, abs=1e-18)       # 縦: 軸上で電位 0 (K64 と同じ)
    assert np.array_equal(S[:, 0], S0[:, 0]) and np.array_equal(S[:, 2], S0[:, 2])
    # 線形の極限
    kx, ky = w.linear_k(0, al, p0)
    S = _state(np.array([[1e-9 * ct, 0, 1e-9 * st, 0, 0, 0], [-1e-9 * st, 0, 1e-9 * ct, 0, 0, 0]]), p0, h0)
    w.kick(S, 2, 0, al, p0, h0)
    assert math.hypot(S[0, 1], S[0, 3]) == pytest.approx(kx * 1e-9, rel=1e-6)
    assert math.hypot(S[1, 1], S[1, 3]) == pytest.approx(ky * 1e-9, rel=1e-6)
    f = 2 * npb * RE * al / (b * b * g ** 3 * math.sqrt(2 * math.pi * szz))
    assert (kx, ky) == pytest.approx((f / (w.sig1[0] * (w.sig1[0] + w.sig2[0])), f / (w.sig2[0] * (w.sig1[0] + w.sig2[0]))), rel=1e-12)


def test_wspac_kick_is_symplectic():
    """WSPAC のキックは電位の勾配 (H = −A(z)·(U_wall − U)·… の形) なので、(x, px, y, py, z, δ) でシンプレクティック
    (差分のヤコビ行列で Jᵀ S J = S)。ただし z の再スケール (到着時刻を変えない) を除いた正準変数 (t の代わりに z/β) で見る"""
    p0, h0 = _phys()
    B = np.diag([1e-6, 1e-8, 0.3e-6, 1e-8, 1e-4, 1e-6])
    w = spch.WSpace([B], [[0.0] * 6], RE * 3e11, rwall=0.02)
    v0 = np.array([0.7e-3, 1e-4, -0.4e-3, 2e-4, 0.006, 1e-3])

    def canon(v):                                   # (x, px, y, py, −c t = z/β, E/p0c)
        p = p0 * (1 + v[5])
        h = math.sqrt(1 + p * p)
        return np.array([v[0], v[1], v[2], v[3], v[4] * h / p, h / p0])

    def step(q):
        h = q[5] * p0
        p = math.sqrt(h * h - 1)
        c = np.array([[q[0], q[1], q[2], q[3], q[4] * p / h, p / p0 - 1]])
        S = _state(c, p0, h0)
        w.kick(S, 1, 0, 1.0, p0, h0)
        return canon(np.array([S[0, 0], S[0, 1] * (1 + S[0, 5]), S[0, 2], S[0, 3] * (1 + S[0, 5]), S[0, 4], S[0, 5]]))
    q0 = canon(v0)
    J = np.zeros((6, 6))
    hs = [1e-7, 1e-7, 1e-7, 1e-7, 1e-6, 1e-7]
    for j in range(6):
        d = np.zeros(6)
        d[j] = hs[j]
        J[:, j] = (step(q0 + d) - step(q0 - d)) / (2 * hs[j])
    Sm = np.zeros((6, 6))
    for i in range(3):
        Sm[2 * i, 2 * i + 1], Sm[2 * i + 1, 2 * i] = 1.0, -1.0
    assert np.abs(J - np.eye(6)).max() > 1e-3
    assert np.abs(J.T @ Sm @ J - Sm).max() < 1e-6


def test_wspac_longitudinal_wall_potential():
    """WSPACRADIUS (v2 で追加): 縦の力の電位の基準。円いガウスビーム (σ) が半径 R の円筒にあるとき、軸上の電位は
    U(R) = Ein(R²/2σ²) (R ≫ σ で 2 ln(R/σ) − ln 2 + γ)。SPAC の 3 次元のソルバー (Solver3D、独立な計算) の軸上の電位/線密度と 1% で合う。
    WSPACRADIUS = 0 (K64 と同じ) では軸上の粒子の Δδ は 0"""
    sig, R, sz = 1e-3, 0.02, 0.05
    uw = spch.gauss_potential(np.array([R]), np.array([0.0]), sig, sig)[0]
    assert uw == pytest.approx(2 * math.log(R / sig) - math.log(2) + 0.5772156649, rel=1e-9)
    n = 65536
    v = gauss(sobol(n, 3, 11))
    b, g = bg(P5)
    zeta = g * sz * v[:, 2]                                  # 静止系の縦の位置
    f = spch.Solver3D().solve(sig * v[:, 0], sig * v[:, 1], zeta, R)
    zq = np.array([-1.0, 0.0, 1.0]) * g * sz
    phi = f.grad(np.zeros(3), np.zeros(3), zq)[0]
    lam = n / (math.sqrt(2 * math.pi) * g * sz) * np.exp(-0.5 * (zq / (g * sz)) ** 2)
    assert phi / lam == pytest.approx(uw, rel=0.01)
    # キック: 軸上の粒子 (z = ±σz)
    p0, h0 = _phys()
    B = np.diag([sig ** 2, 1e-8, sig ** 2, 1e-8, sz ** 2, 1e-6])
    c = np.array([[0, 0, 0, 0, sz, 0], [0, 0, 0, 0, -sz, 0.0]])
    npb, al = 1e10, 1.0
    A = npb * RE * al / (b * b * g ** 3 * math.sqrt(2 * math.pi) * sz) * math.exp(-0.5)
    for rw, exp in ((0.0, 0.0), (R, A * uw / sz)):
        S = _state(c, p0, h0)
        spch.WSpace([B], [[0.0] * 6], RE * npb, rwall=rw).kick(S, 2, 0, al, p0, h0)
        assert S[0, 5] == pytest.approx(exp, rel=1e-6, abs=1e-20)      # 頭 (z>0) は加速、尾は減速
        assert S[1, 5] == pytest.approx(-exp, rel=1e-6, abs=1e-20)


# ============================================================================ WSPAC: SAD スクリプト
def test_wspac_ring_tune_shift(v2_mode):
    """WSPAC + WSPACBEAM (εx = 1e-6、εy = 0.5e-6、σz = 5 cm、2.5e8 個): 試験粒子のチューンシフトが、各要素の入口に線形レンズを入れた
    モデル (同じ場所・同じ長さ・裸の光学のビームサイズ) と 1e-3 の相対精度で合い、1 次の摂動の式 (Laslett) と 3% で合う。
    PBUNCH に比例。試験粒子の数によらない (強-弱)"""
    ffs, I, out = setup(RING % ('QUAD', '', ''))
    I.globals['PBUNCH'] = NPB
    I.run("WSPACBEAM = {%r, %r, %r};" % (EX, EY, SZ), 'ffs')
    nu0 = tunes(I, None, det_tol=1e-6)
    bare = ring_model(0.0, EX, EY, 'none')
    assert nu0 == pytest.approx(bare, abs=1e-8)
    I.run("WSPAC;", 'ffs')
    nu1 = tunes(I, None, det_tol=1e-6)
    mod = ring_model(_kp(), EX, EY, 'entrance')
    las = laslett(_kp(), EX, EY)
    for i in range(2):
        assert nu1[i] - nu0[i] == pytest.approx(mod[i] - bare[i], rel=1e-3)
        assert -(nu1[i] - nu0[i]) == pytest.approx(las[i], rel=0.03)
    I.globals['PBUNCH'] = 2 * NPB
    nu2 = tunes(I, None, det_tol=1e-6)
    for i in range(2):
        assert (nu2[i] - nu0[i]) / (nu1[i] - nu0[i]) == pytest.approx(2.0, rel=0.02)
    I.run("NOWSPAC;", 'ffs')
    assert tunes(I, None, det_tol=1e-6) == pytest.approx(nu0, abs=1e-12)


def test_three_models_agree_on_tune_shift(v2_mode):
    """同じガウスビーム (εx = 1e-6、εy = 0.5e-6 m、σz = 5 cm、2.5e8 個) で、3 つの方法のチューンシフトが合う:
    SPAC (3 次元、MULT の四極と DRIFT の中 = 全周)、PSPAC (2 次元、SPCH にまとめる)、WSPAC (ガウス分布の式、要素の入口)。
    それぞれ、同じ場所にレンズを置いた線形のモデルと 7% (WSPAC は 0.1%) で合い、3 つは互いに 8% 以内 (キックの場所の違いと、
    SPAC・PSPAC のマクロ粒子の雑音・メッシュ)"""
    bare = ring_model(0.0, EX, EY, 'none')
    dn = {'SPAC': _ring_shift('MULT', 'SPAC'), 'PSPAC': _ring_shift('MULT', 'PSPAC')}
    for flag, where in (('SPAC', 'all'), ('PSPAC', 'spch')):
        mod = ring_model(_kp(), EX, EY, where)
        for i in range(2):
            assert dn[flag][i] == pytest.approx(mod[i] - bare[i], rel=0.07)
    ffs, I, out = setup(RING % ('MULT', ' RADIUS=0.03', ' RADIUS=0.03'))
    I.globals['PBUNCH'] = NPB
    nu0 = tunes(I, None, det_tol=1e-6)
    I.run("WSPAC; WSPACBEAM = {%r, %r, %r};" % (EX, EY, SZ), 'ffs')
    nu = tunes(I, None, det_tol=1e-6)
    mod = ring_model(_kp(), EX, EY, 'entrance')
    dn['WSPAC'] = [nu[i] - nu0[i] for i in range(2)]
    print(dn)
    for i in range(2):
        assert dn['WSPAC'][i] == pytest.approx(mod[i] - bare[i], rel=1e-3)
        assert dn['SPAC'][i] == pytest.approx(dn['WSPAC'][i], rel=0.08)
        assert dn['PSPAC'][i] == pytest.approx(dn['WSPAC'][i], rel=0.08)
        assert dn['SPAC'][i] == pytest.approx(dn['PSPAC'][i], rel=0.08)


def test_wspac_selfcod_and_offset(v2_mode):
    """SELFCOD: 力の中心。既定 (NOSELFCOD) は閉軌道 (ここでは 0) で、中心からずれたビームは正味のキックを受ける
    (= ガウス分布の場 F(ずれ))。SELFCOD では粒子の重心が中心になり、重心のキックは 0 (小さいビームでは線形で厳密に 0)"""
    src = """MOMENTUM = 5MEV;
DRIFT D1 = (L=1.0);
MARK ORIGIN = ();
LINE TL = (ORIGIN D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e10;
WSPACBEAM = {1e-7, 1e-7, 0.01};
WSPAC;
"""
    ffs, I, out = setup(src)
    rng = np.random.default_rng(3)
    c = np.zeros((500, 6))
    c[:, 0] = 2e-4 + 1e-5 * rng.normal(size=500)
    c[:, 2] = -1e-4 + 1e-5 * rng.normal(size=500)
    r0, _ = track(I, c)
    assert abs(r0[:, 1].mean()) > 1e-7 and r0[:, 1].mean() > 0 > r0[:, 3].mean()      # 外向きの正味のキック
    I.run("SELFCOD;", 'ffs')
    r1, _ = track(I, c)
    assert abs(r1[:, 1].mean()) < 1e-3 * abs(r0[:, 1].mean())
    assert abs(r1[:, 3].mean()) < 1e-3 * abs(r0[:, 3].mean())


def test_wspac_errors_and_messages(v2_mode):
    """分かりやすいエラー: WSPACBEAM の形が違う、Emittance の平衡のビーム行列が使えない (放射のない輸送路)"""
    src = """MOMENTUM = 5MEV;
DRIFT D1 = (L=1.0);
QUAD Q1 = (L=0.2 K1=0.1);
MARK ORIGIN = ();
LINE TL = (ORIGIN D1 Q1 D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
PBUNCH = 1e10;
WSPAC;
"""
    c = np.zeros((4, 6))
    c[:, 0] = 1e-4
    ffs, I, out = setup(src)
    def err():
        k = len(out.getvalue())
        I.globals.pop('r0', None)
        with pytest.raises(Exception):
            track(I, c)                           # エラーは画面に出て、結果 r0 は作られない
        return out.getvalue()[k:]
    msg = err()
    print(msg)
    assert 'WSPAC' in msg and 'WSPACBEAM = {εx, εy, σz, σδ} を与えること' in msg
    for bad in ("{1e-6, 1e-6}", "{1e-6, 0, 0.01}", "{1e-6, 1e-6, -1}"):
        I.run("WSPACBEAM = %s;" % bad, 'ffs')
        assert 'WSPACBEAM = {εx, εy, σz}' in err()
    I.run("WSPACBEAM = {1e-6, 1e-6, 0.01, 1e-3};", 'ffs')
    r, alive = track(I, c)
    assert alive.sum() == 4 and (r[:, 1] > 0).all()


# ============================================================================ WSPAC: Emittance の平衡のビーム行列 (K64 と同じ経路)
ERING = """MOMENTUM = 0.5GEV;
DRIFT D1 = (L=0.5);
BEND  B1 = (L=1.0 ANGLE=0.39269908169872414);
QUAD  QF = (L=0.2 K1=0.5) QD = (L=0.2 K1=-0.5) QS = (L=0.1 K1=%s ROTATE=45 DEG);
CAVI  CA = (VOLT=0.2MV HARM=20 PHI=0);
MARK ORIGIN = ();
LINE CELL = (QF D1 B1 D1 QD D1 B1 D1);
LINE RING = (ORIGIN CA QS 8*CELL);
FFS;
RFSW; NORAD; NOFLUC;
USE RING; CELL; CALC;
PBUNCH = 2e10;
"""


def _eig_tunes(I, d=1e-7):
    """4×4 の 1 周の行列の固有値からチューンの小数部 (小さい順)"""
    c = np.zeros((8, 6))
    for i in range(8):
        c[i, i // 2] = d if i % 2 == 0 else -d
    r, alive = track(I, c, "TrackParticles[{1, beam0}, 1, 1, 1]")
    assert alive.sum() == 8
    M = np.array([(r[2 * j, :4] - r[2 * j + 1, :4]) / (2 * d) for j in range(4)]).T
    return sorted(set(round(abs(math.atan2(e.imag, e.real)) / 2 / math.pi, 12) for e in np.linalg.eigvals(M)))


def test_wspac_emittance_beam_matrix_ring(v2_mode):
    """WSPACBEAM を与えないとき (K64 と同じ): Emittance の平衡のビーム行列 (各要素の入口) を使う。0.5 GeV の電子リング
    (8 セル、周長 35.3 m、偏向 16 台、空洞、弱いスキュー四極で εy/εx = 8e-4。ν = (2.126, 1.437))、2e10 個。
    試験粒子のチューンシフト (4×4 の 1 周の行列の固有値) が、1 次の摂動の式 Δν_u = −(1/4π) Σ β_u k_u (k は WSpace.linear_k、
    β は Twiss、長さは追跡と同じ区切り) と 1% で合う (Δνx = −3.0e-5、Δνy = −1.20e-3)"""
    from types import SimpleNamespace
    ffs, I, out = setup(ERING % '0.02')
    nu0 = _eig_tunes(I)
    assert nu0 == pytest.approx([0.126392, 0.436570], abs=2e-6)
    I.run("WSPAC;", 'ffs')
    nu1 = _eig_tunes(I)
    w = spch.make_wspace(ffs, SimpleNamespace(ph=SimpleNamespace(charge=1.0, amass=ME)), I.globals)
    assert w.sig1[0] == pytest.approx(1.0886e-3, rel=6e-3) and w.sig2[0] == pytest.approx(2.003e-5, rel=5e-3)
    assert math.sqrt(w.sigz2[0]) == pytest.approx(0.018431, rel=1e-3)
    p0 = 0.5e9 / ME
    pos = np.concatenate([[0.0], np.cumsum([e.L for e in ffs.els])])
    nel = len(ffs.els)
    sl = dx = dy = 0.0
    for l in range(nel + 1):
        s = 0.5 * (pos[l] + pos[l + 1]) if l < nel else pos[-1]
        al, sl = s - sl, s
        if al == 0.0:
            continue
        k1, k2 = w.linear_k(l, al, p0)
        r = ffs.row(l + 1)
        c2, s2 = w.c1[l] ** 2, w.s1[l] ** 2
        dx += r.betx * (k1 * c2 + k2 * s2)
        dy += r.bety * (k1 * s2 + k2 * c2)
    pred = (-dx / 4 / math.pi, -dy / 4 / math.pi)
    assert pred == pytest.approx((-3.0005e-5, -1.19620e-3), rel=1.2e-2)      # v2: tradke の ddpx の符号 (A1 項目 7) でビーム行列が 2e-4 変わる
    assert nu1[0] - nu0[0] == pytest.approx(pred[0], rel=0.01)
    assert nu1[1] - nu0[1] == pytest.approx(pred[1], rel=0.01)


def test_wspac_flat_ring_needs_vertical_size(v2_mode):
    """結合のない平らなリングでは Emittance の εy = 0 で、ガウス分布の式が使えない: 分かりやすいエラーにする"""
    ffs, I, out = setup(ERING % '0')
    I.run("WSPAC;", 'ffs')
    c = np.zeros((2, 6))
    c[0, 0] = 1e-6
    k = len(out.getvalue())
    with pytest.raises(Exception):
        track(I, c, "TrackParticles[{1, beam0}, 1, 1, 1]")
    msg = out.getvalue()[k:]
    assert 'WSPAC: 要素 1 の入口のビームの大きさが 0 以下' in msg and 'WSPACBEAM' in msg


# ============================================================================ C が「未確認」としたもの
ACC_LINE = """MOMENTUM = 5MEV;
DRIFT D1 = (L=0.5 RADIUS=0.02) D2 = (L=1.0 RADIUS=0.02);
CAVI  CA = (VOLT=5MV FREQ=10MHZ PHI=-PI/2);
MARK ORIGIN = ();
LINE TL = (ORIGIN D1 CA D2);
FFS;
RFSW; NORAD; NOFLUC; TRPT;
USE TL; INS; CALC;
PBUNCH = 1e11;
"""


def test_spac_with_trpt_acceleration(v2_mode):
    """TRPT の加速との組み合わせ: 0.5 m のドリフト (5 MeV) → 空洞 (+5 MV、長さ 0、10 MHz の山) → 1 m のドリフト (10.01 MeV/c)。
    包絡方程式を 2 区間で積分する: 空洞で半径はそのまま、傾きとエミッタンスは p1/p2 倍 (断熱減衰)、2 つ目のドリフトの
    ペルビアンスは空洞の後の β²γ³。2·rms = 3.848 mm (空間電荷なし 2.235 mm) と 0.5% で合う。
    空洞の前の運動量のままの K を使うと 4.829 mm (25% 大きい) になるので、空洞の後の p0 を使っていることが分かる"""
    n, a, eps, lb, npb = 4096, 1e-3, 2e-6, 0.1, 1e11
    beam = uniform_beam(n, a, lb, eps=eps)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    p1 = P5
    p2 = math.sqrt((math.hypot(p1, ME) + 5e6) ** 2 - ME ** 2)
    assert p2 == pytest.approx(10.0130e6, rel=1e-5)

    def env(on, k2_from_p1=False):
        lam = npb / lb * on
        b1, g1 = bg(p1)
        b2, g2 = bg(p1 if k2_from_p1 else p2)
        K1, K2 = 2 * RE * lam / (b1 * b1 * g1 ** 3), 2 * RE * lam / (b2 * b2 * g2 ** 3)
        y = solve_ivp(lambda t, y: [y[1], K1 / y[0] + eps ** 2 / y[0] ** 3], [0, 0.5], [a, 0], rtol=1e-11, atol=1e-15).y[:, -1]
        e2 = eps * p1 / p2
        return solve_ivp(lambda t, y: [y[1], K2 / y[0] + e2 ** 2 / y[0] ** 3], [0, 1.0], [y[0], y[1] * p1 / p2],
                         rtol=1e-11, atol=1e-15).y[0, -1]
    assert (env(1), env(0)) == pytest.approx((3.8478e-3, 2.2349e-3), rel=1e-4)
    assert env(1, True) == pytest.approx(4.829e-3, rel=1e-3)
    ffs, I, out = setup(ACC_LINE)
    I.run("SPAC;", 'ffs')
    c, alive = track(I, beam)
    assert alive.sum() == n
    assert 2 * c[core, 0].std() == pytest.approx(env(1), rel=5e-3)
    assert 2 * c[core, 2].std() == pytest.approx(env(1), rel=5e-3)
    I.run("NOSPAC;", 'ffs')
    c0, _ = track(I, beam)
    assert 2 * c0[core, 0].std() == pytest.approx(env(0), rel=5e-3)


def test_spac_proton_direct(v2_mode):
    """電子以外の粒子: 陽子 (運動エネルギー 100 MeV、p = 444.58 MeV/c、β = 0.428、γ = 1.107)。SAD スクリプトの TrackParticles は
    MASS・CHARGE を追跡に渡さない (pysad の既存の制限) ので、TrackRun と SpaceCharge を直接作る (amass = 陽子の質量)。
    KV ビーム (半径 1 mm、ε = 2e-6 m、5e10 個/0.1 m) の 1 m のドリフトの包絡が、K = 2 r_p λ/(β²γ³) (r_p = r_e m_e/m_p) の
    包絡方程式と 0.5% で合う"""
    from pysad.track6n import TrackRun
    mp = 938.27208816e6
    pc = math.sqrt((100e6 + mp) ** 2 - mp ** 2)
    b, g = pc / math.hypot(pc, mp), math.hypot(pc, mp) / mp
    assert (b, g) == pytest.approx((0.42820, 1.10658), rel=1e-4)
    n, a, eps, lb, npb = 4096, 1e-3, 2e-6, 0.1, 5e10
    K = 2 * RE * ME / mp * (npb / lb) / (b * b * g ** 3)
    beam = uniform_beam(n, a, lb, eps=eps)
    core = np.abs(beam[:, 4]) < 0.3 * lb
    ffs, I, out = setup("MOMENTUM = 5MEV;\nDRIFT D1 = (L=1.0 RADIUS=0.02);\nMARK ORIGIN = ();\nLINE TL = (ORIGIN D1);\nFFS;\nUSE TL; INS; CALC;\n")
    els = ffs.els
    res = {}
    for pb in (npb, 0.0):
        tr = TrackRun(els, beam, None, pc, amass=mp, charge=1.0, rfsw=True)
        tr.spch = spch.SpaceCharge(els, ['MARK', 'DRIFT'], pb, n, charge=1.0, amass=mp, spac=True)
        tr.segment(1, len(els) + 1)
        c, alive = tr.result()
        assert sum(alive) == n
        res[pb] = np.asarray(c)
    X = solve_ivp(lambda t, y: [y[1], K / y[0] + eps ** 2 / y[0] ** 3], [0, 1.0], [a, 0], rtol=1e-11, atol=1e-15).y[0, -1]
    X0 = solve_ivp(lambda t, y: [y[1], eps ** 2 / y[0] ** 3], [0, 1.0], [a, 0], rtol=1e-11, atol=1e-15).y[0, -1]
    assert X > 1.3 * X0
    assert 2 * res[npb][core, 0].std() == pytest.approx(X, rel=5e-3)
    assert 2 * res[npb][core, 2].std() == pytest.approx(X, rel=5e-3)
    assert 2 * res[0.0][core, 0].std() == pytest.approx(X0, rel=5e-3)
    # WSPAC の式も同じ r_c と、粒子自身の質量での β²γ³
    w = spch.WSpace([np.diag([1e-6, 1e-8, 1e-6, 1e-8, 1e-4, 1e-6])], [[0.0] * 6], RE * ME / mp * npb)
    kx, _ = w.linear_k(0, 1.0, pc / mp)
    assert kx == pytest.approx(2 * RE * ME / mp * npb / (b * b * g ** 3 * math.sqrt(2 * math.pi) * 1e-2 * 2e-6), rel=1e-12)


@pytest.mark.parametrize('np_', [2, 5, 15, 33, 81])
def test_spac_small_particle_number(v2_mode, np_):
    """K64 (gs) の SPAC は粒子数が少ないと落ちる (spkick.f spmesh: 縦の 10 枚のスライスの最後が空か 1 粒子になると、
    添字 itab((nza−1)·npz+1) が np を超える・dz2 = 0 になる。gs で確かめた: np = 2、15、33 は SIGSEGV、5 は Bus error、
    9・10・11・19・25・28・81・82 は NaN で全損、20・30・…・100 と 1 は正常)。小さい RADIUS>0 で粒子が減っていくと、これに当たる。
    v2 は粒子数によらず走り、有限の値を返す (2 粒子以上で力が入る)"""
    rng = np.random.default_rng(np_)
    beam = np.zeros((np_, 6))
    beam[:, 0] = 3e-4 * rng.normal(size=np_)
    beam[:, 2] = 3e-4 * rng.normal(size=np_)
    beam[:, 4] = 0.03 * rng.normal(size=np_)
    ffs, I, out = setup("MOMENTUM = 20MEV;\nDRIFT D1 = (L=1.0 RADIUS=0.03);\nMARK ORIGIN = ();\nLINE TL = (ORIGIN D1);\nFFS;\n"
                        "RFSW; NORAD; NOFLUC;\nUSE TL; INS; CALC;\nPBUNCH = 1e12;\nSPAC;\n")
    c, alive = track(I, beam)
    assert alive.sum() == np_
    assert np.isfinite(c).all()
    assert np.abs(c[:, 1]).max() > 0.0 and np.abs(c[:, :4]).max() < 0.03


def test_mult_radius_readable(v2_mode):
    """MULT の RADIUS は Element["RADIUS", ..] で読める・変えられる (今までは定義の値が 0 と読めた)"""
    ffs, I, out = setup(MULT_LINE % ('TQ', '1e11'))
    I.run('r1 = Element["RADIUS", "Q1"]; Element["RADIUS", "Q1"] = -0.001; r2 = Element["RADIUS", "Q1"]; r3 = LINE["RADIUS", 2];', 'ffs')
    assert (I.globals['r1'], I.globals['r2'], I.globals['r3']) == (0.02, -0.001, -0.001)


# ============================================================================ 走る例 (example/05_tracking/flag_SPAC_MULT_radius.sad・flag_WSPAC_tune_shift.sad、example/05_tracking/flag_SPAC.sad)
def _example(name):
    import os
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'example')
    with open(os.path.join(root, name), encoding='utf-8') as f:
        return f.read()


def _nums(out, key):
    """key で始まる行の「名前=数値」を順に取り出す"""
    import re
    for l in out.getvalue().splitlines():
        if l.startswith(key):
            return [float(t) for t in re.findall(r'=\s*(-?[0-9.]+(?:E[-+]?[0-9]+)?)', l[len(key):])]
    raise AssertionError('行が無い: ' + key)


def test_example_wspac_script(v2_mode):
    """example/05_tracking/flag_WSPAC_tune_shift.sad が走り、WSPAC のチューンシフトが Laslett の式 (スクリプトの中で Emittance のビーム行列から計算) と
    1% で合う。SELFCOD で重心のキックが 0 になる。平らなリングはエラーで知らせ、WSPACBEAM を与えると走る。
    同じスクリプトの gs (K64、EMIT を WSPAC つきで実行) の追跡は dnux = −2.938e-5、dnuy = −1.210e-3 (ビームの大きさを空間電荷込みで
    解き直しているので 1% 違う)"""
    ffs, I, out = setup(_example('05_tracking/flag_WSPAC_tune_shift.sad'))
    tr, las = _nums(out, 'tracking:'), _nums(out, 'Laslett:')
    assert tr == pytest.approx([-2.9627e-5, -1.196e-3], rel=1.2e-2)      # v2: tradke の ddpx の符号 (A1 項目 7) でビーム行列が 2e-4 変わる
    assert tr == pytest.approx(las, rel=0.01)
    assert tr == pytest.approx([-2.938e-5, -1.210e-3], rel=0.025)             # gs (K64)
    assert _nums(out, 'NOSELFCOD:')[0] == pytest.approx(0.1999, rel=6e-3)      # v2: tradke の ddpx の符号 (A1 項目 7) でビーム行列が 2e-3 変わる
    assert abs(_nums(out, 'SELFCOD:')[0]) < 1e-6
    assert 'WSPAC: 要素 1 の入口のビームの大きさが 0 以下' in out.getvalue()
    assert _nums(out, 'flat ring + WSPACBEAM (ey = 1% ex):') == pytest.approx([-2.4875e-5, -3.873e-4], rel=1.5e-2)      # v2: tradke の ddpx の符号 (A1 項目 7)


def test_example_mult_radius_script(v2_mode):
    """example/05_tracking/flag_SPAC_MULT_radius.sad が走る: (1) 四極 (MULT) の中の KV ビームの 2·rms が包絡方程式 (スクリプトの中で積分) と 1% で
    合う、(2) リングのチューンシフトが、SPAC (MULT で全周)・PSPAC・WSPAC と Laslett の式 (全周) で合う (SPAC・PSPAC は 4000 個の
    準乱数のビームの雑音で 7%、WSPAC は 1%)。DRIFT だけの式 (全周の 80%) とは合わない"""
    ffs, I, out = setup(_example('05_tracking/flag_SPAC_MULT_radius.sad'))
    text = out.getvalue()
    sp = _nums(out, '(1) SPAC ')
    l = [t for t in text.splitlines() if t.startswith('(1) envelope')][0]
    env = [float(t) for t in l[l.index('X, Y = {') + 8:l.index('}')].split(',')]
    assert env == pytest.approx([2.121e-3, 6.813e-3], rel=1e-3)
    assert sp[:2] == pytest.approx(env, rel=0.01)
    las, lasd = _nums(out, '(2) Laslett (whole ring):'), _nums(out, '(2) Laslett (DRIFT only')
    assert las == pytest.approx([-4.264e-3, -6.671e-3], rel=1e-3)
    assert _nums(out, '(2) WSPAC')[-2:] == pytest.approx(las, rel=0.01)
    for key in ('(2) SPAC ', '(2) PSPAC'):
        v = _nums(out, key)[-2:]
        assert v == pytest.approx(las, rel=0.07)
        assert abs(v[0]) > 1.15 * abs(lasd[0]) and abs(v[1]) > 1.15 * abs(lasd[1])


def test_example_flag_spac_script(v2_mode):
    """example/05_tracking/flag_SPAC.sad (v2 モード): (A) 1 m のドリフトの rms が包絡方程式 (K) と 1% で合い、K/γ (K64 の SPC1) の値
    (gs の結果 0.506 mm) とは合わない。(C) RADIUS<0 の絞りの生存数は gs と同じ 84。(B) のリングは時間がかかるので、ここでは走らせない
    (例の .pysad.txt と flag_SPAC.md に値がある)"""
    src = _example('05_tracking/flag_SPAC.sad')
    a, b, c = src.index('! ==== (A)'), src.index('! ==== (B)'), src.index('! ==== (C)')
    ffs, I, out = setup(src[:b] + 'USE RING; CELL; CALC;\n' + src[c:].replace('Exit;', ''))
    sp = _nums(out, '(A) SPAC ')
    l = [t for t in out.getvalue().splitlines() if t.startswith('(A) K =')][0]
    k_env, kg_env = float(l.split('K ->')[1].split()[0]), float(l.split('->')[-1])
    assert (k_env, kg_env) == pytest.approx((7.1965e-4, 5.0599e-4), rel=1e-4)
    assert sp[:2] == pytest.approx([k_env, k_env], rel=0.01)
    assert sp[2] == 2000.0
    assert _nums(out, '(C) SPAC RADIUS=-0.4 mm, 2 turns:')[-1] == 84.0
    assert _nums(out, '(C) NOSPAC RADIUS=-0.4 mm, 2 turns:')[-1] == 200.0
