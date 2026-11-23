# -*- coding: utf-8 -*-
"""v2.4.1 C: MULT (ANGLE=0) の端の縦方向の収束 (docs/30 B292)。

v2 の MULT の端 (AE1・AE2) は、水平方向の運動を厳密な硬い端の写像 (本体の場を面まで通す流れ) で扱い、縦方向の収束は K64 の薄い蹴りの縦成分
Δpy = −tan e·Im((Σ K_n z^n/n!) z)/L を足していた (B224)。厳密な Maxwell 場と比べると 2 つの誤りがあった。
 (1) 面を横切る粒子の傾き θ = atan(px/pz) への依存が無い: 厳密な蹴りは tan(e ± θ) に比例する (面を横切る薄い層の力の積分。入口 +、出口 −)。
     θ の tan θ の分は本体の双極の項 (tsolqu) が既に持つので、足す分は tan(e ± θ) − tan(±θ)。K0 = 0.06・AE = 0.3 で py/y が 1.9% 足りなかった (θ = −0.059)。
 (2) K1 以上の項は Im(z B(z)) (= K64 の cx) でなく、横磁場のスカラーポテンシャル Im ψ (ψ = Σ K_n z^(n+1)/(n+1)!/L) を使う: 水平方向の運動を厳密な流れが担うので
     「面までの余分な経路の寄与」を二重に数えていた (K1 の項が 2 倍)。K1 = 0.4 で py が 4e-7 ずれていた。
独立な基準: 面が z = x tan e の滑らかな階段 f(u) (u = (z − x tan e)/w、5 次の smoothstep) の Maxwell 場 B = ∇(ψ(x, y)·f) (y の 1 次まで厳密。ψ = (K0 y + K1 x y)/L)
の中の運動 (ローレンツ力の常微分方程式、DOP853、rtol 1e-12)。w = 2e-4・1e-4 の線形外挿 (x は w に 1 次で依存する)。ANGLE=0 の BEND の端の鎖 (tbedge、
厳密な硬い端) もこの基準と px 1e-11・py/y 5e-8 で一致する (K0 だけのとき) ので、基準自体の確認になっている。
k64 は K64 の薄い蹴り (追跡) のまま。
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat
from pysad import lattice as lat
from pysad import track6 as t6
from . import _p3_indep as X


def _smooth(s):
    s = min(1.0, max(0.0, s))
    return 10 * s**3 - 15 * s**4 + 6 * s**5, (30 * s**2 - 60 * s**3 + 30 * s**4 if 0.0 < s < 1.0 else 0.0)


def _ref_w(c0, L, K0, K1, e1, e2, w, delta, m=0.03):
    pr = 1.0 + delta
    x, px, y, py = c0
    pz = math.sqrt(pr * pr - px * px - py * py)
    t1, t2 = math.tan(e1), math.tan(e2)
    r0 = np.array([x - px / pz * m, y - py / pz * m, -m, px, py, pz])
    b0, b1 = K0 / L, K1 / L

    def rhs(tau, v):
        X_, Y, Z, PX, PY, PZ = v
        f1, d1 = _smooth((Z - X_ * t1) / w)
        f2, d2 = _smooth((Z - L + X_ * t2) / w)
        g = f1 * (1 - f2)
        gz = d1 / w * (1 - f2) - f1 * d2 / w
        gx = (-t1 * d1 / w) * (1 - f2) - f1 * (t2 * d2 / w)
        beta = b0 + b1 * X_
        Bx = Y * (b1 * g + beta * gx)
        By = beta * g
        Bz = Y * beta * gz
        vx, vy, vz = PX / pr, PY / pr, PZ / pr
        return [vx, vy, vz, vy * Bz - vz * By, vz * Bx - vx * Bz, vx * By - vy * Bx]

    def ev(tau, v):
        return v[2] - (L + m)
    ev.terminal = True
    ev.direction = 1
    sol = solve_ivp(rhs, (0, 1.0), r0, method='DOP853', rtol=1e-12, atol=1e-15, events=ev)
    X_, Y, Z, PX, PY, PZ = sol.y_events[0][0]
    dz = L - Z
    return np.array([X_ + PX / PZ * dz, PX, Y + PY / PZ * dz, PY])


def _ref(c0, L, K0, K1, e1, e2, delta=0.0):
    return 2 * _ref_w(c0, L, K0, K1, e1, e2, 1e-4, delta) - _ref_w(c0, L, K0, K1, e1, e2, 2e-4, delta)


def _elem(src):
    ffs, _ = X.run('MOMENTUM = 3GEV;\n%s\nMARK ST=();\nLINE ONE=(ST SA);\nFFS;\nUSE ONE;\n' % src)
    return [e for e in ffs.lat.elements if e.name == 'SA'][0], ffs


def _calc(mode, c0, L, K0, K1, e1, e2, delta):
    old = compat._mode
    compat.set_mode(mode)
    try:
        el, _f = _elem('MULT SA = (L=%r K0=%r K1=%r AE1=%r AE2=%r);' % (L, K0, K1, e1, e2))
        return lat._el_cod_step(el, np.array(list(c0) + [0.0, delta]), delta)[1][:4]
    finally:
        compat.set_mode(old)


def _track(mode, c0, L, K0, K1, e1, e2, delta):
    old = compat._mode
    compat.set_mode(mode)
    try:
        el, ffs = _elem('MULT SA = (L=%r K0=%r K1=%r AE1=%r AE2=%r);' % (L, K0, K1, e1, e2))
        c = np.zeros((1, 6))
        c[0, :4] = c0
        c[0, 5] = delta
        return np.asarray(t6.track_particles(ffs.lat.elements, c, 1, ffs.momentum(), 0.0, False, False, False))[0][:4]
    finally:
        compat.set_mode(old)


_C0 = (2e-3, 1e-3, 1e-4, 5e-4)
CASES = [(0.2, 0.06, 0.0, 0.3, 0.3), (0.2, 0.06, 0.0, 0.3, -0.2), (0.2, 0.01, 0.0, 0.3, 0.3),
         (0.2, 0.0, 0.4, 0.3, 0.3), (0.2, 0.03, 0.4, 0.3, -0.2)]


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('delta', [0.0, 0.01])
def test_mult_edge_matches_maxwell_reference_v2(case, delta):
    """v2: CALC の要素写像・追跡 (track6) とも独立な Maxwell 場の基準と x・px・y・py が 5e-8 以内 (基準の w 外挿の誤差と同程度。入口 1e-3 の px・0.5 mrad の py、
    y = 0.1 mm の粒子。δ = 0.01 を含む)。追跡と CALC は互いに 1e-11 以内"""
    L, K0, K1, e1, e2 = case
    r = _ref(_C0, L, K0, K1, e1, e2, delta)
    a = _calc('v2', _C0, L, K0, K1, e1, e2, delta)
    b = _track('v2', _C0, L, K0, K1, e1, e2, delta)
    assert np.abs(a - r).max() < 5e-8, (a - r)
    assert np.abs(b - r).max() < 5e-8, (b - r)
    assert np.abs(a - b).max() < 1e-11


def test_k64_mult_edge_keeps_the_thin_kick():
    """k64 は K64 の薄い蹴りのまま (K0 = 0.06・AE = 0.3 の追跡は py が基準と 5e-7、px が 6e-6 ずれる (CALC は 5e-6・4e-4)。v2 は 2e-8 以内)"""
    L, K0, K1, e1, e2 = CASES[0]
    r = _ref(_C0, L, K0, K1, e1, e2)
    k = _track('k64', _C0, L, K0, K1, e1, e2, 0.0)
    assert np.abs(k - r)[3] > 3e-7 and np.abs(k - r)[1] > 3e-6


def test_vertical_focusing_of_uniform_edge_includes_crossing_angle():
    """K0 だけの MULT (AE2 = 0.3、出口の端だけ): py/y は −(K0/L)·tan(e − θ) (θ = 出口の面の内側の傾き ≈ −K0)。K64 の薄い蹴り (θ 依存なし) は −(K0/L)·tan e + 本体の −(K0/L)·tan(−θ)
    を足した値で、厳密な値より 1.9% 小さかった。v2 は基準 (と ANGLE=0 の BEND の端の鎖) と 3e-4 の比 (y = 1e-5 では py が 1e-6 で、基準の w 外挿の誤差) で一致"""
    c0 = (2e-3, 1e-3, 1e-5, 0.0)
    r = _ref(c0, 0.2, 0.06, 0.0, 0.0, 0.3)
    v2 = _track('v2', c0, 0.2, 0.06, 0.0, 0.0, 0.3, 0.0)
    assert abs(v2[3] / r[3] - 1.0) < 1e-3, (v2[3], r[3])
    old = compat._mode
    compat.set_mode('v2')
    try:
        bend = X.one_element_track('BEND SA = (L=0.2 ANGLE=0 K0=0.06 AE2=0.3);', np.array([c0[0], c0[1], c0[2], c0[3], 0.0, 0.0]))
    finally:
        compat.set_mode(old)
    assert abs(bend[3] / r[3] - 1.0) < 1e-3, (bend[3], r[3])      # 基準の確認 (厳密な BEND の鎖)
