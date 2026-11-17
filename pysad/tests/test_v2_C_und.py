# -*- coding: utf-8 -*-
"""v2 の修正 B37: UND の磁場 (FBX・FBY [T])。K64 はキーワードが登録されず磁場なしのドリフトになる。

根拠 (解析式。平面アンジュレータ、KX=QY=0、軸上の磁場 By = F0 cos(kz s)、F0 = FBX·c/p [1/m]):
  - 軸上の軌道: x' = -(F0/kz) sin(kz s)。整数周期の出口で x = x' = 0。
  - 行路長: Δz = -∫x'²/2 ds = -(F0/kz)² L/4。
  - 垂直の収束: y'' = -(F0²/2) y/(1+δ)² (周期平均)。水平は収束なし。KX≠0 では kx = -F0² KX²/(2 kz²)、ky = F0² (kz²+KX²)/(2 kz²)。
  - シンプレクティック性 (正準変数 x, px, y, py の Jacobi 行列)。
  - 放射 (v2 の追加): Δδ = -(2/3) r_e γ³ ∫(1/ρ²) ds = -(2/3) r_e γ³ F0² L/2。
K64 の undulator.f は陰的な刻みの右辺に dpds を二重に掛けていて、収束力が刻み幅に比例して消える (v2 で修正)。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat, emit6d as _e, und
from pysad.elements import Undulator
from pysad.track6n import make_state

P = 3.0e9
RE = 2.8179403e-15


def _run(x, sl=20000, rad=False, **kw):
    ph = _e._Phys(P, 100.0)
    par = dict(FBX=1.0, POLE=20, SLICE=sl)
    par.update(kw)
    el = Undulator(name='U', L=1.0, params=par)
    pr = und.compile_und(el, ph.pgev)
    S = make_state(np.array([x], float), ph)
    rs = np.array([1, 1], dtype=np.int64)
    und.k_und_field(S, 1, *pr, ph.p0, ph.h0, rad, False, rs, ph.anrad, ph.cuc)
    return S[0, :6].copy(), pr, ph


def test_b37_k64_rejects_field(k64_mode):
    ph = _e._Phys(P, 100.0)
    with pytest.raises(NotImplementedError):
        und.compile_und(Undulator(name='U', L=1.0, params=dict(FBX=1.0, POLE=20, SLICE=100)), ph.pgev)


def test_b37_on_axis_orbit_and_path_length(v2_mode):
    o, pr, _ = _run([0, 0, 0, 0, 0, 0])
    f0, kz = pr[1], pr[4]
    assert abs(f0 - 1.0 * 299792458.0 / P) < 1e-15
    assert abs(o[0]) < 1e-12 and abs(o[1]) < 1e-12
    assert abs(o[4] / (-(f0 / kz)**2 / 4.0) - 1.0) < 1e-6


def test_b37_vertical_focusing(v2_mode):
    y0 = 1.0e-4
    o, pr, _ = _run([0, 0, y0, 0, 0, 0])
    r = pr[1] / math.sqrt(2.0)
    assert abs(o[2] - y0 * math.cos(r)) < 2e-6 * y0
    assert abs(o[3] / (-y0 * r * math.sin(r)) - 1.0) < 5e-4
    o, _, _ = _run([0, 0, 0, 1e-5, 0, 0])
    assert abs(o[2] / (1e-5 * math.sin(r) / r) - 1.0) < 1e-5
    # 運動量の依存: 収束力は 1/(1+δ)²
    o, _, _ = _run([0, 0, y0, 0, 0, 1e-2])
    assert abs(o[3] / (-y0 * r * r / 1.01**2) - 1.0) < 2e-3


def test_b37_kx_focusing(v2_mode):
    x0 = y0 = 1.0e-4
    o, pr, _ = _run([x0, 0, y0, 0, 0, 0], KX=30.0)
    f0, kz, kx, ky = pr[1], pr[4], pr[5], pr[6]
    xm = x0 - f0 / kz**2                       # 蛇行の中心 (x = x0 + (F0/kz²)(cos kz s - 1))
    assert abs(o[1] / (f0**2 * kx**2 / (2 * kz**2) * xm) - 1.0) < 5e-3     # 水平は発散
    assert abs(o[3] / (-f0**2 * ky**2 / (2 * kz**2) * y0) - 1.0) < 5e-3


def test_b37_symplectic(v2_mode):
    x0 = np.array([2e-4, 1e-4, -1e-4, 5e-5, 0.0, 0.0])
    h = 1e-6
    J = np.zeros((4, 4))
    for k in range(4):
        d = np.zeros(6)
        d[k] = h
        a, _, _ = _run(x0 + d, sl=2000, KX=20.0, QY=10.0, FBY=0.5, PHI=90.0)
        b, _, _ = _run(x0 - d, sl=2000, KX=20.0, QY=10.0, FBY=0.5, PHI=90.0)
        J[:, k] = (a[:4] - b[:4]) / (2 * h)
    Sm = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], float)
    assert np.abs(J.T @ Sm @ J - Sm).max() < 1e-7


def test_b37_radiation_loss(v2_mode):
    o, pr, ph = _run([0, 0, 0, 0, 0, 0], rad=True)
    exp = -(2.0 / 3.0) * RE * ph.h0**3 * pr[1]**2 / 2.0
    assert abs(o[5] / exp - 1.0) < 1e-4
    assert abs(o[0]) < 1e-9 and abs(o[1]) < 1e-10


SRC = r'''MOMENTUM = 3GEV;
DRIFT DR = (L=1);
UND U1 = (L=1 FBX=1 POLE=20 SLICE=4000);
MARK ST = ();
LINE BL = (ST DR U1 DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
r = TrackParticles[{1, {{0,0},{0,0},{1e-4,0},{0,0},{0,0},{0,0},{1,1}}}, 4, 1, 1];
$FORM = "18.12";
Print["PY ", r[[2,4,1]], " ", r[[2,1,2]]];
Exit;
'''


def test_b37_script_fbx(v2_mode):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('PY')][0].split()
    f0 = 299792458.0 / P
    assert abs(float(line[1]) / (-1e-4 * f0**2 / 2.0) - 1.0) < 2e-3
    assert abs(float(line[2])) < 1e-10


# ---------------------------------------------------------------------------- リングでの光学・エミッタンス (v2)
RING = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
%s
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 8*CELL1 W RFC);
FFS;
USE RING;
CELL; CALC;
$FORM="20.14";
PageWidth=1999;
Print["NU ",Twiss["NX","$$$"]/2/Pi," ",Twiss["NY","$$$"]/2/Pi," ",Twiss["BY","W"]," ",Twiss["AY","W"]];
e=Emittance[];
Print["EM ",(Emittances/.e)[[1]]," ",EnergyLossU0/.e," ",(DampingRate/.e)[[1]]," ",(DampingRate/.e)[[2]]," ",(DampingRate/.e)[[3]]];
Exit;
'''
F0R = 1.8 * 299792458.0 / 3.0e9


def _ring(defn):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(RING % defn, out=buf)
    out = {}
    for l in buf.getvalue().splitlines():
        t = l.split()
        if t and t[0] in ('NU', 'EM'):
            out[t[0]] = [float(v) for v in t[1:] if v != '{}']
    return out


def test_b37_ring_und_without_field_is_a_drift(v2_mode):
    """v2: 光学・Emittance[] でも UND は追跡と同じ写像 (磁場なしなら長さ L のドリフト)。K64 は光学では恒等"""
    a = _ring('DRIFT W = (L=1);')
    b = _ring('UND W = (L=1 POLE=20 SLICE=2000);')
    for k in ('NU', 'EM'):
        for x, y in zip(a[k], b[k]):
            assert abs(x - y) <= 1e-6 * abs(x)


def test_b37_ring_tune_shift_energy_loss_and_emittance(v2_mode):
    d = _ring('DRIFT W = (L=1);')
    u = _ring('UND W = (L=1 FBX=1.8 POLE=20 SLICE=2000);')
    # 垂直チューンのずれ: Δν = (1/4π) ∫β k ds、k = F0²/2、ドリフトの中の β(s) = β0 - 2 α0 s + γ0 s²
    b0, a0 = d['NU'][2], d['NU'][3]
    ib = b0 - a0 + (1 + a0 * a0) / b0 / 3.0
    assert abs((u['NU'][1] - d['NU'][1]) / (F0R**2 / 2 * ib / (4 * math.pi)) - 1.0) < 1e-2
    assert abs(u['NU'][0] - d['NU'][0]) < 1e-6                       # 水平は収束なし
    # 放射損失: ΔU = (2/3) r_e γ³ (F0²/2) L · E
    gam = 3.0e9 / 510998.95
    du = (2.0 / 3.0) * RE * gam**3 * F0R**2 / 2 * 3.0e9
    assert abs((u['EM'][1] - d['EM'][1]) / du - 1.0) < 2e-3
    # Robinson の和則: 減衰率の和 = 2 U0/E
    assert abs(-sum(u['EM'][2:5]) / (2 * u['EM'][1] / 3.0e9) - 1.0) < 1e-3
    # 量子励起: 同じ <1/ρ²> の矩形磁場のウィグラー (BEND の列、pysad の BEND の放射) と比べる。
    # 励起の増分 r = (ε/ε0)(αx/αx0) - 1 は <|1/ρ|³> に比例: 正弦波 4/(3π) F0³、矩形 (F0/√2)³ → 比 1.2004
    h = F0R / math.sqrt(2.0)
    a = h * 0.025
    w = _ring('BEND BH=(L=0.0125 ANGLE=%r E1=0 E2=1) BH2=(L=0.0125 ANGLE=%r E1=1 E2=0) BM=(L=0.025 ANGLE=%r E1=.5 E2=.5) '
              'BP=(L=0.025 ANGLE=%r E1=.5 E2=.5);\nLINE W=(BH 19*(BM BP) BM BH2);' % (a / 2, a / 2, -a, a))
    assert abs(w['NU'][1] - u['NU'][1]) < 2e-5                       # 収束も同じ (矩形磁場の列と一致)
    assert abs(w['EM'][1] / u['EM'][1] - 1.0) < 1e-4
    ru = u['EM'][0] / d['EM'][0] * u['EM'][2] / d['EM'][2] - 1.0
    rw = w['EM'][0] / d['EM'][0] * w['EM'][2] / d['EM'][2] - 1.0
    assert abs(ru / rw / ((4 / (3 * math.pi)) / 2**-1.5) - 1.0) < 2e-2
