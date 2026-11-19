# -*- coding: utf-8 -*-
"""v2next 担当 P3 の試験の補助 (独立な基準と、gs と同じ台本の組み立て)。

- 台本の組み立て: 12 セルの FODO (3 GeV、周長 45.6 m) にステアリング (BEND ANGLE=0 K0≠0) か MULT を入れたリングと、大きなリング A
  (pysad/tests/kek-ls_v03_12.lat、3 GeV、603 m) の先頭に補正磁石 2 台を入れたリング。gs の値は _p3_gsdata.py。
- 独立な基準 (pysad の追跡・光学のコードを使わない):
  * `maxwell_steer`: 補正磁石の端 (F1: 磁場が直線で立ち上がる長さ) の磁場 By = B(s) − y²B''/2、Bs = y B'(s) (x 方向に一様な Maxwell 場) の中の
    厳密な運動 (正準方程式の数値積分)。
  * ロビンソンの定理と厳密な幾何 (円弧と直線) による端の角の Jx: test_v2_E_steer._jx_exact。
"""
import io
import math
import os
import re

import numpy as np
from scipy.integrate import solve_ivp

_HERE = os.path.dirname(os.path.abspath(__file__))


def floats(s):
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', s)]


def tag_values(text, tag):
    """text の中の `tag ` で始まる出力 (複数行にまたがる) の先頭 18 個の数値 (3 粒子 x 6 座標)"""
    i = text.index('\n' + tag + ' ') + 1 if ('\n' + tag + ' ') in text else text.index(tag + ' ')
    seg = text[i + len(tag):]
    m = re.search(r'\n(?=[A-Za-z?])', seg)
    seg = seg[:m.start()] if m else seg
    return floats(seg)[:18]


def run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _ = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


# ---------------------------------------------------------------- 台本 (gs と同じ組み立て)
TRACK_TAIL = '''FFS;
$FORM="14.10";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
USE RING; CELL; CALC; NORFSW; NORAD;
r=TrackParticles[{1,bm},1,1,3];
Print["N3 ", r[[2]]];
RFSW; RAD; NOFLUC;
r=TrackParticles[{1,bm},1,1,3];
Print["R3 ", r[[2]]];
Exit;
'''

_STEER = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) DS = (L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND SA = (L=0.2 ANGLE=0 %s);
BEND SB = (L=0.2 ANGLE=0 %s);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELLA = (QF D2 B D2 QD DS SA D1);
LINE CELLB = (QF D2 B D2 QD DS SB D1);
LINE CC = (CELLA CELLB);
LINE RING = (ST 6*CC RFC);
'''

_MULT = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
MULT QF = (L=0.3 K1=0.7 %s);
MULT QD = (L=0.3 K1=-0.6 %s);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
'''

EMIT_TAIL = '''FFS;
$FORM="14.10";
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
RFSW; RAD; NOFLUC;
e = Emittance[];
Print["Tunes ", Tunes/.e];
Print["U0 ", EnergyLossU0/.e, " alpha ", MomentumCompaction/.e, " EqPos ", EquilibriumPosition/.e];
Print["Emittances ", Emittances/.e];
Print["Damp ", DampingRate/.e];
Exit;
'''

# gs の値と比べる追跡の組 (キー = _p3_gsdata.TRACK)。(kind, 要素 1 の引数, 要素 2 の引数)
TRACK_CASES = {
    's1': ('steer', 'K0=0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01', 'K0=-0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01'),
    'b_c': ('steer', 'K0=0.003 F1=0.04 FRINGE=1', 'K0=-0.003 F1=0.04 FRINGE=2'),
    'b_d': ('steer', 'K0=0.003 F1=0.04 FRINGE=-1', 'K0=-0.003 F1=0.04 FRINGE=-2'),
    'b_e': ('steer', 'K0=0.003 F1=0.04 FRINGE=0', 'K0=-0.003 F1=0.04 FRINGE=0'),
    'b_f': ('steer', 'K0=0.003 AE1=0.2 AE2=0.1', 'K0=-0.003 AE1=0.2 AE2=0.1'),
    'b_g': ('steer', 'K0=0.003 E1=0.01 E2=0.02 DX=1e-4 ROTATE=0.1 FB2=0.03 FRINGE=3', 'K0=-0.003 AE2=0.05 FB2=0.03 FRINGE=2 DISFRIN=1'),
    'm_a': ('mult', 'F1=0.05 FRINGE=3', 'F1=0.05 FRINGE=3'),
    'm_b': ('mult', 'F1=0.05 F2=0.02 FRINGE=3', 'F1=0.04 F2=0.03 FRINGE=3'),
    'm_c': ('mult', 'F1=0.05 FRINGE=1', 'F1=0.05 FRINGE=2'),
    'm_d': ('mult', 'F1=0.05 FRINGE=3 K0=0.0004 FB1=0.03 FB2=0.02', 'F1=0.05 FRINGE=3 K0=-0.0004 FB1=0.03 FB2=0.02'),
    'm_e': ('mult', 'K0=0.001 FRINGE=3 FB1=0.05 FB2=0.04 K2=20', 'K0=-0.001 FRINGE=3 FB1=0.05 FB2=0.04 K2=-20'),
    'm_f': ('mult', 'AE1=0.1 AE2=-0.05 K0=0.001', 'AE1=0.1 AE2=-0.05 K0=-0.001'),
    'm_g': ('mult', 'AE1=0.2 AE2=0.1', 'AE1=-0.2 AE2=0.3'),
    'm_h': ('mult', 'FRINGE=2 K2=30', 'FRINGE=1 K2=-30'),
    'm_i': ('mult', 'F1K1F=0.02 F1K1B=0.04 F1=0.03 FRINGE=3 AE1=0.05 K0=0.0006 FB1=0.03',
            'F1K1F=0.02 F1K1B=0.04 F1=0.03 FRINGE=3 AE2=0.05 K0=-0.0006 FB2=0.03 ROTATE=0.1 DX=1e-4'),
}

# Emittance[]・CALC の gs との比較 (キー = _p3_gsdata.EMIT): MULT の端の角 AE1・AE2
EMIT_CASES = {
    'e_a': ('AE1=0.3 AE2=0.3 K0=0.01', 'AE1=0.3 AE2=0.3 K0=-0.01'),
    'e_b': ('AE1=0.3 K0=0.01', 'AE1=0.3 K0=-0.01'),
    'e_c': ('AE2=0.3 K0=0.01', 'AE2=0.3 K0=-0.01'),
    'e_d': ('AE1=0.1 AE2=-0.1 K0=0.001 K2=5 F1=0.04 FRINGE=3', 'AE1=0.1 AE2=-0.1 K0=-0.001 K2=-5 F1=0.04 FRINGE=3'),
    'e_e': ('AE1=0.2 AE2=0.3', 'AE1=0.3 AE2=-0.2'),
}

# 大きなリング A の先頭に入れる補正磁石 2 台 (キー = ra1..ra3)
RING_A_INSERT = {
    'ra1': 'BEND SA = (L=0.2 ANGLE=0 K0=5e-4 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1);\n'
           'BEND SB = (L=0.2 ANGLE=0 K0=-5e-4 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1);',
    'ra2': 'MULT SA = (L=0.2 K1=0.05 K0=5e-4 F1=0.03 F2=0.01 FB1=0.02 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1);\n'
           'MULT SB = (L=0.2 K1=-0.05 K0=-5e-4 F1=0.03 F2=0.01 FB1=0.02 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1);',
    'ra3': 'BEND SA = (L=0.2 ANGLE=0 K0=5e-4 AE1=0.2 AE2=0.3);\n'
           'BEND SB = (L=0.2 ANGLE=0 K0=-5e-4 AE1=-0.3 AE2=0.2);',
}


def ring_script(kind, a, b, tail):
    return (_STEER if kind == 'steer' else _MULT) % (a, b) + tail


def mult_emit_script(a, b):
    return _MULT % (a, b) + EMIT_TAIL


def ring_a_script(key, tail):
    path = os.path.join(_HERE, 'kek-ls_v03_12.lat')
    with open(path, encoding='utf-8', errors='replace') as fh:
        lat = fh.read()
    lat = lat.replace(' RING      =(MKS', ' RING0     =(MKS', 1)
    return lat + '\n' + RING_A_INSERT[key] + '\nMARK ST0 = ();\nLINE RING = (ST0 SA SB RING0);\n' + tail


RING_A_TAIL = '''FFS;
$FORM="14.10";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
NORFSW; NORAD;
r=TrackParticles[{1,bm},1,1,3];
Print["N3 ", r[[2]]];
RFSW; RAD; NOFLUC;
r=TrackParticles[{1,bm},1,1,3];
Print["R3 ", r[[2]]];
e = Emittance[];
Print["Tunes ", Tunes/.e];
Print["U0 ", EnergyLossU0/.e, " alpha ", MomentumCompaction/.e, " EqPos ", EquilibriumPosition/.e];
Print["Emittances ", Emittances/.e];
Print["Damp ", DampingRate/.e];
Exit;
'''


# ---------------------------------------------------------------- 追跡のチューン (閉軌道のまわりの 4 次元の一周写像) と CALC のチューン
def tracking_and_calc_tunes(src, fast=True):
    """src (FFS; USE RING; CELL; CALC; と Print["CALC tunes ", ...] を含む台本) の追跡 (RF・放射なし) の閉軌道のまわりの一周写像の固有値から
    チューンの小数部を求める (ニュートン法、数値微分)。戻り値: (追跡の x・y のチューンの小数部の組 (昇順)、CALC の x・y のチューンの小数部)"""
    from pysad import track6, track6n
    ffs, out = run(src)
    els = ffs.lat.elements
    p = ffs.momentum()
    tr = track6n.track_particles_fast if fast else track6.track_particles

    def f(v):
        c = np.zeros(6)
        c[:4] = v
        return np.asarray(tr(els, np.array([c]), 1, p, 0.0, False, False, False))[0][:4]
    h = 1e-6

    def jac(v):
        J = np.zeros((4, 4))
        for i in range(4):
            d = np.zeros(4)
            d[i] = h
            J[:, i] = (f(v + d) - f(v - d)) / (2 * h)
        return J
    v = np.zeros(4)
    for _ in range(8):
        v = v - np.linalg.solve(jac(v) - np.eye(4), f(v) - v)
    ev = np.linalg.eigvals(jac(v))
    nu = sorted(set(round(abs(np.angle(z)) / (2 * math.pi), 9) for z in ev))
    m = re.search(r'CALC tunes\s+([-\d.eE+]+)\s+([-\d.eE+]+)', out)
    al = lambda x: min(x % 1.0, 1.0 - x % 1.0)
    return nu, sorted([al(float(m.group(1))), al(float(m.group(2)))])


# ---------------------------------------------------------------- 独立な基準: x 方向に一様な磁場の端 (Maxwell の ramp) の厳密な運動
def maxwell_steer(fb1, fb2, al, K0, y0, delta=0.0):
    """補正磁石 (長さ al、曲げ K0 [rad]、B0 = K0/al) の端の磁場の厳密な運動。磁場は x に依らず、B(s) は入口 [−fb1/2, fb1/2]・出口 [al−fb2/2, al+fb2/2]
    で 0 → B0 → 0 に直線で変わる (積分は B0·al で、F1 を使うときの K64 の模型と同じ)。By = B − y²B''/2、Bs = y B' (∇·B = 0、∇×B = 0)。
    正準方程式 (ベクトルポテンシャル Ax = F(s) − y²B'/2、F' = B) を solve_ivp (DOP853、rtol 1e-13) で積分する。運動量は p0 で規格化。
    s = −fb1/2 で (x, px, y, py) = (0, 0, y0, 0) から s = al + fb2/2 まで。戻り値: (x, px, y, py)、運動量は運動学的な値"""
    B0 = K0 / al
    pts = [-fb1 / 2, fb1 / 2, al - fb2 / 2, al + fb2 / 2]
    segs = [(pts[0], pts[1], 0.0, B0 / fb1), (pts[1], pts[2], B0, 0.0), (pts[2], pts[3], B0, -B0 / fb2)]
    F = 0.0
    state = np.array([0.0, 0.0, y0, 0.0])
    for (sa, sb, b0, bp) in segs:
        def rhs(s, v, F0=F, sa=sa, b0=b0, bp=bp):
            x, px, y, py = v
            Fs = F0 + b0 * (s - sa) + 0.5 * bp * (s - sa)**2
            pi = px - Fs + 0.5 * y * y * bp
            pz = math.sqrt((1 + delta)**2 - pi * pi - py * py)
            return [pi / pz, 0.0, py / pz, -(pi / pz) * y * bp]
        sol = solve_ivp(rhs, (sa, sb), state, method='DOP853', rtol=1e-13, atol=1e-20)
        state = sol.y[:, -1]
        F += b0 * (sb - sa) + 0.5 * bp * (sb - sa)**2
    x, px, y, py = state
    return x, px - F, y, py


def one_element_track(el_def, coords, rad=False, fast=False):
    """1 要素 (MARK のあとの el_def) だけの追跡 (RF なし)。coords: 6 成分。戻り値: 6 成分"""
    from pysad import track6, track6n
    src = 'MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n' % el_def
    ffs, _ = run(src)
    c = np.zeros((1, 6))
    c[0, :] = coords
    f = track6n.track_particles_fast if fast else track6.track_particles
    return np.asarray(f(ffs.lat.elements, c, 1, ffs.momentum(), 0.0, False, rad, False))[0]


# ---------------------------------------------------------------- 独立な基準: 勾配つきの極 (K1≠0) の厳密な幾何 + ロビンソンの定理
def pole_k1(x, phi, d, L, K0, K1, a1, a2):
    """硬い端の磁場 B(x) = (K0 + K1 x)/L の極 (K0 > 0 は −x へ曲がる、K1 > 0 は x に収束)。端の面は実験室系で z = x tan a1 (入口)、z = L − x tan a2 (出口)。
    参照軌道に垂直な面 z = 0 の (x, φ) から出口の面 z = L の (x, φ) と ∫B² dl を返す (運動の常微分方程式を数値積分)"""
    ta1, ta2 = math.tan(a1), math.tan(a2)
    t = x * ta1 / (math.cos(phi) - math.sin(phi) * ta1)
    xp = x + t * math.sin(phi)
    zp = t * math.cos(phi)

    def rhs(l, v):
        xx, _zz, ph, _W = v
        B = (K0 + K1 * xx) / L
        return [math.sin(ph), math.cos(ph), -B / (1 + d), B * B]

    def ev(l, v):
        return v[1] + v[0] * ta2 - L
    ev.terminal = True
    ev.direction = 1
    sol = solve_ivp(rhs, (0.0, 2.0), [xp, zp, phi, 0.0], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    xx, zz, ph, W = sol.y_events[0][0]
    tt = (L - zz) / math.cos(ph)
    return xx + tt * math.sin(ph), ph, W


def jx_exact_k1(K0, K1, a1, a2, LP=0.2, DW=0.05):
    """test_v2_E_steer._jx_exact と同じ 24 セルの FODO の極 3 個の組 (K0・−2K0・K0) に勾配 K1 (積分値、3 個とも同じ符号と大きさ) を足した Jx = 3 − Je、
    Je = 2 + dln(Σ∫B² dl)/dδ (閉軌道は δ ごとの 1 セルの写像の不動点)。K1 = 0 で _jx_exact と 6e-12 で一致"""
    from . import test_v2_E_steer as T
    wB = 1 / T._RHO**2

    def cell(w, d):
        x, phi = T._quad(w[0], w[1], d, 0.3, 0.7)
        x += 0.5 * math.tan(phi)
        x, phi, lb = T._bend(x, phi, d)
        wsum = wB * lb
        x += 0.5 * math.tan(phi)
        x, phi = T._quad(x, phi, d, 0.3, -0.6)
        x += 0.1 * math.tan(phi)
        for i, k in enumerate((K0, -2 * K0, K0)):
            x, phi, W = pole_k1(x, phi, d, LP, k, K1, a1, a2)
            wsum += W
            if i < 2:
                x += DW * math.tan(phi)
        x += 1.1 * math.tan(phi)
        return np.array([x, phi]), wsum

    def wtot(d):
        w = np.zeros(2)
        for _ in range(30):
            r = cell(w, d)[0] - w
            J = np.zeros((2, 2))
            for i in range(2):
                dw = np.zeros(2)
                dw[i] = 1e-7
                J[:, i] = (cell(w + dw, d)[0] - cell(w - dw, d)[0]) / 2e-7
            w = w - np.linalg.solve(J - np.eye(2), r)
        return cell(w, d)[1]
    h = 1e-5
    w0 = wtot(0.0)
    return 3.0 - (2.0 + (wtot(h) - wtot(-h)) / (2 * h) / w0)
