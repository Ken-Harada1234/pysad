#!/usr/bin/env python3
"""
bench_level2_pysad.py
=====================
Level 2 ベンチマーク — pysad 側スクリプト

対応する SAD スクリプト: bench_level2.sad

テスト内容:
  Test A: オフモメンタム Twiss (delta = -0.02, -0.01, 0, +0.01, +0.02)
  Test B: 複合偏向電磁石 (BEND K1 ≠ 0): 放射積分 I1..I5, 自然エミッタンス

K1 単位系の確認 (SAD ↔ pysad):
  SAD   : BEND K1 = 0.08 [m^-1] (integrated over full length L=0.8 m)
  pysad : Bend K1 = 0.10 [m^-2] (per-unit-length; parser が integrated/L を格納)
  → pysad 内 el.K1 = 0.10 m^-2, SAD integrated = el.K1 × L = 0.08 m^-1

事前条件:
  lattice.py の _radiation_integrals_bend の K1 バグ修正済み:
    k = K1      (修正後; K1 は既に per-unit-length [m^-2])
    DKL1 = K1 * DL  (修正後; INTGRL.f: DKL1=RKL1/IRP=K1_per_length×DL)

実行:
  cd "Claude outputs"
  python bench_level2_pysad.py
"""

import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)

from pysad.elements import make_element
from pysad.lattice import Lattice, _transport_6x6

# ---------------------------------------------------------------------------
# ラティスパラメータ (bench_level1 と同一基本セル)
# ---------------------------------------------------------------------------

N_CELLS   = 20           # セル数
ENERGY_EV = 2.5e9        # ビームエネルギー [eV]

L_QF  = 0.3
L_QD  = 0.3
L_D   = 1.2
L_B   = 0.8

K1F   =  1.8             # QF 四極強度 [m^-2] (per-unit-length)
K1D   = -1.8             # QD 四極強度 [m^-2]

B_ANGLE = np.pi / 20     # 偏向角 [rad]
B_E1    = B_ANGLE / 2
B_E2    = B_ANGLE / 2

# Test B: 複合偏向電磁石
# SAD  K1_integrated = 0.08 m^-1
# pysad K1_per_length = 0.08 / 0.8 = 0.10 m^-2
BK_K1_SAD   = 0.08       # SAD 積分 K1 [m^-1]
BK_K1_PYSAD = BK_K1_SAD / L_B   # pysad per-unit-length [m^-2] = 0.10

L_CELL = L_QF + 2*L_D + L_B + L_QD + 2*L_D + L_B
CIRCUMFERENCE = N_CELLS * L_CELL

print(f"ラティス: {N_CELLS} FODO セル, C = {CIRCUMFERENCE:.3f} m, E = {ENERGY_EV*1e-9:.2f} GeV")
print(f"BK: K1_SAD = {BK_K1_SAD} m^-1 (integrated), K1_pysad = {BK_K1_PYSAD} m^-2 (per-length)")
print()

# ---------------------------------------------------------------------------
# 補助関数 (Level 1 と同一)
# ---------------------------------------------------------------------------

def periodic_dispersion(elements, delta=0.0):
    """一周 6×6 行列から周期閉分散 (etax, etapx) を計算。"""
    M = _transport_6x6(elements, delta)
    Mx = M[:2, :2]
    d  = M[:2, 5]
    try:
        eta = np.linalg.solve(np.eye(2) - Mx, d)
    except np.linalg.LinAlgError:
        eta = np.array([0.0, 0.0])
    return eta[0], eta[1]

def compute_phase_table(elements, twiss0_dict, eta0, etap0, delta=0.0):
    """全要素出口で s, betx, bety, alfx, alfy, phix, phiy, etax を計算。"""
    ax = twiss0_dict['ax']; bx = twiss0_dict['bx']
    ay = twiss0_dict['ay']; by = twiss0_dict['by']
    eta = eta0; etap = etap0
    phix = 0.0; phiy = 0.0
    s = 0.0

    rows = [dict(s=s, betx=bx, alfx=ax, bety=by, alfy=ay,
                 phix=phix, phiy=phiy, etax=eta, etapx=etap)]

    gx = (1 + ax**2) / bx
    gy = (1 + ay**2) / by

    for el in elements:
        M = el.transfer_matrix(delta)
        a11 = M[0,0]; a12 = M[0,1]; a21 = M[1,0]; a22 = M[1,1]
        b11 = M[2,2]; b12 = M[2,3]; b21 = M[3,2]; b22 = M[3,3]

        bx_new = a11**2*bx - 2*a11*a12*ax + a12**2*gx
        ax_new = -(a11*a21*bx) + (a11*a22+a12*a21)*ax - a12*a22*gx
        dpx = np.arctan2(a12, a11*bx - a12*ax)

        by_new = b11**2*by - 2*b11*b12*ay + b12**2*gy
        ay_new = -(b11*b21*by) + (b11*b22+b12*b21)*ay - b12*b22*gy
        dpy = np.arctan2(b12, b11*by - b12*ay)

        new_eta  = M[0,0]*eta + M[0,1]*etap + M[0,5]
        new_etap = M[1,0]*eta + M[1,1]*etap + M[1,5]

        gx = (1 + ax_new**2) / bx_new
        gy = (1 + ay_new**2) / by_new
        bx = bx_new; ax = ax_new
        by = by_new; ay = ay_new
        eta = new_eta; etap = new_etap
        phix += dpx; phiy += dpy
        s += el.L

        rows.append(dict(s=s, betx=bx, alfx=ax, bety=by, alfy=ay,
                         phix=phix, phiy=phiy, etax=eta, etapx=etap))
    return rows

# ---------------------------------------------------------------------------
# ラティス構築 (Test A: 純偏向 BEND)
# ---------------------------------------------------------------------------

def make_cell_pure(cell_idx):
    return [
        make_element('QUAD',  f'QF{cell_idx}',  L=L_QF, K1=K1F),
        make_element('DRIFT', f'D1a{cell_idx}', L=L_D),
        make_element('BEND',  f'B1{cell_idx}',  L=L_B, ANGLE=B_ANGLE, E1=B_E1, E2=B_E2),
        make_element('DRIFT', f'D1b{cell_idx}', L=L_D),
        make_element('QUAD',  f'QD{cell_idx}',  L=L_QD, K1=K1D),
        make_element('DRIFT', f'D2a{cell_idx}', L=L_D),
        make_element('BEND',  f'B2{cell_idx}',  L=L_B, ANGLE=B_ANGLE, E1=B_E1, E2=B_E2),
        make_element('DRIFT', f'D2b{cell_idx}', L=L_D),
    ]

all_pure = []
for i in range(N_CELLS):
    all_pure.extend(make_cell_pure(i))

lat_pure = Lattice(all_pure)

# ---------------------------------------------------------------------------
# Test A: オフモメンタム Twiss
# ---------------------------------------------------------------------------
print("=" * 65)
print("Test A: オフモメンタム Twiss")
print("=" * 65)

def print_twiss_at_entrance(delta):
    tw = lat_pure.twiss(delta=delta)
    twiss0 = tw[0]
    tune_x = twiss0['nx'] / (2 * np.pi)
    tune_y = twiss0['ny'] / (2 * np.pi)
    # SAD qtwiss 方式: COD 依存 Jacobian から計算した周期分散を使用
    etax_0  = twiss0['ex']
    etapx_0 = twiss0['epx']
    # Twiss at s=0 (ring entrance)
    print(f"--- dp = {delta:+.2f} ---")
    print(f"  BetaX  = {twiss0['bx']:.6f} m")
    print(f"  BetaY  = {twiss0['by']:.6f} m")
    print(f"  AlphaX = {twiss0['ax']:.6f}")
    print(f"  AlphaY = {twiss0['ay']:.6f}")
    print(f"  EtaX   = {etax_0:.6f} m")
    print(f"  EtaPX  = {etapx_0:.6f}")
    print(f"  Qx     = {tune_x:.6f}")
    print(f"  Qy     = {tune_y:.6f}")
    print()

for dp in [0.00, -0.02, -0.01, 0.01, 0.02]:
    print_twiss_at_entrance(dp)

# ---------------------------------------------------------------------------
# ラティス構築 (Test B: 複合偏向 BEND with K1)
# ---------------------------------------------------------------------------

def make_cell_k1(cell_idx):
    return [
        make_element('QUAD',  f'QF{cell_idx}',  L=L_QF, K1=K1F),
        make_element('DRIFT', f'D1a{cell_idx}', L=L_D),
        make_element('BEND',  f'BK1_{cell_idx}', L=L_B, ANGLE=B_ANGLE,
                     K1=BK_K1_PYSAD, E1=B_E1, E2=B_E2),
        make_element('DRIFT', f'D1b{cell_idx}', L=L_D),
        make_element('QUAD',  f'QD{cell_idx}',  L=L_QD, K1=K1D),
        make_element('DRIFT', f'D2a{cell_idx}', L=L_D),
        make_element('BEND',  f'BK2_{cell_idx}', L=L_B, ANGLE=B_ANGLE,
                     K1=BK_K1_PYSAD, E1=B_E1, E2=B_E2),
        make_element('DRIFT', f'D2b{cell_idx}', L=L_D),
    ]

all_k1 = []
for i in range(N_CELLS):
    all_k1.extend(make_cell_k1(i))

lat_k1 = Lattice(all_k1)

# ---------------------------------------------------------------------------
# Test B: 複合偏向電磁石, 放射積分
# ---------------------------------------------------------------------------
print("=" * 65)
print(f"Test B: 複合偏向 BEND (K1_SAD={BK_K1_SAD} m^-1 integrated, K1_pysad={BK_K1_PYSAD} m^-2)")
print("=" * 65)

tw_k1  = lat_k1.twiss()
tw0_k1 = tw_k1[0]
tune_x = tw0_k1['nx'] / (2 * np.pi)
tune_y = tw0_k1['ny'] / (2 * np.pi)
etax_0, etapx_0 = periodic_dispersion(lat_k1.elements)

print("--- 線形オプティクス (ring entrance) ---")
print(f"  BetaX  = {tw0_k1['bx']:.6f} m")
print(f"  BetaY  = {tw0_k1['by']:.6f} m")
print(f"  AlphaX = {tw0_k1['ax']:.6f}")
print(f"  AlphaY = {tw0_k1['ay']:.6f}")
print(f"  EtaX   = {etax_0:.6f} m")
print(f"  EtaPX  = {etapx_0:.6f}")
print(f"  Qx     = {tune_x:.6f}")
print(f"  Qy     = {tune_y:.6f}")
print()

print("--- 放射積分 (I1..I5) ---")
emit_result = lat_k1.emittance(energy_eV=ENERGY_EV)
print(f"  I1 = {emit_result['I1']:.6e} m")
print(f"  I2 = {emit_result['I2']:.6e} m^-1")
print(f"  I3 = {emit_result['I3']:.6e} m^-2")
print(f"  I4 = {emit_result['I4']:.6e} m^-1")
print(f"  I5 = {emit_result['I5']:.6e} m^-1")
print()
print("--- エミッタンス ---")
print(f"  Emittance_x        = {emit_result['emittance_x']*1e9:.6f} nm·rad")
print(f"  Energy spread      = {emit_result['energy_spread']*1e3:.6f} × 10^-3")
print(f"  Energy loss U0     = {emit_result['energy_loss_eV']*1e-3:.6f} keV/turn")
print(f"  Mom. compaction α  = {emit_result['momentum_compaction']:.6e}")
print()

print("--- Emittance[] (6D: temit/tradke 直訳, pysad.emit6d) ---")
from pysad.emit6d import emittance6d
e6 = emittance6d(lat_k1, ENERGY_EV)
print(f"  Emittance x        = {e6['emittance_x']:.10e} m")
print(f"  MomentumSpread     = {e6['momentum_spread']:.10e}")
print(f"  MomentumCompaction = {e6['momentum_compaction']:.10e}")
print(f"  EnergyLossU0       = {e6['energy_loss_eV']:.6f} eV")
print()
print("=" * 65)
print("完了")
print("=" * 65)
