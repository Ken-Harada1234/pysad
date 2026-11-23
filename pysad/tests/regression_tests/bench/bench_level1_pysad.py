#!/usr/bin/env python3
"""
bench_level1_pysad.py
=====================
Level 1 ベンチマーク — pysad 側スクリプト

対応する SAD スクリプト: bench_level1.sad

ラティス定義 (双方で共通):
  N_CELLS = 20 の FODO リング
  エネルギー: 2.5 GeV (電子)
  CELL = QF - D1 - B - D1 - QD - D1 - B - D1

テスト内容:
  Test 1: 線形オプティクス (Twiss, チューン, 分散)
  Test 2: COD -- QD0 に DX = 1e-4 m 誤差を印加
  Test 3: 放射積分 I1..I5, 自然エミッタンス, エネルギースプレッド

実行:
  cd "Claude outputs"
  python bench_level1_pysad.py

注記:
  pysad の lat.twiss() には既知の2つの留意点がある:
  (a) チューン: tw[0]['nx'] が (tune × 2π) を返す。tw[-1]['nx'] は2倍になるため、
      チューンは tw[0]['nx']/(2π) で求める。
  (b) 分散: periodic_optics が dispersion=0 を初期値とするため、
      このスクリプトでは一周転送行列から直接周期分散を計算して使う。
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
# ラティスパラメータ (SAD スクリプトと完全に一致させること)
# ---------------------------------------------------------------------------

N_CELLS   = 20           # セル数
ENERGY_EV = 2.5e9        # ビームエネルギー [eV]

L_QF  = 0.3              # QF 長さ [m]
L_QD  = 0.3              # QD 長さ [m]
L_D   = 1.2              # ドリフト長さ [m]
L_B   = 0.8              # 偏向電磁石長さ [m]

K1F   =  1.8             # QF 四極強度 [m^-2]
K1D   = -1.8             # QD 四極強度 [m^-2]

B_ANGLE = np.pi / 20     # 偏向角 [rad] (2π / (N_CELLS×2) = π/20)
B_E1    = B_ANGLE / 2    # 端部角 (対称ハードエッジ) [rad]
B_E2    = B_ANGLE / 2

DX_MISALIGN = 1.0e-4     # QD0 据え付け誤差 DX [m]

L_CELL = L_QF + 2*L_D + L_B + L_QD + 2*L_D + L_B
CIRCUMFERENCE = N_CELLS * L_CELL

print(f"ラティス: {N_CELLS} FODO セル, C = {CIRCUMFERENCE:.3f} m, E = {ENERGY_EV*1e-9:.2f} GeV")
print(f"セル長 = {L_CELL:.3f} m, 偏向角/bend = {np.degrees(B_ANGLE):.4f} deg")
print()

# ---------------------------------------------------------------------------
# セル要素の生成
# ---------------------------------------------------------------------------

def make_cell(cell_idx):
    return [
        # 要素名は SAD と同じ (同名 = 同じファミリー; インスタンスは QD.1, QD.2, ...)
        make_element('QUAD',  'QF', L=L_QF, K1=K1F),
        make_element('DRIFT', 'D1', L=L_D),
        make_element('BEND',  'B',  L=L_B, ANGLE=B_ANGLE, E1=B_E1, E2=B_E2),
        make_element('DRIFT', 'D1', L=L_D),
        make_element('QUAD',  'QD', L=L_QD, K1=K1D),
        make_element('DRIFT', 'D1', L=L_D),
        make_element('BEND',  'B',  L=L_B, ANGLE=B_ANGLE, E1=B_E1, E2=B_E2),
        make_element('DRIFT', 'D1', L=L_D),
    ]

all_elements = []
for i in range(N_CELLS):
    all_elements.extend(make_cell(i))

lat = Lattice(all_elements)

IDX_QF0 = 0   # QF.1
IDX_QD0 = 4   # QD.1

# ---------------------------------------------------------------------------
# 補助: 周期分散を一周転送行列から直接計算
# ---------------------------------------------------------------------------
def periodic_dispersion(elements, delta=0.0):
    """
    一周 6×6 行列から周期閉分散 (etax, etapx) を計算する。
    (I - Mx) @ [eta, etap] = [M[0,5], M[1,5]]
    """
    M = _transport_6x6(elements, delta)
    Mx = M[:2, :2]
    d  = M[:2, 5]
    eta = np.linalg.solve(np.eye(2) - Mx, d)
    return eta[0], eta[1]

def propagate_dispersion(elements, eta0, etap0, delta=0.0):
    """
    分散 (eta0, etap0) を全要素出口で伝搬し、各位置の (s, etax, etapx) を返す。
    """
    result = [(0.0, eta0, etap0)]
    s = 0.0
    eta, etap = eta0, etap0
    for el in elements:
        M = el.transfer_matrix(delta)
        new_eta  = M[0,0]*eta + M[0,1]*etap + M[0,5]
        new_etap = M[1,0]*eta + M[1,1]*etap + M[1,5]
        s += el.L
        result.append((s, new_eta, new_etap))
        eta, etap = new_eta, new_etap
    return result

# ---------------------------------------------------------------------------
# 補助: 位相進行を全要素で計算 (twiss0 の beta/alpha から)
# ---------------------------------------------------------------------------
def compute_phase_table(elements, twiss0_dict, eta0, etap0, delta=0.0):
    """
    全要素出口で s, betx, bety, alfx, alfy, phix, phiy, etax を計算して返す。
    phix/phiy は s=0 からの累積位相 [rad]。
    """
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

        # update for next step
        gx = (1 + ax_new**2) / bx_new
        gy = (1 + ay_new**2) / by_new
        bx = bx_new; ax = ax_new
        by = by_new; ay = ay_new
        eta = new_eta; etap = new_etap
        phix += dpx
        phiy += dpy
        s += el.L

        rows.append(dict(s=s, betx=bx, alfx=ax, bety=by, alfy=ay,
                         phix=phix, phiy=phiy, etax=eta, etapx=etap))
    return rows

# ---------------------------------------------------------------------------
# Test 1: 線形オプティクス
# ---------------------------------------------------------------------------
print("=" * 60)
print("Test 1: 線形オプティクス")
print("=" * 60)

# チューン: tw[0]['nx'] = tune×2π (periodic_optics が返す値)
tw = lat.twiss()
twiss0 = tw[0]  # 初期 Twiss (periodic solution at s=0)

tune_x = twiss0['nx'] / (2 * np.pi)
tune_y = twiss0['ny'] / (2 * np.pi)

# 周期分散
etax_0, etapx_0 = periodic_dispersion(lat.elements)

# 全要素 Twiss テーブル
table = compute_phase_table(lat.elements, twiss0, etax_0, etapx_0)

row0  = table[0]
row_qf = table[IDX_QF0]  # QF.1 入口 (SAD Twiss[...,"QF"] は要素入口)
row_qd = table[IDX_QD0]  # QD.1 入口

print(f"--- ラティス入口 (s=0) ---")
print(f"  betx  = {row0['betx']:.6f} m")
print(f"  bety  = {row0['bety']:.6f} m")
print(f"  alfx  = {row0['alfx']:.6f}")
print(f"  alfy  = {row0['alfy']:.6f}")
print(f"  etax  = {row0['etax']:.6f} m")
print(f"  etapx = {row0['etapx']:.6f}")
print(f"  Qx    = {tune_x:.6f}")
print(f"  Qy    = {tune_y:.6f}")
print()

print(f"--- QF.1 入口 (s = {row_qf['s']:.4f} m) ---")
print(f"  betx  = {row_qf['betx']:.6f} m")
print(f"  bety  = {row_qf['bety']:.6f} m")
print(f"  etax  = {row_qf['etax']:.6f} m")
print()

print(f"--- QD.1 入口 (s = {row_qd['s']:.4f} m) ---")
print(f"  betx  = {row_qd['betx']:.6f} m")
print(f"  bety  = {row_qd['bety']:.6f} m")
print(f"  etax  = {row_qd['etax']:.6f} m")
print()

all_bx  = [r['betx'] for r in table]
all_by  = [r['bety'] for r in table]
all_ex  = [r['etax'] for r in table]
print(f"--- 全体集計 ---")
print(f"  max(betx) = {max(all_bx):.6f} m")
print(f"  min(betx) = {min(all_bx):.6f} m")
print(f"  max(bety) = {max(all_by):.6f} m")
print(f"  min(bety) = {min(all_by):.6f} m")
print(f"  max(etax) = {max(all_ex):.6f} m")
print(f"  min(etax) = {min(all_ex):.6f} m")

# ---------------------------------------------------------------------------
# Test 2: COD (QD0 に DX 誤差)
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("Test 2: COD (QD に DX=1e-4 を印加: SAD の LINE[\"DX\",\"QD.1\"] + CELL CALC と同じ)")
print("=" * 60)

# SAD 側は LINE["DX","QD.1"]=1e-4 → CELL; CALC で設定している。
# QD.1 はファミリー QD の代表なので、CALC 時に要素値として全 QD (20 台) に展開される
# (SAD の既定動作)。1 台だけに誤差を入れる場合は FFS["QD.1 DX 1E-4"] /
# pysad では lat.set_error('QD.1', 'DX', v) を使う。
n_set = lat.set_element('QD', 'DX', DX_MISALIGN)
print(f"  QD.DX = {DX_MISALIGN*1e6:.1f} µm  (要素値として全 {n_set} 台に設定 = SAD LINE[] + CELL CALC)")

try:
    tw_cod = lat.twiss_table()
    cod = np.array([tw_cod[0].x_cod, tw_cod[0].px_cod, tw_cod[0].y_cod, tw_cod[0].py_cod])

    print(f"  閉軌道 (リング入口):")
    print(f"    x    = {cod[0]*1e6:.4f} µm")
    print(f"    px   = {cod[1]*1e6:.4f} µrad")
    print(f"    y    = {cod[2]*1e6:.4f} µm")
    print(f"    py   = {cod[3]*1e6:.4f} µrad")
    print()

    row_qd_cod = tw_cod[IDX_QD0]
    print(f"  COD at QD.1 入口 (s={row_qd_cod.s:.4f} m):")
    print(f"    x    = {row_qd_cod.x_cod*1e6:.4f} µm")
    print(f"    px   = {row_qd_cod.px_cod*1e6:.4f} µrad")
    print()

    print(f"  COD 下の Twiss (リング入口):")
    print(f"    BX   = {tw_cod[0].betx:.6f} m")
    print(f"    Qx   = {tw_cod[-1].mux/(2*np.pi):.6f}")
    print()

    x_all = np.array([r.x_cod for r in tw_cod])
    y_all = np.array([r.y_cod for r in tw_cod])
    print(f"  COD RMS:")
    print(f"    x_rms = {np.sqrt(np.mean(x_all**2))*1e6:.4f} µm")
    print(f"    y_rms = {np.sqrt(np.mean(y_all**2))*1e6:.4f} µm")
except RuntimeError as e:
    print(f"  COD 計算失敗: {e}")

lat.set_element('QD', 'DX', 0.0)  # 誤差リセット (SAD: LINE["DX","QD.1"]=0)
lat.twiss()                       # SAD: CELL; CALC (閉軌道を設計軌道に戻す)

# ---------------------------------------------------------------------------
# Test 3: 放射積分 / エミッタンス
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("Test 3: 放射積分 / エミッタンス")
print("=" * 60)

emit_result = lat.emittance(energy_eV=ENERGY_EV)

print(f"  I1 = {emit_result['I1']:.6e} m")
print(f"  I2 = {emit_result['I2']:.6e} m^-1")
print(f"  I3 = {emit_result['I3']:.6e} m^-2")
print(f"  I4 = {emit_result['I4']:.6e} m^-1")
print(f"  I5 = {emit_result['I5']:.6e} m^-1")
print()
print(f"  Emittance_x   = {emit_result['emittance_x']*1e9:.6f} nm·rad")
print(f"  Energy spread = {emit_result['energy_spread']*1e3:.6f} x10^-3")
print(f"  Energy loss   = {emit_result['energy_loss_eV']*1e-3:.6f} keV/turn")
print(f"  Mom. compact. = {emit_result['momentum_compaction']:.6e}")
print()

print("--- Emittance[] (6D: temit/tradke 直訳, pysad.emit6d) ---")
from pysad.emit6d import emittance6d
e6 = emittance6d(lat, ENERGY_EV)
print(f"  Emittance x        = {e6['emittance_x']:.10e} m")
print(f"  MomentumSpread     = {e6['momentum_spread']:.10e}")
print(f"  MomentumCompaction = {e6['momentum_compaction']:.10e}")
print(f"  EnergyLossU0       = {e6['energy_loss_eV']:.6f} eV")
print()
print("=" * 60)
print("完了")
print("=" * 60)
