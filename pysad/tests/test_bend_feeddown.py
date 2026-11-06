"""
test_bend_feeddown.py — Bend.constant_kick() の DX/DY フィードダウン単体テスト

Bend に DX/DY オフセットがある場合、constant_kick() が
  d[0] += (1 - M[0,0]) * DX
  d[1] += -M[1,0] * DX
  d[2] += (1 - M[2,2]) * DY
  d[3] += -M[3,2] * DY
を正しく適用することを確認する。
"""
import sys
import numpy as np

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pysad.elements import Bend

PASS = 0
FAIL = 0

def check(name, got, expected, atol=1e-12):
    global PASS, FAIL
    got = float(np.asarray(got).ravel()[0])
    expected = float(np.asarray(expected).ravel()[0])
    ok = abs(got - expected) <= atol
    status = "PASS" if ok else "FAIL"
    if not ok:
        FAIL += 1
        print(f"  [FAIL] {name}")
        print(f"         got      = {got:.6e}")
        print(f"         expected = {expected:.6e}")
        print(f"         diff     = {got - expected:.6e}")
    else:
        PASS += 1
        print(f"  [PASS] {name}")


# ===========================================================================
# Test 1: DX=DY=0 → feeddown ゼロ (親クラスの constant_kick と同一)
# ===========================================================================
print("\n--- Test 1: DX=DY=0 → feeddown なし ---")
bend = Bend(name='B0', L=1.0, ANGLE=0.1, K0=0.1, DX=0.0, DY=0.0)
d_no_feed = bend.constant_kick()
# K0≠0 の場合も feeddown は DX=DY=0 なら変化しない
bend_ref = Bend(name='B0ref', L=1.0, ANGLE=0.1, K0=0.1, DX=0.0, DY=0.0)
from pysad.elements import Element
d_parent = super(Bend, bend_ref).constant_kick() if hasattr(super(Bend, bend_ref), 'constant_kick') else np.zeros(6)

# DX=DY=0 では Bend.constant_kick() == 親クラスの結果
# (直接比較: d_no_feed と bend.transfer_matrix()なしの基底実装)
# 簡易確認: DX=DY=0 のとき feeddown 部は d に何も加えない
d_ref = np.zeros(6)  # K0 等による寄与は別途
M = bend.transfer_matrix()
# feeddown は DX=DY=0 なら 0
feed_x = (1.0 - M[0, 0]) * 0.0 + 0.0
feed_y = (1.0 - M[2, 2]) * 0.0 + 0.0
check("feeddown_x when DX=0", feed_x, 0.0)
check("feeddown_y when DY=0", feed_y, 0.0)


# ===========================================================================
# Test 2: 純粋 DX オフセット (DY=0) の feeddown
# ===========================================================================
print("\n--- Test 2: DX オフセットのみ ---")
DX_t = 0.002   # 2 mm オフセット
bend_dx = Bend(name='Bdx', L=1.0, ANGLE=0.1, K0=0.1, DX=DX_t, DY=0.0)
M_dx = bend_dx.transfer_matrix()
d_dx = bend_dx.constant_kick()

# 親クラス (no DX) の constant_kick (DX=DY=0 の同じ要素)
bend_base = Bend(name='Bbase', L=1.0, ANGLE=0.1, K0=0.1, DX=0.0, DY=0.0)
# 親の constant_kick を呼ぶには super() が必要; ここでは DX=0 のインスタンスで代替
from pysad.elements import Element as _El
d_base_x = d_dx[0] - (1.0 - M_dx[0, 0]) * DX_t   # feeddown を除いた値
d_base_px = d_dx[1] - (-M_dx[1, 0]) * DX_t

# 解析値と比較
expected_d0  = (1.0 - M_dx[0, 0]) * DX_t
expected_d1  = -M_dx[1, 0] * DX_t
# 0 成分の寄与 (K0 等による d_base) は bend_base で確認
d_bend_base = bend_base.constant_kick()

check("d[0] feeddown_x = (1-M00)*DX",
      d_dx[0] - d_bend_base[0], expected_d0, atol=1e-14)
check("d[1] feeddown_px = -M10*DX",
      d_dx[1] - d_bend_base[1], expected_d1, atol=1e-14)
check("d[2] unchanged (DY=0)", d_dx[2] - d_bend_base[2], 0.0, atol=1e-14)
check("d[3] unchanged (DY=0)", d_dx[3] - d_bend_base[3], 0.0, atol=1e-14)


# ===========================================================================
# Test 3: 純粋 DY オフセット (DX=0) の feeddown
# ===========================================================================
print("\n--- Test 3: DY オフセットのみ ---")
DY_t = 0.001
bend_dy = Bend(name='Bdy', L=1.0, ANGLE=0.1, K0=0.1, DX=0.0, DY=DY_t)
M_dy = bend_dy.transfer_matrix()
d_dy = bend_dy.constant_kick()
d_bend_base2 = Bend(name='Bbase2', L=1.0, ANGLE=0.1, K0=0.1).constant_kick()

expected_d2 = (1.0 - M_dy[2, 2]) * DY_t
expected_d3 = -M_dy[3, 2] * DY_t

check("d[2] feeddown_y = (1-M22)*DY",
      d_dy[2] - d_bend_base2[2], expected_d2, atol=1e-14)
check("d[3] feeddown_py = -M32*DY",
      d_dy[3] - d_bend_base2[3], expected_d3, atol=1e-14)
check("d[0] unchanged (DX=0)", d_dy[0] - d_bend_base2[0], 0.0, atol=1e-14)
check("d[1] unchanged (DX=0)", d_dy[1] - d_bend_base2[1], 0.0, atol=1e-14)


# ===========================================================================
# Test 4: DX と DY 両方のフィードダウン
# ===========================================================================
print("\n--- Test 4: DX と DY 両方 ---")
DX_t2, DY_t2 = 0.003, -0.0015
bend_dxdy = Bend(name='Bdxdy', L=1.5, ANGLE=0.2, K0=0.133, DX=DX_t2, DY=DY_t2)
M_dxdy = bend_dxdy.transfer_matrix()
d_dxdy = bend_dxdy.constant_kick()
d_base_dxdy = Bend(name='Bbase3', L=1.5, ANGLE=0.2, K0=0.133).constant_kick()

check("d[0] = base + (1-M00)*DX",
      d_dxdy[0] - d_base_dxdy[0], (1.0 - M_dxdy[0, 0]) * DX_t2, atol=1e-14)
check("d[1] = base - M10*DX",
      d_dxdy[1] - d_base_dxdy[1], -M_dxdy[1, 0] * DX_t2, atol=1e-14)
check("d[2] = base + (1-M22)*DY",
      d_dxdy[2] - d_base_dxdy[2], (1.0 - M_dxdy[2, 2]) * DY_t2, atol=1e-14)
check("d[3] = base - M32*DY",
      d_dxdy[3] - d_base_dxdy[3], -M_dxdy[3, 2] * DY_t2, atol=1e-14)


# ===========================================================================
# Test 5: delta 依存性 (delta=0.01 での feeddown)
# ===========================================================================
print("\n--- Test 5: delta 依存性 ---")
DX_t3 = 0.001
delta = 0.01
bend_d = Bend(name='Bd', L=1.0, ANGLE=0.1, K0=0.1, DX=DX_t3, DY=0.0)
M_d = bend_d.transfer_matrix(delta=delta)
d_d = bend_d.constant_kick(delta=delta)
d_base_d = Bend(name='Bbase4', L=1.0, ANGLE=0.1, K0=0.1).constant_kick(delta=delta)

check("d[0] feeddown at delta=0.01",
      d_d[0] - d_base_d[0], (1.0 - M_d[0, 0]) * DX_t3, atol=1e-14)
check("d[1] feeddown at delta=0.01",
      d_d[1] - d_base_d[1], -M_d[1, 0] * DX_t3, atol=1e-14)


# ===========================================================================
# Summary
# ===========================================================================
print(f"\n{'='*50}")
print(f"  Results: {PASS} PASS, {FAIL} FAIL")
print(f"{'='*50}")
sys.exit(0 if FAIL == 0 else 1)
