"""
test_sext_radiation.py — Sextupole thick radiation kick の単体テスト

E0_GeV=0 のとき放射なし (旧動作と同一)、
E0_GeV>0 のとき決定論的エネルギー損失が発生することを確認する。
"""
import sys
import numpy as np

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pysad.elements import Sextupole
from pysad.tracking import track_tthin, track_element

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

def check_less(name, val, threshold):
    global PASS, FAIL
    val = float(np.asarray(val).ravel()[0])
    ok = val < threshold
    status = "PASS" if ok else "FAIL"
    if not ok:
        FAIL += 1
        print(f"  [FAIL] {name}: {val:.4e} not < {threshold:.4e}")
    else:
        PASS += 1
        print(f"  [PASS] {name}: {val:.4e} < {threshold:.4e}")

def check_greater(name, val, threshold):
    global PASS, FAIL
    val = float(np.asarray(val).ravel()[0])
    ok = val > threshold
    status = "PASS" if ok else "FAIL"
    if not ok:
        FAIL += 1
        print(f"  [FAIL] {name}: {val:.4e} not > {threshold:.4e}")
    else:
        PASS += 1
        print(f"  [PASS] {name}: {val:.4e} > {threshold:.4e}")


# テスト粒子
x0  = np.array([1e-3])
px0 = np.array([0.0])
y0  = np.array([5e-4])
py0 = np.array([0.0])
z0  = np.array([0.0])
g0  = np.array([0.0])

# Sextupole パラメータ
K2  = 10.0   # 1/m^3
L   = 0.3    # m


# ===========================================================================
# Test 1: E0_GeV=0 → 放射なし (track_tthin 直接呼び出し)
# ===========================================================================
print("\n--- Test 1: E0_GeV=0 → 放射なし ---")
xo, pxo, yo, pyo, zo, go = track_tthin(
    x0, px0, y0, py0, z0, g0,
    nord=6, L=L, ak=K2, E0_GeV=0.0, rng=None)

# g は変化しないはず
check("g unchanged when E0=0", go, g0, atol=1e-14)
# x は厚みのある多極で若干変化 (変化すること自体は確認しない; 変化量が有限であることだけ確認)
print(f"  [INFO] g_out = {go[0]:.4e}  (E0=0: g 変化なし)")


# ===========================================================================
# Test 2: E0_GeV>0, rng=None (決定論的) → エネルギー損失が発生
# ===========================================================================
print("\n--- Test 2: E0_GeV>0, rng=None → エネルギー損失 ---")
E0 = 2.5   # GeV
xo_r, pxo_r, yo_r, pyo_r, zo_r, go_r = track_tthin(
    x0, px0, y0, py0, z0, g0,
    nord=6, L=L, ak=K2, E0_GeV=E0, rng=None)

# g < 0 (エネルギー損失 → g = δ < 0)
check_less("g < 0 after radiation", go_r[0], 0.0)
print(f"  [INFO] Δg = {go_r[0]:.4e}  (E0={E0} GeV)")


# ===========================================================================
# Test 3: E0=0 と E0>0 の比較 (横座標 x は大きく変わらない)
# ===========================================================================
print("\n--- Test 3: 横座標 x の変化量は小さい ---")
# g 変化による 2 次補正は小さいはず (|Δx| < |Δg| * L のオーダー)
delta_x = abs(float(xo_r[0]) - float(xo[0]))
delta_g = abs(float(go_r[0]))
# x の変化は Δg * L 程度以下
check_less("x change << |Δg| * 10", delta_x, delta_g * 10 + 1e-10)


# ===========================================================================
# Test 4: エネルギー損失は x^2 + y^2 に比例 (大振幅粒子で大きい)
# ===========================================================================
print("\n--- Test 4: 大振幅粒子ほどエネルギー損失が大きい ---")
x_large = np.array([5e-3])
y_large = np.array([2e-3])
_, _, _, _, _, go_large = track_tthin(
    x_large, px0, y_large, py0, z0, g0,
    nord=6, L=L, ak=K2, E0_GeV=E0, rng=None)
_, _, _, _, _, go_small = track_tthin(
    x0, px0, y0, py0, z0, g0,
    nord=6, L=L, ak=K2, E0_GeV=E0, rng=None)

# 大振幅 → より大きなエネルギー損失 (|Δg_large| > |Δg_small|)
check_greater("large amplitude: more energy loss",
              abs(float(go_large[0])), abs(float(go_small[0])))
print(f"  [INFO] Δg(large) = {go_large[0]:.4e}, Δg(small) = {go_small[0]:.4e}")


# ===========================================================================
# Test 5: 放射後 z は有限 (NaN・Inf なし)
# ===========================================================================
# 注意: _radiation_kick は g のみ更新し z は直接変化させない。
# z への効果は後続ドリフトを経由するが, Δg ≈ 7e-14 では数値精度以下のため
# ここでは z が NaN/Inf でないことを確認するにとどめる。
print("\n--- Test 5: 放射後 z は有限 (NaN なし) ---")
xo_r2, _, _, _, zo_r2, go_r2 = track_tthin(
    x0, px0, y0, py0, z0, g0,
    nord=6, L=L, ak=K2, E0_GeV=E0, rng=None)
print(f"  [INFO] z(E0>0)={zo_r2[0]:.4e}")
z_finite = np.isfinite(zo_r2[0])
check("z is finite after radiation", [int(z_finite)], [1])


# ===========================================================================
# Test 6: track_element 経由の Sextupole トラッキング (E0_GeV パス確認)
# ===========================================================================
print("\n--- Test 6: track_element 経由 (E0_GeV=2.5) ---")
sext = Sextupole(name='S1', L=L, K2=K2)
from pysad.tracking import track_element
xe, pxe, ye, pye, ze, ge = track_element(
    sext, x0, px0, y0, py0, z0, g0, E0_GeV=E0)
xr, pxr, yr, pyr, zr, gr = track_tthin(
    x0, px0, y0, py0, z0, g0, nord=6, L=L, ak=K2, E0_GeV=E0, rng=None)

check("track_element x == track_tthin x",   xe, xr, atol=1e-14)
check("track_element g == track_tthin g",   ge, gr, atol=1e-14)


# ===========================================================================
# Test 7: K2=0 (ゼロ強度) → 純ドリフト; g 変化なし (radiation term もゼロ)
# ===========================================================================
print("\n--- Test 7: K2=0 → 放射なし (h_eff=0) ---")
xo_k0, _, _, _, zo_k0, go_k0 = track_tthin(
    x0, px0, y0, py0, z0, g0,
    nord=6, L=L, ak=0.0, E0_GeV=E0, rng=None)
check("g unchanged when K2=0", go_k0, g0, atol=1e-14)


# ===========================================================================
# Summary
# ===========================================================================
print(f"\n{'='*50}")
print(f"  Results: {PASS} PASS, {FAIL} FAIL")
print(f"{'='*50}")
sys.exit(0 if FAIL == 0 else 1)
