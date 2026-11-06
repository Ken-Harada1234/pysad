"""
FFSState error spreading functions test.
FFSState を直接インスタンス化してエラーばらまき関数をテストする。
env を mock し、_elements/_names/_s を直接セット。
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

import numpy as np
from pysad.elements import Drift, Quadrupole, Bend, Sextupole
from pysad.ffs import FFSState

# -----------------------------------------------------------------------
# Mock env (build_latticのみ必要)
# -----------------------------------------------------------------------
class MockEnv:
    def __init__(self, elements):
        self._elements = elements

    def build_lattice(self, name):
        class MockLat:
            def __init__(self, elems):
                self.elements = elems
        return MockLat(self._elements)

# -----------------------------------------------------------------------
# シンプルな FODO セルを作る
# -----------------------------------------------------------------------
def make_fodo_elements():
    """QF - D - B - D - QD - D - B - D の FODO セル"""
    qf = Quadrupole(name="QF", L=0.5, K1=1.2)
    d1 = Drift(name="D1", L=1.0)
    b1 = Bend(name="B1", L=2.0, K0=0.05, ROTATE=0.0)
    d2 = Drift(name="D2", L=1.0)
    qd = Quadrupole(name="QD", L=0.5, K1=-1.2)
    d3 = Drift(name="D3", L=1.0)
    b2 = Bend(name="B2", L=2.0, K0=0.05, ROTATE=0.0)
    d4 = Drift(name="D4", L=1.0)
    sf = Sextupole(name="SF", L=0.3, K2=5.0)
    sd = Sextupole(name="SD", L=0.3, K2=-5.0)
    return [qf, d1, b1, d2, qd, d3, b2, d4, sf, sd]

def make_ffs_from_elements(elements):
    """FFSState を直接構築 (mock env 経由)"""
    env = MockEnv(elements)
    ffs = FFSState(env)
    # use() の代わりに直接セット
    ffs._beamline_name = "TEST"
    ffs._elements = elements
    ffs._names = [e.name for e in elements]
    s = 0.0
    ffs._s = []
    for e in elements:
        s += e.L
        ffs._s.append(s)
    return ffs

# -----------------------------------------------------------------------
# テスト
# -----------------------------------------------------------------------

print("=" * 60)
print("FFSState error spreading test")
print("=" * 60)

elements = make_fodo_elements()
ffs = make_ffs_from_elements(elements)

# ------ Test 1: apply_misalignment ------
print("\n[Test 1] apply_misalignment Q* DX=1e-4 (sigma), DY=2e-4")
# 初期値確認
qf0 = next(e for e in ffs._elements if e.name == "QF")
print(f"  Before: QF DX={qf0.DX}, DY={qf0.DY}")

ffs.apply_misalignment("Q*", dx_sigma=1e-4, dy_sigma=2e-4, gcut=3.0, seed=42)

qf1 = next(e for e in ffs._elements if e.name == "QF")
qd1 = next(e for e in ffs._elements if e.name == "QD")
print(f"  After:  QF DX={qf1.DX:.6e}, DY={qf1.DY:.6e}")
print(f"  After:  QD DX={qd1.DX:.6e}, DY={qd1.DY:.6e}")
assert qf1.DX != 0.0, "QF DX should be nonzero after misalignment"
assert qf1.DY != 0.0, "QF DY should be nonzero after misalignment"
assert qd1.DX != 0.0, "QD DX should be nonzero after misalignment"
print("  PASS")

# constant_kick を確認
# 統合 (10/13) で厚肉の式に変更: d[1] = K0 - M[1,0]*DX。
# 10/11 系の期待値 K0 - K1*L*DX は符号が逆 (磁石が +DX にずれると x=0 の粒子は
# 磁石から見て -DX にいるので集束四極では +x 向きに蹴られる)。厚肉の式は
# Level1 Test2 で SAD の COD と一致することを確認済み。
M_q = qf1.transfer_matrix()
d = qf1.constant_kick()
expected_px = float(qf1.K0) - M_q[1, 0] * float(qf1.DX)
print(f"  QF constant_kick[1]={d[1]:.6e}, expected={expected_px:.6e}")
assert abs(d[1] - expected_px) < 1e-15, "constant_kick mismatch"
print("  constant_kick PASS")

# ------ Test 2: apply_field_error ------
print("\n[Test 2] apply_field_error B* K0 sigma=1e-3 (relative)")
b1_0 = next(e for e in ffs._elements if e.name == "B1")
k0_nom = float(b1_0.K0)
print(f"  Before: B1 K0={k0_nom:.6e}")

ffs.apply_field_error("B*", param="K0", sigma=1e-3, relative=True, gcut=3.0, seed=7)

b1_1 = next(e for e in ffs._elements if e.name == "B1")
k0_new = float(b1_1.K0)
print(f"  After:  B1 K0={k0_new:.6e}")
assert k0_new != k0_nom, "B1 K0 should be changed after field error"
rel_err = abs(k0_new - k0_nom) / abs(k0_nom)
print(f"  Relative error = {rel_err:.4e} (expect ~1e-3)")
assert rel_err < 5e-3, f"Relative error too large: {rel_err}"
print("  PASS")

# 4極の K1 誤差テスト
print("\n[Test 2b] apply_field_error Q* K1 sigma=5e-4 (relative)")
qf2 = next(e for e in ffs._elements if e.name == "QF")
k1_nom = float(qf2.K1)

ffs.apply_field_error("Q*", param="K1", sigma=5e-4, relative=True, gcut=3.0, seed=13)

qf3 = next(e for e in ffs._elements if e.name == "QF")
k1_new = float(qf3.K1)
rel_k1 = abs(k1_new - k1_nom) / abs(k1_nom)
print(f"  QF K1: {k1_nom:.6f} -> {k1_new:.6f}, rel={rel_k1:.4e}")
assert k1_new != k1_nom, "QF K1 should change"
print("  PASS")

# ------ Test 3: apply_rotation ------
print("\n[Test 3] apply_rotation Q* sigma=1e-4 rad")
qf4 = next(e for e in ffs._elements if e.name == "QF")
print(f"  Before: QF ROTATION={qf4.ROTATION}")

ffs.apply_rotation("Q*", sigma=1e-4, gcut=3.0, seed=99)

qf5 = next(e for e in ffs._elements if e.name == "QF")
qd5 = next(e for e in ffs._elements if e.name == "QD")
print(f"  After:  QF ROTATION={qf5.ROTATION:.6e}")
print(f"  After:  QD ROTATION={qd5.ROTATION:.6e}")
assert qf5.ROTATION != 0.0, "QF ROTATION should be nonzero"
print("  PASS")

# ROTATION が transfer_matrix に効いているか確認
M_no_rot = Quadrupole(name="QF_test", L=0.5, K1=1.2, ROTATION=0.0).transfer_matrix()
M_with_rot = Quadrupole(name="QF_test2", L=0.5, K1=1.2, ROTATION=qf5.ROTATION).transfer_matrix()
coupling = abs(M_with_rot[0, 2]) + abs(M_with_rot[2, 0])
print(f"  Coupling terms |M[0,2]|+|M[2,0]| = {coupling:.6e} (expect nonzero with rotation)")
assert coupling > 1e-10, "ROTATION should introduce coupling"
print("  transfer_matrix coupling PASS")

# ------ Test 4: reset_errors ------
print("\n[Test 4] reset_errors Q* (DX, DY, ROTATION)")
ffs.reset_errors("Q*", dx=True, dy=True, rotation=True, field=False)

qf6 = next(e for e in ffs._elements if e.name == "QF")
qd6 = next(e for e in ffs._elements if e.name == "QD")
print(f"  After reset: QF DX={qf6.DX}, DY={qf6.DY}, ROTATION={qf6.ROTATION}")
assert float(qf6.DX) == 0.0, "QF DX should be 0 after reset"
assert float(qf6.DY) == 0.0, "QF DY should be 0 after reset"
assert float(qf6.ROTATION) == 0.0, "QF ROTATION should be 0 after reset"
# K1 should NOT be reset
assert float(qf6.K1) != 1.2, "QF K1 should still have field error (not reset)"
print("  PASS")

print("\n[Test 4b] reset_errors Q* field=['K1'] (K1 → 0.0)")
k1_before = float(qf6.K1)
ffs.reset_errors("Q*", dx=False, dy=False, rotation=False, field=['K1'])
qf7 = next(e for e in ffs._elements if e.name == "QF")
print(f"  QF K1 before field reset: {k1_before:.6f}")
print(f"  QF K1 after field reset:  {float(qf7.K1):.6f}")
assert float(qf7.K1) == 0.0, "QF K1 should be 0.0 after field reset"
print("  PASS (field reset sets K to 0 as documented)")

print("\n" + "=" * 60)
print("ALL TESTS PASSED")
print("=" * 60)
