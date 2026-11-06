"""
pysad 新機能テスト:
  1. MARK 先頭要素の DX/DY/DPX/DPY → ins_calc の初期軌道オフセット
  2. save() / save(all_params=True) / restore_saved()
  3. fitp(N) / dp 属性 / 多点チューンマッチング
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

import numpy as np
from pysad.elements import Drift, Quadrupole, Bend, Sextupole, Marker
from pysad.ffs import FFSState

# -----------------------------------------------------------------------
# Mock env
# -----------------------------------------------------------------------
class MockEnv:
    def __init__(self, elements):
        self._elements = {e.name: {} for e in elements}

    def build_lattice(self, name):
        class MockLat:
            def __init__(self, elems): self.elements = elems
        return MockLat(list(self._elements.keys()))


def make_ffs(elements):
    env = MockEnv(elements)
    ffs = FFSState(env)
    ffs._beamline_name = "TEST"
    ffs._elements = list(elements)
    ffs._names = [e.name for e in elements]
    s = 0.0
    ffs._s = []
    for e in elements:
        s += e.L
        ffs._s.append(s)
    return ffs


def make_stable_fodo(with_mark=False):
    """
    安定な FODO セルを作る。
    K1=0.5, L_q=0.3, L_drift=0.5 → f=6.67m, half-cell=0.65m, L/(2f)=0.05 → 安定
    """
    elems = []
    if with_mark:
        elems.append(Marker(name="MARK0", DX=0.0, DPX=0.0, DY=0.0, DPY=0.0))
    elems += [
        Quadrupole(name="QF", L=0.3, K1=0.5),
        Drift(name="D1", L=0.5),
        Quadrupole(name="QD", L=0.3, K1=-0.5),
        Drift(name="D2", L=0.5),
        Sextupole(name="SF", L=0.1, K2=2.0),
        Sextupole(name="SD", L=0.1, K2=-2.0),
    ]
    return elems


print("=" * 65)
print("pysad 新機能テスト")
print("=" * 65)

# ======================================================================
# Test 1: MARK 軌道オフセット
# ======================================================================
print("\n[Test 1] MARK DX/DY/DPX/DPY → ins_calc 初期軌道")

elems_mark = make_stable_fodo(with_mark=True)
ffs_m = make_ffs(elems_mark)
ffs_m.cell_calc()

# 1a: DX=DY=0 なら ins_calc の dx=dy=0 (MARK 先頭位置)
mk = ffs_m._elements[0]
mk.DX = 0.0; mk.DPX = 0.0; mk.DY = 0.0; mk.DPY = 0.0
ffs_m.ins_calc()
dx0 = ffs_m._twiss[0]['dx']
dy0 = ffs_m._twiss[0]['dy']
print(f"  (a) MARK DX=DY=0 → _twiss[0] dx={dx0:.4e}, dy={dy0:.4e}")
assert abs(dx0) < 1e-12 and abs(dy0) < 1e-12
print("      PASS")

# 1b: DX=1e-3 でオフセット軌道
mk.DX = 1e-3
mk.DPX = 2e-4
mk.DY = 0.5e-3
mk.DPY = 0.0
ffs_m.ins_calc()
dx1 = ffs_m._twiss[0]['dx']
dy1 = ffs_m._twiss[0]['dy']
print(f"  (b) MARK DX=1e-3,DPX=2e-4,DY=5e-4 → _twiss[0] dx={dx1:.4e}, dy={dy1:.4e}")
assert abs(dx1 - 1e-3) < 1e-12, f"dx1={dx1}"
assert abs(dy1 - 0.5e-3) < 1e-12, f"dy1={dy1}"
print("      PASS")

# 1c: 軌道が伝播していること (要素を進んだ位置で dx ≠ 1e-3)
# _twiss[1] はインデックス1の要素 (MARK なら MARK 自身、インデックス0がMARK)
# MARK 後の QF 出口 = _twiss[2]
dx_after_qf = ffs_m._twiss[2]['dx']
print(f"  (c) QF 出口 _twiss[2].dx={dx_after_qf:.4e}  (軌道伝播)")
assert dx_after_qf != 1e-3, "軌道は QF で変化するはず"
print("      PASS")

# 1d: MARK が先頭でない場合は通常の ins_calc (例外なし)
elems_no_mark = make_stable_fodo(with_mark=False)
ffs_nm = make_ffs(elems_no_mark)
ffs_nm.cell_calc()
ffs_nm.ins_calc()
dx_nm0 = ffs_nm._twiss[0]['dx']
print(f"  (d) MARK なし ins_calc _twiss[0].dx={dx_nm0:.4e}  (COD 由来)")
print("      PASS (例外なし)")

# ======================================================================
# Test 2: save() / save(all_params=True) / restore_saved()
# ======================================================================
print("\n[Test 2] save() / restore_saved()")

elems_s = make_stable_fodo(with_mark=False)
ffs_s = make_ffs(elems_s)
ffs_s.cell_calc()

# 2a: save() — 主変数のみ
qf_el = next(e for e in ffs_s._elements if e.name == "QF")
qd_el = next(e for e in ffs_s._elements if e.name == "QD")
k1_qf = float(qf_el.K1)   # 0.5
k1_qd = float(qd_el.K1)   # -0.5
k2_sf = float(next(e for e in ffs_s._elements if e.name == "SF").K2)  # 2.0
print(f"  (a) Before save: QF K1={k1_qf}, QD K1={k1_qd}, SF K2={k2_sf}")

ffs_s.save()
assert 'QF' in ffs_s._saved_params, "QF が保存されていない"
assert 'K1' in ffs_s._saved_params['QF'], "QF.K1 が保存されていない"
assert abs(ffs_s._saved_params['QF']['K1'] - k1_qf) < 1e-15, "QF K1 値が違う"
assert 'K2' in ffs_s._saved_params.get('SF', {}), "SF.K2 が保存されていない"
print(f"      saved['QF']={ffs_s._saved_params['QF']}")
print("      save() PASS")

# 2b: パラメータを変更してから restore_saved()
qf_el.K1 = 9.99
ffs_s.restore_saved()
k1_restored = float(qf_el.K1)
print(f"  (b) After restore_saved: QF K1={k1_restored} (expect {k1_qf})")
assert abs(k1_restored - k1_qf) < 1e-15
print("      restore_saved() PASS")

# 2c: save(all_params=True) — 全パラメータ
ffs_s.save(all_params=True)
qf_saved_all = ffs_s._saved_params['QF']
print(f"  (c) save(all_params=True) QF saved keys: {sorted(qf_saved_all.keys())}")
assert 'L'  in qf_saved_all, "L が保存されていない"
assert 'K1' in qf_saved_all, "K1 が保存されていない"
assert 'DX' in qf_saved_all, "DX が保存されていない"
print("      save(all_params=True) PASS")

# 2d: Drift の主変数は L
drift_el = next(e for e in ffs_s._elements if e.name == "D1")
ffs_s.save()
d1_saved = ffs_s._saved_params.get('D1', {})
print(f"  (d) Drift D1 saved (primary): {d1_saved}")
assert 'L' in d1_saved, f"Drift D1 の主変数 L が保存されていない: {d1_saved}"
assert list(d1_saved.keys()) == ['L'], f"Drift は L だけのはず: {d1_saved}"
print("      Drift primary param PASS")

# ======================================================================
# Test 3: fitp() / dp / 多点マッチング
# ======================================================================
print("\n[Test 3] fitp / dp / _cell_calc_at_delta")

ffs_f = make_ffs(make_stable_fodo(with_mark=False))
ffs_f.cell_calc()

# 3a: 基本属性確認
assert ffs_f._fit_npoints == 1
assert ffs_f._dp == 0.0
ffs_f.fitp(2)
assert ffs_f._fit_npoints == 2
ffs_f.dp = 5e-3
assert ffs_f._dp == 5e-3
print(f"  (a) fitp(2), dp=5e-3 → _fit_npoints={ffs_f._fit_npoints}, _dp={ffs_f._dp}  PASS")

ffs_f.fitp(1)
assert ffs_f._fit_npoints == 1
print("  (b) fitp(1) → _fit_npoints=1  PASS")

# 3b: _cell_calc_at_delta
ffs_f.dp = 5e-3
tw0  = ffs_f._cell_calc_at_delta(0.0)
tw_p = ffs_f._cell_calc_at_delta(+5e-3)
tw_m = ffs_f._cell_calc_at_delta(-5e-3)
assert all(tw is not None for tw in [tw0, tw_p, tw_m])
nx0  = tw0['nx'] / (2 * np.pi)
nx_p = tw_p['nx'] / (2 * np.pi)
nx_m = tw_m['nx'] / (2 * np.pi)
print(f"  (c) Tune NX: delta=0 → {nx0:.6f}, +dp → {nx_p:.6f}, -dp → {nx_m:.6f}")
print(f"      Δtune(+dp)={nx_p - nx0:+.2e}, Δtune(-dp)={nx_m - nx0:+.2e}")
print("      _cell_calc_at_delta PASS")

# ======================================================================
# Test 4: 多点マッチング (npoints=2)
# ======================================================================
print("\n[Test 4] 多点チューンマッチング (fitp=2, dp=5e-3)")

ffs4 = make_ffs(make_stable_fodo(with_mark=False))
ffs4.cell_calc()

nx_init = ffs4.Twiss("NX", "$$$") / (2 * np.pi)
ny_init = ffs4.Twiss("NY", "$$$") / (2 * np.pi)
print(f"  Initial tune: NX={nx_init:.6f}, NY={ny_init:.6f}")

# 目標: 現在のチューンに近い値 (収束の確認)
target_nx = nx_init
target_ny = ny_init

ffs4.fitp(2)
ffs4.dp = 5e-3
ffs4._free_params = [
    (i, e.name, "K1")
    for i, e in enumerate(ffs4._elements)
    if e.name in ("QF", "QD")
]
ffs4._constraints = []
ffs4.fit_add("$$$", "NX", target_nx, npoints=2)
ffs4.fit_add("$$$", "NY", target_ny, npoints=2)

success = ffs4.go(niter=10)
ffs4.cell_calc()
nx_fin = ffs4.Twiss("NX", "$$$") / (2 * np.pi)
ny_fin = ffs4.Twiss("NY", "$$$") / (2 * np.pi)
print(f"  After go (npoints=2): NX={nx_fin:.6f} (target={target_nx:.6f})")
print(f"                        NY={ny_fin:.6f} (target={target_ny:.6f})")
assert abs(nx_fin - target_nx) < 2e-3, f"NX 収束失敗: diff={abs(nx_fin-target_nx)}"
assert abs(ny_fin - target_ny) < 2e-3, f"NY 収束失敗: diff={abs(ny_fin-target_ny)}"
print(f"  収束: {success}  PASS")

# ======================================================================
# Test 5: fit_add npoints=None → fitp() 設定を継承
# ======================================================================
print("\n[Test 5] fit_add npoints=None → fitp() 設定を継承")

ffs5 = make_ffs(make_stable_fodo(with_mark=False))
ffs5.fitp(3)
ffs5.dp = 1e-3
ffs5.fit_add("$$$", "NX", 0.1)
assert ffs5._constraints[0]['npoints'] is None, "npoints=None のはず"
print(f"  constraints[0]['npoints'] = {ffs5._constraints[0]['npoints']}  (None = fitp 継承)")
print("  PASS")

print("\n" + "=" * 65)
print("ALL TESTS PASSED")
print("=" * 65)
