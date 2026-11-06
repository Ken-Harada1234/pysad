"""
test_coord.py — COORD 要素の単体テスト

tcoord.f の各分岐とエッジケースを検証する。
"""
import sys
import os
import numpy as np

# pysad パッケージへのパス
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pysad.elements import Coord
from pysad.tracking import track_coord

PASS = 0
FAIL = 0

def check(name, got, expected, atol=1e-12):
    global PASS, FAIL
    got = np.atleast_1d(np.asarray(got, dtype=float))
    expected = np.atleast_1d(np.asarray(expected, dtype=float))
    ok = np.allclose(got, expected, atol=atol, rtol=0)
    status = "PASS" if ok else "FAIL"
    if not ok:
        FAIL += 1
        print(f"  [{status}] {name}")
        print(f"         got      = {got}")
        print(f"         expected = {expected}")
        print(f"         diff     = {got - expected}")
    else:
        PASS += 1
        print(f"  [{status}] {name}")


# ===========================================================================
# Test 1: 恒等変換 (すべての角度=0, オフセット=0)
# DIR=True で何も変化しないことを確認
# ===========================================================================
print("\n--- Test 1: 恒等変換 (角度・オフセットすべて 0, DIR=True) ---")
el = Coord(name='C0', CHI1=0, CHI2=0, CHI3=0, DX=0, DY=0, DZ=0, DIR=True)
x0  = np.array([0.001])
px0 = np.array([0.002])
y0  = np.array([-0.001])
py0 = np.array([0.003])
z0  = np.array([0.0005])
g0  = np.array([0.01])
xo, pxo, yo, pyo, zo, go = track_coord(el, x0, px0, y0, py0, z0, g0)
check("x unchanged",  xo,  x0)
check("px unchanged", pxo, px0)
check("y unchanged",  yo,  y0)
check("py unchanged", pyo, py0)
check("z unchanged",  zo,  z0)
check("g unchanged",  go,  g0)


# ===========================================================================
# Test 2: 恒等変換 (DIR=False も同様)
# ===========================================================================
print("\n--- Test 2: 恒等変換 (角度・オフセットすべて 0, DIR=False) ---")
el = Coord(name='C0', CHI1=0, CHI2=0, CHI3=0, DX=0, DY=0, DZ=0, DIR=False)
xo, pxo, yo, pyo, zo, go = track_coord(el, x0, px0, y0, py0, z0, g0)
check("x unchanged",  xo,  x0)
check("px unchanged", pxo, px0)
check("y unchanged",  yo,  y0)
check("py unchanged", pyo, py0)
check("z unchanged",  zo,  z0)
check("g unchanged",  go,  g0)


# ===========================================================================
# Test 3: 転送行列との一致 (小角度; 線形近似との比較)
# ===========================================================================
print("\n--- Test 3: 転送行列との一致 (小角度, CHI3=0.01 rad) ---")
theta = 0.01   # rad (小角度)
el = Coord(name='C1', CHI3=theta, DIR=True)
M = el.transfer_matrix()

# 小振幅粒子 (線形近似が成立する範囲)
x0s  = np.array([1e-4])
px0s = np.array([5e-5])
y0s  = np.array([2e-4])
py0s = np.array([3e-5])
z0s  = np.array([0.0])
g0s  = np.array([0.0])

xo, pxo, yo, pyo, zo, go = track_coord(el, x0s, px0s, y0s, py0s, z0s, g0s)
coords_in  = np.array([x0s[0], px0s[0], y0s[0], py0s[0], z0s[0], g0s[0]])
coords_mat = M @ coords_in

check("x  (matrix vs track)",  xo,  [coords_mat[0]], atol=1e-9)
check("px (matrix vs track)",  pxo, [coords_mat[1]], atol=1e-9)
check("y  (matrix vs track)",  yo,  [coords_mat[2]], atol=1e-9)
check("py (matrix vs track)",  pyo, [coords_mat[3]], atol=1e-9)
check("z  (matrix vs track)",  zo,  [coords_mat[4]], atol=1e-9)
check("g  unchanged",          go,  g0s)


# ===========================================================================
# Test 4: DIR=False のオフセット適用確認
# SAD では DIR=False は DIR=True の単純な逆変換ではない。
# 出口側の変位は dx1=DX, dy1=-DY, dz1=-DZ として加算される (tcoord.f)。
# CHI1=CHI2=CHI3=0 でオン軸粒子 (px=py=z=0) の場合:
#   x_out = 0 + DX, y_out = 0 - DY, z_out = -DZ (r33=1 のとき)
# ===========================================================================
print("\n--- Test 4: DIR=False のオフセット適用 (角度=0, オン軸粒子) ---")
DX_t, DY_t, DZ_t = 0.002, 0.001, 0.0005
el_out = Coord(name='Cout', DX=DX_t, DY=DY_t, DZ=DZ_t,
               CHI1=0, CHI2=0, CHI3=0, DIR=False)
# オン軸粒子 (x=y=z=px=py=g=0)
x0_  = np.array([0.0])
px0_ = np.array([0.0])
y0_  = np.array([0.0])
py0_ = np.array([0.0])
z0_  = np.array([0.0])
g0_  = np.array([0.0])
xo, pxo, yo, pyo, zo, go = track_coord(el_out, x0_, px0_, y0_, py0_, z0_, g0_)
# chi1=chi2=chi3=0 → R=I, r33=1
# dx0=dy0=dz0=0 (exit)
# xi=0, yi=0, zi=0 → xf=yf=0, zf=dz1=-DZ
# pzf=1 → z_out=zf/1=-DZ
# x_out=0-0*z_out+dx1=DX,  y_out=0-0*z_out+dy1=-DY
check("x_out = +DX (DIR=False)",  xo,   [DX_t],  atol=1e-14)
check("y_out = -DY (DIR=False)",  yo,   [-DY_t], atol=1e-14)
check("z_out = -DZ (DIR=False)",  zo,   [-DZ_t], atol=1e-14)
check("px unchanged",             pxo,  [0.0],   atol=1e-14)
check("py unchanged",             pyo,  [0.0],   atol=1e-14)
check("g unchanged",              go,   [0.0],   atol=1e-14)


# ===========================================================================
# Test 5: 純粋 Z 軸回転 (CHI3 のみ; x-y 平面の回転)
# オン軸粒子 (z=0, px=py=0) でのオフセット検証
# ===========================================================================
print("\n--- Test 5: 純粋 Z 軸回転 (CHI3=pi/6, オン軸, DIR=True) ---")
from math import pi, cos, sin
angle = pi / 6   # 30度
el = Coord(name='Cz', CHI3=angle, DIR=True)
# z=0, px=py=0 のオン軸粒子: pzi=1, xi=x, yi=y, zi=0
# xf = r11*x + r12*y = (c1*c3)*x + (-c2*s3)*y
# c1=1, s1=0, c2=1, s2=0, c3=cos(angle), s3=sin(angle)
# r11=c3, r12=-s3, r21=s3, r22=c3, r33=1
# x_out = xf - px_out*z_out = xf - r11*px*0 = xf = c3*x - s3*y
# y_out = yf = s3*x + c3*y
x_in = np.array([1.0])
y_in = np.array([0.0])
px_in = np.array([0.0])
py_in = np.array([0.0])
z_in  = np.array([0.0])
g_in  = np.array([0.0])
xo, pxo, yo, pyo, zo, go = track_coord(el, x_in, px_in, y_in, py_in, z_in, g_in)
c3, s3 = cos(angle), sin(angle)
check("x after z-rotation",  xo,  [c3 * 1.0],  atol=1e-14)
check("y after z-rotation",  yo,  [s3 * 1.0],  atol=1e-14)
check("px after z-rotation", pxo, [c3 * 0.0],  atol=1e-14)
check("py after z-rotation", pyo, [s3 * 0.0],  atol=1e-14)
check("z unchanged",         zo,  [0.0],        atol=1e-14)


# ===========================================================================
# Test 6: g (δ) は座標変換で変化しない
# ===========================================================================
print("\n--- Test 6: δ (g) は変化しない ---")
el = Coord(name='Cd', CHI1=0.1, CHI2=0.05, CHI3=0.03, DX=0.001, DY=-0.001, DZ=0.0005, DIR=True)
g_vals = np.array([0.0, 0.01, -0.02, 0.005])
x_v = np.zeros_like(g_vals)
px_v = np.zeros_like(g_vals)
y_v = np.zeros_like(g_vals)
py_v = np.zeros_like(g_vals)
z_v = np.zeros_like(g_vals)
_, _, _, _, _, go = track_coord(el, x_v, px_v, y_v, py_v, z_v, g_vals)
check("g unchanged for all particles", go, g_vals)


# ===========================================================================
# Test 7: transfer_matrix の具体的な要素検証 (CHI3=pi/4)
# on-axis での解析値と一致することを確認。
# R = [[c3,-s3,0],[s3,c3,0],[0,0,1]]  (chi1=chi2=0)
# M[0,0]=r11-r13*r31/r33=c3-0=c3,  M[0,2]=-s3, M[1,1]=c3, M[1,3]=-s3
# M[2,0]=s3,  M[2,2]=c3, M[3,1]=s3,  M[3,3]=c3
# M[4,4]=M[5,5]=1, 他=0
# ===========================================================================
print("\n--- Test 7: transfer_matrix 要素検証 (CHI3=pi/4) ---")
from math import pi, cos, sin
angle = pi / 4
el = Coord(name='Cm', CHI3=angle, DIR=True)
M = el.transfer_matrix()
c3, s3 = cos(angle), sin(angle)
check("M[0,0] = c3",  [M[0, 0]], [c3],  atol=1e-14)
check("M[0,2] = -s3", [M[0, 2]], [-s3], atol=1e-14)
check("M[1,1] = c3",  [M[1, 1]], [c3],  atol=1e-14)
check("M[1,3] = -s3", [M[1, 3]], [-s3], atol=1e-14)
check("M[2,0] = s3",  [M[2, 0]], [s3],  atol=1e-14)
check("M[2,2] = c3",  [M[2, 2]], [c3],  atol=1e-14)
check("M[3,1] = s3",  [M[3, 1]], [s3],  atol=1e-14)
check("M[3,3] = c3",  [M[3, 3]], [c3],  atol=1e-14)
check("M[4,4] = 1",   [M[4, 4]], [1.0], atol=1e-14)
check("M[5,5] = 1",   [M[5, 5]], [1.0], atol=1e-14)
# 残り要素はすべて 0
off_diag_ok = (
    M[0,1]==0 and M[0,3]==0 and M[0,4]==0 and M[0,5]==0 and
    M[1,0]==0 and M[1,2]==0 and M[1,4]==0 and M[1,5]==0 and
    M[2,1]==0 and M[2,3]==0 and M[2,4]==0 and M[2,5]==0 and
    M[3,0]==0 and M[3,2]==0 and M[3,4]==0 and M[3,5]==0 and
    M[4,0]==0 and M[4,1]==0 and M[4,2]==0 and M[4,3]==0 and M[4,5]==0 and
    M[5,0]==0 and M[5,1]==0 and M[5,2]==0 and M[5,3]==0 and M[5,4]==0
)
check("off-diagonal zeros", [int(off_diag_ok)], [1])


# ===========================================================================
# Test 8: 多粒子ベクトル化動作確認
# ===========================================================================
print("\n--- Test 8: 多粒子ベクトル化 ---")
n = 50
rng = np.random.default_rng(42)
x_many  = rng.normal(0, 1e-3, n)
px_many = rng.normal(0, 5e-4, n)
y_many  = rng.normal(0, 1e-3, n)
py_many = rng.normal(0, 5e-4, n)
z_many  = rng.normal(0, 1e-4, n)
g_many  = rng.normal(0, 1e-2, n)

el = Coord(name='Cmulti', CHI1=0.1, CHI2=0.05, CHI3=0.02, DIR=True)
xo, pxo, yo, pyo, zo, go = track_coord(el, x_many, px_many, y_many, py_many, z_many, g_many)
check("output shape x",  [xo.shape[0]],  [n])
check("output shape px", [pxo.shape[0]], [n])
check("g unchanged (multi)", go, g_many)


# ===========================================================================
# Test 9: ELEMENT_TYPES に 'COORD' が登録されているか
# ===========================================================================
print("\n--- Test 9: ELEMENT_TYPES 登録確認 ---")
from pysad.elements import ELEMENT_TYPES
check("'COORD' in ELEMENT_TYPES", ['COORD' in ELEMENT_TYPES], [True])


# ===========================================================================
# Test 10: track_element() 経由のディスパッチ
# ===========================================================================
print("\n--- Test 10: track_element() ディスパッチ ---")
from pysad.tracking import track_element
el = Coord(name='Cd2', CHI3=0.05, DIR=True)
xo2, pxo2, yo2, pyo2, zo2, go2 = track_element(el, x0, px0, y0, py0, z0, g0)
xo_ref, pxo_ref, yo_ref, pyo_ref, zo_ref, go_ref = track_coord(el, x0, px0, y0, py0, z0, g0)
check("track_element x",  xo2,  xo_ref)
check("track_element px", pxo2, pxo_ref)
check("track_element y",  yo2,  yo_ref)
check("track_element py", pyo2, pyo_ref)
check("track_element z",  zo2,  zo_ref)


# ===========================================================================
# Summary
# ===========================================================================
print(f"\n{'='*50}")
print(f"  Results: {PASS} PASS, {FAIL} FAIL")
print(f"{'='*50}")
sys.exit(0 if FAIL == 0 else 1)
