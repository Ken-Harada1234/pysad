"""
Numerically solve for what QF1 and QD2H K1 values would reproduce
the SAD-expected transfer matrix at MKBE.

This will tell us: is there a systematic scaling factor between our K1 convention
and what SAD uses?
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import fsolve, minimize

# SAD Twiss
MKS  = dict(BX=10.0, AX=0.0, BY=3.0, AY=0.0, PSIX=0.0, PSIY=0.0)
MKBE = dict(BX=2.6395796774794973, AX=1.7441349553395291,
            BY=13.80022977033686, AY=-0.12051188458546798,
            PSIX=0.45207394807955115, PSIY=0.8393662807518893)

def M_from_twiss(bx1, ax1, bx2, ax2, dpsi_x, by1, ay1, by2, ay2, dpsi_y):
    cx, sx = np.cos(dpsi_x), np.sin(dpsi_x)
    Mxx = np.array([
        [np.sqrt(bx2/bx1)*(cx + ax1*sx),  np.sqrt(bx1*bx2)*sx],
        [-(1+ax1*ax2)*sx/np.sqrt(bx1*bx2) + (ax1-ax2)*cx/np.sqrt(bx1*bx2),
          np.sqrt(bx1/bx2)*(cx - ax2*sx)]
    ])
    cy, sy = np.cos(dpsi_y), np.sin(dpsi_y)
    Myy = np.array([
        [np.sqrt(by2/by1)*(cy + ay1*sy),  np.sqrt(by1*by2)*sy],
        [-(1+ay1*ay2)*sy/np.sqrt(by1*by2) + (ay1-ay2)*cy/np.sqrt(by1*by2),
          np.sqrt(by1/by2)*(cy - ay2*sy)]
    ])
    return Mxx, Myy

Mx_sad, My_sad = M_from_twiss(
    MKS['BX'], MKS['AX'], MKBE['BX'], MKBE['AX'], MKBE['PSIX'],
    MKS['BY'], MKS['AY'], MKBE['BY'], MKBE['AY'], MKBE['PSIY']
)
print("SAD-expected at MKBE:")
print(f"  Mx_sad = {Mx_sad}")
print(f"  My_sad = {My_sad}")

# ─────────────────────────────────────────────────────────────────────────────
# Build MKS→MKBE transfer matrix as a function of (k1_QF1_eff, k1_QD2H_eff)
# where k_eff is per-unit-length [m^-2]
# ─────────────────────────────────────────────────────────────────────────────
DRID1_L = 2.669741905755901
DR05_L  = 0.05
DR40_L  = 0.40
DR18_L  = 0.18
QF1_L   = 0.42
QD2H_L  = 0.16

def quad_matrix_x_focus(k_eff, L):
    """x-focusing quad (k_eff > 0)"""
    kx = np.sqrt(k_eff)
    return np.array([[np.cos(kx*L), np.sin(kx*L)/kx],
                     [-kx*np.sin(kx*L), np.cos(kx*L)]])

def quad_matrix_y_defocus(k_eff, L):
    """y-defocusing quad (k_eff > 0 → x-focusing → y-defocusing)"""
    ky = np.sqrt(k_eff)
    return np.array([[np.cosh(ky*L), np.sinh(ky*L)/ky],
                     [ky*np.sinh(ky*L), np.cosh(ky*L)]])

def quad_matrix_x_defocus(k_eff, L):
    """x-defocusing quad (k_eff > 0 → this is |K1|/L)"""
    kx = np.sqrt(k_eff)
    return np.array([[np.cosh(kx*L), np.sinh(kx*L)/kx],
                     [kx*np.sinh(kx*L), np.cosh(kx*L)]])

def quad_matrix_y_focus(k_eff, L):
    """y-focusing quad (k_eff > 0 → this is |K1|/L)"""
    ky = np.sqrt(k_eff)
    return np.array([[np.cos(ky*L), np.sin(ky*L)/ky],
                     [-ky*np.sin(ky*L), np.cos(ky*L)]])

def drift_2x2(L):
    return np.array([[1.0, L], [0.0, 1.0]])

def compute_MKS_to_MKBE(k1_QF1_eff, k1_QD2H_eff):
    """
    Compute 2x2 x and y transfer matrices from MKS to MKBE.
    k1_QF1_eff  > 0: x-focusing, y-defocusing
    k1_QD2H_eff > 0: y-focusing, x-defocusing (use |k1_QD2H|)
    """
    # x-plane
    Mx = drift_2x2(DRID1_L + DR05_L)
    Mx = quad_matrix_x_focus(k1_QF1_eff, QF1_L) @ Mx
    Mx = drift_2x2(DR40_L) @ Mx
    Mx = quad_matrix_x_defocus(k1_QD2H_eff, QD2H_L) @ Mx  # QD2H is x-defocusing
    # VD: L=0 → skip
    Mx = quad_matrix_x_defocus(k1_QD2H_eff, QD2H_L) @ Mx  # second QD2H
    Mx = drift_2x2(DR18_L) @ Mx

    # y-plane
    My = drift_2x2(DRID1_L + DR05_L)
    My = quad_matrix_y_defocus(k1_QF1_eff, QF1_L) @ My   # QF1 y-defocusing
    My = drift_2x2(DR40_L) @ My
    My = quad_matrix_y_focus(k1_QD2H_eff, QD2H_L) @ My   # QD2H y-focusing
    My = quad_matrix_y_focus(k1_QD2H_eff, QD2H_L) @ My   # second QD2H
    My = drift_2x2(DR18_L) @ My

    return Mx, My

# ─────────────────────────────────────────────────────────────────────────────
# Verify our current values
# ─────────────────────────────────────────────────────────────────────────────
# Our convention: k1_eff = K1/L
K1_QF1  =  0.6673251846840295
K1_QD2H = -0.26151095666312235

k1_QF1_eff_ours  = K1_QF1 / QF1_L
k1_QD2H_eff_ours = abs(K1_QD2H) / QD2H_L

print(f"\nOur k1_eff values:")
print(f"  QF1:  K1={K1_QF1:.6f}, L={QF1_L}, k1_eff = K1/L = {k1_QF1_eff_ours:.6f} m^-2")
print(f"  QD2H: K1={K1_QD2H:.6f}, L={QD2H_L}, k1_eff = |K1|/L = {k1_QD2H_eff_ours:.6f} m^-2")

Mx_ours, My_ours = compute_MKS_to_MKBE(k1_QF1_eff_ours, k1_QD2H_eff_ours)
print(f"\nOurs at MKBE:")
print(f"  Mx = {Mx_ours}")
print(f"  My = {My_ours}")

# ─────────────────────────────────────────────────────────────────────────────
# Alternative: if K1 is stored as m^-2 directly (no division by L)
# ─────────────────────────────────────────────────────────────────────────────
k1_QF1_eff_alt  = abs(K1_QF1)    # no division by L
k1_QD2H_eff_alt = abs(K1_QD2H)   # no division by L

print(f"\n'Alt' k1_eff values (NO division by L):")
print(f"  QF1:  k1_eff = {k1_QF1_eff_alt:.6f} m^-2")
print(f"  QD2H: k1_eff = {k1_QD2H_eff_alt:.6f} m^-2")

Mx_alt, My_alt = compute_MKS_to_MKBE(k1_QF1_eff_alt, k1_QD2H_eff_alt)
print(f"\n'Alt' at MKBE:")
print(f"  Mx = {Mx_alt}")
print(f"  My = {My_alt}")

# ─────────────────────────────────────────────────────────────────────────────
# Numerical fit: find k1_eff values that reproduce SAD matrix at MKBE
# ─────────────────────────────────────────────────────────────────────────────
def residuals(params):
    k1_qf1, k1_qd2h = params[0], params[1]
    if k1_qf1 <= 0 or k1_qd2h <= 0:
        return [1e6] * 8
    Mx, My = compute_MKS_to_MKBE(k1_qf1, k1_qd2h)
    res = []
    for i in range(2):
        for j in range(2):
            res.append(Mx[i,j] - Mx_sad[i,j])
            res.append(My[i,j] - My_sad[i,j])
    return res

x0 = [k1_QF1_eff_ours, k1_QD2H_eff_ours]
sol = fsolve(residuals, x0, full_output=True)
k1_fit = sol[0]
info = sol[1]
print(f"\nNumerical fit (minimize both Mx and My residuals at MKBE):")
print(f"  k1_QF1_eff_fit  = {k1_fit[0]:.8f} m^-2")
print(f"  k1_QD2H_eff_fit = {k1_fit[1]:.8f} m^-2")
print(f"  Ratio QF1:  fit/ours = {k1_fit[0]/k1_QF1_eff_ours:.5f}")
print(f"  Ratio QD2H: fit/ours = {k1_fit[1]/k1_QD2H_eff_ours:.5f}")
print(f"  Ratio QF1:  fit vs K1 (no L) = {k1_fit[0]/k1_QF1_eff_alt:.5f}")
print(f"  Ratio QD2H: fit vs |K1| (no L) = {k1_fit[1]/k1_QD2H_eff_alt:.5f}")

# Check whether there's a common factor
Mx_fit, My_fit = compute_MKS_to_MKBE(k1_fit[0], k1_fit[1])
print(f"\nFit at MKBE:")
print(f"  Mx_fit = {Mx_fit}")
print(f"  My_fit = {My_fit}")
print(f"\nResiduals at fit:")
for i in range(2):
    for j in range(2):
        print(f"  Mx[{i},{j}]: fit={Mx_fit[i,j]:.6f} sad={Mx_sad[i,j]:.6f} diff={Mx_fit[i,j]-Mx_sad[i,j]:.6e}")
        print(f"  My[{i},{j}]: fit={My_fit[i,j]:.6f} sad={My_sad[i,j]:.6f} diff={My_fit[i,j]-My_sad[i,j]:.6e}")

# ─────────────────────────────────────────────────────────────────────────────
# More general: fit only x-plane to see if k1_QF1 changes
# ─────────────────────────────────────────────────────────────────────────────
def residuals_x(params):
    k1_qf1, k1_qd2h = params[0], params[1]
    if k1_qf1 <= 0 or k1_qd2h <= 0:
        return [1e6] * 4
    Mx, _ = compute_MKS_to_MKBE(k1_qf1, k1_qd2h)
    return [Mx[i,j] - Mx_sad[i,j] for i in range(2) for j in range(2)]

sol_x = fsolve(residuals_x, x0)
print(f"\nFit using x-plane only:")
print(f"  k1_QF1_eff  = {sol_x[0]:.8f} (ratio to ours: {sol_x[0]/k1_QF1_eff_ours:.5f})")
print(f"  k1_QD2H_eff = {sol_x[1]:.8f} (ratio to ours: {sol_x[1]/k1_QD2H_eff_ours:.5f})")

def residuals_y(params):
    k1_qf1, k1_qd2h = params[0], params[1]
    if k1_qf1 <= 0 or k1_qd2h <= 0:
        return [1e6] * 4
    _, My = compute_MKS_to_MKBE(k1_qf1, k1_qd2h)
    return [My[i,j] - My_sad[i,j] for i in range(2) for j in range(2)]

sol_y = fsolve(residuals_y, x0)
print(f"\nFit using y-plane only:")
print(f"  k1_QF1_eff  = {sol_y[0]:.8f} (ratio to ours: {sol_y[0]/k1_QF1_eff_ours:.5f})")
print(f"  k1_QD2H_eff = {sol_y[1]:.8f} (ratio to ours: {sol_y[1]/k1_QD2H_eff_ours:.5f})")

# ─────────────────────────────────────────────────────────────────────────────
# Check: what if there's a common multiplier to K1 for ALL elements?
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Scanning common scale factor for k1_eff ===")
print("If all K1/L are multiplied by factor f, what f minimizes the error at MKBE?\n")

best_f = None
best_err = 1e10
for f in np.linspace(0.5, 2.0, 200):
    Mx, My = compute_MKS_to_MKBE(f*k1_QF1_eff_ours, f*k1_QD2H_eff_ours)
    err = np.sum((Mx - Mx_sad)**2) + np.sum((My - My_sad)**2)
    if err < best_err:
        best_err = err
        best_f = f

print(f"  Best common scale factor: f = {best_f:.4f}")
Mx_best, My_best = compute_MKS_to_MKBE(best_f*k1_QF1_eff_ours, best_f*k1_QD2H_eff_ours)
print(f"  Mx at best f: {Mx_best}")
print(f"  My at best f: {My_best}")
print(f"  This corresponds to:")
print(f"    QF1  k1_eff_new = {best_f*k1_QF1_eff_ours:.6f} m^-2 (vs K1={K1_QF1:.6f}, K1/L={k1_QF1_eff_ours:.6f})")
print(f"    QD2H k1_eff_new = {best_f*k1_QD2H_eff_ours:.6f} m^-2 (vs |K1|={abs(K1_QD2H):.6f}, |K1|/L={k1_QD2H_eff_ours:.6f})")

# ─────────────────────────────────────────────────────────────────────────────
# Twiss propagation from M
# ─────────────────────────────────────────────────────────────────────────────
def beta_from_M(M2x2, bx1, ax1):
    """Compute beta at end point from transfer matrix and initial Twiss."""
    b = M2x2[0,0]**2 * bx1 - 2*M2x2[0,0]*M2x2[0,1]*ax1 + M2x2[0,1]**2 * (1+ax1**2)/bx1
    a = -(M2x2[0,0]*M2x2[1,0]*bx1 - (M2x2[0,0]*M2x2[1,1]+M2x2[0,1]*M2x2[1,0])*ax1 + M2x2[0,1]*M2x2[1,1]*(1+ax1**2)/bx1)
    return b, a

bx_mkbe_ours, ax_mkbe_ours = beta_from_M(Mx_ours, 10.0, 0.0)
by_mkbe_ours, ay_mkbe_ours = beta_from_M(My_ours, 3.0, 0.0)

print(f"\nBeta at MKBE from our M:")
print(f"  BX_ours = {bx_mkbe_ours:.4f} (SAD: {MKBE['BX']:.4f})")
print(f"  AX_ours = {ax_mkbe_ours:.4f} (SAD: {MKBE['AX']:.4f})")
print(f"  BY_ours = {by_mkbe_ours:.4f} (SAD: {MKBE['BY']:.4f})")
print(f"  AY_ours = {ay_mkbe_ours:.4f} (SAD: {MKBE['AY']:.4f})")
