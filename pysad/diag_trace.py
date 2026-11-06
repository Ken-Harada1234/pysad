"""
Diagnostic: trace transfer matrix element-by-element through ARCAH (MKS→MKSFA1)
and compare against SAD-expected matrices derived from marker Twiss data.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pysad.elements import Drift, Quadrupole, Bend, Sextupole, Marker, Monitor

# ─────────────────────────────────────────────────────────────────────────────
# Final-group element parameters (last definition in the .n file)
# ─────────────────────────────────────────────────────────────────────────────
DRID1_L = 2.669741905755901
DR05_L  = 0.05
DR40_L  = 0.40
DR18_L  = 0.18
DR45_L  = 0.45
DR20_L  = 0.20
DR26_L  = 0.26
DRID2_L = 0.791

QF1_L   = 0.42;   QF1_K1  =  0.6673251846840295
QD2H_L  = 0.16;   QD2H_K1 = -0.26151095666312235
QD3_L   = 0.32;   QD3_K1  = -0.3255825249794379
QF4_L   = 0.32;   QF4_K1  =  0.4912055094107217

SDAH_L  = 0.07
SFH_L   = 0.07

BID5_L  = 0.32;  BID5_A  = 0.016820289
BID4_L  = 0.32;  BID4_A  = 0.010019129
BID3H_L = 0.16;  BID3H_A = 0.003960072
BID2_L  = 0.32;  BID2_A  = 0.006469649
BID1_L  = 0.32;  BID1_A  = 0.00553479

# ─────────────────────────────────────────────────────────────────────────────
# SAD Twiss at markers (from the .n file)
# ─────────────────────────────────────────────────────────────────────────────
markers_sad = {
    'MKS':    {'BX':10.0,               'AX':0.0,                'BY':3.0,               'AY':0.0,
                'PSIX':0.0,             'PSIY':0.0},
    'MKBE':   {'BX':2.6395796774794973, 'AX':1.7441349553395291,  'BY':13.80022977033686,  'AY':-0.12051188458546798,
                'PSIX':0.45207394807955115, 'PSIY':0.8393662807518893},
    'MKS2':   {'BX':2.4961379130629306, 'AX':-2.3573671436328194, 'BY':11.86922990172713,  'AY':4.793902621384401,
                'PSIX':2.476129782556569,  'PSIY':0.9911015116980837},
    'MKSDA1': {'BX':2.8390412707249384, 'AX':-2.5412522515387135, 'BY':11.207983870513745, 'AY':4.652469253092507,
                'PSIX':2.502428099776169,  'PSIY':0.9971706311667013},
    'MKSFA1': {'BX':7.460481652558058,  'AX':-0.31509699385440815,'BY':6.083154671685557,  'AY':-0.12330304435387954,
                'PSIX':2.693915701151687,  'PSIY':1.1341150283341883},
    'MKC':    {'BX':2.9999999962355566, 'AX':0.0,                'BY':1.4999994619563635,  'AY':0.0,
                'PSIX':7.406304681648944,  'PSIY':2.704911317850623},
}

# ─────────────────────────────────────────────────────────────────────────────
# Build transfer matrix from Courant-Snyder parameters
# ─────────────────────────────────────────────────────────────────────────────
def M_from_twiss(bx1, ax1, bx2, ax2, dpsi_x, by1, ay1, by2, ay2, dpsi_y):
    """Reconstruct 2×2 transfer matrices in x and y from Courant-Snyder."""
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

def sad_M(mk1, mk2):
    """Expected transfer matrix M(mk1→mk2) from SAD Twiss data."""
    d1 = markers_sad[mk1]; d2 = markers_sad[mk2]
    dpsi_x = d2['PSIX'] - d1['PSIX']
    dpsi_y = d2['PSIY'] - d1['PSIY']
    return M_from_twiss(d1['BX'], d1['AX'], d2['BX'], d2['AX'], dpsi_x,
                        d1['BY'], d1['AY'], d2['BY'], d2['AY'], dpsi_y)

# ─────────────────────────────────────────────────────────────────────────────
# Print SAD-expected matrices
# ─────────────────────────────────────────────────────────────────────────────
print("=== SAD-expected transfer matrices (from Courant-Snyder Twiss) ===\n")
for mk in ['MKBE', 'MKS2', 'MKSDA1', 'MKSFA1']:
    Mxx, Myy = sad_M('MKS', mk)
    print(f"MKS → {mk}:")
    print(f"  Mx_SAD = [[{Mxx[0,0]:.5f}, {Mxx[0,1]:.5f}], [{Mxx[1,0]:.5f}, {Mxx[1,1]:.5f}]]")
    print(f"  My_SAD = [[{Myy[0,0]:.5f}, {Myy[0,1]:.5f}], [{Myy[1,0]:.5f}, {Myy[1,1]:.5f}]]")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Build ARCAH elements from MKS to MKSFA1, tracking at each marker
# ─────────────────────────────────────────────────────────────────────────────
sequence = [
    ('MKS',    Marker(name='MKS')),
    ('DRID1',  Drift(name='DRID1',  L=DRID1_L)),
    ('PM',     Monitor(name='PM',   L=0.0)),
    ('DR05',   Drift(name='DR05',   L=DR05_L)),
    ('QF1',    Quadrupole(name='QF1',  L=QF1_L,  K1=QF1_K1)),
    ('DR40',   Drift(name='DR40',   L=DR40_L)),
    ('QD2H',   Quadrupole(name='QD2H', L=QD2H_L, K1=QD2H_K1)),
    ('VD',     Drift(name='VD',     L=0.0)),   # ANGLE=0, ROTATE → drift
    ('QD2H_2', Quadrupole(name='QD2H_2', L=QD2H_L, K1=QD2H_K1)),
    ('DR18',   Drift(name='DR18',   L=DR18_L)),
    ('MKBE',   Marker(name='MKBE')),
    ('BID5',   Bend(name='BID5',   L=BID5_L, ANGLE=BID5_A)),
    ('BID4',   Bend(name='BID4',   L=BID4_L, ANGLE=BID4_A)),
    ('BID3H',  Bend(name='BID3H',  L=BID3H_L, ANGLE=BID3H_A)),
    ('HD',     Drift(name='HD',    L=0.0)),    # ANGLE=0 → drift
    ('BID3H_2',Bend(name='BID3H_2',L=BID3H_L, ANGLE=BID3H_A)),
    ('BID2',   Bend(name='BID2',   L=BID2_L, ANGLE=BID2_A)),
    ('BID1',   Bend(name='BID1',   L=BID1_L, ANGLE=BID1_A)),
    ('DR18_2', Drift(name='DR18_2',L=DR18_L)),
    ('QD3',    Quadrupole(name='QD3',  L=QD3_L,  K1=QD3_K1)),
    ('DR18_3', Drift(name='DR18_3',L=DR18_L)),
    ('MKS2',   Marker(name='MKS2')),
    ('SDAH',   Sextupole(name='SDAH', L=SDAH_L)),   # linear → drift
    ('MKSDA1', Marker(name='MKSDA1')),
    ('SDAH_2', Sextupole(name='SDAH_2', L=SDAH_L)),
    ('DR05_2', Drift(name='DR05_2',L=DR05_L)),
    ('PM2',    Monitor(name='PM2', L=0.0)),
    ('DR45',   Drift(name='DR45',  L=DR45_L)),
    ('QF4',    Quadrupole(name='QF4',  L=QF4_L,  K1=QF4_K1)),
    ('DR20',   Drift(name='DR20',  L=DR20_L)),
    ('SFH',    Sextupole(name='SFH',  L=SFH_L)),    # linear → drift
    ('VD2',    Drift(name='VD2',   L=0.0)),
    ('MKSFA1', Marker(name='MKSFA1')),
]

print("=== Element-by-element trace ===\n")
M_cum = np.eye(6)
s = 0.0
marker_M = {}

for name, elem in sequence:
    Mel = elem.transfer_matrix()
    M_cum = Mel @ M_cum
    s += elem.L
    if name in markers_sad:
        marker_M[name] = M_cum.copy()
        Mx = M_cum[0:2, 0:2]
        My = M_cum[2:4, 2:4]
        print(f"  {name:10s} (s={s:.4f} m):")
        print(f"    Mx = [[{Mx[0,0]:.5f}, {Mx[0,1]:.5f}], [{Mx[1,0]:.5f}, {Mx[1,1]:.5f}]]")
        print(f"    My = [[{My[0,0]:.5f}, {My[0,1]:.5f}], [{My[1,0]:.5f}, {My[1,1]:.5f}]]")

# ─────────────────────────────────────────────────────────────────────────────
# Compare each marker: ours vs SAD
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== Ratio analysis: ours/SAD for cumulative M (MKS→marker) ===\n")
for mk in ['MKBE', 'MKS2', 'MKSDA1', 'MKSFA1']:
    if mk not in marker_M:
        continue
    Mx_sad, My_sad = sad_M('MKS', mk)
    Mx_ours = marker_M[mk][0:2, 0:2]
    My_ours = marker_M[mk][2:4, 2:4]

    print(f"=== MKS → {mk} ===")
    for i in range(2):
        for j in range(2):
            if abs(Mx_sad[i,j]) > 0.005:
                ratio = Mx_ours[i,j] / Mx_sad[i,j]
                diff = Mx_ours[i,j] - Mx_sad[i,j]
                print(f"  Mx[{i},{j}]: ours={Mx_ours[i,j]:9.5f}  SAD={Mx_sad[i,j]:9.5f}  ratio={ratio:.4f}  diff={diff:+.5f}")
    for i in range(2):
        for j in range(2):
            if abs(My_sad[i,j]) > 0.005:
                ratio = My_ours[i,j] / My_sad[i,j]
                diff = My_ours[i,j] - My_sad[i,j]
                print(f"  My[{i},{j}]: ours={My_ours[i,j]:9.5f}  SAD={My_sad[i,j]:9.5f}  ratio={ratio:.4f}  diff={diff:+.5f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sub-segment analysis: check each marker interval
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Sub-segment ratios (each marker interval) ===\n")
for (mk1, mk2) in [('MKS','MKBE'), ('MKBE','MKS2'), ('MKS2','MKSDA1'), ('MKSDA1','MKSFA1')]:
    M1 = marker_M.get(mk1, np.eye(6))
    M2 = marker_M.get(mk2, np.eye(6))
    M_seg = M2 @ np.linalg.inv(M1)
    Mx_ours = M_seg[0:2, 0:2]
    My_ours = M_seg[2:4, 2:4]
    Mx_sad, My_sad = sad_M(mk1, mk2)

    print(f"=== {mk1} → {mk2} ===")
    for i in range(2):
        for j in range(2):
            if abs(Mx_sad[i,j]) > 0.005:
                ratio = Mx_ours[i,j] / Mx_sad[i,j]
                print(f"  Mx[{i},{j}]: ours={Mx_ours[i,j]:9.5f}  SAD={Mx_sad[i,j]:9.5f}  ratio={ratio:.4f}")
    for i in range(2):
        for j in range(2):
            if abs(My_sad[i,j]) > 0.005:
                ratio = My_ours[i,j] / My_sad[i,j]
                print(f"  My[{i},{j}]: ours={My_ours[i,j]:9.5f}  SAD={My_sad[i,j]:9.5f}  ratio={ratio:.4f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Fine-grained trace: after each element within MKS→MKBE segment
# Check where x-scale factor first deviates from 1.0
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Fine trace: MKS→MKBE element by element ===\n")
print("Phase advances Δψ_x and Δψ_y at each step (from ours vs SAD expected at MKBE)\n")

Mx_mkbe_sad, My_mkbe_sad = sad_M('MKS', 'MKBE')

# We'll track M and at each step compute what the 'beta function would be' by M[0,1] = sqrt(b1*b2)*sin(psi)
# More directly: track M_cum[0,0], M_cum[0,1] and compare at MKBE
bx1, ax1 = markers_sad['MKS']['BX'], markers_sad['MKS']['AX']
by1, ay1 = markers_sad['MKS']['BY'], markers_sad['MKS']['AY']

print(f"Initial: BX={bx1}, AX={ax1}, BY={by1}, AY={ay1}")
print()

# Re-run MKS→MKBE sub-sequence
mks_to_mkbe = [
    ('DRID1',  Drift(name='DRID1',  L=DRID1_L)),
    ('PM',     Monitor(name='PM',   L=0.0)),
    ('DR05',   Drift(name='DR05',   L=DR05_L)),
    ('QF1',    Quadrupole(name='QF1',  L=QF1_L,  K1=QF1_K1)),
    ('DR40',   Drift(name='DR40',   L=DR40_L)),
    ('QD2H',   Quadrupole(name='QD2H', L=QD2H_L, K1=QD2H_K1)),
    ('VD',     Drift(name='VD',     L=0.0)),
    ('QD2H_2', Quadrupole(name='QD2H_2', L=QD2H_L, K1=QD2H_K1)),
    ('DR18',   Drift(name='DR18',   L=DR18_L)),
]

M_cur = np.eye(6)
s = 0.0
for name, elem in mks_to_mkbe:
    M_cur = elem.transfer_matrix() @ M_cur
    s += elem.L
    Mx = M_cur[0:2, 0:2]
    My = M_cur[2:4, 2:4]
    # Compute current beta from M[0,1] = sqrt(b1*b2)*sin(psi) and M[0,0] = sqrt(b2/b1)*(cos+ax1*sin)
    # Phase: psi_x = atan2(Mx[0,1]/sqrt(bx1), Mx[0,0] - ax1*Mx[0,1]/bx1)
    # Actually easier: just show M[0,0] and M[0,1]
    print(f"  After {name:8s} (s={s:6.3f}m): Mx=[[{Mx[0,0]:7.4f},{Mx[0,1]:7.4f}],[{Mx[1,0]:7.4f},{Mx[1,1]:7.4f}]]"
          f"  My=[[{My[0,0]:7.4f},{My[0,1]:7.4f}],[{My[1,0]:7.4f},{My[1,1]:7.4f}]]")

print(f"\n  Expected at MKBE:")
print(f"    Mx_SAD = [[{Mx_mkbe_sad[0,0]:.4f},{Mx_mkbe_sad[0,1]:.4f}],[{Mx_mkbe_sad[1,0]:.4f},{Mx_mkbe_sad[1,1]:.4f}]]")
print(f"    My_SAD = [[{My_mkbe_sad[0,0]:.4f},{My_mkbe_sad[0,1]:.4f}],[{My_mkbe_sad[1,0]:.4f},{My_mkbe_sad[1,1]:.4f}]]")
print(f"\n  Ours at MKBE:")
Mx_mkbe_ours = marker_M['MKBE'][0:2, 0:2]
My_mkbe_ours = marker_M['MKBE'][2:4, 2:4]
print(f"    Mx_ours= [[{Mx_mkbe_ours[0,0]:.4f},{Mx_mkbe_ours[0,1]:.4f}],[{Mx_mkbe_ours[1,0]:.4f},{Mx_mkbe_ours[1,1]:.4f}]]")
print(f"    My_ours= [[{My_mkbe_ours[0,0]:.4f},{My_mkbe_ours[0,1]:.4f}],[{My_mkbe_ours[1,0]:.4f},{My_mkbe_ours[1,1]:.4f}]]")

print("\n\n=== K1 convention check for QF1 ===")
print(f"QF1: L={QF1_L}, K1_stored={QF1_K1}")
print(f"  K1_eff (per-unit-length) = K1_stored / L = {QF1_K1/QF1_L:.6f} m^-2")
print(f"  kx = sqrt(K1_eff) = {np.sqrt(QF1_K1/QF1_L):.6f} m^-1")
print(f"  kx*L = {np.sqrt(QF1_K1/QF1_L)*QF1_L:.6f} rad")
print(f"  QF1 Mx = [[cos(kx*L), sin(kx*L)/kx], [-kx*sin(kx*L), cos(kx*L)]]")
kx = np.sqrt(QF1_K1/QF1_L)
L = QF1_L
print(f"         = [[{np.cos(kx*L):.5f}, {np.sin(kx*L)/kx:.5f}], [{-kx*np.sin(kx*L):.5f}, {np.cos(kx*L):.5f}]]")

print(f"\n  If K1 were NOT divided by L (wrong convention):")
k1_wrong = QF1_K1
kx_wrong = np.sqrt(abs(k1_wrong)) if k1_wrong > 0 else None
print(f"  K1_eff = {k1_wrong:.6f} → kx = {np.sqrt(k1_wrong):.6f}")
print(f"  kx*L = {np.sqrt(k1_wrong)*L:.6f} rad")
