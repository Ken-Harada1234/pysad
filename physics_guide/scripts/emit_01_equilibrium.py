# -*- coding: utf-8 -*-
"""平衡エミッタンス (14_emittance.md): CALC6D / tecalc の結果 (6 次元の Σ = MΣM^T + Δ) を、放射積分の標準の式と比べる。

リング: 8 セルの FODO (3 GeV、扇形偏向 16 台、QF = 0.25・QD = -0.30、周長 41.6 m)、CAVI (VOLT = 3 MV、HARM = 64) 1 台 (ラインの先頭: 励起の行列 Δ は RF のあとの放射だけ)。
(a) calc6d の結果: チューン、減衰率、エミッタンス、エネルギー広がり、1 周のエネルギー損失。
(b) 放射積分の式: U0 = C_γ E⁴ I2/2π、J_x = 1 - I4/I2、J_z = 2 + I4/I2、ε_x = C_q γ² I5/(J_x I2)、σ_δ² = C_q γ² I3/(J_z I2)、α_p = I1/C、シンクロトロンチューン。
(c) 励起の行列 Δ = Σ - MΣM^T の δδ 成分と、偏向磁石 16 台の N_γ <u²> の和 (13 章)。Δ の半正定値性。
(d) Σ の固有構造: 正規座標への変換で Σ が対角 (エミッタンス) になるか、ほか。
"""
import io
import math
import numpy as np
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad.calc6d import calc6d

print('pysad:', C.pysad_src())
mc2 = 0.51099895e6
re = 2.8179403262e-15
alpha = 1 / 137.035999084
lam = re / alpha
src = '''MOMENTUM = 3GEV;
DRIFT D1=(L=0.8) D2=(L=0.3);
BEND B=(L=1.2 ANGLE=%.15g);
QUAD QF=(L=0.3 K1=0.25) QD=(L=0.3 K1=-0.30);
CAVI RF=(L=0 VOLT=3e6 HARM=64);
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=(ST RF CELL CELL CELL CELL CELL CELL CELL CELL);
FFS;
USE RING;
''' % (2 * math.pi / 16)
ffs, _ = run_sad_string(src, out=io.StringIO())
lat = ffs.lat
circ = sum(e.L for e in lat.elements)
res = calc6d(lat, ffs.momentum(), rfsw=True, radcod=True)
em = res.emit
E = 3e9
g = E / mc2
print('周長 %.3f m、チューン (νx, νy, νs) = (%.4f, %.4f, %.4f) (符号は固有値の向き)' % (circ, *[abs(t) for t in em['tunes']]))

print('\n== (a) calc6d (tecalc) の結果 ==')
dr = -np.asarray(em['damping_rate'])
print('1 周のエネルギー損失 U0 = %.4f MeV' % (em['energy_loss_eV'] / 1e6))
print('減衰率 (1 周あたり) (x, y, z) = (%.5e, %.5e, %.5e)、和 = %.5e' % (*dr, dr.sum()))
print('Jx, Jy, Jz = 減衰率 x 4 / 和: (%.4f, %.4f, %.4f)' % tuple(dr * 4 / dr.sum()))
print('エミッタンス εx = %.4e m rad、εy = %.2e、σδ = %.4e、σz = %.4e m、α_p = %.5f' % (
    em['emittance_x'], em['emittance_y'], em['momentum_spread'], em['bunch_length'], em['momentum_compaction']))

print('\n== (b) 放射積分の式 ==')
ri = lat.radiation_integrals()
I1, I2, I3, I4, I5 = [ri[k] for k in ('I1', 'I2', 'I3', 'I4', 'I5')]
Cg = 4 * math.pi / 3 * re / (mc2 / 1e9) ** 3
Cq = 55 / (32 * math.sqrt(3)) * lam
U0 = Cg * 3.0**4 * I2 / (2 * math.pi) * 1e9
Jx, Jz = 1 - I4 / I2, 2 + I4 / I2
print('I1..I5 = %.4f %.4f %.4f %.4f %.4f' % (I1, I2, I3, I4, I5))
print('U0 = C_γ E⁴ I2/(2π) = %.4f MeV (calc6d %.4f)' % (U0 / 1e6, em['energy_loss_eV'] / 1e6))
print('Jx = 1 - I4/I2 = %.4f、Jz = 2 + I4/I2 = %.4f、(Jx, Jy, Jz) の和 = %.4f' % (Jx, Jz, Jx + 1 + Jz))
print('εx = C_q γ² I5/(Jx I2) = %.4e m rad (calc6d %.4e、比 %.4f)' % (Cq * g**2 * I5 / (Jx * I2), em['emittance_x'], em['emittance_x'] / (Cq * g**2 * I5 / (Jx * I2))))
sd = math.sqrt(Cq * g**2 * I3 / (Jz * I2))
print('σδ = √(C_q γ² I3/(Jz I2)) = %.4e (calc6d %.4e、比 %.4f)' % (sd, em['momentum_spread'], em['momentum_spread'] / sd))
ap = I1 / circ
print('α_p = I1/C = %.5f (calc6d %.5f)' % (ap, em['momentum_compaction']))
V, h = 3e6, 64
phis = math.asin(U0 / V)
nus = math.sqrt(h * ap * V * math.cos(phis) / (2 * math.pi * E))
print('同期位相 sin φs = U0/V = %.4f (φs = %.4f rad)。シンクロトロンチューン νs = √(h α_p V cos φs/(2π E)) = %.4f (calc6d %.4f。符号の向きにより α_p の符号は逆。大きさで比べる)' % (
    U0 / V, phis, nus, abs(em['tunes'][2])))

print('\n== (c) 励起の行列 Δ = Σ - MΣM^T ==')
Sg = np.asarray(em['beam_matrix'])
M = np.asarray(em['one_turn_matrix'])
Dl = Sg - M @ Sg @ M.T
print('Σ の対称性 |Σ - Σ^T| = %.1e、Δ の対称性 %.1e、Δ の固有値の最小・最大 = %.3e、%.3e' % (
    np.abs(Sg - Sg.T).max(), np.abs(Dl - Dl.T).max(), np.linalg.eigvalsh((Dl + Dl.T) / 2).min(), np.linalg.eigvalsh((Dl + Dl.T) / 2).max()))
th = 2 * math.pi / 16
rho = 1.2 / th
Nph = 5 / (2 * math.sqrt(3)) * alpha * g * th
uc = 1.5 * lam * g**3 / rho * mc2 / E
dee = 16 * Nph * 11 / 27 * uc**2
print('Δδδ = %.4e、16 台の N_γ <u²> の和 (13 章) = %.4e (比 %.4f)' % (Dl[5, 5], dee, Dl[5, 5] / dee))
print('Σδδ = %.4e = σδ² (σδ = %.4e)' % (Sg[5, 5], math.sqrt(Sg[5, 5])))

print('\n== (d) 正規座標 ==')
r = np.asarray(em['normal_r'])
ri_ = np.asarray(em['normal_ri'])
Sn = ri_ @ Sg @ ri_.T
print('正規座標の Σ_n (対角: 各モードのエミッタンス)。対角成分:', np.array2string(np.diag(Sn), precision=4))
print('モード内の (X,Px) の対角の和/2: εx = %.4e、εy = %.2e、εz = %.4e' % ((Sn[0, 0] + Sn[1, 1]) / 2, (Sn[2, 2] + Sn[3, 3]) / 2, (Sn[4, 4] + Sn[5, 5]) / 2))
print('非対角成分の最大 (モード間の相関、0 なら完全に分離) = %.2e (対角の最大 %.2e)' % (
    max(np.abs(Sn[0:2, 2:]).max(), np.abs(Sn[2:4, 4:]).max(), np.abs(Sn[0:2, 4:]).max()), np.abs(np.diag(Sn)).max()))
print('R^{-1} M R の x ブロック (回転 + 減衰): 対角 %.6f、非対角 %.6f (期待: e^{-damp} cos μ, sin μ)' % (
    (ri_ @ M @ r)[0, 0], (ri_ @ M @ r)[0, 1]))

print('\n== (e) 各モードの平衡: ε_k = D_k/(1 - e^{-2 a_k}) ≈ D_k/(2 a_k) (D_k = (Δ_n,XX + Δ_n,PP)/2) ==')
Dn = ri_ @ Dl @ ri_.T
a_k = dr                                     # 1 周あたりの減衰 (振幅)
for k, nm in ((0, 'x'), (2, 'y'), (4, 'z')):
    trD = 0.5 * (Dn[k, k] + Dn[k + 1, k + 1])
    pred = trD / (1 - math.exp(-2 * a_k[k // 2]))
    print('モード %s: Δ_n の (X,Px) 平均 = %.4e、減衰 a = %.5e、予想 ε = %.4e、近似 D/(2a) = %.4e、Σ_n の ε = %.4e' % (
        nm, trD, a_k[k // 2], pred, trD / (2 * a_k[k // 2]), 0.5 * (Sn[k, k] + Sn[k + 1, k + 1])))
print('減衰率の和 Jx + Jy + Jz = %.4f (Robinson の定理: 4)' % (dr.sum() * 4 / dr.sum()))
