# -*- coding: utf-8 -*-
"""六極の組ごとの解析量の表 (R7): 1 次の RDT、δ の項、4 次の RDT (正規形)、ADTS、クロマティシティ (解析式と CALC)、目的関数の内訳。
使い方 (どこからでも): py -3.12 sextupole_rdt_table.py [--cand results/candidates.json] [--no-nform] [--ring RING]
出力: results/rdt_table.json、results/rdt_table.md、results/rdt_table.out.txt、fig_rdt_walk.png (要素ごとの寄与の累積)、fig_rdt_sets.png (組ごとの比較)
条件: 運動量 5 GeV (SAD の K は規格化済みで、2.5 GeV の版も光学は同じ)、リング先頭 (長直線部の中心) が位相の基準、偏向磁石の 3 次の項・八極 (K3 = 0)・運動学・フリンジを入れた
pysad の解析式 (`resdrive.resonance_terms`・`rdt2nd.rdt4`・`resdrive.adts`、多項式モデル `sxtopt.SxtModel` が 1e-12 で一致)。クロマティシティの CALC は 4 次元の閉軌道の δ 依存 (order 3)。
nform: 追跡の 1 周の写像の数値の正規形 (nform.adts_at_delta・adts_delta_derivative(method='nform')) で、全部の組の ADTS (δ=0、解析式の検証) と ∂(ADTS)/∂δ を出す (1 組に数秒)。
解析式の δ 微分は yy などで外れる (docs/49) ので使わない。
"""
import argparse
import io
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import nform
from pysad import rdt2nd as Q
from pysad import resdrive as R
from pysad import sxtopt as S

ap = argparse.ArgumentParser()
ap.add_argument('--cand', default=os.path.join(C.OUT, 'candidates.json'))
ap.add_argument('--no-nform', action='store_true')
ap.add_argument('--ring', default='RING')
ap.add_argument('--no-fig', action='store_true')
args = ap.parse_args()
LOG = io.StringIO()


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')


CAND = json.load(open(args.cand, encoding='utf-8'))
src, lat = C.load_lattice(args.ring)
m = S.SxtModel(lat, C.P0)
assert CAND['family_names'] == m.names
obj = S.SxtObjective(m)
kf = np.array(CAND['candidates']['final']['k'])
f_final = obj.value(kf)
rows_all = {}
P('リング %s (%d 素子)、チューン (%.4f, %.4f)、六極 %d ファミリー' % (args.ring, len(lat.elements), m.qx, m.qy, m.F))
t0 = time.time()
for name, c in CAND['candidates'].items():
    k = np.array(c['k'])
    h1 = m.h1(k)
    h4 = m.h4nf_fast(k)
    a = m.adts(k)
    lat2 = S.with_k2l(lat, m.names, k)
    ch = R.chromaticity_numeric(lat2, order=3, dmax=2e-3)
    p = obj.parts(k)
    nd = np.abs(h4[obj.i4])
    row = dict(name=name, group=c['group'], desc=c['desc'], k_max=float(np.max(np.abs(k))),
               xi_an=[float(v) for v in m.xi(k)], xi_calc=list(ch['d1']), xi2_calc=list(ch['d2']),
               h3=[abs(h1[t]) for t in S.CUBIC5], h3_names=list(S.CUBIC5), h3_sum=float(sum(abs(h1[t]) for t in S.CUBIC5)),
               hchrom=[abs(h1[t]) for t in S.CHROM3], h4_nondiag_sum=float(np.sum(nd)), h4_nondiag_max=float(np.max(nd)),
               adts=[float(a['xx']), float(a['xy']), float(a['yy'])], adts_norm=float(np.sqrt(sum(float(v) ** 2 for v in a.values()))),
               parts=p, obj=float(obj.value(k)), obj_rel=float(obj.value(k) / f_final))
    # 基準の振幅 J_ref = 1e-6 m rad (Jx = Jy) の 1 周あたりのチューンのずれ
    row['dnu_ref'] = [float((a['xx'] + a['xy']) * obj.J), float((a['xy'] + a['yy']) * obj.J)]
    rows_all[name] = row
P('解析量の計算 %.1f s' % (time.time() - t0))

# 追跡の写像の数値の正規形 (nform): δ = 0 の ADTS (解析式の検証) と ∂(ADTS)/∂δ (δ = ±3e-3 の中心差分)。--no-nform で省く (1 組に数秒)
if not args.no_nform:
    t1 = time.time()
    for name in rows_all:
        k = np.array(CAND['candidates'][name]['k'])
        lat2 = S.with_k2l(lat, m.names, k)
        a0 = nform.adts_at_delta(lat2, C.P0, 0.0, amp3=1e-5)[0]
        d = nform.adts_delta_derivative(lat2, momentum_eV=C.P0, delta0=3e-3, method='nform', amp3=1e-5)
        rows_all[name]['adts_nform'] = [float(a0[c]) for c in ('xx', 'xy', 'yy')]
        rows_all[name]['dadts_ddelta'] = [float(d[c]) for c in ('xx', 'xy', 'yy')]
        rows_all[name]['dadts_norm'] = float(np.sqrt(sum(float(d[c]) ** 2 for c in ('xx', 'xy', 'yy'))))
    P('追跡の写像の数値の正規形 (δ=0 の ADTS と ∂ADTS/∂δ): %.0f s' % (time.time() - t1))
    P('\n解析式の ADTS と、追跡の写像の数値の正規形の ADTS (δ=0)、および ∂(ADTS)/∂δ [1/(m rad)]')
    for name, r in rows_all.items():
        if r['group'] != 'random' or name.endswith('_0'):
            P('  %-14s 解析 %s | 正規形 %s | ∂/∂δ %s' % (name, np.round(r['adts'], 0).tolist(), np.round(r['adts_nform'], 0).tolist(), np.round(r['dadts_ddelta'], 0).tolist()))
    dif = [max(abs(r['adts_nform'][i] - r['adts'][i]) for i in range(3)) / max(r['adts_norm'], 1.0) for r in rows_all.values()]
    P('  解析式と正規形の差 / |ADTS| の最大 %.2e、中央値 %.2e (%d 組)' % (max(dif), float(np.median(dif)), len(dif)))

os.makedirs(C.OUT, exist_ok=True)
json.dump(rows_all, open(os.path.join(C.OUT, 'rdt_table.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

# 表
hdr = '| 名前 | 組の種類 | max|K2L| | ξx,ξy (CALC) | ξ2x,ξ2y (CALC) | Σ|h3| (5 項) | h30000 | h10020 | h21000 | 4 次 非対角 Σ | ADTS xx | xy | yy | 目的/最終 |'
lines = [hdr, '|' + '---|' * 15]
for name, r in rows_all.items():
    h = dict(zip(r['h3_names'], r['h3']))
    lines.append('| %s | %s | %.1f | %.3f, %.3f | %.0f, %.0f | %.1f | %.1f | %.1f | %.1f | %.3g | %.3g | %.3g | %.3g | %.3g |' % (
        name, r['group'], r['k_max'], *r['xi_calc'], *r['xi2_calc'], r['h3_sum'], h['h30000'], h['h10020'], h['h21000'], r['h4_nondiag_sum'] * 1e0,
        *r['adts'], r['obj_rel']))
open(os.path.join(C.OUT, 'rdt_table.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
P('\n' + '\n'.join(lines))

# 目的関数の重みへの感度: J_ref・δ_ref を変えて最適化し直し、最終の組の目的関数が最小値の何倍か・解の距離
P('\n[重み (J_ref、δ_ref) への感度] 全項の最適化 (kmax = 最終の最大、10 初期値)。値 = 最終の組の目的関数 / 最小値、解 1 の K2L の基準解 (J_ref=1e-6、δ_ref=1%) からの距離')
xi0, M = m.xi_matrix()
xi_t = xi0 + M @ kf
kmax = float(np.max(np.abs(kf)))
ref_sol = None
for Jr in (1e-6, 3e-7, 3e-6, 1e-5):
    for dr in (0.01, 0.005, 0.02):
        o = S.SxtObjective(m, J_ref=Jr, delta_ref=dr)
        best = None
        for sd in range(10):
            rng = np.random.default_rng(sd)
            k0 = kf.copy() if sd == 0 else S.project_to_xi(m, np.clip(kf * (1 + 0.3 * rng.standard_normal(m.F)), -kmax, kmax), xi_t)
            k, r = S.optimize(o, k0, xi_t, kmax, maxiter=600)
            if best is None or o.value(k) < best[0]:
                best = (o.value(k), k)
        if ref_sol is None:
            ref_sol = best[1]
        parts = o.parts(kf)
        P('  J_ref=%.0e δ_ref=%.3f: 最終/最小 = %.2f、解の距離 %.2f、最小の解の ADTS = %s、(最終の組の内訳 w3 %.1e / w4 %.1e / adts %.1e / chrom %.1e)' % (
            Jr, dr, o.value(kf) / best[0], np.linalg.norm(best[1] - ref_sol), np.round([float(v) for v in m.adts(best[1]).values()], 0).tolist(),
            parts['w3'], parts['w4'], parts['adts'], parts['chrom']))

# ADTS の配置への依存 (最終の組、3 つの RING の並び)
P('\n[配置への依存] 最終の組を 3 つの並び (RINGSP・RINGINJ・RINGSYM) に置いたときの解析量 (ADTS は並びによらず 1e-9 で同じ、1 次の RDT は並びで変わる)')
for ring in ('RINGSP', 'RINGINJ', 'RINGSYM'):
    _, lt = C.load_lattice(ring, src=src)
    mm = S.SxtModel(lt, C.P0, names=m.names, quad_forms=False)
    h = mm.h1(kf)
    a = mm.adts(kf)
    P('  %-8s 鏡面の軸 %s  |h|(5 項) = %s  ADTS = %s' % (ring, C.mirror_axes(lt)[:2], np.round([abs(h[t]) for t in S.CUBIC5], 2).tolist(),
                                                    np.round([float(a[c]) for c in ('xx', 'xy', 'yy')], 3).tolist()))

# 鏡面対称: リングを鏡面の軸から始めると、1 次の RDT は e^{iπΔν} (Δ = (j−k)νx + (l−m)νy) を掛けると実数になる (複素 1 個 = 実数 1 個の条件)
ax = C.mirror_axes(lat)
if ax:
    import copy
    P('\n[鏡面対称] リングを鏡面の軸 (素子の添字 %d) から始めたときの 1 次の RDT (最終の組)。h e^{-iπΔν} の虚部 / |h| (鏡面対称なら 0)' % ax[0])
    l2 = copy.deepcopy(lat)
    l2.elements = list(lat.elements[ax[0]:]) + list(lat.elements[:ax[0]])
    r2 = l2.twiss_table()[-1]
    q2 = (r2.mux / (2 * math.pi), r2.muy / (2 * math.pi))
    rr = R._resonance_terms_cubic(l2, S.CUBIC5 + S.CHROM3, bend_cubic=True, momentum_eV=C.P0)
    for t, v in rr.total.items():
        j, k_, l_, m_ = [int(c) for c in t[1:5]]
        z = v * np.exp(-1j * math.pi * ((j - k_) * q2[0] + (l_ - m_) * q2[1]))
        P('  %-7s |h| = %8.4f   Im/|h| = %.1e' % (t, abs(v), abs(z.imag) / max(abs(z), 1e-30)))
    rr0 = R._resonance_terms_cubic(lat, S.CUBIC5, bend_cubic=True, momentum_eV=C.P0)
    P('  (リング先頭から始めたときは |h| = %s で、虚部は 0 にならない)' % np.round([abs(v) for v in rr0.total.values()], 2).tolist())

# 図
if not args.no_fig:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    show = [n for n in ('v07_12', 'retune1', 'final', 'opt_full_1') if n in CAND['candidates']]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    for ia, term in enumerate(('h21000', 'h10020', 'h30000')):
        for n in show:
            lat2 = S.with_k2l(lat, m.names, np.array(CAND['candidates'][n]['k']))
            r = R._resonance_terms_cubic(lat2, (term,), bend_cubic=True, momentum_eV=C.P0)
            z = np.concatenate([[0], np.cumsum(r.per_element[:, 0])])
            axs[ia].plot(z.real, z.imag, '-', lw=0.9, label='%s (|h|=%.1f)' % (n, abs(z[-1])))
            axs[ia].plot(z[-1].real, z[-1].imag, 'o', ms=4)
        axs[ia].set_aspect('equal')
        axs[ia].set_title('%s: 要素ごとの寄与の累積 (複素平面、s の順)' % term, fontsize=9)
        axs[ia].legend(fontsize=7)
        axs[ia].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_rdt_walk.png'), dpi=130)
    names = [n for n in rows_all if rows_all[n]['group'] != 'random'] + []
    fig, axs = plt.subplots(1, 4, figsize=(15, 4.2))
    for ia, (key, lab) in enumerate((('w3', '1 次の RDT (W3)'), ('w4', '4 次の非対角 (W4)'), ('adts', 'ADTS'), ('chrom', 'δ の 1 次 (chrom)'))):
        vals = [rows_all[n]['parts'][key] for n in names]
        axs[ia].bar(range(len(names)), vals, color=['C0' if rows_all[n]['group'] == 'historic' else 'C1' for n in names])
        axs[ia].set_yscale('log')
        axs[ia].set_xticks(range(len(names)))
        axs[ia].set_xticklabels(names, rotation=90, fontsize=7)
        axs[ia].set_title('%s の目的関数の項 (青: 歴史的、橙: 最適化)' % lab, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_rdt_sets.png'), dpi=130)
open(os.path.join(C.OUT, 'rdt_table.out.txt'), 'w', encoding='utf-8').write(LOG.getvalue())
