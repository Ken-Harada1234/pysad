# -*- coding: utf-8 -*-
"""[試作] 解析式 (RDT・ADTS・クロマティシティ) を目的関数にした六極 K2L (20 ファミリー) の最適化 (R7)。使い方 (どこからでも):
  py -3.12 sextupole_optimize.py [--nstart 40] [--ring RING]
変数: 20 ファミリーの K2L [1/m²] (鏡面対称は配置が保つ。ファミリー = セルの種類 × S1〜S5 なので鏡像の六極は同じファミリー)。
制約: 解析式のクロマティシティ (ξx, ξy) = 最終の組の値 (等式)、|K2L| ≤ kmax (既定は最終の組の最大 |K2L| = 32.65、別に 1.25 倍も)。
目的: pysad/sxtopt.py の SxtObjective (物理的な重みだけ。J_ref = 1e-6 m·rad、δ_ref = 1%、各項の重み 1。追跡に合わせた調整はしない)。
実装: モデルは解析式と 1e-12 で一致する多項式 (SxtModel)、勾配は解析、SLSQP。複数の初期値 (最終の組の周りの摂動・一様乱数) から走らせ、離れた局所解を分ける。
出力: results/candidates.json (追跡の判定 sextupole_tracking_check.py が読む候補: 歴史的な組・最適化の結果・最終の周りの乱数の組)、
      results/optimize.out.txt、fig_opt_convergence.png
解析式と CALC の差 (K2² のクロマティシティ) は、CALC のクロマティシティを最終の組の値に合わせる微調整 (2 ファミリーの K2L をニュートン法で動かす) で補正する (候補ごとに記録)。
"""
import argparse
import io
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import resdrive as R
from pysad import sxtopt as S

ap = argparse.ArgumentParser()
ap.add_argument('--nstart', type=int, default=40)
ap.add_argument('--ring', default='RING')
ap.add_argument('--nrnd', type=int, default=8, help='最終の周りの乱数の組の数 (σ=10% と 30% で各 nrnd 個)')
args = ap.parse_args()
os.makedirs(C.OUT, exist_ok=True)
LOG = io.StringIO()


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')


src, lat = C.load_lattice(args.ring)
sets = C.load_sets()
m = S.SxtModel(lat, C.P0)
assert sets['family_names'] == m.names
kf = np.array(sets['sets']['final']['k'])
xi0, M = m.xi_matrix()
xi_an_final = xi0 + M @ kf
calc_final = R.chromaticity_numeric(S.with_k2l(lat, m.names, kf), order=3, dmax=2e-3)
xi_calc_final = np.array(calc_final['d1'])
kmax = float(np.max(np.abs(kf)))
P('リング %s、運動量 %.1f GeV、六極 %d ファミリー、鏡面対称の軸 (素子の添字) %s' % (args.ring, C.P0 / 1e9, m.F, C.mirror_axes(lat)[:2]))
P('最終の組: 解析式 ξ = (%.4f, %.4f)、CALC ξ = (%.4f, %.4f)、max|K2L| = %.3f' % (*xi_an_final, *xi_calc_final, kmax))
obj = S.SxtObjective(m)
f_final = obj.value(kf)
P('目的関数 (最終の組) = %.4e  内訳 %s' % (f_final, {k: '%.3e' % v for k, v in obj.parts(kf).items()}))

# ξ の微調整に使う 2 ファミリー (K2L の相対変化あたりの ξ の効きが独立なもの)
best, pair = 0, (0, 1)
for i in range(m.F):
    for j in range(i + 1, m.F):
        d = abs(np.linalg.det(M[:, [i, j]] * np.abs(kf[[i, j]])))
        if d > best:
            best, pair = d, (i, j)
P('クロマティシティの微調整に使うファミリー:', [m.names[i] for i in pair])


def trim_xi(k, nit=4):
    """CALC の ξ を最終の組の CALC の値に合わせる (2 ファミリーの K2L をニュートン法で動かす。ヤコビ行列は解析式の M)。戻り (k, 残差)"""
    k = k.copy()
    Mp = M[:, list(pair)]
    for _ in range(nit):
        xi = np.array(R.chromaticity_numeric(S.with_k2l(lat, m.names, k), order=3, dmax=2e-3)['d1'])
        r = xi_calc_final - xi
        if np.max(np.abs(r)) < 1e-4:
            break
        k[list(pair)] += np.linalg.solve(Mp, r)
    xi = np.array(R.chromaticity_numeric(S.with_k2l(lat, m.names, k), order=3, dmax=2e-3)['d1'])
    return k, xi - xi_calc_final


def multistart(obj, kmax, n, seed0=0):
    sols = []
    for sd in range(n):
        rng = np.random.default_rng(seed0 + sd)
        if sd == 0:
            k0 = kf.copy()
        elif sd % 2 == 1:
            k0 = S.project_to_xi(m, np.clip(kf * (1 + 0.3 * rng.standard_normal(m.F)), -kmax, kmax), xi_an_final)
        else:
            k0 = S.project_to_xi(m, rng.uniform(-kmax, kmax, m.F) * 0.6, xi_an_final)
        try:
            k, r = S.optimize(obj, k0, xi_an_final, kmax, maxiter=600)
        except Exception:
            continue
        sols.append((obj.value(k), k, r.success, r.nit, k0))
    sols.sort(key=lambda t: t[0])
    return sols


def distinct(sols, tol=0.5, nmax=4):
    out = []
    for v, k, ok, nit, k0 in sols:
        if all(np.linalg.norm(k - o[1]) > tol for o in out):
            out.append((v, k, ok, nit, k0))
        if len(out) >= nmax:
            break
    return out


def history(obj, k0, kmax):
    """収束曲線 (SLSQP の反復ごとの目的関数 / 最終の組の値)"""
    from scipy.optimize import minimize
    xi0_, M_ = obj.m.xi_matrix()
    hist = []
    cons = [{'type': 'eq', 'fun': lambda k: xi0_ + M_ @ k - xi_an_final, 'jac': lambda k: M_}]
    sc = obj.value(k0)
    minimize(lambda k: obj.value(k) / sc, k0, jac=lambda k: obj.grad(k) / sc, method='SLSQP', bounds=[(-kmax, kmax)] * obj.m.F, constraints=cons,
             options=dict(maxiter=600, ftol=1e-14), callback=lambda k: hist.append(obj.value(k) / f_final))
    return [obj.value(k0) / f_final] + hist


cand = {}
labels = {'v07_11': 'v07_11 (初期の 3 個 (S1=S2=0、S3=S5))', 'v07_12': 'v07_12 (初期の 3 個 (S1=S2=0、S3=S5))', 'retune1': '六極の再調整 第 1 版', 'final': '最終 (20231205 再調整版)',
          'rev9': 'v07_16 sx_rev9 (ξ が違う)'}
for key, s in sets['sets'].items():
    cand[key] = dict(k=s['k'], desc=labels.get(key, s['label']), group='historic', trim=None)


def add(name, k, desc, group, trim=True):
    tr = None
    if trim:
        k, res = trim_xi(k)
        tr = [float(x) for x in res]
    cand[name] = dict(k=[float(x) for x in k], desc=desc, group=group, trim=tr)
    return k


# 1. 目的関数の全体の最適化 (kmax = 最終の最大値) と、離れた局所解
t0 = time.time()
sols = multistart(obj, kmax, args.nstart)
P('\n[全項の最適化 kmax=%.2f] %d 回の初期値から %d 個が走った (%.1f s)。目的関数 / 最終の組の値 の分布: %s' % (
    kmax, args.nstart, len(sols), time.time() - t0, np.round([s[0] / f_final for s in sols][:12], 4).tolist()))
dist = distinct(sols)
P('離れた局所解 (K2L の距離 > 0.5):')
for i, (v, k, ok, nit, k0) in enumerate(dist):
    P('  解 %d: 目的 %.4f (最終比)、success=%s、反復 %d、max|K2L| = %.2f、ADTS = %s' % (
        i + 1, v / f_final, ok, nit, np.max(np.abs(k)), np.round([float(x) for x in m.adts(k).values()], 0).tolist()))
    add('opt_full_%d' % (i + 1), k, '全項の最適化 (kmax=%.2f) の局所解 %d' % (kmax, i + 1), 'optimized')

# 2. 項ごとの最適化 (どの項が追跡と相関するかを分けるための設計された組)
for tag, terms, desc in (('opt_rdt1', ('w3', 'chrom'), '1 次の RDT (3 次とδの項) だけを最小化'), ('opt_adts', ('adts',), 'ADTS だけを最小化'),
                         ('opt_rdt4', ('w4',), '4 次の非対角の項 (正規形) だけを最小化')):
    o2 = S.SxtObjective(m, terms=terms)
    s2 = multistart(o2, kmax, max(8, args.nstart // 4), seed0=100)
    v, k, ok, nit, k0 = s2[0]
    P('[%s] %s: 目的 (その項) %.3e (最終の組 %.3e)、全項の目的 %.4f (最終比)' % (tag, desc, v, o2.value(kf), obj.value(k) / f_final))
    add(tag, k, desc + ' (kmax=%.2f)' % kmax, 'optimized')

# 3. 上限を 1.25 倍にした全項の最適化
s3 = multistart(obj, 1.25 * kmax, max(8, args.nstart // 4), seed0=200)
v, k, ok, nit, k0 = s3[0]
P('[opt_full_k125] kmax=%.2f: 目的 %.4f (最終比)、max|K2L| = %.2f' % (1.25 * kmax, v / f_final, np.max(np.abs(k))))
add('opt_full_k125', k, '全項の最適化 (kmax=%.2f = 最終の最大の 1.25 倍)' % (1.25 * kmax), 'optimized')

# 4. 最終の周りの乱数の組 (ξ を解析式で最終の値に射影)。K2L の相対ずれ σ = 10%、30%
for sig, tag in ((0.10, 'rnd10'), (0.30, 'rnd30')):
    for i in range(args.nrnd):
        rng = np.random.default_rng(1000 + int(sig * 100) * 10 + i)
        k = S.project_to_xi(m, kf * (1 + sig * rng.standard_normal(m.F)), xi_an_final)
        add('%s_%d' % (tag, i), k, '最終の組の周りの乱数 (K2L の相対ずれ σ=%d%%)' % (sig * 100), 'random')

json.dump(dict(family_names=m.names, kmax=kmax, xi_calc_final=list(map(float, xi_calc_final)), xi_an_final=list(map(float, xi_an_final)), f_final=f_final,
               ring=args.ring, candidates=cand), open(os.path.join(C.OUT, 'candidates.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
P('\n候補 %d 個を results/candidates.json に書いた (CALC ξ の微調整後の残差の最大 %.1e)' % (
    len(cand), max([max(abs(np.array(c['trim']))) for c in cand.values() if c['trim']] + [0])))
P('\n名前            目的/最終  w3        chrom     w4        adts     ADTS(xx,xy,yy)                max|K2L|  ξ_an')
for n, c in cand.items():
    k = np.array(c['k'])
    p = obj.parts(k)
    P('%-14s %8.4f  %.2e  %.2e  %.2e  %.2e  %-28s  %6.2f  (%.3f, %.3f)' % (
        n, obj.value(k) / f_final, p['w3'], p['chrom'], p['w4'], p['adts'], np.round([float(x) for x in m.adts(k).values()], 0).tolist(),
        np.max(np.abs(k)), *m.xi(k)))

# 図
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, (v, k, ok, nit, k0) in enumerate(dist[:3]):
        h = history(obj, k0, kmax)
        ax[0].semilogy(h, label='解 %d の初期値 (反復 %d)' % (i + 1, len(h) - 1))
    k0 = S.project_to_xi(m, np.clip(kf * (1 + 0.3 * np.random.default_rng(5).standard_normal(m.F)), -kmax, kmax), xi_an_final)
    ax[0].semilogy(history(obj, k0, kmax), label='最終の周りの摂動 30%')
    ax[0].axhline(1.0, color='k', lw=0.8, ls='--')
    ax[0].text(0.98, 1.1, '最終の組 = 1', ha='right', transform=ax[0].get_yaxis_transform())
    ax[0].set_xlabel('SLSQP の反復')
    ax[0].set_ylabel('目的関数 / 最終の組の値')
    ax[0].legend(fontsize=8)
    ax[0].set_title('最適化の収束 (全項、kmax=最終の最大)')
    x = np.arange(m.F)
    ax[1].bar(x - 0.27, kf, 0.27, label='最終の組')
    ax[1].bar(x, np.array(cand['opt_full_1']['k']), 0.27, label='全項の最適化 解 1')
    ax[1].bar(x + 0.27, np.array(cand['opt_rdt1']['k']), 0.27, label='1 次の RDT だけ')
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(m.names, rotation=90, fontsize=7)
    ax[1].set_ylabel('K2L [1/m²]')
    ax[1].legend(fontsize=8)
    ax[1].set_title('ファミリーごとの K2L')
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_opt_convergence.png'), dpi=130)
except Exception as ex:
    P('図の作成に失敗:', ex)
open(os.path.join(C.OUT, 'optimize.out.txt'), 'w', encoding='utf-8').write(LOG.getvalue())
