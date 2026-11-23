# -*- coding: utf-8 -*-
"""入門 §4: 共鳴駆動項 h の「6 極 1 本ずつの寄与の和」と、厚い 6 極 → 薄いレンズの極限。FODO リング (3 GeV、105.6 m、SF・SD の K2L = 4, −6、各 48 本)。
  左: h30000 (3νx の共鳴) と h21000 (νx の共鳴) の 6 極ごとの寄与 (複素数) を s の順に足したときの軌跡 (複素平面)。6 極の位置の位相 e^{3iμx} が回るので、
      個々の寄与の大きさの和 Σ|h_i| に比べて、和 |Σ h_i| はずっと小さい (位相が回って打ち消す)。
  右: 厚い 6 極 (L = 0.2 m) の解析式 (厳密) と、6 極を n 個の薄レンズに分けた列 (n = 1 … 64) の h21000・h30000 の差。n を 2 倍にすると差は 1/4。
  (薄い 6 極の式との照合は ADTS_summary_rdt_check.py)
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/primer/primer_rdt.py   出力: primer_rdt.png、標準出力
"""
import math
import numpy as np
import primer_common as C
plt = C.setup_plot()
from pysad import resdrive as R

lat, rows, tw, nu, circ = C.fodo()
r = R.resonance_terms(lat, terms=['h21000', 'h30000', 'h10110', 'h10020', 'h10200'], method='analytic', bend_cubic=False)      # 6 極だけ (既定は偏向磁石の項つき)
print('FODO リング (3 GeV、105.6 m)、νx = %.6f、νy = %.6f、位相の基準 = リング先頭、厚い解析式' % nu)
print('3νx = %.6f、νx = %.6f、νx − 2νy = %.6f、νx + 2νy = %.6f' % (3 * nu[0], nu[0], nu[0] - 2 * nu[1], nu[0] + 2 * nu[1]))
idx = [i for i, nm in enumerate(r.elem_names) if nm in ('SF', 'SD')]
print('6 極の数 = %d' % len(idx))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for a, term in zip(axs[:2], ('h30000', 'h21000')):
    col = r.names.index(term)
    per = r.per_element[idx, col]
    order = np.argsort(np.array(r.elem_s)[idx])
    path = np.r_[0, np.cumsum(per[order])]
    print('%s: Σ|h_i| = %.4f、|Σ h_i| = %.4f、Σ h_i = %.6f%+.6fi (total = %.6f%+.6fi)' % (
        term, np.abs(per).sum(), abs(per.sum()), per.sum().real, per.sum().imag, r.total[term].real, r.total[term].imag))
    a.plot(path.real, path.imag, '-o', ms=2.5, lw=0.8)
    a.plot([0, path[-1].real], [0, path[-1].imag], 'r-', lw=1.5, label='和 |Σh| = %.3f' % abs(per.sum()))
    a.set_aspect('equal', adjustable='datalim')
    a.set_title('%s: 6 極 %d 本の寄与を s の順に足す (Σ|h_i| = %.2f)' % (term, len(idx), np.abs(per).sum()), fontsize=9)
    a.set_xlabel('Re h'); a.set_ylabel('Im h'); a.legend(fontsize=8); a.grid(alpha=0.3)
# 厚い → 薄い
ref = R.resonance_terms(lat, terms=['h21000', 'h30000'], method='analytic', bend_cubic=False)
a = axs[2]
ns = [1, 2, 4, 8, 16, 32, 64]
err = {t: [] for t in ('h21000', 'h30000')}
for n in ns:
    th = R.slice_thin(lat, n)
    rt = R.resonance_terms(th, terms=['h21000', 'h30000'], method='analytic', bend_cubic=False)
    for t in err:
        err[t].append(abs(rt.total[t] - ref.total[t]) / abs(ref.total[t]))
print('厚い 6 極 (L = 0.2 m) と n 分割の薄レンズ列の |Δh|/|h|:')
print('n      ' + '  '.join('%10d' % n for n in ns))
for t in err:
    print('%s ' % t + '  '.join('%10.3e' % e for e in err[t]))
    print('   比 (n → 2n): ' + '  '.join('%6.2f' % (err[t][i] / err[t][i + 1]) for i in range(len(ns) - 1)))
    a.loglog(ns, err[t], 'o-', label=t)
a.loglog(ns, err['h21000'][0] / np.array(ns) ** 2, 'k:', lw=0.8, label=r'$1/n^2$')
a.set_xlabel('分割数 n (薄レンズの数)'); a.set_ylabel('|h(薄レンズ n 個) − h(厚い)| / |h(厚い)|')
a.set_title('厚い 6 極 → 薄レンズ列: 差は 1/n²', fontsize=9)
a.legend(fontsize=8); a.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(C.os.path.join(C.HERE, 'primer_rdt.png'), dpi=110)
