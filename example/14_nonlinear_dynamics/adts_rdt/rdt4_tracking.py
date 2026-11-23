# -*- coding: utf-8 -*-
"""4 次の共鳴駆動項 h_jklm (6 極の 2 次を含む) の解析式 (rdt2nd.rdt4、resonance_terms の 4 次の名前) を、追跡の 1 周の写像の 3 次のテイラー展開から
取り出した 4 次の生成関数 f4 (Dragt-Finn。rdttrack.ring_map3・rdt4_from_map) と比べる。

使い方 (リポジトリ直下から)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_tracking.py table fodo|fodo_rect|fodo_asym|astrid|astrid_mpw|kekls   λ = 0, 0.5, 1, 2 の走査 (解析対追跡の表、内訳、λ の 2 次の依存)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_tracking.py adts  fodo|astrid|astrid_mpw|kekls                      h22000 などの対角部分 + 正規形の補正 = adts() と、追跡 + NAFF (direct_adts) との比較
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_tracking.py plot                                                   rdt4_table.png (table の出力から)
条件 (すべて): pysad k64 の追跡 (RF 空洞止め、放射なし、結合なし、閉軌道のまわり、δ=0)、リング先頭が位相の基準 (μ=0)。3 次のテイラー展開: 正規化座標、方向 48 本、
  奇関数の差分 (a, 2a, 3a)、a = amp3 (FODO・ASTRID2 1e-4 √m、大きなリング A 1e-5 √m。大きなリングは 1e-4 だと 7 次の項で 1e-3 外れる)、2 次のテイラー係数の振幅 1e-5 √m。
  解析: bend_nslice=16、厚い 6 極は閉形式。偏向磁石の 4 次の項は要素ごとの写像から (docs/25 §10.3)。6 極の強さは K2 を λ 倍 (8 極は変えない)。
"""
import json
import math
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, rdttrack as T, rdt2nd as Q, compat
compat.set_mode('k64')
import adts_bend as M

AMP3 = {'fodo': 1e-4, 'fodo_rect': 1e-4, 'fodo_asym': 1e-4, 'astrid': 1e-4, 'astrid_mpw': 1e-4, 'kekls': 1e-5}
SHOW = ['h22000', 'h11110', 'h00220', 'h31000', 'h40000', 'h20020', 'h20200', 'h11200', 'h11020', 'h20110', 'h02110', 'h00310', 'h00400', 'h13000', 'h04000']
OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(s)


def save(name, ring, data=None):
    base = os.path.join(HERE, 'rdt4_%s_%s.pysad' % (name, ring))
    with open(base + '.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')
    if data is not None:
        with open(base + '.json', 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else str(o))


def cx(v):
    return '%+.5e%+.5ej' % (v.real, v.imag)


def one(key, lam, amp3):
    S = M.setup(key, lam)
    rm = T.ring_map3(S['lat'].elements, S['p'], S['rows'][0], amp3=amp3)
    trk, info = T.rdt4_from_map(rm)
    info['fit'] = rm.fit_res3
    t0 = time.time()
    r = Q.rdt4(S['lat'], rows=S['rows'])
    return S, trk, info, r, time.time() - t0


def cmd_table(key):
    amp3 = AMP3[key]
    lams = (0.0, 0.5, 1.0, 2.0)
    data = {'ring': key, 'amp3': amp3, 'lam': {}}
    res = {}
    for lam in lams:
        S, trk, info, r, ta = one(key, lam, amp3)
        res[lam] = (trk, r)
        if lam == lams[0]:
            P('R5 4 次の共鳴駆動項: 解析 (rdt4: 8 極 1 次 + 運動学 + 端のフリンジ + 偏向磁石の 4 次 + 3 次の 2 次) 対 追跡の 3 次のテイラー展開 (Dragt-Finn の f4): %s (%.2f GeV、周長 %.1f m、νx=%.6f νy=%.6f)' % (key, S['p'] / 1e9, S['circ'], S['nu'][0], S['nu'][1]))
            P('条件: pysad k64 の追跡、RF・放射なし、閉軌道のまわり、δ=0、リング先頭が位相の基準。3 次のテイラー展開 amp3=%.0e √m (方向 48 本、a・2a・3a の奇関数の差分)、bend_nslice=16' % amp3)
            for w in r.warnings[:4]:
                P('  (解析の注意) ' + w)
        mx = max(abs(trk[n]) for n in r.names)
        d = [abs(r.total[n] - trk[n]) / abs(trk[n]) for n in r.names if abs(trk[n]) > 1e-3 * mx]
        P('')
        P('--- 6 極 ×%g: 解析 %.1f s、追跡の f4 のシンプレクティック残差 %.1e、フィット残差 %.1e。|h| > 1e-3 max の項 %d 個 (35 項中) の最大 相対差 %.2e (中央値 %.2e) ---' % (lam, ta, info['symp'], info['fit'], len(d), max(d), float(np.median(d))))
        if lam in (1.0,):
            P('項       解析                         追跡                         相対差')
            for n in SHOW:
                P('%-8s %s  %s  %.2e' % (n, cx(r.total[n]), cx(trk[n]), abs(r.total[n] - trk[n]) / abs(trk[n])))
            P('内訳 (対角 h22000 / h11110 / h00220、解析):')
            for k, v in r.parts.items():
                P('  %-9s %+.4f %+.4f %+.4f' % (k, v[Q._IDX4[(2, 2, 0, 0)]].real, v[Q._IDX4[(1, 1, 1, 1)]].real, v[Q._IDX4[(0, 0, 2, 2)]].real))
        data['lam'][str(lam)] = {'max_rel': max(d), 'median_rel': float(np.median(d)), 'symp': info['symp'], 'terms': {n: [r.total[n], trk[n]] for n in SHOW}}
    P('')
    P('λ 依存: h(λ) = h0 + aλ + bλ² を λ = 0, 1, 2 の値から決め、λ = 0.5 を予測 (6 極の 1 次 × 他の 3 次の項は λ の 1 次、6 極どうしは 2 次)。|予測 − 実測|/|h(0.5)|:')
    mx = max(abs(res[1.0][0][m]) for m in res[1.0][1].names)
    for lab, ix in (('追跡', 0), ('解析', 1)):
        dd = []
        for n in res[1.0][1].names:
            v = [res[l][0][n] if ix == 0 else res[l][1].total[n] for l in (0.0, 1.0, 2.0, 0.5)]
            if abs(v[3]) < 1e-3 * mx:
                continue
            b = (v[2] - 2 * v[1] + v[0]) / 2.0
            a = v[1] - v[0] - b
            dd.append(abs(v[0] + a * 0.5 + b * 0.25 - v[3]) / abs(v[3]))
        P('  %s: 最大 %.2e (中央値 %.2e)' % (lab, max(dd), float(np.median(dd))))
    P('')
    P('各 λ の h22000 / h11110 / h00220 (解析 | 追跡)')
    for lam in lams:
        trk, r = res[lam]
        P('  λ=%-4g %s' % (lam, '  '.join('%+.4f|%+.4f' % (r.total[n].real, trk[n].real) for n in ('h22000', 'h11110', 'h00220'))))
    save('table', key, data)


def cmd_adts(key):
    from pysad import adtsdirect as D
    S = M.setup(key, 1.0)
    P('R5 ADTS との整合: %s (%.2f GeV、周長 %.1f m)。対角の h22000・h11110・h00220 (写像) から ∂νx/∂Jx = 4 h22000/π、∂νx/∂Jy = 2 h11110/π、∂νy/∂Jy = 4 h00220/π、' % (key, S['p'] / 1e9, S['circ']))
    P('正規形の補正 (3 次の項の非共鳴の摂動、nf_correction) を足すと adts() の 6 極 2 次 + 偏向磁石 3 次 + 交差項と一致するはず。[1/(m rad)]')
    A = R.adts(S['lat'], rows=S['rows'])
    r0 = Q.rdt4(S['lat'], rows=S['rows'], kinematic=False, fringe=False, octupole=False, bend4=False)
    m0, c0 = r0.adts(), Q.nf_correction(r0.h3, r0.qx, r0.qy)
    ex = {k: A.sext[k] + A.bend[k] + A.cross[k] for k in ('xx', 'xy', 'yy')}
    P('(1) 6 極 2 次 + 偏向磁石の 3 次 + 交差項 (4 次の 1 次の項なし):')
    P('    h からの写像の部分  %s' % ' '.join('%+.5f' % m0[k] for k in m0))
    P('    正規形の補正        %s' % ' '.join('%+.5f' % c0[k] for k in m0))
    P('    合計                %s' % ' '.join('%+.5f' % (m0[k] + c0[k]) for k in m0))
    P('    adts() (R3)         %s   相対差 %s' % (' '.join('%+.5f' % ex[k] for k in m0), ' '.join('%.1e' % abs((m0[k] + c0[k] - ex[k]) / ex[k]) for k in m0)))
    a4, r4 = Q.adts_from_rdt4(S['lat'], rows=S['rows'])
    P('(2) 全部入り (8 極 + 運動学 + フリンジ + 偏向磁石の 4 次の項まで h に入れる)。直接の追跡 (追跡 + NAFF + フィット、direct_adts) と比べる:')
    jb = S['jb']
    tr = D.direct_adts(S['lat'], S['p'], jb, order=2, nturn=1024)
    t = {'xx': tr['xx'], 'xy': tr['xy_mean'], 'yy': tr['yy']}
    P('    条件: Jmax=%.1e、70 粒子 (±位相の対)、1024 周、NAFF Hann p=1、order=2' % jb)
    P('    追跡 + NAFF        %s' % ' '.join('%+.4f' % t[k] for k in t))
    P('    adts().total (R3)  %s   相対差 %s' % (' '.join('%+.4f' % A.total[k] for k in t), ' '.join('%+.2e' % ((A.total[k] - t[k]) / t[k]) for k in t)))
    P('    rdt4 + 補正 (R5)   %s   相対差 %s' % (' '.join('%+.4f' % a4[k] for k in t), ' '.join('%+.2e' % ((a4[k] - t[k]) / t[k]) for k in t)))
    save('adts', key, {'ring': key, 'track': t, 'adts_R3': dict(A.total), 'rdt4_nf': a4})


def cmd_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    rings = ('fodo', 'astrid', 'astrid_mpw', 'kekls')
    w = 0.2
    for i, key in enumerate(rings):
        fn = os.path.join(HERE, 'rdt4_table_%s.pysad.json' % key)
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn, encoding='utf-8'))
        lams = sorted(d['lam'], key=float)
        ax[0].bar(np.arange(len(lams)) + (i - 1.5) * w, [d['lam'][l]['max_rel'] for l in lams], w, label=key)
    ax[0].set_xticks(range(4))
    ax[0].set_xticklabels(['x0', 'x0.5', 'x1', 'x2'])
    ax[0].set_yscale('log')
    ax[0].set_ylabel('max |analytic - tracking| / |tracking|')
    ax[0].set_title('4th-order RDT, 35 terms')
    ax[0].legend()
    for key in rings:
        fn = os.path.join(HERE, 'rdt4_table_%s.pysad.json' % key)
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn, encoding='utf-8'))['lam']['1.0']['terms']
        ax[1].plot([abs(complex(*v[1])) for v in d.values()], [abs(complex(*v[0]) - complex(*v[1])) / abs(complex(*v[1])) for v in d.values()], 'o', label=key)
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('|h| (tracking)')
    ax[1].set_ylabel('relative difference')
    ax[1].set_title('x1, 15 terms per ring')
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'rdt4_table.png'), dpi=110)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'table':
        cmd_table(sys.argv[2])
    elif cmd == 'adts':
        cmd_adts(sys.argv[2])
    elif cmd == 'plot':
        cmd_plot()
