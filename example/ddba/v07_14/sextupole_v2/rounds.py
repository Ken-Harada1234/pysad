# -*- coding: utf-8 -*-
"""試行の各ラウンドで人が選んだ項の指定 (SxtSpec の辞書)。README の経過と同じ。
  py -3.12 rounds.py            各候補の指定から K2L を計算し直して results/cands.json の値と比べる (再現の確認。書き換えない)
  py -3.12 rounds.py r6a        1 つだけ
この関数は項と重みを入力して K2L を返すだけ (pysad.sxtopt.design_k2l)。どの項をどの重みで選ぶかは人が周波数マップ・RDT・ADTS の表を見て決めた (README)。
r5・r7 の走査は round_scan.py、揺らぎの確認 (jit_*) は README の §4。
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

S = V.S
src, lat, m, ks, XT = V.load()
KF = ks['final']
_main = ['w3', 'chrom', 'adts', 'h20200', 'h20020']
BASE = {r['name']: r['value'] for r in S.SxtSpecObjective(m, 'w3=objective, chrom=objective, adts=objective, h20200=objective, h20020=objective', warn=False).table(KF, XT)}
A_FINAL = {c: float(v) for c, v in m.adts(KF).items()}


def near(movers, tol, k0=None):
    """movers の項を (係数 × 完成版の値) へ動かし (目的)、他の主な項 (3 次の RDT 5 項・δ の 3 項・ADTS 3 成分・h20200・h20020) は完成版の値の ±tol の制約にする。"""
    spec = {}
    for n, v in BASE.items():
        if n in movers:
            spec[n] = ('objective', 1.0, movers[n] * v)
        else:
            spec[n] = ('constraint', 1.0, v, tol * abs(v) + (1e-3 if n.startswith('h') and abs(v) < 1 else 0.0))
    return spec


ADTS_KEEP = {('adts_' + c): ('constraint', 1.0, A_FINAL[c], 500.0) for c in ('xx', 'xy', 'yy')}

# (名前: (spec, 設計の呼び出しの引数))
ROUNDS = {
    # ラウンド 1: ADTS だけ
    'r1a': ('adts=objective', dict()),
    # ラウンド 2: 誤差感度 + 3 次の RDT + ADTS (+ 4 次)
    'r2a': ('err3=objective, w3=objective, adts=objective', dict()),
    'r2b': ('err3=objective, w3=objective, adts=objective, w4=objective', dict()),
    # ラウンド 3: 完成版の h20200 だけを 0 へ、3 次の RDT は ±30%、ADTS は ±500
    'r3a': (dict([('h20200', ('objective', 1.0, 0.0))] + [(t, ('constraint', 1.0, complex(m.h1(KF)[t]), 0.3 * abs(m.h1(KF)[t]) + 0.3)) for t in S.CUBIC5]
                 + [('adts_' + c, ('constraint', 1.0, A_FINAL[c], 500.0)) for c in ('xx', 'xy', 'yy')]), dict(dk=6)),
    # ラウンド 4: ADTS は完成版の値 ±500 の制約、3 次の RDT を 0 へ
    'r4a': (dict(ADTS_KEEP, w3=('objective', 1.0)), dict()),
    'r4b': (dict(ADTS_KEEP, w3=('objective', 1.0), h20200=('objective', 1.0)), dict()),
    'r4c': (dict(ADTS_KEEP, w3=('objective', 1.0), w4=('objective', 1.0)), dict()),
    # ラウンド 6: 感度の走査 (round_scan.py の r5) で効いた h00201 と h20200 を組み合わせる
    'r6a': (near({'h00201': 0.0, 'h20200': 0.3}, 0.10), dict(dk=4)),
    'r6b': (near({'h00201': 0.0}, 0.10), dict(dk=4)),
    'r6c': (near({'h00201': 0.0, 'h20200': 0.3}, 0.25), dict(dk=4)),
    # ラウンド 9: r6a の変更 (ξ の目標・ゆるい制約・項の追加)
    'r9a': (near({'h00201': 0.0, 'h20200': 0.3}, 0.10), dict(dk=4, xi_target=np.array([1.0, 1.0]))),
    'r9b': (near({'h00201': 0.0, 'h20200': 0.3}, 0.25), dict(dk=10)),
    'r9c': (near({'h00201': 0.0, 'h20200': 0.3, 'h10020': 0.3, 'h10110': 0.3}, 0.10), dict(dk=6)),
}


def compute(name):
    spec, kw = ROUNDS[name]
    kw = dict(kw)
    xi = kw.pop('xi_target', XT)
    return S.design_k2l(m, spec, k0=KF, xi_target=xi, kmax=V.KMAX, warn=False, **kw)


if __name__ == '__main__':
    cands = json.load(open(os.path.join(V.OUT, 'cands.json'), encoding='utf-8'))
    names = sys.argv[1:] or list(ROUNDS)
    for n in names:
        r = compute(n)
        d = np.max(np.abs(r['k'] - np.array(cands[n]['k'])))
        print('%-4s 最適化 成功=%s 制約違反 %.2g | cands.json との max|ΔK2L| = %.2g' % (n, r['success'], r['constraint_violation'], d))
        print('     ' + V.show_k(m, n, r['k']))
