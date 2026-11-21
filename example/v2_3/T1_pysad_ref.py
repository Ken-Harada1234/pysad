# -*- coding: utf-8 -*-
"""v2.3.0 T1: pysad 側の基準値 (光学、IBS 込みのビーム、Touschek 寿命) を作る。elegant と同じ条件を揃えるための入力にもなる。

使い方 (リポジトリ直下から。PYTHONUTF8=1 を付ける):
  PYTHONUTF8=1 py -3.12 example/v2_3/T1_pysad_ref.py astrid2|kekls
出力: example/v2_3/T1_pysad_<名前>.ref.json (elegant の touschekLifetime に渡す値と、pysad の寿命)、標準出力 (T1_pysad_<名前>.ref.out.txt に保存した)。
内容: 光学 (チューン、色収差 (∂ν/∂δ の数値微分)、β・η の最大・最小、αc、U0、自然エミッタンス・エネルギー拡がり・バンチ長 (RFSW・RAD ON の Emittance[]))、
  INTRA ON・PBUNCH = 1e9 粒子/バンチの Emittance[] の収束値 (EMITXC・EMITYC・EMITZC、SIGE、SIGZ)、一様な受容 δ での Touschek 寿命 (TouschekTable の損失率の和)。
局所の LMA の寿命は T1_lma_compare.py が作る。v2 モード。
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

PB = 1e9
UNIFORM = {'astrid2': [0.004, 0.005, 0.006, 0.010, 0.016, 0.0161756586929], 'kekls': [0.004, 0.005, 0.006, 0.010, 0.020, 0.0299397749]}


def main(name):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.touschek import touschek_lifetime
    import T1_sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.010)
    g = ffs.I.globals
    cond = st.conditions()
    out = dict(name=name, momentum_GeV=cond['momentum_GeV'], circumference=cond['circumference_m'], U0_eV=cond['U0_eV'], Vrf=cond['Vrf_eV'],
               nu_s=cond['nu_s'], alpha_c=cond['momentum_compaction'], bucket_height=cond['bucket_height'], phi_s=cond['phi_s_rad'])
    from pysad.trksetup import find_periodicity
    per = find_periodicity(st.els, st.rows)
    out['cell_length'] = float(per['cell_length'])
    out['ncell'] = int(per['nper'])
    out['natural'] = {k: float(g[k]) for k in ('EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ')}
    tabs = {d: ffs.lat.twiss_table(delta=d) for d in (0.0, 1e-3, -1e-3)}
    nx = lambda d: tabs[d][-1].mux / (2 * math.pi)
    ny = lambda d: tabs[d][-1].muy / (2 * math.pi)
    t = tabs[0.0]
    out['tune'] = [nx(0.0), ny(0.0)]
    out['chrom'] = [(nx(1e-3) - nx(-1e-3)) / 2e-3, (ny(1e-3) - ny(-1e-3)) / 2e-3]
    out['beta_x'] = [min(r.betx for r in t), max(r.betx for r in t)]
    out['beta_y'] = [min(r.bety for r in t), max(r.bety for r in t)]
    out['eta_x_max'] = max(r.etx for r in t)
    ffs.flags['INTRA'] = True
    g['PBUNCH'] = PB
    ffs.flags['RFSW'] = True
    ffs.flags['RAD'] = True
    ffs.run_emit(quiet=True)
    out['ibs'] = {k: float(g[k]) for k in ('EMITXC', 'EMITYC', 'EMITZC', 'SIGE', 'SIGZ')}
    tab = ffs.touschek['table']
    sige = float(g['SIGE'])
    out['uniform_a'] = UNIFORM[name]
    out['pysad_uniform'] = {}
    for a in UNIFORM[name]:
        out['pysad_uniform']['%.10g' % a] = float(touschek_lifetime(tab, None, None, a / sige, sige))
    # elegant に渡す値 (トップレベルにも置く)
    for k in ('EMITXC', 'EMITYC', 'EMITZC', 'SIGE', 'SIGZ'):
        out[k] = out['ibs'][k]
    out['PBUNCH'] = PB
    out['pysad_table'] = 'Q1_lma_%s.table.txt' % name if name == 'astrid2' else 'T1_pysad_%s.lma.txt' % name
    # pysad の局所の LMA (表があれば) の寿命: lma.touschek_lifetime_local (1 セルを周期的に延ばす)
    tp = os.path.join(HERE, out['pysad_table'])
    if os.path.exists(tp):
        from types import SimpleNamespace
        from pysad.lma import touschek_lifetime_local
        rows = [ln.split() for ln in open(tp, encoding='utf-8') if ln.strip() and not ln.startswith('#')]
        cs = [r for r in rows if float(r[1]) < per['cell_length'] + 1e-9]
        lm = SimpleNamespace(s=np.array([float(r[1]) for r in cs]), dp_plus=np.array([float(r[3]) for r in cs]),
                             dp_minus=np.array([float(r[4]) for r in cs]), setup=st)
        res = touschek_lifetime_local(tab, lm, setup=st, period=per['cell_length'], s0=per['s_starts'][0])
        out['pysad_local'] = float(res['lifetime'])
    json.dump(out, open(os.path.join(HERE, 'T1_pysad_%s.ref.json' % name), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
