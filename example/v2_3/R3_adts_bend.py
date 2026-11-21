# -*- coding: utf-8 -*-
"""v2.3.0 R3: 偏向磁石の 3 次の項 (曲率 (h/2)x(px²+py²)・硬い端 ±(h/2)px y²・複合機能の本体の K1 h x³/6・E1/E2≠0 の端) を入れた解析の ADTS
(resdrive.adts(bend_cubic=True)) を、実際の強さの追跡 + NAFF + 多項式フィット (adtsdirect.direct_adts、pysad k64。gs `NORFSW;` で追加確認) と比べる。

使い方 (リポジトリ直下から)
  py -3.12 example/v2_3/R3_adts_bend.py scan  fodo|astrid|astrid_mpw|kekls [--lams 0,0.5,1,2] [--jmax J] [--turns 1024] [--nslice 16]
  py -3.12 example/v2_3/R3_adts_bend.py conv  fodo|astrid|astrid_mpw|kekls            (偏向磁石の本体の分割数 nslice への収束。解析のみ)
  py -3.12 example/v2_3/R3_adts_bend.py gs    fodo|astrid|kekls --lam 1               (WSL の gs。NORFSW;。timeout 1800 s。1 本ずつ)
  py -3.12 example/v2_3/R3_adts_bend.py plot                                          (R3_adts_bend.png: scan の出力から)
条件: 3 GeV 等はラティスのまま (FODO は SF・SD の K2L = 4, −6)、RF 空洞は止める (pysad は rfsw=False、gs は NORFSW)、放射なし、結合なし、δ=0 (閉軌道は 0)。
6 極の強さは K2 (数値) の文字列を λ 倍する (8 極 K3 は変えない)。J の格子: Jx, Jy ∈ {0, 1/5, ..., 1}·Jmax の 35 点 × (φ, φ+π) の対 = 70 粒子、位相 φx = φy = 0。
NAFF: Laskar の Hann 窓 (p=1)、探索は線形のチューン ± 0.02、1024 周。フィット: νx (Jx>0)・νy (Jy>0) を全次数 ≤ 2 の多項式で。"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, adtsdirect as D, compat
from pysad.parser import load_sad_string
compat.set_mode(os.environ.get('PYSAD_COMPAT', 'k64'))
import Q2_adts_gs as G

RINGS = {
    'fodo': (os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), 'RING', 3e9, {'SF': 4.0, 'SD': -6.0}, 1e-7),
    'fodo_rect': (os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), 'RING', 3e9, {'SF': 4.0, 'SD': -6.0}, 1e-7),      # 長方形磁石 E1=E2=ANGLE/2
    'fodo_asym': (os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), 'RING', 3e9, {'SF': 4.0, 'SD': -6.0}, 1e-7),      # E1=0.8·ANGLE、E2=−0.3·ANGLE
    'astrid': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING', 0.58e9, None, 2e-8),
    'astrid_mpw': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), 'RINGMPW', 0.58e9, None, 2e-8),
    'kekls': (os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), 'RING', 3e9, None, 2e-10),
}
NSTEP = 5
KEYS = ('xx', 'xy', 'yy')


def setup(key, lam=1.0, mu=1.0):
    path, line, p, k2over, jb = RINGS[key]
    src = open(path, encoding='utf-8').read()
    if k2over:
        src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k2over['SF'], k2over['SD']))
    E_OLD = 'E1=0.1308996938995747 E2=0.1308996938995747'
    if key == 'fodo_rect':
        src = src.replace(E_OLD, 'E1=0.5 E2=0.5')       # SAD の E1・E2 は ANGLE に対する比 (元のファイルは 0.1309 倍 = 0.0343 rad)
    elif key == 'fodo_asym':
        src = src.replace(E_OLD, 'E1=0.8 E2=-0.3')
    if lam != 1.0 or mu != 1.0:
        src = G.scale_text(src, lam, mu)
    lat = load_sad_string(src).build_lattice(line)
    rows = lat.twiss_table()
    r0 = rows[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    nu = (float(rows[-1].mux) / (2 * math.pi), float(rows[-1].muy) / (2 * math.pi))
    return dict(key=key, src=src, line=line, p=p, lat=lat, rows=rows, tw=tw, nu=nu, jb=jb, circ=float(rows[-1].s), lam=lam)


def f3(d):
    return 'xx=% .5f  xy=% .5f  yy=% .5f' % (d['xx'], d['xy'], d['yy'])


def rel(a, b):
    return {k: (a[k] - b[k]) / abs(b[k]) for k in KEYS}


def fr3(d):
    return '%+.2e %+.2e %+.2e' % (d['xx'], d['xy'], d['yy'])


def analytic(S, nslice=16):
    A0 = R.adts(S['lat'], rows=S['rows'], bend_cubic=False)
    A1 = R.adts(S['lat'], rows=S['rows'], bend_cubic=True, bend_nslice=nslice)
    return A0, A1


def cmd_scan(key, lams, jmax, nturn, nslice):
    out = []
    res = {'ring': key, 'jmax': jmax, 'nturn': nturn, 'nslice': nslice, 'rows': []}
    for lam in lams:
        S = setup(key, lam)
        t0 = time.time()
        A0, A1 = analytic(S, nslice)
        ta = time.time() - t0
        tr = D.direct_adts(S['lat'], S['p'], jmax, order=2, nturn=nturn)
        t = {'xx': tr['xx'], 'xy': tr['xy_mean'], 'yy': tr['yy']}
        if lam == lams[0]:
            out.append('リング %s: p=%.4g GeV/c、周長 %.3f m、νx=%.6f νy=%.6f、k64 の追跡 (rfsw=False、放射なし、δ=0、閉軌道 0)、J の格子 Jmax=%.1e (6×6、±位相の対 = 70 粒子)、%d 周、NAFF Hann p=1、order=2、偏向磁石の分割 nslice=%d' % (
                key, S['p'] / 1e9, S['circ'], S['nu'][0], S['nu'][1], jmax, nturn, nslice))
            for w in A1.warnings:
                out.append('(解析の注意) ' + w)
        out.append('')
        out.append('--- 6 極の強さ ×%g (解析 %.1f s) ---' % (lam, ta))
        out.append('  内訳 [1/(m rad)]')
        for nm, d in (('6 極 2 次', A1.sext), ('8 極 1 次', A1.octu), ('運動学的 p^4/8', A1.kin), ('K1 の端のフリンジ', A1.fringe),
                      ('偏向磁石の 3 次どうし (R3)', A1.bend), ('6 極 × 偏向磁石の交差項 (R3)', A1.cross)):
            out.append('    %-28s %s' % (nm, f3(d)))
        t0_ = dict(A0.total)
        t1_ = dict(A1.total)
        out.append('  解析の合計 (R2、偏向磁石の 3 次なし)   %s' % f3(t0_))
        out.append('  解析の合計 (R3、偏向磁石の 3 次あり)   %s' % f3(t1_))
        out.append('  追跡 (pysad k64、直接のフィット)       %s    (フィットの rms %.1e %.1e、∂νx/∂Jy と ∂νy/∂Jx の差 %.1e)' % (
            f3(t), tr['rms'][0], tr['rms'][1], abs(tr['xy'] - tr['yx']) / max(abs(t['xy']), 1e-30)))
        out.append('  相対差 (解析 − 追跡)/|追跡|  R2 (なし)  %s' % fr3(rel(t0_, t)))
        out.append('                              R3 (あり)  %s' % fr3(rel(t1_, t)))
        res['rows'].append({'lam': lam, 'track': t, 'R2': t0_, 'R3': t1_, 'parts': {'sext': dict(A1.sext), 'octu': dict(A1.octu), 'kin': dict(A1.kin), 'fringe': dict(A1.fringe),
                            'bend': dict(A1.bend), 'cross': dict(A1.cross)}, 'rel_R2': rel(t0_, t), 'rel_R3': rel(t1_, t), 'rms': list(tr['rms'])})
        print('\n'.join(out[-14:]), flush=True)
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'R3_scan_%s.pysad.txt' % key), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(res, open(os.path.join(HERE, 'R3_scan_%s.pysad.json' % key), 'w'), indent=1)


def cmd_conv(key, nslices=(1, 2, 4, 8, 16, 32, 64)):
    S = setup(key)
    out = ['リング %s: 偏向磁石の本体の分割数 nslice への収束 (解析のみ、6 極の強さ ×1)' % key,
           '  nslice   偏向磁石の 3 次どうし (xx xy yy)                 6 極 × 偏向磁石 (xx xy yy)                    合計 (xx xy yy)   [最大の nslice との相対差]']
    ref = None
    rows = []
    for n in nslices:
        t0 = time.time()
        A = R.adts(S['lat'], rows=S['rows'], bend_cubic=True, bend_nslice=n)
        rows.append((n, dict(A.bend), dict(A.cross), dict(A.total), time.time() - t0))
    ref = rows[-1][3]
    for n, b, c, tot, tm in rows:
        out.append('  %5d    %s   %s   %s   [%s]  (%.1f s)' % (n, '%.4f %.4f %.4f' % (b['xx'], b['xy'], b['yy']), '%.4f %.4f %.4f' % (c['xx'], c['xy'], c['yy']),
                                                               '%.3f %.3f %.3f' % (tot['xx'], tot['xy'], tot['yy']), fr3(rel(tot, ref)), tm))
    txt = '\n'.join(out)
    print(txt)
    open(os.path.join(HERE, 'R3_conv_%s.pysad.txt' % key), 'w', encoding='utf-8').write(txt + '\n')


def cmd_gs(key, lam, jmax, nturn, nslice):
    S = setup(key, lam)
    A0, A1 = analytic(S, nslice)
    jx, jy = D.grid_actions(jmax, NSTEP)
    JX, JY = np.r_[jx, jx], np.r_[jy, jy]
    PX = np.r_[np.zeros(len(jx)), np.full(len(jx), math.pi)]
    c = D.initial_coords(S['tw'], JX, JY, PX, PX)
    # 位相は φx = φy = 0 (対は π)。script は SAD の源 (K2 は S['src'] で λ 倍済み) をそのまま使う
    nps = len(c)
    scr = G.make_script(S['src'], S['line'], S['p'], 1.0, 1.0, c, nturn, {}, {}, by_text=True)
    t0 = time.time()
    txt, f = G.run_gs(scr, timeout=1800)
    el = time.time() - t0
    rec4 = G.parse(txt, nps, nturn)
    left = G.check_gs_left()
    if np.isnan(rec4).any():
        print('gs の出力が不完全', txt[-600:], 'gs の残り:', left)
        return
    rec = np.zeros((nturn, nps, 6))
    rec[:, :, :4] = rec4
    g = (float(D.fr(S['nu'][0])), float(D.fr(S['nu'][1])))
    nu = D.record_tunes(rec, S['tw'], window=1, guess=g, halfwidth=0.02, orbit=None)
    nu = D.pair_average(nu)
    r = D.adts_from_tunes(jx, jy, nu, order=2, nu_lin=g)
    fg = {'xx': r['xx'], 'xy': 0.5 * (r['xy'] + r['yx']), 'yy': r['yy']}
    trp = D.direct_adts(S['lat'], S['p'], jmax, order=2, nturn=nturn)
    fp = {'xx': trp['xx'], 'xy': trp['xy_mean'], 'yy': trp['yy']}
    out = ['リング %s: gs (SAD K64、WSL、NORFSW;)、6 極 ×%g、%d 周、Jmax=%.1e、%d 粒子、gs の経過 %.0f s、gs の残り: %s' % (key, lam, nturn, jmax, nps, el, left or 'なし'),
           '  解析 R2 (偏向磁石の 3 次なし)  %s' % f3(A0.total), '  解析 R3 (あり)                 %s' % f3(A1.total),
           '  gs の直接のフィット            %s' % f3(fg), '  pysad k64 (同じ粒子)           %s' % f3(fp),
           '  gs と pysad の相対差           %s' % fr3(rel(fg, fp)),
           '  解析 R3 − gs の相対差          %s    (R2: %s)' % (fr3(rel(A1.total, fg)), fr3(rel(A0.total, fg)))]
    s = '\n'.join(out)
    print(s)
    open(os.path.join(HERE, 'R3_adts_gs_%s_lam%g.gs.txt' % (key, lam)), 'w', encoding='utf-8').write(s + '\n')


def cmd_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'Meiryo', 'MS Gothic', 'DejaVu Sans']
    names = [('fodo', 'FODO'), ('fodo_rect', 'FODO 長方形'), ('fodo_asym', 'FODO E1≠E2'), ('astrid', 'ASTRID2'), ('astrid_mpw', 'ASTRID2 MPW'), ('kekls', 'KEK-LS')]
    data = {}
    for k, nm in names:
        fn = os.path.join(HERE, 'R3_scan_%s.pysad.json' % k)
        if os.path.exists(fn):
            data[k] = json.load(open(fn))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), sharey=False)
    for ic, comp in enumerate(KEYS):
        a = ax[ic]
        labels, v2, v3 = [], [], []
        for k, nm in names:
            if k not in data:
                continue
            for row in data[k]['rows']:
                if abs(row['lam'] - 1.0) < 1e-9:
                    labels.append(nm)
                    v2.append(100 * row['rel_R2'][comp])
                    v3.append(100 * row['rel_R3'][comp])
        x = np.arange(len(labels))
        a.bar(x - 0.2, v2, 0.4, color='#9ecae1', label='R2 (偏向磁石の 3 次なし)')
        a.bar(x + 0.2, v3, 0.4, color='#d62728', label='R3 (あり)')
        a.set_xticks(x)
        a.set_xticklabels(labels, fontsize=8)
        a.axhline(0, color='k', lw=0.5)
        a.set_title('∂ν%s/∂J%s: (解析 − 追跡)/追跡 [%%]' % (('x', 'x') if comp == 'xx' else ('x', 'y') if comp == 'xy' else ('y', 'y')))
        if ic == 0:
            a.legend(fontsize=8)
    fig.suptitle('R3: 解析の ADTS と追跡 (pysad k64 + NAFF + フィット) の相対差、6 極の強さ ×1', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'R3_adts_bend.png'), dpi=110)
    # 6 極の強さ走査
    fig, ax = plt.subplots(len(data), 3, figsize=(13, 2.7 * len(data)), squeeze=False)
    for ir, (k, nm) in enumerate([(k, nm) for k, nm in names if k in data]):
        lam = [r['lam'] for r in data[k]['rows']]
        for ic, comp in enumerate(KEYS):
            a = ax[ir, ic]
            a.plot(lam, [100 * r['rel_R2'][comp] for r in data[k]['rows']], 'o--', color='#6baed6', label='R2')
            a.plot(lam, [100 * r['rel_R3'][comp] for r in data[k]['rows']], 's-', color='#d62728', label='R3')
            a.axhline(0, color='k', lw=0.5)
            a.set_title('%s  %s' % (nm, comp), fontsize=9)
            a.set_xlabel('6 極の強さ λ')
            a.set_ylabel('相対差 [%]')
            if ir == 0 and ic == 0:
                a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'R3_adts_scan.png'), dpi=110)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'plot':
        cmd_plot()
        sys.exit(0)
    key = sys.argv[2]
    jm = float(sys.argv[sys.argv.index('--jmax') + 1]) if '--jmax' in sys.argv else RINGS[key][4]
    nt = int(sys.argv[sys.argv.index('--turns') + 1]) if '--turns' in sys.argv else 1024
    ns = int(sys.argv[sys.argv.index('--nslice') + 1]) if '--nslice' in sys.argv else 16
    if cmd == 'scan':
        lams = [float(v) for v in sys.argv[sys.argv.index('--lams') + 1].split(',')] if '--lams' in sys.argv else [0.0, 0.5, 1.0, 2.0]
        cmd_scan(key, lams, jm, nt, ns)
    elif cmd == 'conv':
        cmd_conv(key)
    elif cmd == 'gs':
        lam = float(sys.argv[sys.argv.index('--lam') + 1]) if '--lam' in sys.argv else 1.0
        cmd_gs(key, lam, jm, nt, ns)
