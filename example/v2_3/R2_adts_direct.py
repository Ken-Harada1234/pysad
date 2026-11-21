# -*- coding: utf-8 -*-
"""v2.3.0 R2: 実際の強さの追跡 + NAFF のチューンを J の多項式でフィットして ADTS を直接求め、resdrive の解析値と比べる (pysad の追跡、gs の追跡)。
Q2 の比較 (6 極を ×(+1,-1,0) にした偶数次の差) と違い、強さは変えない。解析側は adts() の 4 つの和 (6 極の 2 次 + 8 極の 1 次 + 運動学的な項 + 端部のフリンジ)。

使い方 (リポジトリ直下から)
  py -3.12 example/v2_3/R2_adts_direct.py pysad fodo|astrid|astrid_mpw|kekls [--jmax J] [--turns 1024] [--study]
  py -3.12 example/v2_3/R2_adts_direct.py gs    fodo|astrid|astrid_mpw|kekls [--jmax J] [--turns 1024]        (WSL の gs。NORFSW;。timeout 1800 s。1 本ずつ)
  py -3.12 example/v2_3/R2_adts_direct.py plot  fodo|astrid|kekls                                              (R2_fit_<ring>.png)
条件: 3 GeV 等はラティスのまま、RF 空洞は止める (pysad は rfsw=False、gs は NORFSW)、放射なし、結合なし、δ=0 (閉軌道は 0)。
J の格子: Jx, Jy ∈ {0, 1/5, ..., 1}·Jmax の 35 点 (Jx=Jy=0 を除く) × (φ, φ+π) の対 = 70 粒子 (x だけ、y だけ、両方を含む)。位相 φx = φy = 0 (対は π)。
NAFF: Laskar の Hann 窓 (p=1)、探索は線形のチューン ± 0.02、既定 1024 周。フィット: νx (Jx>0 の点)・νy (Jy>0 の点) を全次数 ≤ order の多項式で (既定 order=2)。"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, adtsdirect as D, compat
from pysad.parser import load_sad_file, load_sad_string
compat.set_mode(os.environ.get('PYSAD_COMPAT', 'k64'))
import Q2_adts_gs as G

# リング: (ファイル, ライン, 運動量 eV, FODO の 6 極の上書き, 既定の Jmax)
RINGS = {
    'fodo': (os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), 'RING', 3e9, {'SF': 4.0, 'SD': -6.0}, 1e-7),
    'astrid': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING', 0.58e9, None, 2e-8),
    'astrid_mpw': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), 'RINGMPW', 0.58e9, None, 2e-8),
    'kekls': (os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), 'RING', 3e9, None, 2e-10),
}
NSTEP = 5


def setup(key):
    path, line, p, k2over, jb = RINGS[key]
    src = open(path, encoding='utf-8').read()
    if k2over:
        src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k2over['SF'], k2over['SD']))
    lat = load_sad_string(src).build_lattice(line)
    rows = lat.twiss_table()
    r0 = rows[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    nu = (float(rows[-1].mux) / (2 * math.pi), float(rows[-1].muy) / (2 * math.pi))
    g = (float(D.fr(nu[0])), float(D.fr(nu[1])))
    return dict(key=key, src=src, line=line, p=p, lat=lat, rows=rows, tw=tw, nu=nu, g=g, jb=jb, circ=float(rows[-1].s))


def grid(S, jmax, phix=0.0, phiy=0.0):
    jx, jy = D.grid_actions(jmax, NSTEP)
    JX, JY = np.r_[jx, jx], np.r_[jy, jy]
    PX = np.r_[np.full(len(jx), phix), np.full(len(jx), phix + math.pi)]
    PY = np.r_[np.full(len(jx), phiy), np.full(len(jx), phiy + math.pi)]
    return jx, jy, D.initial_coords(S['tw'], JX, JY, PX, PY)


def analytic(S):
    A = R.adts(S['lat'], rows=S['rows'], bend_cubic=False)
    return A, {'sext': dict(A.sext), 'octu': dict(A.octu), 'kin': dict(A.kin), 'fringe': dict(A.fringe), 'total': dict(A.total)}


def fit_pair(S, jx, jy, rec, nturn=1024, window=1, order=2, skip=0, halfwidth=0.02, pair=True):
    nu = D.record_tunes(rec, S['tw'], window=window, guess=S['g'], halfwidth=halfwidth, nturn=nturn, skip=skip)
    if pair:
        nu = D.pair_average(nu)
    else:
        nu = nu[:len(jx)]
    return D.adts_from_tunes(jx, jy, nu, order=order, nu_lin=S['g']), nu


def rel(a, b):
    return (a - b) / abs(b)


def line3(name, d):
    return '%-34s xx=% .5f  xy=% .5f  yy=% .5f' % (name, d['xx'], d['xy'], d['yy'])


def fmt_r(r):
    return {'xx': r['xx'], 'xy': 0.5 * (r['xy'] + r['yx']), 'yy': r['yy']}


def report(S, A, ad, r, label, out):
    f = fmt_r(r)
    out.append(line3(label, f))
    out.append('   解析の全部入り に対する相対差: xx=%.2e xy=%.2e yy=%.2e (xy は ∂νx/∂Jy と ∂νy/∂Jx の平均、両者の差 %.1e、フィットの rms %.1e %.1e)' % (
        rel(f['xx'], ad['total']['xx']), rel(f['xy'], ad['total']['xy']), rel(f['yy'], ad['total']['yy']),
        abs(r['xy'] - r['yx']) / abs(f['xy']), r['rms'][0], r['rms'][1]))


def cmd_pysad(key, jmax, nturn, study):
    S = setup(key)
    out = []
    A, ad = analytic(S)
    out.append('リング %s: p=%.4g GeV/c、周長 %.3f m、νx=%.6f νy=%.6f、k64 の追跡 (rfsw=False、放射なし、δ=0)、J の格子 Jmax=%.1e (%d×%d、±位相の対 = %d 粒子)、%d 周、NAFF Hann p=1、order=2' % (
        key, S['p'] / 1e9, S['circ'], S['nu'][0], S['nu'][1], jmax, NSTEP + 1, NSTEP + 1, 2 * (len(D.grid_actions(jmax, NSTEP)[0])), nturn))
    for w in A.warnings:
        out.append('(解析の注意) ' + w)
    for nm in ('sext', 'octu', 'kin', 'fringe', 'total'):
        out.append(line3('解析 ' + {'sext': '6 極 2 次', 'octu': '8 極 1 次', 'kin': '運動学的 p^4/8', 'fringe': '端部のフリンジ', 'total': '合計'}[nm], ad[nm]))
    t0 = time.time()
    jx, jy, c = grid(S, jmax)
    nmax = max(nturn, 4096 if study else nturn)
    rec = D.track_record(S['lat'].elements, S['p'], c, nmax)
    out.append('(全部入りの追跡 %d 周 %.0f s)' % (nmax, time.time() - t0))
    r, nu = fit_pair(S, jx, jy, rec, nturn=nturn)
    report(S, A, ad, r, '追跡 (pysad k64) の直接のフィット', out)
    f = fmt_r(r)
    out.append('   解析との差 (追跡 − 合計) = xx %.4f xy %.4f yy %.4f  (入れていない項: 偏向磁石の曲率・端部、6 極×(その他) の交差項、J の高次)' % (
        f['xx'] - ad['total']['xx'], f['xy'] - ad['total']['xy'], f['yy'] - ad['total']['yy']))
    out.append('   ν0 のフィット (端数) = %.9f %.9f、線形のチューン = %.9f %.9f' % (r['nu0'][0], r['nu0'][1], S['g'][0], S['g'][1]))
    res = {'ring': key, 'p': S['p'], 'circ': S['circ'], 'nu': S['nu'], 'jmax': jmax, 'nturn': nturn, 'analytic': ad,
           'pysad': {'xx': f['xx'], 'xy': f['xy'], 'yy': f['yy'], 'xy_yx': r['yx']}, 'jx': list(map(float, jx)), 'jy': list(map(float, jy)),
           'nu_pair': nu.tolist(), 'g': list(S['g'])}
    # パラキシャルな追跡 (解析式と同じ模型) との照合
    try:
        nsl = 40 if key in ('fodo', 'astrid', 'astrid_mpw') else 6
        tr = D.paraxial_tracker(S['lat'], nsl)
        recp = D.track_record(S['lat'].elements, S['p'], c, nturn, tracker=tr)
        rp, _ = fit_pair(S, jx, jy, recp, nturn=nturn)
        fp = fmt_r(rp)
        sx = {k: ad['sext'][k] + ad['octu'][k] for k in ('xx', 'xy', 'yy')}
        out.append(line3('追跡 (線形行列 + 薄い 6 極・8 極、パラキシャル)', fp))
        out.append('   解析の 6 極 2 次 + 8 極 1 次 に対する相対差: xx=%.2e xy=%.2e yy=%.2e (薄レンズ %d 分割の誤差を含む)' % (
            rel(fp['xx'], sx['xx']), rel(fp['xy'], sx['xy']), rel(fp['yy'], sx['yy']), nsl))
        res['paraxial'] = fp
    except NotImplementedError as e:
        out.append('(パラキシャルな追跡は省略: %s)' % e)
    if study:
        out.append('')
        out.append('--- フィットの次数と J の範囲 (全部入りの追跡、解析の合計に対する相対差) ---')
        out.append('  Jmax        order   xx          xy          yy         rms(νx)  rms(νy)')
        res['range'] = []
        for fac in (0.01, 0.1, 1.0, 10.0, 100.0):
            jm = jmax * fac
            jx2, jy2, c2 = grid(S, jm)
            try:
                rec2 = D.track_record(S['lat'].elements, S['p'], c2, nturn)
            except FloatingPointError:
                out.append('  %.1e  発散' % jm)
                continue
            for order in (1, 2, 3, 4):
                r2, _ = fit_pair(S, jx2, jy2, rec2, nturn=nturn, order=order)
                f2 = fmt_r(r2)
                out.append('  %.1e   %d     %+.2e  %+.2e  %+.2e   %.1e  %.1e' % (
                    jm, order, rel(f2['xx'], ad['total']['xx']), rel(f2['xy'], ad['total']['xy']), rel(f2['yy'], ad['total']['yy']), r2['rms'][0], r2['rms'][1]))
                res['range'].append({'jmax': jm, 'order': order, 'xx': f2['xx'], 'xy': f2['xy'], 'yy': f2['yy']})
        out.append('')
        out.append('--- NAFF の設定 (Jmax=%.1e、order=2、全部入りの追跡、解析の合計に対する相対差) ---' % jmax)
        out.append('  設定                                 xx          xy          yy')
        res['naff'] = []
        for nt in (128, 256, 512, 1024, 2048, 4096):
            if nt > nmax:
                continue
            r2, _ = fit_pair(S, jx, jy, rec, nturn=nt)
            f2 = fmt_r(r2)
            out.append('  周数 %-5d (p=1, 位相 0)             %+.2e  %+.2e  %+.2e' % (nt, rel(f2['xx'], ad['total']['xx']), rel(f2['xy'], ad['total']['xy']), rel(f2['yy'], ad['total']['yy'])))
            res['naff'].append({'what': 'turns', 'v': nt, 'xx': f2['xx'], 'xy': f2['xy'], 'yy': f2['yy']})
        for wp in (1, 2, 3, 4):
            r2, _ = fit_pair(S, jx, jy, rec, nturn=1024, window=wp)
            f2 = fmt_r(r2)
            out.append('  窓の次数 p=%d (1024 周)                %+.2e  %+.2e  %+.2e' % (wp, rel(f2['xx'], ad['total']['xx']), rel(f2['xy'], ad['total']['xy']), rel(f2['yy'], ad['total']['yy'])))
            res['naff'].append({'what': 'window', 'v': wp, 'xx': f2['xx'], 'xy': f2['xy'], 'yy': f2['yy']})
        for sk in (0, 1, 7, 100):
            r2, _ = fit_pair(S, jx, jy, rec, nturn=1024, skip=sk)
            f2 = fmt_r(r2)
            out.append('  先頭の %-3d 周を捨てる (1024 周)        %+.2e  %+.2e  %+.2e' % (sk, rel(f2['xx'], ad['total']['xx']), rel(f2['xy'], ad['total']['xy']), rel(f2['yy'], ad['total']['yy'])))
            res['naff'].append({'what': 'skip', 'v': sk, 'xx': f2['xx'], 'xy': f2['xy'], 'yy': f2['yy']})
        for (px, py, pr) in ((math.pi / 2, math.pi / 2, True), (0.7, 2.1, True), (0.0, 0.0, False), (1.3, 0.4, False)):
            jx3, jy3, c3 = grid(S, jmax, px, py)
            rec3 = D.track_record(S['lat'].elements, S['p'], c3, 1024)
            r3, _ = fit_pair(S, jx3, jy3, rec3, nturn=1024, pair=pr)
            f3 = fmt_r(r3)
            out.append('  初期位相 (%.2f, %.2f) %s              %+.2e  %+.2e  %+.2e' % (px, py, '±対' if pr else '単独', rel(f3['xx'], ad['total']['xx']), rel(f3['xy'], ad['total']['xy']), rel(f3['yy'], ad['total']['yy'])))
            res['naff'].append({'what': 'phase%s' % ('pair' if pr else 'single'), 'v': [px, py], 'xx': f3['xx'], 'xy': f3['xy'], 'yy': f3['yy']})
        for hw in (0.005, 0.02, 0.1):
            r2, _ = fit_pair(S, jx, jy, rec, nturn=1024, halfwidth=hw)
            f2 = fmt_r(r2)
            out.append('  探索幅 ±%-5g (1024 周)               %+.2e  %+.2e  %+.2e' % (hw, rel(f2['xx'], ad['total']['xx']), rel(f2['xy'], ad['total']['xy']), rel(f2['yy'], ad['total']['yy'])))
    txt = '\n'.join(out)
    print(txt)
    open(os.path.join(HERE, 'R2_adts_direct_%s.pysad.txt' % key), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(res, open(os.path.join(HERE, 'R2_adts_direct_%s.pysad.json' % key), 'w'), indent=1)


def cmd_gs(key, jmax, nturn):
    S = setup(key)
    out = []
    A, ad = analytic(S)
    jx, jy, c = grid(S, jmax)
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
    r, nu = fit_pair(S, jx, jy, rec, nturn=nturn)
    # 同じ粒子の pysad の追跡と、チューンの粒子ごとの差
    recp = D.track_record(S['lat'].elements, S['p'], c, nturn)
    rp, nup = fit_pair(S, jx, jy, recp, nturn=nturn)
    dnu = np.abs(D.fr(nu - nup)).max()
    out.append('リング %s: gs (SAD K64、WSL、NORFSW;)、%d 周、Jmax=%.1e、%d 粒子、gs の経過 %.0f s、gs の残り: %s' % (key, nturn, jmax, nps, el, left or 'なし'))
    out.append(line3('解析 合計', ad['total']))
    report(S, A, ad, r, 'gs の直接のフィット', out)
    f_gs, f_py = fmt_r(r), fmt_r(rp)
    out.append(line3('pysad (同じ粒子)', f_py))
    out.append('   gs と pysad の相対差: xx=%.2e xy=%.2e yy=%.2e;  粒子ごとのチューンの最大差 %.2e' % (
        rel(f_gs['xx'], f_py['xx']), rel(f_gs['xy'], f_py['xy']), rel(f_gs['yy'], f_py['yy']), dnu))
    s = '\n'.join(out)
    print(s)
    open(os.path.join(HERE, 'R2_adts_direct_%s.gs.txt' % key), 'w', encoding='utf-8').write(s + '\n')
    json.dump({'ring': key, 'jmax': jmax, 'nturn': nturn, 'gs': f_gs, 'pysad': f_py, 'analytic': ad, 'jx': list(map(float, jx)), 'jy': list(map(float, jy)),
               'nu_gs': nu.tolist(), 'nu_pysad': nup.tolist(), 'g': list(S['g'])}, open(os.path.join(HERE, 'R2_adts_direct_%s.gs.json' % key), 'w'), indent=1)


def cmd_plot(key):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'Meiryo', 'MS Gothic', 'DejaVu Sans']
    d = json.load(open(os.path.join(HERE, 'R2_adts_direct_%s.pysad.json' % key)))
    jx, jy = np.array(d['jx']), np.array(d['jy'])
    nu = np.array(d['nu_pair'])
    g = d['g']
    nu = nu.copy()
    nu[:, 0] = D.fr(nu[:, 0] - g[0])
    nu[:, 1] = D.fr(nu[:, 1] - g[1])
    r = D.adts_from_tunes(jx, jy, nu, order=3)
    fx, fy = r['fitx'], r['fity']
    tot = d['analytic']['total']
    jm = d['jmax']

    def poly(fit, a, b):
        v = np.zeros_like(a)
        for (p_, q_), cc in fit['c'].items():
            v = v + cc * a ** p_ * b ** q_
        return v
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.5), gridspec_kw={'height_ratios': [3, 1.3]})
    cases = [('Δνx vs Jx (Jy=0)', 0, jy == 0, jx, 'xx'), ('Δνy vs Jy (Jx=0)', 1, jx == 0, jy, 'yy'), ('Δνx vs Jy (Jx=Jmax)', 0, np.isclose(jx, jm, rtol=1e-6, atol=0), jy, 'xy')]
    for k, (ttl, comp, m, xs, key2) in enumerate(cases):
        a = ax[0, k]
        o = np.argsort(xs[m])
        xx = xs[m][o]
        yy = nu[m, comp][o]
        fit = fx if comp == 0 else fy
        xf = np.linspace(0, xx.max(), 100)
        if k == 0:
            pf = poly(fit, xf, 0 * xf)
            pd = poly(fit, xx, 0 * xx)
            ref = nu[m, comp][o]
            lin = d['analytic']['total']['xx'] * xf
        elif k == 1:
            pf = poly(fit, 0 * xf, xf)
            pd = poly(fit, 0 * xx, xx)
            lin = d['analytic']['total']['yy'] * xf
        else:
            pf = poly(fit, jm + 0 * xf, xf)
            pd = poly(fit, jm + 0 * xx, xx)
            lin = d['analytic']['total']['xx'] * jm + d['analytic']['total']['xy'] * xf
        a.plot(xx, yy, 'o', color='#1f77b4', label='追跡 + NAFF (pysad)')
        a.plot(xf, pf, '-', color='#d62728', label='3 次のフィット')
        a.plot(xf, lin, '--', color='#2ca02c', label='解析の接線 (全部入り)')
        a.set_title(ttl)
        a.set_ylabel('Δν (線形のチューンからの差)')
        a.legend(fontsize=7)
        ax[1, k].plot(xx, yy - pd, 'o-', color='#555')
        ax[1, k].set_ylabel('残差')
        ax[1, k].set_xlabel('J [m rad]')
        ax[1, k].axhline(0, color='k', lw=0.5)
    fig.suptitle('R2: %s  ν(J) の追跡 + NAFF と多項式フィット (Jmax=%.1e、%d 周、全部入りの解析 xx=%.2f xy=%.2f yy=%.2f 1/(m rad))' % (
        key, jm, d['nturn'], tot['xx'], tot['xy'], tot['yy']), fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'R2_fit_%s.png' % key), dpi=110)


if __name__ == '__main__':
    cmd, key = sys.argv[1], sys.argv[2]
    jm = float(sys.argv[sys.argv.index('--jmax') + 1]) if '--jmax' in sys.argv else RINGS[key][4]
    nt = int(sys.argv[sys.argv.index('--turns') + 1]) if '--turns' in sys.argv else 1024
    if cmd == 'pysad':
        cmd_pysad(key, jm, nt, '--study' in sys.argv)
    elif cmd == 'gs':
        cmd_gs(key, jm, nt)
    elif cmd == 'plot':
        cmd_plot(key)
