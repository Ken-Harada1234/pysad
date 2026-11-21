# -*- coding: utf-8 -*-
"""v2.3.0 R2: ADTS・クロマティシティの「未確認の部分」を追跡 (pysad k64 + NAFF) で調べる。節ごとに実行 (出力は R2_<節>.pysad.txt)。
  py -3.12 example/v2_3/R2_unverified.py octupole   # 8 極の強さ μ 倍の走査: 1 次の ADTS が μ に比例するか、2 次 (J の 2 次) の係数の μ 依存
  py -3.12 example/v2_3/R2_unverified.py ablate     # 端部・運動学的な項・複合機能偏向磁石を切った追跡との比較 (DISFRIN=1、E1=E2=0、K1=0)
  py -3.12 example/v2_3/R2_unverified.py delta      # チューンの δ 依存 (∂ν/∂δ、∂²ν/∂δ²) とオフモメンタムの ADTS (δ = ±δ0)
  py -3.12 example/v2_3/R2_unverified.py resonance  # 3νx 共鳴線の近くの不安定固定点の振幅: h30000 + ADTS の予測と追跡 (Newton 法)
  py -3.12 example/v2_3/R2_unverified.py coupling   # 結合 (四極の回転)・垂直分散・DX/DY を無視した場合の ADTS の誤差
条件はいずれも RF 空洞は止める (rfsw=False)、放射なし、k64、1024 周、NAFF Hann p=1、J の格子 6×6 × ±位相の対 (adtsdirect.direct_adts)。"""
import copy, dataclasses, json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, adtsdirect as D, compat
compat.set_mode('k64')
import R2_adts_direct as M

OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def save(name):
    open(os.path.join(HERE, 'R2_%s.pysad.txt' % name), 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def f3(d):
    return 'xx=% .4f xy=% .4f yy=% .4f' % (d['xx'], d['xy'], d['yy'])


def dres(r):
    return {'xx': r['xx'], 'xy': r['xy_mean'], 'yy': r['yy']}


def relv(a, b):
    return {k: (a[k] - b[k]) / abs(b[k]) if b[k] else float('nan') for k in ('xx', 'xy', 'yy')}


def fmt_rel(d):
    return 'xx=%+.2e xy=%+.2e yy=%+.2e' % (d['xx'], d['xy'], d['yy'])


def octupole(key='kekls'):
    S = M.setup(key)
    lat = S['lat']
    oc = [e for e in lat.elements if type(e).__name__ == 'Octupole']
    k3 = [e.K3 for e in oc]
    P('節 octupole: リング %s (p=%.3g GeV/c、周長 %.1f m)。8 極 %d 台の K3 を μ 倍。直接のフィット (order=3) の 1 次の係数 a = ∂ν/∂J と 2 次の係数 b = (1/2)∂²ν/∂J² (νx = ν0 + a_xx Jx + a_xy Jy + b_xx Jx² + b_xy JxJy + b_yy Jy² の係数)' % (key, S['p'] / 1e9, S['circ'], len(oc)))
    P('条件: Jmax=%.1e、order=3、1024 周、k64、RF 止め' % (S['jb'] * 10))
    res = []
    for mu in (0.0, 1.0, 2.0, 4.0, 8.0):
        for e, k in zip(oc, k3):
            e.K3 = k * mu
        A = R.adts(lat, bend_cubic=False)
        r = D.direct_adts(lat, S['p'], S['jb'] * 10, order=3, nturn=1024)
        fx = r['fitx']['c']
        P('μ=%g: 解析の合計 %s' % (mu, f3(A.total)))
        P('      追跡の a            %s   相対差 %s' % (f3(dres(r)), fmt_rel(relv(dres(r), A.total))))
        P('      追跡の b (νx)       b_xx=%.4e b_xy=%.4e b_yy=%.4e   (νy) b_xx=%.4e b_xy=%.4e b_yy=%.4e' % (
            fx[(2, 0)], fx[(1, 1)], fx[(0, 2)], r['fity']['c'][(2, 0)], r['fity']['c'][(1, 1)], r['fity']['c'][(0, 2)]))
        res.append((mu, dres(r), A.total, fx[(2, 0)], fx[(1, 1)], fx[(0, 2)]))
    for e, k in zip(oc, k3):
        e.K3 = k
    mus = np.array([r[0] for r in res])
    P('')
    P('b_xx(μ) の 2 次多項式フィット b = b0 + b1 μ + b2 μ²: (b2 = 8 極の 2 次、b1 = 6 極×8 極、b0 = 6 極の 4 次ほか)')
    for nm, idx in (('b_xx', 3), ('b_xy', 4), ('b_yy', 5)):
        v = np.array([r[idx] for r in res])
        co = np.polynomial.polynomial.polyfit(mus, v, 2)
        P('  %s: b0=%.4e b1=%.4e b2=%.4e (最大残差 %.1e)' % (nm, co[0], co[1], co[2], np.abs(v - np.polynomial.polynomial.polyval(mus, co)).max()))
    for nm in ('xx', 'xy', 'yy'):
        a = np.array([r[1][nm] for r in res])
        co = np.polynomial.polynomial.polyfit(mus, a, 1)
        P('  a_%s(μ) の直線フィット: a0=%.4f 傾き=%.4f (解析の 8 極 1 次の傾き %.4f)' % (nm, co[0], co[1], (res[1][2][nm] - res[0][2][nm])))
    save('octupole_' + key)


def ablate(key):
    S = M.setup(key)
    lat0 = S['lat']
    P('節 ablate: リング %s (p=%.3g GeV/c、周長 %.1f m)。解析の合計 (6 極 2 次 + 8 極 1 次 + 運動学的 + 端部) と追跡を、切った場合と比べる' % (key, S['p'] / 1e9, S['circ']))
    A = R.adts(lat0, bend_cubic=False)
    P('  解析: 6 極 %s' % f3(A.sext))
    P('        8 極 %s' % f3(A.octu))
    P('        運動学 %s' % f3(A.kin))
    P('        端部   %s' % f3(A.fringe))
    P('        合計   %s' % f3(A.total))
    jm = S['jb'] * 5
    variants = {}

    def variant(name, fn):
        lat = copy.deepcopy(lat0)
        fn(lat)
        variants[name] = lat

    def nofringe(lat):
        for e in lat.elements:
            if type(e).__name__ in ('Quadrupole', 'Bend'):
                e.DISFRIN = 1.0

    def noedge(lat):
        for e in lat.elements:
            if type(e).__name__ == 'Bend':
                e.E1 = 0.0
                e.E2 = 0.0

    def nosext(lat):
        for e in lat.elements:
            if type(e).__name__ == 'Sextupole':
                e.K2 = 0.0

    def nosextoct(lat):
        for e in lat.elements:
            if type(e).__name__ == 'Sextupole':
                e.K2 = 0.0
            if type(e).__name__ == 'Octupole':
                e.K3 = 0.0
    variant('全部入り (K64 の既定)', lambda l: None)
    variant('端部のフリンジ切り (DISFRIN=1: 四極・偏向磁石)', nofringe)
    variant('6 極・8 極なし (運動学的 + 端部 + 偏向磁石)', nosextoct)
    variant('6 極・8 極・端部のフリンジなし (運動学的 + 偏向磁石の曲率)', lambda l: (nosextoct(l), nofringe(l)))
    if any(type(e).__name__ == 'Bend' and (e.E1 or e.E2) for e in lat0.elements):
        variant('偏向磁石の端 E1=E2=0 (セクター型)', noedge)
    P('')
    P('直接のフィット (Jmax=%.1e、order=2、1024 周):' % jm)
    for name, lat in variants.items():
        r = D.direct_adts(lat, S['p'], jm, order=2)
        A2 = R.adts(lat, bend_cubic=False)
        t = dres(r)
        P('%s' % name)
        P('   追跡   %s' % f3(t))
        if '端部のフリンジ切り' in name or 'フリンジなし' in name:
            an = {k: A2.sext[k] + A2.octu[k] + A2.kin[k] for k in ('xx', 'xy', 'yy')}
            P('   解析   %s (6 極 2 次 + 8 極 1 次 + 運動学的、端部なし)' % f3(an))
        else:
            an = A2.total
            P('   解析   %s (全部入り)' % f3(an))
        P('   差 (追跡 − 解析) %s   相対 %s' % (f3({k: t[k] - an[k] for k in t}), fmt_rel(relv(t, an))))
    save('ablate_' + key)


def delta_sec(key, d0, jmax):
    S = M.setup(key)
    lat = S['lat']
    P('節 delta: リング %s (p=%.3g GeV/c、周長 %.1f m)。RF 止め、δ は一定の運動量偏差、閉軌道から出す' % (key, S['p'] / 1e9, S['circ']))
    P('(1) チューンの δ 依存: 微小振幅 (J=1e-12) の NAFF、δ = ±%g の Chebyshev 節点 11 点、次数 5 の多項式。CALC 数値微分 (chromaticity_numeric) と 6 極の 1 次の解析式 ξ = h/π を比較' % d0)
    nodes = d0 * np.cos(np.pi * (np.arange(11) + 0.5) / 11)
    nut, nuref = D.tune_vs_delta(lat, S['p'], nodes)
    for p_ in (0, 1):
        v = D.fr(nut[:, p_] - nuref[:, p_]) + (nuref[:, p_] - np.round(nuref[:, p_]))
    out = {}
    cn = R.chromaticity_numeric(lat, order=3, dmax=d0)
    sx = R.chromaticity_sextupole(lat)
    for p_, pn in enumerate('xy'):
        val = D.fr(nut[:, p_] - nuref[:, p_]) + D.fr(nuref[:, p_])
        # 追跡のチューン (符号つき端数) の連続化: 線形のチューンの端数との差 + CALC の端数
        fit = D.fit_delta_curve(nodes, val, 5)
        cal = D.fit_delta_curve(nodes, D.fr(nuref[:, p_]), 5)
        out[pn] = (fit, cal)
        P('  ν%s: 追跡 ν0=%.7f ∂ν/∂δ=%.5f ∂²ν/∂δ²=%.3f ∂³ν/∂δ³=%.1f (rms %.1e) | CALC のチューン ν0=%.7f ∂ν/∂δ=%.5f ∂²ν/∂δ²=%.3f ∂³=%.1f | chromaticity_numeric d1=%.5f d2=%.3f d3=%.1f' % (
            pn, fit['nu0'], fit['d1'], fit['d2'], fit['d3'], fit['rms'], cal['nu0'], cal['d1'], cal['d2'], cal['d3'], cn['d1'][p_], cn['d2'][p_], cn['d3'][p_]))
    P('  6 極の 1 次の寄与 (ξx, ξy) = (%.5f, %.5f)' % sx)
    P('  追跡の ∂ν/∂δ と CALC の差: x %.2e y %.2e (規格化せず絶対値)' % (out['x'][0]['d1'] - out['x'][1]['d1'], out['y'][0]['d1'] - out['y'][1]['d1']))
    P('  追跡の ∂²ν/∂δ² と CALC の差: x %.2e y %.2e' % (out['x'][0]['d2'] - out['x'][1]['d2'], out['y'][0]['d2'] - out['y'][1]['d2']))
    # 6 極なしの線形のクロマティシティ
    lat_ns = copy.deepcopy(lat)
    for e in lat_ns.elements:
        if type(e).__name__ == 'Sextupole':
            e.K2 = 0.0
        elif type(e).__name__ == 'Octupole':
            e.K3 = 0.0
    cn0 = R.chromaticity_numeric(lat_ns, order=3, dmax=d0)
    P('  6 極・8 極なしの CALC: ∂ν/∂δ = (%.5f, %.5f)、∂²ν/∂δ² = (%.3f, %.3f) (全部入りの差 = 6 極の寄与: ∂ν/∂δ %.5f, %.5f; ∂²ν/∂δ² %.3f, %.3f)' % (
        cn0['d1'][0], cn0['d1'][1], cn0['d2'][0], cn0['d2'][1], cn['d1'][0] - cn0['d1'][0], cn['d1'][1] - cn0['d1'][1], cn['d2'][0] - cn0['d2'][0], cn['d2'][1] - cn0['d2'][1]))
    P('')
    P('(2) オフモメンタムの ADTS: δ = 0, ±%g で adts(delta=δ) の合計と追跡の直接のフィット (Jmax=%.1e、order=2)。差分 ∂(∂ν/∂J)/∂δ = (a(+δ0) − a(−δ0))/(2δ0)' % (d0, jmax))
    tr, an = {}, {}
    for d in (-d0, 0.0, d0):
        r = D.direct_adts(lat, S['p'], jmax, delta=d, order=2)
        A = R.adts(lat, delta=d, bend_cubic=False)
        tr[d], an[d] = dres(r), dict(A.total)
        P('  δ=%+g 追跡 %s | 解析 合計 %s | 6 極 2 次 %s | 相対差 %s' % (d, f3(tr[d]), f3(an[d]), f3(A.sext), fmt_rel(relv(tr[d], an[d]))))
    dtr = {k: (tr[d0][k] - tr[-d0][k]) / (2 * d0) for k in ('xx', 'xy', 'yy')}
    dan = {k: (an[d0][k] - an[-d0][k]) / (2 * d0) for k in ('xx', 'xy', 'yy')}
    P('  ∂a/∂δ  追跡 %s' % f3(dtr))
    P('  ∂a/∂δ  解析 %s   相対差 %s' % (f3(dan), fmt_rel(relv(dtr, dan))))
    save('delta_' + key)


def resonance(eps_list=(3e-4, 6e-4, 1e-3, 2e-3)):
    """FODO の νx を 3νx = 20 (νx = 20/3) の近くに動かし (QF の K1 の調整)、y=0 の x の運動の 3 周写像の周期点 (Newton 法、位相 24 通りから) の振幅 J の最小・最大を、
    共鳴駆動項 h30000 と ADTS から予測する。1 周の生成関数 G = 2πνJ + 2π(a/2)J² + 2|h|(2J)^{3/2} cos(3φ+θ): 周期点の条件 2π(ε + aJ) = ∓6|h|√(2J)。
    ADTS の a は (1) 解析の全部入り、(2) 追跡の直接のフィット、の 2 通り (周期点の J は共鳴線 ν(J) = 20/3 の J = −ε/a から ±(3|h|√(2J)/π a) だけずれる)"""
    S = M.setup('fodo')
    src = S['src']
    from pysad.parser import load_sad_string
    from scipy.optimize import brentq, fsolve
    from pysad.track6n import track_particles_fast
    P('節 resonance: FODO (3 GeV、105.6 m)、6 極 K2L = 4, −6、QF の K1 を調整して νx = 20/3 + ε (3νx = 20 + 3ε)。y=0 の運動の 3 周写像の周期点 (24 通りの位相から Newton) の J の最小・最大')
    P('予測: 2π(ε + a J) = ∓6 |h30000| √(2J)。a = ADTS (xx)。J_res = −ε/a は共鳴線 (ADTS だけで ν(J) = 20/3 になる J)')
    for eps in eps_list:
        for sgn in (+1, -1):
            e0 = sgn * eps

            def nux_of(k1):
                s = src.replace('QF = (L=0.3 K1=0.7)', 'QF = (L=0.3 K1=%r)' % k1)
                lat = load_sad_string(s).build_lattice('RING')
                return lat, float(lat.twiss_table()[-1].mux) / (2 * math.pi)
            k1 = brentq(lambda k: nux_of(k)[1] - 20.0 / 3.0 - e0, 0.68, 0.72, xtol=1e-12)
            lat, nux = nux_of(k1)
            rows = lat.twiss_table()
            tw = ((float(rows[0].alfx), float(rows[0].betx)), (float(rows[0].alfy), float(rows[0].bety)))
            A = R.adts(lat, rows=rows, bend_cubic=False)
            h = R.resonance_terms(lat, ('h30000',), rows=rows, bend_cubic=False).total['h30000']      # R2 の値 (6 極のみ。R4 は R4_rdt_tracking.py res)
            ee = nux - 20.0 / 3.0
            a_an = A.total['xx']
            a_tr = D.direct_adts(lat, S['p'], 3e-7, order=2, nturn=512, nstep=3)['xx']

            def predict(a):
                out = []
                for sig in (+1, -1):
                    f = lambda u: 2 * math.pi * (ee + a * u * u / 2) + 6 * abs(h) * u * sig
                    try:
                        if f(1e-7) * f(0.3) < 0:
                            u = brentq(f, 1e-7, 0.3)
                            out.append(u * u / 2)
                    except ValueError:
                        pass
                return sorted(out)

            def T3(v):
                c = np.zeros((1, 6))
                c[0, 0], c[0, 1] = v
                for _ in range(3):
                    c = track_particles_fast(lat.elements, c, 1, S['p'], 0.0, False, False, False)
                return np.array([c[0, 0], c[0, 1]]) - v
            Jres = -ee / a_tr
            Jst = abs(ee) / a_tr
            Js = []
            for ph in np.arange(24) * math.pi / 12:
                al, be = tw[0]
                x0 = np.array([math.sqrt(2 * Jst * be) * math.cos(ph), -math.sqrt(2 * Jst / be) * (math.sin(ph) + al * math.cos(ph))])
                try:
                    sol, info, ier, msg = fsolve(T3, x0, full_output=True, xtol=1e-10)
                except Exception:
                    continue
                if np.abs(T3(sol)).max() < 1e-12:
                    x, px = sol
                    Js.append((x * x + (al * x + be * px) ** 2) / (2 * be))
            pa, pt = predict(a_an), predict(a_tr)
            P('ε=%+.2e (K1=%.7f νx=%.6f) |h30000|=%.5f a(解析)=%.2f a(追跡)=%.2f J_res=%.4e: 予測(解析 a) %s | 予測(追跡 a) %s | 追跡の周期点 J の最小・最大 %s' % (
                ee, k1, nux, abs(h), a_an, a_tr, Jres, ['%.4e' % v for v in pa], ['%.4e' % v for v in pt],
                ('%.4e %.4e (%d 個)' % (min(Js), max(Js), len(Js))) if Js else 'なし'))
    save('resonance')


def coupling():
    S = M.setup('astrid')
    lat0 = S['lat']
    P('節 coupling: ASTRID2 (0.58 GeV、45.7 m)。結合・垂直分散・DX/DY を入れて、解析 (結合なしの式) と追跡の直接のフィット (J は結合なしの Twiss から) を比べる')
    A0 = R.adts(lat0, bend_cubic=False)
    r0 = D.direct_adts(lat0, S['p'], 2e-8, order=2)
    P('基準 (誤差なし): 追跡 %s、解析 %s' % (f3(dres(r0)), f3(A0.total)))
    rng = np.random.default_rng(12345)
    quads = [e for e in lat0.elements if type(e).__name__ == 'Quadrupole']
    sexts = [e for e in lat0.elements if type(e).__name__ == 'Sextupole']
    for kind, amp in (('roll', 3e-4), ('roll', 1e-3), ('roll', 3e-3), ('dxdy_sext', 1e-4), ('dxdy_sext', 5e-4), ('dxdy_all', 1e-4)):
        lat = copy.deepcopy(lat0)
        qs = [e for e in lat.elements if type(e).__name__ == 'Quadrupole']
        ss = [e for e in lat.elements if type(e).__name__ == 'Sextupole']
        if kind == 'roll':
            for e in qs:
                e.ROTATION = float(rng.normal(0, amp))
        elif kind == 'dxdy_sext':
            for e in ss:
                e.DX, e.DY = float(rng.normal(0, amp)), float(rng.normal(0, amp))
        else:
            for e in qs + ss:
                e.DX, e.DY = float(rng.normal(0, amp)), float(rng.normal(0, amp))
        try:
            rows = lat.twiss_table()
            A = R.adts(lat, rows=rows, bend_cubic=False)
            tw, orb, g, nu = D.lattice_start(lat, 0.0, rows)
            rmax = max(max(abs(float(r_.r1)), abs(float(r_.r2)), abs(float(r_.r3)), abs(float(r_.r4))) for r_ in rows)
            ey = max(abs(float(r_.ety)) for r_ in rows)
            cod = max(abs(float(r_.x_cod)) for r_ in rows), max(abs(float(r_.y_cod)) for r_ in rows)
            r = D.direct_adts(lat, S['p'], 2e-8, order=2, rows=rows)
            t = dres(r)
            P('%s rms %.1e: |R|max=%.2e ety max=%.2e m COD max (%.1e, %.1e) m | 追跡 %s | 解析 %s | 相対差 %s | 解析の警告 %d 件' % (
                kind, amp, rmax, ey, cod[0], cod[1], f3(t), f3(A.total), fmt_rel(relv(t, A.total)), len(A.warnings)))
        except Exception as ex:
            P('%s rms %.1e: 失敗 %r' % (kind, amp, ex))
    save('coupling')


if __name__ == '__main__':
    sec = sys.argv[1]
    key = sys.argv[2] if len(sys.argv) > 2 else None
    if sec == 'octupole':
        octupole(key or 'kekls')
    elif sec == 'ablate':
        ablate(key)
    elif sec == 'delta':
        delta_sec(key, float(sys.argv[3]) if len(sys.argv) > 3 else 2e-3, float(sys.argv[4]) if len(sys.argv) > 4 else M.RINGS[key][4])
    elif sec == 'resonance':
        resonance()
    elif sec == 'coupling':
        coupling()
