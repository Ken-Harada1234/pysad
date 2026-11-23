# -*- coding: utf-8 -*-
"""共鳴駆動項 (RDT) h_jklmp の解析式 (resdrive.resonance_terms) を、追跡から取り出した値と比べる。

使い方 (リポジトリ直下から)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_tracking.py table  fodo|fodo_rect|fodo_asym|astrid|astrid_mpw|kekls   (a) 1 周の写像の 2 次のテイラー展開 (p=0 の 3 次の項、δ の項)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_tracking.py chrom  fodo|astrid|astrid_mpw|kekls                       δ の 1 次の項と自然クロマティシティ: 解析 (自然 + 6 極 + 偏向磁石) と CALC・追跡
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_tracking.py res    fodo|astrid|astrid_mpw|kekls                       (b) 3νx = 整数の近くの 3 周写像の周期点: |h30000| と位相の検証
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_tracking.py thin                                                      薄い 6 極の δ の項の閉じた式との照合 (導出は docs/25 §9)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_tracking.py plot                                                      rdt_table.png (table の出力から)
条件 (すべて): pysad k64 の追跡 (RF 空洞止め `rfsw=False`、放射なし、結合なし、閉軌道のまわり)、ring_map の振幅は正規化座標で 1e-5 √m (J = 5e-11 m·rad)、
  リング先頭を位相の基準 (μ = 0)。6 極の強さは K2 を λ 倍 (λ = 1 と 0 の差が 6 極の 1 次、λ = 0 が偏向磁石など 6 極以外)。δ の項は中心差分 δ = ±1e-4、±2e-4 の Richardson。
  エネルギー: FODO 3 GeV、ASTRID2 0.58 GeV、大きなリング A (kek-ls_v03_12.lat) 3 GeV。
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
from pysad import resdrive as R, rdttrack as T, compat
compat.set_mode(os.environ.get('PYSAD_COMPAT', 'k64'))
import adts_bend as M

CUBIC5 = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
DELTA5 = ('h11001', 'h00111', 'h20001', 'h00201', 'h10002')
OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(s)


def save(name, ring, data=None):
    base = os.path.join(HERE, 'rdt_%s_%s.pysad' % (name, ring))
    with open(base + '.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')
    if data is not None:
        with open(base + '.json', 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else str(o))


def cx(v):
    return '%+.6f%+.6fj' % (v.real, v.imag)


def relerr(a, t):
    return abs(a - t) / max(abs(t), 1e-300)


def cmd_table(key, amp=1e-5):
    S1, S0 = M.setup(key, 1.0), M.setup(key, 0.0)
    el1, el0, p = S1['lat'].elements, S0['lat'].elements, S1['p']
    row0 = S1['rows'][0]
    nu = S1['nu']
    P('R4 (a) 1 周の写像の 2 次のテイラー展開から取り出した RDT 対 解析 resonance_terms: %s (%.2f GeV、周長 %.1f m、νx=%.6f νy=%.6f)' % (key, p / 1e9, S1['circ'], nu[0], nu[1]))
    P('条件: pysad k64 の追跡、RF・放射なし、閉軌道のまわり、ring_map の振幅 %.0e √m (正規化座標、Richardson)、リング先頭が位相の基準、δ は ±1e-4・±2e-4 の Richardson' % amp)
    t0 = time.time()
    tr1 = T.tracked_rdt_table(el1, p, row0, amp=amp)
    tr0 = T.tracked_rdt_table(el0, p, row0, amp=amp)
    P('追跡からの取り出し %.1f s (6 極 ×1 と ×0)' % (time.time() - t0))
    cub = list(T.CUBIC_MONOS)
    dl = list(T.DELTA1_MONOS) + [(1, 0, 0, 0, 2)]
    t0 = time.time()
    rb = R.resonance_terms(S1['lat'], cub + dl, rows=S1['rows'], bend_cubic=True)
    ra = R.resonance_terms(S1['lat'], cub + dl, rows=S1['rows'], bend_cubic=False)
    P('解析 %.1f s (偏向磁石の 3 次の項は追跡の写像から取り出し、素子の種類ごとに 1 回)。警告 %d 件' % (time.time() - t0, len(rb.warnings)))
    for w in rb.warnings[:5]:
        P('  (注意) ' + w)
    data = {'ring': key, 'amp': amp, 'rows': []}
    P('')
    P('[p=0 の 3 次の項] h = 全周の和 (H の係数、(2J)^((j+k)/2)(2J)^((l+m)/2) の前の係数、複素数)')
    P('項      追跡 (6 極 ×1)                    解析 R3 (6 極のみ)                   解析 R4 (6 極 + 偏向磁石)         差 R3     差 R4     偏向磁石 追跡(λ=0)  解析     (差)/|h|')
    zero = 0.0
    for q in T.CUBIC_MONOS:
        nm = T.mono_name(q)
        if q[0] < q[1] or (q[0] == q[1] and q[2] < q[3]):
            continue
        t1, t0_ = tr1['cubic'][nm], tr0['cubic'][nm]
        if nm not in CUBIC5:
            zero = max(zero, abs(t1), abs(ra.total[nm]), abs(rb.total[nm]))
            continue
        a, b = ra.total[nm], rb.total[nm]
        bb = rb.total_bend[nm]
        P('%-7s %s  %s  %s  %.2e  %.2e  %s %s  %.1e' % (nm, cx(t1), cx(a), cx(b), relerr(a, t1), relerr(b, t1), cx(t0_), cx(bb), abs(t0_ - bb) / abs(b)))
        data['rows'].append({'term': nm, 'track': t1, 'sext': a, 'total': b, 'track_bend': t0_, 'an_bend': bb})
    P('他の 3 次の項 (h20010・h11010・h02100 など、0 のはず) の最大 |h| = %.2e (追跡・解析とも)' % zero)
    P('')
    P('[6 極だけの分: 追跡の ×1 − ×0 と 解析 R3 (6 極のみ)] (6 極の強さに比例する 1 次の部分)')
    for nm in CUBIC5:
        d = tr1['cubic'][nm] - tr0['cubic'][nm]
        P('%-7s 追跡 %s  解析 %s  相対差 %.2e' % (nm, cx(d), cx(ra.total[nm]), relerr(ra.total[nm], d)))
        data['rows'][[r['term'] for r in data['rows']].index(nm)]['sext_rel'] = relerr(ra.total[nm], d)
    P('')
    P('[δ の項 (p ≥ 1)] 6 極の分: 追跡 (×1 − ×0、閉軌道のヤコビ行列の δ 微分) と 解析 (6 極のみ)')
    for nm in DELTA5:
        if nm == 'h10002':
            d = tr1['h10002'] - tr0['h10002']
        else:
            d = tr1['delta1'][nm] - tr0['delta1'][nm]
        a = ra.total[nm]
        P('%-7s 追跡 %s  解析 %s  相対差 %.2e' % (nm, cx(d), cx(a), relerr(a, d)))
        data['rows'].append({'term': nm, 'track_sext': d, 'sext': a, 'total': rb.total[nm], 'an_bend': rb.total_bend[nm], 'sext_rel': relerr(a, d)})
    P('偏向磁石の分 (解析、分散による軌道のずれの展開): ' + '、'.join('%s %s' % (nm, cx(rb.total_bend[nm])) for nm in DELTA5))
    P('  (追跡の ×0 は自然なクロマティックな項を含み、偏向磁石の分だけを取り出せない。h11001・h00111 の偏向磁石の分は `chrom` で CALC と照合)')
    save('table', key, data)


def beta_nux_scan(S, target_nux, qfscale=None):
    """K1>0 の四極 (QF 系) の強さを f 倍して νx = target にする。戻り (f, lat, rows)"""
    import copy
    from scipy.optimize import brentq
    lat0 = S['lat']
    idx = [i for i, e in enumerate(lat0.elements) if type(e).__name__ == 'Quadrupole' and float(e.K1) > 0]
    k0 = {i: float(lat0.elements[i].K1) for i in idx}

    def build(f):
        lat = copy.deepcopy(lat0)
        for i in idx:
            lat.elements[i].K1 = k0[i] * f
        return lat

    def nux(f):
        return float(build(f).twiss_table()[-1].mux) / (2 * math.pi)
    f0 = 1.0
    v0 = nux(f0) - target_nux
    step = 5e-4 if v0 < 0 else -5e-4
    br = None
    for _ in range(200):
        f1 = f0 + step
        try:
            v1 = nux(f1) - target_nux
        except Exception:
            break
        if v0 * v1 < 0:
            br = (min(f0, f1), max(f0, f1))
            break
        f0, v0 = f1, v1
    f = brentq(lambda f: nux(f) - target_nux, br[0], br[1], xtol=1e-13)
    lat = build(f)
    return f, lat, lat.twiss_table()


def cmd_res(key, eps_list=None):
    S = M.setup(key, 1.0)
    p = S['p']
    nux0 = S['nu'][0]
    n3 = round(3 * nux0)
    base = n3 / 3.0
    P('R4 (b) 共鳴線の近くの周期点: %s (%.2f GeV、νx=%.6f → 3νx = %d + 3ε に動かす。K1>0 の四極を f 倍)' % (key, p / 1e9, nux0, n3))
    P('1 周の生成関数 G = 2π[νJ + a J²/2] + 2|h|(2J)^{3/2} cos(3φ + θ)、θ = arg h30000。3 周写像の周期点 (y=0、リング先頭の φ: x = √(2Jβ) cos φ) の条件:')
    P('  2π(ε + a J) + 6|h| √(2J) cos(3φ + θ) = 0 (接線方向) と sin(3φ + θ) = 0 (動径方向、φ が 6 点に決まる)。a は共鳴に依らない ∂νx/∂Jx = 解析の合計 − 共鳴の特異部分 (−12|h30000|²/(π² ε))')
    P('  予測 J は cos = ±1 の 2 つ (J_lo、J_hi。周期点の J の最小・最大)。|h30000| の解析は 6 極 + 偏向磁石 (R4)、「6 極のみ」(R3) の予測も併記。')
    P('  注意: リングが k 回対称 (FODO は 24 セル) で 3νx = 整数が k の倍数でないとき、共鳴は駆動されず (h30000 ∝ ε)、3 周写像の不動点は 6 点でなく曲線になる (J(φ) が J_lo〜J_hi を動く)。J の最小・最大は同じ式で予測できる (R2 と同じ検証)。')
    eps_list = eps_list or {'fodo': (-2e-3, -1e-3, -5e-4), 'astrid': (-2e-3, -1e-3, -5e-4),
                            'astrid_mpw': (-2e-3, -1e-3, -5e-4), 'kekls': (-1e-3, -5e-4, -2e-4)}[key]
    data = {'ring': key, 'rows': []}
    for eps in eps_list:
        f, lat, rows = beta_nux_scan(S, base + eps)
        nux = float(rows[-1].mux) / (2 * math.pi)
        e_ = nux - base
        hs = R.resonance_terms(lat, ('h30000',), rows=rows, bend_cubic=False).total['h30000']
        hb = R.resonance_terms(lat, ('h30000',), rows=rows, bend_cubic=True).total['h30000']
        A = R.adts(lat, rows=rows)
        sing = -12.0 * abs(hb) ** 2 / (math.pi ** 2 * e_)         # 共鳴の項 (k=±(3,0)、n = 3νx に最も近い整数) の 2 次の寄与 (特異部分)
        a_nr = A.total['xx'] - sing
        th = math.atan2(hb.imag, hb.real)
        pred = T.predict_periodic_J(e_, a_nr, abs(hb))
        pred_s = T.predict_periodic_J(e_, a_nr, abs(hs))
        row0 = rows[0]
        pts = []
        for jg in (pred or [abs(e_) / abs(a_nr)]):
            pts += T.periodic_points(lat.elements, p, row0, 3, jg, nphase=24)
        uniq = []
        for J, ph in pts:
            if not any(abs(J - u[0]) < 1e-6 * J and abs(math.remainder(ph - u[1], 2 * math.pi)) < 1e-3 for u in uniq):
                uniq.append((J, ph))
        P('')
        P('ε=%+.3e (f=%.6f νx=%.7f): |h30000| 6 極のみ %.5f、+偏向磁石 %.5f (θ=%.4f rad)、a=%.2f (合計 %.2f、特異部分 %.2f) [1/(m rad)]' % (
            e_, f, nux, abs(hs), abs(hb), th, a_nr, A.total['xx'], sing))
        P('  予測 J (+偏向磁石): %s | (6 極のみ): %s' % (['%.5e' % v for v in pred], ['%.5e' % v for v in pred_s]))
        if not uniq:
            P('  追跡の周期点: なし')
            continue
        Js = [u[0] for u in uniq]
        for J, ph in uniq:
            data['rows'].append({'eps': e_, 'J': J, 'phi': ph, 'h_sext': abs(hs), 'h_total': abs(hb)})
        msg = '  追跡の周期点 %d 個: J の最小 %.5e 最大 %.5e' % (len(uniq), min(Js), max(Js))
        if len(pred) == 2:
            msg += ' | 予測との差 (最小 %+.2f%%、最大 %+.2f%%)、6 極のみの予測との差 (%+.2f%%、%+.2f%%)' % (
                100 * (min(Js) / pred[0] - 1), 100 * (max(Js) / pred[1] - 1), 100 * (min(Js) / pred_s[0] - 1) if len(pred_s) == 2 else float('nan'),
                100 * (max(Js) / pred_s[1] - 1) if len(pred_s) == 2 else float('nan'))
        P(msg)
        if len(uniq) >= 20:
            # k 回対称のリング: 不動点は曲線 (J(φ))。J(φ) の高調波フィット J = c0 + Σ_{k=1..6} (a_k cos kφ + b_k sin kφ) の 3φ 成分から |h30000| と θ を求める。
            # 正規形 (1 次のリー変換 W: (e^{i3μ} − 1) W = h、|W| = |h|/(6π|ε|)) では、不変曲線 I = 一定の J(φ) の 3φ 成分は A cos(3φ + θ)、A = 6|W|(2J0)^{3/2}/(2π...) = 2|h| u0/(π a)
            ph_ = np.array([u[1] for u in uniq])
            J_ = np.array([u[0] for u in uniq])
            cols = [np.ones_like(ph_)] + [f_(k * ph_) for k in range(1, 7) for f_ in (np.cos, np.sin)]
            coef = np.linalg.lstsq(np.array(cols).T, J_, rcond=None)[0]
            a3, b3 = coef[1 + 2 * 2], coef[2 + 2 * 2]
            A3 = math.hypot(a3, b3)
            u0 = math.sqrt(2 * coef[0])
            h_inf = math.pi * abs(a_nr) * A3 / (2 * u0)
            th_inf = math.atan2(-b3, a3)               # a3 cos 3φ + b3 sin 3φ = A cos(3φ + θ) → a3 = A cos θ、b3 = −A sin θ
            resid = float(np.abs(np.array(cols).T @ coef - J_).max() / coef[0])
            P('  J(φ) の 3φ 成分から: |h30000|_推定 %.5f (解析 +偏向磁石 %.5f: %+.2f%%、6 極のみ %.5f: %+.2f%%)、arg_推定 %.4f rad (解析 %.4f、差 %+.4f)、フィットの残り %.1e (J の相対)' % (
                h_inf, abs(hb), 100 * (h_inf / abs(hb) - 1), abs(hs), 100 * (h_inf / abs(hs) - 1), th_inf, th, math.remainder(th_inf - th, 2 * math.pi), resid))
            data['rows'][-1]['h_inferred'] = h_inf
            data['rows'][-1]['theta_inferred'] = th_inf
            data['rows'][-1]['theta_pred'] = th
            data['rows'][-1]['h_err_total'] = h_inf / abs(hb) - 1
            data['rows'][-1]['h_err_sext'] = h_inf / abs(hs) - 1
        if len(pred) == 2 and len(pred_s) == 2:
            sp_t, sp_b, sp_s = max(Js) - min(Js), pred[1] - pred[0], pred_s[1] - pred_s[0]
            P('  J の幅 (最大 − 最小): 追跡 %.4e、予測 (+偏向磁石) %.4e (追跡/予測 − 1 = %+.2f%%)、予測 (6 極のみ) %.4e (%+.2f%%)' % (sp_t, sp_b, 100 * (sp_t / sp_b - 1), sp_s, 100 * (sp_t / sp_s - 1)))
            data['rows'][-1]['spread_err'] = sp_t / sp_b - 1
            data['rows'][-1]['spread_err_sext'] = sp_t / sp_s - 1
        data['rows'][-1]['jmin_err'] = (min(Js) / pred[0] - 1) if len(pred) == 2 else None
        data['rows'][-1]['jmax_err'] = (max(Js) / pred[1] - 1) if len(pred) == 2 else None
        data['rows'][-1]['jmin_err_sext'] = (min(Js) / pred_s[0] - 1) if len(pred_s) == 2 else None
        data['rows'][-1]['jmax_err_sext'] = (max(Js) / pred_s[1] - 1) if len(pred_s) == 2 else None
    save('res', key, data)


def cmd_chrom(key):
    P('R4 δ の 1 次の項と自然クロマティシティ: %s' % key)
    P('CALC の δ 依存 (chromaticity_numeric、order=2) と、解析 = 自然 (chromaticity_natural: −(1/4π)∮(K1+h²)β + 端) + 偏向磁石 (chromaticity_bend: h11001/π の偏向磁石の分) + 6 極 (chromaticity_sextupole)')
    data = {'ring': key, 'rows': []}
    for lam in (0.0, 1.0):
        S = M.setup(key, lam)
        ch = R.chromaticity_numeric(S['lat'], order=2)['d1']
        n = R.chromaticity_natural(S['lat'])
        b = R.chromaticity_bend(S['lat'])
        s = R.chromaticity_sextupole(S['lat'])
        tot = [n[i] + b[i] + (s[i] if lam else 0.0) for i in (0, 1)]
        P('λ=%g: CALC ξ = (%.5f, %.5f)  解析 = (%.5f, %.5f) [自然 (%.5f, %.5f) + 偏向磁石 (%.5f, %.5f) + 6 極 (%.5f, %.5f)]  差 (%+.2e, %+.2e)' % (
            lam, ch[0], ch[1], tot[0], tot[1], n[0], n[1], b[0], b[1], (s[0] if lam else 0), (s[1] if lam else 0), ch[0] - tot[0], ch[1] - tot[1]))
        data['rows'].append({'lam': lam, 'calc': ch, 'natural': n, 'bend': b, 'sext': s if lam else (0, 0), 'diff': [ch[0] - tot[0], ch[1] - tot[1]]})
    P('λ=0 (6 極なし) の差は偏向磁石の 3 次の項を入れた解析と CALC の一致 (自然クロマティシティの式 + 偏向磁石の h11001/π)。λ=1 の差は 6 極の 2 次 (K2² の δ 依存) 以上')
    save('chrom', key, data)


def cmd_thin():
    """薄い 6 極 (FODO の SF・SD を slice_thin で薄レンズ列にしたもの) で、δ の項の閉じた式と resonance_terms (解析・数値積分) を比べる"""
    import math
    S = M.setup('fodo', 1.0)
    lat = R.slice_thin(S['lat'], 4)
    rows = lat.twiss_table()
    els = lat.elements
    sums = {k: 0j for k in DELTA5}
    for i, e in enumerate(els):
        if type(e).__name__ != 'Multipole' or not (getattr(e, 'KnL', {}) or {}).get(2):
            continue
        k2l = float(e.KnL[2])
        r = rows[i]
        bx, by, eta = float(r.betx), float(r.bety), float(r.etx)
        mx, my = float(r.mux), float(r.muy)
        sums['h11001'] += k2l * eta * bx / 4
        sums['h00111'] += -k2l * eta * by / 4
        sums['h20001'] += k2l * eta * bx * np.exp(2j * mx) / 8
        sums['h00201'] += -k2l * eta * by * np.exp(2j * my) / 8
        sums['h10002'] += k2l * eta * eta * math.sqrt(bx) * np.exp(1j * mx) / 4
    r1 = R.resonance_terms(lat, DELTA5, rows=rows, bend_cubic=False)
    r2 = R.resonance_terms(lat, DELTA5, rows=rows, method='numeric', nslice=4, bend_cubic=False)
    P('R4 薄い 6 極 (FODO の SF・SD を 4 分割した薄レンズ列、3 GeV) の δ の項: 閉じた式 (H = (K2L/6)(x³ − 3xy²)、x → x + ηδ から導出) と resonance_terms の差')
    P('閉じた式: h11001 = ¼ΣK2L η βx、h00111 = −¼ΣK2L η βy、h20001 = ⅛ΣK2L η βx e^{2iμx}、h00201 = −⅛ΣK2L η βy e^{2iμy}、h10002 = ¼ΣK2L η² √βx e^{iμx}')
    for k in DELTA5:
        P('%-7s 閉じた式 %s  resonance_terms(解析) %s  差 %.1e  (数値積分 nslice=4 との差 %.1e)' % (k, cx(sums[k]), cx(r1.total[k]), abs(sums[k] - r1.total[k]), abs(r2.total[k] - r1.total[k])))
    save('thin', 'fodo')


def cmd_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False
    rings = ['fodo', 'fodo_rect', 'fodo_asym', 'astrid', 'astrid_mpw', 'kekls']
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    w = 0.38
    for k, ring in enumerate(rings):
        fn = os.path.join(HERE, 'rdt_table_%s.pysad.json' % ring)
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn, encoding='utf-8'))
        rows = [r for r in d['rows'] if r['term'] in CUBIC5 and 'track' in r]
        for j, r in enumerate(rows):
            t = complex(*r['track'])
            s_, tot = complex(*r['sext']), complex(*r['total'])
            ax[0].plot(k + (j - 2) * 0.12 - 0.0, max(abs(s_ - t) / abs(t), 1e-10), 'o', color='tab:red', ms=5, label='6 極のみ (R3)' if (k == 0 and j == 0) else None)
            ax[0].plot(k + (j - 2) * 0.12 + 0.05, max(abs(tot - t) / abs(t), 1e-10), 's', color='tab:blue', ms=5, label='6 極 + 偏向磁石 (R4)' if (k == 0 and j == 0) else None)
    ax[0].set_yscale('log')
    ax[0].set_xticks(range(len(rings)))
    ax[0].set_xticklabels(['FODO', 'FODO\n長方形', 'FODO\n非対称', 'ASTRID2', 'ASTRID2\nMPW', '大きな\nリング A'], fontsize=8)
    ax[0].set_ylabel('|解析 − 追跡| / |追跡|  (5 項)')
    ax[0].set_title('p=0 の 3 次の共鳴駆動項: 1 周の写像の 2 次の展開との差')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    for ring, c, mk in (('fodo', 'tab:blue', 'o'), ('astrid', 'tab:green', 's'), ('kekls', 'tab:red', '^')):
        fn = os.path.join(HERE, 'rdt_res_%s.pysad.json' % ring)
        if not os.path.exists(fn):
            continue
        rows = [r for r in json.load(open(fn, encoding='utf-8'))['rows'] if 'h_err_total' in r]
        if not rows:
            continue
        ax[1].plot([r['eps'] for r in rows], [100 * r['h_err_total'] for r in rows], mk, color=c, label=ring + ' (6 極 + 偏向磁石)')
        ax[1].plot([r['eps'] for r in rows], [100 * r['h_err_sext'] for r in rows], mk, mfc='none', color=c, label=ring + ' (6 極のみ)')
    ax[1].axhline(0, color='k', lw=0.6)
    ax[1].set_xlabel('ε = νx − 整数/3')
    ax[1].set_ylabel('|h30000| (周期点の J(φ) の 3φ 成分から) / 解析 − 1 [%]')
    ax[1].set_title('3νx 共鳴線の近くの不動点から求めた |h30000|')
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'rdt_table.png'), dpi=130)
    print('saved rdt_table.png')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'table':
        cmd_table(sys.argv[2])
    elif cmd == 'res':
        cmd_res(sys.argv[2])
    elif cmd == 'chrom':
        cmd_chrom(sys.argv[2])
    elif cmd == 'thin':
        cmd_thin()
    elif cmd == 'plot':
        cmd_plot()
    else:
        print(__doc__)
