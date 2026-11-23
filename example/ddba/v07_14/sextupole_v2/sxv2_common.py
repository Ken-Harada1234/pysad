# -*- coding: utf-8 -*-
"""sextupole_v2 の共通部分: v07_14 のラティスと 6 極の組 (../sextupole/sextupole_sets.json の final = ラティス埋込みの完成版)、K2L のランダム誤差つきの追跡の設定、
周波数マップの共鳴線の判定。どこからでも実行できる。環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages) を指定する。

試行錯誤は人手 (このフォルダの round*.py と README.md の経過)。ここの関数は、データを出す・見る・計算する道具だけで、反復も探索も自動化していない。
"""
import json
import math
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
OLD = os.path.join(HERE, '..', 'sextupole')               # 既存の例 (読むだけ。書き換えない)
sys.path.insert(0, ROOT)
sys.path.insert(0, OLD)
OUT = os.path.join(HERE, 'results')
os.makedirs(OUT, exist_ok=True)
P0 = 5e9                  # 運動量 [eV]
KMAX = 32.65              # |K2L| の上限 [1/m²] = 完成版の最大 |K2L|
warnings.simplefilter('ignore', RuntimeWarning)

import sxt_common as C                      # noqa: E402  (既存の例の共通部分: ラティスの読み込み・追跡の設定)

C.setup_pysad()
from pysad import sxtopt as S               # noqa: E402


def load():
    """(ソース, Lattice, SxtModel, 組の辞書 {名前: K2L の配列}, xi_target (完成版の解析式の ξ))"""
    src, lat = C.load_lattice('RING')
    m = S.SxtModel(lat, P0)
    sets = C.load_sets()
    assert sets['family_names'] == m.names
    ks = {n: np.array(v['k']) for n, v in sets['sets'].items()}
    return src, lat, m, ks, np.array(m.xi(ks['final']))


def sext_error_vector(lat, rms, seed):
    """600 台の 6 極の K2L の相対誤差 (rms、正規分布) の配列 (素子の順)。"""
    n = sum(1 for e in lat.elements if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0)
    return rms * np.random.default_rng(seed).standard_normal(n)


def make_setup(src, m, k, eps=None):
    """K2L = k (ファミリー) の追跡の設定 (RF・放射なし、LOSSAMPL 50 mm)。eps があれば K2≠0 の 6 極 (素子の順) の K2 に (1 + eps) を掛ける。"""
    s = S.replace_k2_in_source(src, m.names, k)
    st = C.make_track_setup(s, rad=False, rfsw=False, alost=0.05)
    if eps is not None:
        j = 0
        for e in st.els:
            if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0:
                e.K2 = float(e.K2) * (1.0 + float(eps[j]))
                j += 1
        assert j == len(eps)
    return st


def fmap_run(st, nx=16, ny=10, xmax=7e-3, ymax=5e-3, nturn=1024, loc=1):
    from pysad.freqmap import frequency_map
    x = np.linspace(xmax / nx, xmax, nx)
    y = np.linspace(ymax / ny, ymax, ny)
    return frequency_map(st, loc, x=x, y=y, nturn=nturn, workers=1, warn=False)


def lines_of_order(order):
    """|m|+|n| ≤ order の共鳴 (m, n) (m ≥ 0、m = 0 なら n > 0)。"""
    from pysad.freqmap import _mn_list
    return _mn_list(order)


def line_scores(res, qint=(47, 16), order=6, width=0.004, dthr=-4.0, min_pts=3):
    """周波数マップの各点 (生き残った点とその端数のチューン) を、共鳴線 m νx + n νy = p のそばの点 (チューン面の距離 < width) に分け、
    線のそばの点の中の「D が dthr より大きい (拡散が大きい) か、失われた点」の割合 (bad) と、そのそばの点の数を出す。
    戻り [(bad 割合、点の数、m, n, p, 次数)] (bad 割合の大きい順。点が min_pts 未満の線は除く)。p は整数部を含む共鳴 mνx+nνy=p。
    失われた点は、前半の周回で失われたものを含め、チューンが取れた点 (nux1・nuy1 が有限) だけ線に割り当てる。"""
    nux = np.asarray(res.nux1) + qint[0]
    nuy = np.asarray(res.nuy1) + qint[1]
    lost = np.asarray(res.lost_turn) > 0
    D = np.asarray(res.D)
    fin = np.isfinite(nux) & np.isfinite(nuy)
    bad = lost | ~np.isfinite(D) | (D > dthr)
    out = []
    for (m, n) in lines_of_order(order):
        v = m * nux + n * nuy
        p = np.round(v)
        dist = np.abs(v - p) / math.hypot(m, n)
        for pp in np.unique(p[fin & (dist < width)]):
            sel = fin & (dist < width) & (p == pp)
            if sel.sum() < min_pts:
                continue
            out.append((float(bad[sel].mean()), int(sel.sum()), m, n, int(pp), abs(m) + abs(n)))
    out.sort(key=lambda t: (-t[0], -t[1], t[5]))
    return out


def fmap_summary(res):
    D = np.asarray(res.D)
    lost = np.asarray(res.lost_turn) > 0
    ok = np.isfinite(D)
    bad = lost | ~ok | (D > -4.0)
    return dict(n=len(D), frac_lost=float(lost.mean()), frac_bad=float(bad.mean()), frac_chaotic_m3=float((lost | ~ok | (D > -3.0)).mean()),
                meanD=float(np.nanmean(D[ok])) if ok.any() else float('nan'))


def save_fmap(path, res):
    np.savez(path, x=res.x, y=res.y, nux1=res.nux1, nuy1=res.nuy1, nux2=res.nux2, nuy2=res.nuy2, D=res.D, lost=res.lost_turn)


class Saved:
    """save_fmap で保存した周波数マップを line_scores などに渡す形で読む。"""

    def __init__(self, path):
        d = np.load(path)
        self.x, self.y, self.nux1, self.nuy1, self.nux2, self.nuy2, self.D, self.lost_turn = (d[k] for k in ('x', 'y', 'nux1', 'nuy1', 'nux2', 'nuy2', 'D', 'lost'))


def lower_priority():
    C.lower_priority()


def error_rdt_rms(m, k, rms):
    """K2L のランダム誤差 (台ごとに独立、相対 rms) が作る 1 次の項 (CUBIC5・CHROM3) の rms の期待値 (統計平均: 台ごとの寄与の二乗和の平方根)。
    δh_t = Σ_i ε_i k_f(i) a_i,t (a_i,t = 1 台あたり・K2L 1 あたりの寄与、m.r1.per_element から) なので rms(δh_t)² = rms² Σ_f k_f² Σ_{i∈f} |a_i,t|²。
    位相の打ち消しがない (台ごとに独立) ので、ファミリーの強さの組み合わせでは打ち消せず、k_f の二乗の重みで効く。戻り {項: rms}"""
    nm_idx = {}
    for ir, nm in enumerate(m.r1.elem_names):
        nm_idx.setdefault(nm, []).append(ir)
    out = {}
    for t in S.CUBIC5 + S.CHROM3:
        it = m.terms1.index(t)
        s2 = 0.0
        for f, n in enumerate(m.names):
            a = m.r1.per_element[nm_idx[n], it] / m.k_ref[f]
            s2 += k[f] ** 2 * float(np.sum(np.abs(a) ** 2))
        out[t] = rms * math.sqrt(s2)
    return out


def add_cand(name, k, desc, spec=None, path=None):
    """候補を results/cands.json に足す (同名があれば ValueError。上書きしない。版は名前で加算: r1a、r1b …)。spec は記録のために文字列で残す。"""
    path = path or os.path.join(OUT, 'cands.json')
    c = json.load(open(path, encoding='utf-8'))
    if name in c:
        raise ValueError('候補 %s は既にある' % name)
    c[name] = dict(k=[float(x) for x in k], desc=desc, spec=spec)
    json.dump(c, open(path, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)


def show_k(m, tag, k, xt=None):
    """候補の 1 行の要約 (ADTS、1 次の RDT、4 次の主な項、max|K2L|、ξ)。"""
    a = m.adts(k)
    h = m.h1(k)
    h4 = m.h4nf(k)
    ii = S.Q._IDX4
    return ('%-10s ADTS %8.0f %8.0f %8.0f | h21 %.3g h30 %.3g h1011 %.3g h1002 %.3g h1020 %.3g | h20200 %.3g h20020 %.3g | max|k| %.1f xi (%.3f,%.3f)'
            % (tag, a['xx'], a['xy'], a['yy'], abs(h['h21000']), abs(h['h30000']), abs(h['h10110']), abs(h['h10020']), abs(h['h10200']),
               abs(h4[ii[(2, 0, 2, 0)]]), abs(h4[ii[(2, 0, 0, 2)]]), np.max(np.abs(k)), *m.xi(k)))


def lattice_with_errors(lat, m, k, eps):
    """K2L = k (ファミリー) で、K2≠0 の 6 極 (素子の順) の K2 に (1 + eps) を掛けた Lattice の複製 (追跡と同じ誤差の実現を解析式に渡す)。"""
    out = S.with_k2l(lat, m.names, k)
    j = 0
    for e in out.elements:
        if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0:
            e.K2 = float(e.K2) * (1.0 + float(eps[j]))
            j += 1
    assert j == len(eps)
    return out


def analytic_with_errors(lat, m, k, eps=None):
    """誤差つき (または誤差なし) の 6 極の解析式の値 (resdrive・rdt2nd を直接呼ぶ。多項式モデルは使わない: 台ごとの誤差はファミリーの変数にならないため)。
    戻り dict(adts={'xx','xy','yy'}, h1={項: 複素}, h4={名前: 複素}、xi=(ξx, ξy) 解析式)。数秒かかる。"""
    from pysad import resdrive as R
    from pysad import rdt2nd as Q
    lat2 = S.with_k2l(lat, m.names, k) if eps is None else lattice_with_errors(lat, m, k, eps)
    a = R.adts(lat2, momentum_eV=P0, fourth='r3', optics='calc').total
    r1 = R._resonance_terms_cubic(lat2, S.CUBIC5 + S.CHROM3 + S.DELTA_XI, bend_cubic=True, momentum_eV=P0)
    r4 = Q.rdt4(lat2, order2=True, momentum_eV=P0)
    h4 = {S._name4(q): complex(r4.total_arr[i]) for i, q in enumerate(Q.MONOS4)}
    xi_nat = m.xi_nat
    h4nf = np.array(r4.total_arr) + Q.nf_correction_h4(r4.h3, m.qx, m.qy)
    return dict(adts={c: float(a[c]) for c in ('xx', 'xy', 'yy')}, h1={t: complex(r1.total[t]) for t in S.CUBIC5 + S.CHROM3}, h4nf=h4nf,
                h4=h4, xi=(xi_nat[0] + r1.total['h11001'].real / math.pi, xi_nat[1] + r1.total['h00111'].real / math.pi))


def merit_parts(m, r, J_ref=1e-6, delta_ref=0.01):
    """analytic_with_errors の結果 r を、従来の目的関数 (SxtObjective: 物理的な重み。J_ref・δ_ref) の項 w3・chrom・w4・adts に直す (誤差つきの値の目的関数)。"""
    o = S.SxtObjective(m, J_ref=J_ref, delta_ref=delta_ref)
    h1 = np.array([r['h1'][t] for t in S.CUBIC5])
    hc = np.array([r['h1'][t] for t in S.CHROM3])
    a = np.array([r['adts'][c] for c in ('xx', 'xy', 'yy')])
    return dict(w3=float(np.sum(np.abs(h1 * o.s3) ** 2)), chrom=float(np.sum(np.abs(hc * o.sc) ** 2)), w4=float(np.sum(np.abs(r['h4nf'][o.i4] * o.s4) ** 2)),
                adts=float(np.sum((a * o.sA) ** 2)))
