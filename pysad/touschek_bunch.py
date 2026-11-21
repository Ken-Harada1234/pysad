# -*- coding: utf-8 -*-
"""v2.4.0 A3: Touschek の損失率に、バンチ内の δ の分布 (対の重心の δ_c) を入れる。オプション (既定は入れない)。

【原理】 ビームの静止系で 2 粒子の相対運動量はほぼ横方向 (縦は σδ/γ で小さい)。散乱後の 2 粒子の実験室系の δ は
  δ1' = δ_c + Δ、  δ2' = δ_c − Δ,   δ_c = (δ1 + δ2)/2,
Δ = γ q'_s/p0 は相対運動量の散乱 (重心系の角分布) だけで決まる。損失の条件は δ1' > a+ または δ2' < −a−: Δ > a+ − δ_c または Δ > a− + δ_c。
touschek_piwinski (Piwinski 型。K64・elegant も同じ) は δ_c = 0 (kick の基準 δ_old = 0)。ガウス分布のバンチでは δ_c ~ N(0, σδ²/2) で、相対運動量 (δ1 − δ2) と独立
(2 つのガウス変数の和と差は独立)。したがって 1 スライスの損失率は
  R_bunch = E_{δc}[ R_+(a+ − δ_c) + R_−(a− + δ_c) ]   (R_±: 片側の受容の J(2β a) の項)
で、a ≫ σδ では R(a − δ_c) ≈ R(a)(1 + n(n+1)/4 σδ²/a²) (R ∝ a^−n、n ≈ 2〜3)。δ_c の 1 次は ± で消える。a − δ_c の下限は floor × a (既定 0.25)。
【検証】 `example/v2_4/A3_mc_bunch.py`: 実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学、T2_mc と同じ。同じ標本で判定だけ変える)。
 δ_c の効果だけを取り出した比 (絶対値の判定 / δ_c を基準にした kick の判定 − 1) が、理論 (この式) と 2% 以内で一致する:
 ASTRID2 (a/σδ = 9.1・13.7) MC +2.84・+1.32% (±0.3%) に理論 +2.70・+1.28%、KEK-LS (a/σδ = 8.0・13.3) MC +4.38・+1.75% に理論 +4.32・+1.78%。
【入れていない効果 (限界)】
 (1) **相対運動量の縦成分 q_s の効果** (v2.4.0 A10 で `qs_ratio`・`touschek_piwinski_bunch(qs=True)` に入れた。既定は入れない): 通常の kick の判定 (δ_new − 元の自分の δ。K64・elegant・touschek_piwinski が同じ) は、
  δ_c を基準にした判定より損失率が小さい (MC: ASTRID2 で −2.0%、KEK-LS で −0.5%。相対運動量の δ 成分を 0 にすると 0)。この分、絶対値の判定の損失率は δ_c の効果より大きく増える
  (MC の全体の増え: ASTRID2 a/σδ = 9.1 で +4.8%、13.7 で +2.2%、KEK-LS 8.0 で +4.5%、13.3 で +1.8%。δ_c の効果だけの理論は ASTRID2 で +2.6%・+1.3%、KEK-LS で +4.0%・+1.6%)。
  q_s の効果は位置を固定したときの相対 δ の分散 (分散関数との相関) に依る。`qs_ratio` は 3 次元の積分 (T2_check3d の g = n'_s − û_s を、絶対値の判定の g = n'_s に変え、δ_c の平均をとる) で、
  絶対値 / kick の損失率の比を出す: ASTRID2 で +4.81・+2.22% (MC +4.84・+2.22%、スライス 3 つ・a 2 つで 0.1% 以内)。1 スライスに数秒かかるので、`touschek_piwinski_bunch(qs=True)` は
  12 スライスで求めて補正の係数 (3D の比 / δ_c だけの 1 次元の比) をスライスの順番で補間する。
 (2) a − δ_c が 0 に近い粒子は、小角散乱 (Rutherford の 1/K² の発散) で拡散と区別できない (IBS・量子励起の境界の問題)。a − δ_c の下限 floor × a で打ち切る。a ≳ 5σδ では打ち切りへの依存は
  無視できる (δ_c が 0.75a を超える確率 < 1e-8。floor 0.1・0.5 で結果が同じ: A3_mc_bunch の出力)。a ≲ 3σδ では使わない (警告)。
 (3) δ の分布はガウス (シンクロトロン振動の位相の一様分布 + 放射減衰・量子励起の定常分布)。σδ はスライスの実験室系の共分散 (IBS・バンチ長の伸びを含む) から取る。
 (4) 散乱後の粒子の横方向 (ベータトロンの) 振幅の跳び (η Δ) による横方向の受容 (nx, ny) の効果は入れない (LMA の定義と同じ)。受容は位置ごとの LMA δ±(s) (シンクロトロン位相に依らない値)。
 (5) elegant の touschekLifetime には δ の分布を入れる機能が無い (Piwinski の式のまま) ので、elegant との比較は既定の δ_c = 0 の値で行う (T1・T2)。
【T2 の概算との関係】 T2_limits.py の C は元の δ の分散を σδ² としていた (δ_c の σδ²/2 と q_s の効果を一緒にした上限に近い値。ASTRID2 では MC の全体と同程度、KEK-LS では 2 倍近く過大)。
"""
import math
import numpy as np

from .touschek import _piw_J, _R0_M, _CLIGHT_MS


def _slice_parts(sl):
    p0, h0 = sl['p0'], sl['h0']
    beta = p0 / h0
    vol = (4.0 * math.pi) ** 1.5 * math.sqrt(abs(np.linalg.det(sl['xx'])))
    P = np.asarray(sl['pl'])[:2, :2] * p0 ** 2
    ev = np.sort(np.linalg.eigvalsh(0.5 * (P + P.T)))
    l1 = max(ev[1], 1e-300)
    l2 = max(ev[0], 1e-12 * l1)
    return beta, vol, l1, l2, 1.0 / (4.0 * l1), 1.0 / (4.0 * l2)


def slice_sigma_delta(sl):
    """スライスの実験室系の σδ (lab 共分散の δδ 成分)"""
    lab = sl.get('lab')
    if lab is None:
        raise ValueError('スライスに lab 共分散が無い (touschek_table(record_rest=True) の結果が要る)')
    return math.sqrt(abs(float(np.asarray(lab)[5, 5])))


def piwinski_slice_rate_bunch(sl, ap, am, circ, pbunch, sigd=None, floor=0.25, nq=13, npan=40, var_factor=0.5):
    """touschek.piwinski_slice_rate にバンチ内の対の重心の δ_c ~ N(0, var_factor σδ²) (既定 var_factor = 0.5) の平均を入れた 1 スライスの損失率 [1/s]。
    sigd: σδ (省略するとスライスの lab 共分散から)、floor: a ∓ δ_old の下限 = floor × a、nq: Gauss-Hermite の点数、var_factor: δ_c の分散の倍率 (既定 0.5 = (δ1+δ2)/2 の分散。1 は診断用)。σδ = 0 なら piwinski_slice_rate と同じ"""
    p0, h0 = sl['p0'], sl['h0']
    beta, vol, l1, l2, a1, a2 = _slice_parts(sl)
    sd = slice_sigma_delta(sl) if sigd is None else float(sigd)
    pref = _R0_M ** 2 * _CLIGHT_MS * pbunch / (h0 * vol * math.sqrt(l1 * l2))
    if sd <= 0.0:
        t, w = np.array([0.0]), np.array([1.0])
    else:
        t, w = np.polynomial.hermite_e.hermegauss(nq)
        w = w / w.sum()
    dc = sd * math.sqrt(var_factor) * t
    J = 0.0
    for dci, wi in zip(dc, w):
        if ap > 0.0:
            J += wi * _piw_J(2.0 * beta * max(ap - dci, floor * ap), a1, a2, npan)
        if am > 0.0:
            J += wi * _piw_J(2.0 * beta * max(am + dci, floor * am), a1, a2, npan)
    return sl['al'] / circ * pref * J


def touschek_piwinski_bunch(result, a=None, ap=None, am=None, floor=0.25, nq=13, npan=40, sigd=None, var_factor=0.5, qs=False, qs_nsample=12, qs_nu=2048):
    """touschek.touschek_piwinski のバンチ内の δ 分布つき版 (戻り値は同じ形 + 'increase' = (損失率 / δ_c=0 の損失率) − 1、'rate0' = δ_c=0 の損失率、'amin_over_sigma' = min(a/σδ))。
    警告: a/σδ < 3 のスライスがあれば warnings に書く"""
    import warnings
    if 'rest' not in result:
        raise ValueError('touschek_piwinski_bunch: touschek_table(record_rest=True) の結果が必要')
    info = result['rest_info']
    if a is not None:
        ap = am = a
    nlat = info['nlat']

    def _arr(v):
        v = np.asarray(v, dtype=float)
        return np.full(nlat, float(v)) if v.ndim == 0 else v
    apa, ama = _arr(ap), _arr(am)
    from .touschek import piwinski_slice_rate
    rates, rates0, ratio = [], [], []
    for sl in result['rest']:
        l = sl['l'] - 1
        x, y = float(apa[l]), float(abs(ama[l]))
        sd = slice_sigma_delta(sl) if sigd is None else float(sigd)
        rates.append(piwinski_slice_rate_bunch(sl, x, y, info['circ'], info['pbunch'], sigd=sd, floor=floor, nq=nq, npan=npan, var_factor=var_factor))
        rates0.append(piwinski_slice_rate(sl, x, y, info['circ'], info['pbunch'], npan))
        if sd > 0:
            ratio.append(min(x, y) / sd)
    rates, rates0 = np.array(rates), np.array(rates0)
    qs_info = None
    if qs:
        # v2.4.0 A10: q_s の効果。3 次元の積分 (qs_ratio) の絶対値/kick の比を、qs_nsample 個のスライス (損失率の大きい順でなく、順番に等間隔) で求め、
        # δ_c だけの 1 次元の比で割った補正の係数を、スライスの順番で補間してかける (スライスごとに数秒かかるため)
        idx = np.nonzero(rates0 > 0)[0]
        pick = idx[np.unique(np.linspace(0, len(idx) - 1, min(qs_nsample, len(idx))).round().astype(int))]
        fs = []
        for j in pick:
            sl = result['rest'][j]
            l = sl['l'] - 1
            x, y = float(apa[l]), float(abs(ama[l]))
            r3, e3 = qs_ratio(sl, x, y, sigd=(slice_sigma_delta(sl) if sigd is None else float(sigd)), nu=qs_nu, floor=floor)
            fs.append(r3 / (rates[j] / rates0[j]))
        fs = np.array(fs)
        fk = np.interp(np.arange(len(rates)), pick, fs)
        rates = rates * fk
        qs_info = dict(slices=pick.tolist(), factor=fs.tolist())
    rate, rate0 = float(rates.sum()), float(rates0.sum())
    wr = []
    amin = min(ratio) if ratio else math.inf
    if amin < 3.0:
        wr.append('受容が σδ の %.1f 倍 (< 3): バンチ内の δ の分布の効果は展開 (a ≫ σδ) の範囲外。打ち切り (floor) に依存する' % amin)
    return dict(rate=rate, lifetime=(1.0 / rate if rate > 0 else math.inf), rate_per_slice=rates, rate0=rate0, increase=rate / rate0 - 1.0, qs=qs_info,
                l=np.array([sl['l'] for sl in result['rest']]), amin_over_sigma=amin, warnings=wr)


# ---------------------------------------------------------------------------------------------------------------------
# v2.4.0 A10: 相対運動量の縦成分 q_s の効果 (3 次元の積分)。損失の判定を kick (δ_new − 元の自分の δ) から絶対値 (δ_new) に変えたときの損失率の比
# ---------------------------------------------------------------------------------------------------------------------
_R0_M_ = 2.8179403262e-15
_gx, _gw = np.polynomial.legendre.leggauss(48)
_rx, _rw = np.polynomial.legendre.leggauss(40)


def _Kth(th):
    s2 = np.sin(th) ** 2
    return 4 / s2 ** 2 - 3 / s2


def _what(t, ca, absmode):
    """半球 (θ < π/2、û のまわり) の ∫ K(θ) [I(g > t) + I(g < −t)] dΩ'。g = n'_s − û_s (kick の判定、absmode=False)、または g = n'_s (絶対値の判定、True)。
    ca = û_s、t = 2β a/k (配列)。T2_check3d.What の一般化"""
    t = np.asarray(t, float)
    ca = np.asarray(ca, float)
    sa = np.sqrt(np.maximum(1 - ca ** 2, 0))
    th_lo = (0.2 * np.minimum(t, 2.0)) if absmode else 2 * np.arcsin(np.minimum(t / 2, 1.0))
    th_hi = np.full_like(t, math.pi / 2)
    ok = th_lo < th_hi
    lo = np.where(ok, np.log(np.maximum(th_lo, 1e-12)), 0.0)
    hi = np.log(th_hi)
    u = 0.5 * (hi - lo)[:, None] * _gx[None, :] + 0.5 * (hi + lo)[:, None]
    th = np.exp(u)
    wt = 0.5 * (hi - lo)[:, None] * _gw[None, :] * th
    st, ct = np.sin(th), np.cos(th)
    cs, ss, tt = ca[:, None], sa[:, None], t[:, None]
    den = ss * st
    X = (-cs * ct) if absmode else (cs * (1 - ct))
    with np.errstate(divide='ignore', invalid='ignore'):
        y1 = (tt + X) / den
        y2 = (-tt + X) / den
    m1 = np.where(y1 >= 1, 0.0, np.where(y1 <= -1, 2 * np.pi, 2 * np.arccos(np.clip(y1, -1, 1))))
    m2 = np.where(y2 >= 1, 2 * np.pi, np.where(y2 <= -1, 0.0, 2 * (np.pi - np.arccos(np.clip(y2, -1, 1)))))
    small = den < 1e-14
    if small.any():
        g = (cs * ct) if absmode else (-cs * (1 - ct))
        m1 = np.where(small, np.where(g >= tt, 2 * np.pi, 0.0), m1)
        m2 = np.where(small, np.where(g <= -tt, 2 * np.pi, 0.0), m2)
    f = st * _Kth(th) * (m1 + m2) * wt
    f = np.where(ok[:, None], f, 0.0)
    return f.sum(axis=1)


def qs_ratio(sl, ap, am, sigd=None, nu=2048, ngh=7, floor=0.25, seed=1, nk=40):
    """1 スライスの、絶対値の判定 (δ_c ~ N(0, σδ²/2) の平均、相対運動量の縦成分 q_s を含む) の損失率 / kick の判定の損失率 (3 次元の積分、同じ標本)。
    1 次元の式 (touschek_piwinski) は kick の判定の q_s = 0 の値で、3 次元の kick の判定と 1e-3 で一致する (T2) ので、この比を 1 次元の損失率にかけると、
    δ_c と q_s の両方の効果が入る (piwinski_slice_rate_bunch は δ_c だけ)。戻り (比, 比の統計誤差 (標本の相関を使わない上限))。
    sl['pl'] (位置を固定したときの運動量の共分散、3×3、p0 で規格化)・sl['p0']・sl['h0'] を使う"""
    from scipy.stats import qmc
    p0, h0 = sl['p0'], sl['h0']
    beta = p0 / h0
    S = 2 * np.asarray(sl['pl'])[:3, :3] * p0 ** 2
    sd = slice_sigma_delta(sl) if sigd is None else float(sigd)
    Lc = np.linalg.cholesky(S)
    z = qmc.MultivariateNormalQMC(mean=np.zeros(3), cov=np.eye(3), seed=seed).random(nu)
    q = z @ Lc.T
    kk = np.linalg.norm(q, axis=1)
    u = q / kk[:, None]
    A = 0.5 * np.einsum('ni,ij,nj->n', u, np.linalg.inv(S), u)
    ca = u[:, 2]
    tg, wg = np.polynomial.hermite_e.hermegauss(ngh)
    wg = wg / wg.sum()
    dcs = sd / math.sqrt(2.0) * tg

    def radial(c, absmode):
        kmax = np.sqrt(45.0 / A)
        klo = np.full(nu, c / 1.5)
        okk = kmax > klo
        lo = np.log(klo)
        hi = np.log(np.maximum(kmax, klo * 1.0001))
        lk = 0.5 * (hi - lo)[:, None] * _rx[None, :] + 0.5 * (hi + lo)[:, None]
        k = np.exp(lk)
        wk = 0.5 * (hi - lo)[:, None] * _rw[None, :]
        Wv = _what((c / k).ravel(), np.repeat(ca, nk), absmode).reshape(nu, nk)
        return np.where(okk, (np.exp(-k ** 2 * A[:, None]) * Wv * wk).sum(axis=1), 0.0)
    pref = 4 * A ** 1.5 / math.sqrt(math.pi)
    kick = (radial(2 * beta * ap, False) + radial(2 * beta * am, False)) * pref
    ab = np.zeros(nu)
    for dc, wi in zip(dcs, wg):
        ab += wi * (radial(2 * beta * max(ap - dc, floor * ap), True) + radial(2 * beta * max(am + dc, floor * am), True)) * pref
    r = ab.mean() / kick.mean()
    return float(r), float(np.std(ab / kick.mean() - r * kick / kick.mean()) / math.sqrt(nu))
