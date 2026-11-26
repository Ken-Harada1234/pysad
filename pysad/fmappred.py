# -*- coding: utf-8 -*-
"""Rebuilt resonance-line prediction for frequency maps using RDT and non-resonant corrections.
周波数マップの共鳴線の予想の作り直し (v2.3.0 S3、K64にない新機能。互換モードによらず同じ)。

S1 (`fmapres.ResonancePredictor`) は 1 次の RDT (3 次の h) の振り子の式、R5 は 4 次の h (写像の係数) も足した。S3 は R6 の道具で次を変える:
  (a) 駆動項: 3 次の h (`resonance_terms`) + 4 次の h に**非共鳴の 3 次の摂動の補正** (`rdt2nd.nf_correction_h4`、F4' = F4 + ½{R⁻¹W, W}) を足した値。
      補正はチューンに依る (1/(1 − e^{iφ_q}))。`nf_local=True` (既定は False) は補正を初期条件ごとの ν(J) で評価する診断用の選択肢 (作動点 ν0 が定義された値)。
  (b) 作動点からの距離: ν(J) = ν0 + A·J (A = `adts()` の既定 = R6)。J² の項は解析式が未実装 (docs/25 §10.6) で、追跡に合わせた当てはめの項は**足さない** (方針: 摂動論の定義された次数だけ)。
      追跡のチューンを `predict(nu=...)` に渡すのは、ν(J) の誤差の効果を測る診断で、予想ではない。
  (c) 幅: 振り子 H = 2π∫Δ0 dI + V cos(gψ + arg F)、V = 2|F|、G0 = ∂Δ0/∂I (docs/25 §12.3)。島の半幅 (チューン平面、Δ0 の単位) W = √(2|G0|V/π)、位相つきの捕捉条件は S1 と同じ。
      F は共鳴の調和 a ごとに、a の符号の代表の 1 つ (a と −a は共役) で作る。次数 (3, 4) の偶奇が違うので、同じ調和に 3 次と 4 次が同居することはない (3 次どうしの 2 次の効果 = F4' の Magnus の項)。
  (d) 領域: 初期条件 (x, y) の点ごとに、島の中 (位相つき `any`・包絡 `envelope`) と連続の**スコア** (log10(W/|Δ0|) の調和にわたる最大。0 以上なら包絡の中) を返す。
  (e) 信頼度のフラグ (点ごと `reliability` 0/1/2 と、予想全体の `warnings`): `reliability()`・docs/25 §12.5。

使い方:
  pr = ResonancePredictorS3.from_lattice(lat, loc, rows)      # または fmapres.ResonancePredictor.from_lattice(..., model='s3')
  m = pr.predict(r.x, r.y)                                    # {'any','envelope','score','n_env','reliability','by_harmonic'}
  evaluate_prediction(r, m, thr=-5.0)                         # cov・面積・対照・AUC・順位相関・上位 10% の D
"""
import math

import numpy as np

from . import fmapres as FR
from . import rdt2nd as _Q
from . import resdrive as _R

TWO_PI = 2.0 * math.pi
_MONOS3 = _Q.MONOS3
_MONOS4 = _Q.MONOS4


def _mon_exps(q):
    return (q[0] + q[1]) / 2.0, (q[2] + q[3]) / 2.0


class ResonancePredictorS3(FR.ResonancePredictor):
    """S3 resonance-line predictor holding the third- and fourth-order harmonic coefficients, linear ADTS matrix, tunes and Twiss.
    S3 の共鳴線の予想。属性: h3 (MONOS3 の順 20 複素)、h4 (MONOS4 の順 35 複素、写像の係数)、h3sum (正規形の補正に使う 3 次の 1 周の和)、A (2×2 線形の ADTS)、
    nu0 (νx0, νy0)・tw (αx, βx, αy, βy, μx, μy: 位置 loc)・harms ({調和: [(j,k,l,m,h)]}: 3 次と 4 次の写像の係数、補正前)"""

    def __init__(self, h3, h4, h3sum, A, nu0, tw, nf=True, nf_local=False, rowsinfo=None):
        harms = {}
        for q, h in zip(_MONOS3, h3):
            a = (q[0] - q[1], q[2] - q[3])
            if FR._canon(*a) and abs(h) > 0:
                harms.setdefault(a, []).append((q[0], q[1], q[2], q[3], complex(h)))
        for q, h in zip(_MONOS4, h4):
            a = (q[0] - q[1], q[2] - q[3])
            if a == (0, 0) or not FR._canon(*a):
                continue
            harms.setdefault(a, []).append((q[0], q[1], q[2], q[3], complex(h)))
        # 数値の 0 (対称性で消える調和の 1e-8 の丸め) の調和を除く
        hmax = max(abs(t[4]) for ts in harms.values() for t in ts)
        for a in list(harms):
            if max(abs(t[4]) for t in harms[a]) < 1e-6 * hmax:
                del harms[a]
        super().__init__(harms, A, nu0, tw)
        self.h3, self.h4, self.h3sum = np.asarray(h3, dtype=complex), np.asarray(h4, dtype=complex), np.asarray(h3sum, dtype=complex)
        self.nf, self.nf_local = nf, nf_local
        self.rowsinfo = rowsinfo or {}
        self._idx3 = {q: i for i, q in enumerate(_MONOS3)}
        self._idx4 = {q: i for i, q in enumerate(_MONOS4)}

    # ------------------------------------------------------------------
    @classmethod
    def from_lattice(cls, lat, loc, rows=None, delta=0.0, bend_cubic=None, nf=True, nf_local=False, order_adts=True, order2=True, amin2=1e-3, **kw):
        """Build the predictor from a lattice, a location and Twiss rows, using the RDTs (resonance_terms, rdt2nd.rdt4) and ADTS.
        lat: Lattice、loc: 位置 (1 始まり)、rows: lat.twiss_table() (省略可)。RDT は resonance_terms (3 次)・rdt2nd.rdt4 (4 次)、A は adts() の既定 (R6)。
        nf=True: 4 次に非共鳴の 3 次の補正 (非対角の正規形) を入れる。nf_local=True: 補正を初期条件ごとの ν(J) で評価する (False なら作動点 ν0)。δ ≠ 0 は未対応 (RDT は δ = 0)"""
        if delta != 0.0:
            raise NotImplementedError('S3 の予想は δ = 0 だけ (4 次の h が δ = 0 だけ)')
        if rows is None:
            rows = lat.twiss_table(delta=delta)
        h3d, _res = FR.first_order_terms(lat, rows, delta, bend_cubic, amin=1e-12)
        h3 = np.zeros(len(_MONOS3), dtype=complex)
        full3 = _R.resonance_terms(lat, terms=FR._MONOS3, rows=rows, delta=delta, bend_cubic=bend_cubic)
        for nm, q in zip(full3.names, full3.terms):
            h3[_MONOS3.index(tuple(q[:4]))] = full3.total[nm]
        r4 = _Q.rdt4(lat, rows=rows, order2=True)
        a = _R.adts(lat, rows=rows, delta=delta, bend_cubic=bend_cubic, fourth='r3', optics='calc')   # 物理的な項だけ (V2 で adts の既定をこれに戻す予定。既定に左右されないよう明示)
        tot = a.total
        A = [[tot['xx'], tot['xy']], [tot['xy'], tot['yy']]]
        r = rows[loc - 1]
        tw = (float(r.alfx), float(r.betx), float(r.alfy), float(r.bety), float(r.mux), float(r.muy))
        nu0 = (float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI)
        pr = cls(h3, r4.total_arr, r4.h3, A, nu0, tw, nf=nf, nf_local=nf_local)
        pr.warnings = list(r4.warnings)
        pr._info()
        return pr

    # ------------------------------------------------------------------
    def _info(self):
        """作動点の情報: 作動点から最も近い (1 次の RDT の) 3 次の共鳴線までの距離 (正規形の補正の分母の大きさ)"""
        qx, qy = self.nu0
        dist = {}
        for q, h in zip(_MONOS3, self.h3):
            k = (q[0] - q[1], q[2] - q[3])
            if k == (0, 0) or not FR._canon(*k) or abs(h) < 1e-9 * max(abs(self.h3).max(), 1e-300):
                continue
            v = k[0] * qx + k[1] * qy
            dist[k] = abs(v - round(v)) / math.hypot(*k)
        self.rowsinfo['d_first_order'] = dist
        self.rowsinfo['d_first_min'] = min(dist.values()) if dist else float('inf')

    # ------------------------------------------------------------------
    def tunes_model(self, Jx, Jy, quad=False):
        """Return the tune model nu(J) = nu0 + A.J with integer parts, where A is the default ADTS matrix.
        ν(J) = ν0 + A·J (整数部つき。A = adts() の既定)。J² の項は未実装 (摂動論の定義された次数は ADTS の J の 1 次まで)。quad は互換のための引数 (無視)"""
        return self.tunes_lin(Jx, Jy)

    def dnu_dJ(self, Jx, Jy, quad=False):
        """Return the derivative d nu_i / d J_j as an array built from the constant ADTS matrix A.
        ∂νi/∂Jj (定数 A を J の形に広げた (2, 2, ...) 配列)"""
        shp = np.shape(Jx)
        return np.broadcast_to(self.A.reshape((2, 2) + (1,) * len(shp)), (2, 2) + shp).copy()

    # ------------------------------------------------------------------
    def _nf_corr(self, nux, nuy):
        """3 次の項の非共鳴の補正 (35 成分) を、チューン (nux, nuy) (整数部を含む配列、1 点ごと) で評価。戻り (n, 35)"""
        nux = np.atleast_1d(nux).ravel()
        nuy = np.atleast_1d(nuy).ravel()
        out = np.zeros((nux.size, len(_MONOS4)), dtype=complex)
        cache = {}
        for i in range(nux.size):
            key = (round(float(nux[i]), 5), round(float(nuy[i]), 5))
            if key not in cache:
                cache[key] = _Q.nf_correction_h4(self.h3sum, key[0], key[1])
            out[i] = cache[key]
        return out

    def _F(self, a, Jx, Jy, corr):
        """調和 a の F(J) (配列)。corr: (n, 35) の補正 (None なら補正なし)"""
        F = 0j
        sh = np.shape(Jx)
        for (j, k, l, m, h) in self.harms[a]:
            ex, ey = (j + k) / 2.0, (l + m) / 2.0
            term = (2 * Jx) ** ex * (2 * Jy) ** ey
            if j + k + l + m == 4 and corr is not None:
                c = corr[:, self._idx4[(j, k, l, m)]].reshape(sh)
                F = F + (h + c) * term
            else:
                F = F + h * term
        return F

    def predict(self, x, y, nu=None, quad=True, nf=None, nf_local=None):
        """Predict resonance trapping for initial conditions (x, y) [m]; return a dict of scores, envelopes and reliability.
        初期条件 (x, y) [m] (同じ形の配列、px = py = 0) の予想。nu: (νx, νy) (整数部つき、配列) を渡すと ν(J) を実測のチューンにする (G0 は J のモデルの値)。
        (quad は無視)。戻り dict: any (位相つきで島に捕らえられる)、envelope (|Δ0| < W)、score (log10(W/|Δ0|) の調和にわたる最大。≥ 0 で包絡の中)、
        n_env (包絡に入る調和の数)、score_order (最大の調和の次数 |a1|+|a2|)、reliability (0/1/2、`reliability()`)、by_harmonic {a: dict(trapped, envelope, delta0, W, p, G0, g, score)}、
        nu_model (ν(J) の配列 2 つ)"""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        shape = x.shape
        nf = self.nf if nf is None else nf
        nf_local = self.nf_local if nf_local is None else nf_local
        Jx, Jy, phx, phy = self.actions(x, y)
        nux, nuy = (self.tunes_model(Jx, Jy, quad) if nu is None else nu)
        nuxm, nuym = self.tunes_model(Jx, Jy, quad)
        dnu = self.dnu_dJ(Jx, Jy, quad)
        corr = None
        if nf:
            if nf_local:
                corr = self._nf_corr(nux, nuy)
            else:
                corr = np.tile(_Q.nf_correction_h4(self.h3sum, *self.nu0), (Jx.size, 1))
        anyt = np.zeros(shape, dtype=bool)
        anye = np.zeros(shape, dtype=bool)
        nenv = np.zeros(shape, dtype=int)
        best = np.full(shape, -np.inf)
        best_order = np.zeros(shape, dtype=int)
        by = {}
        for a in self.harms:
            g = math.gcd(abs(a[0]), abs(a[1]))
            m0, n0 = a[0] // g, a[1] // g
            G0 = m0 * m0 * dnu[0, 0] + m0 * n0 * (dnu[0, 1] + dnu[1, 0]) + n0 * n0 * dnu[1, 1]
            s = (m0 * nux + n0 * nuy) * g
            p = np.round(s)
            d0 = (s - p) / g
            with np.errstate(divide='ignore', invalid='ignore'):
                dI = np.where(np.abs(G0) > 1e-12, -d0 / G0, np.nan)
            Jxr, Jyr = Jx + dI * m0, Jy + dI * n0
            okr = (Jxr > 0) & (Jyr > 0)
            F = self._F(a, Jx, Jy, corr)
            V = 2.0 * np.abs(F)
            th = a[0] * phx + a[1] * phy + np.angle(F)
            sg = np.where(G0 > 0, 1.0, -1.0)
            with np.errstate(divide='ignore', invalid='ignore'):
                en = math.pi * d0 * d0 / np.abs(G0)
                trapped = okr & (en < V * (1.0 - sg * np.cos(th)))
                W = np.sqrt(2.0 * np.abs(G0) * V / math.pi)
                env = okr & (np.abs(d0) < W)
                sc = np.where(okr, np.log10(np.maximum(W, 1e-300) / np.maximum(np.abs(d0), 1e-300)), -np.inf)
            by[a] = dict(trapped=trapped, envelope=env, delta0=d0, W=W, p=p, G0=G0, g=g, score=sc, V=V)
            anyt |= trapped
            anye |= env
            nenv += env.astype(int)
            upd = sc > best
            best = np.where(upd, sc, best)
            best_order = np.where(upd, abs(a[0]) + abs(a[1]), best_order)
        out = dict(any=anyt, envelope=anye, score=best, n_env=nenv, score_order=best_order, by_harmonic=by, nu_model=(nuxm, nuym), Jx=Jx, Jy=Jy)
        out['reliability'] = self.reliability(out, Jx, Jy)
        return out

    # ------------------------------------------------------------------
    def reliability(self, m, Jx, Jy, dnu_max=0.02, d_first=0.01, dline=0.003):
        """Return the per-point reliability (0 unreliable, 1 caution, 2 reliable) of the prediction.
        点ごとの信頼度 (0 = 信頼できない、1 = 注意、2 = 信頼できる)。判定の基準 (docs/25 §12.5):
          (1) 線形の ν(J) の外れ: |ν_model − ν_lin| (J² の項の寄与) の最大 > 0.005 → 注意、> 0.02 または |ν_lin − ν0| > 0.1 → 信頼できない
              (J² の項を fit していないときは線形の移動 |A·J| で判定: > 0.02 → 注意、> 0.05 → 信頼できない)
          (2) 共鳴の重なり: 包絡に入る調和が 2 つ以上 (カオスの領域) → 信頼できない
          (3) 駆動項が未知の共鳴: 点の ν(J) から次数 5, 6 の共鳴線 (駆動項の解析式が無い) までの距離 < dline → 注意
          (4) 作動点が 1 次の RDT の共鳴線の近く (距離 < d_first): 正規形の補正が使えない → 全点 注意"""
        nux, nuy = m['nu_model']
        shp = np.shape(Jx)
        lin = self.tunes_lin(Jx, Jy)
        shift_lin = np.hypot(lin[0] - self.nu0[0], lin[1] - self.nu0[1])
        rel = np.full(shp, 2, dtype=int)
        rel = np.where(shift_lin > 0.02, np.minimum(rel, 1), rel)
        rel = np.where(shift_lin > 0.05, 0, rel)
        rel = np.where(m['n_env'] >= 2, 0, rel)
        d56 = FR.lines_distance(nux, nuy, _mn_exact_orders((5, 6)))
        rel = np.where(d56 < dline, np.minimum(rel, 1), rel)
        if self.rowsinfo.get('d_first_min', 1.0) < d_first:
            rel = np.minimum(rel, 1)
        return rel

    def warnings_for(self, m, Jx=None, Jy=None):
        """Return the list of warning strings for the whole prediction.
        予想全体の警告 (文字列のリスト)"""
        w = []
        rel = m['reliability']
        n = rel.size
        f0, f1 = float((rel == 0).sum()) / n, float((rel == 1).sum()) / n
        if self.rowsinfo.get('d_first_min', 1.0) < 0.01:
            w.append('作動点が 1 次の RDT の共鳴線の近く (距離 %.4f): 非共鳴の補正 (nf) が使えない' % self.rowsinfo['d_first_min'])
        if f0 > 0.3:
            w.append('初期条件の %.0f%% が「信頼できない」(共鳴の重なり・大きな振幅依存のチューンの動き)' % (100 * f0))
        w.append('ν(J) は線形の ADTS だけ (J² の項は未実装): |ν − ν0| が 0.02 を超える振幅は注意、0.05 を超えると信頼できない')
        if f1 + f0 > 0.9:
            w.append('ほぼ全域が注意以下 (信頼できる点 %.0f%%)' % (100 * (1 - f0 - f1)))
        return w


_MN_CACHE = {}


def _mn_exact_orders(orders):
    """次数がちょうど orders のどれかの (m, n) (既約でなくてもよい。(m, n) と (−m, −n) は 1 つ)"""
    key = tuple(orders)
    if key not in _MN_CACHE:
        from .freqmap import _mn_list
        top = max(orders)
        _MN_CACHE[key] = [mn for mn in _mn_list(top) if abs(mn[0]) + abs(mn[1]) in orders]
    return _MN_CACHE[key]


# ---------------------------------------------------------------------------------------------------------------------
# 判別力の評価 (cov・面積・対照・AUC・順位相関・上位 10%)
# ---------------------------------------------------------------------------------------------------------------------
def _rank(v):
    from math import isnan
    order = np.argsort(v, kind='mergesort')
    ranks = np.empty(len(v), dtype=float)
    ranks[order] = np.arange(1, len(v) + 1)
    # 同順位は平均
    vs = v[order]
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and vs[j + 1] == vs[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = 0.5 * (i + j) + 1
        i = j + 1
    return ranks


def auc_score(score, positive):
    """Compute the AUC of the ROC curve: the probability that a positive score exceeds a negative one (ties count 0.5).
    ROC 曲線の AUC (正例の score が負例より大きい確率。同点は 0.5)。score: 配列、positive: 真偽の配列。正例か負例が無ければ nan"""
    score = np.asarray(score, dtype=float)
    pos = np.asarray(positive, dtype=bool)
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float('nan')
    rk = _rank(score)
    return float((rk[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def spearman(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ra, rb = _rank(a), _rank(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = math.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else float('nan')


def partial_spearman(a, b, c):
    """Compute the partial rank correlation of a and b with the rank of c removed.
    c (例: 振幅) の順位を除いた a と b の順位偏相関"""
    rab, rac, rbc = spearman(a, b), spearman(a, c), spearman(b, c)
    den = math.sqrt(max((1 - rac ** 2) * (1 - rbc ** 2), 1e-300))
    return float((rab - rac * rbc) / den)


def score_stats(D, score, thr=-5.0, top=0.1, amp=None, nboot=0, seed=1):
    """Evaluate scores (larger means stronger resonance) against D: AUC, Spearman rho, and amplitude-based statistics.
    スコア (大きいほど共鳴が強い) の評価。D (nan は除く)。戻り dict: auc (D > thr の点を正例)、rho (D との Spearman)、top の D の中央値・D > thr の割合 (全体の割合と)、
    amp (振幅の指標) があれば amp だけのスコアの auc・rho と、amp の順位を除いた偏相関 rho_partial、amp の層ごと (3 層) の auc"""
    D = np.asarray(D, dtype=float)
    s = np.asarray(score, dtype=float)
    ok = np.isfinite(D) & ~np.isnan(s)
    D, s = D[ok], np.where(np.isinf(s[ok]), -50.0, s[ok])
    pos = D > thr
    out = dict(n=int(ok.sum()), n_pos=int(pos.sum()), auc=auc_score(s, pos), rho=spearman(s, D))
    k = max(int(round(top * len(D))), 1)
    idx = np.argsort(-s, kind='mergesort')[:k]
    out.update(top_n=k, top_Dmed=float(np.median(D[idx])), top_frac_pos=float(pos[idx].mean()), base_frac_pos=float(pos.mean()), base_Dmed=float(np.median(D)))
    if amp is not None:
        a = np.asarray(amp, dtype=float)[ok]
        out.update(amp_auc=auc_score(a, pos), amp_rho=spearman(a, D), rho_partial=partial_spearman(s, D, a))
        qs = np.quantile(a, [1 / 3, 2 / 3])
        layers = np.digitize(a, qs)
        out['auc_layers'] = [auc_score(s[layers == i], pos[layers == i]) for i in range(3)]
        out['rho_layers'] = [spearman(s[layers == i], D[layers == i]) for i in range(3)]
    return out


def evaluate_prediction(r, m, thr=-5.0, key='envelope', nshift=400, amp=None):
    """Evaluate a prediction against an FMapResult r, returning coverage statistics and score statistics.
    FMapResult r に対する予想 m (`ResonancePredictorS3.predict` の戻り) の判別力: coverage_stats (領域 key の cov・面積・対照・p) と score_stats。
    amp 省略時は規格化した振幅 √(Jx/Jxmax + Jy/Jymax) (m['Jx'], m['Jy'] から)。戻り {'cov': coverage_stats の dict, 'score': score_stats の dict}"""
    cs = FR.coverage_stats(r.D, m[key], r.grid, thr=thr, nshift=nshift)
    if amp is None:
        Jx, Jy = np.asarray(m['Jx']).ravel(), np.asarray(m['Jy']).ravel()
        amp = np.sqrt(Jx / max(Jx.max(), 1e-300) + Jy / max(Jy.max(), 1e-300))
    ss = score_stats(r.D, np.asarray(m['score']).ravel(), thr=thr, amp=amp)
    return dict(cov=cs, score=ss)


def line_distance_score(nux, nuy, N=4):
    """Return the comparison score -log10(distance to the nearest resonance line of order N or lower), using proximity only.
    比較用 (予想の強さを使わない) のスコア: −log10(次数 N 以下の共鳴線までの距離)。線への近さだけ"""
    from .freqmap import _mn_list
    d = FR.lines_distance(nux, nuy, _mn_list(N))
    return -np.log10(np.maximum(d, 1e-9))


def nearest_order_stats(D, nux, nuy, thr=-5.0, tol=0.003, orders=(1, 2, 3, 4, 5, 6)):
    """Find the lowest resonance order of lines within tol of each measured tune and return its distribution for D > thr and all points.
    各点の測ったチューン (端数でよい) から、距離 < tol の共鳴線のうち最小の次数 |m|+|n| を求め、D > thr の点と全点の次数の分布を返す。
    戻り dict: orders (配列)、frac_hi・frac_all (最小次数が k の割合、'none' = どの次数 ≤ max にも近くない)、enrich (hi/all)"""
    D = np.asarray(D, dtype=float)
    ok = np.isfinite(D) & np.isfinite(nux)
    hi = ok & (D > thr)
    nmin = np.full(D.shape, 99)
    for k in orders:
        d = FR.lines_distance(nux, nuy, _mn_exact_orders((k,)))
        nmin = np.where((d < tol) & (nmin == 99), k, nmin)
    out = dict(orders=list(orders) + ['none'], n_hi=int(hi.sum()), n_all=int(ok.sum()))
    fh, fa = [], []
    for k in list(orders) + [99]:
        fh.append(float((nmin[hi] == k).sum()) / max(hi.sum(), 1))
        fa.append(float((nmin[ok] == k).sum()) / max(ok.sum(), 1))
    out['frac_hi'], out['frac_all'] = fh, fa
    out['enrich'] = [h / a if a > 0 else float('nan') for h, a in zip(fh, fa)]
    return out
