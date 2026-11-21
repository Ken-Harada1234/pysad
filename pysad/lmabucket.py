# -*- coding: utf-8 -*-
"""RF バケットの限界とアパーチャの限界の区別 (v2.3.0 L2)。局所運動量アクセプタンス (lma.py) の限界が、どちらで決まるかを独立に調べる。

LMA の限界は 2 つの限界の小さい方になる:
  (a) RF バケットの限界 δ_bucket: 縦方向の運動だけで、バケットの外へ出る (位相がすべる) δ。横のアパーチャを無限大にして追跡する
      (`bucket_limits_6d`: 6 次元の追跡)。同じ限界を、縦方向の 1 次元の写像で独立に求める (`bucket_limits_1d`): 解析式の RF の力 (V/E0) sin(φs + k z) と、
      追跡で測った非線形の経路長 Δz = −C (α1 δ + α2 δ²)、放射による δ の変化 (δ の 1 次・2 次) だけを使う。解析式の δ_rf (`rf_bucket_analytic`) は、
      α が一定で U0 が一定の理想化したセパラトリクス (Emittance[] の bucket_height と同じ式) で、放射の減衰・α2 を含まない (それとの差は ±1〜2%)。
  (b) アパーチャの限界 δ_aper: RF・放射なしの 4 次元の追跡 (δ を固定)。位置 loc の閉軌道 (δ = 0) から δ だけずれて出発した粒子が、ベータトロン振動 + 分散の
      ずれで横のアパーチャに当たる最小の δ。縦の運動 (RF のすべり) を含まない。
  `limit_types(setup, result)` は、LMA の結果の位置ごとに (a)・(b) を求めて min(a, b) と LMA の値を比べ、種類 ('rf'・'aper'・'both') を決める。
  LMAResult.reason は「最初に失われた粒子を落とした直接の原因 (|z| か横のアパーチャ)」で、限界の種類とは別: バケットの限界の位置でも、バケットから出た粒子は
  すべって δ が変わり、分散の大きい位置で横のアパーチャに当たって落ちるので 'aper' と記録されうる。
"""
import math

import numpy as np

from . import trksetup as _ts


def rf_bucket_analytic(setup):
    """理想化した RF バケットの高さ (解析式)。α 一定、U0 一定、放射の減衰なし。H = −(kCα/2) δ² + A (cos ψ + ψ sin φs)、A = V/E0 (E0 = p0 c)、ψ = φs + k z。
    不安定固定点 ψ = π − φs を通るセパラトリクスの、安定固定点での δ:
        δ_rf² = (2A/(k C α)) [2 cos φs − (π − 2 φs) sin φs]
    戻り dict: delta_rf、A (= V/E0)、k_rf、C、alpha、phis、h (調和数 kC/2π)、zlimit (バケットの |z| の限界 (π + 2|φs|)/k)、emit_bucket_height (Emittance[] の値)"""
    b = setup.bucket
    if not setup.rfsw or not b or b.get('k_rf', 0.0) <= 0.0:
        return None
    k, C, al, phis = b['k_rf'], setup.circ, b['alpha'], b['phis']
    A = b['Vrf'] / setup.p0
    val = 2.0 * A / (k * C * al) * (2.0 * math.cos(phis) - (math.pi - 2.0 * phis) * math.sin(phis)) if al != 0.0 else float('nan')
    return dict(delta_rf=math.sqrt(abs(val)), A=A, k_rf=k, C=C, alpha=al, phis=phis, h=k * C / (2 * math.pi),
                zlimit=(math.pi + 2.0 * abs(phis)) / k, emit_bucket_height=b['height'])


def _cavity_exit_loc(setup):
    """最後の空洞 (Cavity) の次の位置 (loc: 1 始まり、空洞のすぐ後ろ)。なければ None"""
    idx = [i for i, e in enumerate(setup.els) if type(e).__name__ == 'Cavity']
    if not idx:
        return None
    return idx[-1] + 2


def longitudinal_params(setup, dps=None, loc=None):
    """追跡で測る縦方向の写像のパラメーター (空洞のすぐ後ろの位置から 1 周):
      a1・a2: 1 周の経路長 Δz = z' − z = −C (a1 δ + a2 δ²) (RF ON・放射なし、δ の閉軌道の上、RF の力は周の最後の空洞で 1 回だけなので経路長に入らない)。a1 ≈ αc
      r1・r2・r0: 放射による 1 周の δ の変化 (RF なし。空洞は U0 の定数で補償する SAD の NORFSW の扱いなので、δ の依存だけが出る): Δδ = r0 + r1 δ + r2 δ²
                  (−r1 が 1 周あたりの δ の減衰率、理想的な U ∝ (1 + δ)² なら −r1 = 2U0/E0、r2 = −U0/E0)
      kick_fit: 空洞の力 (V/E0) sin(φs + k z) − U0/E0 を z の格子で追跡して当てはめた (V/E0, φs, U0/E0) と、設計値との最大の差 (解析式の RF の確認)
    δ の範囲 dps (既定: ±バケット高さの 1.3 倍、9 点)。RF なしなら None"""
    b = setup.bucket
    if not setup.rfsw or not b or b.get('k_rf', 0.0) <= 0.0:
        return None
    if loc is None:
        loc = _cavity_exit_loc(setup)
        if loc is None:
            return None
    if dps is None:
        dps = np.linspace(-1.3, 1.3, 9) * b['height']
    dps = np.asarray(dps, dtype=float)
    C = setup.circ
    s4 = setup.copy(rad=False, rfsw=False)
    sR = setup.copy(rad=False, rfsw=True)
    s5 = setup.copy(rad=True, rfsw=False)
    dz, dr, used, skipped = [], [], [], []
    for d in dps:
        try:
            c0 = s4.closed_orbit_dp(loc, float(d))
            big = dict(alost=1e30, zlost=1e30)
            z1 = sR.track(loc, c0[None, :], 1, **big)['final'][0][4] - c0[4]
            r1 = s5.track(loc, c0[None, :], 1, **big)['final'][0][5] - d
        except Exception:       # 大きな δ の閉軌道が求まらない (4 次元の Newton 法が収束しない・粒子が失われる)
            skipped.append(float(d))
            continue
        if not (np.isfinite(z1) and np.isfinite(r1)):
            skipped.append(float(d))
            continue
        dz.append(z1)
        dr.append(r1)
        used.append(float(d))
    if len(used) < 5:
        raise RuntimeError('longitudinal_params: closed orbits found for only %d of %d δ values' % (len(used), len(dps)))
    dps = np.array(used)
    pz = np.polyfit(dps, np.array(dz), 2)
    pr = np.polyfit(dps, np.array(dr), 2)
    # 空洞の力の確認
    k = b['k_rf']
    zs = np.linspace(-3.0, 3.0, 13) / k
    cod = setup.closed_orbit([loc])[0]
    c = np.tile(cod, (len(zs), 1))
    c[:, 4] += zs
    r = setup.track(loc, c, 1, alost=1e30, zlost=1e30)
    kick = r['final'][:, 5] - c[:, 5]
    A, phis, u0 = b['Vrf'] / setup.p0, b['phis'], b['U0'] / setup.p0
    pred = A * np.sin(phis + k * zs) - u0
    # 1 周の始め・終わりの位置が空洞の後ろなので、1 周の経路で δ が変わる (放射) 分は kick に含まれる: 差を見る
    return dict(loc=loc, C=C, a1=float(-pz[1] / C), a2=float(-pz[0] / C), a0=float(pz[2]), r1=float(pr[1]), r2=float(pr[0]), r0=float(pr[2]),
                a1_over_alpha=float(-pz[1] / C / b['alpha']) if b.get('alpha') else float('nan'),
                kick_max_diff=float(np.max(np.abs(kick - pred))), kick_scale=float(A), dps=dps, skipped=skipped)


def _map_limit(par, setup, nturn, sign, dp_top, ngrid=600, nref=48, levels=3, a2=None, r_on=True):
    """縦方向の 1 次元の写像の限界 (1 つの符号)。z* (安定固定点) から δ0 = sign·d で出発し、|z| > バケットの |z| の限界 になる最小の d (外向きに最初に失われる点)。
    戻り (d_limit, lost_turn): 失われる d がなければ (None, 0)"""
    b = setup.bucket
    k, C = b['k_rf'], setup.circ
    A, phis, u0 = b['Vrf'] / setup.p0, b['phis'], b['U0'] / setup.p0
    zb = setup.bucket_zlimit(1.0)
    a1 = par['a1']
    a2 = par['a2'] if a2 is None else a2
    r1, r2 = (par['r1'], par['r2']) if r_on else (0.0, 0.0)

    def run(dgrid):
        d = sign * dgrid.copy()
        z = np.full_like(d, z0)
        alive = np.ones(len(d), bool)
        lost = np.zeros(len(d), dtype=np.int64)
        for n in range(1, nturn + 1):
            d = d + A * np.sin(phis + k * z) - u0 + (r1 * d + r2 * d * d)
            z = z - C * (a1 * d + a2 * d * d)
            out = alive & (np.abs(z) > zb)
            lost[out] = n
            alive &= ~out
            if not alive.any():
                break
        return lost

    # 安定固定点 (δ*, z*): d = d + A sin(φs + kz) − u0 + r(d)、z = z − C a(d)  → a(δ*) = 0 かつ A sin(φs + kz*) = u0 − r(δ*)。δ* ≈ 0、z* ≈ 0
    z0 = (math.asin(max(-1.0, min(1.0, (u0 - par['r0']) / A))) - phis) / k
    grid = np.linspace(0.0, dp_top, ngrid + 1)
    lost = run(grid)
    w = np.nonzero(lost > 0)[0]
    if not len(w) or w[0] == 0:
        return (None, 0) if not len(w) else (0.0, int(lost[0]))
    lo, hi = grid[w[0] - 1], grid[w[0]]
    for _ in range(levels):
        g = np.linspace(lo, hi, nref + 1)
        lost = run(g)
        w = np.nonzero(lost > 0)[0]
        lo, hi = g[w[0] - 1], g[w[0]]
    return float(lo), int(lost[w[0]])


def bucket_limits_1d(setup, nturn, par=None, dp_top=None, **kw):
    """縦方向の 1 次元の写像で求めた RF バケットの限界 (δ+, |δ−|) (nturn 周で |z| がバケットを出ない最大の δ)。LMA の損失判定の |z| と同じ基準。
    par: longitudinal_params の戻り (None なら測る)。戻り dict: dp_plus、dp_minus (負で)、par。a2=0・r_on=False などで項の寄与を切り分けられる"""
    if par is None:
        par = longitudinal_params(setup)
    if par is None:
        return None
    if dp_top is None:
        dp_top = 1.6 * setup.bucket['height']
    p, _ = _map_limit(par, setup, nturn, +1, dp_top, **kw)
    m, _ = _map_limit(par, setup, nturn, -1, dp_top, **kw)
    return dict(dp_plus=p, dp_minus=-m if m is not None else None, par=par)


def bucket_limits_6d(setup, locs, nturn=None, workers=1, **kw):
    """6 次元の追跡の RF バケットの限界: 横のアパーチャを無限大 (alost = 1e30) にして、|z| がバケットを出る (LMA の bucket 判定) δ を位置ごとに求める (momentum_acceptance を使う)"""
    from .lma import momentum_acceptance
    kw.setdefault('check_every', 0)
    kw.setdefault('verify', 0)
    kw.setdefault('turncheck', 0)
    return momentum_acceptance(setup, locs, nturn=nturn, alost=1e30, workers=workers, quiet=True, **kw)


def aperture_limits_4d(setup, locs, nturn=256, dpmax=None, workers=1, check_every=1, **kw):
    """横のアパーチャの限界: RF・放射なしの 4 次元の追跡 (δ を固定)。nturn 周のあいだ |x|,|y| < alost に収まる最大の δ を位置ごとに求める (RF のすべりを含まない)"""
    from .lma import momentum_acceptance
    s4 = setup.copy(rad=False, rfsw=False)
    if dpmax is None:
        dpmax = 1.25 * setup.bucket.get('height', 0.04) if setup.bucket else 0.05
    kw.setdefault('verify', 0)
    kw.setdefault('turncheck', 0)
    return momentum_acceptance(s4, locs, nturn=nturn, dpmax=dpmax, bucket=False, workers=workers, check_every=check_every, quiet=True, **kw)


class LimitTypes:
    """limit_types の結果。配列 (位置の順): loc、s、lma_plus・lma_minus (LMA の値、δ− は負)、rf_plus・rf_minus (6 次元の RF バケットの限界)、
    ap_plus・ap_minus (4 次元のアパーチャの限界、δ− は負)、ap_eff (実効のアパーチャの限界 = min(δ+_aper, |δ−_aper|))、expect_plus・expect_minus (min(rf, ap_eff))、type_plus・type_minus ('rf'・'aper'・'both')、
    reason_plus・reason_minus (LMA の直接の損失の理由)。dev_plus・dev_minus (LMA − expect)。par (1 次元の写像のパラメーター)、
    map1d (1 次元の写像の限界 {dp_plus, dp_minus}、LMA の nturn)、analytic (rf_bucket_analytic)"""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def table(self):
        L = ['# loc  s[m]  name  LMA+  rf+  aper+  type+  reason+ | LMA-  rf-  aper-  type-  reason-  | aper_eff  (rf = RF bucket limit (6D, no aperture), aper = aperture limit (4D, RF/rad off), aper_eff = min(aper+, |aper-|))']
        for i in range(len(self.loc)):
            L.append('%5d %10.4f %-8s %9.5f %9.5f %9.5f %-5s %-5s | %9.5f %9.5f %9.5f %-5s %-5s | %9.5f' % (
                self.loc[i], self.s[i], self.name[i], self.lma_plus[i], self.rf_plus[i], self.ap_plus[i], self.type_plus[i], self.reason_plus[i],
                self.lma_minus[i], self.rf_minus[i], self.ap_minus[i], self.type_minus[i], self.reason_minus[i], self.ap_eff[i]))
        return '\n'.join(L) + '\n'


def limit_types(setup, result, nturn_aper=256, workers=1, tol=None, rf6=None, ap4=None):
    """LMA の結果 (LMAResult) の位置ごとに、RF バケットの限界 (6 次元) とアパーチャの限界 (4 次元) を求め、どちらで決まっているか (type) を返す。
    type: 'rf' = バケットの限界の方が実効のアパーチャの限界より小さい (差 > 2 tol)、'aper' = 実効のアパーチャの限界の方が小さい、'both' = 差が 2 tol 以内。
    LMA の値が min(rf, aper) と tol の数倍以内で合えば、損失の原因は 2 つの限界で説明できる (dev)。rf6・ap4: 計算済みの LMAResult を渡すと再計算しない"""
    c = result.cond
    tol = c.get('tol', 1e-4) if tol is None else tol
    nturn = int(c['nturn'])
    stp = result.setup
    if rf6 is None:
        rf6 = bucket_limits_6d(stp, result.loc, nturn=nturn, workers=workers, tol=tol, nsynch=c.get('nsynch_default', 2.0))
    if ap4 is None:
        ap4 = aperture_limits_4d(stp, result.loc, nturn=nturn_aper, workers=workers, tol=tol)
    par = longitudinal_params(stp)
    m1 = bucket_limits_1d(stp, nturn, par=par) if par is not None else None
    # 6 次元の粒子は、シンクロトロン振動で半周期のあいだに δ の符号が反転するので、反対側のアパーチャの限界にも当たる:
    # 実効のアパーチャの限界 = min(δ+_aper, |δ−_aper|) (両符号に共通)。実際の LMA は、反対側の限界を数% 超えることがある (ベータトロン振幅・位相の違い)
    eff = np.minimum(ap4.dp_plus, -ap4.dp_minus)
    ex_p = np.minimum(rf6.dp_plus, eff)
    ex_m = np.maximum(rf6.dp_minus, -eff)

    def typ(a, b):
        return ['both' if abs(x - y) <= 2 * tol else ('rf' if x < y else 'aper') for x, y in zip(a, b)]
    return LimitTypes(loc=result.loc, s=result.s, name=result.name, lma_plus=result.dp_plus, lma_minus=result.dp_minus,
                      rf_plus=rf6.dp_plus, rf_minus=rf6.dp_minus, ap_plus=ap4.dp_plus, ap_minus=ap4.dp_minus,
                      expect_plus=ex_p, expect_minus=ex_m, type_plus=typ(rf6.dp_plus, eff),
                      type_minus=typ(-rf6.dp_minus, eff), ap_eff=eff, reason_plus=list(result.reason_plus), reason_minus=list(result.reason_minus),
                      dev_plus=result.dp_plus - ex_p, dev_minus=result.dp_minus - ex_m, par=par, map1d=m1, analytic=rf_bucket_analytic(stp),
                      rf6=rf6, ap4=ap4, nturn=nturn, nturn_aper=nturn_aper)
