# -*- coding: utf-8 -*-
"""
RADTAPER (K64 の tcod.f・tturne.f・tturn.f の radtaper): 放射で運動量が周の途中で変わる (エネルギーの鋸歯) のに合わせて、
磁石の強さを閉軌道の局所の運動量に比例させる。フラグ RAD・RADCOD・RADTAPER・RFSW がすべて ON のとき Emittance[]・EMIT と
TrackParticles に効く。

K64 の定義 (収束した状態):
  - 閉軌道の入口の運動量のずれは δ(1) = DP0 に固定し、その代わりに一様な係数 dptaper を未知数にして 6 次元で閉じる (tcod.f)。
  - 要素 l の係数 (tturne.f。v2 は BEND も (δ_l + δ_{l+1})/2): BEND は 1 + (3 δ_l + δ_{l+1})/4 + dptaper、QUAD・薄い多極・MULT は 1 + (δ_l + δ_{l+1})/2 + dptaper
    (δ_l は要素 l の入口の閉軌道の δ。±0.3 で頭打ち)。BEND は (K0 + ANGLE) と K1 に掛ける (幾何の ANGLE は変えない)。
  - 追跡 (tturn.f) は Twiss の DDP (= δ_l) と dptaper から、BEND も含めて 1 + (δ_l + δ_{l+1})/2 + dptaper。
  - Twiss["DDP", l] に δ_l が入る。

実装: K64 の tcod の反復 (1 周ごとに係数を作り直す Newton 法) をそのまま写さず、同じ不動点を外側の反復で求める。
  係数を t_l = 1 + (δ の形 s_l の平均) + u (s_l = δ_l − δ_1、u = DP0 + dptaper) と書き、
  [係数を掛けた要素の列で 6 次元の閉軌道 (既存の closed_orbit6d) → δ の分布を測る → s を更新、u を δ_1 = DP0 になるよう割線法で更新]
  を収束まで繰り返す。CAVI・SOL・ドリフトには掛けない。
限界: CALC (4 次元の光学) には掛けない (K64 の qtwiss の radtaper は未翻訳)。SOL 領域内の要素には掛けない。
"""
from __future__ import annotations

import copy
import math

from . import compat as _compat

import numpy as np

TAPMAX = 0.3
_THIN = {'Sextupole': 'K2', 'Octupole': 'K3', 'Decapole': 'K4', 'Dodecapole': 'K5'}


def _clip(t):
    return min(1.0 + TAPMAX, max(1.0 - TAPMAX, t))


def scaled(el, t):
    """要素 el の磁場を t 倍にした写し (掛ける対象でなければ el そのもの)"""
    if t == 1.0 or hasattr(el, '_solmap'):
        return el
    tn = type(el).__name__
    if tn == 'Bend':
        e = copy.copy(el)
        ang = float(el.ANGLE)
        k0 = float(getattr(el, 'K0', 0.0) or 0.0)
        e.K0 = (k0 + ang) * t - ang
        e.K1 = float(el.K1) * t
        return e
    if tn == 'Quadrupole':
        e = copy.copy(el)
        e.K1 = float(el.K1) * t
        return e
    if tn in _THIN:
        e = copy.copy(el)
        k = _THIN[tn]
        setattr(e, k, float(getattr(el, k)) * t)
        return e
    if tn == 'Multipole':
        e = copy.copy(el)
        knl = {k: v * t for k, v in (getattr(el, 'KnL', {}) or {}).items()}
        ang = float(getattr(el, 'ANGLE', 0.0) or 0.0)
        if ang != 0.0:
            knl[0] = knl.get(0, 0.0) + (t - 1.0) * ang          # tmulta.f: ak(0) = ak0(0)·rtaper + (rtaper − 1)·phi
        if knl:
            e.KnL = knl
        if getattr(el, 'KsL', None):
            e.KsL = {k: v * t for k, v in el.KsL.items()}
        return e
    return el


def factors(els, ddp, u, tracking=False):
    """係数 t_l (要素ごと)。ddp[l] は要素 l (0 始まり) の入口の δ − δ(入口)、ddp[len(els)] は出口"""
    out = []
    for l, el in enumerate(els):
        a, b = ddp[l], ddp[l + 1]
        if type(el).__name__ == 'Bend' and not tracking and not _compat.v2():
            # K64 の tturne の BEND は (3 δ_l + δ_l+1)/4 (放射を刻みの終わりに入れる粗い刻みに合わせた重み)、追跡は (δ_l + δ_l+1)/2。
            # v2 は Emittance[]・EMIT・追跡とも (δ_l + δ_l+1)/2 (連続な放射の極限で、要素の中の損失が線形のとき曲げ角を設計値にする係数は
            # 調和平均 ≈ 算術平均。v2next P2 項目 3)。同じ係数なので Emittance[] の閉軌道は追跡で閉じる
            t = 1.0 + (3.0 * a + b) * 0.25 + u
        else:
            t = 1.0 + (a + b) * 0.5 + u
        out.append(_clip(t))
    return out


def tapered(els, ddp, u, tracking=False):
    return [scaled(el, t) for el, t in zip(els, factors(els, ddp, u, tracking))]


def _profile(els, cod, ph, trf0):
    """閉軌道に沿った δ (要素の入口ごと、最後は出口)。放射込みの COD だけの追跡を要素ごとに"""
    from . import emit6d as _e
    out = [float(cod[5])]
    c = np.array(cod, dtype=float)
    for el in els:
        st = _e.tturne([el], c, ph, trf0=trf0, rfsw=True, radcod=True, want_jac=False)
        c = np.array(st.cod, dtype=float)
        out.append(float(c[5]))
    return np.array(out)


def solve(els, momentum_eV, dp0=0.0, fshift=0.0, charge=1.0, amass=None, tol=1e-10, maxit=30):
    """戻り値 dict(ddp=δ_l (入口ごと + 出口。δ_1 = dp0)、u、dptaper、els=Emittance 用の要素の列、found、iterations)"""
    from . import emit6d as _e
    kw = {} if amass is None else {'amass': amass}
    ph = _e._Phys(momentum_eV, _e._emit_circ(els), charge=charge, **kw)
    ph.dvfs = fshift
    n = len(els)
    s = np.zeros(n + 1)
    u = dp0
    u_prev = err_prev = None
    found = False
    it = 0
    for it in range(1, maxit + 1):
        elt = tapered(els, s, u)
        cod, trf0, fnd, _ = _e.closed_orbit6d(elt, ph, True, True)
        if not fnd or not np.all(np.isfinite(cod)):
            break
        prof = _profile(elt, cod, ph, trf0)
        s_new = prof - prof[0]
        err = prof[0] - dp0
        ds = float(np.max(np.abs(s_new - s)))
        s = s_new
        if ds < tol and abs(err) < tol:
            found = True
            break
        # u の更新: 一様な係数 ε は閉軌道の δ をほぼ ε だけ動かす (傾き ≈ 1)。2 回目からは割線法
        slope = 1.0
        if u_prev is not None and u != u_prev and err != err_prev:
            sl = (err - err_prev) / (u - u_prev)
            if 0.2 < sl < 5.0:
                slope = sl
        u_prev, err_prev = u, err
        u = u - err / slope
    return dict(ddp=s + dp0, shape=s, u=u, dptaper=u - dp0, els=tapered(els, s, u), found=found, iterations=it)
