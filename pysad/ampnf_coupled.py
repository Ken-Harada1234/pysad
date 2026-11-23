# -*- coding: utf-8 -*-
"""v2.4.0 H1: 結合した光学 (ROTATE した磁石) の共鳴駆動項の拡張: 4 次の h、δ≠0、ref≠None、DX・DY・DZ・DROTATE (据え付けのずれ・ロール誤差)。
`ampnf.py` (A1・F1・G の領域) は編集せず、その部品 (`build_ops`・`_compose`・`linear_modes`・`Seg`) を呼ぶ別ファイル。

【方法 (入れている効果と摂動論の次数)】
  模型 (ampnf と同じ: 線形の厳密な転送行列 + 6 極・8 極の kick + 偏向磁石の曲率・端の 3 次の項 + 運動学 + フリンジ) の 1 周の写像を、(x, px, y, py, δ) の 5 変数の
  全次数 D までのべき級数 (TPSA) で作る (物理座標、リング先頭または ref の素子の入口から)。次に
    1. δ = δ0 (固定) として、閉軌道 x* (べき級数の不動点、Newton 法) を求め、その点のまわりに展開し直す (DX・DY・δ0 で閉軌道は 0 でない)、
    2. 展開した点の 1 周の線形写像から 4×4 の固有モード (`ampnf.linear_modes`) を作り、正規化座標 zn (x = x* + Am zn) の写像 R + T (2 次) + T3 (3 次) にする、
    3. その写像から `rdttrack` と同じ手順 (kick の生成関数 → 角度・作用の調和、3 次の Dragt–Finn の f4) で h3 (j+k+l+m = 3)・h4 (j+k+l+m = 4) を出す。
  h4 は「写像 = e^{−:h3 + h4:}」の定義の f4 (非共鳴の項を消す前の写像の係数。docs/25 §10)。含まれるもの: 8 極の 1 次、運動学、フリンジ (fringe=True)、3 次の項どうしの 2 次
  (6 極・偏向磁石の曲率・端 (bend_cubic=True) の Magnus の第 2 項)。**入っていないもの**: 偏向磁石の K1 の 4 次以上、K64 の四極の積分器の台形則 (模型は厳密な線形 + 運動学の流れ)。
  δ0≠0 は、模型の写像の δ のべき級数 (全次数 D の打ち切り: 運動学などの 4 次のハミルトニアンが入るので T3 の δ^p の係数は p ≤ D−4 まで正しい) を δ0 で評価したもの
  = δ0 の閉軌道のまわりの、そのエネルギーの光学 (β・結合・分散による軌道のずれを含む) の h。誤差は O(δ0^(D−3))。h の位相の基準は、リング先頭で x の成分が実で正のモード (`linear_modes`)。
  ref≠None (`ref_phases`・`rotate_to_ref`) は resdrive の ref と同じ定義 = 位相の基準 (角度の原点だけを動かす。ref の素子の入口でのモードの位相進みを h に掛ける。|h| は変わらない)。
  `taylor_map(start=)` は別の量で、1 周の写像の始点を素子 start にする (素子列を循環させる。厳密。始点が違うと非共鳴の部分の h の大きさも変わる)。

【DX・DY・DZ・DROTATE (SAD の tchge)】 `lattice._bend_chg_args` と同じ引数で、tchge の入口・出口をそのまま TPSA にする: 回転 (x' = c x − s y、px・py も) → y のずれ → tbrote
  (x のずれ dx1 = dxr c0 − dz s0、dz1 = dz c0 + dxr s0 のドリフト (tdrife1)、DROTATE≠0 のとき tbgrote の 3 次元の回転 (px・py → pz の平方根の展開を含む))。出口はその逆
  (符号は幾何として正しい A2 の符号に固定。K64 の追跡 track6n の tbshift も同じ符号)。ドリフト dz1 は運動学の項 (kin=True) を含む厳密なドリフト。K64 の `_quirk56` (k64 の行列の癖) は軌道には効かない。
  ずれ・ロールのある素子の中身は、誤差を 0 にした素子として ampnf の `build_ops` に渡す (素子 1 個ずつ)。
【ADTS】 `amplitude_tune_ext` は閉軌道 ≠ 0 のリング (DX・DY・DROTATE) の ν(J, δ) の級数 (ampnf.amplitude_tune の閉軌道版。ampnf の normal_form・tune_series を使う)。`resdrive.adts().rotate` が使う。
【限界】 閉軌道が大きい (数 mm 以上) と、D 次の打ち切りの誤差が増える (D を上げて確かめる)。δ の項 (p ≥ 1、h10002) は ampnf の `rdt_series` (閉軌道 0・δ=0 のとき)。
  運動学・フリンジ・偏向磁石の 4 次のハミルトニアンは TPSA の次数 D に入らないと消える (D=3 では 0)。3 次の写像 T3・J¹ にそれらが入る最小の次数は D=4 (既定)。
  K0≠0 の磁石は `ampnf.build_ops` の対応範囲に従う。K64 の四極の積分器の台形則・偏向磁石の K1 の 4 次以上は模型に入れていない。
"""
import copy
import dataclasses
import math
from types import SimpleNamespace

import numpy as np

from . import ampnf as A
from . import rdttrack as T

_DL = A.RAD * np.exp(2j * np.pi * np.arange(A.NTH) / A.NTH)


def err_args(el):
    """素子の tchge の引数 (dx, dy, dz, theta, dtheta, alg, phig)。誤差・ロールがなければ None"""
    tn = type(el).__name__
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if tn == 'Bend':
        from .lattice import _bend_chg_args
        return _bend_chg_args(el)
    if tn in ('Quadrupole', 'Sextupole', 'Octupole'):
        th = g('ROTATE') + g('ROTATION')                 # DROTATE は SAD の tquade・tmulte(6 極・8 極) に渡らない (効かない)
        if th == 0.0 and g('DX') == 0.0 and g('DY') == 0.0:
            return None
        return (g('DX'), g('DY'), 0.0, th, 0.0, 0.0, 0.0)
    if tn in ('Drift', 'Marker', 'Mark', 'Monitor', 'Moni', 'Cavity'):
        return None
    if any(g(k) for k in ('DX', 'DY', 'ROTATE', 'ROTATION', 'DROTATE')):
        raise NotImplementedError('ampnf_coupled: ずれ・ロールのある %s (%s) は未対応' % (tn, el.name))
    return None


def has_errors(lat):
    """DX・DY (BEND・QUAD・SEXT・OCT)・DZ・DROTATE (ANGLE≠0 の BEND だけ。SAD の tchge が使うのはそこだけ) のあるずれ・ロール誤差の素子の名前 (ROTATE・ROTATION だけの傾きは数えない)"""
    out = []
    for el in getattr(lat, 'elements', []):
        tn = type(el).__name__
        g = lambda k: float(getattr(el, k, 0.0) or 0.0)
        if tn not in ('Bend', 'Quadrupole', 'Sextupole', 'Octupole') or g('L') <= 0:
            continue
        if tn == 'Bend' and not any(g(k) for k in ('ANGLE', 'K1', 'K0')):
            continue
        if g('DX') or g('DY') or (tn == 'Bend' and g('ANGLE') != 0.0 and (g('DZ') or g('DROTATE'))):
            out.append(el.name)
    return out


def _lin_op(M4, d4, D, sd):
    Ms = np.zeros((D + 1, 4, 4), dtype=complex)
    Ms[0] = M4
    ds = np.zeros((D + 1, 4), dtype=complex)
    ds[0] = d4
    return ('lin', A.Seg.from_series(Ms, ds, sd))


def _drift_ops(al, D, sd, kin):
    """tdrife1 (長さ al、負も可) = 線形のドリフト + (kin=True なら) 運動学の厳密な流れ"""
    from .elements import Drift
    M, d, _k = A.element_linear(Drift(name='_d', L=al), _DL)
    ops = [('lin', A.Seg.from_series(A._series_from_values(M, D, _DL), A._series_from_values(d, D, _DL), sd))]
    if kin and al != 0.0:
        tabn = A.get_tab(5, max(D, 3))
        ops.append(('nlphi', A._flow_cached(('kdrift',), tabn, lambda: A._kin_h(tabn), al)))
    return ops


# --- tbgrote の TPSA (DROTATE) ---------------------------------------------------------------
def _tinv(tab, q):
    """TPSA の逆数 (定数項 ≠ 0)"""
    q0 = q[0]
    e = q / q0
    e = e.copy()
    e[0] = 0.0
    out = tab.const(1.0)
    term = tab.const(1.0)
    for _k in range(tab.D):
        term = -tab.mul(term, e)
        if not np.any(term):
            break
        out = out + term
    return out / q0


def _tsqrt1m(tab, a):
    """√(1 − a) (a は定数項 0 の TPSA)"""
    # √(1−a) = Σ binom(1/2, k) (−a)^k
    out = tab.const(1.0)
    pw = tab.const(1.0)
    for k in range(1, tab.D + 1):
        pw = tab.mul(pw, -a)
        if not np.any(pw):
            break
        out = out + A._binom_gen(0.5, k) * pw
    return out


def _tbgrote_frame(rr3, dxs, dzs):
    """tbgrote (tchge.f): dr = rr3 (3×3、回転行列 − I)、位置のずれ (dxs, dzs)。S = [x, px, y, py] の TPSA を書き換える関数 (S, tab, sd) → S"""
    dr = np.asarray(rr3, dtype=float)

    def f(S, tab, sd):
        x, px, y, py = S
        pr = tab.const(1.0) + sd * tab.var(4)
        inv2 = np.array([(-1.0) ** k * (k + 1) * sd ** k for k in range(tab.D + 1)])          # (1+δ)^−2 の δ̃ の級数
        a = tab.mul_dser(inv2, tab.mul(px, px) + tab.mul(py, py))
        pz = tab.mul(pr, _tsqrt1m(tab, a))
        xi = x + tab.const(dxs)
        yi = y
        zi = tab.const(dzs)
        xf = dr[0, 0] * xi + dr[0, 1] * yi + dr[0, 2] * zi + x
        yf = dr[1, 0] * xi + dr[1, 1] * yi + dr[1, 2] * zi + yi
        zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * zi
        pxf = (1.0 + dr[0, 0]) * px + dr[0, 1] * py + dr[0, 2] * pz
        pyf = dr[1, 0] * px + (1.0 + dr[1, 1]) * py + dr[1, 2] * pz
        pzf = dr[2, 0] * px + dr[2, 1] * py + (1.0 + dr[2, 2]) * pz
        ipz = _tinv(tab, pzf)
        q = tab.mul(zf, ipz)
        return [xf - tab.mul(pxf, q), pxf, yf - tab.mul(pyf, q), pyf]

    def fn(v, delta):
        """数値版 (v = (x, px, y, py)、δ 固定)"""
        x, px, y, py = v
        pr = 1.0 + delta
        a = (px / pr) ** 2 + (py / pr) ** 2
        pz = pr * math.sqrt(1.0 - a)
        xi, yi, zi = x + dxs, y, dzs
        xf = dr[0, 0] * xi + dr[0, 1] * yi + dr[0, 2] * zi + x
        yf = dr[1, 0] * xi + dr[1, 1] * yi + dr[1, 2] * zi + yi
        zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * zi
        pxf = (1.0 + dr[0, 0]) * px + dr[0, 1] * py + dr[0, 2] * pz
        pyf = dr[1, 0] * px + (1.0 + dr[1, 1]) * py + dr[1, 2] * pz
        pzf = dr[2, 0] * px + dr[2, 1] * py + (1.0 + dr[2, 2]) * pz
        return np.array([xf - pxf * zf / pzf, pxf, yf - pyf * zf / pzf, pyf])
    return ('frame', f, fn)


def _brote_rr(c0, s0, dtheta, ent):
    st, ct = math.sin(dtheta), math.cos(dtheta)
    dct = ct - 1.0
    rr = np.zeros((3, 3))
    rr[0, 0] = dct * c0 ** 2 + 0.0 * s0 ** 2
    rr[1, 1] = dct
    rr[2, 2] = dct * s0 ** 2
    if ent:
        rr[0, 1] = -c0 * st
        rr[0, 2] = c0 * dct * s0
        rr[1, 0] = c0 * st
        rr[1, 2] = s0 * st
        rr[2, 0] = s0 * c0 * dct
        rr[2, 1] = -s0 * st
    else:
        rr[1, 0] = -c0 * st
        rr[2, 0] = c0 * dct * s0
        rr[0, 1] = c0 * st
        rr[2, 1] = s0 * st
        rr[0, 2] = s0 * c0 * dct
        rr[1, 2] = -s0 * st
    return rr


def _chge_ops(args, ent, D, sd, kin):
    """tchge の入口 (ent=True) または出口の ops"""
    dx0, dy, dz, theta, dtheta, alg, phig = args
    th = theta if phig != 0.0 else theta + dtheta
    if th != 0.0:
        ct, st = math.cos(th), math.sin(th)
        dxr = ct * dx0 - st * dy
        ydis = st * dx0 + ct * dy
    else:
        ct, st = 1.0, 0.0
        dxr = dx0
        ydis = dy
    brote = (dxr != 0.0 or dz != 0.0 or dtheta != 0.0)
    if phig == 0.0:
        c0, s0, dl = 1.0, 0.0, alg
    else:
        c0, s0 = math.cos(phig), math.sin(phig)
        dl = alg * s0 / phig
    ops = []
    v2 = True       # 出口の (dx1, dz1) の符号は幾何として正しい v2 の符号に固定 (K64 の tchge (CALC・Emittance) は符号が逆だが、K64 の追跡 (track6n の tbshift) は v2 と同じ。k64 モードの追跡とも 1e-13 で一致)
    if ent:
        if th != 0.0:
            ops.append(_lin_op(A.rot_matrix(th), np.zeros(4), D, sd))
        if ydis != 0.0:
            ops.append(_lin_op(np.eye(4), np.array([0.0, 0.0, -ydis, 0.0]), D, sd))
        if brote:
            dx1 = dxr * c0 - dz * s0
            dz1 = dz * c0 + dxr * s0
            if dx1 != 0.0:
                ops.append(_lin_op(np.eye(4), np.array([-dx1, 0.0, 0.0, 0.0]), D, sd))
            if dz1 != 0.0:
                ops += _drift_ops(dz1, D, sd, kin)
            if dtheta != 0.0:
                ops.append(_tbgrote_frame(_brote_rr(c0, s0, dtheta, True), dl * s0, -dl * c0))
        return ops
    if brote:
        if dtheta != 0.0:
            ops.append(_tbgrote_frame(_brote_rr(c0, s0, dtheta, False), dl * s0, -dl * c0))
        if v2:
            dx1 = dxr * c0 - dz * s0
            dz1 = dz * c0 + dxr * s0
        else:
            dx1 = dxr * c0 + dz * s0
            dz1 = dz * c0 - dxr * s0
        if dz1 != 0.0:
            ops += _drift_ops(-dz1, D, sd, kin)
        if dx1 != 0.0:
            ops.append(_lin_op(np.eye(4), np.array([dx1, 0.0, 0.0, 0.0]), D, sd))
    if ydis != 0.0:
        ops.append(_lin_op(np.eye(4), np.array([0.0, 0.0, ydis, 0.0]), D, sd))
    if th != 0.0:
        ops.append(_lin_op(A.rot_matrix(-th), np.zeros(4), D, sd))
    return ops


def _quad_k64_ops(el, D, sd, fringe):
    """四極 (K1≠0) の K64 の積分器: nd = 1 + int(L |K1 L| / 0.2) 分割、(厳密な線形写像 ds) の間に運動学の流れ (両端 ds/2、内部 ds) = ∫ p⁴/8 ds の台形則。フリンジは本体の外 (ampnf と同じ)"""
    L, K1 = float(el.L), float(el.K1)
    nd = 1 + int(abs(L * (K1 * L) / 0.2))
    ds = L / nd
    tabn = A.get_tab(5, max(D, 3))
    kin = lambda d: ('nlphi', A._flow_cached(('kdrift',), tabn, lambda: A._kin_h(tabn), d))
    M, d_, _k = A.element_linear(dataclasses.replace(el, L=ds), _DL)
    lin = ('lin', A.Seg.from_series(A._series_from_values(M, D, _DL), A._series_from_values(d_, D, _DL), sd))
    ops = []
    if fringe:
        ops.append(('nlphi', A._flow_cached(('qfr', K1, 1), tabn, lambda: A._fringe_h(tabn, K1, +1.0), 1.0)))
    ops.append(kin(ds / 2))
    for i in range(nd):
        ops.append(lin)
        ops.append(kin(ds if i < nd - 1 else ds / 2))
    if fringe:
        ops.append(('nlphi', A._flow_cached(('qfr', K1, -1), tabn, lambda: A._fringe_h(tabn, K1, -1.0), 1.0)))
    return ops


def build_ops_ext(elements, D, sd, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True, nslice_quad=4, quad_rule='exact'):
    """quad_rule='k64': 四極の運動学の項を K64 の積分器の台形則 (分割数 1 + int(L |K1 L| / 0.2)、運動学の流れを両端に半分・内部の節点に 1 つ) で入れる (K64 の追跡との比較・原因の切り分け用。
    既定 'exact' は nslice_quad 分割の中点則 (nslice_quad → ∞ で厳密な積分)。
    素子列 → ops (ampnf の ('lin', Seg)・('nlphi', 流れ)・('kick', k2l, k3l) と ('frame', 関数))。ずれ・ロール (ROTATE・ROTATION・DX・DY・DZ・DROTATE) は tchge の入口・出口で包む"""
    ops = []
    for el in elements:
        if type(el).__name__ == 'Bend' and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0 and float(getattr(el, 'L', 0.0) or 0.0) > 0.0:
            # v2.4.0 G: K0≠0 の磁石の据え付け (DX・DY・ROTATE・DROTATE・CHI2) は、K64 の追跡と同じ式 (tbshift・tbrot・tent・texit) で ampnf_k0・ampnf._expand_rolls が内部で扱う
            # (tchge の経路 (CALC の式) は追跡と 1e-6 の差があるので使わない)
            sub, _w = A.build_ops(SimpleNamespace(elements=[el]), D, sd, nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=D, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
            ops += sub
            continue
        a = err_args(el)
        if a is None:
            e0 = el
            if float(getattr(el, 'DROTATE', 0.0) or 0.0):          # tchge を呼ばない場合 (ANGLE = K1 = K0 = 0 の BEND など) の DROTATE は SAD でも効かない
                e0 = copy.copy(el)
                e0.DROTATE = 0.0
        else:
            ops += _chge_ops(a, True, D, sd, kin)
            kw = {}
            for k in ('ROTATE', 'ROTATION', 'DX', 'DY'):
                if hasattr(el, k):
                    kw[k] = 0.0
            e0 = dataclasses.replace(el, **kw)
            for k in ('DROTATE', 'DZ'):          # 動的な属性 (dataclass のフィールドではない) は replace で落ちる
                if hasattr(e0, k):
                    setattr(e0, k, 0.0)
        if quad_rule == 'k64' and kin and type(e0).__name__ == 'Quadrupole' and float(e0.L) > 0 and float(e0.K1) != 0.0 and not float(getattr(e0, 'K1S', 0.0) or 0.0)                 and not float(getattr(e0, 'ROTATE', 0.0) or 0.0):
            ops += _quad_k64_ops(e0, D, sd, fringe)
        else:
            sub, _w = A.build_ops(SimpleNamespace(elements=[e0]), D, sd, nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, Dnl=D, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
            ops += sub
        if a is not None:
            ax = a
            if a[6] != 0.0:           # BEND の出口: tchge に alg − L、phig − ANGLE を渡す (elkern.k_bend と同じ)
                ax = a[:5] + (a[5] - float(el.L), a[6] - float(el.ANGLE))
            ops += _chge_ops(ax, False, D, sd, kin)
    out = []
    for o in ops:                               # 隣り合う線形の区間を 1 つにする
        if o[0] == 'lin' and out and out[-1][0] == 'lin':
            out[-1] = ('lin', out[-1][1].then(o[1]))
        else:
            out.append(o)
    return out


def propagate(ops, tab, S, sd, D):
    """ops を TPSA の S (4 成分 + 変数 4 = δ̃) に作用させる (ampnf.one_turn_series の本体と同じ。'frame' を追加)"""
    for o in ops:
        if o[0] == 'lin':
            sg = o[1]
            new = []
            for i in range(4):
                acc = tab.const(0.0)
                for j in range(4):
                    ser = sg.M[:, i, j]
                    if np.any(ser):
                        acc = acc + tab.mul_dser(ser, S[j])
                dser = sg.d[:, i]
                if np.any(dser):
                    acc = acc + tab.mul_dser(dser, tab.const(1.0))
                new.append(acc)
            S = new
        elif o[0] == 'nlphi':
            Q5 = S + [sd * tab.var(4)]
            S = A._compose(tab, o[1], Q5, D)
        elif o[0] == 'frame':
            S = o[1](S, tab, sd)
        elif o[0] == 'tmap':          # v2.4.0 G: K0≠0 の偏向磁石 (ampnf_k0。軌道のまわりの TPSA の写像)
            S = [np.asarray(a) for a in o[1](S, tab, sd)]
        else:
            _, k2l, k3l = o
            x, y = S[0], S[2]
            if k2l:
                x2, y2, xy = tab.mul(x, x), tab.mul(y, y), tab.mul(x, y)
                S[1] = S[1] - 0.5 * k2l * (x2 - y2)
                S[3] = S[3] + k2l * xy
            if k3l:
                x2, y2 = tab.mul(x, x), tab.mul(y, y)
                x3, y3, xy2, x2y = tab.mul(x2, x), tab.mul(y2, y), tab.mul(x, y2), tab.mul(x2, y)
                S[1] = S[1] - k3l / 6.0 * (x3 - 3 * xy2)
                S[3] = S[3] - k3l / 6.0 * (y3 - 3 * x2y)
    return S


def element_index(lat, ref):
    """素子の名前 (最初の一致) または添字 → 添字"""
    if isinstance(ref, str):
        idx = [i for i, e in enumerate(lat.elements) if e.name == ref]
        if not idx:
            raise ValueError('ampnf_coupled: %r の素子がない' % ref)
        return idx[0]
    return int(ref)


def numeric_track(ops, tabn, x4, delta, sd):
    """ops の数値の追跡 (模型の写像。δ 固定、x4 = (x, px, y, py))。閉軌道の探索用 (TPSA と同じ打ち切りの模型)"""
    v = np.array(x4, dtype=float)
    dt = delta / sd
    for o in ops:
        if o[0] == 'lin':
            sg = o[1]
            M = sum(sg.M[p] * dt ** p for p in range(sg.M.shape[0])).real
            d = sum(sg.d[p] * dt ** p for p in range(sg.d.shape[0])).real
            v = M @ v + d
        elif o[0] == 'nlphi':
            full = np.concatenate([v, [delta]])
            mono = np.prod(full[None, :] ** tabn.mons, axis=1)
            v = np.array([np.dot(o[1][j].real, mono) for j in range(4)])
        elif o[0] == 'frame':
            v = o[2](v, delta)
        elif o[0] == 'tmap':          # v2.4.0 G: K0≠0 の偏向磁石は追跡の素子の写像そのもの (数値)
            if len(o) > 3 and o[3]:
                raise NotImplementedError('ampnf_coupled.numeric_track: K0≠0 の磁石の線形だけの写像は数値の追跡が未対応')
            from . import ampnf_k0 as _k0
            v = _k0.numeric_map(o[2], delta)(np.asarray(v)[None, :])[0]
        else:
            _, k2l, k3l = o
            x, y = v[0], v[2]
            if k2l:
                v[1] -= 0.5 * k2l * (x * x - y * y)
                v[3] += k2l * x * y
            if k3l:
                v[1] -= k3l / 6.0 * (x ** 3 - 3 * x * y * y)
                v[3] -= k3l / 6.0 * (y ** 3 - 3 * x * x * y)
    return v


def closed_orbit_numeric(ops, tabn, delta, sd, tol=1e-14, maxit=40, h=1e-7):
    """模型の写像の閉軌道 (Newton 法、数値微分)"""
    z = np.zeros(4)
    for _ in range(maxit):
        f0 = numeric_track(ops, tabn, z, delta, sd) - z
        J = np.zeros((4, 4))
        for j in range(4):
            e = np.zeros(4)
            e[j] = h
            J[:, j] = (numeric_track(ops, tabn, z + e, delta, sd) - numeric_track(ops, tabn, z - e, delta, sd)) / (2 * h)
        dz = np.linalg.solve(J - np.eye(4), -f0)
        z = z + dz
        if np.abs(dz).max() < tol:
            break
    return z, float(np.abs(numeric_track(ops, tabn, z, delta, sd) - z).max())


def taylor_map(lat, delta=0.0, start=None, D=None, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True, nslice_quad=4, sd=1e-2, s=3e-4, ops=None, quad_rule='exact'):
    """模型の 1 周の写像 (リング先頭または start の素子の入口から) を、δ = delta の閉軌道のまわりの 3 次のテイラー展開 (正規化座標、結合したモード) にする。
    閉軌道は模型の数値の追跡の Newton 法で求め、TPSA はその閉軌道から出発する (閉軌道が数 mm でも各素子の展開は閉軌道のまわりで厳密)。
    D は TPSA の次数 (既定 4: 運動学・フリンジ・偏向磁石の 4 次のハミルトニアンの写像への寄与は D−1 次まで正しいので、3 次の写像 T3 には D ≥ 4 が要る。
    δ≠0 では T3 の δ^p の係数が p ≤ D−4 まで正しいので、既定は D = 6 (p ≤ 2、誤差 O(δ³))。
    戻り SimpleNamespace: R (4×4、回転)、T (4,4,4: 出力_i = Σ_ab T[i,a,b] z_a z_b)、T3 (list[dict {(ea,eb,ec,ed): 係数}])、orbit (物理座標の閉軌道)、A (z 物理 = orbit + A zn)、
    mu (モードの位相進み)、coupling (線形の結合の大きさ)、D、resid (閉軌道の残差)、delta"""
    if D is None:
        D = 4 if delta == 0.0 else 6
    D = max(int(D), 4)
    els = list(lat.elements)
    if start is not None:
        start = element_index(lat, start)
        els = els[start:] + els[:start]
    tab = A.get_tab(5, D)
    if ops is None:
        ops = build_ops_ext(els, D, sd, nslice=nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad, quad_rule=quad_rule)
    xs, resid = closed_orbit_numeric(ops, tab, delta, sd)
    S0 = [tab.const(xs[j]) + s * tab.var(j) for j in range(4)]
    Phi = propagate(ops, tab, S0, sd, D)
    dt = delta / sd
    mons, deg = tab.mons, tab.deg
    # δ̃ = δ/sd を代入した 1 周の線形写像 (物理座標のヤコビ行列): v_i の 1 次の係数 × δ̃^p
    Jm = np.zeros((4, 4))
    for m in range(tab.N):
        e = mons[m]
        if e[:4].sum() == 1:
            Jm[:, int(np.argmax(e[:4]))] += np.array([Phi[j][m] for j in range(4)]).real * dt ** int(e[4]) / s
    Am, Ai, mus, _tw, cpl = A.linear_modes(Jm)
    # 変数変換 v = Am u / s、δ̃ = δ/sd (定数)
    Q5 = [sum((Am[j, k] / s) * tab.var(k) for k in range(4) if Am[j, k] != 0.0) for j in range(4)] + [tab.const(dt)]
    Phi2 = A._compose(tab, Phi, Q5, D)
    arr = np.array([Phi2[j] for j in range(4)])
    const = float(np.abs(arr[:, 0] - xs).max())
    arr[:, 0] -= xs
    out = Ai @ arr
    R = np.zeros((4, 4))
    Tt = np.zeros((4, 4, 4))
    T3 = [dict() for _ in range(4)]
    for m in range(tab.N):
        if mons[m, 4] != 0 or deg[m] == 0 or deg[m] > 3:
            continue
        e = [int(x) for x in mons[m, :4]]
        col = out[:, m].real
        if deg[m] == 1:
            R[:, e.index(1)] = col
        elif deg[m] == 2:
            idx = [a for a in range(4) for _ in range(e[a])]
            a, b = idx
            if a == b:
                Tt[:, a, a] += col
            else:
                Tt[:, a, b] += 0.5 * col
                Tt[:, b, a] += 0.5 * col
        else:
            for i in range(4):
                if col[i] != 0.0:
                    T3[i][tuple(e)] = float(col[i])
    return SimpleNamespace(R=R, T=Tt, T3=T3, orbit=xs, A=Am, mu=mus, coupling=float(cpl), D=D, resid=resid, delta=delta, imag=float(np.abs(out.imag).max()), const=const)


def ref_phases(lat, ref, rm, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True, nslice_quad=4, sd=1e-2):
    """位相の基準 ref (素子の添字・名前。その素子の入口) でのモードの位相進み (ψ1, ψ2)。リング先頭の固有モード e_k (rm.A の列: e1 = A[:,0] + i A[:,1]、e2 = A[:,2] + i A[:,3]) を、
    閉軌道 (rm.orbit・rm.delta) のまわりのリング先頭から ref までの線形写像 M で運んだ M e_k の位相: ψ1 = arg (M e1)_x、ψ2 = arg (M e2)_y (結合がなければ Floquet の μx(ref)・μy(ref))"""
    idx = element_index(lat, ref)
    els = list(lat.elements)[:idx]
    D = rm.D
    tab = A.get_tab(5, D)
    ops = build_ops_ext(els, D, sd, nslice=nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    h = 1e-7
    M = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4)
        e[j] = h
        M[:, j] = (numeric_track(ops, tab, rm.orbit + e, rm.delta, sd) - numeric_track(ops, tab, rm.orbit - e, rm.delta, sd)) / (2 * h)
    Am = rm.A
    e1 = Am[:, 0] + 1j * Am[:, 1]
    e2 = Am[:, 2] + 1j * Am[:, 3]
    return float(np.angle((M @ e1)[0])), float(np.angle((M @ e2)[2]))


def rotate_to_ref(h, psi):
    """h の各項 (名前 'hjklmp') に、ref を位相の基準にする因子 exp(−i[(j−k)ψ1 + (l−m)ψ2]) を掛ける (resdrive の ref と同じ定義: h の絶対値は変わらず、角度の原点だけが動く)"""
    out = {}
    for nm, v in h.items():
        j, k, l, m = (int(c) for c in nm[1:5])
        out[nm] = v * np.exp(-1j * ((j - k) * psi[0] + (l - m) * psi[1]))
    return out


def one_turn_series_ext(lat, D=5, nslice=16, bend_cubic=True, nslice_bend=8, kin=True, fringe=True, nslice_quad=4, sd=1e-2, s=3e-4, ops=None):
    """`ampnf.one_turn_series` の閉軌道 ≠ 0 版 (DX・DY・DZ・DROTATE のあるリング): 閉軌道 (δ=0、模型の数値の追跡の Newton 法) から TPSA を出発させ、閉軌道での結合したモードの
    正規化座標の 1 周の写像 res (G = ζ の 4 成分、ampnf の normal_form・tune_series・initial_from_nf・nf_action がそのまま使える dict) を返す。res['orbit'] = 閉軌道 (物理座標)。
    D ≥ 4 (運動学・フリンジのハミルトニアンは写像より 1 次高い)"""
    D = max(int(D), 4)
    tab = A.get_tab(5, D)
    if ops is None:
        ops = build_ops_ext(list(lat.elements), D, sd, nslice=nslice, bend_cubic=bend_cubic, nslice_bend=nslice_bend, kin=kin, fringe=fringe, nslice_quad=nslice_quad)
    xs, resid = closed_orbit_numeric(ops, tab, 0.0, sd)
    h = 1e-7
    Jm = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4)
        e[j] = h
        Jm[:, j] = (numeric_track(ops, tab, xs + e, 0.0, sd) - numeric_track(ops, tab, xs - e, 0.0, sd)) / (2 * h)
    Am, Ai, mus, _tw, cpl = A.linear_modes(Jm)
    z = [tab.var(j) for j in range(5)]
    zn = [0.5 * s * (z[0] + z[1]), 0.5j * s * (z[0] - z[1]), 0.5 * s * (z[2] + z[3]), 0.5j * s * (z[2] - z[3])]
    S = [tab.const(xs[i]) + sum(Am[i, j] * zn[j] for j in range(4) if Am[i, j] != 0.0) for i in range(4)]
    S = propagate(ops, tab, S, sd, D)
    S = [S[i] - tab.const(xs[i]) for i in range(4)]
    xo, po, yo, qo = [sum(Ai[i, j] * S[j] for j in range(4) if Ai[i, j] != 0.0) for i in range(4)]
    G = [(xo - 1j * po) / s, (xo + 1j * po) / s, (yo - 1j * qo) / s, (yo + 1j * qo) / s]
    return dict(G=G, tab=tab, mux=mus[0], muy=mus[1], twiss0=(None,) * 4, Amat=Am, Ainv=Ai, coupling=float(cpl), warnings=[], s=s, sd=sd, D=D, orbit=xs, resid=resid)


def amplitude_tune_ext(lat, D=5, **kw):
    """`ampnf.amplitude_tune` の閉軌道 ≠ 0 版 (据え付けのずれ・DROTATE のあるリング)。ν(Jx, Jy, δ) のべき級数 (J は閉軌道のまわりの結合したモードのアクション)。戻り AmpTuneResult (.res['orbit'])"""
    res = one_turn_series_ext(lat, D=D, **kw)
    nf = A.normal_form(res)
    nux = A.tune_series(res, nf, 0)
    nuy = A.tune_series(res, nf, 1)
    nu0 = (res['mux'] / (2 * math.pi), res['muy'] / (2 * math.pi))
    return A.AmpTuneResult(nux, nuy, nu0, nf, res)


def tracked_adts_fit(lat, r, j=1e-7, mult=(1, 2, 3), nturn=1024, momentum_eV=3e9):
    """検証用: 追跡 + NAFF の J¹ の ADTS の安定なフィット (v2.4.0 H4)。`ampnf.tracked_adts` との違いは **フィットに定数項 ν(0) − ν0 を自由に入れる**こと。
    模型の ν0 と追跡 (K64 の積分器) の線形のチューンは 5e-7〜1e-6 違い、定数項なしのフィットではこの差が ADTS に ~ δν/J として入り、J が小さいほど悪化する
    (J = 3e-8 で 2〜6%)。定数項を入れると J = 3e-8〜1e-7 で 0.04〜0.15%、J = 3e-7 で 0.08〜0.3%。格子は J の倍数 mult × mult (J の 2 次式 + 定数でフィット)。
    戻り ({'xx','xy','yy'}, {'x': ν(0)−ν0, 'y': …})"""
    pts = [(a * j, b * j) for a in mult for b in mult]
    nus = np.array([A.tracked_tunes(lat, r, a, b, nturn=nturn, momentum_eV=momentum_eV) for a, b in pts])
    Am = np.array([[1.0, a, b, a * a, a * b, b * b] for a, b in pts])
    out, off = {}, {}
    for pl, nm in ((0, 'x'), (1, 'y')):
        y = nus[:, pl] - r.nu0[pl] % 1.0
        y = (y + 0.5) % 1.0 - 0.5
        co = np.linalg.lstsq(Am, y, rcond=None)[0]
        off[nm] = co[0]
        out[nm + 'x'], out[nm + 'y'] = co[1], co[2]
    return out, off


def unerrored_lattice(lat):
    """ロール (ROTATE・ROTATION) とずれ・ロール誤差 (DX・DY・DZ・DROTATE) を 0 にした新しい Lattice (BEND・QUAD・SEXT・OCT。元は変えない)"""
    from .lattice import Lattice
    out = []
    for el in lat.elements:
        if type(el).__name__ in ('Bend', 'Quadrupole', 'Sextupole', 'Octupole'):
            if any(float(getattr(el, k, 0.0) or 0.0) for k in ('ROTATE', 'ROTATION', 'DX', 'DY', 'DZ', 'DROTATE')):
                kw = {k: 0.0 for k in ('ROTATE', 'ROTATION', 'DX', 'DY') if hasattr(el, k)}
                e2 = dataclasses.replace(el, **kw)
                for k in ('DROTATE', 'DZ'):
                    if hasattr(e2, k):
                        setattr(e2, k, 0.0)
                el = e2
        out.append(el)
    return Lattice(out, name=getattr(lat, 'name', 'LATTICE'))


def rdt3(rm):
    """3 次の h_jklm (p=0 の 20 項。共役 (j,k,l,m) ↔ (k,j,m,l) も入る) {名前 'hjklm0': 複素数}"""
    return T.rdt_from_map(rm)


def rdt4(rm):
    """4 次の h_jklm (p=0 の 35 項)。戻り ({'hjklm0': 複素数}, info)"""
    return T.rdt4_from_map(rm)
