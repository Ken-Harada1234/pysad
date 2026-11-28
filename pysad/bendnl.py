# -*- coding: utf-8 -*-
"""Extract third-order bend nonlinearities as Lie generators from K64 tracking, for the second-order ADTS perturbation.
v2.3.0 R3: 偏向磁石の 3 次の非線形 (曲率の項・硬い端・複合機能の本体・E1/E2≠0 の端) のリー生成関数を、
K64 の追跡 (track6n) の写像の 2 次のテイラー展開から取り出して、`resdrive.adts()` の 2 次の摂動に渡す。K64 にない新機能。

考え方
  追跡の 1 要素の写像 M(z) = R z + T(z, z) + ... (z = (x, px, y, py)、δ=0、放射なし) の 2 次の部分 T を、(z と −z の平均) の Richardson 外挿
  (振幅 A と 2A) で数値的に取り出す (誤差 ~ A^6)。線形部分 R は中心差分 (Richardson)。
  薄い kick の生成関数 f (3 次の多項式) は、kick 写像 z → z + τ(z) (τ = 2 次の部分) から S∇f = τ (S = [[0,1],[-1,0]] の直和) で決まる:
  ∇f = −S τ、f = (1/3) z·∇f (同次式のオイラーの式)。
  偏向磁石は E1, E2 と本体 (長さ L を n 等分したスライス) に分ける:
    全体 = Edge(E2) ∘ B_n ∘ … ∘ B_1 ∘ Edge(E1)
    K64 の追跡は E=0 の端でも硬い端のフリンジ (Maxwell の項、H = ±(h/2) px y² の流れ) を持つ。長さ L/n・E=0 の偏向磁石 S = A_out ∘ B ∘ A_in
    (A_in: 入口 Δx = +y²/(2ρ)、Δpy = −px·y/ρ、A_out = A_in⁻¹) から B = A_in ∘ S ∘ A_out、端は Edge(E1) = A_in ∘ elem(ε, E1, 0)、
    Edge(E2) = elem(ε, 0, E2) ∘ A_out (ε は 1e-8 m の極薄い偏向磁石)。n 個の部分偏向磁石の連続は 1 つの偏向磁石と 1e-9 (相対) で一致する (確認済み)。
  スライスの生成関数は、スライスの中央の座標 (半分の線形写像 R_h を掛けた座標) で表す (対称な分割で O(ds²))。端の kick は線形部分 R_e の前に置く。
  角度・作用の変数の振幅 A_jklm (H = Σ A (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) e^{i[(j-k)φx+(l-m)φy]}) は、素子内の位置での w = (wx, wpx, wy, wpy)
  (w = √β e^{iμ}、wp = (i−α)/√β e^{iμ}。リングの Twiss から線形写像で伝播) から作る。
"""
import math
import dataclasses
import numpy as np

_S4 = np.zeros((4, 4))
_S4[0, 1] = 1.0
_S4[1, 0] = -1.0
_S4[2, 3] = 1.0
_S4[3, 2] = -1.0
_EPS_EDGE = 1e-8          # 端の取り出しに使う偏向磁石の長さ [m]
_CACHE = {}


# ---------------------------------------------------------------------------------------------
# 2 次のテイラー写像 (R, T): z → R z + T(z, z)。T[i] は対称な 4×4 行列 (τ_i = z^T T[i] z)
# ---------------------------------------------------------------------------------------------
class Taylor2:
    __slots__ = ('R', 'T')

    def __init__(self, R, T):
        self.R = np.asarray(R, dtype=float)
        self.T = np.asarray(T, dtype=float)

    def then(self, other):
        """Compose maps, applying this map first and the other map afterwards.
        self を先に、other を後に掛けた写像 other ∘ self"""
        R = other.R @ self.R
        T = np.einsum('ij,jab->iab', other.R, self.T) + np.einsum('ca,icd,db->iab', self.R, other.T, self.R)
        return Taylor2(R, T)

    def inverse(self):
        Ri = np.linalg.inv(self.R)
        T = -np.einsum('ij,jcd,ca,db->iab', Ri, self.T, Ri, Ri)
        return Taylor2(Ri, T)


def maxwell_in(rho):
    """Return the second-order part of the K64 hard-edge (E=0) entrance map: dx = y^2/(2 rho) and dpy = -px y/rho.
    K64 の硬い端 (E=0) の入口の写像の 2 次の部分: Δx = y²/(2ρ)、Δpy = −px·y/ρ (H = +(h/2) px y² の流れ、h = 1/ρ)。線形部分なし"""
    T = np.zeros((4, 4, 4))
    h = 1.0 / rho
    T[0, 2, 2] = 0.5 * h
    T[3, 1, 2] = T[3, 2, 1] = -0.5 * h
    return Taylor2(np.eye(4), T)


def extract_map(f, amp=5e-3, hlin=1e-4):
    """Extract the Taylor map to second order from a black-box function f of n x 4 coordinates.
    ブラックボックス f(Z (n × 4)) → (n × 4) から Taylor2 を取り出す。R: 中心差分の Richardson、T: ±z の平均の Richardson (振幅 amp と 2·amp)"""
    def central(h):
        e = np.eye(4) * h
        fp, fm = f(e), f(-e)
        return ((fp - fm) / (2 * h)).T
    R = (4.0 * central(hlin) - central(2 * hlin)) / 3.0
    idx = [(a, b) for a in range(4) for b in range(a, 4)]

    def tau_even(A):
        Z = []
        for a, b in idx:
            z = np.zeros(4)
            z[a] += A
            z[b] += A
            Z.append(z)
        Z = np.array(Z)
        E = np.eye(4) * A
        q = 0.5 * (f(Z) + f(-Z))
        qs = 0.5 * (f(E) + f(-E))
        T = np.zeros((4, 4, 4))
        for k, (a, b) in enumerate(idx):
            if a == b:
                T[:, a, a] = q[k] / (4.0 * A * A)
            else:
                c = (q[k] - qs[a] - qs[b]) / (A * A)
                T[:, a, b] = T[:, b, a] = 0.5 * c
        return T
    # c(A) = c + q A², c(2A) = c + 4 q A² → c = (4 c(A) − c(2A))/3
    T = (4.0 * tau_even(amp) - tau_even(2 * amp)) / 3.0
    return Taylor2(R, T)


def cubic_from_kick(Tk):
    """Return the generating function f (cubic homogeneous polynomial) of a kick z to z + tau(z) as a dictionary of coefficients.
    kick z → z + τ(z) (τ_i = z^T Tk[i] z) の生成関数 f (3 次の同次多項式) を {(ea, eb, ec, ed): 係数} で返す (変数 x, px, y, py)。
    ∇f = −S τ、f = (1/3) z·∇f"""
    f = {}
    for j in range(4):
        for i in range(4):
            s = _S4[j, i]
            if s == 0.0:
                continue
            for a in range(4):
                for b in range(4):
                    c = -s * Tk[i, a, b] / 3.0
                    if c == 0.0:
                        continue
                    e = [0, 0, 0, 0]
                    e[j] += 1
                    e[a] += 1
                    e[b] += 1
                    key = tuple(e)
                    f[key] = f.get(key, 0.0) + c
    return f


# ---------------------------------------------------------------------------------------------
# 偏向磁石の分解
# ---------------------------------------------------------------------------------------------
def _tracker_map(elements, momentum, delta=0.0):
    from .track6n import track_particles_fast

    def f(Z):
        c = np.zeros((len(Z), 6))
        c[:, :4] = Z
        c[:, 5] = delta
        o = np.array(track_particles_fast(elements, c, 1, momentum, 0.0, False, False, False))
        return o[:, :4]
    return f


def _mode_key():
    from . import compat
    # v2.3.0 R6: 追跡の写像が変わる設定は全部キーに入れる (U2: BEND の K1 の定義 'log'|'lin'、B161 の dk の修正)
    return (compat.v2(), compat.bend_k1_edge(), compat.bend_k1_def(), compat.bend_dk_fix())


def bend_decomposition(el, nslice=16, momentum=3e9):
    """Decompose a bend (L>0, ANGLE not 0, K0=0) into an entrance edge, nslice body slices, and an exit edge.
    偏向磁石 el (L>0、ANGLE≠0、K0=0) を 入口の端・nslice 個の本体・出口の端 に分ける。戻り dict:
      'R_e1'・'f_e1' (入口の端の線形写像と kick の生成関数)、'R_h' (本体の半スライスの線形写像)、'R_b' (スライス)、'f_b' (スライスの中央の生成関数)、
      'R_e2'・'f_e2'。生成関数は {(ex, epx, ey, epy): 係数} (スライス・端とも 1 つ分の積分値)"""
    dis = float(getattr(el, 'DISFRIN', 0.0) or 0.0)
    key = (el.L, el.ANGLE, el.K1, el.E1, el.E2, dis, nslice, _mode_key())
    if key in _CACHE:
        return _CACHE[key]
    L, ang = float(el.L), float(el.ANGLE)
    rho = L / ang
    n = nslice

    def mk(**kw):
        e2 = dataclasses.replace(el, name='_b', DX=0.0, DY=0.0, ROTATE=0.0, ROTATION=0.0, **kw)
        if dis != 0.0:
            e2.DISFRIN = dis             # dataclass のフィールドではない属性: replace では引き継がれない
        return e2
    ds = L / n
    # K64 の硬い端のフリンジ (Maxwell の項): K1=0 の偏向磁石 (tbendcore) では常に ON、K1≠0 の偏向磁石 (tbendi) では DISFRIN=0 のときだけ
    has_maxwell = (float(el.K1) == 0.0) or dis == 0.0
    A_in = maxwell_in(rho) if has_maxwell else Taylor2(np.eye(4), np.zeros((4, 4, 4)))
    A_out = A_in.inverse()
    # スライス (E=0): S = A_out ∘ B ∘ A_in → B = A_in ∘ S ∘ A_out (A_out を先に)
    S1 = extract_map(_tracker_map([mk(L=ds, ANGLE=ang / n, E1=0.0, E2=0.0)], momentum))
    B = A_out.then(S1).then(A_in)
    Sh = extract_map(_tracker_map([mk(L=ds / 2, ANGLE=ang / 2 / n, E1=0.0, E2=0.0)], momentum))
    Rh = Sh.R
    Rhi = np.linalg.inv(Rh)
    # 中央の座標での kick: z → z + Rh⁻¹ T_B(Rh⁻¹ z, Rh⁻¹ z)
    Tb = np.einsum('ij,jcd,ca,db->iab', Rhi, B.T, Rhi, Rhi)
    # v2.3.0 R4: 分散の生成 (z=0 の粒子の δ への 1 次の応答。半スライス d_h、スライス d_b)。端は δ の 1 次の分散を作らない
    dd = 1e-6
    z0 = np.zeros((1, 4))
    fh = _tracker_map([mk(L=ds / 2, ANGLE=ang / 2 / n, E1=0.0, E2=0.0)], momentum, dd)
    fhm = _tracker_map([mk(L=ds / 2, ANGLE=ang / 2 / n, E1=0.0, E2=0.0)], momentum, -dd)
    fs = _tracker_map([mk(L=ds, ANGLE=ang / n, E1=0.0, E2=0.0)], momentum, dd)
    fsm = _tracker_map([mk(L=ds, ANGLE=ang / n, E1=0.0, E2=0.0)], momentum, -dd)
    d_h = (fh(z0)[0] - fhm(z0)[0]) / (2 * dd)
    d_b = (fs(z0)[0] - fsm(z0)[0]) / (2 * dd)
    out = {'rho': rho, 'd_h': d_h, 'd_b': d_b, 'R_h': Rh, 'R_b': B.R, 'R_b_check': float(np.abs(B.R - Rh @ Rh).max()), 'f_b': cubic_from_kick(Tb), 'ds': ds, 'n': n}
    if el.E1 == 0.0:
        E1m = A_in
    else:
        E1m = extract_map(_tracker_map([mk(L=_EPS_EDGE, ANGLE=ang * _EPS_EDGE / L, E1=el.E1, E2=0.0)], momentum)).then(A_in)
    out['R_e1'] = E1m.R
    out['f_e1'] = cubic_from_kick(np.einsum('ij,jab->iab', np.linalg.inv(E1m.R), E1m.T))
    if el.E2 == 0.0:
        E2m = A_out
    else:
        E2m = A_out.then(extract_map(_tracker_map([mk(L=_EPS_EDGE, ANGLE=ang * _EPS_EDGE / L, E1=0.0, E2=el.E2)], momentum)))
    out['R_e2'] = E2m.R
    out['f_e2'] = cubic_from_kick(np.einsum('ij,jab->iab', np.linalg.inv(E2m.R), E2m.T))
    _CACHE[key] = out
    return out


# ---------------------------------------------------------------------------------------------
# 生成関数 → 調和ごとの振幅 A_jklm
# ---------------------------------------------------------------------------------------------
def _pmul(p, q):
    r = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0] + k2[0], k1[1] + k2[1])
            r[k] = r[k] + c1 * c2 if k in r else c1 * c2
    return r


def harmonic_amplitudes(f, w, monos):
    """Compute the angle-action harmonic amplitudes A_jklm of a cubic generating function f at given positions.
    生成関数 f ({(ex, epx, ey, epy): 係数}、座標 (x, px, y, py) の 3 次式) の、位置 (w = (wx, wpx, wy, wpy)、複素 (4, N) または (4,)) での
    角度・作用の単項式 z^j z̄^k v^l v̄^m (z = √(2Jx) e^{iφx}、v = √(2Jy) e^{iφy}) の係数 A_jklm (monos の順)。
    x = (wx z + 共役)/2、px = (wpx z + 共役)/2、y = (wy v + 共役)/2、py = (wpy v + 共役)/2。戻り: (len(monos), N) の複素配列"""
    w = np.asarray(w, dtype=complex)
    if w.ndim == 1:
        w = w[:, None]
    N = w.shape[1]
    lin = [{(1, 0): 0.5 * w[i], (0, 1): 0.5 * np.conj(w[i])} for i in range(4)]
    out = {}
    for e, c in f.items():
        px = {(0, 0): np.ones(N, dtype=complex)}
        for i in (0, 1):
            for _ in range(e[i]):
                px = _pmul(px, lin[i])
        py = {(0, 0): np.ones(N, dtype=complex)}
        for i in (2, 3):
            for _ in range(e[i]):
                py = _pmul(py, lin[i])
        for (j, k), cx in px.items():
            for (l, m), cy in py.items():
                key = (j, k, l, m)
                out[key] = out[key] + c * cx * cy if key in out else c * cx * cy
    res = np.zeros((len(monos), N), dtype=complex)
    for iq, q in enumerate(monos):
        key = (q[0], q[1], q[2], q[3])
        if key in out:
            res[iq] = out[key]
    return res


def w_from_row(r):
    """Return the complex vector w = (wx, wpx, wy, wpy) from a Twiss row.
    Twiss の行から複素の w = (wx, wpx, wy, wpy) (x = √(2J) Re[wx e^{iφ}]、px = √(2J) Re[wpx e^{iφ}])"""
    bx, ax, mx = float(r.betx), float(r.alfx), float(r.mux)
    by, ay, my = float(r.bety), float(r.alfy), float(r.muy)
    ex, ey = np.exp(1j * mx), np.exp(1j * my)
    return np.array([math.sqrt(bx) * ex, (1j - ax) / math.sqrt(bx) * ex, math.sqrt(by) * ey, (1j - ay) / math.sqrt(by) * ey])


def bend_atoms(el, row_in, row_out, monos, nslice=16, momentum=3e9, check=True):
    """Compute harmonic amplitudes of the thin kicks (entrance edge, slices, exit edge) of a bend from the Twiss row at its entrance.
    偏向磁石 el の入口 (row_in) の Twiss から、入口の端の kick、nslice 個のスライス、出口の端の kick を、線形写像で伝播した位置の w で並べた
    薄い kick の列 (Natom = nslice + 2) の調和の振幅。戻り (a (len(monos), Natom)、info dict)。a は 2^{3/2} を含む
    (resdrive.adts の a_all と同じ規格化: (2J)^{3/2} → J^{3/2} 2^{3/2})。check=True: 出口の w が row_out と一致するか (info['w_err'])"""
    d = bend_decomposition(el, nslice, momentum)
    n = d['n']
    w0 = w_from_row(row_in)
    amps = [harmonic_amplitudes(d['f_e1'], w0, monos)]
    w = d['R_e1'] @ w0
    wc = np.zeros((4, n), dtype=complex)
    for k in range(n):
        wc[:, k] = d['R_h'] @ w
        w = d['R_b'] @ w
    amps.append(harmonic_amplitudes(d['f_b'], wc, monos))
    amps.append(harmonic_amplitudes(d['f_e2'], w, monos))
    w = d['R_e2'] @ w
    info = {'w_err': 0.0, 'R_b_check': d['R_b_check']}
    if check:
        w1 = w_from_row(row_out)
        info['w_err'] = float(np.abs(w - w1).max() / np.abs(w1).max())
    a = np.concatenate(amps, axis=1) * 2.0 ** 1.5
    return a, info


# ---------------------------------------------------------------------------------------------
# v2.3.0 R4: 共鳴駆動項 h_jklmp への寄与 (δ の項を含む)。偏向磁石の 3 次の生成関数 f を、分散による軌道のずれ z → z + e δ で展開する
# ---------------------------------------------------------------------------------------------
def harmonic_amplitudes_shift(f, w, e, monos):
    """Compute the amplitudes of the angle-action monomials for the delta^p coefficients of f(z + e delta).
    f(z + e δ) の δ^p の係数の、角度・作用の単項式 z^j z̄^k v^l v̄^m の振幅 (monos の各要素は (j, k, l, m, p))。
    f: {(ex, epx, ey, epy): 係数} (任意の次数)、w: 複素 (4, N)、e: 実 (4, N) (各位置の分散ベクトル (η, η', η_y, η_y'))。
    x = (wx z + 共役)/2 + ex δ など。p = 0 は harmonic_amplitudes と同じ。戻り (len(monos), N) の複素配列"""
    w = np.asarray(w, dtype=complex)
    e = np.asarray(e, dtype=float)
    if w.ndim == 1:
        w = w[:, None]
    if e.ndim == 1:
        e = np.repeat(e[:, None], w.shape[1], axis=1)
    N = w.shape[1]
    lin = [{(1, 0, 0): 0.5 * w[i], (0, 1, 0): 0.5 * np.conj(w[i]), (0, 0, 1): e[i].astype(complex)} for i in range(4)]

    def mul(p, q):
        r = {}
        for k1, c1 in p.items():
            for k2, c2 in q.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                r[k] = r[k] + c1 * c2 if k in r else c1 * c2
        return r
    out = {}
    for ex_, c in f.items():
        px = {(0, 0, 0): np.ones(N, dtype=complex)}
        for i in (0, 1):
            for _ in range(ex_[i]):
                px = mul(px, lin[i])
        py = {(0, 0, 0): np.ones(N, dtype=complex)}
        for i in (2, 3):
            for _ in range(ex_[i]):
                py = mul(py, lin[i])
        for (j, k, t1), cx in px.items():
            for (l, m, t2), cy in py.items():
                key = (j, k, l, m, t1 + t2)
                out[key] = out[key] + c * cx * cy if key in out else c * cx * cy
    res = np.zeros((len(monos), N), dtype=complex)
    for iq, q in enumerate(monos):
        key = tuple(int(v) for v in q)
        if key in out:
            res[iq] = out[key]
    return res


def bend_atoms_rdt(el, row_in, row_out, monos, nslice=16, momentum=3e9, mu_ref=(0.0, 0.0)):
    """Return the contributions of a bend to the resonance driving terms h_(j,k,l,m,p), per thin kick of edge and slices.
    偏向磁石 el の、共鳴駆動項 h_(j,k,l,m,p) (monos の各要素 (j, k, l, m, p)、j+k+l+m+p = 3 のものだけが 0 でない) への寄与を、
    入口の端・nslice 個のスライス・出口の端の薄い kick ごとに返す (a: (len(monos), nslice + 2)、2^{3/2} を含まない = h の定義そのまま)。
    w は Twiss の行から線形写像で伝播 (位相は mu_ref を引く)、分散ベクトル e = (η, η', 0, 0) は入口の行から線形写像と生成された分散 (d_h・d_b) で伝播する。
    info: w_err (出口の w と row_out の差)、e_err (出口の分散と row_out の差。相対)"""
    d = bend_decomposition(el, nslice, momentum)
    n = d['n']
    w0 = w_from_row(row_in) * np.array([np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[1]), np.exp(-1j * mu_ref[1])])
    e0 = np.array([float(row_in.etx), float(row_in.etpx), 0.0, 0.0])
    ph = np.array([np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[1]), np.exp(-1j * mu_ref[1])])
    amps = [harmonic_amplitudes_shift(d['f_e1'], w0, e0, monos)]
    w = d['R_e1'] @ w0
    e = d['R_e1'] @ e0
    wc = np.zeros((4, n), dtype=complex)
    ec = np.zeros((4, n))
    for k in range(n):
        wc[:, k] = d['R_h'] @ w
        ec[:, k] = d['R_h'] @ e + d['d_h']
        w = d['R_b'] @ w
        e = d['R_b'] @ e + d['d_b']
    amps.append(harmonic_amplitudes_shift(d['f_b'], wc, ec, monos))
    amps.append(harmonic_amplitudes_shift(d['f_e2'], w, e, monos))
    w = d['R_e2'] @ w
    e = d['R_e2'] @ e
    w1 = w_from_row(row_out) * ph
    e1 = np.array([float(row_out.etx), float(row_out.etpx), 0.0, 0.0])
    info = {'w_err': float(np.abs(w - w1).max() / np.abs(w1).max()),
            'e_err': dispersion_err(e, e1, e0, float(el.ANGLE))}
    return np.concatenate(amps, axis=1), info


def dispersion_err(e, e1, e0, angle):
    """Return the relative error of the propagated dispersion vector `e` against the Twiss value `e1` at a bend exit.
    偏向磁石の出口の分散ベクトル e (線形写像と生成分散で伝播) と Twiss の e1 の差の相対値。分母は出口 e1・入口 e0 の大きさと偏向角 (生成される分散 η' の尺度) の最大。
    出口の分散が 0 (アクロマット端。DDBA の鏡像の側の磁石など、分散が 0 に戻る磁石) でも、分母が 1e-12 に落ちて数値の残差 (絶対値 1e-12 m 程度) が相対 1 に見えないようにする"""
    scale = max(float(np.abs(e1).max()), float(np.abs(e0).max()), abs(float(angle)), 1e-12)
    return float(np.abs(e - e1).max() / scale)


# ---------------------------------------------------------------------------------------------
# v2.3.0 R6: 閉軌道・δ のまわりの分解 (K0≠0 の偏向磁石 = 挿入光源、δ≠0)
#   原子 (薄い kick) の生成関数を、各原子の位置の閉軌道 z_k (δ のエネルギー) のまわりのテイラー展開から取り出す。
#   スライスは E=0 の部分偏向磁石 S (長さ L/n、角度・K0 も 1/n) の写像そのもの (硬い端の対は内部に入っていて、S = A⁻¹ B A = exp(ds·Ĥ) の形。
#   Ĥ は定数のハミルトニアンで、スライスの kick の和は ∫Ĥ3 ds の中点則 = 端の項 A を別に置く R3 の分解と同じ積分)。
#   端の薄い kick は、E1≠0 の極薄い偏向磁石 (ε = 1e-8 m) の写像 (その内部の硬い端の対と合わせて、全体の写像と厳密に一致する)。A_in の解析式が要らない。
# ---------------------------------------------------------------------------------------------
_ORB_CACHE = {}


def jac_at(f, z, h=1e-4):
    """Compute the Jacobian matrix of f at z by Richardson-extrapolated central differences.
    f(Z (n × 4)) → (n × 4) の z でのヤコビ行列 (中心差分の Richardson)"""
    z = np.asarray(z, dtype=float)

    def cen(hh):
        e = np.eye(4) * hh
        return ((f(z + e) - f(z - e)) / (2 * hh)).T
    return (4.0 * cen(h) - cen(2 * h)) / 3.0


def _sub_element(el, frac, E1=0.0, E2=0.0, L=None):
    """el の長さ・角度・K0 を frac 倍した E1・E2 つきの複製 (据え付けのずれ・回転なし)"""
    dis = float(getattr(el, 'DISFRIN', 0.0) or 0.0)
    if not hasattr(el, 'ANGLE'):                 # 四極・ドリフトなど (要素全体の写像を取り出すとき)
        return _plain(el)
    kw = dict(name='_b', DX=0.0, DY=0.0, ROTATE=0.0, ROTATION=0.0, L=float(el.L) * frac if L is None else L, ANGLE=float(el.ANGLE) * frac, E1=E1, E2=E2)
    if hasattr(el, 'K0'):
        kw['K0'] = float(el.K0) * frac
    e2 = dataclasses.replace(el, **kw)
    if dis != 0.0:
        e2.DISFRIN = dis
    return e2


def _taylor_about(f, z):
    fz = f(np.asarray(z, dtype=float)[None, :])[0]
    return extract_map(lambda Z: f(Z + z) - fz), fz


def orbit_decomposition(el, z_in, delta=0.0, nslice=16, momentum=3e9):
    """Decompose a bend (L>0, K0 may be nonzero) into entrance edge, slices, and exit edge about the entrance closed orbit and delta.
    偏向磁石 el (L>0。K0≠0 でもよい) の 入口の端・nslice 個の本体・出口の端 への分解を、入口の閉軌道 z_in (物理座標 (x, px, y, py)) と δ のまわりで作る。
    戻り dict: 'units' = [(種類 'e1'|'b'|'e2', R_pre, f, R_post)] (原子 1 個 = 線形写像 R_pre → kick (生成関数 f、R_pre の後の座標で) → 線形写像 R_post。端は R_pre = 1)、
    'z_out' (出口の閉軌道)、'z_mid' (原子の位置の閉軌道)、'n'"""
    z_in = np.asarray(z_in, dtype=float)
    dis = float(getattr(el, 'DISFRIN', 0.0) or 0.0)
    L = float(el.L)
    n = int(nslice)
    key = (type(el).__name__, L, float(el.ANGLE), float(getattr(el, 'K0', 0.0) or 0.0), float(el.K1), float(el.E1), float(el.E2), dis, n, float(delta),
           tuple(np.round(z_in, 13)), momentum, _mode_key())
    if key in _ORB_CACHE:
        return _ORB_CACHE[key]
    f_s = _tracker_map([_sub_element(el, 1.0 / n)], momentum, delta)
    f_h = _tracker_map([_sub_element(el, 0.5 / n)], momentum, delta)
    units = []
    zmid = []
    z = z_in.copy()
    if float(el.E1) != 0.0:
        e1 = _sub_element(el, 0.0, E1=float(el.E1), L=_EPS_EDGE)
        e1.ANGLE, e1.K0 = float(el.ANGLE) * _EPS_EDGE / L, float(getattr(el, 'K0', 0.0) or 0.0) * _EPS_EDGE / L
        f_e = _tracker_map([e1], momentum, delta)
        m, z1 = _taylor_about(f_e, z)
        units.append(('e1', np.eye(4), cubic_from_kick(np.einsum('ij,jab->iab', np.linalg.inv(m.R), m.T)), m.R))
        zmid.append(z.copy())
        z = z1
    for _ in range(n):
        Rh1 = jac_at(f_h, z)
        zm = f_h(z[None, :])[0]
        Rh2 = jac_at(f_h, zm)
        m, z1 = _taylor_about(f_s, z)
        Tb = np.einsum('ij,jcd,ca,db->iab', np.linalg.inv(Rh2), m.T, np.linalg.inv(Rh1), np.linalg.inv(Rh1))
        units.append(('b', Rh1, cubic_from_kick(Tb), Rh2))
        zmid.append(zm.copy())
        z = z1
    if float(el.E2) != 0.0:
        e2 = _sub_element(el, 0.0, E2=float(el.E2), L=_EPS_EDGE)
        e2.ANGLE, e2.K0 = float(el.ANGLE) * _EPS_EDGE / L, float(getattr(el, 'K0', 0.0) or 0.0) * _EPS_EDGE / L
        f_e = _tracker_map([e2], momentum, delta)
        m, z1 = _taylor_about(f_e, z)
        units.append(('e2', np.eye(4), cubic_from_kick(np.einsum('ij,jab->iab', np.linalg.inv(m.R), m.T)), m.R))
        zmid.append(z.copy())
        z = z1
    out = {'units': units, 'z_out': z, 'z_mid': np.array(zmid), 'n': n}
    _ORB_CACHE[key] = out
    return out


def row_orbit(r):
    """Return the closed orbit (x, px, y, py) of a Twiss row.
    Twiss の行の閉軌道 (x, px, y, py)"""
    return np.array([float(r.x_cod), float(r.px_cod), float(r.y_cod), float(r.py_cod)])


def bend_atoms_orbit(el, row_in, row_out, monos, nslice=16, momentum=3e9, delta=0.0, check=True):
    """Closed-orbit and delta version of `bend_atoms`: harmonic amplitudes of the atoms extracted about the entrance closed orbit (R6).
    `bend_atoms` の閉軌道・δ 版 (R6)。入口の閉軌道 (row_in)・δ のまわりで取り出した原子 (端・スライス・端) の調和の振幅 (len(monos), Natom) (2^{3/2} を含む)。
    info: w_err (出口の w と row_out)、z_err (出口の閉軌道と row_out の差、m・rad)"""
    d = orbit_decomposition(el, row_orbit(row_in), delta, nslice, momentum)
    w = w_from_row(row_in)
    amps = []
    for kind, Rpre, f, Rpost in d['units']:
        wc = Rpre @ w
        amps.append(harmonic_amplitudes(f, wc, monos))
        w = Rpost @ wc
    info = {'w_err': 0.0, 'z_err': 0.0}
    if check:
        w1 = w_from_row(row_out)
        info['w_err'] = float(np.abs(w - w1).max() / np.abs(w1).max())
        info['z_err'] = float(np.abs(d['z_out'] - row_orbit(row_out)).max())
    a = np.concatenate(amps, axis=1) * 2.0 ** 1.5
    return a, info


# 閉軌道がこれ以下 (CALC の呼び出しの履歴に由来する 1e-10 程度の雑音) なら原点のまわりの取り出し (R3) を使う。orbit 版と R3 版は K1≠0 の偏向磁石の内部の対の項の分け方が違い、
# fourth='r3' では合計が 1〜2% 変わる (fourth='r6' は取り出した 4 次の項で補うので同じ)。雑音で経路が切り替わらないようにする
ORBIT_TOL = 1e-8


def needs_orbit(el, row_in, delta=0.0):
    """Return whether a bend needs extraction about the closed orbit (K0 not 0, delta not 0, or displaced entrance orbit).
    K0≠0 (挿入光源)・δ≠0・入口の閉軌道がずれている偏向磁石は、閉軌道のまわりの取り出し (bend_atoms_orbit) が要る"""
    return (float(getattr(el, 'K0', 0.0) or 0.0) != 0.0 or float(delta) != 0.0 or float(np.abs(row_orbit(row_in)).max()) > ORBIT_TOL)


def bend_atoms_auto(el, row_in, row_out, monos, nslice=16, momentum=3e9, delta=0.0, check=True, policy='auto'):
    """Choose between `bend_atoms` (R3) and `bend_atoms_orbit` (R6), which return values of the same form.
    `bend_atoms` (閉軌道 0、δ=0、K0=0。R3) と `bend_atoms_orbit` (R6) を選ぶ。どちらも同じ形の戻り値 (a (len(monos), Natom)、info)。
    閉軌道 0・δ=0・K0=0 では両者は 1e-14 で一致 (試験)。policy='k0' は K0≠0 の磁石だけ bend_atoms_orbit (δ・閉軌道のずれは原点のまわりの R3 の分解のまま)"""
    # policy='k0': K0≠0 の磁石だけ閉軌道のまわり (fourth='r3' 用。内部の対の項の分け方を原点のまわりの R3 の分解と合わせる)
    if (float(getattr(el, 'K0', 0.0) or 0.0) != 0.0) if policy == 'k0' else needs_orbit(el, row_in, delta):
        return bend_atoms_orbit(el, row_in, row_out, monos, nslice=nslice, momentum=momentum, delta=delta, check=check)
    return bend_atoms(el, row_in, row_out, monos, nslice=nslice, momentum=momentum, check=check)


def _row_disp(r):
    return np.array([float(r.etx), float(r.etpx), float(getattr(r, 'ety', 0.0) or 0.0), float(getattr(r, 'etpy', 0.0) or 0.0)])


def element_delta_terms(el, row_in, row_out, momentum=3e9, d=1e-4, mu_ref=(0.0, 0.0), slope=True):
    """Extract directly from element tracking the generating functions of the delta terms (p=1 and p=2) of a bend with K0 not 0 (R6).
    K0≠0 の偏向磁石 (挿入光源) 1 台の、δ の項 (p=1: δ·(ζ の 2 次形式)、p=2: δ²·(ζ の 1 次)) の生成関数を、要素の追跡の写像から直接取り出す (R6)。
    展開の中心は 1 次の分散の軌道 z(δ) = z_co + e δ (入口 row_in・出口 row_out の Twiss の行の閉軌道と分散)。
      p=1: ヤコビ行列の δ 微分 X = ∂J(z(δ); δ)/∂δ (中心差分の Richardson、δ = ±d・±2d)。R⁻¹X = S Q (R = δ=0 のヤコビ行列) から g2 = ½ ζᵀQζ。
           **自然なクロマティックな項 (K0 の磁場の 1/(1+δ)・端の収束・傾いた軌道) と分散による軌道のずれの分の全体** (K0≠0 の偏向磁石は自然な項の解析式に入っていない)。
      p=2: 展開の中心が 2 次の閉軌道からずれる分の kick c₂ = ½ ∂²[M(z(δ); δ) − z_out(δ)]/∂δ² (偶部分の Richardson)。u = R⁻¹c₂ = S∇g から g1 = ∇g·ζ。
    戻り: 入口の w での振幅 dict {(j, k, l, m, p): 複素数} (p=1: 2 次形式の単項式、p=2: 1 次の単項式。monos の取り方は harmonic_amplitudes と同じ)"""
    f0 = _tracker_map([_sub_element(el, 1.0)], momentum, 0.0)
    z0, z1 = row_orbit(row_in), row_orbit(row_out)
    e0, e1 = _row_disp(row_in), _row_disp(row_out)
    R0 = jac_at(f0, z0)

    def tmat(dd):
        # スロープ座標 (x, x' = px/(1+δ)) への変換 z_can = T z_slope。追跡は正準運動量 px (p0 で規格化) だが、自然クロマティシティの解析式 (K1/(1+δ)) と
        # 分けて足すには、素子ごとにスロープ座標の δ 依存 (ドリフトの 1/(1+δ) などを含まない) で測る必要がある
        return np.diag([1.0, 1.0 + dd, 1.0, 1.0 + dd]) if slope else np.eye(4)

    def at(dd):
        f = _tracker_map([_sub_element(el, 1.0)], momentum, dd)
        zi = z0 + e0 * dd
        return np.linalg.solve(tmat(dd), jac_at(f, zi) @ tmat(dd)), f(zi[None, :])[0] - (z1 + e1 * dd)
    X = {}
    C = {}
    for dd in (d, 2 * d):
        Jp, cp = at(dd)
        Jm, cm = at(-dd)
        X[dd] = (Jp - Jm) / (2 * dd)
        C[dd] = 0.5 * (cp + cm) / (dd * dd)
    Xr = (4.0 * X[d] - X[2 * d]) / 3.0
    c2 = (4.0 * C[d] - C[2 * d]) / 3.0
    Q = -_S4 @ np.linalg.solve(R0, Xr)
    f2 = {}
    for a in range(4):
        for b in range(4):
            e = [0] * 4
            e[a] += 1
            e[b] += 1
            f2[tuple(e)] = f2.get(tuple(e), 0.0) + 0.5 * Q[a, b]
    u = np.linalg.solve(R0, c2)
    gradg = -_S4 @ u
    f1 = {(1, 0, 0, 0): gradg[0], (0, 1, 0, 0): gradg[1], (0, 0, 1, 0): gradg[2], (0, 0, 0, 1): gradg[3]}
    ph = np.array([np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[0]), np.exp(-1j * mu_ref[1]), np.exp(-1j * mu_ref[1])])
    w = w_from_row(row_in) * ph
    out = {}
    m2 = [(j, k, l, m) for j in range(3) for k in range(3 - j) for l in range(3 - j - k) for m in [2 - j - k - l]]
    a2 = harmonic_amplitudes(f2, w, m2)[:, 0]
    for q, v in zip(m2, a2):
        out[q + (1,)] = complex(v)
    m1 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    a1 = harmonic_amplitudes(f1, w, m1)[:, 0]
    for q, v in zip(m1, a1):
        out[q + (2,)] = complex(v)
    return out


# ---------------------------------------------------------------------------------------------
# v2.3.0 R6: δ≠0・閉軌道がずれた磁石 (四極・ドリフト・K1 の MULT) の 3 次の項。運動学的な項 p⁴/8 や四極端部のフリンジを、傾いた・ずれた軌道のまわりで展開すると
# 3 次の項が出る (px0 (px² + py²) px/2 など。px0 = η'δ)。追跡の写像のテイラー展開から取り出す (6 極・偏向磁石以外のうち、軌道が 0 でないものだけ)。
# ---------------------------------------------------------------------------------------------
def needs_other(el, row_in, delta=0.0):
    tn = type(el).__name__
    if tn not in ('Drift', 'Quadrupole', 'Multipole', 'Sextupole', 'Octupole') or float(getattr(el, 'L', 0.0) or 0.0) <= 0.0:
        return False
    if tn == 'Multipole' and any(int(n) != 1 for n, v in (getattr(el, 'KnL', {}) or {}).items() if float(v) != 0.0):
        return False          # K2 以上を持つ MULT は多極 (collect_nonlinear が扱う)。K1 だけの MULT を四極として扱う
    return float(delta) != 0.0 or float(np.abs(row_orbit(row_in)).max()) > 1e-12


def _plain(el):
    kw = {'name': '_b'}
    for k in ('DX', 'DY', 'ROTATE', 'ROTATION'):
        if hasattr(el, k):
            kw[k] = 0.0
    if type(el).__name__ == 'Sextupole':          # 6 極の kick は解析式 (別の項)。ここはドリフトの運動学的な項だけ
        kw['K2'] = 0.0
    if type(el).__name__ == 'Octupole':
        kw['K3'] = 0.0
    return dataclasses.replace(el, **kw)


def other_atoms(el, row_in, row_out, monos, delta=0.0, nslice_drift=2, momentum=3e9, check=True):
    """Return the harmonic amplitudes of atoms about the closed orbit and delta for quadrupoles, drifts, and K1-only MULT.
    四極・ドリフト・K1 だけの MULT の閉軌道・δ のまわりの原子 (ドリフトは nslice_drift 個の等分、それ以外は要素全体で 1 個。四極は端部のフリンジが
    要素の両端にあるので分割しない)。戻り (a (len(monos), Natom) (2^{3/2} を含む)、info)"""
    z = row_orbit(row_in)
    L = float(el.L)
    n = int(nslice_drift) if type(el).__name__ in ('Drift', 'Sextupole', 'Octupole') else 1
    key = ('other', type(el).__name__, tuple((k, repr(v)) for k, v in sorted(vars(el).items()) if k != 'name'), n, float(delta), tuple(np.round(z, 13)), momentum, _mode_key())
    if key in _ORB_CACHE:
        units, z_out = _ORB_CACHE[key]
    else:
        units = []
        sub = _plain(el) if n == 1 else dataclasses.replace(_plain(el), L=L / n)
        f = _tracker_map([sub], momentum, delta)
        zc = z.copy()
        for _ in range(n):
            m, z1 = _taylor_about(f, zc)
            units.append((np.eye(4), cubic_from_kick(np.einsum('ij,jab->iab', np.linalg.inv(m.R), m.T)), m.R))
            zc = z1
        z_out = zc
        _ORB_CACHE[key] = (units, z_out)
    w = w_from_row(row_in)
    amps = []
    for Rpre, fk, Rpost in units:
        amps.append(harmonic_amplitudes(fk, Rpre @ w, monos))
        w = Rpost @ (Rpre @ w)
    info = {'w_err': 0.0, 'z_err': 0.0}
    if check:
        w1 = w_from_row(row_out)
        info['w_err'] = float(np.abs(w - w1).max() / np.abs(w1).max())
        info['z_err'] = float(np.abs(z_out - row_orbit(row_out)).max())
    return np.concatenate(amps, axis=1) * 2.0 ** 1.5, info


# ---------------------------------------------------------------------------------------------
# v2.3.0 R6: 偏向磁石の 3 次の生成関数の解析式 (閉形式)。追跡から取り出した係数 (bend_decomposition) との比較用で、K64 の追跡の写像の 2 次のテイラー展開から導いた
# ---------------------------------------------------------------------------------------------
def bend_edge_generators(el):
    """Return closed-form generating functions of the edge kicks of a bend (K0=0, E1, E2, curvature h = ANGLE/L).
    偏向磁石 (K0=0、E1・E2・曲率 h = ANGLE/L) の端の kick の生成関数 (bend_decomposition の f_e1・f_e2 と同じ規約: 端 = A_in ∘ (極薄い偏向磁石 E1)、
    (極薄い偏向磁石 E2) ∘ A_out、座標 (x, px, y, py)) の閉形式。t = tan E:
      入口  f_e1 = (h/2) sec²E1 · px y²  + h t² · x y py − (h t²/2) x² px − (h² t³/2) x y² − (h² t³/6) x³
      出口  f_e2 = −(h/2) sec²E2 · px y² − h t² · x y py + (h t²/2) x² px − (h² t/2) x y² + (h² t³/3) x³
    E=0 では ±(h/2) px y² (硬い端、Maxwell) だけ。K64 の tbendcore の端の式を 2 次まで展開して得た係数を、5 組の (E1, E2, h) の追跡から取り出した値と 1e-6 で確認
    (導出は docs/25 §11.3、試験 test_bend_edge_closed_form)。戻り (f_e1, f_e2) = {(ex, epx, ey, epy): 係数}"""
    L, ang = float(el.L), float(el.ANGLE)
    h = ang / L
    out = []
    for E, sg in ((float(el.E1), +1), (float(el.E2), -1)):
        t = math.tan(E)
        s2 = 1.0 + t * t
        if sg > 0:
            f = {(0, 1, 2, 0): 0.5 * h * s2, (1, 0, 1, 1): h * t * t, (2, 1, 0, 0): -0.5 * h * t * t, (1, 0, 2, 0): -0.5 * h * h * t ** 3, (3, 0, 0, 0): -h * h * t ** 3 / 6.0}
        else:
            f = {(0, 1, 2, 0): -0.5 * h * s2, (1, 0, 1, 1): -h * t * t, (2, 1, 0, 0): 0.5 * h * t * t, (1, 0, 2, 0): -0.5 * h * h * t, (3, 0, 0, 0): h * h * t ** 3 / 3.0}
        out.append(f)
    return out[0], out[1]


def bend_analytic_h3(el, row_in, row_out, monos, ngl=16):
    """Compute the cubic harmonic amplitudes of a bend by analytic integrals, for comparison with the sum from `bend_atoms`.
    偏向磁石 (K0=0、ROTATE なし、閉軌道 0・δ=0) の 3 次の項の調和の振幅 (len(monos),) を、解析式の積分で求める (追跡の写像から取り出す `bend_atoms` の和との比較用)。
    本体 (E=0 の部分偏向磁石 = A_in ∘ S ∘ A_out の連続極限): H_B3 = (h/2) x (px² + py²) + (K1 h/6) x³ を、素子内の w(s) (厳密な線形の転送行列、Gauss-Legendre ngl 点) で積分。
    端: `bend_edge_generators` (入口は row_in、出口は row_out の w で評価)。戻り: 振幅 (2^{3/2} を含まない h の定義そのまま)"""
    from . import resdrive as R
    L, ang = float(el.L), float(el.ANGLE)
    h = ang / L
    k1 = float(el.K1)
    nodes, wts = np.polynomial.legendre.leggauss(ngl)
    nodes, wts = 0.5 * (nodes + 1.0), 0.5 * wts
    f = {(1, 2, 0, 0): 0.5 * h, (1, 0, 0, 2): 0.5 * h, (3, 0, 0, 0): k1 * h / 6.0}
    w0 = w_from_row(row_in)

    def edge_lin(E):          # 端の線形写像 (楔の収束): px += h tanE·x、py −= h tanE·y
        t = h * math.tan(E)
        return np.array([[1, 0, 0, 0], [t, 1, 0, 0], [0, 0, 1, 0], [0, 0, -t, 1]], dtype=float)
    wb = edge_lin(float(el.E1)) @ w0
    res = np.zeros(len(monos), dtype=complex)
    for t, wt in zip(nodes, wts):
        sl = R._slice_for_optics(el, L * float(t))
        M = np.asarray(sl.transfer_matrix(0.0), dtype=float)[:4, :4]
        res += L * wt * harmonic_amplitudes(f, M @ wb, monos)[:, 0]
    wend = np.asarray(R._slice_for_optics(el, L).transfer_matrix(0.0), dtype=float)[:4, :4] @ wb
    fe1, fe2 = bend_edge_generators(el)
    res += harmonic_amplitudes(fe1, w0, monos)[:, 0] + harmonic_amplitudes(fe2, wend, monos)[:, 0]
    return res
