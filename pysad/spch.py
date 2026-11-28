# -*- coding: utf-8 -*-
"""Space-charge tracking newly written for v2: SPAC (3D), PSPAC (2D Poisson mesh) and WSPAC (Gaussian beam-matrix) kicks.
spch.py — 空間電荷の追跡 (v2 で新しく実装。K64 の spkick.f / tspch_.c + spch.c の物理を、正しい係数で書き直したもの)。

K64 には空間電荷のモデルが 3 つある (詳細と K64 の誤りは docs/history/v2_changes_C_spch.md)。

- ``SPAC`` フラグ + ``DRIFT`` の ``RADIUS>0`` (tturn.f:288 → spkick.f spdrift_free、SOL 領域は tsol.f:639 → spdrift_solenoid):
  ドリフトを 5 cm ごとに割り、その間に「バンチの静止系での 3 次元のポアソン方程式 (半径 RADIUS の円筒の導体壁) の場」のキックを入れる。
  RADIUS<0 は半径 |RADIUS| の絞りだけ、RADIUS=0 は何もしない。
- ``PSPAC`` フラグ + ``SPCH`` 要素 (tturn.f:521 → tspch_.c tpspac → spch.c): 横 2 次元のポアソン方程式 (PSPACNX×PSPACNY の矩形メッシュ、
  境界で電位 0) の電位 Φ(x,y) と縦の線密度 λ(z) の積で、前の SPCH からの距離の分のキックを SPCH の位置で入れる。
- ``WSPAC`` (tintrb.f twspac): ビーム行列 (Emittance の結果) のガウス分布の力を、すべての要素の入口で入れる (強-弱のモデル)。
- ``SPAC`` + ``MULT`` の ``RADIUS`` (tmulti.f:222): MULT の分割の間に SPAC のキックを挟む (K64 は spkick の引数が合っていない)。

このファイルの実装 (v2 モードで動く。k64 モードでは今までどおり空間電荷を入れない):

**SPAC (3 次元)** `Solver3D`。バンチの静止系 (ζ = γ0 β0 c·(到着時刻の差)) で、マクロ粒子を 2 次の B スプライン (TSC) でメッシュに載せ、
ζ 方向は周期 FFT、横方向はモードごとの自由空間のグリーン関数 2 K0(k d) (k=0 は −2 ln d。セルの平均) の畳み込み (FFT) で解く。
円筒の壁 (半径 R、電位 0) は、方位角の調和関数の級数 −2 Σ_m ε_m K_m(kR)/I_m(kR)·I_m(kr) I_m(kr') cos m(θ−θ') で厳密に入れる。
電位 φ' (単位: 個数/長さ。φ' = Σ 1/r') から、長さ L のキックは

    Δ(px/p0) = −C ∂φ'/∂x,  Δ(py/p0) = −C ∂φ'/∂y,  C = r_c w L / (β0² γ0²)
    Δγ       = −r_c w L ∂φ'/∂ζ

(r_c = q² r_e m_e/m: 粒子の古典半径、w = PBUNCH/NP0: マクロ粒子 1 個の実粒子数、β0 γ0: バンチの平均)。一様な長い円柱 (半径 a、線密度 λ)
では x'' = K x / a²、K = 2 r_c λ /(β² γ³) (ペルビアンス) になる。電位の補間と電荷の配分に同じ形状関数を使い、キックを電位の勾配に
するので、凍結した電位のキックは (x, px, y, py, z, δ) で厳密にシンプレクティック。

**PSPAC (横 2 次元 × 線密度)** `Solver2D`。K64 と同じモデル (矩形の境界、5 点差分のポアソン方程式) を DST で解く。
Δ(px/p0) = −(r_c w L/(β²γ³)) λ(z) ∂Φ/∂x / N、Δγ = −(r_c w L/γ²) λ'(z) Φ / N (Φ = N 個の −2 ln d の和、λ: 個数/m)。

**MULT の RADIUS (SPAC)** `SpaceCharge._mult`。MULT (ANGLE=0、RF なし、SOL 領域の外) を ndiv = max(MULT の分割, |L|/SPACSTEP/EPS) に割り、
本体 (四極 + 双極の厳密解) の半分 → [多極のキック、SPAC のキック (L/ndiv)] → 本体 → … → 本体の半分。場は DRIFT と同じ `Solver3D`。

**WSPAC (ガウス分布の強-弱のモデル)** `WSpace`。要素 l の入口のビーム行列 Σ から、z で条件を付けた横の分散
(Σxx − Σxz²/Σzz など) の主軸と σ1、σ2 を取り、Bassetti-Erskine の場 F (単位の線密度で 2/r) と電位 U (中心で 0) で

    Δ(px/p0) = +A(z) F_x,  Δδ = (z − zc)/σz² · A(z) · (U_wall − U),  A(z) = N r_c L exp(−(z−zc)²/2σz²) / (β²γ³ √(2π) σz)

(L は「前の要素の中点からこの要素の中点まで」。区間の終わりに残りを入れる)。ビーム行列は Emittance の平衡値 (K64 と同じ) か、
`WSPACBEAM = {εx, εy, σz, σδ}` (v2 で追加) と Twiss から作る。

使える範囲・使えない範囲は docs/history/v2_changes_C_spch.md と v2_changes_C2_spch.md。
"""
from __future__ import annotations

import math

import numpy as np
from scipy import fft as _fft
from scipy import special as _sp

from . import quadbody as _qb
from .sadscript import SADError

ELRADI = 2.8179403205e-15        # 電子の古典半径 [m]
ELMASS = 0.51099895e6            # 電子の質量 [eV]
ALSTEP = 0.05                    # spkick.f spdrift_free の alstep [m]
MMAX = 40                        # 壁の像の方位角の調和の上限
NZMAX = 4096                     # ζ 方向の周期メッシュの上限
_QSTEP = math.log(1.04)          # メッシュ幅の量子化 (グリーン関数の表を使い回すため)

_GAUSS6 = np.polynomial.legendre.leggauss(6)
_GAUSS16 = np.polynomial.legendre.leggauss(16)
_EULER = 0.5772156649015329
_CSMOOTH_T = 7.0 / 12.0          # 横: 形状関数のぼかしの分散 (メッシュ幅² の単位。TSC 1/4 ×2 + セルの平均 1/12)
_CSMOOTH_Z = 0.5                 # 縦: TSC 1/4 ×2


class SpaceChargeError(SADError):
    """Exception raised for errors in the space-charge (SPAC/PSPAC) calculation.
    空間電荷 (SPAC/PSPAC) の計算のエラー。"""
    pass


# ============================================================================ 形状関数 (TSC: 2 次の B スプライン)
def _tsc(u):
    """u: ノード番号の単位の座標。戻り: (中央のノード番号 j, 重み (3,n), 重みの u 微分 (3,n))。ノード j−1, j, j+1 に載る"""
    j = np.floor(u + 0.5).astype(np.int64)
    d = u - j
    w = np.empty((3, u.size))
    dw = np.empty((3, u.size))
    w[0] = 0.5 * (0.5 - d) ** 2
    w[1] = 0.75 - d * d
    w[2] = 0.5 * (0.5 + d) ** 2
    dw[0] = d - 0.5
    dw[1] = -2.0 * d
    dw[2] = d + 0.5
    return j, w, dw


def _quant(h):
    """メッシュ幅を 4% 刻みの等比列に切り上げる"""
    return math.exp(math.ceil(math.log(h) / _QSTEP - 1e-9) * _QSTEP)


def _axis(v, nn, hw_floor=0.0):
    """1 方向のメッシュ: 中心 c、幅 h、粒子の座標 (ノード番号の単位、[0.5, nn−1.5] に切る)。
    半幅は min(最大の |v−c|, 4σ)。その外の粒子 (ガウス分布で 6e-5) は端に寄せて扱う"""
    c = float(v.mean())
    dv = v - c
    hw = min(float(np.abs(dv).max()), 4.0 * float(dv.std()))
    hw = max(hw, hw_floor, 1e-300)
    h = _quant(hw * 1.02 / (nn / 2 - 1))
    u = dv / h + 0.5 * (nn - 1)
    np.clip(u, 0.5 + 1e-9, nn - 1.5 - 1e-9, out=u)
    return c, h, u, hw


# ============================================================================ 自由空間のグリーン関数 (セルの平均)
def _avg_mln2(X, Y, hx, hy):
    """−ln(x²+y²) の、(X,Y) を中心とする hx×hy の矩形の上の平均 (解析式)"""
    def F(x, y):
        return x * y * np.log(x * x + y * y) - 3.0 * x * y + x * x * np.arctan(y / x) + y * y * np.arctan(x / y)
    x1, x2 = X - 0.5 * hx, X + 0.5 * hx
    y1, y2 = Y - 0.5 * hy, Y + 0.5 * hy
    return -(F(x2, y2) - F(x1, y2) - F(x2, y1) + F(x1, y1)) / (hx * hy)


def _self_k0(k, a, b):
    """2 K0(k d) の、原点を中心とする 2a×2b の矩形の上の平均。原点の対数特異点は ∫_0^1 K0(κu) u du = (1−κ K1(κ))/κ² で厳密に積分する"""
    gv, gw = _GAUSS16

    def inner(kap):
        out = np.empty_like(kap)
        sm = kap < 1e-3
        ks = kap[sm]
        out[sm] = -0.5 * (np.log(0.5 * ks) + _EULER - 0.5)
        kl = kap[~sm]
        out[~sm] = (1.0 - kl * _sp.k1(kl)) / (kl * kl)
        return out
    c1 = np.sqrt(a * a + b * b * gv * gv)
    c2 = np.sqrt(a * a * gv * gv + b * b)
    i1 = (2.0 * inner(k[:, None] * c1[None, :]) * gw[None, :]).sum(axis=1)
    i2 = (2.0 * inner(k[:, None] * c2[None, :]) * gw[None, :]).sum(axis=1)
    return 0.5 * (i1 + i2)


_GCACHE = {}


def _green_hat(nx, ny, hx, hy, k):
    """自由空間のグリーン関数の表の 2 次元 FFT (モード, 2nx, 2ny)。k[0]=0 (−2 ln d)、k>0 は 2 K0(k d)"""
    key = (nx, ny, hx, hy, k.size, float(k[1]) if k.size > 1 else 0.0)
    g = _GCACHE.get(key)
    if g is not None:
        return g
    ii = np.arange(-(nx - 1), nx)
    jj = np.arange(-(ny - 1), ny)
    X = (ii * hx)[:, None] + 0.0 * jj[None, :]
    Y = 0.0 * ii[:, None] + (jj * hy)[None, :]
    nm = k.size
    G = np.empty((nm, 2 * nx - 1, 2 * ny - 1))
    G[0] = _avg_mln2(X, Y, hx, hy)
    if nm > 1:
        kk = k[1:]
        d = np.hypot(X, Y)
        d[nx - 1, ny - 1] = 1.0
        G[1:] = 2.0 * _sp.k0(kk[:, None, None] * d[None])
        gx, gw = _GAUSS6
        for di in range(-2, 3):
            for dj in range(-2, 3):
                if (di == 0 and dj == 0) or abs(di) > nx - 1 or abs(dj) > ny - 1:
                    continue
                px = di * hx + 0.5 * hx * gx
                py = dj * hy + 0.5 * hy * gx
                r = np.hypot(px[:, None], py[None, :])
                val = 2.0 * _sp.k0(kk[:, None, None] * r[None])
                G[1:, nx - 1 + di, ny - 1 + dj] = 0.25 * np.einsum('nij,i,j->n', val, gw, gw)
        G[1:, nx - 1, ny - 1] = _self_k0(kk, 0.5 * hx, 0.5 * hy)
    # 形状関数 (配分 TSC・補間 TSC・セルの平均) のぼかしを 1 次まで戻す: G ← G − (c/2) δ²G (δ²: 2 階差分、c = 1/4+1/4+1/12)
    D = np.zeros_like(G)
    D[:, 1:-1, :] += G[:, 2:, :] - 2.0 * G[:, 1:-1, :] + G[:, :-2, :]
    D[:, :, 1:-1] += G[:, :, 2:] - 2.0 * G[:, :, 1:-1] + G[:, :, :-2]
    G -= 0.5 * _CSMOOTH_T * D
    W = np.zeros((nm, 2 * nx, 2 * ny))
    W[:, (ii % (2 * nx))[:, None], (jj % (2 * ny))[None, :]] = G
    g = (_fft.fft2(W, axes=(1, 2)).real,        # 実の偶関数
         G[:, nx - 1:nx + 2, ny - 1:ny + 2].copy())       # ずれ 0..2 の値 (自己力を引くのに使う)
    if len(_GCACHE) >= 12:
        _GCACHE.pop(next(iter(_GCACHE)))
    _GCACHE[key] = g
    return g


def _image(Qh, k, X, Y, R):
    """円筒の壁 (半径 R、軸は x=y=0、電位 0) の像の電位 (モード, nx, ny)。Qh: 電荷のモード (モード, nx, ny)"""
    r = np.hypot(X, Y)
    th = np.arctan2(Y, X)
    rc = np.minimum(r, R)
    mask = np.abs(Qh).sum(axis=0) > 0.0
    out = np.zeros_like(Qh)
    if not mask.any():
        return out
    rmax = float(rc[mask].max())
    t = (rmax / R) ** 2
    if t <= 0.0:
        nmh = 1
    elif t >= 1.0:
        nmh = MMAX
    else:
        nmh = int(min(MMAX, max(2, math.ceil(math.log(1e-10) / math.log(t)))))
    rho = rc / R
    out[0] += 2.0 * math.log(R) * Qh[0].sum()
    sel = np.nonzero((k > 0.0) & (k * (R - rmax) < 18.0))[0]
    ks = k[sel]
    for m in range(0, nmh + 1):
        cm = np.cos(m * th)
        sm = np.sin(m * th)
        if m >= 1:
            pw = rho ** m
            C = (Qh[0] * pw * cm).sum()
            S = (Qh[0] * pw * sm).sum()
            out[0] -= (2.0 / m) * pw * (cm * C + sm * S)
        if sel.size == 0:
            continue
        xR = ks * R
        iR = _sp.ive(m, xR)
        kR = _sp.kve(m, xR)
        with np.errstate(all='ignore'):
            coef = kR * iR                              # K_m(kR) I_m(kR)
            B = _sp.ive(m, ks[:, None, None] * rc[None]) / iR[:, None, None] * np.exp(ks[:, None, None] * (rc[None] - R))
        bad = ~(np.isfinite(coef) & (iR > 0.0))
        if m > 0 and bad.any():                                   # kR が小さく I_m が下にあふれる: 小さい引数の極限
            coef[bad] = 0.5 / m
            B[bad] = (rho ** m)[None]
        eps = 1.0 if m == 0 else 2.0
        C = (Qh[sel] * (B * cm[None])).sum(axis=(1, 2))
        S = (Qh[sel] * (B * sm[None])).sum(axis=(1, 2))
        out[sel] -= 2.0 * eps * coef[:, None, None] * B * (cm[None] * C[:, None, None] + sm[None] * S[:, None, None])
    return out


# ============================================================================ 3 次元 (SPAC)
class Field3D:
    """Holds the rest-frame potential phi' on the mesh and its interpolation.
    メッシュの上の静止系の電位 φ' (単位: マクロ粒子の個数/長さ) と、その補間"""

    def __init__(self, phi, xc, yc, hx, hy, hz, nzb):
        self.phi = phi                       # (nx, ny, Nz)
        self.xc, self.yc, self.hx, self.hy, self.hz, self.nzb = xc, yc, hx, hy, hz, nzb
        self.zc = 0.0

    def coords(self, x, y, zeta):
        """Convert the physical coordinates (x, y, zeta) to fractional grid indices (ux, uy, uz), clipped to the interior of the mesh.
        物理座標 (x, y, zeta) を格子の小数の番号 (ux, uy, uz) にする (格子の内側に収まるよう clip する)。"""
        nx, ny, _ = self.phi.shape
        ux = np.clip((x - self.xc) / self.hx + 0.5 * (nx - 1), 0.5 + 1e-9, nx - 1.5 - 1e-9)
        uy = np.clip((y - self.yc) / self.hy + 0.5 * (ny - 1), 0.5 + 1e-9, ny - 1.5 - 1e-9)
        uz = np.clip((zeta - self.zc) / self.hz + 0.5 * (self.nzb - 1), 0.5 + 1e-9, self.nzb - 1.5 - 1e-9)
        return ux, uy, uz

    def grad(self, x, y, zeta, selfsub=False):
        """Return (phi', dphi'/dx, dphi'/dy, dphi'/dzeta) at particle positions; selfsub=True removes the self-force of the macro-particles.
        (φ', ∂φ'/∂x, ∂φ'/∂y, ∂φ'/∂ζ) を粒子の位置で (メッシュの外の粒子は端の値)。
        selfsub=True: この電位を作ったマクロ粒子自身に使うとき。粒子が自分で作った電位の分 (メッシュに対する位置で決まる
        見かけの自己力) を引く。残りは「他の粒子との対のエネルギーの勾配」で、ハミルトン系のまま"""
        ux, uy, uz = self.coords(x, y, zeta)
        u, gx, gy, gz = _interp3(self.phi, ux, uy, uz, self.hx, self.hy, self.hz)
        if selfsub:
            su, sx, sy, sz = _self3(self.kself, ux, uy, uz)
            u -= su
            gx -= sx / self.hx
            gy -= sy / self.hy
            gz -= sz / self.hz
        return u, gx, gy, gz


def _interp3(phi, ux, uy, uz, hx, hy, hz):
    jx, wx, dwx = _tsc(ux)
    jy, wy, dwy = _tsc(uy)
    jz, wz, dwz = _tsc(uz)
    n = ux.size
    u = np.zeros(n)
    gx = np.zeros(n)
    gy = np.zeros(n)
    gz = np.zeros(n)
    for a in range(3):
        ia = jx + (a - 1)
        for b in range(3):
            ib = jy + (b - 1)
            for c in range(3):
                p = phi[ia, ib, jz + (c - 1)]
                u += wx[a] * wy[b] * wz[c] * p
                gx += dwx[a] * wy[b] * wz[c] * p
                gy += wx[a] * dwy[b] * wz[c] * p
                gz += wx[a] * wy[b] * dwz[c] * p
    return u, gx / hx, gy / hy, gz / hz


_OFF3 = np.array([(a, b, c) for a in range(3) for b in range(3) for c in range(3)])
_DIF3 = np.abs(_OFF3[:, None, :] - _OFF3[None, :, :])          # (27, 27, 3)


def _self3(kself, ux, uy, uz):
    """粒子が自分で作った電位の、自分の位置での値と勾配 (ノードの単位): Σ_ab S_a K_ab S_b、Σ_ab ∇S_a K_ab S_b"""
    _, wx, dwx = _tsc(ux)
    _, wy, dwy = _tsc(uy)
    _, wz, dwz = _tsc(uz)
    K = kself[_DIF3[:, :, 0], _DIF3[:, :, 1], _DIF3[:, :, 2]]   # (27, 27)
    a, b, c = _OFF3[:, 0], _OFF3[:, 1], _OFF3[:, 2]
    W = wx[a] * wy[b] * wz[c]                                   # (27, n)
    V = K @ W
    return ((W * V).sum(axis=0), (dwx[a] * wy[b] * wz[c] * V).sum(axis=0),
            (wx[a] * dwy[b] * wz[c] * V).sum(axis=0), (wx[a] * wy[b] * dwz[c] * V).sum(axis=0))


def _deposit3(ux, uy, uz, nx, ny, nz):
    jx, wx, _ = _tsc(ux)
    jy, wy, _ = _tsc(uy)
    jz, wz, _ = _tsc(uz)
    q = np.zeros(nx * ny * nz)
    for a in range(3):
        ia = (jx + (a - 1)) * ny
        for b in range(3):
            iab = (ia + jy + (b - 1)) * nz
            wab = wx[a] * wy[b]
            for c in range(3):
                q += np.bincount(iab + jz + (c - 1), weights=wab * wz[c], minlength=q.size)
    return q.reshape(nx, ny, nz)


class Solver3D:
    """Solver for SPAC; nx, ny are even transverse node counts and nz is the even number of zeta nodes (None: from particle count).
    SPAC のソルバー。nx, ny: 横のノード数 (偶数)、nz: バンチの中の ζ のノード数 (偶数)。None なら粒子数から決める"""

    def __init__(self, nx=None, ny=None, nz=None, warn=None):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.warn = warn
        self._warned = set()

    def _w(self, key, msg):
        if key not in self._warned:
            self._warned.add(key)
            if self.warn is not None:
                self.warn(msg)

    @staticmethod
    def auto_mesh(n):
        """Choose the mesh (nx, nx, nz) automatically for n particles (each between 24 and 64, even).
        粒子数 n から格子 (nx, nx, nz) を自動で決める (各 24〜64、偶数)。"""
        nz = int(min(64, max(24, 2 * round(n ** (1.0 / 3.0)))))
        nx = int(min(64, max(24, 2 * round(math.sqrt(n / nz)))))
        return nx, nx, nz

    def solve(self, x, y, zeta, R):
        """Solve the potential of macro-particles (x, y, zeta) in a cylinder of wall radius R; return Field3D.
        マクロ粒子 (x, y, ζ: 静止系の縦の位置) の電位。R: 円筒の壁の半径 (>0)。戻り: Field3D"""
        n = x.size
        ax, ay, az = self.auto_mesh(n)
        nx = int(self.nx or ax)
        ny = int(self.ny or ay)
        nzb = int(self.nz or az)
        zc = float(zeta.mean())
        if float(np.abs(zeta - zc).max()) <= 0.0:
            raise SpaceChargeError('SPAC: バンチ長が 0 (すべての粒子の z が同じ)。空間電荷の線密度が決まらない。'
                                   '粒子に z の広がりを与え、RFSW を ON にする (NORFSW は区間ごとに z を 0 にする)')
        xc, hx, ux, hwx = _axis(x, nx)
        yc, hy, uy, hwy = _axis(y, ny)
        if hwx < 1e-6 * hwy or hwy < 1e-6 * hwx:            # 片方の幅が 0 (線状のビーム): もう一方の 1e-6 倍を下限にする
            fl = 1e-6 * max(hwx, hwy)
            if fl <= 0.0:
                fl = 1e-9 * R
            xc, hx, ux, hwx = _axis(x, nx, fl)
            yc, hy, uy, hwy = _axis(y, ny, fl)
        _, hz, uz, hwz = _axis(zeta, nzb)
        pad = min(3.0 * R, max(10.0 * max(hwx, hwy), 3.0 * hwz))
        nzt = nzb + int(math.ceil(pad / hz))
        if nzt > NZMAX:
            self._w('pad', f'SPAC: バンチが真空ダクト・ビームの幅に比べて短く、ζ の周期メッシュが上限 {NZMAX} を超える。'
                           f'周期の像までの距離を {pad:.3g} m から {(NZMAX - nzb) * hz:.3g} m (静止系) に縮めた (像の誤差が増える)')
            nzt = NZMAX
        nzt = int(_fft.next_fast_len(nzt, real=True))
        Q = _deposit3(ux, uy, uz, nx, ny, nzt)
        Qh = np.moveaxis(_fft.rfft(Q, axis=2), 2, 0)        # (モード, nx, ny)
        nm = Qh.shape[0]
        k = 2.0 * math.pi * np.arange(nm) / (nzt * hz)
        Gh, Gs = _green_hat(nx, ny, hx, hy, k)
        P = np.zeros((nm, 2 * nx, 2 * ny), dtype=complex)
        P[:, :nx, :ny] = Qh
        ph = _fft.ifft2(_fft.fft2(P, axes=(1, 2)) * Gh, axes=(1, 2))[:, :nx, :ny]
        X = (xc + (np.arange(nx) - 0.5 * (nx - 1)) * hx)[:, None] + np.zeros(ny)[None, :]
        Y = np.zeros(nx)[:, None] + (yc + (np.arange(ny) - 0.5 * (ny - 1)) * hy)[None, :]
        if float(np.hypot(X, Y)[np.abs(Qh).sum(axis=0) > 0.0].max()) > 0.9 * R:
            self._w('wall', 'SPAC: ビームが真空ダクトの半径 RADIUS の 0.9 倍より外まで広がっている。壁の像の級数 '
                            f'(方位角の調和 {MMAX} 次まで) の精度が落ちる')
        ph = ph + _image(Qh, k, X, Y, R)
        sharp = 1.0 + _CSMOOTH_Z * (1.0 - np.cos(k * hz))                       # ζ 方向のぼかしを 1 次まで戻す
        ph *= sharp[:, None, None]
        phi = _fft.irfft(np.moveaxis(ph, 0, 2), n=nzt, axis=2) / hz
        f = Field3D(phi, xc, yc, hx, hy, hz, nzb)
        f.zc = zc
        # ノードの単位電荷が (di, dj, dl) だけ離れたノードに作る電位 (自由空間の部分。0 ≤ di, dj, dl ≤ 2)
        f.kself = _fft.irfft(Gs * sharp[:, None, None], n=nzt, axis=0)[:3].transpose(1, 2, 0) / hz
        return f


# ============================================================================ 横 2 次元 × 線密度 (PSPAC)
_OFF2 = np.array([(a, b) for a in range(3) for b in range(3)])
_DIF2 = np.abs(_OFF2[:, None, :] - _OFF2[None, :, :])


class Field2D:
    """Holds the 2D potential Phi(x, y) (zero at the rectangular boundary) and the line density lambda(z) in the lab frame.
    Φ(x,y) (全マクロ粒子の −2 ln d の和に当たる 2 次元の電位、矩形の境界で 0) と線密度 λ(z) (個数/m、実験室系)"""

    def __init__(self, phi, dx, dy, lam, z0, dz, nxy):
        self.phi, self.dx, self.dy = phi, dx, dy        # phi: (nx+1, ny+1)。ノード i は x = (i − nx/2) dx
        self.nxy = nxy                                  # Φ に載せたマクロ粒子の数
        self.lam, self.z0, self.dz = lam, z0, dz        # lam: (nzn,)。ノード k は z = z0 + k dz

    def eval(self, x, y, z, selfsub=False):
        """Evaluate (Phi, dPhi/dx, dPhi/dy, lambda, dlambda/dz, inside flag); selfsub=True removes the self contribution.
        (Φ, ∂Φ/∂x, ∂Φ/∂y, λ, dλ/dz, 中にあるか)。selfsub=True: Φ を作った粒子自身に使うとき、自分の作った Φ の分を引く"""
        nx = self.phi.shape[0] - 1
        ny = self.phi.shape[1] - 1
        ux = x / self.dx + nx // 2
        uy = y / self.dy + ny // 2
        inside = (ux >= 0.5) & (ux < nx - 0.5) & (uy >= 0.5) & (uy < ny - 0.5)
        ux = np.clip(ux, 0.5 + 1e-9, nx - 0.5 - 1e-9)
        uy = np.clip(uy, 0.5 + 1e-9, ny - 0.5 - 1e-9)
        jx, wx, dwx = _tsc(ux)
        jy, wy, dwy = _tsc(uy)
        n = x.size
        u = np.zeros(n)
        gx = np.zeros(n)
        gy = np.zeros(n)
        for a in range(3):
            for b in range(3):
                p = self.phi[jx + (a - 1), jy + (b - 1)]
                u += wx[a] * wy[b] * p
                gx += dwx[a] * wy[b] * p
                gy += wx[a] * dwy[b] * p
        if selfsub:
            K = self.kself[_DIF2[:, :, 0], _DIF2[:, :, 1]]
            ia, ib = _OFF2[:, 0], _OFF2[:, 1]
            W = wx[ia] * wy[ib]
            V = K @ W
            u -= (W * V).sum(axis=0)
            gx -= (dwx[ia] * wy[ib] * V).sum(axis=0)
            gy -= (wx[ia] * dwy[ib] * V).sum(axis=0)
        nzn = self.lam.size
        uz = (z - self.z0) / self.dz
        inz = (uz >= 0.5) & (uz < nzn - 1.5)
        uzc = np.clip(uz, 0.5 + 1e-9, nzn - 1.5 - 1e-9)
        jz, wz, dwz = _tsc(uzc)
        lam = np.zeros(n)
        dlam = np.zeros(n)
        for c in range(3):
            p = self.lam[jz + (c - 1)]
            lam += wz[c] * p
            dlam += dwz[c] * p
        lam[~inz] = 0.0
        dlam[~inz] = 0.0
        return u, gx / self.dx, gy / self.dy, lam, dlam / self.dz, inside


class Solver2D:
    """Solver for PSPAC with the same model as spch.c in K64; nx, ny are cell counts, dx, dy widths, and Phi = 0 on the boundary.
    PSPAC のソルバー (K64 の spch.c と同じモデル)。nx, ny: セル数 (ノードは nx+1)、dx, dy: 幅。境界 (ノード 0 と nx) で Φ=0。
    nz, dz: 縦のメッシュ (K64: 中心 0、ノード k は z = (k − nz/2) dz)。nz<4 なら、バンチに合わせて自動で 64 ノードを張る"""

    def __init__(self, nx=128, ny=128, nz=1, dx=1e-3, dy=1e-3, dz=1e-3, warn=None):
        self.warn = warn
        self._warned = set()
        if nx < 4 or ny < 4:
            raise SpaceChargeError('PSPAC: PSPACNX, PSPACNY は 4 以上')
        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
        self.dx, self.dy, self.dz = float(dx), float(dy), float(dz)
        px = np.arange(1, self.nx)
        py = np.arange(1, self.ny)
        ex = 2.0 - 2.0 * np.cos(math.pi * px / self.nx)
        ey = 2.0 - 2.0 * np.cos(math.pi * py / self.ny)
        self._eig = (ex / self.dx ** 2)[:, None] + (ey / self.dy ** 2)[None, :]
        # 形状関数 (TSC ×2: 分散 1/2) のぼかしと 5 点差分の誤差 (−1/6) を 1 次まで戻す係数
        self._sharp = 1.0 + (1.0 / 6.0) * (ex[:, None] + ey[None, :])
        # 中央のノードの単位電荷が (di, dj) だけ離れたノードに作る Φ (0 ≤ di, dj ≤ 2。自己力を引くのに使う)
        q = np.zeros((self.nx + 1, self.ny + 1))
        ic, jc = self.nx // 2, self.ny // 2
        q[ic, jc] = 1.0
        self.kself = self._poisson(q)[ic:ic + 3, jc:jc + 3].copy()

    def _poisson(self, q):
        nx, ny = self.nx, self.ny
        rhs = 4.0 * math.pi * q[1:nx, 1:ny] / (self.dx * self.dy)       # ∇²Φ = −4π σ (Φ = Σ −2 ln d)
        phi = np.zeros((nx + 1, ny + 1))
        phi[1:nx, 1:ny] = _fft.idstn(_fft.dstn(rhs, type=1) / self._eig * self._sharp, type=1)
        return phi

    def _w(self, key, msg):
        if key not in self._warned:
            self._warned.add(key)
            if self.warn is not None:
                self.warn(msg)

    def check_resolution(self, x, y, z):
        """Check the mesh resolution against the bunch size, since a coarse mesh underestimates the kick (a source of the gs difference).
        メッシュに対してバンチが粗いと、キックが小さく出る (K64 の gs との差もここから: docs/history/v2next_A1.md 項目 8)。
        一様な円のビームで、半径 a が 1 セル (DX) で −13%、2 セルで −3%、4 セルで −1%、10 セルで −0.6%。縦のメッシュを与えたとき
        (PSPACNZ≥4) は、一様なバンチの長さ lb が 4 セル (DZ) で −12%、10 セルで −4%、20 セルで −2%"""
        sx, sy = float(np.std(x)), float(np.std(y))
        if min(sx / self.dx, sy / self.dy) < 1.5:
            self._w('xy', ' ???-PSPAC-Warning-the transverse mesh is too coarse: rms size x, y = %.3g, %.3g m is only %.1f, %.1f cells of '
                    'PSPACDX, PSPACDY = %.3g, %.3g m (the kick is under-estimated by more than 10%%). Use PSPACDX, PSPACDY <= rms/3.'
                    % (sx, sy, sx / self.dx, sy / self.dy, self.dx, self.dy) + chr(10))
        ins = (np.abs(x) < 0.5 * self.nx * self.dx) & (np.abs(y) < 0.5 * self.ny * self.dy)
        if ins.size and (1.0 - ins.mean()) > 0.01:
            self._w('out', ' ???-PSPAC-Warning-%.1f%% of the particles are outside the mesh (PSPACNX x PSPACDX, PSPACNY x PSPACDY); they get no kick.'
                    % (100.0 * (1.0 - ins.mean())) + chr(10))
        if self.nz >= 4:
            sz = float(np.std(z))
            if sz / self.dz < 1.5:
                self._w('z', ' ???-PSPAC-Warning-the longitudinal mesh is too coarse: rms length %.3g m is only %.1f cells of PSPACDZ = %.3g m '
                        '(the line density is under-estimated). Use PSPACDZ <= rms/3 or PSPACNZ < 4 (automatic mesh).'
                        % (sz, sz / self.dz, self.dz) + chr(10))

    def solve(self, x, y, z):
        """Deposit the particles (x, y, z) on the 2D mesh with the TSC scheme, solve the Poisson equation and return the Field2D with the potential and the line density.
        粒子 (x, y, z) を TSC 法で 2 次元格子に割り振り、ポアソン方程式を解いて、電位と線密度を持つ Field2D を返す。"""
        self.check_resolution(x, y, z)
        nx, ny = self.nx, self.ny
        ux = x / self.dx + nx // 2
        uy = y / self.dy + ny // 2
        ins = (ux >= 0.5) & (ux < nx - 0.5) & (uy >= 0.5) & (uy < ny - 0.5)
        jx, wx, _ = _tsc(ux[ins])
        jy, wy, _ = _tsc(uy[ins])
        q = np.zeros((nx + 1) * (ny + 1))
        for a in range(3):
            ia = (jx + (a - 1)) * (ny + 1)
            for b in range(3):
                q += np.bincount(ia + jy + (b - 1), weights=wx[a] * wy[b], minlength=q.size)
        q = q.reshape(nx + 1, ny + 1)
        phi = self._poisson(q)
        # 線密度
        if self.nz >= 4:
            nzn = self.nz + 1
            dz = self.dz
            z0 = -(self.nz // 2) * dz
        else:
            nzn = 64
            zc = float(z.mean())
            hw = float(np.abs(z - zc).max())
            if hw <= 0.0:
                raise SpaceChargeError('PSPAC: バンチ長が 0 (すべての粒子の z が同じ)。粒子に z の広がりを与え、RFSW を ON にする')
            dz = hw * 1.02 / (nzn / 2 - 1.5)
            z0 = zc - 0.5 * (nzn - 1) * dz
        uz = (z - z0) / dz
        inz = (uz >= 0.5) & (uz < nzn - 1.5)
        jz, wz, _ = _tsc(uz[inz])
        lam = np.zeros(nzn)
        for c in range(3):
            lam += np.bincount(jz + (c - 1), weights=wz[c], minlength=nzn)
        f = Field2D(phi, self.dx, self.dy, lam / dz, z0, dz, float(q.sum()))
        f.kself = self.kself
        return f


# ============================================================================ キック (状態 S: x, px/p, y, py/p, z, δ, dv, …, dvfs)
def _p2h(p):
    return np.sqrt(1.0 + p * p)


def _apply(S, n, dpx, dpy, dgam, p0, h0):
    """正準運動量のキック (dpx, dpy: p0 単位) とエネルギーのキック dgam (γ の変化) を S に入れる。
    到着時刻は変えない (SAD の z = −v t なので z に v'/v を掛ける)"""
    dl = S[:n, 5]
    pr = 1.0 + dl
    p = p0 * pr
    h = _p2h(p)
    h1 = np.maximum(h + dgam, 1.0 + 1e-12)
    p1 = np.sqrt((h1 - 1.0) * (h1 + 1.0))
    dp = (h1 - h) * (h1 + h) / (p1 + p)              # p1 − p
    dl1 = dl + dp / p0
    pr1 = 1.0 + dl1
    S[:n, 1] = (S[:n, 1] * pr + dpx) / pr1
    S[:n, 3] = (S[:n, 3] * pr + dpy) / pr1
    S[:n, 4] = S[:n, 4] * (p1 / h1) / (p / h)
    S[:n, 5] = dl1
    S[:n, 6] = -(p1 + p0) / h1 / (p0 * h1 + p1 * h0) * dl1 + S[:n, 10]


def bunch_frame(S, n, p0):
    """Return the bunch mean beta0 and gamma0 and each particle's rest-frame longitudinal position zeta = gamma0 beta0 z/beta_i.
    バンチの平均の β0、γ0 と、粒子の静止系の縦の位置 ζ = γ0 β0 z/β_i (z/β_i = −c t)"""
    dl = S[:n, 5]
    pm = p0 * (1.0 + float(dl.mean()))
    g0 = math.sqrt(1.0 + pm * pm)
    b0 = pm / g0
    p = p0 * (1.0 + dl)
    bi = p / _p2h(p)
    return b0, g0, g0 * b0 * S[:n, 4] / bi


def kick3d(S, n, field, al, rcw, p0, h0, b0, g0, zeta, selfsub=True):
    """Apply the SPAC kick with rcw = r_c (PBUNCH/NP0) and length al [m]; selfsub is True when S is the same set that made the field.
    SPAC のキック。rcw = r_c·(PBUNCH/NP0)、al: 長さ [m]。selfsub: S の粒子が field を作った粒子そのもののとき True"""
    _, gx, gy, gz = field.grad(S[:n, 0], S[:n, 2], zeta, selfsub=selfsub)
    c = rcw * al / (b0 * g0) ** 2
    _apply(S, n, -c * gx, -c * gy, -rcw * al * gz, p0, h0)


def kick2d(S, n, field, al, rcw, p0, h0, b0, g0, selfsub=True):
    """Apply the PSPAC kick corresponding to the rest-frame potential phi' = (lambda(z)/gamma0) Phi(x, y)/N.
    PSPAC のキック。静止系の電位 φ' = (λ(z)/γ0)·Φ(x,y)/N (N = Φ に載せたマクロ粒子の数) に当たる:
    Δ(px/p0) = −r_c w L/(β²γ³)·λ ∂Φ/∂x / N、Δγ = −r_c w L/γ²·(dλ/dz) Φ / N"""
    dl = S[:n, 5]
    p = p0 * (1.0 + dl)
    bi = p / _p2h(p)
    zb = b0 * S[:n, 4] / bi                           # バンチの速度での位置 (= −β0 c t)
    u, gx, gy, lam, dlam, ins = field.eval(S[:n, 0], S[:n, 2], zb, selfsub=selfsub)
    ntot = max(field.nxy, 1e-300)
    c = rcw * al / (b0 * b0 * g0 ** 3) / ntot
    f = np.where(ins, 1.0, 0.0)
    _apply(S, n, -c * lam * gx * f, -c * lam * gy * f, -(rcw * al / g0 ** 2 / ntot) * dlam * u * f, p0, h0)


# ============================================================================ WSPAC: 2 次元のガウス分布の場と電位
def _ein(u):
    """Ein(u) = ∫_0^u (1 − e^(−t))/t dt = E1(u) + ln u + γ"""
    u = np.asarray(u, dtype=float)
    us = np.where(u < 1e-3, 1.0, u)
    return np.where(u < 1e-3, u * (1.0 - u * (0.25 - u / 18.0)), _sp.exp1(us) + np.log(us) + _EULER)


def gauss_field(x, y, sx, sy):
    """Return the field (Fx, Fy) = grad U of a 2D Gaussian with unit line density by the Bassetti-Erskine formula.
    2 次元のガウス分布 (σx = sx、σy = sy、単位の線密度) の場 (Fx, Fy) = ∇U。遠方で 2 (x, y)/r²、中心で
    Fx = 2x/(σx(σx+σy))。Bassetti-Erskine の式 (複素誤差関数 w)。円 (|σx−σy| ≤ 1e-6 σx) は 2 (1 − e^(−r²/2σ²))/r²"""
    if sx < sy:
        fy, fx = gauss_field(y, x, sy, sx)
        return fx, fy
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if sx - sy <= 1e-6 * sx:
        u = (x * x + y * y) / (2.0 * sx * sy)
        us = np.where(u < 1e-8, 1.0, u)
        f = np.where(u < 1e-8, 1.0 - 0.5 * u, -np.expm1(-us) / us) / (sx * sy)
        return f * x, f * y
    ax, ay = np.abs(x), np.abs(y)
    sq = math.sqrt(2.0 * (sx * sx - sy * sy))
    g = np.exp(-0.5 * ((ax / sx) ** 2 + (ay / sy) ** 2))
    w = (2.0 * math.sqrt(math.pi) / sq) * (_sp.wofz((ax + 1j * ay) / sq) - g * _sp.wofz((ax * (sy / sx) + 1j * ay * (sx / sy)) / sq))
    return np.where(ax == 0.0, 0.0, np.copysign(w.imag, x)), np.where(ay == 0.0, 0.0, np.copysign(w.real, y))


_GAUSS64 = np.polynomial.legendre.leggauss(64)


def gauss_potential(x, y, sx, sy):
    """Return the potential U(x, y) = phi(0,0) - phi(x,y) >= 0 (grad U = gauss_field) by Gauss-Legendre integration.
    U(x, y) = φ(0,0) − φ(x,y) ≥ 0 (∇U = gauss_field)。
    U = ∫_0^∞ [1 − exp(−x²/(2σx²+t) − y²/(2σy²+t))] / sqrt((2σx²+t)(2σy²+t)) dt。
    s = sqrt((2σy²+t)/(2σx²+t)) と置くと U = ∫_r^1 2 (1 − e^(−E))/(1 − s²) ds、E = (x² + y²/s²)(1 − s²)/(2(σx²−σy²))、r = σy/σx。
    さらに 1 − s = (1 − r) e^(−q) と置き、q を Gauss-Legendre (64 点) で積分する (q > Q の残りは 1 次の式)。円は Ein(r²/2σ²)"""
    if sx < sy:
        return gauss_potential(y, x, sy, sx)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if sx - sy <= 1e-6 * sx:
        return _ein((x * x + y * y) / (2.0 * sx * sy))
    r = sy / sx
    d = sx * sx - sy * sy
    x2, y2 = (x * x).ravel(), (y * y).ravel()
    t, wt = _GAUSS64
    t = 0.5 * (t + 1.0)
    sm = max(r, 0.5)                              # s = sm..1: 1 − s = (1 − sm) e^(−q)
    qq = np.log(np.maximum((x2 + y2 / (sm * sm)) * (1.0 - sm) / d, 1.0)) + 21.0
    w = (1.0 - sm) * np.exp(-qq[:, None] * t[None, :])
    s = 1.0 - w
    e = (x2[:, None] + y2[:, None] / (s * s)) * w * (2.0 - w) / (2.0 * d)
    u = 0.5 * qq * ((-np.expm1(-e) * 2.0 / (2.0 - w)) @ wt) + (x2 + y2) * (1.0 - sm) * np.exp(-qq) / d
    if r < sm:                                    # s = r..sm (扁平なビーム): s = r e^p (y²/s² の変化が p で緩やかになる)
        pm = math.log(sm / r)
        s = r * np.exp(pm * t)[None, :]
        e = (x2[:, None] + y2[:, None] / (s * s)) * (1.0 - s * s) / (2.0 * d)
        u = u + 0.5 * pm * ((-np.expm1(-e) * 2.0 * s / (1.0 - s * s)) @ wt)
    return u.reshape(x.shape)


class WSpace:
    """Gaussian beam-matrix space-charge (WSPAC, twspac of tintrb.f); rwall sets the reference of the longitudinal-force potential.
    WSPAC (tintrb.f twspac)。mats: 要素 1..nlat の入口のビーム行列 (6×6、(x, px, y, py, z, δ))、cod: (nlat, 6) の閉軌道。
    rcn = r_c·PBUNCH。rwall: 縦の力の電位の基準 (真空ダクトの半径。0 なら K64 と同じ基準: 軸上で ln((1 + r)/2)、r = σ2/σ1)"""

    def __init__(self, mats, cod, rcn, selfcod=False, rwall=0.0):
        m = len(mats)
        self.sig1 = np.zeros(m)
        self.sig2 = np.zeros(m)
        self.c1 = np.ones(m)
        self.s1 = np.zeros(m)
        self.sigz2 = np.zeros(m)
        self.cod = np.asarray(cod, dtype=float).reshape(m, 6)
        self.rcn = rcn
        self.selfcod = bool(selfcod)
        self.rwall = float(rwall or 0.0)
        self.nkick = 0
        for i, b in enumerate(mats):
            b = np.asarray(b, dtype=float)
            zz = b[4, 4]
            if not (np.isfinite(b).all() and zz > 0.0):
                raise SpaceChargeError('WSPAC: 要素 %d の入口のビーム行列が使えない (σz² = %g)。Emittance の平衡のビーム行列が'
                                       '求まらないとき (放射のないリング、輸送路) は WSPACBEAM = {εx, εy, σz, σδ} を与えること' % (i + 1, zz))
            xx = b[0, 0] - b[0, 4] ** 2 / zz               # z で条件を付けた横の分散 (twspac の xx1、yy1、xy1)
            yy = b[2, 2] - b[2, 4] ** 2 / zz
            xy = b[0, 2] - b[0, 4] * b[2, 4] / zz
            u, v = xx - yy, 2.0 * xy
            a = math.hypot(u, v)
            th = 0.5 * math.atan2(v, u)
            s2 = 0.5 * (xx + yy) - 0.5 * a
            if not (s2 > 1e-12 * (0.5 * (xx + yy) + 0.5 * a)):
                raise SpaceChargeError('WSPAC: 要素 %d の入口のビームの大きさが 0 以下 (σ1² = %g、σ2² = %g)。垂直のエミッタンスが 0 の'
                                       'リングでは、結合を入れるか WSPACBEAM = {εx, εy, σz, σδ} を与えること'
                                       % (i + 1, 0.5 * (xx + yy) + 0.5 * a, s2))
            self.sig1[i] = math.sqrt(0.5 * (xx + yy) + 0.5 * a)
            self.sig2[i] = math.sqrt(s2)
            self.c1[i], self.s1[i] = math.cos(th), math.sin(th)
            self.sigz2[i] = zz

    def linear_k(self, l, al, p0):
        """Return the central linear lens strengths kx, ky (as in twspace), with dpx = kx x1 and dpy = ky y1 in principal-axis coordinates.
        中心の線形のレンズの強さ (twspace と同じ): Δpx = kx·x1、Δpy = ky·y1 (主軸の座標)"""
        cd = self.cod[l]
        p1 = p0 * (1.0 + cd[5])
        f = 2.0 * self.rcn * al / (p1 * p1 * _p2h(p1) * (self.sig1[l] + self.sig2[l]) * math.sqrt(2.0 * math.pi * self.sigz2[l]))
        return f / self.sig1[l], f / self.sig2[l]

    def kick(self, S, n, l, al, p0, h0):
        """Compute the kick of length al using the beam matrix at the entrance of element l (0-based).
        要素 l (0 始まり) の入口のビーム行列で、長さ al の分のキック"""
        if al == 0.0 or n <= 0 or self.rcn == 0.0:
            return
        cd = self.cod[l]
        if self.selfcod:
            xc, yc, zc = float(S[:n, 0].mean()), float(S[:n, 2].mean()), float(S[:n, 4].mean())
        else:
            xc, yc, zc = cd[0], cd[2], cd[4]
        c1, s1, sg1, sg2, sz2 = self.c1[l], self.s1[l], self.sig1[l], self.sig2[l], self.sigz2[l]
        p1 = p0 * (1.0 + cd[5])
        f = self.rcn * al / (p1 * p1 * _p2h(p1) * math.sqrt(2.0 * math.pi * sz2))
        b1 = p1 / _p2h(p1)
        p = p0 * (1.0 + S[:n, 5])
        dx = S[:n, 0] - xc
        dy = S[:n, 2] - yc
        # バンチの速度での位置 (= −β1 c·到着時刻)。K64 は粒子の z (= −β_i c t) をそのまま使う (差は δ/γ² の程度)。
        # こうすると、キックは (x, px, y, py, −c t, E) の正準変数で厳密に電位の勾配になる
        dz = b1 * S[:n, 4] / (p / _p2h(p)) - zc
        x1 = dx * c1 + dy * s1
        y1 = -dx * s1 + dy * c1
        az = f * np.exp(-0.5 * dz * dz / sz2)
        fx, fy = gauss_field(x1, y1, sg1, sg2)
        # 縦の力の電位の基準: K64 (twspfu/twspu) の電位は、軸上で ln((1 + r)/2) (r = σ2/σ1)、遠方で丸いビームの形に合わせた値
        # (gs で確認: 軸上の粒子も縦のキックを受け、−z/σz² az·ln((1+r)/2) と一致)。U (軸上 0) に ln((1+r)/2) を足す。
        # rwall > 0 のときは基準を壁の半径にする (K64 にないオプション)
        uw = -math.log(0.5 + 0.5 * sg2 / sg1)
        if self.rwall > 0.0:
            uw = 0.5 * float(gauss_potential(np.array([self.rwall, 0.0]), np.array([0.0, self.rwall]), sg1, sg2).sum())
        dd = dz / sz2 * az * (uw - gauss_potential(x1, y1, sg1, sg2))           # Δδ (Δγ = p0 β1 Δδ)
        _apply(S, n, az * (fx * c1 - fy * s1), az * (fx * s1 + fy * c1), dd * p0 * b1, p0, h0)
        self.nkick += 1


# ============================================================================ 追跡 (track6n.turn) とのつなぎ
K_DRIFT = 0
K_MULT = 10
K_SOL_DRIFT = 12


class SpaceCharge:
    """Space-charge object called by track6n.turn; els is the full element list and types the SAD type name of each element.
    track6n.turn が呼ぶ。els: 追跡する要素の列 (全体)、types: 要素ごとの SAD の型名 ('DRIFT', 'SPCH', …)"""

    def __init__(self, els, types, pbunch, np0, charge=1.0, amass=ELMASS, spac=False, pspac=False, conf=None, warn=None,
                 wspace=None, k64=False):
        conf = conf or {}
        self.k64 = None                          # k64 モードの PSPAC (K64 の形状関数。spch_k64.K64Pspac)
        self.spac, self.pspac = bool(spac), bool(pspac)
        self.wspace = wspace                     # WSPAC (WSpace)。None なら入れない
        self.wspac = wspace is not None
        self.ctx = (False, False, None, 1.0)     # (rad, fluc, rs, tm): track6n.turn が周の始めに入れる (MULT の中の放射用)
        self.s_wlast = 0.0
        self.rcw = charge ** 2 * ELRADI * ELMASS / amass * pbunch / max(np0, 1)
        self.warn = warn
        nel = len(els)
        self.pos = np.zeros(nel + 1)
        self.radius = np.zeros(nel)
        self.is_spch = np.zeros(nel, dtype=bool)
        self.is_mult = np.zeros(nel, dtype=bool)
        self.mult_skip = []                      # RADIUS≠0 だが空間電荷を入れられない MULT の名前
        for i, (el, tp) in enumerate(zip(els, types)):
            self.pos[i + 1] = self.pos[i] + float(getattr(el, 'L', 0.0) or 0.0)
            if tp == 'DRIFT':
                self.radius[i] = element_radius(el)
            elif tp == 'MULT' and element_radius(el) != 0.0:
                if mult_unsupported(el):
                    if el.name not in self.mult_skip:
                        self.mult_skip.append(el.name)
                else:
                    self.radius[i] = element_radius(el)
                    self.is_mult[i] = True
            elif tp == 'SPCH':
                self.is_spch[i] = True
        self.alstep = float(conf.get('SPACSTEP', ALSTEP) or ALSTEP)

        def iv(key):
            v = conf.get(key)
            return int(v) if v else None
        self.s3 = Solver3D(iv('SPACNX'), iv('SPACNY'), iv('SPACNZ'), warn=warn)
        self.s2 = None
        if self.pspac and k64:
            from .spch_k64 import K64Pspac
            self.k64 = K64Pspac(conf.get('PSPACNX', 128), conf.get('PSPACNY', 128), conf.get('PSPACNZ', 1),
                                conf.get('PSPACDX', 1e-3), conf.get('PSPACDY', 1e-3), conf.get('PSPACDZ', 1e-3), pbunch, amass,
                                nturnend=conf.get('PSPACNTURN', 1) or 1, nturncalc=conf.get('PSPACNTURNCALC', 0) or 0, warn=warn)
        elif self.pspac:
            self.s2 = Solver2D(conf.get('PSPACNX', 128), conf.get('PSPACNY', 128), conf.get('PSPACNZ', 1),
                               conf.get('PSPACDX', 1e-3), conf.get('PSPACDY', 1e-3), conf.get('PSPACDZ', 1e-3), warn=warn)
        self.nturncalc = int(conf.get('PSPACNTURNCALC', 0) or 0)
        self._frozen = {}
        self._turn0 = None
        self.s_last = 0.0
        self.l0 = 0
        self.nkick = 0

    # ---- track6n.turn から
    def begin(self, prog, apert):
        """Start a tracking run with the program prog: check that the element table matches, and initialise the start index l0 and the last kick positions.
        追跡プログラム prog の実行を始める。要素の表が追跡の表と合うか確かめ、開始の番号 l0 と直前のキックの位置を初期化する。"""
        if prog.l0 + len(prog.kinds) > len(self.radius):
            raise SpaceChargeError('SPAC/PSPAC: 要素の表が追跡の表と合わない')
        self.l0 = prog.l0
        self.s_last = self.pos[prog.l0]
        self.s_wlast = self.pos[prog.l0]
        if apert is not None and self._turn0 is None:
            self._turn0 = apert[2]

    def wkick(self, S, n, l, ph):
        """Apply the WSPAC kick at the entrance of element l, from the previous element midpoint to this element midpoint.
        WSPAC: 要素 l の入口で、前の要素の中点 (区間の始め) からこの要素の中点までの長さの分のキック (tturn.f:270-278)"""
        s = 0.5 * (self.pos[l] + self.pos[l + 1])
        if s != self.s_wlast:
            self.wspace.kick(S, n, l, s - self.s_wlast, ph.p0, ph.h0)
            self.s_wlast = s

    def handles(self, kind, l):
        """True if the element l of kind kind is treated by space charge (MULT with SPAC, a drift with a radius, or an element with PSPAC).
        種類 kind の要素 l を空間電荷で扱うなら True (SPAC の MULT、半径のあるドリフト、PSPAC の要素)。"""
        if kind == K_MULT:
            return self.spac and self.is_mult[l]
        if kind == K_DRIFT or kind == K_SOL_DRIFT:
            if self.spac and self.radius[l] != 0.0:
                return True
            if self.pspac and self.is_spch[l]:
                return True
        return False

    def element(self, S, n, kind, pr, l, ph, apert, kdrift, ksol):
        """Track the n particles of S through element l with space-charge kicks (PSPAC 2D kick, MULT, or a drift split into steps with the 3D kick and aperture); returns the new particle count.
        要素 l の中で n 個の粒子 S を空間電荷のキックつきで追跡する (PSPAC の 2 次元キック、MULT、または 3 次元キックと絞りを入れて分割したドリフト)。残った粒子数を返す。"""
        if self.pspac and self.is_spch[l]:
            s = 0.5 * (self.pos[l] + self.pos[l + 1])
            n = self._kick2d(S, n, s - self.s_last, l, ph, apert)
            self.s_last = s
            if n > 0 and pr[0] != 0.0:
                kdrift(S, n, pr[0])
            return n
        R = self.radius[l]
        if kind == K_MULT:
            return self._mult(S, n, pr, R, ph, apert)
        al = pr[0]
        if al == 0.0:
            return self._aperture(S, n, abs(R), apert)
        bz = pr[1] if kind == K_SOL_DRIFT else 0.0
        ndiv = max(1, int(round(abs(al) / self.alstep)))
        if bz != 0.0:
            ndiv = max(ndiv, int(round(abs(bz * al) / 1.5)))
        aln = al / ndiv

        def drift(a):
            if bz != 0.0:
                ksol(S, n, a, bz)
            else:
                kdrift(S, n, a)
        drift(0.5 * aln)
        for i in range(ndiv):
            if i > 0:
                drift(aln)
            n = self._aperture(S, n, abs(R), apert)
            if n <= 0:
                return 0
            if R > 0.0 and n > 1 and self.rcw != 0.0:
                b0, g0, zeta = bunch_frame(S, n, ph.p0)
                f = self.s3.solve(S[:n, 0], S[:n, 2], zeta, R)
                kick3d(S, n, f, aln, self.rcw, ph.p0, ph.h0, b0, g0, zeta)
                self.nkick += 1
        drift(0.5 * aln)
        return self._aperture(S, n, abs(R), apert)

    def _sckick(self, S, n, R, al, ph, apert):
        """半径 |R| の絞りと、R>0 なら長さ al の分の SPAC のキック。戻り: 新しい n"""
        n = self._aperture(S, n, abs(R), apert)
        if n > 1 and R > 0.0 and self.rcw != 0.0:
            b0, g0, zeta = bunch_frame(S, n, ph.p0)
            f = self.s3.solve(S[:n, 0], S[:n, 2], zeta, R)
            kick3d(S, n, f, al, self.rcw, ph.p0, ph.h0, b0, g0, zeta)
            self.nkick += 1
        return n

    def _mult(self, S, n, pr, R, ph, apert):
        """MULT (RADIUS≠0) の中の空間電荷 (tmulti.f の spac1 の分岐の意図)。座標は MULT の中の系 (DX・DY・ROTATE の後) で、
        真空ダクトの軸は磁石の軸。分割は ndiv = max(MULT の分割, nint(|L|/(SPACSTEP·EPS)), 放射の分割)。
        本体の半分 → [多極 (n≥2) のキック、絞り、SPAC のキック (L/ndiv)]×ndiv (間は本体 1 つ分) → 本体の半分 → 絞り。
        K64 は最後の spkick の長さだけ半分 (al*wm) で、合計が (ndiv − 1/2)/ndiv になる。v2 は ndiv 回とも L/ndiv"""
        from . import track6n as _n, track6n_ext as _x
        rad, fluc, rs, tm = self.ctx
        (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad, radok, _rf) = pr
        if al == 0.0:                            # 薄い MULT: 空間電荷は長さ 0。絞りだけ
            _x.k_mult(S, n, al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad, False, fluc, rs,
                      ph.p0, ph.h0 * tm, ph.anrad, ph.cuc)
            return self._aperture(S, n, abs(R), apert)
        krad = bool(rad and radok)
        p0, anrad, cuc = ph.p0, ph.anrad, ph.cuc
        akr1 = float(akr[1].real) if nmmax >= 1 else 0.0
        akr0 = complex(akr[0])
        epsr = eps_rad if eps_rad != 0.0 else 1.0
        ndiv = max(int(ndiv0), int(round(abs(al) / (self.alstep * abs(epsr)))), 1)
        _x.mult_rot_in(S, n, dx, dy, ct, st)
        if krad:
            S[:n, 7] = S[:n, 1]
            S[:n, 8] = S[:n, 3]
            S[:n, 9] = S[:n, 4]
            ndiv = max(ndiv, _n.ndivrad(abs(akr0), akr1, 0.0, eps_rad, ph.h0 * tm))
        ndiv = min(1000, ndiv)
        if fringe:
            _x.mult_fringe_all(S, n, akr, dofr, nmmax, al, 1.0, False)
        wi = 1.0 / ndiv
        akrm = akr * wi
        aln = al * wi

        def body(w):
            _x.mult_qu_all(S, n, al * w, akr1 * w, 0.0, akr0.real * w, akr0.imag * w, _qb.sch_of(rs))
            if krad:
                _n.tradk_all(S, n, al * w, 0.0, fluc, rs, p0, anrad, cuc)
        body(0.5 * wi)
        for m in range(ndiv):
            if m > 0:
                body(wi)
            if nmmax >= 2:
                _x.mult_kick_all(S, n, akrm, nmmax, 2)
            n = self._sckick(S, n, R, aln, ph, apert)
            if n <= 0:
                return 0
        _x.mult_qu_all(S, n, 0.5 * aln, 0.5 * akr1 * wi, 0.0, 0.5 * akr0.real * wi, 0.5 * akr0.imag * wi, _qb.sch_of(rs))
        if fringe:
            _x.mult_fringe_all(S, n, akr, dofr, nmmax, al, -1.0, True)
        if krad:                                 # k_mult と同じ順 (本体 → 出口のフリンジ → 放射)
            _n.tradk_all(S, n, 0.5 * aln, 0.0, fluc, rs, p0, anrad, cuc)
        n = self._aperture(S, n, abs(R), apert)
        if n > 0:
            _x.mult_rot_out(S, n, dx, dy, ct, st)
        return n

    def end(self, S, n, l_end, ph, apert):
        """At the section end, apply the PSPAC kick from the last SPCH, or the WSPAC kick from the last element midpoint, to the endpoint.
        区間の終わり: PSPAC は最後の SPCH から終点までの分のキック、WSPAC は最後の要素の中点から終点までの分"""
        if self.wspac and n > 0:
            s = self.pos[l_end]
            if s != self.s_wlast:
                self.wspace.kick(S, n, l_end, s - self.s_wlast, ph.p0, ph.h0)
                self.s_wlast = s
        if self.pspac and n > 0:
            s = self.pos[l_end]
            if s != self.s_last and self.is_spch[self.l0:l_end].any():
                n = self._kick2d(S, n, s - self.s_last, -1, ph, apert)
            self.s_last = s
        if self.k64 is not None:
            self.k64.end_of_turn()               # tpspac_: nlat == nlatend で isp = 0、nturn++
        return n

    # ---- 中身
    def _kick2d(self, S, n, al, l, ph, apert):
        if self.k64 is not None:                 # k64: K64 の形状関数 (spch_k64)。al = 0 でも電位を計算する (isp が進む)
            return self.k64.kick(S, n, al, ph)
        if al == 0.0 or n <= 1 or self.rcw == 0.0:
            return n
        b0, g0, _ = bunch_frame(S, n, ph.p0)
        f = None
        own = True
        frozen = (self.nturncalc > 0 and apert is not None and self._turn0 is not None
                  and apert[2] - self._turn0 >= self.nturncalc)
        if frozen:
            f = self._frozen.get(l)
            own = f is None
        if f is None:
            p = ph.p0 * (1.0 + S[:n, 5])
            zb = b0 * S[:n, 4] / (p / _p2h(p))
            f = self.s2.solve(S[:n, 0], S[:n, 2], zb)
            if self.nturncalc > 0:
                self._frozen[l] = f
        kick2d(S, n, f, al, self.rcw, ph.p0, ph.h0, b0, g0, selfsub=own)
        self.nkick += 1
        return n

    @staticmethod
    def _aperture(S, n, R, apert):
        """半径 R の円の絞り (spkick.f の tapert(radius, radius))。apert = (kid, lost_turn, kturn, alost, zlost, …)"""
        if apert is None or R <= 0.0 or n <= 0:
            return n
        kid, lost_turn, kturn, _, zlost = apert[:5]
        with np.errstate(all='ignore'):
            ok = ((S[:n, 0] ** 2 + S[:n, 2] ** 2 <= R * R) & (np.abs(S[:n, 4]) <= zlost)
                  & np.isfinite(S[:n, 1]) & np.isfinite(S[:n, 3]) & np.isfinite(S[:n, 5]))       # NaN の粒子も落とす
        if ok.all():
            return n
        drop = ~ok
        for i in np.nonzero(drop)[0]:
            lost_turn[kid[i]] = kturn
        from .track6n import drop_particles
        return drop_particles(S, n, kid, drop)


def element_radius(el):
    """Return the RADIUS of a DRIFT, taking an attribute changed in FFS before the defined value.
    DRIFT の RADIUS (FFS で変えた属性 → 定義の値 の順)"""
    v = getattr(el, 'RADIUS', None)
    if v is None:
        v = (getattr(el, '_raw', None) or {}).get('RADIUS', 0.0)
    try:
        return float(v or 0.0)
    except (TypeError, ValueError):
        return 0.0


def mult_unsupported(el):
    """Return the reason why space charge cannot be applied to a MULT with RADIUS != 0, or '' if it can.
    RADIUS≠0 の MULT に空間電荷を入れられない理由 (入れられるなら '')"""
    def g(k):
        try:
            return float(getattr(el, k, 0.0) or 0.0)
        except (TypeError, ValueError):
            return 0.0
    if g('ANGLE') != 0.0:
        return 'ANGLE≠0'
    if g('VOLT') + g('DVOLT') != 0.0:
        return 'VOLT≠0'
    return ''


def twiss_beam_matrix(r, emx, emy, sigz, sige):
    """Build the 6x6 beam matrix in physical coordinates from a Twiss row and emittances, as in tfbeam.f.
    Twiss の行 r (ffscmd の row) とエミッタンスから、物理座標のビーム行列 (6×6)。tfbeam.f (ffscmd.beam_matrix_twiss) と同じ形:
    ノーマルモードの (β, α, 分散) の行列を R1-R4 で物理座標に移す。x-z、y-z の相関は 0、Σzz = σz²"""
    ax, bx, ay, by = r.alfx, r.betx, r.alfy, r.bety
    ex, epx, ey, epy = r.etx, r.etpx, r.ety, r.etpy
    pfi = sige * sige
    B = np.zeros((6, 6))
    B[0, 0] = bx * emx + ex * ex * pfi
    B[0, 1] = B[1, 0] = -ax * emx + ex * epx * pfi
    B[0, 2] = B[2, 0] = ex * ey * pfi
    B[0, 3] = B[3, 0] = ex * epy * pfi
    B[1, 1] = (1.0 + ax * ax) / bx * emx + epx * epx * pfi
    B[1, 2] = B[2, 1] = epx * ey * pfi
    B[1, 3] = B[3, 1] = epx * epy * pfi
    B[2, 2] = by * emy + ey * ey * pfi
    B[2, 3] = B[3, 2] = -ay * emy + ey * epy * pfi
    B[3, 3] = (1.0 + ay * ay) / by * emy + epy * epy * pfi
    B[0, 5] = B[5, 0] = ex * pfi
    B[1, 5] = B[5, 1] = epx * pfi
    B[2, 5] = B[5, 2] = ey * pfi
    B[3, 5] = B[5, 3] = epy * pfi
    B[5, 5] = pfi
    r1, r2, r3, r4 = r.r1, r.r2, r.r3, r.r4
    if r1 != 0.0 or r2 != 0.0 or r3 != 0.0 or r4 != 0.0:
        cc = math.sqrt(max(0.0, 1.0 - (r1 * r4 - r2 * r3)))
        T = np.eye(6)
        if r.detr < 1.0:
            T[0:4, 0:4] = [[cc, 0.0, r4, -r2], [0.0, cc, -r3, r1], [-r1, -r2, cc, 0.0], [-r3, -r4, 0.0, cc]]
        else:
            T[0:4, 0:4] = [[r4, -r2, cc, 0.0], [-r3, r1, 0.0, cc], [cc, 0.0, -r1, -r2], [-r3, -r4, 0.0, cc]]
        B = T @ B @ T.T
    B[4, 4] = sigz * sigz
    return B


def make_wspace(ffs, tr, g):
    """Set up WSPAC with a beam matrix (from WSPACBEAM or the Emittance equilibrium) and the closed orbit.
    WSPAC の設定: ビーム行列 (WSPACBEAM = {εx, εy, σz, σδ} があれば Twiss から、無ければ Emittance の平衡のビーム行列) と閉軌道"""
    nlat = len(ffs.els) + 1
    wb = g.get('WSPACBEAM')
    rows = [ffs.row(k) for k in range(1, nlat + 1)]
    if isinstance(wb, (list, tuple)):
        try:
            v = [float(t) for t in wb]
        except (TypeError, ValueError):
            v = []
        if len(v) not in (3, 4) or not (v[0] > 0.0 and v[1] > 0.0 and v[2] > 0.0):
            raise SpaceChargeError('WSPAC: WSPACBEAM = {εx, εy, σz} か {εx, εy, σz, σδ} (εx、εy、σz は正の数) で与えること')
        mats = [twiss_beam_matrix(r, v[0], v[1], v[2], v[3] if len(v) == 4 else 0.0) for r in rows]
    else:
        try:
            mats = [ffs.beam_matrix_at(k) for k in range(1, nlat + 1)]
        except SADError:
            raise
        except Exception as e:
            raise SpaceChargeError('WSPAC: Emittance のビーム行列が求まらない (%s)。WSPACBEAM = {εx, εy, σz, σδ} を与えること' % (e,))
    cod = [[r.x_cod, r.px_cod, r.y_cod, r.py_cod, r.z_cod, r.delta_cod] for r in rows]
    rw = g.get('WSPACRADIUS')
    rw = float(rw) if isinstance(rw, (int, float)) and not isinstance(rw, bool) else 0.0
    rcn = tr.ph.charge ** 2 * ELRADI * ELMASS / tr.ph.amass * float(g.get('PBUNCH', 0.0) or 0.0)
    return WSpace(mats, cod, rcn, selfcod=bool(ffs.flags.get('SELFCOD')), rwall=rw)


def make(ffs, tr, npart):
    """Create a SpaceCharge from `TrackParticles` in ffscmd if flag SPAC, PSPAC or WSPAC is ON, else None (v2 mode only).
    ffscmd の TrackParticles から: フラグ SPAC / PSPAC / WSPAC が ON なら SpaceCharge を作って返す (無ければ None)。v2 モードだけ"""
    spac = bool(ffs.flags.get('SPAC'))
    pspac = bool(ffs.flags.get('PSPAC'))
    wspac = bool(ffs.flags.get('WSPAC'))
    if not (spac or pspac or wspac):
        return None
    from . import compat as _cmp
    if _cmp.k64():
        # k64: K64 が追跡で実際に入れるのは PSPAC (SPCH 要素) だけ (SPAC は力が 1/γ 倍に弱い B89、WSPAC の追跡は未対応で今までどおり入れない)。
        # 形状関数は K64 のまま (spch_k64)
        if not pspac or tr is None:
            return None
        return _make_k64(ffs, tr, npart)
    els = ffs.els
    env = getattr(ffs, 'env', None)
    defs = getattr(env, '_elements', {}) if env is not None else {}
    types = []
    for el in els:
        tp = str((defs.get(el.name) or {}).get('TYPE', '')).upper()
        if not tp:
            tp = 'DRIFT' if type(el).__name__ == 'Drift' and not getattr(el, '_apert', False) else type(el).__name__.upper()
        types.append(tp)
    # SOL 領域 (BOUND≠0 の SOL の間) の中の MULT は別のカーネル (K_SOL_MULT) で追うので、空間電荷を入れられない
    insol = False
    solmult = []
    for el, tp in zip(els, types):
        if tp == 'SOL' and float(getattr(el, 'BOUND', 0.0) or 0.0) != 0.0:
            insol = not insol
        elif insol and tp == 'MULT' and element_radius(el) != 0.0 and el.name not in solmult:
            solmult.append(el.name)

    def scr(el, tp):
        return (tp == 'DRIFT' or (tp == 'MULT' and el.name not in solmult and not mult_unsupported(el))) and element_radius(el) != 0.0
    has_r = spac and any(scr(el, tp) for el, tp in zip(els, types))
    has_s = pspac and any(tp == 'SPCH' for tp in types)
    if spac:
        skip = [(nm, 'SOL 領域の中') for nm in solmult]
        for el, tp in zip(els, types):
            why = mult_unsupported(el) if tp == 'MULT' and element_radius(el) != 0.0 and el.name not in solmult else ''
            if why and (el.name, why) not in skip:
                skip.append((el.name, why))
        if skip:
            ffs.write(' SPAC: 次の MULT は RADIUS≠0 だが空間電荷を入れない (ANGLE≠0・VOLT≠0・SOL 領域の中の MULT は未対応): '
                      + ', '.join('%s (%s)' % t for t in skip) + '\n')
    if not (has_r or has_s or wspac):
        return None
    if tr is None:
        raise SpaceChargeError('SPAC/PSPAC/WSPAC: この要素の並びは高速追跡 (TrackRun) で扱えないため、空間電荷を入れられない')
    needs_z = has_s or (spac and any(scr(el, tp) and element_radius(el) > 0.0 for el, tp in zip(els, types)))
    if needs_z and not ffs.flags.get('RFSW'):
        raise SpaceChargeError('SPAC/PSPAC: RFSW を ON にすること。NORFSW は区間 (周) ごとに粒子の z を 0 にするので、'
                               'バンチの線密度が決まらない (空洞が無ければ RFSW ON でも軌道は変わらない)')
    g = ffs.I.globals
    conf = {}
    for key in ('SPACNX', 'SPACNY', 'SPACNZ', 'SPACSTEP', 'PSPACNX', 'PSPACNY', 'PSPACNZ', 'PSPACDX', 'PSPACDY', 'PSPACDZ',
                'PSPACNTURNCALC'):
        v = g.get(key)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            conf[key] = float(v)
    wsp = make_wspace(ffs, tr, g) if wspac else None
    if wspac and (has_r or has_s):
        ffs.write(' WSPAC: SPAC / PSPAC も ON なので、空間電荷の力が二重に入る\n')
    sc = SpaceCharge(els, types, float(g.get('PBUNCH', 0.0) or 0.0), npart, charge=tr.ph.charge, amass=tr.ph.amass,
                     spac=spac, pspac=pspac, conf=conf, warn=lambda m: ffs.write(' ' + m + '\n'), wspace=wsp)
    if spac:
        for i, (el, tp) in enumerate(zip(els, types)):       # SOL 領域の中の MULT は空間電荷なしで追う
            if tp == 'MULT' and el.name in solmult:
                sc.radius[i] = 0.0
                sc.is_mult[i] = False
    return sc


def _make_k64(ffs, tr, npart):
    """k64 モードの PSPAC (SPCH 要素)。SPAC・WSPAC は入れない"""
    els = ffs.els
    env = getattr(ffs, 'env', None)
    defs = getattr(env, '_elements', {}) if env is not None else {}
    types = []
    for el in els:
        tp = str((defs.get(el.name) or {}).get('TYPE', '')).upper()
        if not tp:
            tp = 'DRIFT' if type(el).__name__ == 'Drift' and not getattr(el, '_apert', False) else type(el).__name__.upper()
        types.append(tp)
    if not any(tp == 'SPCH' for tp in types):
        return None
    g = ffs.I.globals
    conf = {}
    for key in ('PSPACNX', 'PSPACNY', 'PSPACNZ', 'PSPACDX', 'PSPACDY', 'PSPACDZ', 'PSPACNTURN', 'PSPACNTURNCALC'):
        v = g.get(key)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            conf[key] = float(v)
    return SpaceCharge(els, types, float(g.get('PBUNCH', 0.0) or 0.0), npart, charge=tr.ph.charge, amass=tr.ph.amass,
                       spac=False, pspac=True, conf=conf, warn=lambda m: ffs.write(' ' + m + chr(10)), k64=True)
