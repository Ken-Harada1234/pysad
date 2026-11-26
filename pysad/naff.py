# -*- coding: utf-8 -*-
"""NAFF-style tune extraction (Laskar numerical frequency analysis) from tracking coordinate sequences (v2next B3).
NAFF 風のチューン抽出 (Laskar の数値的な周波数解析)。追跡の座標列からチューンを取り出す Python の関数 (v2next B3、項目 5)。

K64 の MAIN の TRACK には wtunem (FFT のピークの放物線補間、`maintrack.wtunem`) と、FOURIER フラグの wjfit (wfres・wfbin・wiord、
作用変数の Fourier のフィット) がある。どちらも FFS から呼べず、精度は FFT の分解能 (1/N) の数分の 1 から 1/N² 程度。
ここでは次を実装する (K64 の直訳ではなく、同じ目的の独立な実装)。

  naff(z, nfreq)           複素信号 z_n = Σ a_k exp(2πi ν_k n) の周波数 ν_k と複素振幅 a_k を大きい順に nfreq 個
                           (Hann 窓つき、FFT のピークをゼロ詰めで探し、窓つきの内積 |<z, e^{iνn}>| の最大化 (Brent 法) で精密化、
                           見つけた成分を Gram-Schmidt で引いて次へ)。窓の次数 window=1 で誤差 ~1/N^4 (Laskar)。
  naff_real(x, nfreq)      実数の信号 (x の列だけ): 解析信号 (正の周波数側だけ) にして naff。周波数は 0〜0.5
  tune_xp(x, px, alpha, beta)   (x, px) の列と Twiss から z = x̂ − i p̂ を作り、符号つきのチューン ν (-0.5〜0.5) を返す
  tunes_from_coords(coords, ...)  track_particles の座標列 (ターン × 6) から (νx, νy, νz) を返す。x, y, z 面は (px, py, dp) か
                           位相空間 z−dp を使う。αβ を省くと z = x − i px (αβ の無いリングで近似)

周波数の単位は 1/周回 (チューンの端数)。符号は exp(+2πiνn) の向き: z = x̂ − i p̂ (x̂ = a cos ψ, p̂ = −a sin ψ で ψ = 2πνn) が +ν。
SAD の追跡の座標は px = dx/ds なので、αβ を使って x̂ = x/√β、p̂ = (α x + β px)/√β に直す (tune_xp)。
限界: 周波数が近い成分 (間隔 < 約 2/N) は分離できない。共鳴線の近く (ν の差が 1/N 以下) では振幅が不安定。減衰のある信号 (放射減衰) は
exp(−λ n) の包絡で周波数が 1/(N) の程度ずれる (放射を切って使うこと)。
"""
import math
import numpy as np


def _hann(n, p=1):
    """Laskar の窓: χ_p(t) = 2^p (p!)² / (2p)! (1 + cos πt)^p、t ∈ [−1, 1]。p=1 で (1 + cos) = 2 cos²。和は N に規格化"""
    t = (np.arange(n) + 0.5) / n * 2.0 - 1.0         # n 点の中心 (半整数の点)
    w = (1.0 + np.cos(math.pi * t)) ** p
    return w / w.mean()


def _proj(z, w, nu):
    """窓つきの内積 <z, e^{2πiνn}> = (1/N) Σ w_n z_n e^{−2πiνn}"""
    n = np.arange(len(z))
    return np.sum(w * z * np.exp(-2j * math.pi * nu * n)) / len(z)


def _refine(z, w, nu0, width):
    """|<z, e^{iνn}>| の最大を [nu0 − width, nu0 + width] で Brent 法 (黄金分割 + 放物線) により探す"""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda nu: -abs(_proj(z, w, nu)), bounds=(nu0 - width, nu0 + width), method='bounded',
                          options={'xatol': 1e-14})
    return float(res.x)


def naff(z, nfreq=1, window=1, fmin=-0.5, fmax=0.5, pad=8, tol=1e-14):
    """Compute the NAFF frequencies and complex amplitudes of a complex signal z, returned in order of decreasing amplitude.
    複素信号 z の NAFF。戻り (freqs, amps): freqs は 1/周回 (大きい振幅の順)、amps は exp(2πiνn) の複素振幅。
    fmin・fmax で探す周波数の範囲 (既定は −0.5〜0.5)。amps は、見つけた周波数を固定した窓つきの最小二乗で求める
    (Gram-Schmidt で直交化した基底への射影と同じ)"""
    z = np.asarray(z, dtype=complex)
    nn = len(z)
    if nn < 8:
        raise ValueError('NAFF: 点が少なすぎる')
    w = _hann(nn, window)
    n = np.arange(nn)
    res = z.copy()
    freqs = []
    nfft = int(2 ** math.ceil(math.log2(nn * pad)))
    fgrid = np.fft.fftfreq(nfft)
    ok = (fgrid >= fmin - 1e-15) & (fgrid <= fmax + 1e-15)
    basis = []          # 直交化した基底 (窓つきの内積で)
    for _ in range(nfreq):
        spec = np.fft.fft(w * res, nfft)
        mag = np.abs(spec)
        mag[~ok] = 0.0
        k = int(np.argmax(mag))
        nu0 = float(fgrid[k])
        nu = _refine(res, w, nu0, 1.5 / nn)
        # 範囲外に出たら端に押し戻す
        nu = min(max(nu, fmin), fmax)
        freqs.append(nu)
        e = np.exp(2j * math.pi * nu * n)
        for b in basis:                        # Gram-Schmidt
            e = e - (np.sum(w * e * np.conj(b)) / nn) * b
        e = e / math.sqrt(np.real(np.sum(w * e * np.conj(e))) / nn)
        basis.append(e)
        res = res - (np.sum(w * res * np.conj(e)) / nn) * e
    freqs = np.array(freqs)
    E = np.exp(2j * math.pi * np.outer(n, freqs))
    sw = np.sqrt(w)[:, None]
    amps = np.linalg.lstsq(E * sw, z * sw[:, 0], rcond=None)[0]
    order = np.argsort(-np.abs(amps))
    return freqs[order], amps[order]


def naff_real(x, nfreq=1, window=1, remove_mean=True, **kw):
    """Compute frequencies, real amplitudes, and phases of a real signal by applying `naff` to its analytic signal.
    実数の信号 x_n = Σ A_k cos(2πν_k n + φ_k) の周波数 ν_k (0〜0.5) と振幅 A_k (cos の振幅)、位相 φ_k。
    解析信号 (FFT の正の周波数側だけ) にして naff を呼ぶ。戻り (freqs, amps (実数の振幅), phases)"""
    x = np.asarray(x, dtype=float)
    if remove_mean:
        x = x - x.mean()
    nn = len(x)
    X = np.fft.fft(x)
    h = np.zeros(nn)
    h[0] = 1.0
    if nn % 2 == 0:
        h[nn // 2] = 1.0
        h[1:nn // 2] = 2.0
    else:
        h[1:(nn + 1) // 2] = 2.0
    za = np.fft.ifft(X * h)
    f, a = naff(za, nfreq, window=window, fmin=0.0, fmax=0.5, **kw)
    return f, np.abs(a), np.angle(a)


def twiss_normalize(x, px, alpha=0.0, beta=1.0):
    """x̂ = x/√β、p̂ = (α x + β px)/√β"""
    x = np.asarray(x, dtype=float)
    px = np.asarray(px, dtype=float)
    return x / math.sqrt(beta), (alpha * x + beta * px) / math.sqrt(beta)


def tune_xp(x, px, alpha=0.0, beta=1.0, window=1, nfreq=1):
    """Return the signed tune from an (x, px) sequence, optionally with its amplitude, or several frequencies when nfreq > 1.
    (x, px) の列から符号つきのチューン ν (−0.5〜0.5) を返す。z = x̂ − i p̂。SAD の位相の向き (x̂ = a cos ψ, p̂ = −a sin ψ、ψ = 2πνn) で ν>0。
    戻り (ν, 振幅 |a|) (nfreq=1)、または nfreq>1 のとき (freqs, |amps|)"""
    xh, ph = twiss_normalize(x, px, alpha, beta)
    z = (xh - xh.mean()) - 1j * (ph - ph.mean())
    f, a = naff(z, nfreq, window=window)
    if nfreq == 1:
        return float(f[0]), float(abs(a[0]))
    return f, np.abs(a)


def tunes_from_coords(coords, twiss=None, window=1):
    """Return the signed tunes (nu_x, nu_y, nu_z) and per-plane amplitudes from turn-by-turn six-dimensional coordinates.
    coords: ターン × 6 (x, px, y, py, z, dp)。twiss: ((αx, βx), (αy, βy), (αz, βz)) (省略すると α=0、β=1)。
    戻り (νx, νy, νz) (−0.5〜0.5 の符号つき) と各面の振幅 (3 個)。"""
    c = np.asarray(coords, dtype=float)
    tw = twiss if twiss is not None else ((0.0, 1.0),) * 3
    out = []
    amp = []
    for i in range(3):
        f, a = tune_xp(c[:, 2 * i], c[:, 2 * i + 1], tw[i][0], tw[i][1], window=window)
        out.append(f)
        amp.append(a)
    return tuple(out), tuple(amp)
