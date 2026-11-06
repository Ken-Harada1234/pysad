"""
sadrandom.py — SAD (K64) の既定の乱数プラグイン "SAD" の直訳 (src/random_plugin_sad.c)

- 32 bit 符号付き合同法: seed ← 48828125·seed (mod 2^32, int32 で折り返し)。SeedRandom[n] は
  seed = |n/2|·2+1 (奇数化)。FFS 開始時の種は 17 (SADHelp: 「初期値 17」)。
- Random[] / tran(): [0,1) = (seed + 2^31)/2^32
- GaussRandom[] / tgauss(): Box–Muller。φ = seed·π/2^31、ρ = sqrt(−2 log(seed/2^32 + 0.5))。
  ρ cos φ を返し、ρ sin φ を次回用に保持 (filled/next)。|値| > GCUT は棄却。
- SeedRandom[] の戻り {"SAD", seed, filled, next(2 語)} に対応する dump()。

追跡の FLUC (tradkfn / tdusrn / drndsr) はこの一様乱数を引く。SAD と同じ順で引けば同じ数列になる。
"""
from __future__ import annotations

import math

_A = 48828125
_M32 = 1 << 32
_M31 = 1 << 31


def _lcg(seed: int) -> int:
    s = (_A * seed) % _M32
    return s - _M32 if s >= _M31 else s


class SadRandom:
    def __init__(self, seed: int = 17):
        self.restore(seed)

    # random_plugin_sad.c restore_internal_state (length 1)
    def restore(self, seed: int):
        seed = int(seed)
        if seed < 0:
            self.seed = -(abs(seed // 2 if seed % 2 == 0 else -(-seed // 2)) * 2 + 1)
        else:
            self.seed = abs(seed // 2) * 2 + 1
        # C の seed/2 は 0 方向への切り捨て
        q = abs(seed) // 2
        self.seed = -(q * 2 + 1) if seed < 0 else q * 2 + 1
        self.filled = False
        self.next = 0.0

    def dump(self):
        """SeedRandom[] の戻りの数値部分 {seed, filled, next} 相当"""
        return (self.seed, 1 if self.filled else 0, self.next)

    # RANDOM_UNIFORM0: [0, 1)   (tran(), Random[])
    def uniform0(self) -> float:
        self.seed = _lcg(self.seed)
        return (self.seed + 2147483648.0) * (1.0 / 4294967296.0)

    tran = uniform0

    # RANDOM_UNIFORM: (0, 1)
    def uniform(self) -> float:
        self.seed = _lcg(self.seed)
        return (self.seed + 2147483648.5) * (1.0 / 4294967296.0)

    # RANDOM_GAUSS (GCUT 付き)。n 個まとめて生成 (length の扱いが C と同じになるように)
    def gauss(self, n: int = 1, gcut: float = 1.0e300):
        out = []
        i = 0
        while i < n:
            if not self.filled:
                self.seed = _lcg(self.seed)
                x = self.seed * (math.pi / 2147483648.0)
                self.seed = _lcg(self.seed)
                y = math.sqrt(-2.0 * math.log(self.seed * (1.0 / 4294967296.0) + 0.5))
                out.append(y * math.cos(x))
                i += 1
                if abs(out[-1]) > gcut:
                    out.pop()
                    i -= 1
                if i < n:
                    out.append(y * math.sin(x))
                    i += 1
                    if abs(out[-1]) > gcut:
                        out.pop()
                        i -= 1
                else:
                    self.next = y * math.sin(x)
                    self.filled = True
            else:
                out.append(self.next)
                self.filled = False
                self.next = 0.0
                i += 1
                if abs(out[-1]) > gcut:
                    out.pop()
                    i -= 1
        return out if n != 1 else out[0]

    tgauss = gauss

    # RANDOM_PARABOLA (random_plugin_sad.c, USE_PARABOLA_TRANSFORM なし): x∈[-1,1), y∈[0,1) を x² ≤ y まで引き直す
    def parabola(self, n: int = 1):
        out = []
        for _ in range(n):
            while True:
                self.seed = _lcg(self.seed)
                x = self.seed * (1.0 / 2147483648.0)
                self.seed = _lcg(self.seed)
                y = self.seed * (1.0 / 4294967296.0) + 0.5
                if not (x * x > y):
                    break
            out.append(x)
        return out if n != 1 else out[0]


# ---------------------------------------------------------------------------
# drndsr.f: 放射スペクトルの乱数 (横谷の逆関数近似) と光子数 (tdusrn)
# ---------------------------------------------------------------------------

def drndsr(rng: SadRandom) -> float:
    """drndsr: x = (光子エネルギー)/(臨界エネルギー) の乱数。tran() を 1 回引く"""
    p = rng.tran()
    if p < 0.658:
        p2 = p * p
        return (((0.41841680 * p2 + 0.14175015) * p2 + 0.30528148) * p2 + 0.53520052) * p2 * p
    q = 1.0 - p
    if q > 0.0297:
        return (((-0.32809490 * q + 0.20653610) * q + 0.011920005)
                / (((q + 0.88765651) * q + 0.19270658) * q + 0.0033141692))
    t = -math.log(q)
    return (674.99117 * t + 148.32871) / ((t - 225.45670) * t - 692.23986) + t


_CUTOFF = 0.9999
_NPMAX = 10000


def tdusrn(rng: SadRandom, anp: float):
    """tdusrn: 期待光子数 anp のポアソン分布で光子を出し、
    dph = Σ u_n/uc, r1 = Σ r_n u_n / Σ u_n, r2 = Σ r_n² u_n / Σ r_n u_n / 2, an = 光子数。
    乱数は p (1 回) と光子ごとに drndsr (1 回) + r (1 回)。"""
    p0 = math.exp(-anp)
    pn = p0
    p = rng.tran()
    r1 = r2 = an = dph = 0.0
    x = 1.0
    for _ in range(_NPMAX):
        if p <= pn:
            break
        an += 1.0
        x = x * anp / an
        pn = pn + x * p0
        dph1 = drndsr(rng)
        dph += dph1
        r = rng.tran()
        r1 += r * dph1
        r2 += r * r * dph1
        if pn > _CUTOFF:
            break
    if r1 != 0.0:
        r2 = r2 / r1 * 0.5
        r1 = r1 / dph
    return dph, r1, r2, an
