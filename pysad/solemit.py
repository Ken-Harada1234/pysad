"""
SOL 領域 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで) の Emittance[] / CALC6D / EMIT 用の写像。K64 の tsole1 (tsole.f) 相当。

領域内の要素は solmap.annotate が el._solmap = (kind, params, momentum) を付けておく。ここではその kind ごとに
tsole1 → tdrife0 / tquade (tsolrote → tqfrie → tsolque の分割 → tqfrie → tsolrote) / tsfrie の順に、
写像は追跡のカーネル (track6n_ext の tz_step1 など。同じ解析式を翻訳したもの) の 1 粒子の数値微分 (中心差分 1e-7) で、
放射 (tradke) は bz 付きの emit6d._State.tradke で、K64 と同じ位置に入れる。

領域内の座標は正準運動量 (px + bz/2·y が力学的運動量)。tradke の bzh=bzs/2、tsetr0 の bzh も K64 と同じ。

対応: DRIFT (放射は bz≠0 のとき)、QUAD (K1≠0、ずれ・F1・FRINGE なし)、L=0 の要素 (恒等)、SOL (F1 なし。ずれ・回転なし)。
SOL の端は F1≠0 のときだけ放射する (tsole1 の krad) が、F1≠0 は追跡側の compile が拒否するので扱わない。
"""
from __future__ import annotations

import math

import numpy as np

from . import track6n as _n
from . import track6n_ext as _x
from . import emit6d as _e

EPS = 1e-7
EPSRAD = _e.EPSRAD


def _rows(fn, ph, C):
    """C (m, 6) の物理座標 → fn(S, i) を各行に適用して物理座標で返す"""
    C = np.asarray(C, dtype=float).reshape(-1, 6)
    S = _n.make_state(C, ph)
    for i in range(C.shape[0]):
        fn(S, i)
    return _n.state_out(S, C.shape[0])


def _piece(st, fn, symp_z=True):
    """st.piece: 写像 f は 1 行、ヤコビアンは 12 行をまとめて中心差分"""
    ph = st.ph

    def f(c):
        return _rows(fn, ph, c)[0]

    def jac(c):
        C = np.tile(np.asarray(c, dtype=float), (12, 1))
        for j in range(6):
            C[2 * j, j] += EPS
            C[2 * j + 1, j] -= EPS
        out = _rows(fn, ph, C)
        J = np.empty((6, 6))
        for j in range(6):
            J[:, j] = (out[2 * j] - out[2 * j + 1]) / (2.0 * EPS)
        return J
    st.piece(f, symp_z=symp_z, jac=jac)


def _tomech(S, i, bz):
    """正準 → 力学的運動量 (追跡の tsolqu の先頭と同じ: 規格化した px に bz/pr·y/2 を足す)"""
    bzp = bz / (1.0 + S[i, 5])
    S[i, 1] += bzp * S[i, 2] * 0.5
    S[i, 3] -= bzp * S[i, 0] * 0.5


def _tocan(S, i, bz):
    bzp = bz / (1.0 + S[i, 5])
    S[i, 1] -= bzp * S[i, 2] * 0.5
    S[i, 3] += bzp * S[i, 0] * 0.5


def _itgetqraddiv(ph, c, ak, al, bzh):
    """tsolque.f itgetqraddiv (bz あり)。c は正準運動量のままの cod"""
    xk = abs(ak) * max(1e-6, abs(c[0]) + abs(c[2]))
    xpd = max(1e-6, abs(c[1]) + abs(c[3]))
    a = min(1e-2, xk + xpd)
    b = ph.brhoz * (a * abs(bzh) + xk / abs(al))
    nrad = int(abs(al * ph.crad / EPSRAD * (ph.h0 * b)**2))
    return max(int(ph.emidiv * ph.emidiq * nrad),
               int(abs(ph.h0 * b / ph.brhoz * ph.anrad) / EPSRAD / 1e6 * ph.emidiv * ph.emidib))


def _tsolque(st, aln, akn, bz, eps0, enarad):
    """tsolque.f (bz≠0, ak0x=ak0y=0): 内部分割 ndiv、スライスごとに [tqente + 本体] → tradke (最後以外)、最後に tqente(半分) → tradke"""
    bzh = 0.5 * bz
    eps = 0.1 if eps0 == 0.0 else 0.1 * eps0
    ndiv = min(1000, 1 + int(abs(aln * math.hypot(akn, bz)) / eps))
    aln1 = aln / ndiv
    akk = akn / aln
    for n in range(1, ndiv + 1):
        first = n == 1

        def body(S, i, first=first):
            _tomech(S, i, bz)
            _x.tz_step1(i, S, aln1, akk, bz, first)       # 半 (または全) ドリフト + 本体 (K64 の tqente + 本体)
            _tocan(S, i, bz)
        _piece(st, body)
        if enarad and n != ndiv:
            st.tradke(aln1, 0.0, bzh)

    def last(S, i):
        _tomech(S, i, bz)
        _x.tsolqu_bz_final1(i, S, aln1, bz / (1.0 + S[i, 5]))        # 最後に正準運動量へ戻す所まで入っている
    _piece(st, last)
    if enarad:
        st.tradke(aln1, 0.0, bzh)


def _tquade(st, pr):
    """tquade (bz≠0 か bz=0 の K1≠0。F1=FRINGE=0、ずれなし)"""
    L, ak, ct, stt, fringe, bz, radok = pr
    bzh = 0.5 * bz
    krad = radok and L != 0.0

    def rin(S, i):
        _n.rot1(i, S, ct, stt)
    _piece(st, rin)
    if krad:
        st.setr0(bzh)
    if fringe:
        _piece(st, lambda S, i: _x.ttfrin4_bz1(i, S, ak, L, bz))
    st.setr0(bzh)                              # f1in = 0: tquade の else 側
    if krad:
        nq = min(1000, max(1, _itgetqraddiv(st.ph, st.cod, ak, L, bzh)))
    else:
        nq = 1
    aln = L / nq
    akn = ak / nq
    for _ in range(nq):
        _tsolque(st, aln, akn, bz, 0.0, krad)
    if fringe:
        _piece(st, lambda S, i: _x.ttfrin4_bz1(i, S, -ak, L, bz))

    def rout(S, i):
        _n.rot1(i, S, ct, -stt)
    _piece(st, rout)


def emit(st, el):
    """emit6d.tturne などの要素ループから: el._solmap のある要素 (SOL 領域内) の写像と放射"""
    kind, pr, _mom = el._solmap
    if kind == _n.K_NONE:
        return
    if kind == _x.K_SOL_FRIN:
        dbz = pr[0]
        if dbz != 0.0:
            _piece(st, lambda S, i: _x.tsfrin1(i, S, dbz))
        return
    if kind == _x.K_SOL_DRIFT:
        al, bz, radok = pr
        bzh = 0.5 * bz
        st.setr0(bzh)
        if al == 0.0:
            return
        if bz == 0.0:
            _e._tdrife0(st, al)
            return
        rs = np.array([1], dtype=np.int64)
        ph = st.ph
        _piece(st, lambda S, i: _x.tdrift_solenoid1(i, S, al, bz, False, False, rs, ph.p0, ph.anrad, ph.cuc))
        if radok:
            st.tradke(al, 0.0, bzh)
        return
    if kind == _x.K_SOL_QUAD:
        if pr[5] == 0.0:
            _e._tquade(st, el, krad=True)        # bz=0 の領域内は普通の四極
        else:
            _tquade(st, pr)
        return
    raise NotImplementedError(f'solemit: kind {kind} は未対応')
