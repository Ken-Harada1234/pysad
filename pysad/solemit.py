"""Maps for SOL regions used by `Emittance[]`, CALC6D and EMIT, corresponding to tsole1 (tsole.f) in K64.
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
from . import compat as _compat
from . import spin as _spin

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


def _trk(st, al, bzh, **fld):
    """tradke(al, 0, bzh)。v2 (POL) は spin.tradke_spin に、その刻みの磁場の成分を渡す (E): bz (p0 で規格化。指定がなければ 2·bzh)、
    bx・by (刻みの中点の横磁場)、gb (その (x, y) 勾配 ((∂bx/∂x, ∂bx/∂y), (∂by/∂x, ∂by/∂y)))。端の縦磁場の積分 ψ は st.add_bsi で入れる。
    spin.tradke_spin は運動量の向きの磁場の積分 B∥ を、向きの変化から作る B⊥ と分けて足す (K64 の bsi は z 軸の向きに足す)"""
    if st.calpol and _compat.v2() and _spin.SOL_BPAR:
        h = dict(fld)
        h.setdefault('bz', 2.0 * bzh)
        st.bhint = h
    st.tradke(al, 0.0, bzh)


def _mult_fld(prm):
    """領域内の MULT の横磁場 (p0 で規格化) F(x, y) → (bx, by, ψ) (実験室系の (x, y))。By + i Bx = Σ k_n z^n/n! (k_n = 積分強度/L、z = 局所座標 x + iy:
    ずらして theta2 だけ回した座標。track6._RadHook._field と同じ式)、ψ = Im Σ k_n z^(n+1)/(n+1)! (端の縦磁場の積分。∂ψ/∂x = Bx、∂ψ/∂y = By)。
    磁場のベクトルは実験室系に回し戻す"""
    al, nmmax, akr, _dofr, _nd, dx, dy, ct, stt = prm[:9]
    k = [a / al for a in akr[:nmmax + 1]]

    def F(x, y):
        xl, yl = x - dx, y - dy
        z = complex(ct * xl - stt * yl, stt * xl + ct * yl)
        f = 0j
        g = 0j
        zn = 1.0 + 0j
        fact = 1.0
        for n, kn in enumerate(k):
            f += kn * zn / fact
            zn *= z
            fact *= (n + 1)
            g += kn * zn / fact
        bxl, byl = f.imag, f.real
        return ct * bxl + stt * byl, -stt * bxl + ct * byl, g.imag
    return F


def _mult_trk(st, F, al, bzh):
    """領域内の MULT の 1 刻みの tradke: 刻みの中点の横磁場 (bx, by) とその勾配を渡す (v2 の POL)"""
    if not (st.calpol and _compat.v2() and _spin.SOL_BPAR):
        st.tradke(al, 0.0, bzh)
        return
    xm = 0.5 * (st.codr0[0] + st.cod[0])
    ym = 0.5 * (st.codr0[2] + st.cod[2])
    h = 1e-6
    bx, by, _ = F(xm, ym)
    bxp, byp, _ = F(xm + h, ym)
    bxm, bym, _ = F(xm - h, ym)
    bxq, byq, _ = F(xm, ym + h)
    bxr, byr, _ = F(xm, ym - h)
    gb = (((bxp - bxm) / (2 * h), (bxq - bxr) / (2 * h)), ((byp - bym) / (2 * h), (byq - byr) / (2 * h)))
    _trk(st, al, bzh, bx=bx, by=by, gb=gb)


def _mult_psi(st, F, sgn, at_r0):
    """端の縦磁場の積分 sgn·ψ(x, y) (入口 +、出口 −) を bsir0 に足す (勾配も)"""
    if not (st.calpol and _compat.v2() and _spin.SOL_BPAR):
        return
    x, y = st.cod[0], st.cod[2]
    bx, by, psi = F(x, y)
    st.add_bsi(sgn * psi, sgn * bx, sgn * by, at_r0)


def _spin_rot(st, th, ent):
    """tsolrote の calpol (chi1=chi2=0): スピンを z 軸まわりに th 回す (入口)、戻す (出口)。座標の _n.rot1(ct, ±stt) と同じ向き"""
    if not st.calpol or th == 0.0:
        return
    c, s = math.cos(th), math.sin(th)
    r = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    st.srot = (r if ent else r.T) @ st.srot


def _itgetqraddiv(ph, c, ak, al, bzh):
    """tsolque.f itgetqraddiv (bz あり)。c は正準運動量のままの cod"""
    xk = abs(ak) * max(1e-6, abs(c[0]) + abs(c[2]))
    xpd = max(1e-6, abs(c[1]) + abs(c[3]))
    a = min(1e-2, xk + xpd)
    b = ph.brhoz * (a * abs(bzh) + xk / abs(al))
    nrad = int(abs(al * ph.crad / EPSRAD * (ph.h0 * b)**2))
    return max(int(ph.emidiv * ph.emidiq * nrad),
               int(abs(ph.h0 * b / ph.brhoz * ph.anrad) / EPSRAD / 1e6 * ph.emidiv * ph.emidib))


def _quad_trk(st, kq, al, bzh):
    """領域内の四極の 1 刻みの tradke。v2 (B232) の POL は、刻みの中点の横磁場 (bx, by) = (kq y, kq x) とその勾配を渡して、本体の運動量の向きの磁場
    ∫B·dr = (bx px + by py) al/pz を足す (回した座標系の四極。K64 は螺旋の縦磁場 bz·al だけで本体の横磁場の B∥ が無い)"""
    if kq is None or not (st.calpol and _compat.v2() and _spin.SOL_BPAR and _compat.spin_bpar()):
        _trk(st, al, bzh)
        return
    xm = 0.5 * (st.codr0[0] + st.cod[0])
    ym = 0.5 * (st.codr0[2] + st.cod[2])
    _trk(st, al, bzh, bx=kq * ym, by=kq * xm, gb=((0.0, kq), (kq, 0.0)))


def _tsolque(st, aln, akn, bz, eps0, enarad, psi_out=None, kq=None):
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
            _x.tz_step1(i, S, aln1, akk, bz, first, _compat.v2())       # 半 (または全) ドリフト + 本体 (K64 の tqente + 本体)
            _tocan(S, i, bz)
        _piece(st, body)
        if enarad and n != ndiv:
            _quad_trk(st, kq, aln1, bzh)

    def last(S, i):
        _tomech(S, i, bz)
        _x.tsolqu_bz_final1(i, S, aln1, bz / (1.0 + S[i, 5]))        # 最後に正準運動量へ戻す所まで入っている
    _piece(st, last)
    if enarad:
        if psi_out is not None:
            psi_out()                      # v2 (E): 四極の出口の端の縦磁場の積分 −k x y
        _quad_trk(st, kq, aln1, bzh)


def _tquade(st, pr):
    """tquade (bz≠0 か bz=0 の K1≠0。F1=FRINGE=0、ずれなし)"""
    L, ak, ct, stt, fringe, bz, radok = pr[:7]
    qdx, qdy = (pr[7], pr[8]) if len(pr) > 8 else (0.0, 0.0)      # v2next B3: 四極のずれ (tsolrot の dx, dy。ゲージの分だけ正準運動量もずらす)
    bzh = 0.5 * bz
    krad = radok and L != 0.0
    th2 = math.atan2(stt, ct)

    def rin(S, i):
        if qdx != 0.0 or qdy != 0.0:
            bzp = bz / (1.0 + S[i, 5])
            S[i, 0] -= qdx
            S[i, 2] -= qdy
            S[i, 1] += 0.5 * bzp * qdy
            S[i, 3] -= 0.5 * bzp * qdx
        _n.rot1(i, S, ct, stt)
    _piece(st, rin)
    _spin_rot(st, th2, True)                   # tsolrote (calpol): スピンも回す (E。K64 は srot = rr·srot)
    if krad:
        st.setr0(bzh)
    if fringe:
        _piece(st, lambda S, i: _x.ttfrin4_bz1(i, S, ak, L, bz))
    st.setr0(bzh)                              # f1in = 0: tquade の else 側
    if krad:
        nq = min(1000, max(1, _itgetqraddiv(st.ph, st.cod, ak, L, bzh)))
    else:
        nq = 1
    if krad and st.calpol and _compat.v2():
        # v2 (POL): らせんの回転角が大きいと、1 刻みの放射・スピンの歳差の微分は弦から出るので不正確。bz·L/nq ≤ 0.02 に刻む (領域内の DRIFT・BEND と同じ)
        nq = min(1000, max(nq, 1 + int(abs(bz * L) / 0.02)))
    aln = L / nq
    akn = ak / nq
    kq = ak / L                                # 四極の勾配 (|K1|、回した座標系で By = kq x、Bx = kq y)
    if krad and st.calpol and _compat.v2() and _spin.SOL_BPAR and fringe:
        # v2 (E): 端の縦磁場の積分 ψ = kq x y を入口で +、出口で − (追跡 track6._tsolqu_bz_ens と同じ。P4: DISFRIN=1 では入れない)
        st.add_bsi(kq * st.cod[0] * st.cod[2], kq * st.cod[2], kq * st.cod[0], True)

        def _psi_out():
            st.add_bsi(-kq * st.cod[0] * st.cod[2], -kq * st.cod[2], -kq * st.cod[0], False)
    else:
        _psi_out = None
    for iq in range(nq):
        _tsolque(st, aln, akn, bz, 0.0, krad, psi_out=_psi_out if iq == nq - 1 else None, kq=kq)
    if fringe:
        _piece(st, lambda S, i: _x.ttfrin4_bz1(i, S, -ak, L, bz))

    def rout(S, i):
        _n.rot1(i, S, ct, -stt)
        if qdx != 0.0 or qdy != 0.0:
            bzp = bz / (1.0 + S[i, 5])
            S[i, 0] += qdx
            S[i, 2] += qdy
            S[i, 1] -= 0.5 * bzp * qdy
            S[i, 3] += 0.5 * bzp * qdx
    _piece(st, rout)
    _spin_rot(st, th2, False)


def _cav_fringe(st, el, al, sgn):
    """C2: 領域内の CAVI・RF つき MULT の端の RF のフリンジ (tcavfrie)。追跡 (sol_mult_full の cav_frin_all) と同じ位置に入れる。
    sgn = +1 入口、−1 出口。RFSW OFF なら何もしない (空洞は働かない)"""
    if not st.rfsw or al == 0.0:
        return
    from . import lattice as _lat
    ph = st.ph
    gv = lambda k: float(getattr(el, k, 0.0) or 0.0)
    vc = gv('VOLT') + gv('DVOLT')
    if vc == 0.0:
        return
    harm = gv('HARM') or gv('_harm')
    w = ph.omega0 * harm / _e.CVELOC if harm != 0.0 else 2.0 * math.pi * gv('FREQ') / _e.CVELOC
    v = vc / ph.amass * abs(ph.charge)
    phis = st.trf0 * w
    phic = gv('PHI') * (1.0 if ph.charge >= 0 else -1.0)
    s0 = 0.0 if st.radcod else math.sin(phis)
    cdx, cdy, cth = gv('DX'), gv('DY'), (gv('ROTATION') if type(el).__name__ == 'Cavity' else 0.0)
    sh = type(el).__name__ == 'Cavity' and (cdx != 0.0 or cdy != 0.0 or cth != 0.0)
    if sh:
        st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True),
                 jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True))
    st.piece(_e._tcavfrie_map(al, sgn * v, w, phis - phic, s0, ph.p0), symp_z=False)
    if sh:
        st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False),
                 jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False))


def emit(st, el):
    """Apply the map and radiation of an element with el._solmap (inside a SOL region) from the element loops of emit6d.tturne.
    emit6d.tturne などの要素ループから: el._solmap のある要素 (SOL 領域内) の写像と放射"""
    kind, pr, _mom = el._solmap
    if kind == _n.K_NONE:
        return
    if kind == 41:                       # phsrot.K_PHSROT (v2)
        from . import phsrot as _prt
        _prt.emit(st, pr[0])
        return
    if kind == _n.K_UND:
        from . import und as _un
        _un.emit(st, pr)                 # v2: 磁場のある UND
        return
    if kind == _x.K_SOL_FRIN:
        dbz = pr[0]
        if dbz != 0.0:
            _piece(st, lambda S, i: _x.tsfrin1(i, S, dbz))
        return
    if kind == _x.K_SOL_EDGE:
        # v2 (C2): SOL の端のずれ・回転 (tsconv = trots) と F1 の放射 (tsole1 の icSOL)。
        # 入口: 回転 → tsfrie → tradke(F1, bz/2)。出口: tsfrie → tradke(F1, 0) → 回転。領域内の SOL: tsfrie → tradke(F1, bz1/2)
        rot, r, dx, dy, dz, dbz, bz0, bz1, f1, radok = pr
        if rot == 1:
            _piece(st, lambda S, i: _x.trots1(i, S, r, dx, dy, dz, True))
            if st.calpol:
                st.srot = r @ st.srot              # tsole1: srot = matmul(rr, srot) (calpol)
        if radok:
            st.setr0(0.5 * bz0)
        if dbz != 0.0:
            _piece(st, lambda S, i: _x.tsfrin1(i, S, dbz))
        if radok:
            _trk(st, f1, 0.5 * bz1, bz=0.0)         # v2 (E): 端の放射は縦磁場の積分なし (追跡の _sol_edge は bpar = 0)
        if rot == -1:
            _piece(st, lambda S, i: _x.trots1(i, S, r, dx, dy, dz, False))
            if st.calpol:
                st.srot = r @ st.srot
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
        # v2 (B1・E): らせんの回転角 bz·al が大きいと、1 刻みの放射・スピンの歳差は弦から出るので不正確 (追跡 k_sol_drift・track6._sol_drift_v2 と同じ刻み)
        nd = (1 + int(abs(bz * al) / 0.02)) if (radok and _compat.v2()) else 1
        for _ in range(nd):
            _piece(st, lambda S, i: _x.tdrift_solenoid1(i, S, al / nd, bz, False, False, rs, ph.p0, ph.anrad, ph.cuc))
            if radok:
                _trk(st, al / nd, bzh)
        return
    if kind == _x.K_SOL_QUAD:
        if pr[5] == 0.0:
            _e._tquade(st, el, krad=True)        # bz=0 の領域内は普通の四極
        else:
            _tquade(st, pr)
        return
    if kind == _n.K_APERT:
        return                                   # 絞り: 写像は恒等
    if kind == _x.K_SOL_BEND:
        # v2: 領域内の BEND (一様な双極 + ソレノイド)。刻みごとに写像 → tradke (bz つき)
        al, bz, kx, ky, radok = pr
        bzh = 0.5 * bz
        st.setr0(bzh)
        nd = max(1 + int(abs(bz * al) / 0.02), 1 + int(math.hypot(kx, ky) / 0.005))
        rs = np.array([1], dtype=np.int64)
        for _ in range(nd):
            _piece(st, lambda S, i: _x.tdrift_bz1(i, S, al / nd, bz, kx / nd, ky / nd, False, False, rs, 0.0, 0.0, 0.0))
            if radok:
                _trk(st, al / nd, bzh, bx=ky / al, by=kx / al)       # 一様な双極磁場 (By = kx/al、Bx = ky/al)
        return
    if kind == _x.K_SOL_MULT:
        # v2: 領域内の MULT・CAVI。nd 個に分け (端のフリンジは最初と最後だけ)、分けた要素ごとに写像 → tradke
        prm, bz = pr
        (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, stt, fringe, eps0, radok, rf) = prm
        bzh = 0.5 * bz
        ph = st.ph
        rs = np.array([1], dtype=np.int64)
        st.setr0(bzh)
        if rf is not None and type(el).__name__ != 'Cavity':
            # C2: RF つき MULT (tmulte の RF 部)。nd 個に分け、[多極つきの本体の半分 → 薄い空洞 (emit6d._tcave) → 本体の半分]。
            # 本体は追跡と同じカーネル (ソレノイド + 多極)、空洞の蹴りは CAVI と同じ式 (電圧は MULT の VOLT・FREQ/HARM・PHI)
            from types import SimpleNamespace as _NS
            k1 = abs(akr[1]) if nmmax >= 1 else 0.0
            nd = max(1 + int(abs(bz * al) / 0.02), 1 + int((abs(akr[0]) + k1 * 1e-2) / 0.005), 1 + int(abs(rf[1] * al)))
            gv = lambda k: float(getattr(el, k, 0.0) or 0.0)
            thin = _NS(L=0.0, VOLT=gv('VOLT') / nd, DVOLT=gv('DVOLT') / nd, HARM=gv('HARM'), _harm=gv('_harm'), FREQ=gv('FREQ'),
                       PHI=gv('PHI'), DISFRIN=1.0, FRINGE=0.0, DX=0.0, DY=0.0, ROTATION=0.0, AUTOPHI=gv('AUTOPHI'), name=el.name)
            hl = al / nd * 0.5
            F = _mult_fld(prm)
            if fringe:
                _cav_fringe(st, el, al, 1.0)
            if radok and al != 0.0 and fringe:
                _mult_psi(st, F, 1.0, True)
            for k in range(nd):
                for half in (0, 1):
                    sub = (hl, nmmax, akr / (2 * nd), dofr, ndiv0, dx, dy, ct, stt, fringe, eps0, False, None)
                    fi, fo = (k == 0 and half == 0), (k == nd - 1 and half == 1)

                    def f(S, i, sub=sub, fi=fi, fo=fo):
                        T = S[i:i + 1]
                        _x.sol_mult_full(T, 1, ph, sub, bz, 0.0, False, False, False, False, rs, fr_in=fi, fr_out=fo)
                    _piece(st, f)
                    if radok and al != 0.0:
                        if fo and fringe:
                            _mult_psi(st, F, -1.0, False)
                        _mult_trk(st, F, hl, bzh)
                    if half == 0:
                        _e._tcave(st, thin)
                        st.setr0(bzh)
            if fringe:
                _cav_fringe(st, el, al, -1.0)
            return
        if rf is not None:
            # CAVI: [ソレノイドのドリフト半分 → 薄い空洞 (emit6d._tcave: RF の積算量 vc0 なども更新) → 半分] × nd
            import copy as _copy
            nd = 1 + int(abs(bz * al) / (0.02 if (st.calpol and _compat.v2()) else 0.1)) + int(abs(rf[1] * al))      # v2 (POL): スピンの歳差の刻み
            thin = _copy.copy(el)
            thin.L = 0.0
            thin.VOLT = float(getattr(el, 'VOLT', 0.0) or 0.0) / nd
            if getattr(el, 'DVOLT', None):
                thin.DVOLT = float(el.DVOLT) / nd
            h = al / nd * 0.5
            if fringe:
                _cav_fringe(st, el, al, 1.0)           # C2: 端の RF のフリンジ (追跡と同じ。以前は入っていなかった)
            for _ in range(nd):
                if h != 0.0:
                    _piece(st, lambda S, i: _x.tdrift_solenoid1(i, S, h, bz, False, False, rs, ph.p0, ph.anrad, ph.cuc))
                    _trk(st, h, bzh)
                _e._tcave(st, thin)
                st.setr0(bzh)
                if h != 0.0:
                    _piece(st, lambda S, i: _x.tdrift_solenoid1(i, S, h, bz, False, False, rs, ph.p0, ph.anrad, ph.cuc))
                    _trk(st, h, bzh)
            if fringe:
                _cav_fringe(st, el, al, -1.0)
            return
        if al == 0.0:
            nd = 1
        else:
            k1 = abs(akr[1]) if nmmax >= 1 else 0.0
            nd = max(1 + int(abs(bz * al) / 0.02), 1 + int((abs(akr[0]) + k1 * 1e-2) / 0.005))
        F = _mult_fld(prm)
        if radok and al != 0.0 and fringe:
            _mult_psi(st, F, 1.0, True)                # 入口の端の縦磁場の積分 +ψ (st.setr0 の後。P4: DISFRIN=1 では入れない)
        for k in range(nd):
            sub = (al / nd, nmmax, akr / nd, dofr, ndiv0, dx, dy, ct, stt, fringe, eps0, False, None)
            fi, fo = k == 0, k == nd - 1

            def f(S, i, sub=sub, fi=fi, fo=fo):
                T = S[i:i + 1]
                _x.sol_mult_full(T, 1, ph, sub, bz, 0.0, False, False, False, False, rs, fr_in=fi, fr_out=fo)
            _piece(st, f)
            if radok and al != 0.0:
                if k == nd - 1 and fringe:
                    _mult_psi(st, F, -1.0, False)
                _mult_trk(st, F, al / nd, bzh)
        return
    # C2: bz = 0 の領域の中の CAVI・MULT・DRIFT は、領域の外と同じ写像 (compile_sol_region が普通の kind にする)
    if kind == _n.K_CAV:
        _e._tcave(st, el)
        return
    if kind == _x.K_MULT:
        if pr[-1] is not None:
            raise NotImplementedError('SOL 領域 (bz = 0) の中の RF つき MULT の Emittance[]・EMIT は未対応 (bz≠0 の領域と、追跡は対応)')
        _e._tmulte(st, el)
        return
    if kind == _n.K_DRIFT:
        if pr[0] != 0.0:
            _e._tdrife0(st, pr[0])
        return
    raise NotImplementedError(f'solemit: kind {kind} は未対応')
