"""Translation of the SAD spin (POL flag) code: tradke spin part, spnorm, srequpol, spdepol and srot of tchge.f.
スピン (POL フラグ) の直訳: photons.f の tradke のスピン部、tspin の spnorm / srequpol / spdepol、tchge.f の srot 回転。

状態 (temw / temit の srot(3,9), gintd(3)):
    srot[:, 0:3]  スピン歳差の 1 周行列 (srotinit で単位行列)
    srot[:, 3:9]  (x,px,y,py,z,pz) の初期座標に対する 1 周のスピン微分 (脱偏極ベクトル)
    gintd[3]      sflc·anp·(bt·h1·p/al1)² のスピン軸成分の積算 (Sokolov-Ternov の傾き)。1/sdamp が偏極時間
K64 の癖 (docs/30 B55〜B57) と v2 (compat.v2()) の扱い (docs/history/v2_changes_A.md):
  - tadddsg の呼び出しは dgz に dgx を渡す (dsg は spequ でしか使わず、spequ は呼ばれないので無影響。ここでは dsg を持たない)
  - tbrote の srot 更新は列ごとに行を逐次更新する (rr·srot の行列積ではない: 2 行目が更新後の 1 行目を使う)。v2 は行列積 (B56)
  - serot の c2 = (−s·c1' + c2)/c と、tradke の末尾の srot(3,:) = (−sp·srot(1,:)' + srot(3,:))/cp は、更新後の値 (') を使うが、
    展開すると −s·c1 + c·c2 で、正しい回転をその場で計算する書き方 (誤りではない。v2 でも同じ)
  - 平衡偏極度: K64 の spdepol は NaN (cerfc の桁あふれ)。v2 は桁あふれを直したうえで、線形のスピン軌道結合の脱偏極率
    (slim_depolarization) を使う (B55)
"""
from __future__ import annotations

import math

import numpy as np

GSPIN = 0.001159652181280002
CL = 1.0 + GSPIN
PST = 8.0 * math.sqrt(3.0) / 15.0
CFPD = 5.0 * math.sqrt(3.0) / 8.0
CEZLIM = 18.0
SOL_BPAR = True            # v2 (E): SOL 領域の中の Emittance[] の POL は、運動量の向きの磁場の積分 B∥ を向きの変化から作る B⊥ と分けて足す。False で K64 のまま (bsi を z 軸の向きに足す)
SPIN_ANGLE_EXACT = True    # v2 (E): tradke_spin の向きの変化の大きさを sin θ でなく θ にする (追跡と同じ。False で A の v2 のまま)
DBTZ_PR2 = True          # v2 (B336): dbtz を pr² で割る (試験が False にして K64 の式 (pr) と比べる)
STEP_RATE = True           # v2 (P1): 線形の脱偏極率を、放射の刻みごとのスピン軸の場 N_s (要素の出口の N ではなく) で積分する。False で B3 までの要素の出口の N
KINETIC = True             # v2 (P1): 平衡偏極度の分子に運動論的偏極 b·∂n/∂δ (DK の式の −b·γ∂n/∂γ の項) を入れる。False で n·b だけ (B3 まで)
STEP_ACC = False           # slim_depolarization の中だけ True: tradke_spin が刻みごとの積算量 (kin_A・kin_Y・rt_M・rt_F・rt_G) を作る
ST_EXACT = True            # v2 (P1): 緩和率に Sokolov-Ternov の厳密式 Σ r (1 − 2/9 (n·v)²) を使う。False で A2 の Σ r (b·n)² (傾いた n・垂直の偏向で合わない)
CHROMATIC = True           # v2 (P1): 側帯の式にベータトロンの位相のクロマティックな変調 (チューンの δ 依存 ξ) を入れる。False で B3 の式 (軌道は線形)
CHROM_KMAX = 60000         # 減衰の和の項数 K がこれを超えるとき (放射減衰が弱い大きなリング) は入れない (計算が K×試料数×6 で重い)
CHROM_NEL = 1200           # 要素数がこれを超えるとき (大きなリング A: 3563) は dM/dδ の 2 周も計算しない (未確認の範囲)
CHROM_MIN = 0.03           # 変調指数 (dν/dδ)·σδ/νs がこの値より小さいときは入れない (効果が 1% 未満)
SIDEBAND = True       # v2 (B3): Polarization の脱偏極率にシンクロトロン側帯 (スピンの位相変調) を入れる。False で線形の式 (A の B55) だけ
M_2PI = 2.0 * math.pi
M_SQRTPI = math.sqrt(math.pi)


def sflc():
    from .emit6d import ELRADI, FINEST
    return 0.75 * (ELRADI / FINEST) ** 2


def srotinit():
    s = np.zeros((3, 9))
    s[0, 0] = s[1, 1] = s[2, 2] = 1.0
    return s


def xsincos(x):
    """Compute sin, cos and the accurate form of 1 - cos (s, c, dc), as in xsincos of tfloor.f.
    tfloor.f xsincos: sin, cos と 1−cos の精度の良い形 (s, c, dc)"""
    s = math.sin(x)
    c = math.cos(x)
    dc = -2.0 * math.sin(0.5 * x) ** 2
    return s, c, dc


def _pxy2dpz(px, py):
    a = px * px + py * py
    return -a / (1.0 + math.sqrt(max(0.0, 1.0 - a)))


# ---------------------------------------------------------------------------
# tradke のスピン部 (photons.f 1521-1616)
# ---------------------------------------------------------------------------

def _edge_rot(st, psi, dpsi, pv, dpv, pr, h1, dh1r, transr):
    """v2 (B137): 端の縦磁場の積分 ψ (p0 で規格化。粒子の運動量 (1+δ) で規格化し直す) だけのスピンの歳差を srot (3x9) に掛ける。
    運動量の向き pv (3 成分、単位ベクトル) の成分に (1+a)、垂直な成分に (1+aγ)。dpsi・dpv はローカルな 6 成分座標 (放射の基準点) についての微分 (6 / 3x6)。
    座標系は回さない。srot[:, 3:9] に回転ベクトルの微分を足し (transr 経由)、9 列全部を回す (tradke_spin の本体の歳差と同じ式)"""
    beta = psi / pr
    dbeta = dpsi / pr
    dbeta[5] -= psi / pr ** 2
    b = np.array([0.0, 0.0, beta])
    db = np.zeros((3, 6))
    db[2, :] = dbeta
    bpk = float(b @ pv)
    dbpk = db.T @ pv + dpv.T @ b
    bl = bpk * pv
    dbl = np.outer(pv, dbpk) + bpk * dpv
    bt = b - bl
    dbt = db - dbl
    ct = h1 * GSPIN
    dct = ct * dh1r
    ct += 1.0
    g3 = ct * bt + CL * bl
    dg3 = ct * dbt + CL * dbl
    dg3[:, 5] += dct * bt
    srot = st.srot
    for i in range(3):
        srot[i, 3:9] += dg3[i] @ transr
    g = math.sqrt(float(g3 @ g3))
    if g != 0.0:
        sinu, cosu, dcosu = xsincos(g)
        sinu = sinu / g
        gx, gy, gz = g3
        sx, sy, sz = srot[0].copy(), srot[1].copy(), srot[2].copy()
        sw = -(sx * gx + sy * gy + sz * gz) * dcosu / g ** 2
        sux = sy * gz - sz * gy
        suy = sz * gx - sx * gz
        suz = sx * gy - sy * gx
        srot[0] = cosu * sx + sinu * sux + sw * gx
        srot[1] = cosu * sy + sinu * suy + sw * gy
        srot[2] = cosu * sz + sinu * suz + sw * gz


def tradke_spin(st, c):
    """Update the spin one-turn matrix and depolarization vector (st.srot, st.gintd) from the local quantities of tradke.
    c: _State.tradke のローカル量の dict。st.srot / st.gintd を更新する。transr は更新前の値 (st.transr)"""
    pr = c['pr']
    px, py, pz = c['px'], c['py'], c['pz']
    pxi, pyi, pz00, pz0, pxr0 = c['pxi'], c['pyi'], c['pz00'], c['pz0'], c['pxr0']
    cp, sp = c['cp'], c['sp']
    xpx, xpy, xpz, xpa = c['xpx'], c['xpy'], c['xpz'], c['xpa']
    bzh, bzhr0, al = c['bzh'], c['bzhr0'], c['al']
    h1, p, p0, al1, anp = c['h1'], c['p'], c['p0'], c['al1'], c['anp']
    from .emit6d import _tinv6
    transr = st.transr
    tr2 = st.T @ _tinv6(transr)
    if bzh != 0.0:
        tr2[1, :] = tr2[1, :] + bzh * tr2[2, :]
        tr2[3, :] = tr2[3, :] - bzh * tr2[0, :]
    ddpz = (tr2[5, :] * pr - tr2[1, :] * px - tr2[3, :] * py) / pz
    dpxi = np.array([0.0, 1.0, bzhr0, 0.0, 0.0, 0.0])
    dpyi = np.array([-bzhr0, 0.0, 0.0, 1.0, 0.0, 0.0])
    dpz00 = (-pxi * dpxi - pyi * dpyi) / pz00
    dpz00[5] += pr / pz00
    dpxr0 = cp * dpxi + sp * dpz00
    dpz0 = (-pxr0 * dpxr0 - pyi * dpyi) / pz0
    dpz0[5] += pr / pz0
    dxpx = tr2[3, :] * pz0 + py * dpz0 - ddpz * pyi - pz * dpyi
    dxpy = ddpz * pxr0 + pz * dpxr0 - tr2[1, :] * pz0 - px * dpz0
    dxpz = tr2[1, :] * pyi + px * dpyi - tr2[3, :] * pxr0 - py * dpxr0
    dh1r = p * p0 / h1 ** 2
    from . import compat as _compat0
    ein = getattr(st, 'edge_in', None)
    eout = getattr(st, 'edge_out', None)
    st.edge_in = st.edge_out = None
    if ein is not None and _compat0.spin_edge_split():
        # v2 (B137): 入口の端の ψ は、刻みの本体の前の独立な回転 (放射の基準点の運動量の向き)
        pv_in = np.array([pxi, pyi, pz00]) / pr
        dpv_in = np.array([dpxi, dpyi, dpz00]) / pr
        dpv_in[:, 5] -= pv_in / pr
        _edge_rot(st, ein[0], ein[1].copy(), pv_in, dpv_in, pr, h1, dh1r, transr)

    bsir0 = st.bsir0
    pxh = (pxr0 + px) / pr * 0.5
    pyh = (pyi + py) / pr * 0.5
    dpxh = (tr2[1, :] + dpxr0) / pr * 0.5
    dpxh[5] -= pxh / pr
    dpyh = (tr2[3, :] + dpyi) / pr * 0.5
    dpyh[5] -= pyh / pr
    pzh = 1.0 + _pxy2dpz(pxh, pyh)
    dpzh = -(pxh * dpxh + pyh * dpyh) / pzh
    hint = getattr(st, 'bhint', None)
    st.bhint = None
    from . import compat as _compat
    if hint is not None and _compat.v2():
        # v2 (E): SOL 領域の中。K64 は bsir0 + bz·al = ∫Bs ds (p0 で規格化) を z 軸の向きの磁場として足す。運動量が z 軸から傾いていると
        # (傾いた・ずれたソレノイドの軌道、軸を外れた粒子)、(1) 運動量に垂直な成分が B⊥ に二重に入り (向きの変化 xpx.. はすでに B⊥)、
        # (2) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない。追跡 (track6._sprot の bpar) と同じく、向きの変化から作る B⊥ と、
        # 運動量の向きの磁場の積分 B∥ = [bz·al + (bx px + by py)·al/pz + ψ]/(1+δ) を分けて足す。(bx, by) は刻みの中点での横磁場の成分 (p0 で規格化。
        # hint['gb'] はその (x, y) についての勾配)、ψ は端の縦磁場の積分 (bsir0 と bsi_d0/d1)。向きの変化の大きさは角そのもの (asin。追跡の _sprot と同じ)
        b = np.array([xpx, xpy, xpz]) / pr ** 2
        db = np.array([dxpx, dxpy, dxpz]) / pr ** 2
        db[:, 5] -= 2.0 * b / pr
        ppa = float(np.linalg.norm(b))
        if ppa > 1e-12:
            fs = math.asin(min(1.0, ppa)) / ppa
            b = b * fs
            db = db * fs
        ph3 = np.array([pxh, pyh, pzh])
        dph3 = np.array([dpxh, dpyh, dpzh])
        bpk = float(b @ ph3)
        dbpk = db.T @ ph3 + dph3.T @ b
        bt3 = b - bpk * ph3
        dbt3 = db - np.outer(ph3, dbpk) - bpk * dph3
        bzf, bxf, byf = hint.get('bz', 0.0), hint.get('bx', 0.0), hint.get('by', 0.0)
        pxy = bxf * pxh + byf * pyh
        num = bzf * al + pxy * al / pzh + bsir0
        dnum = (bxf * dpxh + byf * dpyh) * al / pzh - pxy * al / pzh ** 2 * dpzh
        gb = hint.get('gb')
        if gb is not None:
            # 場の (bx, by) は刻みの中点の位置 (x_m, y_m) の関数。∂x_m/∂r0 = (e_x + tr2[0])/2、∂y_m/∂r0 = (e_y + tr2[2])/2
            dxm = 0.5 * tr2[0, :].copy()
            dxm[0] += 0.5
            dym = 0.5 * tr2[2, :].copy()
            dym[2] += 0.5
            dbx = gb[0][0] * dxm + gb[0][1] * dym
            dby = gb[1][0] * dxm + gb[1][1] * dym
            dnum = dnum + (dbx * pxh + dby * pyh) * al / pzh
        bd0 = getattr(st, 'bsi_d0', None)
        if bd0 is not None and (bd0.any() or st.bsi_d1.any()):
            dnum = dnum + bd0 + st.bsi_d1 @ tr2
        bpar = num / pr
        dbpar = dnum / pr
        dbpar[5] -= num / pr ** 2
        bl3 = bpar * ph3
        dbl3 = np.outer(ph3, dbpar) + bpar * dph3
        blx, bly, blz = bl3
        btx, bty, btz = bt3
        dblx, dbly, dblz = dbl3
        dbtx, dbty, dbtz = dbt3
    else:
        if _compat.v2() and SPIN_ANGLE_EXACT:
            # v2 (E、A の N8): 向きの変化 xpx.. = sin(θ) pr² の大きさを、追跡 (track6._sprot の v2) と同じ角 θ = asin(..) にする。K64 は sin θ のまま
            # (刻みの角 θ に対し相対 θ²/6 小さい。0.5 GeV の平面リングでスピンチューンが 2.6e-5 違った)。θ の微分の項 (θ²/6 の微分) は入れない
            _ppa = math.sqrt(xpx * xpx + xpy * xpy + xpz * xpz) / pr ** 2
            if _ppa > 1e-12:
                _fs = math.asin(min(1.0, _ppa)) / _ppa
                xpx, xpy, xpz = xpx * _fs, xpy * _fs, xpz * _fs
                dxpx, dxpy, dxpz = dxpx * _fs, dxpy * _fs, dxpz * _fs
        if _compat.v2() and _compat.spin_bpar():
            # v2 (B242): 端の縦磁場の積分 bsir0 (p0 で規格化した ψ) を粒子の運動量 (1+δ) で規格化し直す。K64 は割らないので、本体の Δψ/(1+δ) (B165) との和が
            # δΔψ だけ 0 にならない。螺旋の縦磁場 bzh·2·al は従来のまま
            xpzb = xpz + bsir0 * pr + bzh * 2.0 * al * pr ** 2
            dxpzb = dxpz.copy()
            dxpzb[5] += bsir0 + 2.0 * bzh * 2.0 * al * pr
            bnorm = pr
        else:
            xpzb = xpz + (bsir0 + bzh * 2.0 * al) * pr ** 2
            dxpzb = dxpz.copy()
            dxpzb[5] += 2.0 * (bsir0 + bzh * 2.0 * al) * pr
            bnorm = pr ** 2
        bd0 = getattr(st, 'bsi_d0', None)
        if bd0 is not None and (bd0.any() or st.bsi_d1.any()):
            # v2 (B58): bsir0 (端の縦磁場の積分 ψ(x, y)) の軌道についての微分。K64 は bsir0 を定数として扱い、
            # スピン軌道結合 (脱偏極ベクトル) に端の縦磁場の y・x 依存が入らない (追跡のスピンの軌道依存と BEND で 6 倍違う)
            dxpzb = dxpzb + (bd0 + st.bsi_d1 @ tr2) * bnorm
        bp = (xpx * pxh + xpy * pyh + xpzb * pzh) / pr ** 2
        dbp = (dxpx * pxh + xpx * dpxh + dxpy * pyh + xpy * dpyh + dxpzb * pzh + xpzb * dpzh) / pr ** 2
        dbp[5] -= 2.0 * bp / pr
        blx, bly, blz = bp * pxh, bp * pyh, bp * pzh
        btx = xpx / pr ** 2 - blx
        bty = xpy / pr ** 2 - bly
        btz = xpzb / pr ** 2 - blz
        dblx = dbp * pxh + bp * dpxh
        dbly = dbp * pyh + bp * dpyh
        dblz = dbp * pzh + bp * dpzh
        dbtx = dxpx / pr ** 2 - dblx
        dbtx[5] -= 2.0 * xpx / pr ** 3
        dbty = dxpy / pr ** 2 - dbly
        dbty[5] -= 2.0 * xpy / pr ** 3
        # v2 (B336): btz = xpzb/pr² − blz の微分は dxpzb/pr² (x・y 成分と同じ)。K64 は pr で割る (誤植。δ ≠ 0 の軌道で効く)
        dbtz = dxpzb / (pr ** 2 if (_compat.v2() and DBTZ_PR2) else pr) - dblz
        dbtz[5] -= 2.0 * xpzb / pr ** 3
    bpot = c.get('bpot')
    if bpot is not None and not (hint is not None and _compat.v2()):
        # v2 (B165): 本体の運動量の向きの磁場の積分 ∫B∥ ds = ψ(終) − ψ(始) (静磁場の ∫B·dr、ψ = ポテンシャル)。粒子の運動量 (1+δ) で規格化して運動量の向きに足す
        from .emit6d import _grad_xy
        cod1 = st.cod
        cod0 = c['codr0']
        bpa = (bpot(cod1) - bpot(cod0)) / pr
        g1 = _grad_xy(bpot, cod1)
        g0 = _grad_xy(bpot, cod0)
        dpsi = g1[0] * tr2[0, :] + g1[1] * tr2[2, :]
        dpsi[0] -= g0[0]
        dpsi[2] -= g0[1]
        dbpa = dpsi / pr
        dbpa[5] -= bpa / pr
        blx = blx + bpa * pxh
        bly = bly + bpa * pyh
        blz = blz + bpa * pzh
        dblx = dblx + dbpa * pxh + bpa * dpxh
        dbly = dbly + dbpa * pyh + bpa * dpyh
        dblz = dblz + dbpa * pzh + bpa * dpzh
    ct = h1 * GSPIN
    dct = ct * dh1r
    ct += 1.0
    gx = ct * btx + CL * blx
    gy = ct * bty + CL * bly
    gz = ct * btz + CL * blz
    dgx = ct * dbtx + CL * dblx
    dgx[5] += dct * btx
    dgy = ct * dbty + CL * dbly
    dgy[5] += dct * bty
    dgz = ct * dbtz + CL * dblz
    dgz[5] += dct * btz
    srot = st.srot
    srot[0, 3:9] += dgx @ transr
    srot[1, 3:9] += dgy @ transr
    srot[2, 3:9] += dgz @ transr
    g = math.sqrt(gx * gx + gy * gy + gz * gz)
    if g != 0.0:
        bt = math.sqrt(btx * btx + bty * bty + btz * btz)
        sinu, cosu, dcosu = xsincos(g)
        sinu = sinu / g
        sx, sy, sz = srot[0].copy(), srot[1].copy(), srot[2].copy()
        sw = (sx * gx + sy * gy + sz * gz) / g
        st.gintd = st.gintd + sw[0:3] * sflc() * anp * (bt * h1 * p / al1) ** 2
        _w = sflc() * anp * (bt * h1 * p / al1) ** 2
        st.gintd2 = st.gintd2 + np.outer(sw[0:3], sw[0:3]) * _w       # v2 (A2): 磁場の向きの符号に依らない緩和率 Σ r (b·n)² 用 (k64 の結果には使わない)
        # v2 (P1): Sokolov-Ternov の厳密な緩和率 Σ r (1 − 2/9 (n·v)²) 用 (gsum = Σ r、gintd2v = Σ r (v̂ v̂ᵀ)。v̂ は刻みの運動量の向き。k64 の結果には使わない)
        _pv = np.array([pxh, pyh, pzh])
        _pv = _pv / np.linalg.norm(_pv)
        _sv = (sx * _pv[0] + sy * _pv[1] + sz * _pv[2])[0:3]
        st.gsum = getattr(st, 'gsum', 0.0) + _w
        st.gintd2v = getattr(st, 'gintd2v', np.zeros((3, 3))) + np.outer(_sv, _sv) * _w
        sw = -sw * dcosu / g
        sux = sy * gz - sz * gy
        suy = sz * gx - sx * gz
        suz = sx * gy - sy * gx
        srot[0] = cosu * sx + sinu * sux + sw * gx
        srot[1] = cosu * sy + sinu * suy + sw * gy
        srot[2] = cosu * sz + sinu * suz + sw * gz
        if STEP_ACC and bt > 0.0:
            _acc_step(st, c, btx / bt, bty / bt, btz / bt, _w)
    srot[0] = cp * srot[0] + sp * srot[2]
    srot[2] = (-sp * srot[0] + srot[2]) / cp       # 更新後の 1 行目を使う (K64 のまま)
    if eout is not None and _compat0.spin_edge_split():
        # v2 (B137): 出口の端の ψ は、刻みの本体と座標系の回転の後の独立な回転 (刻みの終わりの運動量の向き)。ψ の座標の微分は刻みの終わりの座標 → tr2 で基準点に戻す
        pv_out = np.array([px, py, pz]) / pr
        dpv_out = np.array([tr2[1, :], tr2[3, :], ddpz]) / pr
        dpv_out[:, 5] -= pv_out / pr
        _edge_rot(st, eout[0], eout[1] @ tr2, pv_out, dpv_out, pr, h1, dh1r, transr)


def _acc_step(st, c, bhx, bhy, bhz, w):
    """v2 (P1): 1 刻みの積算 (要素の入口の座標・スピンの系で書く)。放射の跳び ξ (刻みの出口の座標、共分散 b_s) が動かすスピン軸は N_s ξ、
    N_s T_e = R_e N_in + [n_e]× D_e (T_e = 入口から刻みの出口までの写像、R_e・D_e = そこまでのスピンの回転とスピン軌道結合、N_in・n_in = 要素の入口の値)。
    R_eᵀ N_s = (N_in + [n_in]× E_s) T_e⁻¹ (E_s = R_eᵀ D_e) なので、脱偏極率 ½ Σ_s tr(N_s⊥ b_s N_s⊥ᵀ) は
        rt_M = Σ M_s,  rt_F = Σ M_s E_sᵀ,  rt_G = Σ E_s M_s E_sᵀ   (M_s = T_e⁻¹ b_s T_e⁻ᵀ)
    から、要素の入口の N_in・n_in を使って後から作れる (slim_depolarization)。運動論的偏極 Σ w b̂·m_s (m_s = N_s の δ 列 = ∂n/∂δ。位置・角度を固定した δ の偏微分) は
        kin_A = Σ w (R_eᵀ b̂) (T_e⁻¹ e_δ)ᵀ,  kin_Y = Σ w R_eᵀ ((D_e T_e⁻¹ e_δ) × b̂)   →   K = Σ A ∘ N_in + n_in·Y
    b̂ は刻みの磁場 (B⊥) の向き、w は gintd と同じ重み。"""
    from .emit6d import CUU
    srot = st.srot
    Re = srot[:, 0:3].copy()
    De = srot[:, 3:9].copy()
    try:
        Tei = np.linalg.inv(st.T + st.D)
    except np.linalg.LinAlgError:
        return
    bh = np.array([bhx, bhy, bhz])
    t = Tei[:, 5]
    rb = Re.T @ bh
    st.kin_A = getattr(st, 'kin_A', 0.0) + w * np.outer(rb, t)
    st.kin_Y = getattr(st, 'kin_Y', 0.0) + w * (Re.T @ np.cross(De @ t, bh))
    px, py, pxi, pyi, bz = c['px'], c['py'], c['pxi'], c['pyi'], c['bzh']
    dee = c['anp'] * c['uc'] ** 2 * CUU
    pxm = pxi + px
    pym = pyi + py
    B = st.B
    bm = np.zeros((6, 6))
    bm[1, 1] = (B[1, 1] + bz * (2.0 * B[1, 2] + bz * B[2, 2]) + (pxm ** 2 + pxi ** 2 + px ** 2) / 6.0) * dee
    bm[3, 1] = bm[1, 3] = (B[3, 1] - bz * (B[0, 1] - B[2, 3] + bz * B[0, 2]) + (pxm * pym + pxi * pyi + px * py) / 6.0) * dee
    bm[3, 3] = (B[3, 3] + bz * (-2.0 * B[0, 3] + bz * B[0, 0]) + (pym ** 2 + pyi ** 2 + py ** 2) / 6.0) * dee
    bm[5, 1] = bm[1, 5] = pxm * dee * 0.5
    bm[5, 3] = bm[3, 5] = pym * dee * 0.5
    bm[5, 5] = dee
    M = Tei @ bm @ Tei.T
    E = Re.T @ De
    st.rt_M = getattr(st, 'rt_M', 0.0) + M
    st.rt_F = getattr(st, 'rt_F', 0.0) + M @ E.T
    st.rt_G = getattr(st, 'rt_G', 0.0) + E @ M @ E.T


# ---------------------------------------------------------------------------
# tchge.f: 座標回転と tbrote の srot 更新
# ---------------------------------------------------------------------------

def tchge_rotate_spin(st, th, ent):
    """Rotate the spin matrix srot at element entrance (inverse rotation at exit), as in tchge.
    tchge: ent → srot(1,:)= c·sx − s·srot(2,:); srot(2,:)= s·sx + c·srot(2,:)。出口は逆回転"""
    if th == 0.0:
        return
    ct, sn = math.cos(th), math.sin(th)
    srot = st.srot
    sx = srot[0].copy()
    if ent:
        srot[0] = ct * sx - sn * srot[1]
        srot[1] = sn * sx + ct * srot[1]
    else:
        srot[0] = ct * sx + sn * srot[1]
        srot[1] = -sn * sx + ct * srot[1]


def tbrote_rr(alg, phig, dtheta, dchi2, ent):
    """Compute the rr rotation matrix of tbrote as the difference from the identity (diagonal in cos-1 form).
    tbrote の rr (単位行列からの差。対角は cos−1 の形)"""
    from .lattice import _xsincos
    if phig == 0.0:
        c0, s0 = 1.0, 0.0
    else:
        c0, s0 = math.cos(phig), math.sin(phig)
    s2, _, c2, dc2 = _xsincos(dchi2)
    st, _, ct, dct = _xsincos(dtheta)
    rr = np.zeros((3, 3))
    rr[0, 0] = dct * c0 ** 2 + dc2 * s0 ** 2 + c0 * s0 * s2 * st
    rr[1, 1] = c2 * dct + dc2
    rr[2, 2] = dc2 * c0 ** 2 + dct * s0 ** 2 - c0 * s0 * s2 * st
    if ent:
        rr[0, 1] = s0 * s2 - c0 * c2 * st
        rr[0, 2] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[1, 0] = -ct * s0 * s2 + c0 * st
        rr[1, 2] = c0 * ct * s2 + s0 * st
        rr[2, 0] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[2, 1] = -c0 * s2 - c2 * s0 * st
    else:
        rr[1, 0] = s0 * s2 - c0 * c2 * st
        rr[2, 0] = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        rr[0, 1] = -ct * s0 * s2 + c0 * st
        rr[2, 1] = c0 * ct * s2 + s0 * st
        rr[0, 2] = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        rr[1, 2] = -c0 * s2 - c2 * s0 * st
    return rr


def tbrote_spin(st, alg, phig, dtheta, dchi2, ent):
    rr = tbrote_rr(alg, phig, dtheta, dchi2, ent)
    srot = st.srot
    from . import compat
    if compat.v2():
        # v2 (B56): 行列積 (I + rr)·srot。回転 (直交行列) なので srot(:,1:3) の直交性が保たれる
        st.srot[:, :] = srot + rr @ srot
        return
    for i in range(9):                            # K64: 行を逐次更新 (2 行目以降は更新後の値を使う)
        srot[0, i] = rr[0, :] @ srot[:, i] + srot[0, i]
        srot[1, i] = rr[1, :] @ srot[:, i] + srot[1, i]
        srot[2, i] = rr[2, :] @ srot[:, i] + srot[2, i]


def tchge_spin(st, dx0, dy, dz, theta, dtheta, alg, phig, ent, dchi2=0.0):
    """Spin (calpol) part of tchge; the cod update is done separately by `lattice._tchge_cod`.
    tchge の calpol 部 (cod の更新は lattice._tchge_cod が別にやる)"""
    th = theta if phig != 0.0 else theta + dtheta
    if ent:
        tchge_rotate_spin(st, th, True)
        if th != 0.0:
            dx = math.cos(th) * dx0 - math.sin(th) * dy
        else:
            dx = dx0
        if dx != 0.0 or dz != 0.0 or dtheta != 0.0 or dchi2 != 0.0:
            tbrote_spin(st, alg, phig, dtheta, dchi2, True)
    else:
        if th != 0.0:
            dx = math.cos(th) * dx0 - math.sin(th) * dy
        else:
            dx = dx0
        if dx != 0.0 or dz != 0.0 or dtheta != 0.0 or dchi2 != 0.0:
            tbrote_spin(st, alg, phig, dtheta, dchi2, False)
        tchge_rotate_spin(st, th, False)


# ---------------------------------------------------------------------------
# spnorm
# ---------------------------------------------------------------------------

def _abs3(a, b, c):
    return math.hypot(a, math.hypot(b, c))


def spnorm(srot, gintd):
    """Orthogonalize srot(:,1:3); return the invariant spin axis sps(:,1), precession angle smu and damping rate sdamp (`spnorm`).
    spnorm: srot(:,1:3) を直交化し、スピン歳差の不変軸 sps(:,1) と歳差角 smu、減衰率 sdamp = n·gintd を返す。srot は更新される"""
    from .emit6d import _tsvd_tsolva
    srot[:, 1] = srot[:, 1] / np.linalg.norm(srot[:, 1])
    s = srot[0, 0] * srot[0, 1] + srot[1, 0] * srot[1, 1] + srot[2, 0] * srot[2, 1]
    srot[:, 0] = srot[:, 0] - s * srot[:, 1]
    s = _abs3(srot[0, 0], srot[1, 0], srot[2, 0])
    srot[:, 0] = srot[:, 0] / s
    srot[:, 2] = np.cross(srot[:, 0], srot[:, 1])
    sps = np.zeros((3, 3))
    sps[0, 0] = srot[1, 2] - srot[2, 1]
    sps[1, 0] = srot[2, 0] - srot[0, 2]
    sps[2, 0] = srot[0, 1] - srot[1, 0]
    s = _abs3(sps[0, 0], sps[1, 0], sps[2, 0])
    smin = 1e-4
    if s < smin:
        # teigen: 実固有値に対応する固有ベクトル
        w, v = np.linalg.eig(srot[:, 0:3])
        for i in range(3):
            if abs(w[i].imag) < 1e-12:
                vv = np.real(v[:, i])
                sps[:, 0] = vv
                s = _abs3(vv[0], vv[1], vv[2])
                break
    sps[:, 0] = sps[:, 0] / s
    dr = sps[:, 0] - srot[:, 0] * sps[0, 0] - srot[:, 1] * sps[1, 0] - srot[:, 2] * sps[2, 0]
    a = srot[:, 0:3].copy()
    a[0, 0] -= 1.0
    a[1, 1] -= 1.0
    a[2, 2] -= 1.0
    dsps = np.asarray(_tsvd_tsolva(a, dr, 1e-8))
    sps[:, 0] = sps[:, 0] + dsps
    sps[:, 0] = sps[:, 0] / np.linalg.norm(sps[:, 0])
    if abs(min(sps[:, 0])) > abs(max(sps[:, 0])):
        sps[:, 0] = -sps[:, 0]
    sps[:, 1] = 0.0
    if abs(sps[0, 0]) > abs(sps[1, 0]):
        sps[1, 1] = 1.0
        s = sps[1, 0]
    else:
        sps[0, 1] = 1.0
        s = sps[0, 0]
    sps[:, 1] = sps[:, 1] - s * sps[:, 0]
    sps[:, 1] = sps[:, 1] / np.linalg.norm(sps[:, 1])
    sps[:, 2] = np.cross(sps[:, 0], sps[:, 1])
    spsa1 = srot[:, 0] * sps[0, 1] + srot[:, 1] * sps[1, 1] + srot[:, 2] * sps[2, 1]
    cm = float(spsa1 @ sps[:, 1])
    sm = float(spsa1 @ sps[:, 2])
    smu = math.atan2(-sm, cm)
    sdamp = float(sps[:, 0] @ gintd)
    return sps, smu, sdamp


# ---------------------------------------------------------------------------
# srequpol / spdepol
# ---------------------------------------------------------------------------

def _serot(xx, xy, yy, c1, c2, d1, d2, e1, e2):
    tx = 0.5 * math.atan2(2.0 * xy, xx - yy)
    c = math.cos(tx)
    s = math.sin(tx)
    em1 = c ** 2 * xx + s ** 2 * yy + 2.0 * c * s * xy
    em2 = c ** 2 * yy + s ** 2 * xx - 2.0 * c * s * xy
    c1 = c * c1 + s * c2
    c2 = (-s * c1 + c2) / c                       # 更新後の c1 (K64 のまま)
    d1 = c * d1 + s * d2
    d2 = (-s * d1 + d2) / c
    e1 = c * e1 + s * e2
    e2 = (-s * e1 + e2) / c
    return c1, c2, d1, d2, e1, e2, em1, em2


def _cexp1(z):
    z = complex(z)
    if abs(z) < 2e-4:
        return z * (1.0 + 0.5 * z * (1.0 + z / 3.0 * (1.0 + z / 4.0 * (1.0 + z / 5.0))))
    return complex(np.exp(z) - 1.0)


def _cerfc(z):
    """gamma.f cerfc: |z|>18 は漸近展開、そうでなければ erfc"""
    z = complex(z)
    if abs(z) > CEZLIM:
        x = 0.5 / (z * z)       # z ** 2 は Python 3.14 では桁あふれで OverflowError (3.12 以前は inf・nan)。z * z は同じ値
        base = 1.0 if z.real == 0.0 else 1.0 - math.copysign(1.0, z.real)
        with np.errstate(all='ignore'):
            return base + (1.0 - x * (1.0 - 3.0 * x * (1.0 - 5.0 * x * (1.0 - 7.0 * x * (1.0 - 9.0 * x)))))\
                * np.exp(-0.5 / x) / z / M_SQRTPI
    from scipy.special import erfc
    return complex(erfc(z))


def _mat(rows):
    return np.array(rows, dtype=complex)


def _indn(i, i1, i2, idx, imx, is_):
    ind = [0] * 13
    ind[2 * i - 2] = i1
    ind[2 * i - 1] = i2
    ind[i + 5] = idx
    ind[i + 8] = imx
    ind[12] = is_
    return tuple(ind)


def _erfcx(z):
    """exp(z²)·erfc(z) (複素数、桁あふれなし) = w(i z) (Faddeeva 関数)"""
    from scipy.special import wofz
    return complex(wofz(1j * complex(z)))


def _spintrm(rm, dx, amx, ams, ssprd):
    from . import compat
    _v2 = compat.v2()
    for ind in list(rm.keys()):
        mat = rm[ind]
        cphi = complex(-(ind[6] * dx[0] + ind[7] * dx[1] + ind[8] * dx[2]),
                       ind[9] * amx[0] + ind[10] * amx[1] + ind[11] * amx[2] + ind[12] * ams)
        with np.errstate(all='ignore'):
            c1 = _cexp1(cphi)
            cm0 = -mat / c1 if c1 != 0 else -mat * complex(float('inf'), 0.0)
            a = math.sqrt(2.0) * abs(ind[12]) * ssprd
            if a == 0.0:
                cm = cm0
            else:
                cb = -_cexp1(-cphi) / a
                if not (math.isfinite(cb.real) and math.isfinite(cb.imag)):
                    # 桁あふれ (不安定なリング)。Python 3.12 以前は inf を含む複素数の割り算が nan+nanj になる。3.14 (C99 Annex G) は inf を残し、
                    # 続く `_cerfc` が 0 割りになるので、3.12 と同じ nan+nanj にそろえる (k64 の値は 3.12 と同じ)
                    cb = complex(float('nan'), float('nan'))
                if _v2:
                    # v2 (B55): exp(cb²)·erfc(cb) を桁あふれしない形 (Faddeeva 関数 w(i·cb)) で計算する。
                    # K64 は |cb|>18 で exp(cb²-cφ) → 0 と erfc の漸近式の exp(-cb²) → ∞ を別々に計算して NaN になる
                    cm = -mat * M_2PI * M_SQRTPI * np.exp(-cphi) * _erfcx(cb) / a
                else:
                    cm = -mat * M_2PI * M_SQRTPI * np.exp(cb * cb - cphi) * _cerfc(cb) / a
        rm[ind] = cm0
        k0 = ind[:6] + (0,) * 7
        if k0 in rm:
            rm[k0] = rm[k0] + cm
        else:
            rm[k0] = cm.copy() if hasattr(cm, 'copy') else cm


def _spdotrm(rma, rmb, rmc):
    for ia, ma in rma.items():
        for ib, mb in rmb.items():
            k = tuple(x + y for x, y in zip(ia, ib))
            t = ma @ mb
            if k in rmc:
                rmc[k] = rmc[k] + t
            else:
                rmc[k] = t


def spdepol(gxr, gxi, gyr, gyi, gzr, gzi, e1, e2, dx, ssprd, amx, ams):
    """Compute the second-order perturbative depolarization matrix (`spdepol`) and return rmd(3,3,1) = rm1.
    spdepol: 2 次までの摂動 (K64 は spcalcres を k=2 だけ呼ぶ)。rmd(3,3,1) = rm1 を返す"""
    rm1 = {}
    cI = 1j
    for i in (1, 2, 3):
        gx = complex(gxr[i - 1], gxi[i - 1])
        gxc = gx.conjugate()
        gy = complex(gyr[i - 1] + gzi[i - 1], gzr[i - 1] - gyi[i - 1])
        gyc = gy.conjugate()
        gz = complex(gyr[i - 1] - gzi[i - 1], gzr[i - 1] + gyi[i - 1])
        gzc = gz.conjugate()
        m1 = 0.25 * gzc * np.array([0, -cI, 1, cI, 0, 0, -1, 0, 0], dtype=complex).reshape((3, 3), order='F')
        m2 = 0.5 * gxc * np.array([0, 0, 0, 0, 0, -1, 0, 1, 0], dtype=complex).reshape((3, 3), order='F')
        m3 = 0.25 * gy * np.array([0, cI, 1, -cI, 0, 0, -1, 0, 0], dtype=complex).reshape((3, 3), order='F')
        rm1[_indn(i, 0, 1, 1, -1, -1)] = m1
        rm1[_indn(i, 0, 1, 1, -1, 0)] = m2
        rm1[_indn(i, 0, 1, 1, -1, 1)] = m3
        rm1[_indn(i, 1, 0, 1, 1, -1)] = m3.conj()
        rm1[_indn(i, 1, 0, 1, 1, 0)] = m2.conj()
        rm1[_indn(i, 1, 0, 1, 1, 1)] = m1.conj()
    rmi1 = {k: v.copy() for k, v in rm1.items()}
    _spintrm(rmi1, dx, amx, ams, ssprd)
    # spcalcres(m=2)
    rm2 = {}
    _spdotrm(rm1, rm1, rm2)
    for k in rm2:
        rm2[k] = rm2[k] / 2.0
    rmi2 = {k: v.copy() for k, v in rm2.items()}
    _spdotrm(rm1, rmi1, rmi2)
    _spintrm(rmi2, dx, amx, ams, ssprd)
    rmd = np.zeros((3, 3))
    z = np.zeros((3, 3), dtype=complex)
    for i in (1, 2, 3):
        a20 = rmi2.get(_indn(i, 2, 0, 0, 0, 0), z)
        a21 = rmi2.get(_indn(i, 1, 1, 0, 0, 0), z)
        se12 = (e1[i - 1] + e2[i - 1]) * 0.5
        de12 = (e1[i - 1] - e2[i - 1]) * 0.5
        rmd = rmd + se12 * np.real(a21) + 2.0 * de12 * np.real(a20)
    return rmd


def srequpol(srot, sps, params, demit, sdamp, r):
    """Compute the equilibrium polarization equpol(3) and depolarization matrix rm1(3,3) (`srequpol`).
    srequpol: 平衡偏極 equpol(3) と脱偏極行列 rm1(3,3)。
    params: dict(emit=(ex,ey,ez), damp=(3), tunes=(3), smu)。demit は規格化モード座標のビーム行列 (21 個の packed)"""
    smu = params['smu']
    drot = sps.T @ srot[:, 3:9] @ r
    c1, c2, c3, c4, c5, c6 = drot[0]
    d1, d2, d3, d4, d5, d6 = drot[1]
    e1, e2, e3, e4, e5, e6 = drot[2]
    c1, c2, d1, d2, e1, e2, dex1, dex2 = _serot(demit[0], demit[1], demit[2], c1, c2, d1, d2, e1, e2)
    c3, c4, d3, d4, e3, e4, dey1, dey2 = _serot(demit[5], demit[8], demit[9], c3, c4, d3, d4, e3, e4)
    c5, c6, d5, d6, e5, e6, dez1, dez2 = _serot(demit[14], demit[19], demit[20], c5, c6, d5, d6, e5, e6)
    emit1 = params['emit']
    damp = np.abs(np.asarray(params['damp'], dtype=float))
    ssprd = math.sqrt(float(np.dot(drot[0] ** 2, [emit1[0], emit1[0], emit1[1], emit1[1], emit1[2], emit1[2]])))
    amu = np.asarray(params['tunes'], dtype=float) * M_2PI
    rm1 = spdepol([c1, c3, c5], [c2, c4, c6], [d1, d3, d5], [d2, d4, d6], [e1, e3, e5], [e2, e4, e6],
                  [dex1, dey1, dez1], [dex2, dey2, dez2], damp, ssprd, amu, smu)
    with np.errstate(all='ignore'):
        equpol = np.array([PST * sdamp / (sdamp - rm1[0, 0]), 0.0, 0.0])
    return equpol, rm1


def _compat_v2():
    from . import compat
    return compat.v2()


def _cross_mat(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def sideband_depolarization(secs, n0, s0, nsamp=1024, tol=1e-7, kmax=300000, modulate=True, seed=20261118, C6=None, us=None):
    """Compute the per-turn depolarization rate including synchrotron sideband (spin phase modulation) effects (v2 Polarization, B3).
    シンクロトロン側帯 (スピンの位相変調) を入れた 1 周あたりの脱偏極率 (v2 の Polarization、B3)。独立な導出 (docs/history/v2next_B3.md §2)。

    線形の式 (slim_depolarization) は、スピン軸の場 n(u) = n0 + N u を軌道 u の 1 次で置き、放射の跳び ξ (共分散 B_i) あたり ½|N ξ|² だけ偏極が減ると
    数える。スピンの 1 周の歳差の位相は 2π ν(u) = 2π ν0 + P·u (P はエネルギー (δ) 依存が主) で、シンクロトロン振動で δ が振動すると位相が変調され
    (変調指数 ~ aγ σδ/νs)、n(u) は u について非線形になる (周波数 ν0 ± kνs の側帯の共鳴)。ここでは変調を厳密に取り (u の 1 次の結合は残す)、
    n 軸の場の勾配 ∂n/∂u を、軌道の分布 (ガウス分布、`s0` = 入口のビーム行列) の各点で求めて平均する。

    導出: スピンのずれを複素数 α = s·(m + i l) (n0 のまわりの回転座標、m×l = n) で表すと、要素 i の (D_i u による余分の回転) で
        α' = e^{i p_i·u} α + g_i·u,   p_i = −n_i·D_i,  g_i = i (m_i + i l_i)·D_i
    (K64 の回転の向き s' = (I − [D u]×) R s)。1 周では α_{n+1} = e^{iμ} (e^{i P·u_n} α_n + G·u_n) (P = Σ p_i Φ_i、G = Σ g_i Φ_i、Φ_i は入口から要素 i への軌道の写像)。
    不変な場 (放射の跳びが無いとき保存される量 α − α*(u))は未来の軌道について和を取って
        α*(u) = −Σ_{k≥0} e^{−iμk} exp(−i S_k) G·M^k u,   S_k = Σ_{t=0}^{k} P·M^t u
    (M は減衰つきの 1 周写像。|eig M| < 1 で収束)。変調を切る (P=0) と N = −G (1 − e^{−iμ}M)^{-1}(線形の周期解) に一致する。
    ξ の跳びでの偏極の減少は ½|∂α*/∂u · ξ|² で、勾配は
        ∂α*/∂u = −Σ_k e^{−iμk} e^{−iS_k} [G M^k − i (G M^k u) Q_k],   Q_k = Σ_{t≤k} P M^t
    要素 i の出口の勾配は、鎖 G_{i+1} = [e^{i p_i u_i} G_i + i α*_i e^{i p_i u_i} p_i + g_i] T_i^{-1}、α*_{i+1} = e^{i p_i u_i} α*_i + g_i u_i で伝える。
    脱偏極率 rate = ½ Σ_i ⟨G_{i+1} B_i G_{i+1}^†⟩ (⟨⟩ は u の分布での平均: Sobol 列の 1024 点)。modulate=False で線形の式 (slim_depolarization) に一致する。
    secs = [(T_i, B_i, R_i, D_i)]、n0 は入口のスピン軸、s0 は入口のビーム行列。
    戻り dict(rate, spread (½⟨|α*_0|²⟩ = 入口での n の広がり), kturns (使った K))
    限界: 結合 b の u の 2 次の項 (変調 × 結合) は入れない (相対 ~aγδ·θ)。|eig M| が 1 に近い (放射減衰が弱い) ときは K が大きく遅い。
    縦の減衰が無い (RF OFF) 場合は呼ばない (和が収束しない)。
    v2 (P1): C6 (1 周の写像の δ 微分 dM/dδ、分散軌道に沿った 6x6) が渡されたとき、ベータトロンの固有モード j の位相がシンクロトロン振動の δ_t で変調される
    (クロマティシティ) ことを入れる: u_{k,j} = c_j λ_j^k exp(i ξ_j Δ_k)、Δ_k = Σ_{t<k} δ_t (δ_t = 線形のシンクロトロン解)、ξ_j = Im[(V⁻¹ C6 V)_jj / λ_j] (rad/δ)。
    スピンの歳差の位相 S_k と n 軸の勾配 ∂α*/∂u も同じ軌道で組み直す (∂Δ_k/∂u を通る項が増える)。スピンチューンがベータトロンの側帯 νx + kνs に近い共鳴
    (4.2 GeV・12 MV の 12 セル FODO: νspin ≈ νx + 2νs) では、スピンと軌道の位相の相対的な変調 (aγ − ξ) が側帯の重みを変える (記録 v2next_P1.md §4)。"""
    N = len(secs)
    m0, l0, ms, ls, ns, mu = _spin_frames(secs, n0)
    p = [-(ns[i] @ secs[i][3]) for i in range(N)]
    g = [1j * ((ms[i] + 1j * ls[i]) @ secs[i][3]) for i in range(N)]
    if not modulate:
        p = [np.zeros(6) for _ in range(N)]
    Phi = np.eye(6)
    Pt = np.zeros(6, dtype=complex)
    Gt = np.zeros(6, dtype=complex)
    for i in range(N):
        Pt = Pt + p[i] @ Phi
        Gt = Gt + g[i] @ Phi
        Phi = secs[i][0] @ Phi
    w, V = np.linalg.eig(Phi)
    lam = float(np.abs(w).max())
    if not lam < 1.0 - 1e-9:
        raise ValueError('側帯の式: 1 周写像に減衰が無い (RF OFF か放射が無い)')
    Vi = np.linalg.inv(V)
    K = int(min(kmax, math.log(tol) / math.log(lam))) + 1
    k = np.arange(K)
    lamk = w[None, :] ** k[:, None]
    GV = Gt @ V
    PV = Pt @ V
    geo = (1.0 - w[None, :] ** (k[:, None] + 1)) / (1.0 - w[None, :])
    GMk = (lamk * GV[None, :]) @ Vi
    Qk = (geo * PV[None, :]) @ Vi
    xi = np.zeros(6)
    if C6 is not None and modulate:
        dlog = np.einsum('ij,jk,ki->i', Vi, np.asarray(C6, dtype=complex), V) / w
        xi = np.imag(dlog)
        ws = np.abs(V[4, :]) ** 2 + np.abs(V[5, :]) ** 2
        for j in np.argsort(-ws)[:2]:
            xi[j] = 0.0                                  # シンクロトロンのモード (RF の非線形・α₂ は入れない)
    ws2 = np.abs(V[4, :]) ** 2 + np.abs(V[5, :]) ** 2
    isy = np.argsort(-ws2)[:2]
    nus = float(np.abs(np.angle(w[isy[0]])) / (2.0 * math.pi))
    mod_spin = float(math.sqrt(max(0.0, float(np.real(Pt) @ s0 @ np.real(Pt)))) / (2.0 * math.pi * nus)) if nus > 0.0 else float('nan')   # スピンの位相の変調指数
    mod_chrom = float(np.abs(xi).max() * math.sqrt(max(0.0, float(s0[5, 5]))) / (2.0 * math.pi * nus)) if nus > 0.0 else float('nan')    # ベータトロンの位相の変調指数
    if C6 is not None and mod_chrom < CHROM_MIN:
        xi = np.zeros(6)                                 # 効果が 1% 未満: 入れない (計算も速い)
    if K > CHROM_KMAX:
        xi = np.zeros(6)
    chrom = bool(np.abs(xi).max() > 0.0)
    from scipy.stats import norm, qmc
    m2 = max(1, int(round(math.log2(nsamp))))
    z = norm.ppf(np.clip(qmc.Sobol(d=6, scramble=True, seed=seed).random_base2(m2), 1e-12, 1 - 1e-12))
    ev, evec = np.linalg.eigh(0.5 * (s0 + s0.T))
    L = evec * np.sqrt(np.clip(ev, 0.0, None))
    us = z @ L.T if us is None else np.asarray(us, dtype=float)       # us: 試料の点を直接渡す (導出の検算用。戻りの G0・a0)
    ns_ = us.shape[0]
    phs = np.exp(-1j * mu * k)
    G0 = np.zeros((ns_, 6), dtype=complex)
    a0 = np.zeros(ns_, dtype=complex)
    batch = max(1, int(4.0e6 // K))                # (batch × K) の複素配列を数個: 4e6 × 16 B × 4 ≈ 250 MB 以下
    if chrom:
        geo0 = (1.0 - w[None, :] ** k[:, None]) / (1.0 - w[None, :])               # Σ_{t<k} λ^t
        dDelta = np.real((geo0 * V[5, :][None, :]) @ Vi)                           # ∂Δ_k/∂u (K, 6)
        batch = max(1, int(4.0e5 // K))                # (batch × K × 6) の複素配列を 6 個前後
    for b in range(0, ns_, batch):
        u = us[b:b + batch]
        c = u @ Vi.T
        if chrom:
            Delta = np.real((c * V[5, :][None, :]) @ geo0.T)                       # Δ_k = Σ_{t<k} δ_t (s, K)
            Eh = lamk[None, :, :] * np.exp(1j * xi[None, None, :] * Delta[:, :, None])
            Ek = Eh * c[:, None, :]                                                # c_j λ_j^k e^{iξ_j Δ_k} (s, K, 6)
            gk = Ek @ GV
            Sk = np.cumsum(np.real(Ek @ PV), axis=1)
            ph = phs[None, :] * np.exp(-1j * Sk)
            a0[b:b + batch] = -(ph * gk).sum(1)
            bg = 1j * (Ek * (GV * xi)[None, None, :]).sum(-1)
            bp = 1j * (Ek * (PV * xi)[None, None, :]).sum(-1)
            dg = (Eh * GV[None, None, :]) @ Vi + bg[:, :, None] * dDelta[None, :, :]
            dS = np.cumsum(np.real((Eh * PV[None, None, :]) @ Vi + bp[:, :, None] * dDelta[None, :, :]), axis=1)
            G0[b:b + batch] = -(ph[:, :, None] * dg).sum(1) + 1j * ((ph * gk)[:, :, None] * dS).sum(1)
            continue
        gk = (c * GV[None, :]) @ lamk.T
        if modulate:
            ph = phs[None, :] * np.exp(-1j * ((c * PV[None, :]) @ geo.T))
            a0[b:b + batch] = -(ph * gk).sum(1)
            G0[b:b + batch] = -(ph @ GMk) + (1j * (ph * gk)) @ Qk
        else:
            ph = np.broadcast_to(phs[None, :], gk.shape)
            a0[b:b + batch] = -(ph * gk).sum(1)
            G0[b:b + batch] = -(ph @ GMk)
    u = us.copy()
    G = G0.copy()
    a = a0.copy()
    rate = 0.0
    for i in range(N):
        T, B, _R, _D = secs[i]
        e = np.exp(1j * (u @ p[i])) if modulate else np.ones(ns_)
        Gn = (e[:, None] * G + ((1j * e * a)[:, None] * p[i][None, :] if modulate else 0.0) + g[i][None, :]) @ np.linalg.inv(T)
        a = e * a + u @ g[i]
        u = u @ T.T
        G = Gn
        rate += 0.5 * float(np.mean(np.einsum('si,ij,sj->s', G, B, G.conj()).real))
    out = dict(rate=rate, spread=0.5 * float(np.mean(np.abs(a0) ** 2)), kturns=K, xi=xi, chromatic=chrom, G0=G0, a0=a0, mod_spin=mod_spin, mod_chrom=mod_chrom, nus=nus)
    # v2 (B341): 和の打ち切り (K が kmax で止まる: |λ|max^K > tol) の補正用に、同じ K・同じ試料の線形の和 (modulate=False) を添える。
    # 放射減衰が弱いリング (ASTRID2: |λ|max = 0.9999927、tol 1e-7 に 2.2e6 周が要るのに kmax = 3e5) では、打ち切った和が e^{-K/τ} ≈ 0.1 の大きさの残りを持ち、
    # 残りの位相 e^{-iμK}λ^K が 1e-6 のチューンの変化で回る (線形の和そのものが厳密な値から 6% ずれた)。厳密な線形の値は呼び出し側にあるので、
    # 側帯の効果は「変調あり − 変調なし (同じ K)」の差として足す (和の残りは変調の有無でほぼ共通で、差では打ち消す)
    out['truncated'] = bool(math.log(tol) / math.log(lam) > kmax)
    if modulate and out['truncated']:
        lin = sideband_depolarization(secs, n0, s0, nsamp=nsamp, tol=tol, kmax=kmax, modulate=False, seed=seed, C6=None, us=us)
        out['rate_lin_sameK'] = lin['rate']
        out['spread_lin_sameK'] = lin['spread']
    return out


def _spin_frames(secs, n0):
    """n0 に垂直な基底 (m0, l0) (m×l = n) を閉軌道のスピンの歳差 R_i で運んだ (m_i, l_i, n_i) と、1 周の歳差角 μ"""
    a = np.array([1.0, 0.0, 0.0]) if abs(n0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    m0 = a - (a @ n0) * n0
    m0 = m0 / np.linalg.norm(m0)
    l0 = np.cross(n0, m0)
    ms, ls, ns_ = [], [], []
    m, l, n = m0, l0, n0
    for _T, _B, R, _D in secs:
        n = R @ n
        m = R @ m
        l = R @ l
        ms.append(m)
        ls.append(l)
        ns_.append(n)
    mu = math.atan2(float(m @ l0), float(m @ m0))
    return m0, l0, ms, ls, ns_, mu


def _chromatic_matrix(els, ph, cod, trf0, rfsw, radcod, secs, eps=1e-4):
    """v2 (P1): 1 周の写像 (T + D) の δ 微分 dM/dδ (6x6)。閉軌道を分散軌道の向き (η_x, η_px, η_y, η_py, 0, 1) に ±eps だけずらした軌道に沿って
    要素ごとのヤコビアンを掛けた 1 周の写像の中心差分。ベータトロンのチューンの δ 依存 (クロマティシティ。自然クロマティシティ + 6 極) を含む"""
    from . import emit6d as _e6
    Phi = np.eye(6)
    for T, _B, _R, _D in secs:
        Phi = T @ Phi
    eta = np.linalg.solve(np.eye(4) - Phi[:4, :4], Phi[:4, 5])
    eta6 = np.array([eta[0], eta[1], eta[2], eta[3], 0.0, 1.0])
    old = _e6.POL_MODE
    _e6.POL_MODE = False
    Ms = []
    try:
        for sgn in (1.0, -1.0):
            c = np.array(cod, dtype=float) + sgn * eps * eta6
            M = np.eye(6)
            for el in els:
                st = _e6.tturne([el], c, ph, trf0=trf0, rfsw=rfsw, radcod=radcod)
                M = (st.T + st.D) @ M
                c = st.cod.copy()
            Ms.append(M)
    finally:
        _e6.POL_MODE = old
    return (Ms[0] - Ms[1]) / (2.0 * eps)


def slim_depolarization(els, ph, cod, trf0=0.0, rfsw=False, radcod=False, beam_matrix=None, sideband=False):
    """Compute the per-turn depolarization rate from linear spin-orbit coupling (SLIM / Derbenev-Kondratenko term) for v2 Polarization.
    線形のスピン軌道結合 (SLIM / Derbenev-Kondratenko の脱偏極項) による 1 周あたりの脱偏極率 (v2 の Polarization が使う)。

    要素 i ごとに emit6d.tturne (POL) で、軌道の行列 T_i (減衰込み)、その要素の放射の拡散行列 B_i (要素の出口での共分散)、
    スピンの回転 R_i、スピン軌道結合 D_i (軌道のずれ u による余分の回転ベクトル D_i u。K64 の回転の向きは ds = s × g なので
    s' = (1 − [D_i u]×) R_i s。追跡 (track6 のスピン) の 1 周の応答と 2〜6% で一致することを確かめた) を取る。
    スピンの周期解 n(u) = n0 + N u (N は 3x6) は  N_{i+1} = (R_i N_i + [n0_{i+1}]× D_i) T_i^{-1}、N_N = N_0 を満たす。
    放射で軌道が ξ (共分散 B_i) だけ跳ぶと、スピンは変わらず n が N ξ だけ動くので、偏極は 1 周あたり
        rate = (1/2) Σ_i tr(N⊥ B_i N⊥ᵀ)   (N⊥ は n0 に垂直な成分)
    だけ減る。これは DK の式の (11/18)|γ ∂n/∂γ|² の項と同じ (B_i の δδ 成分 = (55/24√3) r_e ƛ γ^5/|ρ|³ ds)。
    戻り dict(rate, n0, N0 (入口の ∂n/∂u)、rev (平衡分布での n の広がり (1/2) tr(N⊥ Σ N⊥ᵀ) を呼び出し側が作るための N0))。
    限界: 軌道の振幅について 1 次 (シンクロトロン側帯などの高次の共鳴は入らない)、運動論的偏極 (b·γ∂n/∂γ) の項は入らない"""
    from . import emit6d as _e6
    global STEP_ACC
    old = _e6.POL_MODE
    old_acc = STEP_ACC
    _e6.POL_MODE = True
    STEP_ACC = bool(STEP_RATE or KINETIC) and _compat_v2()
    secs = []
    acc = []
    c = np.array(cod, dtype=float)
    try:
        for el in els:
            st = _e6.tturne([el], c, ph, trf0=trf0, rfsw=rfsw, radcod=radcod)
            secs.append((st.T + st.D, st.B.copy(), st.srot[:, 0:3].copy(), st.srot[:, 3:9].copy()))
            acc.append(None if getattr(st, 'rt_M', None) is None else dict(M=st.rt_M, F=st.rt_F, G=st.rt_G, A=st.kin_A, Y=st.kin_Y))
            c = st.cod.copy()
    finally:
        _e6.POL_MODE = old
        STEP_ACC = old_acc
    rtot = np.eye(3)
    ttot = np.eye(6)
    for _T, _B, R, _D in secs:
        rtot = R @ rtot
        ttot = _T @ ttot
    if _compat_v2() and radcod:
        # v2 (P1): 光学が不安定なリング (1 周の写像 T + D の固有値の大きさが 1 を超える。放射減衰があれば安定なリングは全て < 1) の脱偏極率は意味がない
        # (線形の式の周期解が不安定な軌道の上の値になる: B3 の大きなリング A の縦ずれ 0.2・0.4 mm が x の固有値 1.2 で P = 0.92 を返していた)。NaN にする
        wmax = float(np.abs(np.linalg.eigvals(ttot)).max())
        if wmax > 1.0 + 1e-6:
            raise ValueError('unstable optics: max |eig(one-turn map)| = %.6f > 1' % wmax)
    w, v = np.linalg.eig(rtot)
    k = int(np.argmin(np.abs(w - 1.0)))
    n0 = np.real(v[:, k])
    n0 = n0 / np.linalg.norm(n0)
    if abs(min(n0)) > abs(max(n0)):
        n0 = -n0
    tinv = [np.linalg.inv(T) for T, _B, _R, _D in secs]

    def prop(N):
        n = n0
        out = []
        for (T, _B, R, D), ti in zip(secs, tinv):
            n = R @ n
            N = (R @ N + _cross_mat(n) @ D) @ ti
            out.append((N, n))
        return N, out
    cvec, _ = prop(np.zeros((3, 6)))
    a = np.zeros((18, 18))
    for j in range(18):
        e = np.zeros(18)
        e[j] = 1.0
        a[:, j] = (prop(e.reshape(3, 6))[0] - cvec).ravel()
    try:
        nn0 = np.linalg.solve(np.eye(18) - a, cvec.ravel()).reshape(3, 6)
    except np.linalg.LinAlgError:
        nn0 = np.linalg.lstsq(np.eye(18) - a, cvec.ravel(), rcond=None)[0].reshape(3, 6)
    _, out = prop(nn0)
    rate = 0.0
    for (_T, B, _R, _D), (N, n) in zip(secs, out):
        npp = (np.eye(3) - np.outer(n, n)) @ N
        rate += 0.5 * float(np.trace(npp @ B @ npp.T))
    out_ = dict(rate=rate, n0=n0, N0=nn0, rate_linear=rate, rate_el=rate, kin=0.0)
    if any(a is not None for a in acc):
        # v2 (P1): 放射の刻みごとのスピン軸の場で積分した脱偏極率と、運動論的偏極 K = Σ w b̂·∂n/∂δ (要素の入口の N_in・n_in と刻みごとの積算から)
        rs = 0.0
        kin = 0.0
        N_in, n_in = nn0, n0
        for a, (N_out, n_out) in zip(acc, out):
            if a is not None:
                Nt = (np.eye(3) - np.outer(n_in, n_in)) @ N_in
                cn = _cross_mat(n_in)
                rs += 0.5 * (float(np.trace(Nt @ a['M'] @ Nt.T)) + 2.0 * float(np.trace(Nt @ a['F'] @ cn.T)) + float(np.trace(cn @ a['G'] @ cn.T)))
                kin += float(np.sum(a['A'] * N_in)) + float(n_in @ a['Y'])
            N_in, n_in = N_out, n_out
        out_['rate_step'] = rs
        out_['kin'] = kin
        if STEP_RATE:
            out_.update(rate=rs, rate_linear=rs)
    if sideband and beam_matrix is not None:
        # v2 (B3): シンクロトロン側帯 (スピンの位相変調) を入れた脱偏極率。縦の減衰が無い場合は線形のまま
        try:
            C6 = None
            if CHROMATIC and _compat_v2() and len(els) <= CHROM_NEL:
                C6 = _chromatic_matrix(els, ph, cod, trf0, rfsw, radcod, secs)
            sb = sideband_depolarization([(T, B, R, D) for (T, B, R, D) in secs], n0, np.asarray(beam_matrix, dtype=float), C6=C6)
            # 要素の出口の N での値 sb['rate'] に、刻みごとの積分との差 (rate_step − rate_el: n が傾いたリングだけ 0 でない) を足す
            sb_rate = sb['rate']
            sb_spread = sb['spread']
            if sb.get('truncated') and 'rate_lin_sameK' in sb:
                # v2 (B341): 和が打ち切られている (K = kmax で未収束) ので、sb['rate'] の「線形の厳密な値 (rate_el)」からのずれを、同じ K の線形の和との差で置き換える
                sb_rate = rate + (sb['rate'] - sb['rate_lin_sameK'])
                npp0 = (np.eye(3) - np.outer(n0, n0)) @ nn0
                sb_spread = 0.5 * float(np.trace(npp0 @ np.asarray(beam_matrix, dtype=float) @ npp0.T)) + (sb['spread'] - sb['spread_lin_sameK'])
            out_.update(rate=sb_rate + (out_['rate'] - out_['rate_el']), spread_sb=sb_spread, sideband_kturns=sb['kturns'], sideband_truncated=bool(sb.get('truncated')),
                        chromatic_xi=sb['xi'], chromatic=sb['chromatic'], mod_spin=sb['mod_spin'], mod_chrom=sb['mod_chrom'])
        except Exception as exc:
            out_['sideband_error'] = str(exc)
    return out_


def spin_results(st, em, ph, els=None, cod=None, trf0=0.0, rfsw=False, radcod=False) -> dict:
    """Compute the spin results of tecalc in temit.f (spnorm, srequpol), with polarization time and, in v2, SLIM polarization.
    temit.f tecalc の calpol 部: srot1 = srot(:,1:3) → spnorm → srequpol。偏極時間 = 1/sdamp/revolution frequency [s]。
    v2 (B55) で els・cod が渡されたとき: 平衡偏極度は線形のスピン軌道結合 (slim_depolarization) の脱偏極率から
    P = (8/5√3)·sdamp/(sdamp + rate) とする (K64 の spdepol は NaN、桁あふれを直しても追跡と 3 桁以上合わない)"""
    from . import compat
    srot = st.srot
    srot1 = srot[:, 0:3].copy()
    sps, smu, sdamp = spnorm(srot, st.gintd)
    params = dict(emit=(em['emittance_x'], em['emittance_y'], em['emittance_z']),
                  damp=em['damping_rate'], tunes=em['tunes'], smu=smu)
    equpol, spm = srequpol(srot, sps, params, em['beamn'], sdamp, em['normal_r'])
    frev = ph.omega0 / M_2PI
    with np.errstate(all='ignore'):
        taup = 1.0 / sdamp / frev if sdamp != 0.0 else float('inf')
    res = dict(srot=srot.copy(), srot1=srot1, sps=sps, smu=smu, sdamp=sdamp, nup=smu / M_2PI,
               equpol=equpol, spm=spm, taup=taup, gintd=st.gintd.copy())
    n_ax = sps[:, 0]
    sdamp_rel = float(n_ax @ getattr(st, 'gintd2', np.zeros((3, 3))) @ n_ax)       # Σ r_i (b_i·n)²: ST の緩和率 (符号に依らない)
    res['sdamp_rel'] = sdamp_rel
    if compat.v2() and els is not None and cod is not None:
        res['equpol_k64'] = equpol                      # K64 の spdepol の式 (桁あふれだけ直したもの)。参考値
        try:
            if not np.all(np.isfinite(cod)):
                raise ValueError('closed orbit not found (unstable optics)')        # D1: 要素ごとの tturne を 2000 回回さない
            sl = slim_depolarization(els, ph, cod, trf0=trf0, rfsw=rfsw, radcod=radcod,
                                     beam_matrix=em.get('beam_matrix') if (SIDEBAND and rfsw) else None, sideband=SIDEBAND and rfsw)
            rate = sl['rate']
            # v2 (A2): dP/dt = PST·Σ r_i (b_i·n) − P·(Σ r_i (b_i·n)² + rate)。K64 の sdamp = n·gintd は符号つきの Σ r_i (b_i·n) で、
            # 逆向きの磁場 (挿入光源の交互の極) があると緩和率に使えない (P = PST (一定)、τ が実際より長く出る)。
            # 緩和率は Σ r_i (b_i·n)² (向きの符号に依らない)。一方向の曲げだけのリングでは sdamp と同じ (正の曲げだけ: 完全に一致)
            srel = sdamp_rel if sdamp_rel > 0.0 else sdamp
            gsum = float(getattr(st, 'gsum', 0.0))
            if ST_EXACT and gsum > 0.0:
                # v2 (P1): Sokolov-Ternov の厳密な緩和率。スピンの向き n・磁場の向き b の符号つき・2 乗ではなく、放射 (スピン反転) の率 ∝ ρ⁻³ (1 − 2/9 (n·v)²)
                # (v は速度の向き)。n が b と垂直 (垂直の偏向のある平面リング、水平磁場の挿入光源) でも反転は起きる。n が b と平行 (平面リング) では Σ r (b·n)² と同じ
                res['sdamp_rel_b2'] = sdamp_rel
                srel = gsum - (2.0 / 9.0) * float(n_ax @ getattr(st, 'gintd2v', np.zeros((3, 3))) @ n_ax)
                res['sdamp_rel'] = srel
            kin = 0.0
            if KINETIC and sl.get('kin'):
                # v2 (P1): 運動論的偏極 (DK の式の分子 b·(n − ∂n/∂δ) の第 2 項)。slim の n0 と K64 の n (sps) の向きをそろえる
                kin = float(np.sign(float(sl['n0'] @ n_ax)) * sl['kin'])
            res['kinetic'] = kin
            with np.errstate(all='ignore'):
                pol = PST * (sdamp - kin) / (srel + rate) if (srel + rate) != 0.0 else float('nan')
            if srel > 0.0:
                with np.errstate(all='ignore'):
                    res['taup'] = 1.0 / srel / frev
            npp = (np.eye(3) - np.outer(sl['n0'], sl['n0'])) @ sl['N0']
            res.update(depol_rate_linear=sl.get('rate_linear', rate), depol_sideband=('spread_sb' in sl), spin_spread_sb=sl.get('spread_sb'),
                       depol_chromatic=bool(sl.get('chromatic', False)), mod_spin=sl.get('mod_spin'), mod_chrom=sl.get('mod_chrom'))
            if 'sideband_error' in sl:
                res['sideband_error'] = sl['sideband_error']
            res.update(equpol=np.array([pol, 0.0, 0.0]), depol_rate=rate, spin_orbit=sl['N0'],
                       spin_spread=0.5 * float(np.trace(npp @ em['beam_matrix'] @ npp.T)),
                       taudep=(1.0 / rate / frev if rate > 0.0 else float('inf')))
        except Exception as exc:            # 求まらないときは NaN (K64 と同じ表示) にして理由を残す
            res.update(equpol=np.array([float('nan'), 0.0, 0.0]), depol_error=str(exc))
    return res
