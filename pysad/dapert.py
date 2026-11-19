"""
dapert.py — DynamicApertureSurvey (SAD の Packages/dapert.n → tfdapert.f → trackd.f) の直訳

trackd の流れ (NPARA=1):
  - 軸 Axes ("ZX" 既定 → ivar1=3 (Z), ivar2=1 (X), ivar3=2 (Y))。range[[ivar1]] は ivar1 の値のリスト (a1i)、
    range[[ivar2]] = {a2min,a2max}、range[[ivar3]] = {a3min,a3max}。ivar2/ivar3 は n2p=51 点を同じ j で動かす。
  - 粒子の初期値 (tinip1): 規格化振幅 (x, px, y, py, z, g) に emx = sqrt(|EMITX|+|EMITY|), emz = sqrt(|EMITZ|) (RFSW) か
    |SIGE| を掛け、r (正規モード→物理) を掛けて codin (emit の閉軌道) を足す。
    ivar1 が Z なら z=0, g=a1。ivar2/3 は x = a cosφ, px = −a sinφ (PhaseX 等)。
  - 粒子は SAD と同じ順 (i = 1..n1p、j = 51 → 1) でスロットに置き、全部まとめて追跡する (npmax ≥ n1p·51 のとき)。
    tturn1 の tapert (最初の要素、以後 15 要素ごと、周末) で |x|,|y| < LOSSAMPL, |z| ≤ LOSSDZ を満たさない粒子は
    末尾の生きている粒子と入れ替えて外す (track6n.k_apert = tapert)。この順序で FLUC の乱数も SAD と同じになる。
    ntloss = 損失までに完了した周数 (0 なら 1 周目で損失)。
  - 行ごとに nsc = 最初に損失した j − 1、Score = Σ nsc。表示は f8.2, i3, 1x, 51 文字 (rad62(ntloss/muls)、0 は ' ' か '.')。
  - ResultOfDAPERT = {{{a2min,a2max},{a3min,a3max},{a1min,a1max}}, {{a1i, nsc, {ntloss×51}}, ...}}
    DynamicApertureSurvey は {score, ReplacePart[ResultOfDAPERT, range, 1]} を返す。
未対応: Damping->True (tpdamp)、DAPWIDTH の窓 (npmax < n1p·51 のときの追跡順)。
"""
from __future__ import annotations

import math
import os
import time

import numpy as np

from . import track6 as _t6
from . import track6n as _tn
from .emit6d import _Phys

N2P = 51
MAXTURN = 2**29
_C62 = '0123456789AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz'
_LABEL = {1: 'X', 2: 'Y', 3: 'Z'}


def rad62(i: int) -> str:
    return _C62[i] if 0 <= i <= 61 else '*'


def rad62a(m: int, j: int) -> str:
    if m == 0:
        return '.' if (j % 10 == 1 and j != 1) else ' '
    return rad62(m)


def _muls(nturn: int) -> int:
    muld, mule = 1, 1
    while muld * mule * 62 < nturn:
        if muld == 1:
            muld = 2
        elif muld == 2:
            muld = 5
        else:
            muld = 1
            mule *= 10
    return muld * mule


def dynamic_aperture(elements, range3: list, r1: list, nturn: int, phi, ivar1: int, ivar2: int,
                     momentum_eV: float, emit: dict, emx: float, emz: float, trf0: float, rfsw: bool, rad: bool,
                     fluc: bool = True, alost: float = 1.0, zlost: float = 100.0, out=None, progress: int = 0,
                     fshift: float = 0.0, dapwidth: float = 51.0, np0=None):
    """trackd 相当。range3: [[a2min,a2max],[a3min,a3max]] (ivar2, ivar3 の順)、r1: ivar1 の値。
    戻り: (score, result) result = ResultOfDAPERT"""
    ivar3 = 6 - ivar1 - ivar2
    a2min, a2max = range3[0]
    a3min, a3max = range3[1]
    a1i = list(r1[:256])
    n1p = len(a1i)
    a1min, a1max = min(a1i), max(a1i)
    a2step = (a2max - a2min) / (N2P - 1)
    a3step = (a3max - a3min) / (N2P - 1)
    muls = _muls(nturn)
    w = out if out is not None else (lambda s: None)
    L1, L2, L3 = _LABEL[ivar1], _LABEL[ivar2], _LABEL[ivar3]
    w(f' Range    {L2}min:{a2min:8.3f} {L2}max:{a2max:8.3f}\n')
    w(f'         ({L3}min:{a3min:8.3f} {L3}max:{a3max:8.3f})\n')
    w(f'          {L1}min:{a1min:8.3f} {L1}max:{a1max:8.3f}\n')
    w(f' Display: {muls:5d} turns/character\n')
    w(f'      N{L1}     0----|----1----|----2----|----3----|----4----|----5\n')
    cx, sx = math.cos(phi[0]), math.sin(phi[0])
    cy, sy = math.cos(phi[1]), math.sin(phi[1])
    cz, sz = math.cos(phi[2]), math.sin(phi[2])
    rmat = np.array(emit['normal_r'], dtype=float)
    codin = np.array(emit['closed_orbit'], dtype=float)
    circ = sum(el.L for el in elements)
    ph = _Phys(momentum_eV, circ)
    ph.dvfs = fshift if rfsw else 0.0          # tsetdvfs (trackd.f:86)
    if not progress:
        progress = int(os.environ.get('PYSAD_PROGRESS', '0') or 0)
    # trackd のスロット管理の直訳 (npara = 1): 粒子スロット k (ID) は配列位置 kptbl(k,1) = pos[k] にいて、kid[位置] = ID。
    # 点 (i, j) は、行 i の窓 (損失が見つかった最小の j = nxm(i) から下へ ncons 点) の中だけ追跡する (DAPWIDTH)。
    # 同時に追跡するスロットは npmax = max(1, min(n1p·ncons, np0))。損失・完了でスロットが空くと配置ループへ戻り、
    # 空きがあれば nw = 16 周ごとに戻る (trackd.f)。窓の外 (未起動) の点の ntloss は 2^29 のまま (表示は '*')
    ncons = max(1, int(math.floor(float(dapwidth) + 0.5)))
    if np0 is None:
        np0 = max(2000, n1p, n1p * min(51, ncons))
    npmax = max(1, min(n1p * ncons, int(np0)))

    def _coord(i, j):           # j = 1..51 (Fortran)。ivar1 の値 a1i[i]、ivar2・ivar3 は同じ j で動く
        amp = np.zeros(6)

        def setvar(iv, a, is1=False):
            if iv == 1:
                amp[0], amp[1] = (a, 0.0) if is1 else (a * cx, -a * sx)
            elif iv == 2:
                amp[2], amp[3] = (a, 0.0) if is1 else (a * cy, -a * sy)
            else:
                amp[4], amp[5] = (0.0, a) if is1 else (a * cz, -a * sz)
        setvar(ivar1, a1i[i], True)
        setvar(ivar2, a2step * (j - 1) + a2min)
        setvar(ivar3, a3step * (j - 1) + a3min)
        xa = amp * np.array([emx, emx, emx, emx, emz, emz])
        return rmat @ xa + codin

    prog = _tn.Prog(elements, ph)
    S = _tn.make_state(np.tile(codin, (npmax, 1)), ph)
    kid = np.arange(npmax, dtype=np.int64)          # 配列位置 → ID (kptbl(:,2))
    pos = np.arange(npmax, dtype=np.int64)          # ID → 配列位置 (kptbl(:,1))
    lost_scratch = np.full(npmax, MAXTURN, dtype=np.int64)
    kzx_i = np.full(npmax, -1, dtype=np.int64)       # スロットの (i, j)。i < 0 は空き
    kzx_j = np.zeros(npmax, dtype=np.int64)
    mturn = np.zeros(npmax, dtype=np.int64)
    ntloss = np.full((n1p, N2P + 1), MAXTURN, dtype=np.int64)       # 列 1..51 を使う
    nxm = [N2P + 1] * n1p
    clean = [False] * n1p                            # 走査して起動できる点が無かった行 (ntloss が変わるまで再走査しない。結果は同じ)
    rs = _tn.make_rs(_t6.RNG.seed)          # [seed, v2 の印] (カーネルが rs[1] を読む)
    fl = bool(rad and fluc)
    t0 = time.time()
    nw = 16
    nturns_done = 0
    with np.errstate(all='ignore'):
        while True:                                  # loop_1
            np1 = npmax
            iw = nw
            remain = True
            for k in np.nonzero(kzx_i < 0)[0]:       # LOOP_K (空きスロットを ID 順に)
                ip = int(pos[k])
                if ip >= np1:
                    continue
                if remain:
                    launched = False
                    for i in range(n1p):
                        if clean[i]:
                            continue
                        for j in range(nxm[i] - 1, 0, -1):
                            if j < nxm[i] - ncons:
                                break
                            if ntloss[i, j] < nturn:
                                nxm[i] = j
                            elif ntloss[i, j] == MAXTURN:
                                ntloss[i, j] = nturn
                                kzx_i[k] = i
                                kzx_j[k] = j
                                mturn[k] = 0
                                S[ip] = _tn.make_state(_coord(i, j)[None, :], ph)[0]
                                launched = True
                                break
                        if launched:
                            break
                        clean[i] = True
                    if launched:
                        continue
                remain = False
                jid = int(kid[np1 - 1])              # 位置 ip と np1 の粒子を入れ替え、np1 の座標を ip へ
                pos[jid] = ip
                pos[k] = np1 - 1
                kid[ip] = jid
                kid[np1 - 1] = k
                S[ip] = S[np1 - 1]
                np1 -= 1
                if np1 <= 0:
                    break
            n = np1
            while n > 0:
                n = _tn.turn(prog, S, n, ph, trf0, rfsw, rad, fl, rs, apert=(kid, lost_scratch, nturns_done, alost, zlost))
                nturns_done += 1
                pos[kid] = np.arange(npmax)
                if not rfsw and n > 0:
                    S[:n, 4] = 0.0
                act = np.nonzero(kzx_i >= 0)[0]
                alive = pos[act] < n
                mturn[act[alive]] += 1
                done = act[alive][mturn[act[alive]] >= nturn]
                lost = act[~alive]
                for k in done:
                    ntloss[kzx_i[k], kzx_j[k]] = nturn
                    clean[kzx_i[k]] = False
                    kzx_i[k] = -1
                for k in lost:
                    ntloss[kzx_i[k], kzx_j[k]] = mturn[k]
                    clean[kzx_i[k]] = False
                    kzx_i[k] = -1
                if progress and nturns_done % progress == 0:
                    print(f'    [dapert] pass {nturns_done} alive {n}/{npmax} ({time.time() - t0:.1f} s)', flush=True)
                if len(done) or len(lost):
                    break                            # cycle loop_1
                if np1 < npmax:
                    iw -= 1
                    if iw <= 0:
                        break
            else:
                break                                # np = 0: exit loop_1
    _t6.RNG.seed = int(rs[0])
    ntloss = ntloss[:, 1:]                   # 列 j = 1..51
    rows = []
    nscore = 0
    for i in range(n1p):
        nsc = N2P
        for k in range(N2P):
            if ntloss[i, k] < nturn:
                nsc = k
                break
        nscore += nsc
        w(f' {a1i[i]:8.2f}{nsc:3d} ' + ''.join(rad62a(int(ntloss[i, jj] // muls), jj + 1) for jj in range(N2P)) + '\n')
        rows.append([float(a1i[i]), float(nsc), [float(v) for v in ntloss[i]]])
    w(f'      N{L1}     0----|----1----|----2----|----3----|----4----|----5\n')
    w(f'    Score: {nscore:5d}\n')
    result = [[[a2min, a2max], [a3min, a3max], [a1min, a1max]], rows]
    return float(nscore), result
