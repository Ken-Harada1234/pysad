# -*- coding: utf-8 -*-
"""v2next B3 項目 1: 側帯つきの脱偏極率 (spin.sideband_depolarization) を、線形の軌道 + 厳密なスピン回転の Monte Carlo で検証する。
モデル: 要素ごとに  u <- T_i u + (放射の跳び、共分散 B_i)、  s <- exp(-[D_i u]x) R_i s  (T_i, B_i, R_i, D_i は Emittance[] の要素ごとの写像。
D_i は K64 の回転の向き s' = (I − [D u]x) R s の余分の回転を厳密な回転にしたもの)。1 − <S·n0> の増加の傾きが 1 周あたりの脱偏極率。
追跡 (track6) と違い軌道は線形 (軌道の非線形・スピンの高次を含まない) なので、式の導出 (位相変調) だけを調べられる。
使い方: py -3.12 example/v2_next/B3_spin_mc.py E[GeV] V[MV] 粒子数 周数 種 出力.json [fscale]   (例: 4.2 12 3000 30000 11 out.json)
  fscale = 放射の跳びと初期分布の分散にかける倍率 (1 が実際。0.01 なら変調が小さい線形の領域で、傾きは fscale 倍)。
リング: 12 セル FODO (四極の縦ずれ ±3 mm、RF 調和数 64。test_v2_A_spin.py の _lat と同じ)。1 粒子 1 周 約 10 ms (粒子数 3000)。
"""
import io, json, math, sys, time
import numpy as np
from pysad import compat, emit6d, ffscmd, spin
compat.set_mode('v2')


def lat(E, volt, dy=3e-3):
    return f'''MOMENTUM = {E}GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY={dy}) QDE = (L=0.3 K1=-0.6 DY=-{dy});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={volt}e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''


def build(ffs):
    p = ffs.momentum()
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
    els = ffs.lat.elements
    cod = np.array(e['orbit_entrance'])
    trf0 = e.get('trf0', 0.0)
    ph = emit6d._Phys(p, emit6d._emit_circ(els))
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    secs = []
    c = cod.copy()
    try:
        for el in els:
            st = emit6d.tturne([el], c, ph, trf0=trf0, rfsw=True, radcod=True)
            secs.append((st.T + st.D, st.B.copy(), st.srot[:, :3].copy(), st.srot[:, 3:9].copy()))
            c = st.cod.copy()
    finally:
        emit6d.POL_MODE = old
    return e, els, cod, trf0, ph, secs


def rot(s, w):
    th = np.linalg.norm(w, axis=1)[:, None] + 1e-300
    k = w / th
    c, sn = np.cos(th), np.sin(th)
    return s * c + np.cross(k, s) * sn + k * (np.sum(k * s, axis=1)[:, None]) * (1 - c)


if __name__ == '__main__':
    E, V, N, turns, seed, out = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    FS = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    ffs, _ = ffscmd.run_sad_string(lat(E, V), out=io.StringIO())
    e, els, cod, trf0, ph, secs = build(ffs)
    sl = spin.slim_depolarization(els, ph, cod, trf0=trf0, rfsw=True, radcod=True, beam_matrix=e['beam_matrix'], sideband=True)
    n0, N0 = sl['n0'], sl['N0']
    steps, Tacc = [], np.eye(6)
    for T, B, R, D in secs:                      # 何もしない要素 (ドリフトなど) は次の要素の前に畳む (D u は要素の入口の u を使う)
        if np.abs(D).max() == 0 and np.abs(R - np.eye(3)).max() == 0 and np.abs(B).max() == 0:
            Tacc = T @ Tacc
        else:
            steps.append((Tacc, T, B, R, D))
            Tacc = np.eye(6)
    tail = Tacc
    Ls = []
    for _Tp, T, B, R, D in steps:
        if np.abs(B).max() > 0:
            w, v = np.linalg.eigh(B)
            Ls.append(v * np.sqrt(np.clip(w, 0, None) * FS))
        else:
            Ls.append(None)
    rng = np.random.default_rng(seed)
    u = (np.linalg.cholesky(e['beam_matrix'] * FS) @ rng.standard_normal((6, N))).T
    S = np.tile(n0, (N, 1)) + u @ N0.T
    S /= np.linalg.norm(S, axis=1)[:, None]
    res = dict(E=E, V=V, N=N, seed=seed, fscale=FS, rate_linear=sl['rate_linear'], rate_sideband=sl['rate'], D=[], turns=[])
    t0 = time.time()
    every = 2000
    for t in range(turns + 1):
        if t % every == 0:
            d = 1 - S @ n0
            res['turns'].append(t)
            res['D'].append([float(d.mean()), float(d.std() / math.sqrt(N))])
            json.dump(res, open(out, 'w'))
            print(t, res['D'][-1], time.time() - t0, flush=True)
        for (Tp, T, B, R, D), L in zip(steps, Ls):
            u = u @ Tp.T
            S = S @ R.T
            S = rot(S, -(u @ D.T))
            u = u @ T.T
            if L is not None:
                u = u + rng.standard_normal((N, 6)) @ L.T
        u = u @ tail.T
