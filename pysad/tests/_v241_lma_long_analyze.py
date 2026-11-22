# -*- coding: utf-8 -*-
"""v2.4.1: 長時間追跡 (gs・pysad の出力 `S n i x px y py z dp flag` の行) の解析。
σδ・σz・εx・εy を n ごとに出し、Emittance[] の値 (EMIT の行) から作る予想と比べる。
  FLUC: 閉軌道から出発し ε(n) = ε_eq (1 − exp(−2 α n)) (α = DampingRate の大きさ、平衡の値は Emittances・MomentumSpread。σδ² も同じ形)。
  NOFLUC: 初期の偏りが ε(n) = ε(0) exp(−2 α n) (RF のあるリングでは δ は同期振動をする。δ・z の包絡の減衰 = 振幅² の平均を取る)。
εx = ½<γx̃² + 2αx̃p̃ + βp̃²> (x̃ = x − ηδ)。誤差は N 粒子の標本の統計 (J は指数分布で εの相対誤差 1/√N、σδ² は √(2/N))。
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SP = os.path.join(HERE, 'scripts')
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)


def load(path):
    """戻り: dict(n: (nt,), 座標 (nt, N, 6), flag (nt, N), emit: {...})"""
    txt = open(path, encoding='utf-8').read()
    rows = []
    for ch in re.split(r'(?m)^S ', txt)[1:]:      # gs は 1 行が長いと折り返す (6 個目で改行) ので、チャンク全体の数値を集める
        v = ch.split()[:9]
        rows.append([float(x) for x in v])
    a = np.array(rows)
    ns = np.unique(a[:, 0])
    N = int(a[:, 1].max())
    A = a.reshape(len(ns), N, -1)
    em = {}
    mm = re.search(r'(?m)^EMIT ', txt)
    if mm:
        i = mm.start()
        m2 = re.search(r'(?m)^S ', txt[i:])
        blk = re.sub(r'\s+', '', txt[i:i + m2.start()] if m2 else txt[i:])
        m = re.search(r'DampingRate->\{([^}]*)\}', blk)
        em['alpha'] = [abs(float(v)) for v in m.group(1).split(',')]
        m = re.search(r'Emittances->\{([^}]*)\}', blk)
        em['eps'] = [float(v) for v in m.group(1).split(',')]
        m = re.search(r'MomentumSpread->([-0-9.eE+]*)', blk)
        em['sigd'] = float(m.group(1))
        m = re.search(r'BunchLength->([-0-9.eE+]*)', blk)
        em['sigz'] = float(m.group(1))
    return dict(n=ns, X=A[:, :, 2:8], flag=A[:, :, 8], N=N, emit=em)


def twiss(ring):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    import _v241_lma_tools as T
    compat.set_mode('v2')
    st, _ = TrackSetup.from_source(T.read({'astrid2': T.ASTRID2, 'kekls': T.KEKLS12}[ring]), 'RING', rad=True, rfsw=True, alost=0.1)
    return st.twiss_at(1), st.closed_orbit([1])[0]


def stats(d, tw, cod):
    """各 n で alive の粒子だけ。戻り: dict(n, nalive, sd, sz, ex, ey)"""
    X = d['X'] - cod[None, None, :]
    out = []
    for k in range(len(d['n'])):
        al = d['flag'][k] > 0.5
        x, px, y, py, z, dp = (X[k][al][:, j] for j in range(6))
        xb = x - tw['etx'] * dp
        pb = px - tw['etpx'] * dp
        Jx = (1 + tw['alfx'] ** 2) / tw['betx'] * xb ** 2 + 2 * tw['alfx'] * xb * pb + tw['betx'] * pb ** 2
        Jy = (1 + tw['alfy'] ** 2) / tw['bety'] * y ** 2 + 2 * tw['alfy'] * y * py + tw['bety'] * py ** 2
        out.append((al.sum(), np.sqrt(np.mean(dp ** 2)), np.sqrt(np.mean(z ** 2)), 0.5 * Jx.mean(), 0.5 * Jy.mean()))
    o = np.array(out)
    return dict(n=d['n'], nalive=o[:, 0], sd=o[:, 1], sz=o[:, 2], ex=o[:, 3], ey=o[:, 4])
