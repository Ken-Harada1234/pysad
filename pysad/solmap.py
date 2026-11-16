"""
SOL 領域 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで) の COD 写像と行列 (CALC / INS の Twiss)。K64 の tsole1 (tsole.f) 相当。

領域内の要素の写像は追跡 (track6n_ext の K_SOL_FRIN / K_SOL_DRIFT / K_SOL_QUAD、K64 の tsol.f) と同じ座標 (ソレノイド内の
正準運動量) なので、その numba カーネルを 1 粒子に使い、行列は COD での数値微分にする (tsolque の解析行列ではない)。
領域の外 (BZ=0) の要素は通常の写像。領域内の要素の Twiss は領域内の座標での値。

対応: DRIFT、QUAD (K1≠0、ずれ・F1 なし)、L=0 の要素、SOL (ずれ・回転・F1 なし)。追跡と同じ制限。
未対応: Emittance[]・CALC6D・EMIT (tradke の bz 付き。emit6d 側で拒否する)、BEAMBEAM (K64 が落ちる)、
        MULT・BEND などの領域内の要素 (K64 の tsole1 は対応するが、追跡側 compile が未対応)。
"""
from __future__ import annotations

import numpy as np

from . import track6n as _n
from . import track6n_ext as _x
from . import emit6d as _e


class _Collect:
    def __init__(self):
        self.items = []

    def _add(self, k, p, name, cnt, sm):
        self.items.append((k, p))


def annotate(els, momentum_eV):
    """SOL 領域内の各要素に el._solmap = (kind, params) を付ける (領域外は消す)。領域がなければ False"""
    found = False
    for el in els:
        if hasattr(el, '_solmap'):
            del el._solmap
            el.__dict__.pop('transfer_matrix', None)
            el.__dict__.pop('constant_kick', None)
    ph = _e._Phys(momentum_eV, sum(el.L for el in els))
    i = 0
    while i < len(els):
        el = els[i]
        if type(el).__name__ == 'Solenoid':
            bnd = float(getattr(el, 'BOUND', 0.0) or 0.0)
            bz = float(getattr(el, 'BZ', 0.0) or 0.0)
            if bnd == 0.0:
                if bz != 0.0:
                    raise NotImplementedError('SOL: BOUND=0 で BZ≠0 の SOL (領域の外) の Twiss は未対応')
                i += 1
                continue
            col = _Collect()
            ke = _x.compile_sol_region(col, els, i, ph)
            for j, kp in enumerate(col.items):
                e = els[i + j]
                e._solmap = (kp[0], kp[1], momentum_eV)
                # 線形の経路 (transfer_matrix を直接呼ぶ所) でも同じ写像を使う: 設計軌道 (δ だけ) での行列
                e.transfer_matrix = (lambda delta=0.0, e=e: step(e, np.array([0, 0, 0, 0, 0, delta], dtype=float))[0])
                e.constant_kick = (lambda delta=0.0: np.zeros(6))
            found = True
            i = ke + 1
            continue
        i += 1
    return found


def cod_map(kind, pr, momentum_eV=None):
    """(kind, params) → c (6) → c (6) の写像"""
    ph = _e._Phys(momentum_eV or 1.0e9, 1.0)
    rs = np.array([1], dtype=np.int64)

    def f(c):
        S = _n.make_state(np.array(c, dtype=float).reshape(1, 6), ph)
        _x.step_kind(kind, pr, S, 1, ph, 0.0, False, False, False, rs)
        return _n.state_out(S, 1)[0]          # 物理座標 (px は pr 倍して戻す)
    return f


def step(el, cod):
    """_el_cod_step 用: (M_el, 出口 COD)"""
    from . import lattice as _lat
    kind, pr, momentum_eV = el._solmap
    if kind == _n.K_NONE:
        return np.eye(6), np.array(cod, dtype=float)
    f = cod_map(kind, pr, momentum_eV)
    c = np.array(cod, dtype=float)
    return _lat._num_jacobian(f, c), f(c)
