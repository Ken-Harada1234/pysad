# -*- coding: utf-8 -*-
"""v2.3.0 R5: 追跡の 1 周の多項式の写像から、数値の正規形で ADTS を取り出す (実装の範囲: J の 1 次 = ∂ν/∂J まで。3 次の写像 = f3・f4)。

写像 (正規化座標、リング先頭) M = R e^{−:h3 + h4:}。追跡の 1 周の写像の 2 次・3 次のテイラー展開 (rdttrack.ring_map3) から f3 の係数 h3 (rdttrack.rdt_from_map)・
f4 の係数 h4 (rdttrack.rdt4_from_map、Dragt-Finn) を数値で取り出し、
  ∂νx/∂Jx = 4 h22000/π + (正規形の補正)、∂νx/∂Jy = 2 h11110/π + …、∂νy/∂Jy = 4 h00220/π + …、
補正は 3 次の項 h3 の非共鳴の摂動 `rdt2nd.nf_correction` (解析式の側と同じ式、docs/49 §10.4)。偏向磁石・端部・8 極・運動学的な項・積分器の誤差は追跡の写像に自動で入る。
解析式 (`resdrive.adts`、`rdt2nd.adts_from_rdt4`) との違いは、h3・h4 を追跡から取り出す点だけ (式の枠組みの検証にも使う)。
未実装 (設計は docs/49 §10.6): J² の項 (∂²ν/∂J²。5 次の写像から f5・f6、2 次の正規形の反復)、大振幅の写像の多項式 (7 次以上) と共鳴の島の正規形。
"""
import math
import numpy as np
from . import rdttrack as T
from . import rdt2nd as Q


def adts_from_map(rm):
    """ring_map3 の結果 rm から数値の正規形の ADTS {'xx','xy','yy'} [1/(m rad)] と、使った h3 (MONOS3 の順 20 成分)・h4 の対角 (dict)。
    rm.R の位相進み (チューン) は R の固有値から求める"""
    h3d = T.rdt_from_map(rm)
    h3 = np.array([h3d['h%d%d%d%d0' % q] for q in Q.MONOS3])
    h4d, info = T.rdt4_from_map(rm)
    R = rm.R
    # 正規化座標の R は回転: R[0,0] = cos μ、R[0,1] = sin μ (x_n, p_n) の向きは W_NORM (ζ = x_n − i p_n) に合わせて作る
    mux = math.atan2(R[0, 1], R[0, 0]) % (2 * math.pi)
    muy = math.atan2(R[2, 3], R[2, 2]) % (2 * math.pi)
    qx_, qy_ = mux / (2 * math.pi), muy / (2 * math.pi)
    corr = Q.nf_correction(h3, qx_, qy_)
    diag = {'xx': 4.0 * h4d['h22000'].real / math.pi, 'xy': 2.0 * h4d['h11110'].real / math.pi, 'yy': 4.0 * h4d['h00220'].real / math.pi}
    return {k: diag[k] + corr[k] for k in diag}, dict(h3=h3, diag=diag, corr=corr, tunes=(qx_, qy_), symp=info['symp'])


def adts_at_delta(lat, momentum_eV, delta, amp3=1e-4, amp=1e-5):
    """v2.3.0 R6: δ のエネルギーの ADTS を、追跡の光学 (rdttrack.tracked_rows) と δ の閉軌道のまわりの 1 周の写像 (ring_map3(delta=δ)) の数値の正規形で求める
    (解析式を使わない。直接の追跡 + NAFF と 1e-4 以内: ASTRID2・FODO の δ = ±0.002)。amp3: 3 次のテイラー展開の振幅 (大きなリングは 1e-5)。戻り ({'xx','xy','yy'}, info)"""
    base = lat.twiss_table(delta=delta)
    rows = T.tracked_rows(lat.elements, momentum_eV, delta, base)
    rm = T.ring_map3(lat.elements, momentum_eV, rows[0], delta=delta, amp=amp, amp3=amp3)
    return adts_from_map(rm)


def adts_delta_derivative(lat, momentum_eV=3e9, delta0=2e-3, method='analytic', **kw):
    """v2.3.0 R6: ∂(∂ν/∂J)/∂δ (δ = 0) の中心差分 (δ = ±delta0)。method: 'analytic' = `resdrive.adts(delta=±delta0)` (解析式。偏向磁石・四極・ドリフトの δ の閉軌道のまわりの項。光学は既定で CALC、`optics='tracked'` を渡すと追跡の光学)、
    'nform' = 追跡の写像の数値の正規形 (`adts_at_delta`)。戻り {'xx','xy','yy'} [1/(m rad)]。差分の誤差は O(delta0²) (δ の 3 次の項)"""
    from . import resdrive as R
    out = {}
    for sg in (+1, -1):
        d = sg * delta0
        if method == 'nform':
            out[sg] = adts_at_delta(lat, momentum_eV, d, **kw)[0]
        else:
            out[sg] = {k: float(v) for k, v in R.adts(lat, delta=d, momentum_eV=momentum_eV, **kw).total.items()}
    return {k: (out[1][k] - out[-1][k]) / (2.0 * delta0) for k in ('xx', 'xy', 'yy')}
