# -*- coding: utf-8 -*-
"""バンチ内の δ の分布 (touschek_bunch.py) の独立な検証: 実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学、touschek_formula_mc.py と同じ)。
同じ対の標本で 2 通りの損失の判定を数える:
  (kick) δ_new − δ_old が a+ を超える、または −a− を下回る (touschek_piwinski の基準 = δ_c を 0 とする)
  (abs)  δ_new (実験室系の絶対値) が a+ を超える、または −a− を下回る (元の δ を含める。元の δ の絶対値が 0.75 a を超える対は除く = 受容の近くの粒子の小角散乱を除く)
増え = abs / kick − 1。理論 (touschek_bunch.piwinski_slice_rate_bunch): 左の列 = 分散 σδ² (var_factor=1)、右の列 = 対の重心の δ_c の分散 σδ²/2 (var_factor=0.5、実装の既定)。診断: qs_scale・cm_scale で相対運動量の δ・重心の δ を消した増えの内訳 (A3_decompose)。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_bunch_mc.py astrid2|kekls [a ...]   (出力は標準出力。保存は touschek_bunch_mc_<名前>.out.txt)
"""
import math
import os
import sys

import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from scipy.stats import qmc
import touschek_formula_mc as M


def mc_both(S, al, circ, N, p0, a, npairs=8192, nth=40, nph=32, seed=1, chunk=1024, trunc=0.75, qs_scale=1.0, cm_scale=1.0):
    """(rate_kick, rate_abs, 標本の割合, ...) を同じ標本で返す (ap = am = a)。qs_scale: 相対運動量の δ 成分の倍率 (診断用: 0 で相対の δ を消す)、cm_scale: 対の重心の δ の倍率"""
    h0 = math.sqrt(1 + p0 ** 2)
    xx, pp, P = M.lab_cond(S)
    vol = (4 * np.pi) ** 1.5 * math.sqrt(np.linalg.det(xx))
    Lq = np.linalg.cholesky(2 * P)
    Lc = np.linalg.cholesky(pp / 2)
    z = qmc.MultivariateNormalQMC(mean=np.zeros(6), cov=np.eye(6), seed=seed).random(npairs)
    qrel = z[:, :3] @ Lq.T
    cmom = z[:, 3:] @ Lc.T
    qrel[:, 2] *= qs_scale
    cmom[:, 2] *= cm_scale
    p1 = cmom + qrel / 2
    p2 = cmom - qrel / 2
    f1 = M.mom4(p1[:, 0], p1[:, 1], p1[:, 2], p0)
    f2 = M.mom4(p2[:, 0], p2[:, 1], p2[:, 2], p0)
    tot = f1 + f2
    mtot = np.sqrt(tot[:, 0] ** 2 - np.sum(tot[:, 1:] ** 2, axis=1))
    u = tot / mtot[:, None]
    q1 = M.boost_to(u, f1)
    pm = np.linalg.norm(q1[:, 1:], axis=1)
    k = 2 * pm
    uh = q1[:, 1:] / pm[:, None]
    dot = f1[:, 0] * f2[:, 0] - np.sum(f1[:, 1:] * f2[:, 1:], axis=1)
    vM = np.sqrt(np.maximum(dot ** 2 - 1, 0)) / (f1[:, 0] * f2[:, 0])
    d1 = np.linalg.norm(f1[:, 1:], axis=1) / p0 - 1
    d2 = np.linalg.norm(f2[:, 1:], axis=1) / p0 - 1
    keep = (np.abs(d1) < trunc * a) & (np.abs(d2) < trunc * a)
    ref = np.where(np.abs(uh[:, 0:1]) < 0.9, np.array([[1, 0, 0]]), np.array([[0, 1, 0]]))
    e1 = np.cross(uh, ref)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(uh, e1)
    gx, gw = np.polynomial.legendre.leggauss(nth)
    phis = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    out_k = np.zeros(npairs)
    out_a = np.zeros(npairs)
    out_c = np.zeros(npairs)
    dcen = 0.5 * (d1 + d2)
    for c0 in range(0, npairs, chunk):
        sl = slice(c0, min(npairs, c0 + chunk))
        n = sl.stop - sl.start
        kk = k[sl]
        cmin = (1 - trunc) * a            # abs の判定の最小の kick (δ_old が 0.75a までなので a − 0.75a)
        sinhalf = np.minimum(0.5 * cmin * (p0 / h0) / (2 * pm[sl]), 1.0)
        thlo = 2 * np.arcsin(np.minimum(sinhalf, 1.0))
        thhi = np.full(n, np.pi / 2)
        ok = thlo < thhi
        lo = np.log(np.maximum(thlo, 1e-9))
        hi = np.log(thhi)
        lt = 0.5 * (hi - lo)[:, None] * gx[None, :] + 0.5 * (hi + lo)[:, None]
        th = np.exp(lt)
        wt = 0.5 * (hi - lo)[:, None] * gw[None, :] * th * np.sin(th) * M.Kth(th) * (2 * np.pi / nph)
        st = np.sin(th)[:, :, None, None]
        ct = np.cos(th)[:, :, None, None]
        nd = ct * uh[sl][:, None, None, :] + st * (np.cos(phis)[None, None, :, None] * e1[sl][:, None, None, :] + np.sin(phis)[None, None, :, None] * e2[sl][:, None, None, :])
        q1n = np.empty(nd.shape[:-1] + (4,))
        q1n[..., 0] = q1[sl][:, 0][:, None, None]
        q1n[..., 1:] = pm[sl][:, None, None, None] * nd
        q2n = q1n.copy()
        q2n[..., 1:] *= -1
        un = np.broadcast_to(u[sl][:, None, None, :], q1n.shape)
        l1 = M.boost_from(un, q1n)
        l2 = M.boost_from(un, q2n)
        dn1 = np.linalg.norm(l1[..., 1:], axis=-1) / p0 - 1
        dn2 = np.linalg.norm(l2[..., 1:], axis=-1) / p0 - 1
        K1 = dn1 - d1[sl][:, None, None]
        K2 = dn2 - d2[sl][:, None, None]
        Lk = ((K1 > a) | (K1 < -a)).astype(float) + ((K2 > a) | (K2 < -a)).astype(float)
        La = ((dn1 > a) | (dn1 < -a)).astype(float) + ((dn2 > a) | (dn2 < -a)).astype(float)
        Lc_ = (((dn1 - dcen[sl][:, None, None]) > a) | ((dn1 - dcen[sl][:, None, None]) < -a)).astype(float) + (((dn2 - dcen[sl][:, None, None]) > a) | ((dn2 - dcen[sl][:, None, None]) < -a)).astype(float)
        wc = (Lc_ * wt[:, :, None]).sum(axis=(1, 2))
        wk = (Lk * wt[:, :, None]).sum(axis=(1, 2))
        wa = (La * wt[:, :, None]).sum(axis=(1, 2))
        # kick の判定は θ の下限が abs 用に下げてあるので、kick 側は a 未満の kick を数えない (L=0) → そのまま。
        out_k[sl] = np.where(ok & keep[sl], vM[sl] * (4 * M.R0 ** 2 / kk ** 4) * wk, 0.0)
        out_a[sl] = np.where(ok & keep[sl], vM[sl] * (4 * M.R0 ** 2 / kk ** 4) * wa, 0.0)
        out_c[sl] = np.where(ok & keep[sl], vM[sl] * (4 * M.R0 ** 2 / kk ** 4) * wc, 0.0)
    fac = al / circ * (N / vol) * 0.5 * M.CLIGHT
    return fac * out_k.mean(), fac * out_a.mean(), float(keep.mean()), fac * out_k.std() / math.sqrt(npairs), fac * out_a.std() / math.sqrt(npairs), fac * out_c.mean()


def main():
    name = sys.argv[1]
    amps = [float(v) for v in sys.argv[2:]] or ([0.004, 0.006] if name == 'astrid2' else [0.006, 0.010])
    from touschek_formula_scan import gather
    from pysad.touschek import piwinski_slice_rate
    from pysad.touschek_bunch import piwinski_slice_rate_bunch, slice_sigma_delta
    tt, tt2, ref = gather(name, 1e9, 0.01)
    info = tt2['rest_info']
    rest = tt2['rest']
    n = len(rest)
    picks = [n // 7, n // 3, (2 * n) // 3]
    print('%s: σδ = %.3e、スライス %s' % (name, ref['SIGE'], picks))
    print(' slice  l    a    a/σδ  1D(δ_c=0)   MC kick     MC abs   MC 増え(±)       理論 増え (σδ²)    (σδ²/2: δ_c のみ)   floor 0.1/0.5')
    for j in picks:
        s = rest[j]
        sd = slice_sigma_delta(s)
        for a in amps:
            r0 = piwinski_slice_rate(s, a, a, info['circ'], 1e9)
            th = piwinski_slice_rate_bunch(s, a, a, info['circ'], 1e9, var_factor=1.0) / r0 - 1
            th2 = piwinski_slice_rate_bunch(s, a, a, info['circ'], 1e9, var_factor=0.5) / r0 - 1
            t1 = piwinski_slice_rate_bunch(s, a, a, info['circ'], 1e9, floor=0.1, var_factor=1.0) / r0 - 1
            t5 = piwinski_slice_rate_bunch(s, a, a, info['circ'], 1e9, floor=0.5, var_factor=1.0) / r0 - 1
            rk, ra, fk, ek, ea, rc = mc_both(s['lab'], s['al'], info['circ'], 1e9, s['p0'], a, npairs=16384)
            inc = ra / rk - 1
            print('   [δ_c 基準の kick (q_s の効果を除く) rate_c = %.4e: kick/c = %+.4f、abs/c = %+.4f]' % (rc, rk / rc - 1, ra / rc - 1))
            einc = (ra / rk) * math.hypot(ek / rk, ea / ra) * 0.5     # 同じ標本の相関で実際の誤差はこれより小さい (保守的な上限)
            print(' %5d %4d %.4f %5.1f %10.4e %10.4e %10.4e  %+.4f (%.4f)   %+.4f           %+.4f         %+.4f/%+.4f' % (
                j, s['l'], a, a / sd, r0, rk, ra, inc, einc, th, th2, t1, t5), flush=True)
    s = rest[picks[1]]
    a = amps[0]
    print('内訳 (診断、スライス %d、a = %g): 相対運動量の δ 成分 (qs)・対の重心の δ (cm) を消した増え (abs/kick − 1)' % (picks[1], a))
    for qs, cm in ((1, 1), (0, 1), (1, 0), (0, 0)):
        rk, ra, fk, ek, ea, rc = mc_both(s['lab'], s['al'], info['circ'], 1e9, s['p0'], a, npairs=16384, qs_scale=qs, cm_scale=cm)
        print('  qs_scale %d cm_scale %d: kick/c %+.4f  abs/c %+.4f  abs/kick %+.4f' % (qs, cm, rk / rc - 1, ra / rc - 1, ra / rk - 1), flush=True)


if __name__ == '__main__':
    main()
