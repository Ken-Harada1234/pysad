# -*- coding: utf-8 -*-
"""v2next P2 項目 4: RADCOD ON の 1 周 (または要素) の解析の転送行列 T + D と、写像の数値微分 (4 次の中心差分、刻み 2e-5) の差。
K64 の tradke の微分 (photons.f) は 4 か所が近似・誤りだった (v2 は直した。k64 は K64 のまま):
  (1) 区間の始点の座標に関する微分 tr1 を、D・B の更新 (終点の座標に関する微分を要する) にそのまま使う (スピンの更新は始点の行列 transr を掛けていて正しい)
  (2) dddpx = 0.5((tr2 - dpxr0) dg + ddpx ddg) の 2 項目が dpx·ddg でなく ddpx·ddg (dg 倍だけ小さい)
  (3) x' = x + ddpx c1 の c1 = al1/(3 pr) の微分の項 (ddpx·dc1。さらに δ の項は -ddpx/pr で c1 が欠けている)
  (4) dxpa = (xp·dxp)/(xpa pr²) は pr⁴ で割る (xpa = |xp|/pr²。δ ≠ 0 の軌道で 2δ の相対誤差)
出力は P2_radcod_derivative.out.txt。実行: py -3.12 example/v2_next/P2_radcod_derivative.py"""
import io, math, os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_SAD_PACKAGES', 'D:/work/SAD-master/SAD-master/Packages')
from pysad import ffscmd, compat, emit6d

SRC = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RF = (L=0 VOLT=4e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 24*CELL1 RF);
FFS;
USE RING; CELL; CALC;
'''
C0 = np.array([1e-4, 1e-5, 0, 0, 0, -1e-4])
H = 2e-5
for mode in ('k64', 'v2'):
    compat.set_mode(mode)
    ffs, _ = ffscmd.run_sad_string(SRC, out=io.StringIO())
    els = ffs.lat.elements
    ph = emit6d._Phys(ffs.momentum(), sum(e.L for e in els), amass=emit6d.ELMASS, charge=1.0)

    def f(e_, c0, want=False):
        return emit6d.tturne(e_, np.array(c0, float), ph, trf0=0.0, rfsw=False, radcod=True, want_db=want, want_jac=want)

    def numjac(e_):
        J = np.zeros((6, 6))
        for k in range(6):
            def fe(s):
                c = C0.copy()
                c[k] += s * H
                return f(e_, c).cod
            J[:, k] = (8 * (fe(1) - fe(-1)) - (fe(2) - fe(-2))) / (12 * H)
        return J
    for name, sel in [('QF 1 台', [els[1]]), ('BEND 1 台', [els[3]]), ('QF D BEND', els[1:4]), ('1 セル', els[1:7]), ('24 セルのリング全体', els)]:
        st = f(sel, C0, True)
        d = st.T + st.D - numjac(sel)
        print('%-4s %-22s |T+D - 数値微分| の最大 %.2e (放射の部分 |数値微分 - T| の最大 %.2e)' % (mode, name, np.abs(d).max(), np.abs(numjac(sel) - st.T).max()))
