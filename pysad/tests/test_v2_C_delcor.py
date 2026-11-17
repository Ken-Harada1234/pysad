# -*- coding: utf-8 -*-
"""v2 の修正 B40・B62: DELCOR (相関のあるずれ)。

K64: MARK には DETR・DPX、CAVI には HARM・PHI に (符号も逆に) 書く (添字のずれ)。tvcorr の corlength≠0 の枝は波数 -1 の成分 cv(n) を
設定せず (未初期化)、負の波数の重みが 1 つずれる。
v2: MARK・CAVI も DX・DY に、ほかの磁石と同じ符号で入れる。スペクトルは正負の波数で対称。
根拠: スペクトルの重み w(k) = exp(-(κ x0)²/4)、κ = 2πk/L なので、相関関数は <c(s) c*(s+Δ)> = exp(-Δ²/(2 x0²)) (実数)。
      相関長がリングより十分長ければ、全要素のずれはほぼ同じ値になる (同じ場の値を同じ符号で入れている確認)。
"""
import io
import math

import numpy as np

from pysad.sadrandom import SadRandom


class _Stub:
    def __init__(self, seed):
        self.rng = SadRandom()
        self.rng.restore(seed)

    def gcut(self):
        return 1.0e300


def test_b62_correlation_function_is_gaussian(v2_mode):
    from pysad.ffscmd import FFS
    n, al, x0 = 512, 500.0, 15.0
    st = _Stub(17)
    acc = np.zeros(n, dtype=complex)
    m = 300
    for _ in range(m):
        cv = FFS._tvcorr(st, x0, al, n)
        acc += np.fft.ifft(np.abs(np.fft.fft(cv))**2) / n      # 巡回自己相関 <c_j c*_(j+m)>
    c = acc / m
    c0 = c[0].real
    assert abs(c0 - 2.0) < 0.1                                  # v2 (C2): 実部・虚部それぞれ分散 1 → <|c|²> = 2 (v2_changes_C2_misc.md)
    for k in (1, 5, 10, 20, 30):
        d = k * al / n
        assert abs(c[k].real / c0 - math.exp(-d * d / (2 * x0 * x0))) < 0.05
        assert abs(c[k].imag / c0) < 0.02                       # 対称なスペクトル → 相関関数は実数


SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="20.14";
USE RING;
CELL;
DELCOR 5000 1e-4 2e-4 QF QD RFC ST;
Print["EL ", LINE["DX","QF.1"], " ", LINE["DY","QF.1"], " ", LINE["DX","QD.1"], " ", LINE["DY","QD.1"], " ", LINE["DX","RFC"], " ", LINE["DY","RFC"], " ", LINE["DX","ST"], " ", LINE["DY","ST"], " ", LINE["HARM","RFC"], " ", LINE["DPX","ST"]];
DELCOR ADD 5000 1e-4 2e-4 RFC;
Print["E2 ", LINE["DX","RFC"], " ", LINE["DY","RFC"]];
Exit;
'''


def _run():
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    return {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
            if l.split() and l.split()[0] in ('EL', 'E2')}


def test_b40_cavity_and_mark_get_dx_dy(v2_mode):
    o = _run()
    qx, qy, q2x, q2y, cx, cy, mx, my, harm, dpx = o['EL']
    assert abs(qx) > 1e-7 and abs(qy) > 1e-7
    # 相関長 5 km ≫ 周長 62 m: どの要素もほぼ同じずれ (符号も同じ)
    for v in (q2x, cx, mx):
        assert abs(v / qx - 1.0) < 0.05
    for v in (q2y, cy, my):
        assert abs(v / qy - 1.0) < 0.05
    assert harm == 64.0 and dpx == 0.0
    # ADD: 足される
    assert abs(o['E2'][0] - cx) > 1e-8


def test_b40_k64_mode_writes_harm_and_detr(k64_mode):
    o = _run()
    qx, qy, q2x, q2y, cx, cy, mx, my, harm, dpx = o['EL']
    assert cx == 1e-5 and cy == 2e-5 and mx == 0.0 and my == 0.0     # K64: DX・DY は変わらない
    assert harm != 64.0 and dpx != 0.0
