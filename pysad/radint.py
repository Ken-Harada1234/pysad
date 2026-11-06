"""
radint.py — FFS コマンド RADINT (SAD src/INTGRL.f の直訳)
==========================================================

SAD の RADINT は temit/Emittance[] とは別の、x-y 結合を含む 4D 放射積分 (I1x,I1y,I2,I3,I4u,I4v,I5u,I5v + 色収差)。
直前の CALC の Twiss (twiss(j,0,*) = DP0 の行) を入口条件に、BEND/QUAD/SEXT を 2*IR 個のスライスに分けて
Twiss (R 行列による結合の基準化込み) を伝搬しながら Simpson 則で積分する。

直訳にあたって再現した SAD の癖 (INTGRL.f は 1991 年のコードで、モジュール変数を使い回す):
  1. DL / CKL* / TM / TO / TETA / FUNI などはモジュール変数。SEXT のブロックは DL を設定せず、直前の BEND/QUAD の
     最後のスライス長 DL をそのまま TCAL に使う (FUNI(9)=GU*DL, FUNI(10)=GV*DL が六極でも足される = 色収差の積分に効く)。
     DKL2 は代入されないので六極の CKL0X/CKL1X は常に 0。RLE=0 の六極は TCAL の後で DL=0.2 にしてメッセージを出す。
  2. tbendebody (BEND の 1 スライスの行列と分散ベクトル) は bendib モジュールの rhob/rho0 を使うが INTGRL は設定しない。
     つまり「直前に tbend/tbende/tbendi/tsteee が設定した値」(= 最後に処理した曲率のある BEND の L/ANGLE) が全 BEND に使われる。
     pysad では RadInt.stale_rho (rhob, rho0) として持つ (既定は最後の曲率のある BEND)。
  3. 四極のスライス数は int(|K1L|/0.02) (1〜14)、BEND は int(|ANGLE|/0.002) (1〜14)。ANGLE=0 の BEND、K1=0 の四極は読み飛ばし
     (RLT には長さだけ足される)。
  4. CAVI は VOLT/HARM/FREQ を最後のもので上書き (合計しない)。HARM=0 なら FREQ*RLT/c、FREQ=0 なら HARM*c/RLT。
"""
import math
from typing import List

import numpy as np

from . import emit6d as _E
from .maintrack import fmt_g, fmt_f

NRI = 12
ML = 30
# tffitcode (1 始まり。TWPMG[0] は使わない)
AX, BX, NX, AY, BY, NY, EX, EPX, EY, EPY, R1, R2, R3, R4, DETR, DXC, DPXC, DYC, DPYC = range(1, 20)

# SAD の macphys 定数 (tfstk.f)。elradi, elmass は emit6d と共通
ELMASS = _E.ELMASS
ELRADI = _E.ELRADI
FINEST = _E.FINEST
CVELOC = _E.CVELOC
def _f32(x):
    """Fortran の単精度定数 (INTGRL.f は 1e27, 0.2, 1.E-6, 2.3026 を倍精度の式に単精度のまま書いている)"""
    return float(np.float32(x))


COEFF = (ELRADI / (ELMASS ** 3) / 3.0) * _f32(1e27)       # 1e27 は単精度 (9.99999988e26): ELOSS が相対 1.15e-8 小さくなる
CO4 = math.sqrt(55.0 / 32.0 / math.sqrt(3.0) * ELRADI / FINEST) / ELMASS * 1e9
CO5 = CO4 ** 2


def _symps(z, n):
    """SYMPS.f (z は 0 始まりの配列、Fortran の Z(1)..Z(2N+1))"""
    if n == 0:
        return z[0]
    s2 = 0.0
    s3 = 0.0
    for i in range(1, n):
        s2 = s2 + z[i * 2]
    for i in range(1, n + 1):
        s3 = s3 + z[i * 2 - 1]
    return (z[0] + z[2 * n] + 2.0 * s2 + 4.0 * s3) / 3.0


class RadInt:
    """INTGRL モジュール (module radint) の状態。プロセス中で使い回される (SAD も同じ: 2 回目の RADINT の DL の初期値は前回の最後の値)"""

    def __init__(self):
        self.DL = 0.0
        self.DKL0 = 0.0
        self.DKL1 = 0.0
        self.DKL2 = 0.0
        self.CKL0X = self.CKL0Y = self.CKL1X = self.CKL1Y = 0.0
        self.TETA = 0.0
        self.TM = np.zeros((4, 4))
        self.VM = np.zeros(4)
        self.TO = np.zeros((4, 4))
        self.VO = np.zeros(4)
        self.TU = np.zeros((4, 4))
        self.VU = np.zeros(4)
        self.R0 = np.zeros((4, 4))          # TCAL の save 変数
        self.R0I = np.zeros((4, 4))
        self.FUNI = np.zeros((NRI + 1, ML + 1))   # FUNI[i][k]: i=1..NRI, k=1..ML (1 始まり)
        self.TWPMG = np.zeros(29)
        self.dp0 = 0.0
        self.msgs: List[str] = []

    # ------------------------------------------------------------------ 2x2 / 4x4 の小道具 (DM42, MINV2, CM24 ...)
    @staticmethod
    def _minv2(a):
        det = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]
        if det == 0.0:
            raise ZeroDivisionError('DETERMINANT=0 MATRIX-INV NOT FOUND')
        return np.array([[a[1, 1] / det, -a[0, 1] / det], [-a[1, 0] / det, a[0, 0] / det]])

    @staticmethod
    def _symtr2(a):
        return np.array([[-a[1, 1], a[0, 1]], [a[1, 0], -a[0, 0]]])

    @staticmethod
    def _cm24(a, b, c, d):
        r = np.zeros((4, 4))
        r[0:2, 0:2] = a
        r[0:2, 2:4] = b
        r[2:4, 0:2] = c
        r[2:4, 2:4] = d
        return r

    @staticmethod
    def _r1inv(r0):
        ri = np.zeros((4, 4))
        ri[0:2, 0:2] = r0[0:2, 0:2]
        ri[2:4, 2:4] = r0[2:4, 2:4]
        ri[0:2, 2:4] = -r0[0:2, 2:4]
        ri[2:4, 0:2] = -r0[2:4, 0:2]
        return ri

    @staticmethod
    def _r2inv(r0):
        ri = np.zeros((4, 4))
        ri[0:2, 0:2] = -r0[2:4, 2:4]
        ri[2:4, 2:4] = -r0[0:2, 0:2]
        ri[0:2, 2:4] = r0[0:2, 2:4]
        ri[2:4, 0:2] = r0[2:4, 0:2]
        return ri

    @staticmethod
    def _unit2(x):
        return np.array([[x, 0.0], [0.0, x]])

    # ------------------------------------------------------------------ R0CAL / TMATRS / TWTRANS / TROT
    def r0cal(self):
        tw = self.TWPMG
        rc0 = np.array([[tw[R1], tw[R2]], [tw[R3], tw[R4]]])
        detr = rc0[0, 0] * rc0[1, 1] - rc0[0, 1] * rc0[1, 0]
        if detr <= 1.0:
            rmu0 = math.sqrt(1.0 - detr)
            rb0 = self._symtr2(rc0)
            r0 = self._cm24(self._unit2(rmu0), rb0, rc0, self._unit2(rmu0))
            return r0, self._r1inv(r0)
        # DET(R)>1
        ra0 = np.array([[tw[13] / detr, tw[14] / detr], [tw[11] / detr, tw[12] / detr]])
        rd0 = self._symtr2(ra0)
        rmu0 = math.sqrt(1.0 - 1.0 / detr)
        r0 = self._cm24(ra0, self._unit2(rmu0), self._unit2(rmu0), rd0)
        return r0, self._r2inv(r0)

    def tmatrs(self, R0, R0I):
        """R0, R0I は呼び出し側の配列を書き換える (TCAL は save 変数、TEDGE は局所)"""
        tw = self.TWPMG
        TR = self.TO @ R0I
        A, C = TR[0:2, 0:2], TR[2:4, 0:2]
        AI = self._minv2(A)
        CAI = C @ AI
        detcai = CAI[0, 0] * CAI[1, 1] - CAI[0, 1] * CAI[1, 0]
        if detcai < 1.0:
            smu1 = 1.0 / (1.0 + detcai)
            rmu1 = math.sqrt(smu1)
            RC1 = (-rmu1) * CAI
            RB1 = self._symtr2(RC1)
            RA1 = self._unit2(rmu1)
            RD1 = self._unit2(rmu1)
            R0[:, :] = self._cm24(RA1, RB1, RC1, RD1)
            R0I[:, :] = self._r1inv(R0)
            for j in (1, 2):
                for k in (1, 2):
                    tw[10 + k + 2 * (j - 1)] = RC1[j - 1, k - 1]
        else:
            smu1 = 1.0 / (1.0 + 1.0 / detcai)
            rmu1 = math.sqrt(smu1)
            ACI = self._minv2(CAI)
            RD1 = (-rmu1) * ACI
            RA1 = self._symtr2(RD1)
            RB1 = self._unit2(rmu1)
            RC1 = self._unit2(rmu1)
            R0[:, :] = self._cm24(RA1, RB1, RC1, RD1)
            R0I[:, :] = self._r2inv(R0)
            for j in (1, 2):
                for k in (1, 2):
                    tw[10 + k + 2 * (2 - j)] = -RD1[j - 1, k - 1]
        self.TU = R0 @ TR
        self.VU = R0 @ self.VO

    def twtrans(self):
        tw = self.TWPMG
        TU = self.TU
        ttx = np.array([
            [TU[0, 0] * TU[0, 0], -2.0 * TU[0, 0] * TU[0, 1], TU[0, 1] * TU[0, 1]],
            [-TU[0, 0] * TU[1, 0], 1.0 + 2.0 * TU[0, 1] * TU[1, 0], -TU[0, 1] * TU[1, 1]],
            [TU[1, 0] * TU[1, 0], -2.0 * TU[1, 1] * TU[1, 0], TU[1, 1] * TU[1, 1]]])
        tty = np.array([
            [TU[2, 2] * TU[2, 2], -2.0 * TU[2, 2] * TU[2, 3], TU[2, 3] * TU[2, 3]],
            [-TU[2, 2] * TU[3, 2], 1.0 + 2.0 * TU[2, 3] * TU[3, 2], -TU[2, 3] * TU[3, 3]],
            [TU[3, 2] * TU[3, 2], -2.0 * TU[3, 3] * TU[3, 2], TU[3, 3] * TU[3, 3]]])
        ax, bx, ay, by = tw[AX], tw[BX], tw[AY], tw[BY]
        ex, epx, ey, epy = tw[EX], tw[EPX], tw[EY], tw[EPY]
        gx = (1.0 + ax * ax) / bx
        gy = (1.0 + ay * ay) / by
        bx1 = ttx[0, 0] * bx + ttx[0, 1] * ax + ttx[0, 2] * gx
        ax1 = ttx[1, 0] * bx + ttx[1, 1] * ax + ttx[1, 2] * gx
        by1 = tty[0, 0] * by + tty[0, 1] * ay + tty[0, 2] * gy
        ay1 = tty[1, 0] * by + tty[1, 1] * ay + tty[1, 2] * gy
        VU = self.VU
        ex1 = TU[0, 0] * ex + TU[0, 1] * epx + VU[0]
        epx1 = TU[1, 0] * ex + TU[1, 1] * epx + VU[1]
        ey1 = TU[2, 2] * ey + TU[2, 3] * epy + VU[2]
        epy1 = TU[3, 2] * ey + TU[3, 3] * epy + VU[3]
        tw[AX], tw[AY], tw[BX], tw[BY] = ax1, ay1, bx1, by1
        tw[EX], tw[EPX], tw[EY], tw[EPY] = ex1, epx1, ey1, epy1
        dx1 = self.TO @ tw[DXC:DPYC + 1]
        for i in range(4):
            tw[DXC + i] = dx1[i] + self.VO[i] * self.dp0

    @staticmethod
    def trot(tm, vm, teta):
        ctet = math.cos(teta)
        stet = math.sin(teta)
        rot = np.eye(4) * ctet
        rinv = np.eye(4) * ctet
        rot[0, 2] = -stet
        rot[1, 3] = -stet
        rot[2, 0] = stet
        rot[3, 1] = stet
        rinv[0, 2] = stet
        rinv[1, 3] = stet
        rinv[2, 0] = -stet
        rinv[3, 1] = -stet
        ts = tm @ rot
        return rinv @ ts, rinv @ vm

    # ------------------------------------------------------------------ TCAL
    def tcal(self, i, iflg):
        """FUNI(i,*) を現在の TWPMG で作り、iflg=0 なら次の点へ Twiss を伝搬する。i は 1 始まり"""
        tw = self.TWPMG
        DL = np.float64(self.DL)
        with np.errstate(all='ignore'):
            if DL == 0.0:
                self.msgs.append('  DL=0. IN TCAL ')
            au, bu, av, bv = tw[AX], tw[BX], tw[AY], tw[BY]
            eu, epu, ev, epv = tw[EX], tw[EPX], tw[EY], tw[EPY]
            RK1 = self.CKL1X / DL
            CURV = math.hypot(self.CKL0X, self.CKL0Y) / DL
            gu = (1.0 + au * au) / bu
            gv = (1.0 + av * av) / bv
            if i == 1:
                self.R0[:, :], self.R0I[:, :] = self.r0cal()
            R0I = self.R0I
            CURV2 = CURV * CURV
            CURV3 = CURV2 * abs(CURV)
            CKL0X, CKL0Y, CKL1X, CKL1Y = self.CKL0X, self.CKL0Y, self.CKL1X, self.CKL1Y
            F = self.FUNI
            F[1, i] = CKL0X * (R0I[0, 0] * eu + R0I[0, 1] * epu + R0I[0, 2] * ev + R0I[0, 3] * epv)
            F[2, i] = CURV2 * DL
            F[3, i] = CURV3 * DL
            if CURV > 0:
                F[4, i] = ((CKL0X + 2.0 * CKL1X / CURV) * CURV2 * (R0I[0, 0] * eu + R0I[0, 1] * epu)
                           + (CKL0Y - 2.0 * CKL1Y / CURV) * CURV2 * (R0I[2, 0] * eu + R0I[2, 1] * epu))
            else:
                F[4, i] = 0.0
            F[5, i] = CURV3 * (gu * eu * eu + 2.0 * au * epu * eu + bu * epu * epu) * DL
            F[6, i] = CKL0Y * (R0I[2, 0] * eu + R0I[2, 1] * epu + R0I[2, 2] * ev + R0I[2, 3] * epv)
            if CURV > 0:
                F[7, i] = ((CKL0X + 2.0 * CKL1X / CURV) * CURV2 * (R0I[0, 2] * ev + R0I[0, 3] * epv)
                           + (CKL0Y - 2.0 * CKL1Y / CURV) * CURV2 * (R0I[2, 2] * ev + R0I[2, 3] * epv))
            else:
                F[7, i] = 0.0
            F[8, i] = CURV3 * (gv * ev * ev + 2.0 * av * epv * ev + bv * epv * epv) * DL
            F[9, i] = gu * DL
            F[10, i] = gv * DL
            F[11, i] = CKL0X
            F[12, i] = CKL0Y
            if iflg == 1:
                return
            self.tmatrs(self.R0, self.R0I)
            self.twtrans()

    # ------------------------------------------------------------------ 要素ごとのスライス行列
    def bmag(self, teta1, rhob, rho0):
        """BMAG: tbendebody(trans, cod=0, al=DL, phin=DKL0, ak1=DKL1) の 1 次の行列 (TM) と分散ベクトル (VM)"""
        m = _E._BendIBody(self.DKL1, self.DL, 0.0, rhob, rho0)
        p = 1.0
        tm = np.zeros((4, 4))
        vm = np.zeros(4)
        tm[0, 0] = m.cx
        tm[0, 1] = m.sxkx / p
        vm[0] = m.drhopp * m.dcxkx + m.drhop * m.dcxkxp
        tm[1, 0] = m.aksx * p
        tm[1, 1] = m.cx
        vm[1] = p * (m.drhopp * m.sxkx + m.drhop * m.sxkxp) + m.drhop * m.sxkx
        tm[2, 2] = m.cy
        tm[2, 3] = m.syky / p
        tm[3, 2] = m.aksy * p
        tm[3, 3] = m.cy
        self.TM = tm
        self.VM = vm
        self.CKL0X = math.cos(teta1) * self.DKL0
        self.CKL0Y = -math.sin(teta1) * self.DKL0
        self.CKL1X = math.cos(teta1) * self.DKL1
        self.CKL1Y = -math.cos(teta1) * self.DKL1
        self.TO, self.VO = self.trot(self.TM, self.VM, teta1)

    def quadru(self, teta1, dx, dy):
        self.VM = np.zeros(4)
        tm = np.eye(4)
        DL = self.DL
        rk = self.DKL1 / DL
        sqrk = math.sqrt(rk) if rk >= 0.0 else math.sqrt(-rk)
        p = sqrk * DL
        s1 = math.sin(p)
        s2 = math.cos(p)
        s3 = math.sinh(p)
        s4 = math.cosh(p)
        if rk >= 0.0:
            tm[0, 0] = s2
            tm[0, 1] = s1 / sqrk
            tm[1, 0] = -sqrk * s1
            tm[1, 1] = s2
            tm[2, 2] = s4
            tm[2, 3] = s3 / sqrk
            tm[3, 2] = sqrk * s3
            tm[3, 3] = s4
        else:
            tm[0, 0] = s4
            tm[0, 1] = s3 / sqrk
            tm[1, 0] = sqrk * s3
            tm[1, 1] = s4
            tm[2, 2] = s2
            tm[2, 3] = s1 / sqrk
            tm[3, 2] = -sqrk * s1
            tm[3, 3] = s2
        self.TM = tm
        tet2 = teta1 * 2.0
        self.CKL0Y = -self.DKL1 * (dx * math.sin(tet2) + dy * math.cos(tet2))
        self.CKL0X = self.DKL1 * (dx * math.cos(tet2) - dy * math.sin(tet2))
        dr = math.hypot(dx, dy)
        if dr > 0.0:
            self.CKL1X = self.DKL1 * dx / dr
            self.CKL1Y = -self.DKL1 * dy / dr
        else:
            self.CKL1X = self.DKL1
            self.CKL1Y = -self.DKL1
        self.TO, self.VO = self.trot(self.TM, self.VM, teta1)

    def bedge(self, eang):
        curv = self.DKL0 / self.DL
        self.VM = np.zeros(4)
        tm = np.eye(4)
        tm[1, 0] = curv * math.tan(eang)
        tm[3, 2] = -curv * math.tan(eang)
        self.TM = tm
        self.TO, self.VO = self.trot(self.TM, self.VM, self.TETA)
        self.CKL0X = math.cos(self.TETA) * self.DKL0
        self.CKL0Y = -math.sin(self.TETA) * self.DKL0
        self.CKL1X = 0.0
        self.CKL1Y = 0.0

    def sextu(self, teta1, dx, dy):
        self.TM = np.eye(4)
        self.TO = np.eye(4)
        self.VM = np.zeros(4)
        self.VO = np.zeros(4)
        sint = math.sin(3.0 * teta1)
        cost = math.cos(3.0 * teta1)
        self.CKL0Y = -self.DKL2 * (dx * dy * cost + (dx * dx - dy * dy) / 2.0 * sint)
        self.CKL0X = self.DKL2 * ((dx * dx - dy * dy) / 2.0 * cost - dx * dy * sint)
        self.CKL1X = self.DKL2 * dx
        self.CKL1Y = -self.DKL2 * dy

    def tedge(self, dieg):
        tw = self.TWPMG
        au, bu, av, bv = tw[AX], tw[BX], tw[AY], tw[BY]
        eu, epu, ev, epv = tw[EX], tw[EPX], tw[EY], tw[EPY]
        curvx = self.CKL0X / self.DL
        curvy = self.CKL0Y / self.DL
        R0, R0I = self.r0cal()
        TO, TM = self.TO, self.TM
        dieg[4] = dieg[4] - (TO[1, 0] * curvx * (R0I[0, 0] * eu + R0I[0, 1] * epu)
                             - TO[3, 2] * curvy * (R0I[2, 0] * eu + R0I[2, 1] * epu))
        dieg[7] = dieg[7] - (TO[1, 0] * curvx * (R0I[0, 2] * ev + R0I[0, 3] * epv)
                             - TO[3, 2] * curvy * (R0I[2, 2] * ev + R0I[2, 3] * epv))
        dieg[9] = dieg[9] - TM[1, 0] * bu
        dieg[10] = dieg[10] - TM[3, 2] * bv
        self.tmatrs(R0, R0I)
        self.twtrans()


_STATE = RadInt()
_DEBUG = None        # デバッグ用: list を入れると要素ごとの DI (I=1..10) を溜める (SAD に write を入れたデバッグ版との比較用)


def default_stale_rho(els):
    """INTGRL が使う bendib の rhob/rho0 の既定値: 直前に tbend 系が設定した値 = 最後に処理した曲率のある BEND の (L/ANGLE)"""
    rho = None
    for el in els:
        if type(el).__name__ == 'Bend' and float(getattr(el, 'ANGLE', 0.0) or 0.0) != 0.0 and el.L != 0.0:
            rho = el.L / el.ANGLE
    return (rho, rho) if rho is not None else (1.0, 1.0)


def compute(ffs, rows, stale_rho=None):
    """INTGRL の本体。戻り: 結果の辞書 (RI, 補助量) と出力する行のリスト"""
    g = _STATE
    els = ffs.els
    nlat = ffs.nlat
    g.msgs = []
    g.dp0 = ffs.dp0()
    rhob, rho0 = stale_rho if stale_rho is not None else default_stale_rho(els)
    p0 = ffs.momentum()
    amass = float(ffs.I.get('MASS', ELMASS))            # SAD の amass = 大域変数 MASS (既定 elmass)
    ph = _E._Phys(p0, 1.0, amass=amass)                  # _Phys の第 1 引数は運動量 [eV]
    energy = amass / 1e9 * ph.h0
    RI = np.zeros(NRI + 1)
    rfv = 0.0
    rlt = 0.0
    totang = 0.0
    rharm = 0.0
    freq = 0.0
    out = []
    for j in range(1, nlat):
        el = els[j - 1]
        cls = type(el).__name__
        twpin = rows[j - 1]
        vec = np.zeros(29)
        vec[AX], vec[BX], vec[NX] = twpin.alfx, twpin.betx, twpin.mux
        vec[AY], vec[BY], vec[NY] = twpin.alfy, twpin.bety, twpin.muy
        vec[EX], vec[EPX], vec[EY], vec[EPY] = twpin.etx, twpin.etpx, twpin.ety, twpin.etpy
        vec[R1], vec[R2], vec[R3], vec[R4], vec[DETR] = twpin.r1, twpin.r2, twpin.r3, twpin.r4, twpin.detr
        vec[DXC], vec[DPXC], vec[DYC], vec[DPYC] = twpin.x_cod, twpin.px_cod, twpin.y_cod, twpin.py_cod
        g.TWPMG = vec
        di = np.zeros(NRI + 1)
        dieg = np.zeros(NRI + 1)
        g.CKL0X = g.CKL0Y = g.CKL1X = g.CKL1Y = 0.0
        rle = float(el.L)
        rlt = rlt + rle
        ir = 0
        did = False
        if cls == 'Bend':
            ang = ffs.get_value(el, 'ANGLE')
            e1 = ffs._e_ratio(el, 'E1') * ang
            e2 = ffs._e_ratio(el, 'E2') * ang
            from .lattice import _bend_theta
            g.TETA = _bend_theta(el)
            rkl1 = ffs.get_value(el, 'K1')
            if ang != 0.0:
                ir = min(14, max(1, int(abs(ang) / 0.002)))
                irp = 2 * ir
                g.DL = rle / irp
                g.DKL0 = ang / irp
                g.DKL1 = rkl1 / irp
                totang = totang + ang
                if e1 != 0.0:
                    g.bedge(e1)
                    g.tedge(dieg)
                for i in range(1, irp + 1):
                    g.bmag(g.TETA, rhob, rho0)
                    g.tcal(i, 0)
                g.tcal(irp + 1, 1)
                if e2 != 0.0:
                    g.bedge(e2)
                    g.tedge(dieg)
                did = True
        elif cls == 'Quadrupole':
            rkl1 = ffs.get_value(el, 'K1')
            g.TETA = float(getattr(el, 'ROTATION', 0.0) or 0.0)
            if rkl1 != 0.0:
                ir = int(abs(rkl1) / 0.02)
                if ir == 0:
                    ir = 1
                if ir > 14:
                    ir = 14
                irp = 2 * ir
                g.DL = rle / irp
                g.DKL1 = rkl1 / irp
                for i in range(1, irp + 1):
                    g.quadru(g.TETA, g.TWPMG[DXC], g.TWPMG[DYC])
                    g.tcal(i, 0)
                g.tcal(irp + 1, 1)
                did = True
        elif cls == 'Sextupole':
            g.TETA = float(getattr(el, 'ROTATION', 0.0) or 0.0)
            g.sextu(g.TETA, g.TWPMG[DXC], g.TWPMG[DYC])
            ir = 0
            g.tcal(1, 0)
            if rle == 0.0:
                g.DL = _f32(0.2)
                g.msgs.append(' L(SEXTUPOLE) EQ 0. SET L=0.2')
            did = True
        elif cls == 'Cavity':
            rfv = ffs.get_value(el, 'VOLT')
            rharm = ffs.get_value(el, 'HARM')
            freq = ffs.get_value(el, 'FREQ')
            raw = ffs.env._elements.get(el.name, {})
            if float(raw.get('FREQ', 0.0) or 0.0) == 0.0 and rharm != 0.0:
                freq = 0.0                     # HARM だけの定義 (pysad の parser が FREQ を補った分は SAD の値 0 に戻す)
        if did:
            for i in range(1, NRI + 1):
                di[i] = _symps(g.FUNI[i, 1:], ir) + dieg[i]
                RI[i] = RI[i] + di[i]
            if _DEBUG is not None:
                _DEBUG.append((j, {'Bend': 2, 'Quadrupole': 4, 'Sextupole': 6}.get(cls, 0), [float(di[i]) for i in range(1, 11)]))
    last = rows[nlat - 1]
    qx = last.mux / (2.0 * math.pi)
    qy = last.muy / (2.0 * math.pi)
    res = dict(RI=RI, energy=energy, rlt=rlt, totang=totang, rfv=rfv, rharm=rharm, freq=freq, qx=qx, qy=qy)
    return res, g.msgs


def _fd(x, w, d):
    """Fortran の 1P,Dw.d"""
    if x != x:
        return 'NaN'.rjust(w)
    if math.isinf(x):
        return ('Infinity' if x > 0 else '-Infinity').rjust(w)
    mant, e = f'{x:.{d}E}'.split('E')
    return (mant + f'D{int(e):+03d}').rjust(w)


def report(ffs, rows, stale_rho=None) -> List[str]:
    """INTGRL の出力 (1 行 1 要素のリスト、末尾の改行なし)"""
    res, msgs = compute(ffs, rows, stale_rho)
    RI = res['RI']
    energy = res['energy']
    rlt = res['rlt']
    rfv = res['rfv']
    rharm = res['rharm']
    freq = res['freq']
    dp0 = ffs.dp0()
    pi = math.pi
    L = []
    w = L.append
    stars = '*' * 64
    w(stars)
    w('')
    w('     RADIATION INTEGRALS OF A LATTICE WITH X-Y COUPLING')
    w('')
    w('     CONTRIBUTIONS FROM SOLENOID MAGNETS ARE NOT included ')
    w('              IN THIS VERSION.   1991.2.26 ')
    w('                  Refurbished:  2019.10.10 ')
    w('             Added K1 of BEND:  2022.08.05 ')
    w(stars)
    # DL=0. や L(SEXTUPOLE) のメッセージは計算中に出るので、本来は見出しの直後 (SYNCHROTRON の前) に並ぶ
    L.extend(msgs)
    w('')
    w('')
    w('    SYNCHROTRON RADIATION INTEGRALS')
    if rharm == 0.0:
        rharm = freq * rlt / CVELOC
    if freq == 0.0:
        freq = rharm * CVELOC / rlt
    w('')
    w(' ENERGY (GeV)               = ' + fmt_f(energy, 9, 3))
    w(' DP/P                       = ' + fmt_f(dp0, 9, 3))
    w(' TOTAL LENGTH (m)           = ' + fmt_f(rlt, 9, 3))
    totang = res['totang'] * 180.0 / pi
    w(' TOTAL BENDING ANGLE (DEG)  = ' + fmt_f(totang, 9, 3))
    w(' TOTAL BENDING ANGLE (RAD)  = ' + fmt_g(RI[11], 15, 7))
    rfv = rfv / 1000.0
    w(' RF VOLTAGE (kV)            = ' + fmt_g(rfv, 15, 7))
    freq = freq * _f32(1e-6)
    w(' HARMONIC NUMBER            = ' + (f'{rharm:.0f}.').rjust(6))
    w(' RF FREQUENCY (MHz)         = ' + fmt_f(freq, 10, 4))
    w('')
    w(' RADIATION INTEGRALS  ')
    w('   I1X = ' + fmt_g(RI[1], 15, 7) + '    I1Y = ' + fmt_g(RI[6], 15, 7))
    w('   I2  = ' + fmt_g(RI[2], 15, 7))
    w('   I3  = ' + fmt_g(RI[3], 15, 7))
    w('   I4U = ' + fmt_g(RI[4], 15, 7) + '    I4V = ' + fmt_g(RI[7], 15, 7))
    w('   I5U = ' + fmt_g(RI[5], 15, 7) + '    I5V = ' + fmt_g(RI[8], 15, 7))
    with np.errstate(all='ignore'):
        eloss = 2.0e6 * COEFF * RI[2] * energy ** 4
        alphax = RI[1] / rlt
        alphay = RI[6] / rlt
        dampu = 1.0 - RI[4] / RI[2]
        dampv = 1.0 - RI[7] / RI[2]
        co2 = CVELOC * energy ** 3 * COEFF / rlt
        tauu = 1000.0 / (co2 * (RI[2] - RI[4]))
        tauv = 1000.0 / (co2 * (RI[2] - RI[7]))
        taue = 1000.0 / (co2 * (2 * RI[2] + RI[4] + RI[7]))
        sig = CO4 * energy * math.sqrt(RI[3] / (2.0 * RI[2] + RI[4] + RI[7]))
        emitu = CO5 * energy * energy * RI[5] / (RI[2] - RI[4])
        emitv = CO5 * energy * energy * RI[8] / (RI[2] - RI[7])
        chrx = -RI[9] / 4.0 / pi
        chry = -RI[10] / 4.0 / pi
    w('')
    w(' NU-X-  = ' + fmt_f(res['qx'], 8, 5) + '    NU-Y- = ' + fmt_f(res['qy'], 8, 5))
    w(' MOMENTUM COMPACTION FAC X  = ' + fmt_g(alphax, 12, 4))
    w(' MOMENTUM COMPACTION FAC Y  = ' + fmt_g(alphay, 12, 4))
    w(' ENERGY LOSS PER TURN (kV)  = ' + fmt_g(eloss, 12, 4))
    w(' DU=1-I4U/I2                = ' + fmt_f(dampu, 12, 8))
    w(' DV=1-I4V/I2                = ' + fmt_f(dampv, 12, 8))
    w('')
    w(' DAMPING TIME CONSTANTS')
    w('                TAUU (msec) = ' + fmt_g(tauu, 12, 4))
    w('                TAUV (msec) = ' + fmt_g(tauv, 12, 4))
    w('                TAUE (msec) = ' + fmt_g(taue, 12, 4))
    w(' ENERGY SPREAD (SIG) DE/E   = ' + fmt_g(sig, 12, 4))
    w('            EMITTANCE U     = ' + fmt_g(emitu, 12, 4))
    w('            EMITTANCE V     = ' + fmt_g(emitv, 12, 4))
    if eloss <= rfv:
        siph = eloss / rfv
        coph = math.sqrt(1.0 - siph * siph)
        alpha = alphax + alphay
        qs = 0.00039894228 * math.sqrt(rharm * alpha * rfv * coph / energy)
        sigl = 1.0e6 * sig * rlt * math.sqrt((alpha * energy) / (2 * pi * rharm * rfv * coph))
        overv = rfv / eloss
        overv = math.sqrt(overv * overv - 1)
        rkq = 2.0 * (overv - math.atan(overv))
        epmax = math.sqrt(1.0e-6 * eloss * rkq / (pi * alpha * rharm * energy))
        xi = 0.5 * (epmax / sig) ** 2
        w('')
        w(' NU-S                       = ' + fmt_f(qs, 12, 6))
        w(' BUCKET HIGHT   DEMAX/E     = ' + fmt_f(epmax, 12, 6))
        w(' NATURAL BUNCH LENGTH (mm)  = ' + fmt_f(sigl, 12, 6))
        if xi < 100:
            tauq = taue / 1000.0 * math.exp(xi) / 2.0 / xi
            try:
                tauq = math.exp(_f32(2.3026) * tauq)
            except OverflowError:
                tauq = float('inf')
            w(' QUANTUM LIFETIME(hours)    = ' + _fd(tauq, 12, 4))
        else:
            w(' Q.LIFE IS TOO LONG      XI = ' + _fd(xi, 12, 4))
    else:
        w('  NOT ENOUGH RF VOLTS;  ELOSS,RFV = ' + _lst(eloss, True) + _lst(rfv) + ' ')
    w('')
    w(' LINEAR CHROMATICITY GX     = ' + fmt_f(chrx, 12, 6))
    w(' LINEAR CHROMATICITY GY     = ' + fmt_f(chry, 12, 6))
    return L


def _lst(x, first=False):
    """gfortran のリスト入力出力 (real*8) = G25.17 (0.1<=|x|<1e17 は F21.(17-n) + 4 空白)。項目の間に 1 空白 (文字列の直後の最初の項目には付かない)。
    SAD の実出力 `= ` + 3 空白 + 678.50075327596767 + 8 空白 + 0.0000000000000000 + 5 空白 に一致させた"""
    return ('' if first else ' ') + fmt_g(x, 25, 17)
