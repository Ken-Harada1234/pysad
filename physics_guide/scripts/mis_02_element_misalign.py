# -*- coding: utf-8 -*-
"""要素の取り付け誤差 (11_misalignment_and_coord.md §6〜§8)。

(a) DX・DY・ROTATE: QUAD の追跡 = 「(x-DX, y-DY) を R(+theta) で回し、誤差なしの QUAD を通し、R(-theta) で戻し、DX・DY を足す」と一致するか。
    回す向きの符号 (sigma = +1 か -1) を数値で決める。
(b) 正準性: 4x4 の回転 diag(R, R) は x-px、y-py の対を保つ。シンプレクティックな行列の確認。
(c) BEND の DROTATE・CHI2 (tbrot): 剛体回転 (軸 = 弦の方向のロールと x 軸まわりの CHI2) を、弦の中点を支点に独立に作り、tbrot と比べる。
(d) tsolrot の Bz = 0 の枝 (CHI1・CHI2・DZ。tsol.f の転記) を独立な 3 次元の幾何と比べる。pysad の追跡は直線の MULT の CHI1・CHI2・DZ を未実装。
(e) MARK・APERT・MAP は力学ではない: MARK は恒等、APERT は座標を変えず失う粒子を決めるだけ。
"""
import io
import math
import numpy as np
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6

print('pysad:', C.pysad_src())


def mk(el, name='EL'):
    src = 'MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST %s);\nFFS;\nUSE ONE;\n' % (el, name)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def rot4(c, s):
    """(x, px, y, py, ...) の x-y を R = [[c, -s], [s, c]] で回す (px, py も同じ R)。z, delta は変えない"""
    o = c.copy()
    o[0] = c_ * 0
    return o


def R4(th, cod):
    c, s = math.cos(th), math.sin(th)
    o = np.array(cod, dtype=float)
    o[0], o[2] = c * cod[0] - s * cod[2], s * cod[0] + c * cod[2]
    o[1], o[3] = c * cod[1] - s * cod[3], s * cod[1] + c * cod[3]
    return o


print('== (a) DX・DY・ROTATE の QUAD ==')
c0 = np.array([1e-3, 2e-4, -1e-3, 1e-4, 0.01, 2e-3])
base = mk('QUAD EL=(L=0.3 K1=0.5);')
for dx, dy, th in ((1e-3, 2e-3, 0.0), (0.0, 0.0, 0.1), (1e-3, 2e-3, 0.1)):
    f = mk('QUAD EL=(L=0.3 K1=0.5 DX=%g DY=%g ROTATE=%g);' % (dx, dy, th))
    tr = C.track_one(f, c0)
    best = None
    for sg in (1, -1):
        a = c0.copy()
        a[0] -= dx
        a[2] -= dy
        a = R4(sg * th, a)
        b = C.track_one(base, a)
        b = R4(-sg * th, b)
        b[0] += dx
        b[2] += dy
        e = np.abs(b - tr).max()
        print('DX=%g DY=%g ROTATE=%g: 回す角 sigma*ROTATE の sigma = %+d: 差 %.2e' % (dx, dy, th, sg, e))
    # 順序 (回してから変位) との比較
    a = R4(th, c0)
    a[0] -= dx
    a[2] -= dy
    b = C.track_one(base, a)
    b[0] += dx
    b[2] += dy
    b = R4(-th, b)
    print('   (参考) 先に回し後で DX・DY を引く順序: 差 %.2e' % np.abs(b - tr).max())

print('\n== (b) 回転の 4x4 はシンプレクティック ==')
th = 0.37
c, s = math.cos(th), math.sin(th)
R = np.array([[c, -s], [s, c]])
M = np.zeros((4, 4))
M[0, 0] = M[2, 2] = M[1, 1] = M[3, 3] = c
M[0, 2] = M[1, 3] = -s
M[2, 0] = M[3, 1] = s
S4 = np.zeros((4, 4))
S4[0, 1] = S4[2, 3] = 1
S4[1, 0] = S4[3, 2] = -1
print('|M^T S M - S| = %.1e (x, px, y, py の順。x-y と px-py を同じ R で回す)' % np.abs(M.T @ S4 @ M - S4).max())

print('\n== (c) BEND の tbrot (DROTATE・CHI2): 弦の中点まわりの剛体回転 ==')


def axis_rot(n, ang):
    n = np.asarray(n, dtype=float)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)


def check_tbrot(ent, dtheta, dchi2, alg=0.5, phig=0.1):
    """独立な構成: 支点 (弦の中点) まわりの剛体回転 R = R_n(dchi2) R_t(dtheta)。t = 弦の方向 (ロール)、n = 弦に垂直で水平な軸 (CHI2)。
    t, n の符号と、合成の順序・R か R^T かを総当たりし、最小の差を出す (軸の向きは tbrot.f の読みから決める。数値で確かめる)"""
    ph = _Ph()
    p = track6._Particle(c0, ph)
    track6._tbrot(p, alg, phig, dtheta, dchi2, ent)
    out = p.out()
    pr = 1.0 + c0[5]
    pxs, pys = c0[1], c0[3]
    pzs = math.sqrt(pr * pr - pxs * pxs - pys * pys)
    r = np.array([c0[0], c0[2], 0.0])
    P = np.array([pxs, pys, pzs])
    s0, c0_ = math.sin(phig), math.cos(phig)
    dl = alg if phig == 0.0 else alg * s0 / phig
    piv = np.array([-dl * s0, 0.0, dl * c0_])
    res = {}
    for ts in (1, -1):
        for ns in (1, -1):
            for at in (1, -1):
                for an in (1, -1):
                    t = np.array([ts * s0, 0.0, c0_])
                    n = np.array([c0_, 0.0, ns * s0])
                    Rt = axis_rot(t, at * dtheta)
                    Rn = axis_rot(n, an * dchi2)
                    for nm, Rm in (('Rn Rt', Rn @ Rt), ('Rt Rn', Rt @ Rn)):
                        r1 = piv + Rm @ (r - piv)
                        P1 = Rm @ P
                        u1 = P1 / pr
                        zz = r1[2] / u1[2]
                        X1 = r1 - zz * u1
                        cand = np.array([X1[0], P1[0], X1[1], P1[1], c0[4] + zz, c0[5]])
                        res[(nm, 't=%+d' % ts, 'n=%+d' % ns, 'roll=%+d' % at, 'chi2=%+d' % an)] = np.abs(cand - out).max()
    k = min(res, key=res.get)
    return res[k], k


class _Ph:
    dvfs = 0.0

    def dv(self, d):
        return 0.0


def _phys():
    return _Ph()


try:
    for ent in (True, False):
        for dth, dch in ((0.02, 0.0), (0.0, 0.015), (0.02, 0.015)):
            e, k = check_tbrot(ent, dth, dch)
            print('ent=%s DROTATE=%g CHI2=%g: 最小の差 %.2e (候補 %s)' % (ent, dth, dch, e, k))
except Exception as ex:      # pysad の内部の型が合わないとき
    print('(c) 実行できず:', repr(ex))

print('\n== (d) CHI1・CHI2・DZ (直線の MULT の tsolrot) ==')
ff = mk('MULT EL=(L=1 K1=0.3 DISFRIN=1 CHI1=0.01 CHI2=0.02 DZ=0.02);')
f0 = mk('MULT EL=(L=1 K1=0.3 DISFRIN=1);')
c = np.array([1e-3, 2e-4, -1e-3, 1e-4, 0.01, 2e-3])
print('CHI1・CHI2・DZ を指定した MULT と指定なしの MULT の追跡の差 = %.1e (0 なら pysad の追跡は CHI1・CHI2・DZ を使っていない)' % np.abs(C.track_one(ff, c) - C.track_one(f0, c)).max())

print('\n== (e) MARK・APERT・MAP ==')
f = mk('MARK EL=();')
c = np.array([1e-3, 2e-4, -1e-3, 1e-4, 0.01, 2e-3])
print('MARK: 出力 - 入力 の最大 = %.1e' % np.abs(C.track_one(f, c) - c).max())
f = mk('APERT EL=(AX=0.01 AY=0.01);')
print('APERT (AX = AY = 1 cm): 通過する粒子の出力 - 入力 の最大 = %.1e' % np.abs(C.track_one(f, c) - c).max())
