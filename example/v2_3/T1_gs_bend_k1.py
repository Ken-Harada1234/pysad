# -*- coding: utf-8 -*-
"""v2.3.0 T1: K1≠0 の BEND (複合機能) のオフモメンタムの光学を gs (SAD K64) と pysad で比べ、解析式・elegant と比べる。

リング: 一様な磁場指数 n = 0.4 の扇形偏向磁石 16 台だけ (L = 1 m、ANGLE = 2π/16、K1 = −n/ρ²、ρ = L/ANGLE)、3 GeV、RF なし・放射なし (ここでは磁石の本体だけを見る)。
解析式 (磁場 B_y = B0 + G x が曲線座標の半径方向に厳密に直線 = K1 が一定の勾配): 閉軌道の半径が δ で変わり、磁場指数が n_eff = n (1 + x/R0) / (1 − n x/R0)、x = R0 δ/(1 − n)
 → ν_y² = n_eff、ν_x² = 1 − n_eff。dν_y/dδ = ν_y (1 + n) / (2 (1 − n))、dν_x/dδ = −n (1 + n) / (2 ν_x (1 − n))。n = 0.4 で (−0.6025、+0.7379)。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/T1_gs_bend_k1.py [--no-gs]
手順: 4 次元の閉軌道を pysad で δ = −0.4 … +0.4% で求め、その閉軌道から x・y に 1e-6 ずらした粒子を gs の TrackParticles と pysad で 400 周追跡し、NAFF でチューンの端数を出す。
出力: 標準出力 (T1_gs_bend_k1.out.txt に保存)。gs は WSL (timeout 600 s)、終了後に sad1.exe の残りを確かめる。
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools')); import sadenv
GS = sadenv.gs_wsl()
N, LB, NIDX = 16, 1.0, 0.4
NT = 400


def lattice():
    ang = 2 * math.pi / N
    rho = LB / ang
    k1 = -NIDX / rho ** 2
    return ('MOMENTUM = 3 GEV;\nMARK M0 = ();\nBEND B = (L = %g ANGLE = %.12f K1 = %.12f*%g);\nLINE RING = (M0 %d*B);\n' % (LB, ang, k1, LB, N)), ang, rho, k1


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def main(run_gs=True):
    from pysad import compat, resdrive
    from pysad.trksetup import TrackSetup
    from pysad.naff import naff
    compat.set_mode('v2')
    src, ang, rho, k1 = lattice()
    n = NIDX
    nux, nuy = math.sqrt(1 - n), math.sqrt(n)
    ex_x = -n * (1 + n) / (2 * nux * (1 - n))
    ex_y = nuy * (1 + n) / (2 * (1 - n))
    # pysad が (K64 と同じく) 曲率の因子 (1 + h x) を勾配の力に掛けない場合: n_eff = n (1 + n x/R0) → 係数 (1 + n)/(1 − n) が n/(1 − n) になる
    sad_x = -n * n / (2 * nux * (1 - n))
    sad_y = nuy * n / (2 * (1 - n))
    print('一様な磁場指数 n = %.2f の 16 台の扇形磁石 (3 GeV)。ν = (%.6f, %.6f)' % (n, nux, nuy))
    print('解析式 (K1 = 一定の勾配、曲線座標で厳密): dν/dδ = (%.5f, %.5f)' % (ex_x, ex_y))
    print('曲率の因子 (1 + h x) を勾配の力に掛けない模型: dν/dδ = (%.5f, %.5f)' % (sad_x, sad_y))
    st, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=1.0)
    r = resdrive.chromaticity_numeric(ffs.lat)
    print('pysad (CALC のオフモメンタム、v2): dν/dδ = (%.5f, %.5f)' % tuple(r['d1']))
    dps = [-0.004, -0.002, 0.0, 0.002, 0.004]
    a = 1e-6
    C = []
    for d in dps:
        c0 = st.closed_orbit_dp(1, d)
        c = c0.copy()
        c[0] += a
        c[2] += a
        c[5] = d
        C.append(c)
    C = np.array(C)
    rec = st.track(1, C, NT, record=True)['rec']
    fp = np.array([[abs(naff(rec[:, k, i] - rec[:, k, i].mean(), 1)[0][0]) for i in (0, 2)] for k in range(len(dps))])
    tg = None
    if run_gs:
        def lst(col):
            return '{' + ','.join('%.17g' % v for v in col) + '}'
        beam = '{%s,%s,%s,%s,%s,%s,%s}' % (lst(C[:, 0]), lst(C[:, 1]), lst(C[:, 2]), lst(C[:, 3]), lst(C[:, 4]), lst(C[:, 5]), lst([1] * len(dps)))
        script = ('FFS;\nUSE RING; CELL; CALC; NORFSW; NORAD;\nLOSSAMPL = 1;\nb = %s;\nDo[r = TrackParticles[{1, b}, 1]; b = r[[2]]; '
                  'Print["T ", i, " ", b[[1]], " ", b[[2]], " ", b[[3]], " ", b[[4]], " ", b[[7]]], {i, 1, %d}];\nExit;\n') % (beam, NT)
        open(os.path.join(HERE, 'T1_gs_bend_k1.sad'), 'w', encoding='utf-8').write(src + '\n' + script)
        cmd = 'cd "%s" && timeout 600 %s T1_gs_bend_k1.sad < /dev/null' % (wsl_path(HERE), GS)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=900,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        out = []
        for chunk in re.split(r'(?m)^T ', p.stdout)[1:]:
            lists = re.findall(r'\{([^}]*)\}', ' '.join(chunk.splitlines()))
            out.append([[float(v) for v in s_.split(',')] for s_ in lists])
        tg = np.array(out)
    print('\n追跡 (閉軌道の周りの 1e-6 m の振動の NAFF、%d 周): チューンの端数 |f| (x、y)' % NT)
    for k, d in enumerate(dps):
        line = 'δ = %+.3f %%  pysad %.10f %.10f' % (d * 100, fp[k, 0], fp[k, 1])
        if tg is not None:
            fg = [abs(naff(tg[:, i, k] - tg[:, i, k].mean(), 1)[0][0]) for i in (0, 2)]
            line += '   gs %.10f %.10f   差 %.1e' % (fg[0], fg[1], max(abs(fg[0] - fp[k, 0]), abs(fg[1] - fp[k, 1])))
        print(line)
    # 追跡から dν/dδ (|f| は端数の絶対値。ν_x = 1 − |f_x|、ν_y = 1 − |f_y|)
    print('追跡の dν/dδ: ν_x = 1 − |f_x| なので dν_x/dδ = %.5f、ν_y = 1 − |f_y| なので dν_y/dδ = %.5f' % (-np.polyfit(dps, fp[:, 0], 1)[0], -np.polyfit(dps, fp[:, 1], 1)[0]))


if __name__ == '__main__':
    main(run_gs='--no-gs' not in sys.argv)
