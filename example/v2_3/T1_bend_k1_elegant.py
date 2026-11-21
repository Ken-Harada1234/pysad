# -*- coding: utf-8 -*-
"""v2.3.0 T1: T1_gs_bend_k1.py の一様な磁場指数のリング (n = 0.4、16 台の扇形磁石) を elegant (CSBEND) で計算し、dν/dδ を出す。
K2 = 0 (K1 が一定の勾配: 厳密な曲線座標)、K2 = −h K1 (SAD と同じ)、K2 = −2 h K1 の 3 通り。elegant は WSL。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/T1_bend_k1_elegant.py   出力: 標準出力 (T1_bend_k1_elegant.out.txt に保存した)
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import T1_gs_bend_k1 as g

ENV = dict(os.environ, MSYS_NO_PATHCONV='1')


def main():
    src, ang, rho, k1 = g.lattice()
    h = 1.0 / rho
    n = g.NIDX
    nux, nuy = math.sqrt(1 - n), math.sqrt(n)
    print('n = %.2f、h = %.6f /m、K1 = %.8f /m^2、ν = (%.6f, %.6f)' % (n, h, k1, nux, nuy))
    print('解析式: K1 一定 (厳密) (%.5f, %.5f)、勾配の力に (1 + h x) を掛けない (%.5f, %.5f)' % (
        -n * (1 + n) / (2 * nux * (1 - n)), nuy * (1 + n) / (2 * (1 - n)), -n * n / (2 * nux * (1 - n)), nuy * n / (2 * (1 - n))))
    for lab, k2 in (('K2 = 0 (厳密な勾配)', 0.0), ('K2 = -h K1 (SAD と同じ)', -h * k1), ('K2 = -2 h K1', -2 * h * k1)):
        lte = ('B: CSBEND, L=%g, ANGLE=%.12f, K1=%.12f, K2=%.12f, N_KICKS=32, EDGE_ORDER=2\nRING: LINE=(%s)\n' % (g.LB, ang, k1, k2, ','.join(['B'] * g.N)))
        ele = ('&run_setup lattice=u.lte, use_beamline=RING, p_central_mev=3000, default_order=3, rootname=u, &end\n&run_control &end\n'
               '&twiss_output filename=%s.twi, radiation_integrals=1, higher_order_chromaticity=1 &end\n&bunched_beam n_particles_per_bunch=1 &end\n&track &end\n')
        d = os.path.join(HERE, '..', '..', 'local', 'scratch_T1')
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'u.lte'), 'w', newline='\n').write(lte)
        open(os.path.join(d, 'u.ele'), 'w', newline='\n').write(ele)
        wp = os.path.abspath(d).replace(os.sep, '/')
        wp = '/mnt/' + wp[0].lower() + wp[2:]
        cmd = ('source ~/elegant/env.sh; mkdir -p ~/t1/bendk1; cd ~/t1/bendk1; cp "%s/u.lte" "%s/u.ele" .; elegant u.ele > log.txt 2>&1; '
               'sddsprintout u.twi -par=dnux/dp -par=dnuy/dp') % (wp, wp)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=300, env=ENV)
        import re
        v = re.findall(r'dnu[xy]/dp[^=]*=\s*([-0-9.e+]+)', p.stdout)
        print('elegant CSBEND %-26s dν/dδ = (%.5f, %.5f)' % (lab, float(v[0]), float(v[1])))


if __name__ == '__main__':
    main()
