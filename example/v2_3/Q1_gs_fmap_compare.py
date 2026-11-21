# -*- coding: utf-8 -*-
"""v2.3.0 Q1: 周波数マップの追跡の部分を gs (SAD K64) と比べる。K64 には NAFF がないので、gs には 1 周ずつ TrackParticles を呼んで (x, px, y, py) を周ごとに
印刷させ (位置 2、NORFSW・NORAD、LOSSAMPL = 0.1 m)、その座標列に pysad の NAFF をかけた D を、pysad の frequency_map (と台本の FrequencyMap[]) の D と比べる。
ASTRID2 (0.58 GeV、周長 45.7 m、6 点: 小振幅から大振幅・共鳴の近く)、512 周 (前半・後半 256 周)。gs は WSL (timeout 1800 s)、終了後に sad1.exe の残りを確かめる。
使い方: py -3.12 example/v2_3/Q1_gs_fmap_compare.py [--nogs]
"""
import io
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
NTURN = 512
LOC = 2
PTS = [(1.0e-3, 1.0e-3), (6.0e-3, 3.0e-3), (10.0e-3, 6.0e-3), (14.0e-3, 9.0e-3), (16.25e-3, 0.5e-3), (18.0e-3, 0.95e-3)]


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def parse_tracks(txt):
    """'T i {x...} {px...} {y...} {py...} {flag...}' の行 → (turn, 4, npart)"""
    out = []
    for chunk in re.split(r'(?m)^T ', txt)[1:]:
        lists = re.findall(r'\{([^}]*)\}', ' '.join(chunk.splitlines()))
        out.append([[float(v) for v in s_.split(',')] for s_ in lists])
    return np.array(out)


def main(run_gs=True):
    from pysad import compat, ffscmd, naff
    from pysad.trksetup import TrackSetup
    from pysad.freqmap import frequency_map
    compat.set_mode('v2')
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    n = len(PTS)
    xs = '{' + ','.join('%.17g' % p[0] for p in PTS) + '}'
    ys = '{' + ','.join('%.17g' % p[1] for p in PTS) + '}'
    zero = '{' + ','.join('0' for _ in PTS) + '}'
    one = '{' + ','.join('1' for _ in PTS) + '}'
    script = ('FFS;\nUSE RING; CELL; CALC; NORFSW; NORAD;\nLOSSAMPL = 0.1;\nbeam = {%s, %s, %s, %s, %s, %s, %s};\nb = beam;\n'
              'Do[r = TrackParticles[{%d, b}, %d]; b = r[[2]]; Print["T ", i, " ", b[[1]], " ", b[[2]], " ", b[[3]], " ", b[[4]], " ", b[[7]]], {i, 1, %d}];\nExit;\n') % (
        xs, zero, ys, zero, zero, zero, one, LOC, LOC, NTURN)
    open(os.path.join(HERE, 'Q1_gs_fmap_astrid2.sad'), 'w', encoding='utf-8').write(src + '\n' + script)
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\n' + script, out=buf)
    open(os.path.join(HERE, 'Q1_gs_fmap_astrid2.pysad.txt'), 'w', encoding='utf-8').write(buf.getvalue())
    tp = parse_tracks(buf.getvalue())
    tg = None
    if run_gs:
        cmd = 'cd "%s" && timeout 1700 %s Q1_gs_fmap_astrid2.sad < /dev/null' % (wsl_path(HERE), GS)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=1800,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        open(os.path.join(HERE, 'Q1_gs_fmap_astrid2.gs.txt'), 'w', encoding='utf-8').write(p.stdout)
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        tg = parse_tracks(p.stdout)
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    tw = st.twiss_at(LOC)
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))
    # frequency_map (点ごとに x, y を 1 点の格子で)
    fm = [frequency_map(st, LOC, x=[p[0]], y=[p[1]], nturn=NTURN) for p in PTS]

    def D_of(series):         # series (turn, 4)
        from pysad.freqmap import _tunes_of
        f, _a = _tunes_of(series, twp, 1)
        return f, np.log10(np.hypot(((f[2] - f[0] + .5) % 1) - .5, ((f[3] - f[1] + .5) % 1) - .5))
    print('ASTRID2 0.58 GeV、位置 %d、NORFSW・NORAD、%d 周 (前半・後半 %d 周)、NAFF 窓 Hann' % (LOC, NTURN, NTURN // 2))
    print('%2s %8s %8s | %-9s %-9s %-9s | %s' % ('i', 'x[mm]', 'y[mm]', 'D(API)', 'D(pysad台本座標)', 'D(gs座標)', 'nux(API) nux(gs) 差'))
    worst_x = 0.0
    worst_D = 0.0
    for i, p in enumerate(PTS):
        f_api = [fm[i].nux1[0], fm[i].nuy1[0], fm[i].nux2[0], fm[i].nuy2[0]]
        s_py = np.array([tp[:, 0, i], tp[:, 1, i], tp[:, 2, i], tp[:, 3, i]]).T
        f_py, D_py = D_of(s_py)
        line = '%2d %8.3f %8.3f | %9.4f %9.4f' % (i, p[0] * 1e3, p[1] * 1e3, fm[i].D[0], D_py)
        if tg is not None:
            s_gs = np.array([tg[:, 0, i], tg[:, 1, i], tg[:, 2, i], tg[:, 3, i]]).T
            f_gs, D_gs = D_of(s_gs)
            worst_x = max(worst_x, np.abs(s_gs - s_py).max())
            worst_D = max(worst_D, abs(D_gs - fm[i].D[0]))
            line += ' %9.4f | %.10f %.10f %.1e' % (D_gs, f_api[0], f_gs[0], abs(f_api[0] - f_gs[0]))
        print(line)
    if tg is not None:
        print('座標列の最大の差 (gs と pysad の台本、全周・全粒子・x px y py): %.2e' % worst_x)
        print('|D(gs の座標) − D(frequency_map)| の最大 = %.2e' % worst_D)


if __name__ == '__main__':
    main(run_gs='--nogs' not in sys.argv)
