# -*- coding: utf-8 -*-
"""RFSW ON (RAD OFF) の追跡 (ASTRID2、誤差なし、位置 2) を gs (SAD K64) と比べる。FFS に入ると RFSW は ON なので、gs の既定の台本はこの条件。
6 点 (fmap_ring.py と同じ x, y) を 512 周、1 周ずつ TrackParticles (RFSW・NORAD、LOSSAMPL = 0.1 m) で回し、(x, px, y, py, z, δ) を周ごとに印刷、pysad の同じ追跡
(TrackSetup.track(rfsw=True)) と台本版と比べる。RFSW OFF の同じ追跡 (fmap_ring.py の既定) との差 (z・δ の大きさ) も出す。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.py [--nogs]    gs は WSL (timeout 1700 s)、終了後に sad1.exe の残りを確かめる。出力: fmap_check_gs_rfsw.{sad,gs.txt,pysad.txt}、標準出力
"""
import io
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools')); import sadenv
GS = sadenv.gs_wsl()
NTURN = 512
LOC = 2
PTS = [(1.0e-3, 1.0e-3), (6.0e-3, 3.0e-3), (10.0e-3, 6.0e-3), (14.0e-3, 9.0e-3), (16.25e-3, 0.5e-3), (18.0e-3, 0.95e-3)]


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def parse_tracks(txt):
    out = []
    for chunk in re.split(r'(?m)^T ', txt)[1:]:
        lists = re.findall(r'\{([^}]*)\}', ' '.join(chunk.splitlines()))
        out.append([[float(v) for v in s_.split(',')] for s_ in lists])
    return np.array(out)


def main(run_gs=True):
    from pysad import compat, ffscmd
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    n = len(PTS)
    col = lambda v: '{' + ','.join('%.17g' % t for t in v) + '}'      # noqa
    zero = col([0] * n)
    one = col([1] * n)
    script = ('FFS;\nUSE RING; CELL; CALC; RFSW; NORAD;\nLOSSAMPL = 0.1;\nbeam = {%s, %s, %s, %s, %s, %s, %s};\nb = beam;\n'
              'Do[r = TrackParticles[{%d, b}, %d]; b = r[[2]]; Print["T ", i, " ", b[[1]], " ", b[[2]], " ", b[[3]], " ", b[[4]], " ", b[[5]], " ", b[[6]], " ", b[[7]]], {i, 1, %d}];\nExit;\n') % (
        col([p[0] for p in PTS]), zero, col([p[1] for p in PTS]), zero, zero, zero, one, LOC, LOC, NTURN)
    open(os.path.join(HERE, 'fmap_check_gs_rfsw.sad'), 'w', encoding='utf-8').write(src + '\n' + script)
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\n' + script, out=buf)
    open(os.path.join(HERE, 'fmap_check_gs_rfsw.pysad.txt'), 'w', encoding='utf-8').write(buf.getvalue())
    tp = parse_tracks(buf.getvalue())
    tg = None
    if run_gs:
        cmd = 'cd "%s" && timeout 1700 %s fmap_check_gs_rfsw.sad < /dev/null' % (wsl_path(HERE), GS)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=1800,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        open(os.path.join(HERE, 'fmap_check_gs_rfsw.gs.txt'), 'w', encoding='utf-8').write(p.stdout)
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        tg = parse_tracks(p.stdout)
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=True, alost=0.1)[0]
    coords = np.array([[p[0], 0, p[1], 0, 0, 0] for p in PTS])
    r = st.track(LOC, coords, NTURN, record=True)
    print('位置 %d、RFSW ON・RAD OFF、%d 周。TrackSetup: νs = %.5f、φs = %.4f rad、バケット高さ %.4f' % (
        LOC, NTURN, st.bucket.get('nus', float('nan')), st.bucket.get('phis', float('nan')), st.bucket.get('height', float('nan'))))
    print('%2s %8s %8s | max|z| [mm] max|δ|[1e-4] (TrackSetup) | 座標 (x px y py z δ) の最大の差: 台本−TrackSetup、gs−台本' % ('i', 'x[mm]', 'y[mm]'))
    for i, p in enumerate(PTS):
        s_api = r['rec'][:, i, :6]
        s_py = tp[:, :6, i]
        line = '%2d %8.3f %8.3f | %8.3f %8.3f | %.2e' % (i, p[0] * 1e3, p[1] * 1e3, np.nanmax(np.abs(s_api[:, 4])) * 1e3, np.nanmax(np.abs(s_api[:, 5])) * 1e4, np.nanmax(np.abs(s_py - s_api)))
        if tg is not None:
            line += '  %.2e' % np.abs(tg[:, :6, i] - s_py).max()
        print(line)


if __name__ == '__main__':
    main(run_gs='--nogs' not in sys.argv)
