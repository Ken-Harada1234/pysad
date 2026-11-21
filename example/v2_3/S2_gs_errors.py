# -*- coding: utf-8 -*-
"""v2.3.0 S2 (7): 誤差つき ASTRID2 の追跡を gs (SAD K64) と比べる。
gs の乱数の種は pysad の乱数と同じにできない (SAD の Gauss 乱数は SAD 自身の列で、pysad の周波数マップの誤差は numpy の乱数) ので、乱数の種は渡さず、
**pysad で作った誤差の値 (DX・DY・四極の ROTATE・偏向磁石の DROTATE) を、誤差のある実体ごとに固有の名前の要素定義として gs の台本へ明示的に書く** (error_ring)。
(LINE["DX", 番号] = 値; は同じ名前の全部の実体に効くので使えない。)これで同じ誤差の配置を gs が使う。
初期条件は位置 2 の閉軌道 (pysad の Newton 法) + (x, y) の絶対座標。gs の台本には、CALC のあとの閉軌道 (Twiss["DX", 2]・Twiss["DY", 2]) も印刷させる。
512 周を 1 周ずつ TrackParticles (NORFSW・NORAD、LOSSAMPL = 0.1 m、位置 2) で回し、(x, px, y, py) を周ごとに印刷、pysad の同じ追跡 (TrackSetup.track) と台本版と比べる。
使い方: py -3.12 example/v2_3/S2_gs_errors.py [--nogs] [--seed 1] [--level 1]     gs は WSL (timeout 1700 s)、終了後に sad1.exe の残りを確かめる。
出力: S2_gs_errors.{sad,gs.txt,pysad.txt}、標準出力
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
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools')); import sadenv
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


def error_ring(st0, st1, line='RINGERR'):
    """st0 (誤差なし) と st1 (誤差つき) の要素を比べ、変わった要素ごとに、固有の名前 (E<番号>) の要素定義と、全要素を並べた LINE <line> の SAD の文を返す。
    注意: gs (と pysad の台本) の LINE["DX", 番号] = 値; は、同じ名前の要素 (同じ定義から作られた全ての実体、ASTRID2 の QF3 は 12 個) の全部に効く
    (実体ごとの値は書けない) ので、実体ごとに別の定義を作る。数値は _raw (SAD の値: 四極の K1・6 極の K2・偏向磁石の K1 は積分値) を使う"""
    defs = []
    names = []
    for i, (a, b) in enumerate(zip(st0.els, st1.els)):
        if a is b:
            names.append(a.name)
            continue
        nm = 'E%d' % (i + 1)
        raw = a._raw
        tn = type(b).__name__
        if tn == 'Quadrupole':
            defs.append('QUAD %s = (L=%.17g K1=%.17g DX=%.17g DY=%.17g ROTATE=%.17g);' % (nm, raw['L'], raw['K1'], b.DX, b.DY, b.ROTATION))
        elif tn == 'Sextupole':
            defs.append('SEXT %s = (L=%.17g K2=%.17g DX=%.17g DY=%.17g ROTATE=%.17g);' % (nm, raw['L'], raw['K2'], b.DX, b.DY, b.ROTATION))
        elif tn == 'Bend':
            defs.append('BEND %s = (L=%.17g ANGLE=%.17g K1=%.17g E1=%.17g E2=%.17g DX=%.17g DY=%.17g DROTATE=%.17g);' % (
                nm, raw['L'], raw['ANGLE'], raw['K1'], raw.get('E1', 0), raw.get('E2', 0), b.DX, b.DY, b.DROTATE))
        else:
            raise ValueError('未対応の要素 ' + tn)
        names.append(nm)
    return chr(10).join(defs) + chr(10) + 'LINE %s = (%s);' % (line, (chr(10) + '    ').join(' '.join(names[k:k + 8]) for k in range(0, len(names), 8))) + chr(10), len(defs)


def main(run_gs=True):
    from pysad import compat, ffscmd, naff, fmapres as FR
    from pysad.trksetup import TrackSetup
    from pysad.freqmap import _tunes_of
    from S2_errors_scan import kind_sigmas
    compat.set_mode('v2')
    seed = int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 1
    level = float(sys.argv[sys.argv.index('--level') + 1]) if '--level' in sys.argv else 1.0
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    st0 = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    st1, cnt = FR.apply_alignment_errors(st0, seed=seed, streams='type', **kind_sigmas('astrid2', 'all', level))
    orb = st1.closed_orbit_dp(LOC, 0.0)
    ringdef, ndef = error_ring(st0, st1)
    coords = np.zeros((len(PTS), 6))
    coords[:] = orb
    for i, p in enumerate(PTS):
        coords[i, 0] += p[0]
        coords[i, 2] += p[1]
    col = lambda k: '{' + ','.join('%.17g' % coords[i, k] for i in range(len(PTS))) + '}'      # noqa
    zero = '{' + ','.join('0' for _ in PTS) + '}'
    one = '{' + ','.join('1' for _ in PTS) + '}'
    script = (ringdef + '\nFFS;\nUSE RINGERR; CELL; CALC; NORFSW; NORAD;\nLOSSAMPL = 0.1;\nCALC;\n'
              'Print["ORB ", Twiss["DX", %d], " ", Twiss["DY", %d]];\n'
              'beam = {%s, %s, %s, %s, %s, %s, %s};\nb = beam;\n'
              'Do[r = TrackParticles[{%d, b}, %d]; b = r[[2]]; Print["T ", i, " ", b[[1]], " ", b[[2]], " ", b[[3]], " ", b[[4]], " ", b[[7]]], {i, 1, %d}];\nExit;\n') % (
        LOC, LOC, col(0), col(1), col(2), col(3), zero, zero, one, LOC, LOC, NTURN)
    script = script.replace('b = beam;\n', 'b = beam;\n')
    open(os.path.join(HERE, 'S2_gs_errors.sad'), 'w', encoding='utf-8').write(src + '\n' + script)
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\n' + script, out=buf)
    open(os.path.join(HERE, 'S2_gs_errors.pysad.txt'), 'w', encoding='utf-8').write(buf.getvalue())
    tp = parse_tracks(buf.getvalue())
    m = re.search(r'(?m)^ORB\s+(\S+)\s+(\S+)', buf.getvalue())
    print('誤差を与えた要素 %s (種 %d、強さ ×%g)。位置 %d の閉軌道 x, y [mm]: pysad の Newton 法 (%.5f, %.5f)、pysad の台本 CALC (%s)' % (
        cnt, seed, level, LOC, orb[0] * 1e3, orb[2] * 1e3, ', '.join('%.5f' % (float(v) * 1e3) for v in m.groups()) if m else '?'))
    tg = None
    if run_gs:
        cmd = 'cd "%s" && timeout 1700 %s S2_gs_errors.sad < /dev/null' % (wsl_path(HERE), GS)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=1800,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        open(os.path.join(HERE, 'S2_gs_errors.gs.txt'), 'w', encoding='utf-8').write(p.stdout)
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        tg = parse_tracks(p.stdout)
        mg = re.search(r'(?m)^ORB\s+(\S+)\s+(\S+)', p.stdout)
        print('gs の CALC の閉軌道 x, y [mm]:', ', '.join('%.5f' % (float(v) * 1e3) for v in mg.groups()) if mg else '?(gs の出力に ORB なし)')
    r = st1.track(LOC, coords, NTURN, record=True)
    tw = st1.twiss_at(LOC)
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))

    def D_of(series):
        f, _a = _tunes_of(series, twp, 1)
        return f, np.log10(np.hypot(((f[2] - f[0] + .5) % 1) - .5, ((f[3] - f[1] + .5) % 1) - .5))
    print('%2s %8s %8s | D(TrackSetup) D(pysad 台本) D(gs) | 座標の最大の差 (台本−TrackSetup) (gs−台本)' % ('i', 'x[mm]', 'y[mm]'))
    for i, p in enumerate(PTS):
        s_api = r['rec'][:, i, :4]
        s_py = tp[:, :4, i] if tp.ndim == 3 else None
        d_api = D_of(s_api)[1] if np.all(np.isfinite(s_api)) else float('nan')
        d_py = D_of(s_py)[1]
        line = '%2d %8.3f %8.3f | %10.4f %10.4f' % (i, p[0] * 1e3, p[1] * 1e3, d_api, d_py)
        dd = np.nanmax(np.abs(s_py - s_api))
        if tg is not None:
            s_gs = tg[:, :4, i]
            line += ' %10.4f | %.2e %.2e' % (D_of(s_gs)[1], dd, np.abs(s_gs - s_py).max())
        else:
            line += ' | %.2e' % dd
        print(line)


if __name__ == '__main__':
    main(run_gs='--nogs' not in sys.argv)
