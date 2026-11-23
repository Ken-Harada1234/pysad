# -*- coding: utf-8 -*-
"""v2.4.1: 誤差つき LMA と長時間追跡 (放射・FLUC) の gs・pysad 比較の道具 (台本の生成・gs の実行・出力の解析)。

誤差の与え方: gs の乱数と pysad の乱数は同じ列にできないので、pysad で作った誤差の値 (位置のずれ DX・DY、ロール ROTATE・DROTATE、四極の K1・6 極の K2 の相対誤差)
を、誤差のある実体ごとに固有の名前の要素定義として台本へ明示的に書く (`error_source`)。同じ台本を gs と pysad の台本評価器 (ffscmd.run_sad_string) に渡す。
(`LINE["DX", 番号] = 値;` は同名の全実体に効くので使えない。)
gs は WSL (timeout 長め)。台本名は一意 (`_gs_lma_<番号>_...`)。他のエージェントの sad1.exe は止めない (自分の PID だけ見る)。
"""
import io
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SP = os.path.join(HERE, 'scripts')
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

ASTRID2 = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
KEKLS12 = os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat')

# 誤差の rms (level = 1): 位置ずれは L2・S2 の既定 (四極 50 um・0.2 mrad、6 極 50 um、偏向磁石 50 um・0.1 mrad)、K の相対誤差は四極 1e-3・6 極 1e-3
SIG1 = dict(quad_dxy=50e-6, quad_rot=0.2e-3, sext_dxy=50e-6, bend_dxy=50e-6, bend_rot=0.1e-3)
DK_QUAD = 1e-3
DK_SEXT = 1e-3


def read(p):
    return open(p, encoding='utf-8').read()


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def gs_exe():
    sys.path.insert(0, os.path.join(ROOT, 'example', 'tools'))
    import sadenv
    return sadenv.gs_wsl()


def error_source(src, seed, level=1.0, dk_quad=DK_QUAD, dk_sext=DK_SEXT, line='RING', newline='RINGERR'):
    """src の RING に誤差を与え、誤差のある実体の要素定義と `LINE <newline> = (...)` の文を src の後ろに足した文字列と、情報 (dict) を返す。
    位置ずれ・ロールは fmapres.apply_alignment_errors (streams='type')、K の誤差は種と種類で決まる別の乱数 (numpy)。数値は _raw (SAD の値: 四極の K1・6 極の K2 は積分値)"""
    from pysad import compat, fmapres as FR
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    st0 = TrackSetup.from_source(src, line, rad=False, rfsw=False, alost=0.1)[0]
    sig = {k: v * level for k, v in SIG1.items()}
    st1, cnt = FR.apply_alignment_errors(st0, seed=seed, streams='type', **sig)
    rq = np.random.default_rng([int(seed), 101])
    rs = np.random.default_rng([int(seed), 102])
    defs, names = [], []
    nk = {'quad': 0, 'sext': 0}
    for i, (a, b) in enumerate(zip(st0.els, st1.els)):
        tn = type(a).__name__
        dkq = float(rq.normal(0, 1)) * dk_quad * level if tn == 'Quadrupole' else 0.0
        dks = float(rs.normal(0, 1)) * dk_sext * level if tn == 'Sextupole' else 0.0
        if a is b and dkq == 0.0 and dks == 0.0:
            names.append(a.name)
            continue
        nm = 'E%d' % (i + 1)
        raw = a._raw
        if tn == 'Quadrupole':
            defs.append('QUAD %s = (L=%.17g K1=%.17g DX=%.17g DY=%.17g ROTATE=%.17g);' % (nm, raw['L'], raw['K1'] * (1 + dkq), b.DX, b.DY, b.ROTATION))
            nk['quad'] += 1
        elif tn == 'Sextupole':
            defs.append('SEXT %s = (L=%.17g K2=%.17g DX=%.17g DY=%.17g ROTATE=%.17g);' % (nm, raw['L'], raw['K2'] * (1 + dks), b.DX, b.DY, b.ROTATION))
            nk['sext'] += 1
        elif tn == 'Bend':
            defs.append('BEND %s = (L=%.17g ANGLE=%.17g K1=%.17g E1=%.17g E2=%.17g DX=%.17g DY=%.17g DROTATE=%.17g);' % (
                nm, raw['L'], raw['ANGLE'], raw['K1'], raw.get('E1', 0), raw.get('E2', 0), b.DX, b.DY, b.DROTATE))
        else:
            raise ValueError('未対応の要素 ' + tn)
        names.append(nm)
    text = (src.rstrip() + '\n' + '\n'.join(defs) + '\n' + 'LINE %s = (%s);\n' % (newline, ('\n    ').join(' '.join(names[k:k + 8]) for k in range(0, len(names), 8))))
    return text, dict(count=cnt, kcount=nk, ndef=len(defs))


def run_pysad_text(text, out_path=None, mode='v2'):
    from pysad import compat, ffscmd
    compat.set_mode(mode)
    buf = io.StringIO()
    t0 = time.time()
    ffscmd.run_sad_string(text, out=buf)
    if out_path:
        open(out_path, 'w', encoding='utf-8').write(buf.getvalue())
    return buf.getvalue(), time.time() - t0


def run_gs(sad_path, out_path, timeout=14400):
    """sad_path (一意な名前の台本。_gs_lma_ で始める) を WSL の gs で走らせる。標準出力を out_path に書く。終了後に自分の台本名の sad1.exe の残りを返す"""
    cmd = 'cd "%s" && timeout %d %s %s < /dev/null' % (wsl_path(os.path.dirname(sad_path)), timeout - 60, gs_exe(), os.path.basename(sad_path))
    t0 = time.time()
    p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=timeout,
                       env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
    open(out_path, 'w', encoding='utf-8').write(p.stdout)
    el = time.time() - t0
    ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep | grep %s' % os.path.basename(sad_path)],
                        capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
    return p.stdout, el, ps.stdout.strip()
