# -*- coding: utf-8 -*-
"""pysad の追跡 (trksetup.TrackSetup.track = LMA・周波数マップが使う追跡) と gs (SAD K64) の TrackParticles が、同じ損失判定になることの確認。
ASTRID2 (0.58 GeV、周長 45.7 m、RF 50 kV、RFSW ON・RAD ON・FLUC なし、閉軌道は RADCOD の Emittance[])。位置 loc (要素の入口) から、閉軌道に δ を足した粒子を
NTURN 周追跡し、TrackParticles[{loc, beam}, loc, 1, NTURN] の生き残りの旗 (7 列目) と最後の座標を、gs・pysad の台本・pysad の API の 3 通りで比べる。
損失の限界: LOSSAMPL (±APER)、LOSSDZ (バケットの長さの半分に当たる |z|。LMA の bucket=True と同じ)。δ は、LMA が出した限界の前後。
使い方: py -3.12 example/15_lma_touschek/lma/gs_loss_compare.py [gs を走らせない場合は --nogs]。gs は WSL (timeout 1800 s)、終了後に sad1.exe の残りを確かめる。
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
NTURN = 1312
APER = 0.01
CASES = [(2, [0.0, 0.010, 0.0160, 0.0163, 0.0164, 0.0166, 0.017, 0.020]),      # 長直線部 (RF バケットが限界)。位置の昇順 (台本の出力の順)
         (17, [0.0, 0.007, 0.0082, 0.0083, 0.0084, 0.0086, 0.009, 0.012])]     # 分散の大きい位置 (横のアパーチャが限界)


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def main(run_gs=True):
    from pysad import compat, ffscmd
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    src = open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=APER)
    zb = st.bucket_zlimit()
    print('ASTRID2 0.58 GeV, RFSW ON, RAD ON, FLUC off, nturn=%d, LOSSAMPL=%.3f m, LOSSDZ(bucket)=%.4f m' % (NTURN, APER, zb))
    rows = []
    for loc, dps in CASES:
        cod = st.closed_orbit([loc])[0]
        for sg in (1, -1):
            for d in dps:
                if sg == -1 and d == 0.0:
                    continue
                c = cod.copy()
                c[5] += sg * d
                rows.append((loc, sg * d, c))
    # SAD の台本: 位置ごとに TrackParticles
    lines = ['FFS;', 'USE RING; CELL; CALC; RFSW; RAD; RADCOD; NOFLUC;', 'Emittance[];', 'LOSSAMPL = %.17g; LOSSDZ = %.17g;' % (APER, zb),
             '$FORM = "20.14";', 'NPARA = 1;']
    for loc in sorted(set(r[0] for r in rows)):
        sel = [r for r in rows if r[0] == loc]
        cols = ['{' + ','.join('%.17g' % r[2][k] for r in sel) + '}' for k in range(6)]
        cols.append('{' + ','.join('1' for _ in sel) + '}')
        lines.append('beam = {%s};' % ','.join(cols))
        lines.append('r = TrackParticles[{%d, beam}, %d, 1, %d];' % (loc, loc, NTURN))
        lines.append('Do[Print["P %d ", i, " ", r[[2,7,i]], " ", r[[2,1,i]], " ", r[[2,2,i]], " ", r[[2,3,i]], " ", r[[2,5,i]], " ", r[[2,6,i]]], {i, 1, %d}];' % (loc, len(sel)))
    # 二分法 (純粋な二分、9 回、区間 [0, 0.02]): SAD の台本の中で TrackParticles を 1 粒子ずつ呼ぶ
    for loc in sorted(set(r[0] for r in rows)):
        cod = st.closed_orbit([loc])[0]
        for sg in (1, -1):
            lines.append('cod = {%s}; lo = 0; hi = 0.02;' % ','.join('%.17g' % v for v in cod))
            lines.append('Do[mid = (lo + hi)/2; b1 = {{cod[[1]]}, {cod[[2]]}, {cod[[3]]}, {cod[[4]]}, {cod[[5]]}, {cod[[6]] + (%d)*mid}, {1}};' % sg)
            lines.append('  r = TrackParticles[{%d, b1}, %d, 1, %d]; If[r[[2,7,1]] == 1, lo = mid, hi = mid], {k, 1, 9}];' % (loc, loc, NTURN))
            lines.append('Print["B %d %d ", lo, " ", hi];' % (loc, sg))
    lines.append('Exit;')
    script = '\n'.join(lines) + '\n'
    sad_path = os.path.join(HERE, 'gs_loss_astrid2.sad')
    open(sad_path, 'w', encoding='utf-8').write(src + '\n' + script)
    # pysad の台本
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\n' + script, out=buf)
    py_txt = buf.getvalue()
    open(os.path.join(HERE, 'gs_loss_astrid2.pysad.txt'), 'w', encoding='utf-8').write(py_txt)
    parse = lambda t: [l.split() for l in t.splitlines() if l.startswith('P ')]
    pp = parse(py_txt)
    # pysad の API
    api = []
    for loc in sorted(set(r[0] for r in rows)):
        sel = [r for r in rows if r[0] == loc]
        r = st.track(loc, np.array([x[2] for x in sel]), NTURN, zlost=zb)
        for i in range(len(sel)):
            api.append((loc, sel[i][1], int(r['lost_turn'][i] == 0), r['final'][i]))
    gp = None
    if run_gs:
        gsf = os.path.join(HERE, 'gs_loss_astrid2.gs.txt')
        cmd = 'cd "%s" && timeout 1700 %s gs_loss_astrid2.sad < /dev/null' % (wsl_path(HERE), GS)
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=1800,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        open(gsf, 'w', encoding='utf-8').write(p.stdout)
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        gp = parse(p.stdout)
    print('%4s %10s | %s | %s' % ('loc', 'dp', 'pysad 台本 (生/落)  x', 'API' + ('  gs' if gp else '')))
    nbad = 0
    maxd = 0.0
    for i, (loc, d, c) in enumerate(rows):
        a = pp[i]
        s = 'pysad script: alive=%s x=%s | API: alive=%d' % (a[3], a[4], api[i][2])
        ok = (int(float(a[3])) == api[i][2])
        maxd = max(maxd, abs(float(a[4]) - api[i][3][0]))
        if gp:
            g = gp[i]
            s += ' | gs alive=%s x=%s' % (g[3], g[4])
            ok = ok and int(float(g[3])) == api[i][2]
            maxd = max(maxd, abs(float(g[4]) - api[i][3][0]))
        nbad += (not ok)
        print('%4d %+10.5f | %s %s' % (loc, d, s, '' if ok else '  <-- 不一致'))
    print('不一致 %d / %d、最後の x の最大の差 %.2e m' % (nbad, len(rows), maxd))
    # 二分法: gs の台本の区間 (lo, hi) と、pysad の LMA (lma.momentum_acceptance) の区間が、同じ限界を挟んでいる (区間が重なる)
    from pysad.lma import momentum_acceptance
    locs = sorted(set(r[0] for r in rows))
    res = momentum_acceptance(st, locs, nturn=NTURN, dpmax=0.02, ncoarse=8, nsub=3, tol=2e-5, alost=APER, workers=1)
    pb = lambda t: {(int(l.split()[1]), int(l.split()[2])): (float(l.split()[3]), float(l.split()[4])) for l in t.splitlines() if l.startswith('B ')}
    bpy = pb(py_txt)
    bgs = pb(p.stdout) if gp else {}
    print('二分法 (δ の絶対値): 位置、符号、| pysad の台本 (lo, hi) | LMA (lo, hi) | gs の台本 (lo, hi)')
    nb = 0
    for i, loc in enumerate(res.loc):
        for sg, lo, hi in ((1, res.dp_plus[i], res.hi_plus[i]), (-1, -res.dp_minus[i], -res.hi_minus[i])):
            a = bpy[(loc, sg)]
            line = '%4d %+d | %.6f %.6f | %.6f %.6f' % (loc, sg, a[0], a[1], lo, hi)
            ok = a[0] <= hi and lo <= a[1]
            if gp:
                g = bgs[(loc, sg)]
                line += ' | %.6f %.6f' % g
                ok = ok and g[0] <= hi and lo <= g[1] and g == a
            nb += (not ok)
            print(line, '' if ok else '  <-- 不一致')
    print('二分法の区間の不一致 %d' % nb)


if __name__ == '__main__':
    main(run_gs='--nogs' not in sys.argv)
