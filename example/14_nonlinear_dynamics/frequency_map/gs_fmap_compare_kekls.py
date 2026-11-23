# -*- coding: utf-8 -*-
"""周波数マップの追跡の部分を gs (SAD K64) と比べる (gs_fmap_compare_astrid2.py の一般化。KEK-LS の kek-ls_v03_12.lat でも)。
K64 には NAFF がないので、gs には 1 周ずつ TrackParticles を呼んで (x, px, y, py, 生存の旗) を周ごとに印刷させ (位置 loc、NORFSW・NORAD、LOSSAMPL = 0.1 m)、
その座標列に pysad の NAFF をかけた D を、pysad の frequency_map の D と比べる。失われる点では、落ちた周 (旗が 0 になる最初の周) も比べる。
点は fmap_ring.py の表 (fmap_<ring>.table.txt) から D の値で選ぶ (生き残る 6 点 + 失われる 2 点)。512 周 (前半・後半 256 周)。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/gs_fmap_compare_kekls.py astrid2|kekls [--nogs] [--timeout 1800]
gs は WSL (timeout 付き)、終了後に sad1.exe の残りを確かめる (残りがあれば PID で kill する)。出力: gs_fmap_<ring>.{sad,gs.txt,pysad.txt}、標準出力。
"""
import io
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools')); import sadenv
GS = sadenv.gs_wsl()
NTURN = 512
RINGS = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), loc=2, p=None, table=os.path.join(HERE, 'fmap_astrid2.table.txt')),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), loc=3, p=3e9, table=os.path.join(HERE, 'fmap_kekls.table.txt')),
}


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def parse_tracks(txt):
    """'T i {x...} {px...} {y...} {py...} {flag...}' の行 → (turn, 5, npart)"""
    out = []
    for chunk in re.split(r'(?m)^T ', txt)[1:]:
        lists = re.findall(r'\{([^}]*)\}', ' '.join(chunk.splitlines()))
        out.append([[float(v) for v in s_.split(',')] for s_ in lists])
    return np.array(out)


def choose_points(tab):
    """Q1 の表から、生き残る 6 点 (D が目標 -9.5, -8, -6.5, -5, -4, -3 に近い点) と、失われる 2 点 (落ちた周が最も早い点と、最も遅い点) を選ぶ"""
    from pysad.freqmap import FMapResult
    r = FMapResult.load_table(tab)
    ok = np.isfinite(r.D)
    idx = []
    for tgt in (-9.5, -8.0, -6.5, -5.0, -4.0, -3.0):
        c = np.where(ok, np.abs(r.D - tgt), 1e9)
        for j in idx:
            c[j] = 1e9
        idx.append(int(np.argmin(c)))
    lost = np.nonzero(r.lost_turn > 0)[0]
    if len(lost):
        idx.append(int(lost[np.argmin(r.lost_turn[lost])]))
        idx.append(int(lost[np.argmax(r.lost_turn[lost])]))
    return [(r.x[i], r.y[i], r.D[i], int(r.lost_turn[i])) for i in idx]


def main(ring, run_gs=True, timeout=1800):
    from pysad import compat, ffscmd
    from pysad.trksetup import TrackSetup
    from pysad.freqmap import frequency_map, _tunes_of
    compat.set_mode('v2')
    c = RINGS[ring]
    LOC = c['loc']
    src = open(c['path'], encoding='utf-8').read()
    pts = choose_points(c['table'])
    n = len(pts)
    xs = '{' + ','.join('%.17g' % p[0] for p in pts) + '}'
    ys = '{' + ','.join('%.17g' % p[1] for p in pts) + '}'
    zero = '{' + ','.join('0' for _ in pts) + '}'
    one = '{' + ','.join('1' for _ in pts) + '}'
    mom = ('MOMENTUM = %r;\n' % c['p']) if c['p'] else ''
    script = ('%sFFS;\nUSE RING; CELL; CALC; NORFSW; NORAD;\nLOSSAMPL = 0.1;\nbeam = {%s, %s, %s, %s, %s, %s, %s};\nb = beam;\n'
              'Do[r = TrackParticles[{%d, b}, %d]; b = r[[2]]; Print["T ", i, " ", b[[1]], " ", b[[2]], " ", b[[3]], " ", b[[4]], " ", b[[7]]], {i, 1, %d}];\nExit;\n') % (
        mom, xs, zero, ys, zero, zero, zero, one, LOC, LOC, NTURN)
    fsad = os.path.join(HERE, 'gs_fmap_%s.sad' % ring)
    open(fsad, 'w', encoding='utf-8').write(src + '\n' + script)
    t0 = time.time()
    buf = io.StringIO()
    ffscmd.run_sad_string(src + '\n' + script, out=buf)
    print('pysad の台本 (run_sad_string): %.0f s' % (time.time() - t0))
    open(os.path.join(HERE, 'gs_fmap_%s.pysad.txt' % ring), 'w', encoding='utf-8').write(buf.getvalue())
    tp = parse_tracks(buf.getvalue())
    tg = None
    if run_gs:
        t0 = time.time()
        cmd = 'cd "%s" && timeout %d %s %s < /dev/null' % (wsl_path(HERE), timeout - 100, GS, os.path.basename(fsad))
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=timeout,
                           env=dict(os.environ, MSYS_NO_PATHCONV='1'), encoding='utf-8', errors='replace')
        print('gs: %.0f s' % (time.time() - t0))
        open(os.path.join(HERE, 'gs_fmap_%s.gs.txt' % ring), 'w', encoding='utf-8').write(p.stdout)
        ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                            capture_output=True, text=True, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        print('gs 終了後の sad1.exe の残り:', ps.stdout.strip() or 'なし')
        tg = parse_tracks(p.stdout)
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
    tw = st.twiss_at(LOC)
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))
    fm = [frequency_map(st, LOC, x=[p[0]], y=[p[1]], nturn=NTURN) for p in pts]

    def D_of(series):
        f, _a = _tunes_of(series, twp, 1)
        return f, np.log10(np.hypot(((f[2] - f[0] + .5) % 1) - .5, ((f[3] - f[1] + .5) % 1) - .5))

    def lost_turn(arr, i):
        z = np.nonzero(arr[:, 4, i] == 0)[0]
        return int(z[0]) + 1 if len(z) else 0
    print('%s 位置 %d、NORFSW・NORAD、%d 周 (前半・後半 %d 周)、NAFF 窓 Hann。点は Q1 の表 (1024 周の D) から選んだ。' % (ring, LOC, NTURN, NTURN // 2))
    print('%2s %8s %8s | %-9s %-9s %-9s | %-12s %-12s %s | 落ちた周 (pysad API / pysad 台本 / gs)' % ('i', 'x[mm]', 'y[mm]', 'D(API)', 'D(pysad台本)', 'D(gs座標)', 'nux(API)', 'nux(gs)', '差'))
    worst_x = 0.0
    worst_x_reg = 0.0
    worst_D = 0.0
    worst_D_reg = 0.0
    nlost_match = 0
    nlost_all = 0
    for i, p in enumerate(pts):
        lp = lost_turn(tp, i)
        lg = lost_turn(tg, i) if tg is not None else None
        api = fm[i]
        la = int(api.lost_turn[0])
        if lp > 0 or la > 0 or (lg is not None and lg > 0):
            nlost_all += 1
            nlost_match += int(la == lp and (lg is None or lg == lp))
            print('%2d %8.3f %8.3f | (失われる点: Q1 の 1024 周では落ちた周 %d) 落ちた周: API %d / 台本 %d / gs %s' % (i, p[0] * 1e3, p[1] * 1e3, p[3], la, lp, lg))
            continue
        s_py = np.array([tp[:, 0, i], tp[:, 1, i], tp[:, 2, i], tp[:, 3, i]]).T
        f_py, D_py = D_of(s_py)
        line = '%2d %8.3f %8.3f | %9.4f %9.4f' % (i, p[0] * 1e3, p[1] * 1e3, api.D[0], D_py)
        if tg is not None:
            s_gs = np.array([tg[:, 0, i], tg[:, 1, i], tg[:, 2, i], tg[:, 3, i]]).T
            f_gs, D_gs = D_of(s_gs)
            dx = np.abs(s_gs - s_py).max()
            worst_x = max(worst_x, dx)
            worst_D = max(worst_D, abs(D_gs - api.D[0]))
            if api.D[0] < -3.5:
                worst_x_reg = max(worst_x_reg, dx)
                worst_D_reg = max(worst_D_reg, abs(D_gs - api.D[0]))
            line += ' %9.4f | %.10f %.10f %.1e | 座標の最大の差 %.1e' % (D_gs, api.nux1[0], f_gs[0], abs(api.nux1[0] - f_gs[0]), dx)
        print(line)
    # カオス的な点: 初期の x を 1e-13・1e-12・1e-11 m ずらした pysad の追跡との差 (丸め誤差の増幅の目安)
    base = np.zeros((n, 6))
    base[:, 0] = [p[0] for p in pts]
    base[:, 2] = [p[1] for p in pts]
    r0 = st.track(LOC, base, NTURN, record=True)
    print('初期の x を dx だけずらした pysad の追跡と元の追跡の、512 周の座標 (x px y py) の最大の差 / 落ちた周 (点ごと):')
    for dx in (1e-13, 1e-12, 1e-11):
        pert = base.copy()
        pert[:, 0] += dx
        r1 = st.track(LOC, pert, NTURN, record=True)
        d = np.nanmax(np.abs(r1['rec'][:, :, :4] - r0['rec'][:, :, :4]), axis=(0, 2))
        print('   dx = %.0e: 最大の差 %s、落ちた周 %s (元 %s)' % (dx, ' '.join('%.1e' % v for v in d), [int(v) for v in r1['lost_turn']], [int(v) for v in r0['lost_turn']]))
    if tg is not None:
        print('座標列の最大の差 (gs と pysad の台本、残った点、全周・x px y py): %.2e' % worst_x)
        print('規則的な点 (D < -3.5) の座標列の最大の差: %.2e、|D(gs の座標) − D(frequency_map)| の最大 = %.2e' % (worst_x_reg, worst_D_reg))
        print('カオス的な点 (D > -3.5 または失われる点) を含めた |D(gs の座標) − D(frequency_map)| の最大 = %.2e' % worst_D)
        print('失われる点: 落ちた周が API・台本・gs で一致 %d / %d (一致しない点は、初期値を 1e-13 m ずらすだけで落ちた周が変わるカオス的な点)' % (nlost_match, nlost_all))


if __name__ == '__main__':
    to = 1800
    if '--timeout' in sys.argv:
        to = int(sys.argv[sys.argv.index('--timeout') + 1])
    main(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else 'kekls', '--nogs' not in sys.argv, to)
