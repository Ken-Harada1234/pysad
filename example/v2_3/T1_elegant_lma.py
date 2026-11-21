# -*- coding: utf-8 -*-
"""v2.3.0 T1: elegant の momentum_aperture (追跡、RF・放射つき) で局所運動量アクセプタンスを求める (WSL の Pelegant)。

使い方 (リポジトリ直下から。PYTHONUTF8=1 を付ける):
  PYTHONUTF8=1 py -3.12 example/v2_3/T1_elegant_lma.py astrid2 --mpi [--np 4]   ASTRID2: 全周 (199 要素の出口) を Pelegant 1 回で計算
  PYTHONUTF8=1 py -3.12 example/v2_3/T1_elegant_lma.py kekls [--np 4]        kek-ls_v03_12.lat: 1 セルの位置 (pysad の位置の 4 つおき) + 9 番目のセルの 3 点 (1500 周。下記)
  PYTHONUTF8=1 py -3.12 example/v2_3/T1_elegant_lma.py astrid2 --mpi --exact-k1 --cell   K1≠0 の BEND に K2 = −h K1 を足さない (厳密な曲線座標の勾配)、1 セルだけ (感度の確認)
出力: example/v2_3/T1_elegant_<名前>.lma.txt (loc、s、δ±、損失した周・位置、条件)、T1_<名前>_ma.lte / T1_<名前>_ma.ele (入力)。
elegant は WSL (~/elegant/env.sh)。作業は WSL の ~/t1/<名前>/ で行い、終了後に mpirun・Pelegant の残りを確かめる。

条件 (pysad の LMA と揃える): アパーチャ = 各要素の出口の MAXAMP (矩形 ±10 mm。pysad の check_every = 1 と同じ判定の位置)、
x_initial = y_initial = 0 (閉軌道から。elegant は fiducial 粒子の閉軌道)、周回数 = シンクロトロン振動 2 周期分 (pysad と同じ)、
δ の走査: 刻み delta_step_size で外向きに探し、最初の損失の手前から splits 回 split_step_divisor で割って刻みを細かくする (最終の刻み = step/divisor^splits)。
RF バケットから出た粒子は elegant の momentum_aperture では「損失」ではない (アパーチャに当たるまで追跡が続く)。pysad は |z| がバケットを出た時点で損失とする。
シンクロトロン放射の損失が小さいリング (KEK-LS: 1 周 1e-4) では、バケットから出た δ>0 の粒子が 491 周のうちにアパーチャに届かず、δ+ の限界が見つからない (走査の上限になる)。
kek-ls は既定で 1500 周 (シンクロトロン振動 約 6 周期、pysad は 491 周) にして、最初のセルの位置を 4 つおき (11 点) + 9 番目のセルの 3 点にした (--turns N・--stride k で変える)。
さらに、限界が見つからないとき (走査の上限まで落ちない) の elegant は、soft_failure = 1 でも探索が終わらない (上限を超えて δ を増やし続ける) ことがあった。
ここでは周回数を増やして、スリップした粒子が必ず落ちるようにした。--stride k は位置を k 個おきにする。
"""
import json
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
ENV = dict(os.environ, MSYS_NO_PATHCONV='1')

CASES = {
    'astrid2': dict(pmev=580.0, aper=0.010, step=0.0025, splits=5, div=2, limit=0.03, sel='all'),
    'kekls': dict(pmev=3000.0, aper=0.010, step=0.005, splits=6, div=2, limit=0.05, sel='cell', turns=1500, stride=4, cells=(9,)),
}


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def wsl(cmd, timeout=7200):
    """WSL の bash でコマンドを走らせる。wsl.exe の外側のシェルが $ を展開するので、コマンドをスクリプトファイルに書いて bash に渡す"""
    d = os.path.join(ROOT, 'local', 'scratch_T1')
    os.makedirs(d, exist_ok=True)
    sh = os.path.join(d, '_cmd_%d.sh' % os.getpid())
    with open(sh, 'w', newline=chr(10), encoding='utf-8') as f:
        f.write(cmd + chr(10))
    return subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', wsl_path(sh)], capture_output=True, text=True, timeout=timeout,
                          env=ENV, encoding='utf-8', errors='replace')


def leftovers():
    p = wsl("ps -eo pid,etime,pcpu,cmd | grep -E 'Pelegant|elegant|mpirun' | grep -v grep")
    return p.stdout.strip()


def ele_text(name, lte, pmev, nturn, step, splits, div, limit, tag, pattern='AM'):
    return ("&run_setup\n  lattice = %s, use_beamline = RING,\n  p_central_mev = %.6f, default_order = 3,\n  rootname = %s,\n&end\n"
            "&run_control n_passes = %d &end\n"
            "&momentum_aperture\n  output = %%s.mmap, output_mode = 1, fiducialize = 1, verbosity = 1,\n"
            "  x_initial = 0, y_initial = 0,\n"
            "  delta_negative_start = 0, delta_negative_limit = %.6f,\n  delta_positive_start = 0, delta_positive_limit = %.6f,\n"
            "  delta_step_size = %.8g, splits = %d, split_step_divisor = %d,\n"
            "  include_name_pattern = \"%s\",\n&end\n") % (lte, pmev, tag, nturn, -limit, limit, step, splits, div, pattern)


def parse(txt):
    rows = []
    for ln in txt.splitlines():
        t = ln.split()
        if len(t) < 10:
            continue
        try:
            rows.append(dict(s=float(t[0]), occ=int(t[1]), dir=int(t[2]), delta=float(t[3]), found=int(t[4]), lost_pass=int(t[5]),
                             s_lost=float(t[6]), x_lost=float(t[7]), y_lost=float(t[8]), delta_lost=float(t[9])))
        except ValueError:
            pass
    return rows


def lma_locs(name, st, per, stride=1):
    """比べる位置 (pysad の loc。要素 loc-1 の出口 = loc の入口。loc = 1 はリングの終わり (= s 0) の要素 nlat-1 の出口)。
    astrid2: 全周。kekls: 最初のセルの 0.5 m おきの位置 (pysad の Q1 と同じ select_locations) + 別のセル (3・9・15 番目) の同じ位置の 3 点ずつ"""
    from pysad.trksetup import select_locations
    cs = CASES[name]
    if cs['sel'] == 'all':
        locs = list(range(1, st.nlat))
    else:
        ls, le = per['starts'][0], per['starts'][1] - 1
        locs = select_locations(st, ls, le, 1, 0.5)[::stride]
        pick = [locs[len(locs) // 5], locs[len(locs) // 2], locs[(4 * len(locs)) // 5]]
        for k in cs.get('cells', (3, 9, 15)):
            locs = list(locs) + [per['starts'][k - 1] + (l - per['starts'][0]) for l in pick]
    return sorted(set(int(l) for l in locs))


def main(name, nproc=4, exact=False, cell=False, collect=False, mode='serial', turns=None, stride=1):
    cs = CASES[name]
    sfx = ('_exactk1' if exact else '') + ('_n%d' % turns if turns else '') + ('_s%d' % stride if stride > 1 else '')   # 既定の turns・stride (CASES) では付けない
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity, select_locations
    import T1_sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    st = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=cs['aper'])[0]
    cond = st.conditions()
    per = find_periodicity(st.els, st.rows)
    nturn = turns or cs.get('turns') or st.nturn_for_synchrotron(2.0)
    stride = stride if stride > 1 else cs.get('stride', 1)
    locs = lma_locs(name, st, per, stride)
    if cell:
        locs = [l for l in locs if per['starts'][0] <= l <= per['starts'][1] - 1]
    nel = len(st.els)
    ap_sel = set(((l - 2) % nel) + 1 for l in locs)          # loc l の入口 = 要素 l-1 の出口 (loc = 1 は要素 nel)
    path, _, _ = c2e.convert(name, HERE, cond['U0_eV'], aper=cs['aper'], ap_sel=ap_sel, tag='_ma' + sfx, sad_k1=not exact, am_unique=(mode == 'serial'))
    lte = os.path.basename(path)
    tag = 'T1_%s_ma%s' % (name, sfx)
    jlist = sorted(ap_sel)                                  # 要素の番号 (昇順) = AM の ElementOccurence の順
    occ2loc = {i + 1: ((j % nel) + 1) for i, j in enumerate(jlist)}
    wd = '~/t1/%s_ma%s' % (name, sfx)
    sm = {}
    t0 = time.time()
    if mode == 'serial':
        # 位置ごとに独立な serial の elegant (include_name_pattern = AM<番号> で 1 位置) を nproc 本まで同時に走らせる (xargs -P)。
        # Pelegant の momentum_aperture は、kek-ls の位置の塊で数倍遅くなった (一部のプロセスが終わらない) ので、serial を既定にした
        scratch = os.path.join(ROOT, 'local', 'scratch_T1', 'ma_' + tag)
        os.makedirs(scratch, exist_ok=True)
        for f in os.listdir(scratch):
            os.remove(os.path.join(scratch, f))
        for k in range(len(jlist)):
            ptag = '%s_p%03d' % (tag, k)
            open(os.path.join(scratch, ptag + '.ele'), 'w', newline=chr(10)).write(
                ele_text(name, lte, cs['pmev'], nturn, cs['step'], cs['splits'], cs['div'], cs['limit'], ptag, pattern='AM%03d' % (k + 1)).replace('output_mode = 1', 'output_mode = 0, soft_failure = 1'))
        if not collect:
            cmd = ('source ~/elegant/env.sh; mkdir -p %s; cd %s; rm -f %s_p*; cp "%s" .; cp "%s"/*.ele .; '
                   "ls %s_p*.ele | xargs -P %d -I{} bash -c 'timeout -k 20 1800 elegant {} > {}.log 2>&1'; echo done" % (
                       wd, wd, tag, wsl_path(os.path.join(HERE, lte)), wsl_path(scratch), tag, nproc))
            wsl(cmd, timeout=14400)
        cmd = ('source ~/elegant/env.sh; cd %s; for f in %s_p*.mmap; do sdds2stream -col=ElementName,s,deltaPositive,deltaNegative,'
               'lostOnPassPositive,lostOnPassNegative,sLostPositive,sLostNegative,xLostPositive,xLostNegative,deltaPositiveFound,deltaNegativeFound $f; done' % (wd, tag))
        p = wsl(cmd)
        print('elegant (serial、同時 %d 本): %.0f s' % (nproc, time.time() - t0))
        print('終了後の elegant・mpirun の残り:', leftovers() or 'なし')
        n = 0
        # 各ファイルは 1 位置 = 1 行 (ElementName = AM001、AM002、… が jlist の順)
        for ln in p.stdout.splitlines():
            t = ln.split()
            if len(t) < 12 or not t[0].startswith('AM'):
                continue
            loc = occ2loc.get(int(t[0][2:]))
            sm[loc] = dict(loc=loc, s=float(t[1]), dp_p=float(t[2]), dp_m=float(t[3]), pass_p=int(t[4]), pass_m=int(t[5]), slost_p=float(t[6]),
                           slost_m=float(t[7]), xlost_p=float(t[8]), xlost_m=float(t[9]), found_p=int(t[10]), found_m=int(t[11]))
    else:
        # Pelegant の 1 回の実行 (include_name_pattern = AM で全部の位置)。ASTRID2 の全周はこれで求めた (output_mode = 1: 方向ごとのページ)
        text = ele_text(name, lte, cs['pmev'], nturn, cs['step'], cs['splits'], cs['div'], cs['limit'], tag)
        open(os.path.join(ROOT, 'local', 'scratch_T1', tag + '.ele'), 'w', newline=chr(10)).write(text)
        cp = wsl_path(os.path.join(ROOT, 'local', 'scratch_T1', tag + '.ele'))
        cmd = ('source ~/elegant/env.sh; mkdir -p %s; cd %s; cp "%s" "%s" .; rm -f %s.mmap; '
               'timeout -k 20 3000 mpirun -np %d Pelegant %s.ele > %s.log 2>&1; '
               'sdds2stream -col=s,ElementOccurence,direction,delta,deltaFound,lostOnPass,sLost,xLost,yLost,deltaLost %s.mmap' % (
                   wd, wd, wsl_path(os.path.join(HERE, lte)), cp, tag, nproc, tag, tag, tag))
        if collect:
            cmd = 'source ~/elegant/env.sh; cd %s; sdds2stream -col=s,ElementOccurence,direction,delta,deltaFound,lostOnPass,sLost,xLost,yLost,deltaLost %s.mmap' % (wd, tag)
        p = wsl(cmd)
        for r in parse(p.stdout):
            loc = occ2loc.get(r['occ'])
            d = sm.setdefault(loc, dict(loc=loc, s=r['s']))
            k = 'p' if r['dir'] > 0 else 'm'
            d['dp_' + k] = r['delta']
            d['found_' + k] = r['found']
            d['pass_' + k] = r['lost_pass']
            d['slost_' + k] = r['s_lost']
            d['xlost_' + k] = r['x_lost']
        print('elegant (mpirun -np %d): %.0f s' % (nproc, time.time() - t0))
        print('終了後の elegant・mpirun の残り:', leftovers() or 'なし')
    print('位置', len(sm), '/', len(locs), '、欠けた位置:', sorted(set(locs) - set(sm)) or 'なし')
    p = type('P', (), {'stdout': '', 'stderr': ''})()
    out = os.path.join(HERE, 'T1_elegant_%s%s.lma.txt' % (name, sfx))
    with open(out, 'w', encoding='utf-8', newline='\n') as f:
        f.write('# elegant momentum_aperture (Pelegant -np %d) %s: p = %.1f MeV/c, aperture +-%.1f mm (MAXAMP at every element exit), nturn=%d, '
                'x_initial=y_initial=0, RF+SYNCH_RAD (ISR off)\n' % (nproc, name, cs['pmev'], cs['aper'] * 1e3, nturn))
        f.write('# delta_step_size=%g splits=%d split_step_divisor=%d limit=+-%g (final step %.3g)\n' % (
            cs['step'], cs['splits'], cs['div'], cs['limit'], cs['step'] / cs['div'] ** cs['splits']))
        f.write('# loc  s[m]  dp+  dp-  lostpass+  lostpass-  slost+  slost-  found+ found-\n')
        for loc in sorted(sm):
            d = sm[loc]
            f.write('%5d %10.4f %10.6f %10.6f %6d %6d %9.3f %9.3f %d %d\n' % (
                loc, d['s'], d['dp_p'], d['dp_m'], d['pass_p'], d['pass_m'], d['slost_p'], d['slost_m'], d['found_p'], d['found_m']))
    print('wrote', out)
    open(os.path.join(HERE, 'T1_elegant_%s%s.log.txt' % (name, sfx)), 'w').write(p.stdout[-200:] + '\n' + p.stderr[-2000:])


if __name__ == '__main__':
    nproc = 4
    if '--np' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('--np') + 1])
    main(sys.argv[1], nproc, exact='--exact-k1' in sys.argv, cell='--cell' in sys.argv, collect='--collect' in sys.argv, mode='mpi' if '--mpi' in sys.argv else 'serial',
         turns=int(sys.argv[sys.argv.index('--turns') + 1]) if '--turns' in sys.argv else None,
         stride=int(sys.argv[sys.argv.index('--stride') + 1]) if '--stride' in sys.argv else 1)
