# -*- coding: utf-8 -*-
"""v2.3.0 T2: Touschek 寿命の 3 通りの比較 (pysad = K64 の式 / elegant touschekLifetime / 独立な式 touschek_piwinski)。
同じビーム (pysad の IBS 込みの収束値)・同じ光学・同じ一様な受容。結合 (MINCOUP)・電流 (PBUNCH)・受容を変えて、寿命の表を出す。

使い方 (リポジトリ直下から。PYTHONUTF8=1):
  PYTHONUTF8=1 py -3.12 example/v2_3/T2_scan.py astrid2|kekls [--sets "1e9:0.01,1e9:0.1,..."] [--amps 0.004,0.006,0.010,0.016]
  --sets は PBUNCH:MINCOUP の組。出力は標準出力と T2_scan_<名前>.json。
elegant 側 (WSL、~/t2/<名前>): 光学は elegant の twiss (radiation_integrals=1) を、全要素を 8 分割 (divide_elements) して標本化を細かくする。
 (elegant の touschekLifetime は要素の出口だけで局所の損失率を標本化し台形則で積分するので、粗い要素のリングでは積分の離散化誤差が出る (§3.2))。
 入力 (T1 の lte `T1_<名前>_all.lte` を使う): -emitInput = EMITXC、-emityInput = EMITYC、-deltaInput = SIGE、-length = SIGZ、-particles = PBUNCH。
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
ENV = dict(os.environ, MSYS_NO_PATHCONV='1')
PMEV = {'astrid2': 580.0, 'kekls': 3000.0}


def wsl(cmd, timeout=3600):
    return subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=timeout,
                          env=ENV, encoding='utf-8', errors='replace')


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def gather(name, pb, mc):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.touschek import touschek_table
    import T1_sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)
    g = ffs.I.globals
    ffs.flags['INTRA'] = True
    g['PBUNCH'] = pb
    g['MINCOUP'] = mc
    ffs.flags['RFSW'] = True
    ffs.flags['RAD'] = True
    ffs.run_emit(quiet=True)
    tt = ffs.touschek
    tt2 = touschek_table(ffs.lat, ffs.momentum(), ffs.emit, float(pb), mincoup=float(mc), rfsw=True, radcod=ffs.flags['RADCOD'],
                         verbose=False, intra=True, record_rest=True)
    ref = {k: float(g[k]) for k in ('EMITXC', 'EMITYC', 'EMITZC', 'SIGE', 'SIGZ')}
    return tt, tt2, ref


def elegant_setup(name):
    """全要素を 8 分割した elegant の twiss を ~/t2/<名前>/ に作る (T1 の lte を使う)"""
    wd = '~/t2/%s' % name
    lte = 'T1_%s_all.lte' % name
    ele = ("&divide_elements name=*, type=*, divisions=8, exclude=\"MAXAMP\" &end\n&run_setup\n  lattice = %s, use_beamline = RING, p_central_mev = %.6f, "
           "default_order = 3, rootname = twd, magnets = %%s.mag,\n&end\n&run_control &end\n&twiss_output filename = %%s.twi, radiation_integrals = 1, &end\n"
           "&bunched_beam n_particles_per_bunch = 1 &end\n&track &end\n") % (lte, PMEV[name])
    p = os.path.join(HERE, 'T2_%s_div.ele' % name)
    open(p, 'w', newline='\n').write(ele)
    cmd = ('source ~/elegant/env.sh; mkdir -p %s; cd %s; cp "%s" "%s" .; timeout 900 elegant T2_%s_div.ele > d.log 2>&1; ls twd.twi') % (
        wd, wd, wsl_path(os.path.join(HERE, lte)), wsl_path(p), name)
    return wsl(cmd).stdout.strip()


def elegant_tl(name, a, ref, pb, tag):
    from T1_elegant_touschek import write_sdds
    wd = '~/t2/%s' % name
    fn = os.path.join(ROOT, 'local', 'scratch_T2', 'ap_%s_%.6f.sdds' % (name, a))
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    write_sdds(fn, [(1, 0.0, a, -a)])
    out = 'tl_%s.out' % tag
    cmd = ('source ~/elegant/env.sh; cd %s; cp "%s" ap.sdds; timeout 1800 touschekLifetime %s -twiss=twd.twi -aperture=ap.sdds -particles=%.6g '
           '-emitInput=%.10e -emityInput=%.10e -deltaInput=%.10e -length=%.10g -ignoreMismatch=silently > /dev/null 2>&1; '
           'sdds2stream %s -par=tLifetime; sdds2stream %s -col=tmP | sort -u | head -3') % (wd, wsl_path(fn), out, pb, ref['EMITXC'], ref['EMITYC'], ref['SIGE'],
                                                                                       ref['SIGZ'] * 1e3, out, out)
    r = wsl(cmd).stdout.split()
    return float(r[0]) * 3600.0, [float(x) for x in r[1:]]


def main():
    name = sys.argv[1]
    sets = [(1e9, 0.01)]
    amps = [0.004, 0.006, 0.010, 0.016] if name == 'astrid2' else [0.004, 0.006, 0.010, 0.020]
    if '--sets' in sys.argv:
        sets = [tuple(float(x) for x in s.split(':')) for s in sys.argv[sys.argv.index('--sets') + 1].split(',')]
    if '--amps' in sys.argv:
        amps = [float(x) for x in sys.argv[sys.argv.index('--amps') + 1].split(',')]
    from pysad.touschek import touschek_lifetime, touschek_piwinski
    print('elegant 光学:', elegant_setup(name), flush=True)
    res = []
    for pb, mc in sets:
        tt, tt2, ref = gather(name, pb, mc)
        tab = tt['table']
        print('\nPBUNCH=%g MINCOUP=%g: εx=%.4e εy=%.4e (比 %.4f) σδ=%.4e σz=%.3f mm' % (pb, mc, ref['EMITXC'], ref['EMITYC'], ref['EMITYC'] / ref['EMITXC'],
                                                                                   ref['SIGE'], ref['SIGZ'] * 1e3), flush=True)
        print('   a      pysad(K64) [s]  elegant [s]  touschek_piwinski [s]   K64/独立-1   elegant/独立-1')
        for a in amps:
            k64 = touschek_lifetime(tab, None, None, a / ref['SIGE'], ref['SIGE'])
            piw = touschek_piwinski(tt2, a=a)['lifetime']
            ele, tm = elegant_tl(name, a, ref, pb, 'u%.4f' % a)
            print('  %.4f  %14.4f  %12.4f  %14.4f   %+8.4f   %+8.4f' % (a, k64, ele, piw, k64 / piw - 1, ele / piw - 1), flush=True)
            res.append(dict(pbunch=pb, mincoup=mc, a=a, k64=k64, elegant=ele, piwinski=piw, ref=ref))
    json.dump(res, open(os.path.join(HERE, 'T2_scan_%s.json' % name), 'w'), indent=1)


if __name__ == '__main__':
    main()
