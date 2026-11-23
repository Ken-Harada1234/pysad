# -*- coding: utf-8 -*-
"""Touschek 寿命の式 (pysad = SAD K64 の TouschekTable と、elegant の touschekLifetime = Piwinski の式) の違いを切り分ける。
垂直・水平の運動量の異方性 (結合の大きさ) と、電流 (IBS) を変えて、一様な受容 δ での寿命の比 (elegant / pysad) を出す。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_touschek_scan.py astrid2|kekls [pbunch] [mincoup ...]
 - pbunch (既定 1e6): IBS の影響を小さくするための小さな電流 (寿命は電流に反比例するので、結果は 1e9 粒子/バンチ相当に換算して出す)
 - mincoup: SAD の MINCOUP (垂直エミッタンスの下限 = (εx+εy)·MINCOUP)。大きくすると垂直・水平が等方に近づく。
出力: 標準出力 (elegant_touschek_scan_<名前>.out.txt に保存)。elegant 側の入力は pysad の収束値 (EMITXC・EMITYC・SIGE・SIGZ)。
"""
import json
import os
import sys
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
ENV = dict(os.environ, MSYS_NO_PATHCONV='1')


def wsl(cmd, timeout=3600):
    return subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=timeout,
                          env=ENV, encoding='utf-8', errors='replace')


def tl_elegant(name, apfile, ref, length_mm):
    cmd = ('source ~/elegant/env.sh; cd ~/t1/%s_tl; touschekLifetime scan.out -twiss=tw.twi -aperture=%s -particles=%.6g -emitInput=%.10e '
           '-emityInput=%.10e -deltaInput=%.10e -length=%.10g -ignoreMismatch=silently >/dev/null 2>&1; sddsprintout scan.out -par=tLifetime | grep tLifetime') % (
        name, apfile, ref['PBUNCH'], ref['EMITXC'], ref['EMITYC'], ref['SIGE'], length_mm)
    out = wsl(cmd).stdout
    return float(out.split('=')[1].split()[0]) * 3600.0


def main(name, pb, mincoups):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.touschek import touschek_lifetime
    import sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    rows = []
    amps = {'astrid2': [0.004, 0.006, 0.010, 0.0161756586929], 'kekls': [0.004, 0.006, 0.010, 0.0299397749]}[name]
    for mc in mincoups:
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)
        g = ffs.I.globals
        ffs.flags['INTRA'] = True
        g['PBUNCH'] = pb
        g['MINCOUP'] = mc
        ffs.flags['RFSW'] = True
        ffs.flags['RAD'] = True
        ffs.run_emit(quiet=True)
        ref = dict(PBUNCH=pb, EMITXC=float(g['EMITXC']), EMITYC=float(g['EMITYC']), SIGE=float(g['SIGE']), SIGZ=float(g['SIGZ']))
        tab = ffs.touschek['table']
        print('MINCOUP=%g PBUNCH=%g: EMITXC=%.4e EMITYC=%.4e (ratio %.4f) SIGE=%.4e SIGZ=%.4f mm' % (
            mc, pb, ref['EMITXC'], ref['EMITYC'], ref['EMITYC'] / ref['EMITXC'], ref['SIGE'], ref['SIGZ'] * 1e3), flush=True)
        for a in amps:
            tp = touschek_lifetime(tab, None, None, a / ref['SIGE'], ref['SIGE'])
            fn = '%s.aper_uniform_%.4f.sdds' % (name, a if a < 0.0165 or name == 'kekls' else a)
            # 一様な aperture ファイルがなければ作る
            import elegant_touschek as et
            import numpy as np
            fn = 'scan_u%.6f.sdds' % a
            st0 = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)[0]
            s_el = np.cumsum([0.0] + [float(e.L) for e in st0.els])
            et.write_sdds(os.path.join(HERE, '..', '..', '..', 'local', 'scratch_T1', fn), [(j, s_el[j], a, -a) for j in range(1, len(st0.els) + 1)])
            wsl('source ~/elegant/env.sh; cp "%s" ~/t1/%s_tl/' % (et.wsl_path(os.path.join(HERE, '..', '..', '..', 'local', 'scratch_T1', fn)), name))
            te = tl_elegant(name, fn, ref, ref['SIGZ'] * 1e3)
            print('   a=%.5f  pysad %.4f s  elegant %.4f s  ratio elegant/pysad %.4f   (1e9 粒子/バンチ換算: %.2f / %.2f s)' % (
                a, tp, te, te / tp, tp * pb / 1e9, te * pb / 1e9), flush=True)
            rows.append((mc, a, tp, te))
    return rows


if __name__ == '__main__':
    name = sys.argv[1]
    pb = float(sys.argv[2]) if len(sys.argv) > 2 else 1e6
    mcs = [float(x) for x in sys.argv[3:]] or [0.01, 0.1, 1.0]
    main(name, pb, mcs)
