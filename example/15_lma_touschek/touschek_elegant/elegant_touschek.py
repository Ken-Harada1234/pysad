# -*- coding: utf-8 -*-
"""elegant の touschekLifetime (Piwinski の式。Xiao・Borland の積分) による Touschek 寿命。pysad と同じ入力にして比べる。

使い方 (リポジトリ直下から。PYTHONUTF8=1 を付ける。先に pysad_ref.py で pysad 側の値 pysad_<名前>.ref.json を作る):
  PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_touschek.py astrid2|kekls
出力: elegant_<名前>.touschek.txt (寿命の表)、入力の SDDS (<名前>.aper_*.sdds、ascii の小さなファイル)。
elegant 側の入力: Twiss (radiation_integrals=1。全要素の出口 = MAXAMP 'AM' の位置で標本化)、運動量アクセプタンスのファイル、
  -particles = PBUNCH、-emitInput = 水平エミッタンス (pysad の IBS 込みの収束値 EMITXC)、-emityInput = 垂直エミッタンス (EMITYC)、
  -deltaInput = エネルギー拡がり (IBS 込み)、-length = バンチ長 (IBS 込み)。
受容: (a) 一様な δ (±0.5%、±1.0%、±RF バケットの高さなど)、(b) elegant の momentum_aperture の δ±(s) (ASTRID2 は全周、kek-ls は 1 セルを周期的に延ばす)、
  (c) pysad の LMA の δ±(s) (同じ s の表を elegant の touschekLifetime に渡す: 寿命の式だけの違いを見る)。
"""
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
ENV = dict(os.environ, MSYS_NO_PATHCONV='1')


def wsl_path(p):
    p = os.path.abspath(p).replace(os.sep, '/')
    return '/mnt/' + p[0].lower() + p[2:]


def wsl(cmd, timeout=7200):
    return subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=timeout,
                          env=ENV, encoding='utf-8', errors='replace')


def write_sdds(path, rows):
    """rows: [(occ, s, dpp, dpm)] を ElementName=AM の aperture ファイル (ascii SDDS) に"""
    with open(path, 'w', newline='\n') as f:
        f.write('SDDS1\n&description text="momentum aperture (T1)", contents="momentum aperture" &end\n')
        f.write('&column name=ElementName, type=string &end\n&column name=s, type=double, units=m &end\n')
        f.write('&column name=ElementType, type=string &end\n&column name=ElementOccurence, type=long &end\n')
        f.write('&column name=deltaPositive, type=double &end\n&column name=deltaNegative, type=double &end\n')
        f.write('&data mode=ascii &end\n! page number 1\n%d\n' % len(rows))
        for occ, s, dpp, dpm in rows:
            f.write('AM %.9f MAXAMP %d %.9e %.9e\n' % (s, occ, dpp, dpm))


def run_tl(name, aper, ref, length_mm, extra='', tag='x'):
    wd = '~/t1/%s_tl' % name
    out = 'tl_%s.out' % tag
    cmd = ('source ~/elegant/env.sh; cd %s; touschekLifetime %s -twiss=tw.twi -aperture=%s -particles=%g -emitInput=%.10e -emityInput=%.10e '
           '-deltaInput=%.10e -length=%.10g -ignoreMismatch=silently %s 2>&1 | tail -3; '
           'sddsprintout %s -par=* 2>&1 | grep -i -E "lifetime|Tlife" ') % (wd, out, aper, ref['PBUNCH'], ref['EMITXC'], ref['EMITYC'], ref['SIGE'], length_mm, extra, out)
    p = wsl(cmd)
    return p.stdout


def main(name):
    ref = json.load(open(os.path.join(HERE, 'pysad_%s.ref.json' % name)))
    import sad2ele as c2e
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    wd = '~/t1/%s_tl' % name
    # 全要素の出口を AM にした lte (Twiss の標本化)
    st = TrackSetup.from_source(open(c2e.CASES[name]['path'], encoding='utf-8').read(), 'RING', rad=True, rfsw=True)[0]
    nel = len(st.els)
    path, _, _ = c2e.convert(name, HERE, st.conditions()['U0_eV'], ap_sel=set(range(1, nel + 1)), tag='_all')
    pmev = {'astrid2': 580.0, 'kekls': 3000.0}[name]
    open(os.path.join(HERE, '%s_all.twiss.ele' % name), 'w', newline='\n').write(
        "&run_setup\n  lattice = %s, use_beamline = RING, p_central_mev = %.6f, default_order = 3, rootname = tw, magnets = %%s.mag,\n&end\n"
        "&run_control &end\n&twiss_output filename = %%s.twi, radiation_integrals = 1, &end\n"
        "&bunched_beam n_particles_per_bunch = 1 &end\n&track &end\n" % (os.path.basename(path), pmev))
    p = wsl('source ~/elegant/env.sh; mkdir -p %s; cd %s; cp "%s" "%s" .; elegant %s_all.twiss.ele > tw.log 2>&1; ls tw.twi' % (
        wd, wd, wsl_path(path), wsl_path(os.path.join(HERE, '%s_all.twiss.ele' % name)), name))
    print(p.stdout.strip())
    s_el = np.cumsum([0.0] + [float(e.L) for e in st.els])        # 要素 j (1 始まり) の出口の s = s_el[j]
    circ = s_el[-1]
    # elegant の LMA (cell または全周) → 全要素の出口の δ±
    lma_path = os.path.join(HERE, 'elegant_%s.lma.txt' % name)
    have_ele = os.path.exists(lma_path)
    lma = np.loadtxt(lma_path, comments='#') if have_ele else np.zeros((2, 4))
    L = float(np.atleast_1d(ref.get('cell_length', circ))[0])
    # loc l の s は要素 l-1 の出口 (AM = 要素 j = l-1)
    s_l, dpp_l, dpm_l = lma[:, 1], lma[:, 2], lma[:, 3]
    order = np.argsort(s_l)
    s_l, dpp_l, dpm_l = s_l[order], dpp_l[order], dpm_l[order]
    res = {}

    def aper_rows(kind, arg=None):
        rows = []
        for j in range(1, nel + 1):
            s = s_el[j]
            if kind == 'uniform':
                dpp, dpm = arg, -arg
            elif kind == 'elegant':
                sm = s % L
                if s_l.max() - s_l.min() > 0.9 * circ:              # 全周 (ASTRID2)
                    dpp = np.interp(s, s_l, dpp_l, period=circ)
                    dpm = np.interp(s, s_l, dpm_l, period=circ)
                else:                                                  # 1 セルを周期的に延ばす
                    cm = s_l < L + 1e-9                                      # 最初のセルの位置だけ
                    dpp = np.interp(sm, s_l[cm], dpp_l[cm])
                    dpm = np.interp(sm, s_l[cm], dpm_l[cm])
            elif kind == 'pysad':
                tab = np.loadtxt(os.path.join(HERE, arg), comments='#', usecols=(0, 1, 3, 4))
                sm = s % L
                dpp = np.interp(sm, tab[:, 1], tab[:, 2])
                dpm = np.interp(sm, tab[:, 1], tab[:, 3])
            rows.append((j, s, float(dpp), float(dpm)))
        return rows

    lines = []
    length_mm = ref['SIGZ'] * 1e3
    p = wsl('source ~/elegant/env.sh; cd %s; cp tw.twi tw.twi 2>/dev/null; ls' % wd)
    cases = [('uniform %.4f' % a, 'uniform', a) for a in ref['uniform_a']]
    if have_ele:
        cases += [('elegant LMA', 'elegant', None)]
    if os.path.exists(os.path.join(HERE, ref.get('pysad_table', 'none'))):
        cases += [('pysad LMA', 'pysad', ref['pysad_table'])]
    for label, kind, arg in cases:
        fn = '%s.aper_%s.sdds' % (name, label.replace(' ', '_'))
        write_sdds(os.path.join(HERE, fn), aper_rows(kind, arg))
        tag = label.replace(' ', '_')
        wsl('source ~/elegant/env.sh; cd %s; cp "%s" .' % (wd, wsl_path(os.path.join(HERE, fn))))
        out = run_tl(name, fn, ref, length_mm, tag=tag)
        print(label, '->', out.strip())
        lines.append((label, out.strip()))
    with open(os.path.join(HERE, 'elegant_%s.touschek.txt' % name), 'w', encoding='utf-8') as f:
        f.write('# elegant touschekLifetime, particles=%g emit=(%.5e, %.5e) sigmaDelta=%.5e sigma_z=%.4f mm\n' % (
            ref['PBUNCH'], ref['EMITXC'], ref['EMITYC'], ref['SIGE'], length_mm))
        for lb, o in lines:
            f.write('%s: %s\n' % (lb, o.replace('\n', ' | ')))


if __name__ == '__main__':
    main(sys.argv[1])
