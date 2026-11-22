# -*- coding: utf-8 -*-
"""example/dxf_lattice の台本を作り (mk)、gs と pysad で走らせ、出力の DXF が (1) gs と pysad、(2) pysad の台本と Python 関数
(pysad.dxf_lattice.write_dxf) で一致することを確かめ、図 (PNG) を作る。
使い方: py -3.12 example/dxf_lattice/run_dxf_examples.py [mk] [run] [png]      (既定は 3 つとも)
環境: PYSAD_SAD_PACKAGES (K64 版 SAD の Packages)。gs は example/tools/sadenv.py に従う。"""
import io
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))

# ---- 例の定義: 台本名 → ラティスのファイル、台本の本体、出力 {ファイル名: (USE するライン, Python 関数の引数)}
FODO_SPEC = 'Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1}, "QF"->{Width->.3, PhysicalLength->.26}, "QUAD"->{Width->.36}}'
FODO_SPEC_PY = {'B': {'Width': .5, 'PoleWidth': .15, 'Type': 'Ctype1'}, 'QF': {'Width': .3, 'PhysicalLength': .26}, 'QUAD': {'Width': .36}}
EXAMPLES = {
    'DXF_FODO': {
        'lattice': os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'),
        'outputs': {
            'DXF_FODO_ring.dxf': ('RING', 'DrawMagnetsDXF["DXF_FODO_ring.dxf"]', {}),
            'DXF_FODO_ring_spec.dxf': ('RING', 'DrawMagnetsDXF["DXF_FODO_ring_spec.dxf", %s]' % FODO_SPEC, {'spec': FODO_SPEC_PY}),
            'DXF_FODO_region_mm.dxf': ('RING', 'DrawMagnetsDXF["DXF_FODO_region_mm.dxf", Region->{"QF.1", "QF.4"}, Scale->1000, TextHeight->.1]',
                                       {'region': ('QF.1', 'QF.4'), 'scale': 1000, 'text_height': .1}),
        }},
    'DXF_SOL': {
        'lattice': os.path.join(ROOT, 'example', '00_lattice', 'sol_ring.sad'),
        'outputs': {
            'DXF_SOL_cells.dxf': ('RINGS', 'DrawMagnetsDXF["DXF_SOL_cells.dxf", Region->{1, 25}, TextHeight->.1]',
                                  {'region': (1, 25), 'text_height': .1}),
        }},
    'DXF_ASTRID2': {
        'lattice': os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'),
        'outputs': {
            'DXF_ASTRID2_cell.dxf': ('NCELL', 'DrawMagnetsDXF["DXF_ASTRID2_cell.dxf", TextHeight->.1]', {'text_height': .1}),
            'DXF_ASTRID2_ring.dxf': ('RING', 'DrawMagnetsDXF["DXF_ASTRID2_ring.dxf", TextHeight->.2]', {'text_height': .2}),
        }},
    'DXF_KEKLS': {
        'lattice': os.path.join(ROOT, 'example', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad'),
        'outputs': {
            'DXF_KEKLS_cell.dxf': ('ARCA', 'DrawMagnetsDXF["DXF_KEKLS_cell.dxf", TextHeight->.12]', {'text_height': .12}),
            'DXF_KEKLS_ring.dxf': ('RING', 'DrawMagnetsDXF["DXF_KEKLS_ring.dxf", TextHeight->.3]', {'text_height': .3}),
        }},
}
CLOSE = '''
g1 = LINE["GEO", 1]; g2 = LINE["GEO", LINE["LENGTH"]];
Print["  start (gx, gy, chi1) = ", g1[[1, 1]], " ", g1[[1, 2]], " ", g1[[2, 1]]];
Print["  end   (gx, gy, chi1) = ", g2[[1, 1]], " ", g2[[1, 2]], " ", g2[[2, 1]]];
'''


def footer(ex):
    s = ['', 'FFS;', '$FORM = "14.9";', 'Get["DrawMagnetsDXF.n"];']
    last = None
    for fn, (use, call, _kw) in ex['outputs'].items():
        if use != last:
            s += ['USE %s;' % use, 'Print["USE %s: ", LINE["LENGTH"] - 1, " elements, L = ", LINE["S", LINE["LENGTH"]], " m"];' % use]
            last = use
        s += ['%s;' % call, 'Print["  wrote %s"];' % fn]
        if fn.endswith('_ring.dxf'):
            s += [CLOSE.strip('\n')]
    s += ['Exit;', '']
    return '\n'.join(s)


def mk():
    for name, ex in EXAMPLES.items():
        lat = open(ex['lattice'], encoding='utf-8').read()
        if not lat.endswith('\n'):
            lat += '\n'
        open(os.path.join(HERE, name + '.sad'), 'w', encoding='utf-8', newline='\n').write(lat + footer(ex))
        print('made', name + '.sad')


def run_gs(name, timeout=600):
    import sadenv
    gs = sadenv.gs_wsl()
    tmp = tempfile.gettempdir()
    for fn in EXAMPLES[name]['outputs']:
        if os.path.exists(os.path.join(tmp, fn)):
            os.remove(os.path.join(tmp, fn))
    shutil.copy(os.path.join(HERE, name + '.sad'), os.path.join(tmp, '_dxfex_%s.sad' % name))
    shutil.copy(os.path.join(HERE, 'DrawMagnetsDXF.n'), os.path.join(tmp, 'DrawMagnetsDXF.n'))
    t = tmp.replace(chr(92), '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c',
                        'cd %s && timeout %d %s _dxfex_%s.sad < /dev/null' % (t, timeout, gs, name)], capture_output=True, env=env)
    open(os.path.join(HERE, name + '.gs.txt'), 'wb').write(r.stdout + r.stderr)
    gsdir = os.path.join(tmp, 'dxfex_gs')
    os.makedirs(gsdir, exist_ok=True)
    for fn in EXAMPLES[name]['outputs']:
        shutil.copy(os.path.join(tmp, fn), os.path.join(gsdir, fn))
    return gsdir


def run_pysad(name):
    import sadenv
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8', PYSAD_SAD_PACKAGES=sadenv.packages())
    code = ("import os;os.chdir(%r);from pysad.ffscmd import run_sad_string;run_sad_string(open(%r,encoding='utf-8').read(), cwd=%r)"
            % (HERE, name + '.sad', HERE))
    r = subprocess.run(['py', '-3.12', '-c', code], capture_output=True, env=env)
    open(os.path.join(HERE, name + '.pysad.txt'), 'wb').write(r.stdout + r.stderr)


def lines(path):
    with open(path, encoding='ascii', newline='') as fh:
        return fh.read().replace('\r\n', '\n').split('\n')


def compare(a, b):
    """行ごとの比較。戻り (行数 a, 行数 b, 文字列として違う行の数, 数値として 1e-6 を超えて違う行の数)"""
    la, lb = lines(a), lines(b)
    nd = nn = 0
    for x, y in zip(la, lb):
        if x != y:
            nd += 1
            try:
                if abs(float(x) - float(y)) > 1.5e-6:
                    nn += 1
            except ValueError:
                nn += 1
    return len(la), len(lb), nd, nn


def check(name, gsdir):
    from pysad.dxf_lattice import load_ffs, write_dxf
    ex = EXAMPLES[name]
    ffs_cache = {}
    res = []
    for fn, (use, _call, kw) in ex['outputs'].items():
        if use not in ffs_cache:
            ffs_cache[use] = load_ffs(ex['lattice'], use)
        py = os.path.join(tempfile.gettempdir(), 'dxfex_py_' + fn)
        write_dxf(ffs_cache[use], py, **kw)
        sc = os.path.join(HERE, fn)
        g = os.path.join(gsdir, fn)
        res.append((fn, os.path.getsize(sc), compare(g, sc), compare(py, sc)))
    return res


def png():
    import plot_dxf
    figs = [
        ('DXF_FODO_ring.dxf', 'FODO_ring.png', 'FODO ring (24 cells), plan view', None, False, (8, 8)),
        ('DXF_FODO_ring_spec.dxf', 'FODO_cell_spec.png', 'FODO ring, Spec (B: Ctype1, QF shorter), first cells', (-2, 18, -3, 8), True, (10, 5)),
        ('DXF_FODO_region_mm.dxf', 'FODO_region_mm.png', 'FODO Region {QF.1, QF.4}, Scale 1000 (mm)', None, True, (10, 4)),
        ('DXF_SOL_cells.dxf', 'SOL_cells.png', 'SOL region (SA to SB, BZ = 5 T) in the first cells', (-2, 14, -2, 4), True, (10, 4)),
        ('DXF_ASTRID2_cell.dxf', 'ASTRID2_cell.png', 'ASTRID2 one cell (NCELL)', None, True, (10, 5)),
        ('DXF_ASTRID2_ring.dxf', 'ASTRID2_ring.png', 'ASTRID2 ring', None, False, (8, 8)),
        ('DXF_KEKLS_cell.dxf', 'KEKLS_cell.png', 'KEK-LS one standard cell (ARCA)', None, True, (12, 5)),
        ('DXF_KEKLS_ring.dxf', 'KEKLS_ring.png', 'KEK-LS ring (RING)', None, False, (9, 9)),
    ]
    for dxf, out, title, lim, labels, size in figs:
        plot_dxf.plot(os.path.join(HERE, dxf), os.path.join(HERE, out), title, lim, labels, 130, size)
        print('png', out)


if __name__ == '__main__':
    todo = sys.argv[1:] or ['mk', 'run', 'png']
    if 'mk' in todo:
        mk()
    if 'run' in todo:
        for name in EXAMPLES:
            gsdir = run_gs(name)
            run_pysad(name)
            for r in check(name, gsdir):
                print('%-26s %8d bytes  gs-vs-script (lines, lines, str-diff, num-diff)=%s  python-vs-script=%s' % (r[0], r[1], r[2], r[3]))
    if 'png' in todo:
        png()
