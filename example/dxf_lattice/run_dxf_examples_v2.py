# -*- coding: utf-8 -*-
"""v2 (MagnetSpec[] の読み込み・要素名の重なり消し) の例: KEK-LS の標準セル 1 個。台本 DXF_KEKLS_v2.sad を作り (mk)、gs と pysad で走らせ、
出力の DXF が (1) gs と pysad の台本、(2) pysad の台本と Python 関数 write_dxf (台本の MagnetSpec を自動で読む場合と dict で渡す場合) で
全行一致することを確かめ (run)、図 (PNG、*_v2*.png) を作る (png)。既存の DXF・PNG は上書きしない (出力は名前に _v2 を付けた別ファイルで、
すでにあれば止まる)。
使い方: py -3.12 example/dxf_lattice/run_dxf_examples_v2.py [mk] [run] [png]      (既定は 3 つとも)
環境: PYSAD_SAD_PACKAGES (K64 版 SAD の Packages)。gs は example/tools/sadenv.py に従う (WSL)。"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
sys.path.insert(0, HERE)
from run_dxf_examples import compare   # noqa: E402

NAME = 'DXF_KEKLS_v2'
LATTICE = os.path.join(ROOT, 'example', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad')
# 出力 → (台本の呼び出しのオプション, Python 関数の引数)。TextHeight 0.12 は既存の DXF_KEKLS_cell.dxf と同じ
OUTPUTS = {
    'DXF_KEKLS_cell_v2_remove.dxf': ('TextHeight->.12, UseMagnetSpec->False, Overlap->"remove"',
                                     {'text_height': .12, 'magnet_spec': False, 'overlap': 'remove'}),
    'DXF_KEKLS_cell_v2_spec.dxf': ('TextHeight->.12', {'text_height': .12}),
    'DXF_KEKLS_cell_v2_spec_remove.dxf': ('TextHeight->.12, Overlap->"remove"', {'text_height': .12, 'overlap': 'remove'}),
    # 文字を大きくした (0.3 m) 版: 全部出す版と重なり消しの版
    'DXF_KEKLS_cell_v2_spec_h30.dxf': ('TextHeight->.3', {'text_height': .3}),
    'DXF_KEKLS_cell_v2_spec_remove_h30.dxf': ('TextHeight->.3, Overlap->"remove"', {'text_height': .3, 'overlap': 'remove'}),
}
# MagnetSpec_KEKLS.n と同じ内容を dict で (Python 関数に直接渡す形。台本から読む場合と一致することを確かめる)
MS_DICT = {'BLOW': {'Width': .50, 'Height': .30, 'PoleWidth': .12, 'Type': 'Htype'},
           'BGRAD': {'Width': .56, 'Height': .30, 'PoleWidth': .12, 'Type': 'Ctype1'},
           'QCELL': {'Width': .32, 'Height': .32}, 'SEXT6': {'Width': .24, 'Height': .24}, 'OCFA': {'Width': .20, 'Height': .20}}
_B = ('BID5', 'BID4', 'BID3H', 'BID2', 'BID1', 'BPI1', 'BPI2', 'BPI3H', 'BPI4', 'BPI5')
_Q = ('QFA1H', 'QDA2', 'QDA3', 'QFA4', 'QDA5', 'QFA6H', 'QFA8H')


def mtype(n):
    return 'BLOW' if n in _B else 'BGRAD' if n in ('BPI6H', 'BPI7H') else 'QCELL' if n in _Q else 'SEXT6' if n in ('SDH', 'SFH', 'SHA') else n


def mk():
    lat = open(LATTICE, encoding='utf-8').read()
    if not lat.endswith('\n'):
        lat += '\n'
    s = ['', 'FFS;', '$FORM = "14.9";', 'Get["DrawMagnetsDXF.n"];', 'Get["MagnetSpec_KEKLS.n"];', 'USE ARCA;',
         'Print["USE ARCA: ", LINE["LENGTH"] - 1, " elements, L = ", LINE["S", LINE["LENGTH"]], " m"];']
    for fn, (opts, _kw) in OUTPUTS.items():
        s += ['DrawMagnetsDXF["%s", %s];' % (fn, opts), 'Print["  wrote %s"];' % fn]
    s += ['Exit;', '']
    out = os.path.join(HERE, NAME + '.sad')
    if os.path.exists(out):
        raise SystemExit('%s は既にある (上書きしない)。消してから作り直すこと' % out)
    open(out, 'w', encoding='utf-8', newline='\n').write(lat + '\n'.join(s))
    print('made', NAME + '.sad')


def run_gs(timeout=600):
    import sadenv
    gs = sadenv.gs_wsl()
    tmp = tempfile.gettempdir()
    for fn in OUTPUTS:
        if os.path.exists(os.path.join(tmp, fn)):
            os.remove(os.path.join(tmp, fn))
    shutil.copy(os.path.join(HERE, NAME + '.sad'), os.path.join(tmp, '_dxfex_%s.sad' % NAME))
    for f in ('DrawMagnetsDXF.n', 'MagnetSpec_KEKLS.n'):
        shutil.copy(os.path.join(HERE, f), os.path.join(tmp, f))
    t = tmp.replace(chr(92), '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c',
                        'cd %s && timeout %d %s _dxfex_%s.sad < /dev/null' % (t, timeout, gs, NAME)], capture_output=True, env=env)
    open(os.path.join(HERE, NAME + '.gs.txt'), 'wb').write(r.stdout + r.stderr)
    gsdir = os.path.join(tmp, 'dxfex_gs_v2')
    os.makedirs(gsdir, exist_ok=True)
    for fn in OUTPUTS:
        shutil.copy(os.path.join(tmp, fn), os.path.join(gsdir, fn))
    return gsdir


def run_pysad():
    """pysad で台本を走らせる (出力はこのフォルダ。すでにあれば上書きしないよう止める)"""
    import sadenv
    for fn in OUTPUTS:
        if os.path.exists(os.path.join(HERE, fn)):
            raise SystemExit('%s は既にある (上書きしない)。消してから走らせること' % fn)
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8', PYSAD_SAD_PACKAGES=sadenv.packages())
    code = ("import os;os.chdir(%r);from pysad.ffscmd import run_sad_string;run_sad_string(open(%r,encoding='utf-8').read(), cwd=%r)"
            % (HERE, NAME + '.sad', HERE))
    r = subprocess.run(['py', '-3.12', '-c', code], capture_output=True, env=env)
    open(os.path.join(HERE, NAME + '.pysad.txt'), 'wb').write(r.stdout + r.stderr)


def check(gsdir):
    from pysad.dxf_lattice import load_ffs, read_dxf, write_dxf
    lat = open(LATTICE, encoding='utf-8').read() + '\nFFS;\nGet["MagnetSpec_KEKLS.n"];\n'
    ffs = load_ffs(lat, 'ARCA', cwd=HERE)            # 台本の MagnetSpec を自動で読む
    res = []
    for fn, (_opts, kw) in OUTPUTS.items():
        tmp = tempfile.gettempdir()
        p1 = os.path.join(tmp, 'dxfex_py1_' + fn)
        write_dxf(ffs, p1, **kw)
        kw2 = dict(kw)
        if kw2.get('magnet_spec', None) is None:
            kw2.update(magnet_spec=MS_DICT, magnet_type=mtype)
        p2 = os.path.join(tmp, 'dxfex_py2_' + fn)
        write_dxf(ffs, p2, **kw2)
        sc = os.path.join(HERE, fn)
        nt = sum(1 for e in read_dxf(sc)['entities'] if e['type'] == 'TEXT')
        res.append((fn, os.path.getsize(sc), nt, compare(os.path.join(gsdir, fn), sc), compare(p1, sc), compare(p2, sc)))
    return res


def png():
    import plot_dxf
    figs = [('DXF_KEKLS_cell.dxf', 'KEKLS_cell_v2_default_keep.png', 'KEK-LS ARCA: default widths, all names (existing DXF)'),
            ('DXF_KEKLS_cell_v2_remove.dxf', 'KEKLS_cell_v2_remove.png', 'KEK-LS ARCA: default widths, Overlap="remove"'),
            ('DXF_KEKLS_cell_v2_spec.dxf', 'KEKLS_cell_v2_spec.png', 'KEK-LS ARCA: MagnetSpec widths (illustrative), all names'),
            ('DXF_KEKLS_cell_v2_spec_remove.dxf', 'KEKLS_cell_v2_spec_remove.png', 'KEK-LS ARCA: MagnetSpec widths, Overlap="remove"'),
            ('DXF_KEKLS_cell_v2_spec_h30.dxf', 'KEKLS_cell_v2_spec_h30.png', 'KEK-LS ARCA: MagnetSpec widths, TextHeight 0.3 m, all names'),
            ('DXF_KEKLS_cell_v2_spec_remove_h30.dxf', 'KEKLS_cell_v2_spec_remove_h30.png', 'KEK-LS ARCA: MagnetSpec widths, TextHeight 0.3 m, Overlap="remove"')]
    for dxf, out, title in figs:
        if os.path.exists(os.path.join(HERE, out)):
            print('skip (exists)', out)
            continue
        plot_dxf.plot(os.path.join(HERE, dxf), os.path.join(HERE, out), title, None, True, 130, (12, 5))
        print('png', out)


if __name__ == '__main__':
    todo = sys.argv[1:] or ['mk', 'run', 'png']
    if 'mk' in todo:
        mk()
    if 'run' in todo:
        gsdir = run_gs()
        run_pysad()
        for r in check(gsdir):
            print('%-36s %8d bytes  TEXT=%d  gs-vs-script %s  python(auto)-vs-script %s  python(dict)-vs-script %s' % r)
    if 'png' in todo:
        png()
