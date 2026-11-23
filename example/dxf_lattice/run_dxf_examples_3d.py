# -*- coding: utf-8 -*-
"""3D (磁石の立体外形。Mode->"3D"、Parts->"split") の例: FODO の 1 セル・全周と KEK-LS の標準セル 1 個。台本 DXF_FODO_3d.sad・DXF_KEKLS_3d.sad を作り (mk)、
gs と pysad で走らせて出力の DXF が (1) gs と pysad の台本、(2) pysad の台本と Python 関数 write_dxf(..., mode='3d') で全行一致することを確かめ
(run。あわせて 3DFACE の立体が閉じている・体積が正・body の頂点が 2D の外形と一致するかの検査 `verify`)、mplot3d の図 (PNG、*_3d*.png) を作る (png)。
既存の DXF・PNG は上書きしない (出力は名前に _3d を付けた別ファイルで、すでにあれば止まる)。
使い方: py -3.12 example/dxf_lattice/run_dxf_examples_3d.py [mk] [run] [verify] [png]      (既定は 4 つとも。verify は run の出力を使う)
環境: PYSAD_SAD_PACKAGES (K64 版 SAD の Packages)。gs は example/tools/sadenv.py に従う (WSL。local/wake/gsrun.sh と同じ呼び出し)。"""
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

# 説明用の仮の寸法 (実機の寸法ではない)。FODO は Spec で与える (Height を含む)
FODO_SPEC = ('Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1, Height->.3}, "QUAD"->{Width->.36, Height->.36}, '
             '"SEXT"->{Width->.24, Height->.24}}')
FODO_SPEC_PY = {'B': {'Width': .5, 'PoleWidth': .15, 'Type': 'Ctype1', 'Height': .3}, 'QUAD': {'Width': .36, 'Height': .36},
                'SEXT': {'Width': .24, 'Height': .24}}
EXAMPLES = {
    'DXF_FODO_3d': {
        'lattice': os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'),
        'extra': [],
        'outputs': {
            'DXF_FODO_cell_3d.dxf': ('RING', 'Region->{"SF.1", "SD.1"}, %s, Mode->"3D"' % FODO_SPEC,
                                     {'region': ('SF.1', 'SD.1'), 'spec': FODO_SPEC_PY, 'mode': '3d'}),
            'DXF_FODO_cell_3d_split.dxf': ('RING', 'Region->{"SF.1", "SD.1"}, %s, Mode->"3D", Parts->"split"' % FODO_SPEC,
                                           {'region': ('SF.1', 'SD.1'), 'spec': FODO_SPEC_PY, 'mode': '3d', 'parts': 'split'}),
            'DXF_FODO_ring_3d.dxf': ('RING', 'Mode->"3D"', {'mode': '3d'}),
        }},
    'DXF_KEKLS_3d': {
        'lattice': os.path.join(ROOT, 'example', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad'),
        'extra': ['MagnetSpec_KEKLS.n'],
        'outputs': {
            'DXF_KEKLS_cell_3d.dxf': ('ARCA', 'TextHeight->.12, Mode->"3D"', {'text_height': .12, 'mode': '3d'}),
            'DXF_KEKLS_cell_3d_split.dxf': ('ARCA', 'TextHeight->.12, Mode->"3D", Parts->"split"',
                                            {'text_height': .12, 'mode': '3d', 'parts': 'split'}),
        }},
}


def script_text(name):
    ex = EXAMPLES[name]
    lat = open(ex['lattice'], encoding='utf-8').read()
    if not lat.endswith('\n'):
        lat += '\n'
    s = ['', 'FFS;', '$FORM = "14.9";', 'Get["DrawMagnetsDXF.n"];'] + ['Get["%s"];' % f for f in ex['extra']]
    last = None
    for fn, (use, opts, _kw) in ex['outputs'].items():
        if use != last:
            s += ['USE %s;' % use, 'Print["USE %s: ", LINE["LENGTH"] - 1, " elements, L = ", LINE["S", LINE["LENGTH"]], " m"];' % use]
            last = use
        s += ['DrawMagnetsDXF["%s", %s];' % (fn, opts), 'Print["  wrote %s"];' % fn]
    s += ['Exit;', '']
    return lat + '\n'.join(s)


def mk():
    for name in EXAMPLES:
        out = os.path.join(HERE, name + '.sad')
        if os.path.exists(out):
            raise SystemExit('%s は既にある (上書きしない)。消してから作り直すこと' % out)
        open(out, 'w', encoding='utf-8', newline='\n').write(script_text(name))
        print('made', name + '.sad')


def run_gs(name, timeout=600):
    import sadenv
    gs = sadenv.gs_wsl()
    tmp = tempfile.gettempdir()
    ex = EXAMPLES[name]
    for fn in ex['outputs']:
        if os.path.exists(os.path.join(tmp, fn)):
            os.remove(os.path.join(tmp, fn))
    shutil.copy(os.path.join(HERE, name + '.sad'), os.path.join(tmp, '_dxfex_%s.sad' % name))
    for f in ['DrawMagnetsDXF.n'] + ex['extra']:
        shutil.copy(os.path.join(HERE, f), os.path.join(tmp, f))
    t = tmp.replace(chr(92), '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c',
                        'cd %s && timeout %d %s _dxfex_%s.sad < /dev/null' % (t, timeout, gs, name)], capture_output=True, env=env)
    open(os.path.join(HERE, name + '.gs.txt'), 'wb').write(r.stdout + r.stderr)
    gsdir = os.path.join(tmp, 'dxfex_gs_3d')
    os.makedirs(gsdir, exist_ok=True)
    for fn in ex['outputs']:
        shutil.copy(os.path.join(tmp, fn), os.path.join(gsdir, fn))
    return gsdir


def run_pysad(name):
    import sadenv
    for fn in EXAMPLES[name]['outputs']:
        if os.path.exists(os.path.join(HERE, fn)):
            raise SystemExit('%s は既にある (上書きしない)。消してから走らせること' % fn)
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8', PYSAD_SAD_PACKAGES=sadenv.packages())
    code = ("import os;os.chdir(%r);from pysad.ffscmd import run_sad_string;run_sad_string(open(%r,encoding='utf-8').read(), cwd=%r)"
            % (HERE, name + '.sad', HERE))
    r = subprocess.run(['py', '-3.12', '-c', code], capture_output=True, env=env)
    open(os.path.join(HERE, name + '.pysad.txt'), 'wb').write(r.stdout + r.stderr)


def ffs_of(name, use):
    from pysad.dxf_lattice import load_ffs
    ex = EXAMPLES[name]
    lat = open(ex['lattice'], encoding='utf-8').read() + '\nFFS;\n' + ''.join('Get["%s"];\n' % f for f in ex['extra'])
    return load_ffs(lat, use, cwd=HERE)


def check(name, gsdir):
    from pysad.dxf_lattice import read_dxf, write_dxf
    res = []
    ffs_cache = {}
    for fn, (use, _opts, kw) in EXAMPLES[name]['outputs'].items():
        if use not in ffs_cache:
            ffs_cache[use] = ffs_of(name, use)
        p1 = os.path.join(tempfile.gettempdir(), 'dxfex_py1_' + fn)
        info = write_dxf(ffs_cache[use], p1, **kw)
        sc = os.path.join(HERE, fn)
        nf = sum(1 for e in read_dxf(sc)['entities'] if e['type'] == '3DFACE')
        res.append((fn, os.path.getsize(sc), nf, info['nsolids'], compare(os.path.join(gsdir, fn), sc), compare(p1, sc)))
    return res


def verify(name):
    """3DFACE の立体の検査: (1) 6 枚ずつの立体が閉じている (各有向辺が 1 回、逆向きも 1 回 = 各辺が 2 面で共有され向きがそろう)、(2) 符号つきの体積が正で、
    body は 2D の外形の面積 × Height に一致、(3) body の上面・下面の頂点の (X, Y) が 2D の外形 (同じ引数の平面図) の頂点と一致、
    (4) split は 1 台ごとの部品の体積の和が body 以下で、全部品の bbox が body の bbox の中。戻りは出力ごとの dict"""
    from pysad.dxf_lattice import check_solids, read_dxf, write_dxf
    out = {}
    for fn, (use, _opts, kw) in EXAMPLES[name]['outputs'].items():
        ffs = ffs_of(name, use)
        d3 = read_dxf(os.path.join(HERE, fn))
        sol = check_solids(d3['entities'])
        kw2 = {k: v for k, v in kw.items() if k not in ('mode', 'parts')}
        tmpf = os.path.join(tempfile.gettempdir(), 'dxfex_plan_' + fn)
        write_dxf(ffs, tmpf, **kw2)
        plan = [e for e in read_dxf(tmpf)['entities'] if e['type'] == 'POLYLINE' and e['closed']]
        res = {'nsolids': len(sol), 'all_closed': all(s['closed'] for s in sol), 'min_volume': min(s['volume'] for s in sol)}
        faces = [e for e in d3['entities'] if e['type'] == '3DFACE']
        if kw.get('parts', 'body') == 'body':
            area = lambda P: abs(sum(P[i][0] * P[(i + 1) % 4][1] - P[(i + 1) % 4][0] * P[i][1] for i in range(4))) / 2.0
            ok_xy, dv = len(plan) == len(sol), 0.0
            for i, (s, pl) in enumerate(zip(sol, plan)):
                top = faces[6 * i]['pts']            # 上面 = 6 枚の最初
                bot = faces[6 * i + 1]['pts']        # 下面 = 2 枚目 (逆順)
                ok_xy = ok_xy and [(p[0], p[1]) for p in top] == [tuple(p) for p in pl['pts']]
                ok_xy = ok_xy and [(p[0], p[1]) for p in bot] == [tuple(p) for p in pl['pts'][::-1]]
                h = s['bbox'][5] - s['bbox'][4]
                dv = max(dv, abs(s['volume'] - area(pl['pts']) * h) / (s['volume'] + 1e-30))
            res.update({'vertices_match_2d': ok_xy, 'max_rel_volume_diff_vs_area_x_height': dv})
        else:
            tb = os.path.join(tempfile.gettempdir(), 'dxfex_body_' + fn)
            write_dxf(ffs, tb, **dict(kw, parts='body'))
            isplit = write_dxf(ffs, tmpf, **kw)
            body = check_solids(read_dxf(tb)['entities'])
            # split の部品を 1 台ごとにまとめる (info['solids'] は (要素名, 部品, レイヤー, 高さ) の順。部品の並びは 偏向磁石: yoke… pole pole、四極: yoke×4 pole×4、他: body)
            groups, prev = [], None
            for i, (_en, part, _lay, _h) in enumerate(isplit['solids']):
                if part == 'body' or (part == 'yoke' and prev in (None, 'pole', 'body')):
                    groups.append([])
                groups[-1].append(i)
                prev = part
            ok_in = ok_vol = len(groups) == len(body)
            for g, b in zip(groups, body):
                bb = b['bbox']
                ok_in = ok_in and all(sol[i]['bbox'][0] >= bb[0] - 1e-6 and sol[i]['bbox'][1] <= bb[1] + 1e-6 and sol[i]['bbox'][2] >= bb[2] - 1e-6
                                      and sol[i]['bbox'][3] <= bb[3] + 1e-6 and sol[i]['bbox'][4] >= bb[4] - 1e-6 and sol[i]['bbox'][5] <= bb[5] + 1e-6
                                      for i in g)
                ok_vol = ok_vol and sum(sol[i]['volume'] for i in g) <= b['volume'] * (1 + 1e-9)
            res.update({'parts_total': len(sol), 'bodies': len(body), 'parts_within_body_bbox': ok_in, 'parts_volume_le_body': ok_vol})
        out[fn] = res
    return out


def png():
    import plot_dxf3d
    # (DXF, PNG, タイトル, 3 つの panel: 全体の斜め・平面に近い・拡大)
    figs = [('DXF_FODO_cell_3d.dxf', 'FODO_cell_3d.png', 'FODO, 1 cell (SF.1 to SD.1): solids (Height from Spec), z in m', None),
            ('DXF_FODO_cell_3d_split.dxf', 'FODO_cell_3d_split.png', 'FODO, 1 cell: yoke (layer colors) and pole (gray), Parts="split"', None),
            ('DXF_FODO_ring_3d.dxf', 'FODO_ring_3d.png', 'FODO ring (24 cells), default dimensions, Mode="3D"', None),
            ('DXF_KEKLS_cell_3d.dxf', 'KEKLS_cell_3d.png', 'KEK-LS ARCA cell: solids (illustrative MagnetSpec dimensions)', None),
            ('DXF_KEKLS_cell_3d_split.dxf', 'KEKLS_cell_3d_split.png', 'KEK-LS ARCA cell: yoke and pole (gray), Parts="split"', None)]
    from pysad.dxf_lattice import read_dxf
    for dxf, out, title, _ in figs:
        if os.path.exists(os.path.join(HERE, out)):
            print('skip (exists)', out)
            continue
        d = read_dxf(os.path.join(HERE, dxf))
        faces = [e for e in d['entities'] if e['type'] == '3DFACE']
        xs = [p[0] for f in faces for p in f['pts']]
        ys = [p[1] for f in faces for p in f['pts']]
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        panels = [(38, -70, None)]
        if 'ring' not in dxf:
            # 拡大: 最初の偏向磁石のまわり
            b = [f for f in faces if f['layer'] == 'BEND']
            bx = [p[0] for f in b[:6] for p in f['pts']]
            by = [p[1] for f in b[:6] for p in f['pts']]
            cx, cy, r = (min(bx) + max(bx)) / 2, (min(by) + max(by)) / 2, 1.2 * max(max(bx) - min(bx), max(by) - min(by)) / 2 + 0.3
            panels.append((24, -35, (cx - r, cx + r, cy - r, cy + r)))
        else:
            panels = [(55, -60, None)]
        size = (15, 6) if len(panels) == 2 else (10, 8)
        plot_dxf3d.plot3d(os.path.join(HERE, dxf), os.path.join(HERE, out), title, panels, size, 120)
        print('png', out)


if __name__ == '__main__':
    todo = sys.argv[1:] or ['mk', 'run', 'verify', 'png']
    if 'mk' in todo:
        mk()
    if 'run' in todo:
        for name in EXAMPLES:
            gsdir = run_gs(name)
            run_pysad(name)
            for r in check(name, gsdir):
                print('%-32s %8d bytes  3DFACE=%d  solids=%d  gs-vs-script %s  python-vs-script %s' % r)
    if 'verify' in todo:
        for name in EXAMPLES:
            for fn, r in verify(name).items():
                print(fn, r)
    if 'png' in todo:
        png()
