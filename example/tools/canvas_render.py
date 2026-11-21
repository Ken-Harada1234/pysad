# -*- coding: utf-8 -*-
"""描画の例の図を作る道具 (example/11_graphics と 12_gui_tkinter が使う)。

SAD (gs) も pysad も、描画は Tk の Canvas の項目 (line, rectangle, text, ...) として作る。
gs の窓 (X11) は撮れないので、台本の最後に「Canvas の全項目を Tcl で書き出す」文を足して走らせ、
書き出した項目 (種類・座標・色・文字) を matplotlib で再描画して PNG にする。gs と pysad の両方でこれをして、左右に並べる。

使い方 (Python):
    from canvas_render import make_figure
    make_figure('example/11_graphics/fn_Plot.sad', 'example/11_graphics/fn_Plot_plot.png')
元の .sad は Exit; の直前に書き出しの文を足した一時ファイルで走らせる (元は書き換えない)。
"""
import os, sys, re, subprocess, shutil, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import sadenv
GS = sadenv.gs_wsl()
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
TMP = tempfile.gettempdir()

# Tcl: 全 Canvas の全項目を FN に書く (1 行 1 項目、Tcl のリスト)
DUMP = ('TclEval["update; set f [open FN w]; '
        'proc walk {w} {set r {}; if {[winfo class $w] eq {Canvas}} {lappend r $w}; '
        'foreach c [winfo children $w] {set r [concat $r [walk $c]]}; return $r}; '
        'foreach cv [walk .] {puts $f [list CANVAS $cv [$cv cget -width] [$cv cget -height]]; '
        'foreach i [$cv find all] {set o {}; '
        'foreach k {-fill -outline -width -text -anchor -dash -font -justify -smooth} '
        '{if {![catch {$cv itemcget $i $k} v]} {lappend o $k $v}}; '
        'puts $f [list ITEM $cv $i [$cv type $i] [$cv coords $i] $o [$cv gettags $i]]}}; close $f"];')


def _tmp_sad(path, fn):
    src = open(path, encoding='utf-8').read()
    i = src.rstrip().rfind('Exit;')
    dump = DUMP.replace('FN', fn)
    return src[:i] + dump + '\n' + src[i:]


def run_dump(path, tag, to=120, env_extra=None):
    """path の .sad を gs と pysad で走らせ、Canvas の項目を書き出す。(gs のファイル, pysad のファイル) を返す"""
    path = os.path.abspath(path); name = os.path.basename(path)[:-4]
    out = {}
    # gs
    fn_gs = os.path.join(TMP, '_cv_%s_gs.txt' % name)
    if os.path.exists(fn_gs): os.remove(fn_gs)
    fn_w = fn_gs.replace('\\', '/').replace('C:', '/mnt/c')
    sad = os.path.join(TMP, '_cvg_%s.sad' % name)
    open(sad, 'w', encoding='utf-8', newline='\n').write(_tmp_sad(path, fn_w))
    t = TMP.replace('\\', '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c',
                        'cd %s && timeout %d %s _cvg_%s.sad < /dev/null' % (t, to, GS, name)], capture_output=True, env=env)
    out['gs'] = fn_gs if os.path.exists(fn_gs) else None
    out['gs_log'] = (r.stdout + r.stderr).decode('utf-8', 'replace')
    # pysad
    fn_py = os.path.join(TMP, '_cv_%s_py.txt' % name)
    if os.path.exists(fn_py): os.remove(fn_py)
    sad2 = os.path.join(TMP, '_cvp_%s.sad' % name)
    open(sad2, 'w', encoding='utf-8', newline='\n').write(_tmp_sad(path, fn_py.replace('\\', '/')))
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8', PYSAD_SAD_PACKAGES=sadenv.packages())
    env.update(env_extra or {})
    r = subprocess.run(['py', '-3.12', '-c',
                        "from pysad.ffscmd import run_sad_string;import sys;run_sad_string(open(sys.argv[1],encoding='utf-8').read())", sad2],
                       capture_output=True, env=env, timeout=to * 3)
    out['pysad'] = fn_py if os.path.exists(fn_py) else None
    out['pysad_log'] = (r.stdout + r.stderr).decode('utf-8', 'replace')
    return out


def read_items(fn):
    """書き出しファイルを [(canvas, w, h, [item...])] に。item = dict(type, coords, opt, tags)"""
    import tkinter
    tcl = tkinter.Tcl()
    canv = []
    buf = ''
    for ln in open(fn, encoding='utf-8', errors='replace'):
        buf += ln.rstrip('\n')
        if not buf: continue
        try:
            p = tcl.splitlist(buf)      # 文字項目の中に改行があると 1 項目が複数行になる。閉じるまで足す
        except tkinter.TclError:
            buf += '\n'; continue
        buf = ''
        if p[0] == 'CANVAS':
            canv.append([p[1], float(p[2]), float(p[3]), []])
        elif p[0] == 'ITEM':
            o = tcl.splitlist(p[5]); opt = {o[i]: o[i + 1] for i in range(0, len(o), 2)}
            canv[-1][3].append(dict(id=p[2], type=p[3], coords=[float(x) for x in tcl.splitlist(p[4]) if x.lower() not in ('nan', '-nan', 'inf', '-inf')], opt=opt, tags=tcl.splitlist(p[6])))
    return canv


def _col(c, default=None):
    if not c: return default
    if c.startswith('System'):
        return {'SystemButtonText': 'black', 'SystemWindowText': 'black', 'SystemWindow': 'white', 'SystemButtonFace': '#f0f0f0'}.get(c, 'black')
    c2 = c.replace(' ', '').lower()
    if c2.startswith('#'): return c
    from matplotlib.colors import is_color_like
    if is_color_like(c2): return c2
    c3 = c2.rstrip('0123456789')          # X11 の tomato3 などは番号を落とす
    return c3 if is_color_like(c3) else 'gray'


def draw(ax, items, W, H, title=None):
    from matplotlib.patches import Rectangle, Ellipse, Polygon
    for it in items:
        t, co, o = it['type'], it['coords'], it['opt']
        if t in ('rectangle', 'oval') and len(co) != 4: continue     # NaN を含む座標は描かない
        if len(co) % 2: co = co[:-1]
        fill = _col(o.get('-fill')); outl = _col(o.get('-outline'))
        wd = float(o.get('-width', 1) or 1)
        if t == 'line' and len(co) >= 4:
            xs = co[0::2]; ys = co[1::2]
            ax.plot(xs, ys, color=fill or 'black', lw=min(max(wd * 0.6, 0.3), 4), solid_capstyle='butt')
        elif t == 'rectangle':
            x0, y0, x1, y1 = co
            ax.add_patch(Rectangle((min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0), fc=fill or 'none', ec=outl or 'none', lw=min(wd * 0.6, 3)))
        elif t == 'oval':
            x0, y0, x1, y1 = co
            ax.add_patch(Ellipse(((x0 + x1) / 2, (y0 + y1) / 2), abs(x1 - x0), abs(y1 - y0), fc=fill or 'none', ec=outl or 'none', lw=min(wd * 0.6, 3)))
        elif t == 'polygon':
            pts = list(zip(co[0::2], co[1::2]))
            ax.add_patch(Polygon(pts, closed=True, fc=fill or 'none', ec=outl or 'none', lw=min(wd * 0.6, 3)))
        elif t == 'text':
            txt = o.get('-text', '')
            if txt.startswith('\\@'): txt = '(LaTeX)'
            an = o.get('-anchor', 'center')
            ha = 'left' if 'w' in an else ('right' if 'e' in an else 'center')
            va = 'top' if an.startswith('n') else ('bottom' if an.startswith('s') else 'center')
            m = re.search(r'-(\d+)', o.get('-font', '')) or re.search(r'\s(\d+)\s*$', o.get('-font', ''))
            fs = 5.5
            ax.text(co[0], co[1], txt, ha=ha, va=va, fontsize=fs, color=fill or 'black', clip_on=True)
    pts = [it['coords'] for it in items if len(it['coords']) >= 2]
    xs = [v for c in pts for v in c[0::2]]; ys = [v for c in pts for v in c[1::2]]
    if xs and ys:
        ax.set_xlim(min(xs) - 5, max(xs) + 5); ax.set_ylim(max(ys) + 5, min(ys) - 5)
    else:
        ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect('equal'); ax.axis('off')
    if title: ax.set_title(title, fontsize=9)


def render(files, out_png, labels=('gs', 'pysad'), canvas_index=None):
    """files: [gs の書き出し, pysad の書き出し]。左右に並べて PNG にする。canvas が複数あれば最大の面積のものを使う"""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    panels = []
    for fn, lab in zip(files, labels):
        if not fn:
            panels.append((lab + ' (dump failed)', None)); continue
        cv = read_items(fn)
        if not cv:
            panels.append((lab + ' (no canvas)', None)); continue
        if canvas_index is None:
            k = max(range(len(cv)), key=lambda i: len(cv[i][3]))
        else:
            k = canvas_index
        c = cv[k]
        panels.append(('%s: %d items' % (lab, len(c[3])), c))
    n = len(panels)
    fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 3.9))
    if n == 1: axs = [axs]
    for ax, (ttl, c) in zip(axs, panels):
        if c is None:
            ax.axis('off'); ax.set_title(ttl, fontsize=9); continue
        _, W, H, items = c
        draw(ax, items, W, H, ttl)
    fig.tight_layout(); fig.savefig(out_png, dpi=70); plt.close(fig)
    return panels


def make_figure(sad, png, **kw):
    r = run_dump(sad, 'x', **kw)
    pan = render([r['gs'], r['pysad']], png)
    return r, [(t, (len(c[3]) if c else None)) for t, c in pan]


if __name__ == '__main__':
    reuse = '--reuse' in sys.argv
    for p in [a for a in sys.argv[1:] if not a.startswith('--')]:
        name = os.path.basename(p)[:-4]
        try:
            if reuse:
                pan = render([os.path.join(TMP, '_cv_%s_gs.txt' % name), os.path.join(TMP, '_cv_%s_py.txt' % name)], p[:-4] + '_plot.png')
                info = [(t, (len(c[3]) if c else None)) for t, c in pan]
            else:
                r, info = make_figure(p, p[:-4] + '_plot.png')
        except Exception as e:
            info = 'ERROR %r' % e
        print(p, info, flush=True)
