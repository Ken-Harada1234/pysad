"""
Tk の Canvas の組み込み関数 (src/tfTkInter_.c の TkCanvas*、src/tfcanvasclip.f の CanvasClipLine・CanvasSymbolDirect)。
Packages/Canvas.n (ListPlot・Plots・CSR.n の描画) が使う。K64 は Tk の内部構造を直接たたくが、pysad は tkinter の Tcl に create を送る。

  TkCanvasGetTypes[]                         Tk のキャンバス項目の種類名の表 (CanvasTypeID[名前] = 位置)
  TkCanvasPointer[path]                       キャンバスの識別子 (pysad はパス名のまま)
  TkCanvasCreateItem[p, type, 実数..., 文字列(option 名と値)...]   項目を作って id (実数) を返す
  TkCanvasCreateItemDirect[p, type, 実数 / 実数のリスト / オプション規則...]   `Fill->"red"` などは TkOptionLabel で option 名にする
  TkCanvasLastItem[p]                         最後の項目の id
  TkCanvasFindEnclosed[p, {x1,y1,x2,y2}]      矩形に完全に入る項目の id 表
  CanvasClipLine[{{x..},{y..}}, {{xmin,ymin},{xmax,ymax}}, dash]   枠で切り取った折れ線 (破線の分割つき)
  CanvasSymbolDirect[p, tags, {x,y}|{{x..},{y..}}, sym, size, fill, outline, Canvas$Range, Canvas$Offset, oval, poly, rect, line]   記号を描く
"""
from __future__ import annotations

import os
import math

from .sadscript import Expr, Rule, S, Sym, Null, SADError, is_num

# Tk の標準のキャンバス項目 (ctkCanvGetTypes の順は Tk の登録順。ここは SAD 内で一貫していればよい: 添字は CanvasTypeID で引く)
TYPES = ['arc', 'bitmap', 'image', 'line', 'oval', 'polygon', 'rectangle', 'text', 'window']
_A = 0.866025403784439
_W_MIN = 7.0


def _num(x):
    x = float(x)
    return str(int(x)) if x == int(x) and abs(x) < 1e15 else repr(x)


def install(br, B):
    I = br.I

    def tkcall(*args):
        if br.root is None:
            br.create_interp([])
        return br.root.tk.call(*args)

    def path_of(p):
        if not isinstance(p, str):
            raise SADError('General::wrongtype: Character-string for TkCanvasPointer for #1')
        return p

    def option_label(lhs):
        r = I.eval(Expr(S('TkOptionLabel'), [lhs]))
        if not isinstance(r, str):
            raise SADError('General::wrongtype: TkOptionLabel')
        return r.strip().strip('"').strip()

    def value_string(v):
        if isinstance(v, str):
            return v
        if is_num(v):
            return _num(v)
        if isinstance(v, list):
            return ' '.join(value_string(x) for x in v)
        return I.to_str(v)

    def create(path, typeid, coords, opts):
        tn = TYPES[int(typeid) - 1]
        bm = getattr(br, 'bitmap_files', None)
        if bm:
            opts = list(opts)
            for k in range(len(opts) - 1):
                if opts[k] in ('-stipple', '-bitmap', '-activestipple', '-disabledstipple') and opts[k + 1] in bm:
                    opts[k + 1] = bm[opts[k + 1]]          # TkDefineBitmap した名前は @ファイル にする
        try:
            r = tkcall(path, 'create', tn, *[_num(c) for c in coords], *opts)
        except Exception as e:                      # Tk の create の失敗 (オプション名の誤りなど) は SAD ではメッセージだけで Null
            I.out.write('???Tkinter::tclerror: ' + str(e) + chr(10))
            return Null
        return float(r)

    def create_item(a):
        if len(a) < 3:
            raise SADError('General::narg: Number of arguments is expected to be 3 or more in TkCanvasCreateItem')
        path = path_of(a[0])
        typeid = a[1]
        rest = a[2:]
        coords = []
        k = 0
        if is_num(rest[0]):
            while k < len(rest) and is_num(rest[k]):
                coords.append(rest[k])
                k += 1
        elif isinstance(rest[0], list):
            coords = [x for x in rest[0] if is_num(x)]
            k = 1
        else:
            raise SADError('General::wrongtype: Reals or List of Reals for parameter block')
        opts = []
        if k < len(rest):
            if isinstance(rest[k], str):
                opts = [x for x in rest[k:] if isinstance(x, str)]
            elif isinstance(rest[k], list):
                opts = [x for x in rest[k] if isinstance(x, str)]
        return create(path, typeid, coords, opts)

    def expand_args(x, coords, opts):
        """tfGetCanvasArgStk: 実数・実数の表はそのまま座標、Rule は (option 名, 値の文字列)、リストは入れ子を展開"""
        if is_num(x):
            coords.append(x)
        elif isinstance(x, Rule):
            if x.rhs is Null:
                return                      # Fill->Null などは無視 (gs は既定の値: 線は #000000。GeometryPlot の DRIFT の外形)
            opts.append(option_label(x.lhs))
            opts.append(value_string(x.rhs))
        elif isinstance(x, list):
            for y in x:
                expand_args(y, coords, opts)

    def create_item_direct(a):
        if len(a) < 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 or more in TkCanvasCreateItemDirect')
        coords, opts = [], []
        for x in a[2:]:
            if is_num(x) or isinstance(x, (list, Rule)):
                expand_args(x, coords, opts)
            else:
                break
        if not coords and not opts:
            return Null
        return create(path_of(a[0]), a[1], coords, opts)

    B['TkCanvasGetTypes'] = lambda a: list(TYPES)
    B['TkCanvasPointer'] = lambda a: path_of(a[0])
    B['TkCanvasCreateItem'] = create_item
    B['TkCanvasCreateItemDirect'] = create_item_direct

    _last_id = {}

    def last_item(a):
        # K64 の TkCanvasLastItem は Tk の nextId-1 (最後に作った項目の id)。`delete all` のあとも id は続く。
        # 全部消えても id が戻らないよう、その Canvas (パス名と X のウィンドウ id) で見た最大の id を覚えておく
        pth = path_of(a[0])
        r = tkcall(pth, 'find', 'all')
        ids = [int(x) for x in (r if isinstance(r, (tuple, list)) else str(r).split())]
        try:
            key = (pth, str(tkcall('winfo', 'id', pth)))
        except Exception:
            key = (pth, '')
        m = max(max(ids) if ids else 0, _last_id.get(key, 0))
        _last_id[key] = m
        return float(m)
    B['TkCanvasLastItem'] = last_item

    def find_enclosed(a):
        c = [_num(x) for x in a[1]]
        r = tkcall(path_of(a[0]), 'find', 'enclosed', *c)
        ids = [float(int(x)) for x in (r if isinstance(r, (tuple, list)) else str(r).split())]
        return ids if ids else [Null]
    B['TkCanvasFindEnclosed'] = find_enclosed

    # ---------------------------------------------------------------- CanvasClipLine
    def clip_line(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in CanvasClipLine')
        bl, rng = a[0], a[1]
        dash = a[2] if len(a) == 3 else None
        xm = ym = 0.0
        kd = None
        if dash is not None:
            if isinstance(dash, list) and all(is_num(v) for v in dash):
                if len(dash) > 1:
                    kd = dash
            elif is_num(dash):
                if dash > 0:
                    xm = dash
                elif dash < 0:
                    ym = -dash
        if not (isinstance(rng, list) and len(rng) == 2 and all(isinstance(r, list) for r in rng)):
            raise SADError('General::wrongtype: List of Reals for #2[[1]] and #2[[2]]')
        xmin, ymin, xmax, ymax = (round_half(rng[0][0]), round_half(rng[0][1]), round_half(rng[1][0]), round_half(rng[1][1]))
        if not (isinstance(bl, list) and len(bl) == 2 and isinstance(bl[0], list) and isinstance(bl[1], list)
                and len(bl[0]) == len(bl[1])):
            if os.environ.get('PYSAD_TB'):
                print('CLIP', repr(bl)[:300], flush=True)
            raise SADError('General::wrongtype: {{x1, x2, ..}, {y1, y2, ..}}')
        xs, ys = bl
        out = []

        def pin(x, y):
            return xmin <= x <= xmax and ymin <= y <= ymax

        def put(x, y):
            xl = max(xmin, min(xmax, x))
            yl = max(ymin, min(ymax, y))
            if out and out[-2] == xl and out[-1] == yl:
                return
            out.extend([xl, yl])

        def puta(x0, y0, x1, y1):
            xm_ = [xmin, xmax] if x0 <= x1 else [xmax, xmin]
            ym_ = [ymin, ymax] if y0 <= y1 else [ymax, ymin]
            for i in range(2):
                if (y0 - ym_[i]) * (y1 - ym_[i]) <= 0.0:
                    if x0 <= xmin:
                        put(xmin, ym_[i])
                    elif x0 >= xmax:
                        put(xmax, ym_[i])
                if (x0 - xm_[i]) * (x1 - xm_[i]) <= 0.0:
                    if y0 <= ymin:
                        put(xm_[i], ymin)
                    elif y0 >= ymax:
                        put(xm_[i], ymax)

        x0 = round_half(xs[0])
        y0 = round_half(ys[0])
        in0 = pin(x0, y0)
        put(x0, y0)
        for n in range(1, len(xs)):
            x1 = round_half(xs[n])
            y1 = round_half(ys[n])
            in1 = pin(x1, y1)
            ux, uy = x1 - x0, y1 - y0
            t = ([-1.0, -1.0] if ux == 0.0 else [(xmin - x0) / ux, (xmax - x0) / ux]) + \
                ([-1.0, -1.0] if uy == 0.0 else [(ymin - y0) / uy, (ymax - y0) / uy])
            skip = False
            for i in range(4):
                if 0.0 <= t[i] <= 1.0:
                    xa = round_half(x0 + t[i] * ux)
                    ya = round_half(y0 + t[i] * uy)
                    if pin(xa, ya):
                        both = False
                        for j in range(i + 1, 4):
                            if 0.0 <= t[j] <= 1.0:
                                xb = round_half(x0 + t[j] * ux)
                                yb = round_half(y0 + t[j] * uy)
                                if pin(xb, yb):
                                    if t[i] < t[j]:
                                        put(xa, ya)
                                        put(xb, yb)
                                    else:
                                        put(xb, yb)
                                        put(xa, ya)
                                    both = True
                                    break
                        if both:
                            skip = True
                            break
                        if not in0:
                            put(x0, y0)
                            puta(x0, y0, x1, y1)
                        out.extend([xa, ya])        # K64 は put を通さず直接積む
                        break
            if not skip:
                if in1:
                    put(x1, y1)
                elif not in0:
                    put(x0, y0)
                    puta(x0, y0, x1, y1)
            x0, y0, in0 = x1, y1, in1
        put(x0, y0)
        na = len(out)
        if na == 2:
            out = out + [out[0] + xm, out[1] + ym]
            na = 4
        elif na > 4:
            z0 = out[0:2]
            z1 = out[2:4]
            ignore = 0
            res = out[0:2]
            ispa = 2
            res = list(out)
            ispa = 2                                  # Fortran の ispa = isp0+3 (1 始まりの z1 の x の位置)
            for i in range(4, na - 1, 2):
                z = out[i:i + 2]
                ds = 1.0
                if ignore < 4:
                    if (z[0] - z1[0]) * (z1[0] - z0[0]) + (z[1] - z1[1]) * (z1[1] - z0[1]) >= 0:
                        d1 = (z[0] - z0[0]) * (z1[1] - z0[1]) - (z[1] - z0[1]) * (z1[0] - z0[0])
                        ds = d1 * d1 / (1e-4 + (z[0] - z0[0]) ** 2 + (z[1] - z0[1]) ** 2)
                if ds < 0.16:
                    res[ispa:ispa + 2] = z
                    ignore += 1
                else:
                    res[ispa:ispa + 2] = z1
                    z0 = z1
                    ispa += 2
                    ignore = 0
                z1 = z
            res[ispa:ispa + 2] = z1
            out = res[:ispa + 2]
            na = len(out)
        if kd is not None:
            pieces = []
            if na > 1:
                eps = 1e-6
                cur = []
                dash_ = 0.0
                k = 0
                dmax = max(abs(kd[0]), eps)
                z0 = out[0:2]
                cur.extend(z0)
                i = 2
                while i < na:
                    z1 = out[i:i + 2]
                    while True:
                        dz = [z1[0] - z0[0], z1[1] - z0[1]]
                        s = math.hypot(dz[0], dz[1])
                        if s > 0.0:
                            s1 = s + dash_
                            if dash_ >= 0.0:
                                if s1 <= dmax:
                                    cur.extend(z1)
                                    dash_ = s1
                                    z0 = z1
                                    break
                                z0 = [z0[0] + dz[0] * (dmax - dash_) / s, z0[1] + dz[1] * (dmax - dash_) / s]
                                cur.extend([round_half(z0[0]), round_half(z0[1])])
                                pieces.append([float(v) for v in cur])
                                cur = []
                                k = (k + 1) % len(kd)
                                dash_ = min(-abs(kd[k]), -eps)
                                continue
                            else:
                                if s1 > 0.0:
                                    z0 = [z0[0] - dz[0] * dash_ / s, z0[1] - dz[1] * dash_ / s]
                                    cur.extend([round_half(z0[0]), round_half(z0[1])])
                                    dash_ = 0.0
                                    k = (k + 1) % len(kd)
                                    dmax = max(abs(kd[k]), eps)
                                    continue
                                elif s1 == 0.0:
                                    cur.extend([round_half(z1[0]), round_half(z1[1])])
                                    dash_ = 0.0
                                    k = (k + 1) % len(kd)
                                    dmax = max(abs(kd[k]), eps)
                                    z0 = z1
                                    break
                                else:
                                    dash_ = s1
                                    z0 = z1
                                    break
                        break
                    i += 2
                if dash_ > 0.0 and cur:
                    pieces.append([float(v) for v in cur])
            return pieces if pieces else []
        return [float(v) for v in out] if na > 1 else []

    def round_half(v):
        """Fortran の anint (四捨五入、0.5 は 0 から遠い方)"""
        return float(math.floor(abs(v) + 0.5)) * (1.0 if v >= 0 else -1.0)
    B['CanvasClipLine'] = clip_line

    # ---------------------------------------------------------------- CanvasSymbolDirect
    def symbol_direct(a):
        if len(a) != 13:
            raise SADError('General::narg: Number of arguments is expected to be 13 in CanvasSymbolDirect')
        p, tags, pts, symv, size, fill, outline, rng, off, ovl, pol, rec, lin = a
        sym = '1O'
        if isinstance(symv, str) and len(symv) == 2:
            sym = symv
        elif isinstance(symv, Sym) and symv.name == 'Bar':
            sym = 'BA'
        if not isinstance(pts, list) or not (2 <= len(pts) <= 4):
            raise SADError('General::wrongtype: {x, y} or {{x1,..,xn}, {y1,..,yn}} for #3')
        if all(is_num(v) for v in pts):
            xs, ys = [pts[0]], [pts[1]]
        elif isinstance(pts[0], list) and isinstance(pts[1], list) and len(pts[0]) == len(pts[1]):
            xs, ys = pts[0], pts[1]
        else:
            raise SADError('General::wrongtype: {x, y} or {{x1,..,xn}, {y1,..,yn}} for #3')
        m = len(xs)
        xmin, ymin = rng[0][0], rng[0][1]
        xmax, ymax = rng[1][0], rng[1][1]
        yoff = off[1]
        s = size
        if not is_num(s) or not isinstance(fill, str) or not isinstance(outline, str):
            raise SADError('General::wrongtype: size, fill color or outline for CanvasSymbolDirect')
        kat = None
        if isinstance(tags, list):
            if len(tags) != m:
                raise SADError('General::equalleng: tags (#2) and points (#3)')
            kat = tags
        ol = True
        wbar = None
        if sym in ('BA', 'BR') and s < _W_MIN:
            wbar = _num(round_half(max(1.0, s / _W_MIN * 2)))
            ol = False
        isym = {'1O': 1, 'CI': 1, 'CL': 1, '6O': 2, 'LT': 2, '7O': 3, 'RT': 3, '8O': 4, 'DT': 4, '9O': 5, 'UT': 5,
                'BX': 6, 'SQ': 6, 'RH': 7, 'DM': 7, 'DI': 7, 'PL': 8, 'TI': 9, 'ML': 9, 'MU': 9, 'BA': 99, 'BR': 99}.get(sym, 1)
        path = path_of(p)
        for i in range(m):
            x, y = xs[i], ys[i]
            if x < xmin or x > xmax:
                continue
            if isym != 99:
                if y < ymin or y > ymax:
                    continue
            elif max(y, yoff) < ymin or min(y, yoff) > ymax:
                continue
            if isym == 1:
                typ, c = ovl, [x - s, y - s, x + s, y + s]
            elif isym in (2, 3, 4, 5):
                s1 = s * math.sqrt(math.pi / _A / 1.5)
                sa, sh = s1 * _A, s1 * 0.5
                typ = pol
                c = {2: [x - s1, y, x + sh, y - sa, x + sh, y + sa], 3: [x + s1, y, x - sh, y - sa, x - sh, y + sa],
                     4: [x, y + s1, x - sa, y - sh, x + sa, y - sh], 5: [x, y - s1, x - sa, y + sh, x + sa, y + sh]}[isym]
            elif isym == 6:
                s1 = s * math.sqrt(math.pi / 4.0)
                typ, c = pol, [x + s1, y - s1, x + s1, y + s1, x - s1, y + s1, x - s1, y - s1]
            elif isym == 7:
                s1 = s * math.sqrt(math.pi / 2.0)
                typ, c = pol, [x + s1, y, x, y + s1, x - s1, y, x, y - s1]
            elif isym == 8:
                s1 = s * 1.2
                typ = pol
                c = [x + s1, y - 1, x + s1, y + 1, x + 1, y + 1, x + 1, y + s1, x - 1, y + s1, x - 1, y + 1, x - s1, y + 1,
                     x - s1, y - 1, x - 1, y - 1, x - 1, y - s1, x + 1, y - s1, x + 1, y - 1]
            elif isym == 9:
                typ = pol
                c = [x + s + 1, y + s - 1, x + s - 1, y + s + 1, x, y + 1, x - s + 1, y + s + 1, x - s - 1, y + s - 1, x - 1, y,
                     x - s - 1, y - s + 1, x - s + 1, y - s - 1, x, y - 1, x + s - 1, y - s - 1, x + s + 1, y - s + 1, x + 1, y]
            else:
                yc = min(max(ymin, y), ymax)
                yo = min(max(ymin, yoff), ymax)
                if ol:
                    typ, c = rec, [x - s / _W_MIN, yc, x + s / _W_MIN, yo]
                else:
                    typ, c = lin, [x, yc, x, yo]
            opts = ['-outline', outline] if ol else ['-width', wbar]
            opts += ['-fill', fill]
            if isinstance(tags, str):
                opts += ['-tags', tags]
            elif kat is not None and isinstance(kat[i], str):
                opts += ['-tags', kat[i]]
            create(path, typ, c, opts)
        return Null

    B['CanvasSymbolDirect'] = symbol_direct
