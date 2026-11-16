"""
小さな組み込み関数と Packages/optimize.n の翻訳。

  DownhillSimplex / Amoeba / Amoebar (optimize.n の Simplex クラス): ダウンヒルシンプレックス法。K64 の SADScript をそのまま移植
  StringToStream / TemporaryName / ListRandom / MemoryCheck / TracePrint ($Input と $Output は 77 と -1)

未対応・差: DownhillSimplex の Output (>0) の途中経過は $FORM='10.6' 形式でなく %10.6f で書く。MemoryCheck は SAD のヒープを持たないので
{0,0,0,0} (内容は意味を持たない)。ListRandom は pysad の乱数が SAD の tran() だけなので {"SAD"}。
"""
from __future__ import annotations

import io
import math
import os
import tempfile

from .sadscript import Sym, Expr, Rule, SADError, S, Null, is_num, INF, truthy


def install(I):
    B = I.builtins

    # ------------------------------------------------------------------ 小物
    def _stringtostream(a):
        if not a or not isinstance(a[0], str):
            raise SADError('General::wrongtype: Argument must be Character-string in StringToStream')
        u = I.next_unit
        I.next_unit += 1
        I.files[u] = io.StringIO(a[0])
        return float(u)
    B['StringToStream'] = _stringtostream

    _tmpcount = [0]

    def _temporaryname(a):
        _tmpcount[0] += 1
        d = tempfile.gettempdir()
        return os.path.join(d, 'tmp__%08d.%s' % (os.getpid() + _tmpcount[0], os.urandom(3).hex()[:6])).replace('\\', '/')
    B['TemporaryName'] = _temporaryname
    B['ListRandom'] = lambda a: ['SAD']
    B['MemoryCheck'] = lambda a: [0.0, 0.0, 0.0, 0.0] + ([[]] if a and is_num(a[0]) and a[0] == 2 else [])
    B['TracePrint'] = lambda a: Null
    I.globals.setdefault('$Input', 77.0)
    I.globals.setdefault('$Output', -1.0)

    # ------------------------------------------------------------------ optimize.n: Simplex クラス
    class Simplex:
        def __init__(self):
            self.prolog = None
            self.limit = None
            self.rr = [1, -1]
            self.wexp = None
            self.exp = 2.0
            self.cont = 1
            self.exl = 1
            self.cl = 1
            self.minq = 0.0
            self.per = 1.0
            self.quick = 0.0
            self.f = None

        def fval(self, p):
            v = I.apply_value(self.f, [list(p)], self.f)
            if not is_num(v):
                raise SADError('DownhillSimplex: function value must be Real')
            return float(v)

        def amoebar(self, p, r, kr):
            """Amoebar[p, f, r, kr]: 点 tp (rr の範囲の重み付き平均) から r だけ kr 番目の点を反射/拡大/縮小した点を試す"""
            rr = self.rr
            pts = [q[1] for q in p]
            lo, hi = self.limit
            # Take[pts, rr] (rr={a,b}: a..b。負は末尾から)
            lp = len(pts)
            a0, b0 = rr
            a0 = a0 - 1 if a0 > 0 else lp + a0
            b0 = b0 if b0 > 0 else lp + b0 + 1
            tp_pts = pts[a0:b0]
            nv = len(pts[0])
            tp = [sum(self.wexp[k] * tp_pts[k][j] for k in range(len(tp_pts))) for j in range(nv)]
            kk = kr - 1 if kr > 0 else lp + kr
            pk = p[kk][1]
            pr = [min(max((1.0 - r) * tp[j] + pk[j] * r, lo[j]), hi[j]) for j in range(nv)]
            if self.prolog is not None:
                pr = I.apply_value(self.prolog, [pr], self.prolog)
            fr = self.fval(pr)
            if fr < p[kk][0]:
                q = list(p)
                q[kk] = [fr, pr]
                q.sort(key=lambda t: (t[0], tuple(t[1])))
                return pr, fr, q
            return pr, fr, p

        def amoeba(self, p):
            p0 = list(p)
            lp = len(p)
            fr0 = pr0 = None
            ret = None
            for kr in range(1, self.cont + 1):
                pt, fr, pr = self.amoebar(p, -1.0, -kr)
                if kr == 1:
                    fr0, pr0 = fr, pr
                if fr <= p[self.exl - 1][0]:
                    pt, fr1, pr1 = self.amoebar(p, -self.exp, -kr)
                    return pr1 if fr1 <= fr else pr
            if self.cont != 1:
                fr, pr = fr0, pr0
            if lp > 2 and fr < p[-1 - self.cl][0]:
                return pr
            pt, fr1, pr1 = self.amoebar(pr, 0.5, -1)
            if fr1 < pr[-1 - self.cl][0]:
                return pr1
            if fr < fr1:
                if fr < p0[-1][0]:
                    p0 = pr
            else:
                if fr1 < p0[-1][0]:
                    p0 = pr1
            p0 = [list(t) for t in p0]
            pl = p0[0][1]
            nv = len(pl)
            prr = [[(pl[j] + q[1][j]) / 2.0 for j in range(nv)] for q in pr0[1:]]
            minq = self.minq
            prr = [[((c > pl[j]) * (pl[j] + minq - (pl[j] - minq)) + (pl[j] - minq)) if (pl[j] - minq < c < pl[j] + minq) else c
                    for j, c in enumerate(pt_)] for pt_ in prr]
            for n in range(lp, 1, -1):
                fv = self.fval(prr[n - 2])
                p0[n - 1] = [fv, prr[n - 2]]
                if self.quick and (fv < p0[n - 2][0] or fv > p0[-1][0]):
                    break
            if self.per != 1:
                p0[0] = [p0[0][0] * self.per + self.fval(pl) * (1.0 - self.per), pl]
            p0.sort(key=lambda t: (t[0], tuple(t[1])))
            return p0

    def _opt(opts, name, default):
        return opts.get(name, default)

    def _downhillsimplex(a):
        if len(a) < 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 or more in DownhillSimplex')
        p0, f = a[0], a[1]
        opts = {}
        for o in a[2:]:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym):
                opts[o.lhs.name] = o.rhs
        p = sorted(([float(t[0]), [float(x) for x in t[1]]] for t in p0), key=lambda t: (t[0], tuple(t[1])))
        nv = len(p[0][1])
        sx = Simplex()
        sx.f = f
        lim = opts.get('VariableRange')
        if lim is None:
            sx.limit = ([-INF] * nv, [INF] * nv)
        else:
            sx.limit = ([float(x) for x in lim[0]], [float(x) for x in lim[1]])
        n = int(opts.get('MaxIteration', max(100, 10 * len(p0))))
        lo = int(opts.get('Output', 6))
        tol = float(opts.get('Tolerance', 1e-6))
        sx.per = float(opts.get('Persistence', 1))
        sx.minq = float(opts.get('MinimumSize', 0))
        sx.exp = float(opts.get('ExpansionRatio', 2))
        sx.cont = int(opts.get('TryReflect', 1))
        sx.quick = opts.get('QuickShrink', 0)
        sx.exl = int(opts.get('ExtendLevel', 1))
        sx.cl = int(opts.get('ContractLevel', 1))
        rr = opts.get('ReflectRange', [1, -1])
        sx.rr = [int(x) for x in rr]
        we = float(opts.get('WeightExponent', 0))
        sx.prolog = None if opts.get('Prolog') in (None,) or (isinstance(opts.get('Prolog'), Sym) and opts['Prolog'].name == 'Identity') else opts['Prolog']
        lp = len(p)
        a0, b0 = sx.rr
        a0 = a0 - 1 if a0 > 0 else lp + a0
        b0 = b0 if b0 > 0 else lp + b0 + 1
        lrr = len(p[a0:b0])
        w = [math.exp(-(we / lrr) * (k + 1)) for k in range(lrr)]
        s = sum(w)
        sx.wexp = [x / s for x in w]
        sx.exl = max(1, min(sx.exl, lp))
        sx.cont = max(1, min(sx.cont, lp))
        nc = int((80 - 3) // 22)
        if lo > 0:
            I.write_unit(lo, 'DownhillSimplex {limit, Option}: ...\n')
        I.globals['CurrentSimplex'] = [[t[0], list(t[1])] for t in p]
        for _ in range(n):
            p = sx.amoeba(p)
            I.globals['CurrentSimplex'] = [[t[0], list(t[1])] for t in p]
            if lo > 0:
                vals = [t[0] for t in p]
                if 2 * nc < len(vals):
                    vals = vals[:nc] + vals[len(vals) - nc:]
                I.write_unit(lo, ' '.join('%10.6f' % v for v in vals) + '\n')
            if (p[-1][0] - p[0][0]) / (1e-100 + abs(p[-1][0]) + abs(p[0][0])) < tol:
                break
        I.globals['Status'] = 'Stop'
        return [[t[0], list(t[1])] for t in p]
    B['DownhillSimplex'] = _downhillsimplex


def install_small(I):
    """Short[x, n] (出力: x を書いて Null。n は無視 — gs で確認)、Inequality[a, op1, b, op2, c, …] (連なった比較)"""
    B = I.builtins
    from .sadscript import Expr, S, Sym, truthy

    def short(a):
        if not a:
            return Null
        I.write_unit(6, I.to_str(a[0]) + chr(10))
        return Null
    B['Short'] = short

    def inequality(a):
        if len(a) < 3 or len(a) % 2 == 0:
            return Expr(S('Inequality'), list(a))
        res = 1.0
        for k in range(1, len(a) - 1, 2):
            op = a[k]
            if not isinstance(op, Sym) or op.name not in B:
                return Expr(S('Inequality'), list(a))
            r = B[op.name]([a[k - 1], a[k + 1]])
            if truthy(r) is not True:
                return 0.0
        return res
    B['Inequality'] = inequality

    def messagename(a):
        # gs: MessageName[f, "tag"] の関数形は f を見ず、General::tag の本文 (未定義なら下の文字列) を返す。f::tag の形 (Message の本文) とは別
        if len(a) == 2 and isinstance(a[1], str):
            return 'Undefined MessageName General::"%s"' % a[1]
        return Expr(S('MessageName'), list(a))
    B['MessageName'] = messagename

    # Fork[]: 親では子のプロセス番号 (正の実数)、子では 0。pysad は fork しないので常に親の側 (子のための分岐は実行しない)。
    # 並列の計算は Parallelize (fileops.py) が逐次で行う
    B['Fork'] = lambda a: float(os.getpid())


def install_matrix(I):
    """Packages/Matrix.n の SingularValues と MatrixBalance (SingularValues1 は K64 の組み込み。特異ベクトルの符号は numpy の SVD を、
    各 u の先頭成分が負になるよう揃えて K64 の出力に合わせる)"""
    import numpy as np
    B = I.builtins

    def _opts(rest):
        o = {}
        for r in rest:
            if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                o[r.lhs.name] = r.rhs
        if 'Threshold' in o and 'Tolerance' not in o:
            o['Tolerance'] = o['Threshold']
        return o

    def singular_values(a):
        o = _opts(a[1:])
        e = float(o.get('Tolerance', 1e-8))
        inv = bool(o.get('Inverse', True)) and o.get('Inverse', True) is not False
        srt = bool(o.get('Sort', True)) and o.get('Sort', True) is not False
        m = np.array(a[0], dtype=float)
        u, w, vt = np.linalg.svd(m, full_matrices=False)
        k = len(w)
        uu = [u[:, i] for i in range(k)]
        vv = [vt[i, :] for i in range(k)]
        for i in range(k):
            if uu[i][0] > 0:
                uu[i] = -uu[i]
                vv[i] = -vv[i]
        wmax = float(w.max()) if k else 0.0
        ww = []
        for x in w:
            if x <= e * wmax or x == 0.0:
                ww.append(0.0)
            else:
                ww.append(1.0 / x if inv else float(x))
        trip = list(zip(ww, [list(map(float, x)) for x in uu], [list(map(float, x)) for x in vv]))
        if srt:
            trip.sort(key=lambda t: t[0], reverse=not inv)
        return [[t[1] for t in trip], [t[0] for t in trip], [t[2] for t in trip]]
    B['SingularValues'] = singular_values

    def matrix_balance(a):
        m = np.array(a[0], dtype=float)
        lm, wm = m.shape
        l = min(lm, wm)
        rows = np.max(np.abs(m[:l, :]), axis=1)
        cols = np.max(np.abs(m.T[:l, :]), axis=1)
        r = np.sqrt(np.maximum(1e-100, rows) / np.maximum(1e-100, cols))
        left = 1.0 / np.concatenate([r, np.ones(lm - l)])
        right = np.concatenate([r, np.ones(wm - l)])
        return [[float(x) for x in left], [float(x) for x in right]]
    B['MatrixBalance'] = matrix_balance


def install_system(I):
    """K64 の小さな組み込み (tfinitn.f): Directory SetDirectory Replace FindIndex Attributes ToInputString Context SetEnv GetEnv GetGID
    GarbageCollect ToDate FromDate と $Version。pysad に不要なもの (メモリ・共有メモリ) は含めない"""
    import datetime
    import numpy as np
    from .sadscript import _Lit
    B = I.builtins

    B['Directory'] = lambda a: (I.cwd or os.getcwd()).replace(chr(92), '/')

    def setdir(a):
        p = a[0]
        if not os.path.isabs(p):
            p = os.path.join(I.cwd, p)
        if not os.path.isdir(p):
            raise SADError('General::fileopen: SetDirectory: cannot change directory')
        I.cwd = p
        return p.replace(chr(92), '/')
    B['SetDirectory'] = setdir

    def replace(a):
        x, rules = a[0], a[1]
        rl = [rules] if isinstance(rules, Rule) else list(rules)
        flat = []
        for r in rl:
            if isinstance(r, list):
                flat.extend(r)
            else:
                flat.append(r)
        for r in flat:
            if not isinstance(r, Rule):
                continue
            lhs = r.lhs.v if isinstance(r.lhs, _Lit) else r.lhs
            rhs = r.rhs.v if isinstance(r.rhs, _Lit) else r.rhs
            if I.same(x, lhs):
                return I.eval(rhs) if r.delayed else rhs
            b = {}
            if I._match_one(lhs, x, b):
                return I.eval(I.subst_bindings(rhs, b))
        return x
    B['Replace'] = replace

    def findindex(a):
        xl = np.asarray(a[0], dtype=float)
        n = len(xl)

        def one(v):
            i = int(np.searchsorted(xl, v, side='right'))
            return float(min(max(i, 1), n - 1))
        x = a[1]
        return [one(v) for v in x] if isinstance(x, list) else one(x)
    B['FindIndex'] = findindex

    def attributes(a):
        s = a[0]
        nm = s.name if isinstance(s, Sym) else str(s)
        at = I.attributes.get(nm)
        out = []
        if at:
            out = [S(t) for t in sorted(at)]
        elif nm in B:
            out = [S('HoldNone'), S('Constant'), S('Immediate')] + ([S('OrderLess')] if nm in ('Plus', 'Times') else [])
        return out or [S('None')]
    B['Attributes'] = attributes
    B['ToInputString'] = lambda a: I.to_str(a[0], quote=True)
    B['Context'] = lambda a: 'System`' if not a else '`'
    B['SetEnv'] = lambda a: (os.environ.__setitem__(str(a[0]), str(a[1])), Null)[1]
    B['GetEnv'] = lambda a: os.environ.get(str(a[0]), '')
    B['GetGID'] = lambda a: float(os.getgid()) if hasattr(os, 'getgid') else 0.0
    B['GarbageCollect'] = lambda a: Null
    I.globals.setdefault('$Version', '1.4.4.2k64')

    t0 = datetime.datetime(1900, 1, 1)

    def todate(a):
        d = t0 + datetime.timedelta(seconds=float(a[0]))
        return [float(d.year), float(d.month), float(d.day), float(d.hour), float(d.minute), float(d.second) + d.microsecond * 1e-6]
    B['ToDate'] = todate

    def fromdate(a):
        if not a or not isinstance(a[0], list):
            # FromDate[] / FromDate[Tick->dt]: 現在時刻 (1900-01-01 からの秒、ローカル時刻)。Tkinter.n・KBMainFrame.n が使う
            return (datetime.datetime.now() - t0).total_seconds()
        y, mo, d, h, mi, s = [float(v) for v in a[0]]
        dt = datetime.datetime(int(y), int(mo), int(d), int(h), int(mi), 0) + datetime.timedelta(seconds=s)
        return (dt - t0).total_seconds()
    B['FromDate'] = fromdate


def install_spline(I):
    """K64 の Spline$[{{x1..},{y1..}}, Derivative->{dy1,dy2}] (spline.f: 自然スプライン。戻りは {x, y, y''/6})。
    Spline・SplineData は Packages/D.n (自動読み込み) がこれを使う。Derivative->Periodic (mode 4) は未対応"""
    import numpy as np
    B = I.builtins

    def spline1(a):
        if len(a) not in (1, 2) or not isinstance(a[0], list) or len(a[0]) != 2:
            raise SADError('General::wrongval: Argument must be {{x1, ..}, {y1, ..}}')
        x = np.array(a[0][0], dtype=float)
        y = np.array(a[0][1], dtype=float)
        m = len(x)
        d1 = d2 = None
        if len(a) == 2:
            o = a[1]
            if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'Derivative':
                v = o.rhs
                if isinstance(v, Sym) and v.name in ('Periodic', '`Periodic'):
                    raise SADError('Spline$: Derivative->Periodic is not supported')
                if not isinstance(v, list) or len(v) != 2:
                    raise SADError('General::narg: Number of arguments is expected to be 1 + Derivative->{dy1,dy2}')
                d1 = float(v[0]) if is_num(v[0]) else None
                d2 = float(v[1]) if is_num(v[1]) else None
            else:
                raise SADError('General::narg: Number of arguments is expected to be 1 + Derivative->{dy1,dy2}')
        if m <= 2:
            ddy = np.zeros(max(m, 3))
        else:
            from scipy.interpolate import CubicSpline
            bc = ((1, d1) if d1 is not None else (2, 0.0), (1, d2) if d2 is not None else (2, 0.0))
            ddy = np.zeros(max(m, 3))
            ddy[:m] = CubicSpline(x, y, bc_type=bc)(x, 2) / 6.0
        return [[float(t) for t in x], [float(t) for t in y], [float(t) for t in ddy]]
    B['Spline$'] = spline1
