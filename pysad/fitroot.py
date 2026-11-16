"""
FindRoot / Fit / FitGaussian / FitEmit / NIntegrate。K64 の tffindroot.f (tffindroot, tfnewton, tffit, tffit1)、Packages/Fit.n、Packages/D.n の直訳。

K64 は式を記号微分 (D) して導関数を作り、微分できなかったときだけ数値微分に落ちる。ここでは常に数値微分の経路を使う
(tfnewton の svi=max(|v|*frac, sv)、tffit1 の svi=min(max(svmin,|v|*frac),(vmax-vmin)*factmin) をそのまま移植)。
そのため解析微分で解く gs とは、線形な式で相対 1e-9 程度、非線形でも収束値が数桁 (反復の打ち切り条件に依存) 違いうる。
未対応: 複素数の FindRoot (gs も壊れた値: docs/19 B61)、Fit の Cutoff 以外の位置引数形式 {{x1,x2..},y,dy} (ベクトル独立変数)。FitGaussian の Plot オプションは実装済み。
"""
from __future__ import annotations

import math

import numpy as np

from .sadscript import Sym, Expr, Rule, Frame, SADError, S, Null, is_num, INF
from .tsvd import tsvd

_ALARGE = 1.0e300


def _fmincube(f1, f2, g1, g2, g0, dg):
    """tffindroot.f fmincube"""
    if f2 == 0.0:
        f = -0.5 * f1 * dg / ((g1 - g0) / f1 - dg)
    else:
        a = ((g1 - g0) / f1 ** 2 - (g2 - g0) / f2 ** 2) / (f1 - f2) + dg / f1 / f2
        b = (-f2 * (g1 - g0) / f1 ** 2 + f1 * (g2 - g0) / f2 ** 2) - dg * (f1 + f2) / f1 / f2
        s = math.sqrt(max(0.0, b * b - 3.0 * a * dg))
        f = -dg / (s + b) if b > 0.0 else (s - b) / 3.0 / a
    return min(0.5772156649015329 * f1, max(f1 / 16.0, f))


def _tinvgr(an):
    """tffindroot.f tinvgr: Q(an/2, x/2) = erfc(1/sqrt2) の解 x*2"""
    from scipy.special import gammaincc, gamma
    if an == 1.0:
        return 1.0
    c0 = math.erfc(math.sqrt(0.5))
    anh = an * 0.5
    x0 = anh
    gn = gamma(anh)
    df = gammaincc(anh, x0) - c0
    it = 0
    while abs(df) > 1e-14 and it < 200:
        dfdx = -math.exp(-x0) * x0 ** (anh - 1.0) / gn
        x0 = x0 - df / dfdx
        df = gammaincc(anh, x0) - c0
        it += 1
    return x0 * 2.0


def install(I):
    B = I.builtins

    def ev(expr, env):
        fr = Frame()
        fr.update(env)
        I.frames.append(fr)
        try:
            return I.eval(expr)
        finally:
            I.frames.pop()

    def head_is(e, name):
        return isinstance(e, Expr) and isinstance(e.head, Sym) and e.head.name == name

    def parse_vars(specs):
        """{var, ini} / {var, ini, {min, max}} → (names, v0, vmin, vmax)"""
        names, v0, vmin, vmax = [], [], [], []
        from .sadscript import Seq, _HeldArg, _Lit
        ex = []
        for sp in specs:
            # FitPlot[data_, fun_, x_, par0__] の Fit[data, fun, x, par0]: 配列変数 par0 (Seq) は引数に展開する
            if isinstance(sp, Sym):
                try:
                    v_ = I.lookup(sp.name)
                except KeyError:
                    v_ = None
                if isinstance(v_, Seq):
                    ex.extend(v_.items)
                    continue
            ex.append(sp)
        specs = [(q.node if isinstance(q, _HeldArg) else q.v if isinstance(q, _Lit) else q) for q in ex]
        for sp in specs:
            if isinstance(sp, list):
                sp = Expr(S('List'), list(sp))             # 評価済みの {a, 1} (FitPlot が par0 を Fit に渡す) も同じ
            if not head_is(sp, 'List') or len(sp.args) not in (2, 3) or not isinstance(sp.args[0], Sym):
                import os
                if os.environ.get('PYSAD_TB'):
                    print('FITSPEC', repr(sp)[:200], [repr(q)[:80] for q in specs], flush=True)
                raise SADError('General::wrongtype: Argument must be {var, ini} or {var, ini, {min, max}}')
            names.append(sp.args[0].name)
            ini = I.eval(sp.args[1])
            if not is_num(ini):
                raise SADError('General::wrongval: initial value must be Real')
            if len(sp.args) == 3:
                lim = I.eval(sp.args[2])
                if not isinstance(lim, list) or len(lim) != 2:
                    raise SADError('General::wrongtype: Argument must be {var, ini, {min, max}}')
                lo, hi = float(lim[0]), float(lim[1])
            else:
                lo, hi = -_ALARGE, _ALARGE
            v0.append(min(hi, max(lo, float(ini))))
            vmin.append(lo)
            vmax.append(hi)
        return names, np.array(v0), np.array(vmin), np.array(vmax)

    def get_opts(args):
        """末尾の Rule (評価して) を辞書に。位置引数のオプション (実数 4 つの List) は {MaxIterations, AccuracyGoal, Trace, D}"""
        opts = {}
        rest = []
        for a in args:
            if isinstance(a, Rule) or head_is(a, 'Rule'):
                r = I.eval(a)
                if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                    opts[r.lhs.name] = r.rhs
                continue
            rest.append(a)
        return opts, rest

    # ------------------------------------------------------------------ FindRoot
    def _findroot(args):
        if len(args) < 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 or more in FindRoot')
        eq = args[0]
        opts, rest = get_opts(args[1:])
        names, v0, vmin, vmax = parse_vars(rest)
        maxi = int(opts.get('MaxIterations', 50))
        eps = float(opts.get('AccuracyGoal', 1e-20))
        trace = bool(opts.get('Trace', 0)) and opts.get('Trace') is not False
        # 方程式 f == g か、その List
        if head_is(eq, 'Equal') and len(eq.args) == 2:
            pairs = [(eq.args[0], eq.args[1])]
            lst = False
        elif head_is(eq, 'List') and all(head_is(e, 'Equal') and len(e.args) == 2 for e in eq.args):
            pairs = [(e.args[0], e.args[1]) for e in eq.args]
            lst = True
        else:
            raise SADError('General::wrongtype: Argument must be f == g or List of them')
        neq = len(pairs)
        nvar = len(names)

        def resid(v):
            env = {n: float(x) for n, x in zip(names, v)}
            f = np.zeros(neq)
            am = d = 0.0
            for i, (l, r) in enumerate(pairs):
                a = I.eval_in(l, env) if hasattr(I, 'eval_in') else ev(l, env)
                b = I.eval_in(r, env) if hasattr(I, 'eval_in') else ev(r, env)
                if isinstance(a, (list, Expr)) or isinstance(b, (list, Expr)) or isinstance(a, complex) or isinstance(b, complex):
                    raise SADError('General::wrongval: results must be numbers')
                f[i] = a - b
                d += f[i] ** 2
                am += a * a + b * b
            if trace:
                print(f'FindRoot Residual: {d}', file=I.out)
            return f, am, d

        frac, svmin, svdtol, factmin = 1e-7, 1e-7, 1e-5, 1e-4
        f0, am, d0 = resid(v0)
        goal = am * eps
        sv = (svmin + float(np.sum(np.abs(v0))) * frac) / nvar
        fact = 1.0
        d1 = d0
        fact1 = 0.0
        it = 0
        while True:
            if d0 < goal:
                break
            it += 1
            if it > maxi:
                break
            v = v0.copy()
            a = np.zeros((neq, nvar))
            for i in range(nvar):
                svi = max(abs(v[i]) * frac, sv)
                v[i] = v0[i] + svi
                f, _am, dd = resid(v)
                a[:, i] = (f - f0) / svi
                if dd < d0:
                    v0 = v.copy()
                    f0 = f
                    d0 = dd
                else:
                    v[i] = v0[i]
            a0 = a.copy()
            dv = tsvd(a, -f0, svdtol)
            dg = 0.0
            for i in range(neq):
                dg += f0[i] * float(np.dot(a0[i, :], dv))
            dg *= 2.0
            done = False
            while True:
                v = np.minimum(vmax, np.maximum(vmin, v0 + dv * fact))
                f0n, _am, d = resid(v)
                if d < d0:
                    v0 = v
                    f0 = f0n
                    d0 = d
                    fact = min(1.0, fact * 4.0)
                    fact1 = 0.0
                    d1 = d
                    break
                it += 1
                if it > maxi:
                    done = True
                    break
                d2 = d1
                d1 = d
                fact2 = fact1
                fact1 = fact
                fact = _fmincube(fact1, fact2, d1, d2, d0, dg)
                if fact < factmin:
                    done = True
                    break
            if done:
                break
        out = [Rule(S(n), float(x)) for n, x in zip(names, v0)]
        out.append(Rule(S('Residual'), float(d0)))
        return out
    B['FindRoot'] = _findroot

    # ------------------------------------------------------------------ Fit
    def _fit_core(data, model, xname, names, v0, vmin, vmax, maxi=40, cutoff=0.0, eps=1e-9, xmat=None, xnames=None):
        """data: (m, n) の配列 (n=2: x,y / 3: x,y,dy / 4: x,y,dx,dy)。model は式 (評価器で実行)。戻り: dict"""
        m, n = data.shape
        nvar = len(names)
        xs = data[:, 0]
        ys = data[:, 1]
        sig = data[:, 2] if n == 3 else (data[:, 3] if n == 4 else None)

        def model_at(v, k=None):
            env = {nm: float(x) for nm, x in zip(names, v)}
            out = np.zeros(m)
            for i in range(m):
                if xmat is not None:
                    env.update({nm_: float(x_) for nm_, x_ in zip(xnames, xmat[i])})        # ベクトルの独立変数 (Fit.n の {x1, x2, ..})
                else:
                    env[xname] = float(xs[i])
                val = ev(model, env)
                if isinstance(val, complex):
                    out[i] = _ALARGE
                elif not is_num(val):
                    raise SADError('General::wrongval: result of equation must be Real number')
                else:
                    out[i] = val
            return out

        def evalfit(v, cut):
            df = model_at(v) - ys
            if sig is not None:
                df = df / sig
            if cut != 0.0:
                r = float(np.sum(np.clip(df, -cut, cut) ** 2))
            else:
                r = float(np.sum(df * df))
            return df, r

        frac, factmin, factmin1, svmin, svdeps = 1e-7, 1e-4, 1e-4, 1e-7, 1e-4

        def gammaq(a, x):
            from scipy.special import gammaincc
            return float(gammaincc(a, x))

        v0 = np.minimum(vmax, np.maximum(vmin, v0))
        v00 = v0.copy()
        it = 0
        ajump = 1.0
        dbest = _ALARGE
        acalc = False
        abest = None
        vbest = v0.copy()
        a0 = None
        cut = cutoff
        goto21 = True
        newton = True
        d00 = 0.0
        fact = 1.0
        d1 = 0.0
        fact1 = 0.0
        df0 = None
        d0 = 0.0
        state = 'S21'
        while True:
            if state == 'S21':
                v0 = np.minimum(vmax, np.maximum(vmin, v0))
                df0, d0 = evalfit(v0, cut)
                if d0 < dbest:
                    dbest = d0
                    vbest = v0.copy()
                    if acalc:
                        abest = a0.copy()
                newton = True
                state = 'S11'
            if state == 'S11':
                fact = 1.0
                d1 = d0
                fact1 = 0.0
                d00 = d0
                state = 'S1'
            if state == 'S1':
                it += 1
                if it > maxi:
                    break
                v = v0.copy()
                a = np.zeros((m, nvar))
                for i in range(nvar):
                    svi = min(max(svmin, abs(v0[i]) * frac), (vmax[i] - vmin[i]) * factmin)
                    if vmax[i] - v0[i] >= v0[i] - vmin[i]:
                        svi = min(svi, (vmax[i] - v0[i]) * 0.5)
                    else:
                        svi = max(-svi, (vmin[i] - v0[i]) * 0.5)
                    v[i] = v0[i] + svi
                    df2, db = evalfit(v, 0.0)
                    a[:, i] = (df2 - df0) / svi
                    if db < d0:
                        d0 = db
                        v0[i] = v0[i] + svi
                        df0 = df2
                    else:
                        v[i] = v0[i]
                a0 = a.copy()
                if not acalc:
                    abest = a0.copy()
                    acalc = True
                if newton:
                    dv = tsvd(a, -df0, svdeps)
                else:
                    dv = np.array([-float(np.sum(df0 * a0[:, i])) for i in range(nvar)])
                dg = 2.0 * float(np.sum(df0 * (a0 @ dv)))
                if abs(dg) < d0 * eps:
                    good = gammaq((m - nvar) * 0.5, d0 * 0.5)
                    if good > 0.001 or n == 2:
                        it = maxi
                        continue
                    v0 = np.minimum(vmax, np.maximum(vmin, v00 + ajump * (v00 - v0)))
                    ajump = -ajump * 2.0
                    state = 'S21'
                    continue
                state = 'S2'
            if state == 'S2':
                v = np.minimum(vmax, np.maximum(vmin, v0 + dv * fact))
                df0n, d = evalfit(v, cut)
                if d < d0:
                    df0 = df0n
                    v0 = v.copy()
                    if d < dbest:
                        vbest = v0.copy()
                        abest = a0.copy()
                        dbest = d
                    d0 = d
                    fact = min(1.0, fact * 4.0)
                    fact1 = 0.0
                    d1 = d
                    if not newton:
                        newton = d < 0.99 * d00
                        if newton:
                            state = 'S11'
                            continue
                    state = 'S1'
                    continue
                df0 = df0n if False else df0
                it += 1
                if it <= maxi:
                    d2 = d1
                    d1 = d
                    fact2 = fact1
                    fact1 = fact
                    fact = _fmincube(fact1, fact2, d1, d2, d0, dg)
                    if newton:
                        if fact >= factmin:
                            state = 'S2'
                            continue
                        newton = False
                        state = 'S11'
                        continue
                    if fact >= factmin1:
                        state = 'S2'
                        continue
                    v0 = np.minimum(vmax, np.maximum(vmin, 2.0 * v00 - v0))
                    state = 'S21'
                    continue
                break
        v0 = vbest
        a0 = abest
        d0 = dbest
        if cut != 0.0:
            dfc, _ = evalfit(v0, 0.0)
            for j in range(m):
                if abs(dfc[j]) > cut:
                    a0[j, :] = 0.0
        # tsvdm: A = U S V^T。共分散 = V S^-2 V^T (特異値は svdeps の相対しきい値で打ち切り)
        u, s, vt = np.linalg.svd(a0, full_matrices=False)
        w = np.where(s > svdeps * (s[0] if len(s) else 1.0), 1.0 / np.where(s > 0, s, 1.0), 0.0)
        cv = (vt.T * w ** 2) @ vt
        if n == 2:
            cv = cv * (d0 / max(1, m - nvar))
        chin = _tinvgr(float(nvar))
        ci = [math.sqrt(cv[i, i] * chin) if cv[i, i] >= 0 else float('nan') for i in range(nvar)]
        if n == 2:
            good = gammaq((m - nvar) * 0.5, (m - nvar) * 0.5)
        else:
            good = gammaq((m - nvar) * 0.5, d0 * 0.5)
        return dict(v=v0, chisq=d0, good=good, ci=ci, cov=cv, dm=abest)

    def _fit(args):
        if len(args) < 4:
            raise SADError('General::narg: Number of arguments is expected to be 4 or more in Fit')
        data = I.eval(args[0])
        from .sadscript import _HeldArg, _Lit

        def res(a):
            # 関数の引数 fun_・x_ (FitPlot[data_, fun_, x_, par0__] の Fit[data, fun, x, par0]) が束縛している式・記号に置き換える
            if isinstance(a, Sym):
                try:
                    v_ = I.lookup(a.name)
                except KeyError:
                    return a
                if isinstance(v_, _HeldArg):
                    return v_.node
                if isinstance(v_, _Lit):
                    return v_.v
                if isinstance(v_, (Sym, Expr)) and not (isinstance(v_, Expr) and v_.head == S('$Delayed')):
                    return v_
            return a
        model = res(args[1])
        xs = res(args[2])
        xnames = None
        if head_is(xs, 'List') and xs.args and all(isinstance(t, Sym) for t in xs.args):
            xnames = [t.name for t in xs.args]          # Fit[data, fun, {x1, x2}, ...]: data の行は {{x1, x2}, y [, dy]}
        elif not isinstance(xs, Sym):
            raise SADError('General::wrongtype: Symbol for indep. var as #3')
        rest = args[3:]
        opts = {}
        specs = []
        for a in rest:
            if isinstance(a, Rule) or head_is(a, 'Rule'):
                r = I.eval(a)
                opts[r.lhs.name] = r.rhs
            else:
                specs.append(a)
        names, v0, vmin, vmax = parse_vars(specs)
        xmat = None
        try:
            if xnames is not None:
                xmat = np.array([row[0] for row in data], dtype=float).reshape(len(data), len(xnames))
                arr = np.array([[float(k)] + [float(t) for t in row[1:]] for k, row in enumerate(data)], dtype=float)
            else:
                arr = np.array(data, dtype=float)
        except (TypeError, ValueError, IndexError):
            raise SADError('General::wrongval: data must be List of {x, y}, {x, y, dy}, or {x, y, dx, dy}')
        if arr.ndim != 2 or arr.shape[1] not in (2, 3, 4):
            raise SADError('General::wrongval: data must be List of {x, y}, {x, y, dy}, or {x, y, dx, dy}')
        maxi = int(opts.get('MaxIterations', 40))
        cutoff = abs(float(opts.get('Cutoff', 0.0)))
        cut = cutoff
        xn = xs.name if isinstance(xs, Sym) else 'x$'
        if arr.shape[1] <= 2 and cutoff != 0.0:
            r0 = _fit_core(arr, model, xn, names, v0.copy(), vmin, vmax, maxi, 0.0, xmat=xmat, xnames=xnames)
            cut = math.sqrt(r0['chisq'] / max(1, arr.shape[0] - len(names))) * cutoff
        if arr.shape[1] == 4 and xmat is None:
            # {x, y, dx, dy}: x の誤差を y の誤差に換算して (dy^2 + (dx f')^2)、ChiSquare が 1e-3 で収束するまで繰り返す (Fit.n)
            xs_, ys_, dx_, dy_ = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
            sig = dy_.copy()
            cs0, cs = 1e-100, 1e100
            res = None
            guard = 0
            while abs(cs / cs0 - 1.0) > 1e-3 and guard < 50:
                guard += 1
                cs0 = cs
                res = _fit_core(np.stack([xs_, ys_, sig], axis=1), model, xn, names, v0.copy(), vmin, vmax, maxi, cut)
                cs = res['chisq']
                env = {nm_: float(x_) for nm_, x_ in zip(names, res['v'])}
                df = np.zeros(len(xs_))
                for i_ in range(len(xs_)):
                    h = 1e-6 * max(abs(xs_[i_]), 1.0)
                    e1 = dict(env)
                    e1[xn] = float(xs_[i_]) + h
                    e2 = dict(env)
                    e2[xn] = float(xs_[i_]) - h
                    df[i_] = (float(ev(model, e1)) - float(ev(model, e2))) / (2.0 * h)
                sig = np.sqrt(dy_ ** 2 + (dx_ * df) ** 2)
        else:
            res = _fit_core(arr, model, xn, names, v0.copy(), vmin, vmax, maxi, cut, xmat=xmat, xnames=xnames)
        out = [Rule(S(nm), float(x)) for nm, x in zip(names, res['v'])]
        out.append(Rule(S('ChiSquare'), float(res['chisq'])))
        out.append(Rule(S('GoodnessOfFit'), float(res['good'])))
        out.append(Rule(S('ConfidenceInterval'), [float(x) for x in res['ci']]))
        out.append(Rule(S('CovarianceMatrix'), [[float(x) for x in r] for r in res['cov']]))
        out.append(Rule(S('DesignMatrix'), [[float(x) for x in r] for r in res['dm']]))
        return out
    B['Fit'] = _fit

    # ------------------------------------------------------------------ FitGaussian / FitEmit
    def _createbin(x):
        xs = sorted(float(t) for t in x)
        n = len(xs)
        nb = int(math.floor(math.sqrt(n)))
        rows = [xs[k:k + nb] for k in range(0, n - nb + 1, nb)]
        lxp = [len(r) for r in rows]
        ym = [sum(r) / len(r) for r in rows]
        yms = [math.sqrt(sum((t - ym[i]) ** 2 for t in r)) / lxp[i] for i, r in enumerate(rows)]
        xm = list(np.cumsum(lxp, dtype=float))
        sn = xm[-1]
        xm = [2.0 * (t - xm[0] / 2.0) / sn - 1.0 for t in xm]
        return np.array([xm, ym, yms]).T

    def _fitgaussian(x):
        from scipy.special import erfinv
        d = _createbin(x)
        model = Expr(S('Plus'), [Expr(S('Times'), [S('a$'), Expr(S('InverseErf'), [S('x$')])]), S('b$')])
        v0 = np.array([(d[-1, 1] - d[0, 1]) / 4.0, (d[0, 1] + d[-1, 1]) / 2.0])
        inf = np.full(2, _ALARGE)
        res = _fit_core(d, model, 'x$', ['a$', 'b$'], v0, -inf, inf)
        a, b = res['v']
        ci = res['ci']
        return [a / math.sqrt(2.0), b, [ci[0] / math.sqrt(2.0), ci[1]], res['chisq']]
    def _fitgaussian_opt(a):
        # FitGaussian[d, opt...]: Fit.n と同じく、Plot->True なら {…, Show[Plot (当てはめの曲線), ListPlot (ビンのデータ), opt]} を返す
        f = _fitgaussian(a[0])
        opts = [x for x in a[1:] if isinstance(x, Rule)]
        want = False
        for o in opts:
            if isinstance(o.lhs, Sym) and o.lhs.name == 'Plot':
                want = bool(o.rhs) and o.rhs is not Null and not (isinstance(o.rhs, Sym) and o.rhs.name == 'False')
        if not want:
            return f
        from .sadscript import _Lit
        data = [[float(v) for v in row] for row in _createbin(a[0])]
        curve = Expr(S('Plot'), [
            Expr(S('Plus'), [Expr(S('Times'), [Expr(S('InverseErf'), [S('x')]), f[0] * math.sqrt(2.0)]), f[1]]),
            Expr(S('List'), [S('x'), -0.9999, 0.9999]),
            Rule(S('DisplayFunction'), S('Identity'))])
        g1 = I.eval(curve)
        g2 = I.eval(Expr(S('ListPlot'), [_Lit(data), Rule(S('DisplayFunction'), S('Identity'))]))
        shown = I.eval(Expr(S('Show'), [_Lit(g1), _Lit(g2)] + [_Lit(o) for o in opts]))
        return f + [shown]
    B['FitGaussian'] = _fitgaussian_opt

    def _fitemit(a):
        x, px = list(a[0]), list(a[1])
        sigx, xm, cx, rx = _fitgaussian(x)
        sigpx, pxm, cpx, rpx = _fitgaussian(px)
        x1 = [(t - xm) / sigx for t in x]
        px1 = [(t - pxm) / sigpx for t in px]
        sigx1, xm1, cx1, rx1 = _fitgaussian([p + q for p, q in zip(x1, px1)])
        sigpx1, pxm1, cpx1, rpx1 = _fitgaussian([p - q for p, q in zip(x1, px1)])
        em = sigx1 * sigpx1 * sigx * sigpx / 2.0
        beta = sigx ** 2 / em
        dsigx1 = cx1[0] / sigx1
        dsigpx1 = cpx1[0] / sigpx1
        demr = math.sqrt((cx[0] / sigx) ** 2 + (cpx[0] / sigpx) ** 2 + dsigx1 ** 2 + dsigpx1 ** 2)
        sr = sigpx1 / sigx1
        return [[xm, pxm, (sr - 1 / sr) / 2.0, beta, em],
                [cx[1], cpx[1], (sr + 1 / sr) * math.sqrt(dsigx1 ** 2 + dsigpx1 ** 2) / 4.0, beta * demr, demr * em]]
    B['FitEmit'] = _fitemit

    # ------------------------------------------------------------------ NIntegrate (Packages/D.n)
    def _nintegrate(args):
        if len(args) < 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 or more in NIntegrate')
        f = args[0]
        spec = args[1]
        if not head_is(spec, 'List') or len(spec.args) != 3 or not isinstance(spec.args[0], Sym):
            return Expr(S('NIntegrate'), list(args))
        xn = spec.args[0].name
        x0 = float(I.eval(spec.args[1]))
        x1 = float(I.eval(spec.args[2]))
        opts, _ = get_opts(args[2:])
        eps = float(opts.get('AccuracyGoal', 1e-13))
        n = int(opts.get('InitialPoints', 20))
        spl = opts.get('Spline', 1.0)
        spl = bool(spl) and spl is not False

        def fa(u):
            out = []
            for t in u:
                try:
                    val = ev(f, {xn: float(t)})
                except SADError:
                    val = 0.0
                out.append(float(val) if is_num(val) else 0.0)
            return np.array(out)

        def sp_int(u, v):
            from scipy.interpolate import CubicSpline
            return float(CubicSpline(u, v).integrate(u[0], u[-1]))

        d = (x1 - x0) / n
        u = np.arange(n + 1) * d + x0
        v = fa(u)
        smax = float(np.sum(np.abs(v))) - 0.5 * (abs(v[0]) + abs(v[-1]))
        s = sp_int(u, v) if spl else d * (float(np.sum(v)) - 0.5 * (v[0] + v[-1]))
        while True:
            d *= 0.5
            u1 = u + d
            u1[-1] = x1
            v1 = fa(u1)
            vv = np.empty(2 * len(v))
            vv[0::2] = v
            vv[1::2] = v1
            v = vv[:-1]
            n *= 2
            u = np.arange(n + 1) * d + x0
            if spl:
                s1 = sp_int(u, v)
            else:
                s1 = 0.5 * s + d * float(np.sum(v1[:-1]))
            ds = abs(s1 - s)
            s = s1
            smax += float(np.sum(np.abs(v1)))
            if not ds >= smax * abs(d * eps):
                break
            if n >= 20000:
                print('???NIntegrate::conv:  Poor Convergence ' + repr(float(ds / (smax * abs(d)))), file=I.out)
                break
        return float(s)
    B['NIntegrate'] = _nintegrate
