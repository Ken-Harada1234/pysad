"""
Vector / ReleaseVector (K64 のベクトル評価)。`Vector[{x1,x2,..}]` はベクトル値の引数で、その上の算術・初等関数は成分ごとに働き、
結果も Vector になる。`ReleaseVector[v]` は成分のリストを返す (Plots.n の PlotInterpolate が関数を一度にたくさんの点で評価するのに使う)。

   ReleaseVector[Sin[Vector[{0.,1.,2.}]]]   →  {0, .8414709848078965, .9092974268256817}
   ReleaseVector[Vector[{1,2}]+3]             →  {4, 5}

pysad は、成分ごとの評価を、下の初等関数・算術の組み込みを包んで実現する (Vector を含む呼び出しだけ成分ごとに回す)。
"""
from __future__ import annotations

from .sadscript import Expr, S, Null, is_num

_FUNCS = ('Plus', 'Times', 'Power', 'Divide', 'Subtract', 'Minus', 'Sqrt', 'Exp', 'Log', 'Sin', 'Cos', 'Tan', 'ArcSin', 'ArcCos',
          'ArcTan', 'Sinh', 'Cosh', 'Tanh', 'ArcSinh', 'ArcCosh', 'ArcTanh', 'Abs', 'Re', 'Im', 'Conjugate', 'Arg', 'Sign', 'Round',
          'Floor', 'Ceiling', 'Erf', 'Erfc', 'InverseErf', 'InverseErfc', 'ExpIntegralE', 'ExpIntegralEi', 'ExpIntegralEi', 'Gamma', 'LogGamma', 'BesselJ', 'BesselY', 'BesselI', 'BesselK', 'Sinc', 'Mod', 'Max',
          'Min', 'Less', 'LessEqual', 'Greater', 'GreaterEqual', 'Equal', 'Unequal', 'And', 'Or', 'Not', 'If')


def _isvec(x):
    return isinstance(x, Expr) and x.head == S('Vector') and len(x.args) == 1 and isinstance(x.args[0], list)


def install(I):
    B = I.builtins

    def vector(a):
        if len(a) == 1 and isinstance(a[0], list):
            return Expr(S('Vector'), [a[0]])
        return Expr(S('Vector'), list(a))

    def release(a):
        x = a[0] if a else Null
        return list(x.args[0]) if _isvec(x) else x

    B['Vector'] = vector
    B['ReleaseVector'] = release

    def wrap(name, f):
        def w(args):
            vs = [x for x in args if _isvec(x)]
            if not vs:
                return f(args)
            n = max(len(v.args[0]) for v in vs)
            if any(len(v.args[0]) != n for v in vs):
                return f(args)
            out = []
            for k in range(n):
                out.append(f([x.args[0][k] if _isvec(x) else x for x in args]))
            return Expr(S('Vector'), [out])
        return w

    for nm in _FUNCS:
        f = B.get(nm)
        if f is not None and nm not in ('If', 'And', 'Or'):
            B[nm] = wrap(nm, f)
