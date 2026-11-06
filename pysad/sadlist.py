"""
sadlist.py — SADScript のリスト・式操作関数 (SAD の Mathematica 互換実装を直訳した挙動)

sadscript.Interp._install_builtins の最後で install(interp) が呼ばれ、同名の既存関数を上書きする。
挙動は SAD K64 (bin/gs) で bench_level24a*.sad を実行して得た参照出力に合わせている (docs/history/work_item4.md 参照)。

SAD の定義の出所:
  - 組み込み (src/*.f): MapThread, MapAt, Complement, DeleteCases, ScanThread, Cases, MatchQ, Extract,
    Intersection, Union (Union1), Position, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map など
  - Packages/functions.n: Count[x__]:=Length[Position[x]]、MemberQ=Position[l,p,1,1]<=>{}、
    FreeQ=Position[f,p,{0,Infinity},1]==={}、Nest/NestList/Fold/FoldList (Do ベース)、Through、
    Positive/Negative/IntegerQ、FractionalPart[x]=x-Floor[x]、RotateLeft/Right、Union[x__]=Union1[Join[x]]
  - Packages/Misc.n: FixedPoint (既定の SameTest は NearlySameQ[a,b] = $NearlySameQ[a,b,1e-12,1e-20])
  - src/tfsameq.f: tfnearlysameqf (|a-b| <= (|a|+|b|)*re または <= ae)

値の表現 (sadscript.py と同じ): 実数は float、リストは list、未評価式は Expr(head, args)、文字列は str。
"""
from __future__ import annotations

import functools
import math
import re

import numpy as np

from .sadscript import (Sym, Expr, Rule, Func, Seq, SADError, S, Null, is_num, truthy, _hkey, INF, NAN, _Lit)

_LIST = S('List')


# ---------------------------------------------------------------------------- 基本部品
def _isseq(x):
    """リストまたは複合式 (Expr)"""
    return isinstance(x, (list, Expr))


def _elems(x):
    return x if isinstance(x, list) else x.args


def _rebuild(x, items):
    return list(items) if isinstance(x, list) else Expr(x.head, list(items))


def _head_of(x):
    return _LIST if isinstance(x, list) else x.head


def _int(x, what='number'):
    if not is_num(x) or isinstance(x, complex):
        raise SADError(f'General::wrongtype: Argument must be Real number ({what})')
    return int(math.floor(x)) if x >= 0 else -int(math.floor(-x))


def _isint_num(x):
    return is_num(x) and not isinstance(x, complex)


def _norm_idx(n, length, what='Index'):
    """1 始まり・負は末尾から。範囲外は SADError"""
    k = _int(n)
    if k < 0:
        k = length + k + 1
    if not 1 <= k <= length:
        raise SADError(f'General::index: {what} {int(_int(n))} out of range')
    return k


def _depth(x):
    """Depth: 原子は 1、リスト・複合式は 1+最大の子 (空リストは 2)"""
    if isinstance(x, list):
        return 1 + max((_depth(y) for y in x), default=1)
    if isinstance(x, Expr):
        return 1 + max((_depth(y) for y in x.args), default=1)
    return 1


def _lspec(v, default):
    """レベル指定 → (lo, hi)。n → (1, n)、{n} → (n, n)、{m, n} → (m, n)、Infinity → (1, inf)"""
    if v is None:
        return default
    if isinstance(v, list):
        if len(v) == 1:
            return (_lv(v[0]), _lv(v[0]))
        if len(v) == 2:
            return (_lv(v[0]), _lv(v[1]))
        raise SADError('General::wrongtype: Argument must be level specification')
    return (1, _lv(v))


def _lv(x):
    if is_num(x):
        return x if math.isinf(x) else int(x)
    raise SADError('General::wrongtype: Argument must be level specification')


def _in_level(lo, hi, d, x):
    """深さ d の部分式 x が [lo, hi] に入るか (負は「高さ」で数える)"""
    h = None
    if lo >= 0:
        if d < lo:
            return False
    else:
        h = _depth(x)
        if h > -lo:
            return False
    if hi >= 0:
        if d > hi:
            return False
    else:
        if h is None:
            h = _depth(x)
        if h < -hi:
            return False
    return True


def _walk(x, lo, hi, heads=False, path=None, d=0):
    """部分式を (path, sub, depth) で列挙。子 → 親の順 (SAD の Level/Cases/Position と同じ後順)。
    heads=True なら頭部 (path に 0 を付ける) も、子より先に列挙する (SAD の Position)"""
    path = path or []
    if _isseq(x) and (hi < 0 or d < hi):
        if heads:
            hd = _head_of(x)
            yield from _walk(hd, lo, hi, heads, path + [0.0], d + 1)
        for k, c in enumerate(_elems(x)):
            yield from _walk(c, lo, hi, heads, path + [float(k + 1)], d + 1)
    if _in_level(lo, hi, d, x):
        yield path, x, d


# ---------------------------------------------------------------------------- 標準順序 (Sort / Union / Order)
def _chars_key(s):
    # SAD: a < A < aa < b < B (小文字が同じ大文字より先、先頭文字優先)
    return tuple((c.lower(), c.isupper()) for c in s)


def ckey(x):
    """SAD の標準順序: 数 < 文字列 < シンボル < 複合式 < リスト (リストは長さ優先)"""
    if is_num(x):
        if isinstance(x, complex):
            return (0, x.real, x.imag)
        return (0, x, 0.0)
    if isinstance(x, str):
        return (1, _chars_key(x))
    if isinstance(x, Sym):
        return (2, _chars_key(x.name))
    if isinstance(x, Rule):
        return (3, _chars_key('Rule'), 2, (ckey(x.lhs), ckey(x.rhs)))
    if isinstance(x, Expr):
        hk = ckey(x.head)
        return (3, hk, len(x.args), tuple(ckey(y) for y in x.args))
    if isinstance(x, list):
        return (4, len(x), tuple(ckey(y) for y in x))
    return (5, repr(x))


def _uniq_sorted(items):
    items = sorted(items, key=ckey)
    out = []
    last = None
    for x in items:
        k = _hkey(x)
        if out and k == last:
            continue
        out.append(x)
        last = k
    return out


# ---------------------------------------------------------------------------- 本体
def install(I):
    B = I.builtins

    def apply(f, args):
        return I.apply_value(f, list(args), f)

    def mq(p, x):
        return I._match_one(p, x, {})

    def need_list(x, what='List'):
        if not isinstance(x, list):
            raise SADError(f'General::wrongtype: Argument must be {what}')

    def tf(b):
        return 1.0 if b else 0.0

    # ======================================================================== 位置指定の解釈
    def positions(spec, allow_zero=False):
        """Part の位置指定 → パスのリスト。n → [[n]]、{i,j,...} → [[i,j,...]]、{{..},{..}} → 複数"""
        if is_num(spec):
            return [[spec]]
        if isinstance(spec, list):
            if spec and all(isinstance(s, list) for s in spec):
                return [list(s) for s in spec]
            if all(is_num(s) for s in spec):
                return [list(spec)]
        raise SADError('General::wrongtype: Argument must be Real, List of Reals, or List of Lists')

    def get_at(x, path):
        v = x
        for k in path:
            if not _isseq(v):
                raise SADError('General::index: Part specification out of range')
            el = _elems(v)
            v = el[_norm_idx(k, len(el)) - 1]
        return v

    def modify_at(x, path, fn, insert=False, delete=False):
        """path の位置の要素を fn(old) で置換した新しい式 (insert: その位置の前に fn(None) を挿入、delete: 削除)"""
        if not path:
            return fn(x)
        if not _isseq(x):
            raise SADError('General::index: Part specification out of range')
        el = list(_elems(x))
        n = len(el)
        k = path[0]
        if len(path) == 1 and insert:
            kk = _int(k)
            if kk < 0:
                kk = n + kk + 1
            if not 1 <= kk <= n + 1:
                raise SADError(f'General::index: Index {int(_int(k))} out of range')
            el.insert(kk - 1, fn(None))
            return _rebuild(x, el)
        i = _norm_idx(k, n) - 1
        if len(path) == 1 and delete:
            del el[i]
            return _rebuild(x, el)
        el[i] = modify_at(el[i], path[1:], fn, insert, delete)
        return _rebuild(x, el)

    # ======================================================================== MapThread / ScanThread / MapAt
    def _thread_lists(a, name):
        if len(a) != 2:
            raise SADError(f'General::narg: Number of arguments is expected to be 2 in {name}')
        ls = a[1]
        need_list(ls)
        if not all(isinstance(l, list) for l in ls):
            raise SADError('General::wrongtype: Argument must be List of Lists')
        if ls and len({len(l) for l in ls}) != 1:
            raise SADError(f'General::equalleng: Equal lengths of elements in Thread are expected in {name}')
        return ls

    def _mapthread(a):
        ls = _thread_lists(a, 'MapThread')
        if not ls:
            return []
        return [apply(a[0], [l[k] for l in ls]) for k in range(len(ls[0]))]
    B['MapThread'] = _mapthread

    def _scanthread(a):
        ls = _thread_lists(a, 'ScanThread')
        if ls:
            for k in range(len(ls[0])):
                apply(a[0], [l[k] for l in ls])
        return Null
    B['ScanThread'] = _scanthread

    def _mapat(a):
        if len(a) != 3:
            raise SADError('General::narg: Number of arguments is expected to be 3 in MapAt')
        f, x, spec = a
        for path in positions(spec):
            x = modify_at(x, path, lambda old: apply(f, [old]))
        return x
    B['MapAt'] = _mapat

    # ======================================================================== 集合演算
    def _lists_only(a, name):
        for x in a:
            if not isinstance(x, list):
                raise SADError(f'General::wrongtype: Argument must be List in {name}')

    def _union(a):
        same = None
        if a and isinstance(a[-1], Rule) and isinstance(a[-1].lhs, Sym) and a[-1].lhs.name == 'SameTest':
            same = a[-1].rhs
            a = a[:-1]
        if not a:
            return Expr(S('Union'), [])
        _lists_only(a, 'Union')
        items = sorted([y for x in a for y in x], key=ckey)
        if same is None:
            return _uniq_sorted(items)
        out = []
        for y in items:
            if out and truthy(apply(same, [out[-1], y])) is True:
                continue
            out.append(y)
        return out
    B['Union'] = _union

    def _intersection(a):
        if not a:
            return Expr(S('Intersection'), [])
        _lists_only(a, 'Intersection')
        keys = [{_hkey(y) for y in x} for x in a[1:]]
        return _uniq_sorted([y for y in a[0] if all(_hkey(y) in k for k in keys)])
    B['Intersection'] = _intersection

    def _complement(a):
        if not a:
            return Expr(S('Complement'), [])
        _lists_only(a, 'Complement')
        keys = [{_hkey(y) for y in x} for x in a[1:]]
        return _uniq_sorted([y for y in a[0] if not any(_hkey(y) in k for k in keys)])
    B['Complement'] = _complement

    # ======================================================================== パターン照合系
    def _rule_apply(rule, x):
        """Rule を部分式 x 全体に適用。一致しなければ _NO"""
        b = {}
        if not I._match_one(rule.lhs, x, b):
            return _NO
        rhs = rule.rhs
        return I.eval(I.subst_bindings(rhs, b)) if b else I.eval(rhs)

    _NO = object()

    def _cases(a):
        if len(a) < 2 or len(a) > 4:
            raise SADError('General::narg: Number of arguments is expected to be 2, 3 or 4 in Cases')
        x, p = a[0], a[1]
        lo, hi = _lspec(a[2] if len(a) > 2 else None, (1, 1))
        nmax = int(a[3]) if len(a) > 3 and not math.isinf(a[3]) else None
        out = []
        for _, sub, _ in _walk(x, lo, hi):
            if nmax is not None and len(out) >= nmax:
                break
            if isinstance(p, Rule):
                r = _rule_apply(p, sub)
                if r is not _NO:
                    out.append(r)
            elif mq(p, sub):
                out.append(sub)
        return out
    B['Cases'] = _cases

    def _deletecases(a):
        if len(a) < 2 or len(a) > 4:
            raise SADError('General::narg: Number of arguments is expected to be 2, 3 or 4 in DeleteCases')
        x, p = a[0], a[1]
        lo, hi = _lspec(a[2] if len(a) > 2 else None, (1, 1))
        nmax = int(a[3]) if len(a) > 3 and not math.isinf(a[3]) else None
        left = [nmax if nmax is not None else float('inf')]

        def rec(y, d):
            if not _isseq(y):
                return y
            keep = []
            for c in _elems(y):
                c2 = rec(c, d + 1) if (hi < 0 or d + 1 < hi) else c
                if left[0] > 0 and _in_level(lo, hi, d + 1, c) and mq(p, c):
                    left[0] -= 1
                    continue
                keep.append(c2)
            return _rebuild(y, keep)
        return rec(x, 0)
    B['DeleteCases'] = _deletecases

    def _matchq(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in MatchQ')
        return tf(mq(a[1], a[0]))
    B['MatchQ'] = _matchq

    def _position(a):
        if len(a) < 2 or len(a) > 4:
            raise SADError('General::narg: Number of arguments is expected to be 2, 3 or 4 in Position')
        x, p = a[0], a[1]
        lo, hi = _lspec(a[2] if len(a) > 2 else None, (0, INF))
        nmax = int(a[3]) if len(a) > 3 and not math.isinf(a[3]) else None
        out = []
        for path, sub, _ in _walk(x, lo, hi, heads=True):
            if nmax is not None and len(out) >= nmax:
                break
            if mq(p, sub):
                out.append(list(path))
        return out
    B['Position'] = _position
    B['Count'] = lambda a: float(len(_position(a)))
    B['MemberQ'] = lambda a: tf(len(_position([a[0], a[1], a[2] if len(a) > 2 else 1.0, 1.0])) > 0)
    B['FreeQ'] = lambda a: tf(len(_position([a[0], a[1], a[2] if len(a) > 2 else [0.0, INF], 1.0])) == 0)

    def _extract(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Extract')
        x, spec = a[0], a[1]
        if not isinstance(spec, list):
            raise SADError('General::wrongtype: Argument must be List for #2 in Extract')
        h = a[2] if len(a) > 2 else None

        def one(path):
            v = get_at(x, path)
            return apply(h, [v]) if h is not None else v
        if spec and all(isinstance(s, list) for s in spec):
            return [one(p) for p in spec]
        return one(spec)
    B['Extract'] = _extract

    def _select(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Select')
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Select')
        n = int(a[2]) if len(a) > 2 and not math.isinf(a[2]) else None
        out = []
        for y in _elems(a[0]):
            if n is not None and len(out) >= n:
                break
            if truthy(apply(a[1], [y])) is True:
                out.append(y)
        return out
    B['Select'] = _select

    # ======================================================================== 反復
    def _nest(a):
        f, x, n = a[0], a[1], _int(a[2])
        for _ in range(max(n, 0)):
            x = apply(f, [x])
        return x
    B['Nest'] = _nest

    def _nestlist(a):
        f, x, n = a[0], a[1], _int(a[2])
        out = [x]
        for _ in range(max(n, 0)):
            x = apply(f, [x])
            out.append(x)
        return out
    B['NestList'] = _nestlist

    def _fold_args(a, name):
        if len(a) != 3:
            return None
        return a

    def _fold(a):
        if len(a) != 3:
            return Expr(S('Fold'), a)
        y = a[1]
        for e in _elems(a[2]) if _isseq(a[2]) else []:
            y = apply(a[0], [y, e])
        return y
    B['Fold'] = _fold

    def _foldlist(a):
        if len(a) != 3:
            return Expr(S('FoldList'), a)
        y = a[1]
        out = [y]
        for e in _elems(a[2]) if _isseq(a[2]) else []:
            y = apply(a[0], [y, e])
            out.append(y)
        return out
    B['FoldList'] = _foldlist

    def nearly_same(x, y, re_=1e-12, ae=1e-20):
        """tfsameq.f: tfnearlysameqf"""
        if is_num(x):
            if is_num(y):
                s = abs(x) + abs(y)
                d = abs(y - x)
                return d <= s * re_ or d <= ae
            return False
        if is_num(y):
            return False
        if not isinstance(x, list) or not isinstance(y, list):
            return I.same(x, y) if not (isinstance(x, list) or isinstance(y, list)) else False
        if len(x) != len(y):
            return False
        return all(nearly_same(p, q, re_, ae) for p, q in zip(x, y))

    def _nearlysame(a):
        if len(a) == 2:
            return tf(nearly_same(a[0], a[1]))
        thr, ath = 1e-12, 1e-20
        for o in a[2:]:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym):
                if o.lhs.name == 'Threshold':
                    thr = o.rhs
                elif o.lhs.name == 'AbsoluteThreshold':
                    ath = o.rhs
        return tf(nearly_same(a[0], a[1], thr, ath))
    B['NearlySameQ'] = _nearlysame
    B['$NearlySameQ'] = lambda a: tf(nearly_same(a[0], a[1], a[2], a[3]))

    def _fixedpoint(a):
        """Packages/Misc.n: FixedPoint[f,e,n,opt] = Table[e1=f[e0]; If[NearlySameQ[e1,e0],Break[]]; e0=e1,{n}]; e1"""
        f, e0 = a[0], a[1]
        n = INF
        opts = []
        rest = list(a[2:])
        if rest and is_num(rest[0]):
            n = rest.pop(0)
        opts = rest
        same = None
        thr, ath = 1e-12, 1e-20
        for o in opts:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym):
                if o.lhs.name == 'SameTest':
                    same = o.rhs
                elif o.lhs.name == 'Threshold':
                    thr = o.rhs
                elif o.lhs.name == 'AbsoluteThreshold':
                    ath = o.rhs
        if same is not None and isinstance(same, Sym) and same.name == 'NearlySameQ':
            same = None
        e1 = None
        k = 0
        limit = 10 ** 6 if math.isinf(n) else int(n)
        while k < limit:
            e1 = apply(f, [e0])
            if same is None:
                if nearly_same(e1, e0, thr, ath):
                    break
            elif truthy(apply(same, [e1, e0])) is True:
                break
            e0 = e1
            k += 1
        return e1 if e1 is not None else e0
    B['FixedPoint'] = _fixedpoint

    # ======================================================================== 切り出し・挿入・削除
    def _take_range(x, spec, drop=False):
        el = list(_elems(x))
        n = len(el)
        if isinstance(spec, list):
            if len(spec) == 1:
                i = j = _int(spec[0])
            elif len(spec) == 2:
                i, j = _int(spec[0]), _int(spec[1])
            else:
                raise SADError('General::wrongleng: Length of range-list must be less than 3')
            if i < 0:
                i = n + i + 1
            if j < 0:
                j = n + j + 1
            if drop:
                if len(spec) == 1 and not 1 <= i <= n:
                    raise SADError('General::index: Index out of range in Drop')
                return el[:max(i - 1, 0)] + el[max(j, 0):] if j >= i else el
            if j > n or i < 1 or j < i:
                return []
            return el[i - 1:j]
        k = _int(spec)
        if drop:
            if k >= 0:
                return el[k:]
            return el[:n + k] if -k <= n else []
        if abs(k) > n:
            return []
        return el[:k] if k >= 0 else el[n + k:]

    def _take(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Take')
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Take')
        return _rebuild(a[0], _take_range(a[0], a[1]))
    B['Take'] = _take

    def _drop(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Drop')
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Drop')
        return _rebuild(a[0], _take_range(a[0], a[1], drop=True))
    B['Drop'] = _drop

    def _delete(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Delete')
        x = a[0]
        ps = positions(a[1])
        for p in ps:
            if any(_int(k) == 0 for k in p):
                raise SADError('General::index: Index 0 is not allowed in Delete')
        # 元の位置で指定されるので、後ろ側から消す
        for p in sorted(ps, key=lambda q: _resolve(x, q), reverse=True):
            x = modify_at(x, p, None, delete=True)
        return x

    def _resolve(x, path, ins=False):
        """負の添字を正に直したパス (並べ替え用)"""
        out = []
        v = x
        for k in path:
            el = _elems(v) if _isseq(v) else []
            kk = _int(k)
            kk = len(el) + kk + 1 if kk < 0 else kk
            out.append(kk)
            v = el[kk - 1] if 1 <= kk <= len(el) else None
        return out
    B['Delete'] = _delete

    def _insert(a):
        if len(a) != 3:
            raise SADError('General::narg: Number of arguments is expected to be 3 in Insert')
        x, new = a[0], a[1]
        ps = positions(a[2])
        for p in sorted(ps, key=lambda q: _resolve(x, q, True), reverse=True):
            x = modify_at(x, p, lambda old: new, insert=True)
        return x

    B['Insert'] = _insert

    def _replacepart(a):
        if len(a) == 2:
            rules = a[1] if isinstance(a[1], list) else [a[1]]
            x = a[0]
            for r in rules:
                if not isinstance(r, Rule):
                    raise SADError('General::wrongtype: Argument must be Rule or List of Rules in ReplacePart')
                for p in positions(r.lhs):
                    x = modify_at(x, p, lambda old, r=r: r.rhs)
            return x
        if len(a) != 3:
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in ReplacePart')
        x, new, spec = a
        for p in positions(spec):
            x = modify_at(x, p, lambda old: new)
        return x
    B['ReplacePart'] = _replacepart

    def _append(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Append')
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Append')
        return _rebuild(a[0], list(_elems(a[0])) + [a[1]])
    B['Append'] = _append

    def _prepend(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Prepend')
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Prepend')
        return _rebuild(a[0], [a[1]] + list(_elems(a[0])))
    B['Prepend'] = _prepend

    # ======================================================================== Apply / Map
    def _apply(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Apply')
        f, x = a[0], a[1]
        lo, hi = _lspec(a[2] if len(a) > 2 else None, (0, 0))

        def rec(y, d):
            if not _isseq(y):
                return y
            if hi < 0 or d < hi:
                y = _rebuild(y, [rec(c, d + 1) for c in _elems(y)])
            if _in_level(lo, hi, d, y):
                return apply(f, _elems(y))
            return y
        return rec(x, 0)
    B['Apply'] = _apply

    def _map(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Map')
        f, x = a[0], a[1]
        lo, hi = _lspec(a[2] if len(a) > 2 else None, (1, 1))

        def rec(y, d):
            if _isseq(y) and (hi < 0 or d < hi):
                y = _rebuild(y, [rec(c, d + 1) for c in _elems(y)])
            if _in_level(lo, hi, d, y):
                return apply(f, [y])
            return y
        return rec(x, 0)
    B['Map'] = _map
    B['MapAll'] = lambda a: _map([a[0], a[1], [0.0, INF]])

    def _mapindexed(a):
        f, x = a[0], a[1]
        if not _isseq(x):
            return x
        return _rebuild(x, [apply(f, [c, [float(k + 1)]]) for k, c in enumerate(_elems(x))])
    B['MapIndexed'] = _mapindexed

    def _scan(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Scan')
        f, x = a[0], a[1]
        if _isseq(x):
            lo, hi = _lspec(a[2] if len(a) > 2 else None, (1, 1))
            for _, sub, _ in _walk(x, lo, hi):
                apply(f, [sub])
        return Null
    B['Scan'] = _scan

    def _through(a):
        x = a[0] if a else None
        if isinstance(x, Expr) and isinstance(x.head, Expr) and isinstance(x.head.head, Sym):
            p = x.head.head
            return apply(p, [apply(g, list(x.args)) for g in x.head.args])
        if isinstance(x, Expr) and isinstance(x.head, list):
            raise SADError('General::wrongtype: Argument must be Real, List of Reals, or Null as index')
        return Expr(S('Through'), a)
    B['Through'] = _through

    # ======================================================================== リスト整形
    def _reverse(a):
        if not _isseq(a[0]):
            raise SADError('General::wrongtype: Argument must be List or composition in Reverse')
        return _rebuild(a[0], list(reversed(_elems(a[0]))))
    B['Reverse'] = _reverse

    def _first_last(name, pick):
        def g(a):
            x = a[0]
            if not _isseq(x) or not _elems(x):
                raise SADError(f'General::index: Index out of range in {name}')
            return pick(_elems(x))
        return g
    B['First'] = _first_last('First', lambda e: e[0])
    B['Last'] = _first_last('Last', lambda e: e[-1])

    def _rest(a):
        x = a[0]
        if not _isseq(x) or not _elems(x):
            raise SADError('General::index: Index out of range in Rest')
        return _rebuild(x, _elems(x)[1:])
    B['Rest'] = _rest

    def _join(a):
        if not a:
            raise SADError('General::wrongtype: Argument must be List or composition for all args in Join[]')
        if not all(_isseq(x) for x in a):
            raise SADError('General::wrongtype: Argument must be List or composition for all args in Join')
        h = _head_of(a[0])
        if not all(I.same(_head_of(x), h) for x in a):
            raise SADError('General::samehead: Same heads are expected in Join')
        return _rebuild(a[0], [y for x in a for y in _elems(x)])
    B['Join'] = _join

    def _flatten(a):
        if len(a) not in (1, 2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 1, 2 or 3 in Flatten')
        x = a[0]
        if not _isseq(x):
            raise SADError('General::wrongtype: Argument must be List or composition in Flatten')
        depth = INF
        if len(a) > 1:
            if not _isint_num(a[1]) or a[1] <= 0:
                raise SADError('General::wrongnum: Number must be Positive in Flatten')
            depth = a[1]
        h = a[2] if len(a) > 2 else _head_of(x)

        def rec(y, d):
            out = []
            for c in _elems(y):
                if _isseq(c) and d > 0 and I.same(_head_of(c), h):
                    out.extend(rec(c, d - 1))
                else:
                    out.append(c)
            return out
        return _rebuild(x, rec(x, int(depth) if not math.isinf(depth) else 10 ** 6))
    B['Flatten'] = _flatten

    def _partition(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 or 3 in Partition')
        x = a[0]
        if not _isseq(x):
            raise SADError('General::wrongtype: Argument must be List or composition in Partition')
        n = a[1][0] if isinstance(a[1], list) and len(a[1]) == 1 else a[1]
        n = _int(n)
        d = _int(a[2][0] if isinstance(a[2], list) and len(a[2]) == 1 else a[2]) if len(a) > 2 else n
        if n <= 0 or d <= 0:
            raise SADError('General::wrongnum: Number must be Positive in Partition')
        el = list(_elems(x))
        return [_rebuild(x, el[k:k + n]) for k in range(0, len(el) - n + 1, d)]
    B['Partition'] = _partition

    def _transpose(a):
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in Transpose')
        m = a[0]
        if not isinstance(m, list) or not m or not all(isinstance(r, list) for r in m):
            raise SADError('General::wrongtype: Argument must be Matrix in Transpose')
        if len({len(r) for r in m}) != 1:
            raise SADError('General::equalleng: Equal lengths of rows of matrix are expected in Transpose')
        return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]
    B['Transpose'] = _transpose

    def _rot(x, n):
        el = list(_elems(x))
        if not el:
            return _rebuild(x, el)
        k = _int(n) % len(el)
        return _rebuild(x, el[-k:] + el[:-k] if k else el)

    def _rotright(a, sign=1):
        """RotateRight[a,n] / RotateRight[a,{n1,n2,...}] (functions.n: 先頭の n1 で回してから各要素を残りで回す)"""
        x = a[0]
        if not _isseq(x):
            raise SADError('General::wrongtype: Argument must be List or composition in RotateRight')
        n = a[1] if len(a) > 1 else 1.0
        if isinstance(n, list):
            if not n:
                return x
            y = _rot(x, sign * n[0])
            if len(n) > 1:
                return _rebuild(y, [_rotright([c, list(n[1:])], sign) for c in _elems(y)])
            return y
        return _rot(x, sign * n)
    B['RotateRight'] = lambda a: _rotright(a, 1)
    B['RotateLeft'] = lambda a: _rotright(a, -1)

    def _sort(a):
        if len(a) not in (1, 2):
            raise SADError('General::narg: Number of arguments is expected to be 1 or 2 in Sort')
        x = a[0]
        if not _isseq(x):
            raise SADError('General::wrongtype: Argument must be List or composition in Sort')
        el = list(_elems(x))
        if len(a) == 1:
            return _rebuild(x, sorted(el, key=ckey))
        f = a[1]

        def cmp(p, q):
            r = apply(f, [p, q])
            if truthy(r) is True:
                return -1
            if truthy(r) is None:
                raise SADError('General::wrongval: Real number must be as the result of order-function in Sort')
            return 1 if truthy(apply(f, [q, p])) is True else 0
        return _rebuild(x, sorted(el, key=functools.cmp_to_key(cmp)))
    B['Sort'] = _sort

    def _order(a):
        k1, k2 = ckey(a[0]), ckey(a[1])
        return 0.0 if _hkey(a[0]) == _hkey(a[1]) else (1.0 if k1 < k2 else -1.0)
    B['Order'] = _order

    # ======================================================================== 基本情報
    def _length(a):
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in Length')
        x = a[0]
        return float(len(_elems(x))) if _isseq(x) else 0.0
    B['Length'] = _length

    def _dimensions(a):
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in Dimensions')
        x = a[0]
        if not _isseq(x):
            return []
        dims = [len(_elems(x))]
        cur = list(_elems(x))
        while cur:
            if not all(isinstance(c, list) for c in cur) or len({len(c) for c in cur}) != 1:
                break
            n = len(cur[0])
            dims.append(n)
            cur = [y for c in cur for y in c]
            if n == 0:
                break
        return [float(d) for d in dims]
    B['Dimensions'] = _dimensions
    B['Depth'] = lambda a: float(_depth(a[0]))

    def _level(a):
        if len(a) not in (2, 3):
            raise SADError('General::narg: Number of arguments is expected to be 2 in Level')
        lo, hi = _lspec(a[1], (1, 1))
        return [sub for _, sub, _ in _walk(a[0], lo, hi)]
    B['Level'] = _level

    B['AtomQ'] = lambda a: tf(not _isseq(a[0]))

    def _listwise(name, f):
        def g(a):
            x = a[0]
            if isinstance(x, list):
                return [g([y]) for y in x]
            if _isint_num(x):
                return tf(f(x))
            return Expr(S(name), a)
        return g
    B['IntegerQ'] = _listwise('IntegerQ', lambda x: x == math.floor(x))
    B['Positive'] = _listwise('Positive', lambda x: x > 0)
    B['Negative'] = _listwise('Negative', lambda x: x < 0)
    B['NonNegative'] = _listwise('NonNegative', lambda x: x >= 0)
    B['NonPositive'] = _listwise('NonPositive', lambda x: x <= 0)

    def _evenodd(name, par):
        def g(a):
            x = a[0]
            if _isint_num(x) and x == math.floor(x):
                return tf(int(x) % 2 == par)
            if _isint_num(x):
                return 0.0
            return Expr(S(name), a)
        return g
    B['EvenQ'] = _evenodd('EvenQ', 0)
    B['OddQ'] = _evenodd('OddQ', 1)

    def _binomial(a):
        n, k = a
        if _isint_num(n) and _isint_num(k) and n == math.floor(n) and k == math.floor(k):
            if k < 0 or k > n:
                return 0.0
            return float(math.comb(int(n), int(k)))
        return math.gamma(n + 1) / math.gamma(k + 1) / math.gamma(n - k + 1)
    B['Binomial'] = _binomial

    # FractionalPart[x]=x-Floor[x] (functions.n)
    def _fracpart(a):
        x = a[0]
        if isinstance(x, list):
            return [_fracpart([y]) for y in x]
        if _isint_num(x):
            return x - math.floor(x)
        return Expr(S('FractionalPart'), a)
    B['FractionalPart'] = _fracpart

    # ======================================================================== 数値関数 (非数値は未評価)
    def _flatnums(args):
        out = []
        for x in args:
            if isinstance(x, list):
                out.extend(_flatnums(x))
            else:
                out.append(x)
        return out

    def _allnum(xs):
        return all(_isint_num(x) for x in xs)

    def _minmax(name, fn, empty):
        def g(a):
            xs = _flatnums(a)
            if not _allnum(xs):
                return Expr(S(name), a)
            return fn(xs) if xs else empty
        return g
    B['Max'] = _minmax('Max', max, -INF)
    B['Min'] = _minmax('Min', min, INF)

    def _unary_real(name, fn, nargs=(1,)):
        def g(a):
            if len(a) not in nargs:
                raise SADError(f'General::narg: Number of arguments is expected to be {nargs[0]} in {name}')
            x = a[0]
            if isinstance(x, list):
                if not _allnum(_flatnums([x])):
                    return Expr(S(name), a)
                return _mapnum(x, lambda y: fn(y, a[1:]))
            if _isint_num(x):
                return fn(x, a[1:])
            return Expr(S(name), a)
        return g

    def _mapnum(x, f):
        return [_mapnum(y, f) for y in x] if isinstance(x, list) else f(x)

    def _round(x, rest):
        if rest:
            d = rest[0]
            return _round(x / d, []) * d
        return float(math.floor(x + 0.5)) if x >= 0 else -float(math.floor(-x + 0.5))
    B['Round'] = _unary_real('Round', _round, (1, 2))
    B['Floor'] = _unary_real('Floor', lambda x, r: float(math.floor(x)))
    B['Ceiling'] = _unary_real('Ceiling', lambda x, r: float(math.ceil(x)))

    def _mod(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in Mod')
        x, y = a
        if isinstance(x, list) or isinstance(y, list):
            if not _allnum(_flatnums([x, y])):
                return Expr(S('Mod'), a)

            def m2(p, q):
                if isinstance(p, list) and isinstance(q, list):
                    return [m2(u, v) for u, v in zip(p, q)]
                if isinstance(p, list):
                    return [m2(u, q) for u in p]
                if isinstance(q, list):
                    return [m2(p, v) for v in q]
                return _modn(p, q)
            return m2(x, y)
        if _isint_num(x) and _isint_num(y):
            return _modn(x, y)
        return Expr(S('Mod'), a)

    def _modn(x, y):
        if y == 0:
            raise SADError('General::wrongval: #2 must be nonzero in Mod')
        return x - y * math.floor(x / y)
    B['Mod'] = _mod

    # ======================================================================== 行列
    def _mat(m):
        return np.array(m, dtype=complex if _has_complex(m) else float)

    def _has_complex(m):
        if isinstance(m, list):
            return any(_has_complex(x) for x in m)
        return isinstance(m, complex)

    def _is_numeric_matrix(m):
        return (isinstance(m, list) and m and all(isinstance(r, list) and len(r) == len(m[0]) for r in m)
                and all(is_num(v) for r in m for v in r))

    def _tolist(m):
        if np.iscomplexobj(m):
            m = m.real if np.all(np.abs(m.imag) == 0) else m
        return m.tolist() if m.ndim else float(m)

    def _inverse(a):
        m = a[0]
        if not _is_numeric_matrix(m):
            raise SADError('General::wrongtype: Argument must be Numerical Matrix in Inverse')
        M = _mat(m)
        if not np.any(M):
            return [[NAN] * M.shape[1] for _ in range(M.shape[0])]
        if M.shape[0] == M.shape[1]:
            try:
                if np.linalg.matrix_rank(M) == M.shape[0]:
                    return _tolist(np.linalg.inv(M))
            except np.linalg.LinAlgError:
                pass
        return _tolist(np.linalg.pinv(M).T.T)
    B['Inverse'] = _inverse

    def _det(a):
        m = a[0]
        if not _is_numeric_matrix(m) or len(m) != len(m[0]):
            raise SADError('General::wrongtype: Argument must be Numerical Square Matrix in Det')
        d = np.linalg.det(_mat(m))
        return complex(d) if abs(np.imag(d)) > 0 else float(np.real(d))
    B['Det'] = _det

    def _dot(a):
        if len(a) > 2:
            return []
        x, y = a
        if not (isinstance(x, list) and isinstance(y, list)):
            return Expr(S('Dot'), a)
        try:
            X, Y = _mat(x), _mat(y)
        except (ValueError, TypeError):
            raise SADError('General::wrongtype: Argument must be numerical in Dot')
        if X.shape[-1] != Y.shape[0]:
            raise SADError('General::equalleng: Equal lengths of arguments are expected in Dot')
        return _tolist(np.dot(X, Y))
    B['Dot'] = _dot

    def _identity(a):
        n = _int(a[0])
        if n <= 0:
            raise SADError('General::wrongnum: Number must be Positive in IdentityMatrix')
        return np.eye(n).tolist()
    B['IdentityMatrix'] = _identity

    def _tr(a):
        m = a[0]
        if not (isinstance(m, list) and m and all(isinstance(r, list) for r in m)):
            raise SADError('General::wrongtype: Argument must be Matrix in Tr')
        n = min(len(m), min(len(r) for r in m))
        d = [m[k][k] for k in range(n)]
        if len(a) > 1 and isinstance(a[1], Sym) and a[1].name == 'List':
            return d
        return I.builtins['Plus'](d) if d else 0.0
    B['Tr'] = _tr

    def _outer(a):
        f, ls = a[0], a[1:]
        if not ls:
            return Expr(S('Outer'), a)
        for l in ls:
            need_list(l)

        def rec(k, acc):
            if k == len(ls):
                return apply(f, acc)
            return [rec(k + 1, acc + [x]) for x in ls[k]]
        return rec(0, [])
    B['Outer'] = _outer

    def _thread(a):
        x = a[0]
        if is_num(x) or isinstance(x, str):
            raise SADError('General::wrongtype: Argument must be List or composition for #1 in Thread')
        if isinstance(x, list):
            if x and all(isinstance(y, list) for y in x) and len({len(y) for y in x}) == 1:
                return [list(r) for r in zip(*x)]
            return x
        if isinstance(x, Rule):
            if isinstance(x.lhs, list) and isinstance(x.rhs, list):
                if len(x.lhs) != len(x.rhs):
                    raise SADError('General::equalleng: Equal lengths of elements in Thread are expected')
                return [Rule(l, r, x.delayed) for l, r in zip(x.lhs, x.rhs)]
            if isinstance(x.lhs, list):
                return [Rule(l, x.rhs, x.delayed) for l in x.lhs]
            return x
        if isinstance(x, Expr):
            lists = [y for y in x.args if isinstance(y, list)]
            if not lists:
                return x
            if len({len(l) for l in lists}) != 1:
                raise SADError('General::equalleng: Equal lengths of elements in Thread are expected')
            n = len(lists[0])
            return [I.eval(Expr(x.head, [_Lit(y[k] if isinstance(y, list) else y) for y in x.args])) for k in range(n)]
        return x
    B['Thread'] = _thread

    # ======================================================================== 文字列
    def _string_replace(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in StringReplace')
        s, rules = a
        if not isinstance(s, str):
            raise SADError('General::wrongtype: Argument must be Character-string in StringReplace')
        if isinstance(rules, Rule):
            rules = [rules]
        if not isinstance(rules, list) or not all(isinstance(r, Rule) for r in rules):
            raise SADError('General::wrongtype: Argument must be Rule or List of Rules in StringReplace')
        if not rules:
            return s
        table = {}
        for r in rules:
            table.setdefault(r.lhs, r.rhs)
        pat = re.compile('|'.join(re.escape(r.lhs) for r in rules))
        return pat.sub(lambda m: I.to_str(table[m.group(0)]), s)
    B['StringReplace'] = _string_replace

    def _fromchar(a):
        x = a[0]
        if isinstance(x, list):
            return ''.join(chr(int(c)) for c in x)
        return chr(int(x))
    B['FromCharacterCode'] = _fromchar
    B['ToCharacterCode'] = lambda a: [float(ord(c)) for c in a[0]]

    def _toexpression(a):
        if not isinstance(a[0], str):
            raise SADError('General::wrongtype: Argument must be Character-string in ToExpression')
        t = a[0].rstrip()
        if t and t[-1] in '+-*/^=,([{&|<>~':
            raise SADError('General::incomplete: Incomplete expression ' + a[0])
        return I.run(a[0], 'script')
    B['ToExpression'] = _toexpression

    B['Evaluate'] = lambda a: a[0] if len(a) == 1 else Seq(list(a))
    B['ReleaseHold'] = lambda a: a[0]

    # ======================================================================== 数学関数: 未定義シンボルには未評価のまま返す
    def _wrap_math(name):
        old = B[name]

        def g(a, old=old, name=name):
            if any(isinstance(x, (Sym, Expr)) for x in _flatnums(list(a))):
                return Expr(S(name), list(a))
            return old(a)
        B[name] = g
    for _nm in ('Sqrt', 'Exp', 'Log', 'Sin', 'Cos', 'Tan', 'ArcSin', 'ArcCos', 'ArcTan', 'Sinh', 'Cosh', 'Tanh',
                'Abs', 'Sign', 'Re', 'Im', 'Conjugate', 'Arg', 'Factorial', 'Gamma', 'Erf', 'Erfc'):
        if _nm in B:
            _wrap_math(_nm)
