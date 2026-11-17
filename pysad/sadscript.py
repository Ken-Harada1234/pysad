"""
sadscript.py — SADScript (SAD の Mathematica 風言語) のインタプリタ

SAD の FFS 区間 (FFS; 以降) のスクリプトを実行するための字句解析・構文解析・評価器。
演算子と優先順位は SAD の tfstk.f (opcode, iprior) に合わせる:
  ^ (15, 右結合)  . (10)  * / (40)  + - (50)  // (80, StringJoin)  == <> > >= < <= (100)
  === <=> (120)  ~ (150, Not)  && (160)  || (170)  : (180)  -> :> (190)  /. //. (200)
  & (220, Function)  = := ^= ^:= (250)  =. (900)  ; (1000)  , (9000)
  /@ //@ @@ (9)  @ (5)  [ [[ (6)  ++ -- (7)  += -= *= /= (210)
値: 実数は float (True=1, False=0 も実数)、文字列 str、リスト list、未定義シンボル Sym、
Null、規則 Rule、純関数 Func、未評価式 Expr。
Print/Write/ToString の実数書式は SAD の autofg.f ($FORM) を直訳した fmt_real()。

FFS コマンド (cell calc, fit ..., free ...) は文の先頭語で判別し、ffs_handler に文字列で渡す
(SAD の tffsa.f と同じく「SADScript として解釈できない語はコマンド」)。
"""
from __future__ import annotations

import itertools
import math
import cmath
import os
import re
import functools
import sys
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from . import compat as _compat

# ============================================================================ 値の型
INF = float('inf')
NAN = float('nan')


class Sym:
    __slots__ = ('name',)

    def __init__(self, name):
        self.name = name

    def __eq__(self, o):
        return isinstance(o, Sym) and o.name == self.name

    def __hash__(self):
        return hash(('Sym', self.name))

    def __repr__(self):
        return self.name


class _VarRef(Sym):
    """f[Unevaluated$[n]] で渡された変数 n そのもの (関数の中で n の値として評価され、i++ や代入は n に届く)。
    パターン照合で記号に束縛された値 (f[x_*y_] の y が記号 x のとき) とは区別する"""
    __slots__ = ()


class _Null:
    __slots__ = ()

    def __repr__(self):
        return 'Null'

    def __bool__(self):
        return False


Null = _Null()


class Expr:
    """未評価の式 head[args]"""
    __slots__ = ('head', 'args')

    def __init__(self, head, args):
        self.head = head
        self.args = list(args)

    def __repr__(self):
        return f'{self.head!r}[{", ".join(map(repr, self.args))}]'

    def __eq__(self, o):
        return isinstance(o, Expr) and o.head == self.head and o.args == self.args

    def __hash__(self):
        return hash(('Expr', repr(self)))


class Rule:
    __slots__ = ('lhs', 'rhs', 'delayed')

    def __init__(self, lhs, rhs, delayed=False):
        self.lhs = lhs
        self.rhs = rhs
        self.delayed = delayed

    def __repr__(self):
        return f'{self.lhs!r}{":>" if self.delayed else "->"}{self.rhs!r}'

    def __eq__(self, o):
        return isinstance(o, Rule) and o.lhs == self.lhs and o.rhs == self.rhs

    def __hash__(self):
        return hash(('Rule', repr(self)))


class Func:
    """純関数 (body&) または Function[{x..}, body]"""
    __slots__ = ('params', 'body', 'env')

    def __init__(self, params, body, env=None):
        self.params = params
        self.body = body
        self.env = env

    def __repr__(self):
        return f'Function[{self.params!r},{self.body!r}]'


class Seq:
    """x__ / x___ に束縛された引数列"""
    __slots__ = ('items',)

    def __init__(self, items):
        self.items = list(items)


class _Return(Exception):
    def __init__(self, value):
        self.value = value


class _Break(Exception):
    pass


class _Continue(Exception):
    pass


class _Throw(Exception):
    def __init__(self, value, tag=None):
        self.value = value
        self.tag = tag


class SADError(Exception):
    pass


def S(name):
    return Sym(name)


# ============================================================================ 字句解析
_OPS3 = ['//.', '^:=', '//@', '...', '===', '<=>']
_OPS2 = [':=', '=.', '==', '<>', '<=', '>=', '&&', '||', '//', '->', ':>', '/.', '^=', '/@', '@@',
         '..', '+=', '-=', '*=', '/=', '++', '--', '::', '/:', '##', '[[', '=>', '=<', '><']
_OPS1 = list('+-*/^=<>~[]{}(),;&:?#.|@')
_SYMCH = re.compile(r'[A-Za-z$`][A-Za-z0-9$`]*')
_NUM = re.compile(r'(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?')


class Tok:
    __slots__ = ('kind', 'val', 'line', 'pre')

    def __init__(self, kind, val, line, pre=False):
        self.kind = kind      # 'num', 'str', 'sym', 'op', 'nl', 'eof'
        self.val = val
        self.line = line
        self.pre = pre        # 直前に空白があるか

    def __repr__(self):
        return f'{self.kind}:{self.val!r}'


def strip_comments(src: str) -> str:
    """! から行末、(* ... *) を除く (文字列の中は除かない)。改行は残す"""
    out = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c == '"' or c == "'":
            q = c
            j = i + 1
            while j < n and src[j] != q:
                if src[j] == '\\':
                    j += 1
                j += 1
            out.append(src[i:j + 1])
            i = j + 1
        elif c == '!':
            j = src.find('\n', i)
            i = n if j < 0 else j
        elif c == '(' and src.startswith('(*', i):
            j = src.find('*)', i + 2)
            seg = src[i:j + 2] if j >= 0 else src[i:]
            out.append('\n' * seg.count('\n'))
            i = n if j < 0 else j + 2
        else:
            out.append(c)
            i += 1
    return ''.join(out)


_UNITS = {'RAD': 1.0, 'MRAD': 1e-3, 'DEG': 3.141592653589793 / 180.0, 'M': 1.0, 'CM': 1e-2, 'MM': 1e-3, 'T': 1.0, 'GAUSS': 1e-4,
          'EV': 1.0, 'KEV': 1e3, 'MEV': 1e6, 'GEV': 1e9, 'TEV': 1e12, 'V': 1.0, 'KV': 1e3, 'MV': 1e6, 'GV': 1e9,
          'HZ': 1.0, 'KHZ': 1e3, 'MHZ': 1e6, 'GHZ': 1e9, 'SEC': 1.0, 'MSEC': 1e-3, 'USEC': 1e-6, 'NSEC': 1e-9, 'PSEC': 1e-12}
_UNIT_SUFFIX = re.compile('(?:' + '|'.join(sorted(_UNITS, key=len, reverse=True)) + r')(?![A-Za-z0-9_$])', re.I)


_HEXD = '0123456789abcdefABCDEF'


_CMPSYM = {'Equal': '==', 'Unequal': '<>', 'Greater': '>', 'GreaterEqual': '>=', 'Less': '<', 'LessEqual': '<='}


def _is_const(a):
    """記号を含まない (数・文字列・それらのリスト/式/Blank) 引数"""
    if isinstance(a, (float, int, str, complex)):
        return True
    if isinstance(a, list):
        return all(_is_const(x) for x in a)
    if isinstance(a, Expr):
        if a.head == S('Blank'):
            return True
        return (isinstance(a.head, Sym) and a.head.name in ('List', 'Plus', 'Times', 'Power', 'Blank', 'Alternatives', 'Minus', 'Divide', 'Subtract')
                and all(_is_const(x) for x in a.args))
    return isinstance(a, Sym) and a.name in ('Pi', 'E', 'Degree')


def _unescape_at(src, j, n):
    r"""文字列の中の `\` の次の位置 j から 1 つの文字 (列) を読む。戻りは (文字列, 次の位置)。tfetok.f kxmakestring (gs で確認):
    `\` + 8 進 1〜3 桁 (値は 256 で折り返す)、`\x` / `\X` + 16 進 0〜2 桁 (桁がなければ NUL)、`\u` / `\U` + 16 進 4 / 8 桁まで、
    `\n` `\r` `\t` `\f` `\e` (ESC)、`\` + 改行は行継続 (何も入れない)。それ以外は `\` を落として次の文字をそのまま入れる (`\q` は q)"""
    c = src[j]
    if '0' <= c <= '7':
        k = j + 1
        v = int(c)
        while k < n and k < j + 3 and '0' <= src[k] <= '7':
            v = v * 8 + int(src[k])
            k += 1
        return chr(v % 256), k
    if c in 'xX' and j + 1 < n:
        k = j + 1
        v = 0
        while k < n and k < j + 3 and src[k] in _HEXD:
            v = v * 16 + int(src[k], 16)
            k += 1
        return chr(v), k
    if c in 'uU' and j + 1 < n:
        k = j + 1
        v = 0
        while k < n and k < j + (9 if c == 'U' else 5) and src[k] in _HEXD:
            v = v * 16 + int(src[k], 16)
            k += 1
        if k > j + 1:
            try:
                return chr(v), k
            except (ValueError, OverflowError):
                pass
        return c, j + 1
    if c == '\n':
        return '', j + 1
    return {'n': '\n', 'r': '\r', 't': '\t', 'f': '\f', 'e': '\x1b'}.get(c, c), j + 1


def tokenize(src: str) -> List[Tok]:
    src = strip_comments(src)
    toks = []
    i, n, line = 0, len(src), 1
    pre = False
    while i < n:
        c = src[i]
        if c == '\n':
            toks.append(Tok('nl', '\n', line))
            line += 1
            i += 1
            pre = True
            continue
        if c in ' \t\r':
            i += 1
            pre = True
            continue
        if c == '\\' and i + 1 < n and src[i + 1] == '\n':
            i += 2
            line += 1
            continue
        if c == '"' or c == "'":
            # SAD (tfetok.f) では ' も文字列の区切り (同じ引用符で閉じる)
            q = c
            j = i + 1
            buf = []
            while j < n and src[j] != q:
                if src[j] == '\\' and j + 1 < n:
                    if src[j + 1] == '\n':
                        line += 1
                    ch_, j = _unescape_at(src, j + 1, n)
                    buf.append(ch_)
                else:
                    if src[j] == '\n':
                        line += 1
                    buf.append(src[j])
                    j += 1
            toks.append(Tok('str', ''.join(buf), line, pre))
            i = j + 1
            pre = False
            continue
        if c == '%':
            # % は直前の出力 Out[$Line]、%n は Out[n] (gs)。%% は gs でも記号のまま
            j = i + 1
            while j < n and src[j].isdigit():
                j += 1
            if j < n and src[j] == '%':
                toks.append(Tok('sym', '%%', line, pre))
                i = j + 1
            else:
                toks.append(Tok('sym', 'Out', line, pre))
                toks.append(Tok('op', '[', line, False))
                toks.append(Tok('num', float(src[i + 1:j]), line, False) if j > i + 1 else Tok('sym', '$Line', line, False))
                toks.append(Tok('op', ']', line, False))
                i = j
            pre = False
            continue
        if c == '0' and i + 2 < n and src[i + 1] in 'xX' and src[i + 2] in '0123456789abcdefABCDEF':
            j = i + 2
            while j < n and src[j] in '0123456789abcdefABCDEF':
                j += 1
            toks.append(Tok('num', float(int(src[i + 2:j], 16)), line, pre))       # 0x33 などの 16 進整数
            i = j
            pre = False
            continue
        if c.isdigit() or (c == '.' and i + 1 < n and src[i + 1].isdigit()):
            m = _NUM.match(src, i)
            txt = m.group(0).replace('d', 'e').replace('D', 'e')
            val = float(txt)
            i = m.end()
            # SAD: 数の直後に単位 (1GEV, 2.5MV, 30DEG, 10MM …、空白なし) が付けば掛ける (lattice と同じ単位表)
            mu = _UNIT_SUFFIX.match(src, i)
            if mu:
                val *= _UNITS[mu.group(0).upper()]
                i = mu.end()
            toks.append(Tok('num', val, line, pre))
            pre = False
            continue
        m = _SYMCH.match(src, i)
        if m:
            nm_ = m.group(0)
            # `X・System`X・Global`X は、Begin[…] で文脈を作らない限り大域の X (OpticsPlot.n は `OpticsPlot=System`OpticsPlot=…)
            if '`' in nm_ and nm_.count('`') == 1 and (nm_[0] == '`' or nm_.startswith(('System`', 'Global`'))) and len(nm_) > nm_.index('`') + 1                     and not (_CUR is not None and (_CUR._ctx_cur != '`' or _CUR._ctx_path != _ROOT_PATH)):
                nm_ = nm_[nm_.index('`') + 1:]          # 文脈の中 (BeginPackage / Begin) では _ctx_resolve が大域の記号として解決する
            toks.append(Tok('sym', nm_, line, pre))
            i = m.end()
            pre = False
            continue
        # 演算子: 長いものから。ただし =. /. //. ++ -- の直後が数字なら短い方 (SAD tfetok)
        op = None
        for cand in _OPS3 + _OPS2:
            if src.startswith(cand, i):
                op = cand
                break
        if op in ('=.', '/.', '//.', '++', '--') and i + len(op) < n and src[i + len(op)].isdigit():
            op = {'=.': '=', '/.': '/', '//.': '//', '++': '+', '--': '-'}[op]
        if op is None:
            if c in _OPS1:
                op = c
            else:
                raise SADError(f'line {line}: unexpected character {c!r}')
        toks.append(Tok('op', op, line, pre))
        i += len(op)
        pre = False
    toks.append(Tok('eof', None, line))
    return toks


# ============================================================================ 構文解析 (Pratt)
# 二項演算子: op → (優先度 iprior, 右結合, head)
_BIN = {
    '^': (15, True, 'Power'), '.': (10, False, 'Dot'),
    '*': (40, False, 'Times'), '/': (40, False, 'Divide'),
    '+': (50, False, 'Plus'), '-': (50, False, 'Subtract'),
    '//': (80, False, 'StringJoin'),
    '==': (100, False, 'Equal'), '<>': (100, False, 'Unequal'), '><': (100, False, 'Unequal'),
    '>': (100, False, 'Greater'), '>=': (100, False, 'GreaterEqual'), '=>': (100, False, 'GreaterEqual'),
    '<': (100, False, 'Less'), '<=': (100, False, 'LessEqual'), '=<': (100, False, 'LessEqual'),
    '===': (120, False, 'SameQ'), '<=>': (120, False, 'UnsameQ'),
    '&&': (160, False, 'And'), '||': (170, False, 'Or'),
    ':': (180, True, 'Pattern'), '?': (4, False, 'PatternTest'),
    '->': (190, True, 'Rule'), ':>': (190, True, 'RuleDelayed'),
    '/.': (200, False, 'ReplaceAll'), '//.': (200, False, 'ReplaceRepeated'),
    '=': (250, True, 'Set'), ':=': (250, True, 'SetDelayed'), '^=': (250, True, 'UpSet'), '^:=': (250, True, 'UpSetDelayed'),
    '+=': (210, True, 'AddTo'), '-=': (210, True, 'SubtractFrom'), '*=': (210, True, 'TimesBy'), '/=': (210, True, 'DivideBy'),
    '/@': (9, True, 'Map'), '//@': (9, True, 'MapAll'), '@@': (9, True, 'Apply'), '@': (5, True, 'At'),
    '|': (175, False, 'Alternatives'), '/:': (240, True, 'TagSet'),
}


class Parser:
    def __init__(self, toks: List[Tok]):
        self.toks = toks
        self.i = 0
        self.depth = 0

    def peek(self, k=0):
        return self.toks[min(self.i + k, len(self.toks) - 1)]

    def next(self):
        t = self.toks[self.i]
        self.i += 1
        return t

    def skip_nl(self):
        while self.peek().kind == 'nl':
            self.i += 1

    def at_op(self, *vals):
        t = self.peek()
        return t.kind == 'op' and t.val in vals

    def expect(self, val):
        t = self.next()
        if not (t.kind == 'op' and t.val == val):
            raise SADError(f'line {t.line}: expected {val!r}, got {t.val!r}')
        return t

    # ---- 文の列 (トップレベル / ファイル)
    def parse_statements(self) -> List[Any]:
        """トップレベル: 改行または ; で区切った文のリスト。各要素は AST か ('ffs', text)"""
        stmts = []
        while True:
            self.skip_nl()
            if self.peek().kind == 'eof':
                break
            if self.at_op(';'):
                self.next()
                continue
            e = self.parse_expr(9000)
            stmts.append(e)
            t = self.peek()
            if t.kind == 'op' and t.val == ';':
                self.next()
            elif t.kind in ('nl', 'eof'):
                pass
            else:
                raise SADError(f'line {t.line}: unexpected token {t.val!r}')
        return stmts

    # ---- 式
    def parse_expr(self, rbp=9999):
        """優先度 rbp より強く結合する式を読む (値が小さいほど強い)"""
        if self.depth > 0:
            self.skip_nl()
        left = self.parse_prefix()
        while True:
            t = self.peek()
            nl_skipped = False
            if t.kind == 'nl':
                if self.depth > 0:
                    # 括弧内では改行は空白。ただし次が演算子でなければ、式の区切り (, ; ]) のはず
                    self.skip_nl()
                    t = self.peek()
                    nl_skipped = True
                else:
                    break
            if t.kind == 'eof':
                break
            if t.kind == 'op':
                v = t.val
                if v == '[':
                    if 6 >= rbp:
                        break
                    left = self.parse_call(left)
                    continue
                if v == '[[':
                    if 6 >= rbp:
                        break
                    self.next()
                    self.depth += 1
                    args = self.parse_args(']')
                    self.expect(']')
                    self.depth -= 1
                    left = Expr(S('Part'), [left] + args)
                    continue
                if v == '::':
                    # Sym::tag (メッセージ名)。Sym::tag = "..." でメッセージを定義、Message[Sym::tag, ...] / Off[Sym::tag] で使う
                    if 6 >= rbp:
                        break
                    self.next()
                    tg = self.next()
                    if tg.kind not in ('sym', 'str'):
                        raise SADError(f'line {tg.line}: expected message tag after ::')
                    left = Expr(S('MessageName'), [left, tg.val if tg.kind == 'str' else tg.val])
                    continue
                if v == '&':
                    if 220 >= rbp:
                        break
                    self.next()
                    left = Expr(S('Function'), [left])
                    continue
                if v in ('++', '--'):
                    if 7 >= rbp:
                        break
                    self.next()
                    left = Expr(S('Increment' if v == '++' else 'Decrement'), [left])
                    continue
                if v == '=.':
                    if 900 >= rbp:
                        break
                    self.next()
                    left = Expr(S('Unset'), [left])
                    continue
                if v == ';':
                    if 1000 >= rbp:
                        break
                    self.next()
                    if self.depth > 0:
                        self.skip_nl()
                    nt = self.peek()
                    if nt.kind in ('eof', 'nl') or (nt.kind == 'op' and nt.val in (']', ')', '}', ',', ';')):
                        right = Null
                    else:
                        right = self.parse_expr(1000)
                    if isinstance(left, Expr) and left.head == S('CompoundExpression'):
                        left = Expr(S('CompoundExpression'), left.args + [right])
                    else:
                        left = Expr(S('CompoundExpression'), [left, right])
                    continue
                if v == '..' or v == '...':
                    self.next()
                    left = Expr(S('Repeated' if v == '..' else 'RepeatedNull'), [left])
                    continue
                if v == '@':
                    nt_ = self.peek(1)
                    if nt_ is not None and nt_.kind == 'sym' and 6 < rbp:
                        # x@member (クラスのメンバー参照): 右は記号だけ。x@m[args] は (x@m)[args]、x@a@b は (x@a)@b
                        self.next()
                        nm_ = self.next()
                        left = Expr(S('At'), [left, Sym(nm_.val)])
                        continue
                if v in _BIN:
                    prio, rassoc, head = _BIN[v]
                    if prio > rbp or (prio == rbp and not rassoc):
                        break
                    self.next()
                    if self.depth == 0:
                        self.skip_nl()
                    if v == '@@' and self.at_op('['):
                        # f@@[x, {1}] = Apply[f, x, {1}] (レベル指定付きの Apply。Tkinter.n, KBFrame.n, Canvas.n)
                        self.next()
                        self.depth += 1
                        extra = self.parse_args(']')
                        self.depth -= 1
                        left = Expr(S('Apply'), [left] + extra)
                        continue
                    if v in (':>', '->') and self.peek().kind == 'op' and self.peek().val in (',', ']', ')', '}'):
                        right = Null                       # `Command$:>,` のように右辺が空のルールは Null
                    elif v == '?':
                        right = self.parse_prefix()
                    else:
                        right = self.parse_expr(prio)      # 左結合は同じ優先度で break する
                    left = self.make_bin(head, left, right)
                    continue
                if v not in ('(', '#', '##'):
                    break
            # 数値・シンボル・( が続く: 暗黙の積 (SAD は 2 Pi, a b, 2 (x+1) を Times とする。改行をまたぐ場合は積にしない)
            prev_ = self.toks[self.i - 1] if self.i > 0 else None
            closed_ = prev_ is not None and ((prev_.kind == 'op' and prev_.val in (')', ']', '}')) or (prev_.kind == 'num' and t.kind in ('sym', 'op')))      # f[x]g[y]、(a)(b)、f[#]#2 は積 (D.n)
            if (not nl_skipped and (t.kind in ('num', 'sym') or (t.kind == 'op' and t.val in ('(', '#', '##')))
                    and (t.pre or (t.kind == 'sym' and isinstance(left, float)) or closed_
                         or (t.kind == 'op' and t.val == '(' and isinstance(left, Sym))) and 40 < rbp):      # s2(a|b): 記号の直後の ( も積 (GeometryPlot.n の s2(_Real|_String))
                right = self.parse_expr(40)
                left = Expr(S('Times'), [left, right])
                continue
            break
        return left

    def make_bin(self, head, left, right):
        if head == 'At':
            return Expr(left, [right])
        if head == 'Pattern':
            # x_ : default  → Optional
            if isinstance(left, Expr) and left.head == S('Pattern'):
                return Expr(S('Optional'), [left, right])
            if isinstance(left, Expr) and left.head == S('Blank'):
                return Expr(S('Optional'), [left, right])
            return Expr(S('Pattern'), [left, right])
        if (head in ('SameQ', 'UnsameQ', 'Less', 'LessEqual', 'Greater', 'GreaterEqual', 'Equal') and isinstance(left, Expr)
                and left.head == S(head) and not getattr(left, 'paren', False)):
            return Expr(S(head), left.args + [right])           # a===b===c, a<=>b<=>c は 1 つの n 項
        if head in ('Plus', 'Times'):
            la = left.args if isinstance(left, Expr) and left.head == S(head) else [left]      # Flat: a*b*c は Times[a,b,c] (パターンの b___*Sin[x_]*Cos[x_] を 1 つの Times にする)
            ra = right.args if isinstance(right, Expr) and right.head == S(head) else [right]
            return Expr(S(head), la + ra)
        return Expr(S(head), [left, right])

    def parse_args(self, close):
        args = []
        self.skip_nl()
        if self.at_op(close):
            self.next()
            return args
        while True:
            self.skip_nl()
            if self.at_op(','):
                args.append(Null)
            else:
                args.append(self.parse_expr(9000))
            self.skip_nl()
            t = self.next()
            if t.kind == 'op' and t.val == ',':
                self.skip_nl()
                if self.at_op(close):
                    args.append(Null)
                    self.next()
                    return args
                continue
            if t.kind == 'op' and t.val == close:
                return args
            if close == ']' and t.kind == 'op' and t.val == '[[':
                raise SADError(f'line {t.line}: unexpected [[')
            raise SADError(f'line {t.line}: expected , or {close}, got {t.val!r}')

    def parse_call(self, head):
        self.next()   # [
        self.depth += 1
        args = self.parse_args(']')
        self.depth -= 1
        # f[[i]] が '[' '[' と来た場合
        return Expr(head, args)

    def parse_prefix(self):
        if self.depth > 0:
            self.skip_nl()
        t = self.next()
        if t.kind == 'num':
            return t.val
        if t.kind == 'str':
            return t.val
        if t.kind == 'sym':
            return self.parse_symbol(t.val)
        if t.kind == 'op':
            v = t.val
            if v == '-':
                e = self.parse_expr(20)   # -2^2 = -4, -a.b = -(a.b) (SAD で確認; Power=15, Dot=10 を含める)
                if isinstance(e, float):
                    return -e
                return Expr(S('Minus'), [e])
            if v == '+':
                return self.parse_expr(20)
            if v == '~':
                return Expr(S('Not'), [self.parse_expr(150)])
            if v == '?':
                # SAD: ?FLAG (FFS のフラグの問い合わせ、True/False)。FFS が FFS$FlagQ を登録する
                t2 = self.peek()
                if t2.kind != 'sym':
                    # `X->??` (Tk の bind の %x が イベントに無いときの ?? ) など、? の後に名前が無ければ Null
                    return Null
                self.next()
                return Expr(S('FFS$FlagQ'), [t2.val])
            if v == '++':
                return Expr(S('PreIncrement'), [self.parse_expr(7)])
            if v == '--':
                return Expr(S('PreDecrement'), [self.parse_expr(7)])
            if v == '(':
                self.depth += 1
                self.skip_nl()
                if self.at_op(')'):
                    self.next()
                    self.depth -= 1
                    return Null
                e = self.parse_compound(')')
                self.depth -= 1
                return e
            if v == '{':
                self.depth += 1
                args = self.parse_args('}')
                self.depth -= 1
                return Expr(S('List'), args)
            if v == '#':
                if self.peek().kind == 'num' and not self.peek().pre:
                    return Expr(S('Slot'), [self.next().val])
                return Expr(S('Slot'), [1.0])
            if v == '##':
                if self.peek().kind == 'num' and not self.peek().pre:
                    return Expr(S('SlotSequence'), [self.next().val])
                return Expr(S('SlotSequence'), [1.0])
            if v == '_':
                return self.parse_blank(None)
            if v == '[':
                # SAD: 前置の [x] は括弧 (Map[f,[g[x]]] の [..] は丸括弧と同じ。corutility.n の CorStat)
                self.depth += 1
                args = self.parse_args(']')
                self.depth -= 1
                return args[0] if len(args) == 1 else Expr(S('Sequence'), args)
            if v == '{':
                pass
        if t.kind == 'nl':
            return self.parse_prefix()
        if t.kind == 'eof':
            return Null
        raise SADError(f'line {t.line}: unexpected token {t.val!r}')

    def parse_compound(self, close):
        """( e1; e2; ... ) — ; 区切り。改行も区切りになる"""
        exprs = []
        while True:
            self.skip_nl()
            if self.at_op(close):
                self.next()
                break
            if self.at_op(';'):
                self.next()
                exprs.append(Null)
                continue
            exprs.append(self.parse_expr(9000))
            self.skip_nl()
            if self.at_op(';'):
                self.next()
                self.skip_nl()
                if self.at_op(close):
                    exprs.append(Null)
                    self.next()
                    break
                continue
            if self.at_op(close):
                self.next()
                break
            t = self.peek()
            raise SADError(f'line {t.line}: expected ; or {close}, got {t.val!r}')
        if len(exprs) == 1:
            return exprs[0]
        return Expr(S('CompoundExpression'), exprs)

    def parse_symbol(self, name):
        """シンボル。x_, x__, x___, x_Real, _Real などのパターンも"""
        # 名前の中の '_' はトークナイザで分けていない: SYMCH に '_' を含めないので x_ は sym 'x' + op? '_' は OPS にない
        return S(name)

    def parse_blank(self, name):
        return Expr(S('Blank'), [])


# '_' をパターンとして扱うため、トークナイザの後処理でシンボル名の '_' を処理する
def _retokenize_patterns(toks: List[Tok]) -> List[Tok]:
    """'x_', 'x__', 'x___', 'x_Real', '_', '_Real' を構文木用のトークン列に変換する。
    字句解析では '_' を記号として読めないので、ここでソース中の '_' を含む語を sym 'x' + op '_' + sym 'Real' に分ける"""
    return toks


_PAT = re.compile(r'^([A-Za-z$`][A-Za-z0-9$`]*)?(_{1,3})([A-Za-z$`][A-Za-z0-9$`]*)?$')


def tokenize_full(src: str) -> List[Tok]:
    """'_' をシンボル文字として読み、パターン語を分解する"""
    global _SYMCH
    saved = _SYMCH
    _SYMCH = re.compile(r'[A-Za-z$`_][A-Za-z0-9$`_]*')
    try:
        toks = tokenize(src)
    finally:
        _SYMCH = saved
    out = []
    for t in toks:
        if t.kind == 'sym' and '_' in t.val:
            m = _PAT.match(t.val)
            if not m:
                out.append(t)
                continue
            name, blanks, typ = m.groups()
            out.append(Tok('pat', (name, len(blanks), typ), t.line, t.pre))
        else:
            out.append(t)
    if _CUR is not None:
        out = _CUR._ctx_tokens(out)          # 文脈 (BeginPackage / Begin) の名前解決。利用者のコードを読むときだけ
    return out


# ---- 文脈 (名前空間)。_CUR はいま利用者のコードを読んでいる Interp (SAD の Packages を読む間と初期化中は None)
_CUR = None
_ROOT_PATH = ['`', 'System`']
_SYS_FIXED = ('Real', 'String', 'List', 'Symbol', 'Complex', 'Null', 'True', 'False', 'Infinity', 'INF', 'NaN', 'Pi', 'E', 'I', 'Function',
              'Rule', 'RuleDelayed', 'Constant', 'HoldAll', 'HoldFirst', 'HoldRest', 'HoldNone', 'Undefined', 'Automatic', 'All', 'None',
              'Vector', 'Pattern', 'Blank', 'Integer', 'Number', 'Expression', 'Character', 'Word', 'Byte', 'EndOfFile')
_SYSNAMES = None


def _sysnames():
    """pysad が記号として扱う名前 (pysad のソース中の S('Name')・.name == 'Name' の文字列)。文脈の名前解決で System` の記号とみなす。
    大文字か $ で始まる 2 文字以上の名前だけ (1 文字は I・E・D・N)。最初に文脈を使ったときに 1 回だけ集める"""
    global _SYSNAMES
    if _SYSNAMES is None:
        names = set(_SYS_FIXED)
        rx1 = re.compile(r"\bS\(\s*'([A-Z$][A-Za-z0-9$]*)'\s*\)|\.name\s*(?:==|!=)\s*'([A-Z$][A-Za-z0-9$]*)'")
        rx2 = re.compile(r"\.name\s+(?:not\s+)?in\s+\(([^)]*)\)")
        rx3 = re.compile(r"'([A-Z$][A-Za-z0-9$]*)'")
        here = os.path.dirname(os.path.abspath(__file__))
        try:
            files = [f for f in os.listdir(here) if f.endswith('.py')]
        except OSError:
            files = []
        for f in files:
            try:
                txt = open(os.path.join(here, f), encoding='utf-8', errors='replace').read()
            except OSError:
                continue
            for m in rx1.finditer(txt):
                names.add(m.group(1) or m.group(2))
            for m in rx2.finditer(txt):
                names.update(rx3.findall(m.group(1)))
        _SYSNAMES = {n for n in names if len(n) > 1 or n in ('I', 'E', 'D', 'N')}
    return _SYSNAMES


def _ctx_norm(c):
    """文脈名の正規化: Stat`・`Stat`・"Stat`" → 'Stat`'、`・Global` → '`'"""
    if isinstance(c, Sym):
        c = c.name
    if not isinstance(c, str) or not c:
        raise SADError('General::wrongtype: Argument must be a context name (Name`)')
    if not c.endswith('`'):
        c += '`'
    if len(c) > 1 and c[0] == '`':
        c = c[1:]
    return '`' if c in ('`', 'Global`') else c


def _parse_prefix_pat(self: Parser, t: Tok):
    name, nb, typ = t.val
    head = {1: 'Blank', 2: 'BlankSequence', 3: 'BlankNullSequence'}[nb]
    blank = Expr(S(head), [S(typ)] if typ else [])
    if name:
        return Expr(S('Pattern'), [S(name), blank])
    return blank


_orig_parse_prefix = Parser.parse_prefix


def _parse_prefix2(self):
    if self.depth > 0:
        self.skip_nl()
    if self.peek().kind == 'pat':
        return _parse_prefix_pat(self, self.next())
    return _orig_parse_prefix(self)


Parser.parse_prefix = _parse_prefix2


def parse(src: str) -> List[Any]:
    return Parser(tokenize_full(src)).parse_statements()


# ============================================================================ 実数の書式 (autofg.f の直訳)
_ZERO = '0' * 31


def _strfromd(x: float):
    """autofg.f strfromd: |x| を 18 桁の数字列と指数 (iexp: 10^iexp ≥ |x| の桁数) に"""
    if x == 0.0:
        return _ZERO[:18], 0
    a = abs(x)
    if a < 2.0**31:
        ai = float(int(x))
        if ai == x:
            s = f'{abs(int(x)):010d}'
            for i in range(10):
                if s[i] != '0':
                    iexp = 10 - i
                    return (s[i:i + iexp] + _ZERO)[:18], iexp
    xm = [1e-256, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e4, 1e8, 1e16, 1e32, 1e64, 1e128, 1e256]
    mexp = [-256, -128, -64, -32, -16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
    iae = 9
    xl, xh = 1e8, 1e9
    if a >= xh:
        for i in range(0, 9):          # i = -9..-1
            a1 = a * xm[i]
            if a1 >= xl:
                a = a1
                iae -= mexp[i]
        if a >= xh:
            a *= 0.1
            iae += 1
        if a >= xh:
            return 'INF', 0
    elif a < xl:
        for i in range(18, 9, -1):     # i = 9..1
            a1 = a * xm[i]
            if a1 < xh:
                a = a1
                iae -= mexp[i]
        if a < xl:
            a *= 10.0
            iae -= 1
    iexp = iae
    iam = int(a)
    iaf = int((a - iam) * 1e9 + 0.5)
    if iaf > 999999999:
        iaf = 0
        iam += 1
        if iam > 999999999:
            iam //= 10
            iexp += 1
    return f'{iam:09d}' + (_ZERO[:9] if iaf == 0 else f'{iaf:09d}'), iexp


def _roundnumstr(s: str, icol: int):
    """autofg.f roundnumstr: 文字列 s (s[0] は '0' の桁上がり用) を icol 桁目 (1 始まり) で四捨五入"""
    l = len(s)
    if icol > l:
        return s, False
    ch = s[icol - 1]
    s = list(s[:icol - 1] + _ZERO[:l - icol + 1])
    j = icol - 2
    inc = ch >= '5'
    s[0] = '0'
    while inc:
        s[j] = chr(ord(s[j]) + 1)
        inc = s[j] > '9'
        if inc:
            s[j] = '0'
            j -= 1
    shift = False
    if s[0] == '1':
        shift = True
        for i in range(icol - 1, 0, -1):
            s[i] = s[i - 1]
        s[0] = '0'
    return ''.join(s), shift


def fmt_real(x: float, form: str = '') -> str:
    """SAD の autofg(x, form)。form: '15.10' / 'S15.10' / '' (= 最短表現)"""
    if form == '' or form == 'S' or form.strip() == '':      # K64: $FORM=" " は既定と同じ (最短表現。fm1.sad で確認)
        return _autos(x)
    if form[0] in 'FG':
        try:
            w, d = form[1:].split('.')
            s = f'{x:{int(w)}.{int(d)}f}' if form[0] == 'F' else f'{x:{int(w)}.{int(d)}g}'
            return s if len(s) <= int(w) else '$' * int(w)
        except Exception:
            return repr(x)
    if form[0] == 'E':
        w, d = form[1:].split('.')
        return f'{x:{int(w)}.{int(d) - 1}E}'
    tzero = form[0] in 'SMC'
    f = form[1:] if tzero else form
    if f.strip() == '':
        tzero = True
        lc, lf = 18, 15
    else:
        idot = f.index('.')
        lc = int(f[:idot] or 0)
        lf = int(f[idot + 1:] or 0)
    lc = min(18, lc)
    lf = min(15, lf)
    if isinstance(x, complex):
        return _fmt_complex(x, lambda v: fmt_real(v, form))
    if x != x:
        return 'NaN'   # SAD は $FORM でも NaN を幅なしで出す (確認済み)
    if x == 0.0:
        if tzero:
            return '0'
        if lf == 0:
            return ' ' * (lc - 1) + '0'
        return ' ' * (lc - lf - 2) + ' .' + _ZERO[:lf]
    if x == INF:
        return 'INF' if tzero else ' ' * (lc - 3) + 'INF'
    if x == -INF:
        return '-INF' if tzero else ' ' * (lc - 4) + '-INF'
    ax = abs(x)
    digits, iexp = _strfromd(ax)
    buff = '0' + digits          # buff(1:1)='0', buff(2:19)=digits
    if tzero:
        lc1 = lc if x >= 0 else lc - 1
        lc2 = lc1
    else:
        if x >= 0:
            lc1 = lc - min(lf, 1)
            lc2 = lc - 1
        else:
            lc1 = lc - 1 - min(lf, 1)
            lc2 = lc - 2
    expstr = ''
    while True:
        if iexp >= 101:
            lm = lc2 - 5
        elif iexp >= 11:
            lm = lc2 - 4
        elif iexp >= 1:
            lm = lc2 - 3
        elif iexp >= -8:
            lm = lc2 - 4
        elif iexp >= -98:
            lm = lc2 - 5
        else:
            lm = lc2 - 6
        if lm < 0:
            buf1 = '***' if iexp >= 1 else '0.'
            is_ = len(buf1)
            break
        if iexp <= min(10, lc2) and iexp > -min(3, lf):
            buff, shift = _roundnumstr(buff, min(lc1, iexp + lf) + 2)
            if shift:
                iexp += 1
                continue
            lf1 = min(lf, lc1 - 1)
            li1 = lc1 - lf1 - 1
            is_ = lc1
            if iexp == lc1:
                buf1 = buff[1:lc1 + 1]
            elif iexp >= li1:
                if lf == 0:
                    buf1 = buff[1:iexp + 1]
                else:
                    buf1 = buff[1:iexp + 1] + '.' + buff[iexp + 1:lc1]
            elif iexp >= 1:
                if lf == 0:
                    buf1 = ' ' * (li1 - iexp) + buff[1:iexp + 1]
                else:
                    buf1 = ' ' * (li1 - iexp) + buff[1:iexp + 1] + '.' + buff[iexp + 1:min(lc1, iexp + lf1 + 1)]
            else:
                if iexp == 0:
                    buf1 = ' ' * li1 + '.' + buff[1:lf1 + 1]
                else:
                    buf1 = ' ' * li1 + '.' + _ZERO[:-iexp] + buff[1:lf1 + iexp + 1]
            break
        else:
            buff, shift = _roundnumstr(buff, max(3, lm + 2))
            if shift:
                iexp += 1
                continue
            expstr = 'E' + str(iexp - 1)
            if lm == 0:
                buf1 = buff[1] + expstr
                is_ = 1
            elif lm == 1:
                buf1 = buff[1] + '.' + expstr
                is_ = 2
            else:
                buf1 = buff[1] + '.' + buff[2:lm + 1] + expstr
                is_ = lm + 1
            break
    if tzero and '.' in buf1:
        # 末尾の 0 と '.' を落とす (autofg の tzero 処理)
        core = buf1[:len(buf1) - len(expstr)] if expstr else buf1
        core = core.rstrip('0').rstrip('.')
        buf1 = core + expstr
    if tzero:
        s = buf1 if x >= 0 else '-' + buf1.lstrip()
        return s.strip()
    if x >= 0:
        s = ' ' + buf1
    else:
        s = '  ' + buf1
        i = len(s) - len(s[2:].lstrip())
        s = s[:i - 1] + '-' + s[i:]
    return s


def fmt_g(x: float, w: int = 11, d: int = 4) -> str:
    """Fortran の 1PGw.d (SAD の tffsmatch の表示): 0.1 ≤ |x| < 10^d は F(w-4).(d-k) + 4 空白、他は 1PE 形式"""
    ax = abs(x)
    if x == 0.0:
        return ('0.' + '0' * (d - 1)).rjust(w - 4) + ' ' * 4
    if 0.1 <= ax < 10.0**d:
        k = int(math.floor(math.log10(ax))) + 1

        def _f(nd):
            # Fortran の F7.0 は "1080." (小数点を付ける)。Python の .0f は付けない
            return (f'{x:.0f}' + '.').rjust(w - 4) if nd == 0 else f'{x:{w - 4}.{nd}f}'
        s = _f(d - k)
        if len(s) > w - 4:
            s = _f(max(0, d - k - 1))
        return s + ' ' * 4
    m = f'{x:.{d}E}'
    mant, ex = m.split('E')
    return f'{mant}E{int(ex):+03d}'.rjust(w)


def _fmt_complex(x: complex, f) -> str:
    """SAD の複素数表示: I, (-I), (im I), (re+I), (re-I), (re+im I), (re-im I)。各部は f で書式化 ($FORM が効く)。
    虚部が負の一般の場合は符号は虚部の書式の中 (例 $FORM=15.10 で (   1.0000000000  -2.0000000000 I))"""
    re_, im = x.real, x.imag
    if re_ == 0.0:
        if im == 1.0:
            return 'I'
        if im == -1.0:
            return '(-I)'
        return '(' + f(im) + ' I)'
    if im == 1.0:
        return '(' + f(re_) + '+I)'
    if im == -1.0:
        return '(' + f(re_) + '-I)'
    if im >= 0 or im != im:
        return '(' + f(re_) + '+' + f(im) + ' I)'
    return '(' + f(re_) + f(im) + ' I)'


def _autos(x) -> str:
    """autos (gdtoa g_dfmt: 最短の往復可能表現)"""
    if isinstance(x, complex):
        return _fmt_complex(x, _autos)
    if x != x:
        return 'NaN'
    if x == INF:
        return 'Infinity'
    if x == -INF:
        return '-Infinity'
    # gdtoa g_dfmt: 最短の桁列 (digits, decpt) から、decpt <= -4 または decpt > 桁数+5 なら指数表記 (指数が 10 以上のときだけ + を付ける: 1e6, 1e+10)、
    # それ以外は位取り表記 (K64 で確認。1e6 → 1e6、123456789012345678. → 123456789012345680)
    if x == 0:
        return '0'
    r = repr(abs(float(x)))
    if 'e' in r:
        mant, ex = r.split('e')
        ex = int(ex)
    else:
        mant, ex = r, 0
    ip, _, fp = mant.partition('.')
    digits = (ip + fp)
    decpt = len(ip) + ex
    lead = len(digits) - len(digits.lstrip('0'))
    digits = digits.lstrip('0')
    decpt -= lead
    digits = digits.rstrip('0') or '0'
    sign = '-' if x < 0 else ''
    nd = len(digits)
    if decpt <= -4 or decpt > nd + 5:
        e10 = decpt - 1
        m = digits[0] + ('.' + digits[1:] if nd > 1 else '')
        es = ('-' if e10 < 0 else ('+' if e10 >= 10 else '')) + str(abs(e10))
        return sign + m + 'e' + es
    if decpt <= 0:
        return sign + '.' + '0' * (-decpt) + digits
    if decpt >= nd:
        return sign + digits + '0' * (decpt - nd)
    return sign + digits[:decpt] + '.' + digits[decpt:]


# SAD の Packages/Fit.n (PolynomialFit) をそのまま評価する。xa=0 のとき xl=1 になる SAD の癖もそのまま
# v2 (compat.v2()、B18): xa=0 のときは平行移動が無いので係数は sl そのもの (K64 は xl=1 のまま二項係数の和を取って係数が壊れる)
_PRELUDE = '''
PolynomialFit[d_,nord_]:=Module[{{x,b}=Thread[d],
  s=Table[1,{Length[d]}],a,xa,sl,xl,cm,n,m,sx},
  xa=(Plus@@x)/Length[x];
  scale=1/Max[Abs[x-xa]];
  If[Abs[xa]*scale<0.1,
    xa=0;
    scale=1/Max[Abs[x]]];
  x=(x-xa)*scale;
  a=Transpose[Table[First[{s,s*=x}],{nord+1}]];
  s=LinearSolve[a,b];
  sl=xl=1;
  sl=s*Table[First[{sl,sl*=scale}],{nord+1}];
  xl=If[xa,Table[First[{xl,xl*=-xa}],{nord+1}],
      1];
  sx=sl*xl;
  cm=If[xa==0&&Compat$V2[],sl,Table[
    Take[sx,{m+1,nord+1}].
    Table[Factorial[n]/Factorial[n-m]/Factorial[m],{n,m,nord}],
    {m,0,nord}]/xl];
  {cm,{Residual->Sqrt[(Plus@@((a.s-b)^2))/Length[d]]}}];
'''


def _exit_child(I, code=0):
    """Fork[] の子プロセスを終える (バッファを流して os._exit。親と共有する後始末は走らせない: K64 の exit_without_hooks)"""
    for st in (getattr(I, 'out', None), sys.stdout, sys.stderr):
        try:
            st.flush()
        except Exception:
            pass
    os._exit(code)


# ============================================================================ 評価器
_HOLD = {'Fit', 'FindRoot', 'NIntegrate', 'AccelerateParticles', 'Set', 'SetDelayed', 'UpSet', 'UpSetDelayed', 'Unset', 'Do', 'If', 'While', 'For', 'Table', 'Sum', 'Product',
         'Module', 'With', 'Block', 'Which', 'Switch', 'Function', 'Check', 'Timing', 'AppendTo', 'PrependTo',
         'Increment', 'Decrement', 'PreIncrement', 'PreDecrement', 'AddTo', 'SubtractFrom', 'TimesBy', 'DivideBy',
         'CompoundExpression', 'SetAttributes', 'Protect', 'Unprotect', 'Clear', 'Hold', 'HoldForm', 'RuleDelayed',
         'Catch', 'Pattern', 'Blank', 'BlankSequence', 'BlankNullSequence', 'Optional', 'PatternTest', 'Condition',
         'Literal', 'TagSet', 'Return', 'Break', 'Continue', 'Message', 'Definition', 'Slot', 'SlotSequence',
         'FFS$NumericalDerivative', 'Unevaluated', 'SetOptions', 'DefineOptions', 'Class', 'AutoLoad', 'AutoLoadSet', 'WithPriority', 'Priority', 'And', 'Or', 'At', 'StandardForm', 'Context'}


def is_num(v):
    return isinstance(v, (float, int, complex)) and not isinstance(v, bool)


def truthy(v):
    if is_num(v):
        return v != 0
    return None


class _HeldArg:
    """HoldAll / HoldFirst / HoldRest 属性の関数に評価せずに渡した引数 (本体の中で使われるときに評価する)"""
    __slots__ = ('node',)

    def __init__(self, node):
        self.node = node

    def __repr__(self):
        return 'Held[%r]' % (self.node,)


class Frame(dict):
    boundary = False       # True: ユーザ定義関数の呼び出し (定義の中へ値を埋め込むとき、これより外側の呼び出しの引数は見ない)
    immut = False          # True: 関数の引数・With の変数 (定義の中へ値を埋め込める)


import operator as _op
_FOLD_HEADS = {'Plus', 'Times', 'Power', 'Divide', 'Subtract', 'Minus', 'Sin', 'Cos', 'Tan', 'Exp', 'Log', 'Sqrt', 'ArcTan', 'ArcSin', 'ArcCos', 'Abs', 'Sinh', 'Cosh', 'Tanh'}
_FAST2 = {'Plus': _op.add, 'Times': _op.mul, 'Subtract': _op.sub, 'Divide': _op.truediv}


class Interp:
    def __init__(self, out=None):
        self.globals: Dict[str, Any] = {}
        self.frames: List[Frame] = []
        self.downvalues: Dict[str, List[Tuple[list, Any, bool]]] = {}   # name → [(patterns, body, delayed)]
        self.builtins: Dict[str, Callable] = {}
        self.out = out or sys.stdout
        self.files: Dict[int, Any] = {}
        self.next_unit = 10
        self.ffs_handler: Optional[Callable[[str], Any]] = None
        self.ffs_err = 0                    # termes.f の iffserr (FFS コマンドのエラーが出た)
        self.ffs_words: Callable[[str], bool] = lambda w: False
        self.cwd = os.getcwd()
        self.globals.update({'Pi': math.pi, 'E': math.e, 'True': 1.0, 'False': 0.0, 'INF': INF, 'Infinity': INF,
                             'NaN': NAN, 'I': 1j, '$FORM': '', 'PageWidth': 131.0, 'Degree': math.pi / 180,
                             'GoldenRatio': (1 + 5**0.5) / 2, 'EulerGamma': 0.5772156649015329, 'ElectronMass': 510998.95, 'SpeedOfLight': 2.99792458e8,
                             # K64 の物理定数 (gs で確認。ElectronMass は 510998.95 eV)
                             'ElectronCharge': 1.602176634e-19, 'FineStructureConstant': 0.0072973525692838015, 'ElectronRadius': 2.817940326186004e-15,
                             'ProtonMass': 938272088.16, 'ProtonRadius': 1.5346982671812722e-18, 'PlanckConstant': 6.62607015e-34,
                             'PlanckHbar': 1.0545718176461565e-34, 'BoltzmannConstant': 1.380649e-23, 'ElectronGminus2over2': 0.001159652181280002,
                             'SIMu0': 1.2566370621210481e-6, 'MKSAMu0': 1.2566370621210481e-6, 'SIEpsilon0': 8.854187812793e-12, 'MKSAEpsilon0': 8.854187812793e-12})
        self.protected = set()
        self.out_hist = {}                  # Out[n] (最上位の文の値。run_sad_string が echo_out を立てたとき)
        self.out_line = 0
        self.globals['$Line'] = 0.0
        self.rng = None     # SadRandom (ffs 側で設定)
        self.stopped = False
        self.set_hooks = {}
        self.globals['Null'] = Null
        # 文脈 (BeginPackage / Begin)。_ctx_cur は $Context、_ctx_path は $ContextPath、_ctx_syms は文脈つきの記号 ('Stat`x')、
        # _root_syms は利用者のコードが大域 (`) に作った記号。_pkg_loading > 0 の間 (SAD の Packages の読み込み) は文脈を無視する
        self._ctx_cur = '`'
        self._ctx_path = list(_ROOT_PATH)
        self._ctx_stack = []
        self._pkg_stack = []
        self._ctx_syms = set()
        self._root_syms = set()
        self._pkg_loading = 0
        self._ctx_track = False
        self._install_builtins()
        self.run(_PRELUDE, 'script')
        self._ctx_track = True
        self._ctx_publish()

    # ------------------------------------------------------------------ 文脈 (名前空間)
    def _ctx_publish(self):
        self.globals['$Context'] = Sym(self._ctx_cur)
        self.globals['$ContextPath'] = [Sym(c) for c in self._ctx_path]

    def _ctx_known(self, n):
        return (n in self.builtins or n in self.globals or n in self.downvalues or n in self.autoload
                or n in getattr(self, 'attributes', ()))

    def _ctx_is_system(self, n):
        return n in _sysnames() or (n not in self._root_syms and self._ctx_known(n))

    def _ctx_resolve(self, n):
        """文脈の中で読んだ名前 n を、記号の実際の名前にする (K64: 現在の文脈、$ContextPath の順に探し、無ければ現在の文脈に作る)"""
        if n.endswith('`'):
            return n                                  # 文脈名そのもの (Stat`)
        if '`' in n:
            ctx, _, base = n.rpartition('`')
            ctx = _ctx_norm(ctx + '`')
            if ctx in ('`', 'System`'):
                if ctx == '`' and base not in self._root_syms and not self._ctx_known(base):
                    self._root_syms.add(base)         # `x: 文脈の中から大域 (`) に作る記号
                return base
            self._ctx_syms.add(ctx + base)
            return ctx + base
        cur = self._ctx_cur
        if cur == '`':
            if n in self._root_syms:
                return n
        elif cur + n in self._ctx_syms:
            return cur + n
        for c in self._ctx_path:
            if c == cur:
                continue
            if c == '`':
                if n in self._root_syms:
                    return n
            elif c == 'System`':
                if self._ctx_is_system(n):
                    return n
            elif c + n in self._ctx_syms:
                return c + n
        if cur == '`':
            self._root_syms.add(n)
            return n
        self._ctx_syms.add(cur + n)
        return cur + n

    def _ctx_tokens(self, toks):
        active = self._ctx_cur != '`' or self._ctx_path != _ROOT_PATH
        if not active:
            # 文脈を使っていない間は、利用者のコードが新しく作った大域の記号を覚えるだけ (あとで文脈に入ったときに System` の記号と区別する)
            rs = self._root_syms
            for t in toks:
                if t.kind == 'sym':
                    n = t.val
                    if n not in rs and '`' not in n and not self._ctx_known(n):
                        rs.add(n)
                    elif '`' in n and not n.endswith('`'):
                        self._ctx_syms.add(n)
                elif t.kind == 'pat' and t.val[0]:
                    n = t.val[0]
                    if n not in rs and '`' not in n and not self._ctx_known(n):
                        rs.add(n)
            return toks
        out = []
        prev = None
        for t in toks:
            if t.kind == 'sym' and not (prev is not None and prev.kind == 'op' and prev.val in ('::', '?', '@')):
                n2 = self._ctx_resolve(t.val)
                if n2 != t.val:
                    t = Tok('sym', n2, t.line, t.pre)
            elif t.kind == 'pat':
                name, nb, typ = t.val
                t = Tok('pat', (self._ctx_resolve(name) if name else name, nb, self._ctx_resolve(typ) if typ else typ), t.line, t.pre)
            out.append(t)
            prev = t
        return out

    # ------------------------------------------------------------------ 変数
    def lookup(self, name):
        for fr in reversed(self.frames):
            if name in fr:
                return fr[name]
        if name in self.globals:
            return self.globals[name]
        raise KeyError(name)

    def has(self, name):
        for fr in reversed(self.frames):
            if name in fr:
                return True
        return name in self.globals

    def assign(self, name, value):
        skip = False
        for fr in reversed(self.frames):
            if name in fr:
                cur = fr[name]
                if fr.immut and isinstance(cur, Sym):
                    # 変数そのものを渡された引数 (f[Unevaluated$[n]] の i): 代入は渡された変数 n に届く (同名なら、より外側の n)
                    if cur.name != name:
                        self.assign(cur.name, value)
                        return
                    skip = True
                    continue
                fr[name] = value
                return
        hook = self.set_hooks.get(name)
        if hook:
            value = hook(value)
        self.globals[name] = value

    def get(self, name, default=None):
        try:
            return self.lookup(name)
        except KeyError:
            return default

    # ------------------------------------------------------------------ 実行
    def run(self, src: str, mode: str = 'ffs'):
        """mode 'ffs': 文の先頭語が FFS コマンドなら ffs_handler へ。'script': SADScript のみ"""
        result = Null
        global _CUR
        cur_saved = _CUR
        _CUR = self if (getattr(self, '_ctx_track', False) and not self._pkg_loading) else None
        self._run_depth = getattr(self, '_run_depth', 0) + 1
        stmts = self.split_statements(src)
        semis = list(self._last_semis)
        k = 0
        # tffsa.f の REPEAT/UNTIL (levelr, nrpt, irptp) と tfif.f の IF/ELSEIF/ELSE/ENDIF (iflevel) は
        # 文の列を読む側の制御 (ファイル位置の巻き戻し/読み飛ばし) なので、最上位の FFS 文列でだけ働く
        rep_stack = []
        iflevel = [0]
        try:
            while k < len(stmts):
                stmt_src = stmts[k]
                k += 1
                if self.stopped:
                    break
                try:
                    if mode == 'ffs' and self._run_depth == 1:
                        ctl = self._control_word(stmt_src)
                        if ctl is not None:
                            self._ctl_extra = None
                            k = self._control(ctl, stmts, k, rep_stack, iflevel)
                            if self._ctl_extra is None:
                                continue
                            stmt_src, self._ctl_extra = self._ctl_extra, None
                    # FFS コマンドかどうかは実行時に判定する (USE の前後で要素名が変わるため)
                    if mode == 'ffs' and self._is_ffs_command(stmt_src):
                        if self.ffs_handler is None:
                            raise SADError(f'FFS command without handler: {stmt_src}')
                        self.ffs_handler(stmt_src.strip())
                        continue
                    nodes_ = parse(stmt_src)
                    for node in nodes_:
                        if self.stopped:
                            break
                        result = self.eval(node)
                    if getattr(self, 'echo_out', False) and mode == 'ffs' and self._run_depth == 1 and nodes_                             and k - 1 < len(semis) and not semis[k - 1] and stmts[k - 1] is stmt_src:
                        v_ = result
                        if isinstance(nodes_[-1], Expr) and nodes_[-1].head == S('SetDelayed') and len(nodes_[-1].args) == 2:
                            v_ = nodes_[-1].args[1]                 # gs: f[x_] := body の値は未評価の本体
                        if v_ is not Null and not isinstance(v_, Seq):
                            self.out_line += 1
                            self.out_hist[self.out_line] = v_
                            self.globals['$Line'] = float(self.out_line)
                            self.out.write(self.wrap('Out[%d]:= %s' % (self.out_line, self.to_str(v_, None, True))) + chr(10))
                except SADError as e:
                    # SAD は文の実行中のエラーを表示してその文を打ち切り、次の文に進む (最上位の run だけ)
                    if not getattr(self, 'continue_on_error', False) or self._run_depth > 1:
                        raise
                    self.out.write('???' + str(e) + chr(10))
                    if os.environ.get('PYSAD_TRACE'):
                        sys.stderr.write('[run] ' + str(e)[:160] + chr(10) + '    ' + ' '.join(stmt_src.split())[:200] + chr(10))
                    msg = str(e)
                    if msg.startswith('FFS::error') or msg.startswith('General::'):
                        # SAD は打ち切った文を '... in 式' で示す (式の形は SAD の表示と違い、文の先頭 200 文字)
                        src1 = ' '.join(stmt_src.split())
                        self.out.write('... in ' + (src1[:200] + (' ...' if len(src1) > 200 else '')) + chr(10))
                except (_Return, _Break, _Continue, _Throw):
                    # gs: 関数の外の Return / Break / Continue / Throw はエラーでその文を打ち切る
                    if not getattr(self, 'continue_on_error', False) or self._run_depth > 1:
                        raise
                    self.out.write('???General::unexpbreak:  Unexpected Break, Continue, Goto, Return, Throw' + chr(10))
                except (TypeError, IndexError, ValueError, KeyError, ZeroDivisionError, AttributeError,
                        OverflowError) as e:
                    # 組み込み関数への不正な引数 (SAD では General::wrongtype 等で文を打ち切る)。Python 例外も同様に扱う
                    if not getattr(self, 'continue_on_error', False) or self._run_depth > 1:
                        raise
                    if os.environ.get('PYSAD_TRACE'):
                        import traceback
                        traceback.print_exc()
                    self.out.write('???General::wrongtype: ' + type(e).__name__ + ': ' + str(e) + chr(10))
        finally:
            _CUR = cur_saved
            self._run_depth -= 1
            if self._run_depth == 0 and getattr(self, '_forked_child', False):
                _exit_child(self, 0)            # Fork[] の子は、最上位の文が尽きたら終わる (呼び出し元の Python に戻らない)
        return result

    _CTL = re.compile(r'\s*(REP(?:E(?:A(?:T)?)?)?|UNTIL|IF|ELSEIF|ELSE|ENDIF)(?![A-Za-z0-9$_])(.*)$', re.I | re.S)

    def _control_word(self, stmt):
        """文の先頭が FFS の制御語 (REPEAT / UNTIL / IF / ELSEIF / ELSE / ENDIF) なら (語, 残り)。
        tfif.f の delim と同じく、直後が [ * / ) } ] > < = . & ^ | なら関数呼び出しや式なので制御語ではない"""
        m = self._CTL.match(stmt)
        if not m:
            return None
        rest = m.group(2)
        if rest[:1] and rest[:1] in '[*/)}]><=.&^|':
            return None
        w = m.group(1).upper()
        if w.startswith('REP'):
            w = 'REPEAT'
        return w, rest.strip()

    @staticmethod
    def _split_ws(text):
        """空白区切りの語 ({..} [..] (..) "..." の中は分けない)"""
        out, buf, depth, q = [], '', 0, False
        for c in text:
            if q:
                buf += c
                q = c != '"'
                continue
            if c == '"':
                q = True
                buf += c
            elif c in '{[(':
                depth += 1
                buf += c
            elif c in '}])':
                depth -= 1
                buf += c
            elif c.isspace() and depth <= 0:
                if buf:
                    out.append(buf)
                    buf = ''
            else:
                buf += c
        if buf:
            out.append(buf)
        return out

    def _eval_cond(self, text):
        """getva の代わり: text の先頭の式を評価して実数を返す。FFS の式は空白で終わる (QF K1 0.5 のように語が並ぶため)。
        ただし演算子の前後の空白は式の続き (IF KSTEP<= ISTEP)。(値, 残りの語) — 式が無い/数でなければ (None, text)"""
        text = text.strip()
        if not text:
            return None, ''
        words = self._split_ws(text)
        n = 1
        expr = words[0]
        while n < len(words) and (expr[-1] in '+-*/^<>=&|~,@:' or words[n][0] in '<>=&|*/^@,'):
            expr += ' ' + words[n]
            n += 1
        try:
            v = Null
            for nd in parse(expr):
                v = self.eval(nd)
        except Exception:
            return None, text
        if is_num(v):
            return v, ' '.join(words[n:])
        return None, text

    def _control(self, ctl, stmts, k, rep_stack, iflevel):
        """制御語を実行し、次に読む文の添字を返す"""
        word, rest = ctl
        if word == 'REPEAT':
            v, left = self._eval_cond(rest)
            n = int(v) if v is not None else 65535
            rep_stack.append([k, max(1, n)])
            return k
        if word == 'UNTIL':
            if not rep_stack:
                self.out.write(' ???-FFS-Error-?UNTIL without REPEAT. ' + chr(10))
                return k
            v, left = self._eval_cond(rest)
            if v is not None and v != 0:
                rep_stack.pop()
                return k
            rep_stack[-1][1] -= 1
            if rep_stack[-1][1] > 0:
                return rep_stack[-1][0]
            rep_stack.pop()
            return k
        if word == 'ENDIF':
            if iflevel[0] <= 0:
                self.out.write(' ???-FFS-Error-?ENDIF without IF ' + chr(10))
                iflevel[0] = 0
            else:
                iflevel[0] -= 1
            return k
        if word in ('ELSE', 'ELSEIF'):
            # 実行中の分岐の終わりに来た: 対応する ENDIF まで読み飛ばす
            if iflevel[0] <= 0:
                self.out.write(' ???-FFS-Error-?%s without IF ' % word + chr(10))
                iflevel[0] = 0
                return k
            return self._skip_if(stmts, k, iflevel, False)
        # IF
        iflevel[0] += 1
        return self._if_cond(rest, stmts, k, iflevel)

    def _if_cond(self, rest, stmts, k, iflevel):
        v, left = self._eval_cond(rest)
        if v is None:
            self.out.write(' ???-FFS-Error-?Missing condition for IF ' + chr(10))
            iflevel[0] = 0
            return k
        if v != 0:
            if left.strip():
                self._ctl_extra = left        # IF 条件 の後ろに同じ文で続くコマンド/式 (run が続けて実行する)
            return k
        return self._skip_if(stmts, k, iflevel, True, rest)

    def _skip_if(self, stmts, k, iflevel, allow_else, cond_text=None):
        """偽の分岐 (allow_else=True) は ELSE / 同じ深さの ELSEIF (条件を再評価) まで、真の分岐の終わり (False) は ENDIF まで読み飛ばす"""
        ifl = iflevel[0]
        while k < len(stmts) and iflevel[0] >= ifl:
            ctl = self._control_word(stmts[k])
            k += 1
            if ctl is None:
                continue
            word, rest = ctl
            if word == 'ELSE':
                if allow_else and iflevel[0] == ifl:
                    return k
            elif word == 'IF':
                iflevel[0] += 1
            elif word == 'ENDIF':
                iflevel[0] -= 1
            elif word == 'ELSEIF':
                if allow_else and iflevel[0] == ifl:
                    return self._if_cond(rest, stmts, k, iflevel)
        return k

    def split_statements(self, src: str):
        """ソースを文 (文字列) に分ける: 括弧の外の ; と改行で区切る。演算子で終わる行は次の行に続く"""
        src = strip_comments(src)
        out = []
        semis = self._last_semis = []          # 各文が ; で終わっていたか (最上位の Out[n]:= 表示の判定)
        i, n = 0, len(src)
        depth = 0
        stmt_start = 0
        in_str = False
        while i <= n:
            c = src[i] if i < n else chr(10)
            if in_str:
                if c == chr(92):
                    i += 2
                    continue
                if c == '"':
                    in_str = False
                i += 1
                continue
            if c == '"':
                in_str = True
                i += 1
                continue
            if c in '([{':
                depth += 1
            elif c in ')]}':
                depth -= 1
            if depth <= 0 and (c == ';' or c == chr(10) or i == n):
                depth = max(depth, 0)
                stmt = src[stmt_start:i]
                if stmt.strip():
                    if c == chr(10) and (self._incomplete(stmt) or self._wake_table_open(stmt)):
                        i += 1
                        continue
                    out.append(stmt)
                    semis.append(c == ';')
                stmt_start = i + 1
            i += 1
        return out

    @staticmethod
    def _incomplete(stmt: str) -> bool:
        s = stmt.rstrip()
        if not s:
            return False
        if s[-1] == '*' and len(s) > 1 and (s[-2].isalnum() or s[-2] in '$_.}'):
            return False          # "free Q*" のようなワイルドカード (演算子ではない)
        return s[-1] in '+-*/^=,<>&|:@.~' and not s.endswith('=.')

    @staticmethod
    def _wake_table_open(stmt: str) -> bool:
        """v2: 旧式の `WAKE L|T 要素 s1 w1 … ENDWAKE` (K64 の tfwake.f) の表は行をまたいでよい。ENDWAKE (または ;) が来るまで次の行に続く。
        k64 モードでは従来どおり行で切る (K64 は WAKE を WAKECOMMAND[] に渡すだけで、表を読まない)"""
        s = stmt.lstrip()
        if s[:4].upper() != 'WAKE' or not _compat.v2():
            return False
        m = re.match(r'WAKE\s+([A-Za-z]+)(?:\s|$)', s, re.I)
        return (m is not None and ('LONGITUDINAL'.startswith(m.group(1).upper()) or 'TRANSVERSE'.startswith(m.group(1).upper()))
                and re.search(r'(?<![\w$])ENDWAKE(?![\w$])', s, re.I) is None)

    def _is_ffs_command(self, stmt: str) -> bool:
        s = stmt.strip()
        m = re.match(r'([A-Za-z$][A-Za-z0-9$_.*%{}|<>]*)(.*)$', s, re.S)
        if not m:
            return False
        word, rest = m.group(1), m.group(2)
        r = rest.lstrip()
        if r == '*':
            # "disp *" / "free *" の `*` 1 字はワイルドカード (後ろに何もない `x *` は式にならない)
            return bool(self.ffs_words(word))
        if r.startswith('^^^') and self.ffs_words(word):
            return True           # "disp ^^^ 5": ^^^ は先頭の位置 (べき乗の演算子ではない)
        if r[:1] in ('=', '[', '(', ':', '/', '+', '*', '^', '<', '>', '&', '|', '@', ',', '.', '?') and not r.startswith('->'):
            # "x = 1", "f[..", "a/b" は SADScript。ただし "QF K1 =" 形式は SAD にはない
            if r[:2] in ('=.', '==', '=>', '=<') or r[:1] != '=':
                return False
            return False
        if r[:1] == '-' and not r[1:2].isdigit() and not r[1:2] == '.':
            return False
        return bool(self.ffs_words(word))

    # ------------------------------------------------------------------ eval
    def eval(self, node):
        if isinstance(node, (float, int, str, complex, Rule, Func, _Null)) or node is None:
            return node
        if isinstance(node, Sym):
            try:
                v = self.lookup(node.name)
            except KeyError:
                if node.name in self.autoload and node.name not in self.downvalues and node.name not in self.builtins:
                    self._load_autoload(node.name)
                    return self.eval(node)
                if node.name in self.builtins and node.name not in self.downvalues:
                    return node
                # 引数なしの定義 f := ... は downvalues に () パターンで
                dv = self.downvalues.get(node.name)
                if dv:
                    for pats, body, delayed in dv:
                        if pats is None:
                            return self.eval(body) if delayed else body
                return node
            if isinstance(v, Expr) and v.head == S('$Delayed'):
                return self.eval(v.args[0])
            if isinstance(v, Seq) and any(isinstance(q_, _HeldArg) for q_ in v.items):
                # HoldRest の関数の par0__ (保持された引数の列) は、使われるときに呼び出し側の環境で 1 つずつ評価する (Plots.n の FitPlot の {par0})
                idx = None
                for i_ in range(len(self.frames) - 1, -1, -1):
                    if node.name in self.frames[i_]:
                        idx = i_
                        break
                saved = self.frames
                items_ = []
                try:
                    self.frames = saved[:idx] if idx is not None else saved
                    for q_ in v.items:
                        items_.append(self.eval(q_.node) if isinstance(q_, _HeldArg) else q_)
                finally:
                    self.frames = saved
                return Seq(items_)
            if isinstance(v, _HeldArg):
                # 保持された引数は呼び出し側の環境で評価する (この束縛より内側の frame を外す)
                idx = None
                for i_ in range(len(self.frames) - 1, -1, -1):
                    if node.name in self.frames[i_]:
                        idx = i_
                        break
                if idx is None:
                    return self.eval(v.node)
                saved = self.frames
                # 保持引数は呼び出し側の環境で評価するが、関数の中の With[{x = …}, …] の x は保持引数の中の x にも効く
                # (SAD は引数を本体へ文字どおり代入するので、With の変数が引数の式の中の同名の記号も置き換える: Plots.n の PlotInterpolate)
                self.frames = saved[:idx] + [f_ for f_ in saved[idx + 1:] if getattr(f_, 'is_with', False)]
                try:
                    n_ = v.node
                    if any(getattr(f_, 'is_with', False) for f_ in self.frames[idx:]) and _mentions_hold(n_):
                        n_ = self._close(n_)          # Hold[Sin[y]] の中の y にも With の値を埋め込む (SAD は文字どおり代入する。Plot[{Sin[y], Cos[y]}, …])
                    return self.eval(n_)
                finally:
                    self.frames = saved
            if isinstance(v, _VarRef) and any(node.name in fr_ and fr_.immut for fr_ in self.frames):
                # 変数 n そのものを渡された引数 (f[Unevaluated$[n]] の i): 実数・文字列・リストを持つ変数なら、その値として評価する
                # (引数名と変数名が同じ i のときは、引数の束縛より外側の i を見る)
                idx_ = None
                for k_ in range(len(self.frames) - 1, -1, -1):
                    if node.name in self.frames[k_]:
                        idx_ = k_
                        break
                v2 = None
                found_ = False
                if v.name == node.name:
                    for k_ in range(idx_ - 1, -1, -1):
                        if v.name in self.frames[k_]:
                            v2, found_ = self.frames[k_][v.name], True
                            break
                    if not found_ and v.name in self.globals:
                        v2, found_ = self.globals[v.name], True
                else:
                    try:
                        v2, found_ = self.lookup(v.name), True
                    except KeyError:
                        pass
                if found_ and (is_num(v2) or isinstance(v2, (str, list))):
                    return v2
            if isinstance(v, Expr) and v.head != S('$Delayed') and not isinstance(v.head, Expr) and getattr(self, '_expr_reeval', False):
                # 関数の引数に束縛された未評価の式 (Sin[x]) は、本体の中で評価が続く: SAD は引数を本体へ文字どおり代入してから評価するので、
                # With[{x = Vector[#]}, … Sin[x] …] の x が効く (Plots.n の PlotInterpolate の fxin)。評価は呼び出し環境に戻さず現在の環境で行う
                fr_b = None
                with_in = False
                for f_ in reversed(self.frames):
                    if node.name in f_:
                        fr_b = f_
                        break
                    if getattr(f_, 'is_with', False):
                        with_in = True               # その引数を使う場所の内側に With がある (With の変数が引数の式の中の同名の記号に効く場合だけ)
                if fr_b is not None and fr_b.immut and fr_b.boundary and with_in:
                    g_ = self.__dict__.setdefault('_symre', set())
                    key_ = ('E', node.name, id(v))
                    if key_ not in g_:
                        g_.add(key_)
                        try:
                            return self.eval(v)
                        finally:
                            g_.discard(key_)
            if isinstance(v, Sym) and v.name != node.name and (v.name in self.globals or v.name in self.downvalues):
                # a = b; b = 3 のとき a は 3。パターン変数 (c_) が記号 (Canvas$Canvas) に束縛されたときも同じ (評価は続く)。
                # 記号の値は大域の定義で決まるので、評価は局所の束縛を外した環境で行う (f[x_*y_] の y = x が x_ に化けない)
                g_ = self.__dict__.setdefault('_symre', set())
                if v.name in g_:
                    return v
                g_.add(v.name)
                saved_ = self.frames
                self.frames = []
                try:
                    return self.eval(v)
                finally:
                    self.frames = saved_
                    g_.discard(v.name)
            return v
        if isinstance(node, list):
            return node
        if isinstance(node, Expr):
            return self.eval_expr(node)
        return node

    def _install_autoload(self):
        """AutoLoad[x1, x2, ..., "file"] / AutoLoadSet[...] (init.n): 名前が最初に使われたときに file を Get する"""
        I = self

        def _al(args):
            if not args:
                return Null
            path = I.eval(args[-1])
            if not isinstance(path, str):
                return Null
            for a in args[:-1]:
                if isinstance(a, Sym):
                    I.autoload[a.name] = path
            return Null
        self.builtins['AutoLoad'] = _al
        self.builtins['AutoLoadSet'] = _al
        self.builtins['Library'] = lambda a: Null
        self.builtins['Require'] = lambda a: Null

        def _fbound(a):
            nm = a[0].name if isinstance(a[0], Sym) else None
            return 1.0 if nm and (nm in I.builtins or nm in I.downvalues or nm in I.autoload) else 0.0
        self.builtins['FBoundQ$'] = _fbound

        def _names(a):
            import fnmatch
            pat = a[0] if a and isinstance(a[0], str) else '*'
            if '`' in pat:
                # 文脈つきの名前。K64 は Names["Stat`*"] に {"Stat`"} を返す (文脈名だけ)。v2 はその文脈の記号の名前の一覧
                from . import compat
                hit = sorted(n for n in I._ctx_syms if fnmatch.fnmatchcase(n, pat) or fnmatch.fnmatchcase(n, pat.lstrip('`')))
                if compat.v2():
                    return hit
                return sorted({n[:n.rindex('`') + 1] for n in hit})
            allnames = set(I.builtins) | set(I.downvalues) | set(I.globals) | set(I.autoload)
            return sorted(n for n in allnames if fnmatch.fnmatchcase(n, pat))
        self.builtins['Names'] = _names

        def _baseform(a):
            """BaseForm[n, b]: 整数 n の b 進表記の文字列 (0-9, a-z)。Tkinter.n のウィジェットのパス名 (.5 や .1r) に使われる"""
            if len(a) != 2 or not is_num(a[0]) or not is_num(a[1]):
                raise SADError('General::wrongtype: Argument must be Real number in BaseForm')
            n, base = int(a[0]), int(a[1])
            if base < 2 or base > 36:
                raise SADError('General::wrongtype: Base must be 2..36 in BaseForm')
            digs = '0123456789abcdefghijklmnopqrstuvwxyz'
            if n == 0:
                return '0'
            sign = '-' if n < 0 else ''
            n = abs(n)
            out = ''
            while n:
                n, r = divmod(n, base)
                out = digs[r] + out
            return sign + out
        self.builtins['BaseForm'] = _baseform

    def _load_autoload(self, name):
        path = self.autoload.get(name)
        if path is None:
            return
        for k in [k for k, v in self.autoload.items() if v == path]:
            del self.autoload[k]
        self.load_package(path)

    def load_package(self, path):
        """SAD のパッケージ (.n) を SADScript として 1 文ずつ評価する。文のエラーは (SAD と同じく) 表示して次の文へ進む。
        環境変数 PYSAD_TRACE があれば、エラーになった文の先頭を stderr に出す (パッケージを動かすための調査用)"""
        if not os.path.exists(path):
            return
        src = open(path, encoding='utf-8', errors='replace').read()
        depth = getattr(self, '_run_depth', 0)
        self._run_depth = depth + 1
        saved_frames = self.frames
        self.frames = []      # 関数の途中で自動読み込みされても、呼び出し元の引数・局所変数をパッケージに見せない
        global _CUR
        cur_saved = _CUR
        _CUR = None           # SAD の Packages は文脈を無視して全部大域で読む (従来どおり)
        self._pkg_loading += 1
        try:
            for st in self.split_statements(src):
                try:
                    for node in parse(st):
                        self.eval(node)
                except (SADError, TypeError, IndexError, ValueError, KeyError, ZeroDivisionError, AttributeError, OverflowError, RecursionError) as e:
                    if os.environ.get('PYSAD_TRACE'):
                        sys.stderr.write('[load %s] %s: %s\n    %s\n' % (os.path.basename(path), type(e).__name__, str(e)[:160],
                                                                         ' '.join(st.split())[:200]))
        finally:
            self._pkg_loading -= 1
            _CUR = cur_saved
            self._run_depth = depth
            self.frames = saved_frames

    def eval_args(self, args):
        out = []
        for a in args:
            v = self.eval(a)
            if isinstance(v, Seq):
                out.extend(v.items)
            else:
                out.append(v)
        return out

    def eval_expr(self, node: Expr):
        head = node.head
        if isinstance(head, Sym):
            name = head.name
            if name in _HOLD:
                f = self.builtins.get(name)
                if f is None:
                    raise SADError(f'unsupported hold form {name}')
                return f(node.args)
            if name == 'List':
                return self.eval_args(node.args)
            if self.frames and name not in self.builtins and name not in self.downvalues:
                # 関数の引数・With の変数が Symbol に束縛されていれば、その Symbol が head (a[$Window] の a が w など)
                for fr_ in reversed(self.frames):
                    if name in fr_:
                        v_ = fr_[name]
                        if isinstance(v_, Sym) and fr_.immut:
                            head = v_
                            name = v_.name
                            node = Expr(v_, node.args)
                        break
            if name in self.autoload and name not in self.downvalues and name not in self.builtins:
                self._load_autoload(name)
            # ユーザ定義
            dv = self.downvalues.get(name)
            at = self.attributes.get(name)
            held_inert = bool(at) and not dv and name not in self.builtins and ('HoldAll' in at or 'HoldFirst' in at or 'HoldRest' in at)
            if at and (dv or held_inert):
                hold_all = 'HoldAll' in at
                hold_first = 'HoldFirst' in at
                hold_rest = 'HoldRest' in at
                args = []
                for i, a in enumerate(node.args):
                    if hold_all or (hold_first and i == 0) or (hold_rest and i > 0):
                        h_ = self._close(a) if self._has_immut() else a
                        if isinstance(h_, _Lit) and isinstance(h_.v, Expr):
                            h_ = h_.v
                        if isinstance(h_, Expr) and h_.head == S('Unevaluated') and len(h_.args) == 1:
                            h_ = h_.args[0]                   # 保持引数の位置では Unevaluated[x] の外側を外す (Widget[Unevaluated[Root$Window], Null])
                        elif isinstance(h_, Expr) and h_.head == S('Evaluate') and len(h_.args) == 1:
                            h_ = self.eval(h_.args[0])        # 保持引数の位置の Evaluate[x] は評価する (ParametricPlot[Evaluate[ReleaseHold[a]], …])
                        args.append(_HeldArg(h_))
                    else:
                        v = self.eval(a)
                        if isinstance(v, Seq):
                            args.extend(v.items)
                        else:
                            args.append(v)
            else:
                if dv and any(isinstance(a_, Expr) and a_.head == S('Unevaluated$') and len(a_.args) == 1 for a_ in node.args):
                    # f[Unevaluated$[i]]: i を評価せずに (変数 i そのものを) 渡す。関数の中の i++ や w[i] = … が呼び出し側の i を進める (KBF$ArrangeComponents)
                    args = []
                    for a_ in node.args:
                        if isinstance(a_, Expr) and a_.head == S('Unevaluated$') and len(a_.args) == 1:
                            u_ = a_.args[0]
                            if isinstance(u_, Sym):
                                u_ = _VarRef(u_.name)
                            if isinstance(u_, Expr) and self._has_immut():
                                # f[Unevaluated$[g[y, x]]]: 呼び出し側の関数引数 y x の値を埋め込んで、使われるときに評価する (D.n の ReduceExpression[Unevaluated$[D$[y,x]]])
                                u_ = _HeldArg(self._subst_inert(u_))
                            args.append(u_)
                        else:
                            v_ = self.eval(a_)
                            if isinstance(v_, Seq):
                                args.extend(v_.items)
                            else:
                                args.append(v_)
                else:
                    args = self.eval_args(node.args)
            if dv:
                r = self.apply_downvalues(name, dv, args)
                if r is not _NoMatch:
                    return r
            if len(args) == 2:
                op = _FAST2.get(name)
                if op is not None:
                    a0 = args[0]
                    a1 = args[1]
                    if type(a0) is float and type(a1) is float and not (name == 'Divide' and a1 == 0.0):
                        return op(a0, a1)         # 実数どうしの四則 (lmap2/_fold を経由しない近道。結果は同じ)
            f = self.builtins.get(name)
            if f is not None:
                return f(args)
            # 値を持つシンボル (リスト → Part、純関数 → 適用)
            try:
                hv = self.lookup(name)
            except KeyError:
                if dv and any(d[0] is None for d in dv):
                    hv2 = self.eval(head)             # 遅延した値を持つ変数 (fr := Widget[fr, ...]) を f[x] の f に置く: Widget[fr, ...][x] として適用する
                    if hv2 != head and (isinstance(hv2, (Expr, Func)) or isinstance(hv2, list)):
                        return self.apply_value(hv2, args, head)       # fun := Sqrt[#]& (純関数の遅延した値) も f[x] で適用する
                if at:
                    args = [x.node if isinstance(x, _HeldArg) else x for x in args]    # 保持された引数は式のまま返す (Widget[w, ...] など)
                return Expr(head, args)
            if isinstance(hv, Expr) and hv.head == S('$Delayed'):
                hv = self.eval(hv.args[0])            # Module[{h}, h := Sqrt[#]&; h[4]]: 局所変数の遅延した値
            return self.apply_value(hv, args, head)
        # head が式 (f[a][b], (f@g)[x], Function[...][x])
        if isinstance(head, Expr) and head.head == S('At') and len(head.args) == 2 and isinstance(head.args[1], Sym):
            ms_ = self.resolve_member_head(head)
            if ms_ is not None:
                return self.eval_expr(Expr(ms_, node.args))            # x@method[args]: インスタンス x のメソッド
            return self.eval(Expr(head.args[0], [Expr(head.args[1], node.args)]))     # f@g[x] は f[g[x]]
        if isinstance(head, Expr) and head.head == S('Unevaluated$') and len(head.args) == 1:
            head = head.args[0]                   # Unevaluated$[tf1][Background] は tf1[Background]
            if isinstance(head, Sym):
                return self.eval_expr(Expr(head, node.args))
        hv = self.eval(head)
        args = self.eval_args(node.args)
        return self.apply_value(hv, args, head)

    def apply_value(self, hv, args, head=None):
        if hv is Null:
            return Seq(list(args))                # Hold[x] /. Hold -> Null は Null[x] で、Null[a, b] は a, b を並べる Sequence
        if isinstance(hv, Expr) and hv.head == S('Unevaluated$') and len(hv.args) == 1:
            hv = self.eval(hv.args[0])            # Unevaluated$[tf1][Background] は tf1[Background]
            if isinstance(hv, Sym):
                return self.eval_expr(Expr(hv, [_Lit(x) for x in args]))
        if isinstance(hv, Expr) and isinstance(hv.head, Expr) and isinstance(hv.head.head, Sym) and hv.head.head.name in self.subvalues3:
            for g0, g1, g2, body, delayed in self.subvalues3[hv.head.head.name]:
                b = {}
                if self._match_seq(g0, list(hv.head.args), b) and self._match_seq(g1, list(hv.args), b) and self._match_seq(g2, list(args), b):
                    if not delayed:
                        return self.eval(self.subst_bindings(body, b)) if b else body
                    fr_ = Frame(b)
                    fr_.immut = True
                    fr_.boundary = True
                    self.frames.append(fr_)
                    try:
                        return self.eval(body)
                    except _Return as r:
                        return r.value
                    finally:
                        self.frames.pop()
        if isinstance(hv, Expr) and isinstance(hv.head, Sym) and hv.head.name in self.subvalues:
            for hp, ap, body, delayed in self.subvalues[hv.head.name]:
                b = {}
                if self._match_seq(hp, list(hv.args), b) and self._match_seq(ap, list(args), b):
                    if not delayed:
                        return self.eval(self.subst_bindings(body, b)) if b else body
                    fr_ = Frame(b)
                    fr_.immut = True
                    fr_.boundary = True
                    self.frames.append(fr_)
                    try:
                        return self.eval(body)
                    except _Return as r:
                        return r.value
                    finally:
                        self.frames.pop()
        if isinstance(hv, list):
            return self.part(hv, args)
        if isinstance(hv, str):
            # 文字列の添字 (SAD は文字列の Part を許す)
            if len(args) == 2 and all(is_num(x) for x in args):
                # s[i, j]: 文字列の部分 (1 始まり、負は末尾から。SAD の文字列の Part)
                n = len(hv)
                i, j = int(args[0]), int(args[1])
                i = n + 1 + i if i < 0 else i
                j = n + 1 + j if j < 0 else j
                return hv[max(i, 1) - 1:max(j, 0)]
            if args and any(isinstance(x, str) for x in args):
                return Expr(hv, args)        # "K1"["QD"] (ElementValues のキー) は評価せず式のまま (gs)
            return self.part(list(hv), args) if args else hv
        if isinstance(hv, Func):
            return self.apply_func(hv, args)
        if isinstance(hv, Sym):
            dv = self.downvalues.get(hv.name)
            if dv:
                r = self.apply_downvalues(hv.name, dv, args)
                if r is not _NoMatch:
                    return r
            f = self.builtins.get(hv.name)
            if f is not None:
                return f(args)
            return Expr(hv, args)
        if isinstance(hv, Expr):
            return Expr(hv, args)
        if is_num(hv) and not args:
            return hv
        return Expr(head if head is not None else hv, args)

    # ------------------------------------------------------------------ 純関数・ユーザ定義
    def apply_func(self, fn: Func, args):
        if fn.params is not None:
            fr = Frame()
            fr.immut = True
            for p, a in zip(fn.params, args):
                fr[p] = a
            self.frames.append(fr)
            try:
                return self.eval(fn.body)           # Return は純関数を抜けて、外側のユーザ定義関数まで届く (gs で確認)
            finally:
                self.frames.pop()
        body = self.subst_slots(fn.body, args)
        return self.eval(body)

    def subst_slots(self, node, args):
        if isinstance(node, Expr):
            if node.head == S('Slot'):
                k = int(node.args[0]) if node.args else 1
                if 1 <= k <= len(args):
                    return _Lit(args[k - 1])
                return Null
            if node.head == S('SlotSequence'):
                return Expr(S('Sequence'), [_Lit(a) for a in args])
            if node.head == S('Function') or node.head == S('Hold'):
                return node       # 内側の純関数・Hold の中は置換しない (gs: (Hold[f[#]]&)[3] は Hold[f[#]])
            return Expr(self.subst_slots(node.head, args), [self.subst_slots(a, args) for a in node.args])
        if isinstance(node, list) and any(isinstance(a, (Expr, list)) for a in node):
            return [self.subst_slots(a, args) for a in node]
        return node

    def apply_downvalues(self, name, dv, args):
        for pats, body, delayed in dv:
            if pats is None:
                continue
            b = self.match_args(pats, args)
            if b is None:
                continue
            if not delayed:
                return self.eval(self.subst_bindings(body, b)) if b else body
            fr_ = Frame(b)
            fr_.immut = True
            fr_.boundary = True
            self.frames.append(fr_)
            try:
                return self.eval(body)
            except _Return as r:
                return r.value
            finally:
                self.frames.pop()
        return _NoMatch

    def _subst_inert(self, node):
        """Unevaluated$[g[y, x]] の y x に、いま束縛されている値を埋め込む。値は評価しない (_Lit)。
        SAD は呼び出しのたびに本体の記号を値に置き換えるので、値の中の記号 (y に x*y を入れたときの x) が別の引数名に捕まらない"""
        if isinstance(node, Sym):
            for fr in reversed(self.frames):
                if node.name in fr:
                    if not fr.immut:
                        return node
                    v = fr[node.name]
                    if isinstance(v, _HeldArg):
                        return v.node
                    if isinstance(v, Expr) and v.head == S('$Delayed'):
                        return node
                    return _Lit(v)
            return node
        if isinstance(node, Expr):
            return Expr(self._subst_inert(node.head) if isinstance(node.head, Expr) else node.head, [self._subst_inert(a) for a in node.args])
        return node

    def _has_immut(self):
        return any(fr.immut for fr in self.frames)

    _xs = frozenset()          # _close が置き換えない名前 (内側の With / Module / Function / パターン変数が束縛し直したもの)

    @staticmethod
    def _declared_names(decl):
        names = set()
        items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
        for it in items:
            if isinstance(it, Expr) and it.head == S('CompoundExpression'):
                items = list(items) + list(it.args)
                continue
            tgt = it.args[0] if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')) and it.args else it
            if isinstance(tgt, Sym):
                names.add(tgt.name)
            elif isinstance(tgt, Expr) and tgt.head == S('List'):
                for x in tgt.args:
                    if isinstance(x, Sym):
                        names.add(x.name)
        return names

    @staticmethod
    def _pattern_var_names(x, acc):
        if isinstance(x, Expr):
            if x.head == S('Pattern') and x.args and isinstance(x.args[0], Sym):
                acc.add(x.args[0].name)
            for a in x.args:
                Interp._pattern_var_names(a, acc)
            if isinstance(x.head, Expr):
                Interp._pattern_var_names(x.head, acc)
        elif isinstance(x, Rule):
            Interp._pattern_var_names(x.lhs, acc)
            Interp._pattern_var_names(x.rhs, acc)
        elif isinstance(x, list):
            for a in x:
                Interp._pattern_var_names(a, acc)
        return acc

    def _close(self, node, head=False):
        """定義の中へ、いま見えている関数引数・With の変数の値を埋め込む (head=True なら値が Symbol のときだけ)"""
        if isinstance(node, Sym):
            if node.name in self._xs:
                return node
            for fr in reversed(self.frames):
                if node.name not in fr and fr.boundary:
                    break                           # 別の関数の呼び出しの引数は見ない (f[x_] := … が呼び出し元の x に化けない)
                if node.name in fr:
                    if not fr.immut:
                        return node
                    v = fr[node.name]
                    if isinstance(v, (Expr,)) and v.head == S('$Delayed'):
                        return self._close(v.args[0]) if not head else node     # With[{g:=e/2}, f[e_]:=g*3]: 定義の中へは式 (評価前) を埋め込む
                    if isinstance(v, _HeldArg):
                        return v.node
                    if isinstance(v, (Sym, Expr)):
                        return v                    # 未評価の式 (a[$Menu$KEK] などウィジェットを表す式) はそのまま埋め込む
                    return node if head else _Lit(v)
            return node
        if isinstance(node, Expr):
            if node.head == S('Pattern') and node.args:
                return Expr(node.head, [node.args[0]] + [self._close(a) for a in node.args[1:]])
            hn_ = node.head.name if isinstance(node.head, Sym) else None
            if hn_ in ('With', 'Module', 'Block') and len(node.args) >= 2 and not self._xs.issuperset(self._declared_names(node.args[0])):
                # 内側の宣言が同じ名前を束縛し直す: その名前は本体の中では置き換えない (宣言の右辺は外側の環境)
                decl = node.args[0]
                names = self._declared_names(decl)
                if isinstance(decl, Expr) and decl.head == S('List'):
                    nd = []
                    for it in decl.args:
                        if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')) and len(it.args) == 2:
                            nd.append(Expr(it.head, [it.args[0], self._close(it.args[1])]))
                        else:
                            nd.append(it)
                    decl = Expr(decl.head, nd)
                old_xs = self._xs
                self._xs = old_xs | names
                try:
                    body = [self._close(a) for a in node.args[1:]]
                finally:
                    self._xs = old_xs
                return Expr(node.head, [decl] + body)
            if hn_ in ('Set', 'SetDelayed', 'UpSet', 'UpSetDelayed') and len(node.args) == 2:
                pv = self._pattern_var_names(node.args[0], set())
                if pv and not self._xs.issuperset(pv):
                    lhs = self._close(node.args[0], True) if isinstance(node.args[0], Expr) or isinstance(node.args[0], Sym) else node.args[0]
                    old_xs = self._xs
                    self._xs = old_xs | pv
                    try:
                        rhs = self._close(node.args[1])
                    finally:
                        self._xs = old_xs
                    return Expr(node.head, [lhs, rhs])
            if hn_ == 'Function' and len(node.args) == 2:
                pn = self._declared_names(node.args[0] if isinstance(node.args[0], Expr) and node.args[0].head == S('List') else Expr(S('List'), [node.args[0]]))
                if pn:
                    old_xs = self._xs
                    self._xs = old_xs | pn
                    try:
                        body = self._close(node.args[1])
                    finally:
                        self._xs = old_xs
                    return Expr(node.head, [node.args[0], body])
            new_args = []
            for a_ in node.args:
                c_ = self._close(a_)
                if isinstance(c_, _Lit) and isinstance(c_.v, Seq):
                    for it_ in c_.v.items:           # opt___ の引数列は展開して埋め込む
                        new_args.append(it_.node if isinstance(it_, _HeldArg) else (it_ if isinstance(it_, (Sym, Expr)) else _Lit(it_)))
                else:
                    new_args.append(c_)
            return Expr(self._close(node.head, True) if isinstance(node.head, Sym) else self._close(node.head), new_args)
        if isinstance(node, Rule):
            return Rule(self._close(node.lhs), self._close(node.rhs), node.delayed) if hasattr(node, 'delayed') else node
        return node

    def _mentions(self, node, b):
        """値のリスト node の中に、束縛 b の名前の記号があるか (大きな数表を毎回作り直さないため)"""
        stack = [node]
        while stack:
            x = stack.pop()
            if isinstance(x, Sym):
                if x.name in b:
                    return True
            elif isinstance(x, list):
                stack.extend(y for y in x if not isinstance(y, (float, int, str)))
            elif isinstance(x, Expr):
                stack.append(x.head)
                stack.extend(x.args)
        return False

    def subst_bindings(self, node, b):
        if isinstance(node, Sym) and node.name in b:
            return b[node.name].node if isinstance(b[node.name], _HeldArg) else _Lit(b[node.name])
        if isinstance(node, Expr):
            return Expr(self.subst_bindings(node.head, b), [self.subst_bindings(a, b) for a in node.args])
        if isinstance(node, list):
            # 評価した右辺 (G[x_] = Round[{x}*3]) は値のリストを含む。その中の記号 x も置き換えるので、List 式にして評価し直す
            if not self._mentions(node, b):
                return node
            return Expr(S('List'), [self.subst_bindings(a, b) for a in node])
        return node

    # ---- パターン照合
    def match_args(self, pats, args) -> Optional[dict]:
        b = {}
        if self._match_seq(pats, list(args), b):
            return b
        return None

    def _match_seq(self, pats, args, b) -> bool:
        if not pats:
            return not args
        p = pats[0]
        kind, name, typ, default = self._pat_info(p)
        if kind in ('seq', 'nullseq'):
            # 残りのパターン数だけ後ろに残して貪欲に取る
            rest = len(pats) - 1
            nmin = 1 if kind == 'seq' else 0
            inner_ = None
            q_ = p.args[1] if (isinstance(p, Expr) and p.head == S('Pattern') and len(p.args) == 2) else p
            pinner_ = None
            if isinstance(q_, Expr) and q_.head in (S('Repeated'), S('RepeatedNull')) and q_.args:
                if not _is_pattern(q_.args[0]) and not (isinstance(q_.args[0], Expr) and q_.args[0].head == S('Alternatives')):
                    inner_ = q_.args[0]                      # {0 ..}: 要素がすべて 0 のリスト
                else:
                    pinner_ = q_.args[0]                     # (_Rule|_RuleDelayed)... : 要素がすべてこのパターンに合う
            for take in range(len(args) - rest, nmin - 1, -1):
                if take < nmin:
                    break
                chunk = args[:take]
                if typ and not all(self._type_ok(typ, a) for a in chunk):
                    continue
                if inner_ is not None and not all(self.same(inner_, a) for a in chunk):
                    continue
                if pinner_ is not None and not all(self._match_one(pinner_, x.v if isinstance(x, _Lit) else x, {}) for x in chunk):
                    continue
                b2 = dict(b)
                if name:
                    if name in b2 and not (isinstance(b2[name], Seq) and b2[name].items == chunk):
                        continue
                    b2[name] = Seq([x.v if isinstance(x, _Lit) else x for x in chunk])      # 保持された式の中の _Lit は値に戻して束縛する
                if self._match_seq(pats[1:], args[take:], b2):
                    b.clear()
                    b.update(b2)
                    return True
            return False
        if not args:
            if kind == 'optional':
                b2 = dict(b)
                if name:
                    b2[name] = default
                if self._match_seq(pats[1:], args, b2):
                    b.clear()
                    b.update(b2)
                    return True
            return False
        a = args[0]
        if (isinstance(p, Expr) and (p.head == S('Plus') or p.head == S('Times')) and isinstance(a, Expr) and a.head == p.head):
            # Orderless の引数の入れ替え方のうち、後ろのパターン (D$[x_*y_, x_] の 2 つ目の x_) とも合うものを探す
            for b3 in self._orderless_alts(p.args, [x.v if isinstance(x, _Lit) else x for x in a.args], b):
                if self._match_seq(pats[1:], args[1:], b3):
                    b.clear()
                    b.update(b3)
                    return True
            return False
        b2 = dict(b)
        if self._match_one(p, a, b2) and self._match_seq(pats[1:], args[1:], b2):
            b.clear()
            b.update(b2)
            return True
        if kind == 'optional':
            b2 = dict(b)
            if name:
                b2[name] = default
            if self._match_seq(pats[1:], args, b2):
                b.clear()
                b.update(b2)
                return True
        return False

    def _pat_info(self, p):
        """(kind, name, type, default): kind ∈ blank/seq/nullseq/optional/literal/test"""
        if isinstance(p, Expr):
            h = p.head
            if h == S('Pattern'):
                name = p.args[0].name if isinstance(p.args[0], Sym) else None
                k, _, typ, d = self._pat_info(p.args[1])
                return k, name, typ, d
            if h == S('Blank'):
                return 'blank', None, (p.args[0].name if p.args else None), None
            if h == S('BlankSequence'):
                return 'seq', None, (p.args[0].name if p.args else None), None
            if h == S('BlankNullSequence'):
                return 'nullseq', None, (p.args[0].name if p.args else None), None
            if h in (S('Repeated'), S('RepeatedNull')):
                k, name, typ, d = self._pat_info(p.args[0])
                return ('seq' if h == S('Repeated') else 'nullseq'), name, typ, d
            if h == S('Optional'):
                k, name, typ, _ = self._pat_info(p.args[0])
                return 'optional', name, typ, self.eval(p.args[1]) if len(p.args) > 1 else Null
            if h == S('PatternTest'):
                k, name, typ, d = self._pat_info(p.args[0])
                return k, name, typ, d
        return 'literal', None, None, None

    def _type_ok(self, typ, a):
        if typ is None:
            return True
        if typ == 'Real':
            return is_num(a) and not isinstance(a, complex)
        if typ == 'Complex':
            return isinstance(a, complex)
        if typ == 'String':
            return isinstance(a, str)
        if typ == 'List':
            return isinstance(a, list) or (isinstance(a, Expr) and a.head == S('List'))     # 保持された (評価前の) {Sin[y], Cos[y]} も List
        if typ == 'Symbol':
            return isinstance(a, Sym)
        if typ == 'Rule':
            return isinstance(a, Rule)
        if isinstance(a, Expr) and isinstance(a.head, Sym):
            return a.head.name == typ
        return False

    def _const_fold(self, node):
        """保持された引数が数と定数 (Pi・E・Degree) と組み込みの数値関数だけの式なら、その値 (gs は 2*Pi、Sin[1.]、1+2. を _Real に照合する)。そうでなければ None"""
        def ok(x):
            if is_num(x):
                return True
            if isinstance(x, Sym):
                return x.name in ('Pi', 'E', 'Degree')
            if isinstance(x, Expr) and isinstance(x.head, Sym) and x.args:
                n_ = x.head.name
                return n_ in self.builtins and n_ not in self.downvalues and n_ in _FOLD_HEADS and all(ok(y) for y in x.args)
            return False
        if isinstance(node, Expr) and ok(node):
            try:
                v = self.eval(node)
            except (SADError, ArithmeticError, ValueError, TypeError):
                return None
            return v if is_num(v) else None
        if isinstance(node, Sym) and ok(node):
            return self.eval(node)
        return None

    def _match_held(self, p, a, b) -> bool:
        """a は評価されていない引数 (_HeldArg)。リテラルは式として比べ、構造パターンは式の構造で照合する (束縛される値は保持された式)"""
        node = a.node.v if isinstance(a.node, _Lit) else a.node
        if isinstance(p, Expr) and not is_num(node):
            pt_ = p.args[1] if (p.head == S('Pattern') and len(p.args) == 2) else p
            if isinstance(pt_, Expr) and pt_.head == S('Blank') and pt_.args:
                fv_ = self._const_fold(node)
                if fv_ is not None:
                    a = _HeldArg(fv_)           # x_Real は定数式 (2*Pi など) の値に照合し、値を束縛する
                    node = fv_
        if isinstance(p, _Lit):
            return self.same(p.v, node) or repr(p.v) == repr(node)
        if isinstance(p, Expr):
            h = p.head
            if h == S('Pattern'):
                name = p.args[0].name
                if not self._match_one(p.args[1], a, b):
                    return False
                if name in b and not (isinstance(b[name], _HeldArg) and repr(b[name].node) == repr(node)):
                    return False
                b[name] = a
                return True
            if h in (S('Blank'), S('BlankSequence'), S('BlankNullSequence')):
                return self._type_ok(p.args[0].name if p.args else None, node)
            if h == S('Optional'):
                return self._match_one(p.args[0], a, b)
            if h == S('PatternTest'):
                if not self._match_one(p.args[0], a, b):
                    return False
                try:
                    return truthy(self.eval(Expr(p.args[1], [_Lit(self.eval(node))]))) is True
                except SADError:
                    return False
            if h == S('Alternatives'):
                for q in p.args:
                    b2 = dict(b)
                    if self._match_one(q, a, b2):
                        b.update(b2)               # 一致した選択肢の束縛を残す (Widget[master_,_] | master_)
                        return True
                return False
            if (isinstance(node, Expr) and (node.head == h or (isinstance(h, Expr) and _is_pattern(h) and self._match_one(h, node.head, b))) and (len(node.args) == len(p.args) or any(_is_seq_pat(q) for q in p.args))):
                return self._match_seq(p.args, [_HeldArg(c) for c in node.args], b)
            if isinstance(node, list) and h == S('List') and (len(node) == len(p.args) or any(_is_seq_pat(q) for q in p.args)):
                return self._match_seq(p.args, [_HeldArg(c) for c in node], b)
            return False
        if isinstance(p, Sym):
            return isinstance(node, Sym) and node.name == p.name
        return self.same(p, node)

    def _orderless_alts(self, pargs, aargs, b):
        """Plus・Times の引数の入れ替え (数を先頭にした標準順が先、続いて 6 個までの順列) ごとの束縛 (照合に成功したもの)"""
        import itertools
        cand = sorted(aargs, key=lambda x: 0 if is_num(x) else 1)
        orders = [cand]
        # 列の途中に素の Blank (x_) と列パターン (m__) が両方あるときは標準順だけ (K64 の x_*(m__)+x_*(n__) の共通因子の括り出しは標準順の先頭どうしだけ)
        kinds = [self._pat_info(q)[0] for q in pargs]
        bare = [k_ == 'blank' and (isinstance(q, Expr) and (q.head == S('Blank') or (q.head == S('Pattern') and len(q.args) == 2 and isinstance(q.args[1], Expr) and q.args[1].head == S('Blank') and not q.args[1].args)))
                for k_, q in zip(kinds, pargs)]
        restrict = any(k_ in ('seq', 'nullseq') for k_ in kinds) and any(bare)
        if len(aargs) <= 6 and not restrict:
            seen = {tuple(map(id, cand))}
            for perm in itertools.permutations(cand):
                key = tuple(map(id, perm))
                if key not in seen:
                    seen.add(key)
                    orders.append(list(perm))
        for od in orders:
            b2 = dict(b)
            if self._match_seq(pargs, od, b2):
                yield b2

    def _match_orderless(self, pargs, aargs, b) -> bool:
        """Plus・Times (Orderless) の照合: 数を先頭にした標準順で試し、だめなら引数の順列 (6 個まで) を試す"""
        import itertools
        cand = sorted(aargs, key=lambda x: 0 if is_num(x) else 1)
        b0 = dict(b)
        if self._match_seq(pargs, list(cand), b):
            return True
        b.clear()
        b.update(b0)
        if len(aargs) > 6:
            return False
        seen = {tuple(map(id, cand))}
        for perm in itertools.permutations(cand):
            key = tuple(map(id, perm))
            if key in seen:
                continue
            seen.add(key)
            b2 = dict(b0)
            if self._match_seq(pargs, list(perm), b2):
                b.clear()
                b.update(b2)
                return True
        return False

    def _match_one(self, p, a, b) -> bool:
        if isinstance(a, _HeldArg):
            return self._match_held(p, a, b)
        if isinstance(p, _Lit):
            if self.same(p.v, a):
                return True
            if isinstance(a, Sym) and isinstance(p.v, Expr) and not any(d_[0] is None for d_ in (self.downvalues.get(a.name) or ())):
                # 記号 (Canvas$Canvas) で来た引数は、その値でも照合する (UpValue 内の `m[c] = …` は c を評価した値で保存される)。
                # 遅延した定義 (The$CanvasObject := …) を持つ記号は評価しない (評価すると代入の左辺の定義が走ってしまう)
                try:
                    va = self.eval(a)
                except Exception:
                    va = a
                if not (isinstance(va, Sym) and va.name == a.name):
                    return self.same(p.v, va) or repr(p.v) == repr(va)
            return False
        if isinstance(p, Rule):
            # Rule[_String, _] のように評価されて Rule オブジェクトになったパターン: a->b の両辺を照合する
            if not isinstance(a, Rule) or bool(a.delayed) != bool(p.delayed):
                return False
            return self._match_one(p.lhs, a.lhs, b) and self._match_one(p.rhs, a.rhs, b)
        if isinstance(p, list):
            # リストのパターン {_, _} (評価済みのリスト。要素にパターンを含む)
            if not isinstance(a, list):
                return False
            if len(a) != len(p) and not any(_is_seq_pat(q) for q in p):
                return False
            return self._match_seq(p, list(a), b)
        if isinstance(p, Expr):
            h = p.head
            if h == S('Pattern'):
                name = p.args[0].name
                if not self._match_one(p.args[1], a, b):
                    return False
                if name in b and b[name] != a:
                    return False
                b[name] = a.v if isinstance(a, _Lit) else a
                return True
            if h in (S('Blank'), S('BlankSequence'), S('BlankNullSequence')):
                return self._type_ok(p.args[0].name if p.args else None, a)
            if h == S('Optional'):
                return self._match_one(p.args[0], a, b)
            if h == S('PatternTest'):
                if not self._match_one(p.args[0], a, b):
                    return False
                try:
                    t = self.eval(Expr(p.args[1], [_Lit(a)]))
                except SADError:
                    return False                   # gs: 検査の中のエラー (Second[List] など) は不一致 (メッセージなし)
                return truthy(t) is True
            if h == S('Alternatives'):
                for q in p.args:
                    b2 = dict(b)
                    if self._match_one(q, a, b2):
                        b.update(b2)               # 一致した選択肢の束縛を残す (Widget[master_,_] | master_)
                        return True
                return False
            if h == S('List'):
                if isinstance(a, Expr) and a.head == S('List'):
                    a = list(a.args)                  # 評価されていない {x, y} (a[Protocol] := {"WM_DELETE_WINDOW", cb} の右辺) も List として照合する
                if not isinstance(a, list):
                    return False
                if len(a) != len(p.args) and not any(_is_seq_pat(q) for q in p.args):
                    return False
                return self._match_seq(p.args, list(a), b)
            if isinstance(a, Rule) and len(p.args) == 2:
                # a->b / a:>b は Rule[a, b] / RuleDelayed[a, b] と同じ形として照合する。head は Rule, RuleDelayed, (Rule|RuleDelayed) など
                nm_ = S('RuleDelayed') if a.delayed else S('Rule')
                if h == nm_ or (isinstance(h, Expr) and self._match_one(h, nm_, dict(b))):
                    return self._match_seq(p.args, [a.lhs, a.rhs], b)
                return False
            if (isinstance(a, Expr) and (a.head == h or (isinstance(h, Expr) and _is_pattern(h) and self._match_one(h, a.head, b))) and (len(a.args) == len(p.args) or any(_is_seq_pat(q) for q in p.args))):
                aa_ = [x.v if isinstance(x, _Lit) else x for x in a.args]
                if h == S('Plus') or h == S('Times'):
                    return self._match_orderless(p.args, aa_, b)           # Orderless: 引数の順を入れ替えて照合 (D.n の x_*y_、b___*Sin[x_]*Cos[x_])
                return self._match_seq(p.args, aa_, b)         # wid_[master_, opt___] のように head もパターン (保持引数の _Lit は値に戻す)
            return False
        if isinstance(p, Sym):
            if isinstance(a, Sym):
                return a.name == p.name
            # パターン中のシンボルは値を持てば値と比較
            try:
                return self.same(self.lookup(p.name), a)
            except KeyError:
                return False
        return self.same(p, a)

    @staticmethod
    def same(a, b):
        if is_num(a) and is_num(b):
            return a == b
        if type(a) != type(b):
            return False
        if isinstance(a, list):
            return len(a) == len(b) and all(Interp.same(x, y) for x, y in zip(a, b))
        return a == b

    # ------------------------------------------------------------------ Part
    def part(self, lst, idx):
        v = lst
        for ci_, k in enumerate(idx):
            if k is Null:                         # l[[, 2]]: 空の添字は全要素 (All)
                if isinstance(v, Expr):
                    v = list(v.args)
                if not isinstance(v, list):
                    raise SADError(f'Part: not a list ({v!r})')
                rest_ = list(idx[ci_ + 1:])
                return [self.part(x, rest_) if rest_ else x for x in v]
            if isinstance(k, list):
                rest_ = list(idx[ci_ + 1:])        # l[[{i,j}, k]] = {l[[i,k]], l[[j,k]]} (残りの添字は各要素に)
                return [self.part(v, [j] + rest_) for j in k]
            if k == 0 and isinstance(v, Expr):
                return v.head                 # f[a][[0]] は頭部 f (Plots.n の Extract[?contl, {…,0,1}])
            if isinstance(v, Rule):
                v = [v.lhs, v.rhs]            # Rule[a, b] の第 1/2 成分 (#[[1]] で名前、#[[2]] で値)
            elif isinstance(v, Expr):
                v = list(v.args)              # f[a,b][[2]] = b (BeamLine[...][[i]] など)
            if not isinstance(v, (list, str)):
                raise SADError(f'Part: not a list ({v!r})')
            i = int(k)
            if i == 0:
                return S('List')
            if i < 0:
                i = len(v) + i + 1
            if not 1 <= i <= len(v):
                if os.environ.get('PYSAD_TB'):
                    import traceback
                    print('PART', repr(lst)[:300], idx, flush=True)
                    traceback.print_stack(limit=8)
                raise SADError(f'Part: index {i} out of range (length {len(v)})')
            v = v[i - 1]
        return v

    def set_part(self, name, idx, value):
        lst = self.lookup(name)
        if not isinstance(lst, list):
            raise SADError(f'Part assignment to non-list {name}')
        v = lst
        for k in idx[:-1]:
            i = int(k)
            if i < 0:
                i = len(v) + i + 1
            v = v[i - 1]
        i = int(idx[-1])
        if i < 0:
            i = len(v) + i + 1
        v[i - 1] = value
        return value

    # ------------------------------------------------------------------ 文字列化
    def form(self):
        f = self.get('$FORM', '')
        return f if isinstance(f, str) else ''

    def to_str(self, v, form=None, quote=False) -> str:
        if form is None:
            form = self.form()
        if isinstance(v, str):
            return f'"{v}"' if quote else v
        if isinstance(v, bool):
            return '1' if v else '0'
        if is_num(v):
            return fmt_real(v, form)
        if isinstance(v, list):
            return '{' + ','.join(self.to_str(x, form, True) for x in v) + '}'
        if isinstance(v, _Null):
            # SAD: Print[Null] / ToString[Null] は "Null"。リストや関数の引数の中 (quote=True の再帰呼出し) では空
            return '' if quote else 'Null'
        if isinstance(v, Sym):
            n = v.name
            if '`' in n and self._ctx_path != _ROOT_PATH and not n.endswith('`') and n[:n.rindex('`') + 1] in self._ctx_path:
                return n[n.rindex('`') + 1:]          # gs: $ContextPath にある文脈の記号は、文脈名を付けずに表示する
            return n
        if isinstance(v, Rule):
            return self.to_str(v.lhs, form, True) + (':>' if v.delayed else '->') + self.to_str(v.rhs, form, True)
        if isinstance(v, Expr):
            hn = v.head.name if isinstance(v.head, Sym) else None
            if hn == 'Plus' and len(v.args) >= 2:
                return '+'.join(self.to_str(a, form, True) for a in v.args)
            if hn == 'Times' and len(v.args) >= 2:
                if v.args[0] == -1.0:
                    return '(-' + ' '.join(self._paren_plus(a, form) for a in v.args[1:]) + ')'       # Times[-1, x, y] は (-x y)
                return '(' + ' '.join(self._paren_plus(a, form) for a in v.args) + ')'
            if hn == 'Power' and len(v.args) == 2:
                base_ = self._paren_plus(v.args[0], form)
                if isinstance(v.args[0], Expr) and v.args[0].head == S('Times') and len(v.args[0].args) >= 2:
                    base_ = '(' + base_ + ')'              # K64: Power の底の Times は括弧が 2 重 ((a x))^3
                return base_ + '^' + self._paren_plus(v.args[1], form)
            if hn == 'Blank' and len(v.args) <= 1:
                return '_' + (self.to_str(v.args[0], form, True) if v.args else '')
            if hn in _CMPSYM and len(v.args) == 2:
                return self.to_str(v.args[0], form, True) + _CMPSYM[hn] + self.to_str(v.args[1], form, True)
            if hn == 'Dot' and len(v.args) == 2:
                return self.to_str(v.args[0], form, True) + '.' + self.to_str(v.args[1], form, True)
            if hn == 'Minus' and len(v.args) == 1:
                return '(-' + self._paren_plus(v.args[0], form) + ')'
            return self.to_str(v.head, form, True) + '[' + ','.join(self.to_str(a, form, True) for a in v.args) + ']'
        if isinstance(v, Func):
            return repr(v)
        return str(v)

    def _paren_plus(self, a, form):
        t = self.to_str(a, form, True)
        if isinstance(a, Expr) and a.head == S('Plus'):
            return '(' + t + ')'
        return t

    def wrap(self, s: str) -> str:
        """PageWidth で折り返す (SAD は ', ' の後で改行し、続きの行は 1 空白から)"""
        try:
            pw = int(self.get('PageWidth', 131.0))
        except Exception:
            pw = 131
        if len(s) <= pw or '\n' in s:
            return s
        lines = []
        cur = ''
        parts = s.split(',')
        for k, p in enumerate(parts):
            seg = p + (',' if k < len(parts) - 1 else '')
            if cur and len(cur.rstrip()) + len(seg.rstrip()) > pw:
                lines.append(cur.rstrip())
                cur = ' ' + seg
            else:
                cur += seg
        lines.append(cur)
        return '\n'.join(lines)

    def write_unit(self, unit, text):
        if unit == 6 or unit is None or unit == -1:        # $Output = -1
            self.out.write(text)
            self.out.flush()
        else:
            f = self.files.get(int(unit))
            if f is None:
                raise SADError(f'Write: unit {unit} is not open')
            f.write(text)

    # ------------------------------------------------------------------ 組み込み関数
    def _install_builtins(self):
        B = self.builtins
        I = self

        def num(v, what='argument'):
            if is_num(v):
                return v
            raise SADError(f'{what}: number expected, got {I.to_str(v)}')

        # ---- 代入系 (hold)
        def _set(args, delayed=False):
            lhs, rhs = args[0], args[1] if len(args) > 1 else Null
            if isinstance(lhs, _Lit) and isinstance(lhs.v, (Sym, Expr)):
                lhs = lhs.v                       # 引数を代入した形 (x → cbv[0]) の代入先
            if isinstance(lhs, Expr) and lhs.head == S('Unevaluated$') and len(lhs.args) == 1:
                lhs = lhs.args[0]                 # Unevaluated$[tf1] = … は tf1 = … (With[{wsym = Extract[wv, {1,2}, Unevaluated$]}, wsym = …])
            elif isinstance(lhs, Expr) and isinstance(lhs.head, Expr) and lhs.head.head == S('Unevaluated$') and len(lhs.head.args) == 1:
                lhs = Expr(lhs.head.args[0], lhs.args)    # Unevaluated$[w][Create$Line] = … は w[Create$Line] = … (CanvasW[Create$Line] = {…})
            if isinstance(lhs, Expr) and lhs.head == S('At') and len(lhs.args) == 2 and isinstance(lhs.args[1], Sym):
                ms_ = I.resolve_member_head(lhs)
                if ms_ is not None:
                    lhs = ms_                         # x@v = 3 はインスタンス x のメンバー v への代入
            elif isinstance(lhs, Expr) and isinstance(lhs.head, Expr) and lhs.head.head == S('At'):
                ms_ = I.resolve_member_head(lhs.head)
                if ms_ is not None:
                    lhs = Expr(ms_, lhs.args)         # x@c[Background] = v
            if isinstance(lhs, Expr) and lhs.head == S('TagSet') and len(lhs.args) == 2:
                lhs = lhs.args[1]                 # a/:f[a, x] := … は f[a, x] := … と同じ (タグは定義の置き場なので無視)
            if I._has_immut():
                # 関数の引数・With の変数が Symbol のとき、代入先は (Mathematica の代入と同じく) その Symbol
                if isinstance(lhs, Sym):
                    for fr_ in reversed(I.frames):
                        if lhs.name in fr_:
                            v_ = fr_[lhs.name]
                            if fr_.immut and not isinstance(v_, (Sym, Expr, _HeldArg)) and not delayed:
                                nv_ = _copy(I.eval(rhs))           # 値で渡された引数への代入 (w = Second[wl]) は、その引数 (局所) を書き換える。同名の大域変数には触れない
                                fr_[lhs.name] = nv_
                                return nv_
                            break
                    lhs = I._close(lhs, True)
                elif isinstance(lhs, Expr) and isinstance(lhs.head, Sym):
                    nh = I._close(lhs.head, True)
                    if nh is not lhs.head:
                        lhs = Expr(nh, lhs.args)
                elif isinstance(lhs, Expr) and isinstance(lhs.head, Expr) and not _is_pattern(lhs.head):
                    nh = I._close(lhs.head)               # a[$Window][Title] = … の head (a[$Window]) の中の a も、その Symbol に置き換える
                    if nh is not lhs.head:
                        lhs = Expr(nh, lhs.args)
            if isinstance(lhs, Expr) and lhs.head == S('Unevaluated$') and len(lhs.args) == 1:
                lhs = lhs.args[0]                 # 引数 w が Unevaluated$[tf1] に束縛されている: w = … は tf1 = …
            if delayed and I._has_immut() and isinstance(lhs, (Sym, Expr)):
                # 関数・With の中の遅延定義 (a$ := Widget[a$, w[master$]]): 本体へ現在の引数値を埋め込む。
                # 定義自身のパターン変数 (f[x_] := … の x) は外側の同名の変数に置き換えない (Canvas.n の CoordLimitedX[x_] := Restrict[x, …])
                pn_ = set()

                def _pn(e_):
                    if isinstance(e_, Expr):
                        if e_.head == S('Pattern') and e_.args and isinstance(e_.args[0], Sym):
                            pn_.add(e_.args[0].name)
                        _pn(e_.head)
                        for a_ in e_.args:
                            _pn(a_)
                    elif isinstance(e_, list):
                        for a_ in e_:
                            _pn(a_)
                _pn(lhs)
                old_xs_ = I._xs
                I._xs = old_xs_ | pn_
                try:
                    rhs = I._close(rhs)
                finally:
                    I._xs = old_xs_
            if (not delayed and isinstance(lhs, Expr) and isinstance(lhs.head, Sym) and lhs.head.name == 'Shared' and len(lhs.args) == 1
                    and not _is_pattern(lhs) and 'WriteShared' in I.builtins):
                # shared.n: (Shared[n_]=v_)^:=WriteShared[n,v] (++・+= なども _set を通る)
                return I.builtins['WriteShared']([I.eval(lhs.args[0]), I.eval(rhs)])
            # `(lhs = x_) ^:= body` で定義された UpValue (Tkinter.n のウィジェットの `w[Text]="abc"` や Tk$Echo=0 など)。
            # Set / SetDelayed の downvalue として保持してあるので、通常の代入の前に引く (第 1 引数は評価せずに照合する)
            if (isinstance(lhs, Expr) and isinstance(lhs.head, Sym) and not _is_pattern(lhs) and lhs.args and lhs.head.name not in ('Part', 'List', 'Symbol', 'Set', 'SetDelayed')
                    and not I.attributes.get(lhs.head.name) and I.downvalues.get('SetDelayed' if delayed else 'Set')
                    and any(isinstance(a_, Sym) and I.has(a_.name) for a_ in lhs.args)):
                lhs = Expr(lhs.head, [I.eval(a_) for a_ in lhs.args])        # f[x] = v は x を評価してから UpValue を探す (PkgConfig$Default[pkg, var] = val)
            if isinstance(lhs, Expr) and isinstance(lhs.head, Expr) and not _is_pattern(lhs) and (I.downvalues.get('SetDelayed' if delayed else 'Set')):
                lhs = Expr(lhs.head, [a_.v if isinstance(a_, _Lit) else I.eval(a_) for a_ in lhs.args])      # a$[op$] のように変数を引数に持つ f[x][op] = v: 引数は評価してから UpValue を探す
            if (not delayed and isinstance(lhs, Sym) and lhs.name in I.autoload and lhs.name not in I.downvalues and lhs.name not in I.builtins
                    and not any(lhs.name in fr_ for fr_ in I.frames)):
                I._load_autoload(lhs.name)            # Canvas$Widget = c1: Canvas.n の `(Canvas$Widget=Widget[c_,_Canvas])^:=…` (UpValue) も読み込み後でないと引けない
            if (not delayed and isinstance(rhs, Expr) and isinstance(rhs.head, Sym) and rhs.head.name in I.autoload
                    and rhs.head.name not in I.downvalues and rhs.head.name not in I.builtins
                    and os.path.basename(I.autoload.get(rhs.head.name) or '') in ('correction.n', 'bump.n', 'corutility.n')):
                I._load_autoload(rhs.head.name)       # b = Bump[…]: パッケージの `(a_ = Bump[…]) ^:= …` (UpValue) は読み込み後でないと引けない
            dvh = I.downvalues.get('SetDelayed' if delayed else 'Set')
            if dvh and isinstance(lhs, (Sym, Expr)) and not _is_pattern(lhs):
                if not delayed:
                    rhs = _Lit(I.eval(rhs))
                r = I.apply_downvalues('SetDelayed' if delayed else 'Set', dvh, [lhs, I.eval(rhs) if not delayed else rhs])
                if r is not _NoMatch:
                    return r
            if isinstance(lhs, Sym):
                v = rhs if delayed else I.eval(rhs)
                if delayed and any(lhs.name in fr_ and not fr_.immut for fr_ in I.frames):
                    I.assign(lhs.name, Expr(S('$Delayed'), [rhs]))        # Module の局所変数への h := … は、その局所の遅延した値 (同名の大域の定義は触らない)
                    return Null
                if delayed:
                    I.downvalues.setdefault(lhs.name, [])
                    I.downvalues[lhs.name] = [d for d in I.downvalues[lhs.name] if d[0] is not None]
                    I.downvalues[lhs.name].insert(0, (None, rhs, True))
                    I.globals.pop(lhs.name, None)
                    return Null
                if not any(lhs.name in fr_ for fr_ in I.frames):
                    I.downvalues.pop(lhs.name, None)          # 局所変数 (Module の w など) への代入で、同名の大域変数の定義を消さない
                I.assign(lhs.name, _copy(v))
                return v
            if isinstance(lhs, Expr) and I._has_immut() and _is_pattern(lhs) or (isinstance(lhs, Expr) and isinstance(lhs.head, Sym) and lhs.head.name in ('Set', 'SetDelayed') and I._has_immut()):
                # 関数・With の中で定義される規則 (UpValue を含む): 左辺のパターンと遅延の右辺に現在の値を埋め込む
                lhs = I._close(lhs, True) if isinstance(lhs, Expr) else lhs
                if delayed:
                    rhs = I._close(rhs)
            if isinstance(lhs, Expr):
                h = lhs.head
                if h == S('List') and not delayed:
                    v = I.eval(rhs)
                    if not isinstance(v, list) or len(v) != len(lhs.args):
                        if os.environ.get('PYSAD_TB'):
                            print('LSET', lhs, repr(v)[:300], flush=True)
                        raise SADError('Set: list assignment length mismatch')
                    for l1, v1 in zip(lhs.args, v):
                        _set([l1, _Lit(v1)])
                    return v
                if h == S('Part'):
                    base = lhs.args[0]
                    idx = I.eval_args(lhs.args[1:])
                    v = I.eval(rhs)
                    return I.set_part(base.name, idx, _copy(v))
                if h == S('Symbol') and len(lhs.args) == 1:
                    nm = I.eval(lhs.args[0])
                    if not isinstance(nm, str):
                        raise SADError('General::wrongtype: Argument must be Character-string in Symbol')
                    v = rhs if delayed else I.eval(rhs)
                    I.assign(nm, _copy(v))
                    return v
                if isinstance(h, Sym) and (I.has(h.name) and isinstance(I.get(h.name), list)) and not delayed \
                        and all(not _is_pattern(a) for a in lhs.args):
                    idx = I.eval_args(lhs.args)
                    v = I.eval(rhs)
                    return I.set_part(h.name, idx, _copy(v))
                if isinstance(h, Sym):
                    # 関数定義 f[x_, ...] := body  または f[literal] = value
                    if h.name in ('Set', 'SetDelayed') and lhs.args:
                        # UpValue の定義: 第 1 引数 (代入先) は評価せずリテラルとして照合する
                        pats = [lhs.args[0] if _is_pattern(lhs.args[0]) else _Lit(lhs.args[0])] + [a if _is_pattern(a) else _Lit(I.eval(a)) for a in lhs.args[1:]]
                    else:
                        at_ = I.attributes.get(h.name)
                        if at_ and ('HoldAll' in at_ or 'HoldFirst' in at_ or 'HoldRest' in at_):
                            pats = []
                            for i_, a in enumerate(lhs.args):
                                held = 'HoldAll' in at_ or ('HoldFirst' in at_ and i_ == 0) or ('HoldRest' in at_ and i_ > 0)
                                if held and not _is_pattern(a) and I._has_immut():
                                    a = I._close(a)            # 保持される位置の変数 (w など) は、いま束縛されている式に置き換えて鍵にする
                                pats.append(a if _is_pattern(a) else (_Lit(a.v if isinstance(a, _Lit) else a) if held else _Lit(I.eval(a))))
                        else:
                            pats = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in lhs.args]
                    dv = I.downvalues.setdefault(h.name, [])
                    dv[:] = [d for d in dv if not (d[0] is not None and _same_pats(d[0], pats))]
                    body = rhs if delayed else I.eval(rhs)
                    # 具体的なパターン (リテラル) を先に
                    nlit = _spec_score(pats)
                    pos = 0
                    for pos, d in enumerate(dv):
                        if d[0] is not None and _spec_score(d[0]) < nlit:
                            break
                    else:
                        pos = len(dv)
                    if getattr(I, 'priority', False):
                        pos = 0                       # WithPriority[...] の中の定義は先頭に置く (init.n の $Priority)
                    dv.insert(pos, (pats, body, delayed))
                    I.globals.pop(h.name, None)
                    return Null if delayed else body          # f[x] = v の値は v (WidgetName[x]=... を本体に持つ関数の戻り値)
                if isinstance(h, Expr) and isinstance(h.head, Sym):
                    # SubValues: Widget[b_,_][c_] := …、ReturnToSAD[$Window][Title] = … (f[a][b] の形の定義)
                    nm = h.head.name
                    hp = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in h.args]
                    ap = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in lhs.args]
                    body = rhs if delayed else I.eval(rhs)
                    lst = I.subvalues.setdefault(nm, [])
                    lst[:] = [e for e in lst if not (_same_pats(e[0], hp) and _same_pats(e[1], ap))]
                    nlit = sum(1 for q in hp + ap if isinstance(q, _Lit))
                    pos = len(lst)
                    for k, e in enumerate(lst):
                        if sum(1 for q in e[0] + e[1] if isinstance(q, _Lit)) < nlit:
                            pos = k
                            break
                    lst.insert(pos, (hp, ap, body, delayed))
                    return Null
                if isinstance(h, Expr) and isinstance(h.head, Expr) and isinstance(h.head.head, Sym):
                    # 3 段の SubValues: cf[n][c][Background] = … (KBFrame.n の部品の設定)
                    nm = h.head.head.name
                    g0 = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in h.head.args]
                    g1 = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in h.args]
                    g2 = [a if _is_pattern(a) else _Lit(I.eval(a)) for a in lhs.args]
                    body = rhs if delayed else I.eval(rhs)
                    lst = I.subvalues3.setdefault(nm, [])
                    lst[:] = [e for e in lst if not (_same_pats(e[0], g0) and _same_pats(e[1], g1) and _same_pats(e[2], g2))]
                    lst.insert(0, (g0, g1, g2, body, delayed))
                    return Null
                if isinstance(h, Expr):
                    raise SADError(f'Set: unsupported lhs {lhs!r}')
            raise SADError(f'Set: unsupported lhs {lhs!r}')
        B['Set'] = lambda a: _set(a, False)
        B['SetDelayed'] = lambda a: _set(a, True)
        B['UpSet'] = lambda a: _set(a, False)
        B['UpSetDelayed'] = lambda a: _set(a, True)
        B['TagSet'] = lambda a: _set(a[1:], False)

        def _unset(args):
            lhs = args[0]
            if I._has_immut() and isinstance(lhs, Sym):
                lhs = I._close(lhs, True)
            if isinstance(lhs, Sym) and lhs.name in I.globals and isinstance(I.globals[lhs.name], I._sad_instance_type):
                I.destroy_instance(I.globals[lhs.name])     # x=. はインスタンスの解体 (Destructor を呼び、メンバーを消す)
            if isinstance(lhs, Sym):
                for fr in reversed(I.frames):
                    if lhs.name in fr:
                        del fr[lhs.name]
                        return Null
                I.globals.pop(lhs.name, None)
                I.downvalues.pop(lhs.name, None)
            elif isinstance(lhs, Expr) and isinstance(lhs.head, Sym):
                # f[a, b] =. : その定義だけを消す (以前は f の定義を全部消していた)。widget.n の `a[$Window]=.` など
                nm = lhs.head.name
                dv = I.downvalues.get(nm)
                if dv:
                    at_ = I.attributes.get(nm) or ()
                    pats = []
                    for i_, a in enumerate(lhs.args):
                        held = 'HoldAll' in at_ or ('HoldFirst' in at_ and i_ == 0) or ('HoldRest' in at_ and i_ > 0)
                        if _is_pattern(a):
                            pats.append(a)
                        else:
                            if held and I._has_immut():
                                a = I._close(a)
                            pats.append(_Lit(a.v if isinstance(a, _Lit) else a) if held else _Lit(I.eval(a)))
                    keep = [d for d in dv if not (d[0] is not None and _same_pats(d[0], pats))]
                    if len(keep) == len(dv) and not any(_is_pattern(q) for q in pats):
                        pass                      # 該当なし
                    I.downvalues[nm] = keep
                    if not keep:
                        I.downvalues.pop(nm, None)
            elif isinstance(lhs, Expr) and isinstance(lhs.head, Expr) and isinstance(lhs.head.head, Sym):
                # contl[k][x] =. : 添字つきの値 (SubValues) のその定義だけを消す (Plots.n の等高線の端点の表)
                nm = lhs.head.head.name
                lst = I.subvalues.get(nm)
                if lst:
                    hp = [a if _is_pattern(a) else _Lit(a.v if isinstance(a, _Lit) else I.eval(a)) for a in lhs.head.args]
                    ap = [a if _is_pattern(a) else _Lit(a.v if isinstance(a, _Lit) else I.eval(a)) for a in lhs.args]
                    lst[:] = [e for e in lst if not (_same_pats(e[0], hp) and _same_pats(e[1], ap))]
            return Null
        B['Unset'] = _unset
        B['Clear'] = lambda a: [_unset([x]) for x in a] and Null

        def _compound(args):
            r = Null
            for a in args:
                r = I.eval(a)
            return r
        B['CompoundExpression'] = _compound

        def _inc(args, d, pre):
            lhs = args[0]
            if isinstance(lhs, Sym) and I._has_immut():
                lhs = I._close(lhs, True)         # i++ の i が変数 n を指す引数 (Unevaluated$[n]) なら n を進める
            old = I.eval(lhs)
            new = num(old) + d
            _set([lhs, _Lit(new)])
            return new if pre else old
        B['Increment'] = lambda a: _inc(a, 1.0, False)
        B['Decrement'] = lambda a: _inc(a, -1.0, False)
        B['PreIncrement'] = lambda a: _inc(a, 1.0, True)
        B['PreDecrement'] = lambda a: _inc(a, -1.0, True)

        def _opassign(args, op):
            lhs = args[0]
            old = I.eval(lhs)
            v = I.eval(args[1])
            new = op(old, v)
            _set([lhs, _Lit(new)])
            return new
        B['AddTo'] = lambda a: _opassign(a, lambda x, y: I.builtins['Plus']([x, y]))
        B['SubtractFrom'] = lambda a: _opassign(a, lambda x, y: I.builtins['Subtract']([x, y]))
        B['TimesBy'] = lambda a: _opassign(a, lambda x, y: I.builtins['Times']([x, y]))
        B['DivideBy'] = lambda a: _opassign(a, lambda x, y: I.builtins['Divide']([x, y]))

        def _appendto(args, front=False):
            lhs = args[0]
            cur = I.eval(lhs)
            v = I.eval(args[1])
            if not isinstance(cur, list):
                cur = []
            new = ([v] + cur) if front else (cur + [v])
            _set([lhs, _Lit(new)])
            return new
        B['AppendTo'] = lambda a: _appendto(a)
        B['PrependTo'] = lambda a: _appendto(a, True)

        # ---- 制御 (hold)
        def _iter_values(spec):
            # {i, a, b, d}: 先頭の i は評価しない (Module の局所変数 i があっても反復変数として束縛する)
            if isinstance(spec, Expr) and spec.head == S('List'):
                if len(spec.args) >= 2 and isinstance(spec.args[0], Sym):
                    spec = [spec.args[0]] + [I.eval(x) for x in spec.args[1:]]
                else:
                    spec = [I.eval(x) for x in spec.args]
            else:
                spec = I.eval(spec)
            if not isinstance(spec, list):
                if num(spec) == INF:
                    return None, itertools.repeat(None)
                return None, [None] * int(num(spec))
            if len(spec) == 1:
                if num(spec[0]) == INF:
                    return None, itertools.repeat(None)
                return None, [None] * int(num(spec[0]))
            var = spec[0]
            name = var.name if isinstance(var, Sym) else None
            if len(spec) == 2:
                b = spec[1]
                if isinstance(b, list):
                    # SAD は Table/Do/Sum で {i, {a,b,c}} 形式を受け付けない (General::wrongtype)
                    raise SADError('General::wrongtype: Argument must be {symbol, n}, {symbol, start, stop}, {symbol, start, stop, step}')
                if num(b) == INF:
                    return name, itertools.count(1.0)             # Do[…, {i, Infinity}]・Do[…, {Infinity}]: Break するまで続ける
                return name, [float(k) for k in range(1, int(math.floor(num(b) + 1e-9)) + 1)]
            a = num(spec[1])
            b = num(spec[2])
            d = num(spec[3]) if len(spec) > 3 else 1.0
            if d == 0:
                raise SADError('Do: zero step')
            n = int(math.floor((b - a) / d + 1e-9)) + 1
            return name, [a + k * d for k in range(max(n, 0))]

        def _do(args):
            body = args[0]
            specs = args[1:]

            def loop(k):
                if k == len(specs):
                    I.eval(body)
                    return
                name, vals = _iter_values(specs[k])
                fr = Frame()
                I.frames.append(fr)
                try:
                    for v in vals:
                        if name:
                            fr[name] = v
                        try:
                            loop(k + 1)
                        except _Continue:
                            continue
                finally:
                    I.frames.pop()
            try:
                loop(0)
            except _Break:
                pass
            return Null
        B['Do'] = _do

        def _table(args):
            body = args[0]
            specs = args[1:]
            if not specs:
                return [I.eval(body)]

            def build(k):
                name, vals = _iter_values(specs[k])
                fr = Frame()
                I.frames.append(fr)
                out = []
                try:
                    for v in vals:
                        if name:
                            fr[name] = v
                        r_ = build(k + 1) if k + 1 < len(specs) else _copy(I.eval(body))
                        if isinstance(r_, Seq):
                            out.extend(r_.items)          # Table[If[…, Null[], i], …]: 空の Sequence は取り除き、Null[a,b] は a, b を並べる
                        else:
                            out.append(r_)
                finally:
                    I.frames.pop()
                return out
            return build(0)
        B['Table'] = _table

        def _sum(args, prod=False):
            """Sum/Product[body, {i,..}, {j,..}, ...]: 最初の反復子が外側 (Table と同じ)。項が無ければ 0 / 1"""
            body = args[0]
            specs = args[1:]
            op = I.builtins['Times' if prod else 'Plus']

            def rec(k):
                name, vals = _iter_values(specs[k])
                fr = Frame()
                I.frames.append(fr)
                acc = 1.0 if prod else 0.0
                try:
                    for v in vals:
                        if name:
                            fr[name] = v
                        x = rec(k + 1) if k + 1 < len(specs) else I.eval(body)
                        acc = op([acc, x])
                finally:
                    I.frames.pop()
                return acc
            if not specs:
                raise SADError('General::narg: Sum/Product needs an iterator')
            return rec(0)
        B['Sum'] = _sum
        B['Product'] = lambda a: _sum(a, True)

        def _if(args):
            c = I.eval(args[0])
            t = truthy(c)
            if t is True:
                return I.eval(args[1]) if len(args) > 1 else Null
            if t is False:
                return I.eval(args[2]) if len(args) > 2 else Null
            return I.eval(args[3]) if len(args) > 3 else Null
        B['If'] = _if

        def _while(args):
            try:
                while truthy(I.eval(args[0])) is True:
                    try:
                        if len(args) > 1:
                            I.eval(args[1])
                    except _Continue:
                        pass
            except _Break:
                pass
            return Null
        B['While'] = _while

        def _for(args):
            from . import compat
            v2 = compat.v2()
            if len(args) != 4:
                if v2 and len(args) == 3:
                    args = list(args) + [Null]         # v2: For[ini,cond,incr] は本体なしのループ (K64 は未評価のまま)
                else:
                    return Expr(S('For'), list(args))      # gs: For[ini_,cond_,incr_,body_] は 4 引数だけ (3 引数は未評価のまま)
            I.eval(args[0])
            try:
                while truthy(I.eval(args[1])) is True:
                    try:
                        if v2:                         # K64 (gs) は 4 番目の引数 (本体) を実行せず incr だけを繰り返す (B71)
                            I.eval(args[3])
                    except _Continue:
                        pass
                    I.eval(args[2])
            except _Break:
                pass
            return Null
        B['For'] = _for

        def _compat_v2(args):
            from . import compat
            return 1.0 if compat.v2() else 0.0
        B['Compat$V2'] = _compat_v2

        def _which(args):
            if len(args) % 2 or not args:
                raise SADError('General::narg:  Number of arguments is expected to be even number in Which[%s]' % ','.join(I.to_str(a) for a in args))
            for k in range(0, len(args) - 1, 2):
                if truthy(I.eval(args[k])) is True:
                    return I.eval(args[k + 1])
            return Expr(S('Which'), list(args))        # gs: 一致なしは式をそのまま返す
        B['Which'] = _which

        def _switch(args):
            x = I.eval(args[0])
            for k in range(1, len(args) - 1, 2):
                p = args[k]
                if isinstance(p, Sym):
                    pv_ = I.eval(p)
                    if isinstance(pv_, Expr) and _is_pattern(pv_):
                        p = pv_                       # 変数が持つパターン (Switch[val, types, ...] の types = _Real)
                if (isinstance(p, Expr) and p.head in (S('Blank'),) and not p.args) or I.same(I.eval(p), x) or I._match_one(p, x, {}):
                    return I.eval(args[k + 1])
            if len(args) % 2 == 0:
                raise SADError('General::narg:  Number of arguments is expected to be odd number in Switch[%s]' % ','.join(I.to_str(a) for a in args))
            # gs: 一致なしは Switch[値, 規則…] を返す。引数がすべて定数なら読み込み時に畳まれて Switch[値] だけになる
            if all(_is_const(a) for a in args):
                return Expr(S('Switch'), [x])
            return Expr(S('Switch'), [x] + list(args[1:]))
        B['Switch'] = _switch

        def _resolve(nm):
            # 関数の引数 z が記号 x に束縛されているとき、Block[{z}, …]・Module[{z}, …]・With[{z = …}, …] は x を局所にする
            # (SAD は引数を本体へ文字どおり代入するので宣言の z も x になる)
            for fr0_ in reversed(I.frames):
                if nm in fr0_:
                    v0_ = fr0_[nm]
                    if isinstance(v0_, _HeldArg) and isinstance(v0_.node, Sym):
                        v0_ = v0_.node                      # HoldAll の関数 (Plot) の x_Symbol は Held[y]
                    if fr0_.immut and isinstance(v0_, Sym):
                        return v0_.name
                    break
            return nm

        def _module(args):
            fr = Frame()
            decl = args[0]
            items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
            inits = []
            for it in items:
                if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')):
                    inits.append((it.args[0], it.args[1], it.head == S('SetDelayed')))
                elif isinstance(it, Sym):
                    inits.append((it, None, False))
            I.frames.append(fr)
            try:
                for lhs, init, delayed in inits:
                    if isinstance(lhs, Expr) and lhs.head == S('List'):
                        vals = I.eval(init)
                        for sym, v in zip(lhs.args, vals if isinstance(vals, list) else [vals]):
                            fr[sym.name] = _copy(v)
                        continue
                    name = _resolve(lhs.name)
                    if init is None:
                        fr[name] = Sym(name)          # 未代入のローカル変数は Symbol 自身として評価される
                    elif delayed:
                        fr[name] = Expr(S('$Delayed'), [init])
                    else:
                        fr[name] = _copy(I.eval(init))
                r_ = Null
                for body_ in args[1:]:          # Module[vars, e1, e2, …] は e1, e2, … を順に評価して最後の値 (Plots.n の FitPlot)
                    r_ = I.eval(body_)
                return r_
            finally:
                I.frames.pop()
        _modctr = [0]

        def _rn(node, mp, fl):
            """宣言した局所名を name$N に改名した式を返す (字句スコープ。fl に定義・純関数の有無を記録)"""
            if isinstance(node, Sym):
                n_ = mp.get(node.name)
                return Sym(n_) if n_ is not None else node
            if isinstance(node, Expr):
                h_ = node.head
                if isinstance(h_, Sym):
                    hn_ = h_.name
                    if hn_ == 'Function':
                        fl[0] = True
                    elif hn_ in ('Set', 'SetDelayed', 'UpSet', 'UpSetDelayed') and node.args and isinstance(node.args[0], Expr)                             and not (isinstance(node.args[0].head, Sym) and node.args[0].head.name in ('List', 'Part')):
                        fl[0] = True
                return Expr(_rn(h_, mp, fl), [_rn(a, mp, fl) for a in node.args])
            if isinstance(node, Rule):
                return Rule(_rn(node.lhs, mp, fl), _rn(node.rhs, mp, fl), node.delayed)
            if isinstance(node, list):
                return [_rn(a, mp, fl) for a in node]
            return node

        def _module_lex(args):
            decl = args[0]
            names = Interp._declared_names(decl)
            if not names or any(_resolve(n_) != n_ for n_ in names):
                return _module(args)
            _modctr[0] += 1
            mp = {n_: '%s$%d' % (n_, _modctr[0]) for n_ in names}
            fl = [False]
            items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
            vals = {}
            for it in items:                     # 初期値は外側の環境で評価する (Module[{x = 7, y = x}, y] の y は外側の x)
                if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')):
                    lhs, init = it.args[0], it.args[1]
                    if isinstance(lhs, Expr) and lhs.head == S('List'):
                        v_ = I.eval(init)
                        for sym, x_ in zip(lhs.args, v_ if isinstance(v_, list) else [v_]):
                            vals[mp[sym.name]] = _copy(x_)
                    elif it.head == S('SetDelayed'):
                        vals[mp[lhs.name]] = Expr(S('$Delayed'), [_rn(init, mp, fl)])
                    else:
                        vals[mp[lhs.name]] = _copy(I.eval(init))
                elif isinstance(it, Sym):
                    vals[mp[it.name]] = Sym(mp[it.name])
            fr = Frame()
            fr.update(vals)
            I.frames.append(fr)
            try:
                r_ = Null
                for body_ in args[1:]:
                    r_ = I.eval(_rn(body_, mp, fl))
            finally:
                I.frames.pop()
            if fl[0] or not isinstance(r_, (float, str, int, complex, list)) and r_ is not Null:
                for k_, v_ in fr.items():            # 定義・純関数・記号の結果が局所名を参照しうるので、値を残す
                    if k_ in fr and not (isinstance(v_, Sym) and v_.name == k_):
                        I.globals[k_] = v_
            return r_
        B['Module'] = _module_lex

        def _block(args):
            # Block[{mesh}, …]: 局所にした記号の関数定義 (mesh[x_] := …) は入るとき退避して空にし、出るとき戻す
            # (Plots.n は ListContourPlot の Block の中で、色の目盛りの ListDensityPlot が同じ Block[{nx, ny, mesh}] を使う)
            decl = args[0]
            items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
            names = []
            for it in items:
                t_ = it.args[0] if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')) and it.args else it
                if isinstance(t_, Sym):
                    names.append(_resolve(t_.name))
            saved = {}
            for nm in names:
                if nm in I.builtins:
                    continue
                saved[nm] = (I.downvalues.get(nm), I.subvalues.get(nm))
                I.downvalues.pop(nm, None)
                I.subvalues.pop(nm, None)
            try:
                return _module(args)
            finally:
                for nm, (dv_, sv_) in saved.items():
                    if dv_ is None:
                        I.downvalues.pop(nm, None)
                    else:
                        I.downvalues[nm] = dv_
                    if sv_ is None:
                        I.subvalues.pop(nm, None)
                    else:
                        I.subvalues[nm] = sv_
        B['Block'] = _block

        def _with(args):
            fr = Frame()
            fr.immut = True
            fr.is_with = True
            decl = args[0]
            items = list(decl.args) if isinstance(decl, Expr) and decl.head == S('List') else [decl]
            k_ = 0
            while k_ < len(items):                 # With[{a=1; b=2, c=3}, …]: 宣言の中の ; で区切った代入も全部束縛する (KBMainFrame.n)
                if isinstance(items[k_], Expr) and items[k_].head == S('CompoundExpression'):
                    items[k_:k_ + 1] = list(items[k_].args)
                    continue
                k_ += 1
            for it in items:
                if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')) and not isinstance(it.args[0], Sym):
                    # With[{f[x_] := …}, …]: 関数の定義を宣言に書ける (Canvas.n の TickFun)。定義はそのまま行う (K64 も同じ)
                    I.eval(it)
                    continue
                if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')):
                    nm_ = it.args[0].name
                    for fr0_ in reversed(I.frames):
                        if nm_ not in fr0_ and getattr(fr0_, 'boundary', False):
                            break                     # 別の関数 (規則) の呼び出しの外側の引数は見ない
                        if nm_ in fr0_:
                            if fr0_.immut and isinstance(fr0_[nm_], Sym):
                                nm_ = fr0_[nm_].name      # 関数の引数 z が記号 x に束縛されているとき、With[{z = …}, …] は x を束縛する (SAD は引数を本体に文字どおり代入する)
                            break
                    fr[nm_] = I.eval(it.args[1]) if it.head == S('Set') else Expr(S('$Delayed'), [it.args[1]])
                elif isinstance(it, Sym):
                    try:
                        fr[it.name] = I.eval(it)               # With[{x}, ...] は x の現在値をそのまま埋め込む (Tkinter.n の With[{TkSymbol,WInfoKeys}, …])
                    except SADError:
                        pass
            I.frames.append(fr)
            try:
                return I.eval(args[1]) if len(args) > 1 else Null
            finally:
                I.frames.pop()
        B['With'] = _with

        def _withpriority(args):
            old_ = getattr(I, 'priority', False)
            I.priority = True
            try:
                return _with(args)
            finally:
                I.priority = old_
        B['WithPriority'] = _withpriority

        def _priority(args):
            old_ = getattr(I, 'priority', False)
            I.priority = True
            try:
                return I.eval(args[0]) if args else Null
            finally:
                I.priority = old_
        B['Priority'] = _priority

        def _function(args):
            wf_ = [fr_ for fr_ in I.frames if getattr(fr_, 'is_with', False)]
            if len(args) == 1:
                body_ = args[0]
                if wf_:
                    body_ = _with_subst(body_)             # With[{f}, f&] は f の値を本体へ埋め込む (With は純関数の本体にも効く)
                return Func(None, body_)
            params = args[0]
            names = [p.name for p in (params.args if isinstance(params, Expr) else [params])]
            body_ = args[1]
            if wf_:
                body_ = _with_subst(body_, set(names))
            return Func(names, body_)
        B['Function'] = _function

        def _return(args):
            raise _Return(I.eval(args[0]) if args else Null)
        B['Return'] = _return
        B['Break'] = lambda a: (_ for _ in ()).throw(_Break())
        B['Continue'] = lambda a: (_ for _ in ()).throw(_Continue())

        def _check(args):
            try:
                return I.eval(args[0])
            except (SADError, ZeroDivisionError, ValueError, KeyError, IndexError, TypeError, _Throw):
                return I.eval(args[1]) if len(args) > 1 else Null
        B['Check'] = _check

        def _catch(args):
            try:
                return I.eval(args[0])
            except _Throw as t:
                return t.value
        B['Catch'] = _catch
        B['Throw'] = lambda a: (_ for _ in ()).throw(_Throw(a[0] if a else Null))

        def _timing(args):
            t0 = time.time()
            v = I.eval(args[0])
            return [time.time() - t0, v]
        B['Timing'] = _timing
        B['TimeUsed'] = lambda a: time.process_time()
        for nm in ('Protect', 'Unprotect', 'Hold', 'HoldForm', 'Message', 'Definition',
                   'Literal', 'Unevaluated', 'SetOptions', 'DefineOptions', 'FFS$NumericalDerivative', 'Off', 'On'):
            B[nm] = lambda a, nm=nm: (I.eval(a[0]) if nm in ('HoldForm', 'Unevaluated') and a else Null)
        B['Literal'] = lambda a: a[0] if len(a) == 1 else Null     # Literal[x] は x そのもの (評価も改名もしない記号。Canvas.n の Literal[TickFontScale]->1)
        def _with_subst(x, skip=()):
            """With[{y=…}, Hold[y]] は Hold[…の値] (With は保持された式の中にも値を埋め込む)"""
            if isinstance(x, Sym):
                if x.name in skip:
                    return x
                for fr_ in reversed(I.frames):
                    if x.name not in fr_ and getattr(fr_, 'boundary', False):
                        break                              # 別の関数の呼び出しの内側では、呼び出し元の With の変数は見ない
                    if x.name in fr_:
                        if getattr(fr_, 'is_with', False):
                            v_ = fr_[x.name]
                            if not (isinstance(v_, Expr) and v_.head == S('$Delayed')):
                                return v_
                        return x
                return x
            if isinstance(x, Expr):
                return Expr(_with_subst(x.head, skip), [_with_subst(a_, skip) for a_ in x.args])
            if isinstance(x, list):
                return [_with_subst(a_, skip) for a_ in x]
            if isinstance(x, Rule):
                return Rule(_with_subst(x.lhs, skip), _with_subst(x.rhs, skip), x.delayed)
            return x

        def _hold(a):
            a = [x_.v if isinstance(x_, _Lit) else x_ for x_ in a]
            if any(getattr(fr_, 'is_with', False) for fr_ in I.frames):
                a = [_with_subst(x_) for x_ in a]
            return Expr(S('Hold'), list(a))
        B['Hold'] = _hold                 # Hold[x] は x を評価しないまま Hold[x] を返す (WidgetCopyList = Hold[] に Union で足す)

        def _setattrs(a):
            """SetAttributes[f, attr] / SetAttributes[{f1, f2}, {attr1, attr2}] (引数は保持されたまま来る)"""
            def _names(x):
                if isinstance(x, Sym):
                    return [x.name]
                if isinstance(x, Expr) and x.head == S('List'):
                    return [n for y in x.args for n in _names(y)]
                if isinstance(x, list):
                    return [n for y in x for n in _names(y)]
                return []
            if len(a) >= 2:
                new_ = set(_names(a[1]))
                for nm in _names(a[0]):
                    cur_ = I.attributes.setdefault(nm, set())
                    if 'HoldNone' in new_:
                        cur_ -= {'HoldAll', 'HoldFirst', 'HoldRest'}      # SetAttributes[f, HoldNone] は保持の指定を外す
                    cur_.update(new_ - {'HoldNone'})
            return Null
        B['SetAttributes'] = _setattrs
        def _ruledelayed(a):
            # With[{s1, leng}, {PlotRange :> {s1 - 0.001, s1 + leng}}]: 遅延規則の右辺にも With の値を (評価せずに) 埋め込む。
            # 右辺は規則を使うときに評価されるが、そのとき With の束縛はもうない (gs は規則を作る時に記号を値に置き換える)
            rhs = I._close(a[1]) if I._has_immut() else a[1]
            return Rule(I.eval(a[0]), rhs, True)
        B['RuleDelayed'] = _ruledelayed
        B['Pattern'] = lambda a: Expr(S('Pattern'), a)
        B['Blank'] = lambda a: Expr(S('Blank'), a)
        B['BlankSequence'] = lambda a: Expr(S('BlankSequence'), a)
        B['BlankNullSequence'] = lambda a: Expr(S('BlankNullSequence'), a)
        B['Optional'] = lambda a: Expr(S('Optional'), a)
        B['PatternTest'] = lambda a: Expr(S('PatternTest'), a)
        B['Condition'] = lambda a: Expr(S('Condition'), a)
        B['Slot'] = lambda a: Expr(S('Slot'), a)
        B['SlotSequence'] = lambda a: Expr(S('SlotSequence'), a)

        # ---- 算術 (リスト対応)
        def lmap2(f, a, b, name=None):
            if isinstance(a, list) and isinstance(b, list):
                if len(a) != len(b):
                    raise SADError('list length mismatch')
                return [lmap2(f, x, y, name) for x, y in zip(a, b)]
            if isinstance(a, list):
                return [lmap2(f, x, b, name) for x in a]
            if isinstance(b, list):
                return [lmap2(f, a, y, name) for y in b]
            if is_num(a) and is_num(b):
                return f(a, b)
            return _symb(name, a, b)

        def lmap1(f, a):
            if isinstance(a, list):
                return [lmap1(f, x) for x in a]
            if is_num(a):
                return f(a)
            raise SADError(f'number expected, got {I.to_str(a)}')

        def _plus(args):
            if not args:
                return 0.0
            if len(args) == 1 and _is_pattern(args[0]):
                return Expr(S('Plus'), list(args))          # Plus[x_] は評価しない (D.n の規則 Plus[x_]:>x の左辺)
            r = args[0]
            for b in args[1:]:
                if isinstance(r, str) and isinstance(b, str):
                    r = r + b
                else:
                    r = lmap2(lambda x, y: x + y, r, b, 'Plus')
            return r
        B['Plus'] = _plus
        B['Times'] = lambda a: (Expr(S('Times'), list(a)) if len(a) == 1 and _is_pattern(a[0])
                                else _fold(a, lambda x, y: lmap2(lambda p, q: p * q, x, y, 'Times'), 1.0))
        B['Subtract'] = lambda a: lmap2(lambda x, y: x - y, a[0], a[1], 'Subtract')

        def _div(x, y):
            if y == 0:
                return INF if x > 0 else (-INF if x < 0 else NAN)
            return x / y
        B['Divide'] = lambda a: lmap2(_div, a[0], a[1], 'Divide')
        B['Minus'] = lambda a: _symneg(a[0]) if isinstance(a[0], (Sym, Expr, str)) else lmap1(lambda x: -x, a[0])

        def _pow(x, y):
            try:
                if isinstance(x, complex) or isinstance(y, complex):
                    return x**y
                if x < 0 and y != int(y):
                    return cmath.exp(y * cmath.log(x))      # K64: (-8.)^(1/3) = 1+1.732 I (主値)
                if x == 0 and y < 0:
                    return INF
                return float(x)**y
            except OverflowError:
                return INF
        B['Power'] = lambda a: lmap2(_pow, a[0], a[1], 'Power')

        def _is_rat(x):
            return isinstance(x, Expr) and x.head == S('Rational')

        def _wrap_rational(name, fn):
            def w(a):
                # rational.n を読んだあと: a-b は a+(-b)、a/b は a*(1/b) として Rational の規則 (UpValue) へ回す (gs は Subtract・Divide を Plus・Times に直す)
                if len(a) == 2 and (_is_rat(a[0]) or _is_rat(a[1])) and I.downvalues.get('Rational'):
                    if name == 'Subtract':
                        return I.eval(Expr(S('Plus'), [_Lit(a[0]), Expr(S('Minus'), [_Lit(a[1])])]))
                    if not (name == 'Divide' and a[0] == 1.0):
                        return I.eval(Expr(S('Times'), [_Lit(a[0]), Expr(S('Divide'), [1.0, _Lit(a[1])])]))
                return fn(a)
            return w
        B['Subtract'] = _wrap_rational('Subtract', B['Subtract'])
        B['Divide'] = _wrap_rational('Divide', B['Divide'])

        def _fold(args, f, unit):
            r = unit
            for a in args:
                r = f(r, a)
            return r

        def _cmp(f):
            def g(args):
                if len(args) > 2:
                    # a<x<b などの連鎖: 隣どうしの比較がすべて成り立つとき真 (K64 も同じ)
                    rs = [g([args[k], args[k + 1]]) for k in range(len(args) - 1)]
                    if all(isinstance(r, float) for r in rs):
                        return 1.0 if all(r == 1.0 for r in rs) else 0.0
                    return Expr(S('$Cmp'), args)
                a, b = args[0], args[1]
                if isinstance(a, str) and isinstance(b, str):
                    return 1.0 if f(a, b) else 0.0
                if is_num(a) and is_num(b):
                    return 1.0 if f(a, b) else 0.0
                if isinstance(a, list) or isinstance(b, list):
                    try:
                        return lmap2(lambda x, y: 1.0 if f(x, y) else 0.0, a, b)       # gs: 比較はリストの要素ごと ({1,2}>=2 は {0,1})
                    except Exception:
                        return Expr(S('$Cmp'), args)
                return Expr(S('$Cmp'), args)
            return g
        def _equal2(a):
            return (1.0 if I.same(a[0], a[1]) else 0.0) if not (is_num(a[0]) and is_num(a[1])) else (1.0 if a[0] == a[1] else 0.0)
        def _eqthread(x, y, neg):
            # gs: リストどうしの == は要素ごと ({1,2}=={1,3} は {1,0})、入れ子も再帰。長さが違えば 0。リストとスカラーは未評価
            if isinstance(x, list) and isinstance(y, list):
                if len(x) != len(y):
                    return 1.0 if neg else 0.0
                return [_eqthread(p, q, neg) for p, q in zip(x, y)]
            if isinstance(x, list) or isinstance(y, list):
                return Expr(S('Unequal' if neg else 'Equal'), [x, y])
            if I.same(x, y):
                return 0.0 if neg else 1.0
            if isinstance(x, (Sym, Expr)) or isinstance(y, (Sym, Expr)):
                return Expr(S('Unequal' if neg else 'Equal'), [x, y])
            return 1.0 if neg else 0.0

        def _equal(a):
            if len(a) == 2 and (isinstance(a[0], list) or isinstance(a[1], list)):
                return _eqthread(a[0], a[1], False)
            return _equal2(a) if len(a) == 2 else (1.0 if all(_equal2([a[k], a[k + 1]]) == 1.0 for k in range(len(a) - 1)) else 0.0)
        B['Equal'] = _equal
        B['Unequal'] = lambda a: _eqthread(a[0], a[1], True) if len(a) == 2 and (isinstance(a[0], list) or isinstance(a[1], list)) else 1.0 - _equal(a)
        B['Greater'] = _cmp(lambda x, y: x > y)
        B['GreaterEqual'] = _cmp(lambda x, y: x >= y)
        B['Less'] = _cmp(lambda x, y: x < y)
        B['LessEqual'] = _cmp(lambda x, y: x <= y)
        B['SameQ'] = lambda a: 1.0 if all(I.same(a[k], a[k + 1]) for k in range(len(a) - 1)) else 0.0
        # a<=>b<=>c は a<=>b かつ b<=>c (隣どうしが違う)
        B['UnsameQ'] = lambda a: 1.0 if all(not I.same(a[k], a[k + 1]) for k in range(len(a) - 1)) else 0.0
        def _and(args):
            for x in args:                              # 左から評価し、偽が出たらそこで止める (a<=>{}&&a[[1,2]])
                if truthy(I.eval(x)) is not True:
                    return 0.0
            return 1.0

        def _or(args):
            for x in args:
                if truthy(I.eval(x)) is True:
                    return 1.0
            return 0.0
        B['And'] = _and
        B['Or'] = _or
        B['Not'] = lambda a: (0.0 if truthy(a[0]) else 1.0) if truthy(a[0]) is not None else Expr(S('Not'), a)
        B['Xor'] = lambda a: 1.0 if (truthy(a[0]) is True) != (truthy(a[1]) is True) else 0.0

        # ---- 数学
        def mathf(f, fc=None):
            # 複素数には cmath 版 (fc) を使う (Sqrt が負で複素数を返したあとも計算を続けられるように)
            if fc is None:
                return lambda a: lmap1(f, a[0])
            return lambda a: lmap1(lambda x: fc(x) if isinstance(x, complex) else f(x), a[0])

        def _sqrt(x):
            if isinstance(x, complex):
                return x**0.5
            if x != x:
                return x
            return math.sqrt(x) if x >= 0 else complex(0.0, math.sqrt(-x))
        B['Sqrt'] = mathf(_sqrt)
        def _exp(x):
            try:
                return math.exp(x)
            except OverflowError:
                return INF
        B['Exp'] = mathf(_exp, cmath.exp)
        def _log1(x):
            if isinstance(x, complex):
                return cmath.log(x) if x != 0 else -INF
            if x > 0:
                return math.log(x)
            if x == 0:
                return -INF
            if x != x:
                return NAN
            return complex(math.log(-x), math.pi)       # K64: Log[-2.] = .693 + 3.14159 I
        B['Log'] = lambda a: lmap1(_log1, a[0]) if len(a) == 1 else lmap2(lambda b, x: _log1(x) / _log1(b), a[0], a[1])
        B['Sin'] = mathf(math.sin, cmath.sin)
        B['Cos'] = mathf(math.cos, cmath.cos)
        B['Tan'] = mathf(math.tan, cmath.tan)
        B['ArcSin'] = mathf(lambda x: math.asin(x) if -1 <= x <= 1 else (NAN if x != x else cmath.asin(complex(x, 0.0))), cmath.asin)
        B['ArcCos'] = mathf(lambda x: math.acos(x) if -1 <= x <= 1 else (NAN if x != x else cmath.acos(complex(x, 0.0))), cmath.acos)
        B['ArcTan'] = lambda a: lmap1(lambda x: cmath.atan(x) if isinstance(x, complex) else math.atan(x), a[0]) if len(a) == 1 else lmap2(lambda x, y: math.atan2(y, x), a[0], a[1])
        B['Sinh'] = mathf(math.sinh, cmath.sinh)
        B['Cosh'] = mathf(math.cosh, cmath.cosh)
        B['Tanh'] = mathf(math.tanh, cmath.tanh)
        B['Abs'] = mathf(abs)
        B['Sign'] = mathf(lambda x: (x > 0) - (x < 0) + 0.0)
        B['Floor'] = mathf(lambda x: float(math.floor(x)))
        B['Ceiling'] = mathf(lambda x: float(math.ceil(x)))
        B['Round'] = mathf(lambda x: float(math.floor(x + 0.5)) if x >= 0 else -float(math.floor(-x + 0.5)))
        B['IntegerPart'] = mathf(lambda x: float(int(x)))
        B['FractionalPart'] = mathf(lambda x: x - float(int(x)))
        B['Mod'] = lambda a: lmap2(lambda x, y: x - y * math.floor(x / y), a[0], a[1])
        B['Re'] = mathf(lambda x: x.real if isinstance(x, complex) else x)
        B['Im'] = mathf(lambda x: x.imag if isinstance(x, complex) else 0.0)
        B['Conjugate'] = mathf(lambda x: x.conjugate() if isinstance(x, complex) else x)
        B['Arg'] = mathf(lambda x: math.atan2(x.imag, x.real) if isinstance(x, complex) else (0.0 if x >= 0 else math.pi))
        B['N'] = lambda a: a[0]
        def _factorial(x):
            if x == int(x):
                if x < 0:
                    return NAN
                return float(math.factorial(int(x))) if x < 171 else INF
            try:
                return math.gamma(x + 1.0)
            except OverflowError:
                return INF
        B['Factorial'] = mathf(_factorial)
        B['Gamma'] = mathf(math.gamma)
        def _cerf(f, cf):
            return mathf(f, lambda z: (lambda v: v.real if v.imag == 0 else v)(complex(cf(z))))
        from scipy import special as _sps
        B['Erf'] = _cerf(math.erf, _sps.erf)
        B['Erfc'] = _cerf(math.erfc, _sps.erfc)

        def _gc4(args):
            if len(args) != 4 or not all(is_num(v) for v in args):
                raise SADError('General::wrongtype: "x, y, sigx, sigy"')
            return [float(v) for v in args]

        def _gaussian_coulomb(a):
            from .gausscoul import bbkick1
            z = bbkick1(*_gc4(a))
            return z.real if z.imag == 0.0 else z

        def _gaussian_coulomb_u(a):
            from .gausscoul import twspu
            x, y, sx, sy = _gc4(a)
            return twspu(x, y, sx, sy, 1e-10, 1e-13)

        def _bbbrem1(a):
            from .bbbrem import bbbrem1
            if len(a) != 2:
                raise SADError('General::narg: "2"')
            if not is_num(a[0]):
                raise SADError('General::wrongtype: "Sqrt[s] (eV) for #1"')
            if not is_num(a[1]):
                raise SADError('General::wrongtype: "Min. ds/s for #2"')
            if a[1] <= 0.0 or a[1] >= 1.0:
                raise SADError('General::wrongval: "0 < #2 <1", "ds/s"')
            p, sig = bbbrem1(float(a[0]), float(a[1]), I.rng.uniform0)
            return [list(p), sig]

        B['BBBrem1'] = _bbbrem1
        B['GaussianCoulomb'] = _gaussian_coulomb
        B['GaussianCoulombU'] = _gaussian_coulomb_u

        def _gaussian_coulomb_fitted(a):
            from .gausscoul import twspfu
            x, y, sx, sy = _gc4(a)
            return [float(v) for v in twspfu(x, y, sx, sy)]
        B['GaussianCoulombFitted'] = _gaussian_coulomb_fitted

        def _flatnums(args):
            out = []
            for a in args:
                if isinstance(a, list):
                    out.extend(_flatnums(a))
                else:
                    out.append(a)
            return out
        B['Max'] = lambda a: max(_flatnums(a)) if _flatnums(a) else -INF
        B['Min'] = lambda a: min(_flatnums(a)) if _flatnums(a) else INF
        B['Random'] = lambda a: I.rng.uniform0() if I.rng is not None else float(np.random.random())
        B['GaussRandom'] = lambda a: (I.rng.gauss() if I.rng is not None else float(np.random.normal())) if not a else [I.rng.gauss() for _ in range(int(a[0]))]
        B['SeedRandom'] = lambda a: (I.rng.seed_random(a[0]) if I.rng is not None else np.random.seed(int(a[0]))) or Null

        # ---- リスト
        B['Length'] = lambda a: float(len(a[0])) if isinstance(a[0], (list, str)) else 0.0

        def _range(a):
            if len(a) == 1:
                return [float(k) for k in range(1, int(math.floor(num(a[0]) + 1e-9)) + 1)]
            lo, hi = num(a[0]), num(a[1])
            d = num(a[2]) if len(a) > 2 else 1.0
            n = int(math.floor((hi - lo) / d + 1e-9)) + 1
            return [lo + k * d for k in range(max(n, 0))]
        B['Range'] = _range
        def _part_builtin(a):
            r_ = I.part(a[0], a[1:])
            if isinstance(a[0], Expr) and a[0].head == S('Hold') and isinstance(r_, (list, Expr)):
                r_ = I.eval(Expr(S('List'), list(r_)) if isinstance(r_, list) else r_)            # Hold[…][[1]] は取り出した式を評価し直す (Couple.n: EVCompiled[[1]])
            return r_
        B['Part'] = _part_builtin
        B['First'] = lambda a: a[0][0]
        B['Last'] = lambda a: a[0][-1]
        B['Rest'] = lambda a: a[0][1:]
        B['Most'] = lambda a: a[0][:-1]

        def _take(a):
            lst, n = a[0], a[1]
            if isinstance(n, list):
                i, j = int(n[0]), int(n[1])
                i = len(lst) + i + 1 if i < 0 else i
                j = len(lst) + j + 1 if j < 0 else j
                return lst[i - 1:j]
            n = int(n)
            return lst[:n] if n >= 0 else lst[n:]
        B['Take'] = _take
        B['Drop'] = lambda a: (a[0][int(a[1]):] if int(a[1]) >= 0 else a[0][:int(a[1])]) if not isinstance(a[1], list) else a[0][:int(a[1][0]) - 1] + a[0][int(a[1][1]):]
        B['Append'] = lambda a: list(a[0]) + [a[1]]
        B['Prepend'] = lambda a: [a[1]] + list(a[0])
        B['Join'] = lambda a: _fold(a, lambda x, y: (x + y) if isinstance(x, list) else y, [])

        def _flatten(a):
            out = []
            depth = int(a[1]) if len(a) > 1 else 10**6

            def rec(x, d):
                for y in x:
                    if isinstance(y, list) and d > 0:
                        rec(y, d - 1)
                    else:
                        out.append(y)
            rec(a[0], depth)
            return out
        B['Flatten'] = _flatten
        B['Reverse'] = lambda a: list(reversed(a[0]))

        def _sort(a):
            lst = list(a[0])
            if len(a) > 1:
                import functools
                f = a[1]

                def cmp(x, y):
                    r = I.apply_value(f, [x, y])
                    return -1 if truthy(r) else (1 if truthy(I.apply_value(f, [y, x])) else 0)
                return sorted(lst, key=functools.cmp_to_key(cmp))

            def key(x):
                if is_num(x):
                    return (0, x, '')
                if isinstance(x, str):
                    return (1, 0, x)
                if isinstance(x, list):
                    return (2, 0, tuple(key(y) for y in x))
                return (3, 0, repr(x))
            return sorted(lst, key=key)
        B['Sort'] = _sort
        B['Union'] = lambda a: _sort([list({_hkey(x): x for x in _flatten([a if len(a) > 1 else a[0], 1])}.values())])
        B['Transpose'] = lambda a: [list(r) for r in zip(*a[0])]
        B['Dimensions'] = lambda a: _dims(a[0])
        B['Total'] = lambda a: _fold(a[0], lambda x, y: I.builtins['Plus']([x, y]), 0.0)
        B['Apply'] = lambda a: I.apply_value(a[0], list(a[1]) if isinstance(a[1], list) else [a[1]], a[0])
        B['Map'] = lambda a: [I.apply_value(a[0], [x], a[0]) for x in a[1]] if isinstance(a[1], list) else a[1]
        B['MapAll'] = B['Map']

        def _scan(a):
            """Scan[f, list]: 各要素に f を適用 (結果は捨てる)。戻りは Null"""
            if isinstance(a[1], list):
                for x in a[1]:
                    I.apply_value(a[0], [x], a[0])
            return Null
        B['Scan'] = _scan
        B['MapIndexed'] = lambda a: [I.apply_value(a[0], [x, [float(k + 1)]], a[0]) for k, x in enumerate(a[1])]
        B['Select'] = lambda a: [x for x in a[0] if truthy(I.apply_value(a[1], [x], a[1])) is True]
        B['Position'] = lambda a: [[float(k + 1)] for k, x in enumerate(a[0]) if I.same(x, a[1])]
        B['MemberQ'] = lambda a: 1.0 if any(I.same(x, a[1]) for x in a[0]) else 0.0
        B['FreeQ'] = lambda a: 0.0 if any(I.same(x, a[1]) for x in a[0]) else 1.0
        B['Count'] = lambda a: float(sum(1 for x in a[0] if I.same(x, a[1])))
        B['Delete'] = lambda a: [x for k, x in enumerate(a[0]) if k + 1 != int(a[1])]
        B['Insert'] = lambda a: a[0][:int(a[2]) - 1] + [a[1]] + a[0][int(a[2]) - 1:]
        B['ReplacePart'] = lambda a: _replace_part(a[0], a[1], a[2])
        B['Partition'] = lambda a: [a[0][k:k + int(a[1])] for k in range(0, len(a[0]) - int(a[1]) + 1, int(a[2]) if len(a) > 2 else int(a[1]))]
        B['Thread'] = lambda a: _thread(a[0])
        B['Table$'] = B['Table']
        B['Accumulate'] = lambda a: list(np.cumsum(a[0]))
        B['Fold'] = lambda a: _fold(a[2], lambda x, y: I.apply_value(a[0], [x, y], a[0]), a[1])
        B['Nest'] = lambda a: _fold(range(int(a[2])), lambda x, _: I.apply_value(a[0], [x], a[0]), a[1])
        B['Identity'] = lambda a: a[0]
        B['List'] = lambda a: list(a)
        B['Sequence'] = lambda a: Seq(a)
        B['Null'] = lambda a: Seq(list(a))          # Null[] は引数列に何も足さない (リストや引数の最後に置く SAD の習慣)

        def _thread(x):
            if isinstance(x, list) and x and all(isinstance(y, list) for y in x) and len({len(y) for y in x}) == 1:
                return [list(r) for r in zip(*x)]
            if isinstance(x, Rule):
                if isinstance(x.lhs, list) and isinstance(x.rhs, list):
                    return [Rule(l, r) for l, r in zip(x.lhs, x.rhs)]
                if isinstance(x.lhs, list):
                    return [Rule(l, x.rhs) for l in x.lhs]
            if isinstance(x, Expr):
                lists = [a for a in x.args if isinstance(a, list)]
                if not lists:
                    return x
                n = len(lists[0])
                return [I.eval(Expr(x.head, [_Lit(a[k] if isinstance(a, list) else a) for a in x.args])) for k in range(n)]
            return x

        # ---- 行列
        def _mat(a):
            return np.array(a, dtype=complex if _has_complex(a) else float)

        def _tolist(m):
            if np.iscomplexobj(m):
                m = m.real if np.all(np.abs(m.imag) == 0) else m
            return m.tolist() if m.ndim else float(m)
        B['Dot'] = lambda a: _tolist(np.dot(_mat(a[0]), _mat(a[1])))
        B['Inverse'] = lambda a: _tolist(np.linalg.inv(_mat(a[0])))
        B['Det'] = lambda a: float(np.linalg.det(_mat(a[0])).real)
        B['LinearSolve'] = lambda a: _tolist(np.linalg.lstsq(_mat(a[0]), _mat(a[1]), rcond=None)[0])
        def _diagmat(a):
            v = a[0]
            cx = any(isinstance(x, complex) for x in v) if isinstance(v, list) else False
            m = np.diag(np.array(v, dtype=complex if cx else float)).tolist()
            return [[(x.real if isinstance(x, complex) and x.imag == 0.0 else x) for x in r] for r in m]
        B['DiagonalMatrix'] = _diagmat
        B['IdentityMatrix'] = lambda a: np.eye(int(a[0])).tolist()
        B['Outer'] = lambda a: [[I.apply_value(a[0], [x, y], a[0]) for y in a[2]] for x in a[1]]

        def _eigensystem(a):
            m = _mat(a[0])
            if m.ndim != 2 or m.shape[0] != m.shape[1]:
                raise SADError('General::wrongtype: Argument must be Numerical Square Matrix in Eigensystem')
            if not np.iscomplexobj(m) and m.shape[0] > 8 and np.array_equal(m, m.T):
                # 大きな対称行列 (CodCorr の応答行列 tres.res など): teigen の直訳は縮退した固有値で非直交のベクトルを返し、SAD (-Ofast) とも
                # 丸めで一致しない (160x160 で CodCorr の残差が SAD と 4 倍違った)。対称行列の厳密な固有値分解 (単位直交ベクトル) を使う。
                # 小さな行列 (n<=8) は teigen の直訳で SAD と並び・規格化が一致する
                w_, v_ = np.linalg.eigh(m)
                order = np.argsort(-w_)
                return [[float(w_[k]) for k in order], [[float(x) for x in v_[:, k]] for k in order]]
            if not np.iscomplexobj(m):
                # tfdot.f tfreigen: 実行列は teigen.f (K. Oide) で解く。固有ベクトルは単位ベクトルではない (SAD 固有の規格化)。
                # 行 k = W の列 k。複素固有値の対 (k, k+1) は 行 k = 実部 + i 虚部、行 k+1 = 実部 - i 虚部 (W の列 k, k+1 が実部, 虚部)
                # 固有値の並び・規格化は teigen の直訳と同じ (numpy.linalg.eig は並びも規格化も違い、縮退で虚部の雑音も出る)
                from .teigen import teigen
                ceig, W = teigen(m)
                n = len(ceig)
                vals = [complex(x) if x.imag != 0.0 else float(x.real) for x in ceig]
                vecs = []
                k = 0
                while k < n:
                    if ceig[k].imag != 0.0 and k + 1 < n:
                        vr, vi = W[:, k], W[:, k + 1]
                        vecs.append([complex(x, y) if y != 0.0 else float(x) for x, y in zip(vr, vi)])
                        vecs.append([complex(x, -y) if y != 0.0 else float(x) for x, y in zip(vr, vi)])
                        k += 2
                    else:
                        vecs.append([float(x) for x in W[:, k]])
                        k += 1
                return [vals, vecs]
            w, v = np.linalg.eig(m)
            vals = [complex(x) if abs(x.imag) > 0 else float(x.real) for x in w]
            vecs = [[complex(x) if abs(x.imag) > 0 else float(x.real) for x in v[:, k]] for k in range(len(w))]
            return [vals, vecs]
        B['Eigensystem'] = _eigensystem
        B['Eigenvalues'] = lambda a: _eigensystem(a)[0]
        B['SingularValues'] = lambda a: _tolist(np.linalg.svd(_mat(a[0]), compute_uv=False))

        # ---- 文字列
        B['StringJoin'] = lambda a: ''.join(I.to_str(x) for x in a)         # gs: リストは文字列に連結せず {"a","b"} の形の文字列 (StringJoin[{"a","b"}] = {"a","b"})
        B['ToString'] = lambda a: I.to_str(a[0], a[1] if len(a) > 1 and isinstance(a[1], str) else None)
        B['ToExpression'] = lambda a: I.run(a[0], 'script') if isinstance(a[0], str) else a[0]
        B['StringLength'] = lambda a: float(len(a[0]))
        def _standard_form(a):
            """StandardForm[expr]: expr の評価の間だけ $FORM を標準 (既定) の書式にする。台本が $FORM = "6.2" としていても、
            Tcl へ送る命令の中の数 ("after " // 1001) が 1001. にならない (init.n の StandardForm)"""
            if not a:
                return Null
            had = I.has('$FORM')
            old = I.get('$FORM')
            I.assign('$FORM', '')
            try:
                return I.eval(a[0])
            finally:
                if had:
                    I.assign('$FORM', old)
                else:
                    I.globals.pop('$FORM', None)
        B['StandardForm'] = _standard_form
        def _rgbcolor(a):
            """RGBColor[r, g, b] = "#RRGGBB" (成分は 0..1、nint(256 x) を 0..255 に切る)。RGBColor[{{r,g,b}, …}] は文字列のリスト (tfloor.f の tfrgbcolor)"""
            def one(r, g, b):
                def q(x):
                    return min(255, max(0, int(round(256.0 * x))))
                return '#%02X%02X%02X' % (q(r), q(g), q(b))
            if len(a) == 1 and isinstance(a[0], list):
                if all(isinstance(t, list) and len(t) == 3 and all(is_num(x) for x in t) for t in a[0]):
                    return [one(*t) for t in a[0]]
            elif len(a) == 3 and all(is_num(x) for x in a):
                return one(*a)
            raise SADError('General::wrongtype: Argument must be "r,g,b or {{r1, g1, b1}, ...}" in RGBColor')
        B['RGBColor'] = _rgbcolor
        B['ToUpperCase'] = lambda a: a[0].upper() if isinstance(a[0], str) else a[0]
        B['ToLowerCase'] = lambda a: a[0].lower() if isinstance(a[0], str) else a[0]      # 文字列以外はそのまま返す
        B['Characters'] = lambda a: list(a[0])
        B['StringTake'] = lambda a: _take([list(a[0]), a[1]]) and ''.join(_take([list(a[0]), a[1]]))
        B['StringDrop'] = lambda a: ''.join(B['Drop']([list(a[0]), a[1]]))
        B['StringReplace'] = lambda a: _string_replace(a[0], a[1])
        def _stringposition(a):
            """StringPosition[s, pat] / [s, pat, n]: pat は文字列・文字列のリスト・"a"|"b"。重なる一致も数える (K64 と同じ)。
            選択肢は与えた順に、それぞれの全位置を並べる (同じ選択肢は 1 回だけ)。n があれば先頭の n 個"""
            sv, pat = a[0], a[1]
            if isinstance(pat, str):
                alts = [pat]
            elif isinstance(pat, list):
                alts = list(pat)
            elif isinstance(pat, Expr) and pat.head == S('Alternatives'):
                alts = list(pat.args)
            else:
                raise SADError('General::wrongtype: Argument must be String, List of Strings or Alternatives for #2 in StringPosition')
            seen, out = set(), []
            for q in alts:
                if not isinstance(q, str) or q in seen:
                    continue
                seen.add(q)
                if q == '':
                    continue
                i = sv.find(q)
                while i >= 0:
                    out.append([float(i + 1), float(i + len(q))])
                    i = sv.find(q, i + 1)
            if len(a) > 2 and is_num(a[2]):
                out = out[:int(a[2])]
            return out
        B['StringPosition'] = _stringposition
        def _smq(a):
            if len(a) != 2 or not isinstance(a[0], str) or not isinstance(a[1], str):
                raise SADError('General::wrongtype:  Argument must be [String, Pattern-string] in StringMatchQ[%s]' % ','.join(I.to_str(x, None, True) for x in a))
            if a[0] == '':
                return 0.0                        # gs: 空文字列は何にも合わない
            return 1.0 if _sad_wild_rx(a[1]).fullmatch(a[0]) else 0.0
        B['StringMatchQ'] = _smq
        B['StringForm'] = lambda a: _string_form(a)
        B['StringTrim'] = lambda a: a[0].strip()
        B['StringQ'] = lambda a: 1.0 if isinstance(a[0], str) else 0.0
        B['NumberQ'] = lambda a: 1.0 if is_num(a[0]) else 0.0
        B['RealQ'] = lambda a: 1.0 if (is_num(a[0]) and not isinstance(a[0], complex)) else 0.0
        B['ComplexQ'] = lambda a: 1.0 if isinstance(a[0], complex) else 0.0
        B['ListQ'] = lambda a: 1.0 if isinstance(a[0], list) else 0.0
        def _vectorq(a):
            if not (isinstance(a[0], list) and all(not isinstance(x, list) for x in a[0])):
                return 0.0
            if len(a) > 1:
                return 1.0 if all(truthy(I.apply_value(a[1], [x], a[1])) is True for x in a[0]) else 0.0
            return 1.0
        B['VectorQ'] = _vectorq
        def _matrixq(a):
            if not (isinstance(a[0], list) and a[0] and all(isinstance(x, list) for x in a[0])):
                return 0.0
            if len(a) > 1:
                return 1.0 if all(truthy(I.apply_value(a[1], [y], a[1])) is True for x in a[0] for y in x) else 0.0
            return 1.0
        B['MatrixQ'] = _matrixq
        B['SymbolQ'] = lambda a: 1.0 if isinstance(a[0], Sym) else 0.0
        def _head(a):
            x = a[0]
            if isinstance(x, complex):
                return S('Complex')
            if is_num(x):
                return S('Real')
            if isinstance(x, str):
                return S('String')
            if isinstance(x, list):
                return S('List')
            if x is Null or isinstance(x, Sym):
                return S('Symbol')
            if isinstance(x, Rule):
                return S('RuleDelayed') if getattr(x, 'delayed', False) else S('Rule')
            if isinstance(x, Func):
                return S('Function')
            if isinstance(x, Expr):
                return S('Times') if (isinstance(x.head, Sym) and x.head.name == 'Minus') else x.head
            return S('Symbol')
        B['Head'] = _head

        def _string_replace(s, rules):
            if isinstance(rules, Rule):
                rules = [rules]
            for r in rules:
                s = s.replace(r.lhs, r.rhs)
            return s

        def _string_form(a):
            s = a[0]
            for k, v in enumerate(a[1:]):
                s = s.replace('`' + str(k + 1) + '`', I.to_str(v))
            return s.replace('``', '')

        # ---- 規則
        B['Rule'] = lambda a: Rule(a[0], a[1])
        B['ReplaceAll'] = lambda a: _replace_all(a[0], a[1])
        def _replace_repeated(a):
            x = a[0]
            for _ in range(10000):
                y = _replace_all(x, a[1])
                if os.environ.get('PYSAD_RRTRACE'):
                    sys.stderr.write('[//.] ' + I.to_str(x, quote=True)[:150] + '  =>  ' + I.to_str(y, quote=True)[:150] + chr(10))
                if I.same(x, y) or repr(x) == repr(y):
                    return y
                x = y
            return x
        B['ReplaceRepeated'] = _replace_repeated
        B['Override'] = lambda a: _override(a)
        B['Options'] = lambda a: []

        def _override(a):
            rules = []
            for x in a:
                if isinstance(x, Rule):
                    rules.append(x)
                elif isinstance(x, list):
                    rules.extend(y for y in _flatten([x]) if isinstance(y, Rule))
            seen = {}
            out = []
            for r in rules:
                k = _hkey(r.lhs)
                if k not in seen:
                    seen[k] = True
                    out.append(r)
            if len(out) > 1:
                # gs の Override は結果を左辺の標準順 (文字列、続いて記号を大小無視、同じなら小文字が先) に並べる
                try:
                    ks = I.builtins['Sort']([[ (r.lhs.v if isinstance(r.lhs, _Lit) else r.lhs) for r in out]])
                    ordered = []
                    pool = list(out)
                    for k in ks:
                        for j_, r in enumerate(pool):
                            if I.same(r.lhs.v if isinstance(r.lhs, _Lit) else r.lhs, k):
                                ordered.append(pool.pop(j_))
                                break
                    if len(ordered) == len(out):
                        out = ordered
                except Exception:
                    pass
            return out

        def _replace_all(x, rules):
            if isinstance(rules, Expr) and rules.head == S('Unevaluated$') and len(rules.args) == 1:
                rules = I.eval(rules.args[0])             # x//.Unevaluated$[ReductionRule] (D.n)
            if isinstance(rules, Rule):
                rules = [rules]
            if not isinstance(rules, list):
                return x
            flat = [r for r in _flatten([rules]) if isinstance(r, Rule)]
            flat = [Rule(r.lhs.v if isinstance(r.lhs, _Lit) else r.lhs, r.rhs.v if isinstance(r.rhs, _Lit) else r.rhs, r.delayed) for r in flat]

            def rep(v, ev=True):
                if isinstance(v, _Lit):
                    v = v.v                      # 保持引数 (HoldAll の Image[w, d] など) の中の _Lit は値に戻して照合する
                for r in flat:
                    if I.same(v, r.lhs) or (isinstance(r.lhs, Sym) and isinstance(v, Sym) and r.lhs.name == v.name):
                        return I.eval(r.rhs) if (r.delayed and ev) else r.rhs
                    if _has_pat(r.lhs):
                        b_ = {}
                        if I._match_one(r.lhs, v, b_):          # x_ -> y, {x_->_, y_} :> ... のようなパターンの規則
                            if r.delayed:
                                fr_ = Frame(b_)
                                fr_.immut = True
                                I.frames.append(fr_)
                                try:
                                    return I.eval(r.rhs) if ev else _unlit(I.subst_bindings(r.rhs, b_))
                                finally:
                                    I.frames.pop()
                            return I.eval(I.subst_bindings(r.rhs, b_))
                if isinstance(v, Func):
                    return Func(v.params, _as_expr(rep(v.body, False)), v.env)        # (…&) /. Hold->Null: 純関数の本体にも規則を適用する
                if isinstance(v, list):
                    return _splice([rep(y, ev) for y in v])
                if isinstance(v, Rule):
                    return Rule(rep(v.lhs, ev), rep(v.rhs, ev and not v.delayed), v.delayed)       # a:>b の b は評価しない (Command:>KBF$AboutBox[])
                if isinstance(v, Expr):
                    if not ev or (v.head == S('Hold') and not any(isinstance(r_.lhs, Sym) and r_.lhs.name == 'Hold' for r_ in flat)):          # Hold[…] の中は規則の右辺 (:>) も評価せず埋め込む (Hold を別の頭部に替えるときは評価する)
                        h_ = rep(v.head, False)
                        a_ = _splice([rep(y, False) for y in v.args])
                        if isinstance(h_, _Null) or h_ == S('Null'):
                            return Seq(a_)                 # Hold->Null: Null[x, …] は並びとして外側のリストに展開される
                        return Expr(h_, a_)
                    return I.eval(Expr(rep(v.head), [_Lit(y) for y in _splice([rep(y) for y in v.args])]))
                return v

            def _splice(items):
                """置き換えの結果が Null[] (空の Sequence) の要素は取り除く (`{a, INF, b} /. INF :> Null[]` = {a, b})"""
                if not any(isinstance(y, Seq) for y in items):
                    return items
                out = []
                for y in items:
                    if isinstance(y, Seq):
                        out.extend(y.items)
                    else:
                        out.append(y)
                return out
            return rep(x)

        # ---- 入出力
        def _print(a):
            s = ''.join(I.to_str(x) for x in a)
            I.write_unit(6, I.wrap(s) + '\n')
            return Null
        B['Print'] = _print

        def _write(a):
            if len(a) < 2:
                raise SADError('General::narg:  Number of arguments is expected to be 2 or more in Write[%s]' % ','.join(I.to_str(x) for x in a))
            unit = a[0]
            s = ''.join(I.to_str(x) for x in a[1:])
            I.write_unit(unit if is_num(unit) else 6, I.wrap(s) + '\n')
            return Null
        B['Write'] = _write

        def _writestring(a):
            unit = a[0]
            s = ''.join(I.to_str(x) for x in a[1:])
            I.write_unit(unit if is_num(unit) else 6, s)
            return Null
        B['WriteString'] = _writestring

        def _open(a, mode):
            path = a[0]
            if mode == 'r' and isinstance(path, str) and path.startswith('!'):
                # OpenRead["!コマンド"]: コマンドの出力を読む (CSR.n の FreeMem の `top ... | grep Mem:`)。実行できなければ出力は空
                import io
                import subprocess
                try:
                    out_ = subprocess.run(path[1:], shell=True, capture_output=True, text=True, timeout=30).stdout
                except Exception:
                    out_ = ''
                if not out_.strip() and 'top ' in path and 'Mem:' in path:
                    # top のない環境 (Windows): 空き記憶領域を top の旧形式 (Linux: `Mem: …k total, …k used, …k free, …k buffers`) で作る
                    try:
                        import psutil
                        vm_ = psutil.virtual_memory()
                        tot_, free_, buf_ = vm_.total // 1024, vm_.available // 1024, 0
                    except Exception:
                        tot_, free_, buf_ = 16 * 1024 * 1024, 8 * 1024 * 1024, 0
                    out_ = 'Mem: %dk total, %dk used, %dk free, %dk buffers' % (tot_, tot_ - free_, free_, buf_) + chr(10)
                u = I.next_unit
                I.next_unit += 1
                I.files[u] = io.StringIO(out_)
                return float(u)
            if not os.path.isabs(path):
                path = os.path.join(I.cwd, path)
            u = I.next_unit
            try:
                f_ = open(path, mode, encoding='utf-8', errors='replace')
            except OSError:
                if mode == 'r':
                    # gs: OpenRead が開けないとメッセージを出して、その入力行の評価を打ち切る (h=OpenRead[..] の h は未代入のまま)。OpenWrite・OpenAppend は $Failed を返す
                    raise SADError('General::fileopen:  Open error for file %s in OpenRead[%s]' % (a[0], I.to_input_string(a[0]) if hasattr(I, 'to_input_string') else chr(34) + str(a[0]) + chr(34)))
                I.out.write('???General::fileopen:  Open error for file ' + str(a[0]) + ':' + chr(10))
                return S('$Failed')
            I.next_unit += 1
            I.files[u] = f_
            return float(u)
        B['OpenWrite'] = lambda a: _open(a, 'w')
        B['OpenAppend'] = lambda a: _open(a, 'a')
        B['OpenRead'] = lambda a: _open(a, 'r')

        def _close(a):
            if a and isinstance(a[0], list) and all(is_num(x) for x in a[0]):
                for x_ in a[0]:                       # Close[{f, g}]: 複数のファイルを閉じる
                    _close([x_])
                return Null
            if not a or not is_num(a[0]):
                raise SADError('General::wrongtype:  Argument must be Real number in Close$[%s]' % ','.join(I.to_str(x) for x in a))
            f = I.files.pop(int(a[0]), None) if int(a[0]) > 0 else None
            if f:
                _rst.pop(id(f), None)
                f.close()
            return Null
        B['Close'] = _close

        _rst = {}          # 読み取りの状態: ストリーム → [現在の記録 (行 + 改行)、位置]。tfwrite.f の ipoint / lrecl に相当

        def _rec_next(f, st):
            line = f.readline()
            if line == '':
                return False
            st[0] = line if line.endswith('\n') else line + '\n'
            st[1] = 0
            return True

        def _rd_string(f, st, char1, new, delim, null):
            """tfreadstringf: delim が空なら記録の残り (ReadString)。Word は delim で区切り、Character は 1 文字。EOF は None"""
            while True:
                if st[0] is None or st[1] >= len(st[0]):
                    if not _rec_next(f, st):
                        return None
                rest = st[0][st[1]:]
                nc = len(rest)
                if char1:
                    nw, nxt, is_ = 1, 2, 0
                elif delim:
                    if null:
                        is_ = 0
                        if rest[0] in delim:
                            nw, nxt = 0, 2
                        else:
                            e = next((i for i in range(1, nc) if rest[i] in delim), None)
                            nxt = (e + 2) if e is not None else nc + 1
                            nw = e if e is not None else nc
                    else:
                        is_ = next((i for i in range(nc) if rest[i] not in delim), None)
                        if is_ is None:
                            if new:
                                st[1] = len(st[0])
                                continue
                            nw, nxt, is_ = 0, nc + 1, nc
                        else:
                            e = next((i for i in range(is_ + 1, nc) if rest[i] in delim), None)
                            nxt = (e + 2) if e is not None else nc + 1
                            nw = (e if e is not None else nc) - is_
                else:
                    is_, nw, nxt = 0, nc, nc + 1
                break
            w = rest[is_:is_ + nw]
            if new and nxt > nc:
                st[1] = len(st[0])
            else:
                st[1] += nxt - 1
            if w.endswith('\n'):
                w = w[:-1]
            return w

        def _rd_expr(f, st):
            """tfread1: 現在位置から式を 1 つ読む (隣り合う項で切る。括弧が閉じなければ次の行も使う)。評価して返す"""
            while st[0] is None or not st[0][st[1]:].strip():
                if not _rec_next(f, st):
                    return S('EndOfFile')
            txt = st[0][st[1]:]
            while True:
                depth = 0
                inq = False
                cut = None
                last = ''
                i = 0
                while i < len(txt):
                    c = txt[i]
                    if inq:
                        if c == '\\':
                            i += 1
                        elif c == '"':
                            inq = False
                    elif c == '"':
                        inq = True
                    elif c in '([{':
                        depth += 1
                    elif c in ')]}':
                        depth -= 1
                    elif c in ' \t' and depth == 0 and last:
                        j = i
                        while j < len(txt) and txt[j] in ' \t':
                            j += 1
                        if j < len(txt) and txt[j] != '\n' and (last.isalnum() or last in '$_.")]}') \
                                and (txt[j].isalnum() or txt[j] in '$"({.'):
                            cut = i
                            break
                    if c not in ' \t':
                        last = c
                    i += 1
                if cut is None and (depth > 0 or inq):
                    line2 = f.readline()
                    if line2:
                        txt += line2
                        continue
                break
            if cut is None:
                piece = txt
                st[0], st[1] = txt, len(txt)
            else:
                piece = txt[:cut]
                st[0], st[1] = txt, cut
            return I.run(piece, 'script') if piece.strip() else Null

        def _read_item(f, st, item, opts):
            new, delim, null = opts
            if isinstance(item, Sym):
                t = item.name
                if t == 'Expression':
                    return _rd_expr(f, st)
                if t == 'String':
                    w = _rd_string(f, st, False, new, '', False)
                    return S('EndOfFile') if w is None else w
                if t in ('Word', 'Real', 'Character'):
                    w = _rd_string(f, st, t == 'Character', new, delim, null)
                    if w is None:
                        return S('EndOfFile')
                    if t != 'Real':
                        return w
                    if w == '':
                        return Null
                    try:
                        return float(w)
                    except ValueError:
                        return I.run(w, 'script')
            elif isinstance(item, list):
                out = []
                for it in item:
                    v = _read_item(f, st, it, opts)
                    if isinstance(v, Seq):
                        out.extend(v.items)
                    else:
                        out.append(v)
                return out
            elif isinstance(item, Expr) and item.head == S('Times') and len(item.args) == 2 and is_num(item.args[0]):
                return Seq([_read_item(f, st, item.args[1], opts) for _ in range(int(item.args[0]))])
            raise SADError('General::wrongval:  (n*) BinaryDouble, BinaryInteger, BinaryReal,  BinaryShortInteger, Byte, Character, Expression,  BinaryFloat, Real, String, or Word is for 2nd arg')

        def _read(a):
            if not a or not is_num(a[0]):
                # gs: 装置番号でないもの (開けなかった OpenRead の h など) は wrongtype
                raise SADError('General::wrongtype:  Argument must be Real number or List of two Reals in Read[%s]' % ','.join(I.to_input_string(x) if hasattr(I, 'to_input_string') else I.to_str(x) for x in a))
            f = I.files.get(int(a[0]))
            if f is None:
                return S('EndOfFile')
            st = _rst.setdefault(id(f), [None, 0])
            item = a[1] if len(a) > 1 else S('Expression')
            new, delim, null = True, ' ,\n\t\r', False
            for o in a[2:]:                       # オプション WordSeparators / ReadNewRecord / NullWords (SADHelp の Read)
                if isinstance(o, Rule) and isinstance(o.lhs, Sym):
                    if o.lhs.name == 'WordSeparators' and isinstance(o.rhs, str):
                        delim = o.rhs
                    elif o.lhs.name == 'ReadNewRecord' and is_num(o.rhs):
                        new = o.rhs != 0
                    elif o.lhs.name == 'NullWords' and is_num(o.rhs):
                        null = o.rhs != 0
            return _read_item(f, st, item, (new, delim, null))
        B['Read'] = _read

        def _skip(args):
            # Skip[unit, type, n]: type の項目を n 個読み飛ばし、最後に読んだ項目を返す (gs で確認: Skip[g, Word, 2] → 2 つ目の語、EOF なら EndOfFile)。
            # n=0 や Record などの型は gs も評価しない (引数の誤りのメッセージを出して未評価)
            if len(args) != 3 or not is_num(args[0]) or not is_num(args[2]) or int(args[2]) <= 0:
                return Expr(S('Skip'), list(args))
            r = Null
            for _ in range(int(args[2])):
                r = _read([args[0], args[1]])
                if isinstance(r, Sym) and r.name == 'EndOfFile':
                    break
            return r
        B['Skip'] = _skip
        B['ReadString'] = lambda a: _read([a[0], S('String')] + list(a[1:]))
        B['ReadList'] = lambda a: _readlist(a)

        def _readlist(a):
            path = a[0] if isinstance(a[0], str) else None
            if path is not None:
                if not os.path.isabs(path):
                    path = os.path.join(I.cwd, path)
                lines = open(path, encoding='utf-8', errors='replace').read().split()
            else:
                lines = I.files[int(a[0])].read().split()
            try:
                return [float(x) for x in lines]
            except ValueError:
                return lines

        def _get(a):
            path = a[0]
            if not os.path.isabs(path):
                path = os.path.join(I.cwd, path)
            if not os.path.exists(path):
                I.write_unit(6, f'?Get: file not found {a[0]}\n')
                return S('$Failed')
            src = open(path, encoding='utf-8', errors='replace').read()
            pk = I.globals.get('PackagesDirectory')
            if isinstance(pk, str) and pk and os.path.normcase(os.path.abspath(path)).startswith(os.path.normcase(os.path.abspath(pk))):
                I._pkg_loading += 1            # SAD 自身の Packages: 文脈 (BeginPackage / Begin) は無視して大域で読む (load_package と同じ)
                try:
                    return I.run(src, 'ffs')
                finally:
                    I._pkg_loading -= 1
            return I.run(src, 'ffs')
        B['Get'] = _get

        # ---- 文脈 (名前空間): BeginPackage / Begin / End / EndPackage / Context (init.n の定義と gs の実測に合わせた)
        def _ctxargs(name, a):
            return [_ctx_norm(x) for x in a if x is not Null]

        def _beginpackage(a):
            if I._pkg_loading or not I._ctx_track:
                return Expr(S('BeginPackage'), list(a))
            cs = _ctxargs('BeginPackage', a)
            if not cs:
                raise SADError('General::narg: Number of arguments is expected to be 1 or more in BeginPackage')
            I._pkg_stack.append(list(I._ctx_path))
            path = []
            for c in cs + ['System`']:
                if c not in path:
                    path.append(c)
            I._ctx_path = path
            I._ctx_publish()
            return [Sym(c) for c in path]
        B['BeginPackage'] = _beginpackage

        def _endpackage(a):
            if I._pkg_loading or not I._ctx_track:
                return Expr(S('EndPackage'), list(a))
            I._ctx_path = I._pkg_stack.pop() if I._pkg_stack else list(_ROOT_PATH)
            I._ctx_publish()
            return [Sym(c) for c in I._ctx_path]
        B['EndPackage'] = _endpackage

        def _begin(a):
            if I._pkg_loading or not I._ctx_track:
                return Expr(S('Begin'), list(a))
            cs = _ctxargs('Begin', a)
            if len(cs) != 1:
                raise SADError('General::narg: Number of arguments is expected to be 1 in Begin')
            I._ctx_stack.append(I._ctx_cur)
            I._ctx_cur = cs[0]
            I._ctx_publish()
            return Null
        B['Begin'] = _begin

        def _end(a):
            if I._pkg_loading or not I._ctx_track:
                return Expr(S('End'), list(a))
            I._ctx_cur = I._ctx_stack.pop() if I._ctx_stack else '`'
            I._ctx_publish()
            return Null
        B['End'] = _end
        from . import miscfn as _miscfn
        B['System'] = lambda a: _miscfn.system_command(I, a)
        def _exit(a):
            if getattr(I, '_forked_child', False):
                _exit_child(I, int(a[0]) if a and is_num(a[0]) and not isinstance(a[0], complex) else 0)
            raise SystemExit(0)
        B['Exit'] = _exit
        B['Date'] = lambda a: [float(x) for x in time.localtime()[:6]]
        def _datestring(a):
            # DateString[] は現在、DateString[t] は t (1900-01-01 0:0:0 からの秒、または {Y,m,d[,H[,M[,S[,offset]]]]}) を mm/dd/yyyy hh:mm:ss に (SADHelp、gs で確認)
            import datetime as _dt
            if not a or a[0] is Null:
                if len(a) > 1:
                    raise SADError('General::narg:  Number of arguments is expected to be 0 or 1 (+ option) in DateString$[%s]' % ','.join(I.to_str(x) for x in a))
                return time.strftime('%m/%d/%Y %H:%M:%S')
            if len(a) > 1:
                raise SADError('General::narg:  Number of arguments is expected to be 0 or 1 (+ option) in DateString$[%s]' % ','.join(I.to_str(x) for x in a))
            x = a[0]
            try:
                if is_num(x):
                    d = _dt.datetime(1900, 1, 1) + _dt.timedelta(seconds=float(x))
                elif isinstance(x, list) and 3 <= len(x) <= 7 and all(is_num(v) for v in x):
                    v = [int(t) for t in x[:6]] + [0] * (6 - min(len(x), 6))
                    d = _dt.datetime(v[0], v[1], v[2], v[3], v[4], v[5])
                else:
                    raise ValueError
            except (ValueError, OverflowError):
                raise SADError('General::wrongtype:  Argument must be Null or Real number or {Y,m,d[,H[,M[,S[,offset]]]]} in DateString$[%s]' % I.to_str(x, None, True))
            return d.strftime('%m/%d/%Y %H:%M:%S')
        B['DateString'] = _datestring
        B['Pause'] = lambda a: Null
        B['TkWait'] = lambda a: Null
        B['Update'] = lambda a: Null
        B['FFSFlags'] = lambda a: Null
        def _out(a):
            if len(a) == 1 and is_num(a[0]) and int(a[0]) in I.out_hist and float(a[0]) == int(a[0]):
                return I.out_hist[int(a[0])]
            return Expr(S('Out'), list(a))
        B['Out'] = _out

        # ---- 多項式フィット (SAD の PolynomialFit: {{x,y}...} と次数 → 係数 {c0, c1, ...})
        def _polyfit(a):
            data = a[0]
            n = int(a[1])
            xs = np.array([d[0] for d in data], dtype=float)
            ys = np.array([d[1] for d in data], dtype=float)
            ws = np.array([d[2] if len(d) > 2 else 1.0 for d in data], dtype=float)
            # 最小二乗 (重み w: 残差に w を掛ける)
            A = np.vander(xs, n + 1, increasing=True) * ws[:, None]
            c, *_ = np.linalg.lstsq(A, ys * ws, rcond=None)
            return [float(x) for x in c]
        B['PolynomialFit$numpy'] = _polyfit
        # Spline / SplineData は Packages/D.n (自動読み込み)。Spline$ (自然スプライン y''/6) は miscfn.py

        # ---- リスト・式操作 (SAD 互換の挙動は sadlist.py にまとめた。同名の上の定義を上書きする)
        from . import sadlist
        sadlist.install(I)
        from . import fitroot
        fitroot.install(I)
        from . import miscfn
        miscfn.install(I)
        miscfn.install_matrix(I)
        miscfn.install_small(I)
        miscfn.install_system(I)
        miscfn.install_spline(I)
        from . import fileops
        fileops.install(I)
        fileops.install_misc(I)
        fileops.install_parallel(I)
        from . import radfield
        radfield.install(I)
        from . import csr
        csr.install(I)
        from . import vecfn
        from . import specfn
        specfn.install(I)
        vecfn.install(I)             # specfn の Erf・InverseErf なども Vector に対応させるため specfn の後
        self.autoload = {}
        self.subvalues = {}
        self.subvalues3 = {}
        self.attributes = {}
        self._install_autoload()
        from . import sadclass
        sadclass.install(I)
        from . import tclbridge
        tclbridge.install(I)


def _spec_score(pats):
    """定義の特殊さ (リテラル > 構造パターン > 型つき Blank・選択肢 > 素の Blank)。Mathematica のように、特殊な定義を先に試すための順序"""
    def one(p):
        if isinstance(p, _Lit):
            # リテラルのリスト {Automatic, Automatic} は、一部だけリテラルのパターン {Automatic, r2_} (20 + 100 + 0) より特殊
            return 21 + 100 * len(p.v) if isinstance(p.v, list) else 100
        if isinstance(p, Expr):
            h = p.head
            if h == S('Pattern') and len(p.args) == 2:
                return one(p.args[1])
            if h in (S('Blank'), S('BlankSequence'), S('BlankNullSequence')):
                # 長さの決まらない __・___ は、Blank を並べた定義 (f[_,_,_]) よりわずかに一般的 (OpticsPlot.n の OpticsList[__]={} は後の定義が先)
                return (10 if p.args else 0) - (0 if h == S('Blank') else 1)
            if h == S('Alternatives'):
                # リテラルだけの選択肢 ("" | " ") は個々のリテラルの定義と同じ扱い (K64 は Canvas.n の UnitLabelY[""|" ",__]=Null を先に試す)
                return 105 if p.args and all(isinstance(q, str) for q in p.args) else 5
            if h == S('Optional'):
                return one(p.args[0]) if p.args else 0
            if h in (S('PatternTest'), S('Condition')):
                # `{f_}?test` は検査つきの構造パターン (構造の特殊さ + 検査)。OpticsPlot.n の OpticsList1$[pr_,{f_}?(…)] は
                # OpticsList1$[{pr1_,pr2_},f_,…] より先に試す (gs は定義の順)
                return 8 + (one(p.args[0]) if p.args and isinstance(p.args[0], (Expr, list)) and not (isinstance(p.args[0], Expr) and p.args[0].head in (S('Blank'), S('Pattern'))) else 0)
            if _is_pattern(h):
                return 3 + sum(one(a) for a in p.args)         # f_[y__] のように head もパターン: p_Plus (型つき Blank) より一般的
            return 20 + sum(one(a) for a in p.args)
        if isinstance(p, (Sym, str, float, int)):
            return 100            # パターンの中のリテラル ({Default, y_} の Default) は、素の Blank ({x_, y_}) より特殊 (Plots.n の DP$DefMesh)
        return 0
    score = sum(one(p) for p in pats)
    # 同じ名前のパターン変数が繰り返される定義 (f[x_, x_]) は、型つき Blank (f[_Symbol, _]) より先に試す (K64 の D.n: D$[x_,x_]=1 と D$[_Symbol,_]=0)
    names = Interp._pattern_var_names(list(pats), set()) if False else None
    cnt = {}
    def walk(x):
        if isinstance(x, Expr):
            if x.head == S('Pattern') and x.args and isinstance(x.args[0], Sym):
                cnt[x.args[0].name] = cnt.get(x.args[0].name, 0) + 1
            for a_ in x.args:
                walk(a_)
        elif isinstance(x, list):
            for a_ in x:
                walk(a_)
    walk(list(pats))
    score += sum(12 * (k - 1) for k in cnt.values() if k > 1)
    return score


def _as_expr(x):
    """式の中の Python のリストを List[…] の式に戻す (純関数の本体は評価のたびに要素を評価する)"""
    if isinstance(x, list):
        return Expr(S('List'), [_as_expr(a) for a in x])
    if isinstance(x, Expr):
        return Expr(_as_expr(x.head), [_as_expr(a) for a in x.args])
    return x


def _unlit(x):
    """式の中の _Lit (評価済みの値を埋めるラッパー) を値に戻す"""
    if isinstance(x, _Lit):
        return _unlit(x.v)
    if isinstance(x, Expr):
        if x.head == S('List'):
            return [_unlit(a) for a in x.args]
        return Expr(_unlit(x.head), [_unlit(a) for a in x.args])
    if isinstance(x, list):
        return [_unlit(a) for a in x]
    if isinstance(x, Rule):
        return Rule(_unlit(x.lhs), _unlit(x.rhs), x.delayed)
    return x


class _Lit:
    """評価済みの値を AST に埋め込むためのラッパー (eval で中身をそのまま返す)"""
    __slots__ = ('v',)

    def __init__(self, v):
        self.v = v

    def __repr__(self):
        return f'Lit({self.v!r})'

    def __eq__(self, o):
        return isinstance(o, _Lit) and Interp.same(o.v, self.v)

    def __hash__(self):
        return hash(('Lit', repr(self.v)))


_orig_eval = Interp.eval


def _eval2(self, node):
    if type(node) is float:
        return node
    if isinstance(node, _Lit):
        return node.v
    try:
        return _orig_eval(self, node)
    except (TypeError, AttributeError) as e:          # 組み込み関数に想定外の型が来た: SAD と同じく式のエラーとして扱う (パッケージを動かす途中で止めない)
        if os.environ.get('PYSAD_TB'):
            import traceback
            traceback.print_exc()
        raise SADError('General::wrongtype: %s: %s' % (type(e).__name__, e)) from e


Interp.eval = _eval2
Interp._expr_reeval = True
_NoMatch = object()


def _mentions_hold(node):
    """式の中に Hold[…] があるか"""
    stack = [node]
    while stack:
        x = stack.pop()
        if isinstance(x, Expr):
            if x.head == S('Hold'):
                return True
            stack.append(x.head)
            stack.extend(x.args)
        elif isinstance(x, list):
            stack.extend(x)
    return False


def _is_pattern(a):
    if isinstance(a, Expr):
        if a.head in (S('Pattern'), S('Blank'), S('BlankSequence'), S('BlankNullSequence'), S('Optional'), S('PatternTest'), S('Alternatives')):
            return True
        return any(_is_pattern(x) for x in a.args)
    return False


def _has_pat(x):
    """パターンを (リストや a->b の中も含めて) 含むか"""
    if isinstance(x, list):
        return any(_has_pat(y) for y in x)
    if isinstance(x, Rule):
        return _has_pat(x.lhs) or _has_pat(x.rhs)
    if isinstance(x, Expr):
        return _is_pattern(x) or any(_has_pat(y) for y in x.args)
    return False


def _is_seq_pat(q):
    return isinstance(q, Expr) and (q.head in (S('BlankSequence'), S('BlankNullSequence'), S('Repeated'), S('RepeatedNull'), S('Optional'))
                                    or (q.head == S('Pattern') and _is_seq_pat(q.args[1])))      # Optional (Right[sep_:0]) も引数の数が変わる


def _same_pats(p1, p2):
    return repr(p1) == repr(p2)


def _copy(v):
    if isinstance(v, list):
        return [_copy(x) for x in v]
    return v


def _hkey(v):
    if isinstance(v, list):
        return ('L',) + tuple(_hkey(x) for x in v)
    if isinstance(v, (Sym, Rule, Expr)):
        return repr(v)
    return v


def _dims(v):
    d = []
    while isinstance(v, list):
        d.append(float(len(v)))
        if not v:
            break
        v = v[0]
    return d


def _has_complex(a):
    if isinstance(a, list):
        return any(_has_complex(x) for x in a)
    return isinstance(a, complex)


def _replace_part(lst, value, pos):
    out = _copy(lst)
    idx = pos if isinstance(pos, list) else [pos]
    v = out
    for k in idx[:-1]:
        v = v[int(k) - 1]
    v[int(idx[-1]) - 1] = value
    return out


@functools.lru_cache(maxsize=1024)
def _sad_wild_rx(p):
    """StringMatchQ のワイルドカード (SADHelp の wildcards、gs で確認): * (0 個以上)、% (1 文字)、{abc} (集合。a-c は範囲でなく a - c の 3 文字)、
    {^abc} (否定)、a|b (どちらか)。? は普通の文字。大文字小文字は区別する。閉じていない { は無視する"""
    rx = ''
    i = 0
    while i < len(p):
        c = p[i]
        if c == '*':
            rx += '.*'
        elif c == '%':
            rx += '.'
        elif c == '{':
            j = p.find('}', i)
            if j < 0:
                pass
            else:
                body = p[i + 1:j]
                if body.startswith('^'):
                    rx += '[^' + re.escape(body[1:]) + ']' if body[1:] else '.'
                else:
                    rx += '[' + re.escape(body) + ']' if body else '(?!)'
                i = j
        elif c == '|':
            rx += '|'
        else:
            rx += re.escape(c)
        i += 1
    return re.compile(rx, re.S)


def _wild(p):
    return ''.join('.*' if c == '*' else '.' if c == '?' else re.escape(c) for c in p)


# ----------------------------------------------------------------------------
# 記号 (未定義シンボル) を含む四則演算: SAD は引数を整列し、数の定数項を先頭に畳み込む (x+y+x → x+x+y、1+x+3 → 4+x)
_BEAMLINE = None


def _beamline_neg(b):
    """Packages/beamline.n: -BeamLine[x__] は順序を逆にして各要素の符号を反転 (逆向き通過)"""
    def neg(x):
        if isinstance(x, Expr) and x.head == S('Minus'):
            return x.args[0]
        if isinstance(x, (Sym, str)):
            return Expr(S('Minus'), [x])                 # ビームラインの逆向き要素は Minus[要素] (lattice の USE が読む形)
        return _symneg(x)
    return Expr(S('BeamLine'), [neg(x) for x in reversed(b.args)])


def _beamline_times(n, b):
    """n*BeamLine[x__]: n>=0 なら n 回繰り返し、n<0 なら Abs[n]*(-BeamLine)"""
    if n < 0:
        b = _beamline_neg(b)
        n = -n
    return Expr(S('BeamLine'), [x for _ in range(int(n)) for x in b.args])


def _symneg(b):
    if is_num(b):
        return -b
    if isinstance(b, Expr) and b.head == S('BeamLine'):
        return _beamline_neg(b)
    if isinstance(b, Expr) and b.head == S('Minus'):
        return b.args[0]
    if isinstance(b, Expr) and b.head == S('Times') and b.args:
        if is_num(b.args[0]):
            c = -b.args[0]
            rest = b.args[1:]
            if c == 1.0:
                return rest[0] if len(rest) == 1 else Expr(S('Times'), rest)      # -(-x) = x
            return Expr(S('Times'), [c] + rest)
        return Expr(S('Times'), [-1.0] + list(b.args))      # -(a b) = Times[-1, a, b] (入れ子にしない)
    return Expr(S('Times'), [-1.0, b])                    # K64: -x は Times[-1, x] (Length[-x] は 2)


def _flat_head(name, items):
    out = []
    for x in items:
        if isinstance(x, Expr) and x.head == S(name):
            out.extend(x.args)
        else:
            out.append(x)
    return out


def _sym_assoc(name, items, unit, absorb=None):
    from .sadlist import ckey
    items = _flat_head(name, items)
    const = unit
    rest = []
    for x in items:
        if is_num(x):
            const = const + x if name == 'Plus' else const * x
        else:
            rest.append(x)
    if absorb is not None and is_num(const) and const == absorb:
        return absorb
    rest.sort(key=ckey)
    terms = ([const] if const != unit else []) + rest
    if not terms:
        return unit
    if len(terms) == 1:
        return terms[0]
    return Expr(S(name), terms)


def _symb(name, a, b):
    """lmap2 の非数値オペランド用: 記号式を作る (数値が揃ったら評価時に畳まれる)"""
    if name == 'Plus':
        return _sym_assoc('Plus', [a, b], 0.0)
    if name == 'Times':
        if is_num(a) and isinstance(b, Expr) and b.head == S('BeamLine'):
            return _beamline_times(a, b)
        if is_num(b) and isinstance(a, Expr) and a.head == S('BeamLine'):
            return _beamline_times(b, a)
        return _sym_assoc('Times', [a, b], 1.0, 0.0)
    if name == 'Subtract':
        return _sym_assoc('Plus', [a, _symneg(b)], 0.0)
    if name == 'Divide':
        if is_num(b) and not isinstance(b, complex):
            return _sym_assoc('Times', [a, (1.0 / b) if b != 0 else INF], 1.0, 0.0)       # x/2 は (.5 x)
        return _sym_assoc('Times', [a, Expr(S('Power'), [b, -1.0])], 1.0, 0.0)
    if name == 'Power':
        if is_num(b) and b == 1:
            return a
        if is_num(b) and not isinstance(b, complex) and isinstance(a, Expr) and a.head == S('Times') and a.args and is_num(a.args[0]) \
                and not isinstance(a.args[0], complex) and a.args[0] != 0 and a.args[0] > 0:
            rest = a.args[1:]       # (2 x)^2 は (4 x^2): 数の係数だけ外へ出す (K64 で確認。負の係数は出さない)
            r_ = rest[0] if len(rest) == 1 else Expr(S('Times'), rest)
            return _sym_assoc('Times', [a.args[0] ** b, _symb('Power', r_, b)], 1.0, 0.0)
        if is_num(b) and b == 0:
            return 1.0
        return Expr(S('Power'), [a, b])
    return Expr(S('$Op'), [a, b])


# ----------------------------------------------------------------------------
# パターン照合で Expr の args に _Lit が混ざった場合のための表示
def expr_to_source(v) -> str:
    """デバッグ用"""
    return repr(v)
