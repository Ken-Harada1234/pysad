"""
sadelem.py — SADScript の FFS 専用関数のうち、要素とビームラインに関するもの

  Element[key, element(s)]            要素 (ファミリー) の値の取得・代入 (tfelement)
  SetElement[name, type, opts]        MAIN レベルの要素の作成・変更・参照 (Packages/element.n の SetElement$)
  FLAG[name]                          FFS フラグの値 (tfflag)
  ExtractBeamLine[line]               展開したビームラインを BeamLine[...] で返す (tfextractbeamline)
  BeamLine[...]                       BeamLine オブジェクト (Packages/beamline.n)。Join/Reverse/-/n* が使える
  PrintBeamLine / WriteBeamLine       ビームラインの書き出し (Packages/beamline.n。FFS 形式と MAIN 形式)
  GetMAIN[file]                       MAIN 形式のファイルを読む
  USE / VISIT <BeamLine 式>           式から新しい LINE (名前 L000001$ ...) を作って使う。VISIT/BYE は直前の LINE に戻る

挙動は SAD K64 (gs) で probe9〜12 / bench_level24f.sad を実行して得た出力に合わせている (docs/history/work_item4.md §F)。
SAD の癖: Element は現在のビームライン中のファミリーの一覧 (重複なし、初出順) を相手にする。単一のマッチはスカラー、
リスト指定 ({"QF"}、{2}) はリストで返す。キーワードがその型に無ければ 0、要素が無ければ {}。
"""
from __future__ import annotations

import math
import re
from typing import Any, List

from . import sadkw
from .sadscript import Sym, Expr, Rule, Null, SADError, S, is_num, _Lit


def _fill_string(s: str, n: int) -> str:
    """Packages/beamline.n: Fill$String[s,n] = (s//"                     "[1,n])[1, n*Floor[(StringLength[s]-1)/n]+n]"""
    pad = (' ' * 21)[:n]
    k = n * ((len(s) - 1) // n) + n
    return (s + pad)[:k]


def _remove_paren(s: str) -> str:
    return s[1:-1] if len(s) >= 2 and s[0] == '(' and s[-1] == ')' else s


def install(ffs):
    """FFS インスタンスに関数を登録 (ffscmd.FFS.install から呼ばれる)"""
    from . import ffscmd as _fc
    I = ffs.I
    B = I.builtins

    def tname(el) -> str:
        return _fc._TYPENAME.get(_fc._cls(el), _fc._cls(el).upper())

    def families() -> List[str]:
        seen: List[str] = []
        for el in ffs.els:
            if el.name not in seen:
                seen.append(el.name)
        return seen

    def first_el(fam: str):
        for el in ffs.els:
            if el.name == fam:
                return el
        return None

    def keyword_list(fam: str) -> List[str]:
        return list(sadkw.KEYWORDS.get(tname(first_el(fam)), []))

    def default_key(fam: str) -> str:
        return sadkw.DEFAULT_KEY.get(tname(first_el(fam)), '')

    def value_of(fam: str, K: str):
        el = first_el(fam)
        tn = tname(el)
        if K == 'NAME':
            return fam
        if K == 'POSITION':
            return float(families().index(fam) + 1)
        if K == 'COMPONENT':
            return [float(i + 1) for i, e in enumerate(ffs.els) if e.name == fam]     # LINE 上の位置 (成分) の一覧
        if K == 'KEYWORDS':
            return list(sadkw.KEYWORDS.get(tn, []))
        if K == 'DEFAULT':
            return sadkw.DEFAULT_KEY.get(tn, '')
        if K == 'TYPE':
            return float(sadkw.TYPECODE.get(tn, 0.0))
        if K == 'TYPENAME':
            return tn
        if K == 'VALUE':
            dk = sadkw.DEFAULT_KEY.get(tn, '')
            return value_of(fam, dk) if dk else 0.0
        if K not in sadkw.KEYWORDS.get(tn, []):
            return 0.0
        try:
            return float(ffs.get_value(el, K))
        except (TypeError, ValueError):
            return 0.0

    def _wild(name: str, pat: str) -> bool:
        return _fc._wild_match(name, pat)

    def match_fams(spec) -> List[str]:
        """spec (文字列 / 数) に合うファミリー名のリスト"""
        fams = families()
        if isinstance(spec, Sym):
            spec = spec.name
        if isinstance(spec, str):
            if any(c in spec for c in '*%{|?'):
                return [f for f in fams if _wild(f, spec)]
            return [f for f in fams if f.upper() == spec.upper()]
        if is_num(spec):
            n = len(fams)
            k = int(math.floor(spec))
            if k < 0:
                k = n + k + 1
            if not 1 <= k <= n:
                raise SADError('General::wrongnum: Number must be positive and less than length of beam line in Element$')
            return [fams[k - 1]]
        raise SADError('General::wrongtype: Argument must be element name or position')

    def resolve(spec, K: str):
        if spec is None:
            return [value_of(f, K) for f in families()]
        if isinstance(spec, list):
            return [resolve(s, K) for s in spec]
        ms = match_fams(spec)
        if not ms:
            return []
        if len(ms) == 1:
            return value_of(ms[0], K)
        return [value_of(f, K) for f in ms]

    def _strip_opts(args):
        return [a for a in args if not (isinstance(a, Rule) and isinstance(a.lhs, Sym) and a.lhs.name == 'Saved')]

    def _keyname(key) -> str:
        if isinstance(key, Sym):
            return key.name
        if isinstance(key, str):
            return key
        raise SADError('General::wrongtype: Argument must be Keyword for #1 in Element$[]')

    def fn_element(args):
        args = _strip_opts(list(args))
        if not args:
            raise SADError('General::wrongtype: Argument must be Keyword for #1 in Element$[]')
        key = args[0]
        if isinstance(key, list):
            return [fn_element([k] + args[1:]) for k in key]
        K = _keyname(key).upper()
        if len(args) > 2:
            return Expr(S('Element$'), args)
        if ffs.lat is None:
            raise SADError('No beam line (USE first)')
        if K == 'LENGTH':
            if len(args) != 1:
                raise SADError('General::narg: Number of arguments is expected to be 1 in Element$["LENGTH"]')
            return float(len(families()))
        if K == 'EXPAND':
            return Null
        if K == 'UPDATE':
            return [0.0] * len(families())
        return resolve(args[1] if len(args) > 1 else None, K)
    B['Element'] = fn_element

    def set_element(args, v):
        args = _strip_opts(list(args))
        key = args[0]
        K = _keyname(key).upper()
        spec = args[1] if len(args) > 1 else None
        if spec is None:
            fams = families()
        elif isinstance(spec, list):
            fams = [f for s in spec for f in match_fams(s)]
        else:
            fams = match_fams(spec)
        vals = v if isinstance(v, list) else None
        for i, fam in enumerate(fams):
            x = vals[i] if vals is not None and i < len(vals) else v
            kk = K
            if kk == 'VALUE':
                kk = default_key(fam)
                if not kk:
                    raise SADError('General::wrongtype: element has no default keyword in Element$')
            ffs.set_family(fam, kk, float(x))
        if isinstance(spec, list) and vals is None:
            return [v] * len(spec)        # Element["K1",{"QF","QD"}]=0.45 は {.45,.45} を返す
        return v
    ffs.set_element = set_element
    I.set_call_hooks['Element'] = set_element

    # ------------------------------------------------------------------ FLAG
    def fn_flag(args):
        if len(args) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in FLAG')
        a = args[0]
        if isinstance(a, list):
            return [fn_flag([x]) for x in a]
        name = a.name if isinstance(a, Sym) else a
        if not isinstance(name, str):
            raise SADError('General::wrongtype: Argument must be Character-string in FLAG')
        K = name                    # SAD は FLAG["cell"] を未定義のフラグとする (大文字のみ)
        if K in ffs.flags:
            return 1.0 if ffs.flags[K] else 0.0
        if K in _fc._SINO:
            return 0.0 if ffs.flags.get(_fc._SINO[K]) else 1.0
        raise SADError(f'FFS::undefflag: Undefined flag in FLAG["{name}"]')
    B['FLAG'] = fn_flag

    # ------------------------------------------------------------------ SetElement
    def _flat_rules(xs):
        out = []
        for x in xs:
            if isinstance(x, list):
                out.extend(_flat_rules(x))
            elif isinstance(x, Rule):
                out.append(x)
            else:
                raise SADError('General::wrongtype: Argument must be Rule or List of Rules in SetElement$')
        return out

    def fn_setelement(args):
        if not args:
            raise SADError('General::narg: Number of arguments is expected to be 1 or more in SetElement')
        nm = args[0]
        sym_name = isinstance(nm, Sym)
        if sym_name:
            nm = nm.name
        if not isinstance(nm, str):
            raise SADError('General::wrongtype: Argument must be Character-string in SetElement$')
        NAME = nm.upper()
        ret_name = '' if sym_name else nm        # SAD はシンボルで渡すと戻り値の名前が "" になる (要素は変更される)
        typ = args[1] if len(args) > 1 else None
        opts = _flat_rules(args[2:])
        env = ffs.env
        old = env._elements.get(NAME)
        if isinstance(typ, Sym):
            typ = typ.name
        if is_num(typ):
            inv = {int(c): n for n, c in sadkw.TYPECODE.items()}
            typ = inv.get(int(typ))
        if typ is Null:
            typ = None
        tn = typ.upper() if isinstance(typ, str) else None
        if old is None:
            if tn is None:
                return [ret_name, '']
            if tn not in sadkw.KEYWORDS:
                raise SADError(f'MAIN::wrongtype: Wrong type of symbol {typ} in MAIN in SetElement$["{nm}","{typ}"]')
            env.define_element(NAME, tn, {})
            old = env._elements[NAME]
        elif tn is not None and old.get('TYPE', '').upper() != tn:
            raise SADError(f'FFS::equaltype: Type {tn} must be equal to the previous definition in SetElement$["{nm}","{typ}"]')
        tn = old.get('TYPE', 'DRIFT').upper()
        fam_in_line = first_el(NAME) is not None
        for r in opts:
            k = r.lhs.name if isinstance(r.lhs, Sym) else r.lhs
            if not isinstance(k, str):
                raise SADError('General::wrongtype: Argument must be keyword in SetElement$')
            K = k.upper()
            val = I.eval(r.rhs) if r.delayed else r.rhs
            if not is_num(val):
                raise SADError('General::wrongtype: Argument must be Real number in SetElement$')
            old[K] = float(val)
            if fam_in_line:
                ffs.set_family(NAME, K, float(val))
        keys = sadkw.KEYWORDS.get(tn, [])
        vals = []
        for k in keys:
            if fam_in_line and tn != 'MARK':
                x = ffs.get_value(first_el(NAME), k)
            else:
                x = old.get(k, 0.0)        # MARK のキーワードは MAIN の定義値 (CALC で書き戻される Twiss ではない)
            try:
                x = float(x if x is not None else 0.0)
            except (TypeError, ValueError):
                x = 0.0
            vals.append(Rule(k, x))
        return [ret_name, tn, vals]
    B['SetElement'] = fn_setelement

    # ------------------------------------------------------------------ BeamLine
    BL = S('BeamLine')

    def _neg(x):
        if isinstance(x, Expr) and x.head == S('Minus') and len(x.args) == 1:
            return x.args[0]
        return Expr(S('Minus'), [x])

    def fn_beamline(args):
        out: List[Any] = []
        for a in args:
            if isinstance(a, Expr) and a.head == BL:
                out.extend(a.args)
            elif isinstance(a, Expr) and a.head == S('Times') and len(a.args) == 2 and is_num(a.args[0]):
                n = int(a.args[0])               # BeamLine[2*"DS"] → "DS","DS" (負なら -x を |n| 個)
                x = a.args[1]
                if n < 0:
                    x = _neg(x)
                out.extend([x] * abs(n))
            else:
                out.append(a)
        return Expr(BL, out)
    B['BeamLine'] = fn_beamline

    def _expand_line(name: str) -> List[Any]:
        env = ffs.env
        refs = env._beamlines.get(env.resolve_beamline(name) or name.upper())
        if refs is None:
            raise SADError(f'MAIN::wrongtype: Wrong type of symbol {name} in MAIN in ExtractBeamLine["{name}"]')
        out: list = []
        env._expand_refs(refs, out)
        return [_neg(S(x[1])) if isinstance(x, tuple) else S(x) for x in out]

    def fn_extract(args):
        if not args:
            if ffs.lat is None:
                raise SADError('No beam line (USE first)')
            return Expr(BL, [_neg(S(el.name)) if getattr(el, '_rev', False) else S(el.name) for el in ffs.els])
        a = args[0]
        if isinstance(a, Sym):
            a = a.name
        if not isinstance(a, str):
            raise SADError(f'General::wrongtype: Argument must be Character-string in ExtractBeamLine[{I.to_str(args[0])}]')
        return Expr(BL, _expand_line(a))
    B['ExtractBeamLine'] = fn_extract

    # ---- 書き出し (Packages/beamline.n WriteBeamLine1)
    def _write_unit(unit, text):
        I.write_unit(int(unit) if is_num(unit) else 6, text)

    def _write(unit, text):
        """Write[f, ...]: 末尾に改行 (文字列が改行で終わっていれば足さない)"""
        _write_unit(unit, text if text.endswith('\n') else text + '\n')

    def write_beamline1(f, b, opts, name):
        form = 'FFS'
        for o in opts:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'Format':
                form = o.rhs
        # Print["WriteBeamLine ",{{opt},form,name}] — SAD の Packages に残っているデバッグ出力 (標準出力に出る)
        nm_val = Sym('None') if name is None else name
        I.out.write(I.wrap('WriteBeamLine ' + I.to_str([list(opts), form, nm_val], quote=True)) + chr(10))
        items = list(b.args) if isinstance(b, Expr) else list(b)
        lb = len(items)
        txt = [_remove_paren(I.to_str(x)) for x in items]
        if form == 'MAIN':
            col = 9
            if name is not None:
                _write(f, '\nLINE ' + str(name) + ' = ')
            ch = '       ('
            for t in txt:
                _write_unit(f, _fill_string(ch + t, 9))
                col += 9
                ch = ' '
                if col > 69:
                    _write(f, '\n')
                    col = 9
                    ch = '        '
            _write(f, ')')
            if name is not None:
                _write(f, ';')
        else:
            col = 9
            if name is not None:
                _write(f, '\n' + str(name) + ' = BeamLine[')
            else:
                _write(f, 'BeamLine[')
            ch = '        '
            for i, t in enumerate(txt):
                s = _fill_string(ch + t + ('' if i == lb - 1 else ','), 9)
                _write_unit(f, s)
                col += len(s)
                ch = ''
                if col > 69:
                    _write(f, '')
                    col = 9
                    ch = '        '
            _write(f, '];')

    def fn_writebeamline(args):
        if not args:
            raise SADError('General::narg: Number of arguments is expected to be 1 or more in WriteBeamLine')
        f = args[0]
        rest = list(args[1:])
        opts = [x for x in rest if isinstance(x, Rule)]
        bls = [x for x in rest if not isinstance(x, Rule)]
        if len(bls) == 1 and isinstance(bls[0], list):
            bls = list(bls[0])
        elif len(bls) > 1:
            # WriteBeamLine[f_,b__BeamLine:ExtractBeamLine[],opt___] は Complement[Flatten[{b,opt}],op] で BeamLine を集める:
            # 標準順序に並べ替えられ、同じものは 1 つになる (引数の順に書き出されない)
            from .sadlist import _uniq_sorted
            bls = _uniq_sorted(bls)
        if not bls:
            bls = [fn_extract([])]
        for x in bls:
            if not (isinstance(x, Expr) and x.head == BL):
                raise SADError('General::wrongtype: Argument must be BeamLine in WriteBeamLine')
        names = None
        for o in opts:
            if isinstance(o.lhs, Sym) and o.lhs.name == 'Name':
                names = o.rhs
        if names is None:
            names = [None] * len(bls)
        elif not isinstance(names, list):
            names = [names]
        if len(names) < len(bls):
            I.out.write('WriteBeamLine - Insufficient number of Names' + chr(10))
            return Null
        for b, nm in zip(bls, names):
            write_beamline1(f, b, [Rule(S('Name'), Sym('None') if nm is None else nm)] + opts, nm)
        return Null
    B['WriteBeamLine'] = fn_writebeamline
    B['PrintBeamLine'] = lambda a: fn_writebeamline([6.0] + list(a))

    # ------------------------------------------------------------------ GetMAIN
    def fn_getmain(args):
        if len(args) != 1 or not isinstance(args[0], str):
            raise SADError('General::wrongtype: Argument must be Character-string in GetMAIN')
        import os
        path = args[0]
        if not os.path.isabs(path) and getattr(I, 'cwd', None):
            path = os.path.join(I.cwd, path)
        try:
            src = open(path, encoding='utf-8', errors='replace').read()
        except OSError:
            raise SADError(f'GetMAIN: cannot open {args[0]}')
        from . import parser as _parser
        _parser._Parser(_parser.tokenize(src), ffs.env).parse_program()
        return Null
    B['GetMAIN'] = fn_getmain
    return fn_beamline
