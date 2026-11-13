"""SAD のクラス (Class[...]、x = c[引数]、x@メンバー、This、Constructor / Destructor)。SADHelp の「Class」の節に合わせた実装。

方式: 「クラスのメンバーが出てきたら This@メンバー と読む」という SAD の規則を、メンバーの記号をインスタンスごとに一意な記号
(`c$i12` のようなもの) へ置き換えることで実現する。インスタンスの作成時に
  - インスタンス変数・クラス変数・メソッドの名前を、そのインスタンス (クラス) 専用の記号へ写す対応表を作り、
  - メソッドの本体・パターンの中の記号を写して、写した名前の定義 (DownValue) として登録する
ので、メソッドの本体は普通の関数として評価でき、ウィジェットの生成 (`c = Canvas[...]`) や `coords[0] = …` のように、
記号に定義を付ける処理も、そのままインスタンスごとの記号に対して働く。

対応している範囲: Class[親の列, クラス変数, インスタンス変数, メソッド…]、メソッドの形 (f[x_] := g、With[_, f[x_] := g]、If[_, …, …])、
x = c[引数] (規則 a->1 はメンバーの初期値、他の規則は Options)、Constructor / Destructor、x@メンバー、x@メソッド[引数]、x@メンバー = 値、
This、親クラスの継承、x = c1 の複写 (同じインスタンスを指す)、x=. の解体。
対応していない範囲: `親`メンバー の書き方、h_[f_[arg___], b___] ^:= g 形の上位値のメソッド、Members[] / Parents[] などの問い合わせ。
"""
from __future__ import annotations

from .sadscript import (Expr, Frame, Rule, Sym, S, Null, SADError, _Lit, _is_pattern, _same_pats, _spec_score, _Return)

_LOCAL_SCOPES = ('Module', 'Block', 'With')


class SadClass:
    def __init__(self, cid, supers):
        self.cid = cid
        self.supers = supers            # [SadClass]
        self.cvars = []                 # [(名前, 初期値の式 or None)]
        self.ivars = []                 # [(名前, 初期値の式 or None)]
        self.methods = []               # [(名前, パターンの列, 本体, 遅延か)]
        self.cmap = {}                  # クラス変数の名前 → 記号名 (親の分を含む)

    def __repr__(self):
        return 'SadClass#%d' % self.cid


class SadInstance:
    def __init__(self, cls, iid):
        self.cls = cls
        self.iid = iid
        self.map = {}                   # メンバーの名前 → このインスタンスの記号名
        self.syms = []                  # このインスタンスのために作った記号名 (解体で消す)

    def __repr__(self):
        return '%s@%d' % (self.cls, self.iid)


def install(I):
    B = I.builtins
    state = {'cid': 0, 'iid': 0}

    # ------------------------------------------------------------------ 宣言の読み取り
    def _decl(node):
        """{a, b=1, {c, d}={2, 3}, e=f=0} → [(名前, 初期値の式 or None)]。初期値は先頭の変数に付ける。{c, d}={2,3} は要素ごと"""
        out = []
        items = list(node.args) if isinstance(node, Expr) and node.head == S('List') else [node]
        for it in items:
            if isinstance(it, Sym):
                out.append((it.name, None))
            elif isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')) and len(it.args) == 2:
                names = []
                cur = it
                while isinstance(cur, Expr) and cur.head == S('Set') and len(cur.args) == 2:
                    names.append(cur.args[0])
                    cur = cur.args[1]
                val = cur
                if len(names) == 1 and isinstance(names[0], Expr) and names[0].head == S('List') and isinstance(val, Expr) and val.head == S('List') \
                        and len(names[0].args) == len(val.args):
                    for nm, v in zip(names[0].args, val.args):
                        if isinstance(nm, Sym):
                            out.append((nm.name, v))
                    continue
                if it.head == S('SetDelayed'):
                    nm = it.args[0]
                    if isinstance(nm, Sym):
                        out.append((nm.name, Expr(S('$Delayed'), [it.args[1]])))
                    continue
                for nm in names:
                    if isinstance(nm, Sym):
                        out.append((nm.name, val))
        return out

    # ------------------------------------------------------------------ メソッド定義の収集
    def _collect(node, cls, frames_note):
        """クラスの本体の中のメソッド定義 (f[...] := g、With[_, …]、If[_, …]、a; b) を cls.methods に集める"""
        if node is None or node is Null:
            return
        if isinstance(node, Expr):
            h = node.head
            if h == S('CompoundExpression'):
                for a in node.args:
                    _collect(a, cls, frames_note)
                return
            if h in (S('SetDelayed'), S('Set')) and len(node.args) == 2:
                lhs, rhs = node.args
                if isinstance(lhs, Expr) and isinstance(lhs.head, Sym):
                    pats = [a if _is_pattern(a) else a for a in lhs.args]
                    cls.methods.append((lhs.head.name, pats, rhs, h == S('SetDelayed')))
                elif isinstance(lhs, Sym):
                    cls.methods.append((lhs.name, None, rhs, h == S('SetDelayed')))
                return
            if h == S('With') and len(node.args) == 2:
                # With[{a = 1}, f[x_] := …]: 宣言の値を今の環境で求めて、本体の定義へ埋め込む
                fr = Frame()
                fr.immut = True
                decl = node.args[0]
                items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
                for it in items:
                    if isinstance(it, Expr) and it.head == S('Set') and isinstance(it.args[0], Sym):
                        fr[it.args[0].name] = I.eval(it.args[1])
                I.frames.append(fr)
                try:
                    before = len(cls.methods)
                    _collect(node.args[1], cls, frames_note)
                    for k in range(before, len(cls.methods)):
                        nm, pats, rhs, dl = cls.methods[k]
                        cls.methods[k] = (nm, pats, I._close(rhs), dl)
                finally:
                    I.frames.pop()
                return
            if h == S('If') and len(node.args) >= 2:
                c = I.eval(node.args[0])
                from .sadscript import truthy
                if truthy(c) is True:
                    _collect(node.args[1], cls, frames_note)
                elif len(node.args) >= 3:
                    _collect(node.args[2], cls, frames_note)
                return
            if h in (S('UpSetDelayed'), S('UpSet')):
                return                      # 上位値形式のメソッドは未対応
        # 他の式は無視する

    def _class(args):
        if len(args) < 3:
            raise SADError('General::narg: Number of arguments is expected to be 3 or more in Class')
        state['cid'] += 1
        sup = I.eval(args[0])
        supers = [s for s in (sup if isinstance(sup, list) else [sup]) if isinstance(s, SadClass)]
        cls = SadClass(state['cid'], supers)
        for sc in supers:
            cls.cmap.update(sc.cmap)
        for nm, init in _decl(args[1]):
            sym = '%s$c%d' % (nm, cls.cid)
            cls.cmap[nm] = sym
            cls.cvars.append((nm, init))
        cls.ivars = _decl(args[2])
        for a in args[3:]:
            _collect(a, cls, None)
        # クラス変数の初期化 (クラスのメンバーとして見えるのはインスタンスの対応表を通してだが、初期値は大域で評価する)
        for nm, init in cls.cvars:
            sym = cls.cmap[nm]
            if init is not None:
                I.assign(sym, I.eval(init))
        return cls
    B['Class'] = _class

    # ------------------------------------------------------------------ 記号の置き換え
    def _pat_names(x, acc):
        if isinstance(x, Expr):
            if x.head == S('Pattern') and x.args and isinstance(x.args[0], Sym):
                acc.add(x.args[0].name)
            for a in x.args:
                _pat_names(a, acc)
        elif isinstance(x, list):
            for a in x:
                _pat_names(a, acc)
        return acc

    def _local_names(decl):
        names = set()
        items = decl.args if isinstance(decl, Expr) and decl.head == S('List') else [decl]
        for it in items:
            if isinstance(it, Sym):
                names.add(it.name)
            elif isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')):
                lhs = it.args[0]
                if isinstance(lhs, Sym):
                    names.add(lhs.name)
                elif isinstance(lhs, Expr) and lhs.head == S('List'):
                    for x in lhs.args:
                        if isinstance(x, Sym):
                            names.add(x.name)
        return names

    def _rename(node, mp, excl, inst):
        if isinstance(node, Sym):
            if node.name == 'This':
                return _Lit(inst)
            if node.name in mp and node.name not in excl:
                return Sym(mp[node.name])
            return node
        if isinstance(node, Expr):
            h = node.head
            if isinstance(h, Sym) and h.name in _LOCAL_SCOPES and node.args and isinstance(node.args[0], (Expr, Sym)):
                loc = _local_names(node.args[0])
                decl = node.args[0]
                # 宣言の初期値は外側の環境、本体は局所名を除いた環境
                if isinstance(decl, Expr) and decl.head == S('List'):
                    items = []
                    for it in decl.args:
                        if isinstance(it, Expr) and it.head in (S('Set'), S('SetDelayed')):
                            items.append(Expr(it.head, [it.args[0], _rename(it.args[1], mp, excl, inst)]))
                        else:
                            items.append(it)
                    ndecl = Expr(S('List'), items)
                else:
                    ndecl = decl
                ex2 = excl | loc
                return Expr(h, [ndecl] + [_rename(a, mp, ex2, inst) for a in node.args[1:]])
            if isinstance(h, Sym) and h.name == 'Function' and len(node.args) >= 2 and isinstance(node.args[0], (Sym, Expr)):
                ex2 = excl | _local_names(node.args[0] if isinstance(node.args[0], Expr) and node.args[0].head == S('List') else Expr(S('List'), [node.args[0]]))
                return Expr(h, [node.args[0]] + [_rename(a, mp, ex2, inst) for a in node.args[1:]])
            nh = _rename(h, mp, excl, inst) if not isinstance(h, Sym) or h.name in mp else h
            return Expr(nh, [_rename(a, mp, excl, inst) for a in node.args])
        if isinstance(node, list):
            return [_rename(a, mp, excl, inst) for a in node]
        if isinstance(node, Rule):
            return Rule(_rename(node.lhs, mp, excl, inst), _rename(node.rhs, mp, excl, inst), node.delayed)
        return node

    # ------------------------------------------------------------------ インスタンスの作成
    def _mro(cls):
        out = []
        for sc in cls.supers:
            for c in _mro(sc):
                if c not in out:
                    out.append(c)
        out.append(cls)
        return out

    def _make_instance(cls, args):
        state['iid'] += 1
        iid = state['iid']
        inst = SadInstance(cls, iid)
        mro = _mro(cls)
        mp = {}
        for c in mro:
            mp.update(c.cmap)
        ivars = []
        for c in mro:
            for nm, init in c.ivars:
                ivars.append((nm, init))
        names = set(mp)
        for nm, _ in ivars:
            mp[nm] = '%s$i%d' % (nm, iid)
        mp['Options'] = 'Options$i%d' % iid
        meth_names = []
        for c in mro:
            for nm, _p, _b, _d in c.methods:
                if nm not in meth_names:
                    meth_names.append(nm)
        for nm in meth_names:
            if nm not in mp:
                mp[nm] = '%s$i%d' % (nm, iid)
        inst.map = mp
        inst.syms = [v for k, v in mp.items() if k not in names]
        saved = I.frames
        I.frames = []
        try:
            # インスタンス変数の初期値 (宣言の初期値、次に引数の規則 a->1)
            for nm, init in ivars:
                if init is not None:
                    node = _rename(init, mp, set(), inst)
                    if isinstance(node, Expr) and node.head == S('$Delayed'):
                        I.downvalues.setdefault(mp[nm], []).insert(0, (None, node.args[0], True))
                    else:
                        I.assign(mp[nm], I.eval(node))
            opts = []
            rest = []
            for a in args:
                if isinstance(a, Rule) and isinstance(a.lhs, Sym) and a.lhs.name in mp and a.lhs.name not in meth_names and a.lhs.name != 'Options':
                    if a.delayed:
                        I.downvalues.setdefault(mp[a.lhs.name], []).insert(0, (None, a.rhs, True))
                    else:
                        I.assign(mp[a.lhs.name], a.rhs)
                elif isinstance(a, Rule):
                    opts.append(a)
                    rest.append(a)
                else:
                    rest.append(a)
            I.assign(mp['Options'], opts)
            # メソッドの登録 (子クラスの定義を先に試す)
            for c in mro:                      # 親から順に登録し、子の同じ形の定義が親のものを置き換える
                for nm, pats, body, delayed in c.methods:
                    sym = mp[nm]
                    ex = _pat_names(pats, set()) if pats is not None else set()
                    nb = _rename(body, mp, ex, inst)
                    if pats is None:
                        I.downvalues.setdefault(sym, []).insert(0, (None, nb, True))
                        continue
                    npats = [_rename(p, mp, ex, inst) if _is_pattern(p) else _Lit(I.eval(_rename(p, mp, ex, inst))) for p in pats]
                    dv = I.downvalues.setdefault(sym, [])
                    dv[:] = [d for d in dv if not (d[0] is not None and _same_pats(d[0], npats))]
                    sc_ = _spec_score(npats)
                    pos = len(dv)
                    for k, d in enumerate(dv):
                        if d[0] is not None and _spec_score(d[0]) <= sc_:
                            pos = k
                            break
                    dv.insert(pos, (npats, nb, delayed if delayed else True))
        finally:
            I.frames = saved
        # Constructor
        ctor = mp.get('Constructor')
        if ctor and I.downvalues.get(ctor):
            r = I.apply_downvalues(ctor, I.downvalues[ctor], [a for a in rest if not isinstance(a, Rule)])
        return inst

    # ------------------------------------------------------------------ 評価器との接続
    orig_apply_value = I.apply_value

    def apply_value(hv, args, head=None):
        if isinstance(hv, SadClass):
            return _make_instance(hv, args)
        return orig_apply_value(hv, args, head)
    I.apply_value = apply_value

    def member_sym(obj, name):
        if isinstance(obj, SadInstance):
            sym = obj.map.get(name)
        else:
            sym = obj.cmap.get(name)
        if sym is None:
            raise SADError('Class::nomember: %s is not a member of %r' % (name, obj))
        return Sym(sym)
    I.member_sym = member_sym

    def _at(args):
        if len(args) != 2:
            return Null
        obj = I.eval(args[0])
        if isinstance(obj, (SadInstance, SadClass)) and isinstance(args[1], Sym):
            return I.eval(member_sym(obj, args[1].name))
        return I.eval(Expr(args[0], [args[1]]))
    B['At'] = _at

    def resolve_head(head):
        """Expr(At, [x, m]) → x が SadInstance なら m の記号。そうでなければ None"""
        if isinstance(head, Expr) and head.head == S('At') and len(head.args) == 2 and isinstance(head.args[1], Sym):
            obj = I.eval(head.args[0])
            if isinstance(obj, (SadInstance, SadClass)):
                return member_sym(obj, head.args[1].name)
        return None
    I.resolve_member_head = resolve_head

    def destroy(inst):
        d = inst.map.get('Destructor')
        if d and I.downvalues.get(d):
            I.apply_downvalues(d, I.downvalues[d], [])
        for nm in inst.map.values():
            I.globals.pop(nm, None)
            I.downvalues.pop(nm, None)
    I.destroy_instance = destroy
    I._sad_instance_type = SadInstance
