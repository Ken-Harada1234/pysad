"""
Packages/FileSystem.n のファイル・ディレクトリ関数 (K64 は System["cp ..."] などの外部コマンドで実装。pysad は Python の os・shutil)。

FileQ DirectoryQ FileNames DirectoryName ToFileName CreateDirectory DeleteDirectory RenameDirectory CopyDirectory
CopyFile RenameFile DeleteFile FileDate SetFileDate FileType FileByteCount と $HomeDirectory。

K64 との差: FileByteCount は K64 が `wc -l` (行数) を返す誤り (docs/21 B50) なので、pysad はバイト数を返す。
失敗は $Failed。ファイルの相対パスは Interp.cwd から。
"""
from __future__ import annotations

import datetime
import fnmatch
import os
import shutil

from .sadscript import Rule, Sym, S, Null, truthy

_T0 = datetime.datetime(1900, 1, 1)


def install(I):
    B = I.builtins
    FAILED = S('$Failed')

    def real(p):
        if not os.path.isabs(p):
            p = os.path.join(I.cwd or os.getcwd(), p)
        return p

    def isfile(p):
        return isinstance(p, str) and os.path.isfile(real(p))

    def isdir(p):
        return isinstance(p, str) and os.path.isdir(real(p))

    def mapped(f):
        def w(a):
            if len(a) == 1 and isinstance(a[0], list):
                return [f(x) for x in a[0]]
            if len(a) > 1 and all(isinstance(x, str) for x in a):
                return [f(x) for x in a]
            return f(a[0])
        return w
    B['FileQ'] = mapped(lambda x: 1.0 if isfile(x) else 0.0)
    B['DirectoryQ'] = mapped(lambda x: 1.0 if isdir(x) else 0.0)

    def filenames(a):
        names = sorted(os.listdir(I.cwd or os.getcwd()))
        if not a:
            return names
        if isinstance(a[0], list):
            return [[n for n in names if fnmatch.fnmatchcase(n, f)] for f in a[0]]
        return [n for n in names if fnmatch.fnmatchcase(n, a[0])]
    B['FileNames'] = filenames

    def dirname(a):
        d = a[0]
        for _ in range(int(a[1]) if len(a) > 1 else 1):
            core = d[:-1]
            i = core.rfind('/')
            d = '' if i < 0 else d[:i + 1]
        return d
    B['DirectoryName'] = dirname

    def tofilename(a):
        dirs = a[0] if isinstance(a[0], list) else [a[0]]
        path = [x.strip('/') for x in dirs]
        path = [x for x in path if x != '']
        first = next((x for x in dirs if x != ''), None)
        lead = '/' if first is not None and first.startswith('/') else ''
        base = lead + '/'.join(path) + ('/' if path else '')
        if len(a) > 1:
            return base + (a[1][1:] if a[1].startswith('/') else a[1])
        return base
    B['ToFileName'] = tofilename

    def ok(path):
        return os.path.realpath(real(path)).replace(chr(92), '/')

    def guarded(f):
        def w(a):
            try:
                return f(a)
            except OSError:
                return FAILED
        return w

    def mkdir(a):
        os.mkdir(real(a[0]))
        return a[0]
    B['CreateDirectory'] = guarded(mkdir)

    def rmdir(a):
        contents = any(isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'DeleteContents' and truthy(o.rhs) is True for o in a[1:])
        (shutil.rmtree if contents else os.rmdir)(real(a[0]))
        return Null
    B['DeleteDirectory'] = guarded(rmdir)

    def renamedir(a):
        if not isdir(a[0]) or isdir(a[1]):
            return FAILED
        os.rename(real(a[0]), real(a[1]))
        return ok(a[1])
    B['RenameDirectory'] = guarded(renamedir)

    def copydir(a):
        if not isdir(a[0]) or isdir(a[1]):
            return FAILED
        shutil.copytree(real(a[0]), real(a[1]))
        return ok(a[1])
    B['CopyDirectory'] = guarded(copydir)

    def copyfile(a):
        if not isfile(a[0]) or os.path.exists(real(a[1])):
            return FAILED
        shutil.copy2(real(a[0]), real(a[1]))
        return ok(a[1])
    B['CopyFile'] = guarded(copyfile)

    def renamefile(a):
        if not isfile(a[0]) or os.path.exists(real(a[1])):
            return FAILED
        os.rename(real(a[0]), real(a[1]))
        return ok(a[1])
    B['RenameFile'] = guarded(renamefile)

    def deletefile(a):
        for f in (a[0] if isinstance(a[0], list) else [a[0]]):
            os.remove(real(f))
        return Null
    B['DeleteFile'] = guarded(deletefile)

    def filedate(a):
        p = real(a[0])
        if not os.path.exists(p):
            return FAILED
        d = datetime.datetime.fromtimestamp(os.path.getmtime(p), datetime.timezone.utc)
        return [float(d.year), float(d.month), float(d.day), float(d.hour), float(d.minute), float(d.second)]
    B['FileDate'] = filedate

    def setfiledate(a):
        p = real(a[0])
        if not os.path.exists(p):
            return FAILED
        y, mo, d, h, mi, s = [int(v) for v in a[1]]
        t = datetime.datetime(y, mo, d, h, mi, s, tzinfo=datetime.timezone.utc).timestamp()
        os.utime(p, (t, t))
        return Null
    B['SetFileDate'] = setfiledate
    B['FileType'] = lambda a: S('Directory') if isdir(a[0]) else (S('File') if isfile(a[0]) else S('None'))
    B['FileByteCount'] = lambda a: float(os.path.getsize(real(a[0]))) if isfile(a[0]) else FAILED
    I.globals.setdefault('$HomeDirectory', os.path.expanduser('~').replace(chr(92), '/'))


def install_misc(I):
    """init.n・functions.n の小さな関数: ARRAY Len Put Environment HomeDirectory StringInsert Diagonal Expand SymbolNameRoot HeldPart MkSecureTemp"""
    import tempfile
    from .sadscript import Expr, _Lit
    B = I.builtins
    B['Len'] = lambda a: float(len(a[0])) if isinstance(a[0], str) else Expr(S('Len'), list(a))

    def array(a):
        dims = [int(x) for x in a]

        def build(k):
            return [0.0] * dims[0] if k == len(dims) - 1 else [build(k + 1) for _ in range(dims[k])]
        return [0.0] * dims[0] if len(dims) == 1 else [build(1) for _ in range(dims[0])] if len(dims) == 2 else _nd(dims)
    def _nd(dims):
        if len(dims) == 1:
            return [0.0] * dims[0]
        return [_nd(dims[1:]) for _ in range(dims[0])]
    B['ARRAY'] = lambda a: _nd([int(x) for x in a])

    def put(a):
        f = open(real(a[-1]), 'w', encoding='utf-8')
        for x in a[:-1]:
            f.write(I.wrap(I.to_str(x, quote=True)) + chr(10))
        f.close()
        return Null

    def real(p):
        return p if os.path.isabs(p) else os.path.join(I.cwd or os.getcwd(), p)
    B['Put'] = put
    B['Environment'] = lambda a: (os.environ.__setitem__(str(a[0]), str(a[1])), Null)[1] if len(a) == 2 else os.environ.get(str(a[0]), '')
    B['HomeDirectory'] = lambda a: os.environ.get('HOME', os.path.expanduser('~'))

    def stringinsert(a):
        x, y, n = a
        if isinstance(n, list):
            pos = sorted(int(v) for v in n)
            out = ''
            prev = 0
            for p_ in pos:
                out += x[prev:p_ - 1] + y
                prev = p_ - 1
            return out + x[prev:]
        k = int(n)
        return x[:k - 1] + y + x[k - 1:]
    B['StringInsert'] = stringinsert
    B['Diagonal'] = lambda a: [a[0][i][i] for i in range(min(len(a[0]), len(a[0][0])))]
    B['SymbolNameRoot'] = lambda a: _snr(a[0].name if isinstance(a[0], Sym) else str(a[0]))

    def _snr(n):
        i = n.rfind('$')
        if i >= 0 and (n[i + 1:] == '' or n[i + 1:].isdigit()):
            return n[:i]          # x$12 → x (K64: `$` の後が数字だけなら落とす)
        return n
    B['HeldPart'] = lambda a: Expr(S('Hold'), [I.part(a[0], [float(k) for k in a[1:]])])

    def expand(a):
        x = a[0]
        if not isinstance(x, Expr):
            return x
        if x.head == S('Plus'):
            return I.builtins['Plus']([expand([t]) for t in x.args])
        if x.head == S('Times'):
            terms = [[]]
            for f in x.args:
                if isinstance(f, Expr) and f.head == S('Plus'):
                    ts = [expand([t]) for t in f.args]
                    terms = [t0 + [t1] for t0 in terms for t1 in ts]
                else:
                    terms = [t0 + [f] for t0 in terms]
            prods = [I.builtins['Times'](t) if len(t) > 1 else t[0] for t in terms]
            return I.builtins['Plus'](prods) if len(prods) > 1 else prods[0]
        return x
    B['Expand'] = expand

    def mksecuretemp(a):
        fd, path = tempfile.mkstemp(prefix='mkstemp.')
        os.close(fd)
        return path.replace(chr(92), '/')
    B['MkSecureTemp'] = mksecuretemp


def install_parallel(I):
    """shared.n の Parallelize (逐次実行) と Shared / OpenShared$ / ReadShared / WriteShared (プロセス内の表)。
    K64 は fork して共有メモリで結果を集める。pysad は逐次: Parallelize[f, n, lshare, npara, opts] = {f[1], …, f[n]} (失敗は $Failed)。
    Completion->g があれば結果に g を適用する。Fork は未対応 (プロセスの分岐は移植しない)"""
    from .sadscript import SADError, Rule
    B = I.builtins
    store = {}
    counter = [1024]

    def openshared(a):
        counter[0] += 1
        store[counter[0]] = float('nan') if False else S('Undefined')
        return float(counter[0])
    B['OpenShared'] = openshared
    B['OpenShared$'] = openshared
    B['CloseShared'] = lambda a: (store.pop(int(a[0]), None), Null)[1]
    B['ReadShared'] = lambda a: store.get(int(a[0]), S('Undefined'))
    B['WriteShared'] = lambda a: (store.__setitem__(int(a[0]), a[1]), a[1])[1]
    B['Shared'] = lambda a: store.get(int(a[0]), S('Undefined'))
    B['SharedSize'] = lambda a: float(len(store))

    def parallelize(a):
        if len(a) < 4:
            raise SADError('General::narg: Number of arguments is expected to be 4 or more in Parallelize')
        f, n = a[0], int(a[1])
        comp = None
        for o in a[4:]:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'Completion':
                comp = o.rhs
        out = []
        for k in range(1, n + 1):
            try:
                out.append(I.apply_value(f, [float(k)], f))
            except SADError:
                out.append(S('$Failed'))
        if comp is not None and comp != []:
            I.apply_value(comp, [out], comp)
        return out
    B['Parallelize'] = parallelize
