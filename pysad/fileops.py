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
        from . import compat
        x, y, n = a

        def pos1(v):
            # 負の位置: K64 (gs) は len+v+1 ("abc" の -1 は "abXc")。v2 (B197) は Mathematica の末尾から数える len+v+2 ("abcX")
            v = int(v)
            return v if v >= 0 else len(x) + v + (2 if compat.v2() else 1)
        if isinstance(n, list):
            pos = sorted(pos1(v) for v in n)
            out = ''
            prev = 0
            for p_ in pos:
                out += x[prev:p_ - 1] + y
                prev = p_ - 1
            return out + x[prev:]
        k = pos1(n)
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
    """shared.n の Parallelize (逐次実行) と、共有メモリ OpenShared / Shared / ReadShared / WriteShared / SharedSize、
    プロセスの関数 Fork / Wait / Wait4 / Kill / Sleep (K64 の tfProcess_.c、tfrbuf.f、tfshared.f)。

    - 共有メモリは無名の mmap (POSIX では fork した子と共有される)。値は符号化して書く。OpenShared[n] の n はバイト数で、
      K64 と同じく SharedSize[値] が n を超える値は Shared::toolarge。
    - Fork[] は POSIX (Linux・WSL・macOS) では本当にプロセスを分ける (親には子の pid、子には 0)。子は Exit[] か最上位の文の終わりで
      os._exit する (呼び出し元の Python に戻らない)。Windows には fork が無いので、v2 モードでは分かるエラーにし
      (k64 モードは従来どおり「親の側」として pid を返す)、並列の計算は Parallelize (逐次) を使う。
    - Parallelize[f, n, lshare, npara, opts] は {f[1], …, f[n]} (失敗は $Failed)。POSIX で npara > 1・n > 1 なら fork した子で npara 個まで同時に計算する
      (環境変数 PYSAD_PARALLELIZE=serial で逐次)。Windows (fork が無い) と npara = 1 は逐次。Completion->g があれば結果に g を適用する。"""
    import mmap
    import pickle
    import struct
    import time
    import sys
    from .sadscript import SADError, Rule, Expr, is_num, _exit_child
    from . import compat
    B = I.builtins
    store = {}                       # 番号 → (mmap, n バイト)
    counter = [1024]
    UNDEF = S('Undefined')
    T_NULL, T_SYM, T_RULE, T_EXPR = '#N', '#S', '#R', '#E'

    def enc(v):
        if isinstance(v, (float, int, complex, str)):
            return v
        if isinstance(v, list):
            return [enc(x) for x in v]
        if v is Null:
            return (T_NULL,)
        if isinstance(v, Sym):
            return (T_SYM, v.name)
        if isinstance(v, Rule):
            return (T_RULE, enc(v.lhs), enc(v.rhs), bool(v.delayed))
        if isinstance(v, Expr):
            return (T_EXPR, enc(v.head), [enc(x) for x in v.args])
        raise SADError('General::wrongtype: Shared: value must be numbers, strings, symbols, lists or expressions')

    def dec(v):
        if isinstance(v, list):
            return [dec(x) for x in v]
        if isinstance(v, tuple):
            if v[0] == T_NULL:
                return Null
            if v[0] == T_SYM:
                return S(v[1])
            if v[0] == T_RULE:
                return Rule(dec(v[1]), dec(v[2]), v[3])
            if v[0] == T_EXPR:
                return Expr(dec(v[1]), [dec(x) for x in v[2]])
        return v

    def nwords(v):
        """tfshared.f tfsharedsize: 実数・記号は 1 語、文字列は 3+(len+7)/8 語、リストは 3+len 語 + 中の実体 (それぞれ -1)"""
        if isinstance(v, str):
            return 3 + (len(v.encode('utf-8', 'replace')) + 7) // 8
        if isinstance(v, (list, Expr, Rule)):
            items = v if isinstance(v, list) else ([v.lhs, v.rhs] if isinstance(v, Rule) else v.args)
            n = 3 + len(items)
            for x in items:
                if isinstance(x, (str, list, Expr, Rule)):
                    n += nwords(x) - 1
            return n
        return 1

    def sharedsize(a):
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in SharedSize')
        return float(max(0, nwords(a[0]) - 2) * 8)
    B['SharedSize'] = sharedsize

    def openshared(a):
        if len(a) != 1 or not isinstance(a[0], (int, float)) or isinstance(a[0], bool):
            raise SADError('General::wrongtype: Argument must be Real in OpenShared')
        n = int(a[0])
        if n < 0:
            return S('$Failed')
        counter[0] += 1
        mm = mmap.mmap(-1, 8 + 64 * (n + 8) + 4096)        # 符号化 (pickle) の分の余裕。上限の判定は SharedSize で K64 と同じにする
        mm[0:8] = struct.pack('<q', -1)                       # -1: まだ書かれていない (Undefined)
        store[counter[0]] = (mm, n)
        return float(counter[0])
    B['OpenShared'] = openshared
    B['OpenShared$'] = openshared

    def _slot(x, name):
        if not isinstance(x, (int, float)) or int(x) not in store:
            raise SADError('Shared::notopen: Shared memory is not open in %s' % name)
        return store[int(x)]

    def readshared(a):
        mm, _n = _slot(a[0], 'ReadShared')
        ln = struct.unpack('<q', mm[0:8])[0]
        if ln < 0:
            return UNDEF
        return dec(pickle.loads(mm[8:8 + ln]))

    def writeshared(a):
        mm, n = _slot(a[0], 'WriteShared')
        v = a[1]
        if isinstance(v, (str, list, Expr, Rule)):
            over = max(0, nwords(v) - 2) * 8 - n
            if over > 0:
                raise SADError('Shared::toolarge: Value is larger than the shared area by %d bytes (give a larger size to OpenShared)' % over)
        data = pickle.dumps(enc(v), protocol=4)
        if 8 + len(data) > len(mm):
            raise SADError('Shared::toolarge: Value is larger than the shared area by %d bytes (give a larger size to OpenShared)' % (8 + len(data) - len(mm)))
        mm[0:8] = struct.pack('<q', -1)                       # 書いている間は未定義 (読む側が半端な値を見ない)
        mm[8:8 + len(data)] = data
        mm[0:8] = struct.pack('<q', len(data))
        return v

    def closeshared(a):
        ent = store.pop(int(a[0]), None) if a and isinstance(a[0], (int, float)) else None
        if ent is not None:
            try:
                ent[0].close()
            except (BufferError, ValueError):
                pass
        return Null
    B['CloseShared'] = closeshared
    B['ReadShared'] = readshared
    B['WriteShared'] = writeshared
    B['Shared'] = lambda a: readshared(a) if len(a) == 1 and isinstance(a[0], (int, float)) else Expr(S('Shared'), list(a))
    _close0 = B.get('Close')

    def close(a):
        if len(a) == 1 and isinstance(a[0], (int, float)) and int(a[0]) in store:
            return closeshared(a)
        return _close0(a) if _close0 is not None else Null
    B['Close'] = close

    # ---- プロセス (tfProcess_.c)
    children = set()
    I._forked_child = False

    def fork(a):
        if a and not (len(a) == 1 and a[0] is Null):
            raise SADError('General::narg: Number of arguments is expected to be 0 in Fork')
        if not hasattr(os, 'fork'):
            if compat.v2():
                raise SADError('System::error: Fork[] is not available on this platform (Windows has no fork). '
                               'Use Parallelize[f, n, lshare, npara] (runs serially in pysad), or run pysad on Linux/WSL')
            return float(os.getpid())            # k64 モード: 従来どおり (子は作らず、親の側として pid を返す)
        import sys as _sys
        for st in (I.out, _sys.stdout, _sys.stderr):
            try:
                st.flush()
            except Exception:
                pass
        pid = os.fork()
        if pid == 0:
            I._forked_child = True               # 子: Exit[] と最上位の文の終わりで os._exit (sadscript.Interp.run)
            children.clear()
            return 0.0
        children.add(pid)
        return float(pid)
    B['Fork'] = fork

    def _nochild():
        if compat.v2():
            raise SADError('System::error: no child process (Fork[] is not available on this platform)')

    def wait(a):
        if a and not (len(a) == 1 and a[0] is Null):
            return Expr(S('Wait'), list(a))       # Wait[mem] (CSR.n) などの別の定義のための形
        if not hasattr(os, 'fork'):
            _nochild()
            return [0.0, 0.0]
        try:
            pid, status = os.wait()
        except ChildProcessError:
            raise SADError('System::error: No child processes in Wait[]')
        children.discard(pid)
        return [float(pid), float(status)]
    B['Wait'] = wait

    def wait4(a):
        if not a or (len(a) == 1 and a[0] is Null):
            wpid = -1
            flags = []
        elif isinstance(a[0], (int, float)):
            wpid = int(a[0])
            flags = a[1:]
        else:
            raise SADError('General::wrongtype: Argument must be wpid for #1 in Wait4')
        if not hasattr(os, 'fork'):
            _nochild()
            return [0.0, 0.0]
        opt0 = 0
        for f in flags:
            if not isinstance(f, (int, float)):
                raise SADError('General::wrongtype: Argument must be Wait flag for ##2 in Wait4')
            opt0 |= int(f)
        opt = (os.WNOHANG if opt0 & 1 else 0) | (os.WUNTRACED if opt0 & 2 else 0) | (getattr(os, 'WCONTINUED', 0) if opt0 & 4 else 0)
        try:
            pid, status, ru = os.wait4(wpid, opt)
        except ChildProcessError:
            return [0.0, 0.0]
        except InterruptedError:
            return [-1.0, 0.0]
        if pid == 0:
            return [0.0, 0.0]
        children.discard(pid)
        b = lambda x: 1.0 if x else 0.0
        out = [float(pid), float(os.WEXITSTATUS(status)) if os.WIFEXITED(status) else 0.0]
        if hasattr(os, 'WIFCONTINUED'):
            out.append(Rule(S('Continued'), b(os.WIFCONTINUED(status))))
        out += [Rule(S('Exited'), b(os.WIFEXITED(status))), Rule(S('Signaled'), b(os.WIFSIGNALED(status))),
                Rule(S('Stopped'), b(os.WIFSTOPPED(status)))]
        if os.WIFSIGNALED(status):
            out += [Rule(S('Signal'), float(os.WTERMSIG(status))), Rule(S('Core'), b(os.WCOREDUMP(status)))]
        if os.WIFSTOPPED(status):
            out.append(Rule(S('Signal'), float(os.WSTOPSIG(status))))
        out += [Rule(S('UserTime'), float(ru.ru_utime)), Rule(S('SystemTime'), float(ru.ru_utime))]      # K64 は SystemTime にも utime を入れる
        return out
    B['Wait4'] = wait4

    def kill(a):
        if len(a) not in (1, 2):
            raise SADError('General::narg: Number of arguments is expected to be 1 or 2 in Kill')
        if not isinstance(a[0], (int, float)):
            raise SADError('General::wrongtype: Argument must be Process ID for #1 in Kill')
        import signal
        sig = int(a[1]) if len(a) == 2 else signal.SIGTERM
        pid = int(a[0])
        if pid not in children:
            # 安全のため、この pysad が Fork[] で作った子だけに信号を送る (K64 は任意の pid に送れる)
            raise SADError('System::error: Kill[pid]: pid %d is not a child created by Fork[] in this session' % pid)
        try:
            os.kill(pid, sig)
        except OSError as e:
            raise SADError('System::error: Kill: %s' % e)
        return Null
    B['Kill'] = kill

    def sleep(a):
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in Sleep')
        if not isinstance(a[0], (int, float)):
            raise SADError('General::wrongtype: Argument must be Real number (seconds) in Sleep')
        if a[0] > 0:
            time.sleep(float(a[0]))
        return 1.0
    _sleep0 = B.get('Sleep')
    if _sleep0 is None:
        B['Sleep'] = sleep
    _pause0 = B.get('Pause')

    def pause(a):
        # Pause は K64 では Sleep と同じ。pysad は従来どおり待たない (描画の間などに使われる) が、Fork[] の子がいる間は本当に待つ
        if children and len(a) == 1 and isinstance(a[0], (int, float)):
            return sleep(a)
        return _pause0(a) if _pause0 is not None else Null
    B['Pause'] = pause
    for nm_, v_ in (('WNOHANG', 1.0), ('WUNTRACED', 2.0), ('WCONTINUED', 4.0)):
        I.globals.setdefault(nm_, v_)

    def _parallel_fork(f, n, lshare, npara, tmo):
        """POSIX: f[1..n] を fork した子で同時に npara 個まで計算する (shared.n の Parallelize と同じ分担。結果は子が一時ファイルに書き、親が読む)。
        子は親の複製なので、f が書き換える変数は親に戻らない (K64 と同じ)。子の Print は親と同じ出力先に出る。
        結果が lshare バイト (SharedSize) を超えると、子が Shared::toolarge を出し、その結果は $Failed。Timeout 秒を超えた子は殺して $Failed"""
        import tempfile
        d = tempfile.mkdtemp(prefix='pysad_par_')
        for st in (getattr(I, 'out', None), sys.stdout, sys.stderr):
            try:
                st.flush()
            except Exception:
                pass
        running = {}                                     # pid -> (k, 開始時刻)
        status = {}
        results = [S('$Failed')] * n

        def reap(block):
            while running:
                try:
                    pid, st = os.waitpid(-1, 0 if block else os.WNOHANG)
                except ChildProcessError:
                    running.clear()
                    return
                if pid == 0:
                    return
                if pid in running:
                    k0, _t0 = running.pop(pid)
                    status[k0] = st
                    return
        try:
            for k in range(1, n + 1):
                while len(running) >= npara:
                    now = time.time()
                    for pid, (k0, t0) in list(running.items()):
                        if now - t0 > tmo:
                            try:
                                os.kill(pid, 9)
                            except OSError:
                                pass
                    reap(False)
                    if len(running) >= npara:
                        time.sleep(0.005)
                pid = os.fork()
                if pid == 0:
                    I._forked_child = True
                    code = 0
                    try:
                        try:
                            v = I.apply_value(f, [float(k)], f)
                        except SADError as e:
                            I.write_unit(6, '???' + str(e) + chr(10))
                            v = S('$Failed')
                        if isinstance(v, (str, list, Expr, Rule)) and max(0, nwords(v) - 2) * 8 > lshare:
                            I.write_unit(6, '???Shared::toolarge:  Value is larger than the shared area by %d bytes (Parallelize lshare = %d)'
                                         % (max(0, nwords(v) - 2) * 8 - lshare, lshare) + chr(10))
                            v = S('$Failed')
                        with open(os.path.join(d, '%d' % k), 'wb') as fh:
                            pickle.dump(enc(v), fh, protocol=4)
                    except BaseException:
                        code = 1
                    _exit_child(I, code)
                running[pid] = (k, time.time())
            deadline_check = time.time()
            while running:
                now = time.time()
                for pid, (k0, t0) in list(running.items()):
                    if now - t0 > tmo:
                        try:
                            os.kill(pid, 9)
                        except OSError:
                            pass
                reap(False)
                if running:
                    time.sleep(0.005)
            for k in range(1, n + 1):
                p = os.path.join(d, '%d' % k)
                if os.path.exists(p):
                    try:
                        with open(p, 'rb') as fh:
                            results[k - 1] = dec(pickle.load(fh))
                    except Exception:
                        results[k - 1] = S('$Failed')
            return list(results)
        finally:
            import shutil
            shutil.rmtree(d, ignore_errors=True)

    def parallelize(a):
        if len(a) < 4:
            raise SADError('General::narg: Number of arguments is expected to be 4 or more in Parallelize')
        f, n = a[0], int(a[1])
        comp = None
        for o in a[4:]:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'Completion':
                comp = o.rhs
        npara = int(a[3]) if is_num(a[3]) and not isinstance(a[3], complex) else 1
        lshare = int(a[2]) if is_num(a[2]) and not isinstance(a[2], complex) else 0
        tmo = 2000.0
        for o in a[4:]:
            if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == 'Timeout' and is_num(o.rhs):
                tmo = float(o.rhs)
        out = None
        if (hasattr(os, 'fork') and not getattr(I, '_forked_child', False) and npara > 1 and n > 1
                and os.environ.get('PYSAD_PARALLELIZE', '').lower() != 'serial'):
            out = _parallel_fork(f, n, lshare, npara, tmo)
        if out is None:
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
