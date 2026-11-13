# -*- coding: utf-8 -*-
"""SAD の Tcl/Tk 接点 (src/tfTkInter_.c, src/tfTclArg.f) を Python の tkinter で実装する。

GUI のコードそのもの (Packages/Tkinter.n、KBFrame.n、KBMainFrame.n、TabFrame.n、Canvas.n ほか) は SADScript なので、翻訳せずに pysad のインタプリタで読む。
ここに書いたのは、それらが一番下で呼ぶ C / Fortran の組み込み関数だけ:

  TclCreateInterp[]   Tcl インタプリタ (tkinter.Tk) を作る。SAD の `sadeval` / `sadcommand` コマンドを Tcl に登録する
  Tcl[cmd] / Tcl[cmd, x]  Tcl を評価する (2 引数なら結果の文字列を返す)。エラーは errorInfo を出して Tkinter::tclerror
  TclSetResult[s]     コールバックから戻る Tcl の結果を設定する
  TclArg[...] / TclArg1[...]   SAD の値 (文字列、数、リスト、`オプション->値`、Widget[..]、String[..]) を Tcl の引数の文字列にする (tfTclArg.f)
  TkEventLoop[] / TkQuitEventLoop[]   イベントループ (mainloop / quit)
  そのほか TkCreateFileHandler、TkDefineBitmap、TkPhotoPutBlock$、TkCanvas* は最小限

SAD のパッケージ (`<SAD>/Packages`) の場所は環境変数 PYSAD_SAD_PACKAGES で指定する (リポジトリには含めていない)。指定が無いときは何も有効にならず、従来どおり TkWait[] などは何もしない。
環境変数 PYSAD_TK_NOLOOP=1 なら TkEventLoop[] は 1 回 update するだけで戻る (自動テスト・ウィジェット木の比較用)。
"""
from __future__ import annotations

import os
import re
import sys

from .sadscript import SADError, Expr, Rule, Sym, S, Null, is_num

_AUTOLOAD_FILES = ('TkinterCore.n', 'Tkinter.n', 'TkPhoto.n', 'Canvas.n', 'KBMainFrame.n', 'KBFFileSelectionPanel.n', 'TabFrame.n',
                   'FileDialog.n', 'OptionDialog.n', 'KBFrame.n', 'CursorEntry.n', 'CheckEntryList.n')


def packages_dir():
    p = os.environ.get('PYSAD_SAD_PACKAGES')
    if p and os.path.isdir(p):
        return p.rstrip('/\\') + '/'
    return None


class TclBridge:
    def __init__(self, interp):
        self.I = interp
        self.root = None
        self.pending_result = ''
        self.noloop = bool(os.environ.get('PYSAD_TK_NOLOOP'))
        self.packages = packages_dir()

    # ------------------------------------------------------------------ 組み込みの登録
    def install(self, B):
        B['TclCreateInterp'] = self.create_interp
        B['Tcl'] = self.tcl
        B['TclSetResult'] = self.set_result
        B['TclArg'] = lambda a: self.tclarg(a, False)
        B['TclArg1'] = lambda a: self.tclarg(a, True)
        B['TclUnicode'] = lambda a: a[0]
        B['TkEventLoop'] = self.event_loop
        B['TkQuitEventLoop'] = self.quit_loop
        B['TkCreateFileHandler'] = lambda a: Null
        B['TkDeleteFileHandler'] = lambda a: Null
        B['TkDefineBitmap'] = self.define_bitmap
        B.pop('TkWait', None)               # Tkinter.n の TkWait[] を使う (従来の何もしない版を外す)
        # ---- Xlib 関数 (src/tfXlib_.c)。KBMainFrame.n が X のウィンドウプロパティで他の SAD アプリを探す・ウィンドウを前面に出すのに使う。
        #      Windows の tkinter には X のウィンドウ木が無いので、「他のアプリは見つからない」として答える
        def _toinput(a):
            x = a[0] if a else Null
            if isinstance(x, str):
                return chr(34) + x.replace(chr(92), chr(92) * 2).replace(chr(34), chr(92) + chr(34)) + chr(34)
            return self.I.to_str(x)
        B['ToInputString'] = _toinput
        B['Environment'] = lambda a: os.environ.get(str(a[0]), '') if a and isinstance(a[0], str) else ''
        B['GetPID'] = lambda a: float(os.getpid())
        B['XQueryTree'] = lambda a: [a[0] if a else 0.0, 0.0, []]
        B['XInternAtom'] = lambda a: 1.0
        B['XGetTextProperty'] = lambda a: '***'
        B['XGetWindowProperty'] = lambda a: '***'
        B['XListFonts'] = lambda a: []
        B['XGetFontPath'] = lambda a: []
        B['XServerVendor'] = lambda a: 'pysad'
        B['XGetWindowAttributes'] = lambda a: []
        for nm in ('XFlush', 'XRaiseWindow', 'XLowerWindow', 'XRestackWindows', 'XSetTextProperty', 'XChangeProperty', 'XDeleteProperty'):
            B[nm] = lambda a: Null

    # ------------------------------------------------------------------ Tcl インタプリタ
    def create_interp(self, a):
        import tkinter
        if self.root is None:
            for k_ in range(4):                       # Windows の Tcl は、まれに起動直後の Tk の初期化 (tk.tcl の探索) に失敗する。数回やり直す
                try:
                    self.root = tkinter.Tk()
                    break
                except tkinter.TclError:
                    if k_ == 3:
                        raise
                    import time
                    time.sleep(0.3)
            tk = self.root.tk
            tk.createcommand('sadeval', self._sadeval)
            tk.createcommand('sadcommand', self._sadcommand)
            # Windows の Tk には無い (X11 の Tk の `tk::pkgconfig get fontsystem`)。SAD の TkinterCore.n が使うので X11 として答える
            tk.eval('if {[info commands tk::pkgconfig] eq {}} {proc tk::pkgconfig {args} {return x11}}')
        return 1.0

    def _result(self):
        r, self.pending_result = self.pending_result, ''
        return r

    def _sadeval(self, *argv):
        """Tcl の `sadeval code` / `sadeval f a1 a2 ...` (= f[a1,a2,...]) を SAD で評価する (tclsadeval)"""
        if len(argv) == 1:
            code = argv[0]
        elif len(argv) > 1:
            code = argv[0] + '[' + ','.join(argv[1:]) + ']'
        else:
            return ''
        code = re.sub(r'->\?\?(?=[,}])', '->"??"', code)      # Tk は該当しないイベント欄 (%X など) を ?? にする。SAD の式としては文字列として扱う
        self.pending_result = ''
        try:
            self.I.run(code, 'script')
        except SADError as e:
            self.I.out.write('???' + str(e) + chr(10))
        return self._result()

    def _sadcommand(self, *argv):
        """Tcl の `sadcommand N` → SAD の WidgetCommand[N] (tclsadcommand)"""
        if not argv:
            return ''
        self.pending_result = ''
        try:
            self.I.eval(Expr(S('WidgetCommand'), [float(int(argv[0]))]))
        except SADError as e:
            self.I.out.write('???' + str(e) + chr(10))
        return self._result()

    def tcl(self, a):
        if len(a) not in (1, 2):
            raise SADError('General::narg: Number of arguments is expected to be 1 or 2 in Tcl')
        if not isinstance(a[0], str):
            raise SADError('General::wrongtype: Argument must be Character-string in Tcl')
        if self.root is None:
            self.create_interp([])
        if os.environ.get('PYSAD_TRACE_TCL'):
            sys.stderr.write('Tcl> ' + a[0][:300] + chr(10))
        try:
            r = self.root.tk.eval(a[0])
        except Exception as e:                                     # tkinter.TclError
            try:
                self.I.out.write(str(self.root.tk.eval('set errorInfo')) + chr(10))
            except Exception:
                self.I.out.write(str(e) + chr(10))
            self.I.out.write('???Tkinter::tclerror: 1501' + chr(10))
            return Null
        if len(a) == 2:
            return str(r)
        return Null

    def set_result(self, a):
        if len(a) != 1 or not isinstance(a[0], str):
            raise SADError('General::wrongtype: Argument must be Character-string in TclSetResult')
        self.pending_result = a[0]
        return Null

    # ------------------------------------------------------------------ イベントループ
    def event_loop(self, a):
        if self.root is None:
            return Null
        if self.noloop:
            self.root.update()
            return Null
        self.root.mainloop()
        return Null

    def quit_loop(self, a):
        if self.root is not None:
            try:
                self.root.quit()
            except Exception:
                pass
        return Null

    def define_bitmap(self, a):
        # TkDefineBitmap[name, data, width, height]: ビットマップ名を Tcl に定義する
        name, data, w, h = a[0], a[1], int(a[2]), int(a[3])
        bits = ', '.join('0x%02x' % (int(x) & 255) for x in data) if isinstance(data, list) else ''
        self.tcl(['image create bitmap ' + name + ' -data {#define x_width %d\n#define x_height %d\nstatic char x_bits[] = {\n %s };\n}' % (w, h, bits)])
        return Null

    # ------------------------------------------------------------------ TclArg (tfTclArg.f)
    def _quote(self, s):
        out = []
        for ch in s:
            if ch in '"$\\[]':
                out.append('\\')
            out.append(ch)
        return '"' + ''.join(out) + '"'

    def _real(self, x):
        """実数を Tcl に渡す文字列にする。台本の $FORM に左右されず、整数なら整数 (Tk の -width は整数を要求する)、それ以外は最短の表記"""
        x = float(x)
        if x != x or x in (float('inf'), float('-inf')):
            return str(x)
        if x == int(x) and abs(x) < 1e15:
            return str(int(x))
        return repr(x)

    def _arg(self, items, single, parts):
        I = self.I
        for x in items:
            if x is Null or x is None:
                continue
            if isinstance(x, Expr) and isinstance(x.head, Sym) and x.head.name in ('Rule', 'RuleDelayed') and len(x.args) == 2:
                x = Rule(x.args[0], x.args[1], x.head.name == 'RuleDelayed')       # 保持引数の中で Expr になった a->b
            if isinstance(x, list):
                if all(is_num(y) for y in x):
                    parts.append(' '.join(self._real(y) for y in x))
                else:
                    sub = []
                    self._arg(x, single, sub)
                    parts.append(' '.join(sub))
                continue
            if isinstance(x, Rule):
                lab = I.eval(Expr(S('TkOptionLabel'), [x.lhs]))
                if not isinstance(lab, str):
                    raise SADError('Tkinter::badoption: bad option label')
                if lab.strip() in ('-scriptfont', '-altfont', '-altscriptfont'):
                    continue                  # SAD が手を入れた Tk だけにある上付き・下付き用のフォント指定 (標準の Tk には無い)
                rhs = I.eval(x.rhs) if getattr(x, 'delayed', False) is False else x.rhs
                parts.append(lab)
                if isinstance(rhs, list):
                    sub = []
                    self._arg(rhs, single, sub)
                    parts.append('{' + ' '.join(sub) + '}')
                else:
                    self._arg([rhs], single, parts)
                continue
            if isinstance(x, Expr) and isinstance(x.head, Sym) and x.head.name == 'Widget':
                p = I.eval(Expr(S('TkPathName'), list(x.args[:1])))        # Widget[w, type] のパス名は第 1 引数 w から (tfTclArg.f)
                self._arg([p], single, parts)
                continue
            if isinstance(x, Expr) and isinstance(x.head, Sym) and x.head.name == 'String':
                for y in x.args:
                    parts.append(y if isinstance(y, str) else I.to_str(y))
                continue
            if isinstance(x, str):
                if x:
                    parts.append(x if single else self._quote(x))
                elif not single:
                    parts.append('{}')
                continue
            if is_num(x):
                parts.append(self._real(x))
                continue
            parts.append(I.to_str(x))

    def tclarg(self, a, single):
        parts = []
        self._arg(a, single, parts)
        s = ' '.join(parts)
        return ('{' + s + '} ') if single else (s + ' ')

    # ------------------------------------------------------------------ SAD のパッケージの自動読み込み (init.n の AutoLoad)
    def setup_autoload(self):
        """<SAD>/Packages/init.n から Tk / KBFrame 関係の AutoLoad[…] だけを実行して、最初に使われたときに .n を読み込む"""
        I = self.I
        pk = self.packages
        if not pk:
            return False
        I.assign('PackagesDirectory', pk)
        share = os.path.join(os.path.dirname(pk.rstrip('/' + chr(92))), 'share', 'KBFrame').replace(chr(92), '/')
        if not os.path.isdir(share):                 # ソース配布 (SAD-master) では <SAD>/KBFrame に置かれている
            share = os.path.join(os.path.dirname(pk.rstrip('/' + chr(92))), 'KBFrame').replace(chr(92), '/')
        if os.path.isdir(share):
            os.environ.setdefault('KBFRAMEDIR', share)
            I.assign('KBF$ConfDir', share)           # KBFrame の設定ファイル・画像の置き場 (SAD の share/KBFrame)
        # SAD (X11) の既定フォント。init.n が C 側から受け取る値 (WSL の gs で実測)
        for k_, v_ in (('$DefaultFont', ['helvetica neue', 12.0]), ('$DefaultSansSerifFont', 'helvetica neue'), ('$DefaultSerifFont', 'utopia'),
                       ('$DefaultFixedFont', 'courier'), ('$DefaultTerminalFont', 'fixed'), ('$DefaultSymbolFont', 'symbol'),
                       ('$DefaultFontSize', 12.0), ('$DefaultFontSlant', 'roman'), ('$DefaultFontWeight', 'medium'), ('System$Names', [])):
            if not I.has(k_):
                I.assign(k_, v_)
        I.assign('ReturnToSAD', 0.0)                 # SAD の既定 (False)。Tkinter.n が読み込み時に 'Return to SAD' の窓を作らないように
        I.load_package(pk + 'PkgConfig.n')          # init.n の `Get[PackagesDirectory//"PkgConfig.n"]`
        src = open(os.path.join(pk, 'init.n'), encoding='utf-8', errors='replace').read()
        for st in I.split_statements(src):
            t = st.strip()
            if re.match(r'AutoLoad(Set)?\[', t) and any(f in t for f in _AUTOLOAD_FILES) and 'Ca' not in re.sub(r'Canvas|Cache|Capt', '', t.split('PackagesDirectory')[0]).split('[', 1)[-1][:0]:
                try:
                    I.run(t + ';', 'script')
                except SADError as e:
                    I.out.write('[pysad] autoload: ' + str(e) + chr(10))
        return True


def install(I):
    """Interp に Tcl/Tk の接点を付ける (PYSAD_SAD_PACKAGES が指定されているときだけ)"""
    if packages_dir() is None:
        return None
    if sys.getrecursionlimit() < 30000:
        sys.setrecursionlimit(30000)          # KBFrame.n の入れ子の評価は深い (既定の 1000 では足りない)
    br = TclBridge(I)
    br.install(I.builtins)
    br.setup_autoload()
    I.tclbridge = br
    return br
