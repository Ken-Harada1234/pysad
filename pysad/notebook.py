# -*- coding: utf-8 -*-
"""Run SAD lattice and command text pasted into a Jupyter / IPython notebook (the `%%sad` cell magic).
Jupyter / IPython のセルに貼った SAD のラティスとコマンドを実行する (`%%sad` セルのマジック)。

使い方 (Jupyter・VS Code・JupyterLab の notebook):

    %load_ext pysad.notebook          # 最初に 1 回。`sad` (Session) が名前空間に入る

    %%sad lattice                     # ラティスファイルの中身をそのまま貼る (実行せず保存。同じ名前のセルを再実行すると置き換わる)
    ...

    %%sad                             # コマンドファイルの中身をそのまま貼る。保存したラティス + このセルを、最初から実行する
    FFS; USE RING; CELL; CALC; ...    # (SAD の `gs ラティス コマンド` と同じ)

    %%sad --continue                  # 前の実行の状態の続きに、FFS の文を流す (FIT のあとに Twiss、など)

マジックを使わずに Python から書くとき (`r'''...'''` なら引用符とバックスラッシュをそのまま書ける):

    import pysad
    pysad.run(r'''...ラティスとコマンド...''')

実行のあとの状態は `sad.ffs` (FFS)・`sad.interp` (インタプリタ)。SAD の式の値は `sad.value('Twiss["NX","$$$"]')` で Python の値として取れる。

注意: `--continue` で流せるのは FFS 部の文 (USE・CALC・FIT・GO・Twiss[] など)。要素・LINE の定義は、新しい `%%sad` のセル (最初からの実行) に書く。
`stop;` のあとでも `--continue` は使える (停止の印を解除する)。
IPython の仕様で、セルのマジックの本文は変換されずに渡る (`!` で始まる行の SAD のコメントも、そのまま SAD に届く)。
"""
from __future__ import annotations

import re
import shlex
import sys
from typing import Dict, Optional


class Session:
    """Holds the stored lattice texts and the state of the last run.
    保存したラティスの文字列と、最後の実行の状態 (FFS・インタプリタ) を持つ。"""

    def __init__(self) -> None:
        self.lattices: Dict[str, str] = {}      # 名前 → ラティスの文字列 (挿入順に連結する。同じ名前は置き換える)
        self.ffs = None
        self.interp = None

    def __repr__(self) -> str:
        return '<pysad.notebook.Session lattices=%s, ran=%s>' % (list(self.lattices), self.interp is not None)

    # ------------------------------------------------------------------ ラティス
    def add_lattice(self, text: str, name: str = 'main') -> None:
        """Store (or replace) a lattice text under `name`.
        ラティスの文字列を名前 `name` で保存する。同じ名前なら置き換える (セルを再実行しても重ならない)。"""
        self.lattices[name or 'main'] = text

    def clear_lattices(self) -> None:
        """Forget the stored lattices (the last run is kept).
        保存したラティスだけを捨てる (最後の実行の状態は残す)。"""
        self.lattices.clear()

    def reset(self) -> None:
        """Forget the stored lattices and the last run.
        保存したラティスと、最後の実行の状態を捨てる。"""
        self.clear_lattices()
        self.ffs = None
        self.interp = None

    # ------------------------------------------------------------------ 実行
    def run(self, src: str, cont: bool = False, out=None) -> None:
        """Run `src`: from scratch (stored lattices + src, like `gs lattice command`), or with cont=True in the state of the previous run (FFS statements only; RuntimeError if there is none).
        `src` を実行する。
        cont=False: 保存したラティス + src を連結して、最初から実行する (gs の `ラティス コマンド` と同じ。FFS の文の前に `FFS;` が要る)。
        cont=True: 前の実行の状態に、src の FFS の文を流す。まだ実行していなければ RuntimeError。"""
        out = out if out is not None else sys.stdout
        if cont:
            if self.interp is None:
                raise RuntimeError('--continue: 続きにする実行がまだありません。先に `%%sad` のセルを (--continue なしで) 実行してください')
            self._resume_state(out)
            try:
                self.interp.run(src, 'ffs')
            except SystemExit:
                pass
            return
        from .ffscmd import run_sad_string
        text = ''.join(t.rstrip('\n') + '\n' for t in self.lattices.values()) + src
        ffs, interp = run_sad_string(text, out=out)
        self.ffs, self.interp = ffs, interp

    def _resume_state(self, out) -> None:
        """出力先を今のセルに向け、`stop;` の印を解除する。"""
        I, F = self.interp, self.ffs
        I.out = out
        if F is not None:
            F.out = out
            if hasattr(F, 'stopped'):
                F.stopped = False
        if hasattr(I, 'stopped'):
            I.stopped = False

    def help(self, name: Optional[str] = None, search: Optional[str] = None, full: bool = False, kind: Optional[str] = None, out=None) -> None:
        """Show help for a SAD/FFS/SADScript name (`sad.help('Twiss')`) or search by keyword (`search='tune'`); no argument: overview.
        SAD の名前のヘルプを出す (`sad.help('Twiss')`)。`search=` でキーワード検索、引数なしで概要。`full=True` で手元の SADHelp.html の全文も出す。
        中身は `pysad.sad_help` (pysad/sadhelp.py)。"""
        from .sadhelp import sad_help
        sad_help(name, search, full, kind, out)

    def value(self, expr: str):
        """Evaluate a SAD expression in the last run and return a Python value (float, list, or None for Null).
        最後の実行の状態で SAD の式を評価し、Python の値を返す (数は float、並びは list、Null は None)。"""
        if self.interp is None:
            raise RuntimeError('まだ SAD を実行していません')
        import io
        saved_out, saved_echo = self.interp.out, getattr(self.interp, 'echo_out', False)
        self._resume_state(io.StringIO())
        self.interp.echo_out = False
        try:
            r = self.interp.run(expr.rstrip().rstrip(';') + ';', 'ffs')
        finally:
            self.interp.out = saved_out
            self.interp.echo_out = saved_echo
            if self.ffs is not None:
                self.ffs.out = saved_out
        if r.__class__.__name__ == '_Null':
            return None
        return r


session = Session()


def run(src: str, lattice: Optional[str] = None, cont: bool = False, reset: bool = False, out=None, sess: Optional[Session] = None) -> Session:
    """Run SAD text from Python (a string), without the magic: optional lattice text, cont (continue the previous state), reset (forget stored lattices first).
    文字列の SAD のテキストを実行する (マジックなしの書き方)。
    lattice: ラティスの文字列 (任意。無ければ src にラティスも書く)、cont: 前の状態の続き、reset: 保存を捨ててから。"""
    s = sess or session
    if reset:
        s.reset()
    if lattice:
        s.add_lattice(lattice)
    s.run(src, cont=cont, out=out)
    return s


# ---------------------------------------------------------------------- IPython のマジック
_USAGE = ('使い方: `%%sad` (保存したラティス + このセルを最初から実行) / `%%sad lattice [名前]` (ラティスを保存。実行しない) / '
          '`%%sad --continue` (前の状態の続きに FFS の文を流す) / `%%sad --reset` (保存を捨ててから実行)')


def _parse_line(line: str):
    """`%%sad` の行から (モード, 名前, cont, reset) を取り出す。"""
    cont = reset = False
    words = []
    for t in shlex.split(line or ''):
        if t in ('--continue', '-c'):
            cont = True
        elif t in ('--reset', '-r'):
            reset = True
        elif t.startswith('-'):
            raise ValueError('不明なオプション %s。%s' % (t, _USAGE))
        else:
            words.append(t)
    if words and words[0].lower() == 'lattice':
        if len(words) > 2:
            raise ValueError(_USAGE)
        if cont:
            raise ValueError('`lattice` と --continue は同時に使えません。%s' % _USAGE)
        return 'lattice', (words[1] if len(words) > 1 else 'main'), cont, reset
    if words:
        raise ValueError('不明な指定 %r。%s' % (words[0], _USAGE))
    return 'run', 'main', cont, reset


def sad_magic(line: str, cell: Optional[str] = None):
    """Run SAD lattice/command text pasted into a cell. セルに貼った SAD のラティスとコマンドを実行する。

    %%sad lattice [name]   save the cell as a lattice (not run; the same name replaces). ラティスを保存する (実行しない。同じ名前は置き換える)
    %%sad                  run the stored lattices + this cell from scratch (like `gs lattice command`). 保存したラティス + このセルを最初から実行する
    %%sad --continue (-c)  run FFS statements in the state of the previous run. 前の実行の状態の続きに FFS の文を流す
    %%sad --reset (-r)     forget stored lattices and the last run (then run this cell). 保存したラティスと実行の状態を捨てる (そのあと、このセルを実行)
    %sad --reset           forget them without running anything. 実行せずに捨てる

    The first line of the cell is the magic; the rest is passed to SAD as it is (`!` comments, quotes and backslashes are not changed).
    After a run, `sad.value('Twiss["NX","$$$"]')` returns a Python value; `sad.help('Twiss')` / `%sad_help Twiss` shows help for SAD names.
    セルの 1 行目がマジックで、2 行目以降はそのまま SAD に渡る。実行のあと `sad.value('式')` で Python の値を取り出せる。
    SAD の名前の説明は `sad.help('Twiss')` か `%sad_help Twiss`。FFS 部の文だけが --continue で流せる (要素・LINE の定義は新しい `%%sad` のセルに書く)。
    """
    from IPython.core.error import UsageError
    try:
        mode, name, cont, reset = _parse_line(line)
    except ValueError as e:
        raise UsageError(str(e))
    if cell is None:                    # 行のマジック `%sad --reset` だけで、保存を捨てる
        if reset:
            session.reset()
            print('(pysad) 保存したラティスと実行の状態を捨てました')
            return None
        raise UsageError('`%sad` は `--reset` だけに使えます。セルは `%%sad` で書きます。' + _USAGE)
    if reset:
        session.reset()
    if mode == 'lattice':
        session.add_lattice(cell, name)
        print('(pysad) ラティス %r を保存しました (%d 行。保存中: %s)。コマンドは `%%%%sad` のセルに書きます' % (
            name, cell.count('\n') + 1, ', '.join(session.lattices)))
        return None
    try:
        session.run(cell, cont=cont)
    except RuntimeError as e:
        raise UsageError(str(e))
    return None


def sad_help_magic(line: str):
    """Show help for SAD/FFS/SADScript names: `%sad_help Twiss`, `%sad_help --search tune`, `%sad_help FIT --full`.
    SAD の FFS コマンド・関数・特別変数・フラグ・要素の説明を出す。`%sad_help Twiss`、`%sad_help --search tune`、`%sad_help FIT --full`。

    Options: --search/-s KEYWORD (search names, first lines and notes), --kind/-k KIND (ffs-command, function, special-variable, flag, element, ...),
    --full/-f (also print the full SADHelp text if your own SADHelp.html is found; set SADHELP_HTML to its path).
    オプション: --search/-s キーワード、--kind/-k 種類 (ffs-command・function・special-variable・flag・element など)、
    --full/-f (手元の SADHelp.html があれば全文も出す。環境変数 SADHELP_HTML にそのパス)。引数なしで概要。
    """
    from IPython.core.error import UsageError
    from .sadhelp import sad_help
    srch = None
    full = False
    kind = None
    toks = shlex.split(line or '')
    k = 0
    words = []
    while k < len(toks):
        t = toks[k]
        if t in ('--full', '-f'):
            full = True
        elif t in ('--search', '-s'):
            k += 1
            if k >= len(toks):
                raise UsageError('--search にはキーワードが要ります')
            srch = toks[k]
        elif t in ('--kind', '-k'):
            k += 1
            if k >= len(toks):
                raise UsageError('--kind には種類が要ります')
            kind = toks[k]
        elif t.startswith('-'):
            raise UsageError('不明なオプション %s。`%%sad_help?` で使い方を見られます' % t)
        else:
            words.append(t)
        k += 1
    name = None
    if words:
        name = words[0] if len(words) == 1 else ' '.join(words)       # `BEND AE1` (要素 キーワード) も引ける
    sad_help(name, srch, full, kind)
    return None


_NAMES = None


def _sad_names():
    """タブ補完に使う SAD の名前 (FFS コマンド・関数・特別変数・フラグ・要素)。"""
    global _NAMES
    if _NAMES is None:
        from .sadhelp import entries
        names = set()
        for e in entries():
            if e['kind'] in ('ffs-command', 'function', 'special-variable', 'flag', 'element'):
                n = re.sub(r'\(.*$', '', e['name'])
                if re.match(r'^[A-Za-z$][A-Za-z0-9$]*$', n):
                    names.add(n)
        _NAMES = sorted(names, key=str.lower)
    return _NAMES


def _complete_sad_help(self, event):
    """`%sad_help <Tab>` の補完 (IPython の complete_command フック)。"""
    tok = (event.symbol or '').lower()
    return [n for n in _sad_names() if n.lower().startswith(tok)]


def _make_cell_matcher():
    """`%%sad` のセルの中で、SAD の名前をタブ補完するマッチャー (IPython 8.6 以降の context_matcher)。使えない IPython では None。"""
    try:
        from IPython.core.completer import SimpleCompletion, context_matcher
    except ImportError:
        return None

    @context_matcher(identifier='pysad.sad_names')
    def matcher(ctx):
        none = {'completions': [], 'suppress': False}
        if not ctx.full_text.lstrip().startswith('%%sad') or getattr(ctx, 'cursor_line', 1) == 0:
            return none                      # %%sad のセルの 2 行目以降だけ (1 行目はマジックの指定)
        m = re.search(r'[A-Za-z$][A-Za-z0-9$]*$', ctx.line_with_cursor[:ctx.cursor_position])     # cursor_position はその行の中の位置
        if not m:
            return none
        tok = m.group(0).lower()
        return {'completions': [SimpleCompletion(text=n, type='SAD') for n in _sad_names() if n.lower().startswith(tok)],
                'suppress': False}
    return matcher


_CELL_MATCHER = None


def load_ipython_extension(ip) -> None:
    """`%load_ext pysad.notebook`。`%%sad`・`%sad`・`%sad_help` を登録し、名前空間に `sad` (Session) を入れる。SAD の名前のタブ補完も入れる。
    Register `%%sad`, `%sad` and `%sad_help`, put `sad` (Session) into the namespace, and add tab completion of SAD names."""
    global _CELL_MATCHER
    ip.register_magic_function(sad_magic, magic_kind='line_cell', magic_name='sad')
    ip.register_magic_function(sad_help_magic, magic_kind='line', magic_name='sad_help')
    ip.user_ns.setdefault('sad', session)
    try:
        ip.set_hook('complete_command', _complete_sad_help, str_key='%sad_help')
    except Exception:
        pass
    if _CELL_MATCHER is None:
        _CELL_MATCHER = _make_cell_matcher()
    comp = getattr(ip, 'Completer', None)
    if _CELL_MATCHER is not None and comp is not None and _CELL_MATCHER not in comp.custom_matchers:
        comp.custom_matchers.append(_CELL_MATCHER)


def unload_ipython_extension(ip) -> None:
    """`%unload_ext pysad.notebook`: remove the `%%sad`, `%sad` and `%sad_help` magics and the name completer.
    `%unload_ext pysad.notebook`: マジックと名前の補完を外す。"""
    for kind, name in (('cell', 'sad'), ('line', 'sad'), ('line', 'sad_help')):
        try:
            del ip.magics_manager.magics[kind][name]
        except KeyError:
            pass
    comp = getattr(ip, 'Completer', None)
    if _CELL_MATCHER is not None and comp is not None and _CELL_MATCHER in comp.custom_matchers:
        comp.custom_matchers.remove(_CELL_MATCHER)
