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
        """Store (or replace) a lattice text under `name`. 同じ名前なら置き換える (セルを再実行しても重ならない)。"""
        self.lattices[name or 'main'] = text

    def clear_lattices(self) -> None:
        self.lattices.clear()

    def reset(self) -> None:
        """Forget the stored lattices and the last run. 保存したラティスと、最後の実行の状態を捨てる。"""
        self.clear_lattices()
        self.ffs = None
        self.interp = None

    # ------------------------------------------------------------------ 実行
    def run(self, src: str, cont: bool = False, out=None) -> None:
        """Run `src`.
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
    """Run SAD text from Python (a string). 文字列の SAD のテキストを実行する (マジックなしの書き方)。
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
    """`%%sad` (セル) と `%sad` (1 行) の本体。"""
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


def load_ipython_extension(ip) -> None:
    """`%load_ext pysad.notebook`。`%%sad` と `%sad` を登録し、名前空間に `sad` (Session) を入れる。"""
    ip.register_magic_function(sad_magic, magic_kind='line_cell', magic_name='sad')
    ip.user_ns.setdefault('sad', session)


def unload_ipython_extension(ip) -> None:
    try:
        del ip.magics_manager.magics['cell']['sad']
        del ip.magics_manager.magics['line']['sad']
    except KeyError:
        pass
