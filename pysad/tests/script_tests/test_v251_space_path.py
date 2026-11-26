# -*- coding: utf-8 -*-
"""空白を含むパスでの GUI (Tk) の接点: tclbridge.space_free_dir の試験。

SAD の KBFrame.n は、設定ファイルの場所を引用符なしで Tcl の `option readfile` に渡す。空白を含むパスだと Tcl が空白で切るので、
pysad は空白のないパス (8.3 の短い名前、なければ一時フォルダへのコピー) を渡す。
"""
import os
import sys

from pysad import tclbridge


def _mk(tmp_path):
    d = tmp_path / 'with space' / 'KBFrame'
    d.mkdir(parents=True)
    (d / 'KBFrame.res').write_text('*font: fixed\n')
    (d / 'sub').mkdir()
    (d / 'sub' / 'a.xbm').write_text('x')
    return str(d)


def test_path_without_space_is_returned_as_is(tmp_path):
    d = tmp_path / 'nospace'
    d.mkdir()
    assert tclbridge.space_free_dir(str(d)) == str(d).replace(chr(92), '/')


def test_path_with_space_gives_a_space_free_path_with_the_same_files(tmp_path):
    d = _mk(tmp_path)
    r = tclbridge.space_free_dir(d)
    assert not any(ch.isspace() for ch in r), r
    assert os.path.isfile(os.path.join(r, 'KBFrame.res')) and os.path.isfile(os.path.join(r, 'sub', 'a.xbm'))


def test_copy_fallback_when_short_names_are_unavailable(tmp_path, monkeypatch):
    d = _mk(tmp_path)
    monkeypatch.setattr(sys, 'platform', 'linux')           # 8.3 の短い名前を使わない経路 (一時フォルダへのコピー)
    r = tclbridge.space_free_dir(d)
    if any(ch.isspace() for ch in r):                        # 一時フォルダ自体に空白があるときは元のパスを返す
        assert r == d.replace(chr(92), '/')
    else:
        assert os.path.isfile(os.path.join(r, 'KBFrame.res'))
        assert 'pysad_kbframe_' in r


def test_option_readfile_works_with_the_returned_path(tmp_path):
    import pytest
    tk = pytest.importorskip('tkinter')
    d = _mk(tmp_path)
    r = tclbridge.space_free_dir(d)
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip('Tk を開けない')
    try:
        root.tk.eval('option readfile ' + r + '/KBFrame.res')       # KBFrame.n と同じ引用符なしの書き方
    finally:
        root.destroy()
