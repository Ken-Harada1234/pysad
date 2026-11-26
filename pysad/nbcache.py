# -*- coding: utf-8 -*-
"""Tool that discards stale numba caches (.nbi, .nbc) whose called files have changed (v2next D2).
numba のキャッシュ (`pysad/__pycache__/*.nbi`・`*.nbc`) のうち、呼ぶ先のファイルが変わって古くなったものを捨てる道具 (v2next D2)。

numba (0.68 で確認) の `cache=True` の索引 (.nbi) は、**関数を定義したファイル 1 つの内容のハッシュ**で無効になる (更新時刻ではない)。
別のファイルの njit 関数や定数 (大域の値はコンパイル時に埋め込まれる) を使うカーネルは、呼ばれる側を編集しても呼ぶ側のキャッシュが無効にならず、
古い機械語が残って結果が変わらない (D2 の実験: moda.f が modb.g を呼ぶとき、g を書き換えても f の値は古いまま。moda のファイルを書き換えると直る)。
pysad では `elkern.py` (← `perfkern.py`・`perfkern_bend.py`)、`turnkern.py`・`multa.py`・`tcav.py`・`und.py`・`beambeam.py`・`phsrot.py` (← `track6n.py`・`emit6d.py` の定数) がこれに当たる。

この道具は、`__pycache__/nbdeps.json` に「モジュールごとの、(推移的な) 依存先のハッシュ」を覚えておき、依存先のどれかが変わっていたら、
そのモジュールの .nbi・.nbc を削除する (次の実行で 1 回だけ再コンパイルされる)。依存は、各モジュールの先頭 (桁 0) の `from .X import …` だけを辿る
(関数の中の遅延 import は numba のカーネルが使わないので辿らない)。X が njit の関数を持つモジュールなら X の内容 (改行の違いは無視)、
そうでなければ import した大文字の名前 (定数) の代入行だけのハッシュ。
使い方:  `py -3.12 -m pysad.nbcache` (古いキャッシュを削除して記録を更新)、`--dry-run` は一覧だけ、`--force` は全部捨てる。
記録がない最初の 1 回は、現在の状態を記録するだけで何も削除しない (最初から古いキャッシュがあるかもしれないときは `--force`)。
kernel (`elkern.py`・`perfkern*.py`・`turnkern.py`・`track6n*.py`) を編集したら、試験の前に 1 回走らせる。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys

_PKG = os.path.dirname(os.path.abspath(__file__))
_CACHE = os.path.join(_PKG, '__pycache__')
_RECORD = os.path.join(_CACHE, 'nbdeps.json')
_BS = chr(92)
_IMP = re.compile(r'^from\s+\.(\w+)\s+import\s+(\([^)]*\)|[^\n]*)|^from\s+\.\s+import\s+([\w, ]+?)\s*(?:#.*)?$', re.M)
_CAPS = re.compile(r'^[A-Z][A-Z0-9_]+$')
_NJIT = re.compile(r'\bnjit\b|@jit\b|numba\.jit')


def _read(mod: str) -> str:
    try:
        return open(os.path.join(_PKG, mod + '.py'), encoding='utf-8', errors='replace').read().replace('\r\n', '\n')
    except OSError:
        return ''


def _kernel_module(mod: str) -> bool:
    return bool(_NJIT.search(_read(mod)))


def _imports(mod: str) -> set:
    """mod が先頭 (桁 0) で import する依存先: 'X' (njit を持つモジュール X 全体) または 'X:NAME' (X の大文字の定数 NAME)"""
    out = set()
    for m in _IMP.finditer(_read(mod)):
        if m.group(3):
            for n in m.group(3).split(','):
                n = n.split(' as ')[0].strip()
                if n and _kernel_module(n):
                    out.add(n)
            continue
        name = m.group(1)
        if not os.path.exists(os.path.join(_PKG, name + '.py')):
            continue
        if _kernel_module(name):
            out.add(name)
        else:
            for nm in re.split(r'[,\s()]+', m.group(2)):
                if _CAPS.match(nm):
                    out.add('%s:%s' % (name, nm))
    return out


def _closure(mod: str) -> set:
    seen, todo = set(), [mod]
    while todo:
        cur = todo.pop().split(':')[0]
        for d in _imports(cur):
            if d not in seen:
                seen.add(d)
                if ':' not in d:
                    todo.append(d)
    seen.discard(mod)
    return seen


def _hash(dep: str) -> str:
    if ':' in dep:
        mod, name = dep.split(':')
        lines = re.findall(r'^%s\s*=.*$' % re.escape(name), _read(mod), re.M)
        text = '\n'.join(lines)
    else:
        text = _read(dep)
    return hashlib.sha1(text.encode('utf-8')).hexdigest()


def cached_modules() -> set:
    if not os.path.isdir(_CACHE):
        return set()
    return {f.split('.')[0] for f in os.listdir(_CACHE) if f.endswith('.nbi')}


def _snapshot(mod: str) -> dict:
    return {d: _hash(d) for d in sorted(_closure(mod))}


def stale_modules(record: dict | None = None) -> dict:
    """Return {module name: [changed dependencies]} by comparing current dependency hashes with a recorded set.
    {モジュール名: [変わった依存先]}。record は前に記録したハッシュ ({モジュール: {依存先: ハッシュ}})"""
    if record is None:
        try:
            record = json.load(open(_RECORD))
        except (OSError, ValueError):
            return {}
    out = {}
    for mod in cached_modules():
        old = record.get(mod)
        if old is None:
            continue
        changed = sorted(d for d in _closure(mod) if old.get(d) != _hash(d))
        if changed:
            out[mod] = changed
    return out


def purge(dry_run: bool = False, force: bool = False) -> list:
    """Delete stale numba cache files (.nbi, .nbc) and update the record; return the deleted (or to-be-deleted) file names.
    古いキャッシュ (.nbi・.nbc) を削除し、記録を更新する。戻り値は削除した (または削除する) ファイル名"""
    try:
        record = json.load(open(_RECORD))
    except (OSError, ValueError):
        record = {}
    stale = set(cached_modules()) if force else set(stale_modules(record))
    gone = []
    if os.path.isdir(_CACHE):
        for f in sorted(os.listdir(_CACHE)):
            if f.split('.')[0] in stale and f.endswith(('.nbi', '.nbc')):
                gone.append(f)
                if not dry_run:
                    os.remove(os.path.join(_CACHE, f))
    if not dry_run:
        new = {mod: _snapshot(mod) for mod in cached_modules() | set(record) | stale}
        os.makedirs(_CACHE, exist_ok=True)
        json.dump(new, open(_RECORD, 'w'), indent=0)
    return gone


if __name__ == '__main__':
    dry = '--dry-run' in sys.argv
    force = '--force' in sys.argv
    if not force:
        for m, d in sorted(stale_modules().items()):
            print('%-14s 依存先が変わった: %s' % (m, ', '.join(d)))
    g = purge(dry, force)
    print(('削除する' if dry else '削除した'), len(g), 'ファイル')
