# -*- coding: utf-8 -*-
"""pysad help: look up SAD/FFS commands, SADScript functions, special variables, flags and elements, and see whether pysad supports them.
pysad ヘルプ: SAD の FFS コマンド・SADScript の関数・特別変数・フラグ・要素を引き、説明の最初の行と pysad の状態を表示する。

使い方:
    import pysad
    pysad.sad_help('Twiss')                 # 名前で引く (大文字・小文字は区別しない。`Twiss[]`・`CAL` のような略形・`BEND.AE1` も可)
    pysad.sad_help(search='tune')           # キーワードで探す (名前・説明の最初の行・分類・注記)
    pysad.sad_help('FIT', full=True)        # SADHelp の全文も出す (下の「全文」)
    notebook:  %sad_help Twiss   /   sad.help('Twiss')
    コマンド:  pysad-help Twiss   /   py -m pysad help --search tune

表示するもの (同梱のデータ `pysad/helpdata/sadhelp_index.json`。docs/34・35 から `tools/build_help_data.py` で作る):
  - SADHelp (1.4.4.2k64 Updated: 10/19/2026) の説明の最初の行 (英語)
  - pysad の状態 (実装・未対応・無視 など) と注記 (日本語)、SAD 側の定義場所、詳しい文書 (docs の該当ファイルと GitHub Pages の URL)
  - pysad だけの関数 (FrequencyMap[] など) は、英語と日本語の 1 行ずつ

全文 (`full=True`): SAD 自身のヘルプの本文 `Packages/HelpMessages.n` を読んで出す。これは pysad に同梱している SAD のコピー (pysad/sad_copy) の一部で、
`PYSAD_SAD_PACKAGES` を設定していれば、その Packages の HelpMessages.n を使う。コードは評価せず、構文だけを解析して取り出す。
手元の `SADHelp.html` があれば、それも使える。探す順は、環境変数 `SADHELP_HTML` (ファイルのパス) → 同梱の HelpMessages.n →
`PYSAD_SAD_PACKAGES` の隣の `Documents/SADHelp.html` → pysad の隣の `SAD-master/Documents/SADHelp.html`。

.n の関数: 同梱の SAD の `Packages/*.n` で定義されている関数・クラス・Tkinter のウィジェットも引ける (名前・署名・定義されているファイル。コメントがあれば 1 行)。
名前に `$` を含むものと小文字で始まるものは、利用者が呼ぶものではないとみなして入れていない。説明は付いていない (署名とファイル名だけ)。

限界: 全文を同梱の HelpMessages.n から出せるのは、SADHelp 由来の 538 項目のうち 363 項目。FFS コマンド 55/56、関数 94/98 (pysad だけの 6 つの関数は HelpMessages.n にない)、
特別変数 48/48、フラグ 60/60、要素 11/15、要素のキーワード 45/155、演算子などの式 2/46、dynamics の章 1/13。取れない項目は、HelpMessages.n で実行時に作る部分に
ある (手元の SADHelp.html があれば全文が出る)。.n の関数 (487 件) と Tk のウィジェット (19 件) には、署名・定義されているファイルだけで、全文はない。
pysad の状態の「実装」は、名前の登録があるという意味で、K64 と挙動が一致するかは各要素の文書による (docs/35)。
pysad の Python の関数・クラスの説明は、Python の `help()`・Jupyter の `?` で読む。
"""
from __future__ import annotations

import html as _html
import json
import os
import re
import sys
from typing import Dict, List, Optional

SADHELP_URL = 'https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html'
PAGES_DOCS = 'https://ken-harada1234.github.io/pysad/docs/'
_INDEX_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'helpdata', 'sadhelp_index.json')
_INDEX: Optional[dict] = None
_FULL: Dict[str, Optional[dict]] = {}          # SADHelp.html のパス → 解析結果

KIND_LABEL = {
    'ffs-command': 'FFS command', 'function': 'function', 'special-variable': 'special variable', 'flag': 'flag',
    'element': 'element', 'element-keyword': 'element keyword', 'chapter': 'chapter', 'expression': 'expression',
    'beam-line': 'beam-line', 'character-string': 'character-string', 'data-structure': 'data-structure',
    'dynamics': 'dynamics', 'pattern': 'pattern', 'set-value-of-element': 'set-value-of-element',
    'package-function': 'SAD package (.n) function', 'tk-widget': 'Tk widget',
}


def index() -> dict:
    """Load the bundled help index once (a dict with `entries`).
    同梱のヘルプのデータ (辞書。`entries` を持つ) を 1 回だけ読む。"""
    global _INDEX
    if _INDEX is None:
        with open(_INDEX_PATH, encoding='utf-8') as f:
            _INDEX = json.load(f)
    return _INDEX


def entries() -> List[dict]:
    """All help entries (a list of dicts: name, keys, kind, summary, status, ...).
    すべての項目 (辞書の並び。name・keys・kind・summary・status など)。"""
    return index()['entries']


def _clean(name: str) -> str:
    s = (name or '').strip().strip('`')
    s = re.sub(r'\[\s*\]?$', '', s)             # Twiss[] → Twiss
    s = re.sub(r'\s+', ' ', s)
    return s


def lookup(name: str, kind: Optional[str] = None) -> List[dict]:
    """Return the entries named `name` (case-insensitive; abbreviations such as CAL work; `BEND.AE1` or `BEND AE1` picks the keyword of an element).
    名前 (大文字・小文字は区別しない。略形 CAL も可。`BEND.AE1` で要素のキーワード) に一致する項目の並び。`kind` で種類を絞る。"""
    s = _clean(name)
    if not s:
        return []
    parent = ''
    m = re.match(r'^([^.\s]+)[.\s]+([^.\s]+)$', s)
    if m:
        parent, s = m.group(1).upper(), m.group(2)
    key = s.upper()
    out = []
    for e in entries():
        if kind and e['kind'] != kind:
            continue
        if key in [k.upper() for k in e['keys']] or key == e['name'].upper():
            if parent and e.get('parent', '').upper() != parent:
                continue
            out.append(e)
    return out


def search(text: str, kind: Optional[str] = None) -> List[dict]:
    """Return the entries whose name, first line, category or note contains `text` (case-insensitive).
    名前・説明の最初の行・分類・注記・親の要素に `text` を含む項目 (大文字・小文字は区別しない)。"""
    t = (text or '').strip().lower()
    if not t:
        return []
    out = []
    for e in entries():
        if kind and e['kind'] != kind:
            continue
        hay = ' '.join([e['name'], ' '.join(e['keys']), e.get('summary', ''), e.get('summary_ja', ''), e.get('category', ''),
                        e.get('note', ''), e.get('parent', '')]).lower()
        if t in hay:
            out.append(e)
    return out


# ---------------------------------------------------------------------- 全文 (同梱の HelpMessages.n、手元の SADHelp.html)
def bundled_help_path() -> Optional[str]:
    """Path of HelpMessages.n (SAD's own help text in the bundled SAD copy, or in your PYSAD_SAD_PACKAGES), or None.
    SAD 自身のヘルプの本文 HelpMessages.n (同梱の SAD のコピー、または PYSAD_SAD_PACKAGES の中) のパス。無ければ None。"""
    try:
        from . import pkgload
        d = pkgload.packages_dir()
    except Exception:
        return None
    p = os.path.join(d, 'HelpMessages.n') if d else ''
    return p if p and os.path.isfile(p) else None


def find_sadhelp_html() -> Optional[str]:
    """Find the user's own SADHelp.html (env SADHELP_HTML, then <PYSAD_SAD_PACKAGES>/../Documents, then SAD-master/Documents next to pysad).
    利用者の手元の SADHelp.html のパス。見つからなければ None。"""
    cands = []
    env = os.environ.get('SADHELP_HTML', '').strip()
    if env:
        cands.append(env)
    pk = os.environ.get('PYSAD_SAD_PACKAGES', '').strip()
    if pk:
        cands.append(os.path.join(os.path.dirname(os.path.abspath(pk)), 'Documents', 'SADHelp.html'))
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cands.append(os.path.join(os.path.dirname(repo), 'SAD-master', 'Documents', 'SADHelp.html'))
    cands.append(os.path.join(repo, 'SAD-master', 'Documents', 'SADHelp.html'))
    for c in cands:
        if os.path.isfile(c):
            return c
    return None


_H3 = re.compile(r'(?is)(<UL>)|(</UL>)|<H3>\s*<A NAME=L(\d+)>\s*(.*?)</A>\s*</H3>\s*<PRE>(.*?)</PRE>')


def _parse_sadhelp(path: str) -> dict:
    """SADHelp.html を (深さ, 名前, 親, 本文) の並びに解析する。形式: `<H3><A NAME=Lnn>名前</A></H3><PRE>本文</PRE>` (入れ子は <UL>)。"""
    raw = open(path, encoding='utf-8', errors='replace').read()
    depth = 0
    last: Dict[int, str] = {}                 # 深さ → その深さで最後に出た見出し (親は 1 つ浅い深さの見出し)
    items = []
    for m in _H3.finditer(raw):
        if m.group(1):
            depth += 1
        elif m.group(2):
            depth = max(0, depth - 1)
        else:
            name = _html.unescape(re.sub(r'\s+', ' ', m.group(4))).strip()
            body = _html.unescape(m.group(5)).strip('\n')
            items.append({'depth': depth, 'name': name, 'parent': last.get(depth - 1, '') if depth else '', 'body': body})
            last[depth] = name
    return {'path': path, 'items': items}


def _hname(e) -> str:
    h = e.head
    return h if isinstance(h, str) else getattr(h, 'name', str(h))


def _conv(x, syms):
    """HelpMessages.n の構文木を、文字列と入れ子の並びにする (評価はしない)。`Sort[..]` は中身、`StringJoin` は文字列だけを連結、
    `Set[名前, 文字列]` は文字列を `syms` に覚えて返し、その名前への参照 (`IF$Message` など) は覚えた文字列に置き換える。それ以外は None。"""
    from .sadscript import Expr
    if isinstance(x, str):
        return x
    if isinstance(x, Expr):
        h = _hname(x)
        if h == 'List':
            return [_conv(a, syms) for a in x.args]
        if h == 'Sort' and x.args:
            return _conv(x.args[0], syms)
        if h == 'StringJoin':
            return ''.join(p for p in (_conv(a, syms) for a in x.args) if isinstance(p, str))
        if h == 'Set' and len(x.args) == 2 and not isinstance(x.args[0], Expr):
            v = _conv(x.args[1], syms)
            if isinstance(v, str):
                syms[getattr(x.args[0], 'name', str(x.args[0]))] = v
            return v
        return None
    return syms.get(getattr(x, 'name', str(x)))


def _parse_bundled(path: str) -> dict:
    """HelpMessages.n (`Help$Messages = {{name, text, see-also, {children...}}, ...}`) を、`_parse_sadhelp` と同じ形の項目の並びにする。
    実行時に作る部分 (Table[..]・Module[..] など) は、項目として取れないので飛ばす。"""
    from .sadscript import Expr, parse
    nodes = parse(open(path, encoding='utf-8', errors='replace').read())
    root = None
    for n in nodes:
        if isinstance(n, Expr) and _hname(n) == 'CompoundExpression' and n.args and isinstance(n.args[0], Expr) \
                and _hname(n.args[0]) == 'Set' and str(n.args[0].args[0]) == 'Help$Messages':
            rhs = n.args[0].args[1]
            root = _conv(rhs.args[-1] if isinstance(rhs, Expr) and _hname(rhs) == 'CompoundExpression' else rhs, {})
            break
    items: List[dict] = []

    def walk(node, parent, depth):
        kids = node[3] if len(node) > 3 and isinstance(node[3], list) else []
        for c in kids:
            if isinstance(c, list) and len(c) >= 2 and isinstance(c[0], str) and isinstance(c[1], str):
                body = c[1].strip('\n')
                see = c[2].strip() if len(c) > 2 and isinstance(c[2], str) else ''
                if see:
                    body += '\n\n関連 (See also): ' + see
                items.append({'depth': depth, 'name': c[0], 'parent': parent, 'body': body})
                walk(c, c[0], depth + 1)

    if root and isinstance(root, list) and isinstance(root[0], list):
        walk(root[0], '', 1)
    return {'path': path, 'items': items}


def _load(path: str, kind: str) -> Optional[dict]:
    key = kind + '|' + path
    if key not in _FULL:
        try:
            _FULL[key] = _parse_bundled(path) if kind == 'bundled' else _parse_sadhelp(path)
        except Exception:
            _FULL[key] = None
    return _FULL[key]


def _match(doc: dict, entry: dict) -> Optional[str]:
    keys = {k.upper() for k in entry['keys']} | {entry['name'].upper()}
    par = (entry.get('parent') or '').upper()
    best = None
    for it in doc['items']:
        base = it['name'].upper()
        names = {base, re.sub(r'\(.*$', '', base)}
        m = re.match(r'^([^()\s]+)\(([^)]*)\)$', it['name'])
        if m:
            names |= {m.group(1).upper(), m.group(2).upper()}
        if keys & names:
            if par and it['parent'].upper() != par:
                continue
            if best is None or (not par and not it['parent'] and best['parent']):
                best = it
                if par or not it['parent']:
                    break
    return best['body'] if best else None


def full_text_with_source(entry: dict, path: Optional[str] = None):
    """Return (full text, source label) of `entry`, or (None, '') if not found. Order: `path` (a SADHelp.html) or env SADHELP_HTML, the bundled
    HelpMessages.n, then a SADHelp.html found near the SAD folder. 項目の全文と出所のラベルを返す。見つからなければ (None, '')。"""
    sources = []
    explicit = path or (os.environ.get('SADHELP_HTML', '').strip() or None)
    if explicit and os.path.isfile(explicit):
        sources.append(('html', explicit, 'あなたの SADHelp.html'))
    bp = bundled_help_path()
    if bp:
        sources.append(('bundled', bp, 'SAD の HelpMessages.n (同梱の SAD のコピー)'))
    auto = find_sadhelp_html()
    if auto and auto != explicit:
        sources.append(('html', auto, '手元の SADHelp.html'))
    for kind, p, label in sources:
        doc = _load(p, kind)
        if doc:
            t = _match(doc, entry)
            if t is not None:
                return t, label
    return None, ''


def full_text(entry: dict, path: Optional[str] = None) -> Optional[str]:
    """Return the full SADHelp text of `entry` (from the bundled HelpMessages.n or your own SADHelp.html), or None if not found.
    項目の全文を返す (同梱の HelpMessages.n、または手元の SADHelp.html から)。見つからなければ None。"""
    return full_text_with_source(entry, path)[0]


# ---------------------------------------------------------------------- 表示
def _clipped(summary: str) -> str:
    """The first lines in the index are cut at about 125 characters; add an ellipsis when a long one does not end like a sentence.
    索引の最初の行は 125 文字ほどで切ってある。長くて文の終わりで終わっていないものには、切れたことが分かる … を付ける。"""
    return summary + ' ...' if len(summary) >= 120 and not summary.rstrip().endswith(('.', ':', ']', ')', '!', '?')) else summary


def _doc_lines(e: dict) -> List[str]:
    return ['%s  (%s%s)' % (f, PAGES_DOCS, f[:-3] + '.html') for f in e.get('docs', [])]


def format_entry(e: dict, full: bool = False, path: Optional[str] = None) -> str:
    """Format one entry as text; with full=True append the full help text (bundled HelpMessages.n or your own SADHelp.html).
    項目 1 つを文字列にする。`full=True` なら全文 (同梱の HelpMessages.n、または手元の SADHelp.html) も付ける。"""
    kind = KIND_LABEL.get(e['kind'], e['kind'])
    head = e['name'] if not e.get('parent') else '%s.%s' % (e['parent'], e['name'])
    cat = e.get('category')
    lines = ['%s  [%s%s]' % (head, kind, (' / ' + cat) if cat and cat != kind else '')]
    if e.get('summary'):
        lines.append('  SADHelp : %s' % _clipped(e['summary']))
    if e.get('summary_ja'):
        lines.append('  説明    : %s' % e['summary_ja'])
    if e.get('status'):
        s = '  pysad   : %s' % e['status']
        lines.append(s + ('  (%s)' % e['note'] if e.get('note') else ''))
    if e.get('pysad_type'):
        lines.append('  pysad   : 要素の型 %s、キーワード %s' % (e['pysad_type'], e.get('keywords', '')))
        if e.get('missing_keywords'):
            lines.append('            pysad に無いキーワード: %s' % e['missing_keywords'])
    if e.get('sad_source'):
        lines.append('  SAD 側  : %s' % e['sad_source'])
    if e.get('n_files'):
        lines.append('  .n      : %s' % ', '.join('Packages/' + f for f in e['n_files']))
    for sig in e.get('signature', []):
        lines.append('  署名    : %s' % sig)
    if e.get('tk_widget'):
        lines.append('  Tk      : %s' % e['tk_widget'])
    dl = _doc_lines(e)
    if dl:
        lines.append('  文書    : ' + dl[0])
        lines += ['            ' + x for x in dl[1:]]
    if full:
        t, src = full_text_with_source(e, path)
        if t is None:
            lines.append('  全文    : 見つかりませんでした (要素・キーワードの全文は、同梱の HelpMessages.n にはありません)。手元の SADHelp.html があれば、'
                         '環境変数 SADHELP_HTML にそのパスを設定してください。全文は %s (版が違う場合があります。項目名で探してください)' % SADHELP_URL)
        else:
            lines.append('  全文 (%s):' % src)
            lines += ['    | ' + x for x in t.split('\n')]
    return '\n'.join(lines)


def text(name: Optional[str] = None, search_for: Optional[str] = None, full: bool = False, kind: Optional[str] = None,
         path: Optional[str] = None, limit: int = 60) -> str:
    """Return the help text for a name, a keyword search (`search_for`), or (no arguments) an overview.
    名前・キーワード検索・(引数なしで) 概要のヘルプ文字列を返す。"""
    if name:
        r = lookup(name, kind)
        if not r:
            near = [e for e in search(_clean(name), kind)][:8]
            msg = ['「%s」は見つかりませんでした。' % name]
            if near:
                msg.append('名前の一部に一致する項目:')
                msg += ['  ' + _one_line(e) for e in near]
            else:
                msg.append('`search=` でキーワード検索できます (例: pysad.sad_help(search="tune"))。')
            return '\n'.join(msg)
        return '\n\n'.join(format_entry(e, full, path) for e in r)
    if search_for:
        r = search(search_for, kind)
        if not r:
            return '「%s」を含む項目は見つかりませんでした。' % search_for
        shown = r[:limit]
        out = ['「%s」を含む項目: %d 件%s' % (search_for, len(r), '' if len(r) <= limit else ' (先頭の %d 件)' % limit)]
        out += ['  ' + _one_line(e) for e in shown]
        out.append('詳しくは pysad.sad_help("名前") (全文は full=True)')
        return '\n'.join(out)
    return overview()


def _one_line(e: dict) -> str:
    nm = e['name'] if not e.get('parent') else '%s.%s' % (e['parent'], e['name'])
    summ = (e.get('summary') or e.get('summary_ja') or (e.get('signature') or [''])[0]).replace('\n', ' ')
    return '%-26s [%s] %s%s' % (nm, KIND_LABEL.get(e['kind'], e['kind']), ('(' + e['status'] + ') ') if e.get('status') else '', summ[:70])


def overview() -> str:
    """Overview of the help (kinds and counts, how to ask).
    ヘルプの概要 (種類と件数、引き方)。"""
    idx = index()
    cnt: Dict[str, int] = {}
    for e in idx['entries']:
        cnt[e['kind']] = cnt.get(e['kind'], 0) + 1
    lines = ['pysad ヘルプ (SADHelp %s の説明の最初の行 + pysad の状態)。項目 %d 件:' % (idx['sadhelp_version'], len(idx['entries']))]
    lines += ['  %-18s %4d' % (KIND_LABEL.get(k, k), n) for k, n in sorted(cnt.items(), key=lambda x: -x[1])]
    lines += ['', '引き方:',
              '  pysad.sad_help("Twiss")            名前で引く (CAL のような略形・Twiss[]・BEND.AE1 も可)',
              '  pysad.sad_help(search="tune")      キーワードで探す',
              '  pysad.sad_help("FIT", full=True)   全文も出す (同梱の HelpMessages.n。手元の SADHelp.html は SADHELP_HTML)',
              '  notebook: %sad_help Twiss / sad.help("Twiss")    コマンド: pysad-help Twiss',
              'Python の関数・クラスの説明は、help(pysad.run) や Jupyter の ? で読めます。']
    return '\n'.join(lines)


def sad_help(name: Optional[str] = None, search: Optional[str] = None, full: bool = False, kind: Optional[str] = None,
             out=None, path: Optional[str] = None) -> None:
    """Print help for a SAD/FFS/SADScript name, or search by keyword with `search=`; no argument: overview; full=True adds your own SADHelp.html text.
    SAD の名前のヘルプを表示する。`search=` でキーワード検索、引数なしで概要。`full=True` で、手元の SADHelp.html があれば全文も出す。"""
    f = out if out is not None else sys.stdout
    f.write(text(name, search, full, kind, path) + '\n')


__all__ = ['sad_help', 'text', 'lookup', 'search', 'format_entry', 'full_text', 'full_text_with_source', 'bundled_help_path',
           'find_sadhelp_html', 'entries', 'index', 'overview']
