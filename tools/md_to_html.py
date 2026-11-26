# -*- coding: utf-8 -*-
"""Markdown を、図入り・目次つきの 1 枚の HTML にする (markdown パッケージだけを使う。外部通信なし)。

CLI:
    py -3.12 tools/md_to_html.py README.md                 # README.html を同じフォルダに作る
    py -3.12 tools/md_to_html.py a.md -o out.html --fold-regex "^(ラウンド|反省)" --fold-level 3

出力は md と同じ名前の .html (同じフォルダ)。本ツールが作った .html は上書きし、
そうでない既存の .html は、何もせず警告する。

機能:
- 左固定の目次 (見出し 2〜3 階層)、表 (見出し行を固定・クリックで並べ替え)、コードブロック、明暗の配色。
- 画像記法 ![説明](相対パス) の画像 (png/jpg/gif/svg/webp) を base64 で埋め込む。
  1 枚が MAX_EMBED (3 MB) を超えるときは、埋め込まず相対リンクのままにして警告する。
- 本文中にバッククォートで書かれた図のファイル名 (`fig_xxx.png` など。同じフォルダまたは results/ にあるもの) の
  最初の言及 (コードブロックの外) を含むブロックの直後に、その図を置く。
- 文書間リンク [..](xxx.md#anchor) は、変換対象の md なら .html に書き換える。
- 見出しのアンカーは GitHub 風 (日本語も可。重複は -1, -2)。リンクの # は実在の ID に合わせる。
- 数式 ($...$、$$...$$) は、Node.js の KaTeX で HTML の部品にして埋め込む (mathtex.py。開くときに JavaScript は要らない)。
- コードブロック (python・bash・sh・json・sad、言語指定なしの SAD の台本) は、ビルドのときに Pygments で色分けする (codehl.py。開くときに JavaScript は要らない)。
- root を渡すと、各ページの先頭の注記の行に、検索 (search.html) と入口 (index.html) へのリンクを置く。
- --fold-regex に一致する見出し (--fold-level の階層) の節を <details> にする。既定は折りたたまない。
"""
import argparse
import base64
import html as _html
import mimetypes
import os
import re
import sys
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import quote, unquote

import markdown
import markdown.extensions.toc as _toc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mathtex  # noqa: E402
import codehl  # noqa: E402
import site_index  # noqa: E402

MAX_EMBED = 3 * 1024 * 1024
GENERATOR = 'md_to_html'
IMG_EXT = ('.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp')
MIME = {'.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg', '.gif': 'image/gif',
        '.svg': 'image/svg+xml', '.webp': 'image/webp'}
_warned_node = []      # Node.js が見つからない警告は 1 回だけ出す
_SCHEME = re.compile(r'^(?:[A-Za-z][A-Za-z0-9+.\-]*:|//)')

CSS = """
:root { --fg:#1f2328; --muted:#59636e; --line:#d1d9e0; --bg:#ffffff; --bg2:#f6f8fa; --acc:#0969da; }
@media (prefers-color-scheme: dark) { :root { --fg:#e6edf3; --muted:#9198a1; --line:#3d444d; --bg:#0d1117; --bg2:#151b23; --acc:#4493f8; } }
* { box-sizing: border-box; }
body { margin:0; color:var(--fg); background:var(--bg); font:15px/1.7 "Yu Gothic UI","Meiryo","Hiragino Sans",system-ui,sans-serif; }
nav { position:fixed; top:0; left:0; bottom:0; width:300px; overflow:auto; padding:16px 14px; border-right:1px solid var(--line); background:var(--bg2); font-size:13px; }
nav ul { list-style:none; margin:0; padding-left:12px; } nav > ul, nav .toc > ul { padding-left:0; }
nav a { color:var(--fg); text-decoration:none; display:block; padding:2px 4px; border-radius:4px; } nav a:hover { background:var(--line); }
nav .ttl { font-weight:700; margin-bottom:8px; color:var(--muted); }
main { margin-left:300px; padding:24px 40px 80px; max-width:1180px; }
body.nonav main { margin-left:0; margin-right:auto; }
a { color:var(--acc); }
h1 { font-size:1.7em; border-bottom:1px solid var(--line); padding-bottom:.3em; }
h2 { font-size:1.35em; border-bottom:1px solid var(--line); padding-bottom:.25em; margin-top:2em; }
h3 { font-size:1.12em; margin-top:1.6em; }
code { background:var(--bg2); padding:.15em .35em; border-radius:4px; font-size:.9em; font-family:Consolas,"Cascadia Mono",monospace; }
pre { background:var(--bg2); padding:12px 14px; border-radius:6px; overflow:auto; border:1px solid var(--line); }
pre code { background:none; padding:0; }
blockquote { margin:1em 0; padding:.6em 1em; border-left:4px solid #d29922; background:var(--bg2); color:var(--fg); }
table { border-collapse:collapse; margin:1em 0; display:block; overflow:auto; max-width:100%; font-size:13.5px; }
th, td { border:1px solid var(--line); padding:5px 10px; } th { background:var(--bg2); position:sticky; top:0; cursor:pointer; white-space:nowrap; }
th:hover { color:var(--acc); } tr:nth-child(even) td { background:var(--bg2); }
figure { margin:1.2em 0; padding:10px; border:1px solid var(--line); border-radius:6px; background:#fff; }
figure img { max-width:100%; height:auto; display:block; margin:0 auto; }
figcaption { color:#59636e; font-size:12.5px; margin-top:6px; text-align:center; }
details { margin:1em 0; border:1px solid var(--line); border-radius:6px; padding:0 14px; }
details > summary { cursor:pointer; font-weight:700; padding:10px 0; font-size:1.05em; }
details[open] > summary { border-bottom:1px solid var(--line); margin-bottom:8px; }
.meta { color:var(--muted); font-size:12.5px; margin-bottom:1.5em; } .meta a { color:var(--acc); }
@media (max-width: 900px) { nav { position:static; width:auto; border-right:0; } main { margin-left:0; padding:16px; } }
"""
JS = """
document.querySelectorAll('table').forEach(function (t) {
  t.querySelectorAll('th').forEach(function (th, i) {
    th.title = 'クリックで並べ替え';
    th.addEventListener('click', function () {
      var tb = t.querySelector('tbody'); if (!tb) return;
      var rows = Array.from(tb.querySelectorAll('tr')); var asc = th.dataset.asc !== '1'; th.dataset.asc = asc ? '1' : '0';
      var num = function (s) { var v = parseFloat(String(s).replace(/[,×倍%]/g, '')); return isNaN(v) ? null : v; };
      rows.sort(function (a, b) {
        var x = (a.children[i] || {textContent: ''}).textContent.trim(), y = (b.children[i] || {textContent: ''}).textContent.trim();
        var nx = num(x), ny = num(y);
        var c = (nx !== null && ny !== null) ? nx - ny : x.localeCompare(y, 'ja'); return asc ? c : -c; });
      rows.forEach(function (r) { tb.appendChild(r); });
    });
  });
});
function openTarget() {
  var id = decodeURIComponent(location.hash.slice(1)); if (!id) return;
  var el = document.getElementById(id); if (!el) return;
  for (var p = el.parentElement; p; p = p.parentElement) { if (p.tagName === 'DETAILS') p.open = true; }
  el.scrollIntoView();
}
window.addEventListener('hashchange', openTarget); openTarget();
"""


# ---------------------------------------------------------------- 見出しの ID
_SLUG_MATH = None        # (items, prefix): 見出しの数式の置き換え名 (prefix+番号+END) を、TeX の文字に戻して ID を作る (番号は文書内の式の通し番号で、編集でずれるため)


def gh_slug(value, separator='-'):
    """GitHub 風: 小文字にし、英数字・日本語・_・- と空白以外を除き、空白を - にする。
    見出しの数式は、置き換え名ではなく TeX の文字から作る (`$H_2$ の厳密解` は `h_2-の厳密解`。式の数が変わっても ID は変わらない)。"""
    if _SLUG_MATH and _SLUG_MATH[0]:
        value = mathtex.plain(value, _SLUG_MATH[0], _SLUG_MATH[1])
    v = _html.unescape(value).strip().lower()
    v = re.sub(r'[^\w\- ]', '', v)
    return v.replace(' ', '-')


def _gh_unique(id, ids):
    base = id or 'section'
    n = 0
    while not id or id in ids:
        n += 1
        id = '%s-%d' % (base, n)
    ids.add(id)
    return id


@contextmanager
def _gh_ids():
    old = _toc.unique
    _toc.unique = _gh_unique
    try:
        yield
    finally:
        _toc.unique = old


def _new_md():
    return markdown.Markdown(
        extensions=['tables', 'fenced_code', 'toc', 'sane_lists'],
        extension_configs={'toc': {'toc_depth': '2-3', 'permalink': False, 'slugify': gh_slug}})


# ---------------------------------------------------------------- 前処理
_MARK = re.compile(r'^(\s*)([-*+]|\d+[.)])\s+\S')
_FENCE = re.compile(r'^\s*(```+|~~~+)')


def normalize(text):
    """GitHub では通るが Python-Markdown では崩れる書き方を直す (コードブロックの中は触らない)。
    - 段落の直後に空行なしで始まる箇条書き・表に、空行を入れる。
    - 2 スペースなどの入れ子の箇条書きを、階層ごとに 4 スペースにそろえる。"""
    out = []
    fence = None
    stack = []          # 箇条書きの各階層の元のインデント
    prev = ''
    for line in text.split('\n'):
        fm = _FENCE.match(line)
        if fence is not None:
            if fm and fm.group(1)[0] == fence and line.strip().strip(fence) == '':
                fence = None
            out.append(line)
            prev = line
            continue
        if fm:
            fence = fm.group(1)[0]
            out.append(line)
            prev = line
            continue
        s = line.strip()
        if not s:
            out.append(line)
            prev = line
            continue
        mm = _MARK.match(line)
        indented = line.startswith((' ', '\t'))
        if mm:
            ind = len(mm.group(1).expandtabs(4))
            if not stack and prev.strip() and not prev.lstrip().startswith('|') and not prev.startswith((' ', '\t')):
                out.append('')
            while stack and stack[-1] >= ind:
                stack.pop()
            depth = len(stack)
            stack.append(ind)
            line = ' ' * (4 * depth) + line.lstrip()
        elif s.startswith('|') and prev.strip() and not prev.lstrip().startswith('|'):
            out.append('')
        elif not indented:
            stack = []
        out.append(line)
        prev = line
    return '\n'.join(out)


def _ids_of_text(text):
    md = _new_md()
    text = text.replace('\r\n', '\n')
    global _SLUG_MATH
    if mathtex.find_node():
        mtext, items_, prefix_ = mathtex.extract(text)
    else:
        mtext, items_, prefix_ = text, [], ''
    _SLUG_MATH = (items_, prefix_) if items_ else None
    try:
        with _gh_ids():
            body = md.convert(normalize(mtext))
    finally:
        _SLUG_MATH = None
    ids = set(re.findall(r'<h[1-6] id="([^"]+)"', body))
    return ids | set(re.findall(r'<a\s[^>]*?(?:id|name)="([^"]+)"', text))


_ID_CACHE = {}


def ids_of_md(path):
    path = Path(path)
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _ID_CACHE:
        _ID_CACHE[key] = _ids_of_text(path.read_text(encoding='utf-8-sig'))
    return _ID_CACHE[key]


# ---------------------------------------------------------------- HTML の加工
_BLOCK_TAG = re.compile(r'<(/?)(ul|ol|table|blockquote|pre|div|details|p|h[1-6])(?=[\s>/])[^>]*>')


def _top_blocks(h):
    """最上位のブロック要素の (開始, 終了) の一覧。"""
    spans = []
    depth = 0
    start = 0
    for m in _BLOCK_TAG.finditer(h):
        if m.group(1) == '':
            if depth == 0:
                start = m.start()
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                spans.append((start, m.end()))
            elif depth < 0:
                depth = 0
    return spans


def _embed_src(src, srcdir, embed, res):
    """<img> の src を処理する。埋め込めたら data URI を返す。そうでなければ None。"""
    path = (srcdir / unquote(src)).resolve()
    if not path.is_file():
        res.broken.append(('画像', src))
        return None
    if path.suffix.lower() not in MIME:
        return None
    size = path.stat().st_size
    if size > MAX_EMBED:
        res.warnings.append('画像が %.1f MB で 3 MB を超えるので埋め込まない: %s' % (size / 1048576, src))
        return None
    if not embed:
        res.figures += 1
        return None
    data = base64.b64encode(path.read_bytes()).decode('ascii')
    res.figures += 1
    return 'data:%s;base64,%s' % (MIME[path.suffix.lower()], data)


def place_named_figures(h, srcdir):
    """`fig_xxx.png` のように書かれた図の名前の、最初の言及 (コードブロックの外) を含むブロックの直後に図を置く。"""
    used = set()
    for m in re.finditer(r'<img\s[^>]*?src="([^"]*)"', h):
        s = _html.unescape(m.group(1))
        if not _SCHEME.match(s):
            used.add((srcdir / unquote(s)).resolve())
    pres = [(m.start(), m.end()) for m in re.finditer(r'<pre\b.*?</pre>', h, re.S)]
    blocks = _top_blocks(h)
    seen = set()
    inserts = []
    pat = r'<code>([^<>]+\.(?:%s))</code>' % '|'.join(e[1:] for e in IMG_EXT)
    for m in re.finditer(pat, h, re.I):
        pos = m.start()
        if any(a <= pos < b for a, b in pres):
            continue
        name = _html.unescape(m.group(1)).strip()
        if name in seen:
            continue
        found = None
        for cand in (srcdir / name, srcdir / 'results' / name):
            try:
                if cand.is_file():
                    found = cand.resolve()
                    break
            except OSError:
                pass
        if found is None:
            continue
        seen.add(name)
        if found in used:
            continue
        used.add(found)
        end = next((b for a, b in blocks if a <= pos < b), None)
        if end is None:
            continue
        rel = found.relative_to(srcdir.resolve()).as_posix() if srcdir.resolve() in found.parents else name
        fig = ('\n<figure><img alt="%s" src="%s"><figcaption>%s</figcaption></figure>'
               % (_html.escape(name), _html.escape(quote(rel)), _html.escape(name)))
        inserts.append((end, len(inserts), fig))
    for end, _, fig in sorted(inserts, reverse=True):
        h = h[:end] + fig + h[end:]
    return h


def embed_images(h, srcdir, embed, res):
    def sub(m):
        tag = m.group(0)
        sm = re.search(r'\ssrc="([^"]*)"', tag)
        if not sm:
            return tag
        src = _html.unescape(sm.group(1))
        if _SCHEME.match(src) or not src:
            return tag
        data = _embed_src(src, srcdir, embed, res)
        if data is None:
            return tag
        return tag[:sm.start(1)] + data + tag[sm.end(1):]
    h = re.sub(r'<img\s[^>]*>', sub, h)

    def fig(m):
        img = m.group(1)
        am = re.search(r'\salt="([^"]*)"', img)
        cap = am.group(1) if am else ''
        return '<figure>%s%s</figure>' % (img, '<figcaption>%s</figcaption>' % cap if cap else '')
    return re.sub(r'<p>(<img\s[^>]*>)</p>', fig, h)


def rewrite_links(h, src, is_target, anchor_lookup, res, own_ids):
    srcdir = src.parent

    def fix_anchor(frag, ids, shown):
        if not frag:
            return frag
        a = unquote(frag)
        if ids is None or a in ids:
            return frag
        low = {i.lower(): i for i in ids}
        if a.lower() in low:
            return quote(low[a.lower()])
        res.broken.append(('アンカー', shown))
        return frag

    def sub(m):
        href = _html.unescape(m.group(2))
        if not href or _SCHEME.match(href):
            return m.group(0)
        path, _, frag = href.partition('#')
        path = path.partition('?')[0]
        if not path:
            new = '#' + fix_anchor(frag, own_ids, href)
            return m.group(1) + _html.escape(new) + m.group(3)
        target = (srcdir / unquote(path)).resolve()
        if not target.exists():
            res.broken.append(('リンク', href))
            return m.group(0)
        if target.suffix.lower() != '.md' or not target.is_file():
            return m.group(0)
        ids = anchor_lookup(target) if frag else None
        frag = fix_anchor(frag, ids, href)
        if is_target(target):
            path = path[:-3] + '.html'
        new = path + ('#' + frag if frag else '')
        return m.group(1) + _html.escape(new, quote=False).replace('"', '&quot;') + m.group(3)
    return re.sub(r'(<a\s[^>]*?href=")([^"]*)(")', sub, h)


_HEAD = re.compile(r'<h([1-6]) id="([^"]*)">(.*?)</h\1>', re.S)


def fold_sections(h, fold_re, level):
    """見出しの文字列が fold_re に一致する level 階層の節を <details> にする (同じか上の階層の見出しまで)。"""
    rx = re.compile(fold_re)
    out = []
    pos = 0
    stack = []   # 開いている折りたたみの階層
    for m in _HEAD.finditer(h):
        lv = int(m.group(1))
        out.append(h[pos:m.start()])
        pos = m.end()
        while stack and stack[-1] >= lv:
            stack.pop()
            out.append('</details>')
        text = _html.unescape(re.sub(r'<[^>]+>', '', m.group(3))).strip()
        if lv == level and rx.search(text):
            out.append('<details><summary id="%s">%s</summary>' % (m.group(2), m.group(3)))
            stack.append(lv)
        else:
            out.append(m.group(0))
    out.append(h[pos:])
    out.append('</details>' * len(stack))
    return ''.join(out)


# ---------------------------------------------------------------- 変換
class Result:
    def __init__(self, src):
        self.src = src
        self.html = ''
        self.title = ''
        self.figures = 0
        self.warnings = []
        self.broken = []     # (種類, リンク先)
        self.out = None
        self.skipped = False
        self.math_count = 0      # 検出した数式の数
        self.math_errors = 0     # KaTeX で変換できなかった数式の数
        self.code_blocks = 0     # 色分けしたコードブロックの数
        self.index = None        # 全文検索の索引用の、ページの内容 (site_index.extract_page)


def _flat(toks):
    for t in toks:
        yield t
        yield from _flat(t.get('children', []))


def convert_text(text, src, is_target=None, anchor_lookup=None, fold_regex=None, fold_level=3,
                 embed=True, math=True, root=None):
    """md の本文 (str) を HTML (Result.html) にする。src は元の md のパス (相対リンクの基準)。"""
    src = Path(src).resolve()
    res = Result(src)
    is_target = is_target or (lambda p: Path(p).is_file())
    anchor_lookup = anchor_lookup or ids_of_md
    text = text.replace('\r\n', '\n')
    md = _new_md()
    # 数式: markdown の変換の前にプレースホルダーに置き換え (強調などに誤変換されないように)、変換後に戻す
    items, rendered, prefix, mtext = [], [], '', text
    if math:
        node = mathtex.find_node()
        if node is None:
            if '$' in text and not _warned_node:
                _warned_node.append(1)
                res.warnings.append('Node.js が見つからないので、数式は TeX の文字のまま (PYSAD_NODE、PATH、'
                                    'D:\\work\\tools\\node-*\\node.exe の順に探す)')
        else:
            mtext, items, prefix = mathtex.extract(text)
            if items:
                try:
                    rendered = mathtex.render(items, node)
                except mathtex.MathError as e:
                    res.warnings.append('数式を変換できないので、TeX の文字のまま: %s' % e)
                    mtext, items, rendered, prefix = text, [], [], ''
    res.math_count = len(items)
    for it, (h_, err) in zip(items, rendered):
        if err:
            res.math_errors += 1
            res.warnings.append('数式を変換できない (行 %d): %s  [%s]' % (it['line'] + 1, it['tex'][:70], err[:120]))
    global _SLUG_MATH
    _SLUG_MATH = (items, prefix) if items else None
    try:
        with _gh_ids():
            body = md.convert(normalize(mtext))
    finally:
        _SLUG_MATH = None
    tm = re.search(r'<h1 id="[^"]*">(.*?)</h1>', body, re.S)
    title_src = mathtex.plain(tm.group(1), items, prefix) if tm and items else (tm.group(1) if tm else '')
    res.title = _html.unescape(re.sub(r'<[^>]+>', '', title_src)).strip() if tm else src.name
    body = mathtex.restore(body, items, rendered, prefix)
    body, res.code_blocks = codehl.highlight_body(body)
    own_ids = set(re.findall(r'<h[1-6] id="([^"]+)"', body)) | set(re.findall(r'<a\s[^>]*?(?:id|name)="([^"]+)"', text))
    body = place_named_figures(body, src.parent)
    body = embed_images(body, src.parent, embed, res)
    body = rewrite_links(body, src, is_target, anchor_lookup, res, own_ids)
    res.index = site_index.extract_page(body, res.title)
    if fold_regex:
        body = fold_sections(body, fold_regex, fold_level)
    toc = mathtex.restore(md.toc, items, rendered, prefix)
    has_nav = any(t['level'] >= 2 for t in _flat(md.toc_tokens))   # h2 以下がなければ目次なし
    nav = ('<nav><div class="ttl">目次</div>%s</nav>\n' % toc) if has_nav else ''
    css = CSS
    if any(h_ and 'class="katex' in h_ for h_, _ in rendered):     # 数式のある HTML にだけ、使うフォントを埋めた CSS を入れる
        css += mathtex.css_for(''.join(h_ for h_, _ in rendered if h_))
    if res.code_blocks:
        css += codehl.CSS
    links = ''
    if root is not None:
        up = Path(os.path.relpath(Path(root).resolve(), src.parent)).as_posix()
        up = '' if up == '.' else up + '/'
        links = ' <span class="nav-links">| <a href="%ssearch.html">検索</a> | <a href="%sindex.html">入口</a></span>' % (up, up)
    page = ('<!DOCTYPE html>\n<html lang="ja"><head><meta charset="utf-8">'
            '<meta name="generator" content="%s">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>%s</title><style>%s</style></head><body%s>\n%s<main>\n'
            '<div class="meta">元: %s%s</div>\n%s\n</main><script>%s</script></body></html>\n'
            % (GENERATOR, _html.escape(res.title), css, '' if has_nav else ' class="nonav"', nav,
               _html.escape(src.name), links, body, JS))
    res.html = page
    return res


def is_generated(path):
    try:
        with open(path, encoding='utf-8') as f:
            return ('<meta name="generator" content="%s">' % GENERATOR) in f.read(2000)
    except (OSError, UnicodeDecodeError):
        return False


def convert_file(src, out=None, write=True, **opts):
    """md のファイルを変換して、同じ名前の .html に書く。md 以外の既存の .html は上書きしない。"""
    src = Path(src)
    text = src.read_text(encoding='utf-8-sig')
    res = convert_text(text, src, **opts)
    out = Path(out) if out else src.with_suffix('.html')
    res.out = out
    if not write:
        return res
    if out.exists() and not is_generated(out):
        res.skipped = True
        res.warnings.append('同名の既存の .html は本ツールの生成物でないので、何もしない: %s' % out)
        return res
    with open(out, 'w', encoding='utf-8', newline='\n') as f:
        f.write(res.html)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('md')
    ap.add_argument('-o', '--output', help='出力先 (既定: md と同じ名前の .html)')
    ap.add_argument('--fold-regex', help='この正規表現に一致する見出しの節を <details> にする')
    ap.add_argument('--fold-level', type=int, default=3, help='折りたたむ見出しの階層 (既定 3)')
    ap.add_argument('--require-math', action='store_true', help='Node.js が見つからない、または数式の変換エラーがあるとき、エラー終了にする')
    a = ap.parse_args(argv)
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding='utf-8', errors='replace')
        except (AttributeError, ValueError):
            pass
    if a.require_math and mathtex.find_node() is None:
        print('エラー: Node.js が見つからない (PYSAD_NODE、PATH の node、D:/work/tools/node-*/node.exe の順に探した)')
        return 2
    res = convert_file(a.md, a.output, fold_regex=a.fold_regex, fold_level=a.fold_level)
    for w in res.warnings:
        print('警告:', w)
    for k, t in res.broken:
        print('壊れた%s: %s' % (k, t))
    if not res.skipped:
        print('作成:', res.out, '%.1f KB' % (len(res.html.encode('utf-8')) / 1024), '図 %d 枚' % res.figures)
    return 1 if (res.skipped or (a.require_math and res.math_errors)) else 0


if __name__ == '__main__':
    sys.exit(main())
