# -*- coding: utf-8 -*-
"""全文検索の索引 (search_index.json)・検索の画面 (search.html)・入口ページ (index.html) を作る。

    ページの本文 (変換した HTML。数式は TeX の文字列、コードは文字列として) を、節 (h1〜h3 の見出しごと) に分けて集める。
    索引の形式 (search_index.json):
        {"version": 1,
         "pages":    [{"path": "docs/01_overview.html", "title": "...", "desc": "要約", "headings": ["見出し", ...]}, ...],
         "sections": [{"page": 0, "id": "見出しの ID (先頭は空)", "heading": "...", "text": "節の本文 (空白は 1 つに畳む)"}, ...]}
    検索は search.html の JavaScript (外部ライブラリなし)。語 (空白区切り) の全部を、大文字小文字・全角半角を区別せず、
    部分一致で探す (日本語・数式の TeX・コードにも効く)。

    入口ページは md の見出しから自動生成する (docs の目次の分類は docs/01_overview.md の表から読む)。
"""
import html as _html
import json
import re
import unicodedata
from html.parser import HTMLParser
from pathlib import Path

GENERATOR = 'site_index'
MAX_SECTION_CHARS = 8000       # 1 つの節の本文の上限 (表・出力の長いものを切る)
VOID = {'img', 'br', 'hr', 'input', 'meta', 'link', 'col', 'wbr', 'source', 'area', 'base', 'embed', 'param', 'track'}
BLOCK = {'p', 'li', 'tr', 'td', 'th', 'div', 'pre', 'br', 'blockquote', 'ul', 'ol', 'table', 'figure', 'figcaption',
         'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'details', 'summary'}


def indexed(rel):
    """索引と入口の目次の対象の md か (リポジトリの根からの相対パス)。"""
    parts = rel.split('/')
    if '_deletable' in parts or rel.startswith('docs/history/'):
        return False
    return (rel.startswith('docs/') or rel.startswith('example/') or rel == 'README.md' or rel == 'tools/README.md'
            or (rel.startswith('pysad/tests/') and parts[-1] == 'README.md') or rel.startswith('physics_guide/'))


def _norm(s):
    return re.sub(r'\s+', ' ', unicodedata.normalize('NFKC', s)).strip()


# ---------------------------------------------------------------- HTML から節を取り出す
class _Extractor(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.sections = []          # {'id', 'level', 'heading', 'parts': []}
        self.cur = {'id': '', 'level': 0, 'heading': '', 'parts': []}
        self.sections.append(self.cur)
        self.skip = 0               # KaTeX の中の深さ (>0 のあいだ、文字を捨てる)
        self.mathml = False
        self.annot = False
        self.head = None            # 見出しの中 (level, id, parts)
        self.blocks = []            # 段落・箇条書きの文字列 (要約用)
        self._bstack = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if self.skip:
            if tag not in VOID:
                self.skip += 1
            if tag == 'annotation':
                self.annot = True
            return
        cls = a.get('class') or ''
        if tag == 'span' and 'katex-mathml' in cls.split():
            self.skip, self.mathml = 1, True
            return
        if tag == 'span' and 'katex-html' in cls.split():
            self.skip = 1
            return
        if tag in ('script', 'style'):
            self.skip = 1
            return
        if tag in ('h1', 'h2', 'h3', 'h4', 'h5', 'h6') or (tag == 'summary' and a.get('id')):
            lv = int(tag[1]) if tag != 'summary' else 3
            self.head = (lv, a.get('id', ''), [])
        if tag in ('p', 'li'):
            self._bstack.append([])
        if tag in BLOCK:
            self._put(' ')

    def handle_endtag(self, tag):
        if self.skip:
            if tag == 'annotation':
                self.annot = False
                self._put(' ')
            if tag not in VOID:
                self.skip -= 1
                if not self.skip:
                    self.mathml = False
            return
        if self.head and tag in ('h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'summary'):
            lv, hid, parts = self.head
            self.head = None
            text = _norm(''.join(parts))
            if lv <= 3 and hid:
                self.cur = {'id': hid, 'level': lv, 'heading': text, 'parts': []}
                self.sections.append(self.cur)
            else:
                self.cur['parts'].append(' ' + text + ' ')
        if tag in ('p', 'li') and self._bstack:
            t = _norm(''.join(self._bstack.pop()))
            if t:
                self.blocks.append(t)
        if tag in BLOCK:
            self._put(' ')

    def handle_startendtag(self, tag, attrs):
        if self.skip:
            return
        if tag in BLOCK:
            self._put(' ')

    def _put(self, s):
        if self.head:
            self.head[2].append(s)
        else:
            self.cur['parts'].append(s)
        for b in self._bstack:
            b.append(s)

    def handle_data(self, data):
        if self.skip:
            if self.annot:
                self._put(' ' + data + ' ')
            return
        self._put(data)


def extract_page(body, title=''):
    """変換した本文の HTML から、索引の 1 ページ分 (dict) を作る。"""
    ex = _Extractor()
    ex.feed(body)
    ex.close()
    secs = []
    heads = []
    for s in ex.sections:
        text = _norm(''.join(s['parts']))
        if s['level'] == 0 and not text:
            continue
        if s['level'] == 1 and not title:
            title = s['heading']
        if s['level'] >= 2:
            heads.append(s['heading'])
        if len(text) > MAX_SECTION_CHARS:
            text = text[:MAX_SECTION_CHARS]
        secs.append({'id': s['id'], 'heading': s['heading'] or title, 'text': text})
    desc = ''
    for t in ex.blocks:
        if len(t) >= 12:
            desc = t
            break
    for t in ex.blocks[:12]:
        if t.startswith('この文書で分かること'):
            desc = re.sub(r'^この文書で分かること\s*[:：]?\s*', '', t)
            break
    return {'title': _norm(title), 'desc': _shorten(desc, 120), 'headings': heads, 'sections': secs}


def _shorten(t, n):
    t = _norm(t)
    if len(t) <= n:
        return t
    cut = t[:n]
    k = cut.rfind('。')
    if k >= n // 2:
        return cut[:k + 1]
    return cut.rstrip() + '…'


# ---------------------------------------------------------------- 索引の書き出し
def build_index(pages):
    """pages: {html のパス (根から): extract_page の dict} → 索引の dict。"""
    out_pages, out_secs = [], []
    for i, p in enumerate(sorted(pages)):
        d = pages[p]
        out_pages.append({'path': p, 'title': d['title'], 'desc': d['desc'], 'headings': d['headings']})
        for s in d['sections']:
            out_secs.append({'page': i, 'id': s['id'], 'heading': s['heading'], 'text': s['text']})
    return {'version': 1, 'pages': out_pages, 'sections': out_secs}


def index_to_pages(idx):
    """索引の dict を {path: page dict} に戻す (部分的な作り直しで、変換しなかったページを保つため)。"""
    pages = {}
    for i, p in enumerate(idx.get('pages', [])):
        pages[p['path']] = {'title': p['title'], 'desc': p['desc'], 'headings': p['headings'], 'sections': []}
    for s in idx.get('sections', []):
        pg = idx['pages'][s['page']]['path']
        pages[pg]['sections'].append({'id': s['id'], 'heading': s['heading'], 'text': s['text']})
    return pages


def write_index(root, pages):
    root = Path(root)
    idx = build_index(pages)
    data = json.dumps(idx, ensure_ascii=False, separators=(',', ':'))
    with open(root / 'search_index.json', 'w', encoding='utf-8', newline='\n') as f:
        f.write(data)
    return len(data.encode('utf-8')), len(idx['pages']), len(idx['sections'])


def load_index(root):
    p = Path(root) / 'search_index.json'
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None


# ---------------------------------------------------------------- md から入口ページ用の情報
def _plain_md(t):
    t = re.sub(r'!\[[^\]]*\]\([^)]*\)', '', t)
    t = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', t)
    t = re.sub(r'(\*\*|__|`)', '', t)
    t = re.sub(r'<[^>]+>', '', t)
    return _norm(t)


def md_info(path):
    """md の (h1, 最初の段落 (要約))。コードの囲み・引用・箇条書き・表は段落とみなさない。"""
    try:
        lines = Path(path).read_text(encoding='utf-8-sig').replace('\r\n', '\n').split('\n')
    except OSError:
        return '', ''
    title, i, fence = '', 0, False
    for i, l in enumerate(lines):
        if l.startswith('```'):
            fence = not fence
        if not fence and l.startswith('# '):
            title = _plain_md(l[2:])
            break
    else:
        return '', ''
    para = []
    fence = False
    for l in lines[i + 1:]:
        if l.startswith('```'):
            fence = not fence
            continue
        s = l.strip()
        if fence:
            continue
        if not s:
            if para:
                break
            continue
        if s.startswith(('#', '>', '|', '- ', '* ', '+ ', '![', '<', '---')) or re.match(r'\d+[.)] ', s):
            if para:
                break
            continue
        para.append(s)
    return title, _shorten(_plain_md(' '.join(para)), 140)


def _docs_groups(root):
    """docs/01_overview.md の目次の表 (群 | 文書 | 内容) から、(群, [(ファイル名, 内容)]) の一覧を作る。"""
    p = Path(root) / 'docs' / '01_overview.md'
    groups = []
    if not p.exists():
        return groups
    cur = None
    for l in p.read_text(encoding='utf-8-sig').split('\n'):
        m = re.match(r'^\|\s*([^|]*?)\s*\|\s*\[([^\]]+\.md)\]\(([^)]+\.md)\)\s*\|\s*([^|]*?)\s*\|\s*$', l)
        if not m:
            continue
        if m.group(1):
            cur = (m.group(1), [])
            groups.append(cur)
        if cur is not None:
            cur[1].append((m.group(3), _plain_md(m.group(4))))
    return groups


def _e(s):
    return _html.escape(s, quote=True)


INDEX_CSS = """
main { margin:0 auto; max-width:1000px; padding:24px 28px 80px; }
.bar { display:flex; gap:14px; flex-wrap:wrap; align-items:center; margin:0 0 1.2em; }
form.q { display:flex; gap:8px; margin:1em 0 1.4em; }
form.q input[type=search] { flex:1; font:inherit; padding:8px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
form.q button { font:inherit; padding:8px 16px; border:1px solid var(--line); border-radius:6px; background:var(--bg2); color:var(--fg); cursor:pointer; }
ul.toc { list-style:none; padding-left:0; margin:.4em 0; } ul.toc li { padding:3px 0; } ul.toc ul { list-style:none; padding-left:1.4em; font-size:.93em; }
.d { color:var(--muted); font-size:.9em; }
h3.g { margin:1.3em 0 .3em; font-size:1em; color:var(--muted); border:0; }
"""


def _link(href, text, root, desc='', warn=None):
    if not (Path(root) / href).exists():
        if warn is not None:
            warn.append('入口ページ: %s がない' % href)
    return '<li><a href="%s">%s</a>%s</li>' % (_e(quote_path(href)), _e(text), ' <span class="d">— %s</span>' % _e(desc) if desc else '')


def quote_path(p):
    from urllib.parse import quote
    return quote(p, safe='/#')


def build_index_page(root, rels, warn=None):
    """入口ページ index.html の HTML。rels: 対象の md (根からの相対パス) の集まり。"""
    root = Path(root)
    warn = warn if warn is not None else []
    rels = sorted(r for r in rels if indexed(r))
    html_of = lambda r: r[:-3] + '.html'
    info = {}
    for r in rels:
        info[r] = md_info(root / r)

    def title(r):
        return info[r][0] or Path(r).name

    out = []
    # 目的
    purpose = info.get('docs/01_overview.md', ('', ''))[1] or info.get('README.md', ('', ''))[1]
    out.append('<h1>pysad のドキュメント</h1>')
    if purpose:
        out.append('<p>%s</p>' % _e(purpose))
    out.append('<form class="q" action="search.html" method="get"><input type="search" name="q" placeholder="全文検索 (例: 共鳴、ハミルトニアン、FODO、Emittance)" aria-label="全文検索">'
               '<button type="submit">検索</button></form>')
    out.append('<div class="bar"><a href="search.html">全文検索</a><a href="README.html">README</a>'
               '<a href="docs/01_overview.html">文書の全体像</a></div>')
    # docs
    docs = [r for r in rels if r.startswith('docs/') and r.count('/') == 1]
    shown = set()
    out.append('<h2 id="docs">docs: 解説文書 (%d 本)</h2>' % len(docs))
    for gname, items in _docs_groups(root):
        lis = []
        for f, d in items:
            r = 'docs/' + f
            if r in info:
                shown.add(r)
                lis.append(_link(html_of(r), title(r), root, d, warn))
        if lis:
            out.append('<h3 class="g">%s</h3><ul class="toc">%s</ul>' % (_e(gname), ''.join(lis)))
    rest = [r for r in docs if r not in shown]
    if rest:
        out.append('<h3 class="g">その他</h3><ul class="toc">%s</ul>'
                   % ''.join(_link(html_of(r), title(r), root, info[r][1], warn) for r in rest))
    # example
    ex = [r for r in rels if r.startswith('example/')]
    out.append('<h2 id="example">example: 機能別の例</h2>')
    lis = []
    if 'example/README.md' in info:
        lis.append(_link('example/README.html', title('example/README.md'), root, info['example/README.md'][1], warn))
    tops = [r for r in ex if r.count('/') == 2 and r.endswith('/README.md')]
    for r in tops:
        folder = r.rsplit('/', 1)[0] + '/'
        subs = [s for s in ex if s.startswith(folder) and s != r and s.endswith('.md') and Path(s).name == 'README.md']
        li = _link(html_of(r), title(r), root, info[r][1], warn)
        if subs:
            li = li[:-5] + '<ul>%s</ul></li>' % ''.join(
                _link(html_of(s), title(s), root, '', warn) for s in subs)
        lis.append(li)
    out.append('<ul class="toc">%s</ul>' % ''.join(lis))
    # 試験・道具・物理ガイド
    out.append('<h2 id="tests">試験 (pysad/tests)</h2><ul class="toc">%s</ul>' % ''.join(
        _link(html_of(r), title(r), root, info[r][1], warn) for r in rels if r.startswith('pysad/tests/')))
    out.append('<h2 id="tools">道具 (tools)</h2><ul class="toc">%s</ul>' % ''.join(
        _link(html_of(r), title(r), root, info[r][1], warn) for r in rels if r.startswith('tools/')))
    pg = [r for r in rels if r.startswith('physics_guide/')]
    if pg:
        out.append('<h2 id="physics-guide">物理ガイド (physics_guide)</h2><ul class="toc">%s</ul>' % ''.join(
            _link(html_of(r), title(r), root, info[r][1], warn) for r in pg))
    from md_to_html import CSS
    return ('<!DOCTYPE html>\n<html lang="ja"><head><meta charset="utf-8"><meta name="generator" content="%s">'
            '<meta name="viewport" content="width=device-width,initial-scale=1"><title>pysad のドキュメント</title>'
            '<style>%s%s</style></head><body class="nonav">\n<main>\n%s\n</main></body></html>\n'
            % (GENERATOR, CSS, INDEX_CSS, '\n'.join(out)))


# ---------------------------------------------------------------- 検索の画面
SEARCH_JS = r"""
(function () {
  var box = document.getElementById('q'), out = document.getElementById('out'), stat = document.getElementById('stat');
  var scope = document.getElementById('scope');
  var data = null, norm = null;
  function esc(s) { return s.replace(/[&<>"]/g, function (c) { return {'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]; }); }
  function nfkc(s) { return s.normalize ? s.normalize('NFKC') : s; }
  function load() {
    fetch('search_index.json').then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); }).then(function (d) {
      data = d;
      norm = d.sections.map(function (s) {
        var p = d.pages[s.page];
        return { t: p.title.toLowerCase(), h: s.heading.toLowerCase(), x: s.text.toLowerCase() };
      });
      run();
    }).catch(function () {
      stat.textContent = '索引 (search_index.json) を読めません。ファイルを直接開いている (file://) ときは、リポジトリの直下で py -3.12 -m http.server を実行し、http://localhost:8000/search.html を開いてください。';
    });
  }
  function count(s, w) { var n = 0, i = -1; while ((i = s.indexOf(w, i + 1)) >= 0 && n < 50) n++; return n; }
  function snippet(text, low, words) {
    var pos = -1;
    for (var k = 0; k < words.length; k++) { var i = low.indexOf(words[k]); if (i >= 0 && (pos < 0 || i < pos)) pos = i; }
    if (pos < 0) pos = 0;
    var a = Math.max(0, pos - 40), b = Math.min(text.length, pos + 120);
    var raw = text.slice(a, b), l = low.slice(a, b);
    var marks = [];
    words.forEach(function (w) { var i = -1; while ((i = l.indexOf(w, i + 1)) >= 0) { marks.push([i, i + w.length]); i += w.length - 1; } });
    marks.sort(function (x, y) { return x[0] - y[0]; });
    var res = '', cur = 0;
    marks.forEach(function (m) { if (m[0] < cur) return; res += esc(raw.slice(cur, m[0])) + '<mark>' + esc(raw.slice(m[0], m[1])) + '</mark>'; cur = m[1]; });
    res += esc(raw.slice(cur));
    return (a > 0 ? '…' : '') + res + (b < text.length ? '…' : '');
  }
  function run() {
    if (!data) return;
    var q = nfkc(box.value).trim();
    history.replaceState(null, '', q ? '?q=' + encodeURIComponent(q) : location.pathname);
    if (!q) { out.innerHTML = ''; stat.textContent = data.pages.length + ' ページ・' + data.sections.length + ' 節を索引済み。語を空白で区切ると、すべてを含む節を探します (大文字小文字・全角半角の区別なし、部分一致)。'; return; }
    var words = q.toLowerCase().split(/\s+/).filter(Boolean), sc = scope.value, hits = [];
    for (var i = 0; i < data.sections.length; i++) {
      var s = data.sections[i], p = data.pages[s.page], n = norm[i];
      if (sc && p.path.indexOf(sc) !== 0) continue;
      var score = 0, ok = true;
      for (var k = 0; k < words.length; k++) {
        var w = words[k], c = count(n.x, w), th = n.t.indexOf(w) >= 0, hh = n.h.indexOf(w) >= 0;
        if (!c && !th && !hh) { ok = false; break; }
        score += Math.min(c, 5) + (th ? 6 : 0) + (hh ? 12 : 0);
      }
      if (ok) hits.push([score, i]);
    }
    hits.sort(function (a, b) { return b[0] - a[0] || a[1] - b[1]; });
    stat.textContent = hits.length + ' 件の節が見つかりました' + (hits.length > 60 ? ' (上位 60 件を表示)' : '') + '。';
    out.innerHTML = hits.slice(0, 60).map(function (h) {
      var s = data.sections[h[1]], p = data.pages[s.page];
      var href = p.path + (s.id ? '#' + encodeURIComponent(s.id) : '');
      return '<li><a href="' + esc(href) + '"><b>' + esc(p.title) + '</b>' + (s.id && s.heading !== p.title ? ' › ' + esc(s.heading) : '') + '</a>'
        + '<div class="sn">' + snippet(s.text, norm[h[1]].x, words) + '</div><div class="d">' + esc(p.path) + '</div></li>';
    }).join('');
  }
  var m = /[?&]q=([^&]*)/.exec(location.search);
  if (m) box.value = decodeURIComponent(m[1].replace(/\+/g, ' '));
  box.addEventListener('input', run); scope.addEventListener('change', run);
  load(); box.focus();
})();
"""

SEARCH_CSS = """
main { margin:0 auto; max-width:900px; padding:24px 28px 80px; }
.row { display:flex; gap:8px; margin:1em 0 .6em; }
.row input, .row select { font:inherit; padding:8px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
.row input { flex:1; }
#out { list-style:none; padding:0; } #out li { padding:10px 0; border-bottom:1px solid var(--line); }
#out a { text-decoration:none; } .sn { margin-top:3px; font-size:.93em; overflow-wrap:anywhere; } .d { color:var(--muted); font-size:.82em; }
mark { background:#fde68a; color:#1f2328; padding:0 1px; border-radius:2px; }
@media (prefers-color-scheme: dark) { mark { background:#6b5300; color:#fff; } }
#stat { color:var(--muted); font-size:.9em; }
"""


def build_search_page(rels=None):
    from md_to_html import CSS
    scopes = ['<option value="">すべて</option>', '<option value="docs/">docs (解説文書)</option>',
              '<option value="example/">example (例)</option>', '<option value="pysad/tests/">試験</option>',
              '<option value="tools/">tools</option>', '<option value="physics_guide/">物理ガイド</option>']
    return ('<!DOCTYPE html>\n<html lang="ja"><head><meta charset="utf-8"><meta name="generator" content="%s">'
            '<meta name="viewport" content="width=device-width,initial-scale=1"><title>全文検索</title>'
            '<style>%s%s</style></head><body class="nonav">\n<main>\n'
            '<div class="meta"><a href="index.html">入口</a></div>\n<h1>全文検索</h1>\n'
            '<div class="row"><input id="q" type="search" placeholder="例: 共鳴、ハミルトニアン、FODO、Emittance、\\nu_x" aria-label="検索語">'
            '<select id="scope" aria-label="範囲">%s</select></div>\n<div id="stat">読み込み中…</div>\n<ul id="out"></ul>\n'
            '</main><script>%s</script></body></html>\n' % (GENERATOR, CSS, SEARCH_CSS, ''.join(scopes), SEARCH_JS))


# ---------------------------------------------------------------- 書き出し・検査
def is_generated(path):
    try:
        with open(path, encoding='utf-8') as f:
            return ('<meta name="generator" content="%s">' % GENERATOR) in f.read(2000)
    except (OSError, UnicodeDecodeError):
        return False


def write_site(root, pages, rels):
    """索引・検索の画面・入口ページを書く。(メッセージのリスト, 警告のリスト) を返す。"""
    root = Path(root)
    warns, msgs = [], []
    size, npg, nsec = write_index(root, pages)
    msgs.append('search_index.json: %d ページ・%d 節・%.0f KB' % (npg, nsec, size / 1024))
    for name, html in (('search.html', build_search_page()), ('index.html', build_index_page(root, rels, warns))):
        p = root / name
        if p.exists() and not is_generated(p):
            warns.append('%s は本ツールの生成物でないので、何もしない' % name)
            continue
        with open(p, 'w', encoding='utf-8', newline='\n') as f:
            f.write(html)
        msgs.append('%s: %.1f KB' % (name, len(html.encode('utf-8')) / 1024))
    return msgs, warns


def check_site(root):
    """索引の中のパスと、入口ページ・検索の画面のリンク先の存在を確かめる。壊れているものの一覧を返す。"""
    from urllib.parse import unquote
    root = Path(root)
    broken = []
    idx = load_index(root)
    if idx is None:
        broken.append(('索引', 'search_index.json がない、または読めない'))
    else:
        for p in idx['pages']:
            if not (root / p['path']).is_file():
                broken.append(('索引のパス', p['path']))
        for s in idx['sections']:
            if not 0 <= s['page'] < len(idx['pages']):
                broken.append(('索引の節', str(s['page'])))
                break
    for name in ('index.html', 'search.html'):
        f = root / name
        if not f.is_file():
            broken.append(('入口・検索', name + ' がない'))
            continue
        for h in re.findall(r'href="([^"#]*)', f.read_text(encoding='utf-8')):
            h = _html.unescape(h)
            if h and "'" not in h and not re.match(r'^[A-Za-z][A-Za-z0-9+.\-]*:|^//', h) and not (root / unquote(h)).exists():
                broken.append(('リンク', '%s: %s' % (name, h)))
    return broken
