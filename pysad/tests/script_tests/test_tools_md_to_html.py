# -*- coding: utf-8 -*-
"""tools/md_to_html.py と tools/build_html.py の試験 (一時フォルダの小さな md と画像)。"""
import base64
import re
import sys
from pathlib import Path

import pytest

markdown = pytest.importorskip('markdown')
TOOLS = Path(__file__).resolve().parents[3] / 'tools'
sys.path.insert(0, str(TOOLS))
import md_to_html as m2h  # noqa: E402
import build_html  # noqa: E402

PNG = base64.b64decode(
    'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==')
SVG = b'<svg xmlns="http://www.w3.org/2000/svg" width="4" height="4"><rect width="4" height="4"/></svg>'


@pytest.fixture
def work(tmp_path):
    (tmp_path / 'results').mkdir()
    (tmp_path / 'img.png').write_bytes(PNG)
    (tmp_path / 'pic.svg').write_bytes(SVG)
    (tmp_path / 'results' / 'fig_a.png').write_bytes(PNG)
    (tmp_path / 'results' / 'fig_b.png').write_bytes(PNG)
    (tmp_path / 'big.png').write_bytes(PNG + b'\0' * (m2h.MAX_EMBED + 10))
    (tmp_path / 'other.md').write_text('# 別\n\n## 日本語の 見出し（注）\n\n本文\n', encoding='utf-8')
    (tmp_path / 'plain.md').write_text('# 非対象\n', encoding='utf-8')
    return tmp_path


def conv(work, text, **kw):
    md = work / 'main.md'
    md.write_text(text, encoding='utf-8')
    targets = {(work / 'main.md').resolve(), (work / 'other.md').resolve()}
    return m2h.convert_file(md, is_target=targets.__contains__, **kw)


def test_image_embed_and_limits(work):
    r = conv(work, '# T\n\n## 図\n\n![説明](img.png)\n\n![ベクタ](pic.svg)\n\n![大](big.png)\n\n![無](none.png)\n')
    h = r.html
    assert 'src="data:image/png;base64,' in h
    assert 'src="data:image/svg+xml;base64,' in h
    assert 'src="big.png"' in h
    assert any('3 MB' in w for w in r.warnings)
    assert ('画像', 'none.png') in r.broken
    assert r.figures == 2


def test_named_figure_after_first_mention_only(work):
    text = ('# T\n\n## A\n\n最初に `fig_a.png` を述べる段落。\n\n続きの段落。\n\n- 項目 `fig_b.png` の説明\n- 次の項目\n\n'
            '再び `fig_a.png` に言及。\n\n```\n`fig_a.png` はコードの中\n```\n')
    h = conv(work, text).html
    assert h.count('<figure>') == 2
    i = h.index('<figure>')
    assert h.index('最初に') < i < h.index('続きの段落')
    j = h.index('<figure>', i + 1)
    assert h.index('次の項目') < j      # 箇条書きの後ろ


def test_names_only_in_code_block_ignored(work):
    text = '# T\n\n## A\n\n```\n`fig_a.png`\n```\n\n    `fig_b.png`\n'
    assert '<figure>' not in conv(work, text).html


def test_link_rewrite_and_anchor(work):
    text = ('# T\n\n## 節\n\n[a](other.md) [b](other.md#日本語の-見出し注) [c](plain.md) [d](nothing.md) '
            '[e](#節) [f](other.md#存在しない) [g](https://example.com/x.md)\n')
    r = conv(work, text)
    h = r.html
    assert 'href="other.html"' in h
    assert 'href="other.html#日本語の-見出し注"' in h
    assert 'href="plain.md"' in h
    assert 'href="#節"' in h
    assert 'href="https://example.com/x.md"' in h
    kinds = [b for b in r.broken]
    assert ('リンク', 'nothing.md') in kinds
    assert ('アンカー', 'other.md#存在しない') in kinds
    assert len(kinds) == 2


def test_japanese_heading_ids_github_style(work):
    h = conv(work, '# T\n\n## 節\n\n## 日本語の 見出し（注）\n\n## 節\n').html
    assert '<h2 id="節">' in h
    assert '<h2 id="日本語の-見出し注">' in h
    assert '<h2 id="節-1">' in h
    assert 'href="#節-1"' in h      # 目次


def test_table_sort_and_sticky(work):
    h = conv(work, '# T\n\n## 表\n\n| a | b |\n|---|---|\n| 1 | x |\n| 2 | y |\n').html
    assert '<thead>' in h and '<tbody>' in h
    assert 'position:sticky' in h and "addEventListener('click'" in h


def test_table_and_list_without_blank_line(work):
    h = conv(work, '# T\n\n## S\n\n説明:\n- 一\n  - 入れ子\n- 二\n\n表:\n| a |\n|---|\n| 1 |\n').html
    assert h.count('<li>') >= 3 and '<ul>' in h
    assert '<table>' in h


def test_fold(work):
    text = '# T\n\n## S\n\n### ラウンド 1\n\n内容1\n\n### 概要\n\n内容2\n\n### 反省\n\n内容3\n\n## 次\n\n内容4\n'
    h = conv(work, text, fold_regex='^(ラウンド|反省)').html
    assert h.count('<details>') == 2 and h.count('</details>') == 2
    assert h.index('内容1') < h.index('</details>') < h.index('内容2')
    assert h.index('内容3') < h.rindex('</details>') < h.index('内容4')
    assert '<details>' not in conv(work, text).html      # 既定は折りたたまない


def test_code_block_and_meta(work):
    h = conv(work, '# タイトル\n\n## S\n\n```python\nx = 1 < 2\n```\n').html
    assert '<pre>' in h and '&lt;' in h
    assert '<title>タイトル</title>' in h
    assert '元: main.md' in h
    assert not re.search(r'20\d\d-\d\d-\d\d', h) and not re.search(r'\bv\d{3}\b', h)
    assert 'http://' not in re.sub(r'xmlns="[^"]*"', '', h) and 'cdn' not in h.lower()
    assert (work / 'main.html').exists()


def test_overwrite_policy(work):
    (work / 'main.html').write_text('<html>手書き</html>', encoding='utf-8')
    r = conv(work, '# T\n')
    assert r.skipped and (work / 'main.html').read_text(encoding='utf-8') == '<html>手書き</html>'
    (work / 'main.html').unlink()
    conv(work, '# T\n')
    r = conv(work, '# T2\n')      # 生成物は上書きしてよい
    assert not r.skipped and 'T2' in (work / 'main.html').read_text(encoding='utf-8')


def test_exclusions():
    assert build_html.excluded('docs/history/x.md')
    assert build_html.excluded('CHANGELOG.md') and build_html.excluded('a/LICENSE.md')
    assert not build_html.excluded('docs/01_overview.md')


def test_config_fold_for_sextupole():
    cfg = build_html.load_config()
    c = cfg['example/ddba/v07_14/sextupole_v2/README.md']
    assert re.search(c['fold_regex'], 'ラウンド 3') and re.search(c['fold_regex'], '反省')
    assert re.search(c['fold_regex'], '揺らぎの確認 ab') and not re.search(c['fold_regex'], '概要')


# ---------------------------------------------------------------- 数式 (KaTeX。Node.js がなければ描画の試験は skip)
import mathtex  # noqa: E402

needs_node = pytest.mark.skipif(mathtex.find_node() is None, reason='Node.js が見つからない')


def _tex(text):
    return [(it['tex'], it['display']) for it in mathtex.extract(text)[1]]


def test_math_extract_rules():
    assert _tex(r'a $\nu_x$ b $x$ c $x=1$ d $(x,y)$ e $j+k+l+m=4$ f $0.36$') == [
        (r'\nu_x', False), ('x', False), ('x=1', False), ('(x,y)', False), ('j+k+l+m=4', False), ('0.36', False)]
    assert _tex('$$\n\\frac{a}{b}\n+ c\n$$') == [('\\frac{a}{b}\n+ c', True)]
    assert _tex('文 $$E=mc^2$$ 文') == [('E=mc^2', True)]
    # 数式にしないもの: コード・$$$・金額・単独の $・シェル変数・\$
    for s in ['`$x_1$`', '``$a_1$``', '```\n$x_1$\n$$\ny\n$$\n```', '    $x_2$\n', '$$$ と $$$$', '金額 $5 と $6、$1〜$13',
              '$', '$HOME と $PATH', r'\$5 と \$x_1\$', '$ x $', '$FORM と $Context', '$日本語$']:
        assert _tex(s) == [], s
    # 箇条書きの続きの字下げはコードではない
    assert _tex('- a\n\n    $x_1$ と\n') == [('x_1', False)]
    assert 'ZZMATHZZ' in mathtex.extract('$x_1$')[0] and '&#36;' in mathtex.extract(r'\$5')[0]


@needs_node
def test_math_render_inline_display_and_no_script(work):
    h = conv(work, '# T\n\n## S\n\nインライン $\\nu_x = a_1*b_2*c$ と $x$。\n\n$$\n\\sum_{n=1}^{\\infty}\\frac{1}{n^2}\n'
                   '=\\frac{\\pi^2}{6}\n$$\n\n- 項目 $E=mc^2$\n').html
    assert h.count('<span class="katex">') == 4
    assert 'katex-display' in h and '<math' in h             # 表示数式と、併記の MathML
    assert 'ZZMATHZZ' not in h and '$\\nu' not in h
    assert '<em>' not in h                                    # _ や * が強調にならない
    assert '@font-face' in h and 'data:font/woff2;base64,' in h
    assert 'katex' not in re.search(r'<script>(.*?)</script>', h, re.S).group(1).lower()   # 数式のための JS はない
    assert not re.search(r'url\((?!data:)', h) and 'http://' not in re.sub(r'xmlns="[^"]*"', '', h)
    assert 'https://' not in h


@needs_node
def test_math_in_table_and_pipe(work):
    h = conv(work, '# T\n\n## S\n\n| 量 | 式 |\n|---|---|\n| a | $\\nu_x = a|b$ |\n| b | $e^{-t/\\tau}$ と $\\alpha_c$ |\n').html
    assert h.count('<tr>') == 3 and h.count('<td>') == 4 and 'ZZMATHZZ' not in h
    assert h.count('<span class="katex">') == 3


@needs_node
def test_math_code_and_currency_untouched(work):
    text = '# T\n\n## S\n\n`$x_1$` と $5 と $6、$1〜$13、`$$$`。\\$x_1\\$\n\n```\n$x_1 = 2$\n```\n'
    r = conv(work, text)
    h = r.html
    assert 'katex' not in h and '@font-face' not in h
    assert '<code>$x_1$</code>' in h and '$5 と $6' in h and '$x_1$' in h and '<pre>' in h
    assert r.math_count == 0


@needs_node
def test_math_error_warns_and_keeps_tex(work):
    r = conv(work, '# T\n\n## S\n\nよい $x_1$ と悪い $\\undefinedcommand{x}$。\n\n$$\n\\badmacro\n$$\n')
    assert r.math_count == 3 and r.math_errors == 2
    assert sum('数式を変換できない' in w for w in r.warnings) == 2
    assert any('行 5' in w for w in r.warnings) and any('行 7' in w for w in r.warnings)
    assert 'katex-error' in r.html and '\\undefinedcommand' in r.html


@needs_node
def test_math_css_only_in_math_docs_and_font_subset(work):
    plain = conv(work, '# T\n\n本文だけ。\n').html
    assert '@font-face' not in plain and 'katex' not in plain
    small = conv(work, '# T\n\n$x$\n').html
    big = conv(work, '# T\n\n$\\mathbb{R}\\ \\mathcal{H}\\ \\mathfrak{g}\\ \\mathtt{c}\\ \\mathbf{J}$ $$\\sum_i \\left(\\frac{a}{b}\\right)$$\n').html
    n_small, n_big = small.count('@font-face'), big.count('@font-face')
    assert 1 <= n_small < n_big < 20          # 使われているフォントだけ
    face = 'font-family:KaTeX_AMS;font-style'
    assert face in big and face not in small


def test_math_without_node_keeps_raw(work, monkeypatch):
    monkeypatch.setattr(mathtex, 'find_node', lambda refresh=False: None)
    monkeypatch.setattr(m2h, '_warned_node', [])
    r = conv(work, '# T\n\n$\\nu_x$ と $y_1$\n')
    assert '$\\nu_x$' in r.html and 'katex' not in r.html and r.math_count == 0
    assert sum('Node.js が見つからない' in w for w in r.warnings) == 1
    assert not any('Node.js' in w for w in conv(work, '# T\n\n$x_1$\n').warnings)     # 警告は 1 回だけ


def test_find_node_env_order(monkeypatch, tmp_path):
    fake = tmp_path / 'node.exe'
    fake.write_text('x')
    monkeypatch.setenv('PYSAD_NODE', str(fake))
    assert mathtex.find_node(refresh=True) == fake
    monkeypatch.setenv('PYSAD_NODE', str(tmp_path / 'none'))
    assert mathtex.find_node(refresh=True) != tmp_path / 'none'     # 無効な PYSAD_NODE は飛ばして次を探す
    monkeypatch.undo()
    mathtex.find_node(refresh=True)


def test_config_math_off_for_help_index():
    assert build_html.load_config()['docs/34_sad_help_index.md'] == {'math': False}



# ---------------------------------------------------------------- コードの色分け (Pygments がなければ skip)
import codehl  # noqa: E402
import site_index  # noqa: E402

needs_pygments = pytest.mark.skipif(not codehl.AVAILABLE, reason='Pygments がない')


def _code_texts(h):
    """HTML の <pre><code> の中を、タグを除いて元の文字列に戻したものの一覧。"""
    import html as _h
    return [_h.unescape(re.sub(r'<[^>]+>', '', m)) for m in re.findall(r'<pre><code[^>]*>(.*?)</code></pre>', h, re.S)]


@needs_pygments
def test_highlight_python_bash_json_keep_text(work):
    py = 'import numpy as np\n# 注釈\nx = np.array([1, 2.5]) < "a&b"\n'
    sh = 'cd "$HOME"\npy -3.12 run.py --x 1 # 実行\n'
    js = '{"a": [1, 2], "b": "<c>"}\n'
    h = conv(work, '# T\n\n```python\n%s```\n\n```bash\n%s```\n\n```json\n%s```\n' % (py, sh, js)).html
    assert _code_texts(h) == [py, sh, js]               # コピーすると元の文字列
    assert '<span class="k">import</span>' in h and '<span class="c"># 注釈</span>' in h
    assert '<span class="s">' in h and '<span class="m">' in h
    assert '&lt;' in h and '&amp;' in h                  # 実体参照は保たれる


@needs_pygments
def test_highlight_sad_script_and_not_output(work):
    sad = 'FFS;\nUSE RING; CELL; CALC;   ! 光学\nQUAD QF = (L=0.3 K1=0.7);\nPrint["NX=", Twiss["NX","$$$"]/2/Pi];\n'
    out = '$$$ NX  6.696329\n  BX= 5.025185   EX= .565228;\n'
    table = '| a | b |\n|---|---|\n| 1 | 2 |\n'
    plain = 'Q1 ---- Q2\nx = 1;\n'
    md = '# T\n\n```\n%s```\n\n```\n%s```\n\n```\n%s```\n\n```text\n%s```\n' % (sad, out, table, plain)
    h = conv(work, md).html
    assert _code_texts(h) == [sad, out, table, plain]
    blocks = re.findall(r'<pre><code[^>]*>(.*?)</code></pre>', h, re.S)
    assert '<span' in blocks[0] and 'language-sad' in h
    assert all('<span' not in b for b in blocks[1:])    # 出力の例・表・text は色分けしない
    assert '<span class="c">! 光学</span>' in blocks[0]
    assert '<span class="t">QUAD</span>' in blocks[0] and '<span class="k">USE</span>' in blocks[0]
    assert '<span class="b">Print</span>' in blocks[0] and '<span class="s">"NX="</span>' in blocks[0]


@needs_pygments
def test_highlight_css_only_in_pages_with_colored_code(work):
    assert 'pre .k' not in conv(work, '# T\n\n本文だけ。\n').html
    assert 'pre .k' not in conv(work, '# T\n\n```\n$$$ 1 2 3\n```\n\n```text\nx\n```\n').html
    h = conv(work, '# T\n\n```python\nx = 1\n```\n').html
    assert 'pre .k' in h and 'prefers-color-scheme: dark' in h
    assert h.count('pre .k {') == 1


@needs_pygments
def test_sad_lexer_tokens_and_roundtrip():
    from pygments.token import Comment, Keyword, Name, Number, Operator, Punctuation, String
    import sad_lexer
    code = ('USE RING; ! comment != not\n(* block\n comment *) DRIFT D1 = (L=1.5E-3 K1=-.5);\n'
            'Module[{a}, a = "s\\"t" -> {1, 2}]; x != 1; 1.7MV\n')
    lx = sad_lexer.SADLexer(stripnl=False, ensurenl=False)
    allt = list(lx.get_tokens(code))
    assert ''.join(v for _, v in allt) == code
    toks = [(t, v) for t, v in allt if v.strip()]
    d = {v: t for t, v in toks}
    assert d['USE'] in Keyword and d['DRIFT'] in Keyword.Type and d['Module'] in Keyword
    assert d['! comment != not'] in Comment and any(t in Comment.Multiline for t, _ in toks)
    assert d['1.5E-3'] in Number and d['.5'] in Number and d['1.7'] in Number and d['MV'] in Name
    assert d['"s\\"t"'] in String and d['->'] in Operator and d['!='] in Operator
    assert d['['] in Punctuation and d[';'] in Punctuation and d['D1'] in Name and d['K1'] in Name


@needs_pygments
def test_looks_like_sad_rules():
    assert codehl.looks_like_sad('FFS; USE RING; CELL; CALC;\nPrint[1];\n')
    assert codehl.looks_like_sad('QUAD QF = (L=0.3 K1=0.7);\n! 注\nx = 1;\n')
    for s in ['$$$ 1 2\n', 'py -3.12 run.py\n', '1.0 2.0 3.0\n4.0 5.0 6.0\n', '+----+\n| a; |\n+----+\n',
              'NX=  6.69\nBX= 5.02\n', 'ただの文章。\n', '']:
        assert not codehl.looks_like_sad(s), s


# ---------------------------------------------------------------- 全文検索の索引・入口ページ
@pytest.fixture
def site(tmp_path):
    root = tmp_path
    (root / 'docs' / 'history').mkdir(parents=True)
    (root / 'example' / 'A_first').mkdir(parents=True)
    (root / 'example' / '_deletable').mkdir(parents=True)
    (root / 'physics_guide').mkdir()
    (root / 'tools').mkdir()
    (root / 'docs' / '01_overview.md').write_text(
        '# 01. 全体像\n\nこのプロジェクトの目的は、加速器の計算である。\n\n## 目次\n\n'
        '| 群 | 文書 | 内容 |\n|---|---|---|\n| 入門 | [01_overview.md](01_overview.md) | 本書 |\n| 物理 | [02_phys.md](02_phys.md) | 共鳴の式 |\n',
        encoding='utf-8')
    (root / 'docs' / '02_phys.md').write_text(
        '# 02. 物理\n\n導入の段落。\n\n## 共鳴とハミルトニアン\n\n共鳴は $\\nu_x = 6.8$ で起きる。\n\n```\nFFS; USE RING;\nPrint["Emittance"];\n```\n\n'
        '## FODO\n\n別の節。[概要](01_overview.md)\n', encoding='utf-8')
    (root / 'docs' / '03_extra.md').write_text('# 03. 表に無い文書\n\n本文。\n', encoding='utf-8')
    (root / 'docs' / 'history' / 'old.md').write_text('# 履歴\n', encoding='utf-8')
    (root / 'example' / 'README.md').write_text('# example 全体\n\n例の一覧。\n', encoding='utf-8')
    (root / 'example' / 'A_first' / 'README.md').write_text('# A 最初の例\n\n光学の例を集めた。\n', encoding='utf-8')
    (root / 'example' / '_deletable' / 'README.md').write_text('# 削除してよい\n\n作業用。\n', encoding='utf-8')
    (root / 'physics_guide' / '01_quad.md').write_text('# 四極の導出\n\n本文。\n', encoding='utf-8')
    (root / 'tools' / 'README.md').write_text('# tools\n\n道具。\n', encoding='utf-8')
    (root / 'README.md').write_text('# pysad\n\n説明。\n', encoding='utf-8')
    rels = sorted(p.relative_to(root).as_posix() for p in root.rglob('*.md') if 'history' not in p.parts)
    targets = {(root / r).resolve() for r in rels}
    pages = {}
    for r in rels:
        res = m2h.convert_file(root / r, is_target=targets.__contains__, root=root)
        if site_index.indexed(r):
            pages[r[:-3] + '.html'] = res.index
    srel = {r for r in rels if site_index.indexed(r)}
    msgs, warns = site_index.write_site(root, pages, srel)
    return root, srel, msgs, warns


def _search(idx, *words):
    """search.html と同じ規則 (空白区切りの語を全部含む節。大文字小文字・全角半角を区別しない部分一致)。"""
    import unicodedata
    ws = [unicodedata.normalize('NFKC', w).lower() for w in words]
    hits = []
    for s in idx['sections']:
        p = idx['pages'][s['page']]
        hay = (p['title'] + ' ' + s['heading'] + ' ' + s['text']).lower()
        if all(w in hay for w in ws):
            hits.append((p['path'], s['id']))
    return hits


def _load_idx(root):
    import json
    return json.loads((root / 'search_index.json').read_text(encoding='utf-8'))


def test_search_index_format_and_paths(site):
    root, srel, msgs, warns = site
    assert not warns
    idx = _load_idx(root)
    assert idx['version'] == 1 and set(idx) == {'version', 'pages', 'sections'}
    for p in idx['pages']:
        assert set(p) == {'path', 'title', 'desc', 'headings'} and (root / p['path']).is_file()
    for s in idx['sections']:
        assert set(s) == {'page', 'id', 'heading', 'text'} and 0 <= s['page'] < len(idx['pages'])
    paths = {p['path'] for p in idx['pages']}
    assert 'docs/02_phys.html' in paths and 'example/A_first/README.html' in paths and 'physics_guide/01_quad.html' in paths
    assert not any('_deletable' in p or 'history' in p for p in paths)
    assert site_index.check_site(root) == []


def test_search_finds_japanese_code_and_case(site):
    idx = _load_idx(site[0])
    assert ('docs/02_phys.html', '共鳴とハミルトニアン') in _search(idx, '共鳴')
    assert _search(idx, 'ハミルトニアン')[0][0] == 'docs/02_phys.html'
    assert _search(idx, 'fodo') and _search(idx, 'EMITTANCE')                  # 大文字小文字・コードの中
    assert _search(idx, 'Print[', 'USE')                                       # 複数語は AND
    assert not _search(idx, '存在しない語zzz') and not _search(idx, '共鳴', '作業用')


@needs_node
def test_search_finds_math_tex(site):
    import json
    idx = _load_idx(site[0])
    assert _search(idx, '\\nu_x')
    assert 'katex' not in ''.join(s['text'] for s in idx['sections']).lower() and 'ZZMATH' not in json.dumps(idx)


def test_page_header_links_and_index_page(site):
    root, srel, msgs, warns = site
    h = (root / 'docs' / '02_phys.html').read_text(encoding='utf-8')
    assert 'href="../search.html"' in h and 'href="../index.html"' in h and '元: 02_phys.md' in h
    top = (root / 'README.html').read_text(encoding='utf-8')
    assert 'href="search.html"' in top and 'href="index.html"' in top
    ix = (root / 'index.html').read_text(encoding='utf-8')
    for f in ('docs/01_overview.html', 'docs/02_phys.html', 'docs/03_extra.html', 'example/README.html', 'example/A_first/README.html',
              'tools/README.html', 'physics_guide/01_quad.html', 'README.html', 'search.html'):
        assert 'href="%s"' % f in ix, f
    assert '_deletable' not in ix and 'history' not in ix
    assert '加速器の計算である' in ix                       # 目的 (docs/01 の最初の段落)
    assert ix.index('入門') < ix.index('01. 全体像') < ix.index('物理') < ix.index('02. 物理') < ix.index('その他') < ix.index('03. 表に無い文書')
    assert 'http://' not in ix and ix.count('<script') == 0 and not re.search(r'href="/', ix)
    assert '<script>' in (root / 'search.html').read_text(encoding='utf-8')


def test_site_check_detects_broken(site):
    root = site[0]
    (root / 'docs' / '03_extra.html').unlink()
    kinds = site_index.check_site(root)
    assert ('索引のパス', 'docs/03_extra.html') in kinds and any(k == 'リンク' and '03_extra.html' in t for k, t in kinds)


def test_index_pages_do_not_overwrite_handwritten(site):
    root, srel, msgs, warns = site
    (root / 'index.html').write_text('<html>手書き</html>', encoding='utf-8')
    msgs, warns = site_index.write_site(root, {}, srel)
    assert (root / 'index.html').read_text(encoding='utf-8') == '<html>手書き</html>' and warns


def test_repo_index_lists_all_docs():
    """リポジトリの docs/ の全 md が入口ページの目次に載り、生成した索引・入口のリンク先が存在する。"""
    root = TOOLS.parent
    rels = {p.relative_to(root).as_posix() for p in (root / 'docs').glob('*.md')}
    assert len(rels) >= 36
    ix = site_index.build_index_page(root, rels)
    for r in rels:
        assert 'href="%s"' % (r[:-3] + '.html') in ix, r
    if (root / 'index.html').exists() and (root / 'search_index.json').exists():
        assert site_index.check_site(root) == []
