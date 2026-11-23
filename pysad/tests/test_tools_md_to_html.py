# -*- coding: utf-8 -*-
"""tools/md_to_html.py と tools/build_html.py の試験 (一時フォルダの小さな md と画像)。"""
import base64
import re
import sys
from pathlib import Path

import pytest

markdown = pytest.importorskip('markdown')
TOOLS = Path(__file__).resolve().parents[2] / 'tools'
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
