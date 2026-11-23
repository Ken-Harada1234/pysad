# -*- coding: utf-8 -*-
"""README.md を、図 (results/*.png) を base64 で埋め込んだ 1 枚の HTML にする。

使い方 (このフォルダで):
    py -3.12 make_html.py                      # README.md -> README_ja_v001.html
    py -3.12 make_html.py README.md out.html   # 入力・出力を指定

- 本文中に `fig_xxx.png` のように (バッククォートで) 書かれた図の名前を探し、最初に出てくる段落・箇条書きの直後に図を置く
  (コードブロックの中の名前は無視する)。results/ にない名前は飛ばす。
- 目次は左に固定、表は見出し行を固定し、見出しをクリックすると並べ替えられる。
- 「ラウンド」の見出し (### ラウンド…) は折りたたみ (<details>)。
- 外部ライブラリ・通信は使わない (markdown パッケージだけ要る)。HTML は md から作るので、md を直したら作り直す。
既存のファイルは上書きしない: 出力が既にあるときは、版 (v001 → v002 …) を上げた名前にする。
"""
import base64
import html as _html
import os
import re
import sys
from datetime import date

import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
CSS = """
:root { --fg:#1f2328; --muted:#59636e; --line:#d1d9e0; --bg:#ffffff; --bg2:#f6f8fa; --acc:#0969da; }
@media (prefers-color-scheme: dark) { :root { --fg:#e6edf3; --muted:#9198a1; --line:#3d444d; --bg:#0d1117; --bg2:#151b23; --acc:#4493f8; } }
* { box-sizing: border-box; }
body { margin:0; color:var(--fg); background:var(--bg); font:15px/1.7 "Yu Gothic UI","Meiryo","Hiragino Sans",system-ui,sans-serif; }
nav { position:fixed; top:0; left:0; bottom:0; width:300px; overflow:auto; padding:16px 14px; border-right:1px solid var(--line); background:var(--bg2); font-size:13px; }
nav ul { list-style:none; margin:0; padding-left:12px; } nav > ul { padding-left:0; }
nav a { color:var(--fg); text-decoration:none; display:block; padding:2px 4px; border-radius:4px; } nav a:hover { background:var(--line); }
nav .ttl { font-weight:700; margin-bottom:8px; color:var(--muted); }
main { margin-left:300px; padding:24px 40px 80px; max-width:1180px; }
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
.meta { color:var(--muted); font-size:12.5px; margin-bottom:1.5em; }
@media (max-width: 900px) { nav { position:static; width:auto; border-right:0; } main { margin-left:0; padding:16px; } }
"""
JS = """
document.querySelectorAll('table').forEach(function (t) {
  t.querySelectorAll('th').forEach(function (th, i) {
    th.title = 'クリックで並べ替え';
    th.addEventListener('click', function () {
      var rows = Array.from(t.querySelectorAll('tbody tr')); var asc = th.dataset.asc !== '1'; th.dataset.asc = asc ? '1' : '0';
      var num = function (s) { var v = parseFloat(String(s).replace(/[,×倍%]/g, '')); return isNaN(v) ? null : v; };
      rows.sort(function (a, b) {
        var x = a.children[i].textContent.trim(), y = b.children[i].textContent.trim(), nx = num(x), ny = num(y);
        var c = (nx !== null && ny !== null) ? nx - ny : x.localeCompare(y, 'ja'); return asc ? c : -c; });
      var tb = t.querySelector('tbody'); rows.forEach(function (r) { tb.appendChild(r); });
    });
  });
});
"""


def _in_pre(h, pos):
    a = h.rfind('<pre>', 0, pos)
    return a != -1 and h.find('</pre>', a, pos) == -1


def embed_figures(h, resdir):
    names = []
    for m in re.finditer(r'<code>(fig_[A-Za-z0-9_]+\.png)</code>', h):
        if m.group(1) not in names:
            names.append(m.group(1))
    done = 0
    for name in names:
        path = os.path.join(resdir, name)
        if not os.path.exists(path):
            continue
        tag = '<code>%s</code>' % name
        pos = -1
        start = 0
        while True:
            p = h.find(tag, start)
            if p == -1:
                break
            if not _in_pre(h, p):
                pos = p
                break
            start = p + 1
        if pos == -1:
            continue
        ends = [(h.find('</p>', pos), '</p>'), (h.find('</li>', pos), '</li>')]
        ends = [(e, t) for e, t in ends if e != -1]
        e, t = min(ends)
        ins = e + len(t)
        if t == '</li>':
            cl = [(h.find(x, ins), x) for x in ('</ul>', '</ol>')]
            cl = [(c, x) for c, x in cl if c != -1]
            if cl:
                c, x = min(cl)
                ins = c + len(x)
        data = base64.b64encode(open(path, 'rb').read()).decode('ascii')
        fig = ('\n<figure><img alt="%s" src="data:image/png;base64,%s"><figcaption>%s</figcaption></figure>\n'
               % (_html.escape(name), data, _html.escape(name)))
        h = h[:ins] + fig + h[ins:]
        done += 1
    return h, done


def fold_rounds(h):
    """### ラウンド… の節を <details> にする (次の h2/h3 まで)。"""
    parts = re.split(r'(?=<h[23] id=)', h)
    out = []
    for part in parts:
        m = re.match(r'<h3 id="[^"]*">(ラウンド[^<]*|揺らぎの確認[^<]*|反省[^<]*)</h3>', part)
        if m:
            body = part[m.end():]
            out.append('<details><summary>%s</summary>%s</details>' % (m.group(1), body))
        else:
            out.append(part)
    return ''.join(out)


def add_tbody(h):
    """ markdown の表は <thead> と <tbody> を持つ。持たなければ作る (並べ替え用)。"""
    return h


def main(argv):
    src = argv[1] if len(argv) > 1 else os.path.join(HERE, 'README.md')
    out = argv[2] if len(argv) > 2 else None
    if out is None:
        v = 1
        while True:
            out = os.path.join(HERE, 'README_ja_v%03d.html' % v)
            if not os.path.exists(out):
                break
            v += 1
        if v > 1:
            print('既存の版があるので v%03d にします' % v)
    elif os.path.exists(out):
        print('出力が既にあります (上書きしません):', out)
        return 1
    text = open(src, encoding='utf-8').read().replace('\r\n', '\n')
    md = markdown.Markdown(extensions=['tables', 'fenced_code', 'toc', 'sane_lists'],
                           extension_configs={'toc': {'toc_depth': '2-3', 'permalink': False}})
    body = md.convert(text)
    title = re.search(r'^# (.+)$', text, re.M)
    title = title.group(1) if title else os.path.basename(src)
    body, nfig = embed_figures(body, os.path.join(os.path.dirname(os.path.abspath(src)), 'results'))
    body = fold_rounds(body)
    page = ('<!DOCTYPE html>\n<html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>%s</title><style>%s</style></head><body>\n<nav><div class="ttl">目次</div>%s</nav>\n<main>\n'
            '<div class="meta">%s の README.md から自動で作った HTML (%s)。図は埋め込み済み (%d 枚)。表の見出しをクリックすると並べ替えられる。'
            '元は README.md。</div>\n%s\n</main><script>%s</script></body></html>\n'
            % (_html.escape(title), CSS, md.toc, _html.escape(os.path.basename(os.path.dirname(os.path.abspath(src)))),
               date.today().isoformat(), nfig, body, JS))
    with open(out, 'w', encoding='utf-8', newline='\n') as f:
        f.write(page)
    print('wrote', out, '%.1f KB' % (len(page.encode('utf-8')) / 1024), 'figures:', nfig)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
