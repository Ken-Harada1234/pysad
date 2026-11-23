# -*- coding: utf-8 -*-
"""md の中の TeX の数式 ($...$、$$...$$) を、HTML を作るときに KaTeX (Node.js) で HTML の部品にして埋め込む。

開くときに JavaScript は要らない (数式の HTML・CSS・フォントはすべて HTML の中。外部通信なし)。
Node.js が見つからないときは、数式は従来どおり (生の文字のまま) にする。

使い方 (md_to_html.py から呼ばれる):
    find_node()          Node.js の実行ファイル (PYSAD_NODE → PATH の node → D:\\work\\tools\\node-*\\node.exe)
    extract(text)        数式を取り出してプレースホルダーに置き換える (コードの中は数式にしない)
    render(items, node)  数式をまとめて KaTeX にかける (Node.js を 1 回呼ぶ)
    restore(html, ...)   プレースホルダーを数式の HTML に戻す
    css_for(html)        使われているフォントだけを data URI で埋めた KaTeX の CSS
"""
import base64
import glob
import html as _html
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
KATEX_DIR = HERE / 'katex'
RENDER_JS = HERE / 'katex_render.js'
OUTPUT = 'htmlAndMathml'   # 見た目の HTML に、読み上げ・検索・コピー用の MathML を併記する ('html' だけだと aria-hidden で読み上げられない。数式の部分が約 18% 増える)
NODE_GLOB = r'D:\work\tools\node-*\node.exe'

_FENCE = re.compile(r'^\s*(```+|~~~+)')
_MARK = re.compile(r'^\s*([-*+]|\d+[.)])\s+\S')
_SHORT_VAR = re.compile(r'^[A-Za-z][A-Za-z0-9]{0,2}$')
_NUMBER = re.compile(r'^[+\-]?\d+(\.\d+)?$')
_NODE_CACHE = []


# ---------------------------------------------------------------- Node.js
def _ver_key(p):
    m = re.search(r'node-v(\d+)\.(\d+)\.(\d+)', p)
    return tuple(int(x) for x in m.groups()) if m else (0, 0, 0)


def find_node(refresh=False):
    """PYSAD_NODE → PATH の node → D:\\work\\tools\\node-*\\node.exe (版が最新のもの)。無ければ None。"""
    if _NODE_CACHE and not refresh:
        return _NODE_CACHE[0]
    found = None
    env = os.environ.get('PYSAD_NODE')
    if env and Path(env).is_file():
        found = Path(env)
    if found is None:
        w = shutil.which('node')
        if w:
            found = Path(w)
    if found is None:
        c = sorted(glob.glob(NODE_GLOB), key=_ver_key)
        if c:
            found = Path(c[-1])
    _NODE_CACHE[:] = [found]
    return found


# ---------------------------------------------------------------- 数式の検出
def _code_mask(lines):
    """各行がコードブロック (``` の囲み・4 つ空白) の行かどうか。箇条書きの続きの字下げはコードにしない。"""
    mask = []
    fence = None
    in_list = False
    prev_blank = True
    prev_code = False
    for line in lines:
        fm = _FENCE.match(line)
        if fence is not None:
            mask.append(True)
            if fm and fm.group(1)[0] == fence and line.strip().strip(fence) == '':
                fence = None
            prev_blank, prev_code = False, True
            continue
        if fm:
            fence = fm.group(1)[0]
            mask.append(True)
            prev_blank, prev_code = False, True
            continue
        s = line.strip()
        if not s:
            mask.append(prev_code and False)
            prev_blank = True
            continue
        indented = line.startswith(('    ', '\t'))
        if _MARK.match(line):
            in_list = True
            mask.append(False)
            prev_blank, prev_code = False, False
            continue
        if indented and (prev_blank or prev_code) and not in_list:
            mask.append(True)
            prev_blank, prev_code = False, True
            continue
        if not indented and not line.startswith(' '):
            in_list = False
        mask.append(False)
        prev_blank, prev_code = False, False
    return mask


def _is_math_inline(tex):
    """インラインの数式とみなすか。TeX の命令・_・^・{} を含む、または ASCII だけの短い式 (変数・変数どうしの
    式 `x=1`・`(x,y)`・`j+k+l+m=4`・`|x|>0.05` など)。日本語を含む・先頭が数字の金額らしいものは除く。"""
    if '\\' in tex or '_' in tex or '^' in tex or '{' in tex or '}' in tex:
        return True
    if not tex.isascii() or len(tex) > 60:
        return False
    if _SHORT_VAR.match(tex) or _NUMBER.match(tex) or re.match(r'^[<>=]{1,2}$', tex):
        return True
    return bool(re.search(r'[=+<>()|/,\'\-*!]', tex)) and bool(re.search(r'[A-Za-z0-9]', tex))


def _scan(seg, line0, items, ph):
    """コードでない一続きの本文 seg の数式を取り出し、置き換えた文字列を返す。"""
    out = []
    i = 0
    n = len(seg)
    while i < n:
        c = seg[i]
        if c == '\\' and i + 1 < n:
            if seg[i + 1] == '$':
                out.append('&#36;')
            else:
                out.append(seg[i:i + 2])
            i += 2
            continue
        if c == '`':
            j = i
            while j < n and seg[j] == '`':
                j += 1
            run = seg[i:j]
            end = j
            close = -1
            while True:
                k = seg.find(run, end)
                if k < 0:
                    break
                bl = seg.find('\n\n', j)
                if bl != -1 and bl < k:
                    break
                if (k + len(run) >= n or seg[k + len(run)] != '`') and seg[k - 1] != '`':
                    close = k
                    break
                end = k + 1
                while end < n and seg[end] == '`':
                    end += 1
            if close >= 0:
                out.append(seg[i:close + len(run)])
                i = close + len(run)
            else:
                out.append(run)
                i = j
            continue
        if c == '$':
            j = i
            while j < n and seg[j] == '$':
                j += 1
            r = j - i
            if r >= 3:
                out.append(seg[i:j])
                i = j
                continue
            if r == 2:
                k = j
                close = -1
                bl = seg.find('\n\n', j)
                lim = n if bl < 0 else bl
                while k < lim - 1:
                    if seg[k] == '\\':
                        k += 2
                        continue
                    if seg[k] == '$' and seg[k + 1] == '$' and (k + 2 >= n or seg[k + 2] != '$'):
                        close = k
                        break
                    k += 1
                tex = seg[j:close].strip() if close >= 0 else ''
                if close >= 0 and tex:
                    items.append({'tex': tex, 'display': True, 'line': line0 + seg.count('\n', 0, i)})
                    out.append(ph(len(items) - 1))
                    i = close + 2
                    continue
                out.append('$$')
                i = j
                continue
            # インライン
            nl = seg.find('\n', j)
            lim = n if nl < 0 else nl
            k = j
            close = -1
            if j < n and not seg[j].isspace():
                while k < lim:
                    if seg[k] == '\\':
                        k += 2
                        continue
                    if seg[k] == '`':
                        break
                    if seg[k] == '$':
                        close = k
                        break
                    k += 1
            ok = False
            if close > j:
                tex = seg[j:close]
                after = seg[close + 1] if close + 1 < n else ''
                ok = (not tex[-1].isspace() and not after.isdigit() and after != '$'
                      and _is_math_inline(tex))
            if ok:
                items.append({'tex': tex, 'display': False, 'line': line0 + seg.count('\n', 0, i)})
                out.append(ph(len(items) - 1))
                i = close + 1
            else:
                out.append('$')
                i = j
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def make_prefix(text):
    p = 'ZZMATHZZ'
    while p in text:
        p += 'Z'
    return p


def extract(text):
    """(置き換えた本文, 数式の一覧 [{tex, display, line}], プレースホルダーの接頭辞)。"""
    if '$' not in text:
        return text, [], ''
    prefix = make_prefix(text)

    def ph(k):
        return '%s%dEND' % (prefix, k)

    lines = text.split('\n')
    mask = _code_mask(lines)
    items = []
    out = []
    i = 0
    n = len(lines)
    while i < n:
        j = i
        while j < n and mask[j] == mask[i]:
            j += 1
        chunk = lines[i:j]
        if mask[i]:
            out.extend(chunk)
        else:
            seg = '\n'.join(chunk)
            out.extend(_scan(seg, i, items, ph).split('\n') if '$' in seg else chunk)
        i = j
    return '\n'.join(out), items, prefix


# ---------------------------------------------------------------- KaTeX の呼び出し
class MathError(Exception):
    pass


def render(items, node, output=None):
    """数式の一覧を KaTeX にかけ、[(html, err)] を返す (Node.js を 1 回だけ呼ぶ)。"""
    payload = json.dumps([{'tex': it['tex'], 'display': it['display']} for it in items]).encode('ascii')
    try:
        r = subprocess.run([str(node), str(RENDER_JS), output or OUTPUT], input=payload, capture_output=True,
                           timeout=300)
    except (OSError, subprocess.SubprocessError) as e:
        raise MathError('Node.js を実行できない: %r' % e)
    if r.returncode != 0:
        raise MathError('katex_render.js が失敗 (%d): %s' % (r.returncode, r.stderr.decode('utf-8', 'replace')[:300]))
    res = json.loads(r.stdout.decode('utf-8'))
    return [(x['html'], x['err']) for x in res]


def _fallback(item):
    d = '$$' if item['display'] else '$'
    return '<span class="katex-error" style="color:#cc0000">%s</span>' % _html.escape(d + item['tex'] + d)


def restore(h, items, rendered, prefix):
    """プレースホルダーを数式の HTML に戻す。表示数式だけの段落 <p>…</p> は <p> を外す。"""
    if not items:
        return h

    def body(k):
        html, err = rendered[k]
        return html if html else _fallback(items[k])
    pat = re.compile(r'<p>\s*%s(\d+)END\s*</p>' % re.escape(prefix))
    h = pat.sub(lambda m: body(int(m.group(1))) if items[int(m.group(1))]['display'] else m.group(0), h)
    return re.sub(r'%s(\d+)END' % re.escape(prefix), lambda m: body(int(m.group(1))), h)


def plain(s, items, prefix):
    """プレースホルダーを元の TeX (記号つき) に戻す (タイトルなど、HTML にしないところ用)。"""
    def sub(m):
        it = items[int(m.group(1))]
        d = '$$' if it['display'] else '$'
        return d + it['tex'] + d
    return re.sub(r'%s(\d+)END' % re.escape(prefix), sub, s)


# ---------------------------------------------------------------- CSS とフォント
_FACE = re.compile(r'@font-face\{[^}]*\}')
_FONT_CACHE = {}
_CSS_CACHE = []


def _load_css():
    if _CSS_CACHE:
        return _CSS_CACHE[0]
    css = (KATEX_DIR / 'katex.min.css').read_text(encoding='utf-8')
    faces = {}
    for m in _FACE.finditer(css):
        b = m.group(0)
        fam = re.search(r'font-family:"?(\w+)', b).group(1)
        wt = 'bold' if re.search(r'font-weight:(700|bold)', b) else 'normal'
        st = 'italic' if 'font-style:italic' in b else 'normal'
        f = re.search(r'url\(fonts/([\w\-]+\.woff2)\)', b).group(1)
        faces[(fam, wt, st)] = (f, b)
    rest = _FACE.sub('', css)
    rules = []   # (使うクラスの集合, 家族, 太さ, 字体)
    for m in re.finditer(r'([^{}]+)\{([^{}]*)\}', rest):
        decl = m.group(2)
        fm = re.search(r'font-family:\s*"?(KaTeX_\w+)', decl)
        if not fm:
            continue
        wt = 'bold' if re.search(r'font-weight:\s*(700|bold)', decl) else 'normal'
        st = 'italic' if re.search(r'font-style:\s*italic', decl) else 'normal'
        for sel in m.group(1).split(','):
            cls = set(re.findall(r'\.([\w\-]+)', sel)) - {'katex'}
            rules.append((cls, fm.group(1), wt, st))
    _CSS_CACHE[:] = [(rest, faces, rules)]
    return _CSS_CACHE[0]


def needed_faces(h):
    """KaTeX の HTML の class から、必要なフォントの面 (家族, 太さ, 字体) の集合を判定する。"""
    rest, faces, rules = _load_css()
    used = set()
    for m in re.finditer(r'class="([^"]*)"', h):
        used.update(m.group(1).split())
    need = {('KaTeX_Main', 'normal', 'normal')}
    for cls, fam, wt, st in rules:
        if cls <= used:
            need.add((fam, wt, st))
    # 太字と斜体が同じ家族で両方あるときは、太字の斜体も使われうる
    for fam in ('KaTeX_Main', 'KaTeX_Math'):
        if (fam, 'bold', 'normal') in need and (fam, 'normal', 'italic') in need:
            need.add((fam, 'bold', 'italic'))
    return {k for k in need if k in faces}


def css_for(h):
    """h (KaTeX の HTML 断片をつないだもの) に必要なフォントだけを data URI で埋めた CSS を返す。"""
    rest, faces, rules = _load_css()
    parts = []
    for key in sorted(needed_faces(h)):
        fname, block = faces[key]
        if fname not in _FONT_CACHE:
            _FONT_CACHE[fname] = base64.b64encode((KATEX_DIR / 'fonts' / fname).read_bytes()).decode('ascii')
        uri = 'data:font/woff2;base64,' + _FONT_CACHE[fname]
        parts.append(re.sub(r'src:.*\}$', 'src:url(%s) format("woff2")}' % uri, block))
    extra = ('.katex-display{overflow-x:auto;overflow-y:hidden;padding:.2em 0}'
             '.katex-error{font-family:Consolas,"Cascadia Mono",monospace}'
             '@media (prefers-color-scheme: dark){.katex-error{color:#ff7b72 !important}}')
    return ''.join(parts) + rest + extra
