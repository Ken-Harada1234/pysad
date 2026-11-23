# -*- coding: utf-8 -*-
"""コードブロックの色分け (ビルドのときに <span class="..."> にする。開くときの JavaScript は要らない)。

- 対象: ```python・```bash・```sh・```json・```sad (SAD の台本)。
- 言語指定のないブロックは、SAD の台本と判断できるときだけ SAD として色分けする (looks_like_sad)。
  出力の例 (`$$$ ...` や数値の表)・コマンドの行・図・text は色分けしない。
- 色分けしても文字列は変わらない (タグを除いて HTML の実体参照を戻すと元の文字列。合わなければ色分けしない)。
- Pygments がなければ何もしない (色分けなしで HTML は作れる)。
"""
import html as _html
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    from pygments import lex
    from pygments.lexers import BashLexer, JsonLexer, PythonLexer
    from pygments.token import Token
    import sad_lexer
    AVAILABLE = True
except ImportError:      # pragma: no cover
    AVAILABLE = False

# 明暗の配色 (色覚に配慮して、色相だけに頼らず、コメントは斜体、キーワードは太字にする)
CSS = """
:root { --c-kw:#0550ae; --c-ty:#7a3e9d; --c-st:#8a4b00; --c-nu:#1b7a6e; --c-bi:#1f6f8b; --c-cm:#59636e; --c-fn:#6e4aa5; }
@media (prefers-color-scheme: dark) { :root { --c-kw:#79c0ff; --c-ty:#d2a8ff; --c-st:#e3b341; --c-nu:#7ee0c8; --c-bi:#56d4dd; --c-cm:#9198a1; --c-fn:#d2a8ff; } }
pre .c { color:var(--c-cm); font-style:italic; }
pre .k { color:var(--c-kw); font-weight:600; }
pre .t { color:var(--c-ty); font-weight:600; }
pre .s { color:var(--c-st); }
pre .m { color:var(--c-nu); }
pre .b { color:var(--c-bi); }
pre .f { color:var(--c-fn); }
pre .o { opacity:.8; }
"""

# Pygments のトークン → 短いクラス (c コメント、k キーワード、t 型・要素、s 文字列、m 数、b 組み込み、f 関数・クラスの名前、o 演算子)


def _cls(tok):
    if tok in Token.Comment:
        return 'c'
    if tok in Token.Keyword.Type or tok in Token.Name.Class:
        return 't' if tok in Token.Keyword.Type else 'f'
    if tok in Token.Keyword:
        return 'k'
    if tok in Token.Operator.Word:
        return 'k'
    if tok in Token.Literal.String:
        return 's'
    if tok in Token.Literal.Number:
        return 'm'
    if tok in Token.Name.Builtin or tok in Token.Name.Variable:
        return 'b'
    if tok in Token.Name.Function or tok in Token.Name.Decorator:
        return 'f'
    if tok in Token.Name.Tag or tok in Token.Name.Attribute:
        return 'k'
    if tok in Token.Operator:
        return 'o'
    return ''


_LEXERS = {
    'python': lambda: PythonLexer(stripnl=False, ensurenl=False),
    'py': lambda: PythonLexer(stripnl=False, ensurenl=False),
    'bash': lambda: BashLexer(stripnl=False, ensurenl=False),
    'sh': lambda: BashLexer(stripnl=False, ensurenl=False),
    'shell': lambda: BashLexer(stripnl=False, ensurenl=False),
    'json': lambda: JsonLexer(stripnl=False, ensurenl=False),
    'sad': lambda: sad_lexer.SADLexer(stripnl=False, ensurenl=False),
}

# ---------------------------------------------------------------- 言語指定のないブロックの判定
_STMT_END = re.compile(r';\s*(?:!.*|\(\*.*\*\))?\s*$')
_SAD_START = re.compile(r'^\s*(?:FFS|USE|CALC(?:ULATE)?|EMIT|DRAW|DISPLAY|SHOW|GO|FIT|FREE|FIX|MATCH|RFSW|NORFSW|RAD|NORAD|'
                        r'FLUC|NOFLUC|COD|NOCOD|INS|CELL|PRINT|OUTPUT|CLOSE|APPEND|READ|SAVE|MARK|ORG|VISIT|DRIFT|BEND|QUAD|'
                        r'SEXT|OCT|DECA|DODECA|MULT|SOL|CAVI|TCAVI|APERT|COORD|LINE|UND|MAP|BEAMBEAM|MONI|MOVE|[A-Z][A-Za-z0-9]*'
                        r'\[|[A-Za-z_$][\w$.]*\s*(?::?=|\+=|-=|\*=)|EXPAND|TRPT|WAKE|LWAKE|TWAKE)')
_OUTPUT_MARK = re.compile(r'^\s*\$\$\$|^\s*[-=+|]{5,}\s*$|^\s*\|.*\|\s*$|^\s*[A-Z]{1,3}\s+\S+\s+[-\d.]+\s+\|')
_NUMERIC = re.compile(r'^[\s\d.+\-eEdD|,:%]+$')
_SHELL_START = re.compile(r'^\s*(?:py |py$|python|pip |cd |set |git |wsl |export |QUICK=|PYSAD_|MSYS|PYTHONIOENCODING|\$ |>|mpirun|'
                          r'[A-Za-z]:\\|\./|npm )')


def looks_like_sad(text):
    """言語指定のないコードブロックが、SAD の台本と判断できるか (迷うときは False = 色分けしない)。"""
    lines = [l for l in text.split('\n') if l.strip()]
    if not lines or len(lines) > 400:
        return False
    if any(_OUTPUT_MARK.search(l) for l in lines):
        return False
    if sum(1 for l in lines if _SHELL_START.match(l)) * 2 >= len(lines):
        return False
    ends = [l for l in lines if _STMT_END.search(l)]
    if not ends:
        return False
    scripted = 0
    for l in lines:
        s = l.strip()
        if _STMT_END.search(l) or s.startswith(('!', '(*')) or _SAD_START.match(l):
            scripted += 1
    numeric = sum(1 for l in lines if _NUMERIC.match(l))
    return scripted * 10 >= len(lines) * 7 and numeric * 4 < len(lines) and len(ends) * 3 >= len(lines)


# ---------------------------------------------------------------- 色分け
def highlight_text(code, lang):
    """code (展開済みの文字列) を、<span> つきの HTML にする。元の文字列に戻らなければ None。"""
    mk = _LEXERS.get(lang)
    if mk is None:
        return None
    out = []
    try:
        for tok, val in lex(code, mk()):
            if not val:
                continue
            c = _cls(tok)
            esc = _html.escape(val, quote=False)
            out.append('<span class="%s">%s</span>' % (c, esc) if c and val.strip() else esc)
    except Exception:      # noqa: BLE001
        return None
    h = ''.join(out)
    plain = _html.unescape(re.sub(r'<[^>]+>', '', h))
    if plain != code:
        return None
    return h


_PRE = re.compile(r'<pre><code(?: class="language-([^"]*)")?>(.*?)</code></pre>', re.S)


def highlight_body(body):
    """HTML の本文の <pre><code> を色分けする。(新しい本文, 色分けしたブロックの数) を返す。"""
    if not AVAILABLE:
        return body, 0
    n = [0]

    def sub(m):
        lang = (m.group(1) or '').lower()
        code = _html.unescape(m.group(2))
        if lang in ('', 'text', 'txt', 'console', 'output', 'plain'):
            if lang or not looks_like_sad(code):
                return m.group(0)
            lang = 'sad'
        h = highlight_text(code, lang)
        if h is None or '<span' not in h:
            return m.group(0)
        n[0] += 1
        return '<pre><code class="language-%s">%s</code></pre>' % (lang, h)
    return _PRE.sub(sub, body), n[0]
