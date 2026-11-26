"""Small FFS commands: VARIABLES, ATTRIBUTE, MATRIX, SHOW, VARY, RECOVER, REVERSE, PRINT and EXECUTE (translation of K64 sources).
FFS の小さなコマンド: VARIABLES(VAR/VARS)、ATTRIBUTE(ATTR)、MATRIX(MAT)、SHOW、VARY、RECOVER(REC)、REVERSE(REV)、PRINT(PRI)、EXECUTE(EXEC)。
K64 の tfvars.f、tfattr.f、tfmat.f、tshow.f、tffsa.f の直訳 (出力の書式は autofg ('12.9' など) まで gs と同じ)。

pysad の仕様・未対応:
  - PRINT の `Out[n]:= 値` の n は、USE が `Out[1]:= USE` を出したものとして 2 から数える (pysad は USE でこの行を出さない)。
  - VARY の既定キーは FFS のインスタンスの上書き表に入れる (MULT などキーの多い要素は gs と同じ検査)。
  - STATUS・SPLIT (K64 では obsolete)・RENUMBER・MEASURE は未対応 (従来どおり警告して無視)。
"""
from __future__ import annotations

import math
import re

from .sadscript import SADError, fmt_real, Null

PI2 = 2.0 * math.pi


def _pad(s, n):
    return (s + ' ' * n)[:n]


def _match_name(pat, name):
    """tmatch: * (任意)、% (1 文字)、{a,b} (選択)"""
    rx = ''
    i = 0
    while i < len(pat):
        ch = pat[i]
        if ch == '*':
            rx += '.*'
        elif ch == '%':
            rx += '.'
        elif ch == '{':
            j = pat.index('}', i)
            rx += '(?:' + '|'.join(re.escape(x) for x in pat[i + 1:j].split(',')) + ')'
            i = j
        else:
            rx += re.escape(ch)
        i += 1
    return re.fullmatch(rx, name, re.I) is not None


def _saved_value(ffs, el0, key, default):
    """VAR の Saved 欄: 設計値 (SAVE で更新される。gs で確認: FREE・GO では変わらず、SAVE か定義の値)。"""
    sv = getattr(el0, '_saved', None)
    if sv and key in sv:
        return float(sv[key])
    raw = getattr(ffs.env, '_elements', {}).get(el0.name)
    if raw and key in raw:
        try:
            return float(raw[key])
        except (TypeError, ValueError):
            pass
    return default


def cmd_variables(ffs, words, i):
    m = ffs.matcher_()
    pat = None
    n = 1
    if i + 1 < len(words) and words[i + 1].upper() not in _stop_words():
        pat = words[i + 1]
        n = 2
    ffs.write('!Variable     Keyword    Now          !   Previous    Saved         Minimum      Maximum       Couple      Coefficient\n')
    rows = []
    for v in m.variables:
        if pat is not None and not _match_name(pat, v.get('label', v['family'])):
            continue
        el0 = v['els'][0]
        now_el = float(ffs.get_value(el0, v['key']))
        saved = _saved_value(ffs, el0, v['key'], now_el)
        now = now_el if v.get('noop') else float(v.get('val', now_el))      # 何もしない変数の Now は結合先 (先頭の成分) の値
        prev = float(v.get('val2', now))
        lo, hi = v['vlim'][0], v['vlim'][1]
        # Couple の欄: 既定のキーの変数は、その成分が結合している先 (自分が先頭なら <--)。ほかのキー (ROTATE や BEND の K1) は空 (gs)
        rows.append((v.get('label', v['family']), v['key'], now, prev, saved, lo, hi,
                     v.get('couple') or ('<--' if v['key'] == (ffs.default_key(el0) or '').upper() else ''), 1.0))
    if not rows and m.variables and pat is not None and not any(c in pat for c in '*%{'):
        pass
    for name, key, now, prev, saved, lo, hi, nc, co in rows:
        v1 = fmt_real(now, '15.12')
        v2 = fmt_real(prev, '12.9')
        v3 = fmt_real(saved, '12.9')
        v4 = fmt_real(lo, '12.8')
        v5 = fmt_real(hi, '12.8')
        v6 = fmt_real(co, '12.8')
        ffs.write(' ' + _pad(name, 12) + ' ' + _pad(key, 8) + _pad(v1, 15) + ' !' + _pad(v2, 12) + _pad(v3, 12) + ' ' + _pad(v4, 12)
                  + '  ' + _pad(v5, 12) + ' ' + _pad(nc, 12) + _pad(v6, 12) + '\n')
    return n


def _stop_words():
    from .ffscmd import _COMMANDS, _FUNS, _ERRCMDS
    return set(_COMMANDS) | set(_FUNS) | set(_ERRCMDS)


def cmd_attribute(ffs, words, i):
    pats = []
    n = 1
    stop = _stop_words()
    while i + n < len(words) and words[i + n].upper() not in stop:
        pats += [p for p in words[i + n].split(',') if p]
        n += 1
    ffs.write('Element     Keyword    Value     Mimimum  Maximum  Couple         Coefficient\n')
    seen = set()
    for j in range(1, ffs.nlat):
        el = ffs.els[j - 1]
        fam = el.name
        if fam in seen:
            continue
        name = ffs.elname(j)
        if pats and not any(_match_name(p, name) or _match_name(p, fam) for p in pats):
            continue
        seen.add(fam)
        key = ffs.default_key(el)
        if not key:
            continue
        v = float(ffs.get_value(el, key))
        lim = getattr(el, '_vlim', {}).get(key, [-1e10, 1e10])
        ffs.write(_pad(name, 12) + _pad(key, 8) + _pad(fmt_real(v, '10.6'), 10) + _pad(fmt_real(lim[0], '10.6'), 10)
                  + _pad(fmt_real(lim[1], '10.6'), 10) + ' ' + _pad('<--', 12) + fmt_real(1.0, '10.6') + '\n')
    return n


def cmd_matrix(ffs, words, i):
    """Implement `MATRIX [S|P] [from to]`: show the transfer matrix (4x5, or 6x6 for CALC6D).
    MATRIX [S|P] [from to]: 転送行列 (4×5。CALC6D は 6×6)"""
    import numpy as np
    from .sadscript import S
    n = 1
    symp = True
    title = 'Symplectic'
    if i + n < len(words) and re.fullmatch(r'S(?:Y(?:M(?:P(?:L(?:E(?:C(?:T(?:I(?:C)?)?)?)?)?)?)?)?)?', words[i + n].upper()):
        n += 1
    elif i + n < len(words) and re.fullmatch(r'P(?:H(?:Y(?:S(?:I(?:C(?:A(?:L)?)?)?)?)?)?)?', words[i + n].upper()):
        symp = False
        title = 'Physical'
        n += 1
    l1, l2 = 1, ffs.nlat

    def pos(w):
        # 名前だけ (QF) は NAME.1。ラティスに無い名前は位置にならない (gs: MATRIX PHYSICAL MKSFA1 MKSFA2 は MKSFA2 が無いと to=$$$、エラーなしで読み飛ばす)
        k = ffs._position(w)
        if k is None and '.' not in w and w not in ('$$$', '***', '^^^'):
            k = ffs._position(w + '.1')
        return k
    if i + n < len(words):
        k1 = pos(words[i + n])
        if k1 is not None:
            l1 = int(k1[0] if isinstance(k1, tuple) else k1)
            n += 1
            if i + n < len(words):
                k2 = pos(words[i + n])
                if k2 is not None:
                    l2 = int(k2[0] if isinstance(k2, tuple) else k2)
                    n += 1
                elif re.fullmatch(r'[A-Za-z_][A-Za-z0-9_$]*', words[i + n]):
                    n += 1                               # ラティスに無い名前は gs は読み飛ばす (to=$$$)
    ffs.calc(quiet=True) if getattr(ffs, 'rows', None) is None else None
    B = ffs.I.builtins
    r = 1.0
    if not symp:
        g1, g2 = ffs._gammab_at(l1), ffs._gammab_at(l2)
        r = math.sqrt(g1 / g2)
    ffs.write(' ' + title + ' transfer matrix from ' + ffs.elname(l1) + ' to ' + ffs.elname(l2) + '\n')
    m4 = np.array(B['TransferMatrix']([float(l1), float(l2)]), dtype=float)
    calc6d = bool(ffs.flags.get('CALC6D'))
    if calc6d:
        m6 = m4
        for a in range(6):
            ffs.write(' ' + ''.join(_pad(fmt_real(r * float(m6[a, b]), '11.6'), 12) for b in range(6)) + '\n')
        return n
    # 5 列目: 分散を転送した項 (tftmat の第 5 列) = η(l2) - M η(l1)
    p1 = np.array(ffs._optfn._p28(l1), dtype=float)
    p2 = np.array(ffs._optfn._p28(l2), dtype=float)
    eta1 = p1[6:10]
    eta2 = p2[6:10]
    col5 = eta2 - m4[:4, :4] @ eta1
    for a in range(4):
        vals = [float(m4[a, b]) for b in range(4)] + [float(col5[a])]
        ffs.write(' ' + ''.join(_pad(fmt_real(r * v, '11.6'), 12) for v in vals) + '\n')
    return n


def cmd_show(ffs, words, i):
    ffs.write('!    component1   component2   fun         goal-value     np     scale\n')
    m = ffs.matcher_()
    for c in m.conditions:
        if c['mfitp'] == 0:
            continue
        kp, kp1, kf = c['kp'], c['kp1'], c['kf']
        sc = PI2 if kf in ('NX', 'NY', 'NZ') else (math.pi / 180.0 if kf.startswith('CHI') else 1.0)
        name = _pad(ffs.elname(kp), 12)
        name1 = _pad(ffs.elname(kp1) if kp != kp1 else ' ', 12)
        fun = kf + ('M' if c['mfitp'] < 0 else '')
        if c['fmin'] == c['fmax']:
            vout = fmt_real(float(c['fmin']) / sc, '12.9')
        else:
            vout = '{' + fmt_real(float(c['fmin']) / sc, 'S7.5') + ',' + fmt_real(float(c['fmax']) / sc, 'S7.5') + '}'
        sout = fmt_real(sc, '12.9')
        ffs.write(' FIT ' + name + ' ' + name1 + ' ' + _pad(fun, 8) + ' ' + _pad(vout, 17) + '%3d' % (abs(c['mfitp']) - 1) + ' ! *' + sout + '\n')
    return 1


def cmd_vary(ffs, words, i):
    """Implement `VARY keyword element-pattern`: change the default keyword of elements.
    VARY keyword element-pattern: 要素の既定キーを変える"""
    from .sadkw import KEYWORDS
    if i + 2 >= len(words):
        ffs.write('?Missing keyword or element pattern for VARY\n')
        return len(words) - i
    key = words[i + 1].upper()
    pat = words[i + 2]
    fams = ffs.families(pat)
    done = False
    for fam in fams:
        el = ffs.instances(fam)[0]
        tn = type(el).__name__
        from .ffscmd import _TYPENAME, _cls
        if key in KEYWORDS.get(_TYPENAME.get(_cls(el), ''), ()):
            for e in ffs.instances(fam):
                e._vary_key = key
            done = True
    if not done:
        ffs.write(' ???-FFS-Error-No element for VAR_Y with keyword: ' + words[i + 1] + '\n')
    return 3


def cmd_recover(ffs, words, i):
    m = ffs.matcher_()
    for v in m.variables:
        now = float(v.get('val', 0.0))
        prev = float(v.get('val2', now))
        v['val'], v['val2'] = prev, now
        ffs.set_var(v, prev)
    ffs.invalidate()
    return 1


def cmd_reverse(ffs, words, i):
    """Implement `REVERSE`: flip the signs of AX AY EPX EPY R2 R3 DPX DPY at the beamline entrance.
    REVERSE: ビームラインの入口で AX AY EPX EPY R2 R3 DPX DPY の符号を反転"""
    ini = ffs.I.globals
    for k in ('AX', 'AY', 'EPX', 'EPY', 'R2', 'R3', 'DPX', 'DPY'):
        if k in ini and isinstance(ini[k], float):
            ini[k] = -ini[k]
    for attr in ('axi', 'ayi', 'epxi', 'epyi', 'r2i', 'r3i', 'dpxi', 'dpyi'):
        if hasattr(ffs, attr):
            setattr(ffs, attr, -getattr(ffs, attr))
    ffs.invalidate()
    return 1


def cmd_print(ffs, words, i):
    """Implement `PRINT expr`: evaluate an expression and print it as `Out[n]:= value`.
    PRINT expr: 式を評価して `Out[n]:= 値` を出す"""
    text = ' '.join(words[i + 1:])
    try:
        v = ffs.I.run(text, 'script')
    except SADError as e:
        ffs.write(' ???' + str(e) + '\n')
        return len(words) - i
    ffs._out_n = getattr(ffs, '_out_n', 1) + 1
    ffs.write('Out[%d]:= %s\n' % (ffs._out_n, ffs.I.to_str(v, quote=True)))
    return len(words) - i


def cmd_execute(ffs, words, i):
    text = ' '.join(words[i + 1:]).strip()
    try:
        s = ffs.I.run(text, 'script')
    except SADError as e:
        ffs.write(' ???' + str(e) + '\n')
        return len(words) - i
    if isinstance(s, str):
        ffs.run_text(s)
    return len(words) - i


def _unquote(w):
    w = w.strip()
    if len(w) >= 2 and w[0] in (chr(34), chr(39)) and w[-1] == w[0]:
        return w[1:-1]
    return w


def _open_target(ffs, words, i, mode):
    """ファイル名 (引用符つき) またはファイル番号 → (ファイルオブジェクト, 自分で開いたか)"""
    import os
    if i + 1 >= len(words):
        return None, False
    w = _unquote(words[i + 1])
    try:
        num = int(float(w))
        f = ffs.I.files.get(num)
        if f is not None:
            return f, False
    except ValueError:
        pass
    path = w if os.path.isabs(w) else os.path.join(ffs.I.cwd, w)
    return open(path, mode, encoding='utf-8', errors='replace'), True


def cmd_output(ffs, words, i, append=False):
    """Implement `OUTPUT(OUT) file` / `APPEND(APP) file`: redirect FFS display and Print output to a file.
    OUTPUT(OUT) file / APPEND(APP) file: 出力 (FFS の表示と Print) をファイルへ切り替える。CLOSE OUTPUT で前の出力に戻す"""
    name = 'APP_END' if append else 'OUT_PUT'
    f, own = _open_target(ffs, words, i, 'a' if append else 'w')
    if f is None:
        ffs.write(' ???-FFS-Error-?Missing filename for ' + name + ' ' + chr(10))
        return len(words) - i
    st = getattr(ffs, '_out_stack', None)
    if st is None:
        st = ffs._out_stack = []
    st.append((ffs.out, ffs.I.out, f if own else None))
    ffs.out = f
    ffs.I.out = f
    return 2


def _close_output(ffs):
    st = getattr(ffs, '_out_stack', None)
    if not st:
        return
    out, iout, own = st.pop()
    try:
        ffs.I.out.flush()
    except Exception:
        pass
    ffs.out = out
    ffs.I.out = iout
    if own is not None:
        own.close()


def cmd_close(ffs, words, i):
    """CLOSE [INPUT] / CLOSE OUTPUT / TERMINATE [INPUT|OUTPUT]"""
    n = 1
    what = 'IN'
    if i + 1 < len(words) and re.fullmatch(r'(?i)OUT(?:P(?:U(?:T)?)?)?', words[i + 1]):
        what = 'OUT'
        n = 2
    elif i + 1 < len(words) and re.fullmatch(r'(?i)IN(?:P(?:U(?:T)?)?)?', words[i + 1]):
        n = 2
    if what == 'OUT':
        _close_output(ffs)
    return n


def cmd_input(ffs, words, i, rewind=False):
    """Implement `INPUT(IN) file` / `READ file`: execute FFS and SADScript statements from a file, then return to the original input.
    INPUT(IN) file / READ file: ファイルの FFS・SADScript の文を続けて実行 (末尾で元の入力に戻る)"""
    import os
    if i + 1 >= len(words):
        ffs.write(' ???-FFS-Error-?Missing filename for IN_PUT ' + chr(10))
        return len(words) - i
    w = _unquote(words[i + 1])
    try:
        num = int(float(w))
        f = ffs.I.files.get(num)
        text = f.read() if f is not None else ''
    except ValueError:
        path = w if os.path.isabs(w) else os.path.join(ffs.I.cwd, w)
        try:
            text = open(path, encoding='utf-8', errors='replace').read()
        except OSError:
            ffs.write(' ???-FFS-Error-?Cannot open file: ' + w + '' + chr(10))
            return 2
    ffs.run_text(text)
    return 2



class _Twbuf:
    """twbuf.f: 8 桁のタブ位置に語を並べる行バッファ (irmgn=行末の余白)"""

    def __init__(self, ffs, irmgn):
        self.ffs = ffs
        self.buff = ''
        self.irmgn = irmgn

    def word(self, w, unit=''):
        itab = 8
        ilmgn = 1
        while True:
            ip = len(self.buff)
            ip1 = (max(ip, ilmgn - 1) + itab - 1) // itab * itab
            if ip1 >= self.irmgn:
                self.ffs.write(self.buff + chr(10))
                self.buff = ''
                continue
            w1 = w.rstrip()
            l = len(w1) + len(unit) + 1
            if l > 1:
                if ip1 + l > self.irmgn:
                    self.ffs.write(self.buff + chr(10))
                    self.buff = ''
                    continue
                self.buff = self.buff + ' ' * (ip1 + 1 - len(self.buff)) + w1 + unit
            return

    def flush(self):
        if len(self.buff) > 0:
            self.ffs.write(self.buff + chr(10))
        self.buff = ''


def cmd_status(ffs, words, i):
    """Implement `STATUS`: list flags, variables and special variables (STAT_US in tffsa.f).
    STATUS: フラグ、変数、特別変数の一覧 (tffsa.f の STAT_US)。CTIME・DT は gs と一致しない (実行時間)"""
    import time
    from .ffscmd import _FLAGS_DEFAULT, _SINO
    g = ffs.I.globals
    now = time.process_time()
    t0 = getattr(ffs, '_ctime0', None)
    if t0 is None:
        t0 = ffs._ctime0 = 0.0
    t2 = getattr(ffs, '_ctime2', 0.0)
    ffs._ctime2 = now
    tb = _Twbuf(ffs, 130)
    tb.word('CTIME=' + fmt_real(now - t0, 'S8.3'), ' sec')
    tb.word('DT=' + fmt_real(now - t2, 'S8.3'), ' sec')
    tb.flush()
    inv = {v: k for k, v in _SINO.items()}
    for name in list(_FLAGS_DEFAULT)[:53]:
        if ffs.flags.get(name, False):
            tb.word(name)
        else:
            tb.word(inv.get(name) or ('NO' + name))
    tb.flush()

    def num(k, d=0.0):
        v = g.get(k, d)
        return float(v) if isinstance(v, (int, float)) else d
    tb.word('CHARGE=' + fmt_real(num('CHARGE', 1.0), 'S7.4'))
    tb.word('MASS=' + fmt_real(num('MASS', 0.51099895e6) / 1e6, 'S10.7'), ' MeV')
    tb.word('MOMENTUM=' + fmt_real(ffs.momentum() / 1e9, 'S10.6'), ' GeV')
    tb.word('PBUNCH=' + fmt_real(num('PBUNCH', 1e10), 'S9.6'))
    tb.word('DP=' + fmt_real(ffs.dp(), 'S8.5'))
    tb.word('DP0=' + fmt_real(ffs.dp0(), 'S8.5'))
    tb.word('FSHIFT=' + fmt_real(num('FSHIFT'), 'S9.6'))
    tb.word('MINCOUP=' + fmt_real(num('MINCOUP', 1e-2), 'S9.6'))
    tb.word('EMITX=' + fmt_real(num('EMITX'), 'S9.6'), ' m')
    tb.word('EMITY=' + fmt_real(num('EMITY'), 'S9.6'), ' m')
    tb.word('EMITZ=' + fmt_real(num('EMITZ'), 'S9.6'), ' m')
    tb.word('SIGZ=' + fmt_real(num('SIGZ'), 'S9.6'), ' m')
    tb.word('SIGE=' + fmt_real(num('SIGE'), 'S9.6'))
    tb.word('XIX=' + fmt_real(num('XIX'), 'S7.4'))
    tb.word('XIY=' + fmt_real(num('XIY'), 'S7.4'))
    tb.word('NBUNCH=' + fmt_real(num('NBUNCH', 1.0), ''))
    tb.word('NP=' + fmt_real(num('NP', 1000.0), 'S8.1'))
    tb.word('GCUT=' + fmt_real(ffs.gcut() if hasattr(ffs, 'gcut') else num('GCUT', 1e300), ''))
    tb.word('CODCONV=' + fmt_real(num('CODCONV', 1.0), 'S8.3'))
    from . import compat as _cp
    tb.word('EMITDIV=' + fmt_real(num('EMITDIV', _cp.emitdiv_default()), ''))
    tb.word('EMITDIVB=' + fmt_real(num('EMITDIVB', 1.0), ''))
    tb.word('EMITDIVQ=' + fmt_real(num('EMITDIVQ', 1.0), ''))
    tb.word('EMITDIVS=' + fmt_real(num('EMITDIVS', 1.0), ''))
    tb.word('CONVERGENCE=' + fmt_real(num('CONVERGENCE', 1e-9), ''))
    tb.word('MAXITERATION=' + fmt_real(num('MAXITERATION', 40.0), ''))
    tb.word('NPARA=' + fmt_real(num('NPARA', 4.0), ''))
    tb.flush()
    first = ffs.elname(1)
    last = ffs.elname(ffs.nlat)
    tb.word('FIT ' + last + ' ' + last)
    tb.word('DISP_LAY ' + first + ' ' + last)
    tb.word('ORG ' + first)
    tb.word('MEA_SURE ' + last)
    tb.flush()
    return 1


COMMANDS = {
    'STATUS': cmd_status, 'STAT': cmd_status,
    'OUTPUT': cmd_output, 'OUT': cmd_output, 'OUTP': cmd_output,
    'APPEND': lambda f, w, i: cmd_output(f, w, i, True), 'APP': lambda f, w, i: cmd_output(f, w, i, True),
    'CLOSE': cmd_close, 'CLO': cmd_close, 'TERMINATE': cmd_close, 'TERM': cmd_close,
    'INPUT': cmd_input, 'IN': cmd_input, 'READ': lambda f, w, i: cmd_input(f, w, i, True),
    'SUSPEND': lambda f, w, i: 1, 'SUSP': lambda f, w, i: 1, 'RESUME': lambda f, w, i: 1, 'RES': lambda f, w, i: 1,
    'VARIABLES': cmd_variables, 'VARS': cmd_variables, 'VAR': cmd_variables,
    'ATTRIBUTE': cmd_attribute, 'ATTR': cmd_attribute, 'ATT': cmd_attribute,
    'MATRIX': cmd_matrix, 'MAT': cmd_matrix,
    'SHOW': cmd_show,
    'VARY': cmd_vary,
    'RECOVER': cmd_recover, 'REC': cmd_recover,
    'REVERSE': cmd_reverse, 'REV': cmd_reverse,
    'PRINT': cmd_print, 'PRI': cmd_print,
    'EXECUTE': cmd_execute, 'EXEC': cmd_execute,
}


def run(ffs, W, words, i):
    f = COMMANDS.get(W)
    if f is None:
        return None
    return f(ffs, words, i)
