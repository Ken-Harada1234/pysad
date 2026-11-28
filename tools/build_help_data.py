# -*- coding: utf-8 -*-
"""pysad ヘルプ (pysad/sadhelp.py) の元データ `pysad/helpdata/sadhelp_index.json` を、docs/34 と docs/35 の表から作る。

    py -3.12 tools/build_help_data.py            # pysad/helpdata/sadhelp_index.json を作り直す
    py -3.12 tools/build_help_data.py --check    # 作り直した結果が今のファイルと同じかだけ確かめる (書かない)

元:
- docs/34_sad_help_index.md  SADHelp (1.4.4.2k64 Updated: 10/19/2026) から機械的に抽出した、名前と説明の最初の行。
- docs/35_pysad_vs_sadhelp.md  SAD 側の定義場所と、pysad の状態 (実装・未対応・無視 など) と注記。
- PYSAD_ONLY (下)  SAD にない pysad だけの関数。手で書いた表。
- pysad/sad_copy/Packages/*.n  同梱の SAD の .n で定義されている関数・Tkinter のウィジェット (名前・署名・定義されているファイル。コメントがあれば 1 行)。
  名前が $ を含むものと、小文字で始まるものは、利用者が呼ぶものではないとみなして省く。

SADHelp の全文は同梱しない (再配布してよいかが未確認)。最初の行は、すでに docs/34 としてリポジトリで公開している。
"""
import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'pysad', 'helpdata', 'sadhelp_index.json')
SADHELP_VERSION = '1.4.4.2k64 Updated: 10/19/2026'

# SAD にない pysad だけの関数 (CLAUDE.md の一覧と docs/24〜26)。説明は日本語・英語の 1 行ずつ
PYSAD_ONLY = [
    ('FrequencyMap', 'Frequency map (NAFF tunes over an x-y grid of initial amplitudes).', '周波数マップ (x-y の初期振幅の格子で、NAFF のチューンを求める)。', 'docs/26'),
    ('MomentumAcceptance', 'Local momentum acceptance along the ring (tracking with RF and aperture).', 'リングに沿った局所運動量アクセプタンス (RF とアパーチャを入れた追跡)。', 'docs/26'),
    ('TouschekLifetimeLMA', 'Touschek lifetime using the local momentum acceptance.', '局所運動量アクセプタンスを使った Touschek 寿命。', 'docs/26'),
    ('ResonanceDrivingTerms', 'Resonance driving terms (sextupoles and bends) from the optics.', 'ラティスの光学から求める共鳴駆動項 (6 極と偏向磁石)。', 'docs/25'),
    ('AmplitudeTuneShift', 'Amplitude-dependent tune shift (ADTS) from the optics.', 'ラティスの光学から求める振幅依存チューンシフト (ADTS)。', 'docs/25'),
    ('ChromaticityTerms', 'Chromaticity terms from the optics.', 'ラティスの光学から求めるクロマティシティの項。', 'docs/25'),
    ('NormalFormTune', 'Power series of the tunes in the amplitudes and delta (numerical normal form of the one-turn map).', '振幅と δ に依るチューン ν(Jx, Jy, δ) のべき級数 (1 周の写像の正規形)。', 'docs/25'),
    # 手法の名前 (関数ではない)。5 番目の要素は種別、6 番目は注記
    ('NAFF', 'Numerical Analysis of Fundamental Frequencies (Laskar): the tune of a turn-by-turn signal from the maximum of a windowed inner product.', '基本周波数の数値解析 (Laskar)。周回ごとの信号のチューンを、窓つきの内積の最大から求める手法。pysad/naff.py。', 'docs/24',
     'method', 'pysad の関数ではなく手法の名前 (FrequencyMap[]・ADTS の追跡の検証が使う)'),
]

_KIND_BY_34 = {
    'elements': 'element', 'flags': 'flag', 'expression': 'expression', 'beam-line': 'beam-line', 'character-string': 'character-string',
    'data-structure': 'data-structure', 'dynamics': 'dynamics', 'pattern': 'pattern', 'set-value-of-element': 'set-value-of-element',
}
_DOCS = os.path.join(ROOT, 'docs')


def _cells(line):
    return [c for c in line.strip().strip('|').split('|')]


def _norm(raw):
    """`CALCULATE(CAL)` → (表示名, 鍵の並び)。バッククォートと両端の空白を外し、括弧の中は略形の鍵にする"""
    s = raw.strip().strip('`').strip()
    m = re.match(r'^([^()\s]+)\(([^)]*)\)$', s)
    if m:
        return s, [m.group(1), m.group(2)]
    return s, [s]


def parse34(path):
    ents = []
    sec = sub = ''
    parent = ''
    for line in open(path, encoding='utf-8').read().split('\n'):
        if line.startswith('## '):
            sec, sub, parent = line[3:].strip(), '', ''
        elif line.startswith('### '):
            sub, parent = line[4:].strip(), ''
        elif line.startswith('|') and not line.startswith('|---'):
            c = _cells(line)
            if len(c) < 2 or c[0].strip() in ('名前', 'コマンド', '章', '関数', '変数'):
                continue
            indent = len(c[0]) - len(c[0].lstrip(' ')) - 1       # 表の左端の空白 1 つは数えない
            name, keys = _norm(c[0])
            if not name:
                continue
            if sec.startswith('1.'):
                kind, cat = 'ffs-command', 'FFS commands'
            elif sec.startswith('functions'):
                kind, cat = 'function', re.sub(r'\s*\(\d+\)$', '', sub) or 'functions'
            elif sec.startswith('special-variables'):
                kind, cat = 'special-variable', 'special-variables'
            elif sec.startswith('その他の章'):
                if sub:
                    base = re.sub(r'\s*の下位項目.*$', '', sub)
                    kind, cat = _KIND_BY_34.get(base, base), base
                    if base == 'elements':
                        if indent == 0:
                            parent = name
                        else:
                            kind = 'element-keyword'
                else:
                    kind, cat = 'chapter', 'chapters'
            else:
                continue
            summ = c[-1].strip()
            if kind == 'chapter':
                summ = c[2].strip() if len(c) >= 3 else summ
            ents.append({'name': name, 'keys': keys, 'kind': kind, 'category': cat, 'summary': summ,
                         'parent': parent if kind == 'element-keyword' else ''})
    return ents


def parse35(path):
    """(kind の組, 名前) → {sad_source, status, note}。要素の型の表は keywords・missing も"""
    res = {}
    sec = ''
    for line in open(path, encoding='utf-8').read().split('\n'):
        if line.startswith('## '):
            sec = line[3:].strip()
            continue
        if not line.startswith('| `') or '---' in line[:6]:
            continue
        c = [x.strip() for x in _cells(line)]
        if sec.startswith('5.') and len(c) >= 4:
            name, _ = _norm(c[0])
            res[('element', name.upper())] = {'pysad_type': c[1], 'keywords': c[2], 'missing_keywords': c[3], '_name': name, '_keys': [name]}
            continue
        if len(c) < 3:
            continue
        if sec.startswith('1.'):
            kind = 'ffs-command'
        elif sec.startswith('2.'):
            kind = 'function'
        elif sec.startswith('3.'):
            kind = 'special-variable'
        elif sec.startswith('4.'):
            kind = 'flag'
        else:
            continue
        name, keys = _norm(c[0])
        m = re.match(r'^([^\s(（、。]+)\s*[(（]?(.*?)[)）]?$', c[2]) if c[2] else None
        status = m.group(1) if m else c[2]
        note = m.group(2).strip() if m else ''
        if note.endswith(')') and note.count('(') < note.count(')'):
            note = note[:-1]
        d = {'sad_source': c[1], 'status': status, 'note': note, '_name': name, '_keys': keys}
        for k in keys:
            res.setdefault((kind, k.upper()), d)
    return res


def doc_links(note):
    """注記の中の `docs/NN`・`docs/NN_name` を、docs の実在ファイル名に直す"""
    out = []
    names = sorted(f for f in os.listdir(_DOCS) if f.endswith('.md'))
    for m in re.finditer(r'docs/(\d\d)(?:[_\w.-]*)', note or ''):
        nn = m.group(1)
        for f in names:
            if f.startswith(nn + '_'):
                if f not in out:
                    out.append(f)
                break
    return out


SAD_COPY = os.path.join(ROOT, 'pysad', 'sad_copy')
_SKIP_N = ('HelpMessages.n',)                       # ヘルプの本文 (全文は sadhelp.py が実行時に読む)
_TK_FILES = ('Tkinter.n', 'TkinterCore.n', 'KBFrame.n', 'KBMainFrame.n', 'TabFrame.n', 'TkPhoto.n', 'Canvas.n', 'Canvas3D.n', 'CanvasCM.n',
             'CANVASDRAW.n', 'CheckEntryList.n', 'CursorEntry.n', 'FileDialog.n', 'KBFFileSelectionPanel.n', 'OptionDialog.n', 'SADTerminal.n')


def scan_n():
    """Packages/*.n の最上位 (行頭) の定義 `Name[...] := ...`・`Name = ...` を集める。名前 → [(ファイル, 署名, 直前のコメント)]"""
    defs = {}
    d = os.path.join(SAD_COPY, 'Packages')
    for fn in sorted(os.listdir(d)):
        if not fn.endswith('.n') or fn in _SKIP_N:
            continue
        lines = open(os.path.join(d, fn), encoding='utf-8', errors='replace').read().split('\n')
        for i, l in enumerate(lines):
            m = re.match(r'^([A-Z][A-Za-z0-9]+)\s*(\[|:=|=(?!=))', l)
            if not m:
                continue
            name = m.group(1)
            if m.group(2) == '[':
                depth, sig = 0, ''
                for k in range(len(name), len(l)):
                    if l[k] == '[':
                        depth += 1
                    elif l[k] == ']':
                        depth -= 1
                        if depth == 0:
                            sig = l[:k + 1]
                            break
                rest = l[len(sig):].lstrip()
                if not sig or not (rest.startswith(':=') or rest.startswith('=') or rest.startswith('/;')):
                    continue                      # 定義ではなく、式の先頭が大文字の呼び出し
            else:
                sig = name
            cm = ''
            if i > 0:
                mm = re.match(r'^\(\*\s*(.*?)\s*\*\)\s*$', lines[i - 1])
                if mm:
                    cm = mm.group(1)
            mm = re.search(r'\(\*\s*(.*?)\s*\*\)\s*$', l)
            if mm and not cm:
                cm = mm.group(1)
            defs.setdefault(name, []).append((fn, re.sub(r'\s+', ' ', sig)[:90], cm[:100]))
    return defs


def tk_widgets():
    """Tkinter.n の `Widgets=Join[ {Frame, Window, ...}, (* tkFrame.c *) ...` から、ウィジェットのクラス名と、元の Tk の C ファイル名"""
    path = os.path.join(SAD_COPY, 'Packages', 'Tkinter.n')
    out = {}
    on = False
    for l in open(path, encoding='utf-8', errors='replace').read().split('\n'):
        if l.startswith('Widgets=Join['):
            on = True
            continue
        if on:
            m = re.match(r'^\s*\{([^}]*)\}\s*,?\s*(?:\(\*\s*(.*?)\s*\*\))?', l)
            if not m:
                if l.strip().startswith(']') or l.strip() == '':
                    if l.strip().startswith(']'):
                        break
                    continue
                break
            for n in [x.strip() for x in m.group(1).split(',')]:
                if re.match(r'^[A-Za-z][A-Za-z0-9]*$', n):
                    out[n] = (m.group(2) or '').strip()
    return out


def add_n_packages(entries):
    """同梱の .n の関数・Tkinter のウィジェットを、既存の項目に足し、無い名前は新しい項目にする"""
    by_key = {}
    for e in entries:
        for k in e['keys'] + [e['name']]:
            by_key.setdefault(k.upper(), []).append(e)
    defs = scan_n()
    n_entries = {}
    for name, items in sorted(defs.items(), key=lambda x: x[0].lower()):
        files = []
        sigs = []
        com = ''
        for fn, sig, cm in items:
            if fn not in files:
                files.append(fn)
            if sig and sig not in sigs:
                sigs.append(sig)
            com = com or cm
        tk = any(f in _TK_FILES for f in files)
        info = {'n_files': files, 'signature': sigs[:3]}
        same = [e for e in by_key.get(name.upper(), []) if e['kind'] == 'function']       # SADHelp の関数と同じ名前 (FFS コマンドの FIT は別物)
        if same:
            for e in same:
                e.update(info)
            continue
        e = {'name': name, 'keys': [name], 'kind': 'package-function', 'category': files[0], 'summary': com, 'status': 'SAD の .n',
             'note': 'pysad が SAD 自身の .n を解釈して実行する (同梱のコピーか PYSAD_SAD_PACKAGES)',
             'docs': doc_links('docs/20' if tk else 'docs/06')}
        e.update(info)
        entries.append(e)
        n_entries[name.upper()] = e
    for name, cfile in tk_widgets().items():
        summ = 'Tk widget class' + (' (Tk C source %s)' % cfile if cfile else '')
        e = n_entries.get(name.upper())
        if e is None:
            if name.upper() in by_key:                 # 既存の SADHelp の項目と同じ名前 (Canvas など): 注記だけ足す
                for x in by_key[name.upper()]:
                    x.setdefault('tk_widget', summ)
                continue
            e = {'name': name, 'keys': [name], 'kind': 'tk-widget', 'category': 'Tkinter widgets', 'summary': summ, 'status': 'SAD の .n',
                 'note': 'pysad が SAD 自身の .n を解釈して実行する (同梱のコピーか PYSAD_SAD_PACKAGES)', 'docs': doc_links('docs/20')}
            entries.append(e)
        else:
            e['kind'], e['category'] = 'tk-widget', 'Tkinter widgets'
            e['summary'] = e['summary'] or summ
    return entries


def build():
    e34 = parse34(os.path.join(_DOCS, '34_sad_help_index.md'))
    d35 = parse35(os.path.join(_DOCS, '35_pysad_vs_sadhelp.md'))
    entries = []
    used = set()
    entries = []
    for e in e34:
        info = None
        for key in e['keys']:
            info = d35.get((e['kind'], key.upper()))
            if info:
                break
        if info:
            used.add(id(info))
            e.update({kk: vv for kk, vv in info.items() if vv != '' and not kk.startswith('_')})
            links = doc_links(info.get('note', ''))
            if links:
                e['docs'] = links
        e.setdefault('summary', '')
        entries.append(e)
    left = {}
    for (kind, _k), info in d35.items():           # docs/34 に無い名前 (docs/35 にだけある)
        if id(info) not in used and id(info) not in left:
            left[id(info)] = (kind, info)
    for kind, info in left.values():
        e = {'name': info['_name'], 'keys': info['_keys'], 'kind': kind, 'category': '', 'summary': '', 'parent': ''}
        e.update({kk: vv for kk, vv in info.items() if vv != '' and not kk.startswith('_')})
        links = doc_links(info.get('note', ''))
        if links:
            e['docs'] = links
        entries.append(e)
    for name, en, ja, doc, *more in PYSAD_ONLY:
        kind, note = more if more else ('function', 'SAD (gs) には無い pysad の関数')
        entries.append({'name': name, 'keys': [name], 'kind': kind, 'category': 'pysad only', 'summary': en, 'summary_ja': ja,
                        'parent': '', 'status': 'pysad のみ', 'note': note, 'docs': doc_links(doc)})
    entries = add_n_packages(entries)
    for e in entries:
        for k in [kk for kk in ('parent', 'category') if not e.get(kk)]:
            e.pop(k, None)
    return {'sadhelp_version': SADHELP_VERSION,
            'generated_from': ['docs/34_sad_help_index.md', 'docs/35_pysad_vs_sadhelp.md', 'pysad/sad_copy/Packages/*.n'],
            'entries': entries}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true', help='今のファイルと同じかだけ確かめる')
    a = ap.parse_args(argv)
    data = build()
    text = json.dumps(data, ensure_ascii=False, indent=1) + '\n'
    if a.check:
        cur = open(OUT, encoding='utf-8').read() if os.path.exists(OUT) else ''
        print('same' if cur == text else 'DIFFERENT')
        return 0 if cur == text else 1
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, 'w', encoding='utf-8', newline='\n').write(text)
    kinds = {}
    for e in data['entries']:
        kinds[e['kind']] = kinds.get(e['kind'], 0) + 1
    print('entries', len(data['entries']), kinds, 'bytes', len(text.encode('utf-8')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
