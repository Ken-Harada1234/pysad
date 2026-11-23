# -*- coding: utf-8 -*-
"""リポジトリ内の md を、図入り・目次つきの HTML (md と同じ名前の .html) にまとめて変換する。

    py -3.12 tools/build_html.py                 # 追跡対象の md をすべて変換
    py -3.12 tools/build_html.py docs README.md  # 指定したファイル・フォルダの分だけ
    py -3.12 tools/build_html.py --changed       # git status に出ている md と、HTML より新しい md だけ
    py -3.12 tools/build_html.py --check         # 壊れたリンク (md→md・画像・アンカー) と、索引・入口のリンク先の存在の検査だけ (書き出さない)

終わりに、全文検索の索引 (search_index.json)、検索の画面 (search.html)、入口ページ (index.html) を作る (site_index.py)。
一部の md だけを変換したときは、既存の索引のその分だけを入れ替える。

除外: docs/history/ の下、CHANGELOG.md、LICENSE 類、example/tools/kek_ls_src/ (生成用の断片)。作業ツリーに無い (削除済みの) md は対象にしない。
個別の設定 (折りたたみなど) は tools/html_config.json (リポジトリの根からの相対パス → オプション)。
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import md_to_html as m2h  # noqa: E402
import mathtex  # noqa: E402
import site_index  # noqa: E402

CONFIG = HERE / 'html_config.json'


def _git(*args):
    r = subprocess.run(['git', '-C', str(ROOT)] + list(args), capture_output=True, text=True,
                       encoding='utf-8', check=True)
    return r.stdout


def excluded(rel):
    parts = rel.split('/')
    name = parts[-1].lower()
    return (rel.startswith('docs/history/') or name == 'changelog.md'
            or name.startswith('license') or name.startswith('licence')
            or rel.startswith('example/tools/kek_ls_src/'))


def tracked_md():
    out = _git('-c', 'core.quotepath=off', 'ls-files', '*.md')
    return sorted(l for l in out.splitlines() if l and not excluded(l) and (ROOT / l).is_file())


def changed_md():
    out = _git('-c', 'core.quotepath=off', 'status', '--porcelain', '-uall', '--', '*.md')
    res = set()
    for l in out.splitlines():
        if len(l) < 4 or 'D' in l[:2]:
            continue
        p = l[3:].split(' -> ')[-1].strip('"')
        if p.endswith('.md') and not excluded(p) and (ROOT / p).is_file():
            res.add(p)
    return res


def load_config():
    if CONFIG.exists():
        return {k: v for k, v in json.loads(CONFIG.read_text(encoding='utf-8')).items() if not k.startswith('_')}
    return {}


def main(argv=None):
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding='utf-8', errors='replace')
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('paths', nargs='*', help='md のファイルまたはフォルダ (省略で全件)')
    ap.add_argument('--check', action='store_true', help='壊れた相対リンクの検査だけ (HTML は書かない)')
    ap.add_argument('--changed', action='store_true', help='変更された md (git status と、HTML より新しいもの) だけ')
    ap.add_argument('--require-math', action='store_true',
                    help='Node.js が見つからない、または数式の変換エラーがあるとき、エラー終了にする (既定は警告のみ)')
    ap.add_argument('--include-untracked', action='store_true', help='追跡されていない md も対象に含める')
    a = ap.parse_args(argv)

    if a.require_math and mathtex.find_node() is None:
        print('エラー: Node.js が見つからない (PYSAD_NODE、PATH の node、D:/work/tools/node-*/node.exe の順に探した)')
        return 2
    allmd = tracked_md()
    if a.include_untracked:
        allmd = sorted(set(allmd) | changed_md())
    targets = {(ROOT / p).resolve() for p in allmd}
    cand = allmd
    if a.paths:
        sel = []
        for p in a.paths:
            ap_ = Path(p)
            ap_ = (ap_ if ap_.is_absolute() else Path.cwd() / ap_).resolve()
            try:
                rel = ap_.relative_to(ROOT).as_posix()
            except ValueError:
                print('リポジトリの外なので無視:', p)
                continue
            if ap_.is_dir():
                pre = rel.rstrip('/') + '/' if rel != '.' else ''
                sel += [x for x in allmd if x.startswith(pre)]
            elif ap_.suffix == '.md' and ap_.is_file():
                if excluded(rel):
                    print('除外の対象なので無視:', rel)
                    continue
                sel.append(rel)
                targets.add(ap_)
            else:
                print('md でもフォルダでもないので無視:', p)
        cand = sorted(set(sel))
    if a.changed:
        ch = changed_md()
        cand = [p for p in cand
                if p in ch or not (ROOT / p).with_suffix('.html').exists()
                or (ROOT / p).stat().st_mtime > (ROOT / p).with_suffix('.html').stat().st_mtime]

    cfg = load_config()
    nfile = nfig = nmath = nmerr = ndoc = ncode = 0
    newpages = {}
    warns, brokens, skipped = [], [], []
    for rel in cand:
        opts = dict(cfg.get(rel, {}))
        try:
            res = m2h.convert_file(ROOT / rel, write=not a.check, is_target=targets.__contains__,
                                   embed=not a.check, root=ROOT, **opts)
        except Exception as e:  # noqa: BLE001
            warns.append((rel, '変換に失敗: %r' % e))
            continue
        if res.skipped:
            skipped.append(rel)
        else:
            nfile += 1
            nfig += res.figures
            nmath += res.math_count
            nmerr += res.math_errors
            ndoc += 1 if res.math_count else 0
            ncode += res.code_blocks
            if res.index is not None and site_index.indexed(rel):
                newpages[rel[:-3] + '.html'] = res.index
        warns += [(rel, w) for w in res.warnings]
        brokens += [(rel, k, t) for k, t in dict.fromkeys(res.broken)]
    mode = '検査' if a.check else '変換'
    site_msgs = []
    if a.check:
        for k, t in site_index.check_site(ROOT):
            brokens.append(('索引・入口', k, t))
    else:
        site_rels = {r for r in allmd if site_index.indexed(r)}
        full = not a.paths and not a.changed
        pages = {}
        if not full:
            old = site_index.load_index(ROOT)
            if old:
                pages = {p: d for p, d in site_index.index_to_pages(old).items() if p[:-5] + '.md' in site_rels}
        pages.update(newpages)
        miss = sorted(r for r in site_rels if r[:-3] + '.html' not in pages)
        if miss:
            warns.append(('索引', '索引に入っていない md が %d 件 (例: %s)。全件の変換で入る' % (len(miss), miss[0])))
        msgs, sw = site_index.write_site(ROOT, pages, site_rels)
        site_msgs += msgs
        warns += [('索引・入口', w) for w in sw]
    print('%s: %d 件 (対象 %d 件)、図 %d 枚、数式 %d 個 (%d 件の文書)・変換エラー %d、色分けしたコード %d、警告 %d、壊れたリンク %d、スキップ %d'
          % (mode, nfile, len(cand), nfig, nmath, ndoc, nmerr, ncode, len(warns), len(brokens), len(skipped)))
    for m in site_msgs:
        print(m)
    for rel, w in warns:
        print('警告  %s: %s' % (rel, w))
    for rel, k, t in brokens:
        print('壊れた%s  %s: %s' % (k, rel, t))
    return 1 if ((a.check and brokens) or (a.require_math and nmerr)) else 0


if __name__ == '__main__':
    sys.exit(main())
