# -*- coding: utf-8 -*-
"""Command-line entry points `pysad-run`, `pysad-share`, `pysad-warmup` and `pysad-help` of pysad (also `python -m pysad <command>`).
コマンドライン: `pysad-run` (SAD の台本を pysad で実行する)、`pysad-share` (文書・例のある場所を表示する)、`pysad-warmup` (numba のキャッシュを作る)、
`pysad-help` (SAD の関数・コマンドの説明を引く)。`python -m pysad <command>` でも同じことができる。

pip で入れると、文書 (docs/) と例 (example/) はパッケージとは別に `<prefix>/share/pysad/` に入る
(`pysad-share` がその場所を表示する)。リポジトリを clone して使うときは、リポジトリの直下の docs/ と example/ を使う。
"""
import os
import sys
import sysconfig

_HELP_FLAGS = ('-h', '--help')


def share_dir():
    """Return the folder holding docs/ and example/ (repository root or <prefix>/share/pysad), or None if not found.
    文書 (docs/) と例 (example/) を含むフォルダ。clone したリポジトリならリポジトリの直下。pip で入れたなら <prefix>/share/pysad。
    見つからなければ None。"""
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)
    if os.path.isdir(os.path.join(repo, 'docs')) and os.path.isdir(os.path.join(repo, 'example')):
        return repo
    seen = []
    for key in ('data', 'userbase'):
        try:
            base = sysconfig.get_path('data') if key == 'data' else sysconfig.get_config_var('userbase')
        except Exception:
            base = None
        if not base:
            continue
        cand = os.path.join(base, 'share', 'pysad')
        if cand in seen:
            continue
        seen.append(cand)
        if os.path.isdir(os.path.join(cand, 'docs')):
            return cand
    return None


def run_main(argv=None):
    """Run `pysad-run lattice_file [script_file ...]`, equivalent to passing concatenated files to SAD `gs`.
    `pysad-run lattice_file [script_file ...]`: SAD の `gs` に複数のファイルを連結して渡すのと同じ。"""
    from .ffscmd import run_sad
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in _HELP_FLAGS:
        print('usage: pysad-run lattice_file [script_file ...]   (the same as passing the concatenated files to SAD gs)\n'
              '使い方: pysad-run lattice_file [script_file ...]   (SAD の gs に複数ファイルを連結して渡すのと同じ)')
        return 0 if args else 2
    run_sad(args)
    return 0


def share_main(argv=None):
    """Run `pysad-share`: print where the docs, examples and physics guide are (the folder in a clone, otherwise the GitHub and GitHub Pages URLs).
    `pysad-share`: 文書・例・物理ガイドのある場所を表示する。clone したリポジトリならそのフォルダ。pip で入れたなら GitHub と GitHub Pages の URL。"""
    from .banner import PAGES, URL
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in _HELP_FLAGS:
        print('usage: pysad-share   (print the folder of docs/ and example/ in a clone, otherwise the GitHub and GitHub Pages URLs)\n'
              '使い方: pysad-share   (文書 docs/ と例 example/ のあるフォルダを表示する。pip だけのときは GitHub と GitHub Pages の URL)')
        return 0
    d = share_dir()
    if d is not None:
        print(d)
        return 0
    print('文書 (docs/)・例 (example/)・物理ガイド (physics_guide/)・試験 (pysad/tests/) は、pip には入っていません。')
    print(f'  ソース、例、文書、試験: {URL}')
    print(f'  読みやすい HTML (数式つき): {PAGES}   (docs/、example/、physics_guide/ の README.html)')
    return 0


def warmup_main(argv=None):
    """Run `pysad-warmup [--check | --clear] [--quiet]`: precompile and cache the numba kernels.
    `pysad-warmup [--check | --clear] [--quiet]`: numba のカーネルを事前に機械語へ変換してキャッシュを作る (最初の実行を速くする)。
    `--check` はキャッシュが有効か調べ、`--clear` はキャッシュを消す。詳細は `pysad/warmup.py`。"""
    from .warmup import main
    return main(argv)


def help_main(argv=None):
    """Run `pysad-help [NAME] [--search KEYWORD] [--kind KIND] [--full]`: show help for SAD/FFS/SADScript names (see pysad/sadhelp.py).
    `pysad-help [名前] [--search キーワード] [--kind 種類] [--full]`: SAD の FFS コマンド・関数・特別変数・フラグ・要素の説明を出す。"""
    import argparse
    from . import sadhelp
    ap = argparse.ArgumentParser(
        prog='pysad-help', formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Show help for SAD/FFS/SADScript names and whether pysad supports them.\n'
                    'SAD の FFS コマンド・関数・特別変数・フラグ・要素の説明と、pysad で使えるかを出す (引数なしで概要)。',
        epilog='examples / 例:\n  pysad-help Twiss\n  pysad-help CAL\n  pysad-help BEND.AE1\n  pysad-help --search tune\n  pysad-help FIT --full\n\n'
               '--full needs your own SADHelp.html (set SADHELP_HTML). 全文は手元の SADHelp.html が要る (環境変数 SADHELP_HTML)。')
    ap.add_argument('name', nargs='*', help='name to look up (Twiss, CAL, Twiss[], BEND.AE1 ...). 引く名前')
    ap.add_argument('-s', '--search', metavar='KEYWORD', help='search names, first lines and notes. 名前・最初の行・注記から探す')
    ap.add_argument('-k', '--kind', metavar='KIND', help='ffs-command, function, special-variable, flag, element, ... 種類で絞る')
    ap.add_argument('-f', '--full', action='store_true', help='also print the full SADHelp text (your own SADHelp.html). 手元の SADHelp.html の全文も出す')
    a = ap.parse_args(list(sys.argv[1:] if argv is None else argv))
    name = ' '.join(a.name) if a.name else None
    sadhelp.sad_help(name, a.search, a.full, a.kind)
    return 0


COMMANDS = {
    'help': (help_main, 'SAD の関数・コマンドの説明を引く / look up SAD names'),
    'run': (run_main, 'SAD の台本を実行する / run SAD files'),
    'share': (share_main, '文書・例のある場所を表示する / where docs and examples are'),
    'warmup': (warmup_main, 'numba のキャッシュを作る / precompile numba kernels'),
}


def module_main(argv=None):
    """`python -m pysad <command> ...`: dispatch to help / run / share / warmup (no argument: list the commands).
    `python -m pysad <command> ...`: help・run・share・warmup に振り分ける (引数なしでコマンドの一覧)。"""
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in _HELP_FLAGS or args[0] not in COMMANDS:
        from . import __version__
        unknown = bool(args) and args[0] not in _HELP_FLAGS
        print('pysad %s: python -m pysad <command> [...]\n' % __version__ + 'commands / コマンド:')
        for k, (_f, d) in COMMANDS.items():
            print('  %-7s %s' % (k, d))
        print('例 / example: python -m pysad help Twiss')
        if unknown:
            print('不明なコマンド: %s' % args[0])
        return 2 if (unknown or not args) else 0
    return COMMANDS[args[0]][0](args[1:])
