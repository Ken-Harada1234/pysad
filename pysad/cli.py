# -*- coding: utf-8 -*-
"""コマンドライン: `pysad-run` (SAD の台本を pysad で実行する)、`pysad-share` (文書・例のある場所を表示する)、`pysad-warmup` (numba のカーネルを先に変換する)。

pip で入れると、文書 (docs/) と例 (example/) はパッケージとは別に `<prefix>/share/pysad/` に入る
(`pysad-share` がその場所を表示する)。リポジトリを clone して使うときは、リポジトリの直下の docs/ と example/ を返す。
"""
import os
import sys
import sysconfig


def share_dir():
    """文書 (docs/) と例 (example/) を含むフォルダ。clone したリポジトリならリポジトリの直下。pip で入れたなら <prefix>/share/pysad。
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
    """`pysad-run lattice_file [script_file ...]`: SAD の `gs` に複数のファイルを連結して渡すのと同じ。"""
    from .ffscmd import run_sad
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ('-h', '--help'):
        print('使い方: pysad-run lattice_file [script_file ...]   (SAD の gs に複数ファイルを連結して渡すのと同じ)')
        return 0 if args else 2
    run_sad(args)
    return 0


def share_main(argv=None):
    """`pysad-share`: 文書 (docs/) と例 (example/) のあるフォルダを表示する。"""
    d = share_dir()
    if d is None:
        print('文書・例のフォルダが見つかりません (pip で入れた場合は <prefix>/share/pysad/ にあるはずです)', file=sys.stderr)
        return 1
    print(d)
    return 0


def warmup_main(argv=None):
    """`pysad-warmup [--check | --clear] [--quiet]`: numba のカーネルを事前に機械語へ変換してキャッシュを作る (最初の実行の待ち時間をなくす)。
    `--check` はキャッシュが有効か調べ、`--clear` はキャッシュを消す。詳細は `pysad/warmup.py`。"""
    from .warmup import main
    return main(argv)
