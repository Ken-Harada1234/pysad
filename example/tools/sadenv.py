# -*- coding: utf-8 -*-
"""gs (K64 版 SAD 本体) と SAD の Packages の場所を環境変数から決める (example のスクリプト共通)。

環境変数:
  PYSAD_SAD_PACKAGES  K64 版 SAD の Packages のフォルダ (Windows の表記、例: <SAD のフォルダ>/Packages)。必須。
  PYSAD_GS            WSL から見た gs の絶対パス。省略すると PYSAD_SAD_PACKAGES の隣の bin/gs を WSL の表記 (/mnt/<ドライブ>/...) にしたもの。
"""
import os
import re


def packages():
    p = os.environ.get('PYSAD_SAD_PACKAGES')
    if not p:
        raise SystemExit('環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定してください')
    return p.replace(chr(92), '/')


def gs_wsl():
    g = os.environ.get('PYSAD_GS')
    if g:
        return g
    root = packages().rstrip('/').rsplit('/', 1)[0]
    m = re.match(r'([A-Za-z]):(.*)', root)
    if m:
        root = '/mnt/%s%s' % (m.group(1).lower(), m.group(2))
    return root + '/bin/gs'
