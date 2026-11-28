# -*- coding: utf-8 -*-
"""`python -m pysad <command>`: help / run / share / warmup (see pysad/cli.py).
`python -m pysad <command>`: help (SAD の説明を引く)・run (台本の実行)・share (文書の場所)・warmup (numba のキャッシュ)。"""
import sys

from .cli import module_main

if __name__ == '__main__':
    sys.exit(module_main())
