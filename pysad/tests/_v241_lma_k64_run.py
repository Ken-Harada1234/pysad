# -*- coding: utf-8 -*-
"""v2.4.1: 保存済みの台本 (scripts/v241_lma_*.sad) を pysad の k64 モードで走らせ、`.pysadk64.txt` に書く (gs の再現の確認用)。
使い方: py -3.12 pysad/tests/_v241_lma_k64_run.py 台本のタグ (拡張子なし) [v2 も `--v2`: `.pysad.txt` を作り直す]"""
import io
import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _v241_lma_tools as T
from pysad import compat, ffscmd

tag = sys.argv[1]
mode = 'v2' if '--v2' in sys.argv else 'k64'
compat.set_mode(mode)
txt = open(os.path.join(T.SP, tag + '.sad'), encoding='utf-8').read()
buf = io.StringIO()
t0 = time.time()
ffscmd.run_sad_string(txt, out=buf)
open(os.path.join(T.SP, tag + ('.pysad.txt' if mode == 'v2' else '.pysadk64.txt')), 'w', encoding='utf-8').write(buf.getvalue())
print('%s %s: %.0f s' % (tag, mode, time.time() - t0))
