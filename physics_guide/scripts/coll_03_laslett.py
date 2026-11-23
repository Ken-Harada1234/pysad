# -*- coding: utf-8 -*-
r"""WSPAC のチューンシフトと Laslett の式 ((17.9) の離散和) の再確認 (17_collective.md §4.4)。

例 example/05_tracking/flag_WSPAC_tune_shift.sad の前半 (追跡のチューンシフトと、各要素の入口のレンズの和 -(1/4π) Σ β k) を、そのまま走らせる。
リング: 0.5 GeV の電子、8 セル、周長 35.3 m、PBUNCH = 2e10、弱いスキュー四極 QS で結合を入れて垂直のエミッタンスを 0 でなくしてある。
"""
import io
import os
import re
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
from pysad.ffscmd import run_sad_string

path = os.path.join(C.ROOT, 'example', '05_tracking', 'flag_WSPAC_tune_shift.sad')
src = open(path, encoding='utf-8').read()
src = src[:src.index('! SELFCOD: 力の中心')]
buf = io.StringIO()
run_sad_string(src, out=buf)
txt = buf.getvalue()
for line in txt.splitlines():
    if re.search(r'Twiss:|Emittance:|fractional tune|tracking:|beam size|Laslett', line):
        print(line)
m1 = re.search(r'tracking: dnux=\s*(\S+)\s+dnuy=\s*(\S+)', txt)
m2 = re.search(r'Laslett:\s+dnux=\s*(\S+)\s+dnuy=\s*(\S+)', txt)
if m1 and m2:
    a = [float(m1.group(1)), float(m1.group(2))]
    b = [float(m2.group(1)), float(m2.group(2))]
    print('\n追跡/Laslett: Δνx = %.4e / %.4e (比 %.4f)、Δνy = %.4e / %.4e (比 %.4f)' % (a[0], b[0], a[0] / b[0], a[1], b[1], a[1] / b[1]))
