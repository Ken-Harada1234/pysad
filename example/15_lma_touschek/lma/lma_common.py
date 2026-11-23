# -*- coding: utf-8 -*-
"""LMA の例の共通部分: ラティスの読み込みと位置の選び方 (ASTRID2・KEK-LS・MPW 版 ASTRID2)。リポジトリ直下から実行する。"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), ring='RING', aper=0.010, ds=0.0,
                    label='ASTRID2 (0.58 GeV、45.7 m、6 セル)'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), ring='RING', aper=0.010, ds=0.5,
                  label='kek-ls_v03_12.lat (3 GeV、603 m、20 セル)'),
    'mpw': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), ring='RINGMPW', aper=0.010, ds=0.0,
                label='ASTRID2 + MPW (RINGMPW)'),
}


def load(name, **kw):
    """(TrackSetup, FFS)。v2 モード、RF・放射つき、LOSSAMPL = ±10 mm (例の値)"""
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    return TrackSetup.from_source(src, cs['ring'], rad=True, rfsw=True, alost=kw.pop('alost', cs['aper']), **kw)


def elegant_kekls_locs():
    """elegant_lma.py の elegant (kek-ls、1500 周) の位置 (pysad の loc) と δ± の表: {loc: (s, dp+, dp-)}"""
    out = {}
    for ln in open(os.path.join(HERE, '..', 'touschek_elegant', 'elegant_kekls.lma.txt'), encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        out[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
    return out


def elegant_astrid2_table():
    out = {}
    for ln in open(os.path.join(HERE, '..', 'touschek_elegant', 'elegant_astrid2.lma.txt'), encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        out[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
    return out
