# -*- coding: utf-8 -*-
"""v07_14 の六極の検討 (R7) の共通部分: ラティスの読み込み、六極の組 (K2L) の読み込み、リングの鏡面対称の軸、追跡の設定。
リポジトリのどこから実行してもよい。ラティスは既定で example/ddba/v07_14/lattice/v07_14_5GeV.sad (環境変数 V07_14_LATTICE で変えられる)。
六極の組は sextupole_sets.json (sxt_extract_sets.py が元のラティス定義から作った K2L の表。元のラティスがなくても使える)。
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
sys.path.insert(0, ROOT)
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
LATTICE = os.environ.get('V07_14_LATTICE', os.path.join(ROOT, 'example', 'ddba', 'v07_14', 'lattice', 'v07_14_5GeV.sad'))
SETS_JSON = os.path.join(HERE, 'sextupole_sets.json')
OUT = os.path.join(HERE, 'results')
P0 = 5e9           # 運動量 [eV] (SAD の K は運動量で規格化済み。2.5 GeV の版も光学は同じ)


def setup_pysad():
    from pysad import compat
    compat.set_mode('v2')


def load_source(path=None):
    return open(path or LATTICE, encoding='utf-8', errors='replace').read()


def load_lattice(ring='RING', path=None, src=None):
    """(ソース文字列, Lattice)。ring は LINE の名前 (RING = RINGSP、RINGINJ、RINGSYM)"""
    setup_pysad()
    from pysad.parser import load_sad_string
    src = src if src is not None else load_source(path)
    return src, load_sad_string(src).build_lattice(ring)


def load_sets():
    """{'family_names': [...], 'sets': {名前: {'k': [...], 'label': ..., 'origin': ...}}}"""
    return json.load(open(SETS_JSON, encoding='utf-8'))


def mirror_axes(lat):
    """素子の名前の列 (ドリフト・偏向磁石・四極・六極だけ) が回文になる軸の素子の添字 (リングの鏡面対称の軸)。無ければ []"""
    idx = [i for i, e in enumerate(lat.elements) if type(e).__name__ in ('Sextupole', 'Bend', 'Quadrupole', 'Drift')]
    seq = [lat.elements[i].name for i in idx]
    n = len(seq)
    out = []
    for k in range(n):
        if all(seq[(k + j) % n] == seq[(k - 1 - j) % n] for j in range(n // 2)):
            out.append(idx[k])
    return out


def make_track_setup(src, rad=False, rfsw=False, alost=0.05, ring='RING'):
    setup_pysad()
    from pysad.trksetup import TrackSetup
    st, ffs = TrackSetup.from_source(src, ring, rad=rad, rfsw=rfsw, alost=alost)
    return st


def lower_priority():
    try:
        from pysad.trksetup import lower_priority as lp
        lp()
    except Exception:
        pass
