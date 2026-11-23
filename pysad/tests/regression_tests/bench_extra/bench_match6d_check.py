#!/usr/bin/env python3
"""
bench_match6d_check.py — match6d (SAD の解析微分 tffsqu/qdcell の直訳) と match6 の数値微分の照合
  (1) Level 13B: pysad/tests/kek-ls_v03_12.lat USE RING, fit NX 47.60 NY 17.60, MKC BX 1.2 AX 0 BY 7.5 AY 0, free Q*
  (2) esrf/matching 相当: USE ARCA, matching1 の条件 (FitValue/FitFunction は除く), free Q* BPI6 K1 BPI7 K1
      + EX/EPX/EY/EPY/NX/NY 条件を追加して分散の微分も確認
使い方: PYTHONIOENCODING=utf-8 python -u bench_match6d_check.py
"""
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.match6 import Matcher
from pysad.match6d import compare_derivatives, AnalyticMatcher


def check_13b():
    sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
    lat = sc.build_lattice('RING')
    m = Matcher(lat, out=lambda s: None)
    m.fit_clear()
    m.fit('NX', 47.60)
    m.fit('NY', 17.60)
    m.fit('BX', 1.2, loc='MKC')
    m.fit('AX', 0)
    m.fit('BY', 7.5)
    m.fit('AY', 0)
    m.free('Q*')
    t0 = time.time()
    qa, qn, ids, worst = compare_derivatives(m, verbose=True)
    print(f'[13B] 条件 {len(ids)} × 変数 {len(m.variables)}: 最大相対差 {worst:.3e} ({time.time() - t0:.1f} s)')
    return worst


def check_arca():
    sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
    lat = sc.build_lattice('ARCA')
    m = Matcher(lat, out=lambda s: None)
    m.fit_clear()
    m.fit('BY', 3.6, loc='MKS')
    m.fit('NX', 2.38, loc='$$$')
    m.fit('NY', 0.88)
    m.fit('AX', 0, loc='MKC')
    m.fit('AY', 0)
    m.fit('EPX', 0)
    m.fit('BY', 5)
    m.fit('AY', 1.18, loc='MKSFA1')
    m.fit('BY', 5.4)
    m.fit('EX', 0.1504)
    m.fit('BX', 5)
    m.fit('NX', 1.0, loc='MKSFA1.2')
    m.fit('EY', 0, loc='MKSDC1')
    m.fit('EPY', 0)
    m.free('Q*')
    m.free('BPI6', 'K1')
    m.free('BPI7', 'K1')
    t0 = time.time()
    qa, qn, ids, worst = compare_derivatives(m, verbose=True)
    print(f'[ARCA] 条件 {len(ids)} × 変数 {len(m.variables)}: 最大相対差 {worst:.3e} ({time.time() - t0:.1f} s)')
    return worst


def go_13b():
    """Level 13B のマッチングを AnalyticMatcher で実行し SAD (2 反復, Matched 3.3724E-10) と比べる"""
    sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
    lat = sc.build_lattice('RING')
    m = AnalyticMatcher(lat)
    m.fit_clear()
    m.fit('NX', 47.60)
    m.fit('NY', 17.60)
    m.fit('BX', 1.2, loc='MKC')
    m.fit('AX', 0)
    m.fit('BY', 7.5)
    m.fit('AY', 0)
    m.free('Q*')
    m.go()
    print('13B AnalyticMatcher K1:', {k: round(v, 12) for k, v in m.values().items()})


if __name__ == '__main__':
    w1 = check_13b()
    w2 = check_arca()
    print('合格' if max(w1, w2) < 1e-6 else '不一致あり')
    go_13b()
