# -*- coding: utf-8 -*-
"""KEK-LS (v03_68、3 GeV、570.72 m、RING = 標準セル 19 + 入射セル) に誤差を与え、COD を補正し、補正後の状態で MomentumAcceptance[] を走らせる
(SAD の台本の経路: 補正後の状態を TrackSetup に渡す)。台本の本体は lma_corrected_kekls_v03_68.template.sad (example/kek_ls/errors/KEKLS_correct_cod.sad の誤差・補正の部分に、
補正後の LMA を足したもの。ラティスの定義は先頭に貼ってある)。
使い方: PYSAD_SAD_PACKAGES=... PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_corrected_run.py 種 ctrl niter ae fe te [workers [turns [lat]]]
  種: SeedRandom の種。ctrl: 0 = NOEXPAND のまま LMA (補正後の成分ごとの K0・位置ずれを保つ)、1 = LMA の前に `cell calc expand` (既定の CALEXP: 代表の値に揃える対照)。
  niter: 補正の繰り返し (0 = 補正なし)。ae・fe・te: 位置ずれ [m]・強さの相対誤差・回転 [rad] (標準 4E-5・4E-4・1E-4。誤差の種類と GCUT = 1 は KEKLS_correct_cod と同じ)。
  lat: v03_68 (既定。example/kek_ls、3224 要素・570.72 m) か v03_12 (pysad/tests/kek-ls_v03_12.lat、3563 要素・603 m。Distance = 150 m)。
  rad (第 10 引数): 1 なら LMA の前に `RAD` (放射減衰・励起つきの追跡と、放射を含む閉軌道・同期位相。RFSW は FFS の既定で ON)、既定の 0 は RAD OFF (同期位相 0)。
  dist (第 11 引数): LMA の位置の間隔 [m] (既定 v03_68 143・v03_12 150)。
  workers: LMA の子プロセス数 (既定 4)。turns (省略可、第 8 引数): LMA の周回数 (既定は自動 = シンクロトロン振動 2 周期 = 428)。Distance = 143 m ごとの 4 位置、Aperture ±10 mm、RFSW・RAD ON、FLUC なし、Turns は自動 (シンクロトロン振動 2 周期)。
例: 補正あり lma_corrected_run.py 111 0 10 4E-5 4E-4 1E-4 / 対照 lma_corrected_run.py 111 1 10 4E-5 4E-4 1E-4 / 補正なし lma_corrected_run.py 111 0 0 4E-5 4E-4 1E-4 / 誤差なし lma_corrected_run.py 111 0 0 0 0 0"""
import os
import sys
import time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)


def main():
    from pysad.ffscmd import run_sad_string
    a = sys.argv[1:]
    seed, ctrl, niter, ae, fe, te = a[:6]
    wk = a[6] if len(a) > 6 else '4'
    turns = a[7] if len(a) > 7 else None
    lat = a[8] if len(a) > 8 else 'v03_68'
    rad = (a[9] if len(a) > 9 else '0') == '1'
    dist = a[10] if len(a) > 10 else None
    tpl = 'lma_corrected_kekls_v03_68.template.sad' if lat == 'v03_68' else 'lma_corrected_kekls_v03_12.template.sad'
    src = open(os.path.join(HERE, tpl), encoding='utf-8').read()
    if lat != 'v03_68':
        src = src.replace('Dist = 143;', 'Dist = 150;')          # 603 m を 4 位置
    if dist:
        src = src.replace('Dist = 150;', 'Dist = %s;' % dist).replace('Dist = 143;', 'Dist = %s;' % dist)
    for k, v in (('SEEDVAL', seed), ('CTRLVAL', ctrl), ('NITERVAL', niter), ('AEVAL', ae), ('FEVAL', fe), ('TEVAL', te)):
        src = src.replace(k, v)
    src = src.replace('Wk = 4;', 'Wk = %s;' % wk)
    src = src.replace('RADLINE', 'FFS["RAD"];' if rad else '')          # rad = 1: LMA の前に RAD ON (RFSW は FFS の既定で ON)。既定は RAD OFF (同期位相 0)
    if turns:
        src = src.replace('Aperture->0.01]', 'Aperture->0.01, Turns->%s]' % turns)
    t0 = time.time()
    run_sad_string(src, cwd=os.path.join(ROOT, 'example', 'kek_ls', 'errors'))
    print('WALL', time.time() - t0, flush=True)


if __name__ == '__main__':
    main()
