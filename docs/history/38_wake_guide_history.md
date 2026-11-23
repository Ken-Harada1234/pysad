# 38_wake_guide の履歴 (整理で本文から取り除いた記述)

本文 (`docs/38_wake_guide.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 128 行 (節: 5a. wake の重み (v2 モードで直した K64 の誤り))

> ## 5a. wake の重み (v2 モードで直した K64 の誤り)

## 元 166 行 (節: 9. 限界・未確認)

> - TRPT で最初から LOSSAMPL の外にいる粒子があると、元の gs は落ちる (B33)。直した SAD (`~/sad_dbg`) との出力座標の差は **解消** (2026-11-16 確認)。原因は K64 の `tapert` が落とす印に `kptbl(i,4)=l_track` を書くこと: `l_track` はセッション変数で、セッション最初の TrackParticles では 0 (CALC をしても変わらない)。最初の tapert では印が 0 =「生きている」になり落ちず (座標は 3·LOSSAMPL に頭打ちされて追跡を続け、後の tapert で落ちる)、2 回目の TrackParticles からは `l_track≠0` で最初の要素で落ちる。pysad は以前は常に最初の要素で落としていたので `TrackRun.l_track0`・`k_apert` の lcur==0 で再現し、gs の 1 回目・2 回目の呼び出しと直した SAD (LOSSAMPL=1e-4、3 粒子、NORAD、TRPT、2 周、15 値) と全桁一致 ([`test_small_items.py`](../pysad/tests/test_small_items.py))。元の K64 はここで落ちるので、元の K64 に対する基準は無い。

