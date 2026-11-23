# 17_implementation_summary.md から取り除いた履歴の記述

17_implementation_summary.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 17. pysad 実装まとめ

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット)。詳細は `01_overview.md` 〜 `16_fortran_status.md`、要素・機能別の `22_*.md`〜`37_*.md`、大項目の解説 `38_*.md`〜`45_*.md`、`history/work_*.md` を参照。
````

## ## 1. 概要・所在

整理前の 7 行目:

````text
- 目的: SAD (**K64 fork**, Ver.1.4.4.2k64) を Python に直訳し、SAD の数値を再現する (v1 系)。v2 系 (v2.0.0 から) は K64 のバグの修正・物理の誤りの修正・機能の追加を含み、互換モード ([`pysad/compat.py`](../pysad/compat.py)、既定 v2、`PYSAD_COMPAT=k64` で v1 と同じ動作) で切り替える。以下の「再現する」「直訳」は k64 モードの方針で、v2 の変更の一覧は `47_v2_changes.md`。SAD を物理モデルとして作り直さない (v1)。
````

整理前の 17 行目:

````text
  `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` で **1899 通過・8 skip・失敗 0** (2026-11-21、v2.3.0 の最終。1907 件収集)。以前: v2.1.0 は 1345 通過・8 skip (1353 件)、v2.0.0 は 937 通過・4 skip (941 件)。互換モードは k64 が既定 ([`pysad/tests/conftest.py`](../pysad/tests/conftest.py))。v2 の試験 (`test_v2_*.py` 73 本・1081 件。うち v2.1.0 の追加が 412 件) は `compat.set_mode('v2')` を明示する。`--ignore` する 3 本は sys.exit で終わる旧式のスクリプト。
````

## ## 2. モジュール構成と行数 (`ls pysad/*.py | wc -l` で 74 本、`cat pysad/*.py | wc -l` で計 59797 行。2026-11-18 に実測した。v2.0.0 で `compat.py`・`phsrot.py`・`radtaper.py`・`spch.py` を追加し、v2.1.0 で `naff.py`・`specfn_cx.py`・`nbcache.py` を追加した)

整理前の 20 行目:

````text
## 2. モジュール構成と行数 (`ls pysad/*.py | wc -l` で 74 本、`cat pysad/*.py | wc -l` で計 59797 行。2026-11-18 に実測した。v2.0.0 で `compat.py`・`phsrot.py`・`radtaper.py`・`spch.py` を追加し、v2.1.0 で `naff.py`・`specfn_cx.py`・`nbcache.py` を追加した)
````

整理前の 79 行目:

````text
| naff.py | 135 | NAFF 風のチューン抽出 (v2.1.0。`naff`・`naff_real`・`tune_xp`・`tunes_from_coords`)。docs/12 §16.7、47 |
````

整理前の 80 行目:

````text
| specfn_cx.py | 254 | 複素の Zeta・HurwitzZeta・PolyLog・HypergeometricU (v2.1.0)。docs/33 |
````

整理前の 81 行目:

````text
| nbcache.py | 146 | numba のキャッシュの無効化 (v2.1.0。`py -3.12 -m pysad.nbcache`)。docs/12 §16.7 |
````

整理前の 82 行目:

````text
| compat.py | 39 | 互換モード (v2.0.0 から。`mode()`・`set_mode()`・`v2()`・`k64()`)。docs/12 §16.6、47 |
````

整理前の 90 行目:

````text
テストは [`pysad/tests/`](../pysad/tests/) (`test_*.py` 158 ファイル (2026-11-21 に数えた)。`--ignore` の 3 本を除き、パラメータ化で 1907 件に展開される (収集件数)。うち v2 の試験 `test_v2_*.py` は 73 ファイル・1081 件 (v2.1.0 の追加は 412 件)。補助 `_a2_indep.py`・`_c3_indep.py`・`_e_spinkick.py` は独立な計算)。
````

## ## 6. 未対応・制限 (優先順は `history/00_RESUME.md` に従う。一覧は `46_unimplemented_and_unverified.md`、ファイル別は `16_fortran_status.md` §9。以下は主に k64 モードの状態で、v2 モードで実装したものは書き添えた)

整理前の 141 行目:

````text
- 要素: k64 モードでは SOL 領域内の BEND・MULT・CAVI・APERT (BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げた)、SOL の DX DY CHI・F1 (Twiss・Emittance[]・追跡) が未対応で、v2 モードは実装した (docs/44 §7a)。両モードで未対応: SOL 領域内の光子 (PHOTONS の tphrec。BEND・QUAD・MULT・SEXT・OCT・ドリフト・BEAMBEAM + BEND は v2.1.0 で対応)、TrackParticles の F1 (`tsteer`・`tmulti`)、
````

整理前の 143 行目:

````text
- 未解決の差 (`history/work_unimpl.md`。v2.0.0 の BEND の `F1 FRINGE` + RAD (B113)・L≠0 の CAVI の DX・DY・EMIT のあとの USE の差は v2.1.0 で解決): 
````

整理前の 148 行目:

````text
  複素 FindRoot (B61)・MessageName の関数形 (B70) は v2 モードで、PolyLog の |z|≥150 の非整数次は両モードで実装した 複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で両モードで実装した (`specfn_cx.py`)。
````

## ## 7. SADScript・FFS・パッケージの対応範囲 (docs/33〜35、history/work_unimpl.md)

整理前の 157 行目:

````text
- 試験: 1899 通過・8 skip・失敗 0 (2026-11-21、収集 1907 件。v2.3.0 の最終)。回帰 (`QUICK=1`) は 50 ベンチ FAIL 0。
````

