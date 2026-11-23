# 17_sol の統合の照合

docs/17_sol.md を、旧 docs/44_sol_guide.md (土台) と旧 docs/29_sol_optics.md から統合したときの節の対応。

## 旧 44_sol_guide.md

| 旧の節 | 新 17 での場所 |
|---|---|
| §1 何ができるか (表) | 1 |
| §2 SOL 領域の仕組み・最小例・追跡と放射の例 | 2 (最小例の結果の意味を追加) |
| §3 領域内に置ける要素 | 3 |
| §4 各機能の対応 (追跡・Twiss・Emittance・CALC6D・幾何) | 4.1〜4.5 (CALC6D の確認の数値は validation へ) |
| §5 K64 の癖 | 5.2 |
| §6 gs との一致 (表) | history/17_sol_validation.md。本文は 6 の要約 |
| §7 限界・未確認 (1〜10) | 7 (項目 9 は確認して書き直した。下記) |
| §7a v2 モードの SOL 領域 | 5.3 (5.3.1〜5.3.6) |
| §8 関連文書 | 9 |

## 旧 29_sol_optics.md

| 旧の節・表 | 新 17 での場所 |
|---|---|
| 冒頭 (領域の定義、モードによる違い) | 2・5 |
| 作り方 | 8 (Twiss の実装) |
| Emittance[]・EMIT の K64 と solemit の対応表 | 8 |
| 検証 (Twiss・Emittance[]) | history/17_sol_validation.md。本文は 6 の要約 |
| 限界・未確認 | 7 (44 §7 と重複する項目は 1 つに。29 にだけある GEO の周長の説明は 7 の 11、Twiss の座標は 12) |
| 幾何 (tsgeo) | 4.5 |
| v2 モードで対応する範囲 (表) | 5.1 |

## 食い違いの確認 (旧 44 項目 9 と §1)

- 旧 44 §7 の項目 9: 「k64 モードでは領域内の MULT・BEND・CAVI の tmulte1 と、tsgeo の領域幾何、CALC6D、Touschek は未翻訳」。
- 旧 44 §1 の表: 「CALC6D、RADCOD つきの CALC は対応 (両モード)」「Touschek (INTRA) との組み合わせは対応 (両モード)」。旧 44 §4 の幾何の節も「対応 (sgeo.py)」。
- 確認の方法: (a) ソース。`pysad/sgeo.py` に `tsgeo` の領域幾何、`calc6d.py` (339 行付近)・`touschek.py` (503 行付近) に `_solmap` の領域内の要素を通す分岐、`ffscmd._check_sol_optics`。(b) 実行 (2026-11-23、`example/v2_next/B1_sol_calc6d_intra.sad` を `PYSAD_COMPAT=k64` と既定 (v2) で実行): k64 は CALC 4D・CALC6D・RADCOD つき CALC6D・INTRA の Emittances が完走し、gs の値と 1e-9 で一致 (C4D 2.73475252978794、C6D 2.73495437251039、INTRA εx 2.706114557087e-7)。v2 も完走 (C4D 2.73475252592792、INTRA εx 2.787024683991e-7)。(c) 領域内に BEND を置いた台本は、k64 で NotImplementedError、v2 で完走 (チューン 0.29744、0.10784)。
- 結論: 項目 9 は誤り (`tsgeo` の領域幾何は DRIFT・QUAD・SOL について両モードで翻訳済み、CALC6D・Touschek は両モード対応)。k64 モードで未翻訳なのは領域内の MULT・BEND・CAVI (`tmulte1`) だけ。`tdgeo` の領域内の変数は、k64 では K64 と同じ領域の外の式 (`match6d.tdgeo_term`)、v2 では中心差分 (match6d.py のコメントによる。実行しての確認はしていない)。この結果で 7 の項目 9 と 4.5・5.1 の表を書き直した。

## 意図して削った・移したもの

- 落としたものは無い。gs との一致の表、CALC6D の確認の数値、旧 29 の検証の記述は validation に移した。旧 44 §7a の「記録」(旧 history へのリンク) は 5.3 の本文に残した。
