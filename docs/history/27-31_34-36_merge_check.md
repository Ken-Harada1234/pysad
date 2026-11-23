# 27〜31・34〜36 の照合表 (旧文書 → 新文書)

| 旧文書 | 新文書 | 内容 | 対応・変更 |
|---|---|---|---|
| 10_benchmarks.md (全節) | 27_benchmarks.md | 全節 (§1 方針〜§6) を保持。冒頭に「この表は何を意味するか」(§0) を追加 | §5.1 の件数を最新値 (全試験 2507 通過・13 skip・失敗 0、回帰 QUICK=1 の 50 ベンチ FAIL 0) に更新。旧値 (937 通過・4 skip) は削除。§6 の文書参照を新番号へ。gs との比較が主題の文書なので validation への移動はしていない。「旧版の実行ファイル」などの語りを現在形に |
| 47_v2_vs_k64.md (全節) | 28_v2_vs_k64.md | §1〜§6 の全行を保持 | 「要確認」を確認 (B144 は未調査のまま、B185〜B196 は目録で点ごとの値を照合して記入、B 番号なしの行は由来を記入、WSPAC の行は B149 に対応) |
| 11_sad_quirks.md | 29_sad_quirks.md | 全節 | 番号・リンク・「実装した」の現在形化のみ |
| 21_sad_bug_report.md | 30_sad_bug_report.md | 全節 (内容不変) | 見出し番号と `47_v2_changes.md` へのリンク 2 か所を 28 へ |
| 46_unimplemented_and_unverified.md §1〜§5、§7〜§9 | 31_unimplemented_and_unverified.md §1〜§5、§7〜§9 | 全行を保持 (§9 の詳細も) | 変化の語りを現在形に。`ElementValues` の行は同じ行に新旧の記述が混在していたため 1 つの現在の記述に整理 (事実は変えず) |
| 46 §6 | 31 §6 | 現在形に直した表 | 元の全文は history/31_unimplemented_and_unverified_history.md |
| 16 §9、17 §6〜§8、19 §9〜§11、20 §5〜§6 | 31 | 1 項目ずつ 31 と照合 (約 126 項目、31 に既にあるもの約 85、旧記述が古く 31 が正しいもの約 26) | 31 に無かった項目を §1 (tspac・FFSFlags・EPICS・死んだコード)、§2 (MULT の TAPER・空間電荷の範囲・CALC 解析行列/tdgeo・BIPOL)、§4 (1e199 級の値・旧 API) に追加。旧 16 §10 は docs/33 §10 に残る。旧 16 §10 の「trackb の再数えは未実施」は docs/33 側 |
| 13_sad_help_index.md | 34_sad_help_index.md | 全節 | 冒頭の説明を追加 |
| 14_pysad_vs_sadhelp.md | 35_pysad_vs_sadhelp.md | 全節 | 参照を新番号へ、冒頭の説明、「実装した」を現在形に |
| 15_sadhelp_purpose_policy.md | 36_sadhelp_purpose_policy.md | 全節 | 同上 |
