# 35_ffs_small_commands から取り除いた履歴の記述

`docs/35_ffs_small_commands.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # FFS の小さなコマンド (`pysad/ffsextra.py`)

取り除いた部分: 「最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット。」

元の行:

> 最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット。経緯は `history/` にある)

## ## gs と比べた確認 (2026-11-16 夜)

取り除いた部分: 「(成分ごとの FREE) があれば真 (v2.1.0 から。docs/03 §3.3)。」

元の行:

> - `RENUMBER comp` と `LINE["MULT"]` (tfinimult の直訳。同名 1 台だけは 0)、`USE [NO]EXPAND`・`CALC|GO [NO]EXPAND` は gs と一致 (docs/03)。`NOKEEPEXP` (tffsa.f の `calexp = .not. tfvcomp()`) と RENUMBER の INDEX 属性 (2 以上なら同名の要素の番号をその値から振る。`LINE["INDEX",…]` も読める) は実装した (両モード。gs の実測と一致)。`tfvcomp` は、既定でないキーの成分の変数 (成分ごとの FREE) があれば真 (v2.1.0 から。docs/03 §3.3)。
