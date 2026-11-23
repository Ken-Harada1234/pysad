# 03_ffs_commands の照合記録

旧 03_ffs_commands.md (同名で上書き)。§1〜§9 の構成と表はそのまま。

| 変更 | 内容 |
|---|---|
| 冒頭 | この文書で分かること・前提・読む順・用語を追加 (「最終更新」は削除) |
| §2 BIPOL・ABSW の行、§3.3 成分ごとの FREE、§5 EMIT のあとの USE、§7 SynchroBetaEmittance・TrackParticles の l_track | gs との数値 (全桁一致の値、反復の残差など) を history/03_ffs_commands_validation.md に移し、本文は 1 文に短縮 |
| §8 末尾の「Twiss["DETR"] が CALC の後で常に 0 (N12 未修正)」 | 誤りだったので訂正: INS の CALC のあとも値が入る (pysad で確認、入口の位置 1 は R が全て 0 で DETR=0)。DETR ≥ 1 の INS・MEASURE は未確認のまま |
| §10 検証 | 追加 (K64 との比較の要約 5 行) |
| 参照 | docs/31・35→04、23→13、25・26→19、29→17、45→18、36→16、47→28、11→29、34→06、33・32・34→15・06、「08 の 1.2 節」→docs/10 |

落とした表・数値・注意は無い (移した記述は validation に原文のまま)。§8・§9 に残る gs との照合の記述 (DISP 36 ブロック一致など) は短いので本文に残した。
