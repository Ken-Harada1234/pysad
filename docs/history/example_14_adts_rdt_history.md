# example/14_nonlinear_dynamics/adts_rdt の文書から取り除いた履歴の記述

元は `example/v2_3/` の作業フォルダ (旧パス → 新パスは `example/14_nonlinear_dynamics/adts_rdt/` の各ファイル名の対応で、台本・出力の接頭辞 Q2_・R1_〜R6_ を機能名に替えた)。取り除いた文章を、元の文脈つきでそのまま残す。

## rdt_adts_basic.md (旧 `Q2_README.md`)
- 冒頭: 「(`docs/49_resonance_terms_adts.md`、記録は `docs/history/v23_Q2.md`)。(`example/v2_3/README.md` が索引。)」
- 結果表の gs の列の注: 「ASTRID2: … (NORFSW。RF ON の旧出力は 452.7, 1318.5, 1168.1)」「大きなリング A: … (NORFSW。RF ON の旧出力は 121125, −46333, 22196)」。
- 末尾: 「`NORFSW;` を入れると gs は pysad (k64) と解析に一致する。`R1_README.md`、`docs/history/v23_R1.md`。」

## adts_gs_norfsw.md (旧 `R1_README.md`)
- 冒頭: 「gs の追跡の ADTS が pysad と食い違うときの原因と切り分けは `docs/history/v23_R1.md`。」
- 表: 「`R1_adts_gs_*_norfsw.txt`: gs `NORFSW;`、1024 周 (Q2 の C・B が解析と一致する)」「保存してある `adts_gs_*.txt` はこの条件」。

## adts_direct.md (旧 `R2_README.md`)
- 「記録は `docs/history/v23_R2.md`」「詳細は `docs/history/v23_R2.md` の表」。

## adts_bend_cubic.md (旧 `R3_README.md`)
- 「記録は `docs/history/v23_R3.md`」。

## rdt_tracking.md (旧 `R4_README.md`)
- 「記録は `docs/history/v23_R4.md`」。

## rdt4_normal_form.md (旧 `R5_README.md`)
- 「記録は `docs/history/v23_R5.md`」「`FMAP_summary.md` §12」「(原因は) `FMAP_summary.md` §12.6」。

## adts_exact_terms.md (旧 `R6_README.md`)
- 「記録は `docs/history/v23_R6.md`」。表の列名「精密版の項なし」、(a) の原因欄「K0≠0 の偏向磁石 (内部の軌道が ±0.032 rad 傾く) を入れていなかった。」(現在の言い方: 「K0≠0 の偏向磁石の内部の軌道」)。

## ADTS_summary.md (旧 `ADTS_summary.md`)
- 冒頭: 「記録は `docs/history/v23_Q2.md`・`v23_R1.md`・`v23_R2.md`・`v23_R3.md`」。
- §2.6: 「(J = 2e-9 で 5e-4、`docs/history/v23_Q2.md` の記録)」。
- §6.5: 「`docs/history/v23_R2.md` の記録の値は「比 −1.000000 ± 6e-7、位相も一致」」、§6.6 の表: 「`v23_R2.md` の記録では 6e-7」。
- §7.2: 「記録 `docs/history/v23_R1.md`」、§7.5: 「切り分けの実験 (`v23_R1.md`)」。
- §8.6: 「(`docs/history/v23_R2.md` の記録は (追跡 − 解析)/解析 で書いてあるため …)」「(`v23_R2.md` の記録の値)」、§9.5: 「(`docs/history/v23_R2.md` の記録の「解析/追跡の比」は …)」。
- §12・§13 の冒頭: 「記録は `docs/history/v23_R6.md`」。

## 台本の冒頭コメント (docstring)
作業符号と版の接頭辞 (「v2.3.0 Q2:」「v2.3.0 R2:」〜「v2.3.0 R6:」) を取り除いた。コメントの「(R4 以後の既定は偏向磁石の項つき)」「(v2.3.0 R4 以後、resonance_terms の既定は偏向磁石の 3 次の項つき)」「R2 の値 (6 極のみ。R4 は rdt_tracking.py res)」「Q2 の図は 6 極だけ」「(R1 の対照実験用)」「docs/history/v23_R1.md」を、現在の言い方に直した。
`adts_exact.py` の冒頭: 「「結果が合っていない・原因が未解明」の解消 (a)〜(e) の表を作る。」
出力ファイルの中の表の見出し・JSON のキー (`R2`・`R3`・`rel_R2` など。偏向磁石の 3 次なし = R2、あり = R3 の意味) は、数値ファイルなので書き換えていない。
