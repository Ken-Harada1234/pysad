# example/03_matching の md から取り除いた履歴の記述

例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきで残す (行番号は整理前のもの)。この節は、整理前の md (リポジトリの先頭の版) と現在の md の差分から再構成した。

## README.md

元の 18 行目:
> | [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (ダウンヒルシンプレックス) | ○ (v2 で gs と同じ値・同じ反復。デバッグ表示の 2 行だけ違う) |

## fn_OptimizeOptics.md

元の 32 行目:
> **○ 一致 (2026-11-17、v2.0.0 の再実行)。** 最終値は gs と同じ: K1 (QF, QD) = (0.700221, -0.599636)、`Matched. ( 3.5504E-09)`、NX NY = 6.700009、4.900003 (6 桁)。シンプレックスの反復の表も gs と同じ行数 (最小値に達したあとの 1.000000 の行を含む) が出る。最初の版は、反復を 17 行で打ち切り、最終の `Residual = 1.9E-07`、K1 = (0.700225, -0.599635) と gs と少し違っていた (v1 の差。`BeginPackage` の文脈を扱わなかった影響と思われるが、原因はコードまでは追っていない)。互換モードによらず同じ (`PYSAD_COMPAT=k64` でも gs と同じ値を確認)。表示の差が 2 行だけ残る: デバッグ出力 `{{`,System`},System`}` の引用符 (pysad は `"System`"`)、`DownhillSimplex {limit, Option}:` の行を pysad は `...` で省略する。


# 台本 (`.sad`・`.py`) の冒頭コメントから、機能別フォルダへの移動時に取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `flag_FREE_component.sad` 1 行目:
  - 元: `! v2next B2 項目 1: 成分ごとの FREE (`FREE QF.1 K1;`・`FREE B.3 K1;`)。K64 の tffsfreefix の comp 分岐。`
  - 後: `! 成分ごとの FREE (`FREE QF.1 K1;`・`FREE B.3 K1;`)。K64 の tffsfreefix の comp 分岐。`
- `FFS_FIT_geometry_quad_dx.sad` 1 行目:
  - 元: `! v2next B3 項目 4: SOL 領域 (BOUND=1 の SOL の間) の中の四極の DX・DY を変数にした幾何の fit。`
  - 後: `! SOL 領域 (BOUND=1 の SOL の間) の中の四極の DX・DY を変数にした幾何の fit。`
- `flag_FREE_noop_variable.sad` 1 行目:
  - 元: `! v2next P4 項目 6: 要素ごとの FREE の「何もしない変数」(`FREE QF.3 K1`)。gs (K64) と全桁一致 (k64・v2 とも。K64 の動作の再現)。`
  - 後: `! 要素ごとの FREE の「何もしない変数」(`FREE QF.3 K1`)。gs (K64) と全桁一致 (k64・v2 とも。K64 の動作の再現)。`
- `flag_FITP_beta.sad` 1 行目:
  - 元: `! v2next P5 項目 6: FITP>1 (3 点 δ = -1%, 0, +1%)。NX・NY と BX・BY ($$$ の β)。12 セル FODO (3 GeV、周長 52.8 m)。自由変数は QF QD SF SD。`
  - 後: `! FITP>1 (3 点 δ = -1%, 0, +1%)。NX・NY と BX・BY ($$$ の β)。12 セル FODO (3 GeV、周長 52.8 m)。自由変数は QF QD SF SD。`
- `flag_FITP_five.sad` 1 行目:
  - 元: `! v2next P5 項目 6: FITP>1 (オフモメンタム 5 点 δ = -2%..+2%、NX NY)。12 セル FODO (3 GeV、周長 52.8 m)。gs と pysad k64 の反復の残差・最終の K を比べる。`
  - 後: `! FITP>1 (オフモメンタム 5 点 δ = -2%..+2%、NX NY)。12 セル FODO (3 GeV、周長 52.8 m)。gs と pysad k64 の反復の残差・最終の K を比べる。`
- `flag_FITP_one.sad` 1 行目:
  - 元: `! v2next P5 項目 6: FITP>1 を FITP 1 (オフモメンタムなし) と比べる基準。12 セル FODO、同じ条件で FITP 1。`
  - 後: `! FITP>1 を FITP 1 (オフモメンタムなし) と比べる基準。12 セル FODO、同じ条件で FITP 1。`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- この移動で、旧 `v2_next` の `*.pysad_k64.txt` のうち、このフォルダに移る例の分を削除した (k64 モードの出力は `PYSAD_COMPAT=k64` で再現できる。`flag_PSPAC_k64_shape` だけは k64 オプションの例として残した)。
