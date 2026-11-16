# 03 マッチング

FIT・FREE・FIX・VARY・GO・SHOW・VARIABLES・REJECT・RECOVER/SAVE/RESET、フラグ、フック関数、OptimizeOptics、オフモメンタム、特殊変数。解析解は FODO の 2x2 行列の模型 (`tools/ana_fodo_ring.py`) と薄レンズの式。

各例は 4 ファイル: `NAME.sad` (スクリプト)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明)。
走らせ方: `py -3.12 example/tools/run_example.py example/03_matching/NAME.sad`。凡例: ○ 一致 (表示の桁・書式の差は除く)、△ 差あり (内容を `.md` に書いた)、× gs 失敗、× pysad 未対応。

| ファイル | 内容 | gs と pysad |
|---|---|---|
| [FFS_FIT](FFS_FIT.md) | 1 セル FODO で位相 90 度のマッチング。厚いレンズ・薄レンズの解析解 | ○ (`VARIABLES` の Saved 列のみ違う) |
| [FFS_FREE_FIX_VARY](FFS_FREE_FIX_VARY.md) | FREE・FIX・VARY・VARIABLES | ○ (変数の並び順の表示のみ) |
| [FFS_REJECT](FFS_REJECT.md) | REJECT。条件 1 つ・変数 2 つの最小ノルム解の解析値 | ○ |
| [FFS_RECOVER_SAVE_RESET](FFS_RECOVER_SAVE_RESET.md) | RECOVER・SAVE・RESET | ○ (pysad は RESET で警告を出す) |
| [FFS_FITPOINTS](FFS_FITPOINTS.md) | FITPOINTS (オフモメンタムのマッチング、DP)。六極 K2 の 1 次の解析値 | △ (完走、値は 3e-5 以内。pysad は点ごとの表を出さない) |
| [flag_UNIPOL_BIPOL](flag_UNIPOL_BIPOL.md) | UNIPOL と BIPOL | △ (収束しない B の BIPOL で gs と pysad が別の点で止まる) |
| [flag_ABSW_RELW](flag_ABSW_RELW.md) | ABSW と RELW | △ (収束しない問題。pysad は RELW が効かない) |
| [fn_FitValue_FitWeight_VariableRange_VariableWeight](fn_FitValue_FitWeight_VariableRange_VariableWeight.md) | FitValue・FitWeight・VariableRange・VariableWeight | ○ (FitValue 時の残差の表示のみ) |
| [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (ダウンヒルシンプレックス) | ○ (5e-6 以内、反復の表示が違う) |
| [var_CONVERGENCE](var_CONVERGENCE.md) | CONVERGENCE・MAXITERATION・DP・DP0 など | ○ |

FITP>1 は docs/03 では gs が再起動するとあるが、今回の `FFS_FITPOINTS` では gs は完走した。
