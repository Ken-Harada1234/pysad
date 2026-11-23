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
| [FFS_FITPOINTS](FFS_FITPOINTS.md) | FITPOINTS (オフモメンタムのマッチング、DP)。6 極 K2 の 1 次の解析値 | △ (完走、値は 5e-5 以内。pysad は点ごとの表を出さない) |
| [flag_UNIPOL_BIPOL](flag_UNIPOL_BIPOL.md) | UNIPOL と BIPOL | △ (収束しない B の BIPOL で gs と pysad が別の点で止まる) |
| [flag_ABSW_RELW](flag_ABSW_RELW.md) | ABSW と RELW | △ (収束しない問題。pysad は RELW が効かない) |
| [fn_FitValue_FitWeight_VariableRange_VariableWeight](fn_FitValue_FitWeight_VariableRange_VariableWeight.md) | FitValue・FitWeight・VariableRange・VariableWeight | ○ (FitValue 時の残差の表示のみ) |
| [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (ダウンヒルシンプレックス) | ○ (gs と同じ値・同じ反復。デバッグ表示の 2 行だけ違う) |
| [var_CONVERGENCE](var_CONVERGENCE.md) | CONVERGENCE・MAXITERATION・DP・DP0 など | ○ |
| [flag_FREE_component](flag_FREE_component.sad) | 成分ごとの FREE (`FREE QF.1 K1;`・`FREE B.3 K1;`) | ○ (gs と一致) |
| [flag_FREE_noop_variable](flag_FREE_noop_variable.sad) | 要素ごとの FREE の「何もしない変数」(`FREE QF.3 K1`)。VAR の表・エラー・フィットの変数の数 | ○ (全桁一致。v2・k64 とも) |
| [FFS_FIT_geometry_quad_dx](FFS_FIT_geometry_quad_dx.sad) | SOL 領域内の四極の DX・DY を変数にした幾何の fit | ○ (gs と v2 が 2 回の Newton で収束し、QF DX = 2.788484326e-3、QD DY = 2.667980563e-3 が 1e-12 で一致。k64 は SOL 自体のずれが未対応で動かない) |
| [flag_FITP_five](flag_FITP_five.sad) | FITP>1 (オフモメンタム 5 点 δ = −2%〜+2%、NX NY)。12 セル FODO (3 GeV、52.8 m) | △ (最終の残差が gs と 2e-4 以内。反復の経路は違う) |
| [flag_FITP_beta](flag_FITP_beta.sad) | FITP>1 (3 点 δ = −1%・0・+1%)、NX・NY と BX・BY (`$$$` の β)。自由変数は QF QD SF SD | 同上 |
| [flag_FITP_one](flag_FITP_one.sad) | FITP>1 と比べる基準 (FITP 1、オフモメンタムなし) | 同上 |

上の表のうち、拡張子つきのリンクの行 (`.sad`・`.py`) は説明の `.md` を持たない。内容と結果の要点は表と、台本・スクリプトの先頭のコメントに書いた。`.py` の例は、リポジトリ直下から `py -3.12 example/03_matching/NAME.py` で走らせる (出力は `NAME.out.txt`)。gs と同じ動作で走らせるときは `PYSAD_COMPAT=k64`。

FITP>1 は docs/03 では gs が再起動するとあるが、今回の `FFS_FITPOINTS` では gs は完走した。
