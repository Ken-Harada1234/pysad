# 01 ラティスと要素定義

要素 (DRIFT MARK QUAD BEND SEXT.. MULT CAVI APERT COORD)、LINE の定義と USE、FFS の要素・ライン操作コマンドと関数。ラティスは `00_lattice/FODO_ring.sad` (24 セル、3 GeV、周長 105.6 m) か、要素の効果が分かるように作った短い輸送路・小さなリング (ファイル内に定義)。

各例は 4 ファイル: `NAME.sad` (スクリプト)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明)。
走らせ方: `py -3.12 example/tools/run_example.py example/01_lattice_elements/NAME.sad`。凡例: ○ 一致 (表示の桁・書式の差は除く)、△ 差あり (内容を `.md` に書いた)、× gs 失敗、× pysad 未対応。

| ファイル | 内容 | gs と pysad |
|---|---|---|
| [elem_DRIFT_MARK](elem_DRIFT_MARK.md) | DRIFT と MARK (初期 Twiss)。ドリフトの β(s) | ○ |
| [elem_QUAD](elem_QUAD.md) | QUAD。転送行列と β の変換、薄レンズ | ○ |
| [elem_BEND](elem_BEND.md) | BEND の ANGLE・E1・E2・FRINGE。セクター型と長方形、分散 | ○ |
| [elem_SEXT_OCT_DECA_DODECA](elem_SEXT_OCT_DECA_DODECA.md) | SEXT OCT DECA DODECA。線形には効かない。六極の chromaticity 寄与 | ○ |
| [elem_MULT](elem_MULT.md) | MULT。QUAD+BEND と同じチューン | ○ |
| [elem_CAVI](elem_CAVI.md) | CAVI。Emittance[] のシンクロトロンチューンと νs の式 | ○ (`TYPE RF` の HARM/FREQ 表示のみ違う) |
| [elem_APERT](elem_APERT.md) | APERT。追跡での開口 (flag) | ○ (v2 で対応。`TYPE` の種別表示だけ違う) |
| [elem_COORD](elem_COORD.md) | COORD。軌道と幾何の変換 | ○ |
| [elem_LINE](elem_LINE.md) | LINE の定義 (n*、-、入れ子) | ○ |
| [FFS_USE_EXPAND](FFS_USE_EXPAND.md) | USE・EXPAND・NOEXPAND。成分ごとのずれが残る/消える | ○ |
| [FFS_SPLIT](FFS_SPLIT.md) | SPLIT は gs で obsolete。MARK の OFFSET で要素の中点 | △ (SPLIT はどちらも動かない。`LINE["S",小数]` が不一致) |
| [FFS_RENUMBER](FFS_RENUMBER.md) | RENUMBER (成分番号の振り直し) | ○ |
| [FFS_TYPE_ATTRIBUTE](FFS_TYPE_ATTRIBUTE.md) | TYPE と ATTRIBUTE | ○ (表示の差のみ) |
| [FFS_ORG](FFS_ORG.md) | ORG (幾何の原点)。GX の解析値 | ○ |
| [FFS_REVERSE](FFS_REVERSE.md) | REVERSE。K64 では光学の表示に効かない | ○ (両者とも変化なし) |
| [FFS_VISIT](FFS_VISIT.md) | VISIT と BYE。1 セルのチューン = リング/24 | ○ |
| [fn_Element_SetElement](fn_Element_SetElement.md) | Element と SetElement | ○ |
| [fn_LINE](fn_LINE.md) | LINE 関数 (NAME POSITION S ELEMENT MULT GEO …) | ○ |
| [FFS_wildcards](FFS_wildcards.md) | 名前のパターン (* % {} {^} |) | ○ |

リングの BEND (`E1=E2=0.1309`) は SAD では ANGLE に対する**比**と読まれるので長方形ではない。`elem_BEND` は E1=E2=0.5 で長方形にした。
