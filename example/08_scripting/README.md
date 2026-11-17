# 08_scripting 一覧
SADScript の関数・文法・クラス化・乱数の例。ラティスは使わない (`syntax_Semicolon_Toplevel` だけ `FODO_ring` を使う)。比較は、gs と pysad が同じスクリプトで出す画面表示 (`Print` など) どうしで行う。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| fn_Apply_Thread | fn_Apply_Thread (Apply, @@, Thread) | ○ |
| fn_Cases_DeleteCases_Position | fn_Cases_DeleteCases_Position (Cases, DeleteCases, Position) | ○ |
| fn_Class | fn_Class (Class[], メンバー・メソッド・継承・コンストラクタ) | ○ |
| fn_Do | fn_Do (Do, Break, Continue, Catch, Throw) | ○ |
| fn_FixedPoint | fn_FixedPoint (FixedPoint, FixedPointList, Nest, NestList) | ○ |
| fn_For | fn_For (For、本体を実行しない gs の癖) | △ (v2 は K64 のバグ B71 を直したので gs と違う。k64 モードは一致) |
| fn_Function_Options_Return | fn_Function_Options_Return (オプション・既定値・Return・再帰) | ○ |
| fn_Function_Pattern | fn_Function_Pattern (defining-functions, pattern) | ○ |
| fn_Get_Needs_BeginPackage | fn_Get_Needs_BeginPackage (Get, BeginPackage, Begin, End, EndPackage, 文脈) | △ (`Names["Stat`*"]` だけ違う。k64 モードは一致) |
| fn_If | fn_If (If, 論理演算, FORTRAN 風 IF-ELSEIF-ELSE-ENDIF) | ○ |
| fn_List_functions | fn_List_functions (Take, Drop, Join, Flatten, Sort, ...) | ○ |
| fn_Map_MapThread | fn_Map_MapThread (Map, /@, MapAll, MapThread) | ○ |
| fn_Module_Block_With | fn_Module_Block_With (局所変数のスコープ) | ○ |
| fn_ParabolaRandom_ListRandom | fn_ParabolaRandom_ListRandom (ParabolaRandom, ListRandom, SeedRandom[plugin]) | ○ |
| fn_Random_GaussRandom_SeedRandom | fn_Random_GaussRandom_SeedRandom (Random, GaussRandom, SeedRandom) | ○ |
| fn_Scan_ScanThread | fn_Scan_ScanThread (Scan, ScanThread) | ○ |
| fn_String_functions | fn_String_functions (character-string の項) | ○ |
| fn_Switch | fn_Switch (Switch) | ○ |
| fn_SwitchCases | fn_SwitchCases (SwitchCases, SelectCases) | ○ |
| fn_Table | fn_Table (Table, Sum, Product) | ○ |
| fn_Which | fn_Which (Which) | ○ |
| fn_While | fn_While (While, Break, Continue) | ○ |
| syntax_Brackets | syntax_Brackets (() [] {} [[ ]]) | ○ |
| syntax_Comments | syntax_Comments (! と (* *)) | ○ |
| syntax_Iterator | syntax_Iterator (反復子 {n} {i,n} {i,a,b,step}) | ○ |
| syntax_Operators | syntax_Operators (= := == === など) | ○ |
| syntax_Rules | syntax_Rules (/. -> :> //. Replace) | ○ |
| syntax_Semicolon_Compound | syntax_Semicolon_Compound (括弧・関数の中の ;) | ○ |
| syntax_Semicolon_Toplevel | syntax_Semicolon_Toplevel (最上位の ;) | ○ |
| var_Null_Out | var_Null_Out (Null, Out[n], %, `$Line`) | ○ |

○一致 / △差あり (原因は各 .md) / ×gs 失敗 / ×pysad 未対応。

見つけた差のうち、`Module` の字句スコープ、`Return` の範囲、`Switch`/`Which` の一致なし、`{1,2}=={1,3}`・`Thread`、`Out[n]`・`%`、`StringMatchQ`・8 進エスケープ、`Random[n,m]`・`SeedRandom` は 2026-11-17 に gs に合わせた。残る差は `For` (gs は本体を実行しない K64 の癖。pysad v2 は定義どおり実行、docs/21 B71。k64 モードは gs と同じ) と `BeginPackage`/`Begin` の `Names["Stat`*"]` (v2 は文脈の記号の一覧、gs は `{"Stat`"}`)。文脈そのものは v2 で対応した (2026-11-17 に再実行して gs と一致を確認)。乱数の列は一致した。
