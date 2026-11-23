# example/08_scripting の md から取り除いた履歴の記述
例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきでそのまま残す (行番号は整理前のもの)。取り除いたあとの現在の記述は、各 md を参照。
整理日: 2026-11-23。

## README.md

元の 11 行目:
> | fn_For | fn_For (For、本体を実行しない gs の癖) | △ (v2 は K64 のバグ B71 を直したので gs と違う。k64 モードは一致) |

元の 39 行目:
> 見つけた差のうち、`Module` の字句スコープ、`Return` の範囲、`Switch`/`Which` の一致なし、`{1,2}=={1,3}`・`Thread`、`Out[n]`・`%`、`StringMatchQ`・8 進エスケープ、`Random[n,m]`・`SeedRandom` は 2026-11-17 に gs に合わせた。残る差は `For` (gs は本体を実行しない K64 の癖。pysad v2 は定義どおり実行、docs/21 B71。k64 モードは gs と同じ) と `BeginPackage`/`Begin` の `Names["Stat`*"]` (v2 は文脈の記号の一覧、gs は `{"Stat`"}`)。文脈そのものは v2 で対応した (2026-11-17 に再実行して gs と一致を確認)。乱数の列は一致した。

## fn_Apply_Thread.md

元の 17 行目:
> ○ 一致する (2026-11-17 に修正: リストどうしの `==` を要素ごとにしたので、`Thread[{1,2}=={1,3}]` は gs と同じ {1,0})。

## fn_For.md

元の 18 行目:
> △ **v2 では K64 のバグを直したので gs と違う (docs/47 の B71)。** `PYSAD_COMPAT=k64` で走らせると gs と同じ (`t=0 i=5`、`body` の `Print` も出ない。2026-11-17 に確認)。v2 の既定 (`.pysad.txt` はこちら): pysad は `body` を実行して t=30、`body i=0,1,2` も表示する。gs (K64、2026-10-19 のビルド) は `body` を実行せず、3 番目の `incr` だけを繰り返す。`Packages/functions.n:34` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)` は本体を実行する形で、`src/` に For の組み込みはなく、同じ式を別名で定義すると gs も本体を実行する (`Clear[For]` で再定義しても実行する)。したがって gs の動作は定義とも SAD の意図とも合わない K64 の癖 (読み込み時の名前解決と思われる。原因は未特定。docs/21 B71) と判断し、pysad は定義どおり本体を実行する。3 引数の `For` は gs と同じく未評価で返す。`For` は避けて `Do`・`While` を使うのが安全。

## fn_Function_Options_Return.md

元の 20 行目:
> △ (1) `Return` の範囲は 2026-11-17 に gs に合わせた (純関数・`Scan`・`Do`・`Module` を抜けて外側のユーザー定義関数まで届く。`f3[{1,-1}]` は "negative found")。(2) `Definition[fact]` は gs が `Hold[{fact:=fact, ...}]` を返し、pysad は `Null` (未対応のまま)。それ以外 (オプション、既定値、再帰、HoldAll、個数別定義) は一致。

## fn_Get_Needs_BeginPackage.md

元の 20 行目:
> △ **ほぼ一致 (2026-11-17、v2.0.0 の再実行)。** pysad は `BeginPackage`・`Begin` で文脈を作るようになった (最初の版では文脈を無視して定義が大域に入り、`Stat`Mean2[...]` が評価されなかった。v1 の未対応)。`Stat`Mean2[{1,2,3,4}]` = 2.5、`Stat`Var2` = 1.25、`Stat`Calls[]` = 2、`Stat`Counter` = 2、`Context[Stat`sq]` = `Stat``、文脈なしの `Mean2[{1,2}]` が未評価、`Context[]` が読み込み前後とも `System`` は gs と同じ。文脈の規則は `docs/47_v2_changes.md` の「文脈」。

## fn_Module_Block_With.md

元の 20 行目:
> ○ 一致する (2026-11-17 に修正: `Module` を字句スコープにした。局所変数を `名前$N` に改名して評価するので、`show[]` を `Module[{x = 7}, ...]` の中で呼んでも大域の 10 を返し、`Module[{t}, t] === t` は 0)。局所記号の表示は gs が `t$1073741824` のような番号、pysad は連番 (表示だけの差)。

## fn_ParabolaRandom_ListRandom.md

元の 20 行目:
> ○ ヒストグラムは **完全に一致**。`SeedRandom["SAD"]` (プラグイン名の文字列) も 2026-11-17 の修正で gs と同じ (乱数の種類の選択だけで状態は変わらない)。`ListRandom[]` も一致。

## fn_Random_GaussRandom_SeedRandom.md

元の 21 行目:
> ○ **乱数の列は gs と pysad で一致** (平均 0.503181122013554、`GaussRandom` の最大値 4.03258726603887 まで同じ)。2026-11-17 に 2 点を gs に合わせた: (1) `Random[2,3]` は 2×3 のテンソル ({2,3})。(2) `SeedRandom[状態]` は gs でも復元されない (gs の C 実装はリストの数を読まず、先頭が文字列なら乱数の種類の選択だけで状態は変わらない) ので、pysad も同じ (`a===b` は 0)。`SeedRandom["SAD"]` は状態を変えず現在の状態を返す。

## fn_String_functions.md

元の 19 行目:
> ○ 一致する (2026-11-17 に修正)。`StringMatchQ` は SADHelp (wildcards) の `*`、`%`、`{..}`、`{^..}`、`|` を gs と同じに扱い (`?` は普通の文字、大文字小文字を区別、空文字列は 0)、8 進エスケープ `"\101"` は "A" (65)。

## fn_Switch.md

元の 18 行目:
> ○ 一致する (2026-11-17 に修正: 合うものがないとき `Switch[9]` を返す。引数がすべて定数なら gs が読み込み時に畳んで `Switch[9]` になり、変数を含むと `Switch[5,1,2]` のように第 1 引数だけ評価した式を返す)。

## fn_Which.md

元の 17 行目:
> ○ 一致する (2026-11-17 に修正: どれも真でないときは式をそのまま返す)。引数の個数が合わなければ `General::narg`。

## syntax_Operators.md

元の 20 行目:
> ○ 一致する (2026-11-17 に修正: リストどうしの比較は要素ごと、`{1,2} == {1,3}` は {1,0})。

## syntax_Semicolon_Toplevel.md

元の 22 行目:
> ○ 一致する (2026-11-17 に修正: `;` のない最上位の式の値を `Out[n]:= 値` で表示する)。値の末尾の桁は表示の差 (NX=6.696328940901385 gs、...378 pysad)。pysad の CALCULATE は `Matched.` の行を表示する (gs も表示する。`.gs.txt` で確認)。

## var_Null_Out.md

元の 19 行目:
> ○ `%`・`%n`・`Out[n]`・`$Line` は一致 (2026-11-17 に修正: `%` は `Out[$Line]`、`%%` は gs でも記号のまま、未定義の `Out[-1]` は未評価)。差が残るのは `BoundQ` (gs は `BoundQ[zzz]` を 0、`BoundQ[x]` を 1 と評価、pysad は未評価): 今回の対象外。

## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### syntax_context.sad (旧 v2_next/B2_context.sad)
> (* v2next B2 項目 4: $Context・$ContextPath への代入、WithContext、入れ子の文脈 …
### syntax_symbol_message.sad (旧 v2_next/B2_symbol_message.sad)
> (* v2next B2 項目 6: Symbol["x"] は x の値 …
### FFS_VISIT_BYE.sad (旧 v2_next/B2_visit_bye.sad)
> ! v2next B2 項目 5: VISIT で FitFunction・… を退避し、BYE で戻す …
### syntax_package_context.sad (旧 v2_B/B2_package_context.sad)
書き出すパッケージのファイル名 `v2_B_stat_pkg.n` を `stat_pkg.n` にした。
