# SADScript・入出力を gs の実動作に合わせた作業 (2026-11-17)

担当: pysad/sadscript.py、pysad/ffscmd.py (Random・SeedRandom、`echo_out`)、tests/test_scriptgs.py、docs/04_sadscript.md、example/08_scripting・09_io_system

起点は `example` 08・09 の「gs と pysad の比較」の差。基準は gs (K64 2026-10-19 ビルド) の実動作。各項目を `local/wake/` の短い台本で gs と pysad に同じものを走らせて比べた。

## 直したこと
- (a) **Module の字句スコープ**: 宣言した局所名を呼び出しごとに `名前$N` に改名した式を評価する (`_module_lex`)。別関数から局所変数が見えない、初期値は外側の環境、`Module[{t},t]` は `t$N`、局所の関数定義 `f[x_]:=…` が大域の f と混ざらない。定義・純関数・記号を返す Module は出たあとも値を大域に残す。`Block` は従来の動的スコープ。引数が記号に束縛された名前の宣言 (`_resolve` が別名を返す場合) は従来経路。gs は `t$1073741824` のような番号を出すことがあるが、pysad は 1 からの連番 (表示だけの差)。
- (b) **Return**: gs は純関数・`Do`・`Scan`・`Map`・`Module`・`Block` を抜けて、外側の `:=` のユーザー定義関数まで届く。pysad は純関数が捕まえていたので外した。関数の外の Return/Break/Continue/Throw は `General::unexpbreak` で文を打ち切る。
- (c) **Switch/Which の一致なし**: 式を返す。`Which` は引数そのまま、`Switch[x,…]` は第 1 引数だけ評価した `Switch[5,1,2]`。引数がすべて定数なら gs が読み込み時に畳むので `Switch[9]`。引数の個数の誤りは `General::narg` (Switch は奇数、Which は偶数)。
- (d) **リストの比較**: `==`・`<>`・`>`・`>=`・`<`・`<=` はリストの要素ごと (入れ子も。`{1,2}=={1,3}` は {1,0}、`{1,2}>=2` は {0,1})。長さの違うリストの `==` は 0。リストとスカラーの `==` は未評価。`Thread[{1,2}=={1,3}]` は引数が先にリストになるので通る。
- (e) **Out[n] と %**: run_sad_string で `;` なしの最上位の式文 (値が Null でないもの) を `Out[n]:= 値` と表示し、`Out[n]`・`$Line` に記録。`%` は `Out[$Line]`、`%n` は `Out[n]`、`%%` は記号のまま、未定義の `Out[n]` は未評価。`f[x_]:=body` の値は未評価の本体 (gs)。
- (f) `Close[{f,g}]`、(g) `Write`・`WriteString` の `$Output` (-1) は標準出力、`Write[$Output]` は narg。`Close` は未使用の番号・`$Output` を黙って無視 (gs は `Close[$Output]` で Fortran のエラーで落ちる)。
- (h) `DateString[t]`: 1900-01-01 からの秒、または {Y,m,d[,H[,M[,S[,offset]]]]} (offset は無視、未確認)。
- (i) **Read**: tfwrite.f の tfreadstringf・tfword・tfreadfs を直訳した読み取り器 (記録 = 1 行 + 位置)。Word の区切り (既定は空白・`,`・タブ・改行)、リスト、`n*item` (頭 Null のリスト)、Expression (隣り合う項で切る)、Character、String、ReadString、WordSeparators・ReadNewRecord・NullWords。`fn_Read_ReadString`・`fn_StringToStream` が gs と一致。
- (j) 文字列の `\` エスケープ: tfetok.f kxmakestring の直訳 (8 進 1〜3 桁、`\x`、`\u`、`\e`、`\f`、未知の `\q` は q)。
- (k) `StringMatchQ`: `*` `%` `{…}` `{^…}` `|`、大文字小文字を区別、`?` は普通の文字、空文字列は 0、非文字列は wrongtype。
- (l) `Random[n,m,…]`・`GaussRandom[n,m]`・`ParabolaRandom[n,m]` は多次元配列 (行優先に順に引く)。負は wrongnum、リストは wrongtype。
- (m) `SeedRandom`: 種類の選択 (`"SAD"`、先頭が文字列のリスト) は状態を変えず現在の状態を返す。"SAD" 以外は `System::error`。状態の戻りの next は IEEE 倍精度の下位・上位の 2 語。数だけのリストは種 1 (tfRandom_.c が `vi` を読まない癖: gs では 0 になる。未初期化変数なので確実ではない)。**`SeedRandom[状態]` で状態は復元されない** のが gs の動作で、pysad も同じ。
- 表示: 比較式は `a>b` の中置、`Blank` は `_String`。

## 直さなかったこと
- **For の本体**: gs は実行しない。ソース (`Packages/functions.n:34` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)`) は本体を実行する形で、`src/` に For の組み込みはなく (`While` は tfmap.f の `tfwhile`)、同じ式を `For2` として定義すると gs も本体を実行する。`Clear[For]` して再定義しても実行する。名前解決の癖と思われるが原因は未特定。SAD の意図は「本体を実行する」と判断し、pysad は変えない (docs/21 B71)。3 引数の For は gs と同じく未評価で返すようにした。
- **BeginPackage/Begin の文脈**: gs の動作を調べた (読み込み時に新しい記号が `Ctx`名` になる。`Context[]`・`$Context`・`$ContextPath`・`Names["Stat`*"]` の値が独特で、EndPackage のあとも文脈なしの名前は見えない)。読み込み時の名前解決の模倣が要り、SAD の Packages の読み込みに影響するので未対応とした (docs/04 §6)。
- `Check[Return[1],2]` (gs は 2 を返す)、`Thread[a==b]` (a に値、b が未定義: gs は `3==b`)、記号どうしの `==` の未評価形、未定義の記号だけの文の Out の繰り返し、`"x\\y"` の `\\` に見える 2 重解釈 (結果は pysad も gs と同じだったので触っていない)、`DateString` の offset。
- `System`・`Fork`・`Shared` (09 の差) は今回の課題の外。

## 試験
`pysad/tests/test_scriptgs.py` (16 件、値は gs から)。全試験と回帰の結果は 00_RESUME.md §4。
試験結果: 548 通過・3 skip、回帰 QUICK=1 で 50 ベンチ中 FAIL 0 (125 s)。
