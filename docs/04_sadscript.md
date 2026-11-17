# 04. SADScript 言語と pysad の sadscript.py

最終更新: 2026-11-17 (現状のスナップショット)

対象: `pysad/sadscript.py` (字句解析・Pratt 構文解析・評価器、約 4200 行)、`pysad/sadlist.py` (リスト・文字列・行列の関数)、
`pysad/sadelem.py` (FFS 専用関数)、`pysad/teigen.py` (Eigensystem)、`pysad/sadclass.py` (Class)、`pysad/pkgload.py` (Packages/*.n の読み込み)、`pysad/specfn.py`・`fitroot.py`・`miscfn.py`・`fileops.py` (特殊関数・Fit・ファイル関数など)。SAD 側の根拠は
`src/tfetok.f` (字句)、`src/tfstk.f` (演算子表 `opcode`/`iprior`)、`src/autofg.f` (数値書式)、
`Packages/Fit.n` (PolynomialFit)、`src/tffsa.f` (FFS コマンドの振り分け)、`src/teigen.f`、`src/tmatch.f`。

## 1. SADScript の概要と SAD 内での位置付け

SADScript は SAD に内蔵された Mathematica 風の式言語である。`f[x, y]` の関数呼び出し、`{...}` のリスト、
`x_` のパターン、`->` の規則、`#&` の純関数、`/@` の Map など、表記はほぼ Mathematica と同じだが、
数値は倍精度実数のみ (整数型・有理数・任意精度はない) で、`True`/`False` も 1/0 の実数である。

SAD のスクリプトは MAIN 部 (ラティス定義) と `FFS;` 以降の FFS 部に分かれる。FFS 部では各文を
`tffsa.f` が先頭の語で判定し、`CALC`, `FIT`, `FREE`, `GO`, `USE`, `CELL`, `EMIT` などの FFS コマンド
(`abbrev(word,'DELW_AVE ','_')` のような省略形付きの比較) と照合する。SAD では先頭の `tfprint` (tffsa.f:245) が
`itfpeeko` で先読みし、定義済みシンボルや `name[`, `name=` のように式に見える語はコマンド照合より先に SADScript として
評価する (03_ffs_commands §1)。pysad は判定の向きが逆で (コマンドでなければ式)、`Interp.run(src, 'ffs')` が文ごとに
`_is_ffs_command` (先頭語が `ffs_words` に登録されたコマンドか、かつ続きが `=`, `[`, `(` などで
始まらない) を見て、コマンドなら `ffs_handler` (ffscmd.py の `FFS.execute`) へ文字列のまま渡し、
それ以外は `parse` → `eval` する。

モジュールは `Get["file"]` で読み込む。SAD は起動時に `Packages/init.n` を読み、`functions.n`,
`MessageName.n` を `Get` し、`Fit.n` (`PolynomialFit`, `FindRoot`) や `Matrix.n`, `Optics.n` などは
`AutoLoad[symbol, file]` で最初に使われたときに読む。pysad では `Get` はファイルを読んで
`run(src, 'ffs')` するだけ。Packages は、環境変数 `PYSAD_SAD_PACKAGES` (無ければ pysad の隣の `SAD-master/Packages`) が見つかるときだけ `pkgload.py` が `AutoLoad` と同じく最初に使われたときに読む (CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・D・rational・tableform など。docs/34)。見つからなければ何も登録せず、必要な関数は Python 側の組み込みで動く (PolynomialFit は Fit.n の定義文字列 `_PRELUDE` を起動時に評価する。Python で持つ関数はパッケージの定義で上書きしない `pkgload._protect_builtins`)。Tkinter.n・KBFrame.n・ListPlot.n・Plots.n・Canvas.n も同じ方法で読む (docs/22, 35)。

## 2. 演算子と字句

### 2.1 演算子表 (`tfstk.f` の `opcode`/`iprior`)

`iprior` は優先順位で、**数値が小さいほど強く結合**する。pysad の `_BIN` 表と Pratt パーサの
`parse_expr(rbp)` はこの数値をそのまま使う (「rbp より強く結合する式を読む」)。

| 記号 | iprior | 意味 (head) | 備考 |
|---|---|---|---|
| `@` | 5 | `f@x` = `f[x]` | 右結合 |
| `[`, `[[` | 6 | 呼び出し / `Part` | `a[[i]]` も `a[i]` も Part (§3) |
| `++`, `--` | 7 | Increment / Decrement | 後置・前置 |
| `/@`, `//@`, `@@` | 9 | Map / MapAll / Apply | 右結合 |
| `.` | 10 | Dot (内積・行列積) | 単項 `-` もこの強さ |
| `^` | 15 | Power | 右結合 |
| `*`, `/` | 40 | Times / Divide | |
| `+`, `-` | 50 | Plus / Subtract | |
| `//` | 80 | StringJoin | 数値も文字列化して連結 |
| `==`, `<>`, `>`, `>=`, `<`, `<=` | 100 | Equal / **Unequal** / Greater … | `=>`, `=<`, `><` も同義 (tfetok.f) |
| `===`, `<=>` | 120 | SameQ / UnsameQ | |
| `~` | 150 | Not | 前置 |
| `&&` | 160 | And | |
| `\|\|` | 170 | Or | |
| `\|` | 175 | Alternatives (パターン) | |
| `:` | 180 | Pattern / Optional (`x_:v`) | |
| `->`, `:>` | 190 | Rule / RuleDelayed | |
| `/.`, `//.` | 200 | ReplaceAll / ReplaceRepeated | |
| `+=`, `-=`, `*=`, `/=` | 210 | AddTo … | |
| `&` | 220 | Function (純関数) | 後置 |
| `/:` | 240 | TagSet | |
| `=`, `:=`, `^=`, `^:=` | 250 | Set / SetDelayed / UpSet / UpSetDelayed | 右結合 |
| `=.` | 900 | Unset | |
| `;` | 1000 | CompoundExpression | |
| `,` | 9000 | 引数区切り | |
| `#`, `##` | 3 | Slot / SlotSequence | `#2` は直後の数字 |
| `?` | 4 | PatternTest | |
| `..`, `...` | 172 | Repeated / RepeatedNull | |

`/;` (Condition) は SAD に演算子が無く構文エラーになるので、pysad も実装しない。

### 2.2 字句の規則 (`tfetok.f`)

- 演算子は長いものから照合する (`//.`, `^:=`, `===`, `<=>` → 2 文字 → 1 文字)。ただし `=.`, `/.`, `//.`,
  `++`, `--` の直後が数字なら短い方に戻す (tfetok.f の `case (mtfreplace, mtfunset, …)`)。
  したがって `x=.5` は `x = 0.5`、`a/.5` は `a / 0.5`、`i++1` は `i + +1` である。pysad の `tokenize`
  も同じ分岐を持つ。
- `.` の直後が数字なら数値 (`eval2`)。数値は `1.5`, `.5`, `1e-3`, `1d-3` (Fortran 形式の `d` も可)。
- **暗黙の積**: `2 Pi`, `a b`, `2 (x+1)`, `2Pi` は Times として解釈する (SAD と同じ)。改行をまたぐときは積にしない。
- `!` から行末までがコメント。`(* ... *)` もコメント。文字列は `"..."` (または `'...'`) で、エスケープは tfetok.f の kxmakestring と同じ: `\n \r \t \f \e` (ESC)、`\` + 8 進 1〜3 桁 (`"\101"` は A、256 で折り返す)、`\x` + 16 進 2 桁まで、`\u` + 16 進 4 桁まで、`\` + 改行は行継続、それ以外は `\` を落として次の文字 (`\q` は q)。
  pysad は `strip_comments` で文字列の中を保護しつつ先に落とす。
- `\` + 改行は行継続 (`tremovecont`)。
- シンボルは英字・`$`・`` ` `` で始まり英数字・`$`・`` ` `` が続く。`x_`, `x__`, `x___`, `x_Real`, `_Real`
  はパターン語として一語で読み、pysad では `tokenize_full` が `Pattern[x, Blank[Real]]` 等に分解する。

## 3. 評価の意味論 (pysad の実装)

| 項目 | 実装 |
|---|---|
| 値の型 | 実数 `float` (True=1, False=0)、文字列 `str`、リスト `list`、`Sym` (未定義シンボル)、`Rule`、`Func` (純関数)、`Expr` (未評価式)、`Null`、複素数 (`I`) |
| リスト | 1 始まり。負の添字は末尾から。`a[[i,j]]` も `a[i,j]` も `Part` (値がリストのシンボルに `[ ]` を付けると `apply_value` が Part にする)。`a[[i]]=v` と `a[i]=v` は `set_part`。`a[0]` は head `List` |
| 文字列 | `s[i]` で 1 文字を取り出せる |
| 真偽判定 | `truthy`: 数値なら `!= 0`、それ以外は未定 (`If[c, t, f, u]` の第 4 引数へ) |
| 代入 | `x = v` は値を複製して格納 (`_copy`)、`x := e` は `downvalues` に遅延本体として登録 |
| 関数定義 | `f[x_, y_:2] := body` は `downvalues[f]` に (パターン列, 本体, delayed) で登録。リテラル引数の定義 `Mode["ZX"]={3,1}` はリテラルを `_Lit` として先頭に挿入し、`f[x_]` より先に照合する |
| パターン | `x_` (1 個)、`x__` (1 個以上)、`x___` (0 個以上)、`x_Real`/`_String`/`_List`/`_Symbol`/`_Rule` の型制限、`x_:v` (既定値)、`p?test`、`a\|b`、`{x_, y_}` のリスト分解、リスト `{_,_}`、`f[__]`、`Repeated` (`{_Real..}`)。`x__` は残りのパターン数を残して貪欲に取る |
| 束縛 | `:=` 定義は引数を `Frame` に積んで本体を評価、`=` 定義は本体に `_Lit` を代入してから評価 |
| 記号式 | 未定義シンボルを含む `+ - * / ^` は未評価の式にする。SAD と同じ整列と表示: `x+y+x → x+x+y`、`1+x+3 → 4+x`、`(2 x)`、`(x y)`、`x+(-y)`、`(x y^-1)`、`x^2`、`(-x)`。`Sin[x]` などの数学関数は未定義シンボルを渡すと未評価で返す。規則の右辺 (`Cases[l,{x_,y_}->x+y]`) に式を書く使い方が動く |
| 純関数 | `body&` は `Slot` を引数に置換 (`subst_slots`、内側の `Function` には入らない)。`Function[{x,y}, body]` は名前束縛 |
| 局所変数 | `Module` は字句スコープ (gs と同じ): 宣言した名前を呼び出しごとに `名前$N` へ改名して本体を評価するので、本体から呼ぶ別の関数には見えない (`x=10; show[]:=x; Module[{x=7}, show[]]` は 10)。初期値は外側の環境で評価する (`Module[{x=7, y=x}, y]` の y は外側の x)。値を持たない局所は `t$N` と表示する。定義 (`f[x_]:=…`) や純関数を含む Module は、出た後も局所値を残す。`Block` は動的スコープ (`x=3; f[]:=x; Block[{x=5}, f[]]` は 5)。`{{x,b}=Thread[d]}` のようなリスト分解初期化、`x:=e` の遅延初期化 (`$Delayed`) も可。`With` は値を固定 |
| 反復 | `Do`, `Table`, `Sum`, `Product` の反復子 `{i, a, b, d}` / `{i, list}` / `{n}`。`Sum`/`Product` は複数反復子可。`Table/Do/Sum` の `{i,{a,b,c}}` 形式は SAD でエラー。反復変数は Module の局所名より優先して束縛する |
| 制御 | `If`, `Which`, `Switch`, `While`, `For`, `Return` (例外 `_Return`)、`Break`/`Continue`、`Catch`/`Throw`、`Check[e, fail]`。`Return` は純関数・`Do`・`Scan`・`Module`・`Block` を抜けて、外側のユーザー定義関数 (`:=`) まで届く (gs で確認)。関数の外では `General::unexpbreak` でその文を打ち切る。`Which`/`Switch` は一致なしのとき式を返す (`Which[1>2,"a"]` は未評価のまま。`Switch[a,1,2]` は `Switch[5,1,2]`、引数がすべて定数なら gs が読み込み時に畳むので `Switch[9]`)。引数の個数の誤りは `General::narg`。`For` は 4 引数のときだけ評価し、3 引数は未評価で返す (gs と同じ。本体 (4 番目) は pysad では実行する。§6) |
| 代入演算 | `++`, `--`, `+=` 等は `Set` を経由するので `a[[i]]++` も動く。`AppendTo` |
| 改行と `;` | `split_statements` が括弧の外の `;` と改行で文を切る。行末が演算子 (`+ - * / ^ = , < > & \| : @ . ~`) なら次行へ継続 (`_incomplete`)。ただし `free Q*` のように英数字直後の `*` で終わる行はワイルドカードなので継続しない。括弧の中では改行は空白 |
| 式 head | `f[a][b]`、`(f@g)[x]`、`Function[..][x]` のように head が式でもよい |
| 実行時エラー | 最上位の文での Python 例外 (TypeError など) は SAD の `???General::wrongtype` と同じく、その文を打ち切って続行 (`PYSAD_TRACE=1` で traceback) |

## 4. 入出力と書式

- `Print[a, b, …]` は引数を文字列化して連結し、`PageWidth` (既定 131) で折り返して改行付きで出力。
  `Write[unit, …]` は同じく改行付き、`WriteString[unit, …]` は改行なし。unit 6 は標準出力。
- `OpenWrite[path]`/`OpenAppend[path]`/`OpenRead[path]` はユニット番号 (10 から) を返し、`Close[unit]` か `Close[{u1,u2,…}]`。`$Output` (= -1) への `Write`・`WriteString` は標準出力 (unit 6 と同じ)。`Write[$Output]` (引数 1 つ) は `General::narg`。
  相対パスは `Interp.cwd` (run_sad では最後に渡したスクリプトのディレクトリ)。
  **OpenRead の失敗** (gs に合わせた): 開けなければメッセージを出して入力行を打ち切り、戻りの値に対する `Read[h,…]`・`Close[h]` は wrongtype になる (`test_evalfixes.py`。入力ファイルが無い場合がこれに当たる)。
- `Read[unit, item, opts]` は gs の `tfreadstringf` に合わせた (記録 = 1 行の中の位置を持つ)。`Word` は `WordSeparators` (既定は空白・`,`・タブ・改行) で区切り、語が記録の終わりなら次の記録へ (`ReadNewRecord->False` ならそのまま)、`NullWords->True` で隣り合う区切りの間に `""`。`Real` は語を読んで評価、`String` と `ReadString` は記録の残り、`Character` は 1 文字、`Expression` は現在位置から 1 つの式 (隣り合う項 `x+y  3*4` で切り、括弧が閉じなければ次の行も使う) を評価して返す。`{item,…}` はリスト、`n*item` は頭が Null のリスト (Print で連結される)。EOF は `EndOfFile`。`ReadList`・`Get` は §4 のとおり別。
- `DateString[t]` は t (1900-01-01 0:0:0 からの秒、または `{Y,m,d[,H[,M[,S]]]}`) を `mm/dd/yyyy HH:MM:SS` にする (7 番目の offset は無視、未確認)。`System[cmd]` は実行せず、その旨を表示して 0 を返す (EPICS の `Ca*` も動かさない)。`ToString[x, form]`, `StringJoin`, `StringForm` あり。
- リストの文字列化は `{1,2,3}` のように `,` だけで区切る (空白なし)。文字列はリストの中では `"..."` 付き。
  `PageWidth` を超えると SAD と同じく `,` の後で折り返し、続きの行は 1 空白から始める (`wrap`)。
  入れ子の長いリストの折り返しは SAD (要素の境界で改行、行の上限は pw-1) と pysad (区切りのカンマ、上限 pw) で異なる (PageWidth が十分大きければ一致)。
- `Print[Null]` / `ToString[Null]` は `Null` (リストや関数の引数の中では空)。Expr の表示は `f[a,b]` (カンマ後に空白なし)。
- 複素数の表示: `I`、`(-I)`、`(re+I)`、`(re-I)`、`(re+im I)`、虚部が負なら `(re<虚部の書式> I)`。

### 4.1 `$FORM` と autofg.f (`fmt_real`)

`$FORM` が `""` または `"S"` のときは最短の往復可能表現 (`_autos`: `0.5`, `5.07590922e-16`)。
`"15.10"` のような `幅.小数桁` を指定すると `autofg.f` を直訳した `fmt_real` が使われる。
特徴は、整数部が 0 のとき `0.` を書かず ` .0000000000` とすること、小数桁に収まらない小さな数は
`E` 形式にすること、先頭に 1 文字の符号桁を取ること。`S`/`M`/`C` 前置 (`"S15.10"`) は末尾の 0 と
空白を落とす。`F`/`G`/`E` 前置は Fortran の編集記述子そのまま。

| 値 | `"15.10"` | `"S15.10"` | `""` |
|---|---|---|---|
| 7.4210230634 | `   7.4210230634` | `7.4210230634` | `7.4210230634` |
| 0 | `    .0000000000` | `0` | `0` |
| 5.07590922e-16 | ` 5.07590922E-16` | `5.07590922E-16` | `5.07590922e-16` |
| -0.0012345678901 | `   -.0012345679` | `-.0012345679` | `-0.0012345678901` |

内部では `_strfromd` (18 桁の数字列と 10 の指数)、`_roundnumstr` (指定桁で四捨五入、桁上がりで
指数を 1 増やしてやり直し) を autofg.f と同じ手順で行う。`FFS` 側のマッチング表示 (tffsmatch の
`1PG11.4`) は `fmt_g(x, 11, 4)` で、`0.1 ≤ |x| < 10^4` なら F 形式 + 4 空白、それ以外は `2.5977E-03`
のような 1PE 形式になる。整数に近い値は `1080.` (小数点付き)。

## 5. SAD 独自の関数

### 5.1 数値・行列・乱数

- **PolynomialFit[data, nord]**: `Packages/Fit.n` の定義を `_PRELUDE` としてそのまま評価する。
  データの x を平均 `xa` と `scale=1/Max[Abs[x-xa]]` で正規化し `LinearSolve` で最小二乗解を求めてから、
  係数を元のスケールに戻す。ただし `Abs[xa]*scale < 0.1` なら `xa=0` とし、そのとき
  `xl=If[xa, Table[…], 1]` が**スカラー 1** になる (リストでない) ので、`cm = Table[...]/xl` の除算も
  スカラー除算になる、という SAD の癖も含めて再現する。戻り値は `{{c0, c1, …}, {Residual -> r}}`。
  numpy 版 `PolynomialFit$numpy` も残してあるが既定では使わない。
- **Spline[{{x,y},…}]** は scipy の `CubicSpline` を包んだ純関数を返す (境界条件は確認済み: 自然スプライン `Spline[{{x,y},..}]` と `Spline[..,0,0]` の形・`Sort[list, Greater]`・`Scan` が gs と全桁一致。`Spline[data,0,0]` は gs も pysad も未評価のまま)。
- **複素行列の演算** (gs に合わせてある): `Dot`・`Inverse` の結果は、gs が**要素ごと**に虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻す。`DiagonalMatrix` は複素数の列を受けて動く (`LDSolve.n` が動く)。`LDSolve` の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (`test_evalfixes.py` 2 件)。`Packages/Optics.n` の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo` は gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し pysad は 0)。`PeriodicR` は呼び出し方が不明で比較できなかった (gs も引数エラー)。
- **Eigensystem** (`{値, ベクトル}`): 実行列は teigen.f の直訳 (`pysad/teigen.py`) で解く。固有ベクトルは単位ベクトルではなく
  SAD 固有の規格化、複素固有値の対は 行 k = 実部+i虚部 / 行 k+1 = 実部-i虚部。実行列の対称・縮退でも虚部の雑音が出ない。
  ただし n>8 の実対称行列は `np.linalg.eigh` (SAD の teigen は縮退で丸め依存)。n<=8 や非対称は teigen の直訳。
- **LinearSolve** (最小二乗 `lstsq`)、`Inverse`, `Det`, `Dot`, `Transpose`, `DiagonalMatrix`, `IdentityMatrix`, `SingularValues`。
- **Thread**: リストのリストは転置、`Rule` のリストは要素ごとの規則、式は引数リストで分配。
- **Override[rules…]**: 規則を平坦化し、同じ左辺は最初のものだけ残す (オプション処理用)。
- **Random/GaussRandom/SeedRandom**: ffscmd.py が `sadrandom.SadRandom` (random_plugin_sad.c の
  LCG `seed ← 48828125·seed mod 2^32`、初期種 17、Box–Muller、GCUT 棄却) を割り当てるので、
  追跡の FLUC と同じ乱数列を共有する。`Random[n]` は n 個のリスト (`Random[0]` は `{}`、リスト引数は wrongtype、負は wrongnum)、
  `Random[n,m,…]` は n×m×… の配列 (行優先に順に引く。小数は切り捨て。`GaussRandom`・`ParabolaRandom` も同じ)。`GaussRandom[]` はスカラー。
  `SeedRandom[n]` は種、`SeedRandom[]` は `{"SAD", seed, filled, next の下位語, 上位語}`。`SeedRandom["SAD"]` と、先頭が文字列のリストは乱数の種類の選択だけで
  状態は変わらない (`SeedRandom[状態]` で状態は復元されない: gs の tfRandom_.c はリストの数を読まない。数だけのリストは種 1 になる)。"SAD" 以外は `System::error`。
- **StringMatchQ**: ワイルドカードは `*` (0 個以上)、`%` (1 文字)、`{abc}`・`{^abc}`、`a|b` (SADHelp の wildcards)。`?` は普通の文字、大文字小文字は区別、空文字列は何にも合わない (gs で確認)。

### 5.2 リスト・文字列・行列の関数 (`pysad/sadlist.py`)

SAD の実装 (組み込み + Packages/functions.n, Misc.n) の挙動に合わせてある: `MapThread MapAt ScanThread Complement Union Intersection DeleteCases Cases MatchQ Extract Position Count MemberQ FreeQ Select
Sort Reverse Join Flatten Partition Transpose Take Drop Delete Insert ReplacePart Append Prepend Apply Map MapAll Scan Fold FoldList Nest NestList FixedPoint NearlySameQ Through Tr RotateLeft/Right Depth Level AtomQ IntegerQ
EvenQ OddQ Positive.. Binomial FractionalPart Round Floor Ceiling Mod Max Min Length Dimensions Inverse Det Dot IdentityMatrix Outer Thread StringReplace ToExpression FromCharacterCode ToCharacterCode Order First Last Rest`。

- レベル指定 (n → {1,n}, {n}, {m,n}, Infinity, 負は高さ)、後順の走査。**Position だけ頭部 `{0}` と全体 `{}` も走査し既定のレベルは {0,Infinity}** (Count/MemberQ/FreeQ はこれの上に定義されている)。
- 標準順序は 数 < 文字列 < シンボル < 複合式 < リスト (リストは長さ優先、文字列は `a<A<aa<b<B`)。
- `Take[l,5]`/`Take[l,-5]`/`Take[l,{2,5}]` は範囲外で `{}`、`Drop[l,{2,5}]` は存在する分だけ落とす。
- **SAD に無い関数** (未評価のまま返る): `Total Cross Norm Riffle StringSplit Mean Variance Most Accumulate StringTake IntegerPart Eigenvalues Ordering Tally DeleteDuplicates` など。
  pysad の拡張 (Total, Most, Accumulate, StringTake, IntegerPart, Eigenvalues) は残してあるが、SAD のスクリプトは使えない前提。
- SAD の癖で再現していないもの: `Level[..,{-1}]` のバグ、`DeleteCases` の第 4 引数、非正方・ゼロ行列の `Inverse`、`Insert[l,x,0]`。

### 5.3 FFS 専用関数 (`pysad/sadelem.py`, ffscmd.py)

- `Element` (ファミリー単位、単一マッチはスカラー、リスト指定はリスト、マッチ無しは `{}`、`+= -= *= /=`)、`SetElement` (`{name, TYPE, {keyword->value..}}`、キーワード表は `pysad/sadkw.py`)、`FLAG` (大文字のみ)、
  `ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine` (`-BeamLine` は逆順・反転、`n*BeamLine` は繰り返し、`Fill$String` の整形、デバッグ行 `WriteBeamLine {{..},..}` も SAD と同じ)、`GetMAIN`、`CalculateOptics`。
- `Twiss[key,pattern]` / `LINE[key,pattern]`: 一致 1 件ならスカラー、0 件は `{}`、複数はリスト。
- 名前パターン (tmatch.f): `* % {abc} {^abc} |`。`?` は普通の文字。
- Twiss / LINE / FFS / Emittance / TrackParticles などの仕様は 03_ffs_commands §7。

### 5.4 11/13 以降に足した関数群 (詳細は各文書)

| 分野 | 内容 | 文書 |
|---|---|---|
| 特殊関数・数学関数 | Bessel・超幾何・ガンマ族・ゼータ・多重対数・D.n の特殊関数 (458 式が gs と一致)。K64 の `InverseErf` の誤差も再現。基本関数 (Log・ArcSin・Power・Erf・Factorial など) の複素・範囲外の扱いを gs に合わせた | docs/33 |
| Fit・FindRoot・NIntegrate・FitGaussian・DownhillSimplex | 数値微分で gs と 5e-10〜1e-14。ベクトル独立変数・複素 FindRoot は gs が誤値 (B52, B61) で照合なし | docs/33 |
| ファイル・文字列・システム関数 | FileQ・FileNames・CopyFile・ToDate・GetEnv など (Python の os・shutil で実装)、ARRAY・Len・StringInsert・Diagonal・Expand など | docs/33 |
| 記号微分・パッケージ | `D[f,x]` は D.n を SADScript のまま実行 (59 式が gs と文字列まで一致)。CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・Rational・TableForm | docs/34 |
| Class | `Class[親の列, クラス変数, インスタンス変数, メソッド…]`、`x@メンバー`、`Constructor`/`Destructor`、継承 (`sadclass.py`) | docs/22 |
| Tk・描画 | `TclEval`・`TkEventLoop` など Tcl との接点、`TkCanvasCreateItem`・`CanvasClipLine`・`Vector` (ListPlot・Plot・OpticsPlot の下地) | docs/22, 35 |
| CSR・GaussianCoulomb・BBBrem・RadiationField | 物理計算の関数 | docs/32, 34 |

評価器の側で gs に合わせた主な点 (試験は `test_evalfixes.py`・`test_gui_lang.py`・`test_symbolic.py` ほか): 連鎖した比較 `a<x<b`、`Null[a,b]` (空の Sequence) の展開、UpValue (`^:=`・`a/:f[a]:=…`・SubValue `f[a][b]`)、保持属性 (HoldAll など) の引数への置き換え、`//.` の不動点までの反復、Plus・Times の Flat・Orderless 照合、`Part[l,{i,j},k]`、`With[{f[x_]:=…},…]`、`StandardForm`、数の既定の書式 (`1e6`→`1e6`、`1e10`→`1e+10`)、`Skip`・`Read` の型。`ElectronMass` は 510998.95 eV。

## 6. 対応していない構文・関数

未対応 (parse はできるが未実装、または parse できない):

- `Condition` (`/;`: SAD に演算子が無い。`Condition[1,2]` は未評価のまま返ることを確認)、`Literal`、`Hold` 系の厳密な保持、`Message`、`Definition`、
  `SetAttributes`/`Protect` (無視)、関数に対する `Options[f]` (常に `{}`。Class のインスタンスの `Options` は sadclass が扱う)、`DefineOptions`
- `Class`・`^=` による UpValue・`//.` の反復・`ReplaceAll` の構造パターン・`Module` と `Block` の区別は対応済み (§3、§5.4、docs/22, 32)。ただし `Unevaluated` の扱いは `Unevaluated$[…]` の形の一部 (docs/22 §4)、`^:=` のパターン付き UpValue の一般形、`Hold` 系の厳密な保持は確かめていない (未確認)。
- 表示系 `Format`, `TableForm` (SAD の tableform.n は Packages で動く: docs/34)。`TkWait`・`Update`・`Pause` は sadscript の組み込みでは no-op (tclbridge が有効なときの挙動は未確認)
- `System[]` の実行 (表示して 0 を返す)、EPICS (`CaGet` など)。ソース確認: `CSR.n` の `TrackAnimate` は `System["import -silent -window …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶので Windows の pysad では画像が作れない。`Wait[mem]` は `FreeMem[]` (top) と `Sleep` を使い、gs はこの環境で top が無く `FreeMem` が OpenRead エラーになる。pysad の `FreeMem` は psutil で代替して 1.0 を返す。並列を増やしても速くなるかは未確認 (pysad は fork しない)
- **Tk ウィジェットの値 (差あり、未修正)**: `w=Window[]`・`b=Button[w,…]` の直後、gs では `Head[w]=Widget` (`Widget[w,Window[Root$Window]]`) で `e[Insert]={0,"AB"}`・`e[Get]`・`e[Delete]` が働くが、pysad では w・b・e が裸のシンボルのままで、メソッド呼び出しは何もせず (Tcl 呼び出し 0 回、エラーもなし)、親が toplevel でなく Root になる。`Tkinter.n` の `DefineWidget` 内 `a$:=Widget[a$,w[master$]]` が値にならないところまで絞ったが、gs が再帰しない理由 (属性・評価の順) は未特定。パターン変数への代入 (`f[x_]:=(x=5); f[aa]`) は gs と同じ。深い再帰では pysad は Python の RecursionError で落ちる (gs は `General::deep`)。`Insert` の位置に `"End"` (大文字) を渡すと gs は Tcl エラー 1501 (docs/22・39)
- **`BeginPackage`/`Begin` の文脈 (名前空間): 未対応** (gs の動作は調べた: 読み込み時に新しい記号を ``Ctx`名`` に作り、``Stat`Mean2[..]`` の形で呼べ、`Context[]`・`$Context`・`$ContextPath` の値と ``Names["Stat`*"]`` (= ``{"Stat`"}``) の表示が独特。pysad は文脈を無視して全部大域)。読み込み時の名前解決と `$ContextPath` の模倣が要り、SAD の `.n` パッケージ (文脈つきの名前を持つもの) の読み込みに影響するので見送った。
- **`For` の本体 (4 番目の引数)**: gs は実行しない (`Packages/functions.n` の定義 ``For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)`` どおりなら実行するはずで、同じ式を別名で定義すると gs も実行する。原因は未特定の K64 の読み込み時の癖) が、pysad は定義どおり実行する (差あり、意図的)。
- `Check[Return[1], 2]` は gs では Return を失敗として捕まえて 2 を返す。pysad は Return を通す。`Thread[a==b]` (a が値、b が未定義) は gs が `3==b` を返すが、pysad は `==` が 0 を返すので Thread の引数にならない。リストと記号の ==・数どうしでない `==` の未評価形は一部だけ (`{1,2}==1` は未評価、記号どうしは 0)。
- `Out[n]:=` の表示は、最上位 (run_sad_string) の `;` なしの式文だけ。未定義の記号 1 つだけの文での gs の表示 (前の Out の値を繰り返す) は未対応。`%%` は gs でも記号のまま。
- SAD に無い関数 (`Mean` `Cross` `StringSplit` など) は未評価のまま返る (SAD と同じ。`Total` は pysad の拡張で動く)
- FFS 側: `MCOD/TCOD/MON` の COD 補正は K64 でも式として評価されて何もしない。`Steer`/`Monitor`/`CorrectOrbit` は `Packages/correction.n` を SADScript のまま実行して gs と一致 (docs/34。`PYSAD_SAD_PACKAGES` が必要)。FFS コマンド・関数の対応状況と未対応の一覧は 03_ffs_commands §8

## 7. テスト方法と確認

単体では次のように使う (FFS コマンドなしの `'script'` モード)。

```python
import io
from pysad.sadscript import Interp
o = io.StringIO()
I = Interp(out=o)
I.run('f[x_,y_:2]:=x+y; a={1,2,3}; Print[f[1], " ", a[2], " ", 3<>4, " ", "a"//"b"]', 'script')
print(o.getvalue())      # → 3 2 1 ab
```

ラティス付きのスクリプトは `python run_sad.py lattice.sad script.sad` (SAD の `gs` に複数ファイルを
連結して渡すのと同じ) で、`ffscmd.run_sad` が MAIN 部を parser.py に、FFS 部を `Interp.run(src, 'ffs')`
に渡す。`pysad/tests/` (pytest 578 件を収集、2026-11-17 に 575 通過・3 skip) と `run_regression.py` (50 ベンチをすべて実行して FAIL 0、205 s。QUICK=1 は遅い 4 件 7・12t・13a・16a を SKIP して 46 件で、その 4 件も単独で実行して PASS) で回帰を取る。

SAD のベンチでの確認 (docs/history/00_RESUME.md、<work-olddocs>/sad_level3_4_benchmark_results.md Level 14 節):

- `bench_level13a.sad` / `bench_level13b.sad` (FIT/FREE/GO + Print/LINE/Twiss/MatchingResidual) を
  そのまま実行し、反復行・`Matched` 行・条件表・Print の出力が SAD (`*.result`) と同じ (K1 ≤1e-10、
  残差 5.0763e-16 vs 5.0759e-16)。
- `bench_level14a.sad` (= 20 セルのラティス + 色収差スキャンの台本: Get, OpenWrite, Table, Do,
  DP0 スキャン + `FFS["cell calc;"]`, PolynomialFit) で 21 点のチューンが ≤2e-8、PolynomialFit の係数は
  SAD の式どおり (相対差 ≤1e-5)。
- Element / BeamLine 系は bench_level24f (164 タグ一致)。
