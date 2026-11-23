# 05. SADScript 言語と pysad のインタプリタ

最終更新: 2026-11-23 (現状のスナップショット)

**この文書で分かること**: SADScript (SAD に内蔵された式言語) の文法と評価の規則、pysad のインタプリタ (`sadscript.py`) がどこまで対応しているか、対応していない構文、SAD の `Fork`・文脈 (`BeginPackage`) などの特殊な機能の扱い。
**前提**: Mathematica 風の式言語を読める。FFS (SAD の対話コマンド部) の基本は [03_ffs_commands.md](03_ffs_commands.md) §1 を先に見ると、MAIN 部・FFS 部の区別が分かる。
**読む順**: §1 (全体と使い方) → §2・§3 (文法と意味) → §4 (入出力・書式) → §5 (SAD 独自の関数) → §6 (対応していないもの・特殊な機能)。特殊関数・Fit・SAD の Packages は [06_sadscript_library.md](06_sadscript_library.md)。

対象のコード: [`pysad/sadscript.py`](../pysad/sadscript.py) (字句解析・Pratt 構文解析・評価器、約 4800 行)、[`pysad/sadlist.py`](../pysad/sadlist.py) (リスト・文字列・行列の関数)、
[`pysad/sadelem.py`](../pysad/sadelem.py) (FFS 専用関数)、[`pysad/teigen.py`](../pysad/teigen.py) (Eigensystem)、[`pysad/sadclass.py`](../pysad/sadclass.py) (Class)、[`pysad/pkgload.py`](../pysad/pkgload.py) (Packages/*.n の読み込み)、[`pysad/specfn.py`](../pysad/specfn.py)・`fitroot.py`・`miscfn.py`・`fileops.py` (特殊関数・Fit・ファイル関数など)。SAD 側の根拠は
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
評価する ([03_ffs_commands.md](03_ffs_commands.md) §1)。pysad は判定の向きが逆で (コマンドでなければ式)、`Interp.run(src, 'ffs')` が文ごとに
`_is_ffs_command` (先頭語が `ffs_words` に登録されたコマンドか、かつ続きが `=`, `[`, `(` などで
始まらない) を見て、コマンドなら `ffs_handler` (ffscmd.py の `FFS.execute`) へ文字列のまま渡し、
それ以外は `parse` → `eval` する。

モジュールは `Get["file"]` で読み込む。SAD は起動時に `Packages/init.n` を読み、`functions.n`,
`MessageName.n` を `Get` し、`Fit.n` (`PolynomialFit`, `FindRoot`) や `Matrix.n`, `Optics.n` などは
`AutoLoad[symbol, file]` で最初に使われたときに読む。pysad では `Get` はファイルを読んで
`run(src, 'ffs')` するだけ。Packages は、環境変数 `PYSAD_SAD_PACKAGES` (無ければ pysad の隣の `SAD-master/Packages`) が見つかるときだけ `pkgload.py` が `AutoLoad` と同じく最初に使われたときに読む (CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・D・rational・tableform など。[06_sadscript_library.md](06_sadscript_library.md) §5)。見つからなければ何も登録せず、必要な関数は Python 側の組み込みで動く (PolynomialFit は Fit.n の定義文字列 `_PRELUDE` を起動時に評価する。Python で持つ関数はパッケージの定義で上書きしない `pkgload._protect_builtins`)。Tkinter.n・KBFrame.n・ListPlot.n・Plots.n・Canvas.n も同じ方法で読む ([20_tcltk_gui.md](20_tcltk_gui.md)、[04_ffs_small_commands.md](04_ffs_small_commands.md))。

### 1.1 動かし方 (最小の例)

1. 式だけ試す (FFS コマンドなしの `'script'` モード):

```python
import io
from pysad.sadscript import Interp
o = io.StringIO()
I = Interp(out=o)
I.run('f[x_,y_:2]:=x+y; a={1,2,3}; Print[f[1], " ", a[2], " ", 3<>4, " ", "a"//"b"]', 'script')
print(o.getvalue())      # → 3 2 1 ab
```

2. ラティス付きのスクリプトは `py -3.12 run_sad.py lattice.sad script.sad` (SAD の `gs` に複数ファイルを連結して渡すのと同じ)。`ffscmd.run_sad` が MAIN 部を parser.py に、FFS 部を `Interp.run(src, 'ffs')` に渡す。Packages を使うときは先に環境変数 `PYSAD_SAD_PACKAGES` に SAD の `Packages` のフォルダを指定する。

## 2. 演算子と字句

### 2.1 演算子表 (`tfstk.f` の `opcode`/`iprior`)

次の表は、演算子の結合の強さ (`iprior`: **数値が小さいほど強く結合**) と意味を示す。pysad の `_BIN` 表と Pratt パーサの
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

次の表は、値の型・リスト・代入・関数定義・パターン・局所変数・反復・制御などの評価規則をまとめたもの。

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
| 局所変数 | `Module` は字句スコープ (SAD と同じ): 宣言した名前を呼び出しごとに `名前$N` へ改名して本体を評価するので、本体から呼ぶ別の関数には見えない (`x=10; show[]:=x; Module[{x=7}, show[]]` は 10)。初期値は外側の環境で評価する (`Module[{x=7, y=x}, y]` の y は外側の x)。値を持たない局所は `t$N` と表示する。定義 (`f[x_]:=…`) や純関数を含む Module は、出た後も局所値を残す。`Block` は動的スコープ (`x=3; f[]:=x; Block[{x=5}, f[]]` は 5)。`{{x,b}=Thread[d]}` のようなリスト分解初期化、`x:=e` の遅延初期化 (`$Delayed`) も可。`With` は値を固定 |
| 反復 | `Do`, `Table`, `Sum`, `Product` の反復子 `{i, a, b, d}` / `{i, list}` / `{n}`。`Sum`/`Product` は複数反復子可。`Table/Do/Sum` の `{i,{a,b,c}}` 形式は SAD でエラー。反復変数は Module の局所名より優先して束縛する |
| 制御 | `If`, `Which`, `Switch`, `While`, `For`, `Return` (例外 `_Return`)、`Break`/`Continue`、`Catch`/`Throw`、`Check[e, fail]`。`Return` は純関数・`Do`・`Scan`・`Module`・`Block` を抜けて、外側のユーザー定義関数 (`:=`) まで届く (SAD と同じ)。関数の外では `General::unexpbreak` でその文を打ち切る。`Which`/`Switch` は一致なしのとき式を返す (`Which[1>2,"a"]` は未評価のまま。`Switch[a,1,2]` は `Switch[5,1,2]`、引数がすべて定数なら SAD が読み込み時に畳むので `Switch[9]`)。引数の個数の誤りは `General::narg`。`For` は 4 引数のときだけ評価し、3 引数は未評価で返す (k64 モード。本体 (4 番目) の扱いは v2 と k64 で異なる。§6 の `For`) |
| 代入演算 | `++`, `--`, `+=` 等は `Set` を経由するので `a[[i]]++` も動く。`AppendTo` |
| 改行と `;` | `split_statements` が括弧の外の `;` と改行で文を切る。行末が演算子 (`+ - * / ^ = , < > & \| : @ . ~`) なら次行へ継続 (`_incomplete`)。ただし `free Q*` のように英数字直後の `*` で終わる行はワイルドカードなので継続しない。括弧の中では改行は空白 |
| 式 head | `f[a][b]`、`(f@g)[x]`、`Function[..][x]` のように head が式でもよい |
| 実行時エラー | 最上位の文での Python 例外 (TypeError など) は SAD の `???General::wrongtype` と同じく、その文を打ち切って続行 (`PYSAD_TRACE=1` で traceback) |

## 4. 入出力と書式

- `Print[a, b, …]` は引数を文字列化して連結し、`PageWidth` (既定 131) で折り返して改行付きで出力。
  `Write[unit, …]` は同じく改行付き、`WriteString[unit, …]` は改行なし。unit 6 は標準出力。
- `OpenWrite[path]`/`OpenAppend[path]`/`OpenRead[path]` はユニット番号 (10 から) を返し、`Close[unit]` か `Close[{u1,u2,…}]`。`$Output` (= -1) への `Write`・`WriteString` は標準出力 (unit 6 と同じ)。`Write[$Output]` (引数 1 つ) は `General::narg`。
  相対パスは `Interp.cwd` (run_sad では最後に渡したスクリプトのディレクトリ)。
  **OpenRead の失敗** (SAD と同じ): 開けなければメッセージを出して入力行を打ち切り、戻りの値に対する `Read[h,…]`・`Close[h]` は wrongtype になる ([`test_evalfixes.py`](../pysad/tests/script_tests/test_evalfixes.py)。入力ファイルが無い場合がこれに当たる)。
- `Read[unit, item, opts]` は SAD の `tfreadstringf` と同じ規則 (記録 = 1 行の中の位置を持つ)。`Word` は `WordSeparators` (既定は空白・`,`・タブ・改行) で区切り、語が記録の終わりなら次の記録へ (`ReadNewRecord->False` ならそのまま)、`NullWords->True` で隣り合う区切りの間に `""`。`Real` は語を読んで評価、`String` と `ReadString` は記録の残り、`Character` は 1 文字、`Expression` は現在位置から 1 つの式 (隣り合う項 `x+y  3*4` で切り、括弧が閉じなければ次の行も使う) を評価して返す。`{item,…}` はリスト、`n*item` は頭が Null のリスト (Print で連結される)。EOF は `EndOfFile`。`ReadList`・`Get` は上の §1 のとおり別。
- `DateString[t]` は t (1900-01-01 0:0:0 からの秒、または `{Y,m,d[,H[,M[,S]]]}`) を `mm/dd/yyyy HH:MM:SS` にする (7 番目の offset は無視、未確認)。`System[cmd]` は、既定では実行しない (下の §6。EPICS の `Ca*` も動かさない)。`ToString[x, form]`, `StringJoin`, `StringForm` あり。
- リストの文字列化は `{1,2,3}` のように `,` だけで区切る (空白なし)。文字列はリストの中では `"..."` 付き。
  `PageWidth` を超えると SAD と同じく `,` の後で折り返し、続きの行は 1 空白から始める (`wrap`)。
  入れ子の長いリストの折り返しは SAD (要素の境界で改行、行の上限は pw-1) と pysad (区切りのカンマ、上限 pw) で異なる (PageWidth が十分大きければ一致)。
- `Print[Null]` / `ToString[Null]` は `Null` (リストや関数の引数の中では空)。Expr の表示は `f[a,b]` (カンマ後に空白なし)。
- 複素数の表示: `I`、`(-I)`、`(re+I)`、`(re-I)`、`(re+im I)`、虚部が負なら `(re<虚部の書式> I)`。

### 4.1 数値の書式: `$FORM` と `fmt_real`

`$FORM` (数の出力書式を決める SAD の変数) が `""` または `"S"` のときは最短の往復可能表現 (`_autos`: `0.5`, `5.07590922e-16`)。
`"15.10"` のような `幅.小数桁` を指定すると `autofg.f` を直訳した `fmt_real` が使われる。
特徴は、整数部が 0 のとき `0.` を書かず ` .0000000000` とすること、小数桁に収まらない小さな数は
`E` 形式にすること、先頭に 1 文字の符号桁を取ること。`S`/`M`/`C` 前置 (`"S15.10"`) は末尾の 0 と
空白を落とす。`F`/`G`/`E` 前置は Fortran の編集記述子そのまま。次の表は、同じ値を三つの `$FORM` で出したときの表示。

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

- **PolynomialFit[data, nord]** (B18): `Packages/Fit.n` の定義を `_PRELUDE` としてそのまま評価する。**v2 モードは x の平均が 0 のときの係数の誤りを直す** (下の `xl=1` の癖は k64 モードの動作)。例: `{{-1,1},{0,0},{1,1},{2,4.1},{-2,4}}` の 2 次で k64 は {1.0257, 2.0486, 1.0143}、v2 は {-.0085714, .02, 1.0142857} (numpy.polyfit と 1e-13 で一致。`test_b18_polynomialfit_mean_zero_v2`・`_k64_and_nonzero_mean`)。平均が 0 でないデータは両モードで同じ。v2 の判定は `Compat$V2[]` (互換モードが v2 なら 1) を prelude が使う。以下は k64 モードの説明。
  データの x を平均 `xa` と `scale=1/Max[Abs[x-xa]]` で正規化し `LinearSolve` で最小二乗解を求めてから、
  係数を元のスケールに戻す。ただし `Abs[xa]*scale < 0.1` なら `xa=0` とし、そのとき
  `xl=If[xa, Table[…], 1]` が**スカラー 1** になる (リストでない) ので、`cm = Table[...]/xl` の除算も
  スカラー除算になる、という SAD の癖も含めて再現する。戻り値は `{{c0, c1, …}, {Residual -> r}}`。
  numpy 版 `PolynomialFit$numpy` も残してあるが既定では使わない。
- **Spline[{{x,y},…}]** は scipy の `CubicSpline` を包んだ純関数を返す (自然スプライン)。`Spline[data,0,0]` は SAD も pysad も未評価のまま。
- **複素行列の演算**: `Dot`・`Inverse` の結果は、SAD が**要素ごと**に虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻す。`DiagonalMatrix` は複素数の列を受けて動く (`LDSolve.n` が動く)。`Packages/Optics.n` の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo` も動く (SAD は 1e-16 の雑音の要素を出し、pysad は 0)。`PeriodicR` は `PeriodicOptics` を先に呼ぶ必要がある ([06_sadscript_library.md](06_sadscript_library.md) §5)。
- **Eigensystem** (`{値, ベクトル}`): 実行列は teigen.f の直訳 ([`pysad/teigen.py`](../pysad/teigen.py)) で解く。固有ベクトルは単位ベクトルではなく
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
  状態は変わらない (`SeedRandom[状態]` で状態は復元されない: SAD の tfRandom_.c はリストの数を読まない。数だけのリストは種 1 になる)。"SAD" 以外は `System::error`。
- **StringMatchQ**: ワイルドカードは `*` (0 個以上)、`%` (1 文字)、`{abc}`・`{^abc}`、`a|b` (SADHelp の wildcards)。`?` は普通の文字、大文字小文字は区別、空文字列は何にも合わない。

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

- `Element` (ファミリー単位、単一マッチはスカラー、リスト指定はリスト、マッチ無しは `{}`、`+= -= *= /=`)、`SetElement` (`{name, TYPE, {keyword->value..}}`、キーワード表は [`pysad/sadkw.py`](../pysad/sadkw.py))、`FLAG` (大文字のみ)、
  `ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine` (`-BeamLine` は逆順・反転、`n*BeamLine` は繰り返し、`Fill$String` の整形、デバッグ行 `WriteBeamLine {{..},..}` も SAD と同じ)、`GetMAIN`、`CalculateOptics`。
- `Twiss[key,pattern]` / `LINE[key,pattern]`: 一致 1 件ならスカラー、0 件は `{}`、複数はリスト。
- **`Twiss["DETR", …]`** (結合行列 R の行列式 det R = R1·R4 − R2·R3): CELL の CALC のあとだけでなく **INS の CALC のあとも値が入る** (両モード)。小さな結合リング (四極を 45 度回した 4 要素の FODO 風のリング) で INS の CALC のあと、位置 2 以降は DETR = R1R4 − R2R3 (例: 0.00213…)、入口の位置 1 は R1〜R4 が全て 0 で DETR も 0 (2026-11-23 に pysad で確認)。DETR ≥ 1 (モードが入れ替わる) の INS は未確認 ([07_calc_optics.md](07_calc_optics.md))。
- 名前パターン (tmatch.f): `* % {abc} {^abc} |`。`?` は普通の文字。
- Twiss / LINE / FFS / Emittance / TrackParticles などの仕様は [03_ffs_commands.md](03_ffs_commands.md) §7。

### 5.4 追加の関数群 (詳細は各文書)

次の表は、この文書に書かなかった関数群の一覧と、詳しい説明のある文書への案内。

| 分野 | 内容 | 文書 |
|---|---|---|
| 特殊関数・数学関数 | Bessel・超幾何・ガンマ族・ゼータ・多重対数・D.n の特殊関数。K64 の `InverseErf` の誤差は k64 モードで再現。基本関数 (Log・ArcSin・Power・Erf・Factorial など) の複素・範囲外の扱いも SAD に合わせてある | [06_sadscript_library.md](06_sadscript_library.md) §1 |
| Fit・FindRoot・NIntegrate・FitGaussian・DownhillSimplex | 数値微分で計算。ベクトル独立変数・複素 FindRoot は K64 が誤値 (B52, B61) | 06 §2・§3 |
| ファイル・文字列・システム関数 | FileQ・FileNames・CopyFile・ToDate・GetEnv など (Python の os・shutil で実装)、ARRAY・Len・StringInsert・Diagonal・Expand など | 06 §3 |
| 記号微分・パッケージ | `D[f,x]` は D.n を SADScript のまま実行。CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・Rational・TableForm | 06 §5 |
| Class | `Class[親の列, クラス変数, インスタンス変数, メソッド…]`、`x@メンバー`、`Constructor`/`Destructor`、継承 (`sadclass.py`) | [20_tcltk_gui.md](20_tcltk_gui.md) |
| Tk・描画 | `TclEval`・`TkEventLoop` など Tcl との接点、`TkCanvasCreateItem`・`CanvasClipLine`・`Vector` (ListPlot・Plot・OpticsPlot の下地) | 20、[21_graphics.md](21_graphics.md)、[04_ffs_small_commands.md](04_ffs_small_commands.md) |
| CSR・GaussianCoulomb・BBBrem・RadiationField | 物理計算の関数 | [06_sadscript_library.md](06_sadscript_library.md) §4、[14_beambeam.md](14_beambeam.md)、[15_csr.md](15_csr.md)、[18_radlight.md](18_radlight.md) |

評価器の側で SAD に合わせた主な点 (試験は [`test_evalfixes.py`](../pysad/tests/script_tests/test_evalfixes.py)・[`test_gui_lang.py`](../pysad/tests/script_tests/test_gui_lang.py)・[`test_symbolic.py`](../pysad/tests/script_tests/test_symbolic.py) ほか): 連鎖した比較 `a<x<b`、`Null[a,b]` (空の Sequence) の展開、UpValue (`^:=`・`a/:f[a]:=…`・SubValue `f[a][b]`)、保持属性 (HoldAll など) の引数への置き換え、`//.` の不動点までの反復、Plus・Times の Flat・Orderless 照合、`Part[l,{i,j},k]`、`With[{f[x_]:=…},…]`、`StandardForm`、数の既定の書式 (`1e6`→`1e6`、`1e10`→`1e+10`)、`Skip`・`Read` の型。`ElectronMass` は 510998.95 eV。

## 6. 対応していない構文・関数と特殊な機能

### 6.1 未対応

parse はできるが未実装、または parse できないもの:

- `Condition` (`/;`: SAD に演算子が無い。`Condition[1,2]` は未評価のまま返ることを確認)、`Literal`、`Hold` 系の厳密な保持、`Message`、`Definition`、
  `SetAttributes`/`Protect` (無視)、関数に対する `Options[f]` (常に `{}`。Class のインスタンスの `Options` は sadclass が扱う)、`DefineOptions`
- `Class`・`^=` による UpValue・`//.` の反復・`ReplaceAll` の構造パターン・`Module` と `Block` の区別は対応済み (§3、§5.4、[20_tcltk_gui.md](20_tcltk_gui.md))。ただし `Unevaluated` の扱いは `Unevaluated$[…]` の形の一部 ([20_tcltk_gui.md](20_tcltk_gui.md) §4)、`^:=` のパターン付き UpValue の一般形、`Hold` 系の厳密な保持は確かめていない (未確認)。
- 表示系 `Format`, `TableForm` (SAD の tableform.n は Packages で動く: [06_sadscript_library.md](06_sadscript_library.md) §5)。`TkWait`・`Update`・`Pause` は sadscript の組み込みでは no-op (tclbridge が有効なときの挙動は未確認)
- `System[]`: **既定では実行しない** (任意のシェルコマンドを実行できるため)。v2 モードは 1 行の案内を出して `$Failed` を返し、`[ -e path ]`・`[ -d path ]`・`[ -f path ]`・`[ -r path ]`・`[ -s path ]`・`test …` の形だけはシェルを使わずに判定して 0 (真)・1 (偽) を返す (KBFrame などの `System["[ -d dir ]"]==0` のため)。k64 モードは「実行しません」を出して 0。許可: 環境変数 `PYSAD_ALLOW_SYSTEM=1` (pysad を起動する前に設定)、または Python から `pysad.miscfn.allow_system(True)` (SADScript の中からは変えられない)。許可があるとき: シェルで実行し (作業ディレクトリ `Directory[]`、標準入力は空。Windows は cmd.exe、ほかは /bin/sh、`PYSAD_SYSTEM_SHELL` に sh・bash を書けば `<shell> -c cmd`)、戻り値は v2 が終了コード (成功 0)、k64 が 1 (K64 は起動できれば常に 1)。確かめた範囲は Windows (cmd.exe) と WSL (sh) の `echo … > file`・`exit 3`・標準出力で、SAD の Packages が呼ぶ個々のコマンドは未確認 (`test_system_is_disabled_by_default`、`test_system_runs_only_with_permission`)。EPICS (`CaGet` など) は未対応。ソース確認: `CSR.n` の `TrackAnimate` は `System["import -silent -window …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶので Windows の pysad では画像が作れない。`Wait[mem]` は `FreeMem[]` (top) と `Sleep` を使い、K64 (gs) はこの環境で top が無く `FreeMem` が OpenRead エラーになる。pysad の `FreeMem` は psutil で代替して 1.0 を返す。並列を増やしても速くなるかは未確認 (`Fork` は §6.3。`Parallelize` は §6.2 の項)
- **Tk ウィジェットの値**: `w=Window[]`・`e=Entry[w,…]` の値は SAD と同じ `Widget[…]` (`Head[w]=Widget`) で、`e[Insert]`・`e[Get]`・`e[Delete]` (Entry・ListBox・Scale のメソッド) も働く ([20_tcltk_gui.md](20_tcltk_gui.md))。SAD と違う点は、最初の `Window[]` の親 (SAD は Toplevel、pysad は `.` 自身で、`TextVariable` の名前の番号が違う) と、`OpticsPlot` の text 項目の数 (フォントの測り方、[21_graphics.md](21_graphics.md))。パターン変数への代入 (`f[x_]:=(x=5); f[aa]`) は SAD と同じ。深い再帰では pysad は Python の RecursionError で落ちる (SAD は `General::deep`)。`Insert` の位置に `"End"` (大文字) を渡すと SAD は Tcl エラー 1501 ([20_tcltk_gui.md](20_tcltk_gui.md)・[14_beambeam.md](14_beambeam.md))
- `Check[Return[1], 2]` は SAD では Return を失敗として捕まえて 2 を返す。pysad は Return を `General::unexpbreak` で文ごと打ち切る (2 を返さない)。`Thread[a==b]` (a が値、b が未定義) は SAD が `3==b` を返すが、pysad は `==` が 0 を返すので Thread の引数にならない。リストと記号の ==・数どうしでない `==` の未評価形は一部だけ (`{1,2}==1` は未評価、記号どうしは 0)。
- `Out[n]:=` の表示は、最上位 (run_sad_string) の `;` なしの式文だけ。未定義の記号 1 つだけの文での SAD の表示 (前の Out の値を繰り返す) は未対応。`%%` は SAD でも記号のまま。
- SAD に無い関数 (`Mean` `Cross` `StringSplit` など) は未評価のまま返る (SAD と同じ。`Total` は pysad の拡張で動く)
- FFS 側: `MCOD/TCOD/MON` の COD 補正は K64 でも式として評価されて何もしない。`Steer`/`Monitor`/`CorrectOrbit` は `Packages/correction.n` を SADScript のまま実行する ([06_sadscript_library.md](06_sadscript_library.md) §5。`PYSAD_SAD_PACKAGES` が必要)。FFS コマンド・関数の対応状況と未対応の一覧は [03_ffs_commands.md](03_ffs_commands.md) §8
- **Packages の関数** (`Log2$`・`FunctionDefinition`・`Fill$String`・`Remove$Paren`・`Union` の `SameTest`・`Names`・`Emittance->False` ほか): [06_sadscript_library.md](06_sadscript_library.md) §5。除算 `a/b` は K64 が `a*(1/b)` で計算するので最後の桁が違うことがある (`Log[100]/Log[10]` は K64 1.9999999999999998、pysad 2)。評価器全体は変えていない (全試験の基準値に影響するため。要判断)。
- **複素引数の `PolyLog`・`HypergeometricU`・`Zeta`・`HurwitzZeta`** ([06_sadscript_library.md](06_sadscript_library.md) §1): [`pysad/specfn_cx.py`](../pysad/specfn_cx.py)。

### 6.2 文脈・`For`・`FindRoot`・`MessageName`・`Parallelize`・`VISIT` など (両モードの扱い)

#### 文脈 (`BeginPackage`・`Begin`・`End`・`EndPackage`)

文脈 (名前空間) は、K64 (init.n の定義) の規則を写した。**利用者のコードでは両モードで実装している**。SAD 自身の Packages (`.n`) の読み込みは文脈を無視して大域で読む。名前は「現在の文脈 (`$Context`) → `$ContextPath` の順」に探し、無ければ現在の文脈に作る。`BeginPackage[c1,c2,…]` は `$ContextPath` を `{c1,c2,…,System`}` にし、`Begin[c]` は `$Context` を c にする。`` `x `` は大域の x。文脈の中の記号は ``Stat`名`` という名前になり、外からは ``Stat`名`` で呼ぶ。`$Context`・`$ContextPath`・`Context[sym]` が読める。

- 使い方: [`example/08_scripting/syntax_package_context.sad`](../example/08_scripting/syntax_package_context.sad)。大域に公開する関数は `` `Name[…] := … `` と定義する。大域の変数を使うときは `` `var `` と書くか、`BeginPackage[My`, `]` のように `$ContextPath` に大域を入れる。組み込み関数・FFS の関数 (Twiss、LINE、FFS[…]、DP0 など) はそのまま使える。
- k64 モードと v2 モードの違いは `Names` と、`Begin` に文脈を**文字列**で渡す形 (k64 は K64 と同じく拒否、v2 は受ける。`Context["A`x"]` の文字列も同じ) だけ: ``Names["Stat`*"]`` は v2 がその文脈の記号の一覧、k64 が K64 と同じ ``{"Stat`"}``。
- K64 との意図的な差: (1) pysad の ``System` `` は組み込み関数 + pysad が記号として扱う名前 (`AccuracyGoal`・`Plot` などのオプション名を含む) で K64 より広い (K64 ではパッケージの中の `AccuracyGoal` は ``Stat`AccuracyGoal`` になり、`` `AccuracyGoal `` と書かないとオプションとして効かない。pysad はどちらでも効く)。(2) 表示: K64 は `$ContextPath` に無い文脈の記号に文脈名を付けて表示する (大域の記号が `` `sq `` になるなど)。pysad は文脈つきの名前だけ、`$ContextPath` にある文脈名を省く。(3) ``Symbol["Stat`x"]`` は K64 が値を返し、pysad は記号を返す。
- **`$Context`・`$ContextPath` への代入、`WithContext`、入れ子の文脈** (両モード): `$Context` への文脈の記号 (A の末尾にバッククォートを付けたもの) の代入と、文字列 (`"A`B`"` の形) の代入は init.n の `$SetContext` を呼び、読み取りの文脈が変わる (`End` の対象には積まない)。`$ContextPath` に文脈のリスト (A と System をバッククォートつきで並べたもの) を代入すると `$SetContextPath`。値が文脈でなければ `???General::wrongtype:  Argument must be Context (Symbol ending ...) in $SetContext[…]` (文脈はバッククォートで終わる記号) を出して元のまま。`Begin` は引数を `$SetContext` に渡すので、文字列で渡す `Begin` は拒否 (記号だけ。v2 は文字列も受ける)。`WithContext[{c__},body]` は `BeginPackage[c]; Begin[First[c]]; body; End[]; EndPackage[]` (body は文の読み取りの時点で名前が解決済み)。入れ子 ``A`B`` は一つの名前。`Context[x]` は文脈の記号を返し、文字列を渡す `Context` はエラー (k64)。文脈の記号は文脈が大域なので、大域が `$ContextPath` に無いとバッククォートつきの名前で表示する。**K64 は x・l・a など init.n が大域に作った名前を、文脈の中で読んでも大域の記号に当たる** (文脈 A の中で x=1 としたあとも `Context[x]` は大域): この 51 個の名前を pysad の大域の名前の初期値に入れた。使えない・未確認: 文脈つきの `Names` (K64 は状況で空のリストと文脈名のリストを返す。k64 は文脈名の一覧で規則が違う)、init.n が大域に作る名前の全部 (51 個しか調べていない)、System 文脈に入る名前 (pysad の System 文脈は K64 より広い)。文字列を `ToExpression` するときの文脈は「その時点の現在の文脈」。`BeginPackage` の戻り値 (`{$ContextPath,"`"}`) の K64 との一致は確かめていない (要確認)。SAD の Packages の読み込み (`load_package`、`PackagesDirectory` の下の `Get`) は文脈を扱わず、`PYSAD_SAD_PACKAGES` つきで GUI・描画・パッケージ系の試験が通る。
- 試験: `test_context_package_file`、`test_context_resolution_rules`、`test_context_package_using_ffs_and_builtins`、`test_chromaticity_scan_in_user_package` ([`test_v2_B_script.py`](../pysad/tests/script_tests/test_v2_B_script.py))。

#### `For`、`FindRoot` の複素、`MessageName` の関数形

- **`For[ini,cond,incr,body]`** (B71): K64 (gs) は 4 番目の引数 (本体) を実行せず、incr だけ繰り返す (3 引数は未評価)。`Packages/functions.n` の定義 ``For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)`` どおりなら実行するはずで、同じ式を別名で定義すると K64 も実行する (原因は未特定の K64 の読み込み時の癖)。**k64 モードは K64 と同じ** (本体を実行しない。SAD の Packages は `For` を使っていない)。**v2 モードは本体を実行する** (Break・Continue も効く。3 引数は本体なしのループ。t = Σ i² = 30)。`test_b71_for_body_v2`・`_k64`。
- **複素の初期値の `FindRoot`** (B61): v2 モードは複素の Newton 法で解く (前進差分のヤコビアンと最小二乗の Newton 法、歩幅を半分にする減速つき)。`FindRoot[z^2+1==0,{z,1.+1.*I}]` → `{z->(…+I),Residual->0}`。連立も可 (`{z*w==I,z+w==2}`)。`MaxIterations`・`AccuracyGoal` が効く。確かめた: z² = −1 → i、e^z = 2 → ln 2、z³ = 1、連立 zw = i・z + w = 2 (残差 1e-12 以下)。使えない: 変数の範囲 `{var, ini, {min, max}}` は複素では無視する、正則でない式 (Abs・Conjugate・Re・Im を含む) は収束を保証しない、実数の初期値で複素の値になる式は `General::wrongval`。k64 モードは TypeError。
- **`MessageName[f,"usage"]`** (B70): K64 (gs) は関数形で f を無視して `Undefined MessageName General::"usage"` を返す。v2 モードは f の本文を返す (`f::usage` と同じ。無ければ `General::tag` の本文。`General::mytag` の継承は K64 の `h::mytag` と同じ)。k64 モードの pysad も f を見る (`General::tag` への継承だけ)。`Message[f::usage]` の表示は K64 が `???f::usage:  本文` を出す。**pysad も出す** (両モード): `Message[f::tag, 値…]` は `???f::tag:  本文` (本文の `` `1` `` を値で置換)、本文が無ければ `General::tag`、それも無ければ `Undefined MessageName General::tag`。`Message[9, CA::open]` の先頭の数は重大度。引数の評価で再帰するもの (Tkinter の未定義メソッド) では出さず、`Tkinter::undef` は出さない。`Symbol["x"]` は x の値 (未定義なら記号 x)。`Off[f::tag]` による抑止は未対応 (pysad の `Off` は何もしない)。例 [`example/08_scripting/syntax_symbol_message.sad`](../example/08_scripting/syntax_symbol_message.sad)。
- **`Parallelize[f, n, lshare, npara]`** (両モード。[`pysad/fileops.py`](../pysad/fileops.py)。`Packages/shared.n` の `Fork`・共有メモリ・`Wait4` の動作に合わせた): POSIX (Linux・WSL・macOS) では npara > 1・n > 1 のとき `os.fork` で npara 個まで同時に f[k] を計算する (結果は子が一時ファイルに書き、親が回収)。結果の並びは SAD と同じ。子の代入は親に戻らない (親の変数は変わらない)。子の `Print` は出る (出力先が Python の StringIO のときは親に届かない: Fork と同じ)。lshare を超える結果は子が `Shared::toolarge` を出して `$Failed` (K64 は親が `Timeout` (2000 秒) まで待ち続ける)。Timeout 秒を超えた子は殺して `$Failed`。npara = 1・n = 1 は逐次。環境変数 `PYSAD_PARALLELIZE=serial` で逐次。**Windows は逐次** (`multiprocessing` の spawn は Interp・FFS・ラティスを pickle できず、子の側で作り直す方式は f が参照する大域の状態を再現できないので実装していない)。WSL (Python 3.14) で、FODO の CALC を 8 通りの K1 で逐次と並列の結果がビット単位で同じ。4 つの `Sleep[1]` が 1.02 秒。条件: **親でカーネル (numba) を一度走らせてコンパイルとキャッシュを済ませてから呼ぶ**。**numba と fork の併用は実際には試せていない** (WSL に numba が無く Windows に fork が無い。pysad のカーネルは `parallel=True`・`prange`・スレッドプールを使っていないので fork 後の子でそのまま動く見込みだが、(i) numba を `parallel=True` にする、(ii) 初回のコンパイルを子で行ってキャッシュ (`__pycache__`) を複数の子が同時に書く、(iii) OpenMP 版の BLAS を使う numpy、のどれかなら未確認)。例 [`example/09_io_system/fn_Parallelize.sad`](../example/09_io_system/fn_Parallelize.sad)。試験 [`test_v2_B2_script.py`](../pysad/tests/script_tests/test_v2_B2_script.py)。
- **`VISIT`・`BYE`・`USE` と `ElementValues` ほか** (両モード): `VISIT` は FitFunction・FitValue・ElementValues・OpticsProlog・OpticsEpilog の 5 つの記号を隠し (訪問先では未定義)、`BYE` で戻る (入れ子も)、`USE` は 5 つとも消す (`tffsa.f` の `tfblocksym`・`tfunblocksym`)。`BYE` だけで VISIT が無いと `???-FFS-Error-BYE without VISIT.`。`BYE ALL` は K64 が SIGSEGV で動かない (pysad は一番外側の値に戻す)。台本 `FFS_VISIT_BYE.sad`。
- **`Ffs;`・`ffs use ring;`**: `FFS` と `ffs` だけが FFS の語で、`Ffs;`・`FFs;` は MAIN の構文エラー (SAD と同じ)。LINE 名は書かれたとおり (`USE ring` は ring が無いとエラー)。[03_ffs_commands.md](03_ffs_commands.md) の `FFS;` の項。

### 6.3 `Fork`・`Shared`・`Wait`

- **共有メモリ** (両モード): `s = OpenShared[n]` (n バイト)、`Shared[s] = 値`、`Shared[s]`、`Shared[s]++`、`Shared[s] += v`、`Close[s]`。値は数・文字列・記号・リスト・式・規則 (無名の mmap に符号化して書く)。`SharedSize[値]` が n を超えると `Shared::toolarge` (K64 と同じ式: リストは (1 + 要素数) × 8 + 中身。`SharedSize[{1,2,3}]` = 32)。
- **`Fork[]`**: Linux・WSL・macOS では `os.fork` で本当に分ける (親には子の pid、子には 0。子は `Exit[]` か最上位の文の終わりで `os._exit`)。`Wait[]` は `{pid, status}`、`Wait4[pid, flags…]` は `{pid, 終了コード, Exited->…, …}`、`Kill[pid]` はこの pysad が作った子だけ (安全のための制限)、`Sleep[t]` は本当に待って 1 を返す、`Pause[t]` は子がいる間だけ待つ。**Windows には fork が無い**: v2 は `System::error: Fork[] is not available on this platform …` を出して文を打ち切る (k64 モードは pid を返し子は作らない)。並列の計算は `Parallelize[f, n, lshare, npara]` (§6.2)。
- 確かめた範囲: WSL (Ubuntu、Python 3.14、numba は代用の空モジュール) で例 ([`example/09_io_system/fn_SharedMemory_Parallelize_System_Fork.sad`](../example/09_io_system/fn_SharedMemory_Parallelize_System_Fork.sad)) と 3 個の子の試験。ラティスを使う重い計算を子で走らせる使い方 (numba のスレッドと fork の併用) は未確認。子の `Print` は親と同じ出力先に出る (出力先が Python の StringIO のときは親に届かない)。

## 7. 検証・試験

- 試験は [`pysad/tests/`](../pysad/tests/) (pytest。収集 941 件、2026-11-17 に 937 通過・4 skip。この件数は k64 モードを既定にした実行のもので、現在の既定の v2 での件数は未再確認) と `pysad/tests/regression_tests/run_regression.py` (`QUICK=1` で 50 ベンチ FAIL 0。k64 モード)。
- K64 (SAD の gs) との比較では、FIT/FREE/GO の台本 (`bench_level13a`・`13b`) の出力、色収差スキャンの台本 (`bench_level14a`、チューン 21 点が 2e-8 以内)、Element/BeamLine 系 (`bench_level24f`、164 タグ一致) が一致する。詳細な突き合わせは [27_benchmarks.md](27_benchmarks.md)。
- 未確認・限界は §6 の各項と [31_unimplemented_and_unverified.md](31_unimplemented_and_unverified.md)。

関連文書: [03_ffs_commands.md](03_ffs_commands.md) (FFS の語と関数)、[06_sadscript_library.md](06_sadscript_library.md) (特殊関数・Fit・Packages)、[12_pysad_api.md](12_pysad_api.md) (Python API)、[28_v2_vs_k64.md](28_v2_vs_k64.md) (v2 と k64 の違い)、[29_sad_quirks.md](29_sad_quirks.md) (SAD の癖)、[30_sad_bug_report.md](30_sad_bug_report.md) (B 番号の目録)。
