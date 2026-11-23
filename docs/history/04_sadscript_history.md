# 04_sadscript.md から取り除いた履歴の記述

04_sadscript.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 04. SADScript 言語と pysad の sadscript.py

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)
````

## ### 5.4 11/13 以降に足した関数群 (詳細は各文書)

整理前の 188 行目:

````text
### 5.4 11/13 以降に足した関数群 (詳細は各文書)
````

## ## 6. 対応していない構文・関数

整理前の 212 行目:

````text
- **`BeginPackage`/`Begin` の文脈 (名前空間)**: 利用者のコードでは**両モードで実装した** (下の「文脈」)。SAD 自身の Packages (`.n`) の読み込みは従来どおり文脈を無視して大域で読む。
````

## #### 文脈 (`BeginPackage`・`Begin`・`End`・`EndPackage`)

整理前の 219 行目:

````text
- k64 モードと v2 モードの違いは `Names` と、`Begin` に文脈を**文字列**で渡す形 (k64 は gs と同じく拒否、v2 は受ける。`Context["A`x"]` の文字列も同じ) だけ: ``Names["Stat`*"]`` は v2 がその文脈の記号の一覧、k64 が gs と同じ ``{"Stat`"}``。k64 モードで以前と変わった点: 以前は利用者のコードの `BeginPackage`・`Begin` を無視して全部大域だった。
````

整理前の 221 行目:

````text
- **`$Context`・`$ContextPath` への代入、`WithContext`、入れ子の文脈** (v2.1.0。両モード。gs で実測し、台本 [`example/v2_next/B2_context.sad`](../example/v2_next/B2_context.sad) の 19 行が gs と全行一致): `$Context` への文脈の記号 (A の末尾にバッククォートを付けたもの) の代入と、文字列 (`"A`B`"` の形) の代入は init.n の `$SetContext` を呼び、読み取りの文脈が変わる (`End` の対象には積まない)。`$ContextPath` に文脈のリスト (A と System をバッククォートつきで並べたもの) を代入すると `$SetContextPath`。値が文脈でなければ `???General::wrongtype:  Argument must be Context (Symbol ending ...) in $SetContext[…]` (文脈はバッククォートで終わる記号) を出して元のまま。`Begin` は引数を `$SetContext` に渡すので、文字列で渡す `Begin` は拒否 (記号だけ。v2 は文字列も受ける)。`WithContext[{c__},body]` は `BeginPackage[c]; Begin[First[c]]; body; End[]; EndPackage[]` (body は文の読み取りの時点で名前が解決済み)。入れ子 ``A`B`` は一つの名前。`Context[x]` は文脈の記号を返し、文字列を渡す `Context` はエラー (k64)。文脈の記号は文脈が大域なので、大域が `$ContextPath` に無いとバッククォートつきの名前で表示する。**gs は x・l・a など init.n が大域に作った名前を、文脈の中で読んでも大域の記号に当たる** (文脈 A の中で x=1 としたあとも `Context[x]` は大域): この 51 個の名前 (gs の `Context[名前]` の実測) を pysad の大域の名前の初期値に入れた。使えない・未確認: 文脈つきの `Names` (gs は状況で空のリストと文脈名のリストを返す。k64 は文脈名の一覧で規則が違う)、init.n が大域に作る名前の全部 (gs の `Names` が列挙しないので 51 個しか調べていない)、System 文脈に入る名前 (pysad の System 文脈は gs より広い)。文字列を `ToExpression` するときの文脈は「その時点の現在の文脈」。`BeginPackage` の戻り値 (`{$ContextPath,"`"}`) の gs との一致は v2 の記録では確かめていない (要確認)。SAD の Packages の読み込み (`load_package`、`PackagesDirectory` の下の `Get`) は変えていない (`PYSAD_SAD_PACKAGES` つきで GUI・描画・パッケージ系の試験が通る)。
````

## #### `For`、`FindRoot` の複素、`MessageName` の関数形

整理前の 226 行目:

````text
- **`For[ini,cond,incr,body]`** (B71): gs は 4 番目の引数 (本体) を実行せず、incr だけ繰り返す (3 引数は未評価)。`Packages/functions.n` の定義 ``For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)`` どおりなら実行するはずで、同じ式を別名で定義すると gs も実行する (原因は未特定の K64 の読み込み時の癖)。**k64 モードは gs と同じ**にした (本体を実行しない。以前の pysad は本体を実行していた。SAD の Packages は `For` を使っていない)。**v2 モードは本体を実行する** (Break・Continue も効く。3 引数は本体なしのループ。t = Σ i² = 30)。`test_b71_for_body_v2`・`_k64`。
````

整理前の 227 行目:

````text
- **複素の初期値の `FindRoot`** (B61): v2 モードは複素の Newton 法で解く (前進差分のヤコビアンと最小二乗の Newton 法、歩幅を半分にする減速つき)。`FindRoot[z^2+1==0,{z,1.+1.*I}]` → `{z->(…+I),Residual->0}`。連立も可 (`{z*w==I,z+w==2}`)。`MaxIterations`・`AccuracyGoal` が効く。確かめた: z² = −1 → i、e^z = 2 → ln 2、z³ = 1、連立 zw = i・z + w = 2 (残差 1e-12 以下)。使えない: 変数の範囲 `{var, ini, {min, max}}` は複素では無視する、正則でない式 (Abs・Conjugate・Re・Im を含む) は収束を保証しない、実数の初期値で複素の値になる式は従来どおり `General::wrongval`。k64 モードは従来どおり TypeError。
````

整理前の 228 行目:

````text
- **`MessageName[f,"usage"]`** (B70): K64 (gs) は関数形で f を無視して `Undefined MessageName General::"usage"` を返す。v2 モードは f の本文を返す (`f::usage` と同じ。無ければ `General::tag` の本文。`General::mytag` の継承は gs の `h::mytag` の実測と同じ)。k64 モードの pysad も以前から f を見る (`General::tag` への継承だけ足した)。`Message[f::usage]` の表示は gs が `???f::usage:  本文` を出す。**v2.1.0 で pysad も出す** (両モード): `Message[f::tag, 値…]` は `???f::tag:  本文` (本文の `` `1` `` を値で置換)、本文が無ければ `General::tag`、それも無ければ `Undefined MessageName General::tag`。`Message[9, CA::open]` の先頭の数は重大度。引数の評価で再帰するもの (Tkinter の未定義メソッド) では出さず、`Tkinter::undef` は出さない。`Symbol["x"]` は x の値 (未定義なら記号 x)。`Off[f::tag]` による抑止は未対応 (pysad の `Off` は何もしない)。例 [`example/v2_next/B2_symbol_message.sad`](../example/v2_next/B2_symbol_message.sad) (12 行が gs と全行一致)。
````

整理前の 229 行目:

````text
- **`Parallelize[f, n, lshare, npara]`** (v2.1.0。両モード。[`pysad/fileops.py`](../pysad/fileops.py)。`Packages/shared.n` の `Fork`・共有メモリ・`Wait4` の動作を gs で確認): POSIX (Linux・WSL・macOS) では npara > 1・n > 1 のとき `os.fork` で npara 個まで同時に f[k] を計算する (結果は子が一時ファイルに書き、親が回収)。結果の並びは gs と同じ。子の代入は親に戻らない (親の変数は変わらない)。子の `Print` は出る (出力先が Python の StringIO のときは親に届かない: Fork と同じ)。lshare を超える結果は子が `Shared::toolarge` を出して `$Failed` (gs は親が `Timeout` (2000 秒) まで待ち続ける)。Timeout 秒を超えた子は殺して `$Failed`。npara = 1・n = 1 は逐次。環境変数 `PYSAD_PARALLELIZE=serial` で逐次。**Windows は逐次** (`multiprocessing` の spawn は Interp・FFS・ラティスを pickle できず、子の側で作り直す方式は f が参照する大域の状態を再現できないので実装していない)。WSL (Python 3.14) で、A (5 点)・B (リストの結果)・C (1/0 → Infinity)・E (x は 0 のまま)・F (子の Print) が gs と同じ。FODO の CALC を 8 通りの K1 で逐次と並列の結果がビット単位で同じ。4 つの `Sleep[1]` が 1.02 秒。条件: **親でカーネル (numba) を一度走らせてコンパイルとキャッシュを済ませてから呼ぶ**。**numba と fork の併用は実際には試せていない** (WSL に numba が無く Windows に fork が無い。pysad のカーネルは `parallel=True`・`prange`・スレッドプールを使っていないので fork 後の子でそのまま動く見込みだが、(i) numba を `parallel=True` にする、(ii) 初回のコンパイルを子で行ってキャッシュ (`__pycache__`) を複数の子が同時に書く、(iii) OpenMP 版の BLAS を使う numpy、のどれかなら未確認)。例 [`example/v2_next/B2_parallelize.sad`](../example/v2_next/B2_parallelize.sad)。試験 [`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py)。
````

整理前の 230 行目:

````text
- **`VISIT`・`BYE`・`USE` と `ElementValues` ほか** (v2.1.0。両モード): `VISIT` は FitFunction・FitValue・ElementValues・OpticsProlog・OpticsEpilog の 5 つの記号を隠し (訪問先では未定義)、`BYE` で戻る (入れ子も)、`USE` は 5 つとも消す (`tffsa.f` の `tfblocksym`・`tfunblocksym`)。`BYE` だけで VISIT が無いと `???-FFS-Error-BYE without VISIT.`。`BYE ALL` は gs が SIGSEGV で比べられない (pysad は一番外側の値に戻す)。台本 `B2_visit_bye.sad` の 8 行が gs と全行一致。
````

整理前の 231 行目:

````text
- **`Ffs;`・`ffs use ring;`** (v2.1.0): `FFS` と `ffs` だけが FFS の語で、`Ffs;`・`FFs;` は MAIN の構文エラー (gs と同じ。以前の pysad は読んでいた)。LINE 名は書かれたとおり (`USE ring` は ring が無いとエラー)。docs/03 の `FFS;` の項。
````

整理前の 232 行目:

````text
- **Packages の関数** (v2.1.0。`Log2$`・`FunctionDefinition`・`Fill$String`・`Remove$Paren`・`Union` の `SameTest`・`Names`・`Emittance->False` ほか): 小さな台本を gs と行ごとに突き合わせた (docs/34)。`Element[..]+=`・`LINE[..]+=` の複合代入が `Element`・`LINE` のフックを通らず利用者の規則を作っていた pysad の不具合を直した。除算 `a/b` は gs が `a*(1/b)` で最後の桁が違うことがある (`Log[100]/Log[10]` は gs 1.9999999999999998、pysad 2)。評価器全体は変えていない (全試験の基準値に影響するため。要判断)。
````

整理前の 233 行目:

````text
- **複素引数の `PolyLog`・`HypergeometricU`・`Zeta`・`HurwitzZeta`** (v2.1.0。docs/33): [`pysad/specfn_cx.py`](../pysad/specfn_cx.py)。
````

## #### `Fork`・`Shared`・`Wait`

整理前の 238 行目:

````text
- **`Fork[]`**: Linux・WSL・macOS では `os.fork` で本当に分ける (親には子の pid、子には 0。子は `Exit[]` か最上位の文の終わりで `os._exit`)。`Wait[]` は `{pid, status}`、`Wait4[pid, flags…]` は `{pid, 終了コード, Exited->…, …}`、`Kill[pid]` はこの pysad が作った子だけ (安全のための制限)、`Sleep[t]` は本当に待って 1 を返す、`Pause[t]` は子がいる間だけ待つ。**Windows には fork が無い**: v2 は `System::error: Fork[] is not available on this platform …` を出して文を打ち切る (k64 モードは従来どおり pid を返し子は作らない)。並列の計算は `Parallelize[f, n, lshare, npara]` (下の項)。
````

整理前の 239 行目:

````text
- 確かめた範囲: WSL (Ubuntu、Python 3.14、numba は代用の空モジュール) で例 ([`example/v2_B/B3_shared_system.sad`](../example/v2_B/B3_shared_system.sad)) と 3 個の子の試験。ラティスを使う重い計算を子で走らせる使い方 (numba のスレッドと fork の併用) は未確認。子の `Print` は親と同じ出力先に出る (出力先が Python の StringIO のときは親に届かない)。`Parallelize` は v2.1.0 で本当の並列にした (下の項)。
````

