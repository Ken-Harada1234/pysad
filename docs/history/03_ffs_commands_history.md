# 03_ffs_commands.md から取り除いた履歴の記述

03_ffs_commands.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 03. FFS のコマンド解釈 (src/tffsa.f) と pysad の ffscmd.py

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)
````

## ### 3.3 FREE / FIX (tffsfreefix.f)

整理前の 106 行目:

````text
`FREE pat [KEY] pat …`: `tfgetlineps(mode 0)` で要素位置を集め、`.` を含まない名前は **ファミリー先頭 (klp==k) だけ** を変数にする。`KEY` (次の語) が要素型のキーワード (`tfkwrd`) なら ivvar にそれを、無ければ既定キーワード `nelvx%ival` (tfinit.f: DRIFT は L (ival=1, 下限 vlim(1)=0.1)、BEND/QUAD/SEXT/OCT/DECA/DODECA/SOL/CAVI/TCAV は 2 番目のキーワード (ANGLE/K1/K2/K3/…/BZ/VOLT)、MULT は K1、MARK/MONI は 0 → `Can't use as variable`)。`BXI` 等の `nlistI` は先頭 MARK の初期値を変数にする (it=41)。`QF.2` のように成分を指定すると `ivcomp` 付き (`Info-Component … is coupled to`)。`valvar = 値/errk` で保持。**成分ごとの FREE** (v2.1.0。`tffsfreefix.f`・`tfsetv.f`・`tffsa.f` の `tfvcomp` と gs の実測。pysad は以前 `FREE QF.3 K1` をファミリー全体の変数にしていた不具合を直した。両モード): 語に `.` があれば成分の変数。(a) 既定のキー (QUAD の K1 など) の変数は要素全体の値で、成分は先頭の成分に結合する。先頭の成分 (`QF.1`) は要素全体と同じ (結合先を `???-FFS-Info-Component QF.2 is coupled to QF.1 .` で知らせる)。**先頭でない成分 (`QF.3`) の既定のキーの変数は値を持たず効かない** (K1 は動かない。gs は変数の数に数え、pysad は作らない: 反復の行の最後の数が 1 つ違い、その変数だけで GO すると pysad は `?No variable.`)。(b) 既定でないキー (BEND の K1、ROTATE) の変数は、その成分だけの値 (`tfsetv` の `ievar2`)。(c) 同じ要素の変数は、後から FREE したものが先 (`LOOP_I_2`)。(d) 要素全体の変数がある所へ成分、成分がある所へ要素全体を FREE するとエラー (`A component has been already used as variable: B.6`)。(e) `tfvcomp` は (b) の変数があれば真で、NOKEEPEXP の CALC・GO は展開しない (NOCALEXP)。(f) FIX: `FIX B.3` は成分、`FIX B` は成分に触れない、`FIX B*`・`FIX *` は成分も消す。`VAR` の Couple の欄は既定のキーだけ `<--`、制限値は既定のキー ±1e10、ほかは ±1e99。小さな FODO (3 GeV、周長 18.4 m) で B.3 だけの K1 で NX を 1.30 にすると反復 3 回 (残差 2.1729E-03・4.8421E-08・3.2593E-17)、NX 8.168140905042474 が gs と全桁。大きなリング A (3 GeV、周長 603 m、3563 要素。`FREE BPI6.2 K1; FREE BPI7.5 K1` と `FREE QDA2.1 K1; FREE QDA3.4 K1; FREE QDA4`) でも反復の行・チューン (47.6000342385、17.6000441754)・K1 が gs と一致 (2.8 秒)。未確認: COUPLE で結合した要素の成分、`USE NOEXPAND` の展開済みの成分の値の保存。例 [`example/v2_next/B2_free_component.sad`](../example/v2_next/B2_free_component.sad)。`FIX pat [KEY]` は一致する変数を削除。ワイルドカード (`%*<{|`) があれば複数語を続けて読む。
````

## ## 5. USE / CELL / CALC / GO / EMIT / TYPE / DISP / SAVE / RESET / COUPLE / CHRO

整理前の 144 行目:

````text
- **EMIT のあとの `USE`** (v2.1.0。両モード。K64 の仕様): `EMIT` (または `Emittance[]`・`CALC6D`) のあとに `USE RING;` しても、追跡の RF の位相の基準 `trf0` は EMIT の値のまま (`trf0` は `tmacro` の大域変数。`temit.f` が 0 にしてから `tturne` で計算し、`tffsa.f` の MEASURE (`trf0 = 0`) だけが他に触る。`USE` は触らない)。`CELL; CALC;` のあとも変わらない。以前の pysad は `USE` で `trf0 = 0.0` にしていた。gs と全桁一致 (小リングで EMIT 前 δ = (−3.0771610e-4、6.82902209e-4、−1.023444e-2)、EMIT 後・USE 後 (−4.1421124e-5、9.50018237e-4、−9.97733e-3))。大きなリング A (CELL CALC・RAD・EMIT の 5 周追跡 → USE RING → 5 周追跡で 1.3e-11)・PF リングでも一致。`emit`・Touschek の結果 (Twiss など) は USE で作り直される (捨てる。K64 での残り方は未確認)。例 [`example/v2_next/A1_use_after_emit.sad`](../example/v2_next/A1_use_after_emit.sad)。
````

## ## 8. pysad ffscmd.py の対応表

整理前の 226 行目:

````text
- `NOKEEPEXP` の tfvcomp: 既定でないキーの成分の変数があれば真 (成分ごとの FREE を実装した v2.1.0 から。K64 の定義と同じ)。NOKEEPEXP と RENUMBER の INDEX 属性による番号の初期値は実装した (両モード。gs の実測と一致。下記)。
````

## ### v2 モードの MEASURE・TRCOD と、台本の読み方・小さな修正

整理前の 248 行目:

````text
**v2 モードの MEASURE・TRCOD・LOGMEASURE** ([`pysad/measure.py`](../pysad/measure.py)。k64 モードは従来どおり `FFS command MEASURE は未対応 (無視)`)。K64 の実行ファイル (gs) は MEASURE が無出力、TRCOD が SIGSEGV (B42)、MAIN の TRACK も落ちるので、**gs との数値比較は 1 件もない**。実装は K64 のソース (`tffsa.f`、`trackb.f`、`tracke.f`、`trcoda.f`、`ttinit.f`、`ttstat.f`、`tbdecoup.f`) の読みで、値は解析式と pysad の独立な経路 (CALC の Twiss、EMIT、TrackParticles) で確かめた。表示の書式は gfortran で同じ FORMAT を出して合わせた (gs の実出力との比較ではない)。
````

整理前の 252 行目:

````text
- `TRCOD FIXCOD [位置]` / `MVCOD [位置]` / `STNCOD`: 軌道・分散の「基準」を取っておき、`Twiss["RDX", k]` で基準、`Twiss["DDX", k]` で現在 − 基準 (EX EPX EY EPY DX DPX DY DPY)。基準は REFERENCE と同じ参照光学に置く (K64 の `twiss(*,ndim,*)` を利用者が読む手段が無いため。`REF;` で上書きされるので戻すには `TRCOD STNCOD;`)。`TRCOD ALIGN 要素1 要素2`: 要素1 の DX・DY を、その次から 要素2 までの BEND QUAD SEXT OCT DECA DODECA MULT の DX・DY に書く (DECA DODECA MULT は v2 で足した。回転・縦のずれは揃えない)。`ORBITJITTER xj yj`: Twiss の DX・DY に位置ごとの正規乱数を足す。
````

整理前の 268 行目:

````text
| `FFS;` の書き方 | (gs の動作) | gs に合わせた (モードによらない)。`FFS USE RING;`・`FFS USE=RING;`・行の途中の `FFS;`・`FFS NPARA=3;` (gs と同じエラー)。**`FFS` と `ffs` だけが FFS の語**で、`Ffs;`・`FFs;` は MAIN の構文エラー (`LALR(1) parser: syntax error`) で FFS 部に入れず何も実行されない (v2.1.0。以前の pysad は読んでいた)。**LINE 名は書かれたとおり** (LINE ring と USE RING は別。`USE ring` は ring が無いので `???General::wrongtype:  Argument must be BeamLine[ ... ]:` と `???-FFS-Error-Missing beamline in USE.`) | 台本 15 本 + 6 本を gs と照合。同じ行の `FFS;NPARA=3;` の `NPARA=3;` は FFS 部では実行されない (gs と同じ) |
````

整理前の 269 行目:

````text
| MAIN の一般の代入 `NAME = 式;` と要素定義の式 (`ANGLE=PI/8`) | 無視されていた (要素定義の式は `PI`・関数を文字列として扱い警告して要素を捨てた) | gs と同じ式の文法 (v2.1.0 で約 100 通りを gs で実測): 単位つき数 (`DEG MM CM …` は**数の直後だけ**。`10 DEG/2` は (10 DEG)/2)、`+ - * / ^` (`^` は右結合)、関数 `SIN COS TAN ATAN SQRT EXP LOG (常用) LN SINH COSH TANH ATAN2(y,x)` (関数名は大文字小文字を問わない、丸括弧)、`PI` (大文字だけ)、`MOMENTUM`、前の変数。0 での割り算は ±Infinity・NaN、定義域の外は NaN。gs で構文エラーになる形 (`pi/8`・`Pi/8`・`(1+2) MM`・`SQRT 2`・`ABS`・`2**3`・`2^-1` ・引数の数が違う関数など) は警告 `syntax error in MAIN expression (gs: …)` を出して捨てる | 代入 32 個の式を gs の値と 2e-14、要素定義の式 37 通りの値と 19 通りの構文エラーを gs と全桁 (モードによらない)。大きなリング A・B、PF、ASTRID2 の全要素の値は v2.0.0 の読み取りとビット一致。例 [`example/v2_next/D2_main_expr.sad`](../example/v2_next/D2_main_expr.sad) |
````

整理前の 271 行目:

````text
`FFS;` と MAIN の代入の残した違い: `FFS` (`;` なし) の次の行が USE 以外の文のとき gs はエラーにして行の残りを捨てるが pysad は FFS 部として実行する、`A = 3 B = 4;` は gs では A = 3 (エラーつき) だが pysad は代入しない、gs が構文エラーにするが pysad が受ける形 (`45DEG`・`SIN[0.1]`・未定義の名前 `L=A`)、gs が受けるが pysad が構文エラーにする形 (末尾の `1+`)、MAIN の変数を `E`・`Pi`・`Degree` で代入しても要素定義の式の `E` は SADScript の定数 (gs は利用者の値)、`LINE[{"EMITY",…},1]` は gs が MAIN の値、pysad が MARK に書いた値 (0)。k64 モードでも MAIN で PBUNCH・EMITX・MASS などを代入している TRACK なしの台本は、これまで無視されていた値が効くようになる (gs と同じ向き)。
````

## ## 9. この SAD ビルド (WSL, K64) で分かった注意

整理前の 275 行目:

````text
- **NPARA の既定は 4** (Packages/init.n:219 `NPARA=4;`。initb1.f の `RsetGL('NPARA',1)` を init.n が上書き)。TrackParticles は粒子 9 個以上 (有効な並列数 = min(NPARA, np/9 + 1, ne·np/3000 + 1)。np ≥ 13 で実際に分かれる) で fork するため乱数列が分岐する。pysad の既定は NPARA=1。gs と比較するには `NPARA=1;` をスクリプトに入れるか、**`NPARA=n;` を FFS で設定して同じ分け方を再現する** (v2.1.0。子の乱数の種 = 親の種 + 2k。RAD + FLUC の追跡だけ。gs の NPARA=1 と 4 の rms dp・最終の種と全桁一致: 08 の 1.2 節)。
````

整理前の 276 行目:

````text
- **FITP>1 の GO で SAD が再起動することがある** (`*** Run time Environment` バナーが再び出る)。NPARA=1 でも起きるので、オフモメンタム条件の参照値は得られていない。v2.1.0 の確認 (C): [`example/03_matching/FFS_FITPOINTS.sad`](../example/03_matching/FFS_FITPOINTS.sad) (24 セル FODO、δ = −1%・0・+1% の 3 点で NX = 6.7・NY = 4.9) は今の gs は再起動せず完走し、**最小値の手前で止まった** (K1 (QF, QD) = (0.700211, −0.599628)、K2 (SF, SD) = (1.695226, −2.439058)、残差 6.3321e-6 が独立な最小値 6.33178e-6 より相対 5e-5 大きい。`Overflow & closed orbit not found` の試行点での打ち切り)。pysad の GO の最終解は独立な最小化 (残差の定義 Σ w (f − goal)² scale² の scipy の最小二乗) と 5e-8 以内で一致。
````

整理前の 283 行目:

````text
- **`Skip`・`Read` の型** (Character・Real・Word の読み方、数でない語の Real は式の値) は gs に合わせた ([`test_evalfixes.py`](../pysad/tests/test_evalfixes.py))。`MessageName[f,"usage"]` の関数形は、gs は f を無視して `Undefined MessageName General::"usage"` を出す (B70)。v2 モードは f の本文を返す (`f::usage` と同じ。無ければ `General::tag` の本文)。k64 モードの pysad も以前から f を見る (`General::tag` への継承を足した)。
````

