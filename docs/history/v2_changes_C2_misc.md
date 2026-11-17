# v2 の変更 (担当 C2 の一部): 旧式 WAKE の表・DELCOR の規格化・台本の読み方 2 件 (2026-11-17)

gs は K64 の 2026-10-19 ビルド (`/mnt/d/work/SAD-master/SAD-master/bin/gs`)。試験は `pysad/tests/test_v2_C2_wakecmd.py` (6)、
`test_v2_C2_delcor.py` (7)、`test_v2_C2_script.py` (18)。例は `example/v2_C/C2w_*.sad` (出力 `*.pysad.txt` は v2 モード)。
SADHelp.html に DELCOR・ENDWAKE の記述は無い (`WAKE` の旧式の構文も無い)。根拠は K64 のソースと gs の実動作。

## 1. 旧式の `WAKE` コマンドの表を行をまたいで書く (v2 の機能。k64 モードは従来どおり)

| | |
|---|---|
| K64 (gs) | `WAKE` は `WAKECOMMAND[]` を評価するだけ (tffsa.f:976。配布に WAKECOMMAND.n が無く失敗)。tfwake.f は呼ばれないが、その読み方は「`ENDWAKE` まで、行末なら次の行を読む (skiplnget)」「数でない語は `?Missing ENDWAKE.`」「w が無ければ `?Missing value of wake.`」。 |
| v2 (pysad) | `WAKE L…\|T… 要素` で始まる文は、`ENDWAKE` (または `;`) が来るまで次の行に続く。1 行に何組でも、注釈 `!` 可、`ENDWAKE` の後の `;` は任意、同じ行の続き (`ENDWAKE; CALC;`) も実行する。大文字小文字・省略形 (`wake t RF`) 可。要素名はワイルドカード可 (合う CAVI すべてに登録)。 |
| 根拠 | 同じ表を `WakeFunction[…] = …` で与えた追跡と同じ結果: 6 点・5 点の表で gs の参照値 (test_wake の C3_REF) と 2e-9、400 点の表・2 台の CAVI (RF*) で WakeFunction と 1e-12。 |
| 試験 | `test_multiline_table_equals_wakefunction`、`test_table_of_hundreds_of_points`、`test_wildcard_registers_all_cavities`、`test_errors_are_explained`、`test_statement_after_endwake_on_same_line`、`test_k64_mode_reads_line_by_line` |
| 例 | `example/v2_C/C2w_wake_table.sad` |

書き方:
```
WAKE LONGITUDINAL RF
  -0.1   2e9          ! s [m]  w [V/C] (横 TRANSVERSE は V/C/m)
   0     5e8    0.02 -1e9
ENDWAKE;
```
- 使える範囲と条件: v2 モード。CAVI だけ。追跡は WakeFunction と同じ (TRPT と LWAKE / TWAKE が ON の TrackParticles。`docs/38_wake_guide.md`)。点数の上限は無い (400 点で確認。K64 の tfwake.f は約 8000 組まで)。値は数か式 (`1e9*2`)。
- 使えない範囲: 表の途中の `;` (tfwake.f と同じくエラー)。CAVI 以外の要素。k64 モード (WAKE は WAKECOMMAND[] のまま、文は行ごと)。`ENDWAKE` を書き忘れると、次の `;` までを表として読んでエラーにする (その間の文は実行されない)。
- エラーの文: `?Missing ENDWAKE. (WAKE RF: … N 個の値の後の "語" は数として読めない …)`、`?Missing value of wake. (… 値が N 個 (奇数) …)`、`?WAKE: QF は CAVI ではない (Quadrupole)…`、`?Undefined element 名前 …`、`?No wake data. …`。エラーの後も台本は続く。
- 評価器の変更: `sadscript.py` の `Interp.split_statements` が、改行で文を切る前に `_wake_table_open(stmt)` (新設。v2 で、文が `WAKE L|T …` で始まり `ENDWAKE` をまだ含まないとき真) を見る。先頭 4 文字が WAKE でない文は何も変わらない。`ffscmd.py` の `_wake_cmd_v2` はエラーの文と空の表の検査だけ変更。

## 2. DELCOR の corlength≠0 の規格化 (v2 で修正。k64 モードは K64 のまま)

| | |
|---|---|
| K64 (gs) | tvcorr.f: `s = 1/sqrt(4Σw²+2)`。複素数 c の絶対値の 2 乗平均が 1。tcorr.f は DX = Re(c)·delx、DY = Im(c)·dely と使うので、DX の rms は delx/√2。 |
| v2 (pysad) | `s = 1/sqrt(2Σw²+1)` (実部・虚部それぞれの分散が 1)。DX の rms = delx、DY の rms = dely。 |
| 試験 | `test_rms_of_each_component_is_one`、`test_long_correlation_limit_is_a_common_gaussian_shift`、`test_correlation_function_shape_is_kept`、`test_k64_mode_keeps_k64_normalization`、`test_random_walk_branch_is_per_component_and_same_in_both_modes`、`test_delcor_command_rms_of_dx_is_delx`、`test_delcor_command_k64_rms_is_delx_over_sqrt2` |
| 例 | `example/v2_C/C2w_delcor_rms.sad` |

調べたこと (作者の意図):
- SADHelp.html に DELCOR の項は無い。K64 の配布 (src・Packages・例) で DELCOR を使う・説明する箇所は tffsa.f の呼び出しと tcorr.f の `Syntax: DELCOR corlength delx dely element…` だけ。delx が何の量かを書いた文は無い。
- 式の数え方: 波数 ±k の各成分は `dcmplx(w·g, w·g)` で E|·|² = 2w²、波数 0 は 2。合計 4Σw²+2 は E|c|² そのもの。つまり K64 は「複素数として rms 1」に規格化している。一方、c の実部と虚部は DX 用・DY 用の独立な 2 つの実数の場を入れる器で、|c| はどこにも使われない (delx と dely は別々に指定する)。
- corlength=0 の枝 (同じ tvcorr.f): 1 歩 `dcmplx(w·g, w·g)`、w = sqrt(al/n)。実部だけで 1 歩の分散が al/n、つまり 1 m あたり 1。成分ごとの規格化で、1/2 は入っていない (試験: 中点の実部の分散 = al/4)。同じサブルーチンの中で corlength≠0 だけが 2 成分をまとめて規格化していて、整合しない。
- 極限: 相関長 ≫ 周長では波数 0 の成分だけが残り、全要素が同じ値 delx·Re(g₁+i g₂)·s だけずれる。v2 は s = 1 でずれは delx·g (rms delx)、K64 は delx·g/√2。「全体が delx の rms で一緒に動く」が自然。
- 以上から、「delx は DX の rms」と読むのが自然で、4Σw²+2 は 2 成分を数えた取り違えと判断して直した。ただし作者の文書による裏づけは無い (推定)。K64 の値が必要なら k64 モード、または delx に 1/√2 を掛ける。

根拠の数値 (GCUT なし):
- 場の rms (500 m、512 点、400 標本): x0 = 3、15、5000 m で実部・虚部とも 1 ± (2、4、8)% 以内 (x0 = 15: 0.987、1.009。k64: 0.708、0.703)。
- DELCOR コマンド (12 セルのリング 62.4 m、四極 24 台、corlength 5 m、150 回): rms(DX)/delx = 1.021、rms(DY)/dely = 1.043 (統計誤差 約 2.6%)。k64: 0.713、0.746。
- 相関関数: ⟨c(s) c*(s+Δ)⟩/⟨|c|²⟩ = exp(−Δ²/(2x0²)) が 5% 以内、虚部 2% 以内 (規格化の定数によらない。⟨|c|²⟩ = 2)。
- 相関長 1000 km: リング全体が同じ値 (差 < 1e-9)、その値は rms 1.0 ± 5% のガウス分布。
- corlength = 0 の枝は v2 と k64 でビット単位で同じ (変えていない)。

- 使える範囲と条件: v2 モード、corlength > 0。乱数の列・スペクトルの形 (B62 の対称な置き方)・対象要素 (B40) は前任 C の v2 のまま。
- 使えない範囲・限界: rms は多数回の平均としての値。GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc、nc ≥ 128、INS は周長の 2 倍) より短いと線形補間で rms が下がる (間隔 = 相関長/4 で約 1%)。corlength = 0 の枝は別の量 (分散 = delx²·s(1−s/周長) [m]) で、delx の意味は「1 m あたりの rms」。
- DELCOR の CAVI の符号: **決定済み** (v2 の +delx のまま。変更していない)。
- 既存の試験への影響: 前任 C の `test_v2_C_delcor.py::test_b62_correlation_function_is_gaussian` の 38 行 `assert abs(c0 - 1.0) < 0.05` は旧い規格化 (⟨|c|²⟩ = 1) を確かめているので失敗する (c0 = 1.989)。`abs(c0 - 2.0) < 0.1` に直す必要がある (他人の試験なので触っていない)。`example/v2_C/C6_misc.pysad.txt` の DELCOR の行の値も √2 倍になる。

## 3. `FFS` の文の読み方 (pysad の不具合。モードによらない)

gs で確かめた動作 (小さな台本 15 本):

| 台本 | gs | pysad (修正後) |
|---|---|---|
| `FFS;NPARA=3;` 改行 `USE RING; CELL; CALC;` | FFS に入る。同じ行の `NPARA=3;` は FFS 部では実行されない (NPARA は変わらない。gs は FFS を抜けた後に MAIN の文として読む)。FFS 部は次の行から | 同じ (残りは捨て、捨てたことを 1 行で知らせる) |
| `FFS USE RING;`、`FFS USE=RING;` | FFS に入って RING を使う | 同じ |
| `FFS USE RING; NPARA=3; CELL; CALC; Print[…]` | USE RING だけ。残りは実行されない | 同じ |
| `ffs;`、`  FFS  ;  `、`FFS;   ! comment` | 読める | 同じ |
| `LINE RING = (…); FFS; Print["Z"];` (行の途中) | FFS に入る。Print は実行されない | 同じ |
| `FFS` 改行 `USE RING;` | `FFS USE RING;` と同じ (文は `;` まで) | 同じ結果 |
| `FFS NPARA=3;` | `*** $FFS *** syntax error:undefined variable $NPARA` を出して FFS に入る | 同じ (文の後ろに書ける形の説明つき) |

- 残した違い: (a) `Ffs;` (大文字小文字まじり) は gs では読めない。pysad は従来どおり読む。(b) `FFS` (`;` なし) の次の行が USE 以外の文のとき、gs はその文を FFS の引数として読んでエラーにし、その行の残りを捨てる。pysad は従来どおり次の行から FFS 部として実行する。(c) `ffs use ring;` (小文字の名前) は gs では RING が見つからない。pysad は USE と同じ扱い。(d) gs は捨てた残りを FFS の終了後に MAIN の文として読む (ふつうエラー)。pysad は読まない。
- 使える範囲: `run_sad_string` / `run_sad` で読む台本。試験: `test_ffs_statement_forms` (8 通り)、`test_ffs_in_the_middle_of_a_line`、`test_rest_of_ffs_line_is_not_executed`、`test_ffs_with_unknown_argument_is_reported`、`test_ffs_word_inside_other_text_is_not_the_command`。gs の参照値: `Twiss["BX",1]` = 3.6391435504154614 (pysad 3.6391435504154646)。
- 評価器の変更 (`ffscmd.py`): `_locate_ffs(src)` を新設 (注釈を空白にした文字列を走査し、行の先頭か深さ 0 の `;` の直後の `FFS` を探す)。`split_main_ffs` と `main_part` はこれを使う (戻り値の形は同じ)。従来の「行全体が `FFS;?`」の場合の結果は同じ。`FFS USE 名前;` は FFS 部の先頭に `USE 名前;` を置く。

## 4. MAIN 部の一般の代入 `NAME = 式;` (pysad の不具合。モードによらない)

gs で確かめた動作:
- `A = 1.5;` の後の要素定義 `L=A` は定義の時点の値 (後で `A = 9;` としても要素は 1.5)。FFS に入ると `A` は最後の値 (9)。FFS での `A = 4;` は要素を変えない。
- 名前は大文字小文字を区別する (`aa = 7;` は FFS で `aa` = 7、`AA` は未定義)。
- 式: 数 [単位] (`10 DEG`、`2 GEV`、`5 MM`、`1.5 M + 20 CM`、`2 deg`)、`+ - * /`、`^` (右結合: `2^3^2` = 512)、`( )`、先頭の符号 (`-A*2+D`、`-A^2` = −A²)、関数 `SIN COS TAN ATAN SQRT EXP LOG LN COSH` (丸括弧。`LOG` は常用対数、`LN` は自然対数。小文字 `sin(1)` も可)、`PI`、`MOMENTUM`、前に代入した変数。
- gs で構文エラーになる形: `Sqrt[16]`、`ABS(…)`、`ASIN(…)`、`ARCSIN(…)`、`2*-A`、`2^-1`、`2 A`、`A**2`、`pi`・`Pi`、未定義の名前、先頭が括弧の `(1+2)*3`。gs は変数を未定義のままにするか、前の式の値などを入れる。
- 特別な変数: `PBUNCH EMITX EMITY EMITZ SIGE SIGZ MINCOUP LOSSAMPL LOSSDZ FSHIFT TDXI TDYI TDZI MAXORDER ADDTERMS MASS EMITDIV EMITDIVB EMITDIVQ EMITDIVS FRINGDIV NBUNCH PHSPLOTS` は FFS でも MAIN の値。`NPARA`・`DP` (0.01)・`DP0` (0)・`CHARGE` (1) は FFS の入口で既定値に戻る。先頭の MARK に書いた EMITX などは MAIN の値より優先、MARK に無い (0 の) ものは MAIN の値が残る。MAIN の `E = -A;` は FFS の定数 E を変えない。

pysad (修正後): 上と同じ。gs で構文エラーになる形は代入しない (`(1+2)*3` だけは読む)。関数形・PI を含む代入は、要素定義でも使えるように読み取りの前に数値に書き換える (`F9 = SQRT(A+0.75)*2; DRIFT D3 = (L=F9)`。修正前は ValueError で台本全体が読めなかった)。TRACK のある台本 (MAIN の命令を順に実行する経路) でも同じ。

- 試験: `test_main_assignments_reach_elements_and_ffs`、`test_main_expression_grammar_matches_gs` (32 個の式を gs の値と 2e-14)、`test_main_expressions_gs_rejects_are_not_assigned`、`test_main_special_variables_reach_ffs`、`test_main_assignments_with_track_statement`、`test_parse_main_assignments_units`。例: `example/v2_C/C2w_main_script.sad` (gs の出力と全行同じ。末尾の注釈)。
- 使えない範囲・残した違い: (a) 要素定義の中に直接書いた MAIN の関数 (`L=2*SIN(A)`) は従来の読み取りのまま (未確認)。(b) `A = 3 B = 4;` は gs では A = 3 (エラーつき)、pysad は代入しない。(c) `LINE[{"EMITY",…},1]` (MARK の値の問い合わせ) は gs では MAIN の値を返すが、pysad は MARK に書いた値 (0)。大域変数 EMITY は合う。(d) `SINH`・`ACOS` などここに無い関数は gs で確かめていない (pysad は代入しない)。(e) `_MAIN_SET_NAMES` に無い FFS の変数 (GCUT など) を MAIN で代入した場合は、一般の変数と同じく FFS に代入する (gs での確認はしていない)。
- 評価器の変更: `parser.py` の末尾に `_blank_comments`、`_eval_main_expr`、`parse_main_assignments`、`rewrite_main_assignments` を追加 (既存の関数は変えていない。`parse_main_commands` の戻り値も同じ)。`ffscmd.py` の `run_sad_string`: 読み取りの前に代入を拾い、FFS を作った後に一般の変数を `I.assign` で入れる。TRACK の無い経路では、`_MAIN_SET_NAMES` の値を (EMITDIV 類だけでなく) NPARA・DP・DP0・CHARGE 以外すべて FFS に渡す。`FFS.use` の最初の 1 回の `{EMITX,…} = 先頭 MARK の値` は、MARK に無い (0 の) とき MAIN の値を残す。
- k64 の結果への影響: MAIN で PBUNCH・EMITX・MASS などを代入している TRACK なしの台本は、これまで無視されていた値が効くようになる (gs と同じ向き)。走らせた既存の試験 (下) では変化なし。

## 走らせた既存の試験 (修正後)

`test_main_track`、`test_ffs_gs_match`、`test_ffsextra`、`test_scriptgs`、`test_emit_report`、`test_emit_extra`、`test_emit_unstable`、`test_sol_emit`、`test_v2_A_ffs`、`test_v2_B_ffs`、`test_v2_B_script`、`test_delcor`、`test_v2_C_delcor`、`test_v2_C_wakecmd`、`test_wake`、`test_wake_optics`、`test_example_fixes`、`test_evalfixes`、`test_new_features`、`test_small_items`、`test_pysad`、`test_item4`、`test_fileops`、`test_symbol_names`、`test_sadclass`、`test_apert`、`test_corr`、`test_error_spreading`、`test_synchrobeta`、`test_trpt_acc`: 317 通過・1 skip (ほかに重複して走らせた分あり)、失敗は上の `test_b62_…` の 1 件 (規格化の変更による。§2)。全試験と回帰は走らせていない。
