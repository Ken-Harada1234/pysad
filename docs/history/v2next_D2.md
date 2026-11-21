# v2next の記録: 担当 D2 (パーサ・GUI・パッケージ・整理、2026-11-18)

v2.0.0 の記録と `v2next_B2.md`・`v2next_C.md` の「要判断・未確認」のうち、D2 の 7 項目を調べ、直した。
gs は SAD K64 (1.4.4.2k64、2026-10-19 ビルド、WSL)。`docs/46`・`docs/47`・README・CHANGELOG は触っていない (`docs/21` の B 番号も振っていない: 下の「K64 の誤り」に `D2-1` と仮の名前で書いた)。
走る例: `example/v2_next/D2_*` (台本・gs の出力・pysad v2・pysad k64)。説明は `example/v2_next/D2_README.md`。
試験: `pysad/tests/test_v2_D2_parser.py` (58 件)、`pysad/tests/test_v2_D2_packages.py` (9 件。うち 2 件は `PYSAD_SAD_PACKAGES` が無いと skip)。既存の試験 2 件を直した (`test_tclbridge.py::test_tk_version_and_font_system` と `test_v2_C2_script.py::test_parse_main_assignments_units`)。

## 使ったリング

| 名前 | 条件 | 使った項目 |
|---|---|---|
| 小さなリング | 12 セルの FODO (3 GeV、周長 52.8 m、BPM 24 台、ステアラ各 12 台。`D2_micado.sad` ほか)、24 セルの FODO (3 GeV、周長 105.6 m。`example/13_packages` のリング)、12 セルの FODO + RF (3 GeV、周長 49.2 m。`test_spin.py`) | (3) (4) (5) |
| 大きなリング A | 3 GeV、周長 603 m、3563 要素 (`pysad/tests/kek-ls_v03_12.lat`) | (1) の読み取りが従来と同じ (要素の値が全部ビット一致) |
| PF リング・ASTRID2 | `example/pf/lattice/v9_16_asym.sad`・`v9_16_inj_noSF.sad`、`example/astrid/lattice/ASTRID2_v01_01.sad`・`ASTRID2_v01_02_id.sad` | (1) の読み取りが従来と同じ |
| SP8 (大きなリング B、8 GeV) | 周長 1436 m、4 回対称、矩形の偏向磁石 88 台、RF 17 MV (SPring-8 型、デザインレポートより再現。数値だけ) | (7) |

## 1. 一覧

| 項目 | 状況 (調べたこと) | K64 (gs) の動作 | v2 の動作 | 根拠 | 試験 | k64 モード | 使える範囲と条件 | 使えない範囲 |
|---|---|---|---|---|---|---|---|---|
| (1) MAIN の要素定義の式 (`ANGLE=PI/8`) | 要素定義の値は、数と変数・文字列を読む一般の式の読み取り (`parse_expr`) を通っていて、`PI`・関数 (`SIN(...)`)・`2*PI` を文字列として扱い、`[pysad parser] Warning … unsupported operand type(s) for /: 'str' and 'int'` を出して要素を捨てていた。gs で約 100 通りの式を実測して文法を決めた (下の §2.1) | `PI` (大文字だけ)・`SIN COS TAN ATAN SQRT EXP LOG(常用) LN SINH COSH TANH ATAN2(y,x)` (関数名は大文字小文字を問わない、丸括弧)・`+ - * / ^` (`^` は右結合)・括弧・先頭の符号・単位 `DEG MM CM …` (**数の直後だけ**。`10 DEG/2` は (10 DEG)/2、`10/2 DEG` は 10/(2 DEG))・前に代入した変数。0 での割り算は ±Infinity・NaN、定義域の外は NaN (エラーにならない)。`pi`・`Pi`・`ABS ASIN ACOS LOG10`・`SIN[…]`・`2**3`・`2^-1`・`2*-1`・`--1`・`1.5 PI`・`(1+2) MM`・`SQRT 2`・引数の数が違う関数・未定義の名前は `LALR(1) parser: syntax error` | 同じ (gs と同じ文法・同じ値。pysad 自身の不足なので k64 と同じ)。`pysad/parser.py` の `_main_expr_parse` (一般の代入 `A=PI/8;` の評価と共通) を要素定義の `KEY=式` でも使う。gs が構文エラーにする形は、警告 `syntax error in MAIN expression (gs: …)` を出してその定義を捨てる | gs の 37 通りの値 (`Element["ANGLE","B1"]` など) と 19 通りの構文エラーを `test_v2_D2_parser.py` に写した。`.39269908169872414` など全桁 (rel 1e-15) 一致 | `test_v2_D2_parser.py` (58 件) | 同じ (ゲートしない) | `D2_main_expr.sad` (13 の ANGLE・K1・E1・ROTATE が gs と全桁一致)。大きなリング A、PF、ASTRID2、SP8 (大きなリング B)、24 セル FODO の全要素の値を、v2.0.0 のコミットの読み取りと比べてビット一致 | gs が構文エラーにするが pysad が受ける形: `45DEG` (数と単位の間の空白なし。字句に空白の情報が無い)、`SIN[0.1]` (従来の読み取りへ)、未定義の名前 `L=A` (従来の読み取りで文字列になり、あとで型エラー)。gs が受けるが pysad が構文エラーにする形: 末尾の `1+` (gs は 1)。MAIN の変数を `E`・`Pi`・`Degree` の名前で代入しても、要素定義の式の `E` は SADScript の定数 (gs は利用者の値) |
| (2) Python 3.14 (Tk 9.0) の GUI | `py314` の venv (numpy 2.5.5・numba 0.68.1・scipy 1.18.1) で `test_tclbridge.py`・`test_plot.py`・`test_gui_lang.py` の 25 件と、`example/12_gui_tkinter` の 6 本を実行。ウィジェットの木・幾何・フォント・色を 3.12 (Tk 8.6.13) と 3.14 (Tk 9.0.4) で全ウィジェットについて比べた | — (Tk の版の話) | 失敗は `test_tk_version_and_font_system` の 1 件 (Tk の版が `8.6` で始まることを要求) だけ: `8.6` か `9.0` を受けるように直した (判定の意味は同じ。`$FontSystem == 'X11'` の判定も同じ)。ほか 24 件は 3.14 で通る。6 本の例の出力は Tk の版の行 (TkinterCore) 以外が 3.12 と同じ。見た目の差は §2.2。**`tclbridge.py`**: Windows の Tk 9.0 では SAD の既定の固定幅フォント `courier` が Courier New に置き換わらない (行の高さ 13 px の Courier。Tk 8.6 は Courier New の 16 px) ので、`$DefaultFixedFont` を `Courier New` にする | 幾何の JSON (全ウィジェット) の比較 (§2.2)。Tk 8.6 では `font actual {helvetica 9 bold}` が Arial、9.0 では `helvetica` を返すが、`font metrics`・`font measure` は同じ (15 px・65 px) | `test_tclbridge.py::test_tk_version_and_font_system` (直した)、3.14 での 3 ファイル 25 件 | `tclbridge.py` の変更は Tk 9.0 の Windows だけ (3.12 では何も変わらない) | Python 3.14 (Tk 9.0.4)・Windows。SAD の `.n` の GUI (Window・Button・Entry・ListBox・Scale・Canvas・Menu・KBMainFrame・KBFComponentFrame) | 実際の画面の見た目の撮影 (ImageGrab が画面を取れず、`py314` の venv に PIL が無い)、マウス操作、Tk 9.0 以外 (macOS・Linux) |
| (3) `test_spin.py::test_mult_ring_matches_gs` が 3.14 だけ失敗 | 3.14 と 3.12 で閉軌道の反復を、`_one_turn_cod_matrix` の入出力 (cod0・cod・一周行列の 4x4) と `_solve_cod` の入出力を 16 進で記録して比べた。**反復は収束まで 1 ビットも違わない** (3.12 の記録 21 行 = 一周行列の評価 12 回と solve 9 回が、3.14 と 1 ビットも違わない)。閉軌道は原因ではない。失敗は `emittance6d` の例外を `run_emit` が `closed orbit not found` と表示して握りつぶしたもので、本当の例外は `spin.py` の `_spintrm` (spdepol) で出ていた: 不安定なリングで `cb = −_cexp1(−cφ)/a` が inf を含み、(a) 3.14 の複素数の割り算・掛け算は inf を C99 付属書 G に従って残す (3.12 以前は inf を含む複素数の `/` は nan+nanj)、(b) 3.14 の `complex ** 2` は桁あふれで `OverflowError`、(c) `0.5 / x` が `ZeroDivisionError` | `cb` が inf のとき `_cerfc(cb)` は 0 割りと桁あふれを起こす (Fortran の複素数は IEEE 例外なし) | 同じ (k64 の値は 3.12 と全桁同じ) | (i) `cb` が有限でなければ `nan+nanj` に置き換える (3.12 の結果と同じ)、(ii) `_cerfc` の `z ** 2` を `z * z` に | 3.14 と 3.12 の出力が k64・v2 とも一致 (`PolarizationTime` 133.6117332301696、スピンチューン 0.40379908026906、偏極ベクトルまで)。3.14 と 3.12 で `test_spin.py` 8 件通過 | `test_spin.py` 既存の 8 件 (3.14 と 3.12 の両方で確認) | **変えない** (3.12 の結果を保つ修正) | Python 3.12・3.14 | 3.10 は確認していない (同じ nan+nanj の分岐なので通るはず)。握りつぶした例外の表示が `closed orbit not found` なのは K64 の動作 (`Emittance[]` の NaN の報告) なので変えていない |
| (4) `PeriodicR`・`CorrectOrbit` の Micado・`CorrectLocal` | 小さな台本を gs と pysad (k64) で走らせて比べた。**PeriodicR**: 4x4 の一周行列 6 通り (結合あり・なし・半整数に近い・NaN になる 2 通り) で R・A・B の 12 成分が 3e-16 以内。**CorrectLocal**: 12 セル FODO、領域 {BA3, BA8}・垂直、補正量 12 個と 24 BPM の DY が 1e-14 以内。領域の外の軌道は 1e-11 以内で変わらない (局所バンプが閉じる)。**Micado**: 後述 | PeriodicR は `Optics.n` が読み込まれてから使える (最初に `PeriodicOptics` を呼ぶ。先に呼ぶと未評価。gs も同じ)。**Micado**: `pmicad.f` が `tsvd(a,b,x,eps,svd)` を呼ぶが、新しい tsvd は配列の大きさを `size()` で取る (旧 `tsolva(a,b,x,nx,nx,mdim,eps)` は nx を渡していた)。ハウスホルダー変換したあとの全 n 列・全 m 行を解いて、先頭 nx 個だけを使う (選ばなかった補正量の寄与が落ちる)。24 セル FODO (BPM 48 台、垂直ステアラ 24 台) の補正後の垂直 rms は補正前 9.41e-5 m に対して最小二乗 2.06e-5、Micado n=1 は 1.46e-4 (**悪化**)、n=3 8.2e-5、n=6 1.19e-4、n=8 1.75e-4 (本数を増やしても単調に下がらない)。pysad (k64) は gs と 1e-6 以内で同じ | v2 は上三角の nx x nx だけを解く (`pmicad` の `solvf=True` の分岐)。24 セル: n=1 4.74e-5、2 4.43e-5、3 4.37e-5、4 4.31e-5、6 4.11e-5、8 3.97e-5、12 3.51e-5、24 2.06e-5 (最小二乗と同じ) と**単調に下がる**。`Condition` つき (CorrectLocal) の Micado は、選んだ列だけで `psvdcond` を解く元の分岐なので K64 でも正しい | 独立な計算: 選ばれた nx 列の最小二乗 `numpy.linalg.lstsq` と v2 の解が 1e-9 で一致、選び方は独立な貪欲法 (各段で残差が最小になる列) と同じ [22,0,20,21,17,23]、本数を増やすと残差が減り、全部使えば全体の最小二乗に一致、線形モデルの残差の減り (2.27e-4 → 8.8e-5) と実際の rms の減りが整合 | `test_v2_D2_packages.py`: `test_micado_v2_is_least_squares_on_selected_correctors`、`test_micado_k64_keeps_truncated_full_least_squares`、`test_micado_24cell_ring_k64_matches_gs_and_gets_worse`、`test_micado_24cell_ring_v2_decreases_monotonically` | 従来どおり (gs と一致) | v2 の Micado (`Solver->{"Micado",n}`、`Condition` なし)。24 セルと 12 セルの FODO の垂直 | 水平・2 平面同時の Micado、`CorrectOrbit` の `Bump`・`SVDmode->1` つきの Micado は試していない |
| (5) `CSROYMatrix` の wrongtype (B60)・`tffamsetup` (fam 点) | **CSROYMatrix**: K64 の `tfcsroysetup` (csroy.f) を読み、gs で `CSRSetupOY` の出力の各要素の `Head` を調べた。**4 番目 (zj) の `Head` が実数 (1.2952814481588496e-5)** で、`List` でない。`kazjn = ktavaloc(0,nj*mphi)` の実数リストに `c_f_pointer(c_loc(rlist(kazjn)), zjn, …)` で書き込むとき、先頭の語 `rlist(kazjn)` (リストの頭 = head) から書き始める (要素は `rlist(kazjn+1)` から。ほかの 4 つの配列は `rlist(k+1:…)` で正しい)。head が zjn(1,1) で上書きされるので、`tfcsroymat` の `tfreallistq(las%dbody(4), lazj)` が偽になり、`General::wrongtype` を返す。**手で作った正しい形の実数リストなら gs でも動く** (docs/21 B60 の「手で作った実数リストでも」は誤り) | `CSRSetupOY` の出力を `CSROYMatrix` に渡すと常に wrongtype (zj の head の上書き + B59 の `dpjn` が `zjn` と同じ記憶域)。手で作った入力 (nr=64、nj=3、nphi=2) では 6x6 の行列 (36 要素) を返す | pysad の `CSROYMatrix` は両モードで動く (従来どおり。k64 でも `CSRSetupOY` の出力を受ける: gs の wrongtype は再現しない)。手で作った入力の行列は gs と相対 4e-20 で一致 | gs の 36 要素 (`D2_csroymatrix.gs.txt`) と pysad の比較 | `test_csroymatrix_hand_made_input_matches_gs[k64・v2]` (相対 1e-12) | **変えない** (K64 の wrongtype を再現する修正はしない。k64 は従来どおり動く) | `CSROYMatrix` の行列 (`csroymat` の翻訳) は、手で作った入力 (nr=64、nj=3、nphi=2) で gs と一致。`OYStability` が使うのは `CSRSetupOY` の出力なので、gs では完走しない (K64 の誤り `D2-2`: 下の §3) | `CSRSetupOY` の出力から作った行列の gs との比較 (gs が作れない)、nj・nphi が大きい手作りの入力、`OYStability`・`OYScan` の gs との比較 |
|  | **tffamsetup**: K64 の `tffamsetup.f` (184 行) を全部読み、pysad の `match6.py` の `_fam_define` と式ごとに照合 (dx(1)=c x0、dx(2)=s px0 − α/β dx(1)、tfnormaltophysical の R1〜R4、dp の区間 dpm、添字の並び nfa=merge(−nfa, −nfa+1, nfa≥0)、nfp の上限 12)。gs は FITP>1 の GO で再起動するので値は比べられない。独立な解析の基準を 2 つ作った | — | — | (a) 初期条件の Courant-Snyder 不変量 (x² + (αx+βpx)²)/β = xamp² · em (x と y の各面)。(b) MARK の直後に厚さ 0 の八極 (K3L = 6000) がある線形ラティス (QF D QD D の 8 セル、1 GeV) で、fam 点のチューン (固定軌道のまわりの一周行列) を、八極のキックのヤコビアン J = −K3L (x² − y²)/2 を一周行列に掛けた式 cos μ = cos μ0 + J β sin μ0 / 2 と比べる。x の振幅 3σ の点: νx 0.92315 (式 0.92303)、νy 0.63330 (式 0.63320)、y の振幅の点: νx 0.87218 (式 0.87212)、νy 0.67207 (式 0.67199)。差 1e-4 以内は、八極の前後の厚い四極の非線形項 (px² のオーダー) による | `test_fam_initial_conditions_have_courant_snyder_invariant_of_emittance_times_amplitude`、`test_fam_point_tunes_match_thin_octupole_kick_jacobian` | 同じ | fam 点の初期条件と、そこでの一周行列・チューン (結合なし・八極 1 個の線形ラティス) | 足りない独立な基準: (i) 有限振幅の周期解 (固定点) の条件そのもの (tffscalc の fam ループが FIT の残差に入れる値)、(ii) 振幅依存チューンシフトを FIT の目標にした収束の結果 (gs は FITP>1 で止まる)、(iii) 結合・6 極のあるラティスの fam 点。必要なもの: 追跡の NAFF によるチューンか、6 極の解析的な振幅依存チューンシフト |
| (6) 整理 (改行、numba のキャッシュ) | §4 | — | — | — | `test_nbcache_purges_caller_when_callee_file_changes` | — | — | — |
| (7) SP8 (大きなリング B、8 GeV) | §5。読み込み・CALC・RADINT・Emittance[]・COD が動く。**ラティスの文法で pysad が未対応のものは無かった** (修正なし) | — | — | gs と k64 の比較 (§5) | (数値だけ。試験・例には使わない) | 従来どおり | SP8 (大きなリング B) の CALC・RADINT・Emittance[] (RFSW)・COD (四極の DX・DY・ROTATE) | NORFSW + RAD の Emittance[] で軌道がずれたとき gs と 2〜4% 違う (§5、未調査) |

## 2. 補足

### 2.1 MAIN の式の文法 (gs で実測、2026-11-18)
`BEND B1=(L=1 ANGLE=式)` の `Element["ANGLE","B1"]` を、式約 100 通りで gs に読ませた。
- 値が出る形: `PI/8`=.39269908169872414、`2*PI`、`-PI/8`、`PI/8 DEG`=22.5 (単位は 8 だけにかかる)、`45 DEG`、`1 MM`、`2 CM + 1 MM`=.021、`SIN(0.1)`・`sin(0.1)`、`SQRT(2)/10`、`COS(PI/3)`=.5000000000000001、`TAN(PI/8)`、`ATAN(1)`、`EXP(-1)`、`LOG(100)`=2 (常用)、`LN(10)`、`SINH`・`TANH`・`COSH`、`ATAN2(1+1,2)`・`ATAN2(PI,2 DEG)`、`2^3`、`2^3^2`=512 (右結合)、`-2^2`=−4、`-2*3`、`SIN(10 DEG)`、`SIN(0.1)^2`、`-SIN(0.1)`、`SIN (0.1)`、`2*(1+2)`、`1/0`=Infinity、`-1/0`、`0/0`=NaN、`LN(0)`=−Infinity、`LN(-1)`・`SQRT(-1)`・`LOG(-1)`=NaN、`2^1000000`=Infinity、`0^0`=1、`(-8)^(1/3)`=NaN、`1E400`=Infinity、`1E-400`=0、`TAN(PI/2)`=16331239353195370。
- 構文エラー (`LALR(1) parser: syntax error`): `pi/8`・`Pi/8` (`PI` は大文字だけ)、`45DEG` (数と単位の間に空白が要る)、`(1+2) MM`・`(PI/8) DEG`・`PI DEG`・`1.5 PI` (単位は数の直後だけ)、`SIN[0.1]`、`SQRT 2` (丸括弧が要る)、`ABS`・`ASIN`・`ACOS`・`LOG10`・`ERF`・`GAMMA`・`FLOOR`・`MIN`・`MAX`・`SIGN` など未定義の関数、`2**3`、`2^-1`、`2*-1`、`--1` (演算子の直後の符号)、`ATAN2(1 2)`・`ATAN2(1,2,3)`・`SIN(1,2)`・`LOG(2,3)` (引数の数)、`2 SIN(0.1)`・`SIN(0.1) DEG`、`SIN()`、`(1`、未定義の名前 `A`。
- 例外: `1+` は gs が 1 を返す (構文エラーにならない)。pysad は構文エラーにする。
- MAIN の変数 (`A=PI/8; S=SINH(.1); C=1/0; D=ATAN2(1,2);` のあと `ANGLE=A`・`F+A`・`G` (G=10 DEG)) は値が引き継がれる。変数名は大文字小文字を区別 (`a` は未定義で構文エラー)。
- 従来の (`C2`) 一般の代入の評価 `_eval_main_expr` は、関数が `COSH` までで `SINH TANH ATAN2` が無く、0 での割り算と `SQRT(-1)` を構文エラーにしていた (gs は Infinity・NaN)。同じ文法の評価器 `_main_expr_parse` に統合し、`test_v2_C2_script.py::test_parse_main_assignments_units` の 2 行 (`1/0` と `SQRT(-1)` が ValueError) を gs の値に直した。
- 実装: `_main_expr_parse(toks, vars_, getvar, stream)` が字句 (`num`・`id`・`op`) の並びを読む。`stream=True` は式が終わった所で止まる (要素定義の `ANGLE=PI/8 K1=0` は次のキー `K1` で止まる)。要素定義では `_parse_sad_params` が、値が数だけ (`L=1`) なら従来どおり整数のまま、式なら新しい読み取りに渡す。識別子だけの値 (`L=name` の参照) や文字列・リストは従来の読み取り。MAIN の `/` (`_parse_mul`) も 0 での割り算を ±Infinity・NaN にした。

### 2.2 Tk 8.6.13 と 9.0.4 の差 (全ウィジェットの幾何・オプションの比較)
6 本の例 (`fn_Canvas`・`fn_Entry_ListBox_Scale`・`fn_KBMainFrame_KBFComponentFrame`・`fn_Menu_MenuButton`・`fn_Window_Button_TextLabel`、`fn_TkinterCore` は窓なし) で、ウィジェットごとに class・幅・高さ・位置・要求サイズ・font・bg・fg・relief・bd・padx・pady・text を取って 3.12 と 3.14 で比べた。
- **標準出力**: Tk の版の行 (`fn_TkinterCore`: 9.0.4 と 8.6.13) 以外は全部同じ (`.pysad.txt` と一致)。
- **オプションの値に単位が付く**: Tk 9.0 の `cget` は長さを pt で返すことがある (`-padx` 4p・`-pady` 3p (8.6 は 5・4)、`Scale` の `-width` 11.25p (15))。pysad (SAD の `TclReturn`) は文字列のまま渡す。この例と試験の範囲では数値として読む所は無かった。SAD の `.n` が `cget` の値を数として使うと、Tk 9.0 では文字列になる (**未確認・直していない**)。
- **サイズ**: Checkbutton・Radiobutton の幅が Tk 9.0 で 2〜3 px 大きい (インジケーターの寸法の変更)。Menubutton の padx 5 px → 4 pt (5.3 px)。
- **フォント**: Tk 9.0 の `font actual` は要求した family 名をそのまま返す (`helvetica`。8.6 は Arial)。`font metrics`・`font measure` は同じ (15 px・65 px) なので見た目は同じ。`courier` は 9.0 では Courier New に置き換わらず、行の高さ 13 px の Courier になる (8.6 は Courier New の 16 px)。`$DefaultFixedFont` を `Courier New` にして直した (`tclbridge.py`、Tk 9.0 の Windows だけ)。Windows に無い family (`helvetica neue`・`utopia`・`fixed`) は 8.6・9.0 とも既定のフォントに落ちる。
- **色**: `selectBackground`・`activeBackground` が違って見えるのは、`Tkinter.n` が日時で選ぶ (`Tkinter$BG[]`) 乱数の色で、Tk の版の差ではない。
- ルート窓の画面上の位置は、窓の出る位置の既定の違い (26 px・78 px … の増分)。

## 3. K64 の誤り (D2 が見つけたもの)

| 仮の番号 | 内容 | 原因 | v2 |
|---|---|---|---|
| D2-1 | `CorrectOrbit` の `Solver->{"Micado",n}` が、選んだ nx 本の補正量でなく、全補正量の最小二乗解の先頭 nx 個を使う。補正後の rms が補正前より悪化し、本数を増やしても単調に下がらない | `pmicad.f`: ハウスホルダー変換の後 `call tsvd(a,b,x,eps,svd)`。新しい tsvd は `size()` で n・m を取るので全 n 列を解き、x(nx+1:n)=0 にする (コメントアウトされた旧 `tsolva(a,b,x,nx,nx,mdim,eps)` は nx を渡していた) | 修正 (`corr.py` の `pmicad`、`compat.v2()`) |
| D2-2 | `CSRSetupOY` の出力の 4 番目 (zj) が実数リストでなくなり (head が zjn(1,1) で上書き)、`CSROYMatrix` が wrongtype になる。docs/21 B60 の原因。`OYStability`・`OYScan` は gs で完走しない | `csroy.f` の `tfcsroysetup`: `c_f_pointer(c_loc(rlist(kazjn)), zjn, [mphi,nj])` が `rlist(kazjn+1)` でなく `rlist(kazjn)` (head) から。B59 の `dpjn` と `zjn` が同じ記憶域も同じ箇所 | 直していない (pysad の `CSRSetupOY` は両モードで正しい形の出力を返し、`CSROYMatrix` は動く。k64 で wrongtype を再現しない) |

## 4. 整理 (項目 6)

### 4.1 改行 (LF/CRLF)
`git ls-files --eol` で、HEAD (index) と作業コピーの改行が違うファイルを洗い出した (`git diff --stat` と `git diff -w --stat` の行数が大きく違うもの)。内容は変えず、HEAD の改行に戻した。
- HEAD が CRLF で作業コピーが LF: `pysad/perfkern.py` (761 → 9 行の差)、`pysad/track6n_ext.py` (2814 → 42 行) → CRLF に。
- HEAD が LF で作業コピーが CRLF: `pysad/radtaper.py`、`pysad/tests/test_v2_A_spin.py`、`test_v2_C2_optics.py`、`test_v2_C2_spch.py`、`test_v2_C_track.py` → LF に。
- 直した結果、`git ls-files --eol` の不一致 (i/crlf w/lf、i/lf w/crlf) は 0 件。残る `git diff --stat` と `-w` の差 (elkern・emit6d・ffscmd・fileops・lattice・track6・track6n で 4〜44 行) は改行でなく、行の中の空白の変更 (他の担当の編集)。
- D2 が編集したファイル (`parser.py`・`spin.py`・`corr.py`・`tclbridge.py` は CRLF、`test_v2_C2_script.py` は LF、`test_tclbridge.py` は元から混在) は、元の改行のまま。

### 4.2 numba のキャッシュ
- numba 0.68.1 の `cache=True` の索引 (`.nbi`) は、**関数を定義したファイル 1 つの内容のハッシュ**で無効になる (更新時刻ではない。`numba.core.caching` の `get_source_stamp` が SHA のバイト列を返す)。`touch` だけ (内容が同じ) では捨てられず、内容を 1 行変えると捨てられる。
- **他のファイルの njit 関数や定数を使うカーネルは、呼ばれる側の変更で無効にならない**。実験 (`moda.f` が `modb.g` を呼ぶ): g を書き換えても、f の値は古い (f=202 のまま、g 単独は 1001)。f のファイルを 1 行書き換えると直る。大域の定数 (`emit6d.ELMASS` など) もコンパイル時に埋め込まれるので同じ。
- pysad で当たるもの: `elkern.py` (← `perfkern.py`・`perfkern_bend.py`)、`turnkern.py` (← `elkern.py`・`perfkern*.py`)、`multa.py`・`tcav.py`・`und.py`・`beambeam.py`・`phsrot.py`・`track6n_ext.py` (← `track6n.py`、`emit6d.py` の定数)、`wake.py` (← `emit6d.py` の定数)。`perfkern.py` を編集すると、`elkern`・`turnkern` のキャッシュが古いまま残る。
- 対策: `pysad/nbcache.py` を足した (`py -3.12 -m pysad.nbcache`)。`__pycache__/nbdeps.json` に、モジュールごとの依存先 (njit を持つモジュールの内容のハッシュ (改行の違いは無視)、import した大文字の定数の代入行のハッシュ) を覚え、変わっていたらそのモジュールの `.nbi`・`.nbc` を捨てる。記録が無い最初の 1 回は記録だけ (最初から古いものがあるかもしれないときは `--force`)。**kernel (`elkern.py`・`perfkern*.py`・`turnkern.py`・`track6n*.py`) を編集したら、試験の前に 1 回走らせる**。numba に自動の仕組みは無い (自動で呼ぶ仕組みは入れていない: kernel を持つ共有ファイルを編集せずに済ませるため)。
- 試験 `test_nbcache_purges_caller_when_callee_file_changes`: 呼ぶ側だけ捨てる、使っていない定数の変更では捨てない、使っている定数の変更では捨てる。

## 5. SP8 (大きなリング B、8 GeV)

条件: 8 GeV、周長 1435.9488 m、4 回対称 (1 超周期に偏向磁石 22 台。矩形の偏向磁石 88 台、E1=E2=0.5)、薄い6 極 (K2 のみ)、RF 空洞 32 台 (各超周期に 8 台。合計 17 MV、h=2436)、BPM なし。ラティスは `MOMENTUM=8.0 GEV\n;`、`CAVI CAVI=(…)` (要素名が型名と同じ)、`LINE QUAD=(…)` (ライン名が型名と同じ)、`-NCELLH` (逆順)、`4*QUAD` を使うが、pysad は読める。(1) の変更後の読み取りを、v2.0.0 のコミットの読み取りと比べて、49 要素の値が全部ビット一致。

| 量 | gs (k64) | pysad k64 | pysad v2 |
|---|---|---|---|
| CALC: νx, νy | 47.1199940730654, 15.219999034381738 | 47.11999407306542, 15.219999034381747 | 同じ |
| βx, βy (先頭 MARK) | 19.99998259534472, 15.00002095997227 | 19.999982595346857, 15.00002095997252 | 同じ |
| 周長 [m] | 1435.9488 | 1435.9488 | 同じ |
| RADINT: I1x, I2, I3, I4u, I5u | 0.2096002, 0.1599915, 4.0739433e-3, −6.8003193e-5, 1.1853438e-5 | 同じ (7 桁) | 同じ |
| RADINT: εu, σδ, α | 6.9554e-9, 1.0936e-3, 1.4597e-4 | 同じ | 同じ |
| Emittance[] (RFSW): εx [m] | 6.956073911881967e-9 | 6.956073911881036e-9 | 6.965072e-9 (+0.13%) |
| σδ | 1.09399186e-3 | 1.09399186e-3 | 1.093640e-3 |
| 運動量圧縮 α (RFSW) | 1.459653277613363e-4 | 1.45965327763e-4 | 1.45917e-4 |
| 減衰率 (x, y, z) /周 | −5.77478e-4, −5.76477e-4, −1.15100e-3 | 同じ (1e-11) | −5.76732e-4, −5.76489e-4, −1.15174e-3 |
| COD (四極 QF5 に DX=5e-5、QD4 に DY=−4e-5、QF7 に ROTATE=2e-3) 7 点の x・y | 例: DX(1)=−1.1618633751724545e-3 | 1e-12 以内 | 同じ |
| 同じ誤差で Emittance[] (RFSW): εx, εy | 6.864299835776e-9, 4.438232533e-12 | 6.864299835589e-9, 4.438232534e-12 | 6.87296e-9, 4.43816e-12 |

- 所要: pysad の CALC + RADINT + Emittance[] は約 3.5 秒 (gs 0.5 秒)。
- k64 は gs と全部一致 (1e-10 以内)。v2 の εx が +0.13% 大きいのは、A1・A2・C の v2 の修正 (矩形の BEND の E1・E2 と放射、B07 など) による (sector 磁石 E1=E2=0 に変えると k64 と v2 が 9e-6 以内で一致する)。gs の Emittance[] (εx 6.9561e-9) と RADINT の εu (6.9554e-9) は 1e-4 で整合し、v2 の 6.9651e-9 は RADINT より 1.4e-3 大きい: **v2 の Emittance[] と RADINT の差 (矩形磁石) は、v2 の修正の側の問題かもしれない (未調査。A1 と判断が要る)**。
- **新しい発見 (未調査・直していない)**: `NORFSW; RAD` の `Emittance[]` で、閉軌道がずれる (四極・偏向磁石に DX) と、k64 が gs と 1.8〜4% 違う。小さな FODO (DX=1e-4 の QF): εx gs 2.62400554e-7、k64 2.57658900e-7。DX=2e-4 の B: gs 2.59994014e-7、k64 2.64626192e-7。SP8 (大きなリング B、QF5 に DX=5e-5): gs 7.1529e-9、k64 6.8685e-9。**RFSW が ON なら一致 (1e-10)**。軌道のずれがなければ NORFSW でも一致 (εx 6.956073911967e-9)。v2.0.0 のコミットの pysad でも同じ値で、他の担当の編集による変化ではない。RFSW OFF + RAD は縦方向が閉じない設定で、軌道の dp の扱いが gs と違う可能性。使えるのは RFSW ON。

## 6. 要判断・残した項目

| 項目 | 内容 | 判断・理由 |
|---|---|---|
| `45DEG`・`SIN[0.1]`・未定義の名前の受理 | gs は構文エラー、pysad は受ける (字句に空白の情報が無い、従来の読み取りへ) | 直していない。必要なら `tokenize` に位置を持たせる |
| MAIN の `E`・`Pi` などの代入 | gs は要素定義の式で利用者の値を使う。pysad は SADScript の定数のまま | 頻度が低い。直していない |
| 末尾の演算子 (`1+`) | gs は受理して 1 | 入力ミスなので受けない |
| Tk 9.0 の `cget` の単位つきの値 (`4p`) | SAD の `.n` が数として使うと壊れる | 今の例・試験の範囲では無い。`cget` の値を使う SAD の `.n` が見つかったら `winfo fpixels` で直す |
| D2-2 (CSRSetupOY の head の上書き) | gs の `OYStability` は完走しない | k64 で wrongtype を再現するかは判断が要る (今は再現しない) |
| `NORFSW; RAD` の `Emittance[]` の軌道ずれ (§5) | gs と 2〜4% 違う | A1 (放射・Emittance[]) の担当へ。RFSW ON は一致 |
| v2 の Emittance[] と RADINT の差 (矩形磁石、1.4e-3) | v2 は RADINT と合わない | A1・A2 の判断 |
| numba のキャッシュの自動無効化 | `pysad/nbcache.py` を手で走らせる | 自動で呼ぶ仕組みは kernel を持つ共有ファイルの編集が要るので入れていない |

## 7. 全試験・回帰

2026-11-18 に最後に 1 回だけ走らせた (`PYSAD_SAD_PACKAGES` つき、k64 が既定、Python 3.12):

- 全試験 (`py -3.12 -m pytest pysad/tests -q --ignore=test_bend_feeddown.py --ignore=test_coord.py --ignore=test_sext_radiation.py`): **1281 通過・6 skip・1 失敗** (所要 12 分 44 秒)。失敗は `test_v2_C2_spch.py::test_spac_with_trpt_acceleration`: 試験の台本の CAVI の `PHI=-Pi/2` が gs では構文エラーになる形 (MAIN の `PI` は大文字だけ。(1) 参照) で、(1) の変更で捨てられたため。`PHI=-PI/2` に直して通る (C2 の試験の 1 語の修正。この 1 件を含め `test_v2_C2_spch.py` の 35 件は全部通る)。基準 (937 通過・4 skip) は他の担当の試験が増えて総数が違う。
- Python 3.14 (`py314` の venv): `test_tclbridge.py`・`test_plot.py`・`test_gui_lang.py` 25 件、`test_v2_D2_*`・`test_spin.py`・`test_corr.py`・`test_csr.py`・`test_fam_points.py` 90 件が通る (1 skip)。全試験は 3.14 では走らせていない。
- 回帰 (`QUICK=1 run_regression.py`): **50 ベンチ、FAIL 0** (所要 140 秒)。

## 8. 直したファイル

`pysad/parser.py` (要素定義の式、一般の代入の評価の統合、0 での割り算)、`pysad/spin.py` (3.14 の複素数: `_spintrm`・`_cerfc`)、`pysad/corr.py` (v2 の Micado)、`pysad/tclbridge.py` (Tk 9.0 の固定幅フォント)、`pysad/nbcache.py` (新規)。
試験: `pysad/tests/test_v2_D2_parser.py`・`test_v2_D2_packages.py` (新規)、`test_tclbridge.py`・`test_v2_C2_script.py`・`test_v2_C2_spch.py` (各 1〜2 行)。改行だけ直したファイル: §4.1。
