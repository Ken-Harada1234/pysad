# 作業ログ: item4 (SADScript 関数・CHROMATICITY・REPEAT/UNTIL/IF・CalculateOptics・2 段目フィット) (2026-10-15〜16)

担当: pysad/sadscript.py, pysad/sadlist.py (新規), pysad/tsvd.py (新規), pysad/ffscmd.py, pysad/match6.py, pysad/lattice.py。
方針: 1 項目ずつ SAD (WSL の `gs`) で参照出力を取り、pysad の出力と突き合わせる。各項目の前後で回帰
(`QUICK=1 run_regression.py`、`pytest` 57 本、13a/16a の全実行) を取り、変化がないことを確認した。
照合用の小道具: `cmp_tags.py` (タグ行 `T001 ...` の比較。`--num ABS REL` で数値許容)。SAD の実行は
`wsl -d Ubuntu -- bash -lc 'cd ".../<repo>" && .../bin/gs file.sad < /dev/null > file.result'` (標準入力を閉じないと FFS 終了後に待ち続ける)。

回帰のベースライン (作業前): QUICK=1 で FAIL 1 (3c の EMITY 差 5.3e-16 / 許容 3e-16、既存のノイズ級)、KNOWN 3 (3c EX/EPX、14a dz)。
pytest 57 pass。13a は PASS* (TWB 4.7e-8、既知)、16a は PASS* (DA の損失周数 ±1 周が 1 個、既知)。

---------------------------------------------------------------------------------------------------

## A. SADScript の欠落関数 (bench_level24a.sad)

### SAD の実装箇所
- 組み込み (src/*.f, tfinitn.f の itfunaloc): MapThread, MapAt, Complement, DeleteCases, ScanThread, Cases, MatchQ, Extract, Intersection,
  Union (Union1), Position, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map, Tr, Depth, Level, AtomQ, EvenQ/OddQ ほか。
- Packages/functions.n: `Count[x__]:=Length[Position[x]]`、`MemberQ[l,p]:=Position[l,p,1,1]<=>{}`、`FreeQ[f,p]:=Position[f,p,{0,Infinity},1]==={}`、
  Nest/NestList/Fold/FoldList (Do ベース)、Through、Positive/Negative/NonNegative/IntegerQ、`FractionalPart[x]=x-Floor[x]`、RotateLeft/Right、
  `Union[x__]:=Union1[Join[x]]`。Packages/Misc.n: FixedPoint (既定の SameTest は NearlySameQ = $NearlySameQ[a,b,1e-12,1e-20]、src/tfsameq.f)。
- **SAD に無い関数** (SAD では未評価のまま返る。bench で確認): `Total`, `Cross`, `Norm`, `Riffle`, `StringSplit`, `Mean`, `Variance`, `Most`, `Accumulate`,
  `StringTake`, `IntegerPart`, `Eigenvalues`, `Ordering`, `Tally`, `DeleteDuplicates` など。依頼の一覧にあった Total/Cross/Norm/Riffle/StringSplit は
  実装していない (SAD で動かないスクリプトは存在しないため、pysad も SAD と同じく未評価のままにするのが「SAD の数値の再現」になる)。
  既存の pysad 拡張 (Total, Most, Accumulate, StringTake, IntegerPart, Eigenvalues) は壊さず残した (bench では「既知の拡張」として除外: T118 T191 T195 T220-222 T229 T231 U143 W075)。

### pysad の変更
- `pysad/sadlist.py` (新規): 上の関数群を SAD の挙動で実装し、`_install_builtins` の最後で同名の既存関数を上書き (`sadlist.install(I)`)。
  - レベル指定 (`Cases/DeleteCases/Position/Map/Apply/Level/Scan`): n → {1,n}、{n}、{m,n}、Infinity、負 (高さ)。走査は子 → 親の後順。
    **Position だけ頭部 (位置 `{0}`) も走査し既定のレベルは {0,Infinity}** (Count/MemberQ/FreeQ はこれの上に作られているので `Count[{1,2,3},_]` は 5、`FreeQ[f[g[x]],g]` は 0)。
  - 位置指定 (`MapAt/Extract/Delete/Insert/ReplacePart`): 数 n、`{i,j,..}` (1 つの深い位置)、`{{..},{..}}` (複数)。負の添字 OK。
    `Insert[l,x,-1]` は `{1,2,x,3}` (SAD は Mathematica と違い最後の要素の「前」)。`Delete[l,0]` は SAD 自体が異常終了する (ベンチから除外)。
  - 標準順序 (Sort/Union/Complement/Intersection/Order): 数 < 文字列 < シンボル < 複合式 < リスト、リストは長さ優先、文字列/シンボルは `a<A<aa<b<B` (小文字が先、先頭文字優先)。
  - `Take/Drop` の範囲外: `Take[{1,2,3},5]=Take[l,-5]={}`、`Take[l,{2,5}]={}`、`Drop[l,{2,5}]={1}`。2 つ以上の範囲指定 (`Take[m,2,2]`) は SAD ではエラー。
  - `Flatten[l,n,head]`、`Partition[l,n,d]` (3 引数まで)、`Transpose` (行列のみ、不揃いはエラー)、`Join` (同じ頭部のみ)、`Reverse/First/Last/Rest/Append/Prepend/Select/Map/Apply` は複合式 `f[..]` も受ける
    (`Select` は複合式でもリストを返す)。`Apply[f,l,{lo,hi}]`・`Map[f,l,lev]`・`MapAll`。
  - `Round[x,a]`、`Floor/Ceiling` (2 引数はエラー)、`Mod` (0 での割りはエラー)、`Max/Min` (入れ子・複数引数、空は -INF/INF、非数値を含めば未評価)。
    数学関数 (Sin, Sqrt, Abs …) は未定義シンボルを渡すと未評価で返す (従来は SADError)。
  - `FixedPoint` (NearlySameQ 既定、n 回上限)、`NearlySameQ`、`Fold/FoldList/Nest/NestList` (Fold は 3 引数のみ)、`Through`、`Tr[m]/Tr[m,List]`、`Depth/Level/AtomQ/IntegerQ/EvenQ/OddQ/Positive/…`、
    `Binomial`、`FractionalPart`、`Inverse` (特異・非正方は擬似逆行列)、`Det` (非正方はエラー)、`Dot` (形の不一致はエラー)、`IdentityMatrix[0]` はエラー、
    `StringReplace` (規則は先頭から先に当たったもの、同時置換、文字列 1 つのみ)、`ToExpression` (複数文、未完はエラー)、`FromCharacterCode/ToCharacterCode`、`Outer` (3 リスト以上)、`Thread` (数値・文字列はエラー)。
- `sadscript.py` の変更:
  - **暗黙の積** (`2 Pi`, `a b`, `2 (x+1)`) を実装 (SAD は Times として解釈する。従来は構文エラー)。改行をまたぐ場合は積にしない。
  - **記号式の四則演算** (`x+y`, `2 x`, `x/y`, `-x`, `x^2`) を未定義シンボルで未評価のまま作る (従来は `$Op[..]`)。SAD と同じ並べ替え (`x+y+x → x+x+y`、`1+x+3 → 4+x`) と表示 (`(2 x)`, `x+(-y)`, `(x y^-1)`, `x^2`)。
    `Cases[l,{x_,y_}->x+y]` のように規則の右辺に式を書く使い方が動く。
  - Expr の表示は `f[a,b]` (SAD はカンマ後に空白を入れない。従来は `, `)。`Print[Null]` / `ToString[Null]` は `Null` (リストや関数の引数の中では空)。
  - `Sum/Product` の複数反復子 (`Sum[i*j,{i,3},{j,2}]`)。反復子が `{i,{a,b,c}}` 形式のリストは SAD に合わせてエラー (Table/Do/Sum)。
  - パターン: リスト `{_,_}` と `f[__]` の長さ可変パターン、`Repeated` (`{_Real..}`) を `_match_one`/`_pat_info` に追加。`Condition (/;)` は SAD に演算子が無い (構文エラー) ので実装しない。
  - 組み込み関数に不正な引数を渡したときの Python 例外 (TypeError など) は SAD の `???General::wrongtype` と同じく、その文を打ち切って次の文へ (最上位のみ)。
- `lattice.py`: `TwissRow.betz` の既定を 1.0 に (SAD の 4 次元 CALC/INS では `Twiss["BZ"]` は常に 1。`Twiss["*",1]` の 23 番目)。

### 照合結果
- `bench_level24a.sad` (T/U/V/W の 4 系列、559 タグ、約 600 の Print。SAD の `.result` と `$FORM="15.10"` の桁まで文字列一致):
  `python cmp_tags.py bench_level24a.result bench_level24a.pysad T118 T191 T195 T220 T221 T222 T229 T231 U143 W075` → `tags 559 diff 0 known-extension 10`。
- SAD の癖で再現していないもの: `Level[{1,{2,{3}}},{-1}]` が `{1,3}` (SAD のバグ)、`DeleteCases[..,n]` の第 4 引数、`Inverse` の非正方・ゼロ行列 (SAD は別の疑似逆/NaN)、
  `Insert[l,x,0]`、`Fold[Plus,0,{{1,2},{3}}]` (SAD は `{3}+{1,2}` を未評価で返す)、`N[Sqrt[#]]&` の FixedPoint (SAD は "Too deep")。
- 回帰 (QUICK=1 / 13a / 16a / pytest 57): 作業前と同じ (FAIL 1 = 3c EMITY は既存、KNOWN 同じ、13a/16a の出力は警告の行番号以外 diff なし)。

---------------------------------------------------------------------------------------------------

## B. CHROMATICITY コマンド (bench_level24b.sad, bench_level24b2.sad)

- SAD: src/tfchro.f (51 行)。tffsa.f:609 から `tfchro(latt, twiss(AX,BX,PSI,EX..))` (itwisso の i=0 = 直前の CALC の行)。
  四極・六極・MULT の寄与を 1 要素ずつ `1x,a12,1x,2f17.3` で表示し、最後に `TOTAL`。BEND は数えない。
  - 四極/MULT の K1: `ξx = -(K1·βx + (1+αx²)/βx·L + αx(i+1) - αx(i))/2`、`ξy = -(-K1·βy + (1+αy²)/βy·L + Δαy)/2`
  - 六極 (MULT の K2 も): `ξx = (βx(i)ηx(i)+βx(i+1)ηx(i+1))·K2/2`、`ξy = -(βy(i)ηx(i)+βy(i+1)ηx(i+1))·K2/2` (垂直の分散は使わない)
  - 六極単独の ξy は Fortran の代入 (`xiy=-(..)*v*.5`) なので K2=0 のとき `-0.000` と出る (MULT は足し込みなので出ない)。
- pysad: `ffscmd.py` に `chro_cmd` を追加 (`CHRO`/`CHROMATICITY` の略記を `_abbrev` で受理。従来は「未対応」と警告して無視)。
- 照合: `bench_level24b.sad` (pysad/tests/kek-ls_v03_12.lat 連結、`CELL; CALC; CHROMATICITY;`): 要素 960 行 + ヘッダ + TOTAL の 962 行が**文字列で全行一致** (TOTAL -85.752 / 751.974)。
  `bench_level24b2.sad` (小リング: QUAD/SEXT/MULT(K1,K2)/BEND、DP0=0.002 での 2 回目も): 244 行すべて一致 (MULT 経路の検証)。

---------------------------------------------------------------------------------------------------

## C. REPEAT/UNTIL と IF/ELSEIF/ELSE/ENDIF (bench_level24c.sad)

- SAD: tffsa.f:255-290 (REPEAT: `levelr`, `nrpt(levelr)`, `irptp(levelr)=ipoint` で読み位置を覚え、UNTIL で `cssetp(irptp)` により巻き戻す。
  `REPEAT n` の n が無ければ 65535。`UNTIL 条件` は条件が真なら抜ける、条件が無ければ回数だけ。UNTIL/REPEAT の多重は `levelr` のスタック)。
  tfif.f (IF: 条件が偽なら同じ深さの ELSE/ELSEIF まで読み飛ばす。ELSEIF は ELSE として扱って再評価。ELSE に来たら ENDIF まで読み飛ばす。`iflevel` で入れ子を数える)。
  tfif.f の delim 判定: 語の直後が `[ * / ) } ] > < = . & ^ |` なら関数呼び出しや式 (`If[..]`, `If==1`) なので制御語ではない。
  エラー表示: ` ???-FFS-Error-?UNTIL without REPEAT. ` / `?ELSE without IF` / `?ENDIF without IF` / `?ELSEIF without IF` / `?Missing condition for IF` (IF 後は iflevel=0)。
- pysad: `sadscript.Interp.run` の最上位 FFS の文列を、添字つきの while ループにして `_control_word / _control / _if_cond / _skip_if` を追加。
  条件式は `_eval_cond` (FFS の式は空白で終わる。ただし演算子の前後の空白は式の続き: `IF KSTEP<= ISTEP`)。`IF 条件 コマンド` と同じ行に続くコマンドも実行する。
  FFS[..] の中や SADScript の関数内では働かない (SAD も読み位置の巻き戻しなので文字列の FFS[] では同様に使えない)。
- 照合: `bench_level24c.sad` (REPEAT n と UNTIL 条件で K1 を変えながら CALC、条件なし UNTIL、入れ子、IF/ELSEIF/ELSE の分岐、入れ子 IF、REPEAT の中の IF、IF の中の REPEAT、
  1 行 IF、対応のない UNTIL/ELSE/ENDIF): 19 タグすべて一致 (チューンは 10 桁)。エラーメッセージも probe7.sad で SAD と同形式を確認。

---------------------------------------------------------------------------------------------------

## D. CalculateOptics (bench_level24d.sad, bench_level24d2.sad)

- SAD: src/tfoptics.f → qcell1 (qcell.f)。`CalculateOptics[lb, le, inicond(28 個), cell, idim (,geocal)]` が
  `{ {lb..le の各位置の 28 個の Twiss}, {stabx,staby,stabz}, {tracex,tracey,tracez} }` を返し、大域の Twiss 配列も書き換える (呼出し後の `Twiss[]` に反映)。
  - cell≠0: 区間を周期系として閉軌道 + 周期 Twiss を求める (inicond は閉軌道探索の初期値 [DX..DDP] と位相 PSIX/PSIY のみ)。cell=0: inicond を位置 lb の入口条件として伝搬。
  - tracex/tracey = 区間の転送行列の対角和を位置 lb と le の Twiss から作った `x11+x22` (qcell.f:355-390)。cell では |trace|>2 のとき INTRES/HALFRES が真なら不安定、偽なら 4/trace。tracez=0。
  - 周期解が見つからないとき ` *****qcod---> Overflow & closed orbit not found            1` を出し、stab={0,0,1}・trace={0,0,0}・Twiss は元のまま。
  - **stabz は SAD では未初期化変数** (`ffs_stat optstat` が intent(out) の局所変数): 呼出しごとに 0/1 が不定 (部分区間 cell=0 では 0 になった)。idim=3 の stab/trace も未初期化 (`{0,0,0} {1,0,0}`)。pysad は 1 を返す。bench では stab は先頭 2 要素だけ、idim=3 は stab/trace を除いて比較。
- pysad: `ffscmd.fn_calculate_optics` (従来は「not implemented」)。cell≠0 は区間の部分ラティスの `twiss_table(delta, cod)`、cell=0 は `ins_twiss` を使う。idim=3 は全周 (1, nlat) のみ (CALC6D の経路)。
  - `lattice._one_turn_cod_matrix` の軌道ありのドリフトに **経路長 z の更新 `cod[4] += dpz*ale`** を追加 (qtwiss.f:224)。従来は x, y だけ動かしていたので、
    INS で初期軌道がある場合の `DZ` (Twiss["DZ"]) が SAD と 1.9 倍ずれていた (bench_level24d の D30 で発見)。
- 照合: `bench_level24d.sad` (小リング 181 点: cell=1 全周、cell=0 で BX/AX を変えた伝搬、区間 3..9、名前で範囲指定、R1..R4 つきの初期条件、初期軌道 (DX,DPX) つき、分散つき、DP0≠0、
  不安定 (QF を強く)): 29 タグすべて数値一致 (`cmp_tags.py --num 1e-8 1e-7`。R1.. など 1e-16 級の雑音は SAD の丸め誤差)。
  `bench_level24d2.sad` (ESRF 3564 点、idim=3 の 6 次元周期解と idim=2): 6 タグ一致 (6 次元の Twiss 28 成分 3 点とも 1e-8 以内)。

---------------------------------------------------------------------------------------------------

## E. 2 段目以降のフィットの発散 (<work-pf2>/matching_v30_12_r ほか)

### 症状と調べ方
- 症状 (work_pf2_survey.md): 1 段目 (CELL、NCELL の 2 変数) は SAD と NEWTON 残差が 5 桁一致するが、2 段目 (`USE N2; ins calc; Fit MKDM EX 0 …; free Q*A*; go`、INS で 9 条件 × 11 変数) は
  初期残差が SAD の 1e6〜1e13 倍 (1.9e23 vs 4.0e6) で `optics calculation failed`。現行コードでは失敗はしないが、残差の列が SAD と違い (初回 1.5e5 vs SAD 3.979E+06、最終 0.755 vs 0.832)、別の解に着いていた。
- SAD の内部を見るため、スクリプトに **FitValue / FitWeight / VariableRange / VariableWeight のフック** を定義して Print させ、SAD の各反復の変数値・条件ごとの現在値・重みを取り出した
  (`FitValue[loc_,fn_,{id_,dp_},goal_,now_]:=(Print[..];goal)`、VariableRange は `Null` を返して dlim に影響しないようにする)。pysad でも同じ出力を取って比較した (<work-pf2>/e2.sad)。
- 比較の結果、**初期の条件値 (9 条件) と、SAD の 1 歩目の変数値での条件値は pysad と 10 桁一致** (光学の計算は正しい)。違ったのは (1) 重み、(2) 1 歩目の Newton の解。

### 原因と修正 (4 つ)
1. **USE が EMITX.. を毎回マーク値から読み直していた** (pysad `use()`): SAD は `{EMITX,EMITY,EMITZ,SIGZ,SIGE}=LINE[{..},1]; DP0=LINE["DDP",1]` を **FFS に入った最初の 1 回だけ**
   実行する (tffsa.f:135-139。USE の `chguse=.true.` で 101 に戻る経路ではこの文を飛ばす: `if(.not. chguse)`)。1 段目の `emit;` が閉軌道を見つけられず EMITX=EMITY=NaN になると、
   2 段目の `USE N2` 後もこの NaN が残り、tffscalc の重みは gfortran の MAX/MIN が NaN を無視して emminv 由来になる (EX/EPX の重みが 1.686e6 など)。pysad は 0 に戻していたので重みが 50.5 倍小さかった。
   → `use()` を最初の 1 回だけマークから読むようにした。FitWeight の値 (`1686083.081048` など) が SAD と完全に一致。
2. **tsvd (SAD の SVD ソルバ) が numpy の硬い打ち切りだった** (match6 `_tsolva`): J の特異値は 1.6e7 〜 1.5 (比 1e-7〜1) で、SAD の `tsvd(.., eps=TSTOL/fact=1e-6)` は
   2 重対角化 (Householder/Givens) の後、対角要素の最小値が `anorm*eps` を超えれば **特異値を一つも落とさずに** 2 重対角系を直接解く (quick path)。numpy 版 (特異値 ≤ 1e-6·最大 を落とす最小ノルム解) とは
   この問題で解が大きく違い、SAD の 1 歩目は条件をほぼ完全に満たす (残差 5e9 → 0.53) 大きな歩だった。
   → `pysad/tsvd.py` (新規) に tsvd の svd=.false. の経路 (quick path と暗黙 QR + 滑らかな減衰 w=(s·v/anorm²)²) を直訳。良条件の行列では numpy の pinv と 1e-14 で一致 (n<m, n=m, n>m)。
   pf2 の J で SAD の 1 歩目の K1 変化 (10 桁) を再現する (`pysad/tests/test_item4.py` と `data_tsvd_sad.npz`)。環境変数 `PYSAD_TSVD_NUMPY=1` で従来の numpy 版に戻せる。
3. **INS のフィットが数値微分だった** (`Matcher._analytic_unsupported` が「INS は未検証」で全体を中心差分にしていた): SAD は INS でも解析微分 (tffsqu→qdcell→qdtwis、cell=False で周期解の補正なし) を使う
   (数値微分は `FFS$NumericalDerivative` のときだけ。SAD の INS 数値微分は `qtwiss Overflow (fb)` で壊れるので比較にも使えない)。
   match6d の解析微分は cell=False に対応していたので制限を外した (`PYSAD_INS_NUMERIC=1` で数値微分に戻せる)。tune_sv の 2 段目は 100 反復以上 SAD と同じ列になった (数値微分では 55 行目から 1e-5 ずれていた)。
4. フィットの条件フック: SAD は `FitValue[component, function, {id_,dp_}, goal_, now_]` (2 点条件は `FitValue[comp1, comp2, fn, {id,dp}, goal, now1, now2]`)、`FitWeight[comp, fn, {id_,dp_}, w_]` を呼ぶ。
   pysad は dp を裸の実数で渡していたため `{_,dp_}` のパターンに一致せずフックが呼ばれていなかった → `[0.0, dp]` のリスト (id は通常 0) で渡す。
   **VariableRange[element, keyword, v_]** と **VariableWeight[element, keyword, v_]** (tffsvarfun) を実装 (従来は未実装で無視。ユーザースクリプトの 81 ファイルが使う)。VariableRange が実数を返せば dlim が真になり、
   反復の終わりに範囲に収め直す。0/偽なら変更前の値に戻す。`{min,max}` なら範囲を狭める。
- 副産物: `fmt_g` の 1PG11.4 で整数に近い値を `1080.` と出す (従来は `1080`)、`CALC` で光学が計算できないとき Python の RuntimeError で全体が止まっていたのを、`*****qcod---> Overflow ...` を出してエラーの文として続行するように。

### 結果 (SAD の `*_r.sad.out` との比較: <work-pf2>/*.E3*.out、cmpfit.py)
| スクリプト | 結果 |
|---|---|
| matching_v30_12_r (NCELL → N2 → RING) | NEWTON/Residual 13 行が SAD と一致 (`6.1874E+04` と `6.1875E+04` の末尾 1 桁の丸めのみ)。最終 `Print` の νx=12.0851271700, νy=6.7307624317 一致。238 s |
| matching_v31_23_r | 1 段目 Matched 2.6724E-11、2 段目 1.3019E-12 (SAD 1.3022E-12)、最終 Residual 2.1319E-08 が一致 (従来は 2 段目が別解)。318 s |
| tune_sv_v30_24_1 | 2 段目 (12 変数) は **SAD と反復表が全行一致** (10/16 に解消: 113 反復目の DESCEND の分岐は、モード切替後に重みを計算し直す `wcal=chgmod` の欠落が原因で、丸めではなかった。docs/history/work_descend.md)。最終 Residual 1.8569E-02。 続く `USE RING` の `cell calc` は SAD でも不安定 (Vertical Unstable) だが pysad は例外で全体が止まっていた → 続行するよう修正 |
| chrom_r (FITP、FREE S*) | 173 s で完走、Residual 7.9187E-08 (SAD の 1 反復 7.9e-8 と同じ) |
| matching02_r (NCELL, 変数 1 個の解なし) | 不変 (最終 Residual 4.069 vs SAD 4.072。`fact` の列が違う別件: work_integration.md §2) |
- 回帰: QUICK=1 + 13a + 16a (全 21 ベンチ) FAIL 1 (3c の EMITY、既存)、13a の TWB 4.74e-8 も不変 (tsvd は 13a の良条件の J では numpy と同じ解)。pytest 63 pass (新規 test_item4.py 6 本)。

---------------------------------------------------------------------------------------------------

## F. Element / SetElement / FLAG / ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine / GetMAIN / USE BeamLine (bench_level24f.sad)

実装: `pysad/sadelem.py` (新規)、キーワード表 `pysad/sadkw.py` (SAD で型ごとに `SetElement["Zname",type]` を実行して得た順序)。SAD の Packages/element.n と beamline.n を読み、bench_level24f.sad
(164 タグ。Print/WriteBeamLine の出力は `L001 begin … L001 end` のブロック比較) で照合 (known 1 個 = `SetElement["ORG"]` の MARK の値)。

- **Element[key, spec]** (tfelement): 相手は「現在のビームライン中のファミリーの一覧 (重複なし、初出順)」。`Element["LENGTH"]` はファミリー数。spec 省略 = 全ファミリーのリスト、
  文字列 = 名前またはパターン (SAD の tmatch: `* % {abc} {^abc} |`。**`?` は普通の文字** — pysad の `_wild_match` が `?` を 1 文字ワイルドカードにしていたのを修正)、数 = ファミリー番号
  (負は末尾から、範囲外はエラー)、リスト = 要素ごとの結果のリスト。**単一のマッチはスカラー、`{"QF"}` や `{2}` のようなリスト指定はリスト、マッチ無しは `{}`**。
  key: `NAME POSITION VALUE KEYWORDS DEFAULT TYPE TYPENAME COMPONENT` と要素のキーワード (型に無ければ 0)、key のリストは key ごとに写像。`EXPAND` は Null、`UPDATE` は 0 のリスト。
  `Element[k,spec]=v` は spec のすべてのファミリーに代入 (リスト spec の代入は `{v,v}` を返す)、`VALUE` は既定キーワード。`+= -= *= /=` も使える (既存)。
  Element は MAIN の定義値を変えない (`USE` し直すと元に戻る)。
- **SetElement[name, type, opts]** (SetElement$): 新規なら MAIN に作り、既存なら指定のキーワードを変更し、`{name, TYPENAME, {"L"->v, "K1"->v, …}}` を返す (キーワードは型ごとの全一覧)。
  未定義で型なしは `{"XX",""}`、型が違えば `FFS::equaltype` エラー。**MAIN の定義値 (env) と現在のビームラインの要素の両方を更新**。
  SAD の癖: シンボルで渡すと戻り値の名前が `""` (要素は変更される)、クラッシュする呼び方がある (ファイルへリダイレクトした出力では `SetElement[QF,QUAD,K1->0.3]` で SIGSEGV。参照は pty 経由 `script -qec` で取る)。
- **FLAG[name]**: 1/0。リストはリスト。未定義は `FFS::undefflag` エラー (**大文字小文字を区別**: `FLAG["cell"]` はエラー)。
- **ExtractBeamLine[line]**: 展開したビームラインを `BeamLine[ORG,QF,DS,…]` で返す (引数なし = 現在のライン、逆向きの要素は `(-X)`)。**BeamLine[...]**: n*x の繰り返し展開、`-x` は `(-x)`、
  BeamLine の入れ子は平坦化、LINE の名前は展開しない (`BeamLine[2*HALF]` は `BeamLine[HALF,HALF]`)。`Join/Reverse/Length/Part/Map` が使え、**`-BeamLine[..]` は順序を逆にして各要素を反転、`n*BeamLine[..]` は繰り返し**
  (Packages/beamline.n の UpValue。sadscript の `_symneg/_symb` に実装)。
- **PrintBeamLine / WriteBeamLine[f, b.., Format->"MAIN"|"FFS", Name->{..}]**: Packages/beamline.n の WriteBeamLine1 の直訳。`Fill$String[s,9]` で 9 文字の倍数に詰める出力 (`B D1` のように短い名前の次が詰まる癖もそのまま)、
  col>69 で改行。**SAD は各 WriteBeamLine1 の先頭にデバッグ出力 `WriteBeamLine {{Name->None,Format->"MAIN"},"MAIN",None}` を標準出力に出す** (Packages に残っているもの) ので同じ行を出す。
  複数の BeamLine を渡すと `Complement` で標準順序に並べ替えられ重複が消える (`WriteBeamLine[6,BeamLine[QF,DS,QD],BeamLine[QD,DS,QF],..,Name->{"A1","A2"}]` は A1=(QD DS QF) になる)。
- **USE / VISIT <式>**: `USE BeamLine[ORG,QF,DS,QD]` や `USE Join[a,-a]` を評価して BeamLine なら新しい LINE (名前 `L000001$` …、成功のたびに通し番号) を MAIN に作って使う。
  先頭が MARK でなければ `???FFS::firstmark: The 1st component must be MARK:` (+ 要素名の頭文字と `^`) と ` ???-FFS-Error-Missing beamline in USE. `。**VISIT / BYE** は直前のラインを覚えて戻る (従来 BYE は何もしなかった。FFS の他の状態は復元しない)。
- **GetMAIN[file]**: MAIN 形式のファイルを現在の環境に読む (要素・LINE が増える。現在のラインは変わらない)。
- 副産物の修正: **L を変えても K1.. (積分強度) は保たれる** (`Element["L","QF"]=0.5` で `K1` は変わらない: SAD の仕様)、`Interp.part` が複合式 `f[a,b][[2]]` を受ける。
- 注意 (SAD 側の挙動、ベンチを作るときに踏んだもの): 不安定なラインの `CELL; CALC;` は SAD が (x の) チューンらしい値を返すが pysad は光学を計算できずエラーにする (H005 を除外)。
  LINE の名前が FFS のフラグ名 (`CELL` など) と同じだと `BeamLine[ORG,CELL,CELL]` が SAD でエラー。FLAG の未定義エラーの直後に値が元に戻って見えたのは、テストに混ざっていた `USE RING;` のせいで SAD の挙動ではない。
- 未対応: `Element["PROFILE"…]`、`Element[..,Saved->True]` の保存バッファ (オプションは無視)、`LINE["DK",..]=` との DK の連動、Print の入れ子リストの折り返し (SAD は長い入れ子リストを要素の境界で改行するが
  pysad の `wrap` は区切りのカンマだけ見る。行の上限も SAD は pw-1 文字。PageWidth が十分大きければ一致)、3 引数の `Element$`。

---------------------------------------------------------------------------------------------------

## 変更ファイルの一覧
- 新規: `pysad/sadlist.py`、`pysad/tsvd.py`、`pysad/sadelem.py`、`pysad/sadkw.py`、`pysad/tests/test_item4.py` (+ `data_tsvd_sad.npz`)、`cmp_tags.py`、
  ベンチ `bench_level24a.sad / 24b.sad / 24b2.sad / 24c.sad / 24d.sad / 24d2.sad / 24f.sad` (+ `.result` = SAD の出力、`.pysad` = pysad の出力)。
- 変更: `pysad/sadscript.py` (暗黙の積、記号式、制御構文、パターン、Table/Sum、Print Null、fmt_g、BeamLine の演算)、`pysad/ffscmd.py` (CHROMATICITY、CalculateOptics、USE/VISIT/BYE、フック、`use()`、`_wild_match`、L の保存)、
  `pysad/match6.py` (tsvd、INS 解析微分、VariableRange/VariableWeight)、`pysad/lattice.py` (betz、ドリフトの z)、`run_regression.py` (24a〜24f を追加)。

## 最終回帰 (全作業後, 2026-10-16)
`PYTHONIOENCODING=utf-8 python run_regression.py` 全 28 ベンチ (既存 21 + 新規 24a/24b/24b2/24c/24d/24d2/24f): FAIL 1 (3c の Emittance[] EMITY 5.29e-16 / 許容 3e-16 — 作業前から同じ)、
KNOWN は作業前と同じ (3c EX/EPX、13a TWB、14a dz、16a DA の損失周数 ±1)。13a/14a/16a の出力は作業前と差なし、13b は 1e-14 級の雑音の桁のみ。pytest 63 pass。
参照用ファイルの整理: 調査用の probe スクリプトと SAD の出力は `probes_item4/` に移した (`lat_extra.sad` は bench_level24f が GetMAIN で読むのでルートに残す)。
