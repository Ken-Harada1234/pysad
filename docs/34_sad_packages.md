# SAD の Packages/*.n を SADScript のまま使う (軌道補正・最適化・数学パッケージ)

最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。経緯は `history/` にある)

pysad の SADScript インタプリタが、SAD の `Packages/*.n` を最初に使われたときに読む ([`pysad/pkgload.py`](../pysad/pkgload.py))。Packages ディレクトリは
環境変数 `PYSAD_SAD_PACKAGES`、なければ pysad の隣の `SAD-master/Packages`。見つからなければ何も登録しない (組み込みだけで動く)。

| パッケージ | 提供する関数 | 状況 |
|---|---|---|
| correction.n (+ bump.n、corutility.n) | CorrectOrbit、CorrectLocal、MakeBump、Bump、Monitor、Steer、CorStat、RmsReport、NearestMonitors、ReadMonitorData … | gs と一致 (下記) |
| optimize.n | OptimizeOptics、Amoeba、FindClosedOrbit (DownhillSimplex は Python で書いた `miscfn.py`) | OptimizeOptics が gs と 1e-10 で一致 |
| rational.n、tableform.n、Couple.n、emit.n の Beamstrahlung 系 | Rational、TableForm、ElementValues (Couple.n の SetCoupledElements)、BeamstrahlungUpsilon … | 動く (gs と一致) |
| D.n | Spline (K64 の `Spline$` は miscfn.py)、記号微分 `D`、`Derivative`、`ReduceExpression`、`Integrate`、Spline、特殊関数 (HermiteH など) | `D` を含め動く (59 式が gs と文字列まで一致。下記) |
| LDSolve.n、Optics.n | LDSolve、Optics の TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo | gs と一致 (下記「LDSolve・Optics.n」「v2.1.0 の確認」)。PeriodicR も gs と 3e-16 |

組み込み (Python) で持つ関数 (LinearSolve、Binomial、Spline など) は、パッケージの定義で上書きしない (`pkgload._protect_builtins`)。

## UpValue を定義するパッケージ (MakeBump・Rational・ElementValues)
- `bump.n` の `(a_=Bump[…])^:=…`、`rational.n` の `Plus[Rational[…],…]^:=…`、`Canvas.n` の `(Canvas$Widget=Widget[c_,_Canvas])^:=…` は、組み込みの `Set`・`Plus`・`Times`・`Power`・`Minus`・`Divide`・`Subtract` の定義として持ち込まれる。
  `pkgload._protect_builtins` は、組み込みと同名の定義を捨てるが、この UpValue 由来の定義 (Set・SetDelayed・Plus・Times・Power・Minus・Divide・Subtract) だけは残す。
  `b=Bump[…]` のように UpValue の名前 (Bump) が最初に使われる文でも読み込まれるよう、`(a_=Name[…])^:=` の形から名前を登録する (`_UPDEF`)。`Canvas$Widget=c1` のように左辺の記号が読み込みの引き金のときも読む。
  `MakeBump` は gs と同じバンプ (ybump=.0010000141、使うステアラ 12 本) になる。
- `Rational`: `1/3+1/6` は gs と同じ `Rational[3,6]` (`ReducedRational` が p^2>b で打ち切るため。gs も同じ)。`a-b`・`a/b` は Plus・Times の規則へ回す (`Subtract`・`Divide` は Rational を含むとき `a+(-b)`・`a*(1/b)`)。
- `ElementValues = {"K1"["QD"] -> -0.8*"K1"["QF"]}` は、次の `CALC` の前に `SetCoupledElements` (Couple.n を SADScript のまま実行) で要素値に反映する。`USE` で消える (gs と同じ)。
  `"K1"["QD"]` のような文字列を頭部に持つ式は評価せず式のまま (gs と同じ)。QD の K1 は -0.56、NX=6.788257 (gs と一致)。
  `VISIT` での退避と `BYE` での復元は v2.1.0 で実装した (docs/04)。**未対応**: フィットの変数との連動 (`CouplingMatrix`・`CoupledKeys` を使う Fortran 側 `tffssetupcouple` は呼ばない)。
- Couple.n や `With[{f}, f&]` を動かすため、`With` は保持された式 (`Hold[y]`) と純関数の本体にも値を埋め込む。`Hold[…]` の中の `#` は純関数の置換の対象にしない。`Hold[…]/.規則` は規則の右辺 (:>) を評価せず埋め込む (`Hold` を別の頭部に替えるときは評価)。
  `Hold[…][[1]]` は取り出した式を評価し直す。`HoldAll` の引数位置の `Evaluate[x]` は評価する。`ReplaceAll` は純関数の本体にも働く。
- `DownhillSimplex` の選択肢は、規則のリストを 1 つの引数で受けても読む (`OptimizeOptics` は `Join[{opt},{VariableRange->limit}]` を渡す)。`Output->0` で途中経過を出さず、`MaxIteration` も効く。

## 組み込み (K64 の Fortran) の翻訳 (`pysad/corr.py`)
- **ResponseMatrix** (pgrmat.f): COD 応答行列 (リングは周期解、ライン、`X Y PX PY MX MY MPX MPY`、`Tune`)。
- **LinearSolveConditioned** (pgsolvcond.f、pmicad.f): 最小二乗 (tsvd)、条件付き (psvdcond)、Micado (ベストコレクタ法)。
- **Master[]** (pgmast.f)、**FFS["VAR"]** (tfvars.f): FREE した変数の表。

gs との比較 (`tests/test_corr.py`: 12 セル、24 BPM、24 ステアラ、DX/DY の誤差つき QUAD):
応答行列 `ResponseMatrixG` の行、CorrectOrbit の補正量 (Y、X、Micado 6 本) は gs と全桁一致 (補正量 1e-7 相対)。OptimizeOptics の結果 (K1、チューン) も一致。
なお gs でも pysad でも、`CorrectOrbit` が `LINE["K0"]` に足した値は `CALC NOEXP;` では効かない (SAD は "Element values are not expanded" と表示)。

## LDSolve・Optics.n (2026-11-16 夜に gs と比較)
- Optics.n の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo`: gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し、pysad は 0)。
- LDSolve.n: pysad は `DiagonalMatrix` が複素数の列を受けず失敗していた (`float() argument must be … not 'complex'`) ので直した。あわせて行列の演算 (Dot・Inverse) の結果は、gs が要素ごとに虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻すようにした (以前は全要素が実のときだけ)。
  直したあと、LDSolve の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (試験は [`test_evalfixes.py`](../pysad/tests/test_evalfixes.py) の 2 件)。
- `PeriodicR` は v2.1.0 で比較した (下の「v2.1.0 の確認」。`PeriodicOptics` を先に呼ぶ必要がある)。

## パーサ・評価器の修正 (パッケージを動かすために gs に合わせたもの)
- パターン: Plus・Times は Flat に平らにして Orderless に照合、`f_[y__]` より `p_Plus` を先に、同名の変数を繰り返す定義 (`f[x_,x_]`) を `f[_Symbol,_]` より先に。`Plus[x_]` は評価しない。`{0 ..}` の Repeated。`//.` は不動点まで、`Unevaluated$[rules]` を規則にする、`f/@(a+b)` は評価し直す、`f[Unevaluated$[g[y,x]]]` は呼び出し側の y x の値を埋め込む (評価しない)。`-x` は `Times[-1,x]`、`x/2` は `(.5 x)`、`(2 x)^2` は `(4 x^2)`。
- 前置の `[x]` は括弧 (`Map[f,[g[x]]]`)、`f[x]g[y]`・`(a)(b)`・`2#1`・`x^2 y` は積 (D.n)。
- `Part[l,{i,j},k]` は要素ごとに残りの添字を適用。`Head` (Null と記号は Symbol、複素数は Complex、Function、Rule)。`VectorQ[x,test]`・`MatrixQ[x,test]`。
- `With[{g:=e/2}, f[e_]:=g*3]` は定義の中へ式を埋め込む。
- `CALC NOEXP`・`CALC EXP` を受理 (pysad は元から要素値を展開し直さないので `CALC` と同じ)。
- 数の既定の書式 (`autos`) を gdtoa に合わせた: `1e6` → `1e6`、`1e10` → `1e+10`、`123456789012345678.` → `123456789012345680`。`$FORM=" "` は既定と同じ。
- 物理定数 (`ElectronCharge`、`ElectronMass` = 510998.95 eV、`ElectronRadius`、`PlanckConstant` …) をスクリプトから参照できる。

## v2.1.0 の確認 (小さな台本を gs と行ごとに突き合わせた。k64 モード)

数値は 1e-12 を超える差だけを「差」とした。確認の台本 [`example/v2_next/B2_packages_compare.sad`](../example/v2_next/B2_packages_compare.sad)・`D2_periodicr.sad`・`D2_correctlocal.sad`・`D2_micado.sad`。試験 [`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py)・[`test_v2_D2_packages.py`](../pysad/tests/test_v2_D2_packages.py)。

| ファイル | 比べた関数 | 一致の状態 | 直したもの・残った差 |
|---|---|---|---|
| `functions.n` | `Log2$`、`InverseErfc`、`Names`、`RotateRight`・`RotateLeft`、`Union` (`SameTest`)、`ListQ`、`MatrixQ`、`Diagonal`、`Count`、`MemberQ`・`FreeQ` (レベルつき)、`Nest`・`NestList`・`Fold`・`FoldList`、`Through`、`FractionalPart`、`Positive` 系、`IntegerQ`、`Subtract`・`Divide`、`Xor`、`Append`・`Prepend`、`HeldPart`、`Evaluate`、`StringQ`・`Characters`・`StringDrop`・`StringInsert`、`Timing`、`For` | `Log2$`・`Union` (`SameTest`)・`Names` 以外は一致 | **`Log2$`**: 組み込みが無く未評価だった → 追加 (gs は `Log[x]*(1/Log[b])`)。**`Union` の `SameTest`**: gs は並べ替えも比較もせず結合した先頭の要素だけを返す (K64 の誤り。B121) → k64 は gs と同じ、v2 は真値。**`Names`**: gs は大域に作った記号だけ → k64 は gs と同じ、v2 は組み込みも含める。残った差: `Through[{f,g}[1,2]]` はどちらもエラーだが文言が違う。`FixedPoint[Cos,1.,3]` は最後の桁 (libm) |
| `Misc.n` | `NearlySameQ`、`FixedPoint`、`FixedPointList`、`FunctionDefinition` | `FunctionDefinition` 以外は一致 | **`FunctionDefinition`**: HoldAll でなく定義も返さなかった → `{Hold[左辺],Hold[右辺],Hold[引数を代入した右辺]}` (組み込みや定義が無ければ Null)。パターンの表示・`x//y` の表示・代入後の文字列どうしの連結を gs に合わせた |
| `element.n` | `Element` の読み出し (LENGTH・NAME・POSITION・COMPONENT・DEFAULT・VALUE・TYPE・K1・L・ANGLE・GEO・DK・EXPAND)、代入 (`=`・`+=`・`-=`・`*=`・`/=`・リスト形・`VALUE`)、`PROFILE`、`SetElement` | 読み出しは一致 | **`Element[..]+=`・`LINE[..]+=`**: 複合代入が `Element`・`LINE` のフックを通らず、利用者の規則を作って以後の代入が効かなかった pysad の不具合 → 直した。**`SetElement` の戻り値**は gs は MAIN の定義の値 → 合わせた。残った差: `PROFILE` の代入は gs も `General::invset` で pysad は別の文言のエラー、MARK の Twiss のキー (AX など) は CALC のあとの `SetElement` で gs は MAIN の 0、pysad は CALC が書き戻した値、`Element["K1","QF"]=0.9` の読み戻しが 1 ulp (0.8999999999999999) |
| `beamline.n` | `BeamLine` の演算、`ExtractBeamLine`、`WriteBeamLine`・`PrintBeamLine` (FFS 形式・MAIN 形式・`Name` のリスト)、`Fill$String`、`Remove$Paren` | `BeamLine`・`ExtractBeamLine`・`WriteBeamLine` は全行一致 | **`WriteBeamLine` の `Name`**: gs は `Name->"L2"` (文字列) で `Insufficient number of Names` (K64 の誤り。`{"L2"}` と書く。B122) → k64 は gs と同じ、v2 は文字列も受ける。`Fill$String`・`Remove$Paren` は組み込みが無かった → 追加 |
| `emit.n` | `Emittance[]` (既定、`OneTurnInformation`、`Matrix`、`Orbit`、`Emittance->False`、`Region`、`SaveEMIT`) の戻りの規則の名前・並び・値、`SymplecticJ`・`SymplecticInverse`、`TouschekF`、`TouschekLifetime` の未定義のメッセージ、`Beamstrahlung*` 4 つ | 値は 1e-12 以内 (3 GeV、周長 41.8 m 余りの FODO 8 セル。Emittances の y は gs が 2e-38、pysad は 0 の雑音) | **`Emittance->False`**: pysad は無視していた → gs と同じく DampingRate 以降の規則を返さない。`Touschek` のメッセージの形を `???Touschek::undef:  …` に合わせた。残った差: `SymplecticInverse` の `-0` の表示 |
| `dapert.n` | `DynamicApertureSurvey` (`Output`・`Axes`・`ExpandElementValues`、不正な `Axes`)、`AccelerateParticles`、`TrackParticles` | `AccelerateParticles` は 1e-14、`TrackParticles` も一致、`DynamicApertureSurvey` の Score は一致 | **`TURNS` の既定**: gs は −1、pysad は 1 → −1 に。**残った差**: `DAPWIDTH` の既定 (gs は 7 で x の振幅を上から 7 個ずつ追跡し、続けて生き残れば残りは追跡せず `*` で表示。pysad は 51 で全点を追跡。`Maximum number of particles` が gs 42、pysad 306)、`Axes->"XZ"` で z の範囲が 3 点のとき gs はエラー、pysad は受ける、不正な `Axes` のメッセージ (`DALib::wrongmode`) は未比較 |
| `Optics.n` の `PeriodicR` | 4x4 の一周行列 6 通り (結合あり・なし・半整数に近い・NaN になる 2 通り) | R・A・B の 12 成分が 3e-16 以内 | `PeriodicOptics` を先に呼ぶ必要がある (先に呼ぶと未評価。gs も同じ) |
| `correction.n` の `CorrectLocal` | 12 セル FODO (3 GeV、周長 52.8 m、BPM 24 台)、領域 {BA3, BA8}・垂直 | 補正量 12 個と 24 BPM の DY が 1e-14 以内。領域の外の軌道は 1e-11 以内で変わらない (局所バンプが閉じる) | — |
| `correction.n` の `CorrectOrbit` の Micado | 24 セル FODO (BPM 48 台、垂直ステアラ 24 台)、`Solver->{"Micado",n}` | k64 は gs と 1e-6 以内で一致 | **K64 の誤り (B129)**: `pmicad.f` が `tsvd` に nx を渡さず、全補正量の最小二乗の先頭 nx 個を使う (補正前 9.41e-5 m に対し n=1 で 1.46e-4 に悪化、n=8 で 1.75e-4)。v2 は選んだ nx 本だけを解く (n = 1・2・3・4・6・8・12・24 で 4.74e-5・4.43e-5・4.37e-5・4.31e-5・4.11e-5・3.97e-5・3.51e-5・2.06e-5 m と単調に下がる。最小二乗と同じ 2.06e-5)。`Condition` つき (`CorrectLocal`) の Micado は選んだ列だけで `psvdcond` を解く元の分岐で K64 でも正しい。未確認: 水平・2 平面同時の Micado、`Bump`・`SVDmode->1` つきの Micado |

### `dapert.n` ほか gs との主な差の一覧

| 項目 | gs | pysad (v2.1.0) |
|---|---|---|
| `DynamicApertureSurvey` の `TURNS` の既定 | −1 | −1 (以前は 1) |
| `DAPWIDTH` の既定と窓 | 7。窓の追跡順 | 51。全点を追跡 (§未対応・限界) |
| `Emittance->False` | DampingRate 以降の規則を返さない | 同じ (以前は無視) |
| `Union[…, SameTest->f]` | 先頭の要素だけ (K64 の誤り) | k64 は同じ、v2 は真値 |
| `Names["組み込みの名前"]` | `{}` (大域に作った記号だけ) | k64 は同じ、v2 は組み込みも含める |
| `WriteBeamLine[…, Name->"L2"]` | エラー (K64 の誤り) | k64 は同じ、v2 は受ける |

## SAD 自身の Packages を読むときの注意
- Python 3.14 (Tcl/Tk 9.0) でも [`test_tclbridge.py`](../pysad/tests/test_tclbridge.py)・[`test_plot.py`](../pysad/tests/test_plot.py)・[`test_gui_lang.py`](../pysad/tests/test_gui_lang.py)・[`test_corr.py`](../pysad/tests/test_corr.py) が通る (v2.1.0。Windows の Tk 9.0 は `$DefaultFixedFont` を `Courier New` にする: docs/22)。numba のキャッシュは `py -3.12 -m pysad.nbcache` (docs/12)。

## 未対応・限界
- 記号式の標準形は gs と同じではない (`0*x` は gs が `(0 x)` のまま、`x^0` も `x^0` のまま、pysad は 0 と 1 にする)。D の出力は 59 式で一致するが、別の式では形が違いうる。パターンの Orderless 照合は、数を先頭にした標準順→順列 (6 引数まで)、素の Blank と列パターンが両方あるときは標準順だけ (K64 の挙動に合わせて経験的に決めた)。
- `ResponseMatrix` の TRPT (加速、`gammab` の比) は実装済み (応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。gs と全桁: v2.0.0)。
- `OptimizeOptics` の最初に gs が表示する `{{`,System`},System`}` (BeginPackage の戻り値) は、pysad は `{$ContextPath,"`"}` (文脈を持たないため。差は残る)。
- Matrix.n (`LinearSolve`、`Inverse` の SADScript 版) は読まない。`SingularValues`・`MatrixBalance` は Python (`miscfn.py`) で、特異ベクトルの符号は各 u の先頭成分が負になるよう揃えた (gs と一致を 3 行列で確認)。
