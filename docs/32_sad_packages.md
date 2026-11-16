# SAD の Packages/*.n を SADScript のまま使う (軌道補正・最適化・数学パッケージ)

最終更新: 2026-11-16 (現状のスナップショット。経緯は `history/` にある)

pysad の SADScript インタプリタが、SAD の `Packages/*.n` を最初に使われたときに読む (`pysad/pkgload.py`)。Packages ディレクトリは
環境変数 `PYSAD_SAD_PACKAGES`、なければ pysad の隣の `SAD-master/Packages`。見つからなければ何も登録しない (組み込みだけで動く)。

| パッケージ | 提供する関数 | 状況 |
|---|---|---|
| correction.n (+ bump.n、corutility.n) | CorrectOrbit、CorrectLocal、MakeBump、Bump、Monitor、Steer、CorStat、RmsReport、NearestMonitors、ReadMonitorData … | gs と一致 (下記) |
| optimize.n | OptimizeOptics、Amoeba、FindClosedOrbit (DownhillSimplex は Python で書いた `miscfn.py`) | OptimizeOptics が gs と 1e-10 で一致 |
| rational.n、tableform.n、Couple.n、emit.n の Beamstrahlung 系 | Rational、TableForm、BeamstrahlungUpsilon … | 動く (gs と一致) |
| D.n | Spline (K64 の Spline$ は miscfn.py)、記号微分 `D`、`Derivative`、`ReduceExpression`、`Integrate`、Spline、特殊関数 (HermiteH など) | `D` を含め動く (59 式が gs と文字列まで一致。下記) |
| LDSolve.n、Optics.n | LDSolve、Optics の TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo | gs と一致 (下記「LDSolve・Optics.n」)。PeriodicR は未確認 |

組み込み (Python) で持つ関数 (LinearSolve、Binomial、Spline など) は、パッケージの定義で上書きしない (`pkgload._protect_builtins`)。

## 組み込み (K64 の Fortran) の翻訳 (`pysad/corr.py`)
- **ResponseMatrix** (pgrmat.f): COD 応答行列 (リングは周期解、ライン、`X Y PX PY MX MY MPX MPY`、`Tune`)。
- **LinearSolveConditioned** (pgsolvcond.f、pmicad.f): 最小二乗 (tsvd)、条件付き (psvdcond)、Micado (ベストコレクタ法)。
- **Master[]** (pgmast.f)、**FFS["VAR"]** (tfvars.f): FREE した変数の表。

gs との比較 (`tests/test_corr.py`、`co5_corr.sad`: 12 セル、24 BPM、24 ステアラ、DX/DY の誤差つき QUAD):
応答行列 `ResponseMatrixG` の行、CorrectOrbit の補正量 (Y、X、Micado 6 本) は gs と全桁一致 (補正量 1e-7 相対)。OptimizeOptics の結果 (K1、チューン) も一致。
なお gs でも pysad でも、`CorrectOrbit` が `LINE["K0"]` に足した値は `CALC NOEXP;` では効かない (SAD は "Element values are not expanded" と表示)。

## LDSolve・Optics.n (2026-11-16 夜に gs と比較、`local/wake/op1.sad`・`ld2.sad`・`ev1.sad`)
- Optics.n の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo`: gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し、pysad は 0)。
- LDSolve.n: pysad は `DiagonalMatrix` が複素数の列を受けず失敗していた (`float() argument must be … not 'complex'`) ので直した。あわせて行列の演算 (Dot・Inverse) の結果は、gs が要素ごとに虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻すようにした (以前は全要素が実のときだけ)。
  直したあと、LDSolve の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (試験は `test_evalfixes.py` の 2 件)。
- 未確認: `PeriodicR` は呼び出し方が不明で比較できなかった (gs も引数エラー)。

## パーサ・評価器の修正 (パッケージを動かすために gs に合わせたもの)
- パターン: Plus・Times は Flat に平らにして Orderless に照合、`f_[y__]` より `p_Plus` を先に、同名の変数を繰り返す定義 (`f[x_,x_]`) を `f[_Symbol,_]` より先に。`Plus[x_]` は評価しない。`{0 ..}` の Repeated。`//.` は不動点まで、`Unevaluated$[rules]` を規則にする、`f/@(a+b)` は評価し直す、`f[Unevaluated$[g[y,x]]]` は呼び出し側の y x の値を埋め込む (評価しない)。`-x` は `Times[-1,x]`、`x/2` は `(.5 x)`、`(2 x)^2` は `(4 x^2)`。
- 前置の `[x]` は括弧 (`Map[f,[g[x]]]`)、`f[x]g[y]`・`(a)(b)`・`2#1`・`x^2 y` は積 (D.n)。
- `Part[l,{i,j},k]` は要素ごとに残りの添字を適用。`Head` (Null と記号は Symbol、複素数は Complex、Function、Rule)。`VectorQ[x,test]`・`MatrixQ[x,test]`。
- `With[{g:=e/2}, f[e_]:=g*3]` は定義の中へ式を埋め込む。
- `CALC NOEXP`・`CALC EXP` を受理 (pysad は元から要素値を展開し直さないので `CALC` と同じ)。
- 数の既定の書式 (`autos`) を gdtoa に合わせた: `1e6` → `1e6`、`1e10` → `1e+10`、`123456789012345678.` → `123456789012345680`。`$FORM=" "` は既定と同じ。
- 物理定数 (`ElectronCharge`、`ElectronMass` = 510998.95 eV (旧: 510999.06)、`ElectronRadius`、`PlanckConstant` …) をスクリプトから参照できる。

## 未対応・限界
- 記号式の標準形は gs と同じではない (`0*x` は gs が `(0 x)` のまま、`x^0` も `x^0` のまま、pysad は 0 と 1 にする)。D の出力は 59 式で一致するが、別の式では形が違いうる。パターンの Orderless 照合は、数を先頭にした標準順→順列 (6 引数まで)、素の Blank と列パターンが両方あるときは標準順だけ (K64 の挙動に合わせて経験的に決めた)。
- `ResponseMatrix` の TRPT (加速、`gammab` の比) は未対応。
- Matrix.n (`LinearSolve`、`Inverse` の SADScript 版) は読まない。`SingularValues`・`MatrixBalance` は Python (`miscfn.py`) で、特異ベクトルの符号は各 u の先頭成分が負になるよう揃えた (gs と一致を 3 行列で確認)。
