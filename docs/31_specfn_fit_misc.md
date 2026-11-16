# 特殊関数・数学関数、Fit・FindRoot・NIntegrate、DownhillSimplex など

最終更新: 2026-11-16 (現状のスナップショット。経緯は `history/` にある)

## 1. 特殊関数 (`pysad/specfn.py`)
scipy.special で実数引数の値を返す。K64 の実装 (gamma.f、hg.f、Packages/D.n) の値と、gs で 458 式を突き合わせて一致 (`tests/specfn_ref.json`、`tests/test_specfn.py`)。

| 区分 | 関数 |
|---|---|
| Arc 双曲線・対数ガンマ | ArcSinh ArcCosh ArcTanh LogGamma LogGamma1 (= LogGamma[x+1]) Gamma (複素も) Gamma0 (= E1) XSin |
| ガンマ族 | GammaRegularized(P/Q) PolyGamma Beta Pochhammer Factorial2 Binomial Multinomial HarmonicNumber |
| ゼータ・多重対数 | Zeta HurwitzZeta PolyLog LerchPhi HurwitzLerchPhi BernoulliB |
| 超幾何 | Hypergeometric0F1 1F1 2F1 U PFQ と Regularized |
| Bessel | BesselJ Y I K、BesselJZero |
| D.n の特殊関数 | ExpIntegralE ExpIntegralEi LogIntegral SinIntegral CosIntegral SinhIntegral CoshIntegral WhittakerM WhittakerW HermiteH FresnelS FresnelC AiryAi/Bi/AiPrime/BiPrime LegendreP InverseErf InverseErfc ProductLog |
| 基本関数の修正 (sadscript) | Log (負の引数は複素)、ArcSin/ArcCos/ArcTan (範囲外・複素)、Power (負の底の非整数べき)、Erf/Erfc (複素)、Factorial (非整数)、Exp (オーバーフローは INF)、`EulerGamma` |
| 小物 | Restrict MinMax NaNQ BitwiseAnd/Or/Xor DigitQ LetterQ RealListQ SymbolName Complex Duplicate FlattenAt Inner Fourier InverseFourier InversePower Union1 RotateRight1 Range$ |

- **K64 と値を合わせたもの**: `InverseErf` (B47。相対 1e-8 の誤差まで再現)、`ArcCosh[.5]` の符号 (-1.047 I)、`ArcTanh[2.]` = NaN、`Inner` (gs の `Inner[f,x,y,g]` は `g@@f[x,y]`)、`Fourier` の符号 (gs は e^{+2πi jk/n}/√n)。
- **K64 が誤値を返すので真値にしたもの**: docs/19 B43〜B46。
- **未対応**: 複素引数の Bessel・超幾何・PolyLog、PolyLog の |z|>=150 (非整数次)。PolyLog の整数次 n>=2 の z>1 は実装済み (gs と同じく切断線の下側の複素数。試験 test_polylog_above_one_is_complex_like_gs)、HurwitzZeta の q<=0、D.n の記号微分 (`D`、`Derivative`、`ReduceExpression`)、
  `Short`、`Skip`、`Inequality`、`MessageName`、`ToInputString`。

## 2. Fit・FindRoot・NIntegrate (`pysad/fitroot.py`)
- K64 の tffit1 / tfnewton を移植。K64 は式を記号微分して導関数にするが、ここでは常に数値微分の経路 (K64 が D を作れないときと同じ) を使う。
  gs で `Fit` の導関数が解析的な場合との差は線形な式で相対 1e-9、FindRoot の収束値で 1e-14。FitGaussian の出力は gs と 5e-10 で一致。
- `Fit` のベクトル独立変数 (`{{x1,x2},y,dy}`) は対応するが、gs 側が誤った値を返す (docs/19 B52) ので照合なし。`FindRoot` の複素数は gs が壊れた値を返す (B61) ので実数だけ一致させ、複素は TypeError のまま。`FitGaussian` の Plot オプションは実装済み (下記 §4 の小さな組み込み)。
- NIntegrate は D.n の分割細分 + 三次スプライン積分の移植 (Spline->False も)。scipy の CubicSpline を使うので、収束した値は gs と 1e-14 で合う。被積分関数が複素数 (`Exp[I*x]`) でも、実部と虚部を別々に積分するのと同じ (スプラインは線形) 形で動く (gs と虚部 1e-15 で一致。実部は 1e-16 の雑音)。

- gs と全桁一致を確認 (2026-11-16 夜、`local/wake/v1.sad`): `Spline[{{x,y},…}]` (自然スプライン)、`Sort[list, Greater]`、`Scan`。`Spline[data,0,0]` は gs も pysad も未評価のまま。

## 3. DownhillSimplex (`pysad/miscfn.py`、Packages/optimize.n)
Amoeba の移植。gs と 1e-13 (Rosenbrock、500 反復) で一致。長く縮んだ単体 (VariableRange で拘束し続けるもの) では丸め誤差から 1e-6 ずれることがある。
Output>0 の途中経過の書式は `$FORM='10.6'` ではなく `%10.6f`。`OptimizeOptics` は Python では翻訳せず、`Packages/optimize.n` を SADScript のまま実行して gs と一致する (docs/32)。

## 4. 小さな組み込み
StringToStream (文字列を開く)、TemporaryName、ListRandom (`{"SAD"}`)、MemoryCheck (`{0,0,0,0}`)、`$Input`=77、`$Output`=-1、FixedPointList、SelectCases。

## 5. ファイル・文字列・その他の小さな関数 (`pysad/fileops.py`、`miscfn.py`)
- **ファイル**: FileQ DirectoryQ FileNames DirectoryName ToFileName CreateDirectory DeleteDirectory RenameDirectory CopyDirectory CopyFile RenameFile DeleteFile FileDate SetFileDate FileType FileByteCount `$HomeDirectory` (K64 は外部コマンド。Python の os・shutil。FileByteCount は K64 が行数を返す誤りなのでバイト数。B50)。失敗は `$Failed`。
- **init.n・functions.n**: ARRAY Len Put Environment HomeDirectory StringInsert Diagonal Expand SymbolNameRoot HeldPart MkSecureTemp (ARRAY、Len、StringInsert、Diagonal、Expand は gs と一致)。
- **システム**: Directory SetDirectory Replace FindIndex Attributes ToInputString Context SetEnv GetEnv GetGID GarbageCollect ToDate FromDate、`SingularValues` `MatrixBalance` (gs と一致)。
- `FitGaussian[d, Plot->True]` は `{…, Show[Plot (当てはめの曲線), ListPlot (ビン), opt]}` を返す (Fit.n と同じ)。gs と同じ `List` 5・キャンバス 74 項目 (rectangle 1、line 40、text 13、oval 20)。`Plot[InverseErf[x] …]` のため `Vector` が `InverseErf` などの特殊関数を扱えるようにした (vecfn は specfn の後に入れる)。
- `PolyLog[s, z]` の非整数 s で z>1 (切断線の下側の複素数) と z<-1 (実数)、`HurwitzZeta[s, q]` の q<=0 の非整数 (`(n+q)^-s` は主値の複素数) を実装 (gs と相対 1e-14。`Li_s(z) = Γ(1-s)(-μ)^{s-1} + Σ ζ(s-k) μ^k/k!`、|μ|<2π の範囲 = |z|<150)。試験 `test_polylog_noninteger_outside_unit_disk_and_hurwitz_negative_q_like_gs`。
