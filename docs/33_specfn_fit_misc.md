# 特殊関数・数学関数、Fit・FindRoot・NIntegrate、DownhillSimplex など

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

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
| 小物 | Restrict MinMax NaNQ BitwiseAnd/Or/Xor DigitQ LetterQ RealListQ SymbolName Complex Duplicate FlattenAt Inner Fourier InverseFourier InversePower Union1 RotateRight1 `Range$` |

- **K64 と値を合わせたもの** (k64 モードの動作。v2 モードの違いは下の「v2 モードで変わった値」): `InverseErf` (B47。k64 モードは相対 1e-8 の誤差まで再現)、`ArcCosh[.5]` の符号 (-1.047 I)、`ArcTanh[2.]` = NaN、`Inner` (gs の `Inner[f,x,y,g]` は `g@@f[x,y]`)、`Fourier` の符号 (gs は e^{+2πi jk/n}/√n)。
- **K64 が誤値を返すので真値にしたもの**: docs/21 B43〜B46。
- **未対応**: 複素引数の PolyLog・HypergeometricU・Zeta (scipy に複素の実装が無いか、分岐の向きを確かめる独立な参照が無い。SADError・TypeError)。PolyLog の整数次 n>=2 の z>1 は実装済み (gs と同じく切断線の下側の複素数。試験 test_polylog_above_one_is_complex_like_gs)、HurwitzZeta の q が非正の整数 (`HurwitzZeta[2,-1]` は gs が 2.6449 = 極の項を除いた和、pysad は NaN。2026-11-17 に gs で確認)、D.n の記号微分 (`D`、`Derivative`、`ReduceExpression`)、
  `Short`、`Skip`、`Inequality`、`MessageName`、`ToInputString`。

## 2. Fit・FindRoot・NIntegrate (`pysad/fitroot.py`)
- K64 の tffit1 / tfnewton を移植。K64 は式を記号微分して導関数にするが、ここでは常に数値微分の経路 (K64 が D を作れないときと同じ) を使う。
  gs で `Fit` の導関数が解析的な場合との差は線形な式で相対 1e-9、FindRoot の収束値で 1e-14。FitGaussian の出力は gs と 5e-10 で一致。
- `Fit` のベクトル独立変数 (`{{x1,x2},y,dy}`) は対応する。gs 側は narg を出して最小にならない値を返す (docs/21 B52) ので gs とは照合せず、numpy.linalg.lstsq の解・残差と一致を確かめた (係数 1e-6、ChiSquare 1e-8。`test_b52_fit_vector_independent_variable`)。`FindRoot` の複素数は gs が壊れた値を返す (B61)。実数は gs と一致、**複素は k64 モードでは TypeError のまま、v2 モードは複素の Newton 法で解く** (下の「複素の FindRoot」)。`FitGaussian` の Plot オプションは実装済み (下記 §5 の `FitGaussian`)。
- NIntegrate は D.n の分割細分 + 三次スプライン積分の移植 (Spline->False も)。scipy の CubicSpline を使うので、収束した値は gs と 1e-14 で合う。被積分関数が複素数 (`Exp[I*x]`) でも、実部と虚部を別々に積分するのと同じ (スプラインは線形) 形で動く (gs と虚部 1e-15 で一致。実部は 1e-16 の雑音)。

- gs と全桁一致を確認 (2026-11-16 夜): `Spline[{{x,y},…}]` (自然スプライン)、`Sort[list, Greater]`、`Scan`。`Spline[data,0,0]` は gs も pysad も未評価のまま。

## 3. DownhillSimplex (`pysad/miscfn.py`、Packages/optimize.n)
Amoeba の移植。gs と 1e-13 (Rosenbrock、500 反復) で一致。長く縮んだ単体 (VariableRange で拘束し続けるもの) では丸め誤差から 1e-6 ずれることがある。
Output>0 の途中経過の書式は `$FORM='10.6'` ではなく `%10.6f`。`OptimizeOptics` は Python では翻訳せず、`Packages/optimize.n` を SADScript のまま実行して gs と一致する (docs/34)。

## 4. 小さな組み込み
StringToStream (文字列を開く)、TemporaryName、ListRandom (`{"SAD"}`)、MemoryCheck (`{0,0,0,0}`)、`$Input`=77、`$Output`=-1、FixedPointList、SelectCases。

## 5. ファイル・文字列・その他の小さな関数 (`pysad/fileops.py`、`miscfn.py`)
- **ファイル**: FileQ DirectoryQ FileNames DirectoryName ToFileName CreateDirectory DeleteDirectory RenameDirectory CopyDirectory CopyFile RenameFile DeleteFile FileDate SetFileDate FileType FileByteCount `$HomeDirectory` (K64 は外部コマンド。Python の os・shutil。FileByteCount は K64 が行数を返す誤りなのでバイト数。B50)。失敗は `$Failed`。
- **init.n・functions.n**: ARRAY Len Put Environment HomeDirectory StringInsert Diagonal Expand SymbolNameRoot HeldPart MkSecureTemp (ARRAY、Len、StringInsert、Diagonal、Expand は gs と一致)。
- **システム**: Directory SetDirectory Replace FindIndex Attributes ToInputString Context SetEnv GetEnv GetGID GarbageCollect ToDate FromDate、`SingularValues` `MatrixBalance` (gs と一致)。
- `FitGaussian[d, Plot->True]` は `{…, Show[Plot (当てはめの曲線), ListPlot (ビン), opt]}` を返す (Fit.n と同じ)。gs と同じ `List` 5・キャンバス 74 項目 (rectangle 1、line 40、text 13、oval 20)。`Plot[InverseErf[x] …]` のため `Vector` が `InverseErf` などの特殊関数を扱えるようにした (vecfn は specfn の後に入れる)。
- `PolyLog[s, z]` の非整数 s で z>1 (切断線の下側の複素数) と z<-1 (実数)、`HurwitzZeta[s, q]` の q<=0 の非整数 (`(n+q)^-s` は主値の複素数) を実装 (gs と相対 1e-14。`Li_s(z) = Γ(1-s)(-μ)^{s-1} + Σ ζ(s-k) μ^k/k!`、|μ|<2π の範囲 = |z|<150)。試験 `test_polylog_noninteger_outside_unit_disk_and_hurwitz_negative_q_like_gs`。

## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

根拠は gs ではなく独立な計算 (閉じた式、級数、AGM、積分表示の数値積分、Wronskian などの恒等式、numpy の最小二乗)。mpmath はこの環境に無いので使っていない (scipy.special と、scipy に依らない独立な計算を突き合わせた)。gs は K64 の動作の確認にだけ使った。条件: 実数の倍精度。複素引数は実数のパラメータだけ確かめた。例: `example/v2_B/B1_special_functions.sad`。

**(a) v2 モードだけ値が変わるもの (k64 モードは K64 の値)**

| 関数 | k64 モード | v2 モード | 根拠 |
|---|---|---|---|
| `InverseErf[.3]` (B47) | .272462711811 (Newton の係数が sqrt(π/2)、早い打ち切りで相対 1e-8 の誤差。`FitGaussian` も影響) | .2724627147267543 (scipy の erfinv)。`InverseErfc` は小さい z で桁落ちしない erfcinv | erf(InverseErf(x)) = x が 4e-16 以内、正しい係数 (sqrt(π)/2 e^{y²}) の Newton 法と 1e-13 (`test_b47_inverseerf_v2_true_value`、`test_fitgaussian_uses_true_inverseerf_in_v2`、k64: `test_b47_inverseerf_k64_unchanged`) |
| `PolynomialFit` (B18) | x の平均が 0 のとき係数が壊れる | numpy.polyfit と 1e-13 (docs/04) | `test_b18_polynomialfit_mean_zero_v2` |
| `FindRoot` の複素初期値 (B61) | TypeError | 複素の Newton 法 (docs/04) | `test_b61_findroot_complex_v2` |
| ``Names["Stat`*"]`` | gs と同じ ``{"Stat`"}`` | その文脈の記号の一覧 | docs/04 |

**(b) モードによらず直したもの (K64 が誤値を返す、または pysad の不具合。k64 モードでも同じ)**

| 関数 | K64 (gs) の動作 | pysad の値 | 根拠 |
|---|---|---|---|
| `Hypergeometric2F1[.5,.5,1,.99]`、`BesselJ[.5,30.]` (B43) | 2.348414、-.143924709693 (級数の打ち切り) | 2.352715816780、-.143929653370 | 2F1(½,½;1;m) = 2K(m)/π を AGM で、J_{1/2}(x) = sqrt(2/πx) sin x。相対 1e-13 (pysad は以前から真値) |
| `ExpIntegralEi` (B44・B54) | `Ei[1.5]` = .164-2.985 I、`Ei[-.1]` = -2.2265 | 3.301285449130、-1.822923958419 (実数) | 級数、負の側は -E1(x)、大きい x は漸近級数。1e-10〜2e-12 |
| `HermiteH`・`LegendreP`・`Binomial` (B45) | `HermiteH[3,.7]`・`LegendreP[3,1,.4]`・`LegendreP[2,2,.3]`・`Binomial[-3,2]` が NaN | -5.656、+.274954541697、2.73、6 | H3 = 8x³-12x、P_3^1 = -(1-z²)^½ (15z²-3)/2、P_2^2 = 3(1-z²)、二項係数の積 |
| `HypergeometricU`・`BesselY`・`BesselK`・`BesselJ` の負の次数・`BesselJZero` (B46) | `U[1,2,.5]` = INF、`BesselY[0,.1]` = NaN、`BesselK[1.5,.7]` = INF、他は応答なし | 2、-1.534238651350、1.806573612779、-.557936507910 (`BesselJ[-1,1.5]`)、.145918137967 (`BesselY[1,2.5]`)、2.404825557696 | U(1,2,z) = 1/z、Y_0 の級数、K_{3/2} の閉じた式、Wronskian J_1 Y_0 - J_0 Y_1 = 2/(πx) |
| 大きい引数の `BesselJ/K/Y`・`FresnelS`・`AiryAi` (B53) | `BesselJ[0,33.]` 相対 3e-3、`FresnelS[6.1]` = 1.3e7、`AiryAi[15.5]` = -16 | 真値 (scipy) | 積分表示 (台形則)、Simpson 積分、Ai Bi' - Ai' Bi = 1/π。x = 33・100・1000 の Bessel、x = 6.1・20 の Fresnel |
| `PolyLog` (B54) | 複素数 (虚部 1e-15 の雑音つき) | 実数 (z ≤ 1)。z > 1 は切断線の下側 (Li_2(2) = π²/4 - iπ ln 2)、z < -1 は実数。**非整数次は |z| の上限なし** | Li_2(½) = π²/12 - ln²2/2、Fermi–Dirac 積分、Im = -π ln^{s-1} z/Γ(s)、Li_{s-1} = z dLi_s/dz。s ≤ 7.9、|z| ≤ 1e12 で相対 1e-11 |
| `LegendreP` の奇数 m の符号 | (gs は NaN で参照なし。pysad の誤り) | `LegendreP[3,1,.4]` = +.2749545416973507 (以前は符号が逆) | m = 1 ± 1e-6 の式の値と閉じた式 (scipy の `lpmv` に `(-1)^m` を掛けていた誤りを直した) |
| 複素引数 | TypeError (pysad の従来) | Bessel J/Y/I/K (実数の次数)、Hypergeometric0F1・1F1・2F1 (実数のパラメータ)、Airy が複素引数で使える | 積分表示・級数・Wronskian と 1e-12〜1e-13 (`test_complex_argument_bessel_hypergeometric`) |
| `FixedPointList` (B48)・`Fit` のベクトル独立変数 (B52) | 未評価・narg | Misc.n の定義どおり・最小二乗解 | Newton の反復の再計算、numpy.linalg.lstsq |

**PolyLog の非整数次・大きい |z|**: 従来の μ = ln z の展開は |μ| → 2π で収束が遅く z < -12・z > 60 で精度が落ちていた (z = -140、s = -1.5 で相対誤差 1、|z| ≥ 150 は NaN)。`specfn._polylog_far` を足した (s < 2.2 は反転公式、s ≥ 2.2 は積分表示)。使える範囲: 実数の z、実数の非整数 s。使えない範囲: 複素数の z (SADError)、負の非整数 s が大きい (s < -10 程度) ときの精度は未確認。(上の §5 の「|z|<150」の記述は、この修正前の範囲。)

**複素の FindRoot** (v2): 初期値に複素数があると、式を変数について正則と仮定して、前進差分のヤコビアンと最小二乗 (numpy lstsq) の Newton 法 (歩幅を半分にする減速つき) で解く。z² = -1 → i、e^z = 2 → ln 2、z³ = 1、連立 zw = i・z + w = 2 が残差 1e-12 以下。使い方・使えない範囲は docs/04。

**未対応のまま**: 複素引数の PolyLog・HypergeometricU・Zeta。
