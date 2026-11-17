# 10_math_data 一覧
数値・データ処理の例。当てはめ・補間・積分・最小化、特殊関数、行列関数、定数。ラティスは使わない (`FFS;` の後に直接書く)。式で値を検算できるものは、md に式と数値の照合表がある (`R:` で始まる出力行が表の元)。

○一致 / △差あり (理由は各 .md) / ×gs 失敗 / ×pysad 未対応・失敗

| 例 (NAME.sad・.gs.txt・.pysad.txt・.md) | 内容 | gs と pysad |
|---|---|---|
| [fn_Fit_FitGaussian_FitEmit](fn_Fit_FitGaussian_FitEmit.md) | Fit・FitGaussian・FitEmit | ○ (DesignMatrix の表示だけ pysad が数値微分の値。v2 は `InverseErf` が真値で FitGaussian が 1e-8 の桁で違う: B47) |
| [fn_PolynomialFit](fn_PolynomialFit.md) | PolynomialFit (2 次・1 次) | ○ |
| [fn_Spline](fn_Spline.md) | Spline (自然・傾き指定)、導関数、積分 | ○ |
| [fn_NIntegrate](fn_NIntegrate.md) | NIntegrate (滑らか・特異・複素) | △ (特異点の収束の警告の数値が違う。複素数は 2026-11-17 に対応して一致) |
| [fn_DownhillSimplex](fn_DownhillSimplex.md) | DownhillSimplex (SADHelp の例ほか) | ○ (反復の出力形式が違う。A・B は最小に届かない K64 の動作) |
| [fn_Elementary](fn_Elementary.md) | 初等関数 (三角・指数・対数・双曲線・逆関数) | ○ |
| [fn_Gamma_Erf](fn_Gamma_Erf.md) | Gamma・Erf・GammaRegularized ほか | ○ (LogGamma1 は 2026-11-17 に対応。`InverseErf` は v2 が真値を返し gs の 1e-8 の誤差と違う: B47。k64 モードは一致) |
| [fn_Bessel](fn_Bessel.md) | BesselJ・I・K | △ (gs が J0(100)・I0(0) で誤り) |
| [fn_BesselY_JZero](fn_BesselY_JZero.md) | BesselY・BesselJZero | ×gs 失敗 (NaN または応答なし。pysad は値を返す) |
| [fn_SpecialFunctions_D](fn_SpecialFunctions_D.md) | LegendreP・HermiteH・Ei・Si・Fresnel・Airy・Zeta・PolyLog ほか | △ (gs が HermiteH・Ei で誤り) |
| [fn_NumericalFunctions](fn_NumericalFunctions.md) | Abs・Round・Mod・FractionalPart ほか | ○ |
| [fn_Matrix_LinearAlgebra](fn_Matrix_LinearAlgebra.md) | Det・Inverse・LinearSolve・Transpose・Outer・Inner・Dot | ○ (丸め誤差の差) |
| [fn_Eigensystem](fn_Eigensystem.md) | Eigensystem (対称・回転・1 ターン行列) | ○ |
| [fn_SingularValues](fn_SingularValues.md) | SingularValues (既定は逆数を返す) | ○ |
| [var_PhysicalConstants](var_PhysicalConstants.md) | 物理定数 13 個と定義式の再計算 | ○ |
| [var_Constants](var_Constants.md) | Pi・Degree・GoldenRatio・I・INF・NaN | △ (`I*I` の表示、`NaN==NaN`) |

数値の桁・表示の細かい違いは追っていない (大きな差だけ md に書いた)。
