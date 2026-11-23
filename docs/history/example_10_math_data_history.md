# example/10_math_data の履歴 (取り除いた記述)

md から、版・日付つきの修正・経過の記述を取り除いた。取り除いた元の行 (整理前) を、ファイルごとにそのまま残す。現在の説明は各 md を見る。

## example/10_math_data/README.md

整理前の 14 行目:

> | [fn_Gamma_Erf](fn_Gamma_Erf.md) | Gamma・Erf・GammaRegularized ほか | ○ (LogGamma1 は 2026-11-17 に対応。`InverseErf` は v2 が真値を返し gs の 1e-8 の誤差と違う: B47。k64 モードは一致) |

## example/10_math_data/fn_Bessel.md

整理前の 38 行目:

> $n$ が大きいとき、$x$ が負のときは試していない。複素引数: 最初の版では pysad 未対応と書いたが、v2 で Bessel J/Y/I/K (実数の次数) が複素引数に対応した (`docs/47_v2_changes.md`。`BesselJ[0,1+I]` が値を返すことを 2026-11-17 に確認)。この例の表には入れておらず、gs との比較はしていない。gs の J0(100) が誤る境界の引数 (どの $x$ から誤るか) は調べていない。

## example/10_math_data/fn_Fit_FitGaussian_FitEmit.md

整理前の 54 行目:

> - `FitGaussian`: 全桁一致 (k64 モード)。v2 (既定) は `InverseErf` が真値になる (B47) ので、σ が 8 桁 (.8032255801 と gs の .8032255747)、他の値が 1e-8 の桁で変わる。`FitEmit` も alpha・beta・emittance が 1e-8 以内の差 (gs の alpha -1.0081380609、v2 -1.0081380566)。`.pysad.txt` は v2 の出力。`PYSAD_COMPAT=k64` では保存してあった旧い出力と同じ。

## example/10_math_data/fn_Gamma_Erf.md

整理前の 46 行目:

> - `InverseErf[.5]`: gs は真値 (0.476936276204) と相対 1e-8 ずれる (K64 の近似の誤差。B47)。pysad は **v2 (既定) では真値 .47693627620447 を返す** (K64 のバグを直した。`docs/47_v2_changes.md`)。`PYSAD_COMPAT=k64` では gs と同じ .476936271460083。上の表の pysad 列と相対差は k64 の値 (v1 と同じ)。`.pysad.txt` は v2 の出力なので、この行だけ違う。

## example/10_math_data/fn_NIntegrate.md

整理前の 31 行目:

> 複素数: $\int_0^\pi e^{ix}dx=2i$。gs は `(8.56519716264E-16+1.999999999999996 I)`。pysad は `(8.77770078844E-16+1.999999999999996 I)` で、虚部は gs と同じ (2026-11-17 に複素の被積分関数に対応)。実部は 1e-15 程度の丸めの雑音。

## example/10_math_data/fn_PolynomialFit.md

整理前の 35 行目:

> 重み付き、次数が点数に近い場合、ほぼ共線のデータは試していない。x の平均が 0 のデータでは、K64 が `PolynomialFit` の係数を誤る (B18) ことを v2 で直した (v2 は numpy.polyfit と 1e-13、k64 モードは K64 の値)。この例の x は 0〜5 (平均 2.5) なので影響しない。

## example/10_math_data/fn_SpecialFunctions_D.md

整理前の 46 行目:

> 1 点ずつの確認のみ。引数の範囲 (負、大きい)、次数の高い場合は試していない。複素引数: 最初の版では pysad 未対応と書いたが、v2 で Airy・Bessel・Hypergeometric0F1・1F1・2F1 (実数のパラメータ) が複素引数に対応した (`docs/47_v2_changes.md`。`AiryAi[1+I]` が値を返すことを 2026-11-17 に確認)。複素の PolyLog・HypergeometricU・Zeta は v2 でも未対応 (`Zeta[2+I]` は TypeError で確認)。この例の表には複素引数を入れておらず、gs との比較はしていない。`HermiteH`・`ExpIntegralEi` の gs の誤りの条件 (どの引数で誤るか) は調べていない。

## example/10_math_data/var_PhysicalConstants.md

整理前の 45 行目:

> SADHelp の定義値は CODATA 2018 相当。`ProtonRadius` の定義式は推定。旧版の `ElectronMass` (510999.06) からの変更 (docs/34) は、この例の範囲では gs の値が 510998.95 であることまで確認。


## 追加 (整理の 2 回目)

### example/10_math_data/README.md

整理前の 11 行目:

> | [fn_NIntegrate](fn_NIntegrate.md) | NIntegrate (滑らか・特異・複素) | △ (特異点の収束の警告の数値が違う。複素数は 2026-11-17 に対応して一致) |

### example/10_math_data/fn_Gamma_Erf.md

整理前の 44 行目:

> - `LogGamma1[10.5]` = 16.29200047656724 は gs と一致 (○。2026-11-17 に対応。`LogGamma1[x]=LogGamma[x+1]`、負の引数の複素数の虚部も gs と同じ)。

### example/10_math_data/fn_NIntegrate.md

整理前の 37 行目:

> - 複素数 (`Exp[I*x]` を 0〜π): gs (8.6e-16+1.999999999999996 I)、pysad も虚部 1.999999999999996 で一致 (○。2026-11-17 に対応。実部は 1e-15 程度の丸めの雑音)。


## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### fn_SpecialFunctions_complex.sad (旧 v2_next/B2_special_functions_complex.sad)
> (* v2next B2 項目 2: 複素引数の PolyLog・Zeta・HurwitzZeta・HypergeometricU。…
### fn_SpecialFunctions_complex_scan.sad (旧 v2_next/P5_specfn_cx_scan.sad)
> ! v2next P5 項目 3: 複素引数の特殊関数の走査。… gs (.gs.txt) と pysad k64・v2 の値を比べる。
