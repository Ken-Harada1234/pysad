# 33_specfn_fit_misc から取り除いた履歴の記述

`docs/33_specfn_fit_misc.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # 特殊関数・数学関数、Fit・FindRoot・NIntegrate、DownhillSimplex など

取り除いた部分: 「最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。」

元の行:

> 最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。経緯は `history/` にある)

## ## 1. 特殊関数 (`pysad/specfn.py`)

取り除いた部分: 「HurwitzZeta は v2.1.0 で実装した** (下の」

元の行:

> - **複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で実装した** (下の「複素の特殊関数」)。**未対応**: 複素の次数の PolyLog、Re a ≤ 0 の複素の a・b の HypergeometricU (SADError)、z = 0 の U。PolyLog の整数次 n>=2 の z>1 は実装済み (gs と同じく切断線の下側の複素数。試験 test_polylog_above_one_is_complex_like_gs)、HurwitzZeta の q が非正の整数 (`HurwitzZeta[2,-1]` は gs が 2.6449 = 極の項を除いた和、pysad は NaN。2026-11-17 に gs で確認)、D.n の記号微分 (`D`、`Derivative`、`ReduceExpression`)、

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「**(b) モードによらず直したもの (K64」

元の行:

> **(b) モードによらず直したもの (K64 が誤値を返す、または pysad の不具合。k64 モードでも同じ)**

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「相対 1e-13 (pysad は以前から真値) |」

元の行:

> | `Hypergeometric2F1[.5,.5,1,.99]`、`BesselJ[.5,30.]` (B43) | 2.348414、-.143924709693 (級数の打ち切り) | 2.352715816780、-.143929653370 | 2F1(½,½;1;m) = 2K(m)/π を AGM で、J_{1/2}(x) = sqrt(2/πx) sin x。相対 1e-13 (pysad は以前から真値) |

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「(以前は符号が逆) |」

元の行:

> | `LegendreP` の奇数 m の符号 | (gs は NaN で参照なし。pysad の誤り) | `LegendreP[3,1,.4]` = +.2749545416973507 (以前は符号が逆) | m = 1 ± 1e-6 の式の値と閉じた式 (scipy の `lpmv` に `(-1)^m` を掛けていた誤りを直した) |

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「(scipy の `lpmv` に `(-1)^m` を掛けていた誤りを直した) |」

元の行:

> | `LegendreP` の奇数 m の符号 | (gs は NaN で参照なし。pysad の誤り) | `LegendreP[3,1,.4]` = +.2749545416973507 | m = 1 ± 1e-6 の式の値と閉じた式 (scipy の `lpmv` に `(-1)^m` を掛けていた誤りを直した) |

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「TypeError (pysad の従来) |」

元の行:

> | 複素引数 | TypeError (pysad の従来) | Bessel J/Y/I/K (実数の次数)、Hypergeometric0F1・1F1・2F1 (実数のパラメータ)、Airy が複素引数で使える | 積分表示・級数・Wronskian と 1e-12〜1e-13 (`test_complex_argument_bessel_hypergeometric`) |

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「従来の μ = ln z の展開は |μ| → 2π で収束が遅く z < -12・z > 60 で精度が落ちていた (z = -140、s = -1.5 で相対誤差 1、|z| ≥ 150 は NaN)。`specfn._polylog_far` を足した (s < 2.2 は反転公式、s ≥ 2.2 は積分表示)。」

元の行:

> **PolyLog の非整数次・大きい |z|**: 従来の μ = ln z の展開は |μ| → 2π で収束が遅く z < -12・z > 60 で精度が落ちていた (z = -140、s = -1.5 で相対誤差 1、|z| ≥ 150 は NaN)。`specfn._polylog_far` を足した (s < 2.2 は反転公式、s ≥ 2.2 は積分表示)。使える範囲: 実数の z、実数の非整数 s。使えない範囲: 複素数の z (SADError)、負の非整数 s が大きい (s < -10 程度) ときの精度は未確認。(上の §5 の「|z|<150」の記述は、この修正前の範囲。)

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「(上の §5 の「|z|<150」の記述は、この修正前の範囲。)」

元の行:

> **PolyLog の非整数次・大きい |z|**: μ = ln z の展開は |μ| → 2π で収束が遅く z < -12・z > 60 で精度が落ちる (z = -140、s = -1.5 で相対誤差 1、|z| ≥ 150 は NaN) ので、そこは `specfn._polylog_far` を使う (s < 2.2 は反転公式、s ≥ 2.2 は積分表示)。使える範囲: 実数の z、実数の非整数 s。使えない範囲: 複素数の z (SADError)、負の非整数 s が大きい (s < -10 程度) ときの精度は未確認。(上の §5 の「|z|<150」の記述は、この修正前の範囲。)

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「**複素の特殊関数** (v2.1.0。[`pysad/specfn_cx.py`]」

元の行:

> **複素の特殊関数** (v2.1.0。[`pysad/specfn_cx.py`](../pysad/specfn_cx.py)。両モード。mpmath はこの環境に無いので、級数・Euler–Maclaurin・解析接続・積分表示を自前で実装した。試験 [`test_v2_B2_specfn_cx.py`](../pysad/tests/test_v2_B2_specfn_cx.py)、例 [`example/v2_next/B2_special_functions_complex.sad`](../example/v2_next/B2_special_functions_complex.sad)):

## ## 6. v2 モードで変わった値と、モードによらない修正 (`test_v2_B_specfn.py` 13 件、`test_v2_B_script.py`)

取り除いた部分: 「実軸上の z > 1 は従来どおり切断線の下側 |」

元の行:

> | `PolyLog[s,z]` | \|z\| ≤ 0.5 は冪級数、\|ln z\| ≤ 4 は μ = ln z の展開 (110 項まで)、それより大きいときは倍角の公式の再帰 (\|z\| = 1e8 で 5 段)、s が 0 以下の整数は有理関数。**負の実数の z は上側 (arg z = +π)**。実軸上の z > 1 は従来どおり切断線の下側 | 実数の次数、\|z\| ≤ 1e8 (s = 0.5〜7.3) | 冪級数・`scipy.special.spence`・積分表示・整数次の反転公式と 2e-11 |
