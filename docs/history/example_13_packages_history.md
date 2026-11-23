# example/13_packages の履歴 (取り除いた記述)

md から、版・日付つきの修正・経過の記述を取り除いた。取り除いた元の行 (整理前) を、ファイルごとにそのまま残す。現在の説明は各 md を見る。

## example/13_packages/README.md

整理前の 10 行目:

> | [fn_MakeBump](fn_MakeBump.md) | Bump・MakeBump (垂直の閉バンプ) | ○ (2026-11-17 に修正。目標の 1 mm に一致、12 本) |

整理前の 11 行目:

> | [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (チューンの調整) | ○ (2026-11-17 に Output->0 で反復を出さないよう修正) |

整理前の 13 行目:

> | [fn_Rational](fn_Rational.md) | Rational・ContinuedFraction・SmallPrime | ○ (2026-11-17 に和・積を評価するよう修正) |

整理前の 15 行目:

> | [fn_Couple](fn_Couple.md) | Couple.n (ElementValues) | ○ (2026-11-17 に ElementValues を対応) |

## example/13_packages/fn_CorrectOrbit_Monitor_Steer.md

整理前の 19 行目:

> `R:rmx0`・`R:rmy0` は補正前の水平・垂直の rms、`rmy1` は垂直補正後、`rmx1` は水平補正後、`rmy2` は水平補正後の垂直の rms。`nsteer=24`、`max kick` は最大のキック [rad]。gs は `Element values are not expanded` という注意を出すことがある。pysad も 2026-11-17 の再実行 (v2.0.0) で同じ注意と `NOCALEXP` の表示を出し、gs と同じになった (`NOKEEPEXP`・EXPAND の gs に合わせた対応による。最初の版の `.pysad.txt` は出していなかった。補正の数値は変わらない)。

## example/13_packages/fn_Couple.md

整理前の 3 行目:

> 種別: パッケージ (Packages/Couple.n)・フォルダ 13_packages・ラティス FODO_ring (24 セル、3 GeV、周長 105.6 m)。pysad も動く (2026-11-17 に対応)

整理前の 25 行目:

> - 以前の pysad は `"K1"["QD"]` (文字列が頭部の式) を文字列の添字として評価して `invalid literal for int() with base 10: 'QD'` で落ちた。いまは評価せず式のままにし、CALC の前に Couple.n の `SetCoupledElements` を SADScript のまま実行する (`With` が `Hold`・純関数の本体に値を埋め込む、`Hold[…][[1]]` の評価し直しなどを gs に合わせた。docs/34)。

## example/13_packages/fn_MakeBump.md

整理前の 3 行目:

> 種別: パッケージ (Packages/bump.n、correction.n から読まれる)・フォルダ 13_packages・ラティス errors_ring_pkg (誤差入りの 24 セル、3 GeV、周長 105.6 m)。pysad も動く (2026-11-17 に修正)

整理前の 28 行目:

> - 以前の pysad は `???Set: list assignment length mismatch` で失敗した。`b=Bump[…]` は `bump.n` の `(a_=Bump[…])^:=…` (Set の UpValue) で、pysad は (1) `Bump` を自動読み込みの表に持たず、(2) 読み込んだ UpValue を組み込みの `Set` の定義として捨てていた。直した (docs/34、`work_example_fixes_20261117.md`)。

## example/13_packages/fn_OptimizeOptics.md

整理前の 28 行目:

> - 出力: gs・pysad とも `Output->0` で反復を出さない (2026-11-17 に修正。`OptimizeOptics` は選択肢を `Join[{opt},{VariableRange->limit}]` の規則のリストで `DownhillSimplex` に渡すので、pysad は入れ子のリストを読む。`MaxIteration` も効く)。ただし gs は最初に ``{{`,System`},System`}`` (BeginPackage の戻り値) を 1 行出し、pysad は ``{$ContextPath,"`"}`` を出す (文脈を持たないため。結果には影響しない)。

## example/13_packages/fn_Rational.md

整理前の 20 行目:

> - **○ 演算も一致 (2026-11-17 に修正)**: `Rational[1,3]+Rational[1,6]` は gs・pysad とも `Rational[3,6]` (約分は `ReducedRational` が $p^2>b$ で打ち切るため、和は約分されない。gs も同じ)、`Rational[1,3]*Rational[3,4]` は `Rational[1,4]`。以前の pysad は評価せず式のまま残した。`rational.n` の `Plus[Rational[…],…]^:=…` などの UpValue が、組み込みの Plus・Times・Minus・Divide・Subtract の定義として捨てられていたため。`-a`・`1/a`・`a-b`・`a/b`・`0.5+a` も gs と一致する (13 式、試験 `test_rational_arithmetic_uses_package_upvalues_like_gs`)。


## 追加 (整理の 2 回目)

### example/13_packages/fn_LDSolve.md

整理前の 32 行目:

> `AccuracyGoal` を変えたときの収束、$A$ が $s$ に依存しない場合の速い経路、複素数の $A$ (docs/34 で pysad の修正点に挙がっている) は試していない。


## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### fn_CorrectLocal.sad・fn_CorrectOrbit_Micado.sad・fn_PeriodicR.sad・fn_CSROYMatrix.sad (旧 v2_next/D2_*.sad)
> ! v2next D2: 12 セルの FODO …、v2next D2 項目 4: …、v2next D2 項目 5 (B60): …
### fn_CorrectOrbit_Bump.sad・fn_CorrectOrbit_Micado_variants.sad (旧 v2_next/P5_micado_*.sad)
> ! v2next P5 項目 5: … D2_micado.sad のラティスを使う。gs (K64) と pysad k64 の比較 (.gs.txt と .pysad_k64.txt)。
.pysad_k64.txt は削除した (k64 の出力は `PYSAD_COMPAT=k64` で再現できる)。D2_micado.sad は fn_CorrectOrbit_Micado.sad の名前にした。
### fn_OYStability_pipeline.sad・fn_OYStability_pipeline_unstable.sad (旧 v2_next/P5_oy_pipeline*.sad)
> ! v2next P5 項目 4: OYStability の後半 …
### fn_packages_compare.sad (旧 v2_next/B2_packages_compare.sad)
> ! v2next B2 項目 8: Packages/functions.n・Misc.n・element.n・beamline.n の主な関数を gs と比べた …
