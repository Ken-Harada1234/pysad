# example/kek_ls の履歴 (整理で md から取り除いた記述)

元の文書の文脈つきで、取り除いた文章をそのまま入れてある。

## kek_ls/design/KEKLS_match_cell_tune.md

```
## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
```

## kek_ls/design/KEKLS_match_ring_geo.md

```
## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
```

## kek_ls/design/KEKLS_scale_ncell.md

```
## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
```

## kek_ls/design/KEKLS_match_injcell.md

```
## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
```

## kek_ls/design/KEKLS_scale_mini.md

```
## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で 1e-8 (相対) 以内しか変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
```

## kek_ls/design/KEKLS_match_repeat.md

```
判定: **○** (比べた行 12、最大の相対差 2.7e-05、k64 モードの pysad。v2.4.0 の修正のあとの再生成)。3 回とも同じ値で gs と同じ (A: QFA1H 0.5115137183 / gs 0.5115137397、リングのチューン 48.58012266 / 48.58012266)。残りの差は FitFunction の片側差分の数値雑音。
```

```
- 原因は未確認 (pysad の FIT 状態の持ち越し、または FitFunction の重みの扱い)。`KEKLS_tune_survey_cell`・`KEKLS_scan_target` の gs と pysad の違いはこれが原因とみられる。
- 元の台本 (matching2_01 の `FitFunction` が NaN のとき pysad が動かない) と同じ現象の疑いがある。
```

```
## 更新 (v2.4.0 の修正のあと。`.pysad.txt` は再生成)
原因を確定して直した (docs/46 §4)。(1) `FFS["fit MKSFA1 BX 3."]` の末尾の数 `3.` をパーサが落としていた、(2) gs の tfsolv は FitFunction の行を不等式型 (最初は解から外す) として扱う (k64)、(3) gs は USE のたびに DP を MARK の値に戻す (ChromCor の `DP = 0;` の持ち越しが pysad で 2 回目以降を変えていた)。結果: gs と pysad の A〜D の 12 行の最大の相対差は 2.7e-05 (以前 6.0e-02)、3 回とも同じ値 (gs と同じ)。残りの差は FitFunction の片側差分の数値雑音。以前の判定 (△) と上の「原因は未確認」は、この更新で解消した。
```

## kek_ls/design/KEKLS_tune_survey_cell.md

```
k64 モードの pysad。v2.4.0 の修正のあとの再生成)
```

```
## 更新 (v2.4.0 の修正のあと。`.pysad.txt` は再生成)
`KEKLS_match_repeat` の修正 (FitFunction の行の扱い・末尾の数のパース・USE での DP の戻し) を反映して再計算した。νx = 48.3・48.45 の 10 点は εx が gs と 3e-6 以内・K2(SDH・SFH) が 5e-5 以内で一致し、νx ≥ 48.6 の点は εx が 1e-6〜10%・K2 が 1e-5〜5% 違う (6 極の谷が浅く、反復の履歴の数値雑音で別の谷に落ちる。(48.30, 17.30) の εx は gs 411、pysad 411)。上の判定 (点ごとの一致なし) は、前半の 10 点については解消した。
```

```
違う (6 極の谷が浅く、反復の履歴の数値雑音で別の谷に落ちる)。
```

## kek_ls/display/KEKLS_scan_target.md

```
k64 モードの pysad。v2.4.0 の修正のあとの再生成)。設計値の目標の 3 行は gs と同じく行どうしで同じ値になった (以前の pysad は行ごとに違った)。ξ の値は gs −0.0216・pysad −0.0254 で違う (6 極の補正の谷が浅い)。
```

```
## 更新 (v2.4.0 の修正のあと。`.pysad.txt` は再生成)
`KEKLS_match_repeat` の修正を反映して再計算した。設計値の目標の 3 行 (AY 1.19・BY 6.4・EX 0.14) は、gs と同じく行どうしで同じ値になった (以前の pysad は行ごとに違った)。ξ の値そのものは gs と pysad で違う (gs ξx = -0.0216、pysad -0.0254)。6 極 (SDH・SFH) の補正の谷が浅く、反復の履歴の数値雑音で別の解に落ちるため。上の判定 (△) は変わらないが、行の不一致の原因の一つ (FitFunction の扱いと DP の持ち越し) は解消した。
```

## kek_ls/lattice/README.md

```
pysad の v2 (既定) モードは、K64 のバグを直した `Emittance[]` で
```

