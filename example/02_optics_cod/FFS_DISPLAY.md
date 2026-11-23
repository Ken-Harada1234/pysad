# FFS_DISPLAY (DISPLAY: 光学関数の表示)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`DISP [モード] [パターン] [始点 終点]` (略 `DISP`) は、直前の `CALC` の結果を表にする (自身は計算しない)。
モード語: なし (標準: AX BX NX EX EPX、長さ、値、s、AY BY NY EY EPY、DetR)、`R` (R1〜R4 の x-y 結合行列)、`O` (軌道 DX DPX DY DPY)、`E` (β の極値だけ)、`G` (幾何 Gx Gy Gz とχ)、`GA` (γ関数)、`Z` (縦方向)、`B` (ビームサイズ)、`D` (1 行ダンプ)、`REFERENCE`/`DREFERENCE`。
位置は要素名・パターン・番号・範囲 (`1 8`) で指定する。

## スクリプト
```
DISP 1 8;            DISP QF.1|QD.1;          DISP B* 1 30;      DISP R 1 6;
DISP E 1 40;         DISP G 1 8;              DISP GA 1 6;       DISP Z 1 4;
```
(`B*` は B で始まる要素のパターン。`QF.1|QD.1` は 2 つの成分。範囲 `1 8` は要素番号 1〜8。)

## 結果の読み方
- 標準表: 各行は**その要素の入口**の値 (`ORIGIN` の行がリングの始点 = 周期解の入口 AX=-2.1663、BX=5.02518)。`s(m)` は入口の s、`Value` は既定キーワードの値 (QF なら K1)。
- `R` モード: x-y 結合の R1〜R4 (ここは結合がないので 1e-16 程度の雑音)。`E` モード: β の極値をとる行だけ。
- `G` モード: `Gx Gy Gz`、`s`、`Chi1 Chi2 Chi3` (ORIGIN から DS.1 の入口 0.2 m、QF.1 の入口 0.3 m …)。
- `GA`: `GMX` (= γx = (1+αx²)/βx。ORIGIN で 1.13287) と `GMX*L` (ほぼ「自然クロマティシティの寄与」)。`Z`: CALC6D 用で、ここ (CALC4D) では何も出ない。

## 物理・解析式
`GA` の `GMX` = $(1+\alpha_x^2)/\beta_x$: ORIGIN で $(1+2.1663^2)/5.02518=1.13287$ ✓ (表の値)。

## gs と pysad の比較
**○ 一致 (2026-11-17 確認)。** 標準・`R`・`O`・`E`・`G` に加え `GA`・`Z`・`P`・`A`・`B` (`DISP BEAM`、EMIT 後)・`D` (28 光学関数の 1 行ダンプ)・`REFERENCE`・`DREFERENCE` の表を gs と同じ見出し・桁で出す (表の数値は全桁一致。雑音 1e-16 の列だけ pysad は 0 を出す)。`DISP B REFERENCE` は gs と同じ `???-FFS-Info-REF and DREF not implemented for DISP B.`。

## 未確認・限界
- `DISP B` (ビームサイズ) は EMIT のあと (FFS_BEAMSIZE)、`REFERENCE` モードは FFS_REFERENCE で扱う。66 行ごとの停止の対話 (`ALL` で抑止) は非対話の実行なので確認していない。
