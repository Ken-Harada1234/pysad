# FFS_TYPE_ATTRIBUTE (TYPE と ATTRIBUTE)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。

## 何を見る例か
- `TYPE 要素名…` (略 `T`): 要素の定義を MAIN 部の入力と同じ形で表示する (0 の項目は既定キーワード以外を省く)。
- `ATTRIBUTE 要素名パターン…` (略 `ATTR`): 要素の既定キーワードの現在値、最小・最大 (`VARY` や変数の限界値)、`COUPLE` された親要素と係数を表で出す。

## スクリプト
```
USE RING; CELL; CALC;
TYPE QF QD B RF SF;
ATTRIBUTE QF Q*;
ATTRIBUTE B;
```

## 結果の読み方
- `TYPE`: 各要素 `型 名前 =(キーワード 値 …)`。BEND の `E1`/`E2`、CAVI の `VOLT` `FREQ`、SEXT の `K2 =0` など。
- `ATTRIBUTE QF Q*`: 表の列は `Element Keyword Value Minimum Maximum Couple Coefficient`。QF.1 は K1 = 0.7、限界 ±1e10、Couple の欄の `<--` は「結合していない」の意味。(パターンが 2 つあるので QF.1 は 2 回出る: gs。)
- `ATTRIBUTE B`: 既定キーワードは BEND では `ANGLE` (0.261799)。

## 物理・解析式
なし。

## gs と pysad の比較
値は一致。表示の違い 2 点だけ (桁・形式): (1) `CAVI` の `FREQ` が gs では `5e+08`、pysad は `500000000`。(2) `ATTRIBUTE QF Q*` で QF.1 が重複する行を pysad は 1 回だけ出す (gs は 2 回)。

## 未確認・限界
- `TYPE` の `PROF` などリスト値の表示 (docs/03 の限界) は未確認。`ATTRIBUTE` の Couple 欄は FFS_COUPLE (02) で別に確認した (pysad は結合を表示しない)。
