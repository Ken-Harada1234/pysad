# elem_MULT (多極磁石 MULT)
種別: elem。フォルダ: 01_lattice_elements。ラティスは `FODO_ring.sad` と同じ光学のリングを 2 つ (QUAD+BEND 版と MULT 版。定義はファイル内、RF なし)。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`MULT` は多極を 1 つにまとめた汎用要素 (`K0`, `K1`, `K2`…、スキュー `SK0…`、`ANGLE` で偏向も可)。
`K1` だけを与えれば四極と同じ、`ANGLE` `E1` `E2` を与えれば偏向磁石と同じ線形光学になる。QF・QD・B を MULT に置き換えたリングで、チューンが元のリングと同じことを確かめる。

## スクリプト
```
MULT  MF = (L=0.3 K1=0.7) MD = (L=0.3 K1=-0.6);
MULT  MB = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELLM = (SF DS MF D2 MB D2 MD DS SD D1);   LINE RINGM = (ORIGIN 24*CELLM);
USE RINGQ; CELL; CALC;   ...   USE RINGM; CELL; CALC;   ...   TYPE MF MB;
```

## 結果の読み方
2 行のチューン: 「QUAD/BEND ring」と「MULT ring」。`TYPE MF MB` は MULT の定義 (MB に `K1 =0` が補われる)。

## 物理・解析式
チューンは解析式で書ける (FFS_CALCULATE と同じ FODO リングなので、`tools/ana_fodo_ring.py` の 2x2 行列の積から)。

| 量 | 式 (2x2 行列) | gs (QUAD) | gs (MULT) | pysad (QUAD) | pysad (MULT) |
|---|---|---|---|---|---|
| NX | 6.696329 | 6.696329 | 6.696329 | 6.696329 | 6.696329 |
| NY | 4.906056 | 4.906056 | 4.906056 | 4.906056 | 4.906056 |

## gs と pysad の比較
全桁一致 (MULT 版も QUAD+BEND 版と 7 桁で同じ)。

## 未確認・限界
- K2 以上、スキュー成分、`ANGLE`≠0 で `K1`≠0 の複合機能、`ROTATE`/`DX`/`DY`、放射 (`DISRAD`)、追跡は確認していない。
- docs/28 によれば `MULT` の ANGLE≠0 の Twiss・Emittance・追跡は別の機会に gs と照合済み (この例では再確認していない)。
