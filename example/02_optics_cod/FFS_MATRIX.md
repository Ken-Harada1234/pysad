# FFS_MATRIX (MATRIX: 転送行列)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL、CALC4D (4x5 行列)。

## 何を見る例か
`MATRIX [SYMPLECTIC|PHYSICAL] [from to]` は from から to までの転送行列を表示する (CALC4D のときは 4x5: x px y py と運動量偏差 δ の 5 列目。5 列目はδに対する変位)。
`SYMPLECTIC` (既定) は規格化された座標での行列、`PHYSICAL` は TRPT (輸送路で加速があるとき) の物理座標の行列。リングでは同じ。
関数形 `TransferMatrix[from, to]` の戻り値 (4x4) から、1 セルのトレースと行列式、位相進みも確かめる。

## スクリプト
```
MATRIX SF.1 SF.2;       MATRIX PHYSICAL SF.1 SF.2;       MATRIX QF.1 QD.1;
m = TransferMatrix["SF.1","SF.2"];
Print["trace x = ", m[[1,1]]+m[[2,2]], "  trace y = ", m[[3,3]]+m[[4,4]]];
MATRIX ORIGIN $$$;          ! リング 1 周
```

## 結果の読み方
`Symplectic transfer matrix from SF.1 to SF.2` の 4 行 × 5 列。x-y が分離しているので 2x2 ブロックが対角に並ぶ。`trace x =`、`det x =` の行。
`MATRIX ORIGIN $$$` が 1 周の行列 (安定ならトレースが 2cos(2πν))。

## 物理・解析式
1 セルの周期行列: $\det M=1$、$\cos\mu=\tfrac12{\rm tr}M$。
$$\cos\mu_x=-0.18129,\quad \mu_x/2\pi=0.279014,\quad \nu_x/24=0.279014$$

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| tr Mx (x) | -0.362581 | -0.362581 | -0.362581 |
| tr My (y) | 0.56499 | 0.56499 | 0.56499 |
| det Mx | 1 | 1 | 1 |
| det My | 1 | 1 | 1 |
| 位相進み x (1 セル)/2π | 0.279014 | 0.279014 | 0.279014 |
| 位相進み y (1 セル)/2π | 0.204419 | 0.204419 | 0.204419 |
| arccos(tr/2)/2π (x) | 0.279014 |  |  |


## gs と pysad の比較
全桁一致 (`MATRIX` の表、`PHYSICAL`、`Print`)。1e-16 程度の雑音の要素 (x-y 結合要素の 0) は gs が `-1.479E-31` などを出し、pysad は 0 を出す (表示の差だけ)。

## 未確認・限界
- CALC6D の 6x6 行列は flag_CALC4D_CALC6D。TRPT (加速のある輸送路) の PHYSICAL と SYMPLECTIC の違いは未確認 (リングでは同じ)。
- `MATRIX` の引数が名前だけのとき (`NAME.1` 扱い、docs/33) は MATRIX QF.1 QD.1 の形で使った (一致)。
