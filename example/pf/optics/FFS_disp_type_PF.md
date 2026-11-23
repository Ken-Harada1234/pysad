# FFS disp / type (PF): 表示コマンド
種別: FFS コマンド (KEK の SAD チュートリアルの話題をもとにした)。加速器: PF 2.5 GeV リング (`lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。

## 何を見る例か
`disp 場所` は、その場所のツイス・分散と、最後に必ずリング終点 (`$$$`) の行を出す。`type 要素名` は要素の定義 (長さ・強さ) を出す。

## スクリプト
```
USE RING;  cell calc emit;
disp QF;       ! QF の最初の出現 (QF.1)
type QF;  type B02;  type SFH;
disp INJ;      ! 入射点
```

## 結果の読み方
- `disp QF`: 列は AX BX NX EX EPX 要素名 長さ 強さ s(m) AY BY NY EY EPY DetR 番号。QF.1 (149 番目、s = 38.448 m): βx 5.1495 m、βy 7.9238 m、ηx 0.4731 m、νx 1.88056 (2π で割った位相)。
- 最後の行が `$$$` (リング終点): βx 11.9837 m、βy 4.9956 m、νx 9.59753、νy 5.27896。
- `type QF` → `QUAD QF = (L =.4 K1 =.748319526303)`、`type B02` → `BEND B02 = (L =1.942696 ANGLE =.2243994752564138 E1 =.5 E2 =.5)` (E1=E2=0.5 は矩形磁石)、`type SFH` → `SEXT SFH = (L =.1 K2 =3.9606428...)`。
- `disp INJ`: 入射点 (280 番目、s = 66.642 m) で βx 9.4486 m、βy 1.1552 m、ηx 0.6343 m。

## 物理・解析式
式の照合はなし (表示のみ)。

## gs と pysad の比較
**v2 モードと k64 モード**: `.pysad.txt` は既定の v2 モードの出力。v2 は K64 のバグを含まない計算なので (`docs/28_v2_vs_k64.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

○ 数値は全桁一致。`disp` の最後の 3 列 (EY・EPY・DetR) だけ、gs は 1e-16 の雑音を出し pysad は 0 を出す。

## 未確認・限界
`disp` の列の意味は SADHelp とチュートリアルの説明による。
