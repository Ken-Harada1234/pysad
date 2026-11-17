# fn_Emittance (Emittance[])
種別: 関数 / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz / フラグ: RFSW, RAD。

## 何を見る例か
`Emittance[]` が返す規則のリストから、`/.` で個々の値を取り出す。EMIT コマンドと同じ計算 (temit) を、プログラムから使える形にしたもの。

## ラティスについて
- `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで RF 電圧を 4 MV にした `FODO_ring_V4MV.sad` を追加した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## スクリプト
```
RFSW; RAD;
e = Emittance[];
Print[Emittances /. e, MomentumSpread /. e, BunchLength /. e, EnergyLossU0 /. e,
      MomentumCompaction /. e, DampingRate /. e, Tunes /. e, RfVoltageVc /. e, Stable /. e];
```

## 結果の読み方 (gs)
| 規則 | 値 | 単位・意味 |
|---|---|---|
| `Stable` | 1 | 安定 |
| `Tunes` | {-0.303759, -0.093944, -0.027625} | 小数部 ($\nu_x=6.6963$ の小数部 0.6963-1)。第 3 が $-\nu_s$ |
| `EnergyLossU0` | 1563265.7 | eV |
| `RfVoltageVc` | 4000000.0 | eV |
| `EquilibriumPosition` | 0.038316 | m |
| `MomentumCompaction` | 0.022184 | 6D の値 |
| `BucketHeight` | 0.010035 | $\Delta p/p$ |
| `DampingRate` | {-2.514E-4, -2.606E-4, -5.300E-4} | 一周あたり |
| `Emittances` | {1.7358E-7, 2.825E-38, 1.9069E-5} | m ($\epsilon_z\approx\sigma_z\sigma_\delta$) |
| `MomentumSpread` | 0.001191 | $\sigma_\delta$ |
| `BunchLength` | 0.016060 | m |
| `TuneShiftByRadiation` | {-5.112E-5, 5.4112E-5, -1.795E-7} | |

## 物理・解析式
- 一周あたりの減衰率 $=-\frac{U_0}{2E_0}J_i$。$U_0/E_0=5.2109\times10^{-4}$ から $J_x=0.9651$、$J_y=1.0003$、$J_\epsilon=2.0346$ (三つの和は 4.0000)。
- $\sigma_z\sigma_\delta=0.016060\times0.001191=1.9128\times10^{-5}$ m で、表の $\epsilon_z$ (1.9069E-5) と 0.3% 違う。表示桁の丸め (0.03% 程度) では説明できない。原因は未確認 (`Emittances` の第 3 成分が $\sigma_z\sigma_\delta$ と同じ定義かどうか確かめていない)。
- 他の式 ($U_0,\tau,\epsilon_x,\sigma_\delta,\nu_s$) と放射積分との照合は FFS_EMITTANCE.md の表。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

○ 一致。`Emittances` の第 2 成分は 0 同然の雑音 (gs 2.8E-38、pysad 1.2E-67)。`OrbitAtExit` も 1e-17 の雑音で、値の差はない。

## 未確認・限界
- 平面リング・結合なし・3 GeV・4 MV。`POL`、`INTRA` は使っていない (`PolarizationTime` などは 0)。
- `Tunes` は小数部で、符号付き ($\nu_x=6.6963$ に対し -0.3038)。
