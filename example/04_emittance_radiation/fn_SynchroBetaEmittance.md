# fn_SynchroBetaEmittance
種別: 関数 / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz / フラグ: RFSW, RAD。**平面リング (結合なし) のみ**。

## 何を見る例か
`SynchroBetaEmittance[{nus0, nus1, dnus}, AzimuthalModes->7]` が、シンクロ・ベータ結合 (RF が分散のある場所にあるときや、色収差のあるとき) を含む平衡エミッタンスを、$\nu_s$ を与えて計算する。返り値は `{{nus, emitx, emity, emitxp, emityp, conv}, ...}`。

## ラティスについて
- `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで RF 電圧を 4 MV にした `FODO_ring_V4MV.sad` を追加した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## スクリプト
```
RFSW; RAD;
r = SynchroBetaEmittance[{0.005, 0.025, 0.01}, AzimuthalModes->7];
```

## 結果の読み方 (gs)
| $\nu_s$ | emitx | emity | emitxp | conv |
|---|---|---|---|---|
| 0.005 | 1.3380E-6 | -1.57E-19 | 1.3377E-6 | 3.0706E-9 |
| 0.015 | 1.3285E-6 | -1.10E-19 | 1.3283E-6 | 6.289E-10 |
| 0.025 | 1.3393E-6 | -6.76E-21 | 1.3394E-6 | 5.100E-10 |

conv は正 (収束)。`Emittance[]` の $\epsilon_x=1.7358\times10^{-7}$ に対して **約 7.7 倍**。

## 物理・解析式
このリングは RF 空洞 (ORIGIN の直後) が分散のある場所 ($\eta_x=0.6256$ m、$\beta_x=6.43$ m、$\alpha_x=-2.51$ で $\mathcal H=0.443$ m) にある。空洞でエネルギーが変わると betatron 振幅が変わるため、通常の放射積分の $\epsilon_x=C_q\gamma^2I_5/(J_xI_2)$ (1.7e-7) には入らない励起が加わる。目安として $\mathcal H_{cav}\sigma_\delta^2=0.443\times(1.19\times10^{-3})^2=6.3\times10^{-7}$ m で、差 $1.34\times10^{-6}-1.7\times10^{-7}=1.2\times10^{-6}$ と同じ桁。ただし **この桁の一致だけで、係数の厳密な式の検算はしていない**。$\nu_s$ にほとんど依らない ($0.005\to0.025$ で ±0.4%) のも、この説明と矛盾しない。

## gs と pysad の比較
○ 一致。emitx は pysad も 1.3380E-6, 1.3285E-6, 1.3393E-6 (5 桁)、conv は 3.0718E-9 (gs 3.0706E-9) など。emity は両方 0 同然の雑音 (符号が違う: gs -1.57E-19、pysad +1.68E-19)。

## 未確認・限界
- **結合のあるラティスでは gs が SIGSEGV で落ちる** (docs/03 §175 と docs/10 に記録がある。QT・ずれ入りの 8+4 セルのリングでの確認)。pysad も `ValueError: cannot convert float NaN to integer` で落ちる。この例は平面リングだけ。
- 不安定な $\nu_s$、`emiout`、RFSW の組み合わせは見ていない。
- 7.7 倍の物理的な内訳は上の桁の推定だけで、確かめていない。
