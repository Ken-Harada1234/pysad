# flag CALC6D
種別: フラグ / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz / フラグ: RFSW, RAD, CALC6D。

## 何を見る例か
`CALC6D` を立てて `CALCULATE` すると 6D の光学 (縦の $\beta_z$ とチューン $\nu_z$、$z$ との結合) が得られる。4D の結果 (CALC4D、既定) と比べる。

## ラティスについて
- `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで RF 電圧を 4 MV にした `FODO_ring_V4MV.sad` を追加した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## スクリプト
```
RFSW; RAD;
CALC6D;  CELL;  CALCULATE;
Print[Twiss["NZ","$$$"]/2/Pi, Twiss["BZ","$$$"], Twiss["AZ","$$$"]];
Print[Twiss["EX","ORIGIN"], Twiss["ZX","ORIGIN"], Twiss["ZPX","ORIGIN"]];
```

## 結果の読み方 (gs)
| | 4D (CALC4D) | 6D (CALC6D) |
|---|---|---|
| BX (周回の末尾) | 6.426926 | 6.425423 |
| NX | 6.696329 | 6.696241 |
| AX | -2.506166 | -2.506118 |
| NY, BY | 4.906056, 1.919589 | 同じ |
| NZ/2π | (4D では出ない) | -0.027625 |
| BZ, AZ | | 13.530158, 0.087006 |
| EX (ORIGIN) | 0.625610 (EX は 4D でも出る) | 0.625610 |
| ZX, ZPX (ORIGIN) | | 0.002153, -1.474E-5 |

6D では、横の β と $\nu_x$ が RF と放射でわずかに (1e-4 の桁) 変わる。

## 物理・解析式
縦のチューン $\nu_z=-0.027625$ は、`Emittance[]` の `Tunes` の第 3 成分と一致する (fn_Emittance)。縦の $\beta$ 関数は

$$\beta_z=\frac{|\alpha_c|C}{2\pi\nu_s}=\frac{0.0221837\times105.6}{2\pi\times0.027625}=13.50\ \text{m}$$

(CALC6D の BZ=13.530 と 0.2% 違う。式は近似で、$\nu_s$ に 6D の値、$\alpha_c$ に 6D の $-R_{56}/C$ を使っているため)。バンチ長は $\sigma_z=\sqrt{\epsilon_z\beta_z}=\sqrt{1.9069\times10^{-5}\times13.530}=0.016063$ m で、EMIT の Bunch Length 0.016060 m に一致する。

## gs と pysad の比較
○ 一致 (表示 6 桁)。

## 未確認・限界
- この例は RADCOD OFF (閉軌道に放射損失を入れない) の 6D だけ。RADCOD ON の CALC6D は見ていない。
- `CALC6D;` 単独は再計算しない (SAD の癖。CALCULATE が要る。docs/29 A36)。4D の CALC は縦の列を作り直さない。
- $\beta_z$ の式は近似で、結果は 0.2% 違う。
