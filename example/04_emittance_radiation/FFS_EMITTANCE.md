# FFS EMIT (EMITTANCE)
種別: FFS コマンド / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz / フラグ: RFSW ON, RAD ON。

## 何を見る例か
`EMIT;` は、閉軌道、拡張 Twiss、$U_0$、減衰、平衡エミッタンス、エネルギー拡がり、バンチ長、シンクロトロン振動数を一度に表示する。この値を放射積分 ($I_1\ldots I_5$) の式で検算する。

## ラティスについて
- `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで RF 電圧を 4 MV にした `FODO_ring_V4MV.sad` を追加した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## スクリプト
```
USE RINGQ;  CELL;  CALCULATE;
RFSW;  RAD;
EMIT;
```
全文は `FFS_EMITTANCE.sad`。式の検算は `formula_check.py` が行う (SAD を使わず、3x3 行列で $I_1\ldots I_5$ を自前で積分し、結果を `formula_check.txt` に書く)。

## 結果の読み方 (gs の出力から)
| 行 | 値 | 意味 |
|---|---|---|
| `Energy loss per turn U0` | 1.5632657 MV | 一周のエネルギー損失 |
| `Effective voltage Vc` | 4.0000000 MV | RF 実効電圧 |
| `Equilibrium position dz` | 38.315657 mm | 同期位置 |
| `Momentum compact. alpha` | 0.0221837 | 運動量圧縮 (6D、$-R_{56}/C$) |
| `Synchrotron frequency` | 78426.782 Hz (f0=2838943.6 Hz) | $\nu_s=0.0276253$ |
| `Damping time (sec)` | X 1.401206E-3, Y 1.351875E-3, Z 6.646461E-4 | $\tau_i$ |
| `Damping partition number` | 0.9651, 1.0003, 2.0346 | $J_x,J_y,J_\epsilon$ |
| `Emittance X` | 1.73577E-7 m | $\epsilon_x$ (Y は 1e-38 の雑音) |
| `Energy spread` | 0.00119140 | $\sigma_\delta$ |
| `Bunch Length` | 16.0597178 mm | $\sigma_z$ |

## 物理・解析式
放射積分 ($\rho$ は曲率半径、$\eta$ は分散、$\mathcal H=\gamma\eta^2+2\alpha\eta\eta'+\beta\eta'^2$):

$$I_1=\oint\frac{\eta}{\rho}ds,\quad I_2=\oint\frac{ds}{\rho^2},\quad I_3=\oint\frac{ds}{|\rho|^3},\quad I_4=\oint\frac{\eta}{\rho}\Big(\frac1{\rho^2}+2k\Big)ds+(\text{端の項}),\quad I_5=\oint\frac{\mathcal H}{|\rho|^3}ds$$

$$U_0=\frac{C_\gamma E^4 I_2}{2\pi},\quad \tau_i=\frac{2E_0T_0}{J_iU_0},\quad J_x=1-\frac{I_4}{I_2},\ J_y=1,\ J_\epsilon=2+\frac{I_4}{I_2}$$

$$\epsilon_x=C_q\gamma^2\frac{I_5}{J_xI_2},\quad \sigma_\delta^2=C_q\gamma^2\frac{I_3}{J_\epsilon I_2},\quad \alpha_c=\frac{I_1}{C},\quad \nu_s=\sqrt{\frac{h\,\alpha_c\,eV\cos\phi_s}{2\pi\beta^2E_0}},\ \sin\phi_s=\frac{U_0}{eV},\quad \sigma_z=\frac{\alpha_c c}{2\pi f_s}\sigma_\delta$$

$C_\gamma=\frac{4\pi}{3}\frac{r_e}{(mc^2)^3}=8.846\times10^{-5}$ m/GeV³、$C_q=\frac{55}{32\sqrt3}\frac{\hbar}{mc}=3.832\times10^{-13}$ m、$T_0=C/c$、$h=f_{RF}/f_0$。

`formula_check.py` の結果 (自前で積分した $I$ から式で計算した値):

| 量 | 式 (I4 の端の項 係数 -1) | 式 (係数 -2) | RADINT gs | EMIT gs | EMIT pysad |
|---|---|---|---|---|---|
| $U_0$ [eV] | 1.56327e6 | 1.56327e6 | 1.563e6 | 1.56327e6 | 1.56327e6 |
| $J_x$ | 0.941466 | 0.9641 | 0.941466 | 0.9651 | 0.9651 |
| $\tau_x$ [s] | 1.43601e-3 | 1.4023e-3 | 1.436e-3 | 1.40121e-3 | 1.40121e-3 |
| $\tau_y$ [s] | 1.35195e-3 | 1.35195e-3 | 1.352e-3 | 1.35187e-3 | 1.35187e-3 |
| $\tau_z$ [s] | 6.56755e-4 | 6.64057e-4 | 6.568e-4 | 6.64646e-4 | 6.64646e-4 |
| $\epsilon_x$ [m] | 1.7767e-7 | 1.73499e-7 | 1.7767e-7 | 1.73577e-7 | 1.73577e-7 |
| $\sigma_\delta$ | 1.18311e-3 | 1.18967e-3 | 1.1831e-3 | 1.1914e-3 | 1.1914e-3 |
| $\alpha_c$ | 0.0221368 | 0.0221368 | 0.022137 | 0.0221837 | 0.02218 (表示桁) |
| $\nu_s$ | 0.0275961 | 0.0275961 | 0.027596 | 0.0276253 | 0.0276253 |
| $\sigma_z$ [m] | 0.0159506 | 0.016039 | 0.0159506 | 0.0160597 | 0.0160597 |

- 「係数 -1」の列は、SAD の RADINT と同じ $I_4$ (= 0.0802372。自前の積分も 0.0802373 と一致)。式は RADINT の表示と相対 7e-5 以内で一致した ($\epsilon_x$ は 1.6e-6、$J_x$ は 3e-8)。
- **EMIT (6D の一周写像から求める) の $J_x=0.9651$ は RADINT の 0.9415 と 2.4% 違う。** 原因を切り分けた。RFSW OFF でも同じ (0.96487) なので、RF 空洞の分散位置が原因ではない。$I_4$ の端の項 (BEND の端の角 $\epsilon$ の項 $\eta\tan\epsilon/\rho^2$) の係数を -1 (RADINT) から -2 にすると $J_x=0.9641$ になり、EMIT の 0.9651 に 0.1% まで近づく。つまり RADINT (INTGRL.f) は端の項を片側 1 回分しか入れておらず、EMIT のほうが一周写像に忠実、というのが現時点の見立て (**ソースでは確認していない**。端の点の曲率を 1/(2ρ) とみなす扱いかもしれない)。この差は $\epsilon_x$ に 2.4%、$\sigma_\delta$ に 0.7% の差として出る。
- $\alpha_c$: RADINT の $I_1/C=0.022137$ と EMIT の 6D の $\alpha=0.0221837$ は 0.2% 違う (SAD の MomentumCompaction は 6D の値で、RADINT とは別。docs/10 の level 1)。$\nu_s$ に EMIT の $\alpha$ を入れた式の値は 0.0276253 で EMIT と 6 桁一致。
- $\sigma_z=\alpha_cc\,\sigma_\delta/(2\pi f_s)$ に EMIT の値を入れると 16.079 mm (EMIT の Bunch Length 16.060 mm、0.12% 差)。

## gs と pysad の比較
○ 一致。EMIT 報告の全行が表示桁で一致する。差は 0 同然の量だけ: Emittance Y (gs 2.8e-38、pysad 2.4e-50。どちらも 0 とみなせる)、Beam tilt、`alpha` の表示桁 (gs `.0221837`、pysad `.02218`)。

## 未確認・限界
- 使った条件: 3 GeV、周長 105.6 m、$V_{RF}=4$ MV、$f_{RF}=500$ MHz、RAD ON、FLUC は既定、結合なし (平面リング)、六極は K2=0。
- $J_x$ の RADINT 対 EMIT の 2.4% の差の原因は、端の項の取り扱いと推定した段階 (ソース未確認、他の BEND 形状での確認なし)。どちらが物理的に正しいかは、追跡で測っていない。
- 縦の量 ($\sigma_z$) は線形の式で、非線形性は入れていない。
