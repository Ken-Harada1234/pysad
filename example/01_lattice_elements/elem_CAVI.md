# elem_CAVI (高周波空洞)
種別: elem。フォルダ: 01_lattice_elements。ラティスは `FODO_ring.sad` と同じ光学 (24 セル、周長 105.6 m、3 GeV) に `CAVI RF` を 1 台 (定義はファイル内: `HARM=176`)。フラグ: CELL、RFSW (FFS 進入時に ON)、RAD は OFF (Emittance[] が放射損失 U0 を内部で計算する)。

## 何を見る例か
`CAVI RF = (L=0 VOLT=1.5E6 HARM=176 PHI=0)`: `VOLT` は電圧 [V]、`HARM` は周波数を周回周波数の整数倍で指定 (代わりに `FREQ` [Hz] でも書ける)、`PHI` は位相。
RF は `CALC` の Twiss には効かない (リングではドリフト扱い)。**`Emittance[]` でシンクロトロンチューンとして現れる。**
電圧を 1.5 → 3 → 6 MV と上げ、$\nu_s$ が式どおり変わることを見る。

## スクリプト
```
CAVI RF = (L=0 VOLT=1.5E6 HARM=176 PHI=0);
USE RING; CELL; CALC;   TYPE RF;
em = Emittance[];   Print["Tunes (x,y,s) = ", Tunes/.em];
RF VOLT 3.0E6;   CALC;   em = Emittance[];   ...      RF VOLT 6.0E6;  ...    RF PHI 0.5;  ...
```
`RF VOLT 3.0E6;` は FFS から要素の値を変える。`Tunes/.em` で Emittance[] の戻り (規則のリスト) から値を取り出す。

## 結果の読み方
`Tunes (x,y,s)` は小数部のチューン (x, y は負の符号付きの小数部。`-0.3037` は 6.6963 の 1 に対する補数)、3 番目が $\nu_s$ (SAD は符号付きで返す)。
`Stable`=1 は安定。`U0` は 1 周あたりの放射損失 [eV]、`Vc` は RF 電圧。

## 物理・解析式
1 周の放射損失 $U_0$ を差し引いた同期位相 $\sin\phi_s=U_0/V$ で、
$$\nu_s=\sqrt{\frac{h\,\alpha_c\,eV\cos\phi_s}{2\pi E}},\qquad h=176,\ E=3\ \text{GeV},\ U_0=1563265.7\ \text{eV}$$
(α_c は `MomentumCompaction`。この値は RF 電圧で少し変わる: 3 MV で 0.022169、6 MV で 0.022211。同期粒子の軌道が平衡位置にずれるため。)

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| νs (V=3 MV) | 0.0230219 | 0.023022 | 0.023022 |
| νs (V=6 MV) | 0.0346604 | 0.03466 | 0.03466 |
| νs (V=1.5 MV < U0) | 0 | 0 | 0 |

- **1.5 MV は $U_0$=1.56 MV より小さい**ので同期位相が存在せず (`Vc` < `U0`)、$\nu_s$=0。
- `PHI` を 0.5 にしても結果が変わらない (最後の行): SAD は放射損失とつり合う平衡位相を自分で決める (SADHelp の `DPHI` の項の ts)。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

チューン・U0・α は全桁一致。**`TYPE RF` の表示だけ違う**: gs は `HARM =176`、pysad は `FREQ =499654096.67` (HARM を周波数に直した値で表示)。

## 未確認・限界
- RF のあるリングでの Twiss は RF の影響を受けない (SAD の仕様)。6 次元の Twiss は CALC6D (flag_CALC4D_CALC6D) を参照。
- 空洞のフリンジ (`FRINGE`、`DISFRIN`)、`V1`/`V20` (横方向の電圧勾配)、`DX`/`DY`、`RANVOLT` は未確認。
- 追跡での加速・減速は 05_tracking の例で扱う。
