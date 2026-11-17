# calc_respx (PF): 軌道応答行列と解析式
種別: スクリプト (サンプル集 calc_respx.sad と calc_resp_analy.sad を一つにした)。加速器: PF 2.5 GeV リング (`lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。フラグ: `ON RAD FLUC COD RFSW RADCOD` (RFSW は例の中で切り替える)。

## 何を見る例か
ステアリングに蹴り角 K0 を与え、`cell calc` で閉軌道を求め、ビーム位置モニター (`PM*`、65 台) の位置 `Twiss["DX","PM*"]` を読む (calc_respx)。同じ量を解析式
$$\frac{x_i}{\theta_j} = \frac{\sqrt{\beta_i\beta_j}}{2\sin\pi\nu}\cos\bigl(\pi\nu-|\phi_i-\phi_j|\bigr)$$
で計算し (calc_resp_analy)、比べる。蹴り角は偏向磁石 `B*` の K0 (28 台)、ここでは B16・B25・B07 の 3 台に 1e-5 rad (10 μrad) を与えた (チュートリアルの `bkick` は電流から換算した値)。

## スクリプト
```
dtnx = LINE["NAME","B*"];                      ! 28 台の偏向磁石の名前
LINE["K0", dtnx[j]] = kick;  FFS["cell calc;"];  sx = Twiss["DX","PM*"];
LINE["K0", dtnx[j]] = 0;                        ! 戻す
an = -kick*0.5*Sqrt[bpm[i]*bb]*Cos[tunex*Pi - Abs[ppm[i]-pb2]]/Sin[tunex*Pi];   ! 解析式 (SAD の符号)
```
`bb`・`pb2` は蹴った磁石の中央 (`LINE["POSITION",名前]+0.5`) の βx と位相 (チュートリアルは入口の値)。`bpm`・`ppm` はモニターの βx と位相 (`Twiss["NX",...]` は 2π で割る前の位相)。
`FFS["OFF RFSW"]` で RF を切った場合と `ON RFSW` の場合を回す。

## 結果の読み方
- **SAD では K0 > 0 の蹴りで DX は負**になるので、解析式に -1 を掛けた。
| 蹴った磁石 | PM#1 の DX (SAD) | 解析式 | max abs DX | rms (SAD - 解析式) |
|---|---|---|---|---|
| B16 | -9.5809 μm | -9.6028 μm | 26.25 μm | 0.0211 μm |
| B25 | +10.6160 μm | +10.6382 μm | 34.52 μm | 0.0259 μm |
| B07 | +23.7395 μm | +23.7906 μm | 30.76 μm | 0.0270 μm |

(65 台すべての rms。最大軌道変位の 0.08% 以内で解析式と一致する。)
- RFSW を ON にしても OFF と相対 4e-8 しか違わない (B16 の max abs DX は 26.2514467 → 26.2514476 μm)。10 μrad の蹴りでは、周長変化によるエネルギーのずれ (分散項) は見えない。

## 物理・解析式
上の式。ν は小数部を含めた全チューン (9.5975)。表の「解析式」は、蹴った磁石の中央の β と位相を使った薄い蹴りの近似。位置は PM (モニター) の点。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。保存してあった旧い出力は `PYSAD_COMPAT=k64` の出力と同じ (確認済み)。v2 は K64 のバグを直したので (`docs/47_v2_changes.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

○ 全 3 磁石 × 2 モードで全桁一致 (DX は 9 桁、max |DX| は 10 桁)。

## 未確認・限界
- 蹴りは 3 台だけ (28 台全部の応答行列は作っていない)。垂直方向 (`calc_respy`) は作っていない。
- 解析式は分散項 (周長変化) を含まない。RF ON でも差が出なかったので、この大きさの蹴りでは見えていない。
