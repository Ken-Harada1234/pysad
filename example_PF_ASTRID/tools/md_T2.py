# -*- coding: utf-8 -*-
import os
D = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'T2_samples')
PF = "加速器: PF 2.5 GeV リング (`lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。"


def w(name, text):
    open(os.path.join(D, name + '.md'), 'w', encoding='utf-8', newline='\n').write(text.replace('{PF}', PF).strip('\n') + '\n')


w('calc_respx_PF', r'''
# calc_respx (PF): 軌道応答行列と解析式
種別: スクリプト (サンプル集 calc_respx.sad と calc_resp_analy.sad を一つにした)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD` (RFSW は例の中で切り替える)。

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
○ 全 3 磁石 × 2 モードで全桁一致 (DX は 9 桁、max |DX| は 10 桁)。

## 未確認・限界
- 蹴りは 3 台だけ (28 台全部の応答行列は作っていない)。垂直方向 (`calc_respy`) は作っていない。
- 解析式は分散項 (周長変化) を含まない。RF ON でも差が出なかったので、この大きさの蹴りでは見えていない。
''')

w('calc_dispersion_PF', r'''
# calc_dispersion (PF): 分散関数
種別: スクリプト (サンプル集 calc_dispersion.sad)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
モニター `PM*` の位置の分散 `Twiss["EX","PM*"]` と s を印字する。さらに、(1) 分散の定義 $\eta = dx_{\rm co}/d\delta$ を、`DP0 = ±1e-4` の閉軌道の差分で確かめる、(2) 運動量コンパクション率 $\alpha_c = \frac1C\oint\frac{\eta}{\rho}ds$ を偏向磁石の和で確かめる。

## スクリプト
```
pmex = Twiss["EX","PM*"];  pml = Twiss["S","PM*"];  Do[Print[i, " ", pml[i], " ", pmex[i]], {i,1,Length[pml]}];
FFS["OFF RFSW"];  DP0 = 1E-4; FFS["cell calc"]; dxp = Twiss["DX","PM*"];  DP0 = -1E-4; ... dxm = ...;   ! (dxp-dxm)/2e-4 と EX を比べる
Σ ANGLE * (EX(入口) + 4 EX(中央) + EX(出口))/6 / C                                                       ! α_c
```
RF を切る (`OFF RFSW`) のは、RF があると運動量を変えても平衡エネルギーに戻ってしまうため。

## 結果の読み方
![EX at PM](calc_dispersion_PF_ex.png)
- 65 台のモニターの分散は、標準セル部で 0.28〜0.71 m、ID 直線部の近くで 0.001〜0.002 m まで下がる。
- 照合は次の表。

| 量 | 値 |
|---|---|
| max abs((DX(+δ)-DX(-δ))/2δ - EX) (65 台、δ=1e-4) | 2.1e-6 m |
| α_c = (1/C) Σ θ <η> (偏向磁石 28 台、Simpson) | 0.0065622 |
| α (emit の表) | 0.0065627 |

## 物理・解析式
$\eta_x = dx/d\delta$ と $\alpha_c = \frac1C\oint\frac{\eta_x}{\rho}ds = \frac1C\sum_{\rm bends}\theta\langle\eta\rangle$。どちらも数値と合う (差は 1e-5 の差分誤差と Simpson の誤差)。

## gs と pysad の比較
○ 全桁一致 (差分の max は 2.10852e-6 と 2.10848e-6 で、1e-10 の丸め)。pysad の emit の表は α を `.00656` と 3 桁で表示する (値は同じ)。

## 未確認・限界
- 垂直分散 (EY) は設計ラティスでは 0 なので見ていない。
- 周長に対する α の符号・定義 (`dl/L`) の違いは見ていない。
''')

w('calc_chromaticity_PF', r'''
# calc_chromaticity (PF): 自然色収差と六極の寄与
種別: スクリプト (サンプル集 calc_chromaticity.sad)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
`DP0` を -1%〜+2% (31 点) で動かしてチューンを `PolynomialFit[..., 6]` で当てはめ、c1 を色収差とする (サンプル集の書き方)。六極あり (設計値) と、`S* 0` の六極なし (自然色収差)。
さらに、四極と六極による色収差を、Twiss 関数からの和で計算して比べる。

## スクリプト
```
DP0 = (i-11)*0.001;  FFS["cell calc"];  tunex[i] = {DP0, Twiss["NX","$$$"]/2/Pi}; ...      ! i = 1..31 (-1%〜+2%)
Print[PolynomialFit[tunex,6]];  ... S* 0;  FFS["cell calc"];  同じ当てはめ
ChromSum[]  ! 四極: -ΣK1 <βx>/4π (x)、+ΣK1 <βy>/4π (y)。六極 (SFH,SDH): +ΣK2 ηx βx/4π (x)、-ΣK2 ηx βy/4π (y)。<β> = Simpson (入口・中央・出口)
```
範囲を **-1%〜+2% (平均が 0 でない) にした**のは、SAD の `PolynomialFit` が平均 0 の範囲では最高次以外の係数を壊すため (`T1_intro/fn_PolynomialFit_PF` 参照)。サンプル集の -1%〜+1% のままだと c1 が正しく出ない。

## 結果の読み方
| 量 | x | y |
|---|---|---|
| 四極の和 (式) | -13.4349 | -15.7741 |
| 六極の和 (式、設計値) | +13.3008 | +17.2620 |
| 四極 + 六極 (式) | -0.1341 | +1.4879 |
| ξ 差分 (±1e-4、設計の六極、T1 `FFS_FIT_chrom_PF`) | -0.0004 | -0.0013 |
| ξ 差分 (±1e-4、六極 0) | -13.2986 | -17.2514 |
| 当てはめの c1 (6 次、設計の六極) | -0.000434 | -0.001615 |
| 当てはめの c1 (6 次、六極 0) | -13.2983 | -16.979 (残差 5e-3、当てはまりが悪い) |

- 六極の寄与は 式 (13.3008 / 17.2620) と SAD の差 (六極あり-なしの ξ の差: 13.2979 / 17.2501) が 0.02% / 0.07% で一致する。
- 四極の和 (式) と自然色収差 (差分) の差は x で +0.14、y で -1.48 (四極の和が x は 1.0% 大きく、y は 8.6% 小さく出る)。この差は偏向磁石の寄与 (SAD の `CHRO` が数えないものと同じ) と見られる。
  概算では、偏向磁石の本体の収束 (1/ρ²) が x に 0.12、矩形磁石の端の収束 (tan e / ρ) が y に 2×0.80 の桁 (符号と係数は確かめていない) で、桁は合う。厳密な式は書いていない。
- 6 次当てはめは、自然色収差の y では当てはまりが悪い (残差 5.3e-3)。範囲 -1%〜+2% では 6 次でも足りない (高次の項が大きい)。

## 物理・解析式
$\xi_x = -\frac1{4\pi}\oint\beta_x(K_1-K_2\eta_x)\,ds$、$\xi_y = +\frac1{4\pi}\oint\beta_y(K_1-K_2\eta_x)\,ds$ (K1・K2 は SAD の積分強度で、符号は上の和のとおり。偏向磁石の項を除く)。

## gs と pysad の比較
○ c1〜c6 とも 9 桁一致 (数値は表のとおり)。

## 未確認・限界
- 偏向磁石の寄与は解析式で再現していない (上の概算のみ)。
- 六次の当てはめは、チューンの運動量依存が ±1% を超えると高次が大きく、c1 だけを見るべきで、c2 以上は壊れやすい。
''')

w('FFS_FIT_injbump_PF', r'''
# FFS FIT (PF): 入射バンプ (サンプル集「ローカルバンプ」)
種別: FFS コマンド `INS`・`FIT`・`FREE`・`GO` (サンプル集「ローカルバンプの例」)。加速器: PF 2.5 GeV (`lattice/v9_16_inj_noSF.sad`、入射点付近の六極の一部を外した版)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
入射点 (`INJ`) で、蓄積ビームの軌道を水平 +15 mm、角度 +1 mrad に置き、ライン終点で DX・DPX を 0 に戻す (バンプを閉じる)。キッカー `KK1`〜`KK4` の K0 (蹴り角) を自由変数にして `FIT` で解く。

## スクリプト
```
USE RING;  K* 0;  KK* 0;  cell calc emit;  SAVE ALL;     ! リングの周期解を求め、マークに保存
USE INJL;  ins calc;                                     ! 入射部分 (114 要素) を転送モードで
fit INJ DX 15/1000;  fit INJ DPX 1/1000;  fit DX 0;  fit DPX 0;    ! 場所を省略した fit = ライン終点
free KK1 K0;  free KK2 K0;  free KK3 K0;  free KK4 K0;  go;
```
- `ins` は MARK の値を初期条件にして転送するモード。`SAVE ALL` はベンドの K0 のように、要素の第 1 変数でないものも含めて保存する。
- キッカーの要素番号は 17、46、50、79、`INJ` は 48 番。**`KK3` (50) は INJ の下流**にあるが、バンプを閉じるために必要。

## 結果の読み方
| 量 | 値 |
|---|---|
| INJ の DX, DPX | 14.9999959 mm, 0.9999996 mrad (目標 15, 1) |
| ライン終点の DX, DPX | 1.3e-6 mm, -7.0e-7 mrad (目標 0) |
| KK1 K0 | -2.5211 mrad |
| KK2 K0 | -0.0095 mrad |
| KK3 K0 | -1.9202 mrad |
| KK4 K0 | -3.3750 mrad |

収束は `Matched. (1.4e-10)`。

## 物理・解析式
薄い蹴りが下流で作る軌道は $x(s) = -\theta\sqrt{\beta_k\beta_s}\sin(\psi_s-\psi_k)$ (SAD の K0 の符号が負の向きなので負号)。INS モードの Twiss の位相 `Twiss["NX",...]` を使い、各キッカーの中央の値で和を取った。

| 量 | 式 | SAD (FIT 後) |
|---|---|---|
| INJ の DX (INJ より上流のキッカー KK1・KK2 のみ) | 14.897 mm | 15.000 mm |
| ライン終点の DX (4 台全部) | 0.0047 mm | 0.0000013 mm |

式との差は INJ で 0.7%。キッカー (L=0.3 m) を中央の薄い蹴りに置き換えた近似と、KK2 (46 番) と INJ (48 番) の間の距離による。

## gs と pysad の比較
○ K0 は 9 桁一致 (-2.5210688818 / -2.5210688817 など)、DX・DPX の残りは 1e-10 の雑音の桁で違う。

## 未確認・限界
- 目標値 (15 mm, 1 mrad) とキッカー 4 台はサンプル集のとおり。実際の PF のキッカー値とは対応づけていない。
- バンプの内側の六極 (noSF) の影響は見ていない (ここは線形光学の例)。
''')

w('fn_DynamicApertureSurvey_DAsv_PF', r'''
# fn DynamicApertureSurvey (PF): sample_DAsv
種別: 関数 (サンプル集 sample_DAsv.sad)。{PF}追跡: `OFF RAD FLUC; ON COD RFSW RADCOD`。**2000 周** (サンプル集は 5000 周)。

## 何を見る例か
運動量方向 -2%〜+2% (0.5% おき、9 点)、水平 (y=1σ) と垂直 (x=1σ) の振幅 1〜100σ のサーベイ。`dp = 0.01/SIGE` (1% が何 σE か = 13.63)、σ = 水平エミッタンス × β のルート (垂直も水平エミッタンスを使うのがサンプル集の書き方)。

## スクリプト
```
dp = 0.01/SIGE;  zrange = Range[-2.0*dp, 2.01*dp, dp/2];
range = {{1,100},{1,1},zrange};  res = DynamicApertureSurvey[range, 2000, Output->6];
Print[Table[res[2,2,i,2], {i,1,Length[zrange]}]];      ! 運動量ごとのスコア
```
(サンプル集は `Print[res[2,2,1,2], res[2,2,2,2], ...]` と並べるが、`Table` で書き直した。)

## 結果の読み方 (スコア、(スコア-1)×100σ/50 が振幅 σ)
| dp/p | -2% | -1.5% | -1% | -0.5% | 0 | +0.5% | +1% | +1.5% | +2% |
|---|---|---|---|---|---|---|---|---|---|
| x スコア (gs, pysad) | 0 | 0 | 20 | 19 | 22 | 19 | 19 | 0 | 0 |
| y スコア (gs) | 0 | 0 | 11 | 13 | 14 | 14 | 13 | 0 | 0 |
| y スコア (pysad) | 0 | 0 | 13 | 13 | 14 | 14 | 13 | 0 | 0 |

運動量 0 で x のスコア 22 = (22-1)×2σ = 42σx = 27.0 mm、y のスコア 14 = (14-1)×2σ = 26σy = 10.8 mm (σx = 0.6433 mm、σy = 0.4154 mm)。±1.5% はスコア 0 (RF バケットの外、バケット半高さ 1.18%)。

## 物理・解析式
式の照合はなし。1000 周 (`T1_intro/fn_DynamicApertureSurvey_PF`) と比べると、x のスコアは 22 (同じ)、-1% で 21 → 20 と小さくなる。周回数 1000 → 2000 で DA はほぼ飽和している (放射なしの設計ラティス)。

## gs と pysad の比較
△ x は全点で一致。y の -1% (-13.6σE) の行だけ gs 11、pysad 13。DA の境界に近い粒子で、gs と pysad の周回数の端数の差がスコアに出た (他の 8 行は一致)。画面の記号 (`P`、`K` など) も境界付近で違う。

## 未確認・限界
- 5000 周でなく 2000 周。放射・量子励起は切ってある (サンプル集は `ON RAD FLUC`)。
- gs と pysad を合わせて約 14 秒。
''')

w('fn_DynamicApertureSurvey_da_xy_PF', r'''
# fn DynamicApertureSurvey (PF): da_xy (X-Y 面のダイナミックアパーチャ)
種別: 関数 (サンプル集 da_xy.sad)。{PF}追跡: `OFF RAD FLUC; ON COD RFSW RADCOD`。**24 方向 × 1000 周**、運動量 -2%・0・+2% の 3 点 (サンプル集は 161 方向、5000 周、運動量 ±3%)。

## 何を見る例か
X-Y 平面の原点から、角度 θ の方向に r = 1〜100σ の直線上を `DynamicApertureSurvey` で調べ、θ を 0〜360° で変えて、ダイナミックアパーチャの境界を (x, y) で描く。

## スクリプト
```
nstep = 24;  rmin = 1;  rmax = 100;  zrange = {-dp*2, 0, dp*2};
Do[ ang = (2*Pi)*((i+0.5)/nstep);  xmax = rmax*Cos[ang]; ymax = rmax*Sin[ang]; xmin = rmin*Cos[ang]; ymin = rmin*Sin[ang];
    res = DynamicApertureSurvey[{{xmin,xmax},{ymin,ymax},zrange}, nturn, Output->6];
    score[j] = Append[score[j], {(s-1)/50*xmax, (s-1)/50*ymax}] ... ,{i,0,nstep-1}];
```
- 角度に半ステップ (7.5°) のずらしを入れ、x 軸・y 軸の上の直線 (y=0 や x=0 で運動が平面に閉じる) を避けた。サンプル集は `rmin = 0` で原点を含む (SAD の版により原点を含む計算が止まらなくなる、とサンプル集に注意がある。最初に `rmin = 0`・15° 刻み (x 軸・y 軸上の方向を含む) で試したときは、出力の先頭の運動量の行の境界がすべて 2σ になった。原因の切り分けはしておらず、`rmin = 1` と半ステップのずらしに変えたら普通に出た)。
- サンプル集の「途中で止まると出力が整わない」ので、最後の `Print` で一括して出す形は同じ。

## 結果の読み方
![DA](fn_DynamicApertureSurvey_da_xy_PF_da.png)
- 運動量 0: 境界は水平側で 39.7σx (7.5° の方向)、垂直側で約 26σy (82.5° の方向)。原点からの距離で 24〜40σ。
- 運動量 ±2% (±27.3σE) は RF バケット (半高さ 1.18%) の外で、全方向で半径 2σ (最小の格子点 = 全滅)。

## 物理・解析式
式の照合はなし。

## gs と pysad の比較
○ 全方向の境界の値 (3 運動量 × 24 方向) が一致。

## 未確認・限界
- 24 方向・1000 周に短縮 (gs と pysad を合わせて約 47 秒)。境界の細かい構造 (共鳴のくぼみ) は 15° 刻みでは見えない。
- 放射・量子励起は切ってある。
''')

w('FFS_mom_dep_ts_PF', r'''
# mom_dep_ts (PF): 運動量依存チューンシフト
種別: スクリプト (サンプル集 mom_dep_ts.sad)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
`DP0` を -2%〜+2% (0.2% おき、21 点) で動かして、六極あり (設計値) と六極なし (`S* 0`) のチューンを出す。チュートリアルは 0.1% おき 41 点で、ここでは 21 点に間引いた。

## スクリプト
```
Do[ DP0 = (i-11)*0.002;  FFS["cell calc"];  Print[DP0, " ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi]; ,{i,1,21}];
DP = 0;  S* 0;  (同じループ)
```

## 結果の読み方
![tune vs dp](FFS_mom_dep_ts_PF_tune_vs_dp.png)

| dp/p | νx (六極あり) | νy (六極あり) | νx (六極 0) | νy (六極 0) |
|---|---|---|---|---|
| -2% | 9.63054 | 5.23520 | 9.89891 | 5.64483 |
| -1% | 9.60575 | 5.27126 | 9.73599 | 5.45487 |
| 0 | 9.59753 | 5.27896 | 9.59753 | 5.27896 |
| +1% | 9.60673 | 5.27197 | 9.46828 | 5.10787 |
| +2% | 9.63531 | 5.25003 | 9.34421 | 4.94688 |

- 六極ありでは、チューンは δ=0 の近くで極値になる (ξ ≈ 0) 放物線で、±2% でも νx は +0.03〜+0.04、νy は -0.03〜-0.04 の動き。
- 六極なしでは、ほぼ線形 (自然色収差が負) で、-2% で νx は +0.30、νy は +0.37、+2% で νx は -0.25、νy は -0.33 動く (曲がっている)。±1.6% 付近の νy (4.98 など) は整数 5 の近くで、隣の点と滑らかにつながらない。

## 物理・解析式
中央差分の傾き (±0.2%) は、六極なしで νx: (9.57110-9.62431)/0.004 = -13.30、νy: (5.24455-5.31357)/0.004 = -17.25。これは ±1e-4 の差分の ξ (-13.2986、-17.2514、`T1_intro/FFS_FIT_chrom_PF`) と 0.03% で一致する。四極と六極の寄与の和の式は `calc_chromaticity_PF`。

## gs と pysad の比較
○ 全 42 点のチューンが全桁一致 (表示の 10 桁)。

## 未確認・限界
- δ=0 の点で、六極ありの 1 本目 (9.5975304557) と設計チューン (9.5975304909) に 3.5e-8 の差がある (`disp` と `DP=0` の再設定の順序による収束の違いと見ているが、確認していない)。
- 整数・半整数共鳴の位置は、六極なしの ±1.6% 付近の νy が 5 に近い以外は確認していない (六極ありは ν ≈ 9.6、5.27 で整数・半整数に近づかない)。
''')

w('FFS_amp_dep_ts_PF', r'''
# amp_dep_ts (PF): 振幅依存チューンシフト
種別: スクリプト (サンプル集 amp_dep_ts.sad、アイディアは PF の担当者)。{PF}フラグ: `OFF RAD FLUC; ON COD RFSW RADCOD`。

## 何を見る例か
振幅 A σ のクーラン・シュナイダー不変量をもつ点を、x-x' 平面の楕円上に等間隔で置き (サンプル集は 36 点、ここでは **12 点**)、各点を初期条件にして `ins calc` でチューンを求めて平均する。水平方向 (5〜35σ、5σ おき) と垂直方向 (5〜30σ) を見る。

## スクリプト
```
inidx  = nx*Sqrt[w]*Sqrt[betax]*Cos[theta];   inidpx = nx*Sqrt[w]/Sqrt[betax]*(Sin[theta]-alphax*Cos[theta]);   ! w = EMITX、σ = sqrt(w β)
LINE["DX",ininame] = inidx; LINE["DPX",ininame] = inidpx; ...   ! ライン先頭の MARK に初期座標をセット
FFS["ins calc"];  {Twiss["NX","$$$"]/2/Pi, Twiss["NY","$$$"]/2/Pi}
```
- `FFS["cell calc emit fix"]` で始点の α・β を固定して `SetInitialParameters[]` に取る。
- 振幅は始点 (MKSSES、βx 11.98 m、βy 5.00 m) の σ (σx 0.643 mm、σy 0.415 mm)。
- 他の方向は 0.01σ (y = 0 の厳密な平面運動を避けるため)。

## 結果の読み方
![dnu vs A^2](FFS_amp_dep_ts_PF_dnu_vs_A2.png)

| 水平振幅 | Δνx | Δνy | Δνx / A² |
|---|---|---|---|
| 5σ | +0.000819 | +0.000255 | 3.28e-5 |
| 10σ | +0.003285 | +0.001063 | 3.29e-5 |
| 20σ | +0.013301 | +0.004764 | 3.33e-5 |
| 30σ | +0.031169 | +0.011699 | 3.46e-5 |
| 35σ | +0.044158 | +0.016342 | 3.60e-5 |

| 垂直振幅 | Δνx | Δνy | Δνy / A² |
|---|---|---|---|
| 5σ | +0.000764 | +0.007060 | 2.82e-4 |
| 10σ | +0.003046 | +0.027905 | 2.79e-4 |
| 20σ | +0.011666 | +0.09997 | 2.50e-4 |
| 30σ | +0.025465 | +0.18661 | 2.07e-4 |

- Δν は小振幅で A² に比例 (Δνx/A² がほぼ一定 = 3.3e-5/σ²) し、振幅が大きいと高次の項が見える。垂直振幅の Δνy/A² は 5σ で 2.8e-4、30σ で 2.07e-4 に減る。
- 水平 35σ (≈ 22 mm) は、DA (約 42σ、`T1_intro/fn_DynamicApertureSurvey_PF`) の近くで、チューンが 0.044 動く。垂直振幅のチューンシフトは水平より約 8 倍大きい。
- 試しに 100σ を入れると DA の外で、チューンが gs は 1.6e199、pysad は NaN になった (DA を超えた振幅の `ins calc` は意味がない。この例の最終版には入れていない)。

## 物理・解析式
式の照合は省く (非線形の摂動論の 1 次の式は書いていない)。小振幅の A² 比例は表のとおり。

## gs と pysad の比較
○ 全 13 点 (×2 チューン) が全桁一致。

## 未確認・限界
- 12 点平均 (サンプル集は 36 点)、水平 7 点・垂直 6 点 (サンプル集は 10σ おき 20 点)。
- 位相空間の平均を取るので、円状に近い運動の仮定の平均 (SAD の `ins calc` の初期条件は線形楕円で、非線形では楕円は変形する)。
- 表の Δν の一部の値 (20σ 垂直など) は丸めてある。
''')

w('calc_RadInt_PF', r'''
# calc_RadInt (PF): 放射積分と平衡エミッタンス
種別: スクリプト (サンプル集 calc_RadInt.sad)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
偏向磁石 28 台の放射積分 $I_1\ldots I_5$ を Twiss 関数から数値積分し、$\epsilon_x$・$U_0$・$\sigma_E$・減衰分配数・運動量コンパクション率を、`emit` の出力 (SAD が 6 次元で求める値) と照合する。

## スクリプト
サンプル集は、勾配つき偏向磁石の H 関数の積分の閉じた式 (半分ずつ、中央のツイスから) を使う。PF の偏向磁石は勾配なし (K1 = 0) なので、ここでは各磁石の中を 17 点に分けて `Twiss["EX", n+f]` (f = 0〜1) を読み、合成 Simpson で積分した。
```
sa += w*t[3];  sh += w*(((1+t[2]^2)/t[1])*t[3]^2 + 2*t[2]*t[3]*t[4] + t[1]*t[4]^2);   ! ∫η ds、∫H ds (t = {BX,AX,EX,EPX})
I1 += sa/ρ;  I2 += L/ρ²;  I3 += L/|ρ|³;  I4 += sa/ρ³ - 2(η0 tan e1 + η1 tan e2)/ρ²;  I5 += sh/|ρ|³
```
`I4` の最後の項は矩形磁石の端 (E1 = E2 = 0.5 ANGLE) の収束で、強さ $-\tan e/\rho$ の薄い四極として $2k\eta/\rho$ に入る。

## 結果の読み方
| 量 | 式 (数値積分) | SAD (emit) |
|---|---|---|
| I1 | 1.22742 m | |
| I2 | 0.725766 /m | |
| I3 | 0.0838328 /m² | |
| I4 | -0.0224985 /m | |
| I5 | 0.00281098 /m | |
| α_c = I1/C | 0.0065625 | 0.0065627 |
| U0 = Cγ E⁴ I2/(2π) | 0.399139 MeV | 0.3991454 MV |
| Jx = 1 - I4/I2 | 1.0310 | 1.0313 |
| σE = √(Cq γ² I3/(Jz I2)) | 7.3353e-4 | 7.3366e-4 |
| εx = Cq γ² I5/(Jx I2) | 34.456 nm rad | 34.538 nm rad |

($C_q = 3.832\times10^{-13}$ m、$C_\gamma = 8.846\times10^{-5}$ m/GeV³、γ = 4892.38。)

## 物理・解析式
$$\epsilon_x = C_q\gamma^2\frac{I_5}{J_x I_2},\quad J_x = 1-\frac{I_4}{I_2},\quad \sigma_E^2 = C_q\gamma^2\frac{I_3}{J_z I_2},\quad U_0 = \frac{C_\gamma E^4}{2\pi}I_2,\quad \alpha_c = \frac{I_1}{C}$$
εx は 0.24%、Jx は 0.03%、σE は 0.02%、U0 は 0.002% で SAD の値と一致する。

## gs と pysad の比較
○ 全桁一致。

## 未確認・限界
- サンプル集の閉じた式 (半磁石ごと) を最初に PF に当てはめたところ、I5 が 0.00736 で εx が SAD の約 2.8 倍になった。式の転記か適用の誤りと見て、原因は調べずに、数値積分に置き換えた (サンプル集の式そのものは検証していない)。
- `I4` の端の項の符号と係数 2 は上のとおりに導いたが、Jx が SAD の減衰分配数 (1.0313) と 3 桁一致することで確認した (係数 1 だと 1.0042 で合わない)。
- 勾配つき偏向磁石、垂直偏向、ウィグラーは含まない。
''')

w('phase_PF', r'''
# phase (PF): 位相空間
種別: スクリプト (サンプル集 phase.sad)。{PF}追跡: `OFF RAD FLUC; ON COD RFSW RADCOD`。**60 周** (サンプル集は 1000 周)。

## 何を見る例か
リングの始点 (MKSSES) で 1 粒子を、初期振幅 1σ・10σ・20σ・30σ (x のみ、他の座標は 1e-9) で 60 周追跡し、`TrackParticles[res, 1]` (1→1 = リング 1 周) を繰り返して毎周の (x, x') を取る。

## スクリプト
```
res = {1, {x1, {1E-9},{1E-9},{1E-9},{1E-9},{1E-9},{1}}};
Do[ res = TrackParticles[res, 1];  x = res[2,1,1]; px = res[2,2,1]; Print[...] ,{i,1,60}];
```
`res` は `{位置, 座標}` なのでそのまま次の呼び出しに渡せる。クーラン・シュナイダーの不変量 $\gamma x^2 + 2\alpha x x' + \beta x'^2$ を EMITX で割った値も印字する。

## 結果の読み方
![x-px](phase_PF_xpx.png)
- 周回ごとの点が x-x' 平面の楕円の上を回る (1 周の位相進み 0.5975×2π)。振幅が大きいと歪む。

| 初期振幅 | 不変量/EMITX (1 周目) | 60 周目 |
|---|---|---|
| 1σ | 1.0090 | 1.0033 |
| 10σ | 108.56 | 108.51 |
| 20σ | 462.8 | 349.0 |
| 30σ | 1086 | 1052 |

線形なら振幅² (1, 100, 400, 900) で一定になる。1σ ではほぼ 1 (1.0〜1.01 で保存) で、大振幅ほど (20σ 以上で) 値が変わる。1 周目の値が振幅² から外れる (10σ: 108.6) のは、始点の分散・閉軌道 (6 次元) と非線形項 (六極) による。

## 物理・解析式
線形運動では $\gamma x^2 + 2\alpha x x' + \beta x'^2$ は保存 (Courant-Snyder)、点は楕円 $\gamma x^2 + 2\alpha x x'+\beta x'^2 = \epsilon$ の上を回る。1σ でこれが 0.3% で保存されていることが上の表。

## gs と pysad の比較
○ 全 240 点の (x, x', y, z) が 10 桁一致 (z の値だけ最後の桁が 1 違う)。

## 未確認・限界
- 60 周 4 粒子のみ (サンプル集は 1000 周、粒子 8)。y・z 方向の位相空間 (サンプル集の 3 枚の図) は印字だけで、図にしていない。
- 振幅 20σ・30σ の周回ごとの不変量のゆらぎは、非線形共鳴と見ているが、調べていない。
''')

w('FFS_draw_topdrawer_PF', r'''
# FFS draw (PF): topdrawer 形式の出力
種別: FFS コマンド `draw` (入門マニュアル その 4、サンプル集「topdrawer 形式の図の作り方」)。{PF}フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
チュートリアルでは、`FFS["draw bx & by & ex b*",12]` と書くと装置番号 12 (= `fort.12`) に topdrawer 形式の図 (β、分散、ラティスの絵) が書き出され、`td -b -d postscr fort.12` で PostScript にする、とある。
この環境の SAD (K64 版、WSL) では `draw` が topdrawer の出力として動かないので (次の節)、**topdrawer 形式のファイルを台本から手書きで作る**例にした (`Print` の代わりに `Write[12, ...]` で `fort.12` に書ける)。

## スクリプト
```
Print["NEWFRAME;SET CARD 80;SET FONT DUPLEX"];  Print["SET LIMIT X 0 ", smax];  Print["TITLE 7.3 9.3 CENTER SIZE 2.5 '...'"];
Do[ Print[Twiss["S",i], " ", Sqrt[Twiss["BX",i]]], {i,1,ne}];   Print["JOIN 1 TEXT 'sqrt(BX)'"];  Print["NEW PLOT"]; ... (BY、EX も同様)
```
topdrawer の `NEWFRAME`、`SET LIMIT`、`TITLE`、`JOIN 1 TEXT` は topdrawer のコマンド (データ行の `x y` を並べて `JOIN` で折れ線にする)。

## 結果の読み方
- 画面に `TDBEGIN`〜`TDEND` の間に 2277 行が出る。それを取り出したものが `FFS_draw_topdrawer_PF_fort12.td` (`tools/make_figs.py` が `.gs.txt` から作る)。
- ![hand-written td](FFS_draw_topdrawer_PF_plot.png)
  (上の `.td` の内容を matplotlib で描き直した図。`td` は入っていないので `.td` の描画は試していない。)

## gs と pysad の比較
- 手書きの出力 (`Print`): ○ 全桁一致。
- チュートリアルの書き方の `draw` コマンド:
  - gs: `draw bx & by & ex b*;` は表示のない環境では何も出さずに通る。`$DisplayFunction = TopDrawer;` を付けると戻らなくなる (25 秒待って止まったまま)。`FFS["draw bx",12];` は `FFSCanvasDraw` の内部エラーで落ちた。この K64 版の `draw` は Tk のキャンバス (`FFSCanvasDraw`) に描く実装に置き換わっていて、topdrawer への出力は `OpticsPlot` などの `TopDrawer` 関数 (`Packages/TopDrawer.n`、`fort.9` に書く) が受け持つと見ている。
  - pysad: `draw` を未対応として無視し (メッセージ `FFS command DRAW は未対応`)、何も出さない。
  - 試した小さい台本は `T1_intro/FFS_draw_OpticsPlot_PF.sad` に `draw bx & by & ex b*;` の 1 行として残してある (走る)。

## 未確認・限界
- 実際の `td`/`tdr` コマンドで描いた図は作っていない (入っていない)。topdrawer の命令 (`SET LIMIT`、`JOIN` など) の意味は SADHelp の `TopDrawer` の節と `Packages/TopDrawer.n` の書き方に合わせたが、`td` で描画できるかは確認していない。
- 図の BX・BY は 756 要素の値の折れ線 (要素の中は直線補間)。
''')


def readme():
    txt = r'''
# T2_samples (PF の話題): サンプル集
元: KEK SAD チュートリアル「SAD サンプル集」(Tutorial2) の PF のスクリプト。PF 2.5 GeV リング (`../lattice/v9_16_asym.sad`、入射バンプだけ `v9_16_inj_noSF.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。ASTRID2 の話題 (チューンサーベイ・MPW など) は別の担当の `*_ASTRID2.*` (このフォルダの別のファイル)。

| ファイル (.sad/.gs.txt/.pysad.txt/.md) | サンプル集の節 | 内容 | gs と pysad |
|---|---|---|---|
| `calc_respx_PF` | calc_respx / calc_resp_analy | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| `calc_dispersion_PF` | calc_dispersion | 分散 + 差分と α_c の確認 | ○ |
| `calc_chromaticity_PF` | calc_chromaticity | `PolynomialFit` の 6 次、自然色収差、四極・六極の和の式 | ○ |
| `FFS_FIT_injbump_PF` | ローカルバンプ | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| `fn_DynamicApertureSurvey_DAsv_PF` | sample_DAsv | DA サーベイ (x・y 方向、運動量 9 点、2000 周) | △ (y の 1 行だけ gs 11 / pysad 13) |
| `fn_DynamicApertureSurvey_da_xy_PF` | da_xy | X-Y 面の DA (24 方向、1000 周) と図 | ○ |
| `FFS_mom_dep_ts_PF` | mom_dep_ts | 運動量依存チューンシフト (六極あり・なし) と図 | ○ |
| `FFS_amp_dep_ts_PF` | amp_dep_ts | 振幅依存チューンシフト (12 点平均) と図 | ○ |
| `calc_RadInt_PF` | calc_RadInt | 放射積分、εx・U0・σE・Jx の式と数値の照合 | ○ (εx は SAD と 0.24%) |
| `phase_PF` | phase | 位相空間 (60 周) と図 | ○ |
| `FFS_draw_topdrawer_PF` | topdrawer 形式の図 | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |

- 走らせ方: `py -3.12 tools/runpf.py T2_samples/NAME.sad` (全部を回すときは `*_PF.sad` だけを指定すること。このフォルダは ASTRID2 の担当と共有)。`.sad` は `tools/build_T2.py`、図は `tools/make_figs.py`。
- 長い計算 (DA) は 5000 周 → 1000〜2000 周、方向 161 → 24 などに短縮した (各 `.md` に書いてある)。
- 未対応・書けなかったもの: `calc_respy` (垂直の応答)、`drawmagvw3d` (CAD 出力、CAD に読み込む `.vss` は作っていない)、`modules_*` の GUI 版。
'''
    open(os.path.join(D, 'README_PF.md'), 'w', encoding='utf-8', newline='\n').write(txt.strip('\n') + '\n')


readme()
