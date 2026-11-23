# -*- coding: utf-8 -*-
import os
D = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'T1_intro')


def w(name, text):
    open(os.path.join(D, name + '.md'), 'w', encoding='utf-8', newline='\n').write(text.strip('\n') + '\n')


COMMON = "加速器: PF 2.5 GeV リング (`lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。"

w('FFS_cell_calc_emit_PF', '''
# FFS cell / calc / emit (PF): 実行の仕方と最初の三つのコマンド
種別: FFS コマンド (チュートリアル「SAD 入門マニュアル」その 1・その 2)。''' + COMMON + '''フラグ: `ON LOG; ON RAD FLUC COD RFSW RADCOD;`

## 何を見る例か
ラティス台本の後ろに計算台本を付けて gs に渡す形 (`.sad` 一つ)、`USE`・`cell`・`calc`・`emit`・`SAVE ALL` の意味、周長・チューン・エミッタンスを読み出す。
gs への渡し方は `gs file.sad < /dev/null` (チュートリアルは `cat` でラティスと計算台本を 1 ファイルにして渡す)。トップレベル (要素・ライン定義) と FFS (`FFS;` 以降) の二段で、`stop` 1 回目で FFS が、2 回目で全体が終わる (ここでは `Exit;`)。

## スクリプト
```
ON LOG;                              ! 入力の復唱をしない
ON RAD FLUC COD RFSW RADCOD;         ! 放射・量子励起・COD・RF・放射による軌道
FFS;
$FORM = "15.10"; PageWidth = 1999;   ! 数値の桁と、行の折り返し無し
USE RING;  cell calc emit;           ! cell = 周期解、calc = 軌道と光学、emit = 平衡パラメータ
SAVE ALL;
Print["circumference [m] = ", LINE["S","$$$"]];  ...
```
- `cell` は周期解を求めるモードのフラグ (外すと `INS` = MARK の値から転送するモード)。`calc` は閉軌道のあと周期ツイス、`emit` は 6 次元 (RF 入り) の平衡量。
- コマンド (`cell` など) は大文字小文字を区別しない。関数 (`Print[]`、`Twiss[]`) は区別し、`[]` が付く。
- `$$$` はライン終点、`^^^` は始点。

## 結果の読み方
| 量 | 値 | 出所 |
|---|---|---|
| 周長 | 187.0354880 m | `LINE["S","$$$"]` |
| 要素数 | 756 | `LINE["LENGTH"]` |
| チューン | 9.5975304909 / 5.2789602116 | `Twiss["NX","$$$"]/2/Pi` |
| εx | 3.4538e-8 m rad | `EMITX` (`emit` の表: Emittance X) |
| エネルギー拡がり | 7.3366e-4 | `SIGE` |
| バンチ長 | 9.765 mm | emit 表 |
| U0 | 0.3991 MV | emit 表 |
| 運動量コンパクション率 | 6.563e-3 | emit 表 |
| 有効ハーモニック数 | 311.94 | emit 表 (RF の FREQ から求めた値で整数ではない。`v9_16_inj_noSF` は `HARM = 312`) |
| 減衰時間 | x 7.58 ms、y 7.82 ms、z 3.97 ms | emit 表 |
| 減衰分配 | Jx 1.0313、Jy 0.9999、Jz 1.9688 | emit 表 |

## 物理・解析式
放射積分からの εx・U0・α・Jx の照合は `T2_samples/calc_RadInt_PF` (式と数値の表)。

## gs と pysad の比較
○ 全項目が表示桁内で一致。差は εy の浮動小数点雑音だけ (gs 7e-38、pysad -3e-43 で、どちらも 0)。

## 未確認・限界
- 偏向磁石 28 台・誤差なしの設計ラティス。入射点 `INJ` や ID (`U*`) は MARK で、K 値は 0。
- 垂直エミッタンスは結合なしで 0。
''')

w('FFS_disp_type_PF', '''
# FFS disp / type (PF): 表示コマンド
種別: FFS コマンド (入門マニュアル その 2)。''' + COMMON + '''

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
○ 数値は全桁一致。`disp` の最後の 3 列 (EY・EPY・DetR) だけ、gs は 1e-16 の雑音を出し pysad は 0 を出す。

## 未確認・限界
`disp` の列の意味は SADHelp とチュートリアルの説明による。
''')

w('fn_Twiss_LINE_PF', '''
# fn Twiss / LINE (PF): ツイスと要素の読み書き
種別: 関数 (入門マニュアル その 2)。''' + COMMON + '''

## 何を見る例か
`Twiss[量, 場所]` と `LINE[量, 要素]` の使い方。場所は要素名 (ワイルドカード可)、ライン上の番号、小数 (要素の中)。

## スクリプト
```
Twiss["BX","QF"]                  ! 名前が正確に一致するときは最初の 1 個 (QF.1)
Twiss[{"BX","BY"},"INJ"]
Length[Twiss["BX","QF*"]]         ! QF* で始まる名前は全部 (QFA・QFB を含み 14 個)
Twiss[{"BX","BY"},{n, n+0.5, n+1}]   ! 入口・中央・出口 (n は LINE["POSITION","QF"])
LINE[{"POSITION","NAME","ELEMENT","S","L","K1"},"QF"]
LINE["GEO", LINE["LENGTH"]]       ! 最後の要素の幾何学座標
LINE["K1","QF"] = 0.75;           ! 代入 (再計算は cell calc)
```

## 結果の読み方
- `Twiss["BX","QF"]` = 5.1495067164 m (QF.1 の入口)。`Twiss["BX","QF*"]` の長さは 14。
- QF.1 の入口・中央・出口: βx = 5.1495, 5.7313, 5.4930 m。βy = 7.9238, 6.9165, 6.9968 m。
- `LINE[...,"QF"]` = {149, "QF.1", "QF", 38.44848, 0.4, 0.7483195263}。
- リングを一周した最後の要素の GEO は (x, y) = (-1.2e-14, -9.6e-15) m で、閉じている (幾何学が一周で元に戻る)。

## 物理・解析式
周長は `LINE["S","$$$"]` = Σ 要素長。他に式の照合はなし。

## gs と pysad の比較
○ 全桁一致 (GEO の 1e-15 の雑音の桁だけ違う)。

## 未確認・限界
`LINE` の代入は同名のすべての要素 (`QF` の 14 個ではなく QF という名の全部) に効くと見ているが、個別の要素名 (`QF.1`) での代入との違いは確かめていない。
''')

w('fn_PolynomialFit_PF', '''
# fn PolynomialFit (PF): チューンの運動量依存を多項式で当てはめる
種別: 関数 + 特殊変数 `DP0` (入門マニュアル その 3、サンプル集 calc_chromaticity)。''' + COMMON + '''

## 何を見る例か
`DP0` (設計運動量からのずれ) を -2%〜+2% に 0.1% 刻みで 41 点動かし、チューンを `Table[]` のリストに入れて `PolynomialFit[リスト, 次数]` で当てはめる。
`PolynomialFit` は `{{c0, c1, c2, c3}, {Residual -> r}}` を返す ($\\nu = c_0 + c_1\\delta + c_2\\delta^2 + c_3\\delta^3$)。

## スクリプト
```
tunex = Table[{0,0},{41}];
Do[ DP0 = (i-21)*0.001;  FFS["cell calc"];  tunex[i] = {DP0, Twiss["NX","$$$"]/2/Pi}; ...,{i,1,41,1}];
DP0 = 0;   Print[PolynomialFit[tunex,3]];
```
`Do[]` の中でコマンドを使うときは `FFS["cell calc"]` と包む。`DP0 = 0;` に戻し忘れると以後の計算が狂う。

## 結果の読み方
![tune vs dp](fn_PolynomialFit_PF_tune_vs_dp.png)
- チューンは δ=0 の近くで極値 (6 極補正済みなので ξ≈0)、δ が大きいと高次の依存が見える。
- **対称な範囲 (平均 0) の当てはめでは、最高次以外の係数が壊れる** (SAD の `Packages/Fit.n` の癖。pysad の `docs/05_sadscript.md` にも同じ記述。gs も pysad も同じ値を出す)。
  下の表のように、-2%〜+2% の 3 次当てはめは c3 (νx 257.84、νy 1049.07) だけが numpy と一致し、c0〜c2 は別の値 (c0 = 355.6 は明らかに誤り) になる。平均が 0 でない範囲 (-2%〜+1%、31 点) では全係数が numpy と一致する。

| 当てはめ | 範囲 | c0 | c1 | c2 | c3 |
|---|---|---|---|---|---|
| gs (SAD) νx | -2%〜+2% | 355.600 | 949.82 | 861.66 | 257.84 |
| numpy 3 次 νx | -2%〜+2% | 9.59747 | 0.02658 | 88.136 | 257.84 |
| gs (SAD) νx | -2%〜+1% | 9.597557 | 0.03070 | 85.837 | 129.21 |
| numpy 3 次 νx | -2%〜+1% | 9.597557 | 0.03070 | 85.837 | 129.21 |
| gs (SAD) νy | -2%〜+2% | 966.39 | 2971.4 | 3059.3 | 1049.07 |
| numpy 3 次 νy | -2%〜+2% | 5.27981 | -0.08026 | -87.883 | 1049.07 |
| gs (SAD) νy | -2%〜+1% | 5.278889 | -0.16262 | -67.425 | 2328.9 |
| numpy 3 次 νy | -2%〜+1% | 5.278889 | -0.16262 | -67.425 | 2328.9 |

線形の色収差 ξ は c1 だが、この 3 次の c1 は広い範囲での当てはめで、±1e-4 の差分 (`FFS_FIT_chrom_PF`) の ξ ≈ -4e-4 / -1e-3 とは一致しない (3 次式が広い範囲を平均するため)。この PF は6 極で ξ ≈ 0 に補正済み。自然色収差は `T2_samples/calc_chromaticity_PF`。

## 物理・解析式
$\\xi = d\\nu/d\\delta$ ($\\delta = \\Delta p/p$) は c1。次数が高いほど広い範囲を追えるが、範囲の端で振動する。表の「numpy」は、gs が出した点を `numpy.polyfit` で当てはめた検算。

## gs と pysad の比較
○ 41 点のチューンは全桁一致、係数は 10 桁一致 (壊れる癖も同じ)。

## 未確認・限界
- 癖の原因は `Fit.n` の平均 `xa=0` の分岐 (pysad の説明)。ソースの確認はしていない。チュートリアルの NSLS の例 (c1 = -22) には、この癖は現れていない。
- 次数を上げたときの当てはめの残差は `Residual` の値だけを見た。
''')

w('FFS_draw_OpticsPlot_PF', '''
# FFS draw / OpticsPlot (PF): オプティクスの図
種別: FFS コマンド `draw` と関数 `OpticsPlot[]` (入門マニュアル その 4)。''' + COMMON + '''

## 何を見る例か
- `draw bx & by & ex b*`: topdrawer 形式でβ・分散とラティスの絵を描くコマンド (`&` は別の枠、末尾の `b*` は絵の中に名前を出す要素)。
- `OpticsPlot[{{"BX","BY"},{"EX"}}]`: Tk のキャンバスに同じ図を描く。窓を出す書き方は
  `w = KBMainFrame["sample", namae, Title->"Optics Plot!"]; c1 = Canvas[namae, Width->640, Height->400]; $DisplayFunction = CanvasDrawer; Canvas$Widget = c1; pl = OpticsPlot[...]; TkWait[];`

## スクリプト
窓の出ない環境 (gs を `< /dev/null`、WSL に表示なし) で走るので、`OpticsPlot[..., DisplayFunction->Identity]` で `Graphics` の式だけ受け取り、図の元の数値 (s, βx, βy, ηx) を全要素について印字する (図は matplotlib で描いた)。

## 結果の読み方
![optics](FFS_draw_OpticsPlot_PF_optics.png)
- `OpticsPlot returned: Graphics with 2 parts` (Graphics 式が返る)。
- 756 要素の値の範囲は、βx 0.94〜20.2 m、βy 0.41〜29.0 m、ηx 0〜0.71 m。直線部 (ID 部) で βy が小さく、標準セルでは分散が 0.3〜0.7 m。

## gs と pysad の比較
○ 数値と戻りの式の形は一致。
- `draw` コマンド: gs は表示のない環境では何も出さずに通る。`$DisplayFunction = TopDrawer;` を付けると、この環境では gs が戻らなかった (25 秒待っても止まったまま) ので付けていない。`FFS["draw bx",12]` は K64 版 gs では内部エラー (`FFSCanvasDraw` の中) で落ちた。pysad は `draw` を未対応として無視し、メッセージを出す。
- topdrawer 形式のファイルを手書きで作る例は `T2_samples/FFS_draw_topdrawer_PF`。

## 未確認・限界
- Tk の窓を実際に開くところは確認していない (gs は表示環境なし、pysad の窓の PNG は T5 の担当)。
- 軸ラベルの文字化け (チュートリアルの「既知の問題」) は topdrawer/tdr の話で、ここでは再現できない。
''')

w('FFS_FIT_chrom_PF', '''
# FFS FIT (PF): 色収差補正の6 極フィット
種別: FFS コマンド `FIT`・`FREE`・`GO`・`FITP`・`DP` (入門マニュアル その 5)。''' + COMMON + '''

## 何を見る例か
6 極 (SFH・SDH の 2 ファミリー) を 0 にして自然色収差を出し、±DP の 2 点でチューンが設計値になるように `FREE S*` でフィットし直す。設計の6 極値が戻るかを見る。

## スクリプト
```
S* 0;  cell calc;                  ! 6 極なし (自然色収差)
fit nx nx0 2;  fit ny ny0 2;       ! 場所は省略 = リング終点。2 = DP の 2 点 (-DP, +DP)
DP = 0.0001;                       ! 運動量の幅 (既定 1%)
free S*;  go;
FITP 1;                            ! 点数を 1 に戻す (戻さないと次の fit も 2 点になる)
type SFH;  type SDH;  SAVE S*;
```
`ChromFD[]` は `DP0 = ±1e-4` で `cell calc` し $(\\nu_+ - \\nu_-)/(2\\times10^{-4})$ を返す関数 (ξ の確認用)。

## 結果の読み方
| | ξx | ξy | SFH K2 | SDH K2 |
|---|---|---|---|---|
| 設計の6 極 | -4.2e-4 | -1.26e-3 | 3.960643 | -3.917196 |
| 6 極 0 (自然色収差) | -13.2986 | -17.2514 | 0 | 0 |
| フィット後 | -3.5e-6 | -2.1e-5 | 3.960795 | -3.917392 |

- 自然色収差 (-13.3, -17.3) が 2 ファミリーの6 極で ±1e-4 の 2 点において打ち消され、**元の設計の6 極 (3.9606, -3.9172) とフィット値 (3.9608, -3.9174) が 4 桁で一致**する。このラティスの6 極は「線形の色収差ゼロ」の条件で決まっていると読める。
- 収束は `Matched. (8.2e-12)`。

## 物理・解析式
$\\xi_x = -\\frac1{4\\pi}\\oint\\beta_x K_1 ds + \\frac1{4\\pi}\\oint\\beta_x\\eta_x K_2\\, ds$ (K1・K2 は SAD の積分強度)。四極と6 極の和の照合は `T2_samples/calc_chromaticity_PF`。

## gs と pysad の比較
○ 自然色収差の ξ は全桁一致。フィット後の6 極は gs 3.960795117 / pysad 3.960795576 で 7 桁一致 (FIT の収束の途中経過の違い)。フィット後の ξ の残り (gs -3.5e-6、pysad -1.2e-6) はどちらも 0 とみなせる。

## 未確認・限界
- DP=1e-4 の 2 点だけを合わせるので、高次 (δ² 以上) の色収差は補正されない (`fn_PolynomialFit_PF` の図)。
- FIT の途中表示 (残差) の行数は gs の方が多い (表示のみ)。
''')

w('FFS_FIT_modules_PF', '''
# FFS FIT (PF): チューンのマッチングと色収差補正のモジュール化
種別: FFS コマンド `FIT` + ユーザー定義関数 (入門マニュアル その 6)。''' + COMMON + '''標準セル `NCELL` (B20 を含む 1 セル、MKNCES から)。

## 何を見る例か
`ChromCorrect[]`・`MatchingNcellTune[]` という関数にして、標準セルのチューンを 2x2 の格子点にマッチ → 色収差補正、を `Do` で回す。関数の中の FFS コマンドは `FFS["..."]` で包む。

## スクリプト
```
MatchingNcellTune[]:=( FFS["cell calc;"]; FFS["fit NX tunex;"]; FFS["fit NY tuney;"]; FFS["free QD QF;"]; FFS["go;"]; FFS["cell calc"]; );
ChromCorrect[]:=( FFS["cell calc;"]; {nx0, ny0} = Twiss[{"NX","NY"},"***"]/2/Pi; FFS["FIT;"]; FFS["NX nx0 2;"]; ... FFS["FREE SFH SDH;"]; FFS["GO;"]; FFS["FITP 1;"]; );
Do[ Do[ FFS["USE NCELL;"]; MatchingNcellTune[]; ChromCorrect[]; Print[...]; ,{tuney,0.15,0.20,0.05}], {tunex,0.30,0.35,0.05}];
```
- マッチングのたびに `USE NCELL` する (SAVE しないので、四極・6 極の値が元に戻る。発散しても引きずらない)。
- 引数のない関数は大域変数 (`tunex`・`tuney`) を直接使う。値を返す関数は最後の式に `;` を付けない。

## 結果の読み方
元の標準セルのチューンは (0.27076, 0.06247) (1 セルあたり)。目標へ四極 (QF, QD の K1) と6 極 (SFH, SDH の K2) を動かした結果 (gs):

| 目標 (νx, νy) | 達成 | QF K1 | QD K1 | SFH K2 | SDH K2 |
|---|---|---|---|---|---|
| 0.30, 0.15 | 0.300000, 0.149996 | 0.84545 | -0.60989 | 1.5844 | -1.6968 |
| 0.30, 0.20 | 0.300000, 0.200000 | 0.86373 | -0.69321 | 1.5750 | -1.8364 |
| 0.35, 0.15 | 0.350001, 0.149999 | 0.94534 | -0.63985 | 2.1992 | -2.2731 |
| 0.35, 0.20 | 0.350000, 0.199993 | 0.95801 | -0.71745 | 2.1698 | -2.4234 |

達成チューンの 1e-5 程度のずれは、色収差補正 (±DP のフィット) が四極を固定したまま6 極でチューンを少し変えるため、と見ている (確認はしていない)。

## gs と pysad の比較
○ 四極の K1 は 7〜8 桁、6 極の K2 は 6〜7 桁で一致。達成チューンは 1e-8〜1e-5 のずれ (FIT の収束経路の差)。

## 未確認・限界
- 標準セルの `emit` は RF が無いので NaN になる (`cell calc` だけを使っている)。チュートリアルの `EMITX` の印字は外した。
- 0.30〜0.35 は PF の標準セルの運転点ではなく、台本の書き方を示すための格子。
- チュートリアルの `fit Q2 BX 15` と `FitValue` で上下限を付ける書き方は試していない。
''')

w('fn_DynamicApertureSurvey_PF', '''
# fn DynamicApertureSurvey (PF): ダイナミックアパーチャサーベイ
種別: 関数 (入門マニュアル その 7)。''' + COMMON + '''追跡: `OFF RAD FLUC; ON COD RFSW RADCOD` (決定的な 6 次元、RF あり)。1000 周 (チュートリアルは 1000 周、サンプル集は 5000 周)。

## 何を見る例か
振幅を σ 単位で指定し、水平 (y=1σ) と垂直 (x=1σ) の方向に 1〜100σ を 50 分割して、各点が何周回ったかを運動量偏差 (-1.5%〜+1.5%、0.5% おき) ごとに表示する。

## スクリプト
```
w = EMITX;  sigx = Sqrt[Twiss["BX",1]*w];  sigy = Sqrt[Twiss["BY",1]*w];   ! 分散とエネルギーの分は考えない σ
dp = 0.01/SIGE;                    ! dp/p 1% が何 σE か (= 13.63)
DAPWIDTH = 51;                     ! 51 = 内側を省かず全点計算
zrange = Range[-1.5*dp, 1.6*dp, dp/2];
range  = {{1,100},{1,1},zrange};   ! {{xmin,xmax},{ymin,ymax},{z}}
resx = DynamicApertureSurvey[range, 1000, Output->6];
Print[zrange[i], resx[2,2,i,2], (resx[2,2,i,2]-1)*sigx*1000*2];   ! 運動量、スコア、実振幅 [mm]
```
スタート点は LINE の始点 (`Twiss["BX",1]` = MKSSES、βx = 11.98 m、βy = 5.00 m)。σx = 0.6433 mm、σy = 0.4154 mm。

## 結果の読み方
| dp/p | x スコア | x の DA [mm] | y スコア | y の DA [mm] |
|---|---|---|---|---|
| -1.5% (-20.4σE) | 0 | 0 | 0 | 0 |
| -1% (-13.6σE) | 21 | 25.7 | 13 | 10.0 |
| -0.5% | 19 | 23.2 | 13 | 10.0 |
| 0 | 22 | 27.0 | 18 | 14.1 |
| +0.5% | 19 | 23.2 | 14 | 10.8 |
| +1% | 19 | 23.2 | 14 | 10.8 |
| +1.5% | 0 | 0 | 0 | 0 |

- 運動量 0 で水平 DA は 27 mm (約 42σx)、垂直 14 mm (約 34σy)。(スコア-1)×(100σ/50) が σ 単位の振幅。
- ±1.5% はスコア 0 (振幅によらず早く落ちる)。emit の表の RF バケット高さは dV/P0 = ±1.18% で、±1.5% はバケットの外。

## 物理・解析式
バケット半高さ (emit の表、Bucket height dV/P0 = 0.0118) より外の運動量では、6 次元追跡 (RF あり) で粒子が失われる。これは表の ±1.5% の行 (スコア 0) と合っている。

## gs と pysad の比較
○ 全スコア (x: 0 21 19 22 19 19 0、y: 0 13 13 18 14 14 0) が一致。画面表示の周回数の記号 (`U`、`f`、`K` など) は 1〜2 文字違う (表示のみ。境界の粒子の周回数が 10〜20 周ずれる)。

## 未確認・限界
- 1000 周に短縮 (サンプル集は 5000 周、運動量 ±2%)。DA は周回数で小さくなる。
- 放射 (RAD) と量子励起 (FLUC) は切ってある (乱数に依らない比較のため)。実機の DA には放射の効果がある。
- `DAPWIDTH` を 10 などにすると、連続して生存した外側の点で内側を省略して速くなる (試していない)。
''')

w('fn_TrackParticles_PF', '''
# fn TrackParticles (PF): 1 要素のトラッキングと転送行列
種別: 関数 (入門マニュアル その 8)。''' + COMMON + '''追跡: `OFF RAD FLUC; ON COD RFSW RADCOD`。

## 何を見る例か
`TrackParticles[{始点要素番号, 座標}, 終点要素番号]` で、四極 QF (149 番) を 1 個通し、手計算の転送行列と比べる。座標は
`{{x...},{x'...},{y...},{y'...},{z...},{dP/P...},{生死...}}` の 7 個のリスト (粒子ごとの列)。結果は `{終点番号, 座標}`。

## スクリプト
```
rqk = Sqrt[LINE["K1","QF"]/LINE["L","QF"]];   ! SAD の K1 は積分強度 B'L/(Bρ) なので L で割って平方根
qkmath = {{Cos[rqk*ql], Sin[rqk*ql]/rqk},{-rqk*Sin[rqk*ql], Cos[rqk*ql]}};   ! x 方向 (集束)
qkmatv = {{Cosh[rqk*ql], Sinh[rqk*ql]/rqk},{rqk*Sinh[rqk*ql], Cosh[rqk*ql]}}; ! y 方向 (発散)
after = TrackParticles[{iq, c}, iq+1];       ! c は 1 μm の振幅を 1 つの座標にだけ入れる
```
始点は要素番号で、`iq` (入口) から `iq+1` (出口) まで。1 から始めると 1→1 でリング 1 周分回るので注意。

## 結果の読み方
QF.1 (K1 = 0.74832、L = 0.4 m、$\\sqrt{K} = 1.36777\\ {\\rm m^{-1}}$):

| 要素 | 手計算 | TrackParticles (1 μm) |
|---|---|---|
| x: {{cos, sin/√K},{-√K sin, cos}} | {{0.8540322588, 0.3803413501},{-0.7115421473, 0.8540322588}} | 列 1 = (0.8540322588, -0.7115421473)、列 2 = (0.3803413501, 0.8540322588) |
| y: {{cosh, sinh/√K},{√K sinh, cosh}} | {{1.1534345676, 0.4202559819},{0.7862143932, 1.1534345676}} | 列 3 = (1.1534345676, 0.7862143932)、列 4 = (0.4202559819, 1.1534345676) |
| det | 1.0000000000 (x、y とも) | |

10 桁一致。1 mm の振幅だと x の値が .8540322719 になり (1 μm では .8540322588)、1.3e-8 ずれる。原因は確認していない (振幅に依る高次の項の可能性)。

## 物理・解析式
薄くない四極の行列 $\\begin{pmatrix}\\cos\\sqrt K L & \\sin\\sqrt K L/\\sqrt K\\\\ -\\sqrt K\\sin\\sqrt K L & \\cos\\sqrt K L\\end{pmatrix}$ (x、K>0)。y は cos→cosh、sin→sinh、符号反転。上の表で照合済み。

## gs と pysad の比較
○ 全桁一致。

## 未確認・限界
- 4 粒子を一度に通す形 ($x, x', y, y'$ を別の粒子にする) にはしていない (粒子を 1 つずつ 4 回通した)。
- 四極の端のフリンジ場の扱いは SAD の既定のままで、`FRINGE` フラグは変えていない。
''')

def readme():
    txt = '''
# T1_intro: SAD 入門マニュアルの話題を PF で
元: KEK SAD チュートリアル「SAD の入門マニュアル」(その 1〜8)。NSLS X-ring の例を、PF 2.5 GeV リング (`../lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790) に置き換えた。説明は自分の言葉で書き直してある。

| ファイル (.sad/.gs.txt/.pysad.txt/.md) | 入門マニュアルの節 | 内容 | gs と pysad |
|---|---|---|---|
| `FFS_cell_calc_emit_PF` | その 1・2 | 実行、`cell calc emit`、`SAVE`、周長・チューン・エミッタンス | ○ |
| `FFS_disp_type_PF` | その 2 | `disp`・`type` | ○ |
| `fn_Twiss_LINE_PF` | その 2 | `Twiss[]`・`LINE[]` | ○ |
| `fn_PolynomialFit_PF` | その 3 | チューン対運動量、`DP0`、`PolynomialFit` (SAD の係数の癖あり) | ○ |
| `FFS_draw_OpticsPlot_PF` | その 4 | `draw`、`OpticsPlot` (図の数値と PNG) | ○ (draw は gs 無出力・pysad 未対応) |
| `FFS_FIT_chrom_PF` | その 5 | 色収差補正の6 極フィット (設計値が戻る) | ○ |
| `FFS_FIT_modules_PF` | その 6 | `ChromCorrect[]`・`MatchingNcellTune[]` のモジュール化 | ○ |
| `fn_DynamicApertureSurvey_PF` | その 7 | DA サーベイ (1000 周) | ○ (表示記号に 1〜2 文字の差) |
| `fn_TrackParticles_PF` | その 8 | 四極 1 個の転送行列 | ○ |

- 走らせ方: `py -3.12 tools/runpf.py T1_intro/NAME.sad` (gs は WSL、pysad は `PYSAD_SAD_PACKAGES` つき)。`.sad` は `tools/build_T1.py` が作る。
- 追跡の例 (DA、TrackParticles) は `OFF RAD FLUC; ON COD RFSW RADCOD`。それ以外はチュートリアルと同じ `ON RAD FLUC COD RFSW RADCOD`。
- 入門マニュアルで扱っていない/省いたもの: ファイル読み書き (`OpenWrite`、`Write`、`Close`。出力は画面で見る)、`FitValue` の上下限つきマッチング。
'''
    open(os.path.join(D, 'README.md'), 'w', encoding='utf-8', newline='\n').write(txt.strip('\n') + '\n')

readme()
